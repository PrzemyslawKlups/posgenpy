"""Tests for exporting the final clusters as posgen's .pos files and as .apt.

These files are what a user takes into AP Suite or IVAS, so they have to say
the same as cluster_stats: the same atoms, the same cluster ids (posgen's
offset of 1 included), and nothing of the tip lost between the cluster and
unclustered files.
"""
import numpy as np
import pandas as pd
import pytest

import posgenpy as pgp


def read_apt(path):
    """Minimal .apt reader, walking the file the way NOMAD's
    ifes_apt_tc_data_modeling apt6_reader does: a 540-byte header, then
    148-byte section headers each followed by llByteCount bytes of data -
    "Position" with six bound floats in between."""

    raw = open(path, "rb").read()
    header = np.frombuffer(raw, dtype=pgp.io._APT_FILE_HEADER, count=1)[0]
    offset = int(header["iHeaderSize"])
    sections = {}

    while offset < len(raw):
        section = np.frombuffer(
            raw, dtype=pgp.io._APT_SECTION_HEADER, count=1, offset=offset
        )[0]
        name = "".join(chr(c) for c in section["wcSectionType"] if c)
        offset += 148
        extra = None
        if name == "Position":
            extra = np.frombuffer(raw, dtype="<f4", count=6, offset=offset)
            offset += 24
        n_bytes = int(section["llByteCount"])
        data = np.frombuffer(raw, dtype="<f4", count=n_bytes // 4, offset=offset)
        sections[name] = (section, extra, data)
        offset += n_bytes

    return header, sections


def make_tip(n=40, seed=0):
    rng = np.random.default_rng(seed)
    pos = pd.DataFrame({
        "ionIdx": np.arange(1, n + 1),
        "x": rng.uniform(-5, 5, n).astype(np.float32),
        "y": rng.uniform(-5, 5, n).astype(np.float32),
        "z": rng.uniform(0, 20, n).astype(np.float32),
        "mc": rng.choice([27.97, 29.0, 58.0], n).astype(np.float32),
    })
    # run_clustering keeps ranged atoms only: drop a few "unranged" ones
    clustered = pos.iloc[5:].copy()
    clustered["cluster_id"] = -1
    clustered.loc[clustered.index[:6], "cluster_id"] = 3
    clustered.loc[clustered.index[10:14], "cluster_id"] = 0
    return pos, clustered


def test_write_apt_round_trips_positions_and_masses(tmp_path):
    pos, _ = make_tip()
    pgp.write_apt(pos, tmp_path / "tip.apt")
    header, sections = read_apt(tmp_path / "tip.apt")

    assert bytes(header["cSignature"]) == b"APT"
    assert header["iHeaderSize"] == 540 and header["iHeaderVersion"] == 2
    assert header["llIonCount"] == len(pos)
    assert "".join(chr(c) for c in header["wcFileName"] if c) == "tip.apt"
    assert set(sections) == {"Position", "Mass"}

    section, bounds, data = sections["Position"]
    assert section["iHeaderSize"] == 172 and section["iRecordSize"] == 12
    assert "".join(chr(c) for c in section["wcDataUnit"] if c) == "nm"
    xyz = pos[["x", "y", "z"]].to_numpy()
    np.testing.assert_array_equal(data.reshape(-1, 3), xyz)
    np.testing.assert_array_equal(
        bounds, np.column_stack([xyz.min(axis=0), xyz.max(axis=0)]).ravel()
    )

    section, _, data = sections["Mass"]
    assert section["iHeaderSize"] == 148 and section["llRecordCount"] == len(pos)
    np.testing.assert_array_equal(data, pos["mc"].to_numpy())

    # the creation time is a FILETIME of about now, not 0 or a Unix stamp
    assert header["ftCreationTime"] > 130_000_000_000_000_000


def test_write_apt_refuses_an_empty_file(tmp_path):
    pos, _ = make_tip()
    with pytest.raises(ValueError, match="at least one ion"):
        pgp.write_apt(pos.iloc[:0], tmp_path / "empty.apt")


def test_cluster_files_read_back_as_posgen_output(tmp_path):
    pos, clustered = make_tip()
    written = pgp.export_cluster_pos_files(clustered, pos, tmp_path, "tip", apt=True)

    assert set(written) == {
        "cluster", "clusterID", "unclustered",
        "cluster_apt", "clusterID_apt", "unclustered_apt",
    }

    # the package's own reader of real posgen output takes them as such, with
    # posgen's offset of 1 between the file and cluster_stats' ids
    back = pgp.load_posgen_clustering(written["cluster"], written["clusterID"])
    assert back.attrs["row_aligned"]
    in_clusters = clustered.loc[clustered["cluster_id"] != -1]
    assert sorted(back["cluster_id"] - 1) == sorted(in_clusters["cluster_id"])
    # grouped by cluster
    assert list(back["cluster_id"]) == sorted(back["cluster_id"])

    # cluster + unclustered is the whole tip, unranged atoms included
    unclustered = pgp.read_pos_to_df(written["unclustered"])
    assert len(back) + len(unclustered) == len(pos)
    everything = pd.concat([back[["x", "y", "z", "mc"]], unclustered[["x", "y", "z", "mc"]]])
    pd.testing.assert_frame_equal(
        everything.sort_values(["x", "y", "z"]).reset_index(drop=True),
        pos[["x", "y", "z", "mc"]].sort_values(["x", "y", "z"]).reset_index(drop=True),
    )
    assert set(pos["mc"].iloc[:5]) <= set(unclustered["mc"])

    # and the .apt copies say the same as the .pos files
    _, sections = read_apt(written["clusterID_apt"])
    np.testing.assert_array_equal(
        sections["Mass"][2], pgp.read_pos_to_df(written["clusterID"])["mc"].to_numpy()
    )


def test_cluster_files_refuse_atoms_from_another_tip(tmp_path):
    pos, clustered = make_tip()
    other, _ = make_tip(n=30, seed=1)
    with pytest.raises(ValueError, match="same tip"):
        pgp.export_cluster_pos_files(clustered, other, tmp_path, "tip")


def test_no_clusters_still_writes_the_pos_files(tmp_path):
    pos, clustered = make_tip()
    clustered["cluster_id"] = -1
    written = pgp.export_cluster_pos_files(clustered, pos, tmp_path, "tip", apt=True)

    assert pgp.read_pos_to_df(written["cluster"]).empty
    assert "cluster_apt" not in written and "unclustered_apt" in written
