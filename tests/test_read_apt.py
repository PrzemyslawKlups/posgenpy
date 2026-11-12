"""Tests for reading AP Suite .apt files through read_pos_to_df.

An .apt has to come back as the very frame its .pos would: same columns, same
dtypes, same atoms in the same order. A real AP Suite export carries many more
sections than Position and Mass, in its own order, so those have to be skipped.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import posgenpy as pgp
from posgenpy import io

EXAMPLE_APT = Path(__file__).parent.parent / "examples" / "alloy_C_irradiated.apt"


def make_tip(n=50, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "ionIdx": np.arange(1, n + 1),
        "x": rng.uniform(-30, 30, n).astype(np.float32),
        "y": rng.uniform(-30, 30, n).astype(np.float32),
        "z": rng.uniform(-100, 0, n).astype(np.float32),
        "mc": rng.uniform(0.5, 120, n).astype(np.float32),
    })


def write_ap_suite_like(pos, path, with_mass=True):
    """An .apt laid out like an AP Suite 6 export: other per-ion sections
    around Mass and Position (Mass first), an int section among them, and a
    two-float one."""

    n = len(pos)
    xyz = pos[["x", "y", "z"]].to_numpy(dtype="<f4")
    bounds = np.column_stack([xyz.min(axis=0), xyz.max(axis=0)]).ravel()
    rng = np.random.default_rng(1)

    header = np.zeros(1, dtype=io._APT_FILE_HEADER)
    header["cSignature"] = b"APT\x00"
    header["iHeaderSize"] = io._APT_FILE_HEADER.itemsize
    header["iHeaderVersion"] = 2
    header["llIonCount"] = n

    def section(name, unit, values, data, header_size=148, int_data=False):
        sec = io._apt_section(name, unit, n, values, header_size)
        if int_data:
            sec["eRecordDataType"] = 1
        return sec.tobytes() + data.tobytes()

    parts = [
        header.tobytes(),
        section("Voltage", "", 1, rng.uniform(3000, 7000, n).astype("<f4")),
        section("Delta Pulse", "", 1, rng.integers(0, 9000, n).astype("<u4"), int_data=True),
    ]
    if with_mass:
        parts.append(section("Mass", "Da", 1, pos["mc"].to_numpy(dtype="<f4")))
    parts += [
        section("Multiplicity", "", 1, np.ones(n, dtype="<u4"), int_data=True),
        section("Detector Coordinates", "mm", 2, rng.normal(size=(n, 2)).astype("<f4")),
        # the six bounds sit between the header and the data, counted in
        # the header's size
        section(
            "Position", "nm", 3, np.concatenate([bounds, xyz.ravel()]),
            header_size=148 + 24,
        ),
        section("Epos ToF", "ns", 1, rng.uniform(100, 900, n).astype("<f4")),
    ]
    Path(path).write_bytes(b"".join(parts))


def test_apt_reads_back_as_the_same_frame_as_its_pos(tmp_path):
    pos = make_tip()
    pgp.write_pos(pos, tmp_path / "tip.pos")
    pgp.write_apt(pos, tmp_path / "tip.apt")

    from_pos = pgp.read_pos_to_df(tmp_path / "tip.pos")
    from_apt = pgp.read_pos_to_df(tmp_path / "tip.apt")

    pd.testing.assert_frame_equal(from_apt, from_pos)
    pd.testing.assert_frame_equal(from_apt, pgp.read_apt_to_df(tmp_path / "tip.apt"))
    pd.testing.assert_frame_equal(from_apt, pos, check_dtype=False)


def test_ap_suite_sections_around_position_and_mass_are_skipped(tmp_path):
    pos = make_tip()
    write_ap_suite_like(pos, tmp_path / "export.apt")

    read = pgp.read_pos_to_df(tmp_path / "export.apt")

    assert list(read.columns) == ["ionIdx", "x", "y", "z", "mc"]
    np.testing.assert_array_equal(
        read[["x", "y", "z", "mc"]].to_numpy(), pos[["x", "y", "z", "mc"]].to_numpy()
    )
    np.testing.assert_array_equal(read["ionIdx"], pos["ionIdx"])


def test_upper_case_extension_is_read_as_apt(tmp_path):
    pos = make_tip()
    pgp.write_apt(pos, tmp_path / "TIP.APT")

    read = pgp.read_pos_to_df(tmp_path / "TIP.APT")

    np.testing.assert_array_equal(read["mc"], pos["mc"])


def test_a_pos_file_named_apt_is_refused(tmp_path):
    pgp.write_pos(make_tip(), tmp_path / "renamed.apt")

    with pytest.raises(ValueError, match="not an AP Suite .apt file"):
        pgp.read_pos_to_df(tmp_path / "renamed.apt")


def test_missing_mass_section_is_refused(tmp_path):
    write_ap_suite_like(make_tip(), tmp_path / "no_mass.apt", with_mass=False)

    with pytest.raises(ValueError, match="no 'Mass' section"):
        pgp.read_apt_to_df(tmp_path / "no_mass.apt")


def test_truncated_file_is_refused(tmp_path):
    pgp.write_apt(make_tip(), tmp_path / "tip.apt")
    raw = (tmp_path / "tip.apt").read_bytes()
    (tmp_path / "cut.apt").write_bytes(raw[:-10])

    with pytest.raises(ValueError, match="truncated"):
        pgp.read_apt_to_df(tmp_path / "cut.apt")


@pytest.mark.skipif(not EXAMPLE_APT.exists(), reason="example .apt not present")
def test_example_ap_suite_file_opens_cleanly():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        read = pgp.read_pos_to_df(EXAMPLE_APT)

    assert len(read) == 4_353_119
    assert list(read.columns) == ["ionIdx", "x", "y", "z", "mc"]
    assert (read[["x", "y", "z", "mc"]].dtypes == np.float32).all()
    assert not read[["x", "y", "z", "mc"]].isna().any().any()
    assert read["mc"].min() > 0

    # the six floats in Position's header are the cloud's bounds
    with open(EXAMPLE_APT, "rb") as f:
        sections = io._apt_sections(
            f, EXAMPLE_APT, io._APT_FILE_HEADER.itemsize, EXAMPLE_APT.stat().st_size
        )
        _, data_offset = sections["Position"]
        f.seek(data_offset - 24)
        bounds = np.fromfile(f, dtype="<f4", count=6)
    xyz = read[["x", "y", "z"]].to_numpy()
    np.testing.assert_array_equal(
        bounds, np.column_stack([xyz.min(axis=0), xyz.max(axis=0)]).ravel()
    )
