"""Synthetic tests pinning each clustering rule to posgen's implementation.

Every assertion here corresponds to a specific behaviour read out of
`posgen/operators/clusterOps.cpp` (or libatomprobe, where noted). They run on
hand-built point sets of a few dozen atoms, so the whole file is fast and
needs no external data - unlike `test_posgen_parity.py`, which compares
against real posgen output.
"""
import numpy as np
import pandas as pd
import pytest

from posgenpy import assign_rrng_to_pos, run_clustering


def make_pos(points, ions):
    """Build a pos-like frame from explicit coordinates and ion labels."""
    df = pd.DataFrame(np.asarray(points, dtype=float), columns=["x", "y", "z"])
    df["mc"] = np.arange(1.0, len(df) + 1.0)
    df["ion"] = list(ions)
    return df


def line(n, spacing, axis=0, offset=(0.0, 0.0, 0.0)):
    """n points evenly spaced along one axis."""
    pts = np.zeros((n, 3))
    pts[:, axis] = np.arange(n) * spacing
    return pts + np.asarray(offset)


# ---------------------------------------------------------------------------
# Core classification: posgen's knn counts the point itself
# ---------------------------------------------------------------------------


def test_order_counts_the_point_itself():
    """`order` is posgen's knn, which counts the point itself.

    A tight clump of exactly 6 core atoms clusters at order=6 (each atom sees
    itself plus 5 others) but not at order=7.
    """
    pos = make_pos(line(6, 0.1), ["Cu"] * 6)

    at_6 = run_clustering(pos, ["Cu"], [], 0.5, 6, 1, 0.0, 0.0, verbose=False)
    at_7 = run_clustering(pos, ["Cu"], [], 0.5, 7, 1, 0.0, 0.0, verbose=False)

    assert (at_6["cluster_id"] != -1).sum() == 6
    assert (at_7["cluster_id"] != -1).sum() == 0


def test_deprecated_k_still_means_the_old_thing():
    """`k` counted only the other neighbours, so k=5 must equal order=6.

    Kept working (with a warning) so existing notebooks keep producing the
    numbers they produced before, rather than silently shifting by one.
    """
    pos = make_pos(line(6, 0.1), ["Cu"] * 6)

    with pytest.warns(DeprecationWarning, match="order"):
        via_k = run_clustering(
            pos, ["Cu"], [], 0.5, N_min=1, envelope_radius=0.0,
            d_erosion=0.0, k=5, verbose=False,
        )
    via_order = run_clustering(pos, ["Cu"], [], 0.5, 6, 1, 0.0, 0.0, verbose=False)

    assert (via_k["cluster_id"] == via_order["cluster_id"]).all()

    with pytest.raises(ValueError, match="not both"):
        run_clustering(
            pos, ["Cu"], [], 0.5, 6, 1, 0.0, 0.0, k=5, verbose=False
        )


def test_core_link_distance_is_strict():
    """Linking uses `dist < Dmax`, not `<=` (clusterOps.cpp ptsInSphere).

    Coordinates are exact binary fractions so the gap between the two groups
    is exactly Dmax, with no floating-point slack to hide the boundary.
    """
    left = np.array([[0.0, 0.0, 0.0], [0.0, 0.0625, 0.0], [0.0, 0.0, 0.0625]])
    right = left + np.array([0.5, 0.0, 0.0])
    pos = make_pos(np.vstack([left, right]), ["Cu"] * 6)

    assert np.linalg.norm(right[0] - left[0]) == 0.5

    out = run_clustering(pos, ["Cu"], [], 0.5, 3, 1, 0.0, 0.0, verbose=False)
    clustered = out[out["cluster_id"] != -1]

    # at exactly Dmax the two groups must NOT link into one cluster
    assert len(clustered) == 6
    assert clustered["cluster_id"].nunique() == 2


# ---------------------------------------------------------------------------
# Envelope: seeded from classified core points only
# ---------------------------------------------------------------------------


def test_border_points_do_not_seed_the_envelope():
    """posgen grows the envelope from classified core points only.

    A border atom (core species, attached to a cluster but not itself a core
    point) must not pull in bulk atoms that no core point can reach.
    """
    core_clump = line(6, 0.05)                        # dense -> all core points
    border = np.array([[0.65, 0.0, 0.0]])             # within Dmax of the clump end
    far_bulk = np.array([[0.80, 0.0, 0.0]])           # only within dbulk of `border`

    pos = make_pos(
        np.vstack([core_clump, border, far_bulk]),
        ["Cu"] * 6 + ["Cu"] + ["Fe"],
    )

    seeded_core = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1, 0.25, 0.0,
        envelope_seeds="core", verbose=False,
    )
    seeded_all = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1, 0.25, 0.0,
        envelope_seeds="all", verbose=False,
    )

    fe = pos["ion"] == "Fe"
    assert (seeded_core.loc[fe.values, "cluster_id"] == -1).all()
    assert (seeded_all.loc[fe.values, "cluster_id"] != -1).all()


def test_envelope_only_recruits_listed_species():
    core_clump = line(6, 0.05)
    neighbours = np.array([[0.1, 0.1, 0.0], [0.1, -0.1, 0.0]])
    pos = make_pos(
        np.vstack([core_clump, neighbours]), ["Cu"] * 6 + ["Fe", "Cr"]
    )

    out = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1, 0.25, 0.0,
        include_unranged=False, verbose=False,
    )

    assert out.loc[out["ion"] == "Fe", "cluster_id"].iloc[0] != -1
    assert out.loc[out["ion"] == "Cr", "cluster_id"].iloc[0] == -1


# ---------------------------------------------------------------------------
# Erosion: reference set membership
# ---------------------------------------------------------------------------


def test_unclustered_core_species_never_erodes():
    """posgen's erosion queries the bulk tree only.

    Core-species atoms live in its core/border trees, so an isolated solute
    atom near the envelope must NOT strip enveloped atoms - even though it is
    itself unclustered.
    """
    core_clump = line(6, 0.05)
    enveloped = np.array([[0.20, 0.0, 0.0]])          # Fe, within dbulk of the clump
    lone_solute = np.array([[0.30, 0.0, 0.0]])        # Cu, 0.1 nm from the Fe

    pos = make_pos(
        np.vstack([core_clump, enveloped, lone_solute]),
        ["Cu"] * 6 + ["Fe"] + ["Cu"],
    )

    out = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1, 0.25, 0.25,
        include_unranged=False, verbose=False,
    )

    fe_row = out[out["ion"] == "Fe"].iloc[0]
    assert fe_row["cluster_id"] != -1, "a lone core-species atom must not erode"


def test_unranged_atom_erodes_only_when_included():
    """`<unranged>` decides whether unranged atoms join posgen's bulk tree.

    The unranged atom sits beyond dbulk of every core point (so it is never
    itself swallowed by the envelope) but within derode of the enveloped Fe.
    """
    core_clump = line(6, 1 / 128, axis=1)             # tight clump on the y axis
    enveloped = np.array([[0.25, 0.0, 0.0]])
    unranged = np.array([[0.5, 0.0, 0.0]])

    pos = make_pos(
        np.vstack([core_clump, enveloped, unranged]),
        ["Cu"] * 6 + ["Fe"] + [None],
    )

    kwargs = dict(envelope_radius=0.375, d_erosion=0.375, verbose=False)

    with_unranged = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1, include_unranged=True, **kwargs
    )
    without = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1, include_unranged=False, **kwargs
    )

    assert (with_unranged.loc[with_unranged["ion"] == "Fe", "cluster_id"] == -1).all()
    assert (without.loc[without["ion"] == "Fe", "cluster_id"] != -1).all()


def test_erosion_distance_is_strict():
    """`dist < derode` erodes; an atom at exactly derode survives."""
    core_clump = line(6, 1 / 128, axis=1)
    enveloped = np.array([[0.25, 0.0, 0.0]])
    matrix = np.array([[0.5, 0.0, 0.0]])              # exactly derode away

    pos = make_pos(
        np.vstack([core_clump, enveloped, matrix]), ["Cu"] * 6 + ["Fe", "Fe"]
    )

    assert np.linalg.norm(matrix[0] - enveloped[0]) == 0.25

    out = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1,
        envelope_radius=0.375, d_erosion=0.25,
        include_unranged=False, verbose=False,
    )

    fe = out[out["ion"] == "Fe"].sort_values("x")["cluster_id"].to_numpy()
    assert fe[0] != -1, "at exactly derode the atom must survive"
    assert fe[1] == -1, "the matrix atom is beyond dbulk and stays unclustered"


def test_erosion_is_not_iterative():
    """Erosion is a single pass over a fixed snapshot (clusterOps.cpp:376-378).

    `inner` and `outer` are both enveloped; `matrix` is not. `outer` is within
    derode of `matrix` so it is stripped, and `inner` is within derode of
    `outer` but not of `matrix`. A cascading implementation would strip
    `inner` on a second pass; posgen does not, because eroded atoms stay
    tagged and never rejoin the reference set.
    """
    core_clump = line(6, 1 / 128, axis=1)
    inner = np.array([[0.25, 0.0, 0.0]])
    outer = np.array([[0.5, 0.0, 0.0]])
    matrix = np.array([[0.75, 0.0, 0.0]])

    pos = make_pos(
        np.vstack([core_clump, inner, outer, matrix]), ["Cu"] * 6 + ["Fe"] * 3
    )

    out = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1,
        envelope_radius=0.625, d_erosion=0.375,
        include_unranged=False, verbose=False,
    )

    fe = out[out["ion"] == "Fe"].sort_values("x")["cluster_id"].to_numpy()
    assert fe[0] != -1, "inner atom must survive: erosion must not cascade"
    assert fe[1] == -1, "outer atom should have been eroded by the matrix atom"
    assert fe[2] == -1, "matrix atom was never enveloped"


# ---------------------------------------------------------------------------
# Size clipping
# ---------------------------------------------------------------------------


def test_n_min_counts_core_atoms_by_default():
    """posgen's sizeclip counts core+border ions, not enveloped bulk."""
    core_clump = line(6, 0.05)
    bulk = line(10, 0.05, axis=1, offset=(0.0, 0.1, 0.0))

    pos = make_pos(np.vstack([core_clump, bulk]), ["Cu"] * 6 + ["Fe"] * 10)

    # 6 core atoms, plus enveloped Fe: passes N_min=7 only if bulk is counted
    core_counted = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 5, 7, 0.25, 0.0,
        n_min_count="core", include_unranged=False, verbose=False,
    )
    all_counted = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 5, 7, 0.25, 0.0,
        n_min_count="all", include_unranged=False, verbose=False,
    )

    assert (core_counted["cluster_id"] != -1).sum() == 0
    assert (all_counted["cluster_id"] != -1).sum() > 0


# ---------------------------------------------------------------------------
# Range assignment
# ---------------------------------------------------------------------------


def test_range_assignment_is_first_match_wins():
    """libatomprobe's RangeFile::getIonID returns on its first hit."""
    pos = pd.DataFrame({"x": [0.0], "y": [0.0], "z": [0.0], "mc": [29.0]})
    rrng = pd.DataFrame(
        {
            "mc_begin": [28.9, 28.5],
            "mc_end": [29.1, 29.5],
            "ion": ["Fm", "Ni"],
            "volume": [0.0, 0.01],
            "color": ["008000", "00CC00"],
            "formula": [{"Fm": 1}, {"Ni": 1}],
        }
    )

    assert assign_rrng_to_pos(pos, rrng)["ion"].iloc[0] == "Fm"


def test_range_bounds_are_inclusive():
    pos = pd.DataFrame(
        {"x": [0.0] * 2, "y": [0.0] * 2, "z": [0.0] * 2, "mc": [28.9, 29.1]}
    )
    rrng = pd.DataFrame(
        {
            "mc_begin": [28.9],
            "mc_end": [29.1],
            "ion": ["Fm"],
            "volume": [0.0],
            "color": ["008000"],
            "formula": [{"Fm": 1}],
        }
    )

    assert (assign_rrng_to_pos(pos, rrng)["ion"] == "Fm").all()


# ---------------------------------------------------------------------------
# Bookkeeping
# ---------------------------------------------------------------------------


def test_unranged_atoms_are_excluded_from_output():
    """posgen's splitUnrangedClustered pulls them out before writing stats."""
    core_clump = line(6, 0.05)
    unranged = np.array([[0.1, 0.1, 0.0]])

    pos = make_pos(np.vstack([core_clump, unranged]), ["Cu"] * 6 + [None])

    out = run_clustering(
        pos, ["Cu"], ["Fe"], 0.5, 6, 1, 0.25, 0.0,
        include_unranged=True, verbose=False,
    )

    assert len(out) == 6
    assert out["ion"].notna().all()
    assert out.attrs["n_unranged_clustered"] == 1


def test_rejects_dclassify_that_differs_from_dmax():
    pos = make_pos(line(6, 0.05), ["Cu"] * 6)

    with pytest.raises(ValueError, match="dclassify"):
        run_clustering(
            pos, ["Cu"], [], 0.5, 6, 1, 0.25, 0.0, dclassify=0.3, verbose=False
        )
