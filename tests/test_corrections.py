"""Tests for the post-clustering corrections the final run applies.

The final run is the one whose numbers get reported, and its corrections run in
an order that matters: Fm deconvolution before anything is derived from the
counts, then the size columns, then edge detection. These pin the pieces of that
order which can silently produce plausible-but-wrong numbers.

Synthetic data only, so the file is fast and needs no external datasets.
"""
import numpy as np
import pandas as pd
import pytest

import posgenpy as pgp


def make_atoms(ions, cluster_ids):
    """A clustered-atoms frame of the shape run_clustering returns."""
    return pd.DataFrame(
        {
            "ionIdx": np.arange(1, len(ions) + 1),
            "x": np.arange(len(ions), dtype=float),
            "y": 0.0,
            "z": 0.0,
            "ion": list(ions),
            "cluster_id": list(cluster_ids),
        }
    )


# ---------------------------------------------------------------------------
# summarise_unclustered
# ---------------------------------------------------------------------------


def test_unclustered_size_is_the_unclustered_atom_count():
    """get_cluster_metrics reads `size` off this for the matrix atom count, so
    it has to stay the plain count whatever else is added alongside it."""
    atoms = make_atoms(["Fe", "Fe", "Ni", "Cu"], [-1, 0, -1, 0])

    assert pgp.summarise_unclustered(atoms)["size"].iloc[0] == 2


def test_unclustered_counts_are_credited_per_element():
    """Same rule as analyse_clusters: a C2-ranged atom is two carbons.

    Without this the matrix composition disagrees with the cluster composition
    it is compared against.
    """
    atoms = make_atoms(["C2", "C", "Fe", "Fe"], [-1, -1, -1, 0])
    formula_map = {"C2": {"C": 2}, "C": {"C": 1}, "Fe": {"Fe": 1}}

    summary = pgp.summarise_unclustered(atoms, ion_formula_map=formula_map)

    assert summary["size"].iloc[0] == 4          # three ions, four atoms
    assert summary["C"].iloc[0] == 3             # ...crediting two carbons for C2
    assert summary["Fe"].iloc[0] == 1


def test_unclustered_columns_cover_every_element_in_the_range_file():
    """A fixed schema, so frames from different runs concatenate and an absent
    element reads as 0 rather than going missing."""
    atoms = make_atoms(["Fe", "Fe"], [-1, -1])
    formula_map = {"Fe": {"Fe": 1}, "Ni": {"Ni": 1}, "Fm": {"Fm": 1}}

    summary = pgp.summarise_unclustered(atoms, ion_formula_map=formula_map)

    assert {"Fe", "Ni", "Fm"} <= set(summary.columns)
    assert summary["Ni"].iloc[0] == 0


def test_unclustered_summary_survives_a_frame_with_no_ion_column():
    atoms = pd.DataFrame({"cluster_id": [-1, 0, -1]})

    assert pgp.summarise_unclustered(atoms)["size"].iloc[0] == 2


# ---------------------------------------------------------------------------
# Fm deconvolution
# ---------------------------------------------------------------------------


def test_fm_deconvolution_conserves_atoms_and_empties_the_peak():
    """The 29 Da peak is an Fe/Ni overlap being split, not atoms being created.

    Applied to the matrix as well as the clusters, which is why
    summarise_unclustered has to emit these columns at all.
    """
    stats = pd.DataFrame(
        {"cluster_id": [0, 1], "Fe": [100, 50], "Ni": [20, 5], "Fm": [40, 3]}
    )

    out = pgp.deconvolute_Fm_clusters(stats)

    assert (out["Fm"] == 0).all()
    assert out[["Fe", "Ni"]].to_numpy().sum() == stats[["Fe", "Ni", "Fm"]].to_numpy().sum()
    assert (out["Ni"] >= stats["Ni"]).all()
    assert (out["Fe"] >= stats["Fe"]).all()


def test_fm_deconvolution_caps_ni_at_its_natural_abundance():
    """At most Ni/0.32*0.68 of the peak can be Ni; the rest must go to Fe."""
    stats = pd.DataFrame({"cluster_id": [0], "Fe": [0], "Ni": [10], "Fm": [1000]})

    out = pgp.deconvolute_Fm_clusters(stats)

    assert out["Ni"].iloc[0] == 10 + round(10 / 0.32 * 0.68)
    assert out["Fe"].iloc[0] == 1000 - round(10 / 0.32 * 0.68)


def test_fm_deconvolution_refuses_a_frame_without_the_columns():
    """Silently skipping would leave the Fm atoms uncounted in every downstream
    composition."""
    with pytest.raises(ValueError, match="Ni"):
        pgp.deconvolute_Fm_clusters(pd.DataFrame({"cluster_id": [0], "Fm": [5]}))


# ---------------------------------------------------------------------------
# The size columns the old workflow reported on
# ---------------------------------------------------------------------------


def test_cluster_size_excludes_the_matrix_element():
    stats = pd.DataFrame(
        {"cluster_id": [0], "Cu": [10], "Ni": [5], "Fe": [200]}
    )

    with_matrix = pgp.calculate_cluster_size(stats, ["Cu", "Ni"], ["Fe"], [])
    without = pgp.calculate_cluster_size(stats, ["Cu", "Ni"], ["Fe"], ["Fe"])

    assert with_matrix.iloc[0] == 215
    assert without.iloc[0] == 15


# ---------------------------------------------------------------------------
# `size` counts atoms, on both sides of the cluster/matrix split
# ---------------------------------------------------------------------------

C2_FORMULA_MAP = {"C2": {"C": 2}, "Cu": {"Cu": 1}, "Fe": {"Fe": 1}}


def test_cluster_size_counts_atoms_and_size_core_counts_ions():
    """A C2 ion is two atoms of volume: `size` has to say 4 here, matching the
    element columns. `size_core` stays posgen's ion count."""
    atoms = make_atoms(["C2", "Cu", "Fe"], [0, 0, 0])

    stats = pgp.analyse_clusters(
        atoms, ion_formula_map=C2_FORMULA_MAP, core_species=["Cu", "C2"]
    )

    assert stats["size"].iloc[0] == 4
    assert stats["size"].iloc[0] == stats[["C", "Cu", "Fe"]].sum(axis=1).iloc[0]
    assert stats["size_core"].iloc[0] == 2


def test_cluster_size_without_a_formula_map_counts_ions():
    atoms = make_atoms(["C2", "Cu", "Fe"], [0, 0, 0])

    assert pgp.analyse_clusters(atoms)["size"].iloc[0] == 3


def test_control_run_size_counts_atoms_like_the_real_run():
    """Real and randomised clusters are compared on `size`, so both have to
    count the same unit."""
    from posgenpy.cluster import _summarise_sizes

    atoms = make_atoms(["C2", "Cu", "Fe", "Cu"], [0, 0, 0, -1])

    sizes = _summarise_sizes(atoms, ["Cu", "C2"], ion_formula_map=C2_FORMULA_MAP)
    real = pgp.analyse_clusters(
        atoms, ion_formula_map=C2_FORMULA_MAP, core_species=["Cu", "C2"]
    )

    assert sizes["size"].iloc[0] == real["size"].iloc[0] == 4
    assert sizes["size_core"].iloc[0] == real["size_core"].iloc[0] == 2
    assert _summarise_sizes(atoms, ["Cu", "C2"])["size"].iloc[0] == 3


def test_tip_volume_counts_every_ranged_atom():
    """Cluster and matrix totals add up to the tip's atoms, so a molecular ion
    in the matrix adds all of its atoms to the tip volume."""
    atoms = make_atoms(["C2", "Cu", "Fe", "C2", "Fe", "Cu"], [0, 0, 0, -1, -1, -1])

    stats = pgp.analyse_clusters(atoms, ion_formula_map=C2_FORMULA_MAP)
    matrix = pgp.summarise_unclustered(atoms, ion_formula_map=C2_FORMULA_MAP)
    m = pgp.get_cluster_metrics(
        "r", "loc", "s", stats, matrix, ["Cu"], ["Fe"],
        dmax=0.5, order=5, n_min=1,
        detector_efficiency=0.5, atomic_density_per_nm3=10.0,
    )

    assert m["tip volume nm3"] * 0.5 * 10.0 == pytest.approx(8)    # 6 ions, 8 atoms
    assert m["volume fraction"] == pytest.approx(4 / 8)


# ---------------------------------------------------------------------------
# get_cluster_metrics
# ---------------------------------------------------------------------------


def _metrics_100_clustered_40_fe_900_matrix():
    """100 clustered atoms (40 of them Fe) and 900 in the matrix."""
    clusters = pd.DataFrame(
        {"size": [60, 40], "Fe": [20, 20], "is_edge_cluster": [False, False]}
    )
    matrix = pd.DataFrame({"size": [900]})
    return pgp.get_cluster_metrics(
        "r", "loc", "s", clusters, matrix, ["Ni"], ["Fe"],
        dmax=0.5, order=5, n_min=10, excluded_Fe=True,
    )


def test_volume_fraction_no_fe_divides_by_all_ranged_atoms():
    """Taking the clustered Fe out of the clusters hands it to the matrix; the
    total stays at the 1000 ranged atoms. Adding it to the denominator on top
    of that counts it twice (60 / 1040)."""
    m = _metrics_100_clustered_40_fe_900_matrix()

    assert m["volume fraction no Fe"] == pytest.approx(0.06)
    assert m["volume fraction"] == pytest.approx(0.10)


def test_no_fe_fraction_matches_no_fe_volume_over_tip_volume():
    """The function's own outputs have to agree with each other."""
    m = _metrics_100_clustered_40_fe_900_matrix()

    assert m["volume fraction no Fe"] == pytest.approx(
        m["clustered volume no Fe nm3"] / m["tip volume nm3"]
    )


# ---------------------------------------------------------------------------
# The pipeline as the final run assembles it
# ---------------------------------------------------------------------------


def test_n_min_inside_run_clustering_clips_on_core_atoms_before_the_envelope():
    """The final run passes n_min to run_clustering rather than filtering the
    stats afterwards, because posgen clips core+border counts *before*
    enveloping - a cluster below n_min must not first collect a bulk shell that
    lifts it back over the threshold.
    """
    # four core atoms in a tight line, plus bulk neighbours that the envelope
    # would otherwise add
    positions = np.zeros((10, 3))
    positions[:4, 0] = np.arange(4) * 0.1          # core
    positions[4:, 0] = np.arange(6) * 0.1 + 0.05   # bulk, interleaved

    atoms = pd.DataFrame(positions, columns=["x", "y", "z"])
    atoms["mc"] = 1.0
    atoms["ion"] = ["Cu"] * 4 + ["Fe"] * 6

    kept = pgp.run_clustering(
        atoms, ["Cu"], ["Fe"], 0.5, order=3, N_min=4,
        envelope_radius=0.3, d_erosion=0.0, verbose=False,
    )
    clipped = pgp.run_clustering(
        atoms, ["Cu"], ["Fe"], 0.5, order=3, N_min=5,
        envelope_radius=0.3, d_erosion=0.0, verbose=False,
    )

    # 4 core atoms: survives n_min=4, and the envelope then adds bulk on top
    assert (kept["cluster_id"] != -1).sum() > 4
    # at n_min=5 the cluster is gone before the envelope could pad it out
    assert (clipped["cluster_id"] != -1).sum() == 0
