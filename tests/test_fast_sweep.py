"""Tests for the accelerated (order x dmax) grid sweep.

The claim `fast_sweep` rests on is that the core-classification and linking
steps can be answered for every parameter cell from one kNN query and one
radius query per dataset. These tests pin that claim: the same core atoms, the
same core points, the same partition as `run_clustering` (up to border atoms
two clusters could equally claim), and results that the existing sweep
consumers read unchanged.

Synthetic data only - a few thousand atoms - so the file is fast and needs no
external datasets, like `test_cluster_semantics.py` and unlike
`test_posgen_parity.py`.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.cluster import DBSCAN

import posgenpy as pgp
from posgenpy.fast_sweep import (
    CoreNeighbourhood,
    _apply_n_min,
    _sorted_sizes,
    core_code_from_mc,
    extract_core_atoms,
    load_sweep_size_table,
    relabelled_core_atoms,
    run_fast_grid_sweep,
    select_sweep,
    split_core_code,
    summarise_grid,
)

CORE_IONS = ["Cu", "Ni"]
BULK_IONS = ["Fe"]


def make_rrng():
    """Three non-overlapping ranges, one per ion."""
    return pd.DataFrame(
        {
            "mc_begin": [10.0, 20.0, 30.0],
            "mc_end": [11.0, 21.0, 31.0],
            "ion": ["Cu", "Ni", "Fe"],
            "volume": [0.0118, 0.0118, 0.0118],
            "color": ["FF6600", "00CC00", "FF00FF"],
            "formula": [{"Cu": 1}, {"Ni": 1}, {"Fe": 1}],
        }
    )


def make_tip(n_matrix=4000, n_clusters=12, cluster_size=40, seed=0):
    """A matrix of Fe with planted Cu/Ni blobs, as a pos-like frame.

    Blob spread and spacing are chosen so that at the dmax values tested some
    blobs merge and some border atoms sit between two of them - i.e. the cases
    where implementations are allowed to disagree actually occur.
    """
    rng = np.random.default_rng(seed)

    matrix = rng.uniform(-10, 10, size=(n_matrix, 3))
    centres = rng.uniform(-8, 8, size=(n_clusters, 3))
    blobs = np.concatenate(
        [centre + rng.normal(0, 0.35, size=(cluster_size, 3)) for centre in centres]
    )

    positions = np.concatenate([matrix, blobs])
    mc = np.concatenate(
        [
            np.full(n_matrix, 30.5),  # Fe
            rng.choice([10.5, 20.5], size=len(blobs)),  # Cu / Ni
        ]
    )

    shuffle = rng.permutation(len(positions))
    positions, mc = positions[shuffle], mc[shuffle]

    return pd.DataFrame(
        {
            "ionIdx": np.arange(1, len(positions) + 1),
            "x": positions[:, 0],
            "y": positions[:, 1],
            "z": positions[:, 2],
            "mc": mc,
        }
    )


def make_percolating_tip(
    n_atoms=20000, solute_fraction=0.15, n_clusters=0, cluster_size=80, seed=0
):
    """A uniform solid solution dense enough that the solute network percolates
    at the dmax under test. With the default `n_clusters=0` nothing is planted,
    so every cluster found is background - which is the case `make_tip` cannot
    produce and the one that breaks a count-based real-vs-random comparison.

    `n_clusters` plants blobs in that background as well, giving a tip with
    genuine precipitates *and* a percolating matrix: real signal that the
    cluster counts cannot see past the giant cluster to find.
    """
    rng = np.random.default_rng(seed)

    positions = rng.uniform(-10, 10, size=(n_atoms, 3))
    n_solute = int(n_atoms * solute_fraction)
    mc = np.full(n_atoms, 30.5)  # Fe
    mc[:n_solute] = rng.choice([10.5, 20.5], size=n_solute)  # Cu / Ni
    rng.shuffle(mc)

    if n_clusters:
        centres = rng.uniform(-8, 8, size=(n_clusters, 3))
        blobs = np.concatenate(
            [centre + rng.normal(0, 0.3, size=(cluster_size, 3)) for centre in centres]
        )
        positions = np.concatenate([positions, blobs])
        mc = np.concatenate([mc, rng.choice([10.5, 20.5], size=len(blobs))])

    return pd.DataFrame(
        {
            "ionIdx": np.arange(1, len(positions) + 1),
            "x": positions[:, 0],
            "y": positions[:, 1],
            "z": positions[:, 2],
            "mc": mc,
        }
    )


@pytest.fixture(scope="module")
def tip():
    return make_tip(), make_rrng()


# ---------------------------------------------------------------------------
# Dataset extraction
# ---------------------------------------------------------------------------


def test_core_atoms_match_the_full_ranging_path(tip):
    """`extract_core_atoms` must select exactly what assign_rrng_to_pos does."""
    pos_df, rrng_df = tip

    index, code = extract_core_atoms(pos_df, rrng_df, CORE_IONS)

    labelled = pgp.assign_rrng_to_pos(pos_df, rrng_df)
    expected_index = np.flatnonzero(labelled["ion"].isin(CORE_IONS).to_numpy())

    assert np.array_equal(index, expected_index)
    assert list(np.array(CORE_IONS)[code]) == list(labelled["ion"].to_numpy()[index])


def test_overlapping_ranges_are_first_match_wins(tip):
    """libatomprobe returns on its first hit, so an earlier range wins.

    The shortcut must not quietly change which atoms count as core species.
    """
    pos_df, rrng_df = tip

    overlapping = pd.concat(
        [
            rrng_df,
            # covers Cu's range again, but later in the file
            pd.DataFrame(
                {
                    "mc_begin": [10.0], "mc_end": [11.0], "ion": ["Ni"],
                    "volume": [0.0118], "color": ["00CC00"], "formula": [{"Ni": 1}],
                }
            ),
        ],
        ignore_index=True,
    )

    index, code = extract_core_atoms(pos_df, overlapping, CORE_IONS)
    labelled = pgp.assign_rrng_to_pos(pos_df, overlapping)

    assert np.array_equal(
        index, np.flatnonzero(labelled["ion"].isin(CORE_IONS).to_numpy())
    )
    # every atom in the overlapped range is still Cu, not Ni
    assert list(np.array(CORE_IONS)[code]) == list(labelled["ion"].to_numpy()[index])


def test_relabelling_equals_shuffling_mc_and_reranging(tip):
    """The identity the whole relabelling implementation rests on.

    Ranging is elementwise in mc, so ranging a permuted mc column gives the same
    per-atom labels as permuting the ranged ones: R(mc[pi]) == R(mc)[pi]. That
    is what makes shuffling the code array *the* standard relabelling rather
    than a cheaper stand-in for it - the only thing skipped is walking every
    mass range over every atom a second time.
    """
    pos_df, rrng_df = tip
    mc = pos_df["mc"].to_numpy()

    permutation = np.random.default_rng(11).permutation(len(pos_df))

    # permute the labels
    permuted_labels = core_code_from_mc(mc, rrng_df, CORE_IONS)[permutation]
    # ...against ranging the permuted mc column from scratch
    reranged = core_code_from_mc(mc[permutation], rrng_df, CORE_IONS)

    assert np.array_equal(permuted_labels, reranged)

    # and it agrees with the package's own ranging of a shuffled pos frame
    shuffled = pos_df.copy()
    shuffled["mc"] = mc[permutation]
    labelled = pgp.assign_rrng_to_pos(shuffled, rrng_df)
    index, code = split_core_code(reranged)

    assert np.array_equal(
        index, np.flatnonzero(labelled["ion"].isin(CORE_IONS).to_numpy())
    )
    assert list(np.array(CORE_IONS)[code]) == list(labelled["ion"].to_numpy()[index])


def test_relabelling_preserves_the_core_label_multiset(tip):
    """A control dataset is a permutation of the labels, not a resampling.

    The randomised runs are only a fair baseline if they hold exactly as many
    of each core ion as the real tip, spread over the whole tip.
    """
    pos_df, rrng_df = tip
    full_code = core_code_from_mc(pos_df["mc"].to_numpy(), rrng_df, CORE_IONS)
    _, real_code = split_core_code(full_code)

    index, code = relabelled_core_atoms(full_code, np.random.default_rng(0))

    assert len(index) == len(real_code)
    assert len(np.unique(index)) == len(index)
    assert index.max() < len(pos_df)
    assert np.array_equal(np.bincount(code), np.bincount(real_code))
    # the real tip's codes are not disturbed
    assert np.array_equal(
        full_code, core_code_from_mc(pos_df["mc"].to_numpy(), rrng_df, CORE_IONS)
    )


def test_relabelling_is_reproducible_from_the_seed(tip):
    pos_df, rrng_df = tip
    full_code = core_code_from_mc(pos_df["mc"].to_numpy(), rrng_df, CORE_IONS)

    first = relabelled_core_atoms(full_code, np.random.default_rng(7))
    again = relabelled_core_atoms(full_code, np.random.default_rng(7))

    assert np.array_equal(first[0], again[0])
    assert np.array_equal(first[1], again[1])


# ---------------------------------------------------------------------------
# The shortcut itself
# ---------------------------------------------------------------------------


ORDERS = [3, 5, 7]
DMAX_VALUES = [0.4, 0.6, 0.8]


@pytest.fixture(scope="module")
def neighbourhood(tip):
    pos_df, rrng_df = tip
    index, _ = extract_core_atoms(pos_df, rrng_df, CORE_IONS)
    positions = pos_df[["x", "y", "z"]].to_numpy()[index]
    return CoreNeighbourhood(positions, ORDERS, max(DMAX_VALUES)), positions


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("dmax", DMAX_VALUES)
def test_core_points_match_dbscan_exactly(neighbourhood, order, dmax):
    """The order statistic must reproduce DBSCAN's core-point test at every cell.

    This is the whole shortcut: one kNN query answers `min_samples` for every
    order, and re-thresholding it answers every dmax. If this ever drifts, the
    clusters drift with it.
    """
    nb, positions = neighbourhood

    fitted = DBSCAN(eps=float(np.nextafter(dmax, 0.0)), min_samples=order).fit(positions)
    expected = np.zeros(len(positions), dtype=bool)
    expected[fitted.core_sample_indices_] = True

    assert np.array_equal(nb.core_distance(order) < dmax, expected)


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("dmax", DMAX_VALUES)
def test_same_atoms_clustered_as_dbscan(neighbourhood, order, dmax):
    """Which atoms end up in *a* cluster is not open to interpretation."""
    nb, positions = neighbourhood

    fitted = DBSCAN(eps=float(np.nextafter(dmax, 0.0)), min_samples=order).fit(positions)

    assert np.array_equal(nb.labels(order, dmax) != -1, fitted.labels_ != -1)


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("dmax", DMAX_VALUES)
def test_partition_matches_run_clustering_within_border_ambiguity(
    tip, neighbourhood, order, dmax
):
    """Same clusters as the full engine, bar atoms two clusters could claim.

    A border atom within dmax of core atoms of two different clusters is given
    to the nearest one here and to whichever the traversal reaches first by
    sklearn; nothing else may differ.
    """
    pos_df, rrng_df = tip
    nb, _ = neighbourhood

    labelled = pgp.assign_rrng_to_pos(pos_df, rrng_df)
    full = pgp.run_clustering(
        labelled, CORE_IONS, BULK_IONS, dmax, order, 1,
        envelope_radius=dmax / 2, d_erosion=dmax / 2, verbose=False,
    )

    index, _ = extract_core_atoms(pos_df, rrng_df, CORE_IONS)
    full_labels = pd.Series(-1, index=pos_df.index, dtype=np.int64)
    core_rows = full.loc[full["ion"].isin(CORE_IONS)]
    full_labels.loc[core_rows.index] = core_rows["cluster_id"].to_numpy()

    from posgenpy.fast_sweep import _partition_agreement

    agreement = _partition_agreement(
        full_labels.to_numpy()[index], nb.labels(order, dmax)
    )
    ambiguous = nb.ambiguous_border_atoms(order, dmax)

    assert agreement["n_atoms_full"] == agreement["n_atoms_fast"]
    assert agreement["n_differing_atoms"] <= ambiguous


def test_envelope_and_erosion_cannot_change_the_core_size_distribution(tip):
    """Why skipping the envelope is safe for the sweep stage.

    `optimise_n_min` and the real-vs-random plots threshold posgen's
    `Cluster Size(core)`. The envelope only ever adds bulk-tree atoms and
    erosion only ever removes them, so that number is settled before either
    runs - which is what lets this path skip both.
    """
    pos_df, rrng_df = tip
    labelled = pgp.assign_rrng_to_pos(pos_df, rrng_df)

    with_envelope = pgp.run_clustering(
        labelled, CORE_IONS, BULK_IONS, 0.6, 5, 1,
        envelope_radius=0.3, d_erosion=0.3, verbose=False,
    )
    without_envelope = pgp.run_clustering(
        labelled, CORE_IONS, BULK_IONS, 0.6, 5, 1,
        envelope_radius=0.0, d_erosion=0.0, verbose=False,
    )

    def core_sizes(clustered):
        core_only = clustered.loc[
            (clustered["cluster_id"] != -1) & clustered["ion"].isin(CORE_IONS)
        ]
        return np.sort(core_only.groupby("cluster_id").size().to_numpy())

    assert np.array_equal(core_sizes(with_envelope), core_sizes(without_envelope))


def make_dumbbell(blob_size=60, neck_size=6, separation=3.0, seed=0):
    """Two dense blobs joined by a thin line of atoms.

    At a low order the neck atoms pass the knn test, so the whole thing is one
    cluster. Raise the order and the sparse neck stops qualifying as core, and
    the two blobs are no longer connected.
    """
    rng = np.random.default_rng(seed)

    left = rng.normal(0, 0.3, size=(blob_size, 3))
    right = rng.normal(0, 0.3, size=(blob_size, 3)) + [separation, 0, 0]

    neck = np.zeros((neck_size, 3))
    neck[:, 0] = np.linspace(0.6, separation - 0.6, neck_size)

    return np.concatenate([left, neck, right])


def test_raising_the_order_splits_a_cluster():
    """A cluster held together by a thin neck must come apart at higher order.

    Only the dmax axis is walked monotonically here; orders are recomputed from
    scratch, precisely so that a split is reported as a split. If any
    incremental union crept into the linking step, this is the test that would
    catch it - at order 10 the neck atoms fail the knn test and the two blobs
    stop being connected.
    """
    positions = make_dumbbell()
    nb = CoreNeighbourhood(positions, [3, 10], 0.5)

    joined = nb.labels(3, 0.5)
    split = nb.labels(10, 0.5)

    assert len(np.unique(joined[joined != -1])) == 1
    assert len(np.unique(split[split != -1])) == 2

    # the two halves are the two blobs, not an arbitrary cut
    left_of_neck = positions[:, 0] < 1.5
    for cluster_id in np.unique(split[split != -1]):
        members = split == cluster_id
        assert len(np.unique(left_of_neck[members])) == 1

    # and the full engine says the same
    pos_df = pd.DataFrame(positions, columns=["x", "y", "z"])
    pos_df["mc"] = 10.5
    pos_df["ionIdx"] = np.arange(1, len(pos_df) + 1)

    for order, expected in ((3, 1), (10, 2)):
        full = pgp.run_clustering(
            pgp.assign_rrng_to_pos(pos_df, make_rrng()), CORE_IONS, BULK_IONS,
            0.5, order, 1, envelope_radius=0.25, d_erosion=0.25, verbose=False,
        )
        clustered = full.loc[full["cluster_id"] != -1, "cluster_id"]
        assert clustered.nunique() == expected


def test_a_split_is_not_an_artefact_of_sweeping_order_upwards():
    """Cells are independent: the same cell gives the same answer whichever
    order the grid happens to visit it in."""
    positions = make_dumbbell()
    nb = CoreNeighbourhood(positions, [3, 10], 0.5)

    ascending = [nb.labels(order, 0.5) for order in (3, 10)]
    descending = [nb.labels(order, 0.5) for order in (10, 3)][::-1]

    for first, second in zip(ascending, descending):
        assert np.array_equal(first, second)


def test_dmax_only_ever_adds_atoms_to_clusters(neighbourhood):
    """Raising dmax cannot un-cluster an atom, at a fixed order.

    The property the whole "check the next parameter as you go" idea rests on:
    both the core test and the link test are monotone in dmax.
    """
    nb, _ = neighbourhood

    clustered = [nb.labels(5, dmax) != -1 for dmax in DMAX_VALUES]

    for smaller, larger in zip(clustered, clustered[1:]):
        assert np.all(larger[smaller])


def test_dmax_beyond_the_precomputed_radius_is_refused(neighbourhood):
    """Pairs past max_dmax were never collected - answering anyway would
    silently under-link, so it must raise instead."""
    nb, _ = neighbourhood

    with pytest.raises(ValueError, match="max_dmax"):
        nb.labels(5, max(DMAX_VALUES) + 0.1)


def test_order_counts_the_point_itself():
    """Same rule as the full engine: posgen's Order includes the point itself."""
    positions = np.zeros((6, 3))
    positions[:, 0] = np.arange(6) * 0.1

    nb = CoreNeighbourhood(positions, [6, 7], 0.5)

    assert (nb.labels(6, 0.5) != -1).sum() == 6
    assert (nb.labels(7, 0.5) != -1).sum() == 0


def test_distance_test_is_strict():
    """Two atoms exactly dmax apart do not link, matching posgen's `<`."""
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])

    nb = CoreNeighbourhood(positions, [2], 0.6)

    assert (nb.labels(2, 0.5) != -1).sum() == 0
    assert (nb.labels(2, 0.5000001) != -1).sum() == 2


def test_n_min_drops_small_clusters(neighbourhood):
    nb, _ = neighbourhood

    labels = nb.labels(5, 0.6)
    sizes = _sorted_sizes(labels)

    clipped = _apply_n_min(labels, 10)

    assert np.array_equal(_sorted_sizes(clipped), sizes[sizes >= 10])
    # labels stay a dense 0..n-1 range after the clip
    surviving = clipped[clipped != -1]
    assert set(np.unique(surviving)) == set(range(len(np.unique(surviving))))


def test_no_core_species_gives_empty_clusters():
    nb = CoreNeighbourhood(np.zeros((0, 3)), [5], 0.5)

    assert len(nb.labels(5, 0.5)) == 0
    assert nb.ambiguous_border_atoms(5, 0.5) == 0


# ---------------------------------------------------------------------------
# The sweep, and what reads its output
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def grid(tip, tmp_path_factory):
    pos_df, rrng_df = tip
    output_dir = tmp_path_factory.mktemp("fast_grid")

    return run_fast_grid_sweep(
        pos_df, rrng_df, CORE_IONS,
        orders=[5, 7], dmax_values=[0.4, 0.6],
        bulk_ions=BULK_IONS, n_relabel_runs=3,
        output_dir=output_dir, prefix="test", verbose=False,
    ), output_dir


def test_grid_covers_every_cell(grid):
    manifest, _ = grid

    assert len(manifest) == 4
    assert {(e["order"], e["dmax"]) for e in manifest} == {
        (5, 0.4), (5, 0.6), (7, 0.4), (7, 0.6)
    }
    assert all(e["pipeline"] == "core_only" for e in manifest)


def test_sizes_written_are_core_counts(grid):
    """`size` and `size_core` agree here: with no envelope every atom counted
    is a core-species one. That is what makes size_metric="core" exact."""
    manifest, _ = grid

    stats = pd.read_parquet(manifest[0]["real_cluster_stats_path"])
    assert (stats["size"] == stats["size_core"]).all()
    assert (stats[CORE_IONS].sum(axis=1) == stats["size"]).all()


def test_select_sweep_feeds_the_existing_sweep_functions(grid):
    """The point of the file layout: optimise_n_min and the plots take this
    manifest with no changes."""
    manifest, _ = grid

    dmax_sweep = select_sweep(manifest, order=7)
    assert [e["swept_value"] for e in dmax_sweep] == [0.4, 0.6]
    assert all(e["swept_param"] == "dmax" for e in dmax_sweep)

    n_min, real, random = pgp.optimise_n_min(dmax_sweep, n_min_values=range(1, 30))
    assert len(n_min) == len(real) == len(random) == 2

    order_sweep = select_sweep(manifest, dmax=0.4)
    assert [e["swept_value"] for e in order_sweep] == [5, 7]

    figure, _ = pgp.plot_real_clusters_across_swept_param(dmax_sweep, "dmax", [2, 5])
    assert figure is not None


def test_select_sweep_rejects_an_ambiguous_request(grid):
    manifest, _ = grid

    with pytest.raises(ValueError, match="exactly one"):
        select_sweep(manifest)
    with pytest.raises(ValueError, match="exactly one"):
        select_sweep(manifest, order=7, dmax=0.4)
    with pytest.raises(ValueError, match="covers orders"):
        select_sweep(manifest, order=99)


def test_size_table_carries_the_cluster_id_and_run(grid):
    """What the distribution plots hover on: every point knows its own cluster."""
    manifest, _ = grid

    table = load_sweep_size_table(manifest)

    assert {"order", "dmax", "cluster_id", "size_core", "run", "run_type"} <= set(
        table.columns
    )
    assert set(table["run_type"]) <= {"real", "random"}
    assert set(table.loc[table["run_type"] == "random", "run"]) == {
        "random_000", "random_001", "random_002"
    }

    real_only = load_sweep_size_table(manifest, include_random=False)
    assert (real_only["run_type"] == "real").all()


def test_real_atoms_carry_cluster_ids_matching_the_stats(grid):
    """The round trip the notebook relies on: a cluster_id read off a size
    distribution selects the same atoms in the 3D view."""
    manifest, _ = grid

    entry = manifest[0]
    atoms = pd.read_parquet(entry["real_clustered_atoms_path"])
    stats = pd.read_parquet(entry["real_cluster_stats_path"])

    clustered = atoms.loc[atoms["cluster_id"] != -1]
    assert set(clustered["cluster_id"]) == set(stats["cluster_id"])

    for cluster_id, size in stats.set_index("cluster_id")["size"].head(5).items():
        assert (clustered["cluster_id"] == cluster_id).sum() == size


def test_resume_reuses_finished_cells(tip, grid):
    """A rerun must not recompute - long sweeps get interrupted."""
    pos_df, rrng_df = tip
    manifest, output_dir = grid

    stamps = {
        entry["real_cluster_stats_path"]: __import__("os").path.getmtime(
            entry["real_cluster_stats_path"]
        )
        for entry in manifest
    }

    run_fast_grid_sweep(
        pos_df, rrng_df, CORE_IONS,
        orders=[5, 7], dmax_values=[0.4, 0.6],
        bulk_ions=BULK_IONS, n_relabel_runs=3,
        output_dir=output_dir, prefix="test", verbose=False,
    )

    for path, stamp in stamps.items():
        assert __import__("os").path.getmtime(path) == stamp


def test_control_runs_match_the_full_pipeline(tip, tmp_path):
    """The control runs, not just the real one, must agree with the slow path.

    Both paths are given the same permutations, so they are clustering the
    *same* relabelled datasets rather than merely similar ones - which is what
    makes this a comparison of the two algorithms and not of two random draws.
    """
    pos_df, rrng_df = tip

    # a low order and a wide dmax, so the randomised datasets actually produce
    # clusters to compare rather than a handful of empty runs
    dmax, order, n_runs = 0.8, 3, 3
    rng = np.random.default_rng(5)
    shuffles = [rng.permutation(len(pos_df)) for _ in range(n_runs)]

    fast = run_fast_grid_sweep(
        pos_df, rrng_df, CORE_IONS, orders=[order], dmax_values=[dmax],
        bulk_ions=BULK_IONS, n_relabel_runs=n_runs, cached_shuffles=shuffles,
        output_dir=tmp_path / "fast", prefix="fast", verbose=False,
    )[0]

    slow = pgp.run_relabel_sweep(
        pos_df, rrng_df, CORE_IONS, BULK_IONS,
        dmax=dmax, order=order, n_min=1,
        envelope_radius=dmax / 2, d_erosion=dmax / 2,
        n_relabel_runs=n_runs, run_name="slow", output_dir=tmp_path / "slow",
        cached_shuffles=shuffles,
    )

    full_code = core_code_from_mc(pos_df["mc"].to_numpy(), rrng_df, CORE_IONS)

    for run in range(n_runs):
        fast_sizes = pd.read_parquet(
            Path(fast["random_dir"]) / f"run_{run:03d}.parquet"
        )["size_core"]
        slow_sizes = pd.read_parquet(
            Path(slow["random_dir"]) / f"run_{run:03d}.parquet"
        )["size_core"]

        # identical control dataset on both sides
        index, _ = split_core_code(full_code[shuffles[run]])
        nb = CoreNeighbourhood(
            pos_df[["x", "y", "z"]].to_numpy()[index], [order], dmax
        )
        ambiguous = nb.ambiguous_border_atoms(order, dmax)

        assert len(fast_sizes) > 0, "a control run with no clusters proves nothing"

        if ambiguous == 0:
            # nothing either implementation is free to decide differently, so
            # the size distributions must be equal outright
            assert np.array_equal(
                np.sort(fast_sizes.to_numpy()), np.sort(slow_sizes.to_numpy())
            )
        else:
            assert abs(len(fast_sizes) - len(slow_sizes)) <= ambiguous
            assert abs(fast_sizes.sum() - slow_sizes.sum()) <= ambiguous


def test_cached_shuffles_must_cover_every_control_run(tip, tmp_path):
    pos_df, rrng_df = tip

    with pytest.raises(ValueError, match="cached_shuffles"):
        run_fast_grid_sweep(
            pos_df, rrng_df, CORE_IONS, orders=[5], dmax_values=[0.5],
            n_relabel_runs=3,
            cached_shuffles=[np.arange(len(pos_df))],
            output_dir=tmp_path, prefix="short", verbose=False,
        )


def test_drop_in_sweep_matches_run_parameter_sweep_call_for_call(tip, tmp_path):
    """The migration path: swap the function name, change nothing else.

    An existing notebook calls run_parameter_sweep with these arguments and
    feeds the result to optimise_n_min and the plot_* functions. All of that has
    to keep working verbatim, or the change is not the one-liner it claims.
    """
    pos_df, rrng_df = tip

    order_manifest = pgp.run_fast_parameter_sweep(
        pos_df, rrng_df, CORE_IONS, BULK_IONS,
        swept_param="order", swept_values=[5, 7],
        fixed_dmax=0.5, n_min=2, n_relabel_runs=2,
        output_dir=tmp_path, prefix="tip", verbose=False,
    )
    dmax_manifest = pgp.run_fast_parameter_sweep(
        pos_df, rrng_df, CORE_IONS, BULK_IONS,
        swept_param="dmax", swept_values=[0.4, 0.5],
        fixed_order=5, n_min=2, n_relabel_runs=2,
        output_dir=tmp_path, prefix="tip", verbose=False,
    )

    assert [e["swept_value"] for e in order_manifest] == [5, 7]
    assert all(e["swept_param"] == "order" for e in order_manifest)
    assert [e["swept_value"] for e in dmax_manifest] == [0.4, 0.5]
    assert all(e["swept_param"] == "dmax" for e in dmax_manifest)

    # the two sweeps must not overwrite each other's manifest file
    assert (tmp_path / "tip_order_sweep_manifest.json").exists()
    assert (tmp_path / "tip_dmax_sweep_manifest.json").exists()

    # ...and the cell they share is written once, not twice
    shared = [e for e in order_manifest if e["order"] == 5 and e["dmax"] == 0.5]
    assert len(shared) == 1
    assert shared[0]["run_name"] == "tip_order_5_dmax_0p5"

    # everything downstream takes them unchanged
    for manifest, label in ((order_manifest, "Order"), (dmax_manifest, "dmax")):
        n_min, real, random = pgp.optimise_n_min(
            manifest, n_min_values=range(1, 30), verbose=False
        )
        assert len(n_min) == 2

        figure, _ = pgp.plot_real_clusters_across_swept_param(manifest, label, [2, 5])
        assert figure is not None
        figure, _ = pgp.plot_real_cluster_ratio_across_swept_param(
            manifest, label, [2, 5]
        )
        assert figure is not None
        figure, _ = pgp.plot_cluster_composition_across_swept_param_relative(
            manifest, label, CORE_IONS, n_min
        )
        assert figure is not None


def test_drop_in_sweep_rejects_the_same_calls_run_parameter_sweep_does(tip, tmp_path):
    pos_df, rrng_df = tip

    with pytest.raises(ValueError, match="swept_param"):
        pgp.run_fast_parameter_sweep(
            pos_df, rrng_df, CORE_IONS, BULK_IONS,
            swept_param="nmin", swept_values=[1], output_dir=tmp_path,
        )
    with pytest.raises(ValueError, match="fixed_dmax"):
        pgp.run_fast_parameter_sweep(
            pos_df, rrng_df, CORE_IONS, BULK_IONS,
            swept_param="order", swept_values=[5], output_dir=tmp_path,
        )
    with pytest.raises(ValueError, match="fixed_order"):
        pgp.run_fast_parameter_sweep(
            pos_df, rrng_df, CORE_IONS, BULK_IONS,
            swept_param="dmax", swept_values=[0.5], output_dir=tmp_path,
        )


def test_drop_in_sweep_warns_about_the_arguments_it_cannot_honour(tip, tmp_path):
    """envelope_radius/d_erosion/dclassify describe stages this path skips.
    Accepting them silently would let a notebook believe it enveloped."""
    pos_df, rrng_df = tip

    with pytest.warns(UserWarning, match="no envelope or erosion"):
        pgp.run_fast_parameter_sweep(
            pos_df, rrng_df, CORE_IONS, BULK_IONS,
            swept_param="dmax", swept_values=[0.5], fixed_order=5,
            envelope_radius=0.25, d_erosion=0.25,
            n_relabel_runs=0, output_dir=tmp_path, prefix="warned",
        )


def test_summarise_grid_covers_every_cell(grid):
    """The 2D read of the grid: one row per cell, from the same threshold rule
    optimise_n_min applies to a 1D sweep."""
    manifest, _ = grid

    summary = summarise_grid(manifest, n_min_values=range(1, 60))

    assert len(summary) == len(manifest)
    assert set(zip(summary["order"], summary["dmax"])) == {
        (5, 0.4), (5, 0.6), (7, 0.4), (7, 0.6)
    }
    assert {
        "optimal_n_min", "n_real_clusters", "n_random_clusters",
        "true_positive_ratio", "n_min_reached",
    } <= set(summary.columns)

    # the per-cell answer must be the one the 1D path gives for that cell
    dmax_sweep = select_sweep(manifest, order=7)
    n_min, real, _ = pgp.optimise_n_min(
        dmax_sweep, n_min_values=range(1, 60), verbose=False
    )
    at_order_7 = summary.loc[summary["order"] == 7].sort_values("dmax")

    assert list(at_order_7["optimal_n_min"]) == n_min
    assert list(at_order_7["n_real_clusters"]) == real


def test_grid_heatmap_carries_every_cell_column(grid):
    manifest, _ = grid

    summary = summarise_grid(manifest, n_min_values=range(1, 60))
    figure = pgp.plot_grid_heatmap(summary, value="optimal_n_min")

    assert figure.data[0].z.shape == (2, 2)
    hover = figure.data[0].customdata[0][0]
    assert "n_real_clusters" in hover and "true_positive_ratio" in hover

    with pytest.raises(ValueError, match="not a column"):
        pgp.plot_grid_heatmap(summary, value="nope")


def test_heatmap_blanks_cells_that_never_reached_the_threshold():
    """A cell that never got there is reported at the largest n_min tried.
    Plotted as a value it is an outlier that flattens the scale over every cell
    that did work, so it has to be a hole, not a number."""
    summary = pd.DataFrame(
        {
            "order": [5, 5, 7, 7],
            "dmax": [0.4, 0.6, 0.4, 0.6],
            "optimal_n_min": [4, 6, 8, 199],
            "n_real_clusters": [40, 30, 20, 1],
            "n_min_reached": [True, True, True, False],
        }
    )

    masked = np.asarray(
        pgp.plot_grid_heatmap(summary, value="optimal_n_min").data[0].z, dtype=float
    )
    assert np.isnan(masked).sum() == 1
    assert np.nanmax(masked) == 8

    unmasked = np.asarray(
        pgp.plot_grid_heatmap(summary, value="optimal_n_min", mask_unreached=False)
        .data[0].z,
        dtype=float,
    )
    assert np.nanmax(unmasked) == 199


# ---------------------------------------------------------------------------
# The randomised control percolating, and the sweep noticing
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def percolating_sweep(tmp_path_factory):
    """An order sweep at a dmax where the control's solute network spans the
    tip at low order and is sub-critical at high order.

    No planted clusters at all, so anything either dataset finds is background:
    the whole range of the diagnostic in one manifest.
    """
    pos_df = make_percolating_tip()
    output_dir = tmp_path_factory.mktemp("percolating")

    return run_fast_grid_sweep(
        pos_df, make_rrng(), CORE_IONS,
        orders=[1, 2, 3, 4, 5, 6, 7, 8], dmax_values=[1.4],
        bulk_ions=BULK_IONS, n_relabel_runs=2,
        output_dir=output_dir, prefix="perc", verbose=False,
    )


@pytest.fixture(scope="module")
def percolating_sweep_with_clusters(tmp_path_factory):
    """The same, with genuine precipitates planted in the percolating matrix."""
    pos_df = make_percolating_tip(solute_fraction=0.12, n_clusters=8)
    output_dir = tmp_path_factory.mktemp("percolating_blobs")

    return run_fast_grid_sweep(
        pos_df, make_rrng(), CORE_IONS,
        orders=[1, 2, 3, 4, 5, 6, 7, 8], dmax_values=[1.3],
        bulk_ions=BULK_IONS, n_relabel_runs=2,
        output_dir=output_dir, prefix="blobs", verbose=False,
    )


def test_percolation_report_flags_a_uniform_solid_solution(percolating_sweep):
    """At an order this permissive the control is one tip-spanning cluster.
    Every cluster count downstream is then arithmetic on a broken control, and
    this is the only thing in the package that says so."""

    report = pgp.percolation_report(percolating_sweep, warn=False)
    at_order_3 = report.loc[report["swept_value"] == 3].iloc[0]

    assert at_order_3["verdict"] == "percolating"
    assert at_order_3["giant_fraction"] > 0.5
    assert at_order_3["denominator"] == "n_core_atoms"
    # the mean hides how bad a single run gets; the max is why it is reported
    assert (
        at_order_3["largest_random_cluster_max"]
        >= at_order_3["largest_random_cluster"]
    )


def test_percolation_report_is_ok_for_planted_clusters(grid):
    """The case the sweep is for: a dilute matrix with real precipitates in it.
    A guard that fires here would be worse than no guard."""
    manifest, _ = grid

    for order in (5, 7):
        report = pgp.percolation_report(select_sweep(manifest, dmax=0.4), warn=False)
        assert set(report["verdict"]) == {"ok"}


def test_percolation_report_warns_so_it_survives_verbose_false(percolating_sweep):
    with pytest.warns(UserWarning, match="randomised control"):
        pgp.percolation_report(percolating_sweep)


def test_raising_the_order_walks_the_control_through_the_threshold(percolating_sweep):
    """Pins the mechanism, not just the flag: the order filter is what takes
    the control back below the percolation threshold, and the diagnostic has to
    follow it the whole way down."""

    report = pgp.percolation_report(percolating_sweep, warn=False).sort_values(
        "swept_value"
    )
    fractions = report["giant_fraction"].to_numpy()

    # removing core points can only ever shrink the giant cluster
    assert np.all(np.diff(fractions) <= 1e-12)
    assert fractions[0] > 0.30 and fractions[-1] < 0.10
    assert list(report["verdict"])[0] == "percolating"
    assert list(report["verdict"])[-1] == "ok"


def test_n_min_report_marks_an_unreachable_threshold(percolating_sweep):
    """Real clusters running out is not the search succeeding, but the legacy
    3-tuple reports the two identically."""

    with pytest.warns(UserWarning, match="never reached"):
        report = pgp.n_min_report(percolating_sweep, n_min_values=range(1, 400))

    unreachable = report.loc[report["status"] == "unreachable"]
    assert not unreachable.empty
    assert unreachable["n_min"].isna().all()
    assert (unreachable["n_real_clusters"] == 0).all()

    # ...and the legacy value it would have returned is still there, unchanged
    legacy, _, _ = pgp.optimise_n_min(
        percolating_sweep, n_min_values=range(1, 400), verbose=False
    )
    assert list(report["n_min_returned"]) == legacy


def test_n_min_report_marks_an_exhausted_range(percolating_sweep):
    with pytest.warns(UserWarning, match="never reached"):
        report = pgp.n_min_report(percolating_sweep, n_min_values=range(1, 5))

    exhausted = report.loc[report["status"] == "range_exhausted"]
    assert not exhausted.empty
    assert exhausted["n_min"].isna().all()
    assert (exhausted["n_min_returned"] == 4).all()


def test_n_min_report_is_the_same_search_optimise_n_min_runs(grid):
    """The wrapper must not drift from the report it wraps."""
    dmax_sweep = select_sweep(grid[0], order=7)

    report = pgp.n_min_report(dmax_sweep, n_min_values=range(1, 60), warn=False)
    n_min, real, random = pgp.optimise_n_min(
        dmax_sweep, n_min_values=range(1, 60), verbose=False
    )

    assert list(report["n_min_returned"]) == n_min
    assert list(report["n_real_clusters"]) == real
    assert list(report["n_random_clusters"]) == random


def test_summarise_grid_does_not_report_an_unreachable_cell_as_reached(
    percolating_sweep,
):
    """The regression: `n_min_reached` used to be re-derived from the returned
    integer (`n_min < max(n_min_values)`), which is true for a cell that ran out
    of real clusters. Such a cell was reported as a pass, at a plausible-looking
    optimal_n_min, with nothing but a NaN in a column no one filters on."""

    with pytest.warns(UserWarning):
        summary = summarise_grid(percolating_sweep, n_min_values=range(1, 400))

    unreachable = summary.loc[summary["n_min_status"] == "unreachable"]
    assert not unreachable.empty
    assert (unreachable["optimal_n_min"] < 400).all()  # the old heuristic's "pass"
    assert not unreachable["n_min_reached"].any()

    # the column is now exactly the status, not a re-derivation of it
    assert (summary["n_min_reached"] == (summary["n_min_status"] == "reached")).all()


def test_summarise_grid_carries_the_percolation_verdict(percolating_sweep, grid):
    """Both failure modes have to be visible in the one table a cell is
    actually picked from."""

    with pytest.warns(UserWarning):
        summary = summarise_grid(percolating_sweep, n_min_values=range(1, 400))

    assert {"percolation_verdict", "giant_fraction", "true_positive_ratio_atoms"} <= set(
        summary.columns
    )
    assert "percolating" in set(summary["percolation_verdict"])

    clean = summarise_grid(grid[0], n_min_values=range(1, 60), warn=False)
    assert set(clean["percolation_verdict"]) == {"ok"}


def test_percolation_report_runs_on_a_manifest_without_n_core_atoms(
    percolating_sweep,
):
    """n_core_atoms was added after the sweeps that need this diagnostic most
    were written. Without it the denominator falls back to the solute the
    control actually clustered, which is labelled as such because it inflates
    the fraction at high order."""

    older = [
        {k: v for k, v in entry.items() if k != "n_core_atoms"}
        for entry in percolating_sweep
    ]

    report = pgp.percolation_report(older, warn=False)

    assert set(report["denominator"]) == {"clustered_solute"}
    assert report["giant_fraction"].notna().all()
    # the percolating end is still caught, which is the point of the fallback
    assert report.loc[report["swept_value"] == 1, "verdict"].iloc[0] == "percolating"


def test_the_fallback_denominator_does_not_flag_a_handful_of_atoms():
    """`clustered_solute` is a small number when the control clusters almost
    nothing, so a 12-atom "giant" cluster is 100% of it. The absolute floor is
    what stops that reading as a percolating network."""
    from posgenpy.sweep import _percolation_rows

    data = {
        1: {
            "real_sizes": pd.Series([20, 8]),
            "random_sizes": pd.Series([12, 4, 3, 11, 5, 2]),
            "random_runs": [pd.Series([12, 4, 3]), pd.Series([11, 5, 2])],
            "n_runs": 2,
        }
    }
    rows = _percolation_rows([{"swept_value": 1}], data, 0.10, 0.30)

    assert rows[0]["giant_fraction"] > 0.5  # 11.5 of a ~18-atom denominator
    assert rows[0]["verdict"] == "ok"


def test_the_atom_ratio_sees_what_the_cluster_counts_cannot(
    percolating_sweep_with_clusters,
):
    """The tell that the physics is fine and the count-based metric is not.

    A percolating control contributes exactly one cluster however much solute
    it swallowed, so `true_positive_ratio` can sit at or above the threshold
    where there is no signal at all. Counted in atoms the same comparison
    collapses to zero across the percolating orders and only becomes positive
    once the control is sub-critical and the planted clusters are what is left.
    """

    report = pgp.n_min_report(
        percolating_sweep_with_clusters, n_min_values=range(1, 400), warn=False
    ).sort_values("swept_value")
    percolation = pgp.percolation_report(
        percolating_sweep_with_clusters, warn=False
    ).sort_values("swept_value")

    merged = report.merge(percolation[["swept_value", "verdict"]], on="swept_value")

    percolating = merged.loc[merged["verdict"] == "percolating"]
    assert not percolating.empty
    assert (percolating["true_positive_ratio_atoms"] <= 0.01).all()

    at_the_top = merged.iloc[-1]
    assert at_the_top["verdict"] == "ok"
    assert at_the_top["true_positive_ratio_atoms"] > 0.1


def test_the_sweep_plots_shade_the_percolating_swept_values(percolating_sweep):
    """A reader of the current graph has no way to see that the left of the
    x-axis is nonsense."""

    figure, ax = pgp.plot_real_cluster_ratio_across_swept_param(
        percolating_sweep, "Order", [2, 5]
    )
    assert any("percolat" in str(text.get_text()) for text in ax.get_legend().get_texts())

    _, clean_ax = pgp.plot_real_cluster_ratio_across_swept_param(
        percolating_sweep, "Order", [2, 5], flag_percolation=False
    )
    assert not any(
        "percolat" in str(text.get_text()) for text in clean_ax.get_legend().get_texts()
    )

    figure, ax = pgp.plot_atom_ratio_across_swept_param(percolating_sweep, "Order")
    assert figure is not None
