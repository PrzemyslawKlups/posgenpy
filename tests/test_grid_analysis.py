"""Tests for the grid-reading layer: per-cell cluster tables, cluster
separation, split/merge tracking, and the plots the grid notebook is built on.

`test_fast_sweep.py` pins the sweep itself - that every cell's partition is the
one the full engine would give. These pin what is read *off* those cells: the
numbers the heatmaps colour, the distributions the ridgelines draw, and the
atom-by-atom comparison between neighbouring cells.

Synthetic data only, so the file is fast and needs no external datasets.
"""
import numpy as np
import pandas as pd
import pytest

import posgenpy as pgp
from posgenpy.corrections import _alpha_shape_mesh
from posgenpy.fast_sweep import _clip_small_clusters, _n_min_lookup

CORE_IONS = ["Cu", "Ni"]
BULK_IONS = ["Fe"]


def make_rrng():
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


def make_tip(n_matrix=3000, n_clusters=10, cluster_size=40, seed=0):
    """An Fe matrix with planted Cu/Ni blobs, plus one dumbbell.

    The dumbbell is a pair of blobs joined by a thin neck: it is one cluster at
    a low order and two at a high one, so a grid over this tip has a real split
    in it for `track_cluster_splits` to find.
    """
    rng = np.random.default_rng(seed)

    matrix = rng.uniform(-10, 10, size=(n_matrix, 3))
    centres = rng.uniform(-8, 8, size=(n_clusters, 3))
    blobs = np.concatenate(
        [centre + rng.normal(0, 0.3, size=(cluster_size, 3)) for centre in centres]
    )

    left = rng.normal(0, 0.3, size=(cluster_size, 3)) + [9.0, 9.0, 9.0]
    right = rng.normal(0, 0.3, size=(cluster_size, 3)) + [11.5, 9.0, 9.0]
    neck = np.zeros((6, 3))
    neck[:, 0] = np.linspace(9.5, 11.0, 6)
    neck[:, 1] = 9.0
    neck[:, 2] = 9.0
    dumbbell = np.concatenate([left, neck, right])

    positions = np.concatenate([matrix, blobs, dumbbell])
    mc = np.concatenate(
        [
            np.full(n_matrix, 30.5),                              # Fe
            rng.choice([10.5, 20.5], size=len(blobs) + len(dumbbell)),
        ]
    )

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
def grid(tmp_path_factory):
    pos_df, rrng_df = make_tip(), make_rrng()
    output_dir = tmp_path_factory.mktemp("grid_analysis")

    manifest = pgp.run_fast_grid_sweep(
        pos_df, rrng_df, CORE_IONS,
        orders=[3, 5, 8], dmax_values=[0.4, 0.6],
        bulk_ions=BULK_IONS, n_min=1, n_relabel_runs=3,
        output_dir=output_dir, prefix="grid", verbose=False,
    )
    return manifest


@pytest.fixture(scope="module")
def summary(grid):
    return pgp.summarise_grid(grid, n_min_values=range(1, 40))


# ---------------------------------------------------------------------------
# Cluster separation
# ---------------------------------------------------------------------------


def brute_force_separation(stats, radius_column="equivalent_spherical_radius"):
    """The definitions, written out directly, for the fast paths to match."""
    xyz = stats[["x", "y", "z"]].to_numpy(dtype=float)
    radii = stats[radius_column].to_numpy(dtype=float)

    distance = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    ratio = distance / (radii[:, None] + radii[None, :])
    np.fill_diagonal(distance, np.inf)
    np.fill_diagonal(ratio, np.inf)

    return distance.min(axis=1), ratio.min(axis=1)


def test_normalised_separation_is_exact_not_a_knn_approximation():
    """The k-nearest search is only a first guess: a big cluster further away
    can crowd a small one more than its nearest neighbour does. The refinement
    pass has to find it."""
    rng = np.random.default_rng(0)
    n = 200

    stats = pd.DataFrame(
        {
            "cluster_id": np.arange(n),
            "x": rng.uniform(-20, 20, n),
            "y": rng.uniform(-20, 20, n),
            "z": rng.uniform(-20, 20, n),
            # a wide radius spread is what makes the two answers differ
            "equivalent_spherical_radius": rng.uniform(0.2, 6.0, n),
        }
    )

    result = pgp.find_knn_cluster_distance(stats)
    expected_distance, expected_ratio = brute_force_separation(stats)

    assert np.allclose(result["closest_cluster_distance"], expected_distance)
    assert np.allclose(result["closest_cluster_distance_norm"], expected_ratio)

    # the partner ids must be the ones those values were taken from
    partner = result.set_index("cluster_id")["closest_cluster_id_norm"]
    for cid, other in partner.items():
        centres = stats.set_index("cluster_id")[["x", "y", "z"]]
        radii = stats.set_index("cluster_id")["equivalent_spherical_radius"]
        gap = np.linalg.norm(centres.loc[cid] - centres.loc[other])
        assert np.isclose(
            gap / (radii.loc[cid] + radii.loc[other]),
            result.set_index("cluster_id").loc[cid, "closest_cluster_distance_norm"],
        )


def test_overlapping_clusters_come_out_below_one():
    """The number the heatmap is read on: below 1 the equivalent spheres of the
    two clusters intersect, which is one precipitate reported as two."""
    stats = pd.DataFrame(
        {
            "cluster_id": [0, 1],
            "x": [0.0, 1.0], "y": [0.0, 0.0], "z": [0.0, 0.0],
            "equivalent_spherical_radius": [0.8, 0.8],
        }
    )

    result = pgp.find_knn_cluster_distance(stats)

    assert np.isclose(result["closest_cluster_distance"].min(), 1.0)
    assert np.isclose(result["closest_cluster_distance_norm"].min(), 1 / 1.6)
    assert (result["closest_cluster_distance_norm"] < 1).all()


def test_separation_without_radii_or_neighbours_is_not_invented():
    no_radius = pd.DataFrame(
        {"cluster_id": [0, 1], "x": [0.0, 1.0], "y": [0.0, 0.0], "z": [0.0, 0.0]}
    )
    result = pgp.find_knn_cluster_distance(no_radius)
    assert "closest_cluster_distance_norm" not in result.columns

    alone = pd.DataFrame(
        {"cluster_id": [0], "x": [0.0], "y": [0.0], "z": [0.0],
         "equivalent_spherical_radius": [1.0]}
    )
    result = pgp.find_knn_cluster_distance(alone)
    assert result["closest_cluster_distance"].isna().all()
    assert result["closest_cluster_distance_norm"].isna().all()


# ---------------------------------------------------------------------------
# What the heatmaps are coloured by
# ---------------------------------------------------------------------------


def test_summary_cluster_metrics_are_the_cells_own_numbers(grid, summary):
    """Every added column must be reproducible from that cell's stats file, at
    that cell's n_min - the whole point being that one heatmap cell and one
    parquet file say the same thing."""

    assert {
        "n_clusters_all", "max_cluster_size",
        "min_cluster_distance", "min_cluster_distance_norm",
    } <= set(summary.columns)

    for entry in grid:
        row = summary.loc[
            (summary["order"] == entry["order"])
            & np.isclose(summary["dmax"], entry["dmax"])
        ].squeeze()

        stats = pd.read_parquet(entry["real_cluster_stats_path"])
        assert row["n_clusters_all"] == len(stats)

        kept = stats.loc[stats["size_core"] >= row["optimal_n_min"]]
        assert row["n_real_clusters"] == len(kept)
        assert row["max_cluster_size"] == kept["size_core"].max()

        distance, ratio = brute_force_separation(kept)
        assert np.isclose(row["min_cluster_distance"], distance.min())
        assert np.isclose(row["min_cluster_distance_norm"], ratio.min())


def test_summary_metrics_can_be_switched_off(grid):
    lean = pgp.summarise_grid(grid, n_min_values=range(1, 40), cluster_metrics=False)

    assert "n_clusters_all" not in lean.columns
    assert "optimal_n_min" in lean.columns


# ---------------------------------------------------------------------------
# The per-cluster table
# ---------------------------------------------------------------------------


def test_grid_cluster_table_clips_each_cell_at_its_own_n_min(grid, summary):
    table = pgp.load_grid_cluster_table(grid, n_min=summary)

    assert len(table) == summary["n_real_clusters"].sum()
    assert (table["size_core"] >= table["n_min"]).all()

    for _, row in summary.iterrows():
        cell = table.loc[
            (table["order"] == row["order"]) & np.isclose(table["dmax"], row["dmax"])
        ]
        assert len(cell) == row["n_real_clusters"]
        assert (cell["n_min"] == row["optimal_n_min"]).all()


def test_grid_cluster_table_measures_neighbours_within_the_cell(grid, summary):
    """A cluster's neighbours are different at every (order, dmax); measuring
    them across the whole table would compare a cluster against copies of
    itself from other cells."""
    table = pgp.load_grid_cluster_table(grid, n_min=summary)

    for (order, dmax), cell in table.groupby(["order", "dmax"]):
        if len(cell) < 2:
            continue
        distance, ratio = brute_force_separation(cell)
        assert np.allclose(cell["closest_cluster_distance"], distance)
        assert np.allclose(cell["closest_cluster_distance_norm"], ratio)


def test_n_min_spec_accepts_a_number_a_mapping_or_a_summary(grid, summary):
    everything = pgp.load_grid_cluster_table(grid, n_min=None)

    threshold = int(everything["size_core"].median()) + 1
    clipped = pgp.load_grid_cluster_table(grid, n_min=threshold)
    assert (clipped["size_core"] >= threshold).all()
    assert 0 < len(clipped) < len(everything)

    mapping = {(e["order"], e["dmax"]): 4 for e in grid}
    by_mapping = pgp.load_grid_cluster_table(grid, n_min=mapping)
    assert (by_mapping["n_min"] == 4).all()

    lookup = _n_min_lookup(summary)
    for _, row in summary.iterrows():
        assert lookup(row["order"], row["dmax"]) == row["optimal_n_min"]

    with pytest.raises(ValueError, match="not a column"):
        _n_min_lookup(summary, column="nope")


def test_n_min_spec_accepts_optimise_n_min_output_for_a_1d_sweep(grid):
    """`optimise_n_min` returns a plain list, one per manifest entry - the
    1D notebook passes it straight in rather than building a mapping."""
    order_axis = pgp.select_sweep(grid, dmax=0.4)
    n_mins = [2, 30, 5]

    table = pgp.load_grid_cluster_table(order_axis, n_min=n_mins)
    for entry, n_min in zip(order_axis, n_mins):
        cell = table.loc[table["order"] == entry["order"]]
        assert (cell["n_min"] == n_min).all()
        assert (cell["size_core"] >= n_min).all()

    with pytest.raises(ValueError, match="one per manifest entry"):
        pgp.load_grid_cluster_table(order_axis, n_min=n_mins[:2])


# ---------------------------------------------------------------------------
# Following clusters between cells
# ---------------------------------------------------------------------------


def test_raising_the_order_is_reported_as_a_split(grid):
    """The tip has a dumbbell in it, so somewhere on the order axis a cluster
    must come apart - and it must be counted as a split rather than as one
    cluster dissolving and two appearing."""
    splits = pgp.track_cluster_splits(grid, axis="order", n_min=5)

    assert set(splits["axis"]) == {"order"}
    assert (splits["from_value"] < splits["to_value"]).all()
    assert splits["n_split"].sum() >= 1

    a_split = splits.loc[splits["n_split"] > 0].iloc[0]
    assert a_split["n_new_from_splits"] >= a_split["n_split"]
    assert len(a_split["split_cluster_ids"]) == a_split["n_split"]

    # the ids are the earlier cell's, and are real clusters of it
    parent_cell = next(
        e for e in grid
        if e["order"] == a_split["from_order"]
        and np.isclose(e["dmax"], a_split["from_dmax"])
    )
    parent_stats = pd.read_parquet(parent_cell["real_cluster_stats_path"])
    assert set(a_split["split_cluster_ids"]) <= set(parent_stats["cluster_id"])


def test_raising_dmax_merges_and_never_splits(grid):
    """DBSCAN is monotone in dmax: more links, so clusters can merge but not
    come apart. A split along that axis would mean the sweep had gone wrong."""
    splits = pgp.track_cluster_splits(grid, axis="dmax", n_min=1)

    assert (splits["n_split"] == 0).all()
    assert splits["n_merged"].sum() >= 1


def test_split_counts_add_up(grid):
    """Every cluster on each side is accounted for exactly once: it either
    survives (possibly splitting or merging) or it does not exist there."""
    for axis in ("order", "dmax"):
        splits = pgp.track_cluster_splits(grid, axis=axis, n_min=3)

        for _, row in splits.iterrows():
            # clusters of the earlier cell = those with children + dissolved
            children = row["n_clusters_from"] - row["n_dissolved"]
            # each of those maps to at least one later cluster, and the extras
            # are exactly what the splits created
            assert (
                row["n_clusters_to"] - row["n_created"] + row["n_lost_to_merges"]
                == children + row["n_new_from_splits"]
            )


def test_tracking_needs_the_per_atom_files(grid):
    stripped = [{k: v for k, v in e.items() if k != "real_clustered_atoms_path"}
                for e in grid]

    with pytest.raises(ValueError, match="save_real_atoms"):
        pgp.track_cluster_splits(stripped)

    with pytest.raises(ValueError, match='axis must be'):
        pgp.track_cluster_splits(grid, axis="n_min")


def test_clipping_keeps_the_original_cluster_ids():
    """Renumbering after the clip would break the round trip from a split back
    to the cluster in that cell's stats table."""
    labels = np.array([0, 0, 0, 1, 2, 2, 2, 2, -1])

    clipped = _clip_small_clusters(labels, n_min=3)

    assert set(clipped) == {0, 2, -1}
    assert (clipped[labels == 1] == -1).all()
    assert np.array_equal(_clip_small_clusters(labels, 1), labels)


# ---------------------------------------------------------------------------
# Alpha shapes
# ---------------------------------------------------------------------------


def sphere_of_atoms(n=400, radius=1.0, seed=0):
    rng = np.random.default_rng(seed)
    points = rng.normal(size=(n, 3))
    points /= np.linalg.norm(points, axis=1)[:, None]
    return points * radius * rng.uniform(0.85, 1.0, (n, 1))


def test_alpha_shape_of_a_ball_is_closed_and_about_the_right_size():
    atoms = pd.DataFrame(sphere_of_atoms(), columns=["x", "y", "z"])
    atoms["cluster_id"] = 0

    meshes = pgp.cluster_alpha_shapes(atoms, alpha=1.0)

    assert list(meshes) == [0]
    mesh = meshes[0]
    assert mesh["method"] == "alpha"
    assert mesh["watertight"]
    assert mesh["faces"].max() < len(mesh["vertices"])
    # a unit ball, give or take the surface the atoms actually trace out
    assert 2.5 < mesh["volume"] < 4.2


def test_a_degenerate_alpha_shape_falls_back_to_the_hull():
    """alphashape answers alpha=0 with a 2D polygon and a too-large alpha with
    a shredded surface; neither is a cluster shell, so the hull is used."""
    points = sphere_of_atoms(n=60)

    hull = _alpha_shape_mesh(points, alpha=0)
    assert hull["method"] == "hull"
    assert len(hull["faces"]) >= 4
    assert hull["volume"] > 0

    shredded = _alpha_shape_mesh(points, alpha=500)
    assert shredded["method"] == "hull"


def test_clusters_too_small_to_have_a_surface_are_skipped_not_faked():
    atoms = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0], "y": [0.0, 0.0, 0.0], "z": [0.0, 0.0, 0.0],
         "cluster_id": [0, 0, 0]}
    )

    with pytest.warns(UserWarning, match="no usable surface"):
        meshes = pgp.cluster_alpha_shapes(atoms, alpha=1.0)

    assert meshes == {}


# ---------------------------------------------------------------------------
# The plots
# ---------------------------------------------------------------------------


def test_tooltips_are_curated_by_default(grid, summary):
    table = pgp.load_grid_cluster_table(grid, n_min=summary)

    curated = list(pgp.hover_columns(table))
    assert curated[0] == "cluster_id"
    assert len(curated) < len(table.columns)
    assert "bbox_volume" not in curated

    assert len(pgp.hover_columns(table, "all")) == len(table.columns)
    assert list(pgp.hover_columns(table, ["size", "cluster_id"])) == [
        "size", "cluster_id"
    ]

    with pytest.raises(ValueError, match='"all"'):
        pgp.hover_columns(table, "everything")


def test_heatmap_hover_is_curated_and_takes_a_midpoint(summary):
    figure = pgp.plot_grid_heatmap(
        summary, value="min_cluster_distance_norm", colorscale="RdBu", zmid=1
    )

    assert figure.layout.coloraxis.cmid == 1
    hover = figure.data[0].customdata[0][0]
    assert hover.startswith("min_cluster_distance_norm")
    assert "n_real_clusters" in hover
    assert "run_name" not in hover

    # a column the curated list has never heard of stays out of the tooltip
    extended = summary.assign(scratch_column=1.0)
    curated = pgp.plot_grid_heatmap(extended).data[0].customdata[0][0]
    everything = pgp.plot_grid_heatmap(extended, hover="all").data[0].customdata[0][0]

    assert "scratch_column" not in curated
    assert "scratch_column" in everything


def test_distribution_grid_draws_a_curve_and_its_points_per_cell(grid, summary):
    table = pgp.load_grid_cluster_table(grid, n_min=summary)

    figure = pgp.plot_distribution_grid(
        table, value="size_core", row="order", color="dmax", log_x=True
    )

    curves = [t for t in figure.data if t.fill == "tozeroy"]
    points = [t for t in figure.data if t.mode == "markers"]
    assert len(curves) >= 1
    assert len(points) >= 1

    # one legend entry per colour value, not per (row, colour) pair
    legend = [t.name for t in figure.data if t.showlegend]
    assert legend == [f"dmax = {v}" for v in sorted(table["dmax"].unique())]

    # every point carries its cluster_id, which is the round trip to 3D
    assert points[0].hovertemplate.startswith("<b>Cluster %{customdata[0]}</b>")

    with pytest.raises(ValueError, match="not a column"):
        pgp.plot_distribution_grid(table, value="nope")


def test_split_tracking_plot_ticks_only_the_swept_values(grid):
    splits = pgp.track_cluster_splits(grid, axis="order", n_min=3)

    figure = pgp.plot_split_tracking(splits, value="n_new_from_splits")

    assert list(figure.layout.xaxis.tickvals) == sorted(splits["to_value"].unique())
    assert len(figure.data) == splits["dmax"].nunique()

    with pytest.raises(ValueError, match="not a column"):
        pgp.plot_split_tracking(splits, value="nope")


def test_alpha_shape_plot_colours_edge_clusters_apart():
    atoms = pd.concat(
        [
            pd.DataFrame(sphere_of_atoms(seed=1), columns=["x", "y", "z"]).assign(
                cluster_id=0
            ),
            pd.DataFrame(
                sphere_of_atoms(seed=2) + [5, 0, 0], columns=["x", "y", "z"]
            ).assign(cluster_id=1),
        ],
        ignore_index=True,
    )
    stats = pd.DataFrame(
        {"cluster_id": [0, 1], "size": [400, 400], "is_edge_cluster": [False, True]}
    )

    figure = pgp.plot_cluster_alpha_shapes(
        atoms, stats, alpha=1.0, opacity=0.6, progress=False
    )

    assert len(figure.data) == 2
    assert {t.color for t in figure.data} == {
        pgp.EDGE_CLUSTER_COLORS[False], pgp.EDGE_CLUSTER_COLORS[True]
    }
    assert all(t.opacity == 0.6 for t in figure.data)
    # one legend entry per group, and every mesh names its cluster
    assert sum(t.showlegend for t in figure.data) == 2
    assert all("Cluster" in t.hovertext for t in figure.data)


def test_cluster_evolution_view_has_one_legend_entry_per_swept_value(grid, summary):
    dmax_axis = pgp.select_sweep(grid, order=5)
    table = pgp.load_grid_cluster_table(grid, n_min=summary)
    biggest = (
        table.loc[(table["order"] == 5) & np.isclose(table["dmax"], 0.4)]
        .sort_values("size_core", ascending=False)["cluster_id"].iloc[0]
    )

    figure = pgp.compare_selected_clusters_across_param(
        dmax_axis, "dmax", chosen_value=0.4, selected_clusters=[biggest]
    )

    assert [t.name for t in figure.data] == ["0.4", "0.6"]
    # low values large, high values small and drawn on top, so both are visible
    assert figure.data[0].marker.size > figure.data[1].marker.size

    spread = pgp.compare_selected_clusters_across_param(
        dmax_axis, "dmax", chosen_value=0.4, selected_clusters=[biggest],
        spread=True, ion_colors={"Cu": "#b87333"},
    )
    # one 3D subplot per value, titled by it, all framing the same box
    assert [a.text for a in spread.layout.annotations] == ["dmax = 0.4", "dmax = 0.6"]
    assert {t.scene for t in spread.data} == {"scene", "scene2"}
    assert spread.layout.scene.xaxis.range == spread.layout.scene2.xaxis.range
    # coloured by ion, with one legend entry per ion across both subplots
    assert {t.marker.color for t in spread.data if t.name == "Cu"} == {"#b87333"}
    assert sorted(t.name for t in spread.data if t.showlegend) == sorted(
        {t.name for t in spread.data}
    )

    # the page that links the cameras: one script, not closed early by a
    # tooltip's "</extra>", naming both scenes
    html = pgp.viz._linked_scenes_html(spread)
    assert html.count("</script>") == 1
    assert '"scene2"' in html and "plotly_relayouting" in html

    with pytest.raises(ValueError, match="should match exactly one"):
        pgp.compare_selected_clusters_across_param(
            dmax_axis, "dmax", chosen_value=0.55, selected_clusters=[biggest]
        )


def test_cluster_evolution_view_draws_each_value_at_its_own_n_min(grid, summary):
    dmax_axis = pgp.select_sweep(grid, order=5)
    table = pgp.load_grid_cluster_table(grid, n_min=None)
    at_04 = table.loc[(table["order"] == 5) & np.isclose(table["dmax"], 0.4)]
    biggest = at_04.sort_values("size_core", ascending=False)["cluster_id"].iloc[0]
    everything = int(table["size_core"].max()) + 1

    # nothing survives at 0.6, so its subplot has only unclustered atoms
    spread = pgp.compare_selected_clusters_across_param(
        dmax_axis, "dmax", chosen_value=0.4, selected_clusters=[biggest],
        n_min=[1, everything], include_unclustered=True, spread=True,
    )
    assert [a.text for a in spread.layout.annotations] == [
        "dmax = 0.4, n_min 1", f"dmax = 0.6, n_min {everything}",
    ]
    assert {t.name for t in spread.data if t.scene == "scene2"} == {"unclustered"}
    assert {t.name for t in spread.data if t.scene == "scene"} - {"unclustered"}

    # a cluster clipped at its own value is drawn as unclustered - and said so
    with pytest.warns(UserWarning, match="below n_min"):
        pgp.compare_selected_clusters_across_param(
            dmax_axis, "dmax", chosen_value=0.4, selected_clusters=[biggest],
            n_min=everything, spread=True,
        )


def test_size_distribution_marks_each_values_n_min(grid):
    order_axis = pgp.select_sweep(grid, dmax=0.4)
    sizes = pgp.load_sweep_size_table(order_axis)

    figure = pgp.plot_size_distribution(
        sizes, x="order", y="size_core", n_min=[2, 4, 6]
    )
    tick = next(t for t in figure.data if t.name == "n_min")
    assert list(tick.x) == [3, 5, 8] and list(tick.y) == [2, 4, 6]

    with pytest.raises(ValueError, match="3 values of order"):
        pgp.plot_size_distribution(sizes, x="order", y="size_core", n_min=[2, 4])


def test_reference_line_spans_the_plot_whatever_the_x_values(grid, summary):
    """Placed in paper coordinates: an axis-domain ("x domain") line is read as
    data x = 0..1 by older renderers, which misses an Order axis entirely."""
    table = pgp.load_grid_cluster_table(grid, n_min=summary)

    figure = pgp.plot_size_distribution(
        table, x="order", y="closest_cluster_distance_norm", color=None,
        reference_line=(1, "equivalent spheres touch"),
    )
    (line,) = figure.layout.shapes
    assert (line.xref, line.x0, line.x1) == ("paper", 0, 1)
    assert (line.yref, line.y0, line.y1) == ("y", 1, 1)
    (label,) = figure.layout.annotations
    assert (label.xref, label.x, label.y) == ("paper", 1, 1)
    assert label.text == "equivalent spheres touch"

    unlabelled = pgp.plot_size_distribution(
        table, x="order", y="closest_cluster_distance_norm", reference_line=1
    )
    assert len(unlabelled.layout.shapes) == 1 and not unlabelled.layout.annotations


def make_composition_tables():
    stats = pd.DataFrame(
        {
            "cluster_id": [7, 3, 5, 9],
            "Mn": [10, 0, 4, 0], "Ni": [10, 2, 4, 0], "Si": [20, 2, 0, 0],
            "Cu": [0, 6, 2, 5], "Fe": [50, 10, 10, 10],
            "is_edge_cluster": [False, True, False, False],
        }
    )
    stats["cluster_size"] = stats[["Mn", "Ni", "Si", "Cu"]].sum(axis=1)
    matrix = pd.DataFrame(
        {"size": [10000], "Mn": [140], "Ni": [70], "Si": [50], "Cu": [5],
         "Fe": [9700], "P": [0]}
    )
    return stats, matrix


def test_ternary_places_each_cluster_by_its_three_elements_alone():
    stats, matrix = make_composition_tables()

    figure = pgp.plot_cluster_ternary(
        stats, ions=["Mn", "Ni", "Si"], unclustered_stats=matrix
    )
    traces = {t.name: t for t in figure.data}

    # cluster 7: 10 Mn, 10 Ni, 20 Si; cluster 5 has no Si; Cu and Fe ignored
    interior = traces["interior clusters"]
    assert np.allclose(interior.a, [25, 50]) and np.allclose(interior.c, [50, 0])
    assert traces["edge clusters"].marker.color == pgp.EDGE_CLUSTER_COLORS[True]
    # pooled: every cluster's atoms summed first (14 Mn, 16 Ni, 22 Si)
    assert np.allclose(
        [traces["all clusters"].a[0], traces["all clusters"].c[0]],
        [14 / 52 * 100, 22 / 52 * 100],
    )
    assert np.isclose(traces["matrix"].a[0], 140 / 260 * 100)
    # cluster 9 has none of the three and cannot be placed
    assert "1 clusters with none of the three left out" in figure.layout.title.text

    with pytest.warns(UserWarning, match=r"no \['Si'\] column"):
        pgp.plot_cluster_ternary(stats.drop(columns="Si"))


def test_composition_bars_are_each_100_percent_and_ordered_by_size():
    stats, matrix = make_composition_tables()

    figure = pgp.plot_cluster_composition_bars(
        stats, ions=["Cu", "Ni", "Si", "Mn", "P"], unclustered_stats=matrix
    )
    clusters = [t for t in figure.data if t.xaxis == "x"]
    references = [t for t in figure.data if t.xaxis == "x2"]

    # P has no atoms anywhere, so no legend entry for it
    assert [t.name for t in clusters] == ["Cu", "Ni", "Si", "Mn"]
    # smallest first, ties by id, each column labelled with its size
    assert list(clusters[0].x) == ["Cluster 9", "Cluster 3", "Cluster 5", "Cluster 7"]
    assert list(figure.layout.xaxis.ticktext) == ["5", "10", "10", "40"]
    assert np.allclose(np.sum([t.y for t in clusters], axis=0), 100)

    assert list(references[0].x) == ["all clusters", "matrix"]
    assert np.allclose(np.sum([t.y for t in references], axis=0), 100)
    mn = next(t for t in references if t.name == "Mn")
    assert np.allclose(mn.y, [14 / 65 * 100, 140 / 265 * 100])
