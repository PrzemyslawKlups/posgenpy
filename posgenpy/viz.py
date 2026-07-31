"""Interactive views of a cluster search: size/distance distributions where
every point names its own cluster, and 3D scatters of selected clusters.

Reads the persisted per-atom/per-cluster parquet files directly by cluster_id
(no re-reading/re-labelling pos+range files from scratch, unlike the old
file-based show_clusters_in_pos).

Every plot here carries a full hover layer, so a point in a distribution can be
read off and looked straight back up: hovering a cluster in
`plot_size_distribution` gives its cluster_id, which goes into
`show_clusters_in_pos` unchanged. That round trip is the one the old posgen
notebooks had (`hover_data=df.columns`) and is what makes an outlier in a size
distribution findable in 3D.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

# The old notebooks' real-vs-random pair, kept so these graphs read the same:
# a dark blue and a mid orange, which stay distinguishable under every common
# form of colour-vision deficiency and in greyscale print.
RUN_TYPE_COLORS = {"real": "#063970", "random": "#E28743"}

# columns that describe the run rather than the cluster - shown first in a
# tooltip, because they are what you need to reproduce or look up the point
_LEADING_HOVER_COLUMNS = ["cluster_id", "run", "run_type", "order", "dmax"]


def hover_columns(df, columns=None, skip=()):
    """Plotly `hover_data` spec showing every column, cluster_id first.

    Floats are formatted to 4 significant figures so a tooltip of centroid
    coordinates and shape metrics stays readable.
    """

    if columns is None:
        leading = [c for c in _LEADING_HOVER_COLUMNS if c in df.columns]
        columns = leading + [c for c in df.columns if c not in leading]

    columns = [c for c in columns if c in df.columns and c not in set(skip)]

    return {
        column: (":.4g" if pd.api.types.is_float_dtype(df[column]) else True)
        for column in columns
    }


def _with_cluster_label(df):
    """Copy of `df` with a `cluster` column used as the tooltip's heading."""
    out = df.copy()
    if "cluster_id" in out.columns:
        out["cluster"] = "Cluster " + out["cluster_id"].astype(str)
    return out


def _limit_points(stats, x, max_points, seed=0):
    """Cap how many individual points a distribution plot draws.

    A whole grid's size table is one row per cluster per control run - hundreds
    of thousands of rows - and `points="all"` draws every one, which is enough
    to lock up a browser tab. Points are thinned by sampling within each x
    group, which leaves the shape and relative weight of each violin intact.

    Real clusters get first call on the budget, since those are the ones worth
    hovering and looking up, but never more than half of it while there are
    randomised ones to show - a real-vs-random plot with no random points in it
    would be worse than a thinned one.

    Pass `max_points=None` to plot everything anyway.
    """

    if max_points is None or len(stats) <= max_points:
        return stats

    is_real = (
        stats["run_type"] == "real"
        if "run_type" in stats.columns
        else pd.Series(True, index=stats.index)
    )
    real, rest = stats.loc[is_real], stats.loc[~is_real]
    real_budget = (
        max_points
        if rest.empty
        else max(max_points // 2, max_points - len(rest))
    )

    def thin(frame, budget):
        if len(frame) <= budget:
            return frame
        if x and x in frame.columns and frame[x].nunique() > 1:
            # proportional within each x group, so the violins keep their
            # relative weight as well as their shape
            return frame.groupby(x, observed=True, group_keys=False).sample(
                frac=budget / len(frame), random_state=seed
            )
        return frame.sample(budget, random_state=seed)

    real = thin(real, real_budget)
    rest = thin(rest, max(0, max_points - len(real)))

    kept = pd.concat([real, rest]) if len(rest) else real
    warnings.warn(
        f"plotting {len(kept):,} of {len(stats):,} clusters "
        f"({len(real):,} real, {len(rest):,} randomised): points are thinned "
        f"per {x or 'group'} above max_points={max_points:,} so the browser can "
        "draw them. Pass max_points=None to plot every one."
    )
    return kept


def plot_size_distribution(
    stats,
    x: str = "dmax",
    y: str = "size_core",
    color: str = "run_type",
    hover: list = None,
    title: str = None,
    color_map: dict = None,
    figsize: tuple = (900, 520),
    point_size: int = 4,
    max_points: int = 20000,
):
    """Violin + every individual cluster as a hoverable point.

    Parameters
    ----------
    stats : pandas.DataFrame
        Per-cluster rows - a single run's `analyse_clusters` output, or the
        tidy table `fast_sweep.load_sweep_size_table` builds across a whole
        sweep.

    x : str
        What to split the distributions by: "dmax", "order", or any other
        column. Pass None for a single distribution.

    y : str
        The distributed quantity. "size_core" is posgen's `Cluster Size(core)`
        (core+border atoms) - the same number `optimise_n_min` thresholds;
        "size" counts every atom in the cluster, and
        "closest_cluster_distance" (after `find_knn_cluster_distance`) gives
        the kNN-distance distribution.

    color : str
        Column to colour by; "run_type" separates real from randomised runs,
        as in the old notebooks.

    hover : list[str], optional
        Columns to put in the tooltip. Defaults to all of them, cluster_id
        first - so a point can be identified and passed to
        `show_clusters_in_pos`.

    max_points : int
        Cap on individually drawn points, since a whole grid's worth of
        randomised runs is enough to lock up a browser tab. Real clusters are
        kept in full and the randomised ones thinned per x group (see
        `_limit_points`); `None` plots every one.

    Returns
    -------
    plotly.graph_objects.Figure
    """

    stats = _with_cluster_label(_limit_points(stats, x, max_points))

    if color is not None and color not in stats.columns:
        color = None

    fig = px.violin(
        stats,
        x=x,
        y=y,
        color=color,
        box=True,
        points="all",
        hover_name="cluster" if "cluster" in stats.columns else None,
        hover_data=hover_columns(stats, hover, skip=("cluster",)),
        color_discrete_map=color_map or RUN_TYPE_COLORS,
    )

    fig.update_traces(marker_size=point_size, marker_opacity=0.7, meanline_visible=True)
    fig.update_layout(
        title=title or f"{y} distribution" + (f" vs {x}" if x else ""),
        xaxis_title=x,
        yaxis_title=y,
        violingap=0,
        violinmode="overlay",
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
    )
    fig.update_yaxes(rangemode="nonnegative")

    return fig


def plot_grid_heatmap(
    summary,
    value: str = "optimal_n_min",
    title: str = None,
    figsize: tuple = (760, 460),
    colorscale: str = "Viridis",
    text_format: str = None,
    mask_unreached: bool = True,
):
    """Heatmap of one `fast_sweep.summarise_grid` column over order x dmax.

    The two 1D sweep graphs answer "which dmax at this order" and "which order
    at this dmax". This shows the whole grid at once, which is what choosing the
    pair actually needs - a plateau where `optimal_n_min` and
    `n_real_clusters` are both stable is a parameter region that is not
    balanced on a knife edge.

    Every cell's tooltip carries the rest of its row, so the n_min, the real and
    random cluster counts and the ratio can all be read off one cell without
    replotting.

    Parameters
    ----------
    summary : pandas.DataFrame
        `summarise_grid` output (or anything with order/dmax and `value`).
    value : str
        Column to colour by, e.g. "optimal_n_min", "n_real_clusters",
        "true_positive_ratio".
    text_format : str, optional
        Format for the number printed in each cell, e.g. ".2f". Defaults to
        integers for integer columns and 2 decimals otherwise.
    mask_unreached : bool
        Leave cells that never reached the confidence threshold blank. They are
        reported at the largest n_min tried, which is not a result but a
        did-not-finish - plotted, it is an extreme value that flattens the
        colour scale over every cell that did work. Their numbers are still in
        the tooltip.
    """

    if value not in summary.columns:
        raise ValueError(
            f"{value!r} is not a column of this summary ({list(summary.columns)})"
        )

    grid = summary.pivot(index="order", columns="dmax", values=value)

    if mask_unreached and "n_min_reached" in summary.columns:
        reached = summary.pivot(index="order", columns="dmax", values="n_min_reached")
        grid = grid.astype(float).where(reached.astype(bool))

    # every other column of the cell, laid over the same grid, so one hover
    # answers "why is this cell that colour" without replotting. Masked cells
    # keep their numbers here - that is how you see what they did reach.
    hover = pd.DataFrame("", index=grid.index, columns=grid.columns)
    for column in summary.columns:
        if column in ("order", "dmax", "run_name"):
            continue
        plane = summary.pivot(index="order", columns="dmax", values=column)
        formatted = plane.map(lambda v: f"{v:.4g}" if isinstance(v, float) else str(v))
        hover = hover + f"{column}: " + formatted + "<br>"

    if text_format is None:
        integral = np.isclose(
            np.nan_to_num(grid.to_numpy(dtype=float)) % 1, 0
        ).all()
        text_format = "d" if integral else ".2f"

    fig = px.imshow(
        grid,
        labels={"x": "dmax", "y": "order", "color": value},
        color_continuous_scale=colorscale,
        text_auto=text_format,
        aspect="auto",
    )
    fig.update_traces(
        customdata=hover.to_numpy(),
        hovertemplate="order %{y}, dmax %{x}<br><br>%{customdata}<extra></extra>",
    )
    fig.update_layout(
        title=title or f"{value} across the (order x dmax) grid",
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
    )
    fig.update_xaxes(type="category")
    fig.update_yaxes(type="category")

    return fig


def show_clusters_in_pos(
    clustered_atoms,
    clusters_to_plot: list,
    ion_colors: dict = None,
    hover: list = None,
    marker_size: int = 2,
):
    """
    3D scatter of one or more selected clusters, atom by atom.

    clustered_atoms: a path to a `*_real_clustered_atoms.parquet` file (as
    returned in a sweep/run_relabel_sweep manifest entry), or an
    already-loaded DataFrame with x, y, z, ion, cluster_id columns.

    Each atom's tooltip carries its cluster_id and ion (and ionIdx, so it can
    be traced back to a row of the original pos file), which means a cluster
    picked out of `plot_size_distribution` can be confirmed here, and a cluster
    spotted here can be found in the stats table.
    """

    if isinstance(clustered_atoms, (str, Path)):
        df = pd.read_parquet(clustered_atoms)
    else:
        df = clustered_atoms

    selected = df.loc[df["cluster_id"].isin(clusters_to_plot)]
    labelled = _with_cluster_label(selected)

    fig = px.scatter_3d(
        labelled, x="x", y="y", z="z", color="ion",
        color_discrete_map=ion_colors or {},
        hover_name="cluster" if "cluster" in labelled.columns else None,
        hover_data=hover_columns(labelled, hover, skip=("cluster",)),
    )
    fig.update_traces(marker_size=marker_size)
    fig.update_layout(title=f"Cluster ID(s): {clusters_to_plot}")
    fig.show()

    return selected


def compare_selected_clusters_across_param(
    sweep_manifest: list,
    swept_param_name: str,
    chosen_value,
    selected_clusters: list,
    hover: list = None,
):
    """
    3D scatter comparing the spatial neighbourhood of `selected_clusters`
    (chosen from the swept value `chosen_value`) across every swept value in
    `sweep_manifest`, by filtering every value's clustered atoms to the
    bounding box of the selected clusters.

    Cluster ids are per-run, so the same physical cluster generally has a
    different id at each swept value - the tooltip carries each atom's id at
    its own swept value, which is what tells you how a cluster split or merged.
    """

    chosen_entries = [e for e in sweep_manifest if e["swept_value"] == chosen_value]
    if len(chosen_entries) != 1:
        print("Found 0 or 2< entries matching chosen_value. Can't continue")
        return None

    chosen_df = pd.read_parquet(chosen_entries[0]["real_clustered_atoms_path"])
    chosen_df = chosen_df.loc[chosen_df["cluster_id"].isin(selected_clusters)]

    xmin, xmax = chosen_df["x"].min(), chosen_df["x"].max()
    ymin, ymax = chosen_df["y"].min(), chosen_df["y"].max()
    zmin, zmax = chosen_df["z"].min(), chosen_df["z"].max()

    frames = []
    for entry in sweep_manifest:
        df = pd.read_parquet(entry["real_clustered_atoms_path"])
        mask = (
            (df["x"] >= xmin) & (df["x"] <= xmax)
            & (df["y"] >= ymin) & (df["y"] <= ymax)
            & (df["z"] >= zmin) & (df["z"] <= zmax)
        )
        df = df.loc[mask].copy()
        df[swept_param_name] = str(entry["swept_value"])
        df["cluster_label"] = df["cluster_id"].astype(str)
        frames.append(df)

    df_to_plot = pd.concat(frames, ignore_index=True)

    fig = px.scatter_3d(
        df_to_plot, x="x", y="y", z="z",
        color=swept_param_name, symbol="cluster_label",
        hover_data=hover_columns(df_to_plot, hover, skip=("cluster_label",)),
    )
    fig.update_traces(marker_size=3, opacity=0.8)

    return fig
