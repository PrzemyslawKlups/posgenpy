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
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .corrections import cluster_alpha_shapes

# The old notebooks' real-vs-random pair, kept so these graphs read the same:
# a dark blue and a mid orange, which stay distinguishable under every common
# form of colour-vision deficiency and in greyscale print.
RUN_TYPE_COLORS = {"real": "#063970", "random": "#E28743"}

# the same two hues, for the other binary split in this workflow
EDGE_CLUSTER_COLORS = {False: "#063970", True: "#E28743"}

# columns that describe the run rather than the cluster - shown first in a
# tooltip, because they are what you need to reproduce or look up the point
_LEADING_HOVER_COLUMNS = ["cluster_id", "run", "run_type", "order", "dmax"]

# What a tooltip shows when nothing else is asked for. A cluster stats table
# has 20+ columns and a tooltip of all of them is a wall of numbers nobody
# reads; these are the ones a point is actually identified and judged by -
# cluster_id first, because that is what gets typed into show_clusters_in_pos.
# Pass columns="all" for the lot, or an explicit list.
DEFAULT_HOVER_COLUMNS = [
    "cluster_id", "run", "run_type", "order", "dmax", "n_min",
    "size", "size_core", "cluster_size",
    "closest_cluster_distance", "closest_cluster_id",
    "closest_cluster_distance_norm",
    "equivalent_spherical_radius", "is_edge_cluster",
    "ion", "ionIdx",
]

# ...and the same idea for a heatmap cell, whose "row" is one grid cell
DEFAULT_GRID_HOVER_COLUMNS = [
    "optimal_n_min", "n_real_clusters", "n_random_clusters",
    "true_positive_ratio", "n_min_reached",
    "n_clusters_all", "max_cluster_size",
    "min_cluster_distance", "min_cluster_distance_norm",
    "n_clusters_from", "n_clusters_to", "n_split", "n_new_from_splits",
    "n_merged", "n_lost_to_merges", "n_dissolved", "n_created",
]


def hover_columns(df, columns=None, skip=()):
    """Plotly `hover_data` spec, cluster_id first.

    Parameters
    ----------
    columns : list, "all", or None
        None (the default) shows `DEFAULT_HOVER_COLUMNS` - the ~10 columns a
        cluster is identified and judged by, rather than every column of a
        stats table. `"all"` is the old behaviour: every column, cluster_id
        first. A list is taken as given, in the order given.

    Floats are formatted to 4 significant figures so a tooltip of centroid
    coordinates and shape metrics stays readable.
    """

    leading = [c for c in _LEADING_HOVER_COLUMNS if c in df.columns]

    if columns is None:
        columns = [c for c in DEFAULT_HOVER_COLUMNS if c in df.columns]
        # a frame with none of them (say, a bare x/y/z table) is better off
        # with everything than with an empty tooltip
        if not columns:
            columns = leading + [c for c in df.columns if c not in leading]
    elif isinstance(columns, str):
        if columns != "all":
            raise ValueError('columns must be a list, "all", or None')
        columns = leading + [c for c in df.columns if c not in leading]

    columns = [c for c in columns if c in df.columns and c not in set(skip)]

    return {
        column: (":.4g" if pd.api.types.is_float_dtype(df[column]) else True)
        for column in columns
    }


def _hover_layer(df, columns=None, skip=()):
    """(customdata, hovertemplate) for a hand-built trace.

    `hover_columns` covers the plotly-express path, where px builds the
    template itself. A trace assembled by hand needs the same tooltip without
    px inventing entries for the axes it was not given (an `index=` line for a
    scatter whose y is only a baseline offset).
    """

    spec = hover_columns(df, columns, skip=skip)
    columns = list(spec)

    lines = [
        f"{column}: %{{customdata[{i}]{'' if spec[column] is True else spec[column]}}}"
        for i, column in enumerate(columns)
    ]
    if columns and columns[0] == "cluster_id":
        lines[0] = "<b>Cluster %{customdata[0]}</b>"

    return df[columns].to_numpy(), "<br>".join(lines) + "<extra></extra>"


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
    hover: list = None,
    zmid: float = None,
    colorbar_title: str = None,
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
    colorscale : str
        Any plotly colour scale. Use a **single-hue sequential** one
        ("Blues", "Purples", "Greens", "Oranges", "Teal") for a magnitude, and
        a different hue per quantity so two heatmaps of the same grid are not
        mistaken for each other. Use a **diverging** one ("RdBu") with `zmid`
        for a quantity that has a meaningful midpoint.
    zmid : float, optional
        Value the colour scale is centred on. For
        `min_cluster_distance_norm`, `zmid=1` puts the colour break exactly at
        "the two clusters' equivalent spheres touch", so overlapping cells
        read as one pole and separated ones as the other.
    text_format : str, optional
        Format for the number printed in each cell, e.g. ".2f". Defaults to
        integers for integer columns and 2 decimals otherwise.
    mask_unreached : bool
        Leave cells that never reached the confidence threshold blank. They are
        reported at the largest n_min tried, which is not a result but a
        did-not-finish - plotted, it is an extreme value that flattens the
        colour scale over every cell that did work. Their numbers are still in
        the tooltip.
    hover : list, optional
        Columns in the tooltip. Defaults to the plotted value plus whichever
        of `DEFAULT_GRID_HOVER_COLUMNS` the frame has - the n_min, the counts
        and the ratio - rather than every column it happens to carry. Pass
        `"all"` for everything, or a list of your own.
    """

    if value not in summary.columns:
        raise ValueError(
            f"{value!r} is not a column of this summary ({list(summary.columns)})"
        )

    grid = summary.pivot(index="order", columns="dmax", values=value)

    if mask_unreached and "n_min_reached" in summary.columns:
        reached = summary.pivot(index="order", columns="dmax", values="n_min_reached")
        grid = grid.astype(float).where(reached.astype(bool))

    # the cell's other numbers, laid over the same grid, so one hover answers
    # "why is this cell that colour" without replotting. Masked cells keep
    # their numbers here - that is how you see what they did reach.
    if hover is None:
        hover = [value] + [c for c in DEFAULT_GRID_HOVER_COLUMNS if c != value]
    elif isinstance(hover, str):
        if hover != "all":
            raise ValueError('hover must be a list, "all", or None')
        hover = [value] + [c for c in summary.columns if c != value]

    hover_text = pd.DataFrame("", index=grid.index, columns=grid.columns)
    for column in hover:
        if column in ("order", "dmax", "run_name") or column not in summary.columns:
            continue
        plane = summary.pivot(index="order", columns="dmax", values=column)
        formatted = plane.map(lambda v: f"{v:.4g}" if isinstance(v, float) else str(v))
        hover_text = hover_text + f"{column}: " + formatted + "<br>"

    if text_format is None:
        integral = np.isclose(
            np.nan_to_num(grid.to_numpy(dtype=float)) % 1, 0
        ).all()
        text_format = "d" if integral else ".2f"

    fig = px.imshow(
        grid,
        labels={"x": "dmax", "y": "order", "color": colorbar_title or value},
        color_continuous_scale=colorscale,
        color_continuous_midpoint=zmid,
        text_auto=text_format,
        aspect="auto",
    )
    fig.update_traces(
        customdata=hover_text.to_numpy(),
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


def _colour_steps(values, colorscale, low=0.15, high=0.92):
    """One colour per value, stepped along a sequential scale in value order.

    The swept parameters are ordered magnitudes, not names, so their colours
    are a ramp: light at the low end, dark at the high end. A reader can then
    tell which curve is which dmax without going back to the legend.
    """
    if len(values) == 1:
        return {values[0]: px.colors.sample_colorscale(colorscale, [high])[0]}

    positions = np.linspace(low, high, len(values))
    colours = px.colors.sample_colorscale(colorscale, positions)
    return dict(zip(values, colours))


def _to_rgba(colour, alpha):
    """A plotly colour string as rgba with the given alpha."""
    r, g, b = px.colors.unlabel_rgb(px.colors.convert_colors_to_same_type(colour)[0][0])
    return f"rgba({r:.0f},{g:.0f},{b:.0f},{alpha})"


def plot_distribution_grid(
    table,
    value: str = "size_core",
    row: str = "order",
    color: str = "dmax",
    points: bool = True,
    hover: list = None,
    title: str = None,
    colorscale: str = "Viridis",
    normalise: str = "curve",
    bandwidth=None,
    log_x: bool = False,
    n_grid: int = 256,
    figsize: tuple = None,
    max_points: int = 4000,
    point_size: int = 5,
    fill_opacity: float = 0.35,
):
    """Both grid axes in one distribution: a ridgeline of half violins.

    A violin per cell of an (order x dmax) grid is 30-odd violins nobody can
    read side by side. This instead gives **one subplot row per `row` value**
    and, inside each, **one filled density curve per `color` value** - which
    is a half violin laid on its side, so the whole grid is read as "how does
    this distribution move as dmax rises (the coloured curves), and does that
    change with order (the rows)".

    Every individual cluster is still drawn, as a point under its own curve,
    carrying its `cluster_id` in the tooltip - so an outlier in the tail can
    be identified here and passed straight to `show_clusters_in_pos`.

    Parameters
    ----------
    table : pandas.DataFrame
        One row per cluster, with `row`, `color` and `value` columns -
        `fast_sweep.load_grid_cluster_table` output, or
        `load_sweep_size_table`.
    value : str
        The distributed quantity: "size_core" for cluster size,
        "closest_cluster_distance" for the nearest-neighbour distance,
        "closest_cluster_distance_norm" for that distance in units of the two
        clusters' radii.
    normalise : {"curve", "row", "none"}
        What the curve heights mean. "curve" gives every curve the same peak,
        comparing shapes (what a violin plot does). "row" weights each curve
        by how many clusters it holds, so a parameter value that found twice
        as many draws twice the area. "none" leaves each curve as a plain
        probability density, so a narrow distribution stands taller than a
        wide one. In every mode each row is then scaled so its tallest curve
        fills the row - the y axis carries no numbers, only shape.
    log_x : bool
        Estimate the density in log10(value) and plot on a log axis. Cluster
        sizes span two decades and are far more readable this way; distances
        usually are not.
    bandwidth : float, optional
        Passed to `scipy.stats.gaussian_kde` as its `bw_method` - raise it to
        smooth a noisy curve, lower it to keep detail.

    Returns
    -------
    plotly.graph_objects.Figure
    """

    from scipy.stats import gaussian_kde

    if normalise not in ("curve", "row", "none"):
        raise ValueError('normalise must be "curve", "row" or "none"')

    for column in (value, row, color):
        if column not in table.columns:
            raise ValueError(
                f"{column!r} is not a column of this table ({list(table.columns)})"
            )

    table = table.loc[np.isfinite(table[value].to_numpy(dtype=float))]
    if log_x:
        table = table.loc[table[value] > 0]
    if table.empty:
        raise ValueError(f"no finite {value!r} values to plot")

    row_values = sorted(table[row].unique())
    colour_values = sorted(table[color].unique())
    colours = _colour_steps(colour_values, colorscale)

    scale = np.log10 if log_x else (lambda v: v)
    samples = scale(table[value].to_numpy(dtype=float))
    x_grid = np.linspace(samples.min(), samples.max(), n_grid)
    x_plot = 10 ** x_grid if log_x else x_grid

    fig = make_subplots(
        rows=len(row_values), cols=1, shared_xaxes=True,
        # plotly caps this at 1/(rows-1), which a tall grid would otherwise hit
        vertical_spacing=min(0.06, 0.8 / max(len(row_values), 2)),
        subplot_titles=[f"{row} = {v}" for v in row_values],
    )

    drawn = _limit_points(table, color, max_points) if points else None

    for row_index, row_value in enumerate(row_values, start=1):
        in_row = table[table[row] == row_value]
        row_peak = 0.0
        curves = {}

        for colour_value in colour_values:
            values = scale(
                in_row.loc[in_row[color] == colour_value, value].to_numpy(dtype=float)
            )
            if len(values) < 2 or np.ptp(values) == 0:
                continue

            density = gaussian_kde(values, bw_method=bandwidth)(x_grid)
            if normalise == "curve":
                density = density / density.max()
            elif normalise == "row":
                density = density * len(values)

            curves[colour_value] = density
            row_peak = max(row_peak, density.max())

        row_peak = row_peak or 1.0

        for colour_value, density in curves.items():
            colour = colours[colour_value]
            fig.add_trace(
                go.Scatter(
                    x=x_plot, y=density / row_peak,
                    mode="lines", line=dict(color=colour, width=2),
                    fill="tozeroy", fillcolor=_to_rgba(colour, fill_opacity),
                    name=f"{color} = {colour_value}",
                    legendgroup=str(colour_value),
                    showlegend=row_index == 1,
                    hovertemplate=(
                        f"{row} {row_value}, {color} {colour_value}"
                        f"<br>{value} %{{x:.4g}}<extra></extra>"
                    ),
                ),
                row=row_index, col=1,
            )

        if not points:
            continue

        # every cluster under its own curve, on its own baseline per colour so
        # the rows of points do not pile up on each other - and close enough
        # together that they stay a strip under the curves, not a second plot
        in_row_points = drawn[drawn[row] == row_value]
        spacing = min(0.09, 0.6 / len(colour_values))

        for offset, colour_value in enumerate(colour_values, start=1):
            group = in_row_points[in_row_points[color] == colour_value]
            if group.empty:
                continue

            jitter = np.random.default_rng(offset).uniform(
                -0.2 * spacing, 0.2 * spacing, len(group)
            )
            customdata, hovertemplate = _hover_layer(group, hover)

            fig.add_trace(
                go.Scatter(
                    x=group[value], y=-spacing * offset + jitter,
                    mode="markers",
                    marker=dict(
                        color=colours[colour_value], size=point_size, opacity=0.55
                    ),
                    customdata=customdata,
                    hovertemplate=hovertemplate,
                    legendgroup=str(colour_value),
                    showlegend=False,
                ),
                row=row_index, col=1,
            )

    fig.update_layout(
        title=title or f"{value} across the ({row} x {color}) grid",
        template="simple_white",
        width=(figsize or (900, 0))[0],
        height=(figsize[1] if figsize else 130 * len(row_values) + 190),
        legend_title_text=color,
        hovermode="closest",
    )
    fig.update_yaxes(
        title_text="density", showticklabels=False, showgrid=False, zeroline=False
    )
    fig.update_xaxes(showgrid=False)
    if log_x:
        # decades only: plotly labels every minor tick by default, which on a
        # two-decade size axis is 20 numbers along the bottom
        fig.update_xaxes(type="log", dtick=1)
    fig.update_xaxes(title_text=value, row=len(row_values), col=1)

    return fig


def plot_split_tracking(
    splits,
    value: str = "n_new_from_splits",
    color: str = None,
    title: str = None,
    colorscale: str = "Viridis",
    figsize: tuple = (860, 460),
):
    """How many clusters split (or merged) at each step along a grid axis.

    Takes `fast_sweep.track_cluster_splits` output and plots one line per
    value of the parameter that was held fixed, against the one that moved.
    The step where the line takes off is the elbow: below it the parameter is
    only tidying edges, above it the clustering is being rearranged.

    `value` is any count in that table - "n_new_from_splits" (extra clusters
    created by splits), "n_merged", "n_dissolved", "n_created".
    """

    if splits.empty:
        raise ValueError("no transitions to plot")
    if value not in splits.columns:
        raise ValueError(
            f"{value!r} is not a column of this table ({list(splits.columns)})"
        )

    axis = splits["axis"].iloc[0]
    color = color or ("dmax" if axis == "order" else "order")

    colour_values = sorted(splits[color].unique())
    colours = _colour_steps(colour_values, colorscale)

    fig = px.line(
        splits.sort_values([color, "to_value"]),
        x="to_value", y=value, color=color, markers=True,
        color_discrete_map={v: colours[v] for v in colour_values},
        hover_data=[
            c for c in (
                "from_value", "to_value", "n_clusters_from", "n_clusters_to",
                "n_split", "n_new_from_splits", "n_merged", "n_dissolved",
                "n_created",
            ) if c in splits.columns
        ],
    )
    fig.update_traces(line_width=2, marker_size=8)
    fig.update_layout(
        title=title or f"{value} between consecutive {axis} values",
        xaxis_title=f"{axis} (arriving at)",
        yaxis_title=value,
        legend_title_text=color,
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
    )
    # only the values that were actually swept: an Order axis has no 4.5
    fig.update_xaxes(tickmode="array", tickvals=sorted(splits["to_value"].unique()))
    fig.update_yaxes(rangemode="nonnegative")

    return fig


def plot_cluster_alpha_shapes(
    clustered_atoms,
    cluster_stats=None,
    clusters_to_plot: list = None,
    alpha: float = 1.0,
    opacity: float = 0.6,
    colors: dict = None,
    color_by: str = "is_edge_cluster",
    show_atoms: bool = False,
    hover: list = None,
    title: str = None,
    figsize: tuple = (900, 800),
    progress: bool = True,
):
    """The clusters as surfaces: one alpha shape per cluster, see-through.

    A 3D scatter of atoms shows where the atoms are; this shows what the
    search decided their *shape* was. Every mesh is drawn at `opacity` (0.6 by
    default) so clusters behind other clusters stay visible, and each carries
    its cluster_id in the tooltip.

    Clusters flagged by `corrections.remove_edge_clusters` as touching the
    tip boundary are drawn in the second colour - they are the ones whose
    volume and composition are only partly measured, so seeing where they sit
    says how much of the tip's edge the analysis is leaning on.

    Parameters
    ----------
    clustered_atoms : DataFrame or path
        Per-atom table with x/y/z and cluster_id.
    cluster_stats : DataFrame, optional
        Per-cluster table, for the `color_by` flag and the tooltips.
    alpha : float
        Alpha-shape parameter in inverse nm - the rolled probe sphere has
        radius `1/alpha`. `1/dmax` is a good starting point; too large tears
        holes in the surface, and where that happens the convex hull is used
        instead (`corrections.cluster_alpha_shapes`).
    color_by : str
        Column of `cluster_stats` splitting the clusters into two colours.
        Defaults to the edge-cluster flag; None draws them all alike.
    show_atoms : bool
        Also scatter the atoms inside the shells.
    """

    if isinstance(clustered_atoms, (str, Path)):
        clustered_atoms = pd.read_parquet(clustered_atoms)

    if clusters_to_plot is None:
        clusters_to_plot = clustered_atoms.loc[
            clustered_atoms["cluster_id"] != -1, "cluster_id"
        ].unique()
    clusters_to_plot = list(clusters_to_plot)

    meshes = cluster_alpha_shapes(
        clustered_atoms, clusters_to_plot, alpha=alpha, progress=progress
    )

    colors = colors or EDGE_CLUSTER_COLORS
    flags, tooltips = {}, {}

    if cluster_stats is not None:
        indexed = cluster_stats.set_index("cluster_id")
        if color_by and color_by in indexed.columns:
            flags = indexed[color_by].to_dict()

        shown = list(hover_columns(cluster_stats, hover, skip=("cluster_id",)))
        for cid, row in indexed[shown].iterrows():
            tooltips[cid] = "<br>".join(
                f"{column}: {row[column]:.4g}" if isinstance(row[column], float)
                else f"{column}: {row[column]}"
                for column in shown
            )

    fig = go.Figure()
    labelled = set()

    for cid, mesh in meshes.items():
        flag = bool(flags.get(cid, False))
        name = (
            f"{color_by} = {flag}" if flags else "cluster"
        )
        vertices, faces = mesh["vertices"], mesh["faces"]

        fig.add_trace(
            go.Mesh3d(
                x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                color=colors.get(flag, colors[False]),
                opacity=opacity,
                flatshading=True,
                name=name,
                legendgroup=name,
                showlegend=name not in labelled,
                hoverinfo="text",
                hovertext=(
                    f"<b>Cluster {cid}</b><br>"
                    + (f"{tooltips[cid]}<br>" if cid in tooltips else "")
                    + f"surface: {mesh['method']}"
                    + (
                        f", volume {mesh['volume']:.3g} nm3"
                        if np.isfinite(mesh["volume"]) else ", open surface"
                    )
                ),
            )
        )
        labelled.add(name)

    if show_atoms:
        inside = clustered_atoms.loc[
            clustered_atoms["cluster_id"].isin(list(meshes))
        ]
        fig.add_trace(
            go.Scatter3d(
                x=inside["x"], y=inside["y"], z=inside["z"],
                mode="markers",
                marker=dict(size=1.5, color="#444444", opacity=0.5),
                name="atoms",
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        title=title or f"{len(meshes)} cluster alpha shapes (alpha={alpha})",
        scene=dict(
            xaxis_title="x [nm]", yaxis_title="y [nm]", zaxis_title="z [nm]",
            aspectmode="data",
        ),
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
        margin=dict(l=20, r=20, t=60, b=20),
    )

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
    padding: float = 0.5,
    include_unclustered: bool = False,
    spread: bool = False,
    colorscale: str = "Viridis",
    marker_size: int = 3,
    size_step: float = 1.6,
    title: str = None,
    figsize: tuple = (900, 750),
):
    """One cluster, at every value of a swept parameter, in one 3D view.

    Takes the neighbourhood of `selected_clusters` (as they are at
    `chosen_value`) and draws what every other swept value made of the same
    piece of the tip, **each value in its own colour**. A cluster splitting in
    two as the Order rises shows up as the later colours arriving where the
    earlier ones had a single blob; a merge is the same thing backwards.

    Every value draws the *same atoms* at the *same coordinates*, so drawn flat
    they would hide each other and only the last colour would be visible. The
    markers are therefore graded by size - low parameter values large, high
    values small and drawn on top - which leaves all of them readable: a big
    dot with a small dot inside it is an atom clustered at both values, a bare
    small dot is one that only the higher value picked up, and a bare big dot
    is one the higher value dropped.

    There is one legend entry per swept value, so clicking one hides that
    value - comparing two at a time is how this plot is read. `spread=True`
    lays the values out side by side along x instead of on top of each other,
    which is the clearer way to watch a cluster break in two.

    Cluster ids are per-run, so the same physical cluster generally has a
    different id at each swept value; each atom's tooltip carries its id at
    its own value, which is what says *what* it merged into.

    Parameters
    ----------
    padding : float
        Extra nm around the chosen clusters' bounding box, so a cluster that
        grows at a looser parameter is not cut off at the box edge.
    include_unclustered : bool
        Also draw core-species atoms the search left unclustered in that box.
    spread : bool
        Offset each swept value along x into its own copy, side by side.
    """

    chosen_entries = [
        e for e in sweep_manifest
        if np.isclose(float(e["swept_value"]), float(chosen_value))
    ]
    if len(chosen_entries) != 1:
        raise ValueError(
            f"{len(chosen_entries)} entries match {swept_param_name}="
            f"{chosen_value}; it should match exactly one of "
            f"{[e['swept_value'] for e in sweep_manifest]}"
        )

    chosen_df = pd.read_parquet(chosen_entries[0]["real_clustered_atoms_path"])
    chosen_df = chosen_df.loc[chosen_df["cluster_id"].isin(selected_clusters)]
    if chosen_df.empty:
        raise ValueError(
            f"no atoms in clusters {selected_clusters} at "
            f"{swept_param_name}={chosen_value}"
        )

    low = chosen_df[["x", "y", "z"]].min() - padding
    high = chosen_df[["x", "y", "z"]].max() + padding

    swept_values = sorted(e["swept_value"] for e in sweep_manifest)
    colours = _colour_steps(swept_values, colorscale)
    step = (high["x"] - low["x"]) + 2 * padding

    frames = []
    for entry in sweep_manifest:
        df = pd.read_parquet(entry["real_clustered_atoms_path"])
        inside = (
            (df["x"].between(low["x"], high["x"]))
            & (df["y"].between(low["y"], high["y"]))
            & (df["z"].between(low["z"], high["z"]))
        )
        if not include_unclustered:
            inside &= df["cluster_id"] != -1

        df = df.loc[inside].copy()
        if spread:
            df["x"] = df["x"] + step * swept_values.index(entry["swept_value"])
        df[swept_param_name] = str(entry["swept_value"])
        frames.append(df)

    df_to_plot = pd.concat(frames, ignore_index=True)

    fig = px.scatter_3d(
        df_to_plot, x="x", y="y", z="z",
        color=swept_param_name,
        category_orders={swept_param_name: [str(v) for v in swept_values]},
        color_discrete_map={str(v): c for v, c in colours.items()},
        hover_data=hover_columns(df_to_plot, hover),
    )
    fig.update_traces(marker_size=marker_size, opacity=0.8)

    if not spread:
        # graded sizes, so five copies of the same atom do not reduce to
        # whichever colour happened to be drawn last
        rank = {str(value): i for i, value in enumerate(swept_values)}
        for trace in fig.data:
            trace.marker.size = (
                marker_size + size_step * (len(swept_values) - 1 - rank[trace.name])
            )
    fig.update_layout(
        title=title or (
            f"clusters {list(selected_clusters)} at {swept_param_name}="
            f"{chosen_value}, across every {swept_param_name}"
        ),
        scene=dict(
            xaxis_title="x [nm]", yaxis_title="y [nm]", zaxis_title="z [nm]",
            aspectmode="data",
        ),
        legend_title_text=swept_param_name,
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
        margin=dict(l=20, r=20, t=60, b=20),
    )

    return fig
