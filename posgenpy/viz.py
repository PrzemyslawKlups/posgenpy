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
import uuid
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .corrections import cluster_alpha_shapes
from .fast_sweep import _clip_small_clusters, _n_min_lookup

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
    "true_positive_ratio", "n_min_reached", "n_min_status",
    "true_positive_ratio_atoms", "percolation_verdict", "giant_fraction",
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
    n_min=None,
    reference_line=None,
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

    n_min : dict or list, optional
        Mark each x value's n_min as a black tick across its violin: a
        `{x value: n_min}` mapping, or a list in ascending-x order - which is
        what `optimise_n_min` returns for a sweep manifest. The clusters below
        a tick are the ones that n_min throws away. Only meaningful when `y` is
        the size n_min is applied to (`size_core`).

    reference_line : float or (float, str), optional
        A dotted horizontal line at this y, across every violin, with an
        optional label at its right-hand end - e.g. `(1, "equivalent spheres
        touch")` on `closest_cluster_distance_norm`.

    Returns
    -------
    plotly.graph_objects.Figure
    """

    if n_min is not None and x is None:
        raise ValueError("n_min is marked per x value, so it needs an x")

    if n_min is not None and not isinstance(n_min, dict):
        x_values = sorted(stats[x].unique())
        if len(n_min) != len(x_values):
            raise ValueError(
                f"{len(n_min)} n_min values for {len(x_values)} values of "
                f"{x} ({x_values}); pass a {{{x} value: n_min}} mapping instead"
            )
        n_min = dict(zip(x_values, n_min))

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

    if n_min:
        fig.add_trace(
            go.Scatter(
                x=list(n_min), y=[int(n) for n in n_min.values()],
                mode="markers", name="n_min",
                marker=dict(symbol="line-ew", size=40, line=dict(width=3, color="black")),
                hovertemplate=f"{x} %{{x}}<br>n_min %{{y}}<extra></extra>",
            )
        )

    if reference_line is not None:
        value, label = (
            reference_line if isinstance(reference_line, tuple)
            else (reference_line, None)
        )
        # Placed across the plot in paper coordinates rather than with
        # fig.add_hline, whose "x domain" reference an older plotly.js (1.49,
        # for one) does not know and falls back to data x = 0..1 for: on an
        # Order axis (3-8) that is a stub left of every violin, while on a
        # dmax axis (0.2-0.7) it happens to look right.
        fig.add_shape(
            type="line", xref="paper", x0=0, x1=1, yref="y", y0=value, y1=value,
            line=dict(color="black", dash="dot", width=1.5),
        )
        if label:
            fig.add_annotation(
                xref="paper", x=1, xanchor="right", yref="y", y=value,
                yanchor="bottom", text=label, showarrow=False,
            )

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


# how a swept parameter is written in a subplot title
_PARAM_LABELS = {"order": "Order"}

# plotly.js, from the page: every 3D subplot follows the camera of the one being
# turned. plotly can match 2D axes between subplots but has no equivalent for 3D
# cameras, so this is the only way to get it.
_LINKED_SCENES_JS = """
(function () {
  var figure = __FIGURE__;
  var gd = document.getElementById("__DIV_ID__");

  function cameraKey(update) {
    return Object.keys(update).filter(function (key) {
      return /^scene\\d*\\.camera$/.test(key);
    })[0];
  }

  function copy(camera) {
    return JSON.parse(JSON.stringify(camera));
  }

  function linkCameras() {
    var scenes = Object.keys(figure.layout).filter(function (key) {
      return /^scene\\d*$/.test(key);
    });
    var following = false;

    // Live, while one scene is dragged. A relayout here would redraw every
    // scene from the saved layout - the dragged one too, whose camera is only
    // saved on release - so the others are moved by the same per-scene
    // viewport call plotly's own camera redraw makes.
    gd.on("plotly_relayouting", function (update) {
      var key = cameraKey(update);
      if (!key) return;
      scenes.forEach(function (id) {
        var scene = gd._fullLayout[id];
        if (key === id + ".camera" || !scene || !scene._scene) return;
        scene._scene.setViewport({
          camera: copy(update[key]), aspectratio: scene.aspectratio
        });
      });
    });

    // ...and on release or scroll-zoom, written into the layout, so a redraw
    // leaves them where they were turned to
    gd.on("plotly_relayout", function (update) {
      var key = cameraKey(update);
      if (following || !key) return;
      var follow = {};
      scenes.forEach(function (id) {
        if (key !== id + ".camera") follow[id + ".camera"] = copy(update[key]);
      });
      following = true;
      var done = function () { following = false; };
      Plotly.relayout(gd, follow).then(done, done);
    });
  }

  function draw() {
    Plotly.newPlot(gd, figure.data, figure.layout, figure.config || {})
      .then(linkCameras);
  }

  if (window.Plotly && window.Plotly.version === "__PLOTLY_VERSION__") {
    draw();
  } else {
    var script = document.createElement("script");
    script.src = "https://cdn.plot.ly/plotly-__PLOTLY_VERSION__.min.js";
    script.onload = draw;
    script.onerror = function () {
      gd.textContent = "plotly.js could not be loaded from " + script.src +
        " (offline?). fig.show(renderer=...) draws this figure without " +
        "linked cameras.";
    };
    document.head.appendChild(script);
  }
})();
"""


def _linked_scenes_html(fig):
    """`fig` as an HTML snippet whose 3D subplots turn together."""

    from plotly.offline import get_plotlyjs_version

    div_id = f"posgenpy-{uuid.uuid4().hex}"
    # "</" would let a tooltip's "</extra>" - or anything else - close the
    # script tag early; "<\/" is the same string to JSON
    figure_json = fig.to_json().replace("</", "<\\/")

    script = (
        _LINKED_SCENES_JS
        .replace("__DIV_ID__", div_id)
        .replace("__PLOTLY_VERSION__", get_plotlyjs_version())
        .replace("__FIGURE__", figure_json)
    )
    return f'<div id="{div_id}"></div>\n<script type="text/javascript">{script}</script>'


class LinkedScenesFigure(go.Figure):
    """A figure of 3D subplots whose cameras turn together.

    Shown in a notebook, it is drawn as HTML carrying the few lines of
    JavaScript that do the linking (plotly.js is loaded from its CDN the first
    time, so this needs to be online). Anywhere else, or with an explicit
    renderer - `fig.show(renderer="vscode")` - it is a normal, unlinked
    `go.Figure`, as it is for `write_html`, `to_dict` and every `update_*`.
    """

    def show(self, *args, **kwargs):
        if args or kwargs or not _in_notebook():
            return super().show(*args, **kwargs)

        from IPython.display import HTML, display
        display(HTML(_linked_scenes_html(self)))

    def _ipython_display_(self):
        self.show()


def _in_notebook():
    """True inside a Jupyter kernel (VSCode, JupyterLab), which runs HTML output."""
    try:
        from IPython import get_ipython
    except ImportError:
        return False
    return getattr(get_ipython(), "kernel", None) is not None


def compare_selected_clusters_across_param(
    sweep_manifest: list,
    swept_param_name: str,
    chosen_value,
    selected_clusters: list,
    n_min=None,
    hover: list = None,
    padding: float = 0.5,
    include_unclustered: bool = False,
    spread: bool = False,
    ion_colors: dict = None,
    n_cols: int = 3,
    colorscale: str = "Viridis",
    marker_size: int = 3,
    size_step: float = 1.6,
    title: str = None,
    figsize: tuple = None,
):
    """One cluster, at every value of a swept parameter.

    Takes the neighbourhood of `selected_clusters` (as they are at
    `chosen_value`) and draws what every other swept value made of the same
    piece of the tip. Two ways to lay that out:

    - **`spread=True`**: one 3D subplot per value, `n_cols` to a row, titled
      by the value ("Order = 5") and coloured by ion (`ion_colors`). Every
      subplot frames the same box, and turning or zooming one turns all the
      others with it (see `LinkedScenesFigure`), so a cluster breaking in two
      as the Order rises is one blob in one subplot and two in the next,
      seen from the same side. Clicking an ion in the legend hides it
      everywhere.
    - **`spread=False`** (default): every value in one view, **each value in
      its own colour**. Every value draws the *same atoms* at the *same
      coordinates*, so drawn flat they would hide each other; the markers are
      graded by size instead - low parameter values large, high values small
      and drawn on top - so a big dot with a small dot inside it is an atom
      clustered at both values, a bare small dot is one that only the higher
      value picked up, and a bare big dot is one the higher value dropped.
      One legend entry per value, so clicking one hides it and two can be
      compared directly.

    Cluster ids are per-run, so the same physical cluster generally has a
    different id at each swept value; each atom's tooltip carries its id at
    its own value, which is what says *what* it merged into.

    Parameters
    ----------
    n_min : int, list, dict or None
        Draw each swept value's clusters as they are after *its own* n_min:
        the list `optimise_n_min` returns for this manifest, a single number,
        or anything else `fast_sweep.load_grid_cluster_table` takes. Clusters
        below it count as unclustered, so the 2-3 atom noise clusters a sweep
        keeps are not drawn beside real ones as if they were part of them.
        Sizes are counted in the per-atom file - core+border atoms on the fast
        path, the size n_min is chosen on. None draws every cluster.
    padding : float
        Extra nm around the chosen clusters' bounding box, so a cluster that
        grows at a looser parameter is not cut off at the box edge.
    include_unclustered : bool
        Also draw core-species atoms the search left unclustered in that box
        (grey, under their own legend entry, when `spread=True`).
    spread : bool
        One linked 3D subplot per swept value instead of one overlaid view.
    ion_colors : dict, optional
        Ion -> colour, as `get_ion_colours_from_rrng` returns, for
        `spread=True`. Ions it does not name get plotly's default colours.
    n_cols : int
        Subplots per row when `spread=True`.
    colorscale : str
        The per-value colour ramp when `spread=False`.
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
    lookup = _n_min_lookup(n_min, manifest=sweep_manifest)
    n_mins = {e["swept_value"]: lookup(e["order"], e["dmax"]) for e in sweep_manifest}

    chosen_n_min = n_mins[chosen_entries[0]["swept_value"]]
    sizes = chosen_df["cluster_id"].value_counts()
    below = sorted(sizes.index[sizes < chosen_n_min].tolist())
    if below:
        warnings.warn(
            f"clusters {below} are below n_min={chosen_n_min} at "
            f"{swept_param_name}={chosen_value}, so they are drawn as "
            "unclustered there"
        )

    frames = {}
    for entry in sweep_manifest:
        df = pd.read_parquet(entry["real_clustered_atoms_path"])
        # clipped over the whole cell, before the box: a cluster poking into
        # the box is judged on all of its atoms, not the ones inside
        df["cluster_id"] = _clip_small_clusters(
            df["cluster_id"].to_numpy(), n_mins[entry["swept_value"]]
        )
        inside = (
            (df["x"].between(low["x"], high["x"]))
            & (df["y"].between(low["y"], high["y"]))
            & (df["z"].between(low["z"], high["z"]))
        )
        if not include_unclustered:
            inside &= df["cluster_id"] != -1

        df = df.loc[inside].copy()
        df[swept_param_name] = str(entry["swept_value"])
        frames[entry["swept_value"]] = df

    title = title or (
        f"clusters {list(selected_clusters)} at {swept_param_name}="
        f"{chosen_value}, across every {swept_param_name}"
    )

    if spread:
        label = _PARAM_LABELS.get(swept_param_name, swept_param_name)
        titles = {
            value: f"{label} = {value}"
            + (f", n_min {n_mins[value]}" if n_min is not None else "")
            for value in frames
        }
        return _clusters_in_linked_subplots(
            frames, titles, low, high, ion_colors, n_cols, hover,
            marker_size, title, figsize,
        )

    colours = _colour_steps(swept_values, colorscale)
    df_to_plot = pd.concat(frames.values(), ignore_index=True)

    fig = px.scatter_3d(
        df_to_plot, x="x", y="y", z="z",
        color=swept_param_name,
        category_orders={swept_param_name: [str(v) for v in swept_values]},
        color_discrete_map={str(v): c for v, c in colours.items()},
        hover_data=hover_columns(df_to_plot, hover),
    )
    fig.update_traces(marker_size=marker_size, opacity=0.8)

    # graded sizes, so five copies of the same atom do not reduce to whichever
    # colour happened to be drawn last
    rank = {str(value): i for i, value in enumerate(swept_values)}
    for trace in fig.data:
        trace.marker.size = (
            marker_size + size_step * (len(swept_values) - 1 - rank[trace.name])
        )

    figsize = figsize or (900, 750)
    fig.update_layout(
        title=title,
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


def _clusters_in_linked_subplots(
    frames, titles, low, high, ion_colors, n_cols, hover,
    marker_size, title, figsize,
):
    """`compare_selected_clusters_across_param(spread=True)`'s figure."""

    values = sorted(frames)
    n_cols = max(1, min(n_cols, len(values)))
    n_rows = -(-len(values) // n_cols)

    # one colour per ion across every subplot: px would otherwise number the
    # ions ion_colors does not name afresh in each subplot, in whatever order
    # they turn up there
    colour_key = {
        value: np.where(df["cluster_id"] != -1, df["ion"].astype(str), "unclustered")
        for value, df in frames.items()
    }
    keys = sorted(set(np.concatenate(list(colour_key.values()))) - {"unclustered"})
    fallback = iter(px.colors.qualitative.Plotly * (len(keys) // 10 + 1))
    colours = {
        key: (ion_colors or {}).get(key) or next(fallback) for key in keys
    }
    colours["unclustered"] = "lightgrey"

    specs = [
        [{"type": "scene"} if row * n_cols + col < len(values) else None
         for col in range(n_cols)]
        for row in range(n_rows)
    ]
    fig = LinkedScenesFigure(make_subplots(
        rows=n_rows, cols=n_cols, specs=specs,
        subplot_titles=[titles[value] for value in values],
        horizontal_spacing=0.02, vertical_spacing=min(0.06, 0.5 / max(n_rows, 2)),
    ))

    # every subplot frames the same box at the same proportions, so an atom is
    # in the same place in all of them
    extent = high - low
    box = dict(
        xaxis=dict(range=[low["x"], high["x"]], title="x [nm]"),
        yaxis=dict(range=[low["y"], high["y"]], title="y [nm]"),
        zaxis=dict(range=[low["z"], high["z"]], title="z [nm]"),
        aspectmode="manual",
        aspectratio=(extent / extent.max()).to_dict(),
    )

    in_legend = set()
    for i, value in enumerate(values):
        row, col = divmod(i, n_cols)
        df = _with_cluster_label(frames[value].assign(colour=colour_key[value]))

        if len(df):
            panel = px.scatter_3d(
                df, x="x", y="y", z="z", color="colour",
                color_discrete_map=colours,
                category_orders={"colour": list(colours)},
                labels={"colour": "ion"},
                hover_name="cluster",
                hover_data=hover_columns(df, hover, skip=("cluster", "ion")),
            )
            for trace in panel.data:
                trace.showlegend = trace.name not in in_legend
                in_legend.add(trace.name)
                fig.add_trace(trace, row=row + 1, col=col + 1)

        fig.update_scenes(box, row=row + 1, col=col + 1)

    fig.update_traces(marker_size=marker_size, opacity=0.8)

    figsize = figsize or (330 * n_cols + 120, 330 * n_rows + 90)
    fig.update_layout(
        title=title,
        legend_title_text="ion",
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
        margin=dict(l=10, r=10, t=80, b=10),
    )

    return fig


# the two reference compositions drawn beside the clusters
_POOLED_LABEL = "all clusters"
_MATRIX_LABEL = "matrix"
_REFERENCE_STYLE = {
    _POOLED_LABEL: dict(symbol="star", color="#222222"),
    _MATRIX_LABEL: dict(symbol="diamond", color="#8c8c8c"),
}


def _reference_compositions(cluster_stats, unclustered_stats, ions):
    """{label: counts per ion} for all clusters pooled and for the matrix.

    Pooled means every cluster's atoms summed before normalising, so a big
    cluster weighs in proportion to its atoms - the composition of "the
    clustered material", not the average cluster.
    """

    references = {_POOLED_LABEL: cluster_stats[ions].sum()}
    if unclustered_stats is not None and len(unclustered_stats):
        matrix = unclustered_stats.reindex(columns=ions, fill_value=0)
        references[_MATRIX_LABEL] = matrix.iloc[0]

    return {
        label: counts.astype(float) for label, counts in references.items()
        if counts.sum() > 0
    }


def plot_cluster_ternary(
    cluster_stats,
    ions=("Mn", "Ni", "Si"),
    unclustered_stats=None,
    size_column: str = "cluster_size",
    title: str = None,
    figsize: tuple = (760, 660),
    max_marker_size: float = 26,
):
    """Ternary of three elements: one point per cluster, plus the matrix.

    Each cluster is placed by those three elements alone - their atom counts
    normalised to 100% - and sized by `size_column`. Edge clusters get the
    second of `EDGE_CLUSTER_COLORS`, as in `plot_cluster_alpha_shapes`. Two
    reference points sit among them: **all clusters** pooled (every cluster's
    atoms summed, so a big cluster weighs more) and, with `unclustered_stats`,
    the **matrix**.

    A cluster holding none of the three has no place on the triangle and is
    left out; the title says how many were. An element with no column at all
    (not ranged in this dataset) is counted as zero, with a warning.

    Parameters
    ----------
    cluster_stats : pandas.DataFrame
        One row per cluster with a count column per element -
        `analyse_clusters` output, after the corrections.
    ions : three element names, corners a, b, c of the triangle.
    unclustered_stats : pandas.DataFrame, optional
        `summarise_unclustered` output, for the matrix point.
    size_column : str
        What marker area is proportional to; falls back to `size`.
    """

    ions = list(ions)
    if len(ions) != 3:
        raise ValueError(f"a ternary needs exactly three elements, not {ions}")

    stats = cluster_stats.copy()
    missing = [ion for ion in ions if ion not in stats.columns]
    if missing:
        warnings.warn(
            f"no {missing} column in these cluster stats - counted as zero",
            stacklevel=2,
        )
        for ion in missing:
            stats[ion] = 0
    if size_column not in stats.columns:
        size_column = "size"

    total = stats[ions].sum(axis=1)
    placed = stats.loc[total > 0].copy()
    percent = placed[ions].div(total[total > 0], axis=0) * 100
    n_left_out = int((total == 0).sum())

    sizes = placed[size_column].astype(float)
    marker_size = max_marker_size * np.sqrt(sizes / sizes.max()) if len(sizes) else sizes
    marker_size = np.maximum(marker_size, 4)

    corner_lines = [
        f"{ion}: %{{{axis}:.1f}}% (%{{customdata[{i + 2}]}} atoms)"
        for i, (ion, axis) in enumerate(zip(ions, "abc"))
    ]

    fig = go.Figure()

    is_edge = (
        placed["is_edge_cluster"].astype(bool) if "is_edge_cluster" in placed.columns
        else pd.Series(False, index=placed.index)
    )
    groups = [(False, "interior clusters"), (True, "edge clusters")]
    if not is_edge.any():
        groups = [(False, "clusters")]

    for edge, name in groups:
        rows = is_edge == edge
        if not rows.any():
            continue
        fig.add_trace(go.Scatterternary(
            a=percent.loc[rows, ions[0]], b=percent.loc[rows, ions[1]],
            c=percent.loc[rows, ions[2]],
            mode="markers", name=name,
            marker=dict(
                size=marker_size[rows.to_numpy()], color=EDGE_CLUSTER_COLORS[edge],
                opacity=0.7, line=dict(width=0.5, color="white"),
            ),
            customdata=np.column_stack([
                placed.loc[rows, "cluster_id"], placed.loc[rows, size_column],
                placed.loc[rows, ions].to_numpy(),
            ]),
            hovertemplate=(
                "<b>Cluster %{customdata[0]}</b><br>"
                f"{size_column}: " "%{customdata[1]}<br>"
                + "<br>".join(corner_lines) + "<extra></extra>"
            ),
        ))

    for label, counts in _reference_compositions(stats, unclustered_stats, ions).items():
        share = counts / counts.sum() * 100
        fig.add_trace(go.Scatterternary(
            a=[share[ions[0]]], b=[share[ions[1]]], c=[share[ions[2]]],
            mode="markers", name=label,
            marker=dict(size=16, line=dict(width=1.5, color="white"),
                        **_REFERENCE_STYLE[label]),
            customdata=[[label, int(counts.sum()), *counts.astype(int)]],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                + "<br>".join(corner_lines) + "<extra></extra>"
            ),
        ))

    title = title or f"{' - '.join(ions)} composition of each cluster"
    if n_left_out:
        title += f"  ({n_left_out} clusters with none of the three left out)"

    fig.update_layout(
        title=title,
        ternary=dict(
            sum=100,
            aaxis=dict(title=ions[0], ticksuffix="%"),
            baxis=dict(title=ions[1], ticksuffix="%"),
            caxis=dict(title=ions[2], ticksuffix="%"),
        ),
        legend=dict(itemsizing="constant"),
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
    )

    return fig


def plot_cluster_composition_bars(
    cluster_stats,
    ions: list,
    unclustered_stats=None,
    ion_colors: dict = None,
    size_column: str = "cluster_size",
    title: str = None,
    figsize: tuple = (1000, 480),
    max_ticks: int = 40,
):
    """Composition of every cluster, one 100%-stacked column each.

    One column per cluster, smallest on the left, each normalised to 100%, so
    every column is the same height and only the mix changes - the question
    it answers is "does the composition change with cluster size". Beside
    them, in their own panel: **all clusters** pooled (every cluster's atoms
    summed, so a big cluster weighs more) and, with `unclustered_stats`, the
    **matrix**.

    The percentages are over `ions` only. Pass the elements `size_column`
    counts - for `cluster_size`, the core and bulk ions minus
    `elements_to_exclude` - and a column's ions add up to exactly the size
    printed under it. Elements without a column, or with no atoms anywhere,
    are dropped from the legend.

    Parameters
    ----------
    cluster_stats : pandas.DataFrame
        One row per cluster with a count column per element.
    ions : list
        Elements to show, stacked bottom to top in this order.
    unclustered_stats : pandas.DataFrame, optional
        `summarise_unclustered` output, for the matrix column.
    ion_colors : dict, optional
        Element -> colour, as `get_ion_colours_from_rrng` returns.
    size_column : str
        What the columns are ordered by and labelled with; falls back to
        the sum over `ions`.
    max_ticks : int
        At most this many clusters get a size label, so a tip with hundreds
        of clusters stays readable (every column is still hoverable).
    """

    stats = cluster_stats.copy()
    ions = [ion for ion in ions if ion in stats.columns]
    references = _reference_compositions(stats, unclustered_stats, ions)
    in_matrix = references.get(_MATRIX_LABEL, pd.Series(0.0, index=ions))
    ions = [ion for ion in ions if stats[ion].sum() > 0 or in_matrix[ion] > 0]
    if not ions:
        raise ValueError("none of the given elements has any atoms in these clusters")
    references = _reference_compositions(stats, unclustered_stats, ions)

    stats["_total"] = stats[ions].sum(axis=1)
    if size_column not in stats.columns:
        size_column = "_total"
    stats = stats.loc[stats["_total"] > 0].sort_values([size_column, "cluster_id"])
    percent = stats[ions].div(stats["_total"], axis=0) * 100

    # categories are the ids (unique), the tick text is the size (not)
    labels = "Cluster " + stats["cluster_id"].astype(str)
    step = max(1, int(np.ceil(len(stats) / max_ticks)))

    fig = make_subplots(
        rows=1, cols=2, shared_yaxes=True,
        column_widths=[0.86, 0.14], horizontal_spacing=0.02,
    )
    colours = ion_colors or {}

    for ion in ions:
        fig.add_trace(go.Bar(
            x=labels, y=percent[ion], name=ion, legendgroup=ion,
            marker_color=colours.get(ion),
            customdata=np.column_stack(
                [stats["cluster_id"], stats[size_column], stats[ion]]
            ),
            hovertemplate=(
                "<b>Cluster %{customdata[0]}</b><br>"
                f"{'size' if size_column == '_total' else size_column}: "
                "%{customdata[1]}<br>"
                f"{ion}: " "%{y:.1f}% (%{customdata[2]} atoms)<extra></extra>"
            ),
        ), row=1, col=1)

        if references:
            fig.add_trace(go.Bar(
                x=list(references),
                y=[counts[ion] / counts.sum() * 100 for counts in references.values()],
                name=ion, legendgroup=ion, showlegend=False,
                marker_color=colours.get(ion),
                customdata=[[int(counts[ion])] for counts in references.values()],
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    f"{ion}: " "%{y:.1f}% (%{customdata[0]} atoms)<extra></extra>"
                ),
            ), row=1, col=2)

    fig.update_xaxes(
        tickmode="array", tickvals=labels.iloc[::step],
        ticktext=stats[size_column].iloc[::step].astype(int).astype(str),
        title=(
            f"clusters, smallest to largest (labelled with "
            f"{'size' if size_column == '_total' else size_column})"
        ),
        row=1, col=1,
    )
    fig.update_yaxes(range=[0, 100], title="composition [at.%]", row=1, col=1)
    fig.update_layout(
        title=title or "composition of each cluster, smallest to largest",
        barmode="stack",
        bargap=0.1,
        legend=dict(traceorder="reversed"),
        template="simple_white",
        width=figsize[0],
        height=figsize[1],
    )

    return fig
