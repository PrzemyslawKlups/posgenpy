"""Order/dmax parameter sweeps, n_min optimisation, and sweep-result plotting.

Replaces posgen-XML-file-based lookups (posgenpy.py's xml_files convention)
with an in-memory + persisted-to-disk "sweep manifest": a list of small dicts
pointing at the parquet files run_relabel_sweep wrote for each swept value.
"""
import json
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .cluster import build_position_tree, run_relabel_sweep

__all__ = [
    "seconds_to_hhmmss",
    "run_parameter_sweep",
    "optimise_n_min",
    "plot_real_cluster_ratio_across_swept_param",
    "plot_real_clusters_across_swept_param",
    "plot_cluster_composition_across_swept_param_absolute",
    "plot_cluster_composition_across_swept_param_relative",
]


def seconds_to_hhmmss(seconds: float):
    hh = int(np.floor(seconds / 3600))
    hh_rem = seconds % 3600
    mm = int(np.floor(hh_rem / 60))
    ss = int(np.floor(seconds % 60))
    return f"{hh}h {mm}m {ss}s"


def run_parameter_sweep(
    pos_df,
    rrng_df,
    core_ions: list,
    bulk_ions: list,
    swept_param: str,
    swept_values: list,
    fixed_dmax: float = None,
    fixed_order: int = None,
    n_min: int = 2,
    envelope_radius: float = None,
    d_erosion: float = None,
    n_relabel_runs: int = 20,
    output_dir: str = ".",
    prefix: str = "sweep",
    dclassify=None,
    cached_shuffles=None,
    relabel_seed: int = 0,
    verbose: bool = False,
) -> list:
    """
    Sweep either `order` (at a fixed dmax) or `dmax` (at a fixed order),
    running a full real+relabelled cluster search (run_relabel_sweep) at each
    swept value.

    `order` is posgen's knn / Order throughout - the same number that appears
    in a posgen run's `_dmax0.5_order_N` filename - so a sweep over
    [3, 4, 5, 6, 7, 8] lines up directly with posgen runs of those orders.

    The KDTree over all atom positions is built once here and reused across
    every swept value and every relabel run.

    If envelope_radius / d_erosion are not given, they default to
    (current dmax value) / 2 for each swept value, matching the previous
    posgen-driven workflow's dbulk="dmax/2", derode="dmax/2" convention.

    Returns a sweep manifest (list of dicts, one per swept value) - also
    persisted as JSON to `{output_dir}/{prefix}_{swept_param}_sweep_manifest.json`
    so a long sweep can be resumed (an existing run's parquet files are
    reused, not recomputed) or reloaded later without rerunning anything.
    """

    if swept_param not in ("order", "dmax"):
        raise ValueError('swept_param must be "order" or "dmax"')
    if swept_param == "order" and fixed_dmax is None:
        raise ValueError("fixed_dmax is required when sweeping order")
    if swept_param == "dmax" and fixed_order is None:
        raise ValueError("fixed_order is required when sweeping dmax")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tree, tree_positions = build_position_tree(pos_df)

    manifest = []
    tic = perf_counter()

    for value in swept_values:

        if swept_param == "order":
            this_dmax = fixed_dmax
            this_order = value
        else:
            this_dmax = value
            this_order = fixed_order

        this_envelope_radius = (
            envelope_radius if envelope_radius is not None else this_dmax / 2
        )
        this_d_erosion = (
            d_erosion if d_erosion is not None else this_dmax / 2
        )

        run_name = f"{prefix}_{swept_param}_{value}".replace(".", "p")

        entry = run_relabel_sweep(
            pos_df, rrng_df, core_ions, bulk_ions,
            dmax=this_dmax, order=this_order, n_min=n_min,
            envelope_radius=this_envelope_radius, d_erosion=this_d_erosion,
            n_relabel_runs=n_relabel_runs, run_name=run_name,
            output_dir=output_dir, relabel_seed=relabel_seed,
            tree=tree, tree_positions=tree_positions,
            cached_shuffles=cached_shuffles, dclassify=dclassify,
            verbose=verbose,
        )
        entry["swept_param"] = swept_param
        entry["swept_value"] = value
        manifest.append(entry)

        toc = perf_counter()
        print(f"{swept_param}={value}: done in {seconds_to_hhmmss(toc - tic)}")
        tic = toc

    manifest_path = output_dir / f"{prefix}_{swept_param}_sweep_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    return manifest


SIZE_METRICS = ("core", "total", "total_excluding")


def _cluster_sizes(stats, size_metric, entry, exclude_ions, core_species=None):
    """Reduce a per-cluster stats table to the size series being thresholded.

    posgen reports two different notions of "cluster size" and the old
    workflow used both: its `sizedist.txt` counts core+border atoms
    (`Cluster Size(core)`), which is what the real-cluster-count graphs were
    filtered on, while the old `optimise_n_min` read a `Cluster Size` column
    that the notebook had built as all ions minus an excluded matrix element.
    Thresholding the wrong one shifts every curve upwards, increasingly so at
    larger n_min.
    """
    if stats.empty:
        return pd.Series(dtype=int)

    if size_metric == "total":
        return stats["size"]

    if size_metric == "core":
        if "size_core" not in stats.columns:
            named = core_species or entry.get("core_species", [])
            core = [c for c in named if c in stats.columns]
            if not core:
                raise ValueError(
                    "size_metric='core' needs a size_core column or per-ion "
                    f"columns; {entry.get('run_name', '?')} has neither. It was "
                    "written by an older version - pass core_species=[...], "
                    "re-run that sweep, or use size_metric='total'."
                )
            return stats[core].sum(axis=1)
        return stats["size_core"]

    if size_metric == "total_excluding":
        drop = [c for c in (exclude_ions or []) if c in stats.columns]
        return stats["size"] - stats[drop].sum(axis=1) if drop else stats["size"]

    raise ValueError(f"size_metric must be one of {SIZE_METRICS}")


def _load_sweep_sizes(
    sweep_manifest: list, size_metric: str = "core", exclude_ions: list = None,
    core_species: list = None, need_random: bool = True,
) -> dict:
    """Load each swept value's real + random per-cluster sizes into memory.

    Cheap: rows = clusters (tens-hundreds), not atoms.
    """

    data = {}

    for entry in sweep_manifest:
        real_stats = pd.read_parquet(entry["real_cluster_stats_path"])
        real_sizes = _cluster_sizes(
            real_stats, size_metric, entry, exclude_ions, core_species
        )

        if need_random:
            random_dir = Path(entry["random_dir"])
            frames = [
                _cluster_sizes(
                    pd.read_parquet(p), size_metric, entry, exclude_ions, core_species
                )
                for p in sorted(random_dir.glob("run_*.parquet"))
            ]
            random_sizes = (
                pd.concat(frames, ignore_index=True) if frames else pd.Series(dtype=int)
            )
        else:
            random_sizes = pd.Series(dtype=int)

        data[entry["swept_value"]] = {
            "real_sizes": real_sizes,
            "random_sizes": random_sizes,
            "n_runs": entry["n_relabel_runs"],
        }

    return data


def optimise_n_min(
    sweep_manifest: list,
    n_min_values=np.arange(1, 500, 1),
    confidence_threshold: float = 0.95,
    size_metric: str = "core",
    exclude_ions: list = None,
    core_species: list = None,
    verbose: bool = True,
) -> tuple:
    """
    For each swept value, find the smallest n_min at which the true-positive
    ratio (real clusters - mean random clusters) / real clusters reaches
    confidence_threshold, using each run's saved per-cluster sizes.

    `size_metric` selects what n_min is compared against - "core" (posgen's
    `Cluster Size(core)`, core+border atoms), "total" (every atom), or
    "total_excluding" with `exclude_ions` (the old notebook's `Cluster Size`,
    which excluded the matrix element, typically Fe).

    `verbose` prints the chosen n_min per swept value; pass False when calling
    this per cell of a grid, where the running commentary is noise.

    Returns (optimised_n_min_values, optimised_real_clusters, optimised_random_clusters),
    one entry per swept value in sweep_manifest, in order.
    """

    vprint = print if verbose else lambda *args, **kwargs: None

    data = _load_sweep_sizes(
        sweep_manifest, size_metric, exclude_ions, core_species
    )

    optimised_n_min_values = []
    optimised_real_clusters = []
    optimised_random_clusters = []

    for entry in sweep_manifest:
        d = data[entry["swept_value"]]
        real_sizes, random_sizes, n_runs = d["real_sizes"], d["random_sizes"], d["n_runs"]

        n_min_found = False
        local_real_clusters = 0
        local_random_clusters = 0.0

        for n_min in n_min_values:

            local_real_clusters = int((real_sizes >= n_min).sum())
            local_random_clusters = (
                (random_sizes >= n_min).sum() / n_runs if n_runs else 0.0
            )

            if local_real_clusters == 0:
                vprint(f"TPR cannot be reached for {entry['run_name']}")
                n_min_found = True
                break

            local_tpr = (local_real_clusters - local_random_clusters) / local_real_clusters

            if local_tpr >= confidence_threshold:
                vprint(
                    f"{entry['run_name']} | Nmin:{n_min} | TPR: {round(local_tpr, 2)}"
                    f"\t(real:{local_real_clusters} | random:{round(local_random_clusters)})"
                )
                n_min_found = True
                break

        optimised_n_min_values.append(
            int(n_min) if n_min_found else int(max(n_min_values))
        )
        optimised_real_clusters.append(local_real_clusters)
        optimised_random_clusters.append(local_random_clusters)

    return optimised_n_min_values, optimised_real_clusters, optimised_random_clusters


def plot_real_cluster_ratio_across_swept_param(
    sweep_manifest: list,
    swept_parameter_name: str,
    n_min_values: list,
    threshold: float = 0.95,
    graph_alpha: float = 0.6,
    figsize: tuple = (9, 5),
    size_metric: str = "core",
    exclude_ions: list = None,
    core_species: list = None,
):
    """Plot (N_real - N_random)/N_real vs swept parameter, for each n_min.

    `size_metric` defaults to "core" - posgen's `Cluster Size(core)` - so the
    curves match the equivalent posgen-era graph. See `optimise_n_min`.
    """

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions)
    swept_values = sorted(data.keys())

    colors = plt.cm.viridis(np.linspace(0, 1, len(n_min_values)))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    for i, n in enumerate(n_min_values):
        real_local = np.array(
            [(data[v]["real_sizes"] >= n).sum() for v in swept_values], dtype=float
        )
        random_local = np.array(
            [
                (data[v]["random_sizes"] >= n).sum() / data[v]["n_runs"]
                if data[v]["n_runs"] else 0.0
                for v in swept_values
            ]
        )

        # guard div-by-zero explicitly (real_local == 0 -> NaN, not a warning)
        y_values = np.where(
            real_local == 0, np.nan, (real_local - random_local) / real_local
        )

        ax.plot(
            swept_values, y_values, 'o-', color=colors[i],
            label="$N_{min}=%s$" % n, alpha=graph_alpha
        )

    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("$(N_{det}-N_{rand})/N_{det}$")
    ax.plot(
        [min(swept_values), max(swept_values)], [threshold, threshold], 'r-',
        label=f"{round(threshold * 100)}%"
    )
    ax.legend(bbox_to_anchor=(1.05, 1))
    ax.set_ylim(top=1.1)
    ax.set_title("Fraction of Real Clusters")

    return fig, ax


def plot_real_clusters_across_swept_param(
    sweep_manifest: list,
    swept_parameter_name: str,
    n_min_values: list,
    volume: float = 0,
    graph_alpha: float = 0.7,
    figsize: tuple = (9, 5),
    size_metric: str = "core",
    exclude_ions: list = None,
    core_species: list = None,
):
    """Plot real cluster count (or number density, if volume given) vs swept parameter.

    `size_metric` defaults to "core" - posgen's `Cluster Size(core)`, the
    quantity its `sizedist.txt` reports - so these curves match the
    equivalent posgen-era graph. Passing "total" counts every atom in the
    cluster instead, which keeps more clusters at every n_min.
    """

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions)
    swept_values = sorted(data.keys())
    show_volume = volume != 0

    colors = plt.cm.viridis(np.linspace(0, 1, len(n_min_values)))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    for i, n in enumerate(n_min_values):
        real_local = np.array(
            [(data[v]["real_sizes"] >= n).sum() for v in swept_values], dtype=float
        )
        y_values = real_local / volume if show_volume else real_local
        ax.plot(
            swept_values, y_values, 'o-', color=colors[i],
            label="$N_{min}=%s$" % n, alpha=graph_alpha
        )

    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("Number density [$m^{-3}$]" if show_volume else "Number of identified clusters")
    ax.legend(bbox_to_anchor=(1.3, 1))
    ax.set_title(f"Real Clusters and {swept_parameter_name} Values")

    return fig, ax


# Everything analyse_clusters emits that is NOT a per-element atom count.
# `size_core` belongs here: it is a count of atoms already counted under their
# own element columns, so leaving it out double-counts every core-species atom
# in the composition denominator below.
_CLUSTER_STATS_METADATA_COLUMNS = [
    "cluster_id", "x", "y", "z", "size", "size_core", "radius_of_gyration",
    "aspect_ratio", "principal_axis_1", "principal_axis_2", "principal_axis_3",
    "bbox_x", "bbox_y", "bbox_z", "bbox_volume", "cluster_density",
    "equivalent_spherical_radius",
    # added by corrections.remove_edge_clusters / find_knn_cluster_distance
    "convexhull_volume", "convexhull_area", "is_edge_cluster",
    "closest_cluster_distance", "closest_cluster_id",
]


def _prepare_composition_data(
    sweep_manifest: list, core_ions: list, n_min_values_for_swept_values: list,
    exclude_ions=None,
) -> dict:

    exclude_ions = exclude_ions or []
    corrected = {}

    for entry, n_min in zip(sweep_manifest, n_min_values_for_swept_values):
        stats = pd.read_parquet(entry["real_cluster_stats_path"])
        stats = stats.drop(columns=[c for c in exclude_ions if c in stats.columns])
        stats = stats[stats["size"] >= n_min]

        value = entry["swept_value"]
        corrected[value] = {ion: stats[ion].sum() if ion in stats.columns else 0
                             for ion in core_ions}

        element_cols = [c for c in stats.columns if c not in _CLUSTER_STATS_METADATA_COLUMNS]
        corrected[value]["total clusters"] = stats[element_cols].to_numpy().sum()

    return corrected


def plot_cluster_composition_across_swept_param_absolute(
    sweep_manifest: list,
    core_ions: list,
    n_min_values_for_swept_values: list,
    swept_parameter_name: str,
    exclude_ions=None,
    graph_alpha: float = 0.7,
    ion_colors: dict = None,
    figsize: tuple = (9, 5),
):
    ion_colors = ion_colors or {}
    corrected = _prepare_composition_data(
        sweep_manifest, core_ions, n_min_values_for_swept_values, exclude_ions
    )
    swept_values = [entry["swept_value"] for entry in sweep_manifest]

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    for ion in core_ions:
        y_values = [
            corrected[v][ion] / corrected[v]["total clusters"] * 100
            if corrected[v]["total clusters"] else 0
            for v in swept_values
        ]
        ax.plot(
            swept_values, y_values, 'o-', alpha=graph_alpha, label=ion,
            color=ion_colors.get(ion)
        )

    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("Cluster composition [at.%]")
    ax.legend(bbox_to_anchor=(1.2, 1))
    ax.set_title(f"Absolute Cluster Composition and {swept_parameter_name} Values")

    return fig, ax


def plot_cluster_composition_across_swept_param_relative(
    sweep_manifest: list,
    swept_parameter_name: str,
    core_ions: list,
    n_min_values_for_swept_values: list,
    exclude_ions=None,
    graph_alpha: float = 1.0,
    ion_colors: dict = None,
    figsize: tuple = (9, 5),
):
    ion_colors = ion_colors or {}
    corrected = _prepare_composition_data(
        sweep_manifest, core_ions, n_min_values_for_swept_values, exclude_ions
    )
    swept_values = [entry["swept_value"] for entry in sweep_manifest]

    core_ions_sum = {
        v: sum(corrected[v][ion] for ion in core_ions) for v in swept_values
    }

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    if len(swept_values) > 1:
        column_width = np.abs(np.diff(swept_values)).mean() * 0.75
    else:
        column_width = 0.75

    next_values = None
    for i, ion in enumerate(core_ions):
        y_values = np.array([
            corrected[v][ion] / core_ions_sum[v] * 100 if core_ions_sum[v] else 0
            for v in swept_values
        ])
        y_values_bar = np.ones(len(swept_values)) * 100 if i == 0 else next_values
        ax.bar(
            swept_values, y_values_bar, alpha=graph_alpha, label=ion,
            width=column_width, color=ion_colors.get(ion)
        )
        next_values = y_values_bar - y_values

    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("Corrected relative cluster composition [at.%]")
    ax.legend(bbox_to_anchor=(1.05, 1))
    ax.set_title(f"Relative Cluster Composition and {swept_parameter_name} Values")

    return fig, ax
