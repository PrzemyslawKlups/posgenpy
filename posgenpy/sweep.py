"""Order/dmax parameter sweeps, n_min optimisation, and sweep-result plotting.

Replaces posgen-XML-file-based lookups (posgenpy.py's xml_files convention)
with an in-memory + persisted-to-disk "sweep manifest": a list of small dicts
pointing at the parquet files run_relabel_sweep wrote for each swept value.
"""
import json
import warnings
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
    "n_min_report",
    "percolation_report",
    "plot_real_cluster_ratio_across_swept_param",
    "plot_real_clusters_across_swept_param",
    "plot_atom_ratio_across_swept_param",
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

    `random_sizes` pools every control run - which is what the cluster *counts*
    want, since they are divided by `n_runs` afterwards. `random_runs` keeps the
    same sizes split per run, because a per-run statistic (the largest cluster
    in one control dataset - see `percolation_report`) is not recoverable from
    the pool.
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
            frames = []
            random_sizes = pd.Series(dtype=int)

        data[entry["swept_value"]] = {
            "real_sizes": real_sizes,
            "random_sizes": random_sizes,
            "random_runs": frames,
            "n_runs": entry["n_relabel_runs"],
        }

    return data


def _atom_ratio(d):
    """`(A_real - <A_rand>) / A_real` over every cluster in the saved tables.

    Atoms, not clusters: a percolating control contributes exactly one cluster
    no matter how much solute it swallowed, so the count-based ratio can look
    healthy - can be *pushed past the threshold by the percolation itself* -
    where this one is at or below zero.
    """

    atoms_real = float(d["real_sizes"].sum())
    if not atoms_real:
        return np.nan
    n_runs = d["n_runs"]
    atoms_random = float(d["random_sizes"].sum()) / n_runs if n_runs else 0.0
    return (atoms_real - atoms_random) / atoms_real


N_MIN_STATUSES = ("reached", "unreachable", "range_exhausted")

_N_MIN_REPORT_COLUMNS = [
    "swept_value", "run_name", "n_min", "status",
    "n_real_clusters", "n_random_clusters",
    "true_positive_ratio", "true_positive_ratio_atoms", "n_min_returned",
]


def _n_min_rows(sweep_manifest, data, n_min_values, confidence_threshold) -> list:
    """The n_min search itself, over sizes already in memory.

    One row per swept value. The search has three exits and this is the only
    place that knows which one was taken:

    | status | |
    |---|---|
    | `"reached"` | the true-positive ratio crossed `confidence_threshold` |
    | `"unreachable"` | the real clusters ran out first, so no n_min in *any* range can get there |
    | `"range_exhausted"` | no n_min in `n_min_values` got there |

    `n_min` is NaN unless the status is `"reached"` - a threshold that was never
    crossed has no answer. `n_min_returned` is the integer the legacy
    `optimise_n_min` 3-tuple reports for that cell regardless of status, kept so
    the two agree exactly.
    """

    n_min_values = list(n_min_values)
    if not n_min_values:
        raise ValueError("n_min_values is empty; nothing to search")

    rows = []

    for entry in sweep_manifest:
        d = data[entry["swept_value"]]
        real_sizes, random_sizes, n_runs = d["real_sizes"], d["random_sizes"], d["n_runs"]

        status = "range_exhausted"
        n_min = int(max(n_min_values))
        local_real_clusters = 0
        local_random_clusters = 0.0

        for candidate in n_min_values:

            local_real_clusters = int((real_sizes >= candidate).sum())
            local_random_clusters = float(
                (random_sizes >= candidate).sum() / n_runs if n_runs else 0.0
            )

            if local_real_clusters == 0:
                status, n_min = "unreachable", int(candidate)
                break

            local_tpr = (local_real_clusters - local_random_clusters) / local_real_clusters

            if local_tpr >= confidence_threshold:
                status, n_min = "reached", int(candidate)
                break

        rows.append(
            {
                "swept_value": entry["swept_value"],
                "run_name": entry.get("run_name"),
                "n_min": float(n_min) if status == "reached" else np.nan,
                "status": status,
                "n_real_clusters": local_real_clusters,
                "n_random_clusters": local_random_clusters,
                "true_positive_ratio": (
                    (local_real_clusters - local_random_clusters) / local_real_clusters
                    if local_real_clusters
                    else np.nan
                ),
                "true_positive_ratio_atoms": _atom_ratio(d),
                "n_min_returned": int(n_min),
            }
        )

    return rows


def _warn_about_n_min_statuses(rows, confidence_threshold):
    """Surface the two failure exits even when nothing is printing.

    `optimise_n_min`'s `verbose=False` is exactly how `summarise_grid` calls it,
    so a printed message is no message at all in the case that matters most.
    """

    failed = [r for r in rows if r["status"] != "reached"]
    if not failed:
        return

    detail = "; ".join(
        f"{r['run_name'] or r['swept_value']}: {r['status']}" for r in failed
    )
    warnings.warn(
        f"the true-positive ratio never reached {confidence_threshold} for "
        f"{len(failed)} of {len(rows)} swept values ({detail}). The n_min "
        "reported for those is where the search stopped, not an answer - see "
        "n_min_report() for the status and percolation_report() for the usual "
        "cause (a randomised control whose solute network percolates).",
        stacklevel=3,
    )


def n_min_report(
    sweep_manifest: list,
    n_min_values=np.arange(1, 500, 1),
    confidence_threshold: float = 0.95,
    size_metric: str = "core",
    exclude_ions: list = None,
    core_species: list = None,
    warn: bool = True,
) -> pd.DataFrame:
    """`optimise_n_min` with its reasoning shown: one row per swept value.

    Same search, same arguments; what it adds is `status` - whether the
    threshold was genuinely crossed (`"reached"`), the real clusters ran out
    first (`"unreachable"`), or no n_min in `n_min_values` got there
    (`"range_exhausted"`). Only a `"reached"` row carries an n_min; the other
    two are NaN, because there is no answer to report.

    | column | |
    |---|---|
    | `n_min` | smallest n_min reaching `confidence_threshold`; NaN if it never did |
    | `status` | which of the three exits above was taken |
    | `n_real_clusters`, `n_random_clusters` | counts where the search stopped |
    | `true_positive_ratio` | `(N_real - <N_rand>) / N_real` there |
    | `true_positive_ratio_atoms` | the same comparison on *atoms* - `(A_real - <A_rand>) / A_real` over every cluster in the saved tables, no n_min applied. At or below zero the control clusters as much solute as the tip does and the swept value has **no discriminating power**, whatever the count-based ratio says |
    | `n_min_returned` | what the legacy `optimise_n_min` 3-tuple reports here |

    Raises a warning for every non-`"reached"` row unless `warn=False`.
    """

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions, core_species)
    rows = _n_min_rows(sweep_manifest, data, n_min_values, confidence_threshold)

    if warn:
        _warn_about_n_min_statuses(rows, confidence_threshold)

    return pd.DataFrame(rows, columns=_N_MIN_REPORT_COLUMNS)


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

    **A returned n_min is not evidence the threshold was reached.** Two of the
    three ways the search can end return a plausible-looking integer that means
    "this is where I gave up": the real clusters ran out, or `n_min_values` did.
    Both now raise a warning, and `n_min_report` - which this is a thin wrapper
    over - names which one happened per swept value.
    """

    vprint = print if verbose else lambda *args, **kwargs: None

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions, core_species)
    rows = _n_min_rows(sweep_manifest, data, n_min_values, confidence_threshold)

    for row in rows:
        if row["status"] == "unreachable":
            vprint(f"TPR cannot be reached for {row['run_name']}")
        elif row["status"] == "reached":
            vprint(
                f"{row['run_name']} | Nmin:{row['n_min_returned']} | "
                f"TPR: {round(row['true_positive_ratio'], 2)}"
                f"\t(real:{row['n_real_clusters']} | "
                f"random:{round(row['n_random_clusters'])})"
            )

    _warn_about_n_min_statuses(rows, confidence_threshold)

    return (
        [row["n_min_returned"] for row in rows],
        [row["n_real_clusters"] for row in rows],
        [row["n_random_clusters"] for row in rows],
    )


PERCOLATION_VERDICTS = ("ok", "near-threshold", "percolating", "unknown")

_PERCOLATION_REPORT_COLUMNS = [
    "swept_value", "run_name", "largest_random_cluster",
    "largest_random_cluster_max", "giant_fraction", "denominator",
    "denominator_atoms", "largest_real_cluster", "verdict",
]


def _percolation_rows(sweep_manifest, data, warn_fraction, fail_fraction) -> list:
    """Size of the largest cluster in each randomised control, per swept value.

    See `percolation_report` for what the numbers mean; this is the calculation
    with the loading and the warning stripped out, so `summarise_grid` can run
    it on sizes it has already read.
    """

    if not 0 < warn_fraction <= fail_fraction:
        raise ValueError("need 0 < warn_fraction <= fail_fraction")

    rows = []

    for entry in sweep_manifest:
        d = data[entry["swept_value"]]
        runs, real_sizes = d["random_runs"], d["real_sizes"]

        maxima = [float(sizes.max()) if len(sizes) else 0.0 for sizes in runs]

        # the natural denominator is every core-species atom in the tip, since
        # that is what the control's giant cluster is eating into. Manifests
        # written before n_core_atoms existed fall back to the solute that the
        # control actually clustered, which is smaller and shrinks with order -
        # so the fraction it gives is inflated, and is floored below.
        n_core_atoms = entry.get("n_core_atoms")
        if n_core_atoms:
            denominator = float(n_core_atoms)
            denominator_name = "n_core_atoms"
        else:
            totals = [float(sizes.sum()) for sizes in runs]
            denominator = float(np.mean(totals)) if totals else 0.0
            denominator_name = "clustered_solute"

        largest_mean = float(np.mean(maxima)) if maxima else np.nan
        largest_max = float(np.max(maxima)) if maxima else np.nan
        giant_fraction = largest_mean / denominator if denominator else np.nan

        if not maxima or not np.isfinite(giant_fraction):
            verdict = "unknown"
        elif giant_fraction >= fail_fraction:
            verdict = "percolating"
        elif giant_fraction >= warn_fraction:
            verdict = "near-threshold"
        else:
            verdict = "ok"

        if verdict in ("percolating", "near-threshold") and denominator_name == "clustered_solute":
            pooled = d["random_sizes"]
            median = float(pooled.median()) if len(pooled) else 0.0
            if largest_mean <= max(50.0, 5.0 * median):
                # a "giant" cluster of 12 atoms in a control that clustered 30
                # is the degenerate denominator talking, not a percolating network
                verdict = "ok"

        rows.append(
            {
                "swept_value": entry["swept_value"],
                "run_name": entry.get("run_name"),
                "largest_random_cluster": largest_mean,
                "largest_random_cluster_max": largest_max,
                "giant_fraction": giant_fraction,
                "denominator": denominator_name,
                "denominator_atoms": denominator,
                "largest_real_cluster": (
                    float(real_sizes.max()) if len(real_sizes) else 0.0
                ),
                "verdict": verdict,
            }
        )

    return rows


def percolation_report(
    sweep_manifest: list,
    size_metric: str = "core",
    exclude_ions: list = None,
    core_species: list = None,
    warn_fraction: float = 0.10,
    fail_fraction: float = 0.30,
    warn: bool = True,
) -> pd.DataFrame:
    """Is the randomised control still a control, or has it percolated?

    `optimise_n_min` compares *counts of clusters*, and a percolating network
    contributes exactly **one** cluster however much solute it swallowed. Where
    the control's solute network spans the tip, the comparison silently stops
    measuring anything: `<N_rand>` collapses towards 1, the true-positive ratio
    sails past the threshold, and a confident-looking n_min comes back from a
    parameter region with no discriminating power at all. The usual cause is a
    `dmax` large enough that the mean solute-solute degree in the *randomised*
    tip exceeds the ~2.7 continuum percolation threshold - a parameter choice,
    not a bug in the search, but one nothing used to detect.

    This reads it off the per-cluster size tables the sweep already wrote, so it
    **runs on existing sweeps with nothing re-run**.

    | column | |
    |---|---|
    | `largest_random_cluster` | mean over control runs of the largest cluster |
    | `largest_random_cluster_max` | the worst single run - near the threshold the run-to-run spread is several-fold, and the mean hides it |
    | `giant_fraction` | `largest_random_cluster / denominator_atoms` |
    | `denominator` | `"n_core_atoms"` (the tip's core-species atoms, from the manifest) or `"clustered_solute"` (the fallback, see below) |
    | `largest_real_cluster` | the same for the real tip, for contrast: a genuine precipitate stays put as the order rises, a percolating background does not |
    | `verdict` | `"ok"` / `"near-threshold"` / `"percolating"`, or `"unknown"` where there are no control runs |

    Measured regimes on a percolating dataset separate cleanly - sub-critical at
    or below 0.01, critical around 0.21, percolating 0.65-0.88 - so the default
    0.10 / 0.30 cut sits in empty space on both sides. They are keyword
    arguments; the physics is not universal.

    **Older manifests.** `n_core_atoms` was added after this function; without
    it the denominator falls back to the solute each control run actually
    clustered. That fallback *inflates* the fraction at high order, since the
    denominator shrinks faster than the numerator does, so a non-`"ok"` verdict
    is additionally gated on the giant cluster exceeding
    `max(50, 5 x median cluster size)`. Re-running the sweep records the real
    denominator and removes the guesswork.
    """

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions, core_species)
    rows = _percolation_rows(sweep_manifest, data, warn_fraction, fail_fraction)

    if warn:
        flagged = [r for r in rows if r["verdict"] in ("percolating", "near-threshold")]
        if flagged:
            detail = "; ".join(
                f"{r['run_name'] or r['swept_value']}: {r['verdict']} "
                f"({r['giant_fraction']:.2f})"
                for r in flagged
            )
            warnings.warn(
                f"the randomised control's largest cluster spans a large "
                f"fraction of the tip's solute at {len(flagged)} of {len(rows)} "
                f"swept values ({detail}). Cluster *counts* there - including "
                "the true-positive ratio optimise_n_min thresholds - are not "
                "meaningful; lower dmax or raise the order.",
                stacklevel=2,
            )

    return pd.DataFrame(rows, columns=_PERCOLATION_REPORT_COLUMNS)


_PERCOLATION_SHADING = {
    "percolating": ("tab:red", 0.13, "control percolates"),
    "near-threshold": ("tab:orange", 0.10, "control near threshold"),
}


def _shade_percolating(ax, swept_values, rows):
    """Grey out the swept values whose control percolates.

    Without this a reader has no way to see that the left third of one of these
    x-axes is arithmetic on a network that spans the tip.
    """

    verdicts = {row["swept_value"]: row["verdict"] for row in rows}
    if not any(v in _PERCOLATION_SHADING for v in verdicts.values()):
        return

    x = np.asarray(swept_values, dtype=float)
    if len(x) > 1:
        midpoints = (x[:-1] + x[1:]) / 2
        edges = np.concatenate(
            [[x[0] - (x[1] - x[0]) / 2], midpoints, [x[-1] + (x[-1] - x[-2]) / 2]]
        )
    else:
        edges = np.array([x[0] - 0.5, x[0] + 0.5])

    labelled = set()
    for i, value in enumerate(swept_values):
        shading = _PERCOLATION_SHADING.get(verdicts.get(value))
        if shading is None:
            continue
        color, alpha, label = shading
        ax.axvspan(
            edges[i], edges[i + 1], color=color, alpha=alpha, zorder=0, lw=0,
            label=None if label in labelled else label,
        )
        labelled.add(label)


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
    flag_percolation: bool = True,
):
    """Plot (N_real - N_random)/N_real vs swept parameter, for each n_min.

    `size_metric` defaults to "core" - posgen's `Cluster Size(core)` - so the
    curves match the equivalent posgen-era graph. See `optimise_n_min`.

    `flag_percolation` shades the swept values where `percolation_report` says
    the randomised control's solute network has percolated. The curve is still
    drawn there, because it is still what the arithmetic gives - but a ratio
    computed against a control that is one tip-spanning cluster is not a
    measurement, and it is otherwise indistinguishable from a good one.
    """

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions, core_species)
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

    if flag_percolation:
        _shade_percolating(
            ax, swept_values, _percolation_rows(sweep_manifest, data, 0.10, 0.30)
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
    flag_percolation: bool = True,
):
    """Plot real cluster count (or number density, if volume given) vs swept parameter.

    `size_metric` defaults to "core" - posgen's `Cluster Size(core)`, the
    quantity its `sizedist.txt` reports - so these curves match the
    equivalent posgen-era graph. Passing "total" counts every atom in the
    cluster instead, which keeps more clusters at every n_min.

    `flag_percolation` shades swept values where the randomised control has
    percolated - see `plot_real_cluster_ratio_across_swept_param`.
    """

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions, core_species)
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

    if flag_percolation:
        _shade_percolating(
            ax, swept_values, _percolation_rows(sweep_manifest, data, 0.10, 0.30)
        )

    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("Number density [$m^{-3}$]" if show_volume else "Number of identified clusters")
    ax.legend(bbox_to_anchor=(1.3, 1))
    ax.set_title(f"Real Clusters and {swept_parameter_name} Values")

    return fig, ax


def plot_atom_ratio_across_swept_param(
    sweep_manifest: list,
    swept_parameter_name: str,
    graph_alpha: float = 0.9,
    figsize: tuple = (9, 5),
    size_metric: str = "core",
    exclude_ions: list = None,
    core_species: list = None,
    flag_percolation: bool = True,
):
    """Plot the atom-based true-positive ratio vs swept parameter.

    `(A_real - <A_rand>) / A_real`, where `A` is the solute the search put into
    clusters. The companion to `plot_real_cluster_ratio_across_swept_param`,
    and the one to read first: it counts atoms, so a percolating control cannot
    hide inside it as a single cluster and cannot be *pushed past the threshold
    by the percolation itself*. One curve, no n_min - every cluster in the
    saved tables.

    **At or below zero the control clusters as much solute as the tip does**,
    and no n_min chosen in that region means anything.
    """

    data = _load_sweep_sizes(sweep_manifest, size_metric, exclude_ions, core_species)
    swept_values = sorted(data.keys())

    y_values = [_atom_ratio(data[v]) for v in swept_values]

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    ax.plot(swept_values, y_values, 'o-', color="tab:blue", alpha=graph_alpha,
            label="$(A_{det}-A_{rand})/A_{det}$")
    ax.axhline(0, color="r", lw=1, label="no discriminating power")

    if flag_percolation:
        _shade_percolating(
            ax, swept_values, _percolation_rows(sweep_manifest, data, 0.10, 0.30)
        )

    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("$(A_{det}-A_{rand})/A_{det}$")
    ax.legend(bbox_to_anchor=(1.05, 1))
    ax.set_title("Fraction of Clustered Solute Above the Control")

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
