"""Accelerated (order x dmax) grid sweep of the DBSCAN/kNN stage.

Why this is faster
------------------
`sweep.run_parameter_sweep` redoes *everything* for every (parameter value,
dataset) pair: re-range the pos file, rebuild the species partition, run
DBSCAN, grow the envelope one cluster at a time, then erode it. With P
parameter values and R relabelled control datasets that is P x (R+1) full
pipelines - and the two most expensive stages, envelope and erosion, do not
change which clusters exist at all, only how much bulk gets wrapped around
them.

This module inverts the loops: datasets outside, parameters inside. Per
dataset it computes, **once**:

* the distance from every core-species atom to its k-th nearest core-species
  neighbour, for every k in the swept orders (one kNN query at k = max order);
* every pair of core-species atoms closer than max(dmax) (one radius query).

Those two arrays are all the core-classification and linking steps need, at
*any* (order, dmax) in the grid:

* an atom is a core point at (order, dmax) exactly when its order-th
  neighbour distance is < dmax - that is the same order statistic for every
  dmax, so it is computed once and merely re-thresholded;
* two core points link exactly when their pair distance is < dmax - a prefix
  of the distance-sorted pair list, so raising dmax only ever *adds* edges.

Each parameter cell is then a handful of numpy operations plus one
connected-components call on a graph already in memory. Going from one dmax
to the next answers "would this atom have been reached at the next parameter
value?" by construction, without re-walking the tree.

Orders are not treated incrementally
------------------------------------
Only the dmax axis is walked monotonically, and only in the direction where
that is sound: in DBSCAN, raising dmax adds core points and adds links, so
clusters can merge but never split. Raising the *order* does the opposite - it
takes core points away, and a cluster held together by a thin neck of core
atoms falls into two when that neck fails the knn test.

So no state is carried between cells. `CoreNeighbourhood.labels` recomputes
which atoms are core points and re-runs connected components from scratch every
time, over whatever subgraph that leaves. A split at a higher order therefore
comes out as two clusters, exactly as `run_clustering` reports it; there is no
incremental union that could hold them together (`test_raising_the_order_splits_a_cluster`).
What is shared across cells is only the *evidence* - the neighbour distances and
the candidate pairs - never a partition.

What it deliberately does not do
--------------------------------
No envelope, no erosion, no bulk species. The clusters it reports are
core+border atoms only. That is not a loss for the sweep stage, because it is
already what the sweep stage reads: `optimise_n_min` and the two real-vs-random
plots default to ``size_metric="core"``, i.e. posgen's `Cluster Size(core)`,
and the envelope/erosion steps provably cannot change that number - they only
ever add or remove *bulk-tree* atoms, and core-species atoms are never in the
bulk tree. So this module reproduces the sweep's own outputs, and the results
it writes are drop-in compatible with `sweep.optimise_n_min`,
`sweep.plot_real_clusters_across_swept_param` and
`sweep.plot_real_cluster_ratio_across_swept_param`.

It is *not* a substitute for `cluster.run_clustering` at the final chosen
parameters, where the envelope and erosion decide the cluster compositions and
volumes. Run the full pipeline once, at the end, on the parameters this picks.

The one place the two can disagree is a border atom within dmax of core atoms
of two different clusters: which one claims it depends on traversal order.
This module gives it to the *nearest* core atom; sklearn's DBSCAN gives it to
whichever its BFS reaches first. That is the same class of ambiguity that
already puts this engine 1 atom out of 45,194 apart from posgen itself
(see ../README). `validate_fast_grid` measures it on real data.
"""
import json
import multiprocessing as mp
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .cluster import analyse_clusters
from .io import build_ion_formula_map
from .sweep import optimise_n_min

__all__ = [
    "extract_core_atoms",
    "core_code_from_mc",
    "split_core_code",
    "relabelled_core_atoms",
    "CoreNeighbourhood",
    "run_fast_parameter_sweep",
    "run_fast_grid_sweep",
    "select_sweep",
    "summarise_grid",
    "load_sweep_size_table",
    "validate_fast_grid",
]


# ===========================================================================
# Datasets: the core-species atoms of the real tip, and of each relabelling


def extract_core_atoms(pos_df, rrng_df, core_ions):
    """Indices + ion codes of the atoms ranged as one of `core_ions`.

    Equivalent to ``assign_rrng_to_pos(...)`` followed by an ``isin(core_ions)``
    filter - same first-match-wins range semantics (libatomprobe's
    `RangeFile::getIonID`, inclusive bounds) - but it never materialises the
    per-atom object columns, which on a 20M-atom tip cost more memory and time
    than the whole cluster search.

    Returns
    -------
    index : (n_core,) int64
        Row positions in `pos_df`.
    code : (n_core,) int
        Index into `core_ions` for each of those atoms.
    """

    return split_core_code(
        core_code_from_mc(pos_df["mc"].to_numpy(), rrng_df, list(core_ions))
    )


def core_code_from_mc(mc, rrng_df, core_ions):
    """Per-atom core-species code for every atom: -1, or an index into `core_ions`.

    Equivalent to ``assign_rrng_to_pos(...)`` followed by an ``isin(core_ions)``
    filter - same first-match-wins range semantics (libatomprobe's
    `RangeFile::getIonID`, inclusive bounds) - but as one small integer per atom
    rather than a column of Python strings.

    Keeping the *full-length* array, rather than only the core atoms' codes, is
    what makes relabelling cheap: see `relabelled_core_atoms`.
    """

    core_ions = list(core_ions)
    name_to_code = {name: i for i, name in enumerate(core_ions)}

    # int8 covers the 127 core species nobody has; the array is one per atom,
    # copied once per control dataset, so its width is worth minding
    dtype = np.int8 if len(core_ions) < np.iinfo(np.int8).max else np.int32
    code = np.full(len(mc), -1, dtype=dtype)
    unassigned = np.ones(len(mc), dtype=bool)

    for low, high, ion in zip(rrng_df["mc_begin"], rrng_df["mc_end"], rrng_df["ion"]):
        # first match wins: a later range never overwrites an earlier one
        matched = (mc >= low) & (mc <= high) & unassigned
        this_code = name_to_code.get(ion)
        if this_code is not None:
            code[matched] = this_code
        unassigned &= ~matched

    return code


def split_core_code(full_code):
    """A full-length core code array as (row positions, codes) of its core atoms."""
    index = np.flatnonzero(full_code >= 0)
    return index, full_code[index]


def relabelled_core_atoms(full_code, rng):
    """One randomised control dataset, as core-atom positions + ion codes.

    This *is* the standard relabelling - shuffle the pos file's `mc` column and
    re-range every atom - not a substitute for it. Ranging is elementwise in
    `mc`, so ranging a permuted `mc` gives the same thing as permuting the
    ranged labels:

        R(mc[pi])[i] = R(mc[pi[i]]) = R(mc)[pi[i]] = (R(mc)[pi])[i]

    So shuffling the per-atom code array in place is exactly `mc`-shuffle plus
    re-ranging, for the permutation the shuffle drew - it just skips walking all
    the mass ranges over all the atoms again, which is where that work actually
    went (~64 ms per dataset on a 1.85M-atom tip, against ~20 ms for the
    shuffle).

    `full_code` is left untouched; the caller keeps the real tip's codes.
    """

    shuffled = full_code.copy()
    rng.shuffle(shuffled)
    return split_core_code(shuffled)


# ===========================================================================
# Per-dataset precompute, reused by every parameter cell


class CoreNeighbourhood:
    """Everything the core classification and linking need, for one dataset.

    Built once per dataset; `labels(order, dmax)` then answers any cell of the
    (order x dmax) grid from it.

    Parameters
    ----------
    positions : (n, 3) array
        Positions of this dataset's core-species atoms only.
    orders : sequence[int]
        Every order (posgen's knn) the grid will ask for.
    max_dmax : float
        The largest dmax the grid will ask for. Pairs further apart than this
        can never link at any cell, so they are never materialised.
    workers : int
        Threads for the kNN query (-1 = all cores). Set to 1 when several
        datasets are already being processed in parallel processes.
    """

    def __init__(self, positions, orders, max_dmax, workers=-1):

        self.positions = np.ascontiguousarray(positions, dtype=np.float64)
        self.orders = sorted({int(o) for o in orders})
        self.max_dmax = float(max_dmax)
        self.n = len(self.positions)

        if self.n == 0:
            self._core_distance = np.zeros((len(self.orders), 0))
            self.edge_i = np.zeros(0, dtype=np.int32)
            self.edge_j = np.zeros(0, dtype=np.int32)
            self.edge_d = np.zeros(0)
            return

        tree = cKDTree(self.positions)

        # --- order statistics -------------------------------------------
        # posgen's Order counts the point itself, so the order-th neighbour
        # (including self at distance 0) is column order-1 of the kNN query.
        # A point is a core point at (order, dmax) iff that distance < dmax,
        # which is what makes every dmax free once this is computed.
        k_max = max(self.orders)
        if k_max > 1:
            distances, _ = tree.query(self.positions, k=k_max, workers=workers)
            distances = np.atleast_2d(distances)
        else:
            distances = np.zeros((self.n, 1))

        self._core_distance = np.empty((len(self.orders), self.n))
        for slot, order in enumerate(self.orders):
            # order <= 1: every point is its own order-1 neighbour at distance
            # 0, i.e. always a core point, matching min_samples=1
            self._core_distance[slot] = (
                np.zeros(self.n) if order <= 1 else distances[:, order - 1]
            )
        del distances

        # --- link candidates --------------------------------------------
        # Every pair that could ever link, at any dmax in the grid. Sorted by
        # distance so each cell's active edges are the prefix d < dmax.
        #
        # This is the memory high-water mark and it is set by max(dmax) for the
        # whole grid, so the indices are narrowed to int32 - there is no
        # 2-billion-atom core species - for 16 bytes an edge instead of 24.
        # The distances stay float64: they are tested against dmax with a strict
        # `<`, and rounding them to float32 would move edges across that
        # boundary, which is the one thing this path must not do.
        pairs = tree.query_pairs(r=self.max_dmax, output_type="ndarray")
        edge_i = pairs[:, 0].astype(np.int32)
        edge_j = pairs[:, 1].astype(np.int32)
        del pairs
        edge_d = np.linalg.norm(
            self.positions[edge_i] - self.positions[edge_j], axis=1
        )

        by_distance = np.argsort(edge_d, kind="stable")
        self.edge_i = edge_i[by_distance]
        self.edge_j = edge_j[by_distance]
        self.edge_d = edge_d[by_distance]

    @property
    def n_edges(self):
        return len(self.edge_d)

    def core_distance(self, order):
        """Distance to the `order`-th nearest core atom, counting the atom itself."""
        try:
            slot = self.orders.index(int(order))
        except ValueError:
            raise ValueError(
                f"order={order} was not precomputed; this neighbourhood covers "
                f"{self.orders}"
            ) from None
        return self._core_distance[slot]

    def _cell(self, order, dmax):
        """Core classification, linked components and border candidates at one cell."""

        if dmax > self.max_dmax:
            raise ValueError(
                f"dmax={dmax} exceeds the max_dmax={self.max_dmax} this "
                "neighbourhood was built for; pairs beyond it were never "
                "collected"
            )

        # strict <, matching posgen's `<` on squared distances (cluster._shrink)
        is_core = self.core_distance(order) < dmax

        if self.n == 0 or not is_core.any():
            empty = np.zeros(0, dtype=np.int32)
            return is_core, np.zeros(self.n, dtype=np.int32), empty, empty

        # edges are distance-sorted, so "closer than dmax" is a prefix
        cut = int(np.searchsorted(self.edge_d, dmax, side="left"))
        edge_i = self.edge_i[:cut]
        edge_j = self.edge_j[:cut]

        core_i = is_core[edge_i]
        core_j = is_core[edge_j]

        # --- core-to-core linking ---------------------------------------
        both_core = core_i & core_j
        graph = coo_matrix(
            (
                np.ones(int(both_core.sum()), dtype=np.int8),
                (edge_i[both_core], edge_j[both_core]),
            ),
            shape=(self.n, self.n),
        )
        _, component = connected_components(graph, directed=False, return_labels=True)

        # --- border candidates ------------------------------------------
        # edges from a non-core core-species atom to a core one, still in
        # ascending distance order
        one_core = core_i ^ core_j
        i_is_core = core_i[one_core]
        border = np.where(i_is_core, edge_j[one_core], edge_i[one_core])
        claimer = np.where(i_is_core, edge_i[one_core], edge_j[one_core])

        return is_core, component, border, claimer

    def labels(self, order, dmax):
        """Cluster label per core atom at (order, dmax); -1 = unclustered.

        Reproduces `DBSCAN(eps=Dmax, min_samples=order)` on the core-species
        atoms, including border-point attachment, using only the precomputed
        arrays. Labels are 0..n_clusters-1.
        """

        is_core, component, border, claimer = self._cell(order, dmax)

        labels = np.full(self.n, -1, dtype=np.int32)
        if not is_core.any():
            return labels

        labels[is_core] = component[is_core]

        # A non-core core-species atom within dmax of a core atom joins that
        # atom's cluster. Where several clusters could claim it, the nearest
        # core atom wins: the border list is distance-ascending, so writing the
        # assignments in reverse leaves the closest one in place.
        if len(border):
            labels[border[::-1]] = component[claimer[::-1]]

        return _compact_labels(labels)

    def ambiguous_border_atoms(self, order, dmax):
        """Border atoms that two or more clusters could equally claim.

        These are the only atoms on which any correct DBSCAN implementation is
        free to disagree: they sit within dmax of core atoms of more than one
        cluster, and the standard algorithm gives each to whichever cluster its
        traversal happens to reach first. This engine gives it to the nearest
        core atom instead. Reported by `validate_fast_grid` as the ceiling on
        how far the two paths can differ.
        """

        _, component, border, claimer = self._cell(order, dmax)
        if len(border) == 0:
            return 0

        claimed = pd.Series(component[claimer], index=border)
        return int((claimed.groupby(level=0).nunique() > 1).sum())


def _compact_labels(labels):
    """Renumber the surviving labels to a dense 0..n_clusters-1 range."""
    clustered = labels != -1
    if not clustered.any():
        return np.full(len(labels), -1, dtype=np.int32)

    used, inverse = np.unique(labels[clustered], return_inverse=True)
    out = np.full(len(labels), -1, dtype=np.int32)
    out[clustered] = inverse.astype(np.int32)
    return out


def _apply_n_min(labels, n_min):
    """Drop clusters with fewer than `n_min` core+border atoms.

    posgen's `<sizeclip nmin>` counts core+border atoms and applies the clip
    before the envelope, so on this core-only path it is the whole story.
    """
    if n_min <= 1:
        return labels

    clustered = labels != -1
    if not clustered.any():
        return labels

    sizes = np.bincount(labels[clustered])
    too_small = np.isin(labels, np.flatnonzero(sizes < n_min)) & clustered
    labels = labels.copy()
    labels[too_small] = -1
    return _compact_labels(labels)


# ===========================================================================
# Per-cell result tables


def _sizes_frame(labels, ion_code, ion_names):
    """Per-cluster size table for a control run - the same schema
    `cluster._summarise_sizes` writes, so `sweep._load_sweep_sizes` reads it
    unchanged.

    `size` and `size_core` are equal here: with no envelope there are no bulk
    atoms in a cluster, so every atom counted is a core-species one.
    """

    clustered = labels != -1
    columns = ["cluster_id", "size", "size_core", *ion_names]

    if not clustered.any():
        return pd.DataFrame({c: pd.Series(dtype="int64") for c in columns})

    label = labels[clustered]
    n_clusters = int(label.max()) + 1
    n_ions = len(ion_names)

    size = np.bincount(label, minlength=n_clusters)
    per_ion = np.bincount(
        label.astype(np.int64) * n_ions + ion_code[clustered],
        minlength=n_clusters * n_ions,
    ).reshape(n_clusters, n_ions)

    out = pd.DataFrame(
        {
            "cluster_id": np.arange(n_clusters, dtype=np.int64),
            "size": size,
            "size_core": size,
        }
    )
    for slot, name in enumerate(ion_names):
        out[name] = per_ion[:, slot]

    return out


def _atoms_frame(pos_df, source_index, labels, ion_code, ion_names):
    """Per-atom table for the real run: the core-species atoms with their
    cluster ids, in the same shape `viz.show_clusters_in_pos` expects.

    Unlike the full pipeline's output this holds core+border atoms only - the
    enveloped bulk was never computed - so a 3D view of a cluster from here
    shows its core species, not its bulk shell.
    """

    columns = [c for c in ("ionIdx", "x", "y", "z", "mc") if c in pos_df.columns]
    out = pos_df.iloc[source_index][columns].copy()
    out["ion"] = pd.Categorical.from_codes(ion_code, categories=list(ion_names))
    out["cluster_id"] = labels
    return out.reset_index(drop=True)


# ===========================================================================
# The sweep


def _format_value(value):
    return str(value).replace(".", "p").replace("-", "m")


def _cell_run_name(prefix, order, dmax):
    return f"{prefix}_order_{_format_value(order)}_dmax_{_format_value(dmax)}"


def _cell_paths(output_dir, run_name):
    output_dir = Path(output_dir)
    return {
        "real_cluster_stats_path": output_dir / f"{run_name}_real_cluster_stats.parquet",
        "real_clustered_atoms_path": output_dir
        / f"{run_name}_real_clustered_atoms.parquet",
        "random_dir": output_dir / f"{run_name}_random",
    }


_SHARED = {}


def _random_dataset_job(run_index):
    """One control dataset: build its neighbourhood, write every cell's sizes.

    Module-level (and reading module-level state) so it can run in a forked
    worker process without pickling the tip's positions to it.
    """

    state = _SHARED
    started = perf_counter()

    cached_shuffles = state["cached_shuffles"]
    if cached_shuffles is not None:
        # the caller supplied the permutations, so this control dataset is the
        # same one run_relabel_sweep would build from mc[cached_shuffles[r]]
        index, code = split_core_code(state["full_code"][cached_shuffles[run_index]])
    else:
        rng = np.random.default_rng(state["relabel_seed"] + run_index)
        index, code = relabelled_core_atoms(state["full_code"], rng)

    neighbourhood = CoreNeighbourhood(
        state["positions"][index],
        state["orders"],
        state["max_dmax"],
        workers=state["workers"],
    )

    for order, dmax in state["cells"]:
        run_path = (
            Path(state["random_dirs"][(order, dmax)]) / f"run_{run_index:03d}.parquet"
        )
        if state["resume"] and run_path.exists():
            continue

        labels = _apply_n_min(neighbourhood.labels(order, dmax), state["n_min"])
        _sizes_frame(labels, code, state["core_ions"]).to_parquet(
            run_path, index=False
        )

    return run_index, perf_counter() - started


def run_fast_grid_sweep(
    pos_df,
    rrng_df,
    core_ions: list,
    orders,
    dmax_values,
    bulk_ions: list = None,
    n_relabel_runs: int = 20,
    output_dir: str = ".",
    prefix: str = "sweep",
    n_min: int = 1,
    relabel_seed: int = 0,
    cached_shuffles=None,
    save_real_atoms: bool = True,
    atomic_volume: float = 0.0118,
    n_jobs: int = 1,
    resume: bool = True,
    verbose: bool = True,
) -> list:
    """Sweep the whole (order x dmax) grid, one pass per dataset.

    Runs the core-classification/linking half of the cluster search - no
    envelope, no erosion - at every combination of `orders` x `dmax_values`,
    for the real tip and for `n_relabel_runs` randomised control datasets. Cost
    is one KD-tree + one kNN query + one radius query **per dataset**, not per
    parameter value, so adding parameter values is nearly free (see the module
    docstring for why that is exact rather than approximate).

    Parameters
    ----------
    core_ions : list[str]
        posgen's `<core><typelist>`. Note posgen treats each ranged label as a
        distinct type, so "C2" is not covered by listing "C".
    orders, dmax_values : sequence
        The grid. Both may be a single value; the returned manifest then
        already carries `swept_param`/`swept_value` for the other axis, ready
        for `sweep.optimise_n_min` and the `sweep.plot_*` functions. For a 2D
        grid, pick a row or column with `select_sweep` first.
    bulk_ions : list[str], optional
        Recorded in the manifest for the record only - this path never uses it.
    n_min : int
        Clusters below this many core+border atoms are dropped as they are
        found. Defaults to 1 (keep everything) because `optimise_n_min` sweeps
        n_min downstream over the sizes written here; raising it only throws
        information away.
    cached_shuffles : sequence[array], optional
        Explicit permutations of the atom indices, one per control dataset -
        the same parameter `cluster.run_relabel_sweep` takes. Passing the same
        ones to both makes the two paths run on *identical* control datasets,
        which is how their control runs are compared. Left out, each dataset
        draws its own from `relabel_seed + run_index`.
    n_jobs : int
        Control datasets to process in parallel (-1 = one per core). Uses
        forked processes that share the tip's positions copy-on-write; falls
        back to serial where fork is unavailable.
    resume : bool
        Skip cells whose parquet files already exist, so an interrupted sweep
        picks up where it stopped.

    Returns
    -------
    list[dict]
        One entry per (order, dmax) cell, also written to
        `{output_dir}/{prefix}_fast_grid_manifest.json`. Every entry points at
        the same parquet layout `run_relabel_sweep` writes, and carries
        `"pipeline": "core_only"` to record that its sizes are core+border
        counts with no envelope contribution.
    """

    if cached_shuffles is not None and len(cached_shuffles) < n_relabel_runs:
        raise ValueError(
            f"cached_shuffles has {len(cached_shuffles)} permutations but "
            f"n_relabel_runs={n_relabel_runs} datasets were asked for"
        )

    orders = sorted({int(o) for o in np.atleast_1d(orders)})
    dmax_values = sorted({float(d) for d in np.atleast_1d(dmax_values)})
    cells = [(order, dmax) for order in orders for dmax in dmax_values]
    core_ions = list(core_ions)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}
    for order, dmax in cells:
        run_name = _cell_run_name(prefix, order, dmax)
        cell_paths = _cell_paths(output_dir, run_name)
        cell_paths["run_name"] = run_name
        cell_paths["random_dir"].mkdir(parents=True, exist_ok=True)
        paths[(order, dmax)] = cell_paths

    positions = pos_df[["x", "y", "z"]].to_numpy()
    n_atoms = len(pos_df)
    max_dmax = max(dmax_values)

    vprint = print if verbose else lambda *a, **kw: None

    # ---------------------------------
    # The real tip
    # ---------------------------------
    started = perf_counter()
    full_code = core_code_from_mc(pos_df["mc"].to_numpy(), rrng_df, core_ions)
    real_index, real_code = split_core_code(full_code)
    vprint(
        f"{len(real_index):,} core-species atoms of {n_atoms:,} "
        f"({len(real_index) / n_atoms:.1%})" if n_atoms else "empty pos file"
    )
    if len(real_index) == 0:
        warnings.warn(
            f"no atoms ranged as any of {core_ions}; every cell will be empty"
        )

    missing_real = [
        cell
        for cell in cells
        if not (resume and paths[cell]["real_cluster_stats_path"].exists())
    ]

    if missing_real:
        neighbourhood = CoreNeighbourhood(
            positions[real_index], orders, max_dmax, workers=-1
        )
        vprint(
            f"real dataset: {neighbourhood.n_edges:,} candidate links "
            f"within {max_dmax} nm, precomputed in "
            f"{perf_counter() - started:.1f}s"
        )

        formula_map = build_ion_formula_map(rrng_df)

        for order, dmax in missing_real:
            labels = _apply_n_min(neighbourhood.labels(order, dmax), n_min)
            atoms = _atoms_frame(pos_df, real_index, labels, real_code, core_ions)

            stats = analyse_clusters(
                atoms,
                ion_formula_map=formula_map,
                atomic_volume=atomic_volume,
                core_species=core_ions,
            )
            stats.to_parquet(paths[(order, dmax)]["real_cluster_stats_path"], index=False)

            if save_real_atoms:
                atoms.to_parquet(
                    paths[(order, dmax)]["real_clustered_atoms_path"], index=False
                )

            vprint(
                f"  order={order} dmax={dmax}: {len(stats)} clusters, "
                f"{int((labels != -1).sum())} core+border atoms"
            )

        del neighbourhood

    vprint(f"real dataset done in {perf_counter() - started:.1f}s")

    # ---------------------------------
    # The randomised control datasets
    # ---------------------------------
    todo = []
    for run_index in range(n_relabel_runs):
        if resume and all(
            (paths[cell]["random_dir"] / f"run_{run_index:03d}.parquet").exists()
            for cell in cells
        ):
            continue
        todo.append(run_index)

    if todo:
        started = perf_counter()
        n_workers = _resolve_n_jobs(n_jobs, len(todo))

        _SHARED.clear()
        _SHARED.update(
            positions=positions,
            core_ions=core_ions,
            # one code per atom of the whole tip, not just the core ones:
            # a control dataset is a shuffle of it
            full_code=full_code,
            orders=orders,
            cells=cells,
            max_dmax=max_dmax,
            n_min=n_min,
            relabel_seed=relabel_seed,
            cached_shuffles=cached_shuffles,
            resume=resume,
            random_dirs={cell: str(paths[cell]["random_dir"]) for cell in cells},
            workers=1 if n_workers > 1 else -1,
        )

        try:
            if n_workers > 1:
                # fork, so each worker inherits `positions` copy-on-write
                # instead of being sent its own copy
                with ProcessPoolExecutor(max_workers=n_workers) as pool:
                    for run_index, elapsed in pool.map(_random_dataset_job, todo):
                        vprint(f"  control run {run_index}: {elapsed:.1f}s")
            else:
                for run_index in todo:
                    _, elapsed = _random_dataset_job(run_index)
                    vprint(f"  control run {run_index}: {elapsed:.1f}s")
        finally:
            # never leave the tip's positions alive in module state
            _SHARED.clear()

        vprint(
            f"{len(todo)} control datasets x {len(cells)} cells in "
            f"{perf_counter() - started:.1f}s"
        )

    # ---------------------------------
    # Manifest
    # ---------------------------------
    # A grid with one value on an axis IS a 1D sweep, so it is labelled as one
    # and every sweep.* consumer takes it directly, no select_sweep needed.
    if len(dmax_values) == 1 and len(orders) > 1:
        swept_param = "order"
    elif len(orders) == 1:
        swept_param = "dmax"
    else:
        swept_param = None

    manifest = []
    for order, dmax in cells:
        entry = {
            "run_name": paths[(order, dmax)]["run_name"],
            "dmax": dmax,
            "order": order,
            "core_species": core_ions,
            "envelope_species": list(bulk_ions) if bulk_ions else [],
            "n_min": n_min,
            "real_cluster_stats_path": str(
                paths[(order, dmax)]["real_cluster_stats_path"]
            ),
            "random_dir": str(paths[(order, dmax)]["random_dir"]),
            "n_relabel_runs": n_relabel_runs,
            "pipeline": "core_only",
        }
        if save_real_atoms:
            entry["real_clustered_atoms_path"] = str(
                paths[(order, dmax)]["real_clustered_atoms_path"]
            )
        if swept_param is not None:
            entry["swept_param"] = swept_param
            entry["swept_value"] = order if swept_param == "order" else dmax
        manifest.append(entry)

    # named per swept axis when it is a 1D sweep, so an order sweep and a dmax
    # sweep with the same prefix do not overwrite each other's manifest
    manifest_path = output_dir / (
        f"{prefix}_{swept_param}_sweep_manifest.json"
        if swept_param
        else f"{prefix}_fast_grid_manifest.json"
    )
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    return manifest


def run_fast_parameter_sweep(
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
    n_jobs: int = 1,
    verbose: bool = False,
) -> list:
    """Drop-in replacement for `sweep.run_parameter_sweep`, on the fast path.

    Same arguments, same returned manifest shape, same files on disk - so an
    existing notebook only changes the name of this one call and everything
    downstream (`optimise_n_min`, `plot_real_clusters_across_swept_param`,
    `plot_real_cluster_ratio_across_swept_param`, the composition plots) keeps
    working untouched. Use it when you are sweeping one axis at a time;
    `run_fast_grid_sweep` is the same machinery with both axes at once.

    What differs, and it is worth knowing rather than discovering:

    - **No envelope or erosion.** Clusters are core+border atoms, which is
      exactly what this stage reads (`size_metric="core"`, posgen's
      `Cluster Size(core)`) - see the module docstring for why neither stage can
      change that number. `envelope_radius`, `d_erosion` and `dclassify` are
      therefore not used; passing them warns rather than pretending.
    - `size` in the per-cluster stats counts core+border atoms, not every atom,
      and only core-species ions get a column. `optimise_n_min` and the two
      real-vs-random graphs are unaffected (they threshold `size_core`). The
      *absolute* composition graph changes meaning - its denominator becomes the
      core ions rather than every atom in the cluster.
    - `n_min` defaults to 2 as in `run_parameter_sweep`. Below order 2 a cluster
      cannot have one atom anyway, so this only matters at order 1; `n_min=1`
      keeps more for `optimise_n_min` to work with.
    - `n_jobs` is the one argument `run_parameter_sweep` does not have. It is
      worth adding to the call (`n_jobs=-1`) - roughly halves the wall time on
      8 cores - but it is left at 1 by default because it forks, and CPython
      warns about forking a process that already has threads running.

    Run the full `cluster.run_relabel_sweep` once at the chosen parameters, as
    the example notebook does - composition, volume and volume fraction all
    need the envelope.
    """

    if swept_param not in ("order", "dmax"):
        raise ValueError('swept_param must be "order" or "dmax"')
    if swept_param == "order" and fixed_dmax is None:
        raise ValueError("fixed_dmax is required when sweeping order")
    if swept_param == "dmax" and fixed_order is None:
        raise ValueError("fixed_order is required when sweeping dmax")

    unused = [
        name
        for name, value in (
            ("envelope_radius", envelope_radius),
            ("d_erosion", d_erosion),
            ("dclassify", dclassify),
        )
        if value is not None
    ]
    if unused:
        warnings.warn(
            f"{', '.join(unused)} ignored: this path runs the DBSCAN/kNN stage "
            "only, with no envelope or erosion. That does not change what the "
            "sweep graphs read (core+border counts), but the per-cluster `size` "
            "and non-core ion counts it writes are not comparable with "
            "run_parameter_sweep's. Use cluster.run_relabel_sweep for the final "
            "run at the chosen parameters."
        )

    if swept_param == "order":
        orders, dmax_values = swept_values, [fixed_dmax]
    else:
        orders, dmax_values = [fixed_order], swept_values

    return run_fast_grid_sweep(
        pos_df, rrng_df, core_ions,
        orders=orders, dmax_values=dmax_values,
        bulk_ions=bulk_ions, n_min=n_min,
        n_relabel_runs=n_relabel_runs, output_dir=output_dir, prefix=prefix,
        relabel_seed=relabel_seed, cached_shuffles=cached_shuffles,
        n_jobs=n_jobs, verbose=verbose,
    )


def _resolve_n_jobs(n_jobs, n_tasks):
    if n_jobs in (0, 1):
        return 1
    try:
        available = mp.cpu_count()
    except NotImplementedError:
        available = 1
    requested = available if n_jobs < 0 else min(n_jobs, available)
    n_workers = max(1, min(requested, n_tasks))

    if n_workers > 1 and mp.get_start_method(allow_none=False) != "fork":
        warnings.warn(
            "n_jobs > 1 needs the 'fork' start method to share the tip's "
            "positions with the workers; running serially instead"
        )
        return 1
    return n_workers


# ===========================================================================
# Consuming a grid manifest


def select_sweep(manifest: list, order=None, dmax=None) -> list:
    """Take one row or column out of a grid manifest as a 1D sweep manifest.

    Fixing `order` gives a dmax sweep, fixing `dmax` gives an order sweep -
    each in the shape `sweep.optimise_n_min` and the `sweep.plot_*` functions
    expect.

    >>> dmax_sweep = select_sweep(grid, order=7)
    >>> pgp.optimise_n_min(dmax_sweep, n_min_values=range(1, 200))
    """

    if (order is None) == (dmax is None):
        raise ValueError("fix exactly one of order= or dmax=")

    if order is not None:
        selected = [e for e in manifest if e["order"] == int(order)]
        swept_param = "dmax"
    else:
        selected = [e for e in manifest if np.isclose(e["dmax"], float(dmax))]
        swept_param = "order"

    if not selected:
        available_orders = sorted({e["order"] for e in manifest})
        available_dmax = sorted({e["dmax"] for e in manifest})
        raise ValueError(
            f"nothing in this grid at order={order}, dmax={dmax}; "
            f"it covers orders {available_orders} and dmax {available_dmax}"
        )

    out = [dict(e, swept_param=swept_param, swept_value=e[swept_param]) for e in selected]
    return sorted(out, key=lambda e: e["swept_value"])


def summarise_grid(
    manifest: list,
    n_min_values=range(1, 200),
    confidence_threshold: float = 0.95,
    size_metric: str = "core",
    exclude_ions: list = None,
    core_species: list = None,
) -> pd.DataFrame:
    """One row per (order, dmax) cell: the numbers you pick a parameter pair on.

    Reading the grid as two 1D slices (`select_sweep`) answers "which dmax, at
    this order" and "which order, at this dmax" separately. This answers them
    together, which is the actual question when both are free.

    Columns: `order`, `dmax`, `optimal_n_min` (smallest n_min whose true-positive
    ratio reaches `confidence_threshold`), `n_real_clusters` and
    `n_random_clusters` at that n_min, and the `true_positive_ratio` it reached.
    `n_min_reached` is False where no n_min in `n_min_values` got there - the
    cell is then reported at the largest n_min tried, and should not be chosen.

    The threshold rule itself is `sweep.optimise_n_min`, called once per cell,
    so there is exactly one definition of it in the package.

    >>> summary = summarise_grid(grid)
    >>> summary.pivot(index="order", columns="dmax", values="optimal_n_min")
    """

    n_min_values = list(n_min_values)
    rows = []

    for entry in manifest:
        cell = dict(entry, swept_param="dmax", swept_value=entry["dmax"])

        n_min, real, random = optimise_n_min(
            [cell],
            n_min_values=n_min_values,
            confidence_threshold=confidence_threshold,
            size_metric=size_metric,
            exclude_ions=exclude_ions,
            core_species=core_species,
            verbose=False,
        )

        rows.append(
            {
                "order": entry["order"],
                "dmax": entry["dmax"],
                "optimal_n_min": n_min[0],
                "n_real_clusters": real[0],
                "n_random_clusters": random[0],
                "true_positive_ratio": (
                    (real[0] - random[0]) / real[0] if real[0] else np.nan
                ),
                "n_min_reached": n_min[0] < max(n_min_values),
                "run_name": entry["run_name"],
            }
        )

    return pd.DataFrame(rows).sort_values(["order", "dmax"], ignore_index=True)


def load_sweep_size_table(
    manifest: list, include_random: bool = True, size_metric: str = "size_core"
) -> pd.DataFrame:
    """Every cluster of every cell as one tidy row, for distribution plots.

    Columns: `order`, `dmax`, `run_type` ("real"/"random"), `run` ("real" or
    "random_007"), `cluster_id`, `size`, `size_core`, the per-ion counts, and -
    for real runs - the centroid/shape columns `analyse_clusters` produced. That
    is what lets a point in a size distribution carry its own cluster id, so
    the cluster can be looked up in `viz.show_clusters_in_pos`.

    Random rows are cheap to include (one row per random cluster) but there are
    `n_relabel_runs` times as many of them; pass `include_random=False` for
    real clusters only.
    """

    frames = []

    for entry in manifest:

        cell_frames = [(pd.read_parquet(entry["real_cluster_stats_path"]), "real")]

        if include_random:
            cell_frames += [
                (pd.read_parquet(path), f"random_{path.stem.split('_')[-1]}")
                for path in sorted(Path(entry["random_dir"]).glob("run_*.parquet"))
            ]

        for frame, run in cell_frames:
            frame.insert(0, "order", entry["order"])
            frame.insert(1, "dmax", entry["dmax"])
            frame["run"] = run
            frame["run_type"] = "real" if run == "real" else "random"
            frames.append(frame)

    if not frames:
        return pd.DataFrame()

    table = pd.concat(frames, ignore_index=True)

    if size_metric not in table.columns:
        raise ValueError(
            f"{size_metric!r} is not a column of these results "
            f"({list(table.columns)})"
        )

    return table


# ===========================================================================
# Validation against the full pipeline


def validate_fast_grid(
    pos_df,
    rrng_df,
    core_ions: list,
    bulk_ions: list,
    orders,
    dmax_values,
    n_min: int = 1,
    n_min_values=range(1, 51),
    relative_count_floor: int = 10,
    verbose: bool = True,
) -> pd.DataFrame:
    """Check the fast grid against `cluster.run_clustering`, cell by cell.

    Runs both on the real tip and compares what the sweep stage actually reads.
    Three things are reported per cell:

    `n_differing_atoms` vs `n_ambiguous_border_atoms`
        Core-species atoms that ended up in a different cluster, against the
        number that *any* correct implementation is free to place either way
        (see `CoreNeighbourhood.ambiguous_border_atoms`). The first should
        always be no larger than the second - that is the pass condition, and
        `differences_within_ambiguity` records it. Anything else is a bug.

    `max_cluster_count_difference` / `max_relative_count_difference`
        The largest disagreement, over `n_min_values`, in the number of
        clusters of at least n_min core+border atoms. This is the number
        `optimise_n_min` and the real-vs-random plots are built from, so it is
        the one that decides whether a sweep run through this path would have
        chosen differently. The relative form is taken only over n_min values
        where the full pipeline still finds at least `relative_count_floor`
        clusters - past that the curve is down to a handful of clusters and a
        single border atom reads as tens of percent.

    `speedup`
        Wall-clock ratio for the cluster search itself.

    Returns one row per cell; run it whenever this code or the engine changes,
    or to confirm the shortcut holds on a new tip.
    """
    from .cluster import run_clustering
    from .io import assign_rrng_to_pos

    orders = sorted({int(o) for o in np.atleast_1d(orders)})
    dmax_values = sorted({float(d) for d in np.atleast_1d(dmax_values)})

    labelled = assign_rrng_to_pos(pos_df, rrng_df)
    real_index, _ = split_core_code(
        core_code_from_mc(pos_df["mc"].to_numpy(), rrng_df, list(core_ions))
    )
    neighbourhood = CoreNeighbourhood(
        pos_df[["x", "y", "z"]].to_numpy()[real_index], orders, max(dmax_values)
    )

    rows = []

    for order in orders:
        for dmax in dmax_values:

            started = perf_counter()
            full = run_clustering(
                labelled, core_ions, bulk_ions, dmax, order, n_min,
                envelope_radius=dmax / 2, d_erosion=dmax / 2, verbose=False,
            )
            full_seconds = perf_counter() - started

            started = perf_counter()
            fast_labels = _apply_n_min(neighbourhood.labels(order, dmax), n_min)
            fast_seconds = perf_counter() - started

            # the full pipeline's core+border atoms, in pos_df row order
            full_labels = pd.Series(-1, index=pos_df.index, dtype=np.int64)
            core_rows = full.loc[full["ion"].isin(core_ions)]
            full_labels.loc[core_rows.index] = core_rows["cluster_id"].to_numpy()
            full_core = full_labels.to_numpy()[real_index]

            agreement = _partition_agreement(full_core, fast_labels)
            ambiguous = neighbourhood.ambiguous_border_atoms(order, dmax)

            # what the sweep stage actually reads off these runs
            full_counts = _counts_at_least(full_core, n_min_values)
            count_difference = np.abs(
                full_counts - _counts_at_least(fast_labels, n_min_values)
            )
            # only where the curve still has enough clusters for a ratio to
            # mean anything - at the tail, one atom moving reads as tens of %
            meaningful = full_counts >= relative_count_floor
            relative_difference = np.divide(
                count_difference,
                full_counts,
                out=np.zeros(len(full_counts), dtype=float),
                where=meaningful,
            )

            rows.append(
                {
                    "order": order,
                    "dmax": dmax,
                    **agreement,
                    "n_ambiguous_border_atoms": ambiguous,
                    "differences_within_ambiguity": bool(
                        agreement["n_differing_atoms"] <= ambiguous
                    ),
                    "max_cluster_count_difference": int(count_difference.max()),
                    "max_relative_count_difference": float(relative_difference.max()),
                    "full_seconds": full_seconds,
                    "fast_seconds": fast_seconds,
                    "speedup": full_seconds / fast_seconds if fast_seconds else np.nan,
                }
            )

            if verbose:
                row = rows[-1]
                print(
                    f"order={order} dmax={dmax}: "
                    f"clusters {row['n_clusters_full']} vs {row['n_clusters_fast']}, "
                    f"atoms {row['n_atoms_full']} vs {row['n_atoms_fast']}, "
                    f"{row['n_differing_atoms']} atoms placed differently of "
                    f"{ambiguous} ambiguous, cluster counts differ by at most "
                    f"{row['max_cluster_count_difference']} "
                    f"({row['max_relative_count_difference']:.1%}), "
                    f"{row['full_seconds']:.2f}s vs {row['fast_seconds']:.3f}s "
                    f"({row['speedup']:.0f}x)"
                )

    return pd.DataFrame(rows)


def _counts_at_least(labels, n_min_values):
    """Clusters of at least n_min atoms, for each n_min - the sweep's own axis."""
    sizes = _sorted_sizes(labels)
    return np.array([int((sizes >= n).sum()) for n in n_min_values])


def _partition_agreement(labels_a, labels_b):
    """Compare two clusterings of the same atoms, ignoring cluster numbering.

    Cluster ids are arbitrary (this engine numbers from 0, posgen from 1, and
    connected-components numbers by traversal), so agreement is measured by
    matching each cluster of A to the cluster of B it shares most atoms with
    and counting the atoms that fall outside those matches.
    """

    clustered_a = labels_a != -1
    clustered_b = labels_b != -1

    result = {
        "n_clusters_full": int(len(np.unique(labels_a[clustered_a]))),
        "n_clusters_fast": int(len(np.unique(labels_b[clustered_b]))),
        "n_atoms_full": int(clustered_a.sum()),
        "n_atoms_fast": int(clustered_b.sum()),
    }

    # atoms clustered by one and not the other are always a difference
    disagree = int((clustered_a != clustered_b).sum())

    both = clustered_a & clustered_b
    if both.any():
        # sparse contingency table - dense would be n_clusters^2 cells
        overlap = coo_matrix(
            (
                np.ones(int(both.sum()), dtype=np.int64),
                (labels_a[both], labels_b[both]),
            )
        ).tocsr()
        # atoms outside each A-cluster's best-matching B-cluster
        best_match = np.asarray(overlap.max(axis=1).todense()).ravel()
        disagree += int(overlap.sum() - best_match.sum())

    result["n_differing_atoms"] = disagree
    return result


def _sorted_sizes(labels):
    """The clustering's size distribution: cluster sizes, ascending."""
    clustered = labels[labels != -1]
    if len(clustered) == 0:
        return np.zeros(0, dtype=np.int64)
    _, sizes = np.unique(clustered, return_counts=True)
    return np.sort(sizes)
