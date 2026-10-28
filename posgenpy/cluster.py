"""Pure-Python max-separation cluster search (DBSCAN + envelope + erosion),
plus the relabelling/randomisation machinery used to estimate a real-vs-random
cluster threshold without holding many full per-atom datasets in memory.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.neighbors import KDTree
from tqdm import tqdm

from .io import assign_rrng_to_pos, build_ion_formula_map

# ===========================================================================
# Cluster Search


def build_position_tree(pos_df):
    """Build one KDTree over ALL atom positions in pos_df.

    Positions never change across relabel runs (only the mc/species labels
    do), so this tree can be built once per pos file and reused across every
    swept parameter value and every relabel run, instead of being rebuilt on
    every call.
    """
    positions = pos_df[["x", "y", "z"]].to_numpy()
    return KDTree(positions), positions


def _resolve_order(order, k):
    """Accept either `order` (posgen's knn) or the deprecated `k`.

    `order` counts the point itself, matching posgen's `<dclassify knn>` and
    the `_order_N` in its filenames. `k` counted only the other neighbours,
    which made a run labelled "order 7" need `k=6` - the off-by-one that made
    sweeps line up against the wrong posgen run.
    """
    if order is None and k is None:
        raise ValueError("order is required (posgen's knn / Order)")

    if k is not None:
        if order is not None:
            raise ValueError("pass either order or k, not both")
        warnings.warn(
            "k is deprecated; use order instead. posgen's Order counts the "
            f"point itself, so k={k} is order={k + 1}.",
            DeprecationWarning,
            stacklevel=3,
        )
        return k + 1

    return order


def _summarise_sizes(clustered_df, core_species):
    """Per-cluster size table for a relabelled control run.

    Deliberately cheaper than `analyse_clusters` - no centroids, shape or
    PCA - because random runs are only ever consumed as size distributions.
    Records `size` (every atom), `size_core` (core+border atoms, posgen's
    `Cluster Size(core)`) and a column per ion, so a size threshold can later
    be applied on any of those without re-running the sweep.
    """
    cl = clustered_df.loc[clustered_df["cluster_id"] != -1]

    if cl.empty:
        return pd.DataFrame(columns=["cluster_id", "size", "size_core"])

    counts = pd.crosstab(cl["cluster_id"], cl["ion"])
    core_cols = [c for c in counts.columns if c in set(core_species)]

    out = pd.DataFrame(
        {
            "cluster_id": counts.index,
            "size": counts.sum(axis=1).to_numpy(),
            "size_core": counts[core_cols].sum(axis=1).to_numpy() if core_cols else 0,
        }
    )
    for col in counts.columns:
        out[col] = counts[col].to_numpy()

    return out


def _shrink(radius):
    """Largest float strictly below `radius`.

    posgen tests every distance with a strict `<` on squared float32 values;
    sklearn's `eps` / `query_radius` are inclusive (`<=`). Shrinking the radius
    by one ULP reproduces the strict test everywhere except for points sitting
    at exactly the boundary distance, which real reconstructions essentially
    never produce.
    """
    return float(np.nextafter(radius, 0.0))


def run_clustering(
    pos_df: pd.DataFrame,
    core_species: list,
    envelope_species: list,
    Dmax: float,
    order: int = None,
    N_min: int = 2,
    envelope_radius: float = None,
    d_erosion: float = None,
    dclassify: float = None,
    k: int = None,
    envelope_seeds: str = "core",
    erosion_reference: str = "bulk_tree",
    include_unranged: bool = False,
    n_min_count: str = "core",
    erosion_diagnostics: bool = False,
    tree: KDTree = None,
    tree_positions: np.ndarray = None,
    verbose: bool = True,
):
    """
    Maximum-separation cluster search: core classification, envelope, erosion.

    Reproduces posgen's `maxsep` algorithm (`posgen/operators/clusterOps.cpp`,
    `ClusterOps::runLinkClustering`). The defaults are posgen-faithful; the
    alternative settings preserve this engine's earlier behaviour so the two
    can be compared (see `posgenpy.parity`).

    Parameters
    ----------
    pos_df : pandas.DataFrame
        Must contain x, y, z and ion. Unranged atoms (ion is NA) must be kept
        in the frame - with `include_unranged=True` they take part in the
        envelope and erosion steps, exactly as posgen's `<unranged>` element
        makes them part of its bulk tree.

    core_species, envelope_species : list[str]
        posgen's `<core><typelist>` and `<bulk><typelist>`. Note posgen treats
        each ranged ion label as a distinct type, so a molecular label such as
        "C2" is NOT covered by listing "C".

    Dmax : float
        Maximum separation for linking core points (nm).

    order : int
        posgen's `knn` / "Order" - the same number that appears in
        `<dclassify knn="...">` and in `..._dmax0.5_order_N` filenames. A
        point is a core point when `order` points, counting itself, lie within
        `dclassify`, so this maps straight onto sklearn's `min_samples`.

        Verified against real posgen runs: order 7 on MNSP tip C gives
        216 clusters / 26,728 clustered atoms, and order 6 on
        MNSP tip A _top gives 1,244 clusters / 28,387 core+border atoms,
        both matching those runs atom for atom.

    k : int, optional
        Deprecated spelling of the same thing counting only the OTHER
        neighbours, i.e. `order = k + 1`. Kept so existing calls keep their
        original meaning instead of silently shifting by one; pass `order`
        instead.

    N_min : int
        posgen's `<sizeclip nmin>`. See `n_min_count` for what is counted.

    envelope_radius : float
        posgen's `dbulk` (nm).

    d_erosion : float
        posgen's `derode` (nm).

    dclassify : float, optional
        posgen's separate core-classification distance. Every run this engine
        has been validated against uses either 0 (disabled) or exactly Dmax,
        both of which collapse into the single-radius classification used
        here; any other value is rejected rather than silently ignored.

    envelope_seeds : {"core", "all"}
        Which cluster members seed the envelope. posgen grows it from
        classified core points only - border points (core-species atoms that
        failed the knn test but got attached to a cluster) are merged into the
        cluster only after envelope and erosion have run. "all" seeds from
        every cluster member, which grows a larger envelope.

    erosion_reference : {"bulk_tree", "envelope_species"}
        What counts as the non-clustered material that erodes the envelope.
        posgen measures against untagged points of its bulk tree: bulk-species
        atoms plus, when `include_unranged`, unranged and other non-core
        atoms. Unclustered CORE-species atoms are never included - they live
        in posgen's core/border trees, which erosion never queries.
        "envelope_species" restricts it to unclustered envelope-species atoms.

    include_unranged : bool
        Mirrors posgen's `<unranged>`: puts unranged atoms and ranged atoms in
        neither species list into the bulk tree, making them both eligible for
        capture by the envelope and part of the erosion reference.

        Defaults to False because that is what posgen actually did on the runs
        this engine is validated against, even though their XML carries
        `<unranged foroutput="True" forstats="True"/>` - the capitalised
        "True" that posgenpy's XML writer emitted evidently did not parse as
        true. Empirically: on MNSP tip A _top, False reproduces posgen's
        45,194 clustered atoms exactly, while True gives 44,426.

    n_min_count : {"core", "all"}
        posgen counts cluster size as core+border atoms (bulk only with
        `<countbulk/>`, unranged never) and applies the clip BEFORE the
        envelope, so a discarded cluster's surrounding bulk atoms stay
        untagged and still erode neighbouring clusters. "all" counts every
        atom in the cluster and clips at the end instead.

    erosion_diagnostics : bool
        Also record the exact distance from each enveloped atom to the nearest
        non-clustered atom in an `erosion_distance` column. Requires an extra
        KDTree over the reference set, so it is off by default.

    tree, tree_positions : optional
        A prebuilt KDTree over ALL atoms of `pos_df` (see
        `build_position_tree`) and the positions array it was built from.
        Lets one tree be reused across many calls on the same pos file.

    Returns
    -------
    pandas.DataFrame
        Ranged atoms only, with a `cluster_id` column (-1 = unclustered) and,
        if requested, `erosion_distance`. Unranged atoms take part in the
        algorithm but are excluded from the result, matching posgen's
        `splitUnrangedClustered`; the count of unranged atoms that ended up
        inside a cluster is reported in `df.attrs["n_unranged_clustered"]`.
    """

    vprint = print if verbose else lambda *args, **kwargs: None

    order = _resolve_order(order, k)

    if envelope_radius is None or d_erosion is None:
        raise ValueError("envelope_radius and d_erosion are required")

    if dclassify not in (None, 0) and not np.isclose(dclassify, Dmax):
        raise ValueError(
            "dclassify distinct from Dmax/0 is not supported by this engine; "
            "posgen's dclassify is always folded into Dmax here (pass "
            "dclassify=None, 0, or Dmax)."
        )
    if envelope_seeds not in ("core", "all"):
        raise ValueError('envelope_seeds must be "core" or "all"')
    if erosion_reference not in ("bulk_tree", "envelope_species"):
        raise ValueError(
            'erosion_reference must be "bulk_tree" or "envelope_species"'
        )
    if n_min_count not in ("core", "all"):
        raise ValueError('n_min_count must be "core" or "all"')

    # ---------------------------------
    # Species partition, over ALL atoms
    # ---------------------------------
    # posgen splits every ion into one of three KD trees; the bulk tree is the
    # only one the envelope and erosion steps ever query.
    n_all = len(pos_df)

    ion_arr = pos_df["ion"].to_numpy()
    is_ranged = pd.notna(ion_arr)
    ion_str = np.where(is_ranged, ion_arr.astype(str), "")

    is_core_species = np.isin(ion_str, list(core_species)) & is_ranged
    is_env_species = np.isin(ion_str, list(envelope_species)) & is_ranged

    in_bulk_tree = is_env_species.copy()
    if include_unranged:
        in_bulk_tree |= ~is_core_species & ~is_env_species

    if tree is not None:
        if tree_positions is None:
            raise ValueError("tree_positions must be supplied alongside tree")
        all_positions = tree_positions
    else:
        all_positions = pos_df[["x", "y", "z"]].to_numpy()
        vprint("Initialising KDTree")
        tree = KDTree(all_positions)

    if len(all_positions) != n_all:
        raise ValueError("tree_positions does not match pos_df length")

    # cluster label per atom, in the global (all-atom) index space
    labels = np.full(n_all, -1, dtype=int)

    # ---------------------------------
    # Core classification + linking
    # ---------------------------------
    vprint("Classifying and linking core points")

    core_idx = np.where(is_core_species)[0]

    if len(core_idx) == 0:
        out = pos_df.loc[is_ranged].copy()
        out["cluster_id"] = -1
        out.attrs["n_unranged_clustered"] = 0
        return out

    db = DBSCAN(eps=_shrink(Dmax), min_samples=order)
    core_labels = db.fit_predict(all_positions[core_idx])
    labels[core_idx] = core_labels

    # classified core points; the rest of a cluster's core-species members are
    # posgen's "border" ions, which never seed the envelope
    is_core_point = np.zeros(n_all, dtype=bool)
    is_core_point[core_idx[db.core_sample_indices_]] = True

    cluster_ids = np.unique(labels[labels != -1])

    # ---------------------------------
    # Early size clip (posgen's default, pre-envelope)
    # ---------------------------------
    if n_min_count == "core" and N_min > 1 and len(cluster_ids) > 0:
        vprint("Applying N_min to core+border counts (pre-envelope)")
        sizes = np.bincount(labels[labels != -1], minlength=cluster_ids.max() + 1)
        too_small = np.isin(labels, np.where(sizes < N_min)[0]) & (labels != -1)
        labels[too_small] = -1
        cluster_ids = np.unique(labels[labels != -1])

    n_core_clustered = int((labels != -1).sum())
    vprint(f"Core+border atoms clustered: {n_core_clustered} "
           f"in {len(cluster_ids)} clusters")

    # ---------------------------------
    # Envelope
    # ---------------------------------
    vprint("Enveloping")

    labels_env = labels.copy()
    # posgen's bulk-tree "tag": set once a point is claimed by any cluster and
    # never cleared, including by erosion
    tagged = np.zeros(n_all, dtype=bool)

    seed_pool = is_core_point if envelope_seeds == "core" else (labels != -1)
    env_radius = _shrink(envelope_radius)

    for cid in tqdm(cluster_ids, disable=not verbose):

        seed_idx = np.where((labels == cid) & seed_pool)[0]

        if len(seed_idx) == 0:
            continue

        neighbour_lists = tree.query_radius(all_positions[seed_idx], r=env_radius)
        candidate_idx = np.unique(np.concatenate(neighbour_lists))

        # only bulk-tree points, and only ones no cluster has claimed yet
        candidate_idx = candidate_idx[
            in_bulk_tree[candidate_idx] & ~tagged[candidate_idx]
        ]

        labels_env[candidate_idx] = cid
        tagged[candidate_idx] = True

    enveloped_mask = tagged.copy()

    if verbose:
        print("\n===== ENVELOPE STATS =====")
        print("Atoms added by envelope:", int(enveloped_mask.sum()))
        print("\nSpecies added by envelope:")
        print(pd.Series(ion_str[enveloped_mask]).replace("", "<unranged>").value_counts())

    # ---------------------------------
    # Erosion
    # ---------------------------------
    vprint("Eroding enveloped ions")

    labels_final = labels_env.copy()
    erosion_distance = np.full(n_all, np.nan)

    env_idx = np.where(enveloped_mask)[0]

    if erosion_reference == "bulk_tree":
        is_reference = in_bulk_tree & ~tagged
    else:
        is_reference = is_env_species & ~tagged

    if len(env_idx) > 0 and is_reference.any():

        # A single pass against a fixed snapshot: eroded atoms stay "tagged"
        # and never rejoin the reference set, so the result is order
        # independent (posgen runs this loop with `omp parallel for`).
        neighbour_lists = tree.query_radius(
            all_positions[env_idx], r=_shrink(d_erosion)
        )

        remove_mask = np.fromiter(
            (bool(is_reference[nn].any()) for nn in neighbour_lists),
            dtype=bool,
            count=len(env_idx),
        )

        labels_final[env_idx[remove_mask]] = -1

        if erosion_diagnostics:
            ref_idx = np.where(is_reference)[0]
            ref_tree = KDTree(all_positions[ref_idx])
            dist, _ = ref_tree.query(all_positions[env_idx], k=1)
            erosion_distance[env_idx] = dist[:, 0]

    if verbose:
        removed = enveloped_mask & (labels_final == -1)
        print("\n===== EROSION STATS =====")
        print("Atoms removed by erosion:", int(removed.sum()))
        print("\nSpecies removed by erosion:")
        print(pd.Series(ion_str[removed]).replace("", "<unranged>").value_counts())
        print("Envelope ions surviving:",
              int((enveloped_mask & (labels_final != -1)).sum()))

    # ---------------------------------
    # Late size clip (whole-cluster counting)
    # ---------------------------------
    if n_min_count == "all" and N_min > 1:
        vprint("Applying N_min to whole-cluster counts")
        final_ids = np.unique(labels_final[labels_final != -1])
        if len(final_ids) > 0:
            sizes = np.bincount(
                labels_final[labels_final != -1], minlength=final_ids.max() + 1
            )
            too_small = (
                np.isin(labels_final, np.where(sizes < N_min)[0])
                & (labels_final != -1)
            )
            labels_final[too_small] = -1

    # ---------------------------------
    # Result: ranged atoms only (posgen's splitUnrangedClustered)
    # ---------------------------------
    n_unranged_clustered = int(((labels_final != -1) & ~is_ranged).sum())

    if verbose:
        print("\n===== BALANCE =====")
        print("Core+border atoms:", n_core_clustered)
        print("After envelope:", int((labels_env != -1).sum()))
        print("After erosion:", int((labels_final != -1).sum()))
        print("Unranged atoms inside clusters (excluded from output):",
              n_unranged_clustered)

    out = pos_df.loc[is_ranged].copy()
    out["cluster_id"] = labels_final[is_ranged]
    if erosion_diagnostics:
        out["erosion_distance"] = erosion_distance[is_ranged]
    out.attrs["n_unranged_clustered"] = n_unranged_clustered

    return out



def analyse_clusters(
    clustered_df,
    ion_formula_map: dict = None,
    atomic_volume=0.0118,  # nm³/atom
    core_species: list = None,
):
    """
    Analyse clusters and return one row per cluster.

    Parameters
    ----------
    clustered_df : pd.DataFrame
        Output of run_clustering().

    ion_formula_map : dict[str, dict[str, int]], optional
        Maps each ion label to its {element: count} formula (see
        io.build_ion_formula_map). When given, per-cluster atom counts are
        credited per-element (e.g. an "Fe2O3"-ranged atom credits 2 to Fe and
        3 to O) instead of one count per raw (possibly compound) ion label.

    atomic_volume : float
        Atomic volume (nm³/atom) used to estimate
        equivalent spherical radius.

    core_species : list[str], optional
        When given, a `size_core` column counts each cluster's core-species
        atoms (core + border). That is posgen's `Cluster Size(core)`, the
        quantity its `sizedist.txt` reports and the one the old
        real-cluster-count graphs were filtered on - as opposed to `size`,
        which counts every atom including the enveloped bulk.

    Returns
    -------
    pd.DataFrame
    """

    clusters = clustered_df[
        clustered_df["cluster_id"] != -1
    ]

    if ion_formula_map is None:
        species_list = sorted(
            clustered_df["ion"]
            .dropna()
            .unique()
        )
    else:
        species_list = sorted(
            {element for formula in ion_formula_map.values() for element in formula}
        )

    if len(clusters) == 0:
        return pd.DataFrame(columns=[
            "cluster_id", "x", "y", "z", "size", "size_core", "radius_of_gyration",
            "aspect_ratio", "principal_axis_1", "principal_axis_2",
            "principal_axis_3", "bbox_x", "bbox_y", "bbox_z", "bbox_volume",
            "cluster_density", "equivalent_spherical_radius", *species_list,
        ])

    results = []

    for cid, group in clusters.groupby("cluster_id"):

        coords = group[["x", "y", "z"]].to_numpy()

        n = len(group)

        # --------------------------
        # Centroid
        # --------------------------
        centre = coords.mean(axis=0)

        centred = coords - centre

        # --------------------------
        # Radius of gyration
        # --------------------------
        rg = np.sqrt(
            np.mean(
                np.sum(centred**2, axis=1)
            )
        )

        # --------------------------
        # Bounding box
        # --------------------------
        bbox_x = coords[:, 0].max() - coords[:, 0].min()
        bbox_y = coords[:, 1].max() - coords[:, 1].min()
        bbox_z = coords[:, 2].max() - coords[:, 2].min()

        bbox_volume = bbox_x * bbox_y * bbox_z

        cluster_density = (
            n / bbox_volume
            if bbox_volume > 0
            else np.nan
        )

        # --------------------------
        # PCA shape analysis
        # --------------------------
        if n >= 3:

            cov = np.cov(
                centred.T
            )

            eigvals = np.linalg.eigvalsh(cov)

            eigvals = np.maximum(
                eigvals,
                0
            )

            eigvals = np.sort(eigvals)[::-1]

            axis1 = 2.0 * np.sqrt(eigvals[0])
            axis2 = 2.0 * np.sqrt(eigvals[1])
            axis3 = 2.0 * np.sqrt(eigvals[2])

            if axis3 > 0:
                aspect_ratio = axis1 / axis3
            else:
                aspect_ratio = np.nan

        else:

            axis1 = np.nan
            axis2 = np.nan
            axis3 = np.nan
            aspect_ratio = np.nan

        # --------------------------
        # Equivalent spherical radius
        # --------------------------
        cluster_volume = n * atomic_volume

        equivalent_spherical_radius = (
            (3.0 * cluster_volume / (4.0 * np.pi))
            ** (1.0 / 3.0)
        )

        row = {
            "cluster_id": cid,
            "x": centre[0],
            "y": centre[1],
            "z": centre[2],
            "size": n,
            "size_core": (
                int(group["ion"].isin(core_species).sum())
                if core_species is not None else n
            ),
            "radius_of_gyration": rg,
            "aspect_ratio": aspect_ratio,
            "principal_axis_1": axis1,
            "principal_axis_2": axis2,
            "principal_axis_3": axis3,
            "bbox_x": bbox_x,
            "bbox_y": bbox_y,
            "bbox_z": bbox_z,
            "bbox_volume": bbox_volume,
            "cluster_density": cluster_density,
            "equivalent_spherical_radius": (
                equivalent_spherical_radius
            ),
        }

        # --------------------------
        # Ion counts (credited per-element when a formula map is given)
        # --------------------------
        counts = group["ion"].value_counts()

        if ion_formula_map is None:
            for ion in species_list:
                row[ion] = counts.get(ion, 0)
        else:
            element_counts = {}
            for ion, n_ions in counts.items():
                formula = ion_formula_map.get(ion, {ion: 1})
                for element, per_ion_count in formula.items():
                    element_counts[element] = (
                        element_counts.get(element, 0) + n_ions * per_ion_count
                    )
            for element in species_list:
                row[element] = element_counts.get(element, 0)

        results.append(row)

    return pd.DataFrame(results)


# ===========================================================================
# Relabelling / randomisation sweep


def run_relabel_sweep(
    pos_df: pd.DataFrame,
    rrng_df: pd.DataFrame,
    core_ions: list,
    bulk_ions: list,
    dmax: float,
    order: int = None,
    n_min: int = 2,
    envelope_radius: float = None,
    d_erosion: float = None,
    n_relabel_runs: int = 0,
    run_name: str = "run",
    output_dir=".",
    relabel_seed: int = 0,
    tree: KDTree = None,
    tree_positions: np.ndarray = None,
    cached_shuffles=None,
    dclassify: float = None,
    k: int = None,
    atomic_volume: float = 0.0118,
    verbose: bool = False,
) -> dict:
    """
    Run one real cluster search plus n_relabel_runs randomised (mc-shuffled)
    control searches, at a single (dmax, order) parameter combination.

    `order` is posgen's knn / Order, the number in its `_order_N` filenames.
    `k` is the deprecated spelling with `order = k + 1`.

    Memory discipline: only one relabel run's full per-atom labelled data is
    ever in memory at a time. The real run's full labelled atoms + per-cluster
    stats are persisted (matching posgen's old behaviour of only writing full
    pos-format output for the real run); each relabel run's per-atom data is
    discarded immediately after reducing it to a small per-cluster size table,
    which is what's persisted.

    Returns a small manifest dict (paths + counts + params only - never
    per-atom data), suitable for passing to sweep.optimise_n_min /
    sweep.plot_* and for JSON-serialising into a sweep manifest.
    """

    order = _resolve_order(order, k)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    random_dir = output_dir / f"{run_name}_random"
    random_dir.mkdir(parents=True, exist_ok=True)

    if tree is None:
        tree, tree_positions = build_position_tree(pos_df)

    formula_map = build_ion_formula_map(rrng_df)

    real_atoms_path = output_dir / f"{run_name}_real_clustered_atoms.parquet"
    real_stats_path = output_dir / f"{run_name}_real_cluster_stats.parquet"

    # ranged once here rather than only inside the branch below, because the
    # manifest records n_core_atoms whether or not the real run was cached -
    # it is the denominator sweep.percolation_report needs
    labelled = assign_rrng_to_pos(pos_df, rrng_df)
    n_core_atoms = int(labelled["ion"].isin(core_ions).sum())

    if not (real_atoms_path.exists() and real_stats_path.exists()):
        clustered = run_clustering(
            labelled, core_ions, bulk_ions, dmax, order, n_min,
            envelope_radius, d_erosion,
            dclassify=dclassify, tree=tree, tree_positions=tree_positions,
            verbose=verbose,
        )
        stats = analyse_clusters(
            clustered, ion_formula_map=formula_map, atomic_volume=atomic_volume,
            core_species=core_ions,
        )
        clustered.to_parquet(real_atoms_path, index=False)
        stats.to_parquet(real_stats_path, index=False)
        del clustered, stats

    del labelled

    mc_values = pos_df["mc"].to_numpy()
    random_paths = []

    for r in range(n_relabel_runs):

        run_path = random_dir / f"run_{r:03d}.parquet"
        random_paths.append(run_path)

        if run_path.exists():
            continue

        if cached_shuffles is not None:
            shuffled_mc = mc_values[cached_shuffles[r]]
        else:
            rng = np.random.default_rng(relabel_seed + r)
            shuffled_mc = rng.permutation(mc_values)

        shuffled_pos = pos_df.copy()
        shuffled_pos["mc"] = shuffled_mc

        labelled = assign_rrng_to_pos(shuffled_pos, rrng_df)
        clustered = run_clustering(
            labelled, core_ions, bulk_ions, dmax, order, n_min,
            envelope_radius, d_erosion,
            dclassify=dclassify, tree=tree, tree_positions=tree_positions,
            verbose=False,
        )

        sizes = _summarise_sizes(clustered, core_ions)
        sizes.to_parquet(run_path, index=False)

        del shuffled_pos, labelled, clustered, sizes

    return {
        "run_name": run_name,
        "dmax": dmax,
        "order": order,
        "core_species": list(core_ions),
        "envelope_species": list(bulk_ions),
        "n_min": n_min,
        "real_clustered_atoms_path": str(real_atoms_path),
        "real_cluster_stats_path": str(real_stats_path),
        "random_dir": str(random_dir),
        "n_relabel_runs": n_relabel_runs,
        "n_core_atoms": n_core_atoms,
    }
