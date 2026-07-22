"""Pure-Python max-separation cluster search (DBSCAN + envelope + erosion),
plus the relabelling/randomisation machinery used to estimate a real-vs-random
cluster threshold without holding many full per-atom datasets in memory.
"""
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


def run_clustering(
    pos_df: pd.DataFrame,
    core_species: list,
    envelope_species: list,
    Dmax: float,
    k: int,
    N_min: int,
    envelope_radius: float,
    d_erosion: float,
    dclassify: float = None,
    tree: KDTree = None,
    tree_positions: np.ndarray = None,
    verbose: bool = True,
):
    """
    DBSCAN + Envelope + Erosion clustering for APT data.

    Parameters
    ----------
    pos_df : pandas.DataFrame
        Must contain: x, y, z, ion

    core_species : list[str]
        Species used for DBSCAN clustering.

    envelope_species : list[str]
        Species allowed to be added during envelope.

    Dmax : float
        DBSCAN search radius (nm).

    k : int
        Number of neighbours for DBSCAN ("Order").

    N_min : int
        Minimum allowable cluster size at this stage (a permissive baseline
        such as 2 is typical - the finer n_min sweep used to pick a final
        threshold is applied later, post-hoc, against the saved per-cluster
        sizes; it does not need to match this value).

    envelope_radius : float
        Radius used for envelope step (nm).

    d_erosion : float
        Erosion distance (nm).

    dclassify : float, optional
        posgen's separate core-classification distance. In every real usage
        pattern this engine has been validated against, dclassify is always
        either 0 (disabled) or exactly equal to Dmax, both of which collapse
        into the same single-eps DBSCAN classification used here. Any other
        value would require a genuinely different two-phase implementation
        this engine does not (yet) provide, so it's rejected explicitly
        rather than silently ignored.

    tree, tree_positions : optional
        A pre-built KDTree (from build_position_tree) over ALL atoms in
        pos_df, plus the exact positions array it was built from. Passing
        these in lets a caller reuse one KDTree across many calls (e.g. many
        relabel runs on the same pos file). If omitted, a tree scoped to just
        the ranged subset is built here from scratch (matching a one-off,
        single-call usage).

    Returns
    -------
    pandas.DataFrame
        Ranged ions only, with cluster_id column.
    """

    vprint = print if verbose else lambda *args, **kwargs: None

    if dclassify not in (None, 0) and not np.isclose(dclassify, Dmax):
        raise ValueError(
            "dclassify distinct from Dmax/0 is not supported by this engine; "
            "posgen's dclassify is always folded into Dmax here (pass "
            "dclassify=None, 0, or Dmax)."
        )

    # ---------------------------------
    # Keep only ranged ions
    # ---------------------------------
    vprint("Filtering POS file")
    ranged_mask = pos_df["ion"].notna().to_numpy()
    df = pos_df[ranged_mask].copy()

    positions = df[["x", "y", "z"]].to_numpy()
    species = df["ion"].to_numpy()

    n_atoms = len(df)

    labels = np.full(n_atoms, -1, dtype=int)

    # ---------------------------------
    # DBSCAN on core species
    # ---------------------------------
    vprint("Running DBSCAN on core species")
    core_mask = np.isin(species, core_species)

    core_positions = positions[core_mask]
    core_indices = np.where(core_mask)[0]

    if len(core_positions) == 0:
        df["cluster_id"] = labels
        return df

    db = DBSCAN(
        eps=Dmax,
        min_samples=k + 1
    )

    core_labels = db.fit_predict(core_positions)

    labels[core_indices] = core_labels

    cluster_ids = np.unique(labels)
    cluster_ids = cluster_ids[cluster_ids != -1]

    # ---------------------------------
    # Neighbour search setup
    # ---------------------------------
    vprint("Initialising KDTree")

    if tree is not None:
        if tree_positions is None:
            raise ValueError("tree_positions must be supplied alongside tree")

        # translate indices into the full (all-atom) tree space to indices
        # into this run's ranged-only local arrays
        global_to_local = np.full(len(tree_positions), -1, dtype=int)
        global_to_local[ranged_mask] = np.arange(n_atoms)
        query_tree = tree
    else:
        query_tree = KDTree(positions)
        global_to_local = None

    # ---------------------------------
    # Envelope
    # ---------------------------------

    vprint("Enveloping")

    labels_env = labels.copy()

    # Tracks ions added during envelope
    enveloped_mask = np.zeros(n_atoms, dtype=bool)

    # only these species are eligible
    bulk_mask = np.isin(species, envelope_species)

    for cid in tqdm(cluster_ids, disable=not(verbose)):

        cluster_idx = np.where(labels == cid)[0]

        if len(cluster_idx) == 0:
            continue

        neighbour_lists = query_tree.query_radius(
            positions[cluster_idx],
            r=envelope_radius
        )

        candidate_idx = np.unique(
            np.concatenate(neighbour_lists)
        )

        if global_to_local is not None:
            candidate_idx = global_to_local[candidate_idx]
            candidate_idx = candidate_idx[candidate_idx != -1]

        # only bulk species
        candidate_idx = candidate_idx[
            bulk_mask[candidate_idx]
        ]

        # IMPORTANT:
        # mimic C++ tagging behaviour
        # only assign previously-unassigned ions
        new_idx = candidate_idx[
            labels_env[candidate_idx] == -1
        ]

        labels_env[new_idx] = cid

        enveloped_mask[new_idx] = True

    # ---------------------------------
    # Envelope statistics
    # ---------------------------------

    env_added_mask = (
        (labels == -1) &
        (labels_env != -1)
    )

    if verbose:
        print("\n===== ENVELOPE STATS =====")
        print("Atoms added by envelope:",
            env_added_mask.sum())

        print("\nSpecies added by envelope:")
        print(
            pd.Series(
                species[env_added_mask]
            ).value_counts()
        )

        print("Envelope ions added:",
          enveloped_mask.sum())

        print(
            "Non-core clustered ions after envelope:",
            np.sum(
                (~core_mask) &
                (labels_env != -1)
            )
        )

    # ---------------------------------
    # Erosion
    # ---------------------------------
    vprint("Eroding enveloped ions")

    labels_final = labels_env.copy()

    # store nearest-background distance for diagnostics
    erosion_distance = np.full(n_atoms, np.nan)

    all_env_idx = np.where(enveloped_mask)[0]

    if len(all_env_idx) > 0:

        # envelope-species that never got clustered
        unclustered_bulk_mask = (
            np.isin(species, envelope_species)
            &
            (labels_env == -1)
        )

        unclustered_bulk_idx = np.where(
            unclustered_bulk_mask
        )[0]

        if len(unclustered_bulk_idx) > 0:

            bulk_tree = KDTree(
                positions[unclustered_bulk_idx]
            )

            dist, _ = bulk_tree.query(
                positions[all_env_idx],
                k=1
            )

            erosion_distance[all_env_idx] = dist[:, 0]

            remove_mask = (
                dist[:, 0] < d_erosion
            )

            labels_final[
                all_env_idx[remove_mask]
            ] = -1

    # ---------------------------------
    # Erosion statistics
    # ---------------------------------

    erosion_removed_mask = (
        (labels_env != -1) &
        (labels_final == -1)
    )

    if verbose:
        print("\n===== EROSION STATS =====")
        print(
            "Atoms removed by erosion:",
            erosion_removed_mask.sum()
        )

        print("\nSpecies removed by erosion:")
        print(
            pd.Series(
                species[erosion_removed_mask]
            ).value_counts()
        )

        print(
            "Envelope ions:",
            enveloped_mask.sum()
        )

        print(
            "Envelope ions surviving:",
            np.sum(
                enveloped_mask &
                (labels_final != -1)
            )
        )

        print(
            "Envelope ions removed:",
            np.sum(
                enveloped_mask &
                (labels_final == -1)
            )
        )

    # diagnostic output
    df["erosion_distance"] = erosion_distance

    # ---------------------------------
    # N_min filter (C++ equivalent)
    # ---------------------------------

    vprint("Applying N_min filter")

    final_cluster_ids = np.unique(labels_final)
    final_cluster_ids = final_cluster_ids[
        final_cluster_ids != -1
    ]

    for cid in final_cluster_ids:

        cluster_size = np.sum(
            labels_final == cid
        )

        if cluster_size < N_min:

            labels_final[
                labels_final == cid
            ] = -1

    # ---------------------------------
    # Store results
    # ---------------------------------

    if verbose:
        print("\n===== BALANCE =====")
        print("DBSCAN atoms:", np.sum(labels != -1))
        print("After envelope:", np.sum(labels_env != -1))
        print("After erosion:", np.sum(labels_final != -1))
        print(
            "Net gain:",
            np.sum(labels_final != -1) - np.sum(labels != -1)
        )

    df["cluster_id"] = labels_final

    return df


def analyse_clusters(
    clustered_df,
    ion_formula_map: dict = None,
    atomic_volume=0.0118  # nm³/atom
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
            "cluster_id", "x", "y", "z", "size", "radius_of_gyration",
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
    k: int,
    n_min: int,
    envelope_radius: float,
    d_erosion: float,
    n_relabel_runs: int,
    run_name: str,
    output_dir,
    relabel_seed: int = 0,
    tree: KDTree = None,
    tree_positions: np.ndarray = None,
    cached_shuffles=None,
    dclassify: float = None,
    atomic_volume: float = 0.0118,
    verbose: bool = False,
) -> dict:
    """
    Run one real cluster search plus n_relabel_runs randomised (mc-shuffled)
    control searches, at a single (dmax, k) parameter combination.

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

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    random_dir = output_dir / f"{run_name}_random"
    random_dir.mkdir(parents=True, exist_ok=True)

    if tree is None:
        tree, tree_positions = build_position_tree(pos_df)

    formula_map = build_ion_formula_map(rrng_df)

    real_atoms_path = output_dir / f"{run_name}_real_clustered_atoms.parquet"
    real_stats_path = output_dir / f"{run_name}_real_cluster_stats.parquet"

    if not (real_atoms_path.exists() and real_stats_path.exists()):
        labelled = assign_rrng_to_pos(pos_df, rrng_df)
        clustered = run_clustering(
            labelled, core_ions, bulk_ions, dmax, k, n_min,
            envelope_radius, d_erosion,
            dclassify=dclassify, tree=tree, tree_positions=tree_positions,
            verbose=verbose,
        )
        stats = analyse_clusters(
            clustered, ion_formula_map=formula_map, atomic_volume=atomic_volume
        )
        clustered.to_parquet(real_atoms_path, index=False)
        stats.to_parquet(real_stats_path, index=False)
        del labelled, clustered, stats

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
            labelled, core_ions, bulk_ions, dmax, k, n_min,
            envelope_radius, d_erosion,
            dclassify=dclassify, tree=tree, tree_positions=tree_positions,
            verbose=False,
        )

        sizes = (
            clustered.loc[clustered["cluster_id"] != -1]
            .groupby("cluster_id")
            .size()
            .reset_index(name="size")
        )
        sizes.to_parquet(run_path, index=False)

        del shuffled_pos, labelled, clustered, sizes

    return {
        "run_name": run_name,
        "dmax": dmax,
        "k": k,
        "n_min": n_min,
        "real_clustered_atoms_path": str(real_atoms_path),
        "real_cluster_stats_path": str(real_stats_path),
        "random_dir": str(random_dir),
        "n_relabel_runs": n_relabel_runs,
    }
