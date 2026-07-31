"""Post-clustering corrections: edge-cluster detection, Fm peak deconvolution,
per-tip metrics, and cluster-to-cluster distance.

Operates entirely on in-memory DataFrames produced by cluster.run_clustering /
cluster.analyse_clusters, sharing one consistent `cluster_id` between the
per-atom and per-cluster tables - unlike the old file-based pipeline, there is
no separate reconciliation step needed between posgen's cluster-stats.txt and
clusterID.pos files.
"""
import numpy as np
import pandas as pd
import alphashape
from scipy.spatial import ConvexHull, Delaunay
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm


def calculate_cluster_size(stats_df: pd.DataFrame, core_ions: list, bulk_ions: list,
                            elements_to_exclude: list) -> pd.Series:
    """Sum of atom counts across the given core/bulk elements (minus any excluded,
    e.g. a matrix element like Fe), per cluster.
    """
    elements_to_include = [
        e for e in core_ions + bulk_ions if e not in elements_to_exclude
    ]

    if stats_df.empty:
        return pd.Series(dtype=int, index=stats_df.index)

    columns_present = [e for e in elements_to_include if e in stats_df.columns]
    return stats_df.loc[:, columns_present].sum(axis=1)


def deconvolute_Fm_clusters(stats_df: pd.DataFrame) -> pd.DataFrame:
    """
    Deconvolute the 29 Da peak labelled "Fm" (Fe/Ni mass overlap) into Fe and
    Ni. The max Ni expected in the Fm peak is derived from natural Ni
    abundance (32% at this mass being the ambiguous isotope, 68% resolvable
    elsewhere): at most Ni/0.32*0.68 of the Fm counts can be Ni; the rest is
    assigned to Fe.
    """

    for col in ("Fm", "Ni", "Fe"):
        if col not in stats_df.columns:
            raise ValueError(f"{col!r} not in the given cluster stats DataFrame")

    max_ni_in_fm = stats_df["Ni"] / 0.32 * 0.68

    ni_to_be_added = max_ni_in_fm.copy()
    smaller_mask = stats_df["Fm"] < max_ni_in_fm
    ni_to_be_added[smaller_mask] = stats_df.loc[smaller_mask, "Fm"]
    ni_to_be_added = ni_to_be_added.round().astype(int)

    fe_to_be_added = stats_df["Fm"] - ni_to_be_added

    new_stats_df = stats_df.copy()
    new_stats_df["Fe"] = new_stats_df["Fe"] + fe_to_be_added
    new_stats_df["Ni"] = new_stats_df["Ni"] + ni_to_be_added
    new_stats_df["Fm"] = 0

    if not np.isclose(
        (new_stats_df["Fe"] + new_stats_df["Ni"]).sum(),
        (stats_df["Fe"] + stats_df["Ni"] + stats_df["Fm"]).sum(),
    ):
        print("Warning: total Fe+Ni+Fm atoms changed after Fm deconvolution.")

    return new_stats_df


def in_hull(p, hull):
    """Test if points `p` (NxK) are inside `hull` (a Delaunay or MxK point array)."""
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p) >= 0


def not_in_hull(p, hull):
    """Test if points `p` (NxK) are outside `hull` (a Delaunay or MxK point array)."""
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p) < 0


def remove_edge_clusters(
    ranged_atoms_df: pd.DataFrame,
    clustered_atoms_df: pd.DataFrame,
    cluster_stats_df: pd.DataFrame,
    alpha: float = 0.1,
    sample_frac: float = 0.01,
    random_state: int = 0,
    print_log: bool = False,
) -> pd.DataFrame:
    """
    Flag clusters that touch the tip's outer boundary as edge clusters, and
    attach per-cluster convex hull volume/area.

    The tip boundary is approximated by an alpha-shape over a random sample
    of all ranged atoms (`sample_frac`, seeded via `random_state` so this is
    reproducible run-to-run, unlike the legacy unseeded implementation it
    replaces).

    Returns a copy of cluster_stats_df with is_edge_cluster, convexhull_volume,
    and convexhull_area columns added. No cross-file cluster-ID reconciliation
    is needed: clustered_atoms_df and cluster_stats_df already share the same
    `cluster_id`.
    """

    if print_log:
        print("Calculating alphashape (this can take a while)")

    sample_points = (
        ranged_atoms_df[["x", "y", "z"]]
        .sample(frac=sample_frac, random_state=random_state)
        .to_numpy()
    )
    alpha_shape = alphashape.alphashape(list(map(tuple, sample_points)), alpha)
    hull_vertices = np.asarray(alpha_shape.vertices)

    stats = cluster_stats_df.reset_index(drop=True).copy()
    is_edge = np.zeros(len(stats), dtype=bool)
    hull_volumes = np.zeros(len(stats))
    hull_areas = np.zeros(len(stats))

    if print_log:
        print("Checking each cluster against the tip boundary")

    for i, cid in enumerate(tqdm(stats["cluster_id"].to_numpy(), disable=not print_log)):

        cluster_points = clustered_atoms_df.loc[
            clustered_atoms_df["cluster_id"] == cid, ["x", "y", "z"]
        ].to_numpy()

        # a convex hull needs at least 4 non-coplanar points in 3D
        if len(cluster_points) < 4:
            continue

        hull = ConvexHull(cluster_points)
        hull_volumes[i] = hull.volume
        hull_areas[i] = hull.area

        shell_points = cluster_points[hull.vertices]
        is_edge[i] = bool(np.any(not_in_hull(shell_points, hull_vertices)))

    stats["convexhull_volume"] = hull_volumes
    stats["convexhull_area"] = hull_areas
    stats["is_edge_cluster"] = is_edge

    if print_log:
        n_edge = int(is_edge.sum())
        pct = round(n_edge / len(stats) * 100, 1) if len(stats) else 0.0
        print(f"Edge clusters: {n_edge} out of {len(stats)} ({pct}%)")

    return stats


def summarise_unclustered(
    clustered_atoms_df: pd.DataFrame, ion_formula_map: dict = None
) -> pd.DataFrame:
    """One-row summary of ranged atoms that ended up unclustered (cluster_id == -1),
    replacing posgen's separate `_unclustered-stats.txt` file for
    get_cluster_metrics' matrix atom count.

    `size` is that count. With `ion_formula_map` (see io.build_ion_formula_map)
    a column per element is added too, credited per-element exactly as
    `analyse_clusters` does - which is what the matrix composition is read off,
    and what lets `deconvolute_Fm_clusters` be applied to the matrix as well as
    to the clusters. Without it, the columns are the raw ion labels.
    """

    unclustered = clustered_atoms_df.loc[clustered_atoms_df["cluster_id"] == -1]
    summary = pd.DataFrame({"size": [len(unclustered)]})

    if "ion" not in clustered_atoms_df.columns:
        return summary

    counts = unclustered["ion"].value_counts()

    if ion_formula_map is None:
        for ion, n_ions in counts.items():
            summary[ion] = n_ions
        return summary

    element_counts = {}
    for ion, n_ions in counts.items():
        for element, per_ion_count in ion_formula_map.get(ion, {ion: 1}).items():
            element_counts[element] = (
                element_counts.get(element, 0) + n_ions * per_ion_count
            )

    for element in sorted(
        {e for formula in ion_formula_map.values() for e in formula}
    ):
        summary[element] = element_counts.get(element, 0)

    return summary


def find_knn_cluster_distance(cluster_stats_df: pd.DataFrame) -> pd.DataFrame:
    """Attach each cluster's nearest-neighbour cluster distance/ID (by centroid)."""

    this_df = cluster_stats_df.copy()

    if len(this_df) < 2:
        this_df["closest_cluster_distance"] = np.nan
        this_df["closest_cluster_id"] = np.nan
        return this_df

    xyz = this_df.loc[:, ["x", "y", "z"]].to_numpy()

    nbrs = NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(xyz)
    d, indices = nbrs.kneighbors(xyz)

    cluster_id_list = this_df["cluster_id"].to_numpy()

    this_df["closest_cluster_distance"] = d[:, 1]
    this_df["closest_cluster_id"] = cluster_id_list[indices[:, 1]]

    return this_df


def get_cluster_metrics(
    reconstruction_name: str,
    reconstruction_location: str,
    sample_type: str,
    cluster_stats_df: pd.DataFrame,
    unclustered_stats_df: pd.DataFrame,
    core_ions: list,
    bulk_ions: list,
    dmax: float,
    order: int,
    n_min: int,
    detector_efficiency: float = 0.52,
    atomic_density_per_nm3: float = 85.49,
    excluded_Fe: bool = False,
    cluster_analysis_type: str = "main",  # "main" or "cv01"
) -> dict:
    """
    One-row summary of cluster analysis metrics (volume fraction, number
    density, etc.), considering edge clusters as half-counted.
    """

    ranged_atoms_in_matrix = unclustered_stats_df["size"].sum()
    ranged_atoms_in_clusters = cluster_stats_df["size"].sum()
    ranged_atoms_in_bulk = ranged_atoms_in_clusters + ranged_atoms_in_matrix

    clustered_volume = ranged_atoms_in_clusters / (detector_efficiency * atomic_density_per_nm3)
    volume_fraction = (
        ranged_atoms_in_clusters / ranged_atoms_in_bulk if ranged_atoms_in_bulk else np.nan
    )

    if excluded_Fe:
        Fe_in_clusters = cluster_stats_df["Fe"].sum() if "Fe" in cluster_stats_df.columns else 0
        clustered_volume_no_Fe = (
            (ranged_atoms_in_clusters - Fe_in_clusters)
            / (detector_efficiency * atomic_density_per_nm3)
        )
        denom = ranged_atoms_in_bulk + Fe_in_clusters
        volume_fraction_no_Fe = (
            (ranged_atoms_in_clusters - Fe_in_clusters) / denom if denom else np.nan
        )

    tip_volume = ranged_atoms_in_bulk / (detector_efficiency * atomic_density_per_nm3)

    if "is_edge_cluster" in cluster_stats_df.columns:
        is_edge = cluster_stats_df["is_edge_cluster"].astype(bool)
    else:
        is_edge = pd.Series(False, index=cluster_stats_df.index)

    number_of_edge_clusters = int(is_edge.sum())
    number_of_all_clusters = len(cluster_stats_df)
    effective_clusters = number_of_all_clusters - 0.5 * number_of_edge_clusters

    number_density = effective_clusters / tip_volume if tip_volume else np.nan
    number_density_error = np.sqrt(effective_clusters) / tip_volume if tip_volume else np.nan

    metrics = {
        "reconstruction name": reconstruction_name,
        "sample type": sample_type,
        "cluster analysis type": cluster_analysis_type,
        "volume fraction": volume_fraction,
        "number density per nm3": number_density,
        "number density per nm3 error": number_density_error,
        "clustered volume nm3": clustered_volume,
        "tip volume nm3": tip_volume,
        "number of all clusters": number_of_all_clusters,
        "number of non-edge clusters": number_of_all_clusters - number_of_edge_clusters,
        "number of edge clusters": number_of_edge_clusters,
        "d_max [nm]": dmax,
        "Order": order,
        "N_min": n_min,
        "core ions": "+".join(core_ions),
        "bulk ions": "+".join(bulk_ions),
        "detector efficiency": detector_efficiency,
        "assumed atomic density per nm3": atomic_density_per_nm3,
        "reconstruction_location": reconstruction_location,
        "Fe excluded?": "Excluded" if excluded_Fe else "Included",
    }

    if excluded_Fe:
        metrics["clustered volume no Fe nm3"] = clustered_volume_no_Fe
        metrics["volume fraction no Fe"] = volume_fraction_no_Fe

    return metrics
