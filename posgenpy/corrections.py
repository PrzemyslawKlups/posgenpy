"""Post-clustering corrections: edge-cluster detection, Fm peak deconvolution,
per-tip metrics, and cluster-to-cluster distance.

Operates entirely on in-memory DataFrames produced by cluster.run_clustering /
cluster.analyse_clusters, sharing one consistent `cluster_id` between the
per-atom and per-cluster tables - unlike the old file-based pipeline, there is
no separate reconciliation step needed between posgen's cluster-stats.txt and
clusterID.pos files.
"""
import warnings

import numpy as np
import pandas as pd
import alphashape
from scipy.spatial import ConvexHull, Delaunay, cKDTree
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


def find_knn_cluster_distance(
    cluster_stats_df: pd.DataFrame,
    radius_column: str = "equivalent_spherical_radius",
) -> pd.DataFrame:
    """Attach each cluster's nearest-neighbour cluster distance/ID (by centroid).

    Adds two columns, as before:

    - `closest_cluster_distance` - centre-to-centre distance to the nearest
      other cluster, in nm;
    - `closest_cluster_id` - which cluster that is.

    A raw centre-to-centre distance mixes two effects: how far apart the
    clusters are, and how big they are. Two 4 nm precipitates 3 nm apart are
    touching; two 0.5 nm ones 3 nm apart are nowhere near each other. So when
    `radius_column` is present (`analyse_clusters` writes
    `equivalent_spherical_radius`; `radius_of_gyration` also works), two more
    columns are added:

    - `closest_cluster_distance_norm` - the smallest `d_ij / (R_i + R_j)` over
      every other cluster j. **Below 1 the two clusters' equivalent spheres
      overlap**, i.e. a single precipitate the search has split in two, or two
      that a larger dmax will merge. Above ~2 they are clearly separate.
    - `closest_cluster_id_norm` - the cluster that ratio is against, which is
      not always the nearest one by raw distance: a large cluster further away
      can crowd a small one more than its nearest neighbour does.

    The normalised minimum is exact, not a k-nearest-neighbours approximation:
    a cluster beyond the k-th neighbour can only beat the running best ratio
    `t` if it lies within `t * (R_i + R_max)`, so any cluster whose k-th
    neighbour falls inside that radius gets a second, wider search.
    """

    this_df = cluster_stats_df.copy()
    has_radius = radius_column in this_df.columns

    if len(this_df) < 2:
        this_df["closest_cluster_distance"] = np.nan
        this_df["closest_cluster_id"] = np.nan
        if has_radius:
            this_df["closest_cluster_distance_norm"] = np.nan
            this_df["closest_cluster_id_norm"] = np.nan
        return this_df

    xyz = this_df.loc[:, ["x", "y", "z"]].to_numpy(dtype=float)
    cluster_id_list = this_df["cluster_id"].to_numpy()

    tree = cKDTree(xyz)
    # 32 neighbours is plenty for the raw distance (which needs 1) and covers
    # the normalised search in all but the pathological cases refined below
    k = min(len(this_df), 33)
    distances, indices = tree.query(xyz, k=k)

    this_df["closest_cluster_distance"] = distances[:, 1]
    this_df["closest_cluster_id"] = cluster_id_list[indices[:, 1]]

    if not has_radius:
        return this_df

    radii = np.nan_to_num(this_df[radius_column].to_numpy(dtype=float), nan=0.0)

    neighbour_d = distances[:, 1:]
    neighbour_i = indices[:, 1:]
    radii_sum = radii[:, None] + radii[neighbour_i]
    ratios = np.divide(
        neighbour_d, radii_sum,
        out=np.full_like(neighbour_d, np.inf), where=radii_sum > 0,
    )

    rows = np.arange(len(this_df))
    best = np.argmin(ratios, axis=1)
    best_ratio = ratios[rows, best]
    best_index = neighbour_i[rows, best]

    if k < len(this_df):
        # a cluster past the k-th neighbour can still win, but only if it is
        # big enough: it must sit within best_ratio * (R_i + R_max)
        search_radius = best_ratio * (radii + radii.max())
        for i in np.flatnonzero(search_radius > distances[:, -1]):
            candidates = np.asarray(
                tree.query_ball_point(xyz[i], search_radius[i]), dtype=int
            )
            candidates = candidates[candidates != i]
            if not len(candidates):
                continue

            candidate_sum = radii[i] + radii[candidates]
            candidate_ratio = np.divide(
                np.linalg.norm(xyz[candidates] - xyz[i], axis=1), candidate_sum,
                out=np.full(len(candidates), np.inf), where=candidate_sum > 0,
            )
            closest = int(np.argmin(candidate_ratio))
            if candidate_ratio[closest] < best_ratio[i]:
                best_ratio[i] = candidate_ratio[closest]
                best_index[i] = candidates[closest]

    this_df["closest_cluster_distance_norm"] = np.where(
        np.isfinite(best_ratio), best_ratio, np.nan
    )
    this_df["closest_cluster_id_norm"] = cluster_id_list[best_index]

    return this_df


def cluster_alpha_shapes(
    clustered_atoms_df: pd.DataFrame,
    cluster_ids=None,
    alpha: float = 1.0,
    min_atoms: int = 4,
    progress: bool = False,
) -> dict:
    """Alpha-shape surface mesh of each cluster, for 3D display.

    The alpha shape is the shape the cluster search actually found - the
    convex hull rounds off every concavity, which is exactly where a
    dumbbell-shaped "cluster" that should have been two shows itself.

    Parameters
    ----------
    clustered_atoms_df : pandas.DataFrame
        Per-atom table with x/y/z and cluster_id (`run_clustering` output, or
        a `*_real_clustered_atoms.parquet` file).
    cluster_ids : sequence, optional
        Which clusters to mesh. Defaults to every cluster in the table.
    alpha : float
        The alpha-shape parameter, in **inverse nm**: the probe sphere rolled
        over the atoms has radius `1/alpha`. `1/dmax` is a good starting point
        (a probe the size of the linking distance). Too large and the surface
        tears into holes; `alpha=0` gives the convex hull.
    min_atoms : int
        Clusters with fewer atoms than this are skipped - 4 non-coplanar
        points are the minimum a closed 3D surface can be built from.

    Returns
    -------
    dict[cluster_id, dict]
        Per cluster: `vertices` (n,3), `faces` (m,3), `method` ("alpha" or
        "hull" where the alpha shape came out degenerate), `watertight`,
        `volume` and `area` (nm3/nm2, NaN when the mesh is not watertight).
    """

    if cluster_ids is None:
        cluster_ids = clustered_atoms_df.loc[
            clustered_atoms_df["cluster_id"] != -1, "cluster_id"
        ].unique()
    cluster_ids = list(cluster_ids)

    # only the requested clusters' atoms are grouped: on a real tip the
    # unclustered ones are millions of rows nothing here ever looks at
    selected = clustered_atoms_df.loc[
        clustered_atoms_df["cluster_id"].isin(cluster_ids)
    ]
    grouped = {
        cid: group[["x", "y", "z"]].to_numpy(dtype=float)
        for cid, group in selected.groupby("cluster_id")
    }

    meshes = {}
    skipped = 0

    for cid in tqdm(cluster_ids, disable=not progress, desc="alpha shapes"):
        points = grouped.get(cid)
        if points is None or len(points) < max(4, min_atoms):
            skipped += 1
            continue

        mesh = _alpha_shape_mesh(points, alpha)
        if mesh is None:
            skipped += 1
            continue

        meshes[cid] = mesh

    if skipped:
        warnings.warn(
            f"{skipped} of {len(cluster_ids)} clusters have no usable surface "
            f"at alpha={alpha} (too few atoms, or coplanar) and are not meshed"
        )

    return meshes


def _alpha_shape_mesh(points, alpha):
    """One cluster's alpha-shape mesh, falling back to its convex hull.

    `alphashape` returns a 2D shapely polygon for alpha=0 rather than a 3D
    hull, and a shredded (non-closed) surface when alpha is too large for the
    atom spacing, so both cases are caught here and answered with the convex
    hull - a shape that always exists - rather than with nothing.
    """

    shape = None
    if alpha > 0:
        try:
            shape = alphashape.alphashape(points, float(alpha))
        except Exception:
            shape = None

    faces = np.asarray(getattr(shape, "faces", []), dtype=int)

    if len(faces) >= 4:
        vertices = np.asarray(shape.vertices, dtype=float)
        watertight = bool(getattr(shape, "is_watertight", False))
        return {
            "vertices": vertices,
            "faces": faces,
            "method": "alpha",
            "watertight": watertight,
            "volume": abs(float(shape.volume)) if watertight else np.nan,
            "area": float(shape.area) if watertight else np.nan,
        }

    try:
        hull = ConvexHull(points)
    except Exception:
        return None

    # ConvexHull indexes into the full point array; re-index onto its own
    # vertices so the mesh carries no unused points
    used = np.unique(hull.simplices)
    remap = np.full(len(points), -1, dtype=int)
    remap[used] = np.arange(len(used))

    return {
        "vertices": points[used],
        "faces": remap[hull.simplices],
        "method": "hull",
        "watertight": True,
        "volume": float(hull.volume),
        "area": float(hull.area),
    }


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
