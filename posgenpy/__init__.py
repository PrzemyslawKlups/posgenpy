"""posgenpy: pure-Python APT cluster-search analysis (no posgen/libatomprobe dependency).

`import posgenpy as pgp` keeps working exactly as before; the implementation
is now split across io.py, cluster.py, sweep.py, corrections.py, viz.py.
"""
from .io import (
    read_pos_to_df,
    range_info_from_rrng,
    read_rrng,
    read_apt_range_file,
    build_ion_formula_map,
    get_ion_colours_from_rrng,
    assign_rrng_to_pos,
)
from .cluster import (
    build_position_tree,
    run_clustering,
    analyse_clusters,
    run_relabel_sweep,
)
from .sweep import (
    seconds_to_hhmmss,
    run_parameter_sweep,
    optimise_n_min,
    plot_real_cluster_ratio_across_swept_param,
    plot_real_clusters_across_swept_param,
    plot_cluster_composition_across_swept_param_absolute,
    plot_cluster_composition_across_swept_param_relative,
)
from .corrections import (
    calculate_cluster_size,
    deconvolute_Fm_clusters,
    remove_edge_clusters,
    summarise_unclustered,
    find_knn_cluster_distance,
    get_cluster_metrics,
    in_hull,
    not_in_hull,
)
from .viz import (
    show_clusters_in_pos,
    compare_selected_clusters_across_param,
)
from .parity import (
    load_posgen_clustering,
    compare_clusterings,
    run_engine_variant,
    variant_matrix,
    ParityReport,
)

__all__ = [
    "read_pos_to_df",
    "range_info_from_rrng",
    "read_rrng",
    "read_apt_range_file",
    "build_ion_formula_map",
    "get_ion_colours_from_rrng",
    "assign_rrng_to_pos",
    "build_position_tree",
    "run_clustering",
    "analyse_clusters",
    "run_relabel_sweep",
    "seconds_to_hhmmss",
    "run_parameter_sweep",
    "optimise_n_min",
    "plot_real_cluster_ratio_across_swept_param",
    "plot_real_clusters_across_swept_param",
    "plot_cluster_composition_across_swept_param_absolute",
    "plot_cluster_composition_across_swept_param_relative",
    "calculate_cluster_size",
    "deconvolute_Fm_clusters",
    "remove_edge_clusters",
    "summarise_unclustered",
    "find_knn_cluster_distance",
    "get_cluster_metrics",
    "in_hull",
    "not_in_hull",
    "show_clusters_in_pos",
    "compare_selected_clusters_across_param",
    "load_posgen_clustering",
    "compare_clusterings",
    "run_engine_variant",
    "variant_matrix",
    "ParityReport",
]
