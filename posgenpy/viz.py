"""3D visualisation of selected clusters, reading persisted per-atom parquet
files directly by cluster_id (no re-reading/re-labelling pos+range files
from scratch, unlike the old file-based show_clusters_in_pos).
"""
from pathlib import Path

import pandas as pd
import plotly.express as px


def show_clusters_in_pos(clustered_atoms, clusters_to_plot: list, ion_colors: dict = None):
    """
    3D scatter of one or more selected clusters.

    clustered_atoms: a path to a `*_real_clustered_atoms.parquet` file (as
    returned in a sweep/run_relabel_sweep manifest entry), or an
    already-loaded DataFrame with x, y, z, ion, cluster_id columns.
    """

    if isinstance(clustered_atoms, (str, Path)):
        df = pd.read_parquet(clustered_atoms)
    else:
        df = clustered_atoms

    selected = df.loc[df["cluster_id"].isin(clusters_to_plot)]

    fig = px.scatter_3d(
        selected, x="x", y="y", z="z", color="ion",
        color_discrete_map=ion_colors or {},
    )
    fig.update_traces(marker_size=2)
    fig.update_layout(title=f"Cluster ID(s): {clusters_to_plot}")
    fig.show()

    return selected


def compare_selected_clusters_across_param(
    sweep_manifest: list,
    swept_param_name: str,
    chosen_value,
    selected_clusters: list,
):
    """
    3D scatter comparing the spatial neighbourhood of `selected_clusters`
    (chosen from the swept value `chosen_value`) across every swept value in
    `sweep_manifest`, by filtering every value's clustered atoms to the
    bounding box of the selected clusters.
    """

    chosen_entries = [e for e in sweep_manifest if e["swept_value"] == chosen_value]
    if len(chosen_entries) != 1:
        print("Found 0 or 2< entries matching chosen_value. Can't continue")
        return None

    chosen_df = pd.read_parquet(chosen_entries[0]["real_clustered_atoms_path"])
    chosen_df = chosen_df.loc[chosen_df["cluster_id"].isin(selected_clusters)]

    xmin, xmax = chosen_df["x"].min(), chosen_df["x"].max()
    ymin, ymax = chosen_df["y"].min(), chosen_df["y"].max()
    zmin, zmax = chosen_df["z"].min(), chosen_df["z"].max()

    frames = []
    for entry in sweep_manifest:
        df = pd.read_parquet(entry["real_clustered_atoms_path"])
        mask = (
            (df["x"] >= xmin) & (df["x"] <= xmax)
            & (df["y"] >= ymin) & (df["y"] <= ymax)
            & (df["z"] >= zmin) & (df["z"] <= zmax)
        )
        df = df.loc[mask].copy()
        df[swept_param_name] = str(entry["swept_value"])
        df["cluster_label"] = df["cluster_id"].astype(str)
        frames.append(df)

    df_to_plot = pd.concat(frames, ignore_index=True)

    fig = px.scatter_3d(
        df_to_plot, x="x", y="y", z="z",
        color=swept_param_name, symbol="cluster_label",
    )
    fig.update_traces(marker_size=3, opacity=0.8)

    return fig
