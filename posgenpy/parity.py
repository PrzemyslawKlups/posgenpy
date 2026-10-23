"""Compare this engine's cluster search against real posgen output.

posgen wrote, per run, a `_cluster.pos` (the clustered atoms) and a
`_clusterID.pos` (same atoms, same order, with the cluster id in the mass
column). Together those are the only *raw* per-atom ground truth available -
the `cluster-stats.txt` / `unclustered-stats.txt` files that sit next to them
were rewritten in place by the old notebook's post-processing (a zeroed `Fm`
column is the giveaway of `deconvolute_Fm_clusters`), so they must not be used
as a reference for per-atom comparisons.

Typical use::

    old = load_posgen_clustering(cluster_pos, clusterid_pos, rrng_df)
    new = run_engine_variant(pos_df, rrng_df, params)
    report = compare_clusterings(old, new)
    report.summary()
"""
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .cluster import build_position_tree, run_clustering
from .io import assign_rrng_to_pos, read_pos_to_df

KEY_COLUMNS = ["x", "y", "z", "mc"]


def load_posgen_clustering(cluster_pos_path, clusterid_pos_path, rrng_df=None):
    """Load a posgen run's `_cluster.pos` + `_clusterID.pos` as one frame.

    Returns per-atom `x, y, z, mc, cluster_id` (plus `ion` if `rrng_df` is
    given). posgen writes the two files in the same order, so the cluster id
    is taken positionally - but that is verified against the coordinates
    rather than assumed, and a coordinate join is used as a fallback.
    """
    cluster_pos = read_pos_to_df(cluster_pos_path)
    clusterid_pos = read_pos_to_df(clusterid_pos_path)

    if len(cluster_pos) != len(clusterid_pos):
        raise ValueError(
            f"{Path(cluster_pos_path).name} has {len(cluster_pos)} atoms but "
            f"{Path(clusterid_pos_path).name} has {len(clusterid_pos)}"
        )

    aligned = np.allclose(
        cluster_pos[["x", "y", "z"]].to_numpy(),
        clusterid_pos[["x", "y", "z"]].to_numpy(),
    )

    df = cluster_pos.copy()

    if aligned:
        df["cluster_id"] = clusterid_pos["mc"].to_numpy().round().astype(int)
    else:
        # fall back to matching on coordinates
        ids = clusterid_pos.rename(columns={"mc": "cluster_id"})
        ids["cluster_id"] = ids["cluster_id"].round().astype(int)
        df = df.merge(ids[["x", "y", "z", "cluster_id"]], on=["x", "y", "z"], how="left")
        if df["cluster_id"].isna().any():
            raise ValueError(
                "cluster.pos and clusterID.pos are neither row-aligned nor "
                "joinable on coordinates"
            )
        df["cluster_id"] = df["cluster_id"].astype(int)

    if rrng_df is not None:
        df = assign_rrng_to_pos(df, rrng_df)

    df.attrs["row_aligned"] = aligned
    return df


def _make_key(df, round_digits):
    return pd.MultiIndex.from_arrays(
        [df[c].to_numpy().round(round_digits) for c in KEY_COLUMNS]
    )


@dataclass
class ParityReport:
    """Outcome of comparing a posgen run against an engine run."""

    n_old: int
    n_new: int
    only_in_old: pd.DataFrame
    only_in_new: pd.DataFrame
    shared: pd.DataFrame
    n_clusters_old: int
    n_clusters_new: int
    label: str = ""
    notes: dict = field(default_factory=dict)

    @property
    def only_in_new_fraction(self):
        return len(self.only_in_new) / self.n_new if self.n_new else 0.0

    @property
    def only_in_old_fraction(self):
        return len(self.only_in_old) / self.n_old if self.n_old else 0.0

    @property
    def matches(self):
        return len(self.only_in_old) == 0 and len(self.only_in_new) == 0

    def change_types(self):
        """How shared atoms moved between clusterings."""
        s = self.shared
        old_c, new_c = s["old_cluster"] != -1, s["new_cluster"] != -1
        return pd.Series(
            {
                "same": int((old_c & new_c & (s["old_cluster"] == s["new_cluster"])).sum()),
                "reassigned": int((old_c & new_c & (s["old_cluster"] != s["new_cluster"])).sum()),
                "lost": int((old_c & ~new_c).sum()),
                "gained": int((~old_c & new_c).sum()),
            }
        )

    def partition_agreement(self):
        """How well the two clusterings agree once cluster ids are ignored.

        posgen numbers clusters from 1 in KD-tree traversal order, this engine
        from 0 in DBSCAN's order, so comparing ids directly is meaningless.
        What matters is whether the two produce the same *partition*: returns
        the fraction of jointly-clustered atoms that sit on the dominant
        old->new mapping (1.0 means the clusterings are identical up to
        renaming), plus the reverse fraction and how many clusters split or
        merged.
        """
        s = self.shared[(self.shared["old_cluster"] != -1) & (self.shared["new_cluster"] != -1)]

        if s.empty:
            return {"forward": float("nan"), "reverse": float("nan"),
                    "split": 0, "merged": 0, "is_pure_relabelling": False}

        overlap = pd.crosstab(s["old_cluster"], s["new_cluster"])
        counts = overlap.to_numpy()
        total = counts.sum()

        forward = counts.max(axis=1).sum() / total
        reverse = counts.max(axis=0).sum() / total
        split = int((counts > 0).sum(axis=1).max() > 1 and ((counts > 0).sum(axis=1) > 1).sum())
        merged = int((counts > 0).sum(axis=0).max() > 1 and ((counts > 0).sum(axis=0) > 1).sum())

        return {
            "forward": forward,
            "reverse": reverse,
            "split": split,
            "merged": merged,
            "is_pure_relabelling": bool(
                np.isclose(forward, 1.0) and np.isclose(reverse, 1.0)
            ),
        }

    def cluster_jaccard(self):
        """Best Jaccard overlap for each old cluster against any new cluster."""
        s = self.shared[(self.shared["old_cluster"] != -1) & (self.shared["new_cluster"] != -1)]
        if s.empty:
            return pd.DataFrame(columns=["old_cluster", "size", "best_jaccard"])

        overlap = pd.crosstab(s["old_cluster"], s["new_cluster"])
        old_sizes = overlap.sum(axis=1)
        new_sizes = overlap.sum(axis=0)
        union = old_sizes.to_numpy()[:, None] + new_sizes.to_numpy()[None, :] - overlap.to_numpy()
        jaccard = np.divide(overlap.to_numpy(), union, out=np.zeros_like(union, dtype=float), where=union > 0)

        return pd.DataFrame(
            {
                "old_cluster": overlap.index,
                "size": old_sizes.to_numpy(),
                "best_jaccard": jaccard.max(axis=1),
            }
        )

    def summary(self, show=True):
        lines = []
        add = lines.append

        add(f"=== parity report {self.label} ===")
        add(f"clustered atoms   old {self.n_old:>8,}   new {self.n_new:>8,}   "
            f"diff {self.n_new - self.n_old:+,}")
        add(f"clusters          old {self.n_clusters_old:>8,}   new {self.n_clusters_new:>8,}")
        add(f"only in old       {len(self.only_in_old):>8,}  ({self.only_in_old_fraction:.2%})")
        add(f"only in new       {len(self.only_in_new):>8,}  ({self.only_in_new_fraction:.2%})")

        for name, frame in (("only in old", self.only_in_old), ("only in new", self.only_in_new)):
            if len(frame) and "ion" in frame.columns:
                counts = frame["ion"].fillna("<unranged>").value_counts()
                add(f"\n{name}, by species:")
                add("  " + counts.to_string().replace("\n", "\n  "))
            if len(frame) and "erosion_distance" in frame.columns:
                n_nan = int(frame["erosion_distance"].isna().sum())
                add(f"  atoms with no erosion_distance (i.e. not envelope-added): {n_nan}")

        ct = self.change_types()
        add("\nshared atoms:")
        add(f"  clustered in both      {ct['same'] + ct['reassigned']:>8,}")
        add(f"  clustered only in old  {ct['lost']:>8,}")
        add(f"  clustered only in new  {ct['gained']:>8,}")

        agree = self.partition_agreement()
        add("\ncluster partition (ignoring id numbering, which differs by engine):")
        add(f"  atoms on the dominant old->new mapping  {agree['forward']:.4%}")
        add(f"  atoms on the dominant new->old mapping  {agree['reverse']:.4%}")
        add(f"  old clusters split across new clusters  {agree['split']}")
        add(f"  new clusters merging several old ones   {agree['merged']}")
        if agree["is_pure_relabelling"]:
            add("  -> identical partition, cluster ids are just renumbered")

        for key, value in self.notes.items():
            add(f"\n{key}: {value}")

        text = "\n".join(lines)
        if show:
            print(text)
        return text


def compare_clusterings(old_df, new_df, round_digits=5, label=""):
    """Align two per-atom clusterings and quantify how they differ.

    Both frames need `x, y, z, mc, cluster_id`; `ion` and `erosion_distance`
    are used for diagnostics when present. Only atoms that are *clustered*
    (`cluster_id != -1`) count towards `only_in_old` / `only_in_new`; the
    shared set covers every atom present in both frames.
    """
    old_key = _make_key(old_df, round_digits)
    new_key = _make_key(new_df, round_digits)

    old = old_df.copy()
    new = new_df.copy()
    old.index = old_key
    new.index = new_key

    old_clustered = old[old["cluster_id"] != -1]
    new_clustered = new[new["cluster_id"] != -1]

    only_in_old = old_clustered[~old_clustered.index.isin(new_clustered.index)]
    only_in_new = new_clustered[~new_clustered.index.isin(old_clustered.index)]

    # every atom the engine saw, with the old label where posgen had one
    shared = pd.DataFrame(
        {
            "new_cluster": new["cluster_id"],
            "old_cluster": old["cluster_id"].reindex(new.index).fillna(-1).astype(int),
        }
    )

    notes = {}
    if old.index.duplicated().any() or new.index.duplicated().any():
        notes["warning"] = (
            f"duplicate (x,y,z,mc) keys at round_digits={round_digits}: "
            f"{int(old.index.duplicated().sum())} old, "
            f"{int(new.index.duplicated().sum())} new"
        )
    unmatched = int((~old_clustered.index.isin(new.index)).sum())
    if unmatched:
        notes["old atoms absent from the pos file entirely"] = unmatched
    if "n_unranged_clustered" in new_df.attrs:
        notes["unranged atoms inside new clusters (excluded from output)"] = (
            new_df.attrs["n_unranged_clustered"]
        )

    return ParityReport(
        n_old=len(old_clustered),
        n_new=len(new_clustered),
        only_in_old=only_in_old.reset_index(drop=True),
        only_in_new=only_in_new.reset_index(drop=True),
        shared=shared.reset_index(drop=True),
        n_clusters_old=int(old_clustered["cluster_id"].nunique()),
        n_clusters_new=int(new_clustered["cluster_id"].nunique()),
        label=label,
        notes=notes,
    )


def run_engine_variant(
    pos_df,
    rrng_df,
    params,
    cache_dir=None,
    tree=None,
    tree_positions=None,
    verbose=False,
    **overrides,
):
    """Run `run_clustering` for one parameter set, caching the result.

    `params` holds the run's clustering parameters (`core_species`,
    `envelope_species`, `Dmax`, `k`, `N_min`, `envelope_radius`,
    `d_erosion`); `overrides` supplies the variant flags being explored
    (`envelope_seeds`, `erosion_reference`, `include_unranged`, ...). Results
    are cached under a hash of both, so re-running a variant is free.
    """
    call = {**params, **overrides}

    cache_path = None
    if cache_dir is not None:
        digest = hashlib.sha1(
            json.dumps(call, sort_keys=True, default=str).encode()
        ).hexdigest()[:16]
        cache_path = Path(cache_dir) / f"variant_{digest}.parquet"
        if cache_path.exists():
            return pd.read_parquet(cache_path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)

    if tree is None:
        tree, tree_positions = build_position_tree(pos_df)

    labelled = assign_rrng_to_pos(pos_df, rrng_df)
    result = run_clustering(
        labelled,
        tree=tree,
        tree_positions=tree_positions,
        erosion_diagnostics=True,
        verbose=verbose,
        **call,
    )

    if cache_path is not None:
        result.to_parquet(cache_path, index=False)

    return result


def variant_matrix(
    pos_df,
    rrng_df,
    params,
    old_df,
    variants,
    cache_dir=None,
    verbose=False,
):
    """Run several engine variants against one posgen run and tabulate them.

    `variants` is a list of dicts of overrides passed to `run_engine_variant`.
    Returns a DataFrame with one row per variant plus the reports themselves,
    so the combination that reproduces posgen can be read off directly.
    """
    tree, tree_positions = build_position_tree(pos_df)

    rows = []
    reports = {}

    for variant in variants:
        label = ", ".join(f"{k}={v}" for k, v in variant.items()) or "defaults"
        new_df = run_engine_variant(
            pos_df, rrng_df, params, cache_dir=cache_dir,
            tree=tree, tree_positions=tree_positions, verbose=verbose, **variant
        )
        report = compare_clusterings(old_df, new_df, label=label)
        reports[label] = report

        rows.append(
            {
                "variant": label,
                "n_new": report.n_new,
                "n_old": report.n_old,
                "clusters_new": report.n_clusters_new,
                "only_in_old": len(report.only_in_old),
                "only_in_new": len(report.only_in_new),
                "only_in_new_pct": round(report.only_in_new_fraction * 100, 2),
            }
        )
        print(f"  {label}: new={report.n_new:,} "
              f"(only_in_new={len(report.only_in_new):,}, "
              f"only_in_old={len(report.only_in_old):,})")

    return pd.DataFrame(rows), reports
