"""Compare the engine against real posgen output, atom for atom.

Needs the reference datasets (real .POS files and the posgen runs made from
them), which do not live in this repo - see tests/conftest.py for how to point
the suite at them. Without that config the whole module skips.

Run with:  pytest tests/test_posgen_parity.py -m slow
"""
import functools
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from posgenpy import (
    compare_clusterings,
    load_posgen_clustering,
    read_pos_to_df,
    read_rrng,
    run_engine_variant,
)

pytestmark = pytest.mark.slow

# A handful of atoms can legitimately land differently, and no amount of
# fixing will remove them: when two clusters are both within reach of the same
# border or bulk atom, posgen gives it to whichever cluster its KD-tree
# traversal reaches first, and sklearn's DBSCAN to whichever *its* traversal
# reaches first. Neither order is reproducible from the other.
#
# Observed effect on the four reference tips: 3 of 4 match exactly; the
# `_bottom` tip differs by one atom in 19,097 (0.005%), where posgen's
# traversal claimed four border ions for a cluster our run assigned elsewhere,
# leaving a lone core point that N_min=2 then clipped.
#
# The tolerances below are wide enough for that and far too tight for a real
# regression - the envelope-seeding bug these tests were written for showed up
# as a 16.6% mismatch.
DEFAULT_ATOM_TOLERANCE = 1e-3
DEFAULT_PARTITION_AGREEMENT = 0.99
DEFAULT_SIZEDIST_TOLERANCE = 5e-3


def _tolerance(dataset, key, default):
    return dataset.get("tolerance", {}).get(key, default)


@functools.lru_cache(maxsize=4)
def _load(dataset_key):
    """Load and cluster one dataset once, reused across the assertions below."""
    dataset = dict(dataset_key)
    params = dict(dataset["params"])
    params["core_species"] = list(params["core_species"])
    params["envelope_species"] = list(params["envelope_species"])

    pos = read_pos_to_df(dataset["pos"])
    rrng = read_rrng(dataset["rrng"])

    old = load_posgen_clustering(dataset["cluster_pos"], dataset["clusterid_pos"], rrng)
    new = run_engine_variant(pos, rrng, params, cache_dir=dataset.get("cache_dir"))

    return old, new, compare_clusterings(old, new, label=dataset["name"])


def _freeze(obj):
    """Make the dataset dict hashable so it can key an lru_cache."""
    if isinstance(obj, dict):
        return tuple(sorted((k, _freeze(v)) for k, v in obj.items()))
    if isinstance(obj, list):
        return tuple(_freeze(v) for v in obj)
    return obj


@pytest.fixture
def loaded(parity_dataset):
    return _load(_freeze(parity_dataset))


def test_no_atom_is_clustered_only_by_posgen(loaded, parity_dataset):
    """Nothing posgen found may be missing from our result."""
    _, _, report = loaded
    tol = _tolerance(parity_dataset, "atoms", DEFAULT_ATOM_TOLERANCE)

    assert report.only_in_old_fraction <= tol, report.summary(show=False)


def test_no_atom_is_clustered_only_by_the_engine(loaded, parity_dataset):
    """And we must not invent clustered atoms posgen did not find.

    This is the assertion that was failing at 16.6% before the envelope was
    corrected to seed from classified core points only.
    """
    _, _, report = loaded
    tol = _tolerance(parity_dataset, "atoms", DEFAULT_ATOM_TOLERANCE)

    assert report.only_in_new_fraction <= tol, report.summary(show=False)


def test_only_bulk_atoms_may_differ(loaded, parity_dataset):
    """Whatever does differ must be an envelope-stage atom.

    Core linking is connected-components over classified core points, which is
    order independent - so a core-species atom appearing on only one side is a
    real bug, not a tie-break. (The one known exception is an orphaned core
    point whose border ions were claimed by another cluster; allow it only
    within the atom tolerance.)
    """
    _, _, report = loaded
    core = list(parity_dataset["params"]["core_species"])
    tol = _tolerance(parity_dataset, "atoms", DEFAULT_ATOM_TOLERANCE)

    for frame, total in ((report.only_in_old, report.n_old),
                         (report.only_in_new, report.n_new)):
        if not len(frame):
            continue
        n_core = int(frame["ion"].isin(core).sum())
        assert n_core / total <= tol, (
            f"{n_core} core-species atoms differ - that is a structural "
            f"difference, not an assignment tie-break\n"
            + report.summary(show=False)
        )


def test_clustered_atom_counts_match(loaded, parity_dataset):
    _, _, report = loaded
    tol = _tolerance(parity_dataset, "atoms", DEFAULT_ATOM_TOLERANCE)

    assert abs(report.n_new - report.n_old) / report.n_old <= tol

    expected = parity_dataset.get("expected", {})
    if "n_clustered" in expected:
        assert report.n_new == expected["n_clustered"]


def test_cluster_count_matches(loaded, parity_dataset):
    _, _, report = loaded
    tol = _tolerance(parity_dataset, "atoms", DEFAULT_ATOM_TOLERANCE)

    delta = abs(report.n_clusters_new - report.n_clusters_old)
    assert delta <= max(1, tol * report.n_clusters_old)

    expected = parity_dataset.get("expected", {})
    if "n_clusters" in expected:
        assert report.n_clusters_new == expected["n_clusters"]


def test_core_and_bulk_split_matches_expected(loaded, parity_dataset):
    """Check the core/envelope split, not just the total.

    Two engines could agree on the total while disagreeing about which stage
    contributed what, so pin both halves. The expected numbers come from the
    run's own sizedist.txt (core+border) and cluster.pos.
    """
    old, new, _ = loaded
    params = parity_dataset["params"]

    for frame, who in ((old, "posgen"), (new[new["cluster_id"] != -1], "engine")):
        n_core = int(frame["ion"].isin(params["core_species"]).sum())
        n_bulk = int(frame["ion"].isin(params["envelope_species"]).sum())
        n_other = len(frame) - n_core - n_bulk

        assert n_other == 0, f"{who} clustered species outside both typelists"

        expected = parity_dataset.get("expected", {})
        if "n_core" in expected:
            assert n_core == expected["n_core"], who
        if "n_bulk" in expected:
            assert n_bulk == expected["n_bulk"], who


def test_size_distribution_matches_posgen_sizedist(loaded, parity_dataset):
    """Cluster counts vs n_min must match posgen's own sizedist.txt.

    posgen's `sizedist.txt` tabulates `Cluster Size(core)` - core+border
    atoms, NOT the enveloped bulk - and that is what the old real-cluster
    count graphs were thresholded on. Counting total cluster size instead
    keeps far more clusters at every n_min (on this data, 168 vs 136 at
    n_min=8), which is what made the new graphs sit above the old ones.
    """
    sizedist = parity_dataset.get("sizedist")
    if not sizedist or not Path(sizedist).exists():
        pytest.skip("no sizedist.txt configured for this dataset")

    _, new, _ = loaded
    core = list(parity_dataset["params"]["core_species"])

    sd = pd.read_csv(sizedist, sep="\t")
    sd.columns = [c.strip() for c in sd.columns]
    posgen_core = np.repeat(
        sd["Cluster Size(core)"].to_numpy(), sd["Frequency"].to_numpy()
    )

    clustered = new[new["cluster_id"] != -1]
    engine_core = (
        clustered[clustered["ion"].isin(core)].groupby("cluster_id").size()
    )

    tol = _tolerance(
        parity_dataset, "sizedist", DEFAULT_SIZEDIST_TOLERANCE
    )
    allowance = max(5, tol * len(posgen_core))

    # Just above the order threshold the size distribution is a cliff - most
    # clusters have barely more core atoms than `order` demands - so a single
    # border atom changing hands moves a whole cluster across the n_min line.
    # On the `_top` tip that costs at most 5 clusters of 1,244 (0.4%), and
    # nothing at all for n_min >= 16. The metric bug this test guards against
    # was 24% at n_min=8, so the allowance is ~50x tighter than the failure.
    for n_min in range(2, 22, 2):
        expected = int((posgen_core >= n_min).sum())
        actual = int((engine_core >= n_min).sum())
        assert abs(actual - expected) <= allowance, (
            f"n_min={n_min}: posgen {expected}, engine {actual}"
        )


def test_partition_agrees_up_to_renumbering(loaded, parity_dataset):
    """The two engines must induce essentially the same partition.

    Cluster ids themselves are meaningless across engines (posgen numbers from
    1 in KD-tree traversal order, this engine from 0 in DBSCAN order), so
    agreement is measured through the dominant old->new mapping. It will not
    be exactly 1.0: atoms reachable from two clusters go to whichever the
    traversal reaches first, which differs between implementations.
    """
    _, _, report = loaded
    threshold = _tolerance(
        parity_dataset, "partition_agreement", DEFAULT_PARTITION_AGREEMENT
    )
    agreement = report.partition_agreement()

    assert agreement["forward"] >= threshold and agreement["reverse"] >= threshold, (
        f"forward={agreement['forward']:.4%} reverse={agreement['reverse']:.4%} "
        f"split={agreement['split']} merged={agreement['merged']}"
    )
