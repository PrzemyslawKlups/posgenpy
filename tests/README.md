# Tests

Two tiers.

## `test_cluster_semantics.py` / `test_fast_sweep.py` / `test_corrections.py` / `test_grid_analysis.py` — fast, run anywhere

Synthetic point sets of a few dozen atoms, one assertion per rule read out of
posgen's `clusterOps.cpp`. No external data, runs in ~1.5 s:

```bash
pytest tests/test_cluster_semantics.py
```

These are the regression net for the clustering rules — the off-by-one in `k`,
which atoms seed the envelope, what erodes what, strict `<` on distances, and
first-match-wins range assignment. Coordinates are exact binary fractions
(1/128, 0.25, 0.375, …) so the boundary cases test the boundary rather than
floating-point slack.

`test_fast_sweep.py` pins the claim `fast_sweep.run_fast_grid_sweep` rests on:
that one kNN query and one radius query per dataset answer every parameter
cell. On a synthetic tip of a few thousand atoms it asserts that the shortcut
picks the same core atoms as the full ranging path, classifies the same core
points as `DBSCAN(eps=dmax, min_samples=order)` at every cell, clusters the
same atoms as `run_clustering`, and differs from it only on border atoms two
clusters could equally claim — plus that the envelope and erosion provably
cannot change the core size distribution, which is why skipping them is safe.

Three of its assertions are about the specific ways that shortcut could be
wrong, and are worth keeping if the file is ever reorganised:

- `test_raising_the_order_splits_a_cluster` — a dumbbell that is one cluster at
  order 3 and two at order 10, checked against `run_clustering`. Only the dmax
  axis is walked monotonically; if anything incremental ever crept into the
  order axis, a split would silently come back as one cluster.
- `test_relabelling_equals_shuffling_mc_and_reranging` — that permuting the
  ranged labels equals ranging a permuted `mc` column, which is what makes the
  control datasets the standard relabelling rather than a stand-in for it.
- `test_control_runs_match_the_full_pipeline` — both paths given the same
  `cached_shuffles`, so they cluster the *same* relabelled datasets and their
  control-run size distributions can be compared directly rather than only in
  distribution.

`test_corrections.py` covers the final run — the one whose numbers get
reported. Its corrections run in an order that matters (Fm deconvolution before
anything is derived from the counts, then the size columns, then edge
detection), and each step can produce plausible-but-wrong numbers if it moves:
that Fm deconvolution conserves Fe+Ni+Fm and empties the peak, that the matrix
summary credits compound ions per element the same way `analyse_clusters` does,
and that `n_min` passed to `run_clustering` clips core+border counts *before*
the envelope, so a cluster under the threshold cannot first collect a bulk shell
that lifts it back over.

`test_grid_analysis.py` covers what is read *off* a finished grid, for the
`PosgenPy_Grid_Analysis_Example_PythonOnly.ipynb` workflow: the per-cell numbers
the heatmaps colour, the per-cluster table the distributions draw, and the
atom-by-atom comparison of neighbouring cells. Three of its assertions are the
ones worth keeping:

- `test_normalised_separation_is_exact_not_a_knn_approximation` — the
  radius-normalised cluster separation is checked against the brute-force
  all-pairs answer on a spread of radii, because the k-nearest search it starts
  from is only a first guess: a large cluster past the k-th neighbour can crowd
  a small one more than its nearest neighbour does.
- `test_raising_dmax_merges_and_never_splits` — DBSCAN is monotone in dmax, so
  a split reported along that axis would mean the sweep itself had gone wrong.
  The order axis is checked the other way round, on a tip with a dumbbell in it.
- `test_split_counts_add_up` — every cluster on each side of a step is
  accounted for exactly once, so `n_split` / `n_merged` / `n_dissolved` /
  `n_created` cannot quietly double-count the same cluster.

```bash
pytest tests/test_fast_sweep.py tests/test_grid_analysis.py
```

## `test_posgen_parity.py` — slow, needs the reference datasets

Compares the engine atom-for-atom against real posgen output. The reference
data (the `.POS` files and the posgen runs made from them) does not live in
this repo, so the module skips unless you point it at them:

```bash
export POSGENPY_PARITY_CONFIG=/path/to/parity_datasets.json
pytest tests/test_posgen_parity.py -m slow
```

or drop the file at `tests/parity_datasets.json` (gitignored — it holds
absolute paths). Copy `parity_datasets.example.json` and edit the four paths.

Each entry names a `.POS`, its `.RRNG`, and the `_cluster.pos` /
`_clusterID.pos` pair posgen wrote, plus the parameters from that run's XML.
`order` is posgen's Order directly: a run labelled `order_7`
(`<dclassify knn="7">`) is configured with `"order": 7`.

Optional `expected` values (`n_clustered`, `n_clusters`, `n_core`, `n_bulk`)
pin absolute numbers. For the `_top` tip they come from the run's own
`sizedist.txt` (`Σ size × frequency` = 28,387 core+border atoms in 1,244
clusters) and `cluster.pos` (45,194 atoms total).

Do **not** use the `cluster-stats.txt` / `unclustered-stats.txt` files sitting
next to those runs as a reference — they were rewritten in place by the old
notebook's post-processing (a zeroed `Fm` column is the giveaway). Only the
binary `.pos` files and `sizedist.txt` are raw posgen output.
