# Tests

Two tiers.

## `test_cluster_semantics.py` — fast, runs anywhere

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
