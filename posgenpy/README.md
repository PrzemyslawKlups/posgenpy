# posgenpy (pure-Python engine)

A dependency-free, pure-Python replacement for the `posgen` C++ cluster-search
tool. No external binary, no `libatomprobe`, no XML, no WSL/build toolchain -
just `numpy`/`pandas`/`scikit-learn`/`scipy`.

It implements the same "maximum separation" cluster-search algorithm posgen
used (DBSCAN-style core classification + envelope + erosion), plus the
relabelling/randomisation workflow used to estimate a real-vs-random cluster
threshold and pick `n_min` - all driven from ordinary Jupyter notebooks.

If you're looking for the old posgen/XML-driven workflow (`write_xml_with_relabelling`,
`subprocess.call(["./posgen", ...])`), it's kept, unmodified, under
[`legacy/`](../legacy/) for reference/cross-validation only. It is not imported
by anything in this package.

## Quick start

```python
import posgenpy as pgp

pos_df = pgp.read_pos_to_df("mytip.pos")
rrng_df = pgp.read_rrng("mytip.rrng")

core_ions = ["Ni", "Mn", "Si"]
bulk_ions = ["Fe"]

# the whole (order x dmax) grid, real + 20 relabelled control datasets per cell
grid = pgp.run_fast_grid_sweep(
    pos_df, rrng_df, core_ions,
    orders=[4, 5, 6, 7, 8],            # posgen's Order, as in _order_7
    dmax_values=[0.3, 0.4, 0.5, 0.6],
    bulk_ions=bulk_ions, n_relabel_runs=20,
    output_dir="my_output/", prefix="mytip", n_jobs=-1,
)

# the whole plane at once: optimal n_min and cluster counts per (order, dmax)
summary = pgp.summarise_grid(grid)
pgp.plot_grid_heatmap(summary, value="optimal_n_min").show()

# ...or read the grid along either axis - these take the same manifest shape
# run_parameter_sweep returns, so every sweep function below works unchanged
dmax_sweep = pgp.select_sweep(grid, order=7)

n_min_values = list(range(2, 21, 2))
pgp.plot_real_cluster_ratio_across_swept_param(dmax_sweep, "dmax", n_min_values)

optimal_n_min, real_clusters, random_clusters = pgp.optimise_n_min(dmax_sweep, n_min_values=range(1, 200))

# every cluster of every cell as one row, for hoverable size distributions
sizes = pgp.load_sweep_size_table(dmax_sweep)
pgp.plot_size_distribution(sizes, x="dmax", y="size_core").show()
```

See [`../PosgenPy_Full_Analysis_Example_PythonOnly.ipynb`](../PosgenPy_Full_Analysis_Example_PythonOnly.ipynb)
for the full workflow: order sweep, dmax sweep, n_min optimisation, final
cluster search + corrections (edge-cluster removal, Fm-peak deconvolution,
per-tip metrics), composition/violin/3D plots, and a CV01 parameter-sensitivity
check.

## Installation

```bash
uv pip install -e .
```

(or plain `pip install -e .` - see `pyproject.toml` at the repo root for the
dependency list. `pyarrow` is required for the parquet files this package
reads/writes.)

## Matching posgen

The defaults reproduce posgen's `maxsep` algorithm. Verified atom-for-atom
against five real posgen runs (see [../tests/README.md](../tests/README.md)):

| tip | order | posgen | engine | differing atoms |
|---|---|---|---|---|
| `MNSP tip A _top` | 6 | 45,194 | 45,194 | 0 |
| `MNSP tip A _bottom` | 5 | 19,097 | 19,096 | 1 |
| `MNSP tip B` | 7 | 42,919 | 42,919 | 0 |
| `MNSP tip A _GB` | 10 | 18,948 | 18,948 | 0 |
| `MNSP tip C` | 7 | 26,728 | 26,728 | 0 |

**`order` means what posgen means by Order.** It is the number in
`<dclassify knn="N">` and in a run's `_dmax0.5_order_N` filename, so an
`order_7` posgen run is reproduced with `order=7`. (The older `k` parameter
counted only the *other* neighbours, i.e. `order = k + 1`; it still works but
warns.)

The one differing atom is not a bug that can be fixed: when two clusters can
both reach the same border atom, each implementation gives it to whichever
cluster its KD-tree traversal reaches first, and those orders differ. Cluster
*ids* likewise differ (posgen numbers from 1, this engine from 0), so compare
partitions rather than ids.

### Which "cluster size" n_min is compared against

posgen reports two different sizes, and the old workflow used *both*, which is
easy to get wrong when comparing graphs:

- its `sizedist.txt` counts **core+border atoms** (`Cluster Size(core)`) — this
  is what the old real-cluster-count graphs were thresholded on;
- the old `optimise_n_min` read a `Cluster Size` column that the notebook had
  built as **all ions minus the matrix element** (usually Fe).

`optimise_n_min` and the two real-cluster plots therefore take a `size_metric`:

| `size_metric` | counts | matches |
|---|---|---|
| `"core"` (default) | core+border atoms | posgen's `sizedist.txt` / the old graphs |
| `"total"` | every atom, incl. enveloped bulk | nothing in the old workflow |
| `"total_excluding"` + `exclude_ions=["Fe"]` | all ions minus those listed | the old `optimise_n_min` |

Thresholding `"total"` keeps far more clusters at every n_min and the gap grows
with n_min (on `MNSP tip C` at order 7: 168 clusters vs posgen's 136 at
n_min=8), so a curve plotted that way sits visibly above the posgen-era one.

Sweeps written before this existed have no `size_core` in their random-run
files; `size_metric="core"` then raises rather than silently using the wrong
number. Real-cluster graphs still work on that older output (pass
`core_species=[...]`), but `optimise_n_min` needs the sweep re-running.

Several behaviours are switchable, so a run can be reproduced either way:

| parameter | default (posgen) | alternative |
|---|---|---|
| `envelope_seeds` | `"core"` - only classified core points seed the envelope | `"all"` - every cluster member seeds it, growing a much larger envelope |
| `erosion_reference` | `"bulk_tree"` - untagged bulk-tree atoms | `"envelope_species"` |
| `include_unranged` | `False` | `True` - unranged atoms join the bulk tree |
| `n_min_count` | `"core"` - core+border, clipped before the envelope | `"all"` - whole cluster, clipped at the end |

`include_unranged` defaults to `False` even though the reference runs' XML says
`<unranged foroutput="True" forstats="True"/>`: posgen evidently did not parse
that capitalised `"True"` as true, and `False` is what reproduces its output.

## Modules

- **`io.py`** - reading `.pos`/`.rrng` (and AP Suite-format range) files, and
  labelling atoms with ion species + parsed element formula.
- **`cluster.py`** - `run_clustering` (the core DBSCAN+envelope+erosion
  search), `analyse_clusters` (per-cluster summary: size, shape, composition),
  and `run_relabel_sweep` (real run + N mc-shuffled control runs at one
  parameter combination).
- **`sweep.py`** - `run_parameter_sweep` (the order/dmax outer loop),
  `optimise_n_min`, and the real-vs-random / composition plotting functions.
- **`fast_sweep.py`** - the accelerated sweep (see below):
  `run_fast_parameter_sweep` (drop-in for `run_parameter_sweep`),
  `run_fast_grid_sweep` (both axes at once), plus `select_sweep` (one
  row/column of a grid), `summarise_grid` (the whole plane, one row per cell),
  `load_sweep_size_table` and `validate_fast_grid`.
- **`corrections.py`** - edge-cluster detection (alpha-shape boundary),
  Fm-peak deconvolution, per-cluster nearest-neighbour distance, and
  `get_cluster_metrics` (volume fraction, number density, etc.).
- **`viz.py`** - `plot_size_distribution` (violin + every cluster as its own
  hoverable point) and 3D visualisation of selected clusters. Every plot's
  tooltip leads with the cluster's `cluster_id`, so a cluster picked out of a
  size distribution can be typed straight into `show_clusters_in_pos`.
- **`parity.py`** - compare a run against real posgen output
  (`load_posgen_clustering`, `compare_clusterings`, `variant_matrix`). Use it
  whenever you change the clustering code:

  ```python
  old = pgp.load_posgen_clustering(cluster_pos, clusterid_pos, rrng_df)
  new = pgp.run_engine_variant(pos_df, rrng_df, params)
  pgp.compare_clusterings(old, new).summary()
  ```

## The accelerated sweep

`run_parameter_sweep` runs a full cluster search for every (parameter value,
dataset) pair: re-range the pos file, run DBSCAN, grow the envelope one cluster
at a time, erode it. At the 5-10 parameter values and 20-50 randomised control
datasets a real study needs, that is 100-500 full pipelines, and hours.

### Making an existing notebook fast

`run_fast_parameter_sweep` has `run_parameter_sweep`'s signature, returns the
same manifest and writes the same file layout, so a notebook that already works
changes one name per sweep and nothing else - `optimise_n_min`, both
real-vs-random graphs and the composition graphs all keep taking the result
directly:

```python
manifest = pgp.run_fast_parameter_sweep(      # was pgp.run_parameter_sweep
    pos_df, rrng_df, core_ions, bulk_ions,
    swept_param="dmax", swept_values=[0.3, 0.4, 0.5, 0.6],
    fixed_order=7, n_min=2, n_relabel_runs=20,
    output_dir="my_output/", prefix="mytip",
    n_jobs=-1,                                # the one addition worth making
)
```

Two things to know before relying on the numbers:

- `size` in the per-cluster files counts core+border atoms rather than every
  atom, and only core-species ions get a column. `optimise_n_min` and both
  real-vs-random graphs are unaffected - they threshold `size_core`. The
  *absolute* composition graph changes meaning: its denominator becomes the core
  ions instead of every atom in the cluster.
- `envelope_radius`, `d_erosion` and `dclassify` are accepted for signature
  compatibility but describe stages this path skips, so passing them warns.

`run_fast_grid_sweep` is the same machinery with both axes swept at once; the
1D form above is a grid one value wide, which is why they share everything below.

`run_fast_grid_sweep` turns the loops inside out. **Per dataset** it computes
once:

- the distance from every core-species atom to its k-th nearest core-species
  neighbour, for every k in the swept orders (one kNN query at `k = max(order)`);
- every pair of core-species atoms closer than `max(dmax)` (one radius query),
  kept sorted by distance.

Those answer the core-classification and linking steps at *any* cell of the
grid: an atom is a core point at `(order, dmax)` exactly when its order-th
neighbour distance is `< dmax`, and two core points link exactly when their
pair distance is `< dmax` - a prefix of the sorted pair list. Both are monotone
in dmax, so moving to the next parameter value only ever *adds* links, and
asking "would this atom be reached at the next dmax?" is a re-threshold of
arrays already in memory rather than another walk of the tree. Cost scales with
the number of datasets, not datasets x parameters.

### Orders are not treated incrementally

Only the dmax axis is walked monotonically, and only in the direction where that
is sound: raising dmax adds core points and adds links, so DBSCAN clusters can
merge but never split. Raising the **order** does the opposite - it takes core
points away, and a cluster held together by a thin neck of core atoms falls into
two when that neck stops passing the knn test.

So no partition is carried between cells. `CoreNeighbourhood.labels`
re-derives which atoms are core points and re-runs connected components from
scratch at every cell, over whatever subgraph that leaves. A split at higher
order comes out as two clusters, exactly as `run_clustering` reports it; there
is no incremental union that could hold them together. What is shared across
cells is the *evidence* - neighbour distances and candidate pairs - never a
result. `test_raising_the_order_splits_a_cluster` pins this on a dumbbell that
is one cluster at order 3 and two at order 10, against the full engine.

### Relabelling is the standard one

The control datasets are the ordinary global relabelling - shuffle the pos
file's `mc` column, re-range every atom - not a substitute for it. It is just
computed without walking the mass ranges twice. Ranging is elementwise in `mc`,
so ranging a permuted `mc` column and permuting the ranged labels give the same
array:

```
R(mc[pi])[i] = R(mc[pi[i]]) = R(mc)[pi[i]] = (R(mc)[pi])[i]
```

so a control dataset is a shuffle of the precomputed per-atom label array (one
`int8` per atom - 1.9 MB on this tip). Re-ranging, not shuffling, is where the
cost was. Per control dataset on `MNSP tip C`:

| | per dataset |
|---|---|
| shuffle the label array (what this does) | **21 ms** |
| shuffle `mc`, then re-range | 80 ms |
| copy the frame + `assign_rrng_to_pos` (what the slow path does) | 281 ms |

No statistical argument is needed for any of this - the datasets are identical,
permutation for permutation, which
`test_relabelling_equals_shuffling_mc_and_reranging` asserts directly.

Pass the same `cached_shuffles` to `run_fast_grid_sweep` and to
`cluster.run_relabel_sweep` and both paths cluster the *same* control datasets,
which is how their control runs are compared rather than merely their
distributions (`test_control_runs_match_the_full_pipeline`; on synthetic tips
the size distributions come out equal outright).

### One further shortcut

- **No envelope, no erosion.** Clusters here are core+border atoms only. The
  sweep stage never needed more: `optimise_n_min` and both real-vs-random plots
  threshold posgen's `Cluster Size(core)`, and the envelope only ever *adds*
  bulk-tree atoms while erosion only ever removes them - core-species atoms are
  never in the bulk tree, so neither stage can change that number.

Measured on `MNSP tip C` (1.85M atoms, 93,839 core-species atoms), at
the scale a real study runs at - 10 dmax values x (1 real + 50 randomised)
datasets, i.e. 510 cluster searches:

| | time |
|---|---|
| `run_parameter_sweep` | ~27 min |
| `run_fast_grid_sweep`, `n_jobs=1` | 15 s |
| `run_fast_grid_sweep`, `n_jobs=-1` (8 cores) | **8 s** |

The `run_parameter_sweep` figure is `run_clustering` timed on this tip at each
of those dmax values (1.0 s at dmax 0.3 rising to 7.0 s at 0.75) plus its
per-dataset relabelling, scaled to 51 datasets. Per cluster search the ratio is
250-1200x. Widening it to the full 5-order x 10-dmax grid - 2,550 searches -
takes 55 s, because the extra orders and dmax values cost almost nothing once
each dataset's neighbourhood exists.

### How exact is it

`validate_fast_grid` runs both paths on a real tip and compares them cell by
cell. Across that 25-cell grid, every cell had:

- identical core-species atom sets;
- identical core-point classification (verified directly against
  `DBSCAN(eps=dmax, min_samples=order)`);
- an identical set of clustered atoms;
- every atom that changed cluster confined to the **ambiguous border atoms** -
  those within dmax of core atoms of two different clusters, which any correct
  DBSCAN may give to either. This path gives them to the nearest core atom;
  sklearn gives them to whichever its traversal reaches first. Typically half
  the ambiguous atoms land differently (e.g. 30 of 59 at order 7, dmax 0.5, out
  of 15,573 clustered atoms). This is the same ambiguity that already puts this
  engine 1 atom out of 45,194 apart from posgen itself.

For the number that decides anything - clusters of at least n_min core+border
atoms, which is what `optimise_n_min` and the real-vs-random curves are built
from - the largest disagreement anywhere on the n_min axis is 1.6-8.8%
depending on the cell. Those maxima all sit at **small n_min (4-9)**, where the
tip has thousands of two- and three-atom clusters and one border atom moving
creates or destroys a whole cluster.

At the n_min values that actually get selected - 19 to 62 on this tip - the two
paths differ by 0 to 2 clusters out of 55-65, under 2%; at order 7, dmax 0.5 the
curves are identical at n_min 20 and 30, and 1 cluster apart at 25, 40 and 50.
Whichever path picks the parameters, it picks the same ones.

What it is *not* is a substitute for `run_clustering` at the final chosen
parameters: cluster composition, volume and volume fraction all depend on the
envelope and erosion. Run the full pipeline once, at the end, on the parameters
this picks - which is what the example notebook does.

## Memory discipline for relabelling

Because Python (not posgen's C++ internals) now does the relabelling itself,
`run_relabel_sweep` never holds more than one relabel run's full per-atom
data in memory at a time: only the real run's full labelled atoms are kept,
each of the N randomised runs is reduced to a small per-cluster size table and
written to its own parquet file (`{run_name}_random/run_NNN.parquet`) before
the next run starts. A KDTree over all atom positions is built once per pos
file and reused across every swept parameter value and every relabel run,
since atom positions never change across relabelling - only which atoms count
as which species does.

## Terms

**dmax** - the maximum-separation cluster-linking distance (nm).

**order** - the core-point order statistic, posgen's `knn`: a point is a core
point if `order` points *including itself* (among core species) lie within
`dmax`. Same number as in a posgen run's `_order_N` filename.

**n_min** - the minimum allowable cluster size; swept as a diagnostic axis and
then auto-selected via a real-vs-random cluster-count threshold (95% "true
positive ratio" by default).

**relabelling / randomisation** - shuffling the pos file's `mc` column and
re-clustering, to build a random-cluster-count baseline for the same
dataset and parameters.

## Authors

posgenpy was created by Dan Hayley, Andy London, James Famelton, Hazel
Gardner, and Przemek Klups. This pure-Python engine replaces the original
posgen/libatomprobe-based implementation.
