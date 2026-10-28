# PosgenPy (pure-Python engine)

A dependency-free, pure-Python replacement for the `posgen` C++ cluster-search
tool. Only `numpy`/`pandas`/`scikit-learn`/`scipy`.

It implements the same "maximum separation" cluster-search algorithm posgen
used (DBSCAN-style core classification + envelope + erosion), plus the
relabelling/randomisation workflow used to estimate a real-vs-random cluster
threshold and pick `n_min` - all driven from ordinary Jupyter notebooks.

If you're looking for the old posgen/XML-driven workflow (`write_xml_with_relabelling`,
`subprocess.call(["./posgen", ...])`), it's kept, unmodified, under
[`legacy/`](legacy/) for reference/cross-validation only. It is not imported
by anything in this package. The posgen-driven version described in the paper
is tagged [`v1.0.0`](https://github.com/PrzemyslawKlups/posgenpy/tree/v1.0.0),
and [`v1.1.0`](https://github.com/PrzemyslawKlups/posgenpy/tree/v1.1.0) adds
the cluster corrections.

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

# the whole plane at once: one row per (order, dmax) - optimal n_min, cluster
# counts, the largest cluster, and how close the closest pair of clusters got
summary = pgp.summarise_grid(grid)
pgp.plot_grid_heatmap(summary, value="optimal_n_min", colorscale="Purples").show()

# every real cluster of every cell, clipped at that cell's own n_min, with its
# nearest neighbour measured inside its own cell
grid_clusters = pgp.load_grid_cluster_table(grid, n_min=summary)
pgp.plot_distribution_grid(grid_clusters, value="size_core", log_x=True).show()

# and what happened between neighbouring cells: raising dmax merges clusters,
# raising the order splits them - where either takes off is the elbow
splits = pgp.track_cluster_splits(grid, axis="order", n_min=summary)
pgp.plot_split_tracking(splits, value="n_new_from_splits").show()

# ...or read the grid along either axis - these take the same manifest shape
# run_parameter_sweep returns, so every sweep function below works unchanged
dmax_sweep = pgp.select_sweep(grid, order=7)

n_min_values = list(range(2, 21, 2))
pgp.plot_real_cluster_ratio_across_swept_param(dmax_sweep, "dmax", n_min_values)

optimal_n_min, real_clusters, random_clusters = pgp.optimise_n_min(dmax_sweep, n_min_values=range(1, 200))

# ...and check the answer means anything: the n_min search has three exits, and
# a randomised control whose solute network percolates satisfies the threshold
# without measuring a thing (see "When the randomised control percolates")
pgp.n_min_report(dmax_sweep, n_min_values=range(1, 200))
pgp.percolation_report(dmax_sweep)

# every cluster of every cell as one row, for hoverable size distributions
sizes = pgp.load_sweep_size_table(dmax_sweep)
pgp.plot_size_distribution(sizes, x="dmax", y="size_core").show()
```

Two example notebooks, same engine, different way in:

- [`PosgenPy_Full_Analysis_Example_PythonOnly.ipynb`](PosgenPy_Full_Analysis_Example_PythonOnly.ipynb)
  — one parameter at a time: order sweep, dmax sweep, n_min optimisation, final
  cluster search + corrections (edge-cluster removal, Fm-peak deconvolution,
  per-tip metrics), composition/violin/3D plots, and a CV01
  parameter-sensitivity check.
- [`PosgenPy_Grid_Analysis_Example_PythonOnly.ipynb`](PosgenPy_Grid_Analysis_Example_PythonOnly.ipynb)
  — the whole (order x dmax) plane first: six heatmaps of it, size and
  nearest-neighbour-distance distributions as ridgelines over both axes,
  split/merge tracking between neighbouring cells, individual clusters in 3D
  and across an axis, then the same final run with its clusters drawn as
  see-through alpha shapes.

## Installation

```bash
uv pip install -e .
```

(or plain `pip install -e .` - see `pyproject.toml` at the repo root for the
dependency list. `pyarrow` is required for the parquet files this package
reads/writes.)

## Matching posgen

The defaults reproduce posgen's `maxsep` algorithm. Verified atom-for-atom
against five real posgen runs (see [tests/README.md](tests/README.md)):

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

## When the randomised control percolates

`optimise_n_min` compares **counts of clusters** — the true-positive ratio is
`(N_real − ⟨N_rand⟩) / N_real`. A percolating network contributes exactly
**one** cluster no matter how much solute it swallowed. So at a `dmax` large
enough that the solute network of the *randomised* tip spans the cell,
`⟨N_rand⟩` collapses towards 1, the ratio sails past the threshold, and a
confident-looking `n_min` comes back from a parameter region with no
discriminating power at all. On one MD dump it made a strictly more permissive
order look strictly more discriminating — orders 1–3 reporting a *lower* n_min
than order 5, which cannot be true.

That is a parameter choice, not a bug in the search: the giant cluster is the
correct output of a maximum-separation search at those parameters. What was
missing is anything that noticed.

```python
pgp.percolation_report(order_sweep)     # one row per swept value
pgp.n_min_report(order_sweep, n_min_values=range(1, 200))
```

`percolation_report` reads the largest cluster of each *control* run out of the
per-cluster size tables the sweep already wrote — so it **runs on existing
sweeps with nothing re-run** — and divides it by the tip's core-species atoms
(`n_core_atoms` in the manifest; older manifests fall back to the solute the
control actually clustered, which is labelled as such because it inflates the
fraction at high order). The `verdict` is `"ok"` / `"near-threshold"` /
`"percolating"`, cut at 0.10 / 0.30 by default. Measured regimes separate
cleanly around those — sub-critical ≤ 0.01, critical ≈ 0.21, percolating
0.65–0.88 — but they are keyword arguments, not physics.

`n_min_report` is `optimise_n_min` with its reasoning shown. The search has
three exits and the returned integer used to be the same shape for all three:

| `status` | |
|---|---|
| `"reached"` | the threshold was genuinely crossed — the only one that is an answer |
| `"unreachable"` | the real clusters ran out first; no n_min in *any* range can get there |
| `"range_exhausted"` | no n_min in `n_min_values` got there |

Only a `"reached"` row carries an `n_min`; the other two are `NaN`. It also
carries `true_positive_ratio_atoms` — the same real-vs-control comparison
counted in **atoms** rather than clusters, over every cluster in the saved
tables. A percolating control cannot hide inside it as a single cluster, so
where it sits at or below zero the control clusters as much solute as the tip
does and nothing chosen there means anything.
`plot_atom_ratio_across_swept_param` draws it.

Both failure modes reach a caller that is not reading either report:
`optimise_n_min` now raises a `UserWarning` on the two failure exits (its
`verbose=False` is exactly how `summarise_grid` calls it, so a printed message
was no message at all); `summarise_grid` takes `n_min_reached` from the status
instead of re-deriving it from the returned integer, and adds `n_min_status`,
`percolation_verdict`, `giant_fraction` and `true_positive_ratio_atoms`
columns; and the two real-vs-random sweep plots shade the swept values whose
control has percolated.

Nothing about the cluster search itself changed, and `optimise_n_min` keeps its
signature and 3-tuple return.

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
  `load_sweep_size_table`, `load_grid_cluster_table` (every real cluster of
  every cell, clipped at that cell's n_min), `track_cluster_splits` (what
  became of what between neighbouring cells) and `validate_fast_grid`.
- **`corrections.py`** - edge-cluster detection (alpha-shape boundary),
  Fm-peak deconvolution, per-cluster nearest-neighbour distance (raw and
  normalised by the two clusters' radii), per-cluster alpha-shape meshes, and
  `get_cluster_metrics` (volume fraction, number density, etc.).
- **`viz.py`** - `plot_size_distribution` (violin + every cluster as its own
  hoverable point), `plot_distribution_grid` (the same distributions as
  ridgelines over both grid axes at once), `plot_grid_heatmap`,
  `plot_split_tracking`, `plot_cluster_alpha_shapes` (see-through cluster
  surfaces, edge clusters apart) and 3D visualisation of selected clusters -
  one at a time, or across a whole swept axis. Every plot's tooltip leads with
  the cluster's `cluster_id`, so a cluster picked out of a size distribution
  can be typed straight into `show_clusters_in_pos`; tooltips otherwise show
  `DEFAULT_HOVER_COLUMNS` rather than every column a stats table carries
  (`hover="all"` for the lot, or a list of your own).
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

## Reading a finished grid

Three functions read the plane the sweep wrote, and between them they are what
the grid notebook is built from.

**`summarise_grid`** — one row per cell. Beyond the n_min columns it carries
`n_clusters_all` (before n_min), `max_cluster_size`, and how close the closest
pair of clusters got, both in nm (`min_cluster_distance`) and normalised by the
two clusters' radii (`min_cluster_distance_norm`). All of them are measured
over the clusters that survive **that cell's own optimal n_min**, since the
ones below it are exactly the noise the n_min was chosen to remove. Every
column is a heatmap: `plot_grid_heatmap(summary, value=..., colorscale=...)`.

The normalised distance is `min_j d_ij / (R_i + R_j)`, so **below 1 the two
clusters' equivalent spheres overlap** — one precipitate reported as two, or
two about to merge at the next dmax. It is exact rather than a k-nearest
approximation: a large cluster past the k-th neighbour can crowd a small one
more than its nearest neighbour does, so anything that could still win gets a
second, wider search (`corrections.find_knn_cluster_distance`).

**`load_grid_cluster_table`** — every real cluster of every cell as one row,
clipped at each cell's own n_min (pass the `summarise_grid` frame and it uses
`optimal_n_min` per cell), with the nearest-neighbour distance measured *within
its own cell*, because a cluster's neighbours are different at every parameter
pair. `viz.plot_distribution_grid` draws it as a ridgeline — one subplot row
per order, one filled density curve per dmax, with every cluster still drawn as
its own hoverable point.

**`track_cluster_splits`** — every cell clusters the *same atoms*, so
consecutive cells can be compared atom by atom rather than by cluster counts.
That is the difference between "5 fewer clusters because five precipitates
merged into two" and "because five dissolved below n_min": same count,
different clustering. Per step it reports `n_split`, `n_new_from_splits`,
`n_merged`, `n_lost_to_merges`, `n_dissolved`, `n_created`, and *which*
clusters split (`split_cluster_ids`, ids in the earlier cell, so they can be
plotted).

The two axes do opposite things, and that is a property of the algorithm rather
than of any dataset: **raising dmax adds links**, so clusters merge but never
split; **raising the order removes core points**, so a cluster held together by
a thin neck of core atoms falls into two. Where either count takes off is the
elbow — past it the parameter is not refining the clustering, it is
rearranging it. `viz.plot_split_tracking` plots that, and since each row also
carries the order/dmax of the cell it *arrives at*, `plot_grid_heatmap` takes
the same table.

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

**normalised cluster separation** - the distance between two cluster centres
divided by the sum of their radii (`equivalent_spherical_radius` by default).
Below 1 the two clusters' equivalent spheres overlap; above ~2 they are clearly
separate. Written `closest_cluster_distance_norm` per cluster and
`min_cluster_distance_norm` per grid cell.

**alpha shape** - the surface traced by rolling a probe sphere of radius
`1/alpha` over a cluster's atoms. Unlike a convex hull it keeps concavities, so
a dumbbell that should have been two clusters looks like one.

## Citing

If you use posgenpy, please cite:

> P. Klupś, D. Haley, A. J. London, H. Gardner, J. Famelton, B. M. Jenkins,
> J. M. Hyde, P. A. J. Bagot and M. P. Moody, "PosgenPy: An Automated and
> Reproducible Approach to Assessing the Validity of Cluster Search Parameters
> in Atom Probe Tomography Datasets", *Microscopy and Microanalysis* **28**(4),
> 1066–1075 (2022). [doi:10.1017/S1431927621012368](https://doi.org/10.1017/S1431927621012368)

The same reference is in [`CITATION.cff`](CITATION.cff), which GitHub's "Cite
this repository" button reads.

## Authors

posgenpy was created by Przemek Klups, Dan Haley, Andy London, James Famelton
and Hazel Gardner. This pure-Python engine, co-created with Claude Code, replaces the original posgen/libatomprobe-based implementation.

## Contact

Bug reports and questions are welcome as [GitHub issues](https://github.com/PrzemyslawKlups/posgenpy/issues).
For anything else, email <p.klups@imperial.ac.uk>.
