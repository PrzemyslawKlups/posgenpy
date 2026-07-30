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

# sweep dmax at a fixed order, with 20 relabelled control runs per value
manifest = pgp.run_parameter_sweep(
    pos_df, rrng_df, core_ions, bulk_ions,
    swept_param="dmax", swept_values=[0.3, 0.4, 0.5, 0.6],
    fixed_order=7, n_relabel_runs=20,   # posgen's Order, as in _order_7
    output_dir="my_output/", prefix="mytip",
)

n_min_values = list(range(2, 21, 2))
pgp.plot_real_cluster_ratio_across_swept_param(manifest, "dmax", n_min_values)

optimal_n_min, real_clusters, random_clusters = pgp.optimise_n_min(manifest, n_min_values=range(1, 200))
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
- **`corrections.py`** - edge-cluster detection (alpha-shape boundary),
  Fm-peak deconvolution, per-cluster nearest-neighbour distance, and
  `get_cluster_metrics` (volume fraction, number density, etc.).
- **`viz.py`** - 3D visualisation of selected clusters.
- **`parity.py`** - compare a run against real posgen output
  (`load_posgen_clustering`, `compare_clusterings`, `variant_matrix`). Use it
  whenever you change the clustering code:

  ```python
  old = pgp.load_posgen_clustering(cluster_pos, clusterid_pos, rrng_df)
  new = pgp.run_engine_variant(pos_df, rrng_df, params)
  pgp.compare_clusterings(old, new).summary()
  ```

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
