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
    fixed_order=5, n_relabel_runs=20,
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

**order (k)** - the core-point order statistic: a point counts as a core
point if its k-th nearest neighbour (among core species) is within `dmax`.

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
