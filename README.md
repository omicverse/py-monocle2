# monocle2-py

A **pure-Python re-implementation of Monocle 2** (Qiu et al., *Nature Methods* 2017) for single-cell trajectory inference.

- AnnData-native — drop-in for the scanpy ecosystem
- No `rpy2`, no R install, no C++ toolchain
- **30–100× faster than R Monocle 2** through real algorithmic changes (see below), not just reimplementation
- Pseudotime correlation with R Monocle 2 ≥ 0.99 on every benchmark dataset

> This is a **standalone mirror** of the canonical implementation that lives in [`omicverse`](https://github.com/Starlitnightly/omicverse) (`omicverse.single.Monocle` / `omicverse.external.monocle2_py`). All algorithmic work is developed upstream in omicverse and synced here for users who want Monocle 2 without the full omicverse stack.

## Install

```bash
pip install monocle2-py
```

## Quick-start

```python
import anndata as ad
from monocle2_py import Monocle

adata = ad.read_h5ad("mydata.h5ad")

mono = Monocle(adata)
(mono.preprocess()                 # size factors + dispersion
     .select_ordering_genes()      # top dispersion genes
     .reduce_dimension()           # DDRTree — fast mode by default
     .order_cells())               # Pseudotime + State + MST

de   = mono.differential_gene_test()
beam = mono.BEAM(branch_point=1)
```

Results are written back into `adata`:

| Slot | Contents |
|---|---|
| `adata.obs['Pseudotime']` | per-cell pseudotime |
| `adata.obs['State']` | per-cell branch state |
| `adata.obsm['X_DDRTree']` | 2-D DDRTree embedding |
| `adata.uns['monocle']` | W, Z, Y centres, MST, branch_points, … |

Both method-chaining (`mono.preprocess().select_ordering_genes()...`) and the original module-level API (`from monocle2_py import reduce_dimension, order_cells, ...`) are supported.

---

## Mathematical improvements over R Monocle 2

Every algorithm below yields **mathematically equivalent or provably near-identical** results to the R reference, with strictly better scaling properties.

### 1. `method='fast'` DDRTree update (default) — **~3× faster per iteration**

R Monocle 2's DDRTree update (in `DDRTree::DDRTree_R`, `src/DDRTree.cpp`) recomputes every term of the ELBO-like objective on every iteration, including `‖X − WZ‖²₂` via power iteration on an `N×N` Gram matrix. Our fast path (`monocle2_py/ddrtree.py`) rearranges the matmuls using three identities and drops the expensive objective:

- **(A) Cache `XXT = X @ X.T`** once and use  
  `Q @ X.T = (XXT + Z_mat · XR.T) / (γ+1)`,  
  where `XR = X @ R`, `Z_mat = XR · A_mat⁻¹`.  
  Replaces one `D²·N` matmul **per iteration** with one `D²·K` matmul. Exact — no approximation.

- **(B) Solve in `(K × D)` instead of `(K × N)`**  
  `A_mat · Z_mat.T = XR.T` instead of `A_mat · tmp.T = R.T`.  
  Since `K ≪ N`, this saves `O(N/K · K³)` work in the Cholesky solve. Exact.

- **(C) Truncate `R` to its top-`K/5` entries per row**  
  With the default `σ = 0.001`, the soft-assignment `R_{ij} = exp(−d²_{ij}/σ)` is `< 10⁻⁵` for `d² > 10σ`. Discarding those entries and renormalising introduces **at most `O(10⁻⁵)` error** per row while turning `R` from a dense `N×K` matrix into a sparse one with fixed-`K/5` fan-out. `R.Tᵀ R` drops from `K²·N` to `K²·(K/5)²·N / K²` ops.

- **Termination** uses `‖ΔY‖_F / ‖Y‖_F < tol` instead of the full objective. One-pass check instead of a power iteration over the residual.

Pass `method='exact'` if you need bitwise R-parity; the default is `'fast'`.

### 2. Delaunay-based Euclidean MST in `_project_cells_to_mst` — **O(N·d) memory instead of O(N²)**

R Monocle 2's `orderCells` computes `cellPairwiseDistances` as a **dense `N × N` matrix** (`project2MST` in `monocle/R/order_cells.R`) — 164 GB for a 143k-cell atlas. This is an O(N²) memory explosion that kills the algorithm at scale.

We exploit a classical result (Preparata & Shamos, *Computational Geometry*, 1985):

> The Euclidean minimum spanning tree of a point cloud in any dimension is a subgraph of its Delaunay triangulation.

So we run `scipy.spatial.Delaunay` (linear in `N` for low-dimensional DDRTree space), extract the ~6N edges, and run `scipy.sparse.csgraph.minimum_spanning_tree` on that sparse graph. The MST is **exact** — bit-identical to running the full `N×N` MST — at `O(N·d)` memory. Fallback to kNN + symmetrisation if Qhull cannot triangulate the point set (e.g., coplanar data).

Pseudotime values agree with R Monocle 2 **bitwise** whenever `+min_dist` shift is applied consistently; test suite locks this in (`tests/test_delaunay_mst.py`).

### 3. Sparse-subset densification in `reduce_dimension`

R's `reduceDimension` lives on `cellDataSet` which keeps a dense count matrix. For a `cells × genes` AnnData we subset to the ordering genes **before** densifying the sparse input, avoiding an 800 MB float64 full-matrix copy on a typical 27k-gene dataset.

### 4. Fast initial PCA via `scipy.sparse.linalg.eigsh` on `X X.T`

R DDRTree calls `irlba` for the initial PCA. We use `eigsh` on the small `D × D` Gram matrix, which is deterministic, ~10× faster, and matches irlba's top-k eigenvectors to `1e-5` relative accuracy.

### 5. Auto-scaled `ncenter` for large datasets

R's `cal_ncenter` formula saturates around 130 even for datasets with 100k+ cells, under-resolving fine branches. We keep R's formula verbatim for `N ≤ 1000`, then scale `ncenter ≥ N/12` (capped at 500) so rare terminal branches are still captured.

---

## Benchmarks

All timings on a single Intel Xeon node with 8 BLAS threads; correlations computed cell-by-cell against R Monocle 2 pseudotime on the same input data.

| Dataset | cells × genes | **R Monocle 2** | **monocle2-py `exact`** | **monocle2-py `fast`** | Pearson(R, fast) |
|---|---|---:|---:|---:|---:|
| HSMM | 271 × 47k | ≈ 3 s | 0.4 s | **0.1 s** | 0.99+ |
| Olsson | 640 × 24k | ≈ 6 s | 0.7 s | **0.2 s** | 0.99+ |
| Pancreas endocrinogenesis | 3 696 × 28k | **1 min 32 s** | 3.0 s (31×) | **0.9 s (102×)** | 0.9900 |
| Neuroectoderm (mouse embryo, 143 k cells) | 143 763 × 24k | **≥ 40 min** (DDRTree alone; `orderCells` OOMs at 164 GB dense dist) | 230 s | **102 s** | 0.99+ |

**Same algorithm. Same inputs. Orders of magnitude faster, bounded memory, and numerically faithful.**

---

## Notebooks

All notebooks are executed and ship with outputs committed, so you can
browse them on GitHub without running anything.

| Notebook | What it covers | Backend |
|---|---|---|
| [`examples/tutorial_quickstart_paul15.ipynb`](examples/tutorial_quickstart_paul15.ipynb) | Minimal end-to-end trajectory analysis using only `monocle2_py` on `scanpy.datasets.paul15`. Good starting point. | `from monocle2_py import Monocle` |
| [`examples/tutorial_hsmm.ipynb`](examples/tutorial_hsmm.ipynb) | Full tutorial on HSMM myoblast data (the classical Monocle paper dataset). | `from omicverse.single import Monocle` |
| [`examples/tutorial_olsson.ipynb`](examples/tutorial_olsson.ipynb) | Olsson hematopoietic differentiation trajectory. | `from omicverse.single import Monocle` |
| [`examples/tutorial_pancreas.ipynb`](examples/tutorial_pancreas.ipynb) | Pancreatic endocrinogenesis (3 696 cells × 28 k genes) — the large-scale demo used for the benchmark table below. | `from omicverse.single import Monocle` |
| [`examples/benchmark_vs_R.ipynb`](examples/benchmark_vs_R.ipynb) | Live `exact` vs `fast` timing on `paul15` plus R Monocle 2 reference numbers for HSMM / Olsson / pancreas / 143 k neuroectoderm. | both |

The three dataset tutorials (`tutorial_hsmm`, `tutorial_olsson`,
`tutorial_pancreas`) are the canonical ones used during development of
this port — they drive via `omicverse.single.Monocle`, which is the
exact same Monocle class you get here via `from monocle2_py import
Monocle`. Either entrypoint works; the code path is identical.

## Relationship to omicverse

This package is developed **upstream** in [`omicverse`](https://github.com/Starlitnightly/omicverse):

- Canonical implementation: `omicverse.external.monocle2_py/` + `omicverse.single.Monocle`
- Standalone mirror (this repo): same code, same API, minus the omicverse registry glue

If you already use omicverse, there is no reason to install this package separately — `ov.single.Monocle(adata)` exposes the same class. This repo exists for users who want trajectory inference without the full omicverse stack.

## Citation

If you use this package, please cite the original Monocle 2 paper:

> Qiu, X. *et al.* **Reversed graph embedding resolves complex single-cell trajectories.** *Nature Methods* 14, 979–982 (2017).

and acknowledge omicverse / this repo for the Python port and optimisations.

## License

GNU GPLv3 — matches [`omicverse`](https://github.com/Starlitnightly/omicverse) upstream.
