"""
monocle2_py: Pure-Python Monocle 2 for single-cell trajectory analysis.

A standalone mirror of the canonical implementation that lives in
``omicverse.single.Monocle`` / ``omicverse.external.monocle2_py``. This
repo exists for users who want Monocle 2 without pulling in the full
omicverse stack — all algorithmic work is developed upstream in
omicverse and synced here.

Input: AnnData objects. All analysis state is written to
``adata.obs`` / ``adata.var`` / ``adata.uns['monocle']`` / ``adata.obsm``
so the annotated object stays compatible with the scanpy ecosystem.

>>> from monocle2_py import Monocle
>>> mono = Monocle(adata)
>>> (mono.preprocess()
...      .select_ordering_genes()
...      .reduce_dimension()          # method='fast' by default
...      .order_cells())
>>> de = mono.differential_gene_test()
>>> beam = mono.BEAM(branch_point=1)
"""

from .core import (
    set_ordering_filter,
    detect_genes,
    estimate_size_factors,
    estimate_dispersions,
    dispersion_table,
    estimate_t,
    relative2abs,
)
from .dimension_reduction import reduce_dimension
from .ordering import order_cells
from .differential import (
    differential_gene_test,
    BEAM,
    fit_model,
    gen_smooth_curves,
)
from .clustering import cluster_cells, cluster_genes
from .plotting import (
    plot_cell_trajectory,
    plot_trajectory_overlay,
    plot_genes_in_pseudotime,
    plot_genes_branched_heatmap,
    plot_genes_branched_pseudotime,
    plot_cell_clusters,
    plot_genes_jitter,
    plot_genes_violin,
    plot_ordering_genes,
    plot_pseudotime_heatmap,
    plot_complex_cell_trajectory,
    plot_multiple_branches_pseudotime,
    plot_multiple_branches_heatmap,
    plot_rho_delta,
    plot_pc_variance_explained,
)
from .ddrtree import DDRTree
from .utils import cal_ABCs, cal_ILRs

from .monocle import Monocle
# Don't leak the submodule name into the public namespace — only the class.
del monocle  # noqa: F821 — the name exists because Python binds submodules

__version__ = "0.1.0"
__all__ = [
    "Monocle",
    # core / preprocessing
    "set_ordering_filter",
    "detect_genes",
    "estimate_size_factors",
    "estimate_dispersions",
    "dispersion_table",
    "estimate_t",
    "relative2abs",
    # dimension reduction & ordering
    "reduce_dimension",
    "order_cells",
    # differential expression
    "differential_gene_test",
    "BEAM",
    "fit_model",
    "gen_smooth_curves",
    # clustering
    "cluster_cells",
    "cluster_genes",
    # plotting
    "plot_cell_trajectory",
    "plot_trajectory_overlay",
    "plot_genes_in_pseudotime",
    "plot_genes_branched_heatmap",
    "plot_genes_branched_pseudotime",
    "plot_cell_clusters",
    "plot_genes_jitter",
    "plot_genes_violin",
    "plot_ordering_genes",
    "plot_pseudotime_heatmap",
    "plot_complex_cell_trajectory",
    "plot_multiple_branches_pseudotime",
    "plot_multiple_branches_heatmap",
    "plot_rho_delta",
    "plot_pc_variance_explained",
    # core algorithms & utilities
    "DDRTree",
    "cal_ABCs",
    "cal_ILRs",
]
