"""
DDRTree: Discriminative Dimensionality Reduction via Learning a Tree.

Pure Python/NumPy/SciPy implementation of the DDRTree algorithm,
faithfully ported from the R/C++ DDRTree package.
"""

import numpy as np
from scipy import sparse as _sparse
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.linalg import cho_factor, cho_solve

# Optional numba acceleration for Hartigan-Wong inner loops. The patch must
# work without numba; if it is not installed the @njit decorator becomes a
# no-op and the same NumPy code runs unaccelerated (slower but identical).
try:
    from numba import njit  # type: ignore
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):  # noqa: D401
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def _wrap(f):
            return f
        return _wrap


def _sqdist(a, b):
    """Compute squared distance matrix between columns of a and b.
    a: (D, N), b: (D, M) -> returns (N, M)
    """
    aa = np.sum(a ** 2, axis=0)
    bb = np.sum(b ** 2, axis=0)
    ab = a.T @ b
    dist = np.abs(aa[:, None] + bb[None, :] - 2 * ab)
    return dist


def _pca_projection_irlba_like(C, L):
    """Iterative top-L eigenvectors of symmetric matrix C — matches R's irlba.

    R's DDRTree C++ code calls pca_projection_R (which uses irlba) inside
    every iteration. The approximation introduces a small random-direction
    perturbation that helps DDRTree escape smooth local minima and
    discover branching structure that exact eigendecomposition can miss.
    """
    from scipy.stats import norm as _norm
    D = C.shape[0]
    if L >= min(C.shape):
        eigenvalues, eigenvectors = np.linalg.eigh(C)
        # Sort by |eigenvalue| descending. Algebraic-sort + which='LM' below
        # mismatch if sym_mat becomes indefinite (issue #14). Defensive only —
        # HSMM has sym_mat PSD on every iter (verified Subagent C).
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        return eigenvectors[:, idx[:L]]
    v0 = _norm.ppf(np.arange(1, D + 1) / (D + 1))
    v0 = v0 / np.linalg.norm(v0)
    try:
        eigenvalues, eigenvectors = eigsh(C, k=L, which='LM', v0=v0,
                                          maxiter=1000, tol=1e-5)
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        return eigenvectors[:, idx]
    except Exception:
        eigenvalues, eigenvectors = np.linalg.eigh(C)
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        return eigenvectors[:, idx[:L]]


def _pca_projection(C, L):
    """PCA projection: return top-L eigenvectors of symmetric matrix C.

    Uses exact np.linalg.eigh when C is small enough (more accurate than
    iterative eigsh), otherwise falls back to eigsh with R-style initial
    vector to match R's irlba.

    C: (D, D), L: int -> returns (D, L)
    """
    from scipy.stats import norm as _norm

    D = C.shape[0]

    # Use exact eigendecomposition for matrices that fit in memory
    # Memory: D*D*8 bytes, 4000x4000 = 128MB is fine
    if L >= min(C.shape) or D <= 5000:
        eigenvalues, eigenvectors = np.linalg.eigh(C)
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        return eigenvectors[:, idx[:L]]
    else:
        # Match R: initial_v <- qnorm(1:(ncol(C)+1) / (ncol(C)+1))[1:ncol(C)]
        v0 = _norm.ppf(np.arange(1, D + 1) / (D + 1))
        v0 = v0 / np.linalg.norm(v0)

        eigenvalues, eigenvectors = eigsh(C, k=L, which='LM', v0=v0,
                                          maxiter=1000)
        idx = np.argsort(np.abs(eigenvalues))[::-1]
        return eigenvectors[:, idx]


# --- Hartigan-Wong (AS 136, 1979) port -------------------------------------
#
# Translation of R's ``src/library/stats/src/kmns.f`` (KMNS + OPTRA + QTRAN).
# The cross-referenced Rust port at ``ddrtree-core/src/init.rs:153-815`` was
# used as a clean reference. This implementation is **NumPy + numba-optional**
# and reproduces R `stats::kmeans(..., algorithm="Hartigan-Wong")` to
# accumulated FP error on the same seed centers.
#
# R-parity contract — DDRTree.R:151 calls ``kmeans(t(Z), K, centers=centers)``.
# Because ``centers=`` is matched by name and ``K`` is positional, R's partial
# matching binds ``K`` to the next open formal (``iter.max``). So the HW outer
# iteration cap is **K, not R's default 10**. We must replicate that.

_HW_BIG = 1.0e30
_HW_QTRAN_STEP_MULTIPLIER = 50  # R's iMaxQtr = 50 * N


@njit(cache=True)
def _hw_initial_assignments(z, centers, ic1, ic2):
    """kmns.f lines 41-75: initial nearest / second-nearest assignment."""
    d, n = z.shape
    k = centers.shape[0]
    for i in range(n):
        ic1[i] = 0
        ic2[i] = 1
        dt1 = 0.0
        dt2 = 0.0
        for j in range(d):
            a = z[j, i] - centers[0, j]; dt1 += a * a
            b = z[j, i] - centers[1, j]; dt2 += b * b
        if dt1 > dt2:  # strict > preserves cluster-0 on ties (kmns.f line 51)
            ic1[i] = 1; ic2[i] = 0
            dt1, dt2 = dt2, dt1
        for l in range(2, k):
            db = 0.0
            completed = True
            for j in range(d):
                delta = z[j, i] - centers[l, j]
                db += delta * delta
                if db >= dt2:  # early abandon (kmns.f line 63)
                    completed = False
                    break
            if not completed:
                continue
            if db >= dt1:
                dt2 = db; ic2[i] = l
            else:
                dt2 = dt1; ic2[i] = ic1[i]
                dt1 = db; ic1[i] = l


@njit(cache=True)
def _hw_transfer(z, centers, point, src, dst, nc, an1, an2):
    """Move point ``point`` from cluster ``src`` to ``dst``. kmns.f 297-311."""
    al1 = float(nc[src]); alw = al1 - 1.0
    al2 = float(nc[dst]); alt = al2 + 1.0
    for j in range(z.shape[0]):
        centers[src, j] = (centers[src, j] * al1 - z[j, point]) / alw
        centers[dst, j] = (centers[dst, j] * al2 + z[j, point]) / alt
    nc[src] -= 1
    nc[dst] += 1
    an2[src] = alw / al1
    an1[src] = (alw / (alw - 1.0)) if alw > 1.0 else _HW_BIG
    an1[dst] = alt / al2
    an2[dst] = alt / (alt + 1.0)


@njit(cache=True)
def _hw_optra(z, centers, ic1, ic2, nc, an1, an2, ncp, dists, itran, live, indx_arr):
    """OPTRA pass — kmns.f 204-329. ``indx_arr`` is a length-1 array used as
    a mutable scalar so numba can write through it."""
    d, n = z.shape
    k = centers.shape[0]
    m_clock = n
    for c in range(k):
        if itran[c] == 1:
            live[c] = m_clock + 1
    for i in range(n):
        i_clock = i + 1
        indx_arr[0] += 1
        l1 = ic1[i]; ll = ic2[i]; l2 = ll
        if nc[l1] != 1:
            if ncp[l1] != 0:
                acc = 0.0
                for j in range(d):
                    delta = z[j, i] - centers[l1, j]
                    acc += delta * delta
                dists[i] = acc * an1[l1]
            r2 = 0.0
            for j in range(d):
                delta = z[j, i] - centers[l2, j]
                r2 += delta * delta
            r2 *= an2[l2]
            for l in range(k):
                if (i_clock >= live[l1] and i_clock >= live[l]) or l == l1 or l == ll:
                    continue
                limit = r2 / an2[l]
                dc = 0.0
                completed = True
                for j in range(d):
                    delta = z[j, i] - centers[l, j]
                    dc += delta * delta
                    if dc >= limit:
                        completed = False
                        break
                if completed:
                    r2 = dc * an2[l]
                    l2 = l
            if r2 >= dists[i]:
                ic2[i] = l2
            else:
                indx_arr[0] = 0
                live[l1] = m_clock + i_clock
                live[l2] = m_clock + i_clock
                ncp[l1] = i_clock
                ncp[l2] = i_clock
                _hw_transfer(z, centers, i, l1, l2, nc, an1, an2)
                ic1[i] = l2; ic2[i] = l1
        if indx_arr[0] == m_clock:
            return
    for c in range(k):
        itran[c] = 0
        live[c] -= m_clock


@njit(cache=True)
def _hw_qtran(z, centers, ic1, ic2, nc, an1, an2, ncp, dists, itran,
              indx_arr, max_qtr):
    """QTRAN pass — kmns.f 332-448. Returns 0 if iMaxQtr exhausted, else 1."""
    d, n = z.shape
    m_clock = n
    icoun = 0
    istep = 0
    while True:
        for i in range(n):
            icoun += 1
            istep += 1
            if istep >= max_qtr:
                return 0
            l1 = ic1[i]; l2 = ic2[i]
            if nc[l1] != 1:
                if istep <= ncp[l1]:
                    acc = 0.0
                    for j in range(d):
                        delta = z[j, i] - centers[l1, j]
                        acc += delta * delta
                    dists[i] = acc * an1[l1]
                if istep < ncp[l1] or istep < ncp[l2]:
                    limit = dists[i] / an2[l2]
                    dc = 0.0
                    completed = True
                    for j in range(d):
                        delta = z[j, i] - centers[l2, j]
                        dc += delta * delta
                        if dc >= limit:
                            completed = False
                            break
                    if completed:
                        icoun = 0
                        indx_arr[0] = 0
                        itran[l1] = 1
                        itran[l2] = 1
                        ncp[l1] = istep + m_clock
                        ncp[l2] = istep + m_clock
                        _hw_transfer(z, centers, i, l1, l2, nc, an1, an2)
                        ic1[i] = l2; ic2[i] = l1
            if icoun == m_clock:
                return 1


def _kmeans_hw(z, seed_centers, iter_max):
    """Hartigan-Wong AS 136. Inputs/outputs match R `stats::kmeans` semantics.

    Parameters
    ----------
    z : (d, n) float64 — same orientation as ``t(Z)`` in R, **transposed
        relative to scikit's API**. Columns are points.
    seed_centers : (k, d) float64 — each row is a seed center, matching R's
        ``centers`` argument (which is row-major (k, d)).
    iter_max : int — OPTRA outer cap. **For DDRTree this must be ``K``**, not
        10, because of R's partial-matching trap at DDRTree.R:151.

    Returns
    -------
    centers : (k, d) float64 — post-iteration cluster means.
    """
    z = np.ascontiguousarray(z, dtype=np.float64)
    centers = np.ascontiguousarray(seed_centers, dtype=np.float64).copy()
    d, n = z.shape
    k = centers.shape[0]
    if k == 0:
        return np.zeros((0, d), dtype=np.float64)
    if n == 0 or k == 1:
        if k == 1 and n > 0:
            return z.mean(axis=1, keepdims=True).T.copy()
        return centers

    ic1 = np.zeros(n, dtype=np.int64)
    ic2 = np.zeros(n, dtype=np.int64)
    _hw_initial_assignments(z, centers, ic1, ic2)

    nc = np.zeros(k, dtype=np.int64)
    for i in range(n):
        nc[ic1[i]] += 1
    any_empty = bool((nc == 0).any())

    # Cluster-mean recompute (kmns.f 77-104). Empty clusters keep seed row.
    means = np.zeros((k, d), dtype=np.float64)
    for i in range(n):
        means[ic1[i]] += z[:, i]
    for c in range(k):
        if nc[c] == 0:
            means[c] = centers[c]
        else:
            means[c] /= nc[c]
    centers = means
    if any_empty:
        return centers

    an1 = np.full(k, _HW_BIG, dtype=np.float64)
    an2 = np.zeros(k, dtype=np.float64)
    ncp = np.full(k, -1, dtype=np.int64)
    itran = np.ones(k, dtype=np.int64)
    live = np.zeros(k, dtype=np.int64)
    dists = np.zeros(n, dtype=np.float64)
    for c in range(k):
        cnt = float(nc[c])
        an2[c] = cnt / (cnt + 1.0)
        if cnt > 1.0:
            an1[c] = cnt / (cnt - 1.0)

    indx_arr = np.zeros(1, dtype=np.int64)  # shared across OPTRA invocations
    max_qtr = _HW_QTRAN_STEP_MULTIPLIER * n
    for _ in range(iter_max):
        _hw_optra(z, centers, ic1, ic2, nc, an1, an2, ncp, dists,
                  itran, live, indx_arr)
        if indx_arr[0] == n:
            break
        qtran_ok = _hw_qtran(z, centers, ic1, ic2, nc, an1, an2, ncp, dists,
                             itran, indx_arr, max_qtr)
        if qtran_ok == 0 or k == 2:
            break
        ncp.fill(0)

    # Final centroid recompute (kmns.f 175-198).
    final = np.zeros((k, d), dtype=np.float64)
    for i in range(n):
        final[ic1[i]] += z[:, i]
    for c in range(k):
        if nc[c] == 0:
            final[c] = centers[c]
        else:
            final[c] /= nc[c]
    return final


def _get_major_eigenvalue(C, L, compat=True):
    """Top squared singular value of ``C`` (D×N).

    DDRTree's objective uses this as ``||X - WZ||_2^2``. We use
    power iteration directly on ``C`` (alternating ``C @ v`` and
    ``C.T @ u``) so we never materialise ``C @ C.T`` or ``C.T @ C``.
    Each iteration costs 2·D·N FMA operations; convergence is
    geometric with rate = (σ₂/σ₁)². For typical DDRTree residuals
    the top singular value is well separated, so 25-40 iterations
    suffice to 1e-5 relative accuracy, giving roughly 30× speed-up
    over forming the D×D Gram matrix and running a full
    ``eigvalsh``.

    Parameters
    ----------
    C : ndarray, (D, N)
    L : int
        Number of dimensions (kept for API compatibility; only the
        top singular value is needed for DDRTree's termination check).
    compat : bool, default True
        When True, clamp the returned σ₁ to ≤ 1.0 to mimic R's IRLBA
        ``max(abs(eigen_res$v))`` bound at DDRTree.R:40 (an upstream bug —
        the IRLBA right-singular vector is unit-norm so its max-abs entry
        is ≤ 1, then squared by DDRTree.cpp:349-352). Set False for the
        mathematically correct ‖X − WZ‖₂ residual. Defaults to True so the
        ``objective_vals`` trajectory matches R's published values.
    """
    C = np.ascontiguousarray(C)
    D, N = C.shape
    if D == 0 or N == 0:
        return 0.0

    rng = np.random.default_rng(0)            # deterministic
    # Iterate in the smaller of the two spaces (D or N) — same math,
    # same answer, but fewer memory reads per matvec.
    if D <= N:
        # matvec: u_new = C @ (C.T @ u)
        u = rng.standard_normal(D)
        u /= np.linalg.norm(u)
        prev = 0.0
        for _ in range(60):
            w = C.T @ u                       # (N,)
            u_new = C @ w                     # (D,)
            s = np.linalg.norm(u_new)
            if s < 1e-30:
                return 0.0
            u_new /= s
            if abs(s - prev) <= 1e-6 * s:
                return float(min(s, 1.0)) if compat else float(s)
            prev, u = s, u_new
        return float(min(s, 1.0)) if compat else float(s)
    else:
        v = rng.standard_normal(N)
        v /= np.linalg.norm(v)
        prev = 0.0
        for _ in range(60):
            w = C @ v                         # (D,)
            v_new = C.T @ w                   # (N,)
            s = np.linalg.norm(v_new)
            if s < 1e-30:
                return 0.0
            v_new /= s
            if abs(s - prev) <= 1e-6 * s:
                return float(min(s, 1.0)) if compat else float(s)
            prev, v = s, v_new
        return float(min(s, 1.0)) if compat else float(s)


def _build_center_mst_nn(Y, N):
    """Final weighted MST on the K centers, embedded in an (N, N) sparse matrix.

    Mirrors ``DDRTree.cpp:528-539``: the returned ``stree`` is
    ``(N_cells, N_cells)`` but only the top-left ``(K, K)`` center block is
    populated. The matrix is symmetric; entries are the squared Euclidean
    distances between MST-linked centers (``_sqdist(Y, Y)``).

    Shared by the pure-Python path and the ``fastddrtree`` backend so both
    return a byte-for-byte identical ``stree`` construction for the same ``Y``.
    """
    distsqMU = _sqdist(Y, Y)
    dist_dense = distsqMU.copy()
    np.fill_diagonal(dist_dense, 0)
    mst_KK = minimum_spanning_tree(csr_matrix(dist_dense))  # (K, K)
    mst_KK_sym = mst_KK + mst_KK.T
    mst_KK_coo = mst_KK_sym.tocoo()
    return _sparse.coo_matrix(
        (mst_KK_coo.data, (mst_KK_coo.row, mst_KK_coo.col)),
        shape=(N, N),
    ).tocsr()


def _ddrtree_fastddrtree(X, *, dimensions, maxIter, sigma, lambda_param,
                         ncenter, param_gamma, tol, verbose, random_state,
                         method='fast', pca_method='irlba', initial_method=None):
    """Dispatch DDRTree to the optional Rust ``fastddrtree`` backend.

    The Rust core is R-parity (HSMM seed=42 procrustes ≈ 4e-14 vs R Monocle 2
    on a C-contiguous input), so its output matches the ``'exact'`` Python path.
    ``method`` and ``pca_method`` select Python-only optimisation strategies and
    are therefore ignored here; a custom ``initial_method`` cannot be honoured
    (the Rust core always initialises ``W`` via PCA) and raises.

    Returns the same 5-key dict contract as the pure-Python ``DDRTree``.
    """
    # Validate unsupported config before attempting the optional import, so the
    # error is actionable regardless of whether fastddrtree is installed.
    if initial_method is not None:
        raise ValueError(
            "backend='fastddrtree' does not support a custom initial_method; "
            "the Rust core always initialises W via PCA. Use backend='python' "
            "for custom initialisation."
        )

    try:
        import fastddrtree
    except ImportError as exc:
        raise ImportError(
            "backend='fastddrtree' requires the optional 'fastddrtree' package "
            "(Rust DDRTree core via pyo3). Install a prebuilt wheel with "
            "`pip install fastddrtree`, or build it locally with "
            "`maturin develop --release` from crates/ddrtree-py."
        ) from exc

    # The pyo3 binding reads X through a zero-copy C-contiguous slice and
    # rejects non-C-contiguous (e.g. Fortran-order) input — silently mis-reading
    # a Fortran-order buffer as C-order scrambles the matrix. Guarantee layout.
    X = np.ascontiguousarray(X, dtype=np.float64)

    # Match the pure-Python DDRTree contract exactly:
    #   * ncenter=None means K=N (every cell is a center). The Rust core would
    #     otherwise auto-select K=cal_ncenter(N) for N>=100 (ddrtree-core
    #     init.rs), so pass N explicitly to keep backend='fastddrtree' a
    #     transparent swap for backend='python'.
    #   * random_state=None is allowed by reduce_dimension's documented API;
    #     map it to the shared default seed (2016) instead of crashing on
    #     int(None). The Python path is deterministic and ignores the seed.
    ncenter_eff = X.shape[1] if ncenter is None else ncenter
    seed_eff = 2016 if random_state is None else int(random_state)

    result = fastddrtree.fit(
        X,
        dimensions=dimensions,
        max_iter=maxIter,
        sigma=sigma,
        ncenter=ncenter_eff,
        lambda_param=lambda_param,
        param_gamma=float(param_gamma),
        tol=tol,
        verbose=verbose,
        random_seed=seed_eff,
    )

    Y = np.ascontiguousarray(result['Y'])
    stree = _build_center_mst_nn(Y, X.shape[1])
    return {
        'W': result['W'],
        'Z': result['Z'],
        'Y': Y,
        'stree': stree,
        'objective_vals': list(result['objective_vals']),
    }


def DDRTree(X, dimensions=2, initial_method=None, maxIter=20, sigma=0.001,
            lambda_param=None, ncenter=None, param_gamma=10, tol=0.001,
            verbose=False, pca_method='irlba', random_state=2016,
            method='fast', backend='python'):
    """
    Perform DDRTree dimensionality reduction.

    Parameters
    ----------
    method : {'fast', 'exact'}, default 'fast'
        ``'fast'`` (default) runs a reformulated update that caches
        ``X @ X.T``, truncates the soft-assignment matrix ``R`` to its
        top-``K/5`` entries per row (safe for the default
        ``sigma=0.001``), and uses
        ``||Y_new − Y_old||_F / ||Y_old||_F < tol`` as the termination
        criterion.  Trajectory topology and pseudotime *correlation*
        with the exact result are preserved (typically 0.99+), but
        absolute pseudotime values may differ slightly.  About 3× faster
        per call on typical datasets.

        ``'exact'`` reproduces R Monocle 2's convergence exactly by
        evaluating the full objective (including ``||X - WZ||_2^2``)
        on every iteration.  Use this when bitwise agreement with R
        is required.
    backend : {'python', 'fastddrtree'}, default 'python'
        ``'python'`` runs the pure-Python implementation in this module.
        ``'fastddrtree'`` dispatches to the optional Rust core (pyo3 binding,
        ``pip install fastddrtree``), which is R-parity and typically faster.
        With ``'fastddrtree'`` the ``method``/``pca_method`` options are ignored
        (the Rust core uses its own R-faithful update path) and a custom
        ``initial_method`` raises ``ValueError``.
    X : np.ndarray, shape (D, N)
        Input data matrix (genes x cells), already preprocessed.
    dimensions : int
        Number of dimensions to reduce to.
    initial_method : callable or None
        Custom initialization method. If None, uses PCA.
    maxIter : int
        Maximum number of iterations.
    sigma : float
        Bandwidth parameter for soft assignment.
    lambda_param : float or None
        Regularization parameter. If None, set to 5*N.
    ncenter : int or None
        Number of centers. If None, K=N (all cells are centers).
    param_gamma : float
        Gamma parameter controlling the trade-off.
    tol : float
        Convergence tolerance.
    verbose : bool
        Print progress information.

    Returns
    -------
    dict with keys:
        W : (D, dimensions) projection matrix
        Z : (dimensions, N) reduced coordinates of cells
        Y : (dimensions, K) center positions in reduced space
        stree : (N, N) sparse MST adjacency (weighted, symmetric). Entries are
            populated only in the top-left (K, K) block — matches
            ``DDRTree.cpp:538`` which allocates ``SpMat(N_cells, N_cells)``.
        objective_vals : list of objective function values
    """
    # Dtype + finiteness validation. R/C++ DDRTree runs everything in double;
    # silent float32 promotion in BLAS paths can shift Cholesky and eigsh
    # results by O(1e-5), which compounds across iterations. NaN/Inf in X
    # would cascade through min, exp, Cholesky, and MST with no diagnostic.
    X = np.asarray(X)
    if X.dtype != np.float64:
        X = np.ascontiguousarray(X, dtype=np.float64)
    else:
        X = np.ascontiguousarray(X)
    if not np.all(np.isfinite(X)):
        raise ValueError(
            "X contains non-finite values (NaN/Inf); DDRTree requires a "
            "finite float64 matrix."
        )
    D, N = X.shape

    if backend == 'fastddrtree':
        return _ddrtree_fastddrtree(
            X, dimensions=dimensions, maxIter=maxIter, sigma=sigma,
            lambda_param=lambda_param, ncenter=ncenter, param_gamma=param_gamma,
            tol=tol, verbose=verbose, random_state=random_state,
            method=method, pca_method=pca_method, initial_method=initial_method,
        )
    if backend != 'python':
        raise ValueError(
            f"backend must be 'python' or 'fastddrtree', got {backend!r}"
        )

    # PCA initialization for W — cache XXT so the fast-mode inner loop
    # can reuse it instead of recomputing ``Q @ X.T`` (D²·N) every
    # iteration.  The initial ``X @ X.T`` was already needed anyway.
    XXT = X @ X.T
    W = _pca_projection(XXT, dimensions)

    # Initialize Z
    if initial_method is None:
        Z = W.T @ X
    else:
        tmp = initial_method(X)
        Z = tmp[:, :dimensions].T

    # Initialize Y (centers)
    if ncenter is None:
        K = N
        Y = Z.copy()
    else:
        K = ncenter
        if K > Z.shape[1]:
            raise ValueError("ncenter must be <= number of cells")
        # R: centers = t(Z)[seq(1, ncol(Z), length.out=K), ]  — evenly spaced
        # indices, 1-based. np.linspace with dtype=int truncates (floor for
        # nonneg), matching `as.integer(seq(...))` semantics.
        indices = np.linspace(0, Z.shape[1] - 1, K, dtype=int)
        seed_centers = Z[:, indices].T              # (K, dim)
        # R kmeans(t(Z), K, centers=centers) — `K` positional binds iter.max
        # via partial matching. So HW outer cap = K, not 10. See _kmeans_hw.
        Y_kd = _kmeans_hw(Z, seed_centers, iter_max=K)
        Y = Y_kd.T                                  # (dim, K)

    if lambda_param is None:
        lambda_param = 5.0 * N

    # Main iterative optimization
    objective_vals = []
    B = np.zeros((K, K))
    _prev_Y = None   # for the 'fast' method's cheap convergence check

    for iteration in range(maxIter):
        if verbose:
            print(f"Iteration: {iteration}")

        # Step 1: Compute MST on Y
        distsqMU = _sqdist(Y, Y)

        # Build full distance matrix for MST computation
        # scipy MST works on dense or sparse, expects (K, K) with distances
        # Use scipy's minimum_spanning_tree which expects a CSR matrix
        dist_dense = distsqMU.copy()
        np.fill_diagonal(dist_dense, 0)
        mst_sparse = minimum_spanning_tree(csr_matrix(dist_dense))
        mst_array = mst_sparse.toarray()

        # Build symmetric adjacency B
        B = np.zeros((K, K))
        B[mst_array > 0] = 1
        B = B + B.T
        B[B > 0] = 1

        # Graph Laplacian
        L = np.diag(B.sum(axis=1)) - B

        # Step 2: Soft assignment of cells to centers
        distZY = _sqdist(Z, Y)  # (N, K)

        min_dist = distZY.min(axis=1, keepdims=True)  # (N, 1)
        tmp_distZY = distZY - min_dist  # (N, K)

        tmp_R = np.exp(-tmp_distZY / sigma)  # (N, K)
        R = tmp_R / tmp_R.sum(axis=1, keepdims=True)  # (N, K)

        # Gamma = diag(colSums(R))
        Gamma = np.diag(R.sum(axis=0))  # (K, K)

        # Termination check
        if method == 'fast':
            # Skip the expensive ``||X - WZ||_2^2`` evaluation and use
            # the change in Y — the learned tree centres — as a
            # convergence signal.  Frobenius norm on Y (K×dim) is
            # effectively free.
            #
            # Empirically the Y-change signal plateaus around ≈ tol_Y
            # because of MST topology jitter (edges flipping between
            # near-equivalent candidates), so we relax the threshold
            # vs. the 'exact' objective-change criterion. The default
            # scaling of 20× matches 'exact' convergence on pancreas /
            # HSMM in ~5–8 iterations.
            if iteration >= 1 and _prev_Y is not None:
                dy = float(np.linalg.norm(Y - _prev_Y))
                norm_y = max(float(np.linalg.norm(_prev_Y)), 1e-12)
                delta_y = dy / norm_y
                if verbose:
                    print(f"  Iter {iteration}: delta_Y = {delta_y:.6e}")
                if delta_y < 20.0 * tol:
                    if verbose:
                        print("Converged (fast)!")
                    objective_vals.append(delta_y)
                    break
            objective_vals.append(0.0)     # placeholder, not used downstream
            _prev_Y = Y.copy()
        else:
            x1 = np.log(np.exp(-tmp_distZY / sigma).sum(axis=1))
            obj1 = -sigma * (x1 - min_dist.flatten() / sigma).sum()
            try:
                major_ev = _get_major_eigenvalue(X - W @ Z, dimensions, compat=True)
            except Exception:
                major_ev = np.linalg.norm(X - W @ Z, ord=2)
            obj2 = major_ev ** 2
            obj2 += lambda_param * np.trace(Y @ L @ Y.T) + param_gamma * obj1
            objective_vals.append(obj2)

            if verbose:
                print(f"  Objective: {obj2:.6f}")

            if iteration >= 1:
                delta_obj = abs(objective_vals[-1] - objective_vals[-2])
                delta_obj /= abs(objective_vals[-2]) if objective_vals[-2] != 0 else 1.0
                if verbose:
                    print(f"  delta_obj: {delta_obj:.8f}")
                if delta_obj < tol:
                    if verbose:
                        print("Converged!")
                    break

        # Step 3: Update W, Z, Y
        if method == 'fast':
            # Reformulated update that avoids every ``O(N·K·D)`` dense
            # intermediate the 'exact' path materialises.  Three ideas:
            #
            #   (A) Truncate R to its top-``K_trunc`` entries per row.
            #       For the default ``sigma=0.001`` the soft-assignment
            #       kernel ``exp(-d²/σ)`` is numerically zero beyond a
            #       handful of nearest centres, so the truncation error
            #       is ~1e-5.  We renormalise rows so they still sum to 1.
            #   (B) Cache ``XXT = X @ X.T`` (already needed by the
            #       initial PCA) and use the identity
            #           Q @ X.T = (XXT + Z_mat · XR.T) / (γ+1)
            #       where ``XR = X @ R``, ``Z_mat = XR · A_mat⁻¹``.
            #       This replaces one ``D²·N`` matmul and one ``D·N·K``
            #       matmul with two ``D·K·K`` matmuls.
            #   (C) Solve ``A_mat · Z_mat.T = XR.T`` (``K × D`` RHS)
            #       instead of ``A_mat · tmp.T = R.T`` (``K × N`` RHS) —
            #       since ``K << N`` this alone saves ``N/K × K³`` ops.
            #
            # Net: per-iteration cost drops from ~1.7G to ~0.1G ops on
            # pancreas (3696 cells, K=308, D=300).
            K_trunc = min(K, max(30, K // 5))
            if K_trunc < K:
                # Keep top-K_trunc entries per row, renormalise.
                top_idx = np.argpartition(R, -K_trunc, axis=1)[:, -K_trunc:]
                rows = np.repeat(np.arange(N), K_trunc)
                cols = top_idx.ravel()
                vals = R[rows, cols]
                R_sparse = _sparse.csr_matrix((vals, (rows, cols)),
                                              shape=(N, K))
                row_sums = np.asarray(R_sparse.sum(axis=1)).ravel()
                row_sums[row_sums == 0] = 1.0
                R_sparse = _sparse.diags(1.0 / row_sums) @ R_sparse
                Gamma_diag = np.asarray(R_sparse.sum(axis=0)).ravel()
                Gamma = np.diag(Gamma_diag)
                RtR = (R_sparse.T @ R_sparse).toarray()
            else:
                R_sparse = R       # dense fall-through for very small K
                RtR = R.T @ R

            A_mat = ((param_gamma + 1.0) / param_gamma) * (
                (lambda_param / param_gamma) * L + Gamma
            ) - RtR

            # XR = X @ R  (D × K, exploits sparsity of R)
            if _sparse.issparse(R_sparse):
                XR = np.asarray(X @ R_sparse)
            else:
                XR = X @ R_sparse

            try:
                cho = cho_factor(A_mat)
                Z_mat = cho_solve(cho, XR.T).T        # (D, K)
            except np.linalg.LinAlgError:
                Z_mat = np.linalg.solve(A_mat, XR.T).T

            # sym_mat = (XXT + 0.5 (Z_mat·XR.T + XR·Z_mat.T)) / (γ+1)
            sym_inner = Z_mat @ XR.T                  # (D, D)
            sym_mat = (XXT + 0.5 * (sym_inner + sym_inner.T)) \
                / (param_gamma + 1.0)
            W = _pca_projection_irlba_like(sym_mat, dimensions)

            # Z = W.T @ Q  without materialising Q.
            Wx = W.T @ X                              # (dim, N)
            WZmat = W.T @ Z_mat                       # (dim, K)
            if _sparse.issparse(R_sparse):
                WZmat_Rt = np.asarray(WZmat @ R_sparse.T)
            else:
                WZmat_Rt = WZmat @ R_sparse.T
            Z = (Wx + WZmat_Rt) / (param_gamma + 1.0)

            # Y update: same form as exact, but uses sparse R.
            A_Y = (lambda_param / param_gamma) * L + Gamma
            if _sparse.issparse(R_sparse):
                ZR = np.asarray(Z @ R_sparse)
            else:
                ZR = Z @ R_sparse
            try:
                cho_Y = cho_factor(A_Y)
                Y = cho_solve(cho_Y, ZR.T).T
            except np.linalg.LinAlgError:
                Y = np.linalg.solve(A_Y, ZR.T).T
            continue  # skip the exact-path update block below

        # tmp = solve(((gamma+1)/gamma) * (lambda/gamma * L + Gamma) - R^T R, R^T)
        A_mat = ((param_gamma + 1.0) / param_gamma) * (
            (lambda_param / param_gamma) * L + Gamma
        ) - R.T @ R

        try:
            # Use Cholesky if positive definite
            cho = cho_factor(A_mat)
            tmp_dense = cho_solve(cho, R.T).T  # (N, K)
        except np.linalg.LinAlgError:
            tmp_dense = np.linalg.solve(A_mat, R.T).T  # (N, K)

        # Q = (X + (X @ tmp_dense) @ R^T) / (gamma + 1)
        Q = (X + (X @ tmp_dense) @ R.T) / (param_gamma + 1.0)  # (D, N)

        # C = Q (in the optimized path)
        C = Q

        # W = pca_projection((Q X^T + X Q^T)/2, dimensions)
        # sym_mat = (Q X^T + X Q^T) / 2 is D×D but rank ≤ N
        # Use low-rank trick: if M = Q X^T, then M + M^T has same eigenvectors
        # as the N×N matrix (X^T Q + Q^T X) / 2, mapped back via X or Q.
        # Specifically: top eigvecs of M M^T can be found via N×N problem.
        #
        # Direct approach: form the N×N Gram matrices and solve exactly.
        # Let A = Q (D×N), B = X (D×N). M = A B^T.
        # (M+M^T)/2 v = λv  ⟺  (AB^T + BA^T)/2 v = λv
        # Multiply by B^T: B^T(AB^T + BA^T)v/2 = λ B^T v
        # Let u = B^T v (N×1): (B^T A)(B^T v) + (B^T B)(A^T v) = 2λ u  — coupled.
        #
        # Simpler: since rank(sym_mat) ≤ 2N, form [Q, X] (D×2N) and solve
        # the 2N×2N eigenvalue problem.
        QX = np.hstack([Q, X])  # (D, 2N)
        small = QX.T @ QX       # (2N, 2N)
        # sym_mat = (Q X^T + X Q^T)/2, so sym_mat @ QX = QX @ M_small
        # where M_small = [[Q^T Q, Q^T X], [X^T Q, X^T X]] ... not quite.
        # Actually: sym_mat = (Q X^T + X Q^T)/2
        # sym_mat @ v = λv means we need eigvecs of sym_mat.
        # Since sym_mat has rank ≤ N, its top eigvecs live in span(Q, X).
        # Write v = Q a + X b = [Q X] [a; b]
        # Then [Q X]^T sym_mat [Q X] [a;b] = λ [Q X]^T [Q X] [a;b]
        # This is a generalized eigenvalue problem of size 2N×2N.
        #
        # But even simpler: Q = (X + X @ tmp_dense @ R^T)/(γ+1)
        # So Q is in span(X), meaning sym_mat = (QX^T + XQ^T)/2 is rank ≤ N.
        # v = X c for some c. Then X^T sym_mat X c = λ X^T X c
        # X^T (QX^T + XQ^T)/2 X c = λ (X^T X) c
        # ((X^T Q)(X^T X) + (X^T X)(Q^T X))/2 c = λ (X^T X) c

        # sym_mat = (Q X^T + X Q^T) / 2 — top eigvecs give new W.
        #
        # Factoring: sym_mat is the symmetric part of Q X^T, so we
        # only need ONE D×D matmul (not two). Halves the cost of this
        # step, which was the largest numpy hotspot in the loop.
        if pca_method == 'irlba':
            M_qx = Q @ X.T                           # D×D, single matmul
            sym_mat = 0.5 * (M_qx + M_qx.T)
            W = _pca_projection_irlba_like(sym_mat, dimensions)
        elif D <= N:
            M_qx = Q @ X.T
            sym_mat = 0.5 * (M_qx + M_qx.T)
            evals_all, evecs_all = np.linalg.eigh(sym_mat)
            idx = np.argsort(evals_all)[::-1][:dimensions]
            W = evecs_all[:, idx]
        else:
            # Low-rank trick for D > N
            XtX = X.T @ X
            XtQ = X.T @ Q
            lhs = (XtQ @ XtX + XtX @ XtQ.T) / 2.0
            from scipy.linalg import eigh as eigh_gen
            try:
                reg = 1e-10 * np.trace(XtX) / N
                evals, evecs = eigh_gen(
                    lhs, XtX + reg * np.eye(N),
                    subset_by_index=[N - dimensions, N - 1])
                idx = np.argsort(evals)[::-1]
                evecs = evecs[:, idx]
                W = X @ evecs
                for col in range(dimensions):
                    norm_val = np.linalg.norm(W[:, col])
                    if norm_val > 0:
                        W[:, col] /= norm_val
            except np.linalg.LinAlgError:
                sym_mat = (Q @ X.T + X @ Q.T) / 2.0
                W = _pca_projection(sym_mat, dimensions)

        # Z = W^T @ C
        Z = W.T @ C

        # Y = solve(lambda/gamma * L + Gamma, (Z @ R)^T)^T
        A_Y = (lambda_param / param_gamma) * L + Gamma
        try:
            cho = cho_factor(A_Y)
            Y = cho_solve(cho, (Z @ R).T).T
        except np.linalg.LinAlgError:
            Y = np.linalg.solve(A_Y, (Z @ R).T).T

    # Build final MST with weights. Per DDRTree.cpp:528-539, the returned
    # ``stree`` is a sparse (N, N) — NOT (K, K) — matrix populated only at
    # center-pair positions. Downstream Monocle 2 consumers (e.g. the
    # ``project2MST`` helper) index ``stree`` with cell indices and expect
    # the larger shape.
    stree = _build_center_mst_nn(Y, N)

    return {
        'W': W,
        'Z': Z,
        'Y': Y,
        'stree': stree,
        'objective_vals': objective_vals,
    }
