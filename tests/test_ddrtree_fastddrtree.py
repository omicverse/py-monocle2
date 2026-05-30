"""Tests for the optional Rust ``fastddrtree`` backend of ``DDRTree``.

The backend is dispatched via ``DDRTree(..., backend='fastddrtree')`` and must
return the exact same 5-key dict contract as the pure-Python path.

Tests that actually invoke the Rust core are skipped when the optional
``fastddrtree`` package is not installed (e.g. CI without the wheel). The
rigorous R-parity claim (HSMM seed=42 procrustes ≈ 4e-14 vs R Monocle 2) is
validated separately by the fastDDRTree validation harness on a C-contiguous
input; these tests only assert the Python-side integration contract.
"""
from __future__ import annotations

import importlib.util

import numpy as np
import pytest
from scipy.spatial import procrustes

from monocle2_py.ddrtree import DDRTree, _build_center_mst_nn

HAVE_FASTDDRTREE = importlib.util.find_spec("fastddrtree") is not None
requires_rust = pytest.mark.skipif(
    not HAVE_FASTDDRTREE, reason="optional 'fastddrtree' package not installed"
)


def _synthetic_two_clusters(D=40, n=80, seed=0):
    """(D, N) matrix with two well-separated Gaussian blobs."""
    rng = np.random.default_rng(seed)
    left = rng.normal(0.0, 1.0, (D, n // 2))
    right = rng.normal(6.0, 1.0, (D, n - n // 2))
    return np.hstack([left, right])


# --- backend dispatch / validation (no Rust core required) -------------------

def test_unknown_backend_raises():
    X = _synthetic_two_clusters()
    with pytest.raises(ValueError, match="backend must be 'python' or 'fastddrtree'"):
        DDRTree(X, dimensions=2, ncenter=10, backend="nope")


def test_initial_method_rejected_for_fastddrtree():
    """Custom init is unsupported by the Rust core and must fail fast — even
    when the optional package is absent (the check precedes the import)."""
    X = _synthetic_two_clusters()
    with pytest.raises(ValueError, match="does not support a custom initial_method"):
        DDRTree(X, dimensions=2, ncenter=10, backend="fastddrtree",
                initial_method=lambda x: x.T)


def test_default_backend_is_python():
    """The default path must remain pure-Python (no hard fastddrtree dep)."""
    X = _synthetic_two_clusters()
    res = DDRTree(X, dimensions=2, ncenter=10, maxIter=5)
    assert set(res) == {"W", "Z", "Y", "stree", "objective_vals"}


# --- shared stree helper (behaviour-preserving refactor) ---------------------

def test_build_center_mst_nn_shape_and_symmetry():
    rng = np.random.default_rng(3)
    K, N = 12, 50
    Y = rng.normal(size=(2, K))
    stree = _build_center_mst_nn(Y, N)
    assert stree.shape == (N, N)
    # A spanning tree over K nodes has K-1 edges; symmetric storage => 2(K-1).
    assert stree.nnz == 2 * (K - 1)
    dense = stree.toarray()
    np.testing.assert_allclose(dense, dense.T, rtol=0, atol=0)
    # Only the top-left (K, K) block is populated.
    assert dense[K:, :].sum() == 0
    assert dense[:, K:].sum() == 0


# --- Rust backend: contract (requires the optional package) ------------------

@requires_rust
def test_fastddrtree_returns_contract_shapes():
    X = _synthetic_two_clusters(D=40, n=80)
    D, N = X.shape
    K = 12
    res = DDRTree(X, dimensions=2, ncenter=K, maxIter=10, backend="fastddrtree")

    assert set(res) == {"W", "Z", "Y", "stree", "objective_vals"}
    assert res["W"].shape == (D, 2)
    assert res["Z"].shape == (2, N)
    assert res["Y"].shape == (2, K)
    assert res["stree"].shape == (N, N)
    assert res["stree"].nnz == 2 * (K - 1)
    assert len(res["objective_vals"]) >= 1
    for key in ("W", "Z", "Y"):
        assert np.all(np.isfinite(res[key]))
        assert res[key].dtype == np.float64


@requires_rust
def test_fastddrtree_is_deterministic():
    X = _synthetic_two_clusters(seed=1)
    a = DDRTree(X, dimensions=2, ncenter=10, random_state=7, backend="fastddrtree")
    b = DDRTree(X, dimensions=2, ncenter=10, random_state=7, backend="fastddrtree")
    np.testing.assert_array_equal(a["W"], b["W"])
    np.testing.assert_array_equal(a["Z"], b["Z"])
    np.testing.assert_array_equal(a["Y"], b["Y"])


@requires_rust
def test_fastddrtree_accepts_fortran_order_input():
    """A Fortran-order X must be coerced to C-contiguous before reaching the
    binding (which otherwise rejects / mis-reads non-C-contiguous buffers)."""
    X = np.asfortranarray(_synthetic_two_clusters(seed=2))
    assert not X.flags["C_CONTIGUOUS"]
    res = DDRTree(X, dimensions=2, ncenter=10, maxIter=5, backend="fastddrtree")
    assert res["W"].shape == (X.shape[0], 2)
    assert np.all(np.isfinite(res["Z"]))


@requires_rust
def test_fastddrtree_structurally_agrees_with_python_exact():
    """Sanity check (NOT a parity proof): on a well-separated synthetic the
    Rust Y and the Python 'exact' Y should align closely after Procrustes.
    The hard R-parity number lives in the validation harness."""
    X = _synthetic_two_clusters(D=50, n=100, seed=5)
    K = 16
    rust = DDRTree(X, dimensions=2, ncenter=K, maxIter=20, backend="fastddrtree")
    py = DDRTree(X, dimensions=2, ncenter=K, maxIter=20, method="exact",
                 backend="python")
    # procrustes() returns a disparity in [0, 1]; standardised, alignment-free.
    _, _, disparity = procrustes(rust["Y"].T, py["Y"].T)
    assert disparity < 0.1, f"Y disparity {disparity:.4f} too large"


@requires_rust
def test_fastddrtree_ncenter_none_honours_K_equals_N():
    """ncenter=None must yield K=N to match the pure-Python contract, even for
    N>=100 where the Rust core would otherwise auto-select fewer centers."""
    X = _synthetic_two_clusters(D=30, n=120)  # N=120 >= 100
    N = X.shape[1]
    res = DDRTree(X, dimensions=2, ncenter=None, maxIter=5, backend="fastddrtree")
    assert res["Y"].shape == (2, N)
    assert res["Z"].shape == (2, N)
    assert res["stree"].shape == (N, N)
    # Same default on the Python path — backends must agree on K.
    py = DDRTree(X, dimensions=2, ncenter=None, maxIter=5, backend="python")
    assert res["Y"].shape[1] == py["Y"].shape[1] == N


@requires_rust
def test_fastddrtree_accepts_random_state_none():
    """random_state=None is part of reduce_dimension's documented API and is
    forwarded into DDRTree; the fast backend must map it to the default seed
    rather than crash on int(None)."""
    X = _synthetic_two_clusters(seed=4)
    res = DDRTree(X, dimensions=2, ncenter=10, maxIter=5, random_state=None,
                  backend="fastddrtree")
    assert res["Y"].shape == (2, 10)
    assert np.all(np.isfinite(res["W"]))
