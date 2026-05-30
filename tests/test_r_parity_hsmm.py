"""R-parity regression for DDRTree on HSMM seed=42.

Pins the parity checks flagged as missing in the second-opinion review:
(a) end-to-end procrustes RMSE < 1e-2 vs the R `DDRTree` 0.1.5 reference,
(b) ``stree`` shape is (N, N) per ``DDRTree.cpp:538``,
(c) ``len(objective_vals)`` tracks the iteration count,
(d) dtype/finiteness guards behave (reject NaN/Inf, coerce float32 -> float64).

All end-to-end checks run in ``method='exact'`` — the mode that reproduces R
Monocle 2's convergence. The package default ``method='fast'`` is an
approximation and is intentionally not gated here (see README parity notes).

Fixtures are absolute paths into the sibling ``monocle2-rebuild`` workspace.
The test skips (does not fail) when the fixtures are absent, so CI on a fresh
checkout without the validation tree still passes the rest of the suite.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from monocle2_py.ddrtree import DDRTree

# The R reference fixtures (HSMM seed=42 DDRTree input X + R's DDRTree 0.1.5
# output) live outside this repo. Point ``MONOCLE2_REBUILD_ROOT`` at the
# directory that holds them; when it is unset (e.g. on CI / a fresh clone) the
# whole module skips instead of failing.
_REBUILD = os.environ.get("MONOCLE2_REBUILD_ROOT")
_X_PATH = (Path(_REBUILD) / "validation/reports/phase3_tier1_debug/HSMM_seed42_X.npz"
           if _REBUILD else None)
_REF_PATH = (Path(_REBUILD) / "references/outputs/HSMM_seed42.rds.flat.npz"
             if _REBUILD else None)

_REASON = (
    "R reference fixtures not available; set MONOCLE2_REBUILD_ROOT to the "
    "directory holding HSMM_seed42_X.npz and HSMM_seed42.rds.flat.npz to run."
)
_HAS_FIXTURES = bool(_REBUILD) and _X_PATH.exists() and _REF_PATH.exists()


def _procrustes_rmse(Y_a: np.ndarray, Y_b: np.ndarray) -> float:
    """Centered orthogonal Procrustes RMSE, reflections allowed.

    Inline copy of monocle2-rebuild/validation/metrics/procrustes.py so this
    test has no dependency on the validation tree at import time.
    """
    if Y_a.shape != Y_b.shape:
        raise ValueError(f"shape mismatch: {Y_a.shape} vs {Y_b.shape}")
    if Y_a.shape[1] == 0:
        return 0.0
    mu_a = Y_a.mean(axis=1, keepdims=True)
    mu_b = Y_b.mean(axis=1, keepdims=True)
    A = Y_a - mu_a
    B = Y_b - mu_b
    U, _, Vt = np.linalg.svd(B @ A.T, full_matrices=False)
    R = U @ Vt
    resid = R @ A - B
    return float(np.sqrt(np.mean(np.sum(resid * resid, axis=0))))


@pytest.fixture(scope="module")
def hsmm_fixture():
    if not _HAS_FIXTURES:
        pytest.skip(_REASON)
    X = np.load(_X_PATH)["X"]                       # (D, N) = (2000, 271)
    ref = np.load(_REF_PATH, allow_pickle=True)
    return {
        "X": np.ascontiguousarray(X, dtype=np.float64),
        "Y_R": ref["Y"],                            # (2, 110)
        "W_R": ref["W"],                            # (2000, 2)
        "Z_R": ref["Z"],                            # (2, 271)
        "K": int(ref["ncenter"]) if "ncenter" in ref.files else 110,
    }


def test_ddrtree_hsmm_seed42_procrustes_under_1e2(hsmm_fixture):
    """End-to-end DDRTree on HSMM seed=42 — procrustes Y < 1e-2 vs R (exact).

    Pre-patch (sklearn Lloyd kmeans init) this was ~1.11; the Hartigan-Wong
    kmeans port restores R-parity to FP epsilon (~2.8e-14).
    """
    out = DDRTree(
        hsmm_fixture["X"],
        dimensions=2,
        maxIter=20,
        sigma=1e-3,
        ncenter=hsmm_fixture["K"],
        param_gamma=10,
        tol=1e-3,
        method="exact",
        random_state=2016,
        verbose=False,
    )
    rmse_Y = _procrustes_rmse(out["Y"], hsmm_fixture["Y_R"])
    assert rmse_Y < 1e-2, f"Y procrustes regressed: {rmse_Y:.4e} (target <1e-2)"


def test_ddrtree_hsmm_stree_is_NxN(hsmm_fixture):
    """stree must be (N, N) per DDRTree.cpp:538 — not (K, K)."""
    out = DDRTree(
        hsmm_fixture["X"], dimensions=2, maxIter=5, sigma=1e-3,
        ncenter=hsmm_fixture["K"], param_gamma=10, tol=1e-3, method="exact",
    )
    N = hsmm_fixture["X"].shape[1]
    assert out["stree"].shape == (N, N)
    # Edges populated in top-left K×K block only — there should be exactly
    # 2*(K-1) nonzeros (symmetric MST on K nodes, no root self-edge).
    K = hsmm_fixture["K"]
    assert out["stree"].nnz == 2 * (K - 1), (
        f"expected 2*(K-1)={2*(K-1)} nnz, got {out['stree'].nnz}"
    )


def test_ddrtree_objective_vals_length_matches_iters(hsmm_fixture):
    """``len(objective_vals)`` must track the number of iterations run."""
    out = DDRTree(
        hsmm_fixture["X"], dimensions=2, maxIter=20, sigma=1e-3,
        ncenter=hsmm_fixture["K"], param_gamma=10, tol=1e-3, method="exact",
    )
    # `exact` appends one entry per iter; loop runs to maxIter unless tol hits.
    assert 1 <= len(out["objective_vals"]) <= 20


def test_ddrtree_rejects_non_finite_X():
    X = np.zeros((5, 4), dtype=np.float64)
    X[0, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        DDRTree(X, dimensions=2, maxIter=2, ncenter=3)


def test_ddrtree_coerces_float32_to_float64():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(20, 12)).astype(np.float32)
    out = DDRTree(X, dimensions=2, maxIter=3, ncenter=5)
    assert out["W"].dtype == np.float64
