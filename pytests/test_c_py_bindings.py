"""Tests for C ABI Python binding vs Fortran-consistent C rule entry."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

# Prefer interfaces/PyInterface on path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "interfaces" / "PyInterface"))

# Library must be available
if "GJP_CINTERP" not in os.environ:
    for cand in (ROOT / "bbdir").rglob("libgjp_cinterp.so"):
        os.environ["GJP_CINTERP"] = str(cand)
        break

pytest.importorskip("gauss_jacobi_quad")
from gauss_jacobi_quad import (  # noqa: E402
    GJP_ERR_BOGAERT_AB,
    GaussJacobiError,
    gauss_jacobi,
)


def test_auto_legendre_matches_bogaert():
    x_a, w_a = gauss_jacobi(32, 0.0, 0.0, method=None)
    x_b, w_b = gauss_jacobi(32, 0.0, 0.0, method="bogaert")
    assert x_a.shape == (32,)
    assert np.allclose(x_a, x_b, atol=0.0)
    assert np.allclose(w_a, w_b, atol=0.0)
    assert np.all(w_a > 0)
    assert np.all(np.diff(x_a) > 0)


def test_glr_matches_algo665():
    x_g, w_g = gauss_jacobi(40, 0.5, 0.5, method="glr")
    x_r, w_r = gauss_jacobi(40, 0.5, 0.5, method="algo665")
    assert np.allclose(x_g, x_r, atol=1e-12)
    assert np.allclose(w_g, w_r, atol=1e-11)


def test_sturm_matches_algo665():
    x_s, w_s = gauss_jacobi(40, 0.5, 0.5, method="sturm")
    x_r, w_r = gauss_jacobi(40, 0.5, 0.5, method="algo665")
    assert np.allclose(x_s, x_r, atol=1e-12)
    assert np.allclose(w_s, w_r, atol=1e-11)


def test_bogaert_non_legendre_raises():
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi(32, 0.5, 0.0, method="bogaert")
    assert ei.value.status == GJP_ERR_BOGAERT_AB
