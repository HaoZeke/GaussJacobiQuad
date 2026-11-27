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
import gauss_jacobi_quad as _gjq  # noqa: E402

GJP_ERR_BOGAERT_AB = _gjq.GJP_ERR_BOGAERT_AB
GaussJacobiError = _gjq.GaussJacobiError
gauss_jacobi = _gjq.gauss_jacobi
gauss_jacobi_scale = _gjq.gauss_jacobi_scale
GJP_ERR_BOUNDS = _gjq.GJP_ERR_BOUNDS


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


def test_scale_polynomial_on_interval():
    """∫_0^3 (3-t)^1 t^2 t^3 dt = 729/14."""
    alpha, beta, a, b = 1.0, 2.0, 0.0, 3.0
    x, w = gauss_jacobi(8, alpha, beta, method="rec")
    xs, ws = gauss_jacobi_scale(x, w, alpha, beta, a, b)
    quad = float(np.dot(ws, xs**3))
    assert abs(quad - 729.0 / 14.0) < 1e-12
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, w, alpha, beta, b, a)
    assert ei.value.status == GJP_ERR_BOUNDS
