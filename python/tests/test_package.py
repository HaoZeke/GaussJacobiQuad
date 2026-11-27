# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# -----------------------------------------------------------------------------
# END_HEADER
"""Tests against the CPython extension (not ctypes)."""
from __future__ import annotations

import importlib.util
import math
import sys

import numpy as np
import pytest
from gauss_jacobi_quad import (
    GJP_ERR_ALPHA,
    GJP_ERR_BETA,
    GJP_ERR_BOGAERT_AB,
    GJP_ERR_BOUNDS,
    GJP_OK,
    GaussJacobiError,
    __version__,
    gauss_jacobi,
    gauss_jacobi_scale,
    status_string,
)


def test_is_compiled_extension():
    import gauss_jacobi_quad._core as core

    assert hasattr(core, "rule")
    assert hasattr(core, "GaussJacobiError")
    assert core.__name__ == "gauss_jacobi_quad._core"
    # Must be a real extension module, not a pure-Python / ctypes loader
    spec = importlib.util.find_spec("gauss_jacobi_quad._core")
    assert spec is not None
    assert spec.origin is not None
    assert not spec.origin.endswith(".py")
    assert "ctypes" not in (spec.origin or "").lower()


def test_version():
    assert __version__ == "0.2.6"


def test_extension_so_naming():
    """Wheel path is either Stable ABI (abi3) or free-threaded (cpXXXt)."""
    import importlib.metadata as md
    import sysconfig

    try:
        dist = md.distribution("gauss-jacobi-quad")
    except md.PackageNotFoundError:
        pytest.skip("not installed as a distribution")
    files = [str(f) for f in (dist.files or [])]
    so = [
        f for f in files if "_core" in f and (f.endswith(".so") or f.endswith(".pyd"))
    ]
    if not so:
        pytest.skip("no extension file recorded")
    path = so[0]
    ft = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))
    if ft:
        # Free-threaded wheels use the full ABI tag, not abi3
        assert "abi3" not in path or "t" in path, path
    else:
        if "abi3" not in path and "cpython" in path:
            pytest.skip(f"source/local GIL build: {path}")
        assert "abi3" in path, path


def test_status_string_and_constants():
    assert GJP_OK == 0
    assert GJP_ERR_BOGAERT_AB == 5
    assert "ok" in status_string(0).lower()
    assert "bogaert" in status_string(5).lower()


def test_auto_and_glr():
    x, w = gauss_jacobi(32, 0.0, 0.0)
    assert x.shape == (32,) and w.shape == (32,)
    assert np.all(w > 0) and np.all(np.isfinite(x))
    assert np.all((x > -1.0) & (x < 1.0))
    xg, wg = gauss_jacobi(40, 0.5, 0.5, method="glr")
    xa, wa = gauss_jacobi(40, 0.5, 0.5, method="algo665")
    assert np.allclose(np.sort(xg), np.sort(xa), atol=1e-12)
    assert np.allclose(wg[np.argsort(xg)], wa[np.argsort(xa)], atol=1e-11)


def test_glr_high_ab_matches_algo665():
    """Regression: large α,β used to skip the leftmost root in GLR."""
    n, a, b = 50, 5.0, 5.0
    xg, wg = gauss_jacobi(n, a, b, method="glr")
    xa, wa = gauss_jacobi(n, a, b, method="algo665")
    assert xg.shape == (n,)
    assert np.all(np.diff(np.sort(xg)) > 0)
    assert np.allclose(np.sort(xg), np.sort(xa), atol=1e-12)
    assert np.allclose(wg[np.argsort(xg)], wa[np.argsort(xa)], atol=1e-11)


def test_methods_rec_gw_sturm():
    for meth in ("rec", "gw", "algo665", "sturm"):
        x, w = gauss_jacobi(8, 0.0, 0.0, method=meth)
        assert x.shape == (8,) and np.all(w > 0)


def test_bogaert_ok_and_policy():
    x, w = gauss_jacobi(32, 0.0, 0.0, method="bogaert")
    assert x.shape == (32,) and np.all(w > 0)
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi(32, 0.5, 0.0, method="bogaert")
    assert "bogaert" in str(ei.value).lower() or "status=5" in str(ei.value)


def test_bad_npts():
    with pytest.raises(GaussJacobiError):
        gauss_jacobi(0, 0.0, 0.0)


def test_scale_polynomial_on_interval():
    """∫_0^3 (3-t) t^2 t^3 dt = 729/14, alpha=1, beta=2."""
    alpha, beta, a, b = 1.0, 2.0, 0.0, 3.0
    x, w = gauss_jacobi(8, alpha, beta, method="rec")
    xs, ws = gauss_jacobi_scale(x, w, alpha, beta, a, b)
    quad = float(np.dot(ws, xs**3))
    assert abs(quad - 729.0 / 14.0) < 1e-12
    assert abs(float(np.sum(ws)) - 27.0 / 4.0) < 1e-12
    assert np.all((xs > a) & (xs < b))
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, w, alpha, beta, b, a)
    assert "status=7" in str(ei.value)
    assert GJP_ERR_BOUNDS == 7


def _status_code(exc):
    text = str(exc)
    marker = "status="
    start = text.rfind(marker)
    assert start >= 0, text
    return int(text[start + len(marker) :].split(")", 1)[0].split(",", 1)[0])


def test_scale_huge_and_nonfinite():
    """Huge finite intervals stay finite; inf bounds and NaN parameters are errors."""
    x, w = gauss_jacobi(4, 0.0, 0.0, method="rec")
    for a, b in ((1.0e308, 1.7e308), (-1.0e308, 1.0e308)):
        xs, ws = gauss_jacobi_scale(x, w, 0.0, 0.0, a, b)
        assert np.all(np.isfinite(xs)) and np.all(np.isfinite(ws))
        assert np.all((xs > a) & (xs < b))
        assert np.all(ws > 0)
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, w, 0.0, 0.0, 0.0, np.inf)
    assert _status_code(ei.value) == GJP_ERR_BOUNDS
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, w, np.nan, 0.0, 0.0, 1.0)
    assert _status_code(ei.value) == GJP_ERR_ALPHA
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, w, 0.0, np.nan, 0.0, 1.0)
    assert _status_code(ei.value) == GJP_ERR_BETA


def test_scale_endpoints_and_weight_overflow():
    """Endpoints stay inside [a, b]. An overflowing weight product is status 7."""
    xs, _ws = gauss_jacobi_scale([-1.0], [1.0], 0.0, 0.0, 1.0, 1.0e16)
    assert xs[0] == 1.0
    xn = math.nextafter(-1.0, 0.0)
    xs, _ws = gauss_jacobi_scale([xn], [1.0], 0.0, 0.0, 1.0, 1.0e16)
    assert math.isfinite(float(xs[0])) and 1.0 < float(xs[0]) < 1.0e16
    xs, _ws = gauss_jacobi_scale([-1.0], [1.0], 0.0, 0.0, 1.0e6, 1.0e20)
    assert xs[0] == 1.0e6
    hi = sys.float_info.max
    xs, _ws = gauss_jacobi_scale([1.0], [1.0], 0.0, 0.0, 1.0e308, hi)
    assert math.isfinite(float(xs[0])) and float(xs[0]) == hi
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale([0.0], [2.0], 0.0, 0.0, -1.0e308, 1.0e308)
    assert _status_code(ei.value) == GJP_ERR_BOUNDS
    assert "status=7" in str(ei.value)
    assert "invalid interval or scaled weights overflow" in str(ei.value)


def test_scale_cancelled_endpoints_and_weight_range():
    """Lease images inside (a, b) and a lone overflowing factor are not answers."""
    xs, _ws = gauss_jacobi_scale([-1.0], [0.0], 0.0, 0.0, -1.0, 1.0e16)
    assert xs[0] == -1.0
    xn = math.nextafter(-1.0, 0.0)
    xs, _ws = gauss_jacobi_scale([xn], [0.0], 0.0, 0.0, -1.0, 1.0e16)
    assert xs[0] < 0.0 and abs(float(xs[0]) + 0.44488848768742173) < 1e-15
    xs, _ws = gauss_jacobi_scale([1.0], [0.0], 0.0, 0.0, -1.0e16, 1.0)
    assert xs[0] == 1.0
    xp = math.nextafter(1.0, -1.0)
    xs, _ws = gauss_jacobi_scale([xp], [0.0], 0.0, 0.0, -1.0e16, 1.0)
    assert xs[0] > 0.0 and abs(float(xs[0]) - 0.44488848768742173) < 1e-15
    xs, _ws = gauss_jacobi_scale([-1.0], [0.0], 0.0, 0.0, -1.0e16, -1.0)
    assert xs[0] == -1.0e16
    xs, _ws = gauss_jacobi_scale([1.0], [0.0], 0.0, 0.0, -1.0e16, -1.0)
    assert xs[0] == -1.0
    xs, _ws = gauss_jacobi_scale([xn], [0.0], 0.0, 0.0, 0.0, 3.0)
    assert xs[0] == 1.6653345369377348e-16
    xs, _ws = gauss_jacobi_scale([1.0], [0.0], 0.0, 0.0, -1.0e308, 1.0)
    assert xs[0] == 1.0
    x, w = gauss_jacobi(4, 1022.0, 0.0, method="rec")
    _xs, ws = gauss_jacobi_scale(x, w, 1022.0, 0.0, 0.0, 0.8)
    assert np.all(np.isfinite(ws)) and np.all(ws > 1.0e-120) and np.all(ws < 1.0e-80)
    x, w = gauss_jacobi(8, 0.0, 0.0, method="rec")
    _xs, ws = gauss_jacobi_scale(x, w, 1023.0, 0.0, -2.0, 2.0)
    assert np.all(np.isfinite(ws)) and np.all(ws > 1.0e300)
    tiny = float(np.nextafter(0.0, 1.0))
    _xs, ws = gauss_jacobi_scale([0.0], [1.0], -0.25, -0.25, 0.0, tiny)
    assert math.isfinite(float(ws[0])) and 1.0e-170 < float(ws[0]) < 1.0e-150
    _xs, ws = gauss_jacobi_scale([0.0], [1.0], -0.9, -0.9, 0.0, tiny)
    assert math.isfinite(float(ws[0])) and float(ws[0]) > 1.0e250


def test_scale_alpha1023_unit_weight_overflows():
    """w=1 on [-2, 2] with alpha=1023 is 2**1024 and must be status 7."""
    x = np.array([0.25])
    w = np.array([1.0])
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, w, 1023.0, 0.0, -2.0, 2.0)
    assert _status_code(ei.value) == GJP_ERR_BOUNDS
    pred = math.nextafter(1.0, 0.0)
    _xs, ws = gauss_jacobi_scale([0.25], [pred], 1023.0, 0.0, -2.0, 2.0)
    assert float(ws[0]) == sys.float_info.max
