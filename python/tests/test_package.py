"""Tests against the CPython extension (not ctypes)."""
from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from gauss_jacobi_quad import (
    GJP_ERR_BOGAERT_AB,
    GJP_OK,
    GaussJacobiError,
    __version__,
    gauss_jacobi,
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
    assert __version__ == "0.2.0"


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
