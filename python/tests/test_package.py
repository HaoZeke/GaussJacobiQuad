"""Smoke tests for the installable package (needs libgjp_cinterp)."""
import numpy as np
import pytest

from gauss_jacobi_quad import (
    GJP_ERR_BOGAERT_AB,
    GaussJacobiError,
    __version__,
    gauss_jacobi,
)


def test_version():
    assert __version__ == "0.2.0"


def test_auto_and_glr():
    x, w = gauss_jacobi(32, 0.0, 0.0)
    assert x.shape == (32,) and w.shape == (32,)
    assert np.all(w > 0) and np.all(np.isfinite(x))
    xg, wg = gauss_jacobi(40, 0.5, 0.5, method="glr")
    xa, wa = gauss_jacobi(40, 0.5, 0.5, method="algo665")
    assert np.allclose(np.sort(xg), np.sort(xa), atol=1e-12)


def test_bogaert_policy():
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi(32, 0.5, 0.0, method="bogaert")
    assert ei.value.status == GJP_ERR_BOGAERT_AB
