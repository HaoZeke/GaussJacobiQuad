"""High-level Python bindings for GaussJacobiQuad (ctypes → C ABI).

Loads libgjp_cinterp (meson build product). Set GJP_CINTERP to the shared library
path if not found automatically.

Example::

    from gauss_jacobi_quad import gauss_jacobi
    x, w = gauss_jacobi(32, 0.0, 0.0)           # auto
    x, w = gauss_jacobi(40, 0.5, 0.5, method="glr")
"""

from __future__ import annotations

__version__ = "0.2.4"


import ctypes
import os
import sys
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np

# Status codes (match GaussJacobiQuadCInterp.h)
GJP_OK = 0
GJP_ERR_NPTS = 1
GJP_ERR_ALPHA = 2
GJP_ERR_BETA = 3
GJP_ERR_METHOD = 4
GJP_ERR_BOGAERT_AB = 5
GJP_ERR_BOGAERT_N = 6

_STATUS_MSG = {
    GJP_OK: "ok",
    GJP_ERR_NPTS: "npts must be positive",
    GJP_ERR_ALPHA: "alpha must be greater than -1",
    GJP_ERR_BETA: "beta must be greater than -1",
    GJP_ERR_METHOD: "unknown method",
    GJP_ERR_BOGAERT_AB: "bogaert requires alpha=beta=0 (Gauss-Legendre only)",
    GJP_ERR_BOGAERT_N: "bogaert requires n >= 21",
}


class GaussJacobiError(RuntimeError):
    """Raised when gauss_jacobi_rule_c returns a non-zero status."""

    def __init__(self, status: int, method: Optional[str] = None):
        self.status = status
        msg = _STATUS_MSG.get(status, f"unknown status {status}")
        if method is not None:
            msg = f"{msg} (method={method!r})"
        super().__init__(msg)


def _candidate_libs():
    env = os.environ.get("GJP_CINTERP")
    if env:
        yield Path(env)
    here = Path(__file__).resolve().parent
    root = here.parent.parent
    names = [
        "libgjp_cinterp.so",
        "libgjp_cinterp.dylib",
        "gjp_cinterp.dll",
    ]
    search = [
        root / "bbdir",
        root / "build" / "bbdir",
        Path.cwd() / "bbdir",
        Path.cwd(),
        here,
        root,
    ]
    # meson often nests
    for base in list(search):
        if base.is_dir():
            for p in base.rglob("libgjp_cinterp.so"):
                yield p
    for base in search:
        for name in names:
            yield base / name


def _load_lib() -> ctypes.CDLL:
    last = None
    for path in _candidate_libs():
        if path is None or not path.is_file():
            continue
        try:
            return ctypes.CDLL(str(path))
        except OSError as e:
            last = e
    raise FileNotFoundError(
        "Could not load libgjp_cinterp. Build with meson and set GJP_CINTERP "
        f"to the .so path. Last error: {last}"
    )


_lib = None


def _lib_api() -> ctypes.CDLL:
    global _lib
    if _lib is None:
        _lib = _load_lib()
        _lib.gauss_jacobi_rule_c.argtypes = [
            ctypes.c_int,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_char_p,
        ]
        _lib.gauss_jacobi_rule_c.restype = ctypes.c_int
    return _lib


def gauss_jacobi(
    npts: int,
    alpha: float,
    beta: float,
    method: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Gauss–Jacobi nodes and weights.

    Parameters
    ----------
    npts, alpha, beta
        Rule size and Jacobi parameters (alpha, beta > -1).
    method
        None / "auto" for select_method_auto; or rec, gw, algo665, algo665_dc,
        sturm, glr, bogaert (Legendre only).

    Returns
    -------
    nodes, weights : ndarray shape (npts,)

    Raises
    ------
    GaussJacobiError
        On invalid inputs or bogaert policy (does not abort the process).
    """
    n = int(npts)
    x = np.empty(n, dtype=np.float64)
    w = np.empty(n, dtype=np.float64)
    m = None if method is None else str(method).encode("ascii")
    if m == b"auto":
        m = b""  # empty → Fortran auto policy
    lib = _lib_api()
    st = lib.gauss_jacobi_rule_c(
        n,
        float(alpha),
        float(beta),
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        w.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        m if m is not None else b"",
    )
    if st != GJP_OK:
        raise GaussJacobiError(st, method)
    return x, w


# Back-compat alias
rule = gauss_jacobi
