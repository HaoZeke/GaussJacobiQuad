"""GaussJacobiQuad Python package — CPython extension over ISO_C_BINDING C ABI.

Primary path: extension module ``gauss_jacobi_quad._core`` (multi-phase init,
free-threading aware) linked to ``gauss_jacobi_rule_c``.
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

import numpy as np

from gauss_jacobi_quad import _core

__version__ = getattr(_core, "__version__", "0.2.5")

GJP_OK = _core.GJP_OK
GJP_ERR_NPTS = _core.GJP_ERR_NPTS
GJP_ERR_ALPHA = _core.GJP_ERR_ALPHA
GJP_ERR_BETA = _core.GJP_ERR_BETA
GJP_ERR_METHOD = _core.GJP_ERR_METHOD
GJP_ERR_BOGAERT_AB = _core.GJP_ERR_BOGAERT_AB
GJP_ERR_BOGAERT_N = _core.GJP_ERR_BOGAERT_N

GaussJacobiError = _core.GaussJacobiError

status_string = _core.status_string


def gauss_jacobi(
    npts: int,
    alpha: float,
    beta: float,
    method: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Gauss–Jacobi nodes and weights via the C extension.

    Parameters
    ----------
    npts, alpha, beta
        Rule size and Jacobi parameters (alpha, beta > -1).
    method
        None / \"auto\" for select_method_auto; or rec, gw, algo665, algo665_dc,
        sturm, glr, bogaert (Legendre only).

    Returns
    -------
    nodes, weights : ndarray of shape (npts,)
    """
    x_list, w_list = _core.rule(int(npts), float(alpha), float(beta), method)
    return np.asarray(x_list, dtype=np.float64), np.asarray(w_list, dtype=np.float64)


# Back-compat alias
rule = gauss_jacobi

__all__ = [
    "gauss_jacobi",
    "rule",
    "GaussJacobiError",
    "status_string",
    "GJP_OK",
    "GJP_ERR_NPTS",
    "GJP_ERR_ALPHA",
    "GJP_ERR_BETA",
    "GJP_ERR_METHOD",
    "GJP_ERR_BOGAERT_AB",
    "GJP_ERR_BOGAERT_N",
    "__version__",
]
