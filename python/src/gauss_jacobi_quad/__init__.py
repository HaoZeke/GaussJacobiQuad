# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# -----------------------------------------------------------------------------
# END_HEADER
"""GaussJacobiQuad Python package — CPython extension over ISO_C_BINDING C ABI.

Primary path: extension module ``gauss_jacobi_quad._core`` (multi-phase init,
free-threading aware) linked to ``gauss_jacobi_rule_c`` and
``gauss_jacobi_scale_c``.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from gauss_jacobi_quad import _core

__version__ = getattr(_core, "__version__", "0.2.6")

GJP_OK = _core.GJP_OK
GJP_ERR_NPTS = _core.GJP_ERR_NPTS
GJP_ERR_ALPHA = _core.GJP_ERR_ALPHA
GJP_ERR_BETA = _core.GJP_ERR_BETA
GJP_ERR_METHOD = _core.GJP_ERR_METHOD
GJP_ERR_BOGAERT_AB = _core.GJP_ERR_BOGAERT_AB
GJP_ERR_BOGAERT_N = _core.GJP_ERR_BOGAERT_N
GJP_ERR_BOUNDS = _core.GJP_ERR_BOUNDS

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


def gauss_jacobi_scale(nodes, weights, alpha: float, beta: float, a: float, b: float):
    """Map Gauss–Jacobi nodes and weights from [-1, 1] to [a, b].

    Nodes are one ties-to-even binary64 rounding of the real expression
    ``(b - a) / 2 * x + (a + b) / 2``. Weights are one ties-to-even
    binary64 rounding of ``w * H ** (alpha + beta + 1)``, where ``H`` is
    the binary64 half-width ``fl((b - a) / 2)``. A half-width that flushes
    to zero is treated as ``2**-1075``. ``x == -1`` is exactly ``a`` and
    ``x == 1`` is exactly ``b``. The mapped rule integrates
    ``(b - t) ** alpha * (t - a) ** beta * f(t)`` on ``[a, b]``.
    An invalid interval raises ``GaussJacobiError``.
    A non-finite scaled node or weight raises ``GaussJacobiError``.
    A weight factor that overflows on its own is not an error when the
    product with the weight is finite.
    """
    x_list, w_list = _core.scale(
        nodes, weights, float(alpha), float(beta), float(a), float(b)
    )
    return np.asarray(x_list, dtype=np.float64), np.asarray(w_list, dtype=np.float64)


# Back-compat alias
rule = gauss_jacobi

__all__ = [
    "gauss_jacobi",
    "gauss_jacobi_scale",
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
    "GJP_ERR_BOUNDS",
    "__version__",
]
