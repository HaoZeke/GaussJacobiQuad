"""ASV: wall time and SciPy parity for gauss_jacobi methods.

Runs against the installed CPython extension (ISO_C_BINDING path).
"""
from __future__ import annotations

import numpy as np

try:
    from gauss_jacobi_quad import gauss_jacobi
except ImportError:  # pragma: no cover
    gauss_jacobi = None

try:
    from scipy.special import roots_jacobi
except ImportError:  # pragma: no cover
    roots_jacobi = None


def _skip_if_missing():
    if gauss_jacobi is None:
        raise NotImplementedError("gauss_jacobi_quad not installed")


class TimeMethods:
    """Wall-clock for representative (n, α, β, method) rules."""

    timeout = 60
    repeat = (3, 10, 20.0)
    number = 1
    warmup_time = 0.5
    params = [
        [16, 64, 128],
        ["auto", "algo665", "glr", "sturm", "rec"],
        [(0.0, 0.0), (0.5, 0.5), (5.0, 5.0)],
    ]
    param_names = ["n", "method", "ab"]

    def setup(self, n, method, ab):
        _skip_if_missing()
        self.n = n
        self.method = None if method == "auto" else method
        self.alpha, self.beta = ab
        # bogaert only for Legendre n>=21 — skip invalid via setup raise
        if method == "bogaert" and (ab != (0.0, 0.0) or n < 21):
            raise NotImplementedError("bogaert out of regime")
        # warm one call
        gauss_jacobi(self.n, self.alpha, self.beta, method=self.method)

    def time_rule(self, n, method, ab):
        gauss_jacobi(self.n, self.alpha, self.beta, method=self.method)


class TimeBogaertLegendre:
    """Legendre asymptotics (bogaert) vs auto for large n."""

    timeout = 60
    repeat = (3, 10, 20.0)
    number = 1
    params = [32, 64, 128, 256]
    param_names = ["n"]

    def setup(self, n):
        _skip_if_missing()
        gauss_jacobi(n, 0.0, 0.0, method="bogaert")

    def time_bogaert(self, n):
        gauss_jacobi(n, 0.0, 0.0, method="bogaert")

    def time_auto(self, n):
        gauss_jacobi(n, 0.0, 0.0, method=None)


class TrackAccuracyVsScipy:
    """Track max |Δx| vs scipy.special.roots_jacobi (lower is better).

    Guards the GLR / high-α,β regime so regressions show up on PRs.
    """

    timeout = 60
    repeat = 1
    number = 1
    params = [
        [16, 50, 100],
        ["auto", "algo665", "glr", "sturm", "rec"],
        [(0.0, 0.0), (0.5, 0.5), (5.0, 5.0)],
    ]
    param_names = ["n", "method", "ab"]
    unit = "max|dx|"
    less_is_better = True

    def setup(self, n, method, ab):
        _skip_if_missing()
        if roots_jacobi is None:
            raise NotImplementedError("scipy not installed")
        self.n = n
        self.method = None if method == "auto" else method
        self.alpha, self.beta = ab
        self.ref_x, _ = roots_jacobi(n, ab[0], ab[1])

    def track_max_dx(self, n, method, ab):
        x, _ = gauss_jacobi(self.n, self.alpha, self.beta, method=self.method)
        ia = np.argsort(x)
        ir = np.argsort(self.ref_x)
        return float(np.max(np.abs(x[ia] - self.ref_x[ir])))
