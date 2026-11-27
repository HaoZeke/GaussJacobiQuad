#!/usr/bin/env python3
"""Sweep shipped Gauss-Jacobi rules and record relative node/weight errors.

Drives ``gauss_jacobi`` (ctypes → libgjp_cinterp → Fortran) over a grid of
``npts`` and ``alpha, beta`` in ``(-1, 20]``. The reference is
``scipy.special.roots_jacobi`` when that rule is internally consistent, and a
high-precision mpmath Golub-Welsch rule otherwise. Writes a CSV and, with
``--plot``, the worst relative errors against ``npts``.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.special import eval_jacobi, roots_jacobi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "interfaces" / "PyInterface"))

# Cartesian product. alpha and beta are each drawn from AB, so the open
# interval (-1, 20] is sampled from just above -1 through 20.
NPTS: Tuple[int, ...] = (8, 32, 64)
AB: Tuple[float, ...] = (-0.9, 0.0, 1.0, 5.0, 20.0)
METHODS: Tuple[str, ...] = (
    "rec",
    "gw",
    "algo665",
    "algo665_dc",
    "sturm",
    "glr",
    "bogaert",
    "auto",
)

# |dx|/max(|x_ref|, REL_FLOOR) and the same for weights.
REL_FLOOR = 1e-30
# SciPy is treated as unreliable above either of these.
NEWTON_TOL = 1e-8
SUM_TOL = 1e-8
MP_DPS = 40

_API = None


@dataclass
class Row:
    method: str
    n: int
    alpha: float
    beta: float
    max_rel_node: float
    max_rel_weight: float
    ref: str
    status: str


def mu0_jacobi(alpha: float, beta: float) -> float:
    """Integral of (1-x)^alpha (1+x)^beta on [-1, 1]."""
    return math.exp(
        (alpha + beta + 1.0) * math.log(2.0)
        + math.lgamma(alpha + 1.0)
        + math.lgamma(beta + 1.0)
        - math.lgamma(alpha + beta + 2.0)
    )


def max_rel_err(est: np.ndarray, ref: np.ndarray) -> float:
    scale = np.maximum(np.abs(ref), REL_FLOOR)
    return float(np.max(np.abs(est - ref) / scale))


def scipy_reliable(
    n: int, alpha: float, beta: float, x: np.ndarray, w: np.ndarray
) -> bool:
    """True when a SciPy rule is monotone, positive, and consistent with P_n."""
    if x.shape != (n,) or w.shape != (n,):
        return False
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(w)):
        return False
    if np.any(x <= -1.0) or np.any(x >= 1.0) or np.any(w <= 0.0):
        return False
    if np.any(np.diff(x) <= 0.0):
        return False
    mu = mu0_jacobi(alpha, beta)
    if abs(float(np.sum(w)) - mu) > SUM_TOL * max(1.0, abs(mu)):
        return False
    try:
        poly = eval_jacobi(n, alpha, beta, x)
        deriv = 0.5 * (n + alpha + beta + 1.0) * eval_jacobi(
            n - 1, alpha + 1.0, beta + 1.0, x
        )
    except Exception:
        return False
    if not np.all(np.isfinite(poly)) or not np.all(np.isfinite(deriv)):
        return False
    if np.any(deriv == 0.0):
        return False
    newton = float(np.max(np.abs(poly / deriv)))
    return newton <= NEWTON_TOL


def mpmath_reference(
    n: int, alpha: float, beta: float, dps: int = MP_DPS
) -> Tuple[np.ndarray, np.ndarray]:
    """Golub-Welsch nodes and weights at ``dps`` decimal digits."""
    from mpmath import mp

    mp.dps = dps
    aa = mp.mpf(alpha)
    bb = mp.mpf(beta)
    diag: List = []
    off: List = []
    for k in range(n):
        kk = mp.mpf(k)
        if aa + bb == 0:
            an = (bb - aa) / (2 + aa + bb) if k == 0 else mp.mpf(0)
        elif k == 0:
            an = (bb - aa) / (2 + aa + bb)
        else:
            an = (bb**2 - aa**2) / ((2 * kk + aa + bb) * (2 * kk + aa + bb + 2))
        diag.append(an)
        if k == 0:
            continue
        bn = (
            mp.mpf(2)
            / (2 * kk + aa + bb)
            * mp.sqrt((kk + aa) * (kk + bb) / (2 * kk + aa + bb + 1))
        )
        if k != 1:
            bn *= mp.sqrt(kk * (kk + aa + bb) / (2 * kk + aa + bb - 1))
        off.append(bn)
    matrix = mp.matrix(n)
    for i, value in enumerate(diag):
        matrix[i, i] = value
    for i, value in enumerate(off):
        matrix[i, i + 1] = value
        matrix[i + 1, i] = value
    evals, evecs = mp.eigsy(matrix)
    order = sorted(range(n), key=lambda i: evals[i])
    mu0 = mp.power(2, aa + bb + 1) * mp.gamma(aa + 1) * mp.gamma(bb + 1) / mp.gamma(
        aa + bb + 2
    )
    xs = [evals[i] for i in order]
    ws = [mu0 * evecs[0, i] ** 2 for i in order]
    return (
        np.array([float(v) for v in xs], dtype=np.float64),
        np.array([float(v) for v in ws], dtype=np.float64),
    )


def reference(
    n: int, alpha: float, beta: float
) -> Tuple[np.ndarray, np.ndarray, str]:
    """SciPy nodes/weights, or mpmath when SciPy fails ``scipy_reliable``."""
    try:
        x, w = roots_jacobi(n, float(alpha), float(beta))
        x = np.asarray(x, dtype=np.float64)
        w = np.asarray(w, dtype=np.float64)
        order = np.argsort(x)
        x, w = x[order], w[order]
        if scipy_reliable(n, alpha, beta, x, w):
            return x, w, "scipy"
    except Exception:
        pass
    x, w = mpmath_reference(n, alpha, beta)
    return x, w, "mpmath"


def _library():
    global _API
    if _API is None:
        from gauss_jacobi_quad import GaussJacobiError, gauss_jacobi

        _API = (gauss_jacobi, GaussJacobiError)
    return _API


def library_rule(
    n: int, alpha: float, beta: float, method: str
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], str]:
    """One shipped rule. ``status`` is ``ok``, ``policy``, or ``fail``."""
    gauss_jacobi, err_type = _library()
    meth = None if method == "auto" else method
    try:
        x, w = gauss_jacobi(n, alpha, beta, method=meth)
    except err_type:
        return None, None, "policy"
    except Exception:
        return None, None, "fail"
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    if x.shape != (n,) or w.shape != (n,):
        return x, w, "fail"
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(w)):
        return x, w, "fail"
    order = np.argsort(x)
    return x[order], w[order], "ok"


def run_sweep(
    npts: Sequence[int] = NPTS,
    alphas: Sequence[float] = AB,
    betas: Sequence[float] = AB,
    methods: Sequence[str] = METHODS,
) -> List[Row]:
    cache: Dict[Tuple[int, float, float], Tuple[np.ndarray, np.ndarray, str]] = {}
    rows: List[Row] = []
    for n in npts:
        for alpha in alphas:
            for beta in betas:
                for method in methods:
                    x, w, status = library_rule(
                        int(n), float(alpha), float(beta), method
                    )
                    if status != "ok" or x is None or w is None:
                        rows.append(
                            Row(
                                method,
                                int(n),
                                float(alpha),
                                float(beta),
                                math.nan,
                                math.nan,
                                "",
                                status,
                            )
                        )
                        continue
                    key = (int(n), float(alpha), float(beta))
                    if key not in cache:
                        cache[key] = reference(key[0], key[1], key[2])
                    xref, wref, ref_name = cache[key]
                    rows.append(
                        Row(
                            method,
                            key[0],
                            key[1],
                            key[2],
                            max_rel_err(x, xref),
                            max_rel_err(w, wref),
                            ref_name,
                            "ok",
                        )
                    )
    return rows


def _fmt_err(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.8e}"


def write_csv(rows: Sequence[Row], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "method",
                "n",
                "alpha",
                "beta",
                "max_rel_node",
                "max_rel_weight",
                "ref",
                "status",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.method,
                    row.n,
                    f"{row.alpha:g}",
                    f"{row.beta:g}",
                    _fmt_err(row.max_rel_node),
                    _fmt_err(row.max_rel_weight),
                    row.ref,
                    row.status,
                ]
            )


def format_report(
    rows: Sequence[Row],
    npts: Sequence[int],
    alphas: Sequence[float],
    betas: Sequence[float],
) -> str:
    seen: Dict[Tuple[int, float, float], str] = {}
    for row in rows:
        if row.ref:
            seen[(row.n, row.alpha, row.beta)] = row.ref
    n_scipy = sum(name == "scipy" for name in seen.values())
    n_mp = sum(name == "mpmath" for name in seen.values())
    lines = [
        (
            f"GRID npts={','.join(str(n) for n in npts)} "
            f"alpha={','.join(f'{a:g}' for a in alphas)} "
            f"beta={','.join(f'{b:g}' for b in betas)} "
            f"rel_floor={REL_FLOOR:g} newton_tol={NEWTON_TOL:g} "
            f"sum_tol={SUM_TOL:g} mp_dps={MP_DPS} "
            f"ref_scipy={n_scipy} ref_mpmath={n_mp}"
        )
    ]
    methods: List[str] = []
    for row in rows:
        if row.method not in methods:
            methods.append(row.method)
    for method in methods:
        subset = [row for row in rows if row.method == method]
        ok = [row for row in subset if row.status == "ok"]
        n_fail = sum(row.status == "fail" for row in subset)
        n_policy = sum(row.status == "policy" for row in subset)
        if not ok:
            lines.append(
                f"RANGE method={method} ok=0 fail={n_fail} policy={n_policy} "
                "max_rel_node=nan max_rel_weight=nan"
            )
            continue
        worst_x = max(ok, key=lambda row: row.max_rel_node)
        worst_w = max(ok, key=lambda row: row.max_rel_weight)
        lines.append(
            f"RANGE method={method} ok={len(ok)} fail={n_fail} policy={n_policy} "
            f"max_rel_node={worst_x.max_rel_node:.6e} "
            f"at n={worst_x.n} alpha={worst_x.alpha:g} beta={worst_x.beta:g} "
            f"max_rel_weight={worst_w.max_rel_weight:.6e} "
            f"at n={worst_w.n} alpha={worst_w.alpha:g} beta={worst_w.beta:g}"
        )
    return "\n".join(lines)


def write_plot(rows: Sequence[Row], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    methods: List[str] = []
    for row in rows:
        if row.method not in methods:
            methods.append(row.method)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), sharex=True)
    for method in methods:
        ok = [row for row in rows if row.method == method and row.status == "ok"]
        if not ok:
            continue
        ns = sorted({row.n for row in ok})
        node = []
        weight = []
        for n in ns:
            group = [row for row in ok if row.n == n]
            node.append(max(row.max_rel_node for row in group))
            weight.append(max(row.max_rel_weight for row in group))
        node = [max(value, REL_FLOOR) for value in node]
        weight = [max(value, REL_FLOOR) for value in weight]
        axes[0].semilogy(ns, node, marker="o", label=method)
        axes[1].semilogy(ns, weight, marker="o", label=method)
    axes[0].set_ylabel("max relative error")
    axes[0].set_title("nodes")
    axes[1].set_title("weights")
    for axis in axes:
        axis.set_xlabel("npts")
        axis.grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="CSV path for per-cell errors")
    parser.add_argument("--plot", type=Path, help="Optional matplotlib figure")
    args = parser.parse_args(argv)
    try:
        _library()
    except Exception as exc:
        print(f"library load failed: {exc}", file=sys.stderr)
        return 1
    rows = run_sweep()
    report = format_report(rows, NPTS, AB, AB)
    print(report, flush=True)
    if args.csv is not None:
        write_csv(rows, args.csv)
        print(f"CSV {args.csv} bytes={args.csv.stat().st_size}", flush=True)
    if args.plot is not None:
        write_plot(rows, args.plot)
        print(f"PLOT {args.plot}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
