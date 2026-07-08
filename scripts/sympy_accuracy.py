#!/usr/bin/env python3
"""SymPy accuracy matrix for shipped GaussJacobiQuad methods.

Reference: sympy.integrals.quadrature.gauss_jacobi (high n_digits).
Library path: ctypes gauss_jacobi_quad → libgjp_cinterp → Fortran kernels
(or fpm gjp_quad CLI fallback).

No hardcoded node/weight goldens — each cell compares live library output to
a fresh SymPy reference for that (n, α, β).

Regimes (see docs/METHODS.org auto policy + practical defaults):
  mild    — mid-n moderate α,β
  high_ab — elevated α or β (GW/algo665/sturm regime)
  legendre — α=β=0, n≥21 (bogaert + general methods)
  auto    — method=None / "auto" selection

Out-of-regime (not accuracy-claimed):
  bogaert with α≠0 or β≠0 → expect GaussJacobiError status 5
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "interfaces" / "PyInterface"))

# ---------------------------------------------------------------------------
# Case matrix
# ---------------------------------------------------------------------------

# abs tol for nodes; weights use max(abs_tol, rel_tol * max|w_ref|)
DEFAULT_NODE_TOL = 1e-12
DEFAULT_WEIGHT_TOL_ABS = 1e-11
DEFAULT_WEIGHT_TOL_REL = 1e-10
N_DIGITS = 28  # SymPy multiprecision digits (well above float64)


@dataclass(frozen=True)
class Case:
    name: str
    regime: str
    n: int
    alpha: float
    beta: float
    methods: Tuple[str, ...]
    node_tol: float = DEFAULT_NODE_TOL
    weight_tol_abs: float = DEFAULT_WEIGHT_TOL_ABS
    weight_tol_rel: float = DEFAULT_WEIGHT_TOL_REL
    # If True, method is optional string "auto"/None
    use_auto: bool = False


# All required single-image methods that produce a full rule
ALL_METHODS = (
    "rec",
    "gw",
    "algo665",
    "algo665_dc",
    "sturm",
    "glr",
    "bogaert",
)

MATRIX: Tuple[Case, ...] = (
    # Mild Jacobi mid-n — full set except bogaert (not Legendre)
    Case(
        name="mild_jacobi_n40",
        regime="mild",
        n=40,
        alpha=0.5,
        beta=0.5,
        methods=("rec", "gw", "algo665", "algo665_dc", "sturm", "glr"),
    ),
    # High-β mid-n — rec is known-weak at extreme β; include only robust GW stack
    Case(
        name="high_beta_n24",
        regime="high_ab",
        n=24,
        alpha=0.5,
        beta=12.0,
        methods=("gw", "algo665", "algo665_dc", "sturm"),
        node_tol=5e-12,
        weight_tol_abs=1e-10,
        weight_tol_rel=1e-9,
    ),
    # High-α
    Case(
        name="high_alpha_n24",
        regime="high_ab",
        n=24,
        alpha=10.0,
        beta=0.5,
        methods=("gw", "algo665", "algo665_dc", "sturm"),
        node_tol=5e-12,
        weight_tol_abs=1e-10,
        weight_tol_rel=1e-9,
    ),
    # Legendre n≥21 — includes bogaert
    Case(
        name="legendre_n32",
        regime="legendre",
        n=32,
        alpha=0.0,
        beta=0.0,
        methods=("rec", "gw", "algo665", "algo665_dc", "sturm", "glr", "bogaert"),
    ),
    # Auto Legendre → bogaert
    Case(
        name="auto_legendre_n48",
        regime="auto",
        n=48,
        alpha=0.0,
        beta=0.0,
        methods=("auto",),
        use_auto=True,
    ),
    # Auto high-β small-n → algo665
    Case(
        name="auto_high_beta_n16",
        regime="auto",
        n=16,
        alpha=0.5,
        beta=12.0,
        methods=("auto",),
        use_auto=True,
        node_tol=5e-12,
        weight_tol_abs=1e-10,
        weight_tol_rel=1e-9,
    ),
    # Mild larger n for glr-relevant auto path (n=64 mild → glr policy if not Legendre)
    Case(
        name="mild_jacobi_n64",
        regime="mild",
        n=64,
        alpha=0.5,
        beta=0.5,
        methods=("rec", "glr", "algo665", "sturm"),
        node_tol=2e-12,
        weight_tol_abs=2e-11,
    ),
)


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def sympy_reference(n: int, alpha: float, beta: float, n_digits: int = N_DIGITS):
    """High-precision SymPy reference nodes/weights.

    Important: pass *exact* alpha/beta into SymPy (via nsimplify). Passing raw
    Python floats makes sympy.integrals.quadrature.gauss_jacobi use an unstable
    floating Jacobi path that can return the wrong count of roots and nodes
    outside (-1,1). Integers/Rationals are required for a valid reference.
    """
    from sympy import nsimplify
    from sympy.integrals.quadrature import gauss_jacobi

    a = nsimplify(alpha, rational=True)
    b = nsimplify(beta, rational=True)
    roots, weights = gauss_jacobi(n=int(n), alpha=a, beta=b, n_digits=int(n_digits))
    # evalf before float() for full multiprecision digits
    x = np.array([float(r.evalf(n_digits)) for r in roots], dtype=np.float64)
    w = np.array([float(wt.evalf(n_digits)) for wt in weights], dtype=np.float64)
    if x.size != n or w.size != n:
        raise RuntimeError(
            f"SymPy returned {x.size} points for n={n} (alpha={a}, beta={b}); "
            "check that alpha/beta were exact rationals, not raw floats"
        )
    if np.any(x <= -1.0) or np.any(x >= 1.0):
        raise RuntimeError(f"SymPy reference has nodes outside (-1,1): min={x.min()} max={x.max()}")
    order = np.argsort(x)
    return x[order], w[order]


_lib_backend = None


def _ensure_ctypes():
    global _lib_backend
    if _lib_backend is not None:
        return _lib_backend
    try:
        from gauss_jacobi_quad import GaussJacobiError, gauss_jacobi

        # probe load
        gauss_jacobi(4, 0.0, 0.0, method="algo665")
        _lib_backend = ("ctypes", gauss_jacobi, GaussJacobiError)
        return _lib_backend
    except Exception as e:
        _lib_backend = ("fail", e, None)
        return _lib_backend


def library_rule(
    n: int, alpha: float, beta: float, method: Optional[str]
) -> Tuple[np.ndarray, np.ndarray]:
    """Call shipped library; returns sorted (x, w)."""
    kind, a, b = _ensure_ctypes()
    if kind == "ctypes":
        gauss_jacobi, GaussJacobiError = a, b
        meth = None if method in (None, "auto", "") else method
        x, w = gauss_jacobi(n, alpha, beta, method=meth)
        order = np.argsort(x)
        return x[order], w[order]
    # fpm CLI fallback
    return _fpm_rule(n, alpha, beta, method if method not in (None, "") else "auto")


def _fpm_rule(n: int, alpha: float, beta: float, method: str) -> Tuple[np.ndarray, np.ndarray]:
    cmd = [
        "fpm",
        "run",
        "--flag",
        "-fcoarray=single",
        "gjp_quad",
        "--",
        str(n),
        f"{alpha:.10f}",
        f"{beta:.10f}",
        method if method not in (None, "") else "auto",
    ]
    env = os.environ.copy()
    r = subprocess.run(
        cmd,
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        check=False,
    )
    if r.returncode != 0:
        raise RuntimeError(f"fpm gjp_quad failed: {r.stderr}\n{r.stdout}")
    xs, ws = [], []
    for line in r.stdout.splitlines():
        if "Root:" not in line:
            continue
        # Root:  ... Weight: ...
        parts = line.replace("Root:", " ").replace("Weight:", " ").split()
        if len(parts) >= 2:
            xs.append(float(parts[0]))
            ws.append(float(parts[1]))
    if len(xs) != n:
        raise RuntimeError(f"parsed {len(xs)} nodes, expected {n}: {r.stdout[:200]}")
    x = np.asarray(xs, dtype=np.float64)
    w = np.asarray(ws, dtype=np.float64)
    order = np.argsort(x)
    return x[order], w[order]


def weight_tol(w_ref: np.ndarray, abs_tol: float, rel_tol: float) -> float:
    scale = max(1.0, float(np.max(np.abs(w_ref))))
    return max(abs_tol, rel_tol * scale)


def compare(
    x: np.ndarray,
    w: np.ndarray,
    x_ref: np.ndarray,
    w_ref: np.ndarray,
    node_tol: float,
    weight_tol_abs: float,
    weight_tol_rel: float,
) -> Tuple[float, float, bool, str]:
    if x.shape != x_ref.shape:
        return np.inf, np.inf, False, f"shape {x.shape} vs {x_ref.shape}"
    max_dx = float(np.max(np.abs(x - x_ref)))
    max_dw = float(np.max(np.abs(w - w_ref)))
    wt = weight_tol(w_ref, weight_tol_abs, weight_tol_rel)
    ok = True
    notes = []
    if max_dx > node_tol:
        ok = False
        notes.append(f"node tol {node_tol:g}")
    if max_dw > wt:
        ok = False
        notes.append(f"weight tol {wt:g}")
    if np.any(w <= 0) or not np.all(np.isfinite(w)):
        ok = False
        notes.append("nonpositive/nonfinite weight")
    if np.any(x <= -1.0) or np.any(x >= 1.0) or not np.all(np.isfinite(x)):
        ok = False
        notes.append("node outside (-1,1) or nonfinite")
    if not np.all(np.diff(x) > 0):
        ok = False
        notes.append("nodes not strictly increasing")
    return max_dx, max_dw, ok, ";".join(notes) if notes else "ok"


def mu0_jacobi(alpha: float, beta: float) -> float:
    """∫_{-1}^1 (1-x)^α (1+x)^β dx = 2^{α+β+1} B(α+1,β+1)."""
    from math import exp, lgamma

    return exp(
        (alpha + beta + 1) * np.log(2.0)
        + lgamma(alpha + 1)
        + lgamma(beta + 1)
        - lgamma(alpha + beta + 2)
    )


# ---------------------------------------------------------------------------
# Run matrix
# ---------------------------------------------------------------------------

@dataclass
class Row:
    case: str
    regime: str
    method: str
    n: int
    alpha: float
    beta: float
    max_dx: float
    max_dw: float
    sumw_err: float
    passed: bool
    note: str


def run_matrix(
    cases: Sequence[Case] = MATRIX,
    verbose: bool = True,
) -> List[Row]:
    rows: List[Row] = []
    # Cache SymPy refs per (n,a,b)
    ref_cache = {}

    for case in cases:
        key = (case.n, case.alpha, case.beta)
        if key not in ref_cache:
            if verbose:
                print(f"# SymPy ref n={case.n} a={case.alpha} b={case.beta} digits={N_DIGITS}", flush=True)
            ref_cache[key] = sympy_reference(case.n, case.alpha, case.beta)
        x_ref, w_ref = ref_cache[key]
        mu0 = mu0_jacobi(case.alpha, case.beta)

        for method in case.methods:
            meth_label = method
            try:
                x, w = library_rule(
                    case.n,
                    case.alpha,
                    case.beta,
                    None if case.use_auto or method == "auto" else method,
                )
                max_dx, max_dw, ok, note = compare(
                    x,
                    w,
                    x_ref,
                    w_ref,
                    case.node_tol,
                    case.weight_tol_abs,
                    case.weight_tol_rel,
                )
                sumw_err = float(abs(np.sum(w) - mu0))
                if sumw_err > max(1e-10, 1e-9 * abs(mu0)):
                    ok = False
                    note = (note + ";sum(w)" if note != "ok" else "sum(w) vs mu0")
            except Exception as e:
                max_dx = max_dw = sumw_err = float("nan")
                ok = False
                note = f"EXCEPTION: {e}"
            row = Row(
                case=case.name,
                regime=case.regime,
                method=meth_label,
                n=case.n,
                alpha=case.alpha,
                beta=case.beta,
                max_dx=max_dx,
                max_dw=max_dw,
                sumw_err=sumw_err,
                passed=ok,
                note=note,
            )
            rows.append(row)
            status = "PASS" if ok else "FAIL"
            if verbose:
                print(
                    f"{status} case={case.name} regime={case.regime} method={meth_label} "
                    f"n={case.n} a={case.alpha} b={case.beta} "
                    f"max|dx|={max_dx:.3e} max|dw|={max_dw:.3e} sumw_err={sumw_err:.3e} ({note})",
                    flush=True,
                )
    return rows


def run_out_of_regime(verbose: bool = True) -> Tuple[bool, str]:
    """bogaert with α≠0 must not claim accuracy (error/status)."""
    kind, a, b = _ensure_ctypes()
    lines = []
    ok = True
    if kind == "ctypes":
        gauss_jacobi, GaussJacobiError = a, b
        try:
            gauss_jacobi(32, 0.5, 0.0, method="bogaert")
            ok = False
            lines.append("FAIL bogaert non-Legendre unexpectedly succeeded")
        except Exception as e:
            # expect GaussJacobiError with status 5
            st = getattr(e, "status", None)
            lines.append(f"bogaert non-Legendre raised {type(e).__name__} status={st}: {e}")
            if st is not None and st != 5:
                ok = False
                lines.append("FAIL expected status 5 (GJP_ERR_BOGAERT_AB)")
            else:
                lines.append("PASS out-of-regime bogaert policy")
    else:
        # fpm: expect nonzero exit
        r = subprocess.run(
            [
                "fpm",
                "run",
                "--flag",
                "-fcoarray=single",
                "gjp_quad",
                "--",
                "32",
                "0.5",
                "0.0",
                "bogaert",
            ],
            cwd=str(ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        lines.append(f"fpm bogaert non-Legendre exit={r.returncode}")
        if r.returncode == 0:
            ok = False
            lines.append("FAIL expected nonzero exit")
        else:
            lines.append("PASS out-of-regime bogaert policy (process error-stop)")
    text = "\n".join(lines)
    if verbose:
        print(text, flush=True)
    return ok, text


def format_table(rows: Sequence[Row]) -> str:
    lines = [
        "method | case | regime | n | alpha | beta | max|dx| | max|dw| | sumw_err | pass | note",
        "-------|------|--------|---|-------|------|---------|---------|----------|------|------",
    ]
    for r in rows:
        lines.append(
            f"{r.method} | {r.case} | {r.regime} | {r.n} | {r.alpha} | {r.beta} | "
            f"{r.max_dx:.3e} | {r.max_dw:.3e} | {r.sumw_err:.3e} | "
            f"{'PASS' if r.passed else 'FAIL'} | {r.note}"
        )
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", type=Path, help="Write markdown-ish table to path")
    p.add_argument("--spot", action="store_true", help="Only glr Legendre + bogaert + auto")
    p.add_argument("--out-of-regime", action="store_true", help="Only bogaert policy check")
    args = p.parse_args(argv)

    if args.out_of_regime:
        ok, _ = run_out_of_regime()
        return 0 if ok else 1

    cases = MATRIX
    if args.spot:
        cases = tuple(
            c
            for c in MATRIX
            if c.name in ("legendre_n32", "auto_legendre_n48", "mild_jacobi_n40")
        )
        # restrict methods for spot
        cases = (
            Case(
                name="spot_glr_legendre",
                regime="legendre",
                n=32,
                alpha=0.0,
                beta=0.0,
                methods=("glr", "bogaert"),
            ),
            Case(
                name="spot_auto",
                regime="auto",
                n=48,
                alpha=0.0,
                beta=0.0,
                methods=("auto",),
                use_auto=True,
            ),
        )

    kind, a, b = _ensure_ctypes()
    if kind != "ctypes":
        print(f"WARNING: ctypes unavailable ({a}); trying fpm fallback", file=sys.stderr)

    rows = run_matrix(cases)
    oor_ok, oor_text = run_out_of_regime()
    table = format_table(rows)
    print("\n=== SUMMARY TABLE ===")
    print(table)
    if args.table:
        args.table.write_text(table + "\n# out-of-regime\n" + oor_text + "\n")

    failed = [r for r in rows if not r.passed]
    if failed or not oor_ok:
        print(f"FAILED {len(failed)} accuracy cells; out_of_regime_ok={oor_ok}", file=sys.stderr)
        return 1
    print(f"ALL PASSED ({len(rows)} cells + out-of-regime)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
