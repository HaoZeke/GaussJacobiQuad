"""pytest entry for SymPy accuracy matrix (real library path via ctypes)."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "interfaces" / "PyInterface"))

# Locate lib if meson build present
if "GJP_CINTERP" not in os.environ:
    for cand in (ROOT / "bbdir").rglob("libgjp_cinterp.so"):
        os.environ["GJP_CINTERP"] = str(cand)
        break
    for cand in Path.cwd().joinpath("bbdir").rglob("libgjp_cinterp.so") if Path.cwd().joinpath("bbdir").exists() else []:
        os.environ["GJP_CINTERP"] = str(cand)
        break

pytest.importorskip("sympy")
pytest.importorskip("numpy")

from sympy_accuracy import MATRIX, run_matrix, run_out_of_regime  # noqa: E402


def test_ctypes_available():
    from gauss_jacobi_quad import gauss_jacobi

    x, w = gauss_jacobi(8, 0.0, 0.0, method="algo665")
    assert x.shape == (8,)
    assert w.shape == (8,)


def test_full_sympy_matrix():
    rows = run_matrix(MATRIX, verbose=True)
    failed = [r for r in rows if not r.passed]
    assert not failed, "failures: " + "; ".join(
        f"{r.case}/{r.method} dx={r.max_dx} dw={r.max_dw} {r.note}" for r in failed
    )


def test_bogaert_out_of_regime():
    ok, text = run_out_of_regime(verbose=True)
    assert ok, text
