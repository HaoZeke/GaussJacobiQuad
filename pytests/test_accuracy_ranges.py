"""Smoke tests for the accuracy-range sweep (ctypes library path)."""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "interfaces" / "PyInterface"))

if "GJP_CINTERP" not in os.environ:
    for cand in (ROOT / "bbdir").rglob("libgjp_cinterp.so"):
        os.environ["GJP_CINTERP"] = str(cand)
        break

pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("mpmath")

from accuracy_ranges import (  # noqa: E402
    format_report,
    mpmath_reference,
    reference,
    run_sweep,
    scipy_reliable,
    write_csv,
)


def test_scipy_reliable_rejects_bad_nodes():
    x = np.linspace(-0.5, 0.5, 4)
    w = np.ones(4)
    assert not scipy_reliable(4, 0.0, 0.0, x, w)


def test_mpmath_matches_scipy_on_legendre():
    x_s, w_s, kind = reference(8, 0.0, 0.0)
    assert kind == "scipy"
    x_m, w_m = mpmath_reference(8, 0.0, 0.0)
    assert np.max(np.abs(x_s - x_m)) < 1e-12
    assert np.max(np.abs(w_s - w_m)) < 1e-12


def test_algo665_legendre_is_accurate():
    rows = run_sweep(npts=(8,), alphas=(0.0,), betas=(0.0,), methods=("algo665",))
    assert len(rows) == 1
    assert rows[0].status == "ok"
    assert rows[0].ref == "scipy"
    assert rows[0].max_rel_node < 1e-12
    assert rows[0].max_rel_weight < 1e-12


def test_gw_high_beta_is_finite():
    rows = run_sweep(npts=(8,), alphas=(0.0,), betas=(20.0,), methods=("gw",))
    assert rows[0].status == "ok"
    assert math.isfinite(rows[0].max_rel_node)
    assert math.isfinite(rows[0].max_rel_weight)


def test_bogaert_nonlegendre_is_policy():
    rows = run_sweep(npts=(8,), alphas=(1.0,), betas=(0.0,), methods=("bogaert",))
    assert rows[0].status == "policy"
    assert math.isnan(rows[0].max_rel_node)
    assert rows[0].ref == ""


def test_report_and_csv(tmp_path):
    npts = (8,)
    alphas = (0.0, 1.0)
    betas = (0.0,)
    rows = run_sweep(npts=npts, alphas=alphas, betas=betas, methods=("rec", "bogaert"))
    text = format_report(rows, npts, alphas, betas)
    assert "GRID npts=8 " in text
    assert "RANGE method=rec " in text
    assert "RANGE method=bogaert " in text
    path = tmp_path / "accuracy.csv"
    write_csv(rows, path)
    lines = path.read_text().splitlines()
    assert lines[0] == (
        "method,n,alpha,beta,max_rel_node,max_rel_weight,ref,status"
    )
    assert any(
        line.startswith("bogaert,8,1,0,") and line.endswith(",policy")
        for line in lines
    )
