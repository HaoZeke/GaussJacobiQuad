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
    for pattern in ("libgjp_cinterp.so", "libgjp_cinterp.dylib"):
        found = next((ROOT / "bbdir").rglob(pattern), None)
        if found is not None:
            os.environ["GJP_CINTERP"] = str(found)
            break

pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("mpmath")

from accuracy_ranges import (  # noqa: E402
    REF_REL_TOL,
    REL_FLOOR,
    WEIGHT_TAU,
    format_report,
    library_rule,
    max_rel_err,
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
    assert lines[0] == ("method,n,alpha,beta,max_rel_node,max_rel_weight,ref,status")
    assert any(
        line.startswith("bogaert,8,1,0,") and line.endswith(",policy") for line in lines
    )


def test_weight_metric_floors_at_tau_times_max():
    ref = np.array([1.0, 1e-20])
    est = np.array([1.0, 2e-20])
    scored = max_rel_err(est, ref, weight=True)
    assert scored == pytest.approx(1e-20 / WEIGHT_TAU)
    assert scored < 1e-6


def test_scipy_gap_cells_use_mpmath():
    # SciPy weights sit ~5.6e-11 and ~2.3e-11 off mpmath on these cells.
    for alpha, beta in ((1.0, -0.9), (-0.9, 20.0)):
        _, _, kind = reference(64, alpha, beta)
        assert kind == "mpmath"


def test_gw_tiny_weight_does_not_set_the_score():
    rows = run_sweep(npts=(64,), alphas=(20.0,), betas=(-0.9,), methods=("gw",))
    assert rows[0].status == "ok"
    assert rows[0].ref == "mpmath"
    _x_m, w_m = mpmath_reference(64, 20.0, -0.9, dps=40)
    _x, w, status = library_rule(64, 20.0, -0.9, "gw")
    assert status == "ok"
    scored = max_rel_err(w, w_m, weight=True)
    assert scored == pytest.approx(rows[0].max_rel_weight, rel=1e-6, abs=0.0)
    # The smallest weight (~3e-27) is under WEIGHT_TAU * max|w|. Unfloored,
    # that component was the old headline (1.264288e-03). It must not set
    # the score.
    i_tiny = int(np.argmin(np.abs(w_m)))
    peak = float(np.max(np.abs(w_m)))
    assert abs(float(w_m[i_tiny])) < WEIGHT_TAU * peak
    floor = max(REL_FLOOR, WEIGHT_TAU * peak)
    per = np.abs(w - w_m) / np.maximum(np.abs(w_m), floor)
    assert float(per[i_tiny]) < scored
    # Same score on Linux OpenBLAS 0.3.26 and reference LAPACK 3.12
    # (gfortran 13.3.0), 1.442641896562202e-13, matching the measured CSV.
    # GitHub macOS-latest is 2.1245686687601543e-13. Twice that worst still
    # rejects a 10x inflated error. 1e-6 keeps the unfloored failure out.
    worst = 2.1245686687601543e-13
    bound = 2.0 * worst
    assert 10.0 * scored > bound
    assert scored < bound
    assert rows[0].max_rel_weight < bound
    assert scored < 1e-6


def test_beta20_weights_against_mpmath():
    rows = run_sweep(npts=(64,), alphas=(-0.9,), betas=(20.0,), methods=("algo665_dc",))
    assert rows[0].status == "ok"
    assert rows[0].ref == "mpmath"
    _x_m, w_m = mpmath_reference(64, -0.9, 20.0)
    _x, w, status = library_rule(64, -0.9, 20.0, "algo665_dc")
    assert status == "ok"
    scored = max_rel_err(w, w_m, weight=True)
    assert scored == pytest.approx(rows[0].max_rel_weight, rel=1e-6, abs=0.0)
    # Was 1.284004e+01 on a ~3e-27 weight. The corrected score moves with
    # the compiler and BLAS: 2.81570509e-10 in the measured CSV,
    # 1.778188454172291e-10 on GitHub ubuntu-latest, and
    # 1.7790996226302387e-10 with conda-forge OpenBLAS 0.3.34 and
    # gfortran 15.3.0. The bound is twice the worst of those.
    bound = 2.0 * 2.81570509e-10
    assert scored < bound
    assert rows[0].max_rel_weight < bound


def test_bogaert_n64_weight_within_rtol():
    rows = run_sweep(npts=(64,), alphas=(0.0,), betas=(0.0,), methods=("bogaert",))
    assert rows[0].status == "ok"
    assert rows[0].ref == "mpmath"
    _x_m, w_m = mpmath_reference(64, 0.0, 0.0, dps=40)
    _x, w, status = library_rule(64, 0.0, 0.0, "bogaert")
    assert status == "ok"
    scored = max_rel_err(w, w_m, weight=True)
    assert scored == pytest.approx(rows[0].max_rel_weight, rel=1e-6, abs=0.0)
    # Legendre n=64 stays inside REF_REL_TOL. Linux OpenBLAS 0.3.26 and
    # reference LAPACK 3.12 both score 2.739895740818578e-14 (the measured
    # CSV). GitHub macOS-latest is 1.7987101809954026e-14. Twice the worst
    # is still under REF_REL_TOL, and a 10x inflated error still fails.
    worst = 2.739895740818578e-14
    bound = 2.0 * worst
    assert 10.0 * scored > bound
    assert bound < REF_REL_TOL
    assert scored < bound
    assert rows[0].max_rel_weight < bound


def test_bogaert_legendre_domain():
    rows = run_sweep(npts=(32,), alphas=(1.0, 0.0), betas=(0.0,), methods=("bogaert",))
    by_alpha = {row.alpha: row for row in rows}
    assert by_alpha[1.0].status == "policy"
    assert by_alpha[0.0].status == "ok"


def test_unknown_method_and_bad_alpha_are_fail():
    bad_method = run_sweep(
        npts=(8,), alphas=(0.0,), betas=(0.0,), methods=("not-a-method",)
    )
    assert bad_method[0].status == "fail"
    bad_alpha = run_sweep(npts=(8,), alphas=(-1.0,), betas=(0.0,), methods=("rec",))
    assert bad_alpha[0].status == "fail"
