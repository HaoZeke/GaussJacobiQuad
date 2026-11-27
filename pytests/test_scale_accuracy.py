"""mpmath reference for gauss_jacobi_scale_c on wide and ordinary intervals."""
from __future__ import annotations

import math
import os
import random
import struct
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "interfaces" / "PyInterface"))
if "GJP_CINTERP" not in os.environ:
    for cand in (ROOT / "bbdir").rglob("libgjp_cinterp.so"):
        os.environ["GJP_CINTERP"] = str(cand)
        break

pytest.importorskip("gauss_jacobi_quad")
import gauss_jacobi_quad as _gjq  # noqa: E402

gauss_jacobi = _gjq.gauss_jacobi
gauss_jacobi_scale = _gjq.gauss_jacobi_scale
GaussJacobiError = _gjq.GaussJacobiError
GJP_ERR_BOUNDS = _gjq.GJP_ERR_BOUNDS

mp.mp.prec = 8000


def _ulp(a, b):
    a = float(a)
    b = float(b)
    if math.isnan(a) or math.isnan(b) or math.isinf(a) or math.isinf(b):
        return 0 if a == b else 10**18

    def ordered(x):
        i = int.from_bytes(struct.pack(">d", float(x)), "big")
        if i & (1 << 63):
            return (~i) & ((1 << 64) - 1)
        return i + (1 << 63)

    return abs(ordered(a) - ordered(b))


def _cr_node(a, b, x):
    aa, bb, xx = mp.mpf(float(a)), mp.mpf(float(b)), mp.mpf(float(x))
    return float(aa + (xx + 1) * (bb - aa) / 2)


def _cr_weight(half, p, w):
    if w == 0.0:
        return 0.0
    val = mp.exp(mp.log(mp.mpf(abs(w))) + mp.mpf(p) * mp.log(mp.mpf(half)))
    return float(mp.sign(w) * val)


def _intervals():
    huge = sys.float_info.max
    return [
        (0.0, 3.0),
        (-2.0, 5.0),
        (-1.0, 1.0),
        (-(2.0**52), 1.0),
        (-huge, huge),
        (1.0, 1.0e16),
        (1.0e6, 1.0e20),
        (-1.0e308, 1.0),
        (-(2.0**53), 1.0),
        (1.0e308, huge),
        (1.0, huge),
        (-huge, -1.0),
        (-1.0, 1.0e16),
        (-1.0e16, 1.0),
        (-1.0e16, -1.0),
        (0.1, 0.2),
        (1.0, 2.0),
        (-100.0, 100.0),
        (2.0**17, math.nextafter(2.0**17, math.inf)),
    ]


def test_scaled_nodes_within_one_ulp():
    rng = random.Random(1127)
    edges = [
        -1.0,
        1.0,
        0.0,
        0.3,
        -0.3,
        math.nextafter(-1.0, 0.0),
        math.nextafter(1.0, -1.0),
    ]
    worst = 0
    for a, b in _intervals():
        xs = edges + [rng.uniform(-1.0, 1.0) for _ in range(24)]
        got, _ = gauss_jacobi_scale(xs, [0.0] * len(xs), 0.0, 0.0, a, b)
        for x, t in zip(xs, got):
            u = _ulp(t, _cr_node(a, b, x))
            worst = max(worst, u)
            assert u <= 1, (a, b, x, t, _cr_node(a, b, x), u)
    assert worst <= 1


def test_weight_products_track_exact_power():
    w = 5.2872834730297891e304
    half = (mp.mpf("0.8") - mp.mpf(0)) / 2
    _x, ws = gauss_jacobi_scale([0.0], [w], 1022.0, 0.0, 0.0, 0.8)
    ref = float(mp.exp(mp.log(mp.mpf(w)) + 1023 * mp.log(half)))
    assert _ulp(ws[0], ref) <= 2048
    leg_w = [
        0.10122853629037609,
        0.22238103445337445,
        0.31370664587788744,
        0.3626837833783621,
    ]
    _x, ws = gauss_jacobi_scale(leg_w, leg_w, 1023.0, 0.0, -2.0, 2.0)
    for src, got in zip(leg_w, ws):
        assert math.isfinite(float(got)) and got != 0.0
        assert _ulp(got, _cr_weight(2.0, 1024, src)) <= 2048
    half = mp.ldexp(1, -1075)
    tiny = float(np.nextafter(0.0, 1.0))
    _x, ws = gauss_jacobi_scale([0.0], [1.0], -0.25, -0.25, 0.0, tiny)
    ref = float(mp.exp(mp.mpf("0.5") * mp.log(half)))
    assert _ulp(ws[0], ref) <= 2048
    _x, ws = gauss_jacobi_scale([0.0], [1.0], -0.9, -0.9, 0.0, tiny)
    ref = float(mp.exp(mp.mpf("-0.8") * mp.log(half)))
    assert _ulp(ws[0], ref) <= 2048
    nodes = np.array([0.0])
    weights = np.array([2.0])
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(nodes, weights, 0.0, 0.0, -1.0e308, 1.0e308)
    assert ei.value.status == GJP_ERR_BOUNDS
    assert nodes[0] == 0.0 and weights[0] == 2.0


def _review_cases():
    huge = sys.float_info.max
    xn = math.nextafter(-1.0, 0.0)
    return [
        (-3.3806359032229558e00, 3.4998833071125759e01, xn),
        (
            -4.9557863286209466e-02,
            6.8473170239949954e-01,
            float.fromhex("-0x1.bb093a01f6698p-1"),
        ),
        (-3.0885921864910043e-04, 6.2629754935323245e-02, xn),
        (1e300, huge, xn),
        (-(2.0**31), 2.0**85, xn),
    ]


def _random_intervals(rng, n=48):
    huge = sys.float_info.max
    out = [
        (1e300, huge),
        (1e308, huge),
        (-huge, -1e300),
        (-(2.0**31), 2.0**85),
        (-(2.0**20), 2.0**90),
        (1.0, math.nextafter(1.0, math.inf)),
    ]
    while len(out) < n:
        kind = rng.randrange(5)
        if kind == 0:
            a = rng.uniform(-1e4, 1e4)
            b = a
            for _ in range(rng.randrange(1, 8)):
                b = math.nextafter(b, math.inf)
        elif kind == 1:
            a = math.ldexp(rng.random() + 0.25, rng.randrange(-30, 40))
            b = a * math.ldexp(1.0, rng.randrange(20, 80))
        elif kind == 2:
            a = -rng.uniform(1e-6, 1e10)
            b = rng.uniform(1e-6, 1e10)
        elif kind == 3:
            a = rng.choice([1e300, 1e308, -1e308, math.ldexp(1.0, 900)])
            b = huge
        else:
            a = -math.ldexp(1.0, rng.randrange(8, 400))
            b = math.ldexp(1.0, rng.randrange(8, 900))
        if math.isfinite(a) and math.isfinite(b) and a < b:
            out.append((float(a), float(b)))
    return out


def test_review_regressions_within_one_ulp():
    """Named counterexamples through gauss_jacobi_scale_c."""
    worst = 0
    for a, b, x in _review_cases():
        got, _ = gauss_jacobi_scale([x], [0.0], 0.0, 0.0, a, b)
        u = _ulp(got[0], _cr_node(a, b, x))
        worst = max(worst, u)
        assert u <= 1, (a, b, x, got[0], _cr_node(a, b, x), u)
    ends, _ = gauss_jacobi_scale(
        [-1.0, 1.0], [0.0, 0.0], 0.0, 0.0, -(2.0**31), 2.0**85
    )
    assert ends[0] == -(2.0**31) and ends[1] == 2.0**85
    _xs, ws = gauss_jacobi_scale(
        [0.0], [2.0], sys.float_info.max, sys.float_info.max, 0.0, 2.0
    )
    assert ws[0] == 2.0
    print(f"MAX_ULP_REVIEW {worst}")


def test_random_nodes_within_one_ulp():
    """Seeded intervals, every quadrature node, via gauss_jacobi_scale_c."""
    rng = random.Random(112710)
    nodes, _w = gauss_jacobi(12, 0.0, 0.0, method="rec")
    probes = [float(v) for v in nodes]
    probes += [
        -1.0,
        1.0,
        math.nextafter(-1.0, 0.0),
        math.nextafter(1.0, -1.0),
    ]
    worst = 0
    for a, b in _random_intervals(rng):
        got, _ = gauss_jacobi_scale(probes, [0.0] * len(probes), 0.0, 0.0, a, b)
        for x, t in zip(probes, got):
            u = _ulp(t, _cr_node(a, b, x))
            worst = max(worst, u)
            assert u <= 1, (a, b, x, t, _cr_node(a, b, x), u)
    print(f"MAX_ULP_RANDOM {worst}")
