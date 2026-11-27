"""mpmath reference for gauss_jacobi_scale_c on wide and ordinary intervals."""
from __future__ import annotations

import math
import os
import random
import struct
import sys
from fractions import Fraction
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


def _floor_log2_ratio(n, d):
    e = n.bit_length() - d.bit_length()

    def ge(exp):
        if exp >= 0:
            return n >= (d << exp)
        return (n << (-exp)) >= d

    if not ge(e):
        e -= 1
    elif ge(e + 1):
        e += 1
    return e


def _round_fraction(frac):
    """Binary64 round-ties-to-even, including subnormals and overflow."""
    if frac == 0:
        return 0.0
    sign = -1.0 if frac < 0 else 1.0
    n = abs(frac.numerator)
    d = frac.denominator
    if n >= (d << 1024):
        return math.copysign(math.inf, sign)
    e = _floor_log2_ratio(n, d)
    if e < -1022:
        num = n << 1074
        quot, rem = divmod(num, d)
        if 2 * rem > d or (2 * rem == d and (quot & 1)):
            quot += 1
        return math.copysign(math.ldexp(quot, -1074), sign)
    shift = 52 - e
    if shift >= 0:
        num, den = n << shift, d
    else:
        num, den = n, d << (-shift)
    quot, rem = divmod(num, den)
    if 2 * rem > den or (2 * rem == den and (quot & 1)):
        quot += 1
    if quot == (1 << 53):
        quot >>= 1
        e += 1
    if e > 1023:
        return math.copysign(math.inf, sign)
    return math.copysign(math.ldexp(quot, e - 52), sign)


def _mpf_to_float(val):
    if val == 0:
        return math.copysign(0.0, -1.0 if mp.sign(val) < 0 else 1.0)
    sign, man, exp, _bc = val._mpf_
    frac = Fraction(int(man), 1) * Fraction(2) ** int(exp)
    if sign:
        frac = -frac
    return _round_fraction(frac)


def _cr_node(a, b, x):
    a = float(a)
    b = float(b)
    x = float(x)
    if not (math.isfinite(a) and math.isfinite(b) and math.isfinite(x)):
        return float(a + (x + 1.0) * (b - a) / 2.0)
    # Exact dyadic image. float(mpf) drops a subnormal sticky bit.
    fa, fb, fx = Fraction(a), Fraction(b), Fraction(x)
    return _round_fraction((fa * (1 - fx) + fb * (1 + fx)) / 2)


def _cr_weight(half, p, w):
    w = float(w)
    if w == 0.0 or p == 0:
        return w
    hh = half if isinstance(half, mp.mpf) else mp.mpf(float(half))
    if hh == 1:
        return w
    val = mp.power(mp.mpf(abs(w)), 1) * mp.power(hh, mp.mpf(p))
    if w < 0.0:
        val = -val
    return _mpf_to_float(val)


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
    ref = _cr_weight(float(half), 1023, w)
    assert _ulp(ws[0], ref) <= 1
    leg_w = [
        0.10122853629037609,
        0.22238103445337445,
        0.31370664587788744,
        0.3626837833783621,
    ]
    _x, ws = gauss_jacobi_scale(leg_w, leg_w, 1023.0, 0.0, -2.0, 2.0)
    for src, got in zip(leg_w, ws):
        assert math.isfinite(float(got)) and got != 0.0
        assert _ulp(got, _cr_weight(2.0, 1024, src)) <= 1
    half = mp.ldexp(1, -1075)
    tiny = float(np.nextafter(0.0, 1.0))
    _x, ws = gauss_jacobi_scale([0.0], [1.0], -0.25, -0.25, 0.0, tiny)
    ref = _cr_weight(half, 0.5, 1.0)
    assert _ulp(ws[0], ref) <= 1
    _x, ws = gauss_jacobi_scale([0.0], [1.0], -0.9, -0.9, 0.0, tiny)
    ref = _cr_weight(half, -0.8, 1.0)
    assert _ulp(ws[0], ref) <= 1
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


def _bits(x):
    return struct.unpack(">Q", struct.pack(">d", float(x)))[0]


def _from_bits(b):
    return struct.unpack(">d", int(b).to_bytes(8, "big"))[0]


def _first_nonfinite_weight(half, p, seed):
    lo = _bits(seed) - 4000
    hi = _bits(seed) + 4000
    if not math.isfinite(_cr_weight(half, p, _from_bits(lo))):
        raise AssertionError("seed window starts non-finite")
    if math.isfinite(_cr_weight(half, p, _from_bits(hi))):
        raise AssertionError("seed window stays finite")
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if math.isfinite(_cr_weight(half, p, _from_bits(mid))):
            lo = mid
        else:
            hi = mid
    return hi


def test_subnormal_sticky_node_rounds_up():
    """2**-1075 + 2**-2149 is above the tie and rounds to 2**-1074."""
    a = 0.0
    b = 2.0**-1074
    x = 2.0**-1074
    ref = _cr_node(a, b, x)
    assert ref == 2.0**-1074
    assert float(mp.mpf(a) + (mp.mpf(x) + 1) * (mp.mpf(b) - mp.mpf(a)) / 2) == 0.0
    got, _w = gauss_jacobi_scale([x], [0.0], 0.0, 0.0, a, b)
    assert got[0] == ref


def test_weight_overflow_both_sides():
    """Status 7 iff the correctly rounded product is non-finite."""
    cases = [
        (2.0, 1024, 1023.0, 0.0, -2.0, 2.0, 1.0, -8, 260),
        (2.0, 1023, 1022.0, 0.0, -2.0, 2.0, 2.0, -8, 16),
        (4.0, 512, 511.0, 0.0, -4.0, 4.0, 1.0, -8, 16),
        (
            1.25,
            3181,
            3180.0,
            0.0,
            0.0,
            2.5,
            float.fromhex("0x1.ed709438659dcp-1"),
            -40,
            16,
        ),
        (
            10.0,
            309,
            308.0,
            0.0,
            0.0,
            20.0,
            float.fromhex("0x1.702ae4d1fb5b5p-3"),
            -201,
            16,
        ),
        (
            1.5,
            1751,
            1750.0,
            0.0,
            0.0,
            3.0,
            float.fromhex("0x1.a8ce6778183bbp-1"),
            -8,
            8,
        ),
    ]
    pinned = {
        1.0: None,
        float.fromhex("0x1.0000000000100p+0"): None,
        float.fromhex("0x1.0000000000101p+0"): None,
        math.nextafter(1.0, 0.0): sys.float_info.max,
        float.fromhex("0x1.ed709438659dcp-1"): float.fromhex("0x1.ffffffffffff6p+1023"),
        float.fromhex("0x1.702ae4d1fb5b5p-3"): float.fromhex("0x1.fffffffffffd5p+1023"),
        float.fromhex("0x1.702ae4d1fb5cbp-3"): float.fromhex("0x1.ffffffffffff4p+1023"),
        float.fromhex("0x1.a8ce6778183bbp-1"): sys.float_info.max,
    }
    for half, p, alpha, beta, a, b, seed, lo_off, hi_off in cases:
        edge = _first_nonfinite_weight(half, p, seed)
        for step in range(lo_off, hi_off + 1):
            w = _from_bits(edge + step)
            x = np.array([0.25])
            ww = np.array([w])
            ref = _cr_weight(half, p, w)
            if not math.isfinite(ref):
                with pytest.raises(GaussJacobiError) as ei:
                    gauss_jacobi_scale(x, ww, alpha, beta, a, b)
                assert ei.value.status == GJP_ERR_BOUNDS
                assert x[0] == 0.25 and ww[0] == w
            else:
                got_x, got_w = gauss_jacobi_scale(x, ww, alpha, beta, a, b)
                assert _ulp(got_w[0], ref) <= 1, (
                    alpha,
                    w.hex(),
                    got_w[0].hex(),
                    ref.hex(),
                )
                assert math.isfinite(float(got_x[0]))
            if w in pinned and pinned[w] is not None and math.isfinite(ref):
                assert ref == pinned[w]
    # Named overflows through the ctypes wrapper.
    for w in (
        1.0,
        float.fromhex("0x1.0000000000100p+0"),
        float.fromhex("0x1.0000000000101p+0"),
        -1.0,
    ):
        x = np.array([0.25])
        ww = np.array([w])
        with pytest.raises(GaussJacobiError) as ei:
            gauss_jacobi_scale(x, ww, 1023.0, 0.0, -2.0, 2.0)
        assert ei.value.status == GJP_ERR_BOUNDS
        assert ww[0] == w and x[0] == 0.25
    x = np.array([0.25])
    ww = np.array([2.0])
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, ww, 1022.0, 0.0, -2.0, 2.0)
    assert ei.value.status == GJP_ERR_BOUNDS
    x = np.array([0.25])
    ww = np.array([1.0])
    with pytest.raises(GaussJacobiError) as ei:
        gauss_jacobi_scale(x, ww, 511.0, 0.0, -4.0, 4.0)
    assert ei.value.status == GJP_ERR_BOUNDS
