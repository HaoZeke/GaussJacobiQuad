import subprocess

import numpy as np
import pytest


def extract_values(output_str):
    lines = output_str.strip().split("\n")
    roots = []
    weights = []
    for line in lines:
        _, root, _, weight = line.split()
        roots.append(float(root))
        weights.append(float(weight))
    return np.array(roots), np.array(weights)


def run_fortran(n, alpha, beta, method):
    result = subprocess.run(
        [
            "fpm",
            "run",
            "gjp_quad",
            "--",
            str(n),
            "{:.1f}".format(alpha),
            "{:.1f}".format(beta),
            method,
        ],
        stdout=subprocess.PIPE,
    )
    return extract_values(result.stdout.decode())


def run_python(n, alpha, beta):
    result = subprocess.run(
        [
            "python",
            "scripts/sympy_gauss_jac.py",
            "--npts",
            str(n),
            "--alpha",
            str(alpha),
            "--beta",
            str(beta),
        ],
        stdout=subprocess.PIPE,
    )
    return extract_values(result.stdout.decode())


@pytest.mark.parametrize(
    "n, alpha, beta, method",
    [
        (3, 1, 5, "recurrence"),
        (5, 2, 3, "recurrence"),
        pytest.param(
            10,
            0.0,
            30.0,
            "recurrence",
            marks=pytest.mark.xfail(reason="High beta values diverge"),
        ),
        (3, 1, 5, "gw"),
        (5, 2, 3, "gw"),
        pytest.param(
            10,
            0.0,
            30.0,
            "gw",
        ),
    ],
)
def test_gjp_quad_rec(n, alpha, beta, method):
    fortran_roots, fortran_weights = run_fortran(n, alpha, beta, method)
    python_roots, python_weights = run_python(n, alpha, beta)

    assert np.allclose(fortran_roots, python_roots, atol=1e-14)
    assert np.allclose(fortran_weights, python_weights, atol=1e-14)
