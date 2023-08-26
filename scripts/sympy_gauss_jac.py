# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# Date: 2023-08-26
# Commit: 954667c
# -----------------------------------------------------------------------------
# This code is part of the GaussJacobiQuad library, providing an efficient
# implementation for Gauss-Jacobi quadrature nodes and weights computation.
# -----------------------------------------------------------------------------
# END_HEADER
"""!
@brief This script computes Gauss-Jacobi quadrature roots and weights using SciPy.

@details
The script takes polynomial degree, alpha, beta, and number of digits
(precision) parameters to compute Gauss-Jacobi quadrature roots and weights. It
uses sympy.integrals.quadrature.gauss_jacobi for the calculations.

@author Rohit Goswami
@date 26-08-2023
"""
import argparse

from sympy.integrals.quadrature import gauss_jacobi


def main(n, alpha, beta, n_dig):
    roots, weights = gauss_jacobi(n=n, alpha=alpha, beta=beta, n_digits=n_dig)
    for idx, root in enumerate(roots):
        sign = " " if root >= 0 else ""
        root_str = f"{float(root):23.17E}"
        weight_str = f"{float(weights[idx].evalf()):23.17E}"
        print(f"Root:{sign} {root_str} Weight: {weight_str}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Gauss-Jacobi quadrature.")
    parser.add_argument("--npts", type=int, default=5, help="Degree of the polynomial.")
    parser.add_argument("--alpha", type=float, default=0, help="Alpha parameter.")
    parser.add_argument("--beta", type=float, default=12, help="Beta parameter.")
    parser.add_argument("--n_dig", type=int, default=15, help="Precision.")

    args = parser.parse_args()

    main(args.npts, args.alpha, args.beta, args.n_dig)
