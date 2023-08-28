# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# Date: 2023-08-28
# Commit: c442f77
# -----------------------------------------------------------------------------
# This code is part of the GaussJacobiQuad library, providing an efficient
# implementation for Gauss-Jacobi quadrature nodes and weights computation.
# -----------------------------------------------------------------------------
# To cite this software:
# Rohit Goswami (2023). HaoZeke/GaussJacobiQuad: v0.1.0.
# Zenodo: https://doi.org/10.5281/ZENODO.8285112
# ---------------------------------------------------------------------
# END_HEADER
"""!
@brief This script computes Gauss-Jacobi quadrature roots and weights using SymPy.

@details
The script takes polynomial degree, alpha, and beta parameters to compute
Gauss-Jacobi quadrature roots and weights.  It uses scipy.special.roots_jacobi
for the calculations.

@author Rohit Goswami
@date 26-08-2023
"""
import argparse

from scipy.special import roots_jacobi


def main(n, alpha, beta):
    roots, weights = roots_jacobi(n=n, alpha=alpha, beta=beta)
    for idx, root in enumerate(roots):
        sign = " " if root >= 0 else ""
        root_str = f"{root:23.17E}"
        weight_str = f"{weights[idx]:23.17E}"
        print(f"Root:{sign} {root_str} Weight: {weight_str}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Gauss-Jacobi quadrature.")
    parser.add_argument("--npts", type=int, default=5, help="Degree of the polynomial.")
    parser.add_argument("--alpha", type=float, default=0, help="Alpha parameter.")
    parser.add_argument("--beta", type=float, default=12, help="Beta parameter.")

    args = parser.parse_args()

    main(args.npts, args.alpha, args.beta)
