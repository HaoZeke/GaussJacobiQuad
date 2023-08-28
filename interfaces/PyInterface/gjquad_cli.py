# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# Date: 2023-08-26
# Commit: 5dd0ffe
# -----------------------------------------------------------------------------
# This code is part of the GaussJacobiQuad library, providing an efficient
# implementation for Gauss-Jacobi quadrature nodes and weights computation.
# -----------------------------------------------------------------------------
# END_HEADER
"""!
@brief This script computes Gauss-Jacobi quadrature roots and weights using
GaussJacobiQuadPy.

@details
The script takes polynomial degree, alpha, beta, and method parameters to
compute Gauss-Jacobi quadrature roots and weights. It calls the GaussJacobiQuad
library for the calculations.

@author Rohit Goswami
@date 26-08-2023
"""
import argparse

import gjquadpy


def main(n, alpha, beta, meth):
    if meth == "rec":
        roots, weights = gjquadpy.gaussjacobiquadccompat.gauss_jacobi_rec_c(
            n, alpha, beta
        )
    elif meth == "gw":
        roots, weights = gjquadpy.gaussjacobiquadccompat.gauss_jacobi_gw_c(
            n, alpha, beta
        )
    elif meth == "algo665":
        roots, weights = gjquadpy.gaussjacobiquadccompat.gauss_jacobi_algo665_c(
            n, alpha, beta
        )
    # print(f"method: {meth}, npts: {n}, alpha: {alpha}, beta: {beta}, range: [-1, 1]")
    for idx, root in enumerate(roots):
        sign = " " if root >= 0 else ""
        root_str = f"{root:23.17E}"
        weight_str = f"{weights[idx]:23.17E}"
        print(f"Root:{sign} {root_str} Weight: {weight_str}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute Gauss-Jacobi quadrature via GaussJacobiQuad."
    )
    parser.add_argument("--npts", type=int, default=5, help="Degree of the polynomial.")
    parser.add_argument("--alpha", type=float, default=0, help="Alpha parameter.")
    parser.add_argument("--beta", type=float, default=12, help="Beta parameter.")
    parser.add_argument(
        "--meth",
        type=str,
        default="gw",
        help="Method.",
        choices=["gw", "rec", "algo665"],
    )

    args = parser.parse_args()

    main(args.npts, args.alpha, args.beta, args.meth)
