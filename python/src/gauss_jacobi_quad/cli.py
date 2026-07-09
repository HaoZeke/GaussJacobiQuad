"""CLI for GaussJacobiQuad (extension package)."""
from __future__ import annotations

import argparse
import sys

from gauss_jacobi_quad import GaussJacobiError, gauss_jacobi


def main() -> int:
    p = argparse.ArgumentParser(
        description="Gauss–Jacobi nodes/weights via C extension"
    )
    p.add_argument("--npts", type=int, default=5)
    p.add_argument("--alpha", type=float, default=0.0)
    p.add_argument("--beta", type=float, default=0.0)
    p.add_argument(
        "--meth",
        "--method",
        dest="method",
        default="auto",
        help="auto|rec|gw|algo665|algo665_dc|sturm|glr|bogaert",
    )
    args = p.parse_args()
    meth = None if args.method in (None, "", "auto") else args.method
    try:
        roots, weights = gauss_jacobi(args.npts, args.alpha, args.beta, method=meth)
    except GaussJacobiError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    for root, weight in zip(roots, weights):
        sign = " " if root >= 0 else ""
        print(f"Root:{sign} {root:23.17E} Weight: {weight:23.17E}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
