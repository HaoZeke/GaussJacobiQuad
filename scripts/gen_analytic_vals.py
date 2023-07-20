import argparse

from sympy.integrals.quadrature import gauss_jacobi


def main(n, alpha, beta, n_dig):
    roots, weights = gauss_jacobi(n=n, alpha=alpha, beta=beta, n_digits=n_dig)
    for i in range(n):
        sign = " " if roots[i] >= 0 else ""
        root_str = f"{float(roots[i]):23.17E}"
        weight_str = f"{float(weights[i].evalf()):23.17E}"
        print(f"Root:{sign} {root_str} Weight: {weight_str}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Gauss-Jacobi quadrature.")
    parser.add_argument("--npts", type=int, default=5, help="Degree of the polynomial.")
    parser.add_argument("--alpha", type=float, default=0, help="Alpha parameter.")
    parser.add_argument("--beta", type=float, default=12, help="Beta parameter.")
    parser.add_argument("--n_dig", type=int, default=15, help="Precision.")

    args = parser.parse_args()

    main(args.npts, args.alpha, args.beta, args.n_dig)
