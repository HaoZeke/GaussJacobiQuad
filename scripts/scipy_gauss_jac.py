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
