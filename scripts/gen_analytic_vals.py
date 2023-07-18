from sympy.integrals.quadrature import gauss_jacobi

n = 5  # degree of the polynomial
alpha = 0  # alpha parameter
beta = 12  # beta parameter
n_dig = 15  # Precision
roots, weights = gauss_jacobi(n=n, alpha=alpha, beta=beta, n_digits=n_dig)
for i in range(n):
    sign = " " if roots[i] >= 0 else ""
    root_str = f"{float(roots[i]):23.17E}"
    weight_str = f"{float(weights[i].evalf()):23.17E}"
    print(f"Root:{sign} {root_str} Weight: {weight_str}")
