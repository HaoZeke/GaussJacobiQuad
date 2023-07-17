from sympy.integrals.quadrature import gauss_jacobi
n = 5  # degree of the polynomial
alpha = 1  # alpha parameter
beta = 50  # beta parameter
n_dig = 15 # Precision
roots, weights = gauss_jacobi(n=n, alpha=alpha, beta=beta, n_digits=n_dig)
for i in range(n):
    sign = ' ' if roots[i] >= 0 else ''
    print(f"Root:{sign} {float(roots[i]):23.17E} Weight: {float(weights[i].evalf()):23.17E}")
