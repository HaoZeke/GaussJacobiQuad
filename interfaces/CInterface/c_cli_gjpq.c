// BEGIN_HEADER
// -----------------------------------------------------------------------------
// Gauss-Jacobi Quadrature Implementation
// Authors: Rohit Goswami <rgoswami[at]ieee.org>
// Source: GaussJacobiQuad Library
// License: MIT
// GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
// Date: 2023-08-28
// Commit: c442f77
// -----------------------------------------------------------------------------
// This code is part of the GaussJacobiQuad library, providing an efficient
// implementation for Gauss-Jacobi quadrature nodes and weights computation.
// -----------------------------------------------------------------------------
// To cite this software:
// Rohit Goswami (2023). HaoZeke/GaussJacobiQuad: v0.1.0.
// Zenodo: https://doi.org/10.5281/ZENODO.8285112
// ---------------------------------------------------------------------
// END_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GaussJacobiQuadCInterp.h"

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <n_points> <alpha> <beta> <method>\n", argv[0]);
        fprintf(stderr, "  n_points: Number of quadrature points (integer)\n");
        fprintf(stderr, "  alpha: Parameter alpha for Gauss-Jacobi quadrature (must be > -1)\n");
        fprintf(stderr, "  beta: Parameter beta for Gauss-Jacobi quadrature (must be > -1)\n");
        fprintf(stderr, "  method: Method to use for computation (supported: 'rec', 'gw', 'algo665')\n");
        return EXIT_FAILURE;
    }

    int n_points = atoi(argv[1]);
    double alpha = atof(argv[2]);
    double beta = atof(argv[3]);
    const char *method = argv[4];

    double *x = (double*) malloc((size_t) n_points * sizeof(double));
    double *wts = (double*) malloc((size_t) n_points * sizeof(double));

    gauss_jacobi_c(&n_points, &alpha, &beta, x, wts, method);

    for (int i = 0; i < n_points; ++i) {
        printf("Root: %e, Weight: %e\n", x[i], wts[i]);
    }

    free(x);
    free(wts);

    return EXIT_SUCCESS;
}
