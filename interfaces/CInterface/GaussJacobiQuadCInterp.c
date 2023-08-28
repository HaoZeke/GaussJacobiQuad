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
#include "GaussJacobiQuadCInterp.h"

void gauss_jacobi_c(int* npts, double* alpha, double* beta, double x[], double wts[], const char* method) {
    if (*npts <= 0) {
        fprintf(stderr, "Error: npts must positive\n");
        exit(EXIT_FAILURE);
    }

    if (*alpha <= -1.0f) {
        fprintf(stderr, "Error: alpha must be greater than -1\n");
        exit(EXIT_FAILURE);
    }

    if (*beta <= -1.0f) {
        fprintf(stderr, "Error: beta must be greater than -1\n");
        exit(EXIT_FAILURE);
    }

    if (strcmp(method, "rec") == 0) {
        gauss_jacobi_rec_c(npts, alpha, beta, x, wts);
    } else if (strcmp(method, "gw") == 0) {
        gauss_jacobi_gw_c(npts, alpha, beta, x, wts);
    } else if (strcmp(method, "algo665") == 0) {
        gauss_jacobi_algo665_c(npts, alpha, beta, x, wts);
    } else {
        fprintf(stderr, "Error: Unknown method specified: %s\n", method);
        fprintf(stderr, "Supported methods: 'rec', 'gw', 'algo665'\n");
        exit(EXIT_FAILURE);
    }
}
