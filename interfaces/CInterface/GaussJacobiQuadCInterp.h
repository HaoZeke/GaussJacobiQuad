// BEGIN_HEADER
// -----------------------------------------------------------------------------
// Gauss-Jacobi Quadrature Implementation
// Authors: Rohit Goswami <rgoswami[at]ieee.org>
// Source: GaussJacobiQuad Library
// License: MIT
// GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
// Date: 2023-08-26
// Commit: 5dd0ffe
// -----------------------------------------------------------------------------
// This code is part of the GaussJacobiQuad library, providing an efficient
// implementation for Gauss-Jacobi quadrature nodes and weights computation.
// -----------------------------------------------------------------------------
// END_HEADER
#ifndef GAUSSJACOBIQUADCINTERP_H
#define GAUSSJACOBIQUADCINTERP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GaussJacobiQuadCInterp.h"

extern void gauss_jacobi_rec_c(int* npts, double* alpha, double* beta, double x[], double wts[]);
extern void gauss_jacobi_gw_c(int* npts, double* alpha, double* beta, double x[], double wts[]);

void gauss_jacobi_c(int* npts, double* alpha, double* beta, double x[], double wts[], const char* method) {
    if (*npts <= 0) {
        fprintf(stderr, "Error: npts must positive\n");
        exit(EXIT_FAILURE);
    }

    if (*alpha <= -1.0) {
        fprintf(stderr, "Error: alpha must be greater than -1\n");
        exit(EXIT_FAILURE);
    }

    if (*beta <= -1.0) {
        fprintf(stderr, "Error: beta must be greater than -1\n");
        exit(EXIT_FAILURE);
    }

    if (strcmp(method, "rec") == 0) {
        gauss_jacobi_rec_c(npts, alpha, beta, x, wts);
    } else if (strcmp(method, "gw") == 0) {
        gauss_jacobi_gw_c(npts, alpha, beta, x, wts);
    } else {
        fprintf(stderr, "Error: Unknown method specified: %s\n", method);
        fprintf(stderr, "Supported methods: 'rec', 'gw'\n");
        exit(EXIT_FAILURE);
    }
}

#endif /* GAUSSJACOBIQUADCINTERP_H */
