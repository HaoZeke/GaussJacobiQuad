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
#ifndef GAUSSJACOBIQUADCINTERP_H
#define GAUSSJACOBIQUADCINTERP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void gauss_jacobi_rec_c(int* npts, double* alpha, double* beta, double x[], double wts[]);
extern void gauss_jacobi_gw_c(int* npts, double* alpha, double* beta, double x[], double wts[]);
extern void gauss_jacobi_algo665_c(int* npts, double* alpha, double* beta, double x[], double wts[]);

void gauss_jacobi_c(int *npts, double *alpha, double *beta, double x[],
                    double wts[], const char *method);

#endif /* GAUSSJACOBIQUADCINTERP_H */
