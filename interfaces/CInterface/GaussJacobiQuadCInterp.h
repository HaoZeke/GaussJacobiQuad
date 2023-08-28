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

extern void gauss_jacobi_rec_c(int* npts, double* alpha, double* beta, double x[], double wts[]);
extern void gauss_jacobi_gw_c(int* npts, double* alpha, double* beta, double x[], double wts[]);
extern void gauss_jacobi_algo665_c(int* npts, double* alpha, double* beta, double x[], double wts[]);

void gauss_jacobi_c(int *npts, double *alpha, double *beta, double x[],
                    double wts[], const char *method);

#endif /* GAUSSJACOBIQUADCINTERP_H */
