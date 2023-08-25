#ifndef GAUSSJACOBIQUADCINTERP_H
#define GAUSSJACOBIQUADCINTERP_H

#include <string.h>

extern void gauss_jacobi_rec_c(int* npts, double* alpha, double* beta, double x[], double wts[]);
extern void gauss_jacobi_gw_c(int* npts, double* alpha, double* beta, double x[], double wts[]);

void gauss_jacobi_c(int* npts, double* alpha, double* beta, double x[], double wts[], const char* method) {
    if (strcmp(method, "recurrence") == 0) {
        gauss_jacobi_rec_c(npts, alpha, beta, x, wts);
    } else if (strcmp(method, "gw") == 0) {
        gauss_jacobi_gw_c(npts, alpha, beta, x, wts);
    } else {
        // Handle unknown method here
    }
}

#endif /* GAUSSJACOBIQUADCINTERP_H */
