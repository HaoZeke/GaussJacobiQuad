/**
 * C helpers for GaussJacobiQuad.
 * Public gauss_jacobi_rule_c handles NULL method; Fortran core is gjp_rule_f.
 */
#include "GaussJacobiQuadCInterp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Fortran bind(C, name="gjp_rule_f") */
extern int gjp_rule_f(int npts, double alpha, double beta, double *x, double *wts,
                      const char *method_c);

/* Fortran bind(C, name="gjp_scale_f") */
extern int gjp_scale_f(int npts, double alpha, double beta, double a, double b, double *x,
                       double *wts);

int gauss_jacobi_rule_c(int npts, double alpha, double beta, double *x, double *wts,
                        const char *method) {
    const char *m = (method == NULL) ? "" : method;
    return gjp_rule_f(npts, alpha, beta, x, wts, m);
}

int gauss_jacobi_scale_c(int npts, double alpha, double beta, double a, double b, double *x,
                         double *wts) {
    if (npts <= 0)
        return GJP_ERR_NPTS;
    /* npts > 0 with a null buffer would segfault inside Fortran. */
    if (x == NULL || wts == NULL)
        return GJP_ERR_NPTS;
    if (!isfinite(alpha) || alpha <= -1.0)
        return GJP_ERR_ALPHA;
    if (!isfinite(beta) || beta <= -1.0)
        return GJP_ERR_BETA;
    if (!isfinite(a) || !isfinite(b) || !(a < b))
        return GJP_ERR_BOUNDS;
    return gjp_scale_f(npts, alpha, beta, a, b, x, wts);
}

const char *gjp_status_string(int status) {
    switch (status) {
    case GJP_OK:
        return "ok";
    case GJP_ERR_NPTS:
        return "npts must be positive";
    case GJP_ERR_ALPHA:
        return "alpha must be greater than -1";
    case GJP_ERR_BETA:
        return "beta must be greater than -1";
    case GJP_ERR_METHOD:
        return "unknown method (use auto|rec|gw|algo665|algo665_dc|sturm|glr|bogaert)";
    case GJP_ERR_BOGAERT_AB:
        return "bogaert requires alpha=beta=0 (Gauss-Legendre only)";
    case GJP_ERR_BOGAERT_N:
        return "bogaert requires n >= 21";
    case GJP_ERR_BOUNDS:
        return "require a < b";
    default:
        return "unknown status";
    }
}

void gauss_jacobi_c(int *npts, double *alpha, double *beta, double x[], double wts[],
                    const char *method) {
    int st;
    const char *m = method ? method : "auto";

    if (*npts <= 0) {
        fprintf(stderr, "Error: npts must be positive\n");
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

    st = gauss_jacobi_rule_c(*npts, *alpha, *beta, x, wts, m);
    if (st != GJP_OK) {
        fprintf(stderr, "Error: %s (status=%d, method=%s)\n", gjp_status_string(st), st, m);
        exit(EXIT_FAILURE);
    }
}
