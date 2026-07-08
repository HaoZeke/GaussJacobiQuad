/* In-repo C consumer: drives shipped gauss_jacobi_rule_c vs algo665 reference. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "GaussJacobiQuadCInterp.h"

static int check_case(int n, double alpha, double beta, const char *method, double tol) {
    double *x, *w, *xref, *wref;
    int i, st;
    double max_dx = 0.0, max_dw = 0.0, wsum = 0.0;

    x = (double *)malloc((size_t)n * sizeof(double));
    w = (double *)malloc((size_t)n * sizeof(double));
    xref = (double *)malloc((size_t)n * sizeof(double));
    wref = (double *)malloc((size_t)n * sizeof(double));
    if (!x || !w || !xref || !wref) {
        fprintf(stderr, "malloc failed\n");
        return 1;
    }

    st = gauss_jacobi_rule_c(n, alpha, beta, xref, wref, "algo665");
    if (st != GJP_OK) {
        fprintf(stderr, "ref algo665 failed: %s\n", gjp_status_string(st));
        return 1;
    }
    st = gauss_jacobi_rule_c(n, alpha, beta, x, w, method);
    if (st != GJP_OK) {
        fprintf(stderr, "method %s failed: %s\n", method, gjp_status_string(st));
        return 1;
    }
    for (i = 0; i < n; ++i) {
        double dx = fabs(x[i] - xref[i]);
        double dw = fabs(w[i] - wref[i]);
        if (dx > max_dx) max_dx = dx;
        if (dw > max_dw) max_dw = dw;
        wsum += w[i];
        if (!(w[i] > 0.0) || !(x[i] > -1.0 && x[i] < 1.0) || !isfinite(x[i]) || !isfinite(w[i])) {
            fprintf(stderr, "bad node/weight at %d: x=%g w=%g\n", i, x[i], w[i]);
            return 1;
        }
    }
    printf("C method=%s n=%d a=%.2f b=%.2f max|dx|=%.3e max|dw|=%.3e sum(w)=%.6e\n", method, n,
           alpha, beta, max_dx, max_dw, wsum);
    free(x);
    free(w);
    free(xref);
    free(wref);
    if (max_dx > tol || max_dw > tol * 10.0) {
        fprintf(stderr, "FAIL tol\n");
        return 1;
    }
    return 0;
}

int main(void) {
    int rc = 0;
    /* auto Legendre -> bogaert */
    rc |= check_case(32, 0.0, 0.0, "auto", 1e-10);
    rc |= check_case(32, 0.0, 0.0, "bogaert", 1e-10);
    rc |= check_case(40, 0.5, 0.5, "glr", 1e-10);
    rc |= check_case(40, 0.5, 0.5, "sturm", 1e-10);
    rc |= check_case(24, 0.5, 0.5, "rec", 1e-10);
    /* bogaert non-Legendre must return status, not crash */
    {
        double x[32], w[32];
        int st = gauss_jacobi_rule_c(32, 0.5, 0.0, x, w, "bogaert");
        printf("C bogaert non-Legendre status=%d (%s)\n", st, gjp_status_string(st));
        if (st != GJP_ERR_BOGAERT_AB) {
            fprintf(stderr, "expected GJP_ERR_BOGAERT_AB\n");
            rc = 1;
        }
    }
    if (rc == 0)
        printf("test_c_rule PASSED\n");
    else
        printf("test_c_rule FAILED\n");
    return rc;
}
