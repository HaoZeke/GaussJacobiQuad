/* In-repo C consumer: drives shipped gauss_jacobi_rule_c vs algo665 reference. */
#include <float.h>
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
    /* scaled rule: ∫_0^3 (3-t) t^2 * t^3 dt = 729/14 */
    {
        enum { n = 8 };
        double x[8], w[8];
        double quad = 0.0, exact = 729.0 / 14.0;
        int i, st = gauss_jacobi_rule_c(n, 1.0, 2.0, x, w, "rec");
        if (st != GJP_OK) {
            fprintf(stderr, "scale setup failed: %s\n", gjp_status_string(st));
            rc = 1;
        } else {
            st = gauss_jacobi_scale_c(n, 1.0, 2.0, 0.0, 3.0, x, w);
            if (st != GJP_OK) {
                fprintf(stderr, "scale failed: %s\n", gjp_status_string(st));
                rc = 1;
            } else {
                for (i = 0; i < n; ++i)
                    quad += w[i] * x[i] * x[i] * x[i];
                printf("C scale [0,3] poly err=%.3e\n", fabs(quad - exact));
                if (fabs(quad - exact) > 1e-12) {
                    fprintf(stderr, "FAIL scale integral\n");
                    rc = 1;
                }
            }
        }
        st = gauss_jacobi_scale_c(4, 1.0, 2.0, 3.0, 0.0, x, w);
        if (st != GJP_ERR_BOUNDS) {
            fprintf(stderr, "expected GJP_ERR_BOUNDS got %d\n", st);
            rc = 1;
        }
    }
    /* Huge finite intervals: finite values, or a nonzero status.
       Never status 0 with inf/NaN. Non-finite inputs are nonzero. */
    {
        enum { n = 4 };
        double x[4], w[4], y[4], z[4];
        double a, b;
        int i, st = gauss_jacobi_rule_c(n, 0.0, 0.0, x, w, "rec");
        if (st != GJP_OK) {
            fprintf(stderr, "huge setup failed: %s\n", gjp_status_string(st));
            rc = 1;
        } else {
            for (i = 0; i < n; ++i) {
                y[i] = x[i];
                z[i] = w[i];
            }
            a = 1e308;
            b = 1.7e308;
            st = gauss_jacobi_scale_c(n, 0.0, 0.0, a, b, x, w);
            if (st != GJP_OK) {
                fprintf(stderr, "FAIL [1e308, 1.7e308] status %d\n", st);
                rc = 1;
            } else {
                for (i = 0; i < n; ++i) {
                    if (!isfinite(x[i]) || !isfinite(w[i]) || !(x[i] > a && x[i] < b)) {
                        fprintf(stderr, "FAIL [1e308, 1.7e308] status 0 non-finite\n");
                        rc = 1;
                        break;
                    }
                }
            }
            a = -1e308;
            b = 1e308;
            st = gauss_jacobi_scale_c(n, 0.0, 0.0, a, b, y, z);
            if (st != GJP_OK) {
                fprintf(stderr, "FAIL [-1e308, 1e308] status %d\n", st);
                rc = 1;
            } else {
                for (i = 0; i < n; ++i) {
                    if (!isfinite(y[i]) || !isfinite(z[i]) || !(y[i] > a && y[i] < b)) {
                        fprintf(stderr, "FAIL [-1e308, 1e308] status 0 non-finite\n");
                        rc = 1;
                        break;
                    }
                }
            }
            st = gauss_jacobi_scale_c(n, 0.0, 0.0, 0.0, HUGE_VAL, x, w);
            if (st != GJP_ERR_BOUNDS) {
                fprintf(stderr, "FAIL +inf bound status %d\n", st);
                rc = 1;
            }
            st = gauss_jacobi_scale_c(n, nan(""), 0.0, 0.0, 1.0, x, w);
            if (st != GJP_ERR_ALPHA) {
                fprintf(stderr, "FAIL NaN alpha status %d\n", st);
                rc = 1;
            }
            st = gauss_jacobi_scale_c(n, 0.0, nan(""), 0.0, 1.0, x, w);
            if (st != GJP_ERR_BETA) {
                fprintf(stderr, "FAIL NaN beta status %d\n", st);
                rc = 1;
            }
            st = gauss_jacobi_scale_c(n, 100.0, 100.0, -1e200, 1e200, x, w);
            if (st != GJP_ERR_BOUNDS) {
                fprintf(stderr, "FAIL non-finite scale factor status %d\n", st);
                rc = 1;
            }
        }
    }
    /* Errors must not write. NULL pointers return GJP_ERR_NPTS. */
    {
        double x[4] = {1.0, 2.0, 3.0, 4.0};
        double w[4] = {5.0, 6.0, 7.0, 8.0};
        double xc[4], wc[4];
        int st;
        memcpy(xc, x, sizeof x);
        memcpy(wc, w, sizeof w);
        st = gauss_jacobi_scale_c(4, nan(""), 2.0, 0.0, 3.0, x, w);
        if (st != GJP_ERR_ALPHA || memcmp(x, xc, sizeof x) || memcmp(w, wc, sizeof w)) {
            fprintf(stderr, "FAIL NaN alpha wrote buffers st=%d\n", st);
            rc = 1;
        }
        st = gauss_jacobi_scale_c(4, 2.0, nan(""), 0.0, 3.0, x, w);
        if (st != GJP_ERR_BETA || memcmp(x, xc, sizeof x) || memcmp(w, wc, sizeof w)) {
            fprintf(stderr, "FAIL NaN beta wrote buffers st=%d\n", st);
            rc = 1;
        }
        st = gauss_jacobi_scale_c(4, 0.0, 0.0, -HUGE_VAL, 1.0, x, w);
        if (st != GJP_ERR_BOUNDS || memcmp(x, xc, sizeof x) || memcmp(w, wc, sizeof w)) {
            fprintf(stderr, "FAIL -inf endpoint wrote buffers st=%d\n", st);
            rc = 1;
        }
        st = gauss_jacobi_scale_c(4, 0.0, 0.0, 0.0, 1.0, NULL, w);
        if (st != GJP_ERR_NPTS || memcmp(w, wc, sizeof w)) {
            fprintf(stderr, "FAIL NULL x st=%d\n", st);
            rc = 1;
        }
        st = gauss_jacobi_scale_c(4, 0.0, 0.0, 0.0, 1.0, x, NULL);
        if (st != GJP_ERR_NPTS || memcmp(x, xc, sizeof x)) {
            fprintf(stderr, "FAIL NULL w st=%d\n", st);
            rc = 1;
        }
        st = gauss_jacobi_scale_c(4, 0.0, 0.0, 0.0, 1.0, NULL, NULL);
        if (st != GJP_ERR_NPTS) {
            fprintf(stderr, "FAIL NULL both st=%d\n", st);
            rc = 1;
        }
    }
    /* Lease sum misses endpoints; a finite factor times w can still overflow. */
    {
        double x[1], w[1], xn, xkeep, wkeep;
        int st;
        if (strcmp(gjp_status_string(GJP_ERR_BOUNDS),
                   "invalid interval or scaled weights overflow") != 0) {
            fprintf(stderr, "FAIL bounds string: %s\n", gjp_status_string(GJP_ERR_BOUNDS));
            rc = 1;
        }
        x[0] = -1.0;
        w[0] = 1.0;
        st = gauss_jacobi_scale_c(1, 0.0, 0.0, 1.0, 1e16, x, w);
        printf("C [1,1e16] x=-1 st=%d x=%.17g\n", st, x[0]);
        if (st != GJP_OK || !isfinite(x[0]) || x[0] != 1.0) {
            fprintf(stderr, "FAIL [1,1e16] x=-1\n");
            rc = 1;
        }
        xn = nextafter(-1.0, 0.0);
        x[0] = xn;
        w[0] = 1.0;
        st = gauss_jacobi_scale_c(1, 0.0, 0.0, 1.0, 1e16, x, w);
        printf("C [1,1e16] next st=%d x=%.17g\n", st, x[0]);
        if (st != GJP_OK || !isfinite(x[0]) || !(x[0] > 1.0 && x[0] < 1e16)) {
            fprintf(stderr, "FAIL [1,1e16] nextafter\n");
            rc = 1;
        }
        x[0] = -1.0;
        w[0] = 1.0;
        st = gauss_jacobi_scale_c(1, 0.0, 0.0, 1e6, 1e20, x, w);
        printf("C [1e6,1e20] x=-1 st=%d x=%.17g\n", st, x[0]);
        if (st != GJP_OK || !isfinite(x[0]) || x[0] != 1e6) {
            fprintf(stderr, "FAIL [1e6,1e20] x=-1\n");
            rc = 1;
        }
        x[0] = 1.0;
        w[0] = 1.0;
        st = gauss_jacobi_scale_c(1, 0.0, 0.0, 1e308, DBL_MAX, x, w);
        printf("C [1e308,DBL_MAX] x=1 st=%d x=%.17g\n", st, x[0]);
        if (st != GJP_OK || !isfinite(x[0]) || x[0] != DBL_MAX) {
            fprintf(stderr, "FAIL [1e308,DBL_MAX] x=1\n");
            rc = 1;
        }
        x[0] = 0.0;
        w[0] = 2.0;
        xkeep = x[0];
        wkeep = w[0];
        st = gauss_jacobi_scale_c(1, 0.0, 0.0, -1e308, 1e308, x, w);
        printf("C weight overflow st=%d w=%.17g\n", st, w[0]);
        if (st != GJP_ERR_BOUNDS || x[0] != xkeep || w[0] != wkeep) {
            fprintf(stderr, "FAIL weight overflow wrote or status %d\n", st);
            rc = 1;
        }
    }
    if (rc == 0)
        printf("test_c_rule PASSED\n");
    else
        printf("test_c_rule FAILED\n");
    return rc;
}
