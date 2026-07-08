#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "GaussJacobiQuadCInterp.h"

int main(int argc, char *argv[]) {
    int n_points;
    double alpha, beta;
    const char *method;
    double *x, *wts;
    int i, st;
    double wsum;

    if (argc == 4) {
        method = "auto";
    } else if (argc == 5) {
        method = argv[4];
    } else {
        fprintf(stderr, "Usage: %s <n_points> <alpha> <beta> [method]\n", argv[0]);
        fprintf(stderr, "  method: auto|rec|gw|algo665|algo665_dc|sturm|glr|bogaert\n");
        fprintf(stderr, "  omit method or pass auto for select_method_auto\n");
        return EXIT_FAILURE;
    }

    n_points = atoi(argv[1]);
    alpha = atof(argv[2]);
    beta = atof(argv[3]);

    x = (double *)malloc((size_t)n_points * sizeof(double));
    wts = (double *)malloc((size_t)n_points * sizeof(double));
    if (!x || !wts) {
        fprintf(stderr, "malloc failed\n");
        return EXIT_FAILURE;
    }

    st = gauss_jacobi_rule_c(n_points, alpha, beta, x, wts, method);
    if (st != GJP_OK) {
        fprintf(stderr, "gauss_jacobi_rule_c failed: %s (%d)\n", gjp_status_string(st), st);
        free(x);
        free(wts);
        return EXIT_FAILURE;
    }

    wsum = 0.0;
    for (i = 0; i < n_points; ++i) {
        printf("Root: %e, Weight: %e\n", x[i], wts[i]);
        wsum += wts[i];
    }
    fprintf(stderr, "method=%s n=%d sum(w)=%.17e\n", method, n_points, wsum);

    free(x);
    free(wts);
    return EXIT_SUCCESS;
}
