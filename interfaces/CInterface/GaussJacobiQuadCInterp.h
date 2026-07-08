/**
 * GaussJacobiQuad — C/C++ interface
 *
 * Preferred entry: gauss_jacobi_rule_c (by-value scalars, optional method, status).
 * Legacy per-method wrappers remain for ABI stability.
 *
 * Methods: auto (NULL / "" / "auto"), rec, gw, algo665, algo665_dc, sturm, glr, bogaert
 * (bogaert: alpha=beta=0 and n>=21). Multi-image CAF is not exposed here.
 */
#ifndef GAUSSJACOBIQUADCINTERP_H
#define GAUSSJACOBIQUADCINTERP_H

#ifdef __cplusplus
extern "C" {
#endif

enum gjp_status {
    GJP_OK = 0,
    GJP_ERR_NPTS = 1,
    GJP_ERR_ALPHA = 2,
    GJP_ERR_BETA = 3,
    GJP_ERR_METHOD = 4,
    GJP_ERR_BOGAERT_AB = 5,
    GJP_ERR_BOGAERT_N = 6
};

/**
 * Single dispatch (implemented in Fortran via ISO_C_BINDING).
 * @param method  NULL, "", or "auto" → select_method_auto policy
 * @return GJP_OK or GJP_ERR_* (does not abort on policy errors)
 */
int gauss_jacobi_rule_c(int npts, double alpha, double beta, double *x, double *wts,
                        const char *method);

/** Human-readable status (static string). */
const char *gjp_status_string(int status);

/**
 * Legacy dispatcher: pointer scalars; prints and exits on error (old CLI path).
 */
void gauss_jacobi_c(int *npts, double *alpha, double *beta, double x[], double wts[],
                    const char *method);

void gauss_jacobi_rec_c(int *npts, double *alpha, double *beta, double x[], double wts[]);
void gauss_jacobi_gw_c(int *npts, double *alpha, double *beta, double x[], double wts[]);
void gauss_jacobi_algo665_c(int *npts, double *alpha, double *beta, double x[], double wts[]);

#ifdef __cplusplus
} /* extern "C" */

inline int gauss_jacobi_rule(int npts, double alpha, double beta, double *x, double *wts,
                             const char *method = nullptr) {
    return gauss_jacobi_rule_c(npts, alpha, beta, x, wts, method);
}
#endif

#endif /* GAUSSJACOBIQUADCINTERP_H */
