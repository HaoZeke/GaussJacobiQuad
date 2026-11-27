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
    GJP_ERR_BOGAERT_N = 6,
    GJP_ERR_BOUNDS = 7
};

/**
 * Single dispatch (implemented in Fortran via ISO_C_BINDING).
 * @param method  NULL, "", or "auto" → select_method_auto policy
 * @return GJP_OK or GJP_ERR_* (does not abort on policy errors)
 */
int gauss_jacobi_rule_c(int npts, double alpha, double beta, double *x, double *wts,
                        const char *method);

/**
 * Map nodes and weights from [-1, 1] to [a, b] in place.
 * Nodes: t = (b-a)/2 * x + (a+b)/2.
 * Weights are multiplied by pow((b-a)/2, alpha+beta+1) so the rule integrates
 * (b-t)^alpha * (t-a)^beta * f(t) on [a, b].
 * @return GJP_OK, or GJP_ERR_NPTS / GJP_ERR_ALPHA / GJP_ERR_BETA / GJP_ERR_BOUNDS.
 *         Non-finite alpha or beta is ALPHA/BETA. Non-finite a or b, or a
 *         non-finite half-width or weight factor, is GJP_ERR_BOUNDS.
 *         NULL x or wts with npts > 0 is GJP_ERR_NPTS and does not write.
 */
int gauss_jacobi_scale_c(int npts, double alpha, double beta, double a, double b,
                         double *x, double *wts);

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
