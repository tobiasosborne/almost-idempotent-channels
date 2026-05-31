/* aic_ucp_shim.h — flat-double ccall SHIM over the arb UCP / cb-norm cores
 * (bead aic-m24, increment 2a). The C substrate for the reusable
 * eta_idempotence(kraus) entry point: it exposes a FLAT-double ABI so the
 * (later, increment 2b) Julia package can `ccall` the certified arb Choi-diff
 * and the eig-free bracket WITHOUT any FLINT type (acb_mat_t / slong) crossing
 * the language boundary. Pure marshalling — the math lives in aic_ucp.c /
 * aic_cbnorm.c.
 *
 * WHY A SHIM (the ccall contract). Julia's `ccall` can pass int / double /
 * Ptr{Cdouble} cleanly, but NOT FLINT structs (acb_mat_t is an opaque array
 * type; slong is FLINT's machine-word int). So the public language boundary
 * uses int dimensions and Ptr{Cdouble} arrays only; this shim converts to/from
 * the internal aic_ucp_kraus / acb_mat at the boundary, runs the SAME cores the
 * C tests exercise (aic_ucp_compose + aic_ucp_choi_diff; aic_cbnorm_eigfree_*),
 * and writes the result back into flat double arrays the caller owns.
 *
 * LAYOUT (matches tests/fixtures_d24.inc.h and aic_ucp.h Convention-A):
 *   kraus_re[a*n*n + i*n + j], kraus_im[...]   r ops, each n x n (K_a: H->K, n=n)
 *   choi_re[p*N + q], choi_im[...]             N = n*n, row-major
 * with the Convention-A Choi C[i*n+a, j*n+b] = sum_x conj(K_x[i,a]) K_x[j,b]
 * (aic_ucp.h:23-42; the CONJUGATION is on the FIRST (i,a) factor — the shim must
 * preserve it through the double round-trip, guarded by the complex_qubit test).
 */
#ifndef AIC_UCP_SHIM_H
#define AIC_UCP_SHIM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Build J = Choi(Phi^2 - Phi) (Convention-A, n^2 x n^2) for the UCP self-map
 * Phi given by r Kraus operators, entirely through flat double arrays so the
 * routine is ccall-able from Julia (no acb_mat_t / FLINT types cross the
 * boundary). Internally: build aic_ucp_kraus (acb_set_d_d), Phi^2 =
 * aic_ucp_compose(phi,phi), J = aic_ucp_choi_diff(phi2, phi), then extract J
 * back to the double arrays (arf_get_d on the midpoints of real/imag parts).
 * `prec` is the arb working precision (e.g. 106 for headroom; 53 to cross-check
 * the double path). Returns 0 on success; fail-loud asserts on bad shape
 * (n>=1, r>=1). Caller owns all arrays (choi_* sized N*N, N=n*n). */
int aic_ucp_choi_diff_d(double *choi_re, double *choi_im,
                        const double *kraus_re, const double *kraus_im,
                        int n, int r, int prec);

/* Eig-free certified two-sided bracket on eta = ||Phi^2-Phi||_cb, returned as
 * two doubles: a RIGOROUS lower bound `*lo` and a RIGOROUS upper bound `*hi`
 * (so `*lo <= eta <= *hi`). Wraps aic_cbnorm_eigfree_ball (include/aic_cbnorm.h)
 * and converts the arb balls to doubles with the rounding that keeps the bracket
 * rigorous: `*lo = floor(arb_lower(lo_ball))`, `*hi = ceil(arb_upper(hi_ball))`
 * (arb_get_lbound_arf / arb_get_ubound_arf -> arf_get_d with ARF_RND_FLOOR /
 * ARF_RND_CEIL). No SDP, no eigendecomposition — Julia gets a certified bracket
 * via ccall with no solver. Returns 0 on success; fail-loud asserts on bad shape
 * (n>=1, r>=1). `lo`, `hi` are caller-owned scalars. */
int aic_cbnorm_eigfree_d(double *lo, double *hi,
                         const double *kraus_re, const double *kraus_im,
                         int n, int r, int prec);

/* TIGHT certified two-sided bracket on eta = ||Phi^2-Phi||_cb, returned as two
 * RIGOROUS doubles (*lo <= eta <= *hi). Flat-double ABI over aic_cbnorm_certify
 * (include/aic_cbnorm.h; docs/cbnorm_tight_certifier.md) for the 3c Julia ccall:
 * J and the two MOSEK feasible points (MAX-primal X,P,Q; MIN-dual Y0,Y1) are
 * passed as flat N x N [p*N+q] (re,im) arrays (N=n*n); the routine loads them to
 * acb_mat, runs the arb certifier, and converts the balls with the rigorous
 * rounding *lo = floor(arb_lower), *hi = ceil(arb_upper). The bracket is the
 * MOSEK-tight ball (gap ~ solver tol + arb radius), not the eig-free 2n-wide
 * one; it dispatches to the eig-free fallback near eta=0 or when the MOSEK points
 * cannot be restored. Returns 0 on success; fail-loud asserts on bad shape
 * (n>=1) or a null array. `lo`, `hi` are caller-owned scalars. */
int aic_cbnorm_certify_d(double *lo, double *hi,
                         const double *J_re, const double *J_im,
                         const double *X_re, const double *X_im,
                         const double *P_re, const double *P_im,
                         const double *Q_re, const double *Q_im,
                         const double *Y0_re, const double *Y0_im,
                         const double *Y1_re, const double *Y1_im,
                         int n, int prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_UCP_SHIM_H */
