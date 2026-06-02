/* aic_assoc_shim.h — flat-double ccall SHIM for the associated eps-C* algebra
 * SUMMARY (bead aic-exa.3 [C], shim C2). Mirrors aic_factorize_shim.h /
 * aic_opspace_shim.h: it exposes a FLAT-double ABI so the Julia package
 * (AlmostIdempotentChannels.jl, aic-exa.7 [J4] `associated_algebra`) can `ccall`
 * the associated-algebra summary WITHOUT any FLINT type (acb_mat_t / slong)
 * crossing the language boundary. Pure marshalling — the math lives in
 * aic_assoc / aic_ecstar. Design: docs/research/julia_package_design.md §4.3 (C2),
 * CORRECTED by Appendix B.
 *
 * WHAT C2 EMITS (.tex:2184-2195, th_almost_idemp). For an eta-idempotent UCP
 * self-map Phi given by flat Kraus, build A = Img Phi_tilde (the eps-C* algebra,
 * Phi_tilde = theta(2 Phi - 1) the exact-idempotent regularization) via
 * aic_assoc_ecstar_from_phi, and report:
 *   *out_n      = n (ambient dim, A <= M_n);
 *   *out_dimA   = dim A;
 *   eps_assoc[0]= lo, eps_assoc[1] = hi : the certified ASSOCIATOR-defect eps of A,
 *                 = aic_ecstar_defect_assoc(A), as an arb ball rounded OUTWARD
 *                 (FLOOR/CEIL, exactly the aic_cbnorm_eigfree_d rounding) so the
 *                 double bracket stays RIGOROUS: lo <= eps <= hi.
 *   basis_re/im : the dim_A Frobenius-ORTHONORMAL basis {B_k} of A, each n x n,
 *                 laid out K-MAJOR [k*n*n + i*n + j] (the Kraus layout). A is an
 *                 OBLIQUE subspace of M_n (FINDINGS §C3); {B_k} is a SNAPSHOT for
 *                 the Julia EpsCStarAlgebra object — its product is the Choi-Effros
 *                 star X*Y = Phi_tilde(XY), NOT plain XY (FINDINGS §C2), so Julia
 *                 must NOT multiply these B_k directly. Caller pre-allocates length
 *                 n^4 (dim_A <= n^2; safe upper bound (n^2)*(n*n) = n^4).
 *
 * KRAUS LAYOUT (Convention A, K-major; identical to aic_ucp_shim.h:18-19):
 *   kraus_re[a*n*n + i*n + j], kraus_im[...]   r ops, each n x n (K_a: H->K, n=n).
 *
 * FAILS LOUD out of the prop_P basin (the aic_assoc_regularize / aic_prop_P assert
 * rho(Phi^2-Phi) >= 1/4 aborts; the Julia caller pre-checks via the eta bracket and
 * surfaces a helpful error naming rho). Fail-loud (Rule 4): bad shape (n<1, r<1,
 * prec<2, null array) aborts at the call site. <=200 LOC (Rule 10).
 */
#ifndef AIC_ASSOC_SHIM_H
#define AIC_ASSOC_SHIM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Build A = Img Phi_tilde from Phi's flat Kraus (aic_assoc_ecstar_from_phi) and
 * report (out_n, out_dimA, eps_assoc[lo,hi] FLOOR/CEIL, basis K-major). `eps` is
 * UNUSED here (A's associator defect is DERIVED from A directly; the eps argument
 * is kept for ABI uniformity with the other shims and ignored). Returns 0 on
 * success; fail-loud asserts on bad shape (n>=1, r>=1, prec>=2, non-null). */
int aic_assoc_summary_d(int *out_n, int *out_dimA,
                        double *eps_assoc,          /* length 2 [lo,hi] */
                        double *basis_re, double *basis_im, /* length n^4 each */
                        const double *kraus_re, const double *kraus_im,
                        int n, int r, double eps, int prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_ASSOC_SHIM_H */
