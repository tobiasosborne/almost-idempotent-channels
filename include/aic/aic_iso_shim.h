/* aic_iso_shim.h — flat-double ccall SHIM for the th_main isomorphism SUMMARY
 * (bead aic-exa.3 [C], shim C3). Mirrors aic_opspace_shim.h: a FLAT-double ABI so
 * the Julia package (AlmostIdempotentChannels.jl, aic-exa.7 [J4]
 * `main_isomorphism`) can `ccall` the th_main isomorphism summary WITHOUT any
 * FLINT type crossing the language boundary. Pure marshalling — the math lives in
 * aic_assoc / aic_cstar / aic_opspace_choi / aic_cbnorm. Design:
 * docs/research/julia_package_design.md §4.3 (C3), CORRECTED by Appendix B (B2).
 *
 * WHAT C3 EMITS (.tex:460 th_main, :484 / FINDINGS §D2 dimension-INDEPENDENT const,
 * :1538-1540 th_main_ext cb-norm). Rebuild the O(eps)-isomorphism v: B -> A
 * (B = (+)_l M_{d_l} a GENUINE C* algebra) from Phi's flat Kraus
 * (aic_assoc_ecstar_from_phi -> aic_cstar_build) and report:
 *   *out_nB     = n_B  = B's block-diagonal representative size (Sum_l d_l);
 *   *out_dimB   = dim_B = Sum_l d_l^2 (= dim A, bijective);
 *   *out_m      = num_blocks (number of C* blocks);
 *   blocks[0..m-1] = B's block sizes d_l (caller pre-allocates length n; m <= n_B <= n);
 *   cbfwd[2]    = [lo,hi] certified bracket on ||v||_cb;
 *   cbinv[2]    = [lo,hi] certified bracket on ||v^{-1}||_cb;
 *   isodef[2]   = [lo,hi] FLOOR/CEIL bracket on the isomorphism defect iso_def of v;
 *   out_eps[0]  = eps_used (the eps sentinel actually fed to aic_cstar_build).
 *
 * THE cb BRACKETS ARE SOLVER-FREE + EIG-FREE (Appendix B2; FINDINGS §C12, §D2).
 * The cb-norm of v is reached via the adjoint duality ||v||_cb = ||v*||_⋄
 * (Watrous 2009). We form the GOLDEN-RULE Convention-A Choi of the ADJOINT maps
 * J(v*) (aic_opspace_choi_vstar: dims d_maj=N, d_min=n_B) and J((v^{-1})*)
 * (aic_opspace_choi_vinvstar: d_maj=n_B, d_min=N) DIRECTLY (NOT by transposing
 * J(v) — the O2 GOLDEN RULE), then call the RECTANGULAR eig-free bracket
 * aic_cbnorm_eigfree_ball_choi_rect (the rect generalization of the self-map
 * aic_cbnorm_eigfree_ball_choi). This is the OPERATOR-norm-faithful diamond bracket
 * (NEVER the Frobenius sigma_min ampliation — that is the §C12 "test that cannot
 * fail"). The TIGHT cb VALUE (the Watrous SDP) is the MOSEK extension's job; here
 * the bracket is solver-free, so the Julia core needs no solver. The bracket is
 * loose by design (hi/lo ~ d_in*d_out) — it CERTIFIES a value, not computes it.
 *
 * KRAUS LAYOUT (Convention A, K-major; identical to aic_ucp_shim.h:18-19):
 *   kraus_re[a*n*n + i*n + j], kraus_im[...]   r ops, each n x n (K_a: H->K, n=n).
 *
 * THE eps SENTINEL (reused VERBATIM from aic_factorize_shim.h:78-89):
 *   eps == 0.0  : eta=0 exact-idempotent oracle  -> eps_used = 0.0.
 *   eps == -2.0 : MULTI-BLOCK eta>0              -> eps_used = eta_proxy (§C11).
 *   eps  < 0.0 (!= -2.0) : DERIVE eps_used = aic_ecstar_defect_assoc(A) midpoint.
 *   eps  > 0.0  : eps_used = eps verbatim.
 * out_eps[0] receives eps_used so a downstream certifier rebuilds the IDENTICAL v.
 *
 * Fail-loud (Rule 4): bad shape (n<1, r<1, prec<2, null array) aborts at the call
 * site; the prop_P basin / non-bijective-v asserts abort deeper. <=200 LOC (Rule 10).
 */
#ifndef AIC_ISO_SHIM_H
#define AIC_ISO_SHIM_H

#ifdef __cplusplus
extern "C" {
#endif

int aic_main_iso_summary_d(int *out_nB, int *out_dimB, int *out_m,
                           int *blocks,                 /* length n */
                           double *cbfwd, double *cbinv, double *isodef, /* len 2 each */
                           double *out_eps,             /* length 1 */
                           const double *kraus_re, const double *kraus_im,
                           int n, int r, double eps, int prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_ISO_SHIM_H */
