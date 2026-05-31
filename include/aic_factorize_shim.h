/* aic_factorize_shim.h — flat-double ccall SHIM for the factorize F4.2 tight
 * cb-norm SDP fixture generator (bead aic-tff, increment F4.2). Mirrors
 * aic_opspace_shim.h: it exposes a FLAT-double ABI so the (serial, one-time)
 * Julia+MOSEK fixture generator (tools/gen_fixtures_factorize_f4.jl) can `ccall`
 * into build/libaic.so WITHOUT any FLINT type (acb_mat_t / slong) crossing the
 * language boundary. Pure marshalling — the math lives in aic_assoc / aic_cstar /
 * aic_factorize. Design: docs/research/factorize_f4_design.md §C.2.
 *
 * WHY A SHIM (the ccall contract). Same reason as aic_opspace_shim.h: Julia's
 * `ccall` passes int / double / Ptr{Cdouble} cleanly but NOT FLINT structs. The
 * generator needs (i) the full factorize pipeline rebuilt from a UCP channel's
 * flat Kraus, and (ii) the two Convention-A Choi DIFFERENCE matrices the tight
 * Watrous diamond-norm SDP consumes (the th_factorization headline bounds,
 * approximate_algebras.tex:2730-2739). The shim runs the EXACT C cores the C test
 * exercises (aic_assoc_ecstar_from_phi -> aic_cstar_build -> aic_factorize_build
 * -> aic_factorize_delta_build -> aic_factorize_upsilon_build), assembles the two
 * Choi via aic_factorize_choi_delups / aic_factorize_choi_upsdel, and writes them
 * back to flat double arrays the caller owns.
 *
 * THE TWO OUTPUT CHOI (th_factorization, .tex:2733/2739):
 *   J_DelUps = Choi(Delta Upsilon) - Choi(Phi)  : N^2 x N^2  (DelUps, .tex:2733).
 *   J_UpsDel = Choi(Upsilon Delta - 1_B)        : n_B^2 x n_B^2 (UpsDel2, .tex:2739).
 * The generator feeds each to the tight cb-norm SDP; the raw optval is ||.||_cb.
 *
 * KRAUS LAYOUT (Convention A, K-major; identical to aic_opspace_shim.h:23-26):
 *   kraus_re[a*n*n + i*n + j], kraus_im[...]   r ops, each n x n (K_a: H->K, n=n).
 * The Choi convention (conjugation on the FIRST factor) is preserved through
 * aic_ucp_kraus_init — the shim only acb_set_d_d the raw entries.
 *
 * BUFFER SIZING (§C.2 buffer note). J_DelUps is N^2 x N^2 (length N^4); J_UpsDel
 * is n_B^2 x n_B^2 (length n_B^4). Since n_B <= N, allocate BOTH at N^4 = (n*n)^2;
 * the shim asserts the shapes (Rule 4).
 *
 * Fail-loud (Rule 4): bad shape (n<1, r<1, prec<1, null array) aborts at the call
 * site. The factorize build routines fail-loud via assert on out-of-basin input
 * (eta >= 1/4, degenerate spectrum); the generator pre-filters to in-basin
 * channels, so those asserts are not caught here. <=200 LOC (Rule 10).
 */
#ifndef AIC_FACTORIZE_SHIM_H
#define AIC_FACTORIZE_SHIM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Rebuild the full factorize pipeline from the flat Kraus of an eta-idempotent
 * UCP self-map Phi (r ops, each n x n) and assemble the two Convention-A Choi
 * difference matrices the tight cb-norm SDP consumes:
 *   J_DelUps = Choi(Delta Upsilon) - Choi(Phi)  : N^2 x N^2,   row-major [p*sz_du+q],
 *              sz_du = N*N.
 *   J_UpsDel = Choi(Upsilon Delta - 1_B)        : n_B^2 x n_B^2, [p*sz_ud+q],
 *              sz_ud = n_B*n_B.
 * where N = F.N (= n, A's ambient dim) and n_B = F.n_B (B's block-diagonal rep size).
 *
 * Internally: kraus_from_flat -> aic_assoc_ecstar_from_phi(&ae, phi) ->
 * aic_cstar_build(&B, &v, iso, &ae.A, eps_used) -> aic_factorize_build /
 * _delta_build / _upsilon_build -> aic_factorize_choi_delups / _choi_upsdel.
 * The two Choi are extracted to the flat OUT arrays (arf_get_d on the real/imag
 * midpoints, like aic_opspace_shim's choi_to_flat — feasible-point SEEDS the arb
 * certifier restores, not a rigorous bracket).
 *
 * OUTPUTS (caller pre-allocates EACH double buffer to length n^4 = (n*n)*(n*n),
 *   the SAFE upper bound: since n_B <= N = n, sz_du = N*N = n^2 (sz_du^2 = n^4)
 *   and sz_ud = n_B*n_B <= n^2 (sz_ud^2 <= n^4) both fit):
 *   jdu_re, jdu_im : J_DelUps entries [p*sz_du + q], sz_du = N*N.
 *   jud_re, jud_im : J_UpsDel entries [p*sz_ud + q], sz_ud = n_B*n_B.
 *   *out_N         : N    = F.N      (A's ambient dim = n; the DelUps Choi side).
 *   *out_nB        : n_B  = F.n_B    (B's block-diagonal rep size; the UpsDel side).
 *   *out_dimB      : dim_A = F.dim_A (number of B matrix units = dim of A).
 *   out_eta        : LENGTH-2 caller buffer.
 *                    out_eta[0] = eta_proxy = ||S_Phi^2 - S_Phi||_op (op-norm
 *                                 midpoint route; the §C5 eta proxy, recorded /
 *                                 sanity-checked ~0 for the eta=0 oracle).
 *                    out_eta[1] = the eps ACTUALLY used (eps_used, see the eps
 *                                 sentinel below); the generator records it and the
 *                                 C certifier test rebuilds with this positive value.
 *
 * `eps` is the eps of A fed to aic_cstar_build, via this SENTINEL:
 *   eps == 0.0  : the eta=0 exact-idempotent oracle path -> eps_used = 0.0
 *                 (passed through to aic_cstar_build).
 *   eps == -2.0 : MULTI-BLOCK eta>0 sentinel -> eps_used = eta_proxy (out_eta[0]).
 *                 The FINDINGS §C11/§C13(c) caveat: pass eta, NOT the ~700x smaller
 *                 aic_ecstar_defect_assoc, else the Stage-1 errreduce C0 gate fires
 *                 (10*eps < the true O(eta) inclusion defect) and the build aborts.
 *   eps  < 0.0 (and != -2.0) : DERIVE eps_used = aic_ecstar_defect_assoc(A) midpoint
 *                 (the opspace-shim "derive" sentinel) so a single-block certifier
 *                 test rebuilding the same pipeline gets the IDENTICAL J.
 *   eps  > 0.0  : eps_used = eps verbatim (an explicit caller-chosen eps).
 * The eps_used value is written to out_eta[1].
 *
 * `prec` is the arb working precision (e.g. 256, the test_factorize CPREC). Returns
 * 0 on success; fail-loud asserts on bad shape (n>=1, r>=1, prec>=1, null arrays).
 * (Out-of-basin Kraus aborts deeper in the build routines — the generator
 * pre-filters, so those asserts are not caught here.) */
int aic_factorize_choi_shim_d(double *jdu_re, double *jdu_im,
                              double *jud_re, double *jud_im,
                              int *out_N, int *out_nB, int *out_dimB,
                              double *out_eta,
                              const double *kraus_re, const double *kraus_im,
                              int n, int r, double eps, int prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_FACTORIZE_SHIM_H */
