/* aic_factorize_artifacts_shim.h — flat-double ccall SHIM for the th_factorization
 * ARTIFACTS (bead aic-exa.3 [C], shims C4 + C5). Mirrors aic_factorize_shim.h: a
 * FLAT-double ABI so the Julia package (AlmostIdempotentChannels.jl, aic-exa.7
 * [J4] `factorize`) can `ccall` the full factorization artifacts WITHOUT any FLINT
 * type crossing the language boundary. Pure marshalling — the math lives in
 * aic_assoc / aic_cstar / aic_factorize. Design:
 * docs/research/julia_package_design.md §4.3 (C4/C5), CORRECTED by Appendix B (B1).
 *
 * WHAT C4/C5 REALIZE (.tex:2730-2899 th_factorization, :2159 the dual channels).
 * Rebuild the approximate channel factorization Phi ~ Delta Upsilon through a
 * genuine C* algebra B (the SAME aic_factorize_choi_shim_d build path:
 * aic_assoc_ecstar_from_phi -> aic_cstar_build -> aic_factorize_build/_delta_build/
 * _upsilon_build), then emit the CHANNEL (CPTP, state-space) Kraus + the certified
 * eig-free round-trip brackets.
 *
 * TWO-CALL PROTOCOL (the Dec/Enc Kraus counts are DATA-DEPENDENT, design §4.3).
 *   C4 = aic_factorize_artifacts_sizes_d : rebuild, report SIZES only (N, n_B,
 *        dim_B, num_blocks, block sizes, AND *out_rDec / *out_rEnc — the number of
 *        Dec / Enc Kraus operators, which the Choi->Kraus PSD-eigen extraction
 *        produces data-dependently). The Julia caller sizes its buffers from these,
 *        then calls C5.
 *   C5 = aic_factorize_artifacts_d : rebuild, emit the full artifacts; ASSERTS the
 *        rebuilt rDec/rEnc match the caller-passed counts (fail loud if the caller
 *        mis-sized via C4).
 *
 * THE CHANNEL DIRECTION (Appendix B1, the §C13 dual-direction trap — LOAD-BEARING).
 * Ground truth aic_factorize.h:306-316 + .tex:2159:
 *   Delta ("encode" OBSERVABLE map B -> B(H)), Upsilon ("decode" OBSERVABLE map
 *   B(H) -> B); the STATE channels are the duals:
 *     Dec = Delta*   = aic_factorize_dec_kraus : channel B(H) -> B, "produces a
 *                      state on B"; Kraus n_B x N; layout [a*(nB*N) + i*N + j].
 *     Enc = Upsilon* = aic_factorize_enc_kraus : channel B -> B(H); Kraus N x n_B;
 *                      layout [a*(N*nB) + i*nB + j].
 *   So C5 MUST emit dec_* from dec_kraus and enc_* from enc_kraus (DO NOT cross
 *   them): swapping them breaks the decode-then-encode round-trip / dim assertion
 *   (the [J4] mutation tooth). FINDINGS §C13 (F-ancilla ordering), §C14 (Delta' only
 *   O(eta)-CP -> the Choi->Kraus PSD-cone projection is internal, admits the O(eta^2)
 *   cone defect, aborts a genuine non-CP).
 *
 * THE ROUND-TRIP BRACKETS (.tex:2733 DelUps, :2739 UpsDel). delups/upsdel are the
 * FLOOR/CEIL doubles of aic_factorize_eigfree_delups / _upsdel (the eig-free
 * certified brackets on ||Delta Upsilon - Phi||_cb / ||Upsilon Delta - 1_B||_cb;
 * the cores already exist — C5 only rounds them OUTWARD). The TIGHT cb VALUE is the
 * MOSEK extension's job (aic_factorize_choi_shim_d + the SDP); kept solver-free here.
 *
 * out_eta[2] = [eta_proxy, eps_used] (identical to aic_factorize_choi_shim_d).
 *
 * KRAUS LAYOUT (input Phi, Convention A, K-major; identical to aic_ucp_shim.h:18-19):
 *   kraus_re[a*n*n + i*n + j], kraus_im[...]   r ops, each n x n (K_a: H->K, n=n).
 *
 * THE eps SENTINEL (reused VERBATIM from aic_factorize_shim.h:78-89): eps==0.0 ->
 * eta=0 oracle (eps_used=0); eps==-2.0 -> multi-block eta>0 (eps_used = eta_proxy,
 * §C11); eps<0 (!=-2) -> derive aic_ecstar_defect_assoc; eps>0 verbatim. out_eta[1]
 * receives eps_used so a downstream certifier rebuilds the IDENTICAL artifacts.
 *
 * Fail-loud (Rule 4): bad shape aborts at the call site; the prop_P basin and the
 * factorize-build asserts abort deeper. <=200 LOC per .c (Rule 10; C4/C5 share one
 * build helper in src/aic_factorize_artifacts_shim.c).
 */
#ifndef AIC_FACTORIZE_ARTIFACTS_SHIM_H
#define AIC_FACTORIZE_ARTIFACTS_SHIM_H

#ifdef __cplusplus
extern "C" {
#endif

/* C4 — rebuild the factorize pipeline and report ONLY the SIZES needed to allocate
 * the C5 buffers. blocks[] caller-allocated length n (m <= n_B <= n). *out_rDec /
 * *out_rEnc are the data-dependent Dec / Enc Kraus counts. Returns 0; fail-loud
 * asserts on bad shape (n>=1, r>=1, prec>=2, non-null). */
int aic_factorize_artifacts_sizes_d(int *out_N, int *out_nB, int *out_dimB,
                                    int *out_m, int *blocks,    /* blocks length n */
                                    int *out_rDec, int *out_rEnc,
                                    const double *kraus_re, const double *kraus_im,
                                    int n, int r, double eps, int prec);

/* C5 — rebuild the factorize pipeline and emit the full artifacts. dec_* / enc_*
 * are caller-allocated to rDec*(n_B*N) / rEnc*(N*n_B) (from C4). delups, upsdel,
 * out_eta are caller-allocated. ASSERTS the rebuilt Dec/Enc counts match rDec/rEnc
 * (fail loud if the caller mis-sized). Returns 0; fail-loud asserts on bad shape. */
int aic_factorize_artifacts_d(double *dec_re, double *dec_im, int rDec,
                              double *enc_re, double *enc_im, int rEnc,
                              double *delups, double *upsdel,  /* len 2 each */
                              double *out_eta,                  /* len 2 */
                              const double *kraus_re, const double *kraus_im,
                              int n, int r, double eps, int prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_FACTORIZE_ARTIFACTS_SHIM_H */
