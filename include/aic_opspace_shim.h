/* aic_opspace_shim.h — flat-double ccall SHIM for the opspace O2 cb-norm SDP
 * fixture generator (bead aic-pjr, increment O2.5). Mirrors aic_ucp_shim.h: it
 * exposes a FLAT-double ABI so the (serial, one-time) Julia+MOSEK fixture
 * generator (tools/gen_fixtures_opspace_o2.jl) can `ccall` into build/libaic.so
 * WITHOUT any FLINT type (acb_mat_t / slong) crossing the language boundary. Pure
 * marshalling — the math lives in aic_assoc / aic_cstar / aic_opspace_choi.
 *
 * WHY A SHIM (the ccall contract). Same reason as aic_ucp_shim.h: Julia's `ccall`
 * passes int / double / Ptr{Cdouble} cleanly but NOT FLINT structs. The generator
 * needs (i) the th_main_ext isomorphism v: B -> A rebuilt from a UCP channel's
 * flat Kraus, and (ii) the Convention-A adjoint Choi matrices J(v*) and J((v^-1)*)
 * the Watrous diamond-norm SDP consumes (docs/research/opspace_o2_design.md §0.5
 * PINNED CONVENTION). The shim runs the EXACT C cores the C tests exercise
 * (aic_assoc_ecstar_from_phi -> aic_cstar_build -> aic_opspace_choi_v{star,invstar})
 * and writes the two Choi matrices back to flat double arrays the caller owns.
 *
 * WHAT THE GENERATOR DOES WITH THE OUTPUT (the PINNED dims, §0.5). It feeds
 *   J(v*)       (N*n_B square) to diamond_*_rect(J, d_maj = N,   d_min = n_B)  -> ||v||_cb
 *   J((v^-1)*)  (n_B*N square) to diamond_*_rect(J, d_maj = n_B, d_min = N)    -> ||v^-1||_cb
 * with dual tr_sys = 2 (trace the MINOR factor) and primal rho_on = :major. The
 * raw SDP optval IS the cb-norm (normalization factor 1, NO 2/n).
 *
 * KRAUS LAYOUT (Convention A, K-major; identical to aic_ucp_shim.h:18-19):
 *   kraus_re[a*n*n + i*n + j], kraus_im[...]   r ops, each n x n (K_a: H->K, n=n)
 * The Choi convention C[i*n+a,j*n+b] = sum_x conj(K_x[i,a]) K_x[j,b] (the
 * conjugation on the FIRST (i,a) factor) is preserved through aic_ucp_kraus_init.
 *
 * Fail-loud (Rule 4): bad shape (n<1, r<1, null array, n_B>n) aborts at the call
 * site. <=200 LOC (Rule 10).
 */
#ifndef AIC_OPSPACE_SHIM_H
#define AIC_OPSPACE_SHIM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Rebuild the th_main_ext isomorphism v: B -> A from the flat Kraus of an
 * eta-idempotent UCP self-map Phi (r ops, each n x n) and assemble the two
 * Convention-A adjoint Choi matrices the cb-norm SDP consumes:
 *   J(v*)      = aic_opspace_choi_vstar(v)      : (N*n_B) x (N*n_B), v*: M_N -> M_{n_B}.
 *   J((v^-1)*) = aic_opspace_choi_vinvstar(v)   : (n_B*N) x (n_B*N), (v^-1)*: M_{n_B}->M_N.
 * where N = v.n (A's ambient dim) and n_B = v.B->n_B (B's ambient dim).
 *
 * Internally: kraus_from_flat -> aic_assoc_ecstar_from_phi(&ae, phi) ->
 * aic_cstar_build(&B, &v, &iso, &ae.A, eps) -> the two choi assemblers. The two
 * Choi are then extracted to the flat OUT arrays (arf_get_d on the real/imag
 * midpoints, like aic_ucp_choi_diff_d), ROW-MAJOR [p*sz + q] with sz = N*n_B (the
 * two Choi are SAME size N*n_B since they differ only in factor ORDER).
 *
 * OUTPUTS (caller pre-allocates EACH double buffer to length n^4 = (n*n)*(n*n),
 *   the SAFE upper bound: since v is bijective n_B <= N = n, so sz = N*n_B <= n^2
 *   and sz^2 <= n^4 holds for both Choi):
 *   jvs_re,  jvs_im   : J(v*)      entries [p*sz + q], sz = N*n_B.
 *   jvis_re, jvis_im  : J((v^-1)*) entries [p*sz + q], sz = N*n_B (= n_B*N).
 *   *out_N            : N    = v.n        (A's ambient dim; the FORWARD input dim).
 *   *out_nB           : n_B  = v.B->n_B   (B's ambient dim; the INVERSE input dim).
 *   *out_dimB         : dim_B = v.B->dim_B (number of B matrix units = dim of A).
 *   out_iso           : LENGTH-2 caller buffer.
 *                       out_iso[0] = iso_def of v (the O(eps) isomorphism defect
 *                                    midpoint) — record / sanity-check eta=0.
 *                       out_iso[1] = the eps ACTUALLY used (see the eps sentinel
 *                                    below); the generator records it as eps_used.
 *
 * The PINNED dims the generator uses (NOT chosen here — recorded for the caller):
 *   ||v||_cb   : diamond_*_rect(J(v*),      d_maj = *out_N,  d_min = *out_nB).
 *   ||v^-1||_cb: diamond_*_rect(J((v^-1)*), d_maj = *out_nB, d_min = *out_N).
 *
 * `eps` is the eps of A fed to aic_cstar_build:
 *   eps == 0.0 : the eta=0 exact-idempotent oracle path (passed through).
 *   eps  < 0.0 : SENTINEL — DERIVE eps = aic_ecstar_defect_assoc(A) internally,
 *                EXACTLY as test_opspace.c T2 does for the eta>0 fixtures, so the
 *                certifier test (which rebuilds v from the same Kraus with the same
 *                derived eps) gets the IDENTICAL v / Choi. The derived value is
 *                written back to out_iso[1].
 *   eps  > 0.0 : used verbatim (an explicit caller-chosen eps).
 * `prec` is the arb working precision (e.g. 256, the test_opspace CPREC). Returns 0
 * on success; fail-loud asserts on bad shape (n>=1, r>=1, prec>=2, null arrays) and
 * on n_B > n (the caller's n^4 buffer would not hold the (N*n_B)^2 Choi). */
int aic_opspace_choi_shim_d(double *jvs_re, double *jvs_im,
                            double *jvis_re, double *jvis_im,
                            int *out_N, int *out_nB, int *out_dimB,
                            double *out_iso,
                            const double *kraus_re, const double *kraus_im,
                            int n, int r, double eps, int prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_OPSPACE_SHIM_H */
