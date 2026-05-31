/* aic_factorize_shim.c — flat-double ccall SHIM for the factorize F4.2 tight
 * cb-norm SDP fixture generator (bead aic-tff, increment F4.2). See
 * include/aic_factorize_shim.h for the ABI/layout contract and WHY the shim exists
 * (no FLINT type may cross the Julia ccall boundary). This file is PURE
 * MARSHALLING: flat Kraus in, the EXACT factorize pipeline the C tests exercise,
 * two flat Convention-A Choi DIFFERENCE matrices out. No new math. It is a close
 * MIRROR of src/aic_opspace_shim.c (same kraus_from_flat / choi_to_flat helpers,
 * same eps-sentinel discipline). Design: docs/research/factorize_f4_design.md §C.2.
 *
 * PROVENANCE (Law 1). The pipeline is verbatim the test_factorize.c construction
 * (fixture channel -> aic_assoc_ecstar_from_phi -> aic_cstar_build -> (B,v) ->
 * aic_factorize_build -> _delta_build -> _upsilon_build) plus the two F4.1 Choi
 * assemblers (aic_factorize.h:286,293):
 *   aic_factorize_choi_delups(J,  &F)  -> J_DelUps = Choi(Delta Upsilon)-Choi(Phi),
 *                                         N^2 x N^2 (DelUps, .tex:2733).
 *   aic_factorize_choi_upsdel(J2, &F)  -> J_UpsDel = Choi(Upsilon Delta - 1_B),
 *                                         n_B^2 x n_B^2 (UpsDel2, .tex:2739).
 * The eta_proxy = ||S_Phi^2 - S_Phi||_op (op-norm midpoint route) replicates
 * test_factorize.c:58-79 (the §C5 proxy); it is BOTH a recorded diagnostic and the
 * eps_used value for the eps == -2.0 multi-block sentinel.
 *
 * THE eps SENTINEL (FINDINGS §C11/§C13(c)). eps == 0.0 -> eps_used = 0.0 (eta=0
 * oracle); eps == -2.0 -> eps_used = eta_proxy (the multi-block eta>0 mode: pass
 * eta, NOT the ~700x smaller aic_ecstar_defect_assoc, else the Stage-1 errreduce
 * C0 gate fires and aic_cstar_build aborts); eps < 0 (!= -2.0) -> eps_used =
 * aic_ecstar_defect_assoc(A) midpoint (single-block "derive"); eps > 0 verbatim.
 * eps_used is written to out_eta[1] (the generator emits it; the C test rebuilds
 * with this positive value).
 *
 * Number path: arb (acb_mat) internally; only the boundary is double. The eventual
 * make-test reads the committed fixtures the generator emits, so steady-state test
 * is Julia-free (the gen_fixtures discipline). Fail-loud (Rule 4): bad shape aborts
 * at the call site; the factorize build routines fail-loud on out-of-basin input
 * (the generator pre-filters, so those are not caught here). <=200 LOC (Rule 10).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_assoc.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_factorize.h"
#include "aic/aic_factorize_shim.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

/* Build the internal aic_ucp_kraus self-map Phi (dim_K == dim_H == n) from the
 * flat double Kraus arrays, layout kraus_*[a*n*n + i*n + j] (Convention A,
 * K-major). VERBATIM from aic_opspace_shim.c:47-58 — the conjugation-on-first-
 * factor Choi convention is preserved because we only acb_set_d_d the raw entries. */
static void kraus_from_flat(aic_ucp_kraus *phi, const double *kraus_re,
                            const double *kraus_im, int n, int r)
{
    aic_ucp_kraus_init(phi, n, n, r);
    for (int a = 0; a < r; a++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                int off = a * n * n + i * n + j;
                acb_set_d_d(acb_mat_entry(phi->K[a], i, j),
                            kraus_re[off], kraus_im[off]);
            }
}

/* Extract the (sz x sz) acb_mat J to flat ROW-MAJOR double arrays out_re/out_im,
 * [p*sz + q] = midpoint(J[p,q]). VERBATIM from aic_opspace_shim.c:64-72 (arf_get_d
 * on the real/imag midpoints; these are feasible-point SEEDS the arb certifier
 * later restores, not a rigorous bracket). */
static void choi_to_flat(double *out_re, double *out_im, const acb_mat_t J, slong sz)
{
    for (slong p = 0; p < sz; p++)
        for (slong q = 0; q < sz; q++) {
            acb_srcptr e = acb_mat_entry(J, p, q);
            out_re[p * sz + q] = arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR);
            out_im[p * sz + q] = arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR);
        }
}

/* eta_proxy = ||S_Phi^2 - S_Phi||_op via the §C5 midpoint route (replicates
 * test_factorize.c:58-79). Build S = aic_assoc_superop_from_ucp(phi), S2 = S^2,
 * D = S2 - S, take the entrywise midpoint into M, op-norm of M, midpoint double. */
static double eta_proxy(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S, S2, D, M;
    acb_mat_init(S, nn, nn);
    acb_mat_init(S2, nn, nn);
    acb_mat_init(D, nn, nn);
    acb_mat_init(M, nn, nn);
    aic_assoc_superop_from_ucp(S, phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    for (slong i = 0; i < nn; i++)
        for (slong j = 0; j < nn; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    double o = arf_get_d(arb_midref(e), ARF_RND_NEAR);
    arb_clear(e);
    acb_mat_clear(M); acb_mat_clear(D); acb_mat_clear(S2); acb_mat_clear(S);
    return o;
}

/* aic_factorize_shim.h — rebuild the full factorize pipeline from flat Kraus,
 * assemble J_DelUps and J_UpsDel (the two Convention-A Choi difference matrices the
 * tight cb-norm SDP consumes), write them to the flat OUT arrays plus (N, n_B,
 * dim_B, eta_proxy, eps_used). */
int aic_factorize_choi_shim_d(double *jdu_re, double *jdu_im,
                              double *jud_re, double *jud_im,
                              int *out_N, int *out_nB, int *out_dimB,
                              double *out_eta,
                              const double *kraus_re, const double *kraus_im,
                              int n, int r, double eps, int prec)
{
    assert(n >= 1 && "aic_factorize_choi_shim_d: n >= 1 required (Rule 4)");
    assert(r >= 1 && "aic_factorize_choi_shim_d: r >= 1 required (Rule 4)");
    assert(prec >= 1 && "aic_factorize_choi_shim_d: prec >= 1 required (Rule 4)");
    assert(jdu_re && jdu_im && jud_re && jud_im &&
           "aic_factorize_choi_shim_d: null Choi out array (Rule 4)");
    assert(out_N && out_nB && out_dimB && out_eta &&
           "aic_factorize_choi_shim_d: null scalar out (Rule 4)");
    assert(kraus_re && kraus_im &&
           "aic_factorize_choi_shim_d: null kraus array (Rule 4)");

    /* ---- (1) Phi -> A = Img Phi_tilde, plus the recorded eta proxy ------------ */
    aic_ucp_kraus phi;
    kraus_from_flat(&phi, kraus_re, kraus_im, n, r);

    out_eta[0] = eta_proxy(&phi, (slong) prec); /* ||S_Phi^2 - S_Phi||_op */

    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, (slong) prec); /* asserts the prop_P basin */

    /* ---- (2) resolve eps_used per the SENTINEL (header / FINDINGS §C11) -------- */
    double eps_used;
    if (eps == 0.0) {
        eps_used = 0.0;                 /* eta=0 exact-idempotent oracle */
    } else if (eps == -2.0) {
        eps_used = out_eta[0];          /* MULTI-BLOCK eta>0: use eta, NOT eps_assoc */
    } else if (eps < 0.0) {
        arb_t ea;                       /* single-block "derive" (eps_assoc midpoint) */
        arb_init(ea);
        aic_ecstar_defect_assoc(ea, &ae.A, (slong) prec);
        eps_used = arf_get_d(arb_midref(ea), ARF_RND_NEAR);
        arb_clear(ea);
    } else {
        eps_used = eps;                 /* explicit caller-chosen eps */
    }
    out_eta[1] = eps_used;

    /* ---- (3) v: B -> A, then the full factorize handle (Delta + Upsilon) ------ */
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, eps_used, (slong) prec); /* th_main, v: B->A */

    aic_factorize F;
    aic_factorize_build(&F, &v, &ae, &phi, (slong) prec);
    aic_factorize_delta_build(&F, (slong) prec);          /* F2 encode  Delta   */
    aic_factorize_upsilon_build(&F, 0.3, (slong) prec);   /* F3 decode  Upsilon, tol_sigma=0.3 */

    *out_N    = (int) F.N;
    *out_nB   = (int) F.n_B;
    *out_dimB = (int) F.dim_A;

    /* ---- (4) assemble the two Choi differences, extract to flat row-major ----- */
    acb_mat_t J;
    acb_mat_init(J, F.N * F.N, F.N * F.N);
    aic_factorize_choi_delups(J, &F, (slong) prec);       /* N^2 x N^2 (.tex:2733) */
    choi_to_flat(jdu_re, jdu_im, J, F.N * F.N);
    acb_mat_clear(J);

    acb_mat_t J2;
    acb_mat_init(J2, F.n_B * F.n_B, F.n_B * F.n_B);
    aic_factorize_choi_upsdel(J2, &F, (slong) prec);      /* n_B^2 x n_B^2 (.tex:2739) */
    choi_to_flat(jud_re, jud_im, J2, F.n_B * F.n_B);
    acb_mat_clear(J2);

    /* ---- (5) free in reverse-dependency order (v borrows B and ae.A) ---------- */
    aic_factorize_clear(&F);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return 0;
}
