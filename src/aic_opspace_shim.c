/* aic_opspace_shim.c — flat-double ccall SHIM for the opspace O2 cb-norm SDP
 * fixture generator (bead aic-pjr, increment O2.5). See include/aic_opspace_shim.h
 * for the ABI/layout contract and WHY the shim exists (no FLINT type may cross the
 * Julia ccall boundary). This file is PURE MARSHALLING: flat Kraus in, the EXACT
 * cstar_build + opspace adjoint-Choi cores the C tests exercise, two flat Choi out.
 * No new math.
 *
 * PROVENANCE (Law 1). The pipeline is verbatim the test_opspace.c v-construction
 * (tests/test_opspace.c:136-142):
 *   aic_assoc_ecstar_from_phi(&ae, phi, prec)   (aic_assoc.h:177 — A = Img Phi_tilde)
 *   aic_cstar_build(&B, &v, iso, &ae.A, eps, prec)  (aic_cstar.h:489 — th_main, v: B->A)
 * then the two ADJOINT Choi assemblers (aic_opspace.h:217,226;
 * docs/research/opspace_o2_design.md §0.5 PINNED CONVENTION):
 *   aic_opspace_choi_vstar(Jvs, &v, prec)       -> J(v*),      (N*n_B) x (N*n_B)
 *   aic_opspace_choi_vinvstar(Jvis, &v, prec)   -> J((v^-1)*), (n_B*N) x (n_B*N)
 * The Kraus marshalling mirrors aic_ucp_shim.c:43-54 EXACTLY (Convention A,
 * K-major, conjugation on the FIRST (i,a) factor lives inside aic_ucp_kraus ->
 * the Choi; the shim only sets acb_set_d_d, preserving it).
 *
 * Number path: arb (acb_mat) internally; only the boundary is double. The eventual
 * make-test reads the committed fixtures the generator emits, so steady-state test
 * is Julia-free (the aic_ucp_shim.c / gen_fixtures_certify.jl discipline).
 *
 * Fail-loud (Rule 4): bad shape (n<1, r<1, null array) aborts at the call site;
 * n_B > N aborts (the caller's n^4 buffer would not hold the (N*n_B)^2 Choi).
 * <=200 LOC (Rule 10).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic_assoc.h"
#include "aic_cstar.h"
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_opspace.h"
#include "aic_opspace_shim.h"
#include "aic_ucp.h"

/* Build the internal aic_ucp_kraus self-map Phi (dim_K == dim_H == n) from the
 * flat double Kraus arrays, layout kraus_*[a*n*n + i*n + j] (Convention A,
 * K-major). VERBATIM from aic_ucp_shim.c:43-54 — the conjugation-on-first-factor
 * Choi convention is preserved because we only acb_set_d_d the raw entries. */
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
 * [p*sz + q] = midpoint(J[p,q]). Mirrors aic_ucp_shim.c:86-91 (arf_get_d on the
 * real/imag midpoints; no rounding direction matters — these are feasible-point
 * SEEDS the arb certifier later restores, not a rigorous bracket). */
static void choi_to_flat(double *out_re, double *out_im, const acb_mat_t J, slong sz)
{
    for (slong p = 0; p < sz; p++)
        for (slong q = 0; q < sz; q++) {
            acb_srcptr e = acb_mat_entry(J, p, q);
            out_re[p * sz + q] = arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR);
            out_im[p * sz + q] = arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR);
        }
}

/* aic_opspace_shim.h — rebuild v: B -> A from flat Kraus, assemble J(v*) and
 * J((v^-1)*) (the two adjoint Choi the Watrous SDP consumes), write them to the
 * flat OUT arrays plus (N, n_B, dim_B, iso_def). */
int aic_opspace_choi_shim_d(double *jvs_re, double *jvs_im,
                            double *jvis_re, double *jvis_im,
                            int *out_N, int *out_nB, int *out_dimB,
                            double *out_iso,
                            const double *kraus_re, const double *kraus_im,
                            int n, int r, double eps, int prec)
{
    assert(n >= 1 && "aic_opspace_choi_shim_d: n >= 1 required (Rule 4)");
    assert(r >= 1 && "aic_opspace_choi_shim_d: r >= 1 required (Rule 4)");
    assert(prec >= 2 && "aic_opspace_choi_shim_d: prec >= 2 required (Rule 4)");
    assert(jvs_re && jvs_im && jvis_re && jvis_im &&
           "aic_opspace_choi_shim_d: null Choi out array (Rule 4)");
    assert(out_N && out_nB && out_dimB && out_iso &&
           "aic_opspace_choi_shim_d: null scalar out (Rule 4)");
    assert(kraus_re && kraus_im &&
           "aic_opspace_choi_shim_d: null kraus array (Rule 4)");

    /* ---- (1) Phi -> A = Img Phi_tilde -> v: B -> A (the test_opspace pipeline) -- */
    aic_ucp_kraus phi;
    kraus_from_flat(&phi, kraus_re, kraus_im, n, r);

    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, (slong) prec); /* asserts the prop_P basin */

    /* eps fed to aic_cstar_build. SENTINEL: eps < 0 => DERIVE eps from A exactly as
     * test_opspace.c T2 does (eps = aic_ecstar_defect_assoc(A), the op-norm assoc
     * defect of A) so the certifier test, which rebuilds v with the SAME derived
     * eps, gets the IDENTICAL v / Choi. eps == 0 (the eta=0 oracle path) is passed
     * through unchanged. The eps actually used is written back to out_iso[1]. */
    double eps_used = eps;
    if (eps < 0.0) {
        arb_t ea;
        arb_init(ea);
        aic_ecstar_defect_assoc(ea, &ae.A, (slong) prec);
        eps_used = arf_get_d(arb_midref(ea), ARF_RND_NEAR);
        arb_clear(ea);
    }

    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, eps_used, (slong) prec); /* th_main, v: B->A */

    slong N = v.n, nB = v.B->n_B, dimB = v.B->dim_B;
    /* The caller's OUT buffers are length n^4 = (n*n)^2; both Choi are (N*n_B)^2.
     * Bijective v has dim_A == dim_B and (for these fixtures) n_B <= N = n, so
     * sz = N*n_B <= n^2 and sz^2 <= n^4. Assert it (fail loud, not a silent over-
     * run). N == n always (A <= M_n, ambient preserved). */
    assert(N == (slong) n && "aic_opspace_choi_shim_d: N != n (A ambient drift)");
    assert(nB <= N &&
           "aic_opspace_choi_shim_d: n_B > N — caller's n^4 buffer too small");
    slong sz = N * nB;

    /* ---- (2) the two adjoint Choi (PINNED INPUT-major Convention A, §0.5) ----- */
    acb_mat_t Jvs, Jvis;
    acb_mat_init(Jvs, sz, sz);
    acb_mat_init(Jvis, sz, sz);
    aic_opspace_choi_vstar(Jvs, &v, (slong) prec);     /* J(v*),      dims (N,  n_B) */
    aic_opspace_choi_vinvstar(Jvis, &v, (slong) prec); /* J((v^-1)*), dims (n_B, N)  */

    /* ---- (3) extract to flat row-major OUT arrays + the scalar metadata -------- */
    choi_to_flat(jvs_re, jvs_im, Jvs, sz);
    choi_to_flat(jvis_re, jvis_im, Jvis, sz);
    *out_N    = (int) N;
    *out_nB   = (int) nB;
    *out_dimB = (int) dimB;
    out_iso[0] = arf_get_d(arb_midref(iso), ARF_RND_NEAR); /* iso_def of v       */
    out_iso[1] = eps_used;                                  /* eps actually used  */

    /* ---- (4) free everything (reverse build order; v borrows B and ae.A) ------- */
    acb_mat_clear(Jvis);
    acb_mat_clear(Jvs);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    arb_clear(iso);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return 0;
}
