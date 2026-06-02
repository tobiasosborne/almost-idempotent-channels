/* aic_iso_shim.c — flat-double ccall SHIM for the th_main isomorphism SUMMARY
 * (bead aic-exa.3 [C], shim C3). See include/aic/aic_iso_shim.h for the ABI/layout
 * contract and WHY the shim exists. PURE MARSHALLING: flat Kraus in, the EXACT
 * cstar_build + opspace adjoint-Choi + rect eig-free cb-bracket cores the C tests
 * exercise, the isomorphism summary out. No new math (the rect bracket core
 * aic_cbnorm_eigfree_ball_choi_rect lives in aic_cbnorm.c). Close MIRROR of
 * aic_opspace_shim.c. Design: docs/research/julia_package_design.md §4.3 (C3),
 * Appendix B2.
 *
 * PROVENANCE (Law 1). The pipeline is verbatim the test_opspace v-construction:
 *   aic_assoc_ecstar_from_phi(&ae, phi, prec)     (A = Img Phi_tilde)
 *   aic_cstar_build(&B, &v, iso, &ae.A, eps_used, prec)  (th_main, v: B -> A)
 * The cb brackets use the adjoint Choi (aic_opspace.h GOLDEN RULE §0.5):
 *   J(v*)      = aic_opspace_choi_vstar(v)      : (N*n_B)^2, dims (d_maj=N,   d_min=n_B)
 *   J((v^-1)*) = aic_opspace_choi_vinvstar(v)   : (n_B*N)^2, dims (d_maj=n_B, d_min=N)
 * fed to aic_cbnorm_eigfree_ball_choi_rect -> ||v||_cb = ||v*||_⋄, ||v^-1||_cb =
 * ||(v^-1)*||_⋄ (Watrous 2009). These adjoint Choi are built DIRECTLY (NOT by
 * transposing J(v) — the §C12.O2 / §0.5 GOLDEN RULE). The cb-norm is the
 * operator-faithful diamond, NEVER the Frobenius sigma_min ampliation (FINDINGS §C12).
 *
 * THE eps SENTINEL (FINDINGS §C11/§C13(c); identical to aic_factorize_shim.c:22-28).
 * eta_proxy = ||S_Phi^2 - S_Phi||_op (op-norm midpoint route; the §C5 proxy) is the
 * eps_used value for the eps == -2.0 multi-block sentinel.
 *
 * Number path: arb (acb_mat) internally; only the boundary is double. The cb /
 * isodef brackets are rounded OUTWARD (FLOOR/CEIL) so they stay RIGOROUS after the
 * cast. Fail-loud (Rule 4): bad shape aborts at the call site. <=200 LOC (Rule 10).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_assoc.h"
#include "aic/aic_cbnorm.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_iso_shim.h"
#include "aic/aic_mat.h"
#include "aic/aic_opspace.h"
#include "aic/aic_ucp.h"

/* Build Phi from flat Kraus (Convention A, K-major). VERBATIM aic_opspace_shim.c. */
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

/* eta_proxy = ||S_Phi^2 - S_Phi||_op (the §C5 midpoint route; replicates
 * aic_factorize_shim.c:86-107). Used only for the eps == -2.0 sentinel. */
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

/* Two-sided arb bracket -> [lo,hi] doubles: lo FLOOR(arb_lower(loB)), hi
 * CEIL(arb_upper(hiB)), rounded OUTWARD so the double bracket stays RIGOROUS
 * (identical to aic_ucp_shim.c's aic_cbnorm_eigfree_d rounding). The rect cb bracket
 * comes as two separate arb balls; a single ball b is the loB==hiB==b case. */
static void twoball_to_bracket(double *out2, const arb_t loB, const arb_t hiB,
                               slong prec)
{
    arf_t t;
    arf_init(t);
    arb_get_lbound_arf(t, loB, prec);
    out2[0] = arf_get_d(t, ARF_RND_FLOOR);
    arb_get_ubound_arf(t, hiB, prec);
    out2[1] = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
}

/* aic_iso_shim.h — rebuild v: B -> A, report block sizes, iso_def, and the
 * solver-free eig-free cb brackets ||v||_cb / ||v^-1||_cb. */
int aic_main_iso_summary_d(int *out_nB, int *out_dimB, int *out_m,
                           int *blocks,
                           double *cbfwd, double *cbinv, double *isodef,
                           double *out_eps,
                           const double *kraus_re, const double *kraus_im,
                           int n, int r, double eps, int prec)
{
    assert(n >= 1 && "aic_main_iso_summary_d: n >= 1 required (Rule 4)");
    assert(r >= 1 && "aic_main_iso_summary_d: r >= 1 required (Rule 4)");
    assert(prec >= 2 && "aic_main_iso_summary_d: prec >= 2 required (Rule 4)");
    assert(out_nB && out_dimB && out_m && blocks && out_eps &&
           "aic_main_iso_summary_d: null scalar/blocks out (Rule 4)");
    assert(cbfwd && cbinv && isodef &&
           "aic_main_iso_summary_d: null bracket out (Rule 4)");
    assert(kraus_re && kraus_im &&
           "aic_main_iso_summary_d: null kraus array (Rule 4)");

    /* ---- (1) Phi -> A -> resolve eps_used per the SENTINEL --------------------- */
    aic_ucp_kraus phi;
    kraus_from_flat(&phi, kraus_re, kraus_im, n, r);

    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, (slong) prec); /* asserts the prop_P basin */

    double eps_used;
    if (eps == 0.0) {
        eps_used = 0.0;
    } else if (eps == -2.0) {
        eps_used = eta_proxy(&phi, (slong) prec);
    } else if (eps < 0.0) {
        arb_t ea;
        arb_init(ea);
        aic_ecstar_defect_assoc(ea, &ae.A, (slong) prec);
        eps_used = arf_get_d(arb_midref(ea), ARF_RND_NEAR);
        arb_clear(ea);
    } else {
        eps_used = eps;
    }
    out_eps[0] = eps_used;

    /* ---- (2) v: B -> A (th_main) ---------------------------------------------- */
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, eps_used, (slong) prec);

    slong N = v.n, nB = v.B->n_B, m = v.B->num_blocks, dimB = v.B->dim_B;
    *out_nB   = (int) nB;
    *out_dimB = (int) dimB;
    *out_m    = (int) m;
    assert(m <= (slong) n && "aic_main_iso_summary_d: m > n (blocks[] overflow)");
    for (slong l = 0; l < m; l++) blocks[l] = (int) v.B->d[l];

    /* iso_def: it is a (possibly zero-radius) arb ball; round OUTWARD. */
    twoball_to_bracket(isodef, iso, iso, (slong) prec); /* iso single ball */

    /* ---- (3) solver-free eig-free cb brackets via the adjoint Choi (B2) -------- */
    /* ||v||_cb = ||v*||_⋄: J(v*) is (N*n_B)^2, map dims d_maj=N (in), d_min=n_B (out). */
    {
        slong sz = N * nB;
        acb_mat_t Jvs;
        acb_mat_init(Jvs, sz, sz);
        aic_opspace_choi_vstar(Jvs, &v, (slong) prec);
        arb_t loB, hiB;
        arb_init(loB);
        arb_init(hiB);
        aic_cbnorm_eigfree_ball_choi_rect(loB, hiB, Jvs, N, nB, (slong) prec);
        twoball_to_bracket(cbfwd, loB, hiB, (slong) prec);
        arb_clear(hiB);
        arb_clear(loB);
        acb_mat_clear(Jvs);
    }
    /* ||v^-1||_cb = ||(v^-1)*||_⋄: J((v^-1)*) is (n_B*N)^2, dims d_maj=n_B, d_min=N. */
    {
        slong sz = nB * N;
        acb_mat_t Jvis;
        acb_mat_init(Jvis, sz, sz);
        aic_opspace_choi_vinvstar(Jvis, &v, (slong) prec);
        arb_t loB, hiB;
        arb_init(loB);
        arb_init(hiB);
        aic_cbnorm_eigfree_ball_choi_rect(loB, hiB, Jvis, nB, N, (slong) prec);
        twoball_to_bracket(cbinv, loB, hiB, (slong) prec);
        arb_clear(hiB);
        arb_clear(loB);
        acb_mat_clear(Jvis);
    }

    /* ---- (4) free (reverse build order; v borrows B and ae.A) ----------------- */
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    arb_clear(iso);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return 0;
}
