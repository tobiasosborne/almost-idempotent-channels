/* aic_factorize_artifacts_shim.c — flat-double ccall SHIM for the th_factorization
 * ARTIFACTS, shim C4 (aic_factorize_artifacts_sizes_d) + the SHARED C4/C5 pipeline
 * build helpers (aic_factorize_artifacts_internal.h). Shim C5 lives in
 * aic_factorize_artifacts2_shim.c (Rule 10 split). See
 * include/aic/aic_factorize_artifacts_shim.h for the ABI/layout contract, the
 * Appendix-B1 channel-direction binding (Dec=Delta*, Enc=Upsilon*), and WHY the
 * shim exists. PURE MARSHALLING over the SAME pipeline aic_factorize_choi_shim_d
 * rebuilds (do NOT re-derive — design §4.2): no new math. Close MIRROR of
 * aic_factorize_shim.c (same kraus_from_flat / eta_proxy / eps-sentinel discipline).
 *
 * PROVENANCE (Law 1). The build is verbatim aic_factorize_shim.c:130-167:
 *   kraus_from_flat -> aic_assoc_ecstar_from_phi -> (eps sentinel) ->
 *   aic_cstar_build -> aic_factorize_build -> _delta_build -> _upsilon_build(0.3).
 * C4 reports SIZES only, including the DATA-DEPENDENT Dec/Enc Kraus counts (the
 * Choi->Kraus PSD-eigen extraction in aic_factorize_dec_kraus/_enc_kraus produces a
 * data-dependent r). The eps SENTINEL is FINDINGS §C11/§C13(c) (identical to
 * aic_factorize_shim.c:22-28); eta_proxy = ||S_Phi^2 - S_Phi||_op (the §C5 route).
 *
 * Number path: arb (acb_mat) internally; only the boundary is double. Fail-loud
 * (Rule 4): bad shape aborts at the call site; the prop_P basin / factorize-build
 * asserts abort deeper. <=200 LOC (Rule 10).
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
#include "aic/aic_factorize_artifacts_shim.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_factorize_artifacts_internal.h"

/* Build Phi from flat Kraus (Convention A, K-major). VERBATIM aic_factorize_shim.c. */
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

/* eta_proxy = ||S_Phi^2 - S_Phi||_op (the §C5 route; VERBATIM aic_factorize_shim.c:86). */
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

/* Resolve eps_used per the SENTINEL (aic_factorize_shim.c:139-153). */
static double resolve_eps(double eps, double eta_pr, const aic_assoc_ecstar *ae,
                          slong prec)
{
    if (eps == 0.0) return 0.0;
    if (eps == -2.0) return eta_pr;
    if (eps < 0.0) {
        arb_t ea;
        arb_init(ea);
        aic_ecstar_defect_assoc(ea, &ae->A, prec);
        double v = arf_get_d(arb_midref(ea), ARF_RND_NEAR);
        arb_clear(ea);
        return v;
    }
    return eps;
}

/* aic_factorize_artifacts_internal.h — the SHARED C4/C5 pipeline build (+ free). */
void aic_fart_build_all(aic_factorize *F, aic_ucp_kraus *phi, aic_assoc_ecstar *ae,
                        aic_dhom_B *B, aic_dhom_v *v, arb_t iso, double *out_eta,
                        const double *kraus_re, const double *kraus_im,
                        int n, int r, double eps, slong prec)
{
    kraus_from_flat(phi, kraus_re, kraus_im, n, r);
    out_eta[0] = eta_proxy(phi, prec);
    aic_assoc_ecstar_from_phi(ae, phi, prec); /* asserts the prop_P basin */
    double eps_used = resolve_eps(eps, out_eta[0], ae, prec);
    out_eta[1] = eps_used;
    arb_init(iso);
    aic_cstar_build(B, v, iso, &ae->A, eps_used, prec);   /* th_main, v: B -> A */
    aic_factorize_build(F, v, ae, phi, prec);
    aic_factorize_delta_build(F, prec);                   /* F2 encode  Delta   */
    aic_factorize_upsilon_build(F, 0.3, prec);            /* F3 decode  Upsilon */
}

void aic_fart_build_free(aic_factorize *F, aic_ucp_kraus *phi, aic_assoc_ecstar *ae,
                         aic_dhom_B *B, aic_dhom_v *v, arb_t iso)
{
    aic_factorize_clear(F);
    arb_clear(iso);
    aic_dhom_v_clear(v);
    aic_dhom_B_clear(B);
    aic_assoc_ecstar_clear(ae);
    aic_ucp_kraus_clear(phi);
}

/* C4 — sizes only (incl. the data-dependent Dec/Enc Kraus counts rDec/rEnc). */
int aic_factorize_artifacts_sizes_d(int *out_N, int *out_nB, int *out_dimB,
                                    int *out_m, int *blocks,
                                    int *out_rDec, int *out_rEnc,
                                    const double *kraus_re, const double *kraus_im,
                                    int n, int r, double eps, int prec)
{
    assert(n >= 1 && r >= 1 && prec >= 2 &&
           "aic_factorize_artifacts_sizes_d: bad shape (Rule 4)");
    assert(out_N && out_nB && out_dimB && out_m && blocks && out_rDec && out_rEnc &&
           kraus_re && kraus_im &&
           "aic_factorize_artifacts_sizes_d: null out/in (Rule 4)");

    aic_ucp_kraus phi; aic_assoc_ecstar ae; aic_dhom_B B; aic_dhom_v v;
    aic_factorize F; arb_t iso; double et[2];
    aic_fart_build_all(&F, &phi, &ae, &B, &v, iso, et,
                       kraus_re, kraus_im, n, r, eps, (slong) prec);

    *out_N = (int) F.N; *out_nB = (int) F.n_B; *out_dimB = (int) F.dim_A;
    *out_m = (int) v.B->num_blocks;
    assert(v.B->num_blocks <= (slong) n && "C4: m > n (blocks[] overflow)");
    for (slong l = 0; l < v.B->num_blocks; l++) blocks[l] = (int) v.B->d[l];

    aic_ucp_kraus dec, enc;
    aic_factorize_dec_kraus(&dec, &F, (slong) prec);   /* Dec = Delta*  (n_B x N) */
    aic_factorize_enc_kraus(&enc, &F, (slong) prec);   /* Enc = Upsilon* (N x n_B) */
    *out_rDec = (int) dec.r;
    *out_rEnc = (int) enc.r;
    aic_ucp_kraus_clear(&enc);
    aic_ucp_kraus_clear(&dec);

    aic_fart_build_free(&F, &phi, &ae, &B, &v, iso);
    return 0;
}
