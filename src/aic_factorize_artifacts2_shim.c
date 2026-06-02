/* aic_factorize_artifacts2_shim.c — flat-double ccall SHIM for the th_factorization
 * ARTIFACTS, shim C5 (aic_factorize_artifacts_d). Split from
 * aic_factorize_artifacts_shim.c (which holds C4 + the SHARED pipeline build helper
 * aic_fart_build_all/_free, aic_factorize_artifacts_internal.h) per Rule 10
 * (each <=200 LOC). See include/aic/aic_factorize_artifacts_shim.h for the
 * ABI/layout contract and the Appendix-B1 channel-direction binding.
 *
 * PROVENANCE (Law 1). C5 rebuilds the SAME pipeline as C4 (aic_fart_build_all),
 * then the F4.1 dual read-off (aic_factorize.h:306-316) and the F4.1 eig-free
 * brackets (aic_factorize.h:301-304):
 *   Dec = Delta*   = aic_factorize_dec_kraus  (n_B x N, layout [a*(nB*N)+i*N+j])
 *   Enc = Upsilon* = aic_factorize_enc_kraus  (N x n_B, layout [a*(N*nB)+i*nB+j])
 *   ||Delta Upsilon - Phi||_cb : aic_factorize_eigfree_delups (.tex:2733)
 *   ||Upsilon Delta - 1_B||_cb : aic_factorize_eigfree_upsdel (.tex:2739)
 * APPENDIX B1 (LOAD-BEARING): dec_* MUST come from dec_kraus and enc_* from
 * enc_kraus (DO NOT cross — the §C13 dual-direction trap; the [J4] round-trip
 * mutation tooth catches a swap). The brackets are rounded OUTWARD (FLOOR/CEIL) so
 * they stay RIGOROUS after the cast (FINDINGS §C14: the dual read-off projects the
 * O(eta^2) cone defect inside dec_kraus/enc_kraus; multi-block eta>0 needs it).
 *
 * Number path: arb (acb_mat) internally; only the boundary is double. Fail-loud
 * (Rule 4): bad shape aborts; the rebuilt rDec/rEnc must match the caller's C4
 * sizes (else the caller mis-sized — fail loud). <=200 LOC (Rule 10).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_assoc.h"
#include "aic/aic_dhom.h"
#include "aic/aic_factorize.h"
#include "aic/aic_factorize_artifacts_shim.h"
#include "aic/aic_ucp.h"
#include "aic_factorize_artifacts_internal.h"

/* Copy an aic_ucp_kraus (r ops, each dim_K x dim_H) to flat K-major
 * [a*(dim_K*dim_H) + i*dim_H + j] double arrays (midpoints). */
static void kraus_to_flat(double *re, double *im, const aic_ucp_kraus *K)
{
    slong dK = K->dim_K, dH = K->dim_H;
    for (slong a = 0; a < K->r; a++)
        for (slong i = 0; i < dK; i++)
            for (slong j = 0; j < dH; j++) {
                slong off = a * (dK * dH) + i * dH + j;
                acb_srcptr e = acb_mat_entry(K->K[a], i, j);
                re[off] = arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR);
                im[off] = arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR);
            }
}

/* arb ball -> [lo,hi] doubles, FLOOR/CEIL outward (RIGOROUS, like aic_ucp_shim.c). */
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

/* C5 — emit the channel Kraus (Dec from dec_kraus, Enc from enc_kraus — Appendix
 * B1, DO NOT cross) + the FLOOR/CEIL round-trip brackets + out_eta. */
int aic_factorize_artifacts_d(double *dec_re, double *dec_im, int rDec,
                              double *enc_re, double *enc_im, int rEnc,
                              double *delups, double *upsdel,
                              double *out_eta,
                              const double *kraus_re, const double *kraus_im,
                              int n, int r, double eps, int prec)
{
    assert(n >= 1 && r >= 1 && prec >= 2 && rDec >= 1 && rEnc >= 1 &&
           "aic_factorize_artifacts_d: bad shape (Rule 4)");
    assert(dec_re && dec_im && enc_re && enc_im && delups && upsdel && out_eta &&
           kraus_re && kraus_im && "aic_factorize_artifacts_d: null out/in (Rule 4)");

    aic_ucp_kraus phi; aic_assoc_ecstar ae; aic_dhom_B B; aic_dhom_v v;
    aic_factorize F; arb_t iso;
    aic_fart_build_all(&F, &phi, &ae, &B, &v, iso, out_eta,
                       kraus_re, kraus_im, n, r, eps, (slong) prec);

    /* Dec = Delta* (B1): emit dec_* from dec_kraus (n_B x N). */
    aic_ucp_kraus dec, enc;
    aic_factorize_dec_kraus(&dec, &F, (slong) prec);
    aic_factorize_enc_kraus(&enc, &F, (slong) prec);   /* Enc = Upsilon* (N x n_B) */
    assert((int) dec.r == rDec &&
           "aic_factorize_artifacts_d: rebuilt rDec != caller's (mis-sized via C4)");
    assert((int) enc.r == rEnc &&
           "aic_factorize_artifacts_d: rebuilt rEnc != caller's (mis-sized via C4)");
    kraus_to_flat(dec_re, dec_im, &dec);
    kraus_to_flat(enc_re, enc_im, &enc);
    aic_ucp_kraus_clear(&enc);
    aic_ucp_kraus_clear(&dec);

    /* round-trip brackets (.tex:2733/2739), FLOOR/CEIL outward. */
    {
        arb_t loB, hiB;
        arb_init(loB);
        arb_init(hiB);
        aic_factorize_eigfree_delups(loB, hiB, &F, (slong) prec);
        twoball_to_bracket(delups, loB, hiB, (slong) prec);
        aic_factorize_eigfree_upsdel(loB, hiB, &F, (slong) prec);
        twoball_to_bracket(upsdel, loB, hiB, (slong) prec);
        arb_clear(hiB);
        arb_clear(loB);
    }

    aic_fart_build_free(&F, &phi, &ae, &B, &v, iso);
    return 0;
}
