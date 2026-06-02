/* aic_idemp_delta_kraus_certify.c — the certifications of the prop_Delta Kraus
 * extraction (bead aic-ynu, I4): the M-part unitality (Delta_M(1_A) = Pi_M) and THE
 * CRUX TOOTH — that the full Delta-Kraus map equals the STORED d->Delta on the
 * A-basis. Split out of aic_idemp_delta_kraus_ops.c to keep each file near the
 * <=200 LOC soft limit (Rule 10), mirroring the gamma_kraus _certify.c split. See
 * aic_idemp_delta_kraus.c for the narrative (.tex:2098-2103,:2100).
 *
 * Both certifications run IN PRODUCTION (arb, fail-loud, Rule 4): a wrong K_{j,c}
 * build, W_j usage, or Sigma read-off gives an O(1) defect and aborts with a
 * precise message — corrupted output is worse than a crash (CLAUDE.md Rule 4).
 */
#include <flint/flint.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_idemp_delta_internal.h"

/* reshape column k of d->Delta (length n^2, row-major vec) into the n x n B_k. */
static void delta_col_to_mat(acb_mat_t B, const acb_mat_t Delta, slong k, slong n)
{
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(B, a, b), acb_mat_entry(Delta, a * n + b, k));
}

void aic_dk_certify_unital(const aic_idemp_delta *out, const aic_idemp_decomp *d,
                           slong prec)
{
    slong n = out->n;
    /* A_hat = 1_{dL_tot} is the abstract unit 1_A; w(1_A) = 1_M, so the M-part
     * Kraus map should give Delta_M(1_A) = J_M 1_M J_M^dag = Pi_M. */
    acb_mat_t one, Mpart, diff;
    arb_t tol, e;
    acb_mat_init(one, out->dL_tot, out->dL_tot);
    acb_mat_one(one);
    acb_mat_init(Mpart, n, n);
    aic_dk_apply_Mpart(Mpart, out, one, prec);
    acb_mat_init(diff, n, n);
    acb_mat_sub(diff, Mpart, d->Pi_M, prec);
    arb_init(tol);
    arb_init(e);
    arb_set_d(tol, 1e-9);
    aic_mat_opnorm(e, diff, prec);
    if (!arb_le(e, tol)) {
        flint_printf("aic_idemp_delta_kraus: M-part not unital onto M "
                     "(||Delta_M(1_A) - Pi_M||_op too large); the K_{j,c} build or "
                     "W_j usage is wrong -> fail loud (Rule 4, .tex:2100)\n");
        flint_abort();
    }
    arb_clear(e);
    arb_clear(tol);
    acb_mat_clear(diff);
    acb_mat_clear(Mpart);
    acb_mat_clear(one);
}

void aic_dk_certify_match(const aic_idemp_delta *out, const aic_idemp_wedderburn *W,
                          const aic_idemp_decomp *d, slong prec)
{
    slong n = out->n, m = out->dim_M, dA = out->dim_A;
    arb_t tol, e;
    arb_init(tol);
    arb_init(e);
    /* match tol 1e-10 (the OUTPUT-CONTRACT gate). The clean match floors at the
     * double-path Wedderburn/Delta extraction ~1e-13..1e-15; a wrong K_{j,c} or
     * Sigma read-off gives O(1). 1e-10 keeps teeth with ~3 orders slack. */
    arb_set_d(tol, 1e-10);

    acb_mat_t a, wA, Ahat, Mpart, MperpPart, recon, Bk, diff;
    acb_mat_init(a, dA, 1);
    acb_mat_init(wA, m, m);
    acb_mat_init(Ahat, out->dL_tot, out->dL_tot);
    acb_mat_init(Mpart, n, n);
    acb_mat_init(MperpPart, n, n);
    acb_mat_init(recon, n, n);
    acb_mat_init(Bk, n, n);
    acb_mat_init(diff, n, n);
    for (slong k = 0; k < dA; k++) {
        /* A = B_k: A-coordinate a = e_k; wA = w(B_k) = reshape of d->w column k. */
        acb_mat_zero(a);
        acb_set_si(acb_mat_entry(a, k, 0), 1);
        for (slong i = 0; i < m; i++)
            for (slong jj = 0; jj < m; jj++)
                acb_set(acb_mat_entry(wA, i, jj), acb_mat_entry(d->w, i * m + jj, k));

        aic_dk_Ahat_from_wA(Ahat, out, W, wA, prec);
        aic_dk_apply_Mpart(Mpart, out, Ahat, prec);        /* M-part via K_{j,c} */
        aic_dk_apply_Sigma(MperpPart, out, d, a, prec);    /* M^perp via Sigma   */
        acb_mat_add(recon, Mpart, MperpPart, prec);

        delta_col_to_mat(Bk, d->Delta, k, n);              /* stored Delta(B_k) */
        acb_mat_sub(diff, recon, Bk, prec);
        aic_mat_opnorm(e, diff, prec);
        if (!arb_le(e, tol)) {
            flint_printf("aic_idemp_delta_kraus: Delta-Kraus map does NOT match "
                         "d->Delta at B_%wd (||Delta_kraus(B_k) - B_k||_op too "
                         "large); the K_{j,c}, W_j, or Sigma is wrong -> fail loud "
                         "(Rule 4, .tex:2100)\n", k);
            flint_abort();
        }
    }
    acb_mat_clear(diff);
    acb_mat_clear(Bk);
    acb_mat_clear(recon);
    acb_mat_clear(MperpPart);
    acb_mat_clear(Mpart);
    acb_mat_clear(Ahat);
    acb_mat_clear(wA);
    acb_mat_clear(a);
    arb_clear(e);
    arb_clear(tol);
}
