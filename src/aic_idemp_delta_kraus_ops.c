/* aic_idemp_delta_kraus_ops.c — helpers for the prop_Delta Kraus extraction (bead
 * aic-ynu, I4): the gauge-safe A_hat read-off (A_j = Tr_{E_j}(W_j w(A) W_j^dag)/dE),
 * the M-part Kraus application (sum_{j,c} K_{j,c} A_hat K_{j,c}^dag), and the
 * M^perp Sigma read-off + embedding. The unitality + match (crux tooth)
 * certifications are in aic_idemp_delta_kraus_certify.c; all split for the <=200 LOC
 * soft limit (Rule 10). See aic_idemp_delta_kraus.c for the narrative (.tex:2098-2103).
 */
#include <flint/flint.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_idemp_delta_internal.h"

void aic_dk_Ahat_from_wA(acb_mat_t Ahat, const aic_idemp_delta *out,
                         const aic_idemp_wedderburn *W, const acb_mat_t wA,
                         slong prec)
{
    slong m = out->dim_M;
    acb_mat_zero(Ahat);
    for (slong j = 0; j < out->num_blocks; j++) {
        slong dL = out->dim_L[j], dE = out->dim_E[j], db = dL * dE, off = out->off_L[j];
        acb_mat_t Wjd, WX, WXW, Aj;
        acb_mat_init(Wjd, m, db);
        acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
        acb_mat_init(WX, db, m);
        acb_mat_mul(WX, W->W_j[j], wA, prec);          /* W_j wA */
        acb_mat_init(WXW, db, db);
        acb_mat_mul(WXW, WX, Wjd, prec);               /* W_j wA W_j^dag */
        acb_mat_init(Aj, dL, dL);
        aic_mat_partial_trace_right(Aj, WXW, dL, dE, prec);   /* Tr_{E_j}(...) */
        for (slong a = 0; a < dL; a++)
            for (slong b = 0; b < dL; b++) {
                acb_div_si(acb_mat_entry(Aj, a, b), acb_mat_entry(Aj, a, b), dE, prec);
                acb_set(acb_mat_entry(Ahat, off + a, off + b),
                        acb_mat_entry(Aj, a, b));
            }
        acb_mat_clear(Aj);
        acb_mat_clear(WXW);
        acb_mat_clear(WX);
        acb_mat_clear(Wjd);
    }
}

void aic_dk_apply_Mpart(acb_mat_t Mpart, const aic_idemp_delta *out,
                        const acb_mat_t Ahat, slong prec)
{
    slong n = out->n, dL_tot = out->dL_tot;
    acb_mat_zero(Mpart);
    acb_mat_t Kd, KA, KAK;
    acb_mat_init(Kd, dL_tot, n);
    acb_mat_init(KA, n, dL_tot);
    acb_mat_init(KAK, n, n);
    for (slong i = 0; i < out->n_kraus; i++) {
        acb_mat_mul(KA, out->K[i], Ahat, prec);        /* K A_hat */
        acb_mat_conjugate_transpose(Kd, out->K[i]);
        acb_mat_mul(KAK, KA, Kd, prec);                /* K A_hat K^dag */
        acb_mat_add(Mpart, Mpart, KAK, prec);
    }
    acb_mat_clear(KAK);
    acb_mat_clear(KA);
    acb_mat_clear(Kd);
}

/* reshape column k of d->Delta (length n^2, row-major vec) into the n x n B_k. */
static void delta_col_to_mat(acb_mat_t B, const acb_mat_t Delta, slong k, slong n)
{
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(B, a, b), acb_mat_entry(Delta, a * n + b, k));
}

void aic_dk_build_sigma(aic_idemp_delta *out, const aic_idemp_decomp *d, slong prec)
{
    slong n = d->n, dmp = out->dim_Mperp, dA = d->dim_A;
    if (dmp == 0) return;                               /* M = H, no M^perp part */
    acb_mat_init(out->Sigma, dmp * dmp, dA);
    acb_mat_t Jpd, B, JpB, JpBJ;
    acb_mat_init(Jpd, dmp, n);
    acb_mat_conjugate_transpose(Jpd, d->J_Mperp);
    acb_mat_init(B, n, n);
    acb_mat_init(JpB, dmp, n);
    acb_mat_init(JpBJ, dmp, dmp);
    for (slong k = 0; k < dA; k++) {
        delta_col_to_mat(B, d->Delta, k, n);
        acb_mat_mul(JpB, Jpd, B, prec);                /* J_Mperp^dag B_k */
        acb_mat_mul(JpBJ, JpB, d->J_Mperp, prec);      /* J_Mperp^dag B_k J_Mperp */
        for (slong p = 0; p < dmp; p++)
            for (slong q = 0; q < dmp; q++)
                acb_set(acb_mat_entry(out->Sigma, p * dmp + q, k),
                        acb_mat_entry(JpBJ, p, q));
    }
    acb_mat_clear(JpBJ);
    acb_mat_clear(JpB);
    acb_mat_clear(B);
    acb_mat_clear(Jpd);
}

void aic_dk_apply_Sigma(acb_mat_t Mperp_part, const aic_idemp_delta *out,
                        const aic_idemp_decomp *d, const acb_mat_t a, slong prec)
{
    slong n = out->n, dmp = out->dim_Mperp;
    acb_mat_zero(Mperp_part);
    if (dmp == 0) return;
    acb_mat_t s, S, JpS, JpSJ, Jpd;
    acb_mat_init(s, dmp * dmp, 1);
    acb_mat_mul(s, out->Sigma, a, prec);               /* Sigma(a) coords */
    acb_mat_init(S, dmp, dmp);
    for (slong p = 0; p < dmp; p++)
        for (slong q = 0; q < dmp; q++)
            acb_set(acb_mat_entry(S, p, q), acb_mat_entry(s, p * dmp + q, 0));
    acb_mat_init(JpS, n, dmp);
    acb_mat_mul(JpS, d->J_Mperp, S, prec);             /* J_Mperp Sigma(a) */
    acb_mat_init(Jpd, dmp, n);
    acb_mat_conjugate_transpose(Jpd, d->J_Mperp);
    acb_mat_init(JpSJ, n, n);
    acb_mat_mul(JpSJ, JpS, Jpd, prec);                 /* J_Mperp Sigma(a) J_Mperp^dag */
    acb_mat_set(Mperp_part, JpSJ);
    acb_mat_clear(JpSJ);
    acb_mat_clear(Jpd);
    acb_mat_clear(JpS);
    acb_mat_clear(S);
    acb_mat_clear(s);
}
