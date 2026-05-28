/* aic_idemp_cp.c — the Gamma-unital defect and the CP certificate for Gamma, via
 * Psi = C_M Lambda = w Gamma (bead aic-wuh; approximate_algebras.tex:2072, :2084).
 * Split out of aic_idemp_verify.c to keep each file <= 200 LOC (Rule 10).
 *
 * Gamma unital (.tex:2072). The proof shows Lambda(1_M) = Phi(1_H) = 1_H, hence
 * Gamma(1_M) = 1_A (so that Delta Gamma(1_M) = Lambda(1_M) = 1_H). We certify
 * ||Delta Gamma(1_M) - 1_H||_op = 0.
 *
 * Gamma CP — TWO certificates that together have teeth (.tex:2084, :2088).
 *
 *   (i) FACTORIZATION is CP. Psi = C_M Lambda : B(M) -> B(M),
 *         Psi(Y) = J_M^dag Phi(J_M Y J_M^dag) J_M,
 *       is a manifest composition of the CP maps C_M and Lambda. We build its
 *       Choi matrix (Convention A, dim_M^2 square) so the caller certifies CP
 *       with the tested aic_ucp_is_cp_choi. This certifies that the abstract
 *       algebra A's Choi-Effros product is CP. It is, however, BUILT FROM
 *       d->Lambda and d->J_M ALONE — it never reads d->Gamma, so on its own it
 *       says NOTHING about the stored Gamma (zeroing Gamma leaves it green).
 *
 *  (ii) STORED Gamma is tied to the CP object. The paper's logic (.tex:2084):
 *       Psi = C_M Lambda = w Gamma. Gamma : B(M) -> A is UCP because Psi = w
 *       Gamma is UCP and w is an injective *-homomorphism (completely isometric,
 *       .tex:2088): Psi is CP iff Gamma is. We CERTIFY the equation w Gamma =
 *       C_M Lambda as superoperators on B(M): aic_idemp_defect_wG_eq_CML returns
 *       ||w_mat @ Gamma_mat - Psi_CML||_op, which reads BOTH d->w and d->Gamma.
 *       A corrupt Gamma breaks this defect (scaling Gamma by 2 or zeroing a
 *       column makes it O(1)). Together (i)+(ii) certify Gamma is CP.
 *
 * Working through Psi on B(M) is GAUGE-INVARIANT and convention-safe (no abstract
 * A-coordinate Choi). Psi is also idempotent.
 *
 * Convention A Choi (matches include/aic_ucp.h): with matrix units F_kl of B(M),
 *     C_Psi[k*dM + a, l*dM + b] = Psi(F_kl)[a,b].
 * The vec-layout superoperator (for the tie-check) uses the SAME Psi(F_kl) but
 * indexed Psi_CML[a*dM + b, k*dM + l] = Psi(F_kl)[a,b] (row-major vec, matching
 * w_mat @ Gamma_mat).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_idemp.h"
#include "aic_mat.h"
#include "aic_ucp.h"

/* y (dA x 1) = Gamma(Y) coords for Y (dM x dM): Gamma[l,i*dM+j] Y[i,j]. */
static void gamma_apply(acb_mat_t y, const acb_mat_t Gamma, const acb_mat_t Y,
                        slong dA, slong dM, slong prec)
{
    acb_t acc, term;
    acb_init(acc);
    acb_init(term);
    for (slong l = 0; l < dA; l++) {
        acb_zero(acc);
        for (slong i = 0; i < dM; i++)
            for (slong j = 0; j < dM; j++) {
                acb_mul(term, acb_mat_entry(Gamma, l, i * dM + j),
                        acb_mat_entry(Y, i, j), prec);
                acb_add(acc, acc, term, prec);
            }
        acb_set(acb_mat_entry(y, l, 0), acc);
    }
    acb_clear(term);
    acb_clear(acc);
}

/* ||Delta Gamma(1_M) - 1_H||_op (.tex:2072). */
void aic_idemp_defect_gamma_unital(arb_t out, const aic_idemp_decomp *d,
                                   slong prec)
{
    slong n = d->n, dM = d->dim_M, dA = d->dim_A;
    acb_mat_t oneM, coords, v, Z, oneH, diff;
    acb_mat_init(oneM, dM, dM);
    acb_mat_init(coords, dA, 1);
    acb_mat_init(v, n * n, 1);
    acb_mat_init(Z, n, n);
    acb_mat_init(oneH, n, n);
    acb_mat_init(diff, n, n);
    acb_mat_one(oneM);
    gamma_apply(coords, d->Gamma, oneM, dA, dM, prec);   /* Gamma(1_M) coords */
    acb_mat_mul(v, d->Delta, coords, prec);              /* Delta(coords) */
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(Z, a, b), acb_mat_entry(v, a * n + b, 0));
    acb_mat_one(oneH);
    acb_mat_sub(diff, Z, oneH, prec);
    aic_mat_opnorm(out, diff, prec);
    acb_mat_clear(diff);
    acb_mat_clear(oneH);
    acb_mat_clear(Z);
    acb_mat_clear(v);
    acb_mat_clear(coords);
    acb_mat_clear(oneM);
}

/* Psi(F_kl) = J_M^dag Lambda(F_kl) J_M (dM x dM), with Lambda(F_kl) read as
 * column (k*dM+l) of d->Lambda (n x n). The single-column kernel shared by the
 * Choi build and the vec-layout superoperator. */
static void psi_col(acb_mat_t PsiF, const aic_idemp_decomp *d, const acb_mat_t JMd,
                    slong k, slong l, slong prec)
{
    slong n = d->n, dM = d->dim_M;
    acb_mat_t LamF, t1;
    acb_mat_init(LamF, n, n);
    acb_mat_init(t1, dM, n);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(LamF, a, b),
                    acb_mat_entry(d->Lambda, a * n + b, k * dM + l));
    acb_mat_mul(t1, JMd, LamF, prec);     /* J_M^dag Lambda(F) */
    acb_mat_mul(PsiF, t1, d->J_M, prec);  /* J_M^dag Lambda(F) J_M */
    acb_mat_clear(t1);
    acb_mat_clear(LamF);
}

/* Choi matrix of Psi = C_M Lambda : B(M) -> B(M), Convention A. */
void aic_idemp_psi_choi(acb_mat_t C, const aic_idemp_decomp *d, slong prec)
{
    slong n = d->n, dM = d->dim_M;
    assert(acb_mat_nrows(C) == dM * dM && acb_mat_ncols(C) == dM * dM &&
           "aic_idemp_psi_choi: C must be (dim_M*dim_M) square");

    acb_mat_t JMd, PsiF;
    acb_mat_init(JMd, dM, n);
    acb_mat_init(PsiF, dM, dM);
    acb_mat_conjugate_transpose(JMd, d->J_M);

    for (slong k = 0; k < dM; k++)
        for (slong l = 0; l < dM; l++) {
            psi_col(PsiF, d, JMd, k, l, prec);
            for (slong a = 0; a < dM; a++)
                for (slong b = 0; b < dM; b++)
                    acb_set(acb_mat_entry(C, k * dM + a, l * dM + b),
                            acb_mat_entry(PsiF, a, b));
        }

    acb_mat_clear(PsiF);
    acb_mat_clear(JMd);
}

/* ||w Gamma - C_M Lambda||_op (.tex:2084: Psi = C_M Lambda = w Gamma), as
 * superoperators on B(M). TEETH for the Gamma-CP certificate: reads BOTH d->w
 * and d->Gamma and ties them to the CP-certified Psi = C_M Lambda (which reads
 * only d->Lambda, d->J_M). Both are (dM^2 x dM^2) in the row-major vec layout
 * Psi[a*dM+b, k*dM+l] = Psi(F_kl)[a,b]; w_mat @ Gamma_mat is already that layout. */
void aic_idemp_defect_wG_eq_CML(arb_t out, const aic_idemp_decomp *d, slong prec)
{
    slong n = d->n, dM = d->dim_M, dM2 = dM * dM;
    acb_mat_t JMd, PsiF, PsiCML, wG, diff;
    acb_mat_init(JMd, dM, n);
    acb_mat_init(PsiF, dM, dM);
    acb_mat_init(PsiCML, dM2, dM2);
    acb_mat_init(wG, dM2, dM2);
    acb_mat_init(diff, dM2, dM2);
    acb_mat_conjugate_transpose(JMd, d->J_M);

    /* Psi_CML[a*dM+b, k*dM+l] = (J_M^dag Lambda(F_kl) J_M)[a,b]. */
    for (slong k = 0; k < dM; k++)
        for (slong l = 0; l < dM; l++) {
            psi_col(PsiF, d, JMd, k, l, prec);
            for (slong a = 0; a < dM; a++)
                for (slong b = 0; b < dM; b++)
                    acb_set(acb_mat_entry(PsiCML, a * dM + b, k * dM + l),
                            acb_mat_entry(PsiF, a, b));
        }

    acb_mat_mul(wG, d->w, d->Gamma, prec);   /* (dM^2 x dA)(dA x dM^2) = w Gamma */
    acb_mat_sub(diff, wG, PsiCML, prec);
    aic_mat_opnorm(out, diff, prec);

    acb_mat_clear(diff);
    acb_mat_clear(wG);
    acb_mat_clear(PsiCML);
    acb_mat_clear(PsiF);
    acb_mat_clear(JMd);
}
