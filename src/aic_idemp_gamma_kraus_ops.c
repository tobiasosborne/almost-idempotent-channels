/* aic_idemp_gamma_kraus_ops.c — small linear-algebra + density-matrix helpers for
 * the prop_Gamma Kraus extraction (bead aic-ynu, I3): the gauge-safe w-target
 * evaluation (aic_gk_wG), the tiny double-path least-squares solve (aic_gk_lstsq),
 * the Hermitization, and the gamma_j density-matrix certificate (Hermitian / PSD /
 * trace 1). The Kraus L_j construction is in aic_idemp_gamma_kraus_build.c and the
 * match-d->Gamma certification (THE CRUX TOOTH) in aic_idemp_gamma_certify.c — all
 * split for the <=200 LOC soft limit (Rule 10). See aic_idemp_gamma_kraus.c for
 * the narrative (.tex:2106-2122).
 */
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_idemp_gamma_internal.h"

double _Complex aic_gk_d(acb_srcptr z)
{
    double re = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    double im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    return re + im * I;
}

void aic_gk_wG(acb_mat_t out, const aic_idemp_decomp *d, const acb_mat_t X,
               slong prec)
{
    slong m = d->dim_M, dA = d->dim_A;
    acb_mat_t vx, g, wg;
    acb_mat_init(vx, m * m, 1);
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            acb_set(acb_mat_entry(vx, i * m + j, 0), acb_mat_entry(X, i, j));
    acb_mat_init(g, dA, 1);
    acb_mat_mul(g, d->Gamma, vx, prec);          /* A-coordinates of Gamma(X) */
    acb_mat_init(wg, m * m, 1);
    acb_mat_mul(wg, d->w, g, prec);              /* embed via w               */
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            acb_set(acb_mat_entry(out, i, j), acb_mat_entry(wg, i * m + j, 0));
    acb_mat_clear(wg);
    acb_mat_clear(g);
    acb_mat_clear(vx);
}

void aic_gk_lstsq(double _Complex *g, double *resid, const double _Complex *A,
                  const double _Complex *b, slong nR, slong nU)
{
    double _Complex *AhA = flint_calloc((size_t) nU * nU, sizeof(double _Complex));
    double _Complex *Ahb = flint_calloc((size_t) nU, sizeof(double _Complex));
    for (slong p = 0; p < nU; p++) {
        for (slong q = 0; q < nU; q++) {
            double _Complex s = 0;
            for (slong r = 0; r < nR; r++) s += conj(A[r * nU + p]) * A[r * nU + q];
            AhA[p * nU + q] = s;
        }
        double _Complex s = 0;
        for (slong r = 0; r < nR; r++) s += conj(A[r * nU + p]) * b[r];
        Ahb[p] = s;
    }
    /* Gauss-Jordan with partial pivoting on the nU x nU normal system. */
    for (slong c = 0; c < nU; c++) {
        slong piv = c;
        for (slong r = c + 1; r < nU; r++)
            if (cabs(AhA[r * nU + c]) > cabs(AhA[piv * nU + c])) piv = r;
        for (slong k = 0; k < nU; k++) {
            double _Complex t = AhA[c * nU + k];
            AhA[c * nU + k] = AhA[piv * nU + k];
            AhA[piv * nU + k] = t;
        }
        double _Complex t = Ahb[c];
        Ahb[c] = Ahb[piv];
        Ahb[piv] = t;
        double _Complex dpiv = AhA[c * nU + c];
        for (slong k = 0; k < nU; k++) AhA[c * nU + k] /= dpiv;
        Ahb[c] /= dpiv;
        for (slong r = 0; r < nU; r++)
            if (r != c) {
                double _Complex f = AhA[r * nU + c];
                for (slong k = 0; k < nU; k++) AhA[r * nU + k] -= f * AhA[c * nU + k];
                Ahb[r] -= f * Ahb[c];
            }
    }
    for (slong p = 0; p < nU; p++) g[p] = Ahb[p];

    /* consistency residual ||A g - b|| / (||b|| + 1) (~0 iff prop_Gamma form). */
    double rn = 0.0, bn = 0.0;
    for (slong r = 0; r < nR; r++) {
        double _Complex s = 0;
        for (slong p = 0; p < nU; p++) s += A[r * nU + p] * g[p];
        rn += cabs(s - b[r]) * cabs(s - b[r]);
        bn += cabs(b[r]) * cabs(b[r]);
    }
    *resid = sqrt(rn) / (sqrt(bn) + 1.0);
    flint_free(Ahb);
    flint_free(AhA);
}

void aic_gk_hermitize(acb_mat_t gamma, slong prec)
{
    slong dE = acb_mat_nrows(gamma);
    acb_mat_t gd;
    acb_mat_init(gd, dE, dE);
    acb_mat_conjugate_transpose(gd, gamma);
    acb_mat_add(gamma, gamma, gd, prec);
    {
        acb_t half;
        acb_init(half);
        acb_set_d(half, 0.5);
        acb_mat_scalar_mul_acb(gamma, gamma, half, prec);
        acb_clear(half);
    }
    acb_mat_clear(gd);
}

void aic_gk_certify_hermitian(const acb_mat_t gamma, slong j, slong prec)
{
    slong dE = acb_mat_nrows(gamma);
    arb_t tol, e;
    arb_init(tol);
    arb_init(e);
    /* tol 1e-9: gamma_j comes from a double-path least squares (floor ~1e-13);
     * an INVALID gamma (non-Hermitian, negative, mistraced) is O(1e-1)..O(1) off
     * (the diag(0.8,0.2) mutation gives ||match||~0.3). 1e-9 keeps full teeth. */
    arb_set_d(tol, 1e-9);
    acb_mat_t gd, diff;
    acb_mat_init(gd, dE, dE);
    acb_mat_init(diff, dE, dE);
    acb_mat_conjugate_transpose(gd, gamma);
    acb_mat_sub(diff, gamma, gd, prec);
    aic_mat_opnorm(e, diff, prec);                /* ||gamma - gamma^dag||_op */
    if (!arb_le(e, tol)) {
        flint_printf("aic_idemp_gamma_kraus: gamma_%wd not Hermitian "
                     "(||gamma - gamma^dag|| too large); fail loud (Rule 4)\n", j);
        flint_abort();
    }
    acb_mat_clear(diff);
    acb_mat_clear(gd);
    arb_clear(e);
    arb_clear(tol);
}

void aic_gk_certify_density(const acb_mat_t gamma, slong j, slong prec)
{
    slong dE = acb_mat_nrows(gamma);
    arb_t tol, e;
    arb_init(tol);
    arb_init(e);
    arb_set_d(tol, 1e-9);

    /* PSD: min eigenvalue = -max eig(-gamma) >= -tol. */
    acb_mat_t neg;
    acb_mat_init(neg, dE, dE);
    for (slong p = 0; p < dE; p++)
        for (slong q = 0; q < dE; q++)
            acb_neg(acb_mat_entry(neg, p, q), acb_mat_entry(gamma, p, q));
    aic_mat_herm_max_eig(e, neg, prec);           /* = -lambda_min(gamma) */
    if (!arb_le(e, tol)) {
        flint_printf("aic_idemp_gamma_kraus: gamma_%wd not PSD "
                     "(min eigenvalue < -1e-9); fail loud (Rule 4)\n", j);
        flint_abort();
    }
    acb_mat_clear(neg);

    /* trace 1: |Tr gamma - 1| */
    acb_t tr;
    acb_init(tr);
    acb_zero(tr);
    for (slong p = 0; p < dE; p++) acb_add(tr, tr, acb_mat_entry(gamma, p, p), prec);
    acb_sub_si(tr, tr, 1, prec);
    acb_abs(e, tr, prec);
    if (!arb_le(e, tol)) {
        flint_printf("aic_idemp_gamma_kraus: gamma_%wd not trace-1 "
                     "(|Tr gamma - 1| too large); fail loud (Rule 4)\n", j);
        flint_abort();
    }
    acb_clear(tr);
    arb_clear(e);
    arb_clear(tol);
}

