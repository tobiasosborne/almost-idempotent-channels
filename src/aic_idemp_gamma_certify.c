/* aic_idemp_gamma_certify.c — THE CRUX TOOTH of the prop_Gamma Kraus extraction
 * (bead aic-ynu, I3): certify (arb, production, fail-loud) that the assembled
 * prop_Gamma map equals the STORED d->Gamma on a basis of B(M). Split out of
 * aic_idemp_gamma_kraus_ops.c to keep each file near the <=200 LOC soft limit
 * (Rule 10). See aic_idemp_gamma_kraus.c for the narrative (.tex:2106-2122,:2113).
 *
 * The match is checked GAUGE-SAFELY: both A-elements are embedded in B(M) via the
 * *-monomorphism w and compared, ||w(Gamma_kraus(E_ab)) - w(d->Gamma(E_ab))||_op.
 * Gamma_j(E_ab) is computed BOTH via eq Gamma (the Tr_{E_j} form, gamma_j) AND via
 * the Choi/Kraus form (L_j) and the two must AGREE — pinning the gamma_j AND the
 * L_j build at once, and catching a wrong Tr_{E_j}-vs-Tr_{L_j} convention (the
 * forms disagree under a swapped trace; mutation-proven).
 */
#include <complex.h>

#include <flint/flint.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_idemp_gamma_internal.h"

/* Gamma_j(X) via eq Gamma (Tr_{E_j} form): Tr_E(W_j X W_j^dag (1_L (x) gamma_j)). */
static void gamma_j_tr(acb_mat_t Gj, const acb_mat_t Wj, const acb_mat_t gamma,
                       const acb_mat_t X, slong dL, slong dE, slong m, slong prec)
{
    slong db = dL * dE;
    acb_mat_t Wjd, WX, WXW, oneL, Ig, Z;
    acb_mat_init(Wjd, m, db);
    acb_mat_conjugate_transpose(Wjd, Wj);
    acb_mat_init(WX, db, m);
    acb_mat_mul(WX, Wj, X, prec);
    acb_mat_init(WXW, db, db);
    acb_mat_mul(WXW, WX, Wjd, prec);              /* W_j X W_j^dag */
    acb_mat_init(oneL, dL, dL);
    acb_mat_one(oneL);
    acb_mat_init(Ig, db, db);
    aic_mat_kronecker(Ig, oneL, gamma, prec);     /* 1_L (x) gamma_j */
    acb_mat_init(Z, db, db);
    acb_mat_mul(Z, WXW, Ig, prec);
    aic_mat_partial_trace_right(Gj, Z, dL, dE, prec);   /* Tr_E */
    acb_mat_clear(Z);
    acb_mat_clear(Ig);
    acb_mat_clear(oneL);
    acb_mat_clear(WXW);
    acb_mat_clear(WX);
    acb_mat_clear(Wjd);
}

/* Gamma_j(X) via the Choi/Kraus form: L_j^dag (X (x) 1_{F_j}) L_j (F_j = E_j). */
static void gamma_j_kraus(acb_mat_t Gj, const acb_mat_t Lj, const acb_mat_t X,
                          slong dL, slong dE, slong m, slong prec)
{
    acb_mat_t oneF, XxI, Ld, t;
    acb_mat_init(oneF, dE, dE);
    acb_mat_one(oneF);
    acb_mat_init(XxI, m * dE, m * dE);
    aic_mat_kronecker(XxI, X, oneF, prec);        /* X (x) 1_F (left-major) */
    acb_mat_init(Ld, dL, m * dE);
    acb_mat_conjugate_transpose(Ld, Lj);
    acb_mat_init(t, dL, m * dE);
    acb_mat_mul(t, Ld, XxI, prec);
    acb_mat_mul(Gj, t, Lj, prec);                 /* L^dag (X(x)1) L */
    acb_mat_clear(t);
    acb_mat_clear(Ld);
    acb_mat_clear(XxI);
    acb_mat_clear(oneF);
}

void aic_gk_certify_match(const aic_idemp_gamma *out, const aic_idemp_wedderburn *W,
                          const aic_idemp_decomp *d, slong prec)
{
    slong m = W->dim_M;
    arb_t tol, e;
    arb_init(tol);
    arb_init(e);
    /* match tol 1e-10 (the OUTPUT-CONTRACT gate). The clean match floors at the
     * double-path Wedderburn/Gamma extraction ~1e-13..1e-15; a wrong gamma_j or
     * Tr_{E_j} convention gives O(1e-1). 1e-10 keeps teeth with ~3 orders slack. */
    arb_set_d(tol, 1e-10);

    acb_mat_t X, wG, recon, Gj, GjK, oneE, GjXI, Wjd, t2, back, diff;
    acb_mat_init(X, m, m);
    acb_mat_init(wG, m, m);
    acb_mat_init(recon, m, m);
    acb_mat_init(diff, m, m);
    for (slong ia = 0; ia < m; ia++)
        for (slong ib = 0; ib < m; ib++) {
            acb_mat_zero(X);
            acb_set_si(acb_mat_entry(X, ia, ib), 1);
            aic_gk_wG(wG, d, X, prec);            /* stored: w(d->Gamma(E_ab)) */
            acb_mat_zero(recon);
            for (slong j = 0; j < W->num_blocks; j++) {
                slong dL = W->dim_L[j], dE = W->dim_E[j], db = dL * dE;
                /* eq-Gamma Tr_{E_j} form AND Choi L_j form must AGREE. */
                acb_mat_init(Gj, dL, dL);
                acb_mat_init(GjK, dL, dL);
                gamma_j_tr(Gj, W->W_j[j], out->gamma_j[j], X, dL, dE, m, prec);
                gamma_j_kraus(GjK, out->L_j[j], X, dL, dE, m, prec);
                acb_mat_sub(GjK, GjK, Gj, prec);
                aic_mat_opnorm(e, GjK, prec);
                if (!arb_le(e, tol)) {
                    flint_printf("aic_idemp_gamma_kraus: block %wd eq-Gamma vs Choi "
                                 "L_j DISAGREE (||Tr-form - Choi-form|| too large); "
                                 "fail loud (Rule 4)\n", j);
                    flint_abort();
                }
                /* re-embed Gamma_j(X) (x) 1_E via W_j^dag */
                acb_mat_init(oneE, dE, dE);
                acb_mat_one(oneE);
                acb_mat_init(GjXI, db, db);
                aic_mat_kronecker(GjXI, Gj, oneE, prec);
                acb_mat_init(Wjd, m, db);
                acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
                acb_mat_init(t2, m, db);
                acb_mat_mul(t2, Wjd, GjXI, prec);
                acb_mat_init(back, m, m);
                acb_mat_mul(back, t2, W->W_j[j], prec);
                acb_mat_add(recon, recon, back, prec);
                acb_mat_clear(back);
                acb_mat_clear(t2);
                acb_mat_clear(Wjd);
                acb_mat_clear(GjXI);
                acb_mat_clear(oneE);
                acb_mat_clear(GjK);
                acb_mat_clear(Gj);
            }
            acb_mat_sub(diff, recon, wG, prec);   /* prop_Gamma vs d->Gamma */
            aic_mat_opnorm(e, diff, prec);
            if (!arb_le(e, tol)) {
                flint_printf("aic_idemp_gamma_kraus: prop_Gamma Kraus map does NOT "
                             "match d->Gamma at E_{%wd,%wd} "
                             "(||w(Gamma_kraus) - w(d->Gamma)||_op too large); the "
                             "gamma_j or convention is wrong -> fail loud "
                             "(Rule 4, .tex:2113)\n", ia, ib);
                flint_abort();
            }
        }
    acb_mat_clear(diff);
    acb_mat_clear(recon);
    acb_mat_clear(wG);
    acb_mat_clear(X);
    arb_clear(e);
    arb_clear(tol);
}
