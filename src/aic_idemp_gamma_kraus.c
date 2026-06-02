/* aic_idemp_gamma_kraus.c — the explicit prop_Gamma Kraus form of the conditional
 * expectation Gamma : B(M) -> A (bead aic-ynu, I3). Realizes prop_Gamma
 * (approximate_algebras.tex:2106-2122) for the EXACTLY-idempotent (eta=0) case.
 *
 * approximate_algebras.tex:2108-2113 — eq Gamma (verbatim core):
 *   "Gamma_j(X) = Tr_{E_j}( W_j X W_j^dag (1_{L_j} (x) gamma_j) )
 *    for some density matrices gamma_j on E_j."
 * approximate_algebras.tex:2118-2122 — the Choi/Kraus form:
 *   "Gamma_j(X) = L_j^dag (X (x) 1_{F_j}) L_j,  L_j = (W_j^dag (x) 1_{F_j})(1_{L_j}
 *    (x) xi_j),  gamma_j = Tr_{F_j}(xi_j xi_j^dag)."
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (Law 3). The paper PROVES every UCP Gamma
 * with Gamma w = 1_A has eq-Gamma form (a Choi-form + algebraic-elimination
 * existence argument). In finite dim we EXTRACT the gamma_j that reproduces the
 * STORED d->Gamma (which turns out maximally-mixed; paper/FINDINGS.md §C20): eq
 * Gamma is LINEAR in gamma_j, so we feed a spanning set of
 * B(M) inputs, read the per-block target Gamma_j(X) from d->Gamma through the
 * *-monomorphism w, and least-squares solve for gamma_j (unique: W_j is an
 * isometry onto block j, so gamma_j |-> Gamma_j is injective). gamma_j is then
 * CERTIFIED a density matrix (Hermitian / PSD / trace 1) and the assembled Kraus
 * map is CERTIFIED to MATCH d->Gamma (the crux tooth) — both in arb, in
 * production, fail-loud (Rule 4). A genuine mismatch (no valid gamma_j reproduces
 * d->Gamma) is a W_j/convention bug, escalated, not faked (.tex:2113).
 *
 * THE TARGET Gamma_j(X) from d->Gamma (gauge-safe). d->Gamma outputs an A-element
 * as a dim_A coordinate vector g (in the {B_k} basis); w embeds it in B(M) as
 * wG(X) = reshape(w * g) (an m x m matrix, m = dim_M). The block-j component is
 * read by the isometry W_j: W_j wG(X) W_j^dag = Gamma_j(X) (x) 1_{E_j}
 * (W_j W_k^dag = delta_jk 1), so Gamma_j(X) = Tr_{E_j}(W_j wG(X) W_j^dag)/dim E_j.
 * Working through w is gauge-invariant (no abstract A-coordinate convention bites).
 *
 * F_j = E_j PURIFICATION. gamma_j = sum_c p_c |v_c><v_c| (Hermitian eig); xi_j =
 * sum_c sqrt(p_c) |v_c>_E (x) |c>_F gives Tr_{F_j}(xi_j xi_j^dag) = gamma_j, so
 * L_j = (W_j^dag (x) 1_{F_j})(1_{L_j} (x) xi_j) is (dim_M*dim_E) x dim_L.
 *
 * NUMBER PATH. Like aic_idemp / the Wedderburn extraction, the gamma_j solve runs
 * in the LAPACK double path (a tiny dim_E^2 x dim_E^2 normal-equations solve), and
 * the gamma_j / Kraus matrices are written as acb_mat at `prec`; the
 * density-matrix validity, the L_j isometry, and the match-d->Gamma defect are
 * CERTIFIED in arb. See aic_idemp_gamma_kraus_ops.c for the small helpers
 * (the w-target evaluation, the eig-purification, the certification).
 */
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_idemp_gamma_internal.h"

/* gamma_j extraction: least-squares solve of  Tr_{E_j}(W_j X W_j^dag (1(x)gamma_j))
 * = Gamma_j(X)  over X = matrix units of B(M). Returns gamma_j (dE x dE) as a
 * double array (row-major), and the consistency residual (||A gamma - b|| / scale)
 * in *resid — large iff d->Gamma is NOT of prop_Gamma form for this block. */
static void extract_gamma_block(double _Complex *gamma, double *resid,
                                const aic_idemp_wedderburn *W,
                                const aic_idemp_decomp *d, slong j, slong prec)
{
    slong m = W->dim_M, dL = W->dim_L[j], dE = W->dim_E[j], db = dL * dE;
    slong nU = dE * dE;                  /* unknowns gamma[c',c] at col c'*dE+c */
    slong nR = (m * m) * (dL * dL);      /* equations: (input unit) x (output entry) */

    double _Complex *A = flint_calloc((size_t) nR * nU, sizeof(double _Complex));
    double _Complex *b = flint_calloc((size_t) nR, sizeof(double _Complex));

    acb_mat_t Wjd, X, wG, WX, WXW, Tr;
    acb_mat_init(Wjd, m, db);
    acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
    acb_mat_init(X, m, m);
    acb_mat_init(wG, m, m);
    acb_mat_init(WX, db, m);
    acb_mat_init(WXW, db, db);
    acb_mat_init(Tr, dL, dL);

    slong row = 0;
    for (slong ia = 0; ia < m; ia++)
        for (slong ib = 0; ib < m; ib++) {
            acb_mat_zero(X);
            acb_set_si(acb_mat_entry(X, ia, ib), 1);

            /* target: Gamma_j(E_ab) = Tr_E(W_j wG(E_ab) W_j^dag)/dE */
            aic_gk_wG(wG, d, X, prec);
            acb_mat_mul(WX, W->W_j[j], wG, prec);
            acb_mat_mul(WXW, WX, Wjd, prec);
            aic_mat_partial_trace_right(Tr, WXW, dL, dE, prec);

            /* coefficient block C = W_j E_ab W_j^dag (so eq Gamma is linear) */
            acb_mat_mul(WX, W->W_j[j], X, prec);
            acb_mat_mul(WXW, WX, Wjd, prec);

            for (slong a0 = 0; a0 < dL; a0++)
                for (slong a1 = 0; a1 < dL; a1++) {
                    for (slong cp = 0; cp < dE; cp++)
                        for (slong c = 0; c < dE; c++)
                            A[row * nU + cp * dE + c] +=
                                aic_gk_d(acb_mat_entry(WXW, a0 * dE + c, a1 * dE + cp));
                    b[row] = aic_gk_d(acb_mat_entry(Tr, a0, a1)) / (double) dE;
                    row++;
                }
        }
    aic_gk_lstsq(gamma, resid, A, b, nR, nU);

    acb_mat_clear(Tr);
    acb_mat_clear(WXW);
    acb_mat_clear(WX);
    acb_mat_clear(wG);
    acb_mat_clear(X);
    acb_mat_clear(Wjd);
    flint_free(b);
    flint_free(A);
}

void aic_idemp_gamma_kraus(aic_idemp_gamma *out, const aic_idemp_wedderburn *W,
                           const aic_idemp_decomp *d, slong prec)
{
    slong m = W->dim_M, nb = W->num_blocks;
    assert(m == d->dim_M && W->dim_A == d->dim_A && "gamma_kraus: W/d mismatch");

    out->num_blocks = nb;
    out->dim_M = m;
    out->dim_A = W->dim_A;
    out->dim_L = flint_malloc((size_t) nb * sizeof(slong));
    out->dim_E = flint_malloc((size_t) nb * sizeof(slong));
    out->gamma_j = flint_malloc((size_t) nb * sizeof(acb_mat_t));
    out->L_j = flint_malloc((size_t) nb * sizeof(acb_mat_t));

    for (slong j = 0; j < nb; j++) {
        slong dL = W->dim_L[j], dE = W->dim_E[j];
        out->dim_L[j] = dL;
        out->dim_E[j] = dE;

        /* (1) extract gamma_j to MATCH d->Gamma (double-path least squares). */
        double _Complex *g = flint_malloc((size_t) (dE * dE) * sizeof(double _Complex));
        double resid = 0.0;
        extract_gamma_block(g, &resid, W, d, j, prec);

        /* fail loud if NO valid gamma reproduces this block (d->Gamma not prop_Gamma
         * form -> a W_j/convention bug, .tex:2113). The lstsq consistency residual
         * is O(1) when inconsistent, ~1e-13 when consistent. 1e-7 sits in the gulf. */
        if (resid > 1e-7) {
            flint_printf("aic_idemp_gamma_kraus: block %wd of d->Gamma is NOT of "
                         "prop_Gamma form (least-squares residual %.3e > 1e-7); the "
                         "W_j or the Tr_{E_j} convention is wrong -> escalate "
                         "(Rule 4, .tex:2113)\n", j, resid);
            flint_abort();
        }

        acb_mat_init(out->gamma_j[j], dE, dE);
        for (slong p = 0; p < dE; p++)
            for (slong q = 0; q < dE; q++)
                acb_set_d_d(acb_mat_entry(out->gamma_j[j], p, q),
                            creal(g[p * dE + q]), cimag(g[p * dE + q]));
        flint_free(g);

        /* (2) certify gamma_j is a valid density matrix. The RAW least-squares
         * solve has a ~1e-13 non-Hermitian part (double-path noise); we check that
         * asymmetry is below tol (the Hermiticity tooth), THEN symmetrize so the
         * PSD eig routine (which enforces a TIGHT relative Hermiticity tol) accepts
         * it, THEN certify PSD + trace 1 on the symmetrized gamma_j. */
        aic_gk_certify_hermitian(out->gamma_j[j], j, prec);   /* asymmetry tooth */
        aic_gk_hermitize(out->gamma_j[j], prec);
        aic_gk_certify_density(out->gamma_j[j], j, prec);     /* PSD + trace 1 */

        /* (3) build the Kraus L_j = (W_j^dag (x) 1_F)(1_L (x) xi_j), F_j = E_j. */
        acb_mat_init(out->L_j[j], m * dE, dL);
        aic_gk_build_Lj(out->L_j[j], out->gamma_j[j], W->W_j[j], dL, dE, m, j, prec);
    }

    /* (4) THE CRUX TOOTH: certify the prop_Gamma Kraus map == d->Gamma on a basis
     * of B(M) (gauge-safe, embedded via w). Fails loud on any mismatch. */
    aic_gk_certify_match(out, W, d, prec);
}

void aic_idemp_gamma_clear(aic_idemp_gamma *out)
{
    if (out == NULL || out->gamma_j == NULL) return;
    for (slong j = 0; j < out->num_blocks; j++) {
        acb_mat_clear(out->gamma_j[j]);
        acb_mat_clear(out->L_j[j]);
    }
    flint_free(out->gamma_j);
    flint_free(out->L_j);
    flint_free(out->dim_L);
    flint_free(out->dim_E);
    out->gamma_j = NULL;
    out->L_j = NULL;
    out->dim_L = NULL;
    out->dim_E = NULL;
    out->num_blocks = 0;
}
