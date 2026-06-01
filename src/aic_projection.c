/* aic_projection.c — the nontrivial-projection finder (bead aic-mqf), orchestrator
 * + Step 5 (certification). Realizes lem_nontriv_projection
 * (approximate_algebras.tex:931-932). See include/aic_projection.h for the full
 * route, the rem_X2-safety argument, and the audition note (Law 4). Steps 1-2
 * (H-search + gap) live in src/aic_projection_find.c; the UNIT-AWARE gap audition
 * + Steps 3-4 (ambient projector + project-into-A) live in
 * src/aic_projection_audit.c (Rule 10 split).
 *
 * approximate_algebras.tex:935-939 (the reduction, verbatim):
 *   "An element P in A is a projection if and only if X = 2P - I is Hermitian and
 *    satisfies the equation X^2 = I. ... An approximate version of these
 *    conditions ... implies that P = 1/4(2I + U + U^dag) is an O(delta+eps)-
 *    projection."
 * We get X = sgn(Y) AMBIENT (Y = s(H - tI) Hermitian, in M_n where sgn is exact)
 * and P_amb = 1/2(I + X). Then P = Phi_tilde(P_amb): the projection back onto A.
 *
 * approximate_algebras.tex:929 (the nontriviality contract, verbatim):
 *   "A delta-projection P is called nontrivial if both P and I-P are nonvanishing."
 * UNIT-AWARE (bead aic-66n). The relevant unit is the ALGEBRA's unit, NOT the
 * ambient 1_n. For an S_P wrapper the unit is Ptilde_m = Co_P(P) (tex:1082;
 * FINDINGS §C7), so "I-P" must be read as "U_A - P" with U_A = Phi_tilde(1_n).
 * The old internal gate tested the AMBIENT ||1_n - P|| ~ 1 — VACUOUS on a wrapper
 * (FINDINGS §C11) — so a wrapper-collapsing split (P ~ Ptilde_m) slipped through.
 * The fix: aic_projection_audit picks the gap by the UNIT-AWARE ||U_A - P|| > c
 * criterion, auditioning (H_k, gap) candidates largest-gap-first.
 *
 * approximate_algebras.tex:917-919, :926-929 (the output contract):
 *   delta-projection: P^dag = P, ||P^2 - P|| <= delta (P^2 = P * P the STAR);
 *   nonvanishing: |||P|| - 1| <= O(delta+eps); nontrivial: P AND U_A-P nonvanishing.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner_internal.h"   /* aic_corner_gamma_opnorm_ub (aic-qgs dodge) */
#include "aic/aic_ecstar.h"
#include "aic/aic_funcalc.h"       /* aic_sgn */
#include "aic/aic_projection.h"
#include "aic_projection_internal.h"

/* Step 3: ambient spectral projector. Y = s(H - tI), s = 1/max(t-lmin, lmax-t);
 * ASSERT ||Y^2 - I||_op < 1 (tex:516,520 funcalc basin); X = sgn(Y); P_amb =
 * (I + X)/2 (exact ambient idempotent). `Pamb` caller-init'd n x n. */
void aic_projection_ambient(acb_mat_t Pamb, const acb_mat_t H, double t,
                            double lmin, double lmax, slong prec)
{
    slong n = acb_mat_nrows(H);
    double denom = fmax(t - lmin, lmax - t);
    assert(denom > 0.0 && "projection: degenerate scaling (t at spectrum edge)");
    double s = 1.0 / denom;

    acb_mat_t Y, Imat, X, Y2, D;
    arb_t sc, e, one;
    acb_mat_init(Y, n, n);
    acb_mat_init(Imat, n, n);
    acb_mat_init(X, n, n);
    acb_mat_init(Y2, n, n);
    acb_mat_init(D, n, n);
    arb_init(sc); arb_init(e); arb_init(one);

    acb_mat_one(Imat);
    arb_set_d(sc, t);
    acb_mat_scalar_mul_arb(Y, Imat, sc, prec);     /* tI */
    acb_mat_sub(Y, H, Y, prec);                    /* H - tI */
    arb_set_d(sc, s);
    acb_mat_scalar_mul_arb(Y, Y, sc, prec);        /* Y = s(H - tI) */

    /* ASSERT the funcalc basin ||Y^2 - I||_op < 1. Use the certified mid+radius
     * upper bound (aic_corner_gamma_opnorm_ub), NOT aic_mat_opnorm (the Gram
     * false-fail dodge, beads aic-qgs/aic-2yo). */
    acb_mat_sqr(Y2, Y, prec);
    acb_mat_sub(D, Y2, Imat, prec);
    aic_corner_gamma_opnorm_ub(e, NULL, D, prec);   /* certified >= ||Y^2-I||_op */
    arb_set_d(one, 1.0);
    if (!arb_lt(e, one)) {
        fprintf(stderr, "aic_projection: ambient sgn basin VIOLATED: "
                "||Y^2 - I||_op not certainly < 1 (gap too small relative to the "
                "spectral radius, or prec too low). s=%.6g, t=%.6g.\n", s, t);
        abort();
    }

    aic_sgn(X, Y, prec);                           /* X = sgn(Y), AMBIENT */
    acb_mat_add(Pamb, Imat, X, prec);
    acb_mat_scalar_mul_2exp_si(Pamb, Pamb, -1);    /* P_amb = (I + sgn Y)/2 */

    arb_clear(one); arb_clear(e); arb_clear(sc);
    acb_mat_clear(D); acb_mat_clear(Y2); acb_mat_clear(X);
    acb_mat_clear(Imat); acb_mat_clear(Y);
}

/* Step 4: project a matrix W into A via the ALGEBRA'S OWN idempotent Phi_tilde,
 * Phi_tilde(W) = W * I (the star with the unit, .tex:2211; FINDINGS §C3 — A is an
 * OBLIQUE image, so the structure-respecting projection is Phi_tilde, NOT the
 * Frobenius projector). `out` caller-init'd n x n. */
void aic_projection_into_A(acb_mat_t out, const aic_ecstar *A, const acb_mat_t W,
                           slong prec)
{
    slong n = A->n;
    acb_mat_t Imat;
    acb_mat_init(Imat, n, n);
    acb_mat_one(Imat);
    aic_ecstar_star(out, A, W, Imat, prec);   /* W * I = Phi_tilde(W) */
    acb_mat_clear(Imat);
}

/* The nontriviality threshold for the unit-aware audition. Chosen >= the
 * cstar_build wrapper backstop (0.15, src/aic_cstar_build.c:258) so a
 * finder-accepted P always passes that downstream gate: a P with ||P|| > 0.3 and
 * ||U_A - P|| > 0.3 has ||Ptilde_m - P|| = ||U_A - P|| > 0.3 >= 0.15. */
#define AIC_PROJECTION_NONTRIV_C 0.3

/* Step 5a: certified star defect delta = ||P*P - P||_op via aic_ecstar_star +
 * aic_corner_gamma_opnorm_ub (the aic-qgs Gram-false-fail dodge). `delta`
 * caller-init'd arb_t. */
static void star_defect(arb_t delta, const aic_ecstar *A, const acb_mat_t P,
                        slong prec)
{
    slong n = A->n;
    acb_mat_t PP, D;
    acb_mat_init(PP, n, n);
    acb_mat_init(D, n, n);
    aic_ecstar_star(PP, A, P, P, prec);         /* P * P = Phi(P.P) (STAR) */
    acb_mat_sub(D, PP, P, prec);                /* P*P - P */
    aic_corner_gamma_opnorm_ub(delta, NULL, D, prec);   /* certified >= ||D||_op */
    acb_mat_clear(D); acb_mat_clear(PP);
}

void aic_projection_nontrivial(acb_mat_t P, arb_t delta, arb_t Pnorm,
                               arb_t ImPnorm, const aic_ecstar *A,
                               aic_projection_witness *wit, slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(acb_mat_nrows(P) == n && acb_mat_ncols(P) == n &&
           "aic_projection_nontrivial: P must be n x n");
    if (d <= 1) {
        fprintf(stderr, "aic_projection: dim A = %ld <= 1. "
                "lem_nontriv_projection requires 1 < dim A; a 1-dim algebra "
                "(scalars) has NO nontrivial projection (only 0 and I).\n",
                (long) d);
        abort();
    }

    /* The ALGEBRA unit U_A = Phi_tilde(1_n) (= 1_n * 1_n, the star with the unit).
     * = Ptilde_m for an S_P wrapper (tex:1082; FINDINGS §C7), ~1_n at top level.
     * This is the unit the nontriviality test (tex:929) is measured against. */
    acb_mat_t U_A, Imat;
    acb_mat_init(U_A, n, n);
    acb_mat_init(Imat, n, n);
    acb_mat_one(Imat);
    aic_ecstar_star(U_A, A, Imat, Imat, prec);   /* U_A = Phi_tilde(1_n) */

    /* Steps 1-4: the unit-aware audition picks the (H_k, gap) and builds P. */
    aic_projection_witness w;
    aic_projection_audit(P, A, U_A, &w, AIC_PROJECTION_NONTRIV_C, prec);

    /* Step 5: certify. delta = ||P*P-P||_op (star); ||P||, ||U_A - P||. The
     * op-norms use the certified mid+radius upper bound (aic_corner_gamma_opnorm_-
     * ub), NOT aic_mat_opnorm: P's Gram has near-zero off-diagonals whose
     * accumulated radii false-fail aic_mat_opnorm's Gram check (aic-qgs/aic-2yo). */
    star_defect(delta, A, P, prec);
    aic_corner_gamma_opnorm_ub(Pnorm, NULL, P, prec);
    {
        acb_mat_t UmP;
        acb_mat_init(UmP, n, n);
        acb_mat_sub(UmP, U_A, P, prec);          /* U_A - P (the ALGEBRA unit) */
        aic_corner_gamma_opnorm_ub(ImPnorm, NULL, UmP, prec);
        acb_mat_clear(UmP);
    }

    /* Nontriviality (fail loud, Rule 4): both ||P|| and ||U_A - P|| bounded away
     * from 0 (tex:929). The audition already enforced this against U_A, so this is
     * an OUTPUT re-assertion (not the intent). A trivial P (near 0 or near U_A)
     * has one of these ~0. NOTE the comparison is against U_A, the ALGEBRA unit,
     * NOT the ambient 1_n (the old vacuous-on-a-wrapper gate; FINDINGS §C11). */
    double pn = arf_get_d(arb_midref(Pnorm), ARF_RND_NEAR);
    double ipn = arf_get_d(arb_midref(ImPnorm), ARF_RND_NEAR);
    if (pn <= AIC_PROJECTION_NONTRIV_C || ipn <= AIC_PROJECTION_NONTRIV_C) {
        fprintf(stderr, "aic_projection: TRIVIAL projection produced "
                "(||P||=%.4f, ||U_A-P||=%.4f; one is near 0 so P ~ 0 or P ~ the "
                "algebra unit U_A). The unit-aware audition should have prevented "
                "this — investigate (degenerate spectrum / mis-chosen gap).\n",
                pn, ipn);
        abort();
    }

    if (wit) *wit = w;

    acb_mat_clear(Imat);
    acb_mat_clear(U_A);
}
