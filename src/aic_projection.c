/* aic_projection.c — the nontrivial-projection finder (bead aic-mqf), orchestrator
 * + Steps 3-5 of Route A-ambient. Realizes lem_nontriv_projection
 * (approximate_algebras.tex:931-932). See include/aic_projection.h for the full
 * route, the rem_X2-safety argument, and the audition note (Law 4). Steps 1-2
 * (H-search + gap) live in src/aic_projection_find.c.
 *
 * approximate_algebras.tex:935-939 (the reduction, verbatim):
 *   "An element P in A is a projection if and only if X = 2P - I is Hermitian and
 *    satisfies the equation X^2 = I. ... An approximate version of these
 *    conditions ... implies that P = 1/4(2I + U + U^dag) is an O(delta+eps)-
 *    projection."
 * We get X = sgn(Y) AMBIENT (Y = s(H - tI) Hermitian, in M_n where sgn is exact)
 * and P_amb = 1/2(I + X). Then P = Pi_A(P_amb): the projection back onto A.
 *
 * approximate_algebras.tex:516, :520 (the funcalc basin we ASSERT before sgn):
 *   "|X| = (X^2)^{1/2}, sgn(X) = X (X^2)^{-1/2}  ( ||X^2 - x0^2 I|| < x0^2 )";
 *   "||sgn(X) - X|| <= ||X|| O(||X^2 - I||)  if ||X^2 - I|| < 1."
 * With Y's eigenvalues in [-1,1] and |eigval| >= s*g/2 off 0, ||Y^2 - I||_op < 1
 * holds; we certify it and abort otherwise (Rule 4).
 *
 * approximate_algebras.tex:917-919, :926-929 (the output contract):
 *   delta-projection: P^dag = P, ||P^2 - P|| <= delta (P^2 = P * P the STAR);
 *   nonvanishing: |||P|| - 1| <= O(delta+eps); nontrivial: P AND I-P nonvanishing.
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
#include "aic/aic_funcalc.h"
#include "aic/aic_mat.h"
#include "aic/aic_projection.h"
#include "aic_projection_internal.h"

/* Step 3: ambient spectral projector. Y = s(H - tI), s = 1/max(t-lam_min,
 * lam_max-t); ASSERT ||Y^2 - I||_op < 1; X = sgn(Y); P_amb = (I + X)/2 (exact
 * ambient idempotent). `Pamb` caller-init'd n x n. */
static void ambient_projector(acb_mat_t Pamb, const acb_mat_t H, double t,
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

    /* ASSERT the funcalc basin ||Y^2 - I||_op < 1 (tex:516,520). Use the certified
     * mid+radius upper bound (aic_corner_gamma_opnorm_ub), NOT aic_mat_opnorm: the
     * Gram (Y^2-I)^dag(Y^2-I) has near-zero off-diagonals whose accumulated radii
     * false-fail aic_mat_opnorm's Gram Hermiticity check (beads aic-qgs/aic-2yo). */
    acb_mat_sqr(Y2, Y, prec);
    acb_mat_sub(D, Y2, Imat, prec);
    aic_corner_gamma_opnorm_ub(e, NULL, D, prec);   /* certified >= ||Y^2-I||_op */
    arb_set_d(one, 1.0);
    if (!arb_lt(e, one)) {
        fprintf(stderr, "aic_projection: ambient sgn basin VIOLATED: "
                "||Y^2 - I||_op not certainly < 1 (gap too small relative to the "
                "spectral radius, or prec too low). Y = s(H - tI), s=%.6g, "
                "t=%.6g.\n", s, t);
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
 * NOT the Frobenius-orthogonal projector. A = Img Phi_tilde is an OBLIQUE image:
 * Phi_tilde is Hermicity-preserving but in general NOT HS-self-adjoint, so the
 * Frobenius-orthogonal projector Pi_A(W) = sum_k <B_k,W>_F B_k does NOT respect
 * A's star structure. Empirically (bead aic-mqf) the Frobenius route leaves the
 * star defect ||P*P-P|| = O(1) (constant in eta ~ 0.5), whereas Phi_tilde(P_amb)
 * gives the O(eta) defect the lemma promises (C ~ 0.04, flat).
 *
 * Phi_tilde IS available through the public ecstar API as the star with the unit:
 *   Phi_tilde(W) = W * I = Phi_tilde(W . I)   (.tex:2211 "X*I = X = I*X" extends to
 * X*I = Phi_tilde(X) for ANY X). Verified == the superoperator application to
 * ~6e-72 (bead aic-mqf). So this step uses ONLY aic_ecstar_star, no assoc
 * internals. `out` caller-init'd n x n. */
static void project_into_A(acb_mat_t out, const aic_ecstar *A, const acb_mat_t W,
                           slong prec)
{
    slong n = A->n;
    acb_mat_t Imat;
    acb_mat_init(Imat, n, n);
    acb_mat_one(Imat);
    aic_ecstar_star(out, A, W, Imat, prec);   /* W * I = Phi_tilde(W) */
    acb_mat_clear(Imat);
}

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

    /* Steps 1-2: pick H, find the gap + threshold. */
    acb_mat_t H;
    acb_mat_init(H, n, n);
    slong k_chosen = -1;
    aic_projection_pick_H(H, &k_chosen, A, prec);
    double t = 0.0, gap = 0.0, lmin = 0.0, lmax = 0.0;
    slong m = 0;
    aic_projection_gap(&t, &gap, &lmin, &lmax, &m, H, prec);

    /* Step 3: ambient projector P_amb (exact idempotent in M_n). */
    acb_mat_t Pamb;
    acb_mat_init(Pamb, n, n);
    ambient_projector(Pamb, H, t, lmin, lmax, prec);

    /* Step 4: project into A. */
    project_into_A(P, A, Pamb, prec);

    /* Step 5: certify. delta = ||P*P-P||_op (star); ||P||, ||I-P||; membership.
     * The op-norms use the certified mid+radius upper bound (aic_corner_gamma_-
     * opnorm_ub), NOT aic_mat_opnorm: P's Gram has near-zero off-diagonals whose
     * accumulated radii false-fail aic_mat_opnorm's Gram check (aic-qgs/aic-2yo).
     * The +||rad||_F inflation is ~1e-72, negligible for the nontriviality gates. */
    star_defect(delta, A, P, prec);
    aic_corner_gamma_opnorm_ub(Pnorm, NULL, P, prec);
    {
        acb_mat_t Imat, ImP;
        acb_mat_init(Imat, n, n);
        acb_mat_init(ImP, n, n);
        acb_mat_one(Imat);
        acb_mat_sub(ImP, Imat, P, prec);
        aic_corner_gamma_opnorm_ub(ImPnorm, NULL, ImP, prec);
        acb_mat_clear(ImP); acb_mat_clear(Imat);
    }

    /* Nontriviality (fail loud): both ||P|| and ||I-P|| bounded away from 0. A
     * trivial P (near 0 or near I) has one of these ~0; the interior gap (m in
     * [1,n-1]) should prevent that, but we ASSERT the output, not the intent. */
    double pn = arf_get_d(arb_midref(Pnorm), ARF_RND_NEAR);
    double ipn = arf_get_d(arb_midref(ImPnorm), ARF_RND_NEAR);
    if (pn <= 0.3 || ipn <= 0.3) {
        fprintf(stderr, "aic_projection: TRIVIAL projection produced "
                "(||P||=%.4f, ||I-P||=%.4f; one is near 0 so P ~ 0 or P ~ I). "
                "The interior gap (m=%ld of n=%ld) should have prevented this — "
                "investigate (degenerate spectrum / mis-chosen threshold).\n",
                pn, ipn, (long) m, (long) n);
        abort();
    }

    if (wit) {
        wit->k_chosen = k_chosen;
        wit->lam_min = lmin;
        wit->lam_max = lmax;
        wit->t = t;
        wit->gap = gap;
        wit->m = m;
    }

    acb_mat_clear(Pamb);
    acb_mat_clear(H);
}
