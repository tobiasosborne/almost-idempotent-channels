/* aic_corner_hilbert.c — Increment 2b of the corner module: the ONE-DIMENSIONAL-Q
 * machinery of section 7. The lem_PQ_Hilb inner product (the Hilbert structure on
 * S_{P,Q} when Q is one-dimensional), lem_PQR (norm multiplicativity), lem_1d_proj
 * and the 1d-equivalence predicate. The Ha-map (the implicit-equation solve) lives
 * in the sibling src/aic_corner_ha.c (split for the ~200 LOC limit, Rule 10).
 * Builds on Increment 1+2a (aic_corner_Co, aic_corner_apply, aic_corner_extract,
 * aic_corner_cdot, aic_corner_Ptilde). See include/aic_corner.h for the contracts.
 *
 * approximate_algebras.tex:1124-1132 — lem_PQ_Hilb (verbatim, the inner product):
 *   "If Q is one-dimensional, then S_{P,Q} is a Hilbert space with the inner
 *    product <.|.> defined by the equation
 *        Y^dag . X = <Y|X> Qtilde      (X, Y in S_{P,Q}),                (1dQ_ip)
 *    where Y^dag . X = Co_Q(Y^dag X) and Qtilde = Co_Q(Q). Furthermore,
 *        | <X|X> - ||X||^2 | <= O(delta+eps) ||X||^2.                    (1dQ_ip_norm)"
 *
 * THE STAR, NOT THE PLAIN PRODUCT (CLAUDE.md callout, LOAD-BEARING). "Y^dag X" in
 * the paper's eqs is the product in A, i.e. the Choi-Effros STAR Y^dag * X =
 * Phi_tilde(Y^dag X), aic_ecstar_star (NOT acb_mat_mul). So the inner-product
 * metric routes Y^dag * X through the star, then through Co_Q. Dropping the star
 * (plain Y^dag X) is RED on the oblique compress_idemp(4,2) fixture, where the
 * star is load-bearing (tests/test_corner.c T10, non-vacuity gap measured O(1)).
 *
 * WHY THE SCALAR IS <Qtilde,W>_F / <Qtilde,Qtilde>_F. For one-dimensional Q the
 * corner S_Q = S_{Q,Q} = C Qtilde is the COMPLEX line spanned by Qtilde (.tex:1160
 * "S_{Q,Q} = C Qtilde ~ C"). W = Co_Q(Y^dag * X) lies in S_Q, so W = s Qtilde for a
 * unique complex s = <Y|X>; projecting onto the line gives
 * s = <Qtilde, W>_F / <Qtilde, Qtilde>_F (the Frobenius projection coefficient).
 *
 * lem_PQR (.tex:1162-1177): for 1d Q, | ||X.Y|| - ||X|| ||Y|| | <= O(delta+eps)
 * ||X|| ||Y|| (X in S_{P,Q}, Y in S_{Q,R}). The proof uses (XY)^dag(XY) =
 * ||X||^2 Y^dag Y + O(delta+eps), i.e. X^dag X ~ ||X||^2 Qtilde via the 1dQ inner
 * product. We test the inequality directly on the extracted bases (a norm test).
 *
 * lem_1d_proj (.tex:1179-1185): for 1d P AND 1d Q, dim S_{P,Q} <= 1. The
 * equivalence P ~ Q := dim S_{P,Q} = 1 (.tex:1187) is reflexive, symmetric, and
 * transitive (via lem_PQR). aic_corner_equiv_1d/dim_S realize the predicate; the
 * transitivity is a TEST property, not a routine.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_corner.h"
#include "aic_corner_internal.h"  /* aic_corner_gamma_opnorm_ub (aic-2yo workaround) */
#include "aic/aic_ecstar.h"

/* <A, B>_F = Tr(A^dag B) = sum_{a,b} conj(A[a,b]) B[a,b] into `out` (caller-init'd
 * acb_t). A, B are n x n. Kept LOCAL to this TU (the production aic_corner_frob_ip
 * is in the corner_internal split; the Ha sibling has its own copy too) so the
 * Hilbert results are self-contained at the ~200 LOC boundary. */
static void frob_ip(acb_t out, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong n = acb_mat_nrows(B);
    acb_t t;
    acb_init(t);
    acb_zero(out);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++) {
            acb_conj(t, acb_mat_entry(A, a, b));
            acb_mul(t, t, acb_mat_entry(B, a, b), prec);
            acb_add(out, out, t, prec);
        }
    acb_clear(t);
}

void aic_corner_ip_1d(acb_t out, const aic_ecstar *A, const acb_mat_t Qtilde,
                      const acb_mat_t CoQ, const acb_mat_t X, const acb_mat_t Y,
                      slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(acb_mat_nrows(CoQ) == d && acb_mat_ncols(CoQ) == d &&
           "ip_1d: CoQ must be dim_A x dim_A");
    assert(acb_mat_nrows(Qtilde) == n && "ip_1d: Qtilde must be n x n");

    /* W = Co_Q(Y^dag * X): the STAR Y^dag * X (.tex:1128 product-in-A), then Co_Q. */
    acb_mat_t Yd, prod, W;
    acb_mat_init(Yd, n, n);
    acb_mat_init(prod, n, n);
    acb_mat_init(W, n, n);
    acb_mat_conjugate_transpose(Yd, Y);          /* Y^dag in S_{Q,P} */
    aic_ecstar_star(prod, A, Yd, X, prec);       /* Y^dag * X (Choi-Effros star) */
    aic_corner_apply(W, A, CoQ, prod, prec);     /* Co_Q(Y^dag * X) in S_Q */

    /* <Y|X> = <Qtilde, W>_F / <Qtilde, Qtilde>_F (project onto the line C Qtilde). */
    acb_t num, den;
    acb_init(num);
    acb_init(den);
    frob_ip(num, Qtilde, W, prec);
    frob_ip(den, Qtilde, Qtilde, prec);
    /* Fail loud (Rule 4): a 1d delta-projection has ||Qtilde||_F ~ 1, so a
     * near-zero denominator means Q is not a genuine 1d projection in A. */
    {
        arb_t mag, floor;
        arb_init(mag);
        arb_init(floor);
        acb_abs(mag, den, prec);
        arb_set_d(floor, 1e-9);
        assert(arb_gt(mag, floor) &&
               "ip_1d: <Qtilde,Qtilde>_F ~ 0 -- Q not a one-dimensional "
               "delta-projection in A (Qtilde = Co_Q(Q) degenerate)");
        arb_clear(floor);
        arb_clear(mag);
    }
    acb_div(out, num, den, prec);

    acb_clear(den);
    acb_clear(num);
    acb_mat_clear(W);
    acb_mat_clear(prod);
    acb_mat_clear(Yd);
}

slong aic_corner_dim_S(const aic_ecstar *A, const acb_mat_t P, const acb_mat_t Q,
                       slong prec)
{
    slong d = A->dim_A;
    acb_mat_t Co;
    acb_mat_init(Co, d, d);
    aic_corner_Co(Co, A, P, Q, prec);
    /* Reuse aic_corner_extract's tested round(trace)-vs-SVD-gap cross-check
     * (fail loud on mismatch, Rule 4), then discard the basis. */
    acb_mat_t *C;
    slong dim_S;
    aic_corner_extract(&C, &dim_S, Co, A, prec);
    for (slong m = 0; m < dim_S; m++) acb_mat_clear(C[m]);
    if (C) flint_free(C);
    acb_mat_clear(Co);
    return dim_S;
}

int aic_corner_equiv_1d(const aic_ecstar *A, const acb_mat_t P, const acb_mat_t Q,
                        slong prec)
{
    return aic_corner_dim_S(A, P, Q, prec) == 1 ? 1 : 0;
}

void aic_corner_pqr_defect(arb_t out, const aic_ecstar *A, const acb_mat_t P,
                           const acb_mat_t Q, const acb_mat_t R, slong prec)
{
    slong n = A->n, d = A->dim_A;
    /* 1d-Q precondition (lem_PQR hypothesis "Q is one-dimensional"). Fail loud. */
    assert(aic_corner_dim_S(A, Q, Q, prec) == 1 &&
           "pqr_defect: Q not one-dimensional (dim S_Q != 1)");

    acb_mat_t CoPQ, CoQR, CoPR;
    acb_mat_init(CoPQ, d, d);
    acb_mat_init(CoQR, d, d);
    acb_mat_init(CoPR, d, d);
    aic_corner_Co(CoPQ, A, P, Q, prec);          /* extract S_{P,Q} */
    aic_corner_Co(CoQR, A, Q, R, prec);          /* extract S_{Q,R} */
    aic_corner_Co(CoPR, A, P, R, prec);          /* target of X . Y */

    acb_mat_t *CX, *CY;
    slong dPQ, dQR;
    aic_corner_extract(&CX, &dPQ, CoPQ, A, prec);
    aic_corner_extract(&CY, &dQR, CoQR, A, prec);

    arb_t worst, nx, ny, nxy, t;
    arb_init(worst);
    arb_init(nx);
    arb_init(ny);
    arb_init(nxy);
    arb_init(t);
    arb_zero(worst);
    acb_mat_t prod;
    acb_mat_init(prod, n, n);
    for (slong i = 0; i < dPQ; i++)
        for (slong j = 0; j < dQR; j++) {
            /* CERTIFIED mid+radius upper bounds, NOT aic_mat_opnorm's Gram path:
             * the cdot output X.Y carries a tiny accumulated arb radius (~1e-70)
             * whose Gram (X.Y)^dag(X.Y) off-diagonals false-fail the relative-
             * Hermiticity check (bead aic-2yo) and SIGABRT on the dim S_{P,Q}=2,
             * eta>0 mixconj corner -- even though X.Y is a genuine O(1) operator
             * (||X.Y|| ~ 1). aic_corner_gamma_opnorm_ub sidesteps the Gram path;
             * the ~1e-70 inflation is negligible vs the O(delta+eps) defect. */
            aic_corner_gamma_opnorm_ub(nx, NULL, CX[i], prec);
            aic_corner_gamma_opnorm_ub(ny, NULL, CY[j], prec);
            aic_corner_cdot(prod, A, CoPR, CX[i], CY[j], prec);   /* X . Y */
            aic_corner_gamma_opnorm_ub(nxy, NULL, prod, prec);
            /* | ||X.Y|| - ||X|| ||Y|| | */
            arb_mul(t, nx, ny, prec);
            arb_sub(t, nxy, t, prec);
            arb_abs(t, t);
            if (arb_gt(t, worst)) arb_set(worst, t);
        }
    arb_set(out, worst);

    acb_mat_clear(prod);
    arb_clear(t);
    arb_clear(nxy);
    arb_clear(ny);
    arb_clear(nx);
    arb_clear(worst);
    for (slong m = 0; m < dQR; m++) acb_mat_clear(CY[m]);
    if (CY) flint_free(CY);
    for (slong m = 0; m < dPQ; m++) acb_mat_clear(CX[m]);
    if (CX) flint_free(CX);
    acb_mat_clear(CoPR);
    acb_mat_clear(CoQR);
    acb_mat_clear(CoPQ);
}
