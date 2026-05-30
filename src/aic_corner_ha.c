/* aic_corner_ha.c — Increment 2b of the corner module: the Ha-map (the implicit-
 * equation solve of section 7, .tex:1146-1160). The most convention-sensitive
 * routine in the module: a chain of nested COMPRESSED products (the star-then-Co
 * cdot) whose Co indices must each match the target corner, decided by the
 * subspace-membership bookkeeping in include/aic_corner.h. The Hilbert-space
 * results (lem_PQ_Hilb inner product, lem_PQR, lem_1d_proj) live in the sibling
 * src/aic_corner_hilbert.c. Builds on aic_corner_Co/apply/extract, aic_corner_cdot,
 * aic_corner_Ptilde, aic_corner_ip_1d.
 *
 * approximate_algebras.tex:1146-1153 — the Ha-map (verbatim, the definition):
 *   "we define a map Ha^Q_{P,R} : S_{P,R} -> B(S_{R,Q}, S_{P,Q}). ... For each
 *    Z in S_{P,R} and X in S_{R,Q}, the element Ha^Q_{P,R}(Z)(X) in S_{P,Q} is
 *    defined by the condition
 *        (Y^dag . Z) . X + Y^dag . (Z . X) = 2 <Y | Ha^Q_{P,R}(Z)(X)> Qtilde
 *        for all Y in S_{P,Q}.                                          (Ha_def)
 *    This symmetric definition enjoys the exact equality
 *        Ha^Q_{R,P}(Z^dag) = Ha^Q_{P,R}(Z)^dag.                         (Ha_dag)
 *    On the other hand, ||Ha^Q_{P,R}(Z)(X) - Z . X||_Euc <= O(delta+eps)
 *    ||Z|| ||X||_Euc."
 *
 * CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). Ha_def is an IMPLICIT linear equation: it
 * pins Ha(Z)(X) by its inner products against every Y in S_{P,Q}. In finite dim
 * with {C^{PQ}_l} the extracted Frobenius-orthonormal basis of S_{P,Q}, write
 * Ha(Z)(X) = sum_m a_m C^{PQ}_m. Taking Y = C^{PQ}_l in Ha_def and reading off the
 * scalar coefficient of Qtilde on both sides (S_Q = C Qtilde is one-dimensional):
 *      scalar(LHS(C^{PQ}_l)) = 2 sum_m <C^{PQ}_l | C^{PQ}_m> a_m,
 * i.e. b_l := (1/2) scalar(LHS(C^{PQ}_l)) = sum_m G_{lm} a_m, the Gram system
 * G a = b with G_{lm} = <C^{PQ}_l | C^{PQ}_m> the lem_PQ_Hilb inner product (NOT
 * the Frobenius ip). Since |<C|C> - ||C||^2| <= O(delta+eps) and the {C^{PQ}_l}
 * are Frobenius-orthonormal, G = I + O(delta+eps) is well-conditioned; we ASSERT
 * ||G - I|| < 1 (fail loud, Rule 4) and solve G a = b via acb_mat_solve (certified).
 * scalar(W in S_Q) = <Qtilde, W>_F / <Qtilde, Qtilde>_F (the line-projection coeff,
 * same as aic_corner_ip_1d).
 *
 * THE ||G-I||<1 GUARD USES THE CERTIFIED MIDPOINT+RADIUS UPPER BOUND (bead
 * aic-2yo), NOT aic_mat_opnorm: that Gram path's off-diagonal Hermiticity check
 * FALSE-FAILS on the near-zero off-diagonals of an almost-identity G and SIGABRTs
 * on a genuine dim S_{P,Q}=2, eta>0 corner trivially in-basin (||G-I|| ~ 1.3e-5).
 * aic_corner_gamma_opnorm_ub (the alpha module's same-issue workaround) returns a
 * certified ||mid(G-I)||_op + ||rad(G-I)||_F >= ||G-I||_op; the <1 decision is
 * unchanged.
 *
 * THE Co-INDEX BOOKKEEPING (the highest convention-risk in the module). Each cdot
 * computes Co_TARGET(left * right); the TARGET corner is fixed by the operands'
 * membership (see include/aic_corner.h):
 *   Y^dag . Z      : S_{Q,P} x S_{P,R} -> S_{Q,R}   cdot(Co_{Q,R}, Y^dag, Z)
 *   (Y^dag.Z) . X  : S_{Q,R} x S_{R,Q} -> S_Q       cdot(Co_Q,   Y^dag.Z, X)
 *   Z . X          : S_{P,R} x S_{R,Q} -> S_{P,Q}   cdot(Co_{P,Q}, Z, X)
 *   Y^dag . (Z.X)  : S_{Q,P} x S_{P,Q} -> S_Q       cdot(Co_Q,   Y^dag, Z.X)
 * Y^dag = (C^{PQ}_l)^dag via acb_mat_conjugate_transpose (lands in S_{Q,P}). A
 * WRONG Co in this chain (e.g. Co_{P,Q} where Co_{Q,R} is required) corrupts the
 * answer; tests/test_corner.c T11 mutation-proves a wrong-Co swap RED. Z . X does
 * not depend on the test index l, so it is computed once per input column.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner.h"
#include "aic_corner_internal.h"  /* aic_corner_gamma_opnorm_ub (aic-2yo workaround) */
#include "aic_ecstar.h"

/* scalar(W in S_Q) = <Qtilde, W>_F / den, den = <Qtilde, Qtilde>_F (precomputed,
 * caller asserts den != 0). W, Qtilde n x n; `out` caller-init'd acb_t. */
static void scalar_in_SQ(acb_t out, const acb_mat_t Qtilde, const acb_t den,
                         const acb_mat_t W, slong prec)
{
    slong n = acb_mat_nrows(W);
    acb_t num, t;
    acb_init(num);
    acb_init(t);
    acb_zero(num);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++) {
            acb_conj(t, acb_mat_entry(Qtilde, a, b));
            acb_mul(t, t, acb_mat_entry(W, a, b), prec);
            acb_add(num, num, t, prec);
        }
    acb_div(out, num, den, prec);
    acb_clear(t);
    acb_clear(num);
}

void aic_corner_ha(acb_mat_t Ha, const aic_ecstar *A, const acb_mat_t Z,
                   const acb_mat_t P, const acb_mat_t R, const acb_mat_t Q,
                   slong prec)
{
    slong n = A->n, d = A->dim_A;

    /* 1d-Q precondition (.tex:1146 "Q a one-dimensional delta-projection"). */
    assert(aic_corner_dim_S(A, Q, Q, prec) == 1 &&
           "ha: Q not one-dimensional (dim S_Q != 1)");

    /* Co's for the cdot chain + the two extracted corner bases. */
    acb_mat_t CoPQ, CoRQ, CoQR, CoQ;
    acb_mat_init(CoPQ, d, d);
    acb_mat_init(CoRQ, d, d);
    acb_mat_init(CoQR, d, d);
    acb_mat_init(CoQ, d, d);
    aic_corner_Co(CoPQ, A, P, Q, prec);     /* rows of Ha: S_{P,Q}; Z.X target  */
    aic_corner_Co(CoRQ, A, R, Q, prec);     /* cols of Ha: S_{R,Q} (inputs X)   */
    aic_corner_Co(CoQR, A, Q, R, prec);     /* Y^dag . Z target                 */
    aic_corner_Co(CoQ,  A, Q, Q, prec);     /* S_Q target (both scalar terms)   */

    acb_mat_t *CPQ, *CRQ;
    slong dPQ, dRQ;
    aic_corner_extract(&CPQ, &dPQ, CoPQ, A, prec);
    aic_corner_extract(&CRQ, &dRQ, CoRQ, A, prec);
    assert(acb_mat_nrows(Ha) == dPQ && acb_mat_ncols(Ha) == dRQ &&
           "ha: Ha must be dim S_{P,Q} x dim S_{R,Q}");

    acb_mat_t Qtilde;
    acb_mat_init(Qtilde, n, n);
    aic_corner_Ptilde(Qtilde, A, Q, prec);  /* Qtilde = Co_Q(Q), the S_Q generator */
    acb_t den;
    acb_init(den);
    {   /* den = <Qtilde,Qtilde>_F, asserted bounded away from 0 (Rule 4). */
        acb_t t;
        acb_init(t);
        acb_zero(den);
        for (slong a = 0; a < n; a++)
            for (slong b = 0; b < n; b++) {
                acb_conj(t, acb_mat_entry(Qtilde, a, b));
                acb_mul(t, t, acb_mat_entry(Qtilde, a, b), prec);
                acb_add(den, den, t, prec);
            }
        acb_clear(t);
        arb_t mag, floor;
        arb_init(mag);
        arb_init(floor);
        acb_abs(mag, den, prec);
        arb_set_d(floor, 1e-9);
        assert(arb_gt(mag, floor) && "ha: <Qtilde,Qtilde>_F ~ 0 (Q degenerate)");
        arb_clear(floor);
        arb_clear(mag);
    }

    /* Gram G_{lm} = <C^{PQ}_l | C^{PQ}_m> (lem_PQ_Hilb ip); assert ||G-I||<1. */
    acb_mat_t G;
    acb_mat_init(G, dPQ, dPQ);
    for (slong l = 0; l < dPQ; l++)
        for (slong m = 0; m < dPQ; m++)
            aic_corner_ip_1d(acb_mat_entry(G, l, m), A, Qtilde, CoQ,
                             CPQ[m], CPQ[l], prec);   /* ip_1d(X=C_m, Y=C_l) */
    {   /* ||G - I||_op < 1 (well-conditioned; fail loud, Rule 4). Certified
         * mid+radius upper bound, NOT aic_mat_opnorm's Gram path (aic-2yo; see
         * the file docstring). */
        acb_mat_t GmI;
        arb_t e, one;
        acb_mat_init(GmI, dPQ, dPQ);
        arb_init(e);
        arb_init(one);
        acb_mat_one(GmI);
        acb_mat_sub(GmI, G, GmI, prec);
        aic_corner_gamma_opnorm_ub(e, NULL, GmI, prec);
        arb_one(one);
        assert(arb_lt(e, one) &&
               "ha: ||G - I||_op >= 1 -- lem_PQ_Hilb Gram ill-conditioned "
               "(delta+eps too large, or Q not one-dimensional)");
        arb_clear(one);
        arb_clear(e);
        acb_mat_clear(GmI);
    }

    /* Per input column m (X = C^{RQ}_m): b_l = (1/2)[scalar((C_l^dag . Z) . X) +
     * scalar(C_l^dag . (Z . X))]; solve G a = b; store a as column m of Ha. */
    acb_mat_t b, a, Yd, YdZ, term1m, ZX, term2m;
    acb_mat_init(b, dPQ, 1);
    acb_mat_init(a, dPQ, 1);
    acb_mat_init(Yd, n, n);
    acb_mat_init(YdZ, n, n);
    acb_mat_init(term1m, n, n);
    acb_mat_init(ZX, n, n);
    acb_mat_init(term2m, n, n);
    acb_t s1, s2, half;
    acb_init(s1);
    acb_init(s2);
    acb_init(half);
    acb_set_d(half, 0.5);
    for (slong m = 0; m < dRQ; m++) {
        aic_corner_cdot(ZX, A, CoPQ, Z, CRQ[m], prec);      /* Z . X in S_{P,Q} */
        for (slong l = 0; l < dPQ; l++) {
            acb_mat_conjugate_transpose(Yd, CPQ[l]);        /* Y^dag in S_{Q,P} */
            aic_corner_cdot(YdZ, A, CoQR, Yd, Z, prec);     /* Y^dag . Z in S_{Q,R} */
            aic_corner_cdot(term1m, A, CoQ, YdZ, CRQ[m], prec); /* (Y^dag.Z).X in S_Q */
            aic_corner_cdot(term2m, A, CoQ, Yd, ZX, prec);  /* Y^dag.(Z.X) in S_Q */
            scalar_in_SQ(s1, Qtilde, den, term1m, prec);
            scalar_in_SQ(s2, Qtilde, den, term2m, prec);
            acb_add(s1, s1, s2, prec);
            acb_mul(acb_mat_entry(b, l, 0), s1, half, prec);
        }
        int ok = acb_mat_solve(a, G, b, prec);
        assert(ok && "ha: acb_mat_solve(G a = b) failed despite ||G-I||<1");
        for (slong l = 0; l < dPQ; l++)
            acb_set(acb_mat_entry(Ha, l, m), acb_mat_entry(a, l, 0));
    }

    acb_clear(half);
    acb_clear(s2);
    acb_clear(s1);
    acb_mat_clear(term2m);
    acb_mat_clear(ZX);
    acb_mat_clear(term1m);
    acb_mat_clear(YdZ);
    acb_mat_clear(Yd);
    acb_mat_clear(a);
    acb_mat_clear(b);
    acb_mat_clear(G);
    acb_clear(den);
    acb_mat_clear(Qtilde);
    for (slong m = 0; m < dRQ; m++) acb_mat_clear(CRQ[m]);
    if (CRQ) flint_free(CRQ);
    for (slong m = 0; m < dPQ; m++) acb_mat_clear(CPQ[m]);
    if (CPQ) flint_free(CPQ);
    acb_mat_clear(CoQ);
    acb_mat_clear(CoQR);
    acb_mat_clear(CoRQ);
    acb_mat_clear(CoPQ);
}
