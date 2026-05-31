/* aic_ecstar_defect.c — certified-arb axiom-defect estimators for the eps-C*
 * algebra data model (bead aic-knm, module "ecstar"): submult (ax_prodnorm),
 * C* (ax_C*), and unit (ax_eps_unit). Each is a BASIS SWEEP over the stored
 * Frobenius-orthonormal basis {B_k} using the operator norm aic_mat_opnorm,
 * returning a certified arb ball. The trilinear associator (ax_assoc) lives in
 * aic_ecstar_assoc.c, and the involution residual (ax_*) — which routes through
 * a shared generic-apply core so a non-HP teeth can exercise the production code
 * path — lives in aic_ecstar_involution.c (both Rule 10 splits).
 *
 * SCOPE (route decision, bead aic-knm). These are EXACT-ZERO detectors for the
 * eta=0 oracle and basis-sweep LOWER bounds on the true sup-over-unit-ball eps.
 * The involution defect (ax_*) and the unit defect (ax_eps_unit) are ALWAYS-zero
 * invariants for a Hermicity-preserving unital Phi with A = Img Phi (see per-
 * function notes); submult and C* are lower bounds that read ~0 for eta=0. The
 * faithful worst-case HOPM search and the certified SDP upper bound are later
 * cycles (beads aic-0at). See include/aic_ecstar.h and ALGORITHM.md.
 *
 * Convention-safe: every product is the Choi-Effros star X*Y = Phi(XY) via
 * aic_ecstar_star (NOT plain XY — that leaves Img Phi, CLAUDE.md callout 2); the
 * involution is the matrix adjoint; the norm is the operator norm. The n^2 x n^2
 * superoperator is never formed.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_ecstar.h"
#include "aic_ecstar_internal.h"

/* approximate_algebras.tex:410-411 ax_prodnorm:
 *   ||XY|| <= (1+eps) ||X|| ||Y||   (X,Y in A).
 * Basis-sweep LOWER bound: max over pairs (j,k) of the positive part
 *   max(0, ||B_j*B_k||_op - ||B_j||_op ||B_k||_op).
 * =0 for a genuine (submultiplicative) C* algebra. Reads ~0 for eta=0. */
void aic_ecstar_defect_submult(arb_t out, const aic_ecstar *A, slong prec)
{
    slong n = A->n, d = A->dim_A;
    acb_mat_t prod;
    arb_t nj, nk, np, slack;
    acb_mat_init(prod, n, n);
    arb_init(nj);
    arb_init(nk);
    arb_init(np);
    arb_init(slack);
    arb_zero(out);
    for (slong j = 0; j < d; j++) {
        aic_mat_opnorm(nj, A->B[j], prec);
        for (slong k = 0; k < d; k++) {
            aic_mat_opnorm(nk, A->B[k], prec);
            aic_ecstar_star(prod, A, A->B[j], A->B[k], prec);
            aic_mat_opnorm(np, prod, prec);
            arb_mul(slack, nj, nk, prec);     /* ||B_j||_op ||B_k||_op       */
            arb_sub(slack, np, slack, prec);  /* ||B_j*B_k||_op - that       */
            if (arb_gt(slack, out)) arb_set(out, slack); /* positive part    */
        }
    }
    /* out starts at 0 and only grows when slack > out >= 0, so it is the
     * positive part max(0, .) by construction (no separate clamp needed). */
    arb_clear(slack);
    arb_clear(np);
    arb_clear(nk);
    arb_clear(nj);
    acb_mat_clear(prod);
}

/* approximate_algebras.tex:427-428 ax_C*:
 *   ||X^dag X|| >= (1-eps) ||X||^2   (X in A).
 * Basis-sweep LOWER bound: max over k of max(0, 1 - ||B_k^dag * B_k||_op /
 * ||B_k||_op^2). =0 for a genuine C* algebra. FAILS LOUD if any ||B_k||_op is
 * below the floor 1/(2 sqrt n) (a Frobenius-unit op has ||.||_op >= 1/sqrt n;
 * half that is a safe prec-aware floor). Reads ~0 for eta=0. */
void aic_ecstar_defect_cstar(arb_t out, const aic_ecstar *A, slong prec)
{
    slong n = A->n, d = A->dim_A;
    acb_mat_t Bd, prod;
    arb_t nk, npd, floor_tol, slack, sqn;
    acb_mat_init(Bd, n, n);
    acb_mat_init(prod, n, n);
    arb_init(nk);
    arb_init(npd);
    arb_init(floor_tol);
    arb_init(slack);
    arb_init(sqn);
    /* floor = 1/(2 sqrt n): below this an "operator" is numerically zero and the
     * ratio is meaningless. Fail loud (Rule 4). */
    arb_set_si(sqn, n);
    arb_sqrt(sqn, sqn, prec);
    arb_set_d(floor_tol, 0.5);
    arb_div(floor_tol, floor_tol, sqn, prec);
    arb_zero(out);
    for (slong k = 0; k < d; k++) {
        aic_mat_opnorm(nk, A->B[k], prec);
        assert(arb_gt(nk, floor_tol) &&
               "aic_ecstar_defect_cstar: ||B_k||_op below floor (degenerate basis)");
        acb_mat_conjugate_transpose(Bd, A->B[k]);
        aic_ecstar_star(prod, A, Bd, A->B[k], prec);  /* B_k^dag * B_k = Phi(..) */
        aic_mat_opnorm(npd, prod, prec);              /* ||B_k^dag * B_k||_op    */
        arb_mul(slack, nk, nk, prec);                 /* ||B_k||_op^2            */
        arb_div(slack, npd, slack, prec);             /* ratio                  */
        arb_sub_si(slack, slack, 1, prec);            /* ratio - 1              */
        arb_neg(slack, slack);                        /* 1 - ratio              */
        if (arb_gt(slack, out)) arb_set(out, slack);
    }
    /* out starts at 0 and only grows, so it is the positive part max(0, .) by
     * construction; when the C* lower bound holds strictly it stays 0. */
    arb_clear(sqn);
    arb_clear(slack);
    arb_clear(floor_tol);
    arb_clear(npd);
    arb_clear(nk);
    acb_mat_clear(prod);
    acb_mat_clear(Bd);
}

/* approximate_algebras.tex:432-434 ax_eps_unit:
 *   ||XI - X|| <= eps ||X||,  ||IX - X|| <= eps ||X||,  | ||I|| - 1 | <= eps.
 * Returns max of
 *   (a) ||1_n - Pi_A(1_n)||_op   (is the unit IN A? aic_ecstar_proj_residual),
 *   (b) max_k ||B_k * I - B_k||_op,
 *   (c) max_k ||I * B_k - B_k||_op.
 * EXACTLY 0 when each B_k is a Phi fixed point (X*I = Phi(X.1) = Phi(X) = X for
 * X in Img Phi; .tex:2211 "X*I = X = I*X"): an always-zero invariant for unital
 * Phi with A = Img Phi. ( | ||I||-1 | = 0 since ||1_n||_op = 1, not swept.) */
void aic_ecstar_defect_unit(arb_t out, const aic_ecstar *A, slong prec)
{
    slong n = A->n, d = A->dim_A;
    acb_mat_t I, sBI, sIB, diff;
    acb_mat_init(I, n, n);
    acb_mat_init(sBI, n, n);
    acb_mat_init(sIB, n, n);
    acb_mat_init(diff, n, n);
    acb_mat_one(I);

    /* (a) is the unit IN A? */
    aic_ecstar_proj_residual(out, A, I, prec);

    /* (b),(c) */
    for (slong k = 0; k < d; k++) {
        aic_ecstar_star(sBI, A, A->B[k], I, prec);  /* B_k * I */
        acb_mat_sub(diff, sBI, A->B[k], prec);
        aic_ecstar_bump_opnorm(out, diff, prec);
        aic_ecstar_star(sIB, A, I, A->B[k], prec);  /* I * B_k */
        acb_mat_sub(diff, sIB, A->B[k], prec);
        aic_ecstar_bump_opnorm(out, diff, prec);
    }
    acb_mat_clear(diff);
    acb_mat_clear(sIB);
    acb_mat_clear(sBI);
    acb_mat_clear(I);
}
