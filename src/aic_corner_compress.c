/* aic_corner_compress.c — the compression maps L_P, R_Q, Co_{P,Q} and the
 * coordinate<->operator apply for section 7 (approximate_algebras.tex:1052-1066).
 * See include/aic_corner.h for the full why (the star-not-plain-product callout,
 * the two coordinate systems, the prop_P basin). Extraction of S_{P,Q} lives in
 * src/aic_corner_extract.c (double-path SVD).
 *
 * approximate_algebras.tex:1054-1057 — the compression map (verbatim):
 *   "there are two slightly different maps that could be called 'compressions':
 *    L_P R_Q : X |-> P(XQ) and R_Q L_P : X |-> (PX)Q. We will use their symmetric
 *    combination, 1/2(L_P R_Q + R_Q L_P). It is an O(delta+eps)-idempotent ...
 *    can be approximated by an idempotent using Proposition prop_P. The
 *    compression map thus defined,
 *        Co_{P,Q} : A -> A,   Co_{P,Q} = theta(L_P R_Q + R_Q L_P - 1)."
 * approximate_algebras.tex:1064 — the one-sided variants:
 *   "There is a variant of this construction where only L_P or only R_Q is
 *    applied, and one writes Co_{P,1}, S_{P,1} or Co_{1,Q}, S_{1,Q}."
 *
 * Constructive route (CLAUDE.md Law 3). In finite dim an operator on the d-dim
 * space A is a d x d matrix in the orthonormal {B_k} basis. L_P R_Q + R_Q L_P - 1
 * = 2M - I with M = 1/2(L_P R_Q + R_Q L_P), so Co = aic_prop_P(M): prop_P IS
 * theta(2M-1) and asserts the basin ||M^2-M|| < 1/4 (= ||(2M-1)^2-I||/4 < 1/4,
 * the funcalc sgn domain). The d x d products L_P R_Q etc. are ORDINARY matrix
 * products of the coordinate matrices (NOT products in A); only the entries
 * (L_P)_{ij} = <B_i, P*B_j>_F route through the star.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_corner.h"
#include "aic_corner_internal.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_funcalc.h"
#include "aic/aic_mat.h"

/* <B_i, Y>_F = Tr(B_i^dag Y) = sum_{a,b} conj(B_i[a,b]) Y[a,b] into `out`. Shared
 * with src/aic_corner_extract.c (aic_corner_apply) via aic_corner_internal.h. */
void aic_corner_frob_ip(acb_t out, const acb_mat_t Bi, const acb_mat_t Y,
                        slong prec)
{
    slong n = acb_mat_nrows(Y);
    acb_t t;
    acb_init(t);
    acb_zero(out);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++) {
            acb_conj(t, acb_mat_entry(Bi, a, b));
            acb_mul(t, t, acb_mat_entry(Y, a, b), prec);
            acb_add(out, out, t, prec);
        }
    acb_clear(t);
}

/* ASSERT P^dag = P up to a prec-appropriate tol (a delta-projection is Hermitian,
 * .tex:917; fail loud, Rule 4). */
static void assert_hermitian(const acb_mat_t P, slong prec)
{
    slong n = acb_mat_nrows(P);
    assert(acb_mat_ncols(P) == n && "corner: P/Q must be square");
    acb_mat_t Pd, D;
    arb_t e, tol;
    acb_mat_init(Pd, n, n);
    acb_mat_init(D, n, n);
    arb_init(e);
    arb_init(tol);
    acb_mat_conjugate_transpose(Pd, P);
    acb_mat_sub(D, P, Pd, prec);
    aic_mat_frobenius_norm(e, D, prec);
    arb_set_d(tol, 1e-9);
    assert(arb_le(e, tol) &&
           "corner: P/Q not Hermitian (a delta-projection has P^dag=P, .tex:917)");
    arb_clear(tol);
    arb_clear(e);
    acb_mat_clear(D);
    acb_mat_clear(Pd);
}

/* (L_P)_{ij} = <B_i, P * B_j>_F, P * B_j = aic_ecstar_star(A, P, B_j). */
void aic_corner_build_L(acb_mat_t Lp, const aic_ecstar *A, const acb_mat_t P,
                        slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(acb_mat_nrows(Lp) == d && acb_mat_ncols(Lp) == d &&
           "build_L: Lp must be dim_A x dim_A");
    assert(acb_mat_nrows(P) == n && "build_L: P must be n x n");
    assert_hermitian(P, prec);

    acb_mat_t PBj;
    acb_mat_init(PBj, n, n);
    for (slong j = 0; j < d; j++) {
        aic_ecstar_star(PBj, A, P, A->B[j], prec);          /* P * B_j */
        for (slong i = 0; i < d; i++)
            aic_corner_frob_ip(acb_mat_entry(Lp, i, j), A->B[i], PBj, prec);
    }
    acb_mat_clear(PBj);
}

/* (R_Q)_{ij} = <B_i, B_j * Q>_F, B_j * Q = aic_ecstar_star(A, B_j, Q). */
void aic_corner_build_R(acb_mat_t Rq, const aic_ecstar *A, const acb_mat_t Q,
                        slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(acb_mat_nrows(Rq) == d && acb_mat_ncols(Rq) == d &&
           "build_R: Rq must be dim_A x dim_A");
    assert(acb_mat_nrows(Q) == n && "build_R: Q must be n x n");
    assert_hermitian(Q, prec);

    acb_mat_t BjQ;
    acb_mat_init(BjQ, n, n);
    for (slong j = 0; j < d; j++) {
        aic_ecstar_star(BjQ, A, A->B[j], Q, prec);          /* B_j * Q */
        for (slong i = 0; i < d; i++)
            aic_corner_frob_ip(acb_mat_entry(Rq, i, j), A->B[i], BjQ, prec);
    }
    acb_mat_clear(BjQ);
}

/* Co = theta(2M-1) for a d x d near-idempotent coordinate matrix M, via
 * aic_prop_P (asserts the prop_P basin ||M^2-M|| < 1/4). Shared by Co/Co_P1/Co_1Q. */
static void co_from_M(acb_mat_t Co, const acb_mat_t M, slong prec)
{
    aic_prop_P(Co, M, prec);    /* prop_P = theta(2M-1); fail-loud basin inside */
}

void aic_corner_Co(acb_mat_t Co, const aic_ecstar *A, const acb_mat_t P,
                   const acb_mat_t Q, slong prec)
{
    slong d = A->dim_A;
    assert(acb_mat_nrows(Co) == d && acb_mat_ncols(Co) == d &&
           "Co: must be dim_A x dim_A");
    acb_mat_t Lp, Rq, LR, RL, M;
    acb_mat_init(Lp, d, d);
    acb_mat_init(Rq, d, d);
    acb_mat_init(LR, d, d);
    acb_mat_init(RL, d, d);
    acb_mat_init(M, d, d);
    aic_corner_build_L(Lp, A, P, prec);
    aic_corner_build_R(Rq, A, Q, prec);
    acb_mat_mul(LR, Lp, Rq, prec);          /* L_P R_Q */
    acb_mat_mul(RL, Rq, Lp, prec);          /* R_Q L_P */
    acb_mat_add(M, LR, RL, prec);           /* L_P R_Q + R_Q L_P = 2M */
    acb_mat_scalar_mul_2exp_si(M, M, -1);   /* M = 1/2(L_P R_Q + R_Q L_P) */
    co_from_M(Co, M, prec);                 /* theta(2M-1) */
    acb_mat_clear(M);
    acb_mat_clear(RL);
    acb_mat_clear(LR);
    acb_mat_clear(Rq);
    acb_mat_clear(Lp);
}

void aic_corner_Co_P1(acb_mat_t Co, const aic_ecstar *A, const acb_mat_t P,
                      slong prec)
{
    /* R_1 = id, so 1/2(L_P R_1 + R_1 L_P) = L_P; Co_{P,1} = theta(2 L_P - 1)
     * = aic_prop_P(L_P). (.tex:1064) */
    slong d = A->dim_A;
    assert(acb_mat_nrows(Co) == d && acb_mat_ncols(Co) == d &&
           "Co_P1: must be dim_A x dim_A");
    acb_mat_t Lp;
    acb_mat_init(Lp, d, d);
    aic_corner_build_L(Lp, A, P, prec);
    co_from_M(Co, Lp, prec);
    acb_mat_clear(Lp);
}

void aic_corner_Co_1Q(acb_mat_t Co, const aic_ecstar *A, const acb_mat_t Q,
                      slong prec)
{
    /* L_1 = id => 1/2(L_1 R_Q + R_Q L_1) = R_Q; Co_{1,Q} = theta(2 R_Q - 1). */
    slong d = A->dim_A;
    assert(acb_mat_nrows(Co) == d && acb_mat_ncols(Co) == d &&
           "Co_1Q: must be dim_A x dim_A");
    acb_mat_t Rq;
    acb_mat_init(Rq, d, d);
    aic_corner_build_R(Rq, A, Q, prec);
    co_from_M(Co, Rq, prec);
    acb_mat_clear(Rq);
}
