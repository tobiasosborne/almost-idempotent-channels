/* aic_funcalc_proj.c — theta(X) and the prop_P projection-extractor (beads
 * aic-4bg, "T-funcalc").
 *
 * approximate_algebras.tex:527-528 — theta and its use:
 *     Ptilde = theta(2P - I),   where   theta(X) = (1/2)(I + sgn(X)).
 * approximate_algebras.tex:524-533 — Proposition prop_P:
 *     Let P satisfy ||P^2 - P|| <= delta < 1/4, and Ptilde = theta(2P - I).
 *     Then Ptilde^2 = Ptilde, Ptilde commutes with P, and
 *     ||Ptilde - P|| <= ||2P - I|| O(delta).
 *   proof (tex:532): set X = 2P - I; then ||X^2 - I|| = 4||P^2 - P|| <= 4 delta < 1,
 *   so sgn(X) is in its validity domain and theta(X) maps it to an exact
 *   idempotent.
 *
 * Constructive route (CLAUDE.md Law 3). The paper's proof IS the algorithm:
 *   1. X = 2P - I.
 *   2. assert ||P^2 - P||_op < 1/4  (equivalently ||X^2 - I|| = 4||P^2-P|| < 1).
 *   3. S = sgn(X)  (eig-free; aic_sgn -> Newton-Schulz, handles the degenerate
 *      spectrum of a near-projector that the eig path cannot — file docstring of
 *      aic_funcalc_sgn.c, bead aic-w4o.1).
 *   4. Ptilde = (I + S)/2.
 * The precondition guard is delta < 1/4 on P, asserted rigorously (arb_lt on the
 * certified ball); this is STRICTER than aic_sgn's own ||X^2-I||<1 guard and is
 * the prop_P hypothesis, so it is checked here explicitly (fail loud, Rule 4).
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_mat.h"
#include "aic_funcalc.h"
#include "aic_funcalc_internal.h"

void aic_theta(acb_mat_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(acb_mat_ncols(X) == n && "aic_theta: X must be square");
    assert(out != X && "aic_theta: out must not alias X");

    /* S = sgn(X) (asserts ||X^2 - I|| < 1 internally), then (I + S)/2. */
    acb_mat_t S, I;
    acb_mat_init(S, n, n);
    acb_mat_init(I, n, n);
    aic_sgn(S, X, prec);
    acb_mat_one(I);
    acb_mat_add(out, I, S, prec);                /* I + sgn(X) */
    acb_mat_scalar_mul_2exp_si(out, out, -1);    /* * 1/2 */

    acb_mat_clear(I);
    acb_mat_clear(S);
}

void aic_prop_P(acb_mat_t out, const acb_mat_t P, slong prec)
{
    slong n = acb_mat_nrows(P);
    assert(acb_mat_ncols(P) == n && "aic_prop_P: P must be square");
    assert(out != P && "aic_prop_P: out must not alias P");

    /* delta = ||P^2 - P||_op as a certified ball; assert delta < 1/4 (tex:525). */
    acb_mat_t P2, defect, X;
    acb_mat_init(P2, n, n);
    acb_mat_init(defect, n, n);
    acb_mat_init(X, n, n);
    acb_mat_sqr(P2, P, prec);
    acb_mat_sub(defect, P2, P, prec);            /* P^2 - P */

    arb_t delta, quarter;
    arb_init(delta);
    arb_init(quarter);
    aic_mat_opnorm(delta, defect, prec);
    arb_set_si(quarter, 1);
    arb_mul_2exp_si(quarter, quarter, -2);       /* 1/4 */
    if (!arb_lt(delta, quarter)) {
        fprintf(stderr, "aic_prop_P: ||P^2 - P||_op not certainly < 1/4 "
                "(prop_P hypothesis delta < 1/4 violated; tex:525)\n");
        abort();
    }

    /* X = 2P - I, then Ptilde = theta(X). */
    acb_mat_scalar_mul_2exp_si(X, P, 1);         /* 2P */
    acb_mat_one(P2);                             /* reuse P2 as I */
    acb_mat_sub(X, X, P2, prec);                 /* 2P - I */
    aic_theta(out, X, prec);

    arb_clear(quarter);
    arb_clear(delta);
    acb_mat_clear(X);
    acb_mat_clear(defect);
    acb_mat_clear(P2);
}
