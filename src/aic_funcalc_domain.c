/* aic_funcalc_domain.c — the shared validity-domain guard for the funcalc module
 * (beads aic-4bg, "T-funcalc"). See src/aic_funcalc_internal.h.
 *
 * approximate_algebras.tex:516 — the functions are defined only on a ball:
 *     |X| = (X^2)^{1/2},  sgn(X) = X (X^2)^{-1/2}     ( ||X^2 - x0^2 I|| < x0^2 ).
 * approximate_algebras.tex:520-521 — the bound that makes sgn the nearest exact
 * solution of Y^2 = I holds in the same ball with x0 = 1:
 *     ||sgn(X) - X|| <= ||X|| O(||X^2 - I||)   if  ||X^2 - I|| < 1.
 * This module fixes x0 = 1, so the precondition is the certified statement
 * ||X^2 - I||_op < 1. We use the OPERATOR norm (aic_mat_opnorm), not Frobenius:
 * the convergence of Newton-Schulz / the binomial series is governed by the
 * spectral radius of I - X^2, and the operator norm is its tight upper bound
 * (Frobenius would over-count off the boundary and could reject a valid input).
 *
 * Fail-loud (CLAUDE.md Rule 4). The check is RIGOROUS: arb_lt on the certified
 * ball ||X^2 - I|| is nonzero only when every value in the ball is < 1, so a
 * ball that straddles 1 (precision lost, or genuinely on the boundary) aborts
 * rather than proceeding into a divergent iteration.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_mat.h"
#include "aic_funcalc_internal.h"

void aic_funcalc_int_tol(arb_t tol, slong prec)
{
    arb_set_si(tol, 1);
    arb_mul_2exp_si(tol, tol, -(prec - 8));
}

void aic_funcalc_int_def_X2(arb_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(acb_mat_ncols(X) == n && "funcalc input must be square");

    acb_mat_t X2, D;
    acb_mat_init(X2, n, n);
    acb_mat_init(D, n, n);
    acb_mat_sqr(X2, X, prec);
    acb_mat_one(D);                 /* D = I */
    acb_mat_sub(D, X2, D, prec);    /* D = X^2 - I */
    aic_mat_opnorm(out, D, prec);   /* certified ||X^2 - I||_op */

    acb_mat_clear(D);
    acb_mat_clear(X2);
}

void aic_funcalc_int_assert_domain(arb_t bound, const acb_mat_t X, slong prec)
{
    arb_t one;
    arb_init(one);
    arb_one(one);

    aic_funcalc_int_def_X2(bound, X, prec);
    if (!arb_lt(bound, one)) {
        fprintf(stderr,
                "aic_funcalc: input out of validity domain: "
                "||X^2 - I||_op is not certainly < 1 "
                "(tex:516,520; raise prec or shrink ||X^2-I||)\n");
        abort();
    }

    arb_clear(one);
}
