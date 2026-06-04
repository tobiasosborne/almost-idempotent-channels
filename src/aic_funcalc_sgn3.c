/* aic_funcalc_sgn3.c — the CUBIC (order-3) Newton-Schulz sign candidate (bead
 * aic-09a, the sgn Pareto audition; CLAUDE.md Law 4). A third, independent,
 * inverse-free candidate alongside the quadratic Newton-Schulz and Denman-Beavers
 * (src/aic_funcalc_sgn.c), kept in its own file to honor the ~200 LOC rule.
 *
 * approximate_algebras.tex:514-518 — sgn(X) = X (X^2)^{-1/2}, sgn(X)^2 = I.
 * approximate_algebras.tex:520-521 — validity ||X^2 - I|| < 1.
 *
 * THE ITERATION. The Newton-Schulz family computes (X^2)^{-1/2} by approximating
 * (1 - t)^{-1/2} at t = I - Y^2 with a truncated binomial series and lifting:
 *   order 1 (QUADRATIC conv): P = 1 + 1/2 t          => Y(3I - Y^2)/2   [the quadratic NS]
 *   order 2 (CUBIC conv):     P = 1 + 1/2 t + 3/8 t^2 => Y(15I - 10 Y^2 + 3 Y^4)/8
 * The order-2 truncation gives a cubically-convergent iteration:
 *   Y_0 = X,   Y_{k+1} = (1/8) Y_k (15 I - 10 Y_k^2 + 3 Y_k^4).
 * Verified expansion (t = I - Y^2): 1 + 1/2(I - Y^2) + 3/8(I - Y^2)^2
 *   = 1 + 1/2 + 3/8 + (-1/2 - 3/4) Y^2 + 3/8 Y^4 = 15/8 - 10/8 Y^2 + 3/8 Y^4. QED.
 *
 * COST / PARETO. 3 matmuls/step (Y^2, Y^4=(Y^2)^2, Y*poly) vs 2 for the quadratic
 * NS, but cubic vs quadratic convergence => fewer steps to a fixed prec. Whether
 * the net wins depends on the regime (cluster vs near-discontinuity spectrum,
 * size n); that is exactly what bench/bench_funcalc.c auditions (Law 4 — keep the
 * Pareto frontier, dispatch per regime). Inverse-free, so no near-singular stall.
 *
 * RIGOR. Identical discipline to the quadratic NS (bead aic-1vp, FINDINGS §C11):
 * assert the basin ||X^2 - I|| < 1 on the radius-carrying X, iterate on mid(X)
 * (radius stripped once) so an inherited wide radius cannot drag the midpoint off
 * the involution, then restore rigor a-posteriori via aic_funcalc_int_certify_sign
 * (||Y^2 - I|| AND ||YX - XY|| < tol against the ORIGINAL X). Out of basin =>
 * fail loud (Rule 4). On a zero-radius input the strip is a no-op.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/mag.h>

#include "aic/aic_funcalc.h"
#include "aic_funcalc_internal.h"

#define AIC_SGN3_MAX_ITERS 100

void aic_sgn_newton_schulz3(acb_mat_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(out != X && "aic_sgn_newton_schulz3: out must not alias X");

    arb_t bound;
    arb_init(bound);
    aic_funcalc_int_assert_domain(bound, X, prec); /* ||X^2-I|| < 1 (tex:521) */
    arb_clear(bound);

    arb_t tol, step;
    arb_init(tol);
    arb_init(step);
    aic_funcalc_int_tol(tol, prec);

    mag_t in_rad;
    mag_init(in_rad);

    acb_mat_t Y, Ynext, Y2, Y4, poly, scratch;
    acb_mat_init(Y, n, n);
    acb_mat_init(Ynext, n, n);
    acb_mat_init(Y2, n, n);
    acb_mat_init(Y4, n, n);
    acb_mat_init(poly, n, n);
    acb_mat_init(scratch, n, n);
    aic_funcalc_int_to_midpoint(Y, X, in_rad); /* Y_0 = mid(X) */

    int k;
    for (k = 0; k < AIC_SGN3_MAX_ITERS; k++) {
        acb_mat_sqr(Y2, Y, prec);                 /* Y^2 */
        acb_mat_sqr(Y4, Y2, prec);                /* Y^4 = (Y^2)^2 */

        /* poly = 15 I - 10 Y^2 + 3 Y^4 */
        acb_mat_scalar_mul_si(poly, Y4, 3, prec);     /* 3 Y^4 */
        acb_mat_scalar_mul_si(scratch, Y2, 10, prec); /* 10 Y^2 */
        acb_mat_sub(poly, poly, scratch, prec);       /* 3 Y^4 - 10 Y^2 */
        acb_mat_one(scratch);
        acb_mat_scalar_mul_si(scratch, scratch, 15, prec); /* 15 I */
        acb_mat_add(poly, poly, scratch, prec);            /* + 15 I */

        acb_mat_mul(Ynext, Y, poly, prec);            /* Y (15I - 10Y^2 + 3Y^4) */
        acb_mat_scalar_mul_2exp_si(Ynext, Ynext, -3); /* * 1/8 */

        aic_funcalc_int_step_norm(step, Ynext, Y, prec);
        acb_mat_set(Y, Ynext);
        if (arb_lt(step, tol)) { k++; break; }
    }
    if (k >= AIC_SGN3_MAX_ITERS) {
        fprintf(stderr, "aic_sgn_newton_schulz3: no convergence in %d iters "
                "(bug: ||X^2-I||<1 should converge cubically)\n",
                AIC_SGN3_MAX_ITERS);
        abort();
    }
    aic_funcalc_int_certify_sign(Y, X, in_rad, "aic_sgn_newton_schulz3", prec);
    acb_mat_set(out, Y);

    mag_clear(in_rad);
    acb_mat_clear(scratch);
    acb_mat_clear(poly);
    acb_mat_clear(Y4);
    acb_mat_clear(Y2);
    acb_mat_clear(Ynext);
    acb_mat_clear(Y);
    arb_clear(step);
    arb_clear(tol);
}
