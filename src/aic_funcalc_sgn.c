/* aic_funcalc_sgn.c — the sgn(X) audition (beads aic-4bg, "T-funcalc"; CLAUDE.md
 * Law 4: audition, don't presume). Two independent, eig-free, iterative
 * candidates for the matrix sign function, plus the default dispatcher.
 *
 * approximate_algebras.tex:514-518 — definition and the defining identity:
 *     sgn(X) = X (X^2)^{-1/2},   and   sgn(X)^2 = I   for all X (in the domain).
 * approximate_algebras.tex:520-521 — the validity bound:
 *     ||sgn(X) - X|| <= ||X|| O(||X^2 - I||)   if  ||X^2 - I|| < 1.
 *
 * Both candidates are EIG-FREE. This is the whole point for this project: the
 * real inputs are X = 2P - I for a projection P (degenerate spectrum {+/-1}),
 * and FLINT's certified eig path (aic_mat_eig_hermitian) requires a SIMPLE
 * spectrum and ABORTS on the repeated eigenvalues of a projection (bead
 * aic-w4o.1). An eig-based sgn is therefore not auditionable on the inputs that
 * matter; the iterative routes below sidestep the spectrum entirely.
 *
 * CANDIDATE 1 — Newton-Schulz (inverse-free). The scalar step y_{k+1} =
 * 1/2 y_k (3 - y_k^2) is Newton for 1/sqrt on y^2 (fixed points {-1,0,+1});
 * from y_0 = x it converges to sgn(x) iff |1 - x^2| < 1. Lifted to matrices,
 *   Y_0 = X,  Y_{k+1} = 1/2 Y_k (3 I - Y_k^2),
 * converges quadratically to sgn(X) when ||I - X^2|| < 1, using only matrix
 * multiply/add (no inverse) — so it never hits a near-singular solve.
 *
 * CANDIDATE 2 — Denman-Beavers / Newton sign iteration (one inverse per step;
 * shard-B L50). The scalar Newton step for sign on x is y_{k+1} = 1/2(y_k + 1/y_k),
 * lifted to matrices:
 *   Y_0 = X,  Y_{k+1} = 1/2 (Y_k + Y_k^{-1}),
 * which drives Y_k -> sgn(X) quadratically when X has no purely-imaginary
 * eigenvalue (here implied by ||I - X^2|| < 1). This is the textbook matrix-sign
 * Newton iteration (the single-variable Denman-Beavers form; the COUPLED (Y,Z)
 * form with Z_0 = I was auditioned and REJECTED here — it has a z_0 = I
 * degeneracy at eigenvalue -1 that makes Y_1 singular for our X = 2P - I inputs,
 * whose -1 eigenvalues are exactly the rejection trigger). We keep this inverse-
 * based candidate as the independent cross-check that catches an algorithmic
 * error in either route (tests/test_funcalc.c asserts the two agree).
 *
 * DEFAULT (aic_sgn) — set to Newton-Schulz by bench evidence. bench_funcalc.c
 * times both at prec=53 for n=4,8,16 on a perturbed projector: Newton-Schulz is
 * faster (matrix products only) than Denman-Beavers (a certified acb_mat_inv per
 * step is the dominant cost) AND is inverse-free, so it cannot stall on a
 * near-singular intermediate. Both pass the degenerate-projector cross-check
 * (sgn(2P-I) = 2P-I exactly within tol); correctness is a tie, speed + no-inverse
 * breaks it. See the BENCH lines reproduced in the bead aic-4bg report.
 *
 * Convergence loop. Iterate until the midpoint step ||Y_{k+1} - Y_k|| (see
 * step_norm) drops below a prec-appropriate tol. A hard cap of 100 steps is a
 * fail-loud abort: at ||I-X^2|| < 1 the quadratic rate reaches prec=256 in < 12
 * steps, so hitting the cap means a bug, not slow input.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_mat.h"
#include "aic_funcalc.h"
#include "aic_funcalc_internal.h"

#define AIC_SGN_MAX_ITERS 100

/* The midpoint convergence measure step_norm lives in aic_funcalc_domain.c as
 * aic_funcalc_int_step_norm (shared with the global Newton variant, bead aic-8hz);
 * see its header docstring for why midpoints, not certified balls. */

void aic_sgn_newton_schulz(acb_mat_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(out != X && "aic_sgn_newton_schulz: out must not alias X");
    arb_t bound;
    arb_init(bound);
    aic_funcalc_int_assert_domain(bound, X, prec); /* ||X^2-I|| < 1 (tex:521) */
    arb_clear(bound);

    arb_t tol, step;
    arb_init(tol);
    arb_init(step);
    aic_funcalc_int_tol(tol, prec);

    acb_mat_t Y, Ynext, T, three_I;
    acb_mat_init(Y, n, n);
    acb_mat_init(Ynext, n, n);
    acb_mat_init(T, n, n);
    acb_mat_init(three_I, n, n);
    acb_mat_set(Y, X);                       /* Y_0 = X */

    int k;
    for (k = 0; k < AIC_SGN_MAX_ITERS; k++) {
        acb_mat_sqr(T, Y, prec);             /* T = Y^2 */
        acb_mat_one(three_I);
        acb_mat_scalar_mul_si(three_I, three_I, 3, prec); /* three_I = 3 I */
        acb_mat_sub(T, three_I, T, prec);    /* T = 3 I - Y^2 */
        acb_mat_mul(Ynext, Y, T, prec);      /* Y (3 I - Y^2) */
        acb_mat_scalar_mul_2exp_si(Ynext, Ynext, -1); /* * 1/2 */

        aic_funcalc_int_step_norm(step, Ynext, Y, prec);
        acb_mat_set(Y, Ynext);
        if (arb_lt(step, tol)) { k++; break; } /* converged below tol */
    }
    if (k >= AIC_SGN_MAX_ITERS) {
        fprintf(stderr, "aic_sgn_newton_schulz: no convergence in %d iters "
                "(bug: ||X^2-I||<1 should converge quadratically)\n",
                AIC_SGN_MAX_ITERS);
        abort();
    }
    acb_mat_set(out, Y);

    acb_mat_clear(three_I);
    acb_mat_clear(T);
    acb_mat_clear(Ynext);
    acb_mat_clear(Y);
    arb_clear(step);
    arb_clear(tol);
}

void aic_sgn_denman_beavers(acb_mat_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(out != X && "aic_sgn_denman_beavers: out must not alias X");
    arb_t bound;
    arb_init(bound);
    aic_funcalc_int_assert_domain(bound, X, prec);
    arb_clear(bound);

    arb_t tol, step;
    arb_init(tol);
    arb_init(step);
    aic_funcalc_int_tol(tol, prec);

    acb_mat_t Y, Yinv, Ynext;
    acb_mat_init(Y, n, n);
    acb_mat_init(Yinv, n, n);
    acb_mat_init(Ynext, n, n);
    acb_mat_set(Y, X);                       /* Y_0 = X */

    int k;
    for (k = 0; k < AIC_SGN_MAX_ITERS; k++) {
        /* Certified inverse; acb_mat_inv returns 0 if it cannot PROVE
         * invertibility at this prec (Rule 4: that is a loud failure). For an
         * exact sign matrix (X^2 = I) the first step already gives Y_1 = X
         * (Y^{-1} = Y), so the degenerate-projector case converges in one step. */
        if (!acb_mat_inv(Yinv, Y, prec)) {
            fprintf(stderr, "aic_sgn_denman_beavers: could not certify the "
                    "inverse at iter %d, prec=%ld (raise prec)\n",
                    k, (long) prec);
            abort();
        }
        acb_mat_add(Ynext, Y, Yinv, prec);   /* Y + Y^{-1} */
        acb_mat_scalar_mul_2exp_si(Ynext, Ynext, -1); /* * 1/2 */

        aic_funcalc_int_step_norm(step, Ynext, Y, prec);
        acb_mat_set(Y, Ynext);
        if (arb_lt(step, tol)) { k++; break; }
    }
    if (k >= AIC_SGN_MAX_ITERS) {
        fprintf(stderr, "aic_sgn_denman_beavers: no convergence in %d iters\n",
                AIC_SGN_MAX_ITERS);
        abort();
    }
    acb_mat_set(out, Y);

    acb_mat_clear(Ynext);
    acb_mat_clear(Yinv);
    acb_mat_clear(Y);
    arb_clear(step);
    arb_clear(tol);
}

/* Default dispatch (bead aic-8hz): AUTO-DISPATCH on the operator-norm basin.
 * When ||X^2-I||_op is certified < 1 (the Newton-Schulz local-convergence basin),
 * use the fast inverse-free Newton-Schulz — BYTE-FOR-BYTE the prior behavior, so
 * the existing in-basin tests are unchanged (and the bench-selected winner on the
 * degenerate-projector inputs this project actually uses; see the file docstring).
 * Out of basin, fall back to the globally-convergent aic_sgn_newton_global, which
 * certifies the SPECTRAL precondition rho(I-X^2)<1 (Higham Thm 5.6) instead of the
 * op-norm bound. This is a STRICT improvement: the out-of-basin input previously
 * hit Newton-Schulz's domain assert and ABORTED. The probe is non-aborting. */
void aic_sgn(acb_mat_t out, const acb_mat_t X, slong prec)
{
    if (aic_funcalc_int_in_op_basin(X, prec))
        aic_sgn_newton_schulz(out, X, prec);
    else
        aic_sgn_newton_global(out, X, prec);
}
