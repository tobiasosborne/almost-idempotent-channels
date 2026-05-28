/* aic_funcalc_pow.c — x^alpha of a matrix via the binomial Taylor series, and
 * |X| (beads aic-4bg, "T-funcalc"). This is the power-series functional-calculus
 * route (the second sgn audition input would use it for (X^2)^{-1/2}; here it is
 * the closed route for |X| = (X^2)^{1/2}).
 *
 * approximate_algebras.tex:503-511 — Taylor/power-series functional calculus:
 *     f(X) = a_0 I + sum_{n>=1} a_n (X - x0 I)^n              if ||X - x0 I|| < r,
 *     ||f(X) - f(x0 I)|| <= sum_{n>=1} |a_n| ||X - x0 I||^n.
 *   "let f(x) = x^alpha ... Then f(X) is well-defined if ||X - x0 I|| < x0."
 * For f(x) = x^alpha expanded at x0 the coefficients are the binomial series
 *     x^alpha = x0^alpha sum_{n>=0} C(alpha,n) (x/x0 - 1)^n,
 *     C(alpha,n) = prod_{j=0}^{n-1} (alpha - j)/(j+1),
 * convergent for |x/x0 - 1| < 1, i.e. ||A - x0 I|| < x0 (the tex:511 condition).
 *
 * Rigorous truncation (CLAUDE.md Rule 4: certified, not heuristic). Let
 * q = ||A/x0 - I||_op (operator norm; tight upper bound for the series radius).
 * Term n has norm <= |C(alpha,n)| q^n. Once n > |alpha| the ratio of successive
 * term-bounds is rho_n = |alpha - n|/(n+1) * q; for q < 1 and large n,
 * rho_n -> q < 1, so the tail from term N is bounded by the geometric sum
 * term_N * rho_N/(1 - rho_N) when rho_N < 1. We sum until that tail bound is
 * below a prec-appropriate tol, then INFLATE the result's error balls by the
 * tail bound (acb_mat_add_error_mag) so the returned enclosure rigorously
 * contains A^alpha. A hard cap (200 terms) aborts loudly: at q < 1 the tail
 * crosses tol well before then; hitting it means q is at the boundary (the
 * caller violated ||A - x0 I|| < x0).
 *
 * |X| = (X^2)^{1/2} (tex:514): aic_funcalc_xpow(X^2, alpha=1/2, x0=1). The
 * domain ||X^2 - I|| < 1 is exactly this module's normalisation, asserted via
 * the shared guard.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/mag.h>

#include "aic_mat.h"
#include "aic_funcalc.h"
#include "aic_funcalc_internal.h"

#define AIC_XPOW_MAX_TERMS 200

void aic_funcalc_xpow(acb_mat_t out, const acb_mat_t A, double alpha, double x0,
                      slong prec)
{
    slong n = acb_mat_nrows(A);
    assert(acb_mat_ncols(A) == n && "aic_funcalc_xpow: A must be square");
    assert(out != A && "aic_funcalc_xpow: out must not alias A");
    assert(x0 > 0.0 && "aic_funcalc_xpow: x0 must be positive (tex:511)");

    /* B = A/x0 - I, and q = ||B||_op (the series variable and its radius). */
    acb_mat_t B, term, Bpow, scratch;
    acb_mat_init(B, n, n);
    acb_mat_init(term, n, n);   /* C(alpha,k) * B^k */
    acb_mat_init(Bpow, n, n);   /* B^k */
    acb_mat_init(scratch, n, n);

    arb_t inv_x0;
    arb_init(inv_x0);
    arb_set_d(inv_x0, x0);
    arb_inv(inv_x0, inv_x0, prec);
    acb_mat_scalar_mul_arb(B, A, inv_x0, prec); /* A/x0 */
    acb_mat_one(scratch);
    acb_mat_sub(B, B, scratch, prec);           /* A/x0 - I */

    arb_t q;
    arb_init(q);
    aic_mat_opnorm(q, B, prec);
    arb_t one;
    arb_init(one);
    arb_one(one);
    if (!arb_lt(q, one)) {
        fprintf(stderr, "aic_funcalc_xpow: ||A/x0 - I||_op not certainly < 1 "
                "(tex:511 domain violated; x0=%g)\n", x0);
        abort();
    }

    arb_t tol, coeff, qpow, termbound, tailbound, ratio, num;
    arb_init(tol); arb_init(coeff); arb_init(qpow);
    arb_init(termbound); arb_init(tailbound); arb_init(ratio); arb_init(num);
    aic_funcalc_int_tol(tol, prec);

    /* out := C(alpha,0) I = I; Bpow := B^0 = I; coeff := 1; qpow := q^0 = 1. */
    acb_mat_one(out);           /* a_0 I, a_0 = 1 */
    acb_mat_one(Bpow);
    arb_set_si(coeff, 1);
    arb_set_si(qpow, 1);

    int k;
    int converged = 0;
    for (k = 1; k < AIC_XPOW_MAX_TERMS; k++) {
        /* coeff_k = coeff_{k-1} * (alpha - (k-1)) / k  ;  qpow_k = q^k. */
        arb_set_d(num, alpha - (double) (k - 1));
        arb_mul(coeff, coeff, num, prec);
        arb_div_si(coeff, coeff, k, prec);
        arb_mul(qpow, qpow, q, prec);

        /* term_k = C(alpha,k) * B^k. */
        acb_mat_mul(scratch, Bpow, B, prec); /* B^k */
        acb_mat_set(Bpow, scratch);
        acb_mat_scalar_mul_arb(term, Bpow, coeff, prec);
        acb_mat_add(out, out, term, prec);

        /* termbound_k = |coeff_k| q^k. Geometric tail from term k+1:
         *   rho = |alpha - k|/(k+1) * q ;  tail = termbound_k * rho/(1-rho)
         * valid once rho < 1 (k > |alpha| guarantees it for q < 1). */
        arb_abs(termbound, coeff);
        arb_mul(termbound, termbound, qpow, prec);

        arb_set_d(num, alpha - (double) k);
        arb_abs(ratio, num);
        arb_div_si(ratio, ratio, k + 1, prec);
        arb_mul(ratio, ratio, q, prec); /* rho */
        if (arb_lt(ratio, one)) {
            arb_sub(num, one, ratio, prec);       /* 1 - rho */
            arb_div(tailbound, ratio, num, prec);  /* rho/(1-rho) */
            arb_mul(tailbound, tailbound, termbound, prec);
            if (arb_lt(tailbound, tol)) { converged = 1; break; }
        }
    }
    if (!converged) {
        fprintf(stderr, "aic_funcalc_xpow: series tail not below tol in %d "
                "terms (q at boundary?)\n", AIC_XPOW_MAX_TERMS);
        abort();
    }

    /* Inflate every entry by the rigorous tail bound so the enclosure is sound,
     * then scale by x0^alpha. */
    mag_t tailmag;
    mag_init(tailmag);
    arb_get_mag(tailmag, tailbound);
    acb_mat_add_error_mag(out, tailmag);
    mag_clear(tailmag);

    arb_set_d(num, x0);
    arb_set_d(coeff, alpha);
    arb_pow(num, num, coeff, prec); /* x0^alpha */
    acb_mat_scalar_mul_arb(out, out, num, prec);

    arb_clear(num); arb_clear(ratio); arb_clear(tailbound);
    arb_clear(termbound); arb_clear(qpow); arb_clear(coeff); arb_clear(tol);
    arb_clear(one); arb_clear(q); arb_clear(inv_x0);
    acb_mat_clear(scratch);
    acb_mat_clear(Bpow);
    acb_mat_clear(term);
    acb_mat_clear(B);
}

void aic_abs(acb_mat_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(acb_mat_ncols(X) == n && "aic_abs: X must be square");
    assert(out != X && "aic_abs: out must not alias X");

    /* |X| = (X^2)^{1/2}, valid on ||X^2 - I|| < 1 (tex:514,516). The xpow guard
     * re-checks ||X^2/1 - I|| < 1, the same statement. */
    acb_mat_t X2;
    acb_mat_init(X2, n, n);
    acb_mat_sqr(X2, X, prec);
    aic_funcalc_xpow(out, X2, 0.5, 1.0, prec);
    acb_mat_clear(X2);
}
