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
#include <flint/arf.h>
#include <flint/mag.h>

#include "aic/aic_mat.h"
#include "aic_funcalc_internal.h"

/* See the header docstring (aic_funcalc_internal.h). Computed on the ball
 * MIDPOINTS so the convergence test tracks "the iterate has stopped moving",
 * not the (monotonically inflating) certified radius. */
void aic_funcalc_int_step_norm(arb_t out, const acb_mat_t A, const acb_mat_t B,
                               slong prec)
{
    slong n = acb_mat_nrows(A);
    acb_mat_t D;
    acb_mat_init(D, n, n);
    acb_mat_sub(D, A, B, prec);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
            acb_get_mid(acb_mat_entry(D, i, j), acb_mat_entry(D, i, j));
    aic_mat_frobenius_norm(out, D, prec);
    acb_mat_clear(D);
}

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

/* See the header docstring (aic_funcalc_internal.h). bead aic-1vp / FINDINGS §C11. */
void aic_funcalc_int_to_midpoint(acb_mat_t out, const acb_mat_t X, mag_t in_rad)
{
    slong n = acb_mat_nrows(X);
    assert(out != X && "aic_funcalc_int_to_midpoint: out must not alias X");
    assert(acb_mat_ncols(X) == n && "funcalc midpoint: X must be square");
    mag_zero(in_rad);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            const acb_struct *e = acb_mat_entry(X, i, j);
            if (mag_cmp(arb_radref(acb_realref(e)), in_rad) > 0)
                mag_set(in_rad, arb_radref(acb_realref(e)));
            if (mag_cmp(arb_radref(acb_imagref(e)), in_rad) > 0)
                mag_set(in_rad, arb_radref(acb_imagref(e)));
            acb_get_mid(acb_mat_entry(out, i, j), e);
        }
}

/* See the header docstring (aic_funcalc_internal.h). bead aic-1vp / FINDINGS §C11. */
void aic_funcalc_int_certify_sign(const acb_mat_t Y, const acb_mat_t X,
                                  const mag_t in_rad, const char *who, slong prec)
{
    slong n = acb_mat_nrows(Y);
    /* cert tol = max(2^-(prec/2), in_rad) * 8 * sqrt(n) — a SANITY BACKSTOP, not a
     * tight enclosure of sgn(X_true). In the operating regime the prec-floor
     * 2^-(prec/2) DOMINATES (e.g. ~2.9e-39 at prec=256, vs an inherited input radius
     * ~1e-70), so the certificate rejects GROSS failures (a non-converged / wrong Y
     * has ||Y^2-I|| or ||[X,Y]|| = O(1) >> tol) while admitting a converged midpoint
     * sign (||Y^2-I|| ~ the midpoint iteration's machine floor ~1e-74 << tol). The
     * in_rad term only binds in the unusual case in_rad > 2^-(prec/2). 8*sqrt(n)
     * absorbs the O(1) Frobenius constants over n^2 entries (NOTE the binding ||XY-YX||
     * arm IS computed on the radius-carrying X, so it tests sgn(mid X) against the
     * TRUE input ball). The ACTUAL soundness is NOT this tol magnitude — it is (i) the
     * away-from-0 basin/Gelfand precondition (asserted on the radius-carrying X BEFORE
     * to_midpoint: rho(I-X^2)<1 keeps spec(X) away from the imaginary axis, so
     * sgn(mid X)=sgn(X_true)), and (ii) the Y0=mid(X) seeding (Higham Functions of
     * Matrices Thm 5.6: Newton/NS converges to the correct-INERTIA sign). The
     * certificate is the loud-failure backstop on top of those (bead aic-1vp,
     * FINDINGS §C11; the certificate alone cannot distinguish +sgn from -sgn). */
    arb_t cert_tol, cn, sqn;
    arb_init(cert_tol);
    arb_init(cn);
    arb_init(sqn);
    arb_set_si(sqn, n);
    arb_sqrt(sqn, sqn, prec);
    {
        arb_t floor_r, rad_r;
        arf_t rad_a;
        arb_init(floor_r);
        arb_init(rad_r);
        arf_init(rad_a);
        arb_set_si(floor_r, 1);
        arb_mul_2exp_si(floor_r, floor_r, -(prec / 2)); /* 2^-(prec/2) */
        arf_set_mag(rad_a, in_rad);                      /* in_rad as an arf */
        arb_set_arf(rad_r, rad_a);                       /* ... then an arb */
        arb_max(cert_tol, floor_r, rad_r, prec);
        arf_clear(rad_a);
        arb_clear(rad_r);
        arb_clear(floor_r);
    }
    arb_mul_si(cert_tol, cert_tol, 8, prec);
    arb_mul(cert_tol, cert_tol, sqn, prec);

    acb_mat_t Y2, D, XY, YX;
    acb_mat_init(Y2, n, n);
    acb_mat_init(D, n, n);
    acb_mat_init(XY, n, n);
    acb_mat_init(YX, n, n);

    acb_mat_sqr(Y2, Y, prec);
    acb_mat_one(D);
    acb_mat_sub(D, Y2, D, prec);          /* Y^2 - I */
    aic_mat_frobenius_norm(cn, D, prec);
    if (!arb_lt(cn, cert_tol)) {
        fprintf(stderr, "%s: a-posteriori ||Y^2-I|| not < tol "
                "(result is not a certified involution; out-of-basin input or "
                "raise prec)\n", who);
        abort();
    }
    acb_mat_mul(XY, X, Y, prec);
    acb_mat_mul(YX, Y, X, prec);
    acb_mat_sub(D, XY, YX, prec);          /* XY - YX */
    aic_mat_frobenius_norm(cn, D, prec);
    if (!arb_lt(cn, cert_tol)) {
        fprintf(stderr, "%s: a-posteriori ||XY-YX|| not < tol "
                "(result does not commute with X; out-of-basin input or raise "
                "prec)\n", who);
        abort();
    }

    acb_mat_clear(YX);
    acb_mat_clear(XY);
    acb_mat_clear(D);
    acb_mat_clear(Y2);
    arb_clear(sqn);
    arb_clear(cn);
    arb_clear(cert_tol);
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
