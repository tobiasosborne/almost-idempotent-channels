/* test_contraction_sgn.h — the paper cross-check for the contraction solver
 * (beads aic-n3n, "T-contraction"): solve the sgn defining system of
 * approximate_algebras.tex:594-599 and cross-check against aic_sgn.
 *
 * approximate_algebras.tex:596-599 — the defining map and its solution:
 *     f(X,Y1,Y2) = (Y1 Y2 - I,  X Y2 - Y1 X);
 *     f(X,Y1,Y2) = 0  has the solution  Y1 = Y2 = sgn(X).
 * approximate_algebras.tex:606-614 — the Jacobian in (Y1,Y2) and its inverse:
 *     d f / d(Y1,Y2) = [[R_{Y2}, L_{Y1}], [-R_X, L_X]];
 *     at Y1=Y2=sgn(X) all factors commute and the inverse is
 *     [[L_X, -L_{sgn}], [R_X, R_{sgn}]] (L_X R_{sgn} + L_{sgn} R_X)^{-1}.
 *
 * SIMPLER INSTANCE (stated per the task): we take X DIAGONAL with real entries
 * x_i (|x_i^2 - 1| < 1, so sgn(x_i) = +/-1). Then |X|, sgn(X), Y1, Y2 are all
 * diagonal and the operators L_A, R_A reduce to scalar multiplication, so the
 * full operator Jacobian (an operator on B x B) collapses to one 2x2 block per
 * diagonal entry. The Banach element (Y1,Y2) is stored as a 2n x 1 column:
 * rows [0,n) are diag(Y1), rows [n,2n) are diag(Y2). Per entry i:
 *     f: A_i = y1_i y2_i - 1,   B_i = x_i (y2_i - y1_i).
 *     J_i = [[y2_i, y1_i], [-x_i, x_i]];  at y1=y2=s_i: det = 2 s_i x_i = 2|x_i|,
 *     J_i^{-1} = (1/2|x_i|) [[x_i, -s_i], [x_i, s_i]]   (matches tex:612 scalar-
 *     reduced: L_X R_{sgn}+L_{sgn} R_X = 2|x_i|). V := J|_{solution}, fixed.
 * This is the lemma applied to the sgn system, restricted to where the operators
 * are scalars; the recovered Y1=Y2 is cross-checked against aic_sgn(X).
 */
#ifndef TEST_CONTRACTION_SGN_H
#define TEST_CONTRACTION_SGN_H

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_contraction.h"
#include "aic_funcalc.h"
#include "aic_mat.h"
#include "aic_test.h"

typedef struct {
    slong   n;       /* number of diagonal entries           */
    double *x;       /* the diagonal entries x_i (real)       */
    double *s;       /* sgn(x_i) = +/-1                       */
} sgn_ctx;

/* z = x * d, d a double (FLINT 3.0.1 has no acb_mul_d; route via an arb_t). */
static void acb_mul_dbl(acb_t z, const acb_t x, double d, slong prec)
{
    arb_t c;
    arb_init(c);
    arb_set_d(c, d);
    acb_mul_arb(z, x, c, prec);
    arb_clear(c);
}

/* f on the 2n-stacked variable: out rows [0,n) = A_i, rows [n,2n) = B_i. */
static void sgn_f(acb_mat_t out, const acb_mat_t v, void *ctx, slong prec)
{
    sgn_ctx *p = ctx;
    slong n = p->n;
    acb_t y1, y2, t;
    acb_init(y1); acb_init(y2); acb_init(t);
    for (slong i = 0; i < n; i++) {
        acb_set(y1, acb_mat_entry(v, i, 0));
        acb_set(y2, acb_mat_entry(v, n + i, 0));
        acb_mul(t, y1, y2, prec);                       /* y1 y2     */
        acb_sub_si(acb_mat_entry(out, i, 0), t, 1, prec); /* A = y1y2-1 */
        acb_sub(t, y2, y1, prec);                       /* y2 - y1   */
        acb_mul_dbl(acb_mat_entry(out, n + i, 0), t, p->x[i], prec); /* B = x(y2-y1) */
    }
    acb_clear(t); acb_clear(y2); acb_clear(y1);
}

/* V^{-1} = J^{-1} at the solution, block-diagonal: per entry i, with input
 * (a, b) = (r_i, r_{n+i}), output (1/2|x|)[[x,-s],[x,s]](a,b). */
static void sgn_Vinv(acb_mat_t out, const acb_mat_t r, void *ctx, slong prec)
{
    sgn_ctx *p = ctx;
    slong n = p->n;
    acb_t a, b, o1, o2, tmp;
    acb_init(a); acb_init(b); acb_init(o1); acb_init(o2); acb_init(tmp);
    for (slong i = 0; i < n; i++) {
        double x = p->x[i], s = p->s[i];
        double inv = 1.0 / (2.0 * (x * s)); /* 2|x| = 2 x s (s=sgn x) */
        acb_set(a, acb_mat_entry(r, i, 0));
        acb_set(b, acb_mat_entry(r, n + i, 0));
        /* o1 = (x a - s b) * inv */
        acb_mul_dbl(o1, a, x, prec);
        acb_mul_dbl(tmp, b, s, prec);
        acb_sub(o1, o1, tmp, prec);
        acb_mul_dbl(acb_mat_entry(out, i, 0), o1, inv, prec);
        /* o2 = (x a + s b) * inv */
        acb_mul_dbl(o2, a, x, prec);
        acb_mul_dbl(tmp, b, s, prec);
        acb_add(o2, o2, tmp, prec);
        acb_mul_dbl(acb_mat_entry(out, n + i, 0), o2, inv, prec);
    }
    acb_clear(tmp); acb_clear(o2); acb_clear(o1); acb_clear(b); acb_clear(a);
}

static void test_sgn_system(slong prec)
{
    const slong n = 4;
    /* diagonal x_i in (0,sqrt2) and (-sqrt2,0): |x_i^2-1|<1 so sgn=+/-1. */
    double xs[4] = {0.85, 1.15, -0.9, -1.1};
    sgn_ctx p;
    p.n = n; p.x = xs;
    double ss[4];
    for (slong i = 0; i < n; i++) ss[i] = (xs[i] > 0) ? 1.0 : -1.0;
    p.s = ss;

    acb_mat_t v0, y, v;
    acb_mat_init(v0, 2 * n, 1);
    acb_mat_init(y, 2 * n, 1);
    acb_mat_init(v, 2 * n, 1);
    acb_mat_zero(y);                       /* solve f = 0 */
    for (slong i = 0; i < n; i++) {        /* start at Y1=Y2=X (sgn(X)~X) */
        acb_set_d(acb_mat_entry(v0, i, 0), xs[i]);
        acb_set_d(acb_mat_entry(v0, n + i, 0), xs[i]);
    }

    aic_contraction_opts o = {
        .f = sgn_f, .apply_Vinv = sgn_Vinv, .ctx = &p,
        .x0 = v0, .y = y, .r = 0.8, .c = 0.7, .tol = 1e-30,
        .max_iter = 300, .prec = prec
    };
    aic_contraction_stats st;
    aic_contraction_solve(v, &o, &st);

    /* Cross-check: the recovered Y1 (and Y2) equal sgn(X) = diag(s_i). Compare
     * to aic_sgn on the diagonal matrix X (the independent funcalc oracle). */
    acb_mat_t X, sgnX;
    acb_mat_init(X, n, n);
    acb_mat_init(sgnX, n, n);
    acb_mat_zero(X);
    for (slong i = 0; i < n; i++)
        acb_set_d(acb_mat_entry(X, i, i), xs[i]);
    aic_sgn(sgnX, X, prec);

    arb_t tol; arb_init(tol); set_tol(tol, 1e-25);
    for (slong i = 0; i < n; i++) {
        /* Y1_i == Y2_i (the system forces them equal here) */
        AIC_CHECK_ACB_CLOSE(acb_mat_entry(v, i, 0),
                            acb_mat_entry(v, n + i, 0), tol);
        /* Y1_i == sgn(X)_ii (the funcalc oracle) */
        AIC_CHECK_ACB_CLOSE(acb_mat_entry(v, i, 0),
                            acb_mat_entry(sgnX, i, i), tol);
    }
    arb_clear(tol);

    acb_mat_clear(sgnX); acb_mat_clear(X);
    acb_mat_clear(v); acb_mat_clear(y); acb_mat_clear(v0);
}

#endif /* TEST_CONTRACTION_SGN_H */
