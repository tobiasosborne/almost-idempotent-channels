/* aic_contraction_work.c — shared scratch + guards for the contraction solver
 * (beads aic-n3n, "T-contraction"). See src/aic_contraction_internal.h.
 *
 * approximate_algebras.tex:580-591 — the proof's structural facts, enforced here
 * as loud invariants (CLAUDE.md Rule 4):
 *   g_y(x) = x + V^{-1}(y - f(x));   ||d g_y|| <= c < 1;
 *   if ||x-x0|| < r then ||g_y(x)-x0|| < c||x-x0|| + (1-c)r <= r  (self-map);
 *   the fixed point is lim x_n=g_y(x_{n-1}), Cauchy because
 *   ||x_n - x_{n-1}|| < r c^{n-1}.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic_contraction_internal.h"

void aic_contraction_work_init(aic_contraction_work *w, const acb_mat_t x0)
{
    w->nr = acb_mat_nrows(x0);
    w->nc = acb_mat_ncols(x0);
    acb_mat_init(w->fx, w->nr, w->nc);
    acb_mat_init(w->resid, w->nr, w->nc);
    acb_mat_init(w->step, w->nr, w->nc);
    acb_mat_init(w->x0mid, w->nr, w->nc);
    acb_mat_get_mid(w->x0mid, x0);
    arb_init(w->one);
    arb_one(w->one);
}

void aic_contraction_work_clear(aic_contraction_work *w)
{
    arb_clear(w->one);
    acb_mat_clear(w->x0mid);
    acb_mat_clear(w->step);
    acb_mat_clear(w->resid);
    acb_mat_clear(w->fx);
}

void aic_contraction_assert_c(double c)
{
    if (!(c >= 0.0 && c < 1.0)) {
        fprintf(stderr,
                "aic_contraction: hypothesis violated, c = %.17g is not in "
                "[0,1) (tex:569: ||V^{-1} d_x f - 1|| <= c < 1)\n", c);
        abort();
    }
}

void aic_contraction_gstep(aic_contraction_work *w, const acb_mat_t x,
                           const aic_contraction_opts *o)
{
    /* g_y(x) - x = V^{-1}(y - f(x))  (tex:582). */
    o->f(w->fx, x, o->ctx, o->prec);                 /* f(x)        */
    acb_mat_sub(w->resid, o->y, w->fx, o->prec);    /* y - f(x)    */
    o->apply_Vinv(w->step, w->resid, o->ctx, o->prec); /* V^{-1}(.) */
}

void aic_contraction_mid_opnorm(arb_t out, const acb_mat_t A, const acb_mat_t B,
                                slong prec)
{
    slong nr = acb_mat_nrows(A), nc = acb_mat_ncols(A);
    acb_mat_t D;
    acb_mat_init(D, nr, nc);
    acb_mat_sub(D, A, B, prec);
    acb_mat_get_mid(D, D);          /* strip radii: honest "iterate moved" */
    aic_mat_opnorm(out, D, prec);   /* certified op-norm of the midpoint diff */
    acb_mat_clear(D);
}

void aic_contraction_assert_in_ball(aic_contraction_work *w, const acb_mat_t x,
                                    double r, slong prec)
{
    arb_t dist, rad;
    arb_init(dist);
    arb_init(rad);
    aic_contraction_mid_opnorm(dist, x, w->x0mid, prec); /* ||x - x0|| */
    arb_set_d(rad, r);
    if (!arb_le(dist, rad)) {
        fprintf(stderr,
                "aic_contraction: iterate left the ball B_r(x0), r = %.17g "
                "(tex:586-590 self-map property failed: a hypothesis is "
                "wrong, e.g. y too far from f(x0) or c underestimated)\n", r);
        abort();
    }
    arb_clear(rad);
    arb_clear(dist);
}

void aic_contraction_assert_cauchy(const arb_t step_norm, double r, double c,
                                   slong n, slong prec)
{
    /* bound = r * c^(n-1) (tex:591). For n=1 this is r (c^0); the first step is
     * bounded by r because g_y maps into B_r(x0). A small slack factor (1 + 2^-20)
     * absorbs the certified op-norm's own rounding so a step exactly on the
     * geometric law is not spuriously rejected; the bound stays an O(1) envelope
     * of the paper's r c^{n-1}, not a loosened claim. */
    arb_t bound, cc, slack;
    arb_init(bound);
    arb_init(cc);
    arb_init(slack);
    arb_set_d(cc, c);
    arb_pow_ui(cc, cc, (ulong)(n - 1), prec);  /* c^(n-1) */
    arb_set_d(bound, r);
    arb_mul(bound, bound, cc, prec);           /* r c^(n-1) */
    arb_set_d(slack, 1.0 + 9.5367431640625e-07); /* 1 + 2^-20 */
    arb_mul(bound, bound, slack, prec);
    if (!arb_le(step_norm, bound)) {
        char *sn = arb_get_str(step_norm, 10, 0);
        char *bn = arb_get_str(bound, 10, 0);
        fprintf(stderr,
                "aic_contraction: Cauchy bound violated at step n=%ld: "
                "||x_n - x_{n-1}|| = %s > r c^{n-1} = %s (tex:591). "
                "c is underestimated for this f.\n", (long) n, sn, bn);
        flint_free(sn);
        flint_free(bn);
        abort();
    }
    arb_clear(slack);
    arb_clear(cc);
    arb_clear(bound);
}

void aic_contraction_record(aic_contraction_stats *st, const arb_t step,
                            slong n, slong prec)
{
    if (!st) return;
    arf_t u;
    arf_init(u);
    arb_get_ubound_arf(u, step, prec);
    double s = arf_get_d(u, ARF_RND_UP);
    arf_clear(u);
    if (n == 1 || s > st->max_step) st->max_step = s;
    st->last_step = s;
    st->iters = n;
}
