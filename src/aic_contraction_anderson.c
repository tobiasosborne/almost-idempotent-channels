/* aic_contraction_anderson.c — the accelerated audition candidate (beads
 * aic-n3n, "T-contraction"; CLAUDE.md Law 4). Anderson(1) / Type-II (Pulay)
 * mixing on the SAME fixed-point map g_y(x) = x + V^{-1}(y - f(x)) as Picard
 * (approximate_algebras.tex:580-592); see include/aic_contraction.h for the full
 * candidate rationale and the Picard-vs-Anderson selection evidence.
 *
 * With the residual r_k = g_y(x_k) - x_k (the Picard step), Anderson(1) mixes the
 * last two residuals over the real inner product to damp the dominant linear
 * mode (rate c):
 *     theta = <dr, r_k> / <dr, dr>,   dr = r_k - r_{k-1},
 *     x_{k+1} = g_k - theta * (g_k - g_{k-1}),     g_k = x_k + r_k = g_y(x_k).
 * It does NOT obey the r c^{n-1} Cauchy law (the mixing beats it), so the Cauchy
 * guard (Picard-only) is omitted; the c<1 and ball-containment guards (tex:569,
 * tex:586-590) are kept — a c >= 1 or an iterate leaving B_r(x0) is still a
 * fail-loud abort.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic/aic_contraction.h"
#include "aic_contraction_internal.h"

/* Real inner product Re Tr(A^dag B) (midpoints) as a double, used only to form
 * Anderson's scalar mixing coefficient; the coefficient is then applied with
 * certified acb arithmetic, so the iterate remains a rigorous enclosure. */
static double real_dot(const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong nr = acb_mat_nrows(A), nc = acb_mat_ncols(A);
    arb_t acc, t;
    arb_init(acc);
    arb_init(t);
    arb_zero(acc);
    for (slong i = 0; i < nr; i++)
        for (slong j = 0; j < nc; j++) {
            const acb_struct *a = acb_mat_entry(A, i, j);
            const acb_struct *b = acb_mat_entry(B, i, j);
            arb_mul(t, acb_realref(a), acb_realref(b), prec);
            arb_add(acc, acc, t, prec);
            arb_mul(t, acb_imagref(a), acb_imagref(b), prec);
            arb_add(acc, acc, t, prec);
        }
    double d = arf_get_d(arb_midref(acc), ARF_RND_NEAR);
    arb_clear(t);
    arb_clear(acc);
    return d;
}

void aic_contraction_anderson(acb_mat_t x, const aic_contraction_opts *o,
                              aic_contraction_stats *st)
{
    aic_contraction_assert_c(o->c);
    aic_contraction_work w;
    aic_contraction_work_init(&w, o->x0);

    arb_t tol, step;
    arb_init(tol);
    arb_init(step);
    arb_set_d(tol, o->tol);

    acb_mat_t rprev, gprev, dr, dg, zero; /* prev residual/g-value, deltas, 0 */
    acb_mat_init(rprev, w.nr, w.nc);
    acb_mat_init(gprev, w.nr, w.nc);
    acb_mat_init(dr, w.nr, w.nc);
    acb_mat_init(dg, w.nr, w.nc);
    acb_mat_init(zero, w.nr, w.nc);
    acb_mat_set(x, o->x0);
    if (st) { st->iters = 0; st->last_step = 0.0; st->max_step = 0.0; }

    slong n;
    for (n = 1; n <= o->max_iter; n++) {
        aic_contraction_gstep(&w, x, o);            /* r_k = g(x)-x in w.step */
        aic_contraction_mid_opnorm(step, w.step, zero, o->prec);

        acb_mat_t gx;
        acb_mat_init(gx, w.nr, w.nc);
        acb_mat_add(gx, x, w.step, o->prec);        /* g(x_k) = x + r_k */

        if (n == 1) {
            acb_mat_set(x, gx);                     /* first step = Picard */
        } else {
            acb_mat_sub(dr, w.step, rprev, o->prec); /* dr = r_k - r_{k-1} */
            acb_mat_sub(dg, gx, gprev, o->prec);     /* dg = g_k - g_{k-1} */
            double num = real_dot(dr, w.step, o->prec);
            double den = real_dot(dr, dr, o->prec);
            double theta = (den > 0.0) ? num / den : 0.0;
            acb_t th;                                /* x_{k+1} = g_k - theta dg */
            acb_init(th);
            acb_set_d(th, theta);
            acb_mat_scalar_mul_acb(dg, dg, th, o->prec);
            acb_mat_sub(x, gx, dg, o->prec);
            acb_clear(th);
        }
        acb_mat_set(rprev, w.step);
        acb_mat_set(gprev, gx);
        acb_mat_clear(gx);

        aic_contraction_assert_in_ball(&w, x, o->r, o->prec); /* tex:586-590 */
        aic_contraction_record(st, step, n, o->prec);
        if (arb_lt(step, tol)) break;
    }
    if (n > o->max_iter) {
        fprintf(stderr, "aic_contraction_anderson: no convergence in %ld "
                "iters\n", (long) o->max_iter);
        abort();
    }
    acb_mat_clear(zero);
    acb_mat_clear(dg);
    acb_mat_clear(dr);
    acb_mat_clear(gprev);
    acb_mat_clear(rprev);
    arb_clear(step);
    arb_clear(tol);
    aic_contraction_work_clear(&w);
}
