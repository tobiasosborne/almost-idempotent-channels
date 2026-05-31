/* aic_contraction.c — the two audition candidates (Picard, Anderson) and the
 * default dispatcher for the Banach contraction / inverse-function solver
 * (beads aic-n3n, "T-contraction"; CLAUDE.md Law 4: audition, don't presume).
 * See include/aic_contraction.h for the API and full paper provenance
 * (approximate_algebras.tex:564-592, Lemma lem_invfun).
 *
 * approximate_algebras.tex:580-591 — the fixed-point map and its convergence:
 *     g_y(x) = x + V^{-1}(y - f(x));      ||d g_y|| <= c < 1;
 *     fixed point = lim x_n = g_y(x_{n-1}),  ||x_n - x_{n-1}|| < r c^{n-1}.
 *
 * CANDIDATE 1 — Picard (the paper's). x_{n} = g_y(x_{n-1}) verbatim. Each step
 * is one f, one V^{-1}, one matrix add. The r c^{n-1} Cauchy bound (tex:591) is
 * asserted at every step — both the paper's claim and a guard against a
 * mis-estimated c.
 *
 * CANDIDATE 2 — Anderson(1) (Type-II / Pulay mixing) on the SAME g_y, using the
 * residual r_k = g_y(x_k) - x_k (the Picard step):
 *     theta = <dr, r_k> / <dr, dr>,   dr = r_k - r_{k-1},
 *     x_{k+1} = (x_k + r_k) - theta * ((x_k + r_k) - (x_{k-1} + r_{k-1})).
 * It does NOT obey r c^{n-1} by design (the mixing beats it), so the Cauchy
 * guard is Picard-only; Anderson keeps the c<1 and ball-containment guards.
 * (Full rationale in include/aic_contraction.h.)
 *
 * DEFAULT (aic_contraction_solve) — Picard, set by bench evidence, NOT presumed.
 * bench_contraction.c times both on g_y(x) = x + kappa V^{-1}(y - f(x)) for the
 * known-inverse f(x)=Vx+eps*N(x), sweeping the contraction rate via the
 * preconditioner quality kappa (prec=53; iteration counts are deterministic, the
 * ns/op are means over 4 runs on a loaded box, +/-10%):
 *   FAST regime kappa=1, c_true~0.03: Picard 9 iters / Anderson 8. ns/op is a
 *     near-tie and SIZE-DEPENDENT: at n=4 Anderson edges it (~89us vs ~97us), at
 *     n=8 Picard wins (~134us vs ~159us) — Anderson's 1 saved iteration does not
 *     repay its per-step dg + two O(n^2) inner-products as n grows. Within noise.
 *   SLOW regime kappa=0.5, c_true~0.5: Picard 38-40 iters / Anderson 24 —
 *     Anderson cuts ~40% of iterations and ROBUSTLY wins ns/op at both sizes
 *     (n=4: ~280us vs ~390us; n=8: ~460us vs ~630us).
 * Neither dominates: the winner is the contraction rate. The lemma is applied in
 * this project with V the actual Jacobian (the FAST, small-c regime), where the
 * two tie within noise and Picard is faithful to tex:591 and wins at the larger
 * size — so Picard is the default, with Anderson exposed for the stiff
 * (mediocre-preconditioner) regime where it robustly wins. See the BENCH/ITERS
 * lines in the bead aic-n3n report.
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

void aic_contraction_picard(acb_mat_t x, const aic_contraction_opts *o,
                            aic_contraction_stats *st)
{
    aic_contraction_assert_c(o->c);
    aic_contraction_work w;
    aic_contraction_work_init(&w, o->x0);

    arb_t tol, step;
    arb_init(tol);
    arb_init(step);
    arb_set_d(tol, o->tol);
    acb_mat_set(x, o->x0);                  /* x_0 = centre */
    if (st) { st->iters = 0; st->last_step = 0.0; st->max_step = 0.0; }

    acb_mat_t zero;
    acb_mat_init(zero, w.nr, w.nc);          /* fixed 0 for the step-norm */

    slong n;
    for (n = 1; n <= o->max_iter; n++) {
        aic_contraction_gstep(&w, x, o);     /* w.step = g_y(x) - x */
        aic_contraction_mid_opnorm(step, w.step, zero, o->prec); /* ||step|| */
        acb_mat_add(x, x, w.step, o->prec);  /* x_n = g_y(x_{n-1}) */
        aic_contraction_assert_in_ball(&w, x, o->r, o->prec); /* tex:586-590 */
        aic_contraction_assert_cauchy(step, o->r, o->c, n, o->prec); /* tex:591 */
        aic_contraction_record(st, step, n, o->prec);
        if (arb_lt(step, tol)) break;
    }
    if (n > o->max_iter) {
        fprintf(stderr, "aic_contraction_picard: no convergence in %ld iters "
                "(c<1 must converge at rate c; c mis-estimated or ill-posed)\n",
                (long) o->max_iter);
        abort();
    }
    acb_mat_clear(zero);
    arb_clear(step);
    arb_clear(tol);
    aic_contraction_work_clear(&w);
}

/* aic_contraction_anderson lives in src/aic_contraction_anderson.c (kept under
 * the ~200-LOC limit, CLAUDE.md Rule 10). */

/* Default: Picard (the paper's method), chosen on bench evidence (file
 * docstring; bench/bench_contraction.c). Lower per-step cost wins in the common
 * small-c regime; Anderson stays available for c near 1. */
void aic_contraction_solve(acb_mat_t x, const aic_contraction_opts *o,
                           aic_contraction_stats *st)
{
    aic_contraction_picard(x, o, st);
}
