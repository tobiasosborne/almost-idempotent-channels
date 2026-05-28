/* aic_contraction_internal.h — shared internals for the contraction solver
 * (beads aic-n3n, "T-contraction"). NOT public (include/aic_contraction.h).
 *
 * The g_y step, the ball-containment / Cauchy-bound guards (tex:586-591), and
 * the convergence-tol machinery are common to both audition candidates (Picard
 * and Anderson), so they live here and are reused (CLAUDE.md Law 2), keeping
 * each candidate's source file under the 200-LOC limit (Rule 10).
 */
#ifndef AIC_CONTRACTION_INTERNAL_H
#define AIC_CONTRACTION_INTERNAL_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_contraction.h"

/* Scratch matrices reused across iterations (one set per solve), so a candidate
 * does not re-init/clear acb_mat_t every step. All are x0-shaped. */
typedef struct {
    slong       nr, nc;
    acb_mat_t   fx;     /* f(x)          */
    acb_mat_t   resid;  /* y - f(x)      */
    acb_mat_t   step;   /* V^{-1}(y-f(x)) = g_y(x) - x  (the Picard step)     */
    acb_mat_t   x0mid;  /* midpoint of x0 (for the ball-containment check)    */
    arb_t       one;    /* the constant 1 (for the c<1 assertion)            */
} aic_contraction_work;

/* Init/clear scratch from the centre x0 (records its shape). */
void aic_contraction_work_init(aic_contraction_work *w, const acb_mat_t x0);
void aic_contraction_work_clear(aic_contraction_work *w);

/* Assert (fail loud, tex:569) that c is a real number with 0 <= c < 1. */
void aic_contraction_assert_c(double c);

/* Compute the Picard step s = g_y(x) - x = V^{-1}(y - f(x)) into w->step,
 * using callbacks o->f and o->apply_Vinv. Leaves f(x) in w->fx, residual in
 * w->resid. */
void aic_contraction_gstep(aic_contraction_work *w, const acb_mat_t x,
                           const aic_contraction_opts *o);

/* Certified ||A - B||_op of the ball MIDPOINTS (radii stripped) as an arb ball.
 * Midpoints: an interval iteration's radius inflates after the midpoint has
 * stopped moving (the funcalc module's step_norm finding, bead aic-4bg); the
 * midpoint step is the honest "iterate stopped moving" measure. The returned
 * matrix is still a rigorous enclosure. */
void aic_contraction_mid_opnorm(arb_t out, const acb_mat_t A, const acb_mat_t B,
                                slong prec);

/* Fail loud (tex:586-590) unless ||x - x0||_op (midpoints) is certainly <= r:
 * the proof's self-map property B_r(x0) -> B_r(x0). Aborts otherwise. */
void aic_contraction_assert_in_ball(aic_contraction_work *w, const acb_mat_t x,
                                    double r, slong prec);

/* Fail loud (tex:591) unless step_norm is certainly < r * c^(n-1): the Cauchy
 * bound ||x_n - x_{n-1}|| < r c^{n-1} that proves convergence. n is the 1-based
 * iteration index (so the first step uses c^0 = 1). Only enforced for Picard,
 * whose error obeys exactly this geometric law; Anderson's mixing intentionally
 * beats it, so Anderson does not call this. */
void aic_contraction_assert_cauchy(const arb_t step_norm, double r, double c,
                                   slong n, slong prec);

/* Record one step's midpoint norm into *st (no-op if st==NULL): updates the
 * last/max step (the latter is the input to the test's r c^{n-1} check) and the
 * iteration count. Shared by both candidates. */
void aic_contraction_record(aic_contraction_stats *st, const arb_t step,
                            slong n, slong prec);

#endif /* AIC_CONTRACTION_INTERNAL_H */
