/* aic_funcalc_internal.h — shared internals for the funcalc module (beads
 * aic-4bg, "T-funcalc"). NOT part of the public API (include/aic_funcalc.h).
 * The validity-domain precondition ||X^2 - I|| < 1 (tex:516, tex:520-521) is the
 * single guard shared by both sgn candidates, by aic_abs, and by aic_theta, so
 * it is written once here and asserted at every entry point (CLAUDE.md Rule 4:
 * fail fast, fail loud — an out-of-domain functional-calculus call is a bug).
 */
#ifndef AIC_FUNCALC_INTERNAL_H
#define AIC_FUNCALC_INTERNAL_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

/* ||X^2 - I||_op as a certified arb ball into `out` (operator norm = largest
 * singular value, via aic_mat_opnorm). X must be square (asserted). */
void aic_funcalc_int_def_X2(arb_t out, const acb_mat_t X, slong prec);

/* Abort (fail loud) unless ||X^2 - I||_op is CERTAINLY < 1, the validity domain
 * of |X|, sgn(X), theta(X) (tex:516, tex:520). Returns the certified bound on
 * ||X^2 - I|| in *bound (caller-initialised arb_t) for the iteration's tol. */
void aic_funcalc_int_assert_domain(arb_t bound, const acb_mat_t X, slong prec);

/* prec-appropriate convergence tol = 2^-(prec - 8); matches aic_mat_int_tol so
 * the funcalc iterations and the mat substrate use one notion of "resolved". */
void aic_funcalc_int_tol(arb_t tol, slong prec);

#endif /* AIC_FUNCALC_INTERNAL_H */
