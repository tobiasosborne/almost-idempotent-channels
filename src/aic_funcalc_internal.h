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

/* MIDPOINT Frobenius step ||A - B||_F into `out`, on the ball MIDPOINTS (radii
 * stripped). Shared by every sgn iteration (Newton-Schulz, Denman-Beavers, global
 * Newton) as the honest "the iterate has stopped moving" convergence measure: a
 * quadratic interval iteration reaches the fixed point in a handful of steps, but
 * the certified ball RADII keep inflating, so a ball-norm convergence test stalls.
 * The returned matrix stays a rigorous enclosure; we just stop once the midpoint
 * step is below tol. Frobenius (no Gram/eig) so it is robust to non-Hermitian
 * iteration matrices. `out` must be an initialised arb_t. */
void aic_funcalc_int_step_norm(arb_t out, const acb_mat_t A, const acb_mat_t B,
                               slong prec);

/* NON-aborting op-norm-basin probe (bead aic-8hz). Returns 1 iff ||X^2 - I||_op
 * is CERTAINLY < 1 (the Newton-Schulz local-convergence basin; same certified
 * ball as aic_funcalc_int_assert_domain, but a verdict, not an abort). A
 * straddling ball returns 0 (not certified). aic_sgn uses this to dispatch:
 * in-basin -> fast inverse-free Newton-Schulz; out-of-basin -> global Newton. */
int aic_funcalc_int_in_op_basin(const acb_mat_t X, slong prec);

/* NON-aborting Gelfand spectral-radius certifier (bead aic-8hz). Given M,
 * tries k = 1, 2, ..., k_max and accepts the FIRST k whose certified arb bound
 * ||M^k||_F^{1/k} is CERTAINLY < 1 (Gelfand: ||M^k||^{1/k} >= rho(M) for every
 * k and -> rho(M), so a single certified k < 1 proves rho(M) < 1, eig-free).
 * On success returns 1, writes that bound into `rho_ub` and the accepting k into
 * *k_used (if non-NULL). On failure (no k <= k_max certifies; ball straddles 1)
 * returns 0 and `rho_ub` holds the tightest (largest-k) bound tried. Frobenius
 * norm (no SVD); M^k built incrementally; arb_root_ui for the k-th root. The
 * public global-Newton wrapper aborts on a 0 return (Rule 4); the test drives
 * this internal entry directly to assert the not-certified status (T-global-4). */
int aic_funcalc_int_gelfand_rho(arb_t rho_ub, slong *k_used, const acb_mat_t M,
                                slong k_max, slong prec);

#endif /* AIC_FUNCALC_INTERNAL_H */
