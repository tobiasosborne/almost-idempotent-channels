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
#include <flint/mag.h>

/* ||X^2 - I||_op as a certified arb ball into `out` (operator norm = largest
 * singular value, via aic_mat_opnorm). X must be square (asserted). */
void aic_funcalc_int_def_X2(arb_t out, const acb_mat_t X, slong prec);

/* MIDPOINT copy: write mid(X) (every entry's radius stripped) into `out`, and
 * return the LARGEST entrywise radius of X (max over re/im of every entry) in
 * `in_rad` (caller-initialised mag_t). The wide-radius sgn fix (bead aic-1vp,
 * FINDINGS §C11) iterates the Newton-Schulz / global-Newton flow on this
 * midpoint: an inherited input radius (~1e-70 from the upstream corner star-
 * product matmul chain) propagates through acb_mat_sqr/mul and DRAGS THE MIDPOINT
 * off the involution once it dominates, so the iteration's midpoint step never
 * settles below the prec-tol floor (2^-(prec-8) ~ 2.2e-75 at prec=256) and the
 * cap-abort fires on a genuinely in-basin X. mid(X) is well-conditioned for sgn
 * here (spec(X) near +/-1, away from 0), so sgn(mid X) is sgn(X_true) within the
 * propagated radius; the a-posteriori certificate (below) restores rigor.
 * `out` must be n x n, initialised; `out` may NOT alias `X`. */
void aic_funcalc_int_to_midpoint(acb_mat_t out, const acb_mat_t X, mag_t in_rad);

/* A-POSTERIORI SIGN CERTIFICATE (shared by global-Newton and Newton-Schulz, bead
 * aic-1vp). After the midpoint iteration, certify that Y is a genuine sign of X:
 *     ||Y^2 - I||_F < tol   AND   ||YX - XY||_F < tol   (certified balls).
 * tol = max(2^-(prec/2), in_rad) * 8 * sqrt(n) — a SANITY BACKSTOP, not a tight
 * enclosure. In the operating regime the prec-floor 2^-(prec/2) DOMINATES (e.g.
 * ~2.9e-39 at prec=256, vs an inherited input radius ~1e-70), so the certificate
 * rejects GROSS failures (a non-converged / wrong Y has ||Y^2-I|| or ||[X,Y]|| =
 * O(1) >> tol) while admitting a converged midpoint sign (||Y^2-I|| ~ the midpoint
 * iteration's machine floor ~1e-74 << tol); the in_rad term only binds if in_rad >
 * 2^-(prec/2). The ||YX-XY|| arm IS computed on the radius-carrying X, so it tests
 * sgn(mid X) against the TRUE input ball. SOUNDNESS does NOT rest on this tol
 * magnitude: it rests on (i) the away-from-0 basin/Gelfand precondition asserted on
 * the radius-carrying X BEFORE to_midpoint (rho(I-X^2)<1 keeps spec(X) off the
 * imaginary axis, so sgn(mid X)=sgn(X_true)), and (ii) the Y0=mid(X) seeding (Higham
 * Thm 5.6: Newton/NS converges to the correct-INERTIA sign — the certificate alone
 * cannot distinguish +sgn from -sgn). A straddling / too-large ball is a FAIL-LOUD
 * abort (Rule 4): a genuinely out-of-basin input must still abort, NOT silently
 * return the midpoint sign. `who` is the caller name; X carries its radius. */
void aic_funcalc_int_certify_sign(const acb_mat_t Y, const acb_mat_t X,
                                  const mag_t in_rad, const char *who,
                                  slong prec);

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
