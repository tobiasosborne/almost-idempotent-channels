/* aic_adversarial_kahan.c — adversarial instance generator: the Kahan
 * catastrophic-cancellation matrix (bead aic-dbo.2, nla 5b). Split out of
 * aic_adversarial_nla.c because that file is at the ~200 LOC cap (Rule 10); this
 * single generator gets its own translation unit, auto-globbed into
 * aic_testsupport by tests/CMakeLists.txt (the aic_*.c glob). See
 * aic_adversarial.h for the corpus rationale and the knob -> property contract.
 *
 * WHAT THIS MATRIX IS (docs/adversarial/nla.md:619-664, "5b Kahan Matrix"). The
 * Kahan matrix is the textbook example that QR factorization WITH COLUMN PIVOTING
 * is NOT a fail-safe rank-revealer: it is invariant under column-pivoting QR, so
 * the pivoting strategy never reorders it, and the smallest pivot (= the smallest
 * diagonal entry c^{n-1}) OVERESTIMATES the true smallest singular value, which is
 * exponentially smaller. That exponential gap is the catastrophic cancellation
 * that a `double` rank/SVD routine reports wrongly while a certified arb SVD
 * resolves — the arb-justifying instance.
 *
 * WHICH CONVENTION (paper/FINDINGS.md §C21 — a catalog convention bug). The
 * nla.md:621-628 "Math" block writes K[i][i] = s^i, K[i][j] = -c s^{i+j-1}; that
 * sin-graded form does NOT exhibit the pathology (MEASURED: sigma_min ~ smallest
 * diagonal, ratio ~1, no separation). The nla.md PROSE knob ("theta near pi/2,
 * s~0 is lethal") and the classic kappa ~ (2c/s^{n-1})^n both describe the MATLAB
 * gallery('kahan') / Higham convention instead:
 *   c = cos(theta), s = sin(theta),
 *   K[i][i] = c^i,   K[i][j] = -s * c^i  (j > i),   0 below the diagonal.
 * We constructivize the genuinely-pathological gallery form (Law 3: the artifact
 * is the matrix that actually breaks the routine, not a literal transcription of a
 * catalog typo), and record the discrepancy in FINDINGS §C21.
 *
 * KNOB DIRECTION. theta in (0, pi/2). theta -> 0 gives c -> 1, s -> 0, K -> I (the
 * well-conditioned reduction oracle: kappa -> 1, sigma_min -> 1). theta -> pi/2
 * gives c -> 0: the diagonal grading c^i collapses, sigma_min plummets
 * exponentially, kappa explodes — the LETHAL end. So SMALL theta is MILD, LARGE
 * theta (toward pi/2) is LETHAL, matching nla.md's prose.
 *
 * PRECISION / DETERMINISM. c^i is a certified arb power (arb_pow_ui) so the
 * exponentially-small diagonal grading is exact to `prec`, not contaminated by
 * double rounding of c^i for large i. Closed form, no RNG: bit-identical per knob.
 *
 * Measured anchors (prec=256), asserted by tests/test_adversarial.c:
 *   n=6 theta=0.3 (mild)  : sigma_min=0.353, kappa=3.56,   sigma_min/c^5=0.44
 *   n=6 theta=1.4 (lethal): sigma_min=7.96e-6, kappa=3.05e5, sigma_min/c^5=0.056
 *   (the lethal sigma_min is ~18x BELOW the smallest diagonal c^5=1.42e-4 — the
 *    catastrophic cancellation a column-pivoting QR cannot see.)
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"

/* gen-kahan — Kahan matrix, gallery('kahan') convention (see file docstring,
 * FINDINGS §C21). K[i][i] = c^i; K[i][j] = -s c^i for j > i; 0 for j < i. The
 * diagonal c^i is built incrementally (ci *= c) from a certified arb c, and the
 * row's off-diagonal value -s*c^i is the same ci scaled by -s, reused across the
 * row (every super-diagonal entry in row i shares the factor c^i). */
void aic_adv_kahan(acb_mat_t out, slong n, double theta, slong prec)
{
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_adv_kahan: out must be n x n");
    assert(n >= 2 && "aic_adv_kahan: need n >= 2");
    /* pi/2 = 1.5707963... ; explicit literal avoids the M_PI feature-test macro
     * under -std=c11 -Wpedantic (M_PI is not ISO C). */
    assert(theta > 0.0 && theta < 1.5707963267948966 &&
           "aic_adv_kahan: need theta in (0, pi/2) (c=cos, s=sin both positive; "
           "theta=0 is the identity reduction, theta=pi/2 is the rank-1 singular end)");
    acb_mat_zero(out);
    arb_t c, ci, off, negs;
    arb_init(c);
    arb_init(ci);
    arb_init(off);
    arb_init(negs);
    arb_set_d(c, cos(theta));    /* c = cos(theta), certified ball */
    arb_set_d(negs, -sin(theta)); /* -s = -sin(theta) */
    arb_one(ci);                 /* ci = c^0 = 1 */
    for (slong i = 0; i < n; i++) {
        /* K[i][i] = c^i */
        arb_set(acb_realref(acb_mat_entry(out, i, i)), ci);
        /* K[i][j] = -s * c^i  for all j > i (constant across the row) */
        arb_mul(off, ci, negs, prec);
        for (slong j = i + 1; j < n; j++)
            arb_set(acb_realref(acb_mat_entry(out, i, j)), off);
        arb_mul(ci, ci, c, prec); /* advance c^i -> c^{i+1} */
    }
    arb_clear(negs);
    arb_clear(off);
    arb_clear(ci);
    arb_clear(c);
}
