/* test_harness.c — red-green self-test of the audition harness (beads aic-73v,
 * "T-harness"). The harness is itself auditioned (CLAUDE.md Law 4): before any
 * module trusts aic_acb_close / aic_acb_mat_close to certify the temporary
 * acb@53-vs-acb@high-prec ladder (beads aic-5ty), this test proves the checker
 * answers TRUE for balls that agree and FALSE for balls that plainly disagree.
 *
 * It cannot use AIC_CHECK_ACB_CLOSE for the negative case (that macro aborts on
 * a false verdict, which is the whole point); so it calls aic_acb_close
 * directly and asserts the boolean return with AIC_CHECK. That is the red-green
 * proof: the same predicate returns 1 then 0 under controlled inputs.
 *
 * Per CLAUDE.md Rule 5 every check asserts a value. No paper provenance: this
 * tests the harness, not a paper result.
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_test.h"

int main(void)
{
    arb_t tol;
    arb_init(tol);
    arb_set_d(tol, 1e-9); /* "agree within ~1e-9", the ladder's working tol */

    /* --- scalar: equal balls are close, different balls are not --- */
    acb_t a, b, far;
    acb_init(a);
    acb_init(b);
    acb_init(far);

    acb_set_si_si(a, 3, -7); /* 3 - 7i, exact (zero-radius) ball */
    acb_set_si_si(b, 3, -7); /* identical */
    acb_set_si_si(far, 3, -6); /* differs by exactly 1 in the imag part */

    /* GREEN: identical balls agree within tol. */
    AIC_CHECK(aic_acb_close(a, b, tol) == 1);

    /* RED→asserted: a ball 1.0 away is NOT within 1e-9. */
    AIC_CHECK(aic_acb_close(a, far, tol) == 0);

    /* A near-but-not-equal pair just inside / outside tol, to prove the
     * boundary is the tol and not a midpoint heuristic. */
    acb_t a_eps;
    acb_init(a_eps);
    acb_set_d_d(a_eps, 3.0, -7.0 + 1e-12); /* off by 1e-12 << 1e-9 */
    AIC_CHECK(aic_acb_close(a, a_eps, tol) == 1);
    acb_set_d_d(a_eps, 3.0, -7.0 + 1e-6);  /* off by 1e-6 >> 1e-9 */
    AIC_CHECK(aic_acb_close(a, a_eps, tol) == 0);

    /* --- matrix: equal matrices agree, a single perturbed entry fails --- */
    const slong dim = 3;
    acb_mat_t A, B;
    acb_mat_init(A, dim, dim);
    acb_mat_init(B, dim, dim);
    acb_mat_one(A);
    acb_mat_one(B);

    slong fi = -1, fj = -1;
    AIC_CHECK(aic_acb_mat_close(A, B, tol, &fi, &fj) == 1);

    /* Perturb B[1][2] by 1.0; the checker must report FALSE at (1,2). */
    acb_set_d(acb_mat_entry(B, 1, 2), 1.0);
    AIC_CHECK(aic_acb_mat_close(A, B, tol, &fi, &fj) == 0);
    AIC_CHECK_MSG(fi == 1 && fj == 2,
                  "expected first mismatch at (1,2), got (%ld,%ld)",
                  (long) fi, (long) fj);

    /* Exercise the abort-on-failure macros on the passing direction (these
     * would abort the binary if they regressed). */
    AIC_CHECK_ACB_CLOSE(a, b, tol);
    AIC_CHECK_ACB_MAT_CLOSE(A, A, tol);

    acb_mat_clear(B);
    acb_mat_clear(A);
    acb_clear(a_eps);
    acb_clear(far);
    acb_clear(b);
    acb_clear(a);
    arb_clear(tol);

    aic_test_report("test_harness");
    return 0;
}
