/* aic_test.h — header-only test helpers for the audition harness (beads aic-73v,
 * "T-harness"). Every later module's red-green TDD (CLAUDE.md Law 4: audition,
 * don't presume) drives through these assertion macros and arb/acb agreement
 * checkers.
 *
 * Fail-loud (CLAUDE.md Rule 4): a failed check prints file:line + a concrete
 * message to stderr and abort()s. "Runs without errors" is not a pass (Rule 5):
 * AIC_CHECK_* assert a value/bound, and AIC_CHECKS counts them so a test reports
 * how many invariants it actually exercised.
 *
 * The cross-check ladder, temporarily collapsed (beads aic-5ty).
 *   CLAUDE.md §"cross-check ladder" wants, weakest→strongest:
 *     (1) double-path internal sanity,
 *     (2) double vs arb@prec=53 agreement (~1e-12),
 *     (3) eta=0 exact-idempotent oracle (zero C*-defect),
 *     (4) arb bound certification at the hypothesis boundary.
 *   Rung (2) and the eta=0 oracle are DEFERRED here:
 *     - Rung (2) needs a double/LAPACK fast path. LAPACK has no dev lib/headers
 *       on this box (beads aic-5ty), so there is no double path to compare
 *       against. Until one exists the temporary ladder is acb@53 vs
 *       acb@high-prec: run the same arb routine at prec=53 and at a high prec,
 *       and require the two error balls to agree within an explicit tol. That is
 *       exactly what aic_acb_close / aic_acb_mat_close certify.
 *     - The eta=0 oracle (rung 3) needs the funcalc/ucp/idemp_structure modules
 *       (MODULE_PLAN.md "Suggested first vertical slice"); none exist yet. It is
 *       deferred to those beads, not implemented here.
 *   This header builds the REUSABLE checkers only; no oracle, no LAPACK gate.
 *
 * Agreement predicate (rigorous). aic_acb_close(a,b,tol) returns 1 iff the ball
 * |a-b| is CERTAINLY <= tol: it forms a-b (acb_sub), takes its magnitude as an
 * arb ball (acb_abs, whose result encloses every |a-b| value), and asks
 * arb_le(|a-b|, tol) — which by FLINT semantics is nonzero only when every value
 * of the first ball is <= every value of the second. So a "close" verdict is a
 * rigorous enclosure statement, not a midpoint comparison; a ball that has lost
 * precision (wide |a-b|) fails the test loudly rather than passing by accident.
 *
 * Header-only: no test-framework dependency, no global state beyond one counter
 * (CLAUDE.md "no over-build"). Include after the FLINT headers you need.
 */
#ifndef AIC_TEST_H
#define AIC_TEST_H

#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

/* Count of checks exercised in the current binary (Rule 5: report how many
 * invariants ran, not just that nothing crashed). */
static long aic_test_checks = 0;

/* Fail-loud core: print file:line + message to stderr and abort(). Not a
 * silent return (Rule 4). Used by the macros below; takes a printf-style msg. */
#define AIC_FAIL(...)                                                          \
    do {                                                                       \
        fprintf(stderr, "AIC_CHECK FAILED at %s:%d: ", __FILE__, __LINE__);    \
        fprintf(stderr, __VA_ARGS__);                                          \
        fprintf(stderr, "\n");                                                 \
        abort();                                                               \
    } while (0)

/* Boolean invariant. cond must hold; otherwise abort with the source text. */
#define AIC_CHECK(cond)                                                        \
    do {                                                                       \
        aic_test_checks++;                                                     \
        if (!(cond)) AIC_FAIL("expected (%s) to be true", #cond);              \
    } while (0)

/* Same, with a caller-supplied message (printf-style) for context. */
#define AIC_CHECK_MSG(cond, ...)                                               \
    do {                                                                       \
        aic_test_checks++;                                                     \
        if (!(cond)) AIC_FAIL(__VA_ARGS__);                                    \
    } while (0)

/* Assert two acb balls agree within tol (the temporary ladder's rung 2/3
 * surrogate), aborting with the index-free message at file:line on failure. */
#define AIC_CHECK_ACB_CLOSE(a, b, tol)                                         \
    do {                                                                       \
        aic_test_checks++;                                                     \
        if (!aic_acb_close((a), (b), (tol)))                                   \
            AIC_FAIL("acb balls not within tol");                             \
    } while (0)

/* Same for whole matrices; reports the first disagreeing (row,col). */
#define AIC_CHECK_ACB_MAT_CLOSE(A, B, tol)                                     \
    do {                                                                       \
        aic_test_checks++;                                                     \
        slong _ai = -1, _aj = -1;                                              \
        if (!aic_acb_mat_close((A), (B), (tol), &_ai, &_aj))                   \
            AIC_FAIL("acb_mat entries not within tol (first at %ld,%ld)",      \
                     (long) _ai, (long) _aj);                                  \
    } while (0)

/* Print the check count; call at end of main() (greppable, machine-readable). */
static inline void aic_test_report(const char *name)
{
    printf("CHECKS %s n=%ld\n", name, aic_test_checks);
}

/* Return 1 iff the ball |a-b| is certainly <= tol (see header docstring).
 * prec is the working precision for the subtraction/abs; pass a value at least
 * as large as the precision the balls were computed at. */
static inline int aic_acb_close(const acb_t a, const acb_t b, const arb_t tol)
{
    const slong prec = 64; /* tol comparison only; inputs carry their own radii */
    acb_t diff;
    arb_t mag;
    int ok;

    acb_init(diff);
    arb_init(mag);

    acb_sub(diff, a, b, prec);   /* exact-radius difference of the two balls   */
    acb_abs(mag, diff, prec);    /* arb ball enclosing every |a-b|             */
    ok = arb_le(mag, tol);       /* nonzero iff that whole ball is <= tol       */

    arb_clear(mag);
    acb_clear(diff);
    return ok;
}

/* Matrix version: dimensions must match (fail-loud on mismatch), then every
 * entry must satisfy aic_acb_close. On the first failure returns 0 and writes
 * the (row,col) to *fi,*fj (if non-NULL); returns 1 if all entries agree. */
static inline int aic_acb_mat_close(const acb_mat_t A, const acb_mat_t B,
                                    const arb_t tol, slong *fi, slong *fj)
{
    slong r, c, i, j;

    if (acb_mat_nrows(A) != acb_mat_nrows(B) ||
        acb_mat_ncols(A) != acb_mat_ncols(B)) {
        fprintf(stderr,
                "aic_acb_mat_close: shape mismatch %ldx%ld vs %ldx%ld\n",
                (long) acb_mat_nrows(A), (long) acb_mat_ncols(A),
                (long) acb_mat_nrows(B), (long) acb_mat_ncols(B));
        abort();
    }

    r = acb_mat_nrows(A);
    c = acb_mat_ncols(A);
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            if (!aic_acb_close(acb_mat_entry(A, i, j),
                               acb_mat_entry(B, i, j), tol)) {
                if (fi) *fi = i;
                if (fj) *fj = j;
                return 0;
            }
        }
    }
    return 1;
}

#endif /* AIC_TEST_H */
