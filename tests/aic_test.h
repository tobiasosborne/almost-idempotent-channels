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
 * The cross-check ladder (beads aic-5ty; rung 2 RESTORED by aic-w4o.5).
 *   CLAUDE.md §"cross-check ladder" wants, weakest→strongest:
 *     (1) double-path internal sanity,
 *     (2) double vs arb@prec=53 agreement (~1e-12),
 *     (3) eta=0 exact-idempotent oracle (zero C*-defect),
 *     (4) arb bound certification at the hypothesis boundary.
 *   Rung (2) is NOW ACTIVE: the LAPACK double path (module latd,
 *   include/aic_latd.h) exists, so we can compare a double-path result against
 *   an arb@prec=53 result within a 53-bit tol (~1e-12). The helpers for that are
 *   aic_double_close (a double scalar vs an arb ball) and aic_eigset_close (a
 *   double eigenvalue/singular-value SET vs an arb SET, compared SORTED — since
 *   the two solvers may order/phase differently). The pre-existing acb@53 vs
 *   acb@high-prec ladder (aic_acb_close / aic_acb_mat_close) remains: it is the
 *   precision-self-consistency check, orthogonal to the double-vs-arb check.
 *   The eta=0 oracle (rung 3) still needs the funcalc/ucp/idemp_structure
 *   modules (MODULE_PLAN.md "Suggested first vertical slice") and is deferred to
 *   those beads. This header builds the REUSABLE checkers only; no oracle.
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
#include <flint/arf.h>

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

/* Assert a double-path scalar `x` agrees with an arb ball `ref` within tol (the
 * RESTORED rung-2 check: double(LAPACK) vs arb@prec=53). Forms |x - ref| as an
 * arb ball that encloses the difference (x enters as a zero-radius arb, ref
 * keeps its computed radius), and requires |x-ref| <= tol. Aborts on failure. */
#define AIC_CHECK_DOUBLE_CLOSE(x, ref, tol)                                    \
    do {                                                                       \
        aic_test_checks++;                                                     \
        if (!aic_double_close((x), (ref), (tol)))                              \
            AIC_FAIL("double %.17g not within tol of arb ref", (double) (x));  \
    } while (0)

/* Assert two value SETS agree as sorted sets within tol: `d` is a double[n]
 * (double path), `a` is an arb_ptr of n entries (arb path). Both are read
 * sorted ascending (NEITHER input is mutated — d is sorted in a private copy, a
 * via an index permutation) and compared elementwise. This is the correct
 * comparison for EIGENVALUES / SINGULAR VALUES, where the two solvers may return
 * different orderings (CLAUDE.md Rule 6: compare as sets, not by position).
 * Aborts on the first out-of-tol pair. */
#define AIC_CHECK_EIGSET_CLOSE(d, a, n, tol)                                   \
    do {                                                                       \
        aic_test_checks++;                                                     \
        slong _ek = -1;                                                        \
        if (!aic_eigset_close((d), (a), (n), (tol), &_ek))                     \
            AIC_FAIL("eigenvalue/singular-value set differs at sorted index %ld",\
                     (long) _ek);                                              \
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

/* Return 1 iff the double `x` is within tol of the arb ball `ref`: x enters as
 * a zero-radius arb, |x - ref| keeps ref's radius, and we require it <= tol
 * (arb_le, rigorous: nonzero only when the whole |x-ref| ball is <= tol). */
static inline int aic_double_close(double x, const arb_t ref, const arb_t tol)
{
    const slong prec = 64; /* tol comparison only; ref carries its own radius */
    arb_t xb, diff;
    int ok;

    arb_init(xb);
    arb_init(diff);
    arb_set_d(xb, x);          /* zero-radius arb (x is exact as given)        */
    arb_sub(diff, xb, ref, prec);
    arb_abs(diff, diff);       /* ball enclosing |x - ref|                     */
    ok = arb_le(diff, tol);

    arb_clear(diff);
    arb_clear(xb);
    return ok;
}

/* qsort comparator for ascending doubles (used to sort the double set in place). */
static inline int aic_dbl_cmp(const void *p, const void *q)
{
    double a = *(const double *) p, b = *(const double *) q;
    return (a > b) - (a < b);
}

/* Return 1 iff the double set `d[0..n)` matches the arb set `a[0..n)` within tol,
 * compared SORTED ASCENDING (the correct comparison for eigenvalues / singular
 * values, whose ordering is solver-dependent). NEITHER input is mutated: `d` is
 * copied and the copy sorted (so a caller's U/Sigma/Vt ordering is preserved);
 * `a` is read via a sorted index permutation (by midpoint) so its balls are
 * untouched. On the first out-of-tol pair returns 0 and writes the sorted index
 * to *fk. */
static inline int aic_eigset_close(const double *d, arb_srcptr a, slong n,
                                   const arb_t tol, slong *fk)
{
    const slong prec = 64;
    double *ds = (double *) malloc((size_t) n * sizeof(double));
    if (ds == NULL) { fprintf(stderr, "aic_eigset_close: OOM\n"); abort(); }
    for (slong k = 0; k < n; k++) ds[k] = d[k];
    qsort(ds, (size_t) n, sizeof(double), aic_dbl_cmp);

    /* index permutation sorting `a` ascending by midpoint (insertion sort; n is
     * small here). The arb balls themselves are not modified. */
    slong *perm = (slong *) malloc((size_t) n * sizeof(slong));
    if (perm == NULL) { fprintf(stderr, "aic_eigset_close: OOM\n"); abort(); }
    for (slong k = 0; k < n; k++) perm[k] = k;
    for (slong i = 1; i < n; i++) {
        slong p = perm[i];
        slong j = i - 1;
        while (j >= 0 &&
               arf_cmp(arb_midref(a + perm[j]), arb_midref(a + p)) > 0) {
            perm[j + 1] = perm[j];
            j--;
        }
        perm[j + 1] = p;
    }

    int ok = 1;
    arb_t xb, diff;
    arb_init(xb);
    arb_init(diff);
    for (slong k = 0; k < n; k++) {
        arb_set_d(xb, ds[k]);
        arb_sub(diff, xb, a + perm[k], prec);
        arb_abs(diff, diff);
        if (!arb_le(diff, tol)) {
            if (fk) *fk = k;
            ok = 0;
            break;
        }
    }
    arb_clear(diff);
    arb_clear(xb);
    free(perm);
    free(ds);
    return ok;
}

#endif /* AIC_TEST_H */
