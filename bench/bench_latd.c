/* bench_latd.c — head-to-head double(LAPACK) vs arb@prec=53 for the two heavy
 * ops, Hermitian eigendecomposition and operator norm (beads aic-w4o.5).
 * Performance is a first-class audition gate (CLAUDE.md Law 4); these BENCH
 * lines are the SPEED dimension of the double-vs-arb Pareto comparison. The
 * other two dimensions (certification: arb yes / double no; degeneracy:
 * LAPACK yes / arb-simple-eig no) are properties, not timings, and are recorded
 * in the module report and the test suite (tests/test_latd.c) — not here.
 *
 * Input. A Hermitian matrix with a WELL-SEPARATED spectrum (distinct dominant
 * integer diagonal + small conjugate-symmetric off-diagonals). The separation
 * matters for the arb side: aic_mat_eig_hermitian goes through the certified
 * SIMPLE-eig path, which aborts on (near-)degeneracy; a well-separated spectrum
 * lets it isolate at prec=53 so the timing reflects the converged path, not a
 * precision-starved retry. (The degenerate case — where the double path is the
 * ONLY route — is a correctness matter, exercised in the test, not benched: the
 * arb simple-eig side cannot run it at all, so there is no head-to-head to time.)
 *
 * Both eig benches request eigenVECTORS (the heavier job) for a like-for-like
 * comparison. The double-path array conversion is done ONCE outside the timed
 * loop (it is part of neither solver). Sizes n = 4, 8, 16, 32. Concrete numbers
 * only (CLAUDE.md Rule 12).
 */
#include <complex.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_bench.h"
#include "aic/aic_latd.h"
#include "aic/aic_mat.h"

/* n x n Hermitian, well-separated spectrum: diagonal 10*(i+1), small generic
 * conjugate-symmetric off-diagonals. */
static void build_herm(acb_mat_t H, slong n)
{
    for (slong i = 0; i < n; i++) {
        acb_set_si(acb_mat_entry(H, i, i), 10 * (i + 1));
        for (slong j = i + 1; j < n; j++) {
            slong re = (i + 2 * j) % 5 - 2;
            slong im = (3 * i + j) % 5 - 2;
            acb_set_si_si(acb_mat_entry(H, i, j), re, im);
            acb_set_si_si(acb_mat_entry(H, j, i), re, -im);
        }
    }
}

static void bench_size(slong n, long iters)
{
    const slong prec = 53;
    char name[80];

    acb_mat_t H;
    acb_mat_init(H, n, n);
    build_herm(H, n);

    /* double-path input: convert once (not part of either solver's cost). */
    double _Complex *Ha = malloc((size_t) (n * n) * sizeof(double _Complex));
    double _Complex *V = malloc((size_t) (n * n) * sizeof(double _Complex));
    double *evals = malloc((size_t) n * sizeof(double));
    aic_latd_from_acb_mat(Ha, H);

    /* arb-path scratch. */
    arb_ptr ev_arb = _arb_vec_init(n);
    acb_mat_t V_arb;
    acb_mat_init(V_arb, n, n);
    arb_t op_arb;
    arb_init(op_arb);

    /* Warm up both paths outside the timers (page-in, converge the QR seed). */
    aic_latd_eig_hermitian(evals, V, Ha, n);
    aic_mat_eig_hermitian(ev_arb, V_arb, H, prec);
    (void) aic_latd_opnorm(Ha, n, n);
    aic_mat_opnorm(op_arb, H, prec);

    /* Hermitian eig with eigenvectors: double vs arb@53. */
    snprintf(name, sizeof name, "eig_herm_double_n%ld", (long) n);
    AIC_BENCH(name, iters, aic_latd_eig_hermitian(evals, V, Ha, n));
    snprintf(name, sizeof name, "eig_herm_arb53_n%ld", (long) n);
    AIC_BENCH(name, iters, aic_mat_eig_hermitian(ev_arb, V_arb, H, prec));

    /* operator norm: double vs arb@53. */
    snprintf(name, sizeof name, "opnorm_double_n%ld", (long) n);
    AIC_BENCH(name, iters, (void) aic_latd_opnorm(Ha, n, n));
    snprintf(name, sizeof name, "opnorm_arb53_n%ld", (long) n);
    AIC_BENCH(name, iters, aic_mat_opnorm(op_arb, H, prec));

    arb_clear(op_arb);
    acb_mat_clear(V_arb);
    _arb_vec_clear(ev_arb, n);
    free(evals);
    free(V);
    free(Ha);
    acb_mat_clear(H);
}

int main(void)
{
    /* iteration counts shrink with n so each measured loop stays sub-second;
     * the arb side dominates the wall time, so the counts are sized for it. */
    bench_size(4, 5000);
    bench_size(8, 2000);
    bench_size(16, 500);
    bench_size(32, 150);
    return 0;
}
