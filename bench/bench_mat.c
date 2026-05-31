/* bench_mat.c — timings for the heavy Layer-0 matrix ops (beads aic-9zh,
 * "T-mat"). Performance is a first-class audition gate (CLAUDE.md Law 4); these
 * BENCH lines are the baseline against which future routes (power iteration for
 * the operator norm, a degenerate-spectrum eig, a full SVD, ...) are measured.
 *
 * Measured at prec=53 (the cross-check ladder's low rung) for n = 4, 8, 16:
 *   - aic_mat_eig_hermitian: certified Hermitian eigendecomposition WITH
 *     eigenvectors (the two-stage approx_eig_qr -> eig_simple path).
 *   - aic_mat_opnorm: operator norm via the Hermitian PSD Gram + global
 *     enclosure of the largest eigenvalue.
 *
 * The benchmark input is a Hermitian matrix with a deliberately well-separated
 * spectrum (a distinct integer diagonal dominating small generic off-diagonals),
 * so acb_mat_eig_simple isolates the eigenvalues at prec=53 and the eig timing
 * reflects the converged path, not a precision-starved retry. Concrete numbers
 * only (CLAUDE.md Rule 12).
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_bench.h"
#include "aic/aic_mat.h"

/* Build an n x n Hermitian matrix with a well-separated spectrum: real diagonal
 * 10*(i+1) (distinct, dominant) plus small conjugate-symmetric off-diagonals. */
static void build_herm(acb_mat_t H, slong n)
{
    for (slong i = 0; i < n; i++) {
        acb_set_si(acb_mat_entry(H, i, i), 10 * (i + 1));
        for (slong j = i + 1; j < n; j++) {
            slong re = (i + 2 * j) % 5 - 2;
            slong im = (3 * i + j) % 5 - 2;
            acb_set_si_si(acb_mat_entry(H, i, j), re, im);
            acb_set_si_si(acb_mat_entry(H, j, i), re, -im); /* conjugate */
        }
    }
}

static void bench_size(slong n, long iters)
{
    const slong prec = 53;
    char name[64];

    acb_mat_t H;
    acb_mat_init(H, n, n);
    build_herm(H, n);

    arb_ptr ev = _arb_vec_init(n);
    acb_mat_t V;
    acb_mat_init(V, n, n);
    arb_t op;
    arb_init(op);

    /* Warm up both ops (page-in, converge the QR seed) outside the timers. */
    aic_mat_eig_hermitian(ev, V, H, prec);
    aic_mat_opnorm(op, H, prec);

    snprintf(name, sizeof name, "eig_hermitian_n%ld_p53", (long) n);
    AIC_BENCH(name, iters, aic_mat_eig_hermitian(ev, V, H, prec));

    snprintf(name, sizeof name, "opnorm_n%ld_p53", (long) n);
    AIC_BENCH(name, iters, aic_mat_opnorm(op, H, prec));

    arb_clear(op);
    acb_mat_clear(V);
    _arb_vec_clear(ev, n);
    acb_mat_clear(H);
}

int main(void)
{
    bench_size(4, 20000);
    bench_size(8, 8000);
    bench_size(16, 2000);
    return 0;
}
