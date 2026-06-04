/* bench_funcalc.c — the sgn(X) head-to-head (beads aic-4bg, "T-funcalc";
 * CLAUDE.md Law 4: performance is a first-class audition gate). These BENCH lines
 * are the evidence for the default aic_sgn dispatch (Newton-Schulz), recorded in
 * the bead aic-4bg report.
 *
 * Candidates timed at prec=53 (the cross-check ladder's low rung) for n = 4, 8,
 * 16 on a representative input: X = 2P - I for P a PERTURBED rank-floor(n/2)
 * projector. This is the project's real shape (a near-projection), in the sgn
 * validity domain ||X^2 - I|| < 1 with a small margin, so both iterations run
 * their full quadratic descent rather than a one-step exact case.
 *   - aic_sgn_newton_schulz : inverse-free, matrix products only.
 *   - aic_sgn_denman_beavers: one certified acb_mat_inv per step (the dominant
 *     cost), plus the same convergence test.
 * Both produce sgn(X) to the same tol; the BENCH ns/op is the discriminator.
 * Concrete numbers only (CLAUDE.md Rule 12).
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_bench.h"
#include "aic/aic_funcalc.h"

/* X = 2P - I for P a rank-r projector perturbed by a small Hermitian term, so
 * ||X^2 - I|| is O(perturbation) < 1 (in-domain, non-degenerate intermediate).
 * A small FIXED eps = 0.02 keeps ||X^2 - I||_op in [0.1, 0.4] across n=4,8,16:
 * the off-diagonal sum that drives ||X^2 - I||_op grows with n (n=16 at eps=0.05
 * already gives ~1.3 > 1, out of domain), and 0.02 leaves a uniform margin at
 * all three sizes so the iteration counts are comparable. (No 1/sqrt(n) scaling
 * because the build links -lflint only, no -lm; eps is a compile-time constant.) */
static void build_perturbed_projector_sign(acb_mat_t X, slong n, slong prec)
{
    slong r = n / 2;
    const double eps = 0.02;
    acb_mat_zero(X);
    /* diagonal 2*P0 - 1 = +1 for i<r, -1 otherwise */
    for (slong i = 0; i < n; i++)
        acb_set_si(acb_mat_entry(X, i, i), (i < r) ? 1 : -1);
    /* small conjugate-symmetric off-diagonal perturbation */
    for (slong i = 0; i < n; i++)
        for (slong j = i + 1; j < n; j++) {
            double re = eps * (((i + 2 * j) % 5) - 2);
            double im = eps * (((3 * i + j) % 5) - 2);
            acb_t z;
            acb_init(z);
            acb_set_d_d(z, re, im);
            acb_set(acb_mat_entry(X, i, j), z);
            acb_conj(acb_mat_entry(X, j, i), z);
            acb_clear(z);
        }
    (void) prec;
}

static void bench_size(slong n, long iters)
{
    const slong prec = 53;
    char name[64];

    acb_mat_t X, S;
    acb_mat_init(X, n, n);
    acb_mat_init(S, n, n);
    build_perturbed_projector_sign(X, n, prec);

    /* Warm up all candidates (page-in, converge once) outside the timers. */
    aic_sgn_newton_schulz(S, X, prec);
    aic_sgn_newton_schulz3(S, X, prec);
    aic_sgn_denman_beavers(S, X, prec);

    snprintf(name, sizeof name, "sgn_newton_schulz_n%ld_p53", (long) n);
    AIC_BENCH(name, iters, aic_sgn_newton_schulz(S, X, prec));

    /* aic-09a cubic NS candidate: cubic conv, 3 matmuls/step. Auditioned, NOT
     * promoted (quadratic NS wins on large-n wall time + tighter balls). */
    snprintf(name, sizeof name, "sgn_newton_schulz3_n%ld_p53", (long) n);
    AIC_BENCH(name, iters, aic_sgn_newton_schulz3(S, X, prec));

    snprintf(name, sizeof name, "sgn_denman_beavers_n%ld_p53", (long) n);
    AIC_BENCH(name, iters, aic_sgn_denman_beavers(S, X, prec));

    acb_mat_clear(S);
    acb_mat_clear(X);
}

int main(void)
{
    bench_size(4, 5000);
    bench_size(8, 2000);
    bench_size(16, 600);
    return 0;
}
