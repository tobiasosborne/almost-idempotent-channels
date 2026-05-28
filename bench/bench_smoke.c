/* bench_smoke.c — one real example benchmark (beads aic-73v, "T-harness").
 *
 * Proves the benchmark harness works end to end by timing the smoke op: the
 * acb_mat_mul of a small (4x4) complex matrix at prec=53 (the same library
 * operation test_smoke.c exercises for correctness). This is not a math result
 * — it carries no paper provenance — it exists so `make bench` produces an
 * actual measurement and the AIC_BENCH line format is real, not hypothetical.
 *
 * Real later benchmarks (funcalc sign-route timings, projection-route timings,
 * ...) will follow this exact shape: set up acb inputs once, warm up, then call
 * AIC_BENCH on the hot loop. Output is one greppable BENCH line (CLAUDE.md
 * Rule 12: concrete numbers, no adjectives).
 */
#include <flint/acb_mat.h>

#include "aic_bench.h"

int main(void)
{
    const slong dim = 4;
    const slong prec = 53;
    const long n = 200000;

    acb_mat_t a, b, c;
    acb_mat_init(a, dim, dim);
    acb_mat_init(b, dim, dim);
    acb_mat_init(c, dim, dim);

    /* Deterministic, non-trivial fill: a = (i+1)+(j+1)i, b = identity-ish, so
     * the product is real work (not all-zero, not the I*I exact short path). */
    for (slong i = 0; i < dim; i++) {
        for (slong j = 0; j < dim; j++) {
            acb_set_si_si(acb_mat_entry(a, i, j), (slong) (i + 1), (slong) (j + 1));
        }
    }
    acb_mat_one(b);
    for (slong i = 0; i < dim; i++)
        acb_set_si(acb_mat_entry(b, i, i), 2);

    /* Warm up (page-in, branch predictor) outside the measured loop. */
    acb_mat_mul(c, a, b, prec);

    AIC_BENCH("acb_mat_mul_4x4_p53", n, acb_mat_mul(c, a, b, prec));

    acb_mat_clear(c);
    acb_mat_clear(b);
    acb_mat_clear(a);
    return 0;
}
