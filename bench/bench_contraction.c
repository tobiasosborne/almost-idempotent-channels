/* bench_contraction.c — the Picard vs Anderson head-to-head (beads aic-n3n,
 * "T-contraction"; CLAUDE.md Law 4: performance is a first-class audition gate).
 * These BENCH lines are the evidence for the default aic_contraction_solve
 * dispatch, recorded in the bead aic-n3n report.
 *
 * Problem (the test's known-inverse instance, approximate_algebras.tex:564-592):
 *   f(x) = V x + eps N(x),  N(x)_i = x_i - x_i^3/6,  V diagonally dominant.
 * The preconditioner applied is kappa * V^{-1} (kappa in (0,1]), an IMPERFECT
 * preconditioner: g_y(x) = x + kappa V^{-1}(y - f(x)) has derivative
 *   1 - kappa V^{-1} d_x f = (1 - kappa) - kappa eps V^{-1} N',
 * so its norm is c_true ~ (1 - kappa) + O(eps). kappa is the contraction-rate
 * knob; we sweep TWO regimes:
 *   - kappa=1 (FAST, c_true ~ eps small): Picard converges in few steps;
 *     Anderson's per-step inner-product + dg work is pure overhead.
 *   - kappa=0.5 (SLOW, c_true ~ 0.5): Picard crawls at rate 0.5; Anderson's
 *     residual mixing should cut the iteration count.
 * This is the standard "good vs mediocre preconditioner" audition and keeps the
 * problem strictly in-hypothesis: the solution x* is small (||x*|| ~ 0.55) so it
 * sits well inside the guaranteed neighbourhood f(x0)+V(B_{(1-c)r}(0)) (tex:576)
 * for the c we pass. Reported per (n, regime): ns/op for both candidates AND the
 * iteration count each took to the same tol. Concrete numbers (CLAUDE.md Rule 12).
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_bench.h"
#include "aic_contraction.h"

typedef struct { acb_mat_t V, Vinv; double eps; double kappa; } prob_ctx;

static void apply_N(acb_mat_t out, const acb_mat_t x, slong prec)
{
    slong nr = acb_mat_nrows(x);
    acb_t z, z3, six;
    acb_init(z); acb_init(z3); acb_init(six);
    acb_set_si(six, 6);
    for (slong i = 0; i < nr; i++) {
        acb_set(z, acb_mat_entry(x, i, 0));
        acb_pow_ui(z3, z, 3, prec);
        acb_div(z3, z3, six, prec);
        acb_sub(acb_mat_entry(out, i, 0), z, z3, prec);
    }
    acb_clear(six); acb_clear(z3); acb_clear(z);
}

static void prob_f(acb_mat_t out, const acb_mat_t x, void *ctx, slong prec)
{
    prob_ctx *p = ctx;
    acb_mat_t Nx;
    acb_mat_init(Nx, acb_mat_nrows(x), 1);
    apply_N(Nx, x, prec);
    acb_mat_mul(out, p->V, x, prec);
    acb_t e; acb_init(e); acb_set_d(e, p->eps);
    acb_mat_scalar_mul_acb(Nx, Nx, e, prec);
    acb_mat_add(out, out, Nx, prec);
    acb_clear(e);
    acb_mat_clear(Nx);
}

static void prob_Vinv(acb_mat_t out, const acb_mat_t r, void *ctx, slong prec)
{
    prob_ctx *p = ctx;
    acb_mat_mul(out, p->Vinv, r, prec);     /* V^{-1} r */
    acb_t k; acb_init(k); acb_set_d(k, p->kappa);
    acb_mat_scalar_mul_acb(out, out, k, prec); /* kappa V^{-1} r */
    acb_clear(k);
}

static void build_V(prob_ctx *p, slong n, slong prec)
{
    acb_mat_init(p->V, n, n);
    acb_mat_init(p->Vinv, n, n);
    acb_mat_zero(p->V);
    for (slong i = 0; i < n; i++) {
        acb_set_si(acb_mat_entry(p->V, i, i), 4);
        for (slong j = 0; j < n; j++)
            if (j != i)
                acb_set_d(acb_mat_entry(p->V, i, j),
                          0.3 * (((i + 2 * j) % 3) - 1));
    }
    acb_mat_inv(p->Vinv, p->V, prec);
}

/* tag: a short regime label for the BENCH/ITERS lines (kappa has a dot). */
static void bench_case(slong n, double kappa, double c, const char *tag,
                       long iters)
{
    const slong prec = 53;
    char name[80];
    prob_ctx p; p.eps = 0.1; p.kappa = kappa;
    build_V(&p, n, prec);

    acb_mat_t xstar, y, x0, x;
    acb_mat_init(xstar, n, 1); acb_mat_init(y, n, 1);
    acb_mat_init(x0, n, 1); acb_mat_init(x, n, 1);
    for (slong i = 0; i < n; i++)
        acb_set_d(acb_mat_entry(xstar, i, 0), 0.1 * (i + 1)); /* small x* */
    prob_f(y, xstar, &p, prec);
    acb_mat_zero(x0);

    aic_contraction_opts o = {
        .f = prob_f, .apply_Vinv = prob_Vinv, .ctx = &p,
        .x0 = x0, .y = y, .r = 1.5, .c = c, .tol = 1e-12,
        .max_iter = 2000, .prec = prec
    };

    /* Warm up + capture iteration counts (the other half of the audition). */
    aic_contraction_stats sp, sa;
    aic_contraction_picard(x, &o, &sp);
    aic_contraction_anderson(x, &o, &sa);
    printf("ITERS n%ld_%s picard=%ld anderson=%ld\n",
           (long) n, tag, (long) sp.iters, (long) sa.iters);

    snprintf(name, sizeof name, "picard_n%ld_%s_p53", (long) n, tag);
    AIC_BENCH(name, iters, aic_contraction_picard(x, &o, NULL));
    snprintf(name, sizeof name, "anderson_n%ld_%s_p53", (long) n, tag);
    AIC_BENCH(name, iters, aic_contraction_anderson(x, &o, NULL));

    acb_mat_clear(x); acb_mat_clear(x0); acb_mat_clear(y); acb_mat_clear(xstar);
    acb_mat_clear(p.Vinv); acb_mat_clear(p.V);
}

int main(void)
{
    /* FAST regime (kappa=1, c_true ~ eps/4): exact preconditioner; pass c=0.5
     * as a safe envelope. SLOW regime (kappa=0.5, c_true ~ 0.5): mediocre
     * preconditioner; pass c=0.7 (> 0.5 + O(eps), still < 1). */
    bench_case(4, 1.0, 0.5, "fast", 3000);
    bench_case(4, 0.5, 0.7, "slow", 1500);
    bench_case(8, 1.0, 0.5, "fast", 1500);
    bench_case(8, 0.5, 0.7, "slow", 800);
    return 0;
}
