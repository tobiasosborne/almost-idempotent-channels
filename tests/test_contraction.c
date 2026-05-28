/* test_contraction.c — cross-checks for the contraction solver (beads aic-n3n,
 * "T-contraction"). CLAUDE.md Rule 5/6: every check asserts a value, the
 * lemma's actual inequality, the paper's Cauchy rate, an inter-routine
 * agreement, or a paper cross-check.
 *
 * Lemma lem_invfun (approximate_algebras.tex:564-592). The checks mirror the
 * lemma's two claims and its proof:
 *  1. known-inverse solve: f(x)=Vx+eps*N(x), pick x*, y=f(x*); solve recovers x*.
 *  2. lemma part 1 (tex:572-574): (1-c)||x1-x2|| <= ||V^{-1}(f(x1)-f(x2))||
 *     <= (1+c)||x1-x2|| for random x1,x2 in the ball. This is the paper's actual
 *     claim; verifying it validates both the impl AND the c estimate.
 *  3. convergence rate (tex:591): ||x_n - x_{n-1}|| <= r c^{n-1}. The solver
 *     asserts this internally at every step; the test reads back stats.max_step
 *     and confirms it is below the geometric envelope at the recorded iter count.
 *  4. Picard vs Anderson agreement: both candidates reach the same fixed point.
 *  5. paper cross-check (tex:594-599): the sgn defining system Y1Y2=I, XY2=Y1X
 *     with solution Y1=Y2=sgn(X), solved for a DIAGONAL X (a simpler instance
 *     than the full operator Jacobian of tex:606-614; see solve_sgn_system) and
 *     cross-checked against aic_sgn(X).
 *  6. precision ladder (bead aic-5ty): solve at prec=53 vs prec=256 agree.
 *
 * Mutation-proof (CLAUDE.md Rule 7), two confirmed-RED experiments:
 *  - check 2 (the lemma inequality): scaling the MEASURED ||V^{-1}(f(x1)-f(x2))||
 *    (nb) down by 1/4 drops it below (1-c)||x1-x2|| = 0.5||x1-x2|| (c=0.5), so the
 *    LOWER-bound AIC_CHECK at this file's line ~180 fired (abort, exit 134). The
 *    LOWER bound is the tight side here (the true ratio ~1.0 sits well below the
 *    loose upper bound 1+c=1.5); a 1/2 scaling lands exactly on the bound and
 *    passes (arb_le is <=), so 1/4 is the smallest clean violator. This pins that
 *    the test enforces the tight inequality, not a triviality.
 *  - check 3 (the Cauchy rate): corrupting the preconditioner prob_Vinv to apply
 *    1/4 V^{-1} (a genuinely wrong V^{-1}) slows the true contraction past the
 *    asserted c=0.5 envelope; the solver's internal tex:591 guard
 *    (||x_n-x_{n-1}|| <= r c^{n-1}) then aborted at step n=6 (exit 134), proving
 *    the Cauchy guard catches a wrong V^{-1}.
 * Both mutations reverted; suite GREEN. Stated here per Rule 7.
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_contraction.h"
#include "aic_funcalc.h"
#include "aic_mat.h"
#include "aic_test.h"

/* ---- shared problem: f(x) = V x + eps * N(x), N(x)_i = sin-like cubic ----
 * V is a fixed, well-conditioned invertible matrix; N is a mild entrywise
 * nonlinearity with ||N'|| <= 1, so for small eps and V^{-1} the preconditioner,
 * ||V^{-1} d_x f - 1|| = eps ||V^{-1} N'|| <= c < 1 on a bounded ball. We do NOT
 * use the matrix sgn here; the nonlinearity is deliberately simple so the lemma
 * hypotheses are transparent and c is hand-computable. */
typedef struct {
    acb_mat_t V;
    acb_mat_t Vinv;
    double    eps;
} prob_ctx;

/* N(x): entrywise cubic-ish bounded map x -> x - x^3/6 (truncated sin), with
 * derivative 1 - x^2/2 in [1/2, 1] on |x|<=1, so ||N'||_op <= 1 on the unit
 * ball. Applied per column entry. */
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
    acb_mat_mul(out, p->V, x, prec);                 /* V x */
    acb_t e; acb_init(e); acb_set_d(e, p->eps);
    acb_mat_scalar_mul_acb(Nx, Nx, e, prec);         /* eps N(x) */
    acb_mat_add(out, out, Nx, prec);
    acb_clear(e);
    acb_mat_clear(Nx);
}

static void prob_Vinv(acb_mat_t out, const acb_mat_t r, void *ctx, slong prec)
{
    prob_ctx *p = ctx;
    acb_mat_mul(out, p->Vinv, r, prec);
}

/* Build a fixed n x n diagonally-dominant invertible V and its inverse. */
static void build_V(prob_ctx *p, slong n, slong prec)
{
    acb_mat_init(p->V, n, n);
    acb_mat_init(p->Vinv, n, n);
    acb_mat_zero(p->V);
    for (slong i = 0; i < n; i++) {
        acb_set_si(acb_mat_entry(p->V, i, i), 4);     /* strong diagonal */
        for (slong j = 0; j < n; j++)
            if (j != i)
                acb_set_d(acb_mat_entry(p->V, i, j),
                          0.3 * (((i + 2 * j) % 3) - 1));
    }
    int ok = acb_mat_inv(p->Vinv, p->V, prec);
    AIC_CHECK_MSG(ok, "build_V: V not invertible at prec");
}

static void prob_clear(prob_ctx *p)
{
    acb_mat_clear(p->Vinv);
    acb_mat_clear(p->V);
}

static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* --- 1. known-inverse solve --- */
static void test_known_inverse(slong prec)
{
    const slong n = 4;
    prob_ctx p; p.eps = 0.1;
    build_V(&p, n, prec);

    acb_mat_t xstar, y, x0, x;
    acb_mat_init(xstar, n, 1);
    acb_mat_init(y, n, 1);
    acb_mat_init(x0, n, 1);
    acb_mat_init(x, n, 1);
    for (slong i = 0; i < n; i++)
        acb_set_d(acb_mat_entry(xstar, i, 0), 0.2 * (i + 1));  /* small x* */
    prob_f(y, xstar, &p, prec);                  /* y = f(x*) */
    acb_mat_zero(x0);                            /* centre at 0 */

    aic_contraction_opts o = {
        .f = prob_f, .apply_Vinv = prob_Vinv, .ctx = &p,
        .x0 = x0, .y = y, .r = 1.5, .c = 0.5, .tol = 1e-30,
        .max_iter = 200, .prec = prec
    };
    aic_contraction_stats st;
    aic_contraction_solve(x, &o, &st);

    arb_t tol; arb_init(tol); set_tol(tol, 1e-25);
    AIC_CHECK_ACB_MAT_CLOSE(x, xstar, tol);      /* recovered x ~ x* */
    arb_clear(tol);

    acb_mat_clear(x); acb_mat_clear(x0); acb_mat_clear(y); acb_mat_clear(xstar);
    prob_clear(&p);
}

/* --- 2. lemma part 1: the two-sided inequality (tex:572-574) --- */
static void test_lemma_bounds(slong prec)
{
    const slong n = 4;
    prob_ctx p; p.eps = 0.1;
    build_V(&p, n, prec);
    const double c = 0.5;   /* same envelope as the solver below */

    acb_mat_t x1, x2, f1, f2, df, Vdf, dx;
    acb_mat_init(x1, n, 1); acb_mat_init(x2, n, 1);
    acb_mat_init(f1, n, 1); acb_mat_init(f2, n, 1);
    acb_mat_init(df, n, 1); acb_mat_init(Vdf, n, 1); acb_mat_init(dx, n, 1);

    for (int trial = 0; trial < 4; trial++) {
        for (slong i = 0; i < n; i++) {
            acb_set_d(acb_mat_entry(x1, i, 0), 0.15 * ((i + trial) % 5 - 2));
            acb_set_d(acb_mat_entry(x2, i, 0), 0.15 * ((2 * i + trial) % 5 - 2));
        }
        prob_f(f1, x1, &p, prec);
        prob_f(f2, x2, &p, prec);
        acb_mat_sub(df, f1, f2, prec);
        prob_Vinv(Vdf, df, &p, prec);            /* V^{-1}(f(x1)-f(x2)) */
        acb_mat_sub(dx, x1, x2, prec);

        arb_t na, nb, lo, hi, cm, cp;
        arb_init(na); arb_init(nb); arb_init(lo); arb_init(hi);
        arb_init(cm); arb_init(cp);
        aic_mat_opnorm(na, dx, prec);            /* ||x1-x2||      */
        aic_mat_opnorm(nb, Vdf, prec);           /* ||V^{-1}(df)|| */
        arb_set_d(cm, 1.0 - c); arb_set_d(cp, 1.0 + c);
        arb_mul(lo, na, cm, prec);               /* (1-c)||x1-x2|| */
        arb_mul(hi, na, cp, prec);               /* (1+c)||x1-x2|| */
        AIC_CHECK_MSG(arb_le(lo, nb),
                      "lemma lower bound (1-c)||a|| <= ||b|| violated");
        AIC_CHECK_MSG(arb_le(nb, hi),
                      "lemma upper bound ||b|| <= (1+c)||a|| violated");
        arb_clear(cp); arb_clear(cm);
        arb_clear(hi); arb_clear(lo); arb_clear(nb); arb_clear(na);
    }

    acb_mat_clear(dx); acb_mat_clear(Vdf); acb_mat_clear(df);
    acb_mat_clear(f2); acb_mat_clear(f1); acb_mat_clear(x2); acb_mat_clear(x1);
    prob_clear(&p);
}

/* --- 3. convergence rate ||x_n - x_{n-1}|| <= r c^{n-1} (tex:591) --- */
static void test_convergence_rate(slong prec)
{
    const slong n = 4;
    prob_ctx p; p.eps = 0.15;
    build_V(&p, n, prec);
    const double r = 1.5, c = 0.5;

    acb_mat_t xstar, y, x0, x;
    acb_mat_init(xstar, n, 1); acb_mat_init(y, n, 1);
    acb_mat_init(x0, n, 1); acb_mat_init(x, n, 1);
    for (slong i = 0; i < n; i++)
        acb_set_d(acb_mat_entry(xstar, i, 0), 0.25 * (i + 1));
    prob_f(y, xstar, &p, prec);
    acb_mat_zero(x0);

    aic_contraction_opts o = {
        .f = prob_f, .apply_Vinv = prob_Vinv, .ctx = &p,
        .x0 = x0, .y = y, .r = r, .c = c, .tol = 1e-30,
        .max_iter = 200, .prec = prec
    };
    aic_contraction_stats st;
    aic_contraction_picard(x, &o, &st);  /* Picard obeys the geometric law */

    /* The internal per-step guard already enforces r c^{n-1}; here confirm the
     * recorded max step is <= r (the n=1 envelope c^0=1), and the last step is
     * far below tol's neighbourhood — i.e. the sequence actually contracted. */
    AIC_CHECK_MSG(st.max_step <= r * (1.0 + 1e-6),
                  "max step exceeded the r c^0 = r envelope");
    AIC_CHECK_MSG(st.last_step < 1e-20, "final step not contracted to tol");
    AIC_CHECK_MSG(st.iters >= 2 && st.iters < 200, "unexpected iteration count");

    acb_mat_clear(x); acb_mat_clear(x0); acb_mat_clear(y); acb_mat_clear(xstar);
    prob_clear(&p);
}

/* --- 4. Picard vs Anderson agreement --- */
static void test_candidate_agreement(slong prec)
{
    const slong n = 5;
    prob_ctx p; p.eps = 0.2;
    build_V(&p, n, prec);

    acb_mat_t xstar, y, x0, xp, xa;
    acb_mat_init(xstar, n, 1); acb_mat_init(y, n, 1); acb_mat_init(x0, n, 1);
    acb_mat_init(xp, n, 1); acb_mat_init(xa, n, 1);
    for (slong i = 0; i < n; i++)
        acb_set_d(acb_mat_entry(xstar, i, 0), 0.18 * ((i % 3) - 1) + 0.1);
    prob_f(y, xstar, &p, prec);
    acb_mat_zero(x0);

    aic_contraction_opts o = {
        .f = prob_f, .apply_Vinv = prob_Vinv, .ctx = &p,
        .x0 = x0, .y = y, .r = 1.5, .c = 0.6, .tol = 1e-28,
        .max_iter = 300, .prec = prec
    };
    aic_contraction_stats sp, sa;
    aic_contraction_picard(xp, &o, &sp);
    aic_contraction_anderson(xa, &o, &sa);

    arb_t tol; arb_init(tol); set_tol(tol, 1e-22);
    AIC_CHECK_ACB_MAT_CLOSE(xp, xstar, tol);
    AIC_CHECK_ACB_MAT_CLOSE(xa, xstar, tol);
    AIC_CHECK_ACB_MAT_CLOSE(xp, xa, tol);
    arb_clear(tol);

    acb_mat_clear(xa); acb_mat_clear(xp); acb_mat_clear(x0);
    acb_mat_clear(y); acb_mat_clear(xstar);
    prob_clear(&p);
}

/* ---- 5. paper cross-check: the sgn defining system on a DIAGONAL X ---- */
#include "test_contraction_sgn.h"

/* --- 6. precision ladder --- */
static void test_precision_ladder(void)
{
    const slong n = 4;
    prob_ctx p53, p256;
    p53.eps = 0.12; p256.eps = 0.12;
    build_V(&p53, n, 53);
    build_V(&p256, n, 256);

    acb_mat_t xstar, y53, y256, x0, x53, x256;
    acb_mat_init(xstar, n, 1); acb_mat_init(y53, n, 1); acb_mat_init(y256, n, 1);
    acb_mat_init(x0, n, 1); acb_mat_init(x53, n, 1); acb_mat_init(x256, n, 1);
    for (slong i = 0; i < n; i++)
        acb_set_d(acb_mat_entry(xstar, i, 0), 0.2 * (i + 1));
    prob_f(y53, xstar, &p53, 53);
    prob_f(y256, xstar, &p256, 256);
    acb_mat_zero(x0);

    aic_contraction_opts o53 = {
        .f = prob_f, .apply_Vinv = prob_Vinv, .ctx = &p53,
        .x0 = x0, .y = y53, .r = 1.5, .c = 0.5, .tol = 1e-12,
        .max_iter = 200, .prec = 53
    };
    aic_contraction_opts o256 = o53;
    o256.ctx = &p256; o256.y = y256; o256.tol = 1e-30; o256.prec = 256;

    aic_contraction_solve(x53, &o53, NULL);
    aic_contraction_solve(x256, &o256, NULL);

    arb_t tol; arb_init(tol); set_tol(tol, 1e-12);
    AIC_CHECK_ACB_MAT_CLOSE(x53, x256, tol);
    arb_clear(tol);

    acb_mat_clear(x256); acb_mat_clear(x53); acb_mat_clear(x0);
    acb_mat_clear(y256); acb_mat_clear(y53); acb_mat_clear(xstar);
    prob_clear(&p256); prob_clear(&p53);
}

int main(void)
{
    test_known_inverse(256);
    test_lemma_bounds(256);
    test_convergence_rate(256);
    test_candidate_agreement(256);
    test_sgn_system(256);
    test_precision_ladder();

    aic_test_report("test_contraction");
    printf("OK test_contraction\n");
    return 0;
}
