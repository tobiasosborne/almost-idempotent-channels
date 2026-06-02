/* aic_adversarial_nla.c — adversarial instance generators, part 2 of 2 (bead
 * aic-dbo.2, Increment 1). Norm / boundary families: gen2 nonnormal_tri, gen6
 * boundary_x2I, gen7 propP_delta. The structural/spectral families (gen1, gen3,
 * gen4, gen5) live in aic_adversarial.c. See aic_adversarial.h for the corpus
 * rationale and the per-generator knob -> property contracts.
 *
 * Determinism: gen2's "random" off-diagonal pattern is a FIXED LCG seeded at
 * AIC_ADV_SEED, so the instance is bit-identical across calls. gen6/gen7 are exact
 * closed forms.
 *
 * Measured anchors (prec=256), asserted by tests/test_adversarial.c:
 *   gen2 n=4 c=0.1: ||[X,X^dag]||_F=0.069, gap(||I-X^2||op - rho)=0.006
 *        n=4 c=10 : ||[X,X^dag]||_F=172.7, gap=86.5    (both grow with c)
 *   gen6 n=3 s=+0.01: ||X^2-I||op=0.990<1 ; s=-0.01: ||X^2-I||op=1.010>=1
 *   gen7 n=3 delta=0.1:  ||P^2-P||op=0.1 ; delta=0.24: ||P^2-P||op=0.24
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_contraction.h"
#include "aic_adversarial.h"

/* One LCG step -> a fixed rational in [-1, 1] with denominator 2^16. Deterministic
 * (single sequenced update of *state). Gives gen2 a reproducible non-normal pattern
 * without an RNG dependency. */
static double adv_lcg_unit(unsigned long *state)
{
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    long q = (long) ((*state >> 32) & 0xFFFFUL); /* 0 .. 65535 */
    return (double) (q - 32768) / 32768.0;       /* [-1, 1) */
}

/* gen2 — Non-normal upper-triangular (nla 2c, deterministic). Diagonal {(i+1)/n}
 * (the controlled eigenvalues; X is triangular so spec(X) = diag), strict upper
 * entries c * fixed_rational[i][j]. Departure from normality ||[X,X^dag]||_F is 0 at
 * c=0 and grows ~linearly in c; rho(I-X^2) = max_i |1-((i+1)/n)^2| = 1-1/n^2 is fixed
 * by the diagonal, so the op-vs-spectral gap ||I-X^2||_op - rho(I-X^2) opens with c.
 * This is exactly the regime where ||I-X^2||_op > 1 (op-basin Newton-Schulz rejects)
 * while rho(I-X^2) < 1 (the global Newton of aic-8hz still converges). */
void aic_adv_nonnormal_tri(acb_mat_t out, slong n, double c, slong prec)
{
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_adv_nonnormal_tri: out must be n x n");
    assert(n >= 2 && "aic_adv_nonnormal_tri: need n >= 2");
    acb_mat_zero(out);
    unsigned long seed = AIC_ADV_SEED;
    arb_t d;
    arb_init(d);
    for (slong i = 0; i < n; i++) {
        /* diagonal eigenvalue (i+1)/n in (0,1] as a certified arb quotient */
        arb_set_si(d, i + 1);
        arb_div_si(d, d, n, prec);
        arb_set(acb_realref(acb_mat_entry(out, i, i)), d);
        for (slong j = i + 1; j < n; j++) {
            double r = adv_lcg_unit(&seed); /* fixed pattern */
            arb_set_d(acb_realref(acb_mat_entry(out, i, j)), c * r);
        }
    }
    arb_clear(d);
}

/* gen6 — Near-boundary ||X^2-I||=1 straddle (nla 7a, lethal #5). X(s) =
 * diag(sqrt(2-s), 1, ..., 1). X^2 - I = diag((2-s)-1, 0, ..., 0) = diag(1-s,0,...),
 * so ||X^2-I||_op = |1-s|: s>0 => 1-s < 1 (certified inside the funcalc domain),
 * s<0 => 1-s > 1 (outside), s=0 => exactly 1 (the ball straddles and the guard must
 * abort). sqrt(2-s) is a certified arb sqrt so the entry — and the whole edge — is a
 * genuine ball, not a double that lies about which side of 1 it is on. Needs s < 2
 * (so 2-s > 0). */
void aic_adv_boundary_x2I(acb_mat_t out, slong n, double s, slong prec)
{
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_adv_boundary_x2I: out must be n x n");
    assert(n >= 1 && "aic_adv_boundary_x2I: need n >= 1");
    assert(s < 2.0 && "aic_adv_boundary_x2I: need s < 2 (so 2-s > 0)");
    acb_mat_one(out); /* identity: all other diagonal entries are 1 */
    arb_t v;
    arb_init(v);
    arb_set_d(v, 2.0 - s);
    arb_sqrt(v, v, prec); /* sqrt(2-s) as a certified ball */
    arb_set(acb_realref(acb_mat_entry(out, 0, 0)), v);
    arb_zero(acb_imagref(acb_mat_entry(out, 0, 0)));
    arb_clear(v);
}

/* gen7 — prop_P delta near 1/4 (nla 7b). P = diag(p, 0, ..., 0) with p in [0,1/2]
 * chosen so ||P^2-P||_op = |p^2 - p| = p - p^2 = delta. Solving p - p^2 = delta:
 *   p = (1 - sqrt(1 - 4 delta)) / 2,   valid for 0 <= delta <= 1/4.
 * The defect is exactly delta, so this rides the prop_P / assoc_ecsa basin edge
 * delta < 1/4 (tex:525): delta < 1/4 => aic_prop_P proceeds; delta = 1/4 (p=1/2) =>
 * the arb_lt(delta, 1/4) guard fails and prop_P aborts. p is a certified arb root so
 * the defect ball is genuine. */
void aic_adv_propP_delta(acb_mat_t out, slong n, double delta, slong prec)
{
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_adv_propP_delta: out must be n x n");
    assert(n >= 1 && "aic_adv_propP_delta: need n >= 1");
    assert(delta >= 0.0 && delta <= 0.25 &&
           "aic_adv_propP_delta: need 0 <= delta <= 1/4");
    acb_mat_zero(out);
    arb_t p, disc;
    arb_init(p);
    arb_init(disc);
    arb_set_d(disc, 1.0 - 4.0 * delta); /* 1 - 4 delta >= 0 on the valid range */
    arb_sqrt(disc, disc, prec);
    arb_one(p);
    arb_sub(p, p, disc, prec);     /* 1 - sqrt(1-4 delta) */
    arb_mul_2exp_si(p, p, -1);     /* /2  => p in [0, 1/2] */
    arb_set(acb_realref(acb_mat_entry(out, 0, 0)), p);
    arb_clear(disc);
    arb_clear(p);
}

/* ---- gen-7c — scalar slow contraction (nla 7c, docs/adversarial/nla.md
 * "7c Contraction Constant c Near 1"; lem_invfun tex:564-592) ---------------
 * The 1x1 scalar contraction with V=I, x0=1, y=0, r=2 and KNOB c (the lemma's
 * contraction constant, tex:569). The MAP IS f(x) = (1-c)*x, NOT f(x) = c*x:
 * the lemma's contraction constant is the iteration rate of g_y, which is
 *   gamma = ||V^{-1} d_x f - 1|| = |(1-c) - 1| = c   (EXACTLY, tex:569),
 * and the Picard iteration itself runs at that rate
 *   g_y(x) = x + V^{-1}(y - f(x)) = x + (0 - (1-c)x) = c*x,   x_n = c^n -> 0.
 * So the iteration RATIO equals the knob c (NOT 1-c). That is the whole point of
 * the §7c family: c near 1 is SLOW (n ~ 1/(1-c) steps), so a tight max_iter cap
 * fires while the math still converges — the "slow contraction vs no contraction"
 * trap. (The docs/adversarial/nla.md recipe literally writes "f(x)=(1-eps_c)*x",
 * i.e. f(x)=(1-c)*x with eps_c=1-c; its own "c=0.9 -> ~230 iters" expected
 * behaviour is ONLY reproduced by this f, not by f(x)=c*x — see paper/FINDINGS.md
 * §C19. f(x)=c*x would give ratio 1-c, i.e. c=0.9 -> ~13 iters, FAST: the wrong,
 * toothless instance.) The two static callbacks below are the solver's f and
 * V^{-1}; they read the caller-owned aic_adv_cedge_ctx so the same compiled fns
 * serve every knob c. */

/* f(x) = (1-c)*x (scalar, 1x1): out[0][0] = (1 - ctx->c) * x[0][0]. The lemma
 * contraction constant is then |d_x f - 1| = |(1-c)-1| = c (the knob). */
static void cedge_f(acb_mat_t out, const acb_mat_t x, void *ctx, slong prec)
{
    const aic_adv_cedge_ctx *p = ctx;
    acb_t a;
    acb_init(a);
    acb_set_d(a, 1.0 - p->c);   /* slope 1-c, so g_y contracts at rate c */
    acb_mul(acb_mat_entry(out, 0, 0), acb_mat_entry(x, 0, 0), a, prec);
    acb_clear(a);
}

/* V^{-1}(r) = r (V = I, the identity preconditioner): out = r. */
static void cedge_Vinv(acb_mat_t out, const acb_mat_t r, void *ctx, slong prec)
{
    (void) ctx;
    (void) prec;
    acb_mat_set(out, r);
}

void aic_adv_contraction_cedge(aic_contraction_opts *o, aic_adv_cedge_ctx *ctx,
                               acb_mat_t x0_out, acb_mat_t y_out,
                               double c, double tol, slong max_iter, slong prec)
{
    assert(acb_mat_nrows(x0_out) == 1 && acb_mat_ncols(x0_out) == 1 &&
           "aic_adv_contraction_cedge: x0_out must be 1x1");
    assert(acb_mat_nrows(y_out) == 1 && acb_mat_ncols(y_out) == 1 &&
           "aic_adv_contraction_cedge: y_out must be 1x1");
    assert(c >= 0.0 && "aic_adv_contraction_cedge: need c >= 0 (a contraction "
           "constant; c>=1 is the fail-loud case the solver guard catches)");
    ctx->c = c;
    acb_mat_one(x0_out);   /* x0 = [[1]] */
    acb_mat_zero(y_out);   /* y  = [[0]] */
    o->f = cedge_f;
    o->apply_Vinv = cedge_Vinv;
    o->ctx = ctx;
    o->x0 = x0_out;
    o->y = y_out;
    o->r = 2.0;            /* tex:586-590 ball B_2(1): |x*-1|=1 < 2 */
    o->c = c;              /* lemma const ||V^{-1} d_x f - 1|| = |(1-c)-1| = c
                            * (tex:569); NOT ||d_x f|| = 1-c. See FINDINGS C19. */
    o->tol = tol;
    o->max_iter = max_iter;
    o->prec = prec;
}
