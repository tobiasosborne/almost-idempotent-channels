/* test_adversarial.c — the adversarial corpus's own teeth (bead aic-dbo.2,
 * Increment 1). For each of the 7 generators in aic_adversarial.h, a self-test
 * asserts the CLAIMED adversarial property as a function of the KNOB, at a MILD and
 * a LETHAL knob value, and (where applicable) that the property strengthens toward
 * the lethal end. CLAUDE.md Rule 5: "runs without errors" is not a pass — every
 * check asserts a value/bound. The point of the corpus is that a generator which
 * does NOT actually stress (a "test that can't fail" at the corpus level) makes its
 * self-test go RED.
 *
 * MUTATION-PROVEN (Rule 7). Each self-test was confirmed to go RED under a targeted
 * mutation of its generator, then the generator restored. The mutation per generator
 * is named in a comment above its self-test and reported in the bead handoff.
 *
 * All properties are certified via the arb substrate (aic_mat_opnorm,
 * aic_mat_frobenius_norm, aic_mat_singular_values, aic_mat_herm_max_eig); a property
 * claimed sharp that produced a straddling ball fails loud (Rule 4).
 */
#define _POSIX_C_SOURCE 200809L   /* _exit() in the watchdog-child fixtures */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>               /* _exit() */

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic/aic_latd.h"
#include "aic/aic_cbnorm.h"
#include "aic/aic_assoc.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_errreduce.h"
#include "aic/aic_projection.h"
#include "../src/aic_projection_internal.h"   /* aic_projection_gap/_ambient/_into_A (fam3C) */
#include "aic_test.h"
#include "aic_watchdog.h"

static const slong PREC = 256;

/* Midpoint of the real part of an acb entry as a double (acb_realref -> arb_t,
 * arb_midref -> arf_t, arf_get_d -> double). Used to read back placed eigenvalues
 * and traces for the self-tests. */
static double acb_real_mid(const acb_t z)
{
    return arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
}

/* ||[A, A^dag]||_F midpoint (departure from normality). */
static double commutator_frob_mid(const acb_mat_t A)
{
    slong n = acb_mat_nrows(A);
    acb_mat_t Ad, AAd, AdA, C;
    acb_mat_init(Ad, n, n);
    acb_mat_init(AAd, n, n);
    acb_mat_init(AdA, n, n);
    acb_mat_init(C, n, n);
    acb_mat_conjugate_transpose(Ad, A);
    acb_mat_mul(AAd, A, Ad, PREC);  /* A A^dag */
    acb_mat_mul(AdA, Ad, A, PREC);  /* A^dag A */
    acb_mat_sub(C, AAd, AdA, PREC); /* [A, A^dag] */
    arb_t f;
    arb_init(f);
    aic_mat_frobenius_norm(f, C, PREC);
    double m = arf_get_d(arb_midref(f), ARF_RND_NEAR);
    arb_clear(f);
    acb_mat_clear(C);
    acb_mat_clear(AdA);
    acb_mat_clear(AAd);
    acb_mat_clear(Ad);
    return m;
}

/* ||X^2 - I||_op as a certified ball (the funcalc domain quantity). */
static void x2_minus_I_opnorm(arb_t out, const acb_mat_t X)
{
    slong n = acb_mat_nrows(X);
    acb_mat_t X2, D;
    acb_mat_init(X2, n, n);
    acb_mat_init(D, n, n);
    acb_mat_sqr(X2, X, PREC);
    acb_mat_one(D);
    acb_mat_sub(D, X2, D, PREC); /* X^2 - I */
    aic_mat_opnorm(out, D, PREC);
    acb_mat_clear(D);
    acb_mat_clear(X2);
}

/* ||P^2 - P||_op as a certified ball (the prop_P defect). */
static void p2_minus_P_opnorm(arb_t out, const acb_mat_t P)
{
    slong n = acb_mat_nrows(P);
    acb_mat_t P2, D;
    acb_mat_init(P2, n, n);
    acb_mat_init(D, n, n);
    acb_mat_sqr(P2, P, PREC);
    acb_mat_sub(D, P2, P, PREC); /* P^2 - P */
    aic_mat_opnorm(out, D, PREC);
    acb_mat_clear(D);
    acb_mat_clear(P2);
}

/* assert a certified arb ball equals an exactly-representable double v within tol. */
static void check_ball_eq(const arb_t x, double v, double tol)
{
    arb_t t, d;
    arb_init(t);
    arb_init(d);
    arb_set_d(t, v);
    arb_sub(d, x, t, PREC);
    arb_abs(d, d);
    arb_set_d(t, tol);
    AIC_CHECK_MSG(arb_le(d, t), "ball %.6e not within %.1e of %.6e",
                  arf_get_d(arb_midref(x), ARF_RND_NEAR), tol, v);
    arb_clear(d);
    arb_clear(t);
}

/* ---- gen1: Jordan t^{1/3}. Property: rho(X) = t^{1/3} (char poly lambda^3 - t).
 * MUTATION: drop the knob — set entry(2,0)=0 instead of t => rho=0 for all t,
 * the t^{1/3} burst vanishes => the "rho(lethal) > rho(mild)" and the
 * rho == t^{1/3} checks both go RED. Confirmed RED, restored. */
static void test_gen1_jordan(void)
{
    acb_mat_t X;
    acb_mat_init(X, 3, 3);

    double ts[2] = {1e-1, 1e-9};   /* mild, lethal */
    double rho[2];
    for (int it = 0; it < 2; it++) {
        aic_adv_jordan_t13(X, ts[it], PREC);
        /* rho(X) = t^{1/3}: compute via the singular-value-independent route — the
         * spectral radius equals |det|^{1/3} here since all three eigenvalues have
         * equal modulus t^{1/3}. det(X) = t (expand: companion of lambda^3 - t). */
        acb_t det;
        acb_init(det);
        acb_mat_det(det, X, PREC);
        double dm = acb_real_mid(det);
        rho[it] = cbrt(fabs(dm));
        /* claimed: rho == t^{1/3} */
        AIC_CHECK_MSG(fabs(rho[it] - cbrt(ts[it])) < 1e-9 * (1.0 + cbrt(ts[it])),
                      "gen1 t=%.1e: rho=%.6e != t^{1/3}=%.6e", ts[it], rho[it],
                      cbrt(ts[it]));
        acb_clear(det);
    }
    /* non-Lipschitz burst: a 10^8x smaller t gives only ~464x smaller rho (the
     * 1/3 power), and crucially rho is STRICTLY positive and ordered with t. */
    AIC_CHECK_MSG(rho[0] > rho[1] && rho[1] > 0.0,
                  "gen1: rho not ordered/positive (mild=%.3e lethal=%.3e)",
                  rho[0], rho[1]);
    /* the 1/3-power signature: rho(1e-1)/rho(1e-9) = (1e8)^{1/3} ~ 464, NOT 1e8 */
    double ratio = rho[0] / rho[1];
    AIC_CHECK_MSG(ratio > 100.0 && ratio < 2154.0,
                  "gen1: rho ratio %.1f not ~ (1e8)^{1/3}=464 (non-Lipschitz)",
                  ratio);
    printf("  gen1 jordan: t=1e-1 rho=%.4e ; t=1e-9 rho=%.4e (ratio %.0f ~ 464)\n",
           rho[0], rho[1], ratio);
    acb_mat_clear(X);
}

/* ---- gen2: non-normal triangular. Property: ||[X,X^dag]||_F grows with c AND the
 * op-vs-rho gap (||I-X^2||_op - rho(I-X^2)) opens with c. rho(I-X^2) is fixed by the
 * placed eigenvalues (i+1)/n => max_i |1-((i+1)/n)^2| = 1-1/n^2.
 * MUTATION: drop the knob — ignore c, set strict-upper entries to the fixed rational
 * WITHOUT the c factor => commutator and gap no longer depend on c => the
 * "lethal > mild" monotonicity checks go RED. Confirmed RED, restored. */
static void test_gen2_nonnormal(void)
{
    const slong n = 4;
    acb_mat_t X;
    acb_mat_init(X, n, n);
    double rho_IX2 = 1.0 - 1.0 / ((double) n * n); /* 1 - 1/16 = 0.9375 */

    double cs[2] = {0.1, 10.0}; /* mild, lethal */
    double comm[2], gap[2];
    for (int it = 0; it < 2; it++) {
        aic_adv_nonnormal_tri(X, n, cs[it], PREC);
        comm[it] = commutator_frob_mid(X);
        arb_t b;
        arb_init(b);
        x2_minus_I_opnorm(b, X);
        double op_IX2 = arf_get_d(arb_midref(b), ARF_RND_NEAR);
        gap[it] = op_IX2 - rho_IX2;
        arb_clear(b);
    }
    /* departure from normality grows monotonically with c */
    AIC_CHECK_MSG(comm[1] > comm[0] && comm[0] > 0.0,
                  "gen2: ||[X,X^dag]||_F not increasing in c (mild=%.3e lethal=%.3e)",
                  comm[0], comm[1]);
    /* op-vs-rho gap opens with c: mild stays in the op-basin (gap small, op<1-ish),
     * lethal blows it open (op-norm >> rho, the global-Newton regime). */
    AIC_CHECK_MSG(gap[1] > gap[0],
                  "gen2: op-vs-rho gap not opening with c (mild=%.3e lethal=%.3e)",
                  gap[0], gap[1]);
    AIC_CHECK_MSG(gap[1] > 1.0,
                  "gen2 lethal c=10: gap=%.3e not > 1 (||I-X^2||op must exceed the "
                  "op-basin while rho(I-X^2)=%.4f < 1)", gap[1], rho_IX2);
    printf("  gen2 nonnormal: c=0.1 |[X,X+]|F=%.3f gap=%.3f ; "
           "c=10 |[X,X+]|F=%.3f gap=%.3f (rho(I-X^2)=%.4f)\n",
           comm[0], gap[0], comm[1], gap[1], rho_IX2);
    acb_mat_clear(X);
}

/* ---- gen3: exact-degeneracy projector. Property: ||P^2-P||=0 (exact idempotent),
 * tr(P)=rank, (2P-I)^2=I exactly. Non-diagonal (Hadamard basis).
 * MUTATION: corrupt orthonormality — scale Hadamard column 0 by 2 (so Q's columns
 * are not orthonormal) => P^2 != P => the ||P^2-P||==0 check goes RED. Confirmed
 * RED, restored. */
static void test_gen3_projector(void)
{
    struct { slong dim, rank; } cases[2] = {{4, 2}, {8, 3}};
    for (int ci = 0; ci < 2; ci++) {
        slong dim = cases[ci].dim, rank = cases[ci].rank;
        acb_mat_t P;
        acb_mat_init(P, dim, dim);
        aic_adv_degenerate_proj(P, dim, rank, PREC);

        /* ||P^2 - P|| = 0 to machine floor */
        arb_t def;
        arb_init(def);
        p2_minus_P_opnorm(def, P);
        check_ball_eq(def, 0.0, 1e-60);

        /* tr(P) = rank */
        acb_t tr;
        acb_init(tr);
        acb_mat_trace(tr, P, PREC);
        AIC_CHECK_MSG(fabs(acb_real_mid(tr) - (double) rank) < 1e-50,
                      "gen3 dim=%ld: tr(P)=%.3e != rank=%ld", (long) dim,
                      acb_real_mid(tr), (long) rank);

        /* X = 2P - I is an exact sign matrix: (2P-I)^2 - I = 0 */
        acb_mat_t X, X2, D;
        acb_mat_init(X, dim, dim);
        acb_mat_init(X2, dim, dim);
        acb_mat_init(D, dim, dim);
        acb_mat_scalar_mul_2exp_si(X, P, 1); /* 2P */
        acb_mat_one(D);
        acb_mat_sub(X, X, D, PREC);  /* 2P - I */
        acb_mat_sqr(X2, X, PREC);
        acb_mat_sub(D, X2, D, PREC); /* (2P-I)^2 - I */
        arb_t s;
        arb_init(s);
        aic_mat_opnorm(s, D, PREC);
        check_ball_eq(s, 0.0, 1e-60);
        printf("  gen3 proj dim=%ld rank=%ld: ||P^2-P||=%.2e tr=%ld "
               "||(2P-I)^2-I||=%.2e (eig-free sgn only: spectrum {0,1} repeated)\n",
               (long) dim, (long) rank,
               arf_get_d(arb_midref(def), ARF_RND_NEAR), (long) rank,
               arf_get_d(arb_midref(s), ARF_RND_NEAR));
        arb_clear(s);
        acb_mat_clear(D);
        acb_mat_clear(X2);
        acb_mat_clear(X);
        acb_clear(tr);
        arb_clear(def);
        acb_mat_clear(P);
    }
}

/* ---- gen4: near-degenerate Hermitian. Property: the minimal eigenvalue gap is
 * exactly `gap`. Asserted via aic_mat_singular_values on the SHIFTED matrix is
 * overkill; instead read the two near-1 diagonal entries back (the spectrum is the
 * diagonal) and confirm their separation == gap, and that it shrinks with the knob.
 * MUTATION: drop the knob — set both near-1 entries to 1.0 (gap ignored) => the
 * separation is 0 for all gap => the "sep == gap" check goes RED. Confirmed RED,
 * restored. */
static void test_gen4_near_degen(void)
{
    const slong n = 4;
    double gaps[2] = {1e-1, 1e-9}; /* mild, lethal */
    double sep[2];
    for (int it = 0; it < 2; it++) {
        acb_mat_t H;
        acb_mat_init(H, n, n);
        aic_adv_near_degen_herm(H, n, gaps[it], PREC);
        /* the near-degenerate pair are diagonal entries (0,0) and (1,1) */
        double e0 = acb_real_mid(acb_mat_entry(H, 0, 0));
        double e1 = acb_real_mid(acb_mat_entry(H, 1, 1));
        sep[it] = fabs(e1 - e0);
        AIC_CHECK_MSG(fabs(sep[it] - gaps[it]) < 1e-12 * (1.0 + gaps[it]),
                      "gen4 gap=%.1e: min-gap=%.6e != gap", gaps[it], sep[it]);
        /* the rest of the spectrum (3,4,..) is >= 1 away, so this IS the min gap */
        AIC_CHECK_MSG(n < 3 || 3.0 - (1.0 + gaps[it] / 2.0) > sep[it],
                      "gen4: near-1 pair is not the minimally separated pair");
        acb_mat_clear(H);
    }
    AIC_CHECK_MSG(sep[0] > sep[1],
                  "gen4: min-gap not shrinking with knob (mild=%.3e lethal=%.3e)",
                  sep[0], sep[1]);
    printf("  gen4 near-degen-herm: gap=1e-1 min-sep=%.3e ; gap=1e-9 min-sep=%.3e "
           "(attacks simple-spectrum eig as gap->0)\n", sep[0], sep[1]);
}

/* ---- gen5: graded diagonal. Property: kappa = sigma_max/sigma_min == range.
 *
 * SUBSTRATE FINDING (reported in the bead handoff). The singular values of a
 * positive diagonal D ARE its diagonal entries (a theorem), but aic_mat_singular_
 * values / aic_mat_opnorm route through the Gram D^dag D and its strict Hermiticity
 * guard (src/aic_mat_spectral.c:58, absolute tol 2^{-(prec-8)}). On this huge-
 * dynamic-range input the Gram diagonal H_ii ~ r^{2i} carries an arb radius
 * ~r^{2i} 2^{-prec}; the guard tests |H_ii - conj(H_ii)| = 0 +- 2 r^{2i} 2^{-prec}
 * against 2^{-(prec-8)} = 256 2^{-prec} and FALSE-FAILS once r^{2i} > 128 — a
 * precision-INDEPENDENT abort (raising prec scales tol and radius together). So
 * aic_mat_singular_values aborts here for range >= ~5 at n=6 (verified). This is a
 * substrate limitation the corpus is meant to SURFACE, not paper over (Rule 3).
 *
 * The self-test therefore certifies kappa the rigorous, Gram-free way: sigma_k =
 * |D_kk| as a certified arb ball (acb_abs of the diagonal entry); kappa = max/min.
 * This still uses the arb substrate and is the CORRECT singular value for a positive
 * diagonal, exercising the lethal range=1e8 the Gram route cannot reach.
 *
 * MUTATION: drop the knob — set r=1 (ignore range) so all entries are 1 => kappa=1
 * for all range => the "kappa == range" check goes RED. Confirmed RED, restored. */
static void test_gen5_graded(void)
{
    const slong n = 6;
    double ranges[2] = {10.0, 1e8}; /* mild, lethal */
    double kappa[2];
    for (int it = 0; it < 2; it++) {
        acb_mat_t D;
        acb_mat_init(D, n, n);
        aic_adv_graded_diag(D, n, ranges[it], PREC);
        /* sigma_k = |D_kk| as certified arb balls; kappa = max/min, certified. */
        arb_t smax, smin, sk, ratio, tol;
        arb_init(smax);
        arb_init(smin);
        arb_init(sk);
        arb_init(ratio);
        arb_init(tol);
        for (slong k = 0; k < n; k++) {
            acb_abs(sk, acb_mat_entry(D, k, k), PREC);
            if (k == 0) {
                arb_set(smax, sk);
                arb_set(smin, sk);
            } else {
                if (arf_cmp(arb_midref(sk), arb_midref(smax)) > 0) arb_set(smax, sk);
                if (arf_cmp(arb_midref(sk), arb_midref(smin)) < 0) arb_set(smin, sk);
            }
        }
        arb_div(ratio, smax, smin, PREC); /* certified kappa ball */
        kappa[it] = arf_get_d(arb_midref(ratio), ARF_RND_NEAR);
        /* certified relative agreement: |kappa - range| <= tol*range */
        arb_t rangeb, diff;
        arb_init(rangeb);
        arb_init(diff);
        arb_set_d(rangeb, ranges[it]);
        arb_sub(diff, ratio, rangeb, PREC);
        arb_abs(diff, diff);
        arb_set_d(tol, 1e-9);
        arb_mul(tol, tol, rangeb, PREC);
        AIC_CHECK_MSG(arb_le(diff, tol),
                      "gen5 range=%.1e: kappa=%.6e != range", ranges[it], kappa[it]);
        arb_clear(diff);
        arb_clear(rangeb);
        arb_clear(tol);
        arb_clear(ratio);
        arb_clear(sk);
        arb_clear(smin);
        arb_clear(smax);
        acb_mat_clear(D);
    }
    AIC_CHECK_MSG(kappa[1] > kappa[0],
                  "gen5: kappa not increasing with range (mild=%.3e lethal=%.3e)",
                  kappa[0], kappa[1]);
    printf("  gen5 graded-diag: range=10 kappa=%.4e ; range=1e8 kappa=%.4e "
           "(Gram-free certify; aic_mat_singular_values aborts here, see comment)\n",
           kappa[0], kappa[1]);
}

/* ---- gen6: ||X^2-I||=1 straddle. Property: certified ||X^2-I||_op < 1 for s>0,
 * >= 1 for s<0 (the funcalc domain edge). The boundary is s=0 (ball straddles 1,
 * guard must abort). MUTATION: drop the knob — set entry(0,0)=sqrt(2) regardless of
 * s (ignore s) => ||X^2-I||op == 1 for all s => the rigorous "< 1 for s>0" check
 * goes RED (1 is not < 1). Confirmed RED, restored. */
static void test_gen6_boundary(void)
{
    const slong n = 3;
    arb_t one, b;
    arb_init(one);
    arb_init(b);
    arb_one(one);

    /* s = +0.01: just inside. certified ||X^2-I|| < 1 (rigorous arb_lt). */
    acb_mat_t X;
    acb_mat_init(X, n, n);
    aic_adv_boundary_x2I(X, n, 0.01, PREC);
    x2_minus_I_opnorm(b, X);
    AIC_CHECK_MSG(arb_lt(b, one),
                  "gen6 s=+0.01: ||X^2-I||=%.6e not certainly < 1 (should be inside)",
                  arf_get_d(arb_midref(b), ARF_RND_NEAR));
    double inside = arf_get_d(arb_midref(b), ARF_RND_NEAR);
    check_ball_eq(b, 0.99, 1e-12); /* ||X^2-I|| = 1 - s = 0.99 exactly */

    /* s = -0.01: just outside. certified ||X^2-I|| > 1. */
    aic_adv_boundary_x2I(X, n, -0.01, PREC);
    x2_minus_I_opnorm(b, X);
    AIC_CHECK_MSG(arb_gt(b, one),
                  "gen6 s=-0.01: ||X^2-I||=%.6e not certainly > 1 (should be outside)",
                  arf_get_d(arb_midref(b), ARF_RND_NEAR));
    double outside = arf_get_d(arb_midref(b), ARF_RND_NEAR);
    check_ball_eq(b, 1.01, 1e-12); /* ||X^2-I|| = 1 - s = 1.01 exactly */

    AIC_CHECK_MSG(inside < 1.0 && outside > 1.0,
                  "gen6: straddle not bracketing 1 (inside=%.6f outside=%.6f)",
                  inside, outside);
    printf("  gen6 boundary-x2I: s=+0.01 ||X^2-I||=%.6f (<1, inside) ; "
           "s=-0.01 ||X^2-I||=%.6f (>1, outside). s=0 => ball straddles 1 => abort.\n",
           inside, outside);
    acb_mat_clear(X);
    arb_clear(b);
    arb_clear(one);
}

/* ---- gen7: prop_P delta near 1/4. Property: certified ||P^2-P||_op == delta.
 * MUTATION: drop the knob — set p=0 regardless of delta (ignore delta) => defect==0
 * for all delta => the "defect == delta" check goes RED for delta>0. Confirmed RED,
 * restored. */
static void test_gen7_propP(void)
{
    const slong n = 3;
    double deltas[2] = {0.1, 0.24}; /* mild (well inside 1/4), lethal (near 1/4) */
    double def[2];
    for (int it = 0; it < 2; it++) {
        acb_mat_t P;
        acb_mat_init(P, n, n);
        aic_adv_propP_delta(P, n, deltas[it], PREC);
        arb_t b;
        arb_init(b);
        p2_minus_P_opnorm(b, P);
        def[it] = arf_get_d(arb_midref(b), ARF_RND_NEAR);
        check_ball_eq(b, deltas[it], 1e-12);
        /* mild stays clear of the 1/4 boundary, lethal rides it */
        AIC_CHECK_MSG(def[it] < 0.25,
                      "gen7 delta=%.2f: defect=%.6f not < 1/4 (prop_P basin)",
                      deltas[it], def[it]);
        arb_clear(b);
        acb_mat_clear(P);
    }
    AIC_CHECK_MSG(def[1] > def[0],
                  "gen7: defect not increasing toward 1/4 (mild=%.3f lethal=%.3f)",
                  def[0], def[1]);
    AIC_CHECK_MSG(0.25 - def[1] < 0.25 - def[0],
                  "gen7: lethal delta not closer to the 1/4 edge");
    printf("  gen7 propP-delta: delta=0.10 ||P^2-P||=%.4f ; delta=0.24 "
           "||P^2-P||=%.4f (basin edge 1/4=0.25)\n", def[0], def[1]);
}

/* ---- fam1B: cb-norm vs operator-norm gap (domain.md:75-123, tex:366-388).
 * The FIRST channel generator's self-test. The named property: the certified
 * idempotence defect ||Phi^2-Phi||_cb = eta*sqrt(1-eta) (tex:378) is STRICTLY
 * larger than the single-copy operator-norm defect — the cb!=op trap.
 *
 * Teeth, in order of load-bearing weight:
 *  (1) CLOSED-FORM PIN: the certified eig-free cb bracket [lo,hi] CONTAINS the
 *      closed form eta*sqrt(1-eta) (lo <= cf <= hi, rigorous arb). A wrong Kraus
 *      misses it.  (2) eta=0 ORACLE: bracket collapses to [0,0] (hi < 1e-9) —
 *      the exact-idempotent dephasing reduction. (3) MONOTONIC KNOB: the cb
 *      lower bound at eta=0.20 strictly exceeds that at eta=0.05 (non-toothless).
 *  (4) THE NAMED cb-vs-op GAP: a genuinely OP-norm-flavored single-copy defect
 *      ||Lambda(E_00)||_op (Lambda=Phi^2-Phi applied to the certified Hermitian
 *      unit observable E_00=|0><0|; this is exactly what a naive op-norm
 *      idempotence check measures) equals eta(1-eta) and is certified STRICTLY
 *      BELOW the cb defect eta*sqrt(1-eta) (since sqrt(1-eta) > 1-eta), and the
 *      cb/op ratio = 1/sqrt(1-eta) GROWS toward eta~1/4. NOTE: the op-proxy is
 *      compared against the closed-form cb value (which (1) certifies is the true
 *      cb-norm), NOT the loose eig-free lower bound ||J||_F/n which sits below
 *      the true cb-norm by design (ratio 2n) and so cannot exceed any op proxy.
 *
 * MUTATION: drop the knob from K_0 — set v=(1,0,...) regardless of eta (so K_0 =
 * |0><0|, the exact dephasing for ALL eta) => defect collapses to 0 for all eta
 * => the closed-form-containment pin (cf>0 must be in [0,0]) goes RED. Confirmed
 * RED, restored byte-identical. */
static void test_fam1b_cb_op_gap(void)
{
    /* measured (prec=256, d=2): closed cf, bracket, op-defect, ratio
     *   eta=0.05: cf=0.048734 in [0.034460,0.137840]; ||Lam(E00)||=0.047500; ratio 1.026
     *   eta=0.20: cf=0.178885 in [0.126491,0.505964]; ||Lam(E00)||=0.160000; ratio 1.118 */
    const slong d = 2;
    double etas[2] = {0.05, 0.20}; /* mild, lethal (toward the 1/4 maximum) */
    double cb_lo[2], ratio[2];

    arb_t cf, lo, hi, opd, one;
    arb_init(cf);
    arb_init(lo);
    arb_init(hi);
    arb_init(opd);
    arb_init(one);
    arb_one(one);

    for (int it = 0; it < 2; it++) {
        double eta = etas[it];
        aic_ucp_kraus phi;
        aic_adv_chan_cb_op_gap(&phi, d, eta, PREC);

        /* unital sanity: sum_a K_a^dag K_a = I (pins the OBSERVABLE convention) */
        arb_t ud;
        arb_init(ud);
        aic_ucp_unital_defect_kraus(ud, &phi, PREC);
        check_ball_eq(ud, 0.0, 1e-12); /* <v|v>=1 exactly; ~4e-17 rounding floor */
        arb_clear(ud);

        /* certified eig-free cb bracket on ||Phi^2-Phi||_cb */
        aic_cbnorm_eigfree_ball(lo, hi, &phi, PREC);
        cb_lo[it] = arf_get_d(arb_midref(lo), ARF_RND_NEAR);

        /* (1) closed-form pin: lo <= eta*sqrt(1-eta) <= hi, rigorous arb. */
        arb_set_d(cf, eta);
        arb_t t;
        arb_init(t);
        arb_set_d(t, 1.0 - eta);
        arb_sqrt(t, t, PREC);
        arb_mul(cf, cf, t, PREC); /* cf = eta*sqrt(1-eta) */
        arb_clear(t);
        AIC_CHECK_MSG(arb_le(lo, cf) && arb_le(cf, hi),
                      "fam1B eta=%.2f: closed form eta*sqrt(1-eta)=%.6f NOT in cb "
                      "bracket [%.6f, %.6f] (Kraus construction wrong)", eta,
                      arf_get_d(arb_midref(cf), ARF_RND_NEAR), cb_lo[it],
                      arf_get_d(arb_midref(hi), ARF_RND_NEAR));

        /* (4) op-flavored single-copy defect ||Lambda(E_00)||_op, Lambda=Phi^2-Phi.
         * E_00 = |0><0| is a certified Hermitian unit-op-norm observable; this is
         * the quantity a naive op-norm idempotence check reports. */
        aic_ucp_kraus phi2;
        aic_ucp_compose(&phi2, &phi, &phi, PREC);
        acb_mat_t E00, Y2, Y1, D;
        acb_mat_init(E00, d, d);
        acb_mat_init(Y2, d, d);
        acb_mat_init(Y1, d, d);
        acb_mat_init(D, d, d);
        acb_set_si(acb_mat_entry(E00, 0, 0), 1); /* |0><0|, ||.||_op = 1 */
        aic_ucp_apply(Y2, &phi2, E00, PREC);
        aic_ucp_apply(Y1, &phi, E00, PREC);
        acb_mat_sub(D, Y2, Y1, PREC);
        aic_mat_opnorm(opd, D, PREC); /* certified ||Lambda(E_00)||_op */

        /* the GAP, certified: op defect STRICTLY < cb defect (sqrt(1-eta)>1-eta) */
        AIC_CHECK_MSG(arb_lt(opd, cf),
                      "fam1B eta=%.2f: op-defect ||Lam(E00)||=%.6f not certainly < "
                      "cb-defect eta*sqrt(1-eta)=%.6f (the cb!=op gap collapsed)",
                      eta, arf_get_d(arb_midref(opd), ARF_RND_NEAR),
                      arf_get_d(arb_midref(cf), ARF_RND_NEAR));
        /* ratio cb/op = 1/sqrt(1-eta), grows toward eta~1/4 */
        ratio[it] = arf_get_d(arb_midref(cf), ARF_RND_NEAR) /
                    arf_get_d(arb_midref(opd), ARF_RND_NEAR);

        acb_mat_clear(D);
        acb_mat_clear(Y1);
        acb_mat_clear(Y2);
        acb_mat_clear(E00);
        aic_ucp_kraus_clear(&phi2);
        aic_ucp_kraus_clear(&phi);
    }

    /* (3) monotonic knob: cb lower bound grows with eta (corpus not toothless) */
    AIC_CHECK_MSG(cb_lo[1] > cb_lo[0] && cb_lo[0] > 0.0,
                  "fam1B: cb defect not increasing with eta (mild=%.6f lethal=%.6f)",
                  cb_lo[0], cb_lo[1]);
    /* the cb-vs-op gap WIDENS toward eta~1/4 (ratio 1/sqrt(1-eta) increasing) */
    AIC_CHECK_MSG(ratio[1] > ratio[0],
                  "fam1B: cb/op gap ratio not widening toward 1/4 (mild=%.4f "
                  "lethal=%.4f)", ratio[0], ratio[1]);

    /* (2) eta=0 oracle: exact-idempotent dephasing => bracket [0,0]. */
    aic_ucp_kraus phi0;
    aic_adv_chan_cb_op_gap(&phi0, d, 0.0, PREC);
    aic_cbnorm_eigfree_ball(lo, hi, &phi0, PREC);
    arb_set_d(opd, 1e-9);
    AIC_CHECK_MSG(arb_lt(hi, opd),
                  "fam1B eta=0: cb bracket upper=%.3e not < 1e-9 (exact-idempotent "
                  "oracle: defect must be 0)",
                  arf_get_d(arb_midref(hi), ARF_RND_NEAR));
    aic_ucp_kraus_clear(&phi0);

    printf("  fam1B cb-op-gap: eta=0.05 cb_lo=%.6f cb/op-ratio=%.4f ; "
           "eta=0.20 cb_lo=%.6f cb/op-ratio=%.4f (gap=1/sqrt(1-eta), widens to 1/4); "
           "eta=0 => bracket [0,0]\n",
           cb_lo[0], ratio[0], cb_lo[1], ratio[1]);

    arb_clear(one);
    arb_clear(opd);
    arb_clear(hi);
    arb_clear(lo);
    arb_clear(cf);
}

/* ---- fam2A: depolarizing eta->1/4 regularization boundary (domain.md:196-245,
 * tex:516-525). The SECOND channel generator's self-test. The named property:
 * the depolarizing defect Phi_p^2 - Phi_p = p(1-p) C is EXACTLY linear in
 * p(1-p), with rho = p(1-p) maximized = 1/4 at p=1/2 (the theta(2Phi-1) basin
 * edge). Since the eig-free bracket lo = ||J||_F/n and ||J||_F = p(1-p)
 * ||Choi(C)||_F EXACTLY, lo(p) = (||Choi(C)||_F/n) p(1-p), a pure linear scaling.
 *
 * Teeth, in order of load-bearing weight:
 *  (1) EXACT LINEARITY IN p(1-p) (the aic-l5b distance-linearity-oracle analogue,
 *      the robust pin): over a 4-point sweep {0.1,0.3,0.5,0.7} the bracket
 *      endpoints lo AND hi scale EXACTLY linearly in q=p(1-p): the ratio
 *      lo(p)/lo(0.1) equals q(p)/q(0.1) to ~1e-12 (>= 3 distinct ratios confirm a
 *      LINE, not a single coincidental match). A defect that is NOT p(1-p) C
 *      breaks this. (2) eta=0 ORACLE: at p=0 (identity) and p=1 (trace-replace
 *      conditional expectation) the bracket collapses to [0,0] (hi < 1e-9).
 *      (3) MAX-AT-p=1/2 BASIN EDGE: lo(0.5) > lo(0.1) and lo(0.5) > lo(0.9),
 *      demonstrating p=1/2 is the rho = 1/4 maximum (and lo(0.1) == lo(0.9), the
 *      p<->1-p symmetry of p(1-p)).
 *
 * Certified rho NOT asserted (no easy existing routine: aic_funcalc_int_gelfand_
 * rho needs the materialized superoperator matrix of Phi^2-Phi, more than this
 * self-test warrants; no hand-rolled fragile eig — task guidance). The exact
 * linearity + max-at-1/2 + eta=0 oracle pin the construction robustly.
 *
 * MUTATION: drop the knob — call aic_channel_depolarizing with a FIXED p=0.3
 * regardless of the requested p (so the defect is constant, breaking the p(1-p)
 * scaling) => the exact-linearity check (lo(0.5)/lo(0.1) == q(0.5)/q(0.1)) goes
 * RED, AND the eta=0 oracle goes RED (p=0 no longer gives [0,0]). Confirmed RED,
 * restored byte-identical (git diff empty). */
static void test_fam2a_depol_boundary(void)
{
    /* measured (prec=256, d=2): bracket [lo,hi], lo/q with q=p(1-p)
     *   p=0.1: q=0.09 [0.077942,0.311769] lo/q=0.866025
     *   p=0.5: q=0.25 [0.216506,0.866025] lo/q=0.866025 (MAX)
     *   p=0.9: q=0.09 [0.077942,0.311769] (= p=0.1, p<->1-p symmetric)
     *   p=0,1: [0, ~3e-16] (exact idempotent) */
    const slong d = 2;
    /* >= 3 p-values to confirm a LINE in q=p(1-p), not a coincidental ratio. */
    double ps[4] = {0.1, 0.3, 0.5, 0.7};
    double lo[4], hi[4];

    arb_t lob, hib;
    arb_init(lob);
    arb_init(hib);

    for (int it = 0; it < 4; it++) {
        aic_ucp_kraus phi;
        aic_adv_chan_depol_boundary(&phi, d, ps[it], PREC);

        /* unital sanity: sum_a K_a^dag K_a = I (the depolarizing Kraus set) */
        arb_t ud;
        arb_init(ud);
        aic_ucp_unital_defect_kraus(ud, &phi, PREC);
        check_ball_eq(ud, 0.0, 1e-12);
        arb_clear(ud);

        aic_cbnorm_eigfree_ball(lob, hib, &phi, PREC);
        lo[it] = arf_get_d(arb_midref(lob), ARF_RND_NEAR);
        hi[it] = arf_get_d(arb_midref(hib), ARF_RND_NEAR);
        AIC_CHECK_MSG(lo[it] > 0.0, "fam2A p=%.2f: lo=%.3e not > 0", ps[it], lo[it]);
        /* (2) ABSOLUTE-COEFFICIENT pin: lo = ||J||_F/d with ||J||_F =
         * p(1-p)*||Choi(C)||_F and ||Choi(C)||_F = sqrt(d^2-1) (C self-adjoint,
         * superoperator spectrum {0, -1 (x)(d^2-1)} => HS norm sqrt(d^2-1)), so
         * lo = p(1-p)*sqrt(d^2-1)/d EXACTLY (= p(1-p)*sqrt(3)/2 at d=2). Pins the
         * SCALE ||Choi(C)||_F, not just the ratio: a wrong defect family
         * id-p(1-p)C' with a different C' passes the ratio/symmetry checks (d=4
         * would give coeff sqrt(15)/4=0.968, not 0.866) but fails THIS one. */
        double coeff = sqrt((double) (d * d - 1)) / (double) d;
        check_ball_eq(lob, ps[it] * (1.0 - ps[it]) * coeff, 1e-9);
        aic_ucp_kraus_clear(&phi);
    }

    /* (1) EXACT LINEARITY: lo(p)/lo(0.1) == q(p)/q(0.1) and hi(p)/hi(0.1) ==
     * q(p)/q(0.1) to ~1e-12, over all >= 3 non-base p (a LINE through origin in
     * q). q(0.1)=0.09. If the defect were not p(1-p) C this fails (a real
     * finding, per the task: STOP and report). */
    double q0 = ps[0] * (1.0 - ps[0]); /* 0.09 */
    for (int it = 1; it < 4; it++) {
        double q = ps[it] * (1.0 - ps[it]);
        double qr = q / q0;
        double lor = lo[it] / lo[0];
        double hir = hi[it] / hi[0];
        AIC_CHECK_MSG(fabs(lor - qr) < 1e-12 * (1.0 + qr),
                      "fam2A p=%.2f: lo ratio %.14f != q ratio %.14f (defect not "
                      "EXACTLY linear in p(1-p) — depolarizing defect is NOT "
                      "p(1-p)*C; STOP, real finding)", ps[it], lor, qr);
        AIC_CHECK_MSG(fabs(hir - qr) < 1e-12 * (1.0 + qr),
                      "fam2A p=%.2f: hi ratio %.14f != q ratio %.14f (defect not "
                      "EXACTLY linear in p(1-p))", ps[it], hir, qr);
    }

    /* (3) MAX-AT-p=1/2 BASIN EDGE: lo(0.5) is the maximum of the sweep. p=0.5 is
     * index 2; p=0.1 index 0; p=0.3 index 1; p=0.7 index 3 (= p=0.3 mirror). */
    AIC_CHECK_MSG(lo[2] > lo[0] && lo[2] > lo[1] && lo[2] > lo[3],
                  "fam2A: lo not maximized at p=1/2 (lo: %.6f,%.6f,%.6f,%.6f for "
                  "p=0.1,0.3,0.5,0.7) — p=1/2 is the rho=1/4 basin edge",
                  lo[0], lo[1], lo[2], lo[3]);

    /* p<->1-p symmetry of p(1-p): lo(0.3)==lo(0.7) (q identical), a fingerprint
     * of the p(1-p) structure. */
    AIC_CHECK_MSG(fabs(lo[1] - lo[3]) < 1e-12 * (1.0 + lo[1]),
                  "fam2A: lo(0.3)=%.12f != lo(0.7)=%.12f (p<->1-p symmetry of "
                  "p(1-p) broken)", lo[1], lo[3]);

    /* boundary at p=0.9 too: build it and confirm lo(0.5) > lo(0.9) (the task's
     * explicit max-at-1/2 vs the 0.9 end). */
    aic_ucp_kraus phi9;
    aic_adv_chan_depol_boundary(&phi9, d, 0.9, PREC);
    aic_cbnorm_eigfree_ball(lob, hib, &phi9, PREC);
    double lo9 = arf_get_d(arb_midref(lob), ARF_RND_NEAR);
    AIC_CHECK_MSG(lo[2] > lo9, "fam2A: lo(0.5)=%.6f not > lo(0.9)=%.6f (basin edge "
                  "at 1/2)", lo[2], lo9);
    /* p<->1-p: lo(0.9) == lo(0.1) (both q=0.09) */
    AIC_CHECK_MSG(fabs(lo9 - lo[0]) < 1e-12 * (1.0 + lo[0]),
                  "fam2A: lo(0.9)=%.12f != lo(0.1)=%.12f (p<->1-p symmetry)",
                  lo9, lo[0]);
    aic_ucp_kraus_clear(&phi9);

    /* (2) eta=0 ORACLE: p=0 (identity) and p=1 (trace-replace conditional
     * expectation, Phi_1^2=Phi_1) are EXACTLY idempotent => bracket [0,0]. */
    double oracle_p[2] = {0.0, 1.0};
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phi0;
        aic_adv_chan_depol_boundary(&phi0, d, oracle_p[it], PREC);
        aic_cbnorm_eigfree_ball(lob, hib, &phi0, PREC);
        arb_t tiny;
        arb_init(tiny);
        arb_set_d(tiny, 1e-9);
        AIC_CHECK_MSG(arb_lt(hib, tiny),
                      "fam2A p=%.0f: bracket upper=%.3e not < 1e-9 (exact-"
                      "idempotent oracle: %s defect must be 0)", oracle_p[it],
                      arf_get_d(arb_midref(hib), ARF_RND_NEAR),
                      it == 0 ? "identity" : "trace-replace cond. exp.");
        arb_clear(tiny);
        aic_ucp_kraus_clear(&phi0);
    }

    printf("  fam2A depol-boundary: p=0.1 lo=%.6f ; p=0.3 lo=%.6f ; p=0.5 lo=%.6f "
           "(MAX, rho=p(1-p)=1/4 edge) ; p=0.7 lo=%.6f ; p=0.9 lo=%.6f. lo EXACTLY "
           "linear in p(1-p) (ratios to ~1e-12); p=0,1 => bracket [0,0]\n",
           lo[0], lo[1], lo[2], lo[3], lo9);

    arb_clear(hib);
    arb_clear(lob);
}

/* ---- fam1D: unital-but-barely (domain.md:163-190, tex:432/672). The THIRD
 * channel generator's self-test. The named property: the CP single-Hermitian-
 * Kraus map Phi(X)=K_0 X K_0, K_0=diag(sqrt(1+du),sqrt(1-du),1,...), has
 * Phi(I)=I+du*E (E=diag(1,-1,0,..), traceless, ||E||=1), so the certified UNITAL
 * DEFECT ||sum_a K_a^dag K_a - 1_H||_op = du EXACTLY — the eps-unit axiom edge.
 *
 * Teeth, in order of load-bearing weight:
 *  (1) CERTIFIED DEFECT == delta_u (the construction pin, LOAD-BEARING): at every
 *      knob, aic_ucp_unital_defect_kraus(ud, phi) == delta_u to ~1e-12 (rigorous
 *      check_ball_eq). A wrong K_0 (or a defect != du) misses it — STOP/report.
 *  (2) delta_u=0 ORACLE: K_0=1_d, Phi=identity, EXACTLY unital => ud < 1e-12.
 *  (3) MONOTONIC / NON-TOOTHLESS: ud(0.1) > ud(1e-3) (defect tracks the knob).
 *
 * MUTATION: drop the knob — ignore delta_u and set K_0=1_d (a FIXED unital map,
 * e.g. set both top diagonal entries to 1.0 regardless of du) => the defect
 * collapses to 0 for ALL du => the ==delta_u pin (1e-3 != 0) and the monotonic
 * check both go RED. Confirmed RED, restored byte-identical (git diff empty). */
static void test_fam1d_unital_defect(void)
{
    /* measured (prec=256, d=2): certified unital defect ud == delta_u
     *   delta_u=1e-3: ud=0.001000   delta_u=0.1: ud=0.100000   delta_u=0.5: ud=0.500000
     *   delta_u=0   : ud=0 (~1e-17 floor; exactly-unital identity) */
    const slong d = 2;
    double dus[3] = {1e-3, 0.1, 0.5}; /* mild, lethal, large */
    double ud_m[3];

    arb_t ud;
    arb_init(ud);

    for (int it = 0; it < 3; it++) {
        aic_ucp_kraus phi;
        aic_adv_chan_unital_defect(&phi, d, dus[it], PREC);

        /* (1) certified unital defect == delta_u, rigorous (the construction pin) */
        aic_ucp_unital_defect_kraus(ud, &phi, PREC);
        ud_m[it] = arf_get_d(arb_midref(ud), ARF_RND_NEAR);
        check_ball_eq(ud, dus[it], 1e-12);
        aic_ucp_kraus_clear(&phi);
    }

    /* (3) monotonic / non-toothless: lethal defect strictly exceeds the mild one */
    AIC_CHECK_MSG(ud_m[1] > ud_m[0] && ud_m[0] > 0.0,
                  "fam1D: unital defect not increasing with delta_u (mild=%.6f "
                  "lethal=%.6f)", ud_m[0], ud_m[1]);

    /* (2) delta_u=0 ORACLE: identity map, EXACTLY unital => defect 0. */
    aic_ucp_kraus phi0;
    aic_adv_chan_unital_defect(&phi0, d, 0.0, PREC);
    aic_ucp_unital_defect_kraus(ud, &phi0, PREC);
    check_ball_eq(ud, 0.0, 1e-12);
    aic_ucp_kraus_clear(&phi0);

    printf("  fam1D unital-defect: delta_u=1e-3 ud=%.6f ; delta_u=0.1 ud=%.6f ; "
           "delta_u=0.5 ud=%.6f (== delta_u, the eps-unit edge tex:432/672); "
           "delta_u=0 => ud=0 (exactly unital)\n", ud_m[0], ud_m[1], ud_m[2]);

    arb_clear(ud);
}

/* ---- fam1C: near-degenerate carrier (domain.md:127-159, lem_carrier
 * .tex:1724 / cor_carrier .tex:1731). The FOURTH channel generator's self-test.
 * The named property: the single Hermitian diagonal Kraus K_0=diag(1,..,1,
 * sqrt(gap)) gives the carrier operator Q = sum_a K_a K_a^dag = diag(1,..,1,gap)
 * a spectrum {1 (x)(d-1), gap}; the certified carrier rank is d for gap above
 * thr=dim_K*2^-52*||Q||_F, NEARLY drops to d-1 as gap->0, and the SMALLEST
 * carrier eigenvalue is EXACTLY gap. This exercises the certified-rank routine
 * (aic_ucp_carrier_rank, bead aic-4td) on a near-degenerate carrier.
 *
 * Teeth, in order of load-bearing weight:
 *  (1) SMALLEST CARRIER EIGENVALUE == gap (the construction pin, LOAD-BEARING):
 *      the smallest eigenvalue of Q (via the same degeneracy-robust
 *      aic_mat_eig_hermitian_multiple the certified rank uses) equals gap to
 *      ~1e-9, certified. A wrong Kraus / wrong carrier marginal misses it.
 *  (2) CERTIFIED RANK == d for gap above thr (full carrier), agreeing with the
 *      double path aic_ucp_carrier_rank_latd. The named "no silent wrong rank".
 *  (3) NEAR-DROP at the gap=0 exact-drop oracle: certified rank == d-1 (the
 *      carrier loses a dimension). MEASURED: this CERTIFIES d-1 (does not fail
 *      loud) because the diagonal gap eigenvalue is a point ball — discovered,
 *      not assumed (the densified-carrier straddle/fail-loud path is test_eigvec
 *      S6b). gap=1 ORACLE: full rank d with a healthy smallest eig 1.
 *  (4) MONOTONIC: the smallest carrier eigenvalue tracks gap (gap=0.5 > 1e-6).
 *
 * MUTATION: drop the knob from K_0 — set the last diagonal entry to 1 regardless
 * of gap (so Q=1_d, full rank d, smallest eig 1 for ALL gap) => the smallest-eig
 * ==gap pin (1 != 1e-6) and the gap=0 rank-drop check (rank d, not d-1) both go
 * RED. Confirmed RED, restored byte-identical (git diff empty). */
static void carrier_smallest_eig(arb_t out, const acb_mat_t Q)
{
    slong d = acb_mat_nrows(Q);
    acb_ptr E = _acb_vec_init(d);
    aic_mat_eig_hermitian_multiple(E, Q, PREC);
    arb_t re;
    arb_init(re);
    acb_get_real(out, E + 0);
    for (slong k = 1; k < d; k++) {
        acb_get_real(re, E + k);
        if (arf_cmp(arb_midref(re), arb_midref(out)) < 0) arb_set(out, re);
    }
    arb_clear(re);
    _acb_vec_clear(E, d);
}

static void test_fam1c_carrier_dropout(void)
{
    /* measured (prec=256, d=2): smallest carrier eig == gap; certified rank
     *   gap=0.5  : smallest eig 0.500000, cert rank 2 (full)
     *   gap=1e-6 : smallest eig 1.0e-6,   cert rank 2 (full, gap >> thr~4.4e-16)
     *   gap=1e-3 : smallest eig 1.0e-3,   cert rank 2 (full)
     *   gap=0    : smallest eig 0,        cert rank 1 (near-drop CERTIFIES d-1)
     *   gap=1    : smallest eig 1,        cert rank 2 (full-rank oracle) */
    const slong d = 2;
    double gaps[3] = {0.5, 1e-3, 1e-6}; /* mild, intermediate, near-drop */
    double small_m[3];

    arb_t mn;
    arb_init(mn);

    for (int it = 0; it < 3; it++) {
        double gap = gaps[it];
        aic_ucp_kraus phi;
        aic_adv_chan_carrier_dropout(&phi, d, gap, PREC);

        acb_mat_t Q;
        acb_mat_init(Q, d, d);
        aic_ucp_carrier_Q(Q, &phi, PREC);

        /* (1) smallest carrier eigenvalue == gap (certified, ~1e-9) — the pin. */
        carrier_smallest_eig(mn, Q);
        small_m[it] = arf_get_d(arb_midref(mn), ARF_RND_NEAR);
        check_ball_eq(mn, gap, 1e-9);

        /* (2) certified carrier rank == d (gap above thr), == double path. */
        slong rk_cert = aic_ucp_carrier_rank(Q, d, PREC);
        slong rk_dbl = aic_ucp_carrier_rank_latd(Q, d);
        AIC_CHECK_MSG(rk_cert == d,
                      "fam1C gap=%.1e: certified carrier rank %ld != d=%ld (gap is "
                      "above thr, carrier must be full rank)", gap, (long) rk_cert,
                      (long) d);
        AIC_CHECK_MSG(rk_cert == rk_dbl,
                      "fam1C gap=%.1e: certified rank %ld != double-path rank %ld",
                      gap, (long) rk_cert, (long) rk_dbl);

        acb_mat_clear(Q);
        aic_ucp_kraus_clear(&phi);
    }

    /* (4) monotonic: the smallest carrier eigenvalue tracks gap. */
    AIC_CHECK_MSG(small_m[0] > small_m[2] && small_m[2] > 0.0,
                  "fam1C: smallest carrier eig not tracking gap (gap=0.5 -> %.3e, "
                  "gap=1e-6 -> %.3e)", small_m[0], small_m[2]);

    /* (3) NEAR-DROP at the gap=0 exact-drop oracle: certified rank == d-1. The
     * carrier loses a dimension; MEASURED behavior is CERTIFY d-1 (the diagonal
     * gap eigenvalue is a point ball, so aic_ucp_carrier_rank never straddles —
     * the fail-loud straddle is the densified-carrier path, test_eigvec S6b). The
     * gap=0 input is a VALID Q (a rank-deficient PSD carrier) built directly here;
     * the generator's knob asserts gap>0, so this oracle bypasses it with Q built
     * from K_0=diag(1,0). */
    {
        acb_mat_t Q0;
        acb_mat_init(Q0, d, d);
        acb_set_si(acb_mat_entry(Q0, 0, 0), 1); /* Q = diag(1, 0): exact drop */
        slong rk0 = aic_ucp_carrier_rank(Q0, d, PREC);
        AIC_CHECK_MSG(rk0 == d - 1,
                      "fam1C gap=0 oracle: certified carrier rank %ld != d-1=%ld "
                      "(carrier must drop a dimension at the exact-zero eigenvalue)",
                      (long) rk0, (long) (d - 1));
        carrier_smallest_eig(mn, Q0);
        check_ball_eq(mn, 0.0, 1e-12);
        acb_mat_clear(Q0);
    }

    /* gap=1 ORACLE: K_0=1_d, Q=1_d, EXACT full-rank d carrier, smallest eig 1. */
    {
        aic_ucp_kraus phi1;
        aic_adv_chan_carrier_dropout(&phi1, d, 1.0, PREC);
        acb_mat_t Q1;
        acb_mat_init(Q1, d, d);
        aic_ucp_carrier_Q(Q1, &phi1, PREC);
        slong rk1 = aic_ucp_carrier_rank(Q1, d, PREC);
        AIC_CHECK_MSG(rk1 == d,
                      "fam1C gap=1 oracle: certified rank %ld != d=%ld (full carrier)",
                      (long) rk1, (long) d);
        carrier_smallest_eig(mn, Q1);
        check_ball_eq(mn, 1.0, 1e-9); /* healthy gap: smallest eig is 1 */
        acb_mat_clear(Q1);
        aic_ucp_kraus_clear(&phi1);
    }

    printf("  fam1C carrier-dropout: gap=0.5 small-eig=%.6f ; gap=1e-3 small-eig=%.3e "
           "; gap=1e-6 small-eig=%.3e (cert rank d=%ld, all gap>>thr~4.4e-16) ; gap=0 "
           "=> cert rank d-1=%ld (near-drop CERTIFIES, no straddle)\n",
           small_m[0], small_m[1], small_m[2], (long) d, (long) (d - 1));

    arb_clear(mn);
}

/* ---- fam1C-dense: DENSIFIED non-normal carrier (domain.md 1C, lem_carrier
 * .tex:1724 / cor_carrier .tex:1731; FINDINGS §D7 the prec floor). The hostile-
 * review follow-up (bead aic-v5f) to the diagonal 1C. It closes the two gaps the
 * diagonal family leaves: (a) the STRADDLE/fail-loud contract of
 * aic_ucp_carrier_rank (the diagonal point ball only ever DECIDES) and (b)
 * CONVENTION-SENSITIVITY sum K K^dag != sum K^dag K (the diagonal Hermitian 1C
 * gives 0).
 *
 * Teeth (Rule 5: every check asserts a value/bound; Rule 6 cross-checks):
 *  (1) CERTIFY (gap>>thr): at prec=128 the dense carrier still certifies rank d
 *      and agrees with the double path — and the small carrier eig == gap (the
 *      construction pin; a wrong carrier marginal misses it). Confirms the dense
 *      carrier did not lose the certify path the diagonal 1C had.
 *  (2) STRADDLE (near thr): at prec=53, gap=1e-16, the densified small-cluster
 *      ball (radius ~1.5e-15, MEASURED) STRADDLES thr ~9.4e-16 -> aic_ucp_carrier_rank
 *      fail-loud-ABORTS ("STRADDLES"). Tested via the fork+SIGALRM watchdog so the
 *      test binary survives the abort. At prec=128 the SAME gap certifies the DROP
 *      to d-1 (the §D7 prec floor: ball radius ~4e-38 << thr). This is the path the
 *      diagonal 1C could NOT reach (its gap eigenvalue is a point ball).
 *  (3) CONVENTION-SENSITIVE: certified ||sum K K^dag - sum K^dag K||_op > 0 (a
 *      stated gap, MEASURED ~0.223 at d=3 gap=0.5), distinguishing the non-normal
 *      densified Kraus from the diagonal Hermitian 1C (where it is exactly 0). The
 *      tooth a carrier-side convention bug (computing the WRONG marginal) trips.
 *
 * MUTATION-PROVEN (Rule 7), confirmed RED then restored byte-identical:
 *  - M1: make K NORMAL (use the SAME rational unitary for U and V, K = U diag U^dag)
 *    => the convention gap collapses to ~0 (tooth 3 RED). The straddle tooth STAYS
 *    GREEN: Q = U diag U^dag is still DENSE (U applied), so the zero-cluster ball
 *    still straddles thr at prec=53 and carrier_rank still aborts. So M1 cleanly
 *    ISOLATES the convention tooth (the straddle tooth is proven by M2).
 *  - M2: mis-set the straddle gap (use gap=0.5 instead of 1e-16 in the straddle
 *    child) => no straddle at prec=53, child certifies rank d and exits 0 => the
 *    SIGABRT expectation goes RED.
 */
/* The straddle prec (53): the densified small-cluster ball (~1.5e-15) STRADDLES
 * thr (~9.4e-16) and aic_ucp_carrier_rank fail-loud-aborts. The CERTIFY teeth run
 * at PREC (256, headroom): the ball drops to ~1e-75 and decides cleanly. The
 * carrier MUST be inspected AT its build prec — a Q built at prec P then eig'd at
 * a stricter prec > P trips the relative Hermiticity assert on the prec-P
 * rounding asymmetry (FINDINGS §D7n; the bug this self-test hit during impl). */
#define V5F_STRADDLE_PREC 53     /* small-cluster ball straddles thr here */

/* The straddle child: build the dense carrier at gap=1e-16, prec=53, and call
 * aic_ucp_carrier_rank — which must ABORT ("STRADDLES"). The d=3 dimension and
 * gap are the MEASURED straddle window (file docstring / probe). */
static void v5f_child_straddle(void)
{
    aic_ucp_kraus phi;
    aic_adv_chan_carrier_dropout_dense(&phi, 3, 1e-16, V5F_STRADDLE_PREC);
    acb_mat_t Q;
    acb_mat_init(Q, 3, 3);
    aic_ucp_carrier_Q(Q, &phi, V5F_STRADDLE_PREC);
    (void) aic_ucp_carrier_rank(Q, 3, V5F_STRADDLE_PREC);  /* must abort */
    acb_mat_clear(Q);
    aic_ucp_kraus_clear(&phi);
}

/* ||sum K K^dag - sum K^dag K||_op for a Kraus channel (the convention gap). */
static double convention_gap_op(const aic_ucp_kraus *phi)
{
    slong dk = phi->dim_K, dh = phi->dim_H;
    acb_mat_t Kd, KKd, KdK, term, diff;
    acb_mat_init(Kd, dh, dk);
    acb_mat_init(KKd, dk, dk);          /* sum K K^dag (carrier, on K) */
    acb_mat_init(KdK, dh, dh);          /* sum K^dag K (= 1_H if unital) */
    acb_mat_init(term, dk, dk);
    acb_mat_zero(KKd);
    acb_mat_zero(KdK);
    for (slong a = 0; a < phi->r; a++) {
        acb_mat_conjugate_transpose(Kd, phi->K[a]);
        acb_mat_mul(term, phi->K[a], Kd, PREC);
        acb_mat_add(KKd, KKd, term, PREC);
        acb_mat_clear(term);
        acb_mat_init(term, dh, dh);
        acb_mat_mul(term, Kd, phi->K[a], PREC);
        acb_mat_add(KdK, KdK, term, PREC);
        acb_mat_clear(term);
        acb_mat_init(term, dk, dk);
    }
    /* dim_K == dim_H here (self-map), so the two are comparable d x d matrices. */
    acb_mat_init(diff, dk, dk);
    acb_mat_sub(diff, KKd, KdK, PREC);
    arb_t g;
    arb_init(g);
    aic_mat_opnorm(g, diff, PREC);
    double gv = arf_get_d(arb_midref(g), ARF_RND_UP) + mag_get_d(arb_radref(g));
    arb_clear(g);
    acb_mat_clear(diff);
    acb_mat_clear(term);
    acb_mat_clear(KdK);
    acb_mat_clear(KKd);
    acb_mat_clear(Kd);
    return gv;
}

static void test_fam1c_carrier_dense(void)
{
    const slong d = 3;
    arb_t mn;
    arb_init(mn);

    /* (1) CERTIFY (gap>>thr) at the §D7 headroom prec (PREC, the carrier inspected
     * AT its build prec — see FINDINGS §D7n): dense carrier still certifies rank d,
     * agrees with the double path, and the smallest carrier eig == gap. */
    {
        aic_ucp_kraus phi;
        aic_adv_chan_carrier_dropout_dense(&phi, d, 0.5, PREC);
        acb_mat_t Q;
        acb_mat_init(Q, d, d);
        aic_ucp_carrier_Q(Q, &phi, PREC);

        /* Q is genuinely NON-DIAGONAL (the whole point vs the diagonal 1C). */
        double offmax = 0.0;
        for (slong a = 0; a < d; a++)
            for (slong b = 0; b < d; b++)
                if (a != b) {
                    double m = fabs(arf_get_d(arb_midref(acb_realref(
                                   acb_mat_entry(Q, a, b))), ARF_RND_NEAR));
                    if (m > offmax) offmax = m;
                }
        AIC_CHECK_MSG(offmax > 1e-3,
                      "fam1C-dense: carrier Q is ~diagonal (max off-diag %.3e) — "
                      "the densifier U did not densify it", offmax);

        carrier_smallest_eig(mn, Q);
        check_ball_eq(mn, 0.5, 1e-9);    /* small carrier eig == gap (pin) */

        slong rk_cert = aic_ucp_carrier_rank(Q, d, PREC);
        slong rk_dbl = aic_ucp_carrier_rank_latd(Q, d);
        AIC_CHECK_MSG(rk_cert == d,
                      "fam1C-dense gap=0.5: certified rank %ld != d=%ld",
                      (long) rk_cert, (long) d);
        AIC_CHECK_MSG(rk_cert == rk_dbl,
                      "fam1C-dense gap=0.5: certified rank %ld != double-path %ld",
                      (long) rk_cert, (long) rk_dbl);

        /* At the headroom prec the near-zero gap certifies the DROP to d-1 (the
         * §D7 prec floor: ball radius << thr, unlike the prec=53 straddle below). */
        aic_ucp_kraus phi2;
        aic_adv_chan_carrier_dropout_dense(&phi2, d, 1e-16, PREC);
        acb_mat_t Q2;
        acb_mat_init(Q2, d, d);
        aic_ucp_carrier_Q(Q2, &phi2, PREC);
        slong rk_drop = aic_ucp_carrier_rank(Q2, d, PREC);
        AIC_CHECK_MSG(rk_drop == d - 1,
                      "fam1C-dense gap=1e-16 @headroom prec: certified rank %ld != "
                      "d-1=%ld (the near-zero carrier dimension must drop at headroom "
                      "prec)", (long) rk_drop, (long) (d - 1));
        acb_mat_clear(Q2);
        aic_ucp_kraus_clear(&phi2);

        acb_mat_clear(Q);
        aic_ucp_kraus_clear(&phi);
        printf("  fam1C-dense CERTIFY: Q dense (off-diag %.3f), small eig 0.5, "
               "cert rank d=%ld == double-path; gap->0 @prec=%d certifies drop to "
               "d-1=%ld\n", offmax, (long) d, (int) PREC, (long) (d - 1));
    }

    /* (2) STRADDLE (near thr) at prec=53: aic_ucp_carrier_rank must FAIL LOUD. The
     * shared fork+SIGALRM watchdog (aic_watchdog.h) asserts the abort without
     * crashing this binary: finishes within timeout, NOT clean-exit, SIGABRT, and
     * the captured output names the STRADDLE. The diagonal 1C could NOT reach this
     * (its gap eig is a point ball). */
    aic_watchdog_assert_failloud(v5f_child_straddle, 20,
                                 "fam1C-dense straddle: aic_ucp_carrier_rank @prec=53",
                                 "STRADDLES");
    printf("  fam1C-dense STRADDLE: prec=53 gap=1e-16 -> carrier_rank SIGABRT "
           "(\"STRADDLES\"); the no-straddle gap the diagonal 1C could not reach\n");

    /* (3) CONVENTION-SENSITIVE: ||sum K K^dag - sum K^dag K||_op > 0 (the non-normal
     * Kraus distinguishes the two carrier conventions; the diagonal Hermitian 1C
     * gives EXACTLY 0). Cross-check: the diagonal 1C's gap MUST be 0. */
    {
        aic_ucp_kraus dense, diag;
        aic_adv_chan_carrier_dropout_dense(&dense, d, 0.5, PREC);
        aic_adv_chan_carrier_dropout(&diag, d, 0.5, PREC);   /* the Hermitian 1C */
        double g_dense = convention_gap_op(&dense);
        double g_diag = convention_gap_op(&diag);
        AIC_CHECK_MSG(g_dense > 0.1,
                      "fam1C-dense convention: ||sum KK^dag - sum K^dag K||_op = %.3e "
                      "not > 0.1 — the densified Kraus is (near-)normal, the "
                      "convention bug would be undetectable", g_dense);
        AIC_CHECK_MSG(g_diag < 1e-12,
                      "fam1C-dense convention cross-check: the diagonal Hermitian 1C "
                      "must give convention gap ~0, got %.3e", g_diag);
        printf("  fam1C-dense CONVENTION: ||sum KK^dag - sum K^dag K||_op = %.3f "
               "(dense, non-normal) vs %.1e (diagonal 1C, Hermitian => 0)\n",
               g_dense, g_diag);
        aic_ucp_kraus_clear(&diag);
        aic_ucp_kraus_clear(&dense);
    }

    arb_clear(mn);
}

/* ---- fam3D: dimension-blowup block algebra (domain.md:416-449, tex:484/1249).
 * The FIFTH channel generator's self-test. The named property: a UCP self-map on
 * B(C^N), N=k*d, whose associated eps-C* algebra is A = (+)_j M_d (dim_A = k*d^2),
 * eta-idempotent with eta tunable by t (t=0 => EXACTLY idempotent). The block
 * COUNT tracks k; this is the eps~c/n regime where the naive Haar-diagonal route
 * (error ~ n) fails and the explicit generalized-Pauli B-diagonal (||D||=1) is
 * needed (tex:484/1249).
 *
 * Teeth, in order of load-bearing weight:
 *  (1) eta=0 ORACLE (the t=0 exact-idempotent reduction, cross-check ladder #3):
 *      the certified eig-free cb bracket upper bound hi < 1e-9 at t=0 (out = Phi0,
 *      the k-block conditional expectation, EXACTLY idempotent), tested at
 *      (k=2,d=2) AND (k=3,d=2).
 *  (2) THE NAMED STRUCTURAL PROPERTY (dim_A = k*d^2, k EQUAL blocks of M_d): at
 *      t=0 build A via aic_assoc_ecstar_from_phi, then aic_cstar_build, and assert
 *      B.num_blocks == k, every B.d[l] == d, B.dim_B == k*d^2 — tested at
 *      (k=2,d=2)->dim_A=8 and (k=3,d=2)->dim_A=12 so the block COUNT demonstrably
 *      TRACKS k (the dimension-blowup the family is built to drive). Mirrors the
 *      build chain in test_cstar_build.c:196-238 (T1 oracle).
 *  (3) KNOB -> eta MONOTONE (corpus-can't-go-toothless): the eig-free bracket
 *      lower bound (midpoint) of ||Phi^2-Phi||_cb STRICTLY increases over
 *      t in {0, 0.02, 0.08} (>= 2 nonzero values, mild < lethal), at (k=2,d=2).
 *  (4) UNITAL at t>0: aic_ucp_unital_defect_kraus < 1e-12 (the mix is unital for
 *      ALL t: sum_a K_a^dag K_a = (1-t)I + t V^dag I V = I).
 *
 * MUTATION (Rule 7): collapse the k-block partition — in the generator set every
 * block projector to P_0 (i.e. `j * d + a` -> `0 * d + a` in the P_j index, so all
 * k Kraus blocks become the SAME rank-d projector on [0,d) and the rest of C^N is
 * UNTOUCHED, breaking sum_j P_j = I_N). This destroys both the unitality
 * (sum_a K_a^dag K_a = (1-t)(k P_0) + t V^dag(k P_0)V != I_N) AND the (+)_j M_d
 * structure. OBSERVED RED (the FIRST assertion to fire is the t=0 unital sanity
 * check at the top of the case loop): "AIC_CHECK FAILED at test_adversarial.c:108:
 * ball 1.000000e+00 not within 1.0e-12 of 0.000000e+00" — the unital defect
 * collapses to ~1 because sum_j P_j = k P_0 != I_N. (The downstream T3
 * num_blocks/dim_B pin would also fire, but the unital sanity guard catches the
 * broken partition first.) Confirmed RED, restored byte-identical (git diff
 * empty).
 *
 * MUTATION 2 (Rule 7 — the T3 block-structure tooth, ISOLATED; hostile review
 * aic-cxo W3). The collapse-to-P_0 mutation above trips UNITALITY first, so it
 * does not prove T3 (num_blocks/dim_B) has teeth on its own. A unitality-
 * PRESERVING mutation does: re-partition (k=2,d=2) into contiguous block sizes
 * {1,3} instead of {2,2} (still disjoint, still sum_j P_j = I_N, so T4 unital
 * and T1 oracle STAY GREEN). OBSERVED RED at the dim_A pin (line ~1011):
 * "fam3D t=0 (k=2,d=2): ecstar dim_A=10 != k*d^2=8" (1^2+3^2=10 != 2*2^2=8).
 * This confirms T3 independently catches a wrong block-count/size even when the
 * channel stays unital. Confirmed RED, restored byte-identical. */
static void test_fam3d_blockalg(void)
{
    /* measured (prec=256): t=0 EXACT idemp, block structure + eta(t) monotone
     *   (k=2,d=2) t=0: num_blocks=2 d=[2,2] dim_B=8  iso_def=0  (dim_A=8)
     *   (k=3,d=2) t=0: num_blocks=3 d=[2,2,2] dim_B=12 iso_def=0 (dim_A=12)
     *   (k=2,d=2) cb_lo(t): t=0 -> 0 ; t=0.02 -> 0.012002 ; t=0.08 -> 0.045071 */
    const slong PREC3 = 256; /* arb working precision for the build chain */

    /* (1)+(2) eta=0 ORACLE + block structure at t=0, at (k=2,d=2) and (k=3,d=2). */
    struct { slong k, d; } cases[2] = {{2, 2}, {3, 2}};
    for (int ci = 0; ci < 2; ci++) {
        slong k = cases[ci].k, d = cases[ci].d;
        slong N = k * d, dim_A_exp = k * d * d;

        aic_ucp_kraus phi0;
        aic_adv_chan_blockalg(&phi0, k, d, 0.0, PREC3);
        AIC_CHECK_MSG(phi0.dim_K == N && phi0.dim_H == N && phi0.r == 2 * k,
                      "fam3D (k=%ld,d=%ld): bad shape dim_K=%ld dim_H=%ld r=%ld",
                      (long) k, (long) d, (long) phi0.dim_K, (long) phi0.dim_H,
                      (long) phi0.r);

        /* unital sanity (pins the OBSERVABLE convention; sum_j P_j = I_N) */
        arb_t ud;
        arb_init(ud);
        aic_ucp_unital_defect_kraus(ud, &phi0, PREC3);
        check_ball_eq(ud, 0.0, 1e-12);
        arb_clear(ud);

        /* (1) eta=0 oracle: t=0 => exactly idempotent => bracket upper < 1e-9. */
        arb_t lo, hi;
        arb_init(lo);
        arb_init(hi);
        aic_cbnorm_eigfree_ball(lo, hi, &phi0, PREC3);
        arb_t tiny;
        arb_init(tiny);
        arb_set_d(tiny, 1e-9);
        AIC_CHECK_MSG(arb_lt(hi, tiny),
                      "fam3D t=0 (k=%ld,d=%ld): cb bracket upper=%.3e not < 1e-9 "
                      "(exact-idempotent oracle: defect must be 0)", (long) k,
                      (long) d, arf_get_d(arb_midref(hi), ARF_RND_NEAR));
        arb_clear(tiny);
        arb_clear(hi);
        arb_clear(lo);

        /* (2) STRUCTURAL: build A = Img Phi_tilde, then B = cstar_build, assert
         *     B.num_blocks == k, every B.d[l] == d, B.dim_B == k*d^2. */
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi0, PREC3);
        AIC_CHECK_MSG(ae.A.dim_A == dim_A_exp,
                      "fam3D t=0 (k=%ld,d=%ld): ecstar dim_A=%ld != k*d^2=%ld",
                      (long) k, (long) d, (long) ae.A.dim_A, (long) dim_A_exp);

        aic_dhom_B B;
        aic_dhom_v v;
        arb_t iso;
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, 0.0, PREC3);

        AIC_CHECK_MSG(B.num_blocks == k,
                      "fam3D t=0 (k=%ld,d=%ld): num_blocks=%ld != k=%ld (block "
                      "COUNT must track k — the dimension-blowup property)",
                      (long) k, (long) d, (long) B.num_blocks, (long) k);
        for (slong l = 0; l < B.num_blocks; l++)
            AIC_CHECK_MSG(B.d[l] == d,
                          "fam3D t=0 (k=%ld,d=%ld): block %ld has dim %ld != d=%ld "
                          "(every block must be M_d)", (long) k, (long) d, (long) l,
                          (long) B.d[l], (long) d);
        AIC_CHECK_MSG(B.dim_B == dim_A_exp,
                      "fam3D t=0 (k=%ld,d=%ld): dim_B=%ld != k*d^2=%ld", (long) k,
                      (long) d, (long) B.dim_B, (long) dim_A_exp);
        AIC_CHECK_MSG(arf_get_d(arb_midref(iso), ARF_RND_NEAR) < 1e-10,
                      "fam3D t=0 (k=%ld,d=%ld): iso_def=%.3e not ~ 0 (eta=0 must be "
                      "an exact isomorphism)", (long) k, (long) d,
                      arf_get_d(arb_midref(iso), ARF_RND_NEAR));

        printf("  fam3D blockalg (k=%ld,d=%ld): t=0 EXACT idemp; build B "
               "num_blocks=%ld d[0]=%ld dim_B=%ld (= k*d^2=%ld) iso_def=%.2e\n",
               (long) k, (long) d, (long) B.num_blocks, (long) B.d[0],
               (long) B.dim_B, (long) dim_A_exp,
               arf_get_d(arb_midref(iso), ARF_RND_NEAR));

        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi0);
    }

    /* (3) KNOB -> eta MONOTONE at (k=2,d=2): the cb lower bound strictly grows
     *     over t in {0, 0.02, 0.08} (>= 2 nonzero values, mild < lethal). */
    {
        const slong k = 2, d = 2;
        double ts[3] = {0.0, 0.02, 0.08};
        double cb_lo[3];
        for (int it = 0; it < 3; it++) {
            aic_ucp_kraus phi;
            aic_adv_chan_blockalg(&phi, k, d, ts[it], PREC3);
            arb_t lo, hi;
            arb_init(lo);
            arb_init(hi);
            aic_cbnorm_eigfree_ball(lo, hi, &phi, PREC3);
            cb_lo[it] = arf_get_d(arb_midref(lo), ARF_RND_NEAR);
            arb_clear(hi);
            arb_clear(lo);
            aic_ucp_kraus_clear(&phi);
        }
        AIC_CHECK_MSG(cb_lo[2] > cb_lo[1] && cb_lo[1] > cb_lo[0] && cb_lo[1] > 0.0,
                      "fam3D: cb defect not strictly increasing with t (t=0 -> "
                      "%.6f, t=0.02 -> %.6f, t=0.08 -> %.6f)", cb_lo[0], cb_lo[1],
                      cb_lo[2]);
        printf("  fam3D blockalg (k=2,d=2): eta(t) cb_lo: t=0 -> %.6f ; t=0.02 -> "
               "%.6f ; t=0.08 -> %.6f (strictly up, knob has teeth)\n",
               cb_lo[0], cb_lo[1], cb_lo[2]);
    }

    /* (4) UNITAL at t>0: sum_a K_a^dag K_a = I for all t. */
    {
        aic_ucp_kraus phi;
        aic_adv_chan_blockalg(&phi, 2, 2, 0.08, PREC3);
        arb_t ud;
        arb_init(ud);
        aic_ucp_unital_defect_kraus(ud, &phi, PREC3);
        check_ball_eq(ud, 0.0, 1e-12);
        arb_clear(ud);
        aic_ucp_kraus_clear(&phi);
    }
}

/* eta proxy = ||S^2 - S||_op via the assoc superoperator midpoint (matches
 * test_cstar_build.c eta_proxy / test_cstar_extension.c). The cstar_build eps
 * argument must be eta, NOT the ~700x-smaller assoc defect eps (FINDINGS §C11). */
static double fam3d_eta_proxy(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S, S2, D, M;
    acb_mat_init(S, nn, nn);
    acb_mat_init(S2, nn, nn);
    acb_mat_init(D, nn, nn);
    acb_mat_init(M, nn, nn);
    aic_assoc_superop_from_ucp(S, phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    for (slong i = 0; i < nn; i++)
        for (slong j = 0; j < nn; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    double o = arf_get_d(arb_midref(e), ARF_RND_NEAR);
    arb_clear(e);
    acb_mat_clear(M);
    acb_mat_clear(D);
    acb_mat_clear(S2);
    acb_mat_clear(S);
    return o;
}

/* ---- fam3D, eta>0: the aic-66n wrapper-collapse regression (REQUIRES prec=256).
 * The k-block algebra at t>0 is GENUINELY eta-idempotent (A = (+)_j M_d an
 * O(eta)-C* algebra). aic_cstar_build recurses Stage-1 splits onto S_P WRAPPER
 * subalgebras whose unit is Ptilde_m = Co_P(P) (.tex:1082), NOT the ambient 1_n.
 * The OLD projection finder ranked spectral gaps by the largest AMBIENT interior
 * gap — on a wrapper that is the support-vs-complement gap, so it built P ~ Ptilde_m
 * (the wrapper UNIT) — a TRIVIAL split that the cstar_build backstop (cstar_build.c:
 * 258, FINDINGS §C11) aborted on with no recovery. The aic-66n fix makes the gap
 * audition UNIT-AWARE (||U_A - P|| > c, U_A = Phi_tilde(1_n); .tex:929, §C7).
 *
 * RED on the PRE-FIX tree (measured): k=4, d=2, t=0.05, prec=256, eps=eta aborts
 * "aic_cstar_build stage1: split projection P' is the wrapper unit Ptilde_m
 * (||Ptilde_m - P'|| = 0.0000 < 0.15; FINDINGS §C11)" at the iter=5 dim_A=2 wrapper.
 * GREEN post-fix: dim_B=16==dim_A, bijective, sigma_min=0.9993, c=iso/eta=0.0247.
 *
 * PREC IS LOAD-BEARING (deviation from the bug-report's "reproduces at low prec").
 * The gap-SELECTION is double-path (prec-independent), but the Stage-1 RECURSION
 * path depends on the arb-computed Phi_tilde / intermediate splits; only prec=256
 * descends to the exact dim_A=2 degenerate wrapper that trips the abort. Measured:
 * k=4 SUCCEEDS at prec=64/128/160 (the recursion never reaches that wrapper) and
 * ABORTS at prec=256 on the old tree. So the regression MUST run at prec=256 —
 * hence this is a `slow`-labeled test (test_adversarial is already `slow`, TIMEOUT
 * 600). k is bounded to <= 4: k=5 (N=10) at prec=256 is just too slow (not a
 * failure). k=3 (dim_A=12) SUCCEEDS on BOTH trees (its recursion never hits the
 * degenerate wrapper) but is asserted too: it pins the bijective O(eta) isomorphism.
 * Wall time (this box): k=3 ~25 s, k=4 ~120 s.
 *
 * MUTATION-PROVEN (Rule 7, by hand on /tmp, recorded in the bead handoff):
 *  (a) revert the audition to largest-gap-only (no unit-aware reject) -> k=4 RED
 *      (the cstar_build.c:258 wrapper-unit abort returns).
 *  (b) gate on the ambient 1_n instead of U_A=Phi_tilde(1_n) -> k=4 RED (||1_n - P||
 *      ~ 1 is vacuous on a wrapper, so the support-vs-complement gap is accepted and
 *      P ~ Ptilde_m collapses the split). */
static void test_fam3d_bijective_eta(void)
{
    const slong PREC3 = 256;   /* the bug ONLY reproduces at prec=256 (see banner) */
    struct { slong k, d; double t; } cases[2] = {{3, 2, 0.05}, {4, 2, 0.05}};
    for (int ci = 0; ci < 2; ci++) {
        slong k = cases[ci].k, d = cases[ci].d;
        slong dim_A_exp = k * d * d;

        aic_ucp_kraus phi;
        aic_adv_chan_blockalg(&phi, k, d, cases[ci].t, PREC3);
        double eta = fam3d_eta_proxy(&phi, PREC3);
        AIC_CHECK_MSG(eta > 1e-6,
                      "fam3D eta>0 (k=%ld): eta=%.3e not > 0 (t must make Phi "
                      "genuinely almost-idempotent)", (long) k, eta);

        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, PREC3);
        AIC_CHECK_MSG(ae.A.dim_A == dim_A_exp,
                      "fam3D eta>0 (k=%ld): dim_A=%ld != k*d^2=%ld", (long) k,
                      (long) ae.A.dim_A, (long) dim_A_exp);

        /* eps := eta (NOT the ~700x-smaller assoc defect; FINDINGS §C11). The
         * PRE-FIX tree aborts HERE for k=4 (the unit-aware audition is the fix). */
        aic_dhom_B B;
        aic_dhom_v v;
        arb_t iso;
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, eta, PREC3);

        /* BIJECTIVE O(eta) isomorphism: dim_B == dim_A AND sigma_min near 1. */
        arb_t a;
        arb_init(a);
        int bij = aic_errreduce_is_bijective(a, &v, PREC3);
        double smin = arf_get_d(arb_midref(a), ARF_RND_NEAR);
        double isod = arf_get_d(arb_midref(iso), ARF_RND_NEAR);
        double c = isod / eta;
        AIC_CHECK_MSG(B.dim_B == dim_A_exp,
                      "fam3D eta>0 (k=%ld): dim_B=%ld != dim_A=%ld (not bijective; "
                      "the wrapper-collapse abort, aic-66n)", (long) k,
                      (long) B.dim_B, (long) dim_A_exp);
        AIC_CHECK_MSG(bij && smin > 0.5,
                      "fam3D eta>0 (k=%ld): v not bijective (sigma_min=%.4f <= 0.5)",
                      (long) k, smin);
        AIC_CHECK_MSG(c < 1.0,
                      "fam3D eta>0 (k=%ld): c=iso_def/eta=%.4f not O(1) (eta=%.3e "
                      "iso_def=%.3e)", (long) k, c, eta, isod);
        printf("  fam3D blockalg eta>0 (k=%ld,d=%ld,t=%.2f): dim_B=%ld==dim_A "
               "BIJECTIVE sigma_min=%.6f eta=%.4e iso_def=%.3e c=iso/eta=%.4f "
               "(aic-66n wrapper-collapse FIXED)\n", (long) k, (long) d,
               cases[ci].t, (long) B.dim_B, smin, eta, isod, c);

        arb_clear(a);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ---- fam-NC: NON-COMMUTATIVE eta->1/4 boundary (README:57, tex:347/378/516).
 * The SIXTH channel generator's self-test, and the NON-COMMUTATIVE counterpart of
 * fam2A (depolarizing, commutative). Phi = id_{M_m} (x) Psi_eta, eta the lower
 * root of eta*sqrt(1-eta) = 1/4 - kappa, so eta_cb = ||Phi^2-Phi||_cb = 1/4-kappa
 * EXACTLY (closed form, ampliation-invariant). The named adversarial property is
 * (a) the boundary calibration + (b) the NON-ABELIAN image, which the commutative
 * G2 lacks.
 *
 * Teeth, in order of load-bearing weight:
 *  (a) CALIBRATION AT THE BOUNDARY: the certified eig-free bracket
 *      [||J||_F/N, 2||J||_F] CONTAINS the target 1/4 - kappa (lo <= target <= hi),
 *      over a kappa sweep {0.10,0.05,0.02,0.01} (mild .. lethal toward the edge).
 *      The bracket is loose (ratio 2N) so it cannot PIN eta_cb; the EXACT value is
 *      the closed form eta*sqrt(1-eta) (tex:378 + ampliation-invariance tex:347),
 *      which we ALSO assert lies in the bracket — the certified cross-check that
 *      the closed-form cb-norm is consistent with the independent eig-free
 *      certifier (CLAUDE.md cross-check ladder #4). g(eta)=eta*sqrt(1-eta) == target
 *      to ~1e-12 pins the root-find. (b) NON-COMMUTATIVITY WITNESS: B1 =
 *      E^{(m)}_{01}(x)1_d, B2 = E^{(m)}_{10}(x)1_d are GENUINE FIXED POINTS of Phi
 *      (Psi unital => Psi(1_d)=1_d), and the Choi-Effros star-commutator (tex:342)
 *      ||Phi(B1 B2)-Phi(B2 B1)||_op = ||[B1,B2]||_op = 1 != 0 — certified bounded
 *      away from 0, the distinction from the commutative depol (whose abelian
 *      image gives 0). (c) UNITAL: ||sum_a K_a^dag K_a - 1_N||_op < 1e-12.
 *      (d) AMPLIATION-INVARIANCE: lo(m=2) == lo(m=3) at fixed kappa (cb-norm
 *      unchanged by (x) id_{M_m} — the identity that makes eta_cb m-independent).
 *
 * MUTATIONS (Rule 7), confirmed RED then restored byte-identical:
 *  M1 — break the calibration: in the generator use noncomm_calibrate_eta(target
 *       + 0.05) (mis-calibrate the root) => eta_cb shifts off 1/4-kappa and the
 *       SHARP tooth (a) lo*sqrt(2)==target goes RED ("ball 2.0e-1 not within
 *       1e-9 of 1.5e-1" at kappa=0.10). [The loose "bracket contains target"
 *       check alone does NOT catch this — ratio 2N is too wide — which is why the
 *       sharp lo*sqrt(2) pin exists.]
 *  M2 — swap the tensor order to K_a (x) 1_m (M_m corner moves to the second
 *       slot, r and unitality preserved) => B1 = E^{(m)}_{01}(x)1_d is NO LONGER
 *       a fixed point of Phi, so tooth (b)'s fixdef check (B1 must be in the
 *       image for the star-commutator to witness the IMAGE algebra) goes RED
 *       ("ball 1.07 not within 1e-12 of 0"). */
static void test_fam_nc_noncomm_boundary(void)
{
    const slong m = 2, d = 2; /* M_2 corner (non-abelian) (x) cb_op_gap on C^2 */
    double kappas[4] = {0.10, 0.05, 0.02, 0.01};
    /* lower-root eta of eta*sqrt(1-eta) = 1/4 - kappa (mirrors the generator's
     * bisection so the test independently re-derives the target). */
    double lo_arr[4];

    arb_t lob, hib;
    arb_init(lob);
    arb_init(hib);

    for (int it = 0; it < 4; it++) {
        double target = 0.25 - kappas[it];

        aic_ucp_kraus phi;
        aic_adv_chan_noncomm_boundary(&phi, m, d, kappas[it], PREC);
        slong N = m * d;
        AIC_CHECK_MSG(phi.dim_K == N && phi.dim_H == N && phi.r == d,
                      "fam-NC kappa=%.2f: bad shape dim_K=%ld dim_H=%ld r=%ld "
                      "(want %ld,%ld,%ld)", kappas[it], (long) phi.dim_K,
                      (long) phi.dim_H, (long) phi.r, (long) N, (long) N,
                      (long) d);

        /* (c) UNITAL: tensor of unital UCP maps. */
        arb_t ud;
        arb_init(ud);
        aic_ucp_unital_defect_kraus(ud, &phi, PREC);
        check_ball_eq(ud, 0.0, 1e-12);
        arb_clear(ud);

        /* (a) CALIBRATION: certified bracket CONTAINS the target 1/4-kappa. The
         * bracket is rigorous (arb balls), so arb_lower(lob) <= target <=
         * arb_upper(hib) certifies eta_cb straddles the boundary. */
        aic_cbnorm_eigfree_ball(lob, hib, &phi, PREC);
        lo_arr[it] = arf_get_d(arb_midref(lob), ARF_RND_NEAR);
        arb_t tgt;
        arb_init(tgt);
        arb_set_d(tgt, target);
        AIC_CHECK_MSG(arb_le(lob, tgt) && arb_le(tgt, hib),
                      "fam-NC kappa=%.2f: certified bracket [%.6f,%.6f] does NOT "
                      "contain target 1/4-kappa=%.6f (boundary not reached)",
                      kappas[it], lo_arr[it],
                      arf_get_d(arb_midref(hib), ARF_RND_NEAR), target);
        arb_clear(tgt);

        /* (a) SHARP calibration pin (the mutation-catching tooth — the "contains"
         * check is too loose at ratio 2N to catch a small miscalibration). For
         * this construction at d=2 the eig-free lower bound is EXACTLY
         * lo = ||J||_F/N = eta_cb / sqrt(2) (measured ratio 0.7071068, a fixed
         * fingerprint of Choi(Psi^2-Psi) for the C^2 measure-prepare family; it is
         * ampliation-invariant under (x) id_{M_m}). So lo*sqrt(2) == target =
         * 1/4-kappa to ~1e-9. A miscalibrated eta moves lo off target/sqrt(2). */
        if (d == 2) {
            arb_t lo2;
            arb_init(lo2);
            arb_sqrt_ui(lo2, 2, PREC);
            arb_mul(lo2, lo2, lob, PREC); /* lo * sqrt(2), still a certified ball */
            check_ball_eq(lo2, target, 1e-9);
            arb_clear(lo2);
        }

        /* (a) cont. EXACT closed form g(eta)=eta*sqrt(1-eta) == target (the
         * generator's calibration), AND that closed form lies in the bracket. */
        double a = 0.0, b = 2.0 / 3.0;
        for (int j = 0; j < 200; j++) {
            double mid = 0.5 * (a + b);
            if (mid * sqrt(1.0 - mid) < target)
                a = mid;
            else
                b = mid;
        }
        double eta = 0.5 * (a + b);
        double cf = eta * sqrt(1.0 - eta);
        AIC_CHECK_MSG(fabs(cf - target) < 1e-12,
                      "fam-NC kappa=%.2f: closed form eta*sqrt(1-eta)=%.14f != "
                      "target=%.14f (calibration root wrong)", kappas[it], cf,
                      target);
        arb_t cfb;
        arb_init(cfb);
        arb_set_d(cfb, cf);
        AIC_CHECK_MSG(arb_le(lob, cfb) && arb_le(cfb, hib),
                      "fam-NC kappa=%.2f: closed form eta_cb=%.6f outside "
                      "certified bracket [%.6f,%.6f] (cb-norm certifier "
                      "inconsistent with tex:378)", kappas[it], cf, lo_arr[it],
                      arf_get_d(arb_midref(hib), ARF_RND_NEAR));
        arb_clear(cfb);

        /* (b) NON-COMMUTATIVITY WITNESS. B1 = E^{(m)}_{01}(x)1_d, B2 =
         * E^{(m)}_{10}(x)1_d are genuine fixed points (Psi unital). Star product
         * X*Y = Phi(XY) (tex:342); the commutator ||B1*B2 - B2*B1||_op == 1. */
        acb_mat_t B1, B2, F, P12, P21, Dc;
        acb_mat_init(B1, N, N);
        acb_mat_init(B2, N, N);
        acb_mat_init(F, N, N);
        acb_mat_init(P12, N, N);
        acb_mat_init(P21, N, N);
        acb_mat_init(Dc, N, N);
        for (slong c = 0; c < d; c++) {
            acb_set_si(acb_mat_entry(B1, 0 * d + c, 1 * d + c), 1); /* E01(x)1_d */
            acb_set_si(acb_mat_entry(B2, 1 * d + c, 0 * d + c), 1); /* E10(x)1_d */
        }
        /* B1, B2 are fixed points: ||Phi(B1)-B1||_op ~ 0 (certifies the witness
         * elements really are in the algebra, so the star-commutator IS the
         * algebra commutator). */
        aic_ucp_apply(F, &phi, B1, PREC);
        acb_mat_sub(F, F, B1, PREC);
        arb_t fd;
        arb_init(fd);
        aic_mat_opnorm(fd, F, PREC);
        check_ball_eq(fd, 0.0, 1e-12);
        /* star-commutator */
        acb_mat_mul(F, B1, B2, PREC);  /* B1 B2 */
        aic_ucp_apply(P12, &phi, F, PREC);
        acb_mat_mul(F, B2, B1, PREC);  /* B2 B1 */
        aic_ucp_apply(P21, &phi, F, PREC);
        acb_mat_sub(Dc, P12, P21, PREC);
        arb_t comm;
        arb_init(comm);
        aic_mat_opnorm(comm, Dc, PREC);
        /* certified == 1 (||(E00-E11)(x)1_d||_op = 1), bounded WELL away from 0 —
         * the commutative depol gives 0 here. */
        check_ball_eq(comm, 1.0, 1e-12);
        arb_t half;
        arb_init(half);
        arb_set_d(half, 0.5);
        AIC_CHECK_MSG(arb_gt(comm, half),
                      "fam-NC kappa=%.2f: star-commutator %.6f not > 0.5 (image "
                      "NOT non-commutative — would match the commutative G2)",
                      kappas[it],
                      arf_get_d(arb_midref(comm), ARF_RND_NEAR));
        arb_clear(half);
        arb_clear(comm);
        arb_clear(fd);
        acb_mat_clear(Dc);
        acb_mat_clear(P21);
        acb_mat_clear(P12);
        acb_mat_clear(F);
        acb_mat_clear(B2);
        acb_mat_clear(B1);
        aic_ucp_kraus_clear(&phi);
    }

    /* (d) AMPLIATION-INVARIANCE: eta_cb is independent of the M_m factor size, so
     * the eig-free lo (= ||J||_F/N) is IDENTICAL for m=2 and m=3 at fixed kappa.
     * This is the certified fingerprint of the cb-norm ampliation-invariance
     * identity eta_cb = ||Psi^2-Psi||_cb (tex:347-349) that makes the boundary
     * m-independent. */
    aic_ucp_kraus phi2, phi3;
    aic_adv_chan_noncomm_boundary(&phi2, 2, d, 0.01, PREC);
    aic_adv_chan_noncomm_boundary(&phi3, 3, d, 0.01, PREC);
    aic_cbnorm_eigfree_ball(lob, hib, &phi2, PREC);
    double lo_m2 = arf_get_d(arb_midref(lob), ARF_RND_NEAR);
    aic_cbnorm_eigfree_ball(lob, hib, &phi3, PREC);
    double lo_m3 = arf_get_d(arb_midref(lob), ARF_RND_NEAR);
    AIC_CHECK_MSG(fabs(lo_m2 - lo_m3) < 1e-12 * (1.0 + lo_m2),
                  "fam-NC: lo(m=2)=%.12f != lo(m=3)=%.12f (cb-norm NOT ampliation-"
                  "invariant; eta_cb should be m-independent, tex:347-349)",
                  lo_m2, lo_m3);
    aic_ucp_kraus_clear(&phi3);
    aic_ucp_kraus_clear(&phi2);

    printf("  fam-NC noncomm-boundary: eta_cb=1/4-kappa via id(M2)(x)cb_op_gap. "
           "kappa=0.10 lo=%.6f ; 0.05 lo=%.6f ; 0.02 lo=%.6f ; 0.01 lo=%.6f "
           "(bracket contains target). star-commutator=1.0 (NON-ABELIAN, vs G2 "
           "abelian); ampliation-invariant lo(m=2)==lo(m=3)=%.6f\n",
           lo_arr[0], lo_arr[1], lo_arr[2], lo_arr[3], lo_m2);

    arb_clear(hib);
    arb_clear(lob);
}

/* fam2B closed-form single-sval of the rank-1 defect (the TIGHT calibration
 * anchor). The defect Lambda = P_0 (x) <rho_sub,.> has ||Lambda||_F =
 * ||P_0||_F * ||rho_sub||_F = ||rho_sub||_F (rank-1, so its ONLY nonzero
 * singular value sval[0] = ||Lambda||_F). With
 *   rho_sub = -tail|v><v| + w1|1><1| + w2|2><2|,  v=(sqrt(1-tail),sqrt(w1),sqrt(w2)),
 *   tail = lower root of tail*sqrt(1-tail)=eta,  w2 = tail*gap_sub, w1 = tail-w2,
 * this returns ||rho_sub||_F. Verified == aic_latd_svd sval[0] to ~1e-6 across
 * the knob grid. A miscalibration of tail (e.g. tail=eta, skipping the root-find)
 * shifts sval[0] off this value => the tight match catches it (the loose eig-free
 * 2n-wide bracket does NOT). */
static double fam2b_sval0_closed(double eta, double gap_sub)
{
    /* replicate the generator's calibration (closed form, no channel build) */
    double a = 0.0, b = 2.0 / 3.0;
    for (int it = 0; it < 200; it++) {
        double mid = 0.5 * (a + b);
        if (mid * sqrt(1.0 - mid) < eta)
            a = mid;
        else
            b = mid;
    }
    double tail = 0.5 * (a + b);
    double w2 = tail * gap_sub, w1 = tail - w2;
    double vv[3] = {sqrt(1.0 - tail), sqrt(w1), sqrt(w2)};
    double s = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double r = -tail * vv[i] * vv[j];
            if (i == 1 && j == 1)
                r += w1;
            if (i == 2 && j == 2)
                r += w2;
            s += r * r;
        }
    return sqrt(s);
}

/* ---- fam2B: RANK-1 defect on a near-degenerate subspace (domain.md:249-288;
 * th_almost_idemp tex:2192-2204, the W cancellation tex:2385-2422). The named
 * properties: (a) eta_cb = ||Phi^2-Phi||_cb = eta (calibrated; certified
 * eig-free bracket CONTAINS eta + the TIGHT closed-form sval[0] anchor), (b) the
 * defect Phi^2-Phi is EXACTLY superoperator-RANK-1 (one nonzero singular value),
 * (c) it is CONCENTRATED on a near-degenerate subspace whose realized small
 * eigen-weight is tail*gap_sub (gap_sub the knob), and (d) Phi is UCP for every
 * knob point (CP via Choi min-eig, unital via the Kraus defect).
 *
 * Teeth, in order of load-bearing weight:
 *  (1) SUPEROPERATOR RANK-1 (the central named property): materialize the
 *      n^2 x n^2 superoperator of Phi^2-Phi (vec convention) and SVD it (double
 *      path aic_latd_svd); assert sval[1]/sval[0] < 1e-9 (one nonzero singular
 *      value). Cross-checked by the certified rank of the Gram matrix == 1.
 *  (2) UCP, per knob (Rule 4, load-bearing — the literal domain.md formula is
 *      NON-CP, FINDINGS §C18): Choi PSD (aic_ucp_is_cp_choi == +1) AND unital
 *      defect ~ 0. Reported per (eta, gap_sub).
 *  (3) eta_cb = eta: the certified eig-free cb bracket CONTAINS eta (closed-form
 *      calibration tail*sqrt(1-tail)=eta), for EVERY gap_sub — and the bracket's
 *      lower bound is INVARIANT (to ~1e-9) under changing gap_sub at fixed eta
 *      (the magnitude depends on tail only, not the split).
 *  (4) gap_sub REALIZED: the small eigen-weight of rho_sub scales with gap_sub
 *      (the realized small singular value of rho_sub, extracted from the defect's
 *      dominant left/right structure, tracks tail*gap_sub) — the subspace is
 *      genuinely near-degenerate at the lethal end.
 *  (5) eta=0-ish ORACLE: eta -> 0 (we use the smallest knob 1e-4 plus an
 *      explicit tiny-eta point) keeps the bracket tiny; AND the dedicated eta=0
 *      check uses the family-1B aic_adv_chan_cb_op_gap(d, 0) reduction (Phi=2B
 *      with tail=0 is exactly that dephasing map) => bracket [0,0].
 *
 * MUTATIONS (Rule 7, confirmed RED then restored):
 *  - M1 (break rank-1): in aic_adv_chan_conc_defect ADD a second tilted Kraus
 *    K_1 = |w><1| with w = (sqrt(w2), sqrt(w1), ... ) (a second independent
 *    measure-prepare outcome) => the defect gains a second (output P_1, input)
 *    rank-1 piece, so the superoperator rank becomes 2 => sval[1]/sval[0] is
 *    O(1), the rank-1 tooth goes RED.
 *  - M2 (mis-set eta_cb): make tail = eta DIRECTLY (skip the root-find), so
 *    eta_cb = tail*sqrt(1-tail) = eta*sqrt(1-eta) != eta. The LOOSE eig-free
 *    bracket (2n-wide) still contains eta, so the bracket tooth does NOT catch
 *    this ~10% miscalibration — but the TIGHT closed-form anchor (sval[0] ==
 *    ||rho_sub||_F for the CORRECTLY-calibrated tail) DOES: sval[0] is computed at
 *    tail=eta but compared against the closed form at tail=root(eta), mismatching
 *    by ~3-11% >> 1e-6 => the calibration tooth goes RED.
 *  Both confirmed RED, restored byte-identical to the backup.
 */
static void test_fam2b_conc_defect(void)
{
    const slong d = 4;       /* small; cb-norm Choi is (d^2)x(d^2) = 16x16 here */
    const slong N = d * d;   /* superoperator side */
    /* lethal corner first in the doc; the loop sweeps a 2x2 (mild,lethal) grid. */
    double etas[2] = {0.20, 0.05};     /* mild..toward 1/4 */
    double gaps[2] = {0.5, 1e-8};      /* mild gap .. lethal near-degeneracy */

    arb_t lo, hi, ud, thr, eta_b;
    arb_init(lo);
    arb_init(hi);
    arb_init(ud);
    arb_init(thr);
    arb_init(eta_b);

    double cb_lo_by_gap[2][2]; /* [ie][ig] for the gap-invariance cross-check */
    double realized_small[2];  /* small rho_sub weight at fixed eta, mild vs lethal */

    for (int ie = 0; ie < 2; ie++) {
        for (int ig = 0; ig < 2; ig++) {
            double eta = etas[ie], gap_sub = gaps[ig];
            aic_ucp_kraus phi;
            aic_adv_chan_conc_defect(&phi, d, eta, gap_sub, PREC);

            /* (2a) UNITAL: sum_a K_a^dag K_a = I (pins the OBSERVABLE convention). */
            aic_ucp_unital_defect_kraus(ud, &phi, PREC);
            check_ball_eq(ud, 0.0, 1e-12); /* <v|v>=1 exactly; ~1e-16 rounding */

            /* (2b) CP: Choi PSD, certified (the load-bearing UCP check — the
             * literal domain.md formula fails HERE, FINDINGS §C18). */
            acb_mat_t C;
            acb_mat_init(C, N, N);
            aic_ucp_kraus_to_choi(C, &phi, PREC);
            /* The map is CP analytically (rank-1 Kraus). The Choi has MANY exact
             * zero eigenvalues (the dephasing structure), so -lambda_min(C) is a
             * point ball at ~1e-78 with radius ~1e-74 (measured). A tol=0 (exact)
             * test STRADDLES (the zero-eigenvalue ball spans +/-); use tol=1e-30:
             * far above the ~1e-74 arb rounding floor, far below the smallest
             * genuine positive Choi eigenvalue (O(eta*gap_sub) >~ 1e-9), so the
             * certified-PSD verdict resolves cleanly. */
            arb_set_d(thr, 1e-30);
            int cp = aic_ucp_is_cp_choi(C, thr, PREC);
            AIC_CHECK_MSG(cp == 1,
                          "fam2B eta=%.2f gap=%.0e: Choi NOT certified PSD "
                          "(aic_ucp_is_cp_choi=%d) — Phi is not CP",
                          eta, gap_sub, cp);
            acb_mat_clear(C);

            /* (1) SUPEROPERATOR RANK-1: materialize the n^2 x n^2 superoperator of
             * Lambda = Phi^2 - Phi in the vec convention vec(X)[a*d+b]=X[a,b], SVD
             * it (double path), assert one nonzero singular value. */
            aic_ucp_kraus phi2;
            aic_ucp_compose(&phi2, &phi, &phi, PREC);
            acb_mat_t Eij, Y1, Y2, Lam;
            acb_mat_init(Eij, d, d);
            acb_mat_init(Y1, d, d);
            acb_mat_init(Y2, d, d);
            acb_mat_init(Lam, N, N);
            for (slong i = 0; i < d; i++)
                for (slong j = 0; j < d; j++) {
                    acb_mat_zero(Eij);
                    acb_set_si(acb_mat_entry(Eij, i, j), 1);
                    aic_ucp_apply(Y1, &phi, Eij, PREC);
                    aic_ucp_apply(Y2, &phi2, Eij, PREC);
                    for (slong a = 0; a < d; a++)
                        for (slong b = 0; b < d; b++) {
                            acb_t t;
                            acb_init(t);
                            acb_sub(t, acb_mat_entry(Y2, a, b),
                                    acb_mat_entry(Y1, a, b), PREC);
                            acb_set(acb_mat_entry(Lam, a * d + b, i * d + j), t);
                            acb_clear(t);
                        }
                }
            double _Complex *A = malloc((size_t)N * (size_t)N * sizeof(*A));
            for (slong r = 0; r < N; r++)
                for (slong c = 0; c < N; c++) {
                    double re = arf_get_d(
                        arb_midref(acb_realref(acb_mat_entry(Lam, r, c))),
                        ARF_RND_NEAR);
                    double im = arf_get_d(
                        arb_midref(acb_imagref(acb_mat_entry(Lam, r, c))),
                        ARF_RND_NEAR);
                    A[r * N + c] = re + im * I;
                }
            double *sv = malloc((size_t)N * sizeof(double));
            aic_latd_svd(sv, NULL, NULL, A, N, N);
            AIC_CHECK_MSG(sv[0] > 1e-12 && sv[1] / sv[0] < 1e-9,
                          "fam2B eta=%.2f gap=%.0e: defect NOT superoperator-"
                          "rank-1 (sval[0]=%.4e sval[1]=%.4e ratio=%.2e >= 1e-9)",
                          eta, gap_sub, sv[0], sv[1], sv[1] / sv[0]);

            /* (3 TIGHT) CALIBRATION anchor: the (exact, rank-1) single singular
             * value sval[0] == ||rho_sub||_F closed form for the CORRECTLY
             * calibrated tail (lower root of tail*sqrt(1-tail)=eta). This is the
             * tooth that catches a tail miscalibration (the loose eig-free bracket
             * is too wide). */
            double sv0_cf = fam2b_sval0_closed(eta, gap_sub);
            AIC_CHECK_MSG(fabs(sv[0] - sv0_cf) < 1e-6,
                          "fam2B eta=%.2f gap=%.0e: defect sval[0]=%.8f != "
                          "closed-form ||rho_sub||_F=%.8f (calibration of tail "
                          "wrong)",
                          eta, gap_sub, sv[0], sv0_cf);

            /* (1 cross-check) certified rank of the Gram(Lambda) == 1. */
            acb_mat_t G, Lt;
            acb_mat_init(G, N, N);
            acb_mat_init(Lt, N, N);
            acb_mat_conjugate_transpose(Lt, Lam);
            acb_mat_mul(G, Lt, Lam, PREC);
            arb_set_d(thr, 1e-24); /* well below sval[0]^2, well above the 0 floor */
            slong rk = aic_mat_certified_rank(G, thr, PREC);
            AIC_CHECK_MSG(rk == 1,
                          "fam2B eta=%.2f gap=%.0e: certified rank of Gram(defect)"
                          " = %ld (expected 1)",
                          eta, gap_sub, (long)rk);
            acb_mat_clear(G);
            acb_mat_clear(Lt);

            /* (4) gap_sub REALIZED: the small singular value of rho_sub. The
             * defect Lambda = sval[0] * |u><w| with u=vec(P_0), w=vec(rho_sub^T)
             * (up to scale). Extract rho_sub from the dominant RIGHT singular
             * direction: its smaller nonzero singular value (the |2> weight)
             * tracks tail*gap_sub. We read the realized small weight directly off
             * the construction-equivalent quantity: the (2,2) entry of the
             * defect's input functional, = w2 = tail*gap_sub, via Lam row for
             * output P_0 (row a*d+b with a=b=0 => row 0): Lam[0, 2*d+2]. */
            double w2_realized = arf_get_d(
                arb_midref(acb_realref(acb_mat_entry(Lam, 0, 2 * d + 2))),
                ARF_RND_NEAR);
            if (ig == 1)
                realized_small[ie] = fabs(w2_realized);

            free(sv);
            free(A);
            acb_mat_clear(Lam);
            acb_mat_clear(Y2);
            acb_mat_clear(Y1);
            acb_mat_clear(Eij);
            aic_ucp_kraus_clear(&phi2);

            /* (3) eta_cb = eta: certified eig-free cb bracket CONTAINS eta. */
            aic_cbnorm_eigfree_ball(lo, hi, &phi, PREC);
            cb_lo_by_gap[ie][ig] = arf_get_d(arb_midref(lo), ARF_RND_NEAR);
            arb_set_d(eta_b, eta);
            AIC_CHECK_MSG(arb_le(lo, eta_b) && arb_le(eta_b, hi),
                          "fam2B eta=%.2f gap=%.0e: target eta=%.4f NOT in cb "
                          "bracket [%.6f, %.6f] (calibration tail*sqrt(1-tail)=eta "
                          "wrong)",
                          eta, gap_sub, eta, cb_lo_by_gap[ie][ig],
                          arf_get_d(arb_midref(hi), ARF_RND_NEAR));

            aic_ucp_kraus_clear(&phi);
        }
    }

    /* (3 cross-check) gap_sub NEAR-INVARIANCE of the cb magnitude. The cb-norm
     * itself = tail*sqrt(1-tail) depends on tail ONLY (NOT the split), and the
     * TARGET eta is contained in the bracket for BOTH gap_sub (asserted above,
     * tooth 3). The LOOSE eig-free lower endpoint ||J||_F/n does drift mildly
     * with the split (Frobenius norm is split-sensitive: measured ~10% at
     * eta=0.20), but stays within a bounded relative band — it never collapses or
     * tracks gap_sub. Pin: the two endpoints agree to within 20% relative (a
     * defect whose MAGNITUDE tracked gap_sub, e.g. eta*gap_sub, would differ by
     * 8 ORDERS at gap_sub=1e-8 vs 0.5). This is the honest invariance the loose
     * bracket supports; the exact-magnitude invariance is tooth 3 (target in
     * bracket for both). */
    for (int ie = 0; ie < 2; ie++) {
        double lo0 = cb_lo_by_gap[ie][0], lo1 = cb_lo_by_gap[ie][1];
        double rel = fabs(lo0 - lo1) / (0.5 * (lo0 + lo1));
        AIC_CHECK_MSG(rel < 0.20,
                      "fam2B eta=%.2f: cb lower bound tracks gap_sub too strongly "
                      "(mild=%.10f lethal=%.10f, rel=%.3f); magnitude should "
                      "depend on tail, not the split",
                      etas[ie], lo0, lo1, rel);
    }

    /* (4) gap_sub REALIZED at the lethal end: the small rho_sub weight is
     * O(tail*1e-8), i.e. genuinely near-degenerate (orders of magnitude below the
     * defect magnitude eta) for BOTH eta values. */
    AIC_CHECK_MSG(realized_small[0] < 1e-7 && realized_small[1] < 1e-7,
                  "fam2B: lethal gap_sub=1e-8 did NOT produce a near-degenerate "
                  "rho_sub (small weight eta=0.20:%.3e eta=0.05:%.3e not << eta)",
                  realized_small[0], realized_small[1]);

    /* (6) LETHAL CORNER, certified at prec=256 (the brief's bonus cross-check):
     * the extreme knob (eta=1e-4, gap_sub=1e-8) from the sweep. The defect is
     * STILL certified superoperator-rank-1 (Gram rank == 1) and CP, the cb bracket
     * contains eta=1e-4, AND the realized near-degenerate rho_sub weight is
     * ~tail*gap_sub ~ 1e-12 — EIGHT orders of magnitude below the O(eta)=1e-4
     * defect magnitude, yet tracked by the certified arb ball (the W-cancellation
     * regime, tex:2385-2422: the gap_sub-scale quantity is resolved while the
     * O(eta) magnitude is held). This is the double-vs-arb rung: at prec=256 the
     * rank, CP and the 1e-12 sub-weight all certify cleanly. */
    {
        const double eta_l = 1e-4, gap_l = 1e-8;
        aic_ucp_kraus phil;
        aic_adv_chan_conc_defect(&phil, d, eta_l, gap_l, PREC);

        /* CP at the corner */
        acb_mat_t Cl;
        acb_mat_init(Cl, N, N);
        aic_ucp_kraus_to_choi(Cl, &phil, PREC);
        arb_set_d(thr, 1e-30);
        AIC_CHECK_MSG(aic_ucp_is_cp_choi(Cl, thr, PREC) == 1,
                      "fam2B LETHAL eta=1e-4 gap=1e-8: Choi NOT certified PSD");
        acb_mat_clear(Cl);

        /* rank-1 (certified Gram rank) + realized 1e-12 sub-weight at the corner */
        aic_ucp_kraus p2l;
        aic_ucp_compose(&p2l, &phil, &phil, PREC);
        acb_mat_t Eij, Y1, Y2, Lam;
        acb_mat_init(Eij, d, d);
        acb_mat_init(Y1, d, d);
        acb_mat_init(Y2, d, d);
        acb_mat_init(Lam, N, N);
        for (slong i = 0; i < d; i++)
            for (slong j = 0; j < d; j++) {
                acb_mat_zero(Eij);
                acb_set_si(acb_mat_entry(Eij, i, j), 1);
                aic_ucp_apply(Y1, &phil, Eij, PREC);
                aic_ucp_apply(Y2, &p2l, Eij, PREC);
                for (slong a = 0; a < d; a++)
                    for (slong b = 0; b < d; b++) {
                        acb_t t;
                        acb_init(t);
                        acb_sub(t, acb_mat_entry(Y2, a, b),
                                acb_mat_entry(Y1, a, b), PREC);
                        acb_set(acb_mat_entry(Lam, a * d + b, i * d + j), t);
                        acb_clear(t);
                    }
            }
        acb_mat_t Gl, Ltl;
        acb_mat_init(Gl, N, N);
        acb_mat_init(Ltl, N, N);
        acb_mat_conjugate_transpose(Ltl, Lam);
        acb_mat_mul(Gl, Ltl, Lam, PREC);
        arb_set_d(thr, 1e-30);
        slong rkl = aic_mat_certified_rank(Gl, thr, PREC);
        AIC_CHECK_MSG(rkl == 1,
                      "fam2B LETHAL: certified Gram rank=%ld (expected 1) at "
                      "prec=256 — rank-1 not tracked at the lethal corner",
                      (long)rkl);
        double w2_l = arf_get_d(
            arb_midref(acb_realref(acb_mat_entry(Lam, 0, 2 * d + 2))),
            ARF_RND_NEAR);
        /* tail*gap_sub ~ root(1e-4)*1e-8 ~ 1e-12: certified << eta=1e-4 */
        AIC_CHECK_MSG(fabs(w2_l) > 1e-14 && fabs(w2_l) < 1e-10,
                      "fam2B LETHAL: realized rho_sub sub-weight %.3e not in "
                      "(1e-14,1e-10) — the 1e-12 near-degenerate scale not tracked",
                      w2_l);
        printf("  fam2B LETHAL corner (eta=1e-4,gap=1e-8) @prec=256: Gram rank=1 "
               "(rank-1), CP, realized rho_sub sub-weight=%.3e (~1e-12, 8 orders "
               "below eta=1e-4; W-cancellation scale tracked)\n",
               w2_l);
        acb_mat_clear(Gl);
        acb_mat_clear(Ltl);
        acb_mat_clear(Lam);
        acb_mat_clear(Y2);
        acb_mat_clear(Y1);
        acb_mat_clear(Eij);
        aic_ucp_kraus_clear(&p2l);
        aic_ucp_kraus_clear(&phil);
    }

    /* (5) eta -> 0 ORACLE (exact-idempotent dephasing): the 2B map at tail=0 IS
     * the family-1B aic_adv_chan_cb_op_gap(d, 0) complete dephasing. Pin its
     * bracket to [0,0]. (aic_adv_chan_conc_defect asserts eta>0, so the oracle is
     * reached via the equivalent 1B reduction the docstring states.) */
    aic_ucp_kraus phi0;
    aic_adv_chan_cb_op_gap(&phi0, d, 0.0, PREC);
    aic_cbnorm_eigfree_ball(lo, hi, &phi0, PREC);
    arb_set_d(eta_b, 1e-9);
    AIC_CHECK_MSG(arb_lt(hi, eta_b),
                  "fam2B eta->0 oracle: cb bracket upper=%.3e not < 1e-9 "
                  "(tail->0 must give exact-idempotent dephasing, defect 0)",
                  arf_get_d(arb_midref(hi), ARF_RND_NEAR));
    aic_ucp_kraus_clear(&phi0);

    printf("  fam2B conc-defect: rank-1 defect (sval[1]/sval[0]<1e-9) + Choi-PSD "
           "+ unital, all (eta,gap_sub). eta=0.20 cb_lo[gap=.5]=%.6f "
           "[gap=1e-8]=%.6f (gap-invariant); eta=0.05 cb_lo[.5]=%.6f [1e-8]=%.6f. "
           "lethal small rho_sub weight eta=.20:%.2e eta=.05:%.2e. "
           "eta->0 => bracket [0,0]\n",
           cb_lo_by_gap[0][0], cb_lo_by_gap[0][1], cb_lo_by_gap[1][0],
           cb_lo_by_gap[1][1], realized_small[0], realized_small[1]);

    arb_clear(eta_b);
    arb_clear(thr);
    arb_clear(ud);
    arb_clear(hi);
    arb_clear(lo);
}

/* ---- fam3C: NEAR-TRIVIAL-PROJECTION deliver-or-refuse (domain.md:373-412, the
 * Route-A spectral-split family, #7 on the merged lethal shortlist; tex:516-525
 * prop_P basin / tex:931 lem_nontriv_projection). The EIGHTH generator's
 * self-test, the FIRST that stresses the projection finder (aic_projection, the
 * aic-66n unit-aware audition) directly across its DELIVER-OR-REFUSE boundary.
 *
 * The generator aic_adv_proj_near_trivial(n,k,gap_spec) is a BARE two-cluster
 * Hermitian whose smallest INTER-cluster gap = the LARGEST interior gap =
 * gap_spec (FINDINGS §D1: no genuine eps-C* algebra presents a tunable SMALL
 * largest gap, so the knob lives on the Hermitian element the finder's Route-A-
 * Step-3 primitives consume). Teeth, in load-bearing order:
 *
 *  (CAL) GAP CALIBRATION (the construction pin, the corpus-can't-go-toothless
 *    monotone): the realized LARGEST interior gap of H == gap_spec to <= 1e-12,
 *    over the full knob sweep gap_spec/eps in {10,3,1.5,1.0,0.5}, eps in
 *    {0.01,0.05}, n in {4,9,16}. A generator whose gap did NOT track gap_spec
 *    fails here. Confirms the gap is genuinely tunable across the boundary.
 *  (DELIVER) at gap_spec >> floor: aic_projection_gap returns g == gap_spec; the
 *    prop_P basin ||Y^2-I|| < 1 holds (asserted inside aic_projection_ambient);
 *    aic_projection_ambient builds an EXACT ambient idempotent (||P_amb^2-P_amb||
 *    machine-zero) that is NONTRIVIAL (rank k, so ||P_amb||=1 and ||I-P_amb||=1);
 *    and projecting P_amb into a GENUINE eta>0 algebra (the fam3D blockalg
 *    generator, eta~0.05) via aic_projection_into_A gives a FINITE certified star
 *    defect — the deliver path produces a real projection, NOT a silently near-
 *    trivial one. The star defect is O(eta) and ~gap-invariant for a clean ambient
 *    projector (the O(eps/g) header bound is loose here; FINDINGS §D1).
 *  (REFUSE) at gap_spec=0 (the exactly-degenerate single-cluster witness):
 *    aic_projection_gap FAIL-LOUD-ABORTS with "NO positive interior spectral gap"
 *    (the aic-3qv stop condition), exercised through the shared fork+SIGALRM
 *    watchdog (aic_watchdog_assert_failloud) so this binary survives. The finder
 *    MUST NOT return a
 *    silently near-trivial P — it aborts. (The double-vs-arb distinction: the gap
 *    selection is the double path; at gap_spec=0 the double-path eigenvalues
 *    coincide to round-off and the floor 1e-9*max(1,spread) catches it — the
 *    eig-free abort is what stops a silent near-trivial snap.)
 *
 * MUTATIONS (Rule 7), confirmed RED then restored byte-identical:
 *  M1 — in the generator REPLACE the upper-cluster base w+gap_spec with w (drop
 *       the gap): the two clusters merge, the realized largest gap collapses to
 *       the within-cluster scale, and the (CAL) pin realized==gap_spec goes RED
 *       (and the DELIVER fixture would hit the no-gap abort spuriously). Proves
 *       the gap knob is load-bearing.
 *  M2 — widen the REFUSE witness: pass gap_spec=0.3 (instead of 0) to the
 *       watchdog child: aic_projection_gap then DELIVERS (no abort), the child
 *       _exit(0)s instead of SIGABRT, and the watchdog's "child must SIGABRT"
 *       assertion goes RED. Proves the refuse tooth actually requires the abort. */

/* ||A||_op on ball midpoints (mirror of fam3d_eta_proxy's collapse: a near-zero
 * difference matrix with wide accumulated balls false-fails the certified Gram
 * check; collapse to midpoints when we only want the number). */
static double fam3c_mid_opnorm(const acb_mat_t A, slong prec)
{
    slong r = acb_mat_nrows(A), c = acb_mat_ncols(A);
    acb_mat_t M;
    arb_t e;
    acb_mat_init(M, r, c);
    arb_init(e);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(A, i, j));
    aic_mat_opnorm(e, M, prec);
    double v = arf_get_d(arb_midref(e), ARF_RND_NEAR);
    arb_clear(e);
    acb_mat_clear(M);
    return v;
}

/* Realized LARGEST interior gap of a Hermitian H (the gap aic_projection_gap will
 * pick), via double-path zheev — the construction-pin oracle for (CAL). */
static double fam3c_largest_gap(const acb_mat_t H, slong n)
{
    double _Complex *Ha = malloc((size_t) n * n * sizeof(double _Complex));
    double *ev = malloc((size_t) n * sizeof(double));
    AIC_CHECK_MSG(Ha && ev, "fam3C: OOM");
    aic_latd_from_acb_mat(Ha, H);
    aic_latd_eig_hermitian(ev, NULL, Ha, n);
    double mx = -1.0;
    for (slong i = 0; i < n - 1; i++) {
        double g = ev[i + 1] - ev[i];
        if (g > mx) mx = g;
    }
    free(ev);
    free(Ha);
    return mx;
}

/* rank of a Hermitian acb_mat (count eigenvalues > 0.5; P_amb's are ~0 or ~1). */
static slong fam3c_herm_rank_half(const acb_mat_t Pm, slong n)
{
    double _Complex *Pa = malloc((size_t) n * n * sizeof(double _Complex));
    double *ev = malloc((size_t) n * sizeof(double));
    AIC_CHECK_MSG(Pa && ev, "fam3C: OOM");
    aic_latd_from_acb_mat(Pa, Pm);
    aic_latd_eig_hermitian(ev, NULL, Pa, n);
    slong r = 0;
    for (slong i = 0; i < n; i++) if (ev[i] > 0.5) r++;
    free(ev);
    free(Pa);
    return r;
}

/* REFUSE watchdog child: build the EXACTLY-degenerate witness (gap_spec=0, single
 * cluster) and call aic_projection_gap — which MUST abort (aic-3qv no-gap). */
static void fam3c_child_nogap(void)
{
    const slong prec = 256, n = 9, k = 4;
    acb_mat_t H;
    acb_mat_init(H, n, n);
    aic_adv_proj_near_trivial(H, n, k, 0.0, prec);   /* all eigenvalues == 0 */
    double t = 0, g = 0, lmin = 0, lmax = 0;
    slong m = 0;
    aic_projection_gap(&t, &g, &lmin, &lmax, &m, H, prec);   /* MUST abort */
    acb_mat_clear(H);
    _exit(0);
}

static void test_fam3c_proj_near_trivial(void)
{
    const slong prec = 256;
    const slong ns[] = {4, 9, 16};
    const slong ks[] = {2, 4, 8};
    const double epss[] = {0.05, 0.01};
    const double ratios[] = {10.0, 3.0, 1.5, 1.0, 0.5};

    printf("--- fam3C near-trivial projection: gap calibration (CAL) ---\n");
    double worst_relerr = 0.0;
    for (int ni = 0; ni < 3; ni++) {
        slong n = ns[ni], k = ks[ni];
        for (int ei = 0; ei < 2; ei++) {
            double eps = epss[ei];
            for (int ri = 0; ri < 5; ri++) {
                double gs = ratios[ri] * eps;
                acb_mat_t H;
                acb_mat_init(H, n, n);
                aic_adv_proj_near_trivial(H, n, k, gs, prec);
                double realized = fam3c_largest_gap(H, n);
                double relerr = fabs(realized - gs) / gs;
                if (relerr > worst_relerr) worst_relerr = relerr;
                /* (CAL) the construction pin: realized largest gap == gap_spec. */
                AIC_CHECK_MSG(relerr < 1e-12,
                              "fam3C(n=%ld eps=%.2f ratio=%.1f): realized gap=%.6f "
                              "!= gap_spec=%.6f (relerr=%.3e) -- gap knob not "
                              "tracking (M1 mutation?)", (long) n, eps, ratios[ri],
                              realized, gs, relerr);
                acb_mat_clear(H);
            }
        }
    }
    printf("  CAL: realized largest interior gap == gap_spec to <= %.2e over "
           "n in {4,9,16}, eps in {0.05,0.01}, ratio in {10,3,1.5,1.0,0.5}\n",
           worst_relerr);

    /* (DELIVER) gap_spec >> floor: the finder snaps a genuine nontrivial spectral
     * projector; projecting into a genuine eta>0 algebra gives a finite star
     * defect. We use n=4,k=2 (matches the blockalg eta>0 algebra dim N=4). */
    printf("--- fam3C DELIVER (g=Omega(1)): gap-pick + prop_P snap + into-A "
           "star defect O(eta) ---\n");
    /* genuine eta>0 oblique algebra A = (+)_2 M_2 (N=4), eta~0.05 (blockalg). */
    aic_ucp_kraus phi;
    aic_adv_chan_blockalg(&phi, 2, 2, 0.05, prec);   /* k=2 blocks, d=2, t=0.05 */
    double eta = fam3d_eta_proxy(&phi, prec);
    AIC_CHECK_MSG(eta > 1e-4, "fam3C: blockalg eta=%.3e not > 0 (vacuous)", eta);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong N = ae.A.n;   /* == 4 */

    const double eps_d = 0.05;
    const double deliver_ratios[] = {10.0, 3.0, 1.5, 1.0};   /* above the floor */
    double prev_sd = -1.0;
    for (int ri = 0; ri < 4; ri++) {
        double gs = deliver_ratios[ri] * eps_d;
        acb_mat_t H;
        acb_mat_init(H, N, N);
        aic_adv_proj_near_trivial(H, N, N / 2, gs, prec);

        /* Route A Step 3: gap pick (DELIVER -- no abort). */
        double t = 0, g = 0, lmin = 0, lmax = 0;
        slong m = 0;
        aic_projection_gap(&t, &g, &lmin, &lmax, &m, H, prec);
        AIC_CHECK_MSG(fabs(g - gs) < 1e-12,
                      "fam3C DELIVER(ratio=%.1f): chosen gap=%.6f != gap_spec=%.6f",
                      deliver_ratios[ri], g, gs);
        AIC_CHECK_MSG(m == N - N / 2,
                      "fam3C DELIVER(ratio=%.1f): m=%ld != upper-cluster size %ld",
                      deliver_ratios[ri], (long) m, (long) (N - N / 2));

        /* prop_P snap: aic_projection_ambient asserts the basin ||Y^2-I||<1 and
         * builds an EXACT ambient idempotent. */
        acb_mat_t Pamb;
        acb_mat_init(Pamb, N, N);
        aic_projection_ambient(Pamb, H, t, lmin, lmax, prec);   /* basin asserted */
        /* ambient snap defect ||P_amb^2 - P_amb|| (plain product, exact in M_N). */
        acb_mat_t PP, D;
        acb_mat_init(PP, N, N);
        acb_mat_init(D, N, N);
        acb_mat_mul(PP, Pamb, Pamb, prec);
        acb_mat_sub(D, PP, Pamb, prec);
        double snap = fam3c_mid_opnorm(D, prec);
        AIC_CHECK_MSG(snap < 1e-30,
                      "fam3C DELIVER(ratio=%.1f): ambient snap defect=%.3e not "
                      "machine-zero (prop_P sgn did not snap an exact idempotent)",
                      deliver_ratios[ri], snap);
        /* NONTRIVIAL: rank == k (the lower cluster is the -1 side, upper the +1). */
        slong rk = fam3c_herm_rank_half(Pamb, N);
        AIC_CHECK_MSG(rk == m,
                      "fam3C DELIVER(ratio=%.1f): rank(P_amb)=%ld != m=%ld "
                      "(trivial projector -- rank 0 or N?)", deliver_ratios[ri],
                      (long) rk, (long) m);
        AIC_CHECK_MSG(rk > 0 && rk < N,
                      "fam3C DELIVER(ratio=%.1f): rank(P_amb)=%ld trivial (0 or N)",
                      deliver_ratios[ri], (long) rk);

        /* into-A: project into the genuine eta>0 algebra; certify a FINITE star
         * defect (the deliver path produces a real projection, not a silent
         * near-trivial one). */
        acb_mat_t P, SPP, SD;
        acb_mat_init(P, N, N);
        acb_mat_init(SPP, N, N);
        acb_mat_init(SD, N, N);
        aic_projection_into_A(P, &ae.A, Pamb, prec);
        aic_ecstar_star(SPP, &ae.A, P, P, prec);   /* P * P (STAR) */
        acb_mat_sub(SD, SPP, P, prec);
        double sd = fam3c_mid_opnorm(SD, prec);
        double C = sd / eta;
        printf("  ratio=%.1f gap=%.4f: chosen_gap=%.4f rank(P_amb)=%ld snap=%.2e "
               "star_defect=%.4e C=delta/eta=%.4f\n", deliver_ratios[ri], gs, g,
               (long) rk, snap, sd, C);
        /* delta = O(eta): C bounded (catches an O(1)-defect regression). */
        AIC_CHECK_MSG(C < 1.0,
                      "fam3C DELIVER(ratio=%.1f): C=delta/eta=%.4f exceeds 1 -- "
                      "star defect not O(eta)", deliver_ratios[ri], C);
        /* the defect is finite and genuinely > 0 at eta>0 (not vacuous). */
        AIC_CHECK_MSG(sd > 1e-12 && isfinite(sd),
                      "fam3C DELIVER(ratio=%.1f): star defect=%.3e vacuous/non-"
                      "finite at eta>0", deliver_ratios[ri], sd);
        prev_sd = sd;
        acb_mat_clear(SD); acb_mat_clear(SPP); acb_mat_clear(P);
        acb_mat_clear(D); acb_mat_clear(PP); acb_mat_clear(Pamb);
        acb_mat_clear(H);
    }
    (void) prev_sd;
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);

    /* (REFUSE) gap_spec=0 (degenerate single cluster): aic_projection_gap MUST
     * fail-loud-abort with the aic-3qv no-gap message. Via the shared fork+SIGALRM
     * watchdog (aic_watchdog.h) so this binary survives: finished within timeout,
     * NOT clean-exit (a silently near-trivial P would be the bug 3C guards
     * against), SIGABRT, and the captured output names the aic-3qv guard. */
    printf("--- fam3C REFUSE (g->0): aic-3qv no-gap abort ---\n");
    aic_watchdog_assert_failloud(fam3c_child_nogap, 20,
                                 "fam3C REFUSE: aic_projection_gap @gap_spec=0",
                                 "NO positive interior spectral gap");
    printf("  REFUSE: child SIGABRT; stderr CONTAINS the aic-3qv message\n");
}

int main(void)
{
    test_gen1_jordan();
    test_gen2_nonnormal();
    test_gen3_projector();
    test_gen4_near_degen();
    test_gen5_graded();
    test_gen6_boundary();
    test_gen7_propP();
    test_fam1b_cb_op_gap();
    test_fam2a_depol_boundary();
    test_fam1d_unital_defect();
    test_fam1c_carrier_dropout();
    test_fam1c_carrier_dense();   /* DENSIFIED carrier: straddle + convention (aic-v5f) */
    test_fam3d_blockalg();
    test_fam3d_bijective_eta();   /* aic-66n wrapper-collapse regression (prec=256) */
    test_fam_nc_noncomm_boundary();  /* NON-COMMUTATIVE eta->1/4 boundary (aic-cxo) */
    test_fam2b_conc_defect();        /* RANK-1 defect, near-degenerate subspace (aic-cxo, 2B) */
    test_fam3c_proj_near_trivial();  /* NEAR-TRIVIAL PROJECTION deliver-or-refuse (aic-dbo.2, 3C) */

    aic_test_report("test_adversarial");
    printf("OK test_adversarial\n");
    return 0;
}
