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
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic/aic_cbnorm.h"
#include "aic/aic_assoc.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic_test.h"

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
    test_fam3d_blockalg();

    aic_test_report("test_adversarial");
    printf("OK test_adversarial\n");
    return 0;
}
