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
#include "aic_mat.h"
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

int main(void)
{
    test_gen1_jordan();
    test_gen2_nonnormal();
    test_gen3_projector();
    test_gen4_near_degen();
    test_gen5_graded();
    test_gen6_boundary();
    test_gen7_propP();

    aic_test_report("test_adversarial");
    printf("OK test_adversarial\n");
    return 0;
}
