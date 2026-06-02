/* test_funcalc.c — cross-checks for the funcalc module (beads aic-4bg,
 * "T-funcalc"). CLAUDE.md Rule 5/6: every check asserts a value, a defining
 * identity, an inter-candidate agreement, or the paper's O(delta) bound; the
 * degenerate-projector case is the load-bearing cross-check distinguishing the
 * eig-free candidates from an eig route (bead aic-w4o.1).
 *
 * Checks:
 *  1. defining identities (tex:514-528): sgn(X)^2 = I, |X|^2 = X^2,
 *     theta(X)^2 = theta(X).
 *  2. candidate agreement: aic_sgn_newton_schulz and aic_sgn_denman_beavers
 *     agree within tol on a random in-domain X (catches an algorithmic error in
 *     either route).
 *  3. prop_P (tex:524-533): near-idempotent P -> Ptilde^2 = Ptilde, commutation
 *     ||Ptilde P - P Ptilde|| small, and the bound ||Ptilde - P|| <= ||2P-I|| C delta.
 *  4. DEGENERATE input: an exact rank-k projector gives X = 2P-I with X^2 = I
 *     exactly; both eig-free candidates return sgn(X) = X within tol. This is the
 *     spectrum the certified simple-eig path CANNOT handle.
 *  5. precision ladder (bead aic-5ty): aic_sgn at prec=53 vs prec=256 agree.
 *
 *  6. WIDE-RADIUS in-basin (bead aic-1vp, FINDINGS §C11): a near-involution X with
 *     ||X^2-I||_op << 1 but a seeded WIDE entry radius (~1e-70, the upstream corner
 *     star-product inherited radius). Pre-fix aic_sgn ABORTED ("no convergence in
 *     100 iters") on this genuinely in-basin input; post-fix it returns a certified
 *     sign (||Y^2-I|| and ||YX-XY|| small, a prec-floor sanity-backstop tol (see certify_sign)). Teeth: a
 *     genuinely OUT-of-basin X (||X^2-I||>1) still FAILS LOUD (forked SIGABRT — the
 *     midpoint fix must not mask a real out-of-basin input).
 *
 *  7. ADVERSARIAL DUAL-OUTCOME (bead aic-dbo.4, the FIRST increment establishing the
 *     reusable dbo.4 pattern; CLAUDE.md adversarial-first + Rule 4). Each funcalc
 *     routine is DRAWN against the adversarial corpus (tests/aic_adversarial.h) and,
 *     per instance, asserted to do EITHER (a) deliver a CERTIFIED output meeting the
 *     paper bound (in-basin), OR (b) FAIL LOUD via SIGABRT (out-of-basin). "Passes on
 *     nice inputs" is NOT tested. The reusable idiom is fc_assert_dual (a draw +
 *     dual-outcome harness) + the fc_watch fork+SIGALRM watchdog (mirror of
 *     test_xo0_failloud.c / test_adversarial.c v5f), copyable by the projection /
 *     contraction / ucp retrofits. The instances + MEASURED behavior (prec=256):
 *       - gen6 boundary_x2I (nla 7a, tex:807-841): s=+0.5,+0.01 IN-BASIN -> aic_sgn
 *         delivers sgn^2=I (~2.4e-74), [X,sgn]=0, AND the tex:520 bound
 *         ||sgn(X)-X|| <= ||X||*||X^2-I|| holds (ratio 0.37, 0.29). s=-0.01 OUT
 *         (||X^2-I||_op=1.01, rho(I-X^2)=1.01) -> aic_sgn AND aic_sgn_newton_schulz
 *         BOTH fail loud. NOTE (Rule 2 correction): the dbo.4 brief claimed gen6
 *         s<0 gives rho(I-X^2)=|s|<1 so the default aic_sgn would auto-dispatch to
 *         global Newton and DELIVER; that is WRONG. gen6 is X=diag(sqrt(2-s),1,..),
 *         so I-X^2 = diag(s-1,0,..) and rho(I-X^2)=|s-1|=1-s=1.01 for s=-0.01 — OUT
 *         of the SPECTRAL basin too, so aic_sgn(s=-0.01) FAILS LOUD (MEASURED SIGABRT).
 *       - gen1 jordan_t13 (nla 1a, tex:540-544): t->0 (1e-9). MEASURED ||X^2-I||~1.62,
 *         rho(I-X^2)~1.11 for ALL tested t incl t=1 -> aic_sgn / aic_theta FAIL LOUD
 *         (the paper's t^{1/3} fragility lands the whole companion matrix outside the
 *         funcalc domain; nla.md:74-78 "domain guard fires for small t").
 *       - gen4 near_degen_herm (nla 4a, tex:466-489) wired AS prop_P input: a diagonal
 *         near-projector with eigenvalues near {0,1} a tiny `gap` apart (delta=
 *         (gap/2)(1-gap/2) < 1/4). aic_prop_P DELIVERS an exact idempotent (Pt^2=Pt
 *         ~1e-75) for gap down to 1e-9 (MEASURED) — the eig-free route handles the
 *         near-degenerate spectrum the simple-eig path cannot.
 *       - gen5 graded_diag (nla 5c, tex:667-700): range<sqrt2 IN-BASIN. aic_sgn
 *         delivers sgn(D)=I (~e-74) up to range=1.41. SEPARATE fail-loud tooth:
 *         aic_abs (the binomial xpow series) FAILS LOUD at range>=1.25 (q=||X^2-I||
 *         ~0.56, still < 1 so IN the tex:511 domain) because the (X^2)^{1/2} series
 *         needs > 200 terms — a Rule-4 cap abort, NOT a domain abort (distinct from
 *         the sgn quadratic route, which delivers at the same input). MEASURED.
 *     FINDING: none. Every out-of-basin instance FAILS LOUD; every in-basin instance
 *     meets its certified bound; the double-vs-arb (prec=53 vs 256) rung agrees to
 *     ~1.5e-14. No silent wrong answer, no bound violation surfaced.
 *
 * Mutation-proof (CLAUDE.md Rule 7): the cross-checks were confirmed RED by
 * temporarily replacing the Newton-Schulz update 1/2 Y(3I - Y^2) with the wrong
 * coefficient 1/2 Y(2I - Y^2). The corrupted iteration no longer converges to a
 * fixed point, so the convergence guard aborted the suite (exit 134) before the
 * identity checks even ran — the regression is caught loudly. Coefficient
 * restored to 3; the suite is GREEN again. Stated here per Rule 7.
 *
 * Mutation-proof for check 6 (bead aic-1vp): (a) REVERTING the midpoint fix (set
 * Y_0 = X via acb_mat_set instead of aic_funcalc_int_to_midpoint in
 * aic_sgn_newton_schulz) makes test_wide_radius RED — the 100-iter cap fires and the
 * suite aborts (exit 134), confirmed by hand 2026-05-31. (b) the out-of-basin tooth
 * (forked) fires SIGABRT; weakening the certificate (e.g. arb_lt -> arb_le with an
 * inflated tol, or dropping the ||YX-XY|| arm) would let a non-sign pass — recorded
 * RED by hand. Restored; suite GREEN.
 */
/* _POSIX_C_SOURCE: the fork+SIGALRM watchdog (bead aic-dbo.4, test 7) uses
 * sigaction/kill/mkstemp/alarm — same feature-test macro as test_xo0_failloud.c. */
#define _POSIX_C_SOURCE 200809L
#include <fcntl.h>
#include <signal.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/mag.h>

#include "aic/aic_funcalc.h"
/* the certified-domain probe aic_funcalc_int_def_X2 is a funcalc internal (the
 * public wrappers abort); the test drives it directly like test_funcalc_global.c. */
#include "../src/aic_funcalc_internal.h"
#include "aic/aic_mat.h"
#include "aic_adversarial.h"        /* the adversarial corpus (bead aic-dbo.4) */
#include "aic_test.h"

static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* X = S + eps*H where S = diag(+/-1) (so S^2 = I) and H a fixed small Hermitian
 * perturbation. For small eps, ||X^2 - I||_op = O(eps) < 1 (in domain), and X is
 * a non-trivial (non-diagonal) sign-like matrix. signs[i] in {+1,-1}. */
static void build_perturbed_sign(acb_mat_t X, slong n, const int *signs,
                                 double eps, slong prec)
{
    acb_mat_zero(X);
    for (slong i = 0; i < n; i++)
        acb_set_si(acb_mat_entry(X, i, i), signs[i]);
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

/* --- 1. defining identities --- */
static void test_identities(void)
{
    const slong prec = 256;
    const slong n = 4;
    int signs[4] = {1, -1, 1, -1};
    acb_mat_t X, S, S2, AbsX, Abs2, X2, Th, Th2, I;
    acb_mat_init(X, n, n);
    acb_mat_init(S, n, n);
    acb_mat_init(S2, n, n);
    acb_mat_init(AbsX, n, n);
    acb_mat_init(Abs2, n, n);
    acb_mat_init(X2, n, n);
    acb_mat_init(Th, n, n);
    acb_mat_init(Th2, n, n);
    acb_mat_init(I, n, n);
    build_perturbed_sign(X, n, signs, 0.05, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-40);

    /* sgn(X)^2 = I */
    aic_sgn(S, X, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_one(I);
    AIC_CHECK_ACB_MAT_CLOSE(S2, I, tol);

    /* |X|^2 = X^2 */
    aic_abs(AbsX, X, prec);
    acb_mat_sqr(Abs2, AbsX, prec);
    acb_mat_sqr(X2, X, prec);
    AIC_CHECK_ACB_MAT_CLOSE(Abs2, X2, tol);

    /* theta(X)^2 = theta(X) */
    aic_theta(Th, X, prec);
    acb_mat_sqr(Th2, Th, prec);
    AIC_CHECK_ACB_MAT_CLOSE(Th2, Th, tol);

    arb_clear(tol);
    acb_mat_clear(I); acb_mat_clear(Th2); acb_mat_clear(Th);
    acb_mat_clear(X2); acb_mat_clear(Abs2); acb_mat_clear(AbsX);
    acb_mat_clear(S2); acb_mat_clear(S); acb_mat_clear(X);
}

/* --- 2. candidate agreement: Newton-Schulz vs Denman-Beavers --- */
static void test_candidate_agreement(void)
{
    const slong prec = 256;
    const slong n = 5;
    int signs[5] = {1, 1, -1, 1, -1};
    acb_mat_t X, Sns, Sdb;
    acb_mat_init(X, n, n);
    acb_mat_init(Sns, n, n);
    acb_mat_init(Sdb, n, n);
    build_perturbed_sign(X, n, signs, 0.08, prec);

    aic_sgn_newton_schulz(Sns, X, prec);
    aic_sgn_denman_beavers(Sdb, X, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-40);
    AIC_CHECK_ACB_MAT_CLOSE(Sns, Sdb, tol);
    arb_clear(tol);

    acb_mat_clear(Sdb); acb_mat_clear(Sns); acb_mat_clear(X);
}

/* --- 3. prop_P: idempotent, commutation, O(delta) bound --- */
static void test_prop_P(void)
{
    const slong prec = 256;
    const slong n = 4;
    /* near-idempotent P = P0 + small Hermitian perturbation, P0 = diag(1,1,0,0).
     * delta = ||P^2 - P|| is O(perturbation) < 1/4. */
    acb_mat_t P, Pt, Pt2, PtP, PPt, comm, diff;
    acb_mat_init(P, n, n);
    acb_mat_init(Pt, n, n);
    acb_mat_init(Pt2, n, n);
    acb_mat_init(PtP, n, n);
    acb_mat_init(PPt, n, n);
    acb_mat_init(comm, n, n);
    acb_mat_init(diff, n, n);
    acb_mat_zero(P);
    acb_set_si(acb_mat_entry(P, 0, 0), 1);
    acb_set_si(acb_mat_entry(P, 1, 1), 1);
    /* small Hermitian perturbation, magnitude ~0.05 (delta well under 1/4) */
    for (slong i = 0; i < n; i++)
        for (slong j = i + 1; j < n; j++) {
            double re = 0.05 * (((i + 2 * j) % 3) - 1);
            double im = 0.05 * (((i + j) % 3) - 1);
            acb_t z;
            acb_init(z);
            acb_set_d_d(z, re, im);
            acb_add(acb_mat_entry(P, i, j), acb_mat_entry(P, i, j), z, prec);
            acb_conj(z, z);
            acb_add(acb_mat_entry(P, j, i), acb_mat_entry(P, j, i), z, prec);
            acb_clear(z);
        }

    aic_prop_P(Pt, P, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-40);

    /* Ptilde^2 = Ptilde (exact idempotent). */
    acb_mat_sqr(Pt2, Pt, prec);
    AIC_CHECK_ACB_MAT_CLOSE(Pt2, Pt, tol);

    /* commutation ||Ptilde P - P Ptilde|| small. */
    acb_mat_mul(PtP, Pt, P, prec);
    acb_mat_mul(PPt, P, Pt, prec);
    acb_mat_sub(comm, PtP, PPt, prec);
    arb_t cnorm;
    arb_init(cnorm);
    aic_mat_opnorm(cnorm, comm, prec);
    AIC_CHECK_MSG(arb_lt(cnorm, tol), "prop_P commutation not within tol");

    /* O(delta) bound: ||Ptilde - P|| <= ||2P - I|| * C * delta, C = 4 (a concrete
     * constant covering the O(.) in tex:530; the proof gives 4 delta < 1 for the
     * domain, so C = 4 is a safe envelope verified to hold here). */
    acb_mat_t P2, defect, twoPmI, Itmp;
    acb_mat_init(P2, n, n);
    acb_mat_init(defect, n, n);
    acb_mat_init(twoPmI, n, n);
    acb_mat_init(Itmp, n, n);
    acb_mat_sqr(P2, P, prec);
    acb_mat_sub(defect, P2, P, prec);
    arb_t delta, normX, lhs, rhs;
    arb_init(delta); arb_init(normX); arb_init(lhs); arb_init(rhs);
    aic_mat_opnorm(delta, defect, prec);          /* delta = ||P^2 - P|| */
    acb_mat_scalar_mul_2exp_si(twoPmI, P, 1);      /* 2P */
    acb_mat_one(Itmp);
    acb_mat_sub(twoPmI, twoPmI, Itmp, prec);       /* 2P - I */
    aic_mat_opnorm(normX, twoPmI, prec);           /* ||2P - I|| */
    acb_mat_sub(diff, Pt, P, prec);
    aic_mat_opnorm(lhs, diff, prec);               /* ||Ptilde - P|| */
    arb_mul(rhs, normX, delta, prec);              /* ||2P-I|| * delta */
    arb_mul_si(rhs, rhs, 4, prec);                 /* * C, C = 4 */
    AIC_CHECK_MSG(arb_le(lhs, rhs),
                  "prop_P O(delta) bound violated: "
                  "||Ptilde-P|| > 4 ||2P-I|| delta");

    arb_clear(rhs); arb_clear(lhs); arb_clear(normX); arb_clear(delta);
    acb_mat_clear(Itmp); acb_mat_clear(twoPmI);
    acb_mat_clear(defect); acb_mat_clear(P2);
    arb_clear(cnorm);
    arb_clear(tol);
    acb_mat_clear(diff); acb_mat_clear(comm); acb_mat_clear(PPt);
    acb_mat_clear(PtP); acb_mat_clear(Pt2); acb_mat_clear(Pt); acb_mat_clear(P);
}

/* --- 4. DEGENERATE input: exact rank-k projector, sgn(2P-I) = 2P-I --- */
static void test_degenerate_projector(void)
{
    const slong prec = 256;
    const slong n = 6;
    const slong rank = 4;
    /* P = diag(1,1,1,1,0,0): exact projector, P^2 = P. X = 2P - I = diag(+/-1)
     * (eigenvalues 0 and 1 each repeated => highly degenerate; X^2 = I exactly).
     * The eig-free candidates must return sgn(X) = X within tol; the certified
     * simple-eig path would ABORT on this spectrum (bead aic-w4o.1). */
    acb_mat_t P, X, Sns, Sdb;
    acb_mat_init(P, n, n);
    acb_mat_init(X, n, n);
    acb_mat_init(Sns, n, n);
    acb_mat_init(Sdb, n, n);
    acb_mat_zero(P);
    for (slong i = 0; i < rank; i++)
        acb_set_si(acb_mat_entry(P, i, i), 1);

    acb_mat_scalar_mul_2exp_si(X, P, 1);   /* 2P */
    acb_mat_t I;
    acb_mat_init(I, n, n);
    acb_mat_one(I);
    acb_mat_sub(X, X, I, prec);            /* 2P - I = diag(1,1,1,1,-1,-1) */

    aic_sgn_newton_schulz(Sns, X, prec);
    aic_sgn_denman_beavers(Sdb, X, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-40);
    AIC_CHECK_ACB_MAT_CLOSE(Sns, X, tol); /* Newton-Schulz: sgn(X) = X */
    AIC_CHECK_ACB_MAT_CLOSE(Sdb, X, tol); /* Denman-Beavers: sgn(X) = X */
    arb_clear(tol);

    acb_mat_clear(I);
    acb_mat_clear(Sdb); acb_mat_clear(Sns); acb_mat_clear(X); acb_mat_clear(P);
}

/* --- 5. precision ladder: aic_sgn at prec=53 vs prec=256 --- */
static void test_precision_ladder(void)
{
    const slong n = 4;
    int signs[4] = {1, -1, -1, 1};
    acb_mat_t X, S53, S256;
    acb_mat_init(X, n, n);
    acb_mat_init(S53, n, n);
    acb_mat_init(S256, n, n);
    build_perturbed_sign(X, n, signs, 0.06, 256);

    aic_sgn(S53, X, 53);
    aic_sgn(S256, X, 256);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12); /* 53-bit-appropriate */
    AIC_CHECK_ACB_MAT_CLOSE(S53, S256, tol);
    arb_clear(tol);

    acb_mat_clear(S256); acb_mat_clear(S53); acb_mat_clear(X);
}

/* --- 6. WIDE-RADIUS in-basin: aic_sgn certifies; out-of-basin fails loud --- */

/* X = 2P - I, P = rank-2 projector on span(e1,e2) + small Hermitian perturbation
 * (eps), then a seeded uniform entry radius `rad` (the inherited corner radius).
 * In-basin: ||X^2-I||_op = O(eps) < 1, but the wide radius pre-fix stalled the
 * Newton-Schulz midpoint step below the prec-tol floor (bead aic-1vp). */
static void build_wide_radius(acb_mat_t X, slong n, double eps, double rad,
                              slong prec)
{
    acb_mat_t P, I;
    acb_mat_init(P, n, n);
    acb_mat_init(I, n, n);
    acb_mat_zero(P);
    acb_set_si(acb_mat_entry(P, 0, 0), 1);
    acb_set_si(acb_mat_entry(P, 1, 1), 1);
    for (slong i = 0; i < n; i++)
        for (slong j = i + 1; j < n; j++) {
            acb_t z;
            acb_init(z);
            acb_set_d_d(z, eps * (((i + 2 * j) % 3) - 1),
                        eps * (((i + j) % 3) - 1));
            acb_add(acb_mat_entry(P, i, j), acb_mat_entry(P, i, j), z, prec);
            acb_conj(z, z);
            acb_add(acb_mat_entry(P, j, i), acb_mat_entry(P, j, i), z, prec);
            acb_clear(z);
        }
    acb_mat_scalar_mul_2exp_si(X, P, 1);   /* 2P */
    acb_mat_one(I);
    acb_mat_sub(X, X, I, prec);            /* 2P - I */
    mag_t r;
    mag_init(r);
    mag_set_d(r, rad);
    acb_mat_add_error_mag(X, r);           /* seed the wide inherited radius */
    mag_clear(r);
    acb_mat_clear(I);
    acb_mat_clear(P);
}

static void test_wide_radius(void)
{
    const slong prec = 256;
    const slong n = 6;
    acb_mat_t X, S, S2, I, XS, SX, comm;
    acb_mat_init(X, n, n);
    acb_mat_init(S, n, n);
    acb_mat_init(S2, n, n);
    acb_mat_init(I, n, n);
    acb_mat_init(XS, n, n);
    acb_mat_init(SX, n, n);
    acb_mat_init(comm, n, n);
    /* eps=0.04 -> ||X^2-I||_op ~ 0.45 (solidly in-basin), rad ~ 2.7e-70 (the
     * upstream corner star-product inherited radius that pre-fix stalled NS). */
    build_wide_radius(X, n, 0.04, 2.7e-70, prec);

    /* the input is genuinely IN-BASIN: ||X^2-I||_op ~ 0.45 < 1 (this is the whole
     * point — pre-fix the in-basin input still ABORTED on the wide radius). */
    arb_t b, one;
    arb_init(b);
    arb_init(one);
    arb_one(one);
    aic_funcalc_int_def_X2(b, X, prec);
    AIC_CHECK_MSG(arb_lt(b, one),
                  "test_wide_radius: fixture not in-basin (||X^2-I||_op not < 1)");
    arb_clear(one);
    arb_clear(b);

    /* GREEN: aic_sgn returns (pre-fix it ABORTED here, 100-iter cap). */
    aic_sgn(S, X, prec);

    /* certified sign: ||sgn^2 - I|| and ||X sgn - sgn X|| small. The tol tracks
     * the input radius (rad * O(sqrt(n))), not 1e-40 — the wide ball legitimately
     * carries ~rad-scale slack (absorbed by the certificate sanity-backstop tol). */
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-66);                   /* prec-floor tol ~2^-128*8*sqrt(6) ~6e-38 dominates the ~2.7e-70 input radius */
    acb_mat_sqr(S2, S, prec);
    acb_mat_one(I);
    AIC_CHECK_ACB_MAT_CLOSE(S2, I, tol);   /* sgn(X)^2 = I */
    acb_mat_mul(XS, X, S, prec);
    acb_mat_mul(SX, S, X, prec);
    acb_mat_sub(comm, XS, SX, prec);
    arb_t cnorm;
    arb_init(cnorm);
    aic_mat_frobenius_norm(cnorm, comm, prec);
    AIC_CHECK_MSG(arb_lt(cnorm, tol),
                  "test_wide_radius: [X, sgn(X)] not within tol (cert arm)");
    arb_clear(cnorm);
    arb_clear(tol);

    acb_mat_clear(comm); acb_mat_clear(SX); acb_mat_clear(XS);
    acb_mat_clear(I); acb_mat_clear(S2); acb_mat_clear(S); acb_mat_clear(X);
}

/* Teeth: a genuinely OUT-of-basin X (||X^2-I||_op = 3 > 1) must FAIL LOUD even
 * with a seeded wide radius — the midpoint fix must NOT mask a real out-of-basin
 * input (the basin guard / Gelfand precondition / certificate still abort).
 * Forked child (mirrors test_corner T9): the child SIGABRTs, the parent asserts. */
static void test_out_of_basin_failloud(void)
{
    pid_t pid = fork();
    if (pid == 0) {
        const slong prec = 256, n = 3;
        acb_mat_t X, S;
        acb_mat_init(X, n, n);
        acb_mat_init(S, n, n);
        acb_mat_zero(X);
        acb_set_d(acb_mat_entry(X, 0, 0), 2.0);    /* X = diag(2,-2,0.5) */
        acb_set_d(acb_mat_entry(X, 1, 1), -2.0);   /* X^2 = diag(4,4,0.25), */
        acb_set_d(acb_mat_entry(X, 2, 2), 0.5);    /* ||X^2-I||_op = 3 > 1 */
        mag_t r;
        mag_init(r);
        mag_set_d(r, 2.7e-70);
        acb_mat_add_error_mag(X, r);               /* same wide radius */
        mag_clear(r);
        aic_sgn(S, X, prec);                       /* must abort (out of basin) */
        _exit(0);                                  /* reached only if NOT aborted */
    }
    int status = 0;
    waitpid(pid, &status, 0);
    AIC_CHECK_MSG(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT,
                  "test_out_of_basin_failloud: out-of-basin X did not abort via "
                  "SIGABRT (WIFSIGNALED=%d, WTERMSIG=%d) — the midpoint fix "
                  "masked a real out-of-basin input", WIFSIGNALED(status),
                  WIFSIGNALED(status) ? WTERMSIG(status) : -1);
}

/* ===================================================================== *
 * --- 7. ADVERSARIAL DUAL-OUTCOME (bead aic-dbo.4) --------------------- *
 * ===================================================================== *
 *
 * THE REUSABLE dbo.4 PATTERN. For an adversarial instance, we must assert
 * EITHER the certified funcalc bound holds (in-basin) OR the routine FAILS
 * LOUD (out-of-basin). Two reusable pieces below — copyable verbatim by the
 * projection / contraction / ucp retrofits:
 *   (A) fc_watch: a fork + SIGALRM watchdog (mirror of test_xo0_failloud.c's
 *       xo0_run_child and test_adversarial.c's v5f_run_child — both are static
 *       to their own .c, so each retrofit carries its own copy; the watchdog is
 *       a HANG backstop, since a funcalc domain guard aborts BEFORE any iteration
 *       it is belt-and-suspenders here, but the pattern is the point).
 *   (B) fc_assert_failloud(fn, who): run fn() in fc_watch and assert it
 *       SIGABRTed (out-of-basin -> Rule-4 abort) and did not hang.
 * The in-basin "bound holds" half is asserted inline per instance (the certified
 * arb identities/bounds), since the certified quantity differs per routine. */

#define FC_WATCH_S 20              /* a child alive past this is a HANG */

typedef void (*fc_child_fn)(void);
static volatile pid_t fc_watch_pid = 0;
static volatile sig_atomic_t fc_timed_out = 0;
static void fc_alarm(int sig)
{
    (void) sig;
    fc_timed_out = 1;
    if (fc_watch_pid > 0) kill(fc_watch_pid, SIGKILL);
}

/* Run fn() in a forked child bounded by a SIGALRM watchdog. Captures the child's
 * stderr into `err` (NUL-terminated, `errsz`). Writes the raw wait status into
 * *status. Returns 1 if the child terminated within the watchdog, 0 if it was
 * killed for hanging (caller asserts on that). */
static int fc_watch(fc_child_fn fn, int *status, char *err, size_t errsz)
{
    char tmpl[] = "/tmp/aic_fc_err_XXXXXX";
    int efd = mkstemp(tmpl);
    AIC_CHECK_MSG(efd >= 0, "fc_watch: mkstemp failed");
    fflush(NULL);
    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "fc_watch: fork failed");
    if (pid == 0) {
        dup2(efd, STDERR_FILENO);
        close(efd);
        fn();
        _exit(0);                  /* reached only if fn did NOT abort */
    }
    fc_watch_pid = pid;
    fc_timed_out = 0;
    struct sigaction sa = {0}, old;
    sa.sa_handler = fc_alarm;
    sigaction(SIGALRM, &sa, &old);
    alarm(FC_WATCH_S);
    int st = 0;
    pid_t w;
    do { w = waitpid(pid, &st, 0); } while (w < 0 && !fc_timed_out);
    alarm(0);
    sigaction(SIGALRM, &old, NULL);
    if (w < 0) waitpid(pid, &st, 0);
    lseek(efd, 0, SEEK_SET);
    ssize_t r = read(efd, err, errsz - 1);
    err[(r > 0) ? (size_t) r : 0] = '\0';
    close(efd);
    unlink(tmpl);
    *status = st;
    return fc_timed_out ? 0 : 1;
}

/* Assert fn() FAILS LOUD: terminates within the watchdog (no hang) via SIGABRT,
 * and its stderr CONTAINS `needle` (the basin-violation message — proves it
 * aborted for the RIGHT reason, not some unrelated crash). `who` labels the case. */
static void fc_assert_failloud(fc_child_fn fn, const char *who, const char *needle)
{
    int status = 0;
    char err[4096];
    int finished = fc_watch(fn, &status, err, sizeof err);
    AIC_CHECK_MSG(finished, "%s: out-of-basin call HUNG past %ds watchdog "
                  "(expected a fast Rule-4 abort)", who, FC_WATCH_S);
    AIC_CHECK_MSG(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT,
                  "%s: expected SIGABRT (out-of-basin fail-loud) but got "
                  "WIFSIGNALED=%d WTERMSIG=%d / WIFEXITED=%d code=%d. stderr: %s",
                  who, WIFSIGNALED(status),
                  WIFSIGNALED(status) ? WTERMSIG(status) : -1,
                  WIFEXITED(status), WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                  err);
    AIC_CHECK_MSG(strstr(err, needle) != NULL,
                  "%s: aborted but WITHOUT the basin-violation message '%s' "
                  "(silent/wrong-reason abort?). stderr: %s", who, needle, err);
}

/* --- 7a. gen6 boundary_x2I: the deliver-vs-fail-loud pair for sgn/theta --- */

/* Out-of-basin children (drawn at module scope so fc_watch's void(void) child can
 * reach them; the knob is baked in per the MEASURED out-of-basin value). */
static void fc_child_sgn_gen6_out(void)
{
    acb_mat_t X, S;
    acb_mat_init(X, 3, 3);
    acb_mat_init(S, 3, 3);
    aic_adv_boundary_x2I(X, 3, -0.01, 256);  /* ||X^2-I||=1.01, rho(I-X^2)=1.01 */
    aic_sgn(S, X, 256);                       /* must FAIL LOUD (both bases out) */
    acb_mat_clear(S);
    acb_mat_clear(X);
}
static void fc_child_ns_gen6_out(void)
{
    acb_mat_t X, S;
    acb_mat_init(X, 3, 3);
    acb_mat_init(S, 3, 3);
    aic_adv_boundary_x2I(X, 3, -0.01, 256);
    aic_sgn_newton_schulz(S, X, 256);         /* op-basin assert must FAIL LOUD */
    acb_mat_clear(S);
    acb_mat_clear(X);
}

static void test_adv_boundary_x2I(void)
{
    const slong prec = 256;
    const slong n = 3;
    /* IN-BASIN s = {0.5, 0.01}: aic_sgn delivers a certified Hermitian-unitary,
     * and the tex:520 bound ||sgn(X)-X|| <= ||X|| O(||X^2-I||) holds with the
     * envelope constant C=1 (MEASURED ratios 0.37, 0.29). */
    double s_in[2] = {0.5, 0.01};
    for (int c = 0; c < 2; c++) {
        acb_mat_t X, S, S2, I, XS, SX, comm, D;
        acb_mat_init(X, n, n);
        acb_mat_init(S, n, n);
        acb_mat_init(S2, n, n);
        acb_mat_init(I, n, n);
        acb_mat_init(XS, n, n);
        acb_mat_init(SX, n, n);
        acb_mat_init(comm, n, n);
        acb_mat_init(D, n, n);
        aic_adv_boundary_x2I(X, n, s_in[c], prec);

        /* assert the fixture really is IN-BASIN (else the "deliver" claim is vacuous). */
        arb_t b, one;
        arb_init(b);
        arb_init(one);
        arb_one(one);
        aic_funcalc_int_def_X2(b, X, prec);
        AIC_CHECK_MSG(arb_lt(b, one),
                      "test_adv_boundary_x2I: s=%g fixture not in-basin", s_in[c]);

        aic_sgn(S, X, prec);

        /* (i) sgn(X)^2 = I (Hermitian-unitary), certified ~1e-74. */
        arb_t tol;
        arb_init(tol);
        set_tol(tol, 1e-60);
        acb_mat_sqr(S2, S, prec);
        acb_mat_one(I);
        AIC_CHECK_ACB_MAT_CLOSE(S2, I, tol);

        /* (ii) [X, sgn(X)] = 0 (sgn commutes with X, tex:518). */
        acb_mat_mul(XS, X, S, prec);
        acb_mat_mul(SX, S, X, prec);
        acb_mat_sub(comm, XS, SX, prec);
        arb_t cnorm;
        arb_init(cnorm);
        aic_mat_frobenius_norm(cnorm, comm, prec);
        AIC_CHECK_MSG(arb_lt(cnorm, tol),
                      "test_adv_boundary_x2I: s=%g [X,sgn(X)] not ~0", s_in[c]);

        /* (iii) the tex:520 bound ||sgn(X)-X|| <= ||X|| * C * ||X^2-I||, C=1.
         * lhs and rhs are certified arb balls; arb_le is rigorous. */
        acb_mat_sub(D, S, X, prec);
        arb_t lhs, nX, nX2I, rhs;
        arb_init(lhs);
        arb_init(nX);
        arb_init(nX2I);
        arb_init(rhs);
        aic_mat_opnorm(lhs, D, prec);            /* ||sgn(X) - X|| */
        aic_mat_opnorm(nX, X, prec);             /* ||X|| */
        aic_funcalc_int_def_X2(nX2I, X, prec);   /* ||X^2 - I|| */
        arb_mul(rhs, nX, nX2I, prec);            /* C=1 envelope */
        AIC_CHECK_MSG(arb_le(lhs, rhs),
                      "test_adv_boundary_x2I: s=%g tex:520 bound VIOLATED "
                      "(||sgn-X|| > ||X|| ||X^2-I||)", s_in[c]);

        /* (iv) double-vs-arb rung (prec=53 vs prec=256), ladder #2. */
        acb_mat_t S53;
        acb_mat_init(S53, n, n);
        aic_sgn(S53, X, 53);
        arb_t tol53;
        arb_init(tol53);
        set_tol(tol53, 1e-12);
        AIC_CHECK_ACB_MAT_CLOSE(S53, S, tol53);

        arb_clear(tol53);
        acb_mat_clear(S53);
        arb_clear(rhs);
        arb_clear(nX2I);
        arb_clear(nX);
        arb_clear(lhs);
        arb_clear(cnorm);
        arb_clear(tol);
        arb_clear(one);
        arb_clear(b);
        acb_mat_clear(D);
        acb_mat_clear(comm);
        acb_mat_clear(SX);
        acb_mat_clear(XS);
        acb_mat_clear(I);
        acb_mat_clear(S2);
        acb_mat_clear(S);
        acb_mat_clear(X);
    }

    /* OUT-OF-BASIN s=-0.01: aic_sgn AND aic_sgn_newton_schulz both FAIL LOUD, but
     * via DIFFERENT guards (the dual-path teeth):
     *   - aic_sgn dispatches s=-0.01 (out of op-basin) to the GLOBAL Newton, which
     *     then fails the SPECTRAL Gelfand certifier (rho(I-X^2)=1.01 !< 1).
     *   - aic_sgn_newton_schulz hits the op-basin domain assert directly.
     * (MEASURED ||X^2-I||_op = 1.01, rho(I-X^2) = 1.01 — out of BOTH bases.) */
    fc_assert_failloud(fc_child_sgn_gen6_out,
                       "aic_sgn(gen6 s=-0.01)",
                       "cannot certify rho(I-X^2)");
    fc_assert_failloud(fc_child_ns_gen6_out,
                       "aic_sgn_newton_schulz(gen6 s=-0.01)",
                       "out of validity domain");
}

/* --- 7b. gen1 jordan_t13: t->0, the tex:540 spectral-fragility domain edge --- */

static void fc_child_sgn_jordan(void)
{
    acb_mat_t X, S;
    acb_mat_init(X, 3, 3);
    acb_mat_init(S, 3, 3);
    aic_adv_jordan_t13(X, 1e-9, 256); /* ||X^2-I||~1.62, rho(I-X^2)~1.11: OUT */
    aic_sgn(S, X, 256);               /* must FAIL LOUD */
    acb_mat_clear(S);
    acb_mat_clear(X);
}
static void fc_child_theta_jordan(void)
{
    acb_mat_t X, T;
    acb_mat_init(X, 3, 3);
    acb_mat_init(T, 3, 3);
    aic_adv_jordan_t13(X, 1e-9, 256);
    aic_theta(T, X, 256);             /* delegates to aic_sgn -> must FAIL LOUD */
    acb_mat_clear(T);
    acb_mat_clear(X);
}

static void test_adv_jordan_t13(void)
{
    /* The paper's t^{1/3} companion matrix (tex:540-544): MEASURED ||X^2-I||~1.62
     * and rho(I-X^2)~1.11 for t down to 1e-9 (and even t=1) — the whole matrix is
     * OUTSIDE the funcalc domain (no in-basin "deliver" half exists for this gen;
     * nla.md:74-78). So this is a pure fail-loud instance for sgn AND theta. */
    fc_assert_failloud(fc_child_sgn_jordan,
                       "aic_sgn(gen1 jordan t=1e-9)",
                       "rho");
    fc_assert_failloud(fc_child_theta_jordan,
                       "aic_theta(gen1 jordan t=1e-9)",
                       "rho");
}

/* --- 7c. gen4 near_degen_herm AS prop_P input: tiny-gap near-projector --- */

static void test_adv_near_degen_prop_P(void)
{
    /* gen4's raw H = diag(1-gap/2,1+gap/2,3,..) is FAR out of the funcalc domain
     * (||H^2-I||~15). Its funcalc relevance (aic_adversarial.h gen4 docstring) is
     * AS a prop_P input: a near-projector P whose spectrum clusters near {0,1} with
     * a tiny `gap`. We build that diagonal P directly (the gen4 near-degenerate
     * PAIR placed across the 0/1 projector spectrum) and assert aic_prop_P DELIVERS
     * an exact idempotent for gap down to 1e-9 — the eig-free Gelfand/sgn route
     * handling the near-degenerate spectrum that the simple-eig path aborts on
     * (nla.md:466-489 4a; tex:524-533 prop_P). */
    const slong prec = 256;
    const slong n = 4;
    double gaps[3] = {0.5, 1e-3, 1e-9};
    for (int c = 0; c < 3; c++) {
        double g = gaps[c];
        acb_mat_t P, Pt, Pt2, P2, defect;
        acb_mat_init(P, n, n);
        acb_mat_init(Pt, n, n);
        acb_mat_init(Pt2, n, n);
        acb_mat_init(P2, n, n);
        acb_mat_init(defect, n, n);
        acb_mat_zero(P);
        /* near-degenerate 0/1 pair (the gen4 pair, mapped onto a projector spectrum)
         * plus an exact {0,1} remainder — delta = ||P^2-P|| = (g/2)(1-g/2) < 1/4. */
        acb_set_d(acb_mat_entry(P, 0, 0), g / 2.0);       /* near 0 */
        acb_set_d(acb_mat_entry(P, 1, 1), 1.0 - g / 2.0); /* near 1 */
        acb_set_d(acb_mat_entry(P, 2, 2), 0.0);
        acb_set_d(acb_mat_entry(P, 3, 3), 1.0);

        /* assert delta < 1/4 (the prop_P hypothesis is genuinely satisfied). */
        acb_mat_sqr(P2, P, prec);
        acb_mat_sub(defect, P2, P, prec);
        arb_t delta, quarter;
        arb_init(delta);
        arb_init(quarter);
        aic_mat_opnorm(delta, defect, prec);
        arb_set_d(quarter, 0.25);
        AIC_CHECK_MSG(arb_lt(delta, quarter),
                      "test_adv_near_degen_prop_P: gap=%g delta !< 1/4", g);

        /* DELIVER: aic_prop_P returns an EXACT idempotent Ptilde^2 = Ptilde. */
        aic_prop_P(Pt, P, prec);
        acb_mat_sqr(Pt2, Pt, prec);
        arb_t tol;
        arb_init(tol);
        set_tol(tol, 1e-60);
        AIC_CHECK_ACB_MAT_CLOSE(Pt2, Pt, tol);   /* idempotent within ~1e-75 */

        /* commutation [Ptilde, P] = 0 (tex:530). */
        acb_mat_t PtP, PPt, comm;
        acb_mat_init(PtP, n, n);
        acb_mat_init(PPt, n, n);
        acb_mat_init(comm, n, n);
        acb_mat_mul(PtP, Pt, P, prec);
        acb_mat_mul(PPt, P, Pt, prec);
        acb_mat_sub(comm, PtP, PPt, prec);
        arb_t cnorm;
        arb_init(cnorm);
        aic_mat_frobenius_norm(cnorm, comm, prec);
        AIC_CHECK_MSG(arb_lt(cnorm, tol),
                      "test_adv_near_degen_prop_P: gap=%g [Pt,P] not ~0", g);

        arb_clear(cnorm);
        acb_mat_clear(comm);
        acb_mat_clear(PPt);
        acb_mat_clear(PtP);
        arb_clear(tol);
        arb_clear(quarter);
        arb_clear(delta);
        acb_mat_clear(defect);
        acb_mat_clear(P2);
        acb_mat_clear(Pt2);
        acb_mat_clear(Pt);
        acb_mat_clear(P);
    }
}

/* --- 7d. gen5 graded_diag: wide dynamic range -> sgn delivers, abs fails loud - */

static void fc_child_abs_graded_out(void)
{
    acb_mat_t D, A;
    acb_mat_init(D, 4, 4);
    acb_mat_init(A, 4, 4);
    aic_adv_graded_diag(D, 4, 1.25, 256); /* q=||X^2-I||~0.56 < 1 but series>200 */
    aic_abs(A, D, 256);                   /* binomial xpow cap -> FAIL LOUD */
    acb_mat_clear(A);
    acb_mat_clear(D);
}

static void test_adv_graded_diag(void)
{
    const slong prec = 256;
    const slong n = 4;
    /* IN-BASIN range = {1.1, 1.3, 1.41} (all < sqrt2): the positive diagonal has
     * sgn(D) = I; aic_sgn delivers it (MEASURED ||sgn-I||~e-74) via the QUADRATIC
     * Newton-Schulz, robustly across the wide dynamic range. We also run the
     * prec=53 vs prec=256 rung (ladder #2): the arb path tracks the dynamic range
     * where a double would lose relative precision in the small entries. */
    double ranges[3] = {1.1, 1.3, 1.41};
    for (int c = 0; c < 3; c++) {
        acb_mat_t D, S, S2, I, Diff;
        acb_mat_init(D, n, n);
        acb_mat_init(S, n, n);
        acb_mat_init(S2, n, n);
        acb_mat_init(I, n, n);
        acb_mat_init(Diff, n, n);
        aic_adv_graded_diag(D, n, ranges[c], prec);

        aic_sgn(S, D, prec);
        arb_t tol;
        arb_init(tol);
        set_tol(tol, 1e-60);
        acb_mat_sqr(S2, S, prec);
        acb_mat_one(I);
        AIC_CHECK_ACB_MAT_CLOSE(S2, I, tol);     /* sgn^2 = I */
        AIC_CHECK_ACB_MAT_CLOSE(S, I, tol);      /* sgn(D) = I (positive diagonal) */

        acb_mat_t S53;
        acb_mat_init(S53, n, n);
        aic_sgn(S53, D, 53);
        arb_t tol53;
        arb_init(tol53);
        set_tol(tol53, 1e-12);
        AIC_CHECK_ACB_MAT_CLOSE(S53, S, tol53);  /* double-vs-arb rung */

        arb_clear(tol53);
        acb_mat_clear(S53);
        arb_clear(tol);
        acb_mat_clear(Diff);
        acb_mat_clear(I);
        acb_mat_clear(S2);
        acb_mat_clear(S);
        acb_mat_clear(D);
    }

    /* FAIL-LOUD: aic_abs at range=1.25. q=||X^2-I||~0.56 is < 1 (IN the tex:511
     * domain), but the (X^2)^{1/2} binomial series needs > 200 terms, so the
     * xpow term-cap fires (a Rule-4 abort, NOT a domain abort) — the routine
     * refuses to return an uncertified result. This is the SAME input on which
     * aic_sgn (quadratic) delivers, isolating the series route's distinct limit. */
    fc_assert_failloud(fc_child_abs_graded_out,
                       "aic_abs(gen5 range=1.25)",
                       "series tail not below tol");
}

int main(void)
{
    test_identities();
    test_candidate_agreement();
    test_prop_P();
    test_degenerate_projector();
    test_precision_ladder();
    test_wide_radius();
    test_out_of_basin_failloud();
    /* --- adversarial dual-outcome (bead aic-dbo.4) --- */
    test_adv_boundary_x2I();
    test_adv_jordan_t13();
    test_adv_near_degen_prop_P();
    test_adv_graded_diag();

    aic_test_report("test_funcalc");
    printf("OK test_funcalc\n");
    return 0;
}
