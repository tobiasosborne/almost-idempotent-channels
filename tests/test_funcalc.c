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
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/mag.h>

#include "aic_funcalc.h"
/* the certified-domain probe aic_funcalc_int_def_X2 is a funcalc internal (the
 * public wrappers abort); the test drives it directly like test_funcalc_global.c. */
#include "../src/aic_funcalc_internal.h"
#include "aic_mat.h"
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

int main(void)
{
    test_identities();
    test_candidate_agreement();
    test_prop_P();
    test_degenerate_projector();
    test_precision_ladder();
    test_wide_radius();
    test_out_of_basin_failloud();

    aic_test_report("test_funcalc");
    printf("OK test_funcalc\n");
    return 0;
}
