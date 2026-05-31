/* test_mat.c — cross-checks for the Layer-0 matrix substrate (beads aic-9zh,
 * "T-mat"). CLAUDE.md Rule 5/6: every check asserts a value or an inter-routine
 * agreement; unit cases catch typos, the reconstruction/identity/precision-ladder
 * cross-checks catch algorithmic errors.
 *
 * Checks:
 *  1. exact special cases: Frobenius and operator norm of I_n and of an explicit
 *     diagonal matrix equal the hand-computed values.
 *  2. eig reconstruction: a known Hermitian H decomposes to (V, lambda) with
 *     V diag(lambda) V^dag ~ H (eigenvectors normalised to unit 2-norm so V is
 *     unitary for the simple Hermitian spectrum).
 *  3. precision ladder (beads aic-5ty): operator norm and eigenvalues at prec=53
 *     vs prec=256 agree within a 53-bit-appropriate tol.
 *  4. Kronecker / partial-trace identities: Tr_2(A (x) B) = Tr(B) A and
 *     Tr_1(A (x) B) = Tr(A) B for small pseudo-random A, B.
 *  5. singular values of a known matrix equal hand values.
 *  6. aic_mat_herm_max_eig NON-FINITE-eps FALLBACK (bead aic-92f fix pass): on a
 *     near-zero Hermitian Gram whose acb_mat_eig_global_enclosure radius is
 *     non-finite, the eig-free ball [-||H||_F, ||H||_F] must be SOUND (finite,
 *     radius == ||H||_F) and CONTAIN the true lambda_max (cross-checked against
 *     the LAPACK double path). Mutation-proven: arb_zero(out) (unsound point-ball)
 *     and 0.25*||H||_F (too-small radius) both fail to contain lambda_max.
 */
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"   /* evil-matrix corpus (aic-dbo.2): gen5 graded_diag */
#include "aic/aic_assoc.h"
#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"
#include "../src/aic_mat_internal.h" /* aic_mat_int_is_hermitian (aic-2yo teeth test) */

/* tol as an acb (zero imaginary) for the matrix close-checker, and as arb. */
static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* One LCG step -> a small generic integer in [-9, 9]. Advances *state. Defined
 * (single sequenced modification of *state per call). */
static slong lcg_draw(unsigned long *state)
{
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    return (slong) ((*state >> 40) % 19) - 9;
}

/* Assert an arb ball x is within tol of the arb reference ref. Subtraction at
 * 300 bits so the comparison precision is never the bottleneck. */
static void check_arb_close(const arb_t x, const arb_t ref, const arb_t tol)
{
    arb_t d;
    arb_init(d);
    arb_sub(d, x, ref, 300);
    arb_abs(d, d);
    AIC_CHECK_MSG(arb_le(d, tol), "arb value not within tol of reference");
    arb_clear(d);
}

/* Convenience: compare against an EXACTLY-representable double value v (integers
 * and small dyadic rationals are exact as doubles, so arb_set_d is lossless). */
static void check_arb_eq(const arb_t x, double v, const arb_t tol)
{
    arb_t t;
    arb_init(t);
    arb_set_d(t, v);
    check_arb_close(x, t, tol);
    arb_clear(t);
}

/* Normalise each column of V to unit 2-norm in place (so V is unitary when its
 * columns are an orthogonal eigenbasis of a Hermitian matrix). */
static void normalise_columns(acb_mat_t V, slong prec)
{
    slong n = acb_mat_nrows(V);
    arb_t nrm, sq;
    acb_t scal;
    arb_init(nrm);
    arb_init(sq);
    acb_init(scal);
    for (slong c = 0; c < n; c++) {
        arb_zero(nrm);
        for (slong r = 0; r < n; r++) {
            acb_abs(sq, acb_mat_entry(V, r, c), prec);
            arb_mul(sq, sq, sq, prec);
            arb_add(nrm, nrm, sq, prec);
        }
        arb_sqrt(nrm, nrm, prec);
        acb_set_arb(scal, nrm);
        for (slong r = 0; r < n; r++)
            acb_div(acb_mat_entry(V, r, c), acb_mat_entry(V, r, c), scal, prec);
    }
    acb_clear(scal);
    arb_clear(sq);
    arb_clear(nrm);
}

/* --- 1. exact norm special cases --- */
static void test_norms_exact(void)
{
    const slong prec = 256;
    arb_t tol, val;
    arb_init(tol);
    arb_init(val);
    set_tol(tol, 1e-60);

    /* Identity 5x5: ||Id||_F = sqrt(5), ||Id||_op = 1. sqrt(5) is irrational, so
     * the reference is computed in arb at the same prec, not as a double. (Named
     * `Id`, not `I`: aic_latd.h pulls in <complex.h>, whose `I` macro would
     * shadow a matrix called `I`.) */
    acb_mat_t Id;
    acb_mat_init(Id, 5, 5);
    acb_mat_one(Id);
    arb_t sqrt5;
    arb_init(sqrt5);
    arb_sqrt_ui(sqrt5, 5, prec);
    aic_mat_frobenius_norm(val, Id, prec);
    check_arb_close(val, sqrt5, tol);
    arb_clear(sqrt5);
    aic_mat_opnorm(val, Id, prec);
    check_arb_eq(val, 1.0, tol);
    acb_mat_clear(Id);

    /* Diagonal diag(3, -4, 12): ||.||_F = 13, ||.||_op = 12. */
    acb_mat_t D;
    acb_mat_init(D, 3, 3);
    acb_set_si(acb_mat_entry(D, 0, 0), 3);
    acb_set_si(acb_mat_entry(D, 1, 1), -4);
    acb_set_si(acb_mat_entry(D, 2, 2), 12);
    aic_mat_frobenius_norm(val, D, prec);
    check_arb_eq(val, 13.0, tol); /* sqrt(9+16+144) */
    aic_mat_opnorm(val, D, prec);
    check_arb_eq(val, 12.0, tol);
    acb_mat_clear(D);

    arb_clear(val);
    arb_clear(tol);
}

/* Build a fixed 3x3 Hermitian H with known structure (real diagonal, conjugate
 * off-diagonals). */
static void build_hermitian3(acb_mat_t H)
{
    acb_mat_init(H, 3, 3);
    acb_set_si(acb_mat_entry(H, 0, 0), 2);
    acb_set_si(acb_mat_entry(H, 1, 1), 5);
    acb_set_si(acb_mat_entry(H, 2, 2), -1);
    acb_set_si_si(acb_mat_entry(H, 0, 1), 1, 2);  /* 1 + 2i */
    acb_set_si_si(acb_mat_entry(H, 1, 0), 1, -2); /* conj */
    acb_set_si_si(acb_mat_entry(H, 0, 2), 3, -1); /* 3 - i */
    acb_set_si_si(acb_mat_entry(H, 2, 0), 3, 1);  /* conj */
    acb_set_si_si(acb_mat_entry(H, 1, 2), 0, 4);  /* 4i */
    acb_set_si_si(acb_mat_entry(H, 2, 1), 0, -4); /* conj */
}

/* --- 2. eig reconstruction cross-check --- */
static void test_eig_reconstruct(void)
{
    const slong prec = 256;
    const slong n = 3;
    acb_mat_t H;
    build_hermitian3(H);

    arb_ptr ev = _arb_vec_init(n);
    acb_mat_t V, Vd, D, T, R;
    acb_mat_init(V, n, n);
    acb_mat_init(Vd, n, n);
    acb_mat_init(D, n, n);
    acb_mat_init(T, n, n);
    acb_mat_init(R, n, n);

    aic_mat_eig_hermitian(ev, V, H, prec);
    normalise_columns(V, prec);

    /* D = diag(ev) */
    for (slong i = 0; i < n; i++)
        acb_set_arb(acb_mat_entry(D, i, i), ev + i);

    /* R = V D V^dag */
    acb_mat_conjugate_transpose(Vd, V);
    acb_mat_mul(T, V, D, prec);
    acb_mat_mul(R, T, Vd, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-30); /* prec=256 reconstruction is very tight */
    AIC_CHECK_ACB_MAT_CLOSE(R, H, tol);
    arb_clear(tol);

    acb_mat_clear(R);
    acb_mat_clear(T);
    acb_mat_clear(D);
    acb_mat_clear(Vd);
    acb_mat_clear(V);
    _arb_vec_clear(ev, n);
    acb_mat_clear(H);
}

/* --- 3. precision-ladder cross-check (acb@53 vs acb@256) --- */
static void test_precision_ladder(void)
{
    const slong n = 3;
    acb_mat_t H;
    build_hermitian3(H);

    arb_ptr e53 = _arb_vec_init(n);
    arb_ptr e256 = _arb_vec_init(n);
    aic_mat_eig_hermitian(e53, NULL, H, 53);
    aic_mat_eig_hermitian(e256, NULL, H, 256);

    arb_t tol, d;
    arb_init(tol);
    arb_init(d);
    set_tol(tol, 1e-12); /* 53-bit-appropriate */
    for (slong i = 0; i < n; i++) {
        arb_sub(d, e53 + i, e256 + i, 64);
        arb_abs(d, d);
        AIC_CHECK_MSG(arb_le(d, tol),
                      "eig[%ld] disagrees between prec=53 and prec=256",
                      (long) i);
    }

    /* operator norm at both precisions */
    arb_t op53, op256;
    arb_init(op53);
    arb_init(op256);
    aic_mat_opnorm(op53, H, 53);
    aic_mat_opnorm(op256, H, 256);
    arb_sub(d, op53, op256, 64);
    arb_abs(d, d);
    AIC_CHECK_MSG(arb_le(d, tol), "operator norm disagrees 53 vs 256");

    arb_clear(op256);
    arb_clear(op53);
    arb_clear(d);
    arb_clear(tol);
    _arb_vec_clear(e256, n);
    _arb_vec_clear(e53, n);
    acb_mat_clear(H);
}

/* --- 4. Kronecker / partial-trace identities --- */
static void test_kron_ptrace(void)
{
    const slong prec = 200;
    const slong a = 2, b = 3;
    acb_mat_t A, B, K, T2, T1, expA, expB;
    acb_mat_init(A, a, a);
    acb_mat_init(B, b, b);
    acb_mat_init(K, a * b, a * b);
    acb_mat_init(T2, a, a);
    acb_mat_init(T1, b, b);
    acb_mat_init(expA, a, a);
    acb_mat_init(expB, b, b);

    /* Deterministic pseudo-random complex fill via a tiny LCG. Generic entries
     * (no structural coincidence such as a shifted diagonal summing to the
     * trace) so the Tr_2(A(x)B)=Tr(B)A / Tr_1(A(x)B)=Tr(A)B identities are
     * mutation-sensitive — verified by perturbing the partial-trace index. Each
     * draw is sequenced into its own statement (no two LCG steps in one
     * unsequenced expression) to keep the behaviour defined. */
    unsigned long seed = 0x9e3779b9UL;
    for (slong i = 0; i < a; i++)
        for (slong j = 0; j < a; j++) {
            slong re = lcg_draw(&seed), im = lcg_draw(&seed);
            acb_set_si_si(acb_mat_entry(A, i, j), re, im);
        }
    for (slong i = 0; i < b; i++)
        for (slong j = 0; j < b; j++) {
            slong re = lcg_draw(&seed), im = lcg_draw(&seed);
            acb_set_si_si(acb_mat_entry(B, i, j), re, im);
        }

    aic_mat_kronecker(K, A, B, prec);
    aic_mat_partial_trace_right(T2, K, a, b, prec); /* expect Tr(B) * A */
    aic_mat_partial_trace_left(T1, K, a, b, prec);  /* expect Tr(A) * B */

    acb_t trA, trB;
    acb_init(trA);
    acb_init(trB);
    acb_mat_trace(trA, A, prec);
    acb_mat_trace(trB, B, prec);
    acb_mat_scalar_mul_acb(expA, A, trB, prec);
    acb_mat_scalar_mul_acb(expB, B, trA, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-40);
    AIC_CHECK_ACB_MAT_CLOSE(T2, expA, tol);
    AIC_CHECK_ACB_MAT_CLOSE(T1, expB, tol);
    arb_clear(tol);

    acb_clear(trB);
    acb_clear(trA);
    acb_mat_clear(expB);
    acb_mat_clear(expA);
    acb_mat_clear(T1);
    acb_mat_clear(T2);
    acb_mat_clear(K);
    acb_mat_clear(B);
    acb_mat_clear(A);
}

/* --- 5. singular values of a known matrix --- */
static void test_singular_values(void)
{
    const slong prec = 256;
    /* 2x3 real matrix [[3,0,0],[0,4,0]]: singular values 4, 3 (descending). */
    acb_mat_t A;
    acb_mat_init(A, 2, 3);
    acb_set_si(acb_mat_entry(A, 0, 0), 3);
    acb_set_si(acb_mat_entry(A, 1, 1), 4);

    arb_ptr sv = _arb_vec_init(2);
    aic_mat_singular_values(sv, A, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-50);
    check_arb_eq(sv + 0, 4.0, tol);
    check_arb_eq(sv + 1, 3.0, tol);
    arb_clear(tol);

    _arb_vec_clear(sv, 2);
    acb_mat_clear(A);

    /* operator norm of the same A must equal its largest singular value, 4. */
    acb_mat_t A2;
    acb_mat_init(A2, 2, 3);
    acb_set_si(acb_mat_entry(A2, 0, 0), 3);
    acb_set_si(acb_mat_entry(A2, 1, 1), 4);
    arb_t op, tol2;
    arb_init(op);
    arb_init(tol2);
    set_tol(tol2, 1e-50);
    aic_mat_opnorm(op, A2, prec);
    check_arb_eq(op, 4.0, tol2);
    arb_clear(tol2);
    arb_clear(op);
    acb_mat_clear(A2);
}

/* --- 6. aic_mat_herm_max_eig non-finite-eps fallback (bead aic-92f) --------- */
/* The fallback in aic_mat_herm_max_eig (src/aic_mat_spectral.c) returns the
 * rigorous eig-free ball [-||H||_F, ||H||_F] when acb_mat_eig_global_enclosure
 * yields a NON-FINITE radius. That branch is exercised by test_assoc on near-zero
 * Grams of EXACT idempotents but had no direct correctness test; the hostile
 * review showed mutating its body to arb_zero(out) (an UNSOUND point-ball claiming
 * lambda_max==0 exactly) passed the whole suite. This test gives it teeth.
 *
 * TRIGGER (found empirically): the Hermitian Gram H = D^dag D, D = S^2 - S, where
 * S is the superoperator of the EXACT idempotent noiseless_subsystem(2,2) channel.
 * Its entries are ~1e-31, the heuristic-QR-seeded global enclosure overflows, and
 * eps is non-finite at prec=256 -> the fallback fires (verified: mid==0, radius ==
 * ||H||_F == 6.287e-32 exactly). The true lambda_max (LAPACK) is 3.144e-32 =
 * 0.5*||H||_F, strictly inside the fallback ball AND strictly outside 0.25*||H||_F.
 * Plain near-zero diagonals/rank-1 matrices do NOT fire the branch (FLINT returns
 * a finite eps for them) -- the degenerate Gram structure is what overflows.
 *
 * MUTATION-PROVEN (FIX 1 of the aic-92f fix pass):
 *  (a) fallback body -> arb_zero(out): ball [0,0] cannot contain lambda_max>0 -> RED.
 *  (b) fallback radius -> 0.25*||H||_F: ball cannot contain lambda_max=0.5*||H||_F
 *      -> RED.
 *  Both verified RED, then restored.  This test ALSO fails loud (AIC_CHECK at the
 *  "fallback fired" guard) if a future FLINT stops triggering the branch, so it can
 *  never silently stop exercising the fallback. */

/* noiseless_subsystem(dL,dE) Heisenberg Kraus, replicated inline (test_idemp.h is
 * not included here -- its un-attributed static builders would warn under
 * -Wunused-function). Phi(X) = (Tr_E X) (x) (1_E / dE); K_{jk} = (1_L (x)
 * |e_j><e_k|)/sqrt(dE), dE^2 operators, on C^{dL*dE}. */
static void build_noiseless_subsystem(aic_ucp_kraus *phi, slong dL, slong dE)
{
    slong n = dL * dE;
    aic_ucp_kraus_init(phi, n, n, dE * dE);
    double inv = 1.0 / sqrt((double) dE);
    slong op = 0;
    for (slong j = 0; j < dE; j++)
        for (slong k = 0; k < dE; k++) {
            for (slong a = 0; a < dL; a++)
                acb_set_d(acb_mat_entry(phi->K[op], a * dE + j, a * dE + k), inv);
            op++;
        }
}

static void test_herm_max_eig_fallback(void)
{
    const slong prec = 256;

    /* Build the triggering Gram H = D^dag D, D = S^2 - S. */
    aic_ucp_kraus phi;
    build_noiseless_subsystem(&phi, 2, 2);
    slong n = phi.dim_H, N = n * n;
    acb_mat_t S, S2, D, Dd, H;
    acb_mat_init(S, N, N);
    acb_mat_init(S2, N, N);
    acb_mat_init(D, N, N);
    acb_mat_init(Dd, N, N);
    acb_mat_init(H, N, N);
    aic_assoc_superop_from_ucp(S, &phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);            /* ~1e-31 entries */
    acb_mat_conjugate_transpose(Dd, D);
    acb_mat_mul(H, Dd, D, prec);            /* Hermitian PSD Gram */

    arb_t lam, fro;
    arb_init(lam);
    arb_init(fro);
    aic_mat_herm_max_eig(lam, H, prec);
    aic_mat_frobenius_norm(fro, H, prec);

    /* (1) The fallback must actually have FIRED: a sound ball with finite, nonzero
     * radius and midpoint 0 (the eig-free [-||H||_F, ||H||_F]). If a future FLINT
     * returns a finite eps here, this guard fails loud rather than letting the
     * branch go silently untested. */
    AIC_CHECK_MSG(arb_is_finite(lam),
                  "fallback: lambda_max ball must be finite (no nan/inf)");
    double mid = arf_get_d(arb_midref(lam), ARF_RND_NEAR);
    double rad = mag_get_d(arb_radref(lam));
    double frod = arf_get_d(arb_midref(fro), ARF_RND_NEAR);
    AIC_CHECK_MSG(frod > 0.0, "fallback: ||H||_F must be positive (near-zero Gram)");
    AIC_CHECK_MSG(mid == 0.0,
                  "fallback did NOT fire: midpoint=%.3e != 0 (FLINT returned a "
                  "finite eps? the fallback branch is no longer exercised)", mid);
    {
        /* radius == ||H||_F to ~mag precision. The fallback sets the radius via
         * arb_add_error(out, ||H||_F-ball), so the stored radius is the UPPER
         * bound of the ||H||_F ball rounded UP into a mag_t (~30-bit mantissa),
         * hence rad >= frod and rad/frod - 1 ~ 3e-9 (the mag granularity), NOT
         * exactly 1. Assert rad in [||H||_F, 1.01*||H||_F]: the lower bound
         * (rad >= ||H||_F) is the SOUND eig-free bound and fails for both the
         * arb_zero mutation (rad=0) and the 0.25*||H||_F mutation (rad~0.25 frod);
         * the upper bound confirms it is the TIGHT [-||H||_F,||H||_F], not a wide
         * over-claim. */
        double ratio = rad / frod;
        AIC_CHECK_MSG(ratio >= 1.0 && ratio <= 1.01,
                      "fallback radius=%.6e not in [||H||_F, 1.01||H||_F]=%.6e "
                      "(ratio=%.9f); the fallback bound is not [-||H||_F,||H||_F]",
                      rad, frod, ratio);
    }

    /* (2) The ball must CONTAIN the true lambda_max. Compute it via the LAPACK
     * double path (the project's independent double-vs-arb oracle rung), which
     * handles the degenerate near-zero spectrum the arb certifier could not. */
    double _Complex *arr = malloc(sizeof(double _Complex) * (size_t) (N * N));
    double *ev = malloc(sizeof(double) * (size_t) N);
    AIC_CHECK_MSG(arr != NULL && ev != NULL, "fallback test: OOM");
    aic_latd_from_acb_mat(arr, H);
    aic_latd_eig_hermitian(ev, NULL, arr, N);  /* ascending */
    double lam_true = ev[N - 1];
    AIC_CHECK_MSG(lam_true > 0.0,
                  "fallback test: expected lambda_max>0 (PSD Gram), got %.3e",
                  lam_true);

    /* lam_true must lie in [mid-rad, mid+rad]. Build the true value as a
     * zero-radius arb and require the fallback ball to certainly contain it
     * (arb_contains is the rigorous enclosure predicate). */
    arb_t lam_true_arb;
    arb_init(lam_true_arb);
    arb_set_d(lam_true_arb, lam_true);
    AIC_CHECK_MSG(arb_contains(lam, lam_true_arb),
                  "fallback ball does NOT contain true lambda_max=%.6e "
                  "(ball mid=%.6e rad=%.6e) -- UNSOUND (arb_zero point-ball "
                  "mutation, or too-small-radius mutation, would fail here)",
                  lam_true, mid, rad);

    /* Sanity that the test has real teeth on BOTH mutations: lam_true is strictly
     * between 0.25*||H||_F (the too-small mutation) and ||H||_F. Assert it so the
     * mutation analysis is self-documenting and a regression that shifts lam_true
     * out of this band is flagged. */
    AIC_CHECK_MSG(lam_true > 0.25 * frod,
                  "fallback test integrity: lambda_max=%.6e <= 0.25*||H||_F=%.6e, "
                  "the 0.25*||H||_F mutation would NOT be caught", lam_true,
                  0.25 * frod);
    AIC_CHECK_MSG(lam_true < frod,
                  "fallback test integrity: lambda_max=%.6e >= ||H||_F=%.6e",
                  lam_true, frod);

    printf("  test6 herm_max_eig fallback: ||H||_F=%.3e fired (mid=0,rad=||H||_F); "
           "lambda_max(LAPACK)=%.3e in ball, > 0.25*||H||_F=%.3e\n",
           frod, lam_true, 0.25 * frod);

    arb_clear(lam_true_arb);
    free(ev);
    free(arr);
    arb_clear(fro);
    arb_clear(lam);
    acb_mat_clear(H);
    acb_mat_clear(Dd);
    acb_mat_clear(D);
    acb_mat_clear(S2);
    acb_mat_clear(S);
    aic_ucp_kraus_clear(&phi);
}

/* --- 7. assert_hermitian RELATIVE tol on graded/ill-conditioned Gram (aic-2yo) - */
/* The old ABSOLUTE Hermiticity tol 2^-(prec-8) false-failed on a Gram whose entry
 * radii grow with magnitude (the asymmetry it flagged was the arb ball RADIUS, not
 * a real asymmetry), so aic_mat_opnorm / aic_mat_singular_values ABORTED on any
 * ill-conditioned input (cond >~ 1e2, precision-independently). The relative +
 * absolute-floor tol tol*(1+|H_ij|+|H_ji|) fixes it while keeping teeth on a
 * GENUINELY non-Hermitian input. Instance drawn from the adversarial corpus
 * (aic_adv_graded_diag = gen5, the family that surfaced the bug).
 *
 * MUTATION-PROVEN: reverting the relative scale to the bare absolute tol (drop the
 * arb_add/arb_add_si/arb_mul, compare mag <= tol) makes (1) and (2) RED (the graded
 * Gram is rejected and opnorm aborts). Verified RED, then restored. */
static void test_assert_hermitian_relative(void)
{
    const slong prec = 256;
    const slong n = 6;
    const double range = 1e8;

    /* (1) THE FIX: the Gram G = D^dag D of a graded diagonal (kappa=1e8) is a
     * genuine Hermitian; is_hermitian must ACCEPT it (the old absolute tol -> 0). */
    acb_mat_t D, Dd, G;
    acb_mat_init(D, n, n);
    acb_mat_init(Dd, n, n);
    acb_mat_init(G, n, n);
    aic_adv_graded_diag(D, n, range, prec);     /* diag(1, r, ..., range) */
    acb_mat_conjugate_transpose(Dd, D);
    acb_mat_mul(G, Dd, D, prec);                /* G_ii ~ r^{2i}, up to range^2=1e16 */
    AIC_CHECK_MSG(aic_mat_int_is_hermitian(G, prec),
                  "aic-2yo: graded Gram (kappa=%.0e) rejected as non-Hermitian "
                  "(absolute-tol regression)", range);

    /* (2) END-TO-END: aic_mat_opnorm(D) routes through Gram + herm_max_eig +
     * assert_hermitian; it must SUCCEED (no abort) and enclose sigma_max(D)=range
     * (D diagonal). Before the fix this aborted. Containment is the rigorous check;
     * the 1% midpoint band also catches a 0-straddling fallback ball (~5e7). */
    arb_t op, range_arb;
    arb_init(op);
    arb_init(range_arb);
    arb_set_d(range_arb, range);
    aic_mat_opnorm(op, D, prec);
    double opd = arf_get_d(arb_midref(op), ARF_RND_NEAR);
    AIC_CHECK_MSG(arb_is_finite(op) && arb_contains(op, range_arb),
                  "aic-2yo: opnorm(graded_diag) ball does not enclose range=%.0e "
                  "(mid=%.6e)", range, opd);
    AIC_CHECK_MSG(fabs(opd - range) <= 1e-2 * range,
                  "aic-2yo: opnorm(graded_diag)=%.6e far from range=%.0e", opd, range);

    /* (3) TEETH: GENUINELY non-Hermitian inputs must STILL be rejected. Real
     * asymmetric [[1,2],[3,1]] (H_01=2 != conj(H_10)=3) and symmetric-but-not-
     * Hermitian complex [[1,i],[i,1]] (H_01=i, conj(H_10)=-i). */
    acb_mat_t NH;
    acb_mat_init(NH, 2, 2);
    acb_set_si(acb_mat_entry(NH, 0, 0), 1);
    acb_set_si(acb_mat_entry(NH, 1, 1), 1);
    acb_set_si(acb_mat_entry(NH, 0, 1), 2);
    acb_set_si(acb_mat_entry(NH, 1, 0), 3);
    AIC_CHECK_MSG(!aic_mat_int_is_hermitian(NH, prec),
                  "aic-2yo TEETH: real asymmetric [[1,2],[3,1]] wrongly accepted "
                  "(relative tol too loose)");
    acb_onei(acb_mat_entry(NH, 0, 1));          /* H_01 = i */
    acb_onei(acb_mat_entry(NH, 1, 0));          /* H_10 = i  => conj(H_10) = -i */
    AIC_CHECK_MSG(!aic_mat_int_is_hermitian(NH, prec),
                  "aic-2yo TEETH: symmetric-not-Hermitian [[1,i],[i,1]] wrongly "
                  "accepted");

    /* (4) SANITY: a genuine Hermitian (real symmetric) is accepted. */
    acb_set_si(acb_mat_entry(NH, 0, 1), 2);
    acb_set_si(acb_mat_entry(NH, 1, 0), 2);
    AIC_CHECK_MSG(aic_mat_int_is_hermitian(NH, prec),
                  "aic-2yo: real symmetric [[1,2],[2,1]] wrongly rejected");

    printf("  test7 assert_hermitian relative tol (aic-2yo): graded Gram kappa=%.0e "
           "accepted, opnorm=%.4e (=range); [[1,2],[3,1]] & [[1,i],[i,1]] rejected\n",
           range, opd);

    arb_clear(range_arb);
    arb_clear(op);
    acb_mat_clear(NH);
    acb_mat_clear(G);
    acb_mat_clear(Dd);
    acb_mat_clear(D);
}

/* --- 8. aic_mat_opnorm Gram radius-accumulation false-fail (bead aic-qgs) ---- */
/* The Gram path of aic_mat_opnorm / aic_mat_singular_values previously routed the
 * RAW acb_mat_mul output G = A^dag A through the certified-eig Hermiticity
 * predicate. When A is built by a deep matmul chain that CANCELS — the canonical
 * case is D = X^2 - I with X^2 ~ I (so D is SMALL in magnitude ~1e-3) while X
 * carries a LARGE accumulated radius (~1e-70, from the oblique-star products in
 * the corner machinery) — the Gram entries are small-magnitude (~1e-5) yet carry
 * a large accumulated radius (~1e-72). The aic-2yo magnitude-relative Hermiticity
 * floor tol*(1+|G_ij|+|G_ji|) is then ~tol (since |G| << 1), ~1e3x SMALLER than
 * the radius, so the predicate FALSE-FAILS and aic_mat_opnorm SIGABRTs on a
 * matrix that is genuinely Hermitian. This was a BLOCKER for the corner machinery
 * (aic_corner_Co / dim_S) on genuinely-oblique (eta>0) S_P wrappers (FINDINGS
 * §C5/§C10; the §C10 measurement ||(2M-1)^2-I|| ~1e-3, deeply in-basin).
 *
 * FIX (aic_mat_gram, src/aic_mat_norms.c): split the certified Gram into an
 * EXACTLY-Hermitian midpoint matrix Gmid (zero radius -> the predicate passes for
 * the right reason) plus a rigorous Weyl perturbation bound R >= ||G_true - Gmid||;
 * run the certified eig on Gmid and inflate by R. NOT symmetrization alone: the
 * off-diagonal midpoints of A^dag A are ALREADY exact conjugates, so the residual
 * the predicate flags is purely the difference-ball RADIUS, which (G+G^dag)/2 does
 * not remove (the diagnosis in the bead missed this; recorded in the report).
 *
 * This test reproduces the trigger at unit level: X a near-involution with a
 * SEEDED ~1e-70 radius, D = X^2 - I (small magnitude, large radius), then assert
 * aic_mat_opnorm(D) COMPLETES (no abort) and its midpoint matches the independent
 * midpoint-route ||mid(D)||_op (the §C5 workaround value) to tolerance, with a
 * sound ball that contains it.
 *
 * MUTATION-PROVEN (recorded in the increment report): reverting aic_mat_gram to
 * hand the raw Gram to herm_max_eig (or to mere (G+G^dag)/2 symmetrization)
 * re-introduces the SIGABRT here — verified by hand, restored. The pre-fix probe
 * tests/test_qgs_probe.c aborted with exactly this message; this is its permanent
 * unit-level form. */
static double opnorm_mid_route(const acb_mat_t A, slong prec)
{
    /* ||mid(A)||_op: drop the accumulated radii (the §C5 workaround at test layer),
     * an independent value to cross-check the certified opnorm midpoint against. */
    slong r = acb_mat_nrows(A), c = acb_mat_ncols(A);
    acb_mat_t M;
    acb_mat_init(M, r, c);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(A, i, j));
    arb_t op;
    arb_init(op);
    aic_mat_opnorm(op, M, prec);
    double o = arf_get_d(arb_midref(op), ARF_RND_NEAR);
    arb_clear(op);
    acb_mat_clear(M);
    return o;
}

static void test_opnorm_gram_radius_qgs(void)
{
    const slong prec = 256;
    const slong n = 3;

    /* X: a near-involution (X^2 ~ I to ~1e-3), with a SEEDED ~1e-70 absolute
     * radius on every entry (mimics the deep oblique-star radius accumulation the
     * corner machinery feeds in). */
    acb_mat_t X;
    acb_mat_init(X, n, n);
    acb_set_d(acb_mat_entry(X, 0, 0), 1.0);
    acb_set_d(acb_mat_entry(X, 1, 1), -1.0);
    acb_set_d(acb_mat_entry(X, 2, 2), 1.0);
    acb_set_d(acb_mat_entry(X, 0, 1), 1e-3);
    acb_set_d(acb_mat_entry(X, 1, 0), 1e-3);
    acb_set_d(acb_mat_entry(X, 0, 2), 2e-3);
    acb_set_d(acb_mat_entry(X, 2, 0), -2e-3);
    acb_set_d(acb_mat_entry(X, 1, 2), 3e-3);
    acb_set_d(acb_mat_entry(X, 2, 1), 3e-3);
    mag_t r;
    mag_init(r);
    mag_set_d(r, 1e-70);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            arb_add_error_mag(acb_realref(acb_mat_entry(X, i, j)), r);
            arb_add_error_mag(acb_imagref(acb_mat_entry(X, i, j)), r);
        }
    mag_clear(r);

    /* D = X^2 - I: small magnitude (~4e-3), large accumulated radius (~2e-70).
     * Its Gram D^dag D has small-magnitude off-diagonals (~1e-5) with radius
     * ~1e-72 -- the exact aic-qgs false-fail trigger (this IS what
     * aic_funcalc_int_def_X2 forms and op-norms). */
    acb_mat_t X2, D;
    acb_mat_init(X2, n, n);
    acb_mat_init(D, n, n);
    acb_mat_sqr(X2, X, prec);
    acb_mat_one(D);
    acb_mat_sub(D, X2, D, prec); /* D = X^2 - I */

    /* (1) THE FIX: aic_mat_opnorm(D) must COMPLETE (no SIGABRT) and be finite.
     * Before the fix this aborted in the Gram-Hermiticity predicate. */
    arb_t op;
    arb_init(op);
    aic_mat_opnorm(op, D, prec); /* must not abort */
    double opd = arf_get_d(arb_midref(op), ARF_RND_NEAR);
    AIC_CHECK_MSG(arb_is_finite(op),
                  "aic-qgs: opnorm(X^2-I) ball must be finite (no nan/inf)");

    /* (2) VALUE: the certified midpoint matches the independent midpoint-route
     * ||mid(D)||_op (the §C5 workaround value) to tolerance. */
    double op_ref = opnorm_mid_route(D, prec);
    AIC_CHECK_MSG(op_ref > 1e-4,
                  "aic-qgs test integrity: ||mid(D)||_op=%.3e too small (the input "
                  "is not exercising a non-trivial opnorm)", op_ref);
    AIC_CHECK_MSG(fabs(opd - op_ref) <= 1e-9 * (op_ref + 1.0),
                  "aic-qgs: opnorm(X^2-I)=%.6e disagrees with ||mid(D)||_op=%.6e",
                  opd, op_ref);

    /* (3) BALL SANITY: the certified opnorm ball must be TIGHT — the Weyl
     * inflation R (~1e-72 here) must not blow up the result. Assert the ball
     * radius is tiny relative to the midpoint (a wide/0-straddling ball, e.g. from
     * a non-finite-eps fallback or a runaway R, would fail this). */
    double op_rad = mag_get_d(arb_radref(op));
    AIC_CHECK_MSG(op_rad <= 1e-9 * (opd + 1.0),
                  "aic-qgs: opnorm ball radius=%.3e too large relative to mid=%.6e "
                  "(Weyl inflation R should be ~machine level here)", op_rad, opd);

    /* (4) singular_values shares the Gram path: it must also complete and its
     * largest singular value must equal the opnorm (both = sqrt(lambda_max)). */
    arb_ptr sv = _arb_vec_init(n);
    aic_mat_singular_values(sv, D, prec); /* must not abort */
    double sv0 = arf_get_d(arb_midref(sv + 0), ARF_RND_NEAR);
    AIC_CHECK_MSG(arb_is_finite(sv + 0) && fabs(sv0 - opd) <= 1e-9 * (opd + 1.0),
                  "aic-qgs: singular_values[0]=%.6e != opnorm=%.6e", sv0, opd);
    _arb_vec_clear(sv, n);

    printf("  test8 opnorm Gram radius (aic-qgs): opnorm(X^2-I)=%.6e completes "
           "(was SIGABRT); =||mid(D)||_op=%.6e; sigma_max=%.6e\n", opd, op_ref, sv0);

    arb_clear(op);
    acb_mat_clear(D);
    acb_mat_clear(X2);
    acb_mat_clear(X);
}

int main(void)
{
    test_norms_exact();
    test_eig_reconstruct();
    test_precision_ladder();
    test_kron_ptrace();
    test_singular_values();
    test_herm_max_eig_fallback();
    test_assert_hermitian_relative();
    test_opnorm_gram_radius_qgs();

    aic_test_report("test_mat");
    printf("OK test_mat\n");
    return 0;
}
