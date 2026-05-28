/* test_latd.c — cross-checks for the LAPACK double path (beads aic-w4o.5).
 * CLAUDE.md Rule 6: the strongest test here is double(LAPACK) vs arb@prec=53;
 * Rule 5: every check asserts a value or an agreement. Per project policy (bead
 * aic-dbo.4) the inputs are ADVERSARIAL instances from docs/adversarial/nla.md,
 * NOT happy-path matrices.
 *
 * Checks (each cites the adversarial family it draws from):
 *  1. DEGENERATE Hermitian (nla.md Family 4b, exact-degeneracy projector):
 *     X = 2P - I for a rank-k projector has the EXACT repeated spectrum
 *     {-1 (mult n-k), +1 (mult k)}. LAPACKE_zheev returns it; the arb certified
 *     SIMPLE-eig path (aic_mat_eig_hermitian) CANNOT — it aborts on repeated
 *     eigenvalues (bead aic-w4o.1). So this is cross-checked against the
 *     HAND-BUILT spectrum, and it is the fast double path's exclusive territory.
 *  2. NON-NORMAL (nla.md Family 1a, the paper's tex:540 t^{1/3} 3x3; and
 *     Family 2a, a Davies-type bidiagonal): the EIGENVALUES are hypersensitive
 *     (~t^{1/3}) but the SINGULAR VALUES are well conditioned. Operator norm
 *     (largest singular value) and the full singular-value set must agree
 *     double-vs-arb@53 within ~1e-12.
 *  3. NEAR-SINGULAR (nla.md Family 3a/3b): smallest singular value tiny; the
 *     double SVD must report it gracefully (no NaN/abort) and agree with arb@53.
 *  4. random Hermitian double-vs-arb@53 on eigenvalues AND operator norm.
 *  5. round-trip of the acb_mat <-> array conversion (identity within 0).
 *  6. MUTATION PROOF: a deliberately wrong eigenvalue set is REJECTED by the
 *     cross-check macro (proving the check has teeth; CLAUDE.md Rule 7).
 *
 * Convention reminder (include/aic_latd.h): the double arrays are ROW-MAJOR,
 * arr[i*c + j] = entry (i,j).
 */
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_latd.h"
#include "aic_mat.h"
#include "aic_test.h"

static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* Deterministic LCG -> small generic integer in [-4, 4]. */
static slong lcg_draw(unsigned long *state)
{
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    return (slong) ((*state >> 41) % 9) - 4;
}

/* --- 1. degenerate Hermitian: X = 2P - I, P a rank-k projector --- */
/* The arb certified SIMPLE-eig path cannot diagonalise this (repeated
 * eigenvalues {-1,+1}); LAPACKE_zheev can, and that is the point of latd. We
 * cross-check zheev's eigenvalues against the hand-built spectrum. */
static void test_degenerate_hermitian(void)
{
    const slong n = 4, k = 2; /* rank-2 projector in C^4 */

    /* P = diag(1,1,0,0) conjugated by a fixed real orthogonal Q (a single
     * Givens-like mix in the (0,2) and (1,3) planes) so X is NOT diagonal and
     * the degeneracy is genuine, not basis-aligned. Build X = 2P - I directly
     * in arb at high prec, then take midpoints for the double path. */
    acb_mat_t P, Q, Qt, T, X;
    acb_mat_init(P, n, n);
    acb_mat_init(Q, n, n);
    acb_mat_init(Qt, n, n);
    acb_mat_init(T, n, n);
    acb_mat_init(X, n, n);

    acb_mat_zero(P);
    for (slong i = 0; i < k; i++) acb_set_si(acb_mat_entry(P, i, i), 1);

    /* Q: rotation by theta with cos=3/5, sin=4/5 (exact rational, orthogonal)
     * in planes (0,2) and (1,3); identity elsewhere. */
    const slong prec = 200;
    acb_mat_one(Q);
    acb_t c, s, ns;
    acb_init(c);
    acb_init(s);
    acb_init(ns);
    acb_set_si(c, 3);
    acb_div_si(c, c, 5, prec); /* 3/5 */
    acb_set_si(s, 4);
    acb_div_si(s, s, 5, prec); /* 4/5 */
    acb_neg(ns, s);
    /* plane (0,2) */
    acb_set(acb_mat_entry(Q, 0, 0), c);
    acb_set(acb_mat_entry(Q, 2, 2), c);
    acb_set(acb_mat_entry(Q, 0, 2), ns);
    acb_set(acb_mat_entry(Q, 2, 0), s);
    /* plane (1,3) */
    acb_set(acb_mat_entry(Q, 1, 1), c);
    acb_set(acb_mat_entry(Q, 3, 3), c);
    acb_set(acb_mat_entry(Q, 1, 3), ns);
    acb_set(acb_mat_entry(Q, 3, 1), s);

    acb_mat_conjugate_transpose(Qt, Q);
    acb_mat_mul(T, Q, P, prec);
    acb_mat_mul(P, T, Qt, prec); /* P <- Q P Q^dag, still a rank-2 projector */

    /* X = 2P - I */
    acb_mat_scalar_mul_si(X, P, 2, prec);
    for (slong i = 0; i < n; i++)
        acb_sub_si(acb_mat_entry(X, i, i), acb_mat_entry(X, i, i), 1, prec);

    /* double path: zheev on the midpoints. */
    double _Complex *Xa = malloc((size_t) (n * n) * sizeof(double _Complex));
    double evals[4];
    aic_latd_from_acb_mat(Xa, X);
    aic_latd_eig_hermitian(evals, NULL, Xa, n);

    /* hand-built spectrum: {-1, -1, +1, +1} as an arb set. */
    arb_ptr spec = _arb_vec_init(n);
    arb_set_si(spec + 0, -1);
    arb_set_si(spec + 1, -1);
    arb_set_si(spec + 2, 1);
    arb_set_si(spec + 3, 1);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12); /* 53-bit-appropriate */
    AIC_CHECK_EIGSET_CLOSE(evals, spec, n, tol);

    /* also confirm the eigenvectors zheev returns are orthonormal: V^dag V = I
     * (a double-path internal sanity, ladder rung 1). */
    double _Complex *V = malloc((size_t) (n * n) * sizeof(double _Complex));
    double evals2[4];
    aic_latd_eig_hermitian(evals2, V, Xa, n);
    for (slong p = 0; p < n; p++) {
        for (slong q = 0; q < n; q++) {
            double _Complex dot = 0;
            for (slong i = 0; i < n; i++)
                dot += conj(V[i * n + p]) * V[i * n + q];
            double want = (p == q) ? 1.0 : 0.0;
            AIC_CHECK_MSG(cabs(dot - want) < 1e-12,
                          "zheev eigenvectors not orthonormal at (%ld,%ld)",
                          (long) p, (long) q);
        }
    }

    arb_clear(tol);
    _arb_vec_clear(spec, n);
    free(V);
    free(Xa);
    acb_clear(ns);
    acb_clear(s);
    acb_clear(c);
    acb_mat_clear(X);
    acb_mat_clear(T);
    acb_mat_clear(Qt);
    acb_mat_clear(Q);
    acb_mat_clear(P);
}

/* Build the paper's tex:540 companion 3x3: [[0,1,0],[0,0,1],[t,0,0]] into an
 * acb_mat at the given prec (t exact). */
static void build_t13(acb_mat_t A, double t, slong prec)
{
    (void) prec;
    acb_mat_zero(A);
    acb_set_si(acb_mat_entry(A, 0, 1), 1);
    acb_set_si(acb_mat_entry(A, 1, 2), 1);
    acb_set_d(acb_mat_entry(A, 2, 0), t);
}

/* --- 2. non-normal: opnorm + singular values double-vs-arb@53 --- */
static void test_nonnormal_svd(void)
{
    arb_t tol, op_arb;
    arb_init(tol);
    arb_init(op_arb);
    set_tol(tol, 1e-12);

    /* (a) tex:540 t^{1/3} matrix, t = 1e-6: eigenvalues ~1e-2 burst out
     * (1/3-power sensitive, nla.md 1a), but the SINGULAR values are exactly
     * {1, 1, t} (the matrix is t * a permutation on one row + two unit rows),
     * well conditioned. */
    {
        const slong n = 3;
        const double t = 1e-6;
        acb_mat_t A;
        acb_mat_init(A, n, n);
        build_t13(A, t, 200);

        double _Complex *Aa = malloc((size_t) (n * n) * sizeof(double _Complex));
        aic_latd_from_acb_mat(Aa, A);

        /* opnorm: double vs arb@53 */
        double op_d = aic_latd_opnorm(Aa, n, n);
        aic_mat_opnorm(op_arb, A, 53);
        AIC_CHECK_DOUBLE_CLOSE(op_d, op_arb, tol);

        /* singular value SET: double vs arb@53 (arb_mat_singular_values needs a
         * SIMPLE Gram spectrum; for t!=1 the singular values {1,1,t} are NOT
         * distinct, so the arb simple-SVD aborts — exactly the gap latd fills.
         * We therefore cross-check the double SVD set against the hand value
         * {1, 1, t}; the operator-norm double-vs-arb above is the rigorous rung-2
         * anchor for this instance). */
        double sv[3];
        aic_latd_singular_values(sv, Aa, n, n);
        arb_ptr want = _arb_vec_init(n);
        arb_set_si(want + 0, 1);
        arb_set_si(want + 1, 1);
        arb_set_d(want + 2, t);
        AIC_CHECK_EIGSET_CLOSE(sv, want, n, tol);

        _arb_vec_clear(want, n);
        free(Aa);
        acb_mat_clear(A);
    }

    /* (b) Davies-type bidiagonal D_n(c) = diag(1..n) + c*superdiagonal
     * (nla.md 2a): non-normal, eigenvalues {1..n} but eigenvectors
     * ill-conditioned. Singular values are DISTINCT for small c, so the arb
     * simple-SVD path applies and we get a true double-vs-arb@53 set check. */
    {
        const slong n = 5;
        const double cc = 0.3;
        acb_mat_t D;
        acb_mat_init(D, n, n);
        acb_mat_zero(D);
        for (slong i = 0; i < n; i++)
            acb_set_si(acb_mat_entry(D, i, i), i + 1);
        for (slong i = 0; i + 1 < n; i++)
            acb_set_d(acb_mat_entry(D, i, i + 1), cc);

        double _Complex *Da = malloc((size_t) (n * n) * sizeof(double _Complex));
        aic_latd_from_acb_mat(Da, D);

        double op_d = aic_latd_opnorm(Da, n, n);
        aic_mat_opnorm(op_arb, D, 53);
        AIC_CHECK_DOUBLE_CLOSE(op_d, op_arb, tol);

        double sv[5];
        aic_latd_singular_values(sv, Da, n, n);
        arb_ptr sv_arb = _arb_vec_init(n);
        aic_mat_singular_values(sv_arb, D, 53);
        AIC_CHECK_EIGSET_CLOSE(sv, sv_arb, n, tol);

        _arb_vec_clear(sv_arb, n);
        free(Da);
        acb_mat_clear(D);
    }

    arb_clear(op_arb);
    arb_clear(tol);
}

/* --- 3. near-singular: tiny smallest singular value, graceful behavior --- */
static void test_near_singular(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);

    /* X(sigma) = diag(3, 2, 1, sigma) with sigma tiny (nla.md 3a). The full SVD
     * must return svals = {3, 2, 1, sigma} with the tiny value handled
     * gracefully (no NaN/abort). The large values are kept DISTINCT (3,2,1) so
     * the arb simple-SVD reference path applies at prec=53; the degenerate-large
     * variant diag(1,1,1,sigma) is the eig path's territory, not the SVD set
     * reference's, and is exercised in the degenerate-Hermitian test instead. */
    const slong n = 4;
    const double sigma = 1e-9;
    const double diag_large[3] = {3.0, 2.0, 1.0};
    acb_mat_t X;
    acb_mat_init(X, n, n);
    acb_mat_zero(X);
    for (slong i = 0; i + 1 < n; i++)
        acb_set_d(acb_mat_entry(X, i, i), diag_large[i]);
    acb_set_d(acb_mat_entry(X, n - 1, n - 1), sigma);

    double _Complex *Xa = malloc((size_t) (n * n) * sizeof(double _Complex));
    aic_latd_from_acb_mat(Xa, X);

    double sv[4];
    double _Complex U[16], Vt[16];
    aic_latd_svd(sv, U, Vt, Xa, n, n);

    /* smallest singular value must be the tiny sigma, not 0 or NaN. */
    double smin = sv[n - 1];
    AIC_CHECK_MSG(smin > 0.0 && smin == smin, /* >0 and not NaN */
                  "near-singular smallest singular value not graceful: %.3e",
                  smin);
    AIC_CHECK_MSG(fabs(smin - sigma) < 1e-12,
                  "smallest singular value %.3e != sigma %.3e", smin, sigma);

    /* set check vs arb@53 */
    arb_ptr sv_arb = _arb_vec_init(n);
    aic_mat_singular_values(sv_arb, X, 53);
    AIC_CHECK_EIGSET_CLOSE(sv, sv_arb, n, tol);

    /* reconstruction A = U diag(sv) Vt within tol (full-SVD sanity). */
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            double _Complex acc = 0;
            for (slong kk = 0; kk < n; kk++)
                acc += U[i * n + kk] * sv[kk] * Vt[kk * n + j];
            double want = (i == j)
                ? ((i == n - 1) ? sigma : diag_large[i]) : 0.0;
            AIC_CHECK_MSG(cabs(acc - want) < 1e-10,
                          "SVD reconstruction off at (%ld,%ld)",
                          (long) i, (long) j);
        }
    }

    _arb_vec_clear(sv_arb, n);
    free(Xa);
    acb_mat_clear(X);
    arb_clear(tol);
}

/* --- 4. random Hermitian double-vs-arb@53 (eigenvalues + opnorm) --- */
/* A well-separated random Hermitian so the arb SIMPLE-eig path applies; this is
 * the head-to-head double-vs-arb@53 agreement check on a generic input. */
static void test_random_hermitian(void)
{
    const slong n = 6;
    const slong prec = 53;
    acb_mat_t H;
    acb_mat_init(H, n, n);

    /* distinct dominant diagonal 10*(i+1) + small conjugate-symmetric
     * off-diagonals -> simple spectrum, fast arb isolation. */
    unsigned long seed = 0xC0FFEEUL;
    for (slong i = 0; i < n; i++) {
        acb_set_si(acb_mat_entry(H, i, i), 10 * (i + 1));
        for (slong j = i + 1; j < n; j++) {
            slong re = lcg_draw(&seed), im = lcg_draw(&seed);
            acb_set_si_si(acb_mat_entry(H, i, j), re, im);
            acb_set_si_si(acb_mat_entry(H, j, i), re, -im);
        }
    }

    double _Complex *Ha = malloc((size_t) (n * n) * sizeof(double _Complex));
    aic_latd_from_acb_mat(Ha, H);

    /* eigenvalues: double vs arb@53 set */
    double evals[6];
    aic_latd_eig_hermitian(evals, NULL, Ha, n);
    arb_ptr ev_arb = _arb_vec_init(n);
    aic_mat_eig_hermitian(ev_arb, NULL, H, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);
    AIC_CHECK_EIGSET_CLOSE(evals, ev_arb, n, tol);

    /* operator norm: double vs arb@53 */
    double op_d = aic_latd_opnorm(Ha, n, n);
    arb_t op_arb;
    arb_init(op_arb);
    aic_mat_opnorm(op_arb, H, prec);
    AIC_CHECK_DOUBLE_CLOSE(op_d, op_arb, tol);

    arb_clear(op_arb);
    arb_clear(tol);
    _arb_vec_clear(ev_arb, n);
    free(Ha);
    acb_mat_clear(H);
}

/* --- 5. conversion round-trip --- */
static void test_roundtrip(void)
{
    const slong r = 3, c = 4;
    acb_mat_t A, B;
    acb_mat_init(A, r, c);
    acb_mat_init(B, r, c);

    /* fill A with exactly-representable dyadic-rational midpoints so the
     * midpoint->double->zero-radius-acb round-trip is lossless. */
    unsigned long seed = 0x1234UL;
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_set_si_si(acb_mat_entry(A, i, j), lcg_draw(&seed),
                          lcg_draw(&seed));

    double _Complex *buf = malloc((size_t) (r * c) * sizeof(double _Complex));
    aic_latd_from_acb_mat(buf, A);
    aic_latd_to_acb_mat(B, buf, r, c);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 0.0); /* exact: integers are lossless through double */
    AIC_CHECK_ACB_MAT_CLOSE(A, B, tol);
    arb_clear(tol);

    free(buf);
    acb_mat_clear(B);
    acb_mat_clear(A);
}

/* --- 6. mutation proof: a wrong eigenvalue set must be REJECTED --- */
/* CLAUDE.md Rule 7 (mutation-prove): perturb a known-correct result and confirm
 * the cross-check rejects it (returns 0), proving the check has teeth. We call
 * the predicate aic_eigset_close directly (the macro aborts; here we want the
 * boolean) so we can assert it returns FALSE on corrupted input. */
static void test_mutation_proof(void)
{
    const slong n = 4;
    /* correct degenerate spectrum {-1,-1,1,1} */
    arb_ptr spec = _arb_vec_init(n);
    arb_set_si(spec + 0, -1);
    arb_set_si(spec + 1, -1);
    arb_set_si(spec + 2, 1);
    arb_set_si(spec + 3, 1);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);

    /* correct double set passes */
    double good[4] = {-1.0, -1.0, 1.0, 1.0};
    slong fk;
    AIC_CHECK_MSG(aic_eigset_close(good, spec, n, tol, &fk),
                  "mutation proof setup: correct set should pass");

    /* corrupted set (one eigenvalue moved by 1e-3 >> tol) must be REJECTED */
    double bad[4] = {-1.0, -1.0, 1.0, 1.001};
    AIC_CHECK_MSG(!aic_eigset_close(bad, spec, n, tol, &fk),
                  "mutation proof: corrupted set was NOT rejected (check has "
                  "no teeth)");

    arb_clear(tol);
    _arb_vec_clear(spec, n);
}

int main(void)
{
    test_degenerate_hermitian();
    test_nonnormal_svd();
    test_near_singular();
    test_random_hermitian();
    test_roundtrip();
    test_mutation_proof();

    aic_test_report("test_latd");
    printf("OK test_latd\n");
    return 0;
}
