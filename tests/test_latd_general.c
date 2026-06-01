/* test_latd_general.c — cross-checks for the general (non-Hermitian) eig path
 * aic_latd_eig_general + the spectral-gap helpers (bead aic-pvs). This is the
 * FAST/UNCERTIFIED double path (LAPACKE_zgeev); certified enclosures are the
 * separate arb concern (bead aic-w4o.1). CLAUDE.md Rule 5: every check asserts a
 * value/bound; Rule 6: the strongest check here is zgeev-vs-zheev on a Hermitian
 * input (an INDEPENDENT solver agreeing on the same spectrum).
 *
 * Checks:
 *  1. Closed-form spectra: diagonal -> evals == the diagonal (as a SET);
 *     upper-triangular -> evals == the diagonal; the 2x2 rotation [[0,1],[-1,0]]
 *     -> +/- i.
 *  2. THE cross-check (Rule 6): a HERMITIAN input's zgeev eigenvalues (which are
 *     real) match the independent aic_latd_eig_hermitian (zheev) eigenvalues,
 *     sorted, to ~1e-10. This validates zgeev's value path AND the row-major
 *     eigenvector layout (the eigenvector residual below shares the same call).
 *  3. Eigenvector residual ||A v_k - lambda_k v_k|| ~ 0 for each recovered pair,
 *     on a well-conditioned (Hermitian) input.
 *  4. Spectral-gap helpers: distance-to-z and pairwise-separation match hand
 *     values; a degenerate spectrum -> separation ~ 0 (degeneracy detected).
 *  5. Non-normal caveat (the tex:540 t^{1/3} companion): zgeev returns FINITE
 *     eigenvalues of rough magnitude t^{1/3}; we do NOT assert tight accuracy
 *     (uncertified ill-conditioned path).
 *  6. MUTATION PROOF (Rule 7): the zgeev-vs-zheev cross-check has teeth — a
 *     deliberately corrupted eigenvalue set is REJECTED by the set comparator.
 *
 * Convention reminder (include/aic_latd.h): row-major, arr[i*c+j] = entry (i,j);
 * evecs column k (evecs[i*n+k]) is the eigenvector for evals[k].
 */
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <flint/acb_mat.h>

#include "aic/aic_latd.h"
#include "aic_test.h"

/* --- complex eigenvalue-SET comparator (zgeev evals are complex doubles) ---
 * Greedy nearest-match: for each `got` value find an unused `want` within tol.
 * Returns 1 iff every value matches (an O(n^2) bijection check, n small). */
static int cset_close(const double _Complex *got, const double _Complex *want,
                      slong n, double tol)
{
    int *used = calloc((size_t) n, sizeof(int));
    if (!used) { fprintf(stderr, "cset_close: OOM\n"); abort(); }
    int ok = 1;
    for (slong i = 0; i < n && ok; i++) {
        int matched = 0;
        for (slong j = 0; j < n; j++) {
            if (!used[j] && cabs(got[i] - want[j]) <= tol) {
                used[j] = 1;
                matched = 1;
                break;
            }
        }
        if (!matched) ok = 0;
    }
    free(used);
    return ok;
}

/* --- 1. closed-form spectra --- */
static void test_closed_form(void)
{
    /* (a) diagonal: evals == the diagonal as a set. Off-diagonal kept 0. */
    {
        const slong n = 4;
        double _Complex A[16] = {0};
        double _Complex diag[4] = {2.0, -1.5, 3.0 + 1.0 * I, -0.25};
        for (slong i = 0; i < n; i++) A[i * n + i] = diag[i];
        double _Complex ev[4];
        aic_latd_eig_general(ev, NULL, A, n);
        AIC_CHECK_MSG(cset_close(ev, diag, n, 1e-12),
                      "diagonal eigenvalues do not match the diagonal");
    }
    /* (b) upper-triangular: eigenvalues are the diagonal (triangular spectrum is
     * the diagonal regardless of the strictly-upper entries). */
    {
        const slong n = 3;
        double _Complex A[9] = {0};
        double _Complex diag[3] = {1.0, 2.0, 3.0};
        for (slong i = 0; i < n; i++) A[i * n + i] = diag[i];
        A[0 * n + 1] = 5.0;   /* (0,1) */
        A[0 * n + 2] = -2.0;  /* (0,2) */
        A[1 * n + 2] = 7.0;   /* (1,2) */
        double _Complex ev[3];
        aic_latd_eig_general(ev, NULL, A, n);
        AIC_CHECK_MSG(cset_close(ev, diag, n, 1e-12),
                      "upper-triangular eigenvalues do not match the diagonal");
    }
    /* (c) 2x2 rotation [[0,1],[-1,0]] -> eigenvalues +/- i. */
    {
        const slong n = 2;
        double _Complex A[4] = {0.0, 1.0, -1.0, 0.0};
        double _Complex want[2] = {1.0 * I, -1.0 * I};
        double _Complex ev[2];
        aic_latd_eig_general(ev, NULL, A, n);
        AIC_CHECK_MSG(cset_close(ev, want, n, 1e-12),
                      "rotation [[0,1],[-1,0]] eigenvalues are not +/- i");
    }
}

/* Deterministic LCG -> small integer in [-4,4] (matches test_latd.c style). */
static slong lcg_draw(unsigned long *state)
{
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    return (slong) ((*state >> 41) % 9) - 4;
}

/* --- 2 & 3. zgeev-vs-zheev cross-check + eigenvector residual (Hermitian) --- */
static void test_hermitian_crosscheck(void)
{
    const slong n = 6;
    double _Complex H[36] = {0};
    /* distinct dominant diagonal + small conjugate-symmetric off-diagonals:
     * a well-conditioned Hermitian, real spectrum, well-separated. */
    unsigned long seed = 0xBEEFUL;
    for (slong i = 0; i < n; i++) {
        H[i * n + i] = (double) (10 * (i + 1));
        for (slong j = i + 1; j < n; j++) {
            double re = (double) lcg_draw(&seed);
            double im = (double) lcg_draw(&seed);
            H[i * n + j] = re + im * I;
            H[j * n + i] = re - im * I; /* Hermitian */
        }
    }

    /* zheev (independent path): real eigenvalues ascending. */
    double evh[6];
    aic_latd_eig_hermitian(evh, NULL, H, n);

    /* zgeev: complex eigenvalues + right eigenvectors. */
    double _Complex evg[6];
    double _Complex V[36];
    aic_latd_eig_general(evg, V, H, n);

    /* THE cross-check: zgeev's (real, up to ~1e-12 imaginary part) eigenvalues
     * equal zheev's as a set to ~1e-10. Build a complex `want` from the real
     * zheev set and compare with the complex comparator. */
    double _Complex evh_c[6];
    for (slong i = 0; i < n; i++) evh_c[i] = evh[i];
    AIC_CHECK_MSG(cset_close(evg, evh_c, n, 1e-10),
                  "zgeev eigenvalues disagree with zheev (independent paths)");
    /* the imaginary parts of a Hermitian matrix's eigenvalues must be ~0. */
    for (slong i = 0; i < n; i++)
        AIC_CHECK_MSG(fabs(cimag(evg[i])) < 1e-10,
                      "zgeev eigenvalue %ld of a Hermitian matrix is not real "
                      "(imag = %.3e)", (long) i, cimag(evg[i]));

    /* eigenvector residual: ||H v_k - lambda_k v_k|| ~ 0 for each column k. */
    for (slong k = 0; k < n; k++) {
        double res = 0.0;
        for (slong i = 0; i < n; i++) {
            double _Complex Hv = 0.0;
            for (slong j = 0; j < n; j++) Hv += H[i * n + j] * V[j * n + k];
            double _Complex r = Hv - evg[k] * V[i * n + k];
            res += creal(r) * creal(r) + cimag(r) * cimag(r);
        }
        res = sqrt(res);
        AIC_CHECK_MSG(res < 1e-10,
                      "eigenvector residual ||H v - lambda v|| = %.3e too "
                      "large at column %ld", res, (long) k);
    }
}

/* --- 4. spectral-gap helpers --- */
static void test_gap_helpers(void)
{
    /* spectrum {0, 1, 2, 5} (real); distance to z=2.5 is min |lam-2.5| = 0.5
     * (from lam=2); pairwise separation min |lam_i-lam_j| = 1 (the 0-1 / 1-2). */
    {
        const slong n = 4;
        double _Complex ev[4] = {0.0, 1.0, 2.0, 5.0};
        double gap = aic_latd_spectral_gap(ev, n, 2.5 + 0.0 * I);
        AIC_CHECK_MSG(fabs(gap - 0.5) < 1e-14,
                      "spectral_gap to 2.5 = %.6g, expected 0.5", gap);
        double sep = aic_latd_spectral_separation(ev, n);
        AIC_CHECK_MSG(fabs(sep - 1.0) < 1e-14,
                      "spectral_separation = %.6g, expected 1.0", sep);
        /* distance to a spectrum point itself is 0. */
        double on = aic_latd_spectral_gap(ev, n, 1.0 + 0.0 * I);
        AIC_CHECK_MSG(on < 1e-14,
                      "spectral_gap to an eigenvalue should be 0, got %.6g", on);
    }
    /* complex spectrum {+i, -i}: distance to z=0 (the sgn split point / imaginary
     * axis foot) is 1; separation is 2. */
    {
        const slong n = 2;
        double _Complex ev[2] = {1.0 * I, -1.0 * I};
        double gap0 = aic_latd_spectral_gap(ev, n, 0.0 + 0.0 * I);
        AIC_CHECK_MSG(fabs(gap0 - 1.0) < 1e-14,
                      "spectral_gap of {+/-i} to 0 = %.6g, expected 1", gap0);
        double sep = aic_latd_spectral_separation(ev, n);
        AIC_CHECK_MSG(fabs(sep - 2.0) < 1e-14,
                      "spectral_separation of {+/-i} = %.6g, expected 2", sep);
    }
    /* DEGENERATE: a repeated eigenvalue -> separation ~ 0 (degeneracy detected).
     * Use a non-diagonal matrix with eigenvalues {2,2,5}: J = [[2,1,0],[0,2,0],
     * [0,0,5]] (a 2x2 Jordan block at 2 + an isolated 5). */
    {
        const slong n = 3;
        double _Complex A[9] = {0};
        A[0 * n + 0] = 2.0;
        A[0 * n + 1] = 1.0;
        A[1 * n + 1] = 2.0;
        A[2 * n + 2] = 5.0;
        double _Complex ev[3];
        aic_latd_eig_general(ev, NULL, A, n);
        double sep = aic_latd_spectral_separation(ev, n);
        AIC_CHECK_MSG(sep < 1e-8,
                      "degenerate spectrum {2,2,5}: separation %.3e should be "
                      "~0 (degeneracy NOT detected)", sep);
        /* n<2 returns +inf (no pair). */
        double _Complex one = 3.0;
        AIC_CHECK_MSG(isinf(aic_latd_spectral_separation(&one, 1)),
                      "separation of a single eigenvalue should be +inf");
    }
}

/* --- 5. non-normal caveat: tex:540 companion, FINITE + rough magnitude only --- */
static void test_nonnormal_caveat(void)
{
    /* [[0,1,0],[0,0,1],[t,0,0]] has characteristic poly lambda^3 = t, so the
     * eigenvalues are the three cube roots of t, magnitude t^{1/3}. With
     * t=1e-6, |lambda| ~ 1e-2. The SPECTRUM is hypersensitive (1/3-power), so we
     * assert only: (i) all eigenvalues finite (no NaN/Inf), (ii) magnitudes in a
     * generous band around t^{1/3} — NOT tight accuracy (uncertified path). */
    const slong n = 3;
    const double t = 1e-6;
    double _Complex A[9] = {0};
    A[0 * n + 1] = 1.0;
    A[1 * n + 2] = 1.0;
    A[2 * n + 0] = t;
    double _Complex ev[3];
    aic_latd_eig_general(ev, NULL, A, n);

    double expected_mag = cbrt(t); /* 1e-2 */
    for (slong i = 0; i < n; i++) {
        AIC_CHECK_MSG(isfinite(creal(ev[i])) && isfinite(cimag(ev[i])),
                      "non-normal eigenvalue %ld not finite (re=%.3e im=%.3e)",
                      (long) i, creal(ev[i]), cimag(ev[i]));
        double mag = cabs(ev[i]);
        /* a generous half-decade band: this is the uncertified caveat, not a
         * tight check (tex:540 perturbation sensitivity). */
        AIC_CHECK_MSG(mag > 0.3 * expected_mag && mag < 3.0 * expected_mag,
                      "non-normal eigenvalue %ld magnitude %.3e outside rough "
                      "band around t^{1/3}=%.3e (uncertified path; not a tight "
                      "failure but flag it)", (long) i, mag, expected_mag);
    }
    /* product of eigenvalues = det = t (this IS stable — the magnitudes
     * multiply even though individual phases are sensitive). */
    double _Complex prod = ev[0] * ev[1] * ev[2];
    AIC_CHECK_MSG(cabs(prod - t) < 1e-9,
                  "product of eigenvalues %.3e+%.3ei != det = t = %.3e",
                  creal(prod), cimag(prod), t);
}

/* --- 6. mutation proof: the zgeev-vs-zheev set comparator has teeth --- */
static void test_mutation_proof(void)
{
    const slong n = 3;
    double _Complex got[3] = {1.0, 2.0, 3.0};
    double _Complex want_good[3] = {3.0, 1.0, 2.0}; /* same set, permuted */
    AIC_CHECK_MSG(cset_close(got, want_good, n, 1e-10),
                  "mutation setup: a permuted correct set should match");
    /* corrupt one value by 1e-3 >> tol: must be REJECTED. */
    double _Complex want_bad[3] = {3.0, 1.001, 2.0};
    AIC_CHECK_MSG(!cset_close(got, want_bad, n, 1e-10),
                  "mutation proof: corrupted set was NOT rejected (comparator "
                  "has no teeth)");
}

int main(void)
{
    test_closed_form();
    test_hermitian_crosscheck();
    test_gap_helpers();
    test_nonnormal_caveat();
    test_mutation_proof();

    aic_test_report("test_latd_general");
    printf("OK test_latd_general\n");
    return 0;
}
