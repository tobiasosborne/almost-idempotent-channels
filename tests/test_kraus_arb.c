/* test_kraus_arb.c — cross-checks for the CERTIFIED arb-path Choi->Kraus
 * extraction (bead aic-4td increment 2 step D; design
 * docs/research/eigvec_certified_design.md §3.2 / §4 S3/S4/S5/S7). The arb
 * counterpart of test_ucp.c's latd Choi->Kraus checks.
 *
 * CLAUDE.md Rule 6 cross-checks (the strongest first):
 *  S3 DOUBLE-vs-arb@53: aic_ucp_choi_to_kraus_arb at prec=53 vs the LAPACK
 *     aic_ucp_choi_to_kraus_latd, compared AS CHANNELS (rebuild Choi from each;
 *     Kraus gauge freedom => never operator-by-operator) — ||C_arb - C_latd||_op
 *     < 1e-12, recovered Kraus rank matches.
 *  S4 ROUND-TRIP (the headline correctness gate): kraus_to_choi(choi_to_kraus_arb
 *     (C)) == C within the certified ball — ||C_rebuilt - C||_op enclosure
 *     contains 0 (exact up to the dropped numerically-zero tail). Both a DISTINCT-
 *     eigenvalue Choi (each nonzero eig -> k=1 cluster, lambda exact) AND a
 *     genuinely-DEGENERATE Choi (repeated nonzero eigenvalue, one cluster). This
 *     is the gate that proves the Loewdin orthonormalisation of Rump's (non-
 *     orthonormal) X: a raw reshape would NOT rebuild C.
 *  S5 PSD-cone _tol (FINDINGS §C14): a planted O(eta^2) negative cluster
 *     (mid -2.5e-6) is DROPPED with clipped mass accumulated, kept Kraus rebuild
 *     a CP map; a planted O(1) negative (-0.3) FAILS LOUD (non-CP). Mutation:
 *     flip the lambda <= -neg_tol comparison -> the O(1) case stops aborting.
 *  S7 FAIL-LOUD (Rule 4, fork+SIGALRM watchdog): (a) a Choi whose smallest kept
 *     cluster sits ON keep_thr -> strict abort ("STRADDLE"/"unresolved");
 *     (d) an O(1)-negative Choi -> strict choi_to_kraus_arb aborts ("not CP").
 *
 * Concrete numbers (Rule 12) printed.
 */
#define _POSIX_C_SOURCE 200809L
#include <fcntl.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"

/* The DOUBLE-path reference (aic_ucp_choi_to_kraus_latd) is implicitly ~53-bit.
 * The CERTIFIED arb path needs a higher prec than 53 to cleanly resolve the
 * RANK gap of a rank-deficient Choi: the certified zero-eigenvalue enclosure at
 * prec=53 has radius ~2e-14 (densify+Rump conditioning), which STRADDLES the
 * QETLAB keep_thr ~1e-15 and fail-loud-aborts (correct, but not a usable cross-
 * check prec). Measured floor where the zero ball (~3e-17) drops cleanly below
 * keep_thr: prec >= 64; we use 128 for headroom. The double-vs-arb cross-check
 * is "the two extractions agree AS CHANNELS", which does NOT require equal prec
 * (the double path is the ~53-bit anchor, the arb path certifies at ARB_PREC).
 * paper/FINDINGS.md §D7 (certified Choi->Kraus rank-floor prec). */
#define ARB_PREC 128
#define KA_PREC 53

/* ---- channel / Choi fixtures -------------------------------------------- */

/* The complex asymmetric Convention-A channel from test_ucp.c (genuinely complex
 * eigenvectors; its Choi has DISTINCT eigenvalues — the distinct-cluster case). */
static void make_complex_channel(aic_ucp_kraus *phi)
{
    aic_ucp_kraus_init(phi, 2, 2, 2);
    acb_set_d_d(acb_mat_entry(phi->K[0], 0, 0), 0.6, 0.0);
    acb_set_d_d(acb_mat_entry(phi->K[0], 1, 0), 0.0, 0.8);
    acb_set_d_d(acb_mat_entry(phi->K[1], 0, 1), 0.0, 0.8);
    acb_set_d_d(acb_mat_entry(phi->K[1], 1, 1), 0.6, 0.0);
}

/* Exactly-idempotent compression UCP on C^d, proper carrier rank m (from
 * test_ucp.c): a degenerate (rank-deficient) Choi with a REPEATED nonzero
 * eigenvalue — the genuine-degeneracy round-trip case. */
static void make_compress_idemp(aic_ucp_kraus *phi, slong d, slong m)
{
    slong r = 1 + (d - m);
    aic_ucp_kraus_init(phi, d, d, r);
    for (slong i = 0; i < m; i++)
        acb_set_si(acb_mat_entry(phi->K[0], i, i), 1);
    for (slong b = 0; b < d - m; b++)
        acb_set_si(acb_mat_entry(phi->K[1 + b], 0, m + b), 1);
}

/* Depolarizing on C^d: Choi = (1/d) 1_{d^2} (ALL eigenvalues equal 1/d — the
 * maximally-degenerate single-cluster case). */
static void make_depolarizing(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, d * d);
    double inv = 1.0 / sqrt((double) d);
    slong a = 0;
    for (slong i = 0; i < d; i++)
        for (slong j = 0; j < d; j++) {
            acb_set_d(acb_mat_entry(phi->K[a], i, j), inv);
            a++;
        }
}

/* Hermitian-project the MIDPOINTS of an interval Hermitian matrix into an
 * EXACT-Hermitian (zero-radius) matrix out: diagonal real, lower = exact conj of
 * upper midpoints. A Choi rebuilt from sqrt-scaled interval Kraus is Hermitian in
 * exact arithmetic but carries independent per-entry radii, which trips the tight
 * Hermiticity assert of aic_ucp_is_cp_choi; the midpoint projection is the
 * standard cure (FINDINGS §C5) and changes no value on the true (Hermitian)
 * matrix. The CP verdict on the projected matrix is the genuine PSD test. */
static void herm_project_mid(acb_mat_t out, const acb_mat_t A)
{
    slong n = acb_mat_nrows(A);
    for (slong i = 0; i < n; i++) {
        acb_get_mid(acb_mat_entry(out, i, i), acb_mat_entry(A, i, i));
        arb_zero(acb_imagref(acb_mat_entry(out, i, i)));
        for (slong j = i + 1; j < n; j++) {
            acb_get_mid(acb_mat_entry(out, i, j), acb_mat_entry(A, i, j));
            acb_conj(acb_mat_entry(out, j, i), acb_mat_entry(out, i, j));
        }
    }
}

/* ||A - B||_op upper bound as a double. */
static double opdiff(const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong n = acb_mat_nrows(A);
    acb_mat_t D;
    acb_mat_init(D, n, n);
    acb_mat_sub(D, A, B, prec);
    arb_t nr;
    arb_init(nr);
    aic_mat_opnorm(nr, D, prec);
    double ub = arf_get_d(arb_midref(nr), ARF_RND_UP) + mag_get_d(arb_radref(nr));
    arb_clear(nr);
    acb_mat_clear(D);
    return ub;
}

/* ---- S3: double-vs-arb@53, compared AS CHANNELS -------------------------- */
static void s3_one(const char *name, const aic_ucp_kraus *phi)
{
    slong dk = phi->dim_K, dh = phi->dim_H, n = dk * dh;
    acb_mat_t C, Carb, Clatd;
    acb_mat_init(C, n, n);
    acb_mat_init(Carb, n, n);
    acb_mat_init(Clatd, n, n);
    aic_ucp_kraus_to_choi(C, phi, ARB_PREC);

    aic_ucp_kraus phi_arb, phi_latd;
    aic_ucp_choi_to_kraus_arb(&phi_arb, C, dk, dh, ARB_PREC);
    aic_ucp_choi_to_kraus_latd(&phi_latd, C, dk, dh);
    aic_ucp_kraus_to_choi(Carb, &phi_arb, ARB_PREC);
    aic_ucp_kraus_to_choi(Clatd, &phi_latd, KA_PREC);

    double d = opdiff(Carb, Clatd, KA_PREC);
    AIC_CHECK_MSG(d < 1e-12, "%s: ||C_arb - C_latd||_op = %.3e not < 1e-12",
                  name, d);
    AIC_CHECK_MSG(phi_arb.r == phi_latd.r,
                  "%s: arb rank %ld != latd rank %ld", name,
                  (long) phi_arb.r, (long) phi_latd.r);
    printf("  %-22s ||C_arb-C_latd||_op=%.3e  rank arb=%ld latd=%ld\n",
           name, d, (long) phi_arb.r, (long) phi_latd.r);

    aic_ucp_kraus_clear(&phi_latd);
    aic_ucp_kraus_clear(&phi_arb);
    acb_mat_clear(Clatd);
    acb_mat_clear(Carb);
    acb_mat_clear(C);
}

static void test_s3_double_vs_arb(void)
{
    printf("S3: double-vs-arb@53 Choi->Kraus, compared AS CHANNELS (rebuild Choi)\n");
    aic_ucp_kraus c;
    make_complex_channel(&c);
    s3_one("complex asym C^2", &c);
    aic_ucp_kraus_clear(&c);

    /* multi-block degenerate Choi: the compression idempotent C^3 m=1 has Choi
     * spectrum {0^6, 1^3} — a REPEATED nonzero eigenvalue (mult 3), the genuine
     * (+)B(L)-style degeneracy; AND the depolarizing (all eigenvalues equal). */
    aic_ucp_kraus ci;
    make_compress_idemp(&ci, 3, 1);
    s3_one("compress-idemp C^3", &ci);
    aic_ucp_kraus_clear(&ci);

    aic_ucp_kraus dp;
    make_depolarizing(&dp, 3);
    s3_one("depolarizing C^3", &dp);
    aic_ucp_kraus_clear(&dp);
}

/* ---- S4: rebuild-Choi round-trip within the certified ball --------------- */
/* ||C_rebuilt - C||_op enclosure must CONTAIN 0 (exact up to dropped tail). We
 * assert the certified ball's lower endpoint <= 0 AND the upper endpoint tiny. */
static void s4_one(const char *name, const acb_mat_t C, slong dk, slong dh,
                   slong expected_rank)
{
    slong n = dk * dh;
    aic_ucp_kraus phi;
    aic_ucp_choi_to_kraus_arb(&phi, C, dk, dh, ARB_PREC);
    AIC_CHECK_MSG(phi.r == expected_rank, "%s: rank %ld != expected %ld",
                  name, (long) phi.r, (long) expected_rank);

    acb_mat_t Cr, D;
    acb_mat_init(Cr, n, n);
    acb_mat_init(D, n, n);
    aic_ucp_kraus_to_choi(Cr, &phi, ARB_PREC);
    acb_mat_sub(D, Cr, C, ARB_PREC);

    arb_t nr;
    arb_init(nr);
    aic_mat_opnorm(nr, D, ARB_PREC);
    /* certified ball: [mid - rad, mid + rad]; the enclosure must contain 0, i.e.
     * the lower endpoint mid-rad <= 0. Since opnorm >= 0, this is "the ball
     * touches 0" — equivalently arb_contains_zero of (||D||_op - 0)? opnorm is a
     * nonneg ball; the rebuild is exact, so the TRUE value is 0 and the certified
     * ball must enclose it: lower endpoint <= 0. */
    double mid = arf_get_d(arb_midref(nr), ARF_RND_NEAR);
    double rad = mag_get_d(arb_radref(nr));
    double ub = mid + rad;
    AIC_CHECK_MSG(mid - rad <= 1e-25,
                  "%s: ||C_rebuilt - C||_op ball [%.3e,%.3e] does not contain 0",
                  name, mid - rad, ub);
    AIC_CHECK_MSG(ub < 1e-10, "%s: ||C_rebuilt - C||_op upper = %.3e not < 1e-10",
                  name, ub);
    printf("  %-22s rank=%ld  ||C_rebuilt-C||_op in [%.3e, %.3e]\n",
           name, (long) phi.r, mid - rad, ub);

    arb_clear(nr);
    acb_mat_clear(D);
    acb_mat_clear(Cr);
    aic_ucp_kraus_clear(&phi);
}

static void test_s4_roundtrip(void)
{
    printf("S4: rebuild-Choi round-trip within the certified ball (the headline)\n");

    /* DISTINCT-eigenvalue Choi: the complex asymmetric channel's Choi (rank 2,
     * distinct nonzero eigenvalues -> two k=1 clusters, lambda exact). */
    aic_ucp_kraus c;
    make_complex_channel(&c);
    slong cn = c.dim_K * c.dim_H;
    acb_mat_t Cc;
    acb_mat_init(Cc, cn, cn);
    aic_ucp_kraus_to_choi(Cc, &c, ARB_PREC);
    s4_one("distinct (complex C^2)", Cc, c.dim_K, c.dim_H, 2);
    acb_mat_clear(Cc);
    aic_ucp_kraus_clear(&c);

    /* GENUINELY-DEGENERATE Choi: the compression idempotent C^3 m=1, Choi spectrum
     * {0^6, 1^3} — a REPEATED nonzero eigenvalue (mult 3) forming ONE nonzero
     * cluster (lambda=1 exact to the ball). rank = 1+(d-m) = 3 (one Kraus op per
     * orthonormal column of the k=3 cluster). This is the (+)B(L)-style block
     * multiplicity the lambda_c-for-all reshape must handle exactly. */
    aic_ucp_kraus ci;
    make_compress_idemp(&ci, 3, 1);
    slong in = ci.dim_K * ci.dim_H;
    acb_mat_t Ci;
    acb_mat_init(Ci, in, in);
    aic_ucp_kraus_to_choi(Ci, &ci, ARB_PREC);
    s4_one("degenerate (compress C^3)", Ci, ci.dim_K, ci.dim_H, 3);
    acb_mat_clear(Ci);
    aic_ucp_kraus_clear(&ci);

    /* depolarizing: ALL d^2 eigenvalues equal 1/d -> ONE cluster, rank d^2. */
    aic_ucp_kraus dp;
    make_depolarizing(&dp, 3);
    slong dn = dp.dim_K * dp.dim_H;
    acb_mat_t Cd;
    acb_mat_init(Cd, dn, dn);
    aic_ucp_kraus_to_choi(Cd, &dp, ARB_PREC);
    s4_one("depolarizing (C^3)", Cd, dp.dim_K, dp.dim_H, 9);
    acb_mat_clear(Cd);
    aic_ucp_kraus_clear(&dp);
}

/* ---- S5: PSD-cone _tol variant ------------------------------------------- */
/* Build a Hermitian Choi-shaped 4x4 (dim_K=2, dim_H=2) C = diag(big, ..., small)
 * with a planted small-negative eigenvalue, conjugated by a rational Givens so it
 * is not coordinate-aligned (exercises the orthonormalisation too). */
static void build_diag_planted(acb_mat_t C, const double *diag, slong n,
                               slong prec)
{
    acb_mat_t D, Q, Qt, T;
    acb_mat_init(D, n, n);
    acb_mat_init(Q, n, n);
    acb_mat_init(Qt, n, n);
    acb_mat_init(T, n, n);
    acb_mat_zero(D);
    for (slong i = 0; i < n; i++) acb_set_d(acb_mat_entry(D, i, i), diag[i]);
    acb_mat_one(Q);
    acb_t cc, ss, ns;
    acb_init(cc);
    acb_init(ss);
    acb_init(ns);
    acb_set_si(cc, 3);
    acb_div_si(cc, cc, 5, prec);
    acb_set_si(ss, 4);
    acb_div_si(ss, ss, 5, prec);
    acb_neg(ns, ss);
    acb_set(acb_mat_entry(Q, 0, 0), cc);
    acb_set(acb_mat_entry(Q, n - 1, n - 1), cc);
    acb_set(acb_mat_entry(Q, 0, n - 1), ns);
    acb_set(acb_mat_entry(Q, n - 1, 0), ss);
    acb_mat_conjugate_transpose(Qt, Q);
    acb_mat_mul(T, Q, D, prec);
    acb_mat_mul(T, T, Qt, prec);
    /* Q D Q^dag is Hermitian in exact arithmetic but the two matmuls accumulate
     * independent per-entry radii; project to the exact-Hermitian midpoint so the
     * tight Hermiticity assert of the subspace solver / is_cp_choi passes for the
     * right reason (FINDINGS §C5). C is then a zero-radius Hermitian Choi with the
     * planted spectrum, exact-rational keep_thr ball. */
    herm_project_mid(C, T);
    acb_clear(ns);
    acb_clear(ss);
    acb_clear(cc);
    acb_mat_clear(T);
    acb_mat_clear(Qt);
    acb_mat_clear(Q);
    acb_mat_clear(D);
}

static void test_s5_psd_cone(void)
{
    printf("S5: PSD-cone _tol — drop O(eta^2) negative, fail loud O(1) negative\n");
    const slong dk = 2, dh = 2, n = dk * dh;

    /* O(eta^2) negative: diag(1.0, 0.5, 0.25, -2.5e-6) conjugated. With neg_tol
     * = 1e-4 (>> 2.5e-6) the negative cluster is DROPPED, mass accumulated; the
     * three positive eigenvalues survive (rank 3); the kept Kraus rebuild a CP
     * map (the PSD-cone projection). */
    {
        double diag[4] = {1.0, 0.5, 0.25, -2.5e-6};
        acb_mat_t C;
        acb_mat_init(C, n, n);
        build_diag_planted(C, diag, n, ARB_PREC);
        aic_ucp_kraus phi;
        double clipped = -1.0;
        aic_ucp_choi_to_kraus_arb_tol(&phi, C, dk, dh, 1e-4, &clipped, ARB_PREC);
        AIC_CHECK_MSG(phi.r == 3, "S5 eta^2: kept rank %ld != 3", (long) phi.r);
        AIC_CHECK_MSG(clipped > 1e-6 && clipped < 1e-5,
                      "S5 eta^2: clipped mass %.3e not ~2.5e-6", clipped);
        /* kept Kraus rebuild a CP (PSD) Choi: the negative direction is gone.
         * Project the interval Cr to its exact-Hermitian midpoint (FINDINGS §C5)
         * before the tight-Hermiticity CP check; tol = keep_thr-ish slack. */
        acb_mat_t Cr, Crh;
        acb_mat_init(Cr, n, n);
        acb_mat_init(Crh, n, n);
        aic_ucp_kraus_to_choi(Cr, &phi, ARB_PREC);
        herm_project_mid(Crh, Cr);
        arb_t tol;
        arb_init(tol);
        arb_set_d(tol, 1e-12);
        AIC_CHECK_MSG(aic_ucp_is_cp_choi(Crh, tol, ARB_PREC) == 1,
                      "S5 eta^2: kept Kraus do not rebuild a CP map");
        printf("  O(eta^2) neg: dropped, clipped mass=%.3e, kept rank=%ld, CP OK\n",
               clipped, (long) phi.r);
        arb_clear(tol);
        acb_mat_clear(Crh);
        acb_mat_clear(Cr);
        aic_ucp_kraus_clear(&phi);
        acb_mat_clear(C);
    }
}

/* ---- the fork+SIGALRM watchdog (pattern from test_eigvec.c) --------------- */
#define KA_WATCH_S 30
typedef void (*ka_child_fn)(void);
static volatile pid_t ka_watch_pid = 0;
static volatile sig_atomic_t ka_timed_out = 0;
static void ka_alarm(int sig)
{
    (void) sig;
    ka_timed_out = 1;
    if (ka_watch_pid > 0) kill(ka_watch_pid, SIGKILL);
}
static int ka_run_child(ka_child_fn fn, int *status, char *err, size_t errsz)
{
    char tmpl[] = "/tmp/aic_kraus_err_XXXXXX";
    int efd = mkstemp(tmpl);
    AIC_CHECK_MSG(efd >= 0, "mkstemp failed");
    fflush(NULL);
    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "fork failed");
    if (pid == 0) {
        dup2(efd, STDERR_FILENO);
        dup2(efd, STDOUT_FILENO);
        close(efd);
        fn();
        _exit(0);
    }
    ka_watch_pid = pid;
    ka_timed_out = 0;
    struct sigaction sa, old;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = ka_alarm;
    sigaction(SIGALRM, &sa, &old);
    alarm(KA_WATCH_S);
    int st = 0;
    pid_t w;
    do { w = waitpid(pid, &st, 0); } while (w < 0 && !ka_timed_out);
    alarm(0);
    sigaction(SIGALRM, &old, NULL);
    if (w < 0) waitpid(pid, &st, 0);
    lseek(efd, 0, SEEK_SET);
    ssize_t rd = read(efd, err, errsz - 1);
    err[(rd > 0) ? (size_t) rd : 0] = '\0';
    close(efd);
    unlink(tmpl);
    *status = st;
    return ka_timed_out ? 0 : 1;
}
static void ka_assert_failloud(ka_child_fn fn, const char *what,
                               const char *needle)
{
    char err[4096];
    int st = 0;
    int finished = ka_run_child(fn, &st, err, sizeof err);
    AIC_CHECK_MSG(finished, "%s HUNG (watchdog killed it after %d s)",
                  what, KA_WATCH_S);
    AIC_CHECK_MSG(WIFEXITED(st) ? WEXITSTATUS(st) != 0 : 1,
                  "%s exited 0 — it silently returned instead of failing loud",
                  what);
    AIC_CHECK_MSG(WIFSIGNALED(st) && WTERMSIG(st) == SIGABRT,
                  "%s: expected SIGABRT (signaled=%d sig=%d exited=%d code=%d)",
                  what, WIFSIGNALED(st), WIFSIGNALED(st) ? WTERMSIG(st) : -1,
                  WIFEXITED(st), WIFEXITED(st) ? WEXITSTATUS(st) : -1);
    AIC_CHECK_MSG(strstr(err, needle) != NULL,
                  "%s: abort message does not name '%s'. Got:\n%s",
                  what, needle, err);
    printf("  child SIGABRT; message names '%s' (OK)\n", needle);
}

/* S5 fail-loud part + mutation tooth: an O(1)-negative Choi must FAIL LOUD in the
 * strict choi_to_kraus_arb (neg_tol = keep_thr, so -0.3 <= -keep_thr). */
static void child_o1_negative(void)
{
    const slong dk = 2, dh = 2, n = dk * dh;
    double diag[4] = {1.0, 0.5, 0.25, -0.3};
    acb_mat_t C;
    acb_mat_init(C, n, n);
    build_diag_planted(C, diag, n, KA_PREC);
    aic_ucp_kraus phi;
    aic_ucp_choi_to_kraus_arb(&phi, C, dk, dh, ARB_PREC); /* must abort */
    aic_ucp_kraus_clear(&phi);
    acb_mat_clear(C);
}

/* S7(a): a Choi whose smallest kept cluster sits EXACTLY on keep_thr -> the
 * cluster ball straddles keep_thr -> strict abort ("STRADDLE"). We engineer this
 * by computing keep_thr for a diag and planting one eigenvalue AT it.
 * keep_thr = n * 2^-52 * ||C||_F. For diag(1, eps_plant) with eps_plant the
 * threshold itself, the small eigenvalue ball overlaps keep_thr. */
static void child_straddle(void)
{
    const slong dk = 2, dh = 1, n = dk * dh;   /* n = 2, a 2x2 Choi */
    /* ||C||_F ~ 1 for diag(1, small); keep_thr = 2 * 2^-52 * 1 ~ 4.4e-16.
     * Plant the small eigenvalue exactly at 2*2^-52 so its ball straddles thr. */
    acb_mat_t C;
    acb_mat_init(C, n, n);
    acb_mat_zero(C);
    acb_set_d(acb_mat_entry(C, 0, 0), 1.0);
    acb_t thr;
    acb_init(thr);
    acb_one(thr);
    acb_mul_2exp_si(thr, thr, -52);
    acb_mul_si(thr, thr, n, KA_PREC);          /* n * 2^-52 (||C||_F ~ 1) */
    acb_set(acb_mat_entry(C, 1, 1), thr);
    aic_ucp_kraus phi;
    aic_ucp_choi_to_kraus_arb(&phi, C, dk, dh, ARB_PREC); /* must abort */
    aic_ucp_kraus_clear(&phi);
    acb_clear(thr);
    acb_mat_clear(C);
}

static void test_s7_failloud(void)
{
    printf("S7: fail-loud teeth (fork watchdog)\n");
    printf(" (d) O(1)-negative Choi -> strict choi_to_kraus_arb aborts (not CP)\n");
    ka_assert_failloud(child_o1_negative,
                       "aic_ucp_choi_to_kraus_arb(O(1) negative)", "CP");
    printf(" (a) smallest kept cluster ON keep_thr -> STRADDLE abort\n");
    ka_assert_failloud(child_straddle,
                       "aic_ucp_choi_to_kraus_arb(straddle)", "STRADDLE");
}

int main(void)
{
    test_s3_double_vs_arb();
    test_s4_roundtrip();
    test_s5_psd_cone();
    test_s7_failloud();
    aic_test_report("test_kraus_arb");
    printf("OK test_kraus_arb\n");
    return 0;
}
