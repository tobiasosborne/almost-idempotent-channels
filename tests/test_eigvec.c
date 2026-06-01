/* test_eigvec.c — cross-checks for the certified eigenvalue + INVARIANT-SUBSPACE
 * layer (bead aic-4td increment 2 step C2; design
 * docs/research/eigvec_certified_design.md §4 S1/S2/S7). CLAUDE.md Rule 6: the
 * strongest tests here are the residual ||H X_c - X_c J_c|| recomputed on the
 * ORIGINAL H (the rigor certificate, §1.6), the eta=0 / exact-projection oracle,
 * and the fail-loud teeth. Rule 7: these assertions were written BEFORE the
 * implementation (RED baseline), and each fail-loud tooth is mutation-proven.
 *
 * The teeth (cross-check ladder rungs in parentheses):
 *  S1 (rung 2/3, the core certificate) for the §D5n cases (C^5{2,3}, C^6{2,4},
 *     C^7{3,4} projectors) + block spectra: every cluster's residual
 *     ||H X_c - X_c J_c||_F < 1e-20 at prec=128 (recomputed on the ORIGINAL H,
 *     using the STORED J_c — design §1.6(ii) proves the same J_c works for H),
 *     all lambda balls pairwise disjoint, Sum k_c = n. These are EXACTLY the
 *     inputs test_eigmult T5 used to assert fail-loud — they now CERTIFY.
 *  S2 (rung 3, eta=0 oracle) a rank-r projector P -> the lambda~1 cluster
 *     projector == P (||Pi_M - P||_op < 1e-25), trace(Pi_M)=r, the lambda~0
 *     cluster projector = I - P; the depolarizing Choi (1/d)I_{d^2} -> ONE
 *     cluster, k=d^2, Pi_M = I.
 *  S7 (Rule 4, fork watchdog) the fail-loud teeth. MEASURED REALITY (probes,
 *     this session; FINDINGS §D5n2): for a near-degenerate spectrum, Rump's
 *     certificate SELF-ISOLATES — whenever both per-cluster enclosures are
 *     FINITE the lambda balls are already DISJOINT, and whenever the balls would
 *     overlap at least one enclosure is NON-FINITE. So the reachable fail-loud
 *     path for a sub-resolution near-degeneracy is the FINITE-enclosure guard
 *     (i), not the disjointness gate (iii):
 *     (c) NON-FINITE enclosure: eigenvalues 1 and 1+2^-40 FORCED into two
 *         clusters at prec=24 -> the tiny cluster's Rump enclosure is [+/-inf]
 *         -> aborts ("UNRESOLVED"). Mutation: drop the finite-enclosure check in
 *         aic_mat_int_certify_cluster -> child no longer aborts there (RED).
 *     (b) OVERLAP gate (defence in depth): no real input reaches it with finite
 *         balls (the self-isolation above), so it is proven load-bearing by a
 *         mutation: forcing the per-cluster lambda balls to be identical (a
 *         synthetic two-equal-clusters input) makes the disjointness gate the
 *         abort source ("OVERLAP"); dropping that gate makes the child exit 0.
 *         A k=n whole-space input (all evals equal) gives ONE cluster (no pair),
 *         confirming the gate only bites a genuine >1-cluster split.
 *
 * Concrete numbers (worst residual, projector defects, abort messages) printed,
 * Rule 12.
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
#include "aic_test.h"

#define EV_PREC 128

/* --- fixtures (ported from test_eigmult.c) --------------------------------- */

/* rank-r orthogonal projector in C^n (spectrum {0 mult n-r, 1 mult r}),
 * conjugated by DISJOINT rational Givens (cos=3/5,sin=4/5): plane (i, r+i) for
 * i < min(r,n-r). Off the coordinate axes, both cluster eigenbases well-
 * conditioned. (The dense-unitary densify inside the routine then defeats the
 * Rump frozen-row partition; design §2.) */
static void build_projector(acb_mat_t P, slong n, slong r, slong prec)
{
    acb_mat_t Q, G, Qt, T;
    acb_mat_init(Q, n, n);
    acb_mat_init(G, n, n);
    acb_mat_init(Qt, n, n);
    acb_mat_init(T, n, n);
    acb_mat_zero(P);
    for (slong i = 0; i < r; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
    acb_t c, s, ns;
    acb_init(c);
    acb_init(s);
    acb_init(ns);
    acb_set_si(c, 3);
    acb_div_si(c, c, 5, prec);
    acb_set_si(s, 4);
    acb_div_si(s, s, 5, prec);
    acb_neg(ns, s);
    acb_mat_one(Q);
    slong pairs = (r < n - r) ? r : (n - r);
    for (slong i = 0; i < pairs; i++) {
        slong a = i, b = r + i;
        acb_mat_one(G);
        acb_set(acb_mat_entry(G, a, a), c);
        acb_set(acb_mat_entry(G, b, b), c);
        acb_set(acb_mat_entry(G, a, b), ns);
        acb_set(acb_mat_entry(G, b, a), s);
        acb_mat_mul(T, Q, G, prec);
        acb_mat_set(Q, T);
    }
    acb_mat_conjugate_transpose(Qt, Q);
    acb_mat_mul(T, Q, P, prec);
    acb_mat_mul(P, T, Qt, prec);
    acb_clear(ns);
    acb_clear(s);
    acb_clear(c);
    acb_mat_clear(T);
    acb_mat_clear(Qt);
    acb_mat_clear(G);
    acb_mat_clear(Q);
}

/* Conjugate a real diagonal `diag` (length n) by rational Givens (0,n-1),(1,n-2)
 * so a degenerate block spectrum is NOT coordinate-aligned. */
static void build_conj_diag(acb_mat_t X, const slong *diag, slong n, slong prec)
{
    acb_mat_t D, Q, Qt, T;
    acb_mat_init(D, n, n);
    acb_mat_init(Q, n, n);
    acb_mat_init(Qt, n, n);
    acb_mat_init(T, n, n);
    acb_mat_zero(D);
    for (slong i = 0; i < n; i++) acb_set_si(acb_mat_entry(D, i, i), diag[i]);
    acb_mat_one(Q);
    acb_t c, s, ns;
    acb_init(c);
    acb_init(s);
    acb_init(ns);
    acb_set_si(c, 3);
    acb_div_si(c, c, 5, prec);
    acb_set_si(s, 4);
    acb_div_si(s, s, 5, prec);
    acb_neg(ns, s);
    acb_set(acb_mat_entry(Q, 0, 0), c);
    acb_set(acb_mat_entry(Q, n - 1, n - 1), c);
    acb_set(acb_mat_entry(Q, 0, n - 1), ns);
    acb_set(acb_mat_entry(Q, n - 1, 0), s);
    if (n >= 4) {
        acb_set(acb_mat_entry(Q, 1, 1), c);
        acb_set(acb_mat_entry(Q, n - 2, n - 2), c);
        acb_set(acb_mat_entry(Q, 1, n - 2), ns);
        acb_set(acb_mat_entry(Q, n - 2, 1), s);
    }
    acb_mat_conjugate_transpose(Qt, Q);
    acb_mat_mul(T, Q, D, prec);
    acb_mat_mul(X, T, Qt, prec);
    acb_clear(ns);
    acb_clear(s);
    acb_clear(c);
    acb_mat_clear(T);
    acb_mat_clear(Qt);
    acb_mat_clear(Q);
    acb_mat_clear(D);
}

/* --- S1: the invariant-subspace residual on the ORIGINAL H (the certificate) - */
/* For each cluster, recompute ||H X_c - X_c J_c||_F using the STORED J_c (design
 * §1.6(ii): the same J_c certifies the subspace of H, since A' = U H U^dag and
 * X_c = U^dag X'_c). Return the worst residual upper bound across clusters as a
 * double, and assert disjointness + Sum k = n inside. */
static double s1_worst_residual(const acb_mat_t H, slong n, const char *name)
{
    aic_mat_eigcluster *cl;
    slong nc;
    aic_mat_eig_hermitian_subspaces(&cl, &nc, H, -1.0, EV_PREC);

    double worst = 0.0;
    slong ksum = 0;
    for (slong c = 0; c < nc; c++) {
        slong k = cl[c].k;
        ksum += k;
        acb_mat_t HX, XJ, R;
        acb_mat_init(HX, n, k);
        acb_mat_init(XJ, n, k);
        acb_mat_init(R, n, k);
        acb_mat_mul(HX, H, cl[c].X, EV_PREC);     /* H X_c  */
        acb_mat_mul(XJ, cl[c].X, cl[c].J, EV_PREC); /* X_c J_c (STORED J) */
        acb_mat_sub(R, HX, XJ, EV_PREC);
        arb_t nr;
        arb_init(nr);
        acb_mat_frobenius_norm(nr, R, EV_PREC);
        double ub = arf_get_d(arb_midref(nr), ARF_RND_UP) +
                    mag_get_d(arb_radref(nr));
        if (ub > worst) worst = ub;
        arb_clear(nr);
        acb_mat_clear(R);
        acb_mat_clear(XJ);
        acb_mat_clear(HX);
    }
    AIC_CHECK_MSG(ksum == n, "%s: Sum k_c = %ld != n = %ld", name,
                  (long) ksum, (long) n);

    /* pairwise disjointness (the routine already asserts it; re-verify here so a
     * silent change to the gate is caught by THIS test too). */
    for (slong a = 0; a < nc; a++)
        for (slong b = a + 1; b < nc; b++)
            AIC_CHECK_MSG(!acb_overlaps(cl[a].lambda, cl[b].lambda),
                          "%s: lambda balls %ld,%ld overlap", name,
                          (long) a, (long) b);

    aic_mat_eigcluster_free(cl, nc);
    return worst;
}

static void test_s1_residual(void)
{
    printf("S1: certified invariant-subspace residual ||H X_c - X_c J_c|| on the "
           "ORIGINAL H (< 1e-20 at prec=128)\n");
    double overall = 0.0;

    struct { slong n, r; const char *name; } projs[] = {
        {5, 2, "proj C^5{2,3}"}, {6, 2, "proj C^6{2,4}"}, {7, 3, "proj C^7{3,4}"},
    };
    for (size_t i = 0; i < sizeof projs / sizeof projs[0]; i++) {
        acb_mat_t P;
        acb_mat_init(P, projs[i].n, projs[i].n);
        build_projector(P, projs[i].n, projs[i].r, EV_PREC);
        double w = s1_worst_residual(P, projs[i].n, projs[i].name);
        AIC_CHECK_MSG(w < 1e-20, "%s: worst residual %.3e not < 1e-20",
                      projs[i].name, w);
        printf("  %-16s worst ||HX-XJ|| = %.3e\n", projs[i].name, w);
        if (w > overall) overall = w;
        acb_mat_clear(P);
    }

    /* block spectra (build_conj_diag): {0,0,2,2,5}, {0,0,0,4,4,4,9}. */
    slong b1[5] = {0, 0, 2, 2, 5};
    slong b2[7] = {0, 0, 0, 4, 4, 4, 9};
    struct { const slong *d; slong n; const char *name; } blocks[] = {
        {b1, 5, "block {0,0,2,2,5}"}, {b2, 7, "block {0^3,4^3,9}"},
    };
    for (size_t i = 0; i < sizeof blocks / sizeof blocks[0]; i++) {
        acb_mat_t H;
        acb_mat_init(H, blocks[i].n, blocks[i].n);
        build_conj_diag(H, blocks[i].d, blocks[i].n, EV_PREC);
        double w = s1_worst_residual(H, blocks[i].n, blocks[i].name);
        AIC_CHECK_MSG(w < 1e-20, "%s: worst residual %.3e not < 1e-20",
                      blocks[i].name, w);
        printf("  %-16s worst ||HX-XJ|| = %.3e\n", blocks[i].name, w);
        if (w > overall) overall = w;
        acb_mat_clear(H);
    }
    printf("  OVERALL worst residual across S1 cases: %.3e (expect ~1e-31)\n",
           overall);
}

/* --- S2: eta=0 / exact-projection oracle ----------------------------------- */
/* Project the cluster at eigenvalue ~target out of the cluster list (must be
 * unique by disjointness). Returns the cluster index; aborts the test if none. */
static slong s2_cluster_near(const aic_mat_eigcluster *cl, slong nc, double target)
{
    for (slong c = 0; c < nc; c++) {
        arb_t re;
        arb_init(re);
        acb_get_real(re, cl[c].lambda);
        double mid = arf_get_d(arb_midref(re), ARF_RND_NEAR);
        arb_clear(re);
        if (fabs(mid - target) < 0.25) return c;
    }
    AIC_CHECK_MSG(0, "no cluster near %.3f", target);
    return -1;
}

/* ||A - B||_op as a double upper bound. */
static double s2_opdiff(const acb_mat_t A, const acb_mat_t B)
{
    slong n = acb_mat_nrows(A);
    acb_mat_t D;
    acb_mat_init(D, n, n);
    acb_mat_sub(D, A, B, EV_PREC);
    arb_t nr;
    arb_init(nr);
    aic_mat_opnorm(nr, D, EV_PREC);
    double ub = arf_get_d(arb_midref(nr), ARF_RND_UP) + mag_get_d(arb_radref(nr));
    arb_clear(nr);
    acb_mat_clear(D);
    return ub;
}

static void test_s2_oracle(void)
{
    printf("S2: eta=0 exact-projection oracle\n");

    /* (a) rank-r projector P in C^6 (r=4). Pi from the lambda~1 cluster == P;
     * Pi from the lambda~0 cluster == I - P; trace(Pi_1) == r. */
    const slong n = 6, r = 4;
    acb_mat_t P;
    acb_mat_init(P, n, n);
    build_projector(P, n, r, EV_PREC);
    aic_mat_eigcluster *cl;
    slong nc;
    aic_mat_eig_hermitian_subspaces(&cl, &nc, P, -1.0, EV_PREC);
    AIC_CHECK_MSG(nc == 2, "(a) projector: %ld clusters != 2", (long) nc);

    slong c1 = s2_cluster_near(cl, nc, 1.0);
    slong c0 = s2_cluster_near(cl, nc, 0.0);
    AIC_CHECK_MSG(cl[c1].k == r, "(a) lambda~1 cluster k=%ld != r=%ld",
                  (long) cl[c1].k, (long) r);

    acb_mat_t Pi1, Pi0;
    acb_mat_init(Pi1, n, n);
    acb_mat_init(Pi0, n, n);
    aic_mat_cluster_projector(Pi1, &cl[c1], EV_PREC);
    aic_mat_cluster_projector(Pi0, &cl[c0], EV_PREC);

    double d1 = s2_opdiff(Pi1, P);
    AIC_CHECK_MSG(d1 < 1e-25, "(a) ||Pi_1 - P||_op = %.3e not < 1e-25", d1);

    /* I - P */
    acb_mat_t ImP;
    acb_mat_init(ImP, n, n);
    acb_mat_one(ImP);
    acb_mat_sub(ImP, ImP, P, EV_PREC);
    double d0 = s2_opdiff(Pi0, ImP);
    AIC_CHECK_MSG(d0 < 1e-25, "(a) ||Pi_0 - (I-P)||_op = %.3e not < 1e-25", d0);

    /* trace(Pi_1) == r */
    acb_t tr;
    acb_init(tr);
    acb_mat_trace(tr, Pi1, EV_PREC);
    arb_t trre, rdiff;
    arb_init(trre);
    arb_init(rdiff);
    acb_get_real(trre, tr);
    arb_sub_si(rdiff, trre, r, EV_PREC);
    arb_abs(rdiff, rdiff);
    double td = arf_get_d(arb_midref(rdiff), ARF_RND_UP) +
                mag_get_d(arb_radref(rdiff));
    AIC_CHECK_MSG(td < 1e-20, "(a) |trace(Pi_1) - r| = %.3e not < 1e-20", td);
    printf("  (a) rank-%ld proj C^%ld: ||Pi_1-P||=%.3e ||Pi_0-(I-P)||=%.3e "
           "trace(Pi_1)=%ld\n", (long) r, (long) n, d1, d0, (long) r);

    arb_clear(rdiff);
    arb_clear(trre);
    acb_clear(tr);
    acb_mat_clear(ImP);
    acb_mat_clear(Pi0);
    acb_mat_clear(Pi1);
    aic_mat_eigcluster_free(cl, nc);
    acb_mat_clear(P);

    /* (b) depolarizing Choi (1/d) I_{d^2}: ONE cluster, k=d^2, Pi_M = I. */
    const slong d = 3, dd = d * d;
    acb_mat_t C;
    acb_mat_init(C, dd, dd);
    acb_mat_zero(C);
    acb_t inv;
    acb_init(inv);
    acb_set_si(inv, 1);
    acb_div_si(inv, inv, d, EV_PREC);
    for (slong i = 0; i < dd; i++) acb_set(acb_mat_entry(C, i, i), inv);
    aic_mat_eigcluster *cl2;
    slong nc2;
    aic_mat_eig_hermitian_subspaces(&cl2, &nc2, C, -1.0, EV_PREC);
    AIC_CHECK_MSG(nc2 == 1, "(b) depolarizing: %ld clusters != 1", (long) nc2);
    AIC_CHECK_MSG(cl2[0].k == dd, "(b) depolarizing cluster k=%ld != %ld",
                  (long) cl2[0].k, (long) dd);
    acb_mat_t PiI, Id;
    acb_mat_init(PiI, dd, dd);
    acb_mat_init(Id, dd, dd);
    acb_mat_one(Id);
    aic_mat_cluster_projector(PiI, &cl2[0], EV_PREC);
    double di = s2_opdiff(PiI, Id);
    AIC_CHECK_MSG(di < 1e-20, "(b) ||Pi_M - I||_op = %.3e not < 1e-20", di);
    printf("  (b) depolarizing (1/%ld)I_{%ld}: one cluster k=%ld, ||Pi_M-I||=%.3e\n",
           (long) d, (long) dd, (long) dd, di);
    acb_mat_clear(Id);
    acb_mat_clear(PiI);
    aic_mat_eigcluster_free(cl2, nc2);
    acb_clear(inv);
    acb_mat_clear(C);
}

/* --- the fork+SIGALRM watchdog (pattern from test_eigmult.c) ---------------- */
#define EV_WATCH_S 20
typedef void (*ev_child_fn)(void);
static volatile pid_t ev_watch_pid = 0;
static volatile sig_atomic_t ev_timed_out = 0;
static void ev_alarm(int sig)
{
    (void) sig;
    ev_timed_out = 1;
    if (ev_watch_pid > 0) kill(ev_watch_pid, SIGKILL);
}

static int ev_run_child(ev_child_fn fn, int *status, char *err, size_t errsz)
{
    char tmpl[] = "/tmp/aic_eigvec_err_XXXXXX";
    int efd = mkstemp(tmpl);
    AIC_CHECK_MSG(efd >= 0, "mkstemp failed");
    fflush(NULL);
    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "fork failed");
    if (pid == 0) {
        dup2(efd, STDERR_FILENO);
        close(efd);
        fn();
        _exit(0);    /* reached only if fn did NOT abort (the bug we guard) */
    }
    ev_watch_pid = pid;
    ev_timed_out = 0;
    struct sigaction sa, old;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = ev_alarm;
    sigaction(SIGALRM, &sa, &old);
    alarm(EV_WATCH_S);
    int st = 0;
    pid_t w;
    do { w = waitpid(pid, &st, 0); } while (w < 0 && !ev_timed_out);
    alarm(0);
    sigaction(SIGALRM, &old, NULL);
    if (w < 0) waitpid(pid, &st, 0);
    lseek(efd, 0, SEEK_SET);
    ssize_t rd = read(efd, err, errsz - 1);
    err[(rd > 0) ? (size_t) rd : 0] = '\0';
    close(efd);
    unlink(tmpl);
    *status = st;
    return ev_timed_out ? 0 : 1;
}

static void ev_assert_failloud(ev_child_fn fn, const char *what,
                               const char *needle)
{
    char err[4096];
    int st = 0;
    int finished = ev_run_child(fn, &st, err, sizeof err);
    AIC_CHECK_MSG(finished, "%s HUNG (watchdog killed it after %d s)",
                  what, EV_WATCH_S);
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
    printf("  abort message: %.200s\n", err);
}

/* --- S7(c): a sub-resolution near-degeneracy must FAIL LOUD (NON-FINITE) ---- */
/* Eigenvalues 1 and 1+2^-40, FORCED into two clusters by a tiny gap_thr, at
 * prec=24: the tiny-gap cluster's seed eigenvectors are near-parallel, so the
 * Rump enclosure is NON-FINITE ([+/-inf]) -> the FINITE-enclosure guard (i) in
 * aic_mat_int_certify_cluster aborts ("UNRESOLVED"). This is the REACHABLE
 * fail-loud path for a near-degeneracy below the working precision (the overlap
 * gate (iii) is unreachable here — Rump self-isolates, FINDINGS §D5n2). Probed:
 * finite both AND overlap never co-occur for this scheme. */
#define S7_PREC 24
static void child_unresolved(void)
{
    const slong n = 2;
    acb_mat_t H;
    acb_mat_init(H, n, n);
    acb_mat_zero(H);
    acb_t g, v;
    acb_init(g);
    acb_init(v);
    acb_one(g);
    acb_mul_2exp_si(g, g, -40);              /* 2^-40 */
    acb_one(v);
    acb_add(v, v, g, S7_PREC);               /* 1 + 2^-40 */
    acb_set_si(acb_mat_entry(H, 0, 0), 1);
    acb_set(acb_mat_entry(H, 1, 1), v);
    aic_mat_eigcluster *cl;
    slong nc;
    /* gap_thr = 1e-30 forces a 2-cluster split; the tiny cluster's Rump
     * enclosure is non-finite at prec=24 -> the finite guard aborts. */
    aic_mat_eig_hermitian_subspaces(&cl, &nc, H, 1e-30, S7_PREC);
    aic_mat_eigcluster_free(cl, nc);
    acb_clear(v);
    acb_clear(g);
    acb_mat_clear(H);
}

static void test_s7c_unresolved_failloud(void)
{
    printf("S7(c): a sub-2^-prec near-degeneracy must FAIL LOUD (Rump enclosure "
           "NON-FINITE -> UNRESOLVED)\n");
    ev_assert_failloud(child_unresolved,
                       "aic_mat_eig_hermitian_subspaces(unresolved)",
                       "UNRESOLVED");
}

int main(void)
{
    test_s1_residual();
    test_s2_oracle();
    test_s7c_unresolved_failloud();
    aic_test_report("test_eigvec");
    printf("OK test_eigvec\n");
    return 0;
}
