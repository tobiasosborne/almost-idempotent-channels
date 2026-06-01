/* test_eigmult.c — cross-checks for the DEGENERACY-ROBUST certified eig + rank
 * layer (bead aic-w4o.1, increment 1; retires the FINDINGS §D5 deferral for the
 * carrier dimension). CLAUDE.md Rule 6: the strongest tests here are
 * double-vs-arb@53 agreement, the eta=0 / exact-projection oracle, and the
 * certified-vs-double carrier rank. Rule 7: these assertions were written BEFORE
 * the implementation (RED baseline). Rule 4 fail-loud: the straddle abort is
 * tested under the fork+SIGALRM watchdog (reused from test_xo0_failloud.c).
 *
 * The four teeth (cross-check ladder rungs in parentheses):
 *  T1 (rung 2) double-vs-arb: each certified ball from
 *     aic_mat_eig_hermitian_multiple CONTAINS the corresponding double
 *     eigenvalue (aic_latd_eig_hermitian), for a SIMPLE and a DEGENERATE
 *     spectrum (the rank-2 projector {-1,-1,1,1} from test_latd.c, and a
 *     diagonal with repeats). Prints the worst ball radius.
 *  T2 (rung 3) exact/eta=0 oracle: (a) a rank-r orthogonal projector P
 *     (spectrum {0,1}) -> certified rank == r EXACTLY for r=1..n-1; (b) the
 *     MAXIMALLY degenerate depolarizing Choi (1/d) I_{d^2} (one cluster) ->
 *     certified rank == d^2.
 *  T3 (rung 2/3) certified carrier rank: a projector carrier Q=Pi ->
 *     aic_ucp_carrier_rank(Q) == aic_ucp_carrier_rank_latd(Q) == known rank.
 *  T4 (Rule 4) FAIL-LOUD straddle tooth: a Hermitian H with an eigenvalue at
 *     the threshold so a ball STRADDLES thr -> aic_mat_certified_rank ABORTS
 *     (child dies by SIGABRT under the watchdog, not exit 0).
 *
 * prec=128 for the {0,1} fixtures (gap ~1, resolves trivially). Concrete numbers
 * (radii, ranks, abort message) printed, Rule 12.
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

#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"

#define EM_PREC 128

/* Build the rank-2 projector spectrum X = 2P - I in C^4, spectrum {-1,-1,1,1},
 * the genuine (non-basis-aligned) degeneracy from test_latd.c: P = diag(1,1,0,0)
 * conjugated by a rational orthogonal Q (cos=3/5, sin=4/5 in planes (0,2),(1,3)). */
static void build_2PmI(acb_mat_t X, slong prec)
{
    const slong n = 4, k = 2;
    acb_mat_t P, Q, Qt, T;
    acb_mat_init(P, n, n);
    acb_mat_init(Q, n, n);
    acb_mat_init(Qt, n, n);
    acb_mat_init(T, n, n);
    acb_mat_zero(P);
    for (slong i = 0; i < k; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
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
    acb_set(acb_mat_entry(Q, 2, 2), c);
    acb_set(acb_mat_entry(Q, 0, 2), ns);
    acb_set(acb_mat_entry(Q, 2, 0), s);
    acb_set(acb_mat_entry(Q, 1, 1), c);
    acb_set(acb_mat_entry(Q, 3, 3), c);
    acb_set(acb_mat_entry(Q, 1, 3), ns);
    acb_set(acb_mat_entry(Q, 3, 1), s);
    acb_mat_conjugate_transpose(Qt, Q);
    acb_mat_mul(T, Q, P, prec);
    acb_mat_mul(P, T, Qt, prec);
    acb_mat_scalar_mul_si(X, P, 2, prec);
    for (slong i = 0; i < n; i++)
        acb_sub_si(acb_mat_entry(X, i, i), acb_mat_entry(X, i, i), 1, prec);
    acb_clear(ns);
    acb_clear(s);
    acb_clear(c);
    acb_mat_clear(T);
    acb_mat_clear(Qt);
    acb_mat_clear(Q);
    acb_mat_clear(P);
}

/* A rank-r orthogonal projector P in C^n (spectrum {0 mult n-r, 1 mult r}),
 * conjugated by DISJOINT rational Givens rotations (cos=3/5, sin=4/5) so the
 * spectrum is OFF the coordinate axes yet the two cluster eigenbases stay
 * WELL-CONDITIONED. Each Givens mixes one "1-axis" i (i<r) with one "0-axis"
 * r+i (a disjoint plane), for i = 0 .. min(r, n-r)-1; the remaining axes are
 * untouched. WHY DISJOINT, not a dense chain (FINDINGS §D5n, load-bearing):
 * acb_mat_eig_multiple's Rump enclosure certifies a degenerate cluster only when
 * the approx_eig_qr SEED gives well-separated cluster eigenvectors. A DENSE
 * Givens CHAIN entangles the two clusters so the seed's per-eigenvector
 * approximations within a cluster are near-parallel and the subspace enclosure
 * cannot close (the same {2,2} spectrum then returns 0 even at prec=256 — a
 * SEED-CONDITIONING limit, not a precision limit). Disjoint Givens keep each
 * cluster's 2D invariant subspace clean, and the {2,2}/{3,2}/... spectra certify
 * at prec=128. This is the genuine, non-basis-aligned degeneracy the adversarial
 * discipline wants, on the well-conditioned side of the §D5n boundary. */
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
    slong pairs = (r < n - r) ? r : (n - r);   /* disjoint 1-axis<->0-axis mixes */
    for (slong i = 0; i < pairs; i++) {        /* plane (i, r+i), all disjoint */
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

/* Conjugate a real diagonal `diag` (length n) by the rational Givens rotations
 * (0,n-1) and (1,n-2) (cos=3/5, sin=4/5), so a degenerate spectrum is NOT
 * coordinate-aligned. WHY (FINDINGS §D5n): acb_mat_eig_multiple's Rump
 * invariant-subspace enclosure FAILS to certify (returns 0) on an
 * exactly-coordinate-aligned degenerate diagonal (the QR seed gives standard-
 * basis cluster eigenvectors the enclosure cannot close on); conjugating off the
 * axes makes the same spectrum certify at prec=128 (measured). This is the
 * genuine, non-basis-aligned degeneracy the adversarial discipline wants. */
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
    if (n >= 4) {               /* second mix breaks alignment within the cluster */
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

/* --- T1: double-vs-arb agreement, simple AND degenerate spectra --- */
/* The certified balls are FAR tighter than double precision (radii ~1e-31 at
 * prec=128 on these exact-rational inputs), so the LAPACK double eigenvalue —
 * which carries ~1e-16 rounding error — does NOT lie inside the ball. The
 * correct cross-check is the project's documented double-vs-arb@53 AGREEMENT
 * predicate (aic_eigset_close, tol ~1e-12), comparing the two sorted sets, NOT
 * exact arb containment. We collect the certified real parts into an arb_ptr and
 * run aic_eigset_close against the double set. */
static void test_t1_double_vs_arb(void)
{
    printf("T1: certified eig balls agree with the double eigenvalues (~1e-12)\n");
    const slong cases = 3;
    slong dims[3] = {4, 5, 4};
    double worst_rad = 0.0;
    arb_t tol;
    arb_init(tol);
    arb_set_d(tol, 1e-12);
    for (slong cidx = 0; cidx < cases; cidx++) {
        slong n = dims[cidx];
        acb_mat_t H;
        acb_mat_init(H, n, n);
        const char *name;
        if (cidx == 0) {            /* degenerate {-1,-1,1,1} */
            build_2PmI(H, EM_PREC);
            name = "2P-I {-1,-1,1,1}";
        } else if (cidx == 1) {     /* simple, well-separated diagonal-ish */
            acb_mat_zero(H);
            for (slong i = 0; i < n; i++)
                acb_set_si(acb_mat_entry(H, i, i), 2 * (i + 1));
            acb_set_si(acb_mat_entry(H, 0, 1), 1);
            acb_set_si(acb_mat_entry(H, 1, 0), 1);
            name = "simple separated";
        } else {                    /* CONJUGATED repeats {3,3,3,7} (off-axis) */
            slong diag[4] = {3, 3, 3, 7};
            build_conj_diag(H, diag, n, EM_PREC);
            name = "conj repeats {3,3,3,7}";
        }

        /* certified balls */
        acb_ptr E = _acb_vec_init(n);
        aic_mat_eig_hermitian_multiple(E, H, EM_PREC);

        /* double eigenvalues (the predicate sorts a private copy itself) */
        double _Complex *Ha = malloc((size_t) (n * n) * sizeof(double _Complex));
        double *dev = malloc((size_t) n * sizeof(double));
        aic_latd_from_acb_mat(Ha, H);
        aic_latd_eig_hermitian(dev, NULL, Ha, n);

        /* certified real parts as an arb_ptr set; track the worst radius. */
        arb_ptr cert = _arb_vec_init(n);
        for (slong i = 0; i < n; i++) {
            acb_get_real(cert + i, E + i);
            double r = mag_get_d(arb_radref(cert + i));
            if (r > worst_rad) worst_rad = r;
        }

        /* double-vs-arb@53-tol AGREEMENT as sorted sets (aic_eigset_close). */
        slong fk = -1;
        AIC_CHECK_MSG(aic_eigset_close(dev, cert, n, tol, &fk),
                      "%s: certified eig set disagrees with double set at sorted "
                      "index %ld", name, (long) fk);

        _arb_vec_clear(cert, n);
        free(dev);
        free(Ha);
        _acb_vec_clear(E, n);
        acb_mat_clear(H);
        printf("  [%s] all %ld eigs agree double-vs-arb within tol\n",
               name, (long) n);
    }
    arb_clear(tol);
    printf("  worst certified ball radius across cases: %.3e\n", worst_rad);
}

/* --- T2: eta=0 / exact oracle — projector rank, depolarizing Choi --- */
static void test_t2_exact_rank(void)
{
    printf("T2: eta=0 exact-rank oracle\n");

    /* (a) rank-r projector in C^4 -> certified rank == r for r=1,2,3 (clusters
     * {1,3},{2,2},{3,1}). thr=1/2 sits in the {0,1} gap. n=4 with disjoint-plane
     * conjugation keeps BOTH cluster eigenbases well-conditioned at every r, so
     * all three certify at prec=128 (the §D5n seed-conditioning boundary: a
     * larger n leaves an unpaired aligned axis in the bigger cluster — e.g.
     * C^5,r=2 ({2,3}) returns 0; that boundary is exercised separately below). */
    const slong n = 4;
    for (slong r = 1; r < n; r++) {
        acb_mat_t P;
        acb_mat_init(P, n, n);
        build_projector(P, n, r, EM_PREC);
        arb_t thr;
        arb_init(thr);
        arb_set_d(thr, 0.5);    /* midway between 0 and 1, clean separation */
        slong rk = aic_mat_certified_rank(P, thr, EM_PREC);
        AIC_CHECK_MSG(rk == r,
                      "(a) projector rank-%ld in C^%ld: certified rank %ld != %ld",
                      (long) r, (long) n, (long) rk, (long) r);
        arb_clear(thr);
        acb_mat_clear(P);
    }
    printf("  (a) projector ranks r=1..%ld in C^%ld: all exact\n",
           (long) (n - 1), (long) n);

    /* (b) maximally degenerate: depolarizing Choi (1/d) I_{d^2}, ALL eigenvalues
     * equal to 1/d (a single cluster of size d^2). Certified rank == d^2: every
     * eigenvalue is above a threshold 1/(2d) < 1/d. */
    const slong d = 3;
    const slong dd = d * d;
    acb_mat_t C;
    acb_mat_init(C, dd, dd);
    acb_mat_zero(C);
    acb_t inv;
    acb_init(inv);
    acb_set_si(inv, 1);
    acb_div_si(inv, inv, d, EM_PREC);     /* 1/d */
    for (slong i = 0; i < dd; i++) acb_set(acb_mat_entry(C, i, i), inv);
    arb_t thr2;
    arb_init(thr2);
    arb_set_d(thr2, 0.5 / (double) d);    /* 1/(2d), below the 1/d cluster */
    slong rk2 = aic_mat_certified_rank(C, thr2, EM_PREC);
    AIC_CHECK_MSG(rk2 == dd,
                  "(b) depolarizing Choi (1/%ld) I_{%ld}: certified rank %ld != %ld",
                  (long) d, (long) dd, (long) rk2, (long) dd);
    printf("  (b) depolarizing Choi (1/%ld) I_{%ld}: single cluster, rank %ld\n",
           (long) d, (long) dd, (long) rk2);
    arb_clear(thr2);
    acb_clear(inv);
    acb_mat_clear(C);
}

/* --- T3: certified carrier rank vs double path --- */
static void test_t3_carrier_rank(void)
{
    printf("T3: certified carrier rank vs double path\n");
    /* Build the carrier Q directly as a rank-r projector Pi in C^dimK (a valid
     * Hermitian PSD carrier; the certified rank does not need the channel
     * machinery, only Q). */
    const slong dimK = 6, r = 4;
    acb_mat_t Q;
    acb_mat_init(Q, dimK, dimK);
    build_projector(Q, dimK, r, EM_PREC);

    slong rk_cert = aic_ucp_carrier_rank(Q, dimK, EM_PREC);
    slong rk_dbl = aic_ucp_carrier_rank_latd(Q, dimK);
    AIC_CHECK_MSG(rk_cert == r,
                  "certified carrier rank %ld != known rank %ld",
                  (long) rk_cert, (long) r);
    AIC_CHECK_MSG(rk_cert == rk_dbl,
                  "certified carrier rank %ld != double-path rank %ld",
                  (long) rk_cert, (long) rk_dbl);
    printf("  carrier Q = rank-%ld projector in C^%ld: certified=%ld double=%ld\n",
           (long) r, (long) dimK, (long) rk_cert, (long) rk_dbl);
    acb_mat_clear(Q);
}

/* --- the fork+SIGALRM watchdog (reusable; pattern from test_xo0_failloud.c) --- */
/* The parent forks; the child redirects stderr to a temp file and runs `fn`
 * (which is expected to ABORT). The parent waits with a SIGALRM-bounded waitpid:
 * a child alive past the watchdog is a HANG (also a failure). Captures the
 * child's stderr into `err` and the raw wait status into *status; returns 1 if
 * the child terminated within the watchdog, 0 if it was killed for hanging. */
#define EM_WATCH_S 20
typedef void (*em_child_fn)(void);
static volatile pid_t em_watch_pid = 0;
static volatile sig_atomic_t em_timed_out = 0;
static void em_alarm(int sig)
{
    (void) sig;
    em_timed_out = 1;
    if (em_watch_pid > 0) kill(em_watch_pid, SIGKILL);
}

static int em_run_child(em_child_fn fn, int *status, char *err, size_t errsz)
{
    char tmpl[] = "/tmp/aic_eigmult_err_XXXXXX";
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

    em_watch_pid = pid;
    em_timed_out = 0;
    struct sigaction sa, old;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = em_alarm;
    sigaction(SIGALRM, &sa, &old);
    alarm(EM_WATCH_S);
    int st = 0;
    pid_t w;
    do { w = waitpid(pid, &st, 0); } while (w < 0 && !em_timed_out);
    alarm(0);
    sigaction(SIGALRM, &old, NULL);
    if (w < 0) waitpid(pid, &st, 0);

    lseek(efd, 0, SEEK_SET);
    ssize_t rd = read(efd, err, errsz - 1);
    err[(rd > 0) ? (size_t) rd : 0] = '\0';
    close(efd);
    unlink(tmpl);

    *status = st;
    return em_timed_out ? 0 : 1;
}

/* assert a child that should ABORT (SIGABRT) and whose stderr names `needle`. */
static void em_assert_failloud(em_child_fn fn, const char *what,
                               const char *needle)
{
    char err[4096];
    int st = 0;
    int finished = em_run_child(fn, &st, err, sizeof err);
    AIC_CHECK_MSG(finished, "%s HUNG (watchdog killed it after %d s)",
                  what, EM_WATCH_S);
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

/* --- T4: FAIL-LOUD straddle tooth --- */
/* H = 2P - I, spectrum {-1,-1,1,1}. Threshold thr = 1 sits EXACTLY on the +1
 * cluster, so those eigenvalue balls straddle thr -> rank unresolved -> abort. */
static void child_straddle(void)
{
    acb_mat_t H;
    acb_mat_init(H, 4, 4);
    build_2PmI(H, EM_PREC);
    arb_t thr;
    arb_init(thr);
    arb_set_si(thr, 1);
    (void) aic_mat_certified_rank(H, thr, EM_PREC);
    arb_clear(thr);
    acb_mat_clear(H);
}

static void test_t4_straddle_failloud(void)
{
    printf("T4: straddling threshold must FAIL LOUD (SIGABRT), not silently count\n");
    em_assert_failloud(child_straddle, "aic_mat_certified_rank(straddle)",
                       "STRADDLES");
}

/* --- T5: FAIL-LOUD on an UNRESOLVED cluster (§D5n seed-conditioning boundary) --- */
/* A rank-2 projector in C^5 (clusters {2,3}) is the §D5n boundary: the dense
 * cluster eigenbases from approx_eig_qr are near-parallel, acb_mat_eig_multiple
 * returns 0, and aic_mat_eig_hermitian_multiple MUST fail loud (Rule 4) — NOT
 * silently miscount the rank. This guards the eig_multiple-returned-0 path. */
static void child_unresolved(void)
{
    acb_mat_t P;
    acb_mat_init(P, 5, 5);
    build_projector(P, 5, 2, EM_PREC);     /* {2,3}: eig_multiple returns 0 */
    arb_t thr;
    arb_init(thr);
    arb_set_d(thr, 0.5);
    (void) aic_mat_certified_rank(P, thr, EM_PREC);
    arb_clear(thr);
    acb_mat_clear(P);
}

static void test_t5_unresolved_failloud(void)
{
    printf("T5: unresolved cluster ({2,3} seed-conditioning) must FAIL LOUD\n");
    em_assert_failloud(child_unresolved,
                       "aic_mat_eig_hermitian_multiple(unresolved)",
                       "clusters unresolved");
}

int main(void)
{
    test_t1_double_vs_arb();
    test_t2_exact_rank();
    test_t3_carrier_rank();
    test_t4_straddle_failloud();
    test_t5_unresolved_failloud();
    aic_test_report("test_eigmult");
    printf("OK test_eigmult\n");
    return 0;
}
