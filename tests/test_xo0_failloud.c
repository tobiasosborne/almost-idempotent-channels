/* test_xo0_failloud.c — the FAIL-LOUD (not HANG) guard for the almost-idempotent
 * pipeline entry aic_assoc_regularize (bead aic-xo0). CLAUDE.md Rule 4 (fail loud)
 * + Rule 5 ("runs without errors" is not a pass: assert the abort + its message).
 *
 * THE GAP THIS CLOSES. aic_assoc_regularize (src/aic_assoc_regularize.c) realizes
 * Phi_tilde = theta(2 Phi - 1) (approximate_algebras.tex:2171-2174), well-defined
 * only when eta < 1/4 (.tex:2176, "the Taylor expansion in 4(Phi-Phi^2) converges
 * if eta<1/4"; equivalently the prop_P basin ||P^2-P|| <= delta < 1/4, .tex:524-525,
 * i.e. ||(2P-I)^2-I|| = 4 rho(P^2-P) < 1, .tex:516,520). Feeding an OUT-OF-REGIME
 * channel (rho(S_Phi^2-S_Phi) >= 1/4) must FAIL LOUD at the API boundary — NOT hang,
 * NOT silently regularize a non-convergent series. Before bead aic-xo0 the only guard
 * was the DEEP one inside aic_prop_P (a generic "P" message with no channel-level
 * eta<1/4 attribution and no aic_assoc_regularize provenance); there was NO executable
 * test that the pipeline entry aborts (an abort crashes the binary, so it cannot be a
 * plain in-process assertion). This driver introduces the reusable SUBPROCESS pattern
 * (fork + bounded watchdog + captured stderr) so a Rule-4 abort is testable.
 *
 * WHY A UNITARY-CONJUGATION CHANNEL IS THE OUT-OF-REGIME FIXTURE (load-bearing).
 * The CLAUDE.md probe-hygiene note and this bead originally cited make_mixconj with
 * t~0.45 as out-of-regime, via the OP-NORM eta-proxy (~6.5 t). But make_mixconj is a
 * convex mix of an exactly-idempotent compression and its UNITARY conjugate — BOTH
 * summands are (spectrally) idempotent, so the mix stays SPECTRALLY near-idempotent
 * at EVERY t: measured rho(4(S^2-S)) <= 0.66 < 1 for t in [0, 1] (it is 0 at t=0 and
 * t=1). Post bead aic-8hz, aic_prop_P certifies the SPECTRAL basin rho(S^2-S) < 1/4
 * (NOT the op-norm one), so make_mixconj is IN regime for all t and is the WRONG
 * fixture for this guard. The genuine out-of-regime UCP self-map used here is the
 * *-automorphism Phi(X) = U X U^dag with U = diag(1,-1) (a Pauli-Z reflection): its
 * superoperator S = U (x) conj(U) has eigenvalues u_i conj(u_j) in {+1,-1}, so the
 * lambda = -1 eigenvalues give (S^2-S) the eigenvalue (-1)^2-(-1) = 2, i.e.
 * rho(S^2-S) = 2 >= 1/4 and rho(4(S^2-S)) = 8 (measured: Gelfand ub ~ 8.09, NOT
 * certified < 1). It is unital, CP, trace-preserving — a perfectly valid UCP map,
 * just NOT almost-idempotent. This is the bead's real adversarial input.
 *
 * THE POSITIVE NO-FALSE-POSITIVE CASES (guard must not over-reject). Two channels
 * that ARE in regime must pass aic_assoc_regularize without aborting:
 *   (a) eta=0: make_compress_idemp(4,2), exactly idempotent (rho(S^2-S)=0).
 *   (b) the SPECTRAL-RELAXATION regime (aic-8hz; test_assoc2 U6): a non-normal
 *       channel with ||S^2-S||_op >= 1/4 but rho(S^2-S) < 1/4. A cheap op-norm gate
 *       (||S^2-S||_op < 1/4) would WRONGLY reject this; the guard MUST use the
 *       spectral Gelfand certifier (rho), exactly as aic_prop_P does. Fixture:
 *       make_mix(compress_idemp(4,1), dephasing(4), t=0.30): measured ||S^2-S||_op
 *       = 0.420 >= 1/4 yet rho(S^2-S) = 0.210 < 1/4 (Gelfand certifies at k=6).
 *
 * THE WATCHDOG. The parent forks; the child redirects stderr to a temp file and runs
 * the pipeline call. The parent waits with a SIGALRM-bounded waitpid (AIC_XO0_WATCH_S
 * seconds): if the child has not terminated by then it is SIGKILLed and the test
 * FAILS (a hang IS a failure, Rule 4). On a clean termination the parent reads the
 * status + the captured stderr and asserts via the NDEBUG-immune AIC_CHECK harness.
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

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_funcalc.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "test_idemp.h"            /* make_compress_idemp, make_dephasing, PREC */

#define XO0_PREC 256
#define AIC_XO0_WATCH_S 20         /* watchdog: a child alive past this is a HANG */

/* ---- the subprocess watchdog (reusable for any Rule-4 abort guard) -------- */

typedef void (*xo0_child_fn)(void *ctx);

/* SIGALRM handler target: which child to kill if the watchdog fires. */
static volatile pid_t xo0_watch_pid = 0;
static volatile sig_atomic_t xo0_timed_out = 0;

static void xo0_alarm(int sig)
{
    (void) sig;
    xo0_timed_out = 1;
    if (xo0_watch_pid > 0) kill(xo0_watch_pid, SIGKILL);
}

/* Run `fn(ctx)` in a forked child with a bounded watchdog. Captures the child's
 * stderr into `err` (caller buffer of `errsz`, NUL-terminated). Writes the raw
 * wait status into *status. Returns 1 if the child terminated within the watchdog,
 * 0 if it was killed for hanging (caller asserts on this). */
static int xo0_run_child(xo0_child_fn fn, void *ctx, int *status,
                         char *err, size_t errsz)
{
    char tmpl[] = "/tmp/aic_xo0_err_XXXXXX";
    int efd = mkstemp(tmpl);
    AIC_CHECK_MSG(efd >= 0, "xo0_run_child: mkstemp failed");
    fflush(NULL);

    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "xo0_run_child: fork failed");
    if (pid == 0) {
        /* child: stderr -> temp file, run the pipeline call (may abort/hang). */
        dup2(efd, STDERR_FILENO);
        close(efd);
        fn(ctx);
        _exit(0);              /* reached only if fn did NOT abort */
    }

    /* parent: bounded waitpid via SIGALRM. */
    xo0_watch_pid = pid;
    xo0_timed_out = 0;
    struct sigaction sa = {0}, old;
    sa.sa_handler = xo0_alarm;
    sigaction(SIGALRM, &sa, &old);
    alarm(AIC_XO0_WATCH_S);

    int st = 0;
    pid_t w;
    do { w = waitpid(pid, &st, 0); } while (w < 0 && !xo0_timed_out);
    alarm(0);
    sigaction(SIGALRM, &old, NULL);
    if (w < 0) waitpid(pid, &st, 0);   /* reap the SIGKILLed child */

    /* read captured stderr */
    lseek(efd, 0, SEEK_SET);
    ssize_t r = read(efd, err, errsz - 1);
    err[(r > 0) ? (size_t) r : 0] = '\0';
    close(efd);
    unlink(tmpl);

    *status = st;
    return xo0_timed_out ? 0 : 1;
}

/* ---- the child payloads --------------------------------------------------- */

/* ctx = a UCP self-map; the child regularizes it (the pipeline entry). */
static void child_regularize(void *ctx)
{
    aic_ucp_kraus *phi = (aic_ucp_kraus *) ctx;
    slong n = phi->dim_H;
    acb_mat_t St;
    acb_mat_init(St, n * n, n * n);
    aic_assoc_regularize(St, phi, XO0_PREC);
    acb_mat_clear(St);
}

/* ---- fixtures ------------------------------------------------------------- */

/* Phi_t = (1-t) base + t mix as the Kraus union (convex combo of UCP maps is UCP).
 * Copied from test_assoc2.c (local fixture, not in a shared header). */
static void make_mix(aic_ucp_kraus *out, const aic_ucp_kraus *base,
                     const aic_ucp_kraus *mix, double t, slong prec)
{
    slong n = base->dim_H;
    AIC_CHECK_MSG(base->dim_K == n && mix->dim_K == n && mix->dim_H == n,
                  "make_mix: base/mix not self-maps on the same C^n");
    aic_ucp_kraus_init(out, n, n, base->r + mix->r);
    arb_t s;
    arb_init(s);
    arb_set_d(s, sqrt(1.0 - t));
    for (slong a = 0; a < base->r; a++)
        acb_mat_scalar_mul_arb(out->K[a], base->K[a], s, prec);
    arb_set_d(s, sqrt(t));
    for (slong b = 0; b < mix->r; b++)
        acb_mat_scalar_mul_arb(out->K[base->r + b], mix->K[b], s, prec);
    arb_clear(s);
}

/* The OUT-OF-REGIME UCP self-map Phi(X) = U X U^dag, U = diag(1,-1,...) a Pauli-Z
 * reflection on the first two coords (rest +1). superop eigenvalues u_i conj(u_j)
 * in {+1,-1}; the -1 ones give rho(S^2-S) = 2 >= 1/4 (file docstring). `out` init'd
 * here. */
static void make_zconj(aic_ucp_kraus *out, slong n)
{
    aic_ucp_kraus_init(out, n, n, 1);
    acb_set_si(acb_mat_entry(out->K[0], 0, 0), 1);
    acb_set_si(acb_mat_entry(out->K[0], 1, 1), -1);
    for (slong i = 2; i < n; i++) acb_set_si(acb_mat_entry(out->K[0], i, i), 1);
}

/* ---- tests ---------------------------------------------------------------- */

/* NEG: the out-of-regime channel must FAIL LOUD (SIGABRT), NOT hang, and the abort
 * message must attribute the failure to aic_assoc_regularize and name the eta<1/4
 * precondition. RED before the entry guard: the only abort is the DEEP aic_prop_P
 * message, which says neither "aic_assoc_regularize" nor "eta < 1/4". */
static void test_oob_failloud(void)
{
    printf("NEG: out-of-regime Z-conjugation channel must fail loud (not hang)\n");
    aic_ucp_kraus phi;
    make_zconj(&phi, 2);

    char err[4096];
    int status = 0;
    int finished = xo0_run_child(child_regularize, &phi, &status, err, sizeof err);

    /* (1) not a hang */
    AIC_CHECK_MSG(finished,
                  "aic_assoc_regularize HUNG on an out-of-regime channel "
                  "(watchdog killed it after %d s) — Rule 4 fail-loud gap", AIC_XO0_WATCH_S);
    /* (2) it aborted (SIGABRT) or at least exited non-zero — NOT a silent success */
    int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
    int nonzero_exit = WIFEXITED(status) && WEXITSTATUS(status) != 0;
    AIC_CHECK_MSG(aborted || nonzero_exit,
                  "aic_assoc_regularize did NOT fail loud on an out-of-regime channel "
                  "(child exited 0 — it silently regularized a non-convergent series)");
    AIC_CHECK_MSG(aborted,
                  "aic_assoc_regularize terminated non-abort on out-of-regime input "
                  "(expected SIGABRT; signaled=%d sig=%d exited=%d code=%d)",
                  WIFSIGNALED(status), WIFSIGNALED(status) ? WTERMSIG(status) : -1,
                  WIFEXITED(status), WIFEXITED(status) ? WEXITSTATUS(status) : -1);
    /* (3) the entry message names the regularization basin / eta<1/4 + provenance.
     * RED before the fix (the deep prop_P message lacks both). */
    AIC_CHECK_MSG(strstr(err, "aic_assoc_regularize") != NULL,
                  "abort message does not attribute the failure to "
                  "aic_assoc_regularize (the API-boundary entry guard). Got:\n%s", err);
    AIC_CHECK_MSG(strstr(err, "eta < 1/4") != NULL || strstr(err, "eta<1/4") != NULL,
                  "abort message does not name the eta < 1/4 almost-idempotent "
                  "precondition (.tex:2168/2176). Got:\n%s", err);
    printf("  child SIGABRT, message names the eta<1/4 entry precondition (OK)\n");

    aic_ucp_kraus_clear(&phi);
}

/* POS: in-regime channels must pass aic_assoc_regularize WITHOUT aborting. Includes
 * the SPECTRAL-RELAXATION case (||S^2-S||_op >= 1/4 but rho < 1/4) that an op-norm
 * gate would wrongly reject (aic-8hz / U6). Guards against an over-aggressive entry
 * guard. */
static void test_in_regime_passes(void)
{
    printf("POS: in-regime channels must pass without aborting\n");

    /* (a) eta=0 exact idempotent. */
    {
        aic_ucp_kraus phi;
        make_compress_idemp(&phi, 4, 2);
        char err[4096];
        int status = 0;
        int finished = xo0_run_child(child_regularize, &phi, &status, err, sizeof err);
        AIC_CHECK_MSG(finished, "eta=0 compress_idemp(4,2) HUNG in regularize");
        AIC_CHECK_MSG(WIFEXITED(status) && WEXITSTATUS(status) == 0,
                      "eta=0 compress_idemp(4,2) FALSE-rejected by the entry guard "
                      "(child did not exit 0). stderr:\n%s", err);
        printf("  (a) eta=0 compress_idemp(4,2): passed (exit 0)\n");
        aic_ucp_kraus_clear(&phi);
    }

    /* (b) spectral-relaxation: ||S^2-S||_op >= 1/4 but rho(S^2-S) < 1/4 (U6). An
     * op-norm gate would wrongly abort; the spectral Gelfand guard must accept. */
    {
        aic_ucp_kraus base, dep, phi;
        make_compress_idemp(&base, 4, 1);
        make_dephasing(&dep, 4);
        make_mix(&phi, &base, &dep, 0.30, XO0_PREC);
        char err[4096];
        int status = 0;
        int finished = xo0_run_child(child_regularize, &phi, &status, err, sizeof err);
        AIC_CHECK_MSG(finished, "spectral-relaxation mix(4,1)+deph(4) HUNG in regularize");
        AIC_CHECK_MSG(WIFEXITED(status) && WEXITSTATUS(status) == 0,
                      "spectral-relaxation channel (||S^2-S||_op>=1/4 but rho<1/4) "
                      "FALSE-rejected by an over-aggressive (op-norm) entry guard — "
                      "the guard must use the SPECTRAL rho (aic-8hz). stderr:\n%s", err);
        printf("  (b) spectral-relaxation mix(compress(4,1),deph(4),0.30): passed (exit 0)\n");
        aic_ucp_kraus_clear(&phi);
        aic_ucp_kraus_clear(&dep);
        aic_ucp_kraus_clear(&base);
    }
}

int main(void)
{
    test_oob_failloud();
    test_in_regime_passes();
    aic_test_report("test_xo0_failloud");
    return 0;
}
