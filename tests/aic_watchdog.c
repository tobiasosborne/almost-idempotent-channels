#define _POSIX_C_SOURCE 200809L
/* aic_watchdog.c — the canonical fork + SIGALRM fail-loud watchdog (bead aic-8de).
 *
 * _POSIX_C_SOURCE is the VERY FIRST line, before ANY #include: fork, alarm,
 * sigaction, mkstemp, dup2, lseek, and waitpid are POSIX, not C11, so they are
 * only declared once the feature-test macro is set. Keeping the macro HERE (and
 * only here) is the whole point of the header/.c split — a consumer that links
 * aic_watchdog needs no POSIX macro of its own (aic_watchdog.h pulls in only
 * <stddef.h>).
 *
 * WHY this file exists. CLAUDE.md Rule 4 demands routines abort() on a bad input
 * rather than silently return corrupted output. A test that asserts "this call
 * aborts" cannot run the call in-process — the abort() would kill the test
 * driver. So aic_watchdog_run() runs fn() in a forked child and the parent reads
 * the child's wait-status: WIFSIGNALED && WTERMSIG == SIGABRT is the fail-loud
 * success. A SIGALRM watchdog backstops a HANG (a non-converging iteration that
 * never aborts): the alarm handler SIGKILLs the child and we report a hang
 * instead of letting the gate stall forever.
 *
 * This unifies the nine near-identical private copies (fc_watch in test_funcalc,
 * xo0_run_child in test_xo0_failloud, v5f_run_child in test_adversarial, and the
 * EV_/EM_/KA_/UCP_/CHAN_ siblings). Increment 1 (bead aic-8de) migrates the
 * test_funcalc reference consumer; the other eight follow in later increments.
 *
 * Subtleties this code guards (each was load-bearing in the originals):
 *  - _exit(0), NOT exit(0), on the child's normal-return path: exit() would run
 *    the parent's inherited FLINT/stdio atexit handlers in the child and flush
 *    state that is not the child's to flush.
 *  - fflush(NULL) BEFORE fork(): otherwise both processes flush the same buffered
 *    bytes and the captured stderr is doubled.
 *  - mkstemp() BEFORE fork() so the child inherits the open fd; the child does
 *    dup2(fd, STDERR_FILENO) (and STDOUT if asked); the parent reads it back
 *    AFTER the child terminates, then unlinks the temp path.
 *  - SIGALRM via sigaction with the OLD handler saved and restored, so the
 *    watchdog does not leak a handler into later tests in the same driver.
 *  - file-static volatile sig_atomic_t timeout flag (the type the C standard
 *    blesses for handler access) plus a volatile pid_t child pid: kill() needs a
 *    pid_t, and on POSIX pid_t is an int whose aligned load is atomic, with the
 *    parent's write happening-before the handler is ever installed — safe in
 *    practice even though only sig_atomic_t is standard-guaranteed.
 */
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "aic_test.h"      /* AIC_FAIL / AIC_CHECK_MSG — the suite's fail-loud idiom */
#include "aic_watchdog.h"

/* The only state the SIGALRM handler may touch: the child pid to kill and a
 * flag the parent polls. The flag is volatile sig_atomic_t (standard-blessed for
 * handler access); the pid is volatile pid_t because kill() needs a pid_t — safe
 * in practice on POSIX (int, atomic aligned load, written before handler install).
 * File-static so the handler (which takes no argument) can reach them. */
static volatile pid_t aic_wd_child = 0;
static volatile sig_atomic_t aic_wd_timed_out = 0;

static void aic_wd_alarm(int sig)
{
    (void) sig;
    aic_wd_timed_out = 1;
    if (aic_wd_child > 0) kill(aic_wd_child, SIGKILL);
}

int aic_watchdog_run(aic_watchdog_fn fn, int timeout_s, int capture_stdout,
                     int *status_out, char *err, size_t errsz)
{
    AIC_CHECK_MSG(err != NULL && errsz > 0,
                  "aic_watchdog_run: err buffer must be non-empty");
    char tmpl[] = "/tmp/aic_watchdog_err_XXXXXX";
    int efd = mkstemp(tmpl);
    AIC_CHECK_MSG(efd >= 0, "aic_watchdog_run: mkstemp failed");

    fflush(NULL);                  /* flush BEFORE fork: no doubled buffered output */
    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "aic_watchdog_run: fork failed");
    if (pid == 0) {
        /* child: redirect stderr (and optionally stdout) to the inherited fd */
        dup2(efd, STDERR_FILENO);
        if (capture_stdout) dup2(efd, STDOUT_FILENO);
        close(efd);
        fn();
        _exit(0);                  /* reached ONLY if fn did NOT abort (silent) */
    }

    aic_wd_child = pid;
    aic_wd_timed_out = 0;
    struct sigaction sa, old;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = aic_wd_alarm;
    sigaction(SIGALRM, &sa, &old);     /* save the OLD handler */
    alarm((unsigned) (timeout_s > 0 ? timeout_s : 0));

    int st = 0;
    pid_t w;
    do { w = waitpid(pid, &st, 0); }
    while (w < 0 && errno == EINTR && !aic_wd_timed_out);
    alarm(0);
    sigaction(SIGALRM, &old, NULL);    /* RESTORE it: no handler leak */
    if (w < 0) waitpid(pid, &st, 0);   /* reap the SIGKILLed child on timeout */

    lseek(efd, 0, SEEK_SET);
    ssize_t r = read(efd, err, errsz - 1);
    err[(r > 0) ? (size_t) r : 0] = '\0';
    close(efd);
    unlink(tmpl);

    *status_out = st;
    return aic_wd_timed_out ? 0 : 1;
}

void aic_watchdog_assert_failloud(aic_watchdog_fn fn, int timeout_s,
                                  const char *who, const char *needle)
{
    int status = 0;
    char err[4096];
    if (who == NULL) who = "(unnamed)";
    /* Capture BOTH streams: FLINT aborts via flint_printf, which writes the
     * fail-loud message to STDOUT, not stderr. A stderr-only capture would miss
     * every flint_abort()/flint_printf needle (the project's primary abort path),
     * so a fail-loud message check must read both. */
    int finished = aic_watchdog_run(fn, timeout_s, 1, &status, err, sizeof err);

    AIC_CHECK_MSG(finished, "%s: call HUNG past %ds watchdog "
                  "(expected a fast Rule-4 abort, not a hang)", who, timeout_s);

    /* silent-exit check: exit 0 means the guard did NOT fire — corrupted output
     * was returned instead of an abort. That is a TEST FAILURE, not a pass. */
    AIC_CHECK_MSG(!(WIFEXITED(status) && WEXITSTATUS(status) == 0),
                  "%s: exited cleanly (code 0) — the fail-loud guard did NOT "
                  "fire (silent success == wrong)", who);

    AIC_CHECK_MSG(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT,
                  "%s: expected SIGABRT (fail-loud abort) but got "
                  "WIFSIGNALED=%d WTERMSIG=%d / WIFEXITED=%d code=%d. stderr: %s",
                  who, WIFSIGNALED(status),
                  WIFSIGNALED(status) ? WTERMSIG(status) : -1,
                  WIFEXITED(status), WIFEXITED(status) ? WEXITSTATUS(status) : -1,
                  err);

    if (needle != NULL)
        AIC_CHECK_MSG(strstr(err, needle) != NULL,
                      "%s: aborted but WITHOUT the expected message '%s' "
                      "(silent/wrong-reason abort?). stderr: %s",
                      who, needle, err);
}
