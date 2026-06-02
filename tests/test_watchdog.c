#define _POSIX_C_SOURCE 200809L
/* test_watchdog.c — the self-test for the shared fork+SIGALRM watchdog
 * (tests/aic_watchdog.{h,c}, bead aic-8de). CLAUDE.md Rule 5 ("runs without
 * errors is NOT a pass") + Rule 7 (mutation-proven teeth).
 *
 * WHY shared infra must self-test. aic_watchdog is the ONE piece of test
 * machinery every fail-loud guard in the suite leans on, yet its only consumer
 * — aic_watchdog_assert_failloud — always hands it a FAST-ABORTING child. So
 * three load-bearing branches were NEVER exercised in-suite and could rot
 * silently: the SIGALRM hang/timeout path (return 0 + SIGKILL status), the
 * capture_stdout redirect, and the discrimination between SIGABRT (a Rule-4
 * fail-loud abort), SIGSEGV (a REAL crash — must NOT be mistaken for an abort),
 * SIGKILL (a hang the watchdog killed), and a clean _exit(0) (silent success,
 * which assert_failloud must REJECT). This driver pins all four wait-status
 * shapes plus the stdout capture, using aic_watchdog_run directly (it RETURNS
 * the status, it does not abort), so the dead branches become live + regression-
 * guarded. assert_failloud is used ONLY on the happy abort+needle child.
 *
 * Sub-tests (one child payload each):
 *   1. child_abort  — abort() after a stderr marker: assert_failloud PASSES, and
 *      run() reports SIGABRT with the marker captured.
 *   2. child_clean  — empty body (helper _exit(0)): run() reports WIFEXITED code 0
 *      (the silent-success case assert_failloud rejects; we prove it is DISTINCT).
 *   3. child_hang   — for(;;)pause(): run() at a 1 s timeout returns finished==0
 *      and SIGKILL (the watchdog path, otherwise dead in-suite).
 *   4. child_segv   — raise(SIGSEGV): run() reports SIGSEGV, NOT SIGABRT (a real
 *      crash is not laundered into a fail-loud abort).
 *   5. child_stdout — write(2) to stdout then abort: capture_stdout=1 lands the
 *      stdout marker in the captured buffer.
 * No FLINT/arb needed; deps are POSIX signal/wait/unistd + AIC_CHECK.
 */
#include <signal.h>     /* SIGABRT/SIGSEGV/SIGKILL, raise */
#include <stdio.h>
#include <stdlib.h>     /* abort */
#include <string.h>     /* strstr */
#include <sys/wait.h>   /* WIF* / WTERMSIG / WEXITSTATUS */
#include <unistd.h>     /* pause, write, STDOUT_FILENO */

#include "aic_watchdog.h"
#include "aic_test.h"

/* ---- child payloads (each a file-static void(void), run in the fork) ------- */

/* 1: fail-loud abort with an identifiable stderr marker. */
static void child_abort(void)
{
    fputs("wd_self_marker_abc\n", stderr);
    abort();
}

/* 2: returns normally -> the helper reaches its _exit(0) (silent success). */
static void child_clean(void)
{
}

/* 3: hangs forever; the SIGALRM watchdog must SIGKILL it. */
static void child_hang(void)
{
    for (;;) pause();
}

/* 4: a REAL crash. Must be reported as SIGSEGV, NOT confused with SIGABRT. */
static void child_segv(void)
{
    raise(SIGSEGV);
}

/* 5: marker on STDOUT. write(2) (not printf) precisely because the helper ends
 * the child with _exit(), which does NOT flush buffered stdio — a printf marker
 * could be lost, so we write the raw bytes to the fd directly. */
static void child_stdout(void)
{
    ssize_t n = write(STDOUT_FILENO, "wd_stdout_marker", 16);
    (void) n;   /* best-effort; the parent verifies via the captured buffer */
    fputs("se\n", stderr);
    abort();
}

/* ---- tests ---------------------------------------------------------------- */

/* 1: the happy fail-loud path through BOTH entry points. */
static void test_abort_failloud(void)
{
    printf("1: child_abort -> SIGABRT + stderr needle\n");
    /* assert_failloud must NOT abort THIS driver on a genuinely-aborting child. */
    aic_watchdog_assert_failloud(child_abort, 20, "self/abort+needle",
                                 "wd_self_marker_abc");

    int st = 0;
    char err[4096];
    int finished = aic_watchdog_run(child_abort, 20, 0, &st, err, sizeof err);
    AIC_CHECK_MSG(finished == 1, "child_abort: run() reported a hang (finished=0)");
    AIC_CHECK_MSG(WIFSIGNALED(st) && WTERMSIG(st) == SIGABRT,
                  "child_abort: expected SIGABRT, got signaled=%d sig=%d exited=%d",
                  WIFSIGNALED(st), WIFSIGNALED(st) ? WTERMSIG(st) : -1, WIFEXITED(st));
    AIC_CHECK_MSG(strstr(err, "wd_self_marker_abc") != NULL,
                  "child_abort: stderr marker not captured. got:\n%s", err);
}

/* 2: a clean _exit(0) is detectable (distinguishable from an abort/crash). */
static void test_clean_exit(void)
{
    printf("2: child_clean -> WIFEXITED code 0 (silent success, distinguishable)\n");
    int st = 0;
    char err[4096];
    int finished = aic_watchdog_run(child_clean, 20, 0, &st, err, sizeof err);
    AIC_CHECK_MSG(finished == 1, "child_clean: run() reported a hang (finished=0)");
    AIC_CHECK_MSG(WIFEXITED(st) && WEXITSTATUS(st) == 0,
                  "child_clean: expected clean exit 0, got signaled=%d exited=%d code=%d",
                  WIFSIGNALED(st), WIFEXITED(st),
                  WIFEXITED(st) ? WEXITSTATUS(st) : -1);
}

/* 3: the watchdog catches a hang — return 0 AND status is SIGKILL. */
static void test_hang_detected(void)
{
    printf("3: child_hang -> finished==0 + SIGKILL (watchdog path)\n");
    int st = 0;
    char err[4096];
    int finished = aic_watchdog_run(child_hang, 1, 0, &st, err, sizeof err);
    AIC_CHECK_MSG(finished == 0,
                  "child_hang: watchdog did NOT detect the hang (finished=%d)", finished);
    AIC_CHECK_MSG(WIFSIGNALED(st) && WTERMSIG(st) == SIGKILL,
                  "child_hang: expected SIGKILL from the watchdog, got signaled=%d sig=%d",
                  WIFSIGNALED(st), WIFSIGNALED(st) ? WTERMSIG(st) : -1);
}

/* 4: a real SIGSEGV is reported as SIGSEGV, NOT laundered into SIGABRT. */
static void test_segv_not_abort(void)
{
    printf("4: child_segv -> SIGSEGV (not mistaken for a fail-loud abort)\n");
    int st = 0;
    char err[4096];
    int finished = aic_watchdog_run(child_segv, 20, 0, &st, err, sizeof err);
    AIC_CHECK_MSG(finished == 1, "child_segv: run() reported a hang (finished=0)");
    AIC_CHECK_MSG(WIFSIGNALED(st) && WTERMSIG(st) == SIGSEGV,
                  "child_segv: expected SIGSEGV, got signaled=%d sig=%d "
                  "(a crash must not be read as a SIGABRT fail-loud)",
                  WIFSIGNALED(st), WIFSIGNALED(st) ? WTERMSIG(st) : -1);
}

/* 5: capture_stdout=1 redirects the child's stdout into the captured buffer. */
static void test_stdout_capture(void)
{
    printf("5: child_stdout -> stdout marker captured (capture_stdout=1)\n");
    int st = 0;
    char err[4096];
    int finished = aic_watchdog_run(child_stdout, 20, 1, &st, err, sizeof err);
    AIC_CHECK_MSG(finished == 1, "child_stdout: run() reported a hang (finished=0)");
    AIC_CHECK_MSG(strstr(err, "wd_stdout_marker") != NULL,
                  "child_stdout: stdout marker NOT captured under capture_stdout=1. got:\n%s",
                  err);
}

int main(void)
{
    test_abort_failloud();
    test_clean_exit();
    test_hang_detected();
    test_segv_not_abort();
    test_stdout_capture();
    aic_test_report("test_watchdog");
    return 0;
}
