/* aic_watchdog.h — the canonical fork + SIGALRM fail-loud watchdog (bead aic-8de).
 *
 * Declarations only. The ONLY system header needed here is <stddef.h> for
 * size_t; the fork/alarm/sigaction/mkstemp/dup2/waitpid machinery (and the
 * _POSIX_C_SOURCE feature-test macro that exposes it) lives entirely in
 * aic_watchdog.c, so a consumer that only wants the fail-loud check needs
 * nothing more than this header — no feature-test macro, no <sys/wait.h>.
 *
 * WHY this exists: CLAUDE.md Rule 4 (fail fast, fail loud). Many routines must
 * abort() / flint_abort() on a bad input (out-of-basin spectrum, a non-converged
 * iteration, a violated eps-axiom). Asserting "this call aborts" cannot be done
 * in-process — abort() would take the test driver down with it. So we run the
 * routine in a forked child, observe its wait-status from the parent, and require
 * WIFSIGNALED && WTERMSIG == SIGABRT. A SIGALRM watchdog kills (and reports) a
 * child that HANGS instead of aborting, so a non-converging routine fails the
 * test loudly rather than hanging the gate forever.
 *
 * This unifies the nine near-identical private copies that had accreted across
 * the suite (fc_watch/xo0_run_child/v5f_run_child/EV_/EM_/KA_/UCP_/CHAN_ ...).
 */
#ifndef AIC_WATCHDOG_H
#define AIC_WATCHDOG_H

#include <stddef.h>

typedef void (*aic_watchdog_fn)(void);

/* Run fn() in a forked child under a SIGALRM watchdog of timeout_s seconds.
 * Child stderr (and stdout too, if capture_stdout != 0) is redirected to a
 * mkstemp temp file, read back NUL-terminated into err[0..errsz-1].
 * Writes the child wait-status to *status_out.
 * Returns 1 if the child finished within the watchdog, 0 if it was killed for
 * exceeding timeout_s (a hang). */
int aic_watchdog_run(aic_watchdog_fn fn, int timeout_s, int capture_stdout,
                     int *status_out, char *err, size_t errsz);

/* Assert fn() FAILS LOUD: finishes within timeout_s, did NOT exit cleanly
 * (exit 0 == silent success == wrong), died by SIGABRT, and — if needle != NULL
 * — its captured stderr contains needle. `who` is a label for diagnostics. */
void aic_watchdog_assert_failloud(aic_watchdog_fn fn, int timeout_s,
                                  const char *who, const char *needle);

#endif /* AIC_WATCHDOG_H */
