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
 * THE WATCHDOG. The shared aic_watchdog (bead aic-8de) forks; the child redirects
 * stderr to a temp file and runs the pipeline call. The parent waits with a
 * SIGALRM-bounded waitpid (AIC_XO0_WATCH_S seconds): if the child has not
 * terminated by then it is SIGKILLed and the test FAILS (a hang IS a failure,
 * Rule 4). The shared helper's child signature is void(void); this driver needs a
 * UCP self-map per child, so a module-static `xo0_ctx` pointer is set immediately
 * before each call and the void(void) wrapper `child_regularize` reads it.
 */
#include <math.h>
#include <stdio.h>
#include <sys/wait.h>            /* WIFEXITED/WEXITSTATUS for the clean-run checks */

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_funcalc.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "aic_watchdog.h"
#include "test_idemp.h"            /* make_compress_idemp, make_dephasing, PREC */

#define XO0_PREC 256
#define AIC_XO0_WATCH_S 20         /* watchdog: a child alive past this is a HANG */

/* ---- the child payload (ctx-bridge to the void(void) shared watchdog) ------ */

/* The shared aic_watchdog child takes NO args; this driver needs a per-call UCP
 * self-map. Bridge via a module-static pointer set immediately before each run. */
static const aic_ucp_kraus *xo0_ctx = NULL;

/* xo0_ctx = a UCP self-map; the child regularizes it (the pipeline entry). */
static void child_regularize(void)
{
    const aic_ucp_kraus *phi = xo0_ctx;
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

    /* assert_failloud asserts: finished (not a hang), NOT a clean exit-0 (no
     * silent regularization of a non-convergent series), died by SIGABRT, and the
     * captured stderr names the eta<1/4 almost-idempotent precondition. The needle
     * "eta < 1/4" is the load-bearing entry-guard string (.tex:2168/2176); the
     * same flint_printf also names "aic_assoc_regularize" (the API-boundary
     * provenance), so a single needle keeps both teeth. RED before the fix (the
     * deep prop_P message names neither). */
    xo0_ctx = &phi;
    aic_watchdog_assert_failloud(child_regularize, AIC_XO0_WATCH_S,
                                 "aic_assoc_regularize/out-of-basin", "eta < 1/4");
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

    /* (a) eta=0 exact idempotent. A clean-exit case: aic_watchdog_run + inline
     * "finished && exit 0" check (assert_failloud is for the abort cases only). */
    {
        aic_ucp_kraus phi;
        make_compress_idemp(&phi, 4, 2);
        char err[4096];
        int st = 0;
        xo0_ctx = &phi;
        int finished = aic_watchdog_run(child_regularize, AIC_XO0_WATCH_S, 0,
                                        &st, err, sizeof err);
        AIC_CHECK_MSG(finished, "eta=0 compress_idemp(4,2) HUNG in regularize");
        AIC_CHECK_MSG(WIFEXITED(st) && WEXITSTATUS(st) == 0,
                      "eta=0 compress_idemp(4,2): expected a clean in-regime run, "
                      "got FALSE-rejection by the entry guard. stderr:\n%s", err);
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
        int st = 0;
        xo0_ctx = &phi;
        int finished = aic_watchdog_run(child_regularize, AIC_XO0_WATCH_S, 0,
                                        &st, err, sizeof err);
        AIC_CHECK_MSG(finished, "spectral-relaxation mix(4,1)+deph(4) HUNG in regularize");
        AIC_CHECK_MSG(WIFEXITED(st) && WEXITSTATUS(st) == 0,
                      "spectral-relaxation channel (||S^2-S||_op>=1/4 but rho<1/4): "
                      "expected a clean in-regime run, got FALSE-rejection by an "
                      "over-aggressive (op-norm) entry guard — the guard must use "
                      "the SPECTRAL rho (aic-8hz). stderr:\n%s", err);
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
