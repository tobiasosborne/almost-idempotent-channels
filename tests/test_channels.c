/* test_channels.c — invariant TDD for the public channel constructors (bead
 * aic-7hg, module "channels"). CLAUDE.md Rule 5 ("runs without errors" is NOT a
 * pass): every channel is checked against UNITAL + CP + a channel-specific
 * invariant, with the eta=0 idempotence oracle (Rule 6, cross-check #3) for the
 * exactly-idempotent ones, and cross-checked against the static test-fixture
 * oracles (make_dephasing / make_block_cond_exp) at the Choi level.
 *
 * The fail-loud input validation (non-unitary U, unnormalised probs) is exercised
 * via the fork+watchdog subprocess pattern from test_xo0_failloud.c so a Rule-4
 * abort is testable without crashing the test binary.
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

#include "aic/aic_channels.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "test_idemp.h"   /* make_dephasing, make_block_cond_exp, PREC, set_tol */

/* ---- shared invariant checks --------------------------------------------- */

/* UNITAL: ||sum_a K_a^dag K_a - 1||_op ~ 0 (machine). */
static void check_unital(const aic_ucp_kraus *phi, const char *name)
{
    arb_t def, tol;
    arb_init(def);
    arb_init(tol);
    set_tol(tol, 1e-12);
    aic_ucp_unital_defect_kraus(def, phi, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "%s: not unital (defect > 1e-12)", name);
    arb_clear(tol);
    arb_clear(def);
}

/* CP: the Choi matrix is PSD (aic_ucp_is_cp_choi == 1). */
static void check_cp(const aic_ucp_kraus *phi, const char *name)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);
    slong n = phi->dim_K * phi->dim_H;
    acb_mat_t C;
    acb_mat_init(C, n, n);
    aic_ucp_kraus_to_choi(C, phi, PREC);
    AIC_CHECK_MSG(aic_ucp_is_cp_choi(C, tol, PREC) == 1,
                  "%s: Choi matrix not certified PSD (not CP)", name);
    acb_mat_clear(C);
    arb_clear(tol);
}

/* eta=0 oracle: ||Choi(Phi^2 - Phi)||_op ~ 0 (Phi exactly idempotent). */
static void check_idempotent(aic_ucp_kraus *phi, const char *name)
{
    arb_t tol, nrm;
    arb_init(tol);
    arb_init(nrm);
    set_tol(tol, 1e-12);
    slong n = phi->dim_K * phi->dim_H;
    aic_ucp_kraus phi2;
    aic_ucp_compose(&phi2, phi, phi, PREC);
    acb_mat_t L;
    acb_mat_init(L, n, n);
    aic_ucp_choi_diff(L, &phi2, phi, PREC);
    aic_mat_opnorm(nrm, L, PREC);
    AIC_CHECK_MSG(arb_le(nrm, tol),
                  "%s: Choi(Phi^2-Phi) not zero (eta=0 idempotence oracle)", name);
    acb_mat_clear(L);
    aic_ucp_kraus_clear(&phi2);
    arb_clear(nrm);
    arb_clear(tol);
}

/* Choi-equality cross-check vs a static oracle channel (both d x d self-maps). */
static void check_choi_eq(const aic_ucp_kraus *a, const aic_ucp_kraus *b,
                          const char *name)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);
    slong n = a->dim_K * a->dim_H;
    AIC_CHECK_MSG(b->dim_K * b->dim_H == n, "%s: oracle shape mismatch", name);
    acb_mat_t Ca, Cb, D;
    acb_mat_init(Ca, n, n);
    acb_mat_init(Cb, n, n);
    acb_mat_init(D, n, n);
    aic_ucp_kraus_to_choi(Ca, a, PREC);
    aic_ucp_kraus_to_choi(Cb, b, PREC);
    acb_mat_sub(D, Ca, Cb, PREC);
    arb_t nrm;
    arb_init(nrm);
    aic_mat_opnorm(nrm, D, PREC);
    AIC_CHECK_MSG(arb_le(nrm, tol), "%s: Choi mismatch vs static oracle", name);
    arb_clear(nrm);
    acb_mat_clear(D);
    acb_mat_clear(Cb);
    acb_mat_clear(Ca);
    arb_clear(tol);
}

/* ---- the four Pauli unitaries (for the closed-group twirl) ---------------- */

static void set_pauli_unitary(acb_mat_t U, int k)
{
    acb_mat_zero(U);
    switch (k) {
    case 0: acb_set_si(acb_mat_entry(U, 0, 0), 1);
            acb_set_si(acb_mat_entry(U, 1, 1), 1); break;
    case 1: acb_set_si(acb_mat_entry(U, 0, 1), 1);
            acb_set_si(acb_mat_entry(U, 1, 0), 1); break;
    case 2: acb_set_d_d(acb_mat_entry(U, 0, 1), 0.0, -1.0);
            acb_set_d_d(acb_mat_entry(U, 1, 0), 0.0,  1.0); break;
    case 3: acb_set_si(acb_mat_entry(U, 0, 0), 1);
            acb_set_si(acb_mat_entry(U, 1, 1), -1); break;
    }
}

/* ---- per-channel tests ---------------------------------------------------- */

static void test_dephasing(void)
{
    printf("dephasing(4): unital, CP, idempotent, Choi == make_dephasing\n");
    aic_ucp_kraus phi, oracle;
    aic_channel_dephasing(&phi, 4, PREC);
    make_dephasing(&oracle, 4);
    check_unital(&phi, "dephasing");
    check_cp(&phi, "dephasing");
    check_idempotent(&phi, "dephasing");
    check_choi_eq(&phi, &oracle, "dephasing vs make_dephasing");
    aic_ucp_kraus_clear(&oracle);
    aic_ucp_kraus_clear(&phi);
}

static void test_cond_expectation(void)
{
    printf("cond_expectation([2,2]): unital, CP, idempotent, Choi == "
           "make_block_cond_exp(4,2)\n");
    aic_ucp_kraus phi, oracle;
    slong blocks[2] = {2, 2};
    aic_channel_cond_expectation(&phi, blocks, 2, PREC);
    make_block_cond_exp(&oracle, 4, 2);
    check_unital(&phi, "cond_expectation");
    check_cp(&phi, "cond_expectation");
    check_idempotent(&phi, "cond_expectation");
    check_choi_eq(&phi, &oracle, "cond_expectation vs make_block_cond_exp");
    aic_ucp_kraus_clear(&oracle);
    aic_ucp_kraus_clear(&phi);

    /* uneven 3-block partition of C^6 — generalises beyond the 2-block oracle. */
    aic_ucp_kraus phi3;
    slong b3[3] = {1, 2, 3};
    aic_channel_cond_expectation(&phi3, b3, 3, PREC);
    check_unital(&phi3, "cond_expectation[1,2,3]");
    check_cp(&phi3, "cond_expectation[1,2,3]");
    check_idempotent(&phi3, "cond_expectation[1,2,3]");
    aic_ucp_kraus_clear(&phi3);
}

/* depolarizing: idempotent at p=1; closed form at general p on |0><0|. */
static void test_depolarizing(void)
{
    printf("depolarizing(3): p=0 identity, p=1 idempotent, general-p closed "
           "form on |0><0|\n");
    /* p=1: fully depolarizing = conditional expectation onto C*1, idempotent. */
    aic_ucp_kraus phi1;
    aic_channel_depolarizing(&phi1, 3, 1.0, PREC);
    check_unital(&phi1, "depolarizing p=1");
    check_cp(&phi1, "depolarizing p=1");
    check_idempotent(&phi1, "depolarizing p=1");
    aic_ucp_kraus_clear(&phi1);

    /* p=0: identity, idempotent. */
    aic_ucp_kraus phi0;
    aic_channel_depolarizing(&phi0, 3, 0.0, PREC);
    check_unital(&phi0, "depolarizing p=0");
    check_idempotent(&phi0, "depolarizing p=0");
    aic_ucp_kraus_clear(&phi0);

    /* general p: Phi_p(|0><0|) = (1-p)|0><0| + p/d 1_d.  d=3, p=0.3:
     * diagonal (1-p)+p/d = 0.8 at [0,0], p/d = 0.1 at [1,1],[2,2], 0 elsewhere. */
    {
        double p = 0.3;
        slong d = 3;
        aic_ucp_kraus phi;
        aic_channel_depolarizing(&phi, d, p, PREC);
        check_unital(&phi, "depolarizing p=0.3");
        check_cp(&phi, "depolarizing p=0.3");
        acb_mat_t X, Y, Ref, Diff;
        acb_mat_init(X, d, d);
        acb_mat_init(Y, d, d);
        acb_mat_init(Ref, d, d);
        acb_mat_init(Diff, d, d);
        acb_set_si(acb_mat_entry(X, 0, 0), 1);          /* |0><0| */
        aic_ucp_apply(Y, &phi, X, PREC);
        acb_set_d(acb_mat_entry(Ref, 0, 0), (1.0 - p) + p / (double) d);
        acb_set_d(acb_mat_entry(Ref, 1, 1), p / (double) d);
        acb_set_d(acb_mat_entry(Ref, 2, 2), p / (double) d);
        acb_mat_sub(Diff, Y, Ref, PREC);
        arb_t nrm, tol;
        arb_init(nrm);
        arb_init(tol);
        set_tol(tol, 1e-12);
        aic_mat_opnorm(nrm, Diff, PREC);
        AIC_CHECK_MSG(arb_le(nrm, tol),
                      "depolarizing p=0.3: Phi_p(|0><0|) != closed form");
        arb_clear(tol);
        arb_clear(nrm);
        acb_mat_clear(Diff);
        acb_mat_clear(Ref);
        acb_mat_clear(Y);
        acb_mat_clear(X);
        aic_ucp_kraus_clear(&phi);
    }
}

/* pauli: unital + CP always; reduces to identity at (1,0,0,0) and to dephasing
 * at (1/2,0,0,1/2) — cross-checked against make_dephasing(2) at the Choi level. */
static void test_pauli(void)
{
    printf("pauli: unital+CP; (1,0,0,0)=identity idempotent; "
           "(1/2,0,0,1/2)=dephasing\n");
    double pmix[4] = {0.4, 0.3, 0.2, 0.1};
    aic_ucp_kraus phi;
    aic_channel_pauli(&phi, pmix, PREC);
    check_unital(&phi, "pauli mixed");
    check_cp(&phi, "pauli mixed");
    aic_ucp_kraus_clear(&phi);

    double pid[4] = {1.0, 0.0, 0.0, 0.0};
    aic_ucp_kraus phii;
    aic_channel_pauli(&phii, pid, PREC);
    check_idempotent(&phii, "pauli identity");
    aic_ucp_kraus_clear(&phii);

    /* (I + Z)/2 conjugation = qubit dephasing: K = (sigma_0 +- sigma_3)/2 gives
     * the diagonal projectors, so Phi == make_dephasing(2). */
    double pdeph[4] = {0.5, 0.0, 0.0, 0.5};
    aic_ucp_kraus phid, oracle;
    aic_channel_pauli(&phid, pdeph, PREC);
    make_dephasing(&oracle, 2);
    check_idempotent(&phid, "pauli dephasing");
    check_choi_eq(&phid, &oracle, "pauli(1/2,0,0,1/2) vs make_dephasing(2)");
    aic_ucp_kraus_clear(&oracle);
    aic_ucp_kraus_clear(&phid);
}

/* mix_unitaries / group_twirl: the Pauli-group twirl over {I,X,Y,Z} is the
 * maximally depolarizing channel (projects onto C*1) and is EXACTLY idempotent;
 * cross-check its Choi equals depolarizing(2, p=1). A non-group unitary mixture
 * (just {I, X} with weights 1/2) is unital+CP but NOT idempotent. */
static void test_mix_and_twirl(void)
{
    printf("mix_unitaries / group_twirl: Pauli-group twirl = depolarizing(2,1) "
           "(idempotent); non-group mix unital+CP\n");
    acb_mat_t U[4];
    for (int k = 0; k < 4; k++) {
        acb_mat_init(U[k], 2, 2);
        set_pauli_unitary(U[k], k);
    }
    /* closed-group twirl over the Pauli group. */
    aic_ucp_kraus tw, dep1;
    aic_channel_group_twirl(&tw, (const acb_mat_t *) U, 4, 2, PREC);
    aic_channel_depolarizing(&dep1, 2, 1.0, PREC);
    check_unital(&tw, "pauli twirl");
    check_cp(&tw, "pauli twirl");
    check_idempotent(&tw, "pauli twirl");
    check_choi_eq(&tw, &dep1, "pauli twirl vs depolarizing(2,1)");
    aic_ucp_kraus_clear(&dep1);
    aic_ucp_kraus_clear(&tw);

    /* non-idempotent mix {I, X} with UNEQUAL weights {0.7, 0.3}: unital+CP but
     * NOT a projection. (Equal weights {1/2,1/2} WOULD be idempotent — {I,X} is
     * a closed order-2 group and the half-half twirl is the projection onto the
     * X-symmetric part; unequal weights break that: Phi^2 = 0.58 id + 0.42 Ad_X
     * != 0.7 id + 0.3 Ad_X.) U[0],U[1] are contiguous, so &U[0] is a length-2
     * unitary array. */
    double w[2] = {0.7, 0.3};
    aic_ucp_kraus mx;
    aic_channel_mix_unitaries(&mx, (const acb_mat_t *) &U[0], w, 2, 2, PREC);
    check_unital(&mx, "mix {I,X}");
    check_cp(&mx, "mix {I,X}");
    /* assert NOT idempotent: ||Choi(Phi^2-Phi)|| is O(1), not ~0. */
    {
        arb_t nrm, tol;
        arb_init(nrm);
        arb_init(tol);
        set_tol(tol, 1e-6);
        slong n = 4;
        aic_ucp_kraus mx2;
        aic_ucp_compose(&mx2, &mx, &mx, PREC);
        acb_mat_t L;
        acb_mat_init(L, n, n);
        aic_ucp_choi_diff(L, &mx2, &mx, PREC);
        aic_mat_opnorm(nrm, L, PREC);
        AIC_CHECK_MSG(!arb_le(nrm, tol),
                      "mix {I,X} unexpectedly idempotent (should be O(1) defect)");
        acb_mat_clear(L);
        aic_ucp_kraus_clear(&mx2);
        arb_clear(tol);
        arb_clear(nrm);
    }
    aic_ucp_kraus_clear(&mx);

    for (int k = 0; k < 4; k++) acb_mat_clear(U[k]);
}

/* ---- fail-loud validation teeth (fork+watchdog, from test_xo0_failloud.c) -- */

#define CHAN_WATCH_S 15
static volatile pid_t chan_watch_pid = 0;
static volatile sig_atomic_t chan_timed_out = 0;
static void chan_alarm(int sig)
{
    (void) sig;
    chan_timed_out = 1;
    if (chan_watch_pid > 0) kill(chan_watch_pid, SIGKILL);
}

typedef void (*chan_child_fn)(void);
static int chan_run_child(chan_child_fn fn, int *status)
{
    fflush(NULL);
    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "chan_run_child: fork failed");
    if (pid == 0) {
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) { dup2(devnull, STDERR_FILENO); close(devnull); }
        fn();
        _exit(0);
    }
    chan_watch_pid = pid;
    chan_timed_out = 0;
    struct sigaction sa = {0}, old;
    sa.sa_handler = chan_alarm;
    sigaction(SIGALRM, &sa, &old);
    alarm(CHAN_WATCH_S);
    int st = 0;
    pid_t w;
    do { w = waitpid(pid, &st, 0); } while (w < 0 && !chan_timed_out);
    alarm(0);
    sigaction(SIGALRM, &old, NULL);
    if (w < 0) waitpid(pid, &st, 0);
    *status = st;
    return chan_timed_out ? 0 : 1;
}

/* child: probabilities that do NOT sum to 1 -> must abort. */
static void child_bad_probs(void)
{
    double bad[4] = {0.5, 0.5, 0.5, 0.0};   /* sums to 1.5 */
    aic_ucp_kraus phi;
    aic_channel_pauli(&phi, bad, PREC);     /* expected: flint_abort */
    aic_ucp_kraus_clear(&phi);
}

/* child: a non-unitary "U" in mix_unitaries -> must abort. */
static void child_non_unitary(void)
{
    acb_mat_t U;
    acb_mat_init(U, 2, 2);
    acb_set_si(acb_mat_entry(U, 0, 0), 2);  /* 2*I is not unitary */
    acb_set_si(acb_mat_entry(U, 1, 1), 2);
    double w[1] = {1.0};
    aic_ucp_kraus phi;
    /* &U is acb_mat_struct(*)[1]; the param is const acb_mat_t* — same target
     * type, the cast just drops the const-qualifier pedantic mismatch. */
    aic_channel_mix_unitaries(&phi, (const acb_mat_t *) &U, w, 1, 2, PREC);
    aic_ucp_kraus_clear(&phi);
    acb_mat_clear(U);
}

static void test_validation_failloud(void)
{
    printf("validation: bad-probs and non-unitary inputs must fail loud "
           "(SIGABRT, not hang)\n");
    int status = 0;
    int finished = chan_run_child(child_bad_probs, &status);
    AIC_CHECK_MSG(finished, "bad-probs validation HUNG (watchdog killed it)");
    AIC_CHECK_MSG(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT,
                  "bad-probs did NOT fail loud (expected SIGABRT; signaled=%d "
                  "exited=%d code=%d)", WIFSIGNALED(status), WIFEXITED(status),
                  WIFEXITED(status) ? WEXITSTATUS(status) : -1);

    finished = chan_run_child(child_non_unitary, &status);
    AIC_CHECK_MSG(finished, "non-unitary validation HUNG (watchdog killed it)");
    AIC_CHECK_MSG(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT,
                  "non-unitary did NOT fail loud (expected SIGABRT; signaled=%d "
                  "exited=%d code=%d)", WIFSIGNALED(status), WIFEXITED(status),
                  WIFEXITED(status) ? WEXITSTATUS(status) : -1);
    printf("  both invalid inputs SIGABRT (OK)\n");
}

/* The dag-order convention cross-check (hostile review #9, 2026-06-01): every
 * mix_unitaries/group_twirl test above uses HERMITIAN Pauli unitaries (U=U^dag),
 * so a dropped ^dag (the transpose bug aic_channels_unitary.c warns is load-
 * bearing) would pass them ALL. This pins the order with a GENUINELY non-Hermitian
 * U on a non-symmetric Hermitian X: the channel MUST give Phi(X)=U X U^dag, and
 * MUST differ from the transpose U^dag X U. */
static void test_dag_order_nonhermitian(void)
{
    printf("dag-order: mix_unitaries(m=1), NON-Hermitian U gives U X U^dag (not "
           "the transpose)\n");
    const slong d = 2;
    /* U = phase gate diag(1, e^{i 0.7}): unitary, NOT Hermitian. */
    acb_mat_t U;
    acb_mat_init(U, d, d);
    acb_mat_zero(U);
    acb_set_si(acb_mat_entry(U, 0, 0), 1);
    acb_set_d_d(acb_mat_entry(U, 1, 1), cos(0.7), sin(0.7));
    /* X = [[1, 2+i],[2-i, 3]]: Hermitian, non-diagonal, X != X^T. */
    acb_mat_t X;
    acb_mat_init(X, d, d);
    acb_set_si(acb_mat_entry(X, 0, 0), 1);
    acb_set_d_d(acb_mat_entry(X, 0, 1), 2.0, 1.0);
    acb_set_d_d(acb_mat_entry(X, 1, 0), 2.0, -1.0);
    acb_set_si(acb_mat_entry(X, 1, 1), 3);

    double w[1] = {1.0};
    aic_ucp_kraus phi;
    aic_channel_mix_unitaries(&phi, (const acb_mat_t *) &U, w, 1, d, PREC);

    acb_mat_t Y, Ud, T, UXUd, UdXU, Diff;
    acb_mat_init(Y, d, d);
    acb_mat_init(Ud, d, d);
    acb_mat_init(T, d, d);
    acb_mat_init(UXUd, d, d);
    acb_mat_init(UdXU, d, d);
    acb_mat_init(Diff, d, d);

    aic_ucp_apply(Y, &phi, X, PREC);            /* channel output */
    acb_mat_conjugate_transpose(Ud, U);         /* U^dag */
    acb_mat_mul(T, U, X, PREC);
    acb_mat_mul(UXUd, T, Ud, PREC);             /* U X U^dag (expected) */
    acb_mat_mul(T, Ud, X, PREC);
    acb_mat_mul(UdXU, T, U, PREC);              /* U^dag X U (the transpose) */

    arb_t nrm, tol, big;
    arb_init(nrm);
    arb_init(tol);
    arb_init(big);
    set_tol(tol, 1e-12);
    set_tol(big, 0.1);

    acb_mat_sub(Diff, Y, UXUd, PREC);
    aic_mat_opnorm(nrm, Diff, PREC);
    AIC_CHECK_MSG(arb_le(nrm, tol),
                  "dag-order: Phi(X) != U X U^dag — the dag-order is BACKWARDS");

    acb_mat_sub(Diff, Y, UdXU, PREC);
    aic_mat_opnorm(nrm, Diff, PREC);
    AIC_CHECK_MSG(arb_gt(nrm, big),
                  "dag-order: Phi(X) == U^dag X U — U/X too symmetric, tooth VACUOUS");
    printf("  Phi(X) = U X U^dag confirmed (transpose differs by O(1))\n");

    arb_clear(big);
    arb_clear(tol);
    arb_clear(nrm);
    acb_mat_clear(Diff);
    acb_mat_clear(UdXU);
    acb_mat_clear(UXUd);
    acb_mat_clear(T);
    acb_mat_clear(Ud);
    acb_mat_clear(Y);
    aic_ucp_kraus_clear(&phi);
    acb_mat_clear(X);
    acb_mat_clear(U);
}

int main(void)
{
    test_dephasing();
    test_cond_expectation();
    test_depolarizing();
    test_pauli();
    test_mix_and_twirl();
    test_dag_order_nonhermitian();
    test_validation_failloud();
    aic_test_report("test_channels");
    return 0;
}
