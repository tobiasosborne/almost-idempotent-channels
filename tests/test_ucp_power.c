/* test_ucp_power.c — cross-checks for channel power Phi^k + Choi-rank
 * compression (bead aic-xxk). These are THIN wrappers over compose / Choi
 * round-trip, so the tests assert the wrapper INVARIANTS (CLAUDE.md Rule 5/6),
 * not the underlying linear algebra (already covered by test_ucp_d24/test_ucp):
 *
 *   1. power(Phi,1) == Phi (Choi diff ~ 0); power(Phi,2) == compose(Phi,Phi).
 *   2. eta=0 idempotence oracle (the cleanest, Rule 6): for EXACTLY idempotent
 *      Phi (dephasing, cond_expectation, depolarizing(d,1.0)), power(Phi,k) == Phi
 *      for k = 2,3,5 — Choi diff ~ 0.
 *   3. power(Phi,0) == identity channel.
 *   4. Kraus-count BOUND (proves compress-each-step works): r(power(Phi,5)) <=
 *      dim_H^2, NOT phi->r^5 — for dephasing(d) (r=d) the naive Phi^5 is d^5 ops.
 *   5. compress invariants: Choi preserved (~0); r drops to the Choi rank; idem-
 *      potent in rank (compress twice == compress once); a deliberately non-
 *      minimal channel (compose(Phi,Phi), r^2 ops) drops to the Choi rank.
 *   6. NON-idempotent tooth: power(Phi,2) != Phi for depolarizing(3,0.3) (Choi
 *      diff is O(1)) — the power genuinely composes.
 *   7. MUTATION proof (Rule 7): documented below; a manual edit drives a test RED.
 *
 * Convention (aic_ucp.h): Kraus K_a are dim_K x dim_H, Phi(X)=sum K_a^dag X K_a
 * (Heisenberg); compose(out,phi,psi) = Phi o Psi, Kraus {psi.K[b] @ phi.K[a]}.
 * The exactly-idempotent channels come from aic_channels.h (bead aic-7hg). Kraus
 * reps are gauge-non-unique, so EVERY comparison is AS CHANNELS (Choi diff opnorm).
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_channels.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"

#define PREC 53

/* ||Choi(a) - Choi(b)||_op as an arb ball — the AS-CHANNELS equality predicate.
 * Both maps must share (dim_K, dim_H). Returns the opnorm into `out`. */
static void choi_diff_opnorm(arb_t out, const aic_ucp_kraus *a,
                             const aic_ucp_kraus *b)
{
    AIC_CHECK_MSG(a->dim_K == b->dim_K && a->dim_H == b->dim_H,
                  "choi_diff_opnorm: maps must share (dim_K, dim_H)");
    slong n = a->dim_K * a->dim_H;
    acb_mat_t L;
    acb_mat_init(L, n, n);
    aic_ucp_choi_diff(L, a, b, PREC);   /* Choi(a) - Choi(b) */
    aic_mat_opnorm(out, L, PREC);
    acb_mat_clear(L);
}

/* assert ||Choi(a) - Choi(b)||_op ~ 0 (same channel) */
static void check_same_channel(const aic_ucp_kraus *a, const aic_ucp_kraus *b,
                               const char *what)
{
    arb_t nrm, tol;
    arb_init(nrm);
    arb_init(tol);
    arb_set_d(tol, 1e-12);
    choi_diff_opnorm(nrm, a, b);
    AIC_CHECK_MSG(arb_le(nrm, tol), "%s: Choi diff not ~0 (different channel)",
                  what);
    arb_clear(tol);
    arb_clear(nrm);
}

/* ---- test 1: power(Phi,1)==Phi and power(Phi,2)==compose(Phi,Phi) ---- */
static void test_power_basics(void)
{
    aic_ucp_kraus phi;
    aic_channel_dephasing(&phi, 3, PREC);    /* a generic self-map */

    /* power(Phi,1) is a copy of phi (same channel). */
    aic_ucp_kraus p1;
    aic_ucp_power(&p1, &phi, 1, PREC);
    AIC_CHECK_MSG(p1.r == phi.r, "power(Phi,1) should copy phi's Kraus count");
    check_same_channel(&p1, &phi, "power(Phi,1)");
    aic_ucp_kraus_clear(&p1);

    /* power(Phi,2) == compose(Phi,Phi) AS CHANNELS. */
    aic_ucp_kraus p2, c2;
    aic_ucp_power(&p2, &phi, 2, PREC);
    aic_ucp_compose(&c2, &phi, &phi, PREC);
    check_same_channel(&p2, &c2, "power(Phi,2) vs compose(Phi,Phi)");
    aic_ucp_kraus_clear(&c2);
    aic_ucp_kraus_clear(&p2);

    aic_ucp_kraus_clear(&phi);
}

/* ---- test 2: eta=0 idempotence oracle — power(Phi,k)==Phi (Rule 6) ---- */
static void check_idemp_power(aic_ucp_kraus *phi, const char *name)
{
    slong ks[] = {2, 3, 5};
    for (int t = 0; t < 3; t++) {
        aic_ucp_kraus pk;
        aic_ucp_power(&pk, phi, ks[t], PREC);
        char what[96];
        snprintf(what, sizeof what, "%s: power(Phi,%ld)==Phi (eta=0 oracle)",
                 name, (long) ks[t]);
        check_same_channel(&pk, phi, what);
        aic_ucp_kraus_clear(&pk);
    }
}

static void test_power_idempotent_oracle(void)
{
    aic_ucp_kraus phi;

    aic_channel_dephasing(&phi, 4, PREC);
    check_idemp_power(&phi, "dephasing(4)");
    aic_ucp_kraus_clear(&phi);

    slong blocks[] = {2, 1, 3};
    aic_channel_cond_expectation(&phi, blocks, 3, PREC);
    check_idemp_power(&phi, "cond_exp(2,1,3)");
    aic_ucp_kraus_clear(&phi);

    aic_channel_depolarizing(&phi, 3, 1.0, PREC);   /* p=1 -> idempotent */
    check_idemp_power(&phi, "depolarizing(3,1.0)");
    aic_ucp_kraus_clear(&phi);
}

/* ---- test 3: power(Phi,0) == identity channel ---- */
static void test_power_zero_is_identity(void)
{
    aic_ucp_kraus phi;
    aic_channel_dephasing(&phi, 3, PREC);

    aic_ucp_kraus p0;
    aic_ucp_power(&p0, &phi, 0, PREC);
    AIC_CHECK_MSG(p0.dim_K == 3 && p0.dim_H == 3 && p0.r == 1,
                  "power(Phi,0) shape should be 3x3 r=1");

    /* hand-built identity channel: one Kraus op K_0 = 1_3. */
    aic_ucp_kraus id;
    aic_ucp_kraus_init(&id, 3, 3, 1);
    acb_mat_one(id.K[0]);
    check_same_channel(&p0, &id, "power(Phi,0) vs identity");

    /* and power(Phi,0) applied to a matrix is the identity: Phi^0(X) = X. */
    acb_mat_t X, p0X;
    acb_mat_init(X, 3, 3);
    acb_mat_init(p0X, 3, 3);
    for (slong i = 0; i < 3; i++)
        for (slong j = 0; j < 3; j++)
            acb_set_d_d(acb_mat_entry(X, i, j), 0.3 * (i + 1), -0.2 * (j + 1));
    aic_ucp_apply(p0X, &p0, X, PREC);
    arb_t nrm, tol;
    arb_init(nrm);
    arb_init(tol);
    arb_set_d(tol, 1e-12);
    acb_mat_sub(p0X, p0X, X, PREC);
    aic_mat_opnorm(nrm, p0X, PREC);
    AIC_CHECK_MSG(arb_le(nrm, tol), "power(Phi,0)(X) != X");
    arb_clear(tol);
    arb_clear(nrm);

    acb_mat_clear(p0X);
    acb_mat_clear(X);
    aic_ucp_kraus_clear(&id);
    aic_ucp_kraus_clear(&p0);
    aic_ucp_kraus_clear(&phi);
}

/* ---- test 4: Kraus-count bound — compress-each-step keeps r <= dim_H^2 ---- */
static void test_power_kraus_count_bound(void)
{
    /* dephasing(d): r = d, Choi rank = d (d independent |e_i><e_i| Kraus ops).
     * Naive Phi^5 would be d^5 = 4^5 = 1024 ops; compress-each-step keeps it at
     * the Choi rank d = 4 (<= dim_H^2 = 16). Assert the strong bound r == d. */
    const slong d = 4;
    aic_ucp_kraus phi;
    aic_channel_dephasing(&phi, d, PREC);

    aic_ucp_kraus p5;
    aic_ucp_power(&p5, &phi, 5, PREC);
    slong naive = 1;
    for (int i = 0; i < 5; i++) naive *= phi.r;   /* d^5 = 1024 */
    AIC_CHECK_MSG(p5.r <= d * d,
                  "power(dephasing(%ld),5) r=%ld exceeds dim_H^2=%ld "
                  "(compress-each-step failed; naive would be %ld)",
                  (long) d, (long) p5.r, (long) (d * d), (long) naive);
    AIC_CHECK_MSG(p5.r == d,
                  "power(dephasing(%ld),5) r=%ld should equal the Choi rank %ld",
                  (long) d, (long) p5.r, (long) d);
    AIC_CHECK_MSG(p5.r < naive, "compress-each-step did not shrink r below d^5");
    check_same_channel(&p5, &phi, "power(dephasing,5) channel");

    aic_ucp_kraus_clear(&p5);
    aic_ucp_kraus_clear(&phi);
}

/* ---- test 5: compress invariants ---- */
static void test_compress_invariants(void)
{
    /* depolarizing(3,0.3) has r = 1 + d^2 = 10 ops but Choi rank < 10 generally;
     * compress preserves the channel and drops r to the Choi rank. */
    aic_ucp_kraus phi;
    aic_channel_depolarizing(&phi, 3, 0.3, PREC);

    aic_ucp_kraus comp;
    aic_ucp_compress(&comp, &phi, PREC);
    check_same_channel(&comp, &phi, "compress preserves channel");
    AIC_CHECK_MSG(comp.r <= phi.r, "compress should not increase r (%ld -> %ld)",
                  (long) phi.r, (long) comp.r);
    AIC_CHECK_MSG(comp.r <= phi.dim_K * phi.dim_H,
                  "compress r=%ld exceeds Choi-rank bound dim_K*dim_H=%ld",
                  (long) comp.r, (long) (phi.dim_K * phi.dim_H));

    /* idempotent in rank: compress twice == compress once (same r, same chan). */
    aic_ucp_kraus comp2;
    aic_ucp_compress(&comp2, &comp, PREC);
    AIC_CHECK_MSG(comp2.r == comp.r,
                  "compress not rank-idempotent (%ld -> %ld)",
                  (long) comp.r, (long) comp2.r);
    check_same_channel(&comp2, &comp, "compress idempotent (channel)");
    aic_ucp_kraus_clear(&comp2);

    /* deliberately NON-minimal channel: compose(Phi,Phi) has r^2 = 100 ops;
     * compress drops it to the Choi rank of Phi^2 (<= dim_H^2 = 9). */
    aic_ucp_kraus phi2, phi2c;
    aic_ucp_compose(&phi2, &phi, &phi, PREC);
    AIC_CHECK_MSG(phi2.r == phi.r * phi.r, "compose(Phi,Phi) should have r^2 ops");
    aic_ucp_compress(&phi2c, &phi2, PREC);
    AIC_CHECK_MSG(phi2c.r < phi2.r,
                  "compress did not shrink the non-minimal compose (%ld -> %ld)",
                  (long) phi2.r, (long) phi2c.r);
    AIC_CHECK_MSG(phi2c.r <= phi.dim_K * phi.dim_H,
                  "compressed Phi^2 r=%ld exceeds dim_K*dim_H=%ld",
                  (long) phi2c.r, (long) (phi.dim_K * phi.dim_H));
    check_same_channel(&phi2c, &phi2, "compress preserves Phi^2 channel");
    aic_ucp_kraus_clear(&phi2c);
    aic_ucp_kraus_clear(&phi2);

    aic_ucp_kraus_clear(&comp);
    aic_ucp_kraus_clear(&phi);
}

/* ---- test 6: NON-idempotent tooth — power(Phi,2) != Phi ---- */
static void test_power_nonidempotent_tooth(void)
{
    /* depolarizing(3,0.3): Phi^2 = Phi_{p(2-p)} = Phi_{0.51} != Phi_{0.3}, so the
     * power genuinely composes and Choi(Phi^2)-Choi(Phi) is O(1). A tooth that
     * proves power(Phi,2) is NOT a no-op for a non-idempotent channel. */
    aic_ucp_kraus phi;
    aic_channel_depolarizing(&phi, 3, 0.3, PREC);

    aic_ucp_kraus p2;
    aic_ucp_power(&p2, &phi, 2, PREC);

    arb_t nrm, big;
    arb_init(nrm);
    arb_init(big);
    arb_set_d(big, 0.01);
    choi_diff_opnorm(nrm, &p2, &phi);
    AIC_CHECK_MSG(arb_gt(nrm, big),
                  "power(depolarizing(3,0.3),2) should DIFFER from Phi (O(1) "
                  "Choi diff) — power is not composing");
    arb_clear(big);
    arb_clear(nrm);

    aic_ucp_kraus_clear(&p2);
    aic_ucp_kraus_clear(&phi);
}

/* MUTATION PROOF (Rule 7), done as a documented manual procedure since the new
 * src/test files are UNTRACKED (git checkout cannot restore them):
 *
 *   MUT-A (power skips the last compose). In src/aic_ucp_power.c change the loop
 *     bound `step <= k` to `step < k`. Rebuild. test_power_idempotent_oracle still
 *     passes (Phi^k == Phi for idempotent Phi regardless of one fewer compose),
 *     but test_power_nonidempotent_tooth FLIPS: power(Phi,2) now equals Phi (only
 *     the k=2 start, no compose), so the "should DIFFER" check goes RED. Restore.
 *
 *   MUT-B (compress keeps all ops). In aic_ucp_power's loop, replace the
 *     `compress(composed)` with a copy of `composed` (skip compression). Rebuild.
 *     The CHANNEL checks still pass, but test_power_kraus_count_bound goes RED:
 *     r(power(dephasing(4),5)) becomes 4^5 = 1024, blowing the `r <= dim_H^2`
 *     and `r == d` assertions.
 *
 * Both were exercised during development (see the REPORT). This block documents
 * the procedure so a future reader can re-run the mutation proof by hand. */

int main(void)
{
    test_power_basics();
    test_power_idempotent_oracle();
    test_power_zero_is_identity();
    test_power_kraus_count_bound();
    test_compress_invariants();
    test_power_nonidempotent_tooth();
    aic_test_report("test_ucp_power");
    printf("OK test_ucp_power\n");
    return 0;
}
