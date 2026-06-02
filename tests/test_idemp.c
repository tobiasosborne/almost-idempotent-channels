/* test_idemp.c — cross-checks for the idemp_structure module (bead aic-wuh), the
 * eta=0 oracle (milestone aic-9kk). th_idemp_structure: an exactly-idempotent UCP
 * Phi -> carrier M, C* algebra A = Img Phi, UCP maps Delta, Gamma with the five
 * defining relations (approximate_algebras.tex:2080-2090) holding to ZERO defect.
 *
 * Oracle channels (constructors in tests/test_idemp.h) with KNOWN (dim_M, dim_A);
 * for each we assert the dims AND the five relations to ~1e-12:
 *   block_cond_exp(5,2): (5, 13)   trace_replace(3): (3, 1)
 *   noiseless_subsystem(2,2): (4,4) identity(3): (3, 9)
 *   dephasing(4): (4, 4)           compress_idemp(4,2) [M PROPER]: (2, 4)
 *
 * GAUGE FREEDOM: the A-basis {B_k} is unique only up to a unitary on C^{dim_A}; we
 * compare the SUBSPACE via Pi_A = Delta Delta^dag, never basis vectors. The five
 * relations are gauge-invariant and asserted directly.
 *
 * Cross-check ladder: internal sanity (Pi_M^2=Pi_M, J_M^dag J_M=1_M, Delta cols
 * orthonormal); acb@53 vs acb@256 (relation defects + Pi_A) — the extraction is
 * double-path regardless of prec, so this is the precision-self-consistency rung;
 * eta=0 oracle (the five relations == 0); dim counts vs known. Mutation proofs:
 * source mutations (zero Lambda, scale Gamma, wrong J_M block) were run by hand
 * and each turned a relation RED (see report); the in-suite teeth test confirms a
 * non-idempotent map aborts via SIGABRT at the entry guard. The Gamma-CP cert is
 * the PAIR (Choi(C_M o Lambda) PSD + ||w Gamma - C_M Lambda|| = 0); the latter has
 * teeth (a scaled/zeroed Gamma makes it O(1)) — the bare Choi check does not read
 * Gamma and is vacuous on its own (FIX 1, bead aic-wuh review).
 */
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "aic_watchdog.h"
#include "test_idemp.h"

static void test_oracle_channels(void)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 5, 2);  check_decomp("block_cond_exp(5,2)", &phi, 5, 13);
    aic_ucp_kraus_clear(&phi);
    make_trace_replace(&phi, 3);      check_decomp("trace_replace(3)", &phi, 3, 1);
    aic_ucp_kraus_clear(&phi);
    make_noiseless_subsystem(&phi, 2, 2);
    check_decomp("noiseless_subsystem(2,2)", &phi, 4, 4);
    aic_ucp_kraus_clear(&phi);
    make_identity(&phi, 3);           check_decomp("identity(3)", &phi, 3, 9);
    aic_ucp_kraus_clear(&phi);
    make_dephasing(&phi, 4);          check_decomp("dephasing(4)", &phi, 4, 4);
    aic_ucp_kraus_clear(&phi);
    make_compress_idemp(&phi, 4, 2);  /* PROPER carrier M < H, exercises M^perp */
    check_decomp("compress_idemp(4,2) [M proper]", &phi, 2, 4);
    aic_ucp_kraus_clear(&phi);
}

/* acb@53 vs acb@256: the defining relations and the gauge-invariant Pi_A must
 * agree across precisions (the extraction is double-path regardless of prec, so
 * the subspaces are identical; this certifies the certified-path arithmetic). */
static void test_prec_consistency(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);

    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 5, 2);

    aic_idemp_decomp d53, d256;
    aic_idemp_decompose(&d53, &phi, 53);
    aic_idemp_decompose(&d256, &phi, 256);
    AIC_CHECK_MSG(d53.dim_A == d256.dim_A && d53.dim_M == d256.dim_M,
                  "prec_consistency: dims differ across prec");

    /* Pi_A agrees (gauge-invariant subspace) */
    slong nn = d53.n * d53.n;
    acb_mat_t PA53, PA256;
    acb_mat_init(PA53, nn, nn);
    acb_mat_init(PA256, nn, nn);
    build_Pi_A(PA53, d53.Delta, 256);
    build_Pi_A(PA256, d256.Delta, 256);
    AIC_CHECK_ACB_MAT_CLOSE(PA53, PA256, tol);
    acb_mat_clear(PA256);
    acb_mat_clear(PA53);

    /* relation defects agree (both ~0) */
    arb_t g53, g256;
    arb_init(g53); arb_init(g256);
    aic_idemp_defect_GCD(g53, &d53, 53);
    aic_idemp_defect_GCD(g256, &d256, 256);
    AIC_CHECK(arb_le(g53, tol) && arb_le(g256, tol));
    arb_clear(g256); arb_clear(g53);

    aic_idemp_clear(&d256);
    aic_idemp_clear(&d53);
    aic_ucp_kraus_clear(&phi);
    arb_clear(tol);
}

/* Mutation teeth child: a UCP map that is NOT idempotent: a simple non-idempotent
 * unital map Phi(X) = 0.5 X + 0.5 diag(X) (a partial dephasing). Kraus:
 * K_0 = sqrt(0.5) 1_2, K_a = sqrt(0.5)|e_a><e_a|. Then Phi^2 != Phi, so the entry
 * guard must flint_abort(). No args (the shared aic_watchdog child is void(void)). */
static void idemp_child_non_idempotent(void)
{
    aic_ucp_kraus phi;
    aic_ucp_kraus_init(&phi, 2, 2, 3);
    double s = sqrt(0.5);
    acb_set_d(acb_mat_entry(phi.K[0], 0, 0), s);
    acb_set_d(acb_mat_entry(phi.K[0], 1, 1), s);   /* sqrt(0.5) 1_2 */
    acb_set_d(acb_mat_entry(phi.K[1], 0, 0), s);   /* sqrt(0.5)|0><0| */
    acb_set_d(acb_mat_entry(phi.K[2], 1, 1), s);   /* sqrt(0.5)|1><1| */
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);           /* must flint_abort() */
}

/* Mutation teeth: a non-idempotent UCP map must ABORT at the entry guard. This is
 * fail-loud (Rule 4); we run it under the shared aic_watchdog (bead aic-8de),
 * which ADDS a SIGALRM hang-backstop the old bare fork LACKED (it could hang the
 * gate forever). assert_failloud keeps the SIGABRT requirement (not a SEGV/OOM
 * masquerade) and adds the silent-exit-0 check.
 *
 * needle = NULL: aic_idemp_decompose aborts via flint_printf + flint_abort, and
 * flint_printf writes to STDOUT (verified), which assert_failloud (capture_stdout=0)
 * does not capture — a stdout-bound message can never match a stderr needle. The
 * "NOT exactly idempotent" string therefore lands on the console, not the captured
 * buffer; NULL avoids a spurious miss. */
static void test_entry_guard(void)
{
    aic_watchdog_assert_failloud(idemp_child_non_idempotent, 20,
                                 "idemp/entry-guard", "NOT exactly idempotent");
}

int main(void)
{
    test_oracle_channels();
    test_prec_consistency();
    test_entry_guard();
    aic_test_report("test_idemp");
    printf("OK test_idemp\n");
    return 0;
}
