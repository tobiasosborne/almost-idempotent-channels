/* aic_assoc_regularize.c — the exact-idempotent regularization
 * Phi_tilde = theta(2 Phi - 1) (bead aic-92f, module "assoc_ecsa", Increment 1).
 * Arb path; see include/aic_assoc.h for the full why.
 *
 * approximate_algebras.tex:2171-2174 — equation tilde_Phi:
 *     Phi_tilde = theta(2 Phi - 1)
 *               = (1/2)(1 + sgn(2 Phi - 1))
 *               = (1/2)(1 + (2 Phi - 1)(1 - 4(Phi - Phi^2))^{-1/2}).
 *   "The right-hand side involves a Taylor expansion in 4(Phi - Phi^2), which
 *    converges if eta < 1/4."  Properties (.tex:2177-2182):
 *     Phi_tilde^2 = Phi_tilde,  ||Phi_tilde - Phi||_cb <= O(eta),
 *     Phi_tilde(1) = 1,  Phi_tilde(X^dag) = Phi_tilde(X)^dag.
 *
 * This is Proposition prop_P (.tex:524-533) applied to P = Phi at the
 * superoperator level: build S_Phi (aic_assoc_superop_from_ucp), then
 * S_tilde = aic_prop_P(S_Phi) = theta(2 S_Phi - 1). aic_prop_P CERTIFIES the
 * prop_P hypothesis as the SPECTRAL condition rho(S_Phi^2 - S_Phi) < 1/4 (.tex:525)
 * via the eig-free Gelfand bound (bead aic-8hz, LANDED) and aborts loudly at the
 * true rho >= 1/4 boundary (Rule 4) — that IS the fail-loud basin guard. The sign
 * inside theta auto-dispatches to the globally-convergent Newton out of the
 * operator-norm basin (aic-8hz), so the full NON-NORMAL eta < 1/4 regime is now
 * reachable (the earlier op-norm guard ||S_Phi^2 - S_Phi||_op < 1/4 over-rejected
 * oblique channels; test_assoc2 U6 is the payoff witness).
 *
 * ROUTE A1 (default, here): prop_P -> theta -> aic_sgn (Newton-Schulz, eig-free).
 * ROUTE A2 (the independent non-normal cross-check, lives in the test, test T6):
 *   the explicit binomial-series form of the third equality above:
 *     D = 4(Phi - Phi^2) = 4(S - S^2);  M = I - D = 1 - 4(Phi - Phi^2) = (2S-I)^2;
 *     Minvsqrt = aic_funcalc_xpow(M, alpha=-1/2, x0=1) = (X^2)^{-1/2};
 *     S_tilde = (1/2)(I + (2 S - I) Minvsqrt).
 *   NOTE the sign: M = 1 - 4(Phi - Phi^2), so D = 4(S - S^2) (NOT 4(S^2 - S)).
 *   The eta=0 case (D=0) is blind to this sign; the eta>0 non-normal T6 channel
 *   is what pins it (a flipped sign gives sgn^2 != I -> A1/A2 disagree by O(0.1)).
 *   A1 and A2 are algorithmically independent (Newton-Schulz vs binomial Taylor)
 *   and must agree to ~1e-25 at prec=256 on a NON-NORMAL S_Phi — both an
 *   independent check of the regularization AND coverage for funcalc's
 *   never-before-tested non-normal sgn path.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_funcalc.h"
#include "aic/aic_ucp.h"
#include "aic/aic_assoc.h"
#include "aic_funcalc_internal.h"

/* approximate_algebras.tex:2168 — the almost-idempotent hypothesis:
 *     ||Phi^2 - Phi||_cb <= eta,
 * and .tex:2176 (eq. tilde_Phi): "The right-hand side involves a Taylor expansion
 * in 4(Phi - Phi^2), which converges if eta < 1/4."  So Phi_tilde = theta(2 Phi-1)
 * is well-defined ONLY when the channel is inside the eta < 1/4 basin. For the
 * finite-dim SUPEROPERATOR S = S_Phi this is the SPECTRAL condition (NOT the cb- or
 * op-norm one): rho(S^2 - S) < 1/4, equivalently rho(I - (2S-I)^2) = 4 rho(S^2-S)
 * < 1 (.tex:516,520-521; the prop_P basin .tex:524-525). The op-norm guard
 * ||S^2-S||_op < 1/4 is STRICTER and over-rejects valid NON-NORMAL channels (bead
 * aic-8hz; test_assoc2 U6: ||S^2-S||_op = 0.420 >= 1/4 but rho = 0.210 < 1/4), so
 * this gate uses the eig-free certified Gelfand spectral radius, exactly as
 * aic_prop_P does — NOT a cheap op-norm gate (which would FALSE-reject U6).
 *
 * WHY AN ENTRY GUARD when aic_prop_P already certifies the same basin deeper down
 * (and aborts fast — measured ~0.2 s, NOT a hang): the deep aic_prop_P message is
 * GENERIC ("P^2-P", prop_P) with no channel-level attribution. An out-of-regime
 * CHANNEL fed to the pipeline entry should fail loud AT THE API BOUNDARY, naming the
 * map aic_assoc_regularize, the eta < 1/4 almost-idempotent hypothesis, and the
 * measured spectral radius — so the failure is attributable to the input channel,
 * not to an internal functional-calculus call (CLAUDE.md Rule 4; bead aic-xo0;
 * FINDINGS §C15, which also records WHY make_mixconj cannot be the out-of-regime
 * fixture — it is spectrally in-basin at every t). The gate runs on the SAME S the
 * regularization then consumes; the small redundancy
 * with prop_P's own guard is intentional (prop_P's guard stays intact for its other
 * callers). The superop build is O(n^4) but cheap (~20 ms at n=10); the certified
 * Gelfand scan is k_max-bounded (no hang). */
static void aic_assoc_int_assert_eta_basin(const acb_mat_t S, slong n, slong prec)
{
    slong nn = n * n;
    acb_mat_t X, M, I;
    acb_mat_init(X, nn, nn);
    acb_mat_init(M, nn, nn);
    acb_mat_init(I, nn, nn);
    acb_mat_one(I);
    acb_mat_scalar_mul_2exp_si(X, S, 1);     /* 2S */
    acb_mat_sub(X, X, I, prec);              /* X = 2S - I */
    acb_mat_sqr(M, X, prec);                 /* X^2 */
    acb_mat_sub(M, I, M, prec);              /* M = I - X^2 = 4(S - S^2) */

    arb_t rho;
    arb_init(rho);
    slong kused = 0;
    int cert = aic_funcalc_int_gelfand_rho(rho, &kused, M, 32, prec);
    if (!cert) {
        /* rho_ub holds the tightest (largest-k) bound on 4 rho(S^2-S). */
        fprintf(stderr,
                "aic_assoc_regularize: channel is OUT OF the eta < 1/4 "
                "almost-idempotent regime (.tex:2168/2176). The regularization "
                "Phi_tilde = theta(2 Phi - 1) needs rho(Phi^2 - Phi) < 1/4 "
                "(equivalently rho(I-(2 S_Phi-I)^2) = 4 rho(S_Phi^2-S_Phi) < 1; the "
                "prop_P basin .tex:524-525). The certified Gelfand bound on "
                "4 rho(S_Phi^2-S_Phi) did NOT fall < 1 by k=32 (ub ~ %.4g) — the "
                "Taylor series in 4(Phi-Phi^2) does not converge, so Phi is at/over "
                "the eta=1/4 boundary (raise prec only if borderline). Failing loud "
                "at the pipeline entry (CLAUDE.md Rule 4; bead aic-xo0).\n",
                arf_get_d(arb_midref(rho), ARF_RND_NEAR));
        abort();
    }
    arb_clear(rho);
    acb_mat_clear(I);
    acb_mat_clear(M);
    acb_mat_clear(X);
}

void aic_assoc_regularize(acb_mat_t S_tilde, const aic_ucp_kraus *phi,
                          slong prec)
{
    slong n = phi->dim_H;
    assert(phi->dim_K == n &&
           "aic_assoc_regularize: Phi must be a self-map (dim_K==dim_H)");
    assert(acb_mat_nrows(S_tilde) == n * n && acb_mat_ncols(S_tilde) == n * n &&
           "aic_assoc_regularize: S_tilde must be n^2 x n^2");

    /* Build S_Phi, then certify the eta < 1/4 basin AT THIS ENTRY (channel-level
     * fail-loud message, bead aic-xo0), then Phi_tilde = prop_P(S_Phi) =
     * theta(2 S_Phi - 1). prop_P re-certifies the same basin internally (its guard
     * stays intact for its other callers); this entry guard is what attributes an
     * out-of-regime failure to the input channel at the API boundary. */
    acb_mat_t S;
    acb_mat_init(S, n * n, n * n);
    aic_assoc_superop_from_ucp(S, phi, prec);
    aic_assoc_int_assert_eta_basin(S, n, prec);
    aic_prop_P(S_tilde, S, prec);
    acb_mat_clear(S);
}
