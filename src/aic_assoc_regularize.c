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

#include <flint/acb_mat.h>

#include "aic_funcalc.h"
#include "aic_ucp.h"
#include "aic_assoc.h"

void aic_assoc_regularize(acb_mat_t S_tilde, const aic_ucp_kraus *phi,
                          slong prec)
{
    slong n = phi->dim_H;
    assert(phi->dim_K == n &&
           "aic_assoc_regularize: Phi must be a self-map (dim_K==dim_H)");
    assert(acb_mat_nrows(S_tilde) == n * n && acb_mat_ncols(S_tilde) == n * n &&
           "aic_assoc_regularize: S_tilde must be n^2 x n^2");

    /* Build S_Phi, then Phi_tilde = prop_P(S_Phi) = theta(2 S_Phi - 1). The
     * prop_P precondition ||S_Phi^2 - S_Phi||_op < 1/4 (.tex:525, the eta<1/4
     * basin) is asserted inside aic_prop_P (fail loud). */
    acb_mat_t S;
    acb_mat_init(S, n * n, n * n);
    aic_assoc_superop_from_ucp(S, phi, prec);
    aic_prop_P(S_tilde, S, prec);
    acb_mat_clear(S);
}
