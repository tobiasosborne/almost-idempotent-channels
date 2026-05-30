/* aic_projection.h — a constructive nontrivial O(eps)-projection finder for an
 * eps-C* algebra A (bead aic-mqf, module "projection"). Realizes
 * lem_nontriv_projection (approximate_algebras.tex:931-932):
 *
 *   "Any eps-C* algebra A such that 1 < dim A < infinity has a nontrivial
 *    O(eps)-projection."
 *
 * A delta-projection (.tex:917-919) is a Hermitian P with ||P^2 - P|| <= delta,
 * where the product P^2 = P * P is the Choi-Effros STAR (.tex:341, X*Y=Phi(XY)),
 * NOT the plain matrix product. It is NONVANISHING (.tex:926-928) if
 * |||P|| - 1| <= O(delta+eps), and NONTRIVIAL (.tex:929) if both P and I-P are
 * nonvanishing. So the output P is bounded away from both 0 and I by Omega(1).
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The paper proves
 * EXISTENCE non-constructively: it reduces to a fixed point of the inversion map
 * sigma(U)=U^dag on the approximate unitary group U (.tex:934-943), then invokes
 * the Lefschetz-Hopf fixed-point theorem + an H-space cohomology argument
 * (.tex:944-969, prop_H-group) to force a nontrivial fixed point. NONE of that is
 * an algorithm. In finite dimensions we build a projection directly (Route
 * A-AMBIENT, the route reconciled across three research legs, bead aic-mqf):
 *
 *   1. FIND A NON-SCALAR HERMITIAN H in A. Each H_k = (1/2)(B_k + B_k^dag) is
 *      Hermitian and in A (A is dag-closed). Pick the H_k with the largest
 *      INTERIOR SPECTRAL GAP (tie-break by spread), NOT the largest spread:
 *      argmax(spread) != argmax(gap), and the a-posteriori star defect is
 *      delta=O(eps/g), so the largest-gap H gives the tightest certified defect
 *      (and avoids the spurious-abort failure mode where a high-spread H_k's
 *      largest gap fell below the floor). All-near-scalar (spread <= 1e-6) =>
 *      FAIL LOUD (a stop condition: dim A > 1 but no non-scalar Hermitian elt).
 *   2. FIND THE SPECTRAL GAP + THRESHOLD. H is Hermitian => its spectrum is STABLE
 *      (the .tex:540-544 fragility is a NON-NORMAL phenomenon), so the double-path
 *      eigenVALUES (aic_latd_eig_hermitian / LAPACK zheev, degenerate-OK) are
 *      sound for choosing a threshold. Take the largest INTERIOR gap [lam_m,
 *      lam_{m+1}] (1 <= m < n, both sides non-empty => the split is guaranteed
 *      nontrivial REGARDLESS of gap size), t = (lam_m + lam_{m+1})/2,
 *      g = lam_{m+1} - lam_m. FAIL LOUD only if NO positive interior gap exists
 *      (the spectrum is degenerate, all eigenvalues coincide to round-off): the
 *      aic-3qv stop condition. The gap SIZE is NOT gated here; it enters the
 *      a-posteriori star defect delta=O(eps/g), certified downstream in Step 5
 *      and guarded by the T3 universality canary.
 *   3. AMBIENT SPECTRAL PROJECTOR via the EIG-FREE sgn. s = 1/max(t-lam_min,
 *      lam_max-t); Y = s*(H - t*I) is n x n Hermitian with eigenvalues in [-1,1]
 *      and |eigval| >= s*g/2 off 0, so the funcalc basin ||Y^2-I|| < 1 holds
 *      (asserted). X = sgn(Y) is +/-1 on the two clusters (AMBIENT, in M_n where
 *      the C* structure is EXACT — NOT inside A, which rem_X2 (.tex:628-630)
 *      forbids); P_amb = (1/2)(I + X) is an EXACT ambient idempotent.
 *   4. PROJECT INTO A via the ALGEBRA'S OWN idempotent Phi_tilde. A = Img Phi_tilde
 *      is an OBLIQUE image (Phi_tilde is HP but not HS-self-adjoint), so the
 *      structure-respecting projection is Phi_tilde itself, NOT the Frobenius-
 *      orthogonal projector: P = Phi_tilde(P_amb) = P_amb * I (the star with the
 *      unit, .tex:2211; Phi_tilde(X)=X*I). P is in A and Hermitian (Phi_tilde is
 *      HP, P_amb Hermitian). At eta=0, Phi_tilde=id so P=P_amb; at eta>0, the
 *      Frobenius route would leave an O(1) star defect (measured ~0.5, constant
 *      in eta) whereas Phi_tilde(P_amb) gives the O(eta) defect the lemma promises
 *      (measured C=delta/eta ~0.04, flat across eta and dim) — bead aic-mqf.
 *   5. CERTIFY (arb). delta = ||P*P - P||_op (the STAR, via aic_ecstar_star +
 *      the certified aic_corner_gamma_opnorm_ub upper bound, which dodges the
 *      aic_mat_opnorm Gram false-fail on near-zero defects, bead aic-qgs).
 *      Nontriviality: ||P||_op, ||I-P||_op both bounded away from 0 and I.
 *      Membership: aic_ecstar_proj_residual(P) = O(eta). FAIL LOUD on any failure.
 *
 * WHY rem_X2-SAFE. The ONLY functional calculus (aic_sgn) is applied to the
 * AMBIENT Hermitian matrix Y in M_n, a genuine C* algebra where sgn is exact.
 * We never apply sgn/theta INSIDE A (the eps-C* algebra), which Remark rem_X2
 * (.tex:628-630) shows is ill-posed (prop_P does not generalize to the eps-C*
 * setting). The O(eta) gap from P_amb to A is absorbed by Pi_A and CERTIFIED by
 * the star-defect ball, not assumed away.
 *
 * AUDITION (Law 4). Route A-AMBIENT is the PRIMARY candidate. Route B (a
 * sigma-fixed-point search on the approximate unitary group, depends on the
 * unitary-group module aic-q2x) and Route C (minimize ||P*P-P|| over Hermitian
 * P bounded from +/-I, needs a warm start) are Pareto alternatives, deferred
 * until a benchmark motivates them.
 *
 * NUMBER PATH. The gap-finding eigenVALUES are the fast double path (LAPACK,
 * stable for Hermitian); the spread certification (aic_mat_herm_max_eig), the
 * ambient sgn, and ALL the output defects are the certified arb path. Certified
 * eigenvalue ENCLOSURES for a tighter gap floor defer to bead aic-w4o.1 (the
 * gap-finding uses double-path zheev now, the DEFECT is arb-certified).
 */
#ifndef AIC_PROJECTION_H
#define AIC_PROJECTION_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_ecstar.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The chosen Hermitian element H in A and its spectral split, exposed so tests
 * can inspect WHICH H/gap/threshold the finder selected (and assert the gap is
 * the one a hand analysis predicts). All fields are OUTPUT. */
typedef struct {
    slong k_chosen;   /* basis index k whose H_k = (B_k+B_k^dag)/2 was selected   */
    double lam_min;   /* smallest eigenvalue of the chosen H                      */
    double lam_max;   /* largest eigenvalue                                       */
    double t;         /* threshold = (lam_m + lam_{m+1})/2 at the largest gap     */
    double gap;       /* g = lam_{m+1} - lam_m (the largest interior gap)         */
    slong m;          /* gap index: lam_0..lam_{m-1} below t, lam_m..lam_{n-1} above
                       * (so m eigenvalues are >= t; the +1-side cluster size)    */
} aic_projection_witness;

/* Find a nontrivial O(eps)-projection P in A (lem_nontriv_projection). Fills:
 *   P        : n x n, Hermitian, in A, with P*P ~ P (the STAR). Caller-init'd
 *              n x n where n = A->n.
 *   delta    : certified arb ball >= ||P*P - P||_op (the star defect).
 *   Pnorm    : certified arb ball for ||P||_op (a nonvanishing witness, ~1).
 *   ImPnorm  : certified arb ball for ||I - P||_op (~1).
 *   wit      : OPTIONAL (NULL to skip) the chosen H / gap / threshold.
 * ASSERTS (fail loud, Rule 4):
 *   - 1 < dim A (the lemma's hypothesis; dim A == 1 has no nontrivial projection),
 *   - some H_k is non-scalar (spread > a floor); else escalate (stop condition),
 *   - a positive interior spectral gap exists (else the aic-3qv escalation:
 *     a degenerate spectrum with all eigenvalues coincident),
 *   - the ambient sgn basin ||Y^2 - I|| < 1 holds,
 *   - the output is a genuine nontrivial delta-projection in A: delta small,
 *     ||P|| and ||I-P|| both bounded away from 0 (>= 0.3) AND from each other's
 *     trivial value, and proj_residual(P) small.
 * The caller owns P/delta/Pnorm/ImPnorm (init before, clear after). */
void aic_projection_nontrivial(acb_mat_t P, arb_t delta, arb_t Pnorm,
                               arb_t ImPnorm, const aic_ecstar *A,
                               aic_projection_witness *wit, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_PROJECTION_H */
