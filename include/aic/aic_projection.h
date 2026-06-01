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
 *   1+2. ENUMERATE (H_k, gap) CANDIDATES + AUDITION UNIT-AWARE (bead aic-66n,
 *      src/aic_projection_audit.c). Each H_k = (1/2)(B_k + B_k^dag) is Hermitian
 *      and in A (A is dag-closed); its double-path zheev spectrum is STABLE (the
 *      .tex:540-544 fragility is a NON-NORMAL phenomenon). For EVERY non-scalar
 *      H_k, enumerate ALL interior gaps [lam_i, lam_{i+1}] (t = (lam_i+lam_{i+1})/2,
 *      g = lam_{i+1}-lam_i, both clusters non-empty), collected into one candidate
 *      list sorted by g DESCENDING. The OLD route pre-picked ONE H_k by the largest
 *      AMBIENT gap and built P from it — but on an S_P wrapper (unit Ptilde_m =
 *      Co_P(P), .tex:1082; NOT 1_n) the largest gap is the support-vs-complement
 *      gap, so P collapsed to the wrapper UNIT (FINDINGS §C11/§C16). Instead we
 *      AUDITION candidates largest-gap-first (Steps 3-4 build P), ACCEPTING the
 *      first that is NONTRIVIAL VS THE ALGEBRA UNIT U_A = Phi_tilde(1_n): ||P|| > c
 *      AND ||U_A - P|| > c (c = 0.3, .tex:929; FINDINGS §C7). Auditioning across
 *      ALL H_k (not one pre-chosen H) is ROBUST: if the largest-ambient-gap H_k has
 *      its in-support eigenvalues clustered while a different H_k has a usable gap,
 *      we still find it. FAIL LOUD only if NO non-scalar H_k has a positive interior
 *      gap, or NO candidate yields a nontrivial-vs-unit split: the aic-3qv stop
 *      condition. The gap SIZE is NOT gated; it enters the a-posteriori star defect
 *      delta=O(eps/g), certified in Step 5 and guarded by the T3 universality canary.
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
 *      Nontriviality (UNIT-AWARE): ||P||_op, ||U_A - P||_op both bounded away from
 *      0, where U_A = Phi_tilde(1_n) is the ALGEBRA unit (= Ptilde_m for an S_P
 *      wrapper, ~1_n at top level; the old ambient ||1_n - P|| gate was vacuous on
 *      a wrapper, FINDINGS §C11/§C16). FAIL LOUD on any failure.
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

#include "aic/aic_ecstar.h"

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
 *   ImPnorm  : certified arb ball for ||U_A - P||_op (~1; U_A the ALGEBRA unit).
 *   wit      : OPTIONAL (NULL to skip) the chosen H / gap / threshold.
 * ASSERTS (fail loud, Rule 4):
 *   - 1 < dim A (the lemma's hypothesis; dim A == 1 has no nontrivial projection),
 *   - some non-scalar H_k has a positive interior spectral gap; else escalate
 *     (the aic-3qv stop condition: a degenerate algebra),
 *   - the ambient sgn basin ||Y^2 - I|| < 1 holds for the chosen candidate,
 *   - SOME (H_k, gap) candidate yields a nontrivial-vs-unit split; else aic-3qv,
 *   - the output is a genuine nontrivial delta-projection in A: delta small,
 *     ||P|| and ||U_A - P|| both bounded away from 0 (>= 0.3), where U_A =
 *     Phi_tilde(1_n) is the ALGEBRA unit (FINDINGS §C16; NOT the ambient 1_n).
 * The caller owns P/delta/Pnorm/ImPnorm (init before, clear after). */
void aic_projection_nontrivial(acb_mat_t P, arb_t delta, arb_t Pnorm,
                               arb_t ImPnorm, const aic_ecstar *A,
                               aic_projection_witness *wit, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_PROJECTION_H */
