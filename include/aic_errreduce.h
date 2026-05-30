/* aic_errreduce.h — cor_improvement (ERROR REDUCTION), bead aic-t81, module
 * "errreduce". Realizes approximate_algebras.tex:1317-1319:
 *
 *   Corollary (Error reduction, cor_improvement). There exist positive constants
 *   eps_max, delta_max, c_0 such that for all eps < eps_max, if a finite-dim C*
 *   algebra B is delta_max-included into an eps-C* algebra A, there is also a
 *   c_0*eps-inclusion. If the original inclusion is bijective, so is the new one.
 *
 * WHAT THIS MODULE IS. THIN orchestration over the dhom module (§8). The heavy
 * lifting — lem_approx (.tex:1224, the multiplicativity-defect Newton iteration
 * driven by the explicit generalized-Pauli diagonal of the GENUINE C* algebra B)
 * and prop_delta_hominc (.tex:1194, the norm bounds) — lives in src/aic_dhom_*.c.
 * errreduce wires them into the cor_improvement statement: it takes a
 * delta-INCLUSION v (a delta-homomorphism that is ALSO bounded below,
 * (1-delta)||X|| <= ||v(X)||, .tex:451-453), produces vtilde = lem_approx(v) (an
 * O(eps)-HOMOMORPHISM, .tex:1225), and CERTIFIES that vtilde is an O(eps)-
 * INCLUSION via prop_delta_hominc applied to vtilde (whose own multiplicativity
 * defect is now O(eps)).
 *
 * WHY cor_improvement IS THE th_main ENGINE (.tex:490). The master-loop proof of
 * th_main resets the inclusion error at EVERY merging/extension step: a freshly
 * extended v is only a delta-inclusion for delta > delta_0; cor_improvement pulls
 * it back to a fixed delta_0 = c_0*eps that does NOT depend on delta or on dim A.
 *
 * THE UNIVERSALITY CANARY (.tex:484, the headline test). The implicit constant in
 * th_main's O(eps) is DIMENSION-INDEPENDENT (.tex:461). The paper warns (.tex:484)
 * that the naive Haar/cohomology diagonal has error proportional to n = dim A; our
 * route avoids this because the diagonal used inside lem_approx is the EXACT
 * generalized-Pauli diagonal of the GENUINE C* algebra B (not of the eps-
 * associative A). So the MEASURED c_0 must NOT grow with dim A / dim B. test_
 * errreduce T4 is the canary: it sweeps B over increasing dim and asserts the
 * measured c_0 ratio across the sweep is bounded by a constant. A c_0 proportional
 * to dim is the .tex:484 failure mode and a stop condition (escalate aic-1bc).
 *
 * THE c_0 RETURNED IS MEASURED, NOT ANALYTIC. The paper's c_0 is "explicit" in
 * principle but unstated numerically (paper/FINDINGS.md D2; the analytic extraction
 * is bead aic-1bc). This module returns the MEASURED c_0 = max(inclusion-defects of
 * vtilde) / eps for the given instance; that is the empirical witness the canary
 * tracks against the universal-constant claim.
 *
 * THE STAR IS LOAD-BEARING (inherited via dhom, CLAUDE.md / FINDINGS C2). Every
 * A-product inside lem_approx/G_v is the Choi-Effros star X*Y = Phi_tilde(XY)
 * (aic_ecstar_star), NEVER plain matmul. The eta=0 identity-channel oracle has
 * star == plain and is BLIND to a star bug; only the eta>0 fixtures (T3/T4)
 * exercise the star non-trivially.
 *
 * FAIL-LOUD (Rule 4). aic_errreduce ASSERTS (i) the input is a delta-inclusion —
 * its COORDINATE-MATRIX sigma_min (aic_dhom_v_sigma_min, the TRUE Frobenius
 * unit-ball inclusion infimum) is >= 0.5, i.e. v is bounded below / injective and
 * does NOT collapse any direction (a delta-inclusion with delta <= delta_max has
 * a >= 1-delta_max). This is the F1 fix: the old basis-sweep min_i ||vE[i]|| was
 * BLIND to a v bounded below on each basis element that collapses a combination
 * (paper/FINDINGS.md). And (ii) the output is certifiably an O(eps)-inclusion
 * (every inclusion-defect of vtilde — including lower_gap = max(0, 1-sigma_min) —
 * is <= c0_cert * max(eps, tol)). It aborts at the call site otherwise; corrupted
 * output is worse than a crash.
 */
#ifndef AIC_ERRREDUCE_H
#define AIC_ERRREDUCE_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_dhom.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The four INCLUSION-DEFECTS of vtilde (cor_improvement certificate), as certified
 * arb balls. vtilde is a c_0*eps-INCLUSION iff all four are <= c_0*eps:
 *   mult_def : basis-sweep multiplicativity defect max_ij ||G_vtilde(E_i,E_j)||_op
 *              (prop_delta_hominc / lem_approx: O(eps) after the Newton iteration).
 *   norm_excess : max(0, ||vtilde|| - 1)        (prop_delta_hominc (i): <= O(eps)).
 *   lower_gap   : max(0, 1 - a),  a = sigma_min of v's coordinate matrix (the TRUE
 *                 Frobenius unit-ball inclusion infimum, F1 fix — NOT the basis
 *                 sweep; prop_delta_hominc (ii): a >= 1 - O(eps), lower_gap <= O(eps)).
 *   unit_def    : ||vtilde(I_B) - I_A||_op       (prop_delta_hominc (iii): <= O(eps)).
 * Each output arb_t is caller-initialised; pass NULL to skip any one. */
typedef struct {
    arb_t mult_def;
    arb_t norm_excess;
    arb_t lower_gap;
    arb_t unit_def;
} aic_errreduce_defects;

void aic_errreduce_defects_init(aic_errreduce_defects *d);
void aic_errreduce_defects_clear(aic_errreduce_defects *d);

/* The largest of the four inclusion-defects (the certificate magnitude:
 * vtilde is a (max/eps)*eps-inclusion). `out` caller-initialised. */
void aic_errreduce_defects_max(arb_t out, const aic_errreduce_defects *d);

/* cor_improvement (.tex:1317): upgrade a delta-inclusion v_in to a c_0*eps-
 * inclusion vtilde_out.
 *
 *   vtilde_out : OUTPUT delta-hom over the SAME B and A as v_in (a fresh copy; its
 *                vE[] are allocated here — init it with aic_dhom_v_init(.,B,A)
 *                first, or use aic_errreduce_copy_v below). On return its basis
 *                holds vtilde = lem_approx(v_in).
 *   c0_out     : OUTPUT measured c_0 = max(inclusion-defects)/eps (caller-init'd
 *                arb_t). If eps is below `tol_abs` (the eta=0 case), the defects
 *                must be ~machine-zero and c0_out is set to 0 (no divide-by-zero).
 *   defects    : OUTPUT the four certified inclusion-defects (caller-init'd via
 *                aic_errreduce_defects_init; NULL to skip).
 *   v_in       : the input delta-inclusion (BORROWED, unchanged).
 *   eps        : the eps of A (e.g. aic_ecstar_defect_assoc). The O(eps) floor for
 *                the lem_approx termination and the c_0 denominator.
 *   tol_abs    : absolute machine floor (e.g. 1e-13 at prec ~80); added to eps for
 *                the Newton stopping threshold and used as the eps~0 cutoff.
 *   max_steps  : lem_approx Newton cap (fail-loud if not reached, .tex:1312).
 *
 * ASSERTS (fail loud): v_in is a delta-inclusion (coordinate-matrix sigma_min
 * a >= 0.5, i.e. bounded below / injective / no collapsed direction, F1 fix), and
 * vtilde is certifiably an O(eps)-inclusion (max inclusion-defect <=
 * AIC_ERRREDUCE_C0_CERT * max(eps, tol_abs)). */
void aic_errreduce(aic_dhom_v *vtilde_out, arb_t c0_out,
                   aic_errreduce_defects *defects, const aic_dhom_v *v_in,
                   double eps, double tol_abs, slong max_steps, slong prec);

/* The certification ceiling: vtilde's max inclusion-defect must be
 * <= AIC_ERRREDUCE_C0_CERT * max(eps, tol_abs), else aic_errreduce aborts. This is
 * a GENEROUS fail-loud guard (a real O(eps)-inclusion has max-defect a few * eps),
 * NOT the analytic c_0 (aic-1bc); it catches a non-converged or non-inclusion
 * output, not the precise constant. Tightened 50 -> 10 (F3): the worst measured
 * max-defect over the test corpus is ~2.71*eps (T4(A) M_2), so 10 keeps a ~3.7x
 * margin while no longer admitting a 16x-off non-inclusion; the machine-floor
 * cases (T1/T2, eps_eff = tol) clear it with >=30x margin (max-defect/tol ~0.33). */
#define AIC_ERRREDUCE_C0_CERT 10.0

/* The lem_approx termination floor MULTIPLE. lem_approx reduces the multiplicativity
 * defect to O(eps) — concretely to ~c_0*eps with c_0 a few (FINDINGS D2; the
 * achievable floor is a SMALL MULTIPLE of eps, NOT eps itself, because the Newton
 * step cannot push the defect below the algebra's intrinsic O(eps) associativity
 * obstruction). So aic_errreduce stops the iteration at AIC_ERRREDUCE_EPS_FLOOR*eps
 * (+ tol_abs); driving toward bare eps would stall and bounce (Newton can't beat
 * the O(eps) floor — caught fail-loud by aic_dhom_approx's contraction guard). The
 * c_0 DENOMINATOR stays the bare eps, so the reported c_0 is the true defect/eps. */
#define AIC_ERRREDUCE_EPS_FLOOR 4.0

/* Deep-copy v_in's stored basis into v_out (which must be aic_dhom_v_init'd over
 * the SAME B and A). Convenience for setting up vtilde_out. */
void aic_errreduce_copy_v(aic_dhom_v *v_out, const aic_dhom_v *v_in);

/* Bijectivity (cor_improvement last sentence, .tex:1318): a delta-ISOMORPHISM
 * (bijective delta-inclusion) upgrades to a c_0*eps-ISOMORPHISM. In finite dim,
 * vtilde is bijective iff it is injective (coordinate-matrix sigma_min > 0, the
 * SOUND collapse detector — F1 fix; a kernel direction is now caught) AND
 * dim_B == dim_A. Returns nonzero iff both hold; writes the sigma_min lower bound
 * a to `a_out` (caller-init'd arb_t, NULL to skip) so the caller can confirm a > 0.
 * Does NOT abort (a diagnostic, not a guard). */
int aic_errreduce_is_bijective(arb_t a_out, const aic_dhom_v *vtilde,
                               slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_ERRREDUCE_H */
