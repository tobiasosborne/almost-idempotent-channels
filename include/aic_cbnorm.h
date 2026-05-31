/* aic_cbnorm.h — eig-free certified two-sided bracket on the eta-idempotence
 * defect (bead aic-m24, increment 1). Realizes the cb-norm / diamond-norm
 * quantity of approximate_algebras.tex:347-354:
 *
 *   .tex:347-352 (verbatim, the cb-norm):
 *     ||Lambda||_cb = sup_n sup_{X!=0} ||(1_{M_n} (x) Lambda)(X)|| / ||X||,
 *     "efficiently computable using semidefinite programming"; the dual norm
 *     for maps of functionals is the diamond norm ||.||_diamond.
 *   .tex:354 (verbatim): "A UCP map Phi: B(H)->B(H) is called eta-idempotent if
 *     ||Phi^2-Phi||_cb <= eta."
 *
 * For a Hermiticity-preserving map (Lambda = Phi^2 - Phi is such), cb-norm =
 * diamond-norm and the ampliation truncation N = dim(input) = n is rigorous
 * (Watrous 2012); see ALGORITHM.md cb-norm section (bead aic-d24) for the
 * (2/n) normalization of the Watrous SDP relative to the Convention-A Choi.
 *
 * THIS MODULE: the cheapest, ALWAYS-VALID rung of the certified-bound ladder.
 * It uses ONLY the certified Frobenius norm of J = Choi(Lambda) — no SDP solver,
 * no eigendecomposition (so it dodges the degenerate-eig wall, bead aic-w4o.1) —
 * and it NEVER aborts. The bracket [||J||_F/n, 2*||J||_F] has ratio 2n (loose by
 * design); the tight ball is a later increment. This is the fallback the future
 * tight certifier dispatches to near eta=0. See aic_cbnorm.c for the derivation.
 */
#ifndef AIC_CBNORM_H
#define AIC_CBNORM_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Eig-free certified two-sided bracket on eta = ||Phi^2-Phi||_cb, computed from
 * J = Choi(Lambda) (the Convention-A Choi of Lambda = Phi^2 - Phi, n^2 x n^2).
 * On return the arb balls satisfy, rigorously,
 *   arb_lower(lo) <= eta <= arb_upper(hi),  with lo ~ ||J||_F/n, hi ~ 2*||J||_F.
 * Never aborts; eta=0 (J=0) gives [0,0]. nrows(J) == ncols(J) == n*n is asserted
 * (fail loud on a shape mismatch — Rule 4). `lo`, `hi` are caller-initialised
 * arb_t. */
void aic_cbnorm_eigfree_ball_choi(arb_t lo, arb_t hi, const acb_mat_t J,
                                  slong n, slong prec);

/* Same bracket, starting from a UCP self-map Phi (dim_K == dim_H == n,
 * asserted): builds Phi^2 via aic_ucp_compose and J = Choi(Phi^2) - Choi(Phi)
 * via aic_ucp_choi_diff, then calls aic_cbnorm_eigfree_ball_choi. `lo`, `hi` are
 * caller-initialised arb_t. */
void aic_cbnorm_eigfree_ball(arb_t lo, arb_t hi, const aic_ucp_kraus *phi,
                             slong prec);

/* TIGHT certified two-sided bracket on eta = ||Phi^2-Phi||_cb (bead aic-m24,
 * increment 3b). Inputs: J = Choi(Phi^2 - Phi) (n^2 x n^2, Convention A) and TWO
 * SDP feasible points produced by the Julia+MOSEK solves of step 3a, each
 * n^2 x n^2:
 *   - MAX-primal (X,P,Q): max Re tr(J^dag X) s.t. [[P,X],[X^dag,Q]]>=0,
 *     P+Q=I_{n^2}; gives the rigorous LOWER bound eta >= (2/n) Re tr(J^dag X).
 *   - MIN-dual (Y0,Y1): min (1/2)(||Tr_2 Y0||_inf + ||Tr_2 Y1||_inf) s.t.
 *     [[Y0,-J],[-J^dag,Y1]]>=0, Y0,Y1>=0; gives the rigorous UPPER bound
 *     eta <= (1/2)(||Tr_2 Y0||_inf + ||Tr_2 Y1||_inf), Tr_2 = partial_trace_right.
 * Both feasible points are RESTORED to exact feasibility in arb (convex-
 * combination toward the Slater center for the primal; eigenvalue shift for the
 * dual) so the returned ball is RIGOROUS: arb_lower(lo) <= eta <= arb_upper(hi),
 * and TIGHT (gap ~ MOSEK duality gap + arb radius, much tighter than the
 * eig-free 2n-wide bracket). Dispatches to aic_cbnorm_eigfree_ball_choi in the
 * eta=0 regime (||J||_F < ~1e-9*n) or when the MOSEK points are too far off to
 * restore cleanly; fails loud (Rule 4) if even the eig-free fallback straddles.
 * All ONE-SIDED herm_max_eig(-M) PSD tests (no full eig; dodges aic-w4o.1).
 * See docs/cbnorm_tight_certifier.md. `lo`, `hi` are caller-initialised arb_t. */
void aic_cbnorm_certify(arb_t lo, arb_t hi, const acb_mat_t J,
                        const acb_mat_t X, const acb_mat_t P, const acb_mat_t Q,
                        const acb_mat_t Y0, const acb_mat_t Y1,
                        slong n, slong prec);

/* RECTANGULAR certified UPPER bound on a diamond norm, for the opspace O2 cb-norm
 * payoff (bead aic-pjr; src/aic_cbnorm_certify_rect.c). Given the GOLDEN-RULE
 * Convention-A Choi J of a map f: M_{d_maj} -> M_{d_min}
 *   J[s*d_min+i, t*d_min+j] = f(E_st)[i,j]   (INPUT s,t MAJOR; OUTPUT i,j MINOR),
 * each of J, Y0, Y1 (d_maj*d_min) square, and a Watrous MIN-dual feasible point
 * (Y0,Y1), this restores (Y0,Y1) to exact feasibility by a shift and reads off a
 * RIGOROUS upper bound on ||f||_⋄ (= the SDP optval, design §0.5: normalization
 * FACTOR 1, NO 2/n). The dual traces the MINOR/OUTPUT factor (tr_sys = 2, size
 * d_min): Tr_min = aic_mat_partial_trace_right(., d_maj, d_min); the shift term is
 * eps*d_min; hi = (1/2)(lambda_max(Tr_min Y0)+lambda_max(Tr_min Y1)) + eps*d_min.
 * NO eta=0 short-circuit (the trivial value here is 1, not 0 — the complete-
 * isometry oracle goes through the full restoration). On return
 * eta = ||f||_⋄ <= arb_upper(hi), rigorously. The diamond-of-adjoint duality
 *   ||v||_cb = ||v*||_⋄,  ||v^{-1}||_cb = ||(v^{-1})*||_⋄  (Watrous 2009)
 * turns this into the certified cb-norm UPPER bound of th_main_ext (.tex:1538-1540)
 * when J = J(v*) (d_maj=N, d_min=n_B) resp. J = J((v^{-1})*) (d_maj=n_B, d_min=N).
 * All ONE-SIDED herm_max_eig(-M) PSD tests (degeneracy-robust, no full eig).
 * `hi` is a caller-initialised arb_t. ASSERTS the shapes (Rule 4). */
void aic_cbnorm_certify_rect_upper(arb_t hi, const acb_mat_t J,
                                   const acb_mat_t Y0, const acb_mat_t Y1,
                                   slong d_maj, slong d_min, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CBNORM_H */
