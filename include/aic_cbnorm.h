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

#ifdef __cplusplus
}
#endif

#endif /* AIC_CBNORM_H */
