/* aic_cbnorm_distance.c — eig-free certified two-sided bracket on the GENERAL
 * cb-norm / diamond distance ||Phi - Psi||_cb between two UCP self-maps (bead
 * aic-l5b, increment A). See include/aic/aic_cbnorm.h for the .tex provenance
 * (approximate_algebras.tex:347-354) and the public contract.
 *
 * WHY THIS FILE EXISTS / WHY IT IS PURE GLUE. The distance ||Phi - Psi||_cb is
 * the same kind of object as the eta-idempotence defect ||Phi^2 - Phi||_cb that
 * aic_cbnorm.c already certifies: in BOTH cases the argument of the cb-norm is a
 * Hermiticity-preserving (NOT CP) map Lambda. For the defect, Lambda = Phi^2 -
 * Phi; here, Lambda = Phi - Psi, a difference of two Hermiticity-preserving UCP
 * maps, hence itself Hermiticity-preserving. So NONE of the math changes — the
 * whole bracket DERIVATION lives in aic_cbnorm.c's docstring and is NOT re-
 * derived here; we only state why Lambda = Phi - Psi qualifies and reuse it.
 *
 * THE MATH (the spec; the bracket is the one derived in aic_cbnorm.c).
 *   approximate_algebras.tex:347-352 (verbatim, the cb-norm):
 *     ||Lambda||_cb = sup_n sup_{X!=0} ||(1_{M_n} (x) Lambda)(X)|| / ||X||,
 *     "efficiently computable using semidefinite programming"; the dual norm
 *     for maps of functionals is the diamond norm ||.||_diamond.
 *   approximate_algebras.tex:354 (verbatim): "A UCP map Phi: B(H)->B(H) is
 *     called eta-idempotent if ||Phi^2-Phi||_cb <= eta."
 * For a Hermiticity-preserving map cb-norm = diamond-norm and the Watrous 2012
 * ampliation truncation N = n (= input dim) is rigorous. With
 *   J = Choi(Phi) - Choi(Psi) = Choi(Phi - Psi) = Choi(Lambda)
 * (linearity of the Choi map; built by aic_ucp_choi_diff), the SAME derivation
 * recorded in aic_cbnorm.c gives the certified bracket
 *   ||Phi - Psi||_cb  in  [ ||J||_F / n ,  2 * ||J||_F ]   (ratio hi/lo = 2n).
 *
 * IMPLEMENTATION. Build J via the existing aic_ucp_choi_diff (Convention A,
 * n^2 x n^2, Hermitian) and defer to the existing aic_cbnorm_eigfree_ball_choi,
 * which returns the certified [||J||_F/n, 2*||J||_F] arb balls. NO new bound, NO
 * new numerics: this is glue over two already-certified primitives. The only
 * asserts are fail-loud SHAPE/caller-error guards (Rule 4); the numerics never
 * abort. Phi == Psi (J = 0) gives [0,0] (the cleanest oracle).
 *
 * Number path: arb/acb ONLY (NO LAPACK/latd) — this rung must stay certified,
 * matching aic_cbnorm.c.
 */
#include <assert.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_cbnorm.h"
#include "aic/aic_ucp.h"

/* approximate_algebras.tex:347-352 — cb-norm (= diamond-norm for Herm-pres.):
 *   ||Lambda||_cb = sup_n sup_{X!=0} ||(1_{M_n} (x) Lambda)(X)|| / ||X||,
 *   "efficiently computable using semidefinite programming".
 * approximate_algebras.tex:354 — Phi is eta-idempotent if ||Phi^2-Phi||_cb<=eta.
 * Lambda = Phi - Psi is Hermiticity-preserving (difference of two Herm-pres. UCP
 * maps), JUST LIKE Phi^2 - Phi, so the bracket [||J||_F/n, 2*||J||_F] derived in
 * aic_cbnorm.c's docstring applies verbatim with J = Choi(Phi - Psi). */
void aic_cbnorm_eigfree_distance(arb_t lo, arb_t hi, const aic_ucp_kraus *phi,
                                 const aic_ucp_kraus *psi, slong prec)
{
    /* Both maps must be self-maps sharing the same dim (Rule 4, fail loud —
     * caller-error guard; aic_ucp_choi_diff also requires equal (dim_K,dim_H)). */
    assert(phi->dim_K == phi->dim_H); /* phi a self-map */
    assert(psi->dim_K == psi->dim_H); /* psi a self-map */
    assert(phi->dim_H == psi->dim_H); /* same n => all four dims equal */
    slong n = phi->dim_H;
    slong N = n * n;

    acb_mat_t J;
    acb_mat_init(J, N, N);
    aic_ucp_choi_diff(J, phi, psi, prec); /* J = Choi(Phi) - Choi(Psi) = Choi(Lambda) */

    aic_cbnorm_eigfree_ball_choi(lo, hi, J, n, prec);

    acb_mat_clear(J);
}
