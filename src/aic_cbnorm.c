/* aic_cbnorm.c — eig-free certified two-sided bracket on the eta-idempotence
 * defect eta = ||Phi^2 - Phi||_cb = ||Phi^2 - Phi||_diamond (bead aic-m24,
 * increment 1). See include/aic_cbnorm.h for the .tex provenance
 * (approximate_algebras.tex:347-354) and the role of this module (the always-
 * valid fallback rung of the certified-bound ladder; the tight ball is a later
 * increment).
 *
 * THE DERIVATION (Law 1 — recorded in full so the bound survives a rewrite).
 * Let Lambda = Phi^2 - Phi (Hermiticity-preserving, not CP), J = Choi(Lambda)
 * in Convention-A (n^2 x n^2, Hermitian; aic_ucp.h header). The Watrous 2012 SDP
 * (ALGORITHM.md cb-norm section, tools/gen_fixtures_d24.jl) gives
 *     eta = (2/n) * optval,
 *     optval = max Re tr(J^dag X)  s.t.  [[P, X],[X^dag, Q]] >= 0,
 *                                        P + Q = I_{n^2},  P,Q >= 0.
 *
 * RIGOROUS UPPER BOUND (eig-free). On the feasible set P + Q = I_{n^2} with
 * P,Q >= 0 forces 0 <= P,Q <= I (so ||P||_op, ||Q||_op <= 1) and
 * tr(Q) <= tr(P+Q) = tr(I_{n^2}) = n^2. The 2x2-block PSD condition gives the
 * Schur factorization X = P^{1/2} K Q^{1/2} with ||K||_op <= 1. Hence
 *     ||X||_F <= ||P^{1/2}||_op * ||K||_op * ||Q^{1/2}||_F
 *             <= 1 * 1 * sqrt(tr Q) <= sqrt(n^2) = n.
 * (Used ||A B||_F <= ||A||_op ||B||_F and ||Q^{1/2}||_F = sqrt(tr Q).) Then by
 * Cauchy-Schwarz Re tr(J^dag X) <= ||J||_F * ||X||_F <= n * ||J||_F, so
 *     optval <= n * ||J||_F   =>   eta = (2/n)*optval <= 2 * ||J||_F.
 *
 * RIGOROUS LOWER BOUND (eig-free). Convention-A J = (id_n (x) Lambda)(omega),
 * omega = |Omega><Omega|, |Omega> = sum_i |i>|i> (UNNORMALIZED). For the
 * NORMALIZED maximally-entangled state rho = omega/n, (id (x) Lambda)(rho) = J/n.
 * The diamond norm is a sup over input states, and the truncation N = n is
 * rigorous for a Hermiticity-preserving map (Watrous 2012; bead aic-d24), so
 *     eta = ||Lambda||_diamond >= ||(id (x) Lambda)(rho)||_1 = ||J||_1 / n
 *                              >= ||J||_2 / n = ||J||_F / n
 * (last step: Schatten ||.||_1 >= ||.||_2 = ||.||_F).
 *
 * So the certified bracket is  [ ||J||_F / n ,  2 * ||J||_F ].  Ratio hi/lo = 2n
 * (valid but loose — expected; tightening is a later increment, NOT done here).
 *
 * IMPLEMENTATION. fro = ||J||_F via aic_mat_frobenius_norm (already a certified
 * arb ball). lo = fro / n (arb_div_ui), hi = 2*fro (arb_mul_2exp_si by +1).
 * Returning the arb BALLS is sufficient for rigor: the true ||J||_F lies inside
 * the fro ball, so by monotonicity arb_lower(lo) <= ||J||_F/n <= eta and
 * eta <= 2*||J||_F <= arb_upper(hi). NO eigendecomposition, NO SDP, NEVER aborts
 * (the only assert is a fail-loud SHAPE check, Rule 4). eta=0 (J=0) gives [0,0].
 *
 * Number path: arb/acb only (NO LAPACK/latd) — this rung must stay certified.
 */
#include <assert.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_cbnorm.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

/* approximate_algebras.tex:347-352 — cb-norm (= diamond-norm for Herm-pres.):
 *   ||Lambda||_cb = sup_n sup_{X!=0} ||(1_{M_n} (x) Lambda)(X)|| / ||X||,
 *   "efficiently computable using semidefinite programming".
 * approximate_algebras.tex:354 — Phi is eta-idempotent if ||Phi^2-Phi||_cb<=eta.
 * Bracket [||J||_F/n, 2*||J||_F] derived in the file docstring above. */
void aic_cbnorm_eigfree_ball_choi(arb_t lo, arb_t hi, const acb_mat_t J,
                                  slong n, slong prec)
{
    slong N = n * n;
    /* Fail loud on a shape mismatch (Rule 4): J must be the n^2 x n^2 Choi. */
    assert(acb_mat_nrows(J) == N && acb_mat_ncols(J) == N);
    assert(n >= 1);

    arb_t fro;
    arb_init(fro);
    aic_mat_frobenius_norm(fro, J, prec); /* certified ball enclosing ||J||_F */

    /* lo = ||J||_F / n  (lower bound: eta >= ||J||_1/n >= ||J||_F/n). */
    arb_div_ui(lo, fro, (ulong) n, prec);

    /* hi = 2 * ||J||_F  (upper bound: eta = (2/n)*optval <= 2*||J||_F). */
    arb_mul_2exp_si(hi, fro, 1); /* multiply by 2^1 = 2, exact (no rounding) */

    arb_clear(fro);
}

/* Wrapper: bracket straight from a UCP self-map Phi (dim_K == dim_H == n).
 * Builds Phi^2 = compose(Phi,Phi) and J = Choi(Phi^2) - Choi(Phi) = Choi(Lambda)
 * (linearity of Choi), then defers to the choi-level routine. */
void aic_cbnorm_eigfree_ball(arb_t lo, arb_t hi, const aic_ucp_kraus *phi,
                             slong prec)
{
    assert(phi->dim_K == phi->dim_H); /* self-map required (Rule 4) */
    slong n = phi->dim_H;
    slong N = n * n;

    aic_ucp_kraus phi2;
    aic_ucp_compose(&phi2, phi, phi, prec); /* Phi^2, n x n, r^2 Kraus ops */

    acb_mat_t J;
    acb_mat_init(J, N, N);
    aic_ucp_choi_diff(J, &phi2, phi, prec); /* Choi(Phi^2 - Phi) */

    aic_cbnorm_eigfree_ball_choi(lo, hi, J, n, prec);

    acb_mat_clear(J);
    aic_ucp_kraus_clear(&phi2);
}
