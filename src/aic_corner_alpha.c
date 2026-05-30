/* aic_corner_alpha.c — the BLOCK BIJECTION alpha of lem_alpha (.tex:1086-1119),
 * Increment 2a of the corner module. The ASSEMBLY half: the alpha/beta build,
 * gamma = beta alpha - I, the CERTIFIED ||gamma||<1 contraction guard, and the
 * solve alpha^{-1} = (beta alpha)^{-1} beta. The ALLOCATION/DIM half
 * (build_blocks, free_blocks, aic_corner_alpha_dims, and the coords_in /
 * mid+radius certified-opnorm helpers) lives in the sibling
 * src/aic_corner_alpha_build.c (split for the ~200 LOC limit, Rule 10; shared
 * declarations in src/aic_corner_internal.h). This file was itself split from
 * src/aic_corner_product.c (which keeps the compressed product cdot + the unit
 * Ptilde). Builds on Increment 1 (aic_corner_Co, aic_corner_apply,
 * aic_corner_extract); see include/aic_corner.h for the API contracts.
 *
 * approximate_algebras.tex:1092-1119 — lem_alpha (verbatim, the maps):
 *   "alpha_{jk} : X |-> Co_{P,Q}(X) : S_{Pj,Qk} -> S_{P,Q},
 *    alpha = sum_{jk} alpha_{jk} : (+)_{jk} S_{Pj,Qk} -> S_{P,Q}.
 *    Then alpha is a linear bijection satisfying ||alpha|| <= pq + O(pq(d+e)),
 *    ||alpha^{-1}|| <= 1 + O(pq(d+e)) ... beta_{jk} : X |-> Co_{Pj,Qk}(X) :
 *    S_{P,Q} -> S_{Pj,Qk} ... beta alpha = 1 + gamma, ||gamma|| <= pq O(d+e) < 1,
 *    so beta alpha is invertible and (beta alpha)^{-1} beta is a left inverse."
 *
 * .tex TYPO (escalated, bead aic-czm NOTES). tex:1109 PRINTS
 *   "beta_{jk} : X |-> Co_{P_j,Q_j}(X)"  (Q_j),
 * but (a) the codomain is S_{Pj,Qk} and (b) the proof at tex:1114 needs
 *   beta_{jk} alpha_{lm}(X_{lm}) = delta_jl delta_km X_lm + O(d+e),
 * the Kronecker delta_km forcing the SECOND index to be Q_k. We implement
 * beta_{jk} = Co_{Pj,Qk} (Law 1: document, do not silently "fix" the math; the
 * justification is tex:1114). A Q_j build is mutation-tested RED (test_corner T9:
 * the Qj typo makes ||gamma|| jump from 0 to 1.0000, tripping the assert).
 *
 * CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). In finite dim each corner S_{Pj,Qk} and
 * S_{P,Q} is a concrete subspace of A with an extracted Frobenius-orthonormal
 * operator basis ({C^{jk}_m} and {D_l}). alpha/beta are then EXPLICIT matrices in
 * those coordinates: column ((j,k),m) of alpha is <D_l, Co_{P,Q}(C^{jk}_m)>_F;
 * column l of the beta_{jk}-block is <C^{jk}_m, Co_{Pj,Qk}(D_l)>_F. gamma = beta
 * alpha - I; the contraction ||gamma||_op < 1 (asserted, Rule 4) makes beta alpha
 * invertible and alpha^{-1} = (beta alpha)^{-1} beta (acb_mat_solve, certified).
 * The bijection forces N = sum dim S_{Pj,Qk} = d_PQ = dim S_{P,Q} (asserted; the
 * eta=0 dim-count oracle, tex:1124). Used only for p,q <= 2 (tex:1084).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner.h"
#include "aic_corner_internal.h"
#include "aic_ecstar.h"
#include "aic_mat.h"

void aic_corner_alpha(acb_mat_t alpha, acb_mat_t alpha_inv,
                      arb_t norm_alpha, arb_t norm_alpha_inv, arb_t gamma_norm,
                      arb_t gamma_rad,
                      slong *N_out, slong *dPQ_out,
                      const aic_ecstar *A,
                      const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                      const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                      slong prec)
{
    slong n = A->n, d = A->dim_A, pq = p * q;
    acb_mat_t **blocks, *Dbasis;
    slong *dimB, N, dPQ;
    aic_corner_build_blocks(&blocks, &dimB, &N, &Dbasis, &dPQ, A, P_parts, p, P,
                            Q_parts, q, Q, prec);
    *N_out = N;
    *dPQ_out = dPQ;
    /* bijection => N == d_PQ (the dim-count oracle, .tex:1124). Fail loud. */
    assert(N == dPQ && "alpha: dim S_{P,Q} != sum dim S_{Pj,Qk} (not a bijection)");
    assert(acb_mat_nrows(alpha) == dPQ && acb_mat_ncols(alpha) == N &&
           "alpha: alpha must be d_PQ x N");

    acb_mat_t CoPQ, W;
    acb_mat_init(CoPQ, d, d);
    acb_mat_init(W, n, n);
    aic_corner_Co(CoPQ, A, P, Q, prec);
    acb_ptr coordD = _acb_vec_init(dPQ);

    /* alpha[l, col] = <D_l, Co_{P,Q}(C^{jk}_m)>_F. beta[col, l] = <C^{jk}_m,
     * Co_{Pj,Qk}(D_l)>_F. Build alpha column-by-block and beta row-by-block. The
     * column/row offset of block b is the running sum of preceding dimB. */
    acb_mat_t beta;
    acb_mat_init(beta, N, dPQ);
    slong off = 0;
    for (slong j = 0; j < p; j++)
        for (slong k = 0; k < q; k++) {
            slong b = j * q + k;
            acb_mat_t CoBlock;
            acb_mat_init(CoBlock, d, d);
            aic_corner_Co(CoBlock, A, P_parts[j], Q_parts[k], prec); /* Co_{Pj,Qk}
                                                       NOT Co_{Pj,Qj}: tex TYPO */
            for (slong m = 0; m < dimB[b]; m++) {
                /* alpha column off+m: Co_{P,Q}(C^{jk}_m) in {D_l}. */
                aic_corner_apply(W, A, CoPQ, blocks[b][m], prec);
                aic_corner_coords_in(coordD, Dbasis, dPQ, W, prec);
                for (slong l = 0; l < dPQ; l++)
                    acb_set(acb_mat_entry(alpha, l, off + m), coordD + l);
            }
            for (slong l = 0; l < dPQ; l++) {
                /* beta rows off..off+dimB[b]: Co_{Pj,Qk}(D_l) in {C^{jk}_m}. */
                aic_corner_apply(W, A, CoBlock, Dbasis[l], prec);
                acb_ptr coordC = _acb_vec_init(dimB[b]);
                aic_corner_coords_in(coordC, blocks[b], dimB[b], W, prec);
                for (slong m = 0; m < dimB[b]; m++)
                    acb_set(acb_mat_entry(beta, off + m, l), coordC + m);
                _acb_vec_clear(coordC, dimB[b]);
            }
            acb_mat_clear(CoBlock);
            off += dimB[b];
        }

    /* BA = beta alpha (N x N); gamma = BA - I; assert ||gamma||_op < 1 (Rule 4).
     * The contraction decision is made on a CERTIFIED UPPER BOUND
     *   gnorm = ||mid(gamma)||_op + ||rad(gamma)||_F  >=  ||gamma||_op
     * (aic_corner_gamma_opnorm_ub, FIX 2): the bare midpoint discarded gamma's
     * radius, a soundness hole per the arb ladder ("a ball that has lost all
     * precision is a loud failure, not a silent widening"). Adding the Frobenius
     * norm of the entry-radius matrix restores the discarded width as a rigorous
     * inflation, so gnorm < 1 now certifies ||gamma||_op < 1 for EVERY matrix in
     * the ball, not just its midpoint. The radius contribution radc is surfaced
     * via the gamma_rad output (printed by the test; ~1e-72, negligible). The
     * midpoint route (not aic_mat_opnorm's Gram path on gamma) is forced by the
     * Gram off-diagonal false-fail on orthogonal-column matrices (bead aic-2yo);
     * the +||rad||_F makes the SUM certified regardless of that route's looseness. */
    acb_mat_t BA, gamma;
    acb_mat_init(BA, N, N);
    acb_mat_init(gamma, N, N);
    acb_mat_mul(BA, beta, alpha, prec);
    acb_mat_one(gamma);
    acb_mat_sub(gamma, BA, gamma, prec);          /* beta alpha - I_N */
    arb_t gnorm, radc;
    arb_init(gnorm);
    arb_init(radc);
    aic_corner_gamma_opnorm_ub(gnorm, radc, gamma, prec);
    {
        arb_t one;
        arb_init(one);
        arb_one(one);
        assert(arb_lt(gnorm, one) &&
               "alpha: certified ||beta alpha - I||_op upper bound >= 1 -- "
               "lem_alpha hypothesis pq(delta+eps)<const violated (overlapping "
               "projections / delta too large); beta alpha not invertible");
        arb_clear(one);
    }
    if (gamma_norm) arb_set(gamma_norm, gnorm);
    if (gamma_rad) arb_set(gamma_rad, radc);

    /* alpha^{-1} = (beta alpha)^{-1} beta (N x d_PQ), certified solve:
     * solve BA . alpha_inv = beta for alpha_inv. */
    if (alpha_inv) {
        assert(acb_mat_nrows(alpha_inv) == N && acb_mat_ncols(alpha_inv) == dPQ &&
               "alpha: alpha_inv must be N x d_PQ");
        int ok = acb_mat_solve(alpha_inv, BA, beta, prec);
        assert(ok && "alpha: acb_mat_solve(beta alpha) failed -- prec too low to "
                     "certify the inverse though ||gamma||<1");
    }
    /* ||alpha||, ||alpha^{-1}|| are INFORMATIONAL here (they gate nothing — only
     * gamma gates, above). We still report them as CERTIFIED upper bounds
     * (mid+radius) for consistency with the gamma guard; the radius term is
     * discarded (NULL rad_out). */
    if (norm_alpha) aic_corner_gamma_opnorm_ub(norm_alpha, NULL, alpha, prec);
    if (norm_alpha_inv) {
        if (alpha_inv) {
            aic_corner_gamma_opnorm_ub(norm_alpha_inv, NULL, alpha_inv, prec);
        } else {
            acb_mat_t ai;
            acb_mat_init(ai, N, dPQ);
            int ok = acb_mat_solve(ai, BA, beta, prec);
            assert(ok && "alpha: acb_mat_solve for norm_alpha_inv failed");
            aic_corner_gamma_opnorm_ub(norm_alpha_inv, NULL, ai, prec);
            acb_mat_clear(ai);
        }
    }

    arb_clear(radc);
    arb_clear(gnorm);
    acb_mat_clear(gamma);
    acb_mat_clear(BA);
    acb_mat_clear(beta);
    _acb_vec_clear(coordD, dPQ);
    acb_mat_clear(W);
    acb_mat_clear(CoPQ);
    aic_corner_free_blocks(blocks, dimB, pq, Dbasis, dPQ);
}
