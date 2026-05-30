/* aic_corner_alpha_build.c — the ALLOCATION/DIM half of the lem_alpha block
 * bijection (.tex:1086-1119), Increment 2a of the corner module. Split from
 * src/aic_corner_alpha.c (which keeps the alpha/beta assembly + the certified
 * contraction guard + the solve) to honor the ~200 LOC limit (Rule 10). The
 * companion declarations live in src/aic_corner_internal.h.
 *
 * This TU owns:
 *   - aic_corner_build_blocks / aic_corner_free_blocks : extract the per-block
 *     corner bases {C^{jk}} = S_{Pj,Qk} and the target basis {D_l} = S_{P,Q}
 *     (the shared core of both aic_corner_alpha and aic_corner_alpha_dims);
 *   - aic_corner_alpha_dims : report N = sum dim S_{Pj,Qk} and d_PQ = dim S_{P,Q}
 *     without assembling alpha (so the caller can allocate alpha to the right
 *     data-dependent shape, then call aic_corner_alpha);
 *   - aic_corner_coords_in : the Frobenius re-coordination helper;
 *   - aic_corner_gamma_opnorm_ub : the CERTIFIED upper bound on ||M||_op used by
 *     aic_corner_alpha's ||gamma||<1 contraction decision (FIX 2 — see the
 *     internal-header docstring and the inline note below).
 *
 * See src/aic_corner_alpha.c for the full lem_alpha statement, the .tex TYPO note
 * (beta_{jk} = Co_{Pj,Qk}, not the printed Co_{Pj,Qj}), and the constructive
 * route. include/aic_corner.h carries the public API contracts.
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/mag.h>

#include "aic_corner.h"
#include "aic_corner_internal.h"
#include "aic_ecstar.h"
#include "aic_mat.h"

/* CERTIFIED upper bound on ||M||_op = ||mid(M)||_op + ||rad(M)||_F (FIX 2). The
 * Rule-4 ||gamma||<1 contraction decision in aic_corner_alpha was previously made
 * on ||mid(gamma)||_op alone, discarding gamma's certified radius — a soundness
 * hole per the arb ladder ("a ball that has lost all precision is a loud failure,
 * not a silent widening"). We restore the discarded width as a rigorous inflation:
 *   ||M||_op <= ||mid(M)||_op + ||M - mid(M)||_op,
 *   ||M - mid(M)||_op <= ||M - mid(M)||_F <= ||rad(M)||_F,
 * since |M_ij - mid(M)_ij| <= rad(M)_ij (a complex entry's deviation is bounded by
 * hypot(rad Re, rad Im) <= rad Re + rad Im; we use the exact hypot via the squared
 * radii) and ||.||_op <= ||.||_F. ||rad(M)||_F^2 = sum_ij (rad Re_ij)^2 +
 * (rad Im_ij)^2. The midpoint norm uses aic_mat_opnorm on mid(M) (NOT the certified
 * Gram path on M, which false-fails on near-zero off-diagonals, bead aic-2yo); the
 * +||rad||_F term makes the SUM a genuinely certified upper bound regardless. */
void aic_corner_gamma_opnorm_ub(arb_t out, arb_t rad_out, const acb_mat_t M,
                                slong prec)
{
    slong r = acb_mat_nrows(M), c = acb_mat_ncols(M);

    /* term 1: ||mid(M)||_op via aic_mat_opnorm on the zero-radius midpoint. */
    acb_mat_t Mid;
    acb_mat_init(Mid, r, c);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(Mid, i, j), acb_mat_entry(M, i, j));
    arb_t mid_norm;
    arb_init(mid_norm);
    aic_mat_opnorm(mid_norm, Mid, prec);

    /* term 2: ||rad(M)||_F = sqrt(sum_ij rad(Re)_ij^2 + rad(Im)_ij^2). Accumulate
     * the sum-of-squares in a mag_t (rigorous upper-rounded magnitude arithmetic),
     * then convert to an arb upper bound and sqrt it. */
    mag_t sumsq, t;
    mag_init(sumsq);
    mag_init(t);
    mag_zero(sumsq);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++) {
            /* arb_radref already returns the entry radius AS a mag_t (an exact
             * upper bound), so copy it directly — no arf_get_mag needed. */
            mag_set(t, arb_radref(acb_realref(acb_mat_entry(M, i, j))));
            mag_addmul(sumsq, t, t);   /* sumsq += (rad Re)^2 (upper-rounded) */
            mag_set(t, arb_radref(acb_imagref(acb_mat_entry(M, i, j))));
            mag_addmul(sumsq, t, t);   /* sumsq += (rad Im)^2 */
        }
    mag_t fro_mag;
    mag_init(fro_mag);
    mag_sqrt(fro_mag, sumsq);          /* ||rad(M)||_F, upper-rounded (a mag is an
                                        * exact upper bound on the true value) */
    /* a zero-radius arb whose value is the mag's (upper-bound) value: arf_set_mag
     * is exact (mag's mantissa fits arf), so rad_bound = ||rad(M)||_F EXACTLY as a
     * point and stays a rigorous upper bound on ||M - mid(M)||_op. */
    arb_t rad_bound;
    arb_init(rad_bound);
    arf_set_mag(arb_midref(rad_bound), fro_mag);
    mag_zero(arb_radref(rad_bound));

    arb_add(out, mid_norm, rad_bound, prec);   /* certified upper bound */
    if (rad_out) arb_set(rad_out, rad_bound);

    arb_clear(rad_bound);
    mag_clear(fro_mag);
    mag_clear(t);
    mag_clear(sumsq);
    arb_clear(mid_norm);
    acb_mat_clear(Mid);
}

void aic_corner_coords_in(acb_ptr coord, acb_mat_t *Bset, slong m,
                          const acb_mat_t W, slong prec)
{
    for (slong i = 0; i < m; i++)
        aic_corner_frob_ip(coord + i, Bset[i], W, prec);
}

void aic_corner_build_blocks(acb_mat_t ***blocks, slong **dimB, slong *N,
                             acb_mat_t **Dbasis, slong *dPQ, const aic_ecstar *A,
                             const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                             const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                             slong prec)
{
    slong d = A->dim_A, pq = p * q;
    acb_mat_t Co;
    acb_mat_init(Co, d, d);

    /* target basis {D_l} of S_{P,Q}. */
    aic_corner_Co(Co, A, P, Q, prec);
    aic_corner_extract(Dbasis, dPQ, Co, A, prec);

    *blocks = flint_malloc((size_t) pq * sizeof(acb_mat_t *));
    *dimB = flint_malloc((size_t) pq * sizeof(slong));
    *N = 0;
    for (slong j = 0; j < p; j++)
        for (slong k = 0; k < q; k++) {
            slong b = j * q + k;
            aic_corner_Co(Co, A, P_parts[j], Q_parts[k], prec);  /* Co_{Pj,Qk} */
            aic_corner_extract(&(*blocks)[b], &(*dimB)[b], Co, A, prec);
            *N += (*dimB)[b];
        }
    acb_mat_clear(Co);
}

void aic_corner_free_blocks(acb_mat_t **blocks, slong *dimB, slong pq,
                            acb_mat_t *Dbasis, slong dPQ)
{
    for (slong b = 0; b < pq; b++) {
        for (slong m = 0; m < dimB[b]; m++) acb_mat_clear(blocks[b][m]);
        if (blocks[b]) flint_free(blocks[b]);
    }
    flint_free(blocks);
    flint_free(dimB);
    for (slong l = 0; l < dPQ; l++) acb_mat_clear(Dbasis[l]);
    if (Dbasis) flint_free(Dbasis);
}

void aic_corner_alpha_dims(slong *N_out, slong *dPQ_out, const aic_ecstar *A,
                           const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                           const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                           slong prec)
{
    acb_mat_t **blocks, *Dbasis;
    slong *dimB, N, dPQ;
    aic_corner_build_blocks(&blocks, &dimB, &N, &Dbasis, &dPQ, A, P_parts, p, P,
                            Q_parts, q, Q, prec);
    *N_out = N;
    *dPQ_out = dPQ;
    aic_corner_free_blocks(blocks, dimB, p * q, Dbasis, dPQ);
}
