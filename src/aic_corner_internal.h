/* aic_corner_internal.h — shared internal helpers for the corner module
 * (src/aic_corner_compress.c, src/aic_corner_extract.c, src/aic_corner_alpha.c,
 * src/aic_corner_alpha_build.c). NOT a public API.
 *
 * The Frobenius inner product <B_i, Y>_F = Tr(B_i^dag Y) = sum_{a,b}
 * conj(B_i[a,b]) Y[a,b] is used in BOTH translation units (the L_P/R_Q coordinate
 * builds in compress.c, and the operator->coords step of aic_corner_apply in
 * extract.c). Declared here so the single definition (in compress.c) is shared
 * rather than duplicated.
 *
 * The lem_alpha block bijection (.tex:1086-1119) is split across two TUs for the
 * ~200 LOC limit (Rule 10): the ALLOCATION/DIM half (build_blocks, free_blocks,
 * aic_corner_alpha_dims and the coords_in / mid_opnorm / gamma_opnorm_ub helpers)
 * lives in src/aic_corner_alpha_build.c; the ASSEMBLY half (aic_corner_alpha:
 * the alpha/beta build, gamma, the certified ||gamma||<1 guard and the solve)
 * stays in src/aic_corner_alpha.c. The helpers shared across the split are
 * declared below. */
#ifndef AIC_CORNER_INTERNAL_H
#define AIC_CORNER_INTERNAL_H

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_ecstar.h"

/* <Bi, Y>_F into `out` (caller-initialised acb_t). Bi, Y are n x n; n taken from
 * Y (Bi must match). */
void aic_corner_frob_ip(acb_t out, const acb_mat_t Bi, const acb_mat_t Y,
                        slong prec);

/* Frobenius coordinates of W (n x n) in the Frobenius-orthonormal operator basis
 * {Bset[0..m-1]}: coord[i] = <Bset[i], W>_F, into the length-m acb_ptr `coord`
 * (caller-allocated). Shared by the alpha/beta builds. */
void aic_corner_coords_in(acb_ptr coord, acb_mat_t *Bset, slong m,
                          const acb_mat_t W, slong prec);

/* CERTIFIED UPPER BOUND on ||M||_op (FIX 2, bead aic — the soundness item from the
 * Increment-2a hostile review): ub = ||mid(M)||_op + ||rad(M)||_F, where rad(M) is
 * the matrix of entry radii. Both terms are themselves rigorous upper bounds:
 *   ||M||_op <= ||mid(M)||_op + ||M - mid(M)||_op,  and  ||M - mid(M)||_op <=
 *   ||M - mid(M)||_F <= ||rad(M)||_F  (each |M_ij - mid(M)_ij| <= rad(M_ij);
 *   ||.||_op <= ||.||_F always).
 * So `out` >= ||M||_op for every matrix in the ball M, a genuinely certified
 * guard for the ||gamma||<1 contraction decision (the bare midpoint discarded the
 * radius — sound only because the radius was ~1e-72 in the tested regime). The
 * midpoint route is used (not aic_mat_opnorm's certified Gram path) because the
 * Gram off-diagonal relative-Hermiticity check false-fails on the near-zero
 * off-diagonals of an orthogonal-column M (bead aic-2yo); the +||rad||_F term
 * restores the discarded width as a rigorous inflation. Writes the radius
 * contribution to `rad_out` (caller-init'd arb_t, NULL to skip) so the caller can
 * report it (should be ~1e-72, negligible). `out`, `rad_out` caller-initialised. */
void aic_corner_gamma_opnorm_ub(arb_t out, arb_t rad_out, const acb_mat_t M,
                                slong prec);

/* Build the per-block corner bases {C^{jk}} (block index b = j*q + k) and the
 * target basis {D_l} of S_{P,Q}. Allocates *blocks (pq arrays, block b holding
 * dimB[b] operators), *dimB (slong[pq]), *Dbasis (dPQ operators); writes *N (=
 * sum dimB[b]) and *dPQ. Caller frees via aic_corner_free_blocks. Shared core of
 * aic_corner_alpha and aic_corner_alpha_dims. */
void aic_corner_build_blocks(acb_mat_t ***blocks, slong **dimB, slong *N,
                             acb_mat_t **Dbasis, slong *dPQ, const aic_ecstar *A,
                             const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                             const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                             slong prec);

void aic_corner_free_blocks(acb_mat_t **blocks, slong *dimB, slong pq,
                            acb_mat_t *Dbasis, slong dPQ);

#endif /* AIC_CORNER_INTERNAL_H */
