/* aic_opspace_map_internal.h — internal interface for the v-specific opmap
 * builders of the opspace module (bead aic-zwo). NOT a public header (lives in
 * src/): shares the forward/inverse opmap builders + the coordinate-matrix
 * inverter between aic_opspace_map.c (the builders) and aic_opspace_entry.c (the
 * public entry points + the factorize-interface adds), so each file stays under
 * the ~200 LOC core limit (CLAUDE.md Rule 10). See include/aic_opspace.h for the
 * public contract and FINDINGS §C12 (operator-norm not Frobenius).
 */
#ifndef AIC_OPSPACE_MAP_INTERNAL_H
#define AIC_OPSPACE_MAP_INTERNAL_H

#include "aic/aic_dhom.h"
#include "aic_opspace_hopm.h"

/* Forward opmap f = v: B -> A (U = B matrix units, W = {B_k}, F = M_1). */
void aic_opspace_opmap_forward(opmap *m, const aic_dhom_v *v, slong prec);

/* Inverse opmap f = v^{-1}: A -> B (U = {B_k}, W = B matrix units, F = M_1^{-1}).
 * ASSERTS v bijective (dim_A == dim_B). */
void aic_opspace_opmap_inverse(opmap *m, const aic_dhom_v *v, slong prec);

/* F = M_1^{-1} (d x d, row-major double) into `out` (caller-allocated d*d):
 * build M_1, invert via the FLINT acb solver, assert invertibility (Rule 4).
 * The single inversion point shared by opmap_inverse and aic_opspace_build_vinv. */
void aic_opspace_build_M1inv(cx *out, const aic_dhom_v *v, slong d, slong prec);

/* Matrix-unit linear index i -> global (row,col) in the n_B x n_B block-diagonal
 * representative of B (block offsets are the n_B-representative offsets, NOT the
 * dim_B matrix-unit offsets). Shared with the Choi assembler (aic_opspace_choi.c,
 * bead aic-pjr) so the (i -> (r,c)) map is single-sourced (assert on out-of-range). */
void aic_opspace_b_unit_pos(const aic_dhom_B *B, slong i, slong *row, slong *col);

#endif /* AIC_OPSPACE_MAP_INTERNAL_H */
