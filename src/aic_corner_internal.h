/* aic_corner_internal.h — shared internal helper for the corner module
 * (src/aic_corner_compress.c, src/aic_corner_extract.c). NOT a public API.
 *
 * The Frobenius inner product <B_i, Y>_F = Tr(B_i^dag Y) = sum_{a,b}
 * conj(B_i[a,b]) Y[a,b] is used in BOTH translation units (the L_P/R_Q coordinate
 * builds in compress.c, and the operator->coords step of aic_corner_apply in
 * extract.c). Declared here so the single definition (in compress.c) is shared
 * rather than duplicated. */
#ifndef AIC_CORNER_INTERNAL_H
#define AIC_CORNER_INTERNAL_H

#include <flint/acb.h>
#include <flint/acb_mat.h>

/* <Bi, Y>_F into `out` (caller-initialised acb_t). Bi, Y are n x n; n taken from
 * Y (Bi must match). */
void aic_corner_frob_ip(acb_t out, const acb_mat_t Bi, const acb_mat_t Y,
                        slong prec);

#endif /* AIC_CORNER_INTERNAL_H */
