/* aic_idemp_internal.h — private helpers shared between the idemp_structure
 * source files (bead aic-wuh). Not part of the public API (include/aic_idemp.h).
 * The double-path extraction (carrier split, image SVD) lives in
 * aic_idemp_carrier.c / aic_idemp_image.c and is wired together by
 * aic_idemp_maps.c; the verify routines are in aic_idemp_verify.c.
 *
 * Vec convention (see include/aic_idemp.h): vec(X)[i*n+j] = X[i,j] (row-major).
 */
#ifndef AIC_IDEMP_INTERNAL_H
#define AIC_IDEMP_INTERNAL_H

#include <flint/acb_mat.h>

#include "aic_idemp.h"
#include "aic_ucp.h"

/* Step 1 (.tex:2056, lem_carrier). Split the carrier operator Q's spectrum and
 * return dim_M. Because the rank is computed INSIDE, this routine acb_mat_INITs
 * its three outputs to the right shape (the caller passes uninitialised
 * acb_mat_t and must clear them): J_M (n x dim_M), J_Mperp (n x (n-dim_M), which
 * may be n x 0 when M = H), Pi_M (n x n). Double path (LAPACKE_zheev); fails
 * loud on a non-PSD or unresolved-rank Q. */
slong aic_idemp_carrier_split(acb_mat_t J_M, acb_mat_t J_Mperp, acb_mat_t Pi_M,
                              const acb_mat_t Q, slong n, slong prec);

/* Step 2 (.tex:2056, A = Img Phi). Extract an ONB {B_k} of A = Img Phi via the
 * column-image SVD: apply Phi to each matrix unit E_{ij}, stack vec(Phi(E_ij)) as
 * the columns of P (n^2 x n^2), SVD, keep left singular vectors with sigma > thr.
 * INITs Delta (n^2 x dim_A) with columns vec(B_k) (orthonormal) and returns
 * dim_A (caller passes uninitialised Delta and must clear it). Double path
 * (LAPACKE_zgesvd). Fails loud on an unresolved rank. */
slong aic_idemp_image_basis(acb_mat_t Delta, const aic_ucp_kraus *phi,
                            slong n, slong prec);

/* Steps 5-7 (.tex:2061-2084). Given out->n, out->dim_M, out->dim_A, out->J_M and
 * out->Delta already filled, INIT and fill out->Lambda, out->Gamma, out->w. */
void aic_idemp_build_maps(aic_idemp_decomp *out, const aic_ucp_kraus *phi,
                          slong prec);

#endif /* AIC_IDEMP_INTERNAL_H */
