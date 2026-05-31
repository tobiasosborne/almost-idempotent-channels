/* aic_projection_internal.h — shared internals for the projection module (bead
 * aic-mqf), split across src/aic_projection.c (orchestrator + Steps 3-5) and
 * src/aic_projection_find.c (Steps 1-2: H-search + gap). NOT a public API
 * (include/aic_projection.h). Split for the ~200 LOC/file limit (Rule 10).
 */
#ifndef AIC_PROJECTION_INTERNAL_H
#define AIC_PROJECTION_INTERNAL_H

#include <flint/acb_mat.h>

#include "aic/aic_ecstar.h"

/* H_k = (B_k + B_k^dag)/2 into `out` (caller-init'd n x n). */
void aic_projection_herm_part(acb_mat_t out, const acb_mat_t Bk, slong prec);

/* Scan A's basis for the H_k of largest INTERIOR SPECTRAL GAP (tie-break by
 * spread); write the chosen H into Hout (caller-init'd n x n) and its index to
 * *k_out. FAILS LOUD if all near-scalar. Returns the chosen H's largest gap. */
double aic_projection_pick_H(acb_mat_t Hout, slong *k_out, const aic_ecstar *A,
                             slong prec);

/* Largest interior spectral gap of Hermitian H (double-path zheev eigenvalues):
 * threshold t -> t_out, gap -> g_out, lam_min and lam_max -> lmin, lmax, the
 * +1-cluster size m -> m_out. FAILS LOUD if no interior gap >= g_min. */
void aic_projection_gap(double *t_out, double *g_out, double *lmin, double *lmax,
                        slong *m_out, const acb_mat_t H, slong prec);

#endif /* AIC_PROJECTION_INTERNAL_H */
