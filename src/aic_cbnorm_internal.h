/* aic_cbnorm_internal.h — internal split of the TIGHT cb-norm certifier (bead
 * aic-m24, increment 3b). The public entry point aic_cbnorm_certify lives in
 * src/aic_cbnorm_certify.c (dispatch + fail-loud); the two feasibility-
 * restoration recipes (the bulk of the math) live in
 * src/aic_cbnorm_certify_restore.c. Split to keep each file <=200 LOC
 * (CLAUDE.md Rule 10). NOT a public header — the API is include/aic_cbnorm.h.
 */
#ifndef AIC_CBNORM_INTERNAL_H
#define AIC_CBNORM_INTERNAL_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

/* LOWER bound (convex-combination feasibility restoration; design doc, LOWER).
 * Writes the rigorous lower-bound ball into `lo` (arb_lower(lo) <= eta) and the
 * restoration defect into `delta`. Inputs: J and the MAX-primal point (X,P,Q),
 * each n^2 x n^2. */
void aic_cbnorm_int_lower(arb_t lo, arb_t delta, const acb_mat_t J,
                          const acb_mat_t X, const acb_mat_t P,
                          const acb_mat_t Q, slong n, slong prec);

/* UPPER bound (eigenvalue-shift feasibility restoration; design doc, UPPER).
 * Writes the rigorous upper-bound ball into `hi` (eta <= arb_upper(hi)) and the
 * restoration shift into `eps`. Inputs: J and the MIN-dual point (Y0,Y1), each
 * n^2 x n^2. */
void aic_cbnorm_int_upper(arb_t hi, arb_t eps, const acb_mat_t J,
                          const acb_mat_t Y0, const acb_mat_t Y1,
                          slong n, slong prec);

#endif /* AIC_CBNORM_INTERNAL_H */
