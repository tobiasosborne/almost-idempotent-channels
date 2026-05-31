/* aic_ecstar_setup.h — internal interface for the HOPM double-snapshot + scratch
 * + defect-map term thunks (bead aic-knm, Cycle 2). Implemented in
 * aic_ecstar_setup.c; used by the public driver aic_ecstar_certify.c. Not public. */
#ifndef AIC_ECSTAR_SETUP_H
#define AIC_ECSTAR_SETUP_H

#include "aic/aic_ecstar.h"
#include "aic_ecstar_search.h"

/* double snapshot of an aic_ecstar (Kraus + basis -> double buffers, midpoints). */
void aic_ehk_snapshot(aic_ehk *h, const aic_ecstar *A);
void aic_ehk_snapshot_free(aic_ehk *h);

/* the six defect-map term thunks (assoc X/Y/Z, submult X/Y, cstar). */
void aic_ehk_term_assoc_X(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k);
void aic_ehk_term_assoc_Y(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k);
void aic_ehk_term_assoc_Z(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k);
void aic_ehk_term_sub_X(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k);
void aic_ehk_term_sub_Y(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k);
void aic_ehk_term_cstar(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k);

#endif /* AIC_ECSTAR_SETUP_H */
