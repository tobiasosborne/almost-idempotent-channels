/* aic_ecstar_internal.h — shared helper for the ecstar defect estimators (bead
 * aic-knm). The estimators are split across aic_ecstar_defect.c (submult, C*,
 * involution, unit) and aic_ecstar_assoc.c (the trilinear associator) to keep
 * each file <= 200 LOC (Rule 10); this header holds the one helper both use.
 * static inline => no extra TU, no linker symbol. */
#ifndef AIC_ECSTAR_INTERNAL_H
#define AIC_ECSTAR_INTERNAL_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_mat.h"

/* worst = max(worst, ||M||_op) (certified arb). */
static inline void aic_ecstar_bump_opnorm(arb_t worst, const acb_mat_t M,
                                          slong prec)
{
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    if (arb_gt(e, worst)) arb_set(worst, e);
    arb_clear(e);
}

#endif /* AIC_ECSTAR_INTERNAL_H */
