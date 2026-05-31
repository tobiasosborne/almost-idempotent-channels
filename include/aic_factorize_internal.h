/* aic_factorize_internal.h — PRIVATE internal seam shared between the F3 source
 * files aic_factorize_upsilon2.c (the build) and aic_factorize_upsilon3.c (the
 * apply maps + the D5 pin). NOT part of the public API; do not include from tests
 * or downstream modules. Exists only because build_Lj is needed by BOTH the
 * upsilon build (caches F->L[j]) and the D5 pin (rebuilds L_j with the alternate
 * (x)1_F ordering), and Rule 10 (<=200 LOC/file) forces them into separate TUs.
 */
#ifndef AIC_FACTORIZE_INTERNAL_H
#define AIC_FACTORIZE_INTERNAL_H

#include <flint/acb_mat.h>

#include "aic_factorize.h"

#ifdef __cplusplus
extern "C" {
#endif

/* L_j = sum_s p_{js} ((x)1_F-factor of Delta(U_{js}^dag)) V W_j^dag (U_{js}(x)xi_j)
 * (.tex:2861, D5). `f_left` selects the (x)1_F ordering: 1 = F-LEFT (CORRECT,
 * matches V's F-major Stinespring layout, 1_F (x) Delta(U^dag)); 0 = F-RIGHT (the
 * WRONG ordering, used ONLY by the D5 pin). `Lj` caller-init'd (N*r) x d_j. */
void aic_factorize_int_build_Lj(acb_mat_t Lj, const aic_factorize *F, slong j,
                                const acb_mat_t Wj, const acb_mat_t xi, slong e_j,
                                int f_left, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_FACTORIZE_INTERNAL_H */
