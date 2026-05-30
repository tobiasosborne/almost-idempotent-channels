/* aic_dhom_internal.h — shared internal helpers for the dhom module
 * (src/aic_dhom_B.c, src/aic_dhom_diag.c, src/aic_dhom_defect.c,
 * src/aic_dhom_approx.c). NOT a public API.
 *
 * The block-row-offset of block l in the n_B x n_B block-diagonal representative
 * of B = (+)_l M_{d_l} (= sum_{l'<l} d_{l'}) is used by both the B-algebra TU
 * (B_mul/B_unit/B_matunit) and the Pauli-diagonal TU (embed_block). The module is
 * split across those two TUs for the ~200 LOC limit (Rule 10), so the single
 * definition (in aic_dhom_B.c) is shared via this header rather than duplicated.
 */
#ifndef AIC_DHOM_INTERNAL_H
#define AIC_DHOM_INTERNAL_H

#include "aic_dhom.h"

/* Row-offset of block l in the n_B x n_B block-diagonal representative:
 * sum_{l'<l} d[l']. */
slong aic_dhom_block_row_offset(const aic_dhom_B *B, slong l);

#endif /* AIC_DHOM_INTERNAL_H */
