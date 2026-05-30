/* aic_dhom_B.c — the genuine finite-dim C* algebra B = (+)_l M_{d_l}: its
 * matrix-unit basis, block product, and unit (bead aic-c1n, module "dhom"). The
 * EXPLICIT generalized-Pauli diagonal of B (the dimension-independence crux) lives
 * in the sibling TU aic_dhom_diag.c (split for the ~200 LOC limit, Rule 10); the
 * shared block-row-offset helper is declared in aic_dhom_internal.h.
 *
 * B is a GENUINE C* algebra: the block product is EXACT block matmul, the
 * involution is the block adjoint, the unit is I_B = (+)_l I_{d_l}, and the norm is
 * the operator norm of the n_B x n_B (n_B = sum_l d_l) block-diagonal
 * representative. An element X with matrix-unit coords x^{(l)}_{ab} is stored as
 * that block-diagonal matrix; the linear index mu(l,a,b) = offset_l + a*d_l + b
 * (offset_l = sum_{l'<l} d_{l'}^2) addresses the dim_B = sum_l d_l^2 coordinates.
 *
 * NO eigendecomposition: pure entry fills + block matmul.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_dhom.h"
#include "aic_dhom_internal.h"

void aic_dhom_B_init(aic_dhom_B *B, const slong *dims, slong m)
{
    assert(B != NULL && dims != NULL && m >= 1);
    B->num_blocks = m;
    B->d = flint_malloc((size_t) m * sizeof(slong));
    B->offset = flint_malloc((size_t) m * sizeof(slong));
    slong off = 0, nb = 0;
    for (slong l = 0; l < m; l++) {
        assert(dims[l] >= 1);
        B->d[l] = dims[l];
        B->offset[l] = off;
        off += dims[l] * dims[l];
        nb += dims[l];
    }
    B->dim_B = off;
    B->n_B = nb;
}

void aic_dhom_B_clear(aic_dhom_B *B)
{
    if (B == NULL) return;
    flint_free(B->d);
    flint_free(B->offset);
    B->d = NULL;
    B->offset = NULL;
}

slong aic_dhom_B_index(const aic_dhom_B *B, slong l, slong a, slong b)
{
    assert(l >= 0 && l < B->num_blocks);
    assert(a >= 0 && a < B->d[l] && b >= 0 && b < B->d[l]);
    return B->offset[l] + a * B->d[l] + b;
}

/* Row-offset of block l in the n_B x n_B block-diagonal representative. Shared
 * with aic_dhom_diag.c (declared in aic_dhom_internal.h). */
slong aic_dhom_block_row_offset(const aic_dhom_B *B, slong l)
{
    slong r = 0;
    for (slong t = 0; t < l; t++) r += B->d[t];
    return r;
}

void aic_dhom_B_mul(acb_mat_t Z, const aic_dhom_B *B, const acb_mat_t X,
                    const acb_mat_t Y, slong prec)
{
    slong nb = B->n_B;
    assert(acb_mat_nrows(Z) == nb && acb_mat_ncols(Z) == nb);
    assert(acb_mat_nrows(X) == nb && acb_mat_ncols(X) == nb);
    assert(acb_mat_nrows(Y) == nb && acb_mat_ncols(Y) == nb);
    assert(Z != X && Z != Y);
    acb_mat_zero(Z);
    for (slong l = 0; l < B->num_blocks; l++) {
        slong d = B->d[l], r0 = aic_dhom_block_row_offset(B, l);
        acb_mat_t Xb, Yb, Zb;
        acb_mat_init(Xb, d, d);
        acb_mat_init(Yb, d, d);
        acb_mat_init(Zb, d, d);
        for (slong a = 0; a < d; a++)
            for (slong b = 0; b < d; b++) {
                acb_set(acb_mat_entry(Xb, a, b), acb_mat_entry(X, r0 + a, r0 + b));
                acb_set(acb_mat_entry(Yb, a, b), acb_mat_entry(Y, r0 + a, r0 + b));
            }
        acb_mat_mul(Zb, Xb, Yb, prec);
        for (slong a = 0; a < d; a++)
            for (slong b = 0; b < d; b++)
                acb_set(acb_mat_entry(Z, r0 + a, r0 + b), acb_mat_entry(Zb, a, b));
        acb_mat_clear(Zb);
        acb_mat_clear(Yb);
        acb_mat_clear(Xb);
    }
}

void aic_dhom_B_unit(acb_mat_t out, const aic_dhom_B *B)
{
    assert(acb_mat_nrows(out) == B->n_B && acb_mat_ncols(out) == B->n_B);
    acb_mat_zero(out);
    for (slong l = 0; l < B->num_blocks; l++) {
        slong d = B->d[l], r0 = aic_dhom_block_row_offset(B, l);
        for (slong a = 0; a < d; a++)
            acb_set_si(acb_mat_entry(out, r0 + a, r0 + a), 1);
    }
}

void aic_dhom_B_matunit(acb_mat_t out, const aic_dhom_B *B, slong i)
{
    assert(i >= 0 && i < B->dim_B);
    assert(acb_mat_nrows(out) == B->n_B && acb_mat_ncols(out) == B->n_B);
    /* find block l and (a,b) from linear index i */
    slong l = 0;
    while (l + 1 < B->num_blocks && B->offset[l + 1] <= i) l++;
    slong d = B->d[l], rem = i - B->offset[l];
    slong a = rem / d, b = rem % d, r0 = aic_dhom_block_row_offset(B, l);
    acb_mat_zero(out);
    acb_set_si(acb_mat_entry(out, r0 + a, r0 + b), 1);
}
