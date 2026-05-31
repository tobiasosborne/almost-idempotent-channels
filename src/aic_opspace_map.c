/* aic_opspace_map.c — the v-specific opmap BUILDERS of the opspace operator-norm
 * ampliation (bead aic-zwo, increment O1). Builds the level-1 `opmap`
 * (src/aic_opspace_hopm.h) for the FORWARD map v: B -> A and the INVERSE
 * v^{-1}: A -> B from an aic_dhom_v, plus the coordinate-matrix inverter. The
 * public entry points (the stretch HOPM drivers + the factorize-interface adds —
 * the structural ampliation and the v^{-1} builder) are in aic_opspace_entry.c;
 * the HOPM kernel is in aic_opspace_ampliate.c. Split three ways to keep each file
 * under the ~200 LOC core limit (CLAUDE.md Rule 10). See include/aic_opspace.h for
 * the contract and FINDINGS §C12 (operator-norm not Frobenius); the shared builder
 * interface is src/aic_opspace_map_internal.h.
 *
 * THE TWO MAPS.
 *   forward v: B -> A: U = B (ON basis = matrix units {E_i}, sU = n_B), V = A (ON
 *     basis {B_k}, sV = N), F = M_1 (the dim_A x dim_B coordinate matrix of
 *     aic_dhom_v_sigma_min, M_1[k,i] = <B_k, vE[i]>_F).
 *   inverse v^{-1}: A -> B: U = A, V = B, F = M_1^{-1} (square + invertible in the
 *     bijective case; built by the FLINT acb solver, assert on singularity).
 */
#include <assert.h>
#include <complex.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_latd.h"
#include "aic_opspace.h"
#include "aic_opspace_map_internal.h"

/* Matrix-unit linear index i -> global (row,col) in the n_B x n_B block-diagonal
 * representative of B (block offsets are the n_B-representative offsets, NOT the
 * dim_B matrix-unit offsets). Non-static + namespaced (declared in the internal
 * header) so the Choi assembler (src/aic_opspace_choi.c, bead aic-pjr) reuses the
 * SAME matrix-unit -> ambient-position map rather than reimplementing it. */
void aic_opspace_b_unit_pos(const aic_dhom_B *B, slong i, slong *row, slong *col)
{
    slong boff = 0;
    for (slong l = 0; l < B->num_blocks; l++) {
        slong base = B->offset[l], dl = B->d[l];
        if (i >= base && i < base + dl * dl) {
            slong loc = i - base;
            *row = boff + loc / dl;
            *col = boff + loc % dl;
            return;
        }
        boff += dl;
    }
    assert(0 && "b_unit_pos: index out of range");
}

/* M_1[k,i] = <B_k, vE[i]>_F as a row-major double dim_A x dim_B matrix. */
static void build_M1(cx *M1, const aic_dhom_v *v, slong prec)
{
    const aic_ecstar *A = v->A;
    slong dimA = A->dim_A, dimB = v->B->dim_B, N = v->n;
    acb_t z, acc;
    acb_init(z);
    acb_init(acc);
    for (slong k = 0; k < dimA; k++)
        for (slong i = 0; i < dimB; i++) {
            acb_zero(acc);
            for (slong p = 0; p < N; p++)
                for (slong q = 0; q < N; q++) {
                    acb_conj(z, acb_mat_entry(A->B[k], p, q));
                    acb_mul(z, z, acb_mat_entry(v->vE[i], p, q), prec);
                    acb_add(acc, acc, z, prec);
                }
            M1[k * dimB + i]
                = arf_get_d(arb_midref(acb_realref(acc)), ARF_RND_NEAR)
                + arf_get_d(arb_midref(acb_imagref(acc)), ARF_RND_NEAR) * I;
        }
    acb_clear(acc);
    acb_clear(z);
}

/* Forward opmap f = v: B -> A. U = B (matrix units, sU = n_B), W = {B_k} (sV = N),
 * F = M_1. */
void aic_opspace_opmap_forward(opmap *m, const aic_dhom_v *v, slong prec)
{
    const aic_ecstar *A = v->A;
    slong dimA = A->dim_A, dimB = v->B->dim_B, nB = v->B->n_B, N = v->n;
    m->dU = dimB; m->dV = dimA; m->sU = nB; m->sV = N;
    m->U = flint_malloc((size_t) dimB * sizeof(cx *));
    m->W = flint_malloc((size_t) dimA * sizeof(cx *));
    for (slong i = 0; i < dimB; i++) {
        m->U[i] = aic_opmap_cxalloc(nB * nB);
        for (slong p = 0; p < nB * nB; p++) m->U[i][p] = 0.0;
        slong row, col;
        aic_opspace_b_unit_pos(v->B, i, &row, &col);
        m->U[i][row * nB + col] = 1.0;
    }
    for (slong k = 0; k < dimA; k++) {
        m->W[k] = aic_opmap_cxalloc(N * N);
        aic_latd_from_acb_mat(m->W[k], A->B[k]);
    }
    m->F = aic_opmap_cxalloc(dimA * dimB);
    build_M1(m->F, v, prec);
}

/* F = M_1^{-1} (d x d, row-major double) into `out` (caller-allocated d*d): build
 * M_1, invert via the FLINT acb solver, assert invertibility (fail loud, Rule 4).
 * The single inversion point shared by opmap_inverse and aic_opspace_build_vinv. */
void aic_opspace_build_M1inv(cx *out, const aic_dhom_v *v, slong d, slong prec)
{
    cx *M1 = aic_opmap_cxalloc(d * d);
    build_M1(M1, v, prec);
    acb_mat_t Macb, Minv, Id;
    acb_mat_init(Macb, d, d);
    acb_mat_init(Minv, d, d);
    acb_mat_init(Id, d, d);
    aic_latd_to_acb_mat(Macb, M1, d, d);
    acb_mat_one(Id);
    int ok = acb_mat_solve(Minv, Macb, Id, prec);
    assert(ok && "build_M1inv: coordinate matrix M_1 not invertible (v singular)");
    aic_latd_from_acb_mat(out, Minv);
    acb_mat_clear(Id);
    acb_mat_clear(Minv);
    acb_mat_clear(Macb);
    flint_free(M1);
}

/* Inverse opmap f = v^{-1}: A -> B. U = {B_k} (sU = N), W = B (matrix units,
 * sV = n_B), F = M_1^{-1}. ASSERTS v bijective (dim_A == dim_B, square M_1). */
void aic_opspace_opmap_inverse(opmap *m, const aic_dhom_v *v, slong prec)
{
    const aic_ecstar *A = v->A;
    slong dimA = A->dim_A, dimB = v->B->dim_B, nB = v->B->n_B, N = v->n;
    assert(dimA == dimB && "opmap_inverse: v not bijective (dim_A != dim_B)");
    slong d = dimA;
    m->dU = dimA; m->dV = dimB; m->sU = N; m->sV = nB;
    m->U = flint_malloc((size_t) dimA * sizeof(cx *));
    m->W = flint_malloc((size_t) dimB * sizeof(cx *));
    for (slong k = 0; k < dimA; k++) {
        m->U[k] = aic_opmap_cxalloc(N * N);
        aic_latd_from_acb_mat(m->U[k], A->B[k]);
    }
    for (slong i = 0; i < dimB; i++) {
        m->W[i] = aic_opmap_cxalloc(nB * nB);
        for (slong p = 0; p < nB * nB; p++) m->W[i][p] = 0.0;
        slong row, col;
        aic_opspace_b_unit_pos(v->B, i, &row, &col);
        m->W[i][row * nB + col] = 1.0;
    }
    m->F = aic_opmap_cxalloc(d * d);
    aic_opspace_build_M1inv(m->F, v, d, prec);
}
