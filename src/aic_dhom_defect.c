/* aic_dhom_defect.c — the delta-homomorphism v : B -> A, its multiplicativity
 * defect G_v, the basis-sweep defect, and prop_delta_hominc (bead aic-c1n).
 * Realizes the defect definition (approximate_algebras.tex:1258) and the norm
 * bounds prop_delta_hominc (.tex:1194-1222).
 *
 * THE STAR IS LOAD-BEARING (CLAUDE.md domain callout). G_v(X,Y) = v(XY) - v(X)*v(Y)
 * uses TWO different products: XY is B's BLOCK product (aic_dhom_B_mul, exact,
 * since B is a genuine C* algebra), but v(X)*v(Y) is the Choi-Effros STAR of A
 * (aic_ecstar_star = Phi_tilde(v(X) v(Y))), NEVER plain matmul (.tex:341, :2189).
 * Plain matmul leaves A = Img Phi_tilde. The eta=0 identity-channel oracle has
 * star == plain product and CANNOT catch a star-vs-plain bug; test_dhom T4 mutates
 * the star on a genuine eta>0 A and confirms the bound is violated (RED).
 *
 * v IS LINEAR IN B-COORDS. v(X) = sum_i x_i v(E_i) where x_i are the matrix-unit
 * entries of the block-diagonal X (read block by block, .tex data model). So
 * aic_dhom_v_apply reads the n_B x n_B entries of X at the matrix-unit linear
 * indices and accumulates the stored vE[i].
 *
 * prop_delta_hominc — three predicates (.tex:1194). The norms a = inf, b = ||v||
 * are the BASIS-SWEEP estimates (max/min over the matrix-unit basis), a documented
 * LOWER/UPPER-on-the-basis surrogate for the true inf/sup over the operator-norm
 * unit ball (the faithful HOPM search is a later cycle; cf. aic_ecstar). The unit
 * defect ||v(I_B) - I_A||_op is computed exactly (no surrogate).
 *
 * No eigendecomposition: products + a CERTIFIED operator-norm UPPER BOUND
 * (aic_corner_gamma_opnorm_ub, the mid+radius bound), used instead of
 * aic_mat_opnorm because the latter's Gram-route Hermiticity check false-fails on
 * the near-zero off-diagonals of the orthogonal-column-ish defect matrices here
 * (bead aic-qgs/aic-2yo). NOTE: this makes norm_lb a min-of-upper-bounds, an even
 * looser surrogate for the true inf a; documented in aic_dhom_prop_bounds.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner_internal.h" /* aic_corner_gamma_opnorm_ub (certified UB) */
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_mat.h"

/* CERTIFIED operator-norm UPPER BOUND, dodging the aic_mat_opnorm Gram-path
 * Hermiticity false-fail on near-zero off-diagonals (bead aic-qgs/aic-2yo). Used
 * for all star-defect / difference norms here, where the Gram of an
 * orthogonal-column-ish matrix trips the relative-Hermiticity check. The bound is
 * ||mid(M)||_op + ||rad(M)||_F >= ||M||_op (see aic_corner_internal.h). */
static void opnorm_ub(arb_t out, const acb_mat_t M, slong prec)
{
    aic_corner_gamma_opnorm_ub(out, NULL, M, prec);
}

void aic_dhom_v_init(aic_dhom_v *v, const aic_dhom_B *B, const aic_ecstar *A)
{
    assert(v != NULL && B != NULL && A != NULL);
    v->B = B;
    v->A = A;
    v->n = A->n;
    v->vE = flint_malloc((size_t) B->dim_B * sizeof(acb_mat_t));
    for (slong i = 0; i < B->dim_B; i++) {
        acb_mat_init(v->vE[i], v->n, v->n);
        acb_mat_zero(v->vE[i]);
    }
}

void aic_dhom_v_clear(aic_dhom_v *v)
{
    if (v == NULL || v->vE == NULL) return;
    for (slong i = 0; i < v->B->dim_B; i++) acb_mat_clear(v->vE[i]);
    flint_free(v->vE);
    v->vE = NULL;
}

/* block row-offset of block l in the n_B x n_B representative. */
static slong block_row_offset(const aic_dhom_B *B, slong l)
{
    slong r = 0;
    for (slong t = 0; t < l; t++) r += B->d[t];
    return r;
}

void aic_dhom_v_apply(acb_mat_t out, const aic_dhom_v *v, const acb_mat_t X,
                      slong prec)
{
    const aic_dhom_B *B = v->B;
    assert(acb_mat_nrows(X) == B->n_B && acb_mat_ncols(X) == B->n_B);
    assert(acb_mat_nrows(out) == v->n && acb_mat_ncols(out) == v->n);
    acb_mat_zero(out);
    acb_mat_t scaled;
    acb_mat_init(scaled, v->n, v->n);
    for (slong l = 0; l < B->num_blocks; l++) {
        slong d = B->d[l], r0 = block_row_offset(B, l);
        for (slong a = 0; a < d; a++)
            for (slong b = 0; b < d; b++) {
                slong i = aic_dhom_B_index(B, l, a, b);
                acb_srcptr coef = acb_mat_entry(X, r0 + a, r0 + b);
                if (acb_is_zero(coef)) continue;
                acb_mat_scalar_mul_acb(scaled, v->vE[i], coef, prec);
                acb_mat_add(out, out, scaled, prec);
            }
    }
    acb_mat_clear(scaled);
}

void aic_dhom_Gv(acb_mat_t out, const aic_dhom_v *v, const acb_mat_t X,
                 const acb_mat_t Y, slong prec)
{
    /* approximate_algebras.tex:1258  G_v(X,Y) = v(XY) - v(X) * v(Y).
     * XY = B's block product; v(X)*v(Y) = STAR in A (NEVER plain matmul). */
    slong nB = v->B->n_B, n = v->n;
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n);
    acb_mat_t XY, vXY, vX, vY, vXvY;
    acb_mat_init(XY, nB, nB);
    acb_mat_init(vXY, n, n);
    acb_mat_init(vX, n, n);
    acb_mat_init(vY, n, n);
    acb_mat_init(vXvY, n, n);
    aic_dhom_B_mul(XY, v->B, X, Y, prec);     /* XY in B (block)            */
    aic_dhom_v_apply(vXY, v, XY, prec);       /* v(XY)                      */
    aic_dhom_v_apply(vX, v, X, prec);
    aic_dhom_v_apply(vY, v, Y, prec);
    aic_ecstar_star(vXvY, v->A, vX, vY, prec); /* v(X) * v(Y) (STAR)        */
    acb_mat_sub(out, vXY, vXvY, prec);
    acb_mat_clear(vXvY);
    acb_mat_clear(vY);
    acb_mat_clear(vX);
    acb_mat_clear(vXY);
    acb_mat_clear(XY);
}

void aic_dhom_defect_sweep(arb_t out, const aic_dhom_v *v, slong prec)
{
    slong dim_B = v->B->dim_B, nB = v->B->n_B, n = v->n;
    acb_mat_t Ei, Ej, g;
    arb_t e;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(Ej, nB, nB);
    acb_mat_init(g, n, n);
    arb_init(e);
    arb_zero(out);
    for (slong i = 0; i < dim_B; i++) {
        aic_dhom_B_matunit(Ei, v->B, i);
        for (slong j = 0; j < dim_B; j++) {
            aic_dhom_B_matunit(Ej, v->B, j);
            aic_dhom_Gv(g, v, Ei, Ej, prec);
            opnorm_ub(e, g, prec);
            if (arb_gt(e, out)) arb_set(out, e);
        }
    }
    arb_clear(e);
    acb_mat_clear(g);
    acb_mat_clear(Ej);
    acb_mat_clear(Ei);
}

void aic_dhom_prop_bounds(arb_t norm_ub, arb_t norm_lb, arb_t unit_def,
                          const aic_dhom_v *v, slong prec)
{
    /* approximate_algebras.tex:1194: (i) ||v|| <= 1+O; (ii) a >= 1-O if a>2delta;
     * (iii) ||v(I_B)-I_A|| <= O. Basis-sweep a = min_i ||vE[i]||, b = max_i. */
    slong dim_B = v->B->dim_B, n = v->n;
    arb_t e;
    arb_init(e);
    if (norm_ub != NULL || norm_lb != NULL) {
        int first = 1;
        for (slong i = 0; i < dim_B; i++) {
            opnorm_ub(e, v->vE[i], prec);
            if (norm_ub != NULL && (first || arb_gt(e, norm_ub)))
                arb_set(norm_ub, e);
            if (norm_lb != NULL && (first || arb_lt(e, norm_lb)))
                arb_set(norm_lb, e);
            first = 0;
        }
    }
    if (unit_def != NULL) {
        /* ||v(I_B) - I_A||_op. I_A = the ambient identity 1_n (A's unit, .tex:2186). */
        acb_mat_t IB, vIB, IA;
        acb_mat_init(IB, v->B->n_B, v->B->n_B);
        acb_mat_init(vIB, n, n);
        acb_mat_init(IA, n, n);
        aic_dhom_B_unit(IB, v->B);
        aic_dhom_v_apply(vIB, v, IB, prec);
        acb_mat_one(IA);
        acb_mat_sub(vIB, vIB, IA, prec);
        opnorm_ub(unit_def, vIB, prec);
        acb_mat_clear(IA);
        acb_mat_clear(vIB);
        acb_mat_clear(IB);
    }
    arb_clear(e);
}
