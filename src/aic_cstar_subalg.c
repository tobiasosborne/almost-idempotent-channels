/* aic_cstar_subalg.c — the S_P-as-aic_ecstar wrapper (bead aic-097, module
 * "cstar_build", Increment 1, piece A). Re-presents a compressed subalgebra S_P of
 * a parent eps-C* algebra A as a DERIVED aic_ecstar so that projection / corner /
 * dhom recurse on it. See include/aic_cstar.h for the full why; design contract in
 * docs/research/cstar_build_design.md §1.
 *
 * approximate_algebras.tex:1077-1082 — the compressed product compr_prod, and the
 * unit of (S_P, .):
 *   (X,Y) |-> X . Y : S_{P,Q} x S_{Q,R} -> S_{P,R},   X . Y = Co_{P,R}(XY)
 *   "If P is a nonvanishing projection, the compressed product turns the Banach
 *    space S_P (closed under X |-> X^dag) into an O(delta+eps)-C* algebra with
 *    unit Ptilde = Co_P(P)."
 * With P = Q = R = P this is the star of S_P: X * Y = Co_P(Phi_tilde(XY)).
 *
 * THE COMPRESSED STAR via the star_phi SEAM (mirrors aic_assoc_ecstar_from_phi).
 * The derived ecstar's star_phi thunk sp_star computes, for the matrix XY passed
 * by aic_ecstar_star,
 *   out = Co_P( Phi_tilde_parent(XY) )
 *       = aic_corner_apply(parent, Co_P, parent_Phi(XY))
 * where parent_Phi(XY) = aic_ecstar_star(parent, XY, I_n) = parent_Phi(XY . I_n)
 * = parent_Phi(XY) (the parent's star with the ambient unit, .tex:2211). This is
 * EXACTLY the compressed product aic_corner_cdot(parent, Co_P, X, Y) when XY is
 * formed from X, Y in S_P; tests/test_cstar.c asserts the thunk-vs-cdot equality.
 *
 * THE STAR, NOT THE PLAIN PRODUCT (CLAUDE.md callout, FINDINGS §C2/§C3,
 * LOAD-BEARING). S_P's product goes through Co_P and the parent's Phi_tilde, NOT
 * plain acb_mat_mul and NOT the Frobenius projector. The eta=0 identity channel
 * has Phi_tilde=id so its star == plain and is BLIND to a star bug; the oblique
 * eta>0 fixture (test T2) is what gives the star teeth.
 *
 * LIFETIME / RELOCATION-SAFETY (mirrors aic_assoc_ecstar's SEPARATE heap ctx).
 * The derived ecstar's star_ctx points at a SEPARATELY HEAP-ALLOCATED block
 * (aic_cstar_sp_ctx, holding the borrowed parent + an owned Co_P), NOT into this
 * owner struct. The old design set star_ctx = &owner (a self-pointer), which
 * DANGLES the moment the owner is copied by value or relocated inside a realloc'd
 * array — exactly what the I5 master loop does (it builds many wrappers in arrays).
 * With the heap ctx, a shallow copy / relocation of the owner leaves star_ctx
 * still pointing at the valid heap block. aic_ecstar_clear frees only the basis;
 * this owner additionally frees the heap ctx (and its Co_P) and Ptilde. The owner
 * is MOVE-ONLY: a by-value copy shares the ctx/Ptilde/basis, so clear it exactly
 * once (see include/aic_cstar.h). The parent is borrowed throughout.
 */
#include <assert.h>

#include <flint/acb_mat.h>

#include "aic_corner.h"
#include "aic_cstar.h"
#include "aic_ecstar.h"

/* The compressed-product star thunk. ctx is the HEAP aic_cstar_sp_ctx* (the
 * parent + its Co_P), NOT the owner struct — so the wrapper survives relocation.
 *   out = Co_P( parent_Phi(XY) ),  XY the ordinary product aic_ecstar_star formed.
 * parent_Phi(XY) = aic_ecstar_star(parent, XY, I_n) (the parent's star with the
 * ambient unit gives Phi_tilde(XY . I_n) = Phi_tilde(XY), .tex:2211). */
static void sp_star(acb_mat_t out, const acb_mat_t XY, void *ctx, slong prec)
{
    const aic_cstar_sp_ctx *c = (const aic_cstar_sp_ctx *) ctx;
    slong n = c->parent->n;
    acb_mat_t Imat, phiXY;
    acb_mat_init(Imat, n, n);
    acb_mat_init(phiXY, n, n);
    acb_mat_one(Imat);
    /* parent_Phi(XY) = Phi_tilde(XY . I_n) = Phi_tilde(XY). */
    aic_ecstar_star(phiXY, c->parent, XY, Imat, prec);
    /* out = Co_P(parent_Phi(XY)) (out != phiXY, the apply requires out != X). */
    aic_corner_apply(out, c->parent, c->Co_P, phiXY, prec);
    acb_mat_clear(phiXY);
    acb_mat_clear(Imat);
}

void aic_cstar_subalg_build(aic_cstar_subalgebra *out, const aic_ecstar *parent,
                            const acb_mat_t P, slong prec)
{
    assert(out != NULL && parent != NULL);
    slong n = parent->n, d = parent->dim_A;
    assert(acb_mat_nrows(P) == n && acb_mat_ncols(P) == n &&
           "aic_cstar_subalg_build: P must be n x n (n = parent->n)");

    out->parent = parent;

    /* Heap ctx (the F2 fix): a SEPARATE block holding the borrowed parent and an
     * OWNED Co_P = aic_corner_Co(parent, P, P) (the d x d compression idempotent).
     * star_ctx will point HERE, not at &out, so the owner survives relocation. */
    aic_cstar_sp_ctx *ctx = flint_malloc(sizeof(aic_cstar_sp_ctx));
    ctx->parent = parent;
    acb_mat_init(ctx->Co_P, d, d);
    aic_corner_Co(ctx->Co_P, parent, P, P, prec);
    out->ctx = ctx;

    /* {C_m}, dim_S = aic_corner_extract(Co_P): the Frobenius-orthonormal n x n
     * basis of S_P (top-dim_S left singular vectors of the oblique Co_P). */
    acb_mat_t *C;
    slong dim_S;
    aic_corner_extract(&C, &dim_S, ctx->Co_P, parent, prec);
    /* fail loud (Rule 4): dim S_P == 0 means P ~ 0 — there is no algebra. */
    assert(dim_S >= 1 &&
           "aic_cstar_subalg_build: dim S_P == 0 (P ~ 0); no subalgebra");

    /* The derived ecstar. aic_ecstar_init requires a non-NULL Kraus phi for the
     * Kraus path; the star_phi superop path does NOT (same as
     * aic_assoc_ecstar_from_phi), so fill the fields directly and adopt the
     * flint_malloc'd C array as the basis (aic_ecstar_clear frees it). */
    out->A.n = n;
    out->A.dim_A = dim_S;
    out->A.phi = NULL;
    out->A.star_phi = sp_star;
    out->A.star_ctx = out->ctx;     /* HEAP ctx (not &out): relocation-safe (F2) */
    out->A.B = C;                   /* adopt the extracted basis array */

    /* Ptilde = Co_P(P) = the unit of (S_P, .) (.tex:1082). */
    acb_mat_init(out->Ptilde, n, n);
    aic_corner_Ptilde(out->Ptilde, parent, P, prec);
}

void aic_cstar_subalg_clear(aic_cstar_subalgebra *out)
{
    if (out == NULL) return;
    aic_ecstar_clear(&out->A);     /* frees the basis only (heap star_ctx) */
    if (out->ctx != NULL) {
        acb_mat_clear(out->ctx->Co_P);
        flint_free(out->ctx);
        out->ctx = NULL;
    }
    acb_mat_clear(out->Ptilde);
    out->parent = NULL;
}
