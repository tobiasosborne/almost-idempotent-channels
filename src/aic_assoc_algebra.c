/* aic_assoc_algebra.c — the associated eps-C* algebra A = Img Phi_tilde with the
 * approximate Choi-Effros star X*Y = Phi_tilde(XY) (bead aic-92f, module
 * "assoc_ecsa", Increment 2). Wires the regularization (Increment 1) + the
 * A-extraction (aic_assoc_extract.c) into the ecstar verification data model
 * (include/aic_ecstar.h) via its generic star_phi seam. See include/aic_assoc.h.
 *
 * approximate_algebras.tex:2184-2195 — th_almost_idemp:
 *     A = Img Phi_tilde = Ker(1 - Phi_tilde);  X * Y = Phi_tilde(XY);
 *     (A, *, ||.||, dag, I) is an extended O(eta)-C* algebra.
 *
 * THE STAR via the SUPEROP. Phi_tilde is an n^2 x n^2 matrix S_tilde (Increment 1),
 * NOT a UCP/Kraus map (.tex:363 it is not even CP). The ecstar star is generalised
 * through star_phi: aic_ecstar_star computes XY then star_phi(out, XY, ctx). Here
 * star_phi = star_apply (below) applies S_tilde via the tested aic_assoc_superop_apply.
 * out->A.phi stays NULL (no Kraus map exists), out->A.star_phi = star_apply.
 *
 * LIFETIME. The ecstar BORROWS star_phi/star_ctx (aic_ecstar_clear frees only the
 * basis), so this owner holds S_tilde + the ctx and frees them in
 * aic_assoc_ecstar_clear -- mirroring ecstar's borrowed-phi rule.
 *
 * THE MIDPOINT-S_tilde DECISION (load-bearing, NOT a bandaid). The regularized
 * S_tilde carries accumulated funcalc-iteration error balls (~5e-75 radius at
 * prec=256, from Newton-Schulz). The certified rigor of the regularization itself
 * is established in Increment 1 (test T3 eta=0 oracle machine-zero; T6 the
 * non-normal A1-vs-A2 cross-check to ~1e-74) -- that ladder rung does NOT live
 * inside the defect estimators. Inside the estimators, the star applies S_tilde
 * to a product XY and the result is fed to aic_mat_opnorm, which forms the Gram
 * M^dag M; squaring a ~5e-75 radius into the Gram exceeds aic_mat_opnorm's strict
 * Hermiticity tol 2^-(prec-8) ~ 2.2e-75 at prec=256 (the asymmetry is an artifact
 * of acb_mat_mul not knowing M^dag M is Hermitian, NOT real non-Hermiticity). So
 * the algebra is built from the MIDPOINT of S_tilde (a zero-radius idempotent
 * superop): the defect estimators then return certified balls of a well-defined
 * exact-idempotent superop, exactly as the basis B_k is a zero-radius double-path
 * extraction. This mirrors Increment 1's midpoint_opnorm for its trend measurements
 * (ALGORITHM.md) and the project's "extraction is double, the DEFECT checks are
 * arb-certified" stance. Tightening this to a radius-aware opnorm is the same
 * deferral as the inflated-radius wall (Increment 1 ALGORITHM.md, bead aic-w4o.1).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic/aic_assoc.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_ucp.h"

/* ctx for the star thunk: a borrowed pointer to the owner's S_tilde. */
typedef struct { const acb_mat_struct *S; } assoc_star_ctx;

/* star_phi: out = Phi_tilde(X) = reshape(S_tilde vec_r(X)) (the Choi-Effros
 * Phi_tilde applied; aic_assoc_superop_apply is the tested superop-apply). */
static void star_apply(acb_mat_t out, const acb_mat_t X, void *ctx, slong prec)
{
    const assoc_star_ctx *c = (const assoc_star_ctx *) ctx;
    aic_assoc_superop_apply(out, c->S, X, prec);
}

void aic_assoc_ecstar_from_phi(aic_assoc_ecstar *out, const aic_ucp_kraus *phi,
                               slong prec)
{
    assert(out != NULL && phi != NULL);
    slong n = phi->dim_H;
    assert(phi->dim_K == n &&
           "aic_assoc_ecstar_from_phi: phi must be a self-map (dim_K==dim_H)");

    /* (a) S_tilde = theta(2 Phi - 1) (asserts the prop_P basin inside, Rule 4),
     * then strip the funcalc-iteration error balls to the midpoint (file
     * docstring "THE MIDPOINT-S_tilde DECISION"). The certified rigor of the
     * regularization is Increment 1's T3/T6; the defect estimators then operate
     * on a clean exact-idempotent superop. */
    acb_mat_init(out->S_tilde, n * n, n * n);
    aic_assoc_regularize(out->S_tilde, phi, prec);
    for (slong i = 0; i < n * n; i++)
        for (slong j = 0; j < n * n; j++)
            acb_get_mid(acb_mat_entry(out->S_tilde, i, j),
                        acb_mat_entry(out->S_tilde, i, j));

    /* (b) A's basis {B_k} = range(S_tilde). */
    acb_mat_t *B;
    slong dim_A;
    aic_assoc_extract_range(&B, &dim_A, out->S_tilde, prec);

    /* (c) the ctx wrapping S_tilde + the ecstar with star_phi = star_apply. */
    assoc_star_ctx *c = flint_malloc(sizeof(assoc_star_ctx));
    c->S = out->S_tilde;
    out->ctx = c;

    /* aic_ecstar_init requires a non-NULL phi for the Kraus path; the superop
     * path does NOT, so we allocate the ecstar's fields directly here (init's
     * phi-assert is for the Kraus path, per the architecture decision). */
    out->A.n = n;
    out->A.dim_A = dim_A;
    out->A.phi = NULL;
    out->A.star_phi = star_apply;
    out->A.star_ctx = out->ctx;
    out->A.B = flint_malloc((size_t) dim_A * sizeof(acb_mat_t));
    for (slong k = 0; k < dim_A; k++) {
        acb_mat_init(out->A.B[k], n, n);
        acb_mat_set(out->A.B[k], B[k]);
        acb_mat_clear(B[k]);
    }
    flint_free(B);
}

void aic_assoc_ecstar_clear(aic_assoc_ecstar *out)
{
    if (out == NULL) return;
    aic_ecstar_clear(&out->A);     /* frees the basis only (borrowed star_phi) */
    if (out->ctx != NULL) flint_free(out->ctx);
    out->ctx = NULL;
    acb_mat_clear(out->S_tilde);
}
