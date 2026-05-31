/* aic_factorize.c — the non-UCP "tilde" maps Delta~, Upsilon~ of th_factorization
 * (bead aic-tff, module "factorize", Increment F1). See include/aic_factorize.h
 * for the full contract and docs/research/factorize_d4_research.md §1/§6.
 *
 * approximate_algebras.tex:2749-2752 — th_factorization (outline, verbatim):
 *   "By Theorem th_main_ext, there exist a finite-dimensional C* algebra B and an
 *    extended O(eta)-isomorphism v: B -> A. Let Delta~ be v followed by the
 *    inclusion A -> B(H), and let Upsilon~ be Phi~ with the target space A,
 *    followed by v^{-1}. These maps are not UCP but meet the other requirements.
 *    Indeed, it is immediate that
 *        Delta~ Upsilon~ = Phi~,   Upsilon~ Delta~ = 1_B,
 *        ||Delta~||_cb <= 1+O(eta),  ||Upsilon~||_cb <= 1+O(eta)."     (tilde_DelUps)
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The paper's v / v^{-1}
 * / Phi~ are abstract existence objects; in finite dim each is an explicit matrix
 * built and certified by a closed module (v: aic_cstar_build; v^{-1}:
 * aic_opspace_build_vinv, the SAME inverse the inverse stretch certifies; Phi~ =
 * S_tilde: aic_assoc_regularize). F1 is pure linear-algebra glue; no new numerical
 * primitive. The paper's exact identities (.tex:2751) are the spec; the tests
 * MEASURE them (the v-inversion + S_tilde idempotence precision is the floor).
 *
 * THE OBLIQUE PROJECTION IS LOAD-BEARING (FINDINGS §C2/§C3). Upsilon~ projects W
 * onto A = Img Phi~ with Phi~ ITSELF (the superop S_tilde), NOT the Frobenius
 * Pi_A(W) = sum_k <B_k,W>_F B_k: Phi~ is Hermicity-preserving but NOT
 * HS-self-adjoint, so Pi_A leaves the star defect O(1) while Phi~ gives the O(eta)
 * defect (§C3). At eta=0 Phi~ = Phi = id-on-A, so a dropped-Phi~ bug is INVISIBLE
 * on the eta=0 oracle (§C2) — only the eta>0 T2 fixture catches it (the recorded
 * mutation tooth). The v^{-1} contraction sum_k <B_k, Phi~(W)>_F vinvB[k] then uses
 * the SAME Frobenius-ON basis {B_k} = Aec->A.B[k] in which v^{-1} was built, so the
 * Frobenius inner product is the correct coordinate read of v^{-1} on A.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_assoc.h"
#include "aic_factorize.h"
#include "aic_opspace.h"

void aic_factorize_build(aic_factorize *out, const aic_dhom_v *v,
                         const aic_assoc_ecstar *Aec, slong prec)
{
    assert(out != NULL && v != NULL && Aec != NULL);
    /* The iso's codomain IS Aec's A (same object), and v must be bijective so the
     * coordinate matrix M_1 is square + invertible (Rule 4, fail loud). */
    assert(v->A == &Aec->A &&
           "aic_factorize_build: v->A must be the SAME object as &Aec->A");
    assert(v->n == Aec->A.n &&
           "aic_factorize_build: v->n != Aec->A.n");
    assert(v->A->dim_A == v->B->dim_B &&
           "aic_factorize_build: v not bijective (dim_A != dim_B)");

    out->N = v->n;
    out->n_B = v->B->n_B;
    out->v = v;
    out->Aec = Aec;

    slong dimA = 0, nB = 0;
    aic_opspace_build_vinv(&out->vinvB, &dimA, &nB, v, prec);
    assert(dimA == v->A->dim_A && nB == out->n_B &&
           "aic_factorize_build: vinvB shape mismatch");
    out->dim_A = dimA;
}

void aic_factorize_clear(aic_factorize *out)
{
    if (out == NULL) return;
    aic_opspace_vinv_clear(out->vinvB, out->dim_A);
    out->vinvB = NULL;
    out->v = NULL;
    out->Aec = NULL;
    out->N = out->n_B = out->dim_A = 0;
}

/* Delta~(X) = iota(v(X)) = v(X) (.tex:2749). iota: A -> B(H) is a no-op on the
 * N x N matrix; v(X) = sum_i x_i vE[i] is exactly aic_dhom_v_apply. */
void aic_factorize_delta_tilde(acb_mat_t out, const aic_factorize *F,
                               const acb_mat_t X, slong prec)
{
    assert(F != NULL);
    assert(acb_mat_nrows(out) == F->N && acb_mat_ncols(out) == F->N &&
           "delta_tilde: out must be N x N");
    assert(acb_mat_nrows(X) == F->n_B && acb_mat_ncols(X) == F->n_B &&
           "delta_tilde: X must be n_B x n_B (B's block-diagonal rep)");
    aic_dhom_v_apply(out, F->v, X, prec);
}

/* Upsilon~(W) = v^{-1}(Phi~(W)) = sum_k <B_k, Phi~(W)>_F vinvB[k] (.tex:2749).
 * Phi~(W) is the OBLIQUE projection onto A via the superop S_tilde (FINDINGS §C3);
 * the contraction reads its v^{-1} coordinates in A's Frobenius-ON basis {B_k}. */
void aic_factorize_upsilon_tilde(acb_mat_t out, const aic_factorize *F,
                                 const acb_mat_t W, slong prec)
{
    assert(F != NULL);
    assert(acb_mat_nrows(out) == F->n_B && acb_mat_ncols(out) == F->n_B &&
           "upsilon_tilde: out must be n_B x n_B (B's block-diagonal rep)");
    assert(acb_mat_nrows(W) == F->N && acb_mat_ncols(W) == F->N &&
           "upsilon_tilde: W must be N x N");

    slong N = F->N;
    const aic_ecstar *A = &F->Aec->A;

    /* Phi~(W): apply the idempotent superoperator S_tilde (the OBLIQUE projector
     * onto A = Img Phi~, NOT the Frobenius Pi_A — §C3). */
    acb_mat_t PW;
    acb_mat_init(PW, N, N);
    aic_assoc_superop_apply(PW, F->Aec->S_tilde, W, prec);

    /* out = sum_k <B_k, Phi~(W)>_F vinvB[k], <B_k,Z>_F = sum_ij conj(B_k[i,j]) Z[i,j]. */
    acb_mat_zero(out);
    acb_t coef, z;
    acb_init(coef);
    acb_init(z);
    acb_mat_t scaled;
    acb_mat_init(scaled, F->n_B, F->n_B);
    for (slong k = 0; k < F->dim_A; k++) {
        acb_zero(coef);
        for (slong i = 0; i < N; i++)
            for (slong j = 0; j < N; j++) {
                acb_conj(z, acb_mat_entry(A->B[k], i, j));
                acb_mul(z, z, acb_mat_entry(PW, i, j), prec);
                acb_add(coef, coef, z, prec);
            }
        acb_mat_scalar_mul_acb(scaled, F->vinvB[k], coef, prec);
        acb_mat_add(out, out, scaled, prec);
    }
    acb_mat_clear(scaled);
    acb_clear(z);
    acb_clear(coef);
    acb_mat_clear(PW);
}
