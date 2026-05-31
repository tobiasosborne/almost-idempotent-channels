/* aic_factorize.h — the approximate factorization of an almost-idempotent UCP
 * channel through a genuine C* algebra (bead aic-tff, module "factorize"),
 * Increment F1. Realizes the NON-UCP "tilde" maps of th_factorization
 * (approximate_algebras.tex:2730-2899), specifically the outline at .tex:2749-2752.
 * Design / feasibility: docs/research/factorize_d4_research.md §1 (Steps 1-3) + §6
 * (the F1-F4 increment skeleton — this file is F1). FINDINGS §D4 (BUILDABLE-MODULO),
 * §C2 (the star), §C3 (A = Img Phi~ is OBLIQUE).
 *
 * WHAT F1 REALIZES (.tex:2749-2752, verbatim in src/aic_factorize.c). By
 * th_main_ext there is a genuine finite-dim C* algebra B and an extended
 * O(eta)-isomorphism v: B -> A = Img Phi~ (aic_cstar_build + opspace). F1 builds
 * the two maps the paper REVERSES out of the u-isomorphism corollary:
 *
 *   Delta~   = iota o v       : B -> B(H)   (.tex:2749, "v followed by A -> B(H)")
 *   Upsilon~ = v^{-1} o Phi~  : B(H) -> B   (.tex:2749, "Phi~ with target A, then v^{-1}")
 *
 * iota: A -> B(H) is the inclusion (a NO-OP on matrices: A <= M_N), so
 *   Delta~(X) = v(X) = sum_i x_i vE[i]    (literally aic_dhom_v_apply(v, X)),
 * viewing the N x N output as living in B(H) rather than A (X is the n_B x n_B
 * block-diagonal representative of an element of B). And
 *   Upsilon~(W) = v^{-1}(Phi~(W)) = sum_k <B_k, Phi~(W)>_F vinvB[k],
 * where Phi~(W) is the OBLIQUE projection of W onto A (the superop S_tilde applied,
 * NOT the Frobenius Pi_A — FINDINGS §C3), {B_k} = A's Frobenius-orthonormal basis
 * (Aec->A.B[k]), <B_k,Z>_F = sum_ij conj(B_k[i,j]) Z[i,j], and
 * vinvB[k] = v^{-1}(B_k) from aic_opspace_build_vinv.
 *
 * THE EXACT-BY-CONSTRUCTION IDENTITIES (.tex:2751, tilde_DelUps):
 *   Delta~ Upsilon~ = Phi~   (as maps M_N -> M_N: v o v^{-1} = id_A, Phi~ -> A),
 *   Upsilon~ Delta~ = 1_B    (as maps B -> B:    Phi~|_A = id, v^{-1} o v = id_B).
 * They hold to the numerical precision of the v-inversion + S_tilde idempotence;
 * tests/test_factorize.c MEASURES them (T1 eta=0 oracle vs the ORIGINAL channel
 * Phi; T2 eta>0 vs Phi~ = S_tilde, the OBLIQUE star path teeth — §C2/§C3).
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The paper's tilde maps
 * are abstract: v is the th_main_ext iso "that exists", Phi~ the cb-algebra functional
 * calculus theta(2 Phi - 1). In finite dim every piece is an explicit matrix object
 * already built and certified by a closed module: v (aic_cstar_build), v^{-1}
 * (aic_opspace_build_vinv, the SAME inverse the inverse stretch certifies), Phi~ =
 * S_tilde (aic_assoc_regularize). F1 is pure linear-algebra glue over them; the
 * paper's exact identities are the spec the tests assert.
 *
 * NON-SCOPE (later increments). F1 is the NON-UCP tilde maps only. The UCP Delta,
 * Upsilon (the 1-design CP-ization / lem_RC) are F2/F3; the certified cb-norm
 * UPPER bound ||Delta~||_cb, ||Upsilon~||_cb <= 1+O(eta) is opspace O2 (structurally
 * ||v||_cb / ||v^{-1}||_cb); the end-to-end DelUps/UpsDel2 + duals are F4.
 */
#ifndef AIC_FACTORIZE_H
#define AIC_FACTORIZE_H

#include <flint/acb_mat.h>

#include "aic_assoc.h"
#include "aic_dhom.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The F1 factorization handle. BORROWS v (B -> A, from aic_cstar_build) and the
 * aic_assoc_ecstar Aec (A = Img Phi~ + S_tilde + A's Frobenius-ON basis B[k]); OWNS
 * the vinvB array (v^{-1}(B_k), from aic_opspace_build_vinv). v and Aec must stay
 * alive for this struct's lifetime. Free with aic_factorize_clear. */
typedef struct {
    slong N;                       /* dim H = v->n = A->n                          */
    slong n_B;                     /* B's block-diagonal representative size        */
    slong dim_A;                   /* dim A = dim B                                 */
    const aic_dhom_v       *v;     /* BORROWED: B -> A (aic_cstar_build)            */
    const aic_assoc_ecstar *Aec;   /* BORROWED: A = Img Phi~ + S_tilde + basis B[k] */
    acb_mat_t *vinvB;              /* OWNED: dim_A operators v^{-1}(B_k), n_B x n_B */
} aic_factorize;

/* Build `out` from a bijective extended O(eta)-isomorphism v: B -> A and the
 * associated eps-C* algebra Aec (A = Img Phi~). Fills N, n_B, dim_A; builds vinvB
 * via aic_opspace_build_vinv. ASSERTS (fail loud, Rule 4): v->A == &Aec->A
 * (same A object — the iso's codomain IS Aec's A), v bijective (dim_A == dim_B),
 * v->n == Aec->A.n. Free with aic_factorize_clear. */
void aic_factorize_build(aic_factorize *out, const aic_dhom_v *v,
                         const aic_assoc_ecstar *Aec, slong prec);

/* Free the OWNED vinvB array and zero the struct. Does NOT touch the borrowed v
 * or Aec. */
void aic_factorize_clear(aic_factorize *out);

/* Delta~(X) = iota(v(X)) = v(X) = sum_i x_i vE[i] (.tex:2749). X is the n_B x n_B
 * block-diagonal representative of an element of B; `out` is the N x N ambient
 * image (caller-init'd N x N). A thin wrapper over aic_dhom_v_apply: the output
 * lives in B(H) rather than A (iota is a no-op on matrices). */
void aic_factorize_delta_tilde(acb_mat_t out, const aic_factorize *F,
                               const acb_mat_t X, slong prec);

/* Upsilon~(W) = v^{-1}(Phi~(W)) = sum_k <B_k, Phi~(W)>_F vinvB[k] (.tex:2749).
 * W is N x N (an operator in B(H)); `out` is the n_B x n_B block-diagonal
 * B-representative (caller-init'd n_B x n_B). Phi~(W) is the OBLIQUE projection of
 * W onto A via the superop S_tilde (aic_assoc_superop_apply over Aec->S_tilde —
 * NOT the Frobenius Pi_A, FINDINGS §C3), then the v^{-1} contraction over A's
 * Frobenius-ON basis {B_k} = Aec->A.B[k]. */
void aic_factorize_upsilon_tilde(acb_mat_t out, const aic_factorize *F,
                                 const acb_mat_t W, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_FACTORIZE_H */
