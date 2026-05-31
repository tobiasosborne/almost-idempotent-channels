/* aic_cstar_merge.c — the two SIMPLER assembly lemmas of the th_main master loop
 * (bead aic-097, module "cstar_build", Increment 2). See include/aic_cstar.h
 * (Increment 2 banner) and docs/research/cstar_build_design.md §4.1 (lem_add_dim),
 * §4.2 (cor_merge_sum), §3 (the aic_dhom_v direct-construction claim).
 *
 * approximate_algebras.tex:1363-1369 — Lemma lem_add_dim:
 *   "Let B and C be commutative C* algebras with projection bases {Pi_1,...,Pi_p}
 *    and {Sigma_1,...,Sigma_q}. ... Let P_j = v(Pi_j), Q_k = w(Sigma_k), and also
 *    P = v(I_B), Q = w(I_C). Then there is a bounded linear bijection between
 *    S_{P,Q} and (+)_{j,k} S_{P_j,Q_k}; in particular,
 *        dim S_{P,Q} = sum_{j,k} dim S_{P_j,Q_k}."
 *
 * approximate_algebras.tex:1352-1358 — Corollary cor_merge_sum (Merging):
 *   "Let B_1 and B_2 be C* algebras, and let P_1, P_2 be delta-projections in an
 *    eps-C* algebra A such that ||P_1+P_2-I|| <= delta. Consider some
 *    delta-inclusions v_j : B_j -> S_{P_j} (for j=1,2). Then the combined map
 *        v : (X_1, X_2) |-> v_1(X_1) + v_2(X_2) : B_1 (+) B_2 -> A
 *    is an O(delta+eps)-inclusion. If v_1, v_2 are bijective and S_{P_1,P_2} = 0,
 *    then v is also bijective."
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3).
 *   lem_add_dim's proof composes lem_alpha bijections; the conclusion we need is
 *   the DIMENSION identity, which aic_corner_alpha_dims already certifies (it
 *   builds the per-block Co's + Co_{P,Q} and sums the extracted dims). So
 *   aic_cstar_lem_add_dim is the dimension query + the fail-loud N == d_PQ assert
 *   that IS the lemma's conclusion. The Stage-3 corollary S_{P_C,P_D}=0 is a thin
 *   aic_corner_dim_S wrapper.
 *
 *   cor_merge_sum's "combined map" is, in the aic_dhom_v model, a pure
 *   CONCATENATION of the two basis-image lists. B's matrix-unit linear index
 *   mu(l,a,b) = offset_l + a*d_l + b is CUMULATIVE across blocks (aic_dhom.h), so
 *   for the merged B = B_1 (+) B_2 the matrix unit E_i of B_1 keeps its index and
 *   the matrix unit E_j of B_2 lands at dim_B1 + j. Hence vE = [v1->vE, v2->vE]
 *   IS v(E_i) for the merged B, i.e. v(X_1,X_2) = v_1(X_1) + v_2(X_2) (the
 *   block-diagonal X has X_1 in block-range [0,dim_B1) and X_2 in the rest, and v
 *   is linear in B-coords). No new map computation; we copy and certify.
 *
 * THE STAR IS LOAD-BEARING (CLAUDE.md domain callout, FINDINGS §C2). The merged
 *   v's multiplicativity defect (aic_dhom_defect_sweep) is taken w.r.t. A's
 *   Choi-Effros STAR X*Y = Phi_tilde(XY), inherited via v_out->A = A. The eta=0
 *   identity channel has Phi_tilde = id so its star == plain product and is BLIND
 *   to a star bug; the oblique eta>0 fixture (make_mixconj) supplies the teeth
 *   (tests/test_cstar_merge.c B2). errreduce (cor_improvement) is applied by the
 *   CALLER after this routine (.tex:1426/1443), not here.
 */
#include <assert.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_corner.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"

/* lem_add_dim (.tex:1363): certify dim S_{P,Q} = sum_{jk} dim S_{P_j,Q_k} and
 * return the common value. aic_corner_alpha_dims computes BOTH N = sum dim
 * S_{P_j,Q_k} and d_PQ = dim S_{P,Q}; the fail-loud N == d_PQ assert IS the
 * lemma's conclusion (Rule 4: a wrong P_part feeding a mismatched sum aborts). */
slong aic_cstar_lem_add_dim(const aic_ecstar *A,
                            const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                            const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                            slong prec)
{
    assert(A != NULL && P_parts != NULL && Q_parts != NULL);
    assert(p >= 1 && q >= 1 && "aic_cstar_lem_add_dim: need p,q >= 1");

    slong N = -1, d_PQ = -1;
    aic_corner_alpha_dims(&N, &d_PQ, A, P_parts, p, P, Q_parts, q, Q, prec);

    /* The lem_add_dim conclusion (.tex:1364): the block sum equals the total. */
    assert(N == d_PQ &&
           "aic_cstar_lem_add_dim: sum_jk dim S_{Pj,Qk} != dim S_{P,Q} "
           "(lem_add_dim conclusion violated)");

    return d_PQ;
}

/* Stage-3 off-diagonal certification (.tex:1443): dim S_{P_C,P_D} == 0 for
 * distinct equivalence classes. Thin wrapper over aic_corner_dim_S. */
int aic_cstar_off_diag_zero(const aic_ecstar *A, const acb_mat_t P_C,
                            const acb_mat_t P_D, slong prec)
{
    assert(A != NULL);
    return aic_corner_dim_S(A, P_C, P_D, prec) == 0;
}

/* cor_merge_sum (.tex:1352): build B_out = B_1 (+) B_2 and v_out with
 * vE = concat(v1->vE, v2->vE); certify the merged multiplicativity defect and the
 * Frobenius inclusion lower bound. v_out->A = A, v_out->B = B_out (the caller
 * keeps B_out alive for v_out's lifetime). vE operators are DEEP COPIES. */
void aic_cstar_merge_sum(aic_dhom_B *B_out, aic_dhom_v *v_out,
                         arb_t mult_def, arb_t sigma_min,
                         const aic_dhom_v *v1, const aic_dhom_v *v2,
                         const aic_ecstar *A, slong prec)
{
    assert(B_out != NULL && v_out != NULL && v1 != NULL && v2 != NULL && A != NULL);
    /* Both inclusions must target the SAME containing algebra A (.tex:1355). */
    assert(v1->A == A && "aic_cstar_merge_sum: v1->A != A");
    assert(v2->A == A && "aic_cstar_merge_sum: v2->A != A");
    assert(v1->n == A->n && v2->n == A->n &&
           "aic_cstar_merge_sum: ambient dim mismatch (v->n != A->n)");

    const aic_dhom_B *B1 = v1->B, *B2 = v2->B;
    slong m1 = B1->num_blocks, m2 = B2->num_blocks;

    /* Merged block-dim list = concat(B_1.d, B_2.d). aic_dhom_B_init copies it. */
    slong m = m1 + m2;
    slong *dims = (slong *) flint_malloc((size_t) m * sizeof(slong));
    assert(dims != NULL);
    for (slong l = 0; l < m1; l++) dims[l] = B1->d[l];
    for (slong l = 0; l < m2; l++) dims[m1 + l] = B2->d[l];
    aic_dhom_B_init(B_out, dims, m);
    flint_free(dims);

    /* The matrix-unit index is cumulative, so dim_B(B_out) == dim_B1 + dim_B2 and
     * the concatenation vE = [v1->vE, v2->vE] is exactly v(E_i) for the merged B
     * (design §3): block l < m1 -> v_1's images; block l >= m1 -> v_2's. */
    assert(B_out->dim_B == B1->dim_B + B2->dim_B &&
           "aic_cstar_merge_sum: merged dim_B != dim_B1 + dim_B2 (index model)");

    aic_dhom_v_init(v_out, B_out, A);
    for (slong i = 0; i < B1->dim_B; i++)
        acb_mat_set(v_out->vE[i], v1->vE[i]);          /* deep copy, not alias */
    for (slong j = 0; j < B2->dim_B; j++)
        acb_mat_set(v_out->vE[B1->dim_B + j], v2->vE[j]);

    /* Certified multiplicativity defect (A's STAR, FINDINGS §C2) + the SOUND
     * Frobenius inclusion lower bound (FINDINGS §C6). Both NULL-skippable. */
    if (mult_def != NULL) aic_dhom_defect_sweep(mult_def, v_out, prec);
    if (sigma_min != NULL) aic_dhom_v_sigma_min(sigma_min, v_out, prec);
}
