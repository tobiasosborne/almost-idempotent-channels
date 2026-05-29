/* aic_ecstar_involution.c — the involution-residual estimator (ax_*) for the
 * eps-C* algebra data model (bead aic-knm, module "ecstar"). Split from
 * aic_ecstar_defect.c (Rule 10) and refactored around a shared GENERIC-APPLY
 * core so the production code path is exercised by a non-HP teeth.
 *
 * approximate_algebras.tex:422-423 ax_*:
 *   ||X^dag|| = ||X||,  (XY)^dag = Y^dag X^dag   (X,Y in A), EXACT equalities.
 * The estimated quantity is
 *   max over pairs (j,k) of ||(B_j*B_k)^dag - B_k^dag * B_j^dag||_op,   X*Y=Phi(XY).
 * Expanding the star:
 *   (B_j*B_k)^dag - B_k^dag*B_j^dag = Phi(B_j B_k)^dag - Phi((B_j B_k)^dag),
 * which is IDENTICALLY 0 for any Hermicity-preserving (HP) Phi (.tex:2211
 * "(X*Y)^dag = Y^dag*X^dag"). Since the aic_ucp_kraus model Phi(X)=sum K^dag X K
 * is HP for EVERY Kraus set, the production estimator is a STRUCTURAL invariant
 * (a consistency check that the star's Phi is HP), eta-independent — NOT a free
 * measurement. (||X^dag||=||X|| is exact by unitary invariance of the operator
 * norm; not swept.)
 *
 * THE DEFECT-CLASS GUARD (Rule 7). Because a correct impl and `return 0` are
 * observationally identical on all aic_ecstar inputs, the residual computation
 * is factored into aic_ecstar_involution_core, parameterised by the linear map
 * (aic_ecstar_apply_fn) that plays the role of Phi. The production estimator
 * passes a thunk wrapping aic_ucp_apply(A->phi, .) (so its behaviour is
 * UNCHANGED — still ~0 on HP maps), and tests/test_ecstar.c teeth_non_hp passes
 * a synthetic NON-HP map T(X)=M X N with M != N^dag, for which the core residual
 * is O(non-HP-ness) and a return-0 / wrong-matrix mutation of the core turns the
 * teeth RED. The core is therefore genuinely protected.
 *
 * Convention-safe: products are ordinary matrix products fed to the generic
 * apply; the n^2 x n^2 superoperator is never formed.
 */
#include <assert.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_ecstar.h"
#include "aic_ecstar_internal.h"
#include "aic_ucp.h"

/* The shared metric, generalised in the map `apply` (= the star's Phi). For each
 * ordered pair (j,k): apply(B_j B_k)^dag - apply(B_k^dag B_j^dag), op-norm, max.
 * HP apply => identically 0; non-HP apply => generally nonzero (the teeth). */
void aic_ecstar_involution_core(arb_t out, const acb_mat_t *B, slong d, slong n,
                                aic_ecstar_apply_fn apply, void *ctx, slong prec)
{
    assert(apply != NULL);
    acb_mat_t prod, lhs, Bjd, Bkd, prodd, rhs, diff;
    acb_mat_init(prod, n, n);   /* B_j B_k (ordinary product)      */
    acb_mat_init(lhs, n, n);    /* apply(B_j B_k)^dag              */
    acb_mat_init(Bjd, n, n);    /* B_j^dag                        */
    acb_mat_init(Bkd, n, n);    /* B_k^dag                        */
    acb_mat_init(prodd, n, n);  /* B_k^dag B_j^dag                */
    acb_mat_init(rhs, n, n);    /* apply(B_k^dag B_j^dag)          */
    acb_mat_init(diff, n, n);
    arb_zero(out);
    for (slong j = 0; j < d; j++)
        for (slong k = 0; k < d; k++) {
            acb_mat_mul(prod, B[j], B[k], prec);       /* B_j B_k             */
            apply(lhs, prod, ctx, prec);               /* apply(B_j B_k)      */
            acb_mat_conjugate_transpose(lhs, lhs);     /* apply(B_j B_k)^dag  */
            acb_mat_conjugate_transpose(Bjd, B[j]);    /* B_j^dag             */
            acb_mat_conjugate_transpose(Bkd, B[k]);    /* B_k^dag             */
            acb_mat_mul(prodd, Bkd, Bjd, prec);        /* B_k^dag B_j^dag     */
            apply(rhs, prodd, ctx, prec);              /* apply(B_k^dag B_j^dag) */
            acb_mat_sub(diff, lhs, rhs, prec);
            aic_ecstar_bump_opnorm(out, diff, prec);
        }
    acb_mat_clear(diff);
    acb_mat_clear(rhs);
    acb_mat_clear(prodd);
    acb_mat_clear(Bkd);
    acb_mat_clear(Bjd);
    acb_mat_clear(lhs);
    acb_mat_clear(prod);
}

/* Production thunk: apply the star's Phi via aic_ucp_apply (HP for any Kraus
 * set). ctx is the borrowed aic_ucp_kraus*. */
static void apply_phi(acb_mat_t out, const acb_mat_t X, void *ctx, slong prec)
{
    aic_ucp_apply(out, (const aic_ucp_kraus *) ctx, X, prec);
}

void aic_ecstar_defect_involution(arb_t out, const aic_ecstar *A, slong prec)
{
    assert(A != NULL && A->phi != NULL);
    aic_ecstar_involution_core(out, (const acb_mat_t *) A->B, A->dim_A, A->n,
                               apply_phi, (void *) A->phi, prec);
}
