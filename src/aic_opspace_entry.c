/* aic_opspace_entry.c — the PUBLIC entry points of the opspace operator-norm
 * ampliation (bead aic-zwo, increment O1): the forward/inverse stretch HOPM
 * drivers, the STRUCTURAL ampliation entry (the FINDINGS §C12 tooth's hook), and
 * the v^{-1} builder + inverse Smith level the factorize/D4 §12 decode reads
 * (docs/research/factorize_d4_research.md §4 adds #1/#3). Split from the opmap
 * builders (aic_opspace_map.c) and the HOPM kernel (aic_opspace_ampliate.c) to
 * keep each file under the ~200 LOC core limit (CLAUDE.md Rule 10). See
 * include/aic_opspace.h for the contract; the shared builder interface is
 * src/aic_opspace_map_internal.h.
 */
#include <assert.h>
#include <complex.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_opspace.h"
#include "aic_opspace_map_internal.h"

void aic_opspace_forward_stretch(arb_t out, const aic_dhom_v *v, slong n,
                                 slong prec)
{
    assert(v != NULL && out != NULL && n >= 1);
    assert(v->A->dim_A == v->B->dim_B && "forward_stretch: v not bijective");
    opmap m;
    aic_opspace_opmap_forward(&m, v, prec);
    double s = aic_opmap_hopm_stretch(&m, n, OPSPACE_NREST, OPSPACE_NITER);
    arb_set_d(out, s);
    aic_opmap_clear(&m);
}

void aic_opspace_inverse_stretch(arb_t out, const aic_dhom_v *v, slong n,
                                 slong prec)
{
    assert(v != NULL && out != NULL && n >= 1);
    assert(v->A->dim_A == v->B->dim_B && "inverse_stretch: v not bijective");
    opmap m;
    aic_opspace_opmap_inverse(&m, v, prec);
    double s = aic_opmap_hopm_stretch(&m, n, OPSPACE_NREST, OPSPACE_NITER);
    arb_set_d(out, s);
    aic_opmap_clear(&m);
}

/* Structural ampliation Y = (1_{M_n} (x) v)(X). Builds the forward opmap and
 * drives the SAME aic_opmap_apply_fn the HOPM kernel uses (FINDINGS §C12 tooth:
 * the structural-ampliation cross-check in test_opspace.c pins this against an
 * independent Kronecker block assembly so a (0,0)-block-only cripple is caught). */
void aic_opspace_ampliate_forward(double _Complex *Y, const aic_dhom_v *v,
                                  const double _Complex *X, slong n, slong prec)
{
    assert(v != NULL && Y != NULL && X != NULL && n >= 1);
    assert(v->A->dim_A == v->B->dim_B && "ampliate_forward: v not bijective");
    opmap m;
    aic_opspace_opmap_forward(&m, v, prec);
    cx *xb = aic_opmap_cxalloc(m.sU * m.sU), *yb = aic_opmap_cxalloc(m.sV * m.sV);
    aic_opmap_apply_fn(Y, &m, X, n, xb, yb);
    flint_free(yb);
    flint_free(xb);
    aic_opmap_clear(&m);
}

/* v^{-1}(B_k) = sum_i (M_1^{-1})[i,k] E_i, each n_B x n_B block-diagonal in B.
 * The SAME M_1^{-1} the inverse stretch certifies (aic_opspace_build_M1inv), so
 * factorize's Upsilon~ = v^{-1} o Phi~ uses the certified-and-used inverse. */
void aic_opspace_build_vinv(acb_mat_t **vinvB, slong *dim_A, slong *n_B,
                            const aic_dhom_v *v, slong prec)
{
    assert(vinvB != NULL && v != NULL);
    slong dimA = v->A->dim_A, dimB = v->B->dim_B, nB = v->B->n_B;
    assert(dimA == dimB && "build_vinv: v not bijective (dim_A != dim_B)");
    slong d = dimA;
    cx *Minv = aic_opmap_cxalloc(d * d);
    aic_opspace_build_M1inv(Minv, v, d, prec); /* Minv[i*d + k] = (M_1^{-1})_{ik} */

    acb_mat_t *out = flint_malloc((size_t) d * sizeof(acb_mat_t));
    acb_mat_t Ei, term;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(term, nB, nB);
    arb_t re, im;
    acb_t coef;
    arb_init(re);
    arb_init(im);
    acb_init(coef);
    for (slong k = 0; k < d; k++) {
        acb_mat_init(out[k], nB, nB);
        acb_mat_zero(out[k]);
        for (slong i = 0; i < dimB; i++) {
            cx c = Minv[i * d + k];
            aic_dhom_B_matunit(Ei, v->B, i);
            arb_set_d(re, creal(c));
            arb_set_d(im, cimag(c));
            acb_set_arb_arb(coef, re, im);
            acb_mat_scalar_mul_acb(term, Ei, coef, prec);
            acb_mat_add(out[k], out[k], term, prec);
        }
    }
    acb_clear(coef);
    arb_clear(im);
    arb_clear(re);
    acb_mat_clear(term);
    acb_mat_clear(Ei);
    flint_free(Minv);

    *vinvB = out;
    if (dim_A) *dim_A = d;
    if (n_B) *n_B = nB;
}

void aic_opspace_vinv_clear(acb_mat_t *vinvB, slong dim_A)
{
    if (vinvB == NULL) return;
    for (slong k = 0; k < dim_A; k++) acb_mat_clear(vinvB[k]);
    flint_free(vinvB);
}
