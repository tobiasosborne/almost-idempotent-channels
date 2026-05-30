/* aic_dhom_sigmin.c — the TRUE coordinate/Frobenius unit-ball lower bound of a
 * delta-homomorphism v : B -> A, the SOUND delta-inclusion COLLAPSE DETECTOR
 * (bead aic-t81, F1 fix; paper/FINDINGS.md). Realizes the delta-inclusion
 * lower-bound hypothesis (approximate_algebras.tex:451-453):
 *
 *     (1-delta)||X|| <= ||v(X)|| <= (1+delta)||X||      (X in B).
 *
 * THE BUG THIS FIXES (F1, a soundness hole). The old guards used
 * aic_dhom_prop_bounds' norm_lb = min_i ||vE[i]||_op, a BASIS sweep. That is NOT
 * the inclusion infimum inf_{X!=0} ||v(X)||/||X||: a v bounded below on every
 * basis element can still COLLAPSE a general combination (e.g. v(E_0)=diag(1,0),
 * v(E_1)=|u><u| at a small angle: each has op-norm 1, so norm_lb = 1 PASSES, but
 * ||v(E_0-E_1)|| = |sin t| ~ 0.0998 while ||E_0-E_1|| = 1, so the true inf is
 * ~0.0998 — a non-inclusion silently certified as a 0-inclusion).
 *
 * THE FIX — sigma_min of the coordinate matrix. v(X) = sum_i x_i vE[i] is LINEAR.
 * Assemble M (dim_A x dim_B) whose column i is the A-coordinates of vE[i] in A's
 * Frobenius-ORTHONORMAL basis {B_k}:
 *     M[k,i] = <B_k, vE[i]>_F = sum_{p,q} conj(B_k[p,q]) vE[i][p,q].
 * B's matrix-unit basis {E_i} is also Frobenius-orthonormal (the E_i are distinct
 * matrix units of the block-diagonal representative, <E_i,E_j>_F = delta_ij). So
 * for X with B-coords x: ||X||_F = ||x||_2 and ||v(X)||_F = ||M x||_2 (M maps
 * B-coords to A-coords, both Frobenius-orthonormal frames). Hence
 *     sigma_min(M) = min_{||x||_2 = 1} ||M x||_2 = inf_{X != 0} ||v(X)||_F/||X||_F,
 * the EXACT Frobenius/coordinate inclusion infimum. It SEES all combinations (not
 * just basis elements): sigma_min(M) = 0 iff M has a nontrivial kernel iff v
 * collapses a direction. So it is a SOUND collapse detector enforcing the
 * delta-inclusion hypothesis, where the basis sweep is blind.
 *
 * FROBENIUS vs OPERATOR (the documented subtlety). sigma_min(M) is the inf in the
 * FROBENIUS/coordinate norm, NOT the operator-norm inclusion inf the .tex states.
 * They differ by <= sqrt(dim) factors (norm equivalence on M_n), so sigma_min is
 * not the precise operator-norm a. But: it is a TRUE Frobenius unit-ball lower
 * bound AND is 0 iff v collapses a direction, so it correctly REJECTS a
 * non-inclusion (the guard's job). The exact operator-norm inclusion inf via a
 * faithful HOPM is a later cycle (cf. aic_ecstar HOPM, bead aic-0at).
 *
 * PRECISION (Law 2: reuse LAPACK; the GATE is a coarse double midpoint). M is
 * assembled in the certified arb path (exact Frobenius inner products), then its
 * midpoints feed the DOUBLE-path SVD aic_latd_singular_values (uncertified — a
 * certified sigma_min enclosure defers to aic-w4o.1/aic-w4o.2 like the other
 * extraction routines). `out` is the double sigma_min as a zero-radius arb: a
 * coarse fail-loud gate midpoint (like the projection-nontriviality gate),
 * acceptable for a 0.5 threshold. NO eigendecomposition of v itself — only the
 * standard SVD of the coordinate matrix.
 */
#include <assert.h>
#include <complex.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_latd.h"

void aic_dhom_v_sigma_min(arb_t out, const aic_dhom_v *v, slong prec)
{
    assert(v != NULL && out != NULL);
    const aic_ecstar *A = v->A;
    slong dimA = A->dim_A, dimB = v->B->dim_B, n = v->n;
    assert(dimA >= 1 && dimB >= 1);

    /* M[k,i] = <B_k, vE[i]>_F = sum_{p,q} conj(B_k[p,q]) vE[i][p,q]. */
    acb_mat_t M;
    acb_mat_init(M, dimA, dimB);
    acb_t acc, z;
    acb_init(acc);
    acb_init(z);
    for (slong k = 0; k < dimA; k++)
        for (slong i = 0; i < dimB; i++) {
            acb_zero(acc);
            for (slong p = 0; p < n; p++)
                for (slong q = 0; q < n; q++) {
                    acb_conj(z, acb_mat_entry(A->B[k], p, q));
                    acb_mul(z, z, acb_mat_entry(v->vE[i], p, q), prec);
                    acb_add(acc, acc, z, prec);
                }
            acb_set(acb_mat_entry(M, k, i), acc);
        }
    acb_clear(z);
    acb_clear(acc);

    /* sigma_min(M) = smallest singular value of the dimA x dimB coordinate matrix.
     * If dimB > dimA, M cannot have full column rank, so sigma_min = 0 (v MUST
     * collapse a direction — the guard rejects). aic_latd_singular_values returns
     * min(m,n) values DESCENDING; the inclusion inf over the dimB-dim coordinate
     * space is the dimB-th singular value, which is 0 when dimB > dimA. */
    double _Complex *Md = flint_malloc((size_t) dimA * (size_t) dimB
                                       * sizeof(double _Complex));
    aic_latd_from_acb_mat(Md, M);               /* midpoints, row-major */
    slong msv = dimA < dimB ? dimA : dimB;
    double *sv = flint_malloc((size_t) msv * sizeof(double));
    aic_latd_singular_values(sv, Md, dimA, dimB);

    double sigma_min;
    if (dimB > dimA)
        sigma_min = 0.0;                        /* rank-deficient: a kernel exists */
    else
        sigma_min = sv[dimB - 1];               /* smallest of the dimB sing. vals */

    arb_set_d(out, sigma_min);                  /* zero-radius gate midpoint */

    flint_free(sv);
    flint_free(Md);
    acb_mat_clear(M);
}
