/* aic_ecstar_assoc.c — the associativity defect estimator (bead aic-knm,
 * module "ecstar"). Split from aic_ecstar_defect.c to keep each file <= 200 LOC
 * (Rule 10); this is the only TRILINEAR sweep (d^3 star products) and the
 * headline cost of the basis-sweep rung.
 *
 * approximate_algebras.tex:412-413 ax_assoc:
 *   ||(XY)Z - X(YZ)|| <= eps ||X|| ||Y|| ||Z||   (X,Y,Z in A).
 * The associator h(X,Y,Z) = (X*Y)*Z - X*(Y*Z) is TRILINEAR in (X,Y,Z) (each star
 * is bilinear and Phi is linear), so it vanishes on every basis triple iff it
 * vanishes identically: sweeping the basis triples is an EXACT zero-detector for
 * ax_assoc (no operator-norm-unit-ball search is needed to certify the ZERO of
 * the eta=0 oracle). It is NOT the faithful sup-over-unit-ball eps (that is the
 * later HOPM cycle, bead aic-0at); but it is the correct quantity to certify the
 * Choi-Effros product associates exactly when eta=0.
 *
 * Convention-safe: products are the Choi-Effros star aic_ecstar_star (Phi(XY)),
 * compared via aic_mat_opnorm. The n^2 x n^2 superoperator is never formed.
 */
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_ecstar.h"
#include "aic_ecstar_internal.h"

/* max over basis triples (i,j,k) of ||(B_i*B_j)*B_k - B_i*(B_j*B_k)||_op. */
void aic_ecstar_defect_assoc(arb_t out, const aic_ecstar *A, slong prec)
{
    slong n = A->n, d = A->dim_A;
    acb_mat_t ij, jk, lhs, rhs, diff;
    acb_mat_init(ij, n, n);
    acb_mat_init(jk, n, n);
    acb_mat_init(lhs, n, n);
    acb_mat_init(rhs, n, n);
    acb_mat_init(diff, n, n);
    arb_zero(out);
    for (slong i = 0; i < d; i++)
        for (slong j = 0; j < d; j++) {
            aic_ecstar_star(ij, A, A->B[i], A->B[j], prec); /* B_i * B_j */
            for (slong k = 0; k < d; k++) {
                aic_ecstar_star(jk, A, A->B[j], A->B[k], prec);  /* B_j * B_k    */
                aic_ecstar_star(lhs, A, ij, A->B[k], prec);      /* (B_i*B_j)*B_k */
                aic_ecstar_star(rhs, A, A->B[i], jk, prec);      /* B_i*(B_j*B_k) */
                acb_mat_sub(diff, lhs, rhs, prec);
                aic_ecstar_bump_opnorm(out, diff, prec);
            }
        }
    acb_mat_clear(diff);
    acb_mat_clear(rhs);
    acb_mat_clear(lhs);
    acb_mat_clear(jk);
    acb_mat_clear(ij);
}
