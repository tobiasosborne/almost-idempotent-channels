/* aic_idemp_verify.c — certified arb defects for the five th_idemp_structure
 * relations (bead aic-wuh; the eta=0 oracle, approximate_algebras.tex:2080-2090).
 * Each defect is 0 (to rounding) for an exactly idempotent Phi.
 *
 * Convention-safe by construction: maps are applied to matrix units / basis
 * coords and compared via aic_mat_opnorm; the n^2 x n^2 superoperator is NEVER
 * formed (CLAUDE.md domain callout: vectorization-convention bugs). All arithmetic
 * is on the certified acb/arb path; the double-path subspaces (J_M, Delta) are
 * lifted as the fixed inputs and the relations are CERTIFIED on top of them.
 *
 * Coordinate maps (see include/aic_idemp.h vec convention). With Delta's columns
 * an ONB of A:
 *   - Delta(coords) = sum_k coords_k B_k : C^{dim_A} -> Img Phi (n x n);
 *   - the A-coordinates of an element Z in A are Delta^dag vec(Z) (orthonormal cols);
 *   - Gamma : B(M) -> A returns A-COORDINATES: Gamma(Y)_l = sum_ij Gamma[l,i*dM+j]Y[i,j];
 *   - w : A -> B(M) takes A-coordinates: w(coords) = sum_k coords_k (reshape col k of w).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_idemp.h"
#include "aic_mat.h"
#include "aic_ucp.h"

/* y (length dA, as a dA x 1 acb_mat) = Gamma applied to Y (dM x dM): the
 * row-major-vec contraction Gamma[l, i*dM+j] Y[i,j]. */
static void gamma_apply(acb_mat_t y, const acb_mat_t Gamma, const acb_mat_t Y,
                        slong dA, slong dM, slong prec)
{
    acb_t acc, term;
    acb_init(acc);
    acb_init(term);
    for (slong l = 0; l < dA; l++) {
        acb_zero(acc);
        for (slong i = 0; i < dM; i++)
            for (slong j = 0; j < dM; j++) {
                acb_mul(term, acb_mat_entry(Gamma, l, i * dM + j),
                        acb_mat_entry(Y, i, j), prec);
                acb_add(acc, acc, term, prec);
            }
        acb_set(acb_mat_entry(y, l, 0), acc);
    }
    acb_clear(term);
    acb_clear(acc);
}

/* C_M(X) = J_M^dag X J_M : out (dM x dM) from X (n x n). */
static void compress(acb_mat_t out, const acb_mat_t J_M, const acb_mat_t X,
                     slong n, slong dM, slong prec)
{
    acb_mat_t JMd, JdX;
    acb_mat_init(JMd, dM, n);
    acb_mat_init(JdX, dM, n);
    acb_mat_conjugate_transpose(JMd, J_M);
    acb_mat_mul(JdX, JMd, X, prec);
    acb_mat_mul(out, JdX, J_M, prec);
    acb_mat_clear(JdX);
    acb_mat_clear(JMd);
}

/* Delta(coords): Z (n x n) = sum_k coords_k B_k = reshape of Delta @ coords. */
static void delta_apply(acb_mat_t Z, const acb_mat_t Delta, const acb_mat_t coords,
                        slong n, slong prec)
{
    acb_mat_t v;          /* (n*n) x 1 */
    acb_mat_init(v, n * n, 1);
    acb_mat_mul(v, Delta, coords, prec);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(Z, a, b), acb_mat_entry(v, a * n + b, 0));
    acb_mat_clear(v);
}

/* reshape column k of Delta into B (n x n). */
static void col_to_mat(acb_mat_t B, const acb_mat_t Delta, slong k, slong n)
{
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(B, a, b), acb_mat_entry(Delta, a * n + b, k));
}

/* ||Gamma C_M Delta - 1_A||_op (.tex:2081). */
void aic_idemp_defect_GCD(arb_t out, const aic_idemp_decomp *d, slong prec)
{
    slong n = d->n, dM = d->dim_M, dA = d->dim_A;
    acb_mat_t B, cB, gB, GCD, diff;
    acb_mat_init(B, n, n);
    acb_mat_init(cB, dM, dM);
    acb_mat_init(gB, dA, 1);
    acb_mat_init(GCD, dA, dA);
    acb_mat_init(diff, dA, dA);
    for (slong k = 0; k < dA; k++) {           /* e_k -> Delta -> C_M -> Gamma */
        col_to_mat(B, d->Delta, k, n);          /* Delta(e_k) = B_k */
        compress(cB, d->J_M, B, n, dM, prec);   /* C_M(B_k) */
        gamma_apply(gB, d->Gamma, cB, dA, dM, prec);
        for (slong l = 0; l < dA; l++)
            acb_set(acb_mat_entry(GCD, l, k), acb_mat_entry(gB, l, 0));
    }
    acb_mat_one(diff);
    acb_mat_sub(diff, GCD, diff, prec);
    aic_mat_opnorm(out, diff, prec);
    acb_mat_clear(diff);
    acb_mat_clear(GCD);
    acb_mat_clear(gB);
    acb_mat_clear(cB);
    acb_mat_clear(B);
}

/* max_ij ||Delta Gamma C_M (E_ij) - Phi(E_ij)||_op (.tex:2072). */
void aic_idemp_defect_DGC(arb_t out, const aic_idemp_decomp *d,
                          const aic_ucp_kraus *phi, slong prec)
{
    slong n = d->n, dM = d->dim_M, dA = d->dim_A;
    acb_mat_t E, cE, gE, Z, PhiE, diff;
    arb_t dd;
    acb_mat_init(E, n, n);
    acb_mat_init(cE, dM, dM);
    acb_mat_init(gE, dA, 1);
    acb_mat_init(Z, n, n);
    acb_mat_init(PhiE, n, n);
    acb_mat_init(diff, n, n);
    arb_init(dd);
    arb_zero(out);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            acb_mat_zero(E);
            acb_set_si(acb_mat_entry(E, i, j), 1);
            compress(cE, d->J_M, E, n, dM, prec);        /* C_M(E_ij) */
            gamma_apply(gE, d->Gamma, cE, dA, dM, prec); /* Gamma(...) coords */
            delta_apply(Z, d->Delta, gE, n, prec);       /* Delta(coords) */
            aic_ucp_apply(PhiE, phi, E, prec);           /* Phi(E_ij) */
            acb_mat_sub(diff, Z, PhiE, prec);
            aic_mat_opnorm(dd, diff, prec);
            if (arb_gt(dd, out)) arb_set(out, dd);
        }
    arb_clear(dd);
    acb_mat_clear(diff);
    acb_mat_clear(PhiE);
    acb_mat_clear(Z);
    acb_mat_clear(gE);
    acb_mat_clear(cE);
    acb_mat_clear(E);
}
