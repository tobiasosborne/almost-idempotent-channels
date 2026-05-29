/* aic_ecstar.c — eps-C* algebra DATA MODEL: lifecycle, construction from an
 * exactly-idempotent decomposition, the Choi-Effros star, and the
 * subspace-projection residual (bead aic-knm, module "ecstar").
 *
 * Realizes the data model of approximate_algebras.tex:407-439 (the eps-C* axiom
 * block) with the multiplication being the Choi-Effros star
 *     X * Y = Phi(X . Y)                       (.tex:341-342, :2187-2189)
 * the involution X |-> X^dag, the unit I = 1_n (.tex:2186-2187), and the operator
 * norm inherited from B(H) (.tex:2192). See include/aic_ecstar.h for the full
 * data-model contract, the borrowed-phi lifetime rule, and the vec convention.
 *
 * The defect ESTIMATORS live in aic_ecstar_defect.c; this file is the substrate
 * they sweep over. Convention-safe by construction: every map is applied to an
 * actual n x n matrix and compared via aic_mat_opnorm; the n^2 x n^2
 * superoperator is NEVER formed (CLAUDE.md vectorization-convention callout).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_ecstar.h"
#include "aic_mat.h"
#include "aic_ucp.h"

void aic_ecstar_init(aic_ecstar *out, slong n, slong dim_A,
                     const aic_ucp_kraus *phi)
{
    assert(out != NULL);
    assert(phi != NULL);
    /* fail loud (Rule 4): the star's Phi must be a self-map on M_n */
    assert(phi->dim_K == n && phi->dim_H == n &&
           "aic_ecstar_init: phi must be a self-map with dim_K==dim_H==n");
    assert(dim_A >= 0);

    out->n = n;
    out->dim_A = dim_A;
    out->phi = phi;
    out->B = (acb_mat_t *) flint_malloc((size_t) dim_A * sizeof(acb_mat_t));
    assert(out->B != NULL || dim_A == 0);
    for (slong k = 0; k < dim_A; k++) {
        acb_mat_init(out->B[k], n, n);
        acb_mat_zero(out->B[k]);
    }
}

void aic_ecstar_clear(aic_ecstar *out)
{
    if (out == NULL) return;
    for (slong k = 0; k < out->dim_A; k++) acb_mat_clear(out->B[k]);
    if (out->B != NULL) flint_free(out->B);
    out->B = NULL;
    out->dim_A = 0;
    /* phi is BORROWED: never freed here. */
}

void aic_ecstar_from_idemp(aic_ecstar *out, const aic_idemp_decomp *d,
                           const aic_ucp_kraus *phi, slong prec)
{
    assert(out != NULL && d != NULL && phi != NULL);
    slong n = out->n;
    /* fail loud (Rule 4): dims must agree across the decomposition and phi */
    assert(d->n == n && "aic_ecstar_from_idemp: d->n != out->n");
    assert(phi->dim_H == n && phi->dim_K == n &&
           "aic_ecstar_from_idemp: phi not a self-map on M_n");
    assert(out->dim_A == d->dim_A &&
           "aic_ecstar_from_idemp: out->dim_A != d->dim_A");
    assert(acb_mat_nrows(d->Delta) == n * n &&
           acb_mat_ncols(d->Delta) == d->dim_A &&
           "aic_ecstar_from_idemp: Delta shape mismatch");

    /* B_k[i,j] = Delta[i*n + j, k] (reshape column k ROW-MAJOR; vec convention). */
    for (slong k = 0; k < d->dim_A; k++)
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
                acb_set(acb_mat_entry(out->B[k], i, j),
                        acb_mat_entry(d->Delta, i * n + j, k));

    /* ASSERT Frobenius-orthonormality of the columns: ||<B_j,B_k>_F - delta_jk||
     * must be below a prec-appropriate floor. <B_j,B_k>_F = sum_pq conj(B_j[p,q])
     * B_k[p,q]. This is the inherited Delta^dag Delta = 1_A relation
     * (include/aic_idemp.h); reshaping preserves it, so a failure means the input
     * decomposition is corrupt. Fail loud (Rule 4). */
    {
        arb_t worst, floor_tol, mag;
        acb_t ip, t;
        arb_init(worst);
        arb_init(floor_tol);
        arb_init(mag);
        acb_init(ip);
        acb_init(t);
        arb_zero(worst);
        for (slong j = 0; j < d->dim_A; j++)
            for (slong k = 0; k < d->dim_A; k++) {
                acb_zero(ip);
                for (slong p = 0; p < n; p++)
                    for (slong q = 0; q < n; q++) {
                        acb_conj(t, acb_mat_entry(out->B[j], p, q));
                        acb_mul(t, t, acb_mat_entry(out->B[k], p, q), prec);
                        acb_add(ip, ip, t, prec);
                    }
                if (j == k) acb_sub_si(ip, ip, 1, prec);
                acb_abs(mag, ip, prec);
                if (arb_gt(mag, worst)) arb_set(worst, mag);
            }
        /* floor: 1e-9 absolute, generous for a double-path-extracted Delta whose
         * orthonormality floors near machine eps but can be ~1e-12 (see
         * aic_idemp wG note). A corrupt basis is O(1). */
        arb_set_d(floor_tol, 1e-9);
        assert(arb_le(worst, floor_tol) &&
               "aic_ecstar_from_idemp: basis columns not Frobenius-orthonormal");
        acb_clear(t);
        acb_clear(ip);
        arb_clear(mag);
        arb_clear(floor_tol);
        arb_clear(worst);
    }
}

void aic_ecstar_star(acb_mat_t out, const aic_ecstar *A, const acb_mat_t X,
                     const acb_mat_t Y, slong prec)
{
    assert(A != NULL && A->phi != NULL);
    slong n = A->n;
    assert(acb_mat_nrows(X) == n && acb_mat_ncols(X) == n);
    assert(acb_mat_nrows(Y) == n && acb_mat_ncols(Y) == n);
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n);

    acb_mat_t XY;
    acb_mat_init(XY, n, n);
    acb_mat_mul(XY, X, Y, prec);          /* ordinary product X . Y         */
    aic_ucp_apply(out, A->phi, XY, prec); /* Phi(XY) = X * Y (.tex:341-342) */
    acb_mat_clear(XY);
}

void aic_ecstar_proj_residual(arb_t out, const aic_ecstar *A, const acb_mat_t X,
                              slong prec)
{
    assert(A != NULL);
    slong n = A->n, d = A->dim_A;
    assert(acb_mat_nrows(X) == n && acb_mat_ncols(X) == n);

    /* R = X - sum_k <B_k,X>_F B_k, with <B_k,X>_F = sum_pq conj(B_k[p,q]) X[p,q]. */
    acb_mat_t R;
    acb_mat_init(R, n, n);
    acb_mat_set(R, X);
    acb_t coeff, t;
    acb_init(coeff);
    acb_init(t);
    for (slong k = 0; k < d; k++) {
        acb_zero(coeff);
        for (slong p = 0; p < n; p++)
            for (slong q = 0; q < n; q++) {
                acb_conj(t, acb_mat_entry(A->B[k], p, q));
                acb_mul(t, t, acb_mat_entry(X, p, q), prec);
                acb_add(coeff, coeff, t, prec);
            }
        for (slong p = 0; p < n; p++)
            for (slong q = 0; q < n; q++) {
                acb_mul(t, coeff, acb_mat_entry(A->B[k], p, q), prec);
                acb_sub(acb_mat_entry(R, p, q), acb_mat_entry(R, p, q), t, prec);
            }
    }
    aic_mat_opnorm(out, R, prec);
    acb_clear(t);
    acb_clear(coeff);
    acb_mat_clear(R);
}
