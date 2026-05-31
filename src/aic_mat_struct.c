/* aic_mat_struct.c — tensor structure over acb_mat_t (beads aic-9zh, "T-mat"):
 * Kronecker product and partial traces. See include/aic_mat.h for the API and
 * the LOAD-BEARING index convention (left factor major, row-major Kronecker).
 *
 * These are the substrate for the UCP module's decode map
 * rho_j = Tr_{E_j}(W_j rho W_j^dag) (approximate_algebras.tex:288) over the
 * tensor factoring M -> L_j (x) E_j (approximate_algebras.tex:263-264). FLINT
 * 3.0.1 has neither a Kronecker product nor a partial trace, so they are built
 * here directly over acb_mat_entry / acb_trace-style accumulation, reusing
 * acb_mul / acb_add rather than reinventing arithmetic (CLAUDE.md Law 2).
 *
 * The convention is fixed so that the two contraction identities hold exactly
 * (up to rounding) for the Kronecker product:
 *   Tr_2(A (x) B) = Tr(B) * A      (trace out the RIGHT/E factor, keep LEFT)
 *   Tr_1(A (x) B) = Tr(A) * B      (trace out the LEFT factor,  keep RIGHT)
 * tests/test_mat.c asserts both on random A, B as the structural cross-check.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic/aic_mat.h"

void aic_mat_kronecker(acb_mat_t out, const acb_mat_t A, const acb_mat_t B,
                       slong prec)
{
    slong m = acb_mat_nrows(A), n = acb_mat_ncols(A);
    slong p = acb_mat_nrows(B), q = acb_mat_ncols(B);

    assert(acb_mat_nrows(out) == m * p && acb_mat_ncols(out) == n * q &&
           "aic_mat_kronecker: out must be (m*p) x (n*q)");

    /* out[i*p + r, j*q + s] = A[i,j] * B[r,s]; LEFT factor (A, indices i,j)
     * major. */
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
            for (slong r = 0; r < p; r++)
                for (slong s = 0; s < q; s++)
                    acb_mul(acb_mat_entry(out, i * p + r, j * q + s),
                            acb_mat_entry(A, i, j),
                            acb_mat_entry(B, r, s), prec);
}

void aic_mat_partial_trace_right(acb_mat_t out, const acb_mat_t M,
                                 slong a, slong b, slong prec)
{
    assert(acb_mat_nrows(M) == a * b && acb_mat_ncols(M) == a * b &&
           "aic_mat_partial_trace_right: M must be (a*b) x (a*b)");
    assert(acb_mat_nrows(out) == a && acb_mat_ncols(out) == a &&
           "aic_mat_partial_trace_right: out must be a x a");

    /* (Tr_2 M)[i,k] = sum_j M[i*b + j, k*b + j]; trace out the RIGHT factor. */
    for (slong i = 0; i < a; i++) {
        for (slong k = 0; k < a; k++) {
            acb_ptr dst = acb_mat_entry(out, i, k);
            acb_zero(dst);
            for (slong j = 0; j < b; j++)
                acb_add(dst, dst, acb_mat_entry(M, i * b + j, k * b + j), prec);
        }
    }
}

void aic_mat_partial_trace_left(acb_mat_t out, const acb_mat_t M,
                                slong a, slong b, slong prec)
{
    assert(acb_mat_nrows(M) == a * b && acb_mat_ncols(M) == a * b &&
           "aic_mat_partial_trace_left: M must be (a*b) x (a*b)");
    assert(acb_mat_nrows(out) == b && acb_mat_ncols(out) == b &&
           "aic_mat_partial_trace_left: out must be b x b");

    /* (Tr_1 M)[j,l] = sum_i M[i*b + j, i*b + l]; trace out the LEFT factor. */
    for (slong j = 0; j < b; j++) {
        for (slong l = 0; l < b; l++) {
            acb_ptr dst = acb_mat_entry(out, j, l);
            acb_zero(dst);
            for (slong i = 0; i < a; i++)
                acb_add(dst, dst, acb_mat_entry(M, i * b + j, i * b + l), prec);
        }
    }
}
