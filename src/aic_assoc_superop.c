/* aic_assoc_superop.c — the superoperator matrix S_Phi of a UCP self-map, its
 * apply, and the Kronecker cross-check (bead aic-92f, module "assoc_ecsa",
 * Increment 1). This is the arb-path build engine for the regularization
 * Phi_tilde = theta(2 Phi - 1); see include/aic_assoc.h for the full why.
 *
 * approximate_algebras.tex:2171-2174 — the regularization to be built on top:
 *     Phi_tilde = theta(2 Phi - 1) = (1/2)(1 + sgn(2 Phi - 1)).
 * The functional calculus needs Phi as an n^2 x n^2 MATRIX; this file produces it.
 *
 * SUPEROPERATOR / vec CONVENTION (row-major, matches aic_idemp.h / aic_ecstar.h).
 *   vec_r(X)[i*n+j] = X[i,j]; S[a*n+b, p*n+q] = T(E_{pq})[a,b]; so column (p*n+q)
 *   of S is vec_r(T(E_{pq})). The PRIMARY build is column-by-column via the tested
 *   aic_ucp_apply on the matrix units E_{pq} (lowest vec-bug risk, the idemp/
 *   ecstar precedent). The SECONDARY build (aic_assoc_superop_kron) uses the
 *   index-derived formula S_Phi = sum_a K_a^dag (x) K_a^T, pinning the convention
 *   a second, independent way (test T2 asserts the two agree).
 *
 * WHY the Kronecker formula has THAT shape. Phi(X) = sum_a K_a^dag X K_a, so
 *   Phi(X)[a,b] = sum_a sum_{p,q} conj(K_a[p,a]) X[p,q] K_a[q,b]
 *               = sum_a sum_{p,q} (K_a^dag)[a,p] X[p,q] (K_a^T)[q,b].
 *   With vec_r, S[a*n+b, p*n+q] = sum_a (K_a^dag)[a,p] (K_a^T)[b,q]
 *   = sum_a (K_a^dag (x) K_a^T)[a*n+b, p*n+q] in the left-major Kronecker
 *   (A (x) C)[i*n+j, p*n+q] = A[i,p] C[j,q] of aic_mat_kronecker. Hence
 *   S_Phi = sum_a K_a^dag (x) K_a^T. (Oracle-tested anyway, T1/T2 — the project
 *   hard rule on the vec convention.)
 */
#include <assert.h>

#include <flint/acb_mat.h>

#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic/aic_assoc.h"

/* Column (p*n+q) of S = vec_r(Phi(E_{pq})), built via aic_ucp_apply. */
void aic_assoc_superop_from_ucp(acb_mat_t S, const aic_ucp_kraus *phi,
                                slong prec)
{
    slong n = phi->dim_H;
    assert(phi->dim_K == n &&
           "aic_assoc_superop_from_ucp: Phi must be a self-map (dim_K==dim_H)");
    assert(acb_mat_nrows(S) == n * n && acb_mat_ncols(S) == n * n &&
           "aic_assoc_superop_from_ucp: S must be n^2 x n^2");

    acb_mat_t E, img;
    acb_mat_init(E, n, n);      /* matrix unit E_{pq} */
    acb_mat_init(img, n, n);    /* Phi(E_{pq}) */
    for (slong p = 0; p < n; p++)
        for (slong q = 0; q < n; q++) {
            acb_mat_zero(E);
            acb_set_si(acb_mat_entry(E, p, q), 1);   /* E_{pq} = |e_p><e_q| */
            aic_ucp_apply(img, phi, E, prec);        /* Phi(E_{pq}) */
            /* column (p*n+q): S[a*n+b, p*n+q] = Phi(E_{pq})[a,b] (vec_r). */
            for (slong a = 0; a < n; a++)
                for (slong b = 0; b < n; b++)
                    acb_set(acb_mat_entry(S, a * n + b, p * n + q),
                            acb_mat_entry(img, a, b));
        }
    acb_mat_clear(img);
    acb_mat_clear(E);
}

/* S_Phi = sum_a K_a^dag (x) K_a^T (left-major Kronecker; file docstring). */
void aic_assoc_superop_kron(acb_mat_t S, const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H;
    assert(phi->dim_K == n &&
           "aic_assoc_superop_kron: Phi must be a self-map (dim_K==dim_H)");
    assert(acb_mat_nrows(S) == n * n && acb_mat_ncols(S) == n * n &&
           "aic_assoc_superop_kron: S must be n^2 x n^2");

    acb_mat_t Kd, Kt, term;
    acb_mat_init(Kd, n, n);     /* K_a^dag */
    acb_mat_init(Kt, n, n);     /* K_a^T   */
    acb_mat_init(term, n * n, n * n);
    acb_mat_zero(S);
    for (slong a = 0; a < phi->r; a++) {
        acb_mat_conjugate_transpose(Kd, phi->K[a]);
        acb_mat_transpose(Kt, phi->K[a]);
        aic_mat_kronecker(term, Kd, Kt, prec);   /* K_a^dag (x) K_a^T */
        acb_mat_add(S, S, term, prec);
    }
    acb_mat_clear(term);
    acb_mat_clear(Kt);
    acb_mat_clear(Kd);
}

/* out = reshape(S vec_r(X)): out[a,b] = sum_{p,q} S[a*n+b, p*n+q] X[p,q]. */
void aic_assoc_superop_apply(acb_mat_t out, const acb_mat_t S, const acb_mat_t X,
                             slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(acb_mat_ncols(X) == n && "aic_assoc_superop_apply: X must be square");
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_assoc_superop_apply: out must be n x n");
    assert(acb_mat_nrows(S) == n * n && acb_mat_ncols(S) == n * n &&
           "aic_assoc_superop_apply: S must be n^2 x n^2");
    assert(out != X && "aic_assoc_superop_apply: out must not alias X");

    acb_t acc, term;
    acb_init(acc);
    acb_init(term);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++) {
            acb_zero(acc);
            for (slong p = 0; p < n; p++)
                for (slong q = 0; q < n; q++) {
                    acb_mul(term, acb_mat_entry(S, a * n + b, p * n + q),
                            acb_mat_entry(X, p, q), prec);
                    acb_add(acc, acc, term, prec);
                }
            acb_set(acb_mat_entry(out, a, b), acc);
        }
    acb_clear(term);
    acb_clear(acc);
}
