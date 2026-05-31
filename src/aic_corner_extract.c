/* aic_corner_extract.c — extraction of the corner subspace S_{P,Q} = Img Co_{P,Q}
 * from the idempotent coordinate matrix Co (d x d, d = dim_A). The COORDINATE-SPACE
 * analogue of aic_assoc_extract_range (which extracts range of an n^2 x n^2 superop):
 * here Co is a d x d matrix, its top-dim_S LEFT singular vectors u_m are d-vectors,
 * and the corner basis is C_m = sum_k u_m[k] B_k (n x n OPERATORS). See
 * include/aic_corner.h for the full why.
 *
 * approximate_algebras.tex:1064 — S_{P,Q}:
 *   "The image of this map, S_{P,Q} = Img Co_{P,Q} = Ker(1 - Co_{P,Q}), is a
 *    closed linear subspace of A."
 *
 * WHY round(trace) == rank for an idempotent (mirrors assoc_extract). Co^2 = Co
 * (after the theta-snap, exact to machine precision) => its eigenvalues are 0 or
 * 1, so Tr Co = #(unit eigenvalues) = rank = dim S_{P,Q}. The trace is rounded;
 * the SVD-gap count (singular values > 0.5) is an INDEPENDENT rank witness. They
 * MUST agree (fail loud, Rule 4): a mismatch means Co is not cleanly idempotent.
 *
 * THE OBLIQUE-PROJECTOR SUBTLETY (same as assoc_extract). Co is idempotent but in
 * general NOT Hermitian (L_P, R_Q need not commute / be self-adjoint in coords).
 * For an idempotent the nonzero singular values are >= 1 (== 1 iff orthogonal
 * projector, > 1 iff oblique) and the zero ones ~0, so 0.5 cleanly separates
 * range from kernel, and the top-dim_S LEFT singular vectors span Img Co either
 * way. Extraction is the fast DOUBLE path (LAPACK zgesvd via aic_latd_svd; the
 * certified degenerate-SVD wall aic-w4o.1 is deferred); the basis is returned as
 * acb_mat at `prec`.
 *
 * COORDINATE -> OPERATOR (the bridge, .tex / corner.h). u_m is a d-vector in A's
 * {B_k} basis; C_m = sum_k u_m[k] B_k is the corresponding n x n operator. Since
 * {B_k} is Frobenius-orthonormal, <C_m, C_p>_F = sum_k conj(u_m[k]) u_p[k] =
 * <u_m, u_p> = delta_{mp} (left singular vectors orthonormal), so {C_m} is
 * Frobenius-orthonormal (ASSERTED, fail loud).
 */
#include <assert.h>
#include <complex.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic/aic_corner.h"
#include "aic_corner_internal.h"
#include "aic/aic_latd.h"

/* Apply a d x d coordinate-space operator T to X in A (the operator->coords->
 * operator bridge; lives here, the coordinate<->operator TU, to keep
 * src/aic_corner_compress.c under the ~200 LOC limit, Rule 10):
 *   out = sum_m (T x)_m B_m,   x_k = <B_k, X>_F  (aic_corner_frob_ip). */
void aic_corner_apply(acb_mat_t out, const aic_ecstar *A, const acb_mat_t T,
                      const acb_mat_t X, slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(out != X && "corner_apply: out must not alias X");
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "corner_apply: out must be n x n");
    assert(acb_mat_nrows(X) == n && "corner_apply: X must be n x n");
    assert(acb_mat_nrows(T) == d && acb_mat_ncols(T) == d &&
           "corner_apply: T must be dim_A x dim_A");

    /* x_k = <B_k, X>_F (coords of X in A). */
    acb_ptr x = _acb_vec_init(d);
    for (slong k = 0; k < d; k++)
        aic_corner_frob_ip(x + k, A->B[k], X, prec);

    /* y = T x, then out = sum_m y_m B_m. */
    acb_mat_zero(out);
    acb_t ym, t;
    acb_init(ym);
    acb_init(t);
    for (slong m = 0; m < d; m++) {
        acb_zero(ym);
        for (slong k = 0; k < d; k++) {
            acb_mul(t, acb_mat_entry(T, m, k), x + k, prec);
            acb_add(ym, ym, t, prec);
        }
        /* out += y_m B_m */
        for (slong a = 0; a < n; a++)
            for (slong b = 0; b < n; b++) {
                acb_mul(t, ym, acb_mat_entry(A->B[m], a, b), prec);
                acb_add(acb_mat_entry(out, a, b), acb_mat_entry(out, a, b), t, prec);
            }
    }
    acb_clear(t);
    acb_clear(ym);
    _acb_vec_clear(x, d);
}

/* round(Re Tr Co) as a slong (the idempotent rank = dim S_{P,Q}). Accumulates the
 * diagonal at the CALL's `prec` (not a hardcoded 64) for precision consistency. */
static slong trace_rank(const acb_mat_t Co, slong prec)
{
    slong d = acb_mat_nrows(Co);
    arb_t tr;
    arb_init(tr);
    arb_zero(tr);
    for (slong r = 0; r < d; r++)
        arb_add(tr, tr, acb_realref(acb_mat_entry(Co, r, r)), prec);
    double trd = arf_get_d(arb_midref(tr), ARF_RND_NEAR);
    arb_clear(tr);
    return (slong) llround(trd);
}

/* Count singular values above gap_thr (an idempotent's nonzero sigmas are >= 1,
 * zeros ~0, so 0.5 is the clean separator). */
static slong svd_rank(const double *sv, slong d, double gap_thr)
{
    slong r = 0;
    for (slong i = 0; i < d; i++) if (sv[i] > gap_thr) r++;
    return r;
}

void aic_corner_extract(acb_mat_t **C_out, slong *dim_S_out, const acb_mat_t Co,
                        const aic_ecstar *A, slong prec)
{
    slong d = acb_mat_nrows(Co);
    assert(acb_mat_ncols(Co) == d && "corner_extract: Co must be square");
    assert(d == A->dim_A && "corner_extract: Co dim != A->dim_A");
    slong n = A->n;

    /* double-path thin SVD of the d x d Co (oblique idempotent; LAPACK zgesvd). */
    double _Complex *Cd = flint_malloc((size_t) (d * d) * sizeof(double _Complex));
    double _Complex *U = flint_malloc((size_t) (d * d) * sizeof(double _Complex));
    double *sv = flint_malloc((size_t) d * sizeof(double));
    aic_latd_from_acb_mat(Cd, Co);                 /* midpoints, row-major */
    aic_latd_svd(sv, U, NULL, Cd, d, d);           /* U cols = left vectors */

    /* dim_S = round(trace), cross-checked against the SVD gap (fail loud). */
    slong dim_S = trace_rank(Co, prec);
    slong r_svd = svd_rank(sv, d, 0.5);
    assert(dim_S == r_svd &&
           "corner_extract: round(trace) != SVD-gap rank -- Co not cleanly "
           "idempotent (out of basin?) or gap threshold wrong");
    assert(dim_S >= 0 && dim_S <= d && "corner_extract: dim_S out of range");
    /* the gap must be genuine (else the rank is ambiguous). dim_S==0 (P~0) and
     * dim_S==d (P~I) are the legitimate degenerate ends, .tex:1066. */
    assert((dim_S == 0 || sv[dim_S - 1] > 0.5) &&
           (dim_S == d || sv[dim_S] < 0.5) &&
           "corner_extract: singular-value gap not at 0.5 (ill-conditioned rank)");

    if (dim_S == 0) {
        *C_out = NULL;
        *dim_S_out = 0;
        flint_free(sv);
        flint_free(U);
        flint_free(Cd);
        return;
    }

    /* C_m = sum_k u_m[k] B_k, u_m = column m of U (U[k*d + m]). */
    acb_mat_t *C = flint_malloc((size_t) dim_S * sizeof(acb_mat_t));
    acb_t coeff, t;
    acb_init(coeff);
    acb_init(t);
    for (slong m = 0; m < dim_S; m++) {
        acb_mat_init(C[m], n, n);
        acb_mat_zero(C[m]);
        for (slong k = 0; k < d; k++) {
            double _Complex z = U[k * d + m];
            acb_set_d_d(coeff, creal(z), cimag(z));
            /* C[m] += coeff * B_k */
            for (slong a = 0; a < n; a++)
                for (slong b = 0; b < n; b++) {
                    acb_mul(t, coeff, acb_mat_entry(A->B[k], a, b), prec);
                    acb_add(acb_mat_entry(C[m], a, b),
                            acb_mat_entry(C[m], a, b), t, prec);
                }
        }
    }
    acb_clear(t);
    acb_clear(coeff);

    /* ASSERT Frobenius-orthonormality: <C_j,C_m>_F = delta_{jm} to tol. {B_k}
     * orthonormal + u_m orthonormal => the C_m are orthonormal; a failure means a
     * corrupt SVD or a non-orthonormal {B_k}. Fail loud (Rule 4). */
    {
        arb_t worst, tol, mag;
        acb_t ip, tt;
        arb_init(worst); arb_init(tol); arb_init(mag);
        acb_init(ip); acb_init(tt);
        arb_zero(worst);
        for (slong j = 0; j < dim_S; j++)
            for (slong m = 0; m < dim_S; m++) {
                acb_zero(ip);
                for (slong a = 0; a < n; a++)
                    for (slong b = 0; b < n; b++) {
                        acb_conj(tt, acb_mat_entry(C[j], a, b));
                        acb_mul(tt, tt, acb_mat_entry(C[m], a, b), prec);
                        acb_add(ip, ip, tt, prec);
                    }
                if (j == m) acb_sub_si(ip, ip, 1, prec);
                acb_abs(mag, ip, prec);
                if (arb_gt(mag, worst)) arb_set(worst, mag);
            }
        arb_set_d(tol, 1e-9);   /* double-path SVD floors near 1e-12 */
        assert(arb_le(worst, tol) &&
               "corner_extract: C_m not Frobenius-orthonormal (corrupt SVD?)");
        acb_clear(tt); acb_clear(ip);
        arb_clear(mag); arb_clear(tol); arb_clear(worst);
    }

    *C_out = C;
    *dim_S_out = dim_S;
    flint_free(sv);
    flint_free(U);
    flint_free(Cd);
}
