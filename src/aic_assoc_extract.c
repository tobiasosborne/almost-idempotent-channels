/* aic_assoc_extract.c — A-extraction: the Frobenius-orthonormal operator basis
 * {B_k} of A = range(Phi_tilde) from the idempotent superop S_tilde (bead aic-92f,
 * module "assoc_ecsa", Increment 2). Route-B1 (research): dim_A = round(trace)
 * cross-checked against the thin-SVD gap; B_k = reshape(top-k left singular
 * vector). See include/aic_assoc.h for the full why (the oblique-projector subtlety).
 *
 * approximate_algebras.tex:2185-2186 — A:
 *     A = Img Phi_tilde = Ker(1 - Phi_tilde) <= B(H),  a closed subspace.
 *
 * WHY round(trace) == rank for an idempotent. S_tilde^2 = S_tilde => its
 * eigenvalues are 0 or 1, so Tr S_tilde = #(unit eigenvalues) = rank = dim A (an
 * integer). The trace is computed and rounded; the SVD gap is an INDEPENDENT
 * rank witness (count of singular values above a gap threshold). They MUST agree
 * (fail loud, Rule 4) -- a disagreement means S_tilde is not cleanly idempotent
 * (eta too large, out of basin) or the gap threshold is wrong for this instance.
 *
 * THE PROJECTOR-SVD SUBTLETY. S_tilde is idempotent (S_tilde^2 = S_tilde) but in
 * general NOT Hermitian (Phi_tilde is HP but not always HS-self-adjoint). For an
 * idempotent the NONZERO singular values are >= 1, with EQUALITY (all == 1, an
 * ORTHOGONAL projector, SVD spectrum exactly {1,...,1,0,...,0}) when Phi_tilde IS
 * HS-self-adjoint -- the self-dual channels: dephasing, AND the dep(d)+conj(dep(d))
 * eta>0 family (a convex mix of HS-self-adjoint maps is HS-self-adjoint, so its
 * S_tilde is ORTHOGONAL, sigma_max == 1; measured, test U5). The sigmas are
 * STRICTLY > 1 only when Phi_tilde is genuinely OBLIQUE (non-HS-self-adjoint),
 * e.g. the compression channel compress_idemp(4,2), whose S_tilde inflates
 * sigma_max to sqrt(3) ~ 1.732 (measured, U5). Either way range(S_tilde) = span of
 * the top-dim_A LEFT singular vectors, and rank = trace still holds, so the 0.5 gap
 * (nonzero sigmas >= 1, zeros ~0) is a clean separator regardless. Extraction is
 * the fast DOUBLE path (LAPACK zgesvd via aic_latd_svd; certified degenerate
 * eig/SVD is the deferred aic-w4o.1 wall); the basis is returned as acb_mat at
 * `prec` and its Frobenius-orthonormality is asserted (left singular vectors are
 * orthonormal in C^{n^2} = Frobenius-orthonormal as operators).
 *
 * vec CONVENTION (row-major, matches aic_assoc.h / aic_idemp.h). vec_r(X)[i*n+j]
 * = X[i,j]; so the n^2-vector that is the m-th left singular vector reshapes to
 * B[i,j] = u_m[i*n+j].
 */
#include <assert.h>
#include <complex.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_assoc.h"
#include "aic_latd.h"

/* round(Re Tr S_tilde) as a slong (the idempotent rank). */
static slong trace_rank(const acb_mat_t S)
{
    slong nn = acb_mat_nrows(S);
    arb_t tr;
    arb_init(tr);
    arb_zero(tr);
    for (slong r = 0; r < nn; r++)
        arb_add(tr, tr, acb_realref(acb_mat_entry(S, r, r)), 64);
    double trd = arf_get_d(arb_midref(tr), ARF_RND_NEAR);
    arb_clear(tr);
    return (slong) llround(trd);
}

/* Count singular values above gap_thr (0.5: an idempotent's nonzero sigmas are
 * >= 1 -- == 1 for an orthogonal projector, > 1 for an oblique one -- and the zero
 * ones ~0, so 0.5 is the clean gap either way). */
static slong svd_rank(const double *sv, slong nn, double gap_thr)
{
    slong r = 0;
    for (slong i = 0; i < nn; i++) if (sv[i] > gap_thr) r++;
    return r;
}

void aic_assoc_extract_range(acb_mat_t **B_out, slong *dim_A_out,
                             const acb_mat_t S_tilde, slong prec)
{
    slong nn = acb_mat_nrows(S_tilde);
    assert(acb_mat_ncols(S_tilde) == nn && "extract_range: S_tilde must be square");
    /* n = sqrt(nn) (S_tilde is n^2 x n^2). */
    slong n = (slong) llround(sqrt((double) nn));
    assert(n * n == nn && "extract_range: S_tilde dim must be a perfect square");

    /* double-path thin SVD of S_tilde (oblique projector; LAPACK zgesvd). */
    double _Complex *Sd = flint_malloc((size_t) (nn * nn) * sizeof(double _Complex));
    double _Complex *U = flint_malloc((size_t) (nn * nn) * sizeof(double _Complex));
    double *sv = flint_malloc((size_t) nn * sizeof(double));
    aic_latd_from_acb_mat(Sd, S_tilde);             /* midpoints, row-major */
    aic_latd_svd(sv, U, NULL, Sd, nn, nn);          /* U cols = left vectors */

    /* dim_A = round(trace), cross-checked against the SVD gap (fail loud). */
    slong dim_A = trace_rank(S_tilde);
    slong r_svd = svd_rank(sv, nn, 0.5);
    assert(dim_A == r_svd &&
           "extract_range: round(trace) != SVD-gap rank -- S_tilde not cleanly "
           "idempotent (eta out of basin?) or gap threshold wrong");
    assert(dim_A >= 1 && dim_A <= nn && "extract_range: dim_A out of range");
    /* the gap must be genuine: sigma[dim_A-1] >> sigma[dim_A] (else the rank is
     * ambiguous; an idempotent's nonzero sigmas are >= 1 (== 1 orthogonal, > 1
     * oblique), zeros ~0). */
    assert(sv[dim_A - 1] > 0.5 &&
           (dim_A == nn || sv[dim_A] < 0.5) &&
           "extract_range: singular-value gap not at 0.5 (ill-conditioned rank)");

    /* B_k = reshape_row-major(column k of U): B_k[i,j] = U[(i*n+j)*nn + k]. */
    acb_mat_t *B = flint_malloc((size_t) dim_A * sizeof(acb_mat_t));
    for (slong k = 0; k < dim_A; k++) {
        acb_mat_init(B[k], n, n);
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++) {
                double _Complex z = U[(i * n + j) * nn + k];
                acb_set_d_d(acb_mat_entry(B[k], i, j), creal(z), cimag(z));
            }
    }

    /* ASSERT Frobenius-orthonormality of {B_k}: <B_j,B_k>_F = delta_jk to tol.
     * Left singular vectors are orthonormal in C^{nn}; reshaping preserves the
     * inner product, so a failure means a corrupt SVD. Fail loud (Rule 4). */
    {
        arb_t worst, tol, mag;
        acb_t ip, t;
        arb_init(worst); arb_init(tol); arb_init(mag);
        acb_init(ip); acb_init(t);
        arb_zero(worst);
        for (slong j = 0; j < dim_A; j++)
            for (slong k = 0; k < dim_A; k++) {
                acb_zero(ip);
                for (slong p = 0; p < n; p++)
                    for (slong q = 0; q < n; q++) {
                        acb_conj(t, acb_mat_entry(B[j], p, q));
                        acb_mul(t, t, acb_mat_entry(B[k], p, q), prec);
                        acb_add(ip, ip, t, prec);
                    }
                if (j == k) acb_sub_si(ip, ip, 1, prec);
                acb_abs(mag, ip, prec);
                if (arb_gt(mag, worst)) arb_set(worst, mag);
            }
        arb_set_d(tol, 1e-9);   /* double-path SVD floors near 1e-12; corrupt O(1) */
        assert(arb_le(worst, tol) &&
               "extract_range: B_k not Frobenius-orthonormal (corrupt SVD)");
        acb_clear(t); acb_clear(ip);
        arb_clear(mag); arb_clear(tol); arb_clear(worst);
    }

    *B_out = B;
    *dim_A_out = dim_A;
    flint_free(sv);
    flint_free(U);
    flint_free(Sd);
}
