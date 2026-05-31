/* aic_latd_svd.c — singular values, operator norm, and full SVD on the LAPACK
 * double path via LAPACKE_zgesvd (beads aic-w4o.5; double-path coverage of the
 * SVD that aic-w4o.2 wants). See include/aic_latd.h for the contract.
 *
 * Route (Law 2: wrap LAPACK). LAPACKE_zgesvd is the LAPACK divide-by-QR SVD
 * driver: A = U * Sigma * V^dag with Sigma's diagonal (the singular values) in
 * DESCENDING order, U (m x m) the left singular vectors as columns, and Vt
 * (n x n) = V^dag whose ROWS are the conjugated right singular vectors. The
 * SVD is well-conditioned even when the eigenproblem of A is hypersensitive
 * (the non-normal t^{1/3} adversary, nla.md Family 1a/2): singular values are
 * the eigenvalues of the Hermitian PSD A^dag A and inherit its Weyl stability,
 * so the double SVD agrees with arb@53 on operator norm even where the SPECTRUM
 * does not — exactly the cross-check test_latd.c exercises.
 *
 * Signature (grepped from /usr/include/lapacke.h:906):
 *   lapack_int LAPACKE_zgesvd(int matrix_layout, char jobu, char jobvt,
 *       lapack_int m, lapack_int n, lapack_complex_double *a, lapack_int lda,
 *       double *s, lapack_complex_double *u, lapack_int ldu,
 *       lapack_complex_double *vt, lapack_int ldvt, double *superb);
 * matrix_layout = LAPACK_ROW_MAJOR; `a` is destroyed on exit so we copy first;
 * `superb` (length min(m,n)-1) returns the unconverged superdiagonal on info>0.
 * For row-major, lda = n (cols of A), ldu = m, ldvt = n.
 *
 * The operator norm is sigma_max = s[0] (descending order), computed with
 * jobu='N', jobvt='N' (values only — no vectors allocated).
 *
 * Fail-loud (Rule 4): info != 0 (bad arg or non-convergence) aborts.
 */
#define LAPACK_COMPLEX_C99 /* pin lapack_complex_double = C99 double _Complex */
#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include <lapacke.h>

#include "aic/aic_latd.h"

/* Shared core: copy A (zgesvd destroys its input), call zgesvd with the given
 * job chars and (possibly NULL) U/Vt buffers, assert info==0. */
static void latd_svd_core(double *svals, char jobu, char jobvt,
                          double _Complex *U, double _Complex *Vt,
                          const double _Complex *A_arr, slong m, slong n)
{
    assert(m >= 1 && n >= 1 && "SVD: dimensions must be positive");
    assert(svals != NULL);

    slong mn = (m < n) ? m : n;

    double _Complex *a = malloc((size_t) (m * n) * sizeof(double _Complex));
    double *superb = malloc((size_t) (mn > 1 ? mn - 1 : 1) * sizeof(double));
    assert(a != NULL && superb != NULL && "SVD: out of memory");
    for (slong k = 0; k < m * n; k++) a[k] = A_arr[k];

    /* Row-major leading dims: lda = n, ldu = m, ldvt = n. When a vector job is
     * 'N' LAPACKE ignores the corresponding pointer/ld, so passing NULL/1 is
     * safe; we pass the real ld when the buffer is present. */
    lapack_int ldu = (jobu == 'N') ? 1 : (lapack_int) m;
    lapack_int ldvt = (jobvt == 'N') ? 1 : (lapack_int) n;

    lapack_int info = LAPACKE_zgesvd(LAPACK_ROW_MAJOR, jobu, jobvt,
                                     (lapack_int) m, (lapack_int) n,
                                     a, (lapack_int) n, svals,
                                     U, ldu, Vt, ldvt, superb);
    if (info != 0) {
        fprintf(stderr,
                "latd_svd_core: LAPACKE_zgesvd failed, info=%d "
                "(info<0: bad arg; info>0: %d superdiagonals did not "
                "converge)\n", (int) info, (int) info);
        abort();
    }

    free(superb);
    free(a);
}

double aic_latd_opnorm(const double _Complex *A_arr, slong m, slong n)
{
    slong mn = (m < n) ? m : n;
    double *svals = malloc((size_t) mn * sizeof(double));
    assert(svals != NULL && "opnorm: out of memory");

    latd_svd_core(svals, 'N', 'N', NULL, NULL, A_arr, m, n);
    double sigma_max = svals[0]; /* descending order: largest is first */

    free(svals);
    return sigma_max;
}

void aic_latd_singular_values(double *svals, const double _Complex *A_arr,
                              slong m, slong n)
{
    latd_svd_core(svals, 'N', 'N', NULL, NULL, A_arr, m, n);
}

void aic_latd_svd(double *svals, double _Complex *U, double _Complex *Vt,
                  const double _Complex *A_arr, slong m, slong n)
{
    char jobu = (U != NULL) ? 'A' : 'N';
    char jobvt = (Vt != NULL) ? 'A' : 'N';
    latd_svd_core(svals, jobu, jobvt, U, Vt, A_arr, m, n);
}
