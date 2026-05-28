/* aic_latd_eig.c — Hermitian eigendecomposition on the LAPACK double path via
 * LAPACKE_zheev (beads aic-w4o.5, addressing the degenerate-spectrum gap of
 * aic-w4o.1). See include/aic_latd.h for the contract.
 *
 * Route (Law 2: wrap LAPACK, do not reinvent). LAPACKE_zheev is the LAPACK
 * driver for the Hermitian eigenproblem: it reduces H to real tridiagonal form
 * (Householder) and applies implicit-QL/QR, returning REAL eigenvalues ascending
 * and (with jobz='V') an ORTHONORMAL eigenbasis. Unlike FLINT's certified
 * acb_mat_eig_simple, zheev does NOT require a simple spectrum — it diagonalises
 * a rank-k projector (spectrum {0,1}, multiplicities k and n-k) just as happily
 * as a generic matrix. That degeneracy-robustness is the whole point of this
 * path (the arb certified simple-eig path aborts on repeated eigenvalues).
 *
 * Signature (grepped from /usr/include/lapacke.h:1832):
 *   lapack_int LAPACKE_zheev(int matrix_layout, char jobz, char uplo,
 *       lapack_int n, lapack_complex_double *a, lapack_int lda, double *w);
 * matrix_layout = LAPACK_ROW_MAJOR (our convention); uplo = 'U' (read the upper
 * triangle — for an exactly Hermitian input the lower triangle is redundant);
 * jobz = 'V' if eigenvectors are wanted else 'N'. On exit, with jobz='V', `a` is
 * OVERWRITTEN by the eigenvectors (column k = eigenvector for w[k]); `w` holds
 * the eigenvalues ascending. We copy into a scratch buffer so the caller's input
 * is not destroyed, and (row-major) read column k as a[i*n + k].
 *
 * Fail-loud (Rule 4): info != 0 means zheev failed to converge (info > 0: that
 * many off-diagonal elements did not converge) or had a bad argument (info < 0).
 * Either is a hard error here — abort, never a silent partial result.
 */
#define LAPACK_COMPLEX_C99 /* pin lapack_complex_double = C99 double _Complex */
#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include <lapacke.h>

#include "aic_latd.h"

void aic_latd_eig_hermitian(double *evals, double _Complex *evecs,
                            const double _Complex *H_arr, slong n)
{
    assert(n >= 1 && "Hermitian eig: dimension must be positive");
    assert(evals != NULL);

    /* zheev overwrites `a` (with eigenvectors if jobz='V', or junk if 'N'), so
     * work on a copy to leave H_arr intact. */
    double _Complex *a = malloc((size_t) (n * n) * sizeof(double _Complex));
    assert(a != NULL && "Hermitian eig: out of memory");
    for (slong k = 0; k < n * n; k++) a[k] = H_arr[k];

    char jobz = (evecs != NULL) ? 'V' : 'N';
    lapack_int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, 'U',
                                    (lapack_int) n, a, (lapack_int) n, evals);
    if (info != 0) {
        fprintf(stderr,
                "aic_latd_eig_hermitian: LAPACKE_zheev failed, info=%d "
                "(info<0: bad arg; info>0: %d off-diagonals did not converge)\n",
                (int) info, (int) info);
        abort();
    }

    if (evecs != NULL) {
        /* a is now the row-major eigenvector matrix: column k is the
         * eigenvector for evals[k]. Copy straight across (same layout). */
        for (slong k = 0; k < n * n; k++) evecs[k] = a[k];
    }

    free(a);
}
