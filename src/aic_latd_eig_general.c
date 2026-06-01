/* aic_latd_eig_general.c — general (non-Hermitian) eigendecomposition and two
 * spectral-gap helpers on the LAPACK double path via LAPACKE_zgeev (bead
 * aic-pvs). See include/aic_latd.h for the contract.
 *
 * Route (Law 2: wrap LAPACK, do not reinvent). LAPACKE_zgeev is the LAPACK
 * driver for the GENERAL (unsymmetric/non-Hermitian) eigenproblem: it balances
 * A, reduces to upper Hessenberg, and runs the implicit-shift QR algorithm,
 * returning the n complex eigenvalues w and (with jobvr='V') the right
 * eigenvectors vr (A v_k = w_k v_k), each Euclidean-normalized with its largest
 * component real. It is the non-Hermitian counterpart of the zheev path
 * (aic_latd_eig.c) and mirrors that file's signature/style/allocation contract.
 *
 * UNCERTIFIED + ILL-CONDITIONED — read this before trusting a value. The whole
 * reason the project prefers Hermitian/normal constructions and certified arb
 * balls is that for a NON-NORMAL A the spectrum is perturbation-sensitive: a
 * tiny entry can move an eigenvalue by O(t^{1/3}) (the paper's tex:540 3x3
 * companion, CLAUDE.md "Spectra are perturbation-sensitive" callout). zgeev in
 * double STILL RETURNS FINITE eigenvalues for such A, but their accuracy can be
 * only a few digits — this is the fast/uncertified anchor, NOT a certified
 * result. Certified eigenvalue ENCLOSURES are the separate arb concern (bead
 * aic-w4o.1, acb_mat_eig_*). The tests therefore assert TIGHT accuracy only
 * where the spectrum is well conditioned (Hermitian input, exact triangular
 * spectra), and merely finiteness + rough magnitude on the non-normal adversary.
 *
 * Signature (grepped from /usr/include/lapacke.h:538):
 *   lapack_int LAPACKE_zgeev(int matrix_layout, char jobvl, char jobvr,
 *       lapack_int n, lapack_complex_double *a, lapack_int lda,
 *       lapack_complex_double *w, lapack_complex_double *vl, lapack_int ldvl,
 *       lapack_complex_double *vr, lapack_int ldvr);
 * matrix_layout = LAPACK_ROW_MAJOR (our convention); jobvl='N' (no LEFT vectors
 * — we never need them); jobvr='V' iff right eigenvectors are wanted else 'N'.
 *
 * ROW-MAJOR HANDLING (load-bearing; the zheev cross-check is its proof). `a` is
 * destroyed by zgeev, so we copy first. The eigenvalues `w` are a plain vector,
 * layout-independent. For the eigenvectors, LAPACKE_zgeev with
 * LAPACK_ROW_MAJOR returns `vr` ALREADY in our row-major (i,j) layout: column k
 * is the k-th right eigenvector, vr[i*n + k] = (v_k)_i, paired with w[k] — the
 * IDENTICAL column-is-eigenvector convention the Hermitian path documents
 * (aic_latd.h). LAPACKE performs the transpose to/from Fortran column-major
 * internally; we pass ldvr = n (cols) and read column k as a[i*n+k], with NO
 * manual transpose. The Hermitian-vs-zheev test confirms this: for a Hermitian
 * input, zgeev's eigenpairs (lambda_k, vr column k) satisfy A v = lambda v and
 * the eigenvalue SET matches the independent zheev path — which it could not if
 * the layout were transposed.
 *
 * Fail-loud (Rule 4): info != 0 (info<0 bad arg; info>0: the QR algorithm failed
 * to compute all eigenvalues, elements info+1..n of w did converge) aborts — no
 * silent partial spectrum.
 */
#define LAPACK_COMPLEX_C99 /* pin lapack_complex_double = C99 double _Complex */
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lapacke.h>

#include "aic/aic_latd.h"

void aic_latd_eig_general(double _Complex *evals, double _Complex *evecs,
                          const double _Complex *A_arr, slong n)
{
    assert(n >= 1 && "general eig: dimension must be positive");
    assert(evals != NULL);

    /* zgeev destroys `a`, so work on a copy to leave A_arr intact. */
    double _Complex *a = malloc((size_t) (n * n) * sizeof(double _Complex));
    assert(a != NULL && "general eig: out of memory");
    for (slong k = 0; k < n * n; k++) a[k] = A_arr[k];

    char jobvr = (evecs != NULL) ? 'V' : 'N';
    /* vr ld must be valid even when unused; ldvr=1 when jobvr='N' (LAPACKE
     * ignores the pointer/ld then), real ld=n when the buffer is present. */
    lapack_int ldvr = (jobvr == 'N') ? 1 : (lapack_int) n;

    lapack_int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', jobvr,
                                    (lapack_int) n, a, (lapack_int) n, evals,
                                    NULL, 1, evecs, ldvr);
    if (info != 0) {
        fprintf(stderr,
                "aic_latd_eig_general: LAPACKE_zgeev failed, info=%d "
                "(info<0: bad arg; info>0: QR failed, only eigenvalues "
                "%d+1..%ld converged)\n",
                (int) info, (int) info, (long) n);
        abort();
    }

    free(a);
}

double aic_latd_spectral_gap(const double _Complex *evals, slong n,
                             double _Complex z)
{
    assert(n >= 1 && "spectral gap: need at least one eigenvalue");
    assert(evals != NULL);
    double best = cabs(evals[0] - z);
    for (slong i = 1; i < n; i++) {
        double d = cabs(evals[i] - z);
        if (d < best) best = d;
    }
    return best;
}

double aic_latd_spectral_separation(const double _Complex *evals, slong n)
{
    assert(evals != NULL || n < 1);
    if (n < 2) return INFINITY; /* no pair: separation is vacuously infinite */
    double best = INFINITY;
    for (slong i = 0; i < n; i++) {
        for (slong j = i + 1; j < n; j++) {
            double d = cabs(evals[i] - evals[j]);
            if (d < best) best = d;
        }
    }
    return best;
}
