/* aic_latd.h — the LAPACK double path ("latd"), an UNCERTIFIED fast anchor
 * (beads aic-w4o.5). This is the project's second number path (CLAUDE.md
 * "Architecture": double path = fast prec=53 anchor, arb path = certified error
 * balls). It is an AUDITION CANDIDATE (Law 4), NOT a replacement for the arb
 * path: every routine here has an arb counterpart in include/aic_mat.h, and the
 * cross-check that justifies this path is double(LAPACK) vs arb@prec=53 (the
 * su2-fft anchor; CLAUDE.md cross-check ladder rung 2).
 *
 * WHY a double path at all (Law 2: reuse LAPACK, do not hand-roll). FLINT's
 * certified eig path needs a SIMPLE (non-degenerate) spectrum: acb_mat_eig_simple
 * PROVES n isolated eigenvalues and aborts on a repeated eigenvalue (bead
 * aic-w4o.1; aic_mat_eig_hermitian's contract). But the project's real workload
 * — projections P, multiplicities, block-diagonal algebras — is exactly
 * degenerate (a rank-k projector has spectrum {0,1} with multiplicities). For
 * those, LAPACKE_zheev is the fast route: it diagonalises ANY Hermitian matrix,
 * degenerate or not, returning real eigenvalues and an orthonormal eigenbasis.
 * The price is no certified error ball — hence "uncertified". The two paths are
 * complementary: arb for certified O(eps) bounds, LAPACK for fast/degenerate.
 *
 * Number type. LAPACKE's complex scalar is C99 `double _Complex`
 * (lapack_complex_double); src/aic_latd.c pins this with LAPACK_COMPLEX_C99
 * before <lapacke.h> so the ABI is the plain C complex, not the struct form.
 *
 * STORAGE CONVENTION (load-bearing — pinned here). The double path uses
 * LAPACK_ROW_MAJOR throughout: a complex matrix is a flat `double _Complex`
 * array `A` of length rows*cols with
 *     A[i*cols + j] = entry (i, j).
 * This matches the C-natural row-major layout and acb_mat's (i,j) indexing, so
 * the acb_mat <-> array conversions below are a straight transcription with no
 * transpose. LAPACKE_zheev / LAPACKE_zgesvd are all called with
 * LAPACK_ROW_MAJOR and the leading dimension = number of columns.
 *
 * Midpoints only. acb_mat_t -> array takes the BALL MIDPOINT of each entry
 * (the radius is discarded; this is the uncertified path). array -> acb_mat_t
 * writes each double as a ZERO-RADIUS acb (the double value is exact-as-given;
 * no spurious error ball is invented).
 *
 * Fail-loud (Rule 4). Dimensions are asserted on every conversion, and every
 * LAPACKE call's `info` return is asserted == 0 (a nonzero info — failed
 * convergence or a bad argument — aborts, never a silent wrong answer).
 */
#ifndef AIC_LATD_H
#define AIC_LATD_H

#include <complex.h>
#include <flint/acb_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --- conversion helpers (row-major, midpoints; see file docstring) --- */

/* acb_mat_t A (r x c) -> caller-allocated row-major `out` of length r*c,
 * out[i*c + j] = midpoint of A(i,j). `out` must point to >= r*c elements
 * (NOT checked — it is a raw buffer); r, c are taken from A. */
void aic_latd_from_acb_mat(double _Complex *out, const acb_mat_t A);

/* row-major `in` of length r*c (in[i*c+j] = entry (i,j)) -> acb_mat_t out.
 * `out` must be initialised to exactly r x c (asserted). Each entry becomes a
 * zero-radius acb (the double is taken as exact). */
void aic_latd_to_acb_mat(acb_mat_t out, const double _Complex *in,
                         slong r, slong c);

/* --- Hermitian eigendecomposition via LAPACKE_zheev (degeneracy-robust) --- */

/* Eigenvalues (and optionally eigenvectors) of an n x n Hermitian H, double
 * path. HANDLES DEGENERATE SPECTRA (repeated eigenvalues) — the whole reason
 * this path exists (the arb simple-eig path cannot, bead aic-w4o.1).
 *   H_arr : row-major n*n input (only the upper triangle is read by zheev with
 *           uplo='U'; the rest is ignored, so a Hermitian array need not be
 *           perfectly symmetric in its stored lower half). DESTROYED-free: a
 *           local copy is made, H_arr is not modified.
 *   evals : caller-allocated double[n], the real eigenvalues ASCENDING (zheev's
 *           native order).
 *   evecs : if non-NULL, caller-allocated row-major n*n; column k (evecs[i*n+k])
 *           is the orthonormal eigenvector for evals[k]. If NULL, jobz='N'
 *           (eigenvalues only, faster).
 * Asserts n >= 1 and the LAPACKE info == 0 (fail loud). */
void aic_latd_eig_hermitian(double *evals, double _Complex *evecs,
                            const double _Complex *H_arr, slong n);

/* --- operator norm / SVD via LAPACKE_zgesvd --- */

/* Operator norm = largest singular value of an m x n matrix A (row-major),
 * double path. Calls the full SVD internally with jobu='N', jobvt='N' (singular
 * values only). Returns sigma_max as a double. Asserts m,n >= 1 and info == 0. */
double aic_latd_opnorm(const double _Complex *A_arr, slong m, slong n);

/* Singular values of an m x n matrix A (row-major), DESCENDING (zgesvd's native
 * order). svals : caller-allocated double[min(m,n)]. Asserts info == 0. */
void aic_latd_singular_values(double *svals, const double _Complex *A_arr,
                              slong m, slong n);

/* Full SVD A = U * Sigma * V^dag of an m x n matrix A (row-major). Outputs
 * (all caller-allocated, any of U/Vt may be NULL to skip):
 *   svals : double[min(m,n)], singular values DESCENDING.
 *   U     : row-major m*m, the left singular vectors as COLUMNS (jobu='A').
 *   Vt    : row-major n*n, V^dag — its ROWS are the conjugated right singular
 *           vectors (jobvt='A'), i.e. Vt(k, :) = v_k^dag, matching LAPACK.
 * So A(i,j) = sum_k U(i,k) * svals[k] * Vt(k,j). If both U and Vt are NULL this
 * reduces to singular values only. Asserts info == 0 (fail loud). */
void aic_latd_svd(double *svals, double _Complex *U, double _Complex *Vt,
                  const double _Complex *A_arr, slong m, slong n);

#ifdef __cplusplus
}
#endif

#endif /* AIC_LATD_H */
