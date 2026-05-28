/* aic_mat_internal.h — shared internals for the mat module (beads aic-9zh,
 * "T-mat"). NOT part of the public API (include/aic_mat.h); these helpers are
 * shared between aic_mat_spectral.c (degeneracy-robust max-eig) and
 * aic_mat_eig.c (full Hermitian eigendecomposition) so the Hermiticity
 * precondition and the certified two-stage eig path are written once.
 */
#ifndef AIC_MAT_INTERNAL_H
#define AIC_MAT_INTERNAL_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

/* prec-appropriate absolute tol = 2^-(prec - 8); the 8-bit guard leaves headroom
 * so a faithfully-rounded Hermitian input (radii ~2^-prec) is not rejected. */
void aic_mat_int_tol(arb_t tol, slong prec);

/* Abort (fail loud, CLAUDE.md Rule 4) unless |H[i,j] - conj(H[j,i])| is
 * certainly <= tol for all i,j. Asserts squareness too. */
void aic_mat_int_assert_hermitian(const acb_mat_t H, slong prec);

/* Certified eig of a GENERAL n x n matrix A: heuristic acb_mat_approx_eig_qr to
 * seed, then acb_mat_eig_simple to PROVE isolating enclosures. Writes n
 * certified eigenvalue enclosures to E (caller _acb_vec_init'd, length n) and
 * right-eigenvector columns to R (caller acb_mat_init'd n x n). Aborts if the
 * spectrum cannot be isolated at this prec (a simple-spectrum requirement). */
void aic_mat_int_eig_certified(acb_ptr E, acb_mat_t R, const acb_mat_t A,
                               slong prec);

#endif /* AIC_MAT_INTERNAL_H */
