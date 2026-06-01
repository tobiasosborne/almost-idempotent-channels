/* aic_mat_internal.h — shared internals for the mat module (beads aic-9zh,
 * "T-mat"). NOT part of the public API (include/aic_mat.h); these helpers are
 * shared between aic_mat_spectral.c (degeneracy-robust max-eig) and
 * aic_mat_eig.c (full Hermitian eigendecomposition) so the Hermiticity
 * precondition and the certified two-stage eig path are written once.
 */
#ifndef AIC_MAT_INTERNAL_H
#define AIC_MAT_INTERNAL_H

#include <complex.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"      /* aic_mat_eigcluster (subspace helpers below) */

/* prec-appropriate absolute tol = 2^-(prec - 8); the 8-bit guard leaves headroom
 * so a faithfully-rounded Hermitian input (radii ~2^-prec) is not rejected. */
void aic_mat_int_tol(arb_t tol, slong prec);

/* Non-aborting Hermiticity predicate: 1 iff |H[i,j] - conj(H[j,i])| is certainly
 * within the RELATIVE + absolute-floor tol tol*(1 + |H_ij| + |H_ji|) for all i,j
 * (bead aic-2yo: a bare ABSOLUTE tol false-fails on graded/ill-conditioned Gram
 * matrices whose entry radii grow with magnitude). Asserts squareness. */
int aic_mat_int_is_hermitian(const acb_mat_t H, slong prec);

/* Abort (fail loud, CLAUDE.md Rule 4) unless aic_mat_int_is_hermitian(H, prec). */
void aic_mat_int_assert_hermitian(const acb_mat_t H, slong prec);

/* Certified eig of a GENERAL n x n matrix A: heuristic acb_mat_approx_eig_qr to
 * seed, then acb_mat_eig_simple to PROVE isolating enclosures. Writes n
 * certified eigenvalue enclosures to E (caller _acb_vec_init'd, length n) and
 * right-eigenvector columns to R (caller acb_mat_init'd n x n). Aborts if the
 * spectrum cannot be isolated at this prec (a simple-spectrum requirement). */
void aic_mat_int_eig_certified(acb_ptr E, acb_mat_t R, const acb_mat_t A,
                               slong prec);

/* The dense-unitary DENSIFIER (docs/research/eigvec_certified_design.md §1.3).
 * Builds U = product over ALL planes (a,b), a<b, of the rational Givens rotation
 * cos = 3/5, sin = 4/5 (the exact rationals already used in test_eigmult.c's
 * fixtures). U is unitary far below the working precision (||U U† - I||_F
 * certified ~1.3e-37 at n=4, prec=128 — measured). Conjugating A' = U A U†
 * SPREADS every invariant subspace across all n coordinates, defeating FLINT
 * Rump's frozen-row partition failure on ROW-SPARSE subspaces (FINDINGS §D5n;
 * design §2). The spectrum is conjugation-invariant, so eigenvalue balls of A'
 * equal those of A. `U` must be init'd n x n by the caller; the caller forms
 * A' = U A U† (two acb_mat_mul with U† = acb_mat_conjugate_transpose). Shared
 * between aic_mat_eig_multiple.c (eigenvalue retry) and aic_mat_eigvec.c
 * (subspace densify, increment-2 step C2). */
void aic_mat_dense_unitary(acb_mat_t U, slong n, slong prec);

/* Subspace-layer fail-loud helpers shared between aic_mat_eigvec.c (the
 * orchestration) and aic_mat_eigvec_seed.c (the primitives), design
 * docs/research/eigvec_certified_design.md §1.2/§1.6. NOT public.
 *
 * aic_mat_int_assert_densify_unitary: ABORT (Rule 4) unless ||U U^dag - I||_F is
 * certified < 2^-(prec-8); guards the densify conjugation A' = U H U^dag (a loose
 * U silently moves the spectrum). */
void aic_mat_int_assert_densify_unitary(const acb_mat_t U, const acb_mat_t Ud,
                                        slong n, slong prec);

/* aic_mat_int_certify_cluster: certify one cluster [s0, s0+k) on the densified A'.
 * Builds the Rump seed Xa from the zheev eigenvector columns Vd (n x n row-major
 * double, column index = eigenvalue index) for the cluster, lambda_approx =
 * mean(ev[s0..s0+k)), runs acb_mat_eig_enclosure_rump on A1, ASSERTS a FINITE +
 * REAL enclosure (else fail loud "UNRESOLVED"/"non-real"), and back-maps
 * X_c = Ud X'_c. Writes lambda/X/J/k into *out (all init'd here; free via
 * aic_mat_eigcluster_free). cidx is the cluster index (for the abort message). */
void aic_mat_int_certify_cluster(aic_mat_eigcluster *out, const acb_mat_t A1,
                                 const acb_mat_t Ud, const double _Complex *Vd,
                                 const double *ev, slong n, slong s0, slong k,
                                 slong cidx, slong prec);

#endif /* AIC_MAT_INTERNAL_H */
