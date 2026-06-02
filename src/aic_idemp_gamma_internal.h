/* aic_idemp_gamma_internal.h — private helpers for the prop_Gamma Kraus-form
 * extraction (bead aic-ynu, I3). Not part of the public API
 * (include/aic/aic_idemp.h). The driver + gamma_j least-squares extraction live
 * in aic_idemp_gamma_kraus.c; the small linear-algebra + certification helpers in
 * aic_idemp_gamma_kraus_ops.c. See aic_idemp_gamma_kraus.c for the narrative
 * (.tex:2106-2122). Vec convention: vec(X)[i*m+j]=X[i,j] (row-major).
 */
#ifndef AIC_IDEMP_GAMMA_INTERNAL_H
#define AIC_IDEMP_GAMMA_INTERNAL_H

#include <complex.h>

#include <flint/acb_mat.h>

#include "aic/aic_idemp.h"

/* midpoint of an acb entry as a C double _Complex (double-path solve scratch). */
double _Complex aic_gk_d(acb_srcptr z);

/* wG(X) = w( d->Gamma(X) ) in B(M): apply the stored conditional expectation
 * d->Gamma to X in B(M) (giving an A coordinate vector), then embed in B(M) via
 * the *-monomorphism w. out is m x m (m = dim_M). Gauge-safe target for the
 * per-block Gamma_j read-off and the match certification. */
void aic_gk_wG(acb_mat_t out, const aic_idemp_decomp *d, const acb_mat_t X,
               slong prec);

/* Solve the (overdetermined) complex least-squares A g = b (A is nR x nU,
 * row-major; b length nR) via the normal equations A^dag A g = A^dag b and a tiny
 * Gauss-Jordan solve. Writes the nU-vector g (row-major) and the scaled residual
 * ||A g - b|| / (||b|| + 1) into *resid (the consistency check: ~0 iff b is in the
 * column space, i.e. d->Gamma IS of prop_Gamma form for this block). Double path. */
void aic_gk_lstsq(double _Complex *g, double *resid, const double _Complex *A,
                  const double _Complex *b, slong nR, slong nU);

/* In-place Hermitian symmetrization gamma <- (gamma + gamma^dag)/2. Applied AFTER
 * the density certificate (which proves the asymmetry is below tol) so the
 * downstream tight-Hermiticity routines (eig subspaces) accept gamma. */
void aic_gk_hermitize(acb_mat_t gamma, slong prec);

/* Certify (arb, fail-loud) that the RAW extracted gamma is Hermitian:
 * ||gamma - gamma^dag||_op <= tol (caught BEFORE aic_gk_hermitize symmetrizes).
 * `j` is the block index for the error message. */
void aic_gk_certify_hermitian(const acb_mat_t gamma, slong j, slong prec);

/* Certify (arb, fail-loud) that gamma (dE x dE, already Hermitized) is a valid
 * density matrix: PSD (min eigenvalue >= -tol) and trace 1 (|Tr - 1| ~ 0).
 * `j` is the block index for the error message. */
void aic_gk_certify_density(const acb_mat_t gamma, slong j, slong prec);

/* Build the Kraus isometry L_j = (W_j^dag (x) 1_{F_j})(1_{L_j} (x) xi_j), F_j=E_j,
 * with xi_j = sum_c sqrt(p_c) |v_c>_E (x) |c>_F the purification of gamma_j =
 * sum_c p_c |v_c><v_c| (Hermitian eig). L is (m*dE) x dL. Tensor layout is
 * LEFT-major throughout (aic_mat_kronecker): the M index is major over F_j on the
 * output, the L_j index major over E_j inside W_j. Certifies L_j^dag L_j = 1_{L_j}
 * (arb, fail-loud). `j` is the block index for the error message. */
void aic_gk_build_Lj(acb_mat_t L, const acb_mat_t gamma, const acb_mat_t Wj,
                     slong dL, slong dE, slong m, slong j, slong prec);

/* THE CRUX TOOTH (.tex:2113). Certify (arb, fail-loud) that the prop_Gamma Kraus
 * map equals the stored d->Gamma on a basis E_ab of B(M): both A-elements are
 * embedded in B(M) via w (gauge-invariant) and compared,
 *   ||w(Gamma_kraus(E_ab)) - w(d->Gamma(E_ab))||_op <= tol  for all a,b.
 * Gamma_kraus(E_ab) re-embedded as sum_j W_j^dag (Gamma_j(E_ab) (x) 1_{E_j}) W_j,
 * with Gamma_j(E_ab) computed BOTH ways (eq Gamma Tr_{E_j} form AND the Choi L_j
 * form) — both must match d->Gamma, pinning the gamma_j AND the L_j build. */
void aic_gk_certify_match(const aic_idemp_gamma *out, const aic_idemp_wedderburn *W,
                          const aic_idemp_decomp *d, slong prec);

#endif /* AIC_IDEMP_GAMMA_INTERNAL_H */
