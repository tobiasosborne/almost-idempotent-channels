/* aic_ucp_kraus_arb_internal.h — helpers shared between the certified Choi->Kraus
 * orchestration (aic_ucp_kraus_arb.c) and its threshold + orthonormalisation
 * primitives (aic_ucp_kraus_arb_orth.c). Split out for Rule 10 (~200 LOC/file).
 * NOT public (no entry in include/aic/aic_ucp.h). Bead aic-4td increment 2 step D.
 */
#ifndef AIC_UCP_KRAUS_ARB_INTERNAL_H
#define AIC_UCP_KRAUS_ARB_INTERNAL_H

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

/* keep_thr = (dim_K*dim_H) * 2^-52 * ||C||_F as a certified arb ball (mirrors the
 * _latd DBL_EPSILON=2^-52 QETLAB gate; ||C||_F is the certified Frobenius norm).
 * n = dim_K*dim_H. `thr` must be an initialised arb_t. */
void aic_ucp_kraus_arb_int_keep_thr(arb_t thr, const acb_mat_t C, slong n,
                                    slong prec);

/* Orthonormalise a cluster's NON-orthonormal invariant-subspace enclosure X
 * (n x k) into V (n x k, caller-init'd) by LOEWDIN symmetric orthonormalisation
 * V = X (X^dag X)^{-1/2}: V^dag V = I_k, range(V) = range(X). k=1 normalises the
 * single certified eigenvector. See aic_ucp_kraus_arb.c's docstring (R2) for why
 * a raw reshape of X would not rebuild C, and why x0 = ||X||_op^2 (NOT
 * herm_max_eig of the raw interval Gram). */
void aic_ucp_kraus_arb_int_loewdin(acb_mat_t V, const acb_mat_t X, slong prec);

/* FINDING-2 FIX: per-EIGENVALUE orthonormal eigenbasis of a kept cluster (resolves a
 * lumped-distinct cluster, design §3.2 Option B). For an outer cluster's NON-orthonormal
 * X (n x k) and the original Choi C, diagonalises the small dense compression
 * M = V^dag C V (V = Loewdin(X)) and assembles the orthonormal eigenbasis VW_out (n x k,
 * caller-init'd, sub-eigenvalue order) with per-column sqrt(eigenvalue) scales (caller
 * arb[k]). EXACT even when the cluster lumps distinct eigenvalues; reduces to the single
 * sqrt(lambda_c) for a genuine degeneracy. k>=1 (k=1 returns the single normalised
 * eigenvector with sqrt(lambda)). See aic_ucp_kraus_arb_orth.c for the full rationale. */
void aic_ucp_kraus_arb_int_cluster_basis(acb_mat_t VW_out, arb_ptr scales,
                                         const acb_mat_t X, const acb_mat_t C,
                                         slong prec);

#endif /* AIC_UCP_KRAUS_ARB_INTERNAL_H */
