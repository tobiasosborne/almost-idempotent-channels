/* aic_idemp_delta_internal.h — private helpers for the prop_Delta Kraus-form
 * extraction (bead aic-ynu, I4). Not part of the public API
 * (include/aic/aic_idemp.h). The driver + K_{j,c} build live in
 * aic_idemp_delta_kraus.c; the A_hat / Sigma read-offs and the unitality + match
 * certifications in aic_idemp_delta_kraus_ops.c. See aic_idemp_delta_kraus.c for
 * the narrative (.tex:2098-2103). Vec convention: vec(X)[i*c+j]=X[i,j] (row-major).
 */
#ifndef AIC_IDEMP_DELTA_INTERNAL_H
#define AIC_IDEMP_DELTA_INTERNAL_H

#include <flint/acb_mat.h>

#include "aic/aic_idemp.h"

/* Read A_hat = (+)_j A_j (the abstract A-element as an operator on L = (+)_j L_j,
 * a dL_tot x dL_tot block-diagonal matrix) from the B(M)-element wA = w(A) in B(M):
 * the block component A_j = Tr_{E_j}(W_j wA W_j^dag) / dim_E[j] is placed at the
 * block offset off_L[j]. Gauge-safe (works through the *-monomorphism w). */
void aic_dk_Ahat_from_wA(acb_mat_t Ahat, const aic_idemp_delta *out,
                         const aic_idemp_wedderburn *W, const acb_mat_t wA,
                         slong prec);

/* Apply the M-part Kraus map to a block-diagonal A_hat (dL_tot x dL_tot):
 * out_M = sum_{j,c} K_{j,c} A_hat K_{j,c}^dag, an n x n operator on H (supported on
 * M). `Mpart` must be init'd n x n. */
void aic_dk_apply_Mpart(acb_mat_t Mpart, const aic_idemp_delta *out,
                        const acb_mat_t Ahat, slong prec);

/* Read the M^perp Sigma matrix off the stored d->Delta: Sigma column k =
 * vec(J_Mperp^dag B_k J_Mperp), B_k = column k of d->Delta reshaped n x n. Stored
 * (dim_Mperp^2) x dim_A (row-major vec). No-op (Sigma stays 0x0) when dim_Mperp=0
 * (M = H). Allocates out->Sigma. */
void aic_dk_build_sigma(aic_idemp_delta *out, const aic_idemp_decomp *d, slong prec);

/* Embed the M^perp Sigma component of A-coordinate vector a (dim_A) back into B(H):
 * Mperp_part = J_Mperp Sigma(a) J_Mperp^dag, where Sigma(a) = reshape(Sigma * a).
 * `Mperp_part` must be init'd n x n; zeroed when dim_Mperp = 0. */
void aic_dk_apply_Sigma(acb_mat_t Mperp_part, const aic_idemp_delta *out,
                        const aic_idemp_decomp *d, const acb_mat_t a, slong prec);

/* Certify (arb, fail-loud) the M-part unitality Delta_M(1_A) = Pi_M: feed
 * A_hat = 1_{dL_tot} (= the abstract unit 1_A) to the Kraus M-part and compare to
 * d->Pi_M. .tex:2100 (Delta(1_A)|_M = w(1_A) = Pi_M, lem_idemp). */
void aic_dk_certify_unital(const aic_idemp_delta *out, const aic_idemp_decomp *d,
                           slong prec);

/* THE CRUX TOOTH (.tex:2100). Certify (arb, fail-loud) that the full Delta-Kraus
 * map equals the stored d->Delta on the A-basis B_k: for each k, read A_hat from
 * wA = w(B_k) (= reshape of d->w column k), reconstruct Delta_kraus(B_k) = M-part
 * (Kraus K_{j,c}) + M^perp-part (Sigma), and compare to d->Delta column k (reshaped
 * n x n), ||Delta_kraus(B_k) - B_k||_op <= tol. Pins the K_{j,c} build + W_j usage
 * + Sigma read-off at once. */
void aic_dk_certify_match(const aic_idemp_delta *out, const aic_idemp_wedderburn *W,
                          const aic_idemp_decomp *d, slong prec);

#endif /* AIC_IDEMP_DELTA_INTERNAL_H */
