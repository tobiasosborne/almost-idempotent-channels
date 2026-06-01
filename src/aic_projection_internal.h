/* aic_projection_internal.h — shared internals for the projection module (bead
 * aic-mqf), split across src/aic_projection.c (orchestrator + Steps 3-5) and
 * src/aic_projection_find.c (Steps 1-2: H-search + gap). NOT a public API
 * (include/aic_projection.h). Split for the ~200 LOC/file limit (Rule 10).
 */
#ifndef AIC_PROJECTION_INTERNAL_H
#define AIC_PROJECTION_INTERNAL_H

#include <flint/acb_mat.h>

#include "aic/aic_ecstar.h"
#include "aic/aic_projection.h"   /* aic_projection_witness */

/* H_k = (B_k + B_k^dag)/2 into `out` (caller-init'd n x n). */
void aic_projection_herm_part(acb_mat_t out, const acb_mat_t Bk, slong prec);

/* Largest interior spectral gap of Hermitian H (double-path zheev eigenvalues):
 * threshold t -> t_out, gap -> g_out, lam_min and lam_max -> lmin, lmax, the
 * +1-cluster size m -> m_out. FAILS LOUD if no interior gap >= g_min. */
void aic_projection_gap(double *t_out, double *g_out, double *lmin, double *lmax,
                        slong *m_out, const acb_mat_t H, slong prec);

/* Step 3: ambient spectral projector P_amb = (I + sgn(s(H-tI)))/2 (exact ambient
 * idempotent), s = 1/max(t-lmin, lmax-t). ASSERTS the funcalc basin ||Y^2-I||<1
 * (tex:516,520). `Pamb` caller-init'd n x n. */
void aic_projection_ambient(acb_mat_t Pamb, const acb_mat_t H, double t,
                            double lmin, double lmax, slong prec);

/* Step 4: P = Phi_tilde(W) = W * I (the star with the unit; FINDINGS §C3 — A is an
 * oblique image, so Phi_tilde is the structure-respecting projector). `out` n x n. */
void aic_projection_into_A(acb_mat_t out, const aic_ecstar *A, const acb_mat_t W,
                           slong prec);

/* ONE interior-gap candidate of a Hermitian H_k (a (H_k, gap) pair the audition
 * may build P from). Doubles are the DOUBLE-path zheev spectrum (prec-independent
 * gap-selection); the defect is certified downstream in the arb path. */
typedef struct {
    slong k;        /* basis index k whose H_k = (B_k+B_k^dag)/2 supplies this gap */
    double gap;     /* g = lam_{i+1} - lam_i (the interior gap)                    */
    double t;       /* threshold = (lam_i + lam_{i+1})/2                           */
    double lmin;    /* smallest eigenvalue of H_k                                  */
    double lmax;    /* largest eigenvalue of H_k                                   */
    slong m;        /* # eigenvalues >= t (the +1-cluster size)                    */
} aic_projection_cand;

/* Enumerate ALL interior-gap candidates over ALL non-scalar H_k in A's basis,
 * sorted by `gap` DESCENDING (the unit-aware audition tries the largest first).
 * Allocates *cands (caller frees via free()) and writes the count to *n_cands.
 * FAILS LOUD (the aic-3qv stop condition) if NO non-scalar H_k yields a positive
 * interior gap (a genuinely degenerate algebra). */
void aic_projection_enum_cands(aic_projection_cand **cands, slong *n_cands,
                               const aic_ecstar *A, slong prec);

/* The UNIT-AWARE gap audition (the aic-66n fix). Computes U_A = Phi_tilde(1_n)
 * once, then auditions the (H_k, gap) candidates largest-gap-first: for each it
 * builds P_amb (ambient spectral projector of H_k at the gap) and P = Phi_tilde
 * (P_amb), ACCEPTING the first with ||P||_op > c AND ||U_A - P||_op > c. Fills the
 * accepted P (caller-init'd n x n) and *wit (the chosen H/gap). FAILS LOUD (the
 * aic-3qv stop condition) if NO candidate yields a nontrivial-vs-unit split.
 * `c` is the nontriviality threshold (>= the cstar_build wrapper backstop 0.15). */
void aic_projection_audit(acb_mat_t P, const aic_ecstar *A, const acb_mat_t U_A,
                          aic_projection_witness *wit, double c, slong prec);

#endif /* AIC_PROJECTION_INTERNAL_H */
