/* aic_factorize_artifacts_internal.h — PRIVATE shared build helpers for the
 * factorize-artifacts shims C4 (aic_factorize_artifacts_shim.c) and C5
 * (aic_factorize_artifacts2_shim.c), split per Rule 10 (each <=200 LOC). NOT a
 * public ABI header (no aic/ install path): it only shares the SAME pipeline build
 * aic_factorize_choi_shim_d uses so C4 and C5 stay byte-identical and do NOT
 * re-derive the pipeline (design §4.2). Internal symbols (no aic_ prefix collision
 * with the public API): aic_fart_build_all/_free are DEFINED ONCE in the C4 TU
 * (src/aic_factorize_artifacts_shim.c — nm shows `T`) and REFERENCED by the C5 TU
 * (src/aic_factorize_artifacts2_shim.c — nm shows `U`); the static lib resolves the
 * reference at link. No ODR issue (single definition).
 */
#ifndef AIC_FACTORIZE_ARTIFACTS_INTERNAL_H
#define AIC_FACTORIZE_ARTIFACTS_INTERNAL_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_dhom.h"
#include "aic/aic_factorize.h"
#include "aic/aic_ucp.h"

/* Build the full factorize handle from flat Kraus (the SHARED C4/C5 pipeline,
 * verbatim aic_factorize_shim.c:130-167). Fills *F (delta_ready + upsilon_ready),
 * *phi, *ae, *B, *v, *iso. Free in reverse order via aic_fart_build_free.
 * out_eta[0]=eta_proxy, out_eta[1]=eps_used (the eps SENTINEL resolved). */
void aic_fart_build_all(aic_factorize *F, aic_ucp_kraus *phi, aic_assoc_ecstar *ae,
                        aic_dhom_B *B, aic_dhom_v *v, arb_t iso, double *out_eta,
                        const double *kraus_re, const double *kraus_im,
                        int n, int r, double eps, slong prec);

void aic_fart_build_free(aic_factorize *F, aic_ucp_kraus *phi, aic_assoc_ecstar *ae,
                         aic_dhom_B *B, aic_dhom_v *v, arb_t iso);

#endif /* AIC_FACTORIZE_ARTIFACTS_INTERNAL_H */
