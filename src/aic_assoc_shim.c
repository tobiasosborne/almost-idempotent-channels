/* aic_assoc_shim.c — flat-double ccall SHIM for the associated eps-C* algebra
 * SUMMARY (bead aic-exa.3 [C], shim C2). See include/aic/aic_assoc_shim.h for the
 * ABI/layout contract and WHY the shim exists (no FLINT type may cross the Julia
 * ccall boundary). This file is PURE MARSHALLING: flat Kraus in, the EXACT cores
 * the C tests exercise (aic_assoc_ecstar_from_phi + aic_ecstar_defect_assoc), the
 * algebra summary + basis snapshot out. No new math. Close MIRROR of
 * aic_opspace_shim.c (same kraus_from_flat helper, same FLOOR/CEIL outward rounding
 * as aic_ucp_shim.c). Design: docs/research/julia_package_design.md §4.3 (C2).
 *
 * PROVENANCE (Law 1). The pipeline is verbatim the test_assoc2 / test_opspace
 * A-construction:
 *   aic_assoc_ecstar_from_phi(&ae, phi, prec)   (aic_assoc.h:177 — A = Img Phi_tilde,
 *                                                ASSERTS the prop_P basin, Rule 4)
 *   aic_ecstar_defect_assoc(eps, &ae.A, prec)   (aic_ecstar.h:136 — the EXACT-zero
 *                                                associator defect, the eps of A)
 * The basis {B_k} = ae.A.B[k] (aic_ecstar.h:92, Frobenius-ON, each n x n) is copied
 * out K-major. A = Img Phi_tilde is OBLIQUE (FINDINGS §C3) and its product is the
 * Choi-Effros star (FINDINGS §C2) — the basis is a read-once snapshot only.
 *
 * Number path: arb (acb_mat) internally; only the boundary is double. The eps
 * bracket is rounded OUTWARD (FLOOR/CEIL) so it stays RIGOROUS after the cast —
 * identical to aic_ucp_shim.c's aic_cbnorm_eigfree_d. Fail-loud (Rule 4): bad shape
 * aborts at the call site; the prop_P basin assert aborts deeper in
 * aic_assoc_ecstar_from_phi (the Julia caller pre-checks). <=200 LOC (Rule 10).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_assoc.h"
#include "aic/aic_assoc_shim.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_ucp.h"

/* Build the internal aic_ucp_kraus self-map Phi (dim_K == dim_H == n) from the
 * flat double Kraus arrays, layout kraus_*[a*n*n + i*n + j] (Convention A,
 * K-major). VERBATIM from aic_opspace_shim.c:47-58. */
static void kraus_from_flat(aic_ucp_kraus *phi, const double *kraus_re,
                            const double *kraus_im, int n, int r)
{
    aic_ucp_kraus_init(phi, n, n, r);
    for (int a = 0; a < r; a++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                int off = a * n * n + i * n + j;
                acb_set_d_d(acb_mat_entry(phi->K[a], i, j),
                            kraus_re[off], kraus_im[off]);
            }
}

/* aic_assoc_shim.h — build A = Img Phi_tilde, report (n, dim_A, eps_assoc[lo,hi]
 * FLOOR/CEIL, basis K-major). */
int aic_assoc_summary_d(int *out_n, int *out_dimA,
                        double *eps_assoc,
                        double *basis_re, double *basis_im,
                        const double *kraus_re, const double *kraus_im,
                        int n, int r, double eps, int prec)
{
    (void) eps; /* A's associator defect is derived from A directly (ABI uniform) */
    assert(n >= 1 && "aic_assoc_summary_d: n >= 1 required (Rule 4)");
    assert(r >= 1 && "aic_assoc_summary_d: r >= 1 required (Rule 4)");
    assert(prec >= 2 && "aic_assoc_summary_d: prec >= 2 required (Rule 4)");
    assert(out_n && out_dimA && eps_assoc &&
           "aic_assoc_summary_d: null scalar out (Rule 4)");
    assert(basis_re && basis_im &&
           "aic_assoc_summary_d: null basis out (Rule 4)");
    assert(kraus_re && kraus_im &&
           "aic_assoc_summary_d: null kraus array (Rule 4)");

    /* ---- (1) Phi -> A = Img Phi_tilde (asserts the prop_P basin) -------------- */
    aic_ucp_kraus phi;
    kraus_from_flat(&phi, kraus_re, kraus_im, n, r);

    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, (slong) prec);

    slong nn = ae.A.n, dimA = ae.A.dim_A;
    assert(nn == (slong) n && "aic_assoc_summary_d: A ambient n drift");
    *out_n    = (int) nn;
    *out_dimA = (int) dimA;

    /* ---- (2) eps = associator defect, FLOOR/CEIL outward (RIGOROUS) ----------- */
    arb_t e;
    arb_init(e);
    aic_ecstar_defect_assoc(e, &ae.A, (slong) prec);
    arf_t t;
    arf_init(t);
    arb_get_lbound_arf(t, e, (slong) prec);
    eps_assoc[0] = arf_get_d(t, ARF_RND_FLOOR);
    arb_get_ubound_arf(t, e, (slong) prec);
    eps_assoc[1] = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
    arb_clear(e);

    /* ---- (3) copy {B_k} out, K-major [k*n*n + i*n + j] ------------------------ */
    for (slong k = 0; k < dimA; k++)
        for (slong i = 0; i < nn; i++)
            for (slong j = 0; j < nn; j++) {
                slong off = k * nn * nn + i * nn + j;
                acb_srcptr ent = acb_mat_entry(ae.A.B[k], i, j);
                basis_re[off] = arf_get_d(arb_midref(acb_realref(ent)), ARF_RND_NEAR);
                basis_im[off] = arf_get_d(arb_midref(acb_imagref(ent)), ARF_RND_NEAR);
            }

    /* ---- (4) free (ae owns the superop S_tilde + ctx + basis) ----------------- */
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return 0;
}
