/* aic_idemp_maps.c — the maps Lambda, Gamma, w of th_idemp_structure and the
 * top-level orchestrator aic_idemp_decompose (bead aic-wuh). Steps 3-7 of the
 * construction (approximate_algebras.tex:2056-2090). The proof IS the algorithm.
 *
 *   Delta : A -> B(H)  = inclusion  = the (n^2 x dim_A) image basis [vec B_k]
 *           (.tex:2056; Delta^dag Delta = 1_A since the columns are orthonormal).
 *   C_M   : X |-> J_M^dag X J_M : B(H) -> B(M)   (compression, .tex:2059).
 *   Lambda: Y |-> Phi(J_M Y J_M^dag) : B(M) -> B(H)   (.tex:2061), CP + unital.
 *   Gamma : B(M) -> A  via  Lambda = Delta Gamma   (.tex:2065). Delta^dag Delta =
 *           1_A => Gamma_mat = Delta^dag Lambda_mat (NO pseudoinverse; the paper
 *           guarantees this Gamma is UCP, .tex:2088).
 *   w     : A -> B(M),  w = C_M Delta  (.tex:2084),  column k = vec(J_M^dag B_k J_M).
 *
 * Maps are stored as their MATRIX in the row-major vec convention vec(X)[i*n+j]=
 * X[i,j] (see include/aic_idemp.h): a map B(C^p)->B(C^q) is a (q*q)x(p*p) matrix.
 *
 * Exact-idempotence guard (Rule 4). This module is ONLY for exactly-idempotent
 * Phi (the eta=0 oracle; eta>0 is assoc_ecsa's job via Phi_tilde). At ENTRY we
 * compute max_ij ||Phi(Phi(E_ij)) - Phi(E_ij)||_op and ABORT if it exceeds a
 * prec-appropriate tol — a non-idempotent Phi has no carrier/image structure
 * theorem and the downstream relations would be silently wrong.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic_idemp_internal.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

/* max_ij ||Phi(Phi(E_ij)) - Phi(E_ij)||_op as a certified ball (idempotence). */
static void idemp_defect(arb_t out, const aic_ucp_kraus *phi, slong n, slong prec)
{
    acb_mat_t E, PhiE, PhiPhiE, diff;
    arb_t d;
    acb_mat_init(E, n, n);
    acb_mat_init(PhiE, n, n);
    acb_mat_init(PhiPhiE, n, n);
    acb_mat_init(diff, n, n);
    arb_init(d);
    arb_zero(out);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            acb_mat_zero(E);
            acb_set_si(acb_mat_entry(E, i, j), 1);
            aic_ucp_apply(PhiE, phi, E, prec);
            aic_ucp_apply(PhiPhiE, phi, PhiE, prec);
            acb_mat_sub(diff, PhiPhiE, PhiE, prec);
            aic_mat_opnorm(d, diff, prec);
            if (arb_gt(d, out)) arb_set(out, d);
        }
    arb_clear(d);
    acb_mat_clear(diff);
    acb_mat_clear(PhiPhiE);
    acb_mat_clear(PhiE);
    acb_mat_clear(E);
}

/* reshape column k of Delta (length n^2, row-major vec) into the n x n matrix B. */
static void delta_col_to_mat(acb_mat_t B, const acb_mat_t Delta, slong k, slong n)
{
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(B, a, b), acb_mat_entry(Delta, a * n + b, k));
}

void aic_idemp_decompose(aic_idemp_decomp *out, const aic_ucp_kraus *phi,
                         slong prec)
{
    slong n = phi->dim_H;
    assert(phi->dim_K == n && "aic_idemp_decompose: Phi must be a self-map B(H)->B(H)");

    /* --- entry guard: Phi must be EXACTLY idempotent (Rule 4) --- */
    {
        arb_t def, tol;
        arb_init(def);
        arb_init(tol);
        idemp_defect(def, phi, n, prec);
        /* This guard measures max_ij ||Phi^2(E_ij) - Phi(E_ij)||_op: a max over the
         * matrix units E_ij of the OP-NORM defect, which is a LOWER BOUND on the
         * paper's cb-norm eta = ||Phi^2-Phi||_cb (cb-norm != op-norm, CLAUDE.md
         * callout: the cb-norm is a sup over ampliations 1_n (x) Phi). For the
         * exact (eta=0) oracle both vanish, so this lower bound is SUFFICIENT to
         * gate this module: an exactly-idempotent Phi passes and a grossly
         * non-idempotent one (op-defect >> 1e-9) is rejected. It does NOT certify
         * almost-idempotence at small eta — a ~1e-5 Kraus perturbation can have
         * op-defect ~1e-10 < 1e-9 and would be wrongly accepted here; the proper
         * cb-norm (diamond-norm) gate is bead aic-d24, out of scope for this oracle.
         * prec=53 resolves an O(1) difference to ~1e-15; 1e-9 is generous slack. */
        arb_set_d(tol, 1e-9);
        if (!arb_le(def, tol)) {
            flint_printf("aic_idemp_decompose: Phi is NOT exactly idempotent "
                         "(max_ij ||Phi^2(E_ij)-Phi(E_ij)||_op exceeds 1e-9); this "
                         "module is the eta=0 oracle only (Rule 4: fail loud)\n");
            flint_abort();
        }
        arb_clear(tol);
        arb_clear(def);
    }

    out->n = n;

    /* --- step 1: carrier M (.tex:2056) --- */
    acb_mat_t Q;
    acb_mat_init(Q, n, n);
    aic_ucp_carrier_Q(Q, phi, prec);
    out->dim_M = aic_idemp_carrier_split(out->J_M, out->J_Mperp, out->Pi_M, Q, n, prec);
    acb_mat_clear(Q);

    /* --- step 2: A = Img Phi, ONB columns of Delta (.tex:2056) --- */
    out->dim_A = aic_idemp_image_basis(out->Delta, phi, n, prec);
    slong dA = out->dim_A;

    /* verify Phi(B_k) = B_k (fixed points): part of step-2 contract. */
    {
        acb_mat_t B, PhiB, diff;
        arb_t fdef, tol;
        acb_mat_init(B, n, n);
        acb_mat_init(PhiB, n, n);
        acb_mat_init(diff, n, n);
        arb_init(fdef);
        arb_init(tol);
        arb_set_d(tol, 1e-9);
        for (slong k = 0; k < dA; k++) {
            delta_col_to_mat(B, out->Delta, k, n);
            aic_ucp_apply(PhiB, phi, B, prec);
            acb_mat_sub(diff, PhiB, B, prec);
            aic_mat_opnorm(fdef, diff, prec);
            if (!arb_le(fdef, tol)) {
                flint_printf("aic_idemp_decompose: basis element B_%wd is not a "
                             "fixed point of Phi (Rule 4)\n", k);
                flint_abort();
            }
        }
        arb_clear(tol);
        arb_clear(fdef);
        acb_mat_clear(diff);
        acb_mat_clear(PhiB);
        acb_mat_clear(B);
    }

    /* --- steps 5-7: Lambda, Gamma, w (.tex:2061-2084; aic_idemp_build.c) --- */
    aic_idemp_build_maps(out, phi, prec);
}

void aic_idemp_clear(aic_idemp_decomp *out)
{
    /* Requires a FULLY-decomposed struct: clears all 7 acb_mat fields
     * unconditionally. Safe today because aic_idemp_decompose either fills every
     * field or flint_abort()s (terminating the process), so a partially-built
     * struct is never handed to clear. */
    if (out == NULL) return;
    acb_mat_clear(out->J_M);
    acb_mat_clear(out->J_Mperp);
    acb_mat_clear(out->Pi_M);
    acb_mat_clear(out->Delta);
    acb_mat_clear(out->Lambda);
    acb_mat_clear(out->Gamma);
    acb_mat_clear(out->w);
}
