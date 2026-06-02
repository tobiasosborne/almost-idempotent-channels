/* aic_idemp_delta_kraus.c — the explicit prop_Delta Kraus form of the inclusion
 * map Delta : A = (+)_j B(L_j) -> B(H) (bead aic-ynu, I4). Realizes
 * Delta_structure (approximate_algebras.tex:2098-2103) for the EXACTLY-idempotent
 * (eta=0) case, the th_idemp_structure oracle.
 *
 * approximate_algebras.tex:2098-2103 — eq Delta_structure (verbatim core):
 *   "Delta(A) = sum_{j=1}^m J_M W_j^dag (A_j (x) 1_{E_j}) W_j J_M^dag
 *             + J_{M^perp} Sigma(A) J_{M^perp}^dag      (A in A),"
 * where Delta(A)|_M = w(A) (eq Aw, :2095) and Sigma : A -> B(M^perp) is UCP.
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (Law 3). The paper merely ASSERTS the
 * existence of the UCP Sigma realizing the M^perp-part (lem_idemp block-diagonality
 * gives the M-part as the *-hom w; the M^perp-part is "obtained from A by some UCP
 * map Sigma"). In finite dim BOTH parts are computable:
 *   - the M-part Kraus operators K_{j,c} are built DIRECTLY from W_j + J_M (below);
 *   - Sigma is READ OFF the stored, specific inclusion d->Delta (the
 *     th_idemp_structure Delta = [vec B_k] columns): Sigma(B_k) = J_Mperp^dag B_k
 *     J_Mperp. We do not need to re-derive Sigma; we match the stored Delta exactly.
 *
 * THE M-PART KRAUS FORM. Expanding A_j (x) 1_{E_j} = sum_{c in E_j}
 * (1_{L_j} (x) e_c) A_j (1_{L_j} (x) e_c)^dag and J_M w(A) J_M^dag gives
 *     Delta_M(A) = sum_{j,c} K_{j,c} A_hat K_{j,c}^dag,
 *     K_{j,c} = J_M W_j^dag (1_{L_j} (x) e_c) P_j : L=(+)_k L_k -> H,
 * an n x dL_tot operator (dL_tot = sum_j dim_L[j]); P_j : L -> L_j the block
 * projection, (1_{L_j} (x) e_c) : L_j -> L_j (x) E_j the c-th multiplicity
 * injection. A_hat = (+)_j A_j is A as an operator on its carrier L = (+)_j L_j.
 * Because A_hat is block-diagonal, sum_{j,c} K_{j,c} A_hat K_{j,c}^dag picks the
 * diagonal blocks A_j only, reproducing J_M w(A) J_M^dag exactly.
 *
 * READING A_hat FROM A-COORDINATES (gauge-safe, through w). The decomp stores an
 * A-element as A-coordinates a_k (a = (a_k), A = sum_k a_k B_k). Working through w
 * avoids the abstract A-basis convention: w(A) = sum_k a_k w(B_k) (w(B_k) = reshape
 * of d->w column k), and W_j w(A) W_j^dag = A_j (x) 1_{E_j}, so
 *     A_j = Tr_{E_j}(W_j w(A) W_j^dag) / dim_E[j]   (W_j isometry onto block j).
 * This mirrors the gamma_kraus aic_gk_wG idiom. The match certification feeds the
 * A-basis B_k (a = e_k) and compares Delta_kraus(B_k) to d->Delta column k.
 *
 * THE M^perp / Sigma PART. Sigma is the M^perp-block of the SPECIFIC inclusion
 * d->Delta (not free): Sigma(B_k) = J_Mperp^dag B_k J_Mperp (B_k = column k of
 * d->Delta reshaped). Stored as a (dmp^2) x dim_A matrix (row-major vec), empty
 * when dim_Mperp = 0 (M = H, the trace-preserving channels).
 *
 * THE GATE (.tex:2100, the crux tooth). Delta_kraus(B_k) = J_M (M-part) J_M^dag +
 * J_Mperp Sigma_k J_Mperp^dag MUST equal d->Delta column k (reshaped) to ~1e-10,
 * for every A-basis element B_k. Plus: Delta_M(1_A) = Pi_M (M-part unitality). Both
 * arb-certified in production, fail-loud (Rule 4). A genuine mismatch is a
 * W_j/convention bug, escalated, not faked (.tex:2098). See
 * aic_idemp_delta_kraus_ops.c for the helpers (A_hat read-off, the M^perp Sigma
 * read-off, the unitality + match certification).
 */
#include <assert.h>

#include <flint/flint.h>
#include <flint/acb_mat.h>

#include "aic/aic_idemp.h"
#include "aic_idemp_delta_internal.h"

/* Build the M-part Kraus operator K_{j,c} = J_M W_j^dag (1_{L_j} (x) e_c) P_j as an
 * n x dL_tot matrix. P_j embeds L_j at column offset off_L[j] in L = (+)_k L_k;
 * (1_{L_j} (x) e_c) maps L_j -> L_j (x) E_j with W_j-row index a*dim_E+c (the
 * L-major tensor layout of the Wedderburn W_j, include/aic_idemp.h). So
 *   K_{j,c}[:, off_L[j] + a] = (J_M W_j^dag)[:, a*dim_E + c],  a in [0,dim_L). */
static void build_Kjc(acb_mat_t K, const acb_mat_t JMWjd, slong dL, slong dE,
                      slong c, slong off, slong n)
{
    acb_mat_zero(K);
    for (slong a = 0; a < dL; a++)
        for (slong r = 0; r < n; r++)
            acb_set(acb_mat_entry(K, r, off + a),
                    acb_mat_entry(JMWjd, r, a * dE + c));
}

void aic_idemp_delta_kraus(aic_idemp_delta *out, const aic_idemp_wedderburn *W,
                           const aic_idemp_decomp *d, slong prec)
{
    slong n = d->n, m = W->dim_M, nb = W->num_blocks;
    assert(m == d->dim_M && W->dim_A == d->dim_A && "delta_kraus: W/d mismatch");

    out->n = n;
    out->dim_M = m;
    out->dim_Mperp = n - m;
    out->dim_A = W->dim_A;
    out->num_blocks = nb;
    out->dim_L = flint_malloc((size_t) nb * sizeof(slong));
    out->dim_E = flint_malloc((size_t) nb * sizeof(slong));
    out->off_L = flint_malloc((size_t) nb * sizeof(slong));

    /* block offsets in L = (+)_j L_j and the total carrier dim dL_tot. */
    slong dL_tot = 0, n_kraus = 0;
    for (slong j = 0; j < nb; j++) {
        out->dim_L[j] = W->dim_L[j];
        out->dim_E[j] = W->dim_E[j];
        out->off_L[j] = dL_tot;
        dL_tot += W->dim_L[j];
        n_kraus += W->dim_E[j];
    }
    out->dL_tot = dL_tot;
    out->n_kraus = n_kraus;
    out->K = flint_malloc((size_t) n_kraus * sizeof(acb_mat_t));

    /* (1) build the M-part Kraus operators K_{j,c} (one per block j, mult index c). */
    slong kk = 0;
    for (slong j = 0; j < nb; j++) {
        slong dL = W->dim_L[j], dE = W->dim_E[j], db = dL * dE;
        /* J_M W_j^dag : n x (dL*dE) (apply J_M to W_j^dag once per block). */
        acb_mat_t Wjd, JMWjd;
        acb_mat_init(Wjd, m, db);
        acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
        acb_mat_init(JMWjd, n, db);
        acb_mat_mul(JMWjd, d->J_M, Wjd, prec);
        for (slong c = 0; c < dE; c++) {
            acb_mat_init(out->K[kk], n, dL_tot);
            build_Kjc(out->K[kk], JMWjd, dL, dE, c, out->off_L[j], n);
            kk++;
        }
        acb_mat_clear(JMWjd);
        acb_mat_clear(Wjd);
    }

    /* (2) read the M^perp Sigma off the stored d->Delta (empty if M = H). */
    aic_dk_build_sigma(out, d, prec);

    /* (3) certify M-part unitality Delta_M(1_A) = Pi_M (fail-loud). */
    aic_dk_certify_unital(out, d, prec);

    /* (4) THE CRUX TOOTH: the full Delta-Kraus map == the stored d->Delta on the
     * A-basis B_k (M-part via K_{j,c} + M^perp-part via Sigma). Fail-loud on any
     * mismatch (.tex:2100). */
    aic_dk_certify_match(out, W, d, prec);
}

void aic_idemp_delta_clear(aic_idemp_delta *out)
{
    if (out == NULL || out->K == NULL) return;
    for (slong i = 0; i < out->n_kraus; i++) acb_mat_clear(out->K[i]);
    flint_free(out->K);
    if (out->dim_Mperp > 0) acb_mat_clear(out->Sigma);
    flint_free(out->dim_L);
    flint_free(out->dim_E);
    flint_free(out->off_L);
    out->K = NULL;
    out->dim_L = NULL;
    out->dim_E = NULL;
    out->off_L = NULL;
    out->n_kraus = 0;
    out->num_blocks = 0;
}
