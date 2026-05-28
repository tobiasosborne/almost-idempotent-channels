/* aic_ucp_compose.c — UCP-map composition and Choi-difference (bead aic-d24,
 * the first increment of the eta-idempotence defect ||Phi^2-Phi||_cb,
 * approximate_algebras.tex:347-354). Split out of aic_ucp_core.c (near the
 * 200-LOC budget) per CLAUDE.md Rule 10.
 *
 * WHAT THIS COMPUTES, AND WHY IT STOPS HERE.
 *   The paper calls Phi `eta-idempotent` if ||Phi^2-Phi||_cb <= eta
 *   (approximate_algebras.tex:354, verbatim:
 *     "A UCP map Phi: B(H)->B(H) is called eta-idempotent if
 *      ||Phi^2-Phi||_cb <= eta.").
 *   Lambda = Phi^2 - Phi is Hermiticity-preserving but NOT completely positive,
 *   so the CP closed form ||Phi||_cb = ||Phi(1)||_op (.tex:1717-1719, used by
 *   aic_ucp_cbnorm_cp) does NOT apply. ||Lambda||_cb = ||Lambda||_diamond is
 *   computed by a Watrous SDP (.tex:352, verbatim: "The completely bounded norm
 *   in finite dimensions is efficiently computable using semidefinite
 *   programming"). No arb-native SDP solver exists, so per the bead-aic-d24
 *   Route B the SDP lives in Julia+MOSEK; this C core supplies the SDP its input
 *   — the Choi matrix of Lambda — and nothing more. The diamond-norm VALUE is the
 *   Julia golden master (validated against analytic anchors), cross-checked here
 *   only at the Choi level (test_ucp_d24.c).
 *
 * COMPOSITION (Heisenberg picture). Derived in aic_ucp.h above the declaration:
 *   (Phi o Psi)(X) = sum_a K_a^dag (sum_b L_b^dag X L_b) K_a
 *                  = sum_{a,b} (L_b K_a)^dag X (L_b K_a),
 *   so Kraus(Phi o Psi) = { L_b K_a } with the matrix product L_b K_a =
 *   psi.K[b] @ phi.K[a]. Requires phi->dim_K == psi->dim_H. This is PURE LINEAR
 *   ALGEBRA (r1*r2 matrix products) — constructive, exact up to rounding, both
 *   number paths. Phi^2 = compose(Phi,Phi). The cross-check in the tests is the
 *   definition itself: aic_ucp_apply(compose(Phi,Psi), X) == Phi(Psi(X)).
 *
 * CHOI DIFFERENCE. Choi is linear in the map (the convention-A formula
 *   C[i*h+a, j*h+b] = sum_x conj(K_x[i,a]) K_x[j,b] is linear in the rank-one
 *   outer products), so Choi(Phi^2 - Phi) = Choi(Phi^2) - Choi(Phi). We build
 *   each Choi via the existing aic_ucp_kraus_to_choi and subtract. The result is
 *   Hermitian (difference of Hermitian Choi matrices) but generically indefinite.
 *
 * Number paths. Both routines are arb/acb here (acb_mat_mul / aic_ucp_kraus_to_
 * choi at an explicit prec). The double-vs-arb@53 cross-check (CLAUDE.md ladder
 * rung 2) is driven from tests/test_ucp_d24.c: the composed Kraus operators and
 * the Choi-diff entries are recomputed on the LAPACK double path (aic_latd matrix
 * products / direct double Choi formula) and compared within a 53-bit tol.
 *
 * Fail-loud (Rule 4). dim_K/dim_H compatibility is asserted with a message at
 * the call site; a mismatch aborts rather than producing a silently wrong map.
 */
#include <assert.h>

#include <flint/acb_mat.h>

#include "aic_ucp.h"

/* aic_ucp.h: Kraus(Phi o Psi) = { L_b K_a } = { psi.K[b] @ phi.K[a] },
 * requires phi->dim_K == psi->dim_H; out shape (psi->dim_K, phi->dim_H),
 * r = phi->r * psi->r. */
void aic_ucp_compose(aic_ucp_kraus *out, const aic_ucp_kraus *phi,
                     const aic_ucp_kraus *psi, slong prec)
{
    assert(out != NULL && phi != NULL && psi != NULL &&
           "aic_ucp_compose: null argument");
    /* Type check: Psi outputs in B(K_psi), Phi accepts B(K_phi), so the dims
     * must chain — phi's domain dimension equals psi's codomain dimension. */
    assert(phi->dim_K == psi->dim_H &&
           "aic_ucp_compose: phi->dim_K must equal psi->dim_H (domains must "
           "chain for Phi o Psi)");

    slong out_dim_K = psi->dim_K;
    slong out_dim_H = phi->dim_H;
    slong out_r = phi->r * psi->r;
    aic_ucp_kraus_init(out, out_dim_K, out_dim_H, out_r);

    /* L_b K_a : (psi->dim_K x psi->dim_H) * (phi->dim_K x phi->dim_H), with
     * psi->dim_H == phi->dim_K (asserted), giving psi->dim_K x phi->dim_H. The
     * index ordering c = a*psi->r + b is arbitrary (Kraus ops are a set up to
     * gauge); we keep it fixed so the count is exactly phi->r * psi->r. */
    slong c = 0;
    for (slong a = 0; a < phi->r; a++)
        for (slong b = 0; b < psi->r; b++) {
            acb_mat_mul(out->K[c], psi->K[b], phi->K[a], prec); /* L_b K_a */
            c++;
        }
    assert(c == out_r && "aic_ucp_compose: Kraus count mismatch");
}

/* aic_ucp.h: C = Choi(phi1) - Choi(phi2), Convention A; phi1, phi2 must share
 * (dim_K, dim_H); C initialised (dim_K*dim_H) square. */
void aic_ucp_choi_diff(acb_mat_t C, const aic_ucp_kraus *phi1,
                       const aic_ucp_kraus *phi2, slong prec)
{
    assert(phi1 != NULL && phi2 != NULL && "aic_ucp_choi_diff: null argument");
    assert(phi1->dim_K == phi2->dim_K && phi1->dim_H == phi2->dim_H &&
           "aic_ucp_choi_diff: phi1 and phi2 must share (dim_K, dim_H)");
    slong n = phi1->dim_K * phi1->dim_H;
    assert(acb_mat_nrows(C) == n && acb_mat_ncols(C) == n &&
           "aic_ucp_choi_diff: C must be (dim_K*dim_H) square");

    acb_mat_t C1, C2;
    acb_mat_init(C1, n, n);
    acb_mat_init(C2, n, n);
    aic_ucp_kraus_to_choi(C1, phi1, prec);
    aic_ucp_kraus_to_choi(C2, phi2, prec);
    acb_mat_sub(C, C1, C2, prec); /* Choi(phi1) - Choi(phi2), Hermitian */
    acb_mat_clear(C2);
    acb_mat_clear(C1);
}
