/* aic_factorize_delta.c — the UCP ENCODE map Delta of th_factorization (bead
 * aic-tff, module "factorize", Increment F2). EXTENDS aic_factorize.c (F1). See
 * include/aic_factorize.h (the F2 banner) and docs/research/factorize_d4_research.md
 * §1 Step 4, §2 (D4-a).
 *
 * approximate_algebras.tex:2786-2801 — th_factorization Step 4 (verbatim):
 *   "we define a map Delta': B -> B(H) by the equation
 *      Delta'(X) = sum_s p_s Phi( Delta~(X U_s^dag) Delta~(U_s) )           (.tex:2788)
 *    It is evident that Delta' commutes with the involution. The complete
 *    positivity is shown as follows: if X in 1_M_n (x) B is positive, X = Y^dag Y,
 *    and Delta'_n(Y^dag Y) = sum_s p_s Phi_n( Delta~_n(Y^dag(I (x) U_s^dag))
 *    Delta~_n((I (x) U_s) Y) ) >= 0.                                  (.tex:2791-2795)
 *    ... Delta: X |-> (Delta'(I_B))^{-1/2} Delta'(X) (Delta'(I_B))^{-1/2}    (.tex:2799)
 *    is a UCP map such that ||Delta - Delta~||_cb <= O(eta)."
 *
 * THE 1-DESIGN (.tex:2771-2784). Per block j, B(L_j) = M_{d_j}: the design is the
 * generalized Pauli (Heisenberg-Weyl) set U_{j,(a,b)} = S_ab = X^a Z^b with weight
 * p_{j,(a,b)} = d_j^{-2} (d_j^2 terms; sum = 1; centrality diag_j2 .tex:2776). The
 * diagonal of the WHOLE B is the CARTESIAN PRODUCT (.tex:2780-2784): index
 * s = (s_1,...,s_m), weight p_s = prod_l d_l^{-2}, JOINT unitary
 *   U_s = U_{1,s_1} (+) ... (+) U_{m,s_m}   (the n_B x n_B block-diagonal JOIN).
 *
 * THE §A2 DISTINCTION (FINDINGS §A2 — the landmine). F2's U_s is the GENUINE
 * per-block Pauli S_ab (aic_dhom_pauli) joined block-diagonally — a real unitary
 * in B (.tex:2782-2783). This is NOT aic_dhom_diag_build's output (the cross-term-
 * free EMBEDDED partial-isometry sum D = sum_l D_l, for lem_approx's error
 * reduction — a different purpose). We build U_s here ourselves and ASSERT each is
 * unitary (the §A2 sanity gate). A broken design (a single Pauli instead of the
 * full Cartesian sum) breaks CP-ness of Delta' / makes ||Delta'-Delta~|| = O(1).
 *
 * WHY Delta' IS CP BUT Delta~ IS NOT (the F2 payoff, §C2-style teeth). Delta~ = v
 * (an exact *-hom into A) is positive only w.r.t. A's star, not as a B(H) map. The
 * 1-design average routes the ORDINARY M_N product Delta~(.)Delta~(.) through the
 * ORIGINAL UCP Phi (aic_ucp_apply, NOT the superop Phi~ — Phi being genuinely CP
 * is what makes Delta'(Y^dag Y) = sum_s p_s Phi(Z_s^dag Z_s) >= 0, .tex:2791-2795).
 * The CP teeth: the per-block Choi C_j must be PSD; this catches both a broken
 * design and substituting Phi~ (not CP) for Phi.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_dhom.h"
#include "aic/aic_factorize.h"
#include "aic/aic_funcalc.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

slong aic_factorize_design_nterms(const aic_factorize *F)
{
    assert(F != NULL);
    const aic_dhom_B *B = F->v->B;
    slong nterms = 1;
    for (slong l = 0; l < B->num_blocks; l++)
        nterms *= B->d[l] * B->d[l];   /* prod_l d_l^2 */
    return nterms;
}

/* U_s = U_{1,s_1} (+) ... (+) U_{m,s_m), the block-diagonal join of per-block
 * Paulis. Decompose the flat index s (mixed-radix, block l radix d_l^2; lowest-
 * order block FIRST) into a per-block Pauli (a,b) = (q / d_l, q % d_l). Each
 * S_ab is placed in sector l of the n_B x n_B representative (offset r_off =
 * sum_{l'<l} d_{l'}). Off-block entries stay zero. */
void aic_factorize_design_unitary(acb_mat_t out, const aic_factorize *F,
                                  slong s, slong prec)
{
    assert(F != NULL);
    const aic_dhom_B *B = F->v->B;
    assert(acb_mat_nrows(out) == B->n_B && acb_mat_ncols(out) == B->n_B &&
           "design_unitary: out must be n_B x n_B");
    assert(s >= 0 && s < aic_factorize_design_nterms(F) &&
           "design_unitary: index s out of range");

    acb_mat_zero(out);
    slong rem = s, r_off = 0;
    for (slong l = 0; l < B->num_blocks; l++) {
        slong d = B->d[l], dd = d * d;
        slong q = rem % dd;           /* this block's Pauli index in [0, d^2) */
        rem /= dd;
        slong a = q / d, b = q % d;   /* S_ab = X^a Z^b on C^d */
        acb_mat_t S;
        acb_mat_init(S, d, d);
        aic_dhom_pauli(S, d, a, b, prec);
        for (slong i = 0; i < d; i++)
            for (slong k = 0; k < d; k++)
                acb_set(acb_mat_entry(out, r_off + i, r_off + k),
                        acb_mat_entry(S, i, k));
        acb_mat_clear(S);
        r_off += d;
    }
    assert(rem == 0 && "design_unitary: mixed-radix decomposition overflow");

    /* §A2 sanity gate (Rule 4, fail loud): U_s^dag U_s = I_{n_B}. */
    {
        slong nB = B->n_B;
        acb_mat_t Ud, G, I;
        acb_mat_init(Ud, nB, nB);
        acb_mat_init(G, nB, nB);
        acb_mat_init(I, nB, nB);
        acb_mat_conjugate_transpose(Ud, out);
        acb_mat_mul(G, Ud, out, prec);
        acb_mat_one(I);
        acb_mat_sub(G, G, I, prec);
        arb_t e, tol;
        arb_init(e);
        arb_init(tol);
        aic_mat_opnorm(e, G, prec);
        arb_set_d(tol, 1e-9);
        assert(arb_lt(e, tol) &&
               "design_unitary: U_s not unitary (broken 1-design? §A2)");
        arb_clear(tol);
        arb_clear(e);
        acb_mat_clear(I); acb_mat_clear(G); acb_mat_clear(Ud);
    }
}

void aic_factorize_delta_prime(acb_mat_t out, const aic_factorize *F,
                               const acb_mat_t X, slong prec)
{
    assert(F != NULL);
    slong N = F->N, nB = F->n_B;
    assert(acb_mat_nrows(out) == N && acb_mat_ncols(out) == N &&
           "delta_prime: out must be N x N");
    assert(acb_mat_nrows(X) == nB && acb_mat_ncols(X) == nB &&
           "delta_prime: X must be n_B x n_B (B's block-diagonal rep)");

    slong nterms = aic_factorize_design_nterms(F);

    /* p_s = prod_l d_l^{-2} (uniform over all terms). */
    arb_t ps;
    arb_init(ps);
    arb_set_si(ps, 1);
    arb_div_si(ps, ps, nterms, prec);    /* 1 / prod_l d_l^2 */

    acb_mat_t Us, Usd, XUsd, dXUsd, dUs, prod, phiprod, term;
    acb_mat_init(Us, nB, nB);
    acb_mat_init(Usd, nB, nB);
    acb_mat_init(XUsd, nB, nB);
    acb_mat_init(dXUsd, N, N);
    acb_mat_init(dUs, N, N);
    acb_mat_init(prod, N, N);
    acb_mat_init(phiprod, N, N);
    acb_mat_init(term, N, N);

    acb_mat_zero(out);
    for (slong s = 0; s < nterms; s++) {
        aic_factorize_design_unitary(Us, F, s, prec);
        acb_mat_conjugate_transpose(Usd, Us);
        aic_dhom_B_mul(XUsd, F->v->B, X, Usd, prec);     /* X U_s^dag (in B) */
        aic_factorize_delta_tilde(dXUsd, F, XUsd, prec); /* Delta~(X U_s^dag) */
        aic_factorize_delta_tilde(dUs, F, Us, prec);     /* Delta~(U_s) */
        acb_mat_mul(prod, dXUsd, dUs, prec);             /* ORDINARY M_N product */
        aic_ucp_apply(phiprod, F->phi, prod, prec);      /* Phi(...) — ORIGINAL */
        acb_mat_scalar_mul_arb(term, phiprod, ps, prec);
        acb_mat_add(out, out, term, prec);
    }

    acb_mat_clear(term); acb_mat_clear(phiprod); acb_mat_clear(prod);
    acb_mat_clear(dUs); acb_mat_clear(dXUsd); acb_mat_clear(XUsd);
    acb_mat_clear(Usd); acb_mat_clear(Us);
    arb_clear(ps);
}

void aic_factorize_delta_block_choi(acb_mat_t Cj, const aic_factorize *F,
                                    slong j, slong prec)
{
    assert(F != NULL);
    const aic_dhom_B *B = F->v->B;
    assert(j >= 0 && j < B->num_blocks && "block_choi: block index out of range");
    slong N = F->N, dj = B->d[j];
    assert(acb_mat_nrows(Cj) == dj * N && acb_mat_ncols(Cj) == dj * N &&
           "block_choi: Cj must be (d_j*N) x (d_j*N)");

    /* C_j = sum_{a,b in [0,d_j)} E_ab (x) Delta'(E^{(j)}_ab), Convention A
     * (left/major = the d_j source factor): Cj[a*N+p, b*N+q] = Delta'(E_ab)[p,q]. */
    acb_mat_zero(Cj);
    acb_mat_t Eab, DE;
    acb_mat_init(Eab, B->n_B, B->n_B);
    acb_mat_init(DE, N, N);
    for (slong a = 0; a < dj; a++)
        for (slong b = 0; b < dj; b++) {
            slong idx = aic_dhom_B_index(B, j, a, b);
            aic_dhom_B_matunit(Eab, B, idx);             /* E^{(j)}_ab in B */
            aic_factorize_delta_prime(DE, F, Eab, prec); /* Delta'(E^{(j)}_ab) */
            for (slong p = 0; p < N; p++)
                for (slong q = 0; q < N; q++)
                    acb_set(acb_mat_entry(Cj, a * N + p, b * N + q),
                            acb_mat_entry(DE, p, q));
        }
    acb_mat_clear(DE); acb_mat_clear(Eab);
}
