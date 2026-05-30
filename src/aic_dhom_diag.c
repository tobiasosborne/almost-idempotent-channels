/* aic_dhom_diag.c — the EXPLICIT generalized-Pauli diagonal of the genuine
 * finite-dim C* algebra B = (+)_l M_{d_l} (bead aic-c1n, module "dhom"). Realizes
 * the diagonal construction of approximate_algebras.tex:1245-1254 (inside the
 * proof of lem_approx :1224). Split from aic_dhom_B.c for the ~200 LOC limit
 * (Rule 10): the B-algebra primitives live in aic_dhom_B.c; the diagonal here.
 *
 * WHY THE PAULI DIAGONAL (NOT a Haar average). The paper's diagonal of a finite-
 * dim C* algebra is D = integral dU (U^dag (x) U) (Haar, .tex:1245). That is
 * non-constructive as written. The constructive route (CLAUDE.md Law 3) is the
 * EXPLICIT closed form .tex:1250-1254 via generalized Pauli (shift-clock /
 * Heisenberg-Weyl) operators S_{jk} = X^j Z^k:
 *     D_l = sum_{j,k=0}^{d-1} d^-2 S_{jk}^dag (x) S_{jk}.
 * This is EXACT: pi(D_l) = sum d^-2 S_{jk}^dag S_{jk} = I_d, X D_l = D_l X for all
 * X in M_d (centrality, since {S_{jk}} is an operator basis / irreducible
 * projective rep), and ||D_l||_proj = sum d^-2 ||S_{jk}^dag|| ||S_{jk}|| = 1 (each
 * S_{jk} unitary). Exactness is WHAT MAKES THE th_main CONSTANT DIMENSION-
 * INDEPENDENT (.tex:484): a naive eps-associative-averaged diagonal would carry an
 * O(n) error; the genuine-C*-algebra Pauli diagonal carries none. The diagonal of
 * a DIRECT SUM (+)_l M_{d_l} is the cross-term-free EMBEDDED SUM D = sum_l D_l (see
 * the .tex:1254 DISCREPANCY box below — the paper's Cartesian-product/joint-unitary
 * formula is NON-central for the finite Pauli 1-design and is NOT what we build).
 * sum_l d_l^2 terms total (a SUM over blocks, not a product over them).
 *
 * NO eigendecomposition: pure complex-exponential fills + block embedding. The
 * diagonal identities (centrality, pi(D)=I, ||D||_proj) are checked in the tests
 * via the certified path (test_dhom T1).
 *
 * ============================== .tex:1254 DISCREPANCY (escalate) ==============
 * The paper's tex:1254 direct-sum formula is the CARTESIAN PRODUCT over per-block
 * Pauli tuples with the JOINT block-diagonal unitary U_{j_1..j_m}=U_{1j_1}(+)..(+)U_{mj_m}
 * and weight prod_l p_{lj_l}. WE DO NOT IMPLEMENT THAT FORMULA: for the finite
 * Pauli 1-design it is NOT a diagonal — it FAILS centrality (X D != D X, .tex:1241).
 *
 * ROOT CAUSE. A diagonal is D = integral dU (U^dag (x) U) over the Haar measure on
 * B's unitary group = prod_l U(d_l) (block-diagonal unitaries). Expanding the joint
 * unitary U=(+)_l U_l gives U^dag (x) U = sum_{a,b} (U_a^dag in sector a) (x)
 * (U_b in sector b). For the HAAR integral the cross terms a!=b VANISH because the
 * unitary FIRST moment is zero (integral_{U(d)} U dU = 0). The implementer's
 * substitution of the joint Cartesian product of finite Pauli 1-designs does NOT
 * reproduce the Haar integral: the Pauli FIRST moment sum_{jk} S_{jk} is a NONZERO
 * rank-1 matrix (unlike Haar), so the joint product leaves SPURIOUS cross-sector
 * terms -> non-central (measured ||XD-DX|| ~ 0.5-1.0 for M_2(+)M_2). The per-block
 * Pauli sum sum_{jk} d^-2 S_{jk}^dag (x) S_{jk} DOES reproduce the SINGLE-block Haar
 * second moment (single block IS correct), but the JOINT product does not.
 *
 * THE CORRECT FINITE CENTRAL DIAGONAL (what we build):
 *     D = sum_l D_l = sum_l sum_{j,k} d_l^-2 (Shat^{(l)}_{jk})^dag (x) Shat^{(l)}_{jk},
 * where Shat^{(l)}_{jk} is the n_B x n_B matrix equal to S^{(l)}_{jk} on block l and
 * ZERO elsewhere (a PARTIAL ISOMETRY embedded in sector l — NOT the joint unitary,
 * NOT (+)-with-identity). This is the Haar diagonal D=integral dU (U^dag (x) U): the
 * cross terms are absent BY CONSTRUCTION (each term lives in one sector). It is
 * central (X D = D X: for X=(+)_l X_l, X.Shat^{(l)} has support only in block l, and
 * within block l the per-block centrality holds), and pi(D)=sum_l sum_{jk} d_l^-2
 * (Shat)^dag Shat = sum_l (block-l identity) = I_B. nterms = SUM_l d_l^2 (NOT the
 * product prod_l d_l^2). The single-block case (m=1) falls out unchanged.
 *
 * CONSEQUENCE for ||D||_proj: each block contributes sum_{jk} d_l^-2 ||Shat^dag||
 * ||Shat|| = sum_{jk} d_l^-2 = 1, so ||D||_proj = sum_l 1 = m (the NUMBER OF
 * BLOCKS), NOT 1. (The norm-1 diagonal of .tex:1247-1249 is unique; in THIS
 * cross-term-free product representation the projective norm is m because the
 * cross-sector terms that would compress it to 1 are exactly the non-central
 * spurious terms we removed. The element D itself is still THE diagonal.) This
 * pushes the lem_approx w' bound (.tex:1281, ||w'(X)|| <= sum_t ||A_t|| ||B_t|| *
 * O(delta) = m * O(delta)) to O(m*delta) — a NUMBER-OF-BLOCKS dependence to watch
 * against the universal-constant claim .tex:484/460 (tracked by test_dhom T5; a
 * .tex-formula concern flagged for escalation, NOT silently papered over).
 * =============================================================================
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_dhom.h"
#include "aic_dhom_internal.h"

void aic_dhom_pauli(acb_mat_t out, slong d, slong j, slong k, slong prec)
{
    /* approximate_algebras.tex:1252:
     *   <e_l | S_{jk} e_m> = exp(2 pi i k m / d) * delta_{l, (m+j) mod d}.
     * So column m has a single nonzero entry exp(2 pi i k m / d) at row (m+j)%d.
     */
    assert(acb_mat_nrows(out) == d && acb_mat_ncols(out) == d);
    acb_mat_zero(out);
    arb_t two_pi, ang, c, s;
    arb_init(two_pi);
    arb_init(ang);
    arb_init(c);
    arb_init(s);
    arb_const_pi(two_pi, prec);
    arb_mul_si(two_pi, two_pi, 2, prec); /* 2 pi */
    for (slong m = 0; m < d; m++) {
        slong row = (m + j) % d;
        /* angle = 2 pi k m / d */
        arb_mul_si(ang, two_pi, k * m, prec);
        arb_div_si(ang, ang, d, prec);
        arb_sin_cos(s, c, ang, prec); /* s=sin, c=cos */
        acb_set_arb_arb(acb_mat_entry(out, row, m), c, s); /* cos + i sin */
    }
    arb_clear(s);
    arb_clear(c);
    arb_clear(ang);
    arb_clear(two_pi);
}

/* Embed a d x d block (block index l) into an n_B x n_B block-diagonal matrix. */
static void embed_block(acb_mat_t out, const aic_dhom_B *B, slong l,
                        const acb_mat_t blk)
{
    slong d = B->d[l], r0 = aic_dhom_block_row_offset(B, l);
    for (slong a = 0; a < d; a++)
        for (slong b = 0; b < d; b++)
            acb_set(acb_mat_entry(out, r0 + a, r0 + b), acb_mat_entry(blk, a, b));
}

void aic_dhom_diag_build(aic_dhom_diag *out, const aic_dhom_B *B, slong prec)
{
    /* The cross-term-free EMBEDDED SUM D = sum_l D_l (see the .tex:1254 DISCREPANCY
     * box at the top of this file — NOT the paper's joint Cartesian product, which
     * is non-central for the finite Pauli design). Each term is one (l, j, k):
     * U_t = Shat^{(l)}_{jk} (the d_l x d_l Pauli S^{(l)}_{jk} embedded in sector l,
     * zero elsewhere — a partial isometry), with weight w_t = d_l^-2. Then
     *   D = sum_l sum_{j,k} d_l^-2 (Shat^{(l)}_{jk})^dag (x) Shat^{(l)}_{jk}
     *     = integral dU (U^dag (x) U)   (the Haar diagonal; cross terms absent).
     * nterms = SUM_l d_l^2 (single block m=1 -> d^2, matching .tex:1251). */
    slong m = B->num_blocks;
    slong nterms = 0;
    for (slong l = 0; l < m; l++) nterms += B->d[l] * B->d[l];
    out->nterms = nterms;
    out->w = _arb_vec_init(nterms);
    out->U = flint_malloc((size_t) nterms * sizeof(acb_mat_t));

    arb_t weight, dl2;
    acb_mat_t S;
    arb_init(weight);
    arb_init(dl2);
    slong t = 0;
    for (slong l = 0; l < m; l++) {
        slong d = B->d[l];
        /* per-block weight d_l^-2 (.tex:1251), constant across this block's terms. */
        arb_set_si(dl2, d * d);
        arb_inv(weight, dl2, prec);
        acb_mat_init(S, d, d);
        for (slong jj = 0; jj < d; jj++)
            for (slong kk = 0; kk < d; kk++) {
                aic_dhom_pauli(S, d, jj, kk, prec); /* S^{(l)}_{jk} (d x d) */
                arb_set(out->w + t, weight);
                acb_mat_init(out->U[t], B->n_B, B->n_B);
                acb_mat_zero(out->U[t]);
                embed_block(out->U[t], B, l, S); /* Shat: S in sector l, 0 else */
                t++;
            }
        acb_mat_clear(S);
    }
    assert(t == nterms);
    arb_clear(dl2);
    arb_clear(weight);
}

void aic_dhom_diag_clear(aic_dhom_diag *out)
{
    if (out == NULL) return;
    for (slong t = 0; t < out->nterms; t++) acb_mat_clear(out->U[t]);
    flint_free(out->U);
    _arb_vec_clear(out->w, out->nterms);
    out->U = NULL;
}
