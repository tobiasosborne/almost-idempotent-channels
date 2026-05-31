/* aic_factorize_upsilon3.c — the UCP DECODE map Upsilon, PART 3: the apply maps
 * Upsilon'_j / Upsilon' / Upsilon and the D5 F-ordering pin (bead aic-tff, module
 * "factorize", Increment F3). Split from aic_factorize_upsilon2.c (Rule 10,
 * <=200 LOC). See include/aic_factorize.h (the F3 banner) + the LOCKED decision
 * D5 (the F-ancilla ordering) in docs/research/factorize_f3_synthesis.md.
 *
 * approximate_algebras.tex:2866-2899 (verbatim):
 *   Upsilon'_j(X) = L_j^dag (Phi(X) (x) 1_F) L_j : B(H) -> B(L_j) (.tex:2869),
 *   written F-LEFT (1_F (x) Phi(X)) to match V's F-major Stinespring (D5).
 *   Upsilon' = (Upsilon'_1, ..., Upsilon'_m) (.tex:2867), the block-diagonal join.
 *   Upsilon: X |-> (Upsilon'(1_H))^{-1/2} Upsilon'(X) (Upsilon'(1_H))^{-1/2}
 *                                                      (.tex:2897), the UCP map.
 *
 * D5 — THE F-ANCILLA ORDERING (Designer C's "highest-risk", all three designers).
 * DEVIATES from .tex:2869's literal (Phi(X)(x)1_F) side — paper/FINDINGS.md §C13(a)
 * records the deviation + the r>1 evidence (Law 1: a .tex-formula deviation cites a
 * FINDINGS entry). aic_ucp_kraus_to_stinespring packs V F-MAJOR: V[a*N+i, j] = K_a[i,j] (a in F,
 * i in H/K, the row index is F-major, H-minor). With THIS layout
 *   Phi(X) = V^dag (1_F (x) X) V    (NOT V^dag (X (x) 1_F) V):
 *   [V^dag(1_F(x)X)V]_{j,l} = sum_{a,i,k} conj(K_a[i,j]) X[i,k] K_a[k,l]
 *                           = sum_a (K_a^dag X K_a)[j,l] = Phi(X)[j,l]. (checked.)
 * So EVERY (x) 1_F factor MUST be written F-LEFT (1_F (x) .) to match V's layout:
 * Upsilon'_j uses (1_F (x) Phi(X)); L_j uses (1_F (x) Delta(U_{js}^dag)). eta=0
 * with r=1 is BLIND (F is 1-dim, 1_F (x) M == M (x) 1_F); the eta>0 r>1 D5 pin
 * (aic_factorize_upsilon_d5_pin) compares both orderings — only F-LEFT reconstructs
 * Upsilon'_j Delta ~ 1_B (MEASURED: F-LEFT 4.4e-2 ~ O(eta) vs F-RIGHT 7.6e-1 ~ O(1),
 * mixconj(4,2,0.03) r=6; tests/test_factorize.c T6).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_dhom.h"
#include "aic/aic_factorize.h"
#include "aic/aic_factorize_internal.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

void aic_factorize_upsilon_prime_block(acb_mat_t out, const aic_factorize *F,
                                       slong j, const acb_mat_t X, slong prec)
{
    assert(F != NULL && F->upsilon_ready &&
           "upsilon_prime_block: aic_factorize_upsilon_build must run first");
    const aic_dhom_B *B = F->v->B;
    slong N = F->N, r = F->r, dj = B->d[j];
    assert(acb_mat_nrows(out) == dj && acb_mat_ncols(out) == dj &&
           "upsilon_prime_block: out must be d_j x d_j");
    assert(acb_mat_nrows(X) == N && acb_mat_ncols(X) == N && "upsilon_prime: X N x N");

    /* Upsilon'_j(X) = L_j^dag (1_F (x) Phi(X)) L_j (F-LEFT, D5). */
    acb_mat_t PhiX, PhiX_F, IF, Ld, tmp;
    acb_mat_init(PhiX, N, N);
    acb_mat_init(PhiX_F, r * N, r * N);
    acb_mat_init(IF, r, r);
    acb_mat_init(Ld, dj, r * N);
    acb_mat_init(tmp, dj, r * N);
    aic_ucp_apply(PhiX, F->phi, X, prec);             /* Phi(X) */
    acb_mat_one(IF);
    aic_mat_kronecker(PhiX_F, IF, PhiX, prec);        /* 1_F (x) Phi(X) (F-LEFT) */
    acb_mat_conjugate_transpose(Ld, F->L[j]);
    acb_mat_mul(tmp, Ld, PhiX_F, prec);
    acb_mat_mul(out, tmp, F->L[j], prec);
    acb_mat_clear(tmp); acb_mat_clear(Ld); acb_mat_clear(IF);
    acb_mat_clear(PhiX_F); acb_mat_clear(PhiX);
}

void aic_factorize_upsilon_prime(acb_mat_t out, const aic_factorize *F,
                                 const acb_mat_t X, slong prec)
{
    assert(F != NULL && F->upsilon_ready);
    const aic_dhom_B *B = F->v->B;
    assert(acb_mat_nrows(out) == F->n_B && acb_mat_ncols(out) == F->n_B &&
           "upsilon_prime: out must be n_B x n_B");
    acb_mat_zero(out);
    slong r_off = 0;
    for (slong j = 0; j < B->num_blocks; j++) {
        slong dj = B->d[j];
        acb_mat_t blk;
        acb_mat_init(blk, dj, dj);
        aic_factorize_upsilon_prime_block(blk, F, j, X, prec);
        for (slong i = 0; i < dj; i++)
            for (slong k = 0; k < dj; k++)
                acb_set(acb_mat_entry(out, r_off + i, r_off + k),
                        acb_mat_entry(blk, i, k));
        acb_mat_clear(blk);
        r_off += dj;
    }
}

void aic_factorize_upsilon(acb_mat_t out, const aic_factorize *F,
                           const acb_mat_t X, slong prec)
{
    assert(F != NULL && F->upsilon_ready);
    slong nB = F->n_B;
    assert(acb_mat_nrows(out) == nB && acb_mat_ncols(out) == nB &&
           "upsilon: out must be n_B x n_B");
    acb_mat_t UP, tmp;
    acb_mat_init(UP, nB, nB);
    acb_mat_init(tmp, nB, nB);
    aic_factorize_upsilon_prime(UP, F, X, prec);          /* Upsilon'(X) */
    acb_mat_mul(tmp, F->upsI_invsqrt, UP, prec);          /* M^{-1/2} Upsilon'(X) */
    acb_mat_mul(out, tmp, F->upsI_invsqrt, prec);         /* ... M^{-1/2} */
    acb_mat_clear(tmp); acb_mat_clear(UP);
}

/* D5 PIN (test/diagnostic): worst over B matrix units E_i of
 * ||Upsilon'_j(Delta(E_i))_j - (E_i)_j||_op, with the (x)1_F ordering chosen by
 * `f_left` for BOTH L_j and Upsilon'_j (V stays F-major). At eta>0 with r>1:
 * f_left=1 (CORRECT) gives O(eta); f_left=0 (F-RIGHT, WRONG) gives O(1) — the
 * decisive D5 distinguisher. Builds L_j fresh per call (does NOT touch the cached
 * F->L[j]). Returns the worst defect over all (j, i). REQUIRES F->upsilon_ready. */
double aic_factorize_upsilon_d5_pin(const aic_factorize *F, int f_left, slong prec)
{
    const aic_dhom_B *B = F->v->B;
    slong N = F->N, r = F->r, nB = F->n_B;
    double worst = 0.0;
    acb_mat_t *Lalt = flint_malloc(sizeof(acb_mat_t) * B->num_blocks);
    for (slong j = 0; j < B->num_blocks; j++) {
        acb_mat_init(Lalt[j], N * r, B->d[j]);
        aic_factorize_int_build_Lj(Lalt[j], F, j, F->W[j], F->xi[j], F->e[j],
                                   f_left, prec);
    }
    acb_mat_t Ei, DE, PhiDE, PhiDE_F, IF, Ld, tmp, blk, Eij, Mmid;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(DE, N, N);
    acb_mat_init(PhiDE, N, N);
    acb_mat_init(PhiDE_F, r * N, r * N);
    acb_mat_init(IF, r, r);
    acb_mat_one(IF);
    for (slong i = 0; i < B->dim_B; i++) {
        aic_dhom_B_matunit(Ei, B, i);
        aic_factorize_delta(DE, F, Ei, prec);             /* Delta(E_i) */
        aic_ucp_apply(PhiDE, F->phi, DE, prec);           /* Phi(Delta(E_i)) */
        if (f_left) aic_mat_kronecker(PhiDE_F, IF, PhiDE, prec);
        else        aic_mat_kronecker(PhiDE_F, PhiDE, IF, prec);
        slong r_off = 0;
        for (slong j = 0; j < B->num_blocks; j++) {
            slong dj = B->d[j];
            acb_mat_init(Ld, dj, r * N);
            acb_mat_init(tmp, dj, r * N);
            acb_mat_init(blk, dj, dj);
            acb_mat_init(Eij, dj, dj);
            acb_mat_init(Mmid, dj, dj);
            acb_mat_conjugate_transpose(Ld, Lalt[j]);
            acb_mat_mul(tmp, Ld, PhiDE_F, prec);
            acb_mat_mul(blk, tmp, Lalt[j], prec);          /* Upsilon'_j(Phi(Delta(E_i))) */
            for (slong a = 0; a < dj; a++)
                for (slong b = 0; b < dj; b++)
                    acb_set(acb_mat_entry(Eij, a, b),
                            acb_mat_entry(Ei, r_off + a, r_off + b));   /* (E_i)_j */
            acb_mat_sub(blk, blk, Eij, prec);
            for (slong a = 0; a < dj; a++)
                for (slong b = 0; b < dj; b++)
                    acb_get_mid(acb_mat_entry(Mmid, a, b), acb_mat_entry(blk, a, b));
            arb_t e;
            arb_init(e);
            aic_mat_opnorm(e, Mmid, prec);
            double d = arf_get_d(arb_midref(e), ARF_RND_NEAR);
            if (d > worst) worst = d;
            arb_clear(e);
            acb_mat_clear(Mmid); acb_mat_clear(Eij); acb_mat_clear(blk);
            acb_mat_clear(tmp); acb_mat_clear(Ld);
            r_off += dj;
        }
    }
    acb_mat_clear(IF); acb_mat_clear(PhiDE_F); acb_mat_clear(PhiDE);
    acb_mat_clear(DE); acb_mat_clear(Ei);
    for (slong j = 0; j < B->num_blocks; j++) acb_mat_clear(Lalt[j]);
    flint_free(Lalt);
    return worst;
}
