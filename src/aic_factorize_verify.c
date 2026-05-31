/* aic_factorize_verify.c — the end-to-end verification objects of th_factorization
 * (bead aic-tff, module "factorize", Increment F4.1). EXTENDS aic_factorize.c
 * (F1) / aic_factorize_delta.c (F2) / aic_factorize_upsilon*.c (F3). See
 * include/aic_factorize.h (the F4 banner) and docs/research/factorize_f4_design.md
 * §A (the composed-map Choi), §B (the certified headline), §D (the route-(i)
 * probe). FINDINGS §C12 (cb != Frobenius), §C13(c) (the m>=2 ∧ eta>0 debt).
 *
 * approximate_algebras.tex:2730-2739 — th_factorization (verbatim, the headline):
 *   "there are a finite-dimensional C* algebra B and UCP maps Delta: B -> B(H)
 *    and Upsilon: B(H) -> B such that
 *        ||Delta Upsilon - Phi||_cb <= O(eta),                        (.tex:2733 DelUps)
 *        ||Upsilon_n(Delta_n(X)Delta_n(Y)) - XY|| <= O(eta)||X|| ||Y||  (.tex:2736 UpsDel2)
 *    (In particular, (UpsDel2) implies that ||Upsilon Delta - 1_B||_cb <= O(eta).)" (.tex:2739)
 *
 * WHAT F4.1 BUILDS. The two Hermiticity-preserving (HP) DEFECT maps of the
 * theorem, as Convention-A Choi matrices, plus the always-valid eig-free certified
 * UPPER bound on their cb-norms:
 *
 *   J_DelUps = Choi(Delta Upsilon) - Choi(Phi)   (N^2 x N^2; .tex:2733)
 *   J_UpsDel = Choi(Upsilon Delta - 1_B)         (n_B^2 x n_B^2; .tex:2739)
 *
 * Delta Upsilon (W) = Delta(Upsilon(W)) is a self-map on B(H) = M_N; Upsilon Delta
 * (X) = Upsilon(Delta(X)) is a self-map on B = (+)_l M_{d_l}. Neither is an
 * aic_ucp_kraus object (each is a functional composition Delta o Upsilon resp.
 * Upsilon o Delta), so we build their Choi DIRECTLY by the E_pq loop (the F2
 * aic_factorize_delta_block_choi pattern, GENERALIZED to the FULL ambient matrix-
 * unit basis of M_N resp. M_{n_B}) rather than routing through aic_ucp_compose.
 *
 * THE AMBIENT E_pq UNITS (DRIFT #1, design §A.2 route (i)). The matrix units here
 * are the ambient M_N (resp. M_{n_B}) standard units E_pq (a zero matrix with a 1
 * at (p,q)) — NOT aic_dhom_B_matunit (which enumerates only the dim_B = Sum d_l^2
 * IN-B block-diagonal units). Route (i): view Upsilon Delta as a map M_{n_B} ->
 * M_{n_B} on the block-diagonal representative and build its Choi over all n_B^2
 * ambient units. Delta reads only the block-diagonal coordinates (aic_dhom_v_apply),
 * so the off-block-diagonal input columns of J_UpsDel are zero and the ambient
 * cb-norm equals the in-B cb-norm (design §A.2; the §D probe, T8, confirms this at
 * eta>0 since the eta=0 oracle is BLIND to it — Upsilon Delta = 1_B exactly).
 *
 * Convention A (aic_ucp.h:23-43): for a map f on M_n,
 *   Choi(f)[p*n + a, q*n + b] = f(E_pq)[a,b]   (INPUT p,q MAJOR/left, OUTPUT a,b MINOR).
 * For "-1_B" on input E_pq: 1_B acts as the identity MAP on M_{n_B}, so the term is
 * just -E_pq (subtract the ambient unit from Upsilon Delta(E_pq)).
 *
 * THE CERTIFIED BOUND (design §B.2). aic_cbnorm_eigfree_ball_choi gives the always-
 * valid bracket [||J||_F/n, 2||J||_F]; hi = 2||J||_F is a RIGOROUS per-instance
 * upper bound on the cb-norm (= diamond norm) of the HP defect map (aic_cbnorm.h:39).
 * At eta=0 the oracle gives J=0 -> [0,0]. The hi/eta ratio is a finite per-instance
 * constant but is NOT dimension-faithful (its 2N looseness is O(N)); the tight SDP
 * cb-norm and the dimension-independence canary are F4.2 (Julia+MOSEK), NOT here.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_cbnorm.h"
#include "aic_dhom.h"
#include "aic_factorize.h"
#include "aic_mat.h"
#include "aic_ucp.h"

/* Write the N_out x N_out output block of a map applied to ambient unit E_pq into
 * Convention-A position (p,q): J[p*N_out + a, q*N_out + b] = blk[a,b]. */
static void write_choi_block(acb_mat_t J, const acb_mat_t blk, slong p, slong q,
                             slong N_out)
{
    for (slong a = 0; a < N_out; a++)
        for (slong b = 0; b < N_out; b++)
            acb_set(acb_mat_entry(J, p * N_out + a, q * N_out + b),
                    acb_mat_entry(blk, a, b));
}

void aic_factorize_choi_delups(acb_mat_t J, const aic_factorize *F, slong prec)
{
    assert(F != NULL);
    assert(F->delta_ready && F->upsilon_ready &&
           "choi_delups: requires delta_ready && upsilon_ready (Rule 4)");
    slong N = F->N, nB = F->n_B;
    assert(acb_mat_nrows(J) == N * N && acb_mat_ncols(J) == N * N &&
           "choi_delups: J must be N^2 x N^2");

    /* Choi(Delta Upsilon) = Sum_{p,q} E_pq (x) Delta(Upsilon(E_pq)), over the N^2
     * ambient M_N units (DRIFT #1: NOT aic_dhom_B_matunit). */
    acb_mat_t Epq, U, DU;
    acb_mat_init(Epq, N, N);
    acb_mat_init(U, nB, nB);          /* Upsilon(E_pq): B(H) -> B  */
    acb_mat_init(DU, N, N);           /* Delta(Upsilon(E_pq)): B -> B(H) */
    acb_mat_zero(J);
    for (slong p = 0; p < N; p++)
        for (slong q = 0; q < N; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_upsilon(U, F, Epq, prec);   /* Upsilon (UCP decode) */
            aic_factorize_delta(DU, F, U, prec);      /* Delta   (UCP encode) */
            write_choi_block(J, DU, p, q, N);
        }
    acb_mat_clear(DU); acb_mat_clear(U); acb_mat_clear(Epq);

    /* Subtract Choi(Phi) (Convention A, N^2 x N^2). */
    acb_mat_t Cphi;
    acb_mat_init(Cphi, N * N, N * N);
    aic_ucp_kraus_to_choi(Cphi, F->phi, prec);
    acb_mat_sub(J, J, Cphi, prec);
    acb_mat_clear(Cphi);
}

void aic_factorize_choi_upsdel(acb_mat_t J, const aic_factorize *F, slong prec)
{
    assert(F != NULL);
    assert(F->delta_ready && F->upsilon_ready &&
           "choi_upsdel: requires delta_ready && upsilon_ready (Rule 4)");
    slong N = F->N, nB = F->n_B;
    const aic_dhom_B *B = F->v->B;
    assert(acb_mat_nrows(J) == nB * nB && acb_mat_ncols(J) == nB * nB &&
           "choi_upsdel: J must be n_B^2 x n_B^2 (route (i), ambient M_{n_B})");

    /* block_of[i] = the B-block containing ambient index i (B = (+)_l M_{d_l},
     * block-diagonal of size n_B = Sum_l d_l). */
    slong *block_of = (slong *) flint_malloc((size_t) nB * sizeof(slong));
    {
        slong off = 0;
        for (slong l = 0; l < B->num_blocks; l++)
            for (slong a = 0; a < B->d[l]; a++) block_of[off++] = l;
    }

    /* Choi(Upsilon Delta - 1_B) = Sum_{p,q} E_pq (x) (Upsilon(Delta(E_pq)) - 1_B(E_pq)),
     * over the n_B^2 ambient M_{n_B} units (DRIFT #1, route (i)). THE 1_B TERM is the
     * CONDITIONAL EXPECTATION P_B onto B (the block-diagonal projection), NOT the full
     * ambient identity: 1_B is the unit of B, and Upsilon Delta lands in B, so the
     * cb-norm object Upsilon Delta - 1_B as a self-map on M_{n_B} subtracts P_B(E_pq)
     * = E_pq when (p,q) is IN a single block, 0 OFF-block. (The design's "-E_pq for
     * every E_pq" wording is loose: it would leave -E_pq on every off-block input
     * column, making ||J_UpsDel||_F = sqrt(#off-block units) != 0 at eta=0 even though
     * Upsilon Delta = 1_B exactly — a "test that cannot pass". The correct 1_B = P_B
     * makes J_UpsDel = 0 at eta=0 AND keeps the off-block input columns = Upsilon
     * Delta(off-block) ~ 0, the §D probe. See REPORT / FINDINGS.) */
    acb_mat_t Epq, D, UD, blk;
    acb_mat_init(Epq, nB, nB);
    acb_mat_init(D, N, N);            /* Delta(E_pq): B -> B(H) */
    acb_mat_init(UD, nB, nB);         /* Upsilon(Delta(E_pq)): B(H) -> B */
    acb_mat_init(blk, nB, nB);
    acb_mat_zero(J);
    for (slong p = 0; p < nB; p++)
        for (slong q = 0; q < nB; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_delta(D, F, Epq, prec);     /* Delta (UCP encode) */
            aic_factorize_upsilon(UD, F, D, prec);    /* Upsilon (UCP decode) */
            acb_mat_set(blk, UD);
            if (block_of[p] == block_of[q])           /* P_B(E_pq) = E_pq in-block */
                acb_sub_si(acb_mat_entry(blk, p, q), acb_mat_entry(blk, p, q), 1, prec);
            write_choi_block(J, blk, p, q, nB);
        }
    acb_mat_clear(blk); acb_mat_clear(UD); acb_mat_clear(D); acb_mat_clear(Epq);
    flint_free(block_of);
}

void aic_factorize_eigfree_delups(arb_t lo, arb_t hi, const aic_factorize *F,
                                  slong prec)
{
    assert(F != NULL);
    slong N = F->N;
    acb_mat_t J;
    acb_mat_init(J, N * N, N * N);
    aic_factorize_choi_delups(J, F, prec);
    aic_cbnorm_eigfree_ball_choi(lo, hi, J, N, prec);   /* hi = 2||J||_F >= ||DelUps||_cb */
    acb_mat_clear(J);
}

void aic_factorize_eigfree_upsdel(arb_t lo, arb_t hi, const aic_factorize *F,
                                  slong prec)
{
    assert(F != NULL);
    slong nB = F->n_B;
    acb_mat_t J;
    acb_mat_init(J, nB * nB, nB * nB);
    aic_factorize_choi_upsdel(J, F, prec);
    aic_cbnorm_eigfree_ball_choi(lo, hi, J, nB, prec);  /* hi = 2||J||_F >= ||UpsDel||_cb */
    acb_mat_clear(J);
}
