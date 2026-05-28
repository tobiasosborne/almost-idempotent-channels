/* aic_ucp_carrier.c — the carrier of a UCP map on the certified arb path
 * (bead aic-c7n): the Stinespring isometry V, the carrier operator Q, and the
 * checks for cor_carrier (.tex:1731) and PhiX_M (.tex:1740).
 *
 * Paper result vs constructive route.
 *  - Stinespring V (.tex:1573, V:H->K(x)F isometry): the column-stack of the
 *    Kraus ops, V[a*dim_K + i, j] = K_a[i,j]. Unitality V^dag V = 1_H (.tex:1571)
 *    is verified by the tests via aic_mat_opnorm on (V^dag V - I).
 *  - Carrier (lem_carrier, .tex:1724): "the smallest M with M (x) F >= Im V".
 *    The paper's proof is via the full-rank-state support of Phi^*(rho_0) and is
 *    an existence statement. CONSTRUCTIVE finite-dim route: the carrier is the
 *    K-marginal of Im V, i.e. the support (range) of Q = Tr_F(V V^dag). Working
 *    that out in Kraus coords, (V V^dag)[(i,a),(j,b)] = (K_a K_b^dag)[i,j], so the
 *    F-partial-trace (a=b) gives
 *        Q = sum_a K_a K_a^dag          (dim_K x dim_K, Hermitian PSD).
 *    This routine returns Q; its RANGE (the carrier itself) needs a degenerate
 *    eigensolver and lives on the double path (aic_ucp_carrier_rank_latd). The
 *    shard's first carrier formula was sum_a K_a^dag K_a (= 1_H by unitality, the
 *    WRONG marginal); the shard self-corrects and this Q = sum_a K_a K_a^dag is
 *    the right one (CLAUDE.md Rule 2: be skeptical of shard summaries).
 *  - cor_carrier (.tex:1731-1732): X annihilates M iff (X (x) 1_F) V = 0. In
 *    Kraus coords (X (x) 1_F) V is the column-stack of (X K_a), so the defect is
 *    ||stack_a(X K_a)||_op; this routine returns it (0 iff Ker X >= M).
 *  - PhiX_M (.tex:1740): Phi(X) = Phi(Pi_M X Pi_M) for ANY UCP map (proof
 *    .tex:1742 inserts Pi_M (x) 1_F on both sides of V, using (Pi_M (x) 1_F)V=V).
 *    This routine returns ||Phi(X) - Phi(Pi_M X Pi_M)||_op (0 up to rounding).
 *
 * Precision: Q and the defects are O(1) differences; prec=53 resolves them to
 * ~1e-15. The certified carrier RANK is gap-dependent (bead aic-w4o.1) and is
 * NOT attempted here.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_mat.h"
#include "aic_ucp.h"

/* .tex:1573 — V[a*dim_K + i, j] = K_a[i, j] (column-stack), V: H -> K (x) F. */
void aic_ucp_kraus_to_stinespring(acb_mat_t V, const aic_ucp_kraus *phi,
                                  slong prec)
{
    (void) prec;
    slong dk = phi->dim_K, dh = phi->dim_H, r = phi->r;
    assert(acb_mat_nrows(V) == dk * r && acb_mat_ncols(V) == dh &&
           "aic_ucp_kraus_to_stinespring: V must be (dim_K*r) x dim_H");
    for (slong a = 0; a < r; a++)
        for (slong i = 0; i < dk; i++)
            for (slong j = 0; j < dh; j++)
                acb_set(acb_mat_entry(V, a * dk + i, j),
                        acb_mat_entry(phi->K[a], i, j));
}

/* lem_carrier (.tex:1724) — Q = sum_a K_a K_a^dag = Tr_F(V V^dag) on K. */
void aic_ucp_carrier_Q(acb_mat_t Q, const aic_ucp_kraus *phi, slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H;
    assert(acb_mat_nrows(Q) == dk && acb_mat_ncols(Q) == dk &&
           "aic_ucp_carrier_Q: Q must be dim_K x dim_K");
    acb_mat_t Kd, term;
    acb_mat_init(Kd, dh, dk);   /* K_a^dag : K -> H */
    acb_mat_init(term, dk, dk); /* K_a K_a^dag : K -> K */
    acb_mat_zero(Q);
    for (slong a = 0; a < phi->r; a++) {
        acb_mat_conjugate_transpose(Kd, phi->K[a]);
        acb_mat_mul(term, phi->K[a], Kd, prec);
        acb_mat_add(Q, Q, term, prec);
    }
    acb_mat_clear(term);
    acb_mat_clear(Kd);
}

/* cor_carrier (.tex:1731) — ||(X (x) 1_F) V||_op = ||stack_a(X K_a)||_op. */
void aic_ucp_carrier_annihilate_defect(arb_t out, const aic_ucp_kraus *phi,
                                       const acb_mat_t X, slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H, r = phi->r;
    assert(acb_mat_nrows(X) == dk && acb_mat_ncols(X) == dk &&
           "aic_ucp_carrier_annihilate_defect: X must be dim_K x dim_K");
    acb_mat_t stack, blk;
    acb_mat_init(stack, dk * r, dh); /* rows: a-block of (X K_a) */
    acb_mat_init(blk, dk, dh);
    for (slong a = 0; a < r; a++) {
        acb_mat_mul(blk, X, phi->K[a], prec); /* X K_a : H -> K */
        for (slong i = 0; i < dk; i++)
            for (slong j = 0; j < dh; j++)
                acb_set(acb_mat_entry(stack, a * dk + i, j),
                        acb_mat_entry(blk, i, j));
    }
    aic_mat_opnorm(out, stack, prec);
    acb_mat_clear(blk);
    acb_mat_clear(stack);
}

/* PhiX_M (.tex:1740) — ||Phi(X) - Phi(Pi_M X Pi_M)||_op. */
void aic_ucp_phiX_M_defect(arb_t out, const aic_ucp_kraus *phi,
                           const acb_mat_t X, const acb_mat_t Pi_M, slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H;
    assert(acb_mat_nrows(X) == dk && acb_mat_ncols(X) == dk &&
           "aic_ucp_phiX_M_defect: X must be dim_K x dim_K");
    assert(acb_mat_nrows(Pi_M) == dk && acb_mat_ncols(Pi_M) == dk &&
           "aic_ucp_phiX_M_defect: Pi_M must be dim_K x dim_K");
    acb_mat_t PiX, comp, lhs, rhs, diff;
    acb_mat_init(PiX, dk, dk);
    acb_mat_init(comp, dk, dk); /* Pi_M X Pi_M */
    acb_mat_init(lhs, dh, dh);
    acb_mat_init(rhs, dh, dh);
    acb_mat_init(diff, dh, dh);
    acb_mat_mul(PiX, Pi_M, X, prec);
    acb_mat_mul(comp, PiX, Pi_M, prec);
    aic_ucp_apply(lhs, phi, X, prec);
    aic_ucp_apply(rhs, phi, comp, prec);
    acb_mat_sub(diff, lhs, rhs, prec);
    aic_mat_opnorm(out, diff, prec);
    acb_mat_clear(diff);
    acb_mat_clear(rhs);
    acb_mat_clear(lhs);
    acb_mat_clear(comp);
    acb_mat_clear(PiX);
}
