/* aic_channels.h — PUBLIC physical-channel constructors (bead aic-7hg, module
 * "channels"). Kills application friction F7: until now every standard channel
 * (dephasing, depolarizing, Pauli, twirl, conditional expectation) had to be
 * hand-rolled from raw Kraus operators at the call site. This module exports the
 * standard library so applications + the Julia/Python bindings build channels by
 * name. Each constructor INITIALISES the supplied aic_ucp_kraus (calling
 * aic_ucp_kraus_init internally); the CALLER clears it with aic_ucp_kraus_clear.
 *
 * KRAUS CONVENTION (load-bearing; matches include/aic/aic_ucp.h exactly).
 *   These are OBSERVABLE / Heisenberg (UCP) maps Phi : B(C^d) -> B(C^d) stored as
 *   r Kraus operators K_a (each d x d) acting as
 *
 *       Phi(X) = sum_a K_a^dag X K_a,        unital iff  sum_a K_a^dag K_a = 1.
 *
 *   This is NOT the Schrodinger/state dual T = Phi^*(X) = sum_a K_a X K_a^dag
 *   (CLAUDE.md "UCP vs CPTP" callout). Getting the dag-order backwards silently
 *   transposes the channel. The two are equal only for Hermitian Kraus (e.g. the
 *   diagonal dephasing / Pauli builders, where K_a = K_a^dag); the unitary-mixture
 *   builders below carry K_j = sqrt(p_j) U_j^dag PRECISELY so that the OBSERVABLE
 *   action is the textbook Phi(X) = sum_j p_j U_j X U_j^dag.
 *
 * INPUT VALIDATION (CLAUDE.md Rule 4, fail loud). Every constructor asserts its
 * preconditions and ABORTS (flint_abort) on violation, never silently returns a
 * malformed channel:
 *   - probabilities: each >= -tol and sum within tol of 1  (tol ~ 1e-12);
 *   - unitaries: ||U^dag U - 1|| <= tol  (genuinely unitary);
 *   - block dims: each >= 1 and sum == d.
 * Validation runs on the double-precision midpoints at the prec given; it is a
 * caller-error guard, not a certified bound.
 *
 * IDEMPOTENCE PROPERTIES (the cleanest cross-check, CLAUDE.md Rule 6). Several of
 * these are EXACTLY idempotent (eta=0, Phi^2 = Phi to machine) and serve as the
 * eta=0 oracle: dephasing, depolarizing at p=1, cond_expectation, and a twirl over
 * a CLOSED group. depolarizing at general p in (0,1) is almost-idempotent. The
 * per-constructor docstrings state which.
 */
#ifndef AIC_CHANNELS_H
#define AIC_CHANNELS_H

#include <flint/acb_mat.h>

#include "aic/aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Dephasing on C^d: Phi(X) = diag(X) (kills off-diagonals). Kraus K_i =
 * |e_i><e_i| (d Hermitian ops), so Phi(X) = sum_i |e_i><e_i| X |e_i><e_i|.
 * Unital, trace-preserving, EXACTLY idempotent (eta=0). Img Phi = diagonal
 * matrices (commutative), so the associated algebra is C^d. Matches the
 * test-fixture oracle make_dephasing. `phi` is init'd here (caller clears). */
void aic_channel_dephasing(aic_ucp_kraus *phi, slong d, slong prec);

/* Depolarizing on C^d (the UNITAL / observable convention):
 *     Phi_p(X) = (1 - p) X + p (tr(X)/d) 1_d,    p in [0,1].
 * This IS unital: Phi_p(1) = (1-p) 1 + p (d/d) 1 = 1. It is its own observable
 * dual (self-adjoint as a superoperator), so the same closed form holds in either
 * picture. Kraus form: K_0 = sqrt(1 - p + p/d) 1_d for the identity weight plus
 * the trace-replace block { sqrt(p/d) (1/sqrt d) ... } — see the implementation;
 * net (1 + d^2) ops. p=0 is the identity; p=1 is the trace-replace conditional
 * expectation onto C*1 (EXACTLY idempotent); general p in (0,1) is
 * almost-idempotent with Phi_p^2 = Phi_{p(2-p)}. `phi` init'd here. */
void aic_channel_depolarizing(aic_ucp_kraus *phi, slong d, double p, slong prec);

/* Qubit (d=2) Pauli channel: Phi(X) = sum_{k=0}^{3} p_k sigma_k X sigma_k with
 * (sigma_0,sigma_1,sigma_2,sigma_3) = (1, X, Y, Z). The Pauli matrices are
 * Hermitian so K_k = sqrt(p_k) sigma_k (observable = dual here). Unital and
 * trace-preserving. probs[4] must be >= 0 and sum to 1 (asserted). probs =
 * (1,0,0,0) reduces to the identity; (1/2,0,0,1/2) reduces to qubit dephasing.
 * `phi` init'd here (caller clears). */
void aic_channel_pauli(aic_ucp_kraus *phi, const double probs[4], slong prec);

/* Unital mixed-unitary channel on C^d:
 *     Phi(X) = sum_{j=0}^{m-1} p_j U_j X U_j^dag.
 * The general unital-channel constructor the twirl reuses. Kraus operators are
 * K_j = sqrt(p_j) U_j^dag (so the OBSERVABLE action is the stated sum; see the
 * header convention note). Validation: each U_j unitary (||U_j^dag U_j - 1|| small),
 * probs >= 0 and sum to 1 (asserted, fail loud). Unital because sum_j p_j U_j U_j^dag
 * = sum_j p_j 1 = 1. Generally NOT idempotent. `U` is an array of m initialised
 * d x d acb_mat_t; `probs` is a length-m double array. `phi` init'd here. */
void aic_channel_mix_unitaries(aic_ucp_kraus *phi, const acb_mat_t *U,
                               const double *probs, slong m, slong d, slong prec);

/* Group twirl on C^d: Phi(X) = (1/m) sum_{g=0}^{m-1} U_g X U_g^dag — the uniform
 * mixed-unitary channel over the provided representation matrices {U_g}. This is
 * aic_channel_mix_unitaries with weights 1/m. IDEMPOTENCE: if {U_g} form a GROUP
 * (closed under multiplication, m = |G|, a unitary representation) the twirl is
 * the orthogonal projection onto the commutant of the representation and is
 * EXACTLY idempotent (Phi^2 = Phi). If {U_g} is NOT a closed group the map is
 * still a valid unital channel but generally not idempotent — the constructor does
 * NOT verify closure (that is the caller's contract; documented here). `U` is m
 * initialised d x d unitaries. `phi` init'd here (caller clears). */
void aic_channel_group_twirl(aic_ucp_kraus *phi, const acb_mat_t *U, slong m,
                             slong d, slong prec);

/* Conditional expectation onto a block-diagonal subalgebra (pinching map) on
 * C^d = (+)_i C^{block_dims[i]}: Phi(X) = sum_i P_i X P_i, where P_i is the
 * orthogonal projector onto block i (the contiguous index range for that block).
 * Kraus K_i = P_i (nblocks Hermitian projectors with sum_i P_i = 1_d). Unital,
 * trace-preserving, EXACTLY idempotent (eta=0): it is the conditional expectation
 * onto (+)_i M_{block_dims[i]}. Generalises the test-fixture make_block_cond_exp
 * (which is the two-block case). `block_dims` is a length-nblocks array of block
 * sizes (each >= 1; sum must equal d, asserted). `phi` init'd here. */
void aic_channel_cond_expectation(aic_ucp_kraus *phi, const slong *block_dims,
                                  slong nblocks, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CHANNELS_H */
