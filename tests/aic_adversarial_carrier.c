/* aic_adversarial_carrier.c — adversarial CHANNEL generator targeting the
 * CARRIER computation (bead aic-dbo.2, Increment 2 tranche 2). The shortlist's
 * "most likely to reveal a silent wrong-rank answer": a UCP map whose carrier
 * operator Q has one eigenvalue at `gap` (near zero) while the rest are O(1), so
 * the carrier dimension NEARLY drops as gap -> 0. Catalogued at
 * docs/adversarial/domain.md:127-159 (family 1C, "near-degenerate carrier —
 * carrier rank drop"). Sibling of aic_adversarial_domain.c (the cb-op-gap /
 * depol / unital-defect channel families); kept in its own file (Rule 10).
 *
 * WHY THIS FAMILY. lem_carrier (.tex:1724) and cor_carrier (.tex:1731) define the
 * carrier M of a UCP map Phi: B(K) -> B(H) as the support (range) of the carrier
 * operator. Our certified-rank routine aic_ucp_carrier_rank (bead aic-4td) must
 * certify dim M = rank Q for Q-eigenvalues above its threshold and either certify
 * the drop or FAIL LOUD near the boundary (the aic_mat_certified_rank
 * certify-or-straddle-abort contract, Rule 4) — never silently return the wrong
 * rank. In idemp_structure the sub-block structure depends on this rank EXACTLY
 * (domain.md:147, .tex:1917); a silent off-by-one carrier dimension corrupts the
 * whole decomposition downstream. This generator is the instance that drives the
 * carrier rank to the d-1 / d boundary.
 *
 * THE CARRIER-API CONVENTION (pinned from include/aic/aic_ucp.h + the existing
 * callers test_eigmult.c:T3 / test_eigvec.c:S6, NOT guessed — Rule 2). A UCP map
 * is r Kraus operators K_a: H -> K (each dim_K x dim_H), acting in the OBSERVABLE
 * convention Phi(X) = sum_a K_a^dag X K_a. The carrier operator that
 * aic_ucp_carrier_Q / aic_ucp_carrier_rank inspect is
 *     Q = sum_a K_a K_a^dag        (dim_K x dim_K, Hermitian PSD)   [.tex:1724]
 * i.e. the OUTPUT (K) marginal Tr_F(V V^dag) of the Stinespring isometry — NOT
 * sum_a K_a^dag K_a (= 1_H by unitality, the wrong marginal; aic_ucp_carrier.c
 * documents that self-correction). The carrier is the OUTPUT support, range(Q) in
 * K. The certified rank counts Q-eigenvalues certified ABOVE
 *     thr = dim_K * 2^-52 * ||Q||_F        (DBL_EPSILON = 2^-52)
 * an arb ball matching aic_ucp_carrier_rank_latd (aic_ucp_carrier.c:140-147).
 *
 * THE CONSTRUCTION (single Hermitian diagonal Kraus). A self-map on B(C^d),
 * dim_K = dim_H = d, with ONE Kraus operator
 *     K_0 = diag(1, 1, ..., 1, sqrt(gap))    on C^d,   gap in (0, 1].
 * One Kraus op is CP trivially. Then the carrier operator is
 *     Q = K_0 K_0^dag = diag(1, 1, ..., 1, gap),
 * with spectrum {1 (x)(d-1), gap}: the last OUTPUT dimension is barely populated
 * (amplitude sqrt(gap)), so the carrier rank is d for gap above thr and NEARLY
 * drops to d-1 as gap -> 0. The smallest carrier eigenvalue is EXACTLY gap (the
 * K-marginal squares the sqrt(gap) amplitude back to gap), the named adversarial
 * quantity the self-test pins. This is NOT unital for gap != 1 (sum_a K_a^dag K_a
 * = diag(1,...,1,gap) != 1_H) — deliberately: the task allows a non-unital
 * adversarial instance, and a diagonal carrier is what lets the SMALL eigenvalue
 * land cleanly in Q = sum K K^dag (a unital padding would refill the d-th output
 * dimension and destroy the rank near-drop). gap = 1 gives K_0 = 1_d, Q = 1_d,
 * the EXACT full-rank d carrier (the eta=0-style oracle).
 *
 * NEAR-DROP BEHAVIOR IS *CERTIFY d-1*, NOT FAIL-LOUD (the measured finding, Rule
 * 2 / Rule 6 — discovered, not assumed). Because K_0 (hence Q) is DIAGONAL with
 * an exactly-representable small entry, the certified eigensolver returns a POINT
 * ball for the gap eigenvalue (zero radius), and ||Q||_F is tight, so thr is a
 * tight ball too. The decision arb_gt(gap, thr) / arb_lt(gap, thr) is therefore
 * NEVER undecided: the routine CERTIFIES rank d when gap > thr and CERTIFIES rank
 * d-1 when gap < thr (measured d=2 transition at thr ~= 4.44e-16: gap=4.44e-16 ->
 * rank 1, gap=4.5e-16 -> rank 2; gap=0 -> rank 1 exactly). The STRADDLE/fail-loud
 * path needs an eigenvalue ball WITH radius overlapping thr — produced by a
 * DENSIFIED (non-diagonal) rank-deficient carrier, already exercised by
 * test_eigvec.c:S6b (the FINDINGS §D7 zero-cluster straddle). This clean diagonal
 * family deliberately stays on the CERTIFIABLE side so the self-test can assert
 * the EXACT rank and the smallest-eigenvalue == gap pin; the gap sweep at d=2
 * therefore reports rank d for every gap in the catalogue sweep {0.5, 1e-4, 1e-8,
 * ...} (all >> thr) and rank d-1 only at the gap=0 (exact-drop) oracle.
 *
 * Determinism: exact closed form (arb sqrt of gap). No RNG. `out` is
 * aic_ucp_kraus_init'd HERE (dim_K = dim_H = d, r = 1; caller clears with
 * aic_ucp_kraus_clear), matching the aic_adv_chan_* convention.
 *
 * Measured anchors (prec=256, d=2), asserted by tests/test_adversarial.c:
 *   gap=0.5  : Q=diag(1,0.5),  smallest carrier eig = 0.5,  certified rank 2
 *   gap=1e-3 : Q=diag(1,1e-3), smallest carrier eig = 1e-3, certified rank 2
 *   gap=1e-6 : Q=diag(1,1e-6), smallest carrier eig = 1e-6, certified rank 2
 *   gap=0    : Q=diag(1,0),    smallest carrier eig = 0,     certified rank 1 (drop)
 *   gap=1    : Q=diag(1,1),    smallest carrier eig = 1,     full rank 2 (oracle)
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_ucp.h"

/* fam1C — near-degenerate carrier (domain.md:127-159, lem_carrier .tex:1724 /
 * cor_carrier .tex:1731). Single Hermitian diagonal Kraus
 *   K_0 = diag(1, ..., 1, sqrt(gap))   on C^d,
 * so the carrier operator Q = sum_a K_a K_a^dag = K_0 K_0^dag = diag(1,...,1,gap)
 * has spectrum {1 (x)(d-1), gap}; the carrier rank NEARLY drops to d-1 as gap->0.
 * KNOB gap in (0, 1]; asserts gap > 0 and d >= 2 (Rule 4, fail loud). `out` is
 * aic_ucp_kraus_init'd HERE (dim_K=dim_H=d, r=1; caller aic_ucp_kraus_clears it).
 * NOT unital for gap != 1 (deliberate; see file docstring). */
void aic_adv_chan_carrier_dropout(aic_ucp_kraus *out, slong d, double gap,
                                  slong prec)
{
    assert(d >= 2 && "aic_adv_chan_carrier_dropout: need d >= 2");
    assert(gap > 0.0 && gap <= 1.0 &&
           "aic_adv_chan_carrier_dropout: need 0 < gap <= 1");

    /* dim_K = dim_H = d (self-map on B(C^d)); single Hermitian Kraus op (r=1). */
    aic_ucp_kraus_init(out, d, d, 1);

    /* K_0 = diag(1, ..., 1, sqrt(gap)) — init zeroes everything; set the diagonal.
     * Q = K_0 K_0^dag = diag(1, ..., 1, gap) gets the small eigenvalue gap. */
    for (slong j = 0; j < d - 1; j++)
        acb_set_si(acb_mat_entry(out->K[0], j, j), 1);
    arb_t s;
    arb_init(s);
    arb_set_d(s, gap);
    arb_sqrt(s, s, prec); /* sqrt(gap) as a certified ball */
    arb_set(acb_realref(acb_mat_entry(out->K[0], d - 1, d - 1)), s);
    arb_clear(s);
}
