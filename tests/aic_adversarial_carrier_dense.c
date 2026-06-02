/* aic_adversarial_carrier_dense.c — the DENSIFIED near-degenerate carrier CHANNEL
 * generator (bead aic-v5f), the hostile-review follow-up to the diagonal 1C family
 * (aic_adversarial_carrier.c, aic_adv_chan_carrier_dropout). It closes the two gaps
 * the diagonal 1C deliberately leaves open, catalogued at docs/adversarial/domain.md
 * 1C and lem_carrier (.tex:1724) / cor_carrier (.tex:1731). Kept in its own file
 * (Rule 10), sibling of aic_adversarial_carrier.c.
 *
 * WHY A SECOND CARRIER GENERATOR (the two gaps the diagonal 1C cannot reach).
 *  (i) STRADDLE / fail-loud. The diagonal 1C carrier Q = diag(1,...,1,gap) is
 *      coordinate-aligned with an exactly-representable small entry, so the
 *      certified eigensolver returns a POINT ball (zero radius) for the gap
 *      eigenvalue: aic_ucp_carrier_rank ALWAYS decides (arb_gt / arb_lt), it
 *      never exercises the certify-OR-straddle-abort contract of
 *      aic_mat_certified_rank (aic_adversarial_carrier.c:50-64, FINDINGS §D7). A
 *      DENSIFIED carrier Q = U diag(1,...,1,gap) U^dag (U a fixed non-trivial
 *      rational unitary) has a NON-axis-aligned eigenvector, so FLINT Rump's
 *      enclosure of the small/zero cluster carries a finite radius (~the
 *      Krawczyk conditioning floor). MEASURED (probe, prec=53, d=3): the
 *      small-cluster ball radius is ~1.5e-15 while thr = dim_K*2^-52*||Q||_F is
 *      ~9.4e-16, so for gap below ~1e-15 the ball STRADDLES thr and
 *      aic_ucp_carrier_rank fail-loud-ABORTS ("STRADDLES") — the contract the
 *      diagonal family could not trigger. At prec>=64 the ball radius drops to
 *      ~7e-19 (prec=64) / ~4e-38 (prec=128) and the routine certifies cleanly
 *      (rank d for gap>>thr, rank d-1 for the near-zero gap). This is the exact
 *      prec floor of FINDINGS §D7 (zero ball ~2e-14 STRADDLES keep_thr ~1e-15 at
 *      prec=53; clears at prec>=64), here at the channel-carrier level.
 *  (ii) CONVENTION-SENSITIVITY. The diagonal 1C Kraus K_0 is HERMITIAN, so
 *      K_0 K_0^dag = K_0^dag K_0 = K_0^2: the two carrier conventions coincide
 *      and a routine that wrongly computes sum K^dag K (= the H-marginal 1_H by
 *      unitality, the WRONG marginal — aic_ucp_carrier.c:18) instead of the
 *      correct sum K K^dag (the K-marginal carrier, .tex:1724) is UNDETECTABLE.
 *      A NON-NORMAL Kraus K = U diag(...) V^dag with U != V makes the marginals
 *      DIFFER: K K^dag = U diag U^dag (the carrier) vs K^dag K = V diag V^dag,
 *      so ||sum K K^dag - sum K^dag K||_op > 0 (MEASURED: 0.223 at d=3 gap=0.5,
 *      up to 0.605 at d=4 gap->0). A convention bug in any carrier-side routine
 *      is then catchable.
 *
 * THE CONSTRUCTION (single non-normal Kraus, two fixed rational unitaries).
 * A self-map on B(C^d), dim_K = dim_H = d, with ONE Kraus operator
 *     K_0 = U diag(1, ..., 1, sqrt(gap)) V^dag       on C^d,   gap in (0, 1],
 * where U and V are DIFFERENT fixed rational-Givens chain unitaries (no RNG):
 *     U : cos=3/5,  sin=4/5   (the aic_mat_dense_unitary densifier triple)
 *     V : cos=5/13, sin=12/13 (a second Pythagorean triple => U != V)
 * each the product over ALL planes (a,b), a<b, of its Givens rotation (so both
 * spread every coordinate; mirror of src/aic_mat_densify.c, built SELF-CONTAINED
 * here because the generator needs TWO distinct unitaries, not the single shared
 * densifier). One Kraus op is CP trivially. Then:
 *     Q := sum_a K_a K_a^dag = K_0 K_0^dag = U diag(1,...,1,gap) U^dag,
 * a NON-DIAGONAL Hermitian PSD carrier with spectrum {1 (x)(d-1), gap} (identical
 * to the diagonal 1C SPECTRUM, but the eigenvectors are now the dense columns of
 * U). The carrier rank is d for gap above thr and NEARLY drops to d-1 as gap->0;
 * unlike the diagonal 1C, the densified small-cluster ball straddles thr at
 * prec=53. The OTHER convention sum_a K_a^dag K_a = V diag(1,...,1,gap) V^dag is a
 * DIFFERENT dense matrix, the convention-sensitivity above.
 *
 * VALIDITY: a carrier-only fixture, NOT unital for gap != 1 (sum K^dag K =
 * V diag(1,...,1,gap) V^dag != 1_H), DELIBERATELY — exactly like the diagonal 1C
 * (aic_adversarial_carrier.c:43-48): a unital padding would refill the d-th output
 * dimension and destroy the rank near-drop. It is CP (single Kraus) and a valid
 * UCP-shaped aic_ucp_kraus; it punishes the CARRIER routines (aic_ucp_carrier_Q /
 * aic_ucp_carrier_rank), not the full almost-idempotent pipeline. gap = 1 gives
 * K_0 = U V^dag (a unitary), Q = 1_d, the EXACT full-rank d carrier oracle, AND it
 * IS unital then (K_0^dag K_0 = V V^dag = 1_d).
 *
 * Determinism: exact closed form (rational Givens chains + arb sqrt of gap). No
 * RNG. `out` is aic_ucp_kraus_init'd HERE (dim_K = dim_H = d, r = 1; caller clears
 * with aic_ucp_kraus_clear), matching the aic_adv_chan_* convention.
 *
 * Measured anchors (prec as noted, asserted by tests/test_adversarial.c
 * test_fam1c_carrier_dense):
 *   gap=0.5, prec=128 : Q dense, smallest carrier eig = 0.5, certified rank d
 *   gap=1e-16, prec=53 : small-cluster ball straddles thr -> carrier_rank ABORTS
 *   gap=1e-16, prec=128: certified rank d-1 (the drop certifies at headroom prec)
 *   convention gap ||sum K K^dag - sum K^dag K||_op = 0.223 (d=3,gap=0.5) > 0
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_ucp.h"

/* A fixed rational-Givens CHAIN unitary on C^n: product over ALL planes (a,b),
 * a<b, of the rotation with the Pythagorean (cos=cn/cd, sin=sn/sd). Mirror of
 * src/aic_mat_densify.c's aic_mat_dense_unitary (which hard-codes 3/5,4/5); kept
 * local + parametrised because the non-normal Kraus needs TWO DISTINCT unitaries
 * U != V. cn^2/cd^2 + sn^2/sd^2 must = 1 (asserted by the caller's choice of
 * Pythagorean triples). Unitary far below the working prec (the chain of n(n-1)/2
 * exact-rational rotations; ||U U^dag - I||_F ~ n^2 * 2^-prec). */
static void carrier_dense_givens_chain(acb_mat_t U, slong n, long cn, long cd,
                                       long sn, long sd, slong prec)
{
    acb_t c, s, ns;
    acb_init(c);
    acb_init(s);
    acb_init(ns);
    acb_set_si(c, cn);
    acb_div_si(c, c, cd, prec);          /* cos = cn/cd */
    acb_set_si(s, sn);
    acb_div_si(s, s, sd, prec);          /* sin = sn/sd */
    acb_neg(ns, s);

    acb_mat_one(U);
    acb_mat_t G, T;
    acb_mat_init(G, n, n);
    acb_mat_init(T, n, n);
    for (slong a = 0; a < n; a++) {
        for (slong b = a + 1; b < n; b++) {  /* every plane (a,b), a<b */
            acb_mat_one(G);
            acb_set(acb_mat_entry(G, a, a), c);
            acb_set(acb_mat_entry(G, b, b), c);
            acb_set(acb_mat_entry(G, a, b), ns);
            acb_set(acb_mat_entry(G, b, a), s);
            acb_mat_mul(T, U, G, prec);
            acb_mat_set(U, T);
        }
    }
    acb_mat_clear(T);
    acb_mat_clear(G);
    acb_clear(ns);
    acb_clear(s);
    acb_clear(c);
}

/* fam1C-dense — DENSIFIED near-degenerate carrier (domain.md 1C; lem_carrier
 * .tex:1724 / cor_carrier .tex:1731; FINDINGS §D7 prec floor). Single NON-NORMAL
 * Kraus K_0 = U diag(1,...,1,sqrt(gap)) V^dag, U(3/5,4/5) != V(5/13,12/13) fixed
 * rational unitaries, so the carrier Q = sum K K^dag = U diag(1,...,1,gap) U^dag is
 * NON-DIAGONAL (spectrum {1(x)(d-1),gap}, straddles thr at prec=53 for gap->0) and
 * sum K^dag K = V diag(1,...,1,gap) V^dag DIFFERS (convention-sensitive). KNOB gap
 * in (0,1]; asserts gap>0 and d>=2 (Rule 4, fail loud). NOT unital for gap!=1
 * (deliberate carrier-only fixture; see file docstring). `out` is
 * aic_ucp_kraus_init'd HERE (dim_K=dim_H=d, r=1; caller aic_ucp_kraus_clears it). */
void aic_adv_chan_carrier_dropout_dense(aic_ucp_kraus *out, slong d, double gap,
                                        slong prec)
{
    assert(d >= 2 && "aic_adv_chan_carrier_dropout_dense: need d >= 2");
    assert(gap > 0.0 && gap <= 1.0 &&
           "aic_adv_chan_carrier_dropout_dense: need 0 < gap <= 1");

    /* dim_K = dim_H = d (self-map on B(C^d)); single non-normal Kraus op (r=1). */
    aic_ucp_kraus_init(out, d, d, 1);

    acb_mat_t U, V, Vd, D, T;
    acb_mat_init(U, d, d);
    acb_mat_init(V, d, d);
    acb_mat_init(Vd, d, d);
    acb_mat_init(D, d, d);
    acb_mat_init(T, d, d);

    /* Two DISTINCT fixed rational unitaries (Pythagorean triples 3-4-5, 5-12-13).
     * U != V => K_0 = U diag V^dag is non-normal => K K^dag != K^dag K. */
    carrier_dense_givens_chain(U, d, 3, 5, 4, 5, prec);
    carrier_dense_givens_chain(V, d, 5, 13, 12, 13, prec);
    acb_mat_conjugate_transpose(Vd, V);

    /* D = diag(1,...,1,sqrt(gap)); init zeroed D, set the diagonal. The singular
     * values of K_0 are {1 (x)(d-1), sqrt(gap)}, so Q = K_0 K_0^dag has spectrum
     * {1 (x)(d-1), gap}: the small carrier eigenvalue is EXACTLY gap. */
    acb_mat_zero(D);
    for (slong j = 0; j < d - 1; j++)
        acb_set_si(acb_mat_entry(D, j, j), 1);
    arb_t sg;
    arb_init(sg);
    arb_set_d(sg, gap);
    arb_sqrt(sg, sg, prec);              /* sqrt(gap) as a certified ball */
    arb_set(acb_realref(acb_mat_entry(D, d - 1, d - 1)), sg);
    arb_clear(sg);

    /* K_0 = U D V^dag (SVD form with U != V => non-normal). */
    acb_mat_mul(T, U, D, prec);
    acb_mat_mul(out->K[0], T, Vd, prec);

    acb_mat_clear(T);
    acb_mat_clear(D);
    acb_mat_clear(Vd);
    acb_mat_clear(V);
    acb_mat_clear(U);
}
