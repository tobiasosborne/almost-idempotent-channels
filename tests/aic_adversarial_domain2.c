/* aic_adversarial_domain2.c — adversarial CHANNEL generator: an eta-idempotent
 * UCP map whose defect Phi^2-Phi is EXACTLY RANK-1 (superoperator sense) and
 * CONCENTRATED on a NEAR-DEGENERATE subspace (knob gap_sub) — family 2B (bead
 * aic-cxo; docs/adversarial/domain.md:249-288). The stressor for the
 * O(sqrt eta) Phi_assoc1 cancellation (the W quantity, tex:2385-2422) that
 * powers the proof of th_almost_idemp (tex:2192-2204). Kept in its own file
 * (Rule 10; aic_adversarial_domain.c is already ~220 LOC). See aic_adversarial.h
 * for the corpus rationale (CLAUDE.md adversarial-first directive).
 *
 * ───────────────────────────────────────────────────────────────────────────
 * WHY 2B IS EVIL (domain.md:263-273, tex:2385-2422). The proof of
 * th_almost_idemp (tex:2192-2204: the equations Phi_assoc1/Phi_assoc2) hinges
 * on the intermediate quantity W (tex:2385-2422), a Choi diagram with two
 * (1-Pi) insertions (Pi = VV^dag, tex:2271). Expanding both rectangles as
 * 1-VV^dag, the O(1) Pi-terms CANCEL, leaving ||W|| = O(eta) by the
 * O(sqrt eta) bound of tex:2296. When the defect Phi^2-Phi is concentrated on a
 * NEAR-DEGENERATE subspace (eigenvalue gap gap_sub -> 0), each surviving term is
 * O(gap_sub) and the O(eta) result emerges only via a cancellation of two
 * O(1) quantities — catastrophic in double when gap_sub << eta, but TRACKED by
 * the certified arb ball. This generator supplies exactly that concentrated,
 * near-degenerate, rank-1 defect.
 *
 * ───────────────────────────────────────────────────────────────────────────
 * THE LITERAL domain.md FORMULA IS NOT UCP (a load-bearing FINDINGS §C18 result,
 * verified, NOT papered over). domain.md:255-261 proposes
 *   Phi = Phi_idemp + eta*delta_channel,
 *   delta_channel(X) = |v><v| Tr(rho_sub X) - Phi_idemp(|v><v| Tr(rho_sub X)).
 * Writing P = Phi_idemp (idempotent self-adjoint superop), D(X) = |v><v|Tr(rho X),
 * this is Phi = P + eta*(1-P)D. Its Choi has a CERTIFIED-NEGATIVE minimum
 * eigenvalue ~ -eta*c/2 < 0 for EVERY eta>0 and EVERY (v, rho_sub) tried (probe
 * grid eta in {0.20..1e-4}, gap_sub in {0.5..1e-8}: minEig(Choi) ranged
 * -2.7e-2..-1.7e-3, all NON-CP) — the subtracted CP piece -eta*PD always pushes
 * Phi out of the PSD cone. So the literal generator is NOT a channel; shipping it
 * as one would violate Rule 4. We instead build a GENUINELY UCP map with the SAME
 * named properties (rank-1 superoperator defect, concentrated on a near-degenerate
 * subspace, magnitude eta), via the measure-prepare route below.
 *
 * ───────────────────────────────────────────────────────────────────────────
 * THE UCP CONSTRUCTION (a measure-prepare map; the paper's OWN 1B route extended
 * to a 3D measured vector). On B(C^d), d>=3, in the OBSERVABLE / Heisenberg
 * convention Phi(X) = sum_a K_a^dag X K_a (aic_ucp.h):
 *   v = (sqrt(1-tail), sqrt(w1), sqrt(w2), 0, ..., 0)   (unit, supported on 0,1,2)
 *   K_0 = |v><0|  (d x d, column 0 = v),   K_j = |j><j|  (j = 1 .. d-1).
 * UNITAL: sum_a K_a^dag K_a = |0><0|<v|v> + sum_{j>=1}|j><j| = P_0 + (I-P_0) = I.
 * CP: every K_a is rank-1 (entanglement-breaking) => CP trivially. So Phi is UCP
 * for ALL (tail, w1, w2) — NO CP boundary (probe grid: minEig(Choi) ~ -1e-78,
 * certified PSD, every knob point). This is exactly the family-1B Kraus pattern
 * (aic_adv_chan_cb_op_gap), with the measured vector spread over THREE basis
 * directions instead of two so the defect's input functional rho_sub lives on a
 * 3D subspace whose near-degeneracy we tune.
 *
 * THE DEFECT IS EXACTLY RANK-1, OUTPUT P_0 (derived; the named property).
 * For a measure-prepare map Phi(X) = sum_j P_j Tr(gamma_j X) with P_j = |j><j|
 * and only gamma_0 = |v><v| tilted (gamma_j = P_j for j>=1), one computes
 *   Phi^2 - Phi = P_0 (x) <rho_sub, .>,   [superoperator outer product, RANK 1]
 *   rho_sub = -(sum_{k>=1}|v_k|^2) |v><v| + sum_{k>=1} |v_k|^2 |k><k|
 *           = -tail |v><v| + w1 |1><1| + w2 |2><2|,   tail = w1 + w2 = 1 - |v_0|^2.
 * (Verified: the n^2 x n^2 superoperator of Phi^2-Phi has sval[1] = 0.0 EXACTLY
 * across the whole knob grid — superoperator rank 1.) The OUTPUT direction P_0
 * is FIXED; the INPUT functional rho_sub is the "density on the near-degenerate
 * subspace" of domain.md:251-253 (here a Hermitian functional, the defect's
 * left/input covector).
 *
 * THE CALIBRATION (knob = the TARGET defect magnitude eta; closed form +
 * root-find, mirroring aic_adv_chan_noncomm_boundary). We set the TOTAL tail
 * weight via the family-1B closed form (tex:378):
 *   tail = the lower root of  g(tail) = tail*sqrt(1-tail) = eta   on [0, 2/3]
 * (bisection; g monotone increasing 0 -> peak 0.3849 at tail=2/3). The split is
 * then w2 = tail*gap_sub, w1 = tail - w2 (so the SMALL eigen-weight of rho_sub is
 * tail*gap_sub).
 *   CAVEAT ON THE NORM (hostile-review correction, FINDINGS C18): the ACTUAL
 * cb-norm of this rank-1 defect is ||Phi^2-Phi||_cb = ||rho_sub||_op (P_0 has
 * unit op-norm), which equals tail*sqrt(1-tail) = eta ONLY in the symmetric
 * gap_sub -> 0 limit; at a finite split it is eta + O(gap_sub*eta) (measured
 * ~0.9% ABOVE eta at the mild gap_sub=0.5, vanishing to <1e-9 at the lethal
 * gap_sub=1e-8 corner). The certified eig-free bracket merely CONTAINS eta (it is
 * 2n-wide, cannot pin the cb-norm); the TIGHT self-test anchor sval[0] =
 * ||rho_sub||_F pins TAIL (the FROBENIUS magnitude of the rank-1 defect), NOT the
 * cb-norm directly. So "eta_cb = eta" is exact at the adversarial near-degenerate
 * corner where the value lives, and a <1% over-estimate at the mild knob. cb != op
 * != F (the CLAUDE.md callout) — keep them distinct.
 *
 * THE gap_sub KNOB (the near-degenerate subspace). gap_sub in (0, 1] is the
 * RELATIVE weight of the smaller measured direction |2> in rho_sub:
 *   w2 / tail = gap_sub.
 * As gap_sub -> 0, the |2> weight w2 = tail*gap_sub -> 0, so rho_sub collapses
 * toward a 2D operator on span{0,1} — the |2> eigen-direction becomes
 * NEAR-DEGENERATE with the zero (kernel) subspace, the realized small eigenvalue
 * (~ tail*gap_sub) approaching 0. That vanishing-eigenvalue subspace is the
 * "near-degenerate subspace" the W cancellation (tex:2385-2422) is stressed on.
 * The defect magnitude is set by tail; the cb-norm = ||rho_sub||_op = eta exactly
 * as gap_sub -> 0 and within O(gap_sub*eta) (<1%) at finite gap_sub — so gap_sub
 * tunes the SUBSPACE near-degeneracy with only a sub-percent effect on magnitude.
 *
 * eta -> 0 / gap_sub fixed REDUCTION (exact-idempotent oracle, cross-check
 * ladder #3). As eta -> 0, tail -> 0, v -> |0>, K_0 -> |0><0|, Phi -> complete
 * dephasing in the standard basis — EXACTLY idempotent (defect 0); the cb bracket
 * collapses to [0,0]. The self-test pins this.
 *
 * Determinism: closed-form bisection (fixed iteration count) + exact arb sqrt;
 * no RNG. Same (d, eta, gap_sub, prec) -> bit-identical channel.
 *
 * Measured anchors (prec=256, d=4), asserted by tests/test_adversarial.c:
 *   eta=0.20 gap_sub=1e-8: tail=0.22756, sval[0]=0.282843 sval[1]=0.0 (rank-1),
 *     realized small rho_sub weight w2 = 2.2756e-9, minEig(Choi) ~ -5e-156 (CP),
 *     cb bracket contains 0.20, unital defect ~ 1e-16.
 *   eta=0.05 gap_sub=0.5: tail=0.05134, sval[0]=0.061787 sval[1]=0.0,
 *     w2 = 2.567e-2, cb bracket contains 0.05.
 *   eta=0 (any gap_sub): bracket [0,0] (exact-idempotent dephasing).
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_ucp.h"

/* g(tail) = tail*sqrt(1-tail), the family-1B closed-form cb defect (tex:378), a
 * function of the TOTAL tail weight only. On [0, 2/3] g is continuous and
 * strictly increasing from g(0)=0 to the peak g(2/3) = 0.3849001794597505. */
static double conc_g(double tail) { return tail * sqrt(1.0 - tail); }

/* Lower-root calibration: bisection for tail in [0, 2/3] with g(tail) = target.
 * Caller has asserted 0 <= target <= g(2/3); g monotone => bisection converges.
 * 200 iterations drive the bracket below 2^-200*(2/3), under the arb prec floor. */
static double conc_calibrate_tail(double target)
{
    double a = 0.0, b = 2.0 / 3.0; /* rising branch [0, argmax] */
    for (int it = 0; it < 200; it++) {
        double mid = 0.5 * (a + b);
        if (conc_g(mid) < target)
            a = mid;
        else
            b = mid;
    }
    return 0.5 * (a + b);
}

/* fam2B — concentrated rank-1 defect on a near-degenerate subspace
 * (domain.md:249-288, tex:2192-2204/2296/2385-2422). The measure-prepare map
 * K_0 = |v><0|, K_j = |j><j| (j>=1), with v = (sqrt(1-tail), sqrt(w1),
 * sqrt(w2), 0...) supported on span{0,1,2}. tail is the lower root of
 * tail*sqrt(1-tail) = eta (so eta_cb = ||Phi^2-Phi||_cb = eta), and the split is
 * w2 = tail*gap_sub, w1 = tail-w2 (so the near-degenerate eigen-direction |2> of
 * rho_sub carries the small weight tail*gap_sub). UCP for ALL knobs (1B Kraus
 * pattern), defect EXACTLY superoperator-rank-1 (output P_0, input rho_sub).
 * Asserts d>=3 (need three measured directions for the near-degenerate triple),
 * 0 < eta < g(2/3)=0.3849 (the eta<1/4-ish hypothesis; a target above the peak
 * has NO calibration root — fail loud, Rule 4), 0 < gap_sub <= 1. `out` is
 * aic_ucp_kraus_init'd HERE (caller clears with aic_ucp_kraus_clear). See the
 * file docstring for the full derivation, the non-CP-of-the-literal-formula
 * finding (FINDINGS §C18), the rank-1 proof, and the calibration. */
void aic_adv_chan_conc_defect(aic_ucp_kraus *out, slong d, double eta,
                              double gap_sub, slong prec)
{
    assert(d >= 3 && "aic_adv_chan_conc_defect: need d >= 3 (three measured "
                     "directions 0,1,2 for the near-degenerate rho_sub triple)");
    assert(eta > 0.0 && eta <= conc_g(2.0 / 3.0) &&
           "aic_adv_chan_conc_defect: need 0 < eta <= 0.3849 (eta = target cb "
           "defect; above the peak tail*sqrt(1-tail) there is NO calibration "
           "root — the construction cannot reach that magnitude)");
    assert(gap_sub > 0.0 && gap_sub <= 1.0 &&
           "aic_adv_chan_conc_defect: need 0 < gap_sub <= 1 (relative weight of "
           "the near-degenerate measured direction |2>)");

    /* Calibrate tail so ||Phi^2-Phi||_cb = tail*sqrt(1-tail) = eta (lower root,
     * the well-conditioned branch). Split the tail by gap_sub. */
    double tail = conc_calibrate_tail(eta);
    double w2 = tail * gap_sub; /* small eigen-weight of rho_sub: the near-degen one */
    double w1 = tail - w2;      /* dominant tail direction |1> */

    /* dim_K = dim_H = d (self-map on B(C^d)); Kraus rank r = d (rank-1 ops). */
    aic_ucp_kraus_init(out, d, d, d);

    /* K_0 = |v><0|: column 0 holds v = (sqrt(1-tail), sqrt(w1), sqrt(w2), 0...).
     * Certified arb sqrt of the (nonneg) weights. */
    arb_t s;
    arb_init(s);
    arb_set_d(s, 1.0 - tail);
    arb_sqrt(s, s, prec); /* v_0 = sqrt(1-tail) */
    arb_set(acb_realref(acb_mat_entry(out->K[0], 0, 0)), s);
    arb_set_d(s, w1);
    arb_sqrt(s, s, prec); /* v_1 = sqrt(w1) */
    arb_set(acb_realref(acb_mat_entry(out->K[0], 1, 0)), s);
    arb_set_d(s, w2);
    arb_sqrt(s, s, prec); /* v_2 = sqrt(w2) */
    arb_set(acb_realref(acb_mat_entry(out->K[0], 2, 0)), s);
    arb_clear(s);

    /* K_j = |j><j| for j = 1 .. d-1 (exact dephasing blocks; the gamma_j = P_j
     * outcomes that keep the map idempotent away from the tilted gamma_0). */
    for (slong j = 1; j < d; j++)
        acb_set_si(acb_mat_entry(out->K[j], j, j), 1);
}
