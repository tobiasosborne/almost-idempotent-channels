/* aic_adversarial_domain.c — adversarial CHANNEL generators (bead aic-dbo.2,
 * Increment 2). The FIRST channel/UCP-map family in the evil-instance corpus:
 * every prior generator (aic_adversarial.c, aic_adversarial_nla.c) produced a
 * bare acb_mat_t; this file produces an aic_ucp_kraus — a full UCP map — so the
 * corpus can punish the channel routines (ucp, cbnorm, assoc_ecsa, opspace), not
 * just the numerical-linear-algebra cores. See aic_adversarial.h for the corpus
 * rationale (CLAUDE.md adversarial-first directive, docs/adversarial/README.md);
 * THIS file's family is catalogued at docs/adversarial/domain.md:75-123 (family
 * 1B, "cb-norm vs operator-norm gap").
 *
 * WHY A CHANNEL GENERATOR. The single most load-bearing CLAUDE.md hallucination
 * callout is "cb-norm != operator norm": eta-idempotence is
 * ||Phi^2-Phi||_cb <= eta (tex:347, tex:354), a sup over ALL ampliations
 * 1_{M_n} (x) Phi, NOT the plain operator norm of the single-copy defect. A
 * routine that estimates the defect in operator norm underreports eta by up to
 * the ampliation factor n (the trap tex:484 warns about). This generator is the
 * canonical witness of that gap: the paper's OWN motivating measure-prepare
 * example (tex:366-388), whose defect is provably ||Phi^2-Phi||_cb =
 * eta*sqrt(1-eta) (tex:378) while the single-copy operator-norm defect is
 * strictly smaller. Feeding it to a cb-norm certifier and to an op-norm proxy
 * exposes the gap as a measured, certified number.
 *
 * tex:366-377 (verbatim), the example, on B(C^2) (extended to C^d here):
 *   Phi(X) = P_0 Tr(gamma_0 X) + P_1 Tr(gamma_1 X),
 *   P_0 = [[1,0],[0,0]], P_1 = [[0,0],[0,1]],
 *   gamma_0 = [[1-eta, sqrt(eta(1-eta))],[sqrt(eta(1-eta)), eta]],
 *   gamma_1 = [[0,0],[0,1]].
 * tex:378 (verbatim): "One can check that ||Phi^2-Phi||_cb = eta sqrt(1-eta)".
 *
 * KRAUS DERIVATION (the OBSERVABLE/Heisenberg convention Phi(X)=sum_a K_a^dag X K_a,
 * aic_ucp.h). Take the unit vector v = (sqrt(1-eta), sqrt(eta), 0, ..., 0) in C^d
 * (<v|v> = (1-eta)+eta = 1), so gamma_0 = |v><v|: indeed
 *   |v><v| = [[1-eta, sqrt((1-eta)eta)], [sqrt(eta(1-eta)), eta]] = gamma_0.
 * Set K_0 = |v><0| (d x d: column 0 equals v, the rest zero). Then for any X,
 *   (K_0^dag X K_0)[a,b] = sum_{p,q} conj(v[p]delta_{a0}) X[p,q] v[q]delta_{b0}
 *                        = delta_{a0}delta_{b0} <v|X|v>,
 * so K_0^dag X K_0 = <v|X|v> |0><0| = Tr(gamma_0 X) P_0. Likewise K_j = |j><j|
 * (j >= 1) gives K_j^dag X K_j = <j|X|j> |j><j| = Tr(gamma_j X) P_j with
 * gamma_j = |j><j| = P_j. So
 *   Phi(X) = Tr(gamma_0 X) P_0 + sum_{j>=1} Tr(P_j X) P_j,
 * which for d=2 is EXACTLY the paper's map (gamma_1 = |1><1| = P_1), and for d>2
 * pads with exactly-idempotent dephasing blocks |j><j| (j=2..d-1) that contribute
 * ZERO to the defect (Phi restricted to span{|j>} is the identity dephasing, an
 * exact idempotent). Kraus rank d, every K_a rank 1 -> entanglement-breaking.
 * UNITAL: sum_a K_a^dag K_a = |0><0|<v|v> + sum_{j>=1}|j><j| = I.
 *
 * eta=0 REDUCTION (the exact-idempotent oracle, CLAUDE.md cross-check ladder #3).
 * At eta=0, v = (1,0,...,0) = |0>, so K_0 = |0><0| and Phi(X) = sum_j <j|X|j> |j><j|
 * is the complete dephasing channel in the standard basis — EXACTLY idempotent
 * (Phi^2 = Phi), defect 0. The cb bracket collapses to [0,0]. The self-test pins
 * this.
 *
 * Determinism: exact closed form (arb sqrt of 1-eta / eta). No RNG.
 *
 * Measured anchors (prec=256, d=2), asserted by tests/test_adversarial.c:
 *   eta=0.05: eta*sqrt(1-eta)=0.048734, in bracket [0.034460, 0.137840]
 *   eta=0.20: eta*sqrt(1-eta)=0.178885, in bracket [0.126491, 0.505964]
 *   eta=0.00: bracket [0,0]  (exact idempotent)
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_channels.h"
#include "aic/aic_ucp.h"

/* fam1B — cb-norm vs operator-norm gap (domain.md:75-123, tex:366-388). The
 * paper's measure-prepare example on C^d in the OBSERVABLE convention
 * Phi(X)=sum_a K_a^dag X K_a:
 *   v = (sqrt(1-eta), sqrt(eta), 0, ..., 0)   (unit, <v|v>=1)
 *   K_0 = |v><0|  (column 0 = v),   K_j = |j><j|  (j = 1 .. d-1).
 * Defect ||Phi^2-Phi||_cb = eta*sqrt(1-eta) (tex:378), STRICTLY larger than the
 * single-copy operator-norm defect: the named adversarial property. KNOB eta in
 * [0,1] (the gap is maximized near eta~1/4; catalogue sweep {0.05,..,0.20}).
 * `out` is aic_ucp_kraus_init'd HERE (caller clears with aic_ucp_kraus_clear),
 * matching the aic_channels.h channel-constructor convention. Asserts d>=2 and
 * 0<=eta<=1 (Rule 4, fail loud). */
void aic_adv_chan_cb_op_gap(aic_ucp_kraus *out, slong d, double eta, slong prec)
{
    assert(d >= 2 && "aic_adv_chan_cb_op_gap: need d >= 2");
    assert(eta >= 0.0 && eta <= 1.0 &&
           "aic_adv_chan_cb_op_gap: need 0 <= eta <= 1");

    /* dim_K = dim_H = d (self-map on B(C^d)); Kraus rank r = d (rank-1 ops). */
    aic_ucp_kraus_init(out, d, d, d);

    /* K_0 = |v><0|: column 0 holds v = (sqrt(1-eta), sqrt(eta), 0,...,0). */
    arb_t s;
    arb_init(s);
    arb_set_d(s, 1.0 - eta);
    arb_sqrt(s, s, prec); /* sqrt(1-eta) as a certified ball */
    arb_set(acb_realref(acb_mat_entry(out->K[0], 0, 0)), s);
    arb_set_d(s, eta);
    arb_sqrt(s, s, prec); /* sqrt(eta) */
    arb_set(acb_realref(acb_mat_entry(out->K[0], 1, 0)), s);
    arb_clear(s);

    /* K_j = |j><j| for j = 1 .. d-1 (exact dephasing blocks; |1><1| reproduces
     * gamma_1 = P_1, the rest pad with zero-defect idempotent blocks). */
    for (slong j = 1; j < d; j++)
        acb_set_si(acb_mat_entry(out->K[j], j, j), 1);
}

/* fam2A — depolarizing eta->1/4 regularization boundary (domain.md:196-245,
 * "family 2A", tex:516-525 the theta(2Phi-1) basin rho(Phi^2-Phi) < 1/4). The
 * cleanest tunable instance sitting AT the regularization basin edge: the
 * standard depolarizing channel, whose eta-idempotence defect has an EXACTLY
 * computable structure.
 *
 * THE EXACT DEFECT (the value the self-test pins; cite this). The unital
 * depolarizing map (aic_channels.h) is
 *   Phi_p(X) = (1-p) X + p (tr X / d) 1_d,    p in [0,1].
 * Its square is Phi_p^2 = Phi_{p(2-p)} (aic_channels.h docstring), so the
 * idempotence defect is
 *   Phi_p^2 - Phi_p = (Phi_{p(2-p)} - Phi_p) = (p(2-p) - p) C0 = p(1-p) C,
 *   where  C(X) = (tr X / d) 1_d - X     (i.e. Phi_q = id - q C, so
 *   Phi_q - Phi_p = -(q-p) C, and q - p = p(2-p) - p = p(1-p)).
 * C is the trace-orthogonal projector defect: as a SUPEROPERATOR it has
 * eigenvalue 0 on X = 1_d (tr/d = 1, so C(1)=0) and eigenvalue -1 on each of the
 * d^2-1 TRACELESS directions (tr X = 0 => C(X) = -X). Hence
 *   - SPECTRAL RADIUS  rho(Phi_p^2 - Phi_p) = p(1-p),  d-INDEPENDENT,
 *     maximized = 1/4 at p = 1/2  (EXACTLY the theta(2Phi-1) basin edge rho<1/4).
 *   - cb-norm DEFECT  ||Phi_p^2 - Phi_p||_cb = p(1-p) ||C||_cb, ||C||_cb a fixed
 *     constant of d only, so the defect scales EXACTLY LINEARLY in p(1-p). The
 *     eig-free bracket inherits this exactly: ||J||_F = p(1-p) ||Choi(C)||_F,
 *     so lo = ||J||_F / n = (||Choi(C)||_F / n) p(1-p)  (measured at d=2:
 *     lo = (sqrt(3)/2) p(1-p) = 0.8660254 p(1-p), constant lo/[p(1-p)]).
 *
 * KNOB p in [0,1]; the basin boundary is approached as p -> 1/2 (rho = p(1-p)
 * -> 1/4). The catalogue 2A sweep targets eta near 1/4: at d=2, p=0.789 gives
 * rho ~= 0.167; p=1/2 gives rho = 1/4 the exact edge.
 *
 * eta=0 REDUCTION (exact-idempotent oracle, CLAUDE.md cross-check ladder #3).
 * p=0 is the identity (Phi^2=Phi trivially); p=1 is the trace-replace
 * conditional expectation onto C*1 (Phi_1^2 = Phi_{1*(2-1)} = Phi_1, exactly
 * idempotent). Both give defect 0; the bracket collapses to [0, ~3e-16].
 *
 * Determinism: closed form inside aic_channel_depolarizing (sqrt of doubles). No
 * RNG. `out` is aic_ucp_kraus_init'd by aic_channel_depolarizing (caller clears
 * with aic_ucp_kraus_clear); this generator's VALUE is the named 2A instance +
 * its certified defect properties, so it DELEGATES rather than reimplementing.
 *
 * Measured anchors (prec=256, d=2), asserted by tests/test_adversarial.c:
 *   p=0.1: p(1-p)=0.09, bracket [0.077942, 0.311769]  (lo/q = 0.866025)
 *   p=0.5: p(1-p)=0.25, bracket [0.216506, 0.866025]  (lo/q = 0.866025, MAX)
 *   p=0.9: p(1-p)=0.09, bracket [0.077942, 0.311769]  (= p=0.1, symmetric)
 *   p=0 and p=1: bracket [0, ~3e-16]  (exact idempotent)
 */
void aic_adv_chan_depol_boundary(aic_ucp_kraus *out, slong d, double p,
                                 slong prec)
{
    assert(d >= 2 && "aic_adv_chan_depol_boundary: need d >= 2");
    assert(p >= 0.0 && p <= 1.0 &&
           "aic_adv_chan_depol_boundary: need 0 <= p <= 1");

    /* The named 2A instance IS the standard depolarizing channel; delegate
     * (it self-inits out, caller clears). Do not reimplement (aic_channels.h). */
    aic_channel_depolarizing(out, d, p, prec);
}
