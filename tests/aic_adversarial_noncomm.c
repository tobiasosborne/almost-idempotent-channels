/* aic_adversarial_noncomm.c — adversarial CHANNEL generator: a NON-COMMUTATIVE
 * UCP map CALIBRATED to the eta < 1/4 boundary eta_cb = ||Phi^2-Phi||_cb =
 * 1/4 - kappa (bead aic-cxo, the MED "non-commutative calibrated eta->1/4" item;
 * docs/adversarial/README.md:57, the tranche-3 gap). Sibling of the commutative
 * boundary generator aic_adv_chan_depol_boundary (family 2A,
 * aic_adversarial_domain.c): SAME boundary eta_cb -> 1/4, but a GENUINELY
 * NON-COMMUTATIVE image algebra, so it stresses the regularization basin
 * (theta(2Phi-1), tex:516-525) AND the non-commutative machinery (cstar_build,
 * ecstar, dhom) where the symmetric depolarizing G2 (image = scalars) never
 * exercises the off-diagonal structure. Kept in its own file (Rule 10). See
 * aic_adversarial.h for the corpus rationale (CLAUDE.md adversarial-first
 * directive).
 *
 * WHY NON-COMMUTATIVE AT THE BOUNDARY IS MORE ADVERSARIAL (the motivation).
 * The depolarizing 2A instance reaches rho(Phi^2-Phi) = p(1-p) = 1/4 at p=1/2,
 * but its image algebra is the SCALARS C*1 (abelian): the off-diagonal /
 * non-abelian dimension of cstar_build, errreduce, dhom, ecstar is NEVER
 * exercised at the boundary. A non-commutative image (a genuine M_m corner, m>=2)
 * forces those routines to run on a non-abelian algebra exactly where the basin
 * theta(2Phi-1) nearly fails (||(2Phi-1)^2-1||_cb = 4*eta_cb = 1 - 4*kappa -> 1)
 * and the th_main constant threatens to blow up.
 *
 * THE CONSTRUCTION (tensor a NON-COMMUTATIVE exact idempotent with a tunable
 * COMMUTATIVE defect, exploiting cb-norm AMPLIATION-INVARIANCE).
 *   Phi = id_{M_m} (x) Psi_eta   on   B(C^m (x) C^d) = B(C^N),  N = m*d,
 * where:
 *   - id_{M_m} is the IDENTITY channel on B(C^m) — EXACTLY idempotent and the
 *     source of NON-COMMUTATIVITY: its fixed-point algebra is the FULL M_m
 *     (non-abelian for m >= 2). Kraus operator: 1_m.
 *   - Psi_eta is the measure-prepare cb-norm example on B(C^d) (the family-1B
 *     generator aic_adv_chan_cb_op_gap, tex:366-388), whose idempotence defect is
 *     KNOWN IN CLOSED FORM (tex:378):
 *         ||Psi_eta^2 - Psi_eta||_cb = eta * sqrt(1 - eta).
 *     Psi_eta is UNITAL (sum_a K_a^dag K_a = 1_d), so 1_d is a fixed point.
 *   The Kraus operators of Phi are the Kronecker products  1_m (x) K_a^{Psi}
 *   (aic_mat_kronecker, LEFT factor major). UNITAL: sum_a (1_m (x) K_a)^dag
 *   (1_m (x) K_a) = 1_m (x) sum_a K_a^dag K_a = 1_m (x) 1_d = 1_N. UCP (a tensor
 *   of UCP maps). Kraus rank = r(Psi) = d.
 *
 * THE EXACT cb-NORM (the calibration target; cite this). cb-norm is AMPLIATION-
 * INVARIANT: for ANY map L, ||id_p (x) L||_cb = ||L||_cb for all p (the sup over
 * n in tex:347-349 already ranges over all ampliations; Paulsen, Completely
 * Bounded Maps, ch. 3). Since
 *   Phi^2 - Phi = (id_{M_m} (x) Psi)^2 - (id_{M_m} (x) Psi)
 *               = id_{M_m} (x) (Psi^2 - Psi),
 * we get the EXACT, m-INDEPENDENT identity
 *   eta_cb := ||Phi^2 - Phi||_cb = ||Psi^2 - Psi||_cb = eta * sqrt(1 - eta)
 *                                  (tex:347-349 ampliation-invariance + tex:378).
 * MEASURED CONFIRMATION (bead aic-cxo feasibility probe, prec=256): the eig-free
 * lower bound lo = ||J||_F/N is INVARIANT under the (x) id_{M_m} ampliation (it
 * stays at the m=1 value 0.272166 at eta=2/3 for m=1,2,3), and the certified
 * bracket [||J||_F/N, 2||J||_F] always CONTAINS eta*sqrt(1-eta). So the cb-norm
 * (NOT just the spectral radius rho — see FINDINGS C15/C17) DOES reach the 1/4
 * boundary, unlike make_mixconj whose rho stays ~0.165.
 *
 * THE CALIBRATION (knob = kappa; root-find eta). Given the target eta_cb = 1/4 -
 * kappa, we BISECTION-root-find eta in [0, 2/3] (the RISING branch of
 * g(eta) = eta*sqrt(1-eta), which increases monotonically from 0 to its peak
 * g(2/3) = (2/3)*sqrt(1/3) = 0.3849002 at eta = 2/3) such that g(eta) = 1/4-kappa.
 * The lower root (eta < 2/3) is the well-conditioned choice (smaller eta => the
 * measure-prepare vector v = (sqrt(1-eta), sqrt(eta), 0...) is closer to a clean
 * basis vector). g is monotone INCREASING and CONTINUOUS on [0,2/3], so bisection
 * is GUARANTEED to bracket the unique root when 0 <= target <= 0.3849002. FAIL
 * LOUD (Rule 4) if the target is unreachable: kappa <= 0 (target >= 1/4, violates
 * the eta < 1/4 hypothesis) or target > peak (no root, kappa < 0.25-0.3849<0,
 * impossible for kappa>0 but guarded) — assert, do NOT silently clamp.
 *
 * NON-COMMUTATIVITY WITNESS (distinguishes it from the commutative G2). The
 * elements B1 = E^{(m)}_{01} (x) 1_d and B2 = E^{(m)}_{10} (x) 1_d are GENUINE
 * FIXED POINTS of Phi (Psi unital => Psi(1_d)=1_d, and id_{M_m} fixes the M_m
 * factor), so under the Choi-Effros star product X*Y = Phi(XY) (tex:342) they
 * multiply as ordinary matrices:  B1*B2 - B2*B1 = Phi([B1,B2]) = [B1,B2] =
 * (E^{(m)}_{00} - E^{(m)}_{11}) (x) 1_d != 0, with ||B1*B2 - B2*B1||_op = 1
 * EXACTLY. For the commutative depolarizing 2A instance the image is abelian and
 * the analogous star-commutator on its fixed points VANISHES — the named
 * adversarial distinction. The self-test (test_adversarial.c) certifies this
 * commutator is bounded away from 0.
 *
 * kappa -> 1/4^- REDUCTION: as kappa -> 1/4, target -> 0, the root eta -> 0,
 * v -> |0>, Psi -> complete dephasing (EXACTLY idempotent), Phi -> id_{M_m} (x)
 * dephasing (exactly idempotent), defect -> 0. The bracket collapses toward [0,0]
 * (consistent with eta_cb = 1/4 - kappa -> 0). NOT used as the oracle here (the
 * point is the BOUNDARY kappa -> 0), but it confirms the construction degrades
 * gracefully.
 *
 * Determinism: closed-form bisection (fixed iteration count) + exact arb sqrt in
 * cb_op_gap; no RNG. The same (m,d,kappa,prec) yields a bit-identical channel.
 *
 * Measured anchors (prec=256), asserted by tests/test_adversarial.c:
 *   eta_cb target 1/4-kappa, lower-root eta, eig-free bracket [lo,hi] (N=m*d):
 *     kappa=0.10 target=0.150 eta=0.16406042 bracket=[0.106066,0.848528]
 *     kappa=0.05 target=0.200 eta=0.22756104 bracket=[0.141421,1.131371]
 *     kappa=0.02 target=0.230 eta=0.26901280 bracket=[0.162635,1.301076]
 *     kappa=0.01 target=0.240 eta=0.28354076 bracket=[0.169706,1.357645]
 *   (each bracket CONTAINS the target 1/4-kappa; lo<target<hi). NON-COMMUTATIVITY
 *   witness ||B1*B2-B2*B1||_op = 1.0 (m=2). UNITAL defect < 1e-15.
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

/* g(eta) = eta*sqrt(1-eta), the EXACT cb-norm defect of Psi_eta (tex:378). On
 * [0, 2/3] g is continuous and strictly increasing from g(0)=0 to the peak
 * g(2/3) = (2/3)*sqrt(1/3) = 0.3849001794597505. */
static double noncomm_g(double eta) { return eta * sqrt(1.0 - eta); }

/* Lower-root calibration: bisection for eta in [0, 2/3] with g(eta) = target.
 * Returns the unique root on the rising branch. Caller has asserted
 * 0 <= target <= g(2/3); g is monotone so bisection converges. 200 iterations
 * drive the bracket below 2^-200 * (2/3) — far under the arb prec floor. */
static double noncomm_calibrate_eta(double target)
{
    double a = 0.0, b = 2.0 / 3.0; /* rising branch [0, argmax] */
    /* invariant g(a) <= target <= g(b) holds: g(0)=0<=target, g(2/3)=peak>=target */
    for (int it = 0; it < 200; it++) {
        double mid = 0.5 * (a + b);
        if (noncomm_g(mid) < target)
            a = mid;
        else
            b = mid;
    }
    return 0.5 * (a + b);
}

/* fam-NC — NON-COMMUTATIVE eta->1/4 boundary (docs/adversarial/README.md:57,
 * tex:347-349 cb ampliation-invariance + tex:366-388/378 the measure-prepare
 * closed form + tex:516-525 the theta(2Phi-1) basin). Phi = id_{M_m} (x)
 * Psi_eta on B(C^N), N = m*d, with eta the lower root of eta*sqrt(1-eta) =
 * 1/4 - kappa, so eta_cb = ||Phi^2-Phi||_cb = 1/4 - kappa EXACTLY (closed form,
 * ampliation-invariant). Kraus: 1_m (x) K_a^{Psi} (Kronecker, LEFT major). UNITAL
 * for all kappa; NON-COMMUTATIVE image (the M_m corner, m>=2). Asserts m>=2 (the
 * non-commutativity needs a >=2-dim matrix corner), d>=2, 0 < kappa < 1/4 (the
 * eta<1/4 hypothesis; Rule 4, fail loud — a target >= 1/4 is OUT of hypothesis,
 * a target > peak 0.3849 has NO root). `out` is aic_ucp_kraus_init'd HERE (caller
 * clears with aic_ucp_kraus_clear). See the file docstring for the full
 * derivation, the cb-norm identity, the calibration, and the witness. */
void aic_adv_chan_noncomm_boundary(aic_ucp_kraus *out, slong m, slong d,
                                   double kappa, slong prec)
{
    assert(m >= 2 && "aic_adv_chan_noncomm_boundary: need m >= 2 (M_m "
                     "non-commutative corner)");
    assert(d >= 2 && "aic_adv_chan_noncomm_boundary: need d >= 2 (cb_op_gap "
                     "inner dim)");
    assert(kappa > 0.0 && kappa < 0.25 &&
           "aic_adv_chan_noncomm_boundary: need 0 < kappa < 1/4 (eta<1/4 "
           "hypothesis; target = 1/4 - kappa)");

    double target = 0.25 - kappa; /* eta_cb to hit, in (0, 1/4) */
    /* FAIL LOUD (Rule 4) if the target is above the reachable peak: there is then
     * NO root of eta*sqrt(1-eta) = target, the construction genuinely cannot
     * reach the boundary, and silently clamping would fake it. (For kappa>0 this
     * cannot trigger — target<1/4<0.3849 — but the guard documents the obstruction
     * the §C15/§C17 feasibility question is about, and protects a future caller.) */
    assert(target <= noncomm_g(2.0 / 3.0) &&
           "aic_adv_chan_noncomm_boundary: target 1/4-kappa exceeds the peak "
           "eta*sqrt(1-eta) = 0.3849; no calibration root (cannot reach boundary)");

    double eta = noncomm_calibrate_eta(target);

    slong N = m * d;
    /* The inner measure-prepare channel Psi_eta on B(C^d) (family-1B generator;
     * reuse, do not reimplement — Law 2). It self-inits psi; we clear it below. */
    aic_ucp_kraus psi;
    aic_adv_chan_cb_op_gap(&psi, d, eta, prec);

    /* Phi = id_{M_m} (x) Psi: dim_K = dim_H = N, Kraus rank r = r(Psi) = d, each
     * Kraus op K_a^{Phi} = 1_m (x) K_a^{Psi} (Kronecker, LEFT factor major). */
    aic_ucp_kraus_init(out, N, N, psi.r);
    acb_mat_t Im;
    acb_mat_init(Im, m, m);
    acb_mat_one(Im); /* 1_m */
    for (slong a = 0; a < psi.r; a++)
        aic_mat_kronecker(out->K[a], Im, psi.K[a], prec);

    acb_mat_clear(Im);
    aic_ucp_kraus_clear(&psi);
}
