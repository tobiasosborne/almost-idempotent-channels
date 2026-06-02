/* aic_adversarial_projsplit.c — adversarial NLA generator for the NEAR-TRIVIAL-
 * PROJECTION / Route-A spectral-split-failure family (domain.md:373-412, family
 * 3C; #7 on the merged lethal shortlist, bead aic-dbo.2). A bare Hermitian
 * acb_mat_t H whose spectrum is two clusters near {0, 1} with a SMALLEST INTERIOR
 * SEPARATION tunable to exactly `gap_spec`, the knob gap_spec/eps sweeping the
 * deliver-or-refuse boundary of the nontrivial-projection finder (aic_projection).
 *
 * WHY THIS FAMILY (domain.md:384-390, tex:516-525). lem_nontriv_projection
 * (tex:931) Route A — the only constructive route — needs (Step 3) a threshold t
 * with no eigenvalue in (t-g, t+g) for a gap g>0, then snaps a sign matrix via
 * prop_P (tex:524-533): theta(2P-I) = (1/2)(I + sgn(2P-I)), whose funcalc basin is
 * ||X^2 - I|| < 1 (tex:516,520). The finder (src/aic_projection_find.c) realizes
 * Step 3 as: pick the LARGEST interior gap of a Hermitian element H, set
 * t=(lam_i+lam_{i+1})/2, s=1/max(t-lam_min, lam_max-t), Y=s(H-tI), and X=sgn(Y);
 * the basin is then ||Y^2-I|| = 1 - (s*g/2)^2 < 1 for ANY g>0. So the SIZE of the
 * gap is NOT a hard precondition — it enters the A-POSTERIORI star defect
 * delta=O(eps/g) (include/aic_projection.h:39) when the ambient projector is
 * projected into a genuine eps-C* algebra. The ONLY fail-loud wall is a GENUINELY
 * DEGENERATE spectrum (no positive interior gap above the round-off floor
 * 1e-9*max(1,spread)): the aic-3qv "NO positive interior spectral gap" abort
 * (src/aic_projection_find.c:142-151 / src/aic_projection_audit.c:98-106). This
 * generator is the witness that drives a Hermitian element's largest interior gap
 * smoothly from g=Omega(1) (DELIVER: the finder snaps a genuine nontrivial
 * spectral projector) down to g->0 (REFUSE: the aic-3qv abort fires), so the
 * deliver-or-refuse contract of the aic-66n unit-aware finder can be exercised
 * across the boundary.
 *
 * THE OUTPUT TYPE (a bare Hermitian, NOT a channel). Section 7 of the brief asked
 * whether 3C should produce an aic_ecstar (an engineered algebra) or an
 * aic_ucp_kraus (a channel -> A via aic_assoc_ecstar_from_phi). We MEASURED
 * (the 2026-06-02 probes recorded in FINDINGS §D1) that NEITHER algebra route can
 * present a tunable SMALL largest-interior-gap to the finder: any genuine 2+-dim
 * eps-C* algebra's Frobenius-ORTHONORMAL SVD basis (the {B_k} of aic_ecstar)
 * contains an element whose Hermitian part H_k = (B_k+B_k^dag)/2 is a near-
 * projector with an O(1) interior gap (a gauge artifact of orthonormalization:
 * e.g. the 2-block conditional expectation on C^n always yields an H_k with
 * spectrum {0,...,0,1}, gap 1.0). The audition ranges over ALL H_k and picks the
 * LARGEST gap, so a genuine algebra ALWAYS DELIVERS at an O(1) gap — the small-gap
 * regime simply does not exist at the algebra level except at dim_A=1 (the T4
 * dim-1 refuse) or a degenerate spectrum. The tunable gap therefore lives on the
 * Hermitian element the finder's Route-A-Step-3 primitives (aic_projection_gap,
 * aic_projection_ambient) consume DIRECTLY, which is exactly the op domain.md:392
 * names ("projection: Route A Step 3; arb gap certification" and "funcalc: prop_P
 * precondition ||X^2-I||<1"). So 3C is an NLA generator (a caller-managed
 * acb_mat_t H, the gen3/gen4/gen6/gen7 convention), NOT a channel. The self-test
 * (tests/test_adversarial.c test_fam3c_*) then drives H through aic_projection_gap
 * + aic_projection_ambient (DELIVER) and the no-gap abort (REFUSE), and through
 * aic_projection_into_A on a genuine eta>0 algebra to certify the O(eta) star
 * defect.
 *
 * THE CONSTRUCTION (a diagonal two-cluster Hermitian). On C^n, a diagonal real
 * (hence Hermitian) matrix H = diag(lower cluster, upper cluster):
 *   - LOWER cluster: k eigenvalues in [0, w], a deterministic ramp
 *       lam_i = w * i/(k-1),   i = 0 .. k-1   (or {0} if k==1),
 *   - UPPER cluster: (n-k) eigenvalues in [w+gap_spec, 2w+gap_spec],
 *       lam_{k+j} = (w + gap_spec) + w * j/(n-k-1),   j = 0 .. n-k-1,
 *   where w = WITHIN-CLUSTER WIDTH = gap_spec/4 (each cluster is 4x TIGHTER than
 *   the separation). So:
 *     * the SMALLEST INTERIOR gap BETWEEN the clusters is EXACTLY gap_spec
 *       (lam_k - lam_{k-1} = (w+gap_spec) - w = gap_spec),
 *     * EVERY within-cluster gap is <= w/(min(k,n-k)-1) <= w = gap_spec/4 < gap_spec,
 *       so the cluster SEPARATION gap_spec is unambiguously the LARGEST interior
 *       gap — the one aic_projection_gap selects AND certifies. The named
 *       adversarial quantity (the smallest INTER-cluster separation) IS the gap
 *       the finder uses. As gap_spec -> 0 with w=gap_spec/4 -> 0 too, the WHOLE
 *       spectrum collapses to a point and the finder hits the aic-3qv no-gap
 *       abort (degenerate spectrum) — the REFUSE wall.
 *   The spectrum is clustered near {0, ~1} when gap_spec ~ 1 (e.g. gap_spec=0.5,
 *   w=0.125: lower in [0,0.125], upper in [0.625,0.75]); for the small-gap knobs
 *   the two clusters sit near {0, gap_spec} (the {0,1}-near regime of domain.md is
 *   the gap_spec=Omega(1) end; the named EVIL regime is gap_spec small, where the
 *   clusters merge). The eps argument is recorded only as the SCALE the knob ratio
 *   gap_spec/eps references; the construction depends only on gap_spec (eps fixes
 *   the sweep point gap_spec = ratio * eps). w=gap_spec/4 is INDEPENDENT of eps so
 *   the cluster separation is ALWAYS the largest interior gap regardless of the
 *   ratio — the realized largest gap is gap_spec to machine precision (MEASURED
 *   rel_err <= 2.8e-16 across n in {4,9,16}, ratio in {10,3,1.5,1.0,0.5}; the
 *   self-test's construction pin).
 *
 * DETERMINISM. Closed form (a real diagonal of rational ramps of gap_spec). No
 * RNG. The same (n, k, gap_spec) yields a bit-identical H.
 *
 * Asserts (Rule 4, fail loud): n >= 2, 0 < k < n (an INTERIOR split needs both
 * clusters non-empty), gap_spec >= 0. gap_spec == 0 is allowed (the exactly-
 * degenerate refuse witness: H collapses to a single cluster, no interior gap);
 * the caller drives the aic-3qv abort with it. `out` must be init'd n x n
 * (asserted), matching the NLA-generator convention (caller clears it).
 *
 * Measured anchors (prec=256), asserted by tests/test_adversarial.c test_fam3c_*:
 *   gap calibration: realized largest interior gap == gap_spec to <= 2.8e-16
 *     across n in {4,9,16}, ratio in {10,3,1.5,1.0,0.5}, eps in {0.01,0.05}.
 *   DELIVER (ratio 10,3): aic_projection_gap returns g==gap_spec; the prop_P basin
 *     ||Y^2-I|| = 1-(s*gap/4)^2... here = 0.5556 < 1 (w=gap/4 makes it gap-
 *     invariant); aic_projection_ambient builds an EXACT ambient idempotent
 *     (||P_amb^2-P_amb|| ~ 1e-74, rank k); projecting into a genuine eta>0 algebra
 *     (mixconj(5,3,0.02), eta~0.046) gives a finite certified star defect ~3.6e-3
 *     = O(eta), gap-invariant for a clean ambient projector (FINDINGS §D1).
 *   REFUSE (gap_spec=0, single cluster): aic_projection_gap fail-loud-ABORTS with
 *     "NO positive interior spectral gap" (the aic-3qv stop condition), via the
 *     fork+SIGALRM watchdog so the test binary survives.
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_adversarial.h"

/* fam3C — near-trivial-projection two-cluster Hermitian (domain.md:373-412,
 * tex:516-525/931). Diagonal real H on C^n: k eigenvalues in [0, gap_spec/4]
 * (lower cluster) and (n-k) in [gap_spec/4 + gap_spec, gap_spec/2 + gap_spec]
 * (upper cluster), so the SMALLEST INTER-cluster gap = the LARGEST interior gap =
 * EXACTLY gap_spec (within-cluster gaps <= gap_spec/4 < gap_spec). KNOB gap_spec
 * (>= 0; sweep gap_spec = ratio*eps, ratio in {10,3,1.5,1.0,0.5}). Asserts n>=2,
 * 0<k<n, gap_spec>=0 (Rule 4, fail loud). `out` must be init'd n x n (caller
 * clears it). See the file docstring for the full derivation, the output-type
 * decision (bare Hermitian, NOT a channel — FINDINGS §D1), and the deliver-or-
 * refuse boundary. */
void aic_adv_proj_near_trivial(acb_mat_t out, slong n, slong k, double gap_spec,
                               slong prec)
{
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_adv_proj_near_trivial: out must be n x n");
    assert(n >= 2 && "aic_adv_proj_near_trivial: need n >= 2");
    assert(k > 0 && k < n &&
           "aic_adv_proj_near_trivial: need 0 < k < n (interior split)");
    assert(gap_spec >= 0.0 &&
           "aic_adv_proj_near_trivial: need gap_spec >= 0");
    (void) prec;   /* the eigenvalues are exact doubles set into the diagonal */

    acb_mat_zero(out);

    /* within-cluster width w = gap_spec/4: each cluster spans <= w, so EVERY
     * within-cluster gap <= w/(cluster_size-1) <= w = gap_spec/4 < gap_spec, and
     * the inter-cluster gap is EXACTLY gap_spec — the largest interior gap. */
    double w = 0.25 * gap_spec;

    /* LOWER cluster: k eigenvalues, lam_i = w * i/(k-1) in [0, w] (or {0}). */
    for (slong i = 0; i < k; i++) {
        double v = (k > 1) ? w * (double) i / (double) (k - 1) : 0.0;
        acb_set_d(acb_mat_entry(out, i, i), v);
    }
    /* UPPER cluster: (n-k) eigenvalues, base = w + gap_spec, spread upward by w. */
    double base = w + gap_spec;
    slong u = n - k;
    for (slong j = 0; j < u; j++) {
        double v = base + ((u > 1) ? w * (double) j / (double) (u - 1) : 0.0);
        acb_set_d(acb_mat_entry(out, k + j, k + j), v);
    }
}
