/* test_assoc2.c — cross-checks for assoc_ecsa Increment 2 (bead aic-92f): the
 * associated eps-C* algebra A = Img Phi_tilde, the approximate Choi-Effros star
 * X*Y = Phi_tilde(XY), and the certification that (A, *, ||.||, dag, I) is an
 * extended O(eta)-C* algebra (th_almost_idemp, approximate_algebras.tex:2184-2237).
 * Reuses the exact-idempotent channels in tests/test_idemp.h and the ecstar
 * axiom-defect estimators (aic_ecstar_*) as the verification tooling.
 *
 * Cross-check ladder:
 *   U1 eta=0 CONSISTENCY ORACLE: for each exact-idempotent channel, Phi_tilde=Phi,
 *      so A=Img Phi_tilde must MATCH ecstar's A=Img Phi (aic_ecstar_from_idemp):
 *      same dim_A; same SUBSPACE (Pi_A projectors agree, A has unitary gauge
 *      freedom so compared as subspaces); ALL FIVE axiom defects machine-zero
 *      (basis-sweep AND assoc/cstar HOPM). Ties Increment 2 to the trusted oracle.
 *   U2 REAL eta>0 (bead aic-4c7): on dep(d)+conj(dep(d)), the FIVE defects are
 *      <= C*eta, -> 0 as t->0, and NONZERO for t>0 (genuinely approximate, not
 *      exact -- teeth). Uses HOPM for the faithful worst case.
 *   U3 DIMENSION-SWEEP UNIVERSALITY CANARY (bead aic-dbo.3): at FIXED small t
 *      across a WIDE sweep dim_A=4,6,8,10,12,16, HOPM_assoc/eta does NOT grow with
 *      dim_A (flat; real spread 1.42, threshold 2.0). The naive Haar route fails
 *      ~ dim (.tex:484); the threshold is mutation-proven to catch BOTH an injected
 *      linear ~dim factor (spread ~4.67) AND a milder ~sqrt(dim) factor (spread
 *      ~2.33). INTEGRITY guard: dim_A varies (dA_max>dA_min) so the canary is not
 *      vacuously flat-by-construction.
 *   U4 Pi_A ACCEPT-GUARD on a NON-POLAR-CLOSED A (bead aic-3qq): the eta>0 assoc A
 *      is NOT polar-closed, so the HOPM witness exercises the Pi_A(polar(C)) +
 *      monotone-accept correctness half. Verify the witness is IN A (proj_residual
 *      ~0, certified arb) and lo <= ratio@witness.
 *
 * MUTATION-PROVEN new teeth (RED observed, restored): (a) BOTTOM-dim_A SVD vectors
 * in aic_assoc_extract_range -> U1 subspace match RED; (b) star_phi applying
 * S_tilde to X instead of XY -> U1 eta=0 assoc oracle RED; (c) a dim-TRIVIAL
 * canary family -> U3 integrity guard RED. See the comments at each test.
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
#include "aic_ecstar.h"
#include "aic_funcalc.h"
#include "aic_idemp.h"
#include "aic_latd.h"
#include "aic_mat.h"
#include "aic_test.h"
#include "aic_ucp.h"
#include "test_idemp.h"

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Phi_t = (1-t) base + t mix as the Kraus union {sqrt(1-t) Ka} U {sqrt(t) Lb}.
 * Convex combo of UCP maps is UCP; at t=0 it is `base`. `out` init'd here. */
static void make_mix(aic_ucp_kraus *out, const aic_ucp_kraus *base,
                     const aic_ucp_kraus *mix, double t, slong prec)
{
    slong n = base->dim_H;
    AIC_CHECK_MSG(base->dim_K == n && mix->dim_K == n && mix->dim_H == n,
                  "make_mix: base/mix not self-maps on the same C^n");
    aic_ucp_kraus_init(out, n, n, base->r + mix->r);
    arb_t s;
    arb_init(s);
    arb_set_d(s, sqrt(1.0 - t));
    for (slong a = 0; a < base->r; a++)
        acb_mat_scalar_mul_arb(out->K[a], base->K[a], s, prec);
    arb_set_d(s, sqrt(t));
    for (slong b = 0; b < mix->r; b++)
        acb_mat_scalar_mul_arb(out->K[base->r + b], mix->K[b], s, prec);
    arb_clear(s);
}

/* The dep(d)+conj(dep(d)) eta>0 family: dephasing(d) mixed with a complex-unitary
 * conjugate of dephasing(d). Two DIFFERENT commutative algebras (diagonal, and a
 * unitarily-rotated diagonal); Phi_tilde's image is genuinely a NON-algebra, so
 * the associator is nonzero O(eta). dim_A = d (scales with dimension). `out`
 * init'd here. */
static void make_eta_family(aic_ucp_kraus *out, slong d, double t, slong prec)
{
    aic_ucp_kraus base, m0, mix;
    make_dephasing(&base, d);
    make_dephasing(&m0, d);
    make_conjugated(&mix, &m0, prec);
    make_mix(out, &base, &mix, t, prec);
    aic_ucp_kraus_clear(&mix);
    aic_ucp_kraus_clear(&m0);
    aic_ucp_kraus_clear(&base);
}

/* op-norm of the midpoint of S_Phi^2 - S_Phi (the rigorous eta proxy: the
 * idempotence defect of the superop, == ||Phi^2-Phi||_op <= ||.||_cb = eta; the
 * certified cb eta needs the Watrous SDP, not here -- Increment 1 T7 precedent). */
static double eta_proxy(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S, S2, D, M;
    arb_t e;
    acb_mat_init(S, nn, nn);
    acb_mat_init(S2, nn, nn);
    acb_mat_init(D, nn, nn);
    acb_mat_init(M, nn, nn);
    arb_init(e);
    aic_assoc_superop_from_ucp(S, phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    for (slong i = 0; i < nn; i++)
        for (slong j = 0; j < nn; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    aic_mat_opnorm(e, M, prec);
    double o = arf_get_d(arb_midref(e), ARF_RND_NEAR);
    arb_clear(e);
    acb_mat_clear(M); acb_mat_clear(D); acb_mat_clear(S2); acb_mat_clear(S);
    return o;
}

/* Build the subspace projector Pi_A = sum_k vec(B_k) vec(B_k)^dag (n^2 x n^2,
 * vec_r). Gauge-invariant: two bases of the SAME subspace give the same Pi_A
 * (HANDOFF: compare A as a SUBSPACE, never operator-by-operator). */
static void assoc_build_Pi_A(acb_mat_t Pi, const acb_mat_t *B, slong dim_A, slong n,
                       slong prec)
{
    slong nn = n * n;
    acb_mat_zero(Pi);
    acb_t t;
    acb_init(t);
    for (slong k = 0; k < dim_A; k++)
        for (slong a = 0; a < nn; a++)        /* row a = (i*n+j) */
            for (slong b = 0; b < nn; b++) {  /* col b = (p*n+q) */
                /* vec(B_k)[a] conj(vec(B_k)[b]) */
                acb_conj(t, acb_mat_entry(B[k], b / n, b % n));
                acb_mul(t, acb_mat_entry(B[k], a / n, a % n), t, prec);
                acb_add(acb_mat_entry(Pi, a, b), acb_mat_entry(Pi, a, b), t, prec);
            }
    acb_clear(t);
}

/* Top singular value sigma_max of the (already-midpoint, exactly-idempotent)
 * superop S_tilde via the double-path SVD (the same engine the extractor uses).
 * For an idempotent the nonzero sigmas are >= 1: == 1 iff S_tilde is an orthogonal
 * projector (Phi_tilde HS-self-adjoint), > 1 iff genuinely oblique (Fix C). */
static double sigma_max_S_tilde(const acb_mat_t S)
{
    slong nn = acb_mat_nrows(S);
    double _Complex *Sd = flint_malloc((size_t) (nn * nn) * sizeof(double _Complex));
    double *sv = flint_malloc((size_t) nn * sizeof(double));
    aic_latd_from_acb_mat(Sd, S);
    aic_latd_svd(sv, NULL, NULL, Sd, nn, nn);  /* singular values only, descending */
    double smax = sv[0];
    flint_free(sv);
    flint_free(Sd);
    return smax;
}

/* ============================ U1: eta=0 oracle =========================== */

/* For an EXACT idempotent phi: build A both ways (assoc Phi_tilde and ecstar
 * Img Phi), assert same dim_A, same SUBSPACE (||Pi_assoc - Pi_idemp||_op ~0), and
 * all five basis-sweep defects + assoc/cstar HOPM machine-zero on the assoc A. */
static void u1_channel(const char *name, aic_ucp_kraus *phi, slong prec)
{
    /* ecstar A = Img Phi (the trusted oracle). */
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, phi, prec);
    aic_ecstar Ai;
    aic_ecstar_init(&Ai, d.n, d.dim_A, phi);
    aic_ecstar_from_idemp(&Ai, &d, phi, prec);

    /* assoc A = Img Phi_tilde. */
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, prec);

    slong n = d.n, nn = n * n;
    AIC_CHECK_MSG(ae.A.dim_A == d.dim_A,
                  "%s: assoc dim_A=%ld != idemp dim_A=%ld", name,
                  (long) ae.A.dim_A, (long) d.dim_A);

    /* (ii) SUBSPACE match: Pi_assoc == Pi_idemp (gauge-invariant). */
    acb_mat_t Pa, Pi, Diff;
    acb_mat_init(Pa, nn, nn);
    acb_mat_init(Pi, nn, nn);
    acb_mat_init(Diff, nn, nn);
    assoc_build_Pi_A(Pa, (const acb_mat_t *) ae.A.B, ae.A.dim_A, n, prec);
    assoc_build_Pi_A(Pi, (const acb_mat_t *) Ai.B, Ai.dim_A, n, prec);
    acb_mat_sub(Diff, Pa, Pi, prec);
    arb_t e, tol;
    arb_init(e);
    arb_init(tol);
    arb_set_d(tol, 1e-9);  /* both bases are double-path SVD extractions */
    aic_mat_opnorm(e, Diff, prec);
    AIC_CHECK_MSG(arb_le(e, tol), "%s: assoc/idemp subspaces differ ||dPi||=%.3e",
                  name, dd(e));
    double subdiff = dd(e);

    /* (iii) all five basis-sweep defects machine-zero on the assoc A. */
    arb_t da, ds, dc, di, du;
    arb_init(da); arb_init(ds); arb_init(dc); arb_init(di); arb_init(du);
    aic_ecstar_defect_assoc(da, &ae.A, prec);
    aic_ecstar_defect_submult(ds, &ae.A, prec);
    aic_ecstar_defect_cstar(dc, &ae.A, prec);
    aic_ecstar_defect_involution(di, &ae.A, prec);
    aic_ecstar_defect_unit(du, &ae.A, prec);
    /* assoc/cstar HOPM (faithful worst case) also ~0. */
    arb_t ha, hc;
    arb_init(ha); arb_init(hc);
    aic_ecstar_defect_assoc_hopm(ha, &ae.A, 4, 8, prec);
    aic_ecstar_defect_cstar_hopm(hc, &ae.A, 4, 8, prec);

    arb_t z;
    arb_init(z);
    arb_set_d(z, 1e-7);
    AIC_CHECK_MSG(arb_le(da, z), "%s: U1 assoc defect not ~0 (%.3e)", name, dd(da));
    AIC_CHECK_MSG(arb_le(ds, z), "%s: U1 submult defect not ~0 (%.3e)", name, dd(ds));
    AIC_CHECK_MSG(arb_le(dc, z), "%s: U1 cstar defect not ~0 (%.3e)", name, dd(dc));
    AIC_CHECK_MSG(arb_le(di, z), "%s: U1 involution defect not ~0 (%.3e)", name, dd(di));
    AIC_CHECK_MSG(arb_le(du, z), "%s: U1 unit defect not ~0 (%.3e)", name, dd(du));
    AIC_CHECK_MSG(dd(ha) < 1e-7, "%s: U1 HOPM assoc not ~0 (%.3e)", name, dd(ha));
    AIC_CHECK_MSG(dd(hc) < 1e-7, "%s: U1 HOPM cstar not ~0 (%.3e)", name, dd(hc));

    printf("U1 %-26s dim_A=%2ld subdiff=%.2e | assoc=%.2e submult=%.2e cstar=%.2e "
           "invol=%.2e unit=%.2e | HOPM a=%.2e c=%.2e\n",
           name, (long) ae.A.dim_A, subdiff, dd(da), dd(ds), dd(dc), dd(di),
           dd(du), dd(ha), dd(hc));

    arb_clear(z); arb_clear(hc); arb_clear(ha);
    arb_clear(du); arb_clear(di); arb_clear(dc); arb_clear(ds); arb_clear(da);
    arb_clear(tol); arb_clear(e);
    acb_mat_clear(Diff); acb_mat_clear(Pi); acb_mat_clear(Pa);
    aic_assoc_ecstar_clear(&ae);
    aic_ecstar_clear(&Ai);
    aic_idemp_clear(&d);
}

static void test_u1(void)
{
    const slong P = 256;
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 5, 2); u1_channel("block_cond_exp(5,2)", &phi, P);
    aic_ucp_kraus_clear(&phi);
    make_block_cond_exp(&phi, 3, 1); u1_channel("block_cond_exp(3,1)", &phi, P);
    aic_ucp_kraus_clear(&phi);
    make_trace_replace(&phi, 3);     u1_channel("trace_replace(3)", &phi, P);
    aic_ucp_kraus_clear(&phi);
    make_noiseless_subsystem(&phi, 2, 2);
    u1_channel("noiseless_subsystem(2,2)", &phi, P);
    aic_ucp_kraus_clear(&phi);
    make_identity(&phi, 3);          u1_channel("identity(3)", &phi, P);
    aic_ucp_kraus_clear(&phi);
    make_dephasing(&phi, 4);         u1_channel("dephasing(4)", &phi, P);
    aic_ucp_kraus_clear(&phi);
    make_compress_idemp(&phi, 4, 2); u1_channel("compress_idemp(4,2)", &phi, P);
    aic_ucp_kraus_clear(&phi);
    /* asymmetric conjugated channel (transpose-non-closed image -- the reshape
     * tooth; MUTATION: bottom-dim_A SVD vectors / transposed reshape make the
     * subspace match RED here). */
    {
        aic_ucp_kraus base, conj;
        make_compress_idemp(&base, 4, 2);
        make_conjugated(&conj, &base, P);
        u1_channel("conj_compress_idemp(4,2)", &conj, P);
        aic_ucp_kraus_clear(&conj);
        aic_ucp_kraus_clear(&base);
    }
}

/* ============================ U2: real eta>0 ============================= */

/* On dep(4)+conj(dep(4)) across a t-sweep: the five defects are <= C*eta, vanish
 * as t->0, and are NONZERO for t>0 (genuinely approximate). HOPM for assoc/cstar/
 * submult (faithful), basis-sweep for involution/unit.
 *
 * BOTH the associator AND the C*-identity are ASSERTED at eta>0 (FIX A, hostile
 * review R6): the C*-axiom (.tex:2233-2236, ||X^dag * X|| >= (1-O(eta))||X||^2) is
 * one of the three O(eta) axioms th_almost_idemp is about, and the cstar HOPM
 * defect c = sup_X (1 - ||X^dag*X||_op/||X||_op^2) is a GENUINE nonzero tooth at
 * eta>0 (Phi_tilde is HP but NOT positive, .tex:363, so the C* lower bound is
 * genuinely violated by O(eta)). Both ratios use the SAME eta proxy (eta_proxy:
 * ||S_Phi^2-S_Phi||_op) so c/eta and assoc/eta are apples-to-apples. */
static void test_u2(void)
{
    const slong P = 256;
    const slong d = 4;
    const double ts[] = {0.0, 0.01, 0.03, 0.06};
    const int nt = 4;
    double aprev = -1.0, cprev = -1.0;
    double ratio_lo = 1e300, ratio_hi = 0.0;
    double cratio_lo = 1e300, cratio_hi = 0.0;
    printf("U2 dep(4)+conj(dep(4)) eta>0 sweep:\n");
    for (int it = 0; it < nt; it++) {
        aic_ucp_kraus phi;
        make_eta_family(&phi, d, ts[it], P);
        double eta = eta_proxy(&phi, P);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, P);
        arb_t ha, hc, hs, di, du;
        arb_init(ha); arb_init(hc); arb_init(hs); arb_init(di); arb_init(du);
        aic_ecstar_defect_assoc_hopm(ha, &ae.A, 12, 16, P);
        aic_ecstar_defect_cstar_hopm(hc, &ae.A, 10, 14, P);
        aic_ecstar_defect_submult_hopm(hs, &ae.A, 10, 14, P);
        aic_ecstar_defect_involution(di, &ae.A, P);
        aic_ecstar_defect_unit(du, &ae.A, P);
        double a = dd(ha), c = dd(hc), s = dd(hs), iv = dd(di), u = dd(du);
        printf("  t=%.2f eta=%.4e dim_A=%ld: assoc=%.4e cstar=%.4e submult=%.4e "
               "invol=%.2e unit=%.2e", ts[it], eta, (long) ae.A.dim_A, a, c, s, iv, u);
        if (ts[it] == 0.0) {
            /* eta=0: everything machine-zero (eta_proxy is ~0 too -> skip ratio) */
            printf("  (eta=0 baseline)\n");
            AIC_CHECK_MSG(a < 1e-7 && fabs(c) < 1e-7 && fabs(s) < 1e-7 &&
                          iv < 1e-9 && u < 1e-9,
                          "U2: defects nonzero at t=0 (a=%.2e c=%.2e s=%.2e)", a, c, s);
        } else {
            double ra = a / eta, rc = c / eta;
            if (ra < ratio_lo) ratio_lo = ra;
            if (ra > ratio_hi) ratio_hi = ra;
            if (rc < cratio_lo) cratio_lo = rc;
            if (rc > cratio_hi) cratio_hi = rc;
            printf("  assoc/eta=%.4f cstar/eta=%.4f\n", ra, rc);
            /* O(eta): assoc/eta bounded by a chosen constant (universal, dim-indep).*/
            AIC_CHECK_MSG(ra <= 1.0,
                          "U2: assoc/eta=%.3f exceeds 1 (t=%.2f) -- not O(eta)?",
                          ra, ts[it]);
            /* GENUINELY APPROXIMATE: assoc nonzero for t>0 (teeth). */
            AIC_CHECK_MSG(a > 1e-6,
                          "U2: assoc ~0 at t=%.2f -- A is exactly associative? "
                          "(family not eta>0 / search failed)", ts[it]);
            /* C*-IDENTITY (.tex:2233-2236, FIX A): the cstar HOPM defect c is a
             * GENUINE nonzero tooth at eta>0 -- Phi_tilde is NOT positive (.tex:363)
             * so ||X^dag*X|| drops O(eta) below ||X||^2 for the worst-case X in A.
             * Mutating aic_ecstar_defect_cstar_hopm to return 0 must fire this. */
            AIC_CHECK_MSG(c > 1e-6,
                          "U2: cstar ~0 at t=%.2f -- C*-identity exact? "
                          "(Phi_tilde positive / search failed)", ts[it]);
            /* O(eta): cstar/eta bounded by a modest constant (universal, dim-indep;
             * SAME eta proxy as assoc, apples-to-apples). Observed c/eta ~0.75-0.82
             * across t (and ~constant over dim in U3); 1.5 is a modest margin. */
            AIC_CHECK_MSG(rc <= 1.5,
                          "U2: cstar/eta=%.3f exceeds 1.5 (t=%.2f) -- not O(eta)?",
                          rc, ts[it]);
            /* submult stays ~0: any star with ||Phi_tilde||_cb<=1+O(eta) is sub-
             * multiplicative (.tex:2216-2219); the upward slack stays tiny. */
            AIC_CHECK_MSG(s < 1e-2,
                          "U2: submult slack %.3e too large at t=%.2f", s, ts[it]);
            /* involution stays ~0: Phi_tilde is HP (.tex:2182). */
            AIC_CHECK_MSG(iv < 1e-9, "U2: involution moved (Phi_tilde HP)? %.3e", iv);
            /* assoc grows with t (-> 0 as t->0). */
            AIC_CHECK_MSG(a > aprev,
                          "U2: assoc not increasing with t (%.3e <= prev %.3e)",
                          a, aprev);
            /* cstar grows with t (-> 0 as t->0), mirroring assoc. */
            AIC_CHECK_MSG(c > cprev,
                          "U2: cstar not increasing with t (%.3e <= prev %.3e)",
                          c, cprev);
        }
        aprev = a;
        cprev = c;
        arb_clear(du); arb_clear(di); arb_clear(hs); arb_clear(hc); arb_clear(ha);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("  -> assoc/eta in [%.4f,%.4f], cstar/eta in [%.4f,%.4f] across t "
           "(bounded, vanishing at t->0)\n",
           ratio_lo, ratio_hi, cratio_lo, cratio_hi);
}

/* ===================== U3: dimension-sweep universality canary =========== */

/* At FIXED t across a WIDE dim sweep dim_A = 4,6,8,10,12,16 (dim_A == d for the
 * dep(d)+conj(dep(d)) family): HOPM_assoc/eta must NOT grow with dim_A. This is the
 * project's single highest-value test (aic-dbo.3, .tex:484: the naive Haar route
 * fails with error ~dim). The sweep is WIDE (FIX B, hostile review R4) so a genuine
 * dim-dependent constant -- linear (~dim) OR sublinear (~sqrt(dim)) -- is AMPLIFIED
 * enough to exceed the threshold. PRINTS per-dim eta so the reader sees eta is
 * ~constant across dims (apples-to-apples: the ratio's dim-dependence is the
 * defect's, not eta's). INTEGRITY: dim_A genuinely varies AND the measured defect
 * varies (so the canary is not vacuously flat). MUTATION: a dim-trivial family
 * (constant dim_A) fires the integrity guard; an injected linear OR sqrt(dim) factor
 * exceeds the tightened threshold. The upper end is dim_A=16: at prec=128 the seven
 * HOPM searches run in a few seconds, and dim_A=16 (a 256x256 superop) amplifies a
 * sqrt(dim) factor by sqrt(16/4)=2.0 -- comfortably over the tightened threshold. */
static void test_u3(void)
{
    const slong P = 128;       /* the ~1e-3 defects resolve fine at prec=128; lower
                                * prec destabilizes the HOPM worst-case search and
                                * INFLATES the measured spread (a false canary trip),
                                * so prec=128 is kept for a clean, reproducible spread */
    const double t = 0.05;
    /* dim_A == d for this family; the sweep spans dim_A=4..16 (a factor-4 range, so
     * a sqrt(dim) pathology is amplified by sqrt(16/4)=2.0). The DOMINANT cost is the
     * Newton-Schulz regularization of the n^2 x n^2 superop (~68s at dim_A=16, prec
     * 128; it scales steeply in n), NOT the HOPM search -- so dim_A=16 is the chosen
     * upper end (the R4 requirement >=16) and the intermediate dim 14 is dropped to
     * bound total runtime; six points still trace a genuine sweep. */
    const slong dims[] = {4, 6, 8, 10, 12, 16};
    /* HOPM restarts/iters: GENEROUS and roughly FLAT across dim (the search, unlike
     * the regularization, is cheap relative to the n^2 Newton-Schulz build). More
     * restarts pin the worst-case sup more tightly, so the real per-dim profile is
     * flatter and more reproducible -- which both lowers the real spread AND lets a
     * sqrt(dim) injection (multiplying the true O(eta) defect) amplify cleanly. */
    const int rs[] = {24, 24, 20, 18, 16, 14}, its[] = {18, 18, 16, 14, 13, 12};
    const int nd = 6;
    double ratio[6];
    slong dimA[6];
    double rmin = 1e300, rmax = 0.0;
    slong dA_min = 1 << 30, dA_max = 0;
    printf("U3 UNIVERSALITY CANARY dep(d)+conj(dep(d)), assoc/eta @ t=%.2f:\n", t);
    for (int i = 0; i < nd; i++) {
        slong d = dims[i];
        aic_ucp_kraus phi;
        make_eta_family(&phi, d, t, P);
        double eta = eta_proxy(&phi, P);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, P);
        arb_t ha;
        arb_init(ha);
        aic_ecstar_defect_assoc_hopm(ha, &ae.A, rs[i], its[i], P);
        double a = dd(ha);
        ratio[i] = a / eta;
        dimA[i] = ae.A.dim_A;
        if (ratio[i] < rmin) rmin = ratio[i];
        if (ratio[i] > rmax) rmax = ratio[i];
        if (dimA[i] < dA_min) dA_min = dimA[i];
        if (dimA[i] > dA_max) dA_max = dimA[i];
        printf("  d=%2ld dim_A=%2ld: eta=%.4e HOPM_assoc=%.4e  assoc/eta=%.4f\n",
               (long) d, (long) ae.A.dim_A, eta, a, ratio[i]);
        AIC_CHECK_MSG(a > 1e-6, "U3: assoc ~0 at d=%ld (search failed / not eta>0)",
                      (long) d);
        arb_clear(ha);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    /* FLATNESS: bounded above (does NOT grow ~ dim) and below (search finds it). */
    AIC_CHECK_MSG(rmin > 1e-3, "U3: assoc/eta fell to ~0 at some dim (%.3e)", rmin);
    /* THRESHOLD 2.0 (FIX B). Measured REAL spread over dim_A=4..16 is 1.42 (flat,
     * dimension-independent; eta is ~constant 3.52e-2..3.76e-2 across the sweep, so
     * the ratio's flatness is the DEFECT's). 2.0 is a modest margin above 1.42 and
     * is mutation-proven to CATCH both pathologies the canary exists to catch
     * (.tex:484, the naive Haar route fails ~dim): injecting a LINEAR factor (dim/4,
     * normalized to 1 at dim_A=4) drives the spread to ~4.67 (RED); a milder
     * sqrt(dim/4) factor drives it to ~2.33 (RED). The factor-4 dim span (4..16) is
     * what gives the sqrt pathology a sqrt(4)=2.0 amplification, enough to clear 2.0
     * given the flat ~1.42 baseline. */
    AIC_CHECK_MSG(rmax / rmin < 2.0,
                  "U3: assoc/eta GREW with dim (max/min=%.2f >= 2.0) -- the ~dim or "
                  "~sqrt(dim) Haar trap (.tex:484 universality FAIL)", rmax / rmin);
    /* INTEGRITY: the canary must be GENUINELY dimension-nontrivial, else it is the
     * toothless flat-by-construction canary the review warned about. The STRUCTURAL
     * witness (which search-count noise CANNOT fake) is that dim_A genuinely SPANS a
     * range: dA_max must strictly exceed dA_min, i.e. the algebra dimension actually
     * grows across the sweep (here dim_A == d, so 4..16). A dim-trivial family
     * (always d=4) makes dA_min==dA_max and fires this. Mutation-proven: forcing
     * make_eta_family to d=4 -> dA range collapses -> RED (the ratio-noise alone
     * would NOT catch it, so this dim_A-span guard is the load-bearing teeth). */
    AIC_CHECK_MSG(dA_max > dA_min,
                  "U3 INTEGRITY: dim_A does NOT vary across the sweep (range [%ld,%ld])"
                  " -- vacuously flat canary", (long) dA_min, (long) dA_max);
    /* and the measured ratio is genuine data (not byte-identical -- a weaker, second
     * witness; the dim_A-span guard above is the primary one). */
    AIC_CHECK_MSG(fabs(ratio[nd - 1] - ratio[0]) > 1e-5 || rmax > rmin + 1e-5,
                  "U3 INTEGRITY: assoc/eta byte-constant across dim (d4=%.4f d16=%.4f)"
                  " -- suspiciously flat", ratio[0], ratio[nd - 1]);
    printf("  -> assoc/eta in [%.4f,%.4f], spread %.2f (flat = dim-independent; "
           "dim_A spans %ld..%ld)\n", rmin, rmax, rmax / rmin,
           (long) dA_min, (long) dA_max);
}

/* ===================== U4: Pi_A accept-guard on non-polar-closed A ======== */

/* The eta>0 assoc A is genuinely NOT polar-closed (only eps-close), so the HOPM
 * block step's Pi_A(polar(C)) + accept-guard CORRECTNESS half is finally exercised
 * (it was structurally untestable on a polar-closed eta=0 A; bead aic-3qq). Verify
 * the assoc-HOPM witness is IN A (proj_residual ~0, certified arb) and
 * lo <= ratio@witness. Also probe that polar(C) of a generic in-A gradient is NOT
 * in A (the algebra is non-polar-closed), so Pi_A is doing real work. */
static void test_u4(void)
{
    const slong P = 256;
    aic_ucp_kraus phi;
    make_eta_family(&phi, 4, 0.06, P);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, P);
    slong n = ae.A.n;

    acb_mat_t wX, wY, wZ;
    acb_mat_init(wX, n, n); acb_mat_init(wY, n, n); acb_mat_init(wZ, n, n);
    arb_t lo;
    arb_init(lo);
    aic_ecstar_defect_assoc_hopm_witness(lo, wX, wY, wZ, &ae.A, 16, 18, P);

    /* (i) witnesses in A (the accept-guard kept every iterate in A). */
    arb_t rX, rY, rZ;
    arb_init(rX); arb_init(rY); arb_init(rZ);
    aic_ecstar_proj_residual(rX, &ae.A, wX, P);
    aic_ecstar_proj_residual(rY, &ae.A, wY, P);
    aic_ecstar_proj_residual(rZ, &ae.A, wZ, P);
    AIC_CHECK_MSG(dd(rX) < 1e-9 && dd(rY) < 1e-9 && dd(rZ) < 1e-9,
                  "U4: witness not in A (residuals %.2e %.2e %.2e)",
                  dd(rX), dd(rY), dd(rZ));

    /* (ii) lo <= ratio at the witness (rigorous lower bound). */
    acb_mat_t xy, lhs, yz, rhs, hh;
    acb_mat_init(xy, n, n); acb_mat_init(lhs, n, n); acb_mat_init(yz, n, n);
    acb_mat_init(rhs, n, n); acb_mat_init(hh, n, n);
    aic_ecstar_star(xy, &ae.A, wX, wY, P);
    aic_ecstar_star(lhs, &ae.A, xy, wZ, P);
    aic_ecstar_star(yz, &ae.A, wY, wZ, P);
    aic_ecstar_star(rhs, &ae.A, wX, yz, P);
    acb_mat_sub(hh, lhs, rhs, P);
    arb_t nh, nx, ny, nz, ratio;
    arb_init(nh); arb_init(nx); arb_init(ny); arb_init(nz); arb_init(ratio);
    aic_mat_opnorm(nh, hh, P);
    aic_mat_opnorm(nx, wX, P);
    aic_mat_opnorm(ny, wY, P);
    aic_mat_opnorm(nz, wZ, P);
    arb_mul(ratio, nx, ny, P);
    arb_mul(ratio, ratio, nz, P);
    arb_div(ratio, nh, ratio, P);
    AIC_CHECK_MSG(arb_le(lo, ratio) || dd(lo) <= dd(ratio) + 1e-9,
                  "U4: lo %.6e > ratio@witness %.6e (NOT a lower bound)",
                  dd(lo), dd(ratio));
    AIC_CHECK_MSG(dd(lo) > 1e-6,
                  "U4: witness ratio ~0 -- the search found no defect (A exactly "
                  "associative? then the accept-guard half is not exercised)");

    /* (iii) NON-POLAR-CLOSED witness: polar(W) of a witness block is generally NOT
     * in A (a real C* algebra is polar-closed; this A is only eps-close). So the
     * Pi_A projection in the block step is doing real work -- WITHOUT it the
     * iterate would leave A. We exhibit this on wX: ||polar(wX) - Pi_A(polar(wX))||
     * is O(eta) > 0. (polar(wX) via the Gram^{-1/2}: U = wX (wX^dag wX)^{-1/2}.) */
    {
        acb_mat_t G, Gisqrt, U, PiU, R;
        acb_mat_init(G, n, n); acb_mat_init(Gisqrt, n, n);
        acb_mat_init(U, n, n); acb_mat_init(PiU, n, n); acb_mat_init(R, n, n);
        acb_mat_t Xd;
        acb_mat_init(Xd, n, n);
        acb_mat_conjugate_transpose(Xd, wX);
        acb_mat_mul(G, Xd, wX, P);                 /* wX^dag wX */
        aic_funcalc_xpow(Gisqrt, G, -0.5, 1.0, P); /* G^{-1/2} */
        acb_mat_mul(U, wX, Gisqrt, P);             /* polar factor U = wX G^{-1/2} */
        /* Pi_A(U) = sum_k <B_k,U>_F B_k */
        acb_mat_zero(PiU);
        acb_t coeff, tt;
        acb_init(coeff); acb_init(tt);
        for (slong k = 0; k < ae.A.dim_A; k++) {
            acb_zero(coeff);
            for (slong p = 0; p < n; p++)
                for (slong q = 0; q < n; q++) {
                    acb_conj(tt, acb_mat_entry(ae.A.B[k], p, q));
                    acb_mul(tt, tt, acb_mat_entry(U, p, q), P);
                    acb_add(coeff, coeff, tt, P);
                }
            for (slong p = 0; p < n; p++)
                for (slong q = 0; q < n; q++) {
                    acb_mul(tt, coeff, acb_mat_entry(ae.A.B[k], p, q), P);
                    acb_add(acb_mat_entry(PiU, p, q), acb_mat_entry(PiU, p, q), tt, P);
                }
        }
        acb_mat_sub(R, U, PiU, P);
        arb_t resid;
        arb_init(resid);
        aic_mat_opnorm(resid, R, P);
        printf("U4 dep(4)+conj(dep4) t=.06: dim_A=%ld lo=%.6e ratio@wit=%.6e "
               "resid(wX,wY,wZ)<%.1e | ||polar(wX)-Pi_A polar(wX)||=%.4e (>0 => "
               "A NON-polar-closed, Pi_A doing real work)\n",
               (long) ae.A.dim_A, dd(lo), dd(ratio),
               fmax(fmax(dd(rX), dd(rY)), dd(rZ)) + 1e-30, dd(resid));
        AIC_CHECK_MSG(dd(resid) > 1e-6,
                      "U4: polar(wX) already in A (resid %.3e) -- A is polar-closed, "
                      "the Pi_A correctness half is NOT exercised (aic-3qq blind spot "
                      "persists)", dd(resid));
        arb_clear(resid);
        acb_clear(tt); acb_clear(coeff);
        acb_mat_clear(Xd);
        acb_mat_clear(R); acb_mat_clear(PiU); acb_mat_clear(U);
        acb_mat_clear(Gisqrt); acb_mat_clear(G);
    }

    arb_clear(ratio); arb_clear(nz); arb_clear(ny); arb_clear(nx); arb_clear(nh);
    acb_mat_clear(hh); acb_mat_clear(rhs); acb_mat_clear(yz);
    acb_mat_clear(lhs); acb_mat_clear(xy);
    arb_clear(rZ); arb_clear(rY); arb_clear(rX);
    arb_clear(lo);
    acb_mat_clear(wZ); acb_mat_clear(wY); acb_mat_clear(wX);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== U5: the projector-SVD spectrum (Fix C) ============= */

/* CONFIRMS the two extraction regimes the docstrings describe (Fix C, hostile
 * review R7). For an idempotent S_tilde the nonzero singular values are >= 1:
 *   - SELF-DUAL channel (HS-self-adjoint Phi_tilde): S_tilde is an ORTHOGONAL
 *     projector, sigma_max == 1 (SVD spectrum {1,..,1,0,..}). dephasing(4) AND --
 *     a FINDING, see below -- the dep(d)+conj(dep(d)) eta>0 family BOTH fall here;
 *   - NON-SELF-DUAL channel: S_tilde is genuinely OBLIQUE, sigma_max STRICTLY > 1,
 *     exercising the >1 extraction path. compress_idemp(4,2) (a non-self-adjoint
 *     compression) gives sigma_max = sqrt(3) ~ 1.732.
 *
 * THE FINDING (corrects a prior docstring claim). The dep(d)+conj(dep(d)) family
 * is a convex mix of dephasing and a UNITARY-conjugate of dephasing; BOTH are
 * HS-self-adjoint, so the mix is HS-self-adjoint, so theta(2 Phi - 1) is too --
 * its S_tilde is an ORTHOGONAL projector (sigma_max == 1), NOT oblique. So that
 * family does NOT exercise the strictly->1 path (it is the U2/U3/U4 NON-associative
 * eta>0 witness, a SEPARATE property from obliqueness). The genuinely-oblique path
 * is exercised by the compression channels (compress_idemp / conj_compress_idemp,
 * which are NOT HS-self-adjoint). Measured sigma_max: dephasing 1.0, dep+conj 1.0,
 * compress_idemp 1.732. */
static void test_u5(void)
{
    const slong P = 256;
    printf("U5 projector-SVD spectrum of S_tilde (Fix C):\n");

    /* (a) self-dual dephasing(4): orthogonal projector, sigma_max == 1. */
    {
        aic_ucp_kraus phi;
        make_dephasing(&phi, 4);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, P);
        double smax = sigma_max_S_tilde(ae.S_tilde);
        printf("  dephasing(4)        sigma_max(S_tilde)=%.6f (orthogonal => ~1)\n",
               smax);
        AIC_CHECK_MSG(fabs(smax - 1.0) < 1e-6,
                      "U5: dephasing S_tilde sigma_max=%.6f != 1 (should be an "
                      "ORTHOGONAL projector, Phi_tilde HS-self-adjoint)", smax);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }

    /* (b) self-dual dep(4)+conj(dep(4)) at eta>0: ALSO an orthogonal projector,
     * sigma_max == 1 (the finding above) -- the U2/U3/U4 non-associativity witness
     * is NOT the obliqueness witness. */
    {
        aic_ucp_kraus phi;
        make_eta_family(&phi, 4, 0.06, P);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, P);
        double smax = sigma_max_S_tilde(ae.S_tilde);
        printf("  dep(4)+conj(dep4)   sigma_max(S_tilde)=%.6f (HS-self-adjoint mix "
               "=> orthogonal => ~1)\n", smax);
        AIC_CHECK_MSG(fabs(smax - 1.0) < 1e-6,
                      "U5: dep+conj S_tilde sigma_max=%.6f != 1 (a mix of HS-self-"
                      "adjoint maps is HS-self-adjoint => orthogonal projector)",
                      smax);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }

    /* (c) NON-self-dual compress_idemp(4,2): genuinely OBLIQUE, sigma_max > 1 --
     * this is the channel that EXERCISES the strictly->1 extraction path. */
    {
        aic_ucp_kraus phi;
        make_compress_idemp(&phi, 4, 2);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, P);
        double smax = sigma_max_S_tilde(ae.S_tilde);
        printf("  compress_idemp(4,2) sigma_max(S_tilde)=%.6f (oblique => >1)\n",
               smax);
        AIC_CHECK_MSG(smax > 1.0 + 1e-4,
                      "U5: compress_idemp S_tilde sigma_max=%.6f not > 1 -- the "
                      "genuinely-OBLIQUE extraction path (non-HS-self-adjoint Phi) "
                      "is NOT exercised", smax);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

int main(void)
{
    test_u1();
    test_u2();
    test_u3();
    test_u4();
    test_u5();
    aic_test_report("test_assoc2");
    printf("OK test_assoc2\n");
    return 0;
}
