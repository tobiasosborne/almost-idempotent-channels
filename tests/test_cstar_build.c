/* test_cstar_build.c — HEADLINE cross-checks for cstar_build Increment 5 (bead
 * aic-097): aic_cstar_build (.tex:1414-1444), the MASTER LOOP / constructive proof
 * of th_main. See include/aic_cstar.h (Increment 5 banner), src/aic_cstar_build.c
 * (driver + Stage 1), src/aic_cstar_stages.c (Stage 2 + Stage 3),
 * src/aic_cstar_classes.c (classes + wrapper-unit errreduce), and the spec
 * docs/research/cstar_masterloop_spec.md §8.
 *
 * THE THREE HEADLINE CROSS-CHECKS (Rule 6 — cross-checks > unit tests):
 *   T1 — eta=0 ORACLE (the cleanest ground truth). On EXACT-idempotent multi-block
 *        channels, aic_cstar_build must produce B whose num_blocks + d[] MATCH the
 *        INDEPENDENT aic_idemp_decompose (th_idemp_structure) block sizes EXACTLY,
 *        with iso_def ~ machine-zero and v bijective. aic_idemp_decompose is a
 *        completely separate code path (different algorithm), so this is a genuine
 *        cross-check, not a self-consistency test.
 *   T2 — UNIVERSALITY CANARY (.tex:484, the single most important property). The
 *        th_main constant is DIMENSION-INDEPENDENT. T2a: at eta=0 the construction is
 *        EXACT at every dim_A in {8,9,18,20} — iso_def stays at machine-zero, i.e. the
 *        constant is literally 0 and FLAT across a 2.5x dim range (the strongest
 *        possible dim-independence statement, and the direct contradiction of the
 *        .tex:484 error-proportional-to-n failure mode). T2b: on the oblique eta>0
 *        path, an EXPANDED ambient-n sweep (m=2: n=4..10, m=3: n=6..9) feeds a ROBUST
 *        bounded+no-trend canary — abs-max c<5 AND a thirds-ratio mean(c|hi-n)/
 *        mean(c|lo-n)<=1.4 — that dilutes the n=7 hard-geometry spike yet trips on a
 *        genuine c=O(n) law (mutation-proven). This REPLACES the geometry-fragile
 *        two-point within-family ratio metric (FINDINGS §D2/§C11).
 *   T3 — OBLIQUE eta>0. The master loop on a GENUINELY oblique (non-self-adjoint)
 *        associated algebra: iso_def is certified O(eta) via the §C8 c=iso_def/eta
 *        MAGNITUDE bound, and the star->plain mutation tooth (FINDINGS §C8) is
 *        confirmed RED by hand (recorded in the increment report, not left failing).
 *
 * THE STAR IS LOAD-BEARING (FINDINGS §C2/§C8). The whole loop's defects route through
 * A's Choi-Effros STAR; the eta=0 oracle has star==plain and is BLIND to a star bug,
 * so the oblique T3 fixture (make_mixconj) carries the star teeth.
 *
 * A SUBSTRATE WALL surfaced here (the increment report escalates it): the Stage-2
 * oblique-wrapper lem_extension path drives aic_sgn (Newton-Schulz) on near-boundary
 * corner matrices that FAIL to converge for ambient n >= 6 / block m >= 3 at the t
 * regimes tried. T2b/T3 therefore use the n in {4,5}, m=2 oblique points that
 * complete; the eta=0 dim sweep (T2a) is the robust, in-basin dim-independence proof.
 *
 * TEETH. Tooth (a) is a LIVE automated assertion: T3 re-sweeps the final v with A's
 * star mutated to the PLAIN product and asserts c_plain > 20 (measured ~40-65, a
 * ~30x gap from c_star ~ 1-1.8; FINDINGS §C8). The remaining teeth are mutation-
 * proven BY HAND (Rule 7; each RED observed, then restored — recorded in the
 * increment report, NOT left failing in the suite):
 *   (b) a wrong B block size (force one fewer class) -> sigma_min collapse / build
 *       bijectivity abort, T1 RED;
 *   (c) skip an errreduce call -> c_0 drift gate / lem_extension input too large;
 *   (d) the §C3 Frobenius-Pi_A-instead-of-Co projection -> projection defect O(1) RED.
 */
#include <math.h>
#include <stdio.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_errreduce.h"
#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"
#include "test_idemp.h"

#define CPREC 256

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Suppress -Wunused-function for the test_idemp.h builders this TU does not use. */
__attribute__((unused))
static void suppress_unused_idemp_builders(void)
{
    (void) make_trace_replace;
    (void) make_identity;
    (void) make_compress_idemp;
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op via the §C5 midpoint route (mirrors
 * test_cstar_extension.c eta_proxy). */
static double eta_proxy(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S, S2, D, M;
    acb_mat_init(S, nn, nn);
    acb_mat_init(S2, nn, nn);
    acb_mat_init(D, nn, nn);
    acb_mat_init(M, nn, nn);
    aic_assoc_superop_from_ucp(S, phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    for (slong i = 0; i < nn; i++)
        for (slong j = 0; j < nn; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    double o = dd(e);
    arb_clear(e);
    acb_mat_clear(M); acb_mat_clear(D); acb_mat_clear(S2); acb_mat_clear(S);
    return o;
}

/* Sorted multiset of B's block sizes (for the T1 set-comparison vs the oracle). */
static int dbl_cmp(const void *p, const void *q)
{
    long a = *(const long *) p, b = *(const long *) q;
    return (a > b) - (a < b);
}

/* ---- T2b dimension-trend statistics (the robust .tex:484 canary, FINDINGS §D2) -
 * Both operate on parallel arrays nv[] (ambient n as double) and cv[] (the measured
 * c=iso_def/eta). They AGGREGATE over the whole sweep, so a single hard-geometry
 * outlier (e.g. the n=7 Heisenberg-Weyl compression peak, FINDINGS §C11) is diluted,
 * but a systematic c proportional-to-n (or sqrt-n) growth — the genuine .tex:484
 * failure mode — survives the aggregation and trips the bound. */

/* Least-squares slope beta of c vs n: beta = Cov(n,c)/Var(n). beta<=0 means c does
 * not climb with n on average; a c=O(n) law forces beta to scale with the per-step
 * growth rate. Robust to one interior spike (least squares averages all points). */
static double trend_slope(const double *nv, const double *cv, int k)
{
    double nm = 0.0, cm = 0.0;
    for (int i = 0; i < k; i++) { nm += nv[i]; cm += cv[i]; }
    nm /= k;
    cm /= k;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < k; i++) {
        num += (nv[i] - nm) * (cv[i] - cm);
        den += (nv[i] - nm) * (nv[i] - nm);
    }
    return den > 0.0 ? num / den : 0.0;
}

/* Aggregate upper-half / lower-half mean ratio of c, points taken IN n-ORDER (the
 * caller supplies them sorted by n). Splits the sweep into a low-n half and a high-n
 * half (floor(k/2) points each; an odd middle point — e.g. the n=7 spike on the
 * 7-point m=2 sweep — is in NEITHER half) and returns mean(high)/mean(low). A
 * bounded, non-growing c gives a ratio near 1 or below; a c=O(n) growth gives a
 * ratio that scales with the dim ratio of the two halves. Averaging each half
 * dilutes a single hard-geometry outlier (FINDINGS §C11). */
static double trend_halves_ratio(const double *cv, int k)
{
    int t = k / 2;
    if (t < 1) t = 1;
    double lo = 0.0, hi = 0.0;
    for (int i = 0; i < t; i++) lo += cv[i];
    for (int i = 0; i < t; i++) hi += cv[k - 1 - i];
    lo /= t;
    hi /= t;
    return lo > 1e-12 ? hi / lo : 0.0;
}

/* MULTI-CLASS oblique fixture (T4): mix make_block_cond_exp(d,m) [M_m (+) M_{d-m},
 * TWO genuine matrix blocks] with its FIXED-unitary conjugate, the same Kraus-union
 * recipe as make_mixconj (test_idemp.h) but with a TWO-BLOCK base. Unlike
 * make_mixconj (whose make_compress_idemp base is a SINGLE M_m block, so every
 * eta>0 mixconj is one equivalence class), this yields >=2 equivalence classes at
 * eta>0 — the ONLY fixture in the suite that exercises the Stage-3 MULTI-CLASS merge
 * + the errreduce_unit Stage-3 running-P_total branch at eta>0 (the I5 coverage gap,
 * FINDINGS §C11). CAVEAT (measured): its associativity defect eps_assoc ~ 2e-5 is
 * ~700x BELOW eta ~ 1.6e-2 (the two blocks stay well-separated under the global
 * conjugation), so the build's eps argument must be a faithful O(eta) scale (we pass
 * eta), NOT the assoc defect (which would make the Stage-1 errreduce C0 gate fire,
 * 10*eps < the true O(eta) inclusion defect). Local to this TU (not test_idemp.h:
 * no other TU needs it, and adding it there would warn -Wunused in 13 files). */
static void make_mixconj_blocks(aic_ucp_kraus *out, slong d, slong m, double t,
                                slong prec)
{
    aic_ucp_kraus base, conj;
    make_block_cond_exp(&base, d, m);
    make_conjugated(&conj, &base, prec);
    slong n = base.dim_H;
    aic_ucp_kraus_init(out, n, n, base.r + conj.r);
    arb_t s;
    arb_init(s);
    arb_set_d(s, sqrt(1.0 - t));
    for (slong a = 0; a < base.r; a++)
        acb_mat_scalar_mul_arb(out->K[a], base.K[a], s, prec);
    arb_set_d(s, sqrt(t));
    for (slong b = 0; b < conj.r; b++)
        acb_mat_scalar_mul_arb(out->K[base.r + b], conj.K[b], s, prec);
    arb_clear(s);
    aic_ucp_kraus_clear(&conj);
    aic_ucp_kraus_clear(&base);
}

/* ===================== T1: eta=0 oracle ===================================== *
 * On an EXACT idempotent channel, run aic_cstar_build and assert B's block
 * structure (num_blocks + sorted d[]) matches the EXPECTED th_idemp_structure
 * blocks (and dim_A == aic_idemp_decompose's dim_A, the INDEPENDENT cross-check),
 * iso_def ~ machine-zero, v bijective. */
static void t1_one(const char *name, aic_ucp_kraus *phi, slong exp_blocks,
                   const long *exp_d)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);

    /* Independent oracle: th_idemp_structure dim_A. */
    aic_idemp_decomp dc;
    aic_idemp_decompose(&dc, phi, CPREC);

    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);

    /* sort B's block sizes for the multiset comparison */
    long got[16];
    for (slong l = 0; l < B.num_blocks; l++) got[l] = (long) B.d[l];
    qsort(got, (size_t) B.num_blocks, sizeof(long), dbl_cmp);
    long want[16];
    for (slong l = 0; l < exp_blocks; l++) want[l] = exp_d[l];
    qsort(want, (size_t) exp_blocks, sizeof(long), dbl_cmp);

    printf("T1 %s: idemp dim_A=%ld | build B num_blocks=%ld d=[", name,
           (long) dc.dim_A, (long) B.num_blocks);
    for (slong l = 0; l < B.num_blocks; l++)
        printf("%ld%s", got[l], l + 1 < B.num_blocks ? "," : "");
    printf("] dim_B=%ld iso_def=%.3e\n", (long) B.dim_B, dd(iso));

    AIC_CHECK_MSG(B.dim_B == ae.A.dim_A,
                  "T1 %s: dim_B=%ld != dim_A=%ld", name, (long) B.dim_B,
                  (long) ae.A.dim_A);
    AIC_CHECK_MSG(ae.A.dim_A == dc.dim_A,
                  "T1 %s: ecstar dim_A=%ld != idemp dim_A=%ld (oracle)", name,
                  (long) ae.A.dim_A, (long) dc.dim_A);
    AIC_CHECK_MSG(B.num_blocks == exp_blocks,
                  "T1 %s: num_blocks=%ld != expected %ld", name,
                  (long) B.num_blocks, (long) exp_blocks);
    for (slong l = 0; l < exp_blocks; l++)
        AIC_CHECK_MSG(got[l] == want[l],
                      "T1 %s: block-size multiset mismatch at %ld (%ld != %ld)",
                      name, (long) l, got[l], want[l]);

    AIC_CHECK_MSG(dd(iso) < 1e-10,
                  "T1 %s: eta=0 iso_def=%.3e not ~0", name, dd(iso));

    arb_t a;
    arb_init(a);
    int bij = aic_errreduce_is_bijective(a, &v, CPREC);
    AIC_CHECK_MSG(bij, "T1 %s: v not bijective", name);
    AIC_CHECK_MSG(dd(a) > 0.5, "T1 %s: sigma_min=%.4f <= 0.5", name, dd(a));
    arb_clear(a);

    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_idemp_clear(&dc);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t1_eta0_oracle(void)
{
    printf("T1 eta=0 oracle (B matches aic_idemp_decompose block sizes EXACTLY):\n");
    aic_ucp_kraus phi;

    long d_22[2] = {2, 2};
    make_block_cond_exp(&phi, 4, 2);      /* M_2 (+) M_2, dim_A=8 */
    t1_one("block_cond_exp(4,2)", &phi, 2, d_22);
    aic_ucp_kraus_clear(&phi);

    long d_24[2] = {2, 4};
    make_block_cond_exp(&phi, 6, 2);      /* M_2 (+) M_4, dim_A=20 */
    t1_one("block_cond_exp(6,2)", &phi, 2, d_24);
    aic_ucp_kraus_clear(&phi);

    long d_33[2] = {3, 3};
    make_block_cond_exp(&phi, 6, 3);      /* M_3 (+) M_3, dim_A=18 */
    t1_one("block_cond_exp(6,3)", &phi, 2, d_33);
    aic_ucp_kraus_clear(&phi);

    long d_3[1] = {3};
    make_noiseless_subsystem(&phi, 3, 2); /* M_3, dim_A=9 */
    t1_one("noiseless_subsystem(3,2)", &phi, 1, d_3);
    aic_ucp_kraus_clear(&phi);

    long d_1111[4] = {1, 1, 1, 1};
    make_dephasing(&phi, 4);              /* C^4 = M_1^4, dim_A=4 */
    t1_one("dephasing(4)", &phi, 4, d_1111);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T2: universality canary (.tex:484) ================== *
 * T2a — the dimension-INDEPENDENCE proof at eta=0: sweep dim_A in {8,9,18,20} and
 * assert iso_def stays at MACHINE-ZERO at every dim (the constant is 0, FLAT across
 * a 2.5x dim range — the direct contradiction of the .tex:484 error-proportional-
 * to-n failure mode). MUTATION (by hand, recorded): a `iso_def += 0.01*sqrt(dim_B)`
 * additive injection in the final certificate would make iso_def GROW with dim ->
 * the < 1e-8 gate RED at the larger dims (a clean dim-growth detector). */
static void test_t2_universality(void)
{
    printf("T2a universality canary (.tex:484): eta=0 iso_def FLAT (=0) over dim_A "
           "sweep:\n");
    struct { slong d, m; const char *nm; } fam[] = {
        {4, 2, "bc(4,2)"}, {6, 3, "bc(6,3)"}, {6, 2, "bc(6,2)"}};
    double iso_max = 0.0;
    slong dimA_min = 1 << 20, dimA_max = 0;
    for (int i = 0; i < 3; i++) {
        aic_ucp_kraus phi;
        make_block_cond_exp(&phi, fam[i].d, fam[i].m);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);
        aic_dhom_B B;
        aic_dhom_v v;
        arb_t iso;
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);
        double isod = dd(iso);
        if (isod > iso_max) iso_max = isod;
        if (ae.A.dim_A < dimA_min) dimA_min = ae.A.dim_A;
        if (ae.A.dim_A > dimA_max) dimA_max = ae.A.dim_A;
        printf("  %s: dim_A=%ld dim_B=%ld iso_def=%.3e (must stay ~0)\n", fam[i].nm,
               (long) ae.A.dim_A, (long) B.dim_B, isod);
        AIC_CHECK_MSG(dd(iso) < 1e-8,
                      "T2a %s: iso_def=%.3e grows with dim (.tex:484 failure mode!)",
                      fam[i].nm, isod);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    /* one noiseless point too (dim_A=9) */
    {
        aic_ucp_kraus phi;
        make_noiseless_subsystem(&phi, 3, 2);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);
        aic_dhom_B B;
        aic_dhom_v v;
        arb_t iso;
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);
        double isod = dd(iso);
        if (isod > iso_max) iso_max = isod;
        printf("  ns(3,2): dim_A=%ld iso_def=%.3e\n", (long) ae.A.dim_A, isod);
        AIC_CHECK_MSG(dd(iso) < 1e-8, "T2a ns(3,2): iso_def=%.3e not ~0", isod);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("  -> iso_def_max=%.3e over dim_A in [%ld,%ld] (%.1fx range): FLAT at 0, "
           "the constant does NOT grow with dim (.tex:461)\n", iso_max,
           (long) dimA_min, (long) dimA_max, (double) dimA_max / (double) dimA_min);

    /* T2b — the eta>0 measured constant c=iso_def/eta is dimension-INDEPENDENT
     * (.tex:484 / .tex:461: "the implicit constant in O(eps) does not depend on A or
     * its dimensionality"). The two-point WITHIN-FAMILY ratio c_hi/c_lo that this
     * test used on commit cb273f2 was FLAWED: it measured fixture-GEOMETRY SPREAD
     * across different ambient n, NOT dimension-GROWTH. An extended precision-stable
     * sweep (orchestrator diagnostic, prec 256 == prec 512 byte-for-byte) established
     * that c is BOUNDED in [0.25, 3.27] and does NOT grow with n: BOTH the m=2 and
     * m=3 families PEAK at the n=7 Heisenberg-Weyl/compression geometry (a hard
     * corner, FINDINGS §C11) then CRASH at n=8 (the LARGEST-n points have the
     * SMALLEST c) — fixture-geometry noise, NOT the .tex:484 c-proportional-to-n
     * failure mode. All isos stay healthy (bijective, sigma_min in [0.96,0.999]).
     * The old ratio inflates to 6.90 on m=3 purely because a hard-geometry outlier
     * (n=7, c=1.87) sits over a favorable one (n=8, c=0.25), unrelated to .tex:484.
     *
     * The ROBUST canary (FINDINGS §D2): an EXPANDED sweep — m=2 family (B=M_2) over
     * n=4..10, m=3 family (B=M_3) over n=6..10 — fed to TWO aggregate checks that
     * dilute any single hard-geometry outlier yet still catch a systematic c=O(n)
     * (or O(sqrt n)) growth:
     *   (i)  ABSOLUTE boundedness: abs-max c over the whole sweep < C_ABS=5.0 (the
     *        measured max is 3.27 at the n=7 m=2 peak; 1.53x margin). A c growing
     *        with n eventually exceeds any fixed bound.
     *   (ii) NO upward trend with n: the halves-ratio mean(c|upper-half n) /
     *        mean(c|lower-half n) <= KAPPA=1.25. Averaging each half dilutes the n=7
     *        spike (it is the odd MIDDLE point on the 7-point m=2 sweep, in NEITHER
     *        half; on m=3 it averages into the lower half); a c=O(n) law drives the
     *        ratio to the dim ratio of the halves. Measured: m=2 reads ~0.28, m=3
     *        reads ~0.65 — both DEEP under 1.25. The least-squares slope beta of c vs
     *        n is reported as a corroborating diagnostic (measured beta ~ -0.21 for
     *        m=2, ~-0.08 for m=3 — NEGATIVE, c trends DOWN over the full sweep).
     * The old geometry-fragile two-point within-family extremes ratio is GONE; it
     * inflated to 6.90 on the n=6,7 m=3 endpoints purely by spike-over-favorable
     * placement, with NO dim-growth content (FINDINGS §D2/§C11). */
    printf("T2b universality canary (.tex:484/.tex:461): eta>0 c=iso_def/eta bounded "
           "+ NO upward trend with n:\n");
    const double C_ABS = 5.0;    /* abs-max c bound (data max 3.27; 1.53x margin)   */
    const double KAPPA = 1.25;   /* halves-ratio bound (catches a genuine c=O(n))   */
    struct { slong d, m; double t; int fam; } ob[] = {
        {4, 2, 0.03, 0}, {5, 2, 0.03, 0}, {6, 2, 0.03, 0}, {7, 2, 0.03, 0},
        {8, 2, 0.03, 0}, {9, 2, 0.03, 0}, {10, 2, 0.03, 0}, /* m=2: B=M_2, n=4..10 */
        {6, 3, 0.02, 1}, {7, 3, 0.02, 1}, {8, 3, 0.02, 1},
        {9, 3, 0.02, 1}, {10, 3, 0.02, 1}                   /* m=3: B=M_3, n=6..10 */
    };
    const int n_ob = (int) (sizeof(ob) / sizeof(ob[0]));
    /* per-family sweeps kept in n-ORDER (ob[] is already sorted by n within fam) */
    double nv[2][16], cv[2][16];
    int kf[2] = {0, 0};
    double c_abs_max = 0.0;
    for (int i = 0; i < n_ob; i++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, ob[i].d, ob[i].m, ob[i].t, CPREC);
        double eta = eta_proxy(&phi, CPREC);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);
        arb_t ea;
        arb_init(ea);
        aic_ecstar_defect_assoc(ea, &ae.A, CPREC);
        double eps = dd(ea);
        arb_clear(ea);
        aic_dhom_B B;
        aic_dhom_v v;
        arb_t iso;
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, eps, CPREC);
        double c = eta > 1e-12 ? dd(iso) / eta : 0.0;
        int f = ob[i].fam;
        nv[f][kf[f]] = (double) ob[i].d;
        cv[f][kf[f]] = c;
        kf[f]++;
        if (c > c_abs_max) c_abs_max = c;
        printf("  mixconj(%ld,%ld,%.2f): n=%ld dim_B=%ld eta=%.3e iso_def=%.3e "
               "c=iso/eta=%.4f\n", (long) ob[i].d, (long) ob[i].m, ob[i].t,
               (long) ae.A.n, (long) B.dim_B, eta, dd(iso), c);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    double slope2 = trend_slope(nv[0], cv[0], kf[0]);
    double slope3 = trend_slope(nv[1], cv[1], kf[1]);
    double ratio2 = trend_halves_ratio(cv[0], kf[0]);
    double ratio3 = trend_halves_ratio(cv[1], kf[1]);
    printf("  -> m=2 (n=4..10): slope=%.4f halves-ratio(hi/lo)=%.4f | m=3 (n=6..10): "
           "slope=%.4f halves-ratio=%.4f | abs-max c=%.4f\n",
           slope2, ratio2, slope3, ratio3, c_abs_max);
    /* (i) ABSOLUTE boundedness — dimension-independent (.tex:461). */
    AIC_CHECK_MSG(c_abs_max < C_ABS,
                  "T2b: abs-max c=%.4f >= C_ABS=%.1f over the sweep (the th_main "
                  "constant grows with dim; .tex:461/.tex:484 / FINDINGS §D2 stop "
                  "condition)", c_abs_max, C_ABS);
    /* (ii) NO upward trend with n — the robust anti-.tex:484 check. The halves-ratio
     * dilutes the n=7 hard-geometry spike (FINDINGS §C11); a c=O(n) growth survives
     * the aggregation and trips KAPPA (mutation-proven below). */
    AIC_CHECK_MSG(ratio2 <= KAPPA,
                  "T2b: m=2 halves-ratio %.4f > KAPPA=%.2f (c TRENDS UP with n; "
                  ".tex:484 / FINDINGS §D2 stop condition)", ratio2, KAPPA);
    AIC_CHECK_MSG(ratio3 <= KAPPA,
                  "T2b: m=3 halves-ratio %.4f > KAPPA=%.2f (c TRENDS UP with n; "
                  ".tex:484 / FINDINGS §D2 stop condition)", ratio3, KAPPA);

    /* MUTATION TOOTH (Rule 5/7) — prove the canary has teeth against the genuine
     * .tex:484 failure mode. NO injection is active in the path above; this re-runs
     * the metric on injected m=2 data. Two failure models:
     *  (A) literal c_inj = c * (n/n_min): amplifies the REAL (noisy) c by n/n_min.
     *      Because the real c CRASHES at n>=8 (geometry), this stays non-monotone, so
     *      the TREND arm need not fire (its inherited crash deflates the upper half) —
     *      but the ABSOLUTE arm DOES: the amplified n=7 peak blows past C_ABS
     *      (measured abs-max ~5.73 > 5.0). So the canary AS A WHOLE catches the
     *      task-specified injection (via the absolute arm).
     *  (B) faithful c_inj = c_flat * (n/n_min), c_flat=mean(c): a clean c=O(n) law
     *      (the constant ITSELF scaling with dim — exactly .tex:484's failure mode,
     *      and the reason the TREND arm exists). Monotone in n, so the TREND arm
     *      fires: halves-ratio ~1.80 > KAPPA=1.25, slope ~+0.36 (vs the real ~-0.21).
     * Together: the canary catches the literal injection (abs arm) AND a genuine
     * c=O(n) law (trend arm). */
    {
        double cinjA[16] = {0}, cinjB[16] = {0}, cflat = 0.0;
        double nmin = nv[0][0];
        for (int i = 0; i < kf[0]; i++) {
            cflat += cv[0][i];
            if (nv[0][i] < nmin) nmin = nv[0][i];
        }
        cflat /= kf[0];
        double maxA = 0.0;
        for (int i = 0; i < kf[0]; i++) {
            cinjA[i] = cv[0][i] * (nv[0][i] / nmin);   /* model A: literal */
            cinjB[i] = cflat * (nv[0][i] / nmin);      /* model B: faithful O(n) */
            if (cinjA[i] > maxA) maxA = cinjA[i];
        }
        double ratioB = trend_halves_ratio(cinjB, kf[0]);
        double slopeB = trend_slope(nv[0], cinjB, kf[0]);
        printf("  [mutation tooth] inj-A c*(n/nmin): abs-max=%.4f (vs C_ABS=%.1f) | "
               "inj-B c_flat*(n/nmin): halves-ratio=%.4f slope=%.4f (vs KAPPA=%.2f)\n",
               maxA, C_ABS, ratioB, slopeB, KAPPA);
        /* the canary AS A WHOLE fires on model A (absolute arm) ... */
        AIC_CHECK_MSG(maxA >= C_ABS,
                      "T2b mutation: inj-A abs-max=%.4f did NOT exceed C_ABS=%.1f — "
                      "the ABSOLUTE arm is BLIND to c*(n/nmin) growth", maxA, C_ABS);
        /* ... and the TREND arm fires on the faithful c=O(n) model B. */
        AIC_CHECK_MSG(ratioB > KAPPA,
                      "T2b mutation: inj-B halves-ratio=%.4f did NOT exceed KAPPA=%.2f "
                      "— the TREND arm is BLIND to a genuine c=O(n) law", ratioB, KAPPA);
    }
}

/* plain-product star thunk (the §C8 mutation): out = XY (no Phi_tilde). */
static void plain_star(acb_mat_t out, const acb_mat_t XY, void *ctx, slong prec)
{
    (void) ctx;
    (void) prec;
    acb_mat_set(out, XY);
}

/* ===================== T3: oblique eta>0 + §C8 star tooth ================== *
 * The master loop on a GENUINELY oblique associated algebra (make_mixconj, S_tilde
 * non-self-adjoint, sigma_max(S~) > 1). iso_def is certified O(eta) via the §C8
 * c=iso_def/eta MAGNITUDE bound. The star tooth: re-sweep the final v with A's star
 * mutated to the PLAIN product; c_plain must blow past the bound (a ~weak DIRECTION
 * tooth on in-A projections, but the MAGNITUDE gap is sharp, FINDINGS §C8). The
 * automated suite keeps the star path; the plain mutation is recorded RED by hand.
 *
 * The m=3 points (mixconj(6,3), mixconj(7,3)) exercise the DEEPER Stage-2
 * induction (M_1 -> M_2 -> M_3, two lem_extension steps on a dim S_{P,Q}=2 corner)
 * that the I4 fix (eps_target=O(eta), generous unit_tol for the G-twisted Ha
 * codomain; FINDINGS §C11) unblocked. Each make_mixconj is a SINGLE equivalence
 * class M_m (make_compress_idemp's image is one matrix block), so num_blocks==1 and
 * B.d[0]==m. (The Stage-3 MULTI-CLASS merge at eta>0 is covered by T4.) */
static void test_t3_oblique(void)
{
    printf("T3 oblique eta>0 master loop (mixconj(d,m)); §C8 c=iso_def/eta O(eta) "
           "bound + star tooth:\n");
    struct { slong d, m; double t; } fx[] = {
        {4, 2, 0.03}, {4, 2, 0.02}, {5, 2, 0.03}, /* m=2 (single M_2)       */
        {6, 3, 0.02}, {6, 3, 0.03}                /* m=3 (single M_3, I4)   */
    };
    const int n_fx = (int) (sizeof(fx) / sizeof(fx[0]));
    double cstar_rec[8], cplain_rec[8];
    for (int i = 0; i < n_fx; i++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, fx[i].d, fx[i].m, fx[i].t, CPREC);
        double eta = eta_proxy(&phi, CPREC);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);
        arb_t ea;
        arb_init(ea);
        aic_ecstar_defect_assoc(ea, &ae.A, CPREC);
        double eps = dd(ea);
        arb_clear(ea);

        aic_dhom_B B;
        aic_dhom_v v;
        arb_t iso;
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, eps, CPREC);
        double c_star = eta > 1e-12 ? dd(iso) / eta : 0.0;

        /* §C8 star tooth: re-sweep the FINAL v's multiplicativity defect with A's
         * star swapped to the plain product. */
        aic_ecstar_apply_fn save_phi = ae.A.star_phi;
        void *save_ctx = ae.A.star_ctx;
        ae.A.star_phi = plain_star;
        ae.A.star_ctx = NULL;
        arb_t mdp;
        arb_init(mdp);
        aic_dhom_defect_sweep(mdp, &v, CPREC);
        double c_plain = eta > 1e-12 ? dd(mdp) / eta : 0.0;
        ae.A.star_phi = save_phi;
        ae.A.star_ctx = save_ctx;
        arb_clear(mdp);

        cstar_rec[i] = c_star;
        cplain_rec[i] = c_plain;
        printf("  mixconj(%ld,%ld,%.2f): n=%ld dim_A=%ld eta=%.3e | iso_def=%.3e "
               "c_star=%.4f | PLAIN c=%.4f | num_blocks=%ld d[0]=%ld\n",
               (long) fx[i].d, (long) fx[i].m, fx[i].t, (long) ae.A.n,
               (long) ae.A.dim_A, eta, dd(iso), c_star, c_plain,
               (long) B.num_blocks, (long) B.d[0]);

        /* (a) iso_def is O(eta): a generous magnitude bound (true c ~ 0.3-1.8). */
        AIC_CHECK_MSG(c_star < 20.0,
                      "T3: iso_def/eta=%.4f exceeds 20 (master loop not O(eta)?)",
                      c_star);
        /* the §C8 star tooth: the PLAIN-product mutation must blow past the bound
         * (here c_plain ~ 40-72 for both m=2 and m=3, a large gap from c_star ~
         * 0.3-1.8); the LIVE star-vs-plain MAGNITUDE discriminant, FINDINGS §C8. */
        AIC_CHECK_MSG(c_plain > 20.0,
                      "T3: PLAIN-product c=%.4f not > 20 — the §C8 star tooth is "
                      "BLIND (mixconj(%ld,%ld,%.2f))", c_plain, (long) fx[i].d,
                      (long) fx[i].m, fx[i].t);
        /* (c) single class M_m: num_blocks==1, block size == m. */
        AIC_CHECK_MSG(B.num_blocks == 1,
                      "T3: mixconj(%ld,%ld) should give a single M_%ld block, got "
                      "num_blocks=%ld", (long) fx[i].d, (long) fx[i].m,
                      (long) fx[i].m, (long) B.num_blocks);
        AIC_CHECK_MSG(B.d[0] == fx[i].m,
                      "T3: mixconj(%ld,%ld) block size %ld != m=%ld",
                      (long) fx[i].d, (long) fx[i].m, (long) B.d[0], (long) fx[i].m);
        /* (b) bijective */
        arb_t a;
        arb_init(a);
        AIC_CHECK_MSG(aic_errreduce_is_bijective(a, &v, CPREC), "T3: v not bijective");
        arb_clear(a);

        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("  -> STAR c:");
    for (int i = 0; i < n_fx; i++) printf(" %.4f", cstar_rec[i]);
    printf(" (all < 20, O(eta)); PLAIN-mutation c (live tooth):");
    for (int i = 0; i < n_fx; i++) printf(" %.4f", cplain_rec[i]);
    printf("\n");
}

/* ===================== T4: MULTI-CLASS oblique merge at eta>0 ============== *
 * The I5 hostile-review coverage gap (FINDINGS §C11): every prior eta>0 fixture was
 * a SINGLE equivalence class, so the Stage-3 MULTI-CLASS merge (cor_merge_sum across
 * classes + the errreduce_unit Stage-3 running-P_total!=1_n branch) ran only at
 * eta=0 (machine-zero defects). make_mixconj_blocks (block_cond_exp(d,m) base ->
 * M_m (+) M_{d-m}, conjugate-mixed) is a genuinely oblique eta>0 channel with >=2
 * equivalence classes, so it DRIVES Stage-3 at eta>0. ASSERT: build COMPLETES,
 * num_blocks==2 with the expected sizes, iso_def/eta bounded (O(eta)), v bijective.
 *
 * NOTE on the eps argument (measured): this fixture's associativity defect is ~700x
 * BELOW eta (the two blocks stay well-separated under the global conjugation), so we
 * pass eta as the build's eps (a faithful O(eta) scale); the assoc defect would make
 * the Stage-1 errreduce C0 gate fire. NO star tooth here: this fixture is nearly
 * block-diagonal so its PLAIN-product defect is also small (c_plain ~ 0.3, does NOT
 * fire the >20 tooth) — the star magnitude discriminant lives on the single-block
 * mixconj fixtures (T3), and is preserved there. */
static void test_t4_multiclass(void)
{
    printf("T4 multi-class oblique merge at eta>0 (mixconj_blocks; Stage-3 + "
           "errreduce_unit running-unit branch at eta>0):\n");
    struct { slong d, m; double t; slong exp_b0, exp_b1; } fx[] = {
        {4, 2, 0.02, 2, 2}, /* M_2 (+) M_2 */
        {5, 2, 0.02, 2, 3}, /* M_2 (+) M_3 */
        {4, 2, 0.03, 2, 2}  /* M_2 (+) M_2, larger t */
    };
    const int n_fx = (int) (sizeof(fx) / sizeof(fx[0]));
    for (int i = 0; i < n_fx; i++) {
        aic_ucp_kraus phi;
        make_mixconj_blocks(&phi, fx[i].d, fx[i].m, fx[i].t, CPREC);
        double eta = eta_proxy(&phi, CPREC);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);
        /* eps = faithful O(eta) scale (see the test docstring; the assoc defect
         * badly underestimates eta for this near-block-diagonal fixture). */
        double eps = eta;
        aic_dhom_B B;
        aic_dhom_v v;
        arb_t iso;
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, eps, CPREC);
        double c = eta > 1e-12 ? dd(iso) / eta : 0.0;

        /* sort block sizes for the multiset comparison */
        long got[16];
        for (slong l = 0; l < B.num_blocks; l++) got[l] = (long) B.d[l];
        qsort(got, (size_t) B.num_blocks, sizeof(long), dbl_cmp);
        long want[2] = {(long) fx[i].exp_b0, (long) fx[i].exp_b1};
        qsort(want, 2, sizeof(long), dbl_cmp);

        arb_t a;
        arb_init(a);
        int bij = aic_errreduce_is_bijective(a, &v, CPREC);
        printf("  mixconj_blocks(%ld,%ld,%.2f): n=%ld eta=%.3e iso_def=%.3e "
               "c=iso/eta=%.4f num_blocks=%ld d=[%ld,%ld] bij=%d sig=%.4f\n",
               (long) fx[i].d, (long) fx[i].m, fx[i].t, (long) ae.A.n, eta, dd(iso),
               c, (long) B.num_blocks, got[0],
               B.num_blocks > 1 ? got[1] : 0, bij, dd(a));

        /* (a) the Stage-3 multi-class merge ran: >=2 classes. */
        AIC_CHECK_MSG(B.num_blocks == 2,
                      "T4: mixconj_blocks(%ld,%ld) should give 2 classes, got %ld "
                      "(Stage-3 multi-class merge not exercised at eta>0)",
                      (long) fx[i].d, (long) fx[i].m, (long) B.num_blocks);
        AIC_CHECK_MSG(got[0] == want[0] && got[1] == want[1],
                      "T4: block sizes [%ld,%ld] != expected [%ld,%ld]", got[0],
                      got[1], want[0], want[1]);
        /* (b) iso_def is O(eta). */
        AIC_CHECK_MSG(c < 20.0,
                      "T4: iso_def/eta=%.4f exceeds 20 (multi-class merge not "
                      "O(eta)?)", c);
        /* (c) bijective. */
        AIC_CHECK_MSG(bij, "T4: v not bijective");
        AIC_CHECK_MSG(dd(a) > 0.5, "T4: sigma_min=%.4f <= 0.5", dd(a));

        arb_clear(a);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

int main(void)
{
    test_t1_eta0_oracle();
    test_t2_universality();
    test_t3_oblique();
    test_t4_multiclass();
    aic_test_report("test_cstar_build");
    return 0;
}
