/* test_opspace.c — cross-checks for the opspace module O1 (bead aic-zwo):
 * aic_opspace_* (approximate_algebras.tex:1447-1561, §10, th_main_ext). These
 * tests are all on the OPERATOR-NORM ampliated stretch (FINDINGS §C12, the binding
 * landmine): NOT the vacuous Frobenius sigma_min(I_{n^2} (x) M_1). See
 * include/aic_opspace.h and src/aic_opspace_ampliate.c for the contract.
 *
 * THE CHECKS (Rule 6 — cross-checks > unit tests; Rule 5 — every check asserts an
 * invariant, never "runs without error"):
 *   T1 — eta=0 ORACLE (the cleanest faithful operator-norm check). On EXACT-
 *        idempotent channels v is an exact *-iso => a COMPLETE ISOMETRY, so the
 *        ampliated stretch ||v_n||_op == ||v_n^{-1}||_op == 1 (hence a_n == 1) to
 *        machine tol for n = 1,2,4,... up to N_max. The HOPM lower bound EQUALS the
 *        exact value because EVERY unit X maps to a unit (FINDINGS §C12 (a)).
 *   T2 — prop_inc_ext DOUBLING (.tex:1493-1503) on the MEASURED operator-norm a_n:
 *        a_{2n} >= a_n/2 for n = 1,2,4,..., plus the a_1 > 2*delta precondition
 *        (.tex:1505). Exercised through aic_opspace_certify (which guards both).
 *   T3 — operator-norm UNIVERSALITY CANARY (the §C12-faithful .tex:484 check): over
 *        a sweep in dim A AND the ampliation level n, the operator-norm a_n is
 *        BOUNDED and does NOT trend up. Robust bounded + no-upward-trend metric
 *        (learned from cstar_build T2b, FINDINGS §D2/§C11), mutation-proven.
 *   T-struct — STRUCTURAL AMPLIATION tooth (the headline §C12 tooth, the
 *        mathematical content of th_main_ext). Cross-checks the production
 *        ampliation aic_opspace_ampliate_forward (which drives the SAME
 *        aic_opmap_apply_fn the HOPM kernel uses) against an INDEPENDENT explicit
 *        Kronecker block assembly of 1_{M_n}(x)v at n=2 with GENUINELY NONZERO
 *        off-diagonal blocks: the reference's block (k,l) = aic_dhom_v_apply(v,
 *        X_{kl}) (the production v-apply, a path with no opmap coordinate matrix).
 *        They agree to ~1e-12. WHY THIS EXISTS (FINDINGS §C12, 2026-05-31 hostile
 *        review): the reviewer crippled apply_fn to process ONLY the (0,0) block
 *        and ALL prior checks stayed GREEN — the real fixtures give a level-
 *        INDEPENDENT a_n ~ 1.0005, so a structural bug preserving the stretch RATIO
 *        was invisible. MUTATION-D (the (0,0)-only cripple) now makes T-struct RED.
 *   §C12 NON-VACUITY tooth — demonstrate the operator-norm a_n is NOT the trivial
 *        Frobenius sigma_min: a mutation scaling one vE[i] up INFLATES the
 *        operator-norm forward stretch while leaving the Frobenius
 *        sigma_min(I_{n^2}(x)M_1)=sigma_min(M_1) BYTE-FOR-BYTE unchanged (it is the
 *        smallest singular value, blind to scaling one column up, and
 *        level-independent for ANY v). This is exactly the "test that cannot fail"
 *        the Frobenius route would be (CLAUDE.md Rule 5). The n=2 arm is tied to
 *        the STRUCTURAL ampliation (the level-2 off-diagonal block) so MUTATION-D
 *        makes the n=2 arm itself RED — it tests the level-2 ampliation, not just
 *        "scaling a matrix raises its norm" (which n=1 already does).
 *   T-vinv — v^{-1} round-trip (the I1 builder): M_1 . M_1^{-1} = I and
 *        v^{-1}(v(E_i)) = E_i to tol, the explicit standalone inverse check.
 *   T-cross — double vs arb@prec=53 agreement of the stretch (cross-check ladder
 *        rung 2). HONESTY (header aic_opspace.h, finding 4): the HOPM kernel is
 *        PURE DOUBLE — `out` is arb_set_d(double), a zero-radius wrap, NOT an arb
 *        re-evaluation. T-cross therefore tests the prec of the M_1 / M_1^{-1}
 *        ASSEMBLY + dot-product rounding (~1e-15 between prec 53 and 256), a COARSE
 *        gate consistent with aic_dhom_v_sigma_min's posture — NOT a genuine
 *        arb-vs-double ALGORITHM cross-check (which would need an arb HOPM).
 *
 * HONESTY (header): the HOPM is a LOWER bound on the op-norm stretch (a good
 * witness, maybe suboptimal), so O1 certifies the eta=0 complete-isometry, the
 * doubling structure, and the universality canary — NOT ||v||_cb <= 1+O(eps),
 * which needs the SDP UPPER bound (the separate O2 increment).
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_latd.h"
#include "aic/aic_opspace.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"
#include "test_idemp.h"

#define CPREC 256

/* row-major double _Complex buffer element type (matches the opspace public API's
 * `double _Complex *`; the internal aic_opspace_hopm.h `cx` is not included here). */
typedef double _Complex cx;

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Suppress -Wunused for the test_idemp.h builders this TU does not use. */
__attribute__((unused))
static void suppress_unused(void)
{
    (void) make_trace_replace;
    (void) make_identity;
    (void) make_compress_idemp;
    (void) make_dephasing;
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op via the §C5 midpoint route (mirrors
 * test_cstar_build.c). */
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

/* The ampliation levels to sweep: 1, 2, 4, ..., then n_B exactly (the Smith
 * truncation, .tex:1505 / FINDINGS §D3). Fills levels[] (capacity 16), returns the
 * count. Matches aic_opspace_certify's level list so the tests and the pipeline
 * cover the same n. */
static int build_levels(slong levels[16], slong nB)
{
    int k = 0;
    for (slong n = 1; n <= nB && k < 16; n *= 2) levels[k++] = n;
    if (k >= 1 && levels[k - 1] != nB && k < 16) levels[k++] = nB;
    return k;
}

/* ===================== T1: eta=0 oracle ===================================== *
 * v an exact *-iso => a COMPLETE ISOMETRY: ||v_n||_op == ||v_n^{-1}||_op == 1
 * exactly, for n = 1,2,4,... up to the Smith truncation. */
static void t1_one(const char *name, aic_ucp_kraus *phi)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);
    AIC_CHECK_MSG(dd(iso) < 1e-10, "T1 %s: eta=0 iso_def=%.3e not ~0", name, dd(iso));

    arb_t f, g;
    arb_init(f);
    arb_init(g);
    printf("T1 %s: eta=0 COMPLETE ISOMETRY (N=%ld n_B=%ld):\n", name,
           (long) v.n, (long) B.n_B);
    slong lv[16];
    int nlv = build_levels(lv, B.n_B);
    for (int li = 0; li < nlv; li++) {
        slong n = lv[li];
        aic_opspace_forward_stretch(f, &v, n, CPREC);
        aic_opspace_inverse_stretch(g, &v, n, CPREC);
        double a_n = dd(g) > 1e-12 ? 1.0 / dd(g) : 0.0;
        printf("  n=%ld  ||v_n||_op=%.12f  ||v_n^{-1}||_op=%.12f  a_n=%.12f\n",
               (long) n, dd(f), dd(g), a_n);
        AIC_CHECK_MSG(fabs(dd(f) - 1.0) < 1e-9,
                      "T1 %s n=%ld: ||v_n||_op=%.12f != 1 (not a complete isometry)",
                      name, (long) n, dd(f));
        AIC_CHECK_MSG(fabs(dd(g) - 1.0) < 1e-9,
                      "T1 %s n=%ld: ||v_n^{-1}||_op=%.12f != 1", name, (long) n,
                      dd(g));
        AIC_CHECK_MSG(fabs(a_n - 1.0) < 1e-9,
                      "T1 %s n=%ld: a_n=%.12f != 1", name, (long) n, a_n);
    }
    /* certify at eta=0: a_cb == 1, cb_forward == 1, flat == 1. */
    aic_opspace_result r;
    aic_opspace_certify(&r, &v, 0.0, CPREC);
    printf("  certify: a_cb=%.12f cb_forward=%.12f flat=%.9f N_max=%ld\n",
           dd(r.a_cb), dd(r.cb_forward), r.a_cb_flat, (long) r.N_max);
    AIC_CHECK_MSG(fabs(dd(r.a_cb) - 1.0) < 1e-9, "T1 %s: a_cb=%.12f != 1", name,
                  dd(r.a_cb));
    AIC_CHECK_MSG(fabs(dd(r.cb_forward) - 1.0) < 1e-9,
                  "T1 %s: cb_forward=%.12f != 1", name, dd(r.cb_forward));
    AIC_CHECK_MSG(fabs(r.a_cb_flat - 1.0) < 1e-9, "T1 %s: a_cb_flat=%.9f != 1",
                  name, r.a_cb_flat);
    aic_opspace_result_clear(&r);

    arb_clear(g);
    arb_clear(f);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t1_eta0_oracle(void)
{
    printf("\n=== T1 eta=0 oracle: complete isometry (a_n == 1 exact) ===\n");
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 4, 2);     /* M_2 (+) M_2, n_B=4, levels {1,2,4} */
    t1_one("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_block_cond_exp(&phi, 6, 2);     /* M_2 (+) M_4, n_B=6, levels {1,2,4,6} */
    t1_one("block_cond_exp(6,2)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_noiseless_subsystem(&phi, 3, 2); /* M_3, n_B=3, levels {1,2,3} */
    t1_one("noiseless_subsystem(3,2)", &phi);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T2: prop_inc_ext doubling cross-check =============== *
 * a_{2n} >= a_n/2 (.tex:1503) on the MEASURED operator-norm a_n, plus the a_1 >
 * 2*delta precondition (.tex:1505) — both guarded inside aic_opspace_certify (it
 * aborts on violation). T2 confirms certify COMPLETES on genuine eta>0 fixtures
 * (the guards hold) and re-measures the doubling directly for the printed record.
 * Mutation (by hand, recorded): injecting a_2 := a_1/4 trips the doubling abort. */
static void test_t2_doubling(void)
{
    printf("\n=== T2 prop_inc_ext doubling a_{2n} >= a_n/2 (.tex:1503) ===\n");
    struct { slong d, m; double t; } fx[] = {
        {4, 2, 0.03}, {5, 2, 0.03}, {6, 3, 0.02}};
    const int nfx = (int) (sizeof(fx) / sizeof(fx[0]));
    for (int i = 0; i < nfx; i++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, fx[i].d, fx[i].m, fx[i].t, CPREC);
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

        /* direct re-measure of the doubling for the record. */
        arb_t g;
        arb_init(g);
        double a_prev = 0.0;
        printf("  mixconj(%ld,%ld,%.2f): n_B=%ld a_n:", (long) fx[i].d,
               (long) fx[i].m, fx[i].t, (long) B.n_B);
        for (slong n = 1; n <= B.n_B; n *= 2) {
            aic_opspace_inverse_stretch(g, &v, n, CPREC);
            double a_n = dd(g) > 1e-12 ? 1.0 / dd(g) : 0.0;
            printf(" a_%ld=%.6f", (long) n, a_n);
            if (n > 1)
                AIC_CHECK_MSG(a_n >= a_prev / 2.0 - 1e-9,
                              "T2: doubling a_{2n}=%.6f < a_n/2=%.6f at n=%ld",
                              a_n, a_prev / 2.0, (long) n);
            a_prev = a_n;
        }
        printf("\n");
        arb_clear(g);

        /* the certify guards (a_1 > 2*delta AND the doubling) must hold — completes. */
        aic_opspace_result r;
        aic_opspace_certify(&r, &v, dd(iso), CPREC);
        AIC_CHECK_MSG(dd(r.a_cb) > 0.5,
                      "T2: a_cb=%.6f <= 0.5 (lower isometry collapsed)", dd(r.a_cb));
        aic_opspace_result_clear(&r);

        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ===================== T3: operator-norm universality canary ============== *
 * The §C12-faithful .tex:484 / .tex:461 check on the OPERATOR-norm a_n, on BOTH
 * axes (the orchestrator's T3 spec), with the robust bounded + no-upward-trend
 * metric LEARNED from cstar_build T2b (commit f9d2229, FINDINGS §D2/§C11) — NOT a
 * fragile extremes ratio:
 *
 *   AXIS L — AMPLIATION-LEVEL universality (the genuinely NEW th_main_ext content,
 *     §C12 / .tex:484 extended to the level n): for EACH fixture, a_n must NOT
 *     degrade as the ampliation level n grows, i.e. a_cb_flat = max_n a_n/min_n a_n
 *     over n=1,2,4,...,n_B stays ~1. This is the prop_inc_ext / Smith content
 *     measured on the operator norm; mutation-proven (a level-dependent
 *     degradation injection blows a_cb_flat up).
 *
 *   AXIS D — dim A / ambient-N universality (inherited from cstar_build, .tex:461):
 *     the cb-defect c = worst_n(1-a_n)/eta is BOUNDED (abs-max < C_ABS) across a
 *     dim-A sweep. ONLY THE BOUNDEDNESS ARM (mutation-proven: a c=O(dim) law trips
 *     C_ABS). The DISTINCT no-upward-TREND-with-dim-A check is NOT done here: the
 *     opspace fixture array MIXES dim_A=4 (m=2) with dim_A=9 (m=3) AND spans the
 *     d=7 N-geometry spike (FINDINGS §C11), so any position-split "halves-ratio"
 *     here is confounded by ambient-N geometry, not by dim_A — the 2026-05-31
 *     hostile review verified a genuine c∝dim_A law gives halves-ratio 0.206 (does
 *     NOT trip KAPPA=1.6) while the small measured ratio is the n=7 spike landing in
 *     the "lo" half by POSITION. A trend arm here would be a "test that cannot fail"
 *     (CLAUDE.md Rule 5), so it was DROPPED (finding 3a). The RIGOROUS dim-A trend
 *     sweep lives in cstar_build T2b (commit f9d2229, FINDINGS §D2) on the SAME v
 *     (opspace inherits it): bounded abs-max + halves-ratio over a dim-segregated
 *     ambient-n sweep, mutation-proven there (a c·(n/n_min) injection trips it).
 *
 * KEY DIAGNOSTIC (the measured truth): a_cb_flat is ~1.0005 (m=2) / ~1.0001 (m=3)
 * uniformly — the operator-norm a_n is essentially LEVEL-INDEPENDENT; the
 * cross-fixture c variation is pure cstar-geometry (the ambient N), NOT an
 * ampliation effect. So AXIS L is the clean opspace canary; AXIS D re-confirms the
 * inherited dim-independence on the operator norm (boundedness only; the dim-A
 * TREND lives in cstar_build T2b, finding 3a). */
static void test_t3_universality(void)
{
    printf("\n=== T3 operator-norm universality canary (.tex:484/.tex:461, §C12) "
           "===\n");
    const double FLAT_MAX = 1.05; /* AXIS L: a_n does not degrade with level n      */
    const double C_ABS = 2.5;     /* AXIS D: abs-max c=(1-a_n)/eta bounded (data 1.4)*/
    /* Sweep over dim A in {4 (m=2), 9 (m=3)}, spanning the d=7 N-geometry peak and
     * the d=8 crash (FINDINGS §C11). NOTE (finding 3a): NO position-split trend arm
     * here — the array mixes dim_A=4/9 and the d=7 spike confounds it; the rigorous
     * dim-A TREND sweep is cstar_build T2b on the SAME v. Here: boundedness only. */
    struct { slong d, m; double t; } ob[] = {
        {4, 2, 0.03}, {6, 2, 0.03}, {7, 2, 0.03}, {8, 2, 0.03}, /* dim_A=4 */
        {6, 3, 0.02}, {8, 3, 0.02}};                            /* dim_A=9 */
    const int nob = (int) (sizeof(ob) / sizeof(ob[0]));
    double dimv[8];             /* dim A in sweep order (for the abs-max mutation)   */
    double c_abs_max = 0.0;     /* AXIS D: abs-max c=(1-a_n)/eta over the sweep      */
    double flat_max = 0.0;      /* AXIS L: the worst a_cb_flat over fixtures         */
    double flat_worst_an[16];   /* the a_n of the worst-flat fixture (for mutation)  */
    slong flat_worst_lv[16];    /* its levels (for the 1/n mutation, level-aware)    */
    int flat_worst_k = 0;
    double def_abs_max = 0.0;
    arb_t g;
    arb_init(g);
    for (int i = 0; i < nob; i++) {
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
        /* the level-n sweep (1,2,4,...,n_B): collect a_n, track worst cb-defect AND
         * the flatness. */
        slong lv[16];
        int nlv = build_levels(lv, B.n_B);
        double an[16];
        int k = 0;
        double def_worst = 0.0, amax = 0.0, amin = 1e300;
        for (int li = 0; li < nlv && k < 16; li++) {
            slong n = lv[li];
            aic_opspace_inverse_stretch(g, &v, n, CPREC);
            double a_n = dd(g) > 1e-12 ? 1.0 / dd(g) : 0.0;
            an[k++] = a_n;
            if (1.0 - a_n > def_worst) def_worst = 1.0 - a_n;
            if (a_n > amax) amax = a_n;
            if (a_n < amin) amin = a_n;
        }
        double flat = amin > 1e-12 ? amax / amin : 0.0;
        if (flat > flat_max) {
            flat_max = flat;
            flat_worst_k = k;
            for (int j = 0; j < k; j++) {
                flat_worst_an[j] = an[j];
                flat_worst_lv[j] = lv[j];
            }
        }
        if (def_worst > def_abs_max) def_abs_max = def_worst;
        double c = eta > 1e-12 ? def_worst / eta : 0.0;
        dimv[i] = (double) ae.A.dim_A;
        if (c > c_abs_max) c_abs_max = c;
        printf("  mixconj(%ld,%ld): dim_A=%ld N=%ld n_B=%ld eta=%.3e | a_cb_flat="
               "%.6f (AXIS L) worst(1-a_n)=%.3e c=%.4f (AXIS D)\n",
               (long) ob[i].d, (long) ob[i].m, (long) ae.A.dim_A, (long) v.n,
               (long) B.n_B, eta, flat, def_worst, c);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    arb_clear(g);
    printf("  -> AXIS L: worst a_cb_flat=%.6f (<= %.2f); AXIS D: abs-max c=%.4f "
           "(< %.2f) [boundedness only; dim-A trend is cstar_build T2b, finding 3a]"
           "\n", flat_max, FLAT_MAX, c_abs_max, C_ABS);

    /* AXIS L (the opspace-specific canary): a_n does NOT degrade with the
     * ampliation level n — a_cb_flat ~ 1 for every fixture (prop_inc_ext / Smith,
     * §C12). This is on the OPERATOR-norm a_n. */
    AIC_CHECK_MSG(flat_max <= FLAT_MAX,
                  "T3 AXIS L: worst a_cb_flat=%.6f > %.2f — the operator-norm a_n "
                  "DEGRADES with the ampliation level n (th_main_ext / §C12 failure "
                  "mode)", flat_max, FLAT_MAX);
    /* AXIS D (inherited dim-independence, on the operator norm): BOUNDEDNESS ONLY.
     * (i) the lower isometry does not collapse (raw defect small): */
    AIC_CHECK_MSG(def_abs_max < 0.2,
                  "T3 AXIS D: abs-max(1-a_n)=%.4e >= 0.2 — the operator-norm lower "
                  "isometry collapsed", def_abs_max);
    /* (ii) the cb-defect CONSTANT c=(1-a_n)/eta is bounded (dim-independent). The
     * dim-A TREND check is NOT here (finding 3a): the fixture array mixes dim_A=4/9
     * and the d=7 N-spike confounds a position-split ratio (a "test that cannot
     * fail"); the rigorous trend sweep is cstar_build T2b on the SAME v. */
    AIC_CHECK_MSG(c_abs_max < C_ABS,
                  "T3 AXIS D: abs-max c=(1-a_n)/eta=%.4f >= %.2f — the cb-defect "
                  "constant is unbounded (.tex:461, §C12)", c_abs_max, C_ABS);

    /* MUTATION TOOTH (Rule 5/7), AXIS L: a genuine level-dependent degradation
     * a_n := a_n^{(real)} * (n_min/n) (the a_n falling off as 1/n — exactly the
     * "constant depends on the ampliation level" failure th_main_ext forbids)
     * blows a_cb_flat up past FLAT_MAX. Re-run the metric on injected data. */
    {
        double amax = 0.0, amin = 1e300;
        for (int j = 0; j < flat_worst_k; j++) {
            double n = (double) flat_worst_lv[j];  /* the actual level at index j     */
            double a_inj = flat_worst_an[j] / n;   /* 1/n degradation                 */
            if (a_inj > amax) amax = a_inj;
            if (a_inj < amin) amin = a_inj;
        }
        double flat_inj = amin > 1e-12 ? amax / amin : 0.0;
        printf("  [mutation tooth AXIS L] a_n := a_n*(n_min/n): a_cb_flat=%.4f (vs "
               "FLAT_MAX=%.2f)\n", flat_inj, FLAT_MAX);
        AIC_CHECK_MSG(flat_inj > FLAT_MAX,
                      "T3 mutation: a 1/n level-degradation gave a_cb_flat=%.4f, did "
                      "NOT exceed FLAT_MAX=%.2f — AXIS L is BLIND to level growth",
                      flat_inj, FLAT_MAX);
    }
    /* MUTATION TOOTH (Rule 5/7), AXIS D: a genuine c=O(dim) law — the worst-case
     * constant SCALING with dim_A (exactly the .tex:461/.tex:484 failure mode) —
     * would reach c_abs_max * (dim_max/dim_min), tripping the absolute arm. The
     * REAL data (c_abs_max ~ 1.35, fixture geometry) stays under C_ABS; a genuine
     * dim-scaling law clears it. */
    {
        double dmin = dimv[0], dmax = dimv[0];
        for (int i = 0; i < nob; i++) {
            if (dimv[i] < dmin) dmin = dimv[i];
            if (dimv[i] > dmax) dmax = dimv[i];
        }
        double inj_max = c_abs_max * (dmax / dmin);
        printf("  [mutation tooth AXIS D] c_abs_max*(dim_max/dim_min)=%.4f (vs "
               "C_ABS=%.2f; real abs-max c=%.4f stays under)\n", inj_max, C_ABS,
               c_abs_max);
        AIC_CHECK_MSG(inj_max >= C_ABS,
                      "T3 mutation: a c=O(dim) law abs-max=%.4f did NOT exceed "
                      "C_ABS=%.2f — AXIS D is BLIND to dim growth", inj_max, C_ABS);
    }
}

/* ===================== T-struct: structural ampliation tooth ============== *
 * The headline §C12 tooth (the MATHEMATICAL CONTENT of th_main_ext): the block
 * ampliation structure of 1_{M_n}(x)v. Cross-checks aic_opspace_ampliate_forward
 * (the production path, driving the SAME aic_opmap_apply_fn the HOPM kernel uses)
 * against an INDEPENDENT explicit block assembly: the reference block (k,l) is
 * aic_dhom_v_apply(v, X_{kl}) (the production v-apply — a path with NO opmap
 * coordinate matrix, so independent of apply_fn). At n=2 with GENUINELY NONZERO
 * off-diagonal blocks, the two must agree to ~1e-12. The 2026-05-31 hostile review
 * crippled apply_fn to (0,0)-block-only and ALL prior checks stayed GREEN (the real
 * fixtures give level-independent a_n); MUTATION-D now makes THIS RED.
 *
 * COORDINATE NOTE. The opmap forward domain block X_{kl} is an n_B x n_B element of
 * B (B's matrix-unit basis), and aic_opmap_apply_fn computes f(X_{kl}) = projection
 * onto A of v(X_{kl}) = v(X_{kl}) (since v lands in A). aic_dhom_v_apply computes
 * v(X_{kl}) directly from vE. So the two paths agree exactly (up to projection
 * roundoff) — and a (0,0)-only apply_fn zeros the off-diagonal v(X_{kl}) blocks the
 * reference keeps, the mismatch this tooth catches. */

/* Fill the n_B x n_B block-diagonal B-element Xb = sum_i x_i E_i with a determin-
 * istic random combination (LCG seed; no wall-clock RNG, project rule), and record
 * the SAME entries into the (n*n_B)x(n*n_B) row-major buffer X at block (k,l). */
static void fill_B_block(double _Complex *X, acb_mat_t Xb, const aic_dhom_B *B,
                         slong nB, slong n, slong k, slong l, unsigned long *s)
{
    slong cU = n * nB;
    acb_mat_zero(Xb);
    acb_mat_t Ei;
    acb_mat_init(Ei, nB, nB);
    arb_t re, im;
    acb_t coef, t;
    arb_init(re); arb_init(im); acb_init(coef); acb_init(t);
    for (slong i = 0; i < B->dim_B; i++) {
        *s = *s * 6364136223846793005UL + 1442695040888963407UL;
        double xr = ((double) ((*s >> 33) & 0xffff) / 32768.0) - 1.0;
        *s = *s * 6364136223846793005UL + 1442695040888963407UL;
        double xi = ((double) ((*s >> 33) & 0xffff) / 32768.0) - 1.0;
        aic_dhom_B_matunit(Ei, B, i);
        arb_set_d(re, xr); arb_set_d(im, xi);
        acb_set_arb_arb(coef, re, im);
        for (slong a = 0; a < nB; a++)
            for (slong b = 0; b < nB; b++) {
                acb_mul(t, acb_mat_entry(Ei, a, b), coef, CPREC);
                acb_add(acb_mat_entry(Xb, a, b), acb_mat_entry(Xb, a, b), t, CPREC);
            }
    }
    /* copy Xb into block (k,l) of the row-major X buffer */
    for (slong a = 0; a < nB; a++)
        for (slong b = 0; b < nB; b++) {
            cx z = arf_get_d(arb_midref(acb_realref(acb_mat_entry(Xb, a, b))),
                             ARF_RND_NEAR)
                 + arf_get_d(arb_midref(acb_imagref(acb_mat_entry(Xb, a, b))),
                             ARF_RND_NEAR) * I;
            X[(k * nB + a) * cU + (l * nB + b)] = z;
        }
    acb_clear(t); acb_clear(coef); arb_clear(im); arb_clear(re);
    acb_mat_clear(Ei);
}

/* max |Y_test - Y_ref| over the (n*N)x(n*N) ampliation, where Y_test is built by
 * aic_opspace_ampliate_forward and Y_ref is the INDEPENDENT block assembly (each
 * block = aic_dhom_v_apply(v, X_{kl})). Returns the max abs entry difference and,
 * via *off, the max abs entry of Y_ref in the OFF-diagonal blocks (k!=l) so the
 * caller can confirm the off-diagonal content is genuinely nonzero. */
static double struct_ampliation_maxdiff(const aic_dhom_v *v, slong n,
                                        unsigned long seed, double *off)
{
    slong nB = v->B->n_B, N = v->n, cU = n * nB, cV = n * N;
    cx *X = flint_malloc((size_t) cU * cU * sizeof(cx));
    cx *Yt = flint_malloc((size_t) cV * cV * sizeof(cx));
    for (slong p = 0; p < cU * cU; p++) X[p] = 0.0;
    acb_mat_t Xb;
    acb_mat_init(Xb, nB, nB);
    /* build X with NONZERO off-diagonal blocks */
    unsigned long s = seed;
    acb_mat_t *Xkl = flint_malloc((size_t) (n * n) * sizeof(acb_mat_t));
    for (slong k = 0; k < n; k++)
        for (slong l = 0; l < n; l++) {
            acb_mat_init(Xkl[k * n + l], nB, nB);
            fill_B_block(X, Xb, v->B, nB, n, k, l, &s);
            acb_mat_set(Xkl[k * n + l], Xb);
        }
    /* test path (production ampliation) */
    aic_opspace_ampliate_forward(Yt, v, X, n, CPREC);
    /* independent reference: block (k,l) = aic_dhom_v_apply(v, X_{kl}) */
    acb_mat_t vXkl;
    acb_mat_init(vXkl, N, N);
    double maxdiff = 0.0, offmax = 0.0;
    for (slong k = 0; k < n; k++)
        for (slong l = 0; l < n; l++) {
            aic_dhom_v_apply(vXkl, v, Xkl[k * n + l], CPREC);
            for (slong a = 0; a < N; a++)
                for (slong b = 0; b < N; b++) {
                    cx ref = arf_get_d(arb_midref(acb_realref(
                                 acb_mat_entry(vXkl, a, b))), ARF_RND_NEAR)
                           + arf_get_d(arb_midref(acb_imagref(
                                 acb_mat_entry(vXkl, a, b))), ARF_RND_NEAR) * I;
                    cx test = Yt[(k * N + a) * cV + (l * N + b)];
                    double d = cabs(test - ref);
                    if (d > maxdiff) maxdiff = d;
                    if (k != l && cabs(ref) > offmax) offmax = cabs(ref);
                }
        }
    acb_mat_clear(vXkl);
    for (slong p = 0; p < n * n; p++) acb_mat_clear(Xkl[p]);
    flint_free(Xkl);
    acb_mat_clear(Xb);
    flint_free(Yt);
    flint_free(X);
    if (off) *off = offmax;
    return maxdiff;
}

static void test_struct_ampliation(void)
{
    printf("\n=== T-struct structural ampliation 1_{M_n}(x)v vs Kronecker ref "
           "(§C12 headline) ===\n");
    struct { slong d, m; double t; } fx[] = {{4, 2, 0.0}, {6, 3, 0.02}};
    const int nfx = (int) (sizeof(fx) / sizeof(fx[0]));
    for (int i = 0; i < nfx; i++) {
        aic_ucp_kraus phi;
        aic_assoc_ecstar ae;
        arb_t iso, ea;
        aic_dhom_B B;
        aic_dhom_v v;
        double eps;
        if (fx[i].t == 0.0) {
            make_block_cond_exp(&phi, fx[i].d, fx[i].m); /* eta=0 oracle fixture */
            eps = 0.0;
        } else {
            make_mixconj(&phi, fx[i].d, fx[i].m, fx[i].t, CPREC);
        }
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);
        if (fx[i].t != 0.0) {
            arb_init(ea);
            aic_ecstar_defect_assoc(ea, &ae.A, CPREC);
            eps = dd(ea);
            arb_clear(ea);
        }
        arb_init(iso);
        aic_cstar_build(&B, &v, iso, &ae.A, eps, CPREC);
        /* n=2 AND n=3: genuinely multi-block ampliations with nonzero off-diag. */
        for (slong n = 2; n <= 3; n++) {
            double off = 0.0;
            double md = struct_ampliation_maxdiff(&v, n, 0x9e3779b9UL + 7 * n, &off);
            printf("  fixture(%ld,%ld,t=%.2f) n=%ld: max|ampliate - block_ref|=%.3e "
                   "(off-diag ref content=%.4f)\n", (long) fx[i].d, (long) fx[i].m,
                   fx[i].t, (long) n, md, off);
            AIC_CHECK_MSG(off > 0.1,
                          "T-struct: off-diagonal reference content=%.4f ~ 0 — the "
                          "test X has no off-diagonal blocks (tooth is toothless)",
                          off);
            AIC_CHECK_MSG(md < 1e-9,
                          "T-struct fixture(%ld,%ld) n=%ld: ampliation 1_{M_n}(x)v "
                          "DISAGREES with the independent block assembly by %.3e "
                          "(th_main_ext block structure broken; e.g. a (0,0)-only "
                          "apply_fn — MUTATION-D)", (long) fx[i].d, (long) fx[i].m,
                          (long) n, md);
        }
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ===================== §C12 non-vacuity tooth ============================== *
 * The operator-norm a_n is NOT the trivial Frobenius sigma_min. Scaling one vE[i]
 * UP inflates the OPERATOR-norm forward stretch ||v_n||_op (the HOPM catches it),
 * while the Frobenius sigma_min(I_{n^2}(x)M_1)=sigma_min(M_1) — the smallest
 * singular value, blind to scaling one column up, and LEVEL-INDEPENDENT for ANY v
 * — is UNCHANGED. A Frobenius-sigma_min ampliation "canary" for th_main_ext would
 * be BLIND to this operator-norm-only effect (CLAUDE.md Rule 5 "test that cannot
 * fail"); the operator-norm route is faithful. */
static void test_c12_nonvacuity(void)
{
    printf("\n=== §C12 non-vacuity: operator-norm a_n != Frobenius sigma_min ===\n");
    aic_ucp_kraus phi;
    make_mixconj(&phi, 6, 3, 0.02, CPREC);
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

    arb_t f, sm;
    arb_init(f);
    arb_init(sm);

    /* PRE: forward stretch at n=1 and n=2, and the Frobenius sigma_min. */
    aic_opspace_forward_stretch(f, &v, 1, CPREC);
    double fwd1_pre = dd(f);
    aic_opspace_forward_stretch(f, &v, 2, CPREC);
    double fwd2_pre = dd(f);
    aic_dhom_v_sigma_min(sm, &v, CPREC);
    double sm_pre = dd(sm);

    /* MUTATION: scale vE[0] by 1.6 (inflates the operator-norm forward stretch). */
    arb_t sc;
    arb_init(sc);
    arb_set_d(sc, 1.6);
    acb_mat_scalar_mul_arb(v.vE[0], v.vE[0], sc, CPREC);
    arb_clear(sc);

    aic_opspace_forward_stretch(f, &v, 1, CPREC);
    double fwd1_post = dd(f);
    aic_opspace_forward_stretch(f, &v, 2, CPREC);
    double fwd2_post = dd(f);
    aic_dhom_v_sigma_min(sm, &v, CPREC);
    double sm_post = dd(sm);

    /* LEVEL-2-SPECIFIC arm (finding 2): the inflated vE[0] in an OFF-DIAGONAL block
     * of a level-2 X. The production n=2 ampliation must STILL match the independent
     * block reference (both see the mutated vE[0]); a (0,0)-only apply_fn (MUTATION-D)
     * zeros that off-diagonal block while the reference keeps the inflated ~1.6*vE[0],
     * so this struct-diff blows up — making the §C12 n=2 arm itself RED at level 2
     * (NOT just "scaling a matrix raises its norm", which the n=1 stretch already is). */
    double off2 = 0.0;
    double struct2 = struct_ampliation_maxdiff(&v, 2, 0xC12C12C1UL, &off2);

    printf("  PRE : ||v_1||_op=%.6f ||v_2||_op=%.6f | Frobenius sigma_min=%.9f\n",
           fwd1_pre, fwd2_pre, sm_pre);
    printf("  POST: ||v_1||_op=%.6f ||v_2||_op=%.6f | Frobenius sigma_min=%.9f\n",
           fwd1_post, fwd2_post, sm_post);
    printf("  operator-norm stretch MOVED: d(n=1)=%.4f d(n=2)=%.4f | Frobenius "
           "sigma_min delta=%.2e | level-2 struct-diff=%.3e (off-diag ref=%.4f)\n",
           fwd1_post - fwd1_pre, fwd2_post - fwd2_pre, fabs(sm_post - sm_pre),
           struct2, off2);

    /* (a) the operator-norm forward stretch RESPONDED to the mutation (it should:
     * the inflated vE[0] is a witness direction with op-norm > 1). The genuine move
     * is ~+0.60 (1.6x scaling); threshold tightened from the loose 0.3 toward it. */
    AIC_CHECK_MSG(fwd1_post - fwd1_pre > 0.5,
                  "§C12: scaling vE[0] by 1.6 moved the operator-norm forward "
                  "stretch by only %.4f (n=1: %.6f -> %.6f); expected ~+0.60 — the "
                  "op-norm path is wrong", fwd1_post - fwd1_pre, fwd1_pre, fwd1_post);
    AIC_CHECK_MSG(fwd2_post - fwd2_pre > 0.5,
                  "§C12: operator-norm forward stretch at n=2 moved by only %.4f "
                  "(%.6f -> %.6f); expected ~+0.60", fwd2_post - fwd2_pre, fwd2_pre,
                  fwd2_post);
    /* (a') LEVEL-2 STRUCTURE: the n=2 ampliation of the MUTATED v matches the
     * independent block reference (off-diagonal blocks carry the inflated vE[0]).
     * MUTATION-D (a (0,0)-only ampliation) makes this RED — the n=2 arm is now
     * level-2-specific, not a rescaled copy of the n=1 arm. */
    AIC_CHECK_MSG(off2 > 0.1,
                  "§C12 n=2: off-diagonal reference content=%.4f ~ 0 (the inflated "
                  "vE[0] is not exercised off-diagonal — arm is toothless)", off2);
    AIC_CHECK_MSG(struct2 < 1e-9,
                  "§C12 n=2: the level-2 ampliation of the mutated v DISAGREES with "
                  "the independent block assembly by %.3e — a (0,0)-only ampliation "
                  "(MUTATION-D) drops the off-diagonal inflated vE[0] block", struct2);
    /* (b) the Frobenius sigma_min was BLIND to the same mutation (it is the
     * smallest singular value of M_1, unchanged by scaling one column UP; and
     * sigma_min(I_{n^2}(x)M_1)=sigma_min(M_1) is level-independent for ANY v). A
     * Frobenius-sigma_min ampliation canary CANNOT see this operator-norm effect. */
    AIC_CHECK_MSG(fabs(sm_post - sm_pre) < 1e-6,
                  "§C12: the Frobenius sigma_min CHANGED (%.9f -> %.9f) under a "
                  "mutation meant to be sigma_min-blind — the tooth is not sharp",
                  sm_pre, sm_post);

    arb_clear(sm);
    arb_clear(f);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T-vinv: v^{-1} round-trip (I1 public builder) ======= *
 * The factorize/D4 interface add: aic_opspace_build_vinv exposes v^{-1} as dim_A
 * operators vinvB[k] = v^{-1}(B_k) in B. Standalone round-trip (finding 5):
 *   v^{-1}(v(E_i)) = E_i for every matrix unit E_i (i.e. M_1 . M_1^{-1} = I read
 * through the bases). v(E_i) = vE[i] in A; its A-coords are c_k = <B_k, vE[i]>_F;
 * v^{-1}(v(E_i)) = sum_k c_k vinvB[k], which must equal E_i. On eta=0 this is exact
 * to machine prec; on eta>0 the M_1 inverse is still exact (it is a genuine matrix
 * inverse, not an O(eps) approximation), so the round-trip is exact regardless. */
static void t_vinv_one(const char *name, aic_ucp_kraus *phi, double tval)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    double eps = 0.0;
    if (tval != 0.0) {
        arb_t ea;
        arb_init(ea);
        aic_ecstar_defect_assoc(ea, &ae.A, CPREC);
        eps = dd(ea);
        arb_clear(ea);
    }
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, eps, CPREC);

    acb_mat_t *vinvB = NULL;
    slong dimA = 0, nB = 0;
    aic_opspace_build_vinv(&vinvB, &dimA, &nB, &v, CPREC);
    AIC_CHECK_MSG(dimA == v.A->dim_A, "T-vinv %s: build_vinv dim_A=%ld != %ld",
                  name, (long) dimA, (long) v.A->dim_A);
    AIC_CHECK_MSG(nB == B.n_B, "T-vinv %s: build_vinv n_B=%ld != %ld", name,
                  (long) nB, (long) B.n_B);

    slong N = v.n;
    acb_mat_t acc, Ei, term;
    acb_t ck;
    acb_mat_init(acc, nB, nB);
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(term, nB, nB);
    acb_init(ck);
    arb_t err, worst;
    arb_init(err);
    arb_init(worst);
    arb_zero(worst);
    for (slong i = 0; i < B.dim_B; i++) {
        /* v^{-1}(v(E_i)) = sum_k <B_k, vE[i]>_F vinvB[k]. */
        acb_mat_zero(acc);
        for (slong k = 0; k < dimA; k++) {
            acb_zero(ck);
            acb_t z;
            acb_init(z);
            for (slong p = 0; p < N; p++)
                for (slong q = 0; q < N; q++) {
                    acb_conj(z, acb_mat_entry(ae.A.B[k], p, q));
                    acb_mul(z, z, acb_mat_entry(v.vE[i], p, q), CPREC);
                    acb_add(ck, ck, z, CPREC);
                }
            acb_clear(z);
            acb_mat_scalar_mul_acb(term, vinvB[k], ck, CPREC);
            acb_mat_add(acc, acc, term, CPREC);
        }
        aic_dhom_B_matunit(Ei, &B, i);
        acb_mat_sub(acc, acc, Ei, CPREC);
        for (slong p = 0; p < nB; p++)
            for (slong q = 0; q < nB; q++) {
                acb_abs(err, acb_mat_entry(acc, p, q), CPREC);
                if (arb_gt(err, worst)) arb_set(worst, err);
            }
    }
    printf("  T-vinv %s: max|v^{-1}(v(E_i)) - E_i| over %ld units = %.3e\n", name,
           (long) B.dim_B, dd(worst));
    AIC_CHECK_MSG(dd(worst) < 1e-9,
                  "T-vinv %s: v^{-1}(v(E_i)) != E_i (max err %.3e) — the inverse "
                  "coordinate map M_1^{-1} is wrong (M_1.M_1^{-1} != I)", name,
                  dd(worst));

    arb_clear(worst);
    arb_clear(err);
    acb_clear(ck);
    acb_mat_clear(term);
    acb_mat_clear(Ei);
    acb_mat_clear(acc);
    aic_opspace_vinv_clear(vinvB, dimA);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_vinv_roundtrip(void)
{
    printf("\n=== T-vinv v^{-1} round-trip: v^{-1}(v(E_i)) == E_i (I1 builder) ===\n");
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 4, 2);        /* eta=0 oracle */
    t_vinv_one("block_cond_exp(4,2) eta=0", &phi, 0.0);
    aic_ucp_kraus_clear(&phi);
    make_mixconj(&phi, 6, 3, 0.02, CPREC);  /* eta>0 */
    t_vinv_one("mixconj(6,3,0.02)", &phi, 0.02);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T-cross: double vs arb@prec=53 ===================== *
 * Cross-check ladder rung 2: the operator-norm stretch at prec=53 vs prec=256
 * agree to ~1e-12. HONESTY (header aic_opspace.h, finding 4): the HOPM kernel is
 * PURE DOUBLE — `out` is arb_set_d(double), a zero-radius wrap, NOT an arb
 * re-evaluation. So this tests the prec of the M_1 / M_1^{-1} ASSEMBLY + the
 * dot-product rounding stability (~1e-15), a COARSE gate consistent with
 * aic_dhom_v_sigma_min's posture — NOT a genuine arb-vs-double ALGORITHM
 * cross-check (which would need an arb HOPM). */
static void test_cross_prec(void)
{
    printf("\n=== T-cross double vs arb@prec=53 (cross-check ladder rung 2) ===\n");
    aic_ucp_kraus phi;
    make_mixconj(&phi, 5, 2, 0.03, CPREC);
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

    arb_t f53, f256, g53, g256, tol;
    arb_init(f53);
    arb_init(f256);
    arb_init(g53);
    arb_init(g256);
    arb_init(tol);
    arb_set_d(tol, 1e-10);
    for (slong n = 1; n <= B.n_B; n++) {
        aic_opspace_forward_stretch(f53, &v, n, 53);
        aic_opspace_forward_stretch(f256, &v, n, 256);
        aic_opspace_inverse_stretch(g53, &v, n, 53);
        aic_opspace_inverse_stretch(g256, &v, n, 256);
        printf("  n=%ld fwd@53=%.12f fwd@256=%.12f | inv@53=%.12f inv@256=%.12f\n",
               (long) n, dd(f53), dd(f256), dd(g53), dd(g256));
        AIC_CHECK_MSG(fabs(dd(f53) - dd(f256)) < 1e-10,
                      "T-cross n=%ld: fwd 53 vs 256 differ by %.2e", (long) n,
                      fabs(dd(f53) - dd(f256)));
        AIC_CHECK_MSG(fabs(dd(g53) - dd(g256)) < 1e-10,
                      "T-cross n=%ld: inv 53 vs 256 differ by %.2e", (long) n,
                      fabs(dd(g53) - dd(g256)));
    }
    arb_clear(tol);
    arb_clear(g256);
    arb_clear(g53);
    arb_clear(f256);
    arb_clear(f53);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_t1_eta0_oracle();
    test_t2_doubling();
    test_t3_universality();
    test_struct_ampliation();
    test_c12_nonvacuity();
    test_vinv_roundtrip();
    test_cross_prec();
    aic_test_report("test_opspace");
    return 0;
}
