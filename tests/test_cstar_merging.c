/* test_cstar_merging.c — cross-checks for cstar_build Increment 3 (bead aic-097):
 * lem_merging (.tex:1325-1346), the GENERAL 2x2 block-assembly lemma that
 * lem_extension (I4) consumes. cor_merge_sum (I2) is its off-diagonal-zero special
 * case (.tex:1352 = lem_merging with gamma_12 = gamma_21 = 0), so C2 cross-checks
 * the merged v + its certified balls against aic_cstar_merge_sum to machine
 * precision (inheriting I2's oblique STAR teeth and FINDINGS §C8 magnitude check).
 *
 * NO NEW double-vs-arb path: the new code (aic_cstar_merging.c) is STRUCTURAL GLUE
 * (routing + concat + certification) over already-tested numerical routines
 * (aic_dhom_defect_sweep / aic_dhom_v_sigma_min carry the arb path; the merge_cond
 * op-norms use the certified mid+radius UB aic_corner_gamma_opnorm_ub). The
 * cross-checks are (C1) the eta=0 full-assembly identity oracle (exact iso, routing
 * exercised on off-diagonal units), (C2) the cor_merge_sum consistency on the
 * OBLIQUE mixconj(5,3) fixture (the star teeth, FINDINGS §C2/§C8), (C3) the
 * merging-condition input guards (merge_cond_max + sigma_min teeth).
 *
 * MUTATION-PROVEN teeth (Rule 7; each RED observed BY HAND, then restored — NOT
 * left as failing cases in the automated suite). The by-hand runs are recorded in
 * the per-test comments and the increment report.
 */
#include <math.h>
#include <stdio.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_corner.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
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
    (void) make_block_cond_exp;
    (void) make_trace_replace;
    (void) make_noiseless_subsystem;
    (void) make_identity;
    (void) make_dephasing;
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op (mirrors test_cstar_merge). */
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

/* n x n diagonal projector with diagonal entries [lo,hi) set to 1 (rest 0). */
static void diag_proj_range(acb_mat_t P, slong n, slong lo, slong hi)
{
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    for (slong i = lo; i < hi; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
}

/* ||A - B||_op via the certified mid+radius UB-equivalent midpoint route (the same
 * Gram false-fail workaround as production; FINDINGS §C5). */
static double opnorm_diff_d(const acb_mat_t Aop, const acb_mat_t Bop, slong prec)
{
    slong r = acb_mat_nrows(Aop), c = acb_mat_ncols(Aop);
    acb_mat_t D, M;
    acb_mat_init(D, r, c);
    acb_mat_init(M, r, c);
    acb_mat_sub(D, Aop, Bop, prec);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    double o = dd(e);
    arb_clear(e);
    acb_mat_clear(M); acb_mat_clear(D);
    return o;
}

/* ===================== C1: eta=0 full-assembly oracle ===================== *
 * A = genuine M_3 (aic_cstar_matrix_algebra: identity channel, star = plain, basis
 * = matrix units E_{lm} at index l*d+m). B = M_3 with n1=2, n2=1; Pi_1 = E00+E11,
 * Pi_2 = E22; P_1 = diag(1,1,0), P_2 = diag(0,0,1) in A. The NATURAL gamma (each
 * gamma_{jk} the identity inclusion: block-local E_{ab} -> the GLOBAL matrix unit
 * E_{lm} it represents, an operator in A=M_3). Then merged gamma = identity:
 *   mult_def ~ machine-zero (exact iso), sigma_min ~ 1, merge_cond_max ~ 0, and the
 *   merged v IS the identity (v(E_i) == E_i). The routing is exercised on the
 *   off-diagonal units E_02 (block (1,2)) and E_20 (block (2,1)). */
static void test_c1_eta0_oracle(void)
{
    printf("C1 lem_merging eta=0 full-assembly oracle (A=M_3, n1=2 n2=1, natural "
           "identity gamma):\n");
    const slong n = 3, n1 = 2, n2 = 1;
    aic_ecstar A;
    aic_cstar_matrix_algebra(&A, n, CPREC);

    /* P_1 = diag(1,1,0), P_2 = diag(0,0,1). */
    acb_mat_t P1, P2;
    diag_proj_range(P1, n, 0, 2);
    diag_proj_range(P2, n, 2, 3);

    /* Build the four natural-identity gamma_{jk} as contiguous acb_mat_t arrays.
     * gamma_{jk}(E^{block}_{ab}) = E_{lm} (global) with l = a + off_j, m = b + off_k,
     * off_1 = 0, off_2 = n1. Each image is the single-1 n x n matrix at (l,m). */
    slong off[2] = {0, n1};
    slong sz[2]  = {n1, n2};
    acb_mat_t *g[2][2];
    for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++) {
            slong rows = sz[j], cols = sz[k];
            g[j][k] = flint_malloc((size_t) (rows * cols) * sizeof(acb_mat_t));
            for (slong a = 0; a < rows; a++)
                for (slong b = 0; b < cols; b++) {
                    acb_mat_t *E = &g[j][k][a * cols + b];
                    acb_mat_init(*E, n, n);
                    acb_mat_zero(*E);
                    acb_set_si(acb_mat_entry(*E, a + off[j], b + off[k]), 1);
                }
        }

    aic_dhom_B Bm;
    aic_dhom_v vm;
    arb_t mult_def, smin, mcm;
    arb_init(mult_def); arb_init(smin); arb_init(mcm);
    aic_cstar_lem_merging(&Bm, &vm, mult_def, smin, mcm, /*two_block=*/0, n1, n2,
                          (const acb_mat_t *) g[0][0], (const acb_mat_t *) g[0][1],
                          (const acb_mat_t *) g[1][0], (const acb_mat_t *) g[1][1],
                          &A, P1, P2, CPREC);

    double md = dd(mult_def), sm = dd(smin), mc = dd(mcm);
    printf("  merged dim_B=%ld (expect 9), mult_def=%.3e (~0), sigma_min=%.4f (~1), "
           "merge_cond_max=%.3e (~0)\n", (long) Bm.dim_B, md, sm, mc);
    AIC_CHECK_MSG(Bm.dim_B == 9, "C1: merged dim_B=%ld != 9", (long) Bm.dim_B);
    AIC_CHECK_MSG(Bm.num_blocks == 1 && Bm.d[0] == 3,
                  "C1: merged B not a single M_3 block");
    AIC_CHECK_MSG(md < 1e-12, "C1: eta=0 mult_def not ~0 (%.3e)", md);
    AIC_CHECK_MSG(fabs(sm - 1.0) < 1e-9, "C1: eta=0 sigma_min=%.6f != 1", sm);
    AIC_CHECK_MSG(mc < 1e-12, "C1: eta=0 merge_cond_max not ~0 (%.3e)", mc);

    /* The merged v IS the identity: v(E_{lm}) == E_{lm} (the single-1 at (l,m)). */
    double worst = 0.0;
    for (slong l = 0; l < n; l++)
        for (slong m = 0; m < n; m++) {
            acb_mat_t Elm;
            acb_mat_init(Elm, n, n);
            acb_mat_zero(Elm);
            acb_set_si(acb_mat_entry(Elm, l, m), 1);
            double e = opnorm_diff_d(vm.vE[l * n + m], Elm, CPREC);
            if (e > worst) worst = e;
            acb_mat_clear(Elm);
        }
    printf("  max||v(E_lm) - E_lm||_op = %.3e (merged v == identity)\n", worst);
    AIC_CHECK_MSG(worst < 1e-12, "C1: merged v != identity (%.3e)", worst);

    /* Routing spot-check: E_02 (l=0<n1, m=2>=n1) is in block (1,2), local (0,0), so
     * v(E_02) == gamma_12[0] == E_02. E_20 (l=2>=n1, m=0<n1) is in block (2,1),
     * local (0,0), v(E_20) == gamma_21[0] == E_20. (Already covered by the identity
     * sweep above, but assert the off-diagonal routes explicitly.) */
    double r02 = opnorm_diff_d(vm.vE[0 * n + 2], g[0][1][0], CPREC); /* block (1,2) */
    double r20 = opnorm_diff_d(vm.vE[2 * n + 0], g[1][0][0], CPREC); /* block (2,1) */
    printf("  routing: ||v(E_02) - gamma_12[0]||=%.3e, ||v(E_20) - gamma_21[0]||="
           "%.3e\n", r02, r20);
    AIC_CHECK_MSG(r02 < 1e-14, "C1: E_02 mis-routed (block (1,2), %.3e)", r02);
    AIC_CHECK_MSG(r20 < 1e-14, "C1: E_20 mis-routed (block (2,1), %.3e)", r20);

    arb_clear(mcm); arb_clear(smin); arb_clear(mult_def);
    aic_dhom_v_clear(&vm); aic_dhom_B_clear(&Bm);
    for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++) {
            slong rows = sz[j], cols = sz[k];
            for (slong i = 0; i < rows * cols; i++) acb_mat_clear(g[j][k][i]);
            flint_free(g[j][k]);
        }
    acb_mat_clear(P2); acb_mat_clear(P1);
    aic_ecstar_clear(&A);
}

/* Build a rank-1 commutative inclusion v_j : M_1 -> S_{P_j}, v_j(1) = Ptilde_j. */
static void build_rank1_incl(aic_dhom_B *B, aic_dhom_v *v, const aic_ecstar *A,
                             const acb_mat_t Pj, slong prec)
{
    slong dims[1] = {1};
    aic_dhom_B_init(B, dims, 1);
    aic_dhom_v_init(v, B, A);
    aic_corner_Ptilde(v->vE[0], A, Pj, prec);
}

/* ===================== C2: cor_merge_sum consistency (oblique) ============= *
 * lem_merging with two_block=1 (B = M_1 (+) M_1, off-diagonal blocks EMPTY) IS
 * cor_merge_sum (.tex:1352; FINDINGS §C9 — the GENUINE reduction is the two-block
 * shape, NOT a single M_2 block with zeroed off-diagonal, which would be an invalid
 * input whose live (E_01,E_10) pair gives an O(1) defect). On make_mixconj(5,3)
 * (oblique eta>0; FINDINGS §C4), with in-A 1d projections P_1=span(e1), P_2=span(e2),
 * build the SAME two commutative inclusions v_1, v_2 (M_1 blocks, v_j(1)=Ptilde_j)
 * as gamma_11, gamma_22 (gamma_12, gamma_21 unused). Assert the merged v + its
 * mult_def + sigma_min MATCH aic_cstar_merge_sum's output (from I2) to machine
 * precision. This inherits I2's oblique STAR teeth and the FINDINGS §C8 c=defect/eta
 * MAGNITUDE check (c < 0.2; mutation: star -> plain pushes c to ~0.43 -> RED, run by
 * hand). */
static void test_c2_merge_sum_consistency(void)
{
    printf("C2 lem_merging(two_block) == cor_merge_sum (mixconj(5,3), "
           "P_1=span(e1), P_2=span(e2)):\n");
    const slong n = 5;
    const double ts[] = {0.06, 0.02};
    double cs[2];
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, 5, 3, ts[it], CPREC);
        double eta = eta_proxy(&phi, CPREC);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

        acb_mat_t P1, P2;
        diag_proj_range(P1, n, 1, 2);          /* span(e1) */
        diag_proj_range(P2, n, 2, 3);          /* span(e2) */

        /* Reference: cor_merge_sum on the two rank-1 inclusions (I2). */
        aic_dhom_B B1, B2;
        aic_dhom_v v1, v2;
        build_rank1_incl(&B1, &v1, &ae.A, P1, CPREC);
        build_rank1_incl(&B2, &v2, &ae.A, P2, CPREC);
        aic_dhom_B Bref;
        aic_dhom_v vref;
        arb_t md_ref, sm_ref;
        arb_init(md_ref); arb_init(sm_ref);
        aic_cstar_merge_sum(&Bref, &vref, md_ref, sm_ref, &v1, &v2, &ae.A, CPREC);

        /* lem_merging two_block: n1=n2=1, gamma_11 = {Ptilde_1}, gamma_22 =
         * {Ptilde_2}; gamma_12/gamma_21 unused (off-diagonal blocks empty). */
        acb_mat_t g11[1], g22[1];
        acb_mat_init(g11[0], n, n);
        acb_mat_init(g22[0], n, n);
        aic_corner_Ptilde(g11[0], &ae.A, P1, CPREC);
        aic_corner_Ptilde(g22[0], &ae.A, P2, CPREC);

        aic_dhom_B Bm;
        aic_dhom_v vm;
        arb_t md_m, sm_m, mcm;
        arb_init(md_m); arb_init(sm_m); arb_init(mcm);
        aic_cstar_lem_merging(&Bm, &vm, md_m, sm_m, mcm, /*two_block=*/1, 1, 1,
                              (const acb_mat_t *) g11, NULL, NULL,
                              (const acb_mat_t *) g22, &ae.A, P1, P2, CPREC);

        /* The merged v must MATCH cor_merge_sum's vref EXACTLY: both are over
         * B = M_1 (+) M_1 (dim_B = 2), vE = [Ptilde_1, Ptilde_2]. */
        double e0 = opnorm_diff_d(vm.vE[0], vref.vE[0], CPREC); /* Ptilde_1 */
        double e1 = opnorm_diff_d(vm.vE[1], vref.vE[1], CPREC); /* Ptilde_2 */

        double md_m_d = dd(md_m), md_ref_d = dd(md_ref);
        double sm_m_d = dd(sm_m), sm_ref_d = dd(sm_ref);
        double mc = dd(mcm);
        double c = eta > 1e-14 ? md_m_d / eta : 0.0;
        cs[it] = c;
        printf("  t=%.2f eta=%.3e | vE match e0=%.2e e1=%.2e | mult_def m=%.3e "
               "ref=%.3e | sigma_min m=%.4f ref=%.4f | mcm=%.3e c=md/eta=%.3f\n",
               ts[it], eta, e0, e1, md_m_d, md_ref_d, sm_m_d, sm_ref_d, mc, c);

        AIC_CHECK_MSG(Bm.dim_B == 2, "C2: merged dim_B=%ld != 2", (long) Bm.dim_B);
        AIC_CHECK_MSG(Bm.num_blocks == 2, "C2: merged B not two blocks (got %ld)",
                      (long) Bm.num_blocks);
        AIC_CHECK_MSG(e0 < 1e-12 && e1 < 1e-12,
                      "C2: merged vE != cor_merge_sum (e0=%.3e e1=%.3e)", e0, e1);
        AIC_CHECK_MSG(fabs(md_m_d - md_ref_d) < 1e-9,
                      "C2: mult_def disagrees with cor_merge_sum (%.3e vs %.3e)",
                      md_m_d, md_ref_d);
        AIC_CHECK_MSG(fabs(sm_m_d - sm_ref_d) < 1e-9,
                      "C2: sigma_min disagrees with cor_merge_sum (%.4f vs %.4f)",
                      sm_m_d, sm_ref_d);
        /* FINDINGS §C8: c = mult_def/eta magnitude is the sharp star tooth (c < 0.2;
         * a star -> plain mutation pushes c to ~0.43 -> RED, run by hand). */
        AIC_CHECK_MSG(c < 0.2, "C2: mult_def/eta = %.3f exceeds 0.2 — star teeth RED "
                      "(plain-product gap ~0.43) (t=%.2f)", c, ts[it]);
        /* merge_cond_max: merging2 ||Ptilde_j - P_j|| = in-A residual O(eta),
         * merging0/3 hold (Ptilde_j Hermitian, ||Ptilde_j||~1). So mcm = O(eta). */
        AIC_CHECK_MSG(mc < 10.0 * eta + 1e-9,
                      "C2: merge_cond_max=%.3e not O(eta=%.3e)", mc, eta);

        arb_clear(mcm); arb_clear(sm_m); arb_clear(md_m);
        acb_mat_clear(g22[0]); acb_mat_clear(g11[0]);
        aic_dhom_v_clear(&vm); aic_dhom_B_clear(&Bm);
        arb_clear(sm_ref); arb_clear(md_ref);
        aic_dhom_v_clear(&vref); aic_dhom_B_clear(&Bref);
        aic_dhom_v_clear(&v2); aic_dhom_B_clear(&B2);
        aic_dhom_v_clear(&v1); aic_dhom_B_clear(&B1);
        acb_mat_clear(P2); acb_mat_clear(P1);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("  -> c=md/eta: t=0.06 %.3f, t=0.02 %.3f (both << 0.2 star threshold)\n",
           cs[0], cs[1]);
}

/* ===================== C3: merging-condition teeth ======================== *
 * The C1 setup (A=M_3, natural identity gamma), but ONE gamma_{jk} perturbed to
 * VIOLATE a merging condition. (a) scale gamma_12's image by 0.5 -> merging3 lower
 * bound violated AND the off-diagonal direction half-collapses -> merge_cond_max
 * grows AND sigma_min drops below 1-eps. (b) set gamma_21(E_20) != gamma_12(E_02)^dag
 * -> merging0 violated -> merge_cond_max grows. Both are mutation-proven (perturb ->
 * RED on the relevant guard -> restore); here we ASSERT the guard FIRES (the
 * perturbed run is a positive teeth check, not a regression). */
static void test_c3_merging_teeth(void)
{
    printf("C3 lem_merging condition teeth (A=M_3, perturb gamma to violate "
           "merging0/3):\n");
    const slong n = 3, n1 = 2, n2 = 1;

    /* Helper builds the four natural gamma arrays into g[2][2] (caller frees). */
    slong off[2] = {0, n1};
    slong sz[2]  = {n1, n2};

    for (int leg = 0; leg < 2; leg++) {
        aic_ecstar A;
        aic_cstar_matrix_algebra(&A, n, CPREC);
        acb_mat_t P1, P2;
        diag_proj_range(P1, n, 0, 2);
        diag_proj_range(P2, n, 2, 3);

        acb_mat_t *g[2][2];
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                slong rows = sz[j], cols = sz[k];
                g[j][k] = flint_malloc((size_t) (rows * cols) * sizeof(acb_mat_t));
                for (slong a = 0; a < rows; a++)
                    for (slong b = 0; b < cols; b++) {
                        acb_mat_t *E = &g[j][k][a * cols + b];
                        acb_mat_init(*E, n, n);
                        acb_mat_zero(*E);
                        acb_set_si(acb_mat_entry(*E, a + off[j], b + off[k]), 1);
                    }
            }

        const char *what;
        if (leg == 0) {
            /* (a) scale gamma_12's only image (E_02) by 0.5 -> merging3 + collapse. */
            what = "merging3 (scale gamma_12 by 0.5)";
            acb_mat_scalar_mul_2exp_si(g[0][1][0], g[0][1][0], -1); /* x 1/2 */
        } else {
            /* (b) gamma_21(E_20) := E_21 (not = gamma_12(E_02)^dag = E_20) ->
             * merging0 broken (and the routing still consistent: still in block
             * (2,1) shape, just a wrong image). */
            what = "merging0 (gamma_21(E_20) != gamma_12(E_02)^dag)";
            acb_mat_zero(g[1][0][0]);
            acb_set_si(acb_mat_entry(g[1][0][0], 2, 1), 1); /* E_21 instead of E_20 */
        }

        aic_dhom_B Bm;
        aic_dhom_v vm;
        arb_t mult_def, smin, mcm;
        arb_init(mult_def); arb_init(smin); arb_init(mcm);
        aic_cstar_lem_merging(&Bm, &vm, mult_def, smin, mcm, /*two_block=*/0, n1, n2,
                              (const acb_mat_t *) g[0][0], (const acb_mat_t *) g[0][1],
                              (const acb_mat_t *) g[1][0], (const acb_mat_t *) g[1][1],
                              &A, P1, P2, CPREC);
        double mc = dd(mcm), sm = dd(smin);
        printf("  leg %d [%s]: merge_cond_max=%.3f (>0.1), sigma_min=%.4f\n",
               leg, what, mc, sm);

        /* Both legs: the perturbation must make merge_cond_max grow well past the
         * eta=0 machine-zero baseline. */
        AIC_CHECK_MSG(mc > 0.1, "C3 leg %d: merge_cond_max=%.3e did not grow under "
                      "perturbation (teeth blind)", leg, mc);
        if (leg == 0) {
            /* (a) the off-diagonal direction half-collapses: sigma_min < 1 - eps. */
            AIC_CHECK_MSG(sm < 0.95, "C3 leg 0: sigma_min=%.4f did not drop under the "
                          "gamma_12 scaling (collapse teeth blind)", sm);
        }

        arb_clear(mcm); arb_clear(smin); arb_clear(mult_def);
        aic_dhom_v_clear(&vm); aic_dhom_B_clear(&Bm);
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                slong rows = sz[j], cols = sz[k];
                for (slong i = 0; i < rows * cols; i++) acb_mat_clear(g[j][k][i]);
                flint_free(g[j][k]);
            }
        acb_mat_clear(P2); acb_mat_clear(P1);
        aic_ecstar_clear(&A);
    }
}

int main(void)
{
    test_c1_eta0_oracle();
    test_c2_merge_sum_consistency();
    test_c3_merging_teeth();
    aic_test_report("test_cstar_merging");
    return 0;
}
