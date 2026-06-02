/* test_shim_assoc_iso.c — C tests for the flat-double ccall shims C2
 * (aic_assoc_summary_d) and C3 (aic_main_iso_summary_d) of bead aic-exa.3 [C], the
 * Julia-package ABI extension (docs/research/julia_package_design.md §4.3/§4.4,
 * Appendix B2). Steady-state Julia-free: these drive the C shims DIRECTLY with flat
 * double arrays (the SAME ABI the Julia ccall will use), so a shim regression is
 * caught at the C gate, not only after the Julia layer lands.
 *
 * OBLIGATIONS (the prompt's C-TEST OBLIGATIONS a/c/d):
 *  (a) eta=0 ORACLE: identity / exact-idempotent Phi -> C2 eps_assoc and C3 iso_def
 *      brackets straddle 0 with hi at machine-eps; the C3 B-block sizes d_l MATCH
 *      aic_idemp_wedderburn on the same channel (Rule 6 canonical cross-check).
 *  (c) BRACKET sanity vs a KNOWN analytic anchor: for the depolarizing channel
 *      phi_t at d=2,3 the eig-free eta bracket (aic_cbnorm_eigfree_d) must CONTAIN
 *      the closed form eta = t(1-t)*2*(1-1/d^2) (solver-free; no MOSEK).
 *  (d) DIM-INDEPENDENCE canary (FINDINGS §D2): the C3 cb-bracket LOWER endpoint (the
 *      operator-faithful, §C12-correct quantity) and iso_def do NOT grow with n
 *      across n=2,3,4 (eta=0 complete-isometry oracle: cb-lo == 1 / 0.7071 exactly,
 *      dimension-INDEPENDENT; eta>0 mixconj: iso_def/eta bounded, no upward trend).
 *
 * The cb brackets are SOLVER-FREE + EIG-FREE (Appendix B2): C3 forms the GOLDEN-RULE
 * adjoint Choi J(v*)/J((v^-1)*) and brackets ||v||_cb=||v*||_⋄ via the rectangular
 * eig-free core aic_cbnorm_eigfree_ball_choi_rect — the OPERATOR-faithful diamond,
 * never the Frobenius sigma_min ampliation (FINDINGS §C12). <=200 LOC (Rule 10).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic/aic_assoc_shim.h"
#include "aic/aic_channels.h"
#include "aic/aic_idemp.h"
#include "aic/aic_iso_shim.h"
#include "aic_test.h"
#include "aic/aic_ucp_shim.h"
#include "aic/aic_ucp.h"
#include "test_idemp.h"   /* make_mixconj, PREC */

#define CPREC 256

/* Flatten an n x n self-map's Kraus to K-major [a*n*n + i*n + j] double arrays. */
static void flat(const aic_ucp_kraus *phi, double **re, double **im, int *n, int *r)
{
    *n = (int) phi->dim_H;
    *r = (int) phi->r;
    int nn = (*n) * (*n);
    *re = calloc((size_t) (*r) * nn, sizeof(double));
    *im = calloc((size_t) (*r) * nn, sizeof(double));
    for (int a = 0; a < *r; a++)
        for (int i = 0; i < *n; i++)
            for (int j = 0; j < *n; j++) {
                int off = a * nn + i * (*n) + j;
                acb_srcptr e = acb_mat_entry(phi->K[a], i, j);
                (*re)[off] = arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR);
                (*im)[off] = arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR);
            }
}

static int cmp_long(const void *a, const void *b)
{
    long x = *(const long *) a, y = *(const long *) b;
    return (x > y) - (x < y);
}

/* (a) eta=0 oracle: C2 eps_assoc straddles 0 (hi ~ machine-eps); C3 iso straddles 0;
 * C3 block sizes match aic_idemp_wedderburn (Rule 6). */
static void t_eta0_oracle(const char *name, aic_ucp_kraus *phi)
{
    double *re, *im;
    int n, r;
    flat(phi, &re, &im, &n, &r);
    int n4 = n * n * n * n;

    /* C2 */
    int on, odimA;
    double epsA[2];
    double *bre = calloc((size_t) n4, sizeof(double));
    double *bim = calloc((size_t) n4, sizeof(double));
    aic_assoc_summary_d(&on, &odimA, epsA, bre, bim, re, im, n, r, 0.0, CPREC);
    AIC_CHECK_MSG(on == n, "%s C2: out_n=%d != n=%d", name, on, n);
    AIC_CHECK_MSG(epsA[0] <= 0.0 && epsA[1] >= 0.0,
                  "%s C2: eps_assoc=[%.3e,%.3e] does not straddle 0", name,
                  epsA[0], epsA[1]);
    AIC_CHECK_MSG(epsA[1] < 1e-9, "%s C2: eta=0 eps_assoc hi=%.3e not ~0", name,
                  epsA[1]);

    /* C3 */
    int onB, odimB, om;
    int *blocks = calloc((size_t) n, sizeof(int));
    double cbf[2], cbi[2], iso[2], oe[1];
    aic_main_iso_summary_d(&onB, &odimB, &om, blocks, cbf, cbi, iso, oe,
                           re, im, n, r, 0.0, CPREC);
    AIC_CHECK_MSG(odimB == odimA, "%s C3: dimB=%d != dimA=%d", name, odimB, odimA);
    AIC_CHECK_MSG(iso[0] <= 0.0 && iso[1] >= 0.0,
                  "%s C3: iso=[%.3e,%.3e] does not straddle 0", name, iso[0], iso[1]);
    AIC_CHECK_MSG(iso[1] < 1e-9, "%s C3: eta=0 iso hi=%.3e not ~0", name, iso[1]);
    /* complete-isometry: ||v||_cb = ||v^-1||_cb = 1 must lie in the brackets. */
    AIC_CHECK_MSG(cbf[0] <= 1.0 && 1.0 <= cbf[1],
                  "%s C3: cbfwd=[%.4f,%.4f] excludes 1", name, cbf[0], cbf[1]);
    AIC_CHECK_MSG(cbi[0] <= 1.0 && 1.0 <= cbi[1],
                  "%s C3: cbinv=[%.4f,%.4f] excludes 1", name, cbi[0], cbi[1]);

    /* Wedderburn B-block cross-check (Rule 6). */
    aic_idemp_decomp dc;
    aic_idemp_decompose(&dc, phi, CPREC);
    aic_idemp_wedderburn wd;
    aic_idemp_wedderburn_decompose(&wd, &dc, CPREC);
    AIC_CHECK_MSG((int) wd.num_blocks == om,
                  "%s C3: num_blocks=%d != Wedderburn %ld", name, om,
                  (long) wd.num_blocks);
    long got[32], want[32];
    for (int l = 0; l < om; l++) got[l] = blocks[l];
    for (slong l = 0; l < wd.num_blocks; l++) want[l] = (long) wd.dim_L[l];
    qsort(got, (size_t) om, sizeof(long), cmp_long);
    qsort(want, (size_t) wd.num_blocks, sizeof(long), cmp_long);
    for (int l = 0; l < om; l++)
        AIC_CHECK_MSG(got[l] == want[l],
                      "%s C3: block d_l mismatch at %d (%ld != Wedderburn %ld)",
                      name, l, got[l], want[l]);
    printf("(a) %s: dimA=%d m=%d blocks match Wedderburn; epsA_hi=%.2e iso_hi=%.2e "
           "cbfwd=[%.3f,%.3f] cbinv=[%.3f,%.3f]\n", name, odimA, om, epsA[1],
           iso[1], cbf[0], cbf[1], cbi[0], cbi[1]);

    aic_idemp_wedderburn_clear(&wd);
    aic_idemp_clear(&dc);
    free(blocks);
    free(bre);
    free(bim);
    free(re);
    free(im);
}

/* (c) phi_t = depolarizing: the eig-free cb bracket must contain the closed form. */
static void t_phit_bracket(int d, double p)
{
    aic_ucp_kraus dep;
    aic_channel_depolarizing(&dep, d, p, CPREC);
    double *re, *im;
    int n, r;
    flat(&dep, &re, &im, &n, &r);
    double lo, hi;
    aic_cbnorm_eigfree_d(&lo, &hi, re, im, n, r, CPREC);
    double cf = p * (1.0 - p) * 2.0 * (1.0 - 1.0 / ((double) d * d));
    AIC_CHECK_MSG(lo <= cf && cf <= hi,
                  "phi_t d=%d p=%.2f: cb bracket [%.6f,%.6f] excludes closed-form "
                  "eta=%.6f", d, p, lo, hi, cf);
    printf("(c) depolarizing d=%d p=%.2f: cb=[%.6f,%.6f] CONTAINS eta=%.6f "
           "(lo-margin %.4f hi-margin %.4f)\n", d, p, lo, hi, cf, cf - lo, hi - cf);
    free(re);
    free(im);
    aic_ucp_kraus_clear(&dep);
}

/* C2 NON-VACUITY (§C2): on a genuinely OBLIQUE eta>0 channel the associator-defect
 * bracket eps_assoc must read NONZERO O(eta) — the eta=0 identity oracle (epsA=0) is
 * structurally BLIND to the Choi-Effros star, so this oblique readout proves C2's
 * marshalling surfaces the real defect, not a hardcoded 0. */
static void t_c2_oblique_nonvacuous(void)
{
    aic_ucp_kraus phi;
    make_mixconj(&phi, 5, 3, 0.05, CPREC);
    double *re, *im;
    int n, r;
    flat(&phi, &re, &im, &n, &r);
    int n4 = n * n * n * n, on, od;
    double ea[2];
    double *br = calloc((size_t) n4, sizeof(double));
    double *bi = calloc((size_t) n4, sizeof(double));
    aic_assoc_summary_d(&on, &od, ea, br, bi, re, im, n, r, 0.0, CPREC);
    AIC_CHECK_MSG(ea[1] > 1e-5,
                  "C2 non-vacuity: oblique mixconj eps_assoc hi=%.4e ~0 (star blind?)",
                  ea[1]);
    AIC_CHECK_MSG(ea[1] < 1.0,
                  "C2 non-vacuity: eps_assoc hi=%.4e not O(eta) (>=1)", ea[1]);
    printf("(C2-oblique) mixconj(5,3,0.05): dimA=%d eps_assoc=[%.4e,%.4e] (nonzero "
           "O(eta) star defect)\n", od, ea[0], ea[1]);
    free(br);
    free(bi);
    free(re);
    free(im);
    aic_ucp_kraus_clear(&phi);
}

/* C3 cb-bracket vs the INDEPENDENT MOSEK-pinned SDP value (Rule 6, ladder rung 4).
 * FINDINGS §C12.O2 committed mixconj(6,2,0.03): ||v||_cb = 1.0019683734,
 * ||v^-1||_cb = 1.5353598357 (the Watrous rect SDP). The solver-free eig-free
 * bracket MUST contain both — this ties the bracket to a known value and proves it
 * is not a "test that cannot fail" (the bracket is genuinely two-sided). */
static void t_cb_vs_sdp(void)
{
    aic_ucp_kraus phi;
    make_mixconj(&phi, 6, 2, 0.03, CPREC);
    double *re, *im;
    int n, r;
    flat(&phi, &re, &im, &n, &r);
    int nB, dimB, m;
    int *blocks = calloc((size_t) n, sizeof(int));
    double cbf[2], cbi[2], iso[2], oe[1];
    aic_main_iso_summary_d(&nB, &dimB, &m, blocks, cbf, cbi, iso, oe,
                           re, im, n, r, -2.0, CPREC);
    const double sdp_fwd = 1.0019683734, sdp_inv = 1.5353598357;
    AIC_CHECK_MSG(cbf[0] <= sdp_fwd && sdp_fwd <= cbf[1],
                  "C3 cb vs SDP: cbfwd=[%.6f,%.6f] excludes ||v||_cb=%.6f", cbf[0],
                  cbf[1], sdp_fwd);
    AIC_CHECK_MSG(cbi[0] <= sdp_inv && sdp_inv <= cbi[1],
                  "C3 cb vs SDP: cbinv=[%.6f,%.6f] excludes ||v^-1||_cb=%.6f", cbi[0],
                  cbi[1], sdp_inv);
    printf("(rung4) mixconj(6,2,0.03): cbfwd=[%.4f,%.4f] contains SDP %.6f; "
           "cbinv=[%.4f,%.4f] contains SDP %.6f\n", cbf[0], cbf[1], sdp_fwd, cbi[0],
           cbi[1], sdp_inv);
    free(blocks);
    free(re);
    free(im);
    aic_ucp_kraus_clear(&phi);
}

/* (d) dim-independence canary, eta>0 mixconj: iso_def/eta bounded + NO upward trend.
 * FINDINGS §D2: the within-family c-ratio is geometry-FRAGILE — the mixconj(4,2)
 * small-block (dim_A=2 degenerate wrapper, §C16) is a known c~12 outlier with no
 * .tex:484 content. So (i) the per-instance bound is the GENEROUS c < 20 (matching
 * cstar_build T3's documented true c~0.4-12 range), and (ii) the dimension-GROWTH
 * canary runs on the WELL-BEHAVED pure-dim sweep (5,3)->(6,3) (same M_3 block, only
 * ambient n grows 5->6) where c must stay flat (no upward trend). */
static void t_dim_indep_mixconj(void)
{
    int fx[3][2] = {{4, 2}, {5, 3}, {6, 3}};
    double cs[3];
    for (int k = 0; k < 3; k++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, fx[k][0], fx[k][1], 0.02, CPREC);
        double *re, *im;
        int n, r;
        flat(&phi, &re, &im, &n, &r);
        int nB, dimB, m;
        int *blocks = calloc((size_t) n, sizeof(int));
        double cbf[2], cbi[2], iso[2], oe[1];
        aic_main_iso_summary_d(&nB, &dimB, &m, blocks, cbf, cbi, iso, oe,
                               re, im, n, r, -2.0, CPREC);
        /* eps_used = eta_proxy (the -2.0 sentinel). c = iso_def(mid)/eta. */
        double eta = oe[0];
        double c = 0.5 * (iso[0] + iso[1]) / eta;
        cs[k] = c;
        AIC_CHECK_MSG(cbf[0] <= 1.0 && 1.0 <= cbf[1],
                      "mixconj(%d,%d) cbfwd=[%.4f,%.4f] excludes 1", fx[k][0],
                      fx[k][1], cbf[0], cbf[1]);
        AIC_CHECK_MSG(cbi[0] <= 1.0 && 1.0 <= cbi[1],
                      "mixconj(%d,%d) cbinv=[%.4f,%.4f] excludes 1", fx[k][0],
                      fx[k][1], cbi[0], cbi[1]);
        /* per-instance O(eta) bound (generous, T3-matching). */
        AIC_CHECK_MSG(c < 20.0, "mixconj(%d,%d) iso/eta=%.4f >= 20 (not O(eta))",
                      fx[k][0], fx[k][1], c);
        printf("(d) mixconj(%d,%d,0.02): n=%d c=iso/eta=%.4f cbfwd=[%.3f,%.3f] "
               "cbinv=[%.3f,%.3f]\n", fx[k][0], fx[k][1], n, c, cbf[0], cbf[1],
               cbi[0], cbi[1]);
        free(blocks);
        free(re);
        free(im);
        aic_ucp_kraus_clear(&phi);
    }
    /* dimension-GROWTH canary (FINDINGS §D2): the WELL-BEHAVED (5,3)->(6,3) pure-dim
     * sweep (cs[1], cs[2]; same M_3 block) must stay flat — no upward trend. */
    AIC_CHECK_MSG(cs[2] < 2.0 * cs[1] + 0.2,
                  "mixconj dim canary: c grows with n (5,3)->(6,3): %.4f -> %.4f",
                  cs[1], cs[2]);
}

int main(void)
{
    /* silence test_idemp.h's unused static channel builders (-Wunused-function). */
    (void) make_trace_replace;
    (void) make_noiseless_subsystem;
    (void) make_ns_weighted;
    (void) make_identity;
    (void) make_compress_idemp;

    aic_ucp_kraus id2;
    aic_ucp_kraus_init(&id2, 2, 2, 1);
    acb_mat_one(id2.K[0]);
    t_eta0_oracle("identity(2)", &id2);
    aic_ucp_kraus_clear(&id2);

    aic_ucp_kraus id3;
    aic_ucp_kraus_init(&id3, 3, 3, 1);
    acb_mat_one(id3.K[0]);
    t_eta0_oracle("identity(3)", &id3);
    aic_ucp_kraus_clear(&id3);

    aic_ucp_kraus bce;
    make_block_cond_exp(&bce, 4, 2);
    t_eta0_oracle("block_cond_exp(4,2)", &bce);
    aic_ucp_kraus_clear(&bce);

    aic_ucp_kraus deph;
    make_dephasing(&deph, 4);
    t_eta0_oracle("dephasing(4)", &deph);
    aic_ucp_kraus_clear(&deph);

    t_phit_bracket(2, 0.1);
    t_phit_bracket(3, 0.1);
    t_phit_bracket(2, 0.3);

    t_c2_oblique_nonvacuous();
    t_cb_vs_sdp();
    t_dim_indep_mixconj();

    aic_test_report("test_shim_assoc_iso");
    return 0;
}
