/* test_shim_factorize_artifacts.c — C tests for the flat-double ccall shims C4
 * (aic_factorize_artifacts_sizes_d) and C5 (aic_factorize_artifacts_d) of bead
 * aic-exa.3 [C] (docs/research/julia_package_design.md §4.3/§4.4, Appendix B1).
 * Steady-state Julia-free: drives the C shims DIRECTLY with flat double arrays (the
 * SAME ABI the Julia `factorize` ccall uses).
 *
 * OBLIGATIONS (the prompt's C-TEST OBLIGATIONS a/b/d):
 *  (a) eta=0 ORACLE: identity / exact-idempotent Phi -> C5 delups/upsdel brackets
 *      straddle 0 with hi at machine-eps; the C4 B-block sizes d_l MATCH
 *      aic_idemp_wedderburn on the same channel (Rule 6 canonical cross-check).
 *  (b) eta>0 OBLIQUE fixture (make_mixconj; the eta=0 identity oracle is
 *      STRUCTURALLY BLIND to star/dual bugs — §C2/§C13): exercise C5, then
 *      MUTATION-PROVE B1 by SWAPPING dec_kraus and enc_kraus so the
 *      decode-then-encode round-trip defect jumps O(1) -> RED, then restore. This is
 *      the test that proves the Appendix-B1 channel-direction binding (Dec=Delta*,
 *      Enc=Upsilon*; aic_factorize.h:306-316).
 *  (d) DIM-INDEPENDENCE canary (FINDINGS §D2): the round-trip / cb brackets do not
 *      grow with n across n in {4,5,6} mixconj fixtures (c = upsdel/eta bounded).
 *
 * THE DUAL CHANNELS (Appendix B1, the §C13 dual-direction trap). C5 emits dec_*
 * from aic_factorize_dec_kraus (Dec=Delta*, channel B(H)->B, n_B x N) and enc_*
 * from aic_factorize_enc_kraus (Enc=Upsilon*, channel B->B(H), N x n_B). The state
 * round-trip decode-then-encode = Dec o Enc on B-states ~ 1_B; its observable dual
 * is ||Upsilon Delta - 1_B||_cb = upsdel (.tex:2739). <=200 LOC (Rule 10).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_factorize_artifacts_shim.h"
#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"
#include "test_idemp.h"   /* make_mixconj, make_block_cond_exp, make_dephasing */

#define CPREC 256

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

/* Load rK Kraus (each row x col) from a flat K-major [a*(row*col)+i*col+j] array. */
static acb_mat_t *load_kraus(const double *re, const double *im, int rK, int row,
                             int col)
{
    acb_mat_t *K = malloc((size_t) rK * sizeof(acb_mat_t));
    for (int a = 0; a < rK; a++) {
        acb_mat_init(K[a], row, col);
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++) {
                int o = a * (row * col) + i * col + j;
                acb_set_d_d(acb_mat_entry(K[a], i, j), re[o], im[o]);
            }
    }
    return K;
}

static void free_kraus(acb_mat_t *K, int rK)
{
    for (int a = 0; a < rK; a++) acb_mat_clear(K[a]);
    free(K);
}

/* decode-then-encode on B-states: Enc (rE Kraus, each N x n_B) then Dec (rD Kraus,
 * each n_B x N); worst over the IN-BLOCK units of ||(Dec o Enc)(E_pq) - E_pq||_op.
 *
 * THE §C14(b) TARGET IS THE CONDITIONAL EXPECTATION P_B, NOT THE FULL IDENTITY.
 * Dec o Enc = (Upsilon Delta)* (dual; aic_factorize_dual.c) and Upsilon Delta = P_B
 * (block-diagonal conditional expectation), NOT 1_{M_{n_B}}: Delta reads only B's
 * block-diagonal coordinates, so (Dec o Enc)(E_pq) = 0 for OFF-block E_pq while
 * 1_B(E_pq) = E_pq there -> comparing to the full identity is the §C14(b) "test that
 * cannot pass" (reads 1.0 on every off-block unit even at eta=0). So we sweep only
 * IN-BLOCK units (p,q in the same B-block), the regime where Upsilon Delta = id_B.
 * `blocks`/`m` give the B-block sizes; pass m=1 (single block) to sweep all units.
 * Explicit per-set (row,col) so the (b) SWAP can reinterpret the flat arrays. */
static double roundtrip(const double *de_re, const double *de_im, int rD, int drow,
                        int dcol, const double *en_re, const double *en_im, int rE,
                        int erow, int ecol, int nB, int N, const int *blocks, int m)
{
    acb_mat_t *D = load_kraus(de_re, de_im, rD, drow, dcol);
    acb_mat_t *E = load_kraus(en_re, en_im, rE, erow, ecol);
    /* block index of each row/col coordinate (block-diagonal layout). */
    int *blkof = calloc((size_t) (nB), sizeof(int));
    {
        int off = 0;
        for (int l = 0; l < m; l++) {
            for (int i = 0; i < blocks[l] && off < nB; i++) blkof[off++] = l;
        }
        for (; off < nB; off++) blkof[off] = m; /* m==1 path: all same "block" 0 */
        if (m == 1) for (int i = 0; i < nB; i++) blkof[i] = 0;
    }
    double worst = 0.0;
    for (int p = 0; p < nB; p++)
        for (int q = 0; q < nB; q++) {
            if (blkof[p] != blkof[q]) continue; /* §C14(b): in-block units only */
            acb_mat_t rho, inner, out, diff;
            acb_mat_init(rho, nB, nB);
            acb_mat_zero(rho);
            acb_set_si(acb_mat_entry(rho, p, q), 1);
            acb_mat_init(inner, N, N);
            acb_mat_zero(inner);
            for (int a = 0; a < rE; a++) { /* E rho E^dag : N x N */
                acb_mat_t t1, t1d, tmp;
                acb_mat_init(t1, N, nB);
                acb_mat_init(t1d, nB, N);
                acb_mat_init(tmp, N, N);
                acb_mat_mul(t1, E[a], rho, CPREC);
                acb_mat_conjugate_transpose(t1d, E[a]);
                acb_mat_mul(tmp, t1, t1d, CPREC);
                acb_mat_add(inner, inner, tmp, CPREC);
                acb_mat_clear(t1);
                acb_mat_clear(t1d);
                acb_mat_clear(tmp);
            }
            acb_mat_init(out, nB, nB);
            acb_mat_zero(out);
            for (int b = 0; b < rD; b++) { /* D inner D^dag : n_B x n_B */
                acb_mat_t t2, t2d, tmp;
                acb_mat_init(t2, nB, N);
                acb_mat_init(t2d, N, nB);
                acb_mat_init(tmp, nB, nB);
                acb_mat_mul(t2, D[b], inner, CPREC);
                acb_mat_conjugate_transpose(t2d, D[b]);
                acb_mat_mul(tmp, t2, t2d, CPREC);
                acb_mat_add(out, out, tmp, CPREC);
                acb_mat_clear(t2);
                acb_mat_clear(t2d);
                acb_mat_clear(tmp);
            }
            acb_mat_init(diff, nB, nB);
            acb_mat_sub(diff, out, rho, CPREC);
            arb_t e;
            arb_init(e);
            aic_mat_opnorm(e, diff, CPREC);
            double v = arf_get_d(arb_midref(e), ARF_RND_NEAR);
            if (v > worst) worst = v;
            arb_clear(e);
            acb_mat_clear(diff);
            acb_mat_clear(out);
            acb_mat_clear(inner);
            acb_mat_clear(rho);
        }
    free(blkof);
    free_kraus(D, rD);
    free_kraus(E, rE);
    return worst;
}

/* (a) eta=0 oracle: C5 delups/upsdel straddle 0; C4 block sizes match Wedderburn. */
static void t_eta0_oracle(const char *name, aic_ucp_kraus *phi)
{
    double *re, *im;
    int n, r;
    flat(phi, &re, &im, &n, &r);
    int N, nB, dimB, m, rDec, rEnc;
    int *blocks = calloc((size_t) n, sizeof(int));
    aic_factorize_artifacts_sizes_d(&N, &nB, &dimB, &m, blocks, &rDec, &rEnc,
                                    re, im, n, r, 0.0, CPREC);

    double *dre = calloc((size_t) rDec * nB * N, sizeof(double));
    double *dim_ = calloc((size_t) rDec * nB * N, sizeof(double));
    double *ere = calloc((size_t) rEnc * N * nB, sizeof(double));
    double *eim = calloc((size_t) rEnc * N * nB, sizeof(double));
    double du[2], ud[2], eta[2];
    aic_factorize_artifacts_d(dre, dim_, rDec, ere, eim, rEnc, du, ud, eta,
                              re, im, n, r, 0.0, CPREC);

    AIC_CHECK_MSG(du[0] <= 0.0 && du[1] >= 0.0,
                  "%s C5: delups=[%.3e,%.3e] does not straddle 0", name, du[0], du[1]);
    AIC_CHECK_MSG(ud[0] <= 0.0 && ud[1] >= 0.0,
                  "%s C5: upsdel=[%.3e,%.3e] does not straddle 0", name, ud[0], ud[1]);
    AIC_CHECK_MSG(du[1] < 1e-9 && ud[1] < 1e-9,
                  "%s C5: eta=0 brackets hi not ~0 (du=%.3e ud=%.3e)", name, du[1],
                  ud[1]);

    /* round-trip with the CORRECT binding ~ 0 at eta=0 (in-block units, §C14(b)). */
    double good = roundtrip(dre, dim_, rDec, nB, N, ere, eim, rEnc, N, nB, nB, N,
                            blocks, m);
    AIC_CHECK_MSG(good < 1e-6, "%s: eta=0 decode-encode roundtrip=%.3e not ~0", name,
                  good);

    /* Wedderburn B-block cross-check (Rule 6). */
    aic_idemp_decomp dc;
    aic_idemp_decompose(&dc, phi, CPREC);
    aic_idemp_wedderburn wd;
    aic_idemp_wedderburn_decompose(&wd, &dc, CPREC);
    AIC_CHECK_MSG((int) wd.num_blocks == m,
                  "%s C4: num_blocks=%d != Wedderburn %ld", name, m,
                  (long) wd.num_blocks);
    long got[32], want[32];
    for (int l = 0; l < m; l++) got[l] = blocks[l];
    for (slong l = 0; l < wd.num_blocks; l++) want[l] = (long) wd.dim_L[l];
    qsort(got, (size_t) m, sizeof(long), cmp_long);
    qsort(want, (size_t) wd.num_blocks, sizeof(long), cmp_long);
    for (int l = 0; l < m; l++)
        AIC_CHECK_MSG(got[l] == want[l],
                      "%s C4: block d_l mismatch at %d (%ld != Wedderburn %ld)",
                      name, l, got[l], want[l]);
    printf("(a) %s: N=%d nB=%d m=%d rDec=%d rEnc=%d blocks match Wedderburn; "
           "delups_hi=%.2e upsdel_hi=%.2e roundtrip=%.2e\n", name, N, nB, m, rDec,
           rEnc, du[1], ud[1], good);

    aic_idemp_wedderburn_clear(&wd);
    aic_idemp_clear(&dc);
    free(dre);
    free(dim_);
    free(ere);
    free(eim);
    free(blocks);
    free(re);
    free(im);
}

/* (b) eta>0 oblique mutation-proof + (d) dim-independence canary. */
static void t_oblique_mutation(int d, int mblk, double t, double *c_out)
{
    aic_ucp_kraus phi;
    make_mixconj(&phi, d, mblk, t, CPREC); /* oblique eta>0, N=d, n_B=mblk */
    double *re, *im;
    int n, r;
    flat(&phi, &re, &im, &n, &r);
    int N, nB, dimB, m, rDec, rEnc;
    int *blocks = calloc((size_t) n, sizeof(int));
    aic_factorize_artifacts_sizes_d(&N, &nB, &dimB, &m, blocks, &rDec, &rEnc,
                                    re, im, n, r, -2.0, CPREC);

    double *dre = calloc((size_t) rDec * nB * N, sizeof(double));
    double *dim_ = calloc((size_t) rDec * nB * N, sizeof(double));
    double *ere = calloc((size_t) rEnc * N * nB, sizeof(double));
    double *eim = calloc((size_t) rEnc * N * nB, sizeof(double));
    double du[2], ud[2], eta[2];
    aic_factorize_artifacts_d(dre, dim_, rDec, ere, eim, rEnc, du, ud, eta,
                              re, im, n, r, -2.0, CPREC);
    double eta_pr = eta[0];

    /* CORRECT binding: decode-then-encode ~ O(eta). mixconj is single-block (m==1)
     * so P_B = full identity on M_{n_B} and all units are in-block. */
    double good = roundtrip(dre, dim_, rDec, nB, N, ere, eim, rEnc, N, nB, nB, N,
                            blocks, m);
    /* MUTATION (B1): SWAP dec<->enc (feed enc as Dec, dec as Enc; both have nB*N
     * entries so the reshape is valid but the CONTENT is wrong). */
    double bad = roundtrip(ere, eim, rEnc, nB, N, dre, dim_, rDec, N, nB, nB, N,
                           blocks, m);

    double c = good / eta_pr; /* round-trip-over-eta constant (O(eta) per instance) */
    *c_out = c;
    /* per-instance O(eta) bound: the documented true c range is ~0.4-12 (the
     * mixconj(4,2) small-block geometry outlier reaches ~12, FINDINGS §D2/§C11;
     * cstar_build T3 uses the same generous c < 20 magnitude bound). */
    AIC_CHECK_MSG(good < 0.5,
                  "(b) mixconj(%d,%d,%.2f): CORRECT roundtrip=%.4e not O(eta) (<0.5)",
                  d, mblk, t, good);
    AIC_CHECK_MSG(bad > 0.5,
                  "(b) mixconj(%d,%d,%.2f): SWAPPED roundtrip=%.4e not O(1) -> the "
                  "B1 dec/enc binding tooth is BLIND (mutation not RED)", d, mblk, t,
                  bad);
    AIC_CHECK_MSG(bad > 20.0 * good,
                  "(b) mixconj(%d,%d,%.2f): SWAP/CORRECT ratio %.1f too small (<20x)",
                  d, mblk, t, bad / good);
    printf("(b) mixconj(%d,%d,%.2f): N=%d nB=%d eta=%.4e CORRECT roundtrip=%.4e "
           "SWAPPED=%.4e ratio=%.1f (B1 mutation RED) | delups=[%.3e,%.3e] "
           "upsdel=[%.3e,%.3e]\n", d, mblk, t, N, nB, eta_pr, good, bad, bad / good,
           du[0], du[1], ud[0], ud[1]);

    free(dre);
    free(dim_);
    free(ere);
    free(eim);
    free(blocks);
    free(re);
    free(im);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    /* silence test_idemp.h's unused static channel builders (-Wunused-function). */
    (void) make_trace_replace;
    (void) make_noiseless_subsystem;
    (void) make_ns_weighted;
    (void) make_identity;
    (void) make_compress_idemp;

    /* (a) eta=0 oracles: identity, block_cond_exp (multi-block), dephasing. */
    aic_ucp_kraus id2;
    aic_ucp_kraus_init(&id2, 2, 2, 1);
    acb_mat_one(id2.K[0]);
    t_eta0_oracle("identity(2)", &id2);
    aic_ucp_kraus_clear(&id2);

    aic_ucp_kraus bce;
    make_block_cond_exp(&bce, 4, 2);
    t_eta0_oracle("block_cond_exp(4,2)", &bce);
    aic_ucp_kraus_clear(&bce);

    aic_ucp_kraus deph;
    make_dephasing(&deph, 3);
    t_eta0_oracle("dephasing(3)", &deph);
    aic_ucp_kraus_clear(&deph);

    /* (b) + (d): eta>0 oblique mutation-proof + dim-independence canary. */
    double cs[3];
    t_oblique_mutation(4, 2, 0.02, &cs[0]);
    t_oblique_mutation(5, 3, 0.02, &cs[1]);
    t_oblique_mutation(6, 3, 0.02, &cs[2]);
    for (int k = 0; k < 3; k++)
        AIC_CHECK_MSG(cs[k] < 5.0, "(d) canary: roundtrip/eta=%.4f at k=%d >= 5.0",
                      cs[k], k);
    AIC_CHECK_MSG(cs[2] < 2.5 * cs[0] + 0.5,
                  "(d) canary: roundtrip/eta grows with n (c[2]=%.4f vs c[0]=%.4f)",
                  cs[2], cs[0]);
    printf("(d) dim-independence: roundtrip/eta = %.4f, %.4f, %.4f (n=4,5,6) — "
           "bounded, no upward trend\n", cs[0], cs[1], cs[2]);

    aic_test_report("test_shim_factorize_artifacts");
    return 0;
}
