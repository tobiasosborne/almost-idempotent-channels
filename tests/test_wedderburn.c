/* test_wedderburn.c — cross-checks for the Artin-Wedderburn / multiplicity
 * decomposition of A = Img Phi (bead aic-ynu, I1; the eta=0 case). Realizes
 * prop_hom_structure (approximate_algebras.tex:259-272) + prop_Gamma (.tex:2106):
 * aic_idemp_wedderburn_decompose recovers num_blocks m, the block (matrix) dims
 * dim L_j, the multiplicity dims dim E_j, and the block-diagonalizing isometries
 * W_j : M -> L_j (x) E_j so that w(B) = sum_j W_j^dag (A_j (x) 1_{E_j}) W_j.
 *
 * WHAT PINS WHAT (the teeth are complementary; no single one pins correctness):
 *   - the EXACT dim identities (sum dim_L*dim_E == dim_M, sum dim_L^2 == dim_A)
 *     constrain the (dim_L, dim_E) split but do NOT determine it alone;
 *   - the INDEPENDENT aic_cstar_build oracle (a wholly different route) pins
 *     num_blocks AND the {dim L_j} multiset — together with the dim identities this
 *     PINS THE SPLIT (a dim_L=4,dim_E=1 vs dim_L=2,dim_E=2 relabel is rejected by
 *     the cstar {dim L_j} match, which the bare w-recon tooth would miss when
 *     dim_E=1 makes the partial trace trivial);
 *   - the w-RECONSTRUCTION tooth then PINS THE TENSOR ALIGNMENT of W_j (an L<->E
 *     factor swap leaves the split sizes alone but breaks the tensor structure).
 * The w-reconstruction tooth, for every basis element B_k of A:
 *   - M_blk(j) = W_j w(B_k) W_j^dag   (db x db, db = dim L_j * dim E_j),
 *   - A_j = Tr_{E_j}(M_blk) / dim E_j (partial trace over the RIGHT/minor factor;
 *     L-MAJOR layout, aic_mat.h), an INDEPENDENT read-off of the block component,
 *   - ASSERT ||M_blk - A_j (x) 1_{E_j}||_op = 0 (the TENSOR-STRUCTURE defect: a
 *     swapped L<->E factor makes this O(1), e.g. 0.707 for the
 *     noiseless_subsystem(2,2) L<->E swap — the mutation proof below),
 *   - ASSERT ||sum_j W_j^dag (A_j (x) 1_{E_j}) W_j - w(B_k)||_op = 0 (round-trip).
 * Plus: W unitary (each W_j W_j^dag = 1, sum_j W_j^dag W_j = 1_M; the production
 * routine self-certifies dims + unitarity + w-reconstruction, so a successful
 * decompose already implies them, but the test re-asserts independently).
 *
 * SCOPE (I1, eta=0): make_block_cond_exp distinct sizes (M_2 (+) M_3) AND equal
 * sizes (M_2 (+) M_2, M_6 (+) M_6 — the random-H_W BLOCKER regression),
 * make_noiseless_subsystem (single block, nontrivial multiplicity), make_identity
 * (single trivial block). The random H_W + re-draw handles equal-size blocks too.
 */
#include <math.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "test_idemp.h"

#define WPREC 128
static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }
static int long_cmp(const void *a, const void *b)
{
    long x = *(const long *) a, y = *(const long *) b;
    return (x > y) - (x < y);
}

/* reshape column k of w (m^2 row-major) into m x m. */
static void w_col(acb_mat_t W, const acb_mat_t w, slong k, slong m)
{
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            acb_set(acb_mat_entry(W, i, j), acb_mat_entry(w, i * m + j, k));
}

/* The w-reconstruction crux tooth (see file banner). Returns the max over k of
 * BOTH the tensor-structure defect and the round-trip residual (caller asserts
 * <= tol). A_j is extracted via partial trace, INDEPENDENT of how W_j was built. */
static void w_reconstruction_defect(arb_t out, const aic_idemp_wedderburn *W,
                                     const aic_idemp_decomp *d, slong prec)
{
    slong m = W->dim_M;
    arb_t e;
    arb_init(e);
    arb_zero(out);
    acb_mat_t wbk, recon, diff;
    acb_mat_init(wbk, m, m);
    acb_mat_init(recon, m, m);
    acb_mat_init(diff, m, m);
    for (slong k = 0; k < W->dim_A; k++) {
        w_col(wbk, d->w, k, m);
        acb_mat_zero(recon);
        for (slong j = 0; j < W->num_blocks; j++) {
            slong dL = W->dim_L[j], dE = W->dim_E[j], db = dL * dE;
            acb_mat_t Wjd, tmp, Mblk, Aj, one_E, AjXI, td, t2, back;
            acb_mat_init(Wjd, m, db);
            acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
            acb_mat_init(tmp, db, m);
            acb_mat_mul(tmp, W->W_j[j], wbk, prec);   /* W_j w(B_k)        */
            acb_mat_init(Mblk, db, db);
            acb_mat_mul(Mblk, tmp, Wjd, prec);        /* W_j w(B_k) W_j^dag */
            /* A_j = Tr_{E_j}(M_blk) / dim E_j (RIGHT factor, L-major). */
            acb_mat_init(Aj, dL, dL);
            aic_mat_partial_trace_right(Aj, Mblk, dL, dE, prec);
            {
                acb_t inv;
                acb_init(inv);
                acb_set_d(inv, 1.0 / (double) dE);
                acb_mat_scalar_mul_acb(Aj, Aj, inv, prec);
                acb_clear(inv);
            }
            acb_mat_init(one_E, dE, dE);
            acb_mat_one(one_E);
            acb_mat_init(AjXI, db, db);
            aic_mat_kronecker(AjXI, Aj, one_E, prec); /* A_j (x) 1_{E_j}   */
            /* tensor-structure defect ||M_blk - A_j (x) 1_E|| */
            acb_mat_init(td, db, db);
            acb_mat_sub(td, Mblk, AjXI, prec);
            aic_mat_opnorm(e, td, prec);
            if (arb_gt(e, out)) arb_set(out, e);
            /* recon += W_j^dag (A_j (x) 1_E) W_j */
            acb_mat_init(t2, m, db);
            acb_mat_mul(t2, Wjd, AjXI, prec);
            acb_mat_init(back, m, m);
            acb_mat_mul(back, t2, W->W_j[j], prec);
            acb_mat_add(recon, recon, back, prec);
            acb_mat_clear(back);
            acb_mat_clear(t2);
            acb_mat_clear(td);
            acb_mat_clear(AjXI);
            acb_mat_clear(one_E);
            acb_mat_clear(Aj);
            acb_mat_clear(Mblk);
            acb_mat_clear(tmp);
            acb_mat_clear(Wjd);
        }
        acb_mat_sub(diff, recon, wbk, prec);
        aic_mat_opnorm(e, diff, prec);            /* round-trip residual */
        if (arb_gt(e, out)) arb_set(out, e);
    }
    acb_mat_clear(diff);
    acb_mat_clear(recon);
    acb_mat_clear(wbk);
    arb_clear(e);
}

/* ||sum_j W_j^dag W_j - 1_M||_op (independent of the production self-cert). */
static void W_unitarity_defect(arb_t out, const aic_idemp_wedderburn *W, slong prec)
{
    slong m = W->dim_M;
    acb_mat_t WtW, one_m, Wjd, prod;
    acb_mat_init(WtW, m, m);
    acb_mat_init(one_m, m, m);
    acb_mat_zero(WtW);
    for (slong j = 0; j < W->num_blocks; j++) {
        slong db = W->dim_L[j] * W->dim_E[j];
        acb_mat_init(Wjd, m, db);
        acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
        acb_mat_init(prod, m, m);
        acb_mat_mul(prod, Wjd, W->W_j[j], prec);
        acb_mat_add(WtW, WtW, prod, prec);
        acb_mat_clear(prod);
        acb_mat_clear(Wjd);
    }
    acb_mat_one(one_m);
    acb_mat_sub(WtW, WtW, one_m, prec);
    aic_mat_opnorm(out, WtW, prec);
    acb_mat_clear(one_m);
    acb_mat_clear(WtW);
}

/* INDEPENDENT ORACLE (aic_cstar_build, th_main): build B = (+)_l M_{|C_l|} from
 * the SAME channel via a wholly different route (the master loop), and assert
 * num_blocks + the {dim L_j} multiset MATCH the Wedderburn block sizes. */
static void cstar_oracle_check(const char *name, aic_ucp_kraus *phi,
                               const aic_idemp_wedderburn *W)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, WPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, 0.0, WPREC);

    AIC_CHECK_MSG(B.num_blocks == W->num_blocks,
                  "%s: cstar_build num_blocks=%ld != Wedderburn %ld", name,
                  (long) B.num_blocks, (long) W->num_blocks);

    long got[32], want[32];
    for (slong l = 0; l < B.num_blocks; l++) got[l] = (long) B.d[l];
    for (slong l = 0; l < W->num_blocks; l++) want[l] = (long) W->dim_L[l];
    qsort(got, (size_t) B.num_blocks, sizeof(long), long_cmp);
    qsort(want, (size_t) W->num_blocks, sizeof(long), long_cmp);
    for (slong l = 0; l < W->num_blocks; l++)
        AIC_CHECK_MSG(got[l] == want[l],
                      "%s: block-size multiset mismatch at %ld (cstar %ld != wedd %ld)",
                      name, (long) l, got[l], want[l]);

    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

/* Decompose, assert the block structure + the w-reconstruction tooth + W unitary
 * + the dim identities + (optionally) the cstar_build oracle. exp_dimL / exp_dimE
 * are SORTED-multiset expectations (length exp_blocks).
 *
 * run_cstar: the independent aic_cstar_build oracle (the main-theorem master loop)
 * is O(1) for small blocks but ~34 s for M_4 (+) M_4 and far worse for M_6 (+) M_6,
 * so it is GATED off for the large equal-size BLOCKER regression (run_cstar=0).
 * That case is still fully pinned WITHOUT it: the EXACT dim identities + the
 * production w-reconstruction self-cert pin the split, and the cstar oracle is
 * exercised on the cheap equal-size case block_cond_exp(4,2) = M_2 (+) M_2. */
static void check_wedderburn(const char *name, aic_ucp_kraus *phi,
                             slong exp_blocks, const long *exp_dimL,
                             const long *exp_dimE, slong exp_dM, slong exp_dA,
                             int run_cstar)
{
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, phi, WPREC);
    aic_idemp_wedderburn W;
    aic_idemp_wedderburn_decompose(&W, &d, WPREC);

    arb_t tol, def, uni;
    arb_init(tol);
    arb_init(def);
    arb_init(uni);
    set_tol(tol, 1e-10);

    AIC_CHECK_MSG(W.num_blocks == exp_blocks, "%s: num_blocks=%ld != %ld", name,
                  (long) W.num_blocks, (long) exp_blocks);
    AIC_CHECK_MSG(W.dim_M == exp_dM, "%s: dim_M=%ld != %ld", name,
                  (long) W.dim_M, (long) exp_dM);
    AIC_CHECK_MSG(W.dim_A == exp_dA, "%s: dim_A=%ld != %ld", name,
                  (long) W.dim_A, (long) exp_dA);

    /* (b) EXACT dim identities */
    slong sumLE = 0, sumL2 = 0;
    long gotL[32], gotE[32];
    for (slong j = 0; j < W.num_blocks; j++) {
        sumLE += W.dim_L[j] * W.dim_E[j];
        sumL2 += W.dim_L[j] * W.dim_L[j];
        gotL[j] = (long) W.dim_L[j];
        gotE[j] = (long) W.dim_E[j];
    }
    AIC_CHECK_MSG(sumLE == exp_dM, "%s: sum dim_L*dim_E=%ld != dim_M=%ld", name,
                  (long) sumLE, (long) exp_dM);
    AIC_CHECK_MSG(sumL2 == exp_dA, "%s: sum dim_L^2=%ld != dim_A=%ld", name,
                  (long) sumL2, (long) exp_dA);

    /* block-size multisets {dim L_j}, {dim E_j} */
    qsort(gotL, (size_t) W.num_blocks, sizeof(long), long_cmp);
    qsort(gotE, (size_t) W.num_blocks, sizeof(long), long_cmp);
    long wantL[32], wantE[32];
    for (slong j = 0; j < exp_blocks; j++) { wantL[j] = exp_dimL[j]; wantE[j] = exp_dimE[j]; }
    qsort(wantL, (size_t) exp_blocks, sizeof(long), long_cmp);
    qsort(wantE, (size_t) exp_blocks, sizeof(long), long_cmp);
    for (slong j = 0; j < exp_blocks; j++) {
        AIC_CHECK_MSG(gotL[j] == wantL[j], "%s: dim_L multiset mismatch at %ld (%ld != %ld)",
                      name, (long) j, gotL[j], wantL[j]);
        AIC_CHECK_MSG(gotE[j] == wantE[j], "%s: dim_E multiset mismatch at %ld (%ld != %ld)",
                      name, (long) j, gotE[j], wantE[j]);
    }

    /* THE CRUX TOOTH: w-reconstruction (tensor structure + round-trip). */
    w_reconstruction_defect(def, &W, &d, WPREC);
    AIC_CHECK_MSG(arb_le(def, tol), "%s: w-reconstruction defect %.3e > 1e-10", name,
                  dd(def));

    /* (a) W unitary */
    W_unitarity_defect(uni, &W, WPREC);
    AIC_CHECK_MSG(arb_le(uni, tol), "%s: ||sum W_j^dag W_j - 1_M|| %.3e > 1e-10",
                  name, dd(uni));

    /* (c) independent aic_cstar_build oracle (gated; see banner) */
    if (run_cstar) cstar_oracle_check(name, phi, &W);

    printf("  %-28s blocks=%ld dimL=[", name, (long) W.num_blocks);
    for (slong j = 0; j < W.num_blocks; j++)
        printf("%ld%s", (long) W.dim_L[j], j + 1 < W.num_blocks ? "," : "");
    printf("] dimE=[");
    for (slong j = 0; j < W.num_blocks; j++)
        printf("%ld%s", (long) W.dim_E[j], j + 1 < W.num_blocks ? "," : "");
    printf("]  w-recon=%.2e  W-unit=%.2e%s\n", dd(def), dd(uni),
           run_cstar ? "" : "  (cstar-oracle gated: cost)");

    arb_clear(uni);
    arb_clear(def);
    arb_clear(tol);
    aic_idemp_wedderburn_clear(&W);
    aic_idemp_clear(&d);
}

static void test_oracles(void)
{
    printf("Wedderburn decomposition (eta=0 oracle, w-reconstruction tooth):\n");
    aic_ucp_kraus phi;

    /* block_cond_exp(5,2): A = M_2 (+) M_3, distinct sizes, E_j = C^1. */
    long L_23[2] = {2, 3}, E_11[2] = {1, 1};
    make_block_cond_exp(&phi, 5, 2);
    check_wedderburn("block_cond_exp(5,2)", &phi, 2, L_23, E_11, 5, 13, 1);
    aic_ucp_kraus_clear(&phi);

    /* EQUAL-SIZE regression oracles (the aic-ynu BLOCKER: a deterministic H_W
     * manufactures a kernel spanning the two equal blocks -> spurious abort). The
     * random H_W + re-draw must give the correct split. block_cond_exp(2k,k) =
     * M_k (+) M_k, E_j = C^1. (Onset of the deterministic >=2-dim kernel is
     * MEASURED at dim_M=12 / M_6 (+) M_6: M_2..M_5 (+) themselves do not collide
     * under the cos-seed, so M_6 (+) M_6 is the load-bearing FIX-1 pin.) */
    long L_22[2] = {2, 2}, L_66[2] = {6, 6};
    /* M_2 (+) M_2 (dim_M=4): cheap equal-size case; cstar oracle ON (0.13 s). */
    make_block_cond_exp(&phi, 4, 2);
    check_wedderburn("block_cond_exp(4,2)", &phi, 2, L_22, E_11, 4, 8, 1);
    aic_ucp_kraus_clear(&phi);
    /* M_6 (+) M_6 (dim_M=12): the EXACT case that broke the deterministic H_W.
     * decompose+arb-certify is 0.18 s; the cstar oracle alone would be many tens of
     * seconds (it is ~34 s already at M_4 (+) M_4), so it is GATED OFF here — the
     * dim identities + the production w-reconstruction self-cert pin the split. */
    make_block_cond_exp(&phi, 12, 6);
    check_wedderburn("block_cond_exp(12,6)", &phi, 2, L_66, E_11, 12, 72, 0);
    aic_ucp_kraus_clear(&phi);

    /* noiseless_subsystem(2,2): A = M_2, single block, E_1 = C^2 (nontrivial). */
    long L_2[1] = {2}, E_2[1] = {2};
    make_noiseless_subsystem(&phi, 2, 2);
    check_wedderburn("noiseless_subsystem(2,2)", &phi, 1, L_2, E_2, 4, 4, 1);
    aic_ucp_kraus_clear(&phi);

    /* identity(3): A = M_3, single trivial block, E_1 = C^1. */
    long L_3[1] = {3}, E_1[1] = {1};
    make_identity(&phi, 3);
    check_wedderburn("identity(3)", &phi, 1, L_3, E_1, 3, 9, 1);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_oracles();
    aic_test_report("test_wedderburn");
    printf("OK test_wedderburn\n");
    return 0;
}
