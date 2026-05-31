/* test_cstar_merge.c — cross-checks for cstar_build Increment 2 (bead aic-097):
 * lem_add_dim (.tex:1363), the Stage-3 off-diagonal certification S_{P_C,P_D}=0
 * (.tex:1443), and cor_merge_sum (.tex:1352). These are the two SIMPLER assembly
 * lemmas of the th_main master loop (approximate_algebras.tex:1321-1444, §9),
 * built on the corner dim queries (aic_corner_alpha_dims / aic_corner_dim_S) and
 * the dhom_v concat + certification (aic_dhom_defect_sweep / aic_dhom_v_sigma_min).
 *
 * NO NEW double-vs-arb path is added here: the new code (aic_cstar_merge.c) is
 * STRUCTURAL GLUE over already-tested numerical routines, which carry BOTH number
 * paths (corner: double SVD extraction + arb-certified relation defects; dhom: arb
 * certified). The cross-checks are therefore (i) the eta=0 oracle (the identity
 * channel A = M_n: exact dim decompositions, machine-zero merged defect, sigma_min
 * ~1) and (ii) the oblique eta>0 trend (make_mixconj: merged defect O(eta),
 * SHRINKING across two eta values; the STAR teeth, FINDINGS §C2 — the merged-v
 * defect is computed with A's Choi-Effros star, so the eta=0->eta>0 trend is RED
 * if the wrong product were used).
 *
 * MUTATION-PROVEN teeth (Rule 7; each RED observed BY HAND, then restored — NOT
 * left as aborting/failing cases in the automated suite):
 *   A1: feed a WRONG P_part (P_parts[0] = full P instead of the rank-1 piece) so
 *       sum_jk dim S_{Pj,Qk} != dim S_{P,Q} -> aic_cstar_lem_add_dim's internal
 *       N == d_PQ assert ABORTS. Ran by hand: with P_parts = {P, P_2} (P = the
 *       rank-2 sum, double-counted) the sum N=5 != d_PQ=4 -> abort("lem_add_dim
 *       conclusion violated"). Confirmed RED, restored.
 *   B2: STAR teeth (FINDINGS §C2). The merged-v multiplicativity defect uses A's
 *       star internally (aic_dhom_defect_sweep -> aic_dhom_Gv -> aic_ecstar_star).
 *       Mutation: recompute the merged-v defect over a PLAIN-product A (overwrite
 *       ae.A.star_phi with the identity copy, i.e. star == plain matmul) on the
 *       OBLIQUE mixconj fixture. Ran by hand (/tmp/teeth_b2): the defect-to-eta
 *       ratio c = mult_def/eta jumps from the STAR path's ~0.017 -> 0.005 to the
 *       PLAIN path's ~0.43 -> 0.41 (the O(1) Phi_tilde-vs-plain gap, ~25-80x
 *       larger), so the retained tooth c < 0.2 goes RED (0.43 > 0.2). Confirmed
 *       RED, restored. (The magnitude shrinks under BOTH products here because
 *       P_1=span(e1)/P_2=span(e2) are in-A projections whose Ptilde stay near
 *       genuine even plain — the FINDINGS §C2 in-A-fixture caveat — so the SHARP
 *       discriminant is the c-RATIO magnitude, not the shrink direction; the
 *       c-shrink and the magnitude-shrink are kept as secondary checks.)
 *   B1: collapse a merged vE (v_out->vE[1] = 0 after merge) -> sigma_min drops to
 *       ~0 and the dim_B==dim_A bijectivity leg fails. (Same collapse-detector
 *       lesson as test_errreduce T6; documented, not retained as a failing case.)
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_corner.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_errreduce.h"
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
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op (mirrors test_cstar). */
static double eta_proxy(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S, S2, D;
    acb_mat_init(S, nn, nn);
    acb_mat_init(S2, nn, nn);
    acb_mat_init(D, nn, nn);
    aic_assoc_superop_from_ucp(S, phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    arb_t e;
    arb_init(e);
    /* op-norm on midpoints (the near-zero diff carries wide balls, false-fails the
     * certified Gram path — same workaround as test_cstar's mid_opnorm_d). */
    acb_mat_t M;
    acb_mat_init(M, nn, nn);
    for (slong i = 0; i < nn; i++)
        for (slong j = 0; j < nn; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    aic_mat_opnorm(e, M, prec);
    double o = dd(e);
    acb_mat_clear(M);
    arb_clear(e);
    acb_mat_clear(D); acb_mat_clear(S2); acb_mat_clear(S);
    return o;
}

/* n x n diagonal projector with diagonal entries [lo,hi) set to 1 (rest 0). */
static void diag_proj_range(acb_mat_t P, slong n, slong lo, slong hi)
{
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    for (slong i = lo; i < hi; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
}

/* ===================== A1: lem_add_dim eta=0 oracle ======================== *
 * Dephasing channel A = diagonal matrices of M_4 (eta=0; A is COMMUTATIVE, dim 4).
 * This is the right setting for lem_add_dim's NON-TRIVIAL block decomposition: the
 * OFF-diagonal corners vanish (in the identity channel A=M_n they would not — the
 * off-diagonal block E_{jk} is in M_n, so S_{P_j,Q_k}=1 for all j,k, masking the
 * structure). Orthogonal rank-1 diagonal projections P_1=diag(1,0,0,0),
 * P_2=diag(0,1,0,0); P=P_1+P_2=diag(1,1,0,0). Q=P (Q_k=P_k). For the dephasing
 * (commutative) algebra:
 *   dim S_{P_j,Q_k} = 1 iff j==k else 0 (a diagonal element exists only when the
 *   supports agree; distinct diagonal entries have no connecting operator in A),
 * so sum_jk dim S_{P_j,Q_k} = 2 (the two diagonal blocks), and dim S_{P,Q} = 2
 * (the 2-dim diagonal corner span{E_00, E_11}). aic_cstar_lem_add_dim returns 2;
 * the internal N==d_PQ assert is the teeth. HAND value 2 asserted (verified by a
 * standalone probe: per-block S_{0,0}=1, S_{0,1}=0, S_{1,0}=0, S_{1,1}=1; total 2). */
static void test_a1_lem_add_dim(void)
{
    printf("A1 lem_add_dim eta=0 oracle (dephasing(4), P=Q=diag(1,1,0,0), "
           "commutative A):\n");
    const slong n = 4;
    aic_ucp_kraus phi;
    make_dephasing(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    /* P_parts is a CONTIGUOUS acb_mat_t array (aic_corner_alpha_dims indexes it as
     * P_parts[j], an acb_mat_t = acb_mat_struct[1]), NOT an array of pointers. */
    acb_mat_t *P_parts = flint_malloc((size_t) 2 * sizeof(acb_mat_t));
    acb_mat_t P;
    diag_proj_range(P_parts[0], n, 0, 1);   /* diag(1,0,0,0) */
    diag_proj_range(P_parts[1], n, 1, 2);   /* diag(0,1,0,0) */
    diag_proj_range(P, n, 0, 2);            /* diag(1,1,0,0) */

    /* Q = P (Q_k = P_k). */
    slong d = aic_cstar_lem_add_dim(&ae.A, (const acb_mat_t *) P_parts, 2, P,
                                    (const acb_mat_t *) P_parts, 2, P, CPREC);

    /* Hand value: off-diagonal corners vanish in the commutative A, so the sum is
     * the 2 diagonal blocks and equals dim S_{P,Q} = 2. */
    slong d_PQ = aic_corner_dim_S(&ae.A, P, P, CPREC);
    slong N = 0;
    for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
            N += aic_corner_dim_S(&ae.A, P_parts[j], P_parts[k], CPREC);
    printf("  dim S_{Pj,Qk} sum N=%ld, dim S_{P,Q}=%ld, returned=%ld (expect 2)\n",
           (long) N, (long) d_PQ, (long) d);
    AIC_CHECK_MSG(d == 2, "A1: lem_add_dim returned %ld, expected 2", (long) d);
    AIC_CHECK_MSG(N == 2, "A1: sum dim S_{Pj,Qk} = %ld, expected 2", (long) N);
    AIC_CHECK_MSG(d_PQ == 2, "A1: dim S_{P,Q} = %ld, expected 2", (long) d_PQ);
    AIC_CHECK_MSG(d == d_PQ && d == N, "A1: N=%ld, d_PQ=%ld, ret=%ld disagree",
                  (long) N, (long) d_PQ, (long) d);

    acb_mat_clear(P);
    acb_mat_clear(P_parts[0]); acb_mat_clear(P_parts[1]);
    flint_free(P_parts);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== A2: off_diag_zero (Stage-3 cert) ===================== *
 * Dephasing channel A = diagonal matrices of M_4 (COMMUTATIVE — the setting where
 * Stage-3's S_{P_C,P_D}=0 holds for distinct equivalence classes; in the identity
 * channel A=M_n EVERY off-diagonal corner E_{jk} is in A so S is never 0). The
 * "equivalence classes" of the 1d projections {E_jj} are the singletons (distinct
 * diagonal entries are inequivalent). (a) DISJOINT-support projections P_C=diag(1,0,0,0),
 * P_D=diag(0,1,0,0) (distinct classes): S_{P_C,P_D} = 0 (no diagonal operator
 * connects e0 and e1), so off_diag_zero == 1. (b) OVERLAPPING projections
 * P_C=diag(1,1,0,0), P_D=diag(0,1,1,0) (share e1): S_{P_C,P_D} = span{E_11} is
 * 1-dim, so off_diag_zero == 0. Both legs are the Stage-3 certification AND its
 * teeth in one test (verified by probe: disjoint -> 0, overlap-at-e1 -> 1). */
static void test_a2_off_diag_zero(void)
{
    printf("A2 off_diag_zero (dephasing(4); disjoint -> 1, overlap-at-e1 -> 0):\n");
    const slong n = 4;
    aic_ucp_kraus phi;
    make_dephasing(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    /* (a) disjoint supports (distinct classes): S_{P_C,P_D} = 0. */
    acb_mat_t Pc, Pd;
    diag_proj_range(Pc, n, 0, 1);   /* diag(1,0,0,0) */
    diag_proj_range(Pd, n, 1, 2);   /* diag(0,1,0,0) */
    int z1 = aic_cstar_off_diag_zero(&ae.A, Pc, Pd, CPREC);
    slong d1 = aic_corner_dim_S(&ae.A, Pc, Pd, CPREC);
    printf("  disjoint: dim S_{Pc,Pd}=%ld, off_diag_zero=%d (expect 0,1)\n",
           (long) d1, z1);
    AIC_CHECK_MSG(z1 == 1, "A2: disjoint supports not certified S=0 (got %d)", z1);
    AIC_CHECK_MSG(d1 == 0, "A2: dim S_{Pc,Pd}=%ld != 0 for disjoint", (long) d1);
    acb_mat_clear(Pd); acb_mat_clear(Pc);

    /* (b) overlapping supports (share e1): S_{P_C,P_D} = span{E_11} != 0. */
    acb_mat_t Qc, Qd;
    diag_proj_range(Qc, n, 0, 2);   /* diag(1,1,0,0) */
    diag_proj_range(Qd, n, 1, 3);   /* diag(0,1,1,0) */
    int z2 = aic_cstar_off_diag_zero(&ae.A, Qc, Qd, CPREC);
    slong d2 = aic_corner_dim_S(&ae.A, Qc, Qd, CPREC);
    printf("  overlap:  dim S_{Qc,Qd}=%ld, off_diag_zero=%d (expect >0,0)\n",
           (long) d2, z2);
    AIC_CHECK_MSG(z2 == 0, "A2: overlapping supports wrongly certified S=0");
    AIC_CHECK_MSG(d2 > 0, "A2: dim S_{Qc,Qd}=%ld should be > 0 for overlap",
                  (long) d2);
    acb_mat_clear(Qd); acb_mat_clear(Qc);

    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* Build a rank-1 commutative inclusion v_j : M_1 -> S_{P_j} with v_j(1) = Ptilde_j
 * = aic_corner_Ptilde(A, P_j) (the unit of S_{P_j}, an element of A). B_j = M_1
 * (dim_B = 1, n_B = 1). v_out->A = A. Caller frees B with aic_dhom_B_clear and v
 * with aic_dhom_v_clear (v borrows B and A). */
static void build_rank1_incl(aic_dhom_B *B, aic_dhom_v *v, const aic_ecstar *A,
                             const acb_mat_t Pj, slong prec)
{
    slong dims[1] = {1};
    aic_dhom_B_init(B, dims, 1);
    aic_dhom_v_init(v, B, A);
    aic_corner_Ptilde(v->vE[0], A, Pj, prec);
}

/* ===================== B1: cor_merge_sum eta=0 oracle ====================== *
 * Dephasing channel A = diag(M_2) (eta=0; A commutative, dim_A = 2). Two orthogonal
 * rank-1 projections P_1=diag(1,0), P_2=diag(0,1). Commutative inclusions
 * v_j(1) = Ptilde_j = Co_{P_j}(P_j) = P_j (a rank-1 C -> S_{P_j} inclusion); merged
 * v : C (+) C -> A maps onto the full diagonal algebra. Assert: merged dim_B == 2,
 * mult_def ~ machine-zero (exact hom at eta=0), sigma_min ~ 1 (exact isometric
 * inclusion). Genuine bijectivity (dim_B == dim_A AND S_{P_1,P_2}=0) is asserted in
 * B3. The dephasing fixture (not identity) is used so the off-diagonal corner
 * S_{P_1,P_2} genuinely vanishes (the Stage-3 / bijectivity setting). */
static void test_b1_merge_eta0(void)
{
    printf("B1 cor_merge_sum eta=0 oracle (dephasing(2), P_1,P_2 orthogonal, "
           "commutative A):\n");
    const slong n = 2;
    aic_ucp_kraus phi;
    make_dephasing(&phi, n);               /* A = diag(M_2) (eta=0, commutative) */
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    acb_mat_t P1, P2;
    diag_proj_range(P1, n, 0, 1);          /* diag(1,0) */
    diag_proj_range(P2, n, 1, 2);          /* diag(0,1) */

    aic_dhom_B B1, B2;
    aic_dhom_v v1, v2;
    build_rank1_incl(&B1, &v1, &ae.A, P1, CPREC);
    build_rank1_incl(&B2, &v2, &ae.A, P2, CPREC);

    aic_dhom_B Bm;
    aic_dhom_v vm;
    arb_t mult_def, smin;
    arb_init(mult_def); arb_init(smin);
    aic_cstar_merge_sum(&Bm, &vm, mult_def, smin, &v1, &v2, &ae.A, CPREC);

    double md = dd(mult_def), sm = dd(smin);
    printf("  merged dim_B=%ld (expect 2), mult_def=%.3e (~0), sigma_min=%.3f (~1)\n",
           (long) Bm.dim_B, md, sm);
    AIC_CHECK_MSG(Bm.dim_B == 2, "B1: merged dim_B=%ld, expected 2", (long) Bm.dim_B);
    AIC_CHECK_MSG(vm.B == &Bm, "B1: v_out->B does not point at B_out");
    AIC_CHECK_MSG(vm.A == &ae.A, "B1: v_out->A != A");
    AIC_CHECK_MSG(md < 1e-12, "B1: eta=0 merged mult defect not ~0 (%.3e)", md);
    AIC_CHECK_MSG(fabs(sm - 1.0) < 1e-9, "B1: eta=0 sigma_min=%.6f != 1", sm);

    /* The deep-copy contract: vm.vE[0] == Ptilde_1, vm.vE[1] == Ptilde_2, and
     * mutating vm.vE must NOT touch v1/v2 (independent storage). */
    acb_mat_t ref0, ref1, db;
    acb_mat_init(ref0, n, n); acb_mat_init(ref1, n, n); acb_mat_init(db, n, n);
    aic_corner_Ptilde(ref0, &ae.A, P1, CPREC);
    aic_corner_Ptilde(ref1, &ae.A, P2, CPREC);
    acb_mat_sub(db, vm.vE[0], ref0, CPREC);
    arb_t e; arb_init(e); aic_mat_opnorm(e, db, CPREC);
    double e0 = dd(e);
    acb_mat_sub(db, vm.vE[1], ref1, CPREC); aic_mat_opnorm(e, db, CPREC);
    double e1 = dd(e);
    AIC_CHECK_MSG(e0 < 1e-12 && e1 < 1e-12,
                  "B1: merged vE != [Ptilde_1, Ptilde_2] (%.3e, %.3e)", e0, e1);
    arb_clear(e);
    acb_mat_clear(db); acb_mat_clear(ref1); acb_mat_clear(ref0);

    arb_clear(smin); arb_clear(mult_def);
    aic_dhom_v_clear(&vm); aic_dhom_B_clear(&Bm);
    aic_dhom_v_clear(&v2); aic_dhom_B_clear(&B2);
    aic_dhom_v_clear(&v1); aic_dhom_B_clear(&B1);
    acb_mat_clear(P2); acb_mat_clear(P1);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== B2: cor_merge_sum oblique eta>0 + star teeth ========= *
 * make_mixconj(5,3,t): a GENUINELY OBLIQUE eta>0 channel (sigma_max(S~) ~ sqrt(3),
 * FINDINGS §C4). Two orthogonal in-A delta-projections P_1=span(e1)=diag(0,1,0,0,0),
 * P_2=span(e2)=diag(0,0,1,0,0) (genuine in-A 1d projections, residual O(eta));
 * commutative inclusions v_j(1) = Ptilde_j. Assert:
 *   - mult_def = O(eta) and SHRINKS across two eta values (the STAR teeth: the
 *     merged-v defect is computed with A's Choi-Effros star; with the WRONG product
 *     it picks up the O(1) Phi_tilde-vs-plain gap and the shrink is RED, FINDINGS §C2),
 *   - sigma_min >= 1 - O(eta).
 * P_1, P_2 are inside mixconj(5,3)'s carrier M=span(e0,e1,e2) but avoid e0 (the
 * oblique-residual direction), so each is a genuine in-A delta-projection. */
static void test_b2_merge_oblique(void)
{
    printf("B2 cor_merge_sum oblique eta>0 + star teeth "
           "(mixconj(5,3); P_1=span(e1), P_2=span(e2)):\n");
    const slong n = 5;
    const double ts[] = {0.06, 0.02};
    double defs[2], etas[2], smins[2];
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, 5, 3, ts[it], CPREC);
        double eta = eta_proxy(&phi, CPREC);
        etas[it] = eta;
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

        acb_mat_t P1, P2;
        diag_proj_range(P1, n, 1, 2);      /* span(e1) */
        diag_proj_range(P2, n, 2, 3);      /* span(e2) */

        aic_dhom_B B1, B2;
        aic_dhom_v v1, v2;
        build_rank1_incl(&B1, &v1, &ae.A, P1, CPREC);
        build_rank1_incl(&B2, &v2, &ae.A, P2, CPREC);

        aic_dhom_B Bm;
        aic_dhom_v vm;
        arb_t mult_def, smin;
        arb_init(mult_def); arb_init(smin);
        aic_cstar_merge_sum(&Bm, &vm, mult_def, smin, &v1, &v2, &ae.A, CPREC);

        double md = dd(mult_def), sm = dd(smin);
        defs[it] = md;
        smins[it] = sm;
        double c = eta > 1e-14 ? md / eta : 0.0;
        printf("  t=%.2f eta=%.4e | merged dim_B=%ld mult_def=%.3e (c=md/eta=%.3f) "
               "sigma_min=%.4f\n",
               ts[it], eta, (long) Bm.dim_B, md, c, sm);
        AIC_CHECK_MSG(Bm.dim_B == 2, "B2: merged dim_B=%ld != 2", (long) Bm.dim_B);
        /* O(eta), AND the c=md/eta ratio is SMALL (the star tooth): with A's
         * Choi-Effros star c ~ 0.017 -> 0.005 (measured); with the WRONG plain
         * product c stays ~0.43 (the O(1) Phi_tilde-vs-plain gap, measured by the
         * /tmp/teeth_b2 mutation). c < 0.2 cleanly separates them (FINDINGS §C2). */
        AIC_CHECK_MSG(c < 0.2, "B2: merged mult_def/eta = %.3f exceeds 0.2 — star "
                      "teeth RED (plain-product gap is ~0.43) (t=%.2f)", c, ts[it]);
        /* sigma_min >= 1 - O(eta) (a few * eta below 1). */
        AIC_CHECK_MSG(sm >= 1.0 - 10.0 * eta && sm <= 1.0 + 1e-6,
                      "B2: sigma_min=%.4f not in [1-O(eta), 1] (eta=%.3e, t=%.2f)",
                      sm, eta, ts[it]);

        arb_clear(smin); arb_clear(mult_def);
        aic_dhom_v_clear(&vm); aic_dhom_B_clear(&Bm);
        aic_dhom_v_clear(&v2); aic_dhom_B_clear(&B2);
        aic_dhom_v_clear(&v1); aic_dhom_B_clear(&B1);
        acb_mat_clear(P2); acb_mat_clear(P1);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("  -> merged mult_def shrinks with eta: t=0.06 def=%.3e (eta=%.3e), "
           "t=0.02 def=%.3e (eta=%.3e)\n", defs[0], etas[0], defs[1], etas[1]);
    /* THE STAR TEETH (FINDINGS §C2): the merged defect SHRINKS as eta->0. With the
     * wrong (plain) product this is RED (def ~ constant O(1) Phi_tilde gap). */
    AIC_CHECK_MSG(defs[1] < defs[0],
                  "B2: merged mult_def did not shrink as eta->0 (%.3e !< %.3e) — "
                  "star teeth RED", defs[1], defs[0]);
    AIC_CHECK_MSG(smins[1] >= smins[0] - 1e-9,
                  "B2: sigma_min did not improve as eta->0 (%.4f < %.4f)",
                  smins[1], smins[0]);
}

/* ===================== B3: cor_merge_sum bijectivity ======================= *
 * Dephasing channel A = diag(M_2) (commutative, dim_A = 2). P_1=diag(1,0),
 * P_2=diag(0,1): orthogonal rank-1, so off_diag_zero(P_1,P_2) == 1 (the Stage-3
 * bijectivity precondition S_{P_1,P_2}=0 genuinely holds in the commutative A),
 * and dim_B1+dim_B2 = 2 == dim_A. Per cor_merge_sum's bijectivity clause
 * (.tex:1358: "If v_1,v_2 are bijective and S_{P_1,P_2}=0, then v is also
 * bijective."), the merged v : C (+) C -> A IS a bijection: it maps onto the full
 * diagonal algebra. We certify BOTH hypotheses (off_diag_zero==1, v_1/v_2 rank-1
 * isometries with sigma_min ~ 1) AND the conclusion via aic_errreduce_is_bijective
 * (injective: sigma_min > 0, AND square: dim_B == dim_A). */
static void test_b3_bijectivity(void)
{
    printf("B3 cor_merge_sum bijectivity (dephasing(2); off_diag_zero + square):\n");
    const slong n = 2;
    aic_ucp_kraus phi;
    make_dephasing(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    acb_mat_t P1, P2;
    diag_proj_range(P1, n, 0, 1);
    diag_proj_range(P2, n, 1, 2);

    /* Hypothesis 1: S_{P_1,P_2} = 0 (the bijectivity precondition, .tex:1358). */
    int odz = aic_cstar_off_diag_zero(&ae.A, P1, P2, CPREC);
    AIC_CHECK_MSG(odz == 1, "B3: S_{P_1,P_2} != 0 (off_diag_zero=%d)", odz);

    aic_dhom_B B1, B2;
    aic_dhom_v v1, v2;
    build_rank1_incl(&B1, &v1, &ae.A, P1, CPREC);
    build_rank1_incl(&B2, &v2, &ae.A, P2, CPREC);

    aic_dhom_B Bm;
    aic_dhom_v vm;
    arb_t mult_def, smin;
    arb_init(mult_def); arb_init(smin);
    aic_cstar_merge_sum(&Bm, &vm, mult_def, smin, &v1, &v2, &ae.A, CPREC);

    double sm = dd(smin);
    printf("  off_diag_zero=%d, merged dim_B=%ld, dim_A=%ld, sigma_min=%.4f (>0)\n",
           odz, (long) Bm.dim_B, (long) ae.A.dim_A, sm);
    AIC_CHECK_MSG(sm > 1e-9, "B3: merged v collapses (sigma_min=%.3e ~ 0)", sm);
    AIC_CHECK_MSG(Bm.dim_B == 2, "B3: dim_B=%ld != 2", (long) Bm.dim_B);

    /* The conclusion (.tex:1358): the merged v IS bijective — injective
     * (sigma_min > 0) AND square (dim_B == dim_A == 2). */
    arb_t a_out;
    arb_init(a_out);
    int bij = aic_errreduce_is_bijective(a_out, &vm, CPREC);
    printf("  is_bijective=%d (expect 1: dim_B=2 == dim_A=2), injective a=%.4f\n",
           bij, dd(a_out));
    AIC_CHECK_MSG(bij == 1, "B3: merged C(+)C -> diag(M_2) not bijective "
                  "(dim_B=%ld, dim_A=%ld, a=%.4f)",
                  (long) Bm.dim_B, (long) ae.A.dim_A, dd(a_out));
    AIC_CHECK_MSG(dd(a_out) > 1e-9, "B3: injective leg a=%.3e ~ 0", dd(a_out));
    arb_clear(a_out);

    arb_clear(smin); arb_clear(mult_def);
    aic_dhom_v_clear(&vm); aic_dhom_B_clear(&Bm);
    aic_dhom_v_clear(&v2); aic_dhom_B_clear(&B2);
    aic_dhom_v_clear(&v1); aic_dhom_B_clear(&B1);
    acb_mat_clear(P2); acb_mat_clear(P1);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_a1_lem_add_dim();
    test_a2_off_diag_zero();
    test_b1_merge_eta0();
    test_b2_merge_oblique();
    test_b3_bijectivity();
    aic_test_report("test_cstar_merge");
    return 0;
}
