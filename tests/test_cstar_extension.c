/* test_cstar_extension.c — cross-checks for cstar_build Increment 4 (bead aic-097):
 * lem_extension (.tex:1378-1411), the inductive Stage-2 growth step v : M_n -> S_P
 * |-> v_+ : M_{n+1} -> A. See include/aic_cstar.h (Increment 4 banner) and
 * src/aic_cstar_extension.c. lem_extension consumes I3's lem_merging (two_block=0).
 *
 * NO NEW double-vs-arb path: the new code is STRUCTURAL GLUE (the Ha-map loop +
 * lem_approx + the SVD-based U_1 extraction + the gamma assembly) over already-tested
 * numerical routines (aic_corner_ha, aic_dhom_approx, aic_latd_svd, aic_cstar_lem
 * _merging). The cross-checks are:
 *   S0 — §10.C smoke test: aic_corner_ha reduces to exact compression on genuine M_2
 *        (Ha^Q_{P,P}(Ptilde) == [1] on the 1-d S_{P,Q}); run FIRST so a corner-vs-
 *        genuine incompatibility surfaces before the full pipeline.
 *   T1 — eta=0 base case n=1: A = genuine M_2, P=diag(1,0), Q=diag(0,1), v: M_1->S_P
 *        the inclusion. v_+ : M_2 -> M_2 is the IDENTITY (mult_def~0, sigma_min~1,
 *        v_+(E_lm)==E_lm). The cleanest ground truth (Rule 6 eta=0 oracle).
 *   T2 — eta=0 inductive n=2 (M_2 -> M_3): A = M_3, P=diag(1,1,0), Q=diag(0,0,1),
 *        v: M_2->S_P the natural upper-left inclusion. v_+ : M_3 -> M_3 the identity.
 *   T4 — eta=0 NON-SYMMETRIC U_1 directional tooth (hostile-review Finding #1):
 *        A = M_3, same P,Q,n=2, but v = conjugation by a 2x2 REAL ROTATION
 *        (theta=0.6) embedded in M_3, so the CORRECT U_1 ~= a non-symmetric unitary.
 *        v_+ is U-conjugated (NOT the identity); the phase-invariant teeth mult_def~0
 *        and sigma_min~1 break under BOTH replacing extract_U1 with identity AND
 *        transposing U1coord (the SVD-convention/directional teeth S0/T1/T2/T3 lack).
 *   T3 — OBLIQUE eta>0. Two legs honoring FINDINGS §C8 (the c=mult_def/eta MAGNITUDE
 *        tooth) and the spec §10.B ||P+Q-unit|| <= O(eta) hypothesis:
 *        (a) eta=0 M_2-WRAPPER-as-parent inductive: proves lem_extension works with a
 *            COMPRESSED S_P as A_parent and P+Q = the wrapper's unit EXACTLY (the
 *            ||P+Q-unit||=0 complementarity §10.B demands; star=plain at eta=0).
 *        (b) genuinely-OBLIQUE end-to-end on raw mixconj(5,3) (n=1, P=span(e1),
 *            Q=span(e2)): the §C8 star tooth c_star < 0.2 << c_plain (~0.5). This is
 *            a CORNER-LOCAL extension (P+Q = span(e1,e2) != 1_5, PRINTED honestly);
 *            the mult_def multiplicativity tooth is the valid non-vacuous content.
 *        The genuinely-oblique S_P-WRAPPER end-to-end is DEFERRED to aic-qgs: the
 *        oblique wrapper's corner Co hits the §C5 Gram-Hermiticity false-fail in the
 *        sgn-basin opnorm (FINDINGS §C10; measured ||(2M-1)^2-I|| ~1e-3, in-basin).
 *
 * MUTATION-PROVEN teeth (Rule 7; each RED observed BY HAND, then restored — recorded
 * in the increment report, NOT left failing in the suite):
 *   (a) feed h_{11} the wrong pair (Q,P,Q) instead of (P,P,Q) -> 1xn into nxn slot,
 *       aic_corner_ha dimension assert RED (+ a dimension-preserving index mutation);
 *   (b) star -> plain in the defect sweep -> c crosses 0.2 (T3b);
 *   (c) gamma_12 -> 0 -> sigma_min -> 0, mult_def -> O(1);
 *   (d) U_1 unnormalized (u0 /= 2) -> merging3 / mult_def blow up.
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
    (void) make_dephasing;
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op (mirrors test_cstar_merging.c:54): take the ball
 * MIDPOINTS then aic_mat_opnorm — the §C5 Gram false-fail workaround at the test
 * layer (midpoints drop the matmul-accumulated radii that trip the Hermiticity check). */
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

/* ||A - B||_op via the §C5 midpoint route (mirrors test_cstar_merging.c:87): subtract,
 * take ball midpoints, aic_mat_opnorm — drops the accumulated radii that false-fail the
 * Gram Hermiticity check on near-zero-off-diagonal star-defect differences. */
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

/* Build v : M_n -> S_P the NATURAL inclusion on a GENUINE matrix algebra A = M_N
 * (N >= base+n), with P projecting onto the first n coords. v(E_lm) = E_lm (the
 * single-1 N x N matrix at (l,m), l,m < n). Caller frees with the dhom clears. */
static void build_natural_incl(aic_dhom_B *B, aic_dhom_v *v, const aic_ecstar *A,
                               slong n)
{
    slong N = A->n;
    slong dims[1] = {n};
    aic_dhom_B_init(B, dims, 1);
    aic_dhom_v_init(v, B, A);
    for (slong l = 0; l < n; l++)
        for (slong m = 0; m < n; m++) {
            acb_mat_zero(v->vE[l * n + m]);
            acb_set_si(acb_mat_entry(v->vE[l * n + m], l, m), 1);
        }
    (void) N;
}

/* Assert v_+ == the identity inclusion M_N -> M_N (v_+(E_lm) == E_lm for all l,m
 * in [0,N)); return the worst ||v_+(E_lm) - E_lm||_op.
 * PHASE-FRAGILE (Finding #6): this literal v_+(E_lm)==E_lm equality holds ONLY
 * because the real-spectrum eta=0 oracles (S0/T1/T2) give the extract_U1 SVD phase
 * phi=0 (u0 = e0); a phased left singular vector would false-RED on CORRECT code
 * (gamma_12/gamma_21 cancel the phase, so the iso property is phase-invariant). Use
 * it ONLY on those real-spectrum fixtures; the phase/permutation-INVARIANT teeth are
 * mult_def (~0) and sigma_min (~1), which T4 (a NON-symmetric, hence non-identity
 * U_1) uses instead of this helper. */
static double check_vplus_identity(const aic_dhom_v *vp, slong N, slong prec)
{
    double worst = 0.0;
    for (slong l = 0; l < N; l++)
        for (slong m = 0; m < N; m++) {
            acb_mat_t Elm;
            acb_mat_init(Elm, N, N);
            acb_mat_zero(Elm);
            acb_set_si(acb_mat_entry(Elm, l, m), 1);
            double e = opnorm_diff_d(vp->vE[l * N + m], Elm, prec);
            if (e > worst) worst = e;
            acb_mat_clear(Elm);
        }
    return worst;
}

/* ===================== S0: §10.C aic_corner_ha smoke test ================== *
 * On genuine M_2, P=diag(1,0), Q=diag(0,1): dim S_{P,Q} = dim S_P = dim S_Q = 1,
 * and h_{11}(Ptilde_P) = Ha^Q_{P,P}(Ptilde_P) == [1] (the identity on the 1-d
 * S_{P,Q}). Surfaces a corner-vs-genuine incompatibility BEFORE the full pipeline. */
static void test_s0_ha_smoke(void)
{
    printf("S0 aic_corner_ha smoke on genuine M_2 (P=diag(1,0), Q=diag(0,1)):\n");
    aic_ecstar A;
    aic_cstar_matrix_algebra(&A, 2, CPREC);
    acb_mat_t P, Q;
    diag_proj_range(P, 2, 0, 1);
    diag_proj_range(Q, 2, 1, 2);

    slong dPQ = aic_corner_dim_S(&A, P, Q, CPREC);
    slong dQ  = aic_corner_dim_S(&A, Q, Q, CPREC);
    printf("  dim S_{P,Q}=%ld, dim S_Q=%ld (both expect 1)\n", (long) dPQ, (long) dQ);
    AIC_CHECK_MSG(dPQ == 1, "S0: dim S_{P,Q}=%ld != 1", (long) dPQ);
    AIC_CHECK_MSG(dQ == 1, "S0: dim S_Q=%ld != 1", (long) dQ);

    acb_mat_t Ha, PtP;
    acb_mat_init(Ha, dPQ, dPQ);
    acb_mat_init(PtP, 2, 2);
    aic_corner_Ptilde(PtP, &A, P, CPREC);
    aic_corner_ha(Ha, &A, PtP, P, P, Q, CPREC);
    acb_t one;
    acb_init(one);
    acb_set_si(one, 1);
    double dha = 0.0;
    {
        acb_t diff;
        acb_init(diff);
        acb_sub(diff, acb_mat_entry(Ha, 0, 0), one, CPREC);
        arb_t mag;
        arb_init(mag);
        acb_abs(mag, diff, CPREC);
        dha = dd(mag);
        arb_clear(mag);
        acb_clear(diff);
    }
    printf("  |Ha^Q_{P,P}(Ptilde_P) - 1| = %.3e (expect ~0; the identity on S_{P,Q})\n",
           dha);
    AIC_CHECK_MSG(dha < 1e-12, "S0: Ha^Q_{P,P}(Ptilde) != [1] (%.3e)", dha);

    acb_clear(one);
    acb_mat_clear(PtP); acb_mat_clear(Ha);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_ecstar_clear(&A);
}

/* ===================== T1: eta=0 base case n=1 (M_1 -> M_2) ================= */
static void test_t1_eta0_base(void)
{
    printf("T1 eta=0 base n=1 (A=M_2, P=diag(1,0), Q=diag(0,1), v(1)=Ptilde_P):\n");
    const slong N = 2, n = 1;
    aic_ecstar A;
    aic_cstar_matrix_algebra(&A, N, CPREC);
    acb_mat_t P, Q;
    diag_proj_range(P, N, 0, 1);
    diag_proj_range(Q, N, 1, 2);

    /* v: M_1 -> S_P, v(1) = Ptilde_P (the unit of the 1-d S_P = span(E_00)). */
    aic_dhom_B B;
    aic_dhom_v v;
    slong dims[1] = {n};
    aic_dhom_B_init(&B, dims, 1);
    aic_dhom_v_init(&v, &B, &A);
    aic_corner_Ptilde(v.vE[0], &A, P, CPREC);

    aic_dhom_B Bp;
    aic_dhom_v vp;
    arb_t md, sm;
    arb_init(md); arb_init(sm);
    aic_cstar_lem_extension(&Bp, &vp, md, sm, &v, P, Q, &A, CPREC);

    double mdv = dd(md), smv = dd(sm);
    double worst = check_vplus_identity(&vp, N, CPREC);
    printf("  v_+ dim_B=%ld (expect 4) | mult_def=%.3e (~0) sigma_min=%.4f (~1) | "
           "max||v_+(E_lm)-E_lm||=%.3e\n", (long) Bp.dim_B, mdv, smv, worst);
    AIC_CHECK_MSG(Bp.dim_B == 4, "T1: v_+ dim_B=%ld != 4", (long) Bp.dim_B);
    AIC_CHECK_MSG(Bp.num_blocks == 1 && Bp.d[0] == 2,
                  "T1: v_+ B not a single M_2 block");
    AIC_CHECK_MSG(mdv < 1e-12, "T1: eta=0 mult_def not ~0 (%.3e)", mdv);
    AIC_CHECK_MSG(fabs(smv - 1.0) < 1e-9, "T1: eta=0 sigma_min=%.6f != 1", smv);
    AIC_CHECK_MSG(worst < 1e-11, "T1: v_+ != identity inclusion (%.3e)", worst);

    arb_clear(sm); arb_clear(md);
    aic_dhom_v_clear(&vp); aic_dhom_B_clear(&Bp);
    aic_dhom_v_clear(&v); aic_dhom_B_clear(&B);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_ecstar_clear(&A);
}

/* ===================== T2: eta=0 inductive n=2 (M_2 -> M_3) ================= */
static void test_t2_eta0_inductive(void)
{
    printf("T2 eta=0 inductive n=2 (A=M_3, P=diag(1,1,0), Q=diag(0,0,1), v natural "
           "M_2 incl):\n");
    const slong N = 3, n = 2;
    aic_ecstar A;
    aic_cstar_matrix_algebra(&A, N, CPREC);
    acb_mat_t P, Q;
    diag_proj_range(P, N, 0, 2);
    diag_proj_range(Q, N, 2, 3);

    aic_dhom_B B;
    aic_dhom_v v;
    build_natural_incl(&B, &v, &A, n);

    aic_dhom_B Bp;
    aic_dhom_v vp;
    arb_t md, sm;
    arb_init(md); arb_init(sm);
    aic_cstar_lem_extension(&Bp, &vp, md, sm, &v, P, Q, &A, CPREC);

    double mdv = dd(md), smv = dd(sm);
    double worst = check_vplus_identity(&vp, N, CPREC);
    printf("  v_+ dim_B=%ld (expect 9) | mult_def=%.3e (~0) sigma_min=%.4f (~1) | "
           "max||v_+(E_lm)-E_lm||=%.3e\n", (long) Bp.dim_B, mdv, smv, worst);
    AIC_CHECK_MSG(Bp.dim_B == 9, "T2: v_+ dim_B=%ld != 9", (long) Bp.dim_B);
    AIC_CHECK_MSG(Bp.num_blocks == 1 && Bp.d[0] == 3,
                  "T2: v_+ B not a single M_3 block");
    AIC_CHECK_MSG(mdv < 1e-12, "T2: eta=0 mult_def not ~0 (%.3e)", mdv);
    AIC_CHECK_MSG(fabs(smv - 1.0) < 1e-9, "T2: eta=0 sigma_min=%.6f != 1", smv);
    AIC_CHECK_MSG(worst < 1e-11, "T2: v_+ != identity inclusion (%.3e)", worst);

    arb_clear(sm); arb_clear(md);
    aic_dhom_v_clear(&vp); aic_dhom_B_clear(&Bp);
    aic_dhom_v_clear(&v); aic_dhom_B_clear(&B);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_ecstar_clear(&A);
}

/* ===================== T4: eta=0 NON-SYMMETRIC U_1 directional tooth ======= *
 * THE extract_U1 directional/SVD-convention tooth (hostile-review Finding #1). All
 * of S0/T1/T2/T3a have the correct U_1 ~= identity (orthogonal projectors on real-
 * spectrum eta=0 oracles, or a 1x1 phase), so neither replacing the whole U1coord
 * fill with acb_mat_one (identity, ignoring mu11) NOR transposing U1coord changes
 * any check — the SVD left-vector reading and the .tex:1404 formula are never
 * distinguished from "return identity." T4 forces the CORRECT U_1 to be a genuine
 * NON-SYMMETRIC unitary so BOTH mutations break it.
 *
 * Setup: A = genuine M_3 (.tex:1371 matrix units), P=diag(1,1,0), Q=diag(0,0,1),
 * n=2. U = a 2x2 REAL ROTATION by theta=0.6 rad (U^T = U(-theta) != U, the
 * NON-symmetry that gives the transpose tooth), embedded in the 3x3 identity acting
 * on coords {0,1} and fixing coord 2 (so U_full is block-diagonal). v : M_2 -> S_P
 * is the inner automorphism v(E_lm) = U_full E_lm U_full^dag (a valid *-isomorphism
 * of the upper-left 2x2 = S_P block; E_lm in the 2-block, l,m in {0,1}, stays in
 * the block under conjugation by block-diagonal U_full). The CORRECT U_1 ~= U.
 *
 * v_+ is the U-conjugated isomorphism, NOT the identity, so check_vplus_identity
 * does NOT apply (spec §7.1). The PHASE/permutation-invariant teeth are the
 * isomorphism invariants: mult_def ~ 0 (multiplicativity through A's star, which
 * the wrong U_1 violates on the cross-block products) and sigma_min ~ 1. A wrong
 * U_1 (identity or U^T) makes gamma_12/gamma_21 inconsistent with gamma_11 = v's
 * U-conjugation, so v_+(E_{l,2}) * v_+(E_{2,m}) != v_+(E_{l,m}) and mult_def -> O(1).
 * MUTATION-PROVEN by hand (Finding #1): identity-substitution and U1coord-transpose
 * each drive T4 mult_def RED while S0/T1/T2/T3a/T3b stay GREEN; numbers in the
 * increment report. */
static void test_t4_nonsym_u1(void)
{
    printf("T4 eta=0 NON-SYMMETRIC U_1 (A=M_3, P=diag(1,1,0), Q=diag(0,0,1), n=2; "
           "v = conj by 2x2 rotation theta=0.6):\n");
    const slong N = 3, n = 2;
    aic_ecstar A;
    aic_cstar_matrix_algebra(&A, N, CPREC);
    acb_mat_t P, Q;
    diag_proj_range(P, N, 0, 2);
    diag_proj_range(Q, N, 2, 3);

    /* U_full = 3x3 identity with a 2x2 rotation by theta on coords {0,1}. */
    const double theta = 0.6;
    acb_mat_t Uf;
    acb_mat_init(Uf, N, N);
    acb_mat_one(Uf);
    acb_set_d(acb_mat_entry(Uf, 0, 0), cos(theta));
    acb_set_d(acb_mat_entry(Uf, 0, 1), -sin(theta));
    acb_set_d(acb_mat_entry(Uf, 1, 0), sin(theta));
    acb_set_d(acb_mat_entry(Uf, 1, 1), cos(theta));

    /* v : M_2 -> S_P, v(E_lm) = U_full * E_lm * U_full^dag (l,m in {0,1}, the 2-block
     * embedded in M_3). Each image is N x N and lands in the upper-left 2x2 = S_P. */
    aic_dhom_B B;
    aic_dhom_v v;
    slong dims[1] = {n};
    aic_dhom_B_init(&B, dims, 1);
    aic_dhom_v_init(&v, &B, &A);
    acb_mat_t Elm, Ufd, tmp;
    acb_mat_init(Elm, N, N);
    acb_mat_init(Ufd, N, N);
    acb_mat_init(tmp, N, N);
    acb_mat_conjugate_transpose(Ufd, Uf);
    for (slong l = 0; l < n; l++)
        for (slong m = 0; m < n; m++) {
            acb_mat_zero(Elm);
            acb_set_si(acb_mat_entry(Elm, l, m), 1);
            acb_mat_mul(tmp, Uf, Elm, CPREC);            /* U_full E_lm */
            acb_mat_mul(v.vE[l * n + m], tmp, Ufd, CPREC); /* ... U_full^dag */
        }
    acb_mat_clear(tmp); acb_mat_clear(Ufd); acb_mat_clear(Elm);

    aic_dhom_B Bp;
    aic_dhom_v vp;
    arb_t md, sm;
    arb_init(md); arb_init(sm);
    aic_cstar_lem_extension(&Bp, &vp, md, sm, &v, P, Q, &A, CPREC);

    double mdv = dd(md), smv = dd(sm);
    printf("  v_+ dim_B=%ld (expect 9) | mult_def=%.3e (~0; the directional tooth) "
           "sigma_min=%.4f (~1) [NOT identity: v_+ is U-conjugated, spec §7.1]\n",
           (long) Bp.dim_B, mdv, smv);
    AIC_CHECK_MSG(Bp.dim_B == 9, "T4: v_+ dim_B=%ld != 9", (long) Bp.dim_B);
    AIC_CHECK_MSG(Bp.num_blocks == 1 && Bp.d[0] == 3,
                  "T4: v_+ B not a single M_3 block");
    /* The phase/permutation-invariant ISOMORPHISM teeth (Finding #1): a wrong U_1
     * (identity-substitution OR U1coord-transpose) violates cross-block
     * multiplicativity -> mult_def O(1), and collapses sigma_min. */
    AIC_CHECK_MSG(mdv < 1e-10, "T4: non-sym U_1 mult_def=%.3e not ~0 (extract_U1 "
                  "directional/SVD-convention bug)", mdv);
    AIC_CHECK_MSG(smv > 0.95, "T4: non-sym U_1 sigma_min=%.4f too small", smv);

    arb_clear(sm); arb_clear(md);
    aic_dhom_v_clear(&vp); aic_dhom_B_clear(&Bp);
    aic_dhom_v_clear(&v); aic_dhom_B_clear(&B);
    acb_mat_clear(Uf);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_ecstar_clear(&A);
}

/* ===================== T3a: eta=0 M_2-WRAPPER-as-parent inductive ========== *
 * A_parent = the S_P=M_2 wrapper over span(e0,e1) of the identity channel M_3. The
 * wrapper is a GENUINE M_2 (dim_A=4) whose UNIT is Ptilde = diag(1,1,0). Inside it:
 * P=span(e0)=diag(1,0,0), Q=span(e1)=diag(0,1,0), so P+Q = Ptilde EXACTLY (the §10.B
 * ||P+Q-unit||=0 complementarity). v: M_1 -> S_P, v(1)=Ptilde_P. v_+ : M_2 -> wrapper
 * is the identity of the M_2 wrapper (mult_def~0, sigma_min~1). Proves lem_extension
 * works with a COMPRESSED S_P as A_parent (the Stage-2-of-th_main usage). */
static void test_t3a_wrapper_inductive(void)
{
    printf("T3a eta=0 M_2-wrapper-as-parent (identity M_3, wrapper over span(e0,e1); "
           "P=span(e0), Q=span(e1)):\n");
    aic_ucp_kraus phi;
    make_identity(&phi, 3);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    acb_mat_t Pg;
    diag_proj_range(Pg, 3, 0, 2);             /* span(e0,e1) -> M_2 wrapper */
    aic_cstar_subalgebra sg;
    aic_cstar_subalg_build(&sg, &ae.A, Pg, CPREC);
    AIC_CHECK_MSG(sg.A.dim_A == 4, "T3a: wrapper dim_A=%ld != 4", (long) sg.A.dim_A);

    acb_mat_t P, Q;
    diag_proj_range(P, 3, 0, 1);              /* span(e0) */
    diag_proj_range(Q, 3, 1, 2);              /* span(e1) */
    /* ||P + Q - unit(A_parent)|| where unit(A_parent) = sg.Ptilde (§10.B). */
    acb_mat_t PQ;
    acb_mat_init(PQ, 3, 3);
    acb_mat_add(PQ, P, Q, CPREC);
    double pq_unit = opnorm_diff_d(PQ, sg.Ptilde, CPREC);
    acb_mat_clear(PQ);
    printf("  ||P+Q-unit(wrapper)|| = %.3e (must be ~0; the complementarity §10.B)\n",
           pq_unit);
    AIC_CHECK_MSG(pq_unit < 1e-12, "T3a: ||P+Q-unit|| = %.3e not ~0 (vacuous input)",
                  pq_unit);

    aic_dhom_B B;
    aic_dhom_v v;
    slong dims[1] = {1};
    aic_dhom_B_init(&B, dims, 1);
    aic_dhom_v_init(&v, &B, &sg.A);
    aic_corner_Ptilde(v.vE[0], &sg.A, P, CPREC);

    aic_dhom_B Bp;
    aic_dhom_v vp;
    arb_t md, sm;
    arb_init(md); arb_init(sm);
    aic_cstar_lem_extension(&Bp, &vp, md, sm, &v, P, Q, &sg.A, CPREC);
    double mdv = dd(md), smv = dd(sm);
    printf("  v_+ dim_B=%ld (expect 4) | mult_def=%.3e (~0) sigma_min=%.4f (~1)\n",
           (long) Bp.dim_B, mdv, smv);
    AIC_CHECK_MSG(Bp.dim_B == 4, "T3a: v_+ dim_B=%ld != 4", (long) Bp.dim_B);
    AIC_CHECK_MSG(mdv < 1e-10, "T3a: wrapper-parent mult_def not ~0 (%.3e)", mdv);
    AIC_CHECK_MSG(smv > 0.95, "T3a: wrapper-parent sigma_min=%.4f too small", smv);

    arb_clear(sm); arb_clear(md);
    aic_dhom_v_clear(&vp); aic_dhom_B_clear(&Bp);
    aic_dhom_v_clear(&v); aic_dhom_B_clear(&B);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_cstar_subalg_clear(&sg); acb_mat_clear(Pg);
    aic_assoc_ecstar_clear(&ae); aic_ucp_kraus_clear(&phi);
}

/* plain-product star thunk (the §C8 mutation): out = XY (no Phi_tilde). */
static void plain_star(acb_mat_t out, const acb_mat_t XY, void *ctx, slong prec)
{
    (void) ctx;
    (void) prec;
    acb_mat_set(out, XY);
}

/* ===================== T3b: oblique eta>0 §C8 magnitude tooth ============== *
 * Genuinely OBLIQUE end-to-end lem_extension on the raw associated algebra of
 * make_mixconj(5,3,t) (sigma_max(S~) ~ sqrt(3); FINDINGS §C4). n=1, P=span(e1),
 * Q=span(e2): two EQUIVALENT 1-d in-A delta-projections (dim S_{P,Q}=1=n, dim S_Q=1).
 * This is a CORNER-LOCAL extension: P+Q = span(e1,e2) != 1_5 (the ambient unit),
 * PRINTED honestly — the genuinely-oblique S_P-wrapper end-to-end (where P+Q would be
 * the wrapper unit) is DEFERRED to aic-qgs (FINDINGS §C10: the oblique wrapper Co hits
 * the §C5 sgn-basin Gram false-fail). The valid non-vacuous content here is the
 * MULTIPLICATIVITY (merging1) star tooth: the FINDINGS §C8 c=mult_def/eta MAGNITUDE
 * gate. mutating A's star -> plain product pushes c from ~0.02 (star) to ~0.5 (plain),
 * a ~30x gap; the c < 0.2 bound is RED for plain. */
static void test_t3b_oblique_c8_tooth(void)
{
    printf("T3b oblique eta>0 end-to-end (mixconj(5,3); n=1, P=span(e1), Q=span(e2)); "
           "§C8 c=mult_def/eta star tooth:\n");
    const slong n = 5;
    const double ts[] = {0.06, 0.02};
    double cstar_s[2], cplain_s[2];
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, 5, 3, ts[it], CPREC);
        double eta = eta_proxy(&phi, CPREC);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

        acb_mat_t P, Q;
        diag_proj_range(P, n, 1, 2);          /* span(e1) */
        diag_proj_range(Q, n, 2, 3);          /* span(e2) */

        /* HONEST report: ||P+Q - 1_5|| (this is a corner-local extension). */
        acb_mat_t PQ, I5;
        acb_mat_init(PQ, n, n);
        acb_mat_init(I5, n, n);
        acb_mat_add(PQ, P, Q, CPREC);
        acb_mat_one(I5);
        double pq_unit = opnorm_diff_d(PQ, I5, CPREC);
        acb_mat_clear(I5); acb_mat_clear(PQ);

        aic_dhom_B B;
        aic_dhom_v v;
        slong dims[1] = {1};
        aic_dhom_B_init(&B, dims, 1);
        aic_dhom_v_init(&v, &B, &ae.A);
        aic_corner_Ptilde(v.vE[0], &ae.A, P, CPREC);

        aic_dhom_B Bp;
        aic_dhom_v vp;
        arb_t md, sm;
        arb_init(md); arb_init(sm);
        aic_cstar_lem_extension(&Bp, &vp, md, sm, &v, P, Q, &ae.A, CPREC);
        double md_star = dd(md), smv = dd(sm);
        double c_star = eta > 1e-14 ? md_star / eta : 0.0;

        /* §C8 mutation: swap A's star -> plain product, re-sweep the SAME v_+. */
        aic_ecstar_apply_fn save_phi = ae.A.star_phi;
        void *save_ctx = ae.A.star_ctx;
        ae.A.star_phi = plain_star;
        ae.A.star_ctx = NULL;
        arb_t mdp;
        arb_init(mdp);
        aic_dhom_defect_sweep(mdp, &vp, CPREC);
        double md_plain = dd(mdp);
        double c_plain = eta > 1e-14 ? md_plain / eta : 0.0;
        ae.A.star_phi = save_phi;
        ae.A.star_ctx = save_ctx;
        arb_clear(mdp);

        cstar_s[it] = c_star;
        cplain_s[it] = c_plain;
        printf("  t=%.2f eta=%.3e | v_+ dim_B=%ld ||P+Q-1_5||=%.3f (corner-local) | "
               "STAR mult_def=%.3e c=%.4f | PLAIN mult_def=%.3e c=%.4f | sigma_min=%.4f\n",
               ts[it], eta, (long) Bp.dim_B, pq_unit, md_star, c_star, md_plain,
               c_plain, smv);

        AIC_CHECK_MSG(Bp.dim_B == 4, "T3b: v_+ dim_B=%ld != 4", (long) Bp.dim_B);
        AIC_CHECK_MSG(smv > 0.5, "T3b: sigma_min=%.4f <= 0.5 (collapse)", smv);
        /* FINDINGS §C8: the c=mult_def/eta MAGNITUDE is the sharp star discriminant. */
        AIC_CHECK_MSG(c_star < 0.2, "T3b: STAR c=mult_def/eta=%.4f exceeds 0.2 "
                      "(t=%.2f)", c_star, ts[it]);
        AIC_CHECK_MSG(c_plain > 0.2, "T3b: PLAIN c=mult_def/eta=%.4f not > 0.2 — the "
                      "§C8 star tooth is BLIND (t=%.2f)", c_plain, ts[it]);

        arb_clear(sm); arb_clear(md);
        aic_dhom_v_clear(&vp); aic_dhom_B_clear(&Bp);
        aic_dhom_v_clear(&v); aic_dhom_B_clear(&B);
        acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae); aic_ucp_kraus_clear(&phi);
    }
    printf("  -> c=md/eta: STAR t=0.06 %.4f t=0.02 %.4f (<< 0.2); PLAIN %.4f %.4f "
           "(> 0.2): ~%.0fx star/plain gap\n", cstar_s[0], cstar_s[1], cplain_s[0],
           cplain_s[1], cplain_s[0] / cstar_s[0]);
}

int main(void)
{
    test_s0_ha_smoke();
    test_t1_eta0_base();
    test_t2_eta0_inductive();
    test_t4_nonsym_u1();
    test_t3a_wrapper_inductive();
    test_t3b_oblique_c8_tooth();
    aic_test_report("test_cstar_extension");
    return 0;
}
