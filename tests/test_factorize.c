/* test_factorize.c — HEADLINE cross-checks for factorize Increment F1 (bead
 * aic-tff): the non-UCP "tilde" maps Delta~ = iota o v and Upsilon~ = v^{-1} o Phi~
 * of th_factorization (approximate_algebras.tex:2749-2752), plus the exact-by-
 * construction identities tilde_DelUps (.tex:2751). See include/aic_factorize.h
 * (the F1 banner), src/aic_factorize.c. Mirrors the test_cstar_build.c setup
 * (fixture channel -> aic_assoc_ecstar_from_phi -> aic_cstar_build -> (B,v) ->
 * factorize struct), the eta_proxy helper, CPREC=256, the aic_test.h macros.
 *
 * THE TWO IDENTITIES (.tex:2751, exact by construction):
 *   Delta~ Upsilon~ = Phi~  : as maps M_N -> M_N (because v o v^{-1} = id_A and
 *     Phi~ maps into A). At eta=0, Phi~ = Phi, so this is a gauge-invariant cross-
 *     check against the ORIGINAL channel (T1); at eta>0 compare to Phi~ = S_tilde
 *     (T2), NOT to Phi.
 *   Upsilon~ Delta~ = 1_B    : as maps B -> B (because Phi~|_A = id and v^{-1} o v
 *     = id_B). Both at eta=0 (T1) and eta>0 (T2).
 *
 * THE STAR / OBLIQUE-PROJECTION IS LOAD-BEARING (FINDINGS §C2/§C3). Upsilon~
 * routes W through Phi~ = the OBLIQUE projector onto A = Img Phi~ (NOT the
 * Frobenius Pi_A). The eta=0 oracle has Phi~ = Phi = id-on-A so it is BLIND to a
 * dropped-Phi~ bug (§C2); only the eta>0 T2 catches it (the mutation tooth below,
 * recorded by hand). T1's Delta~Upsilon~=Phi cross-check is gauge-invariant: it
 * compares MAPS (not the gauge-dependent bases v/idemp differ by a unitary).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
#include "aic_cbnorm.h"
#include "aic_cstar.h"
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_factorize.h"
#include "aic_mat.h"
#include "aic_test.h"
#include "aic_ucp.h"
#include "test_idemp.h"

#define CPREC 256

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Suppress -Wunused-function for the test_idemp.h builders this TU does not use. */
__attribute__((unused))
static void suppress_unused_idemp_builders(void)
{
    (void) make_trace_replace;
    (void) make_identity;
    (void) make_dephasing;
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op via the §C5 midpoint route (mirrors
 * test_cstar_build.c eta_proxy). */
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

/* ||A - B||_op as a double (midpoint route, no Gram false-fail). */
static double opdiff(const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong r = acb_mat_nrows(A), c = acb_mat_ncols(A);
    acb_mat_t D, M;
    acb_mat_init(D, r, c);
    acb_mat_init(M, r, c);
    acb_mat_sub(D, A, B, prec);
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

/* ===================== P-trace: the D1 partial-trace direction pin =========== *
 * THE #1 F3 RISK (blueprint D1, the cbnorm/O2 partial-trace-direction bug class).
 * lem_RC's C_j is
 *   C_j = (1/d_{L_j}) Tr_{L_j} R_j,   R_j = 1_{L_j} (x) C_j  in B(L_j (x) E_j),
 * with L_j the LEFT/MAJOR factor (forced by the U_{js} (x) 1_{E_j} Kroneckers).
 * So C_j = (1/d_{L_j}) aic_mat_partial_trace_LEFT(R_j, a=d_{L_j}, b=d_{E_j})
 * (traces factor 1 = the LEFT/C^a, leaving the b x b RIGHT = E_j). With a wrong
 * partial_trace_RIGHT the result is the WRONG SHAPE (d_{L_j} x d_{L_j}) and the
 * WRONG value — caught instantly by an ASYMMETRIC fixture d_{L_j} != d_{E_j}.
 * Build R = 1_3 (x) C_true, C_true = diag(1, 0.5): correct trace_left -> 2x2
 * = C_true; wrong trace_right -> 3x3 = Tr(C_true)*I_3 = 1.5*I_3 (shape + value
 * both wrong). Bake as a regression BEFORE any L_j code. */
static void test_ptrace_pin(void)
{
    printf("P-trace (D1, the #1 risk): asymmetric d_Lj=3 != d_Ej=2 partial trace\n");
    slong a = 3, b = 2;             /* d_{L_j} = 3 (LEFT), d_{E_j} = 2 (RIGHT) */
    acb_mat_t Ctrue, IL, R, outL, outR;
    acb_mat_init(Ctrue, b, b);
    acb_mat_init(IL, a, a);
    acb_mat_init(R, a * b, a * b);
    acb_mat_init(outL, b, b);       /* trace_left keeps the b x b RIGHT (E_j)  */
    acb_mat_init(outR, a, a);       /* trace_right keeps the a x a LEFT (wrong) */
    acb_set_d(acb_mat_entry(Ctrue, 0, 0), 1.0);
    acb_set_d(acb_mat_entry(Ctrue, 1, 1), 0.5);
    acb_mat_one(IL);
    aic_mat_kronecker(R, IL, Ctrue, CPREC);          /* R = 1_3 (x) C_true */

    /* CORRECT: (1/a) Tr_left(R) over the LEFT factor a -> 2x2 == C_true. */
    aic_mat_partial_trace_left(outL, R, a, b, CPREC);
    {
        acb_t inv;
        acb_init(inv);
        acb_set_d(inv, 1.0 / (double) a);
        acb_mat_scalar_mul_acb(outL, outL, inv, CPREC);
        acb_clear(inv);
    }
    double cdiff = opdiff(outL, Ctrue, CPREC);

    /* WRONG direction: Tr_right(R) over the RIGHT factor b -> 3x3 = Tr(C_true) I_3
     * = 1.5 I_3. Shape (3x3 != 2x2) AND value are both wrong. */
    aic_mat_partial_trace_right(outR, R, a, b, CPREC);

    printf("  trace_LEFT  -> %ldx%ld, ||(1/%ld)Tr_L R - C_true||=%.3e (expect ~0)\n",
           (long) acb_mat_nrows(outL), (long) acb_mat_ncols(outL), (long) a, cdiff);
    printf("  trace_RIGHT -> %ldx%ld (WRONG shape: != %ldx%ld); diag ~ Tr(C_true)=%.2f\n",
           (long) acb_mat_nrows(outR), (long) acb_mat_ncols(outR), (long) b, (long) b,
           dd(acb_realref(acb_mat_entry(outR, 0, 0))));

    AIC_CHECK_MSG(acb_mat_nrows(outL) == b && acb_mat_ncols(outL) == b,
                  "P-trace: trace_left must return %ldx%ld (E_j)", (long) b, (long) b);
    AIC_CHECK_MSG(cdiff < 1e-12,
                  "P-trace: (1/d_Lj) Tr_left(1(x)C) != C_true (=%.3e) — WRONG "
                  "partial-trace direction (D1, the #1 risk)", cdiff);
    AIC_CHECK_MSG(acb_mat_nrows(outR) == a,
                  "P-trace: trace_right returns the WRONG shape (sanity)");

    acb_mat_clear(outR); acb_mat_clear(outL); acb_mat_clear(R);
    acb_mat_clear(IL); acb_mat_clear(Ctrue);
}

/* ===================== P-cent: the 1-design CENTRALITY pin ================== *
 * lem_RC(i) is the FIRST place the genuine 1-design CENTRALITY (diag_j2,
 * .tex:2776) is load-bearing: R_j = 1_{L_j} (x) C_j requires
 *   D_j(X) := Sum_s p_{js} X U_{js}^dag (x) U_{js}
 *           = Sum_s p_{js} U_{js}^dag (x) U_{js} X     for all X in B(L_j).
 * F2's tests (Delta' CP, ||Delta'-Delta~||) are BLIND to this — they hold for any
 * unitary set with Sum p_s = 1. This pin verifies the per-block Pauli design (the
 * d_j^2 Paulis S_{ab} = X^a Z^b, weight d_j^{-2}, aic_dhom_pauli) is central on a
 * NON-scalar X0. Mutation (drop a Pauli / perturb a weight) -> O(1) RED. */
static double design_centrality_defect(slong d, const acb_mat_t X0, slong prec)
{
    /* L = Sum_s p_s X0 U_s^dag (x) U_s ; Rt = Sum_s p_s U_s^dag (x) U_s X0. */
    acb_mat_t L, Rt, U, Ud, X0U, X0Ud, T;
    acb_mat_init(L, d * d, d * d);
    acb_mat_init(Rt, d * d, d * d);
    acb_mat_init(U, d, d);
    acb_mat_init(Ud, d, d);
    acb_mat_init(X0Ud, d, d);   /* X0 U^dag (LEFT factor of the L term)       */
    acb_mat_init(X0U, d, d);    /* U X0    (RIGHT factor of the R term)       */
    acb_mat_init(T, d * d, d * d);
    acb_mat_zero(L);
    acb_mat_zero(Rt);
    arb_t ps;
    arb_init(ps);
    arb_set_si(ps, 1);
    arb_div_si(ps, ps, d * d, prec);     /* p_s = 1/d^2 */
    for (slong qa = 0; qa < d; qa++)
        for (slong qb = 0; qb < d; qb++) {
            aic_dhom_pauli(U, d, qa, qb, prec);          /* U_s = S_{ab} */
            acb_mat_conjugate_transpose(Ud, U);          /* U_s^dag */
            /* L term: (X0 U_s^dag) (x) U_s */
            acb_mat_mul(X0Ud, X0, Ud, prec);
            aic_mat_kronecker(T, X0Ud, U, prec);
            acb_mat_scalar_mul_arb(T, T, ps, prec);
            acb_mat_add(L, L, T, prec);
            /* R term: U_s^dag (x) (U_s X0) */
            acb_mat_mul(X0U, U, X0, prec);               /* U_s X0 */
            aic_mat_kronecker(T, Ud, X0U, prec);
            acb_mat_scalar_mul_arb(T, T, ps, prec);
            acb_mat_add(Rt, Rt, T, prec);
        }
    double diff = opdiff(L, Rt, prec);
    arb_clear(ps);
    acb_mat_clear(T); acb_mat_clear(X0U); acb_mat_clear(X0Ud);
    acb_mat_clear(Ud); acb_mat_clear(U); acb_mat_clear(Rt); acb_mat_clear(L);
    return diff;
}

static void test_centrality_pin(void)
{
    printf("P-cent (lem_RC centrality diag_j2): per-block Pauli design is central\n");
    slong d = 3;
    acb_mat_t X0;
    acb_mat_init(X0, d, d);
    /* a NON-scalar Hermitian X0 = diag(1,2,3) + a real off-diagonal coupling. */
    acb_set_si(acb_mat_entry(X0, 0, 0), 1);
    acb_set_si(acb_mat_entry(X0, 1, 1), 2);
    acb_set_si(acb_mat_entry(X0, 2, 2), 3);
    acb_set_d(acb_mat_entry(X0, 0, 1), 0.5);
    acb_set_d(acb_mat_entry(X0, 1, 0), 0.5);
    double c = design_centrality_defect(d, X0, CPREC);
    printf("  d=%ld non-scalar X0: ||D_j(X0)_left - D_j(X0)_right|| = %.3e\n",
           (long) d, c);
    AIC_CHECK_MSG(c < 1e-9,
                  "P-cent: per-block Pauli design NOT central (=%.3e) — diag_j2 "
                  "broken (dropped Pauli / wrong weight?)", c);
    acb_mat_clear(X0);
}

/* Max over B's matrix units E_i of ||Upsilon~(Delta~(E_i)) - E_i||_op
 * (Upsilon~ Delta~ = 1_B, .tex:2751). */
static double updel_defect(const aic_factorize *F, slong prec)
{
    slong nB = F->n_B, dimB = F->v->B->dim_B;
    acb_mat_t Ei, DX, UDX;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(DX, F->N, F->N);
    acb_mat_init(UDX, nB, nB);
    double worst = 0.0;
    for (slong i = 0; i < dimB; i++) {
        aic_dhom_B_matunit(Ei, F->v->B, i);
        aic_factorize_delta_tilde(DX, F, Ei, prec);
        aic_factorize_upsilon_tilde(UDX, F, DX, prec);
        double d = opdiff(UDX, Ei, prec);
        if (d > worst) worst = d;
    }
    acb_mat_clear(UDX); acb_mat_clear(DX); acb_mat_clear(Ei);
    return worst;
}

/* Max over M_N matrix units E_pq of ||Delta~(Upsilon~(E_pq)) - target(E_pq)||_op,
 * where target is either Phi (eta=0, applied via aic_ucp_apply) or Phi~ = S_tilde
 * (eta>0, applied via aic_assoc_superop_apply). Delta~ Upsilon~ = Phi~ (.tex:2751). */
static double delup_defect_vs_phi(const aic_factorize *F, const aic_ucp_kraus *phi,
                                  slong prec)
{
    slong N = F->N;
    acb_mat_t Epq, UE, DUE, T;
    acb_mat_init(Epq, N, N);
    acb_mat_init(UE, F->n_B, F->n_B);
    acb_mat_init(DUE, N, N);
    acb_mat_init(T, N, N);
    double worst = 0.0;
    for (slong p = 0; p < N; p++)
        for (slong q = 0; q < N; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_upsilon_tilde(UE, F, Epq, prec);
            aic_factorize_delta_tilde(DUE, F, UE, prec);
            aic_ucp_apply(T, phi, Epq, prec);   /* Phi(E_pq), eta=0: Phi~ = Phi */
            double d = opdiff(DUE, T, prec);
            if (d > worst) worst = d;
        }
    acb_mat_clear(T); acb_mat_clear(DUE); acb_mat_clear(UE); acb_mat_clear(Epq);
    return worst;
}

static double delup_defect_vs_phitilde(const aic_factorize *F, slong prec)
{
    slong N = F->N;
    acb_mat_t Epq, UE, DUE, T;
    acb_mat_init(Epq, N, N);
    acb_mat_init(UE, F->n_B, F->n_B);
    acb_mat_init(DUE, N, N);
    acb_mat_init(T, N, N);
    double worst = 0.0;
    for (slong p = 0; p < N; p++)
        for (slong q = 0; q < N; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_upsilon_tilde(UE, F, Epq, prec);
            aic_factorize_delta_tilde(DUE, F, UE, prec);
            aic_assoc_superop_apply(T, F->Aec->S_tilde, Epq, prec); /* Phi~(E_pq) */
            double d = opdiff(DUE, T, prec);
            if (d > worst) worst = d;
        }
    acb_mat_clear(T); acb_mat_clear(DUE); acb_mat_clear(UE); acb_mat_clear(Epq);
    return worst;
}

/* ===================== T1: eta=0 oracle ===================================== *
 * On EXACT-idempotent channels (Phi~ = Phi, the cleanest ground truth, Rule 6):
 *   - Upsilon~ Delta~ = 1_B as a MAP (sweep B matrix units), ~1e-10.
 *   - Delta~ Upsilon~ = Phi as a MAP (gauge-invariant cross-check vs the ORIGINAL
 *     channel; sweep M_N matrix units), ~1e-10. */
static void t1_one(const char *name, aic_ucp_kraus *phi)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);

    aic_factorize F;
    aic_factorize_build(&F, &v, &ae, phi, CPREC);

    double ud = updel_defect(&F, CPREC);
    double du = delup_defect_vs_phi(&F, phi, CPREC);
    printf("T1 %s: N=%ld n_B=%ld dim_A=%ld | UpsDel-1_B=%.3e DelUps-Phi=%.3e "
           "iso_def=%.3e\n", name, (long) F.N, (long) F.n_B, (long) F.dim_A, ud,
           du, dd(iso));

    AIC_CHECK_MSG(ud < 1e-10, "T1 %s: Upsilon~Delta~ - 1_B = %.3e not ~0", name, ud);
    AIC_CHECK_MSG(du < 1e-10, "T1 %s: Delta~Upsilon~ - Phi = %.3e not ~0", name, du);

    aic_factorize_clear(&F);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t1_eta0_oracle(void)
{
    printf("T1 eta=0 oracle (Upsilon~Delta~=1_B, Delta~Upsilon~=Phi as MAPS):\n");
    aic_ucp_kraus phi;

    make_block_cond_exp(&phi, 4, 2);          /* M_2 (+) M_2, dim_A=8 */
    t1_one("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    make_noiseless_subsystem(&phi, 2, 2);     /* M_2, dim_A=4 */
    t1_one("noiseless_subsystem(2,2)", &phi);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T2: eta>0 exact-by-construction ===================== *
 * On a GENUINELY oblique eta>0 channel (make_mixconj): the tilde maps' identities
 * are exact by construction (limited only by the v-inversion + S_tilde idempotence
 * precision). The OBLIQUE star path (§C2/§C3) is exercised here:
 *   - Delta~ Upsilon~ = Phi~ as maps (compare to S_tilde, NOT Phi: at eta>0
 *     Phi~ != Phi). Should be ~1e-10 or better.
 *   - Upsilon~ Delta~ = 1_B as a map. Should be ~machine. */
static void test_t2_eta_gt0(void)
{
    printf("T2 eta>0 exact-by-construction (oblique star path, §C2/§C3):\n");
    struct { slong d, m; double t; } fx[] = {
        {4, 2, 0.03}, {5, 2, 0.03}};
    const int n_fx = (int) (sizeof(fx) / sizeof(fx[0]));
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

        aic_factorize F;
        aic_factorize_build(&F, &v, &ae, &phi, CPREC);

        double ud = updel_defect(&F, CPREC);
        double du = delup_defect_vs_phitilde(&F, CPREC);
        printf("  mixconj(%ld,%ld,%.2f): N=%ld n_B=%ld dim_A=%ld eta=%.3e | "
               "UpsDel-1_B=%.3e DelUps-Phi~=%.3e iso_def=%.3e\n",
               (long) fx[i].d, (long) fx[i].m, fx[i].t, (long) F.N, (long) F.n_B,
               (long) F.dim_A, eta, ud, du, dd(iso));

        /* exact by construction: tied to inversion + S_tilde idempotence prec. */
        AIC_CHECK_MSG(du < 1e-9,
                      "T2 mixconj(%ld,%ld): Delta~Upsilon~ - Phi~ = %.3e not exact "
                      "(oblique star path broken? §C2/§C3)", (long) fx[i].d,
                      (long) fx[i].m, du);
        AIC_CHECK_MSG(ud < 1e-9,
                      "T2 mixconj(%ld,%ld): Upsilon~Delta~ - 1_B = %.3e not exact",
                      (long) fx[i].d, (long) fx[i].m, ud);

        aic_factorize_clear(&F);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ===================== F2 helpers (Delta', Delta, the CP teeth) ============= */

/* Max over B matrix units E_i of ||Delta'(E_i) - Delta~(E_i)||_op. At eta=0 this
 * is ~machine (Delta' = Delta~, the .tex:2797 O(eta) bound collapses); at eta>0 it
 * is O(eta) (tilde_Del2, .tex:2797). */
static double dprime_vs_dtilde(const aic_factorize *F, slong prec)
{
    slong nB = F->n_B, dimB = F->v->B->dim_B;
    acb_mat_t Ei, DP, DT;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(DP, F->N, F->N);
    acb_mat_init(DT, F->N, F->N);
    double worst = 0.0;
    for (slong i = 0; i < dimB; i++) {
        aic_dhom_B_matunit(Ei, F->v->B, i);
        aic_factorize_delta_prime(DP, F, Ei, prec);
        aic_factorize_delta_tilde(DT, F, Ei, prec);
        double d = opdiff(DP, DT, prec);
        if (d > worst) worst = d;
    }
    acb_mat_clear(DT); acb_mat_clear(DP); acb_mat_clear(Ei);
    return worst;
}

/* Max over B matrix units of ||Delta(E_i) - Delta~(E_i)||_op (requires
 * delta_build; .tex:2801 ||Delta - Delta~||_cb <= O(eta)). */
static double delta_vs_dtilde(const aic_factorize *F, slong prec)
{
    slong nB = F->n_B, dimB = F->v->B->dim_B;
    acb_mat_t Ei, D, DT;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(D, F->N, F->N);
    acb_mat_init(DT, F->N, F->N);
    double worst = 0.0;
    for (slong i = 0; i < dimB; i++) {
        aic_dhom_B_matunit(Ei, F->v->B, i);
        aic_factorize_delta(D, F, Ei, prec);
        aic_factorize_delta_tilde(DT, F, Ei, prec);
        double d = opdiff(D, DT, prec);
        if (d > worst) worst = d;
    }
    acb_mat_clear(DT); acb_mat_clear(D); acb_mat_clear(Ei);
    return worst;
}

/* ||Delta'(I_B) - 1_H||_op (.tex:2797: = O(eta); shard H #7: must be < 1). */
static double dprimeI_minus_1(const aic_factorize *F, slong prec)
{
    acb_mat_t IB, DI, I;
    acb_mat_init(IB, F->n_B, F->n_B);
    acb_mat_init(DI, F->N, F->N);
    acb_mat_init(I, F->N, F->N);
    aic_dhom_B_unit(IB, F->v->B);
    aic_factorize_delta_prime(DI, F, IB, prec);
    double d = (acb_mat_one(I), opdiff(DI, I, prec));
    acb_mat_clear(I); acb_mat_clear(DI); acb_mat_clear(IB);
    return d;
}

/* ||Delta(I_B) - 1_H||_op (UCP unitality; should be ~machine by construction). */
static double deltaI_minus_1(const aic_factorize *F, slong prec)
{
    acb_mat_t IB, D, I;
    acb_mat_init(IB, F->n_B, F->n_B);
    acb_mat_init(D, F->N, F->N);
    acb_mat_init(I, F->N, F->N);
    aic_dhom_B_unit(IB, F->v->B);
    aic_factorize_delta(D, F, IB, prec);
    double d = (acb_mat_one(I), opdiff(D, I, prec));
    acb_mat_clear(I); acb_mat_clear(D); acb_mat_clear(IB);
    return d;
}

/* CP-check of Delta' via the per-block Choi C_j (.tex:2791-2795), the §C2-style
 * teeth: Delta' CP iff every block Choi C_j is PSD. Computed on the matmul-
 * accumulated Choi, which the stock aic_ucp_is_cp_choi CANNOT verdict (FINDINGS §C5:
 * its strict relative-Hermiticity guard false-fails — ||C-C^dag||_op ~ 4e-16 from
 * differing arithmetic paths per entry vastly exceeds the strict 2^(-prec/2) ~ 3e-39
 * tol; symmetrizing is insufficient since arb subtraction adds the radii). So we
 * extract the EXACTLY-Hermitian midpoint Cmid = (mid(C)+mid(C)^dag)/2 (zero radius)
 * and a Weyl inflation R >= ||C - Cmid||_op (Frobenius upper-bounds op-norm), giving
 *   mineig(C_j) >= mineig(Cmid) - R = -maxeig(-Cmid) - R   (.tex Weyl).
 *
 * GATE-MIDPOINT POSTURE (NOT formally endpoint-rigorous), and WHY. The verdict reads
 * the MIDPOINT dd(lb) of the lower-bound ball lb = -maxeig(-Cmid) - R, the project's
 * standard coarse fail-loud gate (cf. the projection-nontriviality / sigma_min / O1
 * HOPM gates). It is NOT the ball's lower endpoint: at eta=0 the Choi is DEGENERATE
 * (spectrum {1..1,0..0}), so aic_mat_herm_max_eig's certified enclosure on the
 * near-zero cluster is WIDE (the aic-w4o.1 degenerate-eig wall) — its lower endpoint
 * dips below -tol even though the true mineig is 0, while the midpoint is well-
 * centred (~1e-71). The gate is sound with a VAST margin: at eta>0 the CP signal is
 * +9.4e-6 / +2.2e-5 vs R ~ 3.7e-16 (~1e10x), and the Delta~-NOT-CP witness is
 * -2.1e-3 (~1e7x the radius) — the decision never lives near the radius. A certified
 * endpoint verdict awaits the aic-w4o.1 degenerate-eig enclosure. C_j PSD iff
 * dd(lb) >= -tol. Returns 1 iff every block passes; *worst_mineig = the most-negative
 * gate-midpoint lower bound. */
static int dprime_is_cp(const aic_factorize *F, double tol, double *worst_mineig,
                        slong prec)
{
    const aic_dhom_B *B = F->v->B;
    slong N = F->N;
    int all_cp = 1;
    double wm = 1e300;
    for (slong j = 0; j < B->num_blocks; j++) {
        slong dj = B->d[j], M = dj * N;
        acb_mat_t Cj, Cmid, Cmd, negCmid, resid;
        acb_mat_init(Cj, M, M);
        acb_mat_init(Cmid, M, M);
        acb_mat_init(Cmd, M, M);
        acb_mat_init(negCmid, M, M);
        acb_mat_init(resid, M, M);
        aic_factorize_delta_block_choi(Cj, F, j, prec);

        /* Cmid = exactly-Hermitian midpoint of C_j. */
        for (slong p = 0; p < M; p++)
            for (slong q = 0; q < M; q++)
                acb_get_mid(acb_mat_entry(Cmid, p, q), acb_mat_entry(Cj, p, q));
        acb_mat_conjugate_transpose(Cmd, Cmid);
        acb_mat_add(Cmid, Cmid, Cmd, prec);
        acb_mat_scalar_mul_2exp_si(Cmid, Cmid, -1);     /* (mid + mid^dag)/2 */

        /* R >= ||C_j - Cmid||_op (rigorous; Frobenius >= op-norm). */
        arb_t R, mx, lb;
        arb_init(R);
        arb_init(mx);
        arb_init(lb);
        acb_mat_sub(resid, Cj, Cmid, prec);
        aic_mat_frobenius_norm(R, resid, prec);

        acb_mat_neg(negCmid, Cmid);
        aic_mat_herm_max_eig(mx, negCmid, prec);        /* maxeig(-Cmid) */
        /* mineig(C_j) >= -maxeig(-Cmid) - R */
        arb_neg(lb, mx);
        arb_sub(lb, lb, R, prec);
        double mineig_lb = dd(lb);   /* gate midpoint — see the docstring (aic-w4o.1) */
        if (mineig_lb < wm) wm = mineig_lb;
        if (mineig_lb < -tol) all_cp = 0;

        arb_clear(lb); arb_clear(mx); arb_clear(R);
        acb_mat_clear(resid); acb_mat_clear(negCmid); acb_mat_clear(Cmd);
        acb_mat_clear(Cmid); acb_mat_clear(Cj);
    }
    *worst_mineig = wm;
    return all_cp;
}

/* The NON-VACUITY guard for the CP-ization: the most-negative min-eigenvalue
 * (midpoint estimate) over the per-block Choi matrices of Delta~ ITSELF (using
 * Delta~(E_ab) in place of Delta'(E_ab)). At eta>0 Delta~ = v is a *-hom into A
 * but NOT a CP B(H)-map, so its block Choi has a NEGATIVE eigenvalue — this is
 * what the 1-design average Delta' fixes. If this were >= 0 the "Delta' CP"
 * headline would be vacuous (§C2). Midpoint estimate (no rigorous Weyl needed: we
 * only need to witness a clearly-negative eigenvalue, not certify a boundary). */
static double dtilde_block_choi_worst_mineig(const aic_factorize *F, slong prec)
{
    const aic_dhom_B *B = F->v->B;
    slong N = F->N;
    double wm = 1e300;
    for (slong j = 0; j < B->num_blocks; j++) {
        slong dj = B->d[j], M = dj * N;
        acb_mat_t Cj, Eab, DE, Cmid, Cmd, negCmid;
        acb_mat_init(Cj, M, M);
        acb_mat_init(Eab, B->n_B, B->n_B);
        acb_mat_init(DE, N, N);
        acb_mat_init(Cmid, M, M);
        acb_mat_init(Cmd, M, M);
        acb_mat_init(negCmid, M, M);
        acb_mat_zero(Cj);
        for (slong a = 0; a < dj; a++)
            for (slong b = 0; b < dj; b++) {
                slong idx = aic_dhom_B_index(B, j, a, b);
                aic_dhom_B_matunit(Eab, B, idx);
                aic_factorize_delta_tilde(DE, F, Eab, prec);   /* Delta~, NOT Delta' */
                for (slong p = 0; p < N; p++)
                    for (slong q = 0; q < N; q++)
                        acb_set(acb_mat_entry(Cj, a * N + p, b * N + q),
                                acb_mat_entry(DE, p, q));
            }
        for (slong p = 0; p < M; p++)
            for (slong q = 0; q < M; q++)
                acb_get_mid(acb_mat_entry(Cmid, p, q), acb_mat_entry(Cj, p, q));
        acb_mat_conjugate_transpose(Cmd, Cmid);
        acb_mat_add(Cmid, Cmid, Cmd, prec);
        acb_mat_scalar_mul_2exp_si(Cmid, Cmid, -1);
        arb_t mx;
        arb_init(mx);
        acb_mat_neg(negCmid, Cmid);
        aic_mat_herm_max_eig(mx, negCmid, prec);
        double mineig = -dd(mx);
        if (mineig < wm) wm = mineig;
        arb_clear(mx);
        acb_mat_clear(negCmid); acb_mat_clear(Cmd); acb_mat_clear(Cmid);
        acb_mat_clear(DE); acb_mat_clear(Eab); acb_mat_clear(Cj);
    }
    return wm;
}

/* Delta is a *-hom: max over B matrix-unit pairs of ||Phi(Delta(X)Delta(Y)) -
 * Delta(XY)||_op (PhiDelta2, .tex:2810; at eta=0 -> machine). */
static double delta_phihom_defect(const aic_factorize *F, slong prec)
{
    const aic_dhom_B *B = F->v->B;
    slong nB = F->n_B, dimB = B->dim_B, N = F->N;
    acb_mat_t Ei, Ej, EiEj, DX, DY, prod, phiprod, DXY;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(Ej, nB, nB);
    acb_mat_init(EiEj, nB, nB);
    acb_mat_init(DX, N, N);
    acb_mat_init(DY, N, N);
    acb_mat_init(prod, N, N);
    acb_mat_init(phiprod, N, N);
    acb_mat_init(DXY, N, N);
    double worst = 0.0;
    for (slong i = 0; i < dimB; i++)
        for (slong k = 0; k < dimB; k++) {
            aic_dhom_B_matunit(Ei, B, i);
            aic_dhom_B_matunit(Ej, B, k);
            aic_dhom_B_mul(EiEj, B, Ei, Ej, prec);     /* XY in B */
            aic_factorize_delta(DX, F, Ei, prec);
            aic_factorize_delta(DY, F, Ej, prec);
            acb_mat_mul(prod, DX, DY, prec);
            aic_ucp_apply(phiprod, F->phi, prod, prec); /* Phi(Delta(X)Delta(Y)) */
            aic_factorize_delta(DXY, F, EiEj, prec);    /* Delta(XY) */
            double d = opdiff(phiprod, DXY, prec);
            if (d > worst) worst = d;
        }
    acb_mat_clear(DXY); acb_mat_clear(phiprod); acb_mat_clear(prod);
    acb_mat_clear(DY); acb_mat_clear(DX); acb_mat_clear(EiEj);
    acb_mat_clear(Ej); acb_mat_clear(Ei);
    return worst;
}

/* ===================== T3: eta=0 oracle (Delta reduces to Delta~) =========== *
 * On EXACT-idempotent channels (Phi~ = Phi). Because Delta~ = v is an exact *-hom
 * into A and Phi projects onto A:
 *   Phi(Delta~(XU_s^dag)Delta~(U_s)) = Delta~(XU_s^dag) * Delta~(U_s)
 *     = Delta~(XU_s^dag U_s) = Delta~(X),   sum_s p_s = 1  =>  Delta' = Delta~.
 * So Delta'(E_i) == Delta~(E_i) to machine; ||Delta'(I_B)-1_H|| ~ machine; the
 * per-block Choi C_j is PSD; Delta is UCP and == Delta~. */
static void t3_one(const char *name, aic_ucp_kraus *phi)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);

    aic_factorize F;
    aic_factorize_build(&F, &v, &ae, phi, CPREC);

    double dpdt = dprime_vs_dtilde(&F, CPREC);
    double dpI = dprimeI_minus_1(&F, CPREC);
    double cp_mineig;
    int cp = dprime_is_cp(&F, 1e-9, &cp_mineig, CPREC);

    aic_factorize_delta_build(&F, CPREC);
    double dI = deltaI_minus_1(&F, CPREC);
    double ddt = delta_vs_dtilde(&F, CPREC);
    double hom = delta_phihom_defect(&F, CPREC);

    printf("T3 %s: N=%ld n_B=%ld | Delta'-Delta~=%.3e ||Delta'(I)-1||=%.3e "
           "CP=%d minEig(C_j)=%.3e | Delta(I)-1=%.3e Delta-Delta~=%.3e "
           "PhiHom=%.3e\n", name, (long) F.N, (long) F.n_B, dpdt, dpI, cp,
           cp_mineig, dI, ddt, hom);

    AIC_CHECK_MSG(dpdt < 1e-10, "T3 %s: Delta'-Delta~ = %.3e not ~0", name, dpdt);
    AIC_CHECK_MSG(dpI < 1e-10, "T3 %s: ||Delta'(I_B)-1_H|| = %.3e not ~0", name, dpI);
    AIC_CHECK_MSG(cp == 1, "T3 %s: Delta' not CP (some C_j not PSD)", name);
    AIC_CHECK_MSG(dI < 1e-10, "T3 %s: Delta(I_B)-1_H = %.3e not ~0", name, dI);
    AIC_CHECK_MSG(ddt < 1e-10, "T3 %s: Delta-Delta~ = %.3e not ~0", name, ddt);
    AIC_CHECK_MSG(hom < 1e-9, "T3 %s: Delta not a *-hom (PhiDelta2=%.3e)", name, hom);

    aic_factorize_clear(&F);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t3_eta0_oracle(void)
{
    printf("T3 eta=0 oracle (Delta reduces to Delta~; Delta' CP; Delta UCP):\n");
    aic_ucp_kraus phi;

    make_block_cond_exp(&phi, 4, 2);          /* M_2 (+) M_2, 16-term design */
    t3_one("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    make_noiseless_subsystem(&phi, 2, 2);     /* M_2, 4-term design */
    t3_one("noiseless_subsystem(2,2)", &phi);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T4: eta>0 (the CP-ization payoff) ==================== *
 * On make_mixconj (genuinely oblique eta>0): Delta~ is NOT CP, the 1-design fixes
 * it. The headline: each per-block Choi C_j is PSD (.tex:2791-2795). And
 *   ||Delta' - Delta~|| = O(eta)  (tilde_Del2, .tex:2797),
 *   ||Delta'(I_B) - 1_H|| = O(eta) < 1,
 *   Delta is UCP (Delta(I_B)=1_H, CP), ||Delta - Delta~|| = O(eta). */
static void test_t4_eta_gt0(void)
{
    printf("T4 eta>0 (CP-ization payoff: Delta~ NOT CP, the 1-design fixes it):\n");
    struct { slong d, m; double t; } fx[] = {
        {4, 2, 0.03}, {5, 2, 0.03}};
    const int n_fx = (int) (sizeof(fx) / sizeof(fx[0]));
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

        aic_factorize F;
        aic_factorize_build(&F, &v, &ae, &phi, CPREC);

        double dt_mineig = dtilde_block_choi_worst_mineig(&F, CPREC); /* NON-VACUITY */
        double dpdt = dprime_vs_dtilde(&F, CPREC);
        double dpI = dprimeI_minus_1(&F, CPREC);
        double cp_mineig;
        int cp = dprime_is_cp(&F, 1e-9, &cp_mineig, CPREC);

        aic_factorize_delta_build(&F, CPREC);
        double dI = deltaI_minus_1(&F, CPREC);
        double ddt = delta_vs_dtilde(&F, CPREC);

        printf("  mixconj(%ld,%ld,%.2f): N=%ld n_B=%ld eta=%.3e | Delta~ NOT CP "
               "(minEig=%.3e) -> Delta' CP=%d (minEig(C_j)=%.3e) | "
               "Delta'-Delta~=%.3e (c=%.2f) ||Delta'(I)-1||=%.3e (c=%.2f) | "
               "Delta(I)-1=%.3e Delta-Delta~=%.3e (c=%.2f)\n",
               (long) fx[i].d, (long) fx[i].m, fx[i].t, (long) F.N, (long) F.n_B,
               eta, dt_mineig, cp, cp_mineig, dpdt, dpdt / eta, dpI, dpI / eta, dI,
               ddt, ddt / eta);

        /* NON-VACUITY (§C2): Delta~ must be genuinely NOT CP on this fixture, else
         * "Delta' CP" is a free pass. The 1-design average is the fix. */
        AIC_CHECK_MSG(dt_mineig < -1e-4,
                      "T4 mixconj(%ld,%ld): Delta~ block Choi minEig=%.3e not "
                      "clearly negative — CP-ization is VACUOUS (fixture too mild)",
                      (long) fx[i].d, (long) fx[i].m, dt_mineig);
        AIC_CHECK_MSG(cp == 1,
                      "T4 mixconj(%ld,%ld): Delta' NOT CP (minEig(C_j)=%.3e) — the "
                      "1-design failed to CP-ize", (long) fx[i].d, (long) fx[i].m,
                      cp_mineig);
        AIC_CHECK_MSG(dpdt < 5.0 * eta,
                      "T4 mixconj(%ld,%ld): ||Delta'-Delta~|| = %.3e not O(eta) "
                      "(eta=%.3e)", (long) fx[i].d, (long) fx[i].m, dpdt, eta);
        AIC_CHECK_MSG(dpI < 1.0,
                      "T4 mixconj(%ld,%ld): ||Delta'(I_B)-1_H|| = %.3e >= 1",
                      (long) fx[i].d, (long) fx[i].m, dpI);
        AIC_CHECK_MSG(dI < 1e-10,
                      "T4 mixconj(%ld,%ld): Delta(I_B)-1_H = %.3e not ~0 (not unital)",
                      (long) fx[i].d, (long) fx[i].m, dI);
        AIC_CHECK_MSG(ddt < 5.0 * eta,
                      "T4 mixconj(%ld,%ld): ||Delta-Delta~|| = %.3e not O(eta)",
                      (long) fx[i].d, (long) fx[i].m, ddt);

        aic_factorize_clear(&F);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ===================== F3 helpers (Upsilon via lem_RC) ===================== */

/* max_j ||R_j - 1_{L_j} (x) C_j||_op (lem_RC(i), .tex:2846). Certifies lem_RC(i)
 * + the D1 partial-trace direction + centrality SIMULTANEOUSLY: it RED-fires if
 * the per-block design is not central, the trace is the wrong direction, or W_j
 * is mis-stacked. Rebuilds R_j/C_j from the cached W_j (gauge-invariant). Reports
 * the worst sigma_max(C_j) deviation too. */
static double Rj_minus_1xCj_worst(const aic_factorize *F, double *worst_smax_dev,
                                  slong prec)
{
    const aic_dhom_B *B = F->v->B;
    double worst = 0.0, wdev = 0.0;
    for (slong j = 0; j < B->num_blocks; j++) {
        slong dj = B->d[j], e_j = F->e[j], dim = dj * e_j;
        acb_mat_t Rj, IL, oneCj;
        acb_mat_init(Rj, dim, dim);
        acb_mat_init(IL, dj, dj);
        acb_mat_init(oneCj, dim, dim);
        aic_factorize_upsilon_Rj(Rj, F, j, F->W[j], e_j, prec);
        acb_mat_one(IL);
        aic_mat_kronecker(oneCj, IL, F->C[j], prec);   /* 1_{L_j} (x) C_j */
        double d = opdiff(Rj, oneCj, prec);
        if (d > worst) worst = d;
        /* sigma_max(C_j) via opnorm of C_j. */
        arb_t s;
        arb_init(s);
        {
            acb_mat_t Cmid;
            acb_mat_init(Cmid, e_j, e_j);
            for (slong p = 0; p < e_j; p++)
                for (slong q = 0; q < e_j; q++)
                    acb_get_mid(acb_mat_entry(Cmid, p, q), acb_mat_entry(F->C[j], p, q));
            aic_mat_opnorm(s, Cmid, prec);
            acb_mat_clear(Cmid);
        }
        double dev = dd(s) - 1.0;
        if (dev < 0) dev = -dev;
        if (dev > wdev) wdev = dev;
        arb_clear(s);
        acb_mat_clear(oneCj); acb_mat_clear(IL); acb_mat_clear(Rj);
    }
    if (worst_smax_dev) *worst_smax_dev = wdev;
    return worst;
}

/* max over B matrix units E_i of ||Upsilon(Delta(E_i)) - E_i||_op (Upsilon Delta
 * = 1_B, .tex:2871). Δ = the UCP F2 encode; Upsilon = the UCP F3 decode. */
static double ups_del_defect(const aic_factorize *F, slong prec)
{
    slong nB = F->n_B, dimB = F->v->B->dim_B;
    acb_mat_t Ei, DE, UDE;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(DE, F->N, F->N);
    acb_mat_init(UDE, nB, nB);
    double worst = 0.0;
    for (slong i = 0; i < dimB; i++) {
        aic_dhom_B_matunit(Ei, F->v->B, i);
        aic_factorize_delta(DE, F, Ei, prec);       /* Delta(E_i) (UCP) */
        aic_factorize_upsilon(UDE, F, DE, prec);     /* Upsilon(Delta(E_i)) */
        double d = opdiff(UDE, Ei, prec);
        if (d > worst) worst = d;
    }
    acb_mat_clear(UDE); acb_mat_clear(DE); acb_mat_clear(Ei);
    return worst;
}

/* max over M_N matrix units E_pq of ||Delta(Upsilon(E_pq)) - target(E_pq)||_op,
 * target = Phi (eta=0, via aic_ucp_apply) or Phi~ = S_tilde (eta>0). DelUps = Phi
 * (.tex:2734, the headline). */
static double del_ups_defect(const aic_factorize *F, const aic_ucp_kraus *phi,
                             int use_phitilde, slong prec)
{
    slong N = F->N;
    acb_mat_t Epq, UE, DUE, T;
    acb_mat_init(Epq, N, N);
    acb_mat_init(UE, F->n_B, F->n_B);
    acb_mat_init(DUE, N, N);
    acb_mat_init(T, N, N);
    double worst = 0.0;
    for (slong p = 0; p < N; p++)
        for (slong q = 0; q < N; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_upsilon(UE, F, Epq, prec);    /* Upsilon(E_pq) */
            aic_factorize_delta(DUE, F, UE, prec);       /* Delta(Upsilon(E_pq)) */
            if (use_phitilde)
                aic_assoc_superop_apply(T, F->Aec->S_tilde, Epq, prec);
            else
                aic_ucp_apply(T, phi, Epq, prec);        /* Phi(E_pq), eta=0 */
            double d = opdiff(DUE, T, prec);
            if (d > worst) worst = d;
        }
    acb_mat_clear(T); acb_mat_clear(DUE); acb_mat_clear(UE); acb_mat_clear(Epq);
    return worst;
}

/* Upsilon(1_H) - I_B (UCP unitality of Upsilon; ~machine by construction). */
static double upsI_minus_IB(const aic_factorize *F, slong prec)
{
    acb_mat_t IH, U, IB;
    acb_mat_init(IH, F->N, F->N);
    acb_mat_init(U, F->n_B, F->n_B);
    acb_mat_init(IB, F->n_B, F->n_B);
    acb_mat_one(IH);
    aic_factorize_upsilon(U, F, IH, prec);
    aic_dhom_B_unit(IB, F->v->B);
    double d = opdiff(U, IB, prec);
    acb_mat_clear(IB); acb_mat_clear(U); acb_mat_clear(IH);
    return d;
}

/* CP-verdict of Upsilon' via its per-block Choi C^{up}_j (the §C5 midpoint+Weyl
 * gate, reusing F2's dprime_is_cp pattern). Upsilon'_j: B(H) -> B(L_j); its
 * Convention-A Choi is C^{up}_j[(p)*d_j + a, (q)*d_j + b] = Upsilon'_j(E_pq)[a,b]
 * with the H source factor (dim N) MAJOR, the L_j target (dim d_j) MINOR. CP iff
 * every C^{up}_j is PSD. */
static int upsilon_prime_is_cp(const aic_factorize *F, double tol,
                               double *worst_mineig, slong prec)
{
    const aic_dhom_B *B = F->v->B;
    slong N = F->N;
    int all_cp = 1;
    double wm = 1e300;
    for (slong j = 0; j < B->num_blocks; j++) {
        slong dj = B->d[j], M = N * dj;
        acb_mat_t Cj, Epq, UE, Cmid, Cmd, negCmid, resid;
        acb_mat_init(Cj, M, M);
        acb_mat_init(Epq, N, N);
        acb_mat_init(UE, dj, dj);
        acb_mat_init(Cmid, M, M);
        acb_mat_init(Cmd, M, M);
        acb_mat_init(negCmid, M, M);
        acb_mat_init(resid, M, M);
        acb_mat_zero(Cj);
        for (slong p = 0; p < N; p++)
            for (slong q = 0; q < N; q++) {
                acb_mat_zero(Epq);
                acb_set_si(acb_mat_entry(Epq, p, q), 1);
                aic_factorize_upsilon_prime_block(UE, F, j, Epq, prec);
                for (slong a = 0; a < dj; a++)
                    for (slong b = 0; b < dj; b++)
                        acb_set(acb_mat_entry(Cj, p * dj + a, q * dj + b),
                                acb_mat_entry(UE, a, b));
            }
        for (slong p = 0; p < M; p++)
            for (slong q = 0; q < M; q++)
                acb_get_mid(acb_mat_entry(Cmid, p, q), acb_mat_entry(Cj, p, q));
        acb_mat_conjugate_transpose(Cmd, Cmid);
        acb_mat_add(Cmid, Cmid, Cmd, prec);
        acb_mat_scalar_mul_2exp_si(Cmid, Cmid, -1);
        arb_t R, mx, lb;
        arb_init(R);
        arb_init(mx);
        arb_init(lb);
        acb_mat_sub(resid, Cj, Cmid, prec);
        aic_mat_frobenius_norm(R, resid, prec);
        acb_mat_neg(negCmid, Cmid);
        aic_mat_herm_max_eig(mx, negCmid, prec);
        arb_neg(lb, mx);
        arb_sub(lb, lb, R, prec);
        double mineig_lb = dd(lb);
        if (mineig_lb < wm) wm = mineig_lb;
        if (mineig_lb < -tol) all_cp = 0;
        arb_clear(lb); arb_clear(mx); arb_clear(R);
        acb_mat_clear(resid); acb_mat_clear(negCmid); acb_mat_clear(Cmd);
        acb_mat_clear(Cmid); acb_mat_clear(UE); acb_mat_clear(Epq); acb_mat_clear(Cj);
    }
    *worst_mineig = wm;
    return all_cp;
}

/* ===================== T5: eta=0 oracle for Upsilon ========================= *
 * On EXACT-idempotent channels (the cleanest ground truth, Rule 6). Tests the
 * per-block design across >= 2 blocks (block_cond_exp(4,2): B = M_2 (+) M_2 — the
 * F2 multi-block coverage). The lem_RC structural identity R_j = 1_{L_j} (x) C_j
 * certifies the D1 trace direction + centrality + W_j stacking at once; the map
 * identities Upsilon Delta = 1_B and Delta Upsilon = Phi are the gauge-invariant
 * cross-checks against th_idemp_structure's Gamma o C_M decode. */
static void t5_one(const char *name, aic_ucp_kraus *phi)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);

    aic_factorize F;
    aic_factorize_build(&F, &v, &ae, phi, CPREC);
    aic_factorize_delta_build(&F, CPREC);
    aic_factorize_upsilon_build(&F, 0.3, CPREC);  /* tol_sigma generous (eta=0: =1) */

    double smax_dev;
    double rc = Rj_minus_1xCj_worst(&F, &smax_dev, CPREC);
    double ud = ups_del_defect(&F, CPREC);
    double du = del_ups_defect(&F, phi, /*phitilde=*/0, CPREC);
    double ui = upsI_minus_IB(&F, CPREC);
    double cp_mineig;
    int cp = upsilon_prime_is_cp(&F, 1e-9, &cp_mineig, CPREC);

    printf("T5 %s: N=%ld n_B=%ld m=%ld | ||R_j-1(x)C_j||=%.3e |sigma_max(C_j)-1|=%.3e "
           "| UpsDel-1_B=%.3e DelUps-Phi=%.3e | Ups(I)-I_B=%.3e Ups' CP=%d "
           "minEig=%.3e\n", name, (long) F.N, (long) F.n_B,
           (long) F.v->B->num_blocks, rc, smax_dev, ud, du, ui, cp, cp_mineig);

    AIC_CHECK_MSG(rc < 1e-9, "T5 %s: ||R_j - 1(x)C_j|| = %.3e not ~0 (lem_RC(i) / "
                  "D1 trace dir / centrality broken)", name, rc);
    AIC_CHECK_MSG(smax_dev < 1e-9, "T5 %s: |sigma_max(C_j)-1| = %.3e not ~0 "
                  "(lem_RC ||C_j||=1 at eta=0)", name, smax_dev);
    AIC_CHECK_MSG(ud < 1e-9, "T5 %s: Upsilon Delta - 1_B = %.3e not ~0", name, ud);
    AIC_CHECK_MSG(du < 1e-9, "T5 %s: Delta Upsilon - Phi = %.3e not ~0", name, du);
    AIC_CHECK_MSG(ui < 1e-9, "T5 %s: Upsilon(I_H) - I_B = %.3e not ~0 (not unital)",
                  name, ui);
    AIC_CHECK_MSG(cp == 1, "T5 %s: Upsilon' not CP (minEig=%.3e)", name, cp_mineig);

    aic_factorize_clear(&F);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t5_eta0_oracle(void)
{
    printf("T5 eta=0 oracle (Upsilon via lem_RC; UpsDel=1_B, DelUps=Phi):\n");
    aic_ucp_kraus phi;

    make_block_cond_exp(&phi, 4, 2);          /* M_2 (+) M_2, m=2 (multi-block) */
    t5_one("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    make_noiseless_subsystem(&phi, 2, 2);     /* M_2, m=1 */
    t5_one("noiseless_subsystem(2,2)", &phi);
    aic_ucp_kraus_clear(&phi);
}

/* max over M_N matrix units E_pq of ||Upsilon(E_pq) - Upsilon~(E_pq)||_op
 * (.tex:2899, ||Upsilon - Upsilon~||_cb <= O(eta)). Upsilon~ = F1's tilde decode. */
static double ups_minus_upstilde(const aic_factorize *F, slong prec)
{
    slong N = F->N, nB = F->n_B;
    acb_mat_t Epq, U, Ut;
    acb_mat_init(Epq, N, N);
    acb_mat_init(U, nB, nB);
    acb_mat_init(Ut, nB, nB);
    double worst = 0.0;
    for (slong p = 0; p < N; p++)
        for (slong q = 0; q < N; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_upsilon(U, F, Epq, prec);          /* UCP Upsilon */
            aic_factorize_upsilon_tilde(Ut, F, Epq, prec);   /* tilde Upsilon~ */
            double d = opdiff(U, Ut, prec);
            if (d > worst) worst = d;
        }
    acb_mat_clear(Ut); acb_mat_clear(U); acb_mat_clear(Epq);
    return worst;
}

/* ===================== T6: eta>0 (the headline + the D5 pin) ================ *
 * On make_mixconj (genuinely oblique eta>0, single-block B = M_2, r > 1 — the
 * regime that EXERCISES the D5 F-ordering). Headlines:
 *   - Upsilon UCP: Upsilon(I_H)=I_B (machine), Upsilon' CP (§C5 midpoint+Weyl);
 *   - sigma_max(C_j) >= 1 - O(eta) (lem_RC(ii), report);
 *   - ||Upsilon - Upsilon~|| = O(eta) (report c = defect/eta);
 *   - the D5 F-ordering pin: F-LEFT reconstructs Upsilon'Delta ~ 1_B (O(eta));
 *     F-RIGHT does NOT (O(1)) — the decisive r>1 distinguisher. */
static void test_t6_eta_gt0(void)
{
    printf("T6 eta>0 (Upsilon UCP + ||Ups-Ups~||=O(eta) + the D5 F-ordering pin):\n");
    struct { slong d, m; double t; } fx[] = {
        {4, 2, 0.03}, {5, 2, 0.03}};
    const int n_fx = (int) (sizeof(fx) / sizeof(fx[0]));
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

        aic_factorize F;
        aic_factorize_build(&F, &v, &ae, &phi, CPREC);
        aic_factorize_delta_build(&F, CPREC);
        aic_factorize_upsilon_build(&F, 0.5, CPREC);  /* tol_sigma generous for O(eta) */

        double smax_dev;
        double rc = Rj_minus_1xCj_worst(&F, &smax_dev, CPREC);
        double ui = upsI_minus_IB(&F, CPREC);
        double cp_mineig;
        int cp = upsilon_prime_is_cp(&F, 1e-9, &cp_mineig, CPREC);
        double uut = ups_minus_upstilde(&F, CPREC);
        double d5_left  = aic_factorize_upsilon_d5_pin(&F, /*f_left=*/1, CPREC);
        double d5_right = aic_factorize_upsilon_d5_pin(&F, /*f_left=*/0, CPREC);

        printf("  mixconj(%ld,%ld,%.2f): N=%ld n_B=%ld r=%ld eta=%.3e | "
               "||R_j-1(x)C_j||=%.3e (c=%.2f) sigma_max(C_j) in [1-%.3e,..] | "
               "Ups(I)-I_B=%.3e Ups' CP=%d minEig=%.3e | ||Ups-Ups~||=%.3e (c=%.2f) "
               "| D5: F-LEFT Ups'Del-1_B=%.3e F-RIGHT=%.3e\n", (long) fx[i].d,
               (long) fx[i].m, fx[i].t, (long) F.N, (long) F.n_B, (long) F.r, eta,
               rc, rc / eta, smax_dev, ui, cp, cp_mineig, uut, uut / eta,
               d5_left, d5_right);

        AIC_CHECK_MSG(rc < 8.0 * eta,
                      "T6 mixconj(%ld,%ld): ||R_j-1(x)C_j||=%.3e not O(eta) "
                      "(lem_RC(i) at eta>0)", (long) fx[i].d, (long) fx[i].m, rc);
        AIC_CHECK_MSG(ui < 1e-9, "T6 mixconj(%ld,%ld): Upsilon(I_H)-I_B=%.3e not ~0",
                      (long) fx[i].d, (long) fx[i].m, ui);
        AIC_CHECK_MSG(cp == 1, "T6 mixconj(%ld,%ld): Upsilon' NOT CP (minEig=%.3e)",
                      (long) fx[i].d, (long) fx[i].m, cp_mineig);
        AIC_CHECK_MSG(smax_dev < 5.0 * eta,
                      "T6 mixconj(%ld,%ld): |sigma_max(C_j)-1|=%.3e not O(eta)",
                      (long) fx[i].d, (long) fx[i].m, smax_dev);
        AIC_CHECK_MSG(uut < 8.0 * eta,
                      "T6 mixconj(%ld,%ld): ||Upsilon-Upsilon~||=%.3e not O(eta)",
                      (long) fx[i].d, (long) fx[i].m, uut);
        /* D5: r>1 here, so F-LEFT must reconstruct Ups'Del~1_B (O(eta)) and
         * F-RIGHT must NOT (clearly larger). If F-RIGHT is NOT clearly larger,
         * the F-ordering is not what we think — a FINDING (fail loud). */
        AIC_CHECK_MSG(F.r > 1, "T6 mixconj(%ld,%ld): r=%ld not >1 — D5 pin VACUOUS "
                      "(F is 1-dim, can't distinguish orderings)",
                      (long) fx[i].d, (long) fx[i].m, (long) F.r);
        AIC_CHECK_MSG(d5_left < 8.0 * eta,
                      "T6 mixconj(%ld,%ld): D5 F-LEFT Ups'Del-1_B=%.3e not O(eta)",
                      (long) fx[i].d, (long) fx[i].m, d5_left);
        AIC_CHECK_MSG(d5_right > 10.0 * d5_left + 1e-3,
                      "T6 mixconj(%ld,%ld): D5 F-RIGHT (%.3e) NOT clearly worse than "
                      "F-LEFT (%.3e) — F-ordering is not as believed (FINDING)",
                      (long) fx[i].d, (long) fx[i].m, d5_right, d5_left);

        aic_factorize_clear(&F);
        arb_clear(iso);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* Forward decl: the multi-block eta>0 fixture (defined near T8). */
static void make_mixconj_blocks_f4(aic_ucp_kraus *out, slong d, slong m, double t,
                                   slong prec);

/* ===================== T6b: lem_RC at m>=2 AND eta>0 (§C13(c) debt) ========= *
 * The FIRST exercise of the full F3 decode at MULTI-BLOCK (>=2 equivalence classes)
 * AND eta>0 — previously BLOCKED on the strict CP-extraction gate, now LIVE via the
 * tolerance-aware PSD-cone extraction (FINDINGS §C14). Repays the §C13(c) coverage
 * debt: the lem_RC teeth (R_j = 1_{L_j} (x) C_j EXACTLY at eta>0, since the per-block
 * Heisenberg-Weyl Pauli set is an EXACT 1-design; sigma_max(C_j) >= 1-O(eta)), the
 * UCP-to-O(eta^2) properties, and the round-trip (Upsilon Delta ~ 1_B, Delta Upsilon
 * ~ Phi~ to O(eta)). MUTATION-PROVEN by hand (Rule 7; RED observed then restored):
 * the ||R_j-1(x)C_j|| tooth FIRES at m>=2 eta>0, where it was VACUOUS at eta=0:
 *   - drop one Pauli in aic_factorize_upsilon_Rj (`if (sa==1&&sb==0) continue;`):
 *     ||R_j-1(x)C_j|| 8.6e-78 -> 7.5e-3 (decisive RED, >> the 1e-9 gate);
 *   - U->I (`acb_mat_one(U)` after aic_dhom_pauli, kills the 1-design averaging):
 *     ||R_j-1(x)C_j|| -> 3.0e-2 (RED).
 *   At eta=0 the SAME drop-Pauli mutation leaves the rc check GREEN (~1e-77; the
 *   smax check fires instead) — confirming the rc tooth was previously vacuous and
 *   is now LIVE only at m>=2 eta>0. */
static void t6b_one(const char *name, aic_ucp_kraus *phi)
{
    double eta = eta_proxy(phi, CPREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, eta, CPREC);   /* §C11: pass eta, not eps_assoc */

    aic_factorize F;
    aic_factorize_build(&F, &v, &ae, phi, CPREC);
    aic_factorize_delta_build(&F, CPREC);
    aic_factorize_upsilon_build(&F, 0.5, CPREC);       /* tol_sigma generous for O(eta) */

    double smax_dev;
    double rc = Rj_minus_1xCj_worst(&F, &smax_dev, CPREC);
    double ui = upsI_minus_IB(&F, CPREC);
    double cp_mineig;
    int cp = upsilon_prime_is_cp(&F, 1e-9, &cp_mineig, CPREC);
    double ud = ups_del_defect(&F, CPREC);                       /* Upsilon Delta - 1_B */
    double du = del_ups_defect(&F, phi, /*use_phitilde=*/1, CPREC); /* Delta Upsilon - Phi~ */

    printf("  %s: N=%ld n_B=%ld m=%ld r=%ld eta=%.3e | ||R_j-1(x)C_j||=%.3e "
           "|sigma_max(C_j)-1|=%.3e (c=%.2f) | Ups(I)-I_B=%.3e Ups' CP=%d minEig=%.3e "
           "| UpsDel-1_B=%.3e (c=%.2f) DelUps-Phi~=%.3e (c=%.2f)\n", name, (long) F.N,
           (long) F.n_B, (long) F.v->B->num_blocks, (long) F.r, eta, rc, smax_dev,
           smax_dev / eta, ui, cp, cp_mineig, ud, ud / eta, du, du / eta);

    /* lem_RC(i): R_j = 1_{L_j} (x) C_j is EXACT at eta>0 (the Pauli set is an exact
     * 1-design — §C13(c)); a drop-Pauli mutation drives this to O(1) (mutation-proven
     * by hand). It was VACUOUS at eta=0. */
    AIC_CHECK_MSG(rc < 1e-9, "T6b %s: ||R_j-1(x)C_j||=%.3e not ~0 (lem_RC(i) — the "
                  "Pauli design is NOT an exact 1-design at eta>0)", name, rc);
    /* lem_RC(ii): sigma_max(C_j) >= 1 - O(eta). */
    AIC_CHECK_MSG(smax_dev < 5.0 * eta, "T6b %s: |sigma_max(C_j)-1|=%.3e not O(eta) "
                  "(lem_RC(ii))", name, smax_dev);
    /* Upsilon UCP: unital to machine, Upsilon' CP (after the PSD-cone projection). */
    AIC_CHECK_MSG(ui < 1e-9, "T6b %s: Upsilon(I_H)-I_B=%.3e not ~0 (not unital)",
                  name, ui);
    AIC_CHECK_MSG(cp == 1, "T6b %s: Upsilon' NOT CP (minEig=%.3e)", name, cp_mineig);
    /* Round-trip to O(eta) (the headline th_factorization defects at m>=2 eta>0). */
    AIC_CHECK_MSG(ud < 20.0 * eta, "T6b %s: Upsilon Delta - 1_B=%.3e not O(eta)",
                  name, ud);
    AIC_CHECK_MSG(du < 20.0 * eta, "T6b %s: Delta Upsilon - Phi~=%.3e not O(eta)",
                  name, du);

    aic_factorize_clear(&F);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t6b_multiblock_eta_gt0(void)
{
    printf("T6b m>=2 ∧ eta>0 (§C13(c) debt — lem_RC + UCP + round-trip on the FIRST "
           "multi-block eta>0 fixture; LIVE via the FINDINGS §C14 PSD-cone extraction):\n");
    aic_ucp_kraus phi;

    make_mixconj_blocks_f4(&phi, 4, 2, 0.03, CPREC);   /* (M_2 (+) M_2) mixed-conj */
    t6b_one("mixconj_blocks(4,2,0.03)", &phi);
    aic_ucp_kraus_clear(&phi);

    make_mixconj_blocks_f4(&phi, 6, 3, 0.03, CPREC);   /* (M_3 (+) M_3) mixed-conj */
    t6b_one("mixconj_blocks(6,3,0.03)", &phi);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== F4.1 helpers (end-to-end verify + duals) ============= *
 * The composed-map Choi (Delta Upsilon, Upsilon Delta - 1_B), the eig-free
 * certified upper bound, and the dual channels Dec=Delta*, Enc=Upsilon*. See
 * include/aic_factorize.h (the F4 banner), src/aic_factorize_verify.c +
 * src/aic_factorize_dual.c, docs/research/factorize_f4_design.md §A/§B/§D. */

/* Frobenius norm of an acb_mat as a double (midpoint route, no false-fail). */
static double frob(const acb_mat_t A, slong prec)
{
    arb_t f;
    arb_init(f);
    aic_mat_frobenius_norm(f, A, prec);
    double o = dd(f);
    arb_clear(f);
    return o;
}

/* ||Sum_a K_a^dag K_a - 1_{dim_H}||_op (the dual TP defect: Dec/Enc trace-
 * preserving <=> Delta/Upsilon unital, since Sum K_a^dag K_a = Phi(1) for the
 * Heisenberg map Phi that the duals are read off from). */
static double kraus_unital_defect(const aic_ucp_kraus *phi, slong prec)
{
    slong dH = phi->dim_H;
    acb_mat_t G, Kd, KdK, I;
    acb_mat_init(G, dH, dH);
    acb_mat_init(I, dH, dH);
    acb_mat_zero(G);
    for (slong a = 0; a < phi->r; a++) {
        acb_mat_init(Kd, dH, phi->dim_K);
        acb_mat_init(KdK, dH, dH);
        acb_mat_conjugate_transpose(Kd, phi->K[a]);   /* K_a^dag : dim_K -> dim_H */
        acb_mat_mul(KdK, Kd, phi->K[a], prec);         /* K_a^dag K_a : dim_H x dim_H */
        acb_mat_add(G, G, KdK, prec);
        acb_mat_clear(KdK); acb_mat_clear(Kd);
    }
    acb_mat_one(I);
    acb_mat_sub(G, G, I, prec);
    double o;
    {
        acb_mat_t M;
        acb_mat_init(M, dH, dH);
        for (slong i = 0; i < dH; i++)
            for (slong j = 0; j < dH; j++)
                acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(G, i, j));
        arb_t e;
        arb_init(e);
        aic_mat_opnorm(e, M, prec);
        o = dd(e);
        arb_clear(e);
        acb_mat_clear(M);
    }
    acb_mat_clear(I); acb_mat_clear(G);
    return o;
}

/* Worst ||block||_F over the OFF-block-diagonal INPUT columns of an ambient route-(i)
 * J_UpsDel (the §D probe; defined below, used by T7's multi-block oracle too). */
static double upsdel_offblock_input_mass(const acb_mat_t J, const aic_factorize *F);

/* ===================== T7: eta=0 ORACLE (the F4 spine) ====================== *
 * On EXACT-idempotent channels (Delta Upsilon = Phi and Upsilon Delta = 1_B as
 * MAPS, by the eta=0 reduction of the whole pipeline). The strongest single
 * check: it pins all three CP-izations + both unitalizations + the v-inversion to
 * exact reduction. Assertions:
 *   - ||J_DelUps||_F < 1e-10 (Choi(Delta Upsilon) = Choi(Phi) to machine).
 *   - ||J_UpsDel||_F < 1e-10 (Choi(Upsilon Delta - 1_B) = 0).
 *   - eig-free hi_du, hi_ud < 1e-9 (the 2||J||_F bracket collapses at J=0).
 *   - dual TP-ness: Dec, Enc trace-preserving <=> Delta, Upsilon unital.
 * MUTATION (recorded, mutation-proven by hand): dropping the Choi(Phi) subtraction
 * leaves J_DelUps = Choi(Delta Upsilon), a UCP Choi with ||.||_F = O(1) -> RED. */
static void t7_one(const char *name, aic_ucp_kraus *phi)
{
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    aic_cstar_build(&B, &v, iso, &ae.A, 0.0, CPREC);

    aic_factorize F;
    aic_factorize_build(&F, &v, &ae, phi, CPREC);
    aic_factorize_delta_build(&F, CPREC);
    aic_factorize_upsilon_build(&F, 0.3, CPREC);

    slong N = F.N, nB = F.n_B;
    acb_mat_t Jdu, Jud;
    acb_mat_init(Jdu, N * N, N * N);
    acb_mat_init(Jud, nB * nB, nB * nB);
    aic_factorize_choi_delups(Jdu, &F, CPREC);
    aic_factorize_choi_upsdel(Jud, &F, CPREC);
    double nf_du = frob(Jdu, CPREC), nf_ud = frob(Jud, CPREC);

    arb_t lo, hi_du, hi_ud;
    arb_init(lo);
    arb_init(hi_du);
    arb_init(hi_ud);
    aic_factorize_eigfree_delups(lo, hi_du, &F, CPREC);
    aic_factorize_eigfree_upsdel(lo, hi_ud, &F, CPREC);

    /* §D route-(i)=in-B structural check at multi-block: the OFF-block-diagonal
     * INPUT columns of the ambient J_UpsDel must be ZERO (Delta drops off-block
     * input). Trivially 0 single-block (no off-block pairs); the teeth live at
     * m>=2. (At eta>0 multi-block this is the LOAD-BEARING §D probe, currently
     * BLOCKED on the F3 multi-block CP-extraction wall — see T8 + the report.) */
    double offblk = upsdel_offblock_input_mass(Jud, &F);

    /* Dual channels Dec=Delta*, Enc=Upsilon*; TP <=> Delta, Upsilon unital. */
    aic_ucp_kraus dec, enc;
    aic_factorize_dec_kraus(&dec, &F, CPREC);
    aic_factorize_enc_kraus(&enc, &F, CPREC);
    double tp_dec = kraus_unital_defect(&dec, CPREC);   /* Sum D_a^dag D_a - 1_H */
    double tp_enc = kraus_unital_defect(&enc, CPREC);   /* Sum E_a^dag E_a - 1_B */

    printf("T7 %s: N=%ld n_B=%ld m=%ld | ||J_DelUps||_F=%.3e ||J_UpsDel||_F=%.3e "
           "| hi_du=%.3e hi_ud=%.3e | offblk=%.3e | Dec TP=%.3e (r=%ld) Enc TP=%.3e "
           "(r=%ld)\n", name, (long) N, (long) nB, (long) F.v->B->num_blocks, nf_du,
           nf_ud, dd(hi_du), dd(hi_ud), offblk, tp_dec, (long) dec.r, tp_enc,
           (long) enc.r);

    AIC_CHECK_MSG(nf_du < 1e-10, "T7 %s: ||J_DelUps||_F = %.3e not ~0 (Delta Upsilon "
                  "!= Phi at eta=0)", name, nf_du);
    AIC_CHECK_MSG(nf_ud < 1e-10, "T7 %s: ||J_UpsDel||_F = %.3e not ~0 (Upsilon Delta "
                  "!= 1_B at eta=0)", name, nf_ud);
    AIC_CHECK_MSG(dd(hi_du) < 1e-9, "T7 %s: eig-free hi_du = %.3e not ~0", name,
                  dd(hi_du));
    AIC_CHECK_MSG(dd(hi_ud) < 1e-9, "T7 %s: eig-free hi_ud = %.3e not ~0", name,
                  dd(hi_ud));
    AIC_CHECK_MSG(offblk < 1e-9, "T7 %s: off-block-diagonal input mass = %.3e not ~0 "
                  "(route-(i) ambient J_UpsDel does NOT zero off-block columns — the "
                  "§A.2 ambient=in-B convention is wrong)", name, offblk);
    AIC_CHECK_MSG(tp_dec < 1e-9, "T7 %s: Dec NOT trace-preserving (=%.3e) — Delta "
                  "not unital", name, tp_dec);
    AIC_CHECK_MSG(tp_enc < 1e-9, "T7 %s: Enc NOT trace-preserving (=%.3e) — Upsilon "
                  "not unital", name, tp_enc);

    aic_ucp_kraus_clear(&enc);
    aic_ucp_kraus_clear(&dec);
    arb_clear(hi_ud); arb_clear(hi_du); arb_clear(lo);
    acb_mat_clear(Jud); acb_mat_clear(Jdu);
    aic_factorize_clear(&F);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t7_eta0_oracle(void)
{
    printf("T7 eta=0 oracle (J_DelUps=J_UpsDel=0; duals TP) — the F4 spine:\n");
    aic_ucp_kraus phi;

    make_block_cond_exp(&phi, 4, 2);          /* M_2 (+) M_2, m=2 (multi-block) */
    t7_one("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    make_block_cond_exp(&phi, 6, 3);          /* M_3 (+) M_3, m=2 (multi-block) */
    t7_one("block_cond_exp(6,3)", &phi);
    aic_ucp_kraus_clear(&phi);

    make_noiseless_subsystem(&phi, 2, 2);     /* M_2, m=1 */
    t7_one("noiseless_subsystem(2,2)", &phi);
    aic_ucp_kraus_clear(&phi);
}

/* Multi-block eta>0 fixture (replicated from test_cstar_build.c make_mixconj_blocks,
 * which is static there): block_cond_exp(d,m) [M_m (+) M_{d-m}, TWO genuine blocks]
 * mixed with its fixed-unitary conjugate — the ONLY suite fixture with >=2
 * equivalence classes at eta>0 (the §C13(c) m>=2 ∧ eta>0 coverage debt; make_*
 * helpers from test_idemp.h). Pass `eta` (not the ~700x smaller eps_assoc) as the
 * build's eps per §C11. */
static void make_mixconj_blocks_f4(aic_ucp_kraus *out, slong d, slong m, double t,
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

/* THE §D PROBE, teeth-bearing form (design §D "OR equivalently assert the
 * off-block-diagonal input columns of J^(i) are zero"). The ambient route-(i)
 * Choi J_UpsDel (n_B^2 x n_B^2) is built over ALL n_B^2 ambient units; ambient =
 * in-B holds IFF Delta drops the off-block-diagonal input coordinates, i.e. IFF
 * the (p,q) OUTPUT BLOCK of J_UpsDel — the entries J[p*n_B+a, q*n_B+b], a,b<n_B —
 * is ZERO whenever the ambient unit E_pq is OFF-block-diagonal (p,q in different
 * B-blocks). For such (p,q), Upsilon Delta(E_pq) = 0 AND the "-1_B(E_pq)" term is
 * 0 (E_pq has no in-B diagonal mass), so the whole block must be ~0. Returns the
 * worst ||block||_F over all off-block-diagonal (p,q). This is the LOAD-BEARING
 * route-(i)=route-(ii) check (eta=0 is BLIND to it: Upsilon Delta=1_B exactly):
 * if it is large, Delta reads off-block-diagonal input and route (ii) (the in-B
 * Choi) would be the correct cb-norm — a BLOCKER. A route-(iii) per-block-max
 * would MISS exactly this cross-block leakage. */
static double upsdel_offblock_input_mass(const acb_mat_t J, const aic_factorize *F)
{
    const aic_dhom_B *B = F->v->B;
    slong nB = F->n_B;
    /* block_of[i] = the B-block containing ambient row/col index i. */
    slong *block_of = (slong *) malloc((size_t) nB * sizeof(slong));
    slong off = 0;
    for (slong l = 0; l < B->num_blocks; l++)
        for (slong a = 0; a < B->d[l]; a++) block_of[off++] = l;

    double worst = 0.0;
    for (slong p = 0; p < nB; p++)
        for (slong q = 0; q < nB; q++) {
            if (block_of[p] == block_of[q]) continue; /* in-block input: skip */
            double s = 0.0;
            for (slong a = 0; a < nB; a++)
                for (slong b = 0; b < nB; b++) {
                    acb_srcptr z = acb_mat_entry(J, p * nB + a, q * nB + b);
                    double re = dd(acb_realref(z)), im = dd(acb_imagref(z));
                    s += re * re + im * im;
                }
            s = sqrt(s);
            if (s > worst) worst = s;
        }
    free(block_of);
    return worst;
}

/* ===================== T8: eta>0 per-instance eig-free + route probe ======== *
 * On make_mixconj (single-block) AND make_mixconj_blocks (multi-block, the
 * §C13(c) debt). Per-instance RIGOROUS certified O(eta) upper bound (eig-free,
 * Julia-free):
 *   - hi_du <= C_max*eta and hi_ud <= C_max*eta (C_max=40 generous: the eig-free
 *     2N looseness inflates the constant; a COARSE O(eta)-SCALING gate, NOT tight).
 * THE §D PROBE (route (i) ambient ≡ route (ii) in-B; the one empirical check, eta=0
 * BLIND so it MUST live at eta>0): assert the OFF-block-diagonal INPUT columns of
 * the ambient route-(i) J_UpsDel are ZERO (upsdel_offblock_input_mass < 1e-9). This
 * is WHY ambient = in-B: Delta drops off-block-diagonal input, so feeding it a non-
 * block-diagonal M_{n_B} unit contributes nothing. If it were nonzero, route (ii)
 * (the in-B Choi) would be the correct cb-norm and the §A.2 convention is WRONG
 * (a BLOCKER — fail loud, STOP-and-report). A route-(iii) per-block-max would MISS
 * exactly this cross-block leakage (it sees only in-block inputs), so i != iii.
 *
 * THE FINDINGS §C14 PSD-CONE GUARDS (mutation-proven by hand, Rule 7; RED observed
 * then restored). The multi-block eta>0 build runs ONLY because the tolerance-aware
 * extraction projects the O(eta^2) cone defect away. Its two genuine-bug guards have
 * teeth:
 *   - upsilon_build's aggregate clip-mass ceiling (AIC_FACTORIZE_CONE_MASS_CEIL=1e-2):
 *     lowering it to 1e-7 (< the measured 8.1e-6 clipped mass at eta=2.3e-2) trips the
 *     assert "clipped PSD-cone-defect mass exceeds ceiling" -> abort (RED);
 *   - the per-eigenvalue abort in aic_ucp_choi_to_kraus_latd_tol (neg_tol=1e-3):
 *     injecting an O(1) negative (subtract 0.5*I from Cj_choi) still fails loud
 *     "eigenvalue -1 <= -neg_tol=-0.001; C is not PSD" -> abort (RED). So the
 *     PSD-cone projection admits ONLY the O(eta^2) cone defect, not a genuine non-CP. */
static void t8_one(const char *name, aic_ucp_kraus *phi, double C_max,
                   int do_route_probe)
{
    double eta = eta_proxy(phi, CPREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, phi, CPREC);
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
    arb_init(iso);
    /* §C11: pass eta (not the ~700x smaller eps_assoc) as the build's eps scale. */
    aic_cstar_build(&B, &v, iso, &ae.A, eta, CPREC);

    aic_factorize F;
    aic_factorize_build(&F, &v, &ae, phi, CPREC);
    aic_factorize_delta_build(&F, CPREC);
    aic_factorize_upsilon_build(&F, 0.5, CPREC);

    slong N = F.N, nB = F.n_B;
    arb_t lo, hi_du, hi_ud, hi_ud_ii;
    arb_init(lo);
    arb_init(hi_du);
    arb_init(hi_ud);
    arb_init(hi_ud_ii);
    aic_factorize_eigfree_delups(lo, hi_du, &F, CPREC);

    /* Build J_UpsDel (route (i), ambient n_B^2) ONCE: use it for the eig-free hi_ud,
     * the §D off-block-input probe, AND (multi-block) the route-(ii) comparison. */
    acb_mat_t Jud;
    acb_mat_init(Jud, nB * nB, nB * nB);
    aic_factorize_choi_upsdel(Jud, &F, CPREC);
    aic_cbnorm_eigfree_ball_choi(lo, hi_ud, Jud, nB, CPREC);
    double rdu = dd(hi_du) / eta, rud = dd(hi_ud) / eta;

    double offblk = do_route_probe ? upsdel_offblock_input_mass(Jud, &F) : -1.0;

    /* ROUTE (i) vs ROUTE (ii) (the §D check; multi-block only). Route (ii) = the
     * genuinely-IN-B Choi: zero the OFF-block-diagonal INPUT columns of the ambient
     * route-(i) J_UpsDel (the rows p*nB+a, cols q*nB+b with block_of[p]!=block_of[q])
     * and recompute the eig-free hi. Since Delta drops off-block input (offblk ~ 0),
     * route (ii) ≡ route (i): hi_ud_ii must AGREE with hi_ud. If they diverged, the
     * ambient cb-norm would NOT equal the in-B cb-norm (the §A.2 convention wrong). */
    double rud_agree = -1.0;
    if (do_route_probe) {
        slong *blk = (slong *) malloc((size_t) nB * sizeof(slong));
        slong off = 0;
        for (slong l = 0; l < F.v->B->num_blocks; l++)
            for (slong a = 0; a < F.v->B->d[l]; a++) blk[off++] = l;
        acb_mat_t Jii;
        acb_mat_init(Jii, nB * nB, nB * nB);
        acb_mat_set(Jii, Jud);
        for (slong p = 0; p < nB; p++)
            for (slong q = 0; q < nB; q++)
                if (blk[p] != blk[q])               /* off-block INPUT column -> 0 */
                    for (slong a = 0; a < nB; a++)
                        for (slong b = 0; b < nB; b++)
                            acb_zero(acb_mat_entry(Jii, p * nB + a, q * nB + b));
        aic_cbnorm_eigfree_ball_choi(lo, hi_ud_ii, Jii, nB, CPREC);
        rud_agree = dd(hi_ud_ii) - dd(hi_ud);
        if (rud_agree < 0) rud_agree = -rud_agree;
        acb_mat_clear(Jii);
        free(blk);
    }
    acb_mat_clear(Jud);
    (void) N;

    printf("  %s: N=%ld n_B=%ld m=%ld eta=%.3e | hi_du=%.3e (hi_du/eta=%.2f) "
           "hi_ud=%.3e (hi_ud/eta=%.2f)", name, (long) F.N, (long) nB,
           (long) F.v->B->num_blocks, eta, dd(hi_du), rdu, dd(hi_ud), rud);
    if (do_route_probe)
        printf(" | offblk-input mass=%.3e | route(i)-vs-(ii) hi |diff|=%.3e "
               "(hi_ud_ii=%.3e)", offblk, rud_agree, dd(hi_ud_ii));
    printf("\n");

    AIC_CHECK_MSG(rdu < C_max, "T8 %s: hi_du/eta = %.2f >= C_max=%.0f (not O(eta))",
                  name, rdu, C_max);
    AIC_CHECK_MSG(rud < C_max, "T8 %s: hi_ud/eta = %.2f >= C_max=%.0f (not O(eta))",
                  name, rud, C_max);
    if (do_route_probe) {
        AIC_CHECK_MSG(offblk < 1e-9, "T8 %s: off-block-diagonal input mass of route-(i) "
                      "J_UpsDel = %.3e not ~0 — ambient != in-B, route (ii) is the "
                      "correct cb-norm, the §A.2 convention is WRONG (BLOCKER)", name,
                      offblk);
        AIC_CHECK_MSG(rud_agree < 1e-9, "T8 %s: route(i)-vs-(ii) hi disagree by %.3e — "
                      "ambient cb-norm != in-B cb-norm (route i != ii)", name, rud_agree);
    }

    arb_clear(hi_ud_ii); arb_clear(hi_ud); arb_clear(hi_du); arb_clear(lo);
    aic_factorize_clear(&F);
    arb_clear(iso);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
}

static void test_t8_eta_gt0(void)
{
    printf("T8 eta>0 per-instance eig-free O(eta) bound + the route-(i)-vs-(ii) probe "
           "(single-block: probe degenerate; MULTI-BLOCK: §C13(c) debt now LIVE via the "
           "tolerance-aware PSD-cone extraction, FINDINGS §C14):\n");
    aic_ucp_kraus phi;

    /* single-block (B = M_2): the §D route probe is DEGENERATE here (n_B=2, all
     * indices in one block, no off-block pairs), so the per-instance O(eta) bound
     * is the live check. The multi-block teeth need the fixture below. */
    make_mixconj(&phi, 4, 2, 0.03, CPREC);
    t8_one("mixconj(4,2,0.03)", &phi, 40.0, /*probe=*/0);
    aic_ucp_kraus_clear(&phi);

    make_mixconj(&phi, 5, 2, 0.03, CPREC);
    t8_one("mixconj(5,2,0.03)", &phi, 40.0, /*probe=*/0);
    aic_ucp_kraus_clear(&phi);

    /* MULTI-BLOCK eta>0 (§C13(c) m>=2 ∧ eta>0 debt; make_mixconj_blocks_f4 — B has
     * >=2 equivalence classes). PREVIOUSLY BLOCKED: the full F3 build aborts on the
     * strict choi_to_kraus_latd PSD gate because the UCP Delta block-Choi has an
     * O(eta^2) negative eigenvalue (Delta' is manifestly CP only when Delta~=v is an
     * EXACT *-hom; at eta>0 multi-block v is an O(eta)-iso, so Delta' is CP only to
     * O(eta^2)). NOW LIVE: the tolerance-aware PSD-cone extraction (FINDINGS §C14)
     * projects that O(eta^2) defect away (keeping Delta within O(eta^2) of Delta~,
     * inside th_factorization's O(eta) tolerance), so upsilon_build COMPLETES and the
     * §D route-(i)-vs-(ii) probe + per-instance O(eta) bound run. C_max generous (the
     * eig-free 2N looseness inflates the constant). */
    make_mixconj_blocks_f4(&phi, 4, 2, 0.03, CPREC);
    t8_one("mixconj_blocks(4,2,0.03)", &phi, 60.0, /*probe=*/1);
    aic_ucp_kraus_clear(&phi);

    make_mixconj_blocks_f4(&phi, 6, 3, 0.03, CPREC);
    t8_one("mixconj_blocks(6,3,0.03)", &phi, 60.0, /*probe=*/1);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    /* F3 standalone GO/NO-GO pins (blueprint §4) — run FIRST. */
    test_ptrace_pin();        /* D1 partial-trace direction (#1 risk) */
    test_centrality_pin();    /* lem_RC centrality (diag_j2)          */

    test_t1_eta0_oracle();
    test_t2_eta_gt0();
    test_t3_eta0_oracle();
    test_t4_eta_gt0();
    test_t5_eta0_oracle();
    test_t6_eta_gt0();
    test_t6b_multiblock_eta_gt0(); /* §C13(c): lem_RC at m>=2 ∧ eta>0 (FINDINGS §C14) */
    test_t7_eta0_oracle();    /* F4.1 eta=0 spine */
    test_t8_eta_gt0();        /* F4.1 per-instance + route probe */
    aic_test_report("test_factorize");
    return 0;
}
