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

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
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

int main(void)
{
    test_t1_eta0_oracle();
    test_t2_eta_gt0();
    test_t3_eta0_oracle();
    test_t4_eta_gt0();
    aic_test_report("test_factorize");
    return 0;
}
