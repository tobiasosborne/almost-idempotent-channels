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
    (void) make_block_cond_exp;
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
    aic_factorize_build(&F, &v, &ae, CPREC);

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
        aic_factorize_build(&F, &v, &ae, CPREC);

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

int main(void)
{
    test_t1_eta0_oracle();
    test_t2_eta_gt0();
    aic_test_report("test_factorize");
    return 0;
}
