/* test_gamma_kraus.c — cross-checks for the prop_Gamma explicit Kraus form of the
 * conditional expectation Gamma : B(M) -> A (bead aic-ynu, I3; the eta=0 case).
 * Realizes prop_Gamma (approximate_algebras.tex:2106-2122):
 *   Gamma_j(X) = Tr_{E_j}(W_j X W_j^dag (1_{L_j} (x) gamma_j))   (eq Gamma, :2110)
 *             = L_j^dag (X (x) 1_{F_j}) L_j                       (Choi form, :2120)
 * aic_idemp_gamma_kraus extracts the density matrices gamma_j and the Kraus L_j of
 * the STORED d->Gamma and certifies (in production) that the prop_Gamma map MATCHES
 * d->Gamma. This test re-asserts the gates INDEPENDENTLY and reports the numbers.
 *
 * THE CRUX TOOTH (.tex:2113). The prop_Gamma Kraus map must EQUAL d->Gamma: apply
 * both to a basis E_ab of B(M) and compare the A-element embedded via the
 * *-monomorphism w (gauge-invariant), ||w(Gamma_kraus(E_ab)) - w(d->Gamma(E_ab))||
 * <= 1e-10 for all a,b. Pins the gamma_j extraction + the Tr_{E_j} convention + the
 * W_j usage. The MUTATION PROOF (perturb gamma_j) confirms the tooth bites.
 *
 * Plus: (a) each gamma_j a valid density matrix (Hermitian/PSD/trace 1);
 * (b) Gamma_kraus unital (Gamma_j(1_M) = 1_{L_j}); (c) Gamma.w = 1_A (apply
 * Gamma_kraus to w(B_k) and recover B_k, the .tex:2107 defining property).
 *
 * SCOPE: noiseless_subsystem(2,2) [E_1 = C^2, the NONTRIVIAL gamma_1 = (1/2)1 —
 * THE crux], block_cond_exp(5,2) [M_2 (+) M_3, E_j = C^1], identity(3) [E = C^1],
 * block_cond_exp(4,2) [M_2 (+) M_2, E_j = C^1].
 */
#include <math.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "test_idemp.h"

#define GKPREC 128
static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* wG(X) = w(d->Gamma(X)) in B(M): the stored conditional expectation embedded via
 * the *-monomorphism w. (Same gauge-safe target the production cert uses.) */
static void wG_stored(acb_mat_t out, const aic_idemp_decomp *d, const acb_mat_t X,
                      slong prec)
{
    slong m = d->dim_M, dA = d->dim_A;
    acb_mat_t vx, g, wg;
    acb_mat_init(vx, m * m, 1);
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            acb_set(acb_mat_entry(vx, i * m + j, 0), acb_mat_entry(X, i, j));
    acb_mat_init(g, dA, 1);
    acb_mat_mul(g, d->Gamma, vx, prec);
    acb_mat_init(wg, m * m, 1);
    acb_mat_mul(wg, d->w, g, prec);
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            acb_set(acb_mat_entry(out, i, j), acb_mat_entry(wg, i * m + j, 0));
    acb_mat_clear(wg);
    acb_mat_clear(g);
    acb_mat_clear(vx);
}

/* Gamma_j(X) via eq Gamma (Tr_{E_j} form) for the PUBLIC gamma_j (used in the
 * mutation proof, where gamma_j is replaced by a perturbed matrix). */
static void gamma_j_tr(acb_mat_t Gj, const acb_mat_t Wj, const acb_mat_t gamma,
                       const acb_mat_t X, slong dL, slong dE, slong m, slong prec)
{
    slong db = dL * dE;
    acb_mat_t Wjd, WX, WXW, oneL, Ig, Z;
    acb_mat_init(Wjd, m, db);
    acb_mat_conjugate_transpose(Wjd, Wj);
    acb_mat_init(WX, db, m);
    acb_mat_mul(WX, Wj, X, prec);
    acb_mat_init(WXW, db, db);
    acb_mat_mul(WXW, WX, Wjd, prec);
    acb_mat_init(oneL, dL, dL);
    acb_mat_one(oneL);
    acb_mat_init(Ig, db, db);
    aic_mat_kronecker(Ig, oneL, gamma, prec);
    acb_mat_init(Z, db, db);
    acb_mat_mul(Z, WXW, Ig, prec);
    aic_mat_partial_trace_right(Gj, Z, dL, dE, prec);
    acb_mat_clear(Z);
    acb_mat_clear(Ig);
    acb_mat_clear(oneL);
    acb_mat_clear(WXW);
    acb_mat_clear(WX);
    acb_mat_clear(Wjd);
}

/* wGk(X) = w(Gamma_kraus(X)) in B(M), re-embedded as sum_j W_j^dag (Gamma_j(X) (x)
 * 1_{E_j}) W_j with Gamma_j via the eq-Gamma Tr form (public gamma array). */
static void wGk_prop(acb_mat_t out, const aic_idemp_wedderburn *W,
                     const acb_mat_t *gamma, const acb_mat_t X, slong prec)
{
    slong m = W->dim_M;
    acb_mat_zero(out);
    for (slong j = 0; j < W->num_blocks; j++) {
        slong dL = W->dim_L[j], dE = W->dim_E[j], db = dL * dE;
        acb_mat_t Gj, oneE, GjXI, Wjd, t2, back;
        acb_mat_init(Gj, dL, dL);
        gamma_j_tr(Gj, W->W_j[j], gamma[j], X, dL, dE, m, prec);
        acb_mat_init(oneE, dE, dE);
        acb_mat_one(oneE);
        acb_mat_init(GjXI, db, db);
        aic_mat_kronecker(GjXI, Gj, oneE, prec);
        acb_mat_init(Wjd, m, db);
        acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
        acb_mat_init(t2, m, db);
        acb_mat_mul(t2, Wjd, GjXI, prec);
        acb_mat_init(back, m, m);
        acb_mat_mul(back, t2, W->W_j[j], prec);
        acb_mat_add(out, out, back, prec);
        acb_mat_clear(back);
        acb_mat_clear(t2);
        acb_mat_clear(Wjd);
        acb_mat_clear(GjXI);
        acb_mat_clear(oneE);
        acb_mat_clear(Gj);
    }
}

/* max over E_ab of ||w(Gamma_kraus(E_ab)) - w(d->Gamma(E_ab))||_op (the crux). */
static double match_resid(const aic_idemp_wedderburn *W, const acb_mat_t *gamma,
                          const aic_idemp_decomp *d, slong prec)
{
    slong m = d->dim_M;
    double worst = 0;
    acb_mat_t X, a, b, diff;
    arb_t e;
    acb_mat_init(X, m, m);
    acb_mat_init(a, m, m);
    acb_mat_init(b, m, m);
    acb_mat_init(diff, m, m);
    arb_init(e);
    for (slong ia = 0; ia < m; ia++)
        for (slong ib = 0; ib < m; ib++) {
            acb_mat_zero(X);
            acb_set_si(acb_mat_entry(X, ia, ib), 1);
            wG_stored(a, d, X, prec);
            wGk_prop(b, W, gamma, X, prec);
            acb_mat_sub(diff, a, b, prec);
            aic_mat_opnorm(e, diff, prec);
            double r = dd(e);
            if (r > worst) worst = r;
        }
    arb_clear(e);
    acb_mat_clear(diff);
    acb_mat_clear(b);
    acb_mat_clear(a);
    acb_mat_clear(X);
    return worst;
}

/* density-matrix validity (Hermitian + PSD + |Tr-1|), independent of production. */
static void check_density(const char *name, const acb_mat_t gamma, slong j,
                          const arb_t tol, slong prec)
{
    slong dE = acb_mat_nrows(gamma);
    arb_t e;
    arb_init(e);
    acb_mat_t gd, diff, neg;
    acb_mat_init(gd, dE, dE);
    acb_mat_init(diff, dE, dE);
    acb_mat_conjugate_transpose(gd, gamma);
    acb_mat_sub(diff, gamma, gd, prec);
    aic_mat_opnorm(e, diff, prec);
    AIC_CHECK_MSG(arb_le(e, tol), "%s: gamma_%ld not Hermitian (%.3e)", name,
                  (long) j, dd(e));
    acb_mat_init(neg, dE, dE);
    for (slong p = 0; p < dE; p++)
        for (slong q = 0; q < dE; q++)
            acb_neg(acb_mat_entry(neg, p, q), acb_mat_entry(gamma, p, q));
    aic_mat_herm_max_eig(e, neg, prec);            /* = -lambda_min(gamma) */
    AIC_CHECK_MSG(arb_le(e, tol), "%s: gamma_%ld not PSD (lambda_min=%.3e)", name,
                  (long) j, -dd(e));
    acb_t tr;
    acb_init(tr);
    acb_zero(tr);
    for (slong p = 0; p < dE; p++) acb_add(tr, tr, acb_mat_entry(gamma, p, p), prec);
    acb_sub_si(tr, tr, 1, prec);
    acb_abs(e, tr, prec);
    AIC_CHECK_MSG(arb_le(e, tol), "%s: gamma_%ld trace != 1 (|Tr-1|=%.3e)", name,
                  (long) j, dd(e));
    acb_clear(tr);
    acb_mat_clear(neg);
    acb_mat_clear(diff);
    acb_mat_clear(gd);
    arb_clear(e);
}

/* unitality: Gamma_j(1_M) = 1_{L_j} for each block (=> Gamma(1_M) = 1_A). */
static double unital_defect(const aic_idemp_gamma *G, const aic_idemp_wedderburn *W,
                            slong prec)
{
    slong m = W->dim_M;
    double worst = 0;
    acb_mat_t oneM, Gj, oneL, diff;
    arb_t e;
    arb_init(e);
    acb_mat_init(oneM, m, m);
    acb_mat_one(oneM);
    for (slong j = 0; j < W->num_blocks; j++) {
        slong dL = W->dim_L[j], dE = W->dim_E[j];
        acb_mat_init(Gj, dL, dL);
        gamma_j_tr(Gj, W->W_j[j], G->gamma_j[j], oneM, dL, dE, m, prec);
        acb_mat_init(oneL, dL, dL);
        acb_mat_one(oneL);
        acb_mat_init(diff, dL, dL);
        acb_mat_sub(diff, Gj, oneL, prec);
        aic_mat_opnorm(e, diff, prec);
        if (dd(e) > worst) worst = dd(e);
        acb_mat_clear(diff);
        acb_mat_clear(oneL);
        acb_mat_clear(Gj);
    }
    acb_mat_clear(oneM);
    arb_clear(e);
    return worst;
}

/* Gamma.w = 1_A (.tex:2107): apply Gamma_kraus to w(B_k) (= reshape of d->w col k)
 * and recover B_k. Through w: w(Gamma_kraus(w(B_k))) must equal w(B_k) itself.
 * Returns max_k ||w(Gamma_kraus(w(B_k))) - w(B_k)||_op. */
static double gamma_w_defect(const aic_idemp_gamma *G, const aic_idemp_wedderburn *W,
                             const aic_idemp_decomp *d, slong prec)
{
    slong m = d->dim_M;
    double worst = 0;
    acb_mat_t wbk, gk, diff;
    arb_t e;
    arb_init(e);
    acb_mat_init(wbk, m, m);
    acb_mat_init(gk, m, m);
    acb_mat_init(diff, m, m);
    for (slong k = 0; k < d->dim_A; k++) {
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < m; j++)
                acb_set(acb_mat_entry(wbk, i, j), acb_mat_entry(d->w, i * m + j, k));
        wGk_prop(gk, W, (const acb_mat_t *) G->gamma_j, wbk, prec);
        acb_mat_sub(diff, gk, wbk, prec);
        aic_mat_opnorm(e, diff, prec);
        if (dd(e) > worst) worst = dd(e);
    }
    acb_mat_clear(diff);
    acb_mat_clear(gk);
    acb_mat_clear(wbk);
    arb_clear(e);
    return worst;
}

static void check_gamma(const char *name, aic_ucp_kraus *phi)
{
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, phi, GKPREC);
    aic_idemp_wedderburn W;
    aic_idemp_wedderburn_decompose(&W, &d, GKPREC);
    aic_idemp_gamma G;
    aic_idemp_gamma_kraus(&G, &W, &d, GKPREC);   /* production self-certifies match */

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-10);

    /* (1) the crux: match-d->Gamma (re-measured independently). */
    double mr = match_resid(&W, (const acb_mat_t *) G.gamma_j, &d, GKPREC);
    AIC_CHECK_MSG(mr <= 1e-10, "%s: match-d.Gamma residual %.3e > 1e-10", name, mr);

    /* (2) each gamma_j a valid density matrix. */
    for (slong j = 0; j < G.num_blocks; j++)
        check_density(name, G.gamma_j[j], j, tol, GKPREC);

    /* (3) unitality + Gamma.w = 1_A. */
    double ud = unital_defect(&G, &W, GKPREC);
    AIC_CHECK_MSG(ud <= 1e-10, "%s: Gamma not unital (%.3e)", name, ud);
    double gw = gamma_w_defect(&G, &W, &d, GKPREC);
    AIC_CHECK_MSG(gw <= 1e-10, "%s: Gamma.w != 1_A (%.3e)", name, gw);

    printf("  %-28s blocks=%ld dimE=[", name, (long) G.num_blocks);
    for (slong j = 0; j < G.num_blocks; j++)
        printf("%ld%s", (long) G.dim_E[j], j + 1 < G.num_blocks ? "," : "");
    printf("]  match=%.2e  unital=%.2e  Gamma.w=%.2e", mr, ud, gw);
    if (G.dim_E[0] >= 2)
        printf("  gamma_0[0,0]=%.4f (maximally-mixed=%.4f)",
               arf_get_d(arb_midref(acb_realref(acb_mat_entry(G.gamma_j[0], 0, 0))),
                         ARF_RND_NEAR), 1.0 / (double) G.dim_E[0]);
    printf("\n");

    arb_clear(tol);
    aic_idemp_gamma_clear(&G);
    aic_idemp_wedderburn_clear(&W);
    aic_idemp_clear(&d);
}

/* MUTATION PROOF (Rule 7): perturb gamma_1 on noiseless_subsystem(2,2) and confirm
 * the match-d.Gamma tooth goes RED. d->Gamma is the MAXIMALLY-MIXED conditional
 * expectation (gamma_1 = (1/2)1), so replacing gamma with 1/dE would NOT break it;
 * we instead use a NON-uniform density matrix diag(0.8, 0.2) (still a valid density
 * matrix, but the WRONG conditional expectation) — match-d.Gamma must blow up. */
static void mutation_proof(void)
{
    printf("MUTATION PROOF (noiseless_subsystem(2,2), gamma_1 perturbed):\n");
    aic_ucp_kraus phi;
    make_noiseless_subsystem(&phi, 2, 2);
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, GKPREC);
    aic_idemp_wedderburn W;
    aic_idemp_wedderburn_decompose(&W, &d, GKPREC);
    aic_idemp_gamma G;
    aic_idemp_gamma_kraus(&G, &W, &d, GKPREC);

    /* baseline (true gamma_1 = (1/2)1): GREEN */
    double base = match_resid(&W, (const acb_mat_t *) G.gamma_j, &d, GKPREC);
    AIC_CHECK_MSG(base <= 1e-10, "mutation: baseline match %.3e not green", base);

    /* perturb gamma_1 -> diag(0.8, 0.2) (a valid density matrix, WRONG cond exp) */
    slong dE = G.dim_E[0];
    acb_mat_t saved;
    acb_mat_init(saved, dE, dE);
    acb_mat_set(saved, G.gamma_j[0]);
    acb_mat_zero(G.gamma_j[0]);
    acb_set_d(acb_mat_entry(G.gamma_j[0], 0, 0), 0.8);
    acb_set_d(acb_mat_entry(G.gamma_j[0], 1, 1), 0.2);
    double red = match_resid(&W, (const acb_mat_t *) G.gamma_j, &d, GKPREC);
    printf("  perturbed gamma_1=diag(0.8,0.2): match-d.Gamma = %.3e (RED line)\n", red);
    AIC_CHECK_MSG(red > 1e-3, "mutation: perturbed match %.3e did NOT go red", red);

    /* restore -> GREEN */
    acb_mat_set(G.gamma_j[0], saved);
    double restored = match_resid(&W, (const acb_mat_t *) G.gamma_j, &d, GKPREC);
    printf("  restored gamma_1=(1/2)1:         match-d.Gamma = %.3e (back to green)\n",
           restored);
    AIC_CHECK_MSG(restored <= 1e-10, "mutation: restored match %.3e not green",
                  restored);

    acb_mat_clear(saved);
    aic_idemp_gamma_clear(&G);
    aic_idemp_wedderburn_clear(&W);
    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&phi);
}

static void test_oracles(void)
{
    printf("prop_Gamma Kraus form (eta=0 oracle, match-d.Gamma tooth):\n");
    aic_ucp_kraus phi;

    /* noiseless_subsystem(2,2): E_1 = C^2, the NONTRIVIAL gamma_1 = (1/2)1 — crux. */
    make_noiseless_subsystem(&phi, 2, 2);
    check_gamma("noiseless_subsystem(2,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* block_cond_exp(5,2): M_2 (+) M_3, E_j = C^1, gamma_j = [1]. */
    make_block_cond_exp(&phi, 5, 2);
    check_gamma("block_cond_exp(5,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* identity(3): single trivial block, E = C^1. */
    make_identity(&phi, 3);
    check_gamma("identity(3)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* block_cond_exp(4,2): M_2 (+) M_2, E_j = C^1. */
    make_block_cond_exp(&phi, 4, 2);
    check_gamma("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_oracles();
    mutation_proof();
    aic_test_report("test_gamma_kraus");
    printf("OK test_gamma_kraus\n");
    return 0;
}
