/* test_ucp_d24.c — cross-checks for UCP composition + Choi-difference (bead
 * aic-d24, the eta-idempotence-defect substrate). The diamond-norm SDP itself is
 * Julia+MOSEK (Route B); here we exercise the C substrate it feeds on:
 *
 *   PHASE 1 (this file initially):
 *     1. compose correctness: aic_ucp_apply(compose(Phi,Psi), X) == Phi(Psi(X))
 *        for random X (the DEFINITION of composition, an independent path); and
 *        for Phi^2: apply(compose(Phi,Phi),X) == Phi(Phi(X)).
 *     2. compose preserves unitality (Phi,Psi unital => Phi o Psi unital).
 *     3. choi_diff Hermitian; eta=0 oracle: idempotent Phi => Choi(Phi^2-Phi)==0
 *        to machine zero (the cleanest ground truth, CLAUDE.md cross-check #3).
 *     4. double-vs-arb@53 on the composed Kraus ops and on Choi_diff entries.
 *     5. MUTATION PROOF: wrong product order (K_a L_b instead of L_b K_a) breaks
 *        the compose cross-check on a NON-COMMUTING example.
 *
 *   PHASE 3 (added once fixtures_d24.inc.h exists): for each Julia-generated
 *     fixture, build the channel, compute Phi^2 via compose, Choi(Phi^2-Phi) via
 *     choi_diff, and assert it matches the fixture's reference Choi (C-vs-Julia
 *     cross-check); assert eta=0 fixtures carry reference eta ~ 0.
 *
 * Convention reminders (include/aic_ucp.h): Kraus K_a are dim_K x dim_H,
 * Phi(X)=sum K_a^dag X K_a (Heisenberg); compose Kraus = {psi.K[b] @ phi.K[a]}.
 */
#include <assert.h>
#include <complex.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic_latd.h"
#include "aic_mat.h"
#include "aic_ucp.h"
#include "aic_test.h"

#if __has_include("fixtures_d24.inc.h")
#include "fixtures_d24.inc.h"
#define AIC_HAVE_FIXTURES 1
#endif

#define PREC 53
static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* deterministic complex pseudo-random fill (LCG) for test matrices. */
static void fill_random(acb_mat_t A, unsigned long *s)
{
    for (slong i = 0; i < acb_mat_nrows(A); i++)
        for (slong j = 0; j < acb_mat_ncols(A); j++) {
            *s = *s * 6364136223846793005UL + 1442695040888963407UL;
            double re = (double) ((*s >> 33) % 2000) / 1000.0 - 1.0;
            *s = *s * 6364136223846793005UL + 1442695040888963407UL;
            double im = (double) ((*s >> 33) % 2000) / 1000.0 - 1.0;
            acb_set_d_d(acb_mat_entry(A, i, j), re, im);
        }
}

/* ---- channel constructors ---- */

/* Genuinely complex, asymmetric, unital channel on C^2 (from test_ucp.c):
 *   K_0 = [[0.6,0],[0.8i,0]], K_1 = [[0,0.8i],[0,0.6]]. */
static void make_complex_channel(aic_ucp_kraus *phi)
{
    aic_ucp_kraus_init(phi, 2, 2, 2);
    acb_set_d_d(acb_mat_entry(phi->K[0], 0, 0), 0.6, 0.0);
    acb_set_d_d(acb_mat_entry(phi->K[0], 1, 0), 0.0, 0.8);
    acb_set_d_d(acb_mat_entry(phi->K[1], 0, 1), 0.0, 0.8);
    acb_set_d_d(acb_mat_entry(phi->K[1], 1, 1), 0.6, 0.0);
}

/* Dephasing M_2 on C^2: Phi(X) = diag(X) (kills off-diagonals). UNITAL,
 * idempotent. Kraus K_0 = |0><0|, K_1 = |1><1|: sum K_a^dag X K_a =
 * X[0,0]|0><0| + X[1,1]|1><1| = diag(X). Idempotent: diag(diag(X)) = diag(X). */
static void make_dephasing(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, d);
    for (slong a = 0; a < d; a++)
        acb_set_si(acb_mat_entry(phi->K[a], a, a), 1);
}

/* Identity channel on C^d (UCP, idempotent). */
static void make_identity(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, 1);
    acb_mat_one(phi->K[0]);
}

/* Completely depolarizing Dep_d: Phi(X) = Tr(X) 1_d / d, K_{ij}=|e_i><e_j|/sqrt d.
 * UNITAL, idempotent (Dep^2 = Dep: Tr(Tr(X)1/d)1/d = Tr(X)1/d). */
static void make_depolarizing(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, d * d);
    double inv = 1.0 / sqrt((double) d);
    slong a = 0;
    for (slong i = 0; i < d; i++)
        for (slong j = 0; j < d; j++) {
            acb_set_d(acb_mat_entry(phi->K[a], i, j), inv);
            a++;
        }
}

/* Block conditional expectation on C^4 = C^2 (+) C^2: project onto the
 * 2x2 block-diagonal subalgebra. Phi(X) = P_0 X P_0 + P_1 X P_1 with
 * P_0 = diag(1,1,0,0), P_1 = diag(0,0,1,1). Kraus {P_0, P_1}. UNITAL
 * (P_0+P_1=1), idempotent (P_i X P_j cross-terms vanish under reapplication;
 * Phi^2(X) = P_0(P_0 X P_0 + P_1 X P_1)P_0 + ... = P_0 X P_0 + P_1 X P_1). */
static void make_block_cond_exp(aic_ucp_kraus *phi)
{
    aic_ucp_kraus_init(phi, 4, 4, 2);
    acb_set_si(acb_mat_entry(phi->K[0], 0, 0), 1);
    acb_set_si(acb_mat_entry(phi->K[0], 1, 1), 1);
    acb_set_si(acb_mat_entry(phi->K[1], 2, 2), 1);
    acb_set_si(acb_mat_entry(phi->K[1], 3, 3), 1);
}

/* Trace-replace channel on C^d: Phi(X) = Tr(X) |0><0| ... not unital. We instead
 * use a unital "amplitude-replace": Phi(X) = (1/d) sum_i <e_i|X|e_i> 1_d, which
 * is exactly Dep_d, already covered. Use the asymmetric non-self-composition test
 * below for the K!=H chaining instead. */

/* ---- test 1: compose correctness (independent path) ---- */
static void test_compose_correctness(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-13);
    unsigned long s = 0x1234567ULL;

    /* (a) Phi=complex_channel, Psi=dephasing (both C^2 self-maps), non-commuting.
     * apply(compose(Phi,Psi), X) must equal Phi(Psi(X)) for random X. */
    aic_ucp_kraus phi, psi, comp;
    make_complex_channel(&phi);
    make_dephasing(&psi, 2);
    aic_ucp_compose(&comp, &phi, &psi, PREC);
    AIC_CHECK_MSG(comp.r == phi.r * psi.r, "compose r wrong");
    AIC_CHECK_MSG(comp.dim_K == psi.dim_K && comp.dim_H == phi.dim_H,
                  "compose dims wrong");

    acb_mat_t X, psiX, phipsiX, compX;
    acb_mat_init(X, 2, 2);
    acb_mat_init(psiX, 2, 2);
    acb_mat_init(phipsiX, 2, 2);
    acb_mat_init(compX, 2, 2);
    for (int t = 0; t < 4; t++) {
        fill_random(X, &s);
        aic_ucp_apply(psiX, &psi, X, PREC);          /* Psi(X)        */
        aic_ucp_apply(phipsiX, &phi, psiX, PREC);     /* Phi(Psi(X))   */
        aic_ucp_apply(compX, &comp, X, PREC);         /* (PhioPsi)(X)  */
        AIC_CHECK_ACB_MAT_CLOSE(compX, phipsiX, tol);
    }
    acb_mat_clear(compX);
    acb_mat_clear(phipsiX);
    acb_mat_clear(psiX);
    acb_mat_clear(X);
    aic_ucp_kraus_clear(&comp);
    aic_ucp_kraus_clear(&psi);
    aic_ucp_kraus_clear(&phi);

    /* (b) Phi^2 = compose(Phi,Phi): apply(Phi2,X) == Phi(Phi(X)). Use the
     * depolarizing channel (large r => r^2 composed ops, exercises the loop). */
    aic_ucp_kraus dep, dep2;
    make_depolarizing(&dep, 3);
    aic_ucp_compose(&dep2, &dep, &dep, PREC);
    AIC_CHECK_MSG(dep2.r == dep.r * dep.r, "Phi^2 r wrong");
    acb_mat_t Y, phiY, phiphiY, phi2Y;
    acb_mat_init(Y, 3, 3);
    acb_mat_init(phiY, 3, 3);
    acb_mat_init(phiphiY, 3, 3);
    acb_mat_init(phi2Y, 3, 3);
    for (int t = 0; t < 3; t++) {
        fill_random(Y, &s);
        aic_ucp_apply(phiY, &dep, Y, PREC);
        aic_ucp_apply(phiphiY, &dep, phiY, PREC);
        aic_ucp_apply(phi2Y, &dep2, Y, PREC);
        AIC_CHECK_ACB_MAT_CLOSE(phi2Y, phiphiY, tol);
    }
    acb_mat_clear(phi2Y);
    acb_mat_clear(phiphiY);
    acb_mat_clear(phiY);
    acb_mat_clear(Y);
    aic_ucp_kraus_clear(&dep2);
    aic_ucp_kraus_clear(&dep);

    /* (c) NON-self-map chaining K!=H: Psi: B(C^3)->B(C^2), Phi: B(C^2)->B(C^3),
     * compose Phi o Psi : B(C^3)->B(C^3). Requires phi->dim_K(=2)==psi->dim_H(=2).
     * Psi has Kraus L_b: dim_K=3 x dim_H=2; Phi has K_a: dim_K=2 x dim_H=3. */
    {
        aic_ucp_kraus p, q, c;
        aic_ucp_kraus_init(&p, 2, 3, 2); /* Phi: B(C^2)->B(C^3), K_a 2x3 */
        aic_ucp_kraus_init(&q, 3, 2, 2); /* Psi: B(C^3)->B(C^2), L_b 3x2 */
        unsigned long s2 = 0xBEEFULL;
        fill_random(p.K[0], &s2); fill_random(p.K[1], &s2);
        fill_random(q.K[0], &s2); fill_random(q.K[1], &s2);
        aic_ucp_compose(&c, &p, &q, PREC); /* B(C^3)->B(C^3): dim_K=3,dim_H=3 */
        AIC_CHECK_MSG(c.dim_K == 3 && c.dim_H == 3 && c.r == 4,
                      "K!=H compose shape wrong");
        acb_mat_t Z, qZ, pqZ, cZ;
        acb_mat_init(Z, 3, 3);
        acb_mat_init(qZ, 2, 2);
        acb_mat_init(pqZ, 3, 3);
        acb_mat_init(cZ, 3, 3);
        fill_random(Z, &s);
        aic_ucp_apply(qZ, &q, Z, PREC);     /* Psi(Z) in B(C^2)     */
        aic_ucp_apply(pqZ, &p, qZ, PREC);   /* Phi(Psi(Z)) in B(C^3) */
        aic_ucp_apply(cZ, &c, Z, PREC);     /* (PhioPsi)(Z)          */
        AIC_CHECK_ACB_MAT_CLOSE(cZ, pqZ, tol);
        acb_mat_clear(cZ); acb_mat_clear(pqZ); acb_mat_clear(qZ); acb_mat_clear(Z);
        aic_ucp_kraus_clear(&c); aic_ucp_kraus_clear(&q); aic_ucp_kraus_clear(&p);
    }

    arb_clear(tol);
}

/* ---- test 2: compose preserves unitality ---- */
static void test_compose_unital(void)
{
    arb_t tol, def;
    arb_init(tol);
    arb_init(def);
    set_tol(tol, 1e-13);

    /* both factors unital (complex_channel, dephasing) => Phi o Psi unital. */
    aic_ucp_kraus phi, psi, comp;
    make_complex_channel(&phi);
    make_dephasing(&psi, 2);
    aic_ucp_compose(&comp, &phi, &psi, PREC);
    aic_ucp_unital_defect_kraus(def, &comp, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "compose of unital maps not unital");
    aic_ucp_kraus_clear(&comp);
    aic_ucp_kraus_clear(&psi);
    aic_ucp_kraus_clear(&phi);

    /* Dep_3^2 unital. */
    aic_ucp_kraus dep, dep2;
    make_depolarizing(&dep, 3);
    aic_ucp_compose(&dep2, &dep, &dep, PREC);
    aic_ucp_unital_defect_kraus(def, &dep2, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "Dep^2 not unital");
    aic_ucp_kraus_clear(&dep2);
    aic_ucp_kraus_clear(&dep);

    arb_clear(def);
    arb_clear(tol);
}

/* ---- test 3: choi_diff Hermitian + eta=0 oracle (Choi(Phi^2-Phi)=0) ---- */
static void check_idemp_choi_zero(aic_ucp_kraus *phi, const char *name)
{
    arb_t tol, nrm;
    arb_init(tol);
    arb_init(nrm);
    set_tol(tol, 1e-12);
    slong dk = phi->dim_K, dh = phi->dim_H, n = dk * dh;

    aic_ucp_kraus phi2;
    aic_ucp_compose(&phi2, phi, phi, PREC);

    acb_mat_t L;
    acb_mat_init(L, n, n);
    aic_ucp_choi_diff(L, &phi2, phi, PREC); /* Choi(Phi^2 - Phi) */

    /* idempotent => Lambda = 0 => Choi(Lambda) = 0 (the eta=0 oracle at the
     * Choi level: the cleanest ground truth, CLAUDE.md cross-check #3). */
    aic_mat_opnorm(nrm, L, PREC);
    AIC_CHECK_MSG(arb_le(nrm, tol),
                  "%s: Choi(Phi^2-Phi) not zero (eta=0 oracle)", name);

    /* Hermiticity: ||L - L^dag||_op ~ 0. */
    acb_mat_t Ld, asym;
    arb_t hdef;
    acb_mat_init(Ld, n, n);
    acb_mat_init(asym, n, n);
    arb_init(hdef);
    acb_mat_conjugate_transpose(Ld, L);
    acb_mat_sub(asym, L, Ld, PREC);
    aic_mat_opnorm(hdef, asym, PREC);
    AIC_CHECK_MSG(arb_le(hdef, tol), "%s: Choi(Lambda) not Hermitian", name);

    arb_clear(hdef);
    acb_mat_clear(asym);
    acb_mat_clear(Ld);
    acb_mat_clear(L);
    aic_ucp_kraus_clear(&phi2);
    arb_clear(nrm);
    arb_clear(tol);
}

static void test_choi_diff_eta0(void)
{
    /* the four idempotent test channels named in the bead. */
    aic_ucp_kraus phi;

    make_identity(&phi, 3);
    check_idemp_choi_zero(&phi, "identity");
    aic_ucp_kraus_clear(&phi);

    make_depolarizing(&phi, 2);
    check_idemp_choi_zero(&phi, "Dep_2");
    aic_ucp_kraus_clear(&phi);

    make_depolarizing(&phi, 3);
    check_idemp_choi_zero(&phi, "Dep_3");
    aic_ucp_kraus_clear(&phi);

    make_dephasing(&phi, 2);
    check_idemp_choi_zero(&phi, "dephasing_M2");
    aic_ucp_kraus_clear(&phi);

    make_block_cond_exp(&phi);
    check_idemp_choi_zero(&phi, "block_cond_exp");
    aic_ucp_kraus_clear(&phi);

    /* a NON-idempotent control: Choi(Phi^2-Phi) must be NONzero. Phi_t =
     * (1-t)id + t Dep on C^2 with t=0.3 is not idempotent (Phi^2-Phi =
     * t(1-t)(Dep-id) != 0). Build it as a 1+4=5-op Kraus mixture? Simpler:
     * convex-combine the two channels' actions by stacking weighted Kraus ops.
     * Phi_t(X) = (1-t) X + t Tr(X)1/d. Kraus: sqrt(1-t) 1_d (id part) and
     * sqrt(t) |e_i><e_j|/sqrt d (depol part). */
    {
        const slong d = 2;
        const double tt = 0.3;
        aic_ucp_kraus pt;
        aic_ucp_kraus_init(&pt, d, d, 1 + d * d);
        acb_set_d(acb_mat_entry(pt.K[0], 0, 0), sqrt(1.0 - tt));
        acb_set_d(acb_mat_entry(pt.K[0], 1, 1), sqrt(1.0 - tt));
        double sc = sqrt(tt) / sqrt((double) d);
        slong a = 1;
        for (slong i = 0; i < d; i++)
            for (slong j = 0; j < d; j++)
                acb_set_d(acb_mat_entry(pt.K[a++], i, j), sc);

        /* unital sanity */
        arb_t def, tol;
        arb_init(def); arb_init(tol);
        set_tol(tol, 1e-12);
        aic_ucp_unital_defect_kraus(def, &pt, PREC);
        AIC_CHECK_MSG(arb_le(def, tol), "Phi_t not unital");

        aic_ucp_kraus pt2;
        aic_ucp_compose(&pt2, &pt, &pt, PREC);
        slong n = d * d;
        acb_mat_t L;
        acb_mat_init(L, n, n);
        aic_ucp_choi_diff(L, &pt2, &pt, PREC);
        arb_t nrm;
        arb_init(nrm);
        aic_mat_opnorm(nrm, L, PREC);
        arb_t small;
        arb_init(small);
        arb_set_d(small, 0.01);
        AIC_CHECK_MSG(arb_gt(nrm, small),
                      "Phi_t (t=0.3) Choi(Lambda) should be bounded away from 0");
        arb_clear(small);
        arb_clear(nrm);
        acb_mat_clear(L);
        aic_ucp_kraus_clear(&pt2);
        arb_clear(tol);
        arb_clear(def);
        aic_ucp_kraus_clear(&pt);
    }
}

/* ---- test 4: double-vs-arb@53 on composed Kraus + Choi_diff ---- */
static void test_double_vs_arb(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-11);
    unsigned long s = 0xACE1ULL;

    /* Phi = complex_channel, Psi = a random (non-unital, that's fine here) C^2
     * self-map. Composed Kraus L_b K_a computed on the arb path must match a
     * LAPACK-double recomputation of the same products. */
    aic_ucp_kraus phi, psi, comp;
    make_complex_channel(&phi);
    aic_ucp_kraus_init(&psi, 2, 2, 2);
    fill_random(psi.K[0], &s);
    fill_random(psi.K[1], &s);
    aic_ucp_compose(&comp, &phi, &psi, PREC);

    /* recompute each L_b K_a on the double path and compare to comp.K[c]. */
    double _Complex La[4], Ka[4], Md[4];
    slong c = 0;
    for (slong a = 0; a < phi.r; a++)
        for (slong b = 0; b < psi.r; b++) {
            aic_latd_from_acb_mat(La, psi.K[b]); /* 2x2 row-major */
            aic_latd_from_acb_mat(Ka, phi.K[a]);
            /* Md = La @ Ka (2x2) */
            for (int ii = 0; ii < 2; ii++)
                for (int jj = 0; jj < 2; jj++) {
                    double _Complex acc = 0;
                    for (int kk = 0; kk < 2; kk++)
                        acc += La[ii * 2 + kk] * Ka[kk * 2 + jj];
                    Md[ii * 2 + jj] = acc;
                }
            acb_mat_t Mref;
            acb_mat_init(Mref, 2, 2);
            aic_latd_to_acb_mat(Mref, Md, 2, 2);
            AIC_CHECK_ACB_MAT_CLOSE(comp.K[c], Mref, tol);
            acb_mat_clear(Mref);
            c++;
        }

    /* Choi_diff double-vs-arb: Choi(Phi^2-Phi) for an idempotent channel is 0 on
     * both paths; for a non-idempotent one, compare the arb Choi_diff against a
     * direct double recomputation of C_{Phi^2}-C_Phi entrywise. Use dephasing
     * (idempotent, exact 0) AND a non-idempotent depolarizing-mix is overkill —
     * verify entrywise on dephasing's NON-trivial sub-step: instead use the arb
     * Choi_diff of two DISTINCT channels (phi vs psi) so entries are nonzero. */
    {
        slong n = 4;
        acb_mat_t L;
        acb_mat_init(L, n, n);
        aic_ucp_choi_diff(L, &phi, &psi, PREC); /* Choi(phi) - Choi(psi) */
        /* double recomputation: Choi_x[i*2+a, j*2+b] = sum conj(Kx[i,a])Kx[j,b] */
        double _Complex Pa[2][4]; /* phi Kraus row-major 2x2 */
        double _Complex Qa[2][4]; /* psi Kraus */
        for (slong x = 0; x < 2; x++) {
            aic_latd_from_acb_mat(Pa[x], phi.K[x]);
            aic_latd_from_acb_mat(Qa[x], psi.K[x]);
        }
        for (slong i = 0; i < 2; i++)
            for (slong a = 0; a < 2; a++)
                for (slong j = 0; j < 2; j++)
                    for (slong bb = 0; bb < 2; bb++) {
                        double _Complex cp = 0, cq = 0;
                        for (slong x = 0; x < 2; x++) {
                            cp += conj(Pa[x][i * 2 + a]) * Pa[x][j * 2 + bb];
                            cq += conj(Qa[x][i * 2 + a]) * Qa[x][j * 2 + bb];
                        }
                        double _Complex want = cp - cq;
                        acb_t got;
                        acb_init(got);
                        acb_set(got, acb_mat_entry(L, i * 2 + a, j * 2 + bb));
                        double gre = arf_get_d(arb_midref(acb_realref(got)),
                                               ARF_RND_NEAR);
                        double gim = arf_get_d(arb_midref(acb_imagref(got)),
                                               ARF_RND_NEAR);
                        AIC_CHECK_MSG(fabs(gre - creal(want)) < 1e-11 &&
                                      fabs(gim - cimag(want)) < 1e-11,
                                      "Choi_diff double-vs-arb off at (%ld,%ld)",
                                      (long) (i * 2 + a), (long) (j * 2 + bb));
                        acb_clear(got);
                    }
        acb_mat_clear(L);
    }

    aic_ucp_kraus_clear(&comp);
    aic_ucp_kraus_clear(&psi);
    aic_ucp_kraus_clear(&phi);
    arb_clear(tol);
}

/* ---- test 5: mutation proof — wrong product order breaks compose ---- */
static void test_mutation_product_order(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-13);
    unsigned long s = 0xD24ULL;

    /* On a NON-commuting pair, (PhioPsi)(X)=Phi(Psi(X)) requires Kraus L_b K_a;
     * the wrong order K_a L_b gives a different channel. We construct the WRONG
     * composition by hand (swap the product) and confirm it FAILS the cross-check
     * — proving the test has teeth (Rule 7 mutation proof, done in-test). */
    aic_ucp_kraus phi, psi;
    make_complex_channel(&phi);            /* non-commuting with dephasing */
    make_dephasing(&psi, 2);

    /* WRONG-order composition: K_a L_b (= phi.K[a] @ psi.K[b]) instead of
     * L_b K_a. Both are 2x2 here so shapes match; only the value differs. */
    aic_ucp_kraus wrong;
    aic_ucp_kraus_init(&wrong, 2, 2, phi.r * psi.r);
    slong c = 0;
    for (slong a = 0; a < phi.r; a++)
        for (slong b = 0; b < psi.r; b++)
            acb_mat_mul(wrong.K[c++], phi.K[a], psi.K[b], PREC); /* K_a L_b */

    acb_mat_t X, psiX, phipsiX, wrongX;
    acb_mat_init(X, 2, 2);
    acb_mat_init(psiX, 2, 2);
    acb_mat_init(phipsiX, 2, 2);
    acb_mat_init(wrongX, 2, 2);
    int any_mismatch = 0;
    for (int t = 0; t < 4; t++) {
        fill_random(X, &s);
        aic_ucp_apply(psiX, &psi, X, PREC);
        aic_ucp_apply(phipsiX, &phi, psiX, PREC);  /* correct Phi(Psi(X)) */
        aic_ucp_apply(wrongX, &wrong, X, PREC);    /* wrong-order action  */
        slong fi, fj;
        if (!aic_acb_mat_close(wrongX, phipsiX, tol, &fi, &fj)) any_mismatch = 1;
    }
    AIC_CHECK_MSG(any_mismatch,
                  "MUTATION: wrong product order (K_a L_b) was NOT caught — "
                  "compose cross-check lacks teeth on a non-commuting pair");

    acb_mat_clear(wrongX);
    acb_mat_clear(phipsiX);
    acb_mat_clear(psiX);
    acb_mat_clear(X);
    aic_ucp_kraus_clear(&wrong);
    aic_ucp_kraus_clear(&psi);
    aic_ucp_kraus_clear(&phi);
    arb_clear(tol);
}

#ifdef AIC_HAVE_FIXTURES
/* ---- Phase 3: C-vs-Julia Choi cross-check against the generated fixtures ---- */
static void build_kraus_from_fixture(aic_ucp_kraus *phi, const aic_d24_fixture *f)
{
    aic_ucp_kraus_init(phi, f->n, f->n, f->r);
    /* fixture Kraus layout: flat [op][dim_K*dim_H] row-major, re and im split. */
    for (slong a = 0; a < f->r; a++)
        for (slong i = 0; i < f->n; i++)
            for (slong j = 0; j < f->n; j++) {
                slong off = a * f->n * f->n + i * f->n + j;
                acb_set_d_d(acb_mat_entry(phi->K[a], i, j),
                            f->kraus_re[off], f->kraus_im[off]);
            }
}

static void test_fixtures_choi(void)
{
    arb_t tol, etatol;
    arb_init(tol);
    arb_init(etatol);
    set_tol(tol, 1e-11);
    set_tol(etatol, 1e-9);

    for (int fi = 0; fi < AIC_D24_NFIX; fi++) {
        const aic_d24_fixture *f = &aic_d24_fixtures[fi];
        slong n = f->n, N = n * n;
        aic_ucp_kraus phi, phi2;
        build_kraus_from_fixture(&phi, f);
        aic_ucp_compose(&phi2, &phi, &phi, PREC);

        acb_mat_t L, Lref;
        acb_mat_init(L, N, N);
        acb_mat_init(Lref, N, N);
        aic_ucp_choi_diff(L, &phi2, &phi, PREC);
        /* fixture reference Choi(Phi^2-Phi), flat N x N row-major re/im. */
        for (slong p = 0; p < N; p++)
            for (slong q = 0; q < N; q++)
                acb_set_d_d(acb_mat_entry(Lref, p, q),
                            f->choi_re[p * N + q], f->choi_im[p * N + q]);
        slong di, dj;
        AIC_CHECK_MSG(aic_acb_mat_close(L, Lref, tol, &di, &dj),
                      "fixture %s: C Choi(Lambda) != Julia ref at (%ld,%ld)",
                      f->name, (long) di, (long) dj);

        /* eta=0 fixtures must carry reference eta ~ 0 (single source of the
         * tolerance: the etatol arb_t, set to 1e-9 above). */
        if (f->eta_is_zero) {
            arb_t eref;
            arb_init(eref);
            arb_set_d(eref, f->eta_ref);
            AIC_CHECK_MSG(arb_lt(eref, etatol),
                          "fixture %s: eta=0 oracle but ref eta=%.3e",
                          f->name, f->eta_ref);
            arb_clear(eref);
        }

        acb_mat_clear(Lref);
        acb_mat_clear(L);
        aic_ucp_kraus_clear(&phi2);
        aic_ucp_kraus_clear(&phi);
    }
    arb_clear(etatol);
    arb_clear(tol);
}
#endif

int main(void)
{
    test_compose_correctness();
    test_compose_unital();
    test_choi_diff_eta0();
    test_double_vs_arb();
    test_mutation_product_order();
#ifdef AIC_HAVE_FIXTURES
    test_fixtures_choi();
    printf("fixtures: %d cross-checked\n", AIC_D24_NFIX);
#else
    printf("fixtures: none (fixtures_d24.inc.h not generated yet)\n");
#endif
    aic_test_report("test_ucp_d24");
    printf("OK test_ucp_d24\n");
    return 0;
}
