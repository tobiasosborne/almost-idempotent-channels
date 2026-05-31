/* test_projection.c — cross-checks for the nontrivial-projection finder (bead
 * aic-mqf, module "projection"): lem_nontriv_projection (approximate_algebras.tex
 * :931-932). Built on the eps-C* algebra data model (aic_ecstar) produced from a
 * UCP map by the assoc module (aic_assoc_ecstar_from_phi). Reuses the
 * exactly-idempotent + eta>0 channels in tests/test_idemp.h.
 *
 * THE TWO BLINDNESSES this suite must NOT have (the corner lesson, CLAUDE.md
 * #1/#2 bugs): an identity-channel eta=0 fixture is blind to whether the finder
 * uses A's STAR (Phi_tilde=id => star==plain) and to whether Pi_A is load-bearing
 * (P_amb already in A). So:
 *   - T2 runs on a GENUINE eta>0 OBLIQUE channel (make_mixconj) and asserts the
 *     Pi_A step actually MOVES P_amb (||P-P_amb|| > 0), proving Pi_A is exercised
 *     and the star defect is the real (not plain) defect.
 *   - non-vacuity guards everywhere: the maps are genuinely nonzero, eta>0 is
 *     genuinely > 0, the gap is the one predicted.
 *
 * Cross-check ladder:
 *   T1 eta=0 ORACLE (headline, Rule 6): A = M_d from the identity channel. The
 *      finder returns a spectral projector (rank determined by the gap side);
 *      ||P*P-P||=0 machine-zero, ||P||=1, ||I-P||=1, P^dag=P, P in A (residual 0).
 *      Plus a BLOCK fixture (block_cond_exp): A = M_m (+) M_{d-m} has nontrivial
 *      projections of rank 1..m and 1..(d-m) (all valid in M_m (+) M_{d-m}); the
 *      finder returns one of them (for d=4,m=2 it is rank-1, eigenvalues 0,0,0,1).
 *      Verify delta=0, nontrivial, P in A (gauge-free) — rank is not asserted.
 *      MUTATION-A: skip Pi_A (P=P_amb) -> green at eta=0 (P_amb in A), RED at
 *      eta>0 (T2). MUTATION-B: threshold t=lam_max+1 -> all eigvals one side ->
 *      P trivial -> nontriviality assert RED.
 *   T2 eta>0 (make_mixconj, genuinely oblique): delta = ||P*P-P|| <= C*eta
 *      (report C); nontriviality holds; P in A (proj_residual O(eta)). NON-VACUITY:
 *      eta_proxy>0 AND ||P-P_amb||=O(eta)>0 (Pi_A load-bearing). MUTATION-A: skip
 *      Pi_A -> proj_residual O(eta) -> if asserted, RED.
 *   T3 universality canary (aic-dbo.3): A=M_d (identity) for d=2,3,4: delta/eps
 *      (here machine-zero) and the constant do NOT grow with d (flat). The eta>0
 *      delta/eta constant across d=4,5 is also reported flat (tex:484).
 *   T4 fail-loud: dim A == 1 (A = C*I, trace_replace) -> the routine ABORTS (via
 *      fork/expected-SIGABRT). Confirms the no-nontrivial-projection path fires.
 *   T5 double vs arb@prec=53 (cross-check ladder rung 2): the finder at prec=53
 *      and prec=256 agree on delta/||P|| on a fixed eta=0 instance.
 */
/* mkstemp/dup2/fileno are POSIX; under -std=c11 they need the feature-test macro
 * (the Makefile adds it only on bench compiles, so define it here for T4's
 * stderr-capture fork helper). Must precede all system includes. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "../src/aic_projection_internal.h"   /* aic_projection_gap (no-gap test) */
#include "aic/aic_assoc.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_funcalc.h"
#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic/aic_projection.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"
#include "test_idemp.h"

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Suppress -Wunused-function for the test_idemp.h builders this TU does not use. */
__attribute__((unused))
static void suppress_unused_idemp_builders(void)
{
    (void) make_noiseless_subsystem;
    (void) make_dephasing;
    (void) make_compress_idemp;
}

/* ||A||_op on ball midpoints (mirrors test_corner: a near-zero difference matrix
 * with wide accumulated balls false-fails the certified Gram check, aic-2yo;
 * collapse to midpoints when we only want the number). */
static double mid_opnorm_d(const acb_mat_t Aop, slong prec)
{
    slong r = acb_mat_nrows(Aop), c = acb_mat_ncols(Aop);
    acb_mat_t M;
    arb_t e;
    acb_mat_init(M, r, c);
    arb_init(e);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(Aop, i, j));
    aic_mat_opnorm(e, M, prec);
    double v = dd(e);
    arb_clear(e);
    acb_mat_clear(M);
    return v;
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op (superop idempotence defect, ~ ||Phi^2-Phi||
 * <= eta). Mirrors test_corner eta_proxy. */
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
    double o = mid_opnorm_d(D, prec);
    acb_mat_clear(D); acb_mat_clear(S2); acb_mat_clear(S);
    return o;
}

/* ||P^dag - P||_op (Hermiticity defect of the returned P). */
static double herm_defect(const acb_mat_t P, slong prec)
{
    slong n = acb_mat_nrows(P);
    acb_mat_t Pd, D;
    acb_mat_init(Pd, n, n);
    acb_mat_init(D, n, n);
    acb_mat_conjugate_transpose(Pd, P);
    acb_mat_sub(D, Pd, P, prec);
    double v = mid_opnorm_d(D, prec);
    acb_mat_clear(D); acb_mat_clear(Pd);
    return v;
}

/* Frobenius in-A residual ||X - Pi_A(X)||_F (test-local; avoids the
 * aic_ecstar_proj_residual Gram-false-fail on an in-A input, beads aic-qgs/2yo).
 * Pi_A(X) = sum_k <B_k,X>_F B_k. 0 iff X in A. */
static double in_A_resid(const aic_ecstar *A, const acb_mat_t X, slong prec)
{
    slong n = A->n, d = A->dim_A;
    acb_mat_t R;
    acb_mat_init(R, n, n);
    acb_mat_set(R, X);
    acb_t c, t;
    acb_init(c); acb_init(t);
    for (slong k = 0; k < d; k++) {
        acb_zero(c);
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++) {
                acb_conj(t, acb_mat_entry(A->B[k], i, j));
                acb_mul(t, t, acb_mat_entry(X, i, j), prec);
                acb_add(c, c, t, prec);
            }
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++) {
                acb_mul(t, c, acb_mat_entry(A->B[k], i, j), prec);
                acb_sub(acb_mat_entry(R, i, j), acb_mat_entry(R, i, j), t, prec);
            }
    }
    arb_t e;
    arb_init(e);
    acb_mat_frobenius_norm(e, R, prec);
    double v = dd(e);
    arb_clear(e);
    acb_clear(t); acb_clear(c);
    acb_mat_clear(R);
    return v;
}

/* Common assertions for a returned projection: Hermitian, nontrivial (both
 * ||P|| and ||I-P|| bounded away from 0), and in A. `name` for messages. */
static void check_nontrivial(const char *name, const acb_mat_t P,
                             const arb_t Pnorm, const arb_t ImPnorm,
                             const aic_ecstar *A, double resid_tol, slong prec)
{
    double pn = dd(Pnorm), ipn = dd(ImPnorm);
    AIC_CHECK_MSG(pn > 0.3, "%s: ||P||=%.4f <= 0.3 (P near 0, trivial)", name, pn);
    AIC_CHECK_MSG(ipn > 0.3, "%s: ||I-P||=%.4f <= 0.3 (P near I, trivial)",
                  name, ipn);
    double hd = herm_defect(P, prec);
    AIC_CHECK_MSG(hd < 1e-9, "%s: P not Hermitian (%.3e)", name, hd);
    double rr = in_A_resid(A, P, prec);
    AIC_CHECK_MSG(rr <= resid_tol, "%s: P not in A (resid=%.3e > %.3e)",
                  name, rr, resid_tol);
}

/* =========================== T1: eta=0 ORACLE ============================ */
static void t1_identity(slong d, slong prec)
{
    aic_ucp_kraus phi;
    make_identity(&phi, d);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    AIC_CHECK_MSG(ae.A.dim_A == d * d, "T1(d=%ld): identity dim_A=%ld != d^2",
                  (long) d, (long) ae.A.dim_A);

    acb_mat_t P;
    arb_t delta, Pn, IPn;
    acb_mat_init(P, d, d);
    arb_init(delta); arb_init(Pn); arb_init(IPn);
    aic_projection_witness w;
    aic_projection_nontrivial(P, delta, Pn, IPn, &ae.A, &w, prec);

    double del = dd(delta);
    printf("T1 identity d=%ld: k=%ld spread=[%.3f,%.3f] t=%.3f gap=%.3f m=%ld "
           "delta=%.3e ||P||=%.6f ||I-P||=%.6f\n", (long) d, (long) w.k_chosen,
           w.lam_min, w.lam_max, w.t, w.gap, (long) w.m, del, dd(Pn), dd(IPn));
    /* eta=0: P_amb in A, so delta = ||P*P-P|| is machine-zero. */
    AIC_CHECK_MSG(del < 1e-9, "T1(d=%ld): delta=%.3e not machine-zero at eta=0",
                  (long) d, del);
    /* nontrivial + Hermitian + in A (residual machine-zero at eta=0). */
    check_nontrivial("T1 identity", P, Pn, IPn, &ae.A, 1e-9, prec);
    /* eta=0 identity: P is an EXACT idempotent (plain == star), ||P||=||I-P||=1. */
    AIC_CHECK_MSG(fabs(dd(Pn) - 1.0) < 1e-6, "T1(d=%ld): ||P||=%.6f != 1",
                  (long) d, dd(Pn));
    AIC_CHECK_MSG(fabs(dd(IPn) - 1.0) < 1e-6, "T1(d=%ld): ||I-P||=%.6f != 1",
                  (long) d, dd(IPn));
    /* the gap split is genuinely interior: 1 <= m <= d-1. */
    AIC_CHECK_MSG(w.m >= 1 && w.m <= d - 1, "T1(d=%ld): m=%ld not interior",
                  (long) d, (long) w.m);

    arb_clear(IPn); arb_clear(Pn); arb_clear(delta);
    acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* Block conditional expectation on C^d, blocks [0,m),[m,d): A = M_m (+) M_{d-m},
 * dim_A = m^2+(d-m)^2. The finder returns a valid nontrivial projection in A; its
 * rank may be 1..m or 1..(d-m) (for d=4,m=2 it is rank-1) -- all are valid here,
 * so the test asserts nontriviality/in-A/delta=0, NOT a specific rank. */
static void t1_block(slong d, slong m, slong prec)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, d, m);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong expect_dA = m * m + (d - m) * (d - m);
    AIC_CHECK_MSG(ae.A.dim_A == expect_dA, "T1 block(d=%ld,m=%ld): dim_A=%ld != %ld",
                  (long) d, (long) m, (long) ae.A.dim_A, (long) expect_dA);

    acb_mat_t P;
    arb_t delta, Pn, IPn;
    acb_mat_init(P, d, d);
    arb_init(delta); arb_init(Pn); arb_init(IPn);
    aic_projection_nontrivial(P, delta, Pn, IPn, &ae.A, NULL, prec);
    double del = dd(delta);
    printf("T1 block d=%ld m=%ld: delta=%.3e ||P||=%.6f ||I-P||=%.6f\n",
           (long) d, (long) m, del, dd(Pn), dd(IPn));
    AIC_CHECK_MSG(del < 1e-9, "T1 block: delta=%.3e not machine-zero at eta=0", del);
    check_nontrivial("T1 block", P, Pn, IPn, &ae.A, 1e-9, prec);

    arb_clear(IPn); arb_clear(Pn); arb_clear(delta);
    acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

static void test_t1(void)
{
    t1_identity(2, 256);
    t1_identity(3, 256);
    t1_block(4, 2, 256);
}

/* =========================== T2: eta>0 (oblique) ========================= */
static void test_t2(void)
{
    const slong prec = 256;
    /* make_mixconj(5,3,0.02): genuinely eta>0, NON-COMMUTATIVE, OBLIQUE. */
    aic_ucp_kraus phi;
    make_mixconj(&phi, 5, 3, 0.02, prec);
    double eta = eta_proxy(&phi, prec);
    AIC_CHECK_MSG(eta > 1e-4, "T2: eta_proxy=%.3e not genuinely > 0 (vacuous)", eta);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong n = ae.A.n;

    acb_mat_t P;
    arb_t delta, Pn, IPn;
    acb_mat_init(P, n, n);
    arb_init(delta); arb_init(Pn); arb_init(IPn);
    aic_projection_witness w;
    aic_projection_nontrivial(P, delta, Pn, IPn, &ae.A, &w, prec);

    double del = dd(delta);
    double C = del / eta;
    printf("T2 mixconj(5,3,t=0.02): eta=%.4e delta=%.4e C=delta/eta=%.3f "
           "||P||=%.6f ||I-P||=%.6f gap=%.3f\n", eta, del, C, dd(Pn), dd(IPn),
           w.gap);
    /* delta = ||P*P-P|| <= C*eta with C a bounded constant (report it). C < 1
     * catches the O(1)-defect regression (Frobenius projection had C ~ 15). */
    AIC_CHECK_MSG(C < 1.0, "T2: C=delta/eta=%.3f exceeds 1 -- delta not O(eta) "
                  "(Frobenius-projection regression?)", C);
    AIC_CHECK_MSG(del > 1e-12, "T2: delta=%.3e ~0 at eta>0 -- star defect vacuous "
                  "(is the star == plain? Pi_A not load-bearing?)", del);
    /* nontrivial + in A; at eta>0 the residual is O(eta), not machine-zero. */
    check_nontrivial("T2 mixconj", P, Pn, IPn, &ae.A, 2.0 * eta + 1e-9, prec);

    /* NON-VACUITY: prove the Phi_tilde projection step (project_into_A) is
     * LOAD-BEARING — i.e. the suite is NOT blind to a skip-Pi_A mutation. Two
     * complementary teeth, both needed since the finder ranks H by GAP (FIX 1,
     * aic-3qv) and for this fixture the largest-gap H gives a P_amb already nearly
     * in A (raw defect ~3.7e-3), making a chosen-H-only obliqueness test toothless:
     *
     *   (A) Phi_tilde works: scan all H_k, build each ambient projector, take the
     *       MOST-OBLIQUE one (largest raw star defect, genuinely outside A), and
     *       show Phi_tilde = (.)*I drops its defect by >=10x and moves it. Proves
     *       the star-with-the-unit operation (.tex:2211, Phi_tilde(X)=X*I) reduces
     *       the defect — fixture-robust, independent of the finder's H choice.
     *   (B) the finder USES Phi_tilde: the finder's output P must equal
     *       Phi_tilde(P_amb) for the finder's OWN chosen H (project_into_A applied),
     *       NOT the raw P_amb. Under a skip-Pi_A mutation the finder returns the raw
     *       P_amb, which differs from Phi_tilde(P_amb) (it moved it), so this check
     *       goes RED. This catches the project_into_A mutation regardless of how
     *       oblique the chosen H is. */
    {
        slong dn = n;
        acb_mat_t Imat;
        acb_mat_init(Imat, dn, dn);
        acb_mat_one(Imat);
        double best_raw = -1.0, best_proj = 0.0, best_move = 0.0;
        slong best_k = -1;
        double chosen_pp = -1.0;   /* ||P - Phi_tilde(P_amb_chosen)|| (tooth B) */
        for (slong kk = 0; kk < ae.A.dim_A; kk++) {
            acb_mat_t H, Y, sgnY, Pamb, Bd;
            acb_mat_init(H, dn, dn); acb_mat_init(Y, dn, dn);
            acb_mat_init(sgnY, dn, dn); acb_mat_init(Pamb, dn, dn);
            acb_mat_init(Bd, dn, dn);
            /* H_k = (B_k + B_k^dag)/2 */
            acb_mat_conjugate_transpose(Bd, ae.A.B[kk]);
            acb_mat_add(H, ae.A.B[kk], Bd, prec);
            acb_mat_scalar_mul_2exp_si(H, H, -1);
            /* largest interior gap of H_k (double path, mirrors the finder). */
            double evk[64];
            {
                double _Complex Harr[64 * 64];
                aic_latd_from_acb_mat(Harr, H);
                aic_latd_eig_hermitian(evk, NULL, Harr, dn);
            }
            double bg = -1.0; slong bi = -1;
            for (slong i = 0; i < dn - 1; i++) {
                double g = evk[i + 1] - evk[i];
                if (g > bg) { bg = g; bi = i; }
            }
            double tk = 0.5 * (evk[bi] + evk[bi + 1]);
            double lmn = evk[0], lmx = evk[dn - 1];
            double s = 1.0 / fmax(tk - lmn, lmx - tk);
            arb_t sc; arb_init(sc);
            arb_set_d(sc, tk);
            acb_mat_scalar_mul_arb(Y, Imat, sc, prec);   /* tI */
            acb_mat_sub(Y, H, Y, prec);                  /* H - tI */
            arb_set_d(sc, s);
            acb_mat_scalar_mul_arb(Y, Y, sc, prec);      /* s(H - tI) */
            arb_clear(sc);
            aic_sgn(sgnY, Y, prec);
            acb_mat_add(Pamb, Imat, sgnY, prec);
            acb_mat_scalar_mul_2exp_si(Pamb, Pamb, -1);  /* P_amb = (I+sgn Y)/2 */
            /* Phi_tilde(P_amb) = P_amb * I (the project_into_A step). */
            acb_mat_t Pp;
            acb_mat_init(Pp, dn, dn);
            aic_ecstar_star(Pp, &ae.A, Pamb, Imat, prec);
            /* raw star defect ||P_amb*P_amb - P_amb|| (skip-Pi_A mutation value). */
            acb_mat_t PaPa, Dsd;
            acb_mat_init(PaPa, dn, dn); acb_mat_init(Dsd, dn, dn);
            aic_ecstar_star(PaPa, &ae.A, Pamb, Pamb, prec);
            acb_mat_sub(Dsd, PaPa, Pamb, prec);
            double raw = mid_opnorm_d(Dsd, prec);
            /* tooth B: for the finder's chosen H, ||finder's P - Phi_tilde(P_amb)||. */
            if (kk == w.k_chosen) {
                acb_mat_t Dchk;
                acb_mat_init(Dchk, dn, dn);
                acb_mat_sub(Dchk, P, Pp, prec);
                chosen_pp = mid_opnorm_d(Dchk, prec);
                acb_mat_clear(Dchk);
            }
            if (raw > best_raw) {
                best_raw = raw; best_k = kk;
                /* tooth A: defect of Phi_tilde(P_amb) and the move it produced. */
                acb_mat_t PpPp, Dp, Dmove;
                acb_mat_init(PpPp, dn, dn);
                acb_mat_init(Dp, dn, dn); acb_mat_init(Dmove, dn, dn);
                aic_ecstar_star(PpPp, &ae.A, Pp, Pp, prec);
                acb_mat_sub(Dp, PpPp, Pp, prec);
                best_proj = mid_opnorm_d(Dp, prec);
                acb_mat_sub(Dmove, Pp, Pamb, prec);
                best_move = mid_opnorm_d(Dmove, prec);
                acb_mat_clear(Dmove); acb_mat_clear(Dp); acb_mat_clear(PpPp);
            }
            acb_mat_clear(Dsd); acb_mat_clear(PaPa); acb_mat_clear(Pp);
            acb_mat_clear(Bd); acb_mat_clear(Pamb); acb_mat_clear(sgnY);
            acb_mat_clear(Y); acb_mat_clear(H);
        }
        printf("T2 non-vacuity (A): most-oblique H_k=%ld raw ||P_amb*P_amb-P_amb||"
               "=%.4e (skip-Pi_A value) -> Phi_tilde defect=%.4e (Pi_A reduces "
               "%.0fx); moved ||.||=%.4e. (B) ||P - Phi_tilde(P_amb_chosen)||=%.2e\n",
               (long) best_k, best_raw, best_proj,
               best_proj > 0 ? best_raw / best_proj : 0.0, best_move, chosen_pp);
        /* (A-i) the skip-Pi_A mutation value is O(1) (fixture genuinely oblique). */
        AIC_CHECK_MSG(best_raw > 0.05, "T2: most-oblique raw P_amb star defect=%.3e "
                      "too small -- no oblique H_k, fixture not exercising Pi_A", best_raw);
        /* (A-ii) Phi_tilde reduces the defect by >=10x. */
        AIC_CHECK_MSG(best_proj < 0.1 * best_raw, "T2: Phi_tilde did not reduce the "
                      "star defect (raw=%.3e -> proj=%.3e, factor %.1f < 10) -- Pi_A "
                      "not load-bearing", best_raw, best_proj, best_raw / best_proj);
        /* (A-iii) Phi_tilde MOVED the oblique P_amb. */
        AIC_CHECK_MSG(best_move > 1e-3, "T2: Phi_tilde did not move the oblique "
                      "P_amb (||.||=%.3e ~0) -- skip-Pi_A mutation invisible", best_move);
        /* (B) the FINDER applied Phi_tilde: its P equals Phi_tilde(P_amb_chosen),
         * NOT raw P_amb. Tolerance ties to the chosen-H projection accuracy; tight
         * (< 1e-7) so a skip-Pi_A mutation (P=raw P_amb, differs by ~5e-3) is RED. */
        AIC_CHECK_MSG(chosen_pp >= 0.0 && chosen_pp < 1e-7, "T2: finder's P does not "
                      "match Phi_tilde(P_amb_chosen) (||.||=%.3e) -- project_into_A "
                      "not applied (skip-Pi_A mutation?)", chosen_pp);
        acb_mat_clear(Imat);
    }

    arb_clear(IPn); arb_clear(Pn); arb_clear(delta);
    acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* =========================== T3: universality canary ===================== */
static void test_t3(void)
{
    const slong prec = 256;
    printf("T3 universality canary (delta should NOT grow with d):\n");
    double maxdel = 0.0;
    for (slong d = 2; d <= 4; d++) {
        aic_ucp_kraus phi;
        make_identity(&phi, d);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        acb_mat_t P;
        arb_t delta, Pn, IPn;
        acb_mat_init(P, d, d);
        arb_init(delta); arb_init(Pn); arb_init(IPn);
        aic_projection_nontrivial(P, delta, Pn, IPn, &ae.A, NULL, prec);
        double del = dd(delta);
        if (del > maxdel) maxdel = del;
        printf("  d=%ld: delta=%.3e ||P||=%.6f ||I-P||=%.6f\n", (long) d, del,
               dd(Pn), dd(IPn));
        arb_clear(IPn); arb_clear(Pn); arb_clear(delta);
        acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    AIC_CHECK_MSG(maxdel < 1e-9, "T3: delta grows with d (max=%.3e) -- NOT "
                  "dimension-independent (tex:484)", maxdel);

    /* eta>0 sweep across 5 dims (FIX 3, the highest-value test, tex:484):
     * mixconj(d,m,0.02) for (d,m) = (4,2),(5,3),(6,3),(7,4),(8,4) spans algebra
     * dims dim_A = 4,9,9,16,16 and ambient n = 4..8. C = delta/eta must NOT grow
     * with d. An ABSOLUTE-magnitude gate at 2 points (the old cmax<1 at d=4,5)
     * would miss a slowly d-growing constant C ~ c0*d that stays < 1 until d ~ 20
     * (the naive Haar/cohomology route tex:484 warns against); so we ALSO assert a
     * SLOPE bound C(largest d)/C(smallest d) < FLAT_RATIO. The d-linear failure
     * mode over d=4..8 gives a slope ~ d_max/d_min = 8/4 = 2.0, so FLAT_RATIO=1.5
     * catches it (mutation-proved by injecting C=0.05*d) while the GENUINE data
     * (C fast-DECREASING, observed endpoint slope ~0.05) clears it with ~28x
     * margin -- teeth against d-growth without flakiness. */
    printf("T3 eta>0 delta/eta across 5 dims (should be FLAT, not d-growing):\n");
    const slong ds[] = {4, 5, 6, 7, 8};
    const slong ms[] = {2, 3, 3, 4, 4};
    const int NSW = 5;
    double cmin = 1e30, cmax = 0.0;
    double c_at_dmin = -1.0, c_at_dmax = -1.0;
    slong d_min = ds[0], d_max = ds[0];
    for (int i = 0; i < NSW; i++) {
        if (ds[i] < d_min) d_min = ds[i];
        if (ds[i] > d_max) d_max = ds[i];
    }
    for (int i = 0; i < NSW; i++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, ds[i], ms[i], 0.02, prec);
        double eta = eta_proxy(&phi, prec);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        slong n = ae.A.n;
        acb_mat_t P;
        arb_t delta, Pn, IPn;
        acb_mat_init(P, n, n);
        arb_init(delta); arb_init(Pn); arb_init(IPn);
        aic_projection_nontrivial(P, delta, Pn, IPn, &ae.A, NULL, prec);
        double C = dd(delta) / eta;
        /* MUTATION HOOK (compile-time): REPLACE C with the adversarial d-linear
         * constant C ~ 0.05*d (the exact failure mode FIX 3 guards against: a
         * slowly d-growing constant that passes an absolute cap until d ~ 20).
         * Setting (not scaling) C models it faithfully — scaling the genuine,
         * fast-DECREASING observed C would be masked by that decrease. Confirms
         * the slope assert below goes RED. Off by default. */
#ifdef AIC_T3_INJECT_DGROWTH
        C = 0.05 * (double) ds[i];
#endif
        if (C < cmin) cmin = C;
        if (C > cmax) cmax = C;
        if (ds[i] == d_min) c_at_dmin = C;
        if (ds[i] == d_max) c_at_dmax = C;
        printf("  d=%ld m=%ld dim_A=%ld: eta=%.4e delta=%.4e C=%.4e\n", (long) ds[i],
               (long) ms[i], (long) ae.A.dim_A, eta, dd(delta), C);
        arb_clear(IPn); arb_clear(Pn); arb_clear(delta);
        acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    double slope = (c_at_dmin > 0.0) ? c_at_dmax / c_at_dmin : 1e30;
    printf("  -> eta>0 C in [%.4e, %.4e]; slope C(d=%ld)/C(d=%ld)=%.3f\n",
           cmin, cmax, (long) d_max, (long) d_min, slope);
    /* (1) absolute cap: C < 1 catches the O(1)-defect regression (the Frobenius-
     * projection bug had C ~ 15-60; the Phi_tilde projection gives C ~ 1e-3..6e-2). */
    AIC_CHECK_MSG(cmax < 1.0, "T3: eta>0 C=%.4e exceeds 1 -- delta not O(eta) "
                  "(Frobenius-projection regression? should use Phi_tilde)", cmax);
    /* (2) flatness/slope: C must NOT grow ~d. FLAT_RATIO=1.5 catches the d-linear
     * failure mode (slope ~2.0 over d=4..8) while clearing the genuine fast-
     * decreasing C (observed slope ~0.05) with ~28x margin. */
    AIC_CHECK_MSG(slope < 1.5, "T3: eta>0 C grows with d (C(d=%ld)/C(d=%ld)=%.3f "
                  ">= 1.5) -- NOT dimension-independent, the naive Haar/cohomology "
                  "route the paper warns against (tex:484)", (long) d_max,
                  (long) d_min, slope);
}

/* =========================== T4: fail-loud (message-checked) ============ */
/* Run `child` in a forked process with stderr captured to a temp file; after it
 * exits, assert (1) it died by SIGABRT and (2) the captured stderr CONTAINS
 * `expect`. The message-grep is what gives the fail-loud tests teeth: deleting
 * the INTENDED guard (and falling through to a different abort, or none) changes
 * the message, turning the test RED. SIGABRT-only would stay green if any other
 * guard happened to fire. The temp file (stderr is unbuffered, and glibc abort()
 * flushes) survives the child's abort. */
static void expect_abort_with_msg(void (*child)(void), const char *expect,
                                  const char *name)
{
    char tmpl[] = "/tmp/aic_proj_failXXXXXX";
    int fd = mkstemp(tmpl);
    AIC_CHECK_MSG(fd >= 0, "%s: mkstemp failed", name);
    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "%s: fork failed", name);
    if (pid == 0) {
        dup2(fd, STDERR_FILENO);   /* child's stderr -> temp file */
        close(fd);
        child();
        _exit(0);   /* reached only if NO abort fired (a test failure) */
    }
    int status = 0;
    waitpid(pid, &status, 0);
    close(fd);
    int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
    /* read the captured stderr. */
    char buf[4096];
    size_t got = 0;
    FILE *f = fopen(tmpl, "r");
    if (f) { got = fread(buf, 1, sizeof(buf) - 1, f); fclose(f); }
    buf[got] = '\0';
    unlink(tmpl);
    int has_msg = (strstr(buf, expect) != NULL);
    printf("%s: child %s; stderr %s expected substring \"%s\"\n", name,
           aborted ? "SIGABRT" : "did NOT abort (WRONG)",
           has_msg ? "CONTAINS" : "MISSING", expect);
    AIC_CHECK_MSG(aborted, "%s: child did not abort via SIGABRT", name);
    AIC_CHECK_MSG(has_msg, "%s: captured stderr does NOT contain the expected "
                  "guard message \"%s\" (got: %.200s) -- the INTENDED guard did "
                  "not fire (deleted/wrong guard?)", name, expect, buf);
}

/* trace_replace(k): A = scalars, dim_A = 1. lem_nontriv_projection requires
 * dim A > 1, so the finder must ABORT via the dim-1 guard. */
static void t4_dim1_child(void)
{
    const slong prec = 128;
    aic_ucp_kraus phi;
    make_trace_replace(&phi, 3);   /* A = C*I, dim_A = 1 */
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    AIC_CHECK_MSG(ae.A.dim_A == 1, "t4: trace_replace dim_A=%ld != 1",
                  (long) ae.A.dim_A);
    acb_mat_t P;
    arb_t delta, Pn, IPn;
    acb_mat_init(P, ae.A.n, ae.A.n);
    arb_init(delta); arb_init(Pn); arb_init(IPn);
    aic_projection_nontrivial(P, delta, Pn, IPn, &ae.A, NULL, prec);   /* MUST abort */
    _exit(0);
}

/* aic-3qv no-interior-gap abort: feed aic_projection_gap a CONSTANT-DIAGONAL
 * Hermitian H (all eigenvalues equal) -> no positive interior gap -> the
 * degenerate-spectrum abort (the central documented escalation). After FIX 1
 * this fires ONLY for a genuinely degenerate spectrum, so a constant-diagonal H
 * is the right witness. */
static void t4_nogap_child(void)
{
    const slong prec = 128, nn = 5;
    acb_mat_t H;
    acb_mat_init(H, nn, nn);
    for (slong i = 0; i < nn; i++) acb_set_d(acb_mat_entry(H, i, i), 2.5);  /* 2.5*I */
    double t = 0, g = 0, lmin = 0, lmax = 0;
    slong m = 0;
    aic_projection_gap(&t, &g, &lmin, &lmax, &m, H, prec);   /* MUST abort (no gap) */
    acb_mat_clear(H);
    _exit(0);
}

static void test_t4(void)
{
    /* (a) dim A == 1: the dim-1 guard, message-checked. */
    expect_abort_with_msg(t4_dim1_child, "dim A = 1 <= 1", "T4(a) dim-1");
    /* (b) aic-3qv: no positive interior gap (degenerate spectrum), message-checked. */
    expect_abort_with_msg(t4_nogap_child, "NO positive interior spectral gap",
                          "T4(b) no-gap aic-3qv");
}

/* =========================== T5: double vs arb@53 ======================= */
static void test_t5(void)
{
    aic_ucp_kraus phi;
    make_identity(&phi, 3);
    double dvals[2][3];   /* [prec_idx] = {delta, ||P||, ||I-P||} */
    slong precs[2] = {53, 256};
    for (int pi = 0; pi < 2; pi++) {
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, precs[pi]);
        acb_mat_t P;
        arb_t delta, Pn, IPn;
        acb_mat_init(P, 3, 3);
        arb_init(delta); arb_init(Pn); arb_init(IPn);
        aic_projection_nontrivial(P, delta, Pn, IPn, &ae.A, NULL, precs[pi]);
        dvals[pi][0] = dd(delta);
        dvals[pi][1] = dd(Pn);
        dvals[pi][2] = dd(IPn);
        arb_clear(IPn); arb_clear(Pn); arb_clear(delta);
        acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
    }
    printf("T5 prec=53 vs 256: delta %.3e/%.3e ||P|| %.6f/%.6f ||I-P|| %.6f/%.6f\n",
           dvals[0][0], dvals[1][0], dvals[0][1], dvals[1][1], dvals[0][2],
           dvals[1][2]);
    AIC_CHECK_MSG(fabs(dvals[0][1] - dvals[1][1]) < 1e-9,
                  "T5: ||P|| differs prec=53 vs 256 (%.6f vs %.6f)",
                  dvals[0][1], dvals[1][1]);
    AIC_CHECK_MSG(fabs(dvals[0][2] - dvals[1][2]) < 1e-9,
                  "T5: ||I-P|| differs prec=53 vs 256");
    aic_ucp_kraus_clear(&phi);
}

/* =========================== T6: asymmetric spectrum ==================== */
/* FIX 4(a): the threshold/m consistency on an ASYMMETRIC (non-symmetric-about-0)
 * spectrum, the case the review verified by hand but the suite did not cover.
 * Feed aic_projection_gap a diagonal Hermitian H with chosen eigenvalues and
 * assert (i) m == count(lam >= t) for the largest-gap threshold t, and (ii) the
 * AMBIENT spectral projector P_amb = (I+sgn(s(H-tI)))/2 has rank == m. Guards a
 * future off-by-one/sign in m_out or t. Two asymmetric spectra. */
static slong count_ge(const double *ev, slong n, double t)
{
    slong c = 0;
    for (slong i = 0; i < n; i++) if (ev[i] >= t) c++;
    return c;
}
/* rank of a Hermitian acb_mat (count of eigenvalues > 0.5; P_amb's are ~0 or ~1). */
static slong herm_rank_half(const acb_mat_t Pm, slong n)
{
    double _Complex Parr[64 * 64];
    double ev[64];
    aic_latd_from_acb_mat(Parr, Pm);
    aic_latd_eig_hermitian(ev, NULL, Parr, n);
    slong r = 0;
    for (slong i = 0; i < n; i++) if (ev[i] > 0.5) r++;
    return r;
}
static void t6_one(const double *spec, slong n, const char *name, slong prec)
{
    acb_mat_t H;
    acb_mat_init(H, n, n);
    for (slong i = 0; i < n; i++) acb_set_d(acb_mat_entry(H, i, i), spec[i]);
    double t = 0, g = 0, lmin = 0, lmax = 0;
    slong m = 0;
    aic_projection_gap(&t, &g, &lmin, &lmax, &m, H, prec);
    slong cge = count_ge(spec, n, t);
    /* build P_amb (ambient sgn route, exactly as the finder's ambient_projector). */
    double s = 1.0 / fmax(t - lmin, lmax - t);
    acb_mat_t Y, sgnY, Pamb, Imat;
    acb_mat_init(Y, n, n); acb_mat_init(sgnY, n, n);
    acb_mat_init(Pamb, n, n); acb_mat_init(Imat, n, n);
    acb_mat_one(Imat);
    arb_t sc; arb_init(sc);
    arb_set_d(sc, t);
    acb_mat_scalar_mul_arb(Y, Imat, sc, prec);
    acb_mat_sub(Y, H, Y, prec);
    arb_set_d(sc, s);
    acb_mat_scalar_mul_arb(Y, Y, sc, prec);
    arb_clear(sc);
    aic_sgn(sgnY, Y, prec);
    acb_mat_add(Pamb, Imat, sgnY, prec);
    acb_mat_scalar_mul_2exp_si(Pamb, Pamb, -1);
    slong rk = herm_rank_half(Pamb, n);
    printf("T6 %s: t=%.3f gap=%.3f m=%ld count(lam>=t)=%ld rank(P_amb)=%ld\n",
           name, t, g, (long) m, (long) cge, (long) rk);
    AIC_CHECK_MSG(m == cge, "T6 %s: m=%ld != count(lam>=t)=%ld (m_out off-by-one"
                  "/sign on asymmetric spectrum?)", name, (long) m, (long) cge);
    AIC_CHECK_MSG(rk == m, "T6 %s: rank(P_amb)=%ld != m=%ld (threshold/sgn "
                  "mismatch on asymmetric spectrum?)", name, (long) rk, (long) m);
    acb_mat_clear(Imat); acb_mat_clear(Pamb); acb_mat_clear(sgnY);
    acb_mat_clear(Y); acb_mat_clear(H);
}
static void test_t6(void)
{
    const slong prec = 256;
    /* {-5,0,1,2}: gaps 5,1,1 -> largest [-5,0], t=-2.5, m=count(>=-2.5)=3. */
    const double s1[] = {-5.0, 0.0, 1.0, 2.0};
    t6_one(s1, 4, "{-5,0,1,2}", prec);
    /* {-2,-1,3}: gaps 1,4 -> largest [-1,3], t=1.0, m=count(>=1)=1. */
    const double s2[] = {-2.0, -1.0, 3.0};
    t6_one(s2, 3, "{-2,-1,3}", prec);
}

int main(void)
{
    test_t1();
    test_t2();
    test_t3();
    test_t4();
    test_t5();
    test_t6();
    aic_test_report("test_projection");
    printf("OK test_projection\n");
    return 0;
}
