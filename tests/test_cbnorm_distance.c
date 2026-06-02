/* test_cbnorm_distance.c — hostile cross-checks for the eig-free certified
 * bracket on the GENERAL cb/diamond distance ||Phi - Psi||_cb between two UCP
 * self-maps (aic_cbnorm_eigfree_distance, src/aic_cbnorm_distance.c). The bracket
 *     [ ||J||_F / n ,  2 * ||J||_F ]   with J = Choi(Phi) - Choi(Psi)  (n^2 x n^2)
 * must rigorously contain ||Phi-Psi||_cb (= ||.||_diamond for the Hermiticity-
 * preserving Lambda = Phi - Psi), approximate_algebras.tex:347-354.
 *
 * Assertion families (each NON-VACUOUS — would go RED under a real J/bound bug):
 *   (A) Phi==Psi oracle: J=0 -> bracket [0,0] (sanity floor, blind to most bugs).
 *   (B) Symmetry: J -> -J leaves ||J||_F (hence the bracket) invariant.
 *   (C) Linearity oracle (LOAD-BEARING): D_p - D_q = (q-p) C, C(X)=X-(trX/d)I a
 *       FIXED Herm-pres map, so BOTH endpoints scale LINEARLY in |q-p|. Wrong J
 *       assembly breaks the ratios.
 *   (D) CLOSED-FORM containment (THE correctness check): identity Phi_I vs the
 *       qubit unitary Phi_W, W=diag(1,e^{i theta}), has EXACT cb/diamond distance
 *           ||Phi_I - Phi_W||_cb = 2 sin(theta/2)   (theta in [0,pi]).
 *       Derivation: spec(W)={1,e^{i theta}}; the distance from 0 to conv(spec) is
 *       nu = cos(theta/2); diamond dist = 2 sqrt(1 - nu^2) = 2 sin(theta/2). The
 *       bracket [sqrt2 sin, 4 sqrt2 sin]*sin(theta/2) must contain 2 sin(theta/2)
 *       with margins (lo/true=0.707, hi/true=2.83). A wrong dimension power or a
 *       sign error in J would push lo above true -> RED.
 *   (E) Precision ladder: family (D) at prec=256 still contains, and the 53 vs
 *       256 endpoints agree (the bracket is determined by ||J||_F, prec-stable).
 *   (F) Fail-loud teeth (Rule 4): different-n maps -> the dim-guard assert ABORTS
 *       (shared aic_watchdog, bead aic-8de; asserts are live, build strips -DNDEBUG).
 */
#include <math.h>
#include <stdio.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_cbnorm.h"
#include "aic/aic_channels.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "aic_watchdog.h"

#define PREC 53

/* --- helpers (mirrored from test_cbnorm.c) -------------------------------- */

static int arb_close(const arb_t a, const arb_t b, const arb_t tol)
{
    arb_t d;
    arb_init(d);
    arb_sub(d, a, b, 64);
    arb_abs(d, d);
    int ok = arb_le(d, tol);
    arb_clear(d);
    return ok;
}
static double arb_lo_d(const arb_t x)
{
    arf_t t;
    arf_init(t);
    arb_get_lbound_arf(t, x, 64);
    double v = arf_get_d(t, ARF_RND_FLOOR);
    arf_clear(t);
    return v;
}
static double arb_hi_d(const arb_t x)
{
    arf_t t;
    arf_init(t);
    arb_get_ubound_arf(t, x, 64);
    double v = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
    return v;
}

/* Build W = diag(1, e^{i theta}) as a 2x2 unitary (off-diagonals 0). */
static void build_W(acb_mat_t W, double theta)
{
    acb_set_si(acb_mat_entry(W, 0, 0), 1);
    acb_set_d_d(acb_mat_entry(W, 1, 1), cos(theta), sin(theta));
}

/* --- (A) Phi==Psi oracle: J=0 -> bracket [0,0] --------------------------- */
static void test_self_distance_zero(void)
{
    printf("(A) Phi==Psi: ||Phi-Phi||_cb bracket collapses to [0,0]\n");

    /* depolarizing self-distance */
    aic_ucp_kraus dep;
    aic_channel_depolarizing(&dep, 3, 0.4, PREC);
    arb_t lo, hi;
    arb_init(lo); arb_init(hi);
    aic_cbnorm_eigfree_distance(lo, hi, &dep, &dep, PREC);
    AIC_CHECK_MSG(arb_hi_d(hi) < 1e-9, "(A) depol self: hi=%.3e !< 1e-9",
                  arb_hi_d(hi));
    arb_clear(hi); arb_clear(lo);
    aic_ucp_kraus_clear(&dep);

    /* unitary-channel self-distance (W = diag(1, e^{i*0.7})) */
    acb_mat_t W; acb_mat_init(W, 2, 2);
    build_W(W, 0.7);
    aic_ucp_kraus phiW;
    double pr1[1] = {1.0};
    const acb_mat_t *Warr = (const acb_mat_t *) W;
    aic_channel_mix_unitaries(&phiW, Warr, pr1, 1, 2, PREC);
    arb_t lo2, hi2;
    arb_init(lo2); arb_init(hi2);
    aic_cbnorm_eigfree_distance(lo2, hi2, &phiW, &phiW, PREC);
    AIC_CHECK_MSG(arb_hi_d(hi2) < 1e-9, "(A) unitary self: hi=%.3e !< 1e-9",
                  arb_hi_d(hi2));
    arb_clear(hi2); arb_clear(lo2);
    aic_ucp_kraus_clear(&phiW);
    acb_mat_clear(W);
}

/* --- (B) symmetry: (phi,psi) bracket == (psi,phi) bracket exactly --------- */
static void test_symmetry(void)
{
    printf("(B) ||Phi-Psi||_cb == ||Psi-Phi||_cb (J -> -J, Frobenius invariant)\n");
    aic_ucp_kraus a, b;
    aic_channel_depolarizing(&a, 2, 0.3, PREC);
    aic_channel_depolarizing(&b, 2, 0.6, PREC);
    arb_t lo_ab, hi_ab, lo_ba, hi_ba, tol;
    arb_init(lo_ab); arb_init(hi_ab);
    arb_init(lo_ba); arb_init(hi_ba);
    arb_init(tol); arb_set_d(tol, 1e-14);
    aic_cbnorm_eigfree_distance(lo_ab, hi_ab, &a, &b, PREC);
    aic_cbnorm_eigfree_distance(lo_ba, hi_ba, &b, &a, PREC);
    /* non-vacuous: confirm the brackets are not trivially [0,0] (p!=q => J!=0) */
    AIC_CHECK_MSG(arb_lo_d(lo_ab) > 1e-6,
                  "(B) distinct channels but lo~0 (bracket vacuous): lo=%.3e",
                  arb_lo_d(lo_ab));
    AIC_CHECK_MSG(arb_close(lo_ab, lo_ba, tol), "(B) lo asymmetric");
    AIC_CHECK_MSG(arb_close(hi_ab, hi_ba, tol), "(B) hi asymmetric");
    arb_clear(tol);
    arb_clear(hi_ba); arb_clear(lo_ba);
    arb_clear(hi_ab); arb_clear(lo_ab);
    aic_ucp_kraus_clear(&b);
    aic_ucp_kraus_clear(&a);
}

/* --- (C) linearity oracle: D_p - D_q endpoints scale as |q-p| (LOAD-BEARING)*/
static void test_linearity(void)
{
    printf("(C) linearity: D_p-D_q endpoints scale LINEARLY in |q-p|\n");
    const double p = 0.1;
    const double qs[3] = {0.4, 0.7, 0.55};      /* |q-p| = 0.3, 0.6, 0.45 */
    double los[3], his[3];
    aic_ucp_kraus Dp;
    aic_channel_depolarizing(&Dp, 2, p, PREC);
    for (int i = 0; i < 3; i++) {
        aic_ucp_kraus Dq;
        aic_channel_depolarizing(&Dq, 2, qs[i], PREC);
        arb_t lo, hi;
        arb_init(lo); arb_init(hi);
        aic_cbnorm_eigfree_distance(lo, hi, &Dp, &Dq, PREC);
        /* endpoints are exact ||J||_F/n and 2||J||_F; midpoint == value here */
        los[i] = 0.5 * (arb_lo_d(lo) + arb_hi_d(lo));
        his[i] = 0.5 * (arb_lo_d(hi) + arb_hi_d(hi));
        arb_clear(hi); arb_clear(lo);
        aic_ucp_kraus_clear(&Dq);
    }
    /* non-vacuous: q=0.4 endpoints strictly positive */
    AIC_CHECK_MSG(los[0] > 1e-6 && his[0] > 1e-6,
                  "(C) base bracket vacuous: lo=%.3e hi=%.3e", los[0], his[0]);
    /* ratio q2/q1 = 0.6/0.3 = 2.0 on BOTH endpoints, to ~1e-12 */
    AIC_CHECK_MSG(fabs(his[1] / his[0] - 2.0) < 1e-12,
                  "(C) hi ratio %.15f != 2.0", his[1] / his[0]);
    AIC_CHECK_MSG(fabs(los[1] / los[0] - 2.0) < 1e-12,
                  "(C) lo ratio %.15f != 2.0", los[1] / los[0]);
    /* third point confirms LINEARITY: |q-p|=0.45 -> ratio 0.45/0.3 = 1.5 */
    AIC_CHECK_MSG(fabs(his[2] / his[0] - 1.5) < 1e-12,
                  "(C) hi ratio3 %.15f != 1.5", his[2] / his[0]);
    AIC_CHECK_MSG(fabs(los[2] / los[0] - 1.5) < 1e-12,
                  "(C) lo ratio3 %.15f != 1.5", los[2] / los[0]);
    printf("    los=[%.6e %.6e %.6e] his=[%.6e %.6e %.6e]\n",
           los[0], los[1], los[2], his[0], his[1], his[2]);
    aic_ucp_kraus_clear(&Dp);
}

/* --- (D) closed-form containment oracle: 2 sin(theta/2) in [lo,hi] -------- */
static void run_D_containment(slong prec, double *lo_out, double *hi_out,
                              const double *thetas, int ntheta)
{
    /* identity channel Phi_I via U = I_2 */
    acb_mat_t I2; acb_mat_init(I2, 2, 2);
    acb_set_si(acb_mat_entry(I2, 0, 0), 1);
    acb_set_si(acb_mat_entry(I2, 1, 1), 1);
    double pr1[1] = {1.0};
    aic_ucp_kraus phiI;
    aic_channel_mix_unitaries(&phiI, (const acb_mat_t *) I2, pr1, 1, 2, prec);

    for (int k = 0; k < ntheta; k++) {
        double th = thetas[k];
        acb_mat_t W; acb_mat_init(W, 2, 2);
        build_W(W, th);
        aic_ucp_kraus phiW;
        aic_channel_mix_unitaries(&phiW, (const acb_mat_t *) W, pr1, 1, 2, prec);

        arb_t lo, hi;
        arb_init(lo); arb_init(hi);
        aic_cbnorm_eigfree_distance(lo, hi, &phiI, &phiW, prec);
        double lo_lb = arb_lo_d(lo), hi_ub = arb_hi_d(hi);
        double truev = 2.0 * sin(th / 2.0);
        AIC_CHECK_MSG(lo_lb <= truev,
                      "(D) prec=%ld theta=%.3f: lo=%.6e > true=%.6e (lower bound "
                      "EXCEEDS the exact cb-distance!)", prec, th, lo_lb, truev);
        AIC_CHECK_MSG(truev <= hi_ub,
                      "(D) prec=%ld theta=%.3f: hi=%.6e < true=%.6e (upper bound "
                      "below the exact cb-distance!)", prec, th, hi_ub, truev);
        /* ENDPOINT PINS (pin the 1/n constant and the factor 2 — catches a
         * wrong-n / wrong-constant mutation that mere containment misses).
         * Here n=2, ||J||_F = 2 sqrt2 sin(th/2), so lo = ||J||_F/2 = sqrt2 sin,
         * hi = 2 ||J||_F = 4 sqrt2 sin. */
        double s = sin(th / 2.0);
        double lo_mid = 0.5 * (arb_lo_d(lo) + arb_hi_d(lo));
        double hi_mid = 0.5 * (arb_lo_d(hi) + arb_hi_d(hi));
        AIC_CHECK_MSG(fabs(lo_mid - sqrt(2.0) * s) < 1e-9,
                      "(D) prec=%ld theta=%.3f: lo=%.9e != sqrt2*sin=%.9e (the "
                      "1/n endpoint constant is wrong)", prec, th, lo_mid,
                      sqrt(2.0) * s);
        AIC_CHECK_MSG(fabs(hi_mid - 4.0 * sqrt(2.0) * s) < 1e-9,
                      "(D) prec=%ld theta=%.3f: hi=%.9e != 4 sqrt2*sin=%.9e (the "
                      "factor-2 endpoint constant is wrong)", prec, th, hi_mid,
                      4.0 * sqrt(2.0) * s);
        if (lo_out) lo_out[k] = 0.5 * (arb_lo_d(lo) + arb_hi_d(lo));
        if (hi_out) hi_out[k] = 0.5 * (arb_lo_d(hi) + arb_hi_d(hi));
        if (prec == PREC)
            printf("    theta=%.2f: lo=%.6e <= true=%.6e <= hi=%.6e\n",
                   th, lo_lb, truev, hi_ub);
        arb_clear(hi); arb_clear(lo);
        aic_ucp_kraus_clear(&phiW);
        acb_mat_clear(W);
    }
    aic_ucp_kraus_clear(&phiI);
    acb_mat_clear(I2);
}

static const double D_THETAS[5] = {0.2, 0.5, 1.0, 1.5, 2.0};

static void test_closed_form(void)
{
    printf("(D) identity vs diag(1,e^{i theta}): 2 sin(theta/2) in [lo,hi]\n");
    run_D_containment(PREC, NULL, NULL, D_THETAS, 5);
}

/* --- (E) precision ladder: (D) at prec=256, endpoints agree 53 vs 256 ----- */
static void test_precision_ladder(void)
{
    printf("(E) prec=256 containment + 53-vs-256 endpoint agreement\n");
    double lo53[5], hi53[5], lo256[5], hi256[5];
    run_D_containment(53, lo53, hi53, D_THETAS, 5);
    run_D_containment(256, lo256, hi256, D_THETAS, 5);
    for (int k = 0; k < 5; k++) {
        AIC_CHECK_MSG(fabs(lo53[k] - lo256[k]) < 1e-12,
                      "(E) theta=%.2f lo drift %.3e", D_THETAS[k],
                      fabs(lo53[k] - lo256[k]));
        AIC_CHECK_MSG(fabs(hi53[k] - hi256[k]) < 1e-12,
                      "(E) theta=%.2f hi drift %.3e", D_THETAS[k],
                      fabs(hi53[k] - hi256[k]));
    }
}

/* --- (F) fail-loud teeth: different n -> the dim guard ABORTS -------------- */

/* child: build phi(d=2), psi(d=3); the dim-guard assert must fire. No args (the
 * shared aic_watchdog child signature is void(void)). */
static void cbd_child_dim_guard(void)
{
    aic_ucp_kraus phi, psi;
    aic_channel_depolarizing(&phi, 2, 0.3, PREC);
    aic_channel_depolarizing(&psi, 3, 0.3, PREC);
    arb_t lo, hi;
    arb_init(lo); arb_init(hi);
    aic_cbnorm_eigfree_distance(lo, hi, &phi, &psi, PREC); /* must abort */
    arb_clear(hi); arb_clear(lo);
    aic_ucp_kraus_clear(&psi);
    aic_ucp_kraus_clear(&phi);
}

static void test_dim_guard_failloud(void)
{
    printf("(F) different-n maps must FAIL LOUD (assert in dim guard, not hang)\n");
    /* Strengthening over the old NULL-needle copy (Rule 5): the dim guard is the
     * live assert(phi->dim_H == psi->dim_H), whose stringized expression appears
     * in the glibc abort message. assert_failloud asserts: not a hang, not a
     * silent exit-0, SIGABRT, and that the abort names this specific guard. */
    aic_watchdog_assert_failloud(cbd_child_dim_guard, 20,
                                 "cbnorm_distance/dim-guard",
                                 "phi->dim_H == psi->dim_H");
    printf("    child SIGABRT on dim mismatch (guard has teeth)\n");
}

int main(void)
{
    test_self_distance_zero();
    test_symmetry();
    test_linearity();
    test_closed_form();
    test_precision_ladder();
    test_dim_guard_failloud();
    aic_test_report("test_cbnorm_distance");
    printf("OK test_cbnorm_distance\n");
    return 0;
}
