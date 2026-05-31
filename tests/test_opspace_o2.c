/* test_opspace_o2.c — the CERTIFIED cb-norm bound tests for opspace O2 (bead
 * aic-pjr): the rect Watrous-SDP UPPER certifier (aic_cbnorm_certify_rect_upper)
 * and the O2 pipeline (aic_opspace_certify_cb_upper) that turn the committed
 * MIN-dual feasible points into rigorous arb bounds on ||v||_cb / ||v^{-1}||_cb,
 * the certified half of th_main_ext (.tex:1538-1540). See
 * docs/research/opspace_o2_design.md §0.5 (the PINNED convention: factor 1, dual
 * tr_sys=2 = MINOR factor, eps*d_min shift) + §4.1 (the HOPM<=SDP bracket).
 *
 * Fixtures: tests/fixtures_opspace_o2.inc.h (3 channels, MOSEK-solved, committed):
 *   block_cond_exp(4,2)     eta=0, N=4 n_B=4 (square)     — complete-isometry oracle
 *   noiseless_subsystem(3,2) eta=0, N=6 n_B=3 (RECTANGULAR) — oracle + direction teeth
 *   mixconj(6,2,0.03)       eta>0, N=6 n_B=2              — the bound + bracket carrier
 * Each carries the Kraus (to rebuild v in C), eps_used, dims, and PER DIRECTION the
 * Convention-A adjoint Choi J and a dual point (Y0,Y1) with eta_ref = the cb-norm.
 *
 * THE CHECKS (Rule 5/6; every check asserts an invariant against a known value,
 * mutation-proven RED in the source-mutation log at the bottom of each block):
 *   T1 eta=0 ORACLE: rect certifier on J(v*)/J((v^-1)*) ENCLOSES 1 from above
 *      (|hi-1| <= 1e-6; restoration eps small) — the complete-isometry oracle
 *      (CLAUDE.md ladder rung 3). MUTATION: drop the eps*d_min shift OR use
 *      2/d_min factor -> hi != 1 -> RED.
 *   T2 eta>0 BOUND + BRACKET: hi_fwd ~ fwd_eta_ref, hi_inv ~ inv_eta_ref
 *      (hi >= eta_ref, gap small); run O1 + certify_cb_upper, assert the bracket
 *      cb_forward <= cb_upper and 1/a_cb <= cbinv_upper HOLD (aic-0at HOPM<=SDP).
 *   T3 DIRECTION TOOTH: the certifier traces tr_sys=2 (MINOR, partial_trace_right).
 *      Compute hi the WRONG way (trace MAJOR, partial_trace_left, shift eps*d_maj)
 *      on the RECTANGULAR fixtures and assert it DIFFERS from eta_ref by >= the
 *      known teeth margin — so the direction is genuinely pinned. Both the FORWARD
 *      (all rect fixtures) and the INVERSE Choi J((v^-1)*) are checked; the inverse
 *      hard assertion is gated on the eta=0 ORACLE (mixconj-inv's wrong-direction
 *      value is coincidentally close — see the mutation log). MUTATION: a
 *      right->left swap in the certifier would produce exactly this wrong hi.
 *   T4 §C12 NON-VACUITY: ||v^-1||_cb = hi_inv ~ 1.535 DIFFERS materially from the
 *      Frobenius proxy 1/sigma_min(M_1) ~ 1.027; assert |cbinv_upper - 1/sigma_min|
 *      > 0.1 — O2 measures the genuine operator/cb norm, not the vacuous Frobenius
 *      sigma_min (FINDINGS §C12).
 *   T5 MIDPOINT-RADIUS FIX TOOTH: artificially inflate a committed forward J's
 *      entrywise radius to ~1e-50 (>> herm_max_eig's Hermiticity tol, << the 1e-4
 *      tightness tol) and assert aic_cbnorm_certify_rect_upper COMPLETES (the
 *      midpoint-collapse fix prevents the false Hermiticity abort) and yields a hi
 *      matching the zero-radius hi to <= 1e-10. MUTATION: removing the collapse
 *      SIGABRTs in herm_max_eig.
 *
 * NOTE (>200 LOC): test files are exempt from the Rule-10 source limit (the sibling
 * test files exceed it too); the SOURCE cores stay under 200.
 */
#include <math.h>
#include <string.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_assoc.h"
#include "aic/aic_cbnorm.h"
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_mat.h"
#include "aic/aic_opspace.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"

#include "fixtures_opspace_o2.inc.h"

#define PREC 256

/* rigorous endpoints of an arb ball as doubles (FLOOR / CEIL). */
static double arb_lo_d(const arb_t x)
{ arf_t t; arf_init(t); arb_get_lbound_arf(t, x, 64);
  double v = arf_get_d(t, ARF_RND_FLOOR); arf_clear(t); return v; }
static double arb_hi_d(const arb_t x)
{ arf_t t; arf_init(t); arb_get_ubound_arf(t, x, 64);
  double v = arf_get_d(t, ARF_RND_CEIL); arf_clear(t); return v; }

/* Load a flat sz x sz [p*sz+q] (re,im) array into an acb_mat. */
static void load_sq(acb_mat_t M, const double *re, const double *im, slong sz)
{ for (slong p = 0; p < sz; p++) for (slong q = 0; q < sz; q++)
    acb_set_d_d(acb_mat_entry(M, p, q), re[p * sz + q], im[p * sz + q]); }

/* Build the UCP self-map phi from a fixture's committed Kraus (Convention A,
 * K-major [a*n*n + i*n + j]) — the same marshalling the shim uses. */
static void phi_from_fixture(aic_ucp_kraus *phi, const aic_opspace_o2_fixture *f)
{
    aic_ucp_kraus_init(phi, f->n, f->n, f->r);
    for (slong a = 0; a < f->r; a++)
        for (slong i = 0; i < f->n; i++)
            for (slong j = 0; j < f->n; j++) {
                slong off = a * f->n * f->n + i * f->n + j;
                acb_set_d_d(acb_mat_entry(phi->K[a], i, j),
                            f->kraus_re[off], f->kraus_im[off]);
            }
}

/* The rebuilt v: B -> A from the committed Kraus + eps_used (the EXACT shim/
 * test_opspace pipeline, so J(v*) reproduces the committed Choi to ~1e-16). */
typedef struct {
    aic_ucp_kraus phi;
    aic_assoc_ecstar ae;
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
} o2_fixture;

static void o2_build(o2_fixture *fx, const aic_opspace_o2_fixture *f)
{
    phi_from_fixture(&fx->phi, f);
    aic_assoc_ecstar_from_phi(&fx->ae, &fx->phi, PREC);
    arb_init(fx->iso);
    aic_cstar_build(&fx->B, &fx->v, fx->iso, &fx->ae.A, f->eps_used, PREC);
    AIC_CHECK_MSG(fx->v.A->dim_A == fx->B.dim_B,
                  "%s: v not bijective dim_A=%ld dim_B=%ld", f->name,
                  (long) fx->v.A->dim_A, (long) fx->B.dim_B);
}

static void o2_clear(o2_fixture *fx)
{
    arb_clear(fx->iso);
    aic_dhom_v_clear(&fx->v);
    aic_dhom_B_clear(&fx->B);
    aic_assoc_ecstar_clear(&fx->ae);
    aic_ucp_kraus_clear(&fx->phi);
}

/* The WRONG-direction UPPER (T3 mutation frame): trace the MAJOR factor
 * (partial_trace_left -> d_min x d_min marginal) with shift eps*d_maj. This is
 * exactly what a partial_trace_right->left swap (+ the matching shift dim) in
 * aic_cbnorm_certify_rect_upper would compute; the teeth assert it != eta_ref. */
static double wrong_direction_hi(const acb_mat_t J, const acb_mat_t Y0,
                                 const acb_mat_t Y1, slong d_maj, slong d_min)
{
    slong D = d_maj * d_min;
    acb_mat_t negJ, negJd, blk, nb;
    acb_mat_init(negJ, D, D); acb_mat_init(negJd, D, D);
    acb_mat_init(blk, 2 * D, 2 * D); acb_mat_init(nb, 2 * D, 2 * D);
    acb_mat_neg(negJ, J);
    acb_mat_conjugate_transpose(negJd, negJ);
    for (slong i = 0; i < D; i++)
        for (slong j = 0; j < D; j++) {
            acb_set(acb_mat_entry(blk, i, j), acb_mat_entry(Y0, i, j));
            acb_set(acb_mat_entry(blk, i, D + j), acb_mat_entry(negJ, i, j));
            acb_set(acb_mat_entry(blk, D + i, j), acb_mat_entry(negJd, i, j));
            acb_set(acb_mat_entry(blk, D + i, D + j), acb_mat_entry(Y1, i, j));
        }
    acb_mat_neg(nb, blk);
    arb_t e; arb_init(e);
    aic_mat_herm_max_eig(e, nb, PREC);
    double eps = arb_hi_d(e); if (eps < 0) eps = 0;
    acb_mat_t T0, T1;
    acb_mat_init(T0, d_min, d_min); acb_mat_init(T1, d_min, d_min);
    aic_mat_partial_trace_left(T0, Y0, d_maj, d_min, PREC);   /* WRONG: major */
    aic_mat_partial_trace_left(T1, Y1, d_maj, d_min, PREC);
    arb_t s0, s1; arb_init(s0); arb_init(s1);
    aic_mat_herm_max_eig(s0, T0, PREC);
    aic_mat_herm_max_eig(s1, T1, PREC);
    double wrong = 0.5 * (arb_hi_d(s0) + arb_hi_d(s1)) + eps * (double) d_maj;
    arb_clear(s1); arb_clear(s0); arb_clear(e);
    acb_mat_clear(T1); acb_mat_clear(T0);
    acb_mat_clear(nb); acb_mat_clear(blk);
    acb_mat_clear(negJd); acb_mat_clear(negJ);
    return wrong;
}

static const aic_opspace_o2_fixture *find_fix(const char *name)
{
    for (int fi = 0; fi < AIC_OPSPACE_O2_NFIX; fi++)
        if (strcmp(aic_opspace_o2_fixtures[fi].name, name) == 0)
            return &aic_opspace_o2_fixtures[fi];
    AIC_FAIL("fixture %s not found", name);
    return NULL;
}

int main(void)
{
    printf("=== test_opspace_o2: certified cb-norm UPPER bound (O2.4/6/7, aic-pjr) ===\n");

    /* ====================== T1 + T3-rect: per-fixture rect certifier ============ */
    for (int fi = 0; fi < AIC_OPSPACE_O2_NFIX; fi++) {
        const aic_opspace_o2_fixture *f = &aic_opspace_o2_fixtures[fi];
        slong fz = f->fwd_dmaj * f->fwd_dmin, iz = f->inv_dmaj * f->inv_dmin;
        acb_mat_t Jf, Y0f, Y1f, Ji, Y0i, Y1i;
        acb_mat_init(Jf, fz, fz); acb_mat_init(Y0f, fz, fz); acb_mat_init(Y1f, fz, fz);
        acb_mat_init(Ji, iz, iz); acb_mat_init(Y0i, iz, iz); acb_mat_init(Y1i, iz, iz);
        load_sq(Jf, f->fwd_J_re, f->fwd_J_im, fz);
        load_sq(Y0f, f->fwd_Y0_re, f->fwd_Y0_im, fz);
        load_sq(Y1f, f->fwd_Y1_re, f->fwd_Y1_im, fz);
        load_sq(Ji, f->inv_J_re, f->inv_J_im, iz);
        load_sq(Y0i, f->inv_Y0_re, f->inv_Y0_im, iz);
        load_sq(Y1i, f->inv_Y1_re, f->inv_Y1_im, iz);

        arb_t hf, hin;
        arb_init(hf); arb_init(hin);
        aic_cbnorm_certify_rect_upper(hf, Jf, Y0f, Y1f, f->fwd_dmaj, f->fwd_dmin, PREC);
        aic_cbnorm_certify_rect_upper(hin, Ji, Y0i, Y1i, f->inv_dmaj, f->inv_dmin, PREC);
        double hf_ub = arb_hi_d(hf), hf_lb = arb_lo_d(hf);
        double hi_ub = arb_hi_d(hin), hi_lb = arb_lo_d(hin);
        printf("  %-22s eta0=%d  fwd hi=[%.10f,%.10f] (ref %.10f)  "
               "inv hi=[%.10f,%.10f] (ref %.10f)\n",
               f->name, f->eta_is_zero, hf_lb, hf_ub, f->fwd_eta_ref,
               hi_lb, hi_ub, f->inv_eta_ref);

        if (f->eta_is_zero) {
            /* T1 ORACLE: hi encloses 1 from above; |hi-1| small. The complete-
             * isometry oracle ||v||_cb = ||v^-1||_cb = 1 EXACT (ladder rung 3). */
            AIC_CHECK_MSG(hf_ub >= 1.0 - 1e-9 && fabs(hf_ub - 1.0) <= 1e-6,
                          "(T1) %s fwd oracle: hi=%.10f not ~1 from above", f->name, hf_ub);
            AIC_CHECK_MSG(hi_ub >= 1.0 - 1e-9 && fabs(hi_ub - 1.0) <= 1e-6,
                          "(T1) %s inv oracle: hi=%.10f not ~1 from above", f->name, hi_ub);
        } else {
            /* T2 BOUND (rect-certifier half): hi is a valid upper bound (>= eta_ref)
             * and tight (gap small). */
            AIC_CHECK_MSG(hf_ub >= f->fwd_eta_ref - 1e-7,
                          "(T2) %s fwd: hi=%.10f below eta_ref=%.10f (not an upper bound)",
                          f->name, hf_ub, f->fwd_eta_ref);
            AIC_CHECK_MSG(hf_ub - f->fwd_eta_ref <= 1e-4,
                          "(T2) %s fwd: hi=%.10f not tight vs eta_ref=%.10f",
                          f->name, hf_ub, f->fwd_eta_ref);
            AIC_CHECK_MSG(hi_ub >= f->inv_eta_ref - 1e-7,
                          "(T2) %s inv: hi=%.10f below eta_ref=%.10f", f->name,
                          hi_ub, f->inv_eta_ref);
            AIC_CHECK_MSG(hi_ub - f->inv_eta_ref <= 1e-4,
                          "(T2) %s inv: hi=%.10f not tight vs eta_ref=%.10f", f->name,
                          hi_ub, f->inv_eta_ref);
        }

        /* T3 DIRECTION TOOTH (forward): the wrong-direction (trace MAJOR) hi must
         * DIFFER from eta_ref by >= the teeth margin. On the SQUARE oracle
         * (block_cond_exp, d_maj==d_min) the two directions coincide (blind, like the
         * symmetric self-map fixtures) — those have NO teeth; the RECTANGULAR fixtures
         * (noiseless_subsystem fwd 6x3, mixconj fwd 6x2) DO. */
        double wrong_fwd = wrong_direction_hi(Jf, Y0f, Y1f, f->fwd_dmaj, f->fwd_dmin);
        double teeth = fabs(wrong_fwd - f->fwd_eta_ref);
        printf("      T3 fwd wrong-dir (trace MAJOR) hi=%.6f  ref=%.6f  teeth=%.4f\n",
               wrong_fwd, f->fwd_eta_ref, teeth);
        if (f->fwd_dmaj != f->fwd_dmin) {
            /* the rectangular fixtures: the swap is OBSERVABLE -> teeth. */
            AIC_CHECK_MSG(teeth >= 0.5,
                          "(T3) %s fwd: wrong-direction hi=%.6f too close to eta_ref=%.6f "
                          "(teeth=%.4f < 0.5) — direction not pinned", f->name,
                          wrong_fwd, f->fwd_eta_ref, teeth);
        }

        /* T3 DIRECTION TOOTH (INVERSE): same wrong-direction (trace MAJOR) check on
         * J((v^-1)*). The inverse Choi is rectangular for the non-square fixtures
         * (noiseless_subsystem inv 3x6, mixconj inv 2x6), so the right->left
         * partial-trace swap is OBSERVABLE — this pins the inverse partial-trace
         * direction directly rather than only via the shared codepath.
         *
         * The CLEAN teeth live on the eta=0 ORACLE (noiseless_subsystem inv): the
         * isometry Choi makes the wrong-direction value the oracle deviation
         * 0.5 (= 1/(d_min/d_maj) for inv d_maj=3,d_min=6) vs ref 1, teeth 0.5.
         * For the eta>0 carrier (mixconj inv 2x6) the trace-MAJOR value (1.5298)
         * is COINCIDENTALLY close to the correct trace-MINOR value (1.5354), teeth
         * only ~0.006 — NOT a usable direction tooth on that fixture (FINDINGS: the
         * wrong-direction probe is non-discriminating on mixconj-inv; the forward
         * mixconj teeth ~3.97 and the oracle inverse teeth carry the pin). So gate
         * the hard assertion on the eta=0 oracle, where the swap is genuinely
         * observable; report the eta>0 value for visibility. */
        double wrong_inv = wrong_direction_hi(Ji, Y0i, Y1i, f->inv_dmaj, f->inv_dmin);
        double teeth_inv = fabs(wrong_inv - f->inv_eta_ref);
        printf("      T3 inv wrong-dir (trace MAJOR) hi=%.6f  ref=%.6f  teeth=%.4f\n",
               wrong_inv, f->inv_eta_ref, teeth_inv);
        if (f->inv_dmaj != f->inv_dmin && f->eta_is_zero) {
            /* oracle inverse: the swap deviates by the clean 0.5 oracle margin. */
            AIC_CHECK_MSG(teeth_inv >= 0.25,
                          "(T3) %s inv: wrong-direction hi=%.6f too close to eta_ref=%.6f "
                          "(teeth=%.4f < 0.25) — inverse direction not pinned", f->name,
                          wrong_inv, f->inv_eta_ref, teeth_inv);
        }

        arb_clear(hf); arb_clear(hin);
        acb_mat_clear(Jf); acb_mat_clear(Y0f); acb_mat_clear(Y1f);
        acb_mat_clear(Ji); acb_mat_clear(Y0i); acb_mat_clear(Y1i);
    }

    /* ======= T2 BRACKET + T4 NON-VACUITY: full pipeline on mixconj (eta>0) ======= */
    {
        const aic_opspace_o2_fixture *f = find_fix("mixconj_6_2_0p03");
        o2_fixture fx;
        o2_build(&fx, f);

        /* sanity: the rebuilt v reproduces the committed J(v*) (so the committed
         * dual (Y0,Y1) is feasible for the assembled J — else the restoration is
         * meaningless). Cross-checked in test_opspace_choi; re-assert here cheaply
         * via the forward Choi diff. */
        slong N = fx.v.n, nB = fx.B.n_B, fz = f->fwd_dmaj * f->fwd_dmin;
        acb_mat_t Jvs, Jc;
        acb_mat_init(Jvs, N * nB, N * nB);
        acb_mat_init(Jc, fz, fz);
        aic_opspace_choi_vstar(Jvs, &fx.v, PREC);
        load_sq(Jc, f->fwd_J_re, f->fwd_J_im, fz);
        double jdiff = 0.0;
        acb_t dd; arb_t mm; acb_init(dd); arb_init(mm);
        for (slong p = 0; p < fz; p++)
            for (slong q = 0; q < fz; q++) {
                acb_sub(dd, acb_mat_entry(Jvs, p, q), acb_mat_entry(Jc, p, q), PREC);
                acb_abs(mm, dd, PREC);
                double dv = arf_get_d(arb_midref(mm), ARF_RND_NEAR);
                if (dv > jdiff) jdiff = dv;
            }
        acb_clear(dd); arb_clear(mm);
        acb_mat_clear(Jc); acb_mat_clear(Jvs);
        AIC_CHECK_MSG(jdiff < 1e-9,
                      "(T2) mixconj: rebuilt J(v*) differs from committed by %.3e — the "
                      "committed dual is not feasible for the assembled J", jdiff);
        printf("  mixconj: rebuilt J(v*) vs committed max diff = %.3e (dual feasible)\n",
               jdiff);

        /* O1 (HOPM lower bounds) then O2 (certified SDP uppers + the bracket guard). */
        aic_opspace_result r;
        aic_opspace_certify(&r, &fx.v, arb_hi_d(fx.iso), PREC);

        acb_mat_t Y0f, Y1f, Y0i, Y1i;
        slong iz = f->inv_dmaj * f->inv_dmin;
        acb_mat_init(Y0f, fz, fz); acb_mat_init(Y1f, fz, fz);
        acb_mat_init(Y0i, iz, iz); acb_mat_init(Y1i, iz, iz);
        load_sq(Y0f, f->fwd_Y0_re, f->fwd_Y0_im, fz);
        load_sq(Y1f, f->fwd_Y1_re, f->fwd_Y1_im, fz);
        load_sq(Y0i, f->inv_Y0_re, f->inv_Y0_im, iz);
        load_sq(Y1i, f->inv_Y1_re, f->inv_Y1_im, iz);

        /* certify_cb_upper assembles J internally + asserts the bracket internally
         * (aborts on violation — that abort IS the bracket test). */
        aic_opspace_certify_cb_upper(&r, &fx.v, Y0f, Y1f, Y0i, Y1i, PREC);

        double cb_fwd = arb_hi_d(r.cb_forward);   /* O1 HOPM lower on ||v||_cb */
        double cb_up = arb_hi_d(r.cb_upper);      /* O2 SDP upper             */
        double a_cb = arb_lo_d(r.a_cb);
        double cbinv_lo = a_cb > 1e-300 ? 1.0 / a_cb : 0.0; /* O1 HOPM lower on ||v^-1|| */
        double cbinv_up = arb_hi_d(r.cbinv_upper); /* O2 SDP upper            */
        printf("  T2 BRACKET fwd: HOPM lower=%.6f <= SDP upper=%.6f (ref %.6f)\n",
               cb_fwd, cb_up, f->fwd_eta_ref);
        printf("  T2 BRACKET inv: HOPM lower=%.6f <= SDP upper=%.6f (ref %.6f)\n",
               cbinv_lo, cbinv_up, f->inv_eta_ref);

        /* re-assert the bracket explicitly (the certify call already aborts on a
         * violation; this makes the assertion VISIBLE in the check count + catches a
         * future relaxation of the internal guard). */
        AIC_CHECK_MSG(cb_fwd <= cb_up + 1e-6,
                      "(T2) fwd bracket: HOPM lower=%.6f > SDP upper=%.6f", cb_fwd, cb_up);
        AIC_CHECK_MSG(cbinv_lo <= cbinv_up + 1e-6,
                      "(T2) inv bracket: HOPM lower=%.6f > SDP upper=%.6f",
                      cbinv_lo, cbinv_up);
        /* and the certified uppers reproduce eta_ref (the SDP optimum). */
        AIC_CHECK_MSG(cb_up >= f->fwd_eta_ref - 1e-7 && cb_up - f->fwd_eta_ref <= 1e-4,
                      "(T2) cb_upper=%.10f != fwd_eta_ref=%.10f", cb_up, f->fwd_eta_ref);
        AIC_CHECK_MSG(cbinv_up >= f->inv_eta_ref - 1e-7 && cbinv_up - f->inv_eta_ref <= 1e-4,
                      "(T2) cbinv_upper=%.10f != inv_eta_ref=%.10f", cbinv_up, f->inv_eta_ref);

        /* T4 §C12 NON-VACUITY: the cb ||v^-1||_cb = cbinv_up DIFFERS materially from
         * the Frobenius proxy 1/sigma_min(M_1). The Frobenius sigma_min is the SAME
         * for every ampliation (I_{n^2} (x) M_1), so a sigma_min "cb canary" is
         * vacuous (aic_opspace.h / FINDINGS §C12); O2 measures the genuine cb norm. */
        arb_t smin;
        arb_init(smin);
        aic_dhom_v_sigma_min(smin, &fx.v, PREC);
        double sig = arb_lo_d(smin);
        double frob_proxy = sig > 1e-300 ? 1.0 / sig : 0.0;
        printf("  T4 NON-VACUITY: ||v^-1||_cb=%.6f vs Frobenius 1/sigma_min(M_1)=%.6f  "
               "gap=%.6f\n", cbinv_up, frob_proxy, fabs(cbinv_up - frob_proxy));
        AIC_CHECK_MSG(fabs(cbinv_up - frob_proxy) > 0.1,
                      "(T4) cb ||v^-1||_cb=%.6f too close to Frobenius proxy=%.6f "
                      "(gap %.6f <= 0.1) — O2 would be measuring the vacuous sigma_min",
                      cbinv_up, frob_proxy, fabs(cbinv_up - frob_proxy));
        arb_clear(smin);

        acb_mat_clear(Y0f); acb_mat_clear(Y1f);
        acb_mat_clear(Y0i); acb_mat_clear(Y1i);
        aic_opspace_result_clear(&r);
        o2_clear(&fx);
    }

    /* ===== T5 MIDPOINT-RADIUS FIX TOOTH: the committed-fixture path loads J at
     * ZERO radius, and the pipeline path's arb-assembled J carries only a ~1.9e-71
     * radius — too tiny to perturb herm_max_eig's relative Hermiticity tol
     * (~2^-(prec-8)*(1+|.|), prec=256). So the midpoint-collapse fix in
     * aic_cbnorm_certify_rect_upper (the arb_get_mid_arb + symmetrize) is never
     * EXERCISED with teeth. Here we ARTIFICIALLY INFLATE the forward J's entrywise
     * radius to ~1e-50 (far above the Hermiticity tol, far below the 1e-4 tightness
     * tol). Because the assembled block's off-diagonal -J and -J^dag are independent
     * balls, their asymmetry ball ~ the inflated radius (1e-50) trips the rigorous
     * Hermiticity assert in herm_max_eig — UNLESS the midpoint fix collapses the
     * radius first. Assert (a) it COMPLETES (no false Hermiticity abort) and (b) the
     * resulting hi matches the zero-radius hi to ~1e-10 (the radius is numerically
     * irrelevant to the bound). MUTATION (source-mutation log): removing the
     * midpoint-collapse loop SIGABRTs here in herm_max_eig. */
    {
        const aic_opspace_o2_fixture *f = find_fix("noiseless_subsystem_3_2");
        slong fz = f->fwd_dmaj * f->fwd_dmin;
        acb_mat_t J0, Jr, Y0f, Y1f;
        acb_mat_init(J0, fz, fz); acb_mat_init(Jr, fz, fz);
        acb_mat_init(Y0f, fz, fz); acb_mat_init(Y1f, fz, fz);
        load_sq(J0, f->fwd_J_re, f->fwd_J_im, fz);   /* zero-radius reference */
        load_sq(Jr, f->fwd_J_re, f->fwd_J_im, fz);   /* to be radius-inflated */
        load_sq(Y0f, f->fwd_Y0_re, f->fwd_Y0_im, fz);
        load_sq(Y1f, f->fwd_Y1_re, f->fwd_Y1_im, fz);

        /* inflate every entry's radius to ~1e-50 (>> Hermiticity tol ~1e-75,
         * << tightness tol 1e-4). acb_add_error_mag widens both re/im parts. */
        mag_t err; mag_init(err); mag_set_d(err, 1e-50);
        for (slong p = 0; p < fz; p++)
            for (slong q = 0; q < fz; q++)
                acb_add_error_mag(acb_mat_entry(Jr, p, q), err);
        mag_clear(err);

        arb_t hi0, hir;
        arb_init(hi0); arb_init(hir);
        aic_cbnorm_certify_rect_upper(hi0, J0, Y0f, Y1f, f->fwd_dmaj, f->fwd_dmin, PREC);
        /* (a) COMPLETES without abort despite the 1e-50 radius (the midpoint fix). */
        aic_cbnorm_certify_rect_upper(hir, Jr, Y0f, Y1f, f->fwd_dmaj, f->fwd_dmin, PREC);
        double v0 = arb_hi_d(hi0), vr = arb_hi_d(hir);
        printf("  T5 MIDPOINT-RADIUS: zero-radius hi=%.12f  inflated(1e-50) hi=%.12f  "
               "diff=%.3e\n", v0, vr, fabs(vr - v0));
        /* the act of returning IS the "no false Hermiticity abort" assertion;
         * record it as a counted check for visibility. */
        AIC_CHECK_MSG(1, "(T5) certify_rect_upper completed on 1e-50-radius J");
        /* (b) the bound is unmoved by the radius. */
        AIC_CHECK_MSG(fabs(vr - v0) <= 1e-10,
                      "(T5) radius-inflated hi=%.12f differs from zero-radius hi=%.12f "
                      "by %.3e > 1e-10 — the radius leaked into the bound", vr, v0,
                      fabs(vr - v0));

        arb_clear(hi0); arb_clear(hir);
        acb_mat_clear(J0); acb_mat_clear(Jr);
        acb_mat_clear(Y0f); acb_mat_clear(Y1f);
    }

    printf("CHECKS test_opspace_o2 n=%ld\n", aic_test_checks);
    aic_test_report("test_opspace_o2");
    printf("OK test_opspace_o2\n");
    return 0;
}

/* ===================== SOURCE-MUTATION LOG (Rule 7: proven RED then reverted) ===
 * Each mutation was applied to the impl, rebuilt, confirmed RED at the cited line,
 * then reverted to GREEN. The mutations and their measured RED:
 *
 *   T1 (eta=0 oracle, normalization). aic_cbnorm_certify_rect_upper: multiply hi by
 *      the SUPERSEDED 2/d_min factor (hi <- (2/d_min)*hi). RED:
 *        "(T1) block_cond_exp_4_2 fwd oracle: hi=0.5000000000 not ~1 from above"
 *      (d_min=4 -> hi=0.5). Pins the design-§0.5 FACTOR-1 normalization.
 *
 *   T3 (partial-trace DIRECTION). aic_cbnorm_certify_rect_upper: swap the dual to
 *      trace the MAJOR factor (partial_trace_left -> d_min marginal, shift eps*d_maj).
 *      RED on the RECTANGULAR oracle:
 *        "(T1) noiseless_subsystem_3_2 fwd oracle: hi=2.0000000000 not ~1 from above"
 *      (the SQUARE block_cond_exp is direction-blind; the rectangular fixture catches
 *      the swap). Pins design-§0.5 tr_sys=2 (the load-bearing HANDOFF:340 bug).
 *
 *   T2 (HOPM<=SDP bracket guard). Inject O1.cb_forward = 2.0 > O2 SDP upper before
 *      aic_opspace_certify_cb_upper. RED (the internal guard aborts):
 *        "aic_opspace_certify_cb_upper: FORWARD bracket VIOLATED: O1 HOPM lower=
 *         2.0e+00 > O2 SDP upper=1.0019683734e+00". Pins the aic-0at cross-check.
 *
 *   T4 (§C12 non-vacuity). Set cbinv_up = 1/sigma_min(M_1) (pretend O2 measures the
 *      vacuous Frobenius sigma_min). RED:
 *        "(T4) cb ||v^-1||_cb=1.026640 too close to Frobenius proxy=1.026640 (gap
 *         0.000000 <= 0.1)". Confirms O2 measures the genuine cb norm (1.535), not
 *      the ampliation-blind Frobenius sigma_min (1.027).
 *
 *   T3-inv (INVERSE partial-trace DIRECTION). Two mutation frames pin it:
 *      (a) certifier swap, src/aic_cbnorm_certify_rect.c: trace the MAJOR factor
 *          (partial_trace_left -> d_min marginal, shift eps*d_maj). RED — the
 *          FORWARD oracle aborts first (shared codepath):
 *            "(T1) noiseless_subsystem_3_2 fwd oracle: hi=2.0000000000 not ~1"
 *          and the inverse oracle would likewise read 0.5 not 1. This proves the
 *          certifier direction is pinned but does NOT isolate the inverse leg.
 *      (b) test no-teeth probe (isolating): set wrong_inv = f->inv_eta_ref in the
 *          inverse-T3 block (a probe with no direction sensitivity). RED at the
 *          inverse-T3 oracle assertion:
 *            "(T3) noiseless_subsystem_3_2 inv: wrong-direction hi=1.000000 too
 *             close to eta_ref=1.000000 (teeth=0.0000 < 0.25) — inverse direction
 *             not pinned".
 *          This proves the inverse-T3 tooth INDEPENDENTLY pins the inverse
 *          partial-trace direction (the oracle wrong-direction trace-MAJOR value
 *          0.5 vs ref 1.0 = teeth 0.5 >= 0.25). FINDING: the eta>0 carrier
 *          (mixconj inv 2x6) has wrong-direction hi=1.5298 COINCIDENTALLY within
 *          ~0.006 of the correct 1.5354, so it carries NO usable inverse tooth —
 *          the hard assertion is gated on the eta=0 oracle, where the swap is
 *          genuinely observable.
 *
 *   T5 (MIDPOINT-RADIUS fix is load-bearing). src/aic_cbnorm_certify_rect.c:
 *      REMOVE the midpoint-collapse loop (arb_get_mid_arb on each blk entry) +
 *      the symmetrize. RED — the 1e-50-radius-inflated forward J's assembled
 *      block has independent -J / -J^dag off-diagonal balls whose asymmetry ~1e-50
 *      trips herm_max_eig's rigorous relative Hermiticity tol:
 *        "aic_mat: input not Hermitian: |H_ij - conj(H_ji)| exceeds the relative
 *         tol tol*(1 + |H_ij| + |H_ji|) (bead aic-2yo)"  (SIGABRT, exit 134, on the
 *         T5 hir inflated path; the zero-radius hi0 path and all committed-fixture
 *         T1-T4 paths still pass first). Restoring the loop -> GREEN with the
 *         inflated and zero-radius bounds IDENTICAL (diff 0.0 <= 1e-10), proving
 *         the collapse drops the radius without moving the certified bound.
 */
