/* test_certify_teeth.c — teeth for two rigor-chain steps in the TIGHT cb-norm
 * certifier (bead aic-m24 increment 3b) that the main test_certify.c corpus could
 * NOT exercise (a hostile review found them untested; CLAUDE.md Rule 6/7, "tests
 * that can't fail"). The SOURCE (src/aic_cbnorm_certify_restore.c) is correct;
 * this file only ADDS assertions so that a WRONG variant of the correct code
 * cannot pass silently. Split from test_certify.c (already >200 LOC, Rule 10);
 * the Makefile auto-discovers tests/test_*.c.
 *
 * The two gaps (both in aic_cbnorm_int_lower / _upper, restore.c):
 *
 *   GAP-ii — psd_defect upper-endpoint choice (restore.c:48-71). psd_defect
 *     returns a RIGOROUS UPPER bound on max(0,-lambda_min(M)) via
 *     arb_get_ubound_arf(herm_max_eig(-M), prec) (restore.c:63). Rigor of BOTH
 *     recipes depends on it being the UPPER (not lower/mid) endpoint. On the
 *     golden corpus MOSEK is so accurate the defect-ball radius is ~1e-36, so the
 *     ub vs lb choice is invisible (swapping ub->lb changes nothing measurable).
 *     Teeth here: drive the REAL static psd_defect (through aic_cbnorm_int_upper)
 *     on a matrix with a genuinely negative eigenvalue evaluated at LOW prec=8,
 *     where the herm_max_eig global enclosure has width ~3.7e-2 >> 1e-6, so ub and
 *     lb are 3.9e-2 apart; assert the returned defect == the upper endpoint.
 *
 *   GAP-iv — the E-symmetrization (restore.c:79-88): E=P+Q-I; P'=P-E/2, Q'=Q-E/2.
 *     It makes P'+Q'=I_{N} EXACTLY, which is what the restored primal point
 *     (tX, tP'+(1-t)I/2, tQ'+(1-t)I/2) needs to be FEASIBLE (its marginal sum is
 *     t(P'+Q')+(1-t)I = I iff P'+Q'=I) -> the weak-duality LOWER bound lo<=eta.
 *     On the corpus ||E||~1e-7 (below the 1e-7 slack), so dropping it is invisible.
 *     Teeth here: inject a LARGE PSD equality violation P -> P + h*I (h=0.1, so
 *     ||E||=0.1) and scale X -> c*X (c=2) so the raw objective ~ c*eta exceeds eta;
 *     with the symmetrization the corrected lo stays <= eta (rigor preserved); the
 *     mutation that drops it (block built from raw P,Q) gives lo > eta (rigor
 *     broken). Each test below carries its mutation-proof RED in a source comment.
 *
 * Inert-by-construction NOTES the reviewer asked be recorded (so a future
 * maintainer does NOT assume these paths are covered here):
 *
 *   NIT-A (J-conjugation path). Both recipes form J^dag (acb_mat_conjugate_
 *     transpose, restore.c:107 / via -J at :138). On THIS corpus every Choi-diff
 *     J = Choi(Phi^2-Phi) is HERMITIAN (the fixtures' Phi are Hermitian-preserving,
 *     so Phi^2-Phi is too, and its Choi matrix is Hermitian), hence J^dag = J and a
 *     transpose-vs-conjugate-transpose bug would be SILENT here. The complex_qubit
 *     fixture has complex entries but a Hermitian J, so it does NOT exercise the
 *     conjugation either. A genuinely non-Hermitian-J teeth fixture is future work
 *     (not built; recorded so coverage is not assumed).
 *
 *   NIT-B (distinct-dual-marginal path). aic_cbnorm_int_upper sums the two MAX-eigs
 *     of the partial-trace marginals Tr_R(Y0), Tr_R(Y1) (restore.c:150-151) with
 *     weight 1/2 each. At the SDP optimum the two marginals are EQUAL (the dual is
 *     symmetric under Y0<->Y1 at optimality on this corpus), so a bug that, say,
 *     doubled one weight and zeroed the other would land on the same value and be
 *     SILENT. The (s0+s1)/2 averaging is therefore inert-by-construction here; a
 *     fixture with distinct optimal marginals would be needed to give it teeth
 *     (not built; recorded so coverage is not assumed).
 */
#include <math.h>
#include <string.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic_cbnorm.h"
#include "aic_mat.h"
#include "aic_ucp.h"
#include "aic_test.h"

/* The teeth read the internal restoration defects (delta/eps) and drive the
 * recipes directly; the public aic_cbnorm_certify hides them. The internal header
 * is under src/ (not a public surface), reached by relative path here only. */
#include "../src/aic_cbnorm_internal.h"

#include "fixtures_certify.inc.h"

#define PREC 128

/* lower/upper endpoints of an arb ball as doubles (rigorous: FLOOR / CEIL). */
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

/* Load a flat N x N [p*N+q] (re,im) array into an acb_mat. */
static void load_NN(acb_mat_t M, const double *re, const double *im, slong N)
{
    for (slong p = 0; p < N; p++)
        for (slong q = 0; q < N; q++)
            acb_set_d_d(acb_mat_entry(M, p, q), re[p * N + q], im[p * N + q]);
}

/* Build J = Choi(Phi^2 - Phi) from the fixture's Kraus ops. */
static void build_J(acb_mat_t J, const aic_certify_fixture *f, slong prec)
{
    aic_ucp_kraus phi, phi2;
    aic_ucp_kraus_init(&phi, f->n, f->n, f->r);
    for (slong a = 0; a < f->r; a++)
        for (slong i = 0; i < f->n; i++)
            for (slong j = 0; j < f->n; j++) {
                slong off = a * f->n * f->n + i * f->n + j;
                acb_set_d_d(acb_mat_entry(phi.K[a], i, j),
                            f->kraus_re[off], f->kraus_im[off]);
            }
    aic_ucp_compose(&phi2, &phi, &phi, prec);
    aic_ucp_choi_diff(J, &phi2, &phi, prec);
    aic_ucp_kraus_clear(&phi2);
    aic_ucp_kraus_clear(&phi);
}

static const aic_certify_fixture *find_fixture(const char *name)
{
    for (int fi = 0; fi < AIC_CERTIFY_NFIX; fi++)
        if (strcmp(aic_certify_fixtures[fi].name, name) == 0)
            return &aic_certify_fixtures[fi];
    AIC_FAIL("fixture %s not found", name);
    return NULL;
}

/* Assemble [[A,B],[C,D]] (each m x m) into a 2m x 2m matrix — same layout as the
 * static block2 in restore.c, replicated here so the test can rebuild block_D
 * INDEPENDENTLY and read the upper endpoint psd_defect should have returned. */
static void block2(acb_mat_t out, const acb_mat_t A, const acb_mat_t B,
                   const acb_mat_t C, const acb_mat_t D, slong m)
{
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++) {
            acb_set(acb_mat_entry(out, i, j), acb_mat_entry(A, i, j));
            acb_set(acb_mat_entry(out, i, m + j), acb_mat_entry(B, i, j));
            acb_set(acb_mat_entry(out, m + i, j), acb_mat_entry(C, i, j));
            acb_set(acb_mat_entry(out, m + i, m + j), acb_mat_entry(D, i, j));
        }
}

/* ===================================================================== */
/* GAP-iv: the E-symmetrization (P'=P-E/2, Q'=Q-E/2; restore.c:79-88).    */
/* ===================================================================== */
static void teeth_symmetrization(void)
{
    /* phi_t_0p3 (n=2, eta_ref=0.315). Inject a LARGE PSD equality violation
     * P -> P + h*I (h=0.1 => E=P+Q-I has ||E||_op=0.1, far above the 1e-7 slack)
     * AND scale X -> c*X (c=2). The raw objective (2/n)Re tr(J^dag (cX)) = c*eta
     * exceeds eta; only the t-scaling can pull lo back, and the t-scaling is only
     * RIGOROUS when the restored point is feasible, i.e. when P'+Q'=I — which is
     * exactly what the symmetrization provides.
     *
     * Why the X-scaling is necessary (not gold-plating): with X unperturbed the
     * objective is pinned at eta and t<=1 caps lo<=eta in BOTH the symmetrized and
     * the dropped-symmetrization paths, so the mutation could not go RED. c>1 lifts
     * the objective above eta so the (in)feasibility of the marginal sum becomes
     * load-bearing. (The corpus G.LOWER test in test_certify.c uses the same
     * c*X device.) Adding +h*I to P alone makes the no-sym top-left block corner
     * MORE positive -> SMALLER defect -> LARGER t in the dropped path, so the
     * dropped-symmetrization lo overshoots eta. */
    const aic_certify_fixture *f = find_fixture("phi_t_0p3");
    slong n = f->n, N = n * n;
    const double h = 0.1, c = 2.0;

    acb_mat_t J, X, P, Q, Hh;
    acb_mat_init(J, N, N); acb_mat_init(X, N, N);
    acb_mat_init(P, N, N); acb_mat_init(Q, N, N); acb_mat_init(Hh, N, N);
    build_J(J, f, PREC);
    load_NN(X, f->X_re, f->X_im, N);
    load_NN(P, f->P_re, f->P_im, N);
    load_NN(Q, f->Q_re, f->Q_im, N);

    /* H = I (Hermitian PSD); P -> P + h*I. */
    acb_mat_one(Hh);
    arb_t hb; arb_init(hb); arb_set_d(hb, h);
    for (slong i = 0; i < N; i++)
        for (slong j = 0; j < N; j++) {
            acb_t tmp; acb_init(tmp);
            acb_mul_arb(tmp, acb_mat_entry(Hh, i, j), hb, PREC);
            acb_add(acb_mat_entry(P, i, j), acb_mat_entry(P, i, j), tmp, PREC);
            acb_clear(tmp);
        }
    arb_clear(hb);
    /* X -> c*X (c integer). */
    acb_mat_scalar_mul_si(X, X, (slong) c, PREC);

    /* confirm ||E||_op = ||P+Q-I||_op ~ h = 0.1 (the injected violation is real). */
    acb_mat_t E, I_N;
    acb_mat_init(E, N, N); acb_mat_init(I_N, N, N);
    acb_mat_one(I_N);
    acb_mat_add(E, P, Q, PREC);
    acb_mat_sub(E, E, I_N, PREC);
    arb_t enorm; arb_init(enorm);
    aic_mat_opnorm(enorm, E, PREC);
    double e_op = arb_hi_d(enorm);
    arb_clear(enorm);
    acb_mat_clear(I_N); acb_mat_clear(E);

    /* the REAL certifier path (symmetrization present). */
    arb_t lo, delta;
    arb_init(lo); arb_init(delta);
    aic_cbnorm_int_lower(lo, delta, J, X, P, Q, n, PREC);
    double lo_corr = arb_hi_d(lo);
    double delta_lb = arb_lo_d(delta);

    printf("  GAP-iv (E-symmetrization, restore.c:79-88) on phi_t_0p3:\n");
    printf("        h=%.2f c=%.1f  ||E||_op=%.4f  delta(lb)=%.4e  "
           "corrected lo=%.6f  eta_ref=%.6f\n",
           h, c, e_op, delta_lb, lo_corr, f->eta_ref);

    const double slack = 1e-7;
    /* the injected violation is genuinely large (precondition for teeth). */
    AIC_CHECK_MSG(e_op > 0.05,
                  "(GAP-iv) ||E||_op=%.4e too small: violation not injected", e_op);
    /* the restoration fired (delta>0, so t<1) — otherwise nothing to symmetrize. */
    AIC_CHECK_MSG(delta_lb > 1e-3,
                  "(GAP-iv) delta=%.4e not > 0: restoration inert", delta_lb);
    /* (i) RIGOR PRESERVED: with the symmetrization, lo is a valid lower bound. */
    AIC_CHECK_MSG(lo_corr <= f->eta_ref + slack,
                  "(GAP-iv) rigor BROKEN: corrected lo=%.10e > eta_ref=%.10e",
                  lo_corr, f->eta_ref);

    /* (ii) MUTATION-PROOF (documents the teeth). Dropping the symmetrization means
     * building the block straight from the perturbed (P, Q) with P+Q != I. The
     * dropped-symmetrization lo (computed here as the mutation would) overshoots
     * eta_ref, so the source mutation turns assertion (i) above RED.
     * Source mutation reproduced: in aic_cbnorm_int_lower, skip P'=P-E/2,Q'=Q-E/2
     * and block2(blk, P, X, Xd, Q, N) directly. */
    acb_mat_t Xd, blk, negblk;
    acb_mat_init(Xd, N, N); acb_mat_init(blk, 2 * N, 2 * N);
    acb_mat_init(negblk, 2 * N, 2 * N);
    acb_mat_conjugate_transpose(Xd, X);
    block2(blk, P, X, Xd, Q, N);              /* NO symmetrization: raw P, Q */
    acb_mat_neg(negblk, blk);
    arb_t mx; arb_init(mx);
    aic_mat_herm_max_eig(mx, negblk, PREC);
    arb_t dmut; arb_init(dmut);
    arf_t ub; arf_init(ub);
    arb_get_ubound_arf(ub, mx, PREC);
    arf_set(arb_midref(dmut), ub);
    mag_zero(arb_radref(dmut));
    if (arf_sgn(ub) < 0) arb_zero(dmut);
    arf_clear(ub);
    /* t_mut = 1/(1+2*delta_mut); lo_mut = t_mut * (2/n) Re tr(J^dag X). */
    arb_t tmut; arb_init(tmut);
    arb_mul_2exp_si(tmut, dmut, 1);
    arb_add_ui(tmut, tmut, 1, PREC);
    arb_inv(tmut, tmut, PREC);
    acb_mat_t Jd, Mm;
    acb_mat_init(Jd, N, N); acb_mat_init(Mm, N, N);
    acb_mat_conjugate_transpose(Jd, J);
    acb_mat_mul(Mm, Jd, X, PREC);
    acb_t tr; acb_init(tr);
    acb_mat_trace(tr, Mm, PREC);
    arb_t re; arb_init(re);
    acb_get_real(re, tr);
    arb_mul_2exp_si(re, re, 1);
    arb_div_ui(re, re, (ulong) n, PREC);
    arb_t lo_mut; arb_init(lo_mut);
    arb_mul(lo_mut, tmut, re, PREC);
    double lo_mut_d = arb_hi_d(lo_mut);
    printf("        MUTATION (drop symmetrization): lo_mut=%.6f  "
           "eta_ref=%.6f  -> %s\n", lo_mut_d, f->eta_ref,
           lo_mut_d > f->eta_ref + slack ? "RED (rigor broken)" : "GREEN");
    AIC_CHECK_MSG(lo_mut_d > f->eta_ref + slack,
                  "(GAP-iv) mutation-guard vacuous: lo_mut=%.6e !> eta_ref=%.6e "
                  "(dropping the symmetrization would NOT be caught)",
                  lo_mut_d, f->eta_ref);

    arb_clear(lo_mut); arb_clear(re); acb_clear(tr);
    acb_mat_clear(Mm); acb_mat_clear(Jd);
    arb_clear(tmut); arb_clear(dmut); arb_clear(mx);
    acb_mat_clear(negblk); acb_mat_clear(blk); acb_mat_clear(Xd);
    arb_clear(delta); arb_clear(lo);
    acb_mat_clear(Hh); acb_mat_clear(Q); acb_mat_clear(P);
    acb_mat_clear(X); acb_mat_clear(J);
}

/* ===================================================================== */
/* GAP-ii: psd_defect upper-endpoint choice (restore.c:48-71).           */
/* ===================================================================== */
static void teeth_psd_defect_ubound(void)
{
    /* Drive the REAL static psd_defect through aic_cbnorm_int_upper. Engineer
     * J=0, Y0,Y1 so block_D=[[Y0,-J],[-J^dag,Y1]] = diag(Y0,Y1) with a single
     * eigenvalue -0.3 (the rest +0.7): lambda_min(block_D)=-0.3, defect 0.3.
     * Evaluate at LOW prec=8: the herm_max_eig global enclosure of -block_D is
     * WIDE (radius ~1.8e-2, ub-lb ~3.7e-2 at 64-bit read; the prec-8 ubound and
     * lbound endpoints are 3.9e-2 apart) so the ub vs lb choice is OBSERVABLE.
     * psd_defect must return the UPPER endpoint; the eps out-param of int_upper IS
     * that defect, so we assert eps == arb_get_ubound_arf(herm_max_eig(-blkD), 8). */
    slong n = 2, N = n * n;       /* block_D is 2N x 2N = 8 x 8 */
    const slong lowprec = 8;

    acb_mat_t J, Y0, Y1;
    acb_mat_init(J, N, N); acb_mat_init(Y0, N, N); acb_mat_init(Y1, N, N);
    acb_mat_zero(J);              /* J=0 -> block_D = diag(Y0, Y1) */
    acb_mat_zero(Y0); acb_mat_zero(Y1);
    acb_set_d(acb_mat_entry(Y0, 0, 0), -0.3);          /* the lambda_min = -0.3 */
    for (slong i = 1; i < N; i++) acb_set_d(acb_mat_entry(Y0, i, i), 0.7);
    for (slong i = 0; i < N; i++) acb_set_d(acb_mat_entry(Y1, i, i), 0.7);

    arb_t hi, eps;
    arb_init(hi); arb_init(eps);
    aic_cbnorm_int_upper(hi, eps, J, Y0, Y1, n, lowprec);

    /* Independent block_D and herm_max_eig(-block_D) at the SAME low prec — same
     * matrix, same deterministic enclosure, so the prec-8 endpoints are bit-exact
     * vs what psd_defect saw. */
    acb_mat_t negJ, negJd, blk, negblk;
    acb_mat_init(negJ, N, N); acb_mat_init(negJd, N, N);
    acb_mat_init(blk, 2 * N, 2 * N); acb_mat_init(negblk, 2 * N, 2 * N);
    acb_mat_neg(negJ, J);
    acb_mat_conjugate_transpose(negJd, negJ);
    block2(blk, Y0, negJ, negJd, Y1, N);
    acb_mat_neg(negblk, blk);
    arb_t mx; arb_init(mx);
    aic_mat_herm_max_eig(mx, negblk, lowprec);
    arf_t ub, lb;
    arf_init(ub); arf_init(lb);
    arb_get_ubound_arf(ub, mx, lowprec);   /* what psd_defect extracts */
    arb_get_lbound_arf(lb, mx, lowprec);   /* the WRONG (lb mutation) endpoint */
    double ub_d = arf_get_d(ub, ARF_RND_NEAR);
    double lb_d = arf_get_d(lb, ARF_RND_NEAR);
    double width_64 = arb_hi_d(mx) - arb_lo_d(mx);

    printf("  GAP-ii (psd_defect ub-choice, restore.c:48-71) at prec=%ld:\n",
           (long) lowprec);
    printf("        herm_max_eig(-block_D) enclosure: ub=%.8e lb=%.8e  "
           "ub-lb(prec)=%.6e  width(64-bit)=%.6e\n",
           ub_d, lb_d, ub_d - lb_d, width_64);
    printf("        eps (psd_defect via int_upper) = %.8e\n", arb_hi_d(eps));

    /* The enclosure must be WIDE so ub != lb (else no teeth). */
    AIC_CHECK_MSG(ub_d - lb_d > 1e-6,
                  "(GAP-ii) enclosure ub-lb=%.6e !> 1e-6: ub/lb indistinguishable",
                  ub_d - lb_d);
    /* psd_defect returns the UPPER endpoint exactly (eps is a zero-radius point at
     * the ubound; bit-identical to our independent extraction). */
    arf_t epsmid; arf_init(epsmid);
    arf_set(epsmid, arb_midref(eps));
    AIC_CHECK_MSG(arf_equal(epsmid, ub),
                  "(GAP-ii) psd_defect != upper endpoint: eps=%.10e ub=%.10e",
                  arf_get_d(epsmid, ARF_RND_NEAR), ub_d);
    /* and it is DISTINGUISHABLE from the lower endpoint by > 1e-6 (so the ub->lb
     * source mutation moves eps by ~3.9e-2 and turns the equality above RED). */
    AIC_CHECK_MSG(arb_hi_d(eps) > lb_d + 1e-6,
                  "(GAP-ii) eps=%.10e not > lb+1e-6=%.10e: ub vs lb has no teeth",
                  arb_hi_d(eps), lb_d + 1e-6);
    /* Source mutation reproduced: in psd_defect, arb_get_ubound_arf ->
     * arb_get_lbound_arf (restore.c:63). Then eps = lb_d = 0.2793 and the
     * arf_equal(epsmid, ub) assertion above goes RED. */
    printf("        MUTATION (ub->lb in psd_defect): eps would be %.8e != ub "
           "%.8e -> RED\n", lb_d, ub_d);

    arf_clear(epsmid);
    arf_clear(lb); arf_clear(ub);
    arb_clear(mx);
    acb_mat_clear(negblk); acb_mat_clear(blk);
    acb_mat_clear(negJd); acb_mat_clear(negJ);
    arb_clear(eps); arb_clear(hi);
    acb_mat_clear(Y1); acb_mat_clear(Y0); acb_mat_clear(J);
}

int main(void)
{
    teeth_psd_defect_ubound();   /* GAP-ii */
    teeth_symmetrization();      /* GAP-iv */

    printf("CHECKS test_certify_teeth n=%ld\n", aic_test_checks);
    aic_test_report("test_certify_teeth");
    printf("OK test_certify_teeth\n");
    return 0;
}
