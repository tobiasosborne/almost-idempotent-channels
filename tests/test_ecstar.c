/* test_ecstar.c — cross-checks for the ecstar module (bead aic-knm): the eps-C*
 * algebra data model + certified-arb axiom-defect estimators. Realizes the eta=0
 * oracle rung for the Choi-Effros star (approximate_algebras.tex:341-344,
 * :2187-2237). Reuses the exact-idempotent channels in tests/test_idemp.h.
 *
 * Cross-check ladder coverage:
 *   A. eta=0 ORACLE (rung 3): for each exact-idempotent channel, build A=Img Phi
 *      and assert ALL FIVE defects ~0. involution/unit at 1e-12 (gauge-clean,
 *      always-zero invariants); assoc/submult/cstar at 1e-9 (they route Phi
 *      through the LAPACK-extracted Delta-basis; floor measured, documented).
 *   B. TEETH (MANDATORY; Rule 7): deliberately non-C* instances make each
 *      estimator demonstrably nonzero and O(t). THE TEETH MAP:
 *        - Phi_t = (1-t) Phi + t Dep (UCP, Hermicity-preserving) evaluated as the
 *          star while keeping the exact A=Img Phi basis: B_k is NO LONGER a
 *          Phi_t fixed point, so submult, cstar, assoc, and unit grow O(t);
 *          involution stays ~0 (Hermicity-preserved). 10x ratio t=1e-3 -> 1e-2.
 *        - a basis operator pushed OFF A: the unit projection residual / term (b)
 *          ||Phi(B_k)-B_k|| becomes nonzero.
 *        - a NON-Hermicity-preserving Phi: the involution defect becomes nonzero
 *          (the only thing that moves it).
 *   C. GAUGE invariance: a random unitary gauge on the d basis vectors leaves all
 *      five defects unchanged (catches a basis-dependent estimator).
 *   D. arb@53 vs arb@256: the oracle defects agree across precision (rung 2).
 */
#include <assert.h>
#include <complex.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_ecstar.h"
#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"
#include "test_idemp.h"

/* ---- helpers ------------------------------------------------------------ */

/* All five defects into out[5] = {assoc, submult, cstar, involution, unit}. */
static void all_defects(arb_ptr out, const aic_ecstar *A, slong prec)
{
    aic_ecstar_defect_assoc(out + 0, A, prec);
    aic_ecstar_defect_submult(out + 1, A, prec);
    aic_ecstar_defect_cstar(out + 2, A, prec);
    aic_ecstar_defect_involution(out + 3, A, prec);
    aic_ecstar_defect_unit(out + 4, A, prec);
}

static double defect_dbl(const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

/* out = s * A (scalar double multiply; FLINT has no acb_mat_scalar_mul_d). */
static void mat_scale_d(acb_mat_t out, const acb_mat_t A, double s, slong prec)
{
    acb_t c;
    acb_init(c);
    acb_set_d(c, s);
    acb_mat_scalar_mul_acb(out, A, c, prec);
    acb_clear(c);
}

/* Build a trace-replace ("Dep") channel on C^n: Dep(X) = (Tr X / n) 1_n. */
static void make_dep(aic_ucp_kraus *dep, slong n)
{
    make_trace_replace(dep, n);
}

/* Phi_t = (1-t) Phi + t Dep as a Kraus set: {sqrt(1-t) K_a} U {sqrt(t) L_b}.
 * Both Phi and Dep are UCP + Hermicity-preserving, so Phi_t is too; the unital
 * relation sum K^dag K = (1-t) 1 + t 1 = 1 holds exactly. `out` is init'd here. */
static void make_phi_mix(aic_ucp_kraus *out, const aic_ucp_kraus *phi,
                         const aic_ucp_kraus *dep, double t, slong prec)
{
    slong n = phi->dim_H;
    assert(phi->dim_K == n && dep->dim_K == n && dep->dim_H == n);
    aic_ucp_kraus_init(out, n, n, phi->r + dep->r);
    double s0 = sqrt(1.0 - t), s1 = sqrt(t);
    for (slong a = 0; a < phi->r; a++)
        mat_scale_d(out->K[a], phi->K[a], s0, prec);
    for (slong b = 0; b < dep->r; b++)
        mat_scale_d(out->K[phi->r + b], dep->K[b], s1, prec);
}

/* Phi_scaled = (1+delta) Phi: scale every Kraus op by sqrt(1+delta). This is
 * NON-contractive (||Phi_scaled||_cb = 1+delta > 1, NOT UCP), so it can violate
 * submultiplicativity upward: ||Phi_scaled(XY)|| = (1+delta) ||Phi(XY)||. Used
 * ONLY to give the submult estimator teeth (a UCP star is submultiplicative by
 * construction; see teeth_phi_mix). `out` is init'd here. */
static void make_phi_scale(aic_ucp_kraus *out, const aic_ucp_kraus *phi,
                           double delta, slong prec)
{
    slong n = phi->dim_H;
    aic_ucp_kraus_init(out, n, n, phi->r);
    double s = sqrt(1.0 + delta);
    for (slong a = 0; a < phi->r; a++)
        mat_scale_d(out->K[a], phi->K[a], s, prec);
}

/* ---- A. eta=0 oracle ---------------------------------------------------- */

/* Decompose phi, build A = Img Phi, assert all five defects below the gates;
 * print the measured ceilings (Rule 12: concrete numbers). Returns nothing. */
static void oracle_channel(const char *name, aic_ucp_kraus *phi, slong prec)
{
    arb_t tol_clean, tol_route;
    arb_init(tol_clean);
    arb_init(tol_route);
    arb_set_d(tol_clean, 1e-12); /* involution, unit: gauge-clean invariants  */
    arb_set_d(tol_route, 1e-9);  /* assoc, submult, cstar: route Phi via Delta */

    aic_idemp_decomp d;
    aic_idemp_decompose(&d, phi, prec);

    aic_ecstar A;
    aic_ecstar_init(&A, d.n, d.dim_A, phi);
    aic_ecstar_from_idemp(&A, &d, phi, prec);

    arb_struct def[5];
    for (int i = 0; i < 5; i++) arb_init(&def[i]);
    all_defects(def, &A, prec);

    AIC_CHECK_MSG(arb_le(&def[0], tol_route), "%s: assoc defect too big", name);
    AIC_CHECK_MSG(arb_le(&def[1], tol_route), "%s: submult defect too big", name);
    AIC_CHECK_MSG(arb_le(&def[2], tol_route), "%s: cstar defect too big", name);
    AIC_CHECK_MSG(arb_le(&def[3], tol_clean), "%s: involution defect too big", name);
    AIC_CHECK_MSG(arb_le(&def[4], tol_clean), "%s: unit defect too big", name);

    printf("ORACLE %-26s assoc=%.2e submult=%.2e cstar=%.2e invol=%.2e unit=%.2e\n",
           name, defect_dbl(&def[0]), defect_dbl(&def[1]), defect_dbl(&def[2]),
           defect_dbl(&def[3]), defect_dbl(&def[4]));

    for (int i = 0; i < 5; i++) arb_clear(&def[i]);
    aic_ecstar_clear(&A);
    aic_idemp_clear(&d);
    arb_clear(tol_route);
    arb_clear(tol_clean);
}

static void test_oracle(void)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 5, 2); oracle_channel("block_cond_exp(5,2)", &phi, PREC);
    aic_ucp_kraus_clear(&phi);
    make_block_cond_exp(&phi, 3, 1); oracle_channel("block_cond_exp(3,1)", &phi, PREC);
    aic_ucp_kraus_clear(&phi);
    make_trace_replace(&phi, 3);     oracle_channel("trace_replace(3)", &phi, PREC);
    aic_ucp_kraus_clear(&phi);
    make_noiseless_subsystem(&phi, 2, 2);
    oracle_channel("noiseless_subsystem(2,2)", &phi, PREC);
    aic_ucp_kraus_clear(&phi);
    make_identity(&phi, 3);          oracle_channel("identity(3)", &phi, PREC);
    aic_ucp_kraus_clear(&phi);
    make_dephasing(&phi, 4);         oracle_channel("dephasing(4)", &phi, PREC);
    aic_ucp_kraus_clear(&phi);
    make_compress_idemp(&phi, 4, 2); oracle_channel("compress_idemp(4,2)", &phi, PREC);
    aic_ucp_kraus_clear(&phi);

    /* ASYMMETRIC oracle: conjugate an exactly-idempotent channel by a FIXED
     * complex unitary V (Kraus V^dag K_a V). Phi' stays exactly idempotent and
     * unital, but its image algebra V^dag (Img Phi) V is NOT transpose-closed —
     * giving the from_idemp ROW-MAJOR reshape a real tooth (the other 7 channels
     * all had transpose-closed images, so a TRANSPOSED reshape spanned the same
     * subspace and was undetectable). The base is compress_idemp (NON-self-dual:
     * its Kraus |e_0><e_{m+b}| are not Hermitian, so Phi != Phi*), so the
     * conjugate ALSO catches a Phi<->Phi* dual swap independently — the
     * Phi-direction defense no longer rests on the un-conjugated compress_idemp
     * alone. (block_cond_exp is self-dual: Phi=sum P_a X P_a = Phi*, so its
     * conjugate stays self-dual and could NOT catch the swap.) compress_idemp
     * also has a PROPER carrier M (dim_M=2 < 4). */
    {
        aic_ucp_kraus base, conj;
        make_compress_idemp(&base, 4, 2);   /* dim_A = 4, dim_M = 2 (M (-> H) */
        make_conjugated(&conj, &base, PREC);
        /* Verify the conjugate is STILL exactly idempotent before the oracle
         * runs (the conjugation must preserve Phi^2=Phi; V unitary to ~1e-15). */
        {
            aic_ucp_kraus sq;
            aic_ucp_compose(&sq, &conj, &conj, PREC);
            acb_mat_t Cd;
            acb_mat_init(Cd, 16, 16);
            aic_ucp_choi_diff(Cd, &sq, &conj, PREC); /* Choi(Phi'^2) - Choi(Phi') */
            arb_t idef;
            arb_init(idef);
            aic_mat_opnorm(idef, Cd, PREC);
            printf("ASYM conjugated compress_idemp(4,2): idempotence Choi-defect=%.2e\n",
                   defect_dbl(idef));
            AIC_CHECK_MSG(defect_dbl(idef) < 1e-9,
                          "asym: conjugated channel not exactly idempotent (%.3e)",
                          defect_dbl(idef));
            arb_clear(idef);
            acb_mat_clear(Cd);
            aic_ucp_kraus_clear(&sq);
        }
        oracle_channel("conj_compress_idemp(4,2)", &conj, PREC);
        aic_ucp_kraus_clear(&conj);
        aic_ucp_kraus_clear(&base);
    }
}

/* ---- B. teeth ----------------------------------------------------------- */

/* The Phi_t-mix teeth: build A's basis from the EXACT Phi, but evaluate the star
 * with Phi_t. assoc/submult/cstar/unit grow O(t); involution stays ~0. */
static void teeth_phi_mix(void)
{
    aic_ucp_kraus phi, dep;
    make_block_cond_exp(&phi, 4, 2); /* dim_A = 2^2 + 2^2 = 8 */
    make_dep(&dep, 4);

    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);

    const double ts[3] = {0.0, 1e-3, 1e-2};
    double dassoc[3], dsub[3], dcst[3], dinv[3], dunit[3];
    for (int it = 0; it < 3; it++) {
        aic_ucp_kraus phit;
        make_phi_mix(&phit, &phi, &dep, ts[it], PREC);
        aic_ecstar A;
        aic_ecstar_init(&A, d.n, d.dim_A, &phit);
        /* basis from exact Phi (orthonormal); star uses phit */
        aic_ecstar_from_idemp(&A, &d, &phit, PREC);
        arb_struct def[5];
        for (int i = 0; i < 5; i++) arb_init(&def[i]);
        all_defects(def, &A, PREC);
        dassoc[it] = defect_dbl(&def[0]);
        dsub[it]   = defect_dbl(&def[1]);
        dcst[it]   = defect_dbl(&def[2]);
        dinv[it]   = defect_dbl(&def[3]);
        dunit[it]  = defect_dbl(&def[4]);
        for (int i = 0; i < 5; i++) arb_clear(&def[i]);
        aic_ecstar_clear(&A);
        aic_ucp_kraus_clear(&phit);
    }

    printf("TEETH phi_mix t=0     assoc=%.2e submult=%.2e cstar=%.2e invol=%.2e unit=%.2e\n",
           dassoc[0], dsub[0], dcst[0], dinv[0], dunit[0]);
    printf("TEETH phi_mix t=1e-3  assoc=%.2e submult=%.2e cstar=%.2e invol=%.2e unit=%.2e\n",
           dassoc[1], dsub[1], dcst[1], dinv[1], dunit[1]);
    printf("TEETH phi_mix t=1e-2  assoc=%.2e submult=%.2e cstar=%.2e invol=%.2e unit=%.2e\n",
           dassoc[2], dsub[2], dcst[2], dinv[2], dunit[2]);

    /* t=0: everything ~0 */
    AIC_CHECK(dassoc[0] < 1e-9 && dsub[0] < 1e-9 && dcst[0] < 1e-9 &&
              dunit[0] < 1e-9);
    /* assoc, cstar, unit grow O(t): nonzero at t=1e-3, ~10x at t=1e-2 ([4x,25x]). */
    AIC_CHECK_MSG(dassoc[1] > 1e-6, "phi_mix: assoc did not move at t=1e-3");
    AIC_CHECK_MSG(dassoc[2] / dassoc[1] > 4.0 && dassoc[2] / dassoc[1] < 25.0,
                  "phi_mix: assoc not ~10x (got %.3g)", dassoc[2] / dassoc[1]);
    AIC_CHECK_MSG(dcst[1] > 1e-6, "phi_mix: cstar did not move at t=1e-3");
    AIC_CHECK_MSG(dcst[2] / dcst[1] > 4.0 && dcst[2] / dcst[1] < 25.0,
                  "phi_mix: cstar not ~10x (got %.3g)", dcst[2] / dcst[1]);
    AIC_CHECK_MSG(dunit[1] > 1e-6, "phi_mix: unit did not move at t=1e-3");
    AIC_CHECK_MSG(dunit[2] / dunit[1] > 4.0 && dunit[2] / dunit[1] < 25.0,
                  "phi_mix: unit not ~10x (got %.3g)", dunit[2] / dunit[1]);
    /* involution stays ~0 (Phi_t is Hermicity-preserving). */
    AIC_CHECK_MSG(dinv[2] < 1e-9, "phi_mix: involution moved but Phi_t is HP");
    /* submult stays ~0 under a UCP mix: ANY UCP-defined star is submultiplicative
     * (||Phi(XY)|| <= ||Phi||_cb ||XY|| <= ||X|| ||Y||, ||Phi||_cb=1), so
     * ax_prodnorm holds with eps<=0 for every Phi_t. Mixing UCP maps CANNOT
     * violate it upward. The submult estimator's teeth need a NON-contractive
     * map (teeth_submult below). FINDING: this is a structural property of the
     * Choi-Effros star, not a test gap. */
    AIC_CHECK_MSG(dsub[2] < 1e-9,
                  "phi_mix: submult positive under a UCP mix (should be 0: UCP "
                  "star is submultiplicative)");

    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&dep);
    aic_ucp_kraus_clear(&phi);
}

/* The submult teeth: a NON-contractive map Phi' = (1+delta) Phi makes
 * ||B_j*'B_k|| = (1+delta)||B_j*B_k|| exceed ||B_j|| ||B_k|| where the C* algebra
 * was tight (e.g. a block projector B with ||B*B|| = ||B||^2), so the submult
 * slack becomes positive O(delta). The basis is the EXACT A=Img Phi basis; only
 * the map that defines the star is scaled. */
static void teeth_submult(void)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 4, 2); /* fixed points include block projectors */
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);

    const double ds[3] = {0.0, 1e-3, 1e-2};
    double dsub[3];
    for (int it = 0; it < 3; it++) {
        aic_ucp_kraus phis;
        make_phi_scale(&phis, &phi, ds[it], PREC);
        aic_ecstar A;
        aic_ecstar_init(&A, d.n, d.dim_A, &phis);
        aic_ecstar_from_idemp(&A, &d, &phis, PREC);
        arb_t s;
        arb_init(s);
        aic_ecstar_defect_submult(s, &A, PREC);
        dsub[it] = defect_dbl(s);
        arb_clear(s);
        aic_ecstar_clear(&A);
        aic_ucp_kraus_clear(&phis);
    }
    printf("TEETH submult delta=0 sub=%.2e  delta=1e-3 sub=%.2e  delta=1e-2 sub=%.2e\n",
           dsub[0], dsub[1], dsub[2]);
    AIC_CHECK_MSG(dsub[0] < 1e-9, "submult: nonzero at delta=0");
    AIC_CHECK_MSG(dsub[1] > 1e-6, "submult: did not move at delta=1e-3");
    AIC_CHECK_MSG(dsub[2] / dsub[1] > 4.0 && dsub[2] / dsub[1] < 25.0,
                  "submult: not ~10x (got %.3g)", dsub[2] / dsub[1]);
    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&phi);
}

/* Basis-off-A teeth: push one stored basis operator off Img Phi by t*D where D
 * is a matrix unit (generically not in A); the unit defect term ||Phi(B_k)-B_k||
 * then becomes O(t). Exercises aic_ecstar_defect_unit's fixed-point sweep. */
static void teeth_basis_off_A(void)
{
    aic_ucp_kraus phi;
    make_dephasing(&phi, 4); /* A = diagonal matrices, dim_A = 4 */
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);

    const double ts[3] = {0.0, 1e-3, 1e-2};
    double dunit[3];
    for (int it = 0; it < 3; it++) {
        aic_ecstar A;
        aic_ecstar_init(&A, d.n, d.dim_A, &phi);
        aic_ecstar_from_idemp(&A, &d, &phi, PREC);
        /* add t * E_{0,1} (an off-diagonal unit, NOT in the diagonal algebra) to
         * B[0]; Phi(off-diagonal) = 0 != off-diagonal, so the fixed-point term
         * ||Phi(B_0) - B_0|| picks up O(t). B[0] is diagonal, so entry (0,1) is 0
         * and setting it to t adds exactly t*E_{0,1}. */
        acb_set_d(acb_mat_entry(A.B[0], 0, 1), ts[it]);
        arb_t u;
        arb_init(u);
        aic_ecstar_defect_unit(u, &A, PREC);
        dunit[it] = defect_dbl(u);
        arb_clear(u);
        aic_ecstar_clear(&A);
    }
    printf("TEETH off_A   t=0 unit=%.2e  t=1e-3 unit=%.2e  t=1e-2 unit=%.2e\n",
           dunit[0], dunit[1], dunit[2]);
    AIC_CHECK_MSG(dunit[0] < 1e-9, "off_A: unit nonzero at t=0");
    AIC_CHECK_MSG(dunit[1] > 1e-5, "off_A: unit did not move at t=1e-3");
    AIC_CHECK_MSG(dunit[2] / dunit[1] > 4.0 && dunit[2] / dunit[1] < 25.0,
                  "off_A: unit not ~10x (got %.3g)", dunit[2] / dunit[1]);

    /* Projection-residual teeth for the unit estimator's term (a) "is 1_n in A?":
     * a 1-d algebra spanned by E_00 (Phi=id) does NOT contain the unit 1_2, so
     * aic_ecstar_proj_residual(1_2) is O(1). With B=E_00, ||1_2 - <E_00,1_2>E_00||
     * = ||1_2 - E_00||_op = ||E_11||_op = 1. Proves term (a) is not vacuous. */
    {
        aic_ucp_kraus idp;
        make_identity(&idp, 2);
        aic_ecstar Ar;
        aic_ecstar_init(&Ar, 2, 1, &idp);
        acb_set_si(acb_mat_entry(Ar.B[0], 0, 0), 1); /* B = E_00 (Frobenius-unit) */
        acb_mat_t I2;
        acb_mat_init(I2, 2, 2);
        acb_mat_one(I2);
        arb_t res, u;
        arb_init(res);
        arb_init(u);
        aic_ecstar_proj_residual(res, &Ar, I2, PREC);
        aic_ecstar_defect_unit(u, &Ar, PREC);
        printf("TEETH proj_res 1_2 off span{E_00}: residual=%.3f unit_defect=%.3f\n",
               defect_dbl(res), defect_dbl(u));
        AIC_CHECK_MSG(defect_dbl(res) > 0.9,
                      "proj_residual: 1_2 outside span{E_00} read as in-A (got %.3e)",
                      defect_dbl(res));
        AIC_CHECK_MSG(defect_dbl(u) > 0.9,
                      "unit term (a): vacuous (1_2 not in A but defect %.3e)",
                      defect_dbl(u));
        arb_clear(u);
        arb_clear(res);
        acb_mat_clear(I2);
        aic_ecstar_clear(&Ar);
        aic_ucp_kraus_clear(&idp);
    }

    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&phi);
}

/* A synthetic NON-Hermicity-preserving (non-HP) one-knob map driven THROUGH the
 * production core aic_ecstar_involution_core (so a return-0 / wrong-matrix
 * mutation of the core is detectable):
 *   T_t(X) = X (I + t G),  G a FIXED non-Hermitian matrix.
 * Then T_t(X)^dag = (I + t G^dag) X^dag while T_t(X^dag) = X^dag (I + t G), so
 * the core's per-pair term apply(Y)^dag - apply(Y^dag) = t (G^dag Y^dag - Y^dag G)
 * is O(t) and generically nonzero (T_t is NOT HP for t != 0; T_0 = id is HP).
 * ctx is the t-scaled (I + t G) matrix N; T_t(X) = X N. */
typedef struct { const acb_mat_struct *N; } nonhp_ctx;

static void apply_nonhp(acb_mat_t out, const acb_mat_t X, void *ctx, slong prec)
{
    const nonhp_ctx *c = (const nonhp_ctx *) ctx;
    acb_mat_mul(out, X, c->N, prec);  /* T_t(X) = X (I + t G) */
}

/* Involution teeth (Rule 7, mutation-proof, NON-VACUOUS). The production
 * estimator aic_ecstar_defect_involution is a STRUCTURAL always-zero: with
 *   (B_j*B_k)^dag - B_k^dag*B_j^dag = Phi(B_jB_k)^dag - Phi((B_jB_k)^dag),
 * it is identically 0 for any Hermicity-preserving (HP) Phi, and aic_ucp_kraus
 * CANNOT encode a non-HP map (Phi(X)=sum K^dag X K is HP for ANY Kraus set). A
 * correct impl and `return 0` are observationally identical on all aic_ecstar
 * inputs. So the residual computation is exercised through the SHARED core
 * aic_ecstar_involution_core, which the production estimator wraps unchanged.
 * Here we drive that SAME core with a synthetic non-HP map T_t(X)=X(I+tG):
 *   t=0   : T_0=id is HP  -> core residual ~0;
 *   t=1e-3: O(t)          -> residual moves;
 *   t=1e-2: ~10x t=1e-3   -> residual is genuinely O(non-HP-ness).
 * Mutating the core to `arb_zero(out); return;` (or transposing/dropping an
 * adjoint) makes the t>0 assertions RED, so the production code path is
 * protected. The HP structural-zero is also kept (the production call). */
static void teeth_non_hp(void)
{
    /* A FIXED non-normal, non-Hermitian basis {B_0, B_1} (n=2), so the core's
     * adjoint bookkeeping is non-trivial. (Frobenius-orthonormality is not
     * required by the core; it only reads B as a matrix list.) */
    acb_mat_t B[2];
    for (int b = 0; b < 2; b++) acb_mat_init(B[b], 2, 2);
    acb_set_d(acb_mat_entry(B[0], 0, 0), 1.0 / sqrt(2.0));
    acb_set_d_d(acb_mat_entry(B[0], 0, 1), 0.0, 1.0 / sqrt(2.0)); /* + i/sqrt2 */
    acb_set_d(acb_mat_entry(B[1], 1, 0), 0.6);
    acb_set_d_d(acb_mat_entry(B[1], 0, 1), 0.3, 0.4);

    /* G fixed, non-Hermitian: G[0,1]=1 only (G != G^dag). N_t = I + t G. */
    const double ts[3] = {0.0, 1e-3, 1e-2};
    double res[3];
    for (int it = 0; it < 3; it++) {
        acb_mat_t N;
        acb_mat_init(N, 2, 2);
        acb_mat_one(N);
        acb_set_d(acb_mat_entry(N, 0, 1), ts[it]); /* + t * G, G = E_01 */
        nonhp_ctx c = { N };
        arb_t r;
        arb_init(r);
        aic_ecstar_involution_core(r, (const acb_mat_t *) B, 2, 2, apply_nonhp,
                                   &c, PREC);
        res[it] = defect_dbl(r);
        arb_clear(r);
        acb_mat_clear(N);
    }

    /* Production structural-zero kept: HP Phi (=id) through the WHOLE public
     * estimator reads ~0 (this also confirms the production wrapper path). */
    aic_ucp_kraus phi;
    make_identity(&phi, 2);
    aic_ecstar A;
    aic_ecstar_init(&A, 2, 2, &phi);
    acb_mat_set(A.B[0], B[0]);
    acb_mat_set(A.B[1], B[1]);
    arb_t inv;
    arb_init(inv);
    aic_ecstar_defect_involution(inv, &A, PREC);
    double inv_hp = defect_dbl(inv);

    printf("TEETH non_hp  invol(HP prod)=%.2e  core t=0=%.2e t=1e-3=%.2e t=1e-2=%.2e\n",
           inv_hp, res[0], res[1], res[2]);

    /* HP structural-zero. */
    AIC_CHECK_MSG(inv_hp < 1e-12,
                  "non_hp: production involution must be 0 for HP Phi (got %.3e)",
                  inv_hp);
    AIC_CHECK_MSG(res[0] < 1e-12,
                  "non_hp: core residual nonzero for HP map T_0=id (got %.3e)",
                  res[0]);
    /* The CORE is exercised and non-vacuous: O(t) growth, ~10x t=1e-3 -> 1e-2.
     * A return-0 mutation of the core makes BOTH of these fire. */
    AIC_CHECK_MSG(res[1] > 1e-6,
                  "non_hp: core residual did not move under a non-HP map at "
                  "t=1e-3 (VACUOUS core; got %.3e)", res[1]);
    AIC_CHECK_MSG(res[2] / res[1] > 4.0 && res[2] / res[1] < 25.0,
                  "non_hp: core residual not ~10x t=1e-3 -> 1e-2 (got %.3g)",
                  res[2] / res[1]);

    arb_clear(inv);
    aic_ecstar_clear(&A);
    aic_ucp_kraus_clear(&phi);
    for (int b = 0; b < 2; b++) acb_mat_clear(B[b]);
}

/* ---- C. gauge invariance ------------------------------------------------ */

/* Apply a random d x d unitary gauge to the basis: B'_k = sum_l U[l,k] B_l, and
 * assert all five defects are unchanged. A 2-d real rotation suffices to expose a
 * basis-dependent estimator. */
static void test_gauge(void)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 4, 2); /* dim_A = 8 */
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);

    aic_ecstar A0, A1;
    aic_ecstar_init(&A0, d.n, d.dim_A, &phi);
    aic_ecstar_from_idemp(&A0, &d, &phi, PREC);
    aic_ecstar_init(&A1, d.n, d.dim_A, &phi);
    aic_ecstar_from_idemp(&A1, &d, &phi, PREC);

    slong n = d.n, dd = d.dim_A;
    /* Gauge: a chain of plane rotations on consecutive coordinate pairs (a real
     * orthogonal matrix, hence unitary). B'_k = sum_l U[l,k] B_l. */
    double theta = 0.7;
    double c = cos(theta), s = sin(theta);
    /* rotate pairs (0,1),(2,3),...: U is block-diagonal 2x2 rotations. Build B' in
     * a scratch array, then copy back. */
    acb_mat_struct *Bp = flint_malloc((size_t) dd * sizeof(acb_mat_struct));
    for (slong k = 0; k < dd; k++) {
        acb_mat_init(&Bp[k], n, n);
        acb_mat_set(&Bp[k], A1.B[k]);
    }
    for (slong p = 0; p + 1 < dd; p += 2) {
        /* B'_p = c B_p - s B_{p+1};  B'_{p+1} = s B_p + c B_{p+1} (columns of U) */
        acb_mat_t t0, t1;
        acb_mat_init(t0, n, n);
        acb_mat_init(t1, n, n);
        mat_scale_d(t0, A1.B[p], c, PREC);
        mat_scale_d(t1, A1.B[p + 1], -s, PREC);
        acb_mat_add(&Bp[p], t0, t1, PREC);
        mat_scale_d(t0, A1.B[p], s, PREC);
        mat_scale_d(t1, A1.B[p + 1], c, PREC);
        acb_mat_add(&Bp[p + 1], t0, t1, PREC);
        acb_mat_clear(t1);
        acb_mat_clear(t0);
    }
    for (slong k = 0; k < dd; k++) acb_mat_set(A1.B[k], &Bp[k]);
    for (slong k = 0; k < dd; k++) acb_mat_clear(&Bp[k]);
    flint_free(Bp);

    arb_struct d0[5], d1[5];
    for (int i = 0; i < 5; i++) { arb_init(&d0[i]); arb_init(&d1[i]); }
    all_defects(d0, &A0, PREC);
    all_defects(d1, &A1, PREC);
    arb_t tol;
    arb_init(tol);
    arb_set_d(tol, 1e-9);
    double maxdiff = 0.0;
    for (int i = 0; i < 5; i++) {
        double diff = fabs(defect_dbl(&d0[i]) - defect_dbl(&d1[i]));
        if (diff > maxdiff) maxdiff = diff;
        arb_t db;
        arb_init(db);
        arb_sub(db, &d0[i], &d1[i], PREC);
        arb_abs(db, db);
        AIC_CHECK_MSG(arb_le(db, tol), "gauge: defect %d changed (diff %.3e)", i,
                      diff);
        arb_clear(db);
    }
    printf("GAUGE max defect diff under unitary gauge = %.3e\n", maxdiff);
    arb_clear(tol);
    for (int i = 0; i < 5; i++) { arb_clear(&d0[i]); arb_clear(&d1[i]); }
    aic_ecstar_clear(&A1);
    aic_ecstar_clear(&A0);
    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&phi);
}

/* ---- D. arb@53 vs arb@256 ----------------------------------------------- */

static void test_prec(void)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 5, 2);
    aic_idemp_decomp d53, d256;
    aic_idemp_decompose(&d53, &phi, 53);
    aic_idemp_decompose(&d256, &phi, 256);

    aic_ecstar A53, A256;
    aic_ecstar_init(&A53, d53.n, d53.dim_A, &phi);
    aic_ecstar_from_idemp(&A53, &d53, &phi, 53);
    aic_ecstar_init(&A256, d256.n, d256.dim_A, &phi);
    aic_ecstar_from_idemp(&A256, &d256, &phi, 256);

    arb_struct f53[5], f256[5];
    for (int i = 0; i < 5; i++) { arb_init(&f53[i]); arb_init(&f256[i]); }
    all_defects(f53, &A53, 53);
    all_defects(f256, &A256, 256);
    arb_t tol;
    arb_init(tol);
    arb_set_d(tol, 1e-9);
    double maxdiff = 0.0;
    for (int i = 0; i < 5; i++) {
        double diff = fabs(defect_dbl(&f53[i]) - defect_dbl(&f256[i]));
        if (diff > maxdiff) maxdiff = diff;
        arb_t db;
        arb_init(db);
        arb_sub(db, &f53[i], &f256[i], 256);
        arb_abs(db, db);
        AIC_CHECK_MSG(arb_le(db, tol), "prec: defect %d differs 53 vs 256 (%.3e)",
                      i, diff);
        arb_clear(db);
    }
    printf("PREC arb@53 vs arb@256 max defect diff = %.3e\n", maxdiff);
    arb_clear(tol);
    for (int i = 0; i < 5; i++) { arb_clear(&f53[i]); arb_clear(&f256[i]); }
    aic_ecstar_clear(&A256);
    aic_ecstar_clear(&A53);
    aic_idemp_clear(&d256);
    aic_idemp_clear(&d53);
    aic_ucp_kraus_clear(&phi);
}

/* ---- Cycle 2: FAITHFUL worst-case HOPM search --------------------------- *
 * The basis-sweep estimators above are exact-zero detectors / cheap lower
 * bounds; the HOPM routines (aic_ecstar_defect_*_hopm) are the second audition
 * candidate (Law 4): a scale-invariant alternating-maximization search over the
 * OPERATOR-norm unit ball of A, returning a RIGOROUS LOWER BOUND (the lo end of a
 * certified arb ball). The tests below establish, in order:
 *   E1. eta=0 oracle: all three HOPM lower bounds ~0 on the exact-idempotent
 *       channels (the search seeks a defect that is genuinely 0).
 *   E2. audition payoff: HOPM_assoc >= basis-sweep_assoc on a perturbed instance
 *       (HOPM is the faithful/tighter lower bound for the trilinear associator).
 *   E3. THE UNIVERSALITY CANARY: at a FIXED perturbation t across growing dim,
 *       the HOPM defect/t ratio does NOT grow (dimension-independent), while the
 *       basis-sweep ratio drifts (the d^{3/2} trap). The headline test.
 *   E4. mutation-proof the search SEARCHES (Rule 7): the weakest search
 *       (n_restarts=1, n_iters=0) finds a STRICTLY SMALLER defect than the full
 *       search on a perturbed instance (a "return 0 / never iterate" stub is
 *       caught), per defect.
 *   E5. certified-ball soundness: the returned lo <= the ratio re-evaluated at
 *       the explicit witness, and the witness is in A (proj_residual ~0).
 *
 * Search counts are kept modest to bound suite time (convergence near eta=0 is
 * fast; the global max is found in ~10-25 iters, see the bead notes). */

/* assoc defect/t and (optionally) basis-sweep/t at a FIXED t for a
 * block_cond_exp(d,m) family mixed with `dep`. Returns both ratios via out-params.
 * The basis-sweep assoc is O(dim_A^3) arb star products (=91125 at dim_A=45), so
 * `do_bs` lets the dimension-sweep skip it at the large-dim point (where the
 * basis sweep is not asserted -- only HOPM/t flatness is); the contrast sweep
 * computes it at the small dims where it is cheap. bs_over_t is set to -1 when
 * skipped. */
static void canary_point(double *hopm_over_t, double *bs_over_t, slong d, slong m,
                         const aic_ucp_kraus *dep, double t, int rs, int its,
                         int do_bs)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, d, m);
    aic_idemp_decomp dd;
    aic_idemp_decompose(&dd, &phi, PREC);
    aic_ucp_kraus phit;
    make_phi_mix(&phit, &phi, dep, t, PREC);
    aic_ecstar A;
    aic_ecstar_init(&A, dd.n, dd.dim_A, &phit);
    aic_ecstar_from_idemp(&A, &dd, &phit, PREC);
    arb_t bs, la;
    arb_init(bs);
    arb_init(la);
    *bs_over_t = -1.0;
    if (do_bs) { aic_ecstar_defect_assoc(bs, &A, PREC); *bs_over_t = defect_dbl(bs) / t; }
    aic_ecstar_defect_assoc_hopm(la, &A, rs, its, PREC);
    *hopm_over_t = defect_dbl(la) / t;
    if (do_bs)
        printf("  CANARY bce(%ld,%ld) dim_A=%2ld  basis_sweep/t=%.4e  HOPM/t=%.4e\n",
               (long) d, (long) m, (long) dd.dim_A, *bs_over_t, *hopm_over_t);
    else
        printf("  CANARY bce(%ld,%ld) dim_A=%2ld  basis_sweep/t=(skipped, O(d^3))  HOPM/t=%.4e\n",
               (long) d, (long) m, (long) dd.dim_A, *hopm_over_t);
    arb_clear(la);
    arb_clear(bs);
    aic_ecstar_clear(&A);
    aic_idemp_clear(&dd);
    aic_ucp_kraus_clear(&phit);
    aic_ucp_kraus_clear(&phi);
}

/* E1: eta=0 oracle — all three HOPM lower bounds ~0 on the 7+1 channels. */
static void hopm_oracle_channel(const char *name, aic_ucp_kraus *phi)
{
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, phi, PREC);
    aic_ecstar A;
    aic_ecstar_init(&A, d.n, d.dim_A, phi);
    aic_ecstar_from_idemp(&A, &d, phi, PREC);
    arb_t la, ls, lc;
    arb_init(la);
    arb_init(ls);
    arb_init(lc);
    aic_ecstar_defect_assoc_hopm(la, &A, 4, 8, PREC);
    aic_ecstar_defect_submult_hopm(ls, &A, 4, 8, PREC);
    aic_ecstar_defect_cstar_hopm(lc, &A, 4, 8, PREC);
    double a = defect_dbl(la), s = defect_dbl(ls), c = defect_dbl(lc);
    printf("HOPM-ORACLE %-24s assoc_lo=%.2e submult_lo=%.2e cstar_lo=%.2e\n",
           name, a, s, c);
    /* the genuine C* algebra has zero defect; the search must read ~0 (a tiny
     * negative is fine — submult/cstar lower bounds clamp at >=0 conceptually but
     * the certified ratio can dip slightly below by rounding). */
    AIC_CHECK_MSG(a < 1e-7, "%s: HOPM assoc_lo not ~0 (%.3e)", name, a);
    AIC_CHECK_MSG(s < 1e-7, "%s: HOPM submult_lo not ~0 (%.3e)", name, s);
    AIC_CHECK_MSG(c < 1e-7, "%s: HOPM cstar_lo not ~0 (%.3e)", name, c);
    arb_clear(lc);
    arb_clear(ls);
    arb_clear(la);
    aic_ecstar_clear(&A);
    aic_idemp_clear(&d);
}

static void test_hopm_oracle(void)
{
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 5, 2); hopm_oracle_channel("block_cond_exp(5,2)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_trace_replace(&phi, 3);     hopm_oracle_channel("trace_replace(3)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_noiseless_subsystem(&phi, 2, 2);
    hopm_oracle_channel("noiseless_subsystem(2,2)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_identity(&phi, 3);          hopm_oracle_channel("identity(3)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_dephasing(&phi, 4);         hopm_oracle_channel("dephasing(4)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_compress_idemp(&phi, 4, 2); hopm_oracle_channel("compress_idemp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);
    /* asymmetric conjugated channel (the transpose-closure tooth from Cycle 1). */
    {
        aic_ucp_kraus base, conj;
        make_compress_idemp(&base, 4, 2);
        make_conjugated(&conj, &base, PREC);
        hopm_oracle_channel("conj_compress_idemp(4,2)", &conj);
        aic_ucp_kraus_clear(&conj);
        aic_ucp_kraus_clear(&base);
    }
}

/* E2: HOPM_assoc >= basis-sweep_assoc on a phi_mix-perturbed instance, across t.
 * The associator basis sweep is a genuine LOWER bound on the trilinear sup, so the
 * faithful HOPM must meet or beat it. (For submult/cstar the basis sweep measures
 * a DIFFERENT per-basis-element quantity, not normalized to op-norm-1, so the >=
 * relation need not hold there; this audition payoff is asserted for assoc.) */
static void test_hopm_vs_basis_sweep(void)
{
    aic_ucp_kraus phi, dep;
    make_block_cond_exp(&phi, 4, 2);
    make_dep(&dep, 4);
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);
    const double ts[2] = {1e-3, 1e-2};
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phit;
        make_phi_mix(&phit, &phi, &dep, ts[it], PREC);
        aic_ecstar A;
        aic_ecstar_init(&A, d.n, d.dim_A, &phit);
        aic_ecstar_from_idemp(&A, &d, &phit, PREC);
        arb_t bs, la;
        arb_init(bs);
        arb_init(la);
        aic_ecstar_defect_assoc(bs, &A, PREC);
        aic_ecstar_defect_assoc_hopm(la, &A, 12, 15, PREC);
        double b = defect_dbl(bs), h = defect_dbl(la);
        printf("AUDITION t=%.0e  assoc basis_sweep=%.4e  HOPM_lo=%.4e  (HOPM/BS=%.2f)\n",
               ts[it], b, h, b > 1e-30 ? h / b : 0.0);
        AIC_CHECK_MSG(h >= b - 1e-9,
                      "audition: HOPM assoc %.3e < basis-sweep %.3e (t=%.0e)", h, b,
                      ts[it]);
        AIC_CHECK_MSG(h > b + 1e-9 || ts[it] < 1e-30,
                      "audition: HOPM not strictly tighter than basis-sweep");
        arb_clear(la);
        arb_clear(bs);
        aic_ecstar_clear(&A);
        aic_ucp_kraus_clear(&phit);
    }
    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&dep);
    aic_ucp_kraus_clear(&phi);
}

/* E3: THE UNIVERSALITY CANARY. Fixed t=1e-2 across block_cond_exp at growing dim.
 * The HOPM defect/t ratio must NOT grow with dim (stays within a constant factor)
 * -- that is the dimension-independence the operator-norm-sphere search buys (no
 * d^{3/2} Frobenius inflation). Two sweeps:
 *   (A) DIMENSION SWEEP d=4,6,9 (dim_A=8,18,45) with a DEPHASING perturbation
 *       (d Kraus -> cheap, the 3-point sweep stays fast); asserts HOPM/t flat.
 *   (B) CONTRAST d=4,6 with a TRACE_REPLACE perturbation, where the basis-sweep
 *       ratio visibly DRIFTS DOWN with dim (0.25 -> 0.17) while HOPM stays high
 *       (~1.4-1.5) -- exposing the basis sweep's dimension-dependence the faithful
 *       search avoids. (trace_replace at d=9 has ~90 Kraus and is too slow for the
 *       suite, so the 3-point sweep uses the cheap dephasing perturbation.) */
static void test_hopm_canary(void)
{
    double t = 1e-2;
    double h4, b4, h6, b6, h9, b9;
    printf("UNIVERSALITY CANARY (A) dimension sweep, dephasing-mix, assoc/t, t=%.0e:\n", t);
    {
        aic_ucp_kraus dep; make_dephasing(&dep, 4);
        canary_point(&h4, &b4, 4, 2, &dep, t, 10, 12, 1);
        aic_ucp_kraus_clear(&dep);
    }
    {
        aic_ucp_kraus dep; make_dephasing(&dep, 6);
        canary_point(&h6, &b6, 6, 3, &dep, t, 6, 8, 1);
        aic_ucp_kraus_clear(&dep);
    }
    {
        aic_ucp_kraus dep; make_dephasing(&dep, 9);
        canary_point(&h9, &b9, 9, 3, &dep, t, 3, 6, 0);
        aic_ucp_kraus_clear(&dep);
    }
    (void) b4; (void) b6; (void) b9; /* printed inside canary_point; A asserts only HOPM */
    /* HOPM/t stays within a constant factor: bounded above (does NOT grow with
     * dim) AND bounded below away from 0 (the search genuinely finds the defect at
     * every dim). */
    double hmin = h4, hmax = h4;
    double hs[3] = {h4, h6, h9};
    for (int i = 0; i < 3; i++) { if (hs[i] < hmin) hmin = hs[i]; if (hs[i] > hmax) hmax = hs[i]; }
    AIC_CHECK_MSG(hmin > 0.3,
                  "canary: HOPM/t fell near 0 at some dim (min %.3e) -- search "
                  "failed to find the defect", hmin);
    AIC_CHECK_MSG(hmax / hmin < 4.0,
                  "canary: HOPM/t GREW with dim (max/min=%.2f) -- the d^{3/2} trap "
                  "(universality FAIL)", hmax / hmin);
    printf("  -> HOPM/t in [%.3f, %.3f], spread factor %.2f (flat = dim-independent)\n",
           hmin, hmax, hmax / hmin);

    /* (B) the trace_replace contrast: basis sweep DRIFTS DOWN with dim. */
    printf("UNIVERSALITY CANARY (B) contrast, trace_replace-mix, assoc/t, t=%.0e:\n", t);
    double hc4, bc4, hc6, bc6;
    {
        aic_ucp_kraus dep; make_trace_replace(&dep, 4);
        canary_point(&hc4, &bc4, 4, 2, &dep, t, 10, 12, 1);
        aic_ucp_kraus_clear(&dep);
    }
    {   /* second contrast point at bce(5,2): dim_A=13 (basis-sweep O(d^3) is
         * 2.6s here vs 12s at dim_A=18 -- the contrast drift is already visible
         * from dim_A=8 to 13). */
        aic_ucp_kraus dep; make_trace_replace(&dep, 5);
        canary_point(&hc6, &bc6, 5, 2, &dep, t, 8, 10, 1);
        aic_ucp_kraus_clear(&dep);
    }
    printf("  -> basis_sweep/t DRIFTS %.4f -> %.4f (down %.0f%%); "
           "HOPM/t stays %.4f -> %.4f (flat)\n",
           bc4, bc6, 100.0 * (1.0 - bc6 / bc4), hc4, hc6);
    /* HOPM stays high & flat where the basis sweep drops: the faithful search does
     * NOT inherit the basis sweep's dimension drift. */
    AIC_CHECK_MSG(hc4 > 1.0 && hc6 > 1.0,
                  "canary(B): HOPM/t dropped below 1.0 (%.3f, %.3f)", hc4, hc6);
    AIC_CHECK_MSG(bc6 < bc4,
                  "canary(B): basis_sweep/t did not drift down with dim (%.3f -> %.3f)",
                  bc4, bc6);
}

/* E4: mutation-proof the search SEARCHES (Rule 7). On a perturbed instance the
 * weakest search (1 restart, 0 iters: a single un-iterated init) must find a
 * STRICTLY SMALLER defect than the full search -- so a "return 0 / never iterate"
 * stub is caught for each defect. Also confirms the full search recovers a value
 * bounded below by the hand-checked basis-sweep witness for assoc. */
static void test_hopm_searches(void)
{
    aic_ucp_kraus phi, dep;
    make_block_cond_exp(&phi, 4, 2);
    make_dep(&dep, 4);
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);
    aic_ucp_kraus phit;
    make_phi_mix(&phit, &phi, &dep, 1e-2, PREC);
    aic_ecstar A;
    aic_ecstar_init(&A, d.n, d.dim_A, &phit);
    aic_ecstar_from_idemp(&A, &d, &phit, PREC);

    arb_t a0, a1, s0, s1, c0, c1, bs;
    arb_init(a0); arb_init(a1); arb_init(s0); arb_init(s1);
    arb_init(c0); arb_init(c1); arb_init(bs);
    aic_ecstar_defect_assoc_hopm(a0, &A, 1, 0, PREC);   /* weakest */
    aic_ecstar_defect_assoc_hopm(a1, &A, 12, 15, PREC); /* full    */
    aic_ecstar_defect_submult_hopm(s0, &A, 1, 0, PREC);
    aic_ecstar_defect_submult_hopm(s1, &A, 12, 15, PREC);
    aic_ecstar_defect_cstar_hopm(c0, &A, 1, 0, PREC);
    aic_ecstar_defect_cstar_hopm(c1, &A, 12, 15, PREC);
    aic_ecstar_defect_assoc(bs, &A, PREC);
    printf("SEARCHES assoc  weakest(1,0)=%.4e  full=%.4e  basis_sweep=%.4e\n",
           defect_dbl(a0), defect_dbl(a1), defect_dbl(bs));
    printf("SEARCHES submult weakest(1,0)=%.4e  full=%.4e\n", defect_dbl(s0), defect_dbl(s1));
    printf("SEARCHES cstar  weakest(1,0)=%.4e  full=%.4e\n", defect_dbl(c0), defect_dbl(c1));
    /* the full search STRICTLY beats the weakest -- a non-searching stub is RED */
    AIC_CHECK_MSG(defect_dbl(a1) > defect_dbl(a0) + 1e-5,
                  "searches: assoc full %.3e not > weakest %.3e (search is a stub?)",
                  defect_dbl(a1), defect_dbl(a0));
    AIC_CHECK_MSG(defect_dbl(s1) > defect_dbl(s0) + 1e-5,
                  "searches: submult full %.3e not > weakest %.3e",
                  defect_dbl(s1), defect_dbl(s0));
    AIC_CHECK_MSG(defect_dbl(c1) > defect_dbl(c0) + 1e-5,
                  "searches: cstar full %.3e not > weakest %.3e",
                  defect_dbl(c1), defect_dbl(c0));
    /* full assoc recovers at least the basis-sweep witness value (hand-checkable). */
    AIC_CHECK_MSG(defect_dbl(a1) >= defect_dbl(bs) - 1e-9,
                  "searches: full assoc %.3e below basis-sweep witness %.3e",
                  defect_dbl(a1), defect_dbl(bs));
    arb_clear(bs); arb_clear(c1); arb_clear(c0);
    arb_clear(s1); arb_clear(s0); arb_clear(a1); arb_clear(a0);
    aic_ecstar_clear(&A);
    aic_ucp_kraus_clear(&phit);
    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&dep);
    aic_ucp_kraus_clear(&phi);
}

/* E5: certified-ball soundness. The returned lo must be <= the assoc ratio
 * re-evaluated at the explicit witness (rigorous lower bound), and each witness
 * block must be in A (aic_ecstar_proj_residual ~0).
 *
 * FINDING (subspace-polar accept-guard, Cycle-2 blind spot). The witness-in-A
 * property here is enforced by Pi_A(polar(C)) in the block update. But on EVERY
 * Cycle-2 test instance, A = Img(exact Phi) is a GENUINE C* algebra, which is
 * polar/von-Neumann-CLOSED, so polar(C) already lands in A and dropping Pi_A is
 * observationally identical (verified: a skip-Pi_A mutation leaves this test
 * GREEN). The Pi_A projection's CORRECTNESS contribution (keeping witnesses in an
 * only-eps-polar-closed A) can only be exercised by a genuine eta>0 A = Img Phi~
 * (Phi~=theta(2Phi-1)) -- that is Cycle 3 (bead aic-92f assoc_ecsa), where the
 * associator is also nonzero so candidates actually get accepted. The accept-
 * guard's MONOTONE-ASCENT contribution IS exercised now: a never-iterate mutation
 * of the engine drives the canary/searches RED (the witness ratio collapses).
 * Bead filed for the Pi_A teeth on a non-polar-closed A. */
static void test_hopm_soundness(void)
{
    aic_ucp_kraus phi, dep;
    make_block_cond_exp(&phi, 4, 2);
    make_dep(&dep, 4);
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, PREC);
    aic_ucp_kraus phit;
    make_phi_mix(&phit, &phi, &dep, 1e-2, PREC);
    aic_ecstar A;
    aic_ecstar_init(&A, d.n, d.dim_A, &phit);
    aic_ecstar_from_idemp(&A, &d, &phit, PREC);

    slong n = A.n;
    acb_mat_t wX, wY, wZ, xy, lhs, yz, rhs, hh;
    acb_mat_init(wX, n, n); acb_mat_init(wY, n, n); acb_mat_init(wZ, n, n);
    acb_mat_init(xy, n, n); acb_mat_init(lhs, n, n); acb_mat_init(yz, n, n);
    acb_mat_init(rhs, n, n); acb_mat_init(hh, n, n);
    arb_t lo;
    arb_init(lo);
    aic_ecstar_defect_assoc_hopm_witness(lo, wX, wY, wZ, &A, 20, 20, PREC);

    /* (i) witnesses in A */
    arb_t rX, rY, rZ;
    arb_init(rX); arb_init(rY); arb_init(rZ);
    aic_ecstar_proj_residual(rX, &A, wX, PREC);
    aic_ecstar_proj_residual(rY, &A, wY, PREC);
    aic_ecstar_proj_residual(rZ, &A, wZ, PREC);
    AIC_CHECK_MSG(defect_dbl(rX) < 1e-9 && defect_dbl(rY) < 1e-9 &&
                  defect_dbl(rZ) < 1e-9,
                  "soundness: witness not in A (residuals %.2e %.2e %.2e)",
                  defect_dbl(rX), defect_dbl(rY), defect_dbl(rZ));

    /* (ii) lo <= ratio at the witness. ratio = ||h||_op/(||X||||Y||||Z||). */
    aic_ecstar_star(xy, &A, wX, wY, PREC);
    aic_ecstar_star(lhs, &A, xy, wZ, PREC);
    aic_ecstar_star(yz, &A, wY, wZ, PREC);
    aic_ecstar_star(rhs, &A, wX, yz, PREC);
    acb_mat_sub(hh, lhs, rhs, PREC);
    arb_t nh, nx, ny, nz, ratio;
    arb_init(nh); arb_init(nx); arb_init(ny); arb_init(nz); arb_init(ratio);
    aic_mat_opnorm(nh, hh, PREC);
    aic_mat_opnorm(nx, wX, PREC);
    aic_mat_opnorm(ny, wY, PREC);
    aic_mat_opnorm(nz, wZ, PREC);
    arb_mul(ratio, nx, ny, PREC);
    arb_mul(ratio, ratio, nz, PREC);
    arb_div(ratio, nh, ratio, PREC);
    printf("SOUNDNESS lo=%.6e  ratio@witness=%.6e  residuals<%.1e\n",
           defect_dbl(lo), defect_dbl(ratio),
           fmax(fmax(defect_dbl(rX), defect_dbl(rY)), defect_dbl(rZ)) + 1e-30);
    AIC_CHECK_MSG(arb_le(lo, ratio) || defect_dbl(lo) <= defect_dbl(ratio) + 1e-9,
                  "soundness: lo %.6e > ratio@witness %.6e (NOT a lower bound)",
                  defect_dbl(lo), defect_dbl(ratio));

    arb_clear(ratio); arb_clear(nz); arb_clear(ny); arb_clear(nx); arb_clear(nh);
    arb_clear(rZ); arb_clear(rY); arb_clear(rX);
    arb_clear(lo);
    acb_mat_clear(hh); acb_mat_clear(rhs); acb_mat_clear(yz); acb_mat_clear(lhs);
    acb_mat_clear(xy); acb_mat_clear(wZ); acb_mat_clear(wY); acb_mat_clear(wX);
    aic_ecstar_clear(&A);
    aic_ucp_kraus_clear(&phit);
    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&dep);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_oracle();
    teeth_phi_mix();
    teeth_submult();
    teeth_basis_off_A();
    teeth_non_hp();
    test_gauge();
    test_prec();
    /* Cycle 2: faithful worst-case HOPM search */
    test_hopm_oracle();
    test_hopm_vs_basis_sweep();
    test_hopm_canary();
    test_hopm_searches();
    test_hopm_soundness();
    aic_test_report("test_ecstar");
    printf("OK test_ecstar\n");
    return 0;
}
