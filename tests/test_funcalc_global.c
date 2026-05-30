/* test_funcalc_global.c — cross-checks for the GLOBALLY-convergent non-normal
 * matrix-sign (bead aic-8hz). Split out of test_funcalc.c (which is at the ~200
 * LOC limit, Rule 10). CLAUDE.md Rule 5/6: every check asserts a defining
 * identity, a closed-form oracle, or an inter-candidate agreement.
 *
 * The deliverable relaxes funcalc's sgn precondition from the OPERATOR-NORM basin
 * ||X^2-I||_op < 1 (Newton-Schulz local convergence) to the SPECTRAL condition
 * rho(I-X^2) < 1 (X has no purely-imaginary eigenvalue), certified EIG-FREE via
 * the Gelfand bound. This file exercises the NON-NORMAL regime that the op-basin
 * path rejected.
 *
 * Checks:
 *  T-global-1 (teeth): non-normal X=[[a,c],[0,-b]], a=1.05,b=0.95,c=20. Then
 *    ||I-X^2||_op ~ 2 > 1 (op-basin EXCEEDED) but rho(I-X^2)=0.1025 < 1, and the
 *    EXACT sgn(X)=(2X-(a-b)I)/(a+b)=[[1,2c/(a+b)],[0,-1]]=[[1,20],[0,-1]]. Asserts
 *    the op-basin probe is FALSE, the global Newton returns the closed form to
 *    ~1e-60, aic_sgn (auto-dispatch) returns the same, and sgn^2=I, [X,sgn]=0.
 *  T-global-2 (eta=0 oracle): exact projector P, X=2P-I (X^2=I) -> sgn=X
 *    machine-zero; plus a non-normal near-sign perturbation off an exact projector.
 *  T-global-3 (in-basin agreement): on an in-basin X, newton_global / DB /
 *    Newton-Schulz all agree to ~1e-25, and aic_sgn dispatches to Newton-Schulz
 *    (output identical to the direct NS call).
 *  T-global-4 (precondition fail-loud / not-certified): a matrix with
 *    rho(I-X^2) >= 1 must NOT be certified by the Gelfand certifier (the
 *    non-aborting internal entry returns 0). The public wrapper aborts; we drive
 *    the internal entry to avoid killing the suite.
 *
 * Mutation-proof (Rule 7): see the per-test MUTATION notes; the teeth were
 * confirmed RED by (a) flipping the Newton-step sign to Y-Y^{-1} and (b) dropping
 * the 1/2 factor, restored after RED was observed (recorded in the bead report).
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_funcalc.h"
/* The non-aborting basin probe + Gelfand certifier are funcalc internals (the
 * public wrappers abort); the test drives them directly, like test_certify.c
 * reaches into src/aic_cbnorm_internal.h. */
#include "../src/aic_funcalc_internal.h"
#include "aic_mat.h"
#include "aic_test.h"

static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* X = [[a, c],[0, -b]] (upper-triangular, non-normal for c != 0). */
static void build_teeth(acb_mat_t X, double a, double b, double c)
{
    acb_mat_zero(X);
    acb_set_d(acb_mat_entry(X, 0, 0), a);
    acb_set_d(acb_mat_entry(X, 0, 1), c);
    acb_set_d(acb_mat_entry(X, 1, 1), -b);
}

/* ||I - X^2||_op and rho(I-X^2) reporting helper (rho via Gelfand internal). */
static void report_basin(const char *tag, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    acb_mat_t X2, M;
    arb_t opn, rho;
    slong kused = 0;
    acb_mat_init(X2, n, n);
    acb_mat_init(M, n, n);
    arb_init(opn);
    arb_init(rho);
    acb_mat_sqr(X2, X, prec);
    acb_mat_one(M);
    acb_mat_sub(M, M, X2, prec);          /* M = I - X^2 */
    aic_mat_opnorm(opn, M, prec);
    int cert = aic_funcalc_int_gelfand_rho(rho, &kused, M, 32, prec);
    printf("  %s: ||I-X^2||_op=%.6g  rho_ub(I-X^2)=%.6g (k=%ld, cert=%d)\n",
           tag, arf_get_d(arb_midref(opn), ARF_RND_NEAR),
           arf_get_d(arb_midref(rho), ARF_RND_NEAR), (long) kused, cert);
    arb_clear(rho);
    arb_clear(opn);
    acb_mat_clear(M);
    acb_mat_clear(X2);
}

/* --- T-global-1: the teeth (non-normal, op-basin exceeded) --- */
static void test_teeth(void)
{
    const slong prec = 256;
    const slong n = 2;
    const double a = 1.05, b = 0.95, c = 20.0;
    acb_mat_t X, S, Sauto, Sgn2, I, comm, XS, SX, Sref;
    acb_mat_init(X, n, n);
    acb_mat_init(S, n, n);
    acb_mat_init(Sauto, n, n);
    acb_mat_init(Sgn2, n, n);
    acb_mat_init(I, n, n);
    acb_mat_init(comm, n, n);
    acb_mat_init(XS, n, n);
    acb_mat_init(SX, n, n);
    acb_mat_init(Sref, n, n);
    build_teeth(X, a, b, c);

    printf("T-global-1 (teeth) a=%.2f b=%.2f c=%.1f:\n", a, b, c);
    report_basin("X", X, prec);

    /* (i) the OLD path's basin is EXCEEDED: the op-basin probe must be FALSE. */
    AIC_CHECK_MSG(aic_funcalc_int_in_op_basin(X, prec) == 0,
                  "T-global-1: op-basin probe true, but ||I-X^2||_op > 1 expected");

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-60);

    /* closed form sgn(X) = (2X - (a-b)I)/(a+b) = [[1, 2c/(a+b)],[0, -1]]. */
    acb_mat_zero(Sref);
    acb_set_d(acb_mat_entry(Sref, 0, 0), 1.0);
    acb_set_d(acb_mat_entry(Sref, 0, 1), 2.0 * c / (a + b));
    acb_set_d(acb_mat_entry(Sref, 1, 1), -1.0);

    /* (ii) global Newton returns the closed form. */
    aic_sgn_newton_global(S, X, prec);
    AIC_CHECK_ACB_MAT_CLOSE(S, Sref, tol);

    /* (iii) aic_sgn (auto-dispatch) returns the same. */
    aic_sgn(Sauto, X, prec);
    AIC_CHECK_ACB_MAT_CLOSE(Sauto, Sref, tol);

    /* (iv) sgn^2 = I and [X, sgn] = 0. */
    acb_mat_sqr(Sgn2, S, prec);
    acb_mat_one(I);
    AIC_CHECK_ACB_MAT_CLOSE(Sgn2, I, tol);
    acb_mat_mul(XS, X, S, prec);
    acb_mat_mul(SX, S, X, prec);
    acb_mat_sub(comm, XS, SX, prec);
    arb_t cnorm;
    arb_init(cnorm);
    aic_mat_opnorm(cnorm, comm, prec);
    AIC_CHECK_MSG(arb_lt(cnorm, tol), "T-global-1: [X, sgn(X)] not within tol");
    arb_clear(cnorm);

    arb_clear(tol);
    acb_mat_clear(Sref); acb_mat_clear(SX); acb_mat_clear(XS);
    acb_mat_clear(comm); acb_mat_clear(I); acb_mat_clear(Sgn2);
    acb_mat_clear(Sauto); acb_mat_clear(S); acb_mat_clear(X);
}

/* --- T-global-2: eta=0 oracle (exact sign matrix) + non-normal near-sign --- */
static void test_eta0_oracle(void)
{
    const slong prec = 256;
    const slong n = 5;
    const slong rank = 3;
    acb_mat_t P, X, S, I;
    acb_mat_init(P, n, n);
    acb_mat_init(X, n, n);
    acb_mat_init(S, n, n);
    acb_mat_init(I, n, n);
    acb_mat_zero(P);
    for (slong i = 0; i < rank; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
    acb_mat_scalar_mul_2exp_si(X, P, 1);   /* 2P */
    acb_mat_one(I);
    acb_mat_sub(X, X, I, prec);            /* X = 2P - I, X^2 = I exactly */

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-60);
    aic_sgn_newton_global(S, X, prec);
    AIC_CHECK_ACB_MAT_CLOSE(S, X, tol);   /* sgn(2P-I) = 2P-I */

    /* non-normal near-sign: add a strictly-upper-triangular nilpotent N (commutes
     * with nothing, makes X non-normal) but keep eigenvalues exactly +/-1 so the
     * closed-form sign is the diagonal-of-signs read through the same poly. We use
     * a small upper c on a 2x2 sign block embedded; sgn must still recover the
     * eigen-sign structure. Reuse the teeth closed form on a=b=1 (X^2=I exactly
     * when a=b=1: X=[[1,c],[0,-1]], X^2=[[1, c-c],[0,1]]=I), sgn(X)=X itself. */
    acb_mat_t Y, Sy;
    acb_mat_init(Y, 2, 2);
    acb_mat_init(Sy, 2, 2);
    build_teeth(Y, 1.0, 1.0, 7.0);        /* X^2 = I exactly, sgn(Y)=Y */
    aic_sgn_newton_global(Sy, Y, prec);
    AIC_CHECK_ACB_MAT_CLOSE(Sy, Y, tol);
    acb_mat_clear(Sy); acb_mat_clear(Y);

    arb_clear(tol);
    acb_mat_clear(I); acb_mat_clear(S); acb_mat_clear(X); acb_mat_clear(P);
}

/* --- T-global-3: in-basin agreement + aic_sgn dispatches to Newton-Schulz --- */
static void test_in_basin_agreement(void)
{
    const slong prec = 256;
    const slong n = 4;
    int signs[4] = {1, -1, 1, -1};
    acb_mat_t X, Sg, Sdb, Sns, Sauto;
    acb_mat_init(X, n, n);
    acb_mat_init(Sg, n, n);
    acb_mat_init(Sdb, n, n);
    acb_mat_init(Sns, n, n);
    acb_mat_init(Sauto, n, n);
    /* a Hermitian near-sign with ||X^2-I||_op < 1 (in basin) */
    acb_mat_zero(X);
    for (slong i = 0; i < n; i++) acb_set_si(acb_mat_entry(X, i, i), signs[i]);
    for (slong i = 0; i < n; i++)
        for (slong j = i + 1; j < n; j++) {
            acb_t z;
            acb_init(z);
            acb_set_d_d(z, 0.06 * (((i + 2 * j) % 5) - 2),
                        0.06 * (((3 * i + j) % 5) - 2));
            acb_set(acb_mat_entry(X, i, j), z);
            acb_conj(acb_mat_entry(X, j, i), z);
            acb_clear(z);
        }

    AIC_CHECK_MSG(aic_funcalc_int_in_op_basin(X, prec) == 1,
                  "T-global-3: in-basin probe false on a small Hermitian perturbation");

    aic_sgn_newton_global(Sg, X, prec);
    aic_sgn_denman_beavers(Sdb, X, prec);
    aic_sgn_newton_schulz(Sns, X, prec);
    aic_sgn(Sauto, X, prec);

    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-25);
    AIC_CHECK_ACB_MAT_CLOSE(Sg, Sdb, tol);
    AIC_CHECK_ACB_MAT_CLOSE(Sg, Sns, tol);
    /* dispatch: aic_sgn must take the Newton-Schulz path here (in basin), so its
     * output is BYTE-FOR-BYTE the direct NS call — acb_mat_equal is exact
     * structural equality (midpoint AND radius identical), the right check for
     * "same computation" (a tol of 0 fails even on identical balls, since the
     * subtraction [m±r]-[m±r] = [0±2r] has nonzero width). */
    AIC_CHECK_MSG(acb_mat_equal(Sauto, Sns),
                  "T-global-3: aic_sgn did not dispatch byte-for-byte to "
                  "Newton-Schulz in basin");
    arb_clear(tol);

    acb_mat_clear(Sauto); acb_mat_clear(Sns); acb_mat_clear(Sdb);
    acb_mat_clear(Sg); acb_mat_clear(X);
}

/* --- T-global-4: Gelfand certifier NOT-certified on rho(I-X^2) >= 1 --- */
static void test_gelfand_not_certified(void)
{
    const slong prec = 256;
    const slong n = 2;
    /* X = diag(eps, -eps) with eps tiny: X^2 = diag(eps^2, eps^2) ~ 0, so
     * I - X^2 ~ I, rho(I-X^2) ~ 1 - eps^2 < 1 but the eigenvalue of X sits NEAR
     * the imaginary axis (small |Re|). For rho >= 1 take eps = 0: X = 0, then
     * I - X^2 = I, rho = 1 exactly -> NOT certifiable (Gelfand bound floors at 1).
     */
    acb_mat_t X, X2, M;
    arb_t rho;
    slong kused = -1;
    acb_mat_init(X, n, n);
    acb_mat_init(X2, n, n);
    acb_mat_init(M, n, n);
    arb_init(rho);
    acb_mat_zero(X);                       /* X = 0 -> I - X^2 = I, rho = 1 */
    acb_mat_sqr(X2, X, prec);
    acb_mat_one(M);
    acb_mat_sub(M, M, X2, prec);
    int cert = aic_funcalc_int_gelfand_rho(rho, &kused, M, 32, prec);
    printf("T-global-4: X=0 -> rho_ub(I-X^2)=%.6g cert=%d (expect not-certified)\n",
           arf_get_d(arb_midref(rho), ARF_RND_NEAR), cert);
    AIC_CHECK_MSG(cert == 0,
                  "T-global-4: Gelfand certified rho(I)<1 — must be NOT certified");

    arb_clear(rho);
    acb_mat_clear(M); acb_mat_clear(X2); acb_mat_clear(X);
}

int main(void)
{
    test_teeth();
    test_eta0_oracle();
    test_in_basin_agreement();
    test_gelfand_not_certified();

    aic_test_report("test_funcalc_global");
    printf("OK test_funcalc_global\n");
    return 0;
}
