/* aic_funcalc_global.c — the globally-convergent non-normal matrix-sign (bead
 * aic-8hz). Relaxes the sgn precondition from the OPERATOR-NORM Newton-Schulz
 * basin ||X^2-I||_op < 1 to the SPECTRAL condition rho(I-X^2) < 1 (X has no
 * purely-imaginary eigenvalue), certified EIG-FREE.
 *
 * approximate_algebras.tex:514-516 — definition:
 *     |X|=(X^2)^{1/2},  sgn(X)=X(X^2)^{-1/2}   (||X^2 - x0^2 I|| < x0^2).
 * approximate_algebras.tex:518-521 — sgn is the nearest exact solution of X^2=I:
 *     sgn(X)^2 = I for all X, and ||sgn(X)-X|| <= ||X|| O(||X^2-I||) if ||X^2-I||<1.
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). tex:520-521 states the
 * bound under the NORM hypothesis ||X^2-I||<1; the underlying functional calculus
 * (tex:546-550) only needs spec(X^2) off (-inf,0], i.e. X with no purely-imaginary
 * eigenvalue. The constructive route is Higham, *Functions of Matrices* (SIAM
 * 2008), Ch.5 Thm 5.6: the Newton sign iteration Y_0=X, Y_{k+1}=1/2(Y_k+Y_k^{-1})
 * converges GLOBALLY and quadratically to sgn(X) under exactly that spectral
 * hypothesis (NOT a norm bound). This is the aic_sgn_denman_beavers step; here it
 * is gated by a SPECTRAL precondition instead of the op-norm basin.
 *
 * THE IMPLEMENTED PRECONDITION is the x0=1 normalized basin rho(I-X^2) < 1 — a
 * SUFFICIENT sub-case of Higham's full no-imaginary-axis hypothesis, NOT the whole
 * of it (e.g. X=diag(2,-3) has spec(X^2)={4,9} off (-inf,0], so sgn is Higham-well-
 * defined, yet rho(I-X^2)=8 and this routine ABORTS). That is by design: our inputs
 * are X=2 S_Phi - 1 with spec clustered near +/-1 (rho(I-X^2)=4 rho(S_Phi^2-S_Phi)
 * = 4 eta < 1), exactly the disk this basin certifies, and the Gelfand bound below
 * is cheap and eig-free there. A general no-imaginary-axis sgn would need an
 * eigenvalue/Schur route (aic-w4o.1), out of scope here.
 *
 * GELFAND CERTIFIER (eig-free, rigorous). After this module's x0=1 normalization,
 * rho(I-X^2) < 1 ==> spec(X^2) in the disk |z-1|<1 subset {Re>0}, so X has no
 * purely-imaginary eigenvalue (this one direction is what soundness uses; the
 * converse fails, as the diag(2,-3) example shows). We certify
 * rho(M)<1 for M=I-X^2 via the Gelfand bound: ||M^k||^{1/k} is a rigorous upper
 * bound on rho(M) for every k and decreases to rho(M), so a single k with
 * ||M^k||_F^{1/k} CERTAINLY < 1 proves rho(M)<1 — no eigensolver (sidesteps the
 * degenerate-eig debt aic-w4o.1). k=1 is the old Frobenius check; non-normal M
 * needs a few powers past the transient. We scan k=1..k_max (k_max=32) and accept
 * the FIRST certified k. arb_lt is true only if the WHOLE ball < 1, so a straddle
 * is not accepted (the rigor) — try a larger k, else fail loud.
 *
 * A-POSTERIORI CERTIFICATE. After convergence the result Y is certified an
 * involution commuting with X: ||Y^2-I|| < tol AND ||XY-YX|| < tol (certified
 * balls; a straddle aborts, Rule 4). With the certified precondition and Y_0=X,
 * Higham Thm 5.6 pins Y=sgn(X) (rules out -sgn and wrong-inertia involutions).
 * Factored into the shared aic_funcalc_int_certify_sign (bead aic-1vp), whose
 * uses a prec-floor sanity-backstop tol (see aic_funcalc_int_certify_sign; soundness is the basin precond + mid(X) seed, not tol) so the certificate does not stall
 * below the radius floor on a wide-radius input.
 *
 * WIDE-RADIUS MIDPOINT FIX (bead aic-1vp, FINDINGS §C11). Like Newton-Schulz, the
 * convergence loop runs on mid(X) (radius stripped once): an inherited input radius
 * (~1e-70) propagates through the per-step acb_mat_inv/add, drags the iterate's
 * midpoint off the involution, and the midpoint step never falls below the prec-tol
 * floor — the 100-iter cap then fires on an in-basin X. On mid(X) the flow is clean;
 * rigor is restored by the a-posteriori certificate against the ORIGINAL X. The
 * Gelfand precondition is still certified on the original X (rigorous, radius-aware).
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/mag.h>

#include "aic/aic_mat.h"
#include "aic/aic_funcalc.h"
#include "aic_funcalc_internal.h"

#define AIC_SGN_MAX_ITERS 100
#define AIC_GELFAND_KMAX 32

int aic_funcalc_int_in_op_basin(const acb_mat_t X, slong prec)
{
    arb_t bound, one;
    arb_init(bound);
    arb_init(one);
    arb_one(one);
    aic_funcalc_int_def_X2(bound, X, prec);   /* certified ||X^2 - I||_op */
    int in = arb_lt(bound, one);              /* whole ball < 1 (rigorous) */
    arb_clear(one);
    arb_clear(bound);
    return in;
}

int aic_funcalc_int_gelfand_rho(arb_t rho_ub, slong *k_used, const acb_mat_t M,
                                slong k_max, slong prec)
{
    slong n = acb_mat_nrows(M);
    assert(acb_mat_ncols(M) == n && "gelfand: M must be square");
    assert(k_max >= 1 && "gelfand: k_max must be >= 1");

    acb_mat_t Mk;            /* M^k, built incrementally */
    arb_t nrm, one;
    acb_mat_init(Mk, n, n);
    arb_init(nrm);
    arb_init(one);
    arb_one(one);
    acb_mat_set(Mk, M);      /* Mk = M^1 */

    int certified = 0;
    for (slong k = 1; k <= k_max; k++) {
        aic_mat_frobenius_norm(nrm, Mk, prec);     /* ||M^k||_F (certified) */
        arb_root_ui(rho_ub, nrm, (ulong) k, prec); /* ||M^k||_F^{1/k} >= rho(M) */
        if (arb_lt(rho_ub, one)) {                 /* whole ball < 1 -> rho < 1 */
            if (k_used) *k_used = k;
            certified = 1;
            break;
        }
        if (k < k_max) acb_mat_mul(Mk, Mk, M, prec);  /* M^{k+1} */
        if (k_used) *k_used = k;                       /* tightest k tried */
    }

    arb_clear(one);
    arb_clear(nrm);
    acb_mat_clear(Mk);
    return certified;
}

void aic_sgn_newton_global(acb_mat_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(acb_mat_ncols(X) == n && "aic_sgn_newton_global: X must be square");
    assert(out != X && "aic_sgn_newton_global: out must not alias X");

    /* Precondition: certify rho(I - X^2) < 1 via Gelfand (NOT the op-basin). */
    acb_mat_t X2, M;
    arb_t rho;
    acb_mat_init(X2, n, n);
    acb_mat_init(M, n, n);
    arb_init(rho);
    acb_mat_sqr(X2, X, prec);
    acb_mat_one(M);
    acb_mat_sub(M, M, X2, prec);             /* M = I - X^2 */
    if (!aic_funcalc_int_gelfand_rho(rho, NULL, M, AIC_GELFAND_KMAX, prec)) {
        fprintf(stderr, "aic_sgn_newton_global: cannot certify rho(I-X^2) < 1 at "
                "prec=%ld, k_max=%d (raise prec, or the input is at/over the "
                "spectral boundary — X has a near-imaginary eigenvalue)\n",
                (long) prec, AIC_GELFAND_KMAX);
        abort();
    }
    acb_mat_clear(M);
    acb_mat_clear(X2);

    arb_t tol, step;
    arb_init(tol);
    arb_init(step);
    aic_funcalc_int_tol(tol, prec);

    /* Iterate on mid(X) (bead aic-1vp); strip the inherited radius once at entry. */
    mag_t in_rad;
    mag_init(in_rad);

    acb_mat_t Y, Yinv, Ynext;
    acb_mat_init(Y, n, n);
    acb_mat_init(Yinv, n, n);
    acb_mat_init(Ynext, n, n);
    aic_funcalc_int_to_midpoint(Y, X, in_rad); /* Y_0 = mid(X) */

    int k;
    for (k = 0; k < AIC_SGN_MAX_ITERS; k++) {
        if (!acb_mat_inv(Yinv, Y, prec)) {   /* certified inverse or fail loud */
            fprintf(stderr, "aic_sgn_newton_global: could not certify the inverse "
                    "at iter %d, prec=%ld (near-singular iterate; raise prec)\n",
                    k, (long) prec);
            abort();
        }
        acb_mat_add(Ynext, Y, Yinv, prec);            /* Y + Y^{-1} */
        acb_mat_scalar_mul_2exp_si(Ynext, Ynext, -1); /* * 1/2 */
        aic_funcalc_int_step_norm(step, Ynext, Y, prec);
        acb_mat_set(Y, Ynext);
        if (arb_lt(step, tol)) { k++; break; }
    }
    if (k >= AIC_SGN_MAX_ITERS) {
        fprintf(stderr, "aic_sgn_newton_global: no convergence in %d iters "
                "(rho(I-X^2)<1 certified; should converge quadratically)\n",
                AIC_SGN_MAX_ITERS);
        abort();
    }

    /* A-posteriori certificate (shared, bead aic-1vp): ||Y^2-I|| AND ||XY-YX|| < tol
     * against the ORIGINAL X, the prec-floor sanity-backstop tol. With the certified Gelfand
     * precondition and Y_0=mid(X), Higham Thm 5.6 pins Y=sgn(X). */
    aic_funcalc_int_certify_sign(Y, X, in_rad, "aic_sgn_newton_global", prec);
    acb_mat_set(out, Y);

    mag_clear(in_rad);
    acb_mat_clear(Ynext);
    acb_mat_clear(Yinv);
    acb_mat_clear(Y);
    arb_clear(step);
    arb_clear(tol);
    arb_clear(rho);
}
