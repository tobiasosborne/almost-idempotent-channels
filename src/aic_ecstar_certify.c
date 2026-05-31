/* aic_ecstar_certify.c — public API + arb CERTIFICATION for the Cycle-2 faithful
 * worst-case defect search (bead aic-knm, module "ecstar"):
 * aic_ecstar_defect_{assoc,submult,cstar}_hopm (include/aic_ecstar.h). It runs the
 * fast double engine (aic_ecstar_search.c, set up by aic_ecstar_setup.c) to FIND a
 * worst-case witness, then re-evaluates the explicit witness IN A in the certified
 * arb path (aic_ecstar_star + aic_mat_opnorm) to return a RIGOROUS LOWER BOUND
 * `lo` on the true sup. The lower bound uses the arb LOWER end of the numerator
 * norm and the arb UPPER end of the denominator norms, so the quotient is provably
 * <= the value at the witness <= the sup (web leg Q8). See include/aic_ecstar.h.
 */
#include <assert.h>
#include <complex.h>

#include <flint/arf.h>

#include "aic/aic_ecstar.h"
#include "aic_ecstar_setup.h"
#include "aic/aic_latd.h"
#include "aic/aic_mat.h"

/* load a double witness block (n*n) into an acb_mat (zero-radius). */
static void cx_to_acb(acb_mat_t M, const cx *X, slong n)
{ aic_latd_to_acb_mat(M, X, n, n); }

/* lo/hi end of the certified ||M||_op ball as an arf. */
static void opnorm_lo(arf_t lo, const acb_mat_t M, slong prec)
{ arb_t e; arb_init(e); aic_mat_opnorm(e, M, prec); arb_get_lbound_arf(lo, e, prec); arb_clear(e); }
static void opnorm_hi(arf_t hi, const acb_mat_t M, slong prec)
{ arb_t e; arb_init(e); aic_mat_opnorm(e, M, prec); arb_get_ubound_arf(hi, e, prec); arb_clear(e); }

/* ---- public API ---------------------------------------------------------- */

/* shared assoc body; optionally copies the witness into wX,wY,wZ (any NULL). */
static void certify_assoc(arb_t lo, acb_mat_t wX, acb_mat_t wY, acb_mat_t wZ,
                          const aic_ecstar *A, int n_restarts, int n_iters,
                          slong prec)
{
    assert(n_restarts >= 1 && n_iters >= 0);
    slong n = A->n, nn = n * n;
    aic_ehk h; aic_ehk_snapshot(&h, A);
    cx *X = flint_malloc((size_t) nn * sizeof(cx));
    cx *Y = flint_malloc((size_t) nn * sizeof(cx));
    cx *Z = flint_malloc((size_t) nn * sizeof(cx));
    aic_ehk_run_search(&h, n_restarts, n_iters, 3, 0, aic_ehk_term_assoc_X,
                       aic_ehk_term_assoc_Y, aic_ehk_term_assoc_Z, X, Y, Z);
    acb_mat_t aX, aY, aZ, hh, xy, lhs, yz, rhs;
    acb_mat_init(aX,n,n); acb_mat_init(aY,n,n); acb_mat_init(aZ,n,n); acb_mat_init(hh,n,n);
    acb_mat_init(xy,n,n); acb_mat_init(lhs,n,n); acb_mat_init(yz,n,n); acb_mat_init(rhs,n,n);
    cx_to_acb(aX,X,n); cx_to_acb(aY,Y,n); cx_to_acb(aZ,Z,n);
    if (wX) acb_mat_set(wX, aX);
    if (wY) acb_mat_set(wY, aY);
    if (wZ) acb_mat_set(wZ, aZ);
    aic_ecstar_star(xy, A, aX, aY, prec);
    aic_ecstar_star(lhs, A, xy, aZ, prec);   /* (X*Y)*Z */
    aic_ecstar_star(yz, A, aY, aZ, prec);
    aic_ecstar_star(rhs, A, aX, yz, prec);   /* X*(Y*Z) */
    acb_mat_sub(hh, lhs, rhs, prec);
    arf_t num, dx, dy, dz, den;
    arf_init(num); arf_init(dx); arf_init(dy); arf_init(dz); arf_init(den);
    opnorm_lo(num, hh, prec);
    opnorm_hi(dx, aX, prec); opnorm_hi(dy, aY, prec); opnorm_hi(dz, aZ, prec);
    arf_mul(den, dx, dy, prec, ARF_RND_UP); arf_mul(den, den, dz, prec, ARF_RND_UP);
    arf_div(num, num, den, prec, ARF_RND_DOWN);
    arb_set_arf(lo, num);
    arf_clear(num); arf_clear(dx); arf_clear(dy); arf_clear(dz); arf_clear(den);
    acb_mat_clear(rhs); acb_mat_clear(yz); acb_mat_clear(lhs); acb_mat_clear(xy);
    acb_mat_clear(hh); acb_mat_clear(aZ); acb_mat_clear(aY); acb_mat_clear(aX);
    flint_free(Z); flint_free(Y); flint_free(X); aic_ehk_snapshot_free(&h);
}

void aic_ecstar_defect_assoc_hopm(arb_t lo, const aic_ecstar *A, int n_restarts,
                                  int n_iters, slong prec)
{ certify_assoc(lo, NULL, NULL, NULL, A, n_restarts, n_iters, prec); }

void aic_ecstar_defect_assoc_hopm_witness(arb_t lo, acb_mat_t wX, acb_mat_t wY,
                                          acb_mat_t wZ, const aic_ecstar *A,
                                          int n_restarts, int n_iters, slong prec)
{ certify_assoc(lo, wX, wY, wZ, A, n_restarts, n_iters, prec); }

void aic_ecstar_defect_submult_hopm(arb_t lo, const aic_ecstar *A, int n_restarts,
                                    int n_iters, slong prec)
{
    assert(n_restarts >= 1 && n_iters >= 0);
    slong n = A->n, nn = n * n;
    aic_ehk h; aic_ehk_snapshot(&h, A);
    cx *X = flint_malloc((size_t) nn * sizeof(cx));
    cx *Y = flint_malloc((size_t) nn * sizeof(cx));
    aic_ehk_run_search(&h, n_restarts, n_iters, 2, 0, aic_ehk_term_sub_X,
                       aic_ehk_term_sub_Y, NULL, X, Y, NULL);
    acb_mat_t aX, aY, p;
    acb_mat_init(aX,n,n); acb_mat_init(aY,n,n); acb_mat_init(p,n,n);
    cx_to_acb(aX,X,n); cx_to_acb(aY,Y,n);
    aic_ecstar_star(p, A, aX, aY, prec);
    arf_t num, dx, dy, den;
    arf_init(num); arf_init(dx); arf_init(dy); arf_init(den);
    opnorm_lo(num, p, prec);
    opnorm_hi(dx, aX, prec); opnorm_hi(dy, aY, prec);
    arf_mul(den, dx, dy, prec, ARF_RND_UP);
    arf_div(num, num, den, prec, ARF_RND_DOWN);  /* ratio.lo */
    arf_sub_si(num, num, 1, prec, ARF_RND_DOWN); /* ratio - 1 */
    arb_set_arf(lo, num);
    arf_clear(num); arf_clear(dx); arf_clear(dy); arf_clear(den);
    acb_mat_clear(p); acb_mat_clear(aY); acb_mat_clear(aX);
    flint_free(Y); flint_free(X); aic_ehk_snapshot_free(&h);
}

void aic_ecstar_defect_cstar_hopm(arb_t lo, const aic_ecstar *A, int n_restarts,
                                  int n_iters, slong prec)
{
    assert(n_restarts >= 1 && n_iters >= 0);
    slong n = A->n, nn = n * n;
    aic_ehk h; aic_ehk_snapshot(&h, A);
    cx *X = flint_malloc((size_t) nn * sizeof(cx));
    aic_ehk_run_search(&h, n_restarts, n_iters, 1, 1 /*minimize*/,
                       aic_ehk_term_cstar, NULL, NULL, X, NULL, NULL);
    /* eps_cstar lower bound = 1 - ratio.hi, ratio = ||X^dag*X||_op / ||X||_op^2;
     * ratio.hi = (||X^dag*X||_op).hi / (||X||_op).lo^2 (denom rounded DOWN). */
    acb_mat_t aX, Xd, g;
    acb_mat_init(aX,n,n); acb_mat_init(Xd,n,n); acb_mat_init(g,n,n);
    cx_to_acb(aX,X,n);
    acb_mat_conjugate_transpose(Xd, aX);
    aic_ecstar_star(g, A, Xd, aX, prec);     /* X^dag * X */
    arf_t num, dlo, den, ratio, one;
    arf_init(num); arf_init(dlo); arf_init(den); arf_init(ratio); arf_init(one);
    opnorm_hi(num, g, prec);
    opnorm_lo(dlo, aX, prec);
    arf_mul(den, dlo, dlo, prec, ARF_RND_DOWN);
    arf_div(ratio, num, den, prec, ARF_RND_UP);
    arf_set_si(one, 1);
    arf_sub(ratio, one, ratio, prec, ARF_RND_DOWN); /* 1 - ratio.hi */
    arb_set_arf(lo, ratio);
    arf_clear(one); arf_clear(ratio); arf_clear(den); arf_clear(dlo); arf_clear(num);
    acb_mat_clear(g); acb_mat_clear(Xd); acb_mat_clear(aX);
    flint_free(X); aic_ehk_snapshot_free(&h);
}
