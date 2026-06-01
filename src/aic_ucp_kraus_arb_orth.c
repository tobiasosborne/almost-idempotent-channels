/* aic_ucp_kraus_arb_orth.c — the keep-threshold ball and the per-cluster LOEWDIN
 * orthonormalisation for the certified Choi->Kraus extraction (bead aic-4td
 * increment 2 step D). Split out of aic_ucp_kraus_arb.c for Rule 10 (~200
 * LOC/file). The full rationale (Convention A, the R2 orthonormalisation
 * subtlety, the x0 = ||X||_op^2 / FINDINGS §C5 interval-Gram hazard) is in that
 * file's docstring; this file holds the two numerical primitives.
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_funcalc.h"
#include "aic/aic_mat.h"
#include "aic_ucp_kraus_arb_internal.h"

void aic_ucp_kraus_arb_int_keep_thr(arb_t thr, const acb_mat_t C, slong n,
                                    slong prec)
{
    aic_mat_frobenius_norm(thr, C, prec);
    arb_mul_2exp_si(thr, thr, -52);          /* * 2^-52 (DBL_EPSILON) */
    arb_mul_si(thr, thr, n, prec);           /* * (dim_K*dim_H) */
}

void aic_ucp_kraus_arb_int_loewdin(acb_mat_t V, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X), k = acb_mat_ncols(X);
    if (k == 1) {
        /* V = X / sqrt(X^dag X) (a single column normalisation). */
        arb_t g;
        arb_init(g);
        acb_t acc, t;
        acb_init(acc);
        acb_init(t);
        for (slong i = 0; i < n; i++) {
            acb_conj(t, acb_mat_entry(X, i, 0));
            acb_addmul(acc, t, acb_mat_entry(X, i, 0), prec); /* sum |X_i|^2 */
        }
        acb_get_real(g, acc);                 /* ||X||^2 (real, PD) */
        arb_sqrt(g, g, prec);                 /* ||X|| */
        arb_inv(g, g, prec);
        acb_mat_scalar_mul_arb(V, X, g, prec);
        acb_clear(t);
        acb_clear(acc);
        arb_clear(g);
        return;
    }
    /* G = X^dag X (k x k Hermitian PD); Ginvsqrt = G^{-1/2}; V = X Ginvsqrt. */
    acb_mat_t Xd, G, Ginvsqrt;
    acb_mat_init(Xd, k, n);
    acb_mat_init(G, k, k);
    acb_mat_init(Ginvsqrt, k, k);
    acb_mat_conjugate_transpose(Xd, X);
    acb_mat_mul(G, Xd, X, prec);
    /* x0 = lambda_max(G) = ||X||_op^2 makes ||G/x0 - I||_op = 1 - l_min/l_max < 1,
     * so the xpow binomial domain (tex:511) always holds for PD G. We get
     * lambda_max(G) as ||X||_op^2 via aic_mat_opnorm (which Weyl-Hermitianizes the
     * interval Gram internally, FINDINGS §C5) — NOT aic_mat_herm_max_eig(G), which
     * would ABORT on the raw interval Gram X^dag X (entry (i,j) and conj(j,i) carry
     * independent radii; the tight Hermiticity assert false-fails). */
    arb_t x0;
    arb_init(x0);
    aic_mat_opnorm(x0, X, prec);             /* ||X||_op */
    arb_sqr(x0, x0, prec);                    /* lambda_max(G) = ||X||_op^2 */
    double x0d = arf_get_d(arb_midref(x0), ARF_RND_UP) + mag_get_d(arb_radref(x0));
    aic_funcalc_xpow(Ginvsqrt, G, -0.5, x0d, prec);
    acb_mat_mul(V, X, Ginvsqrt, prec);
    arb_clear(x0);
    acb_mat_clear(Ginvsqrt);
    acb_mat_clear(G);
    acb_mat_clear(Xd);
}
