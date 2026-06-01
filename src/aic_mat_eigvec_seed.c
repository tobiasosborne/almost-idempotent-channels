/* aic_mat_eigvec_seed.c — the densify-assert + per-cluster Rump certification
 * helpers for aic_mat_eig_hermitian_subspaces (bead aic-4td increment 2 step C2;
 * design docs/research/eigvec_certified_design.md §1.2/§1.6). Split out of
 * aic_mat_eigvec.c (Rule 10, ~200 LOC/file); the orchestration (densify, seed,
 * cluster, disjointness gate) stays there, the two fail-loud-bearing primitives
 * live here.
 *
 *   aic_mat_int_assert_densify_unitary — ASSERT ||U U^dag - I|| certified tiny
 *     before trusting A' = U H U^dag (a loose U silently moves the spectrum;
 *     mirrors the C1 pattern in aic_mat_eig_multiple.c).
 *   aic_mat_int_certify_cluster — build the Rump seed Xa from the zheev cluster
 *     columns, run acb_mat_eig_enclosure_rump on A', ASSERT a FINITE + REAL
 *     enclosure (the Krawczyk certificate; design §1.6(i)), back-map
 *     X_c = U^dag X'_c (design §1.6(ii): the SAME J_c works for the original H),
 *     and ASSERT the SELF-CERTIFYING residual ||H X_c - X_c J_c|| on the ORIGINAL
 *     H is small (the honest certificate; design §1.6 names this residual as THE
 *     certificate — via aic_mat_int_assert_subspace_residual in the sibling
 *     aic_mat_eigvec_resid.c, split out for Rule 10 ~200 LOC/file).
 */
#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic_mat_internal.h"

void aic_mat_int_assert_densify_unitary(const acb_mat_t U, const acb_mat_t Ud,
                                        slong n, slong prec)
{
    acb_mat_t T;
    acb_mat_init(T, n, n);
    acb_mat_mul(T, U, Ud, prec);
    for (slong i = 0; i < n; i++)
        acb_sub_si(acb_mat_entry(T, i, i), acb_mat_entry(T, i, i), 1, prec);
    arb_t uu, utol;
    arb_init(uu);
    arb_init(utol);
    acb_mat_frobenius_norm(uu, T, prec);
    /* Tolerance n^2 * 2^-(prec-8). The densifier chains n(n-1)/2 Givens products
     * (aic_mat_densify.c), so arb's certified ||U U^dag - I||_F radius ACCUMULATES
     * ~ n^2 * 2^-prec (measured: 3.5e-38 at n=2 ... 3.4e-35 at n=12, prec=128).
     * The bare 2^-(prec-8) floor is exceeded for n >= 6, so we scale it by n^2 —
     * still ~4 orders below the certified eigenvalue ball radii (~1e-31), so the
     * conjugation genuinely preserves the spectrum; a broken U (defect ~1) still
     * fails loud. (design §1.3; FINDINGS §D5n2 — the C1 retry guard had the same
     * bare floor, a latent fail-loud on legitimate n>=6 inputs, fixed there too.) */
    aic_mat_int_tol(utol, prec);
    arb_mul_si(utol, utol, n * n, prec);
    if (!arb_lt(uu, utol)) {
        char *s = arb_get_str(uu, 8, 0);
        fprintf(stderr,
                "aic_mat_eig_hermitian_subspaces: densifier U not certified "
                "unitary (||U U^dag - I||_F=%s not < n^2*2^-(prec-8)) at prec=%ld "
                "— the conjugation would not preserve the spectrum; bead aic-4td\n",
                s ? s : "?", (long) prec);
        flint_free(s);
        abort();
    }
    arb_clear(utol);
    arb_clear(uu);
    acb_mat_clear(T);
}

void aic_mat_int_certify_cluster(aic_mat_eigcluster *out, const acb_mat_t H,
                                 const acb_mat_t A1, const acb_mat_t Ud,
                                 const double _Complex *Vd, const double *ev,
                                 slong n, slong s0, slong k, slong cidx, slong prec)
{
    acb_mat_t Xa;
    acb_mat_init(Xa, n, k);
    for (slong col = 0; col < k; col++)
        for (slong row = 0; row < n; row++) {
            double _Complex z = Vd[row * n + (s0 + col)];
            acb_set_d_d(acb_mat_entry(Xa, row, col), creal(z), cimag(z));
        }
    double lam_mid = 0.0;
    for (slong t = s0; t < s0 + k; t++) lam_mid += ev[t];
    lam_mid /= (double) k;
    acb_t lam_approx;
    acb_init(lam_approx);
    acb_set_d(lam_approx, lam_mid);

    acb_init(out->lambda);
    acb_mat_init(out->X, n, k);
    acb_mat_init(out->J, k, k);
    out->k = k;

    acb_mat_t Xp;
    acb_mat_init(Xp, n, k);
    acb_mat_eig_enclosure_rump(out->lambda, out->J, Xp, A1, lam_approx, Xa, prec);

    /* (i) FINITE enclosure is the Rump certificate — fail loud otherwise (Rule 4,
     * NDEBUG-immune fprintf+abort). */
    if (!acb_is_finite(out->lambda) || !acb_mat_is_finite(Xp)) {
        fprintf(stderr,
                "aic_mat_eig_hermitian_subspaces: cluster %ld (k=%ld) UNRESOLVED "
                "at prec=%ld — Rump enclosure non-finite even after dense-unitary "
                "densification (raise prec); FINDINGS §D5n, bead aic-4td\n",
                (long) cidx, (long) k, (long) prec);
        abort();
    }
    /* Hermitian => real spectrum: imag(lambda) must enclose 0. */
    arb_t im;
    arb_init(im);
    acb_get_imag(im, out->lambda);
    if (!arb_contains_zero(im)) {
        fprintf(stderr,
                "aic_mat_eig_hermitian_subspaces: cluster %ld certified non-real "
                "(imag enclosure excludes 0) on a Hermitian input\n", (long) cidx);
        abort();
    }
    arb_clear(im);

    acb_mat_mul(out->X, Ud, Xp, prec);        /* back-map X_c = U^dag X'_c */

    /* (iii) the honest, self-certifying residual on the ORIGINAL H (design §1.6;
     * FINDINGS §D5n). Production certifies the Rump enclosure on A'; this closes
     * the previously-uncertified gap by proving H X_c = X_c J_c on H directly. */
    aic_mat_int_assert_subspace_residual(H, out, n, cidx, prec);

    acb_mat_clear(Xp);
    acb_clear(lam_approx);
    acb_mat_clear(Xa);
}
