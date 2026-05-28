/* aic_mat_spectral.c — shared spectral internals + the degeneracy-robust largest
 * Hermitian eigenvalue (beads aic-9zh, "T-mat"). The full Hermitian
 * eigendecomposition (with eigenvectors and sorting) lives in aic_mat_eig.c; the
 * Hermiticity precondition and the certified two-stage eig path are factored
 * here (declared in src/aic_mat_internal.h) so both files share one copy.
 *
 * The certified eig wrinkle (FLINT 3.0.1). There is NO Hermitian-specialised
 * eigensolver; the only certified path is the general two-stage one:
 *   (1) acb_mat_approx_eig_qr  -> heuristic float approximations, NO error
 *       guarantee (acb_mat.rst:645-651);
 *   (2) acb_mat_eig_simple     -> PROVES n simple eigenvalues, returns certified
 *       isolating intervals E and right-eigenvector enclosures R with L = R^{-1},
 *       so L A R = diag(E) (acb_mat.rst:740-756).
 * aic_mat_int_eig_certified packages both; it aborts (Rule 4) if the spectrum
 * cannot be isolated (a simple-spectrum requirement; the cluster fallback
 * acb_mat_eig_multiple is a separate Law-4 audition).
 *
 * aic_mat_herm_max_eig takes a DIFFERENT, degeneracy-robust route because the
 * operator norm must work on matrices with repeated eigenvalues (the identity's
 * Gram matrix is I, a triply-or-more-degenerate spectrum that eig_simple cannot
 * isolate). It uses acb_mat_eig_global_enclosure (acb_mat.rst:653-694), which
 * gives a rigorous radius eps with every eigenvalue within eps of some
 * approximate eigenvalue — valid regardless of multiplicity — and returns
 * [max_k Re(E_k) - eps, max_k Re(E_k) + eps] as a certified ball on lambda_max.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/mag.h>

#include "aic_mat.h"
#include "aic_mat_internal.h"

void aic_mat_int_tol(arb_t tol, slong prec)
{
    arb_set_si(tol, 1);
    arb_mul_2exp_si(tol, tol, -(prec - 8));
}

void aic_mat_int_assert_hermitian(const acb_mat_t H, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(acb_mat_ncols(H) == n && "Hermitian matrix must be square");

    arb_t tol, mag;
    acb_t diff, cji;
    arb_init(tol);
    arb_init(mag);
    acb_init(diff);
    acb_init(cji);
    aic_mat_int_tol(tol, prec);

    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            acb_conj(cji, acb_mat_entry(H, j, i));
            acb_sub(diff, acb_mat_entry(H, i, j), cji, prec);
            acb_abs(mag, diff, prec);
            if (!arb_le(mag, tol)) {
                fprintf(stderr,
                        "aic_mat: input not Hermitian at (%ld,%ld): "
                        "|H_ij - conj(H_ji)| exceeds tol\n",
                        (long) i, (long) j);
                abort();
            }
        }
    }

    acb_clear(cji);
    acb_clear(diff);
    arb_clear(mag);
    arb_clear(tol);
}

void aic_mat_int_eig_certified(acb_ptr E, acb_mat_t R, const acb_mat_t A,
                               slong prec)
{
    slong n = acb_mat_nrows(A);
    acb_ptr E_approx = _acb_vec_init(n);
    acb_mat_t R_approx;
    acb_mat_init(R_approx, n, n);

    /* Stage 1: heuristic approximations (the return value is only a convergence
     * heuristic; the certified stage below is what we trust). */
    acb_mat_approx_eig_qr(E_approx, NULL, R_approx, A, NULL, 0, prec);

    /* Stage 2: prove isolating enclosures. R columns = right eigenvectors. */
    int ok = acb_mat_eig_simple(E, NULL, R, A, E_approx, R_approx, prec);
    if (!ok) {
        fprintf(stderr,
                "aic_mat_int_eig_certified: acb_mat_eig_simple could not "
                "isolate the spectrum at prec=%ld (clustered eigenvalues?); "
                "raise prec or audition acb_mat_eig_multiple\n",
                (long) prec);
        abort();
    }

    acb_mat_clear(R_approx);
    _acb_vec_clear(E_approx, n);
}

void aic_mat_herm_max_eig(arb_t out, const acb_mat_t H, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(n >= 1);
    aic_mat_int_assert_hermitian(H, prec);

    /* Heuristic approximate eigenpairs to seed the rigorous global enclosure. */
    acb_ptr E = _acb_vec_init(n);
    acb_mat_t R;
    acb_mat_init(R, n, n);
    acb_mat_approx_eig_qr(E, NULL, R, H, NULL, 0, prec);

    /* eps: every eigenvalue of H is within eps of some E[k] (acb_mat.rst:653).
     * Valid for repeated eigenvalues (no isolation required). */
    mag_t eps;
    mag_init(eps);
    acb_mat_eig_global_enclosure(eps, H, E, R, prec);

    /* max_k Re(E_k): the largest approximate eigenvalue. */
    arb_t maxre, re;
    arb_init(maxre);
    arb_init(re);
    acb_get_real(maxre, E + 0);
    for (slong k = 1; k < n; k++) {
        acb_get_real(re, E + k);
        arb_max(maxre, maxre, re, prec);
    }

    /* lambda_max in [maxre - eps, maxre + eps]: inflate the midpoint by eps. */
    arb_add_error_mag(maxre, eps);
    arb_set(out, maxre);

    arb_clear(re);
    arb_clear(maxre);
    mag_clear(eps);
    acb_mat_clear(R);
    _acb_vec_clear(E, n);
}
