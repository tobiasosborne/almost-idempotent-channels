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
 *
 * NEAR-ZERO ROBUSTNESS (bead aic-92f finding). For a near-zero Hermitian H —
 * e.g. the Gram matrix of S^2 - S for an EXACTLY idempotent superoperator, whose
 * entries are ~1e-31 (the eta=0 oracle of assoc_ecsa) — acb_mat_eig_global_-
 * enclosure can return a NON-FINITE radius eps (mag overflow on the tiny,
 * effectively-degenerate spectrum). A naive [maxre +/- inf] then sqrt's to
 * nan+/-inf in aic_mat_opnorm, spuriously tripping prop_P's basin guard on the
 * cleanest ground-truth input. Root-cause fix (Rule 3, not a bandaid): when eps
 * is non-finite, fall back to the RIGOROUS, eig-free bound — every eigenvalue of
 * Hermitian H satisfies |lambda| <= ||H||_op <= ||H||_F — and return the sound
 * enclosure [-||H||_F, ||H||_F] (midpoint 0, radius ||H||_F). This never NaNs,
 * is always valid (the heuristic QR midpoint is NOT trusted when the certifier
 * failed), and is tight exactly where it is used: a near-zero H gives a near-zero
 * ||H||_F, so opnorm returns ~0 with a small, finite radius.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/mag.h>

#include "aic/aic_mat.h"
#include "aic_mat_internal.h"

void aic_mat_int_tol(arb_t tol, slong prec)
{
    arb_set_si(tol, 1);
    arb_mul_2exp_si(tol, tol, -(prec - 8));
}

/* Non-aborting Hermiticity predicate. Returns 1 iff |H[i,j] - conj(H[j,i])| is
 * certainly within a RELATIVE + absolute-floor tol tol*(1 + |H_ij| + |H_ji|) for
 * all i,j, else 0 (bead aic-2yo).
 *
 * Why relative, not the bare absolute tol (the bug aic-2yo fixed). For a
 * faithfully-rounded Hermitian H the flagged "asymmetry" |H_ij - conj(H_ji)| is
 * the arb ball RADIUS, which grows with the entry MAGNITUDE. The canonical
 * trigger is a Gram G = D^dag D of a graded/ill-conditioned D: G_ii ~ r^{2i}
 * carries an imaginary-part radius ~r^{2i} 2^-prec (from the products), so the
 * asymmetry is ~r^{2i} 2^-prec while the old ABSOLUTE tol 2^-(prec-8)=256 2^-prec
 * is fixed -> false-fail once r^{2i} > 128, a PRECISION-INDEPENDENT abort (tol and
 * radius both scale as 2^-prec). aic_mat_singular_values / aic_mat_opnorm (which
 * route through the Gram) then aborted on any condition number >~ 1e2. Scaling the
 * tol by (1 + |H_ij| + |H_ji|) tracks that magnitude: the +1 keeps the original
 * absolute floor for O(1)/zero entries, so a GENUINELY non-Hermitian input
 * (asymmetry ~ magnitude >> tol*magnitude) is still rejected (teeth preserved;
 * tested in test_mat.c). */
int aic_mat_int_is_hermitian(const acb_mat_t H, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(acb_mat_ncols(H) == n && "Hermitian matrix must be square");

    arb_t tol, mag, scale, hij, hji;
    acb_t diff, cji;
    arb_init(tol);
    arb_init(mag);
    arb_init(scale);
    arb_init(hij);
    arb_init(hji);
    acb_init(diff);
    acb_init(cji);
    aic_mat_int_tol(tol, prec);

    int herm = 1;
    for (slong i = 0; i < n && herm; i++) {
        for (slong j = 0; j < n; j++) {
            acb_conj(cji, acb_mat_entry(H, j, i));
            acb_sub(diff, acb_mat_entry(H, i, j), cji, prec);
            acb_abs(mag, diff, prec);                 /* |H_ij - conj(H_ji)| */
            acb_abs(hij, acb_mat_entry(H, i, j), prec);
            acb_abs(hji, acb_mat_entry(H, j, i), prec);
            arb_add(scale, hij, hji, prec);
            arb_add_si(scale, scale, 1, prec);        /* 1 + |H_ij| + |H_ji| */
            arb_mul(scale, scale, tol, prec);         /* effective tol */
            if (!arb_le(mag, scale)) {                /* fail iff certainly > tol */
                herm = 0;
                break;
            }
        }
    }

    acb_clear(cji);
    acb_clear(diff);
    arb_clear(hji);
    arb_clear(hij);
    arb_clear(scale);
    arb_clear(mag);
    arb_clear(tol);
    return herm;
}

void aic_mat_int_assert_hermitian(const acb_mat_t H, slong prec)
{
    if (!aic_mat_int_is_hermitian(H, prec)) {
        fprintf(stderr,
                "aic_mat: input not Hermitian: |H_ij - conj(H_ji)| exceeds the "
                "relative tol tol*(1 + |H_ij| + |H_ji|) (bead aic-2yo)\n");
        abort();
    }
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

    arb_t maxre, re;
    arb_init(maxre);
    arb_init(re);

    if (mag_is_finite(eps)) {
        /* max_k Re(E_k): the largest approximate eigenvalue. lambda_max in
         * [maxre - eps, maxre + eps]: inflate the midpoint by eps. */
        acb_get_real(maxre, E + 0);
        for (slong k = 1; k < n; k++) {
            acb_get_real(re, E + k);
            arb_max(maxre, maxre, re, prec);
        }
        arb_add_error_mag(maxre, eps);
        arb_set(out, maxre);
    } else {
        /* Certifier failed (near-zero / degenerate H): rigorous eig-free bound
         * |lambda| <= ||H||_F => lambda_max in [-||H||_F, ||H||_F]. */
        arb_t fro;
        arb_init(fro);
        aic_mat_frobenius_norm(fro, H, prec);
        arb_zero(out);
        arb_add_error(out, fro);   /* [-||H||_F, ||H||_F] */
        arb_clear(fro);
    }

    arb_clear(re);
    arb_clear(maxre);
    mag_clear(eps);
    acb_mat_clear(R);
    _acb_vec_clear(E, n);
}
