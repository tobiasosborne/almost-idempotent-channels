/* aic_mat_eig_multiple.c — certified eigenvalues + numerical rank for a
 * HERMITIAN matrix with a possibly DEGENERATE spectrum (bead aic-w4o.1, the
 * recurring "certified degenerate-eigenvalue wall", paper/FINDINGS.md §D5).
 *
 * WHY THIS EXISTS. FLINT 3.0.1's certified simple-spectrum path
 * (acb_mat_eig_simple, wrapped by aic_mat_eig_hermitian / aic_mat_spectral.c)
 * REQUIRES distinct eigenvalues and returns 0 — so we abort (Rule 4) — on any
 * repeat. But this project is dominated by DEGENERATE spectra: projections have
 * eigenvalues 0/1 with high multiplicity, and the target algebras are direct
 * sums of matrix blocks (+) B(L_j) carrying multiplicities. Certified rank /
 * subspace extraction has therefore been on the LAPACK double path with the
 * certification deferred here (FINDINGS §D5). This file is the eigenvalue+rank
 * layer of that deferred work.
 *
 * THE AUDITION (CLAUDE.md Law 4), FLINT 3.x constraints.
 *  - Route 1 (CHOSEN): acb_mat_approx_eig_qr seed -> acb_mat_eig_multiple. The
 *    latter certifies ALL n eigenvalues WITH multiplicity: it tries the simple
 *    (van der Hoeven-Mourrain) method first and falls back to Rump cluster
 *    enclosures, so it is degeneracy-native. Output balls are either disjoint or
 *    IDENTICAL, identical ones grouped consecutively; a run of k identical balls
 *    is a certified k-cluster that could not be split at this prec but IS
 *    isolated from the other n-k eigenvalues (acb_mat.rst). Every primitive is
 *    present in FLINT 3.0.1.
 *  - Route 2 (DEAD END in FLINT 3.x): an eig-free pivoted Cholesky. There is no
 *    acb_mat_cho; arb_mat_cho / arb_mat_ldl require STRICT positive-definiteness
 *    and return 0 on the semidefinite Choi/carrier we must handle.
 *  - Route 3 (inertia count): viable but needs the same degeneracy-robust
 *    primitive; not better than Route 1, which gives the eigenvalues outright.
 * Route 1 is the only one whose primitives all exist and which is
 * degeneracy-native, so it is the implemented route for this increment.
 *
 * HERMITIAN => REAL SPECTRUM. The inputs are Hermitian, so the certified
 * imaginary part of each ball must enclose 0; we assert that (a ball excluding 0
 * means the solver certified a non-real eigenvalue, contradicting the asserted
 * Hermiticity — a real bug, abort). The seed (approx_eig_qr) return value is
 * only a convergence heuristic and is NOT trusted; eig_multiple's return value
 * is. If eig_multiple returns 0 the clusters are unresolved at this prec and we
 * FAIL LOUD (Rule 4) with a "raise prec" message.
 *
 * CERTIFIED RANK. aic_mat_certified_rank counts eigenvalues certified ABOVE a
 * threshold ball thr. Each ball must be certified-above (arb_gt) OR
 * certified-below (arb_lt) thr; a ball that STRADDLES thr leaves the rank
 * unresolved at this prec and aborts (Rule 4). The fail-loud path uses an
 * fprintf+abort (NDEBUG-immune) so the tooth fires even when asserts are
 * stripped — matching the convention in aic_funcalc_sgn.c / aic_mat_spectral.c.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic_mat_internal.h"

void aic_mat_eig_hermitian_multiple(acb_ptr E, const acb_mat_t H, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(n >= 1 && acb_mat_ncols(H) == n &&
           "aic_mat_eig_hermitian_multiple: H must be square");
    aic_mat_int_assert_hermitian(H, prec);

    acb_ptr E_approx = _acb_vec_init(n);
    acb_mat_t R_approx;
    acb_mat_init(R_approx, n, n);

    /* Seed: heuristic float approximations (no error guarantee; the certified
     * stage below is what we trust). */
    acb_mat_approx_eig_qr(E_approx, NULL, R_approx, H, NULL, 0, prec);

    /* Certify ALL n eigenvalues with multiplicity (Rump cluster fallback). */
    int ok = acb_mat_eig_multiple(E, H, E_approx, R_approx, prec);
    if (!ok) {
        /* eig_multiple's Rump enclosure failed to certify. This is usually a
         * SEED-CONDITIONING limit, not a precision limit (FINDINGS §D5n): when
         * two clusters each have multiplicity >= 2 and the approx_eig_qr seed
         * gives near-parallel cluster eigenvectors, the invariant-subspace
         * enclosure cannot close — and raising prec does NOT help (measured:
         * still 0 at prec=256). Fail loud (Rule 4); raising prec is worth a try
         * only if the spectrum has a true near-degeneracy at this prec. */
        fprintf(stderr,
                "aic_mat_eig_hermitian_multiple: acb_mat_eig_multiple could not "
                "certify the spectrum at prec=%ld (clusters unresolved — likely "
                "a seed-conditioning limit, FINDINGS §D5n; raising prec may not "
                "help); bead aic-w4o.1\n",
                (long) prec);
        abort();
    }

    /* Hermitian => real spectrum: each certified imaginary part must enclose 0.
     * A ball excluding 0 contradicts the asserted Hermiticity (a real bug). */
    arb_t im;
    arb_init(im);
    for (slong k = 0; k < n; k++) {
        acb_get_imag(im, E + k);
        if (!arb_contains_zero(im)) {
            fprintf(stderr,
                    "aic_mat_eig_hermitian_multiple: eigenvalue %ld certified "
                    "non-real (imag enclosure excludes 0) on a Hermitian input\n",
                    (long) k);
            abort();
        }
    }
    arb_clear(im);

    acb_mat_clear(R_approx);
    _acb_vec_clear(E_approx, n);
}

slong aic_mat_certified_rank(const acb_mat_t H, const arb_t thr, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(n >= 1 && acb_mat_ncols(H) == n &&
           "aic_mat_certified_rank: H must be square");

    acb_ptr E = _acb_vec_init(n);
    aic_mat_eig_hermitian_multiple(E, H, prec);

    arb_t re;
    arb_init(re);
    slong rank = 0;
    for (slong k = 0; k < n; k++) {
        acb_get_real(re, E + k);
        if (arb_gt(re, thr)) {
            rank++;                       /* certified above thr */
        } else if (arb_lt(re, thr)) {
            /* certified below thr — not counted */
        } else {
            /* The eigenvalue ball overlaps thr: above/below is UNDECIDED, so the
             * rank is not certified at this prec. Fail loud (Rule 4), NDEBUG-
             * immune fprintf+abort (the assert is the documented invariant). */
            char *re_s = arb_get_str(re, 12, 0);
            char *thr_s = arb_get_str(thr, 12, 0);
            fprintf(stderr,
                    "aic_mat_certified_rank: eigenvalue ball %ld = %s STRADDLES "
                    "the threshold %s (rank unresolved at prec=%ld); raise prec "
                    "or move the threshold out of the cluster (bead aic-w4o.1)\n",
                    (long) k, re_s ? re_s : "?", thr_s ? thr_s : "?",
                    (long) prec);
            flint_free(re_s);
            flint_free(thr_s);
            assert(0 && "aic_mat_certified_rank: eigenvalue ball straddles thr");
            abort();
        }
    }
    arb_clear(re);
    _acb_vec_clear(E, n);
    return rank;
}
