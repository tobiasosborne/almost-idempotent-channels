/* aic_mat_eig.c — full Hermitian eigendecomposition over acb_mat_t (beads
 * aic-9zh, "T-mat"). Eigenvalues AND eigenvectors, via the shared certified
 * two-stage path (aic_mat_int_eig_certified, see aic_mat_spectral.c / the
 * internal header). See include/aic_mat.h for the API contract.
 *
 * Steps:
 *   1. assert the input is Hermitian (fail loud, CLAUDE.md Rule 4);
 *   2. certified general eig -> n isolating eigenvalue enclosures E with
 *      right-eigenvector columns R such that L A R = diag(E), L = R^{-1};
 *   3. EXPLOIT Hermiticity to check correctness: a Hermitian matrix has a real
 *      spectrum, so each CERTIFIED imaginary-part enclosure must CONTAIN 0 (the
 *      rigorous "this eigenvalue is real"). A ball excluding 0 means the solver
 *      certified a non-real eigenvalue, contradicting the asserted Hermiticity —
 *      a real bug, so abort;
 *   4. the eigenvalues come out "in no particular order"; downstream code
 *      (operator norm, projection spectral split) needs them ordered, so sort
 *      ascending by real-part MIDPOINT and carry the eigenvector columns in
 *      lockstep. Sorting by midpoint does not weaken the certified enclosures.
 *
 * Requires a SIMPLE spectrum (acb_mat_eig_simple proves ISOLATED eigenvalues);
 * a repeated eigenvalue aborts in step 2. Degenerate-spectrum eigenvectors (via
 * acb_mat_eig_multiple's invariant-subspace enclosures) are a later audition.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_mat.h"
#include "aic_mat_internal.h"

void aic_mat_eig_hermitian(arb_ptr evals, acb_mat_t evecs,
                           const acb_mat_t H, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(n >= 1);
    aic_mat_int_assert_hermitian(H, prec);

    acb_ptr E = _acb_vec_init(n);
    acb_mat_t R;
    acb_mat_init(R, n, n);

    aic_mat_int_eig_certified(E, R, H, prec);

    /* Hermitian => real spectrum: each certified imaginary-part enclosure must
     * CONTAIN 0. A ball excluding 0 contradicts the asserted Hermiticity. This
     * is the correct test, not an absolute |imag| <= tol: the eig solver's
     * enclosure radius is wider than 2^-prec, so a tol comparison false-fails. */
    arb_t im;
    arb_init(im);
    for (slong k = 0; k < n; k++) {
        acb_get_imag(im, E + k);
        if (!arb_contains_zero(im)) {
            fprintf(stderr,
                    "aic_mat_eig_hermitian: eigenvalue %ld certified non-real "
                    "(imag enclosure excludes 0) on a Hermitian input\n",
                    (long) k);
            abort();
        }
    }

    /* Insertion sort ascending by real-part midpoint, permuting evec columns of
     * R in lockstep. n is small (matrix dims here), so O(n^2) is fine. */
    slong *perm = malloc((size_t) n * sizeof(slong));
    assert(perm != NULL);
    for (slong k = 0; k < n; k++) perm[k] = k;

    arb_t rk;
    arb_init(rk);
    for (slong a2 = 1; a2 < n; a2++) {
        slong p = perm[a2];
        arb_t key;
        arb_init(key);
        acb_get_real(key, E + p);
        slong b2 = a2 - 1;
        while (b2 >= 0) {
            acb_get_real(rk, E + perm[b2]);
            if (arf_cmp(arb_midref(rk), arb_midref(key)) <= 0) break;
            perm[b2 + 1] = perm[b2];
            b2--;
        }
        perm[b2 + 1] = p;
        arb_clear(key);
    }
    arb_clear(rk);

    /* Emit eigenvalues (real parts) and, if requested, eigenvector columns in
     * the sorted order. */
    for (slong k = 0; k < n; k++)
        acb_get_real(evals + k, E + perm[k]);

    if (evecs != NULL) {
        assert(acb_mat_nrows(evecs) == n && acb_mat_ncols(evecs) == n);
        for (slong col = 0; col < n; col++)
            for (slong row = 0; row < n; row++)
                acb_set(acb_mat_entry(evecs, row, col),
                        acb_mat_entry(R, row, perm[col]));
    }

    free(perm);
    arb_clear(im);
    acb_mat_clear(R);
    _acb_vec_clear(E, n);
}
