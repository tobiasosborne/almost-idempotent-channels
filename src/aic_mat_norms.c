/* aic_mat_norms.c — norms and singular values over acb_mat_t (beads aic-9zh,
 * "T-mat"). See include/aic_mat.h for the API contract.
 *
 * Frobenius norm is a thin wrapper over FLINT's certified acb_mat_frobenius_norm
 * (acb_mat.h:198) — provided so callers depend on the aic_mat surface, not the
 * FLINT spelling.
 *
 * Operator norm and singular values both work through the Hermitian PSD Gram
 * matrix G = A^dag A (acb_mat_conjugate_transpose, acb_mat.h:188, then
 * acb_mat_mul, acb_mat.h:213), whose eigenvalues are the squared singular values
 * of A. They differ in WHICH eig route they take, because singular values need
 * per-value enclosures while the operator norm needs only the largest:
 *   - operator norm ||A||_op (the cb-norm ampliation building block,
 *     approximate_algebras.tex:349) = sqrt(largest eigenvalue of G), via the
 *     degeneracy-ROBUST aic_mat_herm_max_eig (global enclosure). This must work
 *     on matrices with repeated singular values (e.g. the identity), so it must
 *     NOT use the simple-spectrum eig path.
 *   - singular values = sqrt of ALL eigenvalues of G, via aic_mat_eig_hermitian
 *     (the simple/isolated-eigenvalue path); it therefore REQUIRES distinct
 *     singular values (degenerate case is a later audition bead).
 * arb_sqrtpos clamps the ball to [0,inf), correct here because G is PSD so the
 * true eigenvalues are >= 0 and any tiny negative excursion of an enclosure is a
 * rounding artefact, not signal.
 *
 * For a rectangular m x n matrix we form the SMALLER Gram matrix — A^dag A
 * (n x n) when n <= m, else A A^dag (m x m) — since the two share the same
 * nonzero spectrum and min(m,n) determines the count of singular values. This is
 * a cost optimisation, not a correctness requirement.
 *
 * Alternative routes NOT implemented here (Law 4 audition slate, to be filed as
 * beads): (a) power iteration on A^dag A for the largest eigenvalue (cheaper than
 * a full eig, but the certified enclosure needs a Rayleigh-quotient residual
 * bound); (b) a full SVD producing U, Sigma, V (FLINT has no SVD; would need a
 * Golub-Kahan bidiagonalisation + implicit-QR, or a one-sided Jacobi sweep);
 * (c) degenerate-spectrum singular values via acb_mat_eig_multiple's cluster
 * enclosures. We implement the Gram + (global-enclosure | simple-eig) routes.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_mat.h"

void aic_mat_frobenius_norm(arb_t out, const acb_mat_t A, slong prec)
{
    acb_mat_frobenius_norm(out, A, prec);
}

/* Build the smaller Gram matrix of A into G (caller-cleared/uninit'd: this
 * inits it) and return its dimension. n<=m: G = A^dag A (n x n); else
 * G = A A^dag (m x m). G is Hermitian PSD by construction. */
static slong aic_mat_gram(acb_mat_t G, const acb_mat_t A, slong prec)
{
    slong m = acb_mat_nrows(A);
    slong n = acb_mat_ncols(A);
    acb_mat_t Ah;
    acb_mat_init(Ah, n, m);
    acb_mat_conjugate_transpose(Ah, A);

    if (n <= m) {
        acb_mat_init(G, n, n);
        acb_mat_mul(G, Ah, A, prec); /* (n x m)(m x n) = n x n */
    } else {
        acb_mat_init(G, m, m);
        acb_mat_mul(G, A, Ah, prec); /* (m x n)(n x m) = m x m */
    }

    acb_mat_clear(Ah);
    return acb_mat_nrows(G);
}

void aic_mat_opnorm(arb_t out, const acb_mat_t A, slong prec)
{
    acb_mat_t G;
    (void) aic_mat_gram(G, A, prec);

    /* largest eigenvalue of the Hermitian PSD Gram = (largest singular value)^2.
     * Use the degeneracy-robust max-eig (repeated singular values, e.g. the
     * identity, must work) rather than the simple-spectrum eig path. */
    arb_t lmax;
    arb_init(lmax);
    aic_mat_herm_max_eig(lmax, G, prec);
    arb_sqrtpos(out, lmax, prec);

    arb_clear(lmax);
    acb_mat_clear(G);
}

void aic_mat_singular_values(arb_ptr svals, const acb_mat_t A, slong prec)
{
    slong m = acb_mat_nrows(A);
    slong n = acb_mat_ncols(A);
    slong k = (m < n) ? m : n; /* number of singular values */

    acb_mat_t G;
    slong g = aic_mat_gram(G, A, prec);
    assert(g == k);

    arb_ptr ev = _arb_vec_init(g);
    aic_mat_eig_hermitian(ev, NULL, G, prec); /* ascending eigenvalues of G */

    /* singular values DESCENDING: sqrt of eigenvalues, reversed. */
    for (slong i = 0; i < g; i++)
        arb_sqrtpos(svals + i, ev + (g - 1 - i), prec);

    _arb_vec_clear(ev, g);
    acb_mat_clear(G);
}
