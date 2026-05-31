/* aic_mat_norms.c — norms and singular values over acb_mat_t (beads aic-9zh,
 * "T-mat"). See include/aic_mat.h for the API contract.
 *
 * Frobenius norm is a thin wrapper over FLINT's certified acb_mat_frobenius_norm
 * (acb_mat.h:198) — provided so callers depend on the aic_mat surface, not the
 * FLINT spelling.
 *
 * Operator norm and singular values both work through the Hermitian PSD Gram
 * matrix G = A^dag A (acb_mat_conjugate_transpose, acb_mat.h:188, then
 * acb_mat_mul, acb_mat.h:213). The eigenvalues of G are the squared singular
 * values of A. The two routines differ in WHICH eig route they take, because
 * singular values need
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
 *
 * THE GRAM HERMITIANIZATION (bead aic-qgs; FINDINGS §C5/§C10). G = A^dag A is
 * Hermitian BY CONSTRUCTION, but acb_mat_mul accumulates an INDEPENDENT arb
 * radius into each entry, so entry (i,j) and conj of (j,i) end up with ~equal
 * midpoints but unrelated radii. The downstream certified-eig path runs the
 * non-aborting Hermiticity predicate aic_mat_int_is_hermitian (via herm_max_eig
 * / eig_hermitian), which tests the radius of the DIFFERENCE ball
 * |G_ij - conj(G_ji)| against tol*(1 + |G_ij| + |G_ji|). When A is built by a
 * deep matmul chain that CANCELS (e.g. the corner machinery's (2M-1)^2 - I on an
 * oblique S_P wrapper: ||(2M-1)^2-I|| ~ 1e-3 but the entry radii are ~1e-70 from
 * the oblique star), the Gram entries are SMALL in magnitude (~1e-6) yet carry a
 * LARGE accumulated radius (~1e-72), so the magnitude-relative floor (~tol, since
 * |G| << 1) is ~1e3x too small and the predicate FALSE-FAILS -> SIGABRT, on a
 * matrix that is genuinely Hermitian. (Measured: max rad(G)/|G| ~ 3e-67 vs
 * 2^-prec ~ 9e-78 at prec=256 -- ~37 bits lost through the chain.)
 *
 * Symmetrizing G <- (G + G^dag)/2 alone does NOT cure this (verified): the
 * off-diagonal midpoints of A^dag A are ALREADY exact conjugates, so the residual
 * the predicate flags is purely the difference-ball RADIUS, which no midpoint
 * manipulation removes (arb subtraction adds the radii). The rigorous cure is to
 * split the certified Gram into (i) an EXACTLY-Hermitian midpoint matrix Gmid
 * (zero radius, so the predicate passes for the right reason) and (ii) a rigorous
 * operator-norm perturbation bound R = ||G_true - Gmid||_op <= ||G_true - Gmid||_F
 * on the dropped radii. By Weyl's inequality for Hermitian matrices, every
 * eigenvalue of the true Gram lies within R of the corresponding eigenvalue of
 * Gmid, so lambda_max(G_true) in [lambda_max(Gmid) - R, lambda_max(Gmid) + R].
 * aic_mat_gram returns Gmid + R; the callers run the certified eig on Gmid and
 * inflate by R before the sqrt. This is the substrate version of the
 * aic_corner_gamma_opnorm_ub mid+radius workaround (now retireable, FINDINGS
 * §C5). It changes NO value on tight inputs (R ~ 2^-prec there, below the
 * existing global-enclosure radius), and does NOT touch the herm_max_eig /
 * is_hermitian guard -- intact for its DIRECT callers (CP-cert herm_max_eig(-C)).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"

void aic_mat_frobenius_norm(arb_t out, const acb_mat_t A, slong prec)
{
    acb_mat_frobenius_norm(out, A, prec);
}

/* Build the EXACTLY-HERMITIAN midpoint Gram of A into Gmid (caller-cleared/
 * uninit'd: this inits it), set *R to a rigorous operator-norm perturbation
 * bound R >= ||G_true - Gmid||_op, and return the dimension. n<=m: the true Gram
 * is A^dag A (n x n); else A A^dag (m x m). See the aic-qgs paragraph in the file
 * docstring for why we do NOT hand the raw acb_mat_mul output to the certified
 * eig path. Steps:
 *   1. G = A^dag A (or A A^dag), then symmetrize G <- (G + G^dag)/2. G_true is
 *      Hermitian, so this is an exact identity on the true value and a valid
 *      enclosure; it equalises the per-entry radii (cosmetic) but the antiherm
 *      MIDPOINT is already 0.
 *   2. Gmid = the EXACT-Hermitian matrix of midpoints of G: diagonal imaginary
 *      part set to 0, lower triangle set to the exact conjugate of the upper.
 *      Gmid has zero radius, so aic_mat_int_is_hermitian passes it for the right
 *      reason (no false-fail on accumulated radius). Gmid IS Hermitian, so its
 *      eigenvalues are real and herm_max_eig / eig_hermitian apply.
 *   3. R = ||G - Gmid||_F (the certified UPPER bound of that Frobenius ball).
 *      Since G_true in G (the ball matrix) and |G_true[i,j] - Gmid[i,j]| <=
 *      rad(G[i,j]), R >= ||G_true - Gmid||_F >= ||G_true - Gmid||_op. By Weyl's
 *      inequality every eigenvalue of G_true is within R of the matching
 *      eigenvalue of Gmid; the callers inflate by R before the sqrt. */
static slong aic_mat_gram(acb_mat_t Gmid, arb_t R, const acb_mat_t A, slong prec)
{
    slong m = acb_mat_nrows(A);
    slong n = acb_mat_ncols(A);
    acb_mat_t Ah, G;
    acb_mat_init(Ah, n, m);
    acb_mat_conjugate_transpose(Ah, A);

    slong g = (n <= m) ? n : m;
    acb_mat_init(G, g, g);
    if (n <= m)
        acb_mat_mul(G, Ah, A, prec); /* (n x m)(m x n) = n x n */
    else
        acb_mat_mul(G, A, Ah, prec); /* (m x n)(n x m) = m x m */
    acb_mat_clear(Ah);

    /* G <- (G + G^dag)/2: equalise the per-entry radii (G_true is Hermitian, so
     * this is exact on the true value). The antiherm midpoint is already 0. */
    acb_mat_t Gt;
    acb_mat_init(Gt, g, g);
    acb_mat_conjugate_transpose(Gt, G);
    acb_mat_add(G, G, Gt, prec);
    acb_mat_scalar_mul_2exp_si(G, G, -1);
    acb_mat_clear(Gt);

    /* Gmid: exact-Hermitian matrix of midpoints (zero radius). Diagonal real,
     * lower = exact conjugate of upper. This is the matrix the certified eig
     * runs on; it passes the Hermiticity predicate by construction. */
    acb_mat_init(Gmid, g, g);
    for (slong i = 0; i < g; i++) {
        acb_get_mid(acb_mat_entry(Gmid, i, i), acb_mat_entry(G, i, i));
        arb_zero(acb_imagref(acb_mat_entry(Gmid, i, i))); /* diag is real */
        for (slong j = i + 1; j < g; j++) {
            acb_get_mid(acb_mat_entry(Gmid, i, j), acb_mat_entry(G, i, j));
            acb_conj(acb_mat_entry(Gmid, j, i), acb_mat_entry(Gmid, i, j));
        }
    }

    /* R = ||G - Gmid||_F (upper bound): the certified perturbation radius. */
    acb_mat_t Diff;
    acb_mat_init(Diff, g, g);
    acb_mat_sub(Diff, G, Gmid, prec);
    arb_t fro;
    arb_init(fro);
    aic_mat_frobenius_norm(fro, Diff, prec);
    arb_get_ubound_arf(arb_midref(R), fro, prec); /* R = upper bound, zero radius */
    mag_zero(arb_radref(R));
    arb_clear(fro);
    acb_mat_clear(Diff);
    acb_mat_clear(G);

    return g;
}

/* lambda <- [lambda - R, lambda + R] (Weyl inflation by the Gram perturbation
 * bound R), clamped at 0 from below since the true Gram is PSD (any negative
 * excursion is the inflation overshooting the PSD floor, not signal). */
static void inflate_by_R(arb_t lambda, const arb_t R, slong prec)
{
    arb_add_error(lambda, R);

    arf_t lo, ub;
    arf_init(lo);
    arf_init(ub);
    arb_get_lbound_arf(lo, lambda, prec);
    if (arf_sgn(lo) < 0) {
        arb_get_ubound_arf(ub, lambda, prec);
        arf_zero(lo);                              /* clamp lower bound to 0 */
        arb_set_interval_arf(lambda, lo, ub, prec); /* sound enclosure [0, ub] */
    }
    arf_clear(ub);
    arf_clear(lo);
}

void aic_mat_opnorm(arb_t out, const acb_mat_t A, slong prec)
{
    acb_mat_t Gmid;
    arb_t R;
    arb_init(R);
    (void) aic_mat_gram(Gmid, R, A, prec);

    /* largest eigenvalue of the exact-Hermitian midpoint Gram, then Weyl-inflate
     * by R to enclose the true Gram's lambda_max (= largest singular value^2).
     * Use the degeneracy-robust max-eig (repeated singular values, e.g. the
     * identity, must work) rather than the simple-spectrum eig path. */
    arb_t lmax;
    arb_init(lmax);
    aic_mat_herm_max_eig(lmax, Gmid, prec);
    inflate_by_R(lmax, R, prec);
    arb_sqrtpos(out, lmax, prec);

    arb_clear(lmax);
    arb_clear(R);
    acb_mat_clear(Gmid);
}

void aic_mat_singular_values(arb_ptr svals, const acb_mat_t A, slong prec)
{
    slong m = acb_mat_nrows(A);
    slong n = acb_mat_ncols(A);
    slong k = (m < n) ? m : n; /* number of singular values */

    acb_mat_t Gmid;
    arb_t R;
    arb_init(R);
    slong g = aic_mat_gram(Gmid, R, A, prec);
    assert(g == k);

    arb_ptr ev = _arb_vec_init(g);
    aic_mat_eig_hermitian(ev, NULL, Gmid, prec); /* ascending eigenvalues of Gmid */

    /* singular values DESCENDING: sqrt of (eigenvalue Weyl-inflated by R),
     * reversed. Every eigenvalue of the true Gram is within R of Gmid's. */
    for (slong i = 0; i < g; i++) {
        arb_t lam;
        arb_init(lam);
        arb_set(lam, ev + (g - 1 - i));
        inflate_by_R(lam, R, prec);
        arb_sqrtpos(svals + i, lam, prec);
        arb_clear(lam);
    }

    _arb_vec_clear(ev, g);
    arb_clear(R);
    acb_mat_clear(Gmid);
}
