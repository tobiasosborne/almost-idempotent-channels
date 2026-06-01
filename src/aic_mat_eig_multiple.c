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
 * THE AUDITION (CLAUDE.md Law 4), FLINT 3.x constraints. Route 1 (CHOSEN):
 * acb_mat_approx_eig_qr seed -> acb_mat_eig_multiple, which certifies ALL n
 * eigenvalues WITH multiplicity (simple van der Hoeven-Mourrain first, Rump
 * cluster enclosures as fallback) — degeneracy-native, all primitives present in
 * FLINT 3.0.1. Output balls are either disjoint or IDENTICAL, identical ones
 * grouped consecutively; a run of k identical balls is a certified k-cluster
 * isolated from the other n-k eigenvalues. (Route 2, eig-free Cholesky: DEAD END
 * — no acb_mat_cho, arb_mat_cho/ldl need strict PD, fail on semidefinite Choi.
 * Route 3, inertia count: needs the same degeneracy-robust primitive, no better.
 * FINDINGS §D5.)
 *
 * HERMITIAN => REAL SPECTRUM. The inputs are Hermitian, so the certified
 * imaginary part of each ball must enclose 0; we assert that (a ball excluding 0
 * means the solver certified a non-real eigenvalue, contradicting the asserted
 * Hermiticity — a real bug, abort). The seed (approx_eig_qr) return value is
 * only a convergence heuristic and is NOT trusted; eig_multiple's return value
 * is.
 *
 * THE §D5n WALL, RESOLVED BY DENSIFICATION (FINDINGS §D5n RESOLVED; design
 * docs/research/eigvec_certified_design.md §2 has the full root-cause + numbers).
 * acb_mat_eig_multiple returns 0 on two-clusters-each-multiplicity-≥2 Hermitian
 * inputs (a C^5 {2,3} projector, diag(2,2,0,0), ...) even at prec=256/1024. The
 * cause is NOT seed conditioning (the original §D5n hypothesis, proven FALSE) and
 * NOT precision: it is FLINT Rump's frozen-row partition (partition_X_sorted)
 * degenerating on a ROW-SPARSE invariant subspace. THE FIX: on eig_multiple(H)==0
 * RETRY on the densified A' = U H U† (U = aic_mat_dense_unitary, §1.3, spreading
 * every eigenvector across all n rows). The spectrum is conjugation-invariant
 * (U U† = I to ~1e-37 ⇒ spec(A')=spec(H)), so the balls of A' ARE those of H and
 * we write them straight to E. Measured: rescues diag(2,2,0,0)/(5,5,2,2)/(1,1,1,0),
 * C^5 {2,3}, all previously-failing carriers (design §2). Only if the DENSIFIED
 * retry ALSO returns 0 — a genuine near-degeneracy unresolvable at this prec, e.g.
 * two mult-2 clusters separated by 2^-10 at prec=24 — do we FAIL LOUD (Rule 4).
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
        /* DENSIFY-RETRY (FINDINGS §D5n RESOLVED; design §2). The 0 is FLINT
         * Rump's frozen-row partition failing on a ROW-SPARSE invariant subspace
         * (NOT seed conditioning, NOT precision). Conjugate by the dense unitary
         * U: A' = U H U†. The spectrum is conjugation-invariant (U U† = I to
         * ~1e-37), so spec(A') = spec(H) and the eigenvalue balls of A' ARE
         * those of H — we write them straight to E. Densification spreads every
         * eigenvector across all n rows, so Rump's partition is well-conditioned
         * on A' even when it was singular on H. */
        acb_mat_t U, Ud, A1, T;
        acb_mat_init(U, n, n);
        acb_mat_init(Ud, n, n);
        acb_mat_init(A1, n, n);
        acb_mat_init(T, n, n);
        aic_mat_dense_unitary(U, n, prec);
        acb_mat_conjugate_transpose(Ud, U);

        /* Assert U is unitary to far below the working prec before trusting the
         * conjugation: ||U U† - I||_F must be certified tiny (Rule 4). A loose U
         * would silently change the spectrum. */
        acb_mat_mul(T, U, Ud, prec);
        for (slong i = 0; i < n; i++)
            acb_sub_si(acb_mat_entry(T, i, i), acb_mat_entry(T, i, i), 1, prec);
        arb_t uu;
        arb_init(uu);
        acb_mat_frobenius_norm(uu, T, prec);
        arb_t utol;
        arb_init(utol);
        /* Tolerance n^2 * 2^-(prec-8): the densifier chains n(n-1)/2 Givens, so
         * the certified ||U U†-I||_F radius accumulates ~ n^2 * 2^-prec and the
         * bare floor is exceeded for n >= 6 (a latent fail-loud on legitimate
         * inputs). FINDINGS §D5n2; mirrors aic_mat_eigvec_seed.c. */
        aic_mat_int_tol(utol, prec);
        arb_mul_si(utol, utol, n * n, prec);
        if (!arb_lt(uu, utol)) {
            char *uu_s = arb_get_str(uu, 8, 0);
            fprintf(stderr,
                    "aic_mat_eig_hermitian_multiple: densifier U not certified "
                    "unitary (||U U†-I||_F=%s not < n^2*2^-(prec-8)) at prec=%ld "
                    "— the conjugation would not preserve the spectrum; bead "
                    "aic-4td\n",
                    uu_s ? uu_s : "?", (long) prec);
            flint_free(uu_s);
            abort();
        }
        arb_clear(utol);
        arb_clear(uu);

        acb_mat_mul(T, U, H, prec);
        acb_mat_mul(A1, T, Ud, prec);       /* A' = U H U† */

        /* Re-seed on A' (the seed eigenvectors differ; the eigenvalues match). */
        acb_mat_approx_eig_qr(E_approx, NULL, R_approx, A1, NULL, 0, prec);
        int ok2 = acb_mat_eig_multiple(E, A1, E_approx, R_approx, prec);

        acb_mat_clear(T);
        acb_mat_clear(A1);
        acb_mat_clear(Ud);
        acb_mat_clear(U);

        if (!ok2) {
            /* Even the densified A' is unresolved — a GENUINE near-degeneracy at
             * this prec (two mult-≥2 clusters separated by far below 2^-prec;
             * measured: {2,2} clusters at gap 2^-10, prec=24). Raising prec
             * resolves a true gap. Fail loud (Rule 4). */
            fprintf(stderr,
                    "aic_mat_eig_hermitian_multiple: acb_mat_eig_multiple could "
                    "not certify the spectrum at prec=%ld even after dense-unitary "
                    "densification (clusters unresolved — a genuine near-degeneracy "
                    "at this prec; raise prec); FINDINGS §D5n, bead aic-4td\n",
                    (long) prec);
            abort();
        }
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
