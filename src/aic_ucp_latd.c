/* aic_ucp_latd.c — UCP-map extractions that need a DEGENERATE-spectrum
 * eigensolver, on the LAPACK double path (bead aic-c7n). The certified arb
 * eigensolver (acb_mat_eig_simple) requires a SIMPLE spectrum and aborts on
 * repeated eigenvalues (bead aic-w4o.1); the Choi matrix of a UCP map and the
 * carrier operator Q generically have degenerate spectra (a rank-r channel has a
 * Choi matrix of rank r << dim, i.e. a zero eigenvalue of high multiplicity).
 * LAPACKE_zheev diagonalises them happily, so the extraction lives here.
 *
 * Choi->Kraus (shard G:57; classical Choi [Cho75]). C_Phi >= 0 is PSD with a
 * spectral decomposition C_Phi = sum_a lambda_a v_a v_a^dag, i.e.
 *   C[i*h+a, j*h+b] = sum_a lambda_a v_a[i*h+a] conj(v_a[j*h+b]).
 * With the header's Convention A index map
 *   C[i*h+a, j*h+b] = sum_x conj(K_x[i,a]) K_x[j,b],
 * matching term-by-term forces conj(K_a[i,c]) = sqrt(lambda_a) v_a[i*h+c], so
 *   K_a[i, c] = sqrt(lambda_a) * conj(v_a[i*h + c])
 * (the CONJUGATE of the reshaped eigenvector). Then kraus_to_choi(extracted) ==
 * C_Phi exactly: the conjugated eigenvectors reshaped (dim_K x dim_H, row index
 * i*dim_H + c -> (i,c)) ARE the Kraus operators up to the sqrt(eigenvalue)
 * scale. We keep eigenpairs with
 *   lambda_a > thr = (dim_K*dim_H) * DBL_EPSILON * ||C||_F     (QETLAB-style)
 * dropping the numerically-zero tail. The recovered set is a valid Kraus rep of
 * the SAME channel, but only up to the unitary gauge freedom of Kraus
 * decompositions — callers must compare AS CHANNELS (rebuild Choi), never
 * operator-by-operator.
 *
 * Carrier rank (lem_carrier, .tex:1724): the carrier dim is the rank of
 * Q = sum_a K_a K_a^dag, i.e. the count of eigenvalues of Q above
 *   thr = dim_K * DBL_EPSILON * ||Q||_F.
 * This is the UNCERTIFIED fast rank; the certified (arb, gap-dependent) rank is
 * blocked on aic-w4o.1 and is a documented module TODO.
 *
 * This file calls LAPACK (via aic_latd) only — no acb eigensolver — per the
 * latd/mat number-path separation. It uses aic_latd's acb_mat<->array helpers
 * for I/O of the certified-path matrices.
 */
#include <assert.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_latd.h"
#include "aic_ucp.h"

/* Frobenius norm of a row-major n x n complex array (used for the threshold). */
static double frob_norm(const double _Complex *A, slong n)
{
    double s = 0.0;
    for (slong k = 0; k < n * n; k++) {
        double m = cabs(A[k]);
        s += m * m;
    }
    return sqrt(s);
}

void aic_ucp_choi_to_kraus_latd(aic_ucp_kraus *phi, const acb_mat_t C,
                                slong dim_K, slong dim_H)
{
    slong n = dim_K * dim_H;
    assert(acb_mat_nrows(C) == n && acb_mat_ncols(C) == n &&
           "choi_to_kraus_latd: C must be (dim_K*dim_H) square");

    double _Complex *Ca = malloc((size_t) (n * n) * sizeof(double _Complex));
    double _Complex *evecs = malloc((size_t) (n * n) * sizeof(double _Complex));
    double *evals = malloc((size_t) n * sizeof(double));
    assert(Ca && evecs && evals && "choi_to_kraus_latd: out of memory");

    aic_latd_from_acb_mat(Ca, C);
    double thr = (double) n * DBL_EPSILON * frob_norm(Ca, n);
    aic_latd_eig_hermitian(evals, evecs, Ca, n); /* ascending */

    /* count kept eigenpairs (lambda > thr); PSD => negatives are noise below thr.
     * A clearly-negative eigenvalue (< -thr) means C is not PSD: fail loud. */
    slong r = 0;
    for (slong a = 0; a < n; a++) {
        /* A clearly-negative eigenvalue (< -thr) means C is not PSD; sqrt'ing it
         * below would produce NaN Kraus ops. Fail loud in EVERY build (not just
         * with assertions enabled) — corrupted output is worse than a crash
         * (CLAUDE.md Rule 4), matching aic_ucp_is_cp_choi. */
        if (evals[a] <= -thr) {
            flint_printf("aic_ucp_choi_to_kraus_latd: eigenvalue %g < -thr=%g; "
                         "C is not PSD (not CP) — cannot extract Kraus ops "
                         "(Rule 4: fail loud)\n", evals[a], -thr);
            flint_abort();
        }
        if (evals[a] > thr) r++;
    }
    assert(r >= 1 && "choi_to_kraus_latd: zero rank (the zero map is not UCP)");

    aic_ucp_kraus_init(phi, dim_K, dim_H, r);
    /* fill from the kept (largest) eigenpairs; column a of evecs is v_a. */
    slong out_a = 0;
    for (slong a = 0; a < n; a++) {
        if (evals[a] <= thr) continue;
        double scale = sqrt(evals[a]);
        for (slong i = 0; i < dim_K; i++)
            for (slong c = 0; c < dim_H; c++) {
                double _Complex v = evecs[(i * dim_H + c) * n + a];
                /* Convention A requires conj on the (i,a) factor, so the reshape
                 * takes the CONJUGATE eigenvector component:
                 *   K_a[i,c] = sqrt(lambda_a) * conj(v_a[i*dim_H + c])
                 * so that kraus_to_choi(extracted) reproduces C. */
                acb_set_d_d(acb_mat_entry(phi->K[out_a], i, c),
                            scale * creal(v), -scale * cimag(v));
            }
        out_a++;
    }
    assert(out_a == r);

    free(evals);
    free(evecs);
    free(Ca);
}

slong aic_ucp_carrier_rank_latd(const acb_mat_t Q, slong dim_K)
{
    assert(acb_mat_nrows(Q) == dim_K && acb_mat_ncols(Q) == dim_K &&
           "carrier_rank_latd: Q must be dim_K x dim_K");
    double _Complex *Qa = malloc((size_t) (dim_K * dim_K) *
                                 sizeof(double _Complex));
    double *evals = malloc((size_t) dim_K * sizeof(double));
    assert(Qa && evals && "carrier_rank_latd: out of memory");

    aic_latd_from_acb_mat(Qa, Q);
    double thr = (double) dim_K * DBL_EPSILON * frob_norm(Qa, dim_K);
    aic_latd_eig_hermitian(evals, NULL, Qa, dim_K);

    slong rank = 0;
    for (slong a = 0; a < dim_K; a++)
        if (evals[a] > thr) rank++;

    free(evals);
    free(Qa);
    return rank;
}
