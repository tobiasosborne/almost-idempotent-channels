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
 * TOLERANCE-AWARE variant (aic_ucp_choi_to_kraus_latd_tol; FINDINGS §C14). The
 * paper's almost-idempotent CP-ization Delta'(X)=sum_s p_s Phi(Dtilde(XU_s^dag)
 * Dtilde(U_s)) (.tex:2786-2789) is manifestly CP (.tex:2791-2796) ONLY when
 * Dtilde is an EXACT unital homomorphism; for the extended O(eta)-isomorphism
 * Dtilde=v the manifest-positive form holds to O(eta), so the resulting UCP Delta
 * has a per-block Choi that is PSD only to O(eta^2) — a small NEGATIVE eigenvalue
 * of order O(eta^2) (measured worst block-Choi minEig = -2.5e-6 at eta=2.3e-2,
 * -1.5e-7 at eta=7.9e-3, -2.9e-8 at eta=2.4e-3: minEig/eta^2 ~ -0.005, NOT
 * O(eta)). The strict gate above (thr ~ 1e-15) rejects that as not-CP and aborts.
 * The _tol variant admits eigenvalues in (-neg_tol, keep_thr] as cone-defect /
 * noise (DROPPING them — the PSD-CONE PROJECTION onto the nearest genuinely-CP
 * map), accumulates the clipped negative mass as the certified cone-defect
 * magnitude, and STILL fails loud for lambda <= -neg_tol (a genuine non-CP input,
 * a real bug). The strict default delegates with neg_tol = keep_thr, so every
 * existing caller is byte-for-byte unchanged.
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

void aic_ucp_choi_to_kraus_latd_tol(aic_ucp_kraus *phi, const acb_mat_t C,
                                    slong dim_K, slong dim_H, double neg_tol,
                                    double *clipped_neg_out)
{
    slong n = dim_K * dim_H;
    assert(acb_mat_nrows(C) == n && acb_mat_ncols(C) == n &&
           "choi_to_kraus_latd: C must be (dim_K*dim_H) square");
    assert(neg_tol >= 0.0 && "choi_to_kraus_latd_tol: neg_tol must be >= 0");

    double _Complex *Ca = malloc((size_t) (n * n) * sizeof(double _Complex));
    double _Complex *evecs = malloc((size_t) (n * n) * sizeof(double _Complex));
    double *evals = malloc((size_t) n * sizeof(double));
    assert(Ca && evecs && evals && "choi_to_kraus_latd: out of memory");

    aic_latd_from_acb_mat(Ca, C);
    double thr = (double) n * DBL_EPSILON * frob_norm(Ca, n); /* keep_thr */
    aic_latd_eig_hermitian(evals, evecs, Ca, n); /* ascending */

    /* THREE eigenvalue regimes (ascending sweep; FINDINGS §C14):
     *   lambda > keep_thr            -> KEEP as Kraus (PSD signal);
     *   lambda in (-neg_tol, keep_thr] -> DROP (numerically-zero tail OR an
     *                                   approximately-CP-cone DEFECT). If lambda<0
     *                                   accumulate |lambda| into the clipped mass
     *                                   (the certified cone-defect magnitude) — this
     *                                   is the PSD-CONE PROJECTION: we keep only the
     *                                   nonnegative part, i.e. the NEAREST genuinely-CP
     *                                   map to C (within the clipped mass of C);
     *   lambda <= -neg_tol           -> a GENUINELY non-PSD (not-CP) input: a real
     *                                   bug or a non-CP argument. Still FAIL LOUD
     *                                   (Rule 4) — sqrt'ing it would make NaN Kraus,
     *                                   and corrupted output is worse than a crash. */
    slong r = 0;
    double clipped = 0.0;
    for (slong a = 0; a < n; a++) {
        if (evals[a] <= -neg_tol) {
            flint_printf("aic_ucp_choi_to_kraus_latd_tol: eigenvalue %g <= "
                         "-neg_tol=%g; C is not PSD (not CP) beyond the admitted "
                         "cone tolerance — cannot extract Kraus ops (Rule 4: fail "
                         "loud)\n", evals[a], -neg_tol);
            flint_abort();
        }
        if (evals[a] > thr)
            r++;
        else if (evals[a] < 0.0)
            clipped += -evals[a];           /* cone-defect mass (clipped to 0) */
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

    if (clipped_neg_out)
        *clipped_neg_out = clipped;

    free(evals);
    free(evecs);
    free(Ca);
}

void aic_ucp_choi_to_kraus_latd(aic_ucp_kraus *phi, const acb_mat_t C,
                                slong dim_K, slong dim_H)
{
    /* STRICT default: delegate with neg_tol = keep_thr (the historical
     * (dim_K*dim_H)*eps_mach*||C||_F gate). For lambda <= -keep_thr the _tol
     * routine fails loud EXACTLY as the old strict gate did; no caller of this
     * symbol sees any behaviour change. */
    slong n = dim_K * dim_H;
    double _Complex *Ca = malloc((size_t) (n * n) * sizeof(double _Complex));
    assert(Ca && "choi_to_kraus_latd: out of memory");
    aic_latd_from_acb_mat(Ca, C);
    double thr = (double) n * DBL_EPSILON * frob_norm(Ca, n);
    free(Ca);
    aic_ucp_choi_to_kraus_latd_tol(phi, C, dim_K, dim_H, thr, NULL);
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
