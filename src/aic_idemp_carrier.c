/* aic_idemp_carrier.c — the carrier M of an exactly-idempotent UCP map and its
 * isometries, on the LAPACK double path (bead aic-wuh). Step 1 of the
 * th_idemp_structure construction (approximate_algebras.tex:2056).
 *
 * Paper result vs constructive route. lem_carrier (.tex:1724) defines M as the
 * support of Phi^*(rho_0) and proves it the smallest M with M (x) F >= Im V — an
 * existence statement. CONSTRUCTIVE finite-dim route (matching aic_ucp_carrier):
 * M = range(Q), Q = sum_a K_a K_a^dag (the K-marginal of Im V). For an exactly
 * idempotent Phi, Q is Hermitian PSD and we split its spectrum at a numerical
 * threshold:
 *   J_M     = eigenvectors of Q with eigenvalue > thr  (n x dim_M),
 *   J_Mperp = eigenvectors with eigenvalue <= thr       (n x (n-dim_M)),
 *   Pi_M    = J_M J_M^dag.
 * thr = n * DBL_EPSILON * ||Q||_F (the QETLAB-style rank threshold used by
 * aic_ucp_carrier_rank_latd, so dim_M here == aic_ucp_carrier_rank_latd(Q)).
 *
 * Degenerate spectrum (why the double path). Q's eigenvalues cluster (an exactly
 * idempotent Phi has carrier eigenvalues that are NOT generally {0,1}, but the
 * kept block sits well above the dropped near-zero block); LAPACKE_zheev handles
 * the multiplicity. The certified arb degenerate eig is blocked on aic-w4o.1.
 *
 * Fail-loud (Rule 4): a clearly-negative Q eigenvalue (Q must be PSD) aborts; a
 * straddling eigenvalue (one within [thr/2, 2*thr], i.e. an unresolved rank gap)
 * aborts — a rank guess would corrupt every downstream map.
 */
#include <assert.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <flint/acb_mat.h>

#include "aic_idemp.h"
#include "aic_idemp_internal.h"
#include "aic_latd.h"
#include "aic_ucp.h"

static double frob_norm_arr(const double _Complex *A, slong n)
{
    double s = 0.0;
    for (slong k = 0; k < n * n; k++) {
        double m = cabs(A[k]);
        s += m * m;
    }
    return sqrt(s);
}

/* Fill J_M (n x dim_M), J_Mperp (n x (n-dim_M)), Pi_M (n x n) from the carrier
 * operator Q. Returns dim_M. evecs columns are eigenvectors, evals ascending, so
 * the KEPT (large) block is the TOP dim_M columns and M^perp is the bottom. */
slong aic_idemp_carrier_split(acb_mat_t J_M, acb_mat_t J_Mperp, acb_mat_t Pi_M,
                              const acb_mat_t Q, slong n, slong prec)
{
    assert(acb_mat_nrows(Q) == n && acb_mat_ncols(Q) == n);

    double _Complex *Qa = malloc((size_t) (n * n) * sizeof(double _Complex));
    double _Complex *evecs = malloc((size_t) (n * n) * sizeof(double _Complex));
    double *evals = malloc((size_t) n * sizeof(double));
    assert(Qa && evecs && evals && "carrier_split: out of memory");

    aic_latd_from_acb_mat(Qa, Q);
    double fro = frob_norm_arr(Qa, n);
    double thr = (double) n * DBL_EPSILON * fro;
    aic_latd_eig_hermitian(evals, evecs, Qa, n); /* ascending */

    slong dim_M = 0;
    for (slong a = 0; a < n; a++) {
        if (evals[a] < -thr) {
            flint_printf("aic_idemp_carrier_split: Q eigenvalue %g < -thr=%g; "
                         "Q = sum K K^dag must be PSD (Rule 4)\n", evals[a], -thr);
            flint_abort();
        }
        /* an eigenvalue in the gap band [thr/2, 2*thr] means the rank is not
         * resolved at this precision — refuse to guess (Rule 4). */
        if (evals[a] > thr / 2.0 && evals[a] < 2.0 * thr) {
            flint_printf("aic_idemp_carrier_split: Q eigenvalue %g straddles the "
                         "rank threshold thr=%g; carrier rank unresolved (Rule 4)\n",
                         evals[a], thr);
            flint_abort();
        }
        if (evals[a] > thr) dim_M++;
    }
    assert(dim_M >= 1 && "carrier_split: empty carrier (the zero map is not UCP)");
    assert(dim_M <= n);

    slong dperp = n - dim_M;
    acb_mat_init(J_M, n, dim_M);
    acb_mat_init(J_Mperp, n, dperp);   /* n x 0 when M = H (valid in FLINT) */
    acb_mat_init(Pi_M, n, n);

    /* KEPT block = top dim_M columns (largest eigenvalues, indices n-dim_M..n-1).
     * column c of J_M = eigenvector at evals index (n - dim_M + c). */
    for (slong c = 0; c < dim_M; c++) {
        slong col = n - dim_M + c;
        for (slong i = 0; i < n; i++)
            acb_set_d_d(acb_mat_entry(J_M, i, c),
                        creal(evecs[i * n + col]), cimag(evecs[i * n + col]));
    }
    /* M^perp = bottom dperp columns (indices 0..dperp-1). */
    for (slong c = 0; c < dperp; c++) {
        for (slong i = 0; i < n; i++)
            acb_set_d_d(acb_mat_entry(J_Mperp, i, c),
                        creal(evecs[i * n + c]), cimag(evecs[i * n + c]));
    }

    /* Pi_M = J_M J_M^dag. */
    acb_mat_t JMd;
    acb_mat_init(JMd, dim_M, n);
    acb_mat_conjugate_transpose(JMd, J_M);
    acb_mat_mul(Pi_M, J_M, JMd, prec);
    acb_mat_clear(JMd);

    free(evals);
    free(evecs);
    free(Qa);
    return dim_M;
}
