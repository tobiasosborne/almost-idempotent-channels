/* aic_idemp_image.c — the algebra A = Img Phi as an orthonormal operator basis,
 * on the LAPACK double path (bead aic-wuh). Step 2 of th_idemp_structure
 * (approximate_algebras.tex:2056: "A = Img Phi").
 *
 * Constructive route (the COLUMN-IMAGE SVD, not the superoperator matrix). Img
 * Phi is the column space of Phi viewed as a linear operator on B(C^n) ~ C^{n^2}.
 * We never form the n^2 x n^2 superoperator (that route risks a vectorization-
 * convention bug); instead we apply the TESTED aic_ucp_apply to each matrix unit
 * E_{ij} (i,j in [0,n)) and stack vec(Phi(E_ij)) (row-major, vec(Y)[a*n+b]=Y[a,b])
 * as the columns of an n^2 x n^2 matrix P. Then
 *     Img Phi = column space of P = span of the left singular vectors of P with
 *               nonzero singular value,
 * so an ORTHONORMAL basis {B_k} of A is the reshape of the left singular vectors
 * u_k (sigma_k > thr) back to n x n matrices: B_k[a,b] = u_k[a*n+b]. dim_A is the
 * numerical rank (count of sigma_k > thr).
 *
 * Why an ONB. Delta = [vec B_1|...|vec B_d] then has Delta^dag Delta = 1_A, which
 * makes the Lambda = Delta Gamma factorization solvable by Gamma = Delta^dag
 * Lambda (no pseudoinverse) and makes Pi_A = Delta Delta^dag the gauge-invariant
 * subspace projector for cross-checks (the basis itself is unique only up to a
 * unitary on C^{dim_A}).
 *
 * Degenerate spectrum (why the double path). P's nonzero singular values are
 * generically repeated (a projection-like image has many equal singular values),
 * so LAPACKE_zgesvd (degeneracy-robust) is used; the certified arb SVD is blocked
 * on aic-w4o.1.
 *
 * Fail-loud (Rule 4): a singular value in the gap band [thr/2, 2*thr] (unresolved
 * rank) aborts.
 */
#include <assert.h>
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <flint/acb_mat.h>

#include "aic_idemp_internal.h"
#include "aic/aic_latd.h"
#include "aic/aic_ucp.h"

slong aic_idemp_image_basis(acb_mat_t Delta, const aic_ucp_kraus *phi,
                            slong n, slong prec)
{
    slong n2 = n * n;
    assert(phi->dim_K == n && phi->dim_H == n && "image_basis: Phi must be a self-map");

    /* Build P (n^2 x n^2): column (i*n+j) = vec(Phi(E_ij)). */
    double _Complex *P = malloc((size_t) (n2 * n2) * sizeof(double _Complex));
    double _Complex *U = malloc((size_t) (n2 * n2) * sizeof(double _Complex));
    double *svals = malloc((size_t) n2 * sizeof(double));
    assert(P && U && svals && "image_basis: out of memory");

    acb_mat_t E, PhiE;
    acb_mat_init(E, n, n);
    acb_mat_init(PhiE, n, n);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            acb_mat_zero(E);
            acb_set_si(acb_mat_entry(E, i, j), 1);       /* E_ij = |e_i><e_j| */
            aic_ucp_apply(PhiE, phi, E, prec);           /* tested path */
            slong col = i * n + j;
            for (slong a = 0; a < n; a++)
                for (slong b = 0; b < n; b++) {
                    /* vec row index a*n+b; P[row*n2 + col] (row-major P). */
                    acb_srcptr e = acb_mat_entry(PhiE, a, b);
                    P[(a * n + b) * n2 + col] =
                        arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR) +
                        arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR) * I;
                }
        }
    acb_mat_clear(PhiE);
    acb_mat_clear(E);

    /* ||P||_F for the threshold. */
    double fro = 0.0;
    for (slong k = 0; k < n2 * n2; k++) { double m = cabs(P[k]); fro += m * m; }
    fro = sqrt(fro);
    double thr = (double) n2 * DBL_EPSILON * fro;

    /* SVD: U columns = left singular vectors, svals DESCENDING. */
    aic_latd_svd(svals, U, NULL, P, n2, n2);

    slong dim_A = 0;
    for (slong k = 0; k < n2; k++) {
        if (svals[k] > thr / 2.0 && svals[k] < 2.0 * thr) {
            flint_printf("aic_idemp_image_basis: singular value %g straddles the "
                         "rank threshold thr=%g; dim A unresolved (Rule 4)\n",
                         svals[k], thr);
            flint_abort();
        }
        if (svals[k] > thr) dim_A++;
    }
    assert(dim_A >= 1 && "image_basis: empty image (Phi(1)=1 forces dim_A>=1)");

    acb_mat_init(Delta, n2, dim_A);
    /* column k of Delta = left singular vector u_k (the TOP dim_A, descending). */
    for (slong k = 0; k < dim_A; k++)
        for (slong row = 0; row < n2; row++)
            acb_set_d_d(acb_mat_entry(Delta, row, k),
                        creal(U[row * n2 + k]), cimag(U[row * n2 + k]));

    free(svals);
    free(U);
    free(P);
    return dim_A;
}
