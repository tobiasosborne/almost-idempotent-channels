/* aic_ecstar_hopm.c — the double-path KERNEL for the Cycle-2 faithful worst-case
 * defect search (bead aic-knm, module "ecstar"). This file is the fast,
 * UNCERTIFIED inner loop (C99 double _Complex row-major arithmetic + LAPACKE SVD
 * via aic_latd) that the certified driver in aic_ecstar_search.c runs to FIND a
 * witness; the driver then re-evaluates the final witness in arb to certify the
 * lower bound. Split from the driver to keep each file <= 200 LOC (Rule 10).
 *
 * The method is scale-invariant HOPM (web leg "Cycle 2 ... method decision",
 * Method 1): see include/aic_ecstar.h for the algorithm contract and the
 * load-bearing SUBSPACE-POLAR ACCEPT-GUARD. The defects realized:
 *   ax_assoc   (.tex:412-413): h(X,Y,Z) = (X*Y)*Z - X*(Y*Z), trilinear;
 *   ax_prodnorm(.tex:410-411): f(X,Y) = X*Y, bilinear (drop the Z block);
 *   ax_C*      (.tex:427-428): g(X) = X^dag * X, MINIMIZED.
 * The star is X*Y = Phi(X.Y) with Phi(X) = sum_a K_a^dag X K_a (Heisenberg,
 * aic_ucp.h); here computed on midpoint doubles for speed.
 *
 * Storage: every matrix is a length-n*n double _Complex buffer, row-major
 * A[i*n+j] = entry(i,j) (the aic_latd convention). No certified balls here.
 */
#define LAPACK_COMPLEX_C99
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "aic_ecstar_hopm.h"
#include "aic_latd.h"

/* C = A*B (n x n, row-major). */
void aic_ehk_matmul(cx *C, const cx *A, const cx *B, slong n)
{
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            cx s = 0.0;
            for (slong k = 0; k < n; k++) s += A[i * n + k] * B[k * n + j];
            C[i * n + j] = s;
        }
}

/* out = Phi(X) = sum_a K_a^dag X K_a (Heisenberg). K[a] is the a-th Kraus op
 * (n x n row-major); r operators. Scratch t1,t2 are caller-provided n*n. */
void aic_ehk_phi(cx *out, const cx *const *K, slong r, const cx *X, slong n,
                 cx *t1, cx *t2)
{
    for (slong p = 0; p < n * n; p++) out[p] = 0.0;
    for (slong a = 0; a < r; a++) {
        /* t1 = K_a^dag X ; then t2 = t1 K_a ; out += t2 */
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++) {
                cx s = 0.0;
                for (slong k = 0; k < n; k++)
                    s += conj(K[a][k * n + i]) * X[k * n + j]; /* (K^dag)[i,k]=conj K[k,i] */
                t1[i * n + j] = s;
            }
        aic_ehk_matmul(t2, t1, K[a], n);
        for (slong p = 0; p < n * n; p++) out[p] += t2[p];
    }
}

/* out = X * Y = Phi(X.Y). Scratch xy,t1,t2 caller-provided n*n. */
void aic_ehk_star(cx *out, const aic_ehk *h, const cx *X, const cx *Y,
                  cx *xy, cx *t1, cx *t2)
{
    aic_ehk_matmul(xy, X, Y, h->n);
    aic_ehk_phi(out, (const cx *const *) h->K, h->r, xy, h->n, t1, t2);
}

/* out = h(X,Y,Z) = (X*Y)*Z - X*(Y*Z). Scratch a,b,xy,t1,t2 caller-provided. */
void aic_ehk_assoc(cx *out, const aic_ehk *h, const cx *X, const cx *Y,
                   const cx *Z, cx *a, cx *b, cx *xy, cx *t1, cx *t2)
{
    slong n = h->n;
    aic_ehk_star(a, h, X, Y, xy, t1, t2);   /* X*Y          */
    aic_ehk_star(b, h, a, Z, xy, t1, t2);   /* (X*Y)*Z -> b */
    aic_ehk_star(a, h, Y, Z, xy, t1, t2);   /* Y*Z -> a     */
    aic_ehk_star(out, h, X, a, xy, t1, t2); /* X*(Y*Z)      */
    for (slong p = 0; p < n * n; p++) out[p] = b[p] - out[p];
}

/* operator norm of an n x n matrix M (largest singular value), double path. */
double aic_ehk_opnorm(const cx *M, slong n)
{
    return aic_latd_opnorm(M, n, n);
}

/* y = M v (n x n times length-n vector). */
static void matvec(cx *y, const cx *M, const cx *v, slong n)
{
    for (slong i = 0; i < n; i++) {
        cx s = 0.0;
        for (slong k = 0; k < n; k++) s += M[i * n + k] * v[k];
        y[i] = s;
    }
}

/* y = M^dag v. */
static void matvec_dag(cx *y, const cx *M, const cx *v, slong n)
{
    for (slong i = 0; i < n; i++) {
        cx s = 0.0;
        for (slong k = 0; k < n; k++) s += conj(M[k * n + i]) * v[k];
        y[i] = s;
    }
}

/* Normalize a length-n vector to unit ell2; returns the old norm. */
static double vnorm(cx *v, slong n)
{
    double s = 0.0;
    for (slong k = 0; k < n; k++) s += creal(v[k] * conj(v[k]));
    s = sqrt(s);
    if (s > 1e-300)
        for (slong k = 0; k < n; k++) v[k] /= s;
    return s;
}

/* The u,v power-step for ||M||_op = max_{||u||=||v||=1} |<u, M v>|: given M and a
 * current v, set u <- Mv/||Mv|| then v <- M^dag u/||M^dag u||. Returns |<u,Mv>|
 * after the update (== the current Rayleigh estimate of ||M||_op). u,v are n. */
double aic_ehk_uv_step(cx *u, cx *v, const cx *M, slong n, cx *tmp)
{
    matvec(tmp, M, v, n);
    for (slong k = 0; k < n; k++) u[k] = tmp[k];
    vnorm(u, n);
    matvec_dag(tmp, M, u, n);
    for (slong k = 0; k < n; k++) v[k] = tmp[k];
    vnorm(v, n);
    /* <u, M v> with the updated v */
    matvec(tmp, M, v, n);
    cx ip = 0.0;
    for (slong k = 0; k < n; k++) ip += conj(u[k]) * tmp[k];
    return cabs(ip);
}

/* Project P = sum_k <B_k, P>_F B_k onto A, returning the coefficient vector in
 * cf[0..d-1] and writing the projected matrix into out (n*n). <B_k,P>_F =
 * sum_pq conj(B_k[p,q]) P[p,q]. */
void aic_ehk_project_A(cx *out, cx *cf, const aic_ehk *h, const cx *P)
{
    slong n = h->n, d = h->d;
    for (slong p = 0; p < n * n; p++) out[p] = 0.0;
    for (slong k = 0; k < d; k++) {
        cx c = 0.0;
        for (slong p = 0; p < n * n; p++) c += conj(h->B[k][p]) * P[p];
        cf[k] = c;
        for (slong p = 0; p < n * n; p++) out[p] += c * h->B[k][p];
    }
}

/* polar factor U V^dag of an n x n matrix C (full SVD via LAPACKE). Writes the
 * polar factor into pol (n*n). Scratch U,Vt are n*n; sv is double[n]. */
void aic_ehk_polar(cx *pol, const cx *C, slong n, cx *U, cx *Vt, double *sv)
{
    aic_latd_svd(sv, U, Vt, C, n, n);
    /* pol = U * Vt  (U columns are left vecs, Vt rows are v_k^dag) */
    aic_ehk_matmul(pol, U, Vt, n);
}
