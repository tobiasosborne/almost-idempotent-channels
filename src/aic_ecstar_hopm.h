/* aic_ecstar_hopm.h — internal double-path kernel interface for the Cycle-2 HOPM
 * search (bead aic-knm). NOT a public header (lives in src/, not include/): the
 * public contract is aic_ecstar_defect_*_hopm in include/aic_ecstar.h. This
 * exposes the row-major double _Complex primitives in aic_ecstar_hopm.c so the
 * driver in aic_ecstar_search.c (and the kernel itself) share them. The block
 * update rules and the polar/project accept-guard live in aic_ecstar_search.c;
 * this header is the matrix arithmetic only.
 */
#ifndef AIC_ECSTAR_HOPM_H
#define AIC_ECSTAR_HOPM_H

#include <complex.h>

#include <flint/flint.h> /* slong */

typedef double _Complex cx;

/* Double-path snapshot of an aic_ecstar: Kraus ops K[a] (n x n row-major) define
 * Phi(X)=sum K_a^dag X K_a, and B[k] (n x n row-major) is the Frobenius-ONB of A.
 * All buffers are BORROWED (owned by the driver). */
typedef struct {
    slong n;            /* A <= M_n                         */
    slong d;            /* dim A                            */
    slong r;            /* number of Kraus ops              */
    cx **K;             /* r Kraus ops, each n*n row-major  */
    cx **B;             /* d basis ops,  each n*n row-major */
} aic_ehk;

/* C = A*B (n x n). */
void aic_ehk_matmul(cx *C, const cx *A, const cx *B, slong n);

/* out = Phi(X) = sum_a K_a^dag X K_a. t1,t2 scratch n*n. */
void aic_ehk_phi(cx *out, const cx *const *K, slong r, const cx *X, slong n,
                 cx *t1, cx *t2);

/* out = X*Y = Phi(X.Y). xy,t1,t2 scratch n*n. */
void aic_ehk_star(cx *out, const aic_ehk *h, const cx *X, const cx *Y,
                  cx *xy, cx *t1, cx *t2);

/* out = h(X,Y,Z) = (X*Y)*Z - X*(Y*Z). a,b,xy,t1,t2 scratch n*n. */
void aic_ehk_assoc(cx *out, const aic_ehk *h, const cx *X, const cx *Y,
                   const cx *Z, cx *a, cx *b, cx *xy, cx *t1, cx *t2);

/* operator norm (largest singular value) of an n x n matrix M. */
double aic_ehk_opnorm(const cx *M, slong n);

/* u,v power-step for ||M||_op; returns the Rayleigh estimate |<u,Mv>|. tmp n. */
double aic_ehk_uv_step(cx *u, cx *v, const cx *M, slong n, cx *tmp);

/* Project P onto A: out = sum_k <B_k,P>_F B_k, coeffs into cf[0..d-1]. */
void aic_ehk_project_A(cx *out, cx *cf, const aic_ehk *h, const cx *P);

/* polar factor U V^dag of C (full SVD). U,Vt scratch n*n; sv double[n]. */
void aic_ehk_polar(cx *pol, const cx *C, slong n, cx *U, cx *Vt, double *sv);

#endif /* AIC_ECSTAR_HOPM_H */
