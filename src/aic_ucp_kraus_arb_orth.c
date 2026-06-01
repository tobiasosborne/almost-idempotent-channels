/* aic_ucp_kraus_arb_orth.c — the keep-threshold ball, the per-cluster LOEWDIN
 * orthonormalisation, and the per-eigenvalue cluster eigenbasis for the certified
 * Choi->Kraus extraction (bead aic-4td increment 2 step D). Split out of
 * aic_ucp_kraus_arb.c for Rule 10 (~200 LOC/file). The full rationale (Convention A,
 * the R2 orthonormalisation subtlety, the x0 = ||X||_op^2 / FINDINGS §C5 interval-Gram
 * hazard) is in that file's docstring; this file holds the numerical primitives.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_funcalc.h"
#include "aic/aic_mat.h"
#include "aic_ucp_kraus_arb_internal.h"

void aic_ucp_kraus_arb_int_keep_thr(arb_t thr, const acb_mat_t C, slong n,
                                    slong prec)
{
    aic_mat_frobenius_norm(thr, C, prec);
    arb_mul_2exp_si(thr, thr, -52);          /* * 2^-52 (DBL_EPSILON) */
    arb_mul_si(thr, thr, n, prec);           /* * (dim_K*dim_H) */
}

void aic_ucp_kraus_arb_int_loewdin(acb_mat_t V, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X), k = acb_mat_ncols(X);
    if (k == 1) {
        /* V = X / sqrt(X^dag X) (a single column normalisation). */
        arb_t g;
        arb_init(g);
        acb_t acc, t;
        acb_init(acc);
        acb_init(t);
        for (slong i = 0; i < n; i++) {
            acb_conj(t, acb_mat_entry(X, i, 0));
            acb_addmul(acc, t, acb_mat_entry(X, i, 0), prec); /* sum |X_i|^2 */
        }
        acb_get_real(g, acc);                 /* ||X||^2 (real, PD) */
        arb_sqrt(g, g, prec);                 /* ||X|| */
        arb_inv(g, g, prec);
        acb_mat_scalar_mul_arb(V, X, g, prec);
        acb_clear(t);
        acb_clear(acc);
        arb_clear(g);
        return;
    }
    /* G = X^dag X (k x k Hermitian PD); Ginvsqrt = G^{-1/2}; V = X Ginvsqrt. */
    acb_mat_t Xd, G, Ginvsqrt;
    acb_mat_init(Xd, k, n);
    acb_mat_init(G, k, k);
    acb_mat_init(Ginvsqrt, k, k);
    acb_mat_conjugate_transpose(Xd, X);
    acb_mat_mul(G, Xd, X, prec);
    /* x0 = lambda_max(G) = ||X||_op^2 makes ||G/x0 - I||_op = 1 - l_min/l_max < 1,
     * so the xpow binomial domain (tex:511) always holds for PD G. We get
     * lambda_max(G) as ||X||_op^2 via aic_mat_opnorm (which Weyl-Hermitianizes the
     * interval Gram internally, FINDINGS §C5) — NOT aic_mat_herm_max_eig(G), which
     * would ABORT on the raw interval Gram X^dag X (entry (i,j) and conj(j,i) carry
     * independent radii; the tight Hermiticity assert false-fails). */
    arb_t x0;
    arb_init(x0);
    aic_mat_opnorm(x0, X, prec);             /* ||X||_op */
    arb_sqr(x0, x0, prec);                    /* lambda_max(G) = ||X||_op^2 */
    double x0d = arf_get_d(arb_midref(x0), ARF_RND_UP) + mag_get_d(arb_radref(x0));
    aic_funcalc_xpow(Ginvsqrt, G, -0.5, x0d, prec);
    acb_mat_mul(V, X, Ginvsqrt, prec);
    arb_clear(x0);
    acb_mat_clear(Ginvsqrt);
    acb_mat_clear(G);
    acb_mat_clear(Xd);
}

/* Hermitian-project the MIDPOINTS of a (genuinely Hermitian) interval matrix A into
 * an exact-Hermitian out (FINDINGS §C5): the two matmuls forming M = V^dag C V carry
 * independent per-entry radii that trip the tight Hermiticity assert of the recursed
 * subspace eig; the midpoint projection is the standard cure and changes no value on
 * the true Hermitian matrix. */
static void herm_project_mid(acb_mat_t out, const acb_mat_t A)
{
    slong n = acb_mat_nrows(A);
    for (slong i = 0; i < n; i++) {
        acb_get_mid(acb_mat_entry(out, i, i), acb_mat_entry(A, i, i));
        arb_zero(acb_imagref(acb_mat_entry(out, i, i)));
        for (slong j = i + 1; j < n; j++) {
            acb_get_mid(acb_mat_entry(out, i, j), acb_mat_entry(A, i, j));
            acb_conj(acb_mat_entry(out, j, i), acb_mat_entry(out, i, j));
        }
    }
}

/* aic_ucp_kraus_arb_int_cluster_basis — FINDING-2 FIX (the lumped-distinct-cluster
 * bug, inc-2 hostile review). The strict route scaled all k orthonormal columns of a
 * cluster by sqrt(lambda_c) (the SINGLE cluster ball). When gap-clustering LUMPS two
 * DISTINCT nonzero eigenvalues (gap < gap_thr) into one k>=2 cluster, lambda_c=mean is
 * WRONG for both -> silently-wrong Kraus (round-trip error ~gap; measured 5.7e-7 at
 * gap 5e-7, NOT contained in the certified ball). FIX (design §3.2 Option B,
 * EXACT even when a cluster lumps distinct eigenvalues): recover the per-eigenvalue
 * structure from the small DENSE k x k compression M = V^dag C V (V = Loewdin basis
 * of the cluster's range, n x k). M is Hermitian PSD with spec(M) = the cluster's
 * eigenvalues; diagonalising it splits a lumped pair. Build:
 *   V = Loewdin(X)                              (n x k, orthonormal range basis)
 *   M = herm_mid(V^dag C V)                     (k x k, exact-Hermitian, FINDINGS §C5)
 *   {mu_i, w_i} = subspace-eig(M, gap_thr=1e-10*(1+||M||_F)) (degeneracy-robust, certified)
 *   per sub-cluster: W = Loewdin(sub.X) (k x k_sub orthonormal); VW = V W (n x k_sub);
 *     each column carries scale sqrt(mu_i).
 * Writes the assembled orthonormal eigenbasis into VW_out (n x k, caller-init'd) in
 * sub-eigenvalue order, and the per-column sqrt(eigenvalue) into scales (caller arb[k]).
 * gap_thr = 1e-10*(1+||M||_F) splits genuinely-distinct sub-eigenvalues while keeping
 * a genuine degeneracy as ONE sub-cluster. The band is wide+measured (probe, prec=128):
 * gap_thr in [1e-14, 1e-7] keeps (1/3)I_9 as one cluster AND splits {0.5,0.5+5e-7}
 * (gap 5e-7) and {0.5,0.5+1e-8}. The LOWER bound matters: the double-precision zheev
 * SEED inside the recursion has eigenvalue jitter ~ n*DBL_EPSILON*||M|| ~ 1e-15*||M||,
 * so a gap_thr BELOW that (e.g. 2^-(prec/2)=5.4e-20 at prec=128) spuriously SPLITS a
 * genuine degeneracy into k=1 sub-clusters whose degenerate-direction Rump enclosure
 * goes NON-FINITE -> the recursion fails loud UNRESOLVED. 1e-10*(1+||M||_F) sits 5
 * orders above the seed floor and >=3 orders below the project's distinct separations
 * (measured: {0.5,0.5}->k=2, {0.7,0.7,0.2}->{0.2},{0.7}, lumped {0.5,0.5+5e-7}->two k=1).
 * For k=1 callers use the single-eigenvalue path directly.
 * Soundness: a kept outer cluster's ball is certified > keep_thr, so every sub-mu_i
 * (>= the ball's lower bound) is also > keep_thr — no sub-column is ever dropped, so
 * the PASS-1 rank count (k per kept cluster) stays exact. The recursed eig fail-loud
 * guards (finite enclosure, disjointness, the on-H residual) carry over to M. */
void aic_ucp_kraus_arb_int_cluster_basis(acb_mat_t VW_out, arb_ptr scales,
                                         const acb_mat_t X, const acb_mat_t C,
                                         slong prec)
{
    slong n = acb_mat_nrows(X), k = acb_mat_ncols(X);
    acb_mat_t V, Vd, CV, M, Mh;
    acb_mat_init(V, n, k);
    acb_mat_init(Vd, k, n);
    acb_mat_init(CV, n, k);
    acb_mat_init(M, k, k);
    acb_mat_init(Mh, k, k);
    aic_ucp_kraus_arb_int_loewdin(V, X, prec);   /* orthonormal range basis (n x k) */
    acb_mat_conjugate_transpose(Vd, V);
    acb_mat_mul(CV, C, V, prec);
    acb_mat_mul(M, Vd, CV, prec);                /* M = V^dag C V (spec = cluster evals) */
    herm_project_mid(Mh, M);                     /* exact-Hermitian (FINDINGS §C5) */

    /* gap_thr = 1e-10*(1+||M||_F): above the ~1e-15*||M|| double-seed jitter (else a
     * genuine degeneracy splits into non-finite k=1 Rump enclosures -> UNRESOLVED),
     * below the project's distinct separations (see the docstring band). */
    arb_t mf;
    arb_init(mf);
    aic_mat_frobenius_norm(mf, Mh, prec);
    double mfd = arf_get_d(arb_midref(mf), ARF_RND_UP) + mag_get_d(arb_radref(mf));
    arb_clear(mf);
    double gap_thr = 1e-10 * (1.0 + mfd);
    aic_mat_eigcluster *sc;
    slong snc;
    aic_mat_eig_hermitian_subspaces(&sc, &snc, Mh, gap_thr, prec);

    slong col = 0;
    for (slong s = 0; s < snc; s++) {
        slong sk = sc[s].k;
        arb_t mu, scl;
        arb_init(mu);
        arb_init(scl);
        acb_get_real(mu, sc[s].lambda);
        arb_sqrt(scl, mu, prec);                 /* sqrt(mu_i) (mu_i > keep_thr > 0) */
        acb_mat_t W, VWs;
        acb_mat_init(W, k, sk);
        acb_mat_init(VWs, n, sk);
        aic_ucp_kraus_arb_int_loewdin(W, sc[s].X, prec);  /* orthonormal sub-basis */
        acb_mat_mul(VWs, V, W, prec);            /* direction in C^n = V W */
        for (slong m = 0; m < sk; m++) {
            for (slong row = 0; row < n; row++)
                acb_set(acb_mat_entry(VW_out, row, col),
                        acb_mat_entry(VWs, row, m));
            arb_set(&scales[col], scl);
            col++;
        }
        acb_mat_clear(VWs);
        acb_mat_clear(W);
        arb_clear(scl);
        arb_clear(mu);
    }
    /* sub-clusters partition the k-dim cluster space (Sum sk == k, the subspace
     * routine's invariant); a mismatch is a silently-wrong assembly -> fail loud. */
    if (col != k) {
        fprintf(stderr,
                "aic_ucp_choi_to_kraus_arb: cluster eigenbasis assembled %ld of %ld "
                "columns (sub-cluster dims do not sum to k) at prec=%ld; bead "
                "aic-4td, finding-2\n",
                (long) col, (long) k, (long) prec);
        abort();
    }
    aic_mat_eigcluster_free(sc, snc);
    acb_mat_clear(Mh);
    acb_mat_clear(M);
    acb_mat_clear(CV);
    acb_mat_clear(Vd);
    acb_mat_clear(V);
}
