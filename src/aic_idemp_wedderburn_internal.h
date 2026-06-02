/* aic_idemp_wedderburn_internal.h — private helpers for the Artin-Wedderburn
 * decomposition (bead aic-ynu, I1). Not part of the public API
 * (include/aic/aic_idemp.h). The double-path geometric extraction (eigen-cluster,
 * block connectivity, W_j assembly) lives in aic_idemp_wedderburn_assemble.c and
 * is driven by aic_idemp_wedderburn.c.
 *
 * Double arrays are ROW-MAJOR complex (the aic_latd convention). The output W_j
 * tensor layout is L-MAJOR: W_j row t = a*dim_E + c (a in L_j, c in E_j), matching
 * aic_mat_kronecker (A_j (x) 1_{E_j}).
 */
#ifndef AIC_IDEMP_WEDDERBURN_INTERNAL_H
#define AIC_IDEMP_WEDDERBURN_INTERNAL_H

#include <complex.h>

#include <flint/acb_mat.h>

#include "aic/aic_idemp.h"

/* ---- small double-path linear-algebra helpers (aic_idemp_wedderburn_ops.c) ---- */

/* <u, v> = sum_i conj(u_i) v_i over length m. */
double _Complex aic_wedd_cdot(const double _Complex *u, const double _Complex *v,
                              slong m);

/* y = A x for a row-major m x m complex A. */
void aic_wedd_matvec(double _Complex *y, const double _Complex *A,
                     const double _Complex *x, slong m);

/* gather column c of a row-major m x m matrix into contiguous f (length m). */
void aic_wedd_getcol(double _Complex *f, const double _Complex *evecs, slong c,
                     slong m);

/* Gap-cluster the ASCENDING eigenvalue array `evals` (length m) into runs of
 * near-equal values. Writes cl_start[a] / cl_size[a] (caller-allocated length m)
 * for cluster a and returns the number of clusters. A new cluster starts when the
 * eigenvalue gap exceeds an absolute+relative threshold. */
slong aic_wedd_cluster_evals(slong *cl_start, slong *cl_size, const double *evals,
                             slong m);

/* max_k ||Q_a w(B_k) Q_b||_op, the W-connectivity of clusters a (columns sa..+ka
 * of `evecs`) and b (sb..+kb). Zero across distinct Wedderburn blocks. */
double aic_wedd_block_link(const double _Complex *evecs, slong sa, slong ka,
                           slong sb, slong kb, double _Complex *const *Wk,
                           slong dA, slong m);

/* Certified-arb check of the FULL correctness spec of the assembled W:
 *   (i)   the dimension identities (sum dim_L*dim_E == dim_M, sum dim_L^2 == dim_A);
 *   (ii)  the stacked W is unitary (each W_j W_j^dag = 1, sum_j W_j^dag W_j = 1_M);
 *   (iii) the w-RECONSTRUCTION (the gauge-invariant correctness spec, NOT implied
 *         by unitarity): for each generator w(B_k), with A_j = Tr_{E_j}(M_blk)/dim
 *         E_j read by partial trace (the INDEPENDENT route),
 *           ||W_j w(B_k) W_j^dag - A_j (x) 1_{E_j}||_op  (tensor alignment) AND
 *           ||sum_j W_j^dag (A_j (x) 1_{E_j}) W_j - w(B_k)||_op  (round-trip)
 *         must be ~0. Wk is the dA double-array generators (same as `decompose`).
 * Fails loud (Rule 4) on any violation. */
void aic_wedd_certify_W(const aic_idemp_wedderburn *out,
                        double _Complex *const *Wk, slong dA, slong prec);

/* Assemble out->{num_blocks, dim_L, dim_E, W_j} from the *-algebra generators
 * {Wk} (dA arrays, each m x m row-major double), the Hermitian-eig onb `evecs`
 * (m x m row-major, columns = eigenvectors of the generic H_W), and the eigenvalue
 * clustering (cl_start/cl_size, n_cl clusters). Groups clusters into Wedderburn
 * blocks by W-connectivity, builds each W_j by aligning a single E_j onb across the
 * dim_L reference vectors via the connecting partial isometries in W, writes W_j as
 * acb_mat at `prec`, then certifies the full spec in arb (aic_wedd_certify_W).
 *
 * RE-DRAW PROTOCOL (the aic-ynu BLOCKER fix). If a draw is degenerate — a cluster
 * spans two would-be blocks, so the per-block cluster multiplicities disagree —
 * the caller can RE-DRAW H_W. When `fail_loud` is 0 such a draw returns nonzero
 * (out left empty, partial state cleaned) so the caller retries; when `fail_loud`
 * is 1 (the last draw in the budget) it aborts (Rule 4). Returns 0 on a clean,
 * certified decomposition. */
int aic_wedd_assemble(aic_idemp_wedderburn *out, double _Complex *const *Wk,
                      slong dA, const double _Complex *evecs,
                      const slong *cl_start, const slong *cl_size, slong n_cl,
                      slong m, slong prec, int fail_loud);

#endif /* AIC_IDEMP_WEDDERBURN_INTERNAL_H */
