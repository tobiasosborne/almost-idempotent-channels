/* aic_idemp_wedderburn_ops.c — small linear-algebra helpers for the Artin-
 * Wedderburn decomposition (bead aic-ynu, I1): the double-path eigenvalue
 * gap-clustering, and the FULL certified-arb correctness check of the assembled
 * W = [W_0; ...; W_{m-1}] — the dim identities, W-unitarity, AND the
 * w-reconstruction spec (FIX 2: a tensor-MISALIGNED-but-unitary W is rejected in
 * production, not only by the tests). Double path for clustering; arb for every
 * certificate. The inter-cluster connectivity probe (aic_wedd_block_link) lives in
 * aic_idemp_wedderburn_assemble.c with the union-find that consumes it. See
 * aic_idemp_wedderburn.c for the algorithm narrative (.tex:259-272, :2106).
 */
#include <assert.h>
#include <complex.h>
#include <math.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic_idemp_wedderburn_internal.h"

/* <u, v> = sum_i conj(u_i) v_i over length m. */
double _Complex aic_wedd_cdot(const double _Complex *u, const double _Complex *v,
                              slong m)
{
    double _Complex s = 0;
    for (slong i = 0; i < m; i++) s += conj(u[i]) * v[i];
    return s;
}

/* y = A x for a row-major m x m complex A. */
void aic_wedd_matvec(double _Complex *y, const double _Complex *A,
                     const double _Complex *x, slong m)
{
    for (slong i = 0; i < m; i++) {
        double _Complex s = 0;
        for (slong j = 0; j < m; j++) s += A[i * m + j] * x[j];
        y[i] = s;
    }
}

/* gather column c of a row-major m x m matrix into contiguous f (length m). */
void aic_wedd_getcol(double _Complex *f, const double _Complex *evecs, slong c,
                     slong m)
{
    for (slong i = 0; i < m; i++) f[i] = evecs[i * m + c];
}

slong aic_wedd_cluster_evals(slong *cl_start, slong *cl_size, const double *evals,
                             slong m)
{
    slong n_cl = 0, run0 = 0;
    double span = fabs(evals[m - 1] - evals[0]);
    /* gap threshold = 1e-9*(1 + span): absolute (1e-9) + relative (1e-9*span) so a
     * degenerate eigenvalue cluster (equal to the last zheev digit, ~1e-14*span)
     * stays merged while genuinely distinct eigenvalues separate. At eta=0 the H_W
     * draw is generic, so within-cluster eigenvalues are EQUAL (a true E_j
     * multiplicity) to ~1e-14*span while distinct eigenvalues differ by O(span);
     * 1e-9 sits in that gulf. Scaled by span (not prec) because the resolvable gap
     * is set by the zheev backward error relative to the spectral spread. */
    double gap_thr = 1e-9 * (1.0 + span);
    for (slong i = 1; i <= m; i++) {
        if (i == m || (evals[i] - evals[i - 1]) > gap_thr) {
            cl_start[n_cl] = run0;
            cl_size[n_cl] = i - run0;
            n_cl++;
            run0 = i;
        }
    }
    return n_cl;
}

/* w(B_k) as an m x m acb_mat from the row-major double array Wk[k]. */
static void wedd_arr_to_acb(acb_mat_t W, const double _Complex *Wk, slong m)
{
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            acb_set_d_d(acb_mat_entry(W, i, j),
                        creal(Wk[i * m + j]), cimag(Wk[i * m + j]));
}

/* Certify (arb) the w-RECONSTRUCTION correctness spec (FIX 2): NOT implied by
 * W-unitarity, so a tensor-MISALIGNED-but-unitary W_j (an L<->E swap) must be
 * caught HERE, in production. For each generator w(B_k):
 *   M_blk(j) = W_j w(B_k) W_j^dag (db x db); A_j = Tr_{E_j}(M_blk)/dim E_j read by
 *   partial trace (the SAME independent route the test uses, .tex:277). Two teeth:
 *   ||M_blk - A_j (x) 1_{E_j}||_op (tensor alignment) and the round-trip residual
 *   ||sum_j W_j^dag (A_j (x) 1_{E_j}) W_j - w(B_k)||_op (.tex:268). Fails loud on
 *   any defect above `tol`. (Spans a basis of A, so all generators are covered.) */
static void wedd_certify_w_recon(const aic_idemp_wedderburn *out,
                                 double _Complex *const *Wk, slong dA,
                                 const arb_t tol, slong prec)
{
    slong nb = out->num_blocks, m = out->dim_M;
    arb_t e;
    arb_init(e);
    acb_mat_t wbk, recon, diff, Wjd, tmp, Mblk, Aj, one_E, AjXI, td, t2, back;
    acb_mat_init(wbk, m, m);
    acb_mat_init(recon, m, m);
    acb_mat_init(diff, m, m);
    for (slong k = 0; k < dA; k++) {
        wedd_arr_to_acb(wbk, Wk[k], m);
        acb_mat_zero(recon);
        for (slong j = 0; j < nb; j++) {
            slong dL = out->dim_L[j], dE = out->dim_E[j], db = dL * dE;
            acb_mat_init(Wjd, m, db);
            acb_mat_conjugate_transpose(Wjd, out->W_j[j]);
            acb_mat_init(tmp, db, m);
            acb_mat_mul(tmp, out->W_j[j], wbk, prec);   /* W_j w(B_k)         */
            acb_mat_init(Mblk, db, db);
            acb_mat_mul(Mblk, tmp, Wjd, prec);          /* W_j w(B_k) W_j^dag */
            acb_mat_init(Aj, dL, dL);
            aic_mat_partial_trace_right(Aj, Mblk, dL, dE, prec);
            { acb_t inv; acb_init(inv); acb_set_d(inv, 1.0 / (double) dE);
              acb_mat_scalar_mul_acb(Aj, Aj, inv, prec); acb_clear(inv); }
            acb_mat_init(one_E, dE, dE);
            acb_mat_one(one_E);
            acb_mat_init(AjXI, db, db);
            aic_mat_kronecker(AjXI, Aj, one_E, prec);    /* A_j (x) 1_{E_j}    */
            acb_mat_init(td, db, db);
            acb_mat_sub(td, Mblk, AjXI, prec);           /* tensor-alignment   */
            aic_mat_opnorm(e, td, prec);
            if (!arb_le(e, tol)) {
                flint_printf("aic_idemp_wedderburn: W_%wd tensor-MISALIGNED "
                             "(||W_j w(B_k) W_j^dag - A_j (x) 1_E||_op too large at "
                             "k=%wd); fail loud (Rule 4)\n", j, k);
                flint_abort();
            }
            acb_mat_init(t2, m, db);
            acb_mat_mul(t2, Wjd, AjXI, prec);
            acb_mat_init(back, m, m);
            acb_mat_mul(back, t2, out->W_j[j], prec);    /* W_j^dag(..)W_j     */
            acb_mat_add(recon, recon, back, prec);
            acb_mat_clear(back); acb_mat_clear(t2); acb_mat_clear(td);
            acb_mat_clear(AjXI); acb_mat_clear(one_E); acb_mat_clear(Aj);
            acb_mat_clear(Mblk); acb_mat_clear(tmp); acb_mat_clear(Wjd);
        }
        acb_mat_sub(diff, recon, wbk, prec);
        aic_mat_opnorm(e, diff, prec);                   /* round-trip         */
        if (!arb_le(e, tol)) {
            flint_printf("aic_idemp_wedderburn: w-reconstruction round-trip failed "
                         "(||sum_j W_j^dag (A_j (x) 1) W_j - w(B_k)||_op too large "
                         "at k=%wd); fail loud (Rule 4)\n", k);
            flint_abort();
        }
    }
    acb_mat_clear(diff);
    acb_mat_clear(recon);
    acb_mat_clear(wbk);
    arb_clear(e);
}

/* Certify (arb) the dimension identities, that the stacked W is unitary, AND the
 * w-reconstruction correctness spec (FIX 2):
 *   sum_j dim_L*dim_E == dim_M,  sum_j dim_L^2 == dim_A (EXACT),
 *   each W_j W_j^dag = 1_{block} and sum_j W_j^dag W_j = 1_M, to a prec tol,
 *   and w(B_k) = sum_j W_j^dag (A_j (x) 1_{E_j}) W_j with A_j independently read by
 *   partial trace (wedd_certify_w_recon above). Wk are the dA double generators. */
void aic_wedd_certify_W(const aic_idemp_wedderburn *out,
                        double _Complex *const *Wk, slong dA, slong prec)
{
    slong nb = out->num_blocks, m = out->dim_M;
    slong sumLE = 0, sumL2 = 0;
    for (slong j = 0; j < nb; j++) {
        sumLE += out->dim_L[j] * out->dim_E[j];
        sumL2 += out->dim_L[j] * out->dim_L[j];
    }
    if (sumLE != m || sumL2 != out->dim_A) {
        flint_printf("aic_idemp_wedderburn: dim identity failed "
                     "(sum dim_L*dim_E=%wd vs dim_M=%wd; sum dim_L^2=%wd vs "
                     "dim_A=%wd); fail loud (Rule 4)\n",
                     sumLE, m, sumL2, out->dim_A);
        flint_abort();
    }

    /* Certificate tol 1e-9: the W_j rows are double-path eigenvectors (zheev,
     * orthonormal to ~1e-14), so the arb residuals (||W_j W_j^dag - 1||,
     * ||sum W_j^dag W_j - 1||, and the w-reconstruction defects) floor at the
     * extraction error ~1e-13, while a WRONG W (misaligned tensor factor, mis-split
     * dims) gives an O(1) defect. 1e-9 keeps full teeth with ~4 orders of headroom
     * over the clean floor; not 2^-prec because the floor is the double extraction,
     * not the arb working precision. */
    arb_t tol, e;
    arb_init(tol);
    arb_init(e);
    arb_set_d(tol, 1e-9);

    acb_mat_t WtW, one_m, Wd, blk, prod, one_b;
    acb_mat_init(WtW, m, m);
    acb_mat_init(one_m, m, m);
    acb_mat_zero(WtW);
    acb_mat_one(one_m);
    for (slong j = 0; j < nb; j++) {
        slong db = out->dim_L[j] * out->dim_E[j];
        acb_mat_init(Wd, m, db);
        acb_mat_conjugate_transpose(Wd, out->W_j[j]);   /* W_j^dag (m x db) */
        /* W_j W_j^dag = 1_block (isometry onto the block) */
        acb_mat_init(blk, db, db);
        acb_mat_init(one_b, db, db);
        acb_mat_mul(blk, out->W_j[j], Wd, prec);
        acb_mat_one(one_b);
        acb_mat_sub(blk, blk, one_b, prec);
        aic_mat_opnorm(e, blk, prec);
        if (!arb_le(e, tol)) {
            flint_printf("aic_idemp_wedderburn: W_%wd not an isometry "
                         "(||W_j W_j^dag - 1|| too large); fail loud (Rule 4)\n", j);
            flint_abort();
        }
        /* accumulate W_j^dag W_j into WtW */
        acb_mat_init(prod, m, m);
        acb_mat_mul(prod, Wd, out->W_j[j], prec);
        acb_mat_add(WtW, WtW, prod, prec);
        acb_mat_clear(prod);
        acb_mat_clear(one_b);
        acb_mat_clear(blk);
        acb_mat_clear(Wd);
    }
    acb_mat_sub(WtW, WtW, one_m, prec);
    aic_mat_opnorm(e, WtW, prec);
    if (!arb_le(e, tol)) {
        flint_printf("aic_idemp_wedderburn: stacked W not unitary "
                     "(||sum_j W_j^dag W_j - 1_M|| too large); fail loud (Rule 4)\n");
        flint_abort();
    }
    acb_mat_clear(one_m);
    acb_mat_clear(WtW);

    /* the correctness spec proper (FIX 2) — unitarity alone does NOT imply it. */
    wedd_certify_w_recon(out, Wk, dA, tol, prec);

    arb_clear(e);
    arb_clear(tol);
}
