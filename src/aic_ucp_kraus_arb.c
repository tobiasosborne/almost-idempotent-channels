/* aic_ucp_kraus_arb.c — CERTIFIED arb-path Choi->Kraus extraction (bead aic-4td
 * increment 2 step D; design docs/research/eigvec_certified_design.md §3.2). The
 * certified counterpart of the LAPACK double-path aic_ucp_choi_to_kraus_latd
 * (src/aic_ucp_latd.c:1-40): SAME Convention-A conjugate reshape, SAME QETLAB
 * keep_thr gate, SAME PSD-cone _tol variant (FINDINGS §C14), but on the certified
 * arb path via the degeneracy-robust eigensolver of step C1/C2
 * (aic_mat_eig_hermitian_subspaces, src/aic_mat_eigvec.c).
 *
 * THE EXTRACTION (Convention A; src/aic_ucp_latd.c:9-24, include/aic/aic_ucp.h
 * :23-37). C >= 0 is PSD: C = sum_a lambda_a v_a v_a^dag, {v_a} ORTHONORMAL.
 * Matching the header Convention A  C[i*h+a, j*h+b] = sum_x conj(K_x[i,a]) K_x[j,b]
 * term-by-term forces the CONJUGATE reshape
 *   K_a[i, c] = sqrt(lambda_a) * conj(v_a[i*dim_H + c]),
 * so kraus_to_choi(extracted) == C. Keep eigenpairs with lambda > keep_thr =
 * (dim_K*dim_H) * 2^-52 * ||C||_F (the _latd DBL_EPSILON=2^-52 QETLAB gate).
 *
 * ORTHONORMALISATION (design §3.2 R2; the S4 gate). Rump's per-cluster X is NOT
 * orthonormal (FLINT applies no normalisation); a RAW reshape would not rebuild C.
 * Per cluster we LOEWDIN-orthonormalise  V = X (X^dag X)^{-1/2}  (V^dag V = I_k,
 * range(V) = range(X)) via aic_funcalc_xpow(alpha=-1/2) at x0 = ||X||_op^2 =
 * lambda_max(G) (G = X^dag X), which makes ||G/x0 - I||_op = 1 - l_min/l_max < 1
 * so the xpow domain (tex:511) always holds for PD G. x0 comes from aic_mat_opnorm
 * (Weyl-Hermitianizes the interval Gram, FINDINGS §C5) — NOT herm_max_eig(G),
 * which aborts on the raw interval Gram X^dag X. k=1: V = X / sqrt(X^dag X).
 *
 * CLUSTER-EIGENVALUE ROUTE. Each kept cluster's SINGLE eigenvalue ball lambda_c is
 * scaled across all k orthonormal columns. EXACT for (i) a SIMPLE eigenvalue
 * (k=1) and (ii) a GENUINE degeneracy ((+)B(L) block multiplicities, all k equal).
 * Gap-clustering gives a distinct nonzero eigenvalue its OWN k=1 cluster, so the
 * round-trip rebuilds C to WITHIN the ball (S4 is the acceptance gate). Were two
 * DISTINCT nonzero eigenvalues ever lumped, lambda_c-for-all would be off by the
 * cluster width; S4 catches it (RED) and we fail loud (no project Choi triggers
 * it — carrier/block spectra are genuine degeneracies, design §1.5/§2).
 *
 * CERTIFIED THREE-WAY DECISION (Rule 4; the aic_mat_certified_rank pattern). The
 * _tol PSD-cone variant (FINDINGS §C14): the almost-idempotent CP-ization Delta
 * has a per-block Choi PSD only to O(eta^2). Per cluster, re = Re(lambda_c) vs the
 * keep_thr / -neg_tol BALLS with rigorous arb_gt/arb_lt:
 *   re > keep_thr            -> KEEP (emit k Kraus);
 *   re in (-neg_tol, keep_thr] -> DROP (cone-defect/noise; clip negative mid into
 *                                the certified mass — the PSD-cone projection);
 *   re < -neg_tol            -> FAIL LOUD (genuine non-CP);
 *   re STRADDLES either bound -> FAIL LOUD "STRADDLE" (rank unresolved; raise prec)
 *                                — NDEBUG-immune fprintf+abort.
 * The strict aic_ucp_choi_to_kraus_arb delegates with neg_tol = keep_thr.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_ucp_kraus_arb_internal.h"

void aic_ucp_choi_to_kraus_arb_tol(aic_ucp_kraus *phi, const acb_mat_t C,
                                   slong dim_K, slong dim_H, double neg_tol,
                                   double *clipped_neg_out, slong prec)
{
    slong n = dim_K * dim_H;
    assert(acb_mat_nrows(C) == n && acb_mat_ncols(C) == n &&
           "choi_to_kraus_arb: C must be (dim_K*dim_H) square");
    assert(neg_tol >= 0.0 && "choi_to_kraus_arb_tol: neg_tol must be >= 0");

    arb_t thr, negtol;
    arb_init(thr);
    arb_init(negtol);
    aic_ucp_kraus_arb_int_keep_thr(thr, C, n, prec);
    arb_set_d(negtol, -neg_tol);             /* the -neg_tol boundary ball */

    /* Certified eigen-cluster decomposition (degeneracy-robust, C1/C2). */
    aic_mat_eigcluster *cl;
    slong nc;
    aic_mat_eig_hermitian_subspaces(&cl, &nc, C, -1.0, prec);

    /* PASS 1: certified three-way decision per cluster -> kept rank r + clipped. */
    arb_t re;
    arb_init(re);
    slong r = 0;
    double clipped = 0.0;
    for (slong c = 0; c < nc; c++) {
        acb_get_real(re, cl[c].lambda);
        if (arb_gt(re, thr)) {
            r += cl[c].k;                     /* KEEP this cluster (k Kraus ops) */
        } else if (arb_lt(re, negtol)) {
            /* certified < -neg_tol: a GENUINELY non-PSD (not-CP) input. */
            char *re_s = arb_get_str(re, 12, 0);
            fprintf(stderr,
                    "aic_ucp_choi_to_kraus_arb_tol: eigenvalue cluster %ld = %s "
                    "certified < -neg_tol=%g; C is not PSD (not CP) beyond the "
                    "admitted cone tolerance — cannot extract Kraus (Rule 4: fail "
                    "loud); bead aic-4td, FINDINGS §C14\n",
                    (long) c, re_s ? re_s : "?", neg_tol);
            flint_free(re_s);
            assert(0 && "choi_to_kraus_arb_tol: cluster eigenvalue < -neg_tol (not CP)");
            abort();
        } else if (arb_lt(re, thr) && arb_gt(re, negtol)) {
            /* certified in (-neg_tol, keep_thr]: DROP (numerically-zero tail OR an
             * approximately-CP-cone DEFECT). If the midpoint is negative,
             * accumulate k*|mid| as the certified cone-defect mass (the PSD-cone
             * projection: keep only the nonnegative part). */
            double mid = arf_get_d(arb_midref(re), ARF_RND_NEAR);
            if (mid < 0.0) clipped += -mid * (double) cl[c].k;
        } else {
            /* The ball overlaps keep_thr or -neg_tol without being certified
             * cleanly above/below: the rank (kept-vs-dropped) is UNRESOLVED at
             * this prec. Fail loud (Rule 4), like aic_mat_certified_rank's
             * straddle (NDEBUG-immune fprintf+abort). */
            char *re_s = arb_get_str(re, 12, 0);
            char *thr_s = arb_get_str(thr, 12, 0);
            fprintf(stderr,
                    "aic_ucp_choi_to_kraus_arb_tol: eigenvalue cluster %ld = %s "
                    "STRADDLES the keep threshold %s (or -neg_tol=%g) — kept-vs-"
                    "dropped rank unresolved at prec=%ld; raise prec or move the "
                    "threshold out of the cluster; bead aic-4td\n",
                    (long) c, re_s ? re_s : "?", thr_s ? thr_s : "?", neg_tol,
                    (long) prec);
            flint_free(re_s);
            flint_free(thr_s);
            assert(0 && "choi_to_kraus_arb_tol: cluster eigenvalue straddles keep_thr");
            abort();
        }
    }
    assert(r >= 1 && "choi_to_kraus_arb: zero rank (the zero map is not UCP)");

    aic_ucp_kraus_init(phi, dim_K, dim_H, r);

    /* PASS 2: per KEPT cluster, orthonormalise X then emit sqrt(lambda)-scaled
     * conjugate-reshape Kraus ops (one per orthonormal column). */
    acb_mat_t V;
    arb_t scale;
    arb_init(scale);
    slong out_a = 0;
    for (slong c = 0; c < nc; c++) {
        acb_get_real(re, cl[c].lambda);
        if (!arb_gt(re, thr)) continue;       /* only the kept (certified-above) */
        slong k = cl[c].k;
        acb_mat_init(V, n, k);
        aic_ucp_kraus_arb_int_loewdin(V, cl[c].X, prec);
        arb_set(scale, re);
        arb_sqrt(scale, scale, prec);         /* sqrt(lambda_c) (re > thr > 0) */
        for (slong m = 0; m < k; m++) {
            /* K_{out_a}[i, ch] = sqrt(lambda_c) * conj(V[i*dim_H + ch, m]). */
            for (slong i = 0; i < dim_K; i++)
                for (slong ch = 0; ch < dim_H; ch++) {
                    acb_t v;
                    acb_init(v);
                    acb_conj(v, acb_mat_entry(V, i * dim_H + ch, m));
                    acb_mul_arb(v, v, scale, prec);
                    acb_set(acb_mat_entry(phi->K[out_a], i, ch), v);
                    acb_clear(v);
                }
            out_a++;
        }
        acb_mat_clear(V);
    }
    assert(out_a == r && "choi_to_kraus_arb: emitted count != kept rank");

    if (clipped_neg_out) *clipped_neg_out = clipped;

    arb_clear(scale);
    arb_clear(re);
    aic_mat_eigcluster_free(cl, nc);
    arb_clear(negtol);
    arb_clear(thr);
}

void aic_ucp_choi_to_kraus_arb(aic_ucp_kraus *phi, const acb_mat_t C,
                               slong dim_K, slong dim_H, slong prec)
{
    /* STRICT default: delegate with neg_tol = keep_thr (the historical
     * (dim_K*dim_H)*2^-52*||C||_F gate). For lambda <= -keep_thr the _tol routine
     * fails loud EXACTLY as the strict QETLAB gate; clipped mass is discarded. */
    slong n = dim_K * dim_H;
    arb_t thr;
    arb_init(thr);
    aic_ucp_kraus_arb_int_keep_thr(thr, C, n, prec);
    double neg_tol = arf_get_d(arb_midref(thr), ARF_RND_UP) +
                     mag_get_d(arb_radref(thr));
    arb_clear(thr);
    aic_ucp_choi_to_kraus_arb_tol(phi, C, dim_K, dim_H, neg_tol, NULL, prec);
}
