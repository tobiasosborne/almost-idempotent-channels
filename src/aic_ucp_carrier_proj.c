/* aic_ucp_carrier_proj.c — the certified RANGE (support) projector Pi_M of the
 * carrier Q = sum_a K_a K_a^dag (lem_carrier, .tex:1724; bead aic-4td increment
 * 2 step E; design docs/research/eigvec_certified_design.md §3.3). Split out of
 * aic_ucp_carrier.c (Rule 10, ~200 LOC/file).
 *
 * PAPER RESULT vs CONSTRUCTIVE ROUTE.
 *   lem_carrier (.tex:1724): the carrier M is the support (range) of the
 *   Hermitian PSD operator Q = Tr_F(V V^dag) = sum_a K_a K_a^dag (see the sibling
 *   aic_ucp_carrier.c docstring for the derivation). The paper's proof is an
 *   existence statement via the full-rank-state support of Phi^*(rho_0); the
 *   CONSTRUCTIVE finite-dim route is the orthogonal projector onto range(Q),
 *   built from the certified eigendecomposition of Q.
 *
 *   Pi_M = sum over the eigen-clusters of Q whose eigenvalue is certified ABOVE
 *          thr = dim_K * 2^-52 * ||Q||_F of the cluster projector
 *          Pi_c = X_c (X_c^dag X_c)^-1 X_c^dag (aic_mat_cluster_projector).
 *
 * WHY THIS THRESHOLD. thr is the SAME arb ball as aic_ucp_carrier_rank (and the
 * double-path aic_ucp_carrier_rank_latd: dim_K * DBL_EPSILON * ||Q||_F), so the
 * count of clusters summed here equals the certified rank exactly — trace(Pi_M)
 * == aic_ucp_carrier_rank(Q) == aic_ucp_carrier_rank_latd(Q) (test_eigvec S6).
 *
 * DEGENERACY-ROBUST EIG. Q is rank-deficient with a large zero cluster and
 * possibly-degenerate nonzero clusters, so the certified SIMPLE-spectrum path
 * cannot be used; aic_mat_eig_hermitian_subspaces (densify + per-cluster Rump,
 * FINDINGS §D5n RESOLVED) is the degeneracy-native route. Q Hermitian/PSD =>
 * left = right invariant subspace, so Rump's right-subspace X IS the range basis
 * (no left/right hazard, FINDINGS §C4).
 *
 * FAIL LOUD (Rule 4). A cluster eigenvalue ball that STRADDLES thr leaves the
 * in-range/in-kernel decision UNDECIDED, so range(Q) is not certified at this
 * prec -> abort (NDEBUG-immune fprintf+abort, EXACTLY mirroring
 * aic_mat_certified_rank's straddle so the projector and the rank agree on the
 * boundary). PRECISION FLOOR: a rank-deficient Q's zero cluster has a certified
 * radius ~2e-14 at prec=53 that straddles thr ~1e-15 (the same densify+Rump
 * conditioning as the certified Choi->Kraus); the zero ball drops cleanly below
 * thr only at prec >= 64 (FINDINGS §D7). Run this at prec=128 for headroom.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

void aic_ucp_carrier_projector(acb_mat_t out, const acb_mat_t Q, slong dim_K,
                               slong prec)
{
    assert(acb_mat_nrows(Q) == dim_K && acb_mat_ncols(Q) == dim_K &&
           acb_mat_nrows(out) == dim_K && acb_mat_ncols(out) == dim_K &&
           "aic_ucp_carrier_projector: Q and out must be dim_K x dim_K");

    /* thr = dim_K * 2^-52 * ||Q||_F (the SAME arb ball as aic_ucp_carrier_rank,
     * so trace(Pi_M) == aic_ucp_carrier_rank(Q)). */
    arb_t fro, thr;
    arb_init(fro);
    arb_init(thr);
    aic_mat_frobenius_norm(fro, Q, prec);     /* certified ||Q||_F */
    arb_mul_si(thr, fro, dim_K, prec);        /* dim_K * ||Q||_F */
    arb_mul_2exp_si(thr, thr, -52);           /* * 2^-52 = DBL_EPSILON */

    /* Certified eigen-clusters of the Hermitian PSD carrier Q (degeneracy-robust,
     * densify + per-cluster Rump; design §1.2). */
    aic_mat_eigcluster *cl;
    slong nc;
    aic_mat_eig_hermitian_subspaces(&cl, &nc, Q, -1.0, prec);

    acb_mat_zero(out);
    arb_t re;
    arb_init(re);
    acb_mat_t Pi;
    acb_mat_init(Pi, dim_K, dim_K);
    for (slong c = 0; c < nc; c++) {
        acb_get_real(re, cl[c].lambda);
        if (arb_gt(re, thr)) {
            /* cluster certified ABOVE thr: in range(Q) -> add its projector. */
            aic_mat_cluster_projector(Pi, &cl[c], prec);
            acb_mat_add(out, out, Pi, prec);
        } else if (arb_lt(re, thr)) {
            /* cluster certified BELOW thr: in the kernel -> skip. */
        } else {
            /* The cluster eigenvalue ball overlaps thr: in-range/in-kernel is
             * UNDECIDED, so range(Q) is not certified at this prec. Fail loud
             * (Rule 4), NDEBUG-immune fprintf+abort EXACTLY as
             * aic_mat_certified_rank (the two must agree on the straddle). */
            char *re_s = arb_get_str(re, 12, 0);
            char *thr_s = arb_get_str(thr, 12, 0);
            fprintf(stderr,
                    "aic_ucp_carrier_projector: carrier eigenvalue cluster %ld = "
                    "%s STRADDLES the threshold %s (carrier range unresolved at "
                    "prec=%ld); raise prec (the rank-deficient zero cluster needs "
                    "prec >= 64, FINDINGS §D7); bead aic-4td\n",
                    (long) c, re_s ? re_s : "?", thr_s ? thr_s : "?",
                    (long) prec);
            flint_free(re_s);
            flint_free(thr_s);
            assert(0 && "aic_ucp_carrier_projector: cluster ball straddles thr");
            abort();
        }
    }

    acb_mat_clear(Pi);
    arb_clear(re);
    aic_mat_eigcluster_free(cl, nc);
    arb_clear(thr);
    arb_clear(fro);
}
