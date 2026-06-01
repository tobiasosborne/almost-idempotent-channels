/* aic_mat_cluster_proj.c — the certified orthogonal RANGE projector of a
 * certified eigen-cluster (bead aic-4td increment 2 step C2; design
 * docs/research/eigvec_certified_design.md §1.4). Split out of aic_mat_eigvec.c
 * (Rule 10, ~200 LOC/file).
 *
 * Rump's per-cluster X (from aic_mat_eig_hermitian_subspaces) is full-column-rank
 * but NOT orthonormal (FLINT applies no normalization). The orthogonal projector
 * onto its column space — the certified invariant subspace — is the standard
 * oblique-free Gram formula
 *
 *     Pi = X (X^dag X)^-1 X^dag        [G = X^dag X is k x k, certified-invertible]
 *
 * Here X (and the cluster's H) is Hermitian/PSD, so left = right singular vectors
 * and Rump's right-invariant-subspace X IS the range basis directly — no
 * left/right hazard (FINDINGS §C4; the §C4 caveat bites only for a NON-Hermitian
 * oblique idempotent, which the Hermiticity precondition of the subspace routine
 * rules out). G is certified well-conditioned because X's columns are a near-
 * orthonormal seed mapped through a unitary; acb_mat_inv succeeds and
 * ||Pi^2-Pi||, ||Pi-Pi^dag|| measure ~1e-31 at prec=128 (design §2.2). A
 * rank-deficient enclosure (acb_mat_inv == 0) fails loud (Rule 4).
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>

#include "aic/aic_mat.h"

void aic_mat_cluster_projector(acb_mat_t out, const aic_mat_eigcluster *cl,
                               slong prec)
{
    slong n = acb_mat_nrows(cl->X), k = cl->k;
    assert(acb_mat_ncols(cl->X) == k && acb_mat_nrows(out) == n &&
           acb_mat_ncols(out) == n &&
           "aic_mat_cluster_projector: shape mismatch");

    /* Pi = X (X^dag X)^-1 X^dag (design §1.4). */
    acb_mat_t Xd, G, Gi, Tmp;
    acb_mat_init(Xd, k, n);
    acb_mat_init(G, k, k);
    acb_mat_init(Gi, k, k);
    acb_mat_init(Tmp, n, k);
    acb_mat_conjugate_transpose(Xd, cl->X);
    acb_mat_mul(G, Xd, cl->X, prec);
    if (!acb_mat_inv(Gi, G, prec)) {
        fprintf(stderr,
                "aic_mat_cluster_projector: X^dag X NOT certified-invertible at "
                "prec=%ld (rank-deficient invariant-subspace enclosure for the "
                "k=%ld cluster; raise prec); bead aic-4td\n",
                (long) prec, (long) k);
        abort();
    }
    acb_mat_mul(Tmp, cl->X, Gi, prec);
    acb_mat_mul(out, Tmp, Xd, prec);

    acb_mat_clear(Tmp);
    acb_mat_clear(Gi);
    acb_mat_clear(G);
    acb_mat_clear(Xd);
}
