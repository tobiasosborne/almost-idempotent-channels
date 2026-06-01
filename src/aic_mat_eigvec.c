/* aic_mat_eigvec.c — certified eigenvalue + INVARIANT-SUBSPACE decomposition of
 * a HERMITIAN matrix (bead aic-4td increment 2 step C2; design
 * docs/research/eigvec_certified_design.md §1.2/§1.6). Builds on the C1
 * eigenvalue layer (aic_mat_eig_multiple.c) and the shared densifier
 * (aic_mat_densify.c). The cluster->projector step (§1.4) lives in the sibling
 * aic_mat_cluster_proj.c (Rule 10, ~200 LOC/file).
 *
 * WHY THIS EXISTS. The C1 layer certifies eigenvalue BALLS but not the invariant
 * SUBSPACES the carrier split (lem_carrier, .tex:1724), certified Choi->Kraus, and
 * the Choi-Effros structure need. FLINT's acb_mat_eig_enclosure_rump returns, per
 * cluster, {lambda ball, n x k invariant-subspace enclosure X, k x k J with
 * spec(J) = the cluster eigenvalues and A X - X J certified to enclose 0}.
 *
 * DENSIFY-RUMP (FINDINGS §D5n RESOLVED, design §2). Rump's frozen-row partition
 * degenerates on a ROW-SPARSE invariant subspace and emits a NON-FINITE ([+/-inf])
 * enclosure. The fix: conjugate A' = U H U^dag by the dense unitary U
 * (aic_mat_dense_unitary), spreading every eigenvector across all n rows. Spectrum
 * is conjugation-invariant, so A''s lambda balls are H's, and the back-mapped
 * X_c = U^dag X'_c is a certified invariant subspace of H carrying the SAME J_c
 * (design §1.6(ii): A X_c = U^dag A' X'_c = U^dag X'_c J_c = X_c J_c). The residual
 * ||H X_c - X_c J_c|| recomputed on the ORIGINAL H is the certificate (~1e-31 at
 * prec=128, design §2.2).
 *
 * SOUNDNESS (design §1.6). (i) a FINITE Rump output is a Krawczyk certificate —
 * ASSERT finiteness (else fail loud "unresolved"). (ii) the conjugation is
 * spectrum/subspace-preserving (asserted ||U U^dag - I|| tiny). (iii) cross-cluster
 * lambda-ball DISJOINTNESS (!acb_overlaps, asserted, else fail loud "overlap")
 * proves distinct eigenvalues, so Sum k_c = n and the subspaces are mutually
 * orthogonal. gap_thr only chooses which approx evals feed ONE Rump call;
 * (i)+(iii) enforce soundness regardless (§1.5) — never silently wrong.
 */
#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic_mat_internal.h"

void aic_mat_eig_hermitian_subspaces(aic_mat_eigcluster **clusters,
                                     slong *n_clusters, const acb_mat_t H,
                                     double gap_thr, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(n >= 1 && acb_mat_ncols(H) == n &&
           "aic_mat_eig_hermitian_subspaces: H must be square");
    aic_mat_int_assert_hermitian(H, prec);

    /* DENSIFY A' = U H U^dag (design §1.2 step 1). */
    acb_mat_t U, Ud, A1, T;
    acb_mat_init(U, n, n);
    acb_mat_init(Ud, n, n);
    acb_mat_init(A1, n, n);
    acb_mat_init(T, n, n);
    aic_mat_dense_unitary(U, n, prec);
    acb_mat_conjugate_transpose(Ud, U);
    aic_mat_int_assert_densify_unitary(U, Ud, n, prec);
    acb_mat_mul(T, U, H, prec);
    acb_mat_mul(A1, T, Ud, prec);

    /* SEED: double-path zheev on midpoint(A') -> ascending evals + orthonormal
     * eigenvectors (seed only, NOT trusted; Rump certifies). */
    double _Complex *A1a = malloc((size_t) (n * n) * sizeof(double _Complex));
    double _Complex *Vd = malloc((size_t) (n * n) * sizeof(double _Complex));
    double *ev = malloc((size_t) n * sizeof(double));
    assert(A1a && Vd && ev && "aic_mat_eig_hermitian_subspaces: OOM");
    aic_latd_from_acb_mat(A1a, A1);
    aic_latd_eig_hermitian(ev, Vd, A1a, n);

    /* CLUSTER: gap-partition ev[]. Default gap_thr (design §1.5): a relative gap
     * tied to the spectral spread, with an absolute floor, so the zero cluster
     * separates from the smallest nonzero eigenvalue. */
    if (gap_thr <= 0.0) {
        double spread = ev[n - 1] - ev[0];
        if (spread < 1.0) spread = 1.0;
        gap_thr = 1e-6 * spread;
    }

    /* First pass: count clusters (a new cluster starts when the gap exceeds
     * gap_thr). */
    slong nc = (n >= 1) ? 1 : 0;
    for (slong j = 1; j < n; j++)
        if (ev[j] - ev[j - 1] > gap_thr) nc++;

    aic_mat_eigcluster *cl = malloc((size_t) nc * sizeof(aic_mat_eigcluster));
    assert(cl && "aic_mat_eig_hermitian_subspaces: OOM");

    /* Second pass: per cluster, Rump on A', back-map (aic_mat_int_certify_cluster). */
    slong c = 0, s0 = 0, ksum = 0;
    for (slong j = 1; j <= n; j++) {
        if (j == n || ev[j] - ev[j - 1] > gap_thr) {
            slong k = j - s0;                /* cluster [s0, j) */
            aic_mat_int_certify_cluster(&cl[c], A1, Ud, Vd, ev, n, s0, k, c,
                                        prec);
            ksum += k;
            s0 = j;
            c++;
        }
    }
    assert(c == nc);

    free(ev);
    free(Vd);
    free(A1a);
    acb_mat_clear(T);
    acb_mat_clear(A1);
    acb_mat_clear(Ud);
    acb_mat_clear(U);

    /* (iii) DISJOINTNESS gate (design §1.6): all pairs, for rigor. DEFENCE IN
     * DEPTH (FINDINGS §D5n2): Rump self-isolates — when both enclosures are finite
     * the balls are already disjoint, so this gate is reached only when the finite
     * guard (i) is bypassed and non-finite [+/-inf] balls flow here (mutation-proven
     * load-bearing in test_eigvec S7). */
    for (slong a = 0; a < nc; a++)
        for (slong b = a + 1; b < nc; b++)
            if (acb_overlaps(cl[a].lambda, cl[b].lambda)) {
                fprintf(stderr,
                        "aic_mat_eig_hermitian_subspaces: clusters %ld and %ld "
                        "OVERLAP at prec=%ld (lambda balls not disjoint — the "
                        "subspace decomposition is not certified; raise prec); "
                        "FINDINGS §D5n, bead aic-4td\n",
                        (long) a, (long) b, (long) prec);
                abort();
            }
    assert(ksum == n && "aic_mat_eig_hermitian_subspaces: Sum k_c != n");

    *clusters = cl;
    *n_clusters = nc;
}

void aic_mat_eigcluster_free(aic_mat_eigcluster *clusters, slong n_clusters)
{
    if (clusters == NULL) return;
    for (slong c = 0; c < n_clusters; c++) {
        acb_clear(clusters[c].lambda);
        acb_mat_clear(clusters[c].X);
        acb_mat_clear(clusters[c].J);
    }
    free(clusters);
}
