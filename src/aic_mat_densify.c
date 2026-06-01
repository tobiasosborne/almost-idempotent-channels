/* aic_mat_densify.c — the dense-unitary DENSIFIER shared across the certified
 * eig layer (bead aic-4td increment 2; design docs/research/eigvec_certified_design.md
 * §1.3). Split out of aic_mat_eig_multiple.c (Rule 10, ~200 LOC/file) because it
 * is REUSED by both the eigenvalue retry (aic_mat_eig_multiple.c) and the
 * certified invariant-subspace layer (aic_mat_eigvec.c, step C2).
 *
 * WHY IT EXISTS (FINDINGS §D5n RESOLVED). FLINT's certified eig family
 * (acb_mat_eig_multiple / acb_mat_eig_enclosure_rump) returns 0 / a non-finite
 * enclosure on a HERMITIAN input whose degenerate cluster's invariant subspace
 * is ROW-SPARSE (supported on few coordinates — e.g. an axis-aligned or
 * disjoint-Givens projector leaves whole rows ~0 across the cluster's columns).
 * The cause, pinned by measured probes (design §2), is Rump's frozen-row
 * partition (partition_X_sorted) having no well-conditioned k-row frozen set —
 * NOT seed conditioning (the original §D5n hypothesis, proven FALSE) and NOT a
 * precision wall (it persists at prec=256/1024). Conjugating A' = U A U† by a
 * DENSE unitary U spreads every eigenvector across all n rows, so the partition
 * is well-conditioned on A'. The spectrum is conjugation-invariant, so the
 * eigenvalue balls of A' are the eigenvalue balls of A.
 *
 * WHY THIS U (Law 3 "Haar diagonal -> explicit finite element"). U is the product
 * over ALL planes (a,b), a<b, of the rational Givens rotation cos=3/5, sin=4/5 —
 * the exact rationals already used in tests/test_eigmult.c's fixtures. It is a
 * FIXED explicit element: reproducible, RNG-free, exact-rational, and unitary far
 * below the working precision (||U U†-I||_F certified ~1.3e-37 at n=4, ~2.6e-37
 * at n=5, prec=128 — measured). The caller forms A' = U A U† (two acb_mat_mul,
 * U† = acb_mat_conjugate_transpose) and ASSERTS ||U U†-I|| certified-tiny before
 * trusting the conjugation (a loose U would silently move the spectrum).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_mat_internal.h"

void aic_mat_dense_unitary(acb_mat_t U, slong n, slong prec)
{
    assert(acb_mat_nrows(U) == n && acb_mat_ncols(U) == n &&
           "aic_mat_dense_unitary: U must be n x n");
    acb_t c, s, ns;
    acb_init(c);
    acb_init(s);
    acb_init(ns);
    acb_set_si(c, 3);
    acb_div_si(c, c, 5, prec);          /* cos = 3/5 */
    acb_set_si(s, 4);
    acb_div_si(s, s, 5, prec);          /* sin = 4/5 */
    acb_neg(ns, s);

    acb_mat_one(U);
    acb_mat_t G, T;
    acb_mat_init(G, n, n);
    acb_mat_init(T, n, n);
    for (slong a = 0; a < n; a++) {
        for (slong b = a + 1; b < n; b++) {     /* every plane (a,b), a<b */
            acb_mat_one(G);
            acb_set(acb_mat_entry(G, a, a), c);
            acb_set(acb_mat_entry(G, b, b), c);
            acb_set(acb_mat_entry(G, a, b), ns);
            acb_set(acb_mat_entry(G, b, a), s);
            acb_mat_mul(T, U, G, prec);
            acb_mat_set(U, T);
        }
    }
    acb_mat_clear(T);
    acb_mat_clear(G);
    acb_clear(ns);
    acb_clear(s);
    acb_clear(c);
}
