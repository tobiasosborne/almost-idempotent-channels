/* aic_mat.h — Layer-0 matrix substrate (beads aic-9zh, "T-mat").
 *
 * The matrix type IS FLINT's `acb_mat_t` (CLAUDE.md Law 2: reuse, do not invent
 * a new struct). FLINT 3.0.1 already gives multiplication, transpose,
 * conjugate-transpose, trace, the certified Frobenius norm, and a certified
 * eigenvalue family; this module adds ONLY the free functions over `acb_mat_t`
 * that the paper's later modules (ucp, idemp_structure, opspace, ...) need and
 * that FLINT does not expose directly: operator norm, a Hermitian-eig wrapper
 * with the Hermiticity precondition asserted, singular values, the Kronecker
 * product, and partial traces.
 *
 * Number path / cross-check ladder (beads aic-5ty). This module is the arb/acb
 * path; every routine takes an explicit `prec`. Its own internal cross-check is
 * "acb@prec=53 vs acb@high-prec": the tests run it at 53 and at 256 and require
 * the two certified balls to agree within a 53-bit-appropriate tol
 * (tests/test_mat.c). (LAPACK was unavailable when this was first written; it is
 * now installed and the double path lives in the separate `latd` module, against
 * which arb is cross-checked double-vs-arb@53 — see include/aic_latd.h.)
 *
 * Certified vs approximate (the FLINT eig wrinkle). FLINT's certified eig path
 * is two-stage: `acb_mat_approx_eig_qr` produces *heuristic* float
 * approximations (no error guarantee), then `acb_mat_eig_simple` *proves*
 * isolating complex-interval enclosures `E` with eigenvector enclosures `R` and
 * `L = R^{-1}` such that `L A R = diag(E)` (acb_mat.rst lines 740-756). There is
 * NO Hermitian-specialised solver in FLINT 3.0.1 — only this general
 * (non-Hermitian) family — so we feed it the general matrix, then EXPLOIT
 * Hermiticity for checking: a Hermitian input must yield real eigenvalues, which
 * we assert (the imaginary part of each enclosure must contain 0; see
 * aic_mat_eig_hermitian). If `acb_mat_eig_simple` cannot isolate the spectrum at
 * the given prec it returns 0 and we abort (CLAUDE.md Rule 4: fail loud), rather
 * than fall back to a cluster enclosure — the multiple-eigenvalue fallback is a
 * separate Law-4 audition.
 *
 * Paper provenance. These are generic linear-algebra primitives, not a single
 * paper result, but they are the substrate for results that ARE cited:
 *   - operator norm / cb-norm ampliation sup, approximate_algebras.tex:349;
 *   - partial trace in the decode map rho_j = Tr_{E_j}(W_j rho W_j^dag),
 *     approximate_algebras.tex:277 and :288, over the tensor factoring
 *     M -> L_j (x) E_j of approximate_algebras.tex:263-264.
 * The partial-trace index convention chosen here (below) is load-bearing for
 * those UCP constructions, so it is pinned now.
 *
 * Tensor / partial-trace index convention (LOAD-BEARING — pinned here).
 *   A matrix on the tensor-factored space C^a (x) C^b is stored as an
 *   (a*b) x (a*b) `acb_mat_t` with the LEFT factor MAJOR (row-major Kronecker,
 *   matching acb_mat_kronecker_product below and NumPy/textbook `kron`):
 *
 *       row index  I = i*b + j   with i in [0,a), j in [0,b)   (i = left, j = right)
 *       col index  J = k*b + l   with k in [0,a), l in [0,b)
 *       M[I,J] = <i,j| M |k,l>.
 *
 *   "Trace out factor 2 (the RIGHT/C^b factor)" sums the diagonal over j=l,
 *   leaving an a x a operator on C^a (the LEFT factor):
 *       (Tr_2 M)[i,k] = sum_{j} M[i*b + j, k*b + j].
 *   "Trace out factor 1 (the LEFT/C^a factor)" sums over i=k, leaving a b x b
 *   operator on C^b (the RIGHT factor):
 *       (Tr_1 M)[j,l] = sum_{i} M[i*b + j, i*b + l].
 *   So Tr_2 keeps the FIRST factor and Tr_1 keeps the SECOND. With this
 *   convention Tr_2(A (x) B) = Tr(B) * A and Tr_1(A (x) B) = Tr(A) * B for the
 *   Kronecker product below; tests/test_mat.c asserts both identities.
 */
#ifndef AIC_MAT_H
#define AIC_MAT_H

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Frobenius norm sqrt(sum |A_ij|^2) as a certified arb ball. Thin wrapper over
 * FLINT's acb_mat_frobenius_norm; provided for API symmetry with aic_mat_opnorm
 * and so callers depend on aic_mat, not on the FLINT spelling. `out` must be an
 * initialised arb_t; A may be rectangular. */
void aic_mat_frobenius_norm(arb_t out, const acb_mat_t A, slong prec);

/* Largest eigenvalue of a Hermitian H, as a certified arb ball, ROBUST TO
 * DEGENERACY. Asserts Hermiticity (fail loud). Route: FLINT's
 * acb_mat_eig_global_enclosure (acb_mat.rst:653-694) gives a rigorous radius eps
 * with every eigenvalue within eps of some approximate eigenvalue; the largest
 * true eigenvalue then lies in [max_k Re(E_k) - eps, max_k Re(E_k) + eps], which
 * is returned as the ball `out`. Unlike aic_mat_eig_hermitian this does NOT
 * require a simple spectrum, so it works on the identity (lambda_max = 1) and
 * any repeated-eigenvalue Gram matrix. `out` must be initialised. */
void aic_mat_herm_max_eig(arb_t out, const acb_mat_t H, slong prec);

/* Operator norm = largest singular value = sqrt(largest eigenvalue of A^dag A),
 * as a certified arb ball. A may be rectangular (m x n); A^dag A is the n x n
 * Hermitian PSD Gram matrix. Route: form the smaller Gram matrix, take its
 * largest Hermitian eigenvalue (aic_mat_herm_max_eig, degeneracy-robust), sqrt
 * it. `out` must be initialised. */
void aic_mat_opnorm(arb_t out, const acb_mat_t A, slong prec);

/* Hermitian eigendecomposition of a Hermitian H (n x n). Asserts Hermiticity of
 * H up to a prec-appropriate tol (CLAUDE.md Rule 4, fail loud) BEFORE solving.
 * Outputs:
 *   evals : caller-allocated arb_t[n], the n real eigenvalues, ASCENDING. The
 *           certified acb enclosures FLINT returns are asserted to have an
 *           imaginary part containing 0 (Hermitian => real spectrum), and the
 *           real part is copied out.
 *   evecs : if non-NULL, an n x n acb_mat_t (caller-initialised) whose COLUMNS
 *           are the certified right-eigenvector enclosures, column k <-> evals[k]
 *           (same ascending order).
 * Reconstruction V diag(evals) V^dag ~ H holds within tol; tests assert it.
 *
 * REQUIRES A SIMPLE (non-degenerate) SPECTRUM. The underlying acb_mat_eig_simple
 * proves n ISOLATED eigenvalues; on a repeated eigenvalue it returns 0 and this
 * routine aborts (CLAUDE.md Rule 4). Degenerate-spectrum eigenvectors (via
 * acb_mat_eig_multiple's invariant-subspace enclosures) are a separate Law-4
 * audition, not implemented here. For the largest eigenvalue alone (no
 * eigenvectors) use aic_mat_herm_max_eig, which IS degeneracy-robust. */
void aic_mat_eig_hermitian(arb_ptr evals, acb_mat_t evecs,
                           const acb_mat_t H, slong prec);

/* Certified eigenvalue enclosures for a HERMITIAN matrix with a possibly
 * DEGENERATE spectrum. Degeneracy-robust counterpart to aic_mat_eig_hermitian
 * (which REQUIRES a simple spectrum and aborts on any repeat). Route (beads
 * aic-w4o.1 / aic-4td, FINDINGS §D5/§D5n): acb_mat_approx_eig_qr seed ->
 * acb_mat_eig_multiple (Rump cluster enclosures, degeneracy-native). E
 * (caller-allocated acb_ptr of length n = nrows(H)) receives n certified
 * eigenvalue balls, GROUPED: a run of k identical balls is a certified k-cluster
 * (could not be split at this prec but is isolated from the other n-k
 * eigenvalues). DENSIFY-RETRY (FINDINGS §D5n RESOLVED): on acb_mat_eig_multiple
 * returning 0 (FLINT Rump's frozen-row partition failing on a ROW-SPARSE
 * invariant subspace — e.g. a C^5 {2,3} projector or diag(2,2,0,0)) it RETRIES on
 * the densified A' = U H U† (U = aic_mat_dense_unitary; spectrum conjugation-
 * invariant, so the balls are those of H written straight to E), asserting U
 * certified-unitary first. ABORTS (Rule 4) only if the DENSIFIED retry ALSO
 * returns 0 (a genuine near-degeneracy unresolvable at this prec -> raise prec).
 * Asserts each ball's imaginary part encloses 0 (Hermitian => real spectrum). */
void aic_mat_eig_hermitian_multiple(acb_ptr E, const acb_mat_t H, slong prec);

/* Certified numerical rank of a Hermitian matrix H relative to threshold thr
 * (an arb ball): the count of eigenvalues certified ABOVE thr. Every eigenvalue
 * ball must be certified-above (arb_gt) or certified-below (arb_lt) thr; if any
 * ball STRADDLES thr the rank is unresolved at this prec and the routine ABORTS
 * (Rule 4) with a clear message (raise prec / move thr out of the cluster).
 * Degeneracy-robust: delegates to aic_mat_eig_hermitian_multiple. This is the
 * certified counterpart to the LAPACK double-path count in
 * aic_ucp_carrier_rank_latd (bead aic-w4o.1). */
slong aic_mat_certified_rank(const acb_mat_t H, const arb_t thr, slong prec);

/* One certified eigen-cluster of a Hermitian matrix: a real eigenvalue ball
 * lambda, an n x k certified invariant-subspace enclosure X (A X subset of
 * span(X)), the cluster multiplicity k, and the k x k interval matrix J from
 * FLINT's Rump enclosure (acb_mat_eig_enclosure_rump). X is NOT orthonormal
 * (Rump applies no normalization); use aic_mat_cluster_projector for the
 * orthogonal range projector. J is stored because the residual certificate
 * ||H X_c - X_c J_c|| on the ORIGINAL H is the soundness proof
 * (docs/research/eigvec_certified_design.md §1.6(ii): the SAME J_c works for H
 * via X_c = U^dag X'_c, since A' = U H U^dag and the eigenvalues of J_c are the
 * cluster's). Storing J lets a caller recompute that residual on H directly
 * without re-deriving J from a Rayleigh quotient. */
typedef struct {
    acb_t lambda;     /* certified eigenvalue ball (imag part asserted to contain 0) */
    acb_mat_t X;      /* n x k certified invariant-subspace enclosure (A X subset span(X)) */
    acb_mat_t J;      /* k x k Rump interval matrix; spec(J) = the cluster eigenvalues */
    slong k;          /* multiplicity (cluster size) */
} aic_mat_eigcluster;

/* Certified eigenvalue + INVARIANT-SUBSPACE decomposition of a Hermitian H,
 * degeneracy-robust via dense-unitary densification + per-cluster Rump
 * (acb_mat_eig_enclosure_rump). Resolves the §D5n two-clusters-each->=2 wall
 * (FINDINGS §D5n RESOLVED; design §1.2/§2). ALGORITHM: densify A' = U H U^dag
 * (U = aic_mat_dense_unitary; assert ||U U^dag - I|| certified tiny); double-path
 * zheev seed on midpoint(A') for ascending evals + orthonormal eigenvectors;
 * gap-cluster the evals (gap_thr; pass <= 0 for the default, §1.5); per cluster
 * run Rump on A' with the cluster's seed columns, back-map X_c = U^dag X'_c.
 * *clusters is allocated HERE to *n_clusters entries (each lambda/X/J init'd),
 * sorted ascending by eigenvalue; free with aic_mat_eigcluster_free. ASSERTS
 * (Rule 4): Hermiticity; densifier certified unitary; every per-cluster Rump
 * enclosure FINITE (else "clusters unresolved at prec" — raise prec); imag(lambda)
 * contains 0; the SELF-CERTIFYING residual ||H X_c - X_c J_c||_F on the ORIGINAL H
 * small per cluster (the honest certificate, design §1.6, FINDINGS §D5n — else
 * "NOT a certified invariant subspace"); cross-cluster lambda balls DISJOINT
 * (!acb_overlaps, else "clusters overlap at prec"); Sum k_c == n. */
void aic_mat_eig_hermitian_subspaces(aic_mat_eigcluster **clusters,
                                     slong *n_clusters, const acb_mat_t H,
                                     double gap_thr, slong prec);

/* Free a cluster array returned by aic_mat_eig_hermitian_subspaces (clears each
 * lambda/X/J then frees the array). No-op on (NULL, 0). */
void aic_mat_eigcluster_free(aic_mat_eigcluster *clusters, slong n_clusters);

/* Orthogonal projector Pi = X (X^dag X)^-1 X^dag onto a certified invariant
 * subspace (the cluster's range), design §1.4. `out` must be initialised
 * n x n (n = nrows(cl->X)). ASSERTS (Rule 4) X^dag X is certified-invertible
 * (acb_mat_inv != 0) — fail loud if the subspace enclosure is rank-deficient at
 * this prec. The result is certified self-adjoint and idempotent to the ball
 * radius (measured ||Pi^2-Pi||, ||Pi-Pi^dag|| ~ 1e-31 at prec=128, design §2.2). */
void aic_mat_cluster_projector(acb_mat_t out, const aic_mat_eigcluster *cl,
                               slong prec);

/* Singular values of a general m x n matrix A, as sqrt of the eigenvalues of the
 * smaller Gram matrix (A^dag A if n<=m, else A A^dag), in DESCENDING order.
 *   svals : caller-allocated arb_t[min(m,n)], the singular values.
 * Goes through aic_mat_eig_hermitian, so it REQUIRES DISTINCT singular values
 * (simple Gram spectrum); a repeated singular value aborts. Full U, Sigma, V SVD
 * and the degenerate-spectrum case are later beads; only the distinct values are
 * produced here. */
void aic_mat_singular_values(arb_ptr svals, const acb_mat_t A, slong prec);

/* Kronecker (tensor) product: out = A (x) B, LEFT factor major (row-major), so
 * for A (m x n) and B (p x q), out is (m*p) x (n*q) with
 *   out[i*p + r, j*q + s] = A[i,j] * B[r,s].
 * `out` must be initialised to the (m*p) x (n*q) shape (asserted). */
void aic_mat_kronecker(acb_mat_t out, const acb_mat_t A, const acb_mat_t B,
                       slong prec);

/* Partial trace over the RIGHT factor (C^b). Input M is (a*b) x (a*b) on
 * C^a (x) C^b (index convention in the file docstring); `out` is the a x a
 * result on C^a. (Tr_2 M)[i,k] = sum_j M[i*b+j, k*b+j]. `out` must be
 * initialised a x a; a*b must equal nrows(M)=ncols(M) (asserted). */
void aic_mat_partial_trace_right(acb_mat_t out, const acb_mat_t M,
                                 slong a, slong b, slong prec);

/* Partial trace over the LEFT factor (C^a). `out` is the b x b result on C^b.
 * (Tr_1 M)[j,l] = sum_i M[i*b+j, i*b+l]. `out` must be initialised b x b;
 * a*b must equal nrows(M)=ncols(M) (asserted). */
void aic_mat_partial_trace_left(acb_mat_t out, const acb_mat_t M,
                                slong a, slong b, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_MAT_H */
