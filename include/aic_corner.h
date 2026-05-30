/* aic_corner.h — the compression maps Co_{P,Q} and the corner subspaces S_{P,Q}
 * of section 7 "Subspaces associated with projections" (approximate_algebras.tex
 * :1052-1076; shard paper/shards/shard-D-projection-subspaces.md). Bead aic
 * "corner", Increment 1. Built ON TOP of the eps-C* algebra data model
 * (include/aic_ecstar.h): an A <= M_n with a Frobenius-orthonormal operator basis
 * {B_k}_{k=0..d-1} and the Choi-Effros STAR X*Y = Phi_tilde(XY).
 *
 * WHAT THIS REALIZES (.tex:1054-1066, verbatim in src/aic_corner_compress.c).
 *   To each pair of delta-projections (P,Q) in an eps-C* algebra A the paper
 *   assigns the corner subspace S_{P,Q} = {P X Q : X in A} (the EXACT-C* picture)
 *   via the COMPRESSION MAP
 *
 *       Co_{P,Q} = theta(L_P R_Q + R_Q L_P - 1)            (.tex:1057)
 *
 *   where L_P : X |-> P*X (left star-multiply by P) and R_Q : X |-> X*Q (right
 *   star-multiply by Q) are operators ON A, and theta(Y) = (1 + sgn(Y))/2 is the
 *   funcalc projection-snapper (prop_P, .tex:524-533). The symmetric combination
 *   1/2(L_P R_Q + R_Q L_P) is an O(delta+eps)-IDEMPOTENT operator on A, and
 *   theta(2M-1) snaps it to an EXACT idempotent (prop_P with X = 2M-1, so
 *   ||X^2-I|| = 4||M^2-M||, basin ||M^2-M|| < 1/4). Co_{P,Q} then satisfies
 *
 *       Co_{P,Q}^2 = Co_{P,Q},   Co_{P,Q}(X)^dag = Co_{Q,P}(X^dag)   (.tex:1061-1062)
 *       ||L_P R_Q - Co_{P,Q}||, ||R_Q L_P - Co_{P,Q}|| <= O(delta+eps)  (.tex:1064)
 *
 *   and S_{P,Q} = Img Co_{P,Q} = Ker(1 - Co_{P,Q}) is a closed subspace of A
 *   (.tex:1064). One-sided variants Co_{P,1} = theta(L_P) (R_1 = id, so the
 *   symmetric combination collapses to L_P) and Co_{1,Q} = theta(R_Q); Co_P :=
 *   Co_{P,P} (.tex:1064-1066). Degenerate dims (.tex:1066): dim S_P = 0 iff P ~ 0,
 *   dim S_P = dim_A iff P ~ I, dim S_P = 1 iff P is a "one-dimensional
 *   delta-projection".
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The paper's Co is an
 * operator on the (possibly infinite-dim) Banach algebra A, snapped by the
 * holomorphic functional calculus of prop_P. In finite dimensions A is a d-dim
 * space with an orthonormal basis {B_k}, so an operator on A IS a d x d matrix.
 * We build L_P, R_Q as d x d matrices in the {B_k} COORDINATE space, form
 * 2M-1 = L_P R_Q + R_Q L_P - I, and call aic_prop_P (theta(2M-1), the eig-free
 * Newton-Schulz sgn) to get the EXACT idempotent Co (a d x d matrix). The
 * paper's O(delta+eps) bound is the spec; tests/test_corner.c asserts it.
 *
 * THE STAR, NOT THE PLAIN PRODUCT (CLAUDE.md domain callout, LOAD-BEARING).
 *   Every "product in A" goes through aic_ecstar_star (X*Y = Phi_tilde(XY)), NOT
 *   acb_mat_mul. Plain operator multiplication leaves Img Phi_tilde. So
 *       (L_P)_{ij} = <B_i, P * B_j>_F,   (R_Q)_{ij} = <B_i, B_j * Q>_F,
 *   with P*B_j = aic_ecstar_star(A, P, B_j) and B_j*Q = aic_ecstar_star(A, B_j, Q).
 *   acb_mat_mul appears ONLY for forming the Frobenius inner products <,>_F and
 *   for the d x d superop matrix products L_P R_Q, R_Q L_P (those ARE ordinary
 *   matrix products of the coordinate matrices, not products in A).
 *
 * TWO COORDINATE SYSTEMS — KEEP THEM STRAIGHT (CLAUDE.md vec callout). The
 *   project's vec is ROW-MAJOR on M_n (vec_r(X)[i*n+j] = X[i,j]), used by the
 *   superoperator S_tilde (n^2-vectors). But Co/L_P/R_Q live in A's {B_k}
 *   COORDINATE space (d-vectors, d = dim_A <= n^2), a DIFFERENT space. The bridge
 *   is x_k = <B_k, X>_F (coords of X in A) and X = sum_k x_k B_k (operator from
 *   coords). aic_corner_apply applies Co to X in A via these maps.
 *
 * THE EPS-C* AXIOMS HOLD ONLY UP TO EPS (CLAUDE.md callout). Co's relations
 *   (Co^2 = Co, involution, the L_P R_Q ~ Co bound) are EXACT only after the
 *   theta-snap pins Co to an exact idempotent; the underlying L_P R_Q ~ R_Q L_P
 *   ~ Co agreement is O(delta+eps). The one-sided / symmetric construction does
 *   NOT assume exact associativity or an exact unit in A.
 *
 * EXTRACTION = DOUBLE PATH (aic-w4o.1 deferral, mirrors assoc_extract). dim_S =
 *   round(Re Tr Co), CROSS-CHECKED against the count of singular values > 0.5 of
 *   Co (fail loud on mismatch, Rule 4). The basis {C_m} of S_{P,Q} comes from the
 *   top-dim_S LEFT singular vectors u_m of Co (d-vectors), assembled as operators
 *   C_m = sum_k u_m[k] B_k. Co is a (generally oblique) idempotent, so its nonzero
 *   singular values are >= 1 and the left singular vectors span its range (same
 *   oblique-projector subtlety as aic_assoc_extract_range). Relation DEFECTS are
 *   arb-certified; the SVD is double.
 */
#ifndef AIC_CORNER_H
#define AIC_CORNER_H

#include <flint/acb_mat.h>

#include "aic_ecstar.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Build the d x d coordinate matrix of L_P : X |-> P * X on A (* = the star),
 *   (L_P)_{ij} = <B_i, P * B_j>_F = Tr(B_i^dag (P * B_j)),
 * where P is an n x n delta-projection (ASSERTED Hermitian, fail loud: a
 * delta-projection has P^dag = P, .tex:917). `Lp` must be initialised d x d
 * (d = A->dim_A); P is n x n. All products through aic_ecstar_star. */
void aic_corner_build_L(acb_mat_t Lp, const aic_ecstar *A, const acb_mat_t P,
                        slong prec);

/* Build the d x d coordinate matrix of R_Q : X |-> X * Q on A,
 *   (R_Q)_{ij} = <B_i, B_j * Q>_F.
 * `Rq` init'd d x d; Q n x n, ASSERTED Hermitian. Products through the star. */
void aic_corner_build_R(acb_mat_t Rq, const aic_ecstar *A, const acb_mat_t Q,
                        slong prec);

/* The compression map Co_{P,Q} = theta(L_P R_Q + R_Q L_P - 1) as a d x d
 * idempotent coordinate matrix (.tex:1057). Builds L_P, R_Q (build_L/build_R),
 * forms 2M-1 = L_P R_Q + R_Q L_P - I (ordinary d x d matrix products), and calls
 * aic_prop_P (theta(2M-1); ASSERTS the basin ||M^2-M|| < 1/4, fail loud Rule 4).
 * `Co` init'd d x d. After this Co^2 = Co to machine precision (test T1). */
void aic_corner_Co(acb_mat_t Co, const aic_ecstar *A, const acb_mat_t P,
                   const acb_mat_t Q, slong prec);

/* One-sided Co_{P,1} = theta(L_P) (R_1 = id => symmetric combination is L_P,
 * .tex:1064). `Co` init'd d x d; P n x n Hermitian. */
void aic_corner_Co_P1(acb_mat_t Co, const aic_ecstar *A, const acb_mat_t P,
                      slong prec);

/* One-sided Co_{1,Q} = theta(R_Q) (.tex:1064). `Co` init'd d x d. */
void aic_corner_Co_1Q(acb_mat_t Co, const aic_ecstar *A, const acb_mat_t Q,
                      slong prec);

/* Apply a d x d coordinate-space operator T (e.g. Co) to X in A:
 *   out = sum_m (T x)_m B_m,   x_k = <B_k, X>_F.
 * `out` and `X` are n x n (out != X asserted); `T` is d x d. This is the
 * operator->coords->operator round trip the corner maps act through. */
void aic_corner_apply(acb_mat_t out, const aic_ecstar *A, const acb_mat_t T,
                      const acb_mat_t X, slong prec);

/* Extract the corner subspace S_{P,Q} = Img Co from the idempotent coordinate
 * matrix `Co` (d x d). Writes dim_S to *dim_S_out and allocates/fills *C_out, a
 * flint_malloc'd array of dim_S acb_mat_t (each n x n OPERATOR, C_m = sum_k
 * u_m[k] B_k for u_m the m-th left singular vector of Co; caller frees each then
 * the array). dim_S = round(Re Tr Co) CROSS-CHECKED against the SVD-gap count at
 * 0.5 (fail loud on mismatch, Rule 4). ASSERTS {C_m} Frobenius-orthonormal (left
 * singular vectors are orthonormal in C^d; the map C_m = sum_k u_m[k] B_k
 * preserves the inner product since {B_k} is Frobenius-orthonormal). DOUBLE-path
 * SVD (LAPACK zgesvd via aic_latd_svd). dim_S may be 0 (P ~ 0): then *C_out is
 * set to NULL. */
void aic_corner_extract(acb_mat_t **C_out, slong *dim_S_out, const acb_mat_t Co,
                        const aic_ecstar *A, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CORNER_H */
