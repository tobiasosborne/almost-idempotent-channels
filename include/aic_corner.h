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

/* ============================ Increment 2a =============================== *
 * The compressed product (compr_prod, .tex:1077-1082) and the block bijection
 * alpha (lem_alpha, .tex:1086-1119). The full why lives in
 * src/aic_corner_product.c.
 */

/* The COMPRESSED PRODUCT X . Y = Co_{P,R}(X * Y) (.tex:1078-1080, eq compr_prod),
 *   (X,Y) |-> X . Y : S_{P,Q} x S_{Q,R} -> S_{P,R}.
 * Computes Z = X * Y (the STAR aic_ecstar_star, n x n) then applies the d x d
 * compression Co_{P,R} (prebuilt by aic_corner_Co) to Z via aic_corner_apply.
 *   `out`   : n x n (n = A->n), caller-initialised, must NOT alias X or Y.
 *   `CoPR`  : d x d, = Co_{P,R} (A->dim_A x A->dim_A).
 *   `X`,`Y` : n x n; X in S_{P,Q}, Y in S_{Q,R} (the caller's responsibility).
 * The paper's spec ||X . Y - X * Y|| <= O(delta+eps) ||X|| ||Y|| (.tex:1081) is
 * the correctness bound; tests/test_corner.c asserts it. The STAR is LOAD-BEARING
 * (CLAUDE.md callout): dropping it (plain X.Y) is RED on the oblique fixture. */
void aic_corner_cdot(acb_mat_t out, const aic_ecstar *A, const acb_mat_t CoPR,
                     const acb_mat_t X, const acb_mat_t Y, slong prec);

/* The unit of (S_P, .): Ptilde = Co_P(P) = Co_{P,P}(P) (.tex:1082). Builds
 * Co_{P,P} (aic_corner_Co with Q=P) and applies it to P (aic_corner_apply).
 *   `out` : n x n, caller-initialised, must NOT alias P.
 *   `P`   : n x n Hermitian delta-projection.
 * Then ||Ptilde . X - X||, ||X . Ptilde - X|| <= O(delta+eps) ||X|| for X in S_P
 * (.tex:1082; machine-zero at eta=0); tests assert it. */
void aic_corner_Ptilde(acb_mat_t out, const aic_ecstar *A, const acb_mat_t P,
                       slong prec);

/* The block bijection alpha of lem_alpha (.tex:1086-1119). Given p Hermitian
 * P_j and q Hermitian Q_k with P = sum_j P_j, Q = sum_k Q_k (the caller passes
 * the SUMS and the PARTS), builds
 *   alpha   = sum_{jk} Co_{P,Q}|_{S_{Pj,Qk}} : (+)_{jk} S_{Pj,Qk} -> S_{P,Q}
 *   beta    = sum_{jk} Co_{Pj,Qk}            : S_{P,Q} -> (+)_{jk} S_{Pj,Qk}
 * as explicit matrices between the direct-sum COORDINATE space (dim N =
 * sum_{jk} dim S_{Pj,Qk}, blocks contiguous) and S_{P,Q}'s coordinate space
 * (dim d_PQ = dim S_{P,Q}). The columns of alpha are alpha_{jk}(C^{jk}_m) =
 * Co_{P,Q}(C^{jk}_m) re-coordinated in S_{P,Q}'s basis {D_l}; the columns of
 * beta are beta_{jk}(D_l) = Co_{Pj,Qk}(D_l) re-coordinated in {C^{jk}_m}.
 *
 * .tex TYPO (escalated, bead aic-czm): tex:1109 prints beta_{jk} = Co_{P_j,Q_j};
 * the proof tex:1114 (delta_jl delta_km) and the codomain S_{Pj,Qk} both require
 * Co_{P_j,Q_k}. We implement Co_{Pj,Qk} (see src/aic_corner_product.c).
 *
 * Computes gamma = beta alpha - I_N and ASSERTS a CERTIFIED UPPER BOUND on
 * ||gamma||_op is < 1 (the lem_alpha hypothesis pq(delta+eps) < const; fail loud,
 * Rule 4). The bound is ||mid(gamma)||_op + ||rad(gamma)||_F >= ||gamma||_op
 * (FIX 2 — the midpoint alone discarded gamma's certified radius, a soundness hole
 * per the arb ladder). Then alpha is a bijection (so N == d_PQ is ASSERTED, the
 * dim-count oracle .tex:1124) with the certified inverse alpha^{-1} =
 * (beta alpha)^{-1} beta (acb_mat_solve, certified).
 *
 * OUTPUTS (all caller-initialised to the shapes below; *N_out, *dPQ_out written):
 *   `alpha`   : d_PQ x N. `alpha_inv` : N x d_PQ (NULL to skip).
 *   `norm_alpha`, `norm_alpha_inv` : initialised arb_t, certified upper bounds on
 *               ||alpha||_op, ||alpha^{-1}||_op (NULL to skip; informational, they
 *               gate nothing).
 *   `gamma_norm` : initialised arb_t, the CERTIFIED UPPER BOUND on ||gamma||_op
 *               (the contraction constant tested < 1; NULL to skip). Always
 *               computed internally for the assert.
 *   `gamma_rad`  : initialised arb_t, the ||rad(gamma)||_F radius contribution
 *               folded into gamma_norm (NULL to skip; ~1e-72 in the tight regime,
 *               so gamma_norm is dominated by ||mid(gamma)||_op). Lets the caller
 *               report how much certified width the guard accounts for.
 * Because the output shapes depend on N (data-dependent), alpha/alpha_inv are
 * caller-allocated via the dims aic_corner_alpha_dims reports first. */
void aic_corner_alpha(acb_mat_t alpha, acb_mat_t alpha_inv,
                      arb_t norm_alpha, arb_t norm_alpha_inv, arb_t gamma_norm,
                      arb_t gamma_rad,
                      slong *N_out, slong *dPQ_out,
                      const aic_ecstar *A,
                      const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                      const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                      slong prec);

/* Report the direct-sum dimension N = sum_{jk} dim S_{Pj,Qk} and d_PQ =
 * dim S_{P,Q} WITHOUT building alpha (so the caller can allocate alpha/alpha_inv
 * to the right shape, then call aic_corner_alpha). Builds the per-block Co's and
 * Co_{P,Q}, sums the extracted dims. */
void aic_corner_alpha_dims(slong *N_out, slong *dPQ_out, const aic_ecstar *A,
                           const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                           const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                           slong prec);

/* ============================ Increment 2b ============================== *
 * One-dimensional-Q machinery: the lem_PQ_Hilb inner product (the Hilbert
 * structure on S_{P,Q} for 1d Q), the Ha-map (the implicit-equation solve that
 * lem_extension depends on), lem_PQR (norm multiplicativity), lem_1d_proj and
 * the 1d-equivalence predicate. The full why lives in src/aic_corner_hilbert.c
 * (the Hilbert-space results) and src/aic_corner_ha.c (the Ha-map).
 * approximate_algebras.tex:1123-1187.
 */

/* lem_PQ_Hilb (.tex:1123-1132). The inner product <Y,X> on S_{P,Q} for a ONE-
 * DIMENSIONAL delta-projection Q, defined by Y^dag . X = <Y,X> Qtilde where
 *   Y^dag . X = Co_Q(Y^dag * X)   (the STAR Y^dag * X = aic_ecstar_star, then Co_Q),
 *   Qtilde   = Co_Q(Q)  (an operator in A; pass it PREBUILT, = aic_corner_Ptilde
 *              on Q, an element of S_Q),
 * and S_Q = C Qtilde is one-dimensional. The scalar is then
 *   <Y,X> = <Qtilde, W>_F / <Qtilde, Qtilde>_F,   W = Co_Q(Y^dag * X).
 * Writes the COMPLEX scalar to `out` (caller-init'd acb_t). `CoQ` = Co_{Q,Q}
 * (d x d, prebuilt by aic_corner_Co with both args Q). X, Y are n x n elements of
 * S_{P,Q}. The 1d-Q precondition (dim S_Q == 1) is the CALLER's to assert (see
 * aic_corner_dim_S); this routine asserts only that <Qtilde,Qtilde>_F is bounded
 * away from 0 (fail loud, Rule 4). The inner product is conjugate-linear in Y,
 * linear in X, Hermitian-symmetric, and |<X,X> - ||X||^2| <= O(delta+eps)||X||^2
 * (.tex:1130, the Hilbert-norm closeness; tests report the constant). At eta=0 on
 * M_n it recovers the GNS inner product Tr(Y^dag X) restricted to S_{P,Q}. */
void aic_corner_ip_1d(acb_t out, const aic_ecstar *A, const acb_mat_t Qtilde,
                      const acb_mat_t CoQ, const acb_mat_t X, const acb_mat_t Y,
                      slong prec);

/* The Ha-map Ha^Q_{P,R}(Z) (.tex:1146-1160, eqs Ha_def/Ha_dag/Ha_prod), a linear
 * map S_{R,Q} -> S_{P,Q} for Z in S_{P,R} and a ONE-DIMENSIONAL Q. Returned as a
 * d_PQ x d_RQ matrix `Ha` in the EXTRACTED corner bases {C^{PQ}_l} of S_{P,Q}
 * (rows) and {C^{RQ}_m} of S_{R,Q} (columns): column m is the coords of
 * Ha^Q_{P,R}(Z)(C^{RQ}_m) in {C^{PQ}_l}. Solves the implicit Ha_def
 *   (Y^dag . Z) . X + Y^dag . (Z . X) = 2 <Y, Ha(Z)(X)> Qtilde   (all Y in S_{P,Q})
 * by the Gram system: with G_{lm} = <C^{PQ}_l, C^{PQ}_m> (the lem_PQ_Hilb inner
 * product, NOT the Frobenius ip; G = I + O(delta+eps), ASSERTED ||G-I||<1 fail
 * loud) and b_l = (1/2)[scalar((C^{PQ}_l^dag . Z) . X) + scalar(C^{PQ}_l^dag .
 * (Z . X))] for input X = C^{RQ}_m, Ha(Z)(C^{RQ}_m) = sum_l (G^{-1} b)_l C^{PQ}_l
 * (acb_mat_solve, certified). scalar(W in S_Q) = <Qtilde,W>_F/<Qtilde,Qtilde>_F.
 *
 * THE Co-INDEX BOOKKEEPING (the module's highest convention-risk; each cdot's
 * target corner is the Co passed):
 *   Y^dag . Z      : S_{Q,P} x S_{P,R} -> S_{Q,R}   via Co_{Q,R}
 *   (Y^dag.Z) . X  : S_{Q,R} x S_{R,Q} -> S_Q       via Co_Q
 *   Z . X          : S_{P,R} x S_{R,Q} -> S_{P,Q}   via Co_{P,Q}
 *   Y^dag . (Z.X)  : S_{Q,P} x S_{P,Q} -> S_Q       via Co_Q
 * Y^dag = C^{PQ}_l^dag (acb_mat_conjugate_transpose, lands in S_{Q,P}).
 *
 * The exact identity Ha^Q_{R,P}(Z^dag) = Ha^Q_{P,R}(Z)^dag (.tex:1153) holds by
 * the symmetric definition; ||Ha(Z)(X) - Z.X||_Euc <= O(delta+eps)||Z||||X||_Euc
 * (.tex:1155); Ha^Q_{P,P} is an O(delta+eps)-HOMOMORPHISM (.tex:1157, the
 * downstream-critical property for lem_extension). Tests assert all three.
 *
 * `Ha` must be init'd d_PQ x d_RQ (d_PQ = dim S_{P,Q}, d_RQ = dim S_{R,Q};
 * query via aic_corner_dim_S). Z is n x n in S_{P,R}. P,R,Q are n x n Hermitian
 * delta-projections; Q one-dimensional (dim S_Q == 1 ASSERTED). All Co's and
 * bases are built INTERNALLY from P,R,Q. */
void aic_corner_ha(acb_mat_t Ha, const aic_ecstar *A, const acb_mat_t Z,
                   const acb_mat_t P, const acb_mat_t R, const acb_mat_t Q,
                   slong prec);

/* lem_PQR defect (.tex:1162-1177). For a ONE-DIMENSIONAL Q, the maximum over the
 * extracted bases X in S_{P,Q}, Y in S_{Q,R} of the norm-multiplicativity defect
 *   | ||X . Y||_op - ||X||_op ||Y||_op |   (NO eigenproblem; pure norm test),
 * returned as a certified arb ball in `out` (caller-init'd). X . Y via
 * aic_corner_cdot with Co_{P,R}. Bounded by O(delta+eps)||X||||Y|| (the spec);
 * machine-zero at eta=0 (||X.Y|| = ||X|| ||Y|| exactly). Builds all Co's/bases
 * from P,Q,R internally; asserts dim S_Q == 1 (fail loud). */
void aic_corner_pqr_defect(arb_t out, const aic_ecstar *A, const acb_mat_t P,
                           const acb_mat_t Q, const acb_mat_t R, slong prec);

/* dim S_{P,Q} (.tex:1064-1066): build Co_{P,Q} and return round(Re Tr Co)
 * CROSS-CHECKED against the SVD-gap count (the same fail-loud cross-check as
 * aic_corner_extract; this is the cheap query that does not allocate the basis).
 * P,Q n x n Hermitian. */
slong aic_corner_dim_S(const aic_ecstar *A, const acb_mat_t P, const acb_mat_t Q,
                       slong prec);

/* 1d-equivalence predicate (lem_1d_proj + .tex:1187). Returns 1 iff
 * dim S_{P,Q} == 1 (the paper's definition of equivalent one-dimensional
 * delta-projections P ~ Q), else 0. For one-dimensional P and Q, dim S_{P,Q} <= 1
 * (lem_1d_proj, .tex:1179); the relation is reflexive, symmetric (via the
 * involution Co_{P,Q}^dag = Co_{Q,P}) and transitive (via lem_PQR). P,Q n x n
 * Hermitian. Thin wrapper over aic_corner_dim_S. */
int aic_corner_equiv_1d(const aic_ecstar *A, const acb_mat_t P, const acb_mat_t Q,
                        slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CORNER_H */
