/* aic_assoc.h — the associated eps-C* algebra of an almost-idempotent UCP map,
 * STEP 1: the exact-idempotent regularization Phi_tilde = theta(2 Phi - 1)
 * (bead aic-92f, module "assoc_ecsa", Increment 1). Realizes the first half of
 * th_almost_idemp (approximate_algebras.tex:2162-2237).
 *
 * WHAT THIS REALIZES.
 *   Given an eta-idempotent UCP map Phi : B(C^n) -> B(C^n) (||Phi^2-Phi||_cb<=eta,
 *   eta<1/4, .tex:2167-2168), the paper regularizes it to an EXACTLY idempotent
 *   superoperator (.tex:2171-2174):
 *
 *     Phi_tilde = theta(2 Phi - 1)
 *               = (1/2)(1 + sgn(2 Phi - 1))
 *               = (1/2)(1 + (2 Phi - 1)(1 - 4(Phi - Phi^2))^{-1/2}).
 *
 *   This is Proposition prop_P (.tex:524-533, "P^2-P <= delta < 1/4 =>
 *   Ptilde = theta(2P-I) is an exact idempotent") applied to P = Phi AT THE
 *   SUPEROPERATOR LEVEL: Phi is a linear map on M_n, with an n^2 x n^2
 *   superoperator matrix S_Phi, and Phi_tilde is the matrix prop_P(S_Phi).
 *   With A = Img Phi_tilde this opens Increment 2 (the Choi-Effros star, the
 *   ecstar axiom defects); Increment 1 stops at S_Phi and S_Phitilde.
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3).
 *   The paper frames Phi_tilde via the functional calculus of the cb-norm
 *   Banach algebra of completely bounded maps on B(H) (.tex:354-359) — a Taylor
 *   expansion in 4(Phi - Phi^2) that converges for eta<1/4. We do NOT transcribe
 *   that infinite-dim cb-algebra calculus. In finite dimensions a UCP self-map on
 *   M_n IS an n^2 x n^2 matrix S_Phi; theta(2 S_Phi - 1) is then an ORDINARY
 *   matrix-functional-calculus call, computed by the eig-free, certified
 *   funcalc engine (aic_prop_P -> aic_theta -> aic_sgn Newton-Schulz). The
 *   paper's bound ||Phi_tilde - Phi||_cb <= O(eta) (.tex:2179) is the algorithm's
 *   spec; the test suite asserts it (and that it vanishes at eta=0).
 *
 * THE NON-NORMALITY SUBTLETY (load-bearing). S_Phi is NON-NORMAL: Phi is
 *   Hermicity-preserving but NOT self-adjoint w.r.t. the Hilbert-Schmidt inner
 *   product, so S_Phi S_Phi^dag != S_Phi^dag S_Phi in general. funcalc's sgn
 *   (Newton-Schulz, eig-free) was auditioned only on Hermitian/near-normal
 *   inputs; assoc_ecsa is its FIRST non-normal customer. Newton-Schulz uses only
 *   matrix multiply/add and converges for ||X^2-I||_op<1 with NO normality
 *   assumption, so it is sound here — but soundness on non-normal input is
 *   CROSS-CHECKED (test T6: prop_P route vs an independent binomial-series route
 *   via aic_funcalc_xpow, agreeing to ~1e-25 at prec=256 on a genuinely
 *   non-normal almost-idempotent channel).
 *
 * THE FAIL-LOUD BASIN (CLAUDE.md Rule 4). theta/sgn require ||(2S-1)^2 - I||_op
 *   = 4 ||S^2 - S||_op < 1 (.tex:516,520-521); aic_prop_P additionally asserts
 *   the prop_P hypothesis ||S^2-S||_op < 1/4 (.tex:525). This holds at eta=0 and
 *   small eta but TRIPS near eta -> 1/4 or for large-n channels whose defect
 *   inflates. The assert SHOULD fire then (do NOT suppress it); the globally-
 *   convergent out-of-basin sgn variant is bead aic-8hz (out of scope here).
 *
 * SUPEROPERATOR / vec CONVENTION (matches aic_idemp.h, aic_ecstar.h — pinned).
 *   An n x n matrix X is vectorized ROW-MAJOR: vec_r(X)[i*n+j] = X[i,j]. The
 *   superoperator S of a linear map T : M_n -> M_n is the n^2 x n^2 matrix with
 *       vec_r(T(X)) = S vec_r(X),  i.e. S[a*n+b, p*n+q] = T(E_{pq})[a,b],
 *   where E_{pq} = |e_p><e_q|. So COLUMN (p*n+q) of S is vec_r(T(E_{pq})). We
 *   build S_Phi this column-by-column route (reusing the tested aic_ucp_apply +
 *   the tested row-major reshape; the explicit Kronecker formula
 *   S_Phi = sum_a K_a^dag (x) K_a^T is a SECONDARY cross-check, test T2). The
 *   vec convention is this project's #1 historical bug (idemp deliberately
 *   avoided forming the superoperator) — every routine here is oracle-tested
 *   (test T1: S_Phi vec_r(X) == vec_r(aic_ucp_apply(Phi,X))).
 */
#ifndef AIC_ASSOC_H
#define AIC_ASSOC_H

#include <flint/acb_mat.h>

#include "aic/aic_ecstar.h"
#include "aic/aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Build the superoperator matrix S_Phi (n^2 x n^2) of the UCP self-map `phi`
 * (dim_K == dim_H == n; ASSERTED, fail loud). Column (p*n+q) of `S` is
 * vec_r(Phi(E_{pq})), the row-major vectorization of Phi applied to the matrix
 * unit |e_p><e_q| (file docstring). Reuses the tested aic_ucp_apply. `S` must be
 * initialised n^2 x n^2 (asserted). After this, S vec_r(X) = vec_r(Phi(X)) for
 * every n x n X (the oracle, test T1). */
void aic_assoc_superop_from_ucp(acb_mat_t S, const aic_ucp_kraus *phi,
                                slong prec);

/* Build S_Phi via the explicit Kronecker formula (SECONDARY route, for the
 * convention cross-check, test T2):
 *     S_Phi = sum_a K_a^dag (x) K_a^T
 * (left-major Kronecker, left factor K_a^dag, right factor K_a^T), which matches
 * the row-major vec convention. `S` must be initialised n^2 x n^2. Must equal
 * aic_assoc_superop_from_ucp to machine precision. */
void aic_assoc_superop_kron(acb_mat_t S, const aic_ucp_kraus *phi, slong prec);

/* Apply a superoperator to a matrix: out = reshape(S vec_r(X)), i.e.
 *     out[a,b] = sum_{p,q} S[a*n+b, p*n+q] X[p,q].
 * `S` is n^2 x n^2; `X` and `out` are n x n (caller-initialised, out != X
 * required and asserted). Needed by the tests now (T1, T5) and by the
 * Choi-Effros star in Increment 2. */
void aic_assoc_superop_apply(acb_mat_t out, const acb_mat_t S, const acb_mat_t X,
                             slong prec);

/* The regularization Phi_tilde = theta(2 Phi - 1) (.tex:2171-2174). Builds
 * S_Phi via aic_assoc_superop_from_ucp, then S_tilde = aic_prop_P(S_Phi) (the
 * funcalc theta(2P-1) engine; route A1). ASSERTS the prop_P basin
 * ||S_Phi^2 - S_Phi||_op < 1/4 inside aic_prop_P (fail loud, Rule 4 — an
 * out-of-basin eta is NOT silently regularized; bead aic-8hz). `S_tilde` must be
 * initialised n^2 x n^2 (asserted). For an EXACTLY idempotent Phi (eta=0),
 * S_tilde == S_Phi to machine precision (the headline eta=0 oracle, test T3). */
void aic_assoc_regularize(acb_mat_t S_tilde, const aic_ucp_kraus *phi,
                          slong prec);

/* ===================== Increment 2: A = Img Phi_tilde + the star ============
 *
 * th_almost_idemp (.tex:2184-2195). With Phi_tilde = theta(2 Phi - 1) exactly
 * idempotent (Increment 1, the superoperator S_tilde), the associated algebra is
 *     A = Img Phi_tilde = Ker(1 - Phi_tilde) <= B(H)        (.tex:2185-2186)
 * a closed subspace containing I = 1_H and invariant under X |-> X^dag, with the
 * APPROXIMATE Choi-Effros product (.tex:2189 Choi-Effros)
 *     X * Y = Phi_tilde(XY)                                  (X,Y in A).
 * Theorem th_almost_idemp: (A, *, ||.||, dag, I) is an extended O(eta)-C* algebra.
 *
 * THE OWNER struct. We build A as an aic_ecstar (the verification data model;
 * include/aic_ecstar.h), whose star is Phi_tilde via the generic star_phi seam
 * (X*Y = Phi_tilde(XY)). The ecstar struct only BORROWS its star_phi/star_ctx, so
 * the assoc layer OWNS the superop matrix S_tilde + a small ctx and frees them.
 * aic_assoc_ecstar holds the aic_ecstar + the owned S_tilde + ctx; build with
 * aic_assoc_ecstar_from_phi, free with aic_assoc_ecstar_clear (which also clears
 * the embedded ecstar's basis). */
typedef struct {
    acb_mat_t S_tilde;   /* owned n^2 x n^2 idempotent superoperator Phi_tilde  */
    void *ctx;           /* owned ctx wrapping S_tilde for the star_phi thunk   */
    aic_ecstar A;        /* the eps-C* data model; A.star_phi = the thunk       */
} aic_assoc_ecstar;

/* A-EXTRACTION (Route-B1). Given the idempotent superop S_tilde (n^2 x n^2),
 * extract A = range(S_tilde) as a Frobenius-ORTHONORMAL operator basis {B_k}:
 *   dim_A = round(Re Tr S_tilde)  (idempotent => trace = rank = #unit eigenvalues),
 *           CROSS-CHECKED against the thin-SVD gap of S_tilde: the count of
 *           singular values above a gap threshold must equal round(trace), else
 *           FAIL LOUD (Rule 4).
 *   B_k    = reshape_row-major(top-k LEFT singular vector of S_tilde) into n x n.
 * Left singular vectors are orthonormal in C^{n^2} = Frobenius-orthonormal as
 * operators (asserted, <B_j,B_k>_F = delta_jk to tol).
 *
 * THE PROJECTOR-SVD SUBTLETY (load-bearing). S_tilde is idempotent
 * (S_tilde^2 = S_tilde) but in general NOT Hermitian (Phi_tilde is HP but not
 * always HS-self-adjoint, aic_assoc.h non-normality note). For an idempotent the
 * nonzero SINGULAR values are >= 1, with EQUALITY (all == 1, an ORTHOGONAL
 * projector, SVD spectrum exactly {1,...,1,0,...}) when Phi_tilde IS
 * HS-self-adjoint -- the self-dual channels: dephasing, AND the dep(d)+conj(dep(d))
 * eta>0 family (a convex mix of HS-self-adjoint maps stays HS-self-adjoint, so its
 * S_tilde is ORTHOGONAL, sigma_max == 1; measured, test U5) -- and STRICTLY > 1
 * only when Phi_tilde is genuinely OBLIQUE (non-HS-self-adjoint), e.g. the
 * compression channel compress_idemp(4,2), whose S_tilde inflates sigma_max to
 * sqrt(3) ~ 1.732 (measured, U5). Either way rank = trace still holds and
 * range(S_tilde) is correctly spanned by the top-dim_A LEFT singular vectors, so
 * the 0.5 gap separates nonzero (>= 1) from zero (~0) cleanly in both cases.
 * Extraction runs in the DOUBLE path (LAPACK zgesvd, fast; the certified-
 * degenerate-eig wall aic-w4o.1 precedent: extraction is double, the DEFECT checks
 * are arb-certified); the basis is returned as acb_mat at `prec`.
 *
 * `B_out` is a caller-allocated array of >= dim_A acb_mat_t pointers? NO -- instead
 * this fills an aic_ecstar via aic_assoc_ecstar_from_phi; the standalone extractor
 * below is exposed for the tests. It WRITES dim_A to *dim_A_out, allocates and
 * fills *B_out (a flint_malloc'd array of dim_A acb_mat_t, each n x n; caller frees
 * each then the array). */
void aic_assoc_extract_range(acb_mat_t **B_out, slong *dim_A_out,
                             const acb_mat_t S_tilde, slong prec);

/* Build the associated eps-C* algebra A = Img Phi_tilde of the eta-idempotent
 * UCP map `phi` (dim_K==dim_H==n; ASSERTED). Steps:
 *   (a) S_tilde = aic_assoc_regularize(phi)  (asserts the prop_P basin, Rule 4);
 *   (b) A's basis {B_k} = aic_assoc_extract_range(S_tilde);
 *   (c) out->A.star_phi = a superop-apply thunk over out->ctx (= S_tilde),
 *       out->A.phi = NULL, n = n, dim_A = dim_A, B = {B_k}.
 * `out` is OUTPUT (every field allocated here). Free with aic_assoc_ecstar_clear.
 * For an EXACTLY idempotent phi (eta=0) this A matches ecstar's A=Img Phi from
 * aic_ecstar_from_idemp as a SUBSPACE (the eta=0 consistency oracle, test U1). */
void aic_assoc_ecstar_from_phi(aic_assoc_ecstar *out, const aic_ucp_kraus *phi,
                               slong prec);

/* Free the owned S_tilde + ctx AND the embedded ecstar's basis (aic_ecstar_clear).
 * Mirrors ecstar's borrowed-phi rule: the ecstar only borrows star_phi/star_ctx,
 * so this owner is responsible for the superop's lifetime. */
void aic_assoc_ecstar_clear(aic_assoc_ecstar *out);

#ifdef __cplusplus
}
#endif

#endif /* AIC_ASSOC_H */
