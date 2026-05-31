/* aic_ecstar.h — the eps-C* algebra DATA MODEL + certified-arb axiom-defect
 * estimators (bead aic-knm, module "ecstar"). Realizes the eps-C* axiom block
 * (approximate_algebras.tex:407-439) and the Choi-Effros star
 * (approximate_algebras.tex:341-344, :2187-2189) that th_almost_idemp
 * (approximate_algebras.tex:2192-2237) certifies is an O(eta)-C* algebra.
 *
 * WHAT AN eps-C* ALGEBRA IS HERE. A subspace A <= M_n (n x n complex), dim A = d,
 * given by a Frobenius-ORTHONORMAL operator basis {B_1,...,B_d} (each n x n,
 * <B_j,B_k>_F = delta_jk), together with a unital CP map Phi : B(C^n) -> B(C^n)
 * (an aic_ucp_kraus with dim_K == dim_H == n) defining the multiplication
 *
 *     X * Y = Phi(X . Y)        (the Choi-Effros STAR; .tex:341, :2189)
 *
 * where X . Y is the ORDINARY matrix product and Phi is applied via
 * aic_ucp_apply. The involution is the matrix adjoint X |-> X^dag; the unit is
 * I = 1_n (the ambient identity, inherited, .tex:2186-2187); the norm is the
 * OPERATOR norm aic_mat_opnorm (inherited from B(H), .tex:2192). These are the
 * paper's choices; the cb-norm appears ONLY in the definition of eta and is NOT
 * used inside the estimators (the estimators consume A with its inherited
 * operator norm).
 *
 * DOMAIN CALLOUTS (CLAUDE.md). (1) The eps-C* axioms hold only UP TO eps
 * (associativity, submultiplicativity, the C* identity, the unit laws are all
 * approximate, .tex:407-439); do NOT assume exact associativity. (2) The product
 * is the Choi-Effros star X*Y = Phi(XY), NOT plain XY (plain multiplication
 * leaves Img Phi). (3) Norm = OPERATOR norm, not cb-norm and not Frobenius.
 *
 * THE ESTIMATORS ARE BASIS-SWEEP LOWER BOUNDS / EXACT-ZERO DETECTORS. Each defect
 * function below is a sweep over the stored orthonormal basis (its triples /
 * pairs / singletons), NOT the faithful sup over the operator-norm unit ball.
 * They serve the eta=0 oracle (bead aic-knm route decision): for an EXACTLY
 * idempotent Phi with A = Img Phi they read ~0, and a deliberate non-C* instance
 * makes them demonstrably nonzero (the module's teeth). The associator and the
 * involution defect are EXACT zero-detectors (multilinear / always-zero
 * invariants); submult and C* are basis-sweep LOWER bounds on the true eps. The
 * faithful worst-case search (HOPM) and the certified SDP upper bound are later
 * cycles (beads aic-0at). See ALGORITHM.md "Module ecstar" for the full analysis.
 *
 * LOAD-BEARING vec CONVENTION (matches aic_idemp.h). An n x n matrix X is
 * vectorized ROW-MAJOR: vec(X)[i*n+j] = X[i,j]. aic_ecstar_from_idemp reshapes
 * column k of d->Delta this way: B_k[i,j] = Delta[i*n+j, k].
 *
 * THE STAR'S MAP: KRAUS Phi -OR- A GENERIC SUPEROPERATOR (Increment 2 of
 * assoc_ecsa, bead aic-92f). The exact-idempotent case (aic_ecstar_from_idemp)
 * has the star X*Y=Phi(XY) with Phi a UCP Kraus map (the `phi` field). The
 * almost-idempotent case A = Img Phi_tilde (Phi_tilde = theta(2 Phi - 1),
 * .tex:2171-2174, :2189) needs the star X*Y = Phi_tilde(XY) where Phi_tilde is a
 * SUPEROPERATOR (an n^2 x n^2 matrix S_tilde, NOT a UCP/Kraus map: Phi_tilde is
 * Hermicity-preserving but NOT completely positive, .tex:363). So the star is
 * generalised through the apply-fn seam this header already ships: two OPTIONAL
 * fields `star_phi` / `star_ctx`. aic_ecstar_star: if star_phi != NULL it computes
 * XY then star_phi(out, XY, star_ctx, prec); ELSE the Kraus path Phi(XY) via
 * aic_ucp_apply on `phi`. aic_ecstar_init sets star_phi = star_ctx = NULL so the
 * Kraus path is byte-identical (existing tests are the regression guard). The
 * superop constructor lives in the assoc module: aic_assoc_ecstar_from_phi
 * (include/aic_assoc.h), which sets phi = NULL, star_phi to a superop-apply thunk.
 */
#ifndef AIC_ECSTAR_H
#define AIC_ECSTAR_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* A generic linear map B(C^n) -> B(C^n): out = T(X). `ctx` carries T's data
 * (e.g. the aic_ecstar* whose phi is applied, or a synthetic non-HP map, or a
 * superoperator matrix). Used as the "star's Phi" by the shared involution-
 * residual core AND (when set on the struct below) by the star itself. */
typedef void (*aic_ecstar_apply_fn)(acb_mat_t out, const acb_mat_t X, void *ctx,
                                    slong prec);

/* An eps-C* algebra A <= M_n with a Frobenius-orthonormal basis and the star's
 * map. TWO mutually-exclusive ways to define the star X*Y:
 *   (Kraus)   star_phi == NULL: X*Y = Phi(XY) with Phi the BORROWED `phi` (the
 *             aic_ucp_kraus the caller owns and MUST keep alive for the lifetime
 *             of this struct; the estimators call aic_ucp_apply through it).
 *   (superop) star_phi != NULL: X*Y = star_phi(XY, star_ctx) with star_phi a
 *             generic apply (e.g. an n^2 x n^2 superoperator Phi_tilde). `phi`
 *             is then NULL. star_phi/star_ctx are BORROWED too; the assoc-module
 *             constructor (aic_assoc_ecstar_from_phi) owns them and supplies an
 *             aic_assoc_ecstar_clear that frees the superop + ctx + the basis.
 * aic_ecstar_clear frees only the B[] basis (never phi, star_phi, or star_ctx). */
typedef struct {
    slong n;                  /* A <= M_n (n x n complex)                       */
    slong dim_A;              /* d = dim A                                      */
    acb_mat_t *B;             /* d operators, each n x n, <B_j,B_k>_F = delta_jk */
    const aic_ucp_kraus *phi; /* BORROWED; the Kraus star's Phi, or NULL        */
    aic_ecstar_apply_fn star_phi; /* OPTIONAL generic star map; NULL => use phi */
    void *star_ctx;           /* BORROWED ctx for star_phi (e.g. the superop)   */
} aic_ecstar;

/* Allocate B[0..dim_A-1], each an n x n zeroed acb_mat_t; set n, dim_A, phi.
 * ASSERTS phi is a self-map with dim_K == dim_H == n (fail loud, Rule 4). Free
 * with aic_ecstar_clear. */
void aic_ecstar_init(aic_ecstar *out, slong n, slong dim_A,
                     const aic_ucp_kraus *phi);

/* Release the B[] basis (every B[k] and the array). Does NOT touch phi (borrowed). */
void aic_ecstar_clear(aic_ecstar *out);

/* Fill the basis from an exactly-idempotent decomposition: B_k[i,j] =
 * d->Delta[i*n+j, k] (reshape column k ROW-MAJOR). This is the eta=0 case
 * A = Img Phi. ASSERTS on entry that d->n == phi->dim_H, dims match, and the
 * resulting columns are Frobenius-orthonormal (||B_j,B_k>_F - delta_jk|| small);
 * fail loud otherwise. `out` must be aic_ecstar_init'd to (d->n, d->dim_A, phi). */
void aic_ecstar_from_idemp(aic_ecstar *out, const aic_idemp_decomp *d,
                           const aic_ucp_kraus *phi, slong prec);

/* The Choi-Effros star: out = Phi(X . Y) (.tex:341, :2189). X, Y, out are all
 * n x n (caller-initialised). Computes the ordinary product X.Y via acb_mat_mul
 * then applies Phi via aic_ucp_apply. */
void aic_ecstar_star(acb_mat_t out, const aic_ecstar *A, const acb_mat_t X,
                     const acb_mat_t Y, slong prec);

/* Projection residual ||X - sum_k <B_k,X>_F B_k||_op as a certified arb ball: how
 * far X is from the subspace A (Frobenius projection). 0 iff X in A. Used by the
 * unit estimator to check 1_n in A. <B_k,X>_F = sum_ij conj(B_k[i,j]) X[i,j].
 * `out` must be an initialised arb_t; X is n x n. */
void aic_ecstar_proj_residual(arb_t out, const aic_ecstar *A, const acb_mat_t X,
                              slong prec);

/* --- axiom-defect estimators (basis sweeps; certified arb balls). Each cites
 * its axiom .tex line above the definition in aic_ecstar_defect.c. `out` must be
 * an initialised arb_t. ---------------------------------------------------- */

/* ax_assoc (.tex:412-413  ||(XY)Z - X(YZ)|| <= eps ||X|| ||Y|| ||Z||):
 * max over basis triples (i,j,k) of ||(B_i*B_j)*B_k - B_i*(B_j*B_k)||_op. The
 * associator is trilinear, so vanishing on all basis triples <=> vanishing
 * identically: this is an EXACT zero-detector. Reads ~0 for eta=0. */
void aic_ecstar_defect_assoc(arb_t out, const aic_ecstar *A, slong prec);

/* ax_prodnorm (.tex:410-411  ||XY|| <= (1+eps) ||X|| ||Y||):
 * max over basis pairs (j,k) of the positive part
 * max(0, ||B_j*B_k||_op - ||B_j||_op ||B_k||_op). =0 for a genuine (sub-
 * multiplicative) C* algebra. A basis-sweep LOWER bound on the true eps_sub.
 * Reads ~0 for eta=0. */
void aic_ecstar_defect_submult(arb_t out, const aic_ecstar *A, slong prec);

/* ax_C* (.tex:427-428  ||X^dag X|| >= (1-eps) ||X||^2):
 * max over basis k of max(0, 1 - ||B_k^dag * B_k||_op / ||B_k||_op^2). =0 for a
 * genuine C* algebra. A basis-sweep LOWER bound on the true eps_C*. FAILS LOUD if
 * any ||B_k||_op is below a prec floor (a Frobenius-unit op has ||.||_op>=1/sqrt n).
 * Reads ~0 for eta=0. */
void aic_ecstar_defect_cstar(arb_t out, const aic_ecstar *A, slong prec);

/* ax_* (.tex:422-423  (XY)^dag = Y^dag X^dag, EXACT):
 * max over basis pairs (j,k) of ||(B_j*B_k)^dag - B_k^dag * B_j^dag||_op.
 *
 * STRUCTURAL INVARIANT, NOT A FREE MEASUREMENT. With X*Y = Phi(XY) one has
 *   (B_j*B_k)^dag - B_k^dag*B_j^dag = Phi(B_j B_k)^dag - Phi((B_j B_k)^dag),
 * which is IDENTICALLY 0 for ANY Hermicity-preserving (HP) Phi (.tex:2211
 * "(X*Y)^dag = Y^dag*X^dag"). The aic_ucp_kraus data model Phi(X)=sum K^dag X K
 * is HP for EVERY Kraus set, so NO input through this public API can make the
 * production estimator nonzero — it is a consistency check that the star's Phi
 * is HP, eta-independent. A correct impl and `return 0` are observationally
 * identical on all aic_ecstar inputs.
 *
 * The production code PATH is exercised — and a return-0 / wrong-matrix mutation
 * is made detectable — by feeding a deliberately NON-HP linear map through the
 * shared core aic_ecstar_involution_core below: see tests/test_ecstar.c
 * teeth_non_hp, which passes a synthetic T(X)=M X N with M != N^dag and asserts
 * the core's residual GROWS with the non-HP-ness (O(t)). That teeth, not this
 * always-zero call, is what protects the residual computation. */
void aic_ecstar_defect_involution(arb_t out, const aic_ecstar *A, slong prec);

/* Shared involution-residual core (the metric of aic_ecstar_defect_involution,
 * generalised in the map). Returns
 *   max over pairs (j,k) of ||apply(B_j B_k)^dag - apply(B_k^dag B_j^dag)||_op,
 * where {B_k} = B[0..d-1] (each n x n) and apply plays the role of the star's
 * Phi (so the (j,k) term is (B_j*B_k)^dag - B_k^dag*B_j^dag when apply=Phi).
 * For an HP apply this is identically 0; for a non-HP apply it is generally
 * nonzero — that is how a return-0 mutation of this core is detected (Rule 7).
 * `out` initialised; `apply`/`ctx` supply the map; B is a d-array of n x n. */
void aic_ecstar_involution_core(arb_t out, const acb_mat_t *B, slong d, slong n,
                                aic_ecstar_apply_fn apply, void *ctx, slong prec);

/* ax_eps_unit (.tex:432-434): the max of
 *   (a) ||1_n - Pi_A(1_n)||_op   (is the unit IN A?  aic_ecstar_proj_residual),
 *   (b) max_k ||B_k * I - B_k||_op,
 *   (c) max_k ||I * B_k - B_k||_op.
 * EXACTLY 0 when each B_k is a Phi fixed point (X*I = Phi(X) = X for X in Img Phi,
 * .tex:2211 "X*I = X = I*X"). An always-zero invariant for unital Phi, A=Img Phi. */
void aic_ecstar_defect_unit(arb_t out, const aic_ecstar *A, slong prec);

/* --- Cycle 2: FAITHFUL worst-case defect search (HOPM) -------------------- *
 * The basis-sweep estimators above are EXACT-ZERO detectors / cheap LOWER
 * bounds; the routines below are the second audition candidate (Law 4): a
 * scale-invariant alternating-maximization (higher-order power method) search
 * for the worst-case axiom defect over the OPERATOR-norm unit ball of the
 * subspace A. Each returns a RIGOROUS LOWER BOUND on the true sup (the `lo`
 * end of a certified arb ball): the search runs in the fast double path for
 * speed, but the FINAL witness it returns (an explicit operator IN A) is
 * evaluated for its defect ratio in the certified arb path (aic_ecstar_star +
 * aic_mat_opnorm), so `lo` provably under-estimates the sup. There is no global
 * guarantee — any in-A witness gives a valid lower bound (web leg Q8).
 *
 * WHY this is dimension-independent (the universality canary). The search works
 * directly on the spectral-norm unit sphere of A: every iterate X = sum_k x_k B_k
 * is in A (so the witness is valid) and is normalized by ||.||_op (not ||.||_F),
 * so it does NOT pay the d^{3/2} Frobenius-to-operator inflation that makes the
 * basis sweep grow with dim (web leg Cycle-2 decision; HANDOFF universality canary).
 *
 * THE SUBSPACE-POLAR SUBTLETY (load-bearing). The block update maximizes the
 * linear functional Re<C, X>_F over {X in A, ||X||_op <= 1} where the Frobenius
 * gradient C = sum_k c_k B_k is IN A (c_k = <u, h(B_k,...) v>). polar(C)=U V^dag
 * (from the SVD of C) is the exact maximizer over the FULL M_n operator-ball, but
 * for an APPROXIMATE C* algebra (eta>0) polar(C) is generally NOT in A (a genuine
 * C* algebra is polar/von-Neumann closed; ours is only eps-close). Using polar(C)
 * directly would yield a witness OUTSIDE A and an INVALID lower bound. So the step
 * is X' = Pi_A(polar(C)) (project the polar factor back onto A), accepted only if
 * the scale-invariant ratio STRICTLY improves (monotone-ascent guard); else X is
 * kept. The accept-guard keeps every accepted iterate exactly in A and makes the
 * approximate-polar step robust.
 *
 * Deterministic restarts (no wall-clock RNG, project rule): each restart's init
 * is seeded from a fixed counter derived from the restart index, so the result is
 * reproducible. Half the restarts warm-start from the Frobenius-relaxation optimum
 * (top singular vectors / random in A), half from Haar-ish random operators in A.
 *
 * Contract: `lo` (initialised arb_t) receives a LOWER bound on the true sup. The
 * routines ASSERT n_restarts >= 1 and n_iters >= 0 (n_iters==0 is the "never
 * iterate" mutation probe: it returns only the init witness value). */

/* ax_assoc (.tex:412-413): lower bound on
 *   eps_assoc = sup_{X,Y,Z in A, !=0}
 *               ||(X*Y)*Z - X*(Y*Z)||_op / (||X||_op ||Y||_op ||Z||_op). */
void aic_ecstar_defect_assoc_hopm(arb_t lo, const aic_ecstar *A, int n_restarts,
                                  int n_iters, slong prec);

/* ax_prodnorm (.tex:410-411): lower bound on
 *   sup_{X,Y in A} ||X*Y||_op / (||X||_op ||Y||_op) - 1
 * (expected ~0: any UCP star is submultiplicative; the search must confirm ~0). */
void aic_ecstar_defect_submult_hopm(arb_t lo, const aic_ecstar *A, int n_restarts,
                                    int n_iters, slong prec);

/* ax_C* (.tex:427-428): lower bound on
 *   eps_cstar = sup_{X in A} (1 - ||X^dag * X||_op / ||X||_op^2)
 *             = 1 - min_{||X||_op=1} ||X^dag * X||_op   (a MINIMIZATION inside). */
void aic_ecstar_defect_cstar_hopm(arb_t lo, const aic_ecstar *A, int n_restarts,
                                  int n_iters, slong prec);

/* assoc HOPM with the EXPLICIT witness exposed (for the certified-ball soundness
 * test, Rule 7). Same as aic_ecstar_defect_assoc_hopm but also writes the worst-
 * case witness operators into wX,wY,wZ (caller-initialised n x n acb_mat), so a
 * test can check (i) each is in A (aic_ecstar_proj_residual ~0) and (ii) the
 * returned `lo` <= the ratio re-evaluated at (wX,wY,wZ). `lo` is the same
 * rigorous lower bound. wX/wY/wZ may be NULL to skip. */
void aic_ecstar_defect_assoc_hopm_witness(arb_t lo, acb_mat_t wX, acb_mat_t wY,
                                          acb_mat_t wZ, const aic_ecstar *A,
                                          int n_restarts, int n_iters, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_ECSTAR_H */
