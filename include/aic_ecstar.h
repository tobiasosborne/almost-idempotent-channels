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
 */
#ifndef AIC_ECSTAR_H
#define AIC_ECSTAR_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_idemp.h"
#include "aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* An eps-C* algebra A <= M_n with a Frobenius-orthonormal basis and the star's
 * Phi. The Phi pointer is BORROWED: the caller owns the aic_ucp_kraus and MUST
 * keep it alive for the entire lifetime of this struct (the estimators call
 * aic_ucp_apply through it). aic_ecstar_clear frees only the B[] basis, never
 * phi. */
typedef struct {
    slong n;                  /* A <= M_n (n x n complex)                       */
    slong dim_A;              /* d = dim A                                      */
    acb_mat_t *B;             /* d operators, each n x n, <B_j,B_k>_F = delta_jk */
    const aic_ucp_kraus *phi; /* BORROWED; the star's Phi (dim_K==dim_H==n)     */
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

/* A generic linear map B(C^n) -> B(C^n): out = T(X). `ctx` carries T's data
 * (e.g. the aic_ecstar* whose phi is applied, or a synthetic non-HP map). Used
 * as the "star's Phi" by the shared involution-residual core so the SAME code
 * path is exercised by both the production (HP) thunk and a non-HP teeth. */
typedef void (*aic_ecstar_apply_fn)(acb_mat_t out, const acb_mat_t X, void *ctx,
                                    slong prec);

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

#ifdef __cplusplus
}
#endif

#endif /* AIC_ECSTAR_H */
