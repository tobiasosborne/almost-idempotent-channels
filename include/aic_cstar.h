/* aic_cstar.h — the constructive proof of th_main (approximate_algebras.tex
 * :1321-1444, section 9): every finite-dimensional eps-C* algebra A is
 * O(eps)-isomorphic to a genuine finite-dimensional C* algebra B = (+)_C M_{|C|},
 * with a dimension-independent constant, REALIZED by an explicit isomorphism.
 * Module bead aic-097 ("cstar_build"); design contract in
 * docs/research/cstar_build_design.md (§7 increment plan).
 *
 * THE MASTER LOOP (the eventual deliverable, .tex:1414-1444). The proof breaks A
 * into pieces and reassembles them: Stage 1 finds a maximal commutative skeleton
 * (a greedy nontrivial-projection loop, lem_nontriv_projection + cor_merge_sum +
 * cor_improvement); Stage 2 grows each equivalence class into a full matrix block
 * (lem_extension, applied inductively M_{r-1} -> M_r); Stage 3 merges the blocks
 * (cor_merge_sum). The inner objects of every stage are the COMPRESSED
 * subalgebras S_P, which must be presented as eps'-C* algebras IN THEIR OWN RIGHT
 * so that projection / extension recurse on them. This header GROWS across the
 * increments I1-I5 of the design doc §7; THIS FILE currently carries I1 only.
 *
 * INCREMENT 1 (I1) — the two pieces of recursion/codomain infrastructure:
 *
 *   (A) aic_cstar_subalgebra: re-present a compressed subalgebra S_P of a parent
 *       eps-C* algebra A as a DERIVED aic_ecstar, whose star is the COMPRESSED
 *       PRODUCT X . Y = Co_P(Phi_tilde(XY)) (.tex:1078-1082, compr_prod) and whose
 *       unit is Ptilde = Co_P(P). This is what lets aic_projection_nontrivial,
 *       aic_corner, aic_dhom recurse on S_P unchanged. The derived ecstar reuses
 *       the parent's star via the star_phi/star_ctx seam (the same seam
 *       aic_assoc_ecstar_from_phi uses for the superoperator star).
 *
 *       THE COMPRESSED STAR, NOT THE PLAIN PRODUCT (CLAUDE.md domain callout,
 *       FINDINGS §C2/§C3, LOAD-BEARING). S_P's product is Co_P(Phi_tilde(XY)),
 *       NOT plain XY and NOT Frobenius-projected. Co_P = aic_corner_Co(A,P,P) is
 *       the OBLIQUE compression idempotent of the corner module; Phi_tilde is the
 *       parent's star map (applied as aic_ecstar_star(parent, XY, I) = parent's
 *       Phi_tilde(XY)). The eta=0 identity channel has Phi_tilde=id so its star is
 *       the plain product and is BLIND to a star bug (FINDINGS §C2) — the wrapper
 *       must therefore be tested on an OBLIQUE eta>0 channel too.
 *
 *   (B) aic_cstar_matrix_algebra: wrap M_d as a GENUINE finite-dim C* algebra in
 *       the aic_ecstar model. Needed later (lem_extension, .tex:1378) as the
 *       codomain B(S_{P,Q}) ~= M_n and as the domain M_n. Its star is the
 *       IDENTITY map (Phi = id), so X*Y = Phi(XY) = XY = plain matrix product, the
 *       genuine product of M_d. Its basis is the d^2 matrix units E_{lm}
 *       (.tex:1371-1374), which are already Frobenius-orthonormal (||E_lm||_F = 1).
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The paper's S_P is an
 * abstract Banach-space compression snapped by the holomorphic functional calculus
 * of prop_P; in finite dim the corner module already realized Co_P as an exact
 * d x d idempotent (eig-free Newton-Schulz theta), so the S_P wrapper here is pure
 * structural glue over those tested numerical routines: it stores Co_P, Ptilde and
 * an extracted Frobenius-orthonormal basis {C_m} of S_P, and routes the star
 * through a thunk computing aic_corner_cdot's compressed product. No new numerical
 * primitive is introduced; the cross-checks are the eta=0 genuine-C* oracle, the
 * thunk-vs-aic_corner_cdot consistency, and the oblique star teeth.
 */
#ifndef AIC_CSTAR_H
#define AIC_CSTAR_H

#include <flint/acb_mat.h>

#include "aic_ecstar.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The star-thunk context for the S_P wrapper. SEPARATELY HEAP-ALLOCATED (NOT
 * inlined into the owner): the derived ecstar's star_ctx points at THIS block,
 * not into the owner struct, so the owner survives being copied/relocated (the
 * I5 master loop builds many wrappers in a relocatable array). Holds the
 * BORROWED parent and an OWNED copy of Co_P (the d x d compression idempotent the
 * thunk applies). */
typedef struct {
    const aic_ecstar *parent;  /* BORROWED parent eps-C* algebra A               */
    acb_mat_t Co_P;            /* OWNED d_parent x d_parent compression idemp Co_P*/
} aic_cstar_sp_ctx;

/* The compressed subalgebra S_P of a parent eps-C* algebra A, presented as a
 * DERIVED aic_ecstar (the field `A`).
 *
 * RELOCATION-SAFE, MOVE-ONLY (the F2 fix). The derived ecstar's star_ctx points
 * at a SEPARATELY HEAP-ALLOCATED ctx block (`ctx`), NOT into this struct — so a
 * by-value copy of the owner, or relocation of the owner inside a `realloc`'d
 * array, leaves star_ctx pointing at the still-valid heap ctx (the old design,
 * star_ctx = &owner, dangled on any move). The owner is therefore SHALLOW-COPYABLE
 * and RELOCATION-SAFE, but it OWNS heap resources (the ctx block, Ptilde, the
 * derived basis): it is MOVE-ONLY. A by-value copy SHARES the ctx/Ptilde/basis
 * with the original; NEVER aic_cstar_subalg_clear two copies (double free).
 *
 * I5 USAGE CONTRACT. Store wrappers in an array (relocation across `realloc` is
 * safe — star_ctx tracks the heap ctx, not the slot address). On each built
 * wrapper call aic_cstar_subalg_clear EXACTLY ONCE; treat a copy as a borrow of
 * the original (do not clear it).
 *
 * This owner holds the heap ctx, Ptilde and the derived basis; aic_cstar_subalg
 * _clear frees all three (+ aic_ecstar_clear on A). The parent is BORROWED (the
 * caller keeps it alive for this struct's lifetime). */
typedef struct {
    const aic_ecstar *parent;  /* BORROWED parent eps-C* algebra A (convenience)  */
    aic_cstar_sp_ctx *ctx;     /* OWNED heap star-thunk ctx (parent + Co_P)       */
    acb_mat_t Ptilde;          /* n x n unit of S_P, = Co_P(P) (.tex:1082)        */
    aic_ecstar A;              /* derived: basis {C_m}, star_phi=thunk, phi=NULL  */
} aic_cstar_subalgebra;

/* Build the S_P wrapper from a parent eps-C* algebra `parent` and a Hermitian
 * delta-projection `P` (n x n, n = parent->n). Fills `out`:
 *   ctx     = heap block holding `parent` and Co_P = aic_corner_Co(parent,P,P)
 *             (d_parent x d_parent); out->A.star_ctx = out->ctx (the heap block,
 *             NOT &out, so the owner survives relocation/copy — F2),
 *   {C_m}   = aic_corner_extract(Co_P, parent)         (dim_S n x n operators),
 *   out->A  : n = parent->n, dim_A = dim_S, B = {C_m}, phi = NULL,
 *             star_phi = the compressed-product thunk,
 *   Ptilde  = aic_corner_Ptilde(parent, P)             (the unit of S_P).
 * The derived star is X*Y = Co_P(Phi_tilde_parent(XY)) (compr_prod, .tex:1080),
 * which EQUALS aic_corner_cdot(parent, Co_P, X, Y) for X, Y in S_P.
 * ASSERTS dim_S >= 1 (fail loud, Rule 4: dim S_P == 0 means P ~ 0, no algebra).
 * Free with aic_cstar_subalg_clear. */
void aic_cstar_subalg_build(aic_cstar_subalgebra *out, const aic_ecstar *parent,
                            const acb_mat_t P, slong prec);

/* Free the heap ctx (incl. its Co_P), Ptilde, and the derived ecstar basis
 * (aic_ecstar_clear). Does NOT touch the borrowed parent. MOVE-ONLY: call this
 * EXACTLY ONCE per built wrapper — a by-value copy shares the ctx/Ptilde/basis,
 * so clearing both copies double-frees. */
void aic_cstar_subalg_clear(aic_cstar_subalgebra *out);

/* Build a GENUINE d x d matrix C* algebra M_d as an aic_ecstar (.tex:1371-1374):
 *   out->n = d, out->dim_A = d*d,
 *   out->B = the d^2 matrix units E_{lm} (each d x d, single 1 at (l,m);
 *            already Frobenius-orthonormal, ||E_lm||_F = 1),
 *   out->phi = NULL, out->star_ctx = NULL,
 *   out->star_phi = a thunk that just COPIES (Phi = identity, so the star is
 *                   X*Y = Phi(XY) = XY = plain matrix multiplication).
 * Every eps-C* axiom defect of this ecstar is machine-zero (a genuine C* algebra).
 * ASSERTS d >= 1. Free with aic_ecstar_clear (the basis is the only owned data;
 * star_phi is a static function, star_ctx is NULL). */
void aic_cstar_matrix_algebra(aic_ecstar *out, slong d, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CSTAR_H */
