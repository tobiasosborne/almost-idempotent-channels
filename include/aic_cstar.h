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
#include <flint/arb.h>

#include "aic_dhom.h"
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

/* ============================ Increment 2 =============================== *
 * The two SIMPLER assembly lemmas of §9 (the master loop builds on these):
 * lem_add_dim (.tex:1363, dimension decomposition) and cor_merge_sum (.tex:1352,
 * merging two delta-inclusions). Plus the Stage-3 off-diagonal certification
 * S_{P_C,P_D}=0 (.tex:1443), a lem_add_dim corollary. The full why lives in
 * src/aic_cstar_merge.c. These are STRUCTURAL GLUE over closed routines: the
 * corner dim queries (aic_corner_alpha_dims / aic_corner_dim_S) and the dhom_v
 * concatenation + certification (aic_dhom_defect_sweep / aic_dhom_v_sigma_min);
 * no new numerical primitive is introduced.
 */

/* lem_add_dim (.tex:1363-1369): the dimension-decomposition certification.
 *   ".tex:1364: dim S_{P,Q} = sum_{j,k} dim S_{P_j,Q_k}."
 * P = sum_j P_parts[j], Q = sum_k Q_parts[k] (the caller passes BOTH the parts
 * AND the sums; same convention as aic_corner_alpha / aic_corner_alpha_dims).
 * Reuses aic_corner_alpha_dims (which builds the per-block Co's + Co_{P,Q}, sums
 * the extracted dims into N = sum_{jk} dim S_{P_j,Q_k}, and reports d_PQ =
 * dim S_{P,Q}); then ASSERTS N == d_PQ (fail loud, Rule 4: the lem_add_dim
 * conclusion). Returns the common value d_PQ.
 *   `A`       : the eps-C* algebra (BORROWED).
 *   `P_parts` : p Hermitian n x n delta-projections; `P` their sum (n x n).
 *   `Q_parts` : q Hermitian n x n delta-projections; `Q` their sum (n x n).
 * On an exact idempotent (eta=0) with orthogonal rank-1 P_j, Q_k this is exact:
 * dim S_{P_j,Q_k} in {0,1} and the sum equals dim S_{P,Q} (an isometric
 * bijection). The internal N == d_PQ assert is the teeth (a wrong P_part feeding
 * a mismatched sum ABORTS). */
slong aic_cstar_lem_add_dim(const aic_ecstar *A,
                            const acb_mat_t *P_parts, slong p, const acb_mat_t P,
                            const acb_mat_t *Q_parts, slong q, const acb_mat_t Q,
                            slong prec);

/* Stage-3 off-diagonal certification (.tex:1443): returns 1 iff
 * dim S_{P_C,P_D} == 0. By lem_add_dim, distinct equivalence classes C != D have
 * S_{P_C,P_D} = 0 (the sum over j in C, k in D of dim S_{P_j,P_k} is 0 since each
 * j !~ k contributes 0); this is the certification that the merged blocks of
 * Stage 3 are genuinely independent (cor_merge_sum's bijectivity hypothesis).
 * Thin wrapper over aic_corner_dim_S(A, P_C, P_D). P_C, P_D n x n Hermitian. */
int aic_cstar_off_diag_zero(const aic_ecstar *A, const acb_mat_t P_C,
                            const acb_mat_t P_D, slong prec);

/* cor_merge_sum (.tex:1352-1358): merge two delta-inclusions
 *   v_1 : B_1 -> A,  v_2 : B_2 -> A     (both into the SAME containing
 * eps-C* algebra A — the parent, OR an S_P wrapper when used inside Stage 1; the
 * paper writes v_j : B_j -> S_{P_j}, and S_{P_j} <= A, so v1->A == v2->A == A)
 * into the combined map
 *   ".tex:1355: v : (X_1, X_2) |-> v_1(X_1) + v_2(X_2) : B_1 (+) B_2 -> A."
 * The merged B has block dims = concat(v1->B->d, v2->B->d); the merged v's
 * vE[] is the DIRECT concatenation vE = [v1->vE, v2->vE] (design §3): B's
 * matrix-unit linear index mu(l,a,b) = offset_l + a*d_l + b is CUMULATIVE across
 * blocks (aic_dhom.h), so the concatenation is exactly v(E_i) for the merged B.
 * The merged v is an O(delta+eps)-inclusion into A (.tex:1357); cor_improvement
 * (errreduce) is applied by the CALLER afterward (.tex:1426/1443), NOT here.
 *
 *   B_out, v_out : OUTPUT, allocated here. v_out->B points at B_out (which the
 *                  caller MUST keep alive for v_out's lifetime); v_out->A = A.
 *                  Free with aic_dhom_B_clear(B_out) and aic_dhom_v_clear(v_out)
 *                  (in either order; v_out borrows B_out and A). The vE operators
 *                  are DEEP COPIES of v1/v2's vE (acb_mat_set, not aliased).
 *   mult_def     : OUTPUT certified arb ball, the merged v's multiplicativity
 *                  defect aic_dhom_defect_sweep(v_out) w.r.t. A's STAR (NOT plain
 *                  product, FINDINGS §C2; NULL to skip).
 *   sigma_min    : OUTPUT certified arb ball, aic_dhom_v_sigma_min(v_out) (the
 *                  SOUND Frobenius inclusion lower bound, FINDINGS §C6; NULL skip).
 *
 * ASSERTS (fail loud, Rule 4): v1->A == A, v2->A == A, v1->n == v2->n == A->n.
 * The ||P_1 + P_2 - I_A|| <= delta hypothesis (.tex:1353) is the CALLER's
 * responsibility (Stage 1/3 of the master loop maintains it); this routine does
 * NOT re-verify it. */
void aic_cstar_merge_sum(aic_dhom_B *B_out, aic_dhom_v *v_out,
                         arb_t mult_def, arb_t sigma_min,
                         const aic_dhom_v *v1, const aic_dhom_v *v2,
                         const aic_ecstar *A, slong prec);

/* ============================ Increment 3 =============================== *
 * lem_merging (.tex:1325-1350): the GENERAL 2x2 block-assembly lemma that
 * lem_extension (I4) consumes. cor_merge_sum (I2) is its TWO-BLOCK / off-diagonal-
 * empty special case (.tex:1352). The full why lives in src/aic_cstar_merging.c.
 * STRUCTURAL GLUE over closed routines (the dhom_v concat + certification
 * aic_dhom_defect_sweep / aic_dhom_v_sigma_min and the certified op-norm upper
 * bound aic_corner_gamma_opnorm_ub); no new numerical primitive is introduced.
 *
 * THE TWO B-SHAPES (FINDINGS §C9, load-bearing). .tex:1325 says only that B is "a
 * C* algebra" with two COMPLEMENTARY projections Pi_1 + Pi_2 = I. Two shapes arise:
 *   (single block, two_block=0)  B = M_{n1+n2}: Pi_1 = proj onto first n1 coords,
 *     Pi_2 onto last n2. The off-diagonal sectors S_{Pi_1,Pi_2}, S_{Pi_2,Pi_1} are
 *     the LIVE rectangular blocks of M_{n1+n2} (E_{lm} with l,m across the n1
 *     boundary EXIST). gamma_12, gamma_21 must be NONZERO inclusions of them. This
 *     is what lem_extension (I4, .tex:1378) needs (v_+ : M_{n+1} -> A is a single
 *     matrix block, gamma_12 = U_1, gamma_21 = U_1(.)^dag).
 *   (two block,  two_block=1)    B = M_{n1} (+) M_{n2}: Pi_1 = I of block 1, Pi_2 =
 *     I of block 2. The off-diagonal sectors are EMPTY (dim 0) — block-diagonal B
 *     has no E_{lm} across blocks — so gamma_12 = gamma_21 = 0 (their domains are
 *     zero-dimensional). This is what cor_merge_sum (.tex:1352) uses; merged v then
 *     equals aic_cstar_merge_sum's concat exactly.
 * Conflating the two is a SILENT BUG: a single M_{n1+n2} block with gamma_12 =
 * gamma_21 = 0 is NOT a valid lem_merging input — the live off-diagonal pair
 * (E_{0,n1}, E_{n1,0}) has E*E = E_00 in M_d but gamma(E)=0, so merging1 is
 * violated and the multiplicativity defect is O(1) (= ||Ptilde_1||, ~1), NOT the
 * cor_merge_sum cross-defect ||Ptilde_1 * Ptilde_2|| ~ O(eps). The prompt's
 * "single-block + zero off-diagonal == cor_merge_sum" reading is mathematically
 * incorrect; the genuine reduction is the TWO-block shape (FINDINGS §C9).
 */

/* lem_merging (.tex:1325-1346): assemble four block maps gamma_{jk} into the
 * combined O(delta+eps)-inclusion gamma : B -> A (.tex:1339-1344). B's shape is
 * chosen by `two_block` (see the §C9 banner above).
 *
 * B is a GENUINE C* algebra with two COMPLEMENTARY projections Pi_1, Pi_2,
 * Pi_1 + Pi_2 = I_B EXACTLY. gamma_{jk} maps the (j,k) block S_{Pi_j,Pi_k} to
 * S_{P_j,P_k} <= A (P_1, P_2 delta-projections in A, ||P_1+P_2-I_A|| <= delta, the
 * CALLER's responsibility, .tex:1326). When two_block=1 the off-diagonal blocks
 * are empty so gamma_12, gamma_21 are IGNORED (pass NULL or anything).
 *
 * INPUT — the four gamma_{jk} as their IMAGES on each block's matrix-unit basis
 * (row-major local index a*cols + b, a = local row, b = local col):
 *   gamma11 : n1*n1 operators, gamma11[a*n1 + b] = gamma_{11}(E^{(11)}_{ab}),
 *   gamma12 : n1*n2 operators, gamma12[a*n2 + b] = gamma_{12}(E^{(12)}_{ab}),
 *   gamma21 : n2*n1 operators, gamma21[a*n1 + b] = gamma_{21}(E^{(21)}_{ab}),
 *   gamma22 : n2*n2 operators, gamma22[a*n2 + b] = gamma_{22}(E^{(22)}_{ab}).
 * Each image is an n x n operator IN A (n = A->n). (Contiguous acb_mat_t arrays,
 * indexed gammaJK[idx] as an acb_mat_t, same convention as aic_corner_alpha_dims.)
 * When two_block=1, gamma12/gamma21 are unused.
 *
 * THE ROUTING (get this EXACTLY right; an off-by-one silently mis-assembles).
 *   two_block=0: B = M_d (d = n1+n2, single block). The matrix unit E_{lm}
 *     (l = global row, m = global col) has linear index mu(0,l,m) = l*d + m
 *     (aic_dhom_B_index). Routed by which sector its row/col fall in:
 *       l <  n1, m <  n1 : block (1,1), (a,b) = (l,      m),      gamma11[a*n1 + b],
 *       l <  n1, m >= n1 : block (1,2), (a,b) = (l,      m - n1), gamma12[a*n2 + b],
 *       l >= n1, m <  n1 : block (2,1), (a,b) = (l - n1, m),      gamma21[a*n1 + b],
 *       l >= n1, m >= n1 : block (2,2), (a,b) = (l - n1, m - n1), gamma22[a*n2 + b].
 *     So v_out->vE[l*d + m] = gamma_{jk}(E_{lm}) (a DEEP COPY of the input image).
 *   two_block=1: B = M_{n1} (+) M_{n2} (two blocks). Block 0's matrix unit E_{ab}
 *     (index a*n1 + b) -> gamma11[a*n1 + b]; block 1's E_{ab} (index n1*n1 +
 *     a*n2 + b) -> gamma22[a*n2 + b]. No off-diagonal units exist.
 *
 *   B_out, v_out : OUTPUT, allocated here. B_out is M_{n1+n2} (two_block=0) or
 *                  M_{n1} (+) M_{n2} (two_block=1). v_out->B points at B_out (keep
 *                  B_out alive for v_out's lifetime); v_out->A = A. Free with
 *                  aic_dhom_B_clear(B_out) and aic_dhom_v_clear(v_out). The vE
 *                  operators are DEEP COPIES.
 *   mult_def     : OUTPUT certified arb ball, aic_dhom_defect_sweep(v_out) — the
 *                  inclusion's multiplicativity defect w.r.t. A's STAR (NOT plain
 *                  product, FINDINGS §C2). The defect_sweep over all matrix-unit
 *                  pairs certifies the COMBINED map gamma:B->A multiplicativity
 *                  (the O(delta)-homomorphism property the proof at .tex:1349 needs)
 *                  — stronger than, and implying, the per-block merging1 (.tex:1331)
 *                  hypothesis G_gamma(E_i,E_j). NULL to skip.
 *   sigma_min    : OUTPUT certified arb ball, aic_dhom_v_sigma_min(v_out) — the
 *                  SOUND Frobenius inclusion lower bound (FINDINGS §C6). NULL skip.
 *   merge_cond_max : OUTPUT certified arb ball, the max violation of the cheap
 *                  per-block INPUT guards merging0 + merging2 + merging3:
 *                    (merging0, .tex:1329) max_{block,(a,b)}
 *                      ||gamma_{kj}(E_{ba}) - gamma_{jk}(E_{ab})^dag||_op,
 *                    (merging2, .tex:1334) ||gamma_{11}(Pi_1) - P_1||_op and
 *                      ||gamma_{22}(Pi_2) - P_2||_op (Pi_j = sum of diagonal units),
 *                    (merging3, .tex:1336) max_{block,(a,b)}
 *                      | ||gamma_{jk}(E_{ab})||_op - 1 |  (||E_{ab}||_op = 1).
 *                  When two_block=1 the off-diagonal blocks are skipped (only the
 *                  diagonal merging0/merging2/merging3 contribute). Op-norms via the
 *                  certified mid+radius upper bound (FINDINGS §C5 Gram false-fail
 *                  workaround). NULL to skip. This is REPORTABLE (does NOT abort);
 *                  the I5 loop will gate on it (fail-loud).
 *
 * ASSERTS (fail loud, Rule 4): n1 >= 1, n2 >= 1, A->n == P_j size, the diagonal
 * gamma arrays non-NULL (and the off-diagonal too when two_block=0). Does NOT
 * re-verify ||P_1+P_2-I_A|| <= delta (CALLER's, .tex:1326).
 *
 * On an exact idempotent (eta=0) with the NATURAL identity-inclusion gamma_{jk}
 * (each M_d matrix unit -> itself in A=M_d, two_block=0), the merged gamma is the
 * identity: mult_def ~ machine-zero, sigma_min ~ 1, merge_cond_max ~ 0 (C1 oracle).
 * With two_block=1 and gamma_12 = gamma_21 unused it REDUCES to cor_merge_sum
 * exactly (C2 cross-check). */
void aic_cstar_lem_merging(aic_dhom_B *B_out, aic_dhom_v *v_out,
                           arb_t mult_def, arb_t sigma_min, arb_t merge_cond_max,
                           int two_block, slong n1, slong n2,
                           const acb_mat_t *gamma11, const acb_mat_t *gamma12,
                           const acb_mat_t *gamma21, const acb_mat_t *gamma22,
                           const aic_ecstar *A,
                           const acb_mat_t P1, const acb_mat_t P2, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CSTAR_H */
