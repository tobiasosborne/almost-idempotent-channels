/* aic_idemp.h — the structure theorem for finite-dimensional EXACTLY idempotent
 * UCP maps (bead aic-wuh, module "idemp_structure"). Realizes th_idemp_structure
 * (approximate_algebras.tex:318 statement, :2055 proof), with lem_idemp (:1916),
 * lem_carrier (:1724), PhiX_M (:1740), prop_Gamma (:2106) for context. This is
 * the eta=0 oracle (milestone aic-9kk) — the cleanest ground truth in the project.
 *
 * THE THEOREM (.tex:318-331). Given an idempotent UCP map Phi : B(H) -> B(H)
 * (dim H = n), there exist a subspace M <= H with inclusion J_M : M -> H and a
 * C* algebra A with UCP maps Delta, Gamma such that
 *
 *     C_M : X |-> J_M^dag X J_M : B(H) -> B(M)             (compression, .tex:325)
 *     Gamma C_M Delta = 1_A,    Delta Gamma C_M = Phi.      (.tex:328 Gamma_C_Delta)
 *
 * w = C_M Delta : A -> B(M) is a C*-homomorphism, and operators in Im Delta are
 * block-diagonal w.r.t. H = M (+) M^perp.
 *
 * THE CONSTRUCTION (the proof IS the algorithm, .tex:2055-2091):
 *   1. Carrier M = range(Q), Q = sum_a K_a K_a^dag (lem_carrier :1724). M = H iff
 *      Phi is trace-preserving (sum_a K_a K_a^dag = 1); M can be a PROPER subspace
 *      when Phi is unital but not TP. J_M = (eigenvectors of Q with eval > thr),
 *      an n x dim_M isometry; Pi_M = J_M J_M^dag.
 *   2. A = Img Phi, basis {B_k} from the column-image of Phi on the matrix units
 *      E_{ij}: stack vec(Phi(E_{ij})) (n^2 columns), SVD, keep left singular
 *      vectors with sigma > thr, reshape to n x n. dim_A = numerical rank.
 *      Phi(B_k) = B_k (fixed points).
 *   3. Delta : A -> B(H) = inclusion = the n^2 x dim_A matrix [vec B_1|...|vec B_d]
 *      with ORTHONORMAL columns (Delta^dag Delta = 1_A). Im Delta = Img Phi.
 *   4. C_M : X |-> J_M^dag X J_M (compression; just a triple product).
 *   5. Lambda : B(M) -> B(H), Lambda(Y) = Phi(J_M Y J_M^dag) (.tex:2061); CP and
 *      unital (Lambda(1_M)=1_H, .tex:2072). Stored as Lambda_mat (n^2 x dim_M^2).
 *   6. Gamma : B(M) -> A via Lambda = Delta Gamma (.tex:2065). Since Delta has
 *      orthonormal columns, Gamma_mat = Delta^dag Lambda_mat (dim_A x dim_M^2).
 *      The paper GUARANTEES this Gamma is UCP (.tex:2088). NOT a w-pseudoinverse.
 *   7. w = C_M Delta : A -> B(M) (.tex:2084), w_mat (dim_M^2 x dim_A), column k =
 *      vec(J_M^dag B_k J_M). A *-homomorphism (lem_idemp :1916).
 *
 * VECTORIZATION CONVENTION (load-bearing — pinned here). An n x n matrix X is
 * vectorized ROW-MAJOR: vec(X)[i*n + j] = X[i,j]. So a linear map T : B(C^p) ->
 * B(C^q) is stored as a (q*q) x (p*p) matrix T_mat with
 *     vec(T(X))[a*q + b] = sum_{i,j} T_mat[a*q+b, i*p+j] vec(X)[i*p+j].
 * This matches aic_ucp.h's E_{ij} indexing and aic_mat.h's left-major Kronecker.
 *
 * NUMBER PATHS. The extraction (carrier eig, image SVD) is the DEGENERATE-
 * spectrum double path (an exactly-idempotent Phi has Q with spectrum {0,1} and
 * a maximally-degenerate image SVD), via aic_latd (LAPACKE_zheev/zgesvd); the
 * certified arb degenerate eig/SVD is blocked on aic-w4o.1. The arb path here
 * does NOT re-extract; it CERTIFIES the five defining relations (.tex:2080-2090)
 * given the double-path subspaces, by applying maps to matrix units and comparing
 * via aic_mat_opnorm (never forming the n^2 x n^2 superoperator — convention-safe).
 *
 * GAUGE FREEDOM. The A-basis {B_k} is unique only up to a unitary on C^{dim_A}.
 * Cross-run / double-vs-arb comparisons must compare the SUBSPACE via the
 * projector Pi_A = Delta Delta^dag (= sum_k vec(B_k) vec(B_k)^dag), NOT individual
 * basis vectors. The relations (Gamma C_M Delta = 1_A, etc.) are gauge-invariant.
 *
 * DEFERRED (not in this module): the Artin-Wedderburn block decomposition of A,
 * prop_Gamma's explicit Tr_{E_j} conditional-expectation form, and the Kraus
 * reps of Delta/Gamma -> bead aic-ynu. Certified degenerate eig/SVD -> aic-w4o.1.
 * ETA > 0 (almost-idempotent, via Phi_tilde) is assoc_ecsa's job, NOT this module:
 * aic_idemp_decompose ASSERTS exact idempotence at entry and aborts otherwise.
 */
#ifndef AIC_IDEMP_H
#define AIC_IDEMP_H

#include <flint/acb_mat.h>

#include "aic/aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The decomposition of an exactly-idempotent UCP self-map Phi : B(C^n) -> B(C^n).
 * All matrices are stored as acb_mat_t (the certified type) at the working prec,
 * but are populated from the LAPACK double path (degenerate eig/SVD). Linear maps
 * are stored as their MATRIX in the row-major vec convention (see file docstring).
 *
 * Lifecycle: zero-initialise the struct, call aic_idemp_decompose, free with
 * aic_idemp_clear. */
typedef struct {
    slong n;        /* dim H */
    slong dim_M;    /* dim of the carrier M = rank(Q) */
    slong dim_A;    /* dim of A = Img Phi (numerical rank of the image SVD) */

    acb_mat_t J_M;       /* n x dim_M isometry, columns = ONB of M (J_M^dag J_M = 1_M) */
    acb_mat_t J_Mperp;   /* n x (n-dim_M) isometry, columns = ONB of M^perp */
    acb_mat_t Pi_M;      /* n x n orthogonal projector onto M (= J_M J_M^dag) */

    acb_mat_t Delta;     /* (n*n) x dim_A, columns = vec(B_k); Delta^dag Delta = 1_A */
    acb_mat_t Lambda;    /* (n*n) x (dim_M*dim_M): Lambda(Y) = Phi(J_M Y J_M^dag) */
    acb_mat_t Gamma;     /* dim_A x (dim_M*dim_M) = Delta^dag Lambda */
    acb_mat_t w;         /* (dim_M*dim_M) x dim_A, column k = vec(J_M^dag B_k J_M) */
} aic_idemp_decomp;

/* Decompose an EXACTLY idempotent UCP self-map. `phi` must have dim_K==dim_H==n.
 * `out` is OUTPUT: every field is allocated/filled here (caller must
 * aic_idemp_clear it). `prec` is the arb working precision for the assembled
 * matrices; the eig/SVD extraction runs in double regardless.
 *
 * Fails loud (Rule 4) if: phi is not a self-map; phi is not exactly idempotent
 * (max_ij ||Phi(Phi(E_ij)) - Phi(E_ij)||_op exceeds a prec-appropriate tol —
 * this module is ONLY for exact idempotents); the extracted Phi(B_k) are not
 * fixed points. */
void aic_idemp_decompose(aic_idemp_decomp *out, const aic_ucp_kraus *phi,
                         slong prec);

void aic_idemp_clear(aic_idemp_decomp *out);

/* --- certified relation defects (the eta=0 oracle; .tex:2080-2090) ---
 * Each returns a certified arb ball that is 0 (to rounding) for an exact
 * idempotent. The maps are applied to matrix units / basis coords and compared
 * via aic_mat_opnorm; the n^2 x n^2 superoperator is NEVER formed. `out` arb_t
 * must be initialised. */

/* ||Gamma C_M Delta - 1_A||_op (.tex:2081). */
void aic_idemp_defect_GCD(arb_t out, const aic_idemp_decomp *d, slong prec);

/* max_ij ||Delta Gamma C_M (E_ij) - Phi(E_ij)||_op (.tex:2072, = Phi). */
void aic_idemp_defect_DGC(arb_t out, const aic_idemp_decomp *d,
                          const aic_ucp_kraus *phi, slong prec);

/* max_{j,k} ||w(B_j star B_k) - w(B_j) w(B_k)||_op where B_j star B_k =
 * Phi(B_j B_k) is the Choi-Effros product (lem_idemp :1916: w is a
 * *-homomorphism). */
void aic_idemp_defect_w_hom(arb_t out, const aic_idemp_decomp *d,
                            const aic_ucp_kraus *phi, slong prec);

/* max_k ( ||(1-Pi_M) B_k Pi_M||_op + ||Pi_M B_k (1-Pi_M)||_op ): operators in
 * Im Delta are block-diagonal w.r.t. M (+) M^perp (lem_idemp :1916, .tex:2090). */
void aic_idemp_defect_blockdiag(arb_t out, const aic_idemp_decomp *d, slong prec);

/* Gamma unital defect: ||Delta Gamma(1_M) - 1_H||_op (Gamma(1_M)=1_A i.e.
 * Delta Gamma(1_M) = Lambda(1_M) = 1_H, .tex:2072). The CP-ness of Gamma is
 * certified by the PAIR (aic_idemp_psi_choi + aic_idemp_defect_wG_eq_CML); see
 * those declarations. `out` initialised. */
void aic_idemp_defect_gamma_unital(arb_t out, const aic_idemp_decomp *d,
                                   slong prec);

/* Build the Choi matrix of Psi = C_M Lambda : B(M) -> B(M) (.tex:2084), a
 * (dim_M*dim_M) x (dim_M*dim_M) Convention-A Choi, so the caller can certify
 * CP via aic_ucp_is_cp_choi. This certifies that the C_M o Lambda FACTORIZATION
 * is CP — i.e. the abstract algebra A's Choi-Effros product is CP. It is built
 * from d->Lambda and d->J_M ALONE (it does NOT read d->Gamma), so it does NOT by
 * itself certify the STORED Gamma; pair it with aic_idemp_defect_wG_eq_CML.
 * Working through Psi on B(M) is gauge-invariant and avoids an abstract
 * A-coordinate Choi. (Psi is also idempotent, .tex:2084.) `C` must be
 * initialised (dim_M*dim_M) x (dim_M*dim_M). */
void aic_idemp_psi_choi(acb_mat_t C, const aic_idemp_decomp *d, slong prec);

/* ||w Gamma - C_M Lambda||_op (.tex:2084: Psi = C_M Lambda = w Gamma). The TEETH
 * for the Gamma-CP certificate: reads BOTH d->w and d->Gamma and ties the stored
 * Gamma to the CP-certified Psi (aic_idemp_psi_choi). Since w is an injective
 * *-monomorphism (completely isometric, .tex:2088) and Psi = w Gamma is CP, this
 * equation holding makes the stored Gamma CP. A corrupt Gamma (scaled, column
 * zeroed) gives an O(1) defect here (measured: Gamma*=2 -> 1.0).
 *
 * PRECISION FLOOR. w Gamma = C_M Delta Delta^dag Lambda = C_M Pi_A Lambda routes
 * through the double-path projector Pi_A = Delta Delta^dag (idempotent only to
 * ~1e-16 at prec=53), so this defect is NOT machine-zero like the gauge-clean
 * relations: measured ceiling 1.6e-12 (noiseless_subsystem(2,2), non-exact
 * 1/sqrt2 extraction), <= 1e-15 elsewhere. Gate at 1e-9, not 1e-12. `out`
 * initialised. */
void aic_idemp_defect_wG_eq_CML(arb_t out, const aic_idemp_decomp *d, slong prec);

/* --- Artin-Wedderburn / multiplicity decomposition of A (bead aic-ynu, I1) ---
 *
 * prop_hom_structure (.tex:259-272) + prop_Gamma (.tex:2106). The *-homomorphism
 * w : A = Img Phi -> B(M) (column k of d->w is vec(w(B_k)), the m x m matrix
 * w(B_k) = J_M^dag B_k J_M, m = dim_M) has image a *-subalgebra W <= B(M). By
 * Artin-Wedderburn it is isomorphic to (+)_j B(L_j) acting on
 *     M = (+)_j L_j (x) E_j,    w(B) = sum_j W_j^dag (A_j (x) 1_{E_j}) W_j,
 * where A_j is B's component in the j-th block and W = (W_1,...,W_m) : M ->
 * (+)_j (L_j (x) E_j) is UNITARY (.tex:259-272). This struct holds that data:
 * num_blocks m, the block (matrix) dims dim L_j, the multiplicity dims dim E_j,
 * and the per-block isometries W_j : M -> L_j (x) E_j (each (dim L_j * dim E_j) x
 * dim_M; sum_j (dim L_j * dim E_j) = dim_M, so W stacked is dim_M x dim_M unitary).
 *
 * TENSOR LAYOUT (load-bearing, matches aic_mat.h left-major Kronecker). Within a
 * block, L_j is the LEFT/MAJOR factor: W_j row index t = a*dim_E_j + c with a in
 * [0,dim L_j) the L_j (matrix) index and c in [0,dim E_j) the E_j (multiplicity)
 * index, so (A_j (x) 1_{E_j})[a*dim_E_j+c, a'*dim_E_j+c'] = A_j[a,a'] delta_{cc'}
 * (aic_mat_kronecker convention). The w-reconstruction invariant (below) pins this
 * order together with the sizes and the W_j entries at once.
 *
 * SCOPE (I1). At eta=0 the random-Hermitian simple-spectrum extraction handles
 * BOTH distinct AND equal block sizes: distinct sizes (block_cond_exp(5,2) =
 * M_2 (+) M_3), EQUAL sizes (block_cond_exp(2k,k) = M_k (+) M_k, including the
 * 12x12 M_6 (+) M_6), single blocks with nontrivial multiplicity
 * (noiseless_subsystem), single trivial blocks (identity), and the fully-diagonal
 * M_1^n (dephasing). The generic random H_W separates equal-size blocks almost
 * surely; the re-draw budget covers the measure-zero collision set. The routine
 * FAILS LOUD (Rule 4) only on a genuinely-unresolvable spectrum (no clean split in
 * the budget) or a failed arb certification — it does NOT return garbage.
 *
 * GAUGE. W_j is unique only up to a unitary on each L_j and on each E_j (and the
 * per-block carrier scaling of w). The w-reconstruction residual is gauge-invariant
 * and is the correctness specification. */
typedef struct {
    slong dim_M;        /* dim of the carrier M (= sum_j dim_L[j]*dim_E[j])       */
    slong dim_A;        /* dim of A (= sum_j dim_L[j]^2)                          */
    slong num_blocks;   /* m, the number of Wedderburn blocks                    */
    slong *dim_L;       /* dim_L[0..m-1], the block (matrix) dimensions          */
    slong *dim_E;       /* dim_E[0..m-1], the multiplicity dimensions            */
    acb_mat_t *W_j;     /* m matrices, W_j[j] is (dim_L[j]*dim_E[j]) x dim_M      */
} aic_idemp_wedderburn;

/* Compute the Artin-Wedderburn decomposition of A = Img Phi from an exactly-
 * idempotent decomposition `d` (the eta=0 case). `out` is OUTPUT: num_blocks,
 * dim_L[], dim_E[] and the W_j[] are all allocated/filled here (caller must
 * aic_idemp_wedderburn_clear it). `prec` is the arb working precision.
 *
 * ROUTE (commutant-eigendecomposition; src/aic_idemp_wedderburn.c). Reshape d->w
 * into the m x m matrices {w(B_k)} (m = dim_M); they span the *-subalgebra
 * W <= B(M). (1) Eigen-cluster a GENERIC RANDOM Hermitian H_W in W (drawn with a
 * FIXED-SEED PRNG for reproducibility, then the double-path zheev
 * aic_latd_eig_hermitian — NOT the certified-arb degenerate-eig machinery; the
 * extraction is double, the defects are arb-certified, matching aic_idemp): each
 * cluster projector lies in W and equals (v v^dag) (x) 1_{E_j} for an L_j
 * eigenvector v, with multiplicity dim E_j. The randomness is load-bearing — a
 * deterministic H_W manufactures a kernel spanning equal-size blocks; a degenerate
 * draw triggers a re-draw (bounded budget). (2) Group clusters into blocks by
 * W-connectivity (||Q_a w(B_k) Q_a'||_op > tol iff a, a' are in the same block);
 * dim L_j = #clusters in the block, dim E_j = the common cluster multiplicity.
 * (3) Build each W_j by aligning the E_j basis across the dim L_j reference
 * subspaces via the connecting partial isometries in W. The output is then
 * CERTIFIED IN ARB (dim identities + W-unitarity + the w-reconstruction
 * correctness spec) before returning. Cites .tex:259, .tex:2106.
 *
 * Fails loud (Rule 4) ONLY for a genuinely-unresolvable spectrum: no draw in the
 * re-draw budget yields a clean block split, or the assembled W fails the arb
 * certification (dim identity, unitarity, or w-reconstruction). */
void aic_idemp_wedderburn_decompose(aic_idemp_wedderburn *out,
                                    const aic_idemp_decomp *d, slong prec);

void aic_idemp_wedderburn_clear(aic_idemp_wedderburn *out);

/* --- prop_Gamma explicit Kraus form of the conditional expectation (aic-ynu, I3) ---
 *
 * prop_Gamma (.tex:2106-2122). Given the Wedderburn data W (the unitary W : M ->
 * (+)_j L_j (x) E_j with components W_j) of the *-homomorphism w : A -> B(M),
 * ANY UCP map Gamma : B(M) -> A with Gamma w = 1_A (the conditional expectation
 * onto A) has the component-wise representation (.tex:2110, eq Gamma)
 *
 *     Gamma_j : B(M) -> B(L_j),
 *     Gamma_j(X) = Tr_{E_j}( W_j X W_j^dag (1_{L_j} (x) gamma_j) )   (eq Gamma)
 *
 * for density matrices gamma_j on E_j (the general form of a finite-dim
 * conditional expectation). Equivalently the Choi/Kraus form (.tex:2118-2122)
 *
 *     Gamma_j(X) = L_j^dag (X (x) 1_{F_j}) L_j,   L_j^dag L_j = 1_{L_j},
 *     L_j = (W_j^dag (x) 1_{F_j})(1_{L_j} (x) xi_j),
 *
 * with xi_j a unit vector in E_j (x) F_j and gamma_j = Tr_{F_j}(xi_j xi_j^dag).
 * The A-element Gamma(X) = (Gamma_1(X),...,Gamma_m(X)) in A = (+)_j B(L_j); under
 * the *-monomorphism w it embeds as w(Gamma(X)) = sum_j W_j^dag (Gamma_j(X) (x)
 * 1_{E_j}) W_j (.tex:2095 eq Aw).
 *
 * WHICH gamma_j (extracted to MATCH the stored d.Gamma, NOT presumed uniform).
 * The stored d.Gamma = Delta^dag Lambda (aic_idemp_build.c) is a SPECIFIC
 * conditional expectation; eq Gamma is LINEAR in gamma_j, so we SOLVE for the
 * gamma_j that reproduces it (least squares over a spanning set of B(M) inputs,
 * the target Gamma_j(X) read block-wise from d.Gamma through w). The solution is
 * unique (W_j isometry onto block j => gamma_j |-> Gamma_j is injective) and is
 * asserted to be a valid density matrix (Hermitian, PSD, trace 1; arb-certified,
 * Rule 4). For the noiseless-subsystem oracle Phi(X)=(Tr_E X)(x)(1_E/dE) the
 * stored d.Gamma turns out to be the MAXIMALLY-MIXED conditional expectation,
 * gamma_1 = (1/dim E_1) 1_{E_1} (measured residual 4.5e-16); whether d.Gamma is
 * uniform is DETERMINED by the solve, not presumed (a non-uniform gamma_j is a
 * valid output and breaks the uniform-gamma assumption — the mutation proof).
 *
 * THE GATE (.tex:2113, the crux tooth). The prop_Gamma Kraus map MUST equal the
 * stored d.Gamma: applying both to a basis E_ab of B(M) and comparing the
 * A-coordinate (embedded via w, gauge-invariant since w is a *-monomorphism),
 * ||w(Gamma_kraus(E_ab)) - w(d.Gamma(E_ab))||_op <= tol for all a,b. This pins
 * the gamma_j extraction + the Tr_{E_j} convention + the W_j usage at once. The
 * decompose routine FAILS LOUD (Rule 4) if the extracted gamma_j is not a valid
 * density matrix or the match-d.Gamma certification fails (do NOT return a Gamma
 * that does not match — that would be silently wrong). Per .tex:2113 EVERY Gamma
 * with Gamma w = 1_A has this form, so a genuine mismatch means the W_j are wrong
 * or the convention is off — a real finding, escalated, not faked.
 *
 * F_j = E_j with the standard purification: gamma_j = sum_c p_c |v_c><v_c|
 * (eigendecomposition), xi_j = sum_c sqrt(p_c) |v_c>_E (x) |c>_F, so
 * Tr_{F_j}(xi_j xi_j^dag) = gamma_j and L_j is dim_M*dim_E[j] x dim_L[j]. */
typedef struct {
    slong num_blocks;   /* = W->num_blocks                                       */
    slong dim_M;        /* = W->dim_M (the carrier dim, for re-embedding via w)  */
    slong dim_A;        /* = W->dim_A                                            */
    slong *dim_L;       /* dim_L[0..m-1] (copied from W, self-contained)         */
    slong *dim_E;       /* dim_E[0..m-1] = dim_F[j] (F_j = E_j)                  */
    acb_mat_t *gamma_j; /* m density matrices, gamma_j[j] is dim_E[j] x dim_E[j] */
    acb_mat_t *L_j;     /* m Kraus isometries, L_j[j] is (dim_M*dim_E[j]) x dim_L[j],
                         * L_j^dag L_j = 1_{L_j}; Gamma_j(X)=L_j^dag(X(x)1_F)L_j  */
} aic_idemp_gamma;

/* Extract the prop_Gamma density matrices gamma_j and Kraus operators L_j of the
 * stored conditional expectation d->Gamma, from the Wedderburn data W (same `d`).
 * `out` is OUTPUT (caller must aic_idemp_gamma_clear it). `prec` is the arb
 * working precision.
 *
 * Fails loud (Rule 4) if: an extracted gamma_j is not a valid density matrix
 * (Hermitian / PSD / trace 1 beyond a prec tol); the prop_Gamma map does NOT
 * match the stored d->Gamma on a basis of B(M) (the crux tooth — a mismatch is a
 * convention/W_j bug, not a silent miss); or the Kraus L_j fail L_j^dag L_j = 1.
 * Cites .tex:2106-2122. */
void aic_idemp_gamma_kraus(aic_idemp_gamma *out, const aic_idemp_wedderburn *W,
                           const aic_idemp_decomp *d, slong prec);

void aic_idemp_gamma_clear(aic_idemp_gamma *out);

#ifdef __cplusplus
}
#endif

#endif /* AIC_IDEMP_H */
