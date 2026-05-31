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

#ifdef __cplusplus
}
#endif

#endif /* AIC_IDEMP_H */
