/* aic_dhom.h — delta-homomorphisms between a finite-dim C* algebra B and an
 * eps-C* algebra A (bead aic-c1n, module "dhom"). Realizes §8 of Kitaev
 * arXiv:2405.02434 (approximate_algebras.tex:1190-1319): prop_delta_hominc
 * (:1194-1222) and lem_approx (:1224-1314). FULLY CONSTRUCTIVE: pure matrix
 * products + certified operator norms; NO eigendecomposition anywhere.
 *
 * WHY THIS MODULE EXISTS. lem_approx is the error-reduction engine that keeps the
 * th_main isomorphism constant DIMENSION-INDEPENDENT (.tex:484 — naive Haar /
 * cohomology constructions of the diagonal have error proportional to n=dim A;
 * the EXPLICIT generalized-Pauli diagonal of the genuine C* algebra B has
 * ||D||_proj = 1 exactly, so the correction it drives does NOT grow with dim).
 * cor_improvement (.tex:1317, the NEXT module errreduce aic-t81) consumes this:
 * it iterates lem_approx + prop_delta_hominc to a c_0*eps inclusion.
 *
 * --------------------------------------------------------------------------
 * DATA MODEL.
 *
 *   B = (+)_{l=0}^{m-1} M_{d_l}(C), a GENUINE finite-dimensional C* algebra.
 *     num_blocks = m; d[l] = block dims; dim_B = sum_l d_l^2; n_B = sum_l d_l
 *     (the size of the block-diagonal n_B x n_B matrix that represents an element
 *     of B). Matrix-unit basis E^{(l)}_{ab}, a,b in [0,d_l). Linear index
 *       mu(l,a,b) = offset_l + a*d_l + b,   offset_l = sum_{l'<l} d_{l'}^2.
 *     B's PRODUCT is BLOCK MATMUL (a genuine C* algebra, exact). Its involution
 *     is the block adjoint; its unit is I_B = (+)_l I_{d_l}; ||.|| = operator
 *     norm of the n_B x n_B block-diagonal representative.
 *
 *   v : B -> A, a delta-homomorphism, stored as dim_B operators in A (each n x n
 *     where n = A.n, the AMBIENT dim of A <= M_n): vE[i] = v(E_i), and for a
 *     general X in B with B-coords x_i (x_i = the matrix-unit entries of X, read
 *     block by block), v(X) = sum_i x_i vE[i] (LINEAR in B-coords).
 *
 *   THE STAR (LOAD-BEARING — CLAUDE.md domain callout). EVERY product taken IN A
 *   is the Choi-Effros star X*Y = Phi_tilde(XY) via aic_ecstar_star, NEVER plain
 *   acb_mat_mul. Plain matmul is used ONLY for B's block products and for the
 *   Pauli matrices. The multiplicativity defect is
 *       G_v(X,Y) = v(XY) - v(X) * v(Y)        (XY in B; * = star in A)
 *   (.tex:1258). The eta=0 identity-channel oracle has star == plain product, so
 *   it CANNOT catch a star-vs-plain bug; only an eta>0 fixture exercises the
 *   star non-trivially (test T4 mutates the star to plain matmul and confirms RED;
 *   the corner-module lesson, bead aic-czm).
 *
 * --------------------------------------------------------------------------
 * THE PAULI DIAGONAL (the dimension-independence crux, .tex:1250-1254).
 *   Per block M_d, the generalized Pauli (shift-clock / Heisenberg-Weyl)
 *     S_{jk} = X^j Z^k,   <e_l|S_{jk}|e_m> = exp(2 pi i k m / d) delta_{l,(m+j) mod d}.
 *   Each S_{jk} is unitary. The single-block diagonal
 *     D_l = sum_{j,k} d^-2 S_{jk}^dag (x) S_{jk}
 *   is the unique norm-1 diagonal of M_d (.tex:1249): XD_l = D_l X for all X
 *   (centrality), pi(D_l) = sum d^-2 S_{jk}^dag S_{jk} = I_d, ||D_l||_proj = 1.
 *   The diagonal of (+)_l M_{d_l} is the cross-term-free EMBEDDED SUM
 *     D = sum_l D_l = sum_l sum_{jk} d_l^-2 (Shat^{(l)}_{jk})^dag (x) Shat^{(l)}_{jk},
 *   where Shat^{(l)}_{jk} is S^{(l)}_{jk} embedded in sector l of the n_B x n_B
 *   block-diagonal representative and ZERO elsewhere (a partial isometry). This is
 *   the Haar diagonal integral dU (U^dag (x) U) and IS central; the paper's
 *   .tex:1254 joint-unitary Cartesian-product formula is NON-central for the finite
 *   Pauli design — see the prominent .tex:1254 DISCREPANCY box in aic_dhom_diag.c.
 *   nterms = sum_l d_l^2 (a SUM over blocks). Stored as a flat list of (weight, U)
 *   pairs (the A_j (x) B_j of .tex:1247 are A_j = weight * U^dag, B_j = U). Each
 *   block contributes sum_{jk} d_l^-2 = 1 to ||D||_proj, so ||D||_proj = m (the
 *   NUMBER OF BLOCKS), NOT 1 — pushing the w' bound to O(m*delta) (.tex:1281; the
 *   number-of-blocks dependence tracked against .tex:484 by test_dhom T5).
 *
 * --------------------------------------------------------------------------
 * lem_approx — the Newton iteration (.tex:1277-1313).
 *   At step s (v^{(0)}=v, delta_0=delta): g = G_{v^{(s)}};
 *     w'(X)  = sum_terms v^{(s)}(A_term) * g(B_term, X)        (.tex:1279, star!)
 *     w''(X) = w'(X^dag)^dag                                    (.tex:1305)
 *     v^{(s+1)} = v^{(s)} + 1/2 (w' + w'')                      (.tex:1311)
 *   Each step drops the defect O(delta) -> O(delta^2 + eps) (.tex:1302). Iterate
 *   until the basis-sweep defect delta_s <= max(tol, c*eps); CAP iterations and
 *   FAIL LOUD if not converging (Rule 4). After each step the unit-preservation
 *   ||v^{(s)}(I_B) - I_A||_op and the involution symmetry
 *   ||v^{(s)}(E_{ba})^dag - v^{(s)}(E_{ab})||_op are asserted small.
 *
 *   CAVEAT (documented limitation). The basis-sweep defect
 *     delta = max over matrix-unit pairs (i,j) of ||G_v(E_i, E_j)||_op
 *   is a LOWER bound on the true sup-over-unit-ball delta (the true value can be
 *   up to a factor dim_B larger). It is used as the TERMINATION criterion only; a
 *   tighter HOPM upper bound is a later cycle (cf. aic_ecstar HOPM). The c_0 of
 *   cor_improvement (.tex:1317) is unstated in the paper — that is the NEXT module
 *   (errreduce aic-t81). Here we iterate to a documented stopping threshold.
 *
 * --------------------------------------------------------------------------
 * prop_delta_hominc — three precondition-guarded predicates (.tex:1194-1222).
 *   (i)   upper bound ||v|| <= 1 + O(delta + eps);
 *   (ii)  lower bound a >= 1 - O(delta + eps) IF the hypothesis a > 2 delta holds;
 *   (iii) unit-preservation ||v(I_B) - I_A|| <= O(delta + eps) (computed directly).
 *   The norms a, ||v|| here are BASIS-SWEEP estimates over the matrix-unit basis
 *   (not the true inf/sup over the unit ball); documented as such.
 */
#ifndef AIC_DHOM_H
#define AIC_DHOM_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_ecstar.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The genuine finite-dim C* algebra B = (+)_l M_{d[l]}. */
typedef struct {
    slong num_blocks; /* m                                                    */
    slong *d;         /* d[0..m-1] block dims                                  */
    slong *offset;    /* offset[l] = sum_{l'<l} d[l']^2  (matrix-unit index)   */
    slong dim_B;      /* sum_l d[l]^2                                          */
    slong n_B;        /* sum_l d[l]  (size of block-diagonal representative)   */
} aic_dhom_B;

/* A delta-homomorphism v : B -> A, stored by its values on the matrix-unit basis.
 *   B  : BORROWED domain C* algebra.
 *   A  : BORROWED codomain eps-C* algebra (its star is the multiplication in A).
 *   n  : A->n (ambient dim; every vE[i] is n x n).
 *   vE : dim_B operators, vE[i] = v(E_i) (the matrix-unit at linear index i). */
typedef struct {
    const aic_dhom_B  *B;
    const aic_ecstar  *A;
    slong              n;
    acb_mat_t         *vE; /* dim_B operators, each n x n */
} aic_dhom_v;

/* The explicit Pauli diagonal of B as a flat list of (weight, U) pairs (the
 * cross-term-free EMBEDDED SUM; see the .tex:1254 DISCREPANCY box in
 * aic_dhom_diag.c — NOT the paper's joint Cartesian product).
 *   nterms = sum_l d[l]^2 (a SUM over blocks; one term per (l,j,k) triple).
 *   w[t]   = d_{l(t)}^-2  (the per-block weight; depends only on the block l of t).
 *   U[t]   = Shat^{(l)}_{jk}, the d_l x d_l Pauli embedded in sector l (0 else) —
 *            a partial isometry, NOT the joint block-diagonal unitary.
 * The diagonal is D = sum_t w[t] U[t]^dag (x) U[t], i.e. A_t = w[t] U[t]^dag,
 * B_t = U[t] in the .tex:1247 notation. sum_t w[t] ||U[t]^dag|| ||U[t]|| = m
 * (the number of blocks), NOT 1. */
typedef struct {
    slong       nterms;
    arb_ptr     w;   /* nterms weights                              */
    acb_mat_t  *U;   /* nterms n_B x n_B unitaries                  */
} aic_dhom_diag;

/* ---- B algebra (src/aic_dhom_B.c) + Pauli diagonal (src/aic_dhom_diag.c) - */

/* Initialise B = (+)_l M_{dims[l]}. Copies dims (length m, each >= 1). Computes
 * offsets, dim_B, n_B. Free with aic_dhom_B_clear. */
void aic_dhom_B_init(aic_dhom_B *B, const slong *dims, slong m);
void aic_dhom_B_clear(aic_dhom_B *B);

/* Matrix-unit linear index mu(l,a,b) = offset[l] + a*d[l] + b. */
slong aic_dhom_B_index(const aic_dhom_B *B, slong l, slong a, slong b);

/* Block matmul Z = X . Y in B: X, Y, Z are n_B x n_B block-diagonal (caller
 * passes block-diagonal matrices; this multiplies block-by-block so the
 * off-block entries stay zero). Z != X, Z != Y required. */
void aic_dhom_B_mul(acb_mat_t Z, const aic_dhom_B *B, const acb_mat_t X,
                    const acb_mat_t Y, slong prec);

/* The block-diagonal n_B x n_B unit I_B = (+)_l I_{d[l]}. `out` init'd n_B x n_B. */
void aic_dhom_B_unit(acb_mat_t out, const aic_dhom_B *B);

/* The matrix-unit E_i (block-diagonal n_B x n_B, single 1 in block l at (a,b))
 * for linear index i. `out` init'd n_B x n_B. */
void aic_dhom_B_matunit(acb_mat_t out, const aic_dhom_B *B, slong i);

/* Generalized Pauli S_{jk} = X^j Z^k on C^d (.tex:1252):
 *   S[l,m] = exp(2 pi i k m / d) if l == (m+j) mod d, else 0.
 * `out` init'd d x d. Unitary. */
void aic_dhom_pauli(acb_mat_t out, slong d, slong j, slong k, slong prec);

/* Build the explicit Pauli diagonal of B (.tex:1250-1254). `out` is OUTPUT
 * (all fields allocated here). Free with aic_dhom_diag_clear. */
void aic_dhom_diag_build(aic_dhom_diag *out, const aic_dhom_B *B, slong prec);
void aic_dhom_diag_clear(aic_dhom_diag *out);

/* ---- v map + multiplicativity defect (src/aic_dhom_defect.c) ------------ */

/* Allocate vE[0..dim_B-1] (each n x n zeroed), set B, A, n. ASSERTS A->n == n.
 * Free with aic_dhom_v_clear. */
void aic_dhom_v_init(aic_dhom_v *v, const aic_dhom_B *B, const aic_ecstar *A);
void aic_dhom_v_clear(aic_dhom_v *v);

/* Apply v: out = v(X) = sum_i <coord_i(X)> vE[i], where coord_i(X) is the
 * matrix-unit entry of the block-diagonal X at linear index i. X is n_B x n_B;
 * out is n x n (caller-init'd). */
void aic_dhom_v_apply(acb_mat_t out, const aic_dhom_v *v, const acb_mat_t X,
                      slong prec);

/* The multiplicativity defect g = G_v(X,Y) = v(XY) - v(X) * v(Y) (.tex:1258),
 * where XY is B's block product and * is the STAR in A (aic_ecstar_star — NEVER
 * plain matmul). X, Y are n_B x n_B (block-diagonal); out is n x n. */
void aic_dhom_Gv(acb_mat_t out, const aic_dhom_v *v, const acb_mat_t X,
                 const acb_mat_t Y, slong prec);

/* Basis-sweep multiplicativity defect: max over matrix-unit pairs (i,j) of
 * ||G_v(E_i, E_j)||_op, as a certified arb ball (the LOWER-bound termination
 * metric; see header CAVEAT). `out` init'd. */
void aic_dhom_defect_sweep(arb_t out, const aic_dhom_v *v, slong prec);

/* prop_delta_hominc (.tex:1194). Computes the three quantities as certified arb
 * balls (any output may be NULL to skip):
 *   norm_ub : basis-sweep ||v|| = max_i ||vE[i]||_op  (predicate (i): <= 1+O).
 *   norm_lb : basis-sweep a    = min_i ||vE[i]||_op  (predicate (ii): >= 1-O when
 *             a > 2 delta; the routine does NOT enforce the hypothesis, it reports
 *             a so the caller can guard).
 *   unit_def: ||v(I_B) - I_A||_op  (predicate (iii): <= O(delta+eps)).
 * The norms are basis-sweep estimates (documented limitation). */
void aic_dhom_prop_bounds(arb_t norm_ub, arb_t norm_lb, arb_t unit_def,
                          const aic_dhom_v *v, slong prec);

/* ---- lem_approx Newton iteration (src/aic_dhom_approx.c) ---------------- */

/* Result of one Newton sweep / the full lem_approx run. */
typedef struct {
    slong  steps;       /* Newton steps taken                                 */
    double delta0;      /* initial basis-sweep defect                         */
    double delta_final; /* final basis-sweep defect                           */
} aic_dhom_approx_stats;

/* One Newton step IN PLACE: v <- v + 1/2 (w' + w'') (.tex:1311), where
 *   w'(X)  = sum_terms v(A_term) * g(B_term, X)  (g = G_v at entry; .tex:1279)
 *   w''(X) = w'(X^dag)^dag                        (.tex:1305)
 * D is the Pauli diagonal of B. Updates each vE[i] += 1/2(w'(E_i)+w''(E_i)).
 * (Exposed for the Newton-quadratic-convergence test; lem_approx loops it.) */
void aic_dhom_newton_step(aic_dhom_v *v, const aic_dhom_diag *D, slong prec);

/* lem_approx (.tex:1224): iterate aic_dhom_newton_step until the basis-sweep
 * defect <= max(tol_abs, eps_target) or max_steps is hit. FAILS LOUD (Rule 4) if
 * the cap is hit without the defect dropping below the threshold, AND if the
 * defect ever INCREASES across a step (Newton must contract). After every step
 * asserts unit-preservation and involution-symmetry are within unit_tol. Builds
 * the Pauli diagonal internally. Fills *st if non-NULL.
 *   eps_target : the O(eps) floor (pass 0 for the eta=0 exact-hom case, where the
 *                defect drops to machine-zero).
 *   tol_abs    : absolute machine floor added to eps_target (e.g. 1e-12 at p=53).
 *   unit_tol   : tolerance for the per-step unit/involution asserts. */
void aic_dhom_approx(aic_dhom_v *v, double eps_target, double tol_abs,
                     double unit_tol, slong max_steps, slong prec,
                     aic_dhom_approx_stats *st);

#ifdef __cplusplus
}
#endif

#endif /* AIC_DHOM_H */
