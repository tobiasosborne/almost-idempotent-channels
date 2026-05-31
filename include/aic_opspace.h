/* aic_opspace.h — constructive proof of th_main_ext (approximate_algebras.tex
 * :1447-1561, §10): the O(eps)-isomorphism v: B -> A from aic_cstar_build is an
 * EXTENDED delta-isomorphism — every ampliation 1_{M_n} (x) v : M_n (x) B ->
 * M_n (x) A is itself a delta-isomorphism, with the SAME delta = O(eps), the
 * constant independent of the level n (th_main_ext, .tex:1538-1540). Module bead
 * aic-zwo (increment O1, the operator-norm ampliation machinery + structural
 * verification). Design contract: docs/research/opspace_design.md (read its
 * ORCHESTRATOR CORRECTION) + opspace_paper_leg.md + opspace_web_leg.md.
 *
 * --------------------------------------------------------------------------
 * THE OPERATOR-NORM cb-INCLUSION, NOT THE FROBENIUS sigma_min (FINDINGS §C12 —
 * the binding landmine of this module). th_main_ext's content is the
 * OPERATOR-norm cb-inclusion lower bound (.tex:1489, def:opspace .tex:1453-1464):
 *
 *     a_n = inf_{X != 0, X in M_n (x) B} ||(1_{M_n} (x) v)(X)||_op / ||X||_op
 *
 * with BOTH norms the OPERATOR norm (largest singular value of the ampliated
 * matrix). The non-trivial prop_inc_ext induction a_{2n} >= a_n/2 (.tex:1493-1503)
 * is a bound on THIS operator-norm quantity, and a_n >= 1 - delta' uniformly in n
 * (.tex:1505) is th_main_ext.
 *
 * The Frobenius COORDINATE sigma_min route (opspace_design.md §3.2, the route
 * aic_dhom_v_sigma_min uses, §C6) is VACUOUS HERE: the ampliated coordinate
 * matrix is I_{n^2} (x) M_1, whose singular values are those of M_1 repeated n^2
 * times, so sigma_min(I_{n^2} (x) M_1) = sigma_min(M_1) for ANY linear v,
 * INDEPENDENT of the ampliation. A Frobenius-sigma_min "universality canary" for
 * th_main_ext cannot fail (CLAUDE.md Rule 5) — it is exactly the "test that
 * cannot fail" this module exists to avoid. So opspace measures the OPERATOR-norm
 * stretch (a max-stretch HOPM over the operator-norm unit ball of the C* algebra
 * M_n (x) B), NOT the Frobenius sigma_min. The non-vacuity tooth in
 * tests/test_opspace.c demonstrates the distinction concretely.
 *
 * --------------------------------------------------------------------------
 * SMITH TRUNCATION (FINDINGS §D3, RESOLVED). The "for all n" quantifier is
 * answered by Smith's lemma (R.R. Smith, J. London Math. Soc. (2) 27 (1983)
 * 157-166; Pisier, Intro to Operator Space Theory (2003) Prop 1.12 p.26; Watrous,
 * TQI (2018) Thm 3.46 + Cor 3.47 via adjoint duality): for ANY linear map
 * u: E -> M_N from an operator space into the matrix algebra M_N (no CP / C*
 * hypothesis), ||u||_cb = ||1_{M_N} (x) u||_op — the cb-norm is ATTAINED at
 * ampliation level N (the codomain ambient dimension) and is not increased by any
 * larger ampliation. Applied here:
 *   forward v: B -> A <= M_N (N = v->n):   ||v||_cb   attained at level N_max = N.
 *   inverse v^{-1}: A -> B <= M_{n_B}:     ||v^{-1}||_cb attained at level n_B.
 * Since v is bijective, a_n = 1/||(1_{M_n} (x) v)^{-1}||_op, and its all-n value
 * = 1/||v^{-1}||_cb is attained at n = n_B = sum_l d_l (the size of B). So finite
 * truncation suffices; no infinite verification problem.
 *
 * --------------------------------------------------------------------------
 * v IS REUSED UNCHANGED (.tex:1542, 1557). th_main_ext does NOT recompute B or v;
 * the §9 master-loop output (aic_cstar_build) is the same object, and §10 is a
 * post-hoc CERTIFICATION argument over it. opspace is a verification module: it
 * takes the existing aic_dhom_v and measures the operator-norm ampliated stretch.
 *
 * --------------------------------------------------------------------------
 * HONESTY ON BOUND DIRECTIONS (O1 vs O2). The HOPM max-stretch is a LOWER bound on
 * the operator-norm max-stretch ||f_n||_op: it finds a good-but-maybe-suboptimal X
 * over the op-norm unit ball and certifies the witness ratio. Therefore:
 *   - the forward HOPM gives a LOWER bound on ||v_n||_op (so it CANNOT by itself
 *     certify ||v||_cb <= 1+O(eps); that needs an UPPER bound -> the Watrous SDP,
 *     the SEPARATE O2 increment, bead filed by the orchestrator).
 *   - a_n = 1/||v_n^{-1}||_op with the HOPM lower-bounding ||v_n^{-1}||_op, so the
 *     reported a_cb is an UPPER bound on the true lower-isometry constant a_n.
 * What O1 DELIVERS FAITHFULLY (all operator-norm, non-vacuous, exactly what §C12
 * demands):
 *   (a) the eta=0 oracle: v an exact *-iso => a COMPLETE ISOMETRY, so the stretch
 *       is EXACTLY 1 for every X and the HOPM lower bound == exact == 1, at every
 *       level n;
 *   (b) the prop_inc_ext doubling structure on the MEASURED operator-norm a_n
 *       (a_{2n} >= a_n/2, .tex:1493-1503);
 *   (c) the operator-norm UNIVERSALITY CANARY: the measured constant does NOT grow
 *       with dim A nor with the ampliation level n.
 * The CERTIFIED cb-bound (the SDP UPPER bound, ||v||_cb <= 1+O(eps)) is O2.
 *
 * --------------------------------------------------------------------------
 * PRECISION POSTURE (mirrors aic_dhom_v_sigma_min, §C6). The HOPM runs in the fast
 * DOUBLE path (LAPACK SVD via aic_latd, the aic_ecstar HOPM kernel pattern); the
 * reported balls are coarse fail-loud GATE midpoints (zero-radius arb of the
 * double value), like the projection-nontriviality and sigma_min gates. The HOPM
 * KERNEL ITSELF IS PURE DOUBLE — `out` is arb_set_d(double), a zero-radius wrap,
 * NOT an arb re-evaluation of the witness ratio. The only prec-dependence is the
 * acb assembly + inversion of the coordinate matrix M_1 / M_1^{-1} (build_M1, the
 * acb_mat_solve). So the double-vs-arb@53 cross-check (test_opspace T-cross,
 * cross-check ladder rung 2) tests the PRECISION OF THAT ASSEMBLY + the
 * dot-product rounding stability (~1e-15 between prec 53 and 256), a COARSE GATE
 * consistent with aic_dhom_v_sigma_min's posture — NOT a genuine arb-vs-double
 * ALGORITHM cross-check (which would need an arb HOPM). A certified arb enclosure
 * of the worst-case stretch defers to the SDP (O2) / aic-0at, exactly as the
 * ecstar HOPM lower bound does.
 */
#ifndef AIC_OPSPACE_H
#define AIC_OPSPACE_H

#include <flint/arb.h>

#include "aic_dhom.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The extended-cb-isomorphism certificate of v: B -> A (from aic_cstar_build).
 *   v          : BORROWED extended delta-isomorphism (must stay alive while this
 *                struct's arb fields are read; the struct does not own v).
 *   a_cb       : the measured operator-norm lower-isometry constant at the inverse
 *                Smith level n_B = sum_l d_l: a_cb = 1/||v_{n_B}^{-1}||_op (HOPM).
 *                An UPPER bound on the true a_n (HOPM lower-bounds ||v^{-1}_n||).
 *                Reported as a zero-radius arb GATE midpoint.
 *   cb_forward : the measured forward operator-norm max-stretch ||v_N||_op at the
 *                forward Smith level N = v->n (HOPM). A LOWER bound on ||v||_cb.
 *                Reported as a zero-radius arb GATE midpoint.
 *   a_cb_flat  : the operator-norm universality canary metric over the levels run:
 *                max_n a_n / min_n a_n (a no-upward-trend check; ~1 for the
 *                I_{n^2}-style isometric ampliation, FINDINGS §C12). Stored double.
 *   N_max      : the highest ampliation level certified (the Smith truncation).
 *   n_B        : the INVERSE Smith level n_B = sum_l d_l (the size of B). The all-n
 *                value of a_cb is attained at this level (Smith, FINDINGS §D3). A
 *                first-class field so factorize (D4, the §12 decode) can size the
 *                per-block 1-designs + the W_j Choi blocks from one struct (see
 *                docs/research/factorize_d4_research.md §4 add #1). Equals v->B->n_B.
 */
typedef struct {
    const aic_dhom_v *v;
    arb_t  a_cb;
    arb_t  cb_forward;
    double a_cb_flat;
    slong  N_max;
    slong  n_B;
    /* O2 (bead aic-pjr): the CERTIFIED cb-norm UPPER bounds from the rect Watrous
     * SDP-restoration certifier (aic_opspace_certify_cb_upper). Together with the
     * O1 HOPM lower bounds (cb_forward, 1/a_cb) they BRACKET the true cb-norms:
     *   cb_forward (HOPM lower)  <= ||v||_cb     <= cb_upper    (SDP, certified),
     *   1/a_cb     (HOPM lower)  <= ||v^{-1}||_cb <= cbinv_upper (SDP, certified).
     * aic_opspace_certify (O1) initialises these to +inf (the sentinel that flags a
     * struct on which the O2 step has NOT been run); aic_opspace_certify_cb_upper
     * fills them with the rigorous bound. */
    arb_t  cb_upper;
    arb_t  cbinv_upper;
} aic_opspace_result;

/* ---- O1 core: ampliated operator-norm max-stretch (src/aic_opspace_ampliate.c) -
 *
 * Operator-norm max-stretch of the ampliated FORWARD map 1_{M_n} (x) v at level n:
 *   ||v_n||_op = sup_{X != 0 in M_n (x) B} ||(1_{M_n} (x) v)(X)||_op / ||X||_op,
 * computed by a scale-invariant HOPM over the operator-norm unit ball of the C*
 * algebra M_n (x) B (the maximizer of Re<G,X> over that ball is the polar factor
 * of the gradient G = v_n^dag(u w^*) taken WITHIN M_n (x) B; mirror ecstar's
 * polar-then-project accept guard). A rigorous LOWER bound on the true sup. `out`
 * (init'd arb_t) receives the double HOPM value as a zero-radius gate midpoint.
 * ASSERTS n >= 1, v bijective (dim_A == dim_B). */
void aic_opspace_forward_stretch(arb_t out, const aic_dhom_v *v, slong n,
                                 slong prec);

/* Operator-norm max-stretch of the ampliated INVERSE map 1_{M_n} (x) v^{-1} at
 * level n: ||v_n^{-1}||_op (HOPM over the op-norm unit ball of M_n (x) A, the
 * eps-C* subspace — so the polar is projected onto M_n (x) A with the same
 * monotone-accept guard). v^{-1} is built by inverting the coordinate matrix
 * M_1[k,i] = <B_k, vE[i]>_F (square + invertible in the bijective case; reuses the
 * assembly of aic_dhom_v_sigma_min then a LAPACK inverse). A rigorous LOWER bound
 * on ||v_n^{-1}||_op, so 1/out is an UPPER bound on the lower-isometry a_n. `out`
 * init'd. ASSERTS n >= 1, v bijective. */
void aic_opspace_inverse_stretch(arb_t out, const aic_dhom_v *v, slong n,
                                 slong prec);

/* ---- structural ampliation 1_{M_n} (x) v (src/aic_opspace_map.c) -------- *
 *
 * Apply the FORWARD ampliation Y = (1_{M_n} (x) v)(X) directly, where X is an
 * (n*n_B) x (n*n_B) block matrix (block (k,l) = X_{kl} in B, in B's matrix-unit
 * coordinates) and Y is the (n*N) x (n*N) result (block (k,l) = v(X_{kl}) in A,
 * ambient M_N). Both are caller-allocated row-major double _Complex buffers
 * (Y length (n*N)^2, X length (n*n_B)^2). This drives the SAME aic_opmap_apply_fn
 * the HOPM kernel uses (aic_opspace_apply.c) — it is the production ampliation
 * code path, exposed so the structural-ampliation cross-check tooth
 * (test_opspace.c) can pin it against an INDEPENDENT Kronecker block assembly at
 * n >= 2 (FINDINGS §C12, the 2026-05-31 hostile-review landmine: a (0,0)-block-only
 * cripple of apply_fn passed every level-independent fixture; this entry + the
 * tooth catch it). ASSERTS n >= 1, v bijective. `prec` is the M_1 assembly prec. */
void aic_opspace_ampliate_forward(double _Complex *Y, const aic_dhom_v *v,
                                  const double _Complex *X, slong n, slong prec);

/* ---- public v^{-1} builder for factorize / D4 (src/aic_opspace_map.c) --- *
 *
 * v^{-1}: A -> B as dim_A operators vinvB[k] = v^{-1}(B_k), each an n_B x n_B
 * block-diagonal matrix in B (built by inverting the coordinate matrix
 * M_1[k,i] = <B_k, vE[i]>_F and reading the columns back through B's matrix-unit
 * basis: v^{-1}(B_k) = sum_i (M_1^{-1})[i,k] E_i). The SAME inverse the inverse
 * stretch certifies, lifted to a reusable object so factorize's Upsilon~ =
 * v^{-1} o Phi~ uses the certified-and-used inverse, not a re-derivation (the
 * highest-value D4 interface add, docs/research/factorize_d4_research.md §4 #3).
 *   dim_A   : OUT, set to dim_A == dim_B (the count of vinvB operators).
 *   n_B     : OUT, set to v->B->n_B (the ambient size of each vinvB[k]).
 *   vinvB   : OUT, *vinvB allocated here as dim_A acb_mat_t (each n_B x n_B);
 *             v^{-1}(B_k) in vinvB[k]. Free with aic_opspace_vinv_clear.
 * ASSERTS v bijective (M_1 square + invertible; fail loud on singularity). */
void aic_opspace_build_vinv(acb_mat_t **vinvB, slong *dim_A, slong *n_B,
                            const aic_dhom_v *v, slong prec);

/* Free the dim_A operators returned by aic_opspace_build_vinv. */
void aic_opspace_vinv_clear(acb_mat_t *vinvB, slong dim_A);

/* ---- O2.1: ADJOINT Choi assemblers for the cb-norm SDP (src/aic_opspace_choi.c) -
 *
 * Build the Convention-A Choi matrices of the ADJOINT maps v* and (v^{-1})* for the
 * Watrous diamond-norm certifier (O2, bead aic-pjr). We certify the cb-SPECTRAL norm
 * of the th_main_ext isomorphism v: B -> A (.tex:1538-1540) via the adjoint duality
 *   ||v||_cb = ||v*||_⋄,    ||v^{-1}||_cb = ||(v^{-1})*||_⋄
 * (docs/research/opspace_o2_design.md §0.5 PINNED CONVENTION — the empirically-pinned
 * recipe that SUPERSEDES §1.4/§2.4/§2.5).
 *
 * THE GOLDEN RULE (design §0.5). For f: M_in -> M_out, the Convention-A Choi is
 *   J(f)[s*out + i, t*out + j] = f(E_st)[i,j]   (INPUT s,t MAJOR/stride out;
 *                                                OUTPUT i,j MINOR),
 * fed to the rect diamond-norm SDP as (d_maj = in, d_min = out). Then raw optval =
 * ||f||_⋄ EXACTLY — normalization factor 1, NO 2/n. The SDP duals trace the MINOR
 * factor (tr_sys = 2, size d_min); the primal places the density on the MAJOR factor
 * (rho_on = :major, size d_maj). The two assemblers below build v* / (v^{-1})*
 * DIRECTLY in this INPUT-MAJOR convention (NOT by transposing J(v): the full
 * transpose keeps the wrong [out,in] block layout — §0.5).
 *
 * BOTH maps are HERMITICITY-PRESERVING (v, v^{-1} are HP), so their adjoint Choi are
 * HERMITIAN (but indefinite — the maps are not CP). */

/* J(v*) for v*: M_N -> M_{n_B}  (d_maj = in = N, d_min = out = n_B).
 *   v*(E_ab^{(N)}) = sum_i conj(vE[i][a,b]) * E_i  (HS-adjoint of v; E_i the n_B-rep
 *   matrix unit, a single 1 at aic_opspace_b_unit_pos(B,i) = (r,c)).
 * So Jvs[a*n_B + r, b*n_B + c] = conj(vE[i][a,b]) for each i (a,b INPUT/major, (r,c)
 * OUTPUT/minor). `Jvs` is caller-init'd (N*n_B) x (N*n_B) (shape ASSERTED). */
void aic_opspace_choi_vstar(acb_mat_t Jvs, const aic_dhom_v *v, slong prec);

/* J((v^{-1})*) for (v^{-1})*: M_{n_B} -> M_N  (d_maj = in = n_B, d_min = out = N).
 *   (v^{-1})*(E_ab^{(n_B)}) = sum_k conj(vinvB[k][a,b]) * B_k, where
 *   vinvB[k] = v^{-1}(B_k) (aic_opspace_build_vinv) and B_k = v->A->B[k] (A's
 *   Frobenius-ON basis, each N x N).
 * So Jvis[a*N + r, b*N + s] += conj(vinvB[k][a,b]) * B_k[r,s] summed over k (a,b
 * INPUT/major in [0,n_B), (r,s) OUTPUT/minor in [0,N)). `Jvis` is caller-init'd
 * (n_B*N) x (n_B*N) (shape ASSERTED). ASSERTS v bijective (dim_A == dim_B). */
void aic_opspace_choi_vinvstar(acb_mat_t Jvis, const aic_dhom_v *v, slong prec);

/* ---- O1 certification pipeline (src/aic_opspace_cert.c) ------------------ *
 *
 * Run the ampliation levels and fill r. Computes:
 *   - the forward stretch ||v_N||_op at the forward Smith level N = v->n,
 *   - the lower isometry a_n = 1/||v_n^{-1}||_op at the doubling levels
 *     n = 1, 2, 4, ..., up to the inverse Smith level n_B = sum_l d_l,
 *   - a_cb = a_{n_B}, the all-n value (Smith), and the universality-canary
 *     flatness a_cb_flat over those levels.
 * GUARDS (fail loud, Rule 4):
 *   - v bijective (dim_A == dim_B);
 *   - the prop_inc_ext precondition a_1 > 2*delta (.tex:1505); abort with the
 *     measured a_1 and delta if violated (a genuine algorithmic failure / stop
 *     condition, not a truncation artifact);
 *   - the doubling step a_{2n} >= a_n/2 (.tex:1503) at each consecutive pair;
 *     abort if violated (the prop_inc_ext mechanism broke).
 * `delta` is the isomorphism defect of v (iso_def from aic_cstar_build / the O(eps)
 * scale). `prec` is the arb working precision. Free with aic_opspace_result_clear. */
void aic_opspace_certify(aic_opspace_result *r, const aic_dhom_v *v,
                         double delta, slong prec);

/* ---- O2 payoff: certified cb-norm UPPER bounds (src/aic_opspace_o2.c) ---- *
 *
 * Fill r->cb_upper and r->cbinv_upper with the RIGOROUS Watrous-SDP upper bounds
 * on ||v||_cb and ||v^{-1}||_cb (the certified half of th_main_ext, .tex:1538-1540).
 * `r` must already have been filled by aic_opspace_certify (the O1 HOPM lower
 * bounds cb_forward / a_cb are read for the bracket guard). The committed MIN-dual
 * feasible points are passed in:
 *   (Y0_fwd, Y1_fwd) for J(v*)        (d_maj = N,   d_min = n_B),
 *   (Y0_inv, Y1_inv) for J((v^{-1})*) (d_maj = n_B, d_min = N).
 * The Choi matrices J(v*), J((v^{-1})*) are assembled internally
 * (aic_opspace_choi_vstar / _vinvstar), then aic_cbnorm_certify_rect_upper restores
 * each dual point and reads off the bound (design §0.5: factor 1, tr_sys=2 minor,
 * eps*d_min shift). GUARD (fail loud, Rule 4): the bracket O1.HOPM_lower <= O2.SDP_
 * upper must hold for both directions (a HOPM > SDP would be a soundness bug);
 * abort with the values if violated (a small slack covers the arb radius + the
 * double-HOPM gate). ASSERTS v bijective. */
void aic_opspace_certify_cb_upper(aic_opspace_result *r, const aic_dhom_v *v,
                                  const acb_mat_t Y0_fwd, const acb_mat_t Y1_fwd,
                                  const acb_mat_t Y0_inv, const acb_mat_t Y1_inv,
                                  slong prec);

void aic_opspace_result_clear(aic_opspace_result *r);

#ifdef __cplusplus
}
#endif

#endif /* AIC_OPSPACE_H */
