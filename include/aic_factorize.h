/* aic_factorize.h — the approximate factorization of an almost-idempotent UCP
 * channel through a genuine C* algebra (bead aic-tff, module "factorize"),
 * Increment F1. Realizes the NON-UCP "tilde" maps of th_factorization
 * (approximate_algebras.tex:2730-2899), specifically the outline at .tex:2749-2752.
 * Design / feasibility: docs/research/factorize_d4_research.md §1 (Steps 1-3) + §6
 * (the F1-F4 increment skeleton — this file is F1). FINDINGS §D4 (BUILDABLE-MODULO),
 * §C2 (the star), §C3 (A = Img Phi~ is OBLIQUE).
 *
 * WHAT F1 REALIZES (.tex:2749-2752, verbatim in src/aic_factorize.c). By
 * th_main_ext there is a genuine finite-dim C* algebra B and an extended
 * O(eta)-isomorphism v: B -> A = Img Phi~ (aic_cstar_build + opspace). F1 builds
 * the two maps the paper REVERSES out of the u-isomorphism corollary:
 *
 *   Delta~   = iota o v       : B -> B(H)   (.tex:2749, "v followed by A -> B(H)")
 *   Upsilon~ = v^{-1} o Phi~  : B(H) -> B   (.tex:2749, "Phi~ with target A, then v^{-1}")
 *
 * iota: A -> B(H) is the inclusion (a NO-OP on matrices: A <= M_N), so
 *   Delta~(X) = v(X) = sum_i x_i vE[i]    (literally aic_dhom_v_apply(v, X)),
 * viewing the N x N output as living in B(H) rather than A (X is the n_B x n_B
 * block-diagonal representative of an element of B). And
 *   Upsilon~(W) = v^{-1}(Phi~(W)) = sum_k <B_k, Phi~(W)>_F vinvB[k],
 * where Phi~(W) is the OBLIQUE projection of W onto A (the superop S_tilde applied,
 * NOT the Frobenius Pi_A — FINDINGS §C3), {B_k} = A's Frobenius-orthonormal basis
 * (Aec->A.B[k]), <B_k,Z>_F = sum_ij conj(B_k[i,j]) Z[i,j], and
 * vinvB[k] = v^{-1}(B_k) from aic_opspace_build_vinv.
 *
 * THE EXACT-BY-CONSTRUCTION IDENTITIES (.tex:2751, tilde_DelUps):
 *   Delta~ Upsilon~ = Phi~   (as maps M_N -> M_N: v o v^{-1} = id_A, Phi~ -> A),
 *   Upsilon~ Delta~ = 1_B    (as maps B -> B:    Phi~|_A = id, v^{-1} o v = id_B).
 * They hold to the numerical precision of the v-inversion + S_tilde idempotence;
 * tests/test_factorize.c MEASURES them (T1 eta=0 oracle vs the ORIGINAL channel
 * Phi; T2 eta>0 vs Phi~ = S_tilde, the OBLIQUE star path teeth — §C2/§C3).
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The paper's tilde maps
 * are abstract: v is the th_main_ext iso "that exists", Phi~ the cb-algebra functional
 * calculus theta(2 Phi - 1). In finite dim every piece is an explicit matrix object
 * already built and certified by a closed module: v (aic_cstar_build), v^{-1}
 * (aic_opspace_build_vinv, the SAME inverse the inverse stretch certifies), Phi~ =
 * S_tilde (aic_assoc_regularize). F1 is pure linear-algebra glue over them; the
 * paper's exact identities are the spec the tests assert.
 *
 * --------------------------------------------------------------------------
 * INCREMENT F2 — the UCP ENCODE map Delta (.tex:2771-2801). EXTENDS this module.
 *
 * F1's Delta~ = iota o v is CP-blind (NOT positive: v(Y^dag Y) need not be PSD in
 * B(H), only star-positive in A). Step 4 of th_factorization CP-izes it by a
 * unitary 1-design average through the ORIGINAL channel Phi (Phi CP is WHAT makes
 * the average manifestly CP, .tex:2791-2795):
 *
 *   Delta'(X) = sum_s p_s Phi( Delta~(X U_s^dag) Delta~(U_s) )    (.tex:2788)
 *
 * where Delta~(.)Delta~(.) is the ORDINARY M_N product (NOT the star), Phi is the
 * original UCP map (aic_ucp_apply over F->phi), and {(p_s, U_s)} is the unitary
 * 1-design DIAGONAL of B = (+)_j M_{d_j}: the CARTESIAN PRODUCT of the per-block
 * generalized-Pauli designs (.tex:2780-2784). Per block j, U_{j,(a,b)} = S_ab =
 * X^a Z^b (a,b in [0,d_j); aic_dhom_pauli), p_{j,(a,b)} = d_j^{-2}; the joint index
 * s = (s_1,...,s_m) has weight p_s = prod_l d_l^{-2} and JOINT unitary
 * U_s = U_{1,s_1} (+) ... (+) U_{m,s_m} (the n_B x n_B block-diagonal join, a
 * genuine unitary in B). Total terms = prod_l d_l^2 (M_2 -> 4; M_2(+)M_2 -> 16).
 *
 * THE §A2 DISTINCTION (FINDINGS §A2 — the landmine). F2's U_s is the GENUINE
 * per-block Pauli S_ab joined BLOCK-DIAGONALLY (.tex:2782-2783) — a real unitary
 * in B. This is DIFFERENT from aic_dhom_diag_build's output (the cross-term-free
 * EMBEDDED partial-isometry sum D = sum_l D_l, for lem_approx's error reduction).
 * F2 does NOT use aic_dhom_diag_build; it builds U_s itself as the Cartesian-
 * product block-diagonal join of aic_dhom_pauli outputs (each U_s asserted unitary).
 *
 * Then unitalize (.tex:2799):  Delta(X) = Delta'(I_B)^{-1/2} Delta'(X) Delta'(I_B)^{-1/2},
 * UCP (unital by construction, CP because congruence preserves CP). The build
 * ASSERTS ||Delta'(I_B) - 1_H||_op < 1 (Rule 4, shard H #7) before the inverse-sqrt.
 *
 * NON-SCOPE (later increments). The UCP DECODE Upsilon (lem_RC, Step 5) is F3;
 * the certified cb-norm UPPER bound ||Delta~||_cb, ||Upsilon~||_cb <= 1+O(eta) is
 * opspace O2 (structurally ||v||_cb / ||v^{-1}||_cb); the end-to-end
 * DelUps/UpsDel2 + duals are F4.
 */
#ifndef AIC_FACTORIZE_H
#define AIC_FACTORIZE_H

#include <flint/acb_mat.h>

#include "aic_assoc.h"
#include "aic_dhom.h"
#include "aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The factorization handle (grows across increments; F1 + F2). BORROWS v (B -> A,
 * from aic_cstar_build), the aic_assoc_ecstar Aec (A = Img Phi~ + S_tilde + A's
 * Frobenius-ON basis B[k]), and phi (the ORIGINAL UCP channel Phi, used by F2's
 * Delta' — Phi being CP is what makes Delta' manifestly CP, .tex:2791). OWNS the
 * vinvB array (v^{-1}(B_k), from aic_opspace_build_vinv) and, after
 * aic_factorize_delta_build, the precomputed deltaI_invsqrt = Delta'(I_B)^{-1/2}
 * (F2). v, Aec, phi must stay alive for this struct's lifetime. Free with
 * aic_factorize_clear. */
typedef struct {
    slong N;                       /* dim H = v->n = A->n                          */
    slong n_B;                     /* B's block-diagonal representative size        */
    slong dim_A;                   /* dim A = dim B                                 */
    const aic_dhom_v       *v;     /* BORROWED: B -> A (aic_cstar_build)            */
    const aic_assoc_ecstar *Aec;   /* BORROWED: A = Img Phi~ + S_tilde + basis B[k] */
    const aic_ucp_kraus    *phi;   /* BORROWED: the ORIGINAL UCP Phi (F2 Delta')    */
    acb_mat_t *vinvB;              /* OWNED: dim_A operators v^{-1}(B_k), n_B x n_B */
    int        delta_ready;        /* 1 after aic_factorize_delta_build ran (F2)    */
    acb_mat_t  deltaI_invsqrt;     /* OWNED (F2): Delta'(I_B)^{-1/2}, N x N          */

    /* ---- F3 (the UCP DECODE map Upsilon via lem_RC; .tex:2840-2899) ---------
     * Filled by aic_factorize_upsilon_build (REQUIRES delta_ready). All OWNED;
     * arrays length m = v->B->num_blocks (one per B-block j). */
    int        upsilon_ready;      /* 1 after aic_factorize_upsilon_build ran       */
    slong      r;                  /* Phi's Stinespring/ancilla F dim (= phi->r)    */
    acb_mat_t  V;                  /* OWNED (F3): Phi's Stinespring V: H->H(x)F,    *
                                    * (N*r) x N, V[a*N+i,j]=K_a[i,j] (aic_ucp.h)    */
    slong     *e;                  /* OWNED: e[j] = E_j (per-block Stinespring rank) */
    acb_mat_t *W;                  /* OWNED: W[j] = W_j: H -> L_j(x)E_j,            *
                                    * (d_j*e_j) x N, L_j-major (D2)                  */
    acb_mat_t *C;                  /* OWNED: C[j] = C_j in B(E_j), e_j x e_j (D1)    */
    acb_mat_t *xi;                 /* OWNED: xi[j] = top RIGHT sing. vec of C_j,    *
                                    * e_j x 1 unit vector (D4)                       */
    acb_mat_t *L;                  /* OWNED: L[j] = L_j: L_j -> H(x)F, (N*r) x d_j  */
    acb_mat_t  upsI_invsqrt;       /* OWNED: Upsilon'(1_H)^{-1/2}, n_B x n_B (D6)    */
} aic_factorize;

/* Build `out` from a bijective extended O(eta)-isomorphism v: B -> A, the
 * associated eps-C* algebra Aec (A = Img Phi~), and the ORIGINAL UCP channel phi
 * (Phi: B(H) -> B(H), a self-map dim_K == dim_H == N; BORROWED, used by F2's
 * Delta'). Fills N, n_B, dim_A; builds vinvB via aic_opspace_build_vinv; leaves
 * delta_ready = 0 (call aic_factorize_delta_build for F2). ASSERTS (fail loud,
 * Rule 4): v->A == &Aec->A (same A object — the iso's codomain IS Aec's A), v
 * bijective (dim_A == dim_B), v->n == Aec->A.n, phi->dim_K == phi->dim_H == N.
 * Free with aic_factorize_clear. */
void aic_factorize_build(aic_factorize *out, const aic_dhom_v *v,
                         const aic_assoc_ecstar *Aec, const aic_ucp_kraus *phi,
                         slong prec);

/* Free the OWNED vinvB array and zero the struct. Does NOT touch the borrowed v
 * or Aec. */
void aic_factorize_clear(aic_factorize *out);

/* Delta~(X) = iota(v(X)) = v(X) = sum_i x_i vE[i] (.tex:2749). X is the n_B x n_B
 * block-diagonal representative of an element of B; `out` is the N x N ambient
 * image (caller-init'd N x N). A thin wrapper over aic_dhom_v_apply: the output
 * lives in B(H) rather than A (iota is a no-op on matrices). */
void aic_factorize_delta_tilde(acb_mat_t out, const aic_factorize *F,
                               const acb_mat_t X, slong prec);

/* Upsilon~(W) = v^{-1}(Phi~(W)) = sum_k <B_k, Phi~(W)>_F vinvB[k] (.tex:2749).
 * W is N x N (an operator in B(H)); `out` is the n_B x n_B block-diagonal
 * B-representative (caller-init'd n_B x n_B). Phi~(W) is the OBLIQUE projection of
 * W onto A via the superop S_tilde (aic_assoc_superop_apply over Aec->S_tilde —
 * NOT the Frobenius Pi_A, FINDINGS §C3), then the v^{-1} contraction over A's
 * Frobenius-ON basis {B_k} = Aec->A.B[k]. */
void aic_factorize_upsilon_tilde(acb_mat_t out, const aic_factorize *F,
                                 const acb_mat_t W, slong prec);

/* ---- Increment F2 — the UCP encode map Delta (src/aic_factorize_delta.c) --- */

/* The joint 1-design unitary U_s (.tex:2782-2783): the CARTESIAN-PRODUCT block-
 * diagonal join U_s = U_{1,s_1} (+) ... (+) U_{m,s_m} of per-block generalized
 * Paulis, where the flat term index `s` in [0, prod_l d_l^2) decomposes (mixed-
 * radix, block l with radix d_l^2; within a block (a,b) = (q / d_l, q % d_l) ->
 * S_ab) into one Pauli per block. `out` (caller-init'd n_B x n_B) is written as a
 * genuine block-diagonal unitary in B. nterms() = prod_l d_l^2. ASSERTS U_s
 * unitary (U_s^dag U_s = I_{n_B}) — fail loud (Rule 4, the §A2 sanity gate). */
slong aic_factorize_design_nterms(const aic_factorize *F);
void  aic_factorize_design_unitary(acb_mat_t out, const aic_factorize *F,
                                   slong s, slong prec);

/* Delta'(X) = sum_s p_s Phi( Delta~(X U_s^dag) Delta~(U_s) ) (.tex:2788), the
 * manifestly-CP 1-design average. X is n_B x n_B (B's block-diagonal rep); `out`
 * is N x N (caller-init'd). Phi is the ORIGINAL UCP map (F->phi); Delta~(.)Delta~(.)
 * is the ORDINARY M_N product. p_s = prod_l d_l^{-2}; sum_s p_s = 1. */
void aic_factorize_delta_prime(acb_mat_t out, const aic_factorize *F,
                               const acb_mat_t X, slong prec);

/* Per-block Choi of Delta' (the CP teeth, §C2-style): for block j of B = (+)_l
 * M_{d_l}, C_j = sum_{a,b in [0,d_j)} E_ab (x) Delta'(E^{(j)}_ab) in M_{d_j} (x) M_N
 * (left/major = the d_j source factor, Convention A of aic_ucp.h). Delta' is CP
 * iff every C_j is PSD (block-decomposable map). `Cj` is caller-init'd
 * (d_j*N) x (d_j*N). The block dim d_j = F->v->B->d[j]. */
void aic_factorize_delta_block_choi(acb_mat_t Cj, const aic_factorize *F,
                                    slong j, slong prec);

/* Precompute Delta'(I_B)^{-1/2} (.tex:2799) and store it in F->deltaI_invsqrt;
 * sets F->delta_ready = 1. ASSERTS ||Delta'(I_B) - 1_H||_op < 1 (Rule 4, shard H
 * #7 — the inverse-sqrt basin; aic_funcalc_xpow needs ||A - I||_op < 1) BEFORE
 * the (-1/2)-power. Idempotent: re-clears a stale deltaI_invsqrt. Must run before
 * aic_factorize_delta. */
void aic_factorize_delta_build(aic_factorize *F, slong prec);

/* Delta(X) = Delta'(I_B)^{-1/2} Delta'(X) Delta'(I_B)^{-1/2} (.tex:2799), the UCP
 * encode map. REQUIRES aic_factorize_delta_build first (asserts F->delta_ready).
 * X is n_B x n_B; `out` is N x N (caller-init'd). Unital by construction
 * (Delta(I_B) = 1_H), CP (congruence preserves the CP of Delta'). */
void aic_factorize_delta(acb_mat_t out, const aic_factorize *F,
                         const acb_mat_t X, slong prec);

/* ---- Increment F3 — the UCP DECODE map Upsilon via lem_RC (.tex:2840-2899) -
 * (src/aic_factorize_upsilon.c: W_j, R_j, C_j, xi_j, the lem_RC asserts;
 *  src/aic_factorize_upsilon2.c: L_j, Upsilon'_j, the Upsilon unitalization). */

/* Per-block W_j: H -> L_j (x) E_j, the Choi/Stinespring of the UCP Delta (F2),
 * built L_j-MAJOR DIRECTLY (D2). Route: build the Convention-A Choi of the j-th
 * CP term Delta_j: M_{d_j} -> B(H) of Delta (reuse the F2 block-Choi LAYOUT with
 * the UCP Delta swapped for Delta'), extract Kraus {D_{j,c}: H->L_j} (each
 * d_j x N) via aic_ucp_choi_to_kraus_latd, then stack
 *   W_j[a*e_j + c, p] = D_{j,c}[a, p]   (a in L_j, c in E_j, p in H),  e_j = E_j.
 * `Wj` is OUTPUT: aic_mat-init'd HERE to (d_j*e_j) x N (caller clears it);
 * *e_j_out receives E_j (the per-block Stinespring rank). REQUIRES F->delta_ready.
 * Do NOT use aic_ucp_kraus_to_stinespring (it is ancilla(E_j)-MAJOR, D2 warning).
 * The Choi->Kraus extraction is TOLERANCE-AWARE (FINDINGS §C14): the UCP Delta is
 * CP only to O(eta^2), so the per-block Choi has a small O(eta^2) negative
 * eigenvalue, projected onto the PSD cone (aborts only on a genuine non-CP
 * eigenvalue <= -1e-3). `clipped_out` (if non-NULL) receives the clipped negative
 * cone-defect mass for this block (the genuine-bug guard asserts it is small). */
void aic_factorize_upsilon_Wj(acb_mat_t Wj, slong *e_j_out, double *clipped_out,
                              const aic_factorize *F, slong j, slong prec);

/* R_j = sum_s p_{js} (U_{js}^dag (x) 1_{E_j}) W_j W_j^dag (U_{js} (x) 1_{E_j})
 * (.tex:2843, R_def), with U_{js} the PER-BLOCK d_j x d_j Paulis S_ab=X^a Z^b
 * (d_j^2 of them, weight d_j^{-2}; aic_dhom_pauli — NOT F2's whole-B join, D3).
 * R_j in B(L_j (x) E_j), L_j the LEFT factor. `Rj` is caller-init'd
 * (d_j*e_j) x (d_j*e_j); `Wj` is the (d_j*e_j) x N block matrix; e_j = E_j. */
void aic_factorize_upsilon_Rj(acb_mat_t Rj, const aic_factorize *F, slong j,
                              const acb_mat_t Wj, slong e_j, slong prec);

/* C_j = (1/d_{L_j}) Tr_{L_j} R_j (D1; lem_RC(i) R_j = 1_{L_j} (x) C_j). Traces
 * the LEFT/MAJOR factor L_j (= aic_mat_partial_trace_left, a=d_j, b=e_j), leaving
 * the e_j x e_j operator on E_j. `Cj` caller-init'd e_j x e_j. */
void aic_factorize_upsilon_Cj(acb_mat_t Cj, const acb_mat_t Rj, slong d_j,
                              slong e_j, slong prec);

/* Build the full F3 decode-map data (V, W_j, C_j, xi_j, L_j, Upsilon'(1_H)^{-1/2})
 * and cache it in F; sets F->upsilon_ready = 1. REQUIRES F->delta_ready (asserts).
 * For each block j: W_j (D2), R_j (D3), C_j (D1), then ASSERTS lem_RC(ii)
 * sigma_max(C_j) >= 1 - tol_sigma (Rule 4; abort on a straddling ball — fail
 * loud); xi_j = top RIGHT singular vector of C_j (D4); L_j (D5); finally
 * Upsilon'(1_H) and ASSERTS ||Upsilon'(1_H) - 1_B||_op < 1 (the inverse-sqrt
 * basin) before the (-1/2)-power -> upsI_invsqrt (D6). Idempotent: re-clears a
 * stale F3 cache. tol_sigma is the O(eta) slack on the sigma_max(C_j) lower bound
 * (e.g. 0.3 — generous; the eta=0 oracle gives exactly 1). */
void aic_factorize_upsilon_build(aic_factorize *F, double tol_sigma, slong prec);

/* Upsilon'_j(X) = L_j^dag (1_F (x) Phi(X)) L_j : B(H) -> B(L_j) (.tex:2869, the
 * F-LEFT ordering D5). X is N x N; `out` is d_j x d_j (caller-init'd). REQUIRES
 * F->upsilon_ready. The manifestly-CP per-block decode component. */
void aic_factorize_upsilon_prime_block(acb_mat_t out, const aic_factorize *F,
                                       slong j, const acb_mat_t X, slong prec);

/* Upsilon'(X) = (Upsilon'_1(X), ..., Upsilon'_m(X)) : B(H) -> B, the block-
 * diagonal join (.tex:2867). X is N x N; `out` is n_B x n_B (caller-init'd,
 * written block-diagonal). REQUIRES F->upsilon_ready. */
void aic_factorize_upsilon_prime(acb_mat_t out, const aic_factorize *F,
                                 const acb_mat_t X, slong prec);

/* Upsilon(X) = Upsilon'(1_H)^{-1/2} Upsilon'(X) Upsilon'(1_H)^{-1/2} (.tex:2897),
 * the UCP decode map B(H) -> B. X is N x N; `out` is n_B x n_B (caller-init'd, the
 * block-diagonal B-representative). REQUIRES F->upsilon_ready. Unital by
 * construction (Upsilon(1_H) = I_B), CP (congruence preserves the CP of Upsilon'). */
void aic_factorize_upsilon(acb_mat_t out, const aic_factorize *F,
                           const acb_mat_t X, slong prec);

/* D5 PIN (test/diagnostic): worst over B matrix units of ||Upsilon'_j(Delta(E_i))_j
 * - (E_i)_j||_op with the (x)1_F ordering chosen by `f_left` (1=F-LEFT correct;
 * 0=F-RIGHT wrong) for BOTH L_j and Upsilon'_j (V stays F-major). At eta>0 with
 * r>1, f_left=1 gives O(eta), f_left=0 gives O(1) — the decisive D5 distinguisher.
 * REQUIRES F->upsilon_ready. */
double aic_factorize_upsilon_d5_pin(const aic_factorize *F, int f_left, slong prec);

/* ---- Increment F4.1 — end-to-end verification + dual channels --------------
 * (src/aic_factorize_verify.c: the composed-map Choi + eig-free bound;
 *  src/aic_factorize_dual.c: the dual-channel Kraus read-off.)
 * Realizes the th_factorization headline bounds (.tex:2733 DelUps, .tex:2739
 * ||Upsilon Delta - 1_B||_cb) and the dual factorization (.tex:2159). Design:
 * docs/research/factorize_f4_design.md §A/§B. */

/* J_DelUps = Choi(Delta Upsilon) - Choi(Phi) (.tex:2733), N^2 x N^2, Hermitian
 * indefinite (a difference of two unital UCP Convention-A Choi). Delta Upsilon(W)
 * = Delta(Upsilon(W)) is a self-map on M_N; its Choi is built by looping E_pq over
 * the N^2 AMBIENT M_N standard units (NOT aic_dhom_B_matunit), applying Upsilon
 * then Delta, writing the N x N output into Convention-A position (p,q); Choi(Phi)
 * = aic_ucp_kraus_to_choi(., F->phi). `J` caller-init'd N^2 x N^2. REQUIRES
 * F->delta_ready && F->upsilon_ready (asserts, Rule 4). */
void aic_factorize_choi_delups(acb_mat_t J, const aic_factorize *F, slong prec);

/* J_UpsDel = Choi(Upsilon Delta - 1_B) (.tex:2739), route (i): the AMBIENT
 * M_{n_B} Convention-A Choi, n_B^2 x n_B^2. Loops E_pq over the n_B^2 ambient
 * M_{n_B} standard units (DRIFT #1: MANUAL, NOT aic_dhom_B_matunit), applies Delta
 * then Upsilon, subtracts 1_B(E_pq)=E_pq (1_B is the identity MAP on M_{n_B}).
 * `J` caller-init'd n_B^2 x n_B^2. REQUIRES F->delta_ready && F->upsilon_ready. */
void aic_factorize_choi_upsdel(acb_mat_t J, const aic_factorize *F, slong prec);

/* Eig-free certified bracket on ||Delta Upsilon - Phi||_cb (resp.
 * ||Upsilon Delta - 1_B||_cb): builds the respective J, calls
 * aic_cbnorm_eigfree_ball_choi with n = N (resp. n_B). On return hi = 2||J||_F is
 * a RIGOROUS per-instance UPPER bound (lo ~ ||J||_F/n); eta=0 gives [0,0]. NOT
 * dimension-faithful (the 2N looseness is O(N); the tight SDP + dim-independence
 * canary are F4.2). `lo`, `hi` caller-init'd arb_t. */
void aic_factorize_eigfree_delups(arb_t lo, arb_t hi, const aic_factorize *F,
                                  slong prec);
void aic_factorize_eigfree_upsdel(arb_t lo, arb_t hi, const aic_factorize *F,
                                  slong prec);

/* Dual channels (.tex:2159: Dec = Delta*, Enc = Upsilon*). Build the Convention-A
 * Choi of Delta (resp. Upsilon) by the E_pq pattern, then extract Kraus via
 * aic_ucp_choi_to_kraus_latd. The returned aic_ucp_kraus stores {K_a}; the dual's
 * Schrodinger action is Phi*(rho) = Sum_a K_a rho K_a^dag. DRIFT #2 dims:
 *   Dec (= Delta*): dim_K = n_B, dim_H = N (Choi(Delta) is (n_B*N)^2).
 *   Enc (= Upsilon*): dim_K = N, dim_H = n_B (Choi(Upsilon) is (N*n_B)^2).
 * The caller OWNS the returned struct and must aic_ucp_kraus_clear it. Dec is TP
 * iff Delta unital; Enc is TP iff Upsilon unital. REQUIRES F->delta_ready (Dec) /
 * F->upsilon_ready (Enc). */
void aic_factorize_dec_kraus(aic_ucp_kraus *dec, const aic_factorize *F, slong prec);
void aic_factorize_enc_kraus(aic_ucp_kraus *enc, const aic_factorize *F, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_FACTORIZE_H */
