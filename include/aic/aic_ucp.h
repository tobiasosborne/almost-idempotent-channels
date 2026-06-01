/* aic_ucp.h — UCP-map representations and the carrier (bead aic-c7n, module
 * "ucp"). Realizes the Choi/Stinespring representations (approximate_algebras.tex
 * :1566-1633, prop_KLHG :1635), the cb-norm bound (:1717-1719), and the carrier
 * results lem_carrier (:1724), cor_carrier (:1731), PhiX_M (:1740).
 *
 * DATA MODEL — Kraus operators, HEISENBERG picture (LOAD-BEARING, pinned here).
 *   A UCP map Phi : B(K) -> B(H) is stored as r Kraus operators K_a : H -> K,
 *   each a dim_K x dim_H acb_mat. The action is the OBSERVABLE (Heisenberg) form
 *
 *       Phi(X) = sum_a  K_a^dag X K_a            (X in B(K)),
 *
 *   which is the .tex:1570 isometry form Phi(X) = V^dag (X (x) 1_F) V (with
 *   V^dag V = 1_H, .tex:1571) expanded in an F-basis, V the column-stack
 *   V[a*dim_K + i, j] = K_a[i, j]  (V : H -> K (x) F, F = C^r). Unitality is
 *
 *       Phi(1_K) = sum_a K_a^dag K_a = 1_H.                          [.tex:1571]
 *
 *   CAREFUL: this is the OBSERVABLE map, NOT the Schrodinger/state dual
 *   T = Phi^* (X) = sum_a K_a X K_a^dag. Getting the dual backwards silently
 *   transposes everything (CLAUDE.md "UCP vs CPTP" callout). The unital test
 *   sum_a K_a^dag K_a = 1_H pins the direction and is asserted by the tests.
 *
 * CHOI MATRIX — Convention A (paper / Watrous / QETLAB / Qiskit). With ONB {e_i}
 * of K and E_{ij} = |e_i><e_j|,
 *
 *       C_Phi = sum_{i,j} E_{ij} (x) Phi(E_{ij})   in B(K) (x) B(H),  [shard G; :1575]
 *
 *   with the K factor LEFT/MAJOR and the H factor RIGHT, matching aic_mat.h's
 *   left-major (row-major) Kronecker. So C_Phi is (dim_K*dim_H) x (dim_K*dim_H).
 *   Hand-deriving Phi(E_{ij})[a,b] from Phi(X) = sum_x K_x^dag X K_x: since
 *   (K_x^dag E_{ij} K_x)[a,b] = sum_{p,q} conj(K_x[p,a]) E_{ij}[p,q] K_x[q,b] and
 *   E_{ij}[p,q] = delta_{pi} delta_{qj}, we get
 *       Phi(E_{ij})[a,b] = sum_x conj(K_x[i,a]) * K_x[j,b],
 *   so (conjugation on the FIRST, (i,a), factor)
 *
 *       C_Phi[i*dim_H + a, j*dim_H + b] = sum_x conj(K_x[i,a]) * K_x[j,b].
 *
 *   Example (dim_K=dim_H=2): row index i*2+a, col index j*2+b, i,j in {0,1}
 *   the K (left) indices and a,b in {0,1} the H (right) indices.
 *   Phi is CP iff C_Phi >= 0 (PSD); Phi is unital iff Tr_K(C_Phi) = 1_H, i.e.
 *   tracing out the LEFT (K) factor leaves the dim_H identity — this uses
 *   aic_mat_partial_trace_left (traces out the LEFT factor, keeps the RIGHT).
 *
 * Number paths. The Kraus<->Choi conversions, unital/CP/cbnorm checks, carrier
 * Q, Stinespring and the cor_carrier/PhiX_M checks live on the CERTIFIED arb
 * path (acb_mat, explicit prec). The Choi->Kraus EXTRACTION (PSD eigendecomp)
 * and the carrier RANK/range need a degenerate-spectrum eigensolver, so they
 * live on the LAPACK DOUBLE path (aic_ucp_*_latd) — the arb certified eig path
 * needs a SIMPLE spectrum and the Choi spectrum is generically degenerate (bead
 * aic-w4o.1). The certified-rank carrier is therefore a documented module TODO.
 */
#ifndef AIC_UCP_H
#define AIC_UCP_H

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/* A UCP map Phi : B(K) -> B(H) as r Kraus operators K_a : H -> K (each
 * dim_K x dim_H). `K[a]` is an initialised acb_mat_t. See file docstring for
 * the Heisenberg convention Phi(X) = sum_a K_a^dag X K_a. */
typedef struct {
    slong dim_K;   /* domain B(K) acts on C^{dim_K} */
    slong dim_H;   /* codomain B(H) acts on C^{dim_H} */
    slong r;       /* number of Kraus operators (Stinespring ancilla dim) */
    acb_mat_t *K;  /* r matrices, each dim_K x dim_H */
} aic_ucp_kraus;

/* Allocate a Kraus container with r operators each dim_K x dim_H, all zeroed.
 * Fields are set; every K[a] is acb_mat_init'd. Free with aic_ucp_kraus_clear. */
void aic_ucp_kraus_init(aic_ucp_kraus *phi, slong dim_K, slong dim_H, slong r);

/* Release all storage owned by `phi` (every K[a] and the array). */
void aic_ucp_kraus_clear(aic_ucp_kraus *phi);

/* Apply Phi: out = sum_a K_a^dag X K_a (the Heisenberg/observable action,
 * .tex:1570). X is dim_K x dim_K; out is dim_H x dim_H (caller-initialised). */
void aic_ucp_apply(acb_mat_t out, const aic_ucp_kraus *phi, const acb_mat_t X,
                   slong prec);

/* Build the Choi matrix C_Phi (Convention A, file docstring) from the Kraus
 * operators. `C` must be initialised (dim_K*dim_H) x (dim_K*dim_H). */
void aic_ucp_kraus_to_choi(acb_mat_t C, const aic_ucp_kraus *phi, slong prec);

/* Stinespring isometry V : H -> K (x) F, the column-stack
 * V[a*dim_K + i, j] = K_a[i, j]. `V` must be initialised
 * (dim_K*r) x dim_H. Unitality is V^dag V = 1_H. */
void aic_ucp_kraus_to_stinespring(acb_mat_t V, const aic_ucp_kraus *phi,
                                  slong prec);

/* Unital defect from the Choi matrix: ||Tr_K(C_Phi) - 1_H||_op as a certified
 * arb ball (0 iff Phi is unital). Traces out the LEFT (K) factor. `out` must be
 * an initialised arb_t. */
void aic_ucp_unital_defect_choi(arb_t out, const acb_mat_t C, slong dim_K,
                                slong dim_H, slong prec);

/* Unital defect straight from Kraus: ||sum_a K_a^dag K_a - 1_H||_op (certified
 * arb ball). Must agree with aic_ucp_unital_defect_choi. NOTE: that cross-check
 * is transpose-invariant, so it validates only the partial-trace DIRECTION
 * (tracing out the LEFT/K factor), NOT the per-entry conjugation convention. The
 * conjugation/index convention is validated by the Choi-convention oracle test
 * (test_ucp.c: Choi block (i,j) == aic_ucp_apply(E_{ij}), with a complex
 * asymmetric channel). `out` must be initialised. */
void aic_ucp_unital_defect_kraus(arb_t out, const aic_ucp_kraus *phi,
                                 slong prec);

/* Complete-positivity certificate via the Choi matrix (.tex: C_Phi >= 0 iff CP).
 * Computes the largest eigenvalue of (-C_Phi) as a certified ball using the
 * degeneracy-robust aic_mat_herm_max_eig (acb_mat_eig_global_enclosure):
 *   returns +1 if that ball is CERTAINLY <= tol  (C_Phi certified PSD),
 *   returns  0 if it is CERTAINLY  > tol         (C_Phi certified NOT PSD),
 *   ABORTS (fail loud, Rule 4) if the ball straddles tol — the spectral gap is
 *   unresolved at this prec and a verdict would be a guess.
 * `tol` is a small nonneg arb slack for the zero boundary (pass 0 for exact). */
int aic_ucp_is_cp_choi(const acb_mat_t C, const arb_t tol, slong prec);

/* cb-norm of a CP map: ||Phi||_cb = ||Phi(1_K)||_op (Paulsen Prop. 3.6;
 * .tex:1717-1719 gives <=1, with equality =1 for UCP per .tex:359). Returns the
 * certified ball ||Phi(1_K)||_op. THIS IS THE ONLY cb-norm in scope; the
 * idempotency defect ||Phi^2-Phi||_cb needs the diamond-norm SDP (bead aic-d24,
 * out of scope). `out` must be initialised. */
void aic_ucp_cbnorm_cp(arb_t out, const aic_ucp_kraus *phi, slong prec);

/* Carrier operator Q = sum_a K_a K_a^dag on K (dim_K x dim_K, Hermitian PSD).
 * lem_carrier (.tex:1724): the carrier M is the support (range) of Q — it is the
 * K-marginal of Im(V), Q = Tr_F(V V^dag). `Q` must be initialised dim_K x dim_K.
 * The RANGE/RANK of Q (the carrier itself) needs a degenerate eigensolver: use
 * the double path aic_ucp_carrier_rank_latd. */
void aic_ucp_carrier_Q(acb_mat_t Q, const aic_ucp_kraus *phi, slong prec);

/* cor_carrier check (.tex:1731-1732): X in B(K) annihilates the carrier M iff
 * (X (x) 1_F) V = 0. Returns ||(X (x) 1_F) V||_op = ||stack_a(X K_a)||_op as a
 * certified ball (0 iff Ker X >= M). X is dim_K x dim_K. `out` initialised. */
void aic_ucp_carrier_annihilate_defect(arb_t out, const aic_ucp_kraus *phi,
                                       const acb_mat_t X, slong prec);

/* PhiX_M defect (.tex:1740): ||Phi(X) - Phi(Pi_M X Pi_M)||_op as a certified
 * ball, where Pi_M is the orthogonal projection onto the carrier M (passed in as
 * the dim_K x dim_K projector `Pi_M`). The identity Phi(X) = Phi(Pi_M X Pi_M)
 * holds for ANY UCP map (.tex:1742), so the defect is 0 up to rounding. X and
 * Pi_M are dim_K x dim_K. `out` initialised. */
void aic_ucp_phiX_M_defect(arb_t out, const aic_ucp_kraus *phi,
                           const acb_mat_t X, const acb_mat_t Pi_M, slong prec);

/* --- composition + Choi-difference (bead aic-d24; aic_ucp_compose.c) --- */

/* Heisenberg COMPOSITION (Phi o Psi)(X) = Phi(Psi(X)). For UCP maps in the
 * Kraus form Phi(Y)=sum_a K_a^dag Y K_a (K_a: dim_K(phi) x dim_H(phi)) and
 * Psi(X)=sum_b L_b^dag X L_b (L_b: dim_K(psi) x dim_H(psi)),
 *   (Phi o Psi)(X) = sum_a K_a^dag (sum_b L_b^dag X L_b) K_a
 *                  = sum_{a,b} (L_b K_a)^dag X (L_b K_a),
 * so the Kraus operators of Phi o Psi are { L_b K_a } = { psi.K[b] @ phi.K[a] }.
 * TYPE CHECK: Psi outputs in B(K_psi) which Phi must accept, so this REQUIRES
 *   phi->dim_K == psi->dim_H            (asserted, fail loud).
 * The product L_b K_a is (dim_K(psi) x dim_H(psi)) * (dim_K(phi) x dim_H(phi))
 * with dim_K(phi)==dim_H(psi), so it is dim_K(psi) x dim_H(phi). RESULT shape:
 *   out->dim_K = psi->dim_K,  out->dim_H = phi->dim_H,  out->r = phi->r * psi->r.
 * For Phi^2 = compose(Phi,Phi) (requires dim_K==dim_H): { K_b K_a }, n x n,
 * r^2 operators. `out` is aic_ucp_kraus_init'd HERE (caller clears it). */
void aic_ucp_compose(aic_ucp_kraus *out, const aic_ucp_kraus *phi,
                     const aic_ucp_kraus *psi, slong prec);

/* Choi DIFFERENCE C = Choi(phi1) - Choi(phi2), Convention A (header). Both maps
 * must share (dim_K, dim_H); `C` must be initialised (dim_K*dim_H) square. The
 * eta-idempotence defect Lambda = Phi^2 - Phi has Choi C_{Phi^2} - C_Phi by
 * linearity (.tex:347-354); pass phi1 = compose(Phi,Phi), phi2 = Phi. The result
 * is Hermitian (a difference of Hermitian Choi matrices) but in general NOT PSD
 * (Phi^2-Phi is not CP), which is exactly why the diamond-norm SDP is needed. */
void aic_ucp_choi_diff(acb_mat_t C, const aic_ucp_kraus *phi1,
                       const aic_ucp_kraus *phi2, slong prec);

/* --- channel power + Choi-rank compression (bead aic-xxk; aic_ucp_power.c) --- */

/* Choi-rank COMPRESSION: minimal-Kraus-rank form of `phi`. Builds the Choi
 * matrix C_Phi (Convention A) and re-extracts its minimal PSD eigenbasis via
 * aic_ucp_choi_to_kraus_latd, which keeps only eigenpairs with lambda > the
 * QETLAB threshold (dim_K*dim_H)*eps_mach*||C||_F. The CHANNEL is unchanged
 * (same Choi matrix); only the Kraus count drops to the Choi rank, so
 *   out->r = (count of nonzero Choi eigenvalues) <= dim_K * dim_H.
 * EXACT up to the latd extraction tolerance (the QETLAB threshold ~ 1e-15||C||_F:
 * numerically-zero eigenpairs are dropped — that drop IS the compression).
 * Kraus reps are unique only up to a unitary gauge, so the recovered operators
 * are NOT the originals: compare AS CHANNELS (rebuild Choi), never op-by-op.
 * `out` is aic_ucp_kraus_init'd HERE (by the latd extractor); caller clears it. */
void aic_ucp_compress(aic_ucp_kraus *out, const aic_ucp_kraus *phi, slong prec);

/* Channel POWER Phi^k (Heisenberg map composition Phi o ... o Phi), k >= 0.
 *   k = 0  -> identity channel on C^{dim_H} (r = 1, K_0 = 1_{dim_H});
 *   k = 1  -> a copy of phi (uncompressed; phi's own Kraus set);
 *   k >= 2 -> repeated compose, COMPRESSING AFTER EACH compose so the running
 *             Kraus count stays at the Choi rank (<= dim_H^2) instead of blowing
 *             up as phi->r^k (compose multiplies operator counts). The compress-
 *             each-step is why this is paired with aic_ucp_compress.
 * REQUIRES dim_K == dim_H for k != 1 (the identity and composition are
 * endomorphisms); asserted fail-loud (Rule 4). For a power both compose factors
 * are powers of the same Phi (which commute), so the compose direction is moot.
 * `out` is aic_ucp_kraus_init'd HERE; caller clears it. */
void aic_ucp_power(aic_ucp_kraus *out, const aic_ucp_kraus *phi, slong k,
                   slong prec);

/* --- double path (LAPACK): degenerate-spectrum extraction (aic_ucp_*_latd) --- */

/* Choi->Kraus extraction (.tex / shard G): PSD eigendecomposition of C_Phi,
 * keep eigenpairs with lambda > threshold = (dim_K*dim_H) * eps_mach * ||C||_F
 * (QETLAB-style), K_a = sqrt(lambda_a) * reshape(evec_a, dim_K x dim_H).
 * Uses LAPACKE_zheev (degeneracy-robust; the arb simple-eig path cannot, bead
 * aic-w4o.1). `phi` is OUTPUT: it is aic_ucp_kraus_init'd here to the recovered
 * rank r (<= dim_K*dim_H) and filled (caller must aic_ucp_kraus_clear it). The
 * extracted set is a VALID Kraus rep of the same channel, but only up to the
 * unitary gauge freedom of Kraus decompositions — compare AS CHANNELS (rebuild
 * Choi), never operator-by-operator. C is the (dim_K*dim_H)^2 Choi matrix. */
void aic_ucp_choi_to_kraus_latd(aic_ucp_kraus *phi, const acb_mat_t C,
                                slong dim_K, slong dim_H);

/* TOLERANCE-AWARE Choi->Kraus with a PSD-cone projection (bead aic-tff, FINDINGS
 * §C14). Same eigendecomposition as the strict variant, with three regimes:
 *   lambda > keep_thr=(dim_K*dim_H)*eps_mach*||C||_F  -> KEEP as Kraus;
 *   lambda in (-neg_tol, keep_thr]                    -> DROP (cone-defect/noise);
 *                                                        |lambda| (if <0) added to
 *                                                        *clipped_neg_out;
 *   lambda <= -neg_tol                                -> FAIL LOUD (genuine non-CP).
 * This extracts the Kraus of the NEAREST genuinely-CP map (C projected onto the PSD
 * cone), needed for the almost-idempotent (O(eta)-CP) Delta of th_factorization
 * whose per-block Choi is PSD only to O(eta^2). `clipped_neg_out` (if non-NULL)
 * receives the TOTAL clipped negative mass = sum of |lambda| over the dropped
 * negative eigenvalues (the certified cone-defect magnitude; the genuine-bug guard
 * checks it is O(eta^2)-small). `neg_tol` must be >= 0. The strict
 * aic_ucp_choi_to_kraus_latd delegates here with neg_tol = keep_thr. `phi` is
 * OUTPUT (init'd here; caller clears). */
void aic_ucp_choi_to_kraus_latd_tol(aic_ucp_kraus *phi, const acb_mat_t C,
                                    slong dim_K, slong dim_H, double neg_tol,
                                    double *clipped_neg_out);

/* Carrier rank: number of eigenvalues of Q = sum_a K_a K_a^dag above
 * threshold = dim_K * eps_mach * ||Q||_F (double path, LAPACKE_zheev). This is
 * the dimension of the carrier M. The CERTIFIED rank (arb) is gap-dependent and
 * blocked on aic-w4o.1 — this is the uncertified fast route. `Q` is the
 * dim_K x dim_K carrier operator from aic_ucp_carrier_Q. */
slong aic_ucp_carrier_rank_latd(const acb_mat_t Q, slong dim_K);

/* Certified carrier rank: rank of the Hermitian PSD carrier Q = sum_a K_a K_a^dag
 * via aic_mat_certified_rank with thr = dim_K * eps_mach * ||Q||_F (an arb ball,
 * matching the double-path aic_ucp_carrier_rank_latd threshold; eps_mach =
 * DBL_EPSILON = 2^-52). The certified counterpart that retires the aic-w4o.1 /
 * FINDINGS §D5 deferral for the carrier dimension. ABORTS (Rule 4) if the rank
 * gap is unresolved at `prec` (a carrier eigenvalue ball straddles thr). */
slong aic_ucp_carrier_rank(const acb_mat_t Q, slong dim_K, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_UCP_H */
