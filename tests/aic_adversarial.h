/* aic_adversarial.h — reusable ADVERSARIAL INSTANCE GENERATORS (bead aic-dbo.2,
 * Increment 1). An "evil-matrix corpus": deterministic generators over acb_mat_t,
 * each producing a matrix that provably exhibits a named adversarial property as a
 * function of a KNOB. Test/bench infrastructure (lives under tests/, not src/),
 * feeding the shared benchmark corpus aic-f9u.1.
 *
 * WHY THIS EXISTS (CLAUDE.md adversarial-first directive, docs/adversarial/README.md).
 * The happy path is the least useful test. Every AIC routine must be punished by an
 * instance engineered to break it; this module is the curated catalogue of those
 * instances, lifted out of ad-hoc per-test builders so they are reusable AND so each
 * generator's adversarial property is asserted ONCE, centrally, by a self-test (the
 * corpus's own teeth, tests/test_adversarial.c). A generator whose output does NOT
 * vary with its knob, or does NOT achieve the claimed property, is the failure mode
 * to avoid — the project's recurring "test that can't fail", here at the corpus
 * level. The self-tests assert knob -> property at >=2 knob values (mild vs lethal)
 * so the corpus cannot silently go toothless.
 *
 * DETERMINISM. Every generator is a fixed closed-form construction (no RNG). Repeated
 * calls with the same knob produce bit-identical matrices, so a test/bench is
 * reproducible. Where a "random" off-diagonal is wanted (gen2), a fixed deterministic
 * pattern (a small rational LCG, fixed seed AIC_ADV_SEED) stands in.
 *
 * CATALOGUE PROVENANCE (Law 1). Each generator cites its docs/adversarial/nla.md or
 * domain.md family/sub-section and, where the catalogue cites the paper, the
 * approximate_algebras.tex line. Increment 1 (aic_adversarial.c / aic_adversarial_nla.c)
 * has seven NLA generators; Increment 2 (aic_adversarial_domain.c) adds the first two
 * channel/UCP-map generators. The full list:
 *
 * NLA generators (bare acb_mat_t; Increment 1):
 *   1 jordan_t13        nla 1a, tex:540-544  : non-Lipschitz t^{1/3} spectral burst
 *   2 nonnormal_tri     nla 2c               : departure-from-normality / op-vs-rho gap
 *   3 degenerate_proj   nla 4b (lethal #1)   : exact {0,1} projector, eig-free sgn only
 *   4 near_degen_herm   nla 4a               : tunable Hermitian gap -> simple-eig kill
 *   5 graded_diag       nla 5c               : dynamic range = condition number kappa
 *   6 boundary_x2I      nla 7a (lethal #5)   : ||X^2-I||=1 straddle, funcalc domain edge
 *   7 propP_delta       nla 7b               : ||P^2-P||=delta, prop_P basin edge 1/4
 *
 * Channel/UCP-map generators (aic_ucp_kraus; Increment 2, tranche 1):
 *   chan_cb_op_gap      domain 1B, tex:366-388 : cb-norm vs operator-norm gap (measure-prepare)
 *   chan_depol_boundary domain 2A, tex:516-525 : eta->1/4 theta(2Phi-1) basin edge
 *   chan_unital_defect  domain 1D, tex:432/672 : unital-but-barely Phi(I)=I+delta_u E
 *   chan_carrier_dropout domain 1C, tex:1724/1731 : near-degenerate carrier, rank near-drop
 *   chan_carrier_dropout_dense 1C-dense, tex:1724/1731 : DENSIFIED non-normal carrier (straddle + convention-sensitive)
 *   chan_blockalg       domain 3D, tex:484/1249 : eps~c/n dimension-blowup block algebra (+)_j M_d
 *   chan_noncomm_boundary README:57, tex:347/378 : NON-COMMUTATIVE eta_cb=1/4-kappa boundary (id(x)cb_op_gap)
 *   chan_conc_defect    domain 2B, tex:2192-2204/2385-2422 : RANK-1 defect on a near-degenerate subspace (gap_sub)
 *
 * API SHAPE (reusable by tests AND benches). Each NLA generator takes a caller-managed
 * acb_mat_t `out` that the caller has init'd to the right shape (asserted), the
 * knob(s), and a prec; outputs are written into `out`; the caller clears it. The
 * channel generators (aic_adv_chan_*) instead take an aic_ucp_kraus `out` that the
 * generator INITIALISES itself (matching the aic_channels.h constructor convention);
 * the caller clears it with aic_ucp_kraus_clear. Names are aic_adv_<family>. No
 * global state.
 */
#ifndef AIC_ADVERSARIAL_H
#define AIC_ADVERSARIAL_H

#include <flint/acb_mat.h>

#include "aic/aic_ucp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Fixed seed for the one generator (gen2) that needs a reproducible "random"
 * off-diagonal pattern. Documented so the instance is regenerable bit-for-bit. */
#define AIC_ADV_SEED 0x9e3779b97f4a7c15UL

/* gen1 — Jordan t^{1/3} (nla 1a, tex:540-544). The paper's 3x3 companion-style
 * matrix X(t) = [[0,1,0],[0,0,1],[t,0,0]] with spec(X) = {t^{1/3} omega^k}. KNOB t.
 * `out` must be 3x3 (asserted). Adversarial property: a size-t entry shifts the
 * spectral radius by t^{1/3} (non-Lipschitz); the self-test checks
 * rho(X) = |det(X)|^{1/3} = t^{1/3} (char poly is lambda^3 - t). */
void aic_adv_jordan_t13(acb_mat_t out, double t, slong prec);

/* gen2 — Non-normal upper-triangular (nla 2c, DETERMINISTIC variant). Diagonal
 * {(i+1)/n} placed so eigenvalues are controlled; strictly-upper entries are
 * c * fixed_rational[i][j] (fixed LCG pattern, no RNG). KNOB c = departure from
 * normality. `out` must be n x n (asserted, n>=2). Adversarial property:
 * ||[X,X^dag]||_F grows monotonically with c, AND the op-vs-spectral-radius gap
 * ||I-X^2||_op - rho(I-X^2) opens with c (the regime aic-8hz's global Newton sgn
 * handles and the op-basin Newton-Schulz rejects). rho(I-X^2) = max_i|1-((i+1)/n)^2|
 * via the placed eigenvalues. */
void aic_adv_nonnormal_tri(acb_mat_t out, slong n, double c, slong prec);

/* gen3 — Exact-degeneracy projector (nla 4b, lethal #1). Orthogonal projector P of
 * rank `rank` on C^dim, P = sum_{k<rank} |q_k><q_k| for an explicit orthonormal set
 * (columns of a normalised real Hadamard-like / DFT unitary so P is NON-diagonal).
 * KNOBS dim, rank. `out` must be dim x dim (asserted; 0 < rank < dim). Adversarial
 * property: ||P^2-P|| = 0 to machine floor (exact idempotent), tr(P) = rank, and
 * X = 2P-I is an exact sign matrix ((2P-I)^2 = I). This is the instance FLINT's
 * simple-spectrum eig cannot isolate (spectrum {0,1} repeated) — the eig-free sgn
 * path's reason to exist. */
void aic_adv_degenerate_proj(acb_mat_t out, slong dim, slong rank, slong prec);

/* gen4 — Near-degenerate Hermitian (nla 4a). Diagonal Hermitian
 * H = diag(1 - gap/2, 1 + gap/2, 3, 4, ..., n) with the two smallest-separated
 * eigenvalues exactly `gap` apart. KNOB gap. `out` must be n x n (asserted, n>=2).
 * Adversarial property: the minimal eigenvalue gap is exactly `gap` (within tol),
 * attacking aic_mat_eig_hermitian's simple-spectrum requirement as gap -> 0. */
void aic_adv_near_degen_herm(acb_mat_t out, slong n, double gap, slong prec);

/* gen5 — Graded / ill-conditioned diagonal (nla 5c). D = diag(1, r, r^2, ...,
 * r^{n-1}) with r = range^{1/(n-1)} so the dynamic range (and 2-norm condition
 * number kappa = sigma_max/sigma_min) is exactly `range`. KNOB range (>= 1).
 * `out` must be n x n (asserted, n>=2). Adversarial property: kappa(D) = range
 * within tol (precision-consumption / wide-dynamic-range attack). */
void aic_adv_graded_diag(acb_mat_t out, slong n, double range, slong prec);

/* gen6 — Near-boundary ||X^2-I||=1 straddle (nla 7a, lethal #5). X(s) =
 * diag(sqrt(2-s), 1, ..., 1) so ||X^2-I||_op = |(2-s)-1| = |1-s| = 1-s for
 * s in (0,2). KNOB s (signed): s>0 => certified ||X^2-I||<1 (just inside the funcalc
 * domain); s<0 => ||X^2-I|| > 1 (just outside). `out` must be n x n (asserted, n>=1).
 * Adversarial property: certified ||X^2-I||_op is < 1 for s>0 and >= 1 for s<0 — the
 * exact edge aic_funcalc_int_assert_domain guards; at s=0 the ball straddles 1 and
 * the guard must fail loud (Rule 4). */
void aic_adv_boundary_x2I(acb_mat_t out, slong n, double s, slong prec);

/* gen7 — prop_P delta near 1/4 (nla 7b). P = diag(p, 0, ..., 0) with p chosen so
 * ||P^2-P||_op = |p^2 - p| = delta exactly: p = (1 - sqrt(1 - 4 delta))/2 in [0,1/2]
 * (so 0 <= delta <= 1/4). KNOB delta. `out` must be n x n (asserted, n>=1,
 * 0 <= delta <= 1/4). Adversarial property: certified ||P^2-P||_op = delta within
 * tol — the prop_P / assoc_ecsa basin edge delta < 1/4 (tex:525). */
void aic_adv_propP_delta(acb_mat_t out, slong n, double delta, slong prec);

/* ---- Domain / channel generators (Increment 2) ---------------------------
 * The FIRST generators producing a UCP map (aic_ucp_kraus), not a bare matrix.
 * They punish the channel routines (ucp, cbnorm, assoc_ecsa, opspace). See
 * aic_adversarial_domain.c for the docstring and Kraus derivations. */

/* fam1B — cb-norm vs operator-norm gap (domain.md:75-123; the paper's own
 * motivating measure-prepare example, tex:366-388). UCP self-map on B(C^d),
 * d>=2, in the OBSERVABLE convention Phi(X)=sum_a K_a^dag X K_a:
 *   v = (sqrt(1-eta), sqrt(eta), 0,...,0)  (unit vector, gamma_0 = |v><v|),
 *   K_0 = |v><0|  (d x d, column 0 = v),   K_j = |j><j|  (j = 1 .. d-1),
 * so Phi(X) = Tr(gamma_0 X) P_0 + sum_{j>=1} Tr(P_j X) P_j (the tex:369 map for
 * d=2; exactly-idempotent dephasing pads for d>2). It is unital (<v|v>=1),
 * entanglement-breaking (rank-1 Kraus), Kraus rank d.
 * KNOB eta in [0,1] (asserted; gap maximized near eta~1/4, sweep {0.05..0.20}).
 * Adversarial property (the NAMED cb!=op trap, CLAUDE.md callout): the certified
 * idempotence defect ||Phi^2-Phi||_cb = eta*sqrt(1-eta) EXACTLY (tex:378), which
 * is STRICTLY larger than the single-copy operator-norm defect — a routine that
 * uses op-norm as eta underreports by up to the ampliation factor (tex:484).
 * eta=0 REDUCTION: v=|0>, Phi = complete dephasing, EXACTLY idempotent (defect 0).
 * `out` is aic_ucp_kraus_init'd HERE (caller aic_ucp_kraus_clears it). */
void aic_adv_chan_cb_op_gap(aic_ucp_kraus *out, slong d, double eta, slong prec);

/* fam2A — depolarizing eta->1/4 regularization boundary (domain.md:196-245; the
 * theta(2Phi-1) basin edge rho(Phi^2-Phi) < 1/4, tex:516-525). The standard
 * unital depolarizing self-map on B(C^d), d>=2:
 *   Phi_p(X) = (1-p) X + p (tr X / d) 1_d,    p in [0,1]   (the named 2A instance;
 * delegates to aic_channel_depolarizing). Its defect is EXACTLY
 *   Phi_p^2 - Phi_p = p(1-p) C,   C(X) = (tr X / d) 1_d - X,
 * with C having superoperator eigenvalues {0 on 1_d, -1 on the d^2-1 traceless
 * directions}. KNOB p in [0,1] (asserted 0<=p<=1, d>=2; Rule 4, fail loud).
 * Adversarial property (the rho->1/4 basin edge): the certified defect scales
 * EXACTLY LINEARLY in p(1-p) — rho(Phi_p^2-Phi_p) = p(1-p) (d-INDEPENDENT,
 * maximized = 1/4 at p=1/2, the theta(2Phi-1) basin boundary), and the eig-free
 * bracket lo = (||Choi(C)||_F/n) p(1-p) inherits the exact linearity. The 2A
 * sweep targets eta near 1/4 (e.g. d=2, p=0.789 gives rho~=0.167; p=1/2 gives the
 * exact edge rho=1/4). eta=0 REDUCTION: p=0 (identity) and p=1 (trace-replace
 * conditional expectation, Phi_1^2=Phi_1) are EXACTLY idempotent (defect 0).
 * `out` is aic_ucp_kraus_init'd HERE (caller aic_ucp_kraus_clears it). */
void aic_adv_chan_depol_boundary(aic_ucp_kraus *out, slong d, double p,
                                 slong prec);

/* fam1D — unital-but-barely (domain.md:163-190; the eps-unit axiom ax_eps_unit
 * ||XI-X|| <= eps||X||, tex:432, and the prop_unit exact-unit fix, tex:672). A CP
 * self-map on B(C^d), d>=2, unital ONLY up to delta_u, built from a SINGLE
 * Hermitian Kraus operator (CP trivially) in the OBSERVABLE convention
 * Phi(X)=sum_a K_a^dag X K_a = K_0 X K_0:
 *   K_0 = diag(sqrt(1+delta_u), sqrt(1-delta_u), 1, ..., 1)  on C^d,
 * so Phi(I) = K_0^2 = I + delta_u*E with E = diag(1,-1,0,...,0) (traceless,
 * ||E||_op=1). KNOB delta_u in [0,1) (asserted; real sqrt; Rule 4 fail loud).
 * Adversarial property (the eps-unit axiom edge): the certified UNITAL DEFECT
 * ||Phi(I)-I||_op = ||sum_a K_a^dag K_a - 1_H||_op = delta_u EXACTLY — the
 * O(delta_u) approximate-unit term ax_eps_unit/prop_unit must absorb or refuse
 * (assoc_ecsa assumes Phi_tilde(1)=1 exactly, tex:2181). delta_u=0 REDUCTION:
 * K_0=1_d, Phi=identity, EXACTLY unital, defect 0. `out` is aic_ucp_kraus_init'd
 * HERE (dim_K=dim_H=d, r=1; caller aic_ucp_kraus_clears it). */
void aic_adv_chan_unital_defect(aic_ucp_kraus *out, slong d, double delta_u,
                                slong prec);

/* fam1C — near-degenerate carrier (domain.md:127-159; lem_carrier .tex:1724,
 * cor_carrier .tex:1731). A self-map on B(C^d), d>=2, with a SINGLE Hermitian
 * diagonal Kraus operator
 *   K_0 = diag(1, ..., 1, sqrt(gap))   on C^d,   gap in (0, 1],
 * so the carrier operator that aic_ucp_carrier_Q / aic_ucp_carrier_rank inspect,
 *   Q = sum_a K_a K_a^dag = K_0 K_0^dag = diag(1, ..., 1, gap),
 * has spectrum {1 (x)(d-1), gap}: the last OUTPUT dimension is barely populated,
 * so the carrier rank is d for gap above the routine's threshold
 * thr = dim_K * 2^-52 * ||Q||_F and NEARLY drops to d-1 as gap -> 0. The smallest
 * carrier eigenvalue is EXACTLY gap (the named adversarial quantity). NOT unital
 * for gap != 1 (deliberate: a unital padding would refill the d-th output
 * dimension and kill the rank near-drop; an adversarial instance need not be
 * unital). KNOB gap in (0, 1] (asserted; d >= 2; Rule 4, fail loud).
 * Adversarial property (the near-degenerate carrier, the silent-wrong-rank trap):
 * aic_ucp_carrier_rank must certify rank d while gap > thr and certify the drop
 * to d-1 (or fail loud) as gap -> 0 below thr, never silently return the wrong
 * rank. MEASURED near-drop behavior (Rule 2, discovered not assumed): because Q
 * is DIAGONAL with an exactly-representable entry, the gap eigenvalue is a POINT
 * ball, so the routine ALWAYS CERTIFIES (never straddles) — rank d for gap > thr,
 * rank d-1 at gap=0 (the exact-drop oracle); the densified-carrier STRADDLE/fail-
 * loud path is exercised separately by test_eigvec.c:S6b (FINDINGS §D7). gap=1
 * REDUCTION: K_0 = 1_d, Q = 1_d, EXACT full-rank d carrier (oracle). `out` is
 * aic_ucp_kraus_init'd HERE (dim_K=dim_H=d, r=1; caller aic_ucp_kraus_clears it). */
void aic_adv_chan_carrier_dropout(aic_ucp_kraus *out, slong d, double gap,
                                  slong prec);

/* fam1C-dense — DENSIFIED near-degenerate carrier (domain.md 1C; lem_carrier
 * .tex:1724, cor_carrier .tex:1731; FINDINGS §D7 the prec floor). The hostile-
 * review follow-up to the diagonal 1C (aic_adv_chan_carrier_dropout), closing the
 * two gaps that family leaves: it exercises the STRADDLE/fail-loud contract of
 * aic_ucp_carrier_rank (the diagonal carrier's point ball only ever DECIDES) and
 * it is CONVENTION-SENSITIVE (the diagonal Hermitian 1C cannot catch a
 * sum K K^dag vs sum K^dag K bug). A self-map on B(C^d), d>=2, with a SINGLE
 * NON-NORMAL Kraus operator
 *   K_0 = U diag(1, ..., 1, sqrt(gap)) V^dag   on C^d,   gap in (0, 1],
 * U(cos=3/5,sin=4/5) != V(cos=5/13,sin=12/13) being two DIFFERENT fixed rational-
 * Givens chain unitaries (no RNG). The carrier
 *   Q = sum_a K_a K_a^dag = K_0 K_0^dag = U diag(1, ..., 1, gap) U^dag
 * is NON-DIAGONAL, spectrum {1 (x)(d-1), gap} (same as diagonal 1C but dense
 * eigenvectors). KNOB gap in (0,1]; asserts gap>0, d>=2 (Rule 4, fail loud). NOT
 * unital for gap!=1 (deliberate carrier-only fixture; a unital padding would kill
 * the rank near-drop, same rationale as the diagonal 1C).
 * Adversarial properties (the two gaps closed):
 *   (STRADDLE) MEASURED prec=53, d=3: the densified small-cluster ball radius
 *     ~1.5e-15 STRADDLES thr = dim_K*2^-52*||Q||_F ~9.4e-16 for gap below ~1e-15,
 *     so aic_ucp_carrier_rank fail-loud-ABORTS ("STRADDLES"); at prec>=64 the ball
 *     radius drops to ~7e-19 and it certifies (rank d for gap>>thr, d-1 for the
 *     near-zero gap). The diagonal 1C's point ball can NEVER reach this (FINDINGS
 *     §D7, aic_adversarial_carrier.c:50-64).
 *   (CONVENTION) certified ||sum K K^dag - sum K^dag K||_op > 0 (MEASURED 0.223 at
 *     d=3 gap=0.5, up to 0.605 at d=4 gap->0), because sum K^dag K =
 *     V diag(1,...,1,gap) V^dag differs from the carrier sum K K^dag =
 *     U diag(1,...,1,gap) U^dag (U!=V). The diagonal Hermitian 1C gives EXACTLY 0.
 * gap=1 REDUCTION: K_0 = U V^dag (unitary), Q = 1_d, EXACT full-rank d carrier
 * oracle (and unital then). `out` is aic_ucp_kraus_init'd HERE (dim_K=dim_H=d,
 * r=1; caller aic_ucp_kraus_clears it). See aic_adversarial_carrier_dense.c. */
void aic_adv_chan_carrier_dropout_dense(aic_ucp_kraus *out, slong d, double gap,
                                        slong prec);

/* fam3D — dimension-blowup block algebra (domain.md:416-449; the eps~c/n regime
 * tex:484, the explicit B-diagonal ||D||=1 vs naive Haar error ~ n, tex:1249). A
 * UCP self-map on B(C^N), N = k*d, whose associated eps-C* algebra is
 *   A = (+)_{j=1}^{k} M_d   (dim_A = k*d^2),
 * eta-idempotent with eta TUNABLE by t (t=0 => EXACTLY idempotent). The k-block
 * conditional expectation Phi0(X) = sum_j P_j X P_j (P_j = sum_{a<d} |jd+a><jd+a|,
 * the k diagonal block projectors, sum_j P_j = I_N) has Kraus K_j = P_j (Hermitian),
 * Img Phi0 = (+)_j M_d, and is EXACTLY idempotent + unital. Mixed with its DFT
 * conjugate Ad_{V^dag} Phi0 Ad_V (V[a,b] = exp(2 pi i a b / N)/sqrt(N), the DFT —
 * unitary, NON-block-diagonal so eta>0 for t>0): the output Kraus UNION is
 *   { sqrt(1-t) P_j : j } U { sqrt(t) V^dag P_j V : j },   rank r = 2k.
 * UNITAL for all t: sum_a K_a^dag K_a = (1-t) sum_j P_j + t V^dag(sum_j P_j)V =
 * (1-t)I + t I = I (the make_mixconj pattern of test_idemp.h generalized to k equal
 * blocks, built SELF-CONTAINED here). KNOB t in [0,1) (asserted; k>=1, d>=1,
 * k*d>=2; Rule 4, fail loud). Use k>=2 for the blowup property: k=1 gives a single
 * block P_0=I_N, so V^dag P_0 V = I_N = base and eta=0 for all t (degenerate, no
 * t-dependence); k=1 is allowed only as a well-defined single M_d.
 * Adversarial property (the dimension-blowup trap, tex:484): as the block count k
 * grows, dim_A = k*d^2 grows; the th_main constant c_0 must stay DIMENSION-
 * INDEPENDENT (the explicit generalized-Pauli B-diagonal route, ||D||=1; NOT the
 * naive Haar route whose error ~ n). At t=0 the recovered structure is EXACTLY
 * (+)_j M_d: aic_cstar_build gives B.num_blocks == k, every B.d[l] == d,
 * B.dim_B == k*d^2 (the named structural pin; the block COUNT tracks k).
 * t=0 REDUCTION (exact-idempotent oracle, cross-check ladder #3): out = Phi0, the
 * k-block conditional expectation, EXACTLY idempotent, defect 0; bracket [0,0].
 * `out` is aic_ucp_kraus_init'd HERE (dim_K=dim_H=N=k*d, r=2k; caller clears it). */
void aic_adv_chan_blockalg(aic_ucp_kraus *out, slong k, slong d, double t,
                           slong prec);

/* fam-NC — NON-COMMUTATIVE eta->1/4 boundary (docs/adversarial/README.md:57, the
 * tranche-3 "non-commutative calibrated eta->1/4" item; tex:347-349 cb-norm
 * ampliation-invariance, tex:366-388/378 the measure-prepare closed form,
 * tex:516-525 the theta(2Phi-1) basin). A UCP self-map on B(C^N), N = m*d, that
 * is CALIBRATED to the eta < 1/4 boundary with a GENUINELY non-commutative image,
 * unlike the COMMUTATIVE depolarizing 2A instance (aic_adv_chan_depol_boundary,
 * image = scalars):
 *   Phi = id_{M_m} (x) Psi_eta,   Kraus  1_m (x) K_a^{Psi},   r = d,
 * where id_{M_m} is the IDENTITY channel on B(C^m) (exactly idempotent, fixed-
 * point algebra = the full NON-ABELIAN M_m for m >= 2) and Psi_eta is the
 * measure-prepare cb-norm example on B(C^d) (the family-1B generator
 * aic_adv_chan_cb_op_gap, tex:366-388) with the CLOSED-FORM defect
 * ||Psi_eta^2-Psi_eta||_cb = eta*sqrt(1-eta) (tex:378). cb-norm is AMPLIATION-
 * INVARIANT (tex:347-349), so
 *   eta_cb := ||Phi^2-Phi||_cb = ||Psi^2-Psi||_cb = eta*sqrt(1-eta)  (EXACT, m-indep).
 * KNOB kappa in (0, 1/4): eta is the LOWER root of eta*sqrt(1-eta) = 1/4 - kappa
 * (bisection on the rising branch [0,2/3]; the peak is 0.3849 so 1/4-kappa is
 * reachable for all kappa>0). UNITAL for all kappa (tensor of unital UCP maps).
 * Adversarial property (the NON-COMMUTATIVE boundary, distinct from G2): eta_cb =
 * 1/4 - kappa EXACTLY at the theta(2Phi-1) basin edge, AND the image is NON-
 * ABELIAN — the genuine fixed points B1 = E^{(m)}_{01}(x)1_d, B2 = E^{(m)}_{10}(x)1_d
 * satisfy the star-product (tex:342) commutator ||B1*B2-B2*B1||_op = 1 != 0
 * (the depolarizing 2A instance's abelian image gives 0). Asserts m>=2, d>=2,
 * 0 < kappa < 1/4 (Rule 4, fail loud: target >= 1/4 is OUT of the eta<1/4
 * hypothesis; target > peak 0.3849 has NO calibration root). `out` is
 * aic_ucp_kraus_init'd HERE (dim_K=dim_H=N=m*d, r=d; caller aic_ucp_kraus_clears
 * it). See aic_adversarial_noncomm.c for the full derivation. */
void aic_adv_chan_noncomm_boundary(aic_ucp_kraus *out, slong m, slong d,
                                   double kappa, slong prec);

/* fam2B — RANK-1 defect CONCENTRATED on a NEAR-DEGENERATE subspace
 * (docs/adversarial/domain.md:249-288; the O(sqrt eta) Phi_assoc1 / W
 * cancellation stressor of th_almost_idemp, tex:2192-2204, tex:2296,
 * tex:2385-2422). A measure-prepare UCP self-map on B(C^N), N=d>=3, in the
 * OBSERVABLE convention Phi(X)=sum_a K_a^dag X K_a:
 *   v = (sqrt(1-tail), sqrt(w1), sqrt(w2), 0,...,0)   (unit, on span{0,1,2}),
 *   K_0 = |v><0|  (column 0 = v),   K_j = |j><j|  (j = 1 .. d-1),
 * UNITAL (<v|v>=1) and CP (rank-1 Kraus) for ALL knobs — the family-1B Kraus
 * pattern (aic_adv_chan_cb_op_gap) with a 3D measured vector. Defect:
 *   Phi^2-Phi = P_0 (x) <rho_sub, .>,   EXACTLY superoperator-RANK-1,
 *   rho_sub = -tail|v><v| + w1|1><1| + w2|2><2|,  tail = w1+w2 = 1-|v_0|^2.
 * CALIBRATION (knob eta = TARGET cb defect; closed form + bisection root-find,
 * mirroring aic_adv_chan_noncomm_boundary): tail is the lower root of
 * tail*sqrt(1-tail) = eta on [0,2/3], so eta_cb := ||Phi^2-Phi||_cb = eta
 * (tex:378, a function of tail only). The split w2 = tail*gap_sub, w1 = tail-w2
 * makes |2> the NEAR-DEGENERATE measured direction: gap_sub in (0,1] is the
 * RELATIVE weight of |2> in rho_sub, and as gap_sub -> 0 the |2> eigen-direction
 * collapses toward the kernel (realized small weight tail*gap_sub -> 0) — the
 * near-degenerate subspace the W cancellation (tex:2385-2422) is stressed on.
 * eta_cb is UNCHANGED by gap_sub (depends on tail only): gap_sub tunes the
 * SUBSPACE, not the magnitude. WHY THE LITERAL domain.md:255-261 FORMULA IS NOT
 * USED: Phi = Phi_idemp + eta*(1-Phi_idemp)D is NOT CP (certified-negative Choi
 * min-eig ~ -eta*c/2 for every eta>0; FINDINGS §C18) — this UCP measure-prepare
 * map realizes the SAME named properties honestly. Asserts d>=3,
 * 0 < eta <= 0.3849 (above the peak there is NO calibration root — fail loud,
 * Rule 4), 0 < gap_sub <= 1. eta -> 0 REDUCTION: v -> |0>, Phi -> complete
 * dephasing, EXACTLY idempotent, defect 0 (bracket [0,0]). `out` is
 * aic_ucp_kraus_init'd HERE (caller clears with aic_ucp_kraus_clear). See
 * aic_adversarial_domain2.c for the full derivation. */
void aic_adv_chan_conc_defect(aic_ucp_kraus *out, slong d, double eta,
                              double gap_sub, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_ADVERSARIAL_H */
