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
 * CATALOGUE PROVENANCE (Law 1). Each generator cites its docs/adversarial/nla.md
 * family/sub-section and, where the catalogue cites the paper, the
 * approximate_algebras.tex line. The seven generators span the distinct NLA attack
 * modes plus the merged lethal shortlist (README.md "Most lethal instances"):
 *   1 jordan_t13        nla 1a, tex:540-544  : non-Lipschitz t^{1/3} spectral burst
 *   2 nonnormal_tri     nla 2c               : departure-from-normality / op-vs-rho gap
 *   3 degenerate_proj   nla 4b (lethal #1)   : exact {0,1} projector, eig-free sgn only
 *   4 near_degen_herm   nla 4a               : tunable Hermitian gap -> simple-eig kill
 *   5 graded_diag       nla 5c               : dynamic range = condition number kappa
 *   6 boundary_x2I      nla 7a (lethal #5)   : ||X^2-I||=1 straddle, funcalc domain edge
 *   7 propP_delta       nla 7b               : ||P^2-P||=delta, prop_P basin edge 1/4
 *
 * API SHAPE (reusable by tests AND benches). Every generator takes a caller-managed
 * acb_mat_t `out` that the caller has init'd to the right shape (asserted), the
 * knob(s), and a prec. Outputs are written into `out`; the caller clears it. Names
 * are aic_adv_<family>. No global state.
 */
#ifndef AIC_ADVERSARIAL_H
#define AIC_ADVERSARIAL_H

#include <flint/acb_mat.h>

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

#ifdef __cplusplus
}
#endif

#endif /* AIC_ADVERSARIAL_H */
