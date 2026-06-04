/* aic_funcalc.h — matrix functional calculus: sgn, abs, theta, x^alpha, and the
 * prop_P projection-extractor (beads aic-4bg, "T-funcalc"). Over FLINT's
 * acb_mat_t, certified arb balls (CLAUDE.md Law 2: the matrix type IS acb_mat_t).
 *
 * Paper provenance (approximate_algebras.tex §3 "Some analytic tools"). The
 * verbatim equations live in each source file's docstring (Law 1). The three
 * workhorse functions are defined at tex:514-516:
 *     |X| = (X^2)^{1/2},   sgn(X) = X (X^2)^{-1/2}     ( ||X^2 - x0^2 I|| < x0^2 )
 * theta is tex:527-528:  theta(X) = (I + sgn(X))/2,  and prop_P (tex:524-533)
 * extracts an exact idempotent Ptilde = theta(2P - I) from a near-idempotent P.
 *
 * Validity domain (tex:516, tex:520, also the CLAUDE.md hallucination callout).
 * |X|, sgn(X), theta(X) need ||X^2 - I|| < 1 in this module's normalisation
 * (x0 = 1); the general holomorphic-calculus condition (tex:550) is
 * spec(X^2) subset C \ (-inf,0], i.e. X has no purely-imaginary eigenvalue.
 * Every routine here ASSERTS its precondition and aborts loudly (CLAUDE.md
 * Rule 4) when it is violated; an out-of-domain call is a bug, not a
 * recoverable error.
 *
 * THE AUDITION (CLAUDE.md Law 4). sgn(X) is the core primitive. We implement two
 * independent, eig-free, iterative candidates and benchmark them head-to-head
 * (bench/bench_funcalc.c); both are kept alive, the default aic_sgn dispatches to
 * the benchmarked winner (see aic_funcalc_sgn.c for the selection comment):
 *   - aic_sgn_newton_schulz  : inverse-free Newton-Schulz, Y_{k+1}=1/2 Y_k(3I-Y_k^2)
 *   - aic_sgn_denman_beavers : Newton sign iter, Y_{k+1}=1/2(Y_k+Y_k^{-1}),
 *                              one matrix inverse per step
 * Both are eig-free, so they handle the project's real case: a DEGENERATE
 * (projection) spectrum, where FLINT's certified simple-eig path (aic_mat,
 * bead aic-w4o.1) cannot isolate the repeated eigenvalues. The Hermitian-eig
 * route is therefore NOT a viable third sgn candidate on the inputs we care
 * about; that is recorded as the audition's finding, not forced into code.
 *
 * All routines take an explicit prec and write certified acb_mat_t balls; the
 * cross-check ladder (bead aic-5ty) runs them at prec=53 vs prec=256
 * (tests/test_funcalc.c).
 */
#ifndef AIC_FUNCALC_H
#define AIC_FUNCALC_H

#include <flint/acb_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/* sgn(X) via the inverse-free Newton-Schulz iteration (tex:518, sgn(X)^2 = I).
 *   Y_0 = X,   Y_{k+1} = (1/2) Y_k (3 I - Y_k^2).
 * Converges quadratically to sgn(X) when ||I - X^2|| < 1 (tex:521). ASSERTS that
 * precondition (fail loud). `out` and `X` are n x n; `out` must be initialised
 * (may alias nothing; out != X is required and asserted). Iterates until the
 * step ||Y_{k+1} - Y_k|| stops shrinking below a prec-appropriate tol or a hard
 * cap is hit (the cap is a fail-loud abort: non-convergence is a bug). */
void aic_sgn_newton_schulz(acb_mat_t out, const acb_mat_t X, slong prec);

/* sgn(X) via the CUBIC (order-3) inverse-free Newton-Schulz iteration — the next
 * member of the NS / Pade-of-(I-t)^{-1/2} family (t = I - Y^2):
 *   Y_0 = X,   Y_{k+1} = (1/8) Y_k (15 I - 10 Y_k^2 + 3 Y_k^4).
 * CUBIC convergence (vs quadratic for aic_sgn_newton_schulz) when ||I - X^2|| < 1
 * (tex:521, same asserted precondition), at the cost of 3 matmuls/step (Y^2, Y^4,
 * Y*poly) vs 2. Auditioned against the quadratic NS and Denman-Beavers (Law 4):
 * cubic trades fewer iterations for more work/iter, a Pareto candidate whose win
 * regime is decided by bench/bench_funcalc.c. Inverse-free. Same midpoint-iterate
 * + a-posteriori sign certificate as the quadratic NS (bead aic-1vp). out != X. */
void aic_sgn_newton_schulz3(acb_mat_t out, const acb_mat_t X, slong prec);

/* sgn(X) via the Denman-Beavers / Newton sign iteration (shard-B L50, single-
 * variable form):
 *   Y_0 = X,   Y_{k+1} = 1/2 (Y_k + Y_k^{-1}).
 * Y_k -> sgn(X) quadratically. One certified acb_mat_inv per step (the dominant
 * cost vs Newton-Schulz). Same precondition ||I - X^2|| < 1, asserted. For an
 * exact sign matrix (X^2 = I) Y_1 = X already, so the degenerate-projector case
 * converges in one step. `out` n x n, initialised. */
void aic_sgn_denman_beavers(acb_mat_t out, const acb_mat_t X, slong prec);

/* sgn(X) via the GLOBALLY-convergent Newton sign iteration (bead aic-8hz):
 *   Y_0 = X,   Y_{k+1} = 1/2 (Y_k + Y_k^{-1}),
 * the same step as aic_sgn_denman_beavers, but the precondition is the SPECTRAL
 * basin rho(I - X^2) < 1 (a SUFFICIENT sub-case of Higham's full no-imaginary-
 * eigenvalue hypothesis — Functions of Matrices, SIAM 2008, Ch.5 Thm 5.6, global +
 * quadratic convergence), NOT the operator-norm basin ||X^2-I||_op < 1. It does
 * NOT cover every Higham-valid X (e.g. diag(2,-3) aborts), by design: our inputs
 * X=2 S_Phi - 1 have spectrum near +/-1, exactly this disk. The condition is
 * certified EIG-FREE and RIGOROUSLY via the Gelfand bound: accept the first k with
 * ||(I-X^2)^k||_F^{1/k} certainly < 1 (fail loud if no k <= k_max certifies).
 * After convergence an a-posteriori certificate (||Y^2-I|| < tol AND
 * ||XY-YX|| < tol) pins Y = sgn(X); a straddling certificate ball aborts (Rule
 * 4). Reaches the full non-normal regime the op-basin path rejected (e.g.
 * X=[[a,c],[0,-b]] with large c, where ||I-X^2||_op > 1 but rho(I-X^2) < 1).
 * `out` n x n, initialised, out != X (asserted). */
void aic_sgn_newton_global(acb_mat_t out, const acb_mat_t X, slong prec);

/* Default sgn(X): AUTO-DISPATCHES (bead aic-8hz) — the fast inverse-free
 * Newton-Schulz when ||X^2-I||_op is certified < 1 (the in-basin path, BYTE-FOR-
 * BYTE the previous behavior), else the globally-convergent aic_sgn_newton_global
 * (a strict improvement: the out-of-basin input previously ABORTED). The probe is
 * aic_funcalc_int_in_op_basin (non-aborting). See src/aic_funcalc_sgn.c. */
void aic_sgn(acb_mat_t out, const acb_mat_t X, slong prec);

/* x^alpha of a matrix A near x0 I, via the binomial Taylor series (tex:503-511):
 *   A^alpha = x0^alpha sum_{n>=0} C(alpha,n) (A/x0 - I)^n,
 * which converges (tex:511) when ||A - x0 I|| < x0 (x0 > 0 real). ASSERTS that
 * precondition. Used with alpha = +1/2 (|X| from X^2) and alpha = -1/2 (sgn from
 * X (X^2)^{-1/2}); a general real alpha is supported. `out`, `A` are n x n;
 * `out` initialised; x0 > 0. The series is truncated when the tail bound
 * |C(alpha,n)| ||A/x0 - I||^n (tex:509) drops below a prec-appropriate tol;
 * stagnation without reaching it is a fail-loud abort. */
void aic_funcalc_xpow(acb_mat_t out, const acb_mat_t A, double alpha, double x0,
                      slong prec);

/* |X| = (X^2)^{1/2} (tex:514). Computed as aic_funcalc_xpow(X^2, 1/2, 1).
 * Precondition ||X^2 - I|| < 1, asserted. `out` n x n, initialised. */
void aic_abs(acb_mat_t out, const acb_mat_t X, slong prec);

/* theta(X) = (I + sgn(X)) / 2 (tex:527-528). Uses the default aic_sgn.
 * Precondition ||I - X^2|| < 1, asserted (inside aic_sgn). `out` n x n. */
void aic_theta(acb_mat_t out, const acb_mat_t X, slong prec);

/* prop_P (tex:524-533): given a near-idempotent P with ||P^2 - P|| <= delta < 1/4,
 * return the exact idempotent Ptilde = theta(2P - I). ASSERTS the precondition
 * ||P^2 - P|| < 1/4 (equivalently ||(2P-I)^2 - I|| = 4||P^2-P|| < 1), fail loud.
 * Then Ptilde^2 = Ptilde, Ptilde commutes with P, ||Ptilde - P|| <= ||2P-I|| O(delta).
 * `out` n x n, initialised. */
void aic_prop_P(acb_mat_t out, const acb_mat_t P, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_FUNCALC_H */
