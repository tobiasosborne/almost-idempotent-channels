/* aic_contraction.h — generic Banach contraction / inverse-function solver
 * (beads aic-n3n, "T-contraction"). Over FLINT's acb_mat_t a Banach-space
 * element is a column vector or matrix (CLAUDE.md Law 2: the element type IS
 * acb_mat_t).
 *
 * Paper provenance (approximate_algebras.tex:564-592, Lemma lem_invfun). The
 * lemma's PROOF IS THE ALGORITHM (CLAUDE.md THE MANDATE, the lem_invfun bullet
 * at CLAUDE.md:69-72). Verbatim statement and proof (tex:564-592):
 *
 *   Let V: A -> B be a linear isomorphism of Banach spaces and
 *   f: B_r^A(x0) -> B continuously differentiable such that
 *       || V^{-1} d_x f(x) - 1 || <= c     (x in B_r^A(x0)),   0 <= c < 1.
 *   Then (1) f is injective and
 *       (1-c)||x1-x2|| <= || V^{-1}(f(x1)-f(x2)) || <= (1+c)||x1-x2||;
 *   (2) the image of f includes f(x0) + V(B_{(1-c)r}^A(0)).
 *
 *   Proof. For each given y, the equation f(x)=y is equivalent to x being a
 *   fixed point of  g_y(x) = x + V^{-1}(y - f(x))  (x in B_r^A(x0)). The norm
 *   of the derivative of g_y is at most c, so ||g_y(x1)-g_y(x2)|| <= c||x1-x2||.
 *   [Part 2:] if ||V^{-1}(y - f(x0))|| < (1-c)r then g_y maps B_r^A(x0) to
 *   itself, since ||g_y(x)-x0|| <= ||g_y(x)-g_y(x0)|| + ||g_y(x0)-x0||
 *   < c||x-x0|| + (1-c)r <= r. The fixed point is the limit of x_n=g_y(x_{n-1}),
 *   a Cauchy sequence because ||x_n - x_{n-1}|| < r c^{n-1}.
 *
 * Constructive route. We implement g_y(x) = x + V^{-1}(y - f(x)) literally and
 * iterate it to a fixed point. f and V^{-1} are caller CALLBACKS over acb_mat_t,
 * so the solver is generic in the nonlinear map and the linear preconditioner;
 * the lemma's c<1 is the convergence rate, and the proof's two structural
 * facts — g_y maps the ball to itself, and ||x_n-x_{n-1}|| < r c^{n-1} — are the
 * loud invariants we assert at every step (CLAUDE.md Rule 4).
 *
 * THE AUDITION (CLAUDE.md Law 4). Two candidates for the fixed-point solve, both
 * kept alive, benchmarked head-to-head (bench/bench_contraction.c):
 *   - aic_contraction_picard   : the paper's plain iteration x_n=g_y(x_{n-1}),
 *     linearly convergent at rate c (faithful to tex:591).
 *   - aic_contraction_anderson : Anderson acceleration depth 1 (Type-II/Pulay)
 *     on the SAME g_y; mixes the last two residuals to damp the linear error,
 *     reaching tol in fewer g_y evaluations when c is not small.
 * The default aic_contraction_solve dispatches to the benchmarked winner (see
 * the selection comment in src/aic_contraction.c). Anderson is NOT presumed
 * better: it is auditioned, and where the contraction is already fast (small c)
 * Picard's lower per-step cost can win.
 *
 * Fail-loud (CLAUDE.md Rule 4). The solver aborts, never silently returns, when:
 *   - the supplied or estimated c is not certainly < 1 (lemma hypothesis
 *     violated; tex:569);
 *   - an iterate leaves the ball B_r(x0) (the proof's self-map property fails,
 *     so a hypothesis is wrong — tex:586-590);
 *   - no convergence within max_iter (a c<1 contraction must converge at rate c;
 *     hitting the cap means c was mis-estimated or the input is ill-posed).
 *
 * All routines take an explicit prec; the cross-check ladder (bead aic-5ty) runs
 * the solver at prec=53 vs prec=256 (tests/test_contraction.c).
 */
#ifndef AIC_CONTRACTION_H
#define AIC_CONTRACTION_H

#include <flint/acb_mat.h>
#include <flint/arb.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The nonlinear map f: applies y_out = f(x). `ctx` carries problem data (V, the
 * nonlinearity scale, ...). Must write a matrix of x's shape into out. */
typedef void (*aic_contraction_f)(acb_mat_t out, const acb_mat_t x,
                                  void *ctx, slong prec);

/* The linear map V^{-1}: applies out = V^{-1}(r). Same shape in and out. This is
 * the lemma's V^{-1}, the fixed preconditioner that makes g_y a contraction. */
typedef void (*aic_contraction_Vinv)(acb_mat_t out, const acb_mat_t r,
                                     void *ctx, slong prec);

/* Solver configuration. The Banach element shape is taken from x0. */
typedef struct {
    aic_contraction_f    f;        /* the nonlinear map f(x)                  */
    aic_contraction_Vinv apply_Vinv; /* the linear map V^{-1}(r)              */
    void                *ctx;      /* opaque data passed to both callbacks    */
    /* x0/y are stored as acb_mat_struct* (the type acb_mat_t decays to), so a
     * caller assigns the matrix NAME directly (o.x0 = x0), not &x0 — taking the
     * address of the array type acb_mat_t and storing it in a const-pointer
     * field trips -Wpedantic's qualifier check pre-C2x. */
    const acb_mat_struct *x0;      /* ball centre (lemma's x0); shape source  */
    const acb_mat_struct *y;       /* target: solve f(x) = y                  */
    double               r;        /* ball radius (lemma's r); > 0            */
    double               c;        /* contraction bound c in [0,1); tex:569   */
    double               tol;      /* stop when midpoint step < tol           */
    slong                max_iter; /* hard cap; exceeding it is a fail-loud   */
    slong                prec;     /* arb working precision                   */
} aic_contraction_opts;

/* Diagnostics returned by every solver (caller-supplied, may be NULL). */
typedef struct {
    slong iters;        /* g_y evaluations performed before tol              */
    double last_step;   /* final midpoint ||x_n - x_{n-1}|| (Frobenius)      */
    double max_step;    /* max over n of ||x_n - x_{n-1}|| (for the Cauchy   */
                        /*   bound check; <= r c^{n-1} at each step)         */
} aic_contraction_stats;

/* Picard iteration x_n = g_y(x_{n-1}), g_y(x) = x + V^{-1}(y - f(x))
 * (tex:580-591, the paper's method). Writes the fixed point into x (initialised
 * to the result shape). Asserts c<1, ball containment, and the r c^{n-1} Cauchy
 * bound at each step; aborts on cap. Fills *st if non-NULL. */
void aic_contraction_picard(acb_mat_t x, const aic_contraction_opts *o,
                            aic_contraction_stats *st);

/* Anderson-accelerated (depth 1) solve of the SAME fixed point of g_y. Mixes the
 * last two residuals r_k = g_y(x_k) - x_k to damp the linear-rate error. Same
 * fail-loud guards (c<1, ball containment, cap). Fills *st if non-NULL. */
void aic_contraction_anderson(acb_mat_t x, const aic_contraction_opts *o,
                              aic_contraction_stats *st);

/* Default solve: dispatches to the benchmarked-best candidate (see the selection
 * comment in src/aic_contraction.c, citing bench/bench_contraction.c). */
void aic_contraction_solve(acb_mat_t x, const aic_contraction_opts *o,
                           aic_contraction_stats *st);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CONTRACTION_H */
