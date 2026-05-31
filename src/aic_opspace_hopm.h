/* aic_opspace_hopm.h — internal interface for the opspace operator-norm ampliated
 * max-stretch HOPM (bead aic-zwo). NOT a public header (lives in src/): the public
 * contract is aic_opspace_forward_stretch / aic_opspace_inverse_stretch in
 * include/aic_opspace.h. This shares the `opmap` (a level-1 linear map f: U -> V in
 * the double path) and the HOPM kernel between aic_opspace_ampliate.c (the kernel)
 * and aic_opspace_map.c (the v-specific opmap builders + public entry points), so
 * neither file exceeds the ~200 LOC core limit (CLAUDE.md Rule 10).
 *
 * See aic_opspace_ampliate.c for the map abstraction, the ampliation, and the
 * scale-invariant polar-then-project HOPM (the operator-norm max-stretch, FINDINGS
 * §C12 — NOT the vacuous Frobenius sigma_min).
 */
#ifndef AIC_OPSPACE_HOPM_H
#define AIC_OPSPACE_HOPM_H

#include <complex.h>

#include <flint/flint.h> /* slong */

typedef double _Complex cx;

/* A level-1 linear map f: U -> V in the double path. Frobenius-ORTHONORMAL bases
 * {U[a]} (a in [0,dU), each sU x sU) of U and {W[b]} (b in [0,dV), each sV x sV)
 * of V, plus the coordinate matrix F[b*dU + a] = <W[b], f(U[a])>_F. All buffers
 * OWNED (freed by aic_opmap_clear). */
typedef struct {
    slong dU, dV;   /* dim U == dim V (bijective)                          */
    slong sU, sV;   /* U <= M_{sU}, V <= M_{sV}                            */
    cx **U;         /* dU domain ON basis ops, each sU*sU row-major         */
    cx **W;         /* dV codomain ON basis ops, each sV*sV row-major       */
    cx *F;          /* dV*dU coordinate matrix, F[b*dU + a]                 */
} opmap;

/* row-major double _Complex buffer of length nn. */
cx *aic_opmap_cxalloc(slong nn);

/* Free the U/W bases and F. */
void aic_opmap_clear(opmap *m);

/* The block-ampliated maps (src/aic_opspace_apply.c). 1_{M_n}⊗f maps the (n*sU)
 * x (n*sU) block matrix X (block (k,l) in U) to the (n*sV) x (n*sV) matrix with
 * block (k,l) = f(X_{kl}); the HS adjoint 1_{M_n}⊗f^dag maps blockwise by f^dag;
 * project_amp projects an ambient (cU x cU) matrix onto M_n⊗U block-by-block.
 * Scratch buffers caller-provided (sizes in the .c docstrings). These are the
 * mathematical core of th_main_ext; the structural-ampliation cross-check tooth
 * (test_opspace.c, FINDINGS §C12) pins aic_opmap_apply_fn against an independent
 * Kronecker reference so a "(0,0)-block-only" cripple is caught (Rule 5). */
void aic_opmap_apply_fn(cx *Y, const opmap *m, const cx *X, slong n, cx *xb,
                        cx *yb);
void aic_opmap_apply_fn_dag(cx *G, const opmap *m, const cx *P, slong n, cx *pb,
                            cx *gb);
void aic_opmap_project_amp(cx *Pol, const opmap *m, slong n, cx *pb, cx *gb);

/* Scale-invariant HOPM max-stretch ||f_n||_op = sup_{X != 0 in M_n (x) U}
 * ||(1_{M_n} (x) f)(X)||_op / ||X||_op, a rigorous LOWER bound on the true sup
 * (honest direction; see include/aic_opspace.h). nrest restarts, niter inner
 * iterations, both deterministic. */
double aic_opmap_hopm_stretch(const opmap *m, slong n, int nrest, int niter);

/* HOPM tuning: enough restarts/iters to reach the isometric optimum on the small
 * ampliated sizes here (validated by the eta=0 oracle == 1 exactly, test T1). */
#define OPSPACE_NREST 6
#define OPSPACE_NITER 40

#endif /* AIC_OPSPACE_HOPM_H */
