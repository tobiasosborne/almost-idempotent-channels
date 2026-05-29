/* aic_ecstar_iterate.h — internal interface for the HOPM restart/block-update
 * engine (bead aic-knm, Cycle 2). Shares the scratch struct and the block-update
 * primitives between aic_ecstar_iterate.c and the driver aic_ecstar_search.c.
 * Not a public header. */
#ifndef AIC_ECSTAR_ITERATE_H
#define AIC_ECSTAR_ITERATE_H

#include <stdint.h>

#include "aic_ecstar_hopm.h"

/* Per-call scratch: all buffers sized for n x n (or n / d / SVD), allocated once
 * by the driver and reused across restarts. */
typedef struct {
    cx *g;     /* n*n: a defect-map value (h, star, or X^dag*X)  */
    cx *C;     /* n*n: the A-gradient sum_k c_k B_k              */
    cx *pol;   /* n*n: polar factor of C                         */
    cx *U;     /* n*n: SVD left vectors                          */
    cx *Vt;    /* n*n: SVD V^dag                                 */
    cx *vec;   /* n:   a working vector (g v, etc.)              */
    cx *cf;    /* d:   projection coefficients                   */
    double *sv;/* n:   singular values                           */
    /* extra scratch the defect-term thunks need (Phi/star intermediates) */
    cx *s1, *s2, *s3, *s4, *s5, *s6;
    /* the two FROZEN blocks for the current active-block update; a thunk reads
     * these and substitutes B_k into the active slot. For submult/cstar the
     * unused pointer is NULL. */
    const cx *frz1, *frz2;
} aic_ehk_work;

/* Build the A-gradient and the candidate block X'=Pi_A(polar(sign*C))/||.||_op;
 * returns 1 if valid. `term(out,w,h,k)` fills `out` with the defect map when the
 * active block is set to B_k. sign=+1 maximizes, -1 minimizes. */
int aic_ehk_block_candidate(cx *cand, aic_ehk_work *w, const aic_ehk *h,
                            const cx *u, const cx *v, double sign,
                            void (*term)(cx *out, aic_ehk_work *w,
                                         const aic_ehk *h, slong k));

/* Restart-init helpers (deterministic LCG; no wall clock). */
void aic_ehk_init_block(cx *X, const aic_ehk *h, uint64_t *seed, int warm,
                        aic_ehk_work *w);
void aic_ehk_init_uv(cx *u, cx *v, slong n, uint64_t *seed);
uint64_t aic_ehk_seed(int restart_index);

#endif /* AIC_ECSTAR_ITERATE_H */
