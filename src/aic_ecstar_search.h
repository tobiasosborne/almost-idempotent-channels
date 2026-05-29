/* aic_ecstar_search.h — internal interface for the HOPM restart engine (bead
 * aic-knm, Cycle 2). Shared by aic_ecstar_search.c (the engine) and
 * aic_ecstar_certify.c (the public API + arb certification). Not public. */
#ifndef AIC_ECSTAR_SEARCH_H
#define AIC_ECSTAR_SEARCH_H

#include "aic_ecstar_iterate.h"

/* a defect-map term thunk: fills `out` with the active defect map when the
 * active block is set to B_k (others read from w->frz1, w->frz2). */
typedef void (*aic_ehk_term)(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k);

/* g = defect map at (X,Y,Z): nblk==3 assoc; ==2 submult X*Y; ==1 cstar X^dag*X. */
void aic_ehk_objval(cx *g, const aic_ehk *h, const cx *X, const cx *Y,
                    const cx *Z, int nblk, aic_ehk_work *w);

/* Run the multi-restart HOPM search; returns the best double objective
 * ||objval||_op and writes the best witness blocks into wX[,wY[,wZ]]. minimize=1
 * for cstar. tY/tZ may be NULL when nblk<3 / nblk<2. */
double aic_ehk_run_search(const aic_ehk *h, int nrest, int niter, int nblk,
                          int minimize, aic_ehk_term tX, aic_ehk_term tY,
                          aic_ehk_term tZ, cx *wX, cx *wY, cx *wZ);

/* scratch lifecycle (in aic_ecstar_certify.c). */
void aic_ehk_work_alloc(aic_ehk_work *w, slong n, slong d);
void aic_ehk_work_free(aic_ehk_work *w);

#endif /* AIC_ECSTAR_SEARCH_H */
