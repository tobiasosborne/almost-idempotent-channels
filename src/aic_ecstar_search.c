/* aic_ecstar_search.c — the double-path restart ENGINE for the Cycle-2 faithful
 * worst-case defect search (bead aic-knm, module "ecstar"). It runs the HOPM
 * inner loop to FIND a worst-case witness; the public API and the arb
 * CERTIFICATION of that witness live in aic_ecstar_certify.c (Rule 10 split). The
 * matrix kernel is aic_ecstar_hopm.c; the block update + restart init is
 * aic_ecstar_iterate.c. See include/aic_ecstar.h for the algorithm contract and
 * the load-bearing subspace-polar accept-guard.
 *
 * objval() computes the active defect map: assoc -> h(X,Y,Z) (nblk==3);
 * submult -> X*Y (nblk==2); cstar -> X^dag*X (nblk==1). The scale-invariant
 * objective tracked across the accept-guard is ||objval||_op (every accepted
 * block is op-norm 1, so the ratio reduces to this norm). For cstar we MINIMIZE.
 */
#include <complex.h>
#include <math.h>

#include "aic_ecstar_search.h"

static void copyn(cx *dst, const cx *src, slong nn)
{
    for (slong p = 0; p < nn; p++) dst[p] = src[p];
}

/* "is `a` better than `b`?" maximize: a>b; minimize: a<b. */
static int better(double a, double b, int minimize)
{
    return minimize ? (a < b) : (a > b);
}

/* g = the defect map at (X,Y,Z): nblk==3 assoc; ==2 submult X*Y; ==1 cstar X^dag X. */
void aic_ehk_objval(cx *g, const aic_ehk *h, const cx *X, const cx *Y,
                    const cx *Z, int nblk, aic_ehk_work *w)
{
    slong n = h->n;
    if (nblk == 3) {
        aic_ehk_assoc(g, h, X, Y, Z, w->s1, w->s2, w->s3, w->s4, w->s5);
    } else if (nblk == 2) {
        aic_ehk_star(g, h, X, Y, w->s1, w->s2, w->s3);
    } else { /* cstar: X^dag * X */
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
                w->s6[i * n + j] = conj(X[j * n + i]); /* X^dag */
        aic_ehk_star(g, h, w->s6, X, w->s1, w->s2, w->s3);
    }
}

double aic_ehk_run_search(const aic_ehk *h, int nrest, int niter, int nblk,
                          int minimize, aic_ehk_term tX, aic_ehk_term tY,
                          aic_ehk_term tZ, cx *wX, cx *wY, cx *wZ)
{
    slong n = h->n, d = h->d, nn = n * n;
    aic_ehk_work w;
    aic_ehk_work_alloc(&w, n, d);
    cx *X = flint_malloc((size_t) nn * sizeof(cx));
    cx *Y = flint_malloc((size_t) nn * sizeof(cx));
    cx *Z = flint_malloc((size_t) nn * sizeof(cx));
    cx *cand = flint_malloc((size_t) nn * sizeof(cx));
    cx *g = flint_malloc((size_t) nn * sizeof(cx));
    cx *u = flint_malloc((size_t) n * sizeof(cx));
    cx *v = flint_malloc((size_t) n * sizeof(cx));
    cx *tmp = flint_malloc((size_t) n * sizeof(cx));
    double sign = minimize ? -1.0 : 1.0;
    double best = minimize ? 1e300 : -1e300;

    for (int ri = 0; ri < nrest; ri++) {
        uint64_t seed = aic_ehk_seed(ri);
        int warm = (ri % 2) == 0;
        aic_ehk_init_block(X, h, &seed, warm, &w);
        aic_ehk_init_block(Y, h, &seed, warm, &w);
        aic_ehk_init_block(Z, h, &seed, warm, &w);
        aic_ehk_init_uv(u, v, n, &seed);

        for (int it = 0; it <= niter; it++) {
            aic_ehk_objval(g, h, X, Y, Z, nblk, &w);
            double cur = aic_ehk_opnorm(g, n);
            if (it == niter) break;            /* last sample taken above */
            aic_ehk_uv_step(u, v, g, n, tmp);

            w.frz1 = (nblk == 1) ? X : (nblk >= 2 ? Y : NULL);
            w.frz2 = (nblk == 3) ? Z : NULL;
            if (aic_ehk_block_candidate(cand, &w, h, u, v, sign, tX)) {
                aic_ehk_objval(g, h, cand, Y, Z, nblk, &w);
                double nc = aic_ehk_opnorm(g, n);
                if (better(nc, cur, minimize)) { copyn(X, cand, nn); cur = nc; }
            }
            if (nblk >= 2) {
                w.frz1 = X; w.frz2 = (nblk == 3) ? Z : NULL;
                if (aic_ehk_block_candidate(cand, &w, h, u, v, sign, tY)) {
                    aic_ehk_objval(g, h, X, cand, Z, nblk, &w);
                    double nc = aic_ehk_opnorm(g, n);
                    if (better(nc, cur, minimize)) { copyn(Y, cand, nn); cur = nc; }
                }
            }
            if (nblk == 3) {
                w.frz1 = X; w.frz2 = Y;
                if (aic_ehk_block_candidate(cand, &w, h, u, v, sign, tZ)) {
                    aic_ehk_objval(g, h, X, Y, cand, nblk, &w);
                    double nc = aic_ehk_opnorm(g, n);
                    if (better(nc, cur, minimize)) { copyn(Z, cand, nn); cur = nc; }
                }
            }
        }
        /* final objective at the converged (X,Y,Z) of this restart */
        aic_ehk_objval(g, h, X, Y, Z, nblk, &w);
        double cur = aic_ehk_opnorm(g, n);
        if (better(cur, best, minimize)) {
            best = cur;
            copyn(wX, X, nn);
            if (nblk >= 2) copyn(wY, Y, nn);
            if (nblk == 3) copyn(wZ, Z, nn);
        }
    }
    flint_free(tmp); flint_free(v); flint_free(u); flint_free(g);
    flint_free(cand); flint_free(Z); flint_free(Y); flint_free(X);
    aic_ehk_work_free(&w);
    return best;
}
