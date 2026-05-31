/* aic_opspace_ampliate.c — the OPERATOR-NORM ampliated max-stretch HOPM KERNEL of
 * the opspace module (bead aic-zwo, increment O1). Realizes the operator-norm
 * cb-inclusion quantities of th_main_ext (approximate_algebras.tex:1447-1561, §10;
 * a_n .tex:1489, def:opspace .tex:1453-1464). See include/aic_opspace.h for the
 * full contract; the binding point is FINDINGS §C12: this measures the
 * OPERATOR-norm stretch, NOT the vacuous Frobenius sigma_min(I_{n^2} (x) M_1). The
 * v-specific opmap builders + public entry points are in aic_opspace_map.c (Rule
 * 10 split); the block-ampliation primitives (1_{M_n}(x)f, its adjoint, the
 * project-onto-M_n(x)U) are in aic_opspace_apply.c (a further Rule 10 split, so
 * each file stays under the ~200 LOC core limit); the shared `opmap` + this
 * kernel's interface are in src/aic_opspace_hopm.h.
 *
 * THE MAP ABSTRACTION. A level-1 linear map f: U -> V between two subspaces of
 * matrix algebras is the `opmap` (src/aic_opspace_hopm.h): Frobenius-ORTHONORMAL
 * bases {U[a]} (sU x sU) of U, {W[b]} (sV x sV) of V, and the coordinate matrix
 * F[b,a] = <W[b], f(U[a])>_F. The ampliation 1_{M_n} (x) f and its HS adjoint are
 * the aic_opmap_apply_fn / aic_opmap_apply_fn_dag of aic_opspace_apply.c.
 *
 * THE HOPM (scale-invariant, mirrors aic_ecstar_search.c). Maximize ||f_n(X)||_op
 * over the operator-norm unit ball of M_n (x) U: (i) output power step Y = f_n(X),
 * top singular triple (u,w) of Y via aic_ehk_uv_step; (ii) input step — the
 * gradient G = f_n^dag(u w^*) lives in M_n (x) U; its maximizer over that op-norm
 * unit ball is the polar factor of G, taken in the ambient M_{n*sU} then PROJECTED
 * back onto M_n (x) U (each block onto {U[a]}), accepted only if ||f_n(.)||_op
 * strictly improves (monotone guard). Every accepted iterate is op-normalized and
 * exactly in M_n (x) U, so the result is a rigorous LOWER bound on ||f_n||_op. When
 * M_n (x) U is a genuine C* algebra (forward, M_n (x) B) the exact polar already
 * lies in it (projection a safety no-op); when it is the eps-C* subspace M_n (x) A
 * (inverse) the projection is load-bearing (eps-C* is not polar-closed; FINDINGS
 * §C2 / ecstar accept-guard). Row-major double (aic_latd / aic_ehk convention);
 * deterministic restarts (a fixed counter, no wall-clock RNG; project rule).
 */
#include <math.h>
#include <stdlib.h>

#include "aic_ecstar_hopm.h"
#include "aic_opspace_hopm.h"

cx *aic_opmap_cxalloc(slong nn)
{
    return flint_malloc((size_t) nn * sizeof(cx));
}

void aic_opmap_clear(opmap *m)
{
    for (slong a = 0; a < m->dU; a++) flint_free(m->U[a]);
    for (slong b = 0; b < m->dV; b++) flint_free(m->W[b]);
    flint_free(m->U);
    flint_free(m->W);
    flint_free(m->F);
}

/* op-norm-normalize X (q x q) IN PLACE; return the old op-norm. */
static double opnorm_normalize(cx *X, slong q)
{
    double s = aic_ehk_opnorm(X, q);
    if (s > 1e-300)
        for (slong p = 0; p < q * q; p++) X[p] /= s;
    return s;
}

double aic_opmap_hopm_stretch(const opmap *m, slong n, int nrest, int niter)
{
    slong sU = m->sU, sV = m->sV, cU = n * sU, cV = n * sV;
    cx *X = aic_opmap_cxalloc(cU * cU), *Cand = aic_opmap_cxalloc(cU * cU);
    cx *Y = aic_opmap_cxalloc(cV * cV), *G = aic_opmap_cxalloc(cU * cU);
    cx *Pol = aic_opmap_cxalloc(cU * cU);
    cx *u = aic_opmap_cxalloc(cV), *w = aic_opmap_cxalloc(cV);
    cx *tmp = aic_opmap_cxalloc(cV > cU ? cV : cU);
    cx *Usvd = aic_opmap_cxalloc(cU * cU), *Vt = aic_opmap_cxalloc(cU * cU);
    double *sv = flint_malloc((size_t) cU * sizeof(double));
    cx *xb = aic_opmap_cxalloc(sU * sU), *yb = aic_opmap_cxalloc(sV * sV);
    cx *pb = aic_opmap_cxalloc(sU * sU > sV * sV ? sU * sU : sV * sV);
    cx *gb = aic_opmap_cxalloc(sU * sU);
    double best = 0.0;

    for (int ri = 0; ri < nrest; ri++) {
        /* deterministic init: a structured-but-varied seed (no wall-clock RNG).
         * Build a random-ish element of M_n (x) U via the basis, op-normalize. */
        unsigned long s = 2862933555777941757UL * (unsigned long) (ri + 1)
                          + 3037000493UL;
        for (slong k = 0; k < n; k++)
            for (slong l = 0; l < n; l++) {
                for (slong p = 0; p < sU * sU; p++) gb[p] = 0.0;
                for (slong a = 0; a < m->dU; a++) {
                    s = s * 6364136223846793005UL + 1442695040888963407UL;
                    double re = ((double) ((s >> 33) & 0xffff) / 32768.0) - 1.0;
                    s = s * 6364136223846793005UL + 1442695040888963407UL;
                    double im = ((double) ((s >> 33) & 0xffff) / 32768.0) - 1.0;
                    cx c = re + im * I;
                    for (slong p = 0; p < sU * sU; p++) gb[p] += c * m->U[a][p];
                }
                for (slong i = 0; i < sU; i++)
                    for (slong j = 0; j < sU; j++)
                        X[(k * sU + i) * cU + (l * sU + j)] = gb[i * sU + j];
            }
        if (opnorm_normalize(X, cU) < 1e-300) {
            for (slong p = 0; p < cU * cU; p++) X[p] = 0.0;
            for (slong p = 0; p < cU; p++) X[p * cU + p] = 1.0;
            opnorm_normalize(X, cU);
        }
        for (slong p = 0; p < cV; p++) w[p] = (p == 0) ? 1.0 : 0.0;

        for (int it = 0; it <= niter; it++) {
            aic_opmap_apply_fn(Y, m, X, n, xb, yb);
            double cur = aic_ehk_opnorm(Y, cV);
            if (cur > best) best = cur;
            if (it == niter) break;
            aic_ehk_uv_step(u, w, Y, cV, tmp);          /* top triple (u,w) of Y  */
            for (slong i = 0; i < cV; i++)
                for (slong j = 0; j < cV; j++)
                    Y[i * cV + j] = u[i] * conj(w[j]);   /* reuse Y as P = u w^*   */
            aic_opmap_apply_fn_dag(G, m, Y, n, pb, gb);
            aic_ehk_polar(Pol, G, cU, Usvd, Vt, sv);     /* ambient polar of G     */
            aic_opmap_project_amp(Pol, m, n, pb, gb);    /* project onto M_n(x)U   */
            double pn = opnorm_normalize(Pol, cU);
            if (pn < 1e-300) continue;
            for (slong p = 0; p < cU * cU; p++) Cand[p] = Pol[p];
            aic_opmap_apply_fn(Y, m, Cand, n, xb, yb);
            double nc = aic_ehk_opnorm(Y, cV);
            if (nc > cur)                                 /* monotone-ascent guard */
                for (slong p = 0; p < cU * cU; p++) X[p] = Cand[p];
        }
    }

    flint_free(gb); flint_free(pb); flint_free(yb); flint_free(xb);
    flint_free(sv); flint_free(Vt); flint_free(Usvd);
    flint_free(tmp); flint_free(w); flint_free(u);
    flint_free(Pol); flint_free(G); flint_free(Y);
    flint_free(Cand); flint_free(X);
    return best;
}
