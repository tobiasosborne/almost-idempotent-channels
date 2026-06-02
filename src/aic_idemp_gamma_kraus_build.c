/* aic_idemp_gamma_kraus_build.c — the Kraus operator L_j of the prop_Gamma form
 * (bead aic-ynu, I3): build L_j = (W_j^dag (x) 1_{F_j})(1_{L_j} (x) xi_j) from the
 * density matrix gamma_j and its purification xi_j (.tex:2118-2122). Split out of
 * aic_idemp_gamma_kraus_ops.c for the <=200 LOC soft limit (Rule 10). See
 * aic_idemp_gamma_kraus.c for the narrative.
 */
#include <complex.h>
#include <math.h>

#include <flint/flint.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_idemp_gamma_internal.h"

void aic_gk_build_Lj(acb_mat_t L, const acb_mat_t gamma, const acb_mat_t Wj,
                     slong dL, slong dE, slong m, slong j, slong prec)
{
    /* The purification: xi_j[e*dE+f] = S[e,f] with S = gamma_j^{1/2} (Hermitian
     * PSD), so Tr_{F_j}(xi_j xi_j^dag) = S S = gamma_j (S Hermitian). S is built
     * degeneracy-robustly as sum_c sqrt(lambda_c) Pi_c over the certified eigen-
     * clusters (gamma_j is OFTEN maximally mixed = degenerate, e.g. (1/2)1 for
     * noiseless_subsystem(2,2), so the simple-spectrum aic_mat_eig_hermitian would
     * abort — aic_mat_eig_hermitian_subspaces is the degeneracy-robust route). */
    acb_mat_t S;
    acb_mat_init(S, dE, dE);
    acb_mat_zero(S);
    {
        aic_mat_eigcluster *cl;
        slong ncl;
        aic_mat_eig_hermitian_subspaces(&cl, &ncl, gamma, -1.0, prec);
        acb_mat_t Pi;
        acb_mat_init(Pi, dE, dE);
        for (slong c = 0; c < ncl; c++) {
            arb_t lam, sq;
            arb_init(lam);
            arb_init(sq);
            arb_set(lam, acb_realref(cl[c].lambda));
            arb_nonnegative_part(sq, lam);        /* clamp tiny negative */
            arb_sqrt(sq, sq, prec);
            aic_mat_cluster_projector(Pi, &cl[c], prec);
            acb_mat_scalar_mul_arb(Pi, Pi, sq, prec);
            acb_mat_add(S, S, Pi, prec);          /* += sqrt(lambda_c) Pi_c */
            arb_clear(sq);
            arb_clear(lam);
        }
        acb_mat_clear(Pi);
        aic_mat_eigcluster_free(cl, ncl);
    }

    /* 1_L (x) xi_j : (dL*dE*dE) x dL, left-major; xi[e*dE+f] = S[e,f]. */
    acb_mat_t Lxi;
    acb_mat_init(Lxi, dL * dE * dE, dL);
    acb_mat_zero(Lxi);
    for (slong a = 0; a < dL; a++)
        for (slong e = 0; e < dE; e++)
            for (slong f = 0; f < dE; f++)
                acb_set(acb_mat_entry(Lxi, (a * dE + e) * dE + f, a),
                        acb_mat_entry(S, e, f));
    acb_mat_clear(S);

    /* W_j^dag (x) 1_F : (m*dE) x (dL*dE*dE). W_j^dag is m x (dL*dE); left-major
     * tensor over F (=E) -> (W_j^dag (x) 1_F)[r*dE+f, t*dE+f'] = W_j^dag[r,t] d_ff'. */
    acb_mat_t Wjd;
    acb_mat_init(Wjd, m, dL * dE);
    acb_mat_conjugate_transpose(Wjd, Wj);
    acb_mat_t WdxI;
    acb_mat_init(WdxI, m * dE, dL * dE * dE);
    acb_mat_zero(WdxI);
    for (slong r = 0; r < m; r++)
        for (slong t = 0; t < dL * dE; t++)
            for (slong f = 0; f < dE; f++)
                acb_set(acb_mat_entry(WdxI, r * dE + f, t * dE + f),
                        acb_mat_entry(Wjd, r, t));

    acb_mat_mul(L, WdxI, Lxi, prec);              /* L_j = (W^dag(x)1)(1(x)xi)  */

    /* certify L_j^dag L_j = 1_{L_j}. */
    acb_mat_t Ld, G, one;
    arb_t tol, eo;
    acb_mat_init(Ld, dL, m * dE);
    acb_mat_init(G, dL, dL);
    acb_mat_init(one, dL, dL);
    arb_init(tol);
    arb_init(eo);
    arb_set_d(tol, 1e-9);
    acb_mat_conjugate_transpose(Ld, L);
    acb_mat_mul(G, Ld, L, prec);
    acb_mat_one(one);
    acb_mat_sub(G, G, one, prec);
    aic_mat_opnorm(eo, G, prec);
    if (!arb_le(eo, tol)) {
        flint_printf("aic_idemp_gamma_kraus: Kraus L_%wd not an isometry "
                     "(||L^dag L - 1_L|| too large); fail loud (Rule 4)\n", j);
        flint_abort();
    }
    arb_clear(eo);
    arb_clear(tol);
    acb_mat_clear(one);
    acb_mat_clear(G);
    acb_mat_clear(Ld);
    acb_mat_clear(WdxI);
    acb_mat_clear(Wjd);
    acb_mat_clear(Lxi);
}
