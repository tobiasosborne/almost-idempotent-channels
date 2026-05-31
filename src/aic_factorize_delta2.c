/* aic_factorize_delta2.c — the UNITALIZATION of the F2 encode map: Delta'(I_B)^{-1/2}
 * and Delta(X) = Delta'(I_B)^{-1/2} Delta'(X) Delta'(I_B)^{-1/2} (bead aic-tff,
 * module "factorize", Increment F2). Split from aic_factorize_delta.c (Rule 10,
 * <=200 LOC). See include/aic_factorize.h (the F2 banner).
 *
 * approximate_algebras.tex:2799 — th_factorization Step 4 (verbatim):
 *   "Delta: X |-> (Delta'(I_B))^{-1/2} Delta'(X) (Delta'(I_B))^{-1/2}
 *    is a UCP map such that ||Delta - Delta~||_cb <= O(eta)."
 *
 * Delta is UCP: UNITAL by construction (Delta(I_B) = M^{-1/2} Delta'(I_B) M^{-1/2}
 * = M^{-1/2} M M^{-1/2} = 1_H with M = Delta'(I_B)), and CP because a congruence
 * X |-> S^dag X S preserves complete positivity (here S = Delta'(I_B)^{-1/2},
 * Hermitian PSD). The inverse-sqrt uses the certified binomial xpow
 * (aic_funcalc_xpow, alpha = -1/2, x0 = 1), whose basin is ||Delta'(I_B) - 1_H||_op
 * < 1 — ASSERTED loud (Rule 4, shard H #7) in aic_factorize_delta_build BEFORE the
 * power, since Delta'(I_B) = Delta~(I_B) + O(eta) ~ 1_H only for small eta.
 */
#include <assert.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_dhom.h"
#include "aic/aic_factorize.h"
#include "aic/aic_funcalc.h"
#include "aic/aic_mat.h"

void aic_factorize_delta_build(aic_factorize *F, slong prec)
{
    assert(F != NULL);
    slong N = F->N;

    /* Delta'(I_B), I_B = the block-diagonal unit (+)_l I_{d_l}. */
    acb_mat_t IB, DI, DImI;
    acb_mat_init(IB, F->n_B, F->n_B);
    acb_mat_init(DI, N, N);
    acb_mat_init(DImI, N, N);
    aic_dhom_B_unit(IB, F->v->B);
    aic_factorize_delta_prime(DI, F, IB, prec);

    /* ASSERT ||Delta'(I_B) - 1_H||_op < 1 (Rule 4, shard H #7 — the inverse-sqrt
     * basin; aic_funcalc_xpow(., -1/2, 1) needs ||A - I||_op < 1). arb_lt is the
     * loud-on-uncertainty form: a straddling ball does NOT silently pass. */
    {
        acb_mat_one(DImI);
        acb_mat_sub(DImI, DI, DImI, prec);
        arb_t e, one;
        arb_init(e);
        arb_init(one);
        aic_mat_opnorm(e, DImI, prec);
        arb_one(one);
        assert(arb_lt(e, one) &&
               "delta_build: ||Delta'(I_B) - 1_H||_op >= 1 (inverse-sqrt out of "
               "basin; shard H #7)");
        arb_clear(one);
        arb_clear(e);
    }

    if (F->delta_ready) acb_mat_clear(F->deltaI_invsqrt);
    acb_mat_init(F->deltaI_invsqrt, N, N);
    aic_funcalc_xpow(F->deltaI_invsqrt, DI, -0.5, 1.0, prec); /* Delta'(I_B)^{-1/2} */
    F->delta_ready = 1;

    acb_mat_clear(DImI); acb_mat_clear(DI); acb_mat_clear(IB);
}

void aic_factorize_delta(acb_mat_t out, const aic_factorize *F,
                         const acb_mat_t X, slong prec)
{
    assert(F != NULL);
    assert(F->delta_ready &&
           "delta: aic_factorize_delta_build must run first");
    slong N = F->N;
    assert(acb_mat_nrows(out) == N && acb_mat_ncols(out) == N &&
           "delta: out must be N x N");

    acb_mat_t DP, tmp;
    acb_mat_init(DP, N, N);
    acb_mat_init(tmp, N, N);
    aic_factorize_delta_prime(DP, F, X, prec);            /* Delta'(X) */
    acb_mat_mul(tmp, F->deltaI_invsqrt, DP, prec);        /* M^{-1/2} Delta'(X) */
    acb_mat_mul(out, tmp, F->deltaI_invsqrt, prec);       /* ... M^{-1/2} */
    acb_mat_clear(tmp); acb_mat_clear(DP);
}
