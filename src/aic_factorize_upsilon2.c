/* aic_factorize_upsilon2.c — the UCP DECODE map Upsilon, PART 2: V, xi_j, L_j and
 * the full F3 BUILD (bead aic-tff, module "factorize", Increment F3). Split from
 * aic_factorize_upsilon.c (Rule 10, <=200 LOC); the apply maps Upsilon'_j /
 * Upsilon'/Upsilon + the D5 pin live in aic_factorize_upsilon3.c. See
 * include/aic_factorize.h (the F3 banner) and the LOCKED decisions D4-D6 of
 * docs/research/factorize_f3_synthesis.md.
 *
 * approximate_algebras.tex:2859-2861 — th_factorization Step 5 (verbatim):
 *   xi_j in E_j a unit vector with |  ||C_j xi_j|| - 1 | <= O(eta) (.tex:2859).
 *   Phi(X) = V^dag (X (x) 1_F) V,  V: H -> H (x) F (.tex:2859, Stinespring).
 *   L_j: L_j -> H (x) F (.tex:2861):
 *     L_j = sum_s p_{js} (Delta(U_{js}^dag) (x) 1_F) V W_j^dag (U_{js} (x) xi_j).
 *
 * D5 — THE F-ANCILLA ORDERING (the F-LEFT (x)1_F; DEVIATES from .tex:2862's
 * (Delta(U)(x)1_F) side — see paper/FINDINGS.md §C13(a) for why; full derivation in the
 * aic_factorize_upsilon3.c apply-map docstring). aic_ucp_kraus_to_stinespring packs
 * V F-MAJOR (V[a*N+i,j]=K_a[i,j]), so Phi(X) = V^dag (1_F (x) X) V and EVERY (x)1_F
 * factor (here in L_j: 1_F (x) Delta(U_{js}^dag)) MUST be F-LEFT to match V. eta=0
 * /r=1 is BLIND (F 1-dim); the r>1 D5 pin distinguishes (MEASURED: F-LEFT 4.4e-2 ~
 * O(eta) vs F-RIGHT 7.6e-1 ~ O(1), mixconj(4,2,0.03) r=6).
 *
 * D4 — xi_j = top RIGHT singular vector of C_j, so ||C_j xi_j|| = sigma_max(C_j).
 * (The LEFT vector gives the wrong norm for a NON-NORMAL C_j; for the current
 * fixtures C_j is Hermitian PSD — R_j = WWdag central-averaged then partial-traced
 * is Hermitian — so left == right and the choice is currently UNOBSERVABLE; it is
 * still the paper-correct choice and matters if a non-normal C_j arises. FINDING.)
 * Via aic_latd_svd: the right singular vectors are the conjugated ROWS of Vt, so
 * v_0[i] = conj(Vt[0,i]).
 *
 * D6 — the unitalization basin: ASSERT ||Upsilon'(1_H) - 1_B||_op < 1 (Rule 4,
 * the aic_funcalc_xpow(.,-1/2,1) basin) BEFORE the inverse-sqrt; B(L_j)-block
 * structure means Upsilon'(1_H) is block-diagonal, so its inverse-sqrt is the
 * block-diagonal inverse-sqrt (aic_funcalc_xpow on the whole n_B x n_B rep).
 */
#include <assert.h>
#include <complex.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_dhom.h"
#include "aic/aic_factorize.h"
#include "aic/aic_factorize_internal.h"
#include "aic/aic_funcalc.h"
#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

/* xi_j = top RIGHT singular vector of C_j (D4): aic_latd_svd on the e_j x e_j
 * midpoint of C_j; v_0[i] = conj(Vt[0,i]). Returns sigma_max(C_j) (the gauge-
 * invariant lem_RC(ii) quantity ||C_j xi_j||). `xi` is caller-init'd e_j x 1. */
static double xi_top_right(acb_mat_t xi, const acb_mat_t Cj, slong e_j)
{
    double _Complex *A = flint_malloc(sizeof(double _Complex) * e_j * e_j);
    double _Complex *Vt = flint_malloc(sizeof(double _Complex) * e_j * e_j);
    double *sv = flint_malloc(sizeof(double) * e_j);
    for (slong i = 0; i < e_j; i++)
        for (slong k = 0; k < e_j; k++) {
            const acb_ptr z = acb_mat_entry(Cj, i, k);
            A[i * e_j + k] = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR)
                           + arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR) * I;
        }
    aic_latd_svd(sv, NULL, Vt, A, e_j, e_j);   /* descending sv; rows of Vt = v_k^dag */
    for (slong i = 0; i < e_j; i++)
        acb_set_d_d(acb_mat_entry(xi, i, 0),
                    creal(conj(Vt[0 * e_j + i])), cimag(conj(Vt[0 * e_j + i])));
    double smax = sv[0];
    flint_free(sv); flint_free(Vt); flint_free(A);
    return smax;
}

/* The block-j-embedded Pauli U_{js} = S_ab as an element of B (S_ab in block j,
 * ZERO elsewhere — D3/B.1), so Delta(U_{js}) = the single j-th Choi_Delta term.
 * `out` caller-init'd n_B x n_B. */
static void embed_block_pauli(acb_mat_t out, const aic_dhom_B *B, slong j,
                              slong a, slong b, slong prec)
{
    slong dj = B->d[j], r_off = 0;
    for (slong l = 0; l < j; l++) r_off += B->d[l];
    acb_mat_t S;
    acb_mat_init(S, dj, dj);
    aic_dhom_pauli(S, dj, a, b, prec);
    acb_mat_zero(out);
    for (slong i = 0; i < dj; i++)
        for (slong k = 0; k < dj; k++)
            acb_set(acb_mat_entry(out, r_off + i, r_off + k), acb_mat_entry(S, i, k));
    acb_mat_clear(S);
}

void aic_factorize_int_build_Lj(acb_mat_t Lj, const aic_factorize *F, slong j,
                                const acb_mat_t Wj, const acb_mat_t xi, slong e_j,
                                int f_left, slong prec)
{
    const aic_dhom_B *B = F->v->B;
    slong N = F->N, r = F->r, dj = B->d[j], dim = dj * e_j;
    acb_mat_t IF, Us, embUsd, dUsd, dUsd_F, Wd, Uxi, WdUxi, VWdUxi, term;
    acb_mat_init(IF, r, r);
    acb_mat_one(IF);
    acb_mat_init(Us, dj, dj);
    acb_mat_init(embUsd, B->n_B, B->n_B);
    acb_mat_init(dUsd, N, N);                 /* Delta(U_{js}^dag) */
    acb_mat_init(dUsd_F, r * N, r * N);        /* (x)1_F factor */
    acb_mat_init(Wd, N, dim);                  /* W_j^dag */
    acb_mat_init(Uxi, dim, dj);                /* U_{js} (x) xi_j */
    acb_mat_init(WdUxi, N, dj);
    acb_mat_init(VWdUxi, r * N, dj);
    acb_mat_init(term, r * N, dj);
    acb_mat_conjugate_transpose(Wd, Wj);

    arb_t ps;
    arb_init(ps);
    arb_set_si(ps, 1);
    arb_div_si(ps, ps, dj * dj, prec);         /* p_{js} = 1/d_j^2 */

    acb_mat_zero(Lj);
    for (slong sa = 0; sa < dj; sa++)
        for (slong sb = 0; sb < dj; sb++) {
            aic_dhom_pauli(Us, dj, sa, sb, prec);          /* U_{js} on L_j */
            /* Delta(U_{js}^dag): embed S_ab in block j, then in-block dag (the dag
             * of a block-diagonal matrix is block-wise, so this is S_ab^dag in
             * block j, zero elsewhere = U_{js}^dag as an element of B). */
            embed_block_pauli(embUsd, B, j, sa, sb, prec); /* embeds S_ab */
            acb_mat_conjugate_transpose(embUsd, embUsd);   /* -> S_ab^dag in block j */
            aic_factorize_delta(dUsd, F, embUsd, prec);    /* Delta(U_{js}^dag) (UCP) */
            if (f_left)
                aic_mat_kronecker(dUsd_F, IF, dUsd, prec); /* 1_F (x) Delta(U^dag) */
            else
                aic_mat_kronecker(dUsd_F, dUsd, IF, prec); /* WRONG: Delta(U^dag)(x)1_F */
            aic_mat_kronecker(Uxi, Us, xi, prec);          /* U_{js} (x) xi_j */
            acb_mat_mul(WdUxi, Wd, Uxi, prec);             /* W_j^dag (U (x) xi) */
            acb_mat_mul(VWdUxi, F->V, WdUxi, prec);        /* V W_j^dag (U (x) xi) */
            acb_mat_mul(term, dUsd_F, VWdUxi, prec);
            acb_mat_scalar_mul_arb(term, term, ps, prec);
            acb_mat_add(Lj, Lj, term, prec);
        }
    arb_clear(ps);
    acb_mat_clear(term); acb_mat_clear(VWdUxi); acb_mat_clear(WdUxi);
    acb_mat_clear(Uxi); acb_mat_clear(Wd); acb_mat_clear(dUsd_F);
    acb_mat_clear(dUsd); acb_mat_clear(embUsd);
    acb_mat_clear(Us); acb_mat_clear(IF);
}

/* The full F3 build: V, W_j, C_j, xi_j, L_j, Upsilon'(1_H)^{-1/2}. */
void aic_factorize_upsilon_build(aic_factorize *F, double tol_sigma, slong prec)
{
    assert(F != NULL);
    assert(F->delta_ready &&
           "upsilon_build: aic_factorize_delta_build (F2) must run first");
    /* Idempotent: free a stale F3 cache (same teardown as aic_factorize_clear). */
    if (F->upsilon_ready) {
        slong mm = F->v->B->num_blocks;
        for (slong j = 0; j < mm; j++) {
            acb_mat_clear(F->W[j]); acb_mat_clear(F->C[j]);
            acb_mat_clear(F->xi[j]); acb_mat_clear(F->L[j]);
        }
        flint_free(F->W); flint_free(F->C); flint_free(F->xi);
        flint_free(F->L); flint_free(F->e);
        acb_mat_clear(F->V); acb_mat_clear(F->upsI_invsqrt);
        F->upsilon_ready = 0;
    }

    const aic_dhom_B *B = F->v->B;
    slong m = B->num_blocks, N = F->N;

    /* Phi's Stinespring V: H -> H (x) F, F-major (aic_ucp.h). */
    F->r = F->phi->r;
    acb_mat_init(F->V, N * F->r, N);
    aic_ucp_kraus_to_stinespring(F->V, F->phi, prec);

    F->e  = flint_malloc(sizeof(slong) * m);
    F->W  = flint_malloc(sizeof(acb_mat_t) * m);
    F->C  = flint_malloc(sizeof(acb_mat_t) * m);
    F->xi = flint_malloc(sizeof(acb_mat_t) * m);
    F->L  = flint_malloc(sizeof(acb_mat_t) * m);

    /* The genuine-bug guard (Rule 4; FINDINGS §C14). Accumulate the PSD-cone defect
     * mass clipped by each W_j extraction. The almost-idempotent UCP Delta is CP
     * only to O(eta^2), so this should be O(eta^2)-small (measured worst per-block
     * minEig -2.5e-6 at eta=2.3e-2, summed over blocks ~ a few 1e-6). If it exceeds
     * AIC_FACTORIZE_CONE_MASS_CEIL = 1e-2 the construction is FAR more broken than
     * the O(eta)-CP theory predicts (a real bug / non-CP input slipped past the
     * per-eigenvalue 1e-3 abort floor by spreading mass over many eigenvalues) —
     * fail loud, do NOT swallow. 1e-2 is ~4000x above the measured worst and below
     * an O(1) cone violation; it is also above the largest eta in scope so an
     * O(eta) (not O(eta^2)) cone defect would trip it. */
#define AIC_FACTORIZE_CONE_MASS_CEIL 1e-2
    double total_clipped = 0.0;

    for (slong j = 0; j < m; j++) {
        slong dj = B->d[j], e_j;
        double clipped_j = 0.0;
        aic_factorize_upsilon_Wj(F->W[j], &e_j, &clipped_j, F, j, prec); /* W_j (inits) */
        total_clipped += clipped_j;
        F->e[j] = e_j;
        acb_mat_t Rj;
        acb_mat_init(Rj, dj * e_j, dj * e_j);
        aic_factorize_upsilon_Rj(Rj, F, j, F->W[j], e_j, prec);
        acb_mat_init(F->C[j], e_j, e_j);
        aic_factorize_upsilon_Cj(F->C[j], Rj, dj, e_j, prec);
        acb_mat_clear(Rj);

        /* xi_j (D4) + lem_RC(ii) sigma_max(C_j) >= 1 - tol_sigma (Rule 4). */
        acb_mat_init(F->xi[j], e_j, 1);
        double smax = xi_top_right(F->xi[j], F->C[j], e_j);
        assert(smax >= 1.0 - tol_sigma &&
               "upsilon_build: sigma_max(C_j) < 1 - tol_sigma (lem_RC(ii) violated)");
        assert(smax <= 1.0 + tol_sigma &&
               "upsilon_build: sigma_max(C_j) > 1 + tol_sigma (||C_j|| <= 1 violated)");

        acb_mat_init(F->L[j], N * F->r, dj);
        aic_factorize_int_build_Lj(F->L[j], F, j, F->W[j], F->xi[j], e_j,
                                   /*f_left=*/1, prec);
    }

    /* The genuine-bug guard: total clipped PSD-cone-defect mass must be O(eta^2)-
     * small. PRINT it (visibility, CLAUDE.md Rule 12 concrete numbers), then
     * ASSERT the ceiling (Rule 4 — fail loud if the construction is more broken
     * than the O(eta)-CP theory predicts; FINDINGS §C14). */
    flint_printf("aic_factorize_upsilon_build: total PSD-cone-defect mass clipped "
                 "across %wd block(s) = %.6e (ceiling %.1e)\n",
                 m, total_clipped, AIC_FACTORIZE_CONE_MASS_CEIL);
    assert(total_clipped <= AIC_FACTORIZE_CONE_MASS_CEIL &&
           "upsilon_build: clipped PSD-cone-defect mass exceeds ceiling — the UCP "
           "Delta is more non-CP than O(eta^2) (FINDINGS §C14; Rule 4)");
#undef AIC_FACTORIZE_CONE_MASS_CEIL

    F->upsilon_ready = 1;

    /* Upsilon'(1_H) and the inverse-sqrt basin assert (D6). */
    acb_mat_t IH, UI, UImI;
    acb_mat_init(IH, N, N);
    acb_mat_init(UI, F->n_B, F->n_B);
    acb_mat_init(UImI, F->n_B, F->n_B);
    acb_mat_one(IH);
    aic_factorize_upsilon_prime(UI, F, IH, prec);     /* Upsilon'(1_H) */
    {
        aic_dhom_B_unit(UImI, B);                      /* I_B */
        acb_mat_sub(UImI, UI, UImI, prec);
        arb_t e, one;
        arb_init(e);
        arb_init(one);
        aic_mat_opnorm(e, UImI, prec);
        arb_one(one);
        assert(arb_lt(e, one) &&
               "upsilon_build: ||Upsilon'(1_H) - 1_B||_op >= 1 (inverse-sqrt out "
               "of basin; D6)");
        arb_clear(one);
        arb_clear(e);
    }
    acb_mat_init(F->upsI_invsqrt, F->n_B, F->n_B);
    aic_funcalc_xpow(F->upsI_invsqrt, UI, -0.5, 1.0, prec);
    acb_mat_clear(UImI); acb_mat_clear(UI); acb_mat_clear(IH);
}
