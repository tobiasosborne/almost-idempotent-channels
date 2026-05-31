/* aic_ucp_core.c — UCP-map core on the certified arb path (bead aic-c7n):
 * Kraus lifecycle, the Heisenberg action, Kraus->Choi, the unital and CP
 * certificates, and the CP cb-norm. See include/aic_ucp.h for the data model
 * and the load-bearing Choi index convention.
 *
 * Paper result vs constructive route.
 *  - Choi representation: the paper states the isometry form
 *    Phi(X)=V^dag(X(x)1_F)V (.tex:1570) with V^dag V=1_H (.tex:1571); we work in
 *    the equivalent Kraus form Phi(X)=sum_a K_a^dag X K_a, the F-basis expansion
 *    of that isometry form (shard G:25). Kraus->Choi is the textbook contraction
 *    C_Phi[i*h+a, j*h+b] = sum_x conj(K_x[i,a]) K_x[j,b]; one matrix loop, no
 *    eigenproblem (constructive and exact up to rounding).
 *  - Unital (V^dag V=1_H, .tex:1571, i.e. sum_a K_a^dag K_a = 1_H): checked TWO
 *    ways — straight from Kraus, and via Tr_K(C_Phi)=1_H — and the two are
 *    cross-checked. That cross-check is transpose-invariant, so it validates the
 *    partial-trace DIRECTION only; the per-entry conjugation convention is
 *    validated by the Choi-convention oracle test (test_choi_convention_oracle).
 *  - CP (C_Phi >= 0 iff CP): the paper builds Phi as UCP by construction; here we
 *    CERTIFY positivity by enclosing the largest eigenvalue of -C_Phi with
 *    aic_mat_herm_max_eig (degeneracy-robust global enclosure), NOT the simple-eig
 *    path which the degenerate Choi spectrum would abort (bead aic-w4o.1).
 *  - cb-norm (.tex:1717-1719, ||Phi||_cb<=1; :359 =1 for UCP): closed form
 *    ||Phi||_cb=||Phi(1_K)||_op (Paulsen Prop 3.6). No SDP needed for CP maps.
 *
 * Fail-loud (Rule 4): dimensions asserted; aic_ucp_is_cp_choi aborts when the
 * eigenvalue ball straddles the boundary (unresolved gap), never guesses.
 * Precision: the unital/CP differences are of O(1) quantities, so prec=53 already
 * resolves them to ~1e-15; the tests cross-check double-vs-arb@53.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"

void aic_ucp_kraus_init(aic_ucp_kraus *phi, slong dim_K, slong dim_H, slong r)
{
    assert(phi != NULL && dim_K >= 1 && dim_H >= 1 && r >= 1 &&
           "aic_ucp_kraus_init: positive dims and at least one Kraus op");
    phi->dim_K = dim_K;
    phi->dim_H = dim_H;
    phi->r = r;
    phi->K = flint_malloc((size_t) r * sizeof(acb_mat_t));
    assert(phi->K != NULL && "aic_ucp_kraus_init: out of memory");
    for (slong a = 0; a < r; a++) {
        acb_mat_init(phi->K[a], dim_K, dim_H);
        acb_mat_zero(phi->K[a]);
    }
}

void aic_ucp_kraus_clear(aic_ucp_kraus *phi)
{
    if (phi == NULL || phi->K == NULL) return;
    for (slong a = 0; a < phi->r; a++) acb_mat_clear(phi->K[a]);
    flint_free(phi->K);
    phi->K = NULL;
}

/* .tex:1570 — Phi(X) = sum_a K_a^dag X K_a (Heisenberg/observable action). */
void aic_ucp_apply(acb_mat_t out, const aic_ucp_kraus *phi, const acb_mat_t X,
                   slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H;
    assert(acb_mat_nrows(X) == dk && acb_mat_ncols(X) == dk &&
           "aic_ucp_apply: X must be dim_K x dim_K");
    assert(acb_mat_nrows(out) == dh && acb_mat_ncols(out) == dh &&
           "aic_ucp_apply: out must be dim_H x dim_H");

    acb_mat_t Kd, T, term;
    acb_mat_init(Kd, dh, dk);   /* K_a^dag : K -> H */
    acb_mat_init(T, dh, dk);    /* K_a^dag X */
    acb_mat_init(term, dh, dh); /* K_a^dag X K_a */
    acb_mat_zero(out);
    for (slong a = 0; a < phi->r; a++) {
        acb_mat_conjugate_transpose(Kd, phi->K[a]);
        acb_mat_mul(T, Kd, X, prec);
        acb_mat_mul(term, T, phi->K[a], prec);
        acb_mat_add(out, out, term, prec);
    }
    acb_mat_clear(term);
    acb_mat_clear(T);
    acb_mat_clear(Kd);
}

/* Choi Convention A (header): C[i*h+a, j*h+b] = sum_x conj(K_x[i,a]) K_x[j,b]. */
void aic_ucp_kraus_to_choi(acb_mat_t C, const aic_ucp_kraus *phi, slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H;
    assert(acb_mat_nrows(C) == dk * dh && acb_mat_ncols(C) == dk * dh &&
           "aic_ucp_kraus_to_choi: C must be (dim_K*dim_H) square");

    acb_t acc, term, cj;
    acb_init(acc);
    acb_init(term);
    acb_init(cj);
    for (slong i = 0; i < dk; i++)
        for (slong a = 0; a < dh; a++)
            for (slong j = 0; j < dk; j++)
                for (slong b = 0; b < dh; b++) {
                    acb_zero(acc);
                    for (slong x = 0; x < phi->r; x++) {
                        acb_conj(cj, acb_mat_entry(phi->K[x], i, a));
                        acb_mul(term, cj, acb_mat_entry(phi->K[x], j, b), prec);
                        acb_add(acc, acc, term, prec);
                    }
                    acb_set(acb_mat_entry(C, i * dh + a, j * dh + b), acc);
                }
    acb_clear(cj);
    acb_clear(term);
    acb_clear(acc);
}

void aic_ucp_unital_defect_choi(arb_t out, const acb_mat_t C, slong dim_K,
                                slong dim_H, slong prec)
{
    assert(acb_mat_nrows(C) == dim_K * dim_H && "unital_defect_choi: shape");
    acb_mat_t tr, diff;
    acb_mat_init(tr, dim_H, dim_H);
    acb_mat_init(diff, dim_H, dim_H);
    /* Tr_K(C): K is the LEFT factor, so trace OUT the left factor (keeps H). */
    aic_mat_partial_trace_left(tr, C, dim_K, dim_H, prec);
    acb_mat_one(diff);                       /* 1_H */
    acb_mat_sub(diff, tr, diff, prec);       /* Tr_K(C) - 1_H */
    aic_mat_opnorm(out, diff, prec);
    acb_mat_clear(diff);
    acb_mat_clear(tr);
}

void aic_ucp_unital_defect_kraus(arb_t out, const aic_ucp_kraus *phi, slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H;
    acb_mat_t Kd, term, acc;
    acb_mat_init(Kd, dh, dk);
    acb_mat_init(term, dh, dh);
    acb_mat_init(acc, dh, dh);
    acb_mat_zero(acc);
    for (slong a = 0; a < phi->r; a++) {            /* sum_a K_a^dag K_a */
        acb_mat_conjugate_transpose(Kd, phi->K[a]);
        acb_mat_mul(term, Kd, phi->K[a], prec);
        acb_mat_add(acc, acc, term, prec);
    }
    acb_mat_t one;
    acb_mat_init(one, dh, dh);
    acb_mat_one(one);
    acb_mat_sub(acc, acc, one, prec);               /* - 1_H */
    aic_mat_opnorm(out, acc, prec);
    acb_mat_clear(one);
    acb_mat_clear(acc);
    acb_mat_clear(term);
    acb_mat_clear(Kd);
}

int aic_ucp_is_cp_choi(const acb_mat_t C, const arb_t tol, slong prec)
{
    acb_mat_t negC;
    arb_t lam;
    int verdict;
    slong n = acb_mat_nrows(C);
    assert(acb_mat_ncols(C) == n && "aic_ucp_is_cp_choi: C must be square");

    /* CP via Choi presumes C Hermitian; lambda_min is meaningless otherwise.
     * Fail loud with a CP-specific message (Rule 4) BEFORE handing C to
     * aic_mat_herm_max_eig, rather than relying on that routine's internal
     * assert. Require ||C - C^dag||_op certainly small: a Hermitian Choi built
     * from Kraus is symmetric to ~2^-prec rounding, so 2^(-prec/2) (~1e-8 at
     * prec=53) is generous for rounding yet catches genuine non-Hermitian input. */
    {
        acb_mat_t Cd, asym;
        arb_t herm_def, herm_tol;
        acb_mat_init(Cd, n, n);
        acb_mat_init(asym, n, n);
        arb_init(herm_def);
        arb_init(herm_tol);
        acb_mat_conjugate_transpose(Cd, C);
        acb_mat_sub(asym, C, Cd, prec);            /* C - C^dag */
        aic_mat_opnorm(herm_def, asym, prec);
        arb_set_si(herm_tol, 1);
        arb_mul_2exp_si(herm_tol, herm_tol, -prec / 2); /* 2^(-prec/2) */
        if (!arb_le(herm_def, herm_tol)) {
            flint_printf("aic_ucp_is_cp_choi: C is not Hermitian "
                         "(||C - C^dag||_op exceeds 2^(-prec/2) at prec=%wd); a CP "
                         "verdict requires a Hermitian Choi (Rule 4: fail loud)\n",
                         prec);
            flint_abort();
        }
        arb_clear(herm_tol);
        arb_clear(herm_def);
        acb_mat_clear(asym);
        acb_mat_clear(Cd);
    }

    acb_mat_init(negC, n, n);
    arb_init(lam);
    acb_mat_neg(negC, C);
    aic_mat_herm_max_eig(lam, negC, prec); /* lambda_max(-C) = -lambda_min(C) */

    if (arb_le(lam, tol)) {            /* whole ball <= tol  => C PSD certified */
        verdict = 1;
    } else if (arb_gt(lam, tol)) {     /* whole ball >  tol  => C not PSD       */
        verdict = 0;
    } else {                           /* ball straddles tol => unresolved gap  */
        flint_printf("aic_ucp_is_cp_choi: -lambda_min(C) ball straddles tol at "
                     "prec=%wd; spectral gap unresolved (Rule 4: fail loud)\n",
                     prec);
        flint_abort();
        verdict = 0; /* unreachable */
    }
    arb_clear(lam);
    acb_mat_clear(negC);
    return verdict;
}

void aic_ucp_cbnorm_cp(arb_t out, const aic_ucp_kraus *phi, slong prec)
{
    /* .tex:1717-1719 / :359 — ||Phi||_cb = ||Phi(1_K)||_op (=1 for UCP). */
    slong dk = phi->dim_K, dh = phi->dim_H;
    acb_mat_t one_K, img;
    acb_mat_init(one_K, dk, dk);
    acb_mat_init(img, dh, dh);
    acb_mat_one(one_K);
    aic_ucp_apply(img, phi, one_K, prec);
    aic_mat_opnorm(out, img, prec);
    acb_mat_clear(img);
    acb_mat_clear(one_K);
}
