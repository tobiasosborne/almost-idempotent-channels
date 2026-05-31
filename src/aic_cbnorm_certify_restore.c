/* aic_cbnorm_certify_restore.c — the two feasibility-restoration recipes for the
 * TIGHT cb-norm certifier (bead aic-m24, increment 3b). The public dispatcher is
 * in src/aic_cbnorm_certify.c; this file holds the bulk of the math (split for
 * the <=200 LOC limit, Rule 10). See include/aic_cbnorm.h for the API and
 * docs/cbnorm_tight_certifier.md for THE CONTRACT.
 *
 * THE TWO RECIPES (Law 1; design doc LOWER/UPPER; .tex:347-354 cb-norm SDP).
 *
 * LOWER (Watrous MAX primal, Convention-A trace-n; eta = (2/n)*optval):
 *   max Re tr(J^dag X) s.t. block=[[P,X],[X^dag,Q]]>=0, P+Q=I_{n^2}.
 *   Restore by a convex combination toward the Slater center (0, I/2, I/2):
 *     E = P+Q-I; P'=P-E/2, Q'=Q-E/2 (so P'+Q'=I exactly); delta = rigorous
 *     UPPER bound on max(0,-lambda_min(block')); t = 1/(1+2*delta); the point
 *     (tX, tP'+(1-t)I/2, tQ'+(1-t)I/2) is FEASIBLE; lo = t*(2/n)*Re tr(J^dag X).
 *
 * UPPER (QETLAB/Watrous MIN dual, trace-n; eta = optval/2):
 *   min (1/2)(||Tr_2 Y0||_inf + ||Tr_2 Y1||_inf) s.t.
 *       block_D=[[Y0,-J],[-J^dag,Y1]]>=0, Y0,Y1>=0.
 *   Restore by a shift: eps = max(0,-lambda_min(block_D)); (Y0+eps I, Y1+eps I)
 *   is dual-feasible; Tr_2 = partial_trace_right (RIGHT/input factor, PINNED 3a);
 *   Tr_2(Y0+eps I_{n^2}) = Tr_2(Y0)+eps*n*I_n, so
 *     hi = (1/2)(lambda_max(Tr_2 Y0)+lambda_max(Tr_2 Y1)) + eps*n.
 *
 * Every PSD test is the ONE-SIDED herm_max_eig(-M) global enclosure (degeneracy-
 * robust; no full eig, dodging aic-w4o.1). arb/acb only.
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic_cbnorm_internal.h"
#include "aic/aic_mat.h"

/* Assemble out = [[A,B],[C,D]], each block m x m, into a 2m x 2m matrix. */
static void block2(acb_mat_t out, const acb_mat_t A, const acb_mat_t B,
                   const acb_mat_t C, const acb_mat_t D, slong m)
{
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++) {
            acb_set(acb_mat_entry(out, i, j), acb_mat_entry(A, i, j));
            acb_set(acb_mat_entry(out, i, m + j), acb_mat_entry(B, i, j));
            acb_set(acb_mat_entry(out, m + i, j), acb_mat_entry(C, i, j));
            acb_set(acb_mat_entry(out, m + i, m + j), acb_mat_entry(D, i, j));
        }
}

/* Rigorous UPPER bound on max(0, -lambda_min(M)) for Hermitian M: take the upper
 * endpoint of herm_max_eig(-M) (= -lambda_min(M)) and clamp at 0. */
static void psd_defect(arb_t out, const acb_mat_t M, slong prec)
{
    slong d = acb_mat_nrows(M);
    acb_mat_t negM;
    acb_mat_init(negM, d, d);
    acb_mat_neg(negM, M);

    arb_t e;
    arb_init(e);
    aic_mat_herm_max_eig(e, negM, prec);   /* ball enclosing -lambda_min(M) */

    arf_t ub;
    arf_init(ub);
    arb_get_ubound_arf(ub, e, prec);       /* rigorous upper endpoint */
    arf_set(arb_midref(out), ub);
    mag_zero(arb_radref(out));
    if (arf_sgn(ub) < 0) arb_zero(out);    /* clamp: defect >= 0 */

    arf_clear(ub);
    arb_clear(e);
    acb_mat_clear(negM);
}

void aic_cbnorm_int_lower(arb_t lo, arb_t delta, const acb_mat_t J,
                          const acb_mat_t X, const acb_mat_t P,
                          const acb_mat_t Q, slong n, slong prec)
{
    slong N = n * n;

    /* 1. symmetrize the equality: E = P+Q-I; P' = P-E/2, Q' = Q-E/2. */
    acb_mat_t Pp, Qp, E, I_N;
    acb_mat_init(Pp, N, N); acb_mat_init(Qp, N, N);
    acb_mat_init(E, N, N); acb_mat_init(I_N, N, N);
    acb_mat_one(I_N);
    acb_mat_add(E, P, Q, prec);
    acb_mat_sub(E, E, I_N, prec);                 /* E = P+Q-I */
    acb_mat_scalar_mul_2exp_si(E, E, -1);         /* E/2 */
    acb_mat_sub(Pp, P, E, prec);                  /* P' = P - E/2 */
    acb_mat_sub(Qp, Q, E, prec);                  /* Q' = Q - E/2 */

    /* 2. block' = [[P',X],[X^dag,Q']]; delta = max(0,-lambda_min(block')). */
    acb_mat_t Xd, blk;
    acb_mat_init(Xd, N, N); acb_mat_init(blk, 2 * N, 2 * N);
    acb_mat_conjugate_transpose(Xd, X);
    block2(blk, Pp, X, Xd, Qp, N);
    psd_defect(delta, blk, prec);

    /* 3. t = 1/(1 + 2*delta). */
    arb_t t;
    arb_init(t);
    arb_mul_2exp_si(t, delta, 1);                 /* 2*delta */
    arb_add_ui(t, t, 1, prec);                    /* 1 + 2*delta */
    arb_inv(t, t, prec);

    /* 4. lo = t * (2/n) * Re tr(J^dag X). */
    acb_mat_t Jd, M;
    acb_mat_init(Jd, N, N); acb_mat_init(M, N, N);
    acb_mat_conjugate_transpose(Jd, J);
    acb_mat_mul(M, Jd, X, prec);
    acb_t tr;
    acb_init(tr);
    acb_mat_trace(tr, M, prec);
    arb_t re;
    arb_init(re);
    acb_get_real(re, tr);                         /* Re tr(J^dag X) */
    arb_mul_2exp_si(re, re, 1);                   /* 2 * Re tr */
    arb_div_ui(re, re, (ulong) n, prec);          /* (2/n) Re tr */
    arb_mul(lo, t, re, prec);

    arb_clear(re); acb_clear(tr);
    acb_mat_clear(M); acb_mat_clear(Jd);
    arb_clear(t);
    acb_mat_clear(blk); acb_mat_clear(Xd);
    acb_mat_clear(I_N); acb_mat_clear(E);
    acb_mat_clear(Qp); acb_mat_clear(Pp);
}

void aic_cbnorm_int_upper(arb_t hi, arb_t eps, const acb_mat_t J,
                          const acb_mat_t Y0, const acb_mat_t Y1,
                          slong n, slong prec)
{
    slong N = n * n;

    /* block_D = [[Y0,-J],[-J^dag,Y1]]; eps = max(0,-lambda_min(block_D)). */
    acb_mat_t negJ, negJd, blk;
    acb_mat_init(negJ, N, N); acb_mat_init(negJd, N, N);
    acb_mat_init(blk, 2 * N, 2 * N);
    acb_mat_neg(negJ, J);
    acb_mat_conjugate_transpose(negJd, negJ);
    block2(blk, Y0, negJ, negJd, Y1, N);
    psd_defect(eps, blk, prec);

    /* Tr_2 = partial_trace_right (RIGHT/input factor, keep LEFT; PINNED 3a).
     * lambda_max of the n x n Hermitian PSD marginal = ||.||_inf. */
    acb_mat_t T0, T1;
    acb_mat_init(T0, n, n); acb_mat_init(T1, n, n);
    aic_mat_partial_trace_right(T0, Y0, n, n, prec);
    aic_mat_partial_trace_right(T1, Y1, n, n, prec);
    arb_t s0, s1;
    arb_init(s0); arb_init(s1);
    aic_mat_herm_max_eig(s0, T0, prec);
    aic_mat_herm_max_eig(s1, T1, prec);

    /* hi = (1/2)(s0 + s1) + eps*n. */
    arb_add(hi, s0, s1, prec);
    arb_mul_2exp_si(hi, hi, -1);                   /* (s0+s1)/2 */
    arb_t shift;
    arb_init(shift);
    arb_mul_ui(shift, eps, (ulong) n, prec);       /* eps*n */
    arb_add(hi, hi, shift, prec);

    arb_clear(shift);
    arb_clear(s1); arb_clear(s0);
    acb_mat_clear(T1); acb_mat_clear(T0);
    acb_mat_clear(blk); acb_mat_clear(negJd); acb_mat_clear(negJ);
}
