/* aic_idemp_verify2.c — the w-homomorphism and block-diagonal defects of
 * th_idemp_structure (bead aic-wuh; approximate_algebras.tex:1916 lem_idemp,
 * :2090). Split out of aic_idemp_verify.c to keep each file <= 200 LOC (Rule 10);
 * see that file's header for the coordinate conventions.
 *
 *  - w is a *-homomorphism (lem_idemp :1916, second part): w carries the
 *    Choi-Effros product X star Y = Phi(XY) on A to the operator product on B(M).
 *    We certify max_{j,k} ||w(B_j star B_k) - w(B_j) w(B_k)||_op = 0.
 *  - operators in Im Delta are block-diagonal w.r.t. H = M (+) M^perp (lem_idemp
 *    first part, :2090): max_k (||(1-Pi_M)B_k Pi_M|| + ||Pi_M B_k(1-Pi_M)||) = 0.
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_idemp.h"
#include "aic_mat.h"
#include "aic_ucp.h"

/* reshape column k of Delta into B (n x n). */
static void col_to_mat(acb_mat_t B, const acb_mat_t Delta, slong k, slong n)
{
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(B, a, b), acb_mat_entry(Delta, a * n + b, k));
}

void aic_idemp_defect_w_hom(arb_t out, const aic_idemp_decomp *d,
                            const aic_ucp_kraus *phi, slong prec)
{
    slong n = d->n, dM = d->dim_M, dA = d->dim_A;
    acb_mat_t Bj, Bk, BjBk, PhiBjBk, vP, coords, Dd;
    acb_mat_t wBj, wBk, wprod, wstar, wcol, diff;
    arb_t dd;
    acb_mat_init(Bj, n, n);
    acb_mat_init(Bk, n, n);
    acb_mat_init(BjBk, n, n);
    acb_mat_init(PhiBjBk, n, n);
    acb_mat_init(vP, n * n, 1);
    acb_mat_init(coords, dA, 1);
    acb_mat_init(Dd, dA, n * n);
    acb_mat_init(wBj, dM, dM);
    acb_mat_init(wBk, dM, dM);
    acb_mat_init(wprod, dM, dM);
    acb_mat_init(wstar, dM, dM);
    acb_mat_init(wcol, dM * dM, 1);
    acb_mat_init(diff, dM, dM);
    arb_init(dd);
    acb_mat_conjugate_transpose(Dd, d->Delta);
    arb_zero(out);

    for (slong j = 0; j < dA; j++) {
        col_to_mat(Bj, d->Delta, j, n);
        for (slong p = 0; p < dM; p++)             /* w(B_j) = reshape col j of w */
            for (slong q = 0; q < dM; q++)
                acb_set(acb_mat_entry(wBj, p, q), acb_mat_entry(d->w, p * dM + q, j));
        for (slong k = 0; k < dA; k++) {
            col_to_mat(Bk, d->Delta, k, n);
            for (slong p = 0; p < dM; p++)
                for (slong q = 0; q < dM; q++)
                    acb_set(acb_mat_entry(wBk, p, q),
                            acb_mat_entry(d->w, p * dM + q, k));
            acb_mat_mul(wprod, wBj, wBk, prec);          /* w(B_j) w(B_k) */

            acb_mat_mul(BjBk, Bj, Bk, prec);             /* B_j B_k */
            aic_ucp_apply(PhiBjBk, phi, BjBk, prec);     /* star = Phi(B_j B_k) */
            for (slong a = 0; a < n; a++)
                for (slong b = 0; b < n; b++)
                    acb_set(acb_mat_entry(vP, a * n + b, 0),
                            acb_mat_entry(PhiBjBk, a, b));
            acb_mat_mul(coords, Dd, vP, prec);           /* A-coords = Delta^dag vec */
            acb_mat_mul(wcol, d->w, coords, prec);       /* w(star) = w_mat @ coords */
            for (slong p = 0; p < dM; p++)
                for (slong q = 0; q < dM; q++)
                    acb_set(acb_mat_entry(wstar, p, q),
                            acb_mat_entry(wcol, p * dM + q, 0));
            acb_mat_sub(diff, wstar, wprod, prec);
            aic_mat_opnorm(dd, diff, prec);
            if (arb_gt(dd, out)) arb_set(out, dd);
        }
    }
    arb_clear(dd);
    acb_mat_clear(diff);
    acb_mat_clear(wcol);
    acb_mat_clear(wstar);
    acb_mat_clear(wprod);
    acb_mat_clear(wBk);
    acb_mat_clear(wBj);
    acb_mat_clear(Dd);
    acb_mat_clear(coords);
    acb_mat_clear(vP);
    acb_mat_clear(PhiBjBk);
    acb_mat_clear(BjBk);
    acb_mat_clear(Bk);
    acb_mat_clear(Bj);
}

void aic_idemp_defect_blockdiag(arb_t out, const aic_idemp_decomp *d, slong prec)
{
    slong n = d->n, dA = d->dim_A;
    acb_mat_t B, t1, off1, t2, off2, one, Qc;
    arb_t a1, a2, sum;
    acb_mat_init(B, n, n);
    acb_mat_init(one, n, n);
    acb_mat_init(Qc, n, n);          /* 1 - Pi_M */
    acb_mat_init(t1, n, n);
    acb_mat_init(off1, n, n);
    acb_mat_init(t2, n, n);
    acb_mat_init(off2, n, n);
    arb_init(a1);
    arb_init(a2);
    arb_init(sum);
    acb_mat_one(one);
    acb_mat_sub(Qc, one, d->Pi_M, prec);   /* 1 - Pi_M */
    arb_zero(out);
    for (slong k = 0; k < dA; k++) {
        col_to_mat(B, d->Delta, k, n);
        acb_mat_mul(t1, Qc, B, prec);              /* (1-Pi_M) B Pi_M */
        acb_mat_mul(off1, t1, d->Pi_M, prec);
        aic_mat_opnorm(a1, off1, prec);
        acb_mat_mul(t2, d->Pi_M, B, prec);         /* Pi_M B (1-Pi_M) */
        acb_mat_mul(off2, t2, Qc, prec);
        aic_mat_opnorm(a2, off2, prec);
        arb_add(sum, a1, a2, prec);
        if (arb_gt(sum, out)) arb_set(out, sum);
    }
    arb_clear(sum);
    arb_clear(a2);
    arb_clear(a1);
    acb_mat_clear(off2);
    acb_mat_clear(t2);
    acb_mat_clear(off1);
    acb_mat_clear(t1);
    acb_mat_clear(Qc);
    acb_mat_clear(one);
    acb_mat_clear(B);
}
