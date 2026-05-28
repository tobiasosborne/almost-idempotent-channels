/* aic_idemp_build.c — steps 5-7 of th_idemp_structure (Lambda, Gamma, w)
 * (bead aic-wuh; approximate_algebras.tex:2061-2084). Split out of
 * aic_idemp_maps.c to keep each file <= 200 LOC (Rule 10).
 *
 *   Lambda(Y) = Phi(J_M Y J_M^dag) : B(M) -> B(H)   (.tex:2061), CP + unital.
 *   Gamma = Delta^dag Lambda : B(M) -> A             (.tex:2065; Delta^dag Delta=1_A).
 *   w     = C_M Delta : A -> B(M),  col k = vec(J_M^dag B_k J_M)  (.tex:2084).
 *
 * Maps are stored row-major-vec (vec(X)[i*n+j]=X[i,j]); a map B(C^p)->B(C^q) is a
 * (q*q)x(p*p) matrix (see include/aic_idemp.h). Maps are applied via direct
 * per-column reshapes / matmuls in this and the verify files; no generic
 * map-matrix-times-matrix helper is exported (Gamma's codomain is abstract
 * A-coordinates, dA x 1, which does not fit a (q*q)x(p*p) signature).
 */
#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_idemp.h"
#include "aic_idemp_internal.h"
#include "aic_ucp.h"

/* reshape column k of Delta (length n^2, row-major vec) into the n x n matrix B. */
static void delta_col_to_mat(acb_mat_t B, const acb_mat_t Delta, slong k, slong n)
{
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(B, a, b), acb_mat_entry(Delta, a * n + b, k));
}

void aic_idemp_build_maps(aic_idemp_decomp *out, const aic_ucp_kraus *phi,
                          slong prec)
{
    slong n = out->n, dM = out->dim_M, dA = out->dim_A, dM2 = dM * dM;

    /* --- step 5: Lambda (n^2 x dim_M^2), Lambda(Y) = Phi(J_M Y J_M^dag) --- */
    acb_mat_init(out->Lambda, n * n, dM2);
    {
        acb_mat_t JMd, F, JF, JFJ, LamF;
        acb_mat_init(JMd, dM, n);
        acb_mat_init(F, dM, dM);
        acb_mat_init(JF, n, dM);
        acb_mat_init(JFJ, n, n);
        acb_mat_init(LamF, n, n);
        acb_mat_conjugate_transpose(JMd, out->J_M);
        for (slong k = 0; k < dM; k++)
            for (slong l = 0; l < dM; l++) {
                acb_mat_zero(F);
                acb_set_si(acb_mat_entry(F, k, l), 1);     /* F_kl on B(M) */
                acb_mat_mul(JF, out->J_M, F, prec);        /* J_M F */
                acb_mat_mul(JFJ, JF, JMd, prec);           /* J_M F J_M^dag */
                aic_ucp_apply(LamF, phi, JFJ, prec);       /* Phi(...) */
                slong col = k * dM + l;
                for (slong a = 0; a < n; a++)
                    for (slong b = 0; b < n; b++)
                        acb_set(acb_mat_entry(out->Lambda, a * n + b, col),
                                acb_mat_entry(LamF, a, b));
            }
        acb_mat_clear(LamF);
        acb_mat_clear(JFJ);
        acb_mat_clear(JF);
        acb_mat_clear(F);
        acb_mat_clear(JMd);
    }

    /* --- step 6: Gamma = Delta^dag Lambda (dim_A x dim_M^2) --- */
    acb_mat_init(out->Gamma, dA, dM2);
    {
        acb_mat_t Dd;
        acb_mat_init(Dd, dA, n * n);
        acb_mat_conjugate_transpose(Dd, out->Delta);
        acb_mat_mul(out->Gamma, Dd, out->Lambda, prec);
        acb_mat_clear(Dd);
    }

    /* --- step 7: w = C_M Delta (dim_M^2 x dim_A), col k = vec(J_M^dag B_k J_M) --- */
    acb_mat_init(out->w, dM2, dA);
    {
        acb_mat_t JMd, B, JdB, JdBJ;
        acb_mat_init(JMd, dM, n);
        acb_mat_init(B, n, n);
        acb_mat_init(JdB, dM, n);
        acb_mat_init(JdBJ, dM, dM);
        acb_mat_conjugate_transpose(JMd, out->J_M);
        for (slong k = 0; k < dA; k++) {
            delta_col_to_mat(B, out->Delta, k, n);
            acb_mat_mul(JdB, JMd, B, prec);          /* J_M^dag B */
            acb_mat_mul(JdBJ, JdB, out->J_M, prec);  /* J_M^dag B J_M */
            for (slong p = 0; p < dM; p++)
                for (slong q = 0; q < dM; q++)
                    acb_set(acb_mat_entry(out->w, p * dM + q, k),
                            acb_mat_entry(JdBJ, p, q));
        }
        acb_mat_clear(JdBJ);
        acb_mat_clear(JdB);
        acb_mat_clear(B);
        acb_mat_clear(JMd);
    }
}
