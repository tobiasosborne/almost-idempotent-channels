/* aic_ecstar_setup.c — double-snapshot + scratch lifecycle + defect-map term
 * thunks for the Cycle-2 HOPM search (bead aic-knm, module "ecstar"). Split from
 * aic_ecstar_certify.c (Rule 10): this is the glue between the public arb-facing
 * aic_ecstar (acb_mat basis + Kraus) and the double-path engine (aic_ehk). The
 * term thunks fill the active defect map when the active block is set to B_k
 * (the other blocks read from w->frz1, w->frz2); they are the only place the
 * three defects (assoc trilinear, submult bilinear, cstar sesquilinear) differ.
 */
#include <complex.h>

#include "aic_ecstar.h"
#include "aic_ecstar_setup.h"
#include "aic_latd.h"

static cx *mat_to_cx(const acb_mat_t M, slong n)
{
    cx *buf = flint_malloc((size_t) (n * n) * sizeof(cx));
    aic_latd_from_acb_mat(buf, M);
    return buf;
}

void aic_ehk_snapshot(aic_ehk *h, const aic_ecstar *A)
{
    slong n = A->n, d = A->dim_A, r = A->phi->r;
    h->n = n; h->d = d; h->r = r;
    h->K = flint_malloc((size_t) r * sizeof(cx *));
    h->B = flint_malloc((size_t) d * sizeof(cx *));
    for (slong a = 0; a < r; a++) h->K[a] = mat_to_cx(A->phi->K[a], n);
    for (slong k = 0; k < d; k++) h->B[k] = mat_to_cx(A->B[k], n);
}

void aic_ehk_snapshot_free(aic_ehk *h)
{
    for (slong a = 0; a < h->r; a++) flint_free(h->K[a]);
    for (slong k = 0; k < h->d; k++) flint_free(h->B[k]);
    flint_free(h->K);
    flint_free(h->B);
}

void aic_ehk_work_alloc(aic_ehk_work *w, slong n, slong d)
{
    size_t nn = (size_t) (n * n);
    cx **slots[] = {&w->g, &w->C, &w->pol, &w->U, &w->Vt,
                    &w->s1, &w->s2, &w->s3, &w->s4, &w->s5, &w->s6};
    for (size_t i = 0; i < sizeof(slots) / sizeof(slots[0]); i++)
        *slots[i] = flint_malloc(nn * sizeof(cx));
    w->vec = flint_malloc((size_t) n * sizeof(cx));
    w->cf = flint_malloc((size_t) d * sizeof(cx));
    w->sv = flint_malloc((size_t) n * sizeof(double));
    w->frz1 = w->frz2 = NULL;
}

void aic_ehk_work_free(aic_ehk_work *w)
{
    cx *slots[] = {w->g, w->C, w->pol, w->U, w->Vt,
                   w->s1, w->s2, w->s3, w->s4, w->s5, w->s6};
    for (size_t i = 0; i < sizeof(slots) / sizeof(slots[0]); i++)
        flint_free(slots[i]);
    flint_free(w->vec); flint_free(w->cf); flint_free(w->sv);
}

/* ---- defect-map term thunks (active slot <- B_k; others = frz1,frz2) ----- */

void aic_ehk_term_assoc_X(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k)
{ aic_ehk_assoc(out, h, h->B[k], w->frz1, w->frz2, w->s1,w->s2,w->s3,w->s4,w->s5); }
void aic_ehk_term_assoc_Y(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k)
{ aic_ehk_assoc(out, h, w->frz1, h->B[k], w->frz2, w->s1,w->s2,w->s3,w->s4,w->s5); }
void aic_ehk_term_assoc_Z(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k)
{ aic_ehk_assoc(out, h, w->frz1, w->frz2, h->B[k], w->s1,w->s2,w->s3,w->s4,w->s5); }
void aic_ehk_term_sub_X(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k)
{ aic_ehk_star(out, h, h->B[k], w->frz1, w->s1, w->s2, w->s3); }
void aic_ehk_term_sub_Y(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k)
{ aic_ehk_star(out, h, w->frz1, h->B[k], w->s1, w->s2, w->s3); }
void aic_ehk_term_cstar(cx *out, aic_ehk_work *w, const aic_ehk *h, slong k)
{   /* gradient of ||X^dag*X||: (B_k^dag * X) + (X^dag * B_k) */
    slong n = h->n;
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) w->s6[i*n+j] = conj(h->B[k][j*n+i]);
    aic_ehk_star(out, h, w->s6, w->frz1, w->s1, w->s2, w->s3);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) w->s6[i*n+j] = conj(w->frz1[j*n+i]);
    aic_ehk_star(w->s4, h, w->s6, h->B[k], w->s1, w->s2, w->s3);
    for (slong p = 0; p < n * n; p++) out[p] += w->s4[p];
}
