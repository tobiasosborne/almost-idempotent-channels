/* aic_projection_audit.c — the UNIT-AWARE gap audition (bead aic-66n), the
 * root-cause fix for the wrapper-collapse abort. Splits Steps 3-4 (ambient
 * projector + project-into-A) and the candidate enumeration out of
 * src/aic_projection.c (Rule 10, ~200 LOC/file). See include/aic_projection.h for
 * the full Route-A-ambient narrative; this file is the audition loop.
 *
 * THE BUG (aic-66n). aic_cstar_build recurses Stage-1 splits onto S_P WRAPPER
 * algebras whose unit is Ptilde_m = Co_P(P), an n x n rank-r ambient projector
 * (r << n), NOT the ambient 1_n (approximate_algebras.tex:1082: "the compressed
 * product turns the Banach space S_P ... into an O(delta+eps)-C* algebra with unit
 * Ptilde = Co_P(P)"). The chosen H = (B_k+B_k^dag)/2 is supported on range(Ptilde_m),
 * so in M_n its spectrum has r genuine eigenvalues + (n-r) near-ZERO complement
 * ones, and the LARGEST interior gap is the support-vs-complement gap. Building P
 * from THAT gap gives P = Phi_tilde(P_amb) ~= Ptilde_m = the wrapper UNIT (a TRIVIAL
 * split). The old gate tested the AMBIENT ||1_n - P|| ~ 1, VACUOUS on a wrapper
 * (FINDINGS §C11), so it missed the collapse.
 *
 * THE FIX. The nontriviality contract (tex:929: "a delta-projection P is nontrivial
 * if both P and I-P are nonvanishing") is measured against the ALGEBRA'S unit, not
 * 1_n. For a wrapper that is Ptilde_m. So we audit gaps largest-first and ACCEPT
 * the first P with ||P|| > c AND ||U_A - P|| > c, U_A = Phi_tilde(1_n) (= Ptilde_m
 * for a wrapper, ~1_n at top level; FINDINGS §C7 established the wrapper unit is
 * Ptilde, NOT 1_n). ROBUSTNESS: the audition ranges over the (H_k, gap) pairs of
 * ALL non-scalar H_k, so if the largest-ambient-gap H_k has no usable in-support
 * gap (its r genuine eigenvalues clustered) while a different H_k does, we still
 * find it — a single pre-chosen H would fail loud spuriously.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner_internal.h"   /* aic_corner_gamma_opnorm_ub (aic-qgs dodge) */
#include "aic/aic_ecstar.h"
#include "aic/aic_latd.h"
#include "aic_projection_internal.h"

/* All interior gaps of a Hermitian H_k (double-path zheev spectrum), appended to
 * cands[*nc..]. Returns the spread (lam_max - lam_min); appends nothing if H_k is
 * near-scalar (spread <= 1e-6). The gap SIZE is not gated (any positive interior
 * gap is a valid candidate; the floor 1e-9*max(1,spread) only excludes round-off
 * noise — same rationale as aic_projection_gap, FINDINGS-cited there). */
static double append_gaps(aic_projection_cand *cands, slong *nc, slong k,
                          const acb_mat_t H)
{
    slong n = acb_mat_nrows(H);
    double _Complex *Harr = malloc((size_t) n * n * sizeof(double _Complex));
    double *ev = malloc((size_t) n * sizeof(double));
    if (!Harr || !ev) { fprintf(stderr, "append_gaps: OOM\n"); abort(); }
    aic_latd_from_acb_mat(Harr, H);
    aic_latd_eig_hermitian(ev, NULL, Harr, n);   /* ascending */
    double spread = ev[n - 1] - ev[0];
    if (spread > 1e-6) {
        double g_min = 1e-9 * fmax(1.0, spread);
        for (slong i = 0; i < n - 1; i++) {
            double g = ev[i + 1] - ev[i];
            if (g < g_min) continue;
            aic_projection_cand *c = &cands[(*nc)++];
            c->k = k; c->gap = g; c->t = 0.5 * (ev[i] + ev[i + 1]);
            c->lmin = ev[0]; c->lmax = ev[n - 1]; c->m = n - 1 - i;
        }
    }
    free(ev); free(Harr);
    return spread;
}

/* qsort comparator: candidates by gap DESCENDING. */
static int cand_cmp(const void *p, const void *q)
{
    double a = ((const aic_projection_cand *) p)->gap;
    double b = ((const aic_projection_cand *) q)->gap;
    return (a < b) - (a > b);
}

void aic_projection_enum_cands(aic_projection_cand **cands_out, slong *n_cands,
                               const aic_ecstar *A, slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(d > 1 && "projection: lem_nontriv_projection needs 1 < dim A");
    /* at most (n-1) interior gaps per H_k, d of them */
    aic_projection_cand *cands = malloc((size_t) d * (size_t) (n - 1) *
                                        sizeof(aic_projection_cand));
    if (!cands) { fprintf(stderr, "enum_cands: OOM\n"); abort(); }
    slong nc = 0;
    double best_spread = -1.0;
    acb_mat_t Hk;
    acb_mat_init(Hk, n, n);
    for (slong k = 0; k < d; k++) {
        aic_projection_herm_part(Hk, A->B[k], prec);
        double sp = append_gaps(cands, &nc, k, Hk);
        if (sp > best_spread) best_spread = sp;
    }
    acb_mat_clear(Hk);
    if (nc == 0) {
        fprintf(stderr, "aic_projection: dim A = %ld > 1 but NO non-scalar "
                "Hermitian basis part H_k has a positive interior spectral gap "
                "(max spread %.3e). The algebra is DEGENERATE; "
                "lem_nontriv_projection's spectral route cannot proceed — escalate "
                "(bead aic-3qv stop condition).\n", (long) d, best_spread);
        free(cands);
        abort();
    }
    qsort(cands, (size_t) nc, sizeof(aic_projection_cand), cand_cmp);
    *cands_out = cands;
    *n_cands = nc;
}

void aic_projection_audit(acb_mat_t P, const aic_ecstar *A, const acb_mat_t U_A,
                          aic_projection_witness *wit, double c, slong prec)
{
    slong n = A->n;
    aic_projection_cand *cands = NULL;
    slong nc = 0;
    aic_projection_enum_cands(&cands, &nc, A, prec);

    acb_mat_t Hk, Pamb, ImP;
    acb_mat_init(Hk, n, n);
    acb_mat_init(Pamb, n, n);
    acb_mat_init(ImP, n, n);
    arb_t pn, un;
    arb_init(pn); arb_init(un);

    slong accepted = -1;
    /* Audition largest-gap-first. ACCEPT the first P=Phi_tilde(P_amb) that is
     * NONTRIVIAL VS THE ALGEBRA UNIT U_A (tex:929): ||P|| > c AND ||U_A - P|| > c.
     * At top level U_A ~ 1_n so this matches the classic test; on an S_P wrapper
     * U_A = Ptilde_m, so the support-vs-complement gap (which yields P ~ Ptilde_m =
     * U_A) is correctly REJECTED and the next (within-support) gap is auditioned. */
    for (slong i = 0; i < nc && accepted < 0; i++) {
        aic_projection_herm_part(Hk, A->B[cands[i].k], prec);
        aic_projection_ambient(Pamb, Hk, cands[i].t, cands[i].lmin, cands[i].lmax,
                               prec);
        aic_projection_into_A(P, A, Pamb, prec);
        aic_corner_gamma_opnorm_ub(pn, NULL, P, prec);          /* ||P||_op */
        acb_mat_sub(ImP, U_A, P, prec);
        aic_corner_gamma_opnorm_ub(un, NULL, ImP, prec);        /* ||U_A - P||_op */
        double pnd = arf_get_d(arb_midref(pn), ARF_RND_NEAR);
        double und = arf_get_d(arb_midref(un), ARF_RND_NEAR);
        if (pnd > c && und > c) {
            accepted = i;
            if (wit) {
                wit->k_chosen = cands[i].k;
                wit->lam_min = cands[i].lmin;
                wit->lam_max = cands[i].lmax;
                wit->t = cands[i].t;
                wit->gap = cands[i].gap;
                wit->m = cands[i].m;
            }
        }
    }

    arb_clear(un); arb_clear(pn);
    acb_mat_clear(ImP); acb_mat_clear(Pamb); acb_mat_clear(Hk);

    if (accepted < 0) {
        fprintf(stderr, "aic_projection: NO (H_k, gap) candidate (%ld auditioned) "
                "yields a nontrivial-vs-unit split (every P has ||P|| <= %.3g or "
                "||U_A - P|| <= %.3g, i.e. P ~ 0 or P ~ the algebra unit U_A). The "
                "algebra is degenerate / has no usable interior gap — escalate "
                "(bead aic-3qv stop condition).\n", (long) nc, c, c);
        free(cands);
        abort();
    }
    free(cands);
}
