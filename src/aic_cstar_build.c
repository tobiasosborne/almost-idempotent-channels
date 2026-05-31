/* aic_cstar_build.c — the MASTER LOOP, constructive proof of th_main (bead
 * aic-097, module "cstar_build", Increment 5). See include/aic_cstar.h
 * (Increment 5 banner) for the full why, the three-stage route, the asserts, and
 * the ownership/lifetime contract; design contract docs/research/
 * cstar_masterloop_spec.md (the signature-verified spec). This file carries the
 * driver, Stage 1 (the greedy commutative skeleton), the equivalence-class
 * partition, and the c_0 / eps' bookkeeping; Stage 2 (per-class induction) and
 * Stage 3 (class merge) live in src/aic_cstar_stages.c (Rule 10 split).
 *
 * approximate_algebras.tex:1417-1426 — Stage 1 (greedy + contradiction, verbatim):
 *   "Let v_comm : B_comm -> A be a c_0*eps-inclusion of a commutative C* algebra of
 *    maximum dimensionality into A. ... We now show that these projections are
 *    one-dimensional. ... Without loss of generality, ... dim S_{P_m} > 1. ... By
 *    Lemma lem_nontriv_projection, there exists an O(eps)-projection P' in S_{P_m}
 *    such that both P' and P'' = Ptilde_m - P' are nonvanishing. ... Merging
 *    v_comm^(1) and v_comm^(2) using Corollary cor_merge_sum, we obtain an
 *    O(eps)-inclusion v_comm^+ ... Due to Corollary cor_improvement, there also
 *    exists a c_0*eps-inclusion ... But this contradicts the maximum dimensionality."
 *
 * approximate_algebras.tex:1428 — the equivalence-class partition (verbatim):
 *   "The indices j of the one-dimensional projections P_j are subdivided into
 *    equivalence classes such that dim S_{P_j,P_k} is 1 if j and k are in the same
 *    class and 0 otherwise."
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The paper's "maximum
 * dimensionality" v_comm is the ENDPOINT of the greedy loop; we CONSTRUCT it by
 * splitting any P_m with dim S_{P_m} > 1 into two nonvanishing projections, exactly
 * the contradiction step, until none can be split. cor_merge_sum is a pure
 * concatenation of basis-image lists; for the FLAT commutative skeleton (every block
 * is a 1x1 C) the concatenation of {P_1,...,P_{m-1},P',P''} IS the C^m inclusion
 * whose vE are those projections, so we rebuild v_comm directly from the projection
 * list (the equivalent of the §2.2 v_part1 (+) v_2 merge) and then cor_improvement
 * (aic_errreduce) restores the c_0*eps inclusion. No new numerical primitive.
 *
 * THE c_0 IS THE MEASURED FIRST-CALL CONSTANT (FINDINGS §D2/§C11). The analytic c_0
 * (bead aic-1bc) is unstated; we fix c0_nominal to the first aic_errreduce's measured
 * c_0 and GATE every subsequent call against 3x it + 5 (the eta~0 noise addend) — the
 * universality canary that catches the .tex:484 dimension-growth failure mode.
 *
 * FAIL LOUD (Rule 4). Stage 1 caps at dim_A + 1 iterations (the contradiction bound);
 * aic_projection_nontrivial aborts on a degenerate >1-dim sub-algebra (the aic-3qv
 * stop condition, FINDINGS §D1 — NOT a silent skip); the S_{P_m} nontriviality is
 * verified against the wrapper unit Ptilde_m (FINDINGS §C11, route (b)).
 */
#include <assert.h>
#include <stdio.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner.h"
#include "aic_cstar.h"
#include "aic_cstar_internal.h"
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_errreduce.h"
#include "aic_mat.h"
#include "aic_projection.h"

static double dmid(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* ||X - Y||_op via the §C5 midpoint route (drop the matmul-accumulated radii that
 * trip the Gram Hermiticity false-fail on near-zero-off-diagonal differences;
 * mirrors test_cstar_extension.c opnorm_diff_d). A coarse fail-loud gate value. */
static double opnorm_diff(const acb_mat_t X, const acb_mat_t Y, slong prec)
{
    slong r = acb_mat_nrows(X), c = acb_mat_ncols(X);
    acb_mat_t D, M;
    acb_mat_init(D, r, c);
    acb_mat_init(M, r, c);
    acb_mat_sub(D, X, Y, prec);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    double o = dmid(e);
    arb_clear(e);
    acb_mat_clear(M);
    acb_mat_clear(D);
    return o;
}

/* ---- projection-list management (the Stage-1 working set) --------------- */

/* Append a deep copy of P (n x n) at index list->m, growing the array. */
void aic_cstar_plist_push(aic_cstar_plist *list, const acb_mat_t P)
{
    slong i = list->m;
    if (i == list->cap) {
        slong newcap = (list->cap == 0) ? 4 : 2 * list->cap;
        list->P = flint_realloc(list->P, (size_t) newcap * sizeof(acb_mat_t));
        list->cap = newcap;
    }
    acb_mat_init(list->P[i], list->n, list->n);
    acb_mat_set(list->P[i], P);
    list->m = i + 1;
}

void aic_cstar_plist_init(aic_cstar_plist *list, slong n)
{
    list->P = NULL;
    list->m = 0;
    list->cap = 0;
    list->n = n;
}

void aic_cstar_plist_clear(aic_cstar_plist *list)
{
    for (slong i = 0; i < list->m; i++) acb_mat_clear(list->P[i]);
    if (list->P) flint_free(list->P);
    list->P = NULL;
    list->m = 0;
    list->cap = 0;
}

/* Replace P[m_idx] with Pp, then insert Ppp right after it (shifting the tail).
 * Net: the list grows by one, P_m becoming the pair {Pp, Ppp}. */
static void plist_split(aic_cstar_plist *list, slong m_idx,
                        const acb_mat_t Pp, const acb_mat_t Ppp)
{
    aic_cstar_plist_push(list, Ppp);   /* grow by one at the tail */
    /* shift entries (m_idx+1 .. m-2) up by one to open a slot at m_idx+1 */
    for (slong i = list->m - 1; i > m_idx + 1; i--)
        acb_mat_swap(list->P[i], list->P[i - 1]);
    acb_mat_set(list->P[m_idx], Pp);
    acb_mat_set(list->P[m_idx + 1], Ppp);
}

/* Build a flat commutative-skeleton inclusion v : C^m -> A from the projection
 * list: B = (+)_{j} M_1 (m one-dim blocks), v->vE[j] = list->P[j]. B_out, v_out
 * are OUTPUT (init'd here); caller frees both. */
static void build_vcomm(aic_dhom_B *B_out, aic_dhom_v *v_out,
                        const aic_cstar_plist *list, const aic_ecstar *A)
{
    slong m = list->m;
    slong *dims = flint_malloc((size_t) m * sizeof(slong));
    for (slong j = 0; j < m; j++) dims[j] = 1;
    aic_dhom_B_init(B_out, dims, m);
    flint_free(dims);
    aic_dhom_v_init(v_out, B_out, A);
    for (slong j = 0; j < m; j++) acb_mat_set(v_out->vE[j], list->P[j]);
}

/* Find the first projection index m with dim S_{P_m} > 1; -1 if all are 1d. */
static slong find_dim_gt1(const aic_ecstar *A, const aic_cstar_plist *list,
                          slong prec)
{
    for (slong m = 0; m < list->m; m++)
        if (aic_corner_dim_S(A, list->P[m], list->P[m], prec) > 1) return m;
    return -1;
}

/* ---- the c_0 gate (FINDINGS §D2 universality canary) -------------------- */

void aic_cstar_c0_gate(double c0_this, double c0_nominal, const char *where)
{
    /* 3x nominal + 5 addend (the eta~0 regime, spec §6). A drift past this is the
     * .tex:484 dimension-growth failure mode (stop condition, FINDINGS §D2). */
    AIC_FAIL_IF(!(c0_this <= 3.0 * c0_nominal + 5.0),
                "aic_cstar_build %s: c_0 drift: step c_0 = %.4f, nominal c_0 = %.4f "
                "(> 3x+5); universality canary fired (.tex:484 / FINDINGS §D2, bead "
                "aic-1bc)", where, c0_this, c0_nominal);
}

/* ---- Stage 1: the greedy commutative skeleton --------------------------- */

/* Fill `list` with the final one-dimensional projections {P_1,...,P_m} and write
 * the nominal c_0 (first cor_improvement) to *c0_nominal. The list's projections
 * are the errreduce'd v_comm.vE diagonal elements (spec §2.2). */
void aic_cstar_stage1_greedy(aic_cstar_plist *list, double *c0_nominal,
                             const aic_ecstar *A, double eps, double tol_abs,
                             slong max_steps, slong prec)
{
    slong n = A->n;
    aic_cstar_plist_init(list, n);

    /* Degenerate case dim A = 1: B = C, v(1) = unit of A (the ambient 1_n is A's
     * unit, .tex:2186). Exact iso, no split (spec §2.1). */
    if (A->dim_A == 1) {
        acb_mat_t I;
        acb_mat_init(I, n, n);
        acb_mat_one(I);
        aic_cstar_plist_push(list, I);
        acb_mat_clear(I);
        *c0_nominal = 0.0;
        return;
    }

    /* Initialization (spec §2.1): the FIRST nontrivial split P_1 = proj, P_2 =
     * 1_A - P_1 (A's unit is the ambient 1_n). */
    acb_mat_t P1, P2, I;
    acb_mat_init(P1, n, n);
    acb_mat_init(P2, n, n);
    acb_mat_init(I, n, n);
    arb_t delta_p, Pn, IPn;
    arb_init(delta_p); arb_init(Pn); arb_init(IPn);
    aic_projection_nontrivial(P1, delta_p, Pn, IPn, A, NULL, prec);
    acb_mat_one(I);
    acb_mat_sub(P2, I, P1, prec);
    aic_cstar_plist_push(list, P1);
    aic_cstar_plist_push(list, P2);
    arb_clear(IPn); arb_clear(Pn); arb_clear(delta_p);
    acb_mat_clear(I); acb_mat_clear(P2); acb_mat_clear(P1);

    int c0_fixed = 0;
    *c0_nominal = 0.0;

    /* The greedy loop. STAGE1_MAX_ITERS = dim_A + 1: each split adds one projection,
     * starting from 2, and at most dim_A 1-dim projections exist (.tex:1417). */
    slong iter = 0;
    const slong max_iters = A->dim_A + 1;
    for (;;) {
        /* (a) cor_improvement on the current flat skeleton (rebuild v_comm from the
         * list, then errreduce). Updates the list to the normalized projections. */
        aic_dhom_B B;
        aic_dhom_v v_raw, v_comm;
        build_vcomm(&B, &v_raw, list, A);
        aic_dhom_v_init(&v_comm, &B, A);
        arb_t c0;
        arb_init(c0);
        aic_errreduce(&v_comm, c0, NULL, &v_raw, eps, tol_abs, max_steps, prec);
        double c0d = dmid(c0);
        if (!c0_fixed) { *c0_nominal = c0d; c0_fixed = 1; }
        else aic_cstar_c0_gate(c0d, *c0_nominal, "stage1");
        /* write back the normalized projections into the list */
        for (slong j = 0; j < list->m; j++) acb_mat_set(list->P[j], v_comm.vE[j]);
        arb_clear(c0);
        aic_dhom_v_clear(&v_comm);
        aic_dhom_v_clear(&v_raw);
        aic_dhom_B_clear(&B);

        /* (b) any P_m with dim S_{P_m} > 1?  If not, the skeleton is maximal. */
        slong m = find_dim_gt1(A, list, prec);
        if (m < 0) break;

        AIC_FAIL_IF(++iter > max_iters,
                    "aic_cstar_build stage1: loop exceeded dim_A + 1 = %ld iterations "
                    "(contradiction: th_main Stage 1 termination violated, .tex:1417; "
                    "bead aic-3qv stop condition)", (long) max_iters);

        /* (c) split P_m: P' in S_{P_m}, P'' = Ptilde_m - P'. */
        aic_cstar_subalgebra sp;
        aic_cstar_subalg_build(&sp, A, list->P[m], prec);
        acb_mat_t Pp, Ppp;
        acb_mat_init(Pp, n, n);
        acb_mat_init(Ppp, n, n);
        arb_t dp, pn, ipn;
        arb_init(dp); arb_init(pn); arb_init(ipn);
        aic_projection_nontrivial(Pp, dp, pn, ipn, &sp.A, NULL, prec);

        /* Nontriviality vs the WRAPPER unit Ptilde_m (FINDINGS §C11 route (b)): the
         * ambient ||1_n - P'|| the projector reports is vacuous on a corner element.
         * Check ||Ptilde_m - P'||_op >= 0.15 (so P' != Ptilde_m) AND ||P'|| >= 0.15
         * (P' nonvanishing). */
        double pdiff = opnorm_diff(sp.Ptilde, Pp, prec);
        AIC_FAIL_IF(!(pdiff >= 0.15),
                    "aic_cstar_build stage1: split projection P' is the wrapper unit "
                    "Ptilde_m (||Ptilde_m - P'|| = %.4f < 0.15; FINDINGS §C11)", pdiff);
        AIC_FAIL_IF(!(dmid(pn) >= 0.15),
                    "aic_cstar_build stage1: split projection P' vanishing "
                    "(||P'|| = %.4f < 0.15)", dmid(pn));

        acb_mat_sub(Ppp, sp.Ptilde, Pp, prec);   /* P'' = Ptilde_m - P' */
        plist_split(list, m, Pp, Ppp);

        arb_clear(ipn); arb_clear(pn); arb_clear(dp);
        acb_mat_clear(Ppp); acb_mat_clear(Pp);
        aic_cstar_subalg_clear(&sp);
    }
}

/* ---- the public driver -------------------------------------------------- */

void aic_cstar_build(aic_dhom_B *B_out, aic_dhom_v *v_out, arb_t iso_def,
                     const aic_ecstar *A, double eps, slong prec)
{
    assert(B_out != NULL && v_out != NULL && A != NULL);
    assert(A->n >= 1 && "aic_cstar_build: A->n must be >= 1");
    assert(A->dim_A >= 1 && "aic_cstar_build: A->dim_A must be >= 1");

    const double tol_abs = 1e-13;
    const slong max_steps = 40;
    /* eps' = O(eps): the S_P sub-algebras are eps'-C* with eps' >= c_0*eps; the
     * floor (eps and the wrapper-defect constant ~2.7) is captured by this factor
     * (spec §4.4). Bare eps would stall the per-class errreduce (the wrapper defect
     * cannot beat its O(eta) obstruction). */
    const double eps_prime = (eps > tol_abs) ? AIC_STAGE2_EPS_PRIME_FACTOR * eps : eps;

    /* Stage 1: the maximal commutative skeleton {P_1,...,P_m}. */
    aic_cstar_plist list;
    double c0_nominal = 0.0;
    aic_cstar_stage1_greedy(&list, &c0_nominal, A, eps, tol_abs, max_steps, prec);
    slong m = list.m;

    /* Equivalence-class partition (.tex:1428). */
    slong *class_id = flint_malloc((size_t) m * sizeof(slong));
    slong num_classes = 0;
    aic_cstar_compute_classes(class_id, &num_classes, A, &list, prec);

    /* Stages 2+3: per-class induction then merge. Fills B_out, v_out, iso_def. */
    aic_cstar_stages_run(B_out, v_out, iso_def, A, &list, class_id, num_classes,
                         c0_nominal, eps, eps_prime, tol_abs, max_steps, prec);

    flint_free(class_id);
    aic_cstar_plist_clear(&list);
}
