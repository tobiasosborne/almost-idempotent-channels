/* aic_cstar_stages.c — Stage 2 (per-class inductive extension) and Stage 3 (class
 * merge) of the th_main master loop, plus the equivalence-class partition (bead
 * aic-097, module "cstar_build", Increment 5). Split from src/aic_cstar_build.c
 * (the driver + Stage 1) under Rule 10. See include/aic_cstar.h (Increment 5
 * banner) and docs/research/cstar_masterloop_spec.md §3-5.
 *
 * approximate_algebras.tex:1430-1441 — Stage 2 (verbatim):
 *   "Each equivalence class C is considered separately. ... We construct c_0*eps'-
 *    isomorphisms v_r : M_r -> A_r for r=1,...,s inductively. To obtain v_r from
 *    v_{r-1}, we first fit v_{r-1}, P_{[1,r-1]}, and P_r into A_r by applying the
 *    compression map Co_{P_{[1,r]}}: v_{r-1}'(X) = Co_{P_{[1,r]}}(v(X)), P_{[1,r-1]}'
 *    = Co_{P_{[1,r]}}(P_{[1,r-1]}), P_r' = Co_{P_{[1,r]}}(P_r). Then we extend the
 *    O(eps)-isomorphism v_{r-1}' to v_{r-1}^+ using Lemma lem_extension. Finally, we
 *    use Corollary cor_improvement to replace v_{r-1}^+ with a c_0*eps'-isomorphism v_r."
 *
 * approximate_algebras.tex:1443 — Stage 3 (verbatim):
 *   "we have constructed c_0*eps'-isomorphisms v_C : M_{|C|} -> S_{P_C} for all
 *    equivalence classes C. Note that by Lemma lem_add_dim, S_{P_C,P_D}=0 if the
 *    classes C and D are distinct. This allows us to successively merge the
 *    isomorphisms v_C ... Each step includes the application of Corollary
 *    cor_merge_sum followed by the use of Corollary cor_improvement to reduce errors."
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The compression
 * Co_{P_{[1,r]}} is aic_corner_apply over the wrapper's prebuilt Co_P; lem_extension
 * is the closed I4 routine (its oblique-wrapper-as-parent path unblocked by aic-qgs,
 * FINDINGS §C10 RESOLVED); cor_merge_sum is the closed I2 concat; cor_improvement is
 * aic_errreduce. The equivalence classes are union-find on dim S_{P_j,P_k}==1. No new
 * numerical primitive is introduced; this is lifetime + bookkeeping glue.
 */
#include <assert.h>
#include <stdio.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_corner.h"
#include "aic/aic_cstar.h"
#include "aic/aic_cstar_internal.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_errreduce.h"

static double smid(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }


/* ---- Stage 2: build v_C : M_{|C|} -> S_{P_C} for one class -------------- */

/* The sorted list of projection indices in class `cid`. Writes idx[0..*s-1]. */
static void class_members(slong *idx, slong *s, const slong *class_id, slong m,
                          slong cid)
{
    slong c = 0;
    for (slong j = 0; j < m; j++) if (class_id[j] == cid) idx[c++] = j;
    *s = c;
}

/* Build v_C (a c_0*eps'-isomorphism M_s -> S_{P_C}) for class `cid`. v_C->A = A
 * (the ambient parent; the wrapper is internal). B_C, v_C are OUTPUT (init'd here,
 * caller frees both). P_C is written (n x n, the class projection sum). */
static void stage2_class(aic_dhom_B *B_C, aic_dhom_v *v_C, acb_mat_t P_C,
                         const aic_ecstar *A, const aic_cstar_plist *list,
                         const slong *class_id, slong cid, double c0_nominal,
                         double eps_prime, double tol_abs, slong max_steps, slong prec)
{
    slong n = A->n, m = list->m;
    slong *idx = flint_malloc((size_t) m * sizeof(slong));
    slong s;
    class_members(idx, &s, class_id, m, cid);
    assert(s >= 1 && "stage2_class: empty class");

    /* Running v_{r-1} : M_{r-1} -> S_{P_sum_{r-1}}, with A_parent = A throughout
     * (the wrapper is built per step; lem_extension's parent is the wrapper). We
     * keep v over the wrapper at each step. Base r=1: v_1 over the wrapper A_1. */
    aic_cstar_subalgebra sp;
    aic_cstar_subalg_build(&sp, A, list->P[idx[0]], prec);  /* A_1 = S_{P_{r_0}} */

    aic_dhom_B B_cur;
    aic_dhom_v v_cur;
    slong dims1[1] = {1};
    aic_dhom_B_init(&B_cur, dims1, 1);
    aic_dhom_v_init(&v_cur, &B_cur, &sp.A);
    aic_corner_Ptilde(v_cur.vE[0], &sp.A, list->P[idx[0]], prec);  /* v_1(1)=Ptilde */

    /* errreduce the base (a genuine 1-d C* algebra; c_0 ~ machine-zero). The
     * codomain is the wrapper sp.A whose unit is sp.Ptilde (NOT 1_n) — use
     * errreduce_unit (FINDINGS §C7). */
    {
        aic_dhom_v vt;
        aic_dhom_v_init(&vt, &B_cur, &sp.A);
        arb_t c0;
        arb_init(c0);
        aic_cstar_errreduce_unit(&vt, c0, &v_cur, sp.Ptilde, eps_prime, tol_abs, max_steps, prec);
        aic_cstar_c0_gate(smid(c0), c0_nominal, "stage2-base");
        aic_dhom_v_clear(&v_cur);
        v_cur = vt;            /* MOVE: vt's vE now owned by v_cur */
        arb_clear(c0);
    }

    acb_mat_t P_sum;           /* running P_{[1,r]} (ambient n x n) */
    acb_mat_init(P_sum, n, n);
    acb_mat_set(P_sum, list->P[idx[0]]);

    for (slong r = 2; r <= s; r++) {
        slong rj = idx[r - 1];               /* the r-th class projection P_r */
        acb_mat_t P_sum_new;
        acb_mat_init(P_sum_new, n, n);
        acb_mat_add(P_sum_new, P_sum, list->P[rj], prec);  /* P_{[1,r]} */

        /* A_r = S_{P_{[1,r]}} wrapper; Co_{P_{[1,r]}} = sp_r.ctx->Co_P. */
        aic_cstar_subalgebra sp_r;
        aic_cstar_subalg_build(&sp_r, A, P_sum_new, prec);

        /* Compress v_{r-1} into A_r (.tex:1437): vr_prime(E_i) = Co_{P_[1,r]}(v(E_i)).
         * v_cur lives over A_{r-1}=sp.A but its vE are n x n operators in A; the
         * compression re-targets them into S_{P_[1,r]}. NOTE the compression Co_P =
         * sp_r.ctx->Co_P is in the PARENT'S coordinate space (parent->dim_A x
         * parent->dim_A), so aic_corner_apply takes the PARENT A, NOT the wrapper
         * sp_r.A (the spec §4.2 wrote &sp_r.A; that mis-sizes T — corrected here,
         * the result is an n x n operator in S_{P_[1,r]} <= A either way). */
        aic_dhom_v vr_prime;
        aic_dhom_v_init(&vr_prime, v_cur.B, &sp_r.A);
        for (slong i = 0; i < v_cur.B->dim_B; i++)
            aic_corner_apply(vr_prime.vE[i], A, sp_r.ctx->Co_P, v_cur.vE[i], prec);

        /* Compress the projections (.tex:1438-1439). */
        acb_mat_t Pm1_prime, Pr_prime;
        acb_mat_init(Pm1_prime, n, n);
        acb_mat_init(Pr_prime, n, n);
        aic_corner_apply(Pm1_prime, A, sp_r.ctx->Co_P, P_sum, prec);
        aic_corner_apply(Pr_prime, A, sp_r.ctx->Co_P, list->P[rj], prec);

        /* lem_extension: M_{r-1} -> M_r over A_r (.tex:1441). */
        aic_dhom_B B_plus;
        aic_dhom_v v_plus;
        arb_t md, sm;
        arb_init(md); arb_init(sm);
        aic_cstar_lem_extension(&B_plus, &v_plus, md, sm, &vr_prime,
                                Pm1_prime, Pr_prime, &sp_r.A, prec);
        arb_clear(sm); arb_clear(md);

        /* cor_improvement -> v_r (c_0*eps'-iso). Codomain wrapper sp_r.A, unit
         * sp_r.Ptilde (FINDINGS §C7). */
        aic_dhom_v v_r;
        aic_dhom_v_init(&v_r, &B_plus, &sp_r.A);
        arb_t c0;
        arb_init(c0);
        aic_cstar_errreduce_unit(&v_r, c0, &v_plus, sp_r.Ptilde, eps_prime, tol_abs, max_steps, prec);
        aic_cstar_c0_gate(smid(c0), c0_nominal, "stage2-step");
        arb_clear(c0);

        /* Advance: v_cur <- v_r (over sp_r.A); B_cur <- B_plus; sp <- sp_r. The
         * struct copies make v_r->B (= &B_plus, loop-local) and v_r->A (= &sp_r.A,
         * loop-local) STALE; re-point them at the persistent &B_cur / &sp.A after
         * the moves (the same move-aliasing fix as Stage 3). */
        aic_dhom_v_clear(&v_plus);
        aic_dhom_v_clear(&vr_prime);
        acb_mat_clear(Pr_prime); acb_mat_clear(Pm1_prime);
        aic_dhom_v_clear(&v_cur);
        aic_dhom_B_clear(&B_cur);
        aic_cstar_subalg_clear(&sp);
        sp = sp_r;
        B_cur = B_plus;
        v_cur = v_r;
        v_cur.B = &B_cur;
        v_cur.A = &sp.A;
        acb_mat_swap(P_sum, P_sum_new);
        acb_mat_clear(P_sum_new);
    }

    /* v_C : M_s -> S_{P_C} <= A. Re-present over the AMBIENT A: the vE are n x n
     * operators in A already (S_{P_C} <= A), so copy them into a fresh v over A. */
    slong dimsC[1] = {s};
    aic_dhom_B_init(B_C, dimsC, 1);
    aic_dhom_v_init(v_C, B_C, A);
    for (slong i = 0; i < B_C->dim_B; i++) acb_mat_set(v_C->vE[i], v_cur.vE[i]);
    acb_mat_set(P_C, P_sum);

    acb_mat_clear(P_sum);
    aic_dhom_v_clear(&v_cur);
    aic_dhom_B_clear(&B_cur);
    aic_cstar_subalg_clear(&sp);
    flint_free(idx);
}

/* ---- Stages 2 + 3 driver ------------------------------------------------ */

void aic_cstar_stages_run(aic_dhom_B *B_out, aic_dhom_v *v_out, arb_t iso_def,
                          const aic_ecstar *A, const aic_cstar_plist *list,
                          const slong *class_id, slong num_classes,
                          double c0_nominal, double eps, double eps_prime,
                          double tol_abs, slong max_steps, slong prec)
{
    slong n = A->n;
    (void) eps;

    /* Stage 2: per-class isomorphisms v_C : M_{|C|} -> S_{P_C}. */
    aic_dhom_B *B_arr = flint_malloc((size_t) num_classes * sizeof(aic_dhom_B));
    aic_dhom_v *v_arr = flint_malloc((size_t) num_classes * sizeof(aic_dhom_v));
    acb_mat_t *P_arr = flint_malloc((size_t) num_classes * sizeof(acb_mat_t));
    for (slong c = 0; c < num_classes; c++) {
        acb_mat_init(P_arr[c], n, n);
        stage2_class(&B_arr[c], &v_arr[c], P_arr[c], A, list, class_id, c,
                     c0_nominal, eps_prime, tol_abs, max_steps, prec);
    }

    /* Stage 3: successively merge classes (.tex:1443). Start with class 0. */
    aic_dhom_B B_merged;
    aic_dhom_v v_merged;
    aic_dhom_B_init(&B_merged, B_arr[0].d, B_arr[0].num_blocks);
    aic_dhom_v_init(&v_merged, &B_merged, A);
    for (slong i = 0; i < B_arr[0].dim_B; i++) acb_mat_set(v_merged.vE[i], v_arr[0].vE[i]);

    acb_mat_t P_total;
    acb_mat_init(P_total, n, n);
    acb_mat_set(P_total, P_arr[0]);

    for (slong c = 1; c < num_classes; c++) {
        /* off-diagonal zero vs the running total (.tex:1443, lem_add_dim). */
        AIC_FAIL_IF(!aic_cstar_off_diag_zero(A, P_total, P_arr[c], prec),
                    "aic_cstar_build stage3: dim S_{P_total,P_C} != 0 merging class "
                    "%ld (cor_merge_sum bijectivity hypothesis, .tex:1443)", (long) c);

        aic_dhom_B B_new;
        aic_dhom_v v_new;
        aic_cstar_merge_sum(&B_new, &v_new, NULL, NULL, &v_merged, &v_arr[c], A, prec);

        /* cor_improvement -> back to c_0*eps'. The merged map includes into the
         * RUNNING sub-algebra S_{P_total + P_C}, whose unit is P_total + P_C (NOT the
         * ambient 1_n until ALL classes are merged) — so certify against that running
         * unit (FINDINGS §C7; the first merge's v(I_B) = sum of two class units != 1_n
         * spuriously trips the stock unit_def gate, exactly as in Stage 2). */
        acb_mat_t unit_run;
        acb_mat_init(unit_run, n, n);
        acb_mat_add(unit_run, P_total, P_arr[c], prec);
        aic_dhom_v v_next;
        aic_dhom_v_init(&v_next, &B_new, A);
        arb_t c0;
        arb_init(c0);
        aic_cstar_errreduce_unit(&v_next, c0, &v_new, unit_run, eps_prime, tol_abs, max_steps, prec);
        aic_cstar_c0_gate(smid(c0), c0_nominal, "stage3-merge");
        arb_clear(c0);
        acb_mat_clear(unit_run);

        /* advance: v_merged <- v_next (over B_new); P_total += P_C. The struct copy
         * B_merged = B_new makes v_merged->B (= &B_new, a loop-local) STALE, so
         * re-point it at the persistent &B_merged. */
        aic_dhom_v_clear(&v_new);
        aic_dhom_v_clear(&v_merged);
        aic_dhom_B_clear(&B_merged);
        B_merged = B_new;
        v_merged = v_next;
        v_merged.B = &B_merged;
        acb_mat_add(P_total, P_total, P_arr[c], prec);
    }
    acb_mat_clear(P_total);

    /* Final v : B -> A. ASSERT bijective (.tex:1318). */
    arb_t a_bij;
    arb_init(a_bij);
    int bij = aic_errreduce_is_bijective(a_bij, &v_merged, prec);
    AIC_FAIL_IF(!bij,
                "aic_cstar_build: final v not bijective (sigma_min = %.4f, dim_B = "
                "%ld, dim_A = %ld); stop condition (FINDINGS §D1/§D2)",
                smid(a_bij), (long) B_merged.dim_B, (long) A->dim_A);
    arb_clear(a_bij);

    /* iso_def = max of the four inclusion-defects of the final v (the certificate). */
    if (iso_def != NULL) {
        aic_dhom_v vt;
        aic_dhom_v_init(&vt, &B_merged, A);
        aic_errreduce_defects defs;
        aic_errreduce_defects_init(&defs);
        arb_t c0;
        arb_init(c0);
        aic_errreduce(&vt, c0, &defs, &v_merged, eps_prime, tol_abs, max_steps, prec);
        aic_errreduce_defects_max(iso_def, &defs);
        arb_clear(c0);
        aic_errreduce_defects_clear(&defs);
        aic_dhom_v_clear(&vt);
    }

    /* hand off B_merged / v_merged as the outputs (MOVE). */
    *B_out = B_merged;
    *v_out = v_merged;
    v_out->B = B_out;   /* re-point at the OWNED B_out (B_merged was the stack copy) */

    /* cleanup the per-class arrays. */
    for (slong c = 0; c < num_classes; c++) {
        aic_dhom_v_clear(&v_arr[c]);
        aic_dhom_B_clear(&B_arr[c]);
        acb_mat_clear(P_arr[c]);
    }
    flint_free(P_arr);
    flint_free(v_arr);
    flint_free(B_arr);
}
