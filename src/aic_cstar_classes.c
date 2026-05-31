/* aic_cstar_classes.c — the equivalence-class partition (.tex:1428) and the
 * wrapper-unit-aware error reduction (cor_improvement for a map into an S_P
 * sub-algebra), shared by Stage 2 and Stage 3 of the th_main master loop (bead
 * aic-097, module "cstar_build", Increment 5). Split from src/aic_cstar_stages.c
 * under Rule 10. See include/aic_cstar.h (Increment 5 banner) and
 * docs/research/cstar_masterloop_spec.md §3-6.
 *
 * approximate_algebras.tex:1428 — the equivalence-class partition (verbatim):
 *   "The indices j of the one-dimensional projections P_j are subdivided into
 *    equivalence classes such that dim S_{P_j,P_k} is 1 if j and k are in the same
 *    class and 0 otherwise."  We realize this as union-find over the m projections
 *    with the relation dim S_{P_j,P_k} == 1 (aic_corner_equiv_1d), then certify the
 *    Stage-3 precondition S_{P_C,P_D}=0 for distinct classes (.tex:1443, lem_add_dim).
 *
 * THE WRAPPER-UNIT ERROR REDUCTION (FINDINGS §C7, the load-bearing subtlety). The
 * stock aic_errreduce certifies unit_def = ||v(I_B) - 1_n||_op against the AMBIENT
 * unit. That is correct for the final v : B -> A (whose image spans A, unit 1_n),
 * but WRONG for the intermediate maps into a compressed sub-algebra S_P (Stage 2:
 * v_r -> S_{P_[1,r]}; Stage 3: the running merge -> S_{P_total}), whose unit is the
 * compressed Ptilde (NOT 1_n) — there ||v(I_B) - 1_n|| ~ 1 spuriously trips the
 * AIC_ERRREDUCE_C0_CERT gate. errreduce_unit takes the codomain's actual unit and
 * certifies against it, REUSING every primitive aic_errreduce uses (no reinvention).
 */
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

static double smid(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* cor_improvement (.tex:1317) for a map into an S_P sub-algebra, unit = `unit`
 * (Ptilde for a wrapper / the running P_total for a Stage-3 merge), NOT 1_n.
 * REUSES aic_dhom_approx (lem_approx Newton), aic_dhom_defect_sweep, aic_dhom_v_
 * sigma_min, aic_dhom_prop_bounds (norm_ub) and certifies the four inclusion-
 * defects with unit_def = ||v(I_B) - unit||_op (the §C5 midpoint route, to dodge
 * the Gram false-fail on the near-zero difference). FINDINGS §C7.
 *   vt_out : OUTPUT (init'd over v_in's B and A by the caller); = lem_approx(v_in).
 *   c0_out : OUTPUT measured c_0 (caller-init'd arb_t).
 *   unit   : the codomain algebra's unit operator (n x n). */
void aic_cstar_errreduce_unit(aic_dhom_v *vt_out, arb_t c0_out,
                              const aic_dhom_v *v_in, const acb_mat_t unit,
                              double eps, double tol_abs, slong max_steps, slong prec)
{
    slong n = v_in->n;

    /* (A) input must be a delta-inclusion: coordinate-matrix sigma_min >= 0.5. */
    arb_t a_in;
    arb_init(a_in);
    aic_dhom_v_sigma_min(a_in, v_in, prec);
    AIC_FAIL_IF(arf_cmp_d(arb_midref(a_in), 0.5) < 0,
                "errreduce_unit: input not a delta-inclusion (sigma_min = %.4e < 0.5)",
                smid(a_in));
    arb_clear(a_in);

    /* (B) vtilde = lem_approx(v_in). unit_tol = 2.0 so the legitimate
     * ||v(I_B) - 1_n|| ~ 1 wrapper drift does NOT trip dhom_approx's BLOWUP guard
     * (the Newton step itself is unit-independent). */
    aic_errreduce_copy_v(vt_out, v_in);
    aic_dhom_approx(vt_out, AIC_ERRREDUCE_EPS_FLOOR * eps, tol_abs, /*unit_tol=*/2.0,
                    max_steps, prec, NULL);

    /* (C) the four inclusion-defects, unit_def against the SUPPLIED unit. */
    arb_t mult, nub, nlb;
    arb_init(mult); arb_init(nub); arb_init(nlb);
    aic_dhom_defect_sweep(mult, vt_out, prec);
    aic_dhom_prop_bounds(nub, NULL, NULL, vt_out, prec);
    aic_dhom_v_sigma_min(nlb, vt_out, prec);

    acb_mat_t IB, vIB, D, M;
    acb_mat_init(IB, vt_out->B->n_B, vt_out->B->n_B);
    acb_mat_init(vIB, n, n);
    acb_mat_init(D, n, n);
    acb_mat_init(M, n, n);
    aic_dhom_B_unit(IB, vt_out->B);
    aic_dhom_v_apply(vIB, vt_out, IB, prec);
    acb_mat_sub(D, vIB, unit, prec);                 /* v(I_B) - unit */
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    arb_t udef;
    arb_init(udef);
    aic_mat_opnorm(udef, M, prec);                   /* §C5 midpoint route */

    /* max of {mult, max(0,||v||-1), max(0,1-sigma_min), unit_def}. */
    arb_t maxdef, tmp, one;
    arb_init(maxdef); arb_init(tmp); arb_init(one);
    arb_one(one);
    arb_set(maxdef, mult);
    arb_sub(tmp, nub, one, prec); if (arb_is_negative(tmp)) arb_zero(tmp);
    if (arb_gt(tmp, maxdef)) arb_set(maxdef, tmp);
    arb_sub(tmp, one, nlb, prec); if (arb_is_negative(tmp)) arb_zero(tmp);
    if (arb_gt(tmp, maxdef)) arb_set(maxdef, tmp);
    if (arb_gt(udef, maxdef)) arb_set(maxdef, udef);

    /* (D) certify O(eps)-inclusion + (E) measured c_0. */
    double eps_eff = eps > tol_abs ? eps : tol_abs;
    AIC_FAIL_IF(!(smid(maxdef) <= AIC_ERRREDUCE_C0_CERT * eps_eff),
                "errreduce_unit: output not an O(eps)-inclusion (max-defect %.4e > "
                "%.1f*max(eps,tol)=%.4e; eps=%.3e)", smid(maxdef),
                AIC_ERRREDUCE_C0_CERT, AIC_ERRREDUCE_C0_CERT * eps_eff, eps);
    if (eps <= tol_abs) arb_zero(c0_out);
    else { arb_set_d(tmp, eps); arb_div(c0_out, maxdef, tmp, prec); }

    arb_clear(one); arb_clear(tmp); arb_clear(maxdef);
    arb_clear(udef);
    acb_mat_clear(M); acb_mat_clear(D); acb_mat_clear(vIB); acb_mat_clear(IB);
    arb_clear(nlb); arb_clear(nub); arb_clear(mult);
}

/* ---- union-find for the equivalence-class partition (.tex:1428) --------- */

static slong uf_find(slong *parent, slong j)
{
    while (parent[j] != j) { parent[j] = parent[parent[j]]; j = parent[j]; }
    return j;
}

static void uf_union(slong *parent, slong j, slong k)
{
    slong rj = uf_find(parent, j), rk = uf_find(parent, k);
    if (rj != rk) parent[rj < rk ? rk : rj] = (rj < rk ? rj : rk);
}

/* Sum the class projections into P_C (n x n), for class index `cid`. */
static void sum_class(acb_mat_t P_C, const aic_cstar_plist *list,
                      const slong *class_id, slong cid, slong prec)
{
    acb_mat_zero(P_C);
    for (slong j = 0; j < list->m; j++)
        if (class_id[j] == cid) acb_mat_add(P_C, P_C, list->P[j], prec);
}

void aic_cstar_compute_classes(slong *class_id, slong *num_classes,
                               const aic_ecstar *A, const aic_cstar_plist *list,
                               slong prec)
{
    slong m = list->m;
    slong *parent = flint_malloc((size_t) m * sizeof(slong));
    for (slong j = 0; j < m; j++) parent[j] = j;

    /* Union j,k iff dim S_{P_j,P_k} == 1 (.tex:1428). */
    for (slong j = 0; j < m; j++)
        for (slong k = j + 1; k < m; k++)
            if (aic_corner_equiv_1d(A, list->P[j], list->P[k], prec))
                uf_union(parent, j, k);

    /* Compact roots into [0, num_classes). */
    slong *root_to_cid = flint_malloc((size_t) m * sizeof(slong));
    for (slong j = 0; j < m; j++) root_to_cid[j] = -1;
    slong nc = 0;
    for (slong j = 0; j < m; j++) {
        slong r = uf_find(parent, j);
        if (root_to_cid[r] < 0) root_to_cid[r] = nc++;
        class_id[j] = root_to_cid[r];
    }
    *num_classes = nc;

    /* Stage-3 precondition (.tex:1443, lem_add_dim): S_{P_C,P_D}=0 for C != D. */
    slong n = A->n;
    acb_mat_t P_C, P_D;
    acb_mat_init(P_C, n, n);
    acb_mat_init(P_D, n, n);
    for (slong c = 0; c < nc; c++)
        for (slong dcl = c + 1; dcl < nc; dcl++) {
            sum_class(P_C, list, class_id, c, prec);
            sum_class(P_D, list, class_id, dcl, prec);
            AIC_FAIL_IF(!aic_cstar_off_diag_zero(A, P_C, P_D, prec),
                        "aic_cstar_build classes: dim S_{P_C,P_D} != 0 for distinct "
                        "classes %ld,%ld (partition inconsistent, .tex:1443)",
                        (long) c, (long) dcl);
        }
    acb_mat_clear(P_D);
    acb_mat_clear(P_C);

    flint_free(root_to_cid);
    flint_free(parent);
}
