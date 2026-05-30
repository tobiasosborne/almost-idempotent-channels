/* aic_errreduce.c — cor_improvement (error reduction), bead aic-t81. Realizes
 * approximate_algebras.tex:1317-1319. THIN orchestration of the dhom module (§8).
 *
 * THE ALGORITHM (.tex:1317, "a straightforward corollary of lem_approx and
 * prop_delta_hominc"):
 *   Input: a delta-INCLUSION v: B -> A (a delta-hom that is ALSO bounded below,
 *          (1-delta)||X|| <= ||v(X)||, .tex:451-453).
 *   1. vtilde = lem_approx(v)  (aic_dhom_approx; an O(eps)-HOMOMORPHISM whose
 *      multiplicativity defect is O(eps), .tex:1225, and ||vtilde - v|| <= O(delta)).
 *   2. CERTIFY vtilde is an O(eps)-INCLUSION via prop_delta_hominc (.tex:1194)
 *      applied to vtilde:
 *        (i)   ||vtilde|| <= 1 + O(eps)                  -> norm_excess <= O(eps)
 *        (ii)  ||vtilde(X)|| >= (1 - O(eps))||X||         -> lower_gap   <= O(eps)
 *              [hypothesis a > 2*delta_vtilde holds: v was a delta-inclusion so
 *               ||v(X)|| >= (1-delta)||X||, and ||vtilde - v|| <= O(delta), so
 *               vtilde inherits a >= (1 - O(delta)) > 2*O(eps)]
 *        (iii) ||vtilde(I_B) - I_A|| <= O(eps)            -> unit_def    <= O(eps)
 *      together with the multiplicativity defect (i'): O(eps).
 *   3. c_0 = MEASURED max(inclusion-defects) / eps (the analytic c_0 is aic-1bc;
 *      paper/FINDINGS.md D2).
 *   4. Bijectivity: dim_B == dim_A AND vtilde injective (lower bound > 0) =>
 *      vtilde is a c_0*eps-ISOMORPHISM (.tex:1318; lem_approx preserves the lower
 *      bound, dim match closes bijectivity in finite dim).
 *
 * WHY THE c_0 IS DIMENSION-INDEPENDENT (.tex:484/461, the universality canary).
 * The diagonal that drives lem_approx is the EXACT generalized-Pauli diagonal of
 * the GENUINE C* algebra B, NOT a naive Haar diagonal of the eps-associative A
 * (whose error is proportional to n = dim A, .tex:484). So the measured c_0 must
 * NOT grow with dim B / dim A; test_errreduce T4 asserts this (a growing c_0 is a
 * stop condition, escalate aic-1bc).
 *
 * FAIL-LOUD (Rule 4). Aborts if (a) the input is NOT a delta-inclusion (its
 * coordinate-matrix sigma_min < 0.5 — a SOUND collapse detector, F1 fix below), or
 * (b) the output is NOT certifiably an O(eps)-inclusion (max inclusion-defect
 * exceeds AIC_ERRREDUCE_C0_CERT * max(eps,tol)).
 *
 * THE LOWER BOUND a IS sigma_min OF THE COORDINATE MATRIX (F1 fix; paper/
 * FINDINGS.md). v(X) = sum_i x_i vE[i] is linear; M[k,i] = <B_k, vE[i]>_F (A's
 * Frobenius-orthonormal basis {B_k}, B's Frobenius-orthonormal matrix units {E_i})
 * gives ||v(X)||_F = ||M x||_2, ||X||_F = ||x||_2, so sigma_min(M) = inf_{X!=0}
 * ||v(X)||_F/||X||_F — the TRUE Frobenius unit-ball inclusion infimum (a SOUND
 * collapse detector: 0 iff v collapses a direction). The OLD basis-sweep lower
 * bound min_i ||vE[i]|| is BLIND to a v bounded below on each basis element that
 * collapses a combination (the F1 BLOCKER). sigma_min is a COORDINATE/Frobenius
 * inf (differs from the exact OPERATOR-norm inclusion inf by <= sqrt(dim) factors;
 * the faithful operator-norm HOPM is a later cycle), but it correctly REJECTS a
 * non-inclusion. ||vtilde|| (norm_excess) and the multiplicativity defect remain
 * basis-sweep estimates (documented dhom limitation; aic_dhom.h CAVEAT).
 *
 * No eigendecomposition: products + certified op-norm upper bounds, all via dhom.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_dhom.h"
#include "aic_errreduce.h"

/* ------------------------------------------------------------------ defects */

void aic_errreduce_defects_init(aic_errreduce_defects *d)
{
    arb_init(d->mult_def);
    arb_init(d->norm_excess);
    arb_init(d->lower_gap);
    arb_init(d->unit_def);
}

void aic_errreduce_defects_clear(aic_errreduce_defects *d)
{
    arb_clear(d->mult_def);
    arb_clear(d->norm_excess);
    arb_clear(d->lower_gap);
    arb_clear(d->unit_def);
}

void aic_errreduce_defects_max(arb_t out, const aic_errreduce_defects *d)
{
    arb_set(out, d->mult_def);
    if (arb_gt(d->norm_excess, out)) arb_set(out, d->norm_excess);
    if (arb_gt(d->lower_gap, out)) arb_set(out, d->lower_gap);
    if (arb_gt(d->unit_def, out)) arb_set(out, d->unit_def);
}

/* ------------------------------------------------------------------- copy v */

void aic_errreduce_copy_v(aic_dhom_v *v_out, const aic_dhom_v *v_in)
{
    assert(v_out->B == v_in->B && v_out->A == v_in->A && v_out->n == v_in->n);
    for (slong i = 0; i < v_in->B->dim_B; i++)
        acb_mat_set(v_out->vE[i], v_in->vE[i]);
}

/* ----------------------------------------------------- positive-part helper */

/* out = max(0, x - 1): the positive part of (x - 1), as a certified arb ball.
 * Used for norm_excess = max(0, ||vtilde|| - 1). */
static void pos_part_sub1(arb_t out, const arb_t x, slong prec)
{
    arb_t one;
    arb_init(one);
    arb_one(one);
    arb_sub(out, x, one, prec);     /* x - 1 */
    if (arb_is_negative(out)) arb_zero(out);
    arb_clear(one);
}

/* out = max(0, 1 - x): the positive part of (1 - x). Used for lower_gap. */
static void pos_part_1sub(arb_t out, const arb_t x, slong prec)
{
    arb_t one;
    arb_init(one);
    arb_one(one);
    arb_sub(out, one, x, prec);     /* 1 - x */
    if (arb_is_negative(out)) arb_zero(out);
    arb_clear(one);
}

/* ----------------------------------------------------------- the main entry */

void aic_errreduce(aic_dhom_v *vtilde_out, arb_t c0_out,
                   aic_errreduce_defects *defects, const aic_dhom_v *v_in,
                   double eps, double tol_abs, slong max_steps, slong prec)
{
    assert(vtilde_out != NULL && v_in != NULL && c0_out != NULL);
    assert(vtilde_out->B == v_in->B && vtilde_out->A == v_in->A);
    assert(eps >= 0.0 && tol_abs > 0.0);

    /* -- (A) the input MUST be a delta-inclusion (bounded below / injective).
     * prop_delta_hominc (ii) needs a > 2*delta; a delta-inclusion has a >= 1-delta.
     * THE GUARD IS sigma_min of v's coordinate matrix (aic_dhom_v_sigma_min), the
     * TRUE Frobenius unit-ball inclusion infimum inf_{X!=0} ||v(X)||_F/||X||_F —
     * NOT the basis sweep min_i ||vE[i]|| (F1 fix: the basis sweep is BLIND to a v
     * bounded below on every basis element that COLLAPSES a general combination;
     * sigma_min = 0 iff v collapses a direction, so it SOUNDLY rejects a
     * non-inclusion). We require a_in = sigma_min >= 0.5 (a coarse fail-loud gate
     * on the double midpoint, like the projection-nontriviality gate; see
     * aic_dhom.h / paper/FINDINGS.md). This is the cor_improvement hypothesis
     * "B is delta_max-INCLUDED into A" (.tex:1318, :451-453). */
    arb_t a_in;
    arb_init(a_in);
    aic_dhom_v_sigma_min(a_in, v_in, prec);
    if (arf_cmp_d(arb_midref(a_in), 0.5) < 0) {
        fprintf(stderr, "aic_errreduce: input is NOT a delta-inclusion "
                "(coordinate-matrix sigma_min = %.6e < 0.5; v collapses a "
                "direction / is not bounded below / not injective). "
                "cor_improvement hypothesis violated.\n",
                arf_get_d(arb_midref(a_in), ARF_RND_NEAR));
        abort();
    }
    arb_clear(a_in);

    /* -- (B) vtilde = lem_approx(v_in): copy then run the Newton iteration. The
     * termination floor is AIC_ERRREDUCE_EPS_FLOOR*eps + tol_abs (NOT bare eps:
     * lem_approx reduces the defect to ~c_0*eps, a small multiple of eps, and
     * cannot beat the algebra's intrinsic O(eps) obstruction — driving toward
     * bare eps would stall and bounce, fail-loud in aic_dhom_approx). lem_approx
     * fails loud if the cap is hit without reaching the floor (.tex:1312).
     * unit_tol generous (the unit defect is only O(delta+eps) per
     * prop_delta_hominc (iii); dhom guards a BLOWUP, not the legitimate drift). */
    aic_errreduce_copy_v(vtilde_out, v_in);
    aic_dhom_approx_stats st;
    aic_dhom_approx(vtilde_out, /*eps_target*/ AIC_ERRREDUCE_EPS_FLOOR * eps,
                    /*tol_abs*/ tol_abs, /*unit_tol*/ 1.0, max_steps, prec, &st);

    /* -- (C) the four inclusion-defects of vtilde (prop_delta_hominc + lem_approx).
     * lower_gap = max(0, 1 - a) uses the SOUND a = sigma_min(coordinate matrix)
     * (F1 fix), NOT the basis-sweep min_i ||vE[i]|| — so the certificate enforces
     * the TRUE inclusion lower bound (a v that collapsed a direction would have
     * sigma_min ~ 0 -> lower_gap ~ 1 -> the certification guard (D) aborts). */
    arb_t mult, nub, nlb, udef;
    arb_init(mult);
    arb_init(nub);
    arb_init(nlb);
    arb_init(udef);
    aic_dhom_defect_sweep(mult, vtilde_out, prec);            /* mult. defect    */
    aic_dhom_prop_bounds(nub, NULL, udef, vtilde_out, prec);  /* ||v||, unit     */
    aic_dhom_v_sigma_min(nlb, vtilde_out, prec);             /* a = sigma_min   */

    aic_errreduce_defects local;
    aic_errreduce_defects *d = defects ? defects : &local;
    if (!defects) aic_errreduce_defects_init(d);
    arb_set(d->mult_def, mult);
    pos_part_sub1(d->norm_excess, nub, prec);   /* max(0, ||vtilde|| - 1) */
    pos_part_1sub(d->lower_gap, nlb, prec);      /* max(0, 1 - a)          */
    arb_set(d->unit_def, udef);

    arb_t maxdef;
    arb_init(maxdef);
    aic_errreduce_defects_max(maxdef, d);

    /* -- (D) CERTIFY vtilde is an O(eps)-inclusion: max-defect <= C0_CERT*eps_eff.
     * eps_eff = max(eps, tol_abs) (the eta=0 case has eps=0 but the defects floor
     * at the machine tol, so guard against a 0 ceiling). Fail loud (Rule 4). */
    double eps_eff = eps > tol_abs ? eps : tol_abs;
    arb_t ceil;
    arb_init(ceil);
    arb_set_d(ceil, AIC_ERRREDUCE_C0_CERT * eps_eff);
    if (!arb_le(maxdef, ceil)) {
        fprintf(stderr, "aic_errreduce: output is NOT certifiably an O(eps)-"
                "inclusion: max inclusion-defect %.6e > %.1f * max(eps,tol)=%.6e "
                "(eps=%.3e, %ld Newton steps).\n",
                arf_get_d(arb_midref(maxdef), ARF_RND_NEAR), AIC_ERRREDUCE_C0_CERT,
                AIC_ERRREDUCE_C0_CERT * eps_eff, eps, (long) st.steps);
        abort();
    }

    /* -- (E) the MEASURED c_0 = max-defect / eps (the empirical witness; the
     * analytic c_0 is aic-1bc, FINDINGS D2). eta=0 case (eps <= tol_abs): the
     * defects are ~machine-zero, report c_0 = 0 (avoid divide-by-zero). */
    if (eps <= tol_abs) {
        arb_zero(c0_out);
    } else {
        arb_set_d(ceil, eps);              /* reuse `ceil` as the denominator */
        arb_div(c0_out, maxdef, ceil, prec);
    }

    arb_clear(ceil);
    arb_clear(maxdef);
    if (!defects) aic_errreduce_defects_clear(&local);
    arb_clear(udef);
    arb_clear(nlb);
    arb_clear(nub);
    arb_clear(mult);
}

/* --------------------------------------------------------------- bijectivity */

int aic_errreduce_is_bijective(arb_t a_out, const aic_dhom_v *vtilde, slong prec)
{
    /* .tex:1318: bijective delta-inclusion -> bijective new inclusion. In finite
     * dim a delta-isomorphism is a bijective delta-inclusion; vtilde stays
     * injective (lower bound a > 0, prop_delta_hominc (ii) preserved by lem_approx)
     * and bijective iff dim_B == dim_A (a bounded-below linear map between equal
     * finite dims is a bijection). INJECTIVITY uses the SOUND a = sigma_min of the
     * coordinate matrix (F1 fix), NOT the basis-sweep min_i ||vE[i]||: a v that
     * collapses a direction (a kernel) is NOT injective and sigma_min ~ 0 catches
     * it, where the basis sweep would wrongly report bounded-below. */
    arb_t a;
    arb_init(a);
    aic_dhom_v_sigma_min(a, vtilde, prec);
    int injective = arf_cmp_d(arb_midref(a), 0.5) >= 0; /* bounded below */
    int dim_match = (vtilde->B->dim_B == vtilde->A->dim_A);
    if (a_out != NULL) arb_set(a_out, a);
    arb_clear(a);
    return injective && dim_match;
}
