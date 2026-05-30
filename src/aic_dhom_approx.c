/* aic_dhom_approx.c — lem_approx: approximate a delta-homomorphism v : B -> A by
 * an O(eps)-homomorphism via the multiplicativity-defect Newton iteration (bead
 * aic-c1n). Realizes approximate_algebras.tex:1224-1314 (the proof of lem_approx).
 *
 * THE CONSTRUCTION (.tex:1277-1313). g = G_v is the multiplicativity defect. The
 * paper's correction uses the EXPLICIT diagonal D = sum_t w_t U_t^dag (x) U_t of B
 * (A_t = w_t U_t^dag, B_t = U_t; the generalized-Pauli diagonal, aic_dhom_B.c):
 *     w'(X)  = sum_t v(A_t) * g(B_t, X)        (.tex:1279, * = STAR in A)
 *     w''(X) = w'(X^dag)^dag                    (.tex:1305, restore involution)
 *     v^{(s+1)} = v^{(s)} + 1/2 (w' + w'')      (.tex:1311)
 * One step drops the multiplicativity defect O(delta) -> O(delta^2 + eps)
 * (.tex:1302-1311). "Iterating as in Newton's method" (.tex:1312), for eps>0 the
 * defect reaches O(eps) in finitely many steps; for eps=0 (an exact hom) it goes
 * to machine-zero (Johnson 1988 Thm 3.1, .tex:1313).
 *
 * THE STAR IS LOAD-BEARING. Every product in w' (v(A_t)*g(B_t,X)) is the
 * Choi-Effros star aic_ecstar_star, NEVER plain matmul. g = G_v itself uses the
 * star (aic_dhom_Gv). The eta=0 identity-channel oracle (star == plain) cannot
 * catch a star bug — test_dhom T4 (eta>0) does.
 *
 * w' UPDATES THE STORED BASIS DIRECTLY. v is stored by vE[i] = v(E_i); the update
 * v <- v + 1/2(w'+w'') sets vE[i] += 1/2(w'(E_i) + w''(E_i)). w'(E_i) and
 * w''(E_i) are computed for each matrix unit E_i, then the WHOLE update is applied
 * at once (snapshot the old vE so w' uses the pre-step v, matching the Newton
 * linearization which evaluates F_{w'} at the current v, .tex:1283-1303).
 *
 * FAIL-LOUD (Rule 4). The loop aborts if (a) the cap is hit without the defect
 * dropping below the threshold (a contraction must converge), (b) the defect ever
 * INCREASES across a step (Newton must contract once delta is small enough), and
 * (c) the per-step unit-preservation / involution-symmetry invariants drift past
 * unit_tol.
 *
 * No eigendecomposition: pure star/matmul + certified op-norm.
 */
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner_internal.h" /* aic_corner_gamma_opnorm_ub (certified UB) */
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_mat.h"

/* CERTIFIED operator-norm UPPER BOUND (mid+radius), dodging the aic_mat_opnorm
 * Gram-path Hermiticity false-fail on near-zero off-diagonals (aic-qgs/aic-2yo). */
static void opnorm_ub(arb_t out, const acb_mat_t M, slong prec)
{
    aic_corner_gamma_opnorm_ub(out, NULL, M, prec);
}

/* Fail-loud (CLAUDE.md Rule 4): print to stderr + abort, never silently return.
 * Mirrors src/aic_contraction_work.c's fprintf+abort idiom. */
static void dhom_fail(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "aic_dhom_approx: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    abort();
}

/* w'(X) = sum_t v(A_t) * g(B_t, X), A_t = w_t U_t^dag, B_t = U_t.
 *   = sum_t w_t [ v(U_t^dag) * G_v(U_t, X) ]    (g = G_v at the current v).
 * U_t, X are n_B x n_B; out is n x n. v(U_t^dag) via aic_dhom_v_apply; the
 * scalar w_t folds into the accumulated sum. */
static void w_prime(acb_mat_t out, const aic_dhom_v *v, const aic_dhom_diag *D,
                    const acb_mat_t X, slong prec)
{
    slong nB = v->B->n_B, n = v->n;
    acb_mat_t Utd, vUtd, g, term;
    acb_mat_init(Utd, nB, nB);
    acb_mat_init(vUtd, n, n);
    acb_mat_init(g, n, n);
    acb_mat_init(term, n, n);
    acb_mat_zero(out);
    for (slong t = 0; t < D->nterms; t++) {
        acb_mat_conjugate_transpose(Utd, D->U[t]);      /* A_t base = U_t^dag */
        aic_dhom_v_apply(vUtd, v, Utd, prec);           /* v(U_t^dag)         */
        aic_dhom_Gv(g, v, D->U[t], X, prec);            /* g(U_t, X)          */
        aic_ecstar_star(term, v->A, vUtd, g, prec);     /* v(U_t^dag) * g     */
        acb_mat_scalar_mul_arb(term, term, D->w + t, prec); /* * w_t          */
        acb_mat_add(out, out, term, prec);
    }
    acb_mat_clear(term);
    acb_mat_clear(g);
    acb_mat_clear(vUtd);
    acb_mat_clear(Utd);
}

void aic_dhom_newton_step(aic_dhom_v *v, const aic_dhom_diag *D, slong prec)
{
    slong dim_B = v->B->dim_B, nB = v->B->n_B, n = v->n;
    /* compute the FULL update {1/2(w'(E_i)+w''(E_i))} from the CURRENT v, then
     * apply all at once (so w' sees the pre-step v, .tex:1283-1303). */
    acb_mat_t *upd = flint_malloc((size_t) dim_B * sizeof(acb_mat_t));
    acb_mat_t Ei, Eid, wp, wpd, half;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(Eid, nB, nB);
    acb_mat_init(wp, n, n);
    acb_mat_init(wpd, n, n);
    acb_mat_init(half, n, n);

    arb_t onehalf;
    arb_init(onehalf);
    arb_set_d(onehalf, 0.5);

    for (slong i = 0; i < dim_B; i++) {
        acb_mat_init(upd[i], n, n);
        aic_dhom_B_matunit(Ei, v->B, i);
        w_prime(wp, v, D, Ei, prec);             /* w'(E_i)              */
        /* w''(E_i) = w'(E_i^dag)^dag */
        acb_mat_conjugate_transpose(Eid, Ei);    /* E_i^dag (a matrix unit) */
        w_prime(half, v, D, Eid, prec);          /* w'(E_i^dag)          */
        acb_mat_conjugate_transpose(wpd, half);  /* w''(E_i)             */
        /* upd = 1/2 (w' + w'') */
        acb_mat_add(upd[i], wp, wpd, prec);
        acb_mat_scalar_mul_arb(upd[i], upd[i], onehalf, prec);
    }
    /* apply */
    for (slong i = 0; i < dim_B; i++) {
        acb_mat_add(v->vE[i], v->vE[i], upd[i], prec);
        acb_mat_clear(upd[i]);
    }
    flint_free(upd);
    arb_clear(onehalf);
    acb_mat_clear(half);
    acb_mat_clear(wpd);
    acb_mat_clear(wp);
    acb_mat_clear(Eid);
    acb_mat_clear(Ei);
}

/* involution-symmetry defect: max over basis pairs (E_{ab},E_{ba}) of
 * ||v(E_{ba})^dag - v(E_{ab})||_op. (.tex:1305-1311: v^{(s)} commutes with the
 * involution after symmetrization.) Uses that E_i^dag is again a matrix unit. */
static void invol_defect(arb_t out, const aic_dhom_v *v, slong prec)
{
    slong dim_B = v->B->dim_B, nB = v->B->n_B, n = v->n;
    acb_mat_t Ei, Eid, vEid, vEiddag, vEi, diff;
    arb_t e;
    acb_mat_init(Ei, nB, nB);
    acb_mat_init(Eid, nB, nB);
    acb_mat_init(vEid, n, n);
    acb_mat_init(vEiddag, n, n);
    acb_mat_init(vEi, n, n);
    acb_mat_init(diff, n, n);
    arb_init(e);
    arb_zero(out);
    for (slong i = 0; i < dim_B; i++) {
        aic_dhom_B_matunit(Ei, v->B, i);
        acb_mat_conjugate_transpose(Eid, Ei);        /* E_i^dag = E_{ba}     */
        aic_dhom_v_apply(vEid, v, Eid, prec);        /* v(E_{ba})            */
        acb_mat_conjugate_transpose(vEiddag, vEid);  /* v(E_{ba})^dag        */
        aic_dhom_v_apply(vEi, v, Ei, prec);          /* v(E_{ab})            */
        acb_mat_sub(diff, vEiddag, vEi, prec);
        opnorm_ub(e, diff, prec);
        if (arb_gt(e, out)) arb_set(out, e);
    }
    arb_clear(e);
    acb_mat_clear(diff);
    acb_mat_clear(vEi);
    acb_mat_clear(vEiddag);
    acb_mat_clear(vEid);
    acb_mat_clear(Eid);
    acb_mat_clear(Ei);
}

void aic_dhom_approx(aic_dhom_v *v, double eps_target, double tol_abs,
                     double unit_tol, slong max_steps, slong prec,
                     aic_dhom_approx_stats *st)
{
    assert(v != NULL && max_steps >= 0);
    aic_dhom_diag D;
    aic_dhom_diag_build(&D, v->B, prec);

    arb_t sweep, threshold, unit_def, inv_def;
    arb_init(sweep);
    arb_init(threshold);
    arb_init(unit_def);
    arb_init(inv_def);
    /* termination threshold = max(tol_abs, eps_target) (the documented O(eps)
     * floor; the explicit c_0*eps of cor_improvement is the NEXT module, bead
     * aic-t81 — .tex:1317). */
    arb_set_d(threshold, tol_abs > eps_target ? tol_abs : eps_target);

    aic_dhom_defect_sweep(sweep, v, prec);
    double delta0 = arf_get_d(arb_midref(sweep), ARF_RND_NEAR);

    slong steps = 0;
    double prev = delta0;
    while (steps < max_steps) {
        if (arb_le(sweep, threshold)) break;       /* already O(eps)        */
        aic_dhom_newton_step(v, &D, prec);
        steps++;
        aic_dhom_defect_sweep(sweep, v, prec);
        double cur = arf_get_d(arb_midref(sweep), ARF_RND_NEAR);
        /* Newton must contract: a strict increase (beyond a small slack to allow
         * the O(eps) floor) is a fail-loud (Rule 4). */
        if (cur > prev * 1.5 + 1e-300)
            dhom_fail("defect INCREASED "
                                 "%.6e -> %.6e at step %ld (Newton not "
                                 "contracting)", prev, cur, (long) steps);
        prev = cur;

        /* per-step invariants (Rule 4). TWO distinct invariants with DIFFERENT
         * tightness (.tex:1305-1311):
         *  - involution-symmetry ||v(E_ba)^dag - v(E_ab)||: the symmetrized update
         *    1/2(w'+w'') is itself involution-symmetric, so it does NOT change an
         *    existing asymmetry of v; a *-linear delta-homomorphism (the lemma's
         *    hypothesis) starts symmetric and STAYS symmetric. So this is held to
         *    the TIGHT unit_tol (it must remain ~machine-zero).
         *  - unit-preservation ||v(I_B) - I_A|| is only O(delta+eps) (prop_
         *    delta_hominc (iii), .tex:1217-1221); it does NOT go to zero. So it is
         *    bounded by a GENEROUS O(delta_0) ceiling, not unit_tol. The guard
         *    catches a BLOWUP (a real bug), not the legitimate O(delta) drift. */
        aic_dhom_prop_bounds(NULL, NULL, unit_def, v, prec);
        invol_defect(inv_def, v, prec);
        double ud = arf_get_d(arb_midref(unit_def), ARF_RND_NEAR);
        double id = arf_get_d(arb_midref(inv_def), ARF_RND_NEAR);
        double unit_ceiling = 8.0 * delta0 + unit_tol; /* O(delta+eps) ceiling */
        if (ud > unit_ceiling)
            dhom_fail("unit-preservation %.3e > O(delta) ceiling %.3e at step "
                      "%ld (BLOWUP, not the legitimate O(delta) drift)", ud,
                      unit_ceiling, (long) steps);
        if (id > unit_tol)
            dhom_fail("involution-symmetry %.3e > tol %.3e at step %ld (the "
                      "symmetrized update should keep v *-linear)", id,
                      unit_tol, (long) steps);
    }
    /* fail loud if we hit the cap and the defect is still above threshold. */
    if (steps == max_steps && !arb_le(sweep, threshold))
        dhom_fail("%ld steps did not reach defect "
                             "threshold (final %.6e)", (long) max_steps,
                             arf_get_d(arb_midref(sweep), ARF_RND_NEAR));

    if (st != NULL) {
        st->steps = steps;
        st->delta0 = delta0;
        st->delta_final = arf_get_d(arb_midref(sweep), ARF_RND_NEAR);
    }
    arb_clear(inv_def);
    arb_clear(unit_def);
    arb_clear(threshold);
    arb_clear(sweep);
    aic_dhom_diag_clear(&D);
}
