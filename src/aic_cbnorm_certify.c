/* aic_cbnorm_certify.c — TIGHT certified two-sided arb ball on the
 * eta-idempotence defect eta = ||Phi^2 - Phi||_cb = ||.||_diamond (bead aic-m24,
 * increment 3b): the PUBLIC dispatcher. The two feasibility-restoration recipes
 * (the bulk of the math) live in src/aic_cbnorm_certify_restore.c
 * (aic_cbnorm_int_lower / aic_cbnorm_int_upper); this file wires them together,
 * applies the dispatch and the fail-loud stop condition. See
 * include/aic_cbnorm.h for the API and docs/research/cbnorm_tight_certifier.md for THE
 * CONTRACT (the normalization table, the two recipes, the dispatch rules).
 *
 * STRATEGY (Law 1; .tex:347-354, the cb-norm = diamond-norm SDP). Two INDEPENDENT
 * MOSEK feasible points bracket eta by weak duality: the MAX-primal (X,P,Q) gives
 * a rigorous LOWER bound, the MIN-dual (Y0,Y1) a rigorous UPPER bound. Each is a
 * double-precision APPROXIMATE solution, so each is RESTORED to EXACT feasibility
 * in arb (convex combination toward the Slater center for the primal; eigenvalue
 * shift for the dual) before its objective is read off; only then is the bound
 * rigorous. The gap hi-lo is ~ the MOSEK duality gap + the arb radius, much
 * tighter than the always-valid eig-free 2n-wide bracket [||J||_F/n, 2||J||_F].
 *
 * DISPATCH (design doc §dispatch; Rule 4).
 *   - eta=0 regime (||J||_F < ~1e-9*n): defer to aic_cbnorm_eigfree_ball_choi
 *     (returns ~[0,0]); skip the feasible-point path.
 *   - restoration off (delta or eps > 0.1: MOSEK badly off / precision wall) OR
 *     bracket inverted (arb_lower(lo) > arb_upper(hi)): defer to the eig-free
 *     bracket as the rigorous fallback; if even THAT straddles, FAIL LOUD with
 *     eta/n/prec (the restoration recipe is wrong; orchestrator must re-derive).
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_cbnorm.h"
#include "aic_cbnorm_internal.h"
#include "aic/aic_mat.h"

void aic_cbnorm_certify(arb_t lo, arb_t hi, const acb_mat_t J,
                        const acb_mat_t X, const acb_mat_t P, const acb_mat_t Q,
                        const acb_mat_t Y0, const acb_mat_t Y1,
                        slong n, slong prec)
{
    slong N = n * n;
    assert(n >= 1 && "aic_cbnorm_certify: n >= 1 (Rule 4)");
    assert(acb_mat_nrows(J) == N && acb_mat_ncols(J) == N);
    assert(acb_mat_nrows(X) == N && acb_mat_ncols(X) == N);
    assert(acb_mat_nrows(P) == N && acb_mat_ncols(P) == N);
    assert(acb_mat_nrows(Q) == N && acb_mat_ncols(Q) == N);
    assert(acb_mat_nrows(Y0) == N && acb_mat_ncols(Y0) == N);
    assert(acb_mat_nrows(Y1) == N && acb_mat_ncols(Y1) == N);

    /* Dispatch: eta=0 regime (||J||_F < ~1e-9*n) -> eig-free ~[0,0]. */
    arb_t fro;
    arb_init(fro);
    aic_mat_frobenius_norm(fro, J, prec);
    arf_t froub, thr;
    arf_init(froub); arf_init(thr);
    arb_get_ubound_arf(froub, fro, prec);
    arf_set_d(thr, 1e-9 * (double) n);
    int eta_is_zero = (arf_cmp(froub, thr) < 0);
    arf_clear(thr); arf_clear(froub);
    arb_clear(fro);

    if (eta_is_zero) {
        aic_cbnorm_eigfree_ball_choi(lo, hi, J, n, prec);
        return;
    }

    arb_t delta, eps;
    arb_init(delta); arb_init(eps);
    aic_cbnorm_int_lower(lo, delta, J, X, P, Q, n, prec);
    aic_cbnorm_int_upper(hi, eps, J, Y0, Y1, n, prec);

    /* Restoration sanity + inversion check (design doc §dispatch). */
    arf_t dub, eub, big, lolb, hiub;
    arf_init(dub); arf_init(eub); arf_init(big);
    arf_init(lolb); arf_init(hiub);
    arb_get_ubound_arf(dub, delta, prec);
    arb_get_ubound_arf(eub, eps, prec);
    arf_set_d(big, 0.1);
    arb_get_lbound_arf(lolb, lo, prec);
    arb_get_ubound_arf(hiub, hi, prec);
    int restoration_off = (arf_cmp(dub, big) > 0) || (arf_cmp(eub, big) > 0);
    int inverted = (arf_cmp(lolb, hiub) > 0);   /* lo > hi: bracket inverted */

    if (restoration_off || inverted) {
        /* Rigorous fallback: the always-valid eig-free bracket. */
        arb_t flo, fhi;
        arb_init(flo); arb_init(fhi);
        aic_cbnorm_eigfree_ball_choi(flo, fhi, J, n, prec);
        arf_t flolb, fhiub;
        arf_init(flolb); arf_init(fhiub);
        arb_get_lbound_arf(flolb, flo, prec);
        arb_get_ubound_arf(fhiub, fhi, prec);
        if (arf_cmp(flolb, fhiub) > 0) {
            /* Even the eig-free bracket straddles -> escalate (Rule 4). */
            fprintf(stderr,
                    "aic_cbnorm_certify: bracket inverted AND eig-free fallback "
                    "straddles; n=%ld prec=%ld delta_ub=%.3e eps_ub=%.3e\n",
                    (long) n, (long) prec, arf_get_d(dub, ARF_RND_NEAR),
                    arf_get_d(eub, ARF_RND_NEAR));
            abort();
        }
        arb_set(lo, flo);
        arb_set(hi, fhi);
        arf_clear(fhiub); arf_clear(flolb);
        arb_clear(fhi); arb_clear(flo);
    }

    arf_clear(hiub); arf_clear(lolb);
    arf_clear(big); arf_clear(eub); arf_clear(dub);
    arb_clear(eps); arb_clear(delta);
}
