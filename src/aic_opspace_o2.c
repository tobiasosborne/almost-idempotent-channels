/* aic_opspace_o2.c — the O2 payoff pipeline of the opspace module (bead aic-pjr):
 * the CERTIFIED cb-norm UPPER bounds of th_main_ext (.tex:1538-1540). It ties the
 * O2.1 adjoint-Choi assemblers (aic_opspace_choi_vstar / _vinvstar) to the O2.4
 * rectangular Watrous-SDP UPPER certifier (aic_cbnorm_certify_rect_upper), turning
 * the committed MIN-dual feasible points into rigorous arb bounds on ||v||_cb and
 * ||v^{-1}||_cb, and asserts the HOPM(O1) <= SDP(O2) bracket. See
 * include/aic_opspace.h for the contract and docs/research/opspace_o2_design.md
 * §0.5 (the pinned convention) + §4.1 (the bracket cross-check, bead aic-0at).
 *
 * THE BRACKET (the load-bearing cross-check, aic-0at; design §4.1). O1's HOPM gives
 * a LOWER bound on the operator-norm max-stretch, which by Smith equals the cb-norm
 * (forward at level N, inverse at n_B):
 *   r->cb_forward = HOPM(||v_N||_op)        <= ||v||_cb     (O1 lower),
 *   1/r->a_cb     = HOPM(||v_{n_B}^{-1}||)  <= ||v^{-1}||_cb (O1 lower).
 * O2 here produces the matching certified UPPERS. The two MUST bracket:
 *   cb_forward <= cb_upper,   1/a_cb <= cbinv_upper.
 * If the HOPM lower EXCEEDS the SDP upper, one of the two routes is unsound — abort
 * (Rule 4), since that breaks weak duality. A small slack covers the arb radius and
 * the double-HOPM gate posture (the HOPM is a coarse double midpoint, aic_opspace.h
 * PRECISION POSTURE).
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_cbnorm.h"
#include "aic/aic_opspace.h"

/* rigorous lower endpoint of an arb ball as a double (FLOOR). */
static double arb_lb_d(const arb_t x)
{
    arf_t t;
    arf_init(t);
    arb_get_lbound_arf(t, x, 64);
    double v = arf_get_d(t, ARF_RND_FLOOR);
    arf_clear(t);
    return v;
}

/* rigorous upper endpoint of an arb ball as a double (CEIL). */
static double arb_ub_d(const arb_t x)
{
    arf_t t;
    arf_init(t);
    arb_get_ubound_arf(t, x, 64);
    double v = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
    return v;
}

void aic_opspace_certify_cb_upper(aic_opspace_result *r, const aic_dhom_v *v,
                                  const acb_mat_t Y0_fwd, const acb_mat_t Y1_fwd,
                                  const acb_mat_t Y0_inv, const acb_mat_t Y1_inv,
                                  slong prec)
{
    assert(r != NULL && v != NULL);
    slong N = v->n, nB = v->B->n_B;
    assert(v->A->dim_A == v->B->dim_B
           && "aic_opspace_certify_cb_upper: v not bijective (dim_A != dim_B)");

    /* FORWARD: ||v||_cb = ||v*||_⋄, J(v*) sized (d_maj = N, d_min = n_B). */
    acb_mat_t Jvs;
    acb_mat_init(Jvs, N * nB, N * nB);
    aic_opspace_choi_vstar(Jvs, v, prec);
    aic_cbnorm_certify_rect_upper(r->cb_upper, Jvs, Y0_fwd, Y1_fwd, N, nB, prec);
    acb_mat_clear(Jvs);

    /* INVERSE: ||v^{-1}||_cb = ||(v^{-1})*||_⋄, J((v^{-1})*) sized
     * (d_maj = n_B, d_min = N). */
    acb_mat_t Jvis;
    acb_mat_init(Jvis, nB * N, nB * N);
    aic_opspace_choi_vinvstar(Jvis, v, prec);
    aic_cbnorm_certify_rect_upper(r->cbinv_upper, Jvis, Y0_inv, Y1_inv, nB, N, prec);
    acb_mat_clear(Jvis);

    /* BRACKET GUARD (Rule 4; design §4.1, aic-0at): O1 HOPM lower <= O2 SDP upper.
     * Use the MOST CONSERVATIVE reading of each ball so a valid bracket can never
     * be false-aborted: take HOPM's SMALLEST value (its ball's lower endpoint) and
     * compare against the SDP upper's LARGEST value (its ball's upper endpoint).
     * Currently cb_forward/a_cb are zero-radius double-HOPM gate midpoints so
     * lb==ub, but the rigorous endpoints are correct and future-proof should O1
     * ever carry a radius. */
    const double slack = 1e-6;
    double hopm_fwd = arb_lb_d(r->cb_forward);     /* O1 lower on ||v||_cb (smallest) */
    double sdp_fwd = arb_ub_d(r->cb_upper);        /* O2 certified upper (largest)    */
    if (hopm_fwd > sdp_fwd + slack) {
        fprintf(stderr,
                "aic_opspace_certify_cb_upper: FORWARD bracket VIOLATED: O1 HOPM "
                "lower=%.10e > O2 SDP upper=%.10e (HOPM > SDP breaks weak duality — "
                "a soundness bug; design §4.1, aic-0at)\n",
                hopm_fwd, sdp_fwd);
        abort();
    }

    /* 1/a_cb is the O1 HOPM lower on ||v^{-1}||_cb (a_cb = 1/HOPM(||v^{-1}||)).
     * The SMALLEST value of this lower bound uses the LARGEST a_cb = the ball's
     * UPPER endpoint (1/x is decreasing), again the most conservative reading so a
     * valid bracket is never false-aborted. */
    double a_cb = arb_ub_d(r->a_cb);               /* largest a_cb -> smallest 1/a_cb */
    double hopm_inv = a_cb > 1e-300 ? 1.0 / a_cb : 0.0;
    double sdp_inv = arb_ub_d(r->cbinv_upper);
    if (hopm_inv > sdp_inv + slack) {
        fprintf(stderr,
                "aic_opspace_certify_cb_upper: INVERSE bracket VIOLATED: O1 HOPM "
                "lower=1/a_cb=%.10e > O2 SDP upper=%.10e (HOPM > SDP; aic-0at)\n",
                hopm_inv, sdp_inv);
        abort();
    }
}
