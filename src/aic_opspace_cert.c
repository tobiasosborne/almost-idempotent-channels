/* aic_opspace_cert.c — the O1 certification pipeline of the opspace module (bead
 * aic-zwo): aic_opspace_certify runs the operator-norm ampliation levels of
 * th_main_ext (approximate_algebras.tex:1447-1561, §10) over the existing
 * v: B -> A and fills aic_opspace_result. See include/aic_opspace.h for the
 * contract and the FINDINGS §C12 (operator-norm not Frobenius) / §D3 (Smith
 * truncation) / honesty (HOPM lower bound, SDP is O2) notes.
 *
 * WHAT THIS RUNS.
 *   forward Smith level N = v->n:  cb_forward = ||v_N||_op (HOPM, LOWER bound on
 *     ||v||_cb).
 *   inverse Smith level n_B = sum_l d_l: a_cb = a_{n_B} = 1/||v_{n_B}^{-1}||_op.
 *   doubling levels n = 1, 2, 4, ..., n_B: the lower-isometry a_n = 1/||v_n^{-1}||,
 *     checked for the prop_inc_ext doubling a_{2n} >= a_n/2 (.tex:1493-1503) and
 *     fed to the universality-canary flatness a_cb_flat = max a_n / min a_n.
 *
 * GUARDS (fail loud, Rule 4):
 *   - v bijective (dim_A == dim_B); else the inverse map is undefined.
 *   - prop_inc_ext precondition a_1 > 2*delta (.tex:1505): if violated, the
 *     induction that lifts a_1 to a_n for all n breaks — a genuine stop condition,
 *     NOT a truncation artifact. Abort with the measured a_1 and delta.
 *   - doubling a_{2n} >= a_n/2 at each consecutive level pair (.tex:1503): if
 *     violated, the prop_inc_ext mechanism itself broke. Abort.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/arb.h>

#include "aic_opspace.h"

static double dmid(const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

void aic_opspace_certify(aic_opspace_result *r, const aic_dhom_v *v,
                         double delta, slong prec)
{
    assert(r != NULL && v != NULL);
    slong dimA = v->A->dim_A, dimB = v->B->dim_B, N = v->n, nB = v->B->n_B;
    /* v bijective: dim_A == dim_B (so M_1 is square + invertible — the case
     * th_main_ext / prop_inc_ext assume; FINDINGS §D3, the a_n = 1/||v^{-1}_n||
     * identity). Fail loud otherwise (stop condition, not a silent skip). */
    if (dimA != dimB) {
        fprintf(stderr,
                "aic_opspace_certify: v not bijective (dim_A=%ld != dim_B=%ld); "
                "the extended-isomorphism a_n = 1/||v^{-1}_n||_op needs a square "
                "invertible coordinate matrix (th_main_ext bijective case)\n",
                (long) dimA, (long) dimB);
        abort();
    }

    r->v = v;
    arb_init(r->a_cb);
    arb_init(r->cb_forward);

    /* the doubling levels 1, 2, 4, ..., up to the inverse Smith truncation n_B. */
    assert(nB >= 1);
    slong levels[16];
    int nlev = 0;
    for (slong n = 1; n <= nB && nlev < 16; n *= 2) levels[nlev++] = n;
    assert(nlev >= 1);
    if (levels[nlev - 1] != nB && nlev < 16) levels[nlev++] = nB; /* ensure n_B */

    /* a_n = 1/||v_n^{-1}||_op at each level. */
    double an[16] = {0};
    arb_t s;
    arb_init(s);
    for (int li = 0; li < nlev; li++) {
        aic_opspace_inverse_stretch(s, v, levels[li], prec);
        double inv = dmid(s);
        an[li] = inv > 1e-300 ? 1.0 / inv : 0.0;
    }
    arb_clear(s);

    /* prop_inc_ext precondition: a_1 > 2*delta (.tex:1505). */
    double a1 = an[0];
    if (!(a1 > 2.0 * delta)) {
        fprintf(stderr,
                "aic_opspace_certify: prop_inc_ext precondition a_1 > 2*delta "
                "VIOLATED: a_1=%.6e delta=%.6e (2*delta=%.6e). The extended "
                "delta-isomorphism induction (.tex:1505) breaks — v is not a good "
                "enough inclusion (stop condition, FINDINGS §C12/§D3)\n",
                a1, delta, 2.0 * delta);
        abort();
    }

    /* prop_inc_ext doubling a_{2n} >= a_n/2 (.tex:1503), with a small numerical
     * slack for the double HOPM gate. */
    for (int li = 0; li + 1 < nlev; li++) {
        if (levels[li + 1] == 2 * levels[li]) {
            if (!(an[li + 1] >= an[li] / 2.0 - 1e-9)) {
                fprintf(stderr,
                        "aic_opspace_certify: prop_inc_ext doubling a_{2n} >= a_n/2 "
                        "VIOLATED at n=%ld: a_n=%.6e a_{2n}=%.6e (.tex:1503 mechanism "
                        "broke)\n",
                        (long) levels[li], an[li], an[li + 1]);
                abort();
            }
        }
    }

    /* a_cb = a_{n_B} (the all-n value by Smith); the universality canary flatness. */
    arb_set_d(r->a_cb, an[nlev - 1]);
    double amax = an[0], amin = an[0];
    for (int li = 1; li < nlev; li++) {
        if (an[li] > amax) amax = an[li];
        if (an[li] < amin) amin = an[li];
    }
    r->a_cb_flat = amin > 1e-12 ? amax / amin : 0.0;

    /* forward stretch at the forward Smith level N = v->n. */
    aic_opspace_forward_stretch(r->cb_forward, v, N, prec);

    r->N_max = nB > N ? nB : N;
    r->n_B = nB; /* the inverse Smith level (I2): factorize/D4 reads it off r */
}

void aic_opspace_result_clear(aic_opspace_result *r)
{
    if (r == NULL) return;
    arb_clear(r->a_cb);
    arb_clear(r->cb_forward);
    r->v = NULL;
}
