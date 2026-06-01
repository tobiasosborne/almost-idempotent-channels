/* aic_channels_diag.c — diagonal / pinching public channel constructors (bead
 * aic-7hg): dephasing, depolarizing, conditional expectation. See
 * include/aic/aic_channels.h for the Heisenberg Kraus convention
 * Phi(X) = sum_a K_a^dag X K_a and the per-constructor properties.
 *
 * All three Kraus sets here are HERMITIAN (diagonal entries / projectors), so the
 * observable map equals its Schrodinger dual and the dag-order is moot — but the
 * code still stores them in the documented observable form.
 *
 *  - dephasing / cond_expectation: K_a are orthogonal projectors P_a with
 *    sum_a P_a = 1 and P_a P_b = delta_{ab} P_a, hence sum K_a^dag K_a = 1
 *    (unital) AND Phi^2 = Phi exactly (P_a X P_a then P_a(.)P_a is idempotent
 *    because the cross terms vanish): these are the eta=0 oracle channels.
 *  - depolarizing Phi_p(X) = (1-p)X + p tr(X)/d 1: built from the identity-weight
 *    Kraus and the trace-replace block. Unital; Phi_p^2 = Phi_{p(2-p)} (so p=1 is
 *    idempotent, the conditional expectation onto C*1; general p almost-idempotent).
 */
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic/aic_channels.h"
#include "aic/aic_ucp.h"
#include "aic_channels_internal.h"

void aic_channel_dephasing(aic_ucp_kraus *phi, slong d, slong prec)
{
    (void) prec;
    if (d < 1) { flint_printf("aic_channel_dephasing: d<1\n"); flint_abort(); }
    /* K_i = |e_i><e_i| (d Hermitian rank-1 projectors). matches make_dephasing. */
    aic_ucp_kraus_init(phi, d, d, d);
    for (slong i = 0; i < d; i++)
        acb_set_si(acb_mat_entry(phi->K[i], i, i), 1);
}

void aic_channel_cond_expectation(aic_ucp_kraus *phi, const slong *block_dims,
                                  slong nblocks, slong prec)
{
    (void) prec;
    if (block_dims == NULL || nblocks < 1) {
        flint_printf("aic_channel_cond_expectation: block_dims NULL or "
                     "nblocks<1 (Rule 4)\n");
        flint_abort();
    }
    slong d = 0;
    for (slong b = 0; b < nblocks; b++) {
        if (block_dims[b] < 1) {
            flint_printf("aic_channel_cond_expectation: block_dims[%wd]=%wd "
                         "< 1 (Rule 4)\n", b, block_dims[b]);
            flint_abort();
        }
        d += block_dims[b];
    }
    /* K_b = P_b = projector onto the contiguous index range of block b; the
     * ranges partition [0,d), so sum_b P_b = 1_d and P_a P_b = delta_{ab} P_a:
     * unital + exactly idempotent. Generalises make_block_cond_exp. */
    aic_ucp_kraus_init(phi, d, d, nblocks);
    slong off = 0;
    for (slong b = 0; b < nblocks; b++) {
        for (slong i = 0; i < block_dims[b]; i++)
            acb_set_si(acb_mat_entry(phi->K[b], off + i, off + i), 1);
        off += block_dims[b];
    }
}

void aic_channel_depolarizing(aic_ucp_kraus *phi, slong d, double p, slong prec)
{
    (void) prec;
    if (d < 1) { flint_printf("aic_channel_depolarizing: d<1\n"); flint_abort(); }
    double pv = p;
    aic_chan_validate_probs((double[]){1.0 - pv, pv}, 2, AIC_CHAN_TOL,
                            "aic_channel_depolarizing (p in [0,1])");

    /* Phi_p(X) = (1-p) X + p tr(X)/d 1_d. Kraus realisation:
     *   identity weight  : K_0 = sqrt(1-p) 1_d                 [ (1-p) X ]
     *   trace-replace    : K_{1+(i*d+j)} = sqrt(p/d) |e_i><e_j| (d^2 ops)
     *                      gives sum K^dag X K = (p/d) sum_{i,j} |e_j><e_i| X
     *                      |e_i><e_j| = (p/d) tr(X) 1_d.
     * Total r = 1 + d^2. Hermitian-superoperator self-dual; unital since
     * sum K^dag K = (1-p) 1 + (p/d) sum_{i,j}|e_j><e_j| = (1-p)1 + p 1 = 1. */
    slong r = 1 + d * d;
    aic_ucp_kraus_init(phi, d, d, r);

    double s_id = sqrt(1.0 - pv);
    for (slong i = 0; i < d; i++)
        acb_set_d(acb_mat_entry(phi->K[0], i, i), s_id);

    double s_tr = sqrt(pv / (double) d);
    slong a = 1;
    for (slong i = 0; i < d; i++)
        for (slong j = 0; j < d; j++) {
            acb_set_d(acb_mat_entry(phi->K[a], i, j), s_tr);
            a++;
        }
}
