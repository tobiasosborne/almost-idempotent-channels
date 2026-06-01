/* aic_channels_validate.c — fail-loud input validators for the public channel
 * constructors (bead aic-7hg). See aic_channels_internal.h. These guard the
 * CLAUDE.md Rule 4 contract: a caller that passes negative/unnormalised
 * probabilities or a non-unitary U gets a loud flint_abort with a concrete
 * message, never a silently malformed channel.
 *
 * The probability check works on the raw doubles. The unitarity check forms
 * U^dag U - 1 in arb and certifies ||.||_op <= tol via the existing
 * aic_mat_opnorm; a U that is unitary to round-off (e.g. exp(iH) from arb)
 * passes, a genuinely non-unitary U aborts.
 */
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic_channels_internal.h"

void aic_chan_validate_probs(const double *probs, slong m, double tol,
                             const char *who)
{
    if (probs == NULL || m < 1) {
        flint_printf("%s: probabilities array NULL or m<1 (Rule 4: fail loud)\n",
                     who);
        flint_abort();
    }
    double sum = 0.0;
    for (slong j = 0; j < m; j++) {
        if (probs[j] < -tol) {
            flint_printf("%s: probs[%wd] = %.17g is negative (< -%g); "
                         "a UCP mixture needs nonnegative weights (Rule 4)\n",
                         who, j, probs[j], tol);
            flint_abort();
        }
        sum += probs[j];
    }
    double d = sum - 1.0;
    if (d < 0) d = -d;
    if (d > tol) {
        flint_printf("%s: probabilities sum to %.17g, not 1 (|sum-1| = %g > %g); "
                     "a unital mixture needs sum_j p_j = 1 (Rule 4: fail loud)\n",
                     who, sum, d, tol);
        flint_abort();
    }
}

void aic_chan_validate_unitary(const acb_mat_t U, slong d, double tol,
                               const char *who, slong idx, slong prec)
{
    if (acb_mat_nrows(U) != d || acb_mat_ncols(U) != d) {
        flint_printf("%s: U[%wd] is %wd x %wd, expected %wd x %wd (Rule 4)\n",
                     who, idx, acb_mat_nrows(U), acb_mat_ncols(U), d, d);
        flint_abort();
    }
    acb_mat_t Ud, G, one;
    arb_t def, tolb;
    acb_mat_init(Ud, d, d);
    acb_mat_init(G, d, d);
    acb_mat_init(one, d, d);
    arb_init(def);
    arb_init(tolb);

    acb_mat_conjugate_transpose(Ud, U);
    acb_mat_mul(G, Ud, U, prec);          /* U^dag U          */
    acb_mat_one(one);
    acb_mat_sub(G, G, one, prec);         /* U^dag U - 1      */
    aic_mat_opnorm(def, G, prec);
    arb_set_d(tolb, tol);

    if (!arb_le(def, tolb)) {
        flint_printf("%s: U[%wd] is not unitary (||U^dag U - 1||_op exceeds %g); "
                     "a mixed-unitary channel needs unitary components (Rule 4)\n",
                     who, idx, tol);
        flint_abort();
    }

    arb_clear(tolb);
    arb_clear(def);
    acb_mat_clear(one);
    acb_mat_clear(G);
    acb_mat_clear(Ud);
}
