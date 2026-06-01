/* aic_channels_internal.h — shared input-validation helpers for the public
 * channel constructors (bead aic-7hg). NOT installed (internal). The validators
 * implement the CLAUDE.md Rule 4 (fail-loud) contract documented in
 * include/aic/aic_channels.h: a malformed channel request ABORTS with a concrete
 * message rather than silently producing a non-channel. They check the
 * double-precision midpoints (caller-error guard, not a certified bound).
 */
#ifndef AIC_CHANNELS_INTERNAL_H
#define AIC_CHANNELS_INTERNAL_H

#include <flint/acb_mat.h>

/* Default validation tolerance (probabilities sum / unitarity). Generous for
 * double round-off in the caller's inputs, tight enough to catch real errors. */
#define AIC_CHAN_TOL 1e-10

/* Abort unless `probs[0..m)` are each >= -tol and sum to 1 within tol. `who`
 * names the calling constructor for the message. */
void aic_chan_validate_probs(const double *probs, slong m, double tol,
                             const char *who);

/* Abort unless U (d x d) satisfies ||U^dag U - 1_d||_op <= tol. `who`/`idx`
 * identify the offending unitary in the message. */
void aic_chan_validate_unitary(const acb_mat_t U, slong d, double tol,
                               const char *who, slong idx, slong prec);

#endif /* AIC_CHANNELS_INTERNAL_H */
