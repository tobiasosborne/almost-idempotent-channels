/* aic.h — public C API for almost-idempotent-channels.
 *
 * This is Layer-0 scaffolding (beads aic-aox, "T-build"): the header exists so
 * later modules (funcalc, contraction, ucp, ... per MODULE_PLAN.md) have a
 * declared place to live, and so the build/test harness has a real symbol to
 * link against. No algorithm APIs are declared yet — those arrive with their
 * owning beads, each carrying its own paper/src/approximate_algebras.tex
 * provenance (CLAUDE.md Law 1). Inventing the full surface now would violate
 * Law 4 (audition, don't presume): the ABI sketch in MODULE_PLAN.md §"C ABI
 * surface" is a candidate, not a commitment.
 *
 * Storage-layout contract (placeholder).
 *   The arb number path represents complex matrices as FLINT `acb_mat_t`
 *   (row-major, each entry an `acb_t` = real+imag arb balls with certified
 *   error radii). The deferred double/fast path (LAPACK; open audition,
 *   beads aic-5ty) will fix a column-major `double _Complex` convention to
 *   match LAPACK's Fortran ABI. The exact layout, ownership, and the
 *   precision-argument convention are pinned here once the first algorithm
 *   module lands and the convention is exercised by a test — not before.
 */
#ifndef AIC_H
#define AIC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Returns a static, non-NULL, non-empty version string ("aic <major.minor.patch>").
 * Lifetime is the program's; the caller must not free it. */
const char *aic_version(void);

#ifdef __cplusplus
}
#endif

#endif /* AIC_H */
