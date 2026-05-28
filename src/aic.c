/* aic.c — Layer-0 implementation (beads aic-aox, "T-build").
 *
 * The only routine here is aic_version(): a single real, linkable symbol whose
 * job is to prove the public header, the source build, and the test harness are
 * wired together end to end before any math lands. It carries no paper
 * provenance because it realizes no result from
 * paper/src/approximate_algebras.tex — when an algorithm module is added it
 * gets its own file (≤200 LOC, CLAUDE.md Rule 10) with its own .tex citation
 * (Law 1).
 *
 * The version is hand-maintained here as the single source of truth; the Julia
 * ccall layer (later) reads it through this same symbol rather than duplicating
 * a constant.
 */
#include "aic.h"

#define AIC_VERSION_STRING "aic 0.0.1"

const char *aic_version(void)
{
    return AIC_VERSION_STRING;
}
