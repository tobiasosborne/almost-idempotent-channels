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
 * Version single-sourcing (bead aic-95g.1, Inc4): CMake's project(AIC VERSION)
 * is now the ONE source of truth. The CMake build injects the real string via
 * a compile definition (-DAIC_VERSION_STRING="aic ${PROJECT_VERSION}", set on
 * the aic_objs OBJECT library in CMakeLists.txt), so the package/SONAME and
 * this symbol can never disagree. The in-file #define below is ONLY a fallback
 * so a non-CMake build (a bare `cc`) still compiles; it is intentionally a
 * loud "-dev" placeholder, never a real release number. The Julia ccall layer
 * (later) reads the version through this same symbol rather than duplicating a
 * constant.
 */
#include "aic/aic.h"

#ifndef AIC_VERSION_STRING
#define AIC_VERSION_STRING "aic 0.0.0-dev"   /* CMake injects the real version */
#endif

const char *aic_version(void)
{
    return AIC_VERSION_STRING;
}
