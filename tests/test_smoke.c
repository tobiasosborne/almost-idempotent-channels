/* test_smoke.c — end-to-end link check (beads aic-aox, "T-build").
 *
 * Not a math test: it proves the toolchain is sound before any algorithm
 * exists. Two assertions:
 *   (a) aic_version() returns a non-NULL, non-empty string (the src + header
 *       link works);
 *   (b) the FLINT/arb path is reachable and computes: build a 2x2 identity as
 *       an acb_mat_t, square it (I*I = I), and assert the exact result equals
 *       identity both structurally (acb_mat_eq) and entrywise (acb_is_one /
 *       acb_is_zero). acb_mat_one yields exact balls (zero radius) and the
 *       product of exact integer entries is exact, so acb_mat_eq — which
 *       requires *definite* equality of exact balls — is the right predicate;
 *       acb_mat_equal (structural) or acb_mat_overlaps (only enclosure) would
 *       be weaker checks. This is the cross-check ladder rung 0 (CLAUDE.md
 *       §"cross-check ladder"): internal sanity that the library answers.
 *
 * Per CLAUDE.md Rule 5 every assertion checks a value, not mere absence of a
 * crash. assert() is the fail-fast mechanism (Rule 4); compiled without
 * NDEBUG.
 */
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <flint/acb_mat.h>

#include "aic/aic.h"

int main(void)
{
    /* (a) version string is real. */
    const char *v = aic_version();
    assert(v != NULL);
    assert(strlen(v) > 0);

    /* (b) FLINT/arb computes: I * I == I for the 2x2 identity. */
    const slong prec = 53;
    acb_mat_t id, prod;
    acb_mat_init(id, 2, 2);
    acb_mat_init(prod, 2, 2);

    acb_mat_one(id); /* exact identity: 1 on the diagonal, 0 off it */
    acb_mat_mul(prod, id, id, prec);

    /* Structural + definite equality of exact balls. */
    assert(acb_mat_eq(prod, id));

    /* Entrywise, to exercise the scalar acb_* predicates too. */
    assert(acb_is_one(acb_mat_entry(prod, 0, 0)));
    assert(acb_is_one(acb_mat_entry(prod, 1, 1)));
    assert(acb_is_zero(acb_mat_entry(prod, 0, 1)));
    assert(acb_is_zero(acb_mat_entry(prod, 1, 0)));

    acb_mat_clear(prod);
    acb_mat_clear(id);

    printf("test_smoke: OK\n");
    return 0;
}
