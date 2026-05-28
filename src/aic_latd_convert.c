/* aic_latd_convert.c — acb_mat_t <-> row-major double-complex array conversions
 * for the LAPACK double path (beads aic-w4o.5). See include/aic_latd.h for the
 * full storage/convention contract; the short version:
 *
 *   - Layout is ROW-MAJOR: arr[i*cols + j] = entry (i, j). This matches both the
 *     C-natural layout and acb_mat's (i,j) indexing, so neither direction needs
 *     a transpose — the conversions are a straight per-entry copy. The LAPACKE
 *     calls in the sibling files all pass LAPACK_ROW_MAJOR, consistent with this.
 *
 *   - acb_mat -> array takes the BALL MIDPOINT (this is the UNCERTIFIED path;
 *     the radius is intentionally discarded). We read the midpoint via
 *     arb_get_d(... ARF_RND_NEAR): arf_get_d rounds the exact arf midpoint to
 *     the nearest double, which is the correct double image of the ball centre.
 *
 *   - array -> acb_mat writes each double as a ZERO-RADIUS acb (acb_set_d_d puts
 *     the exact double value as the midpoint with radius 0): the double is taken
 *     as given, no spurious error ball is invented.
 *
 * Fail-loud (CLAUDE.md Rule 4): the array->matrix direction asserts the target
 * acb_mat is exactly r x c. The matrix->array direction writes r*c entries into
 * a raw caller buffer (size is the caller's contract, documented in the header).
 */
#include <assert.h>
#include <complex.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic_latd.h"

void aic_latd_from_acb_mat(double _Complex *out, const acb_mat_t A)
{
    slong r = acb_mat_nrows(A);
    slong c = acb_mat_ncols(A);
    assert(r >= 0 && c >= 0);

    for (slong i = 0; i < r; i++) {
        for (slong j = 0; j < c; j++) {
            const acb_struct *e = acb_mat_entry(A, i, j);
            /* midpoint of the real/imag arb balls, rounded to nearest double. */
            double re = arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR);
            double im = arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR);
            out[i * c + j] = re + im * _Complex_I;
        }
    }
}

void aic_latd_to_acb_mat(acb_mat_t out, const double _Complex *in,
                         slong r, slong c)
{
    assert(acb_mat_nrows(out) == r && "array->acb_mat: row mismatch");
    assert(acb_mat_ncols(out) == c && "array->acb_mat: col mismatch");

    for (slong i = 0; i < r; i++) {
        for (slong j = 0; j < c; j++) {
            double _Complex z = in[i * c + j];
            /* exact double -> zero-radius acb (midpoint = z, radius = 0). */
            acb_set_d_d(acb_mat_entry(out, i, j), creal(z), cimag(z));
        }
    }
}
