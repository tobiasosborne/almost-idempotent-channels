/* aic_mat_eigvec_resid.c — the SELF-CERTIFYING invariant-subspace residual on the
 * ORIGINAL H (bead aic-4td increment 2 step C2, inc-2 hostile-review finding-1 fix;
 * design docs/research/eigvec_certified_design.md §1.6, FINDINGS §D5n). Split out of
 * aic_mat_eigvec_seed.c (Rule 10, ~200 LOC/file).
 *
 * WHY THIS EXISTS (the finding-1 fix). aic_mat_eig_hermitian_subspaces certifies
 * the Rump enclosure on the DENSIFIED A' = U H U^dag, then back-maps X_c = U^dag X'_c
 * and claims X_c is a certified invariant subspace of the ORIGINAL H (design §1.6(ii):
 * A X_c = U^dag A' X'_c = U^dag X'_c J_c = X_c J_c). Before this fix the residual on H
 * was recomputed ONLY in test S1, leaving an O(||H|| * conjugation-error) gap
 * UNcertified in production. design §1.6 explicitly names ||H X_c - X_c J_c|| as THE
 * certificate; this routine recomputes it on H per cluster and FAILS LOUD (Rule 4) if
 * it is not small — so the production routine is SELF-CERTIFYING on H, not merely on A'.
 *
 * THE TOLERANCE (honest, prec-appropriate — measured, not guessed). The residual is
 * tied to the ENCLOSURE'S OWN certified ball radii, NOT to 2^-(prec/2). MEASURED
 * (probe, prec in {53,128,256}): the densify+Rump residual hits a CONDITIONING floor
 * (~1e-31 at prec=128, ~1e-27 at the n=7 block) that is PREC-INDEPENDENT above 128 —
 * so a 2^-(prec/2) tolerance would FALSE-FAIL at prec>=192 (the residual then exceeds
 * the tol). The residual genuinely tracks the enclosure radii: resid ~ (a few * n) *
 * maxrad(X_c) at every prec (ratio 9.5 at n=5 ... 25.8 at n=12). So
 *   tol = 64 * n * (1 + ||H||_F) * max(maxrad(X_c) + maxrad(J_c), 2^-(prec/2))
 * PASSES the legitimate residual with margin ~1e10..1e14 at prec 53/128/256, while
 * CATCHING a genuine subspace-invariance failure: an O(1) corruption of a column of
 * X_c gives residual ~||H|| (0.8, 4, 7.2 measured) >> tol (~1e-16..1e-25). The
 * 2^-(prec/2) term is a backstop for an EXACT (zero-radius) input; the 64 * n * (1 +
 * ||H||_F) factor absorbs the O(1) Frobenius/matmul constants over the n*k entries
 * (cf. the funcalc a-posteriori certificate aic_funcalc_int_certify_sign, which uses
 * the same max(prec-floor, in_rad) * O(sqrt(n)) shape).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic_mat_internal.h"

/* Largest entry-radius (real or imag) over an acb_mat — the intrinsic certified
 * error scale of an interval matrix. */
static double aic_mat_int_maxrad(const acb_mat_t M)
{
    double mx = 0.0;
    slong r = acb_mat_nrows(M), c = acb_mat_ncols(M);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++) {
            const acb_struct *e = acb_mat_entry(M, i, j);
            double a = mag_get_d(arb_radref(acb_realref(e)));
            double b = mag_get_d(arb_radref(acb_imagref(e)));
            if (a > mx) mx = a;
            if (b > mx) mx = b;
        }
    return mx;
}

void aic_mat_int_assert_subspace_residual(const acb_mat_t H,
                                          const aic_mat_eigcluster *out,
                                          slong n, slong cidx, slong prec)
{
    slong k = out->k;
    acb_mat_t HX, XJ, R;
    acb_mat_init(HX, n, k);
    acb_mat_init(XJ, n, k);
    acb_mat_init(R, n, k);
    acb_mat_mul(HX, H, out->X, prec);           /* H X_c */
    acb_mat_mul(XJ, out->X, out->J, prec);      /* X_c J_c (the STORED Rump J_c) */
    acb_mat_sub(R, HX, XJ, prec);

    arb_t resid, hf;
    arb_init(resid);
    arb_init(hf);
    acb_mat_frobenius_norm(resid, R, prec);     /* certified ||H X_c - X_c J_c||_F */
    aic_mat_frobenius_norm(hf, H, prec);
    double residub = arf_get_d(arb_midref(resid), ARF_RND_UP) +
                     mag_get_d(arb_radref(resid));
    double hfd = arf_get_d(arb_midref(hf), ARF_RND_UP) + mag_get_d(arb_radref(hf));

    double rad = aic_mat_int_maxrad(out->X) + aic_mat_int_maxrad(out->J);
    double floorr = ldexp(1.0, -(int) (prec / 2));   /* 2^-(prec/2) backstop */
    double radscale = (rad > floorr) ? rad : floorr;
    double tol = 64.0 * (double) n * (1.0 + hfd) * radscale;

    if (!(residub < tol)) {
        char *r_s = arb_get_str(resid, 8, 0);
        fprintf(stderr,
                "aic_mat_eig_hermitian_subspaces: cluster %ld (k=%ld) NOT a "
                "certified invariant subspace of the ORIGINAL H: ||H X_c - X_c J_c||"
                "_F=%s (ub %.3e) not < tol=%.3e at prec=%ld — the back-mapped X_c "
                "fails the self-certifying residual (raise prec or the input is not "
                "Hermitian-decomposable here); design §1.6, FINDINGS §D5n, bead "
                "aic-4td\n",
                (long) cidx, (long) k, r_s ? r_s : "?", residub, tol, (long) prec);
        flint_free(r_s);
        abort();
    }
    arb_clear(hf);
    arb_clear(resid);
    acb_mat_clear(R);
    acb_mat_clear(XJ);
    acb_mat_clear(HX);
}
