/* test_shim.c — round-trip cross-check for the flat-double ccall SHIM
 * (bead aic-m24, increment 2a). The shim aic_ucp_choi_diff_d marshals
 * J = Choi(Phi^2 - Phi) (Convention-A, n^2 x n^2; aic_ucp.h) entirely through
 * flat double arrays so the future Julia package can ccall it WITHOUT any
 * FLINT type (acb_mat_t / slong) crossing the language boundary. The bonus
 * aic_cbnorm_eigfree_d does the same for the eig-free certified bracket.
 *
 * Golden master: tests/fixtures_d24.inc.h (17 fixtures; the complex_qubit
 * fixture carries NONZERO imaginary Choi entries and is the conjugation guard).
 *
 * Assertion families:
 *   (RT-golden) the shim's double-out Choi == the fixture's stored golden
 *       choi_re/choi_im, max abs entry diff <= 1e-12 (validates the math AND
 *       the marshalling against the independent Julia/MOSEK golden Choi).
 *   (RT-arb)    the shim's double-out Choi == the pure-arb path (build Kraus as
 *       acb, aic_ucp_compose + aic_ucp_choi_diff, extract to double), <= 1e-12
 *       at prec=53 (the double-vs-arb@53 ladder rung 2; isolates the shim's
 *       OWN marshalling from the underlying arb cores).
 *   (BONUS)     aic_cbnorm_eigfree_d (lo,hi) brackets eta_ref and matches the
 *       C-level aic_cbnorm_eigfree_ball endpoints to ~1e-12.
 *
 * The complex_qubit fixture MUST be exercised (it is, in the all-fixtures loop)
 * — it guards the conjugation surviving the double round-trip; the project has a
 * documented Choi-conjugation bug class (aic_ucp.h:104-107). The mutation proof
 * (drop the conjugation / transpose the index in the shim) is caught HERE.
 */
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic_cbnorm.h"
#include "aic_ucp.h"
#include "aic_ucp_shim.h"
#include "aic_test.h"

#include "fixtures_d24.inc.h"

#define PREC 53

/* Build the channel Phi from a fixture's Kraus ops (aic_ucp layout, identical to
 * test_cbnorm.c / test_ucp_d24.c). */
static void build_kraus_from_fixture(aic_ucp_kraus *phi, const aic_d24_fixture *f)
{
    aic_ucp_kraus_init(phi, f->n, f->n, f->r);
    for (slong a = 0; a < f->r; a++)
        for (slong i = 0; i < f->n; i++)
            for (slong j = 0; j < f->n; j++) {
                slong off = a * f->n * f->n + i * f->n + j;
                acb_set_d_d(acb_mat_entry(phi->K[a], i, j),
                            f->kraus_re[off], f->kraus_im[off]);
            }
}

/* lower/upper endpoints of an arb ball as doubles (rigorous after rounding). */
static double arb_lo_d(const arb_t x)
{
    arf_t t;
    arf_init(t);
    arb_get_lbound_arf(t, x, 64);
    double v = arf_get_d(t, ARF_RND_FLOOR);
    arf_clear(t);
    return v;
}
static double arb_hi_d(const arb_t x)
{
    arf_t t;
    arf_init(t);
    arb_get_ubound_arf(t, x, 64);
    double v = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);
    return v;
}

int main(void)
{
    double worst_golden = 0.0; /* max |shim Choi - golden Choi| over corpus    */
    double worst_arb = 0.0;    /* max |shim Choi - pure-arb Choi| over corpus  */
    int nfix = 0;

    for (int fi = 0; fi < AIC_D24_NFIX; fi++) {
        const aic_d24_fixture *f = &aic_d24_fixtures[fi];
        int n = (int) f->n, r = (int) f->r, N = n * n;
        nfix++;

        /* ---- shim path: flat double in, flat double out ---- */
        double *cr = (double *) malloc((size_t) N * N * sizeof(double));
        double *ci = (double *) malloc((size_t) N * N * sizeof(double));
        AIC_CHECK(cr != NULL && ci != NULL);
        int rc = aic_ucp_choi_diff_d(cr, ci, f->kraus_re, f->kraus_im, n, r, PREC);
        AIC_CHECK_MSG(rc == 0, "fixture %s: shim returned %d", f->name, rc);

        /* ---- (RT-golden): shim Choi vs the fixture's stored golden Choi ---- */
        for (int k = 0; k < N * N; k++) {
            double dr = cr[k] - f->choi_re[k];
            double di = ci[k] - f->choi_im[k];
            double d = sqrt(dr * dr + di * di);
            if (d > worst_golden) worst_golden = d;
        }
        AIC_CHECK_MSG(worst_golden <= 1e-12,
                      "fixture %s: shim Choi vs golden max-diff %.3e > 1e-12",
                      f->name, worst_golden);

        /* ---- pure-arb path: build Kraus as acb, compose + choi_diff, extract ---- */
        aic_ucp_kraus phi, phi2;
        build_kraus_from_fixture(&phi, f);
        aic_ucp_compose(&phi2, &phi, &phi, PREC);
        acb_mat_t J;
        acb_mat_init(J, N, N);
        aic_ucp_choi_diff(J, &phi2, &phi, PREC);

        /* ---- (RT-arb): shim Choi vs the pure-arb Choi (rung 2) ---- */
        double fmax = 0.0;
        for (int p = 0; p < N; p++)
            for (int q = 0; q < N; q++) {
                double are = arf_get_d(arb_midref(acb_realref(acb_mat_entry(J, p, q))),
                                       ARF_RND_NEAR);
                double aim = arf_get_d(arb_midref(acb_imagref(acb_mat_entry(J, p, q))),
                                       ARF_RND_NEAR);
                double dr = cr[p * N + q] - are;
                double di = ci[p * N + q] - aim;
                double d = sqrt(dr * dr + di * di);
                if (d > fmax) fmax = d;
            }
        AIC_CHECK_MSG(fmax <= 1e-12,
                      "fixture %s: shim Choi vs pure-arb max-diff %.3e > 1e-12",
                      f->name, fmax);
        if (fmax > worst_arb) worst_arb = fmax;

        /* ---- (BONUS): aic_cbnorm_eigfree_d brackets eta_ref and matches the
         * C-level aic_cbnorm_eigfree_ball endpoints. ---- */
        {
            double lo = 0.0, hi = 0.0;
            int brc = aic_cbnorm_eigfree_d(&lo, &hi, f->kraus_re, f->kraus_im,
                                           n, r, PREC);
            AIC_CHECK_MSG(brc == 0, "fixture %s: cbnorm shim returned %d",
                          f->name, brc);
            /* rigorous bracket contains eta_ref (slack 1e-7 for SDP tol). */
            const double slack = 1e-7;
            AIC_CHECK_MSG(lo <= f->eta_ref + slack,
                          "fixture %s: cbnorm lo=%.6e > eta_ref=%.6e",
                          f->name, lo, f->eta_ref);
            AIC_CHECK_MSG(f->eta_ref - slack <= hi,
                          "fixture %s: cbnorm hi=%.6e < eta_ref=%.6e",
                          f->name, hi, f->eta_ref);
            AIC_CHECK_MSG(lo <= hi + 1e-15,
                          "fixture %s: cbnorm lo=%.6e > hi=%.6e (bad bracket)",
                          f->name, lo, hi);

            /* matches the C-level arb wrapper endpoints to ~1e-12. The shim's
             * doubles are FLOOR(lo)/CEIL(hi) of the same balls, so they must be
             * <= arb_lo_d(lo) and >= arb_hi_d(hi) respectively (rigorous after
             * rounding) and within 1e-12 of those endpoints. */
            arb_t loB, hiB;
            arb_init(loB);
            arb_init(hiB);
            aic_cbnorm_eigfree_ball(loB, hiB, &phi, PREC);
            double loB_lo = arb_lo_d(loB), hiB_hi = arb_hi_d(hiB);
            AIC_CHECK_MSG(lo <= loB_lo + 1e-12 && loB_lo - 1e-12 <= lo,
                          "fixture %s: cbnorm shim lo=%.17g vs C lo=%.17g",
                          f->name, lo, loB_lo);
            AIC_CHECK_MSG(hi >= hiB_hi - 1e-12 && hiB_hi + 1e-12 >= hi,
                          "fixture %s: cbnorm shim hi=%.17g vs C hi=%.17g",
                          f->name, hi, hiB_hi);
            arb_clear(hiB);
            arb_clear(loB);
        }

        acb_mat_clear(J);
        aic_ucp_kraus_clear(&phi2);
        aic_ucp_kraus_clear(&phi);
        free(ci);
        free(cr);
    }

    printf("fixtures: %d round-tripped\n", nfix);
    printf("worst |shim Choi - golden| = %.3e, worst |shim Choi - pure-arb| = %.3e\n",
           worst_golden, worst_arb);
    aic_test_report("test_shim");
    printf("OK test_shim\n");
    return 0;
}
