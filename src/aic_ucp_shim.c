/* aic_ucp_shim.c — flat-double ccall SHIM over the arb UCP / cb-norm cores
 * (bead aic-m24, increment 2a). See include/aic_ucp_shim.h for the ABI/layout
 * contract and WHY the shim exists (no FLINT type may cross the Julia ccall
 * boundary). This file is PURE MARSHALLING: it converts flat double arrays to
 * the internal aic_ucp_kraus / acb_mat, runs the EXACT cores the C tests
 * exercise, and writes the result back to flat double arrays. No new math.
 *
 * PROVENANCE (Law 1). The Choi convention and the Heisenberg composition the
 * shim feeds are derived verbatim in include/aic_ucp.h:
 *   - Convention-A Choi (aic_ucp.h:23-42):
 *       C[i*n+a, j*n+b] = sum_x conj(K_x[i,a]) * K_x[j,b]
 *     (the K factor LEFT/MAJOR, the conjugation on the FIRST (i,a) factor). The
 *     shim must preserve this conjugation through the double round-trip — the
 *     documented Choi-conjugation bug class (aic_ucp.h:104-107) is guarded by
 *     the complex_qubit fixture in tests/test_shim.c.
 *   - Heisenberg composition (aic_ucp.h:149-165, aic_ucp_compose.c):
 *       Kraus(Phi o Psi) = { L_b K_a }; Phi^2 = compose(Phi,Phi).
 *   - eta-idempotence defect Lambda = Phi^2 - Phi, Choi by linearity
 *     (aic_ucp_choi_diff; .tex:347-354). The eig-free bracket
 *     [||J||_F/n, 2*||J||_F] is aic_cbnorm.c (ALGORITHM.md cbnorm section).
 *
 * Number path: the arb CERTIFIED path internally (acb_mat, explicit prec); only
 * the BOUNDARY is double. The double-vs-arb@53 cross-check (CLAUDE.md ladder
 * rung 2) and the double-vs-golden check are driven from tests/test_shim.c.
 *
 * Fail-loud (Rule 4): bad shape (n<1, r<1, null array) aborts at the call site
 * rather than producing silently wrong output. <=200 LOC (CLAUDE.md Rule 10).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic_cbnorm.h"
#include "aic_ucp.h"
#include "aic_ucp_shim.h"

/* Build the internal aic_ucp_kraus self-map Phi (dim_K == dim_H == n) from the
 * flat double Kraus arrays, layout kraus_*[a*n*n + i*n + j] (aic_ucp.h /
 * fixtures_d24.inc.h). `phi` is aic_ucp_kraus_init'd here; caller clears it. */
static void kraus_from_flat(aic_ucp_kraus *phi, const double *kraus_re,
                            const double *kraus_im, int n, int r)
{
    aic_ucp_kraus_init(phi, n, n, r);
    for (int a = 0; a < r; a++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                int off = a * n * n + i * n + j;
                acb_set_d_d(acb_mat_entry(phi->K[a], i, j),
                            kraus_re[off], kraus_im[off]);
            }
}

/* aic_ucp.h:167-174 — C = Choi(phi1) - Choi(phi2), Convention A; here
 * J = Choi(Phi^2) - Choi(Phi) = Choi(Lambda), Lambda = Phi^2 - Phi (.tex:347).
 * Flat-double ABI: builds Phi, Phi^2 = compose(Phi,Phi), J via choi_diff, then
 * extracts the n^2 x n^2 J back to choi_re/choi_im[p*N + q]. */
int aic_ucp_choi_diff_d(double *choi_re, double *choi_im,
                        const double *kraus_re, const double *kraus_im,
                        int n, int r, int prec)
{
    assert(n >= 1 && "aic_ucp_choi_diff_d: n >= 1 required (Rule 4)");
    assert(r >= 1 && "aic_ucp_choi_diff_d: r >= 1 required (Rule 4)");
    assert(prec >= 2 && "aic_ucp_choi_diff_d: prec >= 2 required (Rule 4)");
    assert(choi_re != NULL && choi_im != NULL &&
           "aic_ucp_choi_diff_d: null choi array (Rule 4)");
    assert(kraus_re != NULL && kraus_im != NULL &&
           "aic_ucp_choi_diff_d: null kraus array (Rule 4)");

    int N = n * n;

    aic_ucp_kraus phi, phi2;
    kraus_from_flat(&phi, kraus_re, kraus_im, n, r);
    aic_ucp_compose(&phi2, &phi, &phi, (slong) prec); /* Phi^2, n x n, r^2 ops */

    acb_mat_t J;
    acb_mat_init(J, N, N);
    aic_ucp_choi_diff(J, &phi2, &phi, (slong) prec);  /* Choi(Phi^2 - Phi)     */

    /* Extract J to the flat double arrays. The conjugation lives INSIDE the
     * Convention-A Choi (aic_ucp_kraus_to_choi); here we only copy real/imag
     * midpoints, so the round-trip preserves it iff we keep (p,q) and the imag
     * SIGN unchanged — the complex_qubit fixture asserts both. */
    for (int p = 0; p < N; p++)
        for (int q = 0; q < N; q++) {
            acb_srcptr e = acb_mat_entry(J, p, q);
            choi_re[p * N + q] = arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR);
            choi_im[p * N + q] = arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR);
        }

    acb_mat_clear(J);
    aic_ucp_kraus_clear(&phi2);
    aic_ucp_kraus_clear(&phi);
    return 0;
}

/* aic_cbnorm.h — eig-free certified bracket [||J||_F/n, 2*||J||_F] on
 * eta = ||Phi^2-Phi||_cb (.tex:347-354). Flat-double ABI: builds Phi, runs
 * aic_cbnorm_eigfree_ball, then converts the arb balls to a RIGOROUS double
 * bracket: *lo = FLOOR(arb_lower(lo_ball)), *hi = CEIL(arb_upper(hi_ball)) —
 * the only rounding directions that keep *lo <= eta <= *hi after the double
 * conversion. */
int aic_cbnorm_eigfree_d(double *lo, double *hi,
                         const double *kraus_re, const double *kraus_im,
                         int n, int r, int prec)
{
    assert(n >= 1 && "aic_cbnorm_eigfree_d: n >= 1 required (Rule 4)");
    assert(r >= 1 && "aic_cbnorm_eigfree_d: r >= 1 required (Rule 4)");
    assert(prec >= 2 && "aic_cbnorm_eigfree_d: prec >= 2 required (Rule 4)");
    assert(lo != NULL && hi != NULL &&
           "aic_cbnorm_eigfree_d: null output scalar (Rule 4)");
    assert(kraus_re != NULL && kraus_im != NULL &&
           "aic_cbnorm_eigfree_d: null kraus array (Rule 4)");

    aic_ucp_kraus phi;
    kraus_from_flat(&phi, kraus_re, kraus_im, n, r);

    arb_t loB, hiB;
    arb_init(loB);
    arb_init(hiB);
    aic_cbnorm_eigfree_ball(loB, hiB, &phi, (slong) prec);

    /* RIGOROUS double bracket: round the lower endpoint DOWN, the upper endpoint
     * UP, so the double bracket still encloses [arb_lower(loB), arb_upper(hiB)]
     * (and hence eta). arb_get_lbound_arf gives the ball's lower endpoint; FLOOR
     * keeps *lo <= it. arb_get_ubound_arf gives the upper; CEIL keeps *hi >= it. */
    arf_t t;
    arf_init(t);
    arb_get_lbound_arf(t, loB, (slong) prec);
    *lo = arf_get_d(t, ARF_RND_FLOOR);
    arb_get_ubound_arf(t, hiB, (slong) prec);
    *hi = arf_get_d(t, ARF_RND_CEIL);
    arf_clear(t);

    arb_clear(hiB);
    arb_clear(loB);
    aic_ucp_kraus_clear(&phi);
    return 0;
}
