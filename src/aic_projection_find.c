/* aic_projection_find.c — the Hermitian-part + single-gap primitives of the
 * Route-A-ambient nontrivial-projection finder (bead aic-mqf). See
 * include/aic_projection.h for the full route and lem_nontriv_projection
 * (approximate_algebras.tex:931).
 *
 * approximate_algebras.tex:935 (the reduction this whole module serves):
 *   "An element P in A is a projection if and only if X = 2P - I is Hermitian and
 *    satisfies X^2 = I. Equivalently, X is both Hermitian and unitary."
 * We realize X as the AMBIENT sign of a gap-shifted Hermitian element of A. The
 * eigenVALUES are the DOUBLE path (LAPACK zheev): H is Hermitian, whose spectrum
 * is STABLE (the .tex:540-544 eigenvalue fragility is a NON-normal phenomenon), so
 * the double-path values are a sound basis for choosing a threshold. The certified
 * DEFECT is the arb path (aic_projection.c); the certified gap ENCLOSURE defers to
 * bead aic-w4o.1.
 *
 * POST-aic-66n. The old single-H pre-selection aic_projection_pick_H was REMOVED:
 * it ranked H_k by the largest AMBIENT interior gap, which on an S_P wrapper is the
 * support-vs-complement gap, collapsing the split to the wrapper unit Ptilde_m
 * (FINDINGS §C11/§C16). The UNIT-AWARE audition (src/aic_projection_audit.c)
 * supersedes it. This file now provides only aic_projection_herm_part (H_k) and
 * aic_projection_gap (the single largest-gap finder, directly tested in
 * tests/test_projection.c: the aic-3qv no-gap abort + the gap-pick teeth).
 *
 * FAIL-LOUD (Rule 4): aic_projection_gap aborts when there is no positive interior
 * gap (a degenerate spectrum, all eigenvalues coincident): the aic-3qv stop
 * condition, certified per-instance.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_ecstar.h"
#include "aic/aic_latd.h"
#include "aic_projection_internal.h"

/* H_k = (B_k + B_k^dag)/2 (Hermitian, in A since A is dag-closed). `out` is
 * caller-init'd n x n. */
void aic_projection_herm_part(acb_mat_t out, const acb_mat_t Bk, slong prec)
{
    slong n = acb_mat_nrows(Bk);
    acb_mat_t Bd;
    acb_mat_init(Bd, n, n);
    acb_mat_conjugate_transpose(Bd, Bk);
    acb_mat_add(out, Bk, Bd, prec);
    acb_mat_scalar_mul_2exp_si(out, out, -1);   /* * 1/2 */
    acb_mat_clear(Bd);
}

/* qsort comparator: ascending doubles. */
static int dbl_cmp(const void *p, const void *q)
{
    double a = *(const double *) p, b = *(const double *) q;
    return (a > b) - (a < b);
}

/* Largest INTERIOR gap of Hermitian H (double-path zheev eigenvalues). Returns
 * the gap; writes the spread (lam_max-lam_min) to *spread_out and, if non-NULL,
 * the sorted eigenvalues into ev_out[0..n) and the gap index to *best_i_out
 * (lam_{best_i}..lam_{best_i+1} straddle the threshold). Every i in 0..n-2 is
 * interior (both clusters non-empty), so the split is ALWAYS nontrivial; we rank
 * candidates by this gap (not spread) because argmax(spread) != argmax(gap).
 * Pure double path (the .tex:540-544 fragility is a non-normal phenomenon; H is
 * Hermitian, so zheev's values are stable). */
static double largest_interior_gap(const acb_mat_t H, double *spread_out,
                                   double *ev_out, slong *best_i_out)
{
    slong n = acb_mat_nrows(H);
    assert(n >= 2 && "projection gap: need n >= 2");
    double _Complex *Harr = (double _Complex *) malloc((size_t) n * n *
                                                       sizeof(double _Complex));
    if (!Harr) { fprintf(stderr, "largest_interior_gap: OOM\n"); abort(); }
    aic_latd_from_acb_mat(Harr, H);
    double *ev = (double *) malloc((size_t) n * sizeof(double));
    if (!ev) { fprintf(stderr, "largest_interior_gap: OOM\n"); abort(); }
    aic_latd_eig_hermitian(ev, NULL, Harr, n);   /* ascending */
    qsort(ev, (size_t) n, sizeof(double), dbl_cmp);   /* defensive (zheev is asc) */

    double best_gap = -1.0;
    slong best_i = -1;
    for (slong i = 0; i < n - 1; i++) {
        double g = ev[i + 1] - ev[i];
        if (g > best_gap) { best_gap = g; best_i = i; }
    }
    if (spread_out) *spread_out = ev[n - 1] - ev[0];
    if (best_i_out) *best_i_out = best_i;
    if (ev_out) for (slong i = 0; i < n; i++) ev_out[i] = ev[i];
    free(ev);
    free(Harr);
    return best_gap;
}

/* NOTE (bead aic-66n): the old single-H pre-selection aic_projection_pick_H was
 * REMOVED. It ranked H_k by the largest AMBIENT interior gap, which on an S_P
 * wrapper is the support-vs-complement gap, so the projection collapsed to the
 * wrapper unit Ptilde_m (FINDINGS §C11). The UNIT-AWARE audition
 * (src/aic_projection_audit.c, aic_projection_audit) supersedes it: it enumerates
 * the (H_k, gap) pairs over ALL H_k and selects by the ||U_A - P|| > c criterion,
 * not a pre-chosen H. aic_projection_gap below remains as the directly-tested
 * single-gap finder (tests/test_projection.c aic-3qv abort / gap-pick teeth). */

/* Largest INTERIOR spectral gap of a Hermitian H (eigenVALUES via LAPACK zheev,
 * degenerate-OK double path). Among the n-1 consecutive gaps lam_{i+1}-lam_i
 * (i = 0..n-2) pick the largest with BOTH sides non-empty (any i is interior:
 * lam_0..lam_i below, lam_{i+1}..lam_{n-1} above), guaranteeing nontriviality.
 * Writes threshold t = (lam_i + lam_{i+1})/2 to t_out, gap to g_out, lam_min and
 * lam_max to lmin, lmax, and m = n-1-i (the count of eigenvalues >= t, the
 * +1-cluster size) to m_out.
 *
 * FAILS LOUD only if there is NO positive interior gap (the spectrum is
 * degenerate, all eigenvalues coincide to round-off): the aic-3qv stop
 * condition. The gap SIZE is NOT gated here — it enters the a-posteriori star
 * defect delta=O(eps/g), certified downstream. spread = lam_max - lam_min. */
void aic_projection_gap(double *t_out, double *g_out, double *lmin, double *lmax,
                        slong *m_out, const acb_mat_t H, slong prec)
{
    slong n = acb_mat_nrows(H);
    assert(n >= 2 && "projection gap: need n >= 2");
    (void) prec;   /* the gap-finding is the double path (zheev, midpoints) */
    double *ev = (double *) malloc((size_t) n * sizeof(double));
    if (!ev) { fprintf(stderr, "aic_projection_gap: OOM\n"); abort(); }
    double spread = 0.0;
    slong best_i = -1;
    double best_gap = largest_interior_gap(H, &spread, ev, &best_i);
    /* g_min: the ONLY requirement is a genuine POSITIVE interior gap. Math: the
     * sgn basin needs ANY g>0 (||Y^2-I||=1-(s*g/2)^2<1, s=1/max(t-lmin,lmax-t));
     * nontriviality is guaranteed by the INTERIOR split (both clusters non-empty),
     * independent of gap SIZE. The gap SIZE only enters the a-posteriori star
     * defect delta=O(eps/g) at eta>0 (research finding) — which is CERTIFIED
     * a-posteriori (star_defect in aic_projection.c) and guarded by the T3
     * universality canary, NOT by this abort. So the abort is reserved for a
     * GENUINELY DEGENERATE spectrum (no positive interior gap), the aic-3qv stop
     * condition. The floor 1e-9*max(1,spread) only ensures the split is REAL
     * (above precision/round-off noise), not that the relative gap is modest;
     * the old 0.05*spread floor spuriously aborted well-conditioned splits at
     * n>=22 (uniform spectrum {0..n-1}: largest gap 1.0 < 0.05*21=1.05). */
    double g_min = 1e-9 * fmax(1.0, spread);
    if (best_gap < g_min) {
        fprintf(stderr, "aic_projection: NO positive interior spectral gap "
                "(largest gap %.3e < floor %.3e, spread %.3e) in the chosen "
                "Hermitian H — the spectrum is DEGENERATE (all eigenvalues "
                "coincide to round-off). lem_nontriv_projection's spectral split "
                "cannot proceed; escalate (bead aic-3qv).\n",
                best_gap, g_min, spread);
        free(ev);
        abort();
    }
    *t_out = 0.5 * (ev[best_i] + ev[best_i + 1]);
    *g_out = best_gap;
    if (lmin) *lmin = ev[0];
    if (lmax) *lmax = ev[n - 1];
    if (m_out) *m_out = n - 1 - best_i;   /* # eigenvalues >= t */
    free(ev);
}
