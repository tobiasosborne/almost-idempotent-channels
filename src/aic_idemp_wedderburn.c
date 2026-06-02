/* aic_idemp_wedderburn.c — the Artin-Wedderburn / multiplicity decomposition of
 * the *-algebra A = Img Phi for an EXACTLY idempotent UCP map (bead aic-ynu, I1).
 * Realizes prop_hom_structure (.tex:259-272) + prop_Gamma (.tex:2106) for eta=0.
 *
 * approximate_algebras.tex:259-264 — prop_hom_structure (verbatim core):
 *   "M = (+)_j L_j (x) E_j, and ... w(B) = sum_j W_j^dag (B_j (x) 1_{E_j}) W_j,
 *    where W_j : M -> L_j (x) E_j and sum_j W_j^dag W_j = 1_M."
 *
 * THE OBJECT. d->w stores the *-homomorphism w : A -> B(M): column k is
 * vec(w(B_k)) (row-major, vec(X)[i*m+j]=X[i,j]), with m = dim_M and w(B_k) =
 * J_M^dag B_k J_M (.tex:2084). The {w(B_k)} span a *-subalgebra W <= B(M); by
 * Artin-Wedderburn W ~ (+)_j B(L_j) acting on M = (+)_j L_j (x) E_j. We recover
 * num_blocks, dim L_j, dim E_j and the block-diagonalizing isometries W_j.
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (Law 3). The paper ASSERTS the Wedderburn
 * structure exists (the structure theorem for finite-dim *-algebras). In finite
 * dim it is a matrix computation: eigendecompose a generic Hermitian H_W in W;
 * its spectral projectors lie in W and equal (v v^dag) (x) 1_{E_j} for an L_j
 * eigenvector v (multiplicity dim E_j). Two such cluster projectors are in the
 * SAME block iff the algebra connects them (||Q_a w(B_k) Q_a'||_op > tol), since
 * every w(B_k) is block-diagonal across j. The blocks are the connected
 * components; dim L_j = #clusters in the block, dim E_j = the common cluster
 * multiplicity. W_j is assembled by aligning a single E_j onb across the dim L_j
 * reference vectors via the connecting partial isometries in W.
 *
 * WHY H_W MUST BE RANDOM (FINDINGS, the aic-ynu BLOCKER). A *deterministic* H_W =
 * sum_k r_k(w(B_k)+w(B_k)^dag) with fixed r_k is NOT generic: on an equal-size
 * block algebra (e.g. M_6 (+) M_6, dim_M=12) the same fixed combination lands the
 * two blocks' spectra in lock-step, manufacturing a >=2-dim kernel that SPANS the
 * two blocks. Gap-clustering then merges two distinct Wedderburn blocks into one
 * cluster and the per-block multiplicity check spuriously aborts a VALID input.
 * The cure is genuine genericity: draw a_k, b_k uniformly in [-1,1] and form the
 * generic Hermitian element  H_W = sum_k a_k(w(B_k)+w(B_k)^dag)
 *                                   + i b_k(w(B_k)-w(B_k)^dag)  of the *-algebra
 * (the i-part is the anti-Hermitian generator i(w-w^dag), so {a_k,b_k} sweep a
 * generic real Hermitian element of W, not just its symmetric real span). Almost
 * every draw gives a simple spectrum within each block AND distinct spectra ACROSS
 * equal-size blocks, so the kernel does not span blocks. We DRAW REPRODUCIBLY (a
 * fixed-seed LCG, AIC_ADV_SEED convention — no wall-clock; same input -> same
 * decomposition every run) and RE-DRAW (up to WEDD_MAX_DRAWS) if a draw yields a
 * cluster that still spans two would-be blocks (the multiplicity-disagreement
 * symptom). Only the genuinely-unresolvable input fails loud after the budget.
 *
 * NUMBER PATH. Like aic_idemp (which extracts in the LAPACK double path and
 * CERTIFIES the defects in arb), the geometric extraction here runs in the double
 * path (aic_latd Hermitian eig + opnorm); W_j is written out as acb_mat at `prec`.
 * The w-reconstruction invariant ||sum_j W_j^dag (A_j (x) 1) W_j - w(B_k)|| is the
 * gauge-invariant correctness certificate; aic_wedd_certify_W certifies it (and
 * the dim identities + W-unitarity) in arb in PRODUCTION (not only in the tests),
 * so a tensor-MISALIGNED W (an L<->E swap that is still unitary) is rejected.
 *
 * SCOPE (I1). At eta=0 the random-Hermitian simple-spectrum extraction handles
 * BOTH distinct AND equal block sizes (block_cond_exp, including the equal-size
 * M_d (+) M_d), single blocks (noiseless_subsystem, identity), and the diagonal
 * M_1^n (dephasing). The ONLY fail-loud (Rule 4) is the genuinely-unresolvable
 * spectrum: no draw in the WEDD_MAX_DRAWS budget yields a clean block split, or
 * the assembled W fails the arb certification.
 */
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <flint/flint.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic_idemp_wedderburn_internal.h"

/* Re-draw budget for the generic Hermitian H_W (see file banner). One generic
 * draw separates equal-size blocks almost surely; a handful covers the measure-
 * zero collision set. Exhausting it is the genuinely-unresolvable fail-loud. */
#define WEDD_MAX_DRAWS 8

/* Fixed PRNG seed (the corpus AIC_ADV_SEED convention, golden-ratio splitmix
 * constant; tests/aic_adversarial.h:75). Reproducible: same input -> same H_W
 * draw stream -> same decomposition every run. NO wall-clock. */
#define WEDD_SEED 0x9e3779b97f4a7c15UL

/* Fixed-seed LCG (the AIC_ADV_SEED corpus convention, no wall-clock) returning a
 * reproducible rational in [-1,1) with denominator 2^16. Same input -> same
 * coefficient stream -> same H_W -> same decomposition on every run. */
static double wedd_lcg_unit(unsigned long *state)
{
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    long q = (long) ((*state >> 32) & 0xFFFFUL);   /* 0 .. 65535 */
    return (double) (q - 32768) / 32768.0;         /* [-1, 1)    */
}

/* reshape column k of d->w (length m^2, row-major) into the m x m double array. */
static void w_col_to_arr(double _Complex *Wk, const acb_mat_t w, slong k, slong m)
{
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++) {
            acb_srcptr e = acb_mat_entry(w, i * m + j, k);
            double re = arf_get_d(arb_midref(acb_realref(e)), ARF_RND_NEAR);
            double im = arf_get_d(arb_midref(acb_imagref(e)), ARF_RND_NEAR);
            Wk[i * m + j] = re + im * I;
        }
}

void aic_idemp_wedderburn_decompose(aic_idemp_wedderburn *out,
                                    const aic_idemp_decomp *d, slong prec)
{
    slong m = d->dim_M, dA = d->dim_A;
    assert(m >= 1 && dA >= 1 && "wedderburn: empty decomposition");

    out->dim_M = m;
    out->dim_A = dA;
    out->W_j = NULL;        /* so aic_..._clear is safe even if we never assemble */
    out->dim_L = NULL;
    out->dim_E = NULL;
    out->num_blocks = 0;

    /* {w(B_k)} as double arrays (the *-subalgebra W <= B(C^m)). */
    double _Complex **Wk = flint_malloc((size_t) dA * sizeof(*Wk));
    for (slong k = 0; k < dA; k++) {
        Wk[k] = flint_malloc((size_t) (m * m) * sizeof(double _Complex));
        w_col_to_arr(Wk[k], d->w, k, m);
    }

    double _Complex *HW = flint_malloc((size_t) (m * m) * sizeof(double _Complex));
    double *evals = flint_malloc((size_t) m * sizeof(double));
    double _Complex *evecs = flint_malloc((size_t) (m * m) * sizeof(double _Complex));
    slong *cl_start = flint_malloc((size_t) m * sizeof(slong));
    slong *cl_size = flint_malloc((size_t) m * sizeof(slong));

    /* Step 1: a RANDOM generic Hermitian H_W in the *-algebra W, with re-draw.
     * H_W = sum_k a_k (w(B_k)+w(B_k)^dag) + i b_k (w(B_k)-w(B_k)^dag); a_k,b_k in
     * [-1,1) from the fixed-seed LCG. The i(w-w^dag) term makes H_W sweep a generic
     * Hermitian element of W (not just the real-symmetric span), so its spectrum is
     * simple within each block and distinct across equal-size blocks a.s. The seed
     * advances each draw, so retries are fresh yet reproducible. */
    unsigned long seed = WEDD_SEED;
    int ok = 0;
    for (int draw = 0; draw < WEDD_MAX_DRAWS; draw++) {
        for (slong i = 0; i < m * m; i++) HW[i] = 0;
        for (slong k = 0; k < dA; k++) {
            double ak = wedd_lcg_unit(&seed);
            double bk = wedd_lcg_unit(&seed);
            for (slong i = 0; i < m; i++)
                for (slong j = 0; j < m; j++) {
                    double _Complex wij = Wk[k][i * m + j];
                    double _Complex wji_d = conj(Wk[k][j * m + i]); /* (w^dag)[i,j] */
                    HW[i * m + j] += ak * (wij + wji_d)
                                   + I * bk * (wij - wji_d);
                }
        }

        /* Step 1b: Hermitian eig (double path): ascending evals + onb evecs. */
        aic_latd_eig_hermitian(evals, evecs, HW, m);

        /* Step 2: gap-cluster -> cluster column ranges of evecs. Each cluster =
         * an L_j-eigenvector tensor E_j subspace (dim = dim E_j). */
        slong n_cl = aic_wedd_cluster_evals(cl_start, cl_size, evals, m);

        /* Step 3: group clusters into blocks + build W_j. On a degenerate draw
         * (a cluster spanning two would-be blocks -> per-block multiplicity
         * disagreement) this returns nonzero; we re-draw. On the LAST draw it
         * is asked to fail loud instead. The arb w-reconstruction + unitarity
         * certificate runs inside aic_wedd_assemble (production self-cert). */
        int last = (draw + 1 == WEDD_MAX_DRAWS);
        int rc = aic_wedd_assemble(out, Wk, dA, evecs, cl_start, cl_size, n_cl,
                                   m, prec, /*fail_loud=*/last);
        if (rc == 0) { ok = 1; break; }
        /* rc != 0: a degenerate draw on a non-final attempt; out left empty. */
    }
    if (!ok) {
        flint_printf("aic_idemp_wedderburn: no clean block split in %d random "
                     "H_W draws (genuinely-unresolvable spectrum); fail loud "
                     "(Rule 4)\n", WEDD_MAX_DRAWS);
        flint_abort();
    }

    /* --- cleanup --- */
    flint_free(cl_size);
    flint_free(cl_start);
    flint_free(evecs);
    flint_free(evals);
    flint_free(HW);
    for (slong k = 0; k < dA; k++) flint_free(Wk[k]);
    flint_free(Wk);
}

void aic_idemp_wedderburn_clear(aic_idemp_wedderburn *out)
{
    if (out == NULL || out->W_j == NULL) return;
    for (slong j = 0; j < out->num_blocks; j++) acb_mat_clear(out->W_j[j]);
    flint_free(out->W_j);
    flint_free(out->dim_L);
    flint_free(out->dim_E);
    out->W_j = NULL;
    out->dim_L = NULL;
    out->dim_E = NULL;
    out->num_blocks = 0;
}
