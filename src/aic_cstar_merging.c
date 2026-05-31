/* aic_cstar_merging.c — lem_merging, the GENERAL 2x2 block-assembly lemma of the
 * th_main master loop (bead aic-097, module "cstar_build", Increment 3). See
 * include/aic_cstar.h (Increment 3 banner) and docs/research/cstar_build_design.md
 * §4.3 (lem_merging) and §4.4 Step 5-6 (how lem_extension's four gamma_{jk} feed
 * lem_merging). lem_extension (I4) is the only consumer.
 *
 * approximate_algebras.tex:1325-1346 — Lemma lem_merging:
 *   "Let Pi_1,Pi_2 be projections in a C* algebra B such that Pi_1+Pi_2=I, and
 *    similarly, let P_1,P_2 be delta-projections in an eps-C* algebra A such that
 *    ||P_1+P_2-I|| <= delta. Consider some linear maps
 *      gamma_{jk} : S_{Pi_j,Pi_k} -> S_{P_j,P_k}  (for j,k in {1,2})
 *    satisfying the following conditions:
 *      (merging0) gamma_{kj}(X^dag) = gamma_{jk}(X)^dag,
 *      (merging1) ||gamma_{jl}(XY) - gamma_{jk}(X) . gamma_{kl}(Y)|| <= delta||X||||Y||,
 *      (merging2) ||gamma_{jj}(Pi_j) - P_j|| <= delta,
 *      (merging3) (1-delta)||X|| <= ||gamma_{jk}(X)|| <= (1+delta)||X||.
 *    (The dot in (merging1) denotes the compressed product (compr_prod).) Then the
 *    combined map gamma : (X_{jk}) |-> sum_{jk} gamma_{jk}(X_{jk}) : B -> A is an
 *    O(delta+eps)-inclusion. If all maps gamma_{jk} are bijective, then gamma is
 *    also bijective."
 *
 * PAPER TECHNIQUE vs CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The proof (.tex:1348)
 *   writes gamma = alpha (sum_{jk} gamma_{jk}) mu with mu the canonical bijection
 *   B -> (+)_{jk} S_{Pi_j,Pi_k} and alpha the lem_alpha block-merge map; merging0-2
 *   give an O(delta)-homomorphism, merging3 + lem_alpha give norm equivalence, and
 *   prop_delta_hominc upgrades to an O(delta+eps)-inclusion. In finite dim B = M_d
 *   (d = n1+n2) is concrete: Pi_1 = proj onto the first n1 coords, Pi_2 the last n2,
 *   the block S_{Pi_j,Pi_k} is the (j,k) rectangular sub-block of M_d, and mu is the
 *   identity bookkeeping that reads a matrix unit's (row,col) sector. So gamma is
 *   ROUTING: each matrix unit E_{lm} of M_d is sent to the gamma_{jk} image of its
 *   block-local unit. We assemble the resulting flat vE[] (the aic_dhom_v model) and
 *   CERTIFY the inclusion with the closed dhom routines:
 *     - merging1 (the compressed-product multiplicativity) <- aic_dhom_defect_sweep,
 *       which is exactly G_gamma(E_i,E_j) over all matrix-unit pairs, w.r.t. A's STAR
 *       (NOT plain product; FINDINGS §C2 — the eta=0 oracle is blind, the oblique
 *       fixture supplies the teeth, and per FINDINGS §C8 the magnitude c=defect/eta
 *       is the sharp discriminant on in-A fixtures);
 *     - the inclusion lower bound (a, b of prop_delta_hominc) <- aic_dhom_v_sigma_min,
 *       the SOUND combination-level Frobenius lower bound (FINDINGS §C6);
 *     - merging0/2/3 are the CHEAP per-block input guards (involution-compat, unit,
 *       norm-equivalence), computed directly into merge_cond_max (reportable, not
 *       aborting; the I5 loop wants a fail-loud version).
 *
 * THE TWO B-SHAPES (FINDINGS §C9, the central subtlety this increment surfaced).
 *   .tex:1325 only says B is "a C* algebra" with Pi_1 + Pi_2 = I. TWO shapes arise:
 *     two_block=0  B = M_{n1+n2} (single block): off-diagonal sectors S_{Pi_1,Pi_2},
 *       S_{Pi_2,Pi_1} are LIVE (E_{lm} across the n1 boundary exist), so gamma_12,
 *       gamma_21 must be NONZERO inclusions — the shape lem_extension (I4) needs.
 *     two_block=1  B = M_{n1} (+) M_{n2} (two blocks): off-diagonal sectors are
 *       EMPTY (block-diagonal B has no cross-block units), so gamma_12 = gamma_21 = 0
 *       (zero-dim domain) — the shape cor_merge_sum (.tex:1352) uses; merged v then
 *       equals aic_cstar_merge_sum's concat exactly.
 *   A single M_{n1+n2} block with zeroed off-diagonal gamma is NOT a valid input:
 *   the pair (E_{0,n1}, E_{n1,0}) has E*E = E_00 but gamma(E)=0, so merging1 is
 *   violated and mult_def is O(1) (= ||Ptilde_1||), NOT the cor_merge_sum cross-
 *   defect ||Ptilde_1 * Ptilde_2|| ~ O(eps). Conflating the two is a silent bug.
 *
 * THE ROUTING (the module's highest convention-risk; an off-by-one mis-assembles
 *   gamma silently and is exactly what a hostile review hunts). For B = M_d (single
 *   block, two_block=0) the matrix unit E_{lm} has linear index mu(0,l,m) = l*d + m
 *   (aic_dhom_B_index). l,m are GLOBAL row/col; the block (j,k) and local (a,b) are
 *   read off the sector boundaries n1, n2; the gamma_{jk} input array is indexed
 *   row-major within its block (a*cols + b, cols = n1 for k=1 else n2). See the
 *   route_unit helper — the single point where the convention lives. For two_block=1
 *   the routing is the trivial block-diagonal one (assemble_two_block).
 */
#include <assert.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner_internal.h" /* aic_corner_gamma_opnorm_ub (§C5 UB workaround) */
#include "aic/aic_cstar.h"
#include "aic/aic_dhom.h"
#include "aic/aic_ecstar.h"

/* Route a global matrix unit E_{lm} of B = M_{n1+n2} to the gamma_{jk} image of
 * its block-local unit. Returns a BORROWED pointer to the input image operator
 * (an acb_mat_struct *, the element of the relevant gammaJK array). This is the
 * single point where the (j,k)/(a,b) bookkeeping lives (.tex:1339-1344). */
static const acb_mat_struct *route_unit(slong l, slong m, slong n1, slong n2,
                                        const acb_mat_t *g11, const acb_mat_t *g12,
                                        const acb_mat_t *g21, const acb_mat_t *g22)
{
    /* Array column counts: g11=n1, g12=n2, g21=n1, g22=n2 (n2 used in the g12/g22
     * index expressions below). The local (a,b) is (l,m) shifted by the sector
     * offsets n1 on whichever of row/col crosses the boundary. */
    if (l < n1 && m < n1)  return g11[l * n1 + m];               /* block (1,1) */
    if (l < n1 && m >= n1) return g12[l * n2 + (m - n1)];        /* block (1,2) */
    if (l >= n1 && m < n1) return g21[(l - n1) * n1 + m];        /* block (2,1) */
    return g22[(l - n1) * n2 + (m - n1)];                        /* block (2,2) */
}

/* ||A - B||_op via the certified mid+radius upper bound (FINDINGS §C5: the bare
 * aic_mat_opnorm Gram path false-fails on near-zero-off-diagonal differences). */
static void opnorm_diff(arb_t out, const acb_mat_t Aop, const acb_mat_t Bop,
                        slong prec)
{
    acb_mat_t D;
    acb_mat_init(D, acb_mat_nrows(Aop), acb_mat_ncols(Aop));
    acb_mat_sub(D, Aop, Bop, prec);
    aic_corner_gamma_opnorm_ub(out, NULL, D, prec);
    acb_mat_clear(D);
}

/* The cheap per-block INPUT guards merging0 + merging2 + merging3 (.tex:1329/1334/
 * 1336), max-reduced into merge_cond_max. Op-norms via the certified UB (§C5). When
 * two_block=1 the off-diagonal blocks are EMPTY, so only the diagonal contributes. */
static void merge_cond_sweep(arb_t out, int two_block, slong n1, slong n2,
                             const acb_mat_t *g11, const acb_mat_t *g12,
                             const acb_mat_t *g21, const acb_mat_t *g22,
                             const aic_ecstar *A, const acb_mat_t P1,
                             const acb_mat_t P2, slong prec)
{
    slong n = A->n;
    arb_t e;
    arb_init(e);
    arb_zero(out);

    /* merging3 (.tex:1336): | ||gamma_{jk}(E_{ab})||_op - 1 | (||E_{ab}||_op = 1).
     * Sweep the diagonal blocks always; the off-diagonal only when they are LIVE. */
    acb_mat_t dag;
    acb_mat_init(dag, n, n);
    int nblk = two_block ? 2 : 4;
    struct { const acb_mat_t *g; slong rows, cols; } blk[4] = {
        {g11, n1, n1}, {g22, n2, n2}, {g12, n1, n2}, {g21, n2, n1}
    };
    for (int B_id = 0; B_id < nblk; B_id++) {
        slong rows = blk[B_id].rows, cols = blk[B_id].cols;
        const acb_mat_t *g = blk[B_id].g;
        for (slong a = 0; a < rows; a++)
            for (slong b = 0; b < cols; b++) {
                /* merging3: ||gamma(E_ab)||_op should be 1. */
                aic_corner_gamma_opnorm_ub(e, NULL, g[a * cols + b], prec);
                arb_sub_si(e, e, 1, prec);
                arb_abs(e, e);
                arb_max(out, out, e, prec);
            }
    }

    /* merging0 (.tex:1329): ||gamma_{kj}(E_{ba}) - gamma_{jk}(E_{ab})^dag||. The
     * diagonal blocks internally (gamma_jj(E_ba) vs gamma_jj(E_ab)^dag); the
     * off-diagonal couples (1,2)<->(2,1) — only when LIVE (two_block=0). */
    const acb_mat_t *gdiag[2] = {g11, g22};
    slong ndiag[2] = {n1, n2};
    for (int d = 0; d < 2; d++)
        for (slong a = 0; a < ndiag[d]; a++)
            for (slong b = 0; b < ndiag[d]; b++) {
                acb_mat_conjugate_transpose(dag, gdiag[d][a * ndiag[d] + b]);
                opnorm_diff(e, gdiag[d][b * ndiag[d] + a], dag, prec);
                arb_max(out, out, e, prec);
            }
    if (!two_block)
        for (slong a = 0; a < n1; a++)
            for (slong b = 0; b < n2; b++) {
                acb_mat_conjugate_transpose(dag, g12[a * n2 + b]); /* g12(E_ab)^dag */
                opnorm_diff(e, g21[b * n1 + a], dag, prec);        /* g21(E_ba)     */
                arb_max(out, out, e, prec);
            }

    /* merging2 (.tex:1334): ||gamma_{jj}(Pi_j) - P_j||, Pi_j = sum_a E^{(jj)}_{aa}.
     * gamma is linear, so gamma_{jj}(Pi_j) = sum_a gamma_{jj}(E_aa). */
    acb_mat_t accP;
    acb_mat_init(accP, n, n);
    acb_mat_zero(accP);
    for (slong a = 0; a < n1; a++)
        acb_mat_add(accP, accP, g11[a * n1 + a], prec);
    opnorm_diff(e, accP, P1, prec);
    arb_max(out, out, e, prec);
    acb_mat_zero(accP);
    for (slong a = 0; a < n2; a++)
        acb_mat_add(accP, accP, g22[a * n2 + a], prec);
    opnorm_diff(e, accP, P2, prec);
    arb_max(out, out, e, prec);

    acb_mat_clear(accP);
    acb_mat_clear(dag);
    arb_clear(e);
}

void aic_cstar_lem_merging(aic_dhom_B *B_out, aic_dhom_v *v_out,
                           arb_t mult_def, arb_t sigma_min, arb_t merge_cond_max,
                           int two_block, slong n1, slong n2,
                           const acb_mat_t *gamma11, const acb_mat_t *gamma12,
                           const acb_mat_t *gamma21, const acb_mat_t *gamma22,
                           const aic_ecstar *A,
                           const acb_mat_t P1, const acb_mat_t P2, slong prec)
{
    assert(B_out != NULL && v_out != NULL && A != NULL);
    assert(gamma11 != NULL && gamma22 != NULL &&
           "aic_cstar_lem_merging: diagonal gamma arrays required");
    assert((two_block || (gamma12 != NULL && gamma21 != NULL)) &&
           "aic_cstar_lem_merging: off-diagonal gamma required when two_block=0");
    assert(n1 >= 1 && n2 >= 1 && "aic_cstar_lem_merging: need n1,n2 >= 1");
    assert(acb_mat_nrows(P1) == A->n && acb_mat_nrows(P2) == A->n &&
           "aic_cstar_lem_merging: P_j size != A->n");

    if (two_block) {
        /* B = M_{n1} (+) M_{n2}: block-diagonal, no cross-block units. Block 0's
         * E_{ab} (index a*n1+b) -> gamma11; block 1's E_{ab} (index n1*n1 + a*n2+b)
         * -> gamma22. (FINDINGS §C9; the genuine cor_merge_sum shape, .tex:1352.) */
        slong dims[2] = {n1, n2};
        aic_dhom_B_init(B_out, dims, 2);
        aic_dhom_v_init(v_out, B_out, A);
        for (slong a = 0; a < n1; a++)
            for (slong b = 0; b < n1; b++)
                acb_mat_set(v_out->vE[a * n1 + b], gamma11[a * n1 + b]);
        slong off1 = n1 * n1;
        for (slong a = 0; a < n2; a++)
            for (slong b = 0; b < n2; b++)
                acb_mat_set(v_out->vE[off1 + a * n2 + b], gamma22[a * n2 + b]);
    } else {
        /* B = M_{n1+n2}: single block, live off-diagonal. Route each matrix unit
         * E_{lm} (index l*d+m) to its gamma_{jk} image (.tex:1339-1344). DEEP COPY
         * (acb_mat_set), never alias the input. */
        slong d = n1 + n2;
        slong dims[1] = {d};
        aic_dhom_B_init(B_out, dims, 1);
        aic_dhom_v_init(v_out, B_out, A);
        for (slong l = 0; l < d; l++)
            for (slong m = 0; m < d; m++) {
                const acb_mat_struct *img =
                    route_unit(l, m, n1, n2, gamma11, gamma12, gamma21, gamma22);
                acb_mat_set(v_out->vE[l * d + m], img);
            }
    }

    /* merging1 (.tex:1331): the multiplicativity defect w.r.t. A's STAR (the
     * blockwise compressed-product condition; FINDINGS §C2). */
    if (mult_def != NULL) aic_dhom_defect_sweep(mult_def, v_out, prec);
    /* the inclusion lower bound a (prop_delta_hominc): SOUND combination-level
     * Frobenius lower bound (FINDINGS §C6). */
    if (sigma_min != NULL) aic_dhom_v_sigma_min(sigma_min, v_out, prec);
    /* merging0 + merging2 + merging3: the cheap per-block input guards (reportable). */
    if (merge_cond_max != NULL)
        merge_cond_sweep(merge_cond_max, two_block, n1, n2, gamma11, gamma12,
                         gamma21, gamma22, A, P1, P2, prec);
}
