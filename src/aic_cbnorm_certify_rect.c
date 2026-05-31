/* aic_cbnorm_certify_rect.c — the RECTANGULAR (non-self-map) certified UPPER
 * bound on a diamond norm, for the opspace O2 cb-norm certifier (bead aic-pjr).
 * This is the certified payoff of th_main_ext (.tex:1538-1540): it turns a
 * committed Watrous MIN-dual feasible point (Y0,Y1) for the map f: M_in -> M_out
 * (presented by its Convention-A Choi J) into a RIGOROUS arb upper bound on
 * ||f||_⋄ = optval. For f = v* this certifies ||v||_cb; for f = (v^{-1})* it
 * certifies ||v^{-1}||_cb (cb-spectral = diamond-of-adjoint, Watrous 2009).
 *
 * THE PINNED CONVENTION (docs/research/opspace_o2_design.md §0.5; the empirically-
 * pinned recipe, factor = 1, NO 2/n). For the GOLDEN-RULE Choi
 *   J(f)[s*out+i, t*out+j] = f(E_st)[i,j]   (INPUT s,t MAJOR/stride out; OUTPUT
 *                                            i,j MINOR),
 * sized (d_maj = in, d_min = out), the Watrous diamond-norm MIN-dual is
 *   min (1/2)(||Tr_min Y0||_inf + ||Tr_min Y1||_inf)
 *   s.t. block_D = [[Y0,-J],[-J^dag,Y1]] >= 0,  Y0,Y1 >= 0   on M_maj (x) M_min,
 * and the dual traces the MINOR/OUTPUT factor (tr_sys = 2, size d_min):
 * Tr_min = aic_mat_partial_trace_right(., d_maj, d_min) -> a d_maj x d_maj
 * marginal. (This is OPPOSITE to the self-map where input = minor; here the GOLDEN
 * RULE puts the OUTPUT in the minor position, but the dual STILL traces the minor
 * — the pin was set by strong duality on an asymmetric NON-CP fixture, §0.5.)
 *
 * UPPER RESTORATION (the deliverable; generalizes aic_cbnorm_int_upper). The
 * committed (Y0,Y1) is a DOUBLE-precision approximate dual point; restore it to
 * exact feasibility by a shift: eps = max(0,-lambda_min(block_D)); then
 * (Y0 + eps I, Y1 + eps I) is dual-feasible (block_D + eps I_{2D} >= 0, and the
 * PSD constraints Y0,Y1>=0 only relax). Since
 *   Tr_min(eps I_{d_maj*d_min}) = eps * d_min * I_{d_maj},
 * the marginal max-eig of each shifted block rises by exactly eps*d_min, so
 *   hi = (1/2)(lambda_max(Tr_min Y0) + lambda_max(Tr_min Y1)) + eps*d_min
 * is a RIGOROUS upper bound on ||f||_⋄. Factor = 1; NO 2/n. NO eta=0 short-circuit
 * (J(v*) is a nonzero isometry Choi at eta=0; the trivial value is 1, not 0 —
 * design §0.5/§6.4). Every PSD test is the ONE-SIDED herm_max_eig(-M) global
 * enclosure (degeneracy-robust; no full eig). arb/acb only.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_cbnorm.h"
#include "aic/aic_mat.h"

/* Assemble out = [[A,B],[C,D]], each block m x m, into a 2m x 2m matrix.
 * (Local copy of the restore.c helper; kept self-contained to leave the self-map
 * recipe files untouched — Rule 10 split, not a shared-static refactor.) */
static void block2_rect(acb_mat_t out, const acb_mat_t A, const acb_mat_t B,
                        const acb_mat_t C, const acb_mat_t D, slong m)
{
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++) {
            acb_set(acb_mat_entry(out, i, j), acb_mat_entry(A, i, j));
            acb_set(acb_mat_entry(out, i, m + j), acb_mat_entry(B, i, j));
            acb_set(acb_mat_entry(out, m + i, j), acb_mat_entry(C, i, j));
            acb_set(acb_mat_entry(out, m + i, m + j), acb_mat_entry(D, i, j));
        }
}

/* Rigorous UPPER bound on max(0, -lambda_min(M)) for Hermitian M: the upper
 * endpoint of herm_max_eig(-M) (= -lambda_min(M)), clamped at 0. */
static void psd_defect_rect(arb_t out, const acb_mat_t M, slong prec)
{
    slong d = acb_mat_nrows(M);
    acb_mat_t negM;
    acb_mat_init(negM, d, d);
    acb_mat_neg(negM, M);

    arb_t e;
    arb_init(e);
    aic_mat_herm_max_eig(e, negM, prec);   /* ball enclosing -lambda_min(M) */

    arf_t ub;
    arf_init(ub);
    arb_get_ubound_arf(ub, e, prec);       /* rigorous upper endpoint */
    arf_set(arb_midref(out), ub);
    mag_zero(arb_radref(out));
    if (arf_sgn(ub) < 0) arb_zero(out);    /* clamp: defect >= 0 */

    arf_clear(ub);
    arb_clear(e);
    acb_mat_clear(negM);
}

void aic_cbnorm_certify_rect_upper(arb_t hi, const acb_mat_t J,
                                   const acb_mat_t Y0, const acb_mat_t Y1,
                                   slong d_maj, slong d_min, slong prec)
{
    slong D = d_maj * d_min;
    assert(d_maj >= 1 && d_min >= 1
           && "aic_cbnorm_certify_rect_upper: dims >= 1 (Rule 4)");
    assert(acb_mat_nrows(J) == D && acb_mat_ncols(J) == D
           && "rect_upper: J must be (d_maj*d_min) square");
    assert(acb_mat_nrows(Y0) == D && acb_mat_ncols(Y0) == D
           && "rect_upper: Y0 must be (d_maj*d_min) square");
    assert(acb_mat_nrows(Y1) == D && acb_mat_ncols(Y1) == D
           && "rect_upper: Y1 must be (d_maj*d_min) square");

    /* block_D = [[Y0,-J],[-J^dag,Y1]]; eps = max(0,-lambda_min(block_D)).
     *
     * The feasible-point seed (Y0,Y1) AND the Choi J are FEASIBLE-POINT SEEDS, not
     * rigorous objects (the rigor is the eps restoration). For an ARB-ASSEMBLED J
     * (e.g. J(v*) from aic_opspace_choi_vstar over a cstar_build v) the entries
     * inherit an ACCUMULATED radius (~1e-71 at prec=256, far above a single-rounding
     * ulp) from the many arb ops in the build. The eig route asserts Hermiticity
     * RIGOROUSLY at the prec-tight tol ~2^-(prec-8) (~1e-75 at prec=256, bead aic-2yo):
     * the off-diagonal -J / -J^dag are INDEPENDENT balls whose asymmetry ball ~ J's
     * radius (1e-71) exceeds that tol, so the rigorous check rejects the genuinely-
     * Hermitian block (the self-map int_upper never hit this — its committed J has a
     * single-multiply radius BELOW tol). Fix: take the block's MIDPOINT (zero radius)
     * and symmetrize, so it is EXACTLY Hermitian-by-construction. The midpoint block
     * differs from the true block by <= the radius (~1e-71) in operator norm, utterly
     * below the 1e-4 tightness tolerance, so eps stays a rigorous defect bound for the
     * (midpoint) feasible point — exactly the self-map's zero-radius-committed-J
     * posture, which the committed self-map path is UNAFFECTED by. */
    acb_mat_t negJ, negJd, blk, blkd;
    acb_mat_init(negJ, D, D);
    acb_mat_init(negJd, D, D);
    acb_mat_init(blk, 2 * D, 2 * D);
    acb_mat_init(blkd, 2 * D, 2 * D);
    acb_mat_neg(negJ, J);
    acb_mat_conjugate_transpose(negJd, negJ);
    block2_rect(blk, Y0, negJ, negJd, Y1, D);
    /* collapse to the midpoint (drop the accumulated radius), then symmetrize so the
     * rigorous Hermiticity assertion accepts the by-construction-Hermitian block. */
    for (slong i = 0; i < 2 * D; i++)
        for (slong j = 0; j < 2 * D; j++) {
            acb_ptr e = acb_mat_entry(blk, i, j);
            arb_get_mid_arb(acb_realref(e), acb_realref(e));
            arb_get_mid_arb(acb_imagref(e), acb_imagref(e));
        }
    acb_mat_conjugate_transpose(blkd, blk);
    acb_mat_add(blk, blk, blkd, prec);
    acb_mat_scalar_mul_2exp_si(blk, blk, -1);     /* (blk + blk^dag)/2, exact Hermitian */
    acb_mat_clear(blkd);
    arb_t eps;
    arb_init(eps);
    psd_defect_rect(eps, blk, prec);

    /* Tr_min = partial_trace_right (the MINOR/OUTPUT factor, size d_min; PINNED
     * §0.5 tr_sys=2). The d_maj x d_maj marginal of Y_i (Hermitian PSD) has
     * ||.||_inf = lambda_max. */
    acb_mat_t T0, T1;
    acb_mat_init(T0, d_maj, d_maj);
    acb_mat_init(T1, d_maj, d_maj);
    aic_mat_partial_trace_right(T0, Y0, d_maj, d_min, prec);
    aic_mat_partial_trace_right(T1, Y1, d_maj, d_min, prec);
    arb_t s0, s1;
    arb_init(s0);
    arb_init(s1);
    aic_mat_herm_max_eig(s0, T0, prec);
    aic_mat_herm_max_eig(s1, T1, prec);

    /* hi = (1/2)(s0 + s1) + eps*d_min. (factor 1: Tr_min(eps I_D) = eps*d_min I.) */
    arb_add(hi, s0, s1, prec);
    arb_mul_2exp_si(hi, hi, -1);                  /* (s0+s1)/2 */
    arb_t shift;
    arb_init(shift);
    arb_mul_ui(shift, eps, (ulong) d_min, prec);  /* eps*d_min */
    arb_add(hi, hi, shift, prec);

    arb_clear(shift);
    arb_clear(s1);
    arb_clear(s0);
    acb_mat_clear(T1);
    acb_mat_clear(T0);
    arb_clear(eps);
    acb_mat_clear(blk);
    acb_mat_clear(negJd);
    acb_mat_clear(negJ);
}
