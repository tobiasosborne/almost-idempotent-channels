/* test_cstar.c — cross-checks for cstar_build Increment 1 (bead aic-097): the
 * S_P-as-aic_ecstar wrapper (compressed subalgebra) and the genuine M_d
 * matrix-algebra ecstar. These realize the recursion + codomain infrastructure of
 * the th_main master loop (approximate_algebras.tex:1321-1444, §9; compr_prod
 * .tex:1078-1082, M_n basis .tex:1371-1374). Built on the corner module
 * (aic_corner_*) and the ecstar axiom-defect estimators (aic_ecstar_*); the parent
 * eps-C* algebras come from the assoc module (aic_assoc_ecstar_from_phi) on the
 * exact-idempotent + oblique channels in tests/test_idemp.h.
 *
 * NO NEW double-vs-arb path is added here: the new code is STRUCTURAL GLUE over
 * already-tested numerical routines (aic_corner_Co/_extract/_Ptilde/_cdot carry
 * both number paths, the ecstar estimators are arb-certified). The cross-checks
 * are therefore (i) the eta=0 genuine-C* oracle (T1, T4: all five axiom defects
 * machine-zero); (ii) the eta=0 INDEPENDENT plain-product oracle (T1: wrapper star
 * == plain C_a.C_b, no Co_P in the RHS, so it catches a Co_P corruption — the F1
 * fix); (iii) the oblique eta>0 star teeth (T2 on mixconj: the S_P star defect is
 * O(eta), shrinking as eta->0, and the M1-RED is O(1) — the F3 fix); (iv) the
 * projection-recursion linchpin (T3: the nontrivial-projection finder runs ON the
 * S_P wrapper and returns an S_P-aware nontrivial projection); and (v) the
 * relocation/copy survival of the wrapper (T5: the heap star_ctx survives a
 * by-value copy + a move — the F2 fix). FINDINGS §C2 CAVEAT: the identity channel
 * has Phi_tilde=id so its star == plain product; T1's genuine-C* and plain-product
 * oracles are blind to a Phi_tilde bug — T2's oblique fixture gives the star teeth.
 *
 * MUTATION-PROVEN teeth (Rule 7; each RED observed by hand, then restored):
 *   T1(F1): scale s->ctx->Co_P by 2 (destroy idempotency) -> the INDEPENDENT
 *           plain-product oracle wrapper_vs_plain goes RED (~1e0) whereas the old
 *           wrapper_vs_cdot stays 0.0 (both sides read the same Co_P, the reviewer's
 *           "test that cannot fail"). Confirmed RED, restored.
 *   T2(M1): sp_star returning the plain XY (skipping Co_P + parent Phi_tilde) ->
 *           wrapper_vs_cdot goes RED on the OBLIQUE mixconj fixture at O(1) (~6e-1,
 *           the F3 fix — vs the old orthogonal dep+conj family's weak ~1.1e-2; the
 *           printed M1_RED column measures exactly this magnitude) AND the eta=0
 *           genuine-C* limit of the wrapper breaks. Confirmed RED, restored.
 *   T5(F2): revert star_ctx = out (the self-pointer into the owner) -> after the
 *           memset of the original bytes the copy/relocated star reads a zeroed
 *           parent pointer -> SIGSEGV / wrong star -> RED. Confirmed RED, restored.
 *   T4(M2): matalg_star using E_{ml} (transpose) instead of copying XY -> the
 *           star(X,Y)==plain X*Y check goes RED (E_{00}*E_{01} would give the
 *           wrong unit). Confirmed RED, restored.
 *   T1(M3): aic_cstar_subalg_build adopting the BOTTOM-dim_S singular vectors in
 *           aic_corner_extract -> dim_S still 4 but the genuine-C* defects of the
 *           wrapper blow up (wrong subspace). (Inherited corner tooth; documented.)
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_assoc.h"
#include "aic/aic_corner.h"
#include "aic/aic_cstar.h"
#include "aic/aic_ecstar.h"
#include "aic/aic_mat.h"
#include "aic/aic_projection.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"
#include "test_idemp.h"

#define CPREC 256

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Suppress -Wunused-function for the test_idemp.h builders this TU does not use.
 * (make_mixconj IS used now — T2/T5 — and make_compress_idemp/make_conjugated are
 * pulled in transitively by make_mixconj, so none of those is listed here.) */
__attribute__((unused))
static void suppress_unused_idemp_builders(void)
{
    (void) make_block_cond_exp;
    (void) make_trace_replace;
    (void) make_noiseless_subsystem;
    (void) make_dephasing;
}

/* ||A||_op on the BALL MIDPOINTS (mirrors test_corner/test_projection: a near-zero
 * difference matrix carrying wide accumulated balls false-fails the certified Gram
 * check, aic-2yo/aic-qgs; collapse to midpoints when we only want the number). */
static double mid_opnorm_d(const acb_mat_t Aop, slong prec)
{
    slong r = acb_mat_nrows(Aop), c = acb_mat_ncols(Aop);
    acb_mat_t M;
    arb_t e;
    acb_mat_init(M, r, c);
    arb_init(e);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(Aop, i, j));
    aic_mat_opnorm(e, M, prec);
    double v = dd(e);
    arb_clear(e);
    acb_mat_clear(M);
    return v;
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op (superop idempotence defect, ~ ||Phi^2-Phi||
 * <= eta; mirrors test_corner/test_projection eta_proxy). */
static double eta_proxy(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S, S2, D;
    acb_mat_init(S, nn, nn);
    acb_mat_init(S2, nn, nn);
    acb_mat_init(D, nn, nn);
    aic_assoc_superop_from_ucp(S, phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    double o = mid_opnorm_d(D, prec);
    acb_mat_clear(D); acb_mat_clear(S2); acb_mat_clear(S);
    return o;
}

/* sigma_max(S_tilde) = ||theta(2 Phi - 1)||_op (the regularized idempotent
 * superop's largest singular value). For an HS-SELF-ADJOINT (orthogonal) Phi_tilde
 * sigma_max == 1; for a GENUINELY OBLIQUE one it is > 1 (FINDINGS §C4): the visible
 * obliqueness witness. For mixconj(4,2) this reads ~sqrt(3)=1.732 (test_idemp.h /
 * aic_assoc.h note), vs 1.0 for the orthogonal dep+conj family. */
static double sigma_max_stilde(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S;
    acb_mat_init(S, nn, nn);
    aic_assoc_regularize(S, phi, prec);
    double s = mid_opnorm_d(S, prec);
    acb_mat_clear(S);
    return s;
}

/* n x n diagonal projector with the first `k` diagonal entries 1 (rest 0). */
static void diag_proj(acb_mat_t P, slong n, slong k)
{
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    for (slong i = 0; i < k; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
}

/* n x n diagonal projector onto span(e_lo .. e_{hi-1}) (entries [lo,hi) set 1).
 * For mixconj's carrier M = span(e0..e_{m-1}), a SUB-projector inside M minus e0
 * (e.g. span(e1,e2)) is a genuine in-A delta-projection (in-A residual O(eta)),
 * whereas one touching e0 or M^perp is genuinely oblique (residual O(1)). */
static void diag_proj_range(acb_mat_t P, slong n, slong lo, slong hi)
{
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    for (slong i = lo; i < hi; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
}

/* Max over basis pairs of ||wrapper_star(C_a,C_b) - aic_corner_cdot(parent,Co_P,
 * C_a,C_b)||_op: the thunk-vs-cdot consistency. SECONDARY M1-localized check only:
 * it DOES catch sp_star->plainXY (cdot keeps the star, so plainXY != cdot), but it
 * is BLIND to a Co_P CORRUPTION (FINDINGS F1, the reviewer's proof): both routes
 * read the SAME s->ctx->Co_P, so scaling Co_P leaves the diff at 0.0. The headline
 * correctness claim at eta=0 is therefore wrapper_vs_plain (below), not this. */
static double wrapper_vs_cdot(const aic_cstar_subalgebra *s, slong prec)
{
    slong n = s->parent->n, ds = s->A.dim_A;
    acb_mat_t sw, cd, db;
    acb_mat_init(sw, n, n); acb_mat_init(cd, n, n); acb_mat_init(db, n, n);
    double worst = 0.0;
    for (slong a = 0; a < ds; a++)
        for (slong b = 0; b < ds; b++) {
            aic_ecstar_star(sw, &s->A, s->A.B[a], s->A.B[b], prec);     /* wrapper * */
            aic_corner_cdot(cd, s->parent, s->ctx->Co_P, s->A.B[a], s->A.B[b], prec);
            acb_mat_sub(db, sw, cd, prec);
            double e = mid_opnorm_d(db, prec);
            if (e > worst) worst = e;
        }
    acb_mat_clear(db); acb_mat_clear(cd); acb_mat_clear(sw);
    return worst;
}

/* INDEPENDENT eta=0 oracle (the F1 fix): max over basis pairs of
 *   ||wrapper_star(C_a, C_b) - C_a . C_b||_op   (PLAIN acb_mat_mul on the RHS).
 * For the identity channel (Phi_tilde = id) and an ORTHOGONAL projector P (here
 * P = diag(1,1,0), so Co_P is the orthogonal compression and {C_m} spans
 * S_P = P M_n P), the compressed product reduces to the ORDINARY matrix product:
 * for C_a, C_b in P M_n P, C_a C_b is already in P M_n P and Co_P is the identity
 * on S_P, so Co_P(Phi_tilde(C_a C_b)) = C_a C_b. The RHS contains NO Co_P, so this
 * oracle is GENUINELY INDEPENDENT of Co_P and CATCHES a Co_P corruption (unlike
 * wrapper_vs_cdot, which reads the same Co_P on both sides). Mutation-proof: scale
 * s->ctx->Co_P by 2 -> this oracle goes RED (~1e0), wrapper_vs_cdot stays ~0.
 * ONLY valid when Phi_tilde = id AND P is an orthogonal projector (T1). */
static double wrapper_vs_plain(const aic_cstar_subalgebra *s, slong prec)
{
    slong n = s->parent->n, ds = s->A.dim_A;
    acb_mat_t sw, pl, db;
    acb_mat_init(sw, n, n); acb_mat_init(pl, n, n); acb_mat_init(db, n, n);
    double worst = 0.0;
    for (slong a = 0; a < ds; a++)
        for (slong b = 0; b < ds; b++) {
            aic_ecstar_star(sw, &s->A, s->A.B[a], s->A.B[b], prec);   /* wrapper * */
            acb_mat_mul(pl, s->A.B[a], s->A.B[b], prec);             /* C_a . C_b  */
            acb_mat_sub(db, sw, pl, prec);
            double e = mid_opnorm_d(db, prec);
            if (e > worst) worst = e;
        }
    acb_mat_clear(db); acb_mat_clear(pl); acb_mat_clear(sw);
    return worst;
}

/* Max of the unit-INDEPENDENT eps-C* axiom defects (assoc, submult, C*,
 * involution) of an ecstar A. EXCLUDES aic_ecstar_defect_unit because that
 * estimator is hardcoded to test the unit laws against the AMBIENT 1_n (its
 * docstring: it checks "1_n in A" and "B_k * I_n"), which is the WRONG unit for a
 * compressed subalgebra S_P (whose unit is Ptilde = Co_P(P), not 1_n; .tex:1082
 * and CLAUDE.md "the exact unit is only available after prop_unit"). The S_P unit
 * law is tested separately via ptilde_unit_defect. See FINDINGS §C7. */
static double max_axiom_defect_no_unit(const aic_ecstar *A, slong prec)
{
    arb_t d;
    arb_init(d);
    double m = 0.0, v;
    aic_ecstar_defect_assoc(d, A, prec);      v = dd(d); if (v > m) m = v;
    aic_ecstar_defect_submult(d, A, prec);    v = dd(d); if (v > m) m = v;
    aic_ecstar_defect_cstar(d, A, prec);      v = dd(d); if (v > m) m = v;
    aic_ecstar_defect_involution(d, A, prec); v = dd(d); if (v > m) m = v;
    arb_clear(d);
    return m;
}

/* Max of all FIVE eps-C* axiom defects (incl. the ambient-unit law) — correct
 * ONLY when the algebra's unit IS the ambient 1_n (a genuine M_d, T4). */
static double max_axiom_defect(const aic_ecstar *A, slong prec)
{
    arb_t d;
    arb_init(d);
    double m = max_axiom_defect_no_unit(A, prec), v;
    aic_ecstar_defect_unit(d, A, prec);       v = dd(d); if (v > m) m = v;
    arb_clear(d);
    return m;
}

/* The S_P unit law against Ptilde (the genuine unit of the compressed subalgebra,
 * .tex:1082): max over the basis {C_m} of max(||Ptilde * C_m - C_m||,
 * ||C_m * Ptilde - C_m||)_op. ~0 at eta=0, O(delta+eps) at eta>0. This replaces
 * aic_ecstar_defect_unit's ambient-1_n test for the wrapper. */
static double ptilde_unit_defect(const aic_cstar_subalgebra *s, slong prec)
{
    slong n = s->parent->n, ds = s->A.dim_A;
    acb_mat_t lu, ru, db;
    acb_mat_init(lu, n, n); acb_mat_init(ru, n, n); acb_mat_init(db, n, n);
    double worst = 0.0;
    for (slong m = 0; m < ds; m++) {
        aic_ecstar_star(lu, &s->A, s->Ptilde, s->A.B[m], prec);  /* Ptilde * C_m */
        aic_ecstar_star(ru, &s->A, s->A.B[m], s->Ptilde, prec);  /* C_m * Ptilde */
        acb_mat_sub(db, lu, s->A.B[m], prec);
        double le = mid_opnorm_d(db, prec);
        acb_mat_sub(db, ru, s->A.B[m], prec);
        double re = mid_opnorm_d(db, prec);
        if (le > worst) worst = le;
        if (re > worst) worst = re;
    }
    acb_mat_clear(db); acb_mat_clear(ru); acb_mat_clear(lu);
    return worst;
}

/* ===================== T1: eta=0 genuine-C* oracle for the S_P wrapper ======= *
 * Identity channel Phi=id on C^3 (eta=0; A = all of M_3). P = diag(1,1,0) (rank-2
 * ambient projector in A). The S_P wrapper -> S_P ~= M_2 (dim 4). Assert:
 *  (a) dim_S == 4 (= rank(P)^2 for the genuine M_2 corner);
 *  (b) HEADLINE CORRECTNESS (F1 fix): the wrapper's star == the PLAIN matrix
 *      product C_a . C_b on basis pairs (wrapper_vs_plain). Independent of Co_P,
 *      so it CATCHES a Co_P corruption (proven below in test_t5_copy_survives's
 *      Co_P-scale mutation; the old wrapper_vs_cdot stays 0.0 there).
 *  (c) (secondary) wrapper star == aic_corner_cdot (M1-localized; catches
 *      sp_star->plainXY but NOT a Co_P corruption);
 *  (d) the wrapper is a GENUINE C* algebra: all five axiom defects ~ machine-zero.
 * CAVEAT (FINDINGS §C2): Phi=id => star==plain, so T1 is the genuine-C* oracle
 * ONLY; it is BLIND to a star bug. T2's oblique fixture supplies the star teeth. */
static void test_t1_identity_oracle(void)
{
    printf("T1 eta=0 genuine-C* oracle (identity(3), P=diag(1,1,0), S_P ~= M_2):\n");
    aic_ucp_kraus phi;
    make_identity(&phi, 3);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    acb_mat_t P;
    diag_proj(P, 3, 2);
    aic_cstar_subalgebra s;
    aic_cstar_subalg_build(&s, &ae.A, P, CPREC);

    printf("  dim_S=%ld (expect 4)\n", (long) s.A.dim_A);
    AIC_CHECK_MSG(s.A.dim_A == 4, "T1: dim S_P=%ld, expected 4 (M_2)",
                  (long) s.A.dim_A);

    /* (b) HEADLINE: independent plain-product oracle (no Co_P in the RHS). */
    double vp = wrapper_vs_plain(&s, CPREC);
    printf("  max||wrapper_star - plain C_a.C_b||=%.3e (independent oracle)\n", vp);
    AIC_CHECK_MSG(vp < 1e-12, "T1: wrapper star != plain product (%.3e)", vp);

    /* (c) secondary thunk-vs-cdot consistency (M1-localized). */
    double vc = wrapper_vs_cdot(&s, CPREC);
    printf("  max||wrapper_star - corner_cdot||=%.3e (secondary)\n", vc);
    AIC_CHECK_MSG(vc < 1e-12, "T1: wrapper star != aic_corner_cdot (%.3e)", vc);

    /* (d) Genuine-C* oracle: the four unit-INDEPENDENT axiom defects machine-zero,
     * AND the S_P unit law against Ptilde machine-zero. (aic_ecstar_defect_unit's
     * ambient-1_n test reads 1.0 here because S_P's unit is Ptilde, not 1_n —
     * FINDINGS §C7; the correct unit test is ptilde_unit_defect.) */
    double md = max_axiom_defect_no_unit(&s.A, CPREC);
    double ud = ptilde_unit_defect(&s, CPREC);
    printf("  max(assoc,submult,C*,involution) defect=%.3e, "
           "Ptilde unit-law defect=%.3e (expect machine-zero)\n", md, ud);
    AIC_CHECK_MSG(md < 1e-12, "T1: S_P wrapper not a genuine C* algebra (%.3e)", md);
    AIC_CHECK_MSG(ud < 1e-12, "T1: Ptilde not the genuine unit of S_P (%.3e)", ud);

    acb_mat_clear(P);
    aic_cstar_subalg_clear(&s);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* max over basis pairs of ||plain C_a C_b - aic_corner_cdot(C_a,C_b)||_op: the
 * EXACT magnitude the wrapper_vs_cdot check would read under MUTATION M1 (sp_star
 * returning the plain XY instead of Co_P(Phi_tilde(XY))). Printed in T2 so the
 * M1 RED is visible (and confirmed O(1) on the oblique mixconj fixture). The plain
 * product leaves S_P / drops Phi_tilde, so on a GENUINELY OBLIQUE channel this is
 * O(1) — the F3 teeth (vs the orthogonal dep+conj family's ~1.1e-2). */
static double m1_red_plainxy_vs_cdot(const aic_cstar_subalgebra *s, slong prec)
{
    slong n = s->parent->n, ds = s->A.dim_A;
    acb_mat_t pl, cd, db;
    acb_mat_init(pl, n, n); acb_mat_init(cd, n, n); acb_mat_init(db, n, n);
    double worst = 0.0;
    for (slong a = 0; a < ds; a++)
        for (slong b = 0; b < ds; b++) {
            acb_mat_mul(pl, s->A.B[a], s->A.B[b], prec);     /* plain (the M1 bug) */
            aic_corner_cdot(cd, s->parent, s->ctx->Co_P, s->A.B[a], s->A.B[b], prec);
            acb_mat_sub(db, pl, cd, prec);
            double e = mid_opnorm_d(db, prec);
            if (e > worst) worst = e;
        }
    acb_mat_clear(db); acb_mat_clear(cd); acb_mat_clear(pl);
    return worst;
}

/* ===================== T2: oblique eta>0 star teeth ========================= *
 * make_mixconj(5, 3, t): a GENUINELY OBLIQUE eta>0 channel — sigma_max(S_tilde) ~
 * sqrt(3) = 1.732 (FINDINGS §C4; an orthogonal/HS-self-adjoint Phi_tilde would read
 * sigma_max == 1). The old dep+conj eta-family was ORTHOGONAL (sigma_max == 1), so
 * the F3 M1 mutation only gave a weak ~1.1e-2 RED; this oblique fixture supplies an
 * O(1) (~5e-1) RED (the F3 fix — the teeth bite). The carrier of compress_idemp(5,3)
 * is M = span(e0,e1,e2), so the choice of P matters (the corner-T8(c) finding,
 * FINDINGS: the S_P unit law ||Ptilde.X - X|| is O(delta+eps) ONLY when P is a
 * genuine delta-PROJECTION IN A; an ambient P far from A has delta = O(1) and the
 * bound is VACUOUS). We therefore split the two concerns across two P's of the SAME
 * oblique channel:
 *  - P_good = span(e1,e2) = diag(0,1,1,0,0): a GENUINE in-A rank-2 delta-projection
 *    (in-A residual O(eta), measured 2.0e-2 -> 6.6e-3 as t: 0.06 -> 0.02), dim S_P
 *    = 4 (an M_2 corner). Used for (a) thunk-vs-cdot consistency and (b) the O(eta)-
 *    SHRINKING S_P unit defect ||Ptilde * Ptilde - Ptilde|| (7.0e-4 -> 7.4e-5).
 *  - P_obl  = span(e0,e1,e2) = diag(1,1,1,0,0): the full carrier, GENUINELY OBLIQUE
 *    (in-A residual O(1) ~0.67), dim S_P = 9. Used for (c) the M1 O(1)-RED teeth
 *    (plain XY vs cdot ~ 0.43; the star path is genuinely non-plain here).
 * Print sigma_max(S_tilde) + the in-A residual of each P so the obliqueness and the
 * good/oblique split are visible. MUTATION M1 (the star teeth): sp_star returning
 * plain XY -> (c)'s plainXY==cdot equality breaks at O(1) on P_obl, AND (b)'s
 * eta=0 genuine-C* limit breaks. */
static void test_t2_oblique_teeth(void)
{
    printf("T2 oblique eta>0 star teeth (mixconj(5,3); P_good=span(e1,e2), "
           "P_obl=span(e0,e1,e2)):\n");
    const double ts[] = {0.06, 0.02};
    double defs[2], etas[2];
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phi;
        make_mixconj(&phi, 5, 3, ts[it], CPREC);
        double eta = eta_proxy(&phi, CPREC);
        double smax = sigma_max_stilde(&phi, CPREC);
        etas[it] = eta;
        AIC_CHECK_MSG(smax > 1.3,
                      "T2: sigma_max(S_tilde)=%.3f not oblique (>1.3) at t=%.2f",
                      smax, ts[it]);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

        /* ---- P_good: the genuine in-A rank-2 delta-projection (a, b) ---- */
        acb_mat_t Pg;
        diag_proj_range(Pg, 5, 1, 3);                /* span(e1,e2) */
        arb_t rg;
        arb_init(rg);
        aic_ecstar_proj_residual(rg, &ae.A, Pg, CPREC);
        double resg = dd(rg);
        arb_clear(rg);
        aic_cstar_subalgebra sg;
        aic_cstar_subalg_build(&sg, &ae.A, Pg, CPREC);
        AIC_CHECK_MSG(sg.A.dim_A >= 2,
                      "T2: dim S_P_good=%ld < 2 (need a >1-dim corner)",
                      (long) sg.A.dim_A);

        /* (a) thunk-vs-cdot consistency (machine-zero at any eta). */
        double vc = wrapper_vs_cdot(&sg, CPREC);
        AIC_CHECK_MSG(vc < 1e-12,
                      "T2(a): wrapper star != aic_corner_cdot (%.3e) at t=%.2f",
                      vc, ts[it]);

        /* (b) S_P unit defect of Ptilde: ||Ptilde * Ptilde - Ptilde|| -> O(eta),
         * SHRINKS as eta->0 (genuine because P_good is a real in-A delta-proj). */
        slong n = 5;
        acb_mat_t pp, db;
        acb_mat_init(pp, n, n); acb_mat_init(db, n, n);
        aic_ecstar_star(pp, &sg.A, sg.Ptilde, sg.Ptilde, CPREC);  /* Ptilde * Ptilde */
        acb_mat_sub(db, pp, sg.Ptilde, CPREC);
        double sdef = mid_opnorm_d(db, CPREC);
        defs[it] = sdef;
        double c = eta > 1e-14 ? sdef / eta : 0.0;
        acb_mat_clear(db); acb_mat_clear(pp);

        /* ---- P_obl: the genuinely oblique full-carrier projection (c) ---- */
        acb_mat_t Po;
        diag_proj_range(Po, 5, 0, 3);                /* span(e0,e1,e2) */
        arb_t ro;
        arb_init(ro);
        aic_ecstar_proj_residual(ro, &ae.A, Po, CPREC);
        double reso = dd(ro);
        arb_clear(ro);
        aic_cstar_subalgebra so;
        aic_cstar_subalg_build(&so, &ae.A, Po, CPREC);
        AIC_CHECK_MSG(so.A.dim_A >= 2, "T2: dim S_P_obl=%ld < 2", (long) so.A.dim_A);
        /* thunk-vs-cdot is machine-zero here too (the wrapper star IS the cdot). */
        double vco = wrapper_vs_cdot(&so, CPREC);
        AIC_CHECK_MSG(vco < 1e-12, "T2: P_obl wrapper star != cdot (%.3e)", vco);
        /* (c) the M1-RED magnitude (plain XY vs cdot): O(1) on this oblique corner.
         * THIS is the F3 tooth — dropping the star to plain XY corrupts the product
         * by O(1) (~0.4), not the weak ~1.1e-2 of the old orthogonal fixture. */
        double m1 = m1_red_plainxy_vs_cdot(&so, CPREC);
        AIC_CHECK_MSG(m1 > 0.1,
                      "T2(c): M1 RED (plainXY vs cdot)=%.3e not O(1) — P_obl not "
                      "oblique enough (t=%.2f)", m1, ts[it]);

        printf("  t=%.2f eta=%.4e sigma_max(S~)=%.4f | P_good: res=%.3e dim_S=%ld "
               "||wrap-cdot||=%.3e ||Ptil*Ptil-Ptil||=%.3e (c=%.3f) | P_obl: res=%.3e "
               "dim_S=%ld M1_RED=%.3e\n",
               ts[it], eta, smax, resg, (long) sg.A.dim_A, vc, sdef, c, reso,
               (long) so.A.dim_A, m1);
        AIC_CHECK_MSG(c < 5.0,
                      "T2(b): S_P unit star defect / eta = %.3f exceeds 5 (t=%.2f)",
                      c, ts[it]);

        aic_cstar_subalg_clear(&so);
        acb_mat_clear(Po);
        aic_cstar_subalg_clear(&sg);
        acb_mat_clear(Pg);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    /* The unit defect (on the genuine in-A P_good) shrinks with eta (the O(eta)
     * bound is genuine, not vacuous). */
    printf("  -> P_good unit defect shrinks: t=0.06 def=%.3e (eta=%.3e), "
           "t=0.02 def=%.3e (eta=%.3e)\n", defs[0], etas[0], defs[1], etas[1]);
    AIC_CHECK_MSG(defs[1] < defs[0],
                  "T2: S_P unit defect did not shrink as eta->0 (%.3e !< %.3e)",
                  defs[1], defs[0]);
}

/* ===================== T3: projection-recursion linchpin =================== *
 * On the T1 identity-channel A = M_3 with P = diag(1,1,0) (S_P ~= M_2, dim 4 > 1),
 * run aic_projection_nontrivial ON the S_P WRAPPER. Assert it returns a nontrivial
 * projection P' with:
 *  - star defect (via the WRAPPER's star) ~ 0 (eta=0 => exact);
 *  - ||P'||_op bounded away from 0;
 *  - ||Ptilde - P'||_op >= 0.3 (the S_P-AWARE nontriviality: the complement in S_P
 *    is Ptilde - P', NOT 1_n - P'; design §2.3 / §G1).
 * This proves the projection finder recurses on the wrapper. */
static void test_t3_projection_recursion(void)
{
    printf("T3 projection recursion (nontrivial proj IN the S_P=M_2 wrapper):\n");
    aic_ucp_kraus phi;
    make_identity(&phi, 3);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    acb_mat_t P;
    diag_proj(P, 3, 2);
    aic_cstar_subalgebra s;
    aic_cstar_subalg_build(&s, &ae.A, P, CPREC);
    AIC_CHECK_MSG(s.A.dim_A > 1, "T3: need dim S_P > 1, got %ld", (long) s.A.dim_A);

    slong n = 3;
    acb_mat_t Pp;
    arb_t delta, Pnorm, ImPnorm;
    acb_mat_init(Pp, n, n);
    arb_init(delta); arb_init(Pnorm); arb_init(ImPnorm);
    /* The finder recurses on the wrapper's ecstar (its basis = the S_P corner
     * basis, its star = the compressed product).
     * F5 (design §G1, FINDINGS §C7): aic_projection_nontrivial's INTERNAL
     * nontriviality assert is hardcoded to the AMBIENT 1_n (it checks ||1_n - P'||),
     * which is VACUOUS on an S_P wrapper (1_n is NOT the unit of S_P; the S_P unit
     * is Ptilde = Co_P(P)). So the wrapper's internal ImPnorm is NOT the S_P
     * complement norm — the S_P-AWARE nontriviality (||Ptilde - P'||) MUST be
     * checked EXTERNALLY here, which is exactly what the (b) leg below does. */
    aic_projection_nontrivial(Pp, delta, Pnorm, ImPnorm, &s.A, NULL, CPREC);

    double del = dd(delta), pn = dd(Pnorm);
    printf("  star defect ||P'*P'-P'||=%.3e (wrapper star), ||P'||=%.3f\n",
           del, pn);
    AIC_CHECK_MSG(del < 1e-9, "T3: S_P projection star defect not ~0 (%.3e)", del);
    /* F4 (design §2.3 / §G1): nontriviality IN S_P needs BOTH legs JOINTLY —
     * (a) ||P'|| > 0.3 (P' bounded away from 0, not the zero projection) AND
     * (b) ||Ptilde - P'|| > 0.3 (P' bounded away from Ptilde, the S_P unit / "full"
     *     projection). Either leg alone is insufficient: a projection equal to
     *     Ptilde passes (a) but is trivial; the zero projection passes (b) but is
     *     trivial. I5's Stage-1 greedy loop MUST keep both. */
    AIC_CHECK_MSG(pn > 0.3, "T3: ||P'|| not bounded away from 0 (%.3f)", pn);

    /* (b) S_P-aware nontriviality: ||Ptilde - P'||_op >= 0.3 (the complement in S_P
     * is Ptilde - P', NOT the ambient 1_n - P'). The F5 external check. */
    acb_mat_t comp;
    acb_mat_init(comp, n, n);
    acb_mat_sub(comp, s.Ptilde, Pp, CPREC);
    double cn = mid_opnorm_d(comp, CPREC);
    printf("  ||Ptilde - P'||=%.3f (S_P-aware nontriviality, expect >= 0.3)\n", cn);
    AIC_CHECK_MSG(cn >= 0.3,
                  "T3: ||Ptilde - P'|| = %.3f < 0.3 (P' trivial IN S_P)", cn);
    acb_mat_clear(comp);

    arb_clear(ImPnorm); arb_clear(Pnorm); arb_clear(delta);
    acb_mat_clear(Pp);
    acb_mat_clear(P);
    aic_cstar_subalg_clear(&s);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T4: genuine M_d matrix algebra ====================== *
 * aic_cstar_matrix_algebra(M_3): assert dim_A == 9, basis Frobenius-orthonormal,
 * ALL five ecstar axiom defects machine-zero (a genuine C* algebra), and the star
 * star(X,Y) == plain X*Y on a couple of basis pairs (E_{lm}E_{m'n} = delta_{m,m'}
 * E_{ln}, .tex:1373).
 * MUTATION M2 (the star teeth): matalg_star transposing -> the star==plain check
 * RED. */
static void test_t4_matrix_algebra(void)
{
    printf("T4 genuine M_d matrix algebra (M_3):\n");
    const slong d = 3;
    aic_ecstar M;
    aic_cstar_matrix_algebra(&M, d, CPREC);

    AIC_CHECK_MSG(M.dim_A == d * d, "T4: dim_A=%ld, expected %ld", (long) M.dim_A,
                  (long) (d * d));

    /* Frobenius-orthonormality: <E_a,E_b>_F = delta_ab (each unit is a single 1). */
    double worst_orth = 0.0;
    for (slong a = 0; a < M.dim_A; a++)
        for (slong b = 0; b < M.dim_A; b++) {
            acb_t ip;
            acb_init(ip);
            /* <E_a,E_b>_F = sum_ij conj(E_a[i,j]) E_b[i,j] = Tr(E_a^dag E_b). */
            acb_mat_t Ead, G;
            acb_mat_init(Ead, d, d);
            acb_mat_init(G, d, d);
            acb_mat_conjugate_transpose(Ead, M.B[a]);
            acb_mat_mul(G, Ead, M.B[b], CPREC);
            acb_mat_trace(ip, G, CPREC);
            arb_t mag;
            arb_init(mag);
            acb_sub_si(ip, ip, a == b ? 1 : 0, CPREC);
            acb_abs(mag, ip, CPREC);
            double e = dd(mag);
            if (e > worst_orth) worst_orth = e;
            arb_clear(mag);
            acb_mat_clear(G); acb_mat_clear(Ead);
            acb_clear(ip);
        }
    printf("  max|<E_a,E_b>_F - delta_ab|=%.3e\n", worst_orth);
    AIC_CHECK_MSG(worst_orth < 1e-14, "T4: matrix units not Frobenius-orthonormal "
                  "(%.3e)", worst_orth);

    double md = max_axiom_defect(&M, CPREC);
    printf("  max eps-C* axiom defect=%.3e (expect machine-zero)\n", md);
    AIC_CHECK_MSG(md < 1e-12, "T4: M_d not a genuine C* algebra (%.3e)", md);

    /* star(E_{lm}, E_{m'n}) == plain E_{lm} E_{m'n}. Check E_{01}*E_{12} = E_{02}
     * (index l*d+m): E_{01}=B[1], E_{12}=B[5], E_{02}=B[2]. And E_{01}*E_{00}=0. */
    {
        acb_mat_t prod, plain, db;
        acb_mat_init(prod, d, d); acb_mat_init(plain, d, d); acb_mat_init(db, d, d);
        /* pair 1: E_{01} * E_{12} (should be E_{02}). */
        aic_ecstar_star(prod, &M, M.B[0 * d + 1], M.B[1 * d + 2], CPREC);
        acb_mat_mul(plain, M.B[0 * d + 1], M.B[1 * d + 2], CPREC);
        acb_mat_sub(db, prod, plain, CPREC);
        double e1 = mid_opnorm_d(db, CPREC);
        acb_mat_sub(db, prod, M.B[0 * d + 2], CPREC);   /* == E_{02}? */
        double e1u = mid_opnorm_d(db, CPREC);
        /* pair 2: E_{01} * E_{00} (delta_{1,0}=0 -> zero). */
        aic_ecstar_star(prod, &M, M.B[0 * d + 1], M.B[0 * d + 0], CPREC);
        acb_mat_mul(plain, M.B[0 * d + 1], M.B[0 * d + 0], CPREC);
        acb_mat_sub(db, prod, plain, CPREC);
        double e2 = mid_opnorm_d(db, CPREC);
        printf("  star==plain: E01*E12 diff=%.3e (==E02 diff=%.3e), "
               "E01*E00 diff=%.3e\n", e1, e1u, e2);
        AIC_CHECK_MSG(e1 < 1e-14 && e2 < 1e-14,
                      "T4: star != plain product (%.3e, %.3e)", e1, e2);
        AIC_CHECK_MSG(e1u < 1e-14, "T4: E01*E12 != E02 (%.3e)", e1u);
        acb_mat_clear(db); acb_mat_clear(plain); acb_mat_clear(prod);
    }

    aic_ecstar_clear(&M);
}

/* ===================== T5: wrapper survives copy / relocation (F2) ========== *
 * The F2 fix: aic_cstar_subalgebra is RELOCATION-SAFE and SHALLOW-COPYABLE because
 * its derived ecstar's star_ctx points at a SEPARATELY HEAP-ALLOCATED ctx block,
 * NOT into the owner struct. The I5 master loop stores many wrappers in a
 * relocatable array, so a wrapper must keep computing the correct star after being
 * copied by value and after the original storage is moved/overwritten.
 *
 * On the OBLIQUE mixconj(4,2) fixture (Phi_tilde != id, so the star genuinely
 * needs ctx->Co_P + ctx->parent; an eta=0 identity fixture would be BLIND here):
 *  1. build wrapper s; capture a reference star ref = corner_cdot(C_0, C_1)
 *     computed directly from the parent + Co_P (independent of the wrapper's seam);
 *  2. by-value copy scopy = s;
 *  3. memcpy s into a heap array slot, then ZERO the original s bytes (simulating
 *     the original stack/array slot being reused after a MOVE);
 *  4. call the wrapper star through BOTH the copy and the relocated slot and assert
 *     each equals ref to machine precision.
 * With the OLD star_ctx = &owner self-pointer, steps 3-4 read the zeroed original
 * (garbage parent pointer) -> wrong star or crash: this test would be RED. With the
 * heap ctx it is unaffected. Clear EXACTLY ONCE (via the slot; move-only). */
static void test_t5_copy_survives(void)
{
    printf("T5 wrapper survives copy/relocation (F2, oblique mixconj(4,2)):\n");
    aic_ucp_kraus phi;
    make_mixconj(&phi, 4, 2, 0.02, CPREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);

    slong n = 4;
    acb_mat_t P;
    diag_proj(P, 4, 2);
    aic_cstar_subalgebra s;
    aic_cstar_subalg_build(&s, &ae.A, P, CPREC);
    AIC_CHECK_MSG(s.A.dim_A >= 2, "T5: dim S_P=%ld < 2", (long) s.A.dim_A);

    /* (1) reference star ref = C_0 . C_1 via the parent + Co_P directly. */
    acb_mat_t ref, got, db;
    acb_mat_init(ref, n, n); acb_mat_init(got, n, n); acb_mat_init(db, n, n);
    aic_corner_cdot(ref, s.parent, s.ctx->Co_P, s.A.B[0], s.A.B[1], CPREC);

    /* (2) by-value copy. */
    aic_cstar_subalgebra scopy = s;

    /* (3) relocate: memcpy s into a heap slot, then zero the original bytes. */
    aic_cstar_subalgebra *slot = flint_malloc(sizeof(aic_cstar_subalgebra));
    memcpy(slot, &s, sizeof(aic_cstar_subalgebra));
    memset(&s, 0, sizeof(aic_cstar_subalgebra));   /* original slot now garbage */

    /* (4) star through the COPY (basis pointers shared, still valid). */
    aic_ecstar_star(got, &scopy.A, scopy.A.B[0], scopy.A.B[1], CPREC);
    acb_mat_sub(db, got, ref, CPREC);
    double ecopy = mid_opnorm_d(db, CPREC);
    /* star through the RELOCATED slot. */
    aic_ecstar_star(got, &slot->A, slot->A.B[0], slot->A.B[1], CPREC);
    acb_mat_sub(db, got, ref, CPREC);
    double eslot = mid_opnorm_d(db, CPREC);
    printf("  ||copy_star - ref||=%.3e, ||relocated_star - ref||=%.3e "
           "(expect machine-zero)\n", ecopy, eslot);
    AIC_CHECK_MSG(ecopy < 1e-12, "T5: copy's star != ref (%.3e); star_ctx dangled?",
                  ecopy);
    AIC_CHECK_MSG(eslot < 1e-12, "T5: relocated star != ref (%.3e); star_ctx dangled?",
                  eslot);

    acb_mat_clear(db); acb_mat_clear(got); acb_mat_clear(ref);
    acb_mat_clear(P);
    /* MOVE-ONLY: clear the live wrapper EXACTLY ONCE (via the slot). scopy shares
     * the same ctx/Ptilde/basis -> must NOT be cleared (would double-free). */
    aic_cstar_subalg_clear(slot);
    flint_free(slot);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_t1_identity_oracle();
    test_t2_oblique_teeth();
    test_t3_projection_recursion();
    test_t4_matrix_algebra();
    test_t5_copy_survives();
    aic_test_report("test_cstar");
    return 0;
}
