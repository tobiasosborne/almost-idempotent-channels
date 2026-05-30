/* test_dhom.c — cross-checks for the dhom module (bead aic-c1n): the
 * delta-homomorphism approximation of §8 (approximate_algebras.tex:1190-1319,
 * prop_delta_hominc + lem_approx). FULLY CONSTRUCTIVE (no eig): the tests are
 * pure products + certified operator norms.
 *
 * Cross-check ladder (CLAUDE.md Rule 6/7; mutation-proven, see per-test comments):
 *   T1 Pauli unit tests — the diagonal is the dimension-independence crux. Each
 *      S_{jk} unitary; pi(D_l)=I; centrality X D_l = D_l X (random X, d=2..4);
 *      ||D_l||_proj = 1. MUTATION: drop the clock phase exp(2 pi i k m/d) (set
 *      k=0) so {S_{jk}} is no longer an operator basis -> centrality RED (defect
 *      0.48 at d=2, abort). A mere index relabel (shift j or a global phase per
 *      S) leaves the FULL sum invariant and is correctly NOT caught; only a
 *      basis-breaking change moves the diagonal. Verified RED 2026-05-30.
 *   T2 eta=0 ORACLE (headline) — A=M_2/M_3 identity channel, v(X)=U X U^dag an
 *      exact *-isomorphism -> G_v = 0 (machine-zero, all basis pairs); lem_approx
 *      returns v unchanged (w'=w''=0). NOTE: identity star == plain product, so
 *      T2 is BLIND to the star bug (only T4 eta>0 catches it).
 *   T3 delta>0, eps=0 (Johnson) — perturb an exact hom by delta*R -> lem_approx
 *      drives the defect to machine-zero (exact hom limit) and ||vtilde - v_pert||
 *      <= O(delta). MUTATION: 0 Newton steps -> defect = delta >> 0 -> RED.
 *   T4 eps>0 GENUINE almost-idempotent A (make_mixconj) — lem_approx reduces the
 *      defect to O(eps); report C = defect/eps. STAR NON-VACUITY (the corner
 *      lesson): MUTATE G_v's star to plain matmul -> bound violated at eps>0 -> RED.
 *   T5 universality canary (.tex:484/460) — defect/eps constant must NOT grow with
 *      EITHER the block dim OR the number of blocks. (A) block-dim sweep
 *      B=M_2,M_3,M_4 (ratio < 2). (B) number-of-blocks sweep m=1,2,3 of M_2
 *      (C does not grow with m; the FIX-3 escalation check — the embedded-sum
 *      diagonal has ||D||_proj = m, so the w' bound is O(m*delta), but the achieved
 *      C stays bounded and is largest at m=1). MUTATION: replace the Pauli diagonal
 *      with a non-diagonal -> canary RED. Proves the Pauli diagonal load-bearing.
 *   T6 prop_delta_hominc bounds — exact hom: ||v||=1, a=1, unit exact; perturbed:
 *      bounds hold.
 */
#include <complex.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_mat.h"
#include "aic_test.h"
#include "aic_ucp.h"
#include "test_idemp.h"

/* test_idemp.h defines PREC 53 for its eta=0 harness; dhom wants a higher working
 * precision for the certified defect balls, so use a distinct name. */
#undef PREC
#define PREC 80
static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Suppress -Wunused-function for the test_idemp.h channel builders this TU does
 * not exercise. References the addresses so the file is warning-clean. */
__attribute__((unused))
static void suppress_unused_idemp_builders(void)
{
    (void) make_block_cond_exp;
    (void) make_trace_replace;
    (void) make_noiseless_subsystem;
    (void) make_dephasing;
    (void) make_compress_idemp;
}

/* ||A||_op as a double, measured on the BALL MIDPOINTS (mirrors test_corner's
 * mid_opnorm_d). For near-zero difference / orthogonal-column matrices carrying
 * wide accumulated radii, the certified aic_mat_opnorm forms the Gram A^dag A
 * whose off-diagonal relative-Hermiticity check false-fails (bead aic-2yo);
 * collapsing to midpoints first gives the number we want (the production code
 * uses the certified mid+radius upper bound aic_corner_gamma_opnorm_ub). */
static double opnorm_d(const acb_mat_t A, slong prec)
{
    slong r = acb_mat_nrows(A), c = acb_mat_ncols(A);
    acb_mat_t M;
    acb_mat_init(M, r, c);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(A, i, j));
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    double v = dd(e);
    arb_clear(e);
    acb_mat_clear(M);
    return v;
}

/* A pseudo-random complex matrix from a fixed integer seed (deterministic, no RNG;
 * project rule). Entries in roughly [-1,1] + i[-1,1]. */
static void fill_random(acb_mat_t M, unsigned seed)
{
    slong r = acb_mat_nrows(M), c = acb_mat_ncols(M);
    unsigned s = seed * 2654435761u + 1u;
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++) {
            s = s * 1103515245u + 12345u;
            double re = ((double) ((s >> 9) & 0xffff) / 32768.0) - 1.0;
            s = s * 1103515245u + 12345u;
            double im = ((double) ((s >> 9) & 0xffff) / 32768.0) - 1.0;
            acb_set_d_d(acb_mat_entry(M, i, j), re, im);
        }
}

/* ====================================================================== T1 */
/* Diagonal identities for an arbitrary B = (+)_l M_{dims[l]} (.tex:1239-1241):
 *   pi(D) = sum_t w_t U_t^dag U_t = I_B,
 *   centrality X D = D X (the tensor identity
 *     sum_t w_t (X U_t^dag) (x) U_t == sum_t w_t U_t^dag (x) (U_t X), random X in B),
 *   and reports ||D||_proj = sum_t w_t ||U_t^dag||_op ||U_t||_op. proj_expect is the
 *   expected ||D||_proj (1 for a single block; m = #blocks for the embedded-sum
 *   diagonal, since each block contributes sum_{jk} d_l^-2 = 1 — see the
 *   .tex:1254 DISCREPANCY note in aic_dhom_diag.c). Returns the measured
 *   centrality defect (so the caller can report it / mutation-prove against it). */
static double check_diag(const slong *dims, slong nblk, double proj_expect,
                         unsigned xseed)
{
    aic_dhom_B B;
    aic_dhom_B_init(&B, dims, nblk);
    aic_dhom_diag Dg;
    aic_dhom_diag_build(&Dg, &B, PREC);
    slong nB = B.n_B, nn = nB * nB;

    /* nterms = sum_l d_l^2 (the cross-term-free embedded sum, NOT prod_l d_l^2). */
    slong expect_nterms = 0;
    for (slong l = 0; l < nblk; l++) expect_nterms += dims[l] * dims[l];
    AIC_CHECK_MSG(Dg.nterms == expect_nterms, "T1: nterms=%ld != sum d_l^2=%ld",
                  (long) Dg.nterms, (long) expect_nterms);

    /* pi(D) = sum_t w_t U_t^dag U_t == I_B ; ||D||_proj = sum_t w_t (each U_t a
     * partial isometry that is unitary on its block -> ||U_t^dag|| ||U_t|| = 1). */
    acb_mat_t piD, Ud, prod, IB;
    acb_mat_init(piD, nB, nB);
    acb_mat_init(Ud, nB, nB);
    acb_mat_init(prod, nB, nB);
    acb_mat_init(IB, nB, nB);
    acb_mat_zero(piD);
    arb_t proj;
    arb_init(proj);
    arb_zero(proj);
    for (slong t = 0; t < Dg.nterms; t++) {
        acb_mat_conjugate_transpose(Ud, Dg.U[t]);
        acb_mat_mul(prod, Ud, Dg.U[t], PREC); /* U^dag U (block-l identity) */
        acb_mat_scalar_mul_arb(prod, prod, Dg.w + t, PREC);
        acb_mat_add(piD, piD, prod, PREC);
        arb_add(proj, proj, Dg.w + t, PREC); /* w_t * 1 (partial isometry norm) */
    }
    aic_dhom_B_unit(IB, &B);
    acb_mat_sub(piD, piD, IB, PREC);
    AIC_CHECK_MSG(opnorm_d(piD, PREC) < 1e-16, "T1: pi(D) != I_B, defect %.3e",
                  opnorm_d(piD, PREC));
    AIC_CHECK_MSG(fabs(dd(proj) - proj_expect) < 1e-15,
                  "T1: ||D||_proj %.6g != expected %.6g", dd(proj), proj_expect);

    /* centrality X D = D X (.tex:1241): build both d_B^2-sized tensor sides and
     * compare. X is a random element of B (block-diagonal n_B x n_B). */
    acb_mat_t X, lhs, rhs, XUd, UtX, kl, kr;
    acb_mat_init(X, nB, nB);
    acb_mat_init(lhs, nn, nn);
    acb_mat_init(rhs, nn, nn);
    acb_mat_init(XUd, nB, nB);
    acb_mat_init(UtX, nB, nB);
    acb_mat_init(kl, nn, nn);
    acb_mat_init(kr, nn, nn);
    acb_mat_zero(X);
    for (slong l = 0; l < nblk; l++) { /* random block-diagonal X in B */
        slong d = dims[l], r0 = 0;
        for (slong ll = 0; ll < l; ll++) r0 += dims[ll];
        acb_mat_t blk;
        acb_mat_init(blk, d, d);
        fill_random(blk, xseed + (unsigned) l);
        for (slong a = 0; a < d; a++)
            for (slong b = 0; b < d; b++)
                acb_set(acb_mat_entry(X, r0 + a, r0 + b), acb_mat_entry(blk, a, b));
        acb_mat_clear(blk);
    }
    acb_mat_zero(lhs);
    acb_mat_zero(rhs);
    for (slong t = 0; t < Dg.nterms; t++) {
        acb_mat_conjugate_transpose(Ud, Dg.U[t]); /* U^dag */
        acb_mat_mul(XUd, X, Ud, PREC);            /* X U^dag */
        acb_mat_mul(UtX, Dg.U[t], X, PREC);       /* U X */
        aic_mat_kronecker(kl, XUd, Dg.U[t], PREC);
        aic_mat_kronecker(kr, Ud, UtX, PREC);
        acb_mat_scalar_mul_arb(kl, kl, Dg.w + t, PREC);
        acb_mat_scalar_mul_arb(kr, kr, Dg.w + t, PREC);
        acb_mat_add(lhs, lhs, kl, PREC);
        acb_mat_add(rhs, rhs, kr, PREC);
    }
    acb_mat_sub(lhs, lhs, rhs, PREC);
    double cdef = opnorm_d(lhs, PREC);

    acb_mat_clear(kr);
    acb_mat_clear(kl);
    acb_mat_clear(UtX);
    acb_mat_clear(XUd);
    acb_mat_clear(rhs);
    acb_mat_clear(lhs);
    acb_mat_clear(X);
    arb_clear(proj);
    acb_mat_clear(IB);
    acb_mat_clear(prod);
    acb_mat_clear(Ud);
    acb_mat_clear(piD);
    aic_dhom_diag_clear(&Dg);
    aic_dhom_B_clear(&B);
    return cdef;
}

/* Pauli unit tests: per-S unitarity; then the diagonal identities pi(D)=I_B,
 * centrality X D = D X, and ||D||_proj for SINGLE blocks (M_2,M_3,M_4) AND
 * DIRECT SUMS (M_2(+)M_2, M_2(+)M_3). The direct-sum cases are the BLOCKER guard:
 * the old joint-unitary/Cartesian-product build (.tex:1254 verbatim) was
 * NON-central there (||XD-DX|| ~ 0.5-1.0); the embedded-sum diagonal is central
 * to ~1e-20. MUTATION (proven 2026-05-30): restoring the joint build makes the
 * direct-sum centrality go RED at ~1.0 while the single-block cases stay green. */
static void test_pauli(void)
{
    printf("-- T1 Pauli diagonal --\n");
    /* each S_{jk} unitary: S^dag S = I */
    for (slong d = 2; d <= 4; d++)
        for (slong j = 0; j < d; j++)
            for (slong k = 0; k < d; k++) {
                acb_mat_t S, Sd, G, one;
                acb_mat_init(S, d, d);
                acb_mat_init(Sd, d, d);
                acb_mat_init(G, d, d);
                acb_mat_init(one, d, d);
                aic_dhom_pauli(S, d, j, k, PREC);
                acb_mat_conjugate_transpose(Sd, S);
                acb_mat_mul(G, Sd, S, PREC);
                acb_mat_one(one);
                acb_mat_sub(G, G, one, PREC);
                AIC_CHECK_MSG(opnorm_d(G, PREC) < 1e-18,
                              "T1: S_{%ld,%ld} (d=%ld) not unitary", (long) j,
                              (long) k, (long) d);
                acb_mat_clear(one);
                acb_mat_clear(G);
                acb_mat_clear(Sd);
                acb_mat_clear(S);
            }

    /* single-block diagonals: ||D||_proj = 1, centrality ~ machine-zero. */
    for (slong d = 2; d <= 4; d++) {
        slong dims[1] = {d};
        double cdef = check_diag(dims, 1, /*proj_expect*/ 1.0, (unsigned) (17 + d));
        AIC_CHECK_MSG(cdef < 1e-16, "T1: single-block centrality (d=%ld) %.3e",
                      (long) d, cdef);
    }

    /* DIRECT-SUM diagonals (the BLOCKER guard). proj_expect = m (#blocks): each
     * block contributes sum_{jk} d_l^-2 = 1 to ||D||_proj (the embedded-sum repr;
     * see the .tex:1254 DISCREPANCY note in aic_dhom_diag.c). Centrality must be
     * ~machine-zero (the old joint build had ~0.5-1.0 here). */
    {
        slong dims2[2] = {2, 2};
        double c22 = check_diag(dims2, 2, /*proj_expect*/ 2.0, 31u);
        printf("   M_2(+)M_2: ||XD-DX|| = %.3e (||D||_proj = 2)\n", c22);
        AIC_CHECK_MSG(c22 < 1e-16, "T1: M_2(+)M_2 centrality %.3e (non-central "
                      "diagonal — the BLOCKER)", c22);
        slong dims23[2] = {2, 3};
        double c23 = check_diag(dims23, 2, /*proj_expect*/ 2.0, 41u);
        printf("   M_2(+)M_3: ||XD-DX|| = %.3e (||D||_proj = 2)\n", c23);
        AIC_CHECK_MSG(c23 < 1e-16, "T1: M_2(+)M_3 centrality %.3e (non-central "
                      "diagonal — the BLOCKER)", c23);
    }
}

/* Set v(E_i) = U E_i U^dag for B = M_n (single block, n_B = n), an EXACT
 * *-isomorphism M_n -> A = M_n. U is n x n unitary. */
static void set_v_conjugation(aic_dhom_v *v, const acb_mat_t U, slong prec)
{
    slong n = v->n;
    acb_mat_t Ud, Ei, tmp;
    acb_mat_init(Ud, n, n);
    acb_mat_init(Ei, n, n);
    acb_mat_init(tmp, n, n);
    acb_mat_conjugate_transpose(Ud, U);
    for (slong i = 0; i < v->B->dim_B; i++) {
        aic_dhom_B_matunit(Ei, v->B, i);     /* matrix unit E_i (n x n)    */
        acb_mat_mul(tmp, U, Ei, prec);       /* U E_i                      */
        acb_mat_mul(v->vE[i], tmp, Ud, prec);/* U E_i U^dag                */
    }
    acb_mat_clear(tmp);
    acb_mat_clear(Ei);
    acb_mat_clear(Ud);
}

/* ====================================================================== T2 */
/* eta=0 oracle: A = M_n (identity channel), v(X) = U X U^dag exact *-iso ->
 * G_v = 0 (machine-zero) on all basis pairs; lem_approx returns v unchanged. */
static void test_eta0_oracle(void)
{
    printf("-- T2 eta=0 oracle --\n");
    for (slong n = 2; n <= 3; n++) {
        aic_ucp_kraus phi;
        make_identity(&phi, n);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
        AIC_CHECK_MSG(ae.A.dim_A == n * n, "T2(n=%ld): dim_A=%ld != n^2",
                      (long) n, (long) ae.A.dim_A);

        slong dims[1] = {n};
        aic_dhom_B B;
        aic_dhom_B_init(&B, dims, 1);
        aic_dhom_v v;
        aic_dhom_v_init(&v, &B, &ae.A);

        acb_mat_t U;
        acb_mat_init(U, n, n);
        make_fixed_unitary(U, n, PREC);
        set_v_conjugation(&v, U, PREC);

        /* G_v = 0 on every basis pair (the headline machine-zero oracle). */
        arb_t sweep;
        arb_init(sweep);
        aic_dhom_defect_sweep(sweep, &v, PREC);
        AIC_CHECK_MSG(dd(sweep) < 1e-20, "T2(n=%ld): G_v sweep %.3e != 0",
                      (long) n, dd(sweep));

        /* snapshot vE for the "unchanged" check after lem_approx */
        acb_mat_t *vE0 = flint_malloc((size_t) B.dim_B * sizeof(acb_mat_t));
        for (slong i = 0; i < B.dim_B; i++) {
            acb_mat_init(vE0[i], n, n);
            acb_mat_set(vE0[i], v.vE[i]);
        }

        aic_dhom_approx_stats st;
        aic_dhom_approx(&v, /*eps_target*/ 0.0, /*tol_abs*/ 1e-18,
                        /*unit_tol*/ 1e-12, /*max_steps*/ 8, PREC, &st);
        AIC_CHECK_MSG(st.steps == 0, "T2(n=%ld): lem_approx took %ld steps on an "
                      "exact hom (expected 0)", (long) n, (long) st.steps);

        /* v unchanged: ||vtilde - v||_op = 0 for every basis op. */
        double maxchg = 0.0;
        acb_mat_t diff;
        acb_mat_init(diff, n, n);
        for (slong i = 0; i < B.dim_B; i++) {
            acb_mat_sub(diff, v.vE[i], vE0[i], PREC);
            double c = opnorm_d(diff, PREC);
            if (c > maxchg) maxchg = c;
        }
        AIC_CHECK_MSG(maxchg < 1e-18, "T2(n=%ld): lem_approx changed v by %.3e",
                      (long) n, maxchg);
        acb_mat_clear(diff);
        for (slong i = 0; i < B.dim_B; i++) acb_mat_clear(vE0[i]);
        flint_free(vE0);

        arb_clear(sweep);
        acb_mat_clear(U);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ====================================================================== T3 */
/* delta>0, eps=0 (Johnson): perturb an exact hom v(X)=U X U^dag by delta*R (||R||
 * <= 1 per basis op). lem_approx (eps=0) drives the defect to machine-zero (exact
 * C* algebra A=M_n; Johnson 1988 Thm 3.1, .tex:1313) and ||vtilde - v_pert|| <=
 * O(delta). */
static void test_delta_eps0(void)
{
    printf("-- T3 delta>0, eps=0 (Johnson) --\n");
    slong n = 2;
    double delta = 0.02;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);

    slong dims[1] = {n};
    aic_dhom_B B;
    aic_dhom_B_init(&B, dims, 1);
    aic_dhom_v v;
    aic_dhom_v_init(&v, &B, &ae.A);

    acb_mat_t U;
    acb_mat_init(U, n, n);
    make_fixed_unitary(U, n, PREC);
    set_v_conjugation(&v, U, PREC);

    /* perturb: vE[i] += delta * R_i, *-SYMMETRICALLY so v stays *-linear (a
     * delta-homomorphism is *-linear by definition, .tex:1311 premise): for the
     * matrix-unit index i = mu(a,b) we set R at the *-conjugate index j = mu(b,a)
     * to R_i^dag, leaving v(E_ba)^dag = v(E_ab) intact. The off-diagonal pair
     * (a!=b) is perturbed by an arbitrary R_i (with the conjugate forced); the
     * diagonal (a==b) is perturbed by a HERMITIAN R_i (R := (R+R^dag)/2). */
    arb_t dlt, inv, ahalf;
    arb_init(dlt);
    arb_init(inv);
    arb_init(ahalf);
    arb_set_d(dlt, delta);
    arb_set_d(ahalf, 0.5);
    acb_mat_t R, Rn, Rd;
    acb_mat_init(R, n, n);
    acb_mat_init(Rn, n, n);
    acb_mat_init(Rd, n, n);
    /* B = M_n single block: linear index i = a*n + b. */
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++) {
            slong i = a * n + b, j = b * n + a;
            if (j < i) continue; /* handle each *-pair once */
            fill_random(R, (unsigned) (101 + i));
            if (a == b) { /* Hermitian perturbation keeps the diagonal *-symmetric */
                acb_mat_conjugate_transpose(Rd, R);
                acb_mat_add(R, R, Rd, PREC);
                acb_mat_scalar_mul_arb(R, R, ahalf, PREC);
            }
            double rn = opnorm_d(R, PREC);
            arb_set_d(inv, 1.0 / rn);
            acb_mat_scalar_mul_arb(Rn, R, inv, PREC);  /* ||Rn||_op = 1      */
            acb_mat_scalar_mul_arb(Rn, Rn, dlt, PREC); /* delta * Rn         */
            acb_mat_add(v.vE[i], v.vE[i], Rn, PREC);
            if (j != i) { /* force the *-conjugate: vE[j] += (delta Rn)^dag */
                acb_mat_conjugate_transpose(Rd, Rn);
                acb_mat_add(v.vE[j], v.vE[j], Rd, PREC);
            }
        }
    acb_mat_clear(Rd);

    /* snapshot v_pert */
    acb_mat_t *vp = flint_malloc((size_t) B.dim_B * sizeof(acb_mat_t));
    for (slong i = 0; i < B.dim_B; i++) {
        acb_mat_init(vp[i], n, n);
        acb_mat_set(vp[i], v.vE[i]);
    }

    arb_t sweep;
    arb_init(sweep);
    aic_dhom_defect_sweep(sweep, &v, PREC);
    double d0 = dd(sweep);
    AIC_CHECK_MSG(d0 > 1e-3, "T3: perturbed defect %.3e too small to be a test",
                  d0);
    printf("   initial defect delta_0 = %.3e\n", d0);

    aic_dhom_approx_stats st;
    aic_dhom_approx(&v, /*eps_target*/ 0.0, /*tol_abs*/ 1e-18,
                    /*unit_tol*/ 1e-6, /*max_steps*/ 40, PREC, &st);
    printf("   steps=%ld  delta_final=%.3e\n", (long) st.steps, st.delta_final);
    AIC_CHECK_MSG(st.delta_final < 1e-15, "T3: eps=0 defect %.3e did not reach "
                  "machine-zero (Johnson)", st.delta_final);
    AIC_CHECK_MSG(st.steps >= 1, "T3: lem_approx took 0 steps on a perturbed hom");

    /* ||vtilde - v_pert|| <= C*delta (report C) */
    double maxchg = 0.0;
    acb_mat_t diff;
    acb_mat_init(diff, n, n);
    for (slong i = 0; i < B.dim_B; i++) {
        acb_mat_sub(diff, v.vE[i], vp[i], PREC);
        double c = opnorm_d(diff, PREC);
        if (c > maxchg) maxchg = c;
    }
    double C = maxchg / delta;
    printf("   ||vtilde - v_pert||_max = %.3e = %.2f * delta\n", maxchg, C);
    AIC_CHECK_MSG(maxchg <= 10.0 * delta, "T3: ||vtilde - v_pert|| = %.3e not "
                  "O(delta=%.3e)", maxchg, delta);
    acb_mat_clear(diff);

    for (slong i = 0; i < B.dim_B; i++) acb_mat_clear(vp[i]);
    flint_free(vp);
    arb_clear(sweep);
    acb_mat_clear(Rn);
    acb_mat_clear(R);
    arb_clear(ahalf);
    arb_clear(inv);
    arb_clear(dlt);
    acb_mat_clear(U);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* Project R (n x n) onto A = span{B_k}: Rp = sum_k <B_k,R>_F B_k. Keeps a
 * perturbation IN A (so the perturbed v stays a map into A). */
static void project_onto_A(acb_mat_t Rp, const aic_ecstar *A, const acb_mat_t R,
                           slong prec)
{
    slong n = A->n;
    acb_mat_zero(Rp);
    acb_t c, z;
    acb_init(c);
    acb_init(z);
    acb_mat_t tmp;
    acb_mat_init(tmp, n, n);
    for (slong k = 0; k < A->dim_A; k++) {
        acb_zero(c);
        for (slong x = 0; x < n; x++)
            for (slong y = 0; y < n; y++) {
                acb_conj(z, acb_mat_entry(A->B[k], x, y));
                acb_mul(z, z, acb_mat_entry(R, x, y), prec);
                acb_add(c, c, z, prec);
            }
        acb_mat_scalar_mul_acb(tmp, A->B[k], c, prec);
        acb_mat_add(Rp, Rp, tmp, prec);
    }
    acb_mat_clear(tmp);
    acb_clear(z);
    acb_clear(c);
}

/* Multi-block conditional expectation Phi(X) = sum_l Pi_l X Pi_l on C^n, where
 * Pi_l projects onto the contiguous coordinate block l (dims[l] coords, offset
 * sum_{l'<l} dims[l']). Img Phi = (+)_l M_{dims[l]} (block-diagonal matrices), an
 * EXACTLY idempotent (eta=0) channel. Kraus: one Pi_l per block. `phi`
 * aic_ucp_kraus_init'd here. Generalizes make_block_cond_exp (which is the m=2,
 * dims={m,d-m} case) to any number of blocks (needed for the T5 m=1,2,3 sweep). */
static void make_multiblock_cond_exp(aic_ucp_kraus *phi, const slong *dims,
                                     slong nblk)
{
    slong n = 0;
    for (slong l = 0; l < nblk; l++) n += dims[l];
    aic_ucp_kraus_init(phi, n, n, nblk);
    slong off = 0;
    for (slong l = 0; l < nblk; l++) {
        for (slong i = 0; i < dims[l]; i++)
            acb_set_si(acb_mat_entry(phi->K[l], off + i, off + i), 1);
        off += dims[l];
    }
}

/* Convex mix of make_multiblock_cond_exp(dims) [Img Phi = (+)_l M_{dims[l]}] with
 * its unitary conjugate (mirrors make_mixconj for the single-block compress
 * channel): a genuine eta>0 almost-idempotent channel whose A ~ (+)_l M_{dims[l]}.
 * `out` aic_ucp_kraus_init'd here. Used by the T5 direct-sum / number-of-blocks
 * canary points. */
static void make_mixconj_bce(aic_ucp_kraus *out, const slong *dims, slong nblk,
                             double t, slong prec)
{
    aic_ucp_kraus base, conj;
    make_multiblock_cond_exp(&base, dims, nblk);
    make_conjugated(&conj, &base, prec);
    slong n = base.dim_H;
    aic_ucp_kraus_init(out, n, n, base.r + conj.r);
    arb_t s;
    arb_init(s);
    arb_set_d(s, sqrt(1.0 - t));
    for (slong a = 0; a < base.r; a++)
        acb_mat_scalar_mul_arb(out->K[a], base.K[a], s, prec);
    arb_set_d(s, sqrt(t));
    for (slong b = 0; b < conj.r; b++)
        acb_mat_scalar_mul_arb(out->K[base.r + b], conj.K[b], s, prec);
    arb_clear(s);
    aic_ucp_kraus_clear(&conj);
    aic_ucp_kraus_clear(&base);
}

/* General lem_approx reduction run. `ae` is a built eta>0 algebra A (~ B =
 * (+)_l M_{dims[l]}); `roff[l]` is the row-offset of block l in M_d (the ambient
 * dim n = ae->A.n). Sets v(E^{(l)}_ab) = Phi_tilde(iota_l(E_ab)) (the natural
 * delta-hom B -> A, iota_l = embedding of block l at roff[l]), *-symmetrically
 * perturbs v by delta (projected into A), runs lem_approx, and returns
 * C = defect_final / eps_assoc. Reports d0, eps, steps via out params. */
static double run_reduce(aic_assoc_ecstar *ae, const slong *dims, slong nblk,
                         const slong *roff, double delta, double *d0_out,
                         double *eps_out, slong *steps_out)
{
    slong prec = PREC, n = ae->A.n;
    arb_t eps;
    arb_init(eps);
    aic_ecstar_defect_assoc(eps, &ae->A, prec); /* A's intrinsic epsilon */
    double eps_v = dd(eps);

    aic_dhom_B B;
    aic_dhom_B_init(&B, dims, nblk);
    AIC_CHECK_MSG(ae->A.dim_A == B.dim_B, "run_reduce: dim_A=%ld != dim_B=%ld",
                  (long) ae->A.dim_A, (long) B.dim_B);
    aic_dhom_v v;
    aic_dhom_v_init(&v, &B, &ae->A);

    acb_mat_t iE, vE;
    acb_mat_init(iE, n, n);
    acb_mat_init(vE, n, n);
    for (slong l = 0; l < nblk; l++) {
        slong dl = dims[l], r0 = roff[l];
        for (slong a = 0; a < dl; a++)
            for (slong b = 0; b < dl; b++) {
                slong i = aic_dhom_B_index(&B, l, a, b);
                acb_mat_zero(iE);
                acb_set_si(acb_mat_entry(iE, r0 + a, r0 + b), 1); /* iota_l(E_ab) */
                aic_assoc_superop_apply(vE, ae->S_tilde, iE, prec); /* Phi_tilde */
                acb_mat_set(v.vE[i], vE);
            }
    }

    /* *-symmetric delta-perturbation, projected into A. Pair up matrix-unit index
     * i = mu(l,a,b) with its *-conjugate j = mu(l,b,a) within each block. */
    arb_t dl_arb, inv, ah;
    arb_init(dl_arb);
    arb_init(inv);
    arb_init(ah);
    arb_set_d(dl_arb, delta);
    arb_set_d(ah, 0.5);
    acb_mat_t R, Rd, Rp;
    acb_mat_init(R, n, n);
    acb_mat_init(Rd, n, n);
    acb_mat_init(Rp, n, n);
    for (slong l = 0; l < nblk; l++)
        for (slong a = 0; a < dims[l]; a++)
            for (slong b = a; b < dims[l]; b++) {
                slong i = aic_dhom_B_index(&B, l, a, b);
                slong j = aic_dhom_B_index(&B, l, b, a);
                fill_random(R, (unsigned) (201 + i));
                if (a == b) {
                    acb_mat_conjugate_transpose(Rd, R);
                    acb_mat_add(R, R, Rd, prec);
                    acb_mat_scalar_mul_arb(R, R, ah, prec);
                }
                project_onto_A(Rp, &ae->A, R, prec);
                double rn = opnorm_d(Rp, prec);
                arb_set_d(inv, 1.0 / rn);
                acb_mat_scalar_mul_arb(Rp, Rp, inv, prec);
                acb_mat_scalar_mul_arb(Rp, Rp, dl_arb, prec);
                acb_mat_add(v.vE[i], v.vE[i], Rp, prec);
                if (j != i) {
                    acb_mat_conjugate_transpose(Rd, Rp);
                    acb_mat_add(v.vE[j], v.vE[j], Rd, prec);
                }
            }

    arb_t sw0;
    arb_init(sw0);
    aic_dhom_defect_sweep(sw0, &v, prec);
    double d0 = dd(sw0);

    aic_dhom_approx_stats st;
    aic_dhom_approx(&v, /*eps_target*/ 4.0 * eps_v, /*tol_abs*/ 1e-13,
                    /*unit_tol*/ 1.0, /*max_steps*/ 30, prec, &st);
    double C = st.delta_final / eps_v;

    if (d0_out) *d0_out = d0;
    if (eps_out) *eps_out = eps_v;
    if (steps_out) *steps_out = st.steps;

    arb_clear(sw0);
    acb_mat_clear(Rp);
    acb_mat_clear(Rd);
    acb_mat_clear(R);
    arb_clear(ah);
    arb_clear(inv);
    arb_clear(dl_arb);
    acb_mat_clear(vE);
    acb_mat_clear(iE);
    arb_clear(eps);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    return C;
}

/* Single-block M_m convenience wrapper over run_reduce: A from make_mixconj(d,m,t),
 * B = M_m embedded top-left (roff = 0). */
static double t4_run(slong d, slong m, double t, double delta, double *d0_out,
                     double *eps_out, slong *steps_out)
{
    aic_ucp_kraus phi;
    make_mixconj(&phi, d, m, t, PREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    slong dims[1] = {m}, roff[1] = {0};
    double C = run_reduce(&ae, dims, 1, roff, delta, d0_out, eps_out, steps_out);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return C;
}

/* Multi-block convenience wrapper: A from make_mixconj_bce(dims) (~ (+)_l M_{dims[l]},
 * blocks laid out contiguously so roff[l] = sum_{l'<l} dims[l']), reduce, return C. */
static double mb_run(const slong *dims, slong nblk, double t, double delta,
                     double *d0_out, double *eps_out, slong *steps_out)
{
    aic_ucp_kraus phi;
    make_mixconj_bce(&phi, dims, nblk, t, PREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    slong roff[8];
    AIC_CHECK_MSG(nblk <= 8, "mb_run: too many blocks %ld", (long) nblk);
    slong o = 0;
    for (slong l = 0; l < nblk; l++) { roff[l] = o; o += dims[l]; }
    double C = run_reduce(&ae, dims, nblk, roff, delta, d0_out, eps_out, steps_out);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return C;
}

/* ====================================================================== T4 */
/* eps>0 GENUINE almost-idempotent A (make_mixconj): a delta-hom v into it,
 * perturbed by delta >> eps; lem_approx reduces the defect to O(eps). Report
 * C = defect/eps. STAR NON-VACUITY: the source mutation replacing G_v's star with
 * plain matmul makes the defect floor at O(1) (not O(eps)) at eta>0 -> the C bound
 * RED (mutation-proven separately; see the module report). */
static void test_eps_pos(void)
{
    printf("-- T4 eps>0 genuine almost-idempotent A --\n");
    double d0, eps_v;
    slong steps;
    double C = t4_run(/*d*/ 5, /*m*/ 3, /*t*/ 0.015, /*delta*/ 0.05, &d0, &eps_v,
                      &steps);
    printf("   eps_assoc=%.3e  delta_0=%.3e  steps=%ld  C=defect/eps=%.3f\n",
           eps_v, d0, (long) steps, C);
    AIC_CHECK_MSG(d0 > 20.0 * eps_v, "T4: delta_0=%.3e not >> eps=%.3e (not "
                  "exercising the reduction)", d0, eps_v);
    AIC_CHECK_MSG(steps >= 1, "T4: lem_approx took 0 steps (defect already small)");
    AIC_CHECK_MSG(C < 6.0, "T4: defect/eps = %.3f not O(1) (lem_approx did not "
                  "reduce to O(eps))", C);
}

/* ====================================================================== T5 */
/* Universality canary (.tex:484/460): the th_main constant must be UNIVERSAL —
 * independent of BOTH the block dimension AND the number of blocks. C = defect/eps
 * is the proxy. Two sweeps:
 *
 *  (A) BLOCK-DIMENSION sweep — single block B = M_2, M_3, M_4 (dim_B = 4,9,16).
 *      C must not grow with the block dim. The EXPLICIT Pauli diagonal (per block
 *      ||D_l||_proj = 1 exactly) keeps it flat; a naive non-Pauli diagonal would
 *      make C grow ~ dim. MUTATION: replace the Pauli diagonal with a non-diagonal
 *      -> C grows -> RED (mutation-proven separately, see the module report).
 *
 *  (B) NUMBER-OF-BLOCKS sweep — B = M_2, M_2(+)M_2, M_2(+)M_2(+)M_2 (m = 1,2,3) at
 *      fixed delta. This is the FIX-3 escalation question: the embedded-sum diagonal
 *      has sum_t ||A_t|| ||B_t|| = m (NOT 1; see the .tex:1254 DISCREPANCY note in
 *      aic_dhom_diag.c), so the w' bound is O(m*delta) and the reduced C COULD carry
 *      an m-factor. We REPORT C(m) and step counts and assert the MEASURED behavior:
 *      C does NOT grow with m (it is in fact LARGEST at m=1; the m-factor in the
 *      worst-case w' bound does not manifest as growth in the achieved C, because
 *      the per-step quadratic contraction reaches the O(eps) floor regardless).
 *      The honest threshold: max over m of C(m) stays below a fixed block-count-
 *      INDEPENDENT ceiling, AND C(3) does not exceed C(1) (no upward block-count
 *      trend). If a future fixture DID show C growing with m, this assert turns RED
 *      and the .tex:484 universal-constant concern (aic-1bc) escalates. */
static void test_canary(void)
{
    printf("-- T5 universality canary --\n");
    double eps_v, d0;
    slong steps;

    /* ---- (A) block-dimension sweep: single block M_2, M_3, M_4 ---- */
    double Cdim[3];
    Cdim[0] = t4_run(4, 2, 0.015, 0.05, &d0, &eps_v, &steps);
    printf("   (A) B=M_2  dim_B=4 : eps=%.3e C=%.3f (steps %ld)\n", eps_v, Cdim[0],
           (long) steps);
    Cdim[1] = t4_run(5, 3, 0.015, 0.05, &d0, &eps_v, &steps);
    printf("   (A) B=M_3  dim_B=9 : eps=%.3e C=%.3f (steps %ld)\n", eps_v, Cdim[1],
           (long) steps);
    Cdim[2] = t4_run(6, 4, 0.015, 0.05, &d0, &eps_v, &steps);
    printf("   (A) B=M_4  dim_B=16: eps=%.3e C=%.3f (steps %ld)\n", eps_v, Cdim[2],
           (long) steps);
    double dmin = Cdim[0], dmax = Cdim[0];
    for (int i = 1; i < 3; i++) {
        if (Cdim[i] < dmin) dmin = Cdim[i];
        if (Cdim[i] > dmax) dmax = Cdim[i];
    }
    printf("   (A) C in [%.3f, %.3f], ratio Cmax/Cmin = %.3f (dim_B 4..16)\n", dmin,
           dmax, dmax / dmin);
    AIC_CHECK_MSG(dmax / dmin < 2.0, "T5(A): C ratio %.3f grows with dim_B (the "
                  "Pauli diagonal should keep it dimension-independent)",
                  dmax / dmin);

    /* ---- (B) number-of-blocks sweep: m = 1,2,3 of M_2 ---- */
    double Cm[3];
    slong sm[3];
    Cm[0] = t4_run(4, 2, 0.015, 0.05, &d0, &eps_v, &sm[0]); /* m=1: M_2 */
    printf("   (B) m=1 B=M_2          dim_B=4 : eps=%.3e C=%.3f (steps %ld)\n",
           eps_v, Cm[0], (long) sm[0]);
    {
        slong dims2[2] = {2, 2};
        Cm[1] = mb_run(dims2, 2, 0.015, 0.05, &d0, &eps_v, &sm[1]); /* m=2 */
        printf("   (B) m=2 B=M_2(+)M_2    dim_B=8 : eps=%.3e C=%.3f (steps %ld)\n",
               eps_v, Cm[1], (long) sm[1]);
        slong dims3[3] = {2, 2, 2};
        Cm[2] = mb_run(dims3, 3, 0.015, 0.05, &d0, &eps_v, &sm[2]); /* m=3 */
        printf("   (B) m=3 B=M_2(+)M_2(+)M_2 dim_B=12: eps=%.3e C=%.3f (steps %ld)\n",
               eps_v, Cm[2], (long) sm[2]);
    }
    double cmmax = Cm[0];
    for (int i = 1; i < 3; i++) if (Cm[i] > cmmax) cmmax = Cm[i];
    printf("   (B) max_m C = %.3f; C(1)=%.3f C(2)=%.3f C(3)=%.3f -> C does %s "
           "with m\n", cmmax, Cm[0], Cm[1], Cm[2],
           Cm[2] > Cm[0] ? "GROW" : "NOT grow");
    /* VERDICT: C does not grow with the number of blocks m. Two guards:
     *  - the achieved C stays below a fixed (block-count-INDEPENDENT) ceiling; and
     *  - no upward trend: C(3) <= C(1) * 1.2 (a small slack for the eps-regime
     *    spread across the three fixtures). The worst-case O(m*delta) w' bound does
     *    NOT show up as C growth (the quadratic Newton contraction reaches the
     *    O(eps) floor at every m). If this RED-fires, the .tex:484 universal-constant
     *    concern (aic-1bc) escalates — do NOT relax the threshold to mask it. */
    AIC_CHECK_MSG(cmmax < 6.0, "T5(B): max_m C = %.3f exceeds block-count-"
                  "independent ceiling (the O(m*delta) w' bound is manifesting)",
                  cmmax);
    AIC_CHECK_MSG(Cm[2] <= 1.2 * Cm[0], "T5(B): C(m=3)=%.3f > C(m=1)=%.3f — C is "
                  "GROWING with the number of blocks m (the .tex:484 universal-"
                  "constant concern, escalate aic-1bc)", Cm[2], Cm[0]);
}

/* ====================================================================== T6 */
/* prop_delta_hominc bounds (.tex:1194). Exact hom v(X)=U X U^dag: ||v||=1, a=1,
 * unit defect = 0. Perturbed (delta>0): ||v|| <= 1+O(delta), a >= 1-O(delta) when
 * a > 2 delta, unit defect <= O(delta). */
static void test_prop_bounds(void)
{
    printf("-- T6 prop_delta_hominc bounds --\n");
    slong n = 3;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    slong dims[1] = {n};
    aic_dhom_B B;
    aic_dhom_B_init(&B, dims, 1);
    aic_dhom_v v;
    aic_dhom_v_init(&v, &B, &ae.A);
    acb_mat_t U;
    acb_mat_init(U, n, n);
    make_fixed_unitary(U, n, PREC);
    set_v_conjugation(&v, U, PREC);

    /* exact hom: matrix units E_ab map to U E_ab U^dag, each ||.||_op = 1 (U
     * unitary). a = ||v|| = 1, unit defect = 0. */
    arb_t nub, nlb, ud;
    arb_init(nub);
    arb_init(nlb);
    arb_init(ud);
    aic_dhom_prop_bounds(nub, nlb, ud, &v, PREC);
    printf("   exact hom: ||v||=%.6f  a=%.6f  unit_def=%.3e\n", dd(nub), dd(nlb),
           dd(ud));
    AIC_CHECK_MSG(fabs(dd(nub) - 1.0) < 1e-12, "T6: exact-hom ||v|| %.6f != 1",
                  dd(nub));
    AIC_CHECK_MSG(fabs(dd(nlb) - 1.0) < 1e-12, "T6: exact-hom a %.6f != 1",
                  dd(nlb));
    AIC_CHECK_MSG(dd(ud) < 1e-12, "T6: exact-hom unit defect %.3e != 0", dd(ud));

    /* perturb diagonally by a small Hermitian delta so v stays *-symmetric; the
     * bounds must hold: ||v|| <= 1+8 delta, a >= 1-8 delta, unit_def <= 8 delta. */
    double delta = 0.01;
    arb_t dl, ah;
    arb_init(dl);
    arb_init(ah);
    arb_set_d(dl, delta);
    arb_set_d(ah, 0.5);
    acb_mat_t R, Rd, Rn;
    acb_mat_init(R, n, n);
    acb_mat_init(Rd, n, n);
    acb_mat_init(Rn, n, n);
    for (slong a = 0; a < n; a++)
        for (slong b = a; b < n; b++) {
            slong i = a * n + b, j = b * n + a;
            fill_random(R, (unsigned) (501 + i));
            if (a == b) {
                acb_mat_conjugate_transpose(Rd, R);
                acb_mat_add(R, R, Rd, PREC);
                acb_mat_scalar_mul_arb(R, R, ah, PREC);
            }
            double rn = opnorm_d(R, PREC);
            arb_t inv;
            arb_init(inv);
            arb_set_d(inv, 1.0 / rn);
            acb_mat_scalar_mul_arb(Rn, R, inv, PREC);
            acb_mat_scalar_mul_arb(Rn, Rn, dl, PREC);
            acb_mat_add(v.vE[i], v.vE[i], Rn, PREC);
            if (j != i) {
                acb_mat_conjugate_transpose(Rd, Rn);
                acb_mat_add(v.vE[j], v.vE[j], Rd, PREC);
            }
            arb_clear(inv);
        }
    aic_dhom_prop_bounds(nub, nlb, ud, &v, PREC);
    printf("   perturbed (delta=%.2f): ||v||=%.6f  a=%.6f  unit_def=%.3e\n",
           delta, dd(nub), dd(nlb), dd(ud));
    AIC_CHECK_MSG(dd(nub) <= 1.0 + 8.0 * delta, "T6: ||v|| %.6f > 1+O(delta)",
                  dd(nub));
    AIC_CHECK_MSG(dd(nlb) >= 1.0 - 8.0 * delta, "T6: a %.6f < 1-O(delta)",
                  dd(nlb));
    AIC_CHECK_MSG(dd(ud) <= 8.0 * delta, "T6: unit defect %.3e > O(delta)",
                  dd(ud));

    acb_mat_clear(Rn);
    acb_mat_clear(Rd);
    acb_mat_clear(R);
    arb_clear(ah);
    arb_clear(dl);
    arb_clear(ud);
    arb_clear(nlb);
    arb_clear(nub);
    acb_mat_clear(U);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    flint_set_num_threads(1);
    test_pauli();
    test_eta0_oracle();
    test_delta_eps0();
    test_eps_pos();
    test_canary();
    test_prop_bounds();
    aic_test_report("test_dhom");
    return 0;
}
