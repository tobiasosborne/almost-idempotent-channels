/* test_idemp.h — exactly-idempotent UCP channel constructors for the
 * idemp_structure cross-checks (bead aic-wuh). Header-only (static); split out of
 * test_idemp.c to keep each test file <= 200 LOC (Rule 10). Heisenberg convention
 * Phi(X) = sum_a K_a^dag X K_a (see include/aic_ucp.h).
 *
 * Each channel has KNOWN (dim_M, dim_A); the values are documented per builder
 * and asserted by check_decomp in test_idemp.c.
 */
#ifndef TEST_IDEMP_H
#define TEST_IDEMP_H

#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_idemp.h"
#include "aic_mat.h"
#include "aic_test.h"
#include "aic_ucp.h"

#define PREC 53
static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* Block-diagonal conditional expectation on C^d, blocks [0,m) and [m,d):
 * Phi(X) = P_1 X P_1 + P_2 X P_2 (K_1=P_1, K_2=P_2, P_1+P_2=1_d). Unital +
 * trace-preserving => M=H (dim_M=d). Img Phi = B(C^m)(+)B(C^{d-m}),
 * dim_A = m^2 + (d-m)^2. */
static void make_block_cond_exp(aic_ucp_kraus *phi, slong d, slong m)
{
    aic_ucp_kraus_init(phi, d, d, 2);
    for (slong i = 0; i < m; i++) acb_set_si(acb_mat_entry(phi->K[0], i, i), 1);
    for (slong i = m; i < d; i++) acb_set_si(acb_mat_entry(phi->K[1], i, i), 1);
}

/* Trace-replace on C^k: Phi(X) = (Tr X / k) 1_k. Kraus K_{ij}=(1/sqrt k)|e_i><e_j|
 * (k^2 ops). Unital + trace-preserving => M=H (dim_M=k); Img Phi = scalars,
 * dim_A=1. */
static void make_trace_replace(aic_ucp_kraus *phi, slong k)
{
    aic_ucp_kraus_init(phi, k, k, k * k);
    double inv = 1.0 / sqrt((double) k);
    slong a = 0;
    for (slong i = 0; i < k; i++)
        for (slong j = 0; j < k; j++) {
            acb_set_d(acb_mat_entry(phi->K[a], i, j), inv);
            a++;
        }
}

/* Noiseless subsystem on H = C^dL (x) C^dE: Phi(X) = (Tr_E X) (x) (1_E / dE).
 * Heisenberg Kraus K_{jk} = (1_L (x) |e_j><e_k|)/sqrt(dE). Img Phi ~ B(C^dL),
 * dim_A = dL^2; M = H (Q=1), dim_M = dL*dE. H row I = a*dE + j (left-major). */
static void make_noiseless_subsystem(aic_ucp_kraus *phi, slong dL, slong dE)
{
    slong n = dL * dE;
    aic_ucp_kraus_init(phi, n, n, dE * dE);
    double inv = 1.0 / sqrt((double) dE);
    slong op = 0;
    for (slong j = 0; j < dE; j++)
        for (slong k = 0; k < dE; k++) {
            for (slong a = 0; a < dL; a++)
                acb_set_d(acb_mat_entry(phi->K[op], a * dE + j, a * dE + k), inv);
            op++;
        }
}

/* Identity Phi = id on C^d: K_0 = 1_d. dim_M=d, dim_A=d^2. */
static void make_identity(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, 1);
    acb_mat_one(phi->K[0]);
}

/* Dephasing Phi(X) = diag(X) on C^d: K_i = |e_i><e_i| (d ops). Img Phi = diagonal
 * matrices (commutative), dim_A=d; M=H, dim_M=d. */
static void make_dephasing(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, d);
    for (slong i = 0; i < d; i++)
        acb_set_si(acb_mat_entry(phi->K[i], i, i), 1);
}

/* Exactly-idempotent, PROPER carrier M=span(e0..e_{m-1}) (unital, NOT
 * trace-preserving). Phi(X) = Pi_M X Pi_M + X[0,0](1-Pi_M). K_0 = Pi_M,
 * K_{1+b} = |e_0><e_{m+b}|. Carrier Q = Pi_M + (d-m)|e0><e0|, support = M
 * (dim_M=m<d). Img Phi = { B + B[0,0](1-Pi_M) : B in B(M) }, dim_A = m^2
 * (the M^perp scalar is determined by B[0,0], not free). */
static void make_compress_idemp(aic_ucp_kraus *phi, slong d, slong m)
{
    slong r = 1 + (d - m);
    aic_ucp_kraus_init(phi, d, d, r);
    for (slong i = 0; i < m; i++) acb_set_si(acb_mat_entry(phi->K[0], i, i), 1);
    for (slong b = 0; b < d - m; b++)
        acb_set_si(acb_mat_entry(phi->K[1 + b], 0, m + b), 1);
}

/* ---- shared cross-check harness ---- */

/* Pi_A = Delta Delta^dag (the gauge-invariant subspace projector for A). */
static void build_Pi_A(acb_mat_t Pi_A, const acb_mat_t Delta, slong prec)
{
    acb_mat_t Dd;
    acb_mat_init(Dd, acb_mat_ncols(Delta), acb_mat_nrows(Delta));
    acb_mat_conjugate_transpose(Dd, Delta);
    acb_mat_mul(Pi_A, Delta, Dd, prec);
    acb_mat_clear(Dd);
}

/* internal sanity: ||M^dag M - 1||_op (orthonormal columns) <= tol. */
static void check_iso(const char *name, const char *what, const acb_mat_t M,
                      const arb_t tol, slong prec)
{
    slong c = acb_mat_ncols(M);
    acb_mat_t Md, G, one;
    arb_t e;
    acb_mat_init(Md, c, acb_mat_nrows(M));
    acb_mat_init(G, c, c);
    acb_mat_init(one, c, c);
    arb_init(e);
    acb_mat_conjugate_transpose(Md, M);
    acb_mat_mul(G, Md, M, prec);
    acb_mat_one(one);
    acb_mat_sub(G, G, one, prec);
    aic_mat_opnorm(e, G, prec);
    AIC_CHECK_MSG(arb_le(e, tol), "%s: %s columns not orthonormal", name, what);
    arb_clear(e);
    acb_mat_clear(one); acb_mat_clear(G); acb_mat_clear(Md);
}

/* Decompose, assert dims, internal sanity, the five relations + Gamma CP. */
static void check_decomp(const char *name, aic_ucp_kraus *phi,
                         slong expect_dM, slong expect_dA)
{
    arb_t tol, def;
    arb_init(tol);
    arb_init(def);
    set_tol(tol, 1e-12);

    aic_idemp_decomp d;
    aic_idemp_decompose(&d, phi, PREC);

    AIC_CHECK_MSG(d.dim_M == expect_dM, "%s: dim_M=%ld, expected %ld", name,
                  (long) d.dim_M, (long) expect_dM);
    AIC_CHECK_MSG(d.dim_A == expect_dA, "%s: dim_A=%ld, expected %ld", name,
                  (long) d.dim_A, (long) expect_dA);

    slong n = d.n, dM = d.dim_M;

    check_iso(name, "J_M", d.J_M, tol, PREC);
    check_iso(name, "Delta", d.Delta, tol, PREC);
    {   /* Pi_M^2 = Pi_M */
        acb_mat_t P2, diff;
        arb_t e;
        acb_mat_init(P2, n, n);
        acb_mat_init(diff, n, n);
        arb_init(e);
        acb_mat_mul(P2, d.Pi_M, d.Pi_M, PREC);
        acb_mat_sub(diff, P2, d.Pi_M, PREC);
        aic_mat_opnorm(e, diff, PREC);
        AIC_CHECK_MSG(arb_le(e, tol), "%s: Pi_M not idempotent", name);
        arb_clear(e);
        acb_mat_clear(diff); acb_mat_clear(P2);
    }

    aic_idemp_defect_GCD(def, &d, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "%s: Gamma C_M Delta != 1_A", name);
    aic_idemp_defect_DGC(def, &d, phi, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "%s: Delta Gamma C_M != Phi", name);
    aic_idemp_defect_w_hom(def, &d, phi, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "%s: w not a *-homomorphism", name);
    aic_idemp_defect_blockdiag(def, &d, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "%s: Im Delta not block-diagonal", name);
    aic_idemp_defect_gamma_unital(def, &d, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "%s: Gamma not unital", name);

    {   /* Gamma CP — PAIR of certificates (see include/aic_idemp.h):
         *  (i) the C_M o Lambda factorization is CP: Choi(Psi) PSD;
         * (ii) the STORED Gamma is tied to that CP object: ||w Gamma - C_M Lambda|| = 0
         *      (reads d->w AND d->Gamma; w a *-monomorphism, .tex:2088). */
        acb_mat_t C;
        acb_mat_init(C, dM * dM, dM * dM);
        aic_idemp_psi_choi(C, &d, PREC);
        AIC_CHECK_MSG(aic_ucp_is_cp_choi(C, tol, PREC) == 1,
                      "%s: C_M o Lambda factorization not certified CP", name);
        acb_mat_clear(C);
        /* Tie-check tol = 1e-9 (NOT the gauge-clean 1e-12): w Gamma routes through
         * Pi_A = Delta Delta^dag (double-path projector, idempotent only to ~1e-16),
         * so the clean defect floors at 1.6e-12 (noiseless_subsystem). A corrupt
         * Gamma gives O(1), so 1e-9 keeps full teeth. See include/aic_idemp.h
         * aic_idemp_defect_wG_eq_CML for the full precision-floor analysis. */
        arb_t tol_tie;
        arb_init(tol_tie);
        set_tol(tol_tie, 1e-9);
        aic_idemp_defect_wG_eq_CML(def, &d, PREC);
        AIC_CHECK_MSG(arb_le(def, tol_tie),
                      "%s: stored Gamma not tied to CP object (w Gamma != C_M Lambda)",
                      name);
        arb_clear(tol_tie);
    }

    aic_idemp_clear(&d);
    arb_clear(def);
    arb_clear(tol);
}

#endif /* TEST_IDEMP_H */
