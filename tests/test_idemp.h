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

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic_test.h"
#include "aic/aic_ucp.h"

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

/* Weighted noiseless subsystem on H = C^dL (x) C^dE: the make_noiseless_subsystem
 * generalization with a NON-UNIFORM multiplicity weight sigma = diag(w_0..w_{dE-1})
 * (a density matrix on E, sum_k w_k = 1):
 *   Phi(X) = Tr_E( (1_L (x) sigma) X ) (x) 1_E.
 * Heisenberg Kraus K_{jk} = sqrt(w_j) (1_L (x) |e_j><e_k|) (dE^2 ops); H row I =
 * a*dE + j (left-major, L the major factor, matching make_noiseless_subsystem and
 * the W_j row convention t = a*dE + c, include/aic/aic_idemp.h:179). EXACTLY
 * idempotent and unital (verified, /tmp probe: ||Phi^2-Phi||~4e-17, ||sum K^dag K
 * - 1||~5e-17): Phi^2=Phi because sum_j w_j (1_L(x)|e_k><e_k|) collapses to a
 * projection-like reset of E, and sum_j w_j = 1 gives unitality. Img Phi ~ B(C^dL),
 * dim_A = dL^2; M = H (Q=1), dim_M = dL*dE.
 *
 * WHY this exists (aic-ynu I3 hostile-review gap): make_noiseless_subsystem mixes E
 * MAXIMALLY (sigma = 1_E/dE) so its extracted gamma_j is UNIFORM ((1/2)1) — every
 * existing prop_Gamma oracle has a uniform gamma_j, so none exercises the data-driven
 * gamma_j solve (a hardcoded "gamma = 1/dE" shortcut would pass them all). With
 * w = (0.8, 0.2) the prop_Gamma extraction returns gamma_0 = diag(0.2, 0.8) (the
 * NON-uniform conditional-expectation density; measured, see test_gamma_kraus.c) —
 * a GREEN test where the correctly-solved gamma_j is genuinely non-uniform. */
static void make_ns_weighted(aic_ucp_kraus *phi, slong dL, slong dE, const double *w)
{
    slong n = dL * dE;
    aic_ucp_kraus_init(phi, n, n, dE * dE);
    slong op = 0;
    for (slong j = 0; j < dE; j++)
        for (slong k = 0; k < dE; k++) {
            double s = sqrt(w[j]);
            for (slong a = 0; a < dL; a++)
                acb_set_d(acb_mat_entry(phi->K[op], a * dE + j, a * dE + k), s);
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

/* Build a FIXED n x n complex unitary V = exp(i H) for a fixed Hermitian H with
 * a real (basis-mixing) part AND a complex (phase) part, so V is non-real and
 * not basis-aligned. exp(i H) is unitary to ~machine precision (arb), enough to
 * keep a conjugated idempotent EXACTLY idempotent at the 1e-12 oracle gate.
 * H[p,q] for p<q: real 0.4 on the (p,p+1) super-diagonal (basis mixing) + a
 * complex 0.3 i on (p, p+2) when it exists; H Hermitian by symmetrising. */
__attribute__((unused))
static void make_fixed_unitary(acb_mat_t V, slong n, slong prec)
{
    acb_mat_t H, iH;
    acb_mat_init(H, n, n);
    acb_mat_init(iH, n, n);
    acb_mat_zero(H);
    for (slong p = 0; p + 1 < n; p++) {
        acb_set_d(acb_mat_entry(H, p, p + 1), 0.4);     /* real coupling   */
        acb_set_d(acb_mat_entry(H, p + 1, p), 0.4);     /* Hermitian       */
    }
    for (slong p = 0; p + 2 < n; p++) {
        acb_set_d_d(acb_mat_entry(H, p, p + 2), 0.0, 0.3);   /* + 0.3 i     */
        acb_set_d_d(acb_mat_entry(H, p + 2, p), 0.0, -0.3);  /* Hermitian   */
    }
    for (slong p = 0; p < n; p++)
        acb_set_d(acb_mat_entry(H, p, p), 0.1 * (double) (p + 1)); /* phases */
    /* iH = i * H */
    {
        acb_t I_unit;
        acb_init(I_unit);
        acb_onei(I_unit);
        acb_mat_scalar_mul_acb(iH, H, I_unit, prec);
        acb_clear(I_unit);
    }
    acb_mat_exp(V, iH, prec); /* V = exp(i H), unitary */
    acb_mat_clear(iH);
    acb_mat_clear(H);
}

/* Conjugate an exactly-idempotent self-map `base` (on C^n) by a FIXED complex
 * unitary V: K'_a = V^dag K_a V, so Phi'(X) = V^dag Phi(V X V^dag) V =
 * Ad_{V^dag} o Phi o Ad_V. Phi' is still EXACTLY idempotent (Ad_V Ad_{V^dag}=id)
 * and unital, but its image algebra V^dag (Img Phi) V is generally NOT transpose-
 * closed for complex V — which gives the from_idemp reshape a real tooth (a
 * TRANSPOSED reshape B_k[i,j]=Delta[j*n+i,k] would span a DIFFERENT subspace,
 * making the gauge-invariant defects O(1); the original 7 channels all had
 * transpose-closed images and could not catch it). `out` is init'd here. */
__attribute__((unused))
static void make_conjugated(aic_ucp_kraus *out, const aic_ucp_kraus *base,
                            slong prec)
{
    slong n = base->dim_H;
    AIC_CHECK_MSG(base->dim_K == n,
                  "make_conjugated: base must be a self-map (dim_K==dim_H)");
    acb_mat_t V, Vd, tmp;
    acb_mat_init(V, n, n);
    acb_mat_init(Vd, n, n);
    acb_mat_init(tmp, n, n);
    make_fixed_unitary(V, n, prec);
    acb_mat_conjugate_transpose(Vd, V);
    aic_ucp_kraus_init(out, n, n, base->r);
    for (slong a = 0; a < base->r; a++) {
        acb_mat_mul(tmp, Vd, base->K[a], prec);   /* V^dag K_a       */
        acb_mat_mul(out->K[a], tmp, V, prec);     /* V^dag K_a V     */
    }
    acb_mat_clear(tmp);
    acb_mat_clear(Vd);
    acb_mat_clear(V);
}

/* Convex mix of make_compress_idemp(d,m) with its unitary conjugate (mirrors the
 * test_corner make_eta_family Kraus-union pattern): Phi_t = (1-t) Phi_compress +
 * t (V^dag Phi_compress V), realised as the Kraus union {sqrt(1-t) K_a} U
 * {sqrt(t) V^dag K_a V}. This is a GENUINELY eta>0 (eta proxy ~ 6.5 t for d=5,m=3),
 * NON-COMMUTATIVE, OBLIQUE almost-idempotent channel whose associated algebra A
 * has genuine 1d delta-projections with a dim S_{P,Q}=2 corner (P = projector onto
 * span(e1,e2), Q = |e1><e1|): the regime lem_extension needs (Ha^Q_{P,P} a
 * homomorphism on a >1-dim corner at eta>0). At t=0.02: eta ~ 1.30e-2, P/Q in-A
 * residuals ~ 9.2e-3 / 4.8e-3, dim S_{P,Q}=2, dim S_Q=1, lem_PQ_Hilb Gram
 * off-diagonal G[0,1] ~ 1.7e-6, ||G-I|| ~ 1.3e-5 (trivially in-basin). The bare
 * aic_mat_opnorm ||G-I|| guard SIGABRTs on this corner (the near-zero off-diagonal
 * false-fails the Gram Hermiticity check, bead aic-2yo); aic_corner_ha now uses the
 * certified mid+radius upper bound, so it runs. `out` is aic_ucp_kraus_init'd here. */
__attribute__((unused))
static void make_mixconj(aic_ucp_kraus *out, slong d, slong m, double t, slong prec)
{
    aic_ucp_kraus base, conj;
    make_compress_idemp(&base, d, m);
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

/* ---- shared cross-check harness ---- */

/* Pi_A = Delta Delta^dag (the gauge-invariant subspace projector for A).
 * __attribute__((unused)): the channel constructors above are reused by
 * tests that do not drive the full idemp-relation harness (e.g. test_ecstar.c);
 * marking the harness-only helpers keeps those TUs warning-clean. */
__attribute__((unused))
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
__attribute__((unused))
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
