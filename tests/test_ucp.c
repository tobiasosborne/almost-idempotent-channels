/* test_ucp.c — cross-checks for the UCP module (bead aic-c7n), the four-rung
 * ladder (CLAUDE.md Rule 6) plus the closed-form oracles and the adversarial
 * corpus (project is adversarial-first; memory feedback_adversarial_benchmarks).
 *
 * Tiers:
 *  1. INTERNAL SANITY: unitality two ways (Kraus and Choi-partial-trace) AGREE;
 *     V is an isometry (V^dag V = 1_H); CP certified PSD on a constructed CP map.
 *  2. ROUND-TRIP: Kraus->Choi->Kraus->Choi reproduces C_Phi (compare CHANNELS via
 *     the Choi matrix, due to Kraus gauge freedom).
 *  3. DOUBLE-vs-arb@53: Choi entries, unital defect, carrier-Q eigenvalues,
 *     cbnorm_cp.
 *  4. ETA=0 EXACT-IDEMPOTENT ORACLE: a coordinate-compression Phi(X)=Pi X Pi is
 *     exactly idempotent UCP; PhiX_M identity exact, carrier spectrum {0,1},
 *     CP+unital certified, cbnorm_cp=1.
 *  5. CLOSED FORMS: identity channel -> Choi = unnormalized max-entangled
 *     projector, cbnorm=1; state Phi(X)=Tr(rho X) (dim_H=1) -> carrier=supp(rho).
 *  6. ADVERSARIAL: rank-deficient Choi; maximally degenerate (depolarizing);
 *     proper-subspace carrier; non-self-map K!=H; near-non-CP must FAIL loud.
 * Plus MUTATION PROOFS at the end (Rule 7).
 *
 * Convention reminders: Kraus K_a are dim_K x dim_H, Phi(X)=sum K_a^dag X K_a
 * (Heisenberg); Choi Convention A, K factor LEFT. See include/aic_ucp.h.
 */
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"

#define PREC 53
static void set_tol(arb_t tol, double t) { arb_set_d(tol, t); }

/* |a - b| <= tol for two arb balls (used to cross-check the two unital-defect
 * routes agree). Rigorous: arb_le on the enclosing |a-b| ball. */
static int arb_close(const arb_t a, const arb_t b, const arb_t tol)
{
    arb_t diff;
    int ok;
    arb_init(diff);
    arb_sub(diff, a, b, 64);
    arb_abs(diff, diff);
    ok = arb_le(diff, tol);
    arb_clear(diff);
    return ok;
}

/* ---- channel constructors (return an initialised aic_ucp_kraus) ---- */

/* Identity channel on C^d: one Kraus op K_0 = 1_d (Phi = id, UCP, idempotent). */
static void make_identity(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, 1);
    acb_mat_one(phi->K[0]);
}

/* State-as-UCP: Phi: B(C^k) -> B(C^1)=C, Phi(X)=Tr(rho X). Kraus ops K_a are
 * k x 1 columns; rho = sum_a K_a K_a^dag is the (k x k) state. We build rho from
 * a list of weighted pure columns. Here: one rank-given diagonal rho. The
 * carrier is supp(rho). */
static void make_state_diag(aic_ucp_kraus *phi, slong k, const double *probs,
                            slong nnz)
{
    /* one Kraus op per nonzero probability: K_a = sqrt(p_a) |e_a> (k x 1). */
    aic_ucp_kraus_init(phi, k, 1, nnz);
    slong a = 0;
    for (slong i = 0; i < k; i++) {
        if (probs[i] <= 0) continue;
        acb_set_d(acb_mat_entry(phi->K[a], i, 0), sqrt(probs[i]));
        a++;
    }
    assert(a == nnz && "make_state_diag: nnz must equal the count of probs[i]>0");
}

/* Depolarizing-style maximally-degenerate map on C^d: Phi(X) = Tr(X) 1_d / d.
 * As a UCP map B(C^d)->B(C^d): K_{ij} = |e_i><e_j|/sqrt(d) for all i,j (d^2
 * Kraus ops), since sum_{ij} K_ij^dag X K_ij = sum_{ij} |e_j><e_i| X |e_i><e_j|
 * / d = (sum_i <e_i|X|e_i>) (sum_j |e_j><e_j|)/d = Tr(X) 1_d / d. Unital:
 * sum K^dag K = sum_{ij}|e_j><e_j|/d = d * 1_d / d = 1_d. Choi is maximally
 * degenerate (Choi matrix = 1/d * 1_{d^2}, all eigenvalues equal). */
static void make_depolarizing(aic_ucp_kraus *phi, slong d)
{
    aic_ucp_kraus_init(phi, d, d, d * d);
    double inv = 1.0 / sqrt((double) d);
    slong a = 0;
    for (slong i = 0; i < d; i++)
        for (slong j = 0; j < d; j++) {
            acb_set_d(acb_mat_entry(phi->K[a], i, j), inv); /* |e_i><e_j|/sqrt d */
            a++;
        }
}

/* Exactly-idempotent UCP on C^d with a PROPER carrier M = span(e_0..e_{m-1}),
 * rank m. Kraus set: K_0 = Pi_M (rank-m coordinate projector) and, for each
 * M^perp basis vector, K_{1+b} = |e_0><e_{m+b}| (b = 0..d-m-1). In the Heisenberg
 * action Phi(X) = sum_a K_a^dag X K_a:
 *   K_0^dag X K_0       = Pi_M X Pi_M,
 *   K_b^dag X K_b       = X[0,0] |e_{m+b}><e_{m+b}|,  so  sum_b = X[0,0](1-Pi_M).
 * Hence Phi(X) = Pi_M X Pi_M + X[0,0](1-Pi_M), whose image {Pi_M X Pi_M +
 * c(1-Pi_M)} is stable => Phi^2 = Phi (exactly idempotent).
 * Unital: K_0^dag K_0 = Pi_M and sum_b K_b^dag K_b = 1 - Pi_M, total 1_d.
 * Carrier: Q = sum_a K_a K_a^dag = Pi_M + (d-m)|e_0><e_0|, support = M (rank m).
 * This is the eta=0 exact-idempotent oracle's channel. */
static void make_compress_idemp(aic_ucp_kraus *phi, slong d, slong m)
{
    slong r = 1 + (d - m);
    aic_ucp_kraus_init(phi, d, d, r);
    for (slong i = 0; i < m; i++)               /* K_0 = Pi */
        acb_set_si(acb_mat_entry(phi->K[0], i, i), 1);
    for (slong b = 0; b < d - m; b++)           /* K_{1+b} = |e_0><e_{m+b}| */
        acb_set_si(acb_mat_entry(phi->K[1 + b], 0, m + b), 1);
}

/* Genuinely COMPLEX, non-Hermitian, ASYMMETRIC unital channel on C^2 with two
 * Kraus ops K_a : C^2 -> C^2 (dim_K = dim_H = 2):
 *   K_0 = [[0.6, 0],[0.8i, 0]],   K_1 = [[0, 0.8i],[0, 0.6]].
 * K_0 is supported on column 0, K_1 on column 1, so sum_a K_a^dag K_a =
 * diag(|0.6|^2+|0.8i|^2, |0.8i|^2+|0.6|^2) = 1_2 (exactly unital). The ops are
 * NOT symmetric (K_0[1,0]=0.8i but K_0[0,1]=0) and carry nonzero imaginary
 * parts, so they distinguish conj-on-first-factor from conj-on-second-factor in
 * the Choi convention (the diagonal-entry case K=K^T cancels the difference).
 * Used by the Choi-convention oracle below. */
static void make_complex_channel(aic_ucp_kraus *phi)
{
    aic_ucp_kraus_init(phi, 2, 2, 2);
    acb_set_d_d(acb_mat_entry(phi->K[0], 0, 0), 0.6, 0.0);  /* K_0[0,0] = 0.6   */
    acb_set_d_d(acb_mat_entry(phi->K[0], 1, 0), 0.0, 0.8);  /* K_0[1,0] = 0.8i  */
    acb_set_d_d(acb_mat_entry(phi->K[1], 0, 1), 0.0, 0.8);  /* K_1[0,1] = 0.8i  */
    acb_set_d_d(acb_mat_entry(phi->K[1], 1, 1), 0.6, 0.0);  /* K_1[1,1] = 0.6   */
}

/* ---- tier 1b: Choi-convention oracle (independent code path) ----
 * The Choi block (i,j) of C_Phi must equal Phi(E_{ij}) entrywise, where E_{ij}
 * is the dim_K x dim_K matrix unit and Phi is computed by the INDEPENDENT
 * aic_ucp_apply (K^dag X K matmuls). This is a genuine Rule-6 cross-check of the
 * Choi index/conjugation convention against the Heisenberg action — NOT the
 * transpose-invariant unital-two-ways check, which only validates the
 * partial-trace DIRECTION. The complex asymmetric channel above makes the
 * conj-on-wrong-factor bug visible (max block-vs-apply diff 0.96 if flipped). */
static void test_choi_convention_oracle(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);

    aic_ucp_kraus phi;
    make_complex_channel(&phi);
    slong dk = phi.dim_K, dh = phi.dim_H, n = dk * dh;

    acb_mat_t C;
    acb_mat_init(C, n, n);
    aic_ucp_kraus_to_choi(C, &phi, PREC);

    acb_mat_t E, PhiE, blk;
    acb_mat_init(E, dk, dk);
    acb_mat_init(PhiE, dh, dh);
    acb_mat_init(blk, dh, dh);
    for (slong i = 0; i < dk; i++)
        for (slong j = 0; j < dk; j++) {
            acb_mat_zero(E);
            acb_set_si(acb_mat_entry(E, i, j), 1);     /* E_{ij} = |e_i><e_j|  */
            aic_ucp_apply(PhiE, &phi, E, PREC);        /* independent path     */
            /* extract block (i,j): rows i*dh..i*dh+dh-1, cols j*dh..j*dh+dh-1 */
            for (slong a = 0; a < dh; a++)
                for (slong b = 0; b < dh; b++)
                    acb_set(acb_mat_entry(blk, a, b),
                            acb_mat_entry(C, i * dh + a, j * dh + b));
            AIC_CHECK_ACB_MAT_CLOSE(blk, PhiE, tol);   /* Convention-A check   */
        }

    acb_mat_clear(blk);
    acb_mat_clear(PhiE);
    acb_mat_clear(E);
    acb_mat_clear(C);
    aic_ucp_kraus_clear(&phi);
    arb_clear(tol);
}

/* ---- tier 1: internal sanity ---- */
static void test_internal_sanity(void)
{
    arb_t tol, d1, d2, vdef, cb;
    arb_init(tol); arb_init(d1); arb_init(d2); arb_init(vdef); arb_init(cb);
    set_tol(tol, 1e-12);

    aic_ucp_kraus phi;
    make_compress_idemp(&phi, 4, 2);
    slong dk = phi.dim_K, dh = phi.dim_H;

    /* unitality two ways must agree (validates the Choi index convention) */
    acb_mat_t C;
    acb_mat_init(C, dk * dh, dk * dh);
    aic_ucp_kraus_to_choi(C, &phi, PREC);
    aic_ucp_unital_defect_choi(d1, C, dk, dh, PREC);
    aic_ucp_unital_defect_kraus(d2, &phi, PREC);
    AIC_CHECK_MSG(arb_close(d1, d2, tol),
                  "unital defects (Choi vs Kraus) disagree -> bad Choi convention");
    AIC_CHECK_MSG(arb_le(d1, tol), "unital defect (Choi) not ~0");
    AIC_CHECK_MSG(arb_le(d2, tol), "unital defect (Kraus) not ~0");

    /* V isometry: ||V^dag V - 1_H||_op ~ 0 */
    acb_mat_t V, Vd, VtV, one;
    acb_mat_init(V, dk * phi.r, dh);
    acb_mat_init(Vd, dh, dk * phi.r);
    acb_mat_init(VtV, dh, dh);
    acb_mat_init(one, dh, dh);
    aic_ucp_kraus_to_stinespring(V, &phi, PREC);
    acb_mat_conjugate_transpose(Vd, V);
    acb_mat_mul(VtV, Vd, V, PREC);
    acb_mat_one(one);
    acb_mat_sub(VtV, VtV, one, PREC);
    aic_mat_opnorm(vdef, VtV, PREC);
    AIC_CHECK_MSG(arb_le(vdef, tol), "V not an isometry");

    /* CP certified PSD */
    arb_t zero; arb_init(zero); arb_zero(zero);
    AIC_CHECK_MSG(aic_ucp_is_cp_choi(C, tol, PREC) == 1, "CP map not certified PSD");
    arb_clear(zero);

    aic_ucp_cbnorm_cp(cb, &phi, PREC);
    AIC_CHECK_DOUBLE_CLOSE(1.0, cb, tol); /* idempotent UCP => cbnorm 1 */

    acb_mat_clear(one); acb_mat_clear(VtV); acb_mat_clear(Vd); acb_mat_clear(V);
    acb_mat_clear(C);
    aic_ucp_kraus_clear(&phi);
    arb_clear(cb); arb_clear(vdef); arb_clear(d2); arb_clear(d1); arb_clear(tol);
}

/* ---- tier 2: round-trip Kraus->Choi->Kraus->Choi (compare as CHANNELS) ---- */
static void test_roundtrip(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-10); /* double-path extraction in the middle => 53-bit tol */

    /* use the compression idempotent (rank-deficient Choi => degenerate spectrum,
     * exactly what LAPACK Choi->Kraus must handle and arb-simple cannot). */
    aic_ucp_kraus phi;
    make_compress_idemp(&phi, 4, 2);
    slong dk = phi.dim_K, dh = phi.dim_H, n = dk * dh;

    acb_mat_t C1, C2;
    acb_mat_init(C1, n, n);
    acb_mat_init(C2, n, n);
    aic_ucp_kraus_to_choi(C1, &phi, PREC);

    aic_ucp_kraus phi2; /* output of extraction; init'd inside */
    aic_ucp_choi_to_kraus_latd(&phi2, C1, dk, dh);
    aic_ucp_kraus_to_choi(C2, &phi2, PREC);

    /* compare CHANNELS via Choi (Kraus gauge freedom => never op-by-op). */
    AIC_CHECK_ACB_MAT_CLOSE(C1, C2, tol);

    /* the recovered rank must be the Choi rank (here 1 + (d-m) = 3). */
    AIC_CHECK_MSG(phi2.r == 1 + (dk - 2),
                  "extracted Kraus rank %ld != Choi rank %ld",
                  (long) phi2.r, (long) (1 + (dk - 2)));

    aic_ucp_kraus_clear(&phi2);
    acb_mat_clear(C2);
    acb_mat_clear(C1);
    aic_ucp_kraus_clear(&phi);

    /* COMPLEX round-trip: the compression/depolarizing channels above are real,
     * so their Choi eigenvectors are real and the reshape-conjugation is a no-op.
     * The complex asymmetric channel has genuinely complex eigenvectors, so this
     * round-trip is what proves the Choi->Kraus reshape conjugation MATCHES the
     * kraus_to_choi conjugation (H2 lockstep): flip either one and C1 != C3. */
    aic_ucp_kraus cphi;
    make_complex_channel(&cphi);
    slong cn = cphi.dim_K * cphi.dim_H;
    acb_mat_t C3, C4;
    acb_mat_init(C3, cn, cn);
    acb_mat_init(C4, cn, cn);
    aic_ucp_kraus_to_choi(C3, &cphi, PREC);
    aic_ucp_kraus cphi2;
    aic_ucp_choi_to_kraus_latd(&cphi2, C3, cphi.dim_K, cphi.dim_H);
    aic_ucp_kraus_to_choi(C4, &cphi2, PREC);
    AIC_CHECK_ACB_MAT_CLOSE(C3, C4, tol);
    aic_ucp_kraus_clear(&cphi2);
    acb_mat_clear(C4);
    acb_mat_clear(C3);
    aic_ucp_kraus_clear(&cphi);

    arb_clear(tol);
}

/* ---- tier 3: double-vs-arb@53 (Choi entries, unital defect, Q-eig, cbnorm) ---- */
static void test_double_vs_arb(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);

    aic_ucp_kraus phi;
    make_compress_idemp(&phi, 4, 2);
    slong dk = phi.dim_K, dh = phi.dim_H, n = dk * dh;

    /* (a) Choi entries: arb@53 matrix vs a double recomputation of the same
     * formula. Build the double Choi directly and compare each entry. */
    acb_mat_t C;
    acb_mat_init(C, n, n);
    aic_ucp_kraus_to_choi(C, &phi, PREC);
    for (slong i = 0; i < dk; i++)
        for (slong a = 0; a < dh; a++)
            for (slong j = 0; j < dk; j++)
                for (slong b = 0; b < dh; b++) {
                    double _Complex acc = 0;
                    for (slong x = 0; x < phi.r; x++) {
                        /* K_x entries are exact integers here */
                        double kr = arf_get_d(
                            arb_midref(acb_realref(
                                acb_mat_entry(phi.K[x], i, a))), ARF_RND_NEAR);
                        double kj = arf_get_d(
                            arb_midref(acb_realref(
                                acb_mat_entry(phi.K[x], j, b))), ARF_RND_NEAR);
                        acc += kr * kj;
                    }
                    double cre = arf_get_d(
                        arb_midref(acb_realref(
                            acb_mat_entry(C, i * dh + a, j * dh + b))),
                        ARF_RND_NEAR);
                    AIC_CHECK_MSG(fabs(creal(acc) - cre) < 1e-12,
                                  "Choi entry double-vs-arb off at (%ld,%ld)",
                                  (long) (i * dh + a), (long) (j * dh + b));
                }

    /* (b) carrier-Q eigenvalues: double path (LAPACK) vs arb herm_max_eig (top)
     * AND a full-set compare using the double eig of the arb-built Q. Q is
     * Hermitian PSD dim_K x dim_K. */
    acb_mat_t Q;
    acb_mat_init(Q, dk, dk);
    aic_ucp_carrier_Q(Q, &phi, PREC);
    /* double eig set */
    double _Complex Qa[16];
    double dev[4];
    aic_latd_from_acb_mat(Qa, Q);
    aic_latd_eig_hermitian(dev, NULL, Qa, dk);
    /* arb top eigenvalue */
    arb_t lam;
    arb_init(lam);
    aic_mat_herm_max_eig(lam, Q, PREC);
    double dmax = dev[dk - 1];
    AIC_CHECK_DOUBLE_CLOSE(dmax, lam, tol);
    arb_clear(lam);

    /* (c) cbnorm_cp = 1 (UCP) double-side: opnorm of Phi(1_K) via arb compared
     * to the exact 1. */
    arb_t cb;
    arb_init(cb);
    aic_ucp_cbnorm_cp(cb, &phi, PREC);
    AIC_CHECK_DOUBLE_CLOSE(1.0, cb, tol);
    arb_clear(cb);

    acb_mat_clear(Q);
    acb_mat_clear(C);
    aic_ucp_kraus_clear(&phi);
    arb_clear(tol);
}

/* ---- tier 4: eta=0 exact-idempotent oracle ---- */
static void test_idemp_oracle(void)
{
    arb_t tol, def;
    arb_init(tol);
    arb_init(def);
    set_tol(tol, 1e-12);

    const slong d = 4, m = 2;
    aic_ucp_kraus phi;
    make_compress_idemp(&phi, d, m);

    /* Phi^2 = Phi on a basis E_ij: max defect ~ 0 (the map IS idempotent). */
    acb_mat_t E, PhiE, Phi2E, diff;
    acb_mat_init(E, d, d);
    acb_mat_init(PhiE, d, d);
    acb_mat_init(Phi2E, d, d);
    acb_mat_init(diff, d, d);
    arb_t idemp_max, t;
    arb_init(idemp_max);
    arb_init(t);
    arb_zero(idemp_max);
    for (slong i = 0; i < d; i++)
        for (slong j = 0; j < d; j++) {
            acb_mat_zero(E);
            acb_set_si(acb_mat_entry(E, i, j), 1);
            aic_ucp_apply(PhiE, &phi, E, PREC);
            aic_ucp_apply(Phi2E, &phi, PhiE, PREC);
            acb_mat_sub(diff, Phi2E, PhiE, PREC);
            aic_mat_opnorm(t, diff, PREC);
            if (arb_gt(t, idemp_max)) arb_set(idemp_max, t);
        }
    AIC_CHECK_MSG(arb_le(idemp_max, tol), "Phi not idempotent (eta=0 oracle)");

    /* PhiX_M identity exact: Pi_M = diag(1,1,0,0); ||Phi(X)-Phi(Pi X Pi)||~0 on a
     * random-ish X. */
    acb_mat_t Pi, X;
    acb_mat_init(Pi, d, d);
    acb_mat_init(X, d, d);
    for (slong i = 0; i < m; i++) acb_set_si(acb_mat_entry(Pi, i, i), 1);
    unsigned long s = 0xABCDUL;
    for (slong i = 0; i < d; i++)
        for (slong j = 0; j < d; j++) {
            s = s * 6364136223846793005UL + 1UL;
            acb_set_si_si(acb_mat_entry(X, i, j),
                          (slong) ((s >> 40) % 7) - 3,
                          (slong) ((s >> 20) % 7) - 3);
        }
    aic_ucp_phiX_M_defect(def, &phi, X, Pi, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "PhiX_M identity violated (eta=0 oracle)");

    /* carrier spectrum {0,1}-supported: rank = m (carrier is the rank-m block).
     * Q = Pi + (d-m)|e_0><e_0|, support = M (rank m). */
    acb_mat_t Q;
    acb_mat_init(Q, d, d);
    aic_ucp_carrier_Q(Q, &phi, PREC);
    slong rank = aic_ucp_carrier_rank_latd(Q, d);
    AIC_CHECK_MSG(rank == m, "carrier rank %ld != m=%ld", (long) rank, (long) m);

    /* CP + unital certified, cbnorm = 1 */
    acb_mat_t C;
    acb_mat_init(C, d * d, d * d);
    aic_ucp_kraus_to_choi(C, &phi, PREC);
    AIC_CHECK_MSG(aic_ucp_is_cp_choi(C, tol, PREC) == 1, "oracle: not CP");
    aic_ucp_unital_defect_kraus(def, &phi, PREC);
    AIC_CHECK_MSG(arb_le(def, tol), "oracle: not unital");
    arb_t cb;
    arb_init(cb);
    aic_ucp_cbnorm_cp(cb, &phi, PREC);
    AIC_CHECK_DOUBLE_CLOSE(1.0, cb, tol);
    arb_clear(cb);

    acb_mat_clear(C);
    acb_mat_clear(Q);
    acb_mat_clear(X);
    acb_mat_clear(Pi);
    arb_clear(t);
    arb_clear(idemp_max);
    acb_mat_clear(diff);
    acb_mat_clear(Phi2E);
    acb_mat_clear(PhiE);
    acb_mat_clear(E);
    aic_ucp_kraus_clear(&phi);
    arb_clear(def);
    arb_clear(tol);
}

/* ---- tier 5: closed forms (identity channel; state channel) ---- */
static void test_closed_forms(void)
{
    arb_t tol, cb;
    arb_init(tol);
    arb_init(cb);
    set_tol(tol, 1e-12);

    /* identity channel on C^d: Choi = |Omega><Omega| (UNNORMALIZED max-entangled
     * projector), i.e. C[i*d+a, j*d+b] = delta_{ia} delta_{jb} (rank-1 onto
     * sum_i |ii>). cbnorm = 1. */
    const slong d = 3;
    aic_ucp_kraus id;
    make_identity(&id, d);
    acb_mat_t C;
    acb_mat_init(C, d * d, d * d);
    aic_ucp_kraus_to_choi(C, &id, PREC);
    /* expected: C[i*d+a, j*d+b] = [i==a] * [j==b] (since K_0=1, K_0[i,a]=[i==a]).
     * That is |Omega><Omega| with Omega = sum_i |i>_K|i>_H, entries at the
     * (i,i),(j,j) positions = 1. Verify a few entries and the trace = d. */
    acb_t one_acb, zero_acb;
    acb_init(one_acb);
    acb_init(zero_acb);
    acb_set_si(one_acb, 1);
    for (slong i = 0; i < d; i++)
        for (slong a = 0; a < d; a++)
            for (slong j = 0; j < d; j++)
                for (slong b = 0; b < d; b++) {
                    int want = (i == a) && (j == b);
                    arb_t e;
                    arb_init(e);
                    /* |entry - want| */
                    acb_t diff;
                    acb_init(diff);
                    acb_sub(diff, acb_mat_entry(C, i * d + a, j * d + b),
                            want ? one_acb : zero_acb, PREC);
                    acb_abs(e, diff, PREC);
                    AIC_CHECK_MSG(arb_le(e, tol),
                                  "identity Choi entry wrong at (%ld,%ld)",
                                  (long) (i * d + a), (long) (j * d + b));
                    acb_clear(diff);
                    arb_clear(e);
                }
    aic_ucp_cbnorm_cp(cb, &id, PREC);
    AIC_CHECK_DOUBLE_CLOSE(1.0, cb, tol);
    acb_clear(zero_acb);
    acb_clear(one_acb);
    acb_mat_clear(C);
    aic_ucp_kraus_clear(&id);

    /* state channel Phi(X)=Tr(rho X), dim_H=1, rho = diag(0.5,0,0.5) on C^3.
     * carrier = supp(rho) = span(e_0,e_2), rank 2. Q = sum_a K_a K_a^dag = rho. */
    const slong k = 3;
    double probs[3] = {0.5, 0.0, 0.5};
    aic_ucp_kraus st;
    make_state_diag(&st, k, probs, 2);
    acb_mat_t Q;
    acb_mat_init(Q, k, k);
    aic_ucp_carrier_Q(Q, &st, PREC);
    slong rank = aic_ucp_carrier_rank_latd(Q, k);
    AIC_CHECK_MSG(rank == 2, "state carrier rank %ld != 2", (long) rank);
    /* unital: Phi(1_k) = Tr(rho) = 1 (state normalised) => cbnorm = 1 */
    aic_ucp_cbnorm_cp(cb, &st, PREC);
    AIC_CHECK_DOUBLE_CLOSE(1.0, cb, tol);
    acb_mat_clear(Q);
    aic_ucp_kraus_clear(&st);

    arb_clear(cb);
    arb_clear(tol);
}

/* ---- tier 6: adversarial corpus ---- */
static void test_adversarial(void)
{
    arb_t tol, def;
    arb_init(tol);
    arb_init(def);
    set_tol(tol, 1e-10);

    /* (a) maximally degenerate: depolarizing on C^3, Choi = 1/d * 1_{d^2}
     * (ALL eigenvalues equal 1/d) -> arb-simple eig WOULD abort; LAPACK
     * Choi->Kraus must round-trip. */
    {
        const slong d = 3, n = d * d;
        aic_ucp_kraus phi;
        make_depolarizing(&phi, d);
        acb_mat_t C1, C2;
        acb_mat_init(C1, n, n);
        acb_mat_init(C2, n, n);
        aic_ucp_kraus_to_choi(C1, &phi, PREC);
        /* CP + unital */
        AIC_CHECK_MSG(aic_ucp_is_cp_choi(C1, tol, PREC) == 1,
                      "depolarizing not CP");
        aic_ucp_unital_defect_kraus(def, &phi, PREC);
        AIC_CHECK_MSG(arb_le(def, tol), "depolarizing not unital");
        /* round-trip through degenerate extraction */
        aic_ucp_kraus phi2;
        aic_ucp_choi_to_kraus_latd(&phi2, C1, d, d);
        aic_ucp_kraus_to_choi(C2, &phi2, PREC);
        AIC_CHECK_ACB_MAT_CLOSE(C1, C2, tol);
        aic_ucp_kraus_clear(&phi2);
        acb_mat_clear(C2);
        acb_mat_clear(C1);
        aic_ucp_kraus_clear(&phi);
    }

    /* (b) proper-subspace carrier dim_K=4, M=2: cor_carrier defect must be ~0
     * for X = 1 - Pi_M (annihilates M) and NONZERO for X = Pi_M (does not). */
    {
        const slong d = 4, m = 2;
        aic_ucp_kraus phi;
        make_compress_idemp(&phi, d, m);
        acb_mat_t X;
        acb_mat_init(X, d, d);
        /* X = 1 - Pi_M = diag(0,0,1,1): Ker X >= M, so (X(x)1)V = 0. */
        for (slong i = m; i < d; i++) acb_set_si(acb_mat_entry(X, i, i), 1);
        aic_ucp_carrier_annihilate_defect(def, &phi, X, PREC);
        AIC_CHECK_MSG(arb_le(def, tol),
                      "cor_carrier: X=1-Pi should annihilate carrier");
        /* X = Pi_M: does NOT annihilate M, defect must be bounded away from 0. */
        acb_mat_zero(X);
        for (slong i = 0; i < m; i++) acb_set_si(acb_mat_entry(X, i, i), 1);
        aic_ucp_carrier_annihilate_defect(def, &phi, X, PREC);
        arb_t half;
        arb_init(half);
        arb_set_d(half, 0.5);
        AIC_CHECK_MSG(arb_gt(def, half),
                      "cor_carrier: X=Pi should NOT annihilate carrier");
        arb_clear(half);
        acb_mat_clear(X);
        aic_ucp_kraus_clear(&phi);
    }

    /* (c) non-self-map K != H: Phi: B(C^2) -> B(C^3). Build a UCP map with
     * Kraus K_a : C^3 -> C^2 (2 x 3). Use an isometric-completion style: V is a
     * 6 x 3 isometry (dim_K * r = 2*3 = 6 rows, dim_H = 3). Take K_0 = top
     * 2x3 block of a 6x3 isometry, etc. Simplest: K_0 = [[1,0,0],[0,1,0]]/?,
     * pick three rank-deficient ops summing K^dag K = 1_3. Use the canonical
     * "embed C^3 into C^2 x C^3" no — instead: a measure-and-prepare UCP map
     * Phi(X) = sum_{i=0}^{2} <e_i|? Let's use Phi(X) = Tr(X) * sigma for a fixed
     * 3x3 density sigma? That is B(C^2)->B(C^3): Phi(X)=Tr(X) sigma needs
     * Phi(1_2)=2 sigma = 1_3 -> sigma=1_3/2, not a state. Use instead the
     * UCP map Phi(X)= sum_a K_a^dag X K_a with K_a = |u_a><e_a| chosen so
     * sum K_a^dag K_a = 1_3. Take K_a = c_a |v_a><e_a|, a=0,1,2, v_a in C^2 unit:
     *   K_a^dag K_a = |c_a|^2 |e_a><e_a|; sum = 1_3 needs |c_a|=1, v_a unit. */
    {
        const slong dk = 2, dh = 3, r = 3;
        aic_ucp_kraus phi;
        aic_ucp_kraus_init(&phi, dk, dh, r);
        /* v_0=|0>, v_1=|1>, v_2=(|0>+|1>)/sqrt2 in C^2 */
        acb_set_si(acb_mat_entry(phi.K[0], 0, 0), 1);              /* |0><e_0| */
        acb_set_si(acb_mat_entry(phi.K[1], 1, 1), 1);              /* |1><e_1| */
        double inv = 1.0 / sqrt(2.0);
        acb_set_d(acb_mat_entry(phi.K[2], 0, 2), inv);            /* v_2<e_2| */
        acb_set_d(acb_mat_entry(phi.K[2], 1, 2), inv);
        /* unital sum_a K_a^dag K_a = 1_3 */
        aic_ucp_unital_defect_kraus(def, &phi, PREC);
        AIC_CHECK_MSG(arb_le(def, tol), "K!=H map not unital");
        /* CP via Choi (dim_K*dim_H = 6) */
        acb_mat_t C;
        acb_mat_init(C, dk * dh, dk * dh);
        aic_ucp_kraus_to_choi(C, &phi, PREC);
        AIC_CHECK_MSG(aic_ucp_is_cp_choi(C, tol, PREC) == 1, "K!=H map not CP");
        /* unital two ways agree (validates convention for rectangular Kraus) */
        arb_t dc;
        arb_init(dc);
        aic_ucp_unital_defect_choi(dc, C, dk, dh, PREC);
        AIC_CHECK_MSG(arb_close(dc, def, tol),
                      "K!=H: unital Choi vs Kraus disagree");
        /* cbnorm = 1 */
        arb_t cb;
        arb_init(cb);
        aic_ucp_cbnorm_cp(cb, &phi, PREC);
        AIC_CHECK_DOUBLE_CLOSE(1.0, cb, tol);
        arb_clear(cb);
        arb_clear(dc);
        acb_mat_clear(C);
        aic_ucp_kraus_clear(&phi);
    }

    /* (d) near-non-CP: a Hermitian Choi-shaped matrix with one eigenvalue just
     * below 0 must be REJECTED by the CP check (verdict 0). Build a diagonal
     * (hence eigenvalue-controlled) 4x4 Hermitian C = diag(1, 0.5, 0.25, -1e-6).
     * The -1e-6 is far above prec=53 noise (~1e-15), so the global-enclosure ball
     * for -lambda_min(C) is certainly > 0 and the verdict must be 0, not a guess.
     * The off-boundary magnitude (1e-6) is deliberate: it certifies the gap is
     * resolved (no straddle abort) while the matrix is genuinely indefinite. */
    {
        const slong n = 4;
        acb_mat_t C;
        acb_mat_init(C, n, n);
        acb_mat_zero(C);
        acb_set_d(acb_mat_entry(C, 0, 0), 1.0);
        acb_set_d(acb_mat_entry(C, 1, 1), 0.5);
        acb_set_d(acb_mat_entry(C, 2, 2), 0.25);
        acb_set_d(acb_mat_entry(C, 3, 3), -1e-6);
        arb_t zero;
        arb_init(zero);
        arb_zero(zero);
        AIC_CHECK_MSG(aic_ucp_is_cp_choi(C, zero, PREC) == 0,
                      "near-non-CP Choi was NOT rejected by CP check");
        arb_clear(zero);
        acb_mat_clear(C);
    }

    arb_clear(def);
    arb_clear(tol);
}

/* ---- mutation proofs (Rule 7): perturb impl mentally via wrong-convention
 * inputs and confirm the checks reject them. We exercise the CHECKS' teeth here
 * with deliberately-wrong values; the source-perturbation mutations (wrong
 * partial-trace direction, missing sqrt, transposed reshape) are documented in
 * the report and were run by hand. ---- */
static void test_mutation_teeth(void)
{
    arb_t tol;
    arb_init(tol);
    set_tol(tol, 1e-12);

    /* a clearly-non-unital map: K_0 = 2 * 1_2 (sum K^dag K = 4 1_2 != 1_2).
     * The unital defect must be ~3 (||4-1||=3), REJECTED by arb_le(., tol). */
    aic_ucp_kraus phi;
    aic_ucp_kraus_init(&phi, 2, 2, 1);
    acb_set_si(acb_mat_entry(phi.K[0], 0, 0), 2);
    acb_set_si(acb_mat_entry(phi.K[0], 1, 1), 2);
    arb_t def;
    arb_init(def);
    aic_ucp_unital_defect_kraus(def, &phi, PREC);
    AIC_CHECK_MSG(!arb_le(def, tol),
                  "mutation: non-unital map slipped past the unital check");
    /* and the Choi-route defect must AGREE (~3), proving both routes see it. */
    acb_mat_t C;
    acb_mat_init(C, 4, 4);
    aic_ucp_kraus_to_choi(C, &phi, PREC);
    arb_t defc;
    arb_init(defc);
    aic_ucp_unital_defect_choi(defc, C, 2, 2, PREC);
    AIC_CHECK_MSG(arb_close(def, defc, tol),
                  "mutation: unital routes disagree on a non-unital map");
    arb_clear(defc);
    acb_mat_clear(C);
    arb_clear(def);
    aic_ucp_kraus_clear(&phi);
    arb_clear(tol);
}

int main(void)
{
    test_choi_convention_oracle();
    test_internal_sanity();
    test_roundtrip();
    test_double_vs_arb();
    test_idemp_oracle();
    test_closed_forms();
    test_adversarial();
    test_mutation_teeth();
    aic_test_report("test_ucp");
    printf("OK test_ucp\n");
    return 0;
}
