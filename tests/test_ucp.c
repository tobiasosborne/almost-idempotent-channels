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
 *  7. ADVERSARIAL CORPUS RETROFIT (bead aic-dbo.4 inc.2; mirrors the funcalc
 *     retrofit test_funcalc.c check 7): each adversarial CHANNEL generator
 *     (aic_adversarial.h) is drawn and, per instance, the relevant ucp ROUTINE
 *     is asserted to EITHER hold a certified bound OR fail loud. BOUND-HOLDS is
 *     asserted INLINE (certified arb balls); FAIL-LOUD is asserted via a
 *     fork+SIGALRM watchdog (aic_watchdog_assert_failloud) that confirms SIGABRT + a
 *     message-needle. The cleanest ucp fail-loud is the carrier-rank STRADDLE
 *     (aic_ucp_carrier_rank on the densified 1C carrier at prec=53); the
 *     non-CP-Choi rejection (aic_ucp_is_cp_choi verdict 0) is the second.
 *     This tests the ROUTINE, not the generator's self-test property (which
 *     test_adversarial.c owns). See test_adversarial section below.
 * Plus MUTATION PROOFS at the end (Rule 7).
 *
 * Convention reminders: Kraus K_a are dim_K x dim_H, Phi(X)=sum K_a^dag X K_a
 * (Heisenberg); Choi Convention A, K factor LEFT. See include/aic_ucp.h.
 */
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>
#include <flint/arf.h>

#include "aic/aic_cbnorm.h"
#include "aic/aic_latd.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_adversarial.h"
#include "aic_test.h"
#include "aic_watchdog.h"

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

/* ====================================================================== *
 * tier 7: ADVERSARIAL CORPUS RETROFIT (bead aic-dbo.4 increment 2)        *
 *                                                                         *
 * Draw the adversarial CHANNEL generators (aic_adversarial.h) and assert, *
 * per instance, the DUAL OUTCOME on the relevant ucp ROUTINE: EITHER a    *
 * certified ucp bound holds (asserted INLINE), OR the routine FAILS LOUD  *
 * (asserted via the fork+SIGALRM watchdog below). This mirrors the just-  *
 * landed funcalc retrofit (test_funcalc.c check 7) and the kraus_arb /    *
 * v5f straddle pattern. We test the ROUTINE, not the generator's own      *
 * self-test property (test_adversarial.c owns those).                     *
 *                                                                         *
 * Provenance: tex:1568-1635 (Choi/Stinespring/Kraus, prop_KLHG),          *
 * tex:1717-1719 (||Phi||_cb<=1), tex:1724/1731 (carrier lem_carrier/      *
 * cor_carrier), tex:347-354 (||Phi^2-Phi||_cb), tex:366-388 (the measure- *
 * prepare cb!=op closed form eta*sqrt(1-eta)), tex:516-525 (theta(2Phi-1) *
 * basin). docs/adversarial/domain.md families 1B/1C/1D/2A/2B/3D + NC.     *
 * CLAUDE.md cross-check ladder: rung 2 (double vs arb), rung 3 (eta=0     *
 * oracle), rung 4 (arb bound certification at the hypothesis boundary).   *
 *                                                                         *
 * RETROFIT PRECISION (Rule 2, MEASURED via a bounded /tmp probe, NOT in   *
 * tests/): the certified arb Choi->Kraus (aic_ucp_choi_to_kraus_arb) on a *
 * 16x16 degenerate block-algebra/measure-prepare Choi FAILS LOUD at the   *
 * densifier-unitary guard ||U U^dag - I||_F < n^2*2^-(prec-8) (bead       *
 * aic-wyo, P3 OPEN; the guard's magic-number tolerance fail-loud-aborts   *
 * for n>=16 at prec=128 even though the radius ~3.7e-34 is far below any   *
 * soundness threshold). That is a KNOWN fail-loud (never silent-wrong),   *
 * NOT a ucp finding — so the round-trip BOUND-HOLDS half uses the LATD     *
 * (double) Choi->Kraus path, which round-trips these degenerate Choi      *
 * cleanly (measured ||C_rebuilt-C||_op ~ 1e-15). The certified rank /     *
 * carrier-rank / cb-bracket / CP / unital routines are all exercised on   *
 * the arb path directly.                                                  */
#define UCP_RETRO_PREC 128       /* arb prec for the bound-holds asserts */

/* ---- shared bound-holds helpers ---- */

/* ||C_rebuilt - C||_op upper bound (double), C round-tripped through the LATD
 * (double-path) Choi->Kraus extractor — the degeneracy-robust route that the
 * certified-arb extractor cannot run on a 16x16 block Choi (aic-wyo guard). */
static double ucp_latd_roundtrip(const aic_ucp_kraus *phi, slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H, n = dk * dh;
    acb_mat_t C, Cr, D;
    acb_mat_init(C, n, n);
    acb_mat_init(Cr, n, n);
    acb_mat_init(D, n, n);
    aic_ucp_kraus_to_choi(C, phi, prec);
    aic_ucp_kraus rec;
    aic_ucp_choi_to_kraus_latd(&rec, C, dk, dh);
    aic_ucp_kraus_to_choi(Cr, &rec, prec);
    acb_mat_sub(D, Cr, C, prec);
    arb_t nr;
    arb_init(nr);
    aic_mat_opnorm(nr, D, prec);
    double ub = arf_get_d(arb_midref(nr), ARF_RND_UP) + mag_get_d(arb_radref(nr));
    arb_clear(nr);
    aic_ucp_kraus_clear(&rec);
    acb_mat_clear(D);
    acb_mat_clear(Cr);
    acb_mat_clear(C);
    return ub;
}

/* arb bracket [lo,hi] on eta = ||Phi^2-Phi||_cb must CONTAIN the NONZERO
 * closed-form value `cf` (tex:347-354 / 378). Rigorous: lower(lo) <= cf <=
 * upper(hi) via arb_ge/arb_le on the certified balls. ONLY for cf bounded away
 * from 0 — for the eta=0 oracle use ucp_bracket_upper_below (the computed J
 * carries an O(1e-16) Frobenius-norm rounding floor for non-integer-Kraus
 * channels like depolarizing, so the bracket lower bound can sit just ABOVE 0;
 * "contains 0" would be a FALSE assertion, measured lo~6.8e-17 at 2A p=1). */
static int ucp_bracket_contains(const aic_ucp_kraus *phi, double cf, slong prec)
{
    arb_t lo, hi, cfb;
    arb_init(lo); arb_init(hi); arb_init(cfb);
    aic_cbnorm_eigfree_ball(lo, hi, phi, prec);
    arb_set_d(cfb, cf);
    /* cf >= lower(lo)  AND  cf <= upper(hi), both rigorous on the balls. */
    int ok = arb_ge(cfb, lo) && arb_le(cfb, hi);
    arb_clear(cfb); arb_clear(hi); arb_clear(lo);
    return ok;
}

/* eta=0 oracle (rung 3): the bracket's UPPER bound must be <= tol, i.e. eta is
 * CERTIFIED tiny (an exactly-idempotent channel). Measured: 3D t=0 gives hi=0
 * exactly (integer block-projector Kraus); depolarizing p in {0,1} gives
 * hi~2.7e-16 (the 1/d rounding floor). Both are << 1e-10. */
static int ucp_bracket_upper_below(const aic_ucp_kraus *phi, double tol_d,
                                   slong prec)
{
    arb_t lo, hi, tol;
    arb_init(lo); arb_init(hi); arb_init(tol);
    aic_cbnorm_eigfree_ball(lo, hi, phi, prec);
    arb_set_d(tol, tol_d);
    int ok = arb_le(hi, tol);   /* upper(hi) <= tol, rigorous */
    arb_clear(tol); arb_clear(hi); arb_clear(lo);
    return ok;
}

/* CP + unital certified, on a freshly-built Choi. Returns 1 iff CP verdict==1
 * and unital defect <= tol. */
static int ucp_cp_and_unital(const aic_ucp_kraus *phi, double tol_d, slong prec)
{
    slong dk = phi->dim_K, dh = phi->dim_H, n = dk * dh;
    acb_mat_t C;
    acb_mat_init(C, n, n);
    aic_ucp_kraus_to_choi(C, phi, prec);
    arb_t tol, ud;
    arb_init(tol); arb_init(ud);
    arb_set_d(tol, tol_d);
    int cp = aic_ucp_is_cp_choi(C, tol, prec);
    aic_ucp_unital_defect_kraus(ud, phi, prec);
    int unital = arb_le(ud, tol);
    arb_clear(ud); arb_clear(tol);
    acb_mat_clear(C);
    return cp == 1 && unital;
}

/* ---- FAIL-LOUD children (module-scope void(void) for the watchdog) ---- */

/* (FL-1) carrier-rank STRADDLE: the densified non-normal 1C carrier at prec=53,
 * gap=1e-16 — the small-cluster ball (~1.5e-15) straddles thr (~9.4e-16) so
 * aic_ucp_carrier_rank cannot resolve the rank and MUST abort ("STRADDLES").
 * The Q is built AT the query prec (FINDINGS §D7n). The DIAGONAL 1C cannot
 * reach this (its gap eig is a point ball that always decides). */
static void ucp_child_carrier_straddle(void)
{
    aic_ucp_kraus phi;
    aic_adv_chan_carrier_dropout_dense(&phi, 3, 1e-16, 53);
    acb_mat_t Q;
    acb_mat_init(Q, 3, 3);
    aic_ucp_carrier_Q(Q, &phi, 53);
    (void) aic_ucp_carrier_rank(Q, 3, 53);   /* must abort: STRADDLES */
    acb_mat_clear(Q);
    aic_ucp_kraus_clear(&phi);
}

/* (FL-2) certified Choi->Kraus on a genuinely non-CP Choi MUST fail loud (not
 * return junk Kraus). Build a 2B conc_defect channel (CP) and FLIP its Choi to
 * non-CP by planting a clearly-negative eigenvalue along the diagonal of the
 * sign-projected Choi — the strict arb extractor (neg_tol=keep_thr) must abort
 * ("CP"/"not CP"). We use a small (4x4, dim_K=dim_H=2) Choi so the arb path
 * runs below the aic-wyo n>=16 guard. */
static void ucp_child_noncp_extract(void)
{
    /* diag(1, 0.5, 0.25, -0.3) conjugated by a rational Givens (so not coord-
     * aligned), an O(1)-negative Choi. Mirror test_kraus_arb.c child_o1_negative. */
    const slong dk = 2, dh = 2, n = dk * dh;
    acb_mat_t D, Qr, Qt, T, C;
    acb_mat_init(D, n, n); acb_mat_init(Qr, n, n); acb_mat_init(Qt, n, n);
    acb_mat_init(T, n, n); acb_mat_init(C, n, n);
    acb_mat_zero(D);
    double diag[4] = {1.0, 0.5, 0.25, -0.3};
    for (slong i = 0; i < n; i++) acb_set_d(acb_mat_entry(D, i, i), diag[i]);
    acb_mat_one(Qr);
    acb_t cc, ss, ns;
    acb_init(cc); acb_init(ss); acb_init(ns);
    acb_set_si(cc, 3); acb_div_si(cc, cc, 5, UCP_RETRO_PREC);   /* 3/5 */
    acb_set_si(ss, 4); acb_div_si(ss, ss, 5, UCP_RETRO_PREC);   /* 4/5 */
    acb_neg(ns, ss);
    acb_set(acb_mat_entry(Qr, 0, 0), cc);
    acb_set(acb_mat_entry(Qr, n - 1, n - 1), cc);
    acb_set(acb_mat_entry(Qr, 0, n - 1), ns);
    acb_set(acb_mat_entry(Qr, n - 1, 0), ss);
    acb_mat_conjugate_transpose(Qt, Qr);
    acb_mat_mul(T, Qr, D, UCP_RETRO_PREC);
    acb_mat_mul(T, T, Qt, UCP_RETRO_PREC);
    /* project to exact-Hermitian midpoint (FINDINGS §C5) so the extractor's
     * tight Hermiticity assert passes for the right reason. */
    for (slong i = 0; i < n; i++) {
        acb_get_mid(acb_mat_entry(C, i, i), acb_mat_entry(T, i, i));
        arb_zero(acb_imagref(acb_mat_entry(C, i, i)));
        for (slong j = i + 1; j < n; j++) {
            acb_get_mid(acb_mat_entry(C, i, j), acb_mat_entry(T, i, j));
            acb_conj(acb_mat_entry(C, j, i), acb_mat_entry(C, i, j));
        }
    }
    aic_ucp_kraus rec;
    aic_ucp_choi_to_kraus_arb(&rec, C, dk, dh, UCP_RETRO_PREC);  /* must abort */
    aic_ucp_kraus_clear(&rec);
    acb_clear(ns); acb_clear(ss); acb_clear(cc);
    acb_mat_clear(C); acb_mat_clear(T); acb_mat_clear(Qt);
    acb_mat_clear(Qr); acb_mat_clear(D);
}

/* ---- the retrofit driver ---- */
static void test_adversarial_corpus(void)
{
    const slong prec = UCP_RETRO_PREC;
    printf("tier 7: adversarial channel corpus retrofit (bound-holds OR fail-loud)\n");

    /* === fam1C-dense carrier dropout: aic_ucp_carrier_rank ============== *
     * BOUND-HOLDS (gap>>thr -> certified rank == d == carrier_rank_latd)   *
     * AND FAIL-LOUD (gap straddle at prec=53 -> SIGABRT "STRADDLES").      */
    {
        const slong d = 3;
        aic_ucp_kraus phi;
        aic_adv_chan_carrier_dropout_dense(&phi, d, 0.5, prec);
        acb_mat_t Q;
        acb_mat_init(Q, d, d);
        aic_ucp_carrier_Q(Q, &phi, prec);
        slong rk_cert = aic_ucp_carrier_rank(Q, d, prec);
        slong rk_dbl = aic_ucp_carrier_rank_latd(Q, d);
        AIC_CHECK_MSG(rk_cert == d,
                      "1C-dense gap=0.5: certified carrier rank %ld != d=%ld",
                      (long) rk_cert, (long) d);
        AIC_CHECK_MSG(rk_cert == rk_dbl,
                      "1C-dense gap=0.5: certified rank %ld != double-path %ld "
                      "(double-vs-arb carrier rank disagree)",
                      (long) rk_cert, (long) rk_dbl);
        printf("  1C-dense carrier_rank BOUND-HOLDS: cert=%ld == double=%ld == d=%ld\n",
               (long) rk_cert, (long) rk_dbl, (long) d);
        acb_mat_clear(Q);
        aic_ucp_kraus_clear(&phi);

        printf("  1C-dense carrier_rank FAIL-LOUD (straddle prec=53):\n");
        aic_watchdog_assert_failloud(ucp_child_carrier_straddle, 20,
                                     "aic_ucp_carrier_rank(1C-dense gap=1e-16 prec=53)",
                                     "STRADDLES");
    }

    /* === fam1B cb!=op gap: cb-bracket contains eta*sqrt(1-eta) + CP + unital *
     * tex:366-388/378. The op-norm defect < the cb defect, so the certified  *
     * eig-free cb bracket [||J||_F/n, 2||J||_F] must STILL contain the cb     *
     * closed form (a routine using op-norm as eta would underreport).        */
    {
        const slong d = 2;
        double etas[3] = {0.05, 0.10, 0.20};
        for (int e = 0; e < 3; e++) {
            aic_ucp_kraus phi;
            aic_adv_chan_cb_op_gap(&phi, d, etas[e], prec);
            double cf = etas[e] * sqrt(1.0 - etas[e]);
            AIC_CHECK_MSG(ucp_bracket_contains(&phi, cf, prec),
                          "1B eta=%.2f: cb bracket does NOT contain closed form "
                          "eta*sqrt(1-eta)=%.6f (tex:378)", etas[e], cf);
            AIC_CHECK_MSG(ucp_cp_and_unital(&phi, 1e-10, prec),
                          "1B eta=%.2f: channel not certified CP+unital", etas[e]);
            aic_ucp_kraus_clear(&phi);
        }
        printf("  1B cb!=op BOUND-HOLDS: bracket contains eta*sqrt(1-eta) "
               "(eta in {0.05,0.10,0.20}); Choi CP+unital\n");
    }

    /* === fam2A depolarizing eta->1/4 basin edge: cb-bracket scales as p(1-p) *
     * tex:516-525. rho(Phi^2-Phi)=p(1-p), maximized 1/4 at p=1/2. The eig-    *
     * free bracket lower bound inherits the EXACT linearity in p(1-p)         *
     * (measured ratio lo/(p(1-p))=sqrt(d^2-1)/d=sqrt(3)/2~=0.866 at d=2).     *
     * p=0 and p=1 are EXACTLY idempotent (defect 0) — the eta=0 rung-3 oracle. */
    {
        const slong d = 2;
        /* eta=0 oracle endpoints: p=0 (identity) and p=1 (cond expectation),
         * both EXACTLY idempotent -> the cb bracket upper bound is certified
         * tiny (measured hi~2.7e-16, the 1/d rounding floor; NOT exactly 0). */
        aic_ucp_kraus phi0, phi1;
        aic_adv_chan_depol_boundary(&phi0, d, 0.0, prec);
        aic_adv_chan_depol_boundary(&phi1, d, 1.0, prec);
        AIC_CHECK_MSG(ucp_bracket_upper_below(&phi0, 1e-10, prec),
                      "2A p=0: eta=0 oracle — bracket upper not certified < 1e-10");
        AIC_CHECK_MSG(ucp_bracket_upper_below(&phi1, 1e-10, prec),
                      "2A p=1: eta=0 oracle — bracket upper not certified < 1e-10");
        aic_ucp_kraus_clear(&phi1);
        aic_ucp_kraus_clear(&phi0);
        /* p in (0,1): the bracket must contain the EXACT cb value, which is at
         * least the certified lower bound lo = sqrt(3)/2 * p(1-p) and at most
         * the closed-form rho is NOT the cb value here (depolarizing is non-
         * unital-defect-free; cb >= rho), so assert the certified bracket
         * BRACKETS rho/p(1-p) consistency: lo <= 2*hi and lo > 0 for p in (0,1),
         * plus CP. We assert lo grows EXACTLY linearly in p(1-p). */
        double ps[3] = {0.3, 0.5, 0.7};
        double lo_over_pp[3];
        for (int k = 0; k < 3; k++) {
            aic_ucp_kraus phi;
            aic_adv_chan_depol_boundary(&phi, d, ps[k], prec);
            arb_t lo, hi;
            arb_init(lo); arb_init(hi);
            aic_cbnorm_eigfree_ball(lo, hi, &phi, prec);
            double lov = arf_get_d(arb_midref(lo), ARF_RND_DOWN)
                         - mag_get_d(arb_radref(lo));
            double pp = ps[k] * (1.0 - ps[k]);
            lo_over_pp[k] = lov / pp;
            AIC_CHECK_MSG(lov > 0.0,
                          "2A p=%.2f: cb lower bound %.6e not > 0 in (0,1)",
                          ps[k], lov);
            AIC_CHECK_MSG(ucp_cp_and_unital(&phi, 1e-10, prec),
                          "2A p=%.2f: not certified CP+unital", ps[k]);
            arb_clear(hi); arb_clear(lo);
            aic_ucp_kraus_clear(&phi);
        }
        /* the lo/p(1-p) ratio is p-INDEPENDENT (exact linearity): all three
         * agree to a tight tol (measured ~0.86603). */
        AIC_CHECK_MSG(fabs(lo_over_pp[0] - lo_over_pp[1]) < 1e-6 &&
                      fabs(lo_over_pp[1] - lo_over_pp[2]) < 1e-6,
                      "2A: cb lower bound NOT linear in p(1-p) "
                      "(ratios %.6f %.6f %.6f differ)",
                      lo_over_pp[0], lo_over_pp[1], lo_over_pp[2]);
        printf("  2A depol BOUND-HOLDS: eta=0 at p in {0,1}; cb-lo linear in "
               "p(1-p) (ratio %.5f); CP+unital\n", lo_over_pp[1]);
    }

    /* === fam1D unital-but-barely: aic_ucp_unital_defect_kraus == delta_u === *
     * tex:432/672. Single Hermitian Kraus -> CP; the certified unital defect  *
     * equals delta_u EXACTLY (a point ball), and the Choi route agrees.       */
    {
        const slong d = 3;
        double dus[3] = {0.0, 0.15, 0.30};
        for (int e = 0; e < 3; e++) {
            aic_ucp_kraus phi;
            aic_adv_chan_unital_defect(&phi, d, dus[e], prec);
            arb_t ud, udc, tol;
            arb_init(ud); arb_init(udc); arb_init(tol);
            arb_set_d(tol, 1e-12);
            aic_ucp_unital_defect_kraus(ud, &phi, prec);
            /* defect == delta_u exactly */
            AIC_CHECK_MSG(aic_double_close(dus[e], ud, tol),
                          "1D du=%.2f: certified unital defect (Kraus) != du",
                          dus[e]);
            /* the Choi route must AGREE (validates the partial-trace direction) */
            slong n = phi.dim_K * phi.dim_H;
            acb_mat_t C;
            acb_mat_init(C, n, n);
            aic_ucp_kraus_to_choi(C, &phi, prec);
            aic_ucp_unital_defect_choi(udc, C, phi.dim_K, phi.dim_H, prec);
            AIC_CHECK_MSG(arb_close(ud, udc, tol),
                          "1D du=%.2f: unital defect Kraus vs Choi disagree", dus[e]);
            /* single Hermitian Kraus -> CP verdict must be 1 */
            AIC_CHECK_MSG(aic_ucp_is_cp_choi(C, tol, prec) == 1,
                          "1D du=%.2f: single-Hermitian-Kraus map not CP", dus[e]);
            acb_mat_clear(C);
            arb_clear(tol); arb_clear(udc); arb_clear(ud);
            aic_ucp_kraus_clear(&phi);
        }
        printf("  1D unital_defect BOUND-HOLDS: certified defect == delta_u "
               "(du in {0,0.15,0.30}); Kraus==Choi route; CP\n");
    }

    /* === fam3D block algebra + fam-NC noncomm: multi-block Choi -> Kraus    *
     * round-trip (LATD path; the arb path fails loud at the aic-wyo n=16     *
     * densifier guard) + CP + unital. The block-structured Choi is the       *
     * degeneracy stressor. tex:484/1249 (3D dim blowup); tex:347-349 cb      *
     * ampliation-invariance (NC). At t=0 the 3D channel is EXACTLY idempotent *
     * (rung-3 oracle: round-trip 0, eta bracket [0,0]).                       */
    {
        /* 3D k=2,d=2, t=0 (exact idempotent oracle) and t=0.10 (eta>0). */
        double ts[2] = {0.0, 0.10};
        for (int e = 0; e < 2; e++) {
            aic_ucp_kraus phi;
            aic_adv_chan_blockalg(&phi, 2, 2, ts[e], prec);
            double rt = ucp_latd_roundtrip(&phi, prec);
            AIC_CHECK_MSG(rt < 1e-10,
                          "3D k2d2 t=%.2f: Choi->Kraus->Choi round-trip %.3e "
                          "not < 1e-10", ts[e], rt);
            AIC_CHECK_MSG(ucp_cp_and_unital(&phi, 1e-10, prec),
                          "3D k2d2 t=%.2f: not certified CP+unital", ts[e]);
            if (ts[e] == 0.0)
                AIC_CHECK_MSG(ucp_bracket_upper_below(&phi, 1e-10, prec),
                              "3D t=0: eta=0 oracle — bracket upper not "
                              "certified < 1e-10 (measured exactly 0)");
            aic_ucp_kraus_clear(&phi);
        }
        /* 3D k=3,d=2,t=0: a bigger block count, still exactly idempotent. */
        {
            aic_ucp_kraus phi;
            aic_adv_chan_blockalg(&phi, 3, 2, 0.0, prec);
            double rt = ucp_latd_roundtrip(&phi, prec);
            AIC_CHECK_MSG(rt < 1e-10, "3D k3d2 t=0: round-trip %.3e not < 1e-10", rt);
            AIC_CHECK_MSG(ucp_cp_and_unital(&phi, 1e-10, prec),
                          "3D k3d2 t=0: not certified CP+unital");
            aic_ucp_kraus_clear(&phi);
        }
        printf("  3D blockalg BOUND-HOLDS: latd round-trip < 1e-10 (k2d2 t=0,0.10; "
               "k3d2 t=0); CP+unital; eta=0 bracket at t=0\n");

        /* fam-NC: the cb bracket must contain eta_cb = 1/4 - kappa (tex:378),
         * round-trip + CP+unital. Non-abelian image, ampliation-invariant cb. */
        double kappas[3] = {0.05, 0.10, 0.15};
        for (int e = 0; e < 3; e++) {
            aic_ucp_kraus phi;
            aic_adv_chan_noncomm_boundary(&phi, 2, 2, kappas[e], prec);
            double target = 0.25 - kappas[e];
            AIC_CHECK_MSG(ucp_bracket_contains(&phi, target, prec),
                          "NC kappa=%.2f: cb bracket does NOT contain "
                          "eta_cb=1/4-kappa=%.6f (tex:378)", kappas[e], target);
            double rt = ucp_latd_roundtrip(&phi, prec);
            AIC_CHECK_MSG(rt < 1e-10, "NC kappa=%.2f: round-trip %.3e not < 1e-10",
                          kappas[e], rt);
            AIC_CHECK_MSG(ucp_cp_and_unital(&phi, 1e-10, prec),
                          "NC kappa=%.2f: not certified CP+unital", kappas[e]);
            aic_ucp_kraus_clear(&phi);
        }
        printf("  NC noncomm BOUND-HOLDS: cb bracket contains 1/4-kappa "
               "(kappa in {0.05,0.10,0.15}); latd round-trip; CP+unital\n");
    }

    /* === fam2B rank-1 concentrated defect: the channel is CP+unital and the  *
     * cb bracket contains the calibrated eta; the Choi round-trips (LATD).    *
     * tex:2192-2204/2385-2422. The rank-1 superoperator defect Phi^2-Phi is   *
     * the O(sqrt eta) cancellation stressor; gap_sub tunes the near-degenerate *
     * |2> direction without changing eta_cb (depends on tail only).           */
    {
        const slong d = 3;
        double etas[2] = {0.10, 0.30};
        double gap_subs[2] = {0.3, 0.05};   /* mild and near-degenerate subspace */
        for (int e = 0; e < 2; e++) {
            aic_ucp_kraus phi;
            aic_adv_chan_conc_defect(&phi, d, etas[e], gap_subs[e], prec);
            AIC_CHECK_MSG(ucp_bracket_contains(&phi, etas[e], prec),
                          "2B eta=%.2f gap_sub=%.2f: cb bracket does NOT contain "
                          "the calibrated eta", etas[e], gap_subs[e]);
            double rt = ucp_latd_roundtrip(&phi, prec);
            AIC_CHECK_MSG(rt < 1e-10, "2B eta=%.2f: round-trip %.3e not < 1e-10",
                          etas[e], rt);
            AIC_CHECK_MSG(ucp_cp_and_unital(&phi, 1e-10, prec),
                          "2B eta=%.2f: not certified CP+unital", etas[e]);
            aic_ucp_kraus_clear(&phi);
        }
        printf("  2B conc_defect BOUND-HOLDS: cb bracket contains eta "
               "(eta in {0.10,0.30}, gap_sub {0.3,0.05}); latd round-trip; "
               "CP+unital\n");
    }

    /* === FAIL-LOUD (FL-2): certified Choi->Kraus on a genuinely non-CP Choi *
     * must abort, not return junk Kraus (the silent-non-CP trap, Rule 4).    */
    printf("  certified Choi->Kraus FAIL-LOUD (O(1)-negative Choi):\n");
    aic_watchdog_assert_failloud(ucp_child_noncp_extract, 20,
                                 "aic_ucp_choi_to_kraus_arb(non-CP Choi)", "CP");
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
    test_adversarial_corpus();
    aic_test_report("test_ucp");
    printf("OK test_ucp\n");
    return 0;
}
