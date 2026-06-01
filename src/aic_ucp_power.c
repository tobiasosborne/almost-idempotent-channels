/* aic_ucp_power.c — channel power Phi^k and Choi-rank compression (bead aic-xxk,
 * module "ucp"). THIN wrappers over already-tested primitives (CLAUDE.md Law 2:
 * reuse, do not reinvent): aic_ucp_compose (Phi o Phi), aic_ucp_kraus_to_choi +
 * aic_ucp_choi_to_kraus_latd (the minimal-Choi-rank re-extraction).
 *
 * WHY THESE TWO LIVE TOGETHER (CLAUDE.md Rule 3, "all bugs are deep").
 *   aic_ucp_compose returns the Kraus SET { L_b K_a }, multiplying the operator
 *   counts: compose(Phi,Phi) of a rank-r channel has r^2 Kraus operators, and
 *   naive k-fold repeated composition blows the running count up as r^k. For
 *   dephasing(d) (r=d) that is d^k operators — d=4,k=5 is 1024 ops, all but
 *   <= dim_H^2 of them linearly dependent (gauge-redundant). That is not a
 *   correctness bug (the CHANNEL is right) but it is an exponential-memory
 *   trap. The fix is COMPRESSION: every UCP map of Choi rank rho has a Kraus rep
 *   with exactly rho operators (rho <= dim_K*dim_H), recoverable by rebuilding
 *   the Choi matrix and re-extracting its minimal PSD eigenbasis. So
 *   aic_ucp_power COMPRESSES AFTER EACH compose, keeping the running Kraus count
 *   at the Choi rank (<= dim_H^2) instead of r^k. The two utilities are paired
 *   precisely so the power can stay polynomial in memory.
 *
 * COMPRESSION (aic_ucp_compress) is exactly the Choi round-trip. The Choi matrix
 *   C_Phi (.tex:1575, Convention A in aic_ucp.h) is a complete invariant of the
 *   channel: two Kraus reps give the same Choi iff they are the same map. So
 *     Phi  --kraus_to_choi-->  C_Phi  --choi_to_kraus_latd-->  Phi'
 *   produces Phi' with the SAME channel (same Choi) but with the MINIMAL number
 *   of Kraus operators — exactly the count of nonzero eigenvalues of C_Phi, i.e.
 *   the Choi rank (aic_ucp_latd.c keeps eigenpairs with lambda > QETLAB
 *   threshold (dim_K*dim_H)*eps_mach*||C||_F). The channel is preserved up to the
 *   latd extraction tolerance (that QETLAB threshold ~ 1e-15 ||C||_F): tiny
 *   numerically-zero eigenpairs are dropped, which is the compression. Kraus reps
 *   are unique only up to a unitary gauge, so the extracted operators are NOT the
 *   originals — callers must compare AS CHANNELS (rebuild Choi), the discipline
 *   the tests follow.
 *
 * POWER (aic_ucp_power). Phi^k in the Heisenberg picture is the k-fold Choi-Effros
 *   composition (the channel is a UCP self-map of B(H), .tex:342 reminds the
 *   product on Img Phi is the Choi-Effros star, but a CHANNEL POWER is plain map
 *   composition Phi o ... o Phi — we are iterating the map, not multiplying inside
 *   its image). aic_ucp_compose(out, phi, psi) returns Phi o Psi (aic_ucp.h:151,
 *   verified: Kraus { psi.K[b] @ phi.K[a] }, action Phi(Psi(X))). For a POWER both
 *   factors are Phi so the direction is moot (Phi o Phi = Phi^2 either way); we
 *   pass acc as BOTH arguments. k=0 is the identity channel (one Kraus op K_0 =
 *   1_{dim_H}); k=1 is a copy of phi; k>=2 is compress(compose(acc, phi)) iterated.
 *   REQUIRES phi to be a self-map (dim_K == dim_H) for k != 1 — composition and the
 *   k=0 identity both need a square endomorphism; asserted fail-loud (Rule 4).
 *
 * eta=0 ORACLE (CLAUDE.md Rule 6, the cleanest cross-check). For an EXACTLY
 *   idempotent Phi (dephasing, cond_expectation, depolarizing at p=1) Phi^2 = Phi
 *   exactly, so Phi^k = Phi for every k>=1 — Choi(power(Phi,k)) == Choi(Phi) to
 *   machine zero. test_ucp_power.c asserts this for k=2,3,5.
 *
 * Number path. Both routines are on the certified arb/acb path EXCEPT the Choi
 *   re-extraction, which routes through aic_ucp_choi_to_kraus_latd (LAPACK double
 *   path) because the Choi spectrum of a rank-deficient channel is degenerate and
 *   the certified arb eigensolver needs a simple spectrum (aic_ucp.h:44-50, bead
 *   aic-w4o.1). The Choi BUILD and every compose stay on acb at the given prec.
 */
#include <assert.h>

#include <flint/acb_mat.h>

#include "aic/aic_ucp.h"

/* aic_ucp.h: minimal-Kraus-rank (Choi-rank) form of phi. Build C_Phi, then
 * re-extract its minimal PSD eigenbasis. `out` is aic_ucp_kraus_init'd HERE by
 * aic_ucp_choi_to_kraus_latd (to the recovered Choi rank); the CALLER clears it.
 * Channel is preserved up to the latd extraction tolerance (QETLAB threshold
 * (dim_K*dim_H)*eps_mach*||C||_F): numerically-zero eigenpairs are dropped, which
 * is the compression. out->r <= dim_K*dim_H (the Choi rank). */
void aic_ucp_compress(aic_ucp_kraus *out, const aic_ucp_kraus *phi, slong prec)
{
    assert(out != NULL && phi != NULL && "aic_ucp_compress: null argument");

    slong dk = phi->dim_K, dh = phi->dim_H;
    slong n = dk * dh;

    acb_mat_t C;
    acb_mat_init(C, n, n);
    aic_ucp_kraus_to_choi(C, phi, prec);          /* C_Phi (Convention A) */
    /* Re-extract the MINIMAL PSD eigenbasis: inits `out` to the Choi rank and
     * fills it. This IS the compression — only lambda > threshold survive. */
    aic_ucp_choi_to_kraus_latd(out, C, dk, dh);
    acb_mat_clear(C);

    assert(out->r <= n &&
           "aic_ucp_compress: recovered rank exceeds dim_K*dim_H (impossible "
           "for a Choi rank — extraction is broken)");
}

/* aic_ucp.h: Phi^k, k >= 0. k=0 -> identity channel (dim_H x dim_H, r=1,
 * K_0 = 1_{dim_H}); k=1 -> a copy of phi; k>=2 -> compress(compose(acc, phi))
 * iterated. `out` is aic_ucp_kraus_init'd HERE (caller clears). REQUIRES
 * dim_K == dim_H for k != 1 (self-map needed for the identity and for compose). */
void aic_ucp_power(aic_ucp_kraus *out, const aic_ucp_kraus *phi, slong k,
                   slong prec)
{
    assert(out != NULL && phi != NULL && "aic_ucp_power: null argument");
    assert(k >= 0 && "aic_ucp_power: k must be >= 0 (fail loud, Rule 4)");

    /* k=1 is the only case that does NOT require a square map: a copy of phi
     * with its own (dim_K, dim_H). Every other case (identity at k=0, the
     * compose chain at k>=2) needs an endomorphism dim_K == dim_H. */
    if (k != 1)
        assert(phi->dim_K == phi->dim_H &&
               "aic_ucp_power: Phi^k for k != 1 needs a self-map (dim_K == "
               "dim_H) — composition and the k=0 identity are endomorphisms");

    if (k == 0) {
        /* Identity channel on C^{dim_H}: one Kraus op K_0 = 1_{dim_H}. The
         * identity is the empty product, the neutral element of composition;
         * Phi^0 = id and id o Phi = Phi recovers k=1 from one more compose. */
        slong d = phi->dim_H;
        aic_ucp_kraus_init(out, d, d, 1);
        acb_mat_one(out->K[0]);                   /* K_0 = 1_d */
        return;
    }

    if (k == 1) {
        /* A copy of phi (same Kraus set, same channel). NOT compressed: k=1 is
         * "return phi unchanged" by contract; compression is the caller's job
         * via aic_ucp_compress if they want the minimal rep. */
        aic_ucp_kraus_init(out, phi->dim_K, phi->dim_H, phi->r);
        for (slong a = 0; a < phi->r; a++)
            acb_mat_set(out->K[a], phi->K[a]);
        return;
    }

    /* k >= 2: start the accumulator at Phi (compressed so the FIRST compose
     * already runs on the minimal rep), then repeatedly acc <- compress(acc o
     * Phi). COMPRESS-EACH-STEP is the whole point (see file docstring): without
     * it the running Kraus count is r^k (exponential); with it it stays at the
     * Choi rank <= dim_H^2 every iteration. */
    aic_ucp_kraus acc;
    aic_ucp_compress(&acc, phi, prec);            /* acc = Phi (minimal rep) */

    for (slong step = 2; step <= k; step++) {
        /* acc o Phi: direction is moot for a power (acc and phi are powers of
         * the same Phi, which commute), but we keep the fixed order compose(acc,
         * phi) = acc o phi for definiteness. */
        aic_ucp_kraus composed;
        aic_ucp_compose(&composed, &acc, phi, prec);  /* r grows to acc.r*phi.r */

        aic_ucp_kraus compressed;
        aic_ucp_compress(&compressed, &composed, prec); /* back to Choi rank */
        aic_ucp_kraus_clear(&composed);

        /* roll the accumulator forward */
        aic_ucp_kraus_clear(&acc);
        acc = compressed;                         /* shallow move (owns storage) */
    }

    /* hand the accumulator's storage to `out` (shallow move; acc is not cleared,
     * its buffers now belong to out, which the caller clears). */
    *out = acc;
}
