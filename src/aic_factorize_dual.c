/* aic_factorize_dual.c — the DUAL channels Dec = Delta*, Enc = Upsilon* of
 * th_factorization (bead aic-tff, module "factorize", Increment F4.1). Split out
 * of aic_factorize_verify.c to respect the <=200 LOC rule (Rule 10). See
 * include/aic_factorize.h (the F4 banner) and docs/research/factorize_f4_design.md
 * §B.3 (the four dual channels).
 *
 * approximate_algebras.tex:2159 — the dual factorization (verbatim, with our
 * Upsilon in the role of the paper's Gamma o Co_M):
 *   "Dec = Delta*,   Enc = (Gamma Co_M)*,   Dec Enc = (Upsilon Delta)* ~ 1_{B*},
 *    Enc Dec = (Delta Upsilon)* ~ Phi*."
 *
 * MIND THE DIRECTION (CLAUDE.md "UCP vs CPTP" callout; design §B.3). Delta: B ->
 * B(H) is the observable-ENCODE; its dual Dec = Delta*: B(H)* -> B* is the state-
 * DECODE. Upsilon: B(H) -> B is the observable-DECODE; its dual Enc = Upsilon*:
 * B* -> B(H)* is the state-ENCODE. The dual SWAPS the words encode/decode relative
 * to the observable maps; .tex:2159 (Dec = Delta*) is the authority. Getting this
 * backwards silently transposes everything.
 *
 * HOW THE DUAL IS READ OFF (constructive, no new construction). For a UCP map in
 * Heisenberg Kraus form Phi(X) = Sum_a K_a^dag X K_a, the dual is the Kraus-adjoint
 * Phi*(rho) = Sum_a K_a rho K_a^dag (the SAME K_a with the dagger moved). So we:
 *   (1) build the Convention-A Choi of Delta / Upsilon by the E_pq loop (the F2
 *       block-Choi pattern over the FULL ambient matrix-unit basis), then
 *   (2) extract Kraus {K_a} via aic_ucp_choi_to_kraus_latd (LAPACK double path).
 * The returned aic_ucp_kraus stores the SAME {K_a}; the dual's Schrodinger action
 * Phi*(rho) = Sum_a K_a rho K_a^dag is the caller's to apply. The verification
 * quantities ||Dec Enc - 1_{B*}||_⋄ = ||Upsilon Delta - 1_B||_cb and
 * ||Enc Dec - Phi*||_⋄ = ||Delta Upsilon - Phi||_cb are ALREADY computed by
 * aic_factorize_eigfree_{upsdel,delups} (the diamond-of-adjoint duality
 * ||Lambda*||_⋄ = ||Lambda||_cb, Watrous), so no separate diamond-SDP is needed.
 *
 * DUAL Kraus DIMS (DRIFT #2). Delta: dim_K = n_B (domain B), dim_H = N (codomain
 * B(H)); Choi(Delta) is (n_B*N)^2; choi_to_kraus_latd(&dec, ., n_B, N) gives Dec's
 * D_a: H -> B, each n_B x N. Upsilon: dim_K = N, dim_H = n_B; Choi(Upsilon) is
 * (N*n_B)^2; choi_to_kraus_latd(&enc, ., N, n_B) gives Enc's E_a: B -> H, each
 * N x n_B. Convention-A index order:
 *   Choi(Delta)[p*N + a, q*N + b]    = Delta(E_pq^{(n_B)})[a,b]   (p,q<n_B, a,b<N)
 *   Choi(Upsilon)[p*n_B + a, q*n_B + b] = Upsilon(E_pq^{(N)})[a,b] (p,q<N, a,b<n_B)
 *
 * The caller OWNS the returned aic_ucp_kraus and must aic_ucp_kraus_clear it.
 * Trace-preservation of Dec/Enc is FREE from the unitality of Delta/Upsilon
 * (Sum_a D_a^dag D_a = Delta(1_B) = 1_H; Sum_a E_a^dag E_a = Upsilon(1_H) = 1_B):
 * Dec is TP iff Delta is unital, Enc is TP iff Upsilon is unital (the dual TP-ness
 * tooth, T7 — a non-unital map would break it).
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_factorize.h"
#include "aic_ucp.h"

/* PSD-cone abort floor for the dual Kraus read-off (FINDINGS §C14). At multi-block
 * eta>0 Choi(Delta)/Choi(Upsilon) inherit the UCP maps' O(eta^2) cone defect (the
 * SAME root cause as the W_j builder), so the strict choi_to_kraus_latd would
 * abort. 1e-3 (absolute) admits the O(eta^2) defect (measured ~1e-6 at eta~1e-2,
 * ~400x below) but aborts a genuine O(1) — or unexpectedly-O(eta) — negative; it
 * matches AIC_FACTORIZE_CONE_NEG_TOL in the W_j builder. */
#define AIC_FACTORIZE_DUAL_NEG_TOL 1e-3

/* Build the Convention-A Choi of Delta (B -> B(H)) over the n_B^2 ambient M_{n_B}
 * units, then extract Dec = Delta* Kraus. `dec` is aic_ucp_kraus_init'd by
 * choi_to_kraus_latd (caller clears). REQUIRES F->delta_ready. */
void aic_factorize_dec_kraus(aic_ucp_kraus *dec, const aic_factorize *F, slong prec)
{
    assert(F != NULL && dec != NULL);
    assert(F->delta_ready && "dec_kraus: requires delta_ready (Rule 4)");
    slong N = F->N, nB = F->n_B, sz = nB * N;

    acb_mat_t C, Epq, DE;
    acb_mat_init(C, sz, sz);
    acb_mat_init(Epq, nB, nB);
    acb_mat_init(DE, N, N);
    acb_mat_zero(C);
    for (slong p = 0; p < nB; p++)
        for (slong q = 0; q < nB; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_delta(DE, F, Epq, prec);   /* Delta(E_pq), N x N */
            for (slong a = 0; a < N; a++)
                for (slong b = 0; b < N; b++)
                    acb_set(acb_mat_entry(C, p * N + a, q * N + b),
                            acb_mat_entry(DE, a, b));
        }
    acb_mat_clear(DE); acb_mat_clear(Epq);

    /* Dec = Delta*: dim_K = n_B (B), dim_H = N (H). TOLERANCE-AWARE: Choi(Delta)
     * inherits Delta's O(eta^2) cone defect (multi-block eta>0; FINDINGS §C14), so
     * the _tol extraction projects onto the PSD cone (neg_tol = 1e-3 absolute
     * floor, matching the W_j builder) and aborts only on a genuine non-CP input. */
    aic_ucp_choi_to_kraus_latd_tol(dec, C, nB, N, AIC_FACTORIZE_DUAL_NEG_TOL, NULL);
    acb_mat_clear(C);
}

/* Build the Convention-A Choi of Upsilon (B(H) -> B) over the N^2 ambient M_N
 * units, then extract Enc = Upsilon* Kraus. `enc` is aic_ucp_kraus_init'd by
 * choi_to_kraus_latd (caller clears). REQUIRES F->upsilon_ready. */
void aic_factorize_enc_kraus(aic_ucp_kraus *enc, const aic_factorize *F, slong prec)
{
    assert(F != NULL && enc != NULL);
    assert(F->upsilon_ready && "enc_kraus: requires upsilon_ready (Rule 4)");
    slong N = F->N, nB = F->n_B, sz = N * nB;

    acb_mat_t C, Epq, UE;
    acb_mat_init(C, sz, sz);
    acb_mat_init(Epq, N, N);
    acb_mat_init(UE, nB, nB);
    acb_mat_zero(C);
    for (slong p = 0; p < N; p++)
        for (slong q = 0; q < N; q++) {
            acb_mat_zero(Epq);
            acb_set_si(acb_mat_entry(Epq, p, q), 1);
            aic_factorize_upsilon(UE, F, Epq, prec); /* Upsilon(E_pq), n_B x n_B */
            for (slong a = 0; a < nB; a++)
                for (slong b = 0; b < nB; b++)
                    acb_set(acb_mat_entry(C, p * nB + a, q * nB + b),
                            acb_mat_entry(UE, a, b));
        }
    acb_mat_clear(UE); acb_mat_clear(Epq);

    /* Enc = Upsilon*: dim_K = N (H), dim_H = n_B (B). TOLERANCE-AWARE: Choi(Upsilon)
     * inherits the O(eta^2) cone defect at multi-block eta>0 (FINDINGS §C14); _tol
     * projects onto the PSD cone, aborting only on a genuine non-CP eigenvalue. */
    aic_ucp_choi_to_kraus_latd_tol(enc, C, N, nB, AIC_FACTORIZE_DUAL_NEG_TOL, NULL);
    acb_mat_clear(C);
}
