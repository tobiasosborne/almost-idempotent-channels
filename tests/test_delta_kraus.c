/* test_delta_kraus.c — cross-checks for the prop_Delta explicit Kraus form of the
 * inclusion map Delta : A -> B(H) (bead aic-ynu, I4; the eta=0 case). Realizes
 * Delta_structure (approximate_algebras.tex:2098-2103):
 *   Delta(A) = sum_j J_M W_j^dag (A_j (x) 1_{E_j}) W_j J_M^dag   (M-part = w(A))
 *            + J_{Mperp} Sigma(A) J_{Mperp}^dag                  (M^perp-part)
 * The M-part Kraus form is sum_{j,c} K_{j,c} A_hat K_{j,c}^dag,
 *   K_{j,c} = J_M W_j^dag (1_{L_j} (x) e_c) P_j  (an n x dL_tot operator).
 * aic_idemp_delta_kraus builds the K_{j,c} + the M^perp Sigma matrix and certifies
 * (in production) that the full Delta-Kraus map MATCHES the stored d->Delta. This
 * test re-asserts the gates INDEPENDENTLY and reports the numbers.
 *
 * THE CRUX TOOTH (.tex:2100). The Delta-Kraus map must EQUAL d->Delta: for each
 * A-basis element B_k, reconstruct Delta_kraus(B_k) (M-part Kraus + M^perp Sigma)
 * and compare to d->Delta column k reshaped, ||Delta_kraus(B_k) - B_k|| <= 1e-10.
 * The MUTATION PROOF (drop a Kraus operator) confirms the tooth bites.
 *
 * Plus: (a) UCP — M-part unital onto M (Delta_M(1_A) = Pi_M); (b) round-trip
 * Delta.Gamma.C_M = Phi (.tex:2080): apply stored Gamma to C_M(E_ij), feed the
 * A-coordinates to Delta_kraus, recover Phi(E_ij).
 *
 * SCOPE: noiseless_subsystem(2,2) [E_1 = C^2, nontrivial multiplicity, M=H],
 * block_cond_exp(5,2) [M_2 (+) M_3, E_j = C^1, M=H], identity(3) [M=H],
 * block_cond_exp(4,2) [M_2 (+) M_2, E_j = C^1, M=H]. (All these channels are
 * trace-preserving so M=H and Sigma is empty; the Sigma read-off path + dim_Mperp>0
 * is exercised by make_compress_idemp(4,2) [PROPER carrier, dim_Mperp=2].)
 */
#include <math.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_idemp.h"
#include "aic/aic_mat.h"
#include "aic/aic_ucp.h"
#include "aic_test.h"
#include "test_idemp.h"

#define DKPREC 128
static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* reshape column k of an (n^2 x cols) row-major-vec matrix into n x n. */
static void col_to_mat(acb_mat_t B, const acb_mat_t M, slong k, slong n)
{
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++)
            acb_set(acb_mat_entry(B, a, b), acb_mat_entry(M, a * n + b, k));
}

/* Apply the Delta-Kraus map to an A-element given by A-coordinates a (dim_A x 1):
 * wA = reshape(d->w * a) = w(A) in B(M); A_hat = (+)_j Tr_{E_j}(W_j wA W_j^dag)/dE;
 * M-part = sum_{j,c} K_{j,c} A_hat K_{j,c}^dag; M^perp-part = J_Mperp Sigma(a)
 * J_Mperp^dag. Returns the n x n operator on H. Independent re-impl of the
 * production reconstruction (uses ONLY the public aic_idemp_delta + W_j). */
static void apply_delta(acb_mat_t out, const aic_idemp_delta *D,
                        const aic_idemp_wedderburn *W, const aic_idemp_decomp *d,
                        const acb_mat_t a, slong prec)
{
    slong n = D->n, m = D->dim_M, dmp = D->dim_Mperp;

    /* wA = w(A) in B(M) */
    acb_mat_t wv, wA;
    acb_mat_init(wv, m * m, 1);
    acb_mat_mul(wv, d->w, a, prec);
    acb_mat_init(wA, m, m);
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            acb_set(acb_mat_entry(wA, i, j), acb_mat_entry(wv, i * m + j, 0));

    /* A_hat = (+)_j Tr_{E_j}(W_j wA W_j^dag)/dE at block offset off_L[j] */
    acb_mat_t Ahat;
    acb_mat_init(Ahat, D->dL_tot, D->dL_tot);
    acb_mat_zero(Ahat);
    for (slong j = 0; j < W->num_blocks; j++) {
        slong dL = W->dim_L[j], dE = W->dim_E[j], db = dL * dE, off = D->off_L[j];
        acb_mat_t Wjd, WX, WXW, Aj;
        acb_mat_init(Wjd, m, db);
        acb_mat_conjugate_transpose(Wjd, W->W_j[j]);
        acb_mat_init(WX, db, m);
        acb_mat_mul(WX, W->W_j[j], wA, prec);
        acb_mat_init(WXW, db, db);
        acb_mat_mul(WXW, WX, Wjd, prec);
        acb_mat_init(Aj, dL, dL);
        aic_mat_partial_trace_right(Aj, WXW, dL, dE, prec);
        for (slong p = 0; p < dL; p++)
            for (slong q = 0; q < dL; q++) {
                acb_div_si(acb_mat_entry(Aj, p, q), acb_mat_entry(Aj, p, q), dE, prec);
                acb_set(acb_mat_entry(Ahat, off + p, off + q),
                        acb_mat_entry(Aj, p, q));
            }
        acb_mat_clear(Aj);
        acb_mat_clear(WXW);
        acb_mat_clear(WX);
        acb_mat_clear(Wjd);
    }

    /* M-part = sum_i K_i A_hat K_i^dag */
    acb_mat_zero(out);
    acb_mat_t Kd, KA, KAK;
    acb_mat_init(Kd, D->dL_tot, n);
    acb_mat_init(KA, n, D->dL_tot);
    acb_mat_init(KAK, n, n);
    for (slong i = 0; i < D->n_kraus; i++) {
        acb_mat_mul(KA, D->K[i], Ahat, prec);
        acb_mat_conjugate_transpose(Kd, D->K[i]);
        acb_mat_mul(KAK, KA, Kd, prec);
        acb_mat_add(out, out, KAK, prec);
    }
    acb_mat_clear(KAK);
    acb_mat_clear(KA);
    acb_mat_clear(Kd);

    /* M^perp-part = J_Mperp Sigma(a) J_Mperp^dag */
    if (dmp > 0) {
        acb_mat_t sv, S, JpS, Jpd, JpSJ;
        acb_mat_init(sv, dmp * dmp, 1);
        acb_mat_mul(sv, D->Sigma, a, prec);
        acb_mat_init(S, dmp, dmp);
        for (slong p = 0; p < dmp; p++)
            for (slong q = 0; q < dmp; q++)
                acb_set(acb_mat_entry(S, p, q), acb_mat_entry(sv, p * dmp + q, 0));
        acb_mat_init(JpS, n, dmp);
        acb_mat_mul(JpS, d->J_Mperp, S, prec);
        acb_mat_init(Jpd, dmp, n);
        acb_mat_conjugate_transpose(Jpd, d->J_Mperp);
        acb_mat_init(JpSJ, n, n);
        acb_mat_mul(JpSJ, JpS, Jpd, prec);
        acb_mat_add(out, out, JpSJ, prec);
        acb_mat_clear(JpSJ);
        acb_mat_clear(Jpd);
        acb_mat_clear(JpS);
        acb_mat_clear(S);
        acb_mat_clear(sv);
    }
    acb_mat_clear(Ahat);
    acb_mat_clear(wA);
    acb_mat_clear(wv);
}

/* max_k ||Delta_kraus(B_k) - B_k||_op (the crux: B_k = d->Delta col k, a = e_k). */
static double match_resid(const aic_idemp_delta *D, const aic_idemp_wedderburn *W,
                          const aic_idemp_decomp *d, slong prec)
{
    slong n = d->n, dA = d->dim_A;
    double worst = 0;
    acb_mat_t a, rec, Bk, diff;
    arb_t e;
    acb_mat_init(a, dA, 1);
    acb_mat_init(rec, n, n);
    acb_mat_init(Bk, n, n);
    acb_mat_init(diff, n, n);
    arb_init(e);
    for (slong k = 0; k < dA; k++) {
        acb_mat_zero(a);
        acb_set_si(acb_mat_entry(a, k, 0), 1);
        apply_delta(rec, D, W, d, a, prec);
        col_to_mat(Bk, d->Delta, k, n);
        acb_mat_sub(diff, rec, Bk, prec);
        aic_mat_opnorm(e, diff, prec);
        if (dd(e) > worst) worst = dd(e);
    }
    arb_clear(e);
    acb_mat_clear(diff);
    acb_mat_clear(Bk);
    acb_mat_clear(rec);
    acb_mat_clear(a);
    return worst;
}

/* M-part unitality: ||Delta_M(1_A) - Pi_M||_op (1_A = abstract unit, A-coords for
 * the unit are NOT e-vectors; instead feed A_hat = 1_{dL_tot} directly). */
static double unital_defect(const aic_idemp_delta *D, const aic_idemp_decomp *d,
                            slong prec)
{
    slong n = D->n;
    acb_mat_t one, Mpart, Kd, KA, KAK, diff;
    arb_t e;
    acb_mat_init(one, D->dL_tot, D->dL_tot);
    acb_mat_one(one);
    acb_mat_init(Mpart, n, n);
    acb_mat_zero(Mpart);
    acb_mat_init(Kd, D->dL_tot, n);
    acb_mat_init(KA, n, D->dL_tot);
    acb_mat_init(KAK, n, n);
    for (slong i = 0; i < D->n_kraus; i++) {
        acb_mat_mul(KA, D->K[i], one, prec);
        acb_mat_conjugate_transpose(Kd, D->K[i]);
        acb_mat_mul(KAK, KA, Kd, prec);
        acb_mat_add(Mpart, Mpart, KAK, prec);
    }
    acb_mat_init(diff, n, n);
    acb_mat_sub(diff, Mpart, d->Pi_M, prec);
    arb_init(e);
    aic_mat_opnorm(e, diff, prec);
    double r = dd(e);
    arb_clear(e);
    acb_mat_clear(diff);
    acb_mat_clear(KAK);
    acb_mat_clear(KA);
    acb_mat_clear(Kd);
    acb_mat_clear(Mpart);
    acb_mat_clear(one);
    return r;
}

/* round-trip Delta.Gamma.C_M = Phi (.tex:2080): for each E_ij, a = Gamma vec(C_M(E_ij)),
 * compare Delta_kraus(a) to Phi(E_ij). Returns max_ij ||...||_op. */
static double roundtrip_defect(const aic_idemp_delta *D, const aic_idemp_wedderburn *W,
                               const aic_idemp_decomp *d, const aic_ucp_kraus *phi,
                               slong prec)
{
    slong n = d->n, m = d->dim_M, dA = d->dim_A;
    double worst = 0;
    acb_mat_t E, JMd, CME, cmv, a, rec, PhiE, diff;
    arb_t e;
    acb_mat_init(JMd, m, n);
    acb_mat_conjugate_transpose(JMd, d->J_M);
    acb_mat_init(E, n, n);
    acb_mat_init(CME, m, m);
    acb_mat_init(cmv, m * m, 1);
    acb_mat_init(a, dA, 1);
    acb_mat_init(rec, n, n);
    acb_mat_init(PhiE, n, n);
    acb_mat_init(diff, n, n);
    arb_init(e);
    acb_mat_t JdE;
    acb_mat_init(JdE, m, n);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            acb_mat_zero(E);
            acb_set_si(acb_mat_entry(E, i, j), 1);
            acb_mat_mul(JdE, JMd, E, prec);
            acb_mat_mul(CME, JdE, d->J_M, prec);       /* C_M(E_ij) = J_M^dag E J_M */
            for (slong p = 0; p < m; p++)
                for (slong q = 0; q < m; q++)
                    acb_set(acb_mat_entry(cmv, p * m + q, 0), acb_mat_entry(CME, p, q));
            acb_mat_mul(a, d->Gamma, cmv, prec);       /* a = Gamma(C_M(E_ij)) */
            apply_delta(rec, D, W, d, a, prec);        /* Delta_kraus(a) */
            aic_ucp_apply(PhiE, phi, E, prec);         /* Phi(E_ij) */
            acb_mat_sub(diff, rec, PhiE, prec);
            aic_mat_opnorm(e, diff, prec);
            if (dd(e) > worst) worst = dd(e);
        }
    acb_mat_clear(JdE);
    arb_clear(e);
    acb_mat_clear(diff);
    acb_mat_clear(PhiE);
    acb_mat_clear(rec);
    acb_mat_clear(a);
    acb_mat_clear(cmv);
    acb_mat_clear(CME);
    acb_mat_clear(E);
    acb_mat_clear(JMd);
    return worst;
}

static void check_delta(const char *name, aic_ucp_kraus *phi)
{
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, phi, DKPREC);
    aic_idemp_wedderburn W;
    aic_idemp_wedderburn_decompose(&W, &d, DKPREC);
    aic_idemp_delta D;
    aic_idemp_delta_kraus(&D, &W, &d, DKPREC);   /* production self-certifies match */

    double mr = match_resid(&D, &W, &d, DKPREC);
    AIC_CHECK_MSG(mr <= 1e-10, "%s: match-d.Delta residual %.3e > 1e-10", name, mr);

    double ud = unital_defect(&D, &d, DKPREC);
    AIC_CHECK_MSG(ud <= 1e-10, "%s: Delta_M not unital onto M (%.3e)", name, ud);

    double rt = roundtrip_defect(&D, &W, &d, phi, DKPREC);
    AIC_CHECK_MSG(rt <= 1e-10, "%s: Delta.Gamma.C_M != Phi (%.3e)", name, rt);

    printf("  %-28s blocks=%ld n_kraus=%ld dimMperp=%ld  match=%.2e  unital=%.2e  "
           "Delta.Gamma=%.2e\n", name, (long) D.num_blocks, (long) D.n_kraus,
           (long) D.dim_Mperp, mr, ud, rt);

    aic_idemp_delta_clear(&D);
    aic_idemp_wedderburn_clear(&W);
    aic_idemp_clear(&d);
}

/* MUTATION PROOF (Rule 7): drop one M-part Kraus operator (zero it) on
 * noiseless_subsystem(2,2) and confirm the match-d.Delta tooth goes RED — a missing
 * K_{j,c} removes one c-slice of the multiplicity sum, so the M-part of Delta(B_k)
 * is wrong. Restore (re-decompose) -> GREEN. */
static void mutation_proof(void)
{
    printf("MUTATION PROOF (noiseless_subsystem(2,2), drop one M-part Kraus op):\n");
    aic_ucp_kraus phi;
    make_noiseless_subsystem(&phi, 2, 2);
    aic_idemp_decomp d;
    aic_idemp_decompose(&d, &phi, DKPREC);
    aic_idemp_wedderburn W;
    aic_idemp_wedderburn_decompose(&W, &d, DKPREC);
    aic_idemp_delta D;
    aic_idemp_delta_kraus(&D, &W, &d, DKPREC);

    double base = match_resid(&D, &W, &d, DKPREC);
    AIC_CHECK_MSG(base <= 1e-10, "mutation: baseline match %.3e not green", base);

    /* drop K[n_kraus-1] (zero it): one multiplicity slice missing -> M-part wrong. */
    slong last = D.n_kraus - 1;
    acb_mat_t saved;
    acb_mat_init(saved, acb_mat_nrows(D.K[last]), acb_mat_ncols(D.K[last]));
    acb_mat_set(saved, D.K[last]);
    acb_mat_zero(D.K[last]);
    double red = match_resid(&D, &W, &d, DKPREC);
    printf("  dropped K[%ld]:                  match-d.Delta = %.3e (RED line)\n",
           (long) last, red);
    AIC_CHECK_MSG(red > 1e-3, "mutation: dropped-Kraus match %.3e did NOT go red", red);

    acb_mat_set(D.K[last], saved);
    double restored = match_resid(&D, &W, &d, DKPREC);
    printf("  restored K[%ld]:                 match-d.Delta = %.3e (back to green)\n",
           (long) last, restored);
    AIC_CHECK_MSG(restored <= 1e-10, "mutation: restored match %.3e not green", restored);

    acb_mat_clear(saved);
    aic_idemp_delta_clear(&D);
    aic_idemp_wedderburn_clear(&W);
    aic_idemp_clear(&d);
    aic_ucp_kraus_clear(&phi);
}

static void test_oracles(void)
{
    printf("prop_Delta Kraus form (eta=0 oracle, match-d.Delta tooth):\n");
    aic_ucp_kraus phi;

    /* noiseless_subsystem(2,2): E_1 = C^2, nontrivial multiplicity, M=H. */
    make_noiseless_subsystem(&phi, 2, 2);
    check_delta("noiseless_subsystem(2,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* block_cond_exp(5,2): M_2 (+) M_3, E_j = C^1, M=H. */
    make_block_cond_exp(&phi, 5, 2);
    check_delta("block_cond_exp(5,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* identity(3): single trivial block, M=H. */
    make_identity(&phi, 3);
    check_delta("identity(3)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* block_cond_exp(4,2): M_2 (+) M_2, E_j = C^1, M=H. */
    make_block_cond_exp(&phi, 4, 2);
    check_delta("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* compress_idemp(4,2): PROPER carrier (dim_M=2 < n=4), dim_Mperp=2 -> exercises
     * the Sigma / M^perp read-off path (the others are all trace-preserving, M=H). */
    make_compress_idemp(&phi, 4, 2);
    check_delta("compress_idemp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_oracles();
    mutation_proof();
    aic_test_report("test_delta_kraus");
    printf("OK test_delta_kraus\n");
    return 0;
}
