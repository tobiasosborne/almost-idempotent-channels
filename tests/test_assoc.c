/* test_assoc.c — cross-checks for assoc_ecsa Increment 1 (bead aic-92f): the
 * superoperator build S_Phi and the regularization Phi_tilde = theta(2 Phi - 1).
 * Realizes/anchors approximate_algebras.tex:2171-2182 (tilde_Phi + its
 * properties) and :524-533 (prop_P) at the superoperator level.
 *
 * Cross-check ladder (funcalc's arb-only precedent: a prec=53-vs-256 rung
 * substitutes for double-vs-arb where there is no double sgn).
 *   T1 SUPEROPERATOR ORACLE  : S_Phi vec_r(X) == vec_r(Phi(X)), machine-zero,
 *                              on >=2 channels incl. an asymmetric/complex one.
 *   T2 KRONECKER cross-check : column-image S_Phi == sum_a K_a^dag (x) K_a^T.
 *   T3 eta=0 ORACLE (headline): S_Phitilde == S_Phi machine-zero on each exact
 *                              idempotent channel (reuses test_idemp.h).
 *   T4 IDEMPOTENCY           : S_Phitilde^2 == S_Phitilde machine-zero.
 *   T5 UNIT + HP             : Phi_tilde(I)=I; Phi_tilde(X^dag)=Phi_tilde(X)^dag.
 *   T6 NON-NORMAL sgn check  : route A1 (prop_P) vs A2 (binomial xpow) agree
 *                              ~1e-25 @256 on TWO structurally-different
 *                              genuinely-NON-NORMAL in-basin channels.
 *   T7 ||Phi_tilde-Phi||=O(eta): (a) linear-in-defect growth + ->0 at t->0 on a
 *                              smooth family; (b) ratio dist/eta stays BOUNDED
 *                              across a DIMENSION-NONTRIVIAL family (d=3,4,6,8).
 *                              PRELIMINARY -- the definitive axiom-defect
 *                              universality canary is Increment 2 (not this).
 *   T8 PRECISION LADDER      : Phi_tilde @53 vs @256 agree ~1e-12.
 *   T9 FAIL-LOUD basin       : documented (aic_prop_P aborts out-of-basin).
 *
 * MUTATION-PROVEN teeth (see the comments at T1 and T4): a transposed reshape /
 * swapped Kronecker factor breaks T1; returning S_Phi unmodified for a
 * non-idempotent channel breaks T4. Verified RED, then restored.
 */
#include <math.h>
#include <stdio.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
#include "aic_funcalc.h"
#include "aic_mat.h"
#include "aic_test.h"
#include "aic_ucp.h"

#include "test_idemp.h"   /* exact-idempotent channel builders + make_conjugated */

/* ---- helpers -------------------------------------------------------------- */

/* A fixed, reproducible complex non-Hermitian n x n matrix (deterministic; no
 * RNG, project rule). Entry (i,j) = (i+1)/(n+i+j+1) + i*((j+1)/(n+2i+j+2)). */
static void fill_fixed_complex(acb_mat_t X, slong n)
{
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            double re = (double) (i + 1) / (double) (n + i + j + 1);
            double im = (double) (j + 1) / (double) (n + 2 * i + j + 2);
            acb_set_d_d(acb_mat_entry(X, i, j), re, im);
        }
}

/* Operator norm of the MIDPOINT of A (radii stripped). The funcalc iterations
 * (Newton-Schulz / xpow) inflate A's error balls to ~5e-75 at prec=256; squaring
 * those into the Gram A^dag A propagates a radius that trips aic_mat_opnorm's
 * strict Hermiticity assert (tol 2^-(prec-8) ~ 2.2e-75) even when the midpoint is
 * exactly Hermitian. The TREND measurements (T7, defect proxy) are honest
 * midpoint magnitudes, not certified balls (rigor lives in the T3 eta=0 oracle
 * and the T6 A1-vs-A2 cross-check), so we strip radii first. */
static double midpoint_opnorm(const acb_mat_t A, slong prec)
{
    slong m = acb_mat_nrows(A), c = acb_mat_ncols(A);
    acb_mat_t M;
    arb_t e;
    double out;
    acb_mat_init(M, m, c);
    arb_init(e);
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(A, i, j));
    aic_mat_opnorm(e, M, prec);
    out = arf_get_d(arb_midref(e), ARF_RND_NEAR);
    arb_clear(e);
    acb_mat_clear(M);
    return out;
}

/* ||[S, S^dag]||_F as a double (non-normality witness; 0 iff S normal). */
static double commutator_normF(const acb_mat_t S, slong prec)
{
    slong m = acb_mat_nrows(S);
    acb_mat_t Sd, A, B, C;
    arb_t f;
    double out;
    acb_mat_init(Sd, m, m);
    acb_mat_init(A, m, m);
    acb_mat_init(B, m, m);
    acb_mat_init(C, m, m);
    arb_init(f);
    acb_mat_conjugate_transpose(Sd, S);
    acb_mat_mul(A, S, Sd, prec);     /* S S^dag */
    acb_mat_mul(B, Sd, S, prec);     /* S^dag S */
    acb_mat_sub(C, A, B, prec);
    aic_mat_frobenius_norm(f, C, prec);
    out = arf_get_d(arb_midref(f), ARF_RND_NEAR);
    arb_clear(f);
    acb_mat_clear(C); acb_mat_clear(B); acb_mat_clear(A); acb_mat_clear(Sd);
    return out;
}

/* op-norm of S^2 - S (the eta proxy used by T7) as a double. */
static double defect_opnorm(const acb_mat_t S, slong prec)
{
    slong m = acb_mat_nrows(S);
    acb_mat_t S2, D;
    double out;
    acb_mat_init(S2, m, m);
    acb_mat_init(D, m, m);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    out = midpoint_opnorm(D, prec);
    acb_mat_clear(D); acb_mat_clear(S2);
    return out;
}

/* Almost-idempotent UCP family: Phi_t = (1-t) base + t mix, as the Heisenberg
 * Kraus union { sqrt(1-t) Kbase_a } U { sqrt(t) Kmix_b }. A convex combination
 * of UCP maps is UCP (unital: (1-t) 1 + t 1 = 1). At t=0 it is `base` (exactly
 * idempotent if base is). `out` is init'd here. */
static void make_mix(aic_ucp_kraus *out, const aic_ucp_kraus *base,
                     const aic_ucp_kraus *mix, double t, slong prec)
{
    slong n = base->dim_H;
    AIC_CHECK_MSG(base->dim_K == n && mix->dim_K == n && mix->dim_H == n,
                  "make_mix: base and mix must be self-maps on the same C^n");
    aic_ucp_kraus_init(out, n, n, base->r + mix->r);
    arb_t s;
    arb_init(s);
    arb_set_d(s, sqrt(1.0 - t));
    for (slong a = 0; a < base->r; a++)
        acb_mat_scalar_mul_arb(out->K[a], base->K[a], s, prec);
    arb_set_d(s, sqrt(t));
    for (slong b = 0; b < mix->r; b++)
        acb_mat_scalar_mul_arb(out->K[base->r + b], mix->K[b], s, prec);
    arb_clear(s);
}

/* Route A2 (independent binomial-series regularization, test-only): the third
 * equality of tilde_Phi (.tex:2174):
 *   D = 4(S^2 - S); M = I - D; Minvsqrt = xpow(M, -1/2, 1);
 *   S_tilde = (1/2)(I + (2 S - I) Minvsqrt). */
static void regularize_A2(acb_mat_t S_tilde, const acb_mat_t S, slong prec)
{
    slong m = acb_mat_nrows(S);
    acb_mat_t S2, D, M, Minv, X2mI, prod, I;
    acb_mat_init(S2, m, m); acb_mat_init(D, m, m); acb_mat_init(M, m, m);
    acb_mat_init(Minv, m, m); acb_mat_init(X2mI, m, m); acb_mat_init(prod, m, m);
    acb_mat_init(I, m, m);
    acb_mat_one(I);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S, S2, prec);                   /* S - S^2 = Phi - Phi^2 */
    acb_mat_scalar_mul_si(D, D, 4, prec);          /* D = 4(Phi - Phi^2) */
    acb_mat_sub(M, I, D, prec);                    /* M = 1 - 4(Phi - Phi^2) = (2S-I)^2 */
    aic_funcalc_xpow(Minv, M, -0.5, 1.0, prec);    /* (1 - D)^{-1/2} = (X^2)^{-1/2} */
    acb_mat_scalar_mul_si(X2mI, S, 2, prec);
    acb_mat_sub(X2mI, X2mI, I, prec);              /* 2 S - I */
    acb_mat_mul(prod, X2mI, Minv, prec);           /* (2S - I)(1-D)^{-1/2} = sgn */
    acb_mat_add(S_tilde, I, prod, prec);
    acb_mat_scalar_mul_2exp_si(S_tilde, S_tilde, -1); /* (1/2)(I + sgn) */
    acb_mat_clear(I); acb_mat_clear(prod); acb_mat_clear(X2mI);
    acb_mat_clear(Minv); acb_mat_clear(M); acb_mat_clear(D); acb_mat_clear(S2);
}

/* ---- T6 helper: NON-NORMAL sgn cross-check on one channel `phi`. ------------
 * Builds S_Phi, asserts the channel is genuinely NON-NORMAL (||[S,S^dag]||_F >
 * 1e-3) and IN-BASIN (4||S^2-S||_op < 1, both printed), then asserts route A1
 * (aic_assoc_regularize = prop_P, Newton-Schulz sgn) equals route A2
 * (regularize_A2 = binomial xpow) to ~1e-25 at prec=256. A1 and A2 are
 * algorithmically independent, so agreement on a non-normal S is both an
 * independent check of the regularization AND coverage for funcalc's
 * never-before-tested non-normal sgn path (.tex:2174 third equality). */
static void T6_nonnormal(const char *name, const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H;
    acb_mat_t S, A1, A2, diff;
    arb_t e, tol;
    acb_mat_init(S, n * n, n * n);
    acb_mat_init(A1, n * n, n * n);
    acb_mat_init(A2, n * n, n * n);
    acb_mat_init(diff, n * n, n * n);
    arb_init(e);
    arb_init(tol);
    arb_set_d(tol, 1e-25);
    aic_assoc_superop_from_ucp(S, phi, prec);

    double cnorm = commutator_normF(S, prec);
    double eta_proxy = defect_opnorm(S, prec);
    AIC_CHECK_MSG(cnorm > 1e-3,
                  "%s: channel must be genuinely NON-NORMAL (||[S,S^dag]||_F=%.3e)",
                  name, cnorm);
    AIC_CHECK_MSG(4.0 * eta_proxy < 1.0,
                  "%s: channel must be IN-BASIN (4||S^2-S||_op=%.3f >= 1)",
                  name, 4.0 * eta_proxy);

    aic_assoc_regularize(A1, phi, prec);   /* route A1: prop_P (Newton-Schulz) */
    regularize_A2(A2, S, prec);            /* route A2: binomial xpow          */
    acb_mat_sub(diff, A1, A2, prec);
    aic_mat_frobenius_norm(e, diff, prec);
    AIC_CHECK_MSG(arb_le(e, tol),
                  "%s: A1 (prop_P) vs A2 (xpow) disagree on non-normal channel "
                  "||A1-A2||_F=%.3e", name, arf_get_d(arb_midref(e), ARF_RND_NEAR));
    printf("  T6 [%s] A1-vs-A2: %.3e; ||[S,S^dag]||_F=%.3f; 4||S^2-S||_op=%.3f\n",
           name, arf_get_d(arb_midref(e), ARF_RND_NEAR), cnorm, 4.0 * eta_proxy);

    arb_clear(tol);
    arb_clear(e);
    acb_mat_clear(diff);
    acb_mat_clear(A2);
    acb_mat_clear(A1);
    acb_mat_clear(S);
}

/* ---- T1: superoperator oracle (mutation-proven) -------------------------- */
/* S_Phi vec_r(X) == vec_r(Phi(X)). MUTATION: a transposed reshape (storing
 * Phi(E_pq)[b,a] in column (p*n+q), or applying with a swapped index in
 * aic_assoc_superop_apply) makes Phi(X)_via_S != aic_ucp_apply for the
 * asymmetric complex channel -> RED. Verified, restored. */
static void T1_oracle(const char *name, const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H;
    acb_mat_t S, X, viaS, viaPhi;
    arb_t tol;
    acb_mat_init(S, n * n, n * n);
    acb_mat_init(X, n, n);
    acb_mat_init(viaS, n, n);
    acb_mat_init(viaPhi, n, n);
    arb_init(tol);
    arb_set_d(tol, (prec >= 200) ? 1e-25 : 1e-12);
    aic_assoc_superop_from_ucp(S, phi, prec);
    for (int trial = 0; trial < 3; trial++) {
        /* three distinct fixed X: complex non-Hermitian, then two shifts. */
        fill_fixed_complex(X, n);
        if (trial == 1)
            for (slong i = 0; i < n; i++)
                acb_add_si(acb_mat_entry(X, i, i), acb_mat_entry(X, i, i), 2, prec);
        if (trial == 2)  /* a Hermitian X: X <- X + X^dag */
        {
            acb_mat_t Xd;
            acb_mat_init(Xd, n, n);
            acb_mat_conjugate_transpose(Xd, X);
            acb_mat_add(X, X, Xd, prec);
            acb_mat_clear(Xd);
        }
        aic_assoc_superop_apply(viaS, S, X, prec);
        aic_ucp_apply(viaPhi, phi, X, prec);
        AIC_CHECK_ACB_MAT_CLOSE(viaS, viaPhi, tol);
    }
    /* T2: Kronecker cross-check (same convention, second way). */
    {
        acb_mat_t Sk;
        acb_mat_init(Sk, n * n, n * n);
        aic_assoc_superop_kron(Sk, phi, prec);
        AIC_CHECK_ACB_MAT_CLOSE(S, Sk, tol);
        acb_mat_clear(Sk);
    }
    (void) name;
    arb_clear(tol);
    acb_mat_clear(viaPhi); acb_mat_clear(viaS); acb_mat_clear(X); acb_mat_clear(S);
}

/* ---- T3/T4/T5: eta=0 oracle, idempotency, unit + HP ---------------------- */
static void T345_idemp_channel(const char *name, aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H;
    acb_mat_t S, St, S2, diff;
    arb_t tol, tol_def, e;
    acb_mat_init(S, n * n, n * n);
    acb_mat_init(St, n * n, n * n);
    acb_mat_init(S2, n * n, n * n);
    acb_mat_init(diff, n * n, n * n);
    arb_init(tol);
    arb_init(tol_def);
    arb_init(e);
    /* Machine-floor tol (for properties prop_P/theta make EXACT regardless of
     * the input: idempotency, Hermicity). prec=256 leaves ~1e-22 headroom. */
    arb_set_d(tol, 1e-22);

    aic_assoc_superop_from_ucp(S, phi, prec);
    aic_assoc_regularize(St, phi, prec);

    /* defect-aware tol for the eta=0 oracle. Some "exact" idempotents carry
     * IRRATIONAL Kraus entries (trace_replace, noiseless_subsystem use 1/sqrt k
     * set from a double), so Phi is idempotent only to the input's own residual
     * delta = ||S^2 - S||_op (~2.7e-16 for trace_replace(3)); integer-Kraus
     * channels are exact (delta ~ 1e-30). prop_P bounds ||S_tilde - S|| <=
     * ||2S-I|| O(delta), so the HONEST oracle bound is C*delta + machine-floor.
     * This is NOT a bandaid: it is the prop_P bound with the input's true delta
     * (the channel reduces to its OWN idempotence precision). */
    {
        acb_mat_t D;
        arb_t delta;
        arf_t ub;
        acb_mat_init(D, n * n, n * n);
        arb_init(delta);
        arf_init(ub);
        acb_mat_sqr(S2, S, prec);
        acb_mat_sub(D, S2, S, prec);
        /* Frobenius defect: exact-radius and a rigorous upper bound on the
         * idempotence defect (||.||_op <= ||.||_F). Take its UPPER bound as a
         * SCALAR (not a ball) so tol_def is a thin threshold, not a wide ball
         * whose arb_le would compare against its near-zero lower end. */
        aic_mat_frobenius_norm(delta, D, prec);
        arb_get_ubound_arf(ub, delta, prec);
        arb_set_arf(tol_def, ub);
        /* x100: a CHOSEN safety multiplier on top of the prop_P bound
         * ||S_tilde - S|| <= ||2S-I|| * O(delta) (which is tight to ~1x in
         * practice -- see T3 prints, ~3.5e-16 vs delta ~2.7e-16), NOT a derived
         * constant. It is headroom against the O(.) hidden factor + accumulated
         * funcalc-iteration rounding, not a claim that prop_P's constant is 100. */
        arb_mul_si(tol_def, tol_def, 100, prec);
        arb_add(tol_def, tol_def, tol, prec);     /* + machine floor */
        arf_clear(ub);
        arb_clear(delta);
        acb_mat_clear(D);
    }

    /* T3: eta=0 oracle -- S_tilde == S_Phi to the input's idempotence precision.
     * Since Phi^2=Phi (up to delta), (2S-1)^2 = I (up to 4delta), so theta gives
     * S back up to O(delta). For integer-Kraus channels delta~1e-30 => machine-
     * zero; for 1/sqrt-k channels => ~1e-16. */
    acb_mat_sub(diff, St, S, prec);
    aic_mat_frobenius_norm(e, diff, prec);
    AIC_CHECK_MSG(arb_le(e, tol_def),
                  "%s: T3 eta=0 oracle ||S_tilde - S_Phi||_F too big", name);
    double t3 = arf_get_d(arb_midref(e), ARF_RND_NEAR);

    /* T4: S_tilde idempotent EXACTLY (machine-floor tol, not defect-aware):
     * prop_P/sgn give sgn^2 = I exactly, so theta^2 = theta exactly, INDEPENDENT
     * of the input defect. MUTATION: aic_assoc_regularize returning S unmodified
     * passes T4 here only because an exact idempotent's S is itself idempotent;
     * the unmodified-S mutation is caught on the NON-idempotent T4 teeth below. */
    acb_mat_sqr(S2, St, prec);
    acb_mat_sub(diff, S2, St, prec);
    aic_mat_frobenius_norm(e, diff, prec);
    AIC_CHECK_MSG(arb_le(e, tol), "%s: T4 S_tilde not idempotent", name);

    /* T5a: Phi_tilde(I) = I. I is a fixed point of Phi (up to delta), so this
     * holds to the input's idempotence precision -> defect-aware tol. */
    {
        acb_mat_t In, out, d2;
        acb_mat_init(In, n, n);
        acb_mat_init(out, n, n);
        acb_mat_init(d2, n, n);
        acb_mat_one(In);
        aic_assoc_superop_apply(out, St, In, prec);
        acb_mat_sub(d2, out, In, prec);
        aic_mat_frobenius_norm(e, d2, prec);
        AIC_CHECK_MSG(arb_le(e, tol_def), "%s: T5 Phi_tilde(I) != I", name);
        acb_mat_clear(d2); acb_mat_clear(out); acb_mat_clear(In);
    }
    /* T5b: Phi_tilde(X^dag) = Phi_tilde(X)^dag (Hermicity preservation). */
    {
        acb_mat_t X, Xd, lhs, rhs, d2;
        acb_mat_init(X, n, n); acb_mat_init(Xd, n, n);
        acb_mat_init(lhs, n, n); acb_mat_init(rhs, n, n); acb_mat_init(d2, n, n);
        fill_fixed_complex(X, n);
        acb_mat_conjugate_transpose(Xd, X);
        aic_assoc_superop_apply(lhs, St, Xd, prec);        /* Phi_tilde(X^dag) */
        aic_assoc_superop_apply(rhs, St, X, prec);
        acb_mat_conjugate_transpose(rhs, rhs);             /* Phi_tilde(X)^dag */
        acb_mat_sub(d2, lhs, rhs, prec);
        aic_mat_frobenius_norm(e, d2, prec);
        AIC_CHECK_MSG(arb_le(e, tol), "%s: T5 Phi_tilde not Hermicity-preserving", name);
        acb_mat_clear(d2); acb_mat_clear(rhs); acb_mat_clear(lhs);
        acb_mat_clear(Xd); acb_mat_clear(X);
    }

    printf("  %s: T3 ||S_tilde-S_Phi||_F=%.2e (eta=0 oracle)\n", name, t3);
    arb_clear(e); arb_clear(tol_def); arb_clear(tol);
    acb_mat_clear(diff); acb_mat_clear(S2); acb_mat_clear(St); acb_mat_clear(S);
}

int main(void)
{
    const slong PREC256 = 256;

    /* ---- T1 (oracle) + T2 (Kronecker) on >=2 channels incl. asymmetric. ---- */
    {
        aic_ucp_kraus phi;
        make_block_cond_exp(&phi, 4, 2);
        T1_oracle("block_cond_exp(4,2)", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
    }
    {   /* asymmetric / complex: conjugate dephasing by a fixed complex unitary */
        aic_ucp_kraus base, conj;
        make_dephasing(&base, 3);
        make_conjugated(&conj, &base, PREC256);
        T1_oracle("conjugated dephasing(3) [complex/asymmetric]", &conj, PREC256);
        aic_ucp_kraus_clear(&conj);
        aic_ucp_kraus_clear(&base);
    }

    /* ---- T3/T4/T5: eta=0 oracle on the exact-idempotent corpus. ---- */
    {
        aic_ucp_kraus phi;
        make_block_cond_exp(&phi, 4, 2);
        T345_idemp_channel("block_cond_exp(4,2)", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
    }
    {
        aic_ucp_kraus phi;
        make_trace_replace(&phi, 3);
        T345_idemp_channel("trace_replace(3)", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
    }
    {
        aic_ucp_kraus phi;
        make_noiseless_subsystem(&phi, 2, 2);
        T345_idemp_channel("noiseless_subsystem(2,2)", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
    }
    {
        aic_ucp_kraus phi;
        make_identity(&phi, 3);
        T345_idemp_channel("identity(3)", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
    }
    {
        aic_ucp_kraus phi;
        make_dephasing(&phi, 4);
        T345_idemp_channel("dephasing(4)", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
    }
    {
        aic_ucp_kraus phi;
        make_compress_idemp(&phi, 4, 2);
        T345_idemp_channel("compress_idemp(4,2)", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
    }
    {   /* the asymmetric/complex idempotent (transpose-non-closed image) */
        aic_ucp_kraus base, conj;
        make_dephasing(&base, 3);
        make_conjugated(&conj, &base, PREC256);
        T345_idemp_channel("conjugated dephasing(3) [asymmetric]", &conj, PREC256);
        aic_ucp_kraus_clear(&conj);
        aic_ucp_kraus_clear(&base);
    }

    /* ---- T4 teeth: a NON-idempotent channel -> S_tilde must DIFFER from S
     * (and still be idempotent). Mutating aic_assoc_regularize to return S
     * unmodified makes ||S_tilde - S||_F ~ 0 here -> RED on the "must differ"
     * check. Verified, restored. ---- */
    {
        aic_ucp_kraus base, mix, phi;
        make_dephasing(&base, 3);
        make_trace_replace(&mix, 3);
        make_mix(&phi, &base, &mix, 0.10, PREC256);  /* small but nonzero eta */
        slong n = phi.dim_H;
        acb_mat_t S, St, S2, diff;
        arb_t e, tol_id, floor;
        acb_mat_init(S, n * n, n * n);
        acb_mat_init(St, n * n, n * n);
        acb_mat_init(S2, n * n, n * n);
        acb_mat_init(diff, n * n, n * n);
        arb_init(e); arb_init(tol_id); arb_init(floor);
        arb_set_d(tol_id, 1e-22);
        arb_set_d(floor, 1e-3);
        aic_assoc_superop_from_ucp(S, &phi, PREC256);
        aic_assoc_regularize(St, &phi, PREC256);
        /* must DIFFER from S (the regularization actually moved it) */
        acb_mat_sub(diff, St, S, PREC256);
        aic_mat_frobenius_norm(e, diff, PREC256);
        AIC_CHECK_MSG(arb_gt(e, floor),
                      "T4 teeth: S_tilde must differ from S on a non-idempotent "
                      "channel (regularize is a no-op?) ||diff||_F=%.2e",
                      arf_get_d(arb_midref(e), ARF_RND_NEAR));
        /* but STILL be exactly idempotent (this is prop_P's guarantee) */
        acb_mat_sqr(S2, St, PREC256);
        acb_mat_sub(diff, S2, St, PREC256);
        aic_mat_frobenius_norm(e, diff, PREC256);
        AIC_CHECK_MSG(arb_le(e, tol_id),
                      "T4 teeth: S_tilde not idempotent on non-idempotent channel");
        printf("  T4 teeth: non-idempotent channel, ||S_tilde^2-S_tilde||_F=%.2e\n",
               arf_get_d(arb_midref(e), ARF_RND_NEAR));
        arb_clear(floor); arb_clear(tol_id); arb_clear(e);
        acb_mat_clear(diff); acb_mat_clear(S2); acb_mat_clear(St); acb_mat_clear(S);
        aic_ucp_kraus_clear(&phi);
        aic_ucp_kraus_clear(&mix);
        aic_ucp_kraus_clear(&base);
    }

    /* ---- T6: NON-NORMAL sgn cross-check A1 (prop_P) vs A2 (binomial xpow), on
     * TWO structurally-different non-normal in-basin channels (each genuinely
     * non-normal + in-basin is asserted live inside T6_nonnormal). ---- */
    {   /* channel 1: compress_idemp(3,1) (+) dephasing(3), t=0.05. A PROPER-
         * carrier idempotent (non-self-adjoint in HS -> non-normal superoperator)
         * mixed with a DIFFERENT idempotent. (||[S,S^dag]||_F~3.13, 4eta~0.33.) */
        aic_ucp_kraus base, mix, phi;
        make_compress_idemp(&base, 3, 1);
        make_dephasing(&mix, 3);
        make_mix(&phi, &base, &mix, 0.05, PREC256);
        T6_nonnormal("compress(3,1)+dephasing(3) t=.05", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
        aic_ucp_kraus_clear(&mix);
        aic_ucp_kraus_clear(&base);
    }
    {   /* channel 2: STRUCTURALLY DIFFERENT -- compress_idemp(4,1) (+)
         * block_cond_exp(4,2), t=0.05. Different base dimension (d=4 vs 3) and a
         * different mix family (a two-block conditional expectation, not
         * dephasing), so the non-normal superoperator has a different structure.
         * Genuinely non-normal (||[S,S^dag]||_F~4.42) and in-basin (4eta~0.38);
         * A1==A2 to ~3e-74 (verified). */
        aic_ucp_kraus base, mix, phi;
        make_compress_idemp(&base, 4, 1);
        make_block_cond_exp(&mix, 4, 2);
        make_mix(&phi, &base, &mix, 0.05, PREC256);
        T6_nonnormal("compress(4,1)+block_cond_exp(4,2) t=.05", &phi, PREC256);
        aic_ucp_kraus_clear(&phi);
        aic_ucp_kraus_clear(&mix);
        aic_ucp_kraus_clear(&base);
    }

    /* ---- T7: ||Phi_tilde - Phi||_op = O(eta). prop_P (.tex:524-533) only gives
     * ||Phi_tilde - Phi|| <= ||2 Phi - 1|| * O(delta) at the SUPEROPERATOR level;
     * the cb-norm bound ||Phi_tilde - Phi||_cb <= O(eta) (.tex:2179) needs the
     * Watrous SDP, not done here. eta proxy = ||S_Phi^2 - S_Phi||_op, a computable
     * stand-in. T7 is a PRELIMINARY sanity check, NOT the paper's universality
     * claim: that claim is about the ALGEBRA's axiom-defect constant being
     * dimension-independent (Increment 2's canary, aic-3qq/aic-4c7), and
     * ||2 Phi - 1|| itself need not be dimension-independent. T7 shows two honest
     * things:
     *   T7a -- linear-in-defect growth + vanishing at eta->0, on the SMOOTH
     *          dephasing(+)trace_replace family (where the op-norm proxy is clean).
     *   T7b -- the ratio dist/eta stays BOUNDED across a genuinely DIMENSION-
     *          NONTRIVIAL family (compress_idemp(d,1)(+)dephasing(d), d=3,4,6,8),
     *          i.e. does not grow with dim. The earlier d-trivial dephasing(+)trace
     *          family gave byte-IDENTICAL ratios at every d (the hostile review
     *          showed it would pass even with a dim-dependent bug); T7b uses a
     *          family whose non-normality ||[S,S^dag]||_F GROWS with d, so the
     *          ratio is genuinely MEASURED, not constant-by-construction. ---- */

    /* T7a: smooth family, ONE dim, linear growth + vanishing. */
    {
        const double ts[] = {0.02, 0.05, 0.10, 0.15, 0.18};
        const int nt = (int) (sizeof(ts) / sizeof(ts[0]));
        double prev_dist = -1.0;       /* monotone-increasing dist witness */
        double ratio_lo = 1e300, ratio_hi = 0.0;
        for (int it = 0; it < nt; it++) {
            aic_ucp_kraus base, mix, phi;
            make_dephasing(&base, 3);
            make_trace_replace(&mix, 3);
            make_mix(&phi, &base, &mix, ts[it], PREC256);
            slong n = phi.dim_H;
            acb_mat_t S, St, diff;
            acb_mat_init(S, n * n, n * n);
            acb_mat_init(St, n * n, n * n);
            acb_mat_init(diff, n * n, n * n);
            aic_assoc_superop_from_ucp(S, &phi, PREC256);
            double eta = defect_opnorm(S, PREC256);
            aic_assoc_regularize(St, &phi, PREC256);
            acb_mat_sub(diff, St, S, PREC256);
            double dist = midpoint_opnorm(diff, PREC256);
            double ratio = dist / eta;
            if (ratio < ratio_lo) ratio_lo = ratio;
            if (ratio > ratio_hi) ratio_hi = ratio;
            printf("  T7a d=3 t=%.2f: eta_proxy=%.3e ||Phi_tilde-Phi||_op=%.3e "
                   "ratio=%.3f\n", ts[it], eta, dist, ratio);
            /* O(eta): ratio bounded by a chosen constant (prop_P gives O(delta);
             * the multiplier is empirical headroom, not a derived bound). */
            AIC_CHECK_MSG(ratio <= 10.0,
                          "T7a: ||Phi_tilde-Phi||_op/eta=%.3f exceeds 10 (t=%.2f)",
                          ratio, ts[it]);
            /* linear-in-defect: dist grows monotonically with t (eta grows with t,
             * dist tracks it). The smallest t gives the smallest dist (-> 0). */
            AIC_CHECK_MSG(dist > prev_dist,
                          "T7a: dist not increasing with t (t=%.2f dist=%.3e "
                          "prev=%.3e) -- not linear-in-defect", ts[it], dist, prev_dist);
            prev_dist = dist;
            acb_mat_clear(diff); acb_mat_clear(St); acb_mat_clear(S);
            aic_ucp_kraus_clear(&phi);
            aic_ucp_kraus_clear(&mix);
            aic_ucp_kraus_clear(&base);
        }
        printf("  T7a (vanishing+linear): ratio in [%.3f,%.3f], dist monotone up\n",
               ratio_lo, ratio_hi);
    }

    /* T7b: DIMENSION-NONTRIVIAL family, ratio stays bounded across d=3,4,6,8. */
    {
        const slong dims[] = {3, 4, 6, 8};
        const double ts[] = {0.02, 0.04};   /* both in-basin at every dim (4eta<0.45) */
        const int nd = 4, nt = 2;
        double ratio_hi = 0.0, ratio_lo = 1e300;
        double cnorm_lo = 1e300, cnorm_hi = 0.0;    /* non-normality vs dim */
        double ratio_d3 = -1.0, ratio_dmax = -1.0;  /* nontriviality witness (t=0.02) */
        for (int di = 0; di < nd; di++) {
            slong d = dims[di];
            for (int it = 0; it < nt; it++) {
                aic_ucp_kraus base, mix, phi;
                make_compress_idemp(&base, d, 1);   /* non-normal, dim-scaling */
                make_dephasing(&mix, d);            /* a DIFFERENT idempotent */
                make_mix(&phi, &base, &mix, ts[it], PREC256);
                slong n = phi.dim_H;
                acb_mat_t S, St, diff;
                acb_mat_init(S, n * n, n * n);
                acb_mat_init(St, n * n, n * n);
                acb_mat_init(diff, n * n, n * n);
                aic_assoc_superop_from_ucp(S, &phi, PREC256);
                double eta = defect_opnorm(S, PREC256);
                double cn = commutator_normF(S, PREC256);
                /* must be in-basin so regularize does not abort (Rule 4 guard). */
                AIC_CHECK_MSG(4.0 * eta < 1.0,
                              "T7b: out of basin 4||S^2-S||_op=%.3f (d=%ld t=%.2f)",
                              4.0 * eta, (long) d, ts[it]);
                aic_assoc_regularize(St, &phi, PREC256);
                acb_mat_sub(diff, St, S, PREC256);
                double dist = midpoint_opnorm(diff, PREC256);
                double ratio = dist / eta;
                if (ratio < ratio_lo) ratio_lo = ratio;
                if (ratio > ratio_hi) ratio_hi = ratio;
                if (cn < cnorm_lo) cnorm_lo = cn;
                if (cn > cnorm_hi) cnorm_hi = cn;
                if (it == 0 && d == 3) ratio_d3 = ratio;
                if (it == 0 && d == 8) ratio_dmax = ratio;
                printf("  T7b d=%ld t=%.2f: eta=%.3e dist=%.3e ratio=%.3f "
                       "|[S,Sd]|_F=%.2f 4eta=%.3f\n",
                       (long) d, ts[it], eta, dist, ratio, cn, 4.0 * eta);
                /* bounded constant across dims (preliminary universality canary):
                 * a constant ~ dim would blow this up at d=8 vs d=3. */
                AIC_CHECK_MSG(ratio <= 5.0,
                              "T7b: ratio=%.3f exceeds 5 (d=%ld t=%.2f) -- constant "
                              "growing with dim?", ratio, (long) d, ts[it]);
                acb_mat_clear(diff); acb_mat_clear(St); acb_mat_clear(S);
                aic_ucp_kraus_clear(&phi);
                aic_ucp_kraus_clear(&mix);
                aic_ucp_kraus_clear(&base);
            }
        }
        /* INTEGRITY: the family must be genuinely dimension-nontrivial, else this
         * is the same toothless canary the review flagged. (i) EVERY dim is
         * genuinely non-normal (cnorm_lo > 1e-3 -- the half that catches a NORMAL
         * family like dephasing(+)trace, whose cnorm is 0 at every d; a bare
         * cnorm_hi >= 2*cnorm_lo would be vacuously true at 0>=0) AND the non-
         * normality GROWS with dim (||[S,S^dag]||_F at d=8 >= 2x its d=3 value,
         * 3.3->10.2); (ii) the ratio is genuinely MEASURED, not constant-by-
         * construction: the d=8 t=0.02 ratio differs from the d=3 t=0.02 ratio
         * (0.52 vs 1.02). Both halves were mutation-proven (swap to dephasing(+)
         * trace -> (i) RED on cnorm_lo and (ii) RED on the byte-identical ratios). */
        AIC_CHECK_MSG(cnorm_lo > 1e-3 && cnorm_hi >= 2.0 * cnorm_lo,
                      "T7b INTEGRITY: family not dim-nontrivial -- ||[S,S^dag]||_F "
                      "band [%.2f,%.2f] (need lo>1e-3 AND hi>=2*lo)", cnorm_lo, cnorm_hi);
        AIC_CHECK_MSG(ratio_d3 > 0.0 && ratio_dmax > 0.0 &&
                          fabs(ratio_dmax - ratio_d3) > 1e-3,
                      "T7b INTEGRITY: ratio is constant-by-construction across dim "
                      "(d=3 %.4f == d=8 %.4f) -- toothless canary", ratio_d3, ratio_dmax);
        /* finite positive band sanity. */
        AIC_CHECK_MSG(ratio_hi < 1e9 && ratio_lo > 0.0,
                      "T7b: ratio band [%.3f,%.3f] not finite positive",
                      ratio_lo, ratio_hi);
        printf("  T7b summary: dist/eta in [%.3f,%.3f] across d=3,4,6,8 (nontrivial: "
               "|[S,Sd]|_F %.2f->%.2f, ratio@t=.02 d3=%.3f d8=%.3f). PRELIMINARY -- "
               "the definitive axiom-defect universality canary is Increment 2.\n",
               ratio_lo, ratio_hi, cnorm_lo, cnorm_hi, ratio_d3, ratio_dmax);
    }

    /* ---- T8: precision ladder, Phi_tilde @53 vs @256. ---- */
    {
        aic_ucp_kraus base, mix, phi;
        make_dephasing(&base, 3);
        make_trace_replace(&mix, 3);
        make_mix(&phi, &base, &mix, 0.10, PREC256);
        slong n = phi.dim_H;
        acb_mat_t S53, S256, diff;
        arb_t e, tol;
        acb_mat_init(S53, n * n, n * n);
        acb_mat_init(S256, n * n, n * n);
        acb_mat_init(diff, n * n, n * n);
        arb_init(e); arb_init(tol);
        arb_set_d(tol, 1e-12);
        aic_assoc_regularize(S53, &phi, 53);
        aic_assoc_regularize(S256, &phi, PREC256);
        acb_mat_sub(diff, S53, S256, PREC256);
        aic_mat_frobenius_norm(e, diff, PREC256);
        AIC_CHECK_MSG(arb_le(e, tol),
                      "T8: prec=53 vs 256 disagree ||diff||_F=%.3e",
                      arf_get_d(arb_midref(e), ARF_RND_NEAR));
        printf("  T8 precision ladder @53 vs @256: %.3e\n",
               arf_get_d(arb_midref(e), ARF_RND_NEAR));
        arb_clear(tol); arb_clear(e);
        acb_mat_clear(diff); acb_mat_clear(S256); acb_mat_clear(S53);
        aic_ucp_kraus_clear(&phi);
        aic_ucp_kraus_clear(&mix);
        aic_ucp_kraus_clear(&base);
    }

    /* T9 FAIL-LOUD basin: aic_prop_P asserts ||S^2-S||_op < 1/4 on a certified
     * ball and aborts otherwise (src/aic_funcalc_proj.c). It is documented and
     * exercised by ALGORITHM.md; not invoked here (an abort would fail the test
     * binary). The in-basin guard 4||S^2-S||_op<1 IS asserted live in T6. */

    aic_test_report("test_assoc");
    return 0;
}
