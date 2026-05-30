/* test_corner.c — cross-checks for the corner module (bead aic "corner",
 * Increment 1): the compression maps Co_{P,Q} and the corner subspaces S_{P,Q}
 * of section 7 (approximate_algebras.tex:1052-1076). Built on the eps-C* algebra
 * data model (aic_ecstar) produced from an eta-idempotent UCP map by the assoc
 * module (aic_assoc_ecstar_from_phi). Reuses the exact-idempotent channels in
 * tests/test_idemp.h (the eta=0 oracle fixtures).
 *
 * Cross-check ladder:
 *   T1 Co^2 = Co (.tex:1061): the theta-snap makes Co an EXACT idempotent --
 *      ||Co^2 - Co||_op machine-zero. MUTATION: replacing aic_prop_P(M) with M
 *      (skip the theta snap) leaves Co = 1/2(L_P R_Q + R_Q L_P), idempotent only
 *      O(delta+eps) -> RED at eta>0; at eta=0 M is already idempotent so use the
 *      eta>0 instance for the mutation. Verified RED on the eta>0 fixture.
 *   T2 involution (.tex:1062): Co_{P,Q}(X)^dag = Co_{Q,P}(X^dag) for each basis
 *      B_k. MUTATION: build Co_{Q,P} as Co_{P,Q} (swap P,Q dropped) -> RED on a
 *      P != Q instance.
 *   T3 (.tex:1064): ||L_P R_Q - Co|| and ||R_Q L_P - Co|| <= c(delta+eps); report
 *      the measured c. (Operator norms of the d x d coordinate matrices.)
 *   T4 almost-containment (.tex:1074): X in S_{P,Q} (= Img Co) is EXACT-fixed by
 *      Co (machine-zero); and the genuine P1 subset P case ||Co_{P,Q}(X)-X|| <=
 *      c(delta+eps)||X|| for X in S_{P1,Q1}.
 *   T5 eta=0 ORACLE (headline, Rule 6): A = M_n from the identity channel
 *      (Phi(X)=X => Phi_tilde=Phi=id, A=M_n, star=plain product). EXACT diagonal
 *      projectors P,Q. (a) apply(Co_{P,Q},X) = P X Q EXACTLY (machine-zero) for
 *      each basis X -- the cleanest oracle; (b) dim S_{P,Q} = rank(P) rank(Q);
 *      (c) S_{P,Q} spans the same subspace as {P E_ab Q} (gauge-free via the
 *      projector Pi_S = sum vec(C_m) vec(C_m)^dag). MUTATION: swap P,Q in the Co
 *      build -> Co gives Q X P, oracle (a) RED.
 *   T6 degenerate dims (.tex:1066): dim S_P = 0 for P~0, = dim_A for P~I, = 1 for
 *      rank-1 P (in the eta=0 M_n algebra).
 *   T7 OBLIQUE oracle (the suite's teeth): T1-T6 are all DEGENERATE for two
 *      correctness questions -- build_L's STAR vs plain product (identity- and
 *      orthogonal-Phi_tilde fixtures make L_P star==plain in coords), and
 *      extract's LEFT vs RIGHT singular vectors (every Co in T1-T6 is an
 *      ORTHOGONAL projector, left==right). compress_idemp(4,2) has a genuinely
 *      OBLIQUE Phi_tilde (S_tilde sigma_max = sqrt(3)) and is EXACTLY idempotent;
 *      with P=|e0+e1><e0+e1|/2, Q=|e0><e0| the compression Co is oblique
 *      (sigma_max = 2/sqrt(3) > 1). (a) independent build_L oracle via
 *      aic_assoc_superop_apply + a non-vacuity gap kills the star mutation;
 *      (b) the oblique-Co fixed-point kills the right-vs-left mutation. Both
 *      mutation-proven RED (see test_t7's docstring).
 *
 * The eta>0 instances (T1 mutation, T3, T4-genuine) use the assoc dep+conj family
 * the assoc module already exercises. The eta=0 oracle uses the identity channel;
 * the oblique oracle (T7) uses compress_idemp.
 */
#include <complex.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <sys/wait.h>
#include <unistd.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
#include "aic_corner.h"
#include "aic_ecstar.h"
#include "aic_mat.h"
#include "aic_test.h"
#include "aic_ucp.h"
#include "test_idemp.h"

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* <Bi, Y>_F = Tr(Bi^dag Y) = sum_{a,b} conj(Bi[a,b]) Y[a,b], computed INLINE in
 * the test so T7's build_L oracle is independent of aic_corner_frob_ip (the
 * production internal). `out` is a caller-init'd acb_t; Bi, Y are n x n. */
static void frob_ip(acb_t out, const acb_mat_t Bi, const acb_mat_t Y, slong prec)
{
    slong n = acb_mat_nrows(Y);
    acb_t t;
    acb_init(t);
    acb_zero(out);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++) {
            acb_conj(t, acb_mat_entry(Bi, a, b));
            acb_mul(t, t, acb_mat_entry(Y, a, b), prec);
            acb_add(out, out, t, prec);
        }
    acb_clear(t);
}

/* Suppress -Wunused-function for the test_idemp.h channel builders this TU does
 * not exercise (it uses make_identity / make_dephasing / make_conjugated and,
 * since T7, make_compress_idemp). Never called; references the addresses so the
 * remaining static helpers are "used". */
__attribute__((unused))
static void suppress_unused_idemp_builders(void)
{
    (void) make_block_cond_exp;
    (void) make_trace_replace;
    (void) make_noiseless_subsystem;
}

/* ||A||_op of an acb_mat (certified ball midpoint). */
static double opnorm_d(const acb_mat_t Aop, slong prec)
{
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, Aop, prec);
    double r = dd(e);
    arb_clear(e);
    return r;
}

/* ||A||_op on the BALL MIDPOINTS (mirrors eta_proxy below). For a NEAR-ZERO
 * difference matrix carrying wide accumulated error balls, the certified
 * aic_mat_opnorm forms the Gram A^dag A whose off-diagonal relative-Hermiticity
 * check can false-fail (bead aic-2yo); collapsing to midpoints first is the
 * test-side measurement (we only want the number, not a certified ball, exactly
 * as eta_proxy does). Used by T8/T9 where cdot outputs can be tiny. */
static double mid_opnorm_d(const acb_mat_t Aop, slong prec)
{
    slong r = acb_mat_nrows(Aop), c = acb_mat_ncols(Aop);
    acb_mat_t M;
    acb_mat_init(M, r, c);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(Aop, i, j));
    double v = opnorm_d(M, prec);
    acb_mat_clear(M);
    return v;
}

/* |z| of a complex acb scalar as a double (ball midpoint). Used by the Increment-2b
 * tests where the inner product / scalar differences are complex numbers. */
static double abs_acb(const acb_t z)
{
    arb_t m;
    arb_init(m);
    acb_abs(m, z, 53);
    double r = dd(m);
    arb_clear(m);
    return r;
}

/* arb ball certainly <= t? (mirrors AIC_CHECK arb_le idiom for a double bound). */
static int arb_le_dbl(const arb_t x, double t)
{
    arb_t tt;
    arb_init(tt);
    arb_set_d(tt, t);
    int r = arb_le(x, tt);
    arb_clear(tt);
    return r;
}

/* n x n diagonal projector with the first `k` diagonal entries 1 (rest 0). */
static void diag_proj(acb_mat_t P, slong n, slong k)
{
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    for (slong i = 0; i < k; i++) acb_set_si(acb_mat_entry(P, i, i), 1);
}

/* a single diagonal projector |e_idx><e_idx| (rank 1). */
static void diag_proj_at(acb_mat_t P, slong n, slong idx)
{
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    acb_set_si(acb_mat_entry(P, idx, idx), 1);
}

/* The eta>0 dep(d)+conj(dep(d)) family (mirrors test_assoc2 make_eta_family),
 * inline to avoid a cross-TU dependency: Phi_t = (1-t) dephasing(d) + t
 * conj(dephasing(d)). Genuinely eta>0 (assoc nonzero). `out` init'd here. */
static void make_eta_family(aic_ucp_kraus *out, slong d, double t, slong prec)
{
    aic_ucp_kraus base, m0, mix;
    make_dephasing(&base, d);
    make_dephasing(&m0, d);
    make_conjugated(&mix, &m0, prec);
    /* convex mix as a Kraus union {sqrt(1-t) Ka} U {sqrt(t) Lb}. */
    slong n = base.dim_H;
    aic_ucp_kraus_init(out, n, n, base.r + mix.r);
    arb_t s;
    arb_init(s);
    arb_set_d(s, sqrt(1.0 - t));
    for (slong a = 0; a < base.r; a++)
        acb_mat_scalar_mul_arb(out->K[a], base.K[a], s, prec);
    arb_set_d(s, sqrt(t));
    for (slong b = 0; b < mix.r; b++)
        acb_mat_scalar_mul_arb(out->K[base.r + b], mix.K[b], s, prec);
    arb_clear(s);
    aic_ucp_kraus_clear(&mix);
    aic_ucp_kraus_clear(&m0);
    aic_ucp_kraus_clear(&base);
}

/* eta proxy = ||S_Phi^2 - S_Phi||_op (the idempotence defect of the superop,
 * == ||Phi^2-Phi||_op <= eta; mirrors test_assoc2 eta_proxy). */
static double eta_proxy(const aic_ucp_kraus *phi, slong prec)
{
    slong n = phi->dim_H, nn = n * n;
    acb_mat_t S, S2, D, M;
    arb_t e;
    acb_mat_init(S, nn, nn);
    acb_mat_init(S2, nn, nn);
    acb_mat_init(D, nn, nn);
    acb_mat_init(M, nn, nn);
    arb_init(e);
    aic_assoc_superop_from_ucp(S, phi, prec);
    acb_mat_sqr(S2, S, prec);
    acb_mat_sub(D, S2, S, prec);
    for (slong i = 0; i < nn; i++)
        for (slong j = 0; j < nn; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(D, i, j));
    aic_mat_opnorm(e, M, prec);
    double o = dd(e);
    arb_clear(e);
    acb_mat_clear(M); acb_mat_clear(D); acb_mat_clear(S2); acb_mat_clear(S);
    return o;
}

/* =========================== T5: eta=0 oracle ============================ *
 * A = M_n from identity(n). P, Q exact diagonal projectors. Co_{P,Q}(X) = PXQ
 * exactly; dim S = rank(P) rank(Q); S spans span{P E_ab Q}. */
static void t5_oracle(slong n, slong rkP, slong rkQ, slong prec)
{
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    AIC_CHECK_MSG(ae.A.dim_A == n * n,
                  "T5(n=%ld): identity dim_A=%ld != n^2", (long) n,
                  (long) ae.A.dim_A);

    acb_mat_t P, Q;
    diag_proj(P, n, rkP);       /* P = diag(1..1,0..) rank rkP */
    /* Q = diag(0..0,1..1) supported on the LAST rkQ coords (disjoint from P when
     * rkP+rkQ<=n; for n=2 P=diag(1,0) Q=diag(0,1)). */
    acb_mat_init(Q, n, n);
    acb_mat_zero(Q);
    for (slong i = n - rkQ; i < n; i++) acb_set_si(acb_mat_entry(Q, i, i), 1);

    acb_mat_t Co;
    acb_mat_init(Co, n * n, n * n);
    aic_corner_Co(Co, &ae.A, P, Q, prec);

    /* (a) apply(Co,X) = P X Q EXACTLY for each basis X. */
    double maxdiff = 0.0;
    acb_mat_t coX, PXQ, XQ, diff;
    acb_mat_init(coX, n, n);
    acb_mat_init(PXQ, n, n);
    acb_mat_init(XQ, n, n);
    acb_mat_init(diff, n, n);
    for (slong k = 0; k < ae.A.dim_A; k++) {
        aic_corner_apply(coX, &ae.A, Co, ae.A.B[k], prec);
        acb_mat_mul(XQ, ae.A.B[k], Q, prec);   /* X Q */
        acb_mat_mul(PXQ, P, XQ, prec);         /* P X Q */
        acb_mat_sub(diff, coX, PXQ, prec);
        double dnorm = opnorm_d(diff, prec);
        if (dnorm > maxdiff) maxdiff = dnorm;
    }
    AIC_CHECK_MSG(maxdiff < 1e-9,
                  "T5(n=%ld): Co_{P,Q}(X) != PXQ (max ||diff||=%.3e)", (long) n,
                  maxdiff);

    /* (b) dim S = rank(P) rank(Q). */
    acb_mat_t *C;
    slong dim_S;
    aic_corner_extract(&C, &dim_S, Co, &ae.A, prec);
    AIC_CHECK_MSG(dim_S == rkP * rkQ,
                  "T5(n=%ld): dim S=%ld != rank(P)rank(Q)=%ld", (long) n,
                  (long) dim_S, (long) (rkP * rkQ));

    /* (c) S spans span{P E_ab Q} (gauge-free). Build Pi_S = sum_m vec(C_m)
     * vec(C_m)^dag (row-major vec) and Pi_ref from an orthonormal basis of
     * {P E_ab Q}. Compare ||Pi_S - Pi_ref||_op ~ 0. The reference set P E_ab Q
     * (a,b run over rows-of-P-support x cols-of-Q-support) is exactly the rkP rkQ
     * matrices with support (row in [0,rkP), col in [n-rkQ,n)), already
     * orthonormal as matrix units. */
    slong nn = n * n;
    acb_mat_t PiS, PiR, Diff;
    acb_mat_init(PiS, nn, nn);
    acb_mat_init(PiR, nn, nn);
    acb_mat_init(Diff, nn, nn);
    acb_mat_zero(PiS);
    acb_mat_zero(PiR);
    acb_t tip;
    acb_init(tip);
    /* Pi_S = sum_m vec(C_m) vec(C_m)^dag, vec_r index = a*n+b. */
    for (slong m = 0; m < dim_S; m++)
        for (slong r = 0; r < nn; r++)
            for (slong c = 0; c < nn; c++) {
                acb_conj(tip, acb_mat_entry(C[m], c / n, c % n));
                acb_mul(tip, acb_mat_entry(C[m], r / n, r % n), tip, prec);
                acb_add(acb_mat_entry(PiS, r, c), acb_mat_entry(PiS, r, c), tip, prec);
            }
    /* Pi_R = sum over the matrix units E_{a,b}, a in [0,rkP), b in [n-rkQ,n):
     * vec(E_ab) is the unit basis vector at index a*n+b, so Pi_R is diagonal with
     * 1s at those indices. */
    for (slong a = 0; a < rkP; a++)
        for (slong b = n - rkQ; b < n; b++)
            acb_set_si(acb_mat_entry(PiR, a * n + b, a * n + b), 1);
    acb_mat_sub(Diff, PiS, PiR, prec);
    double spand = opnorm_d(Diff, prec);
    AIC_CHECK_MSG(spand < 1e-9,
                  "T5(n=%ld): S != span{P E_ab Q} (||Pi_S - Pi_ref||=%.3e)",
                  (long) n, spand);

    printf("T5 oracle n=%ld P(rk%ld) Q(rk%ld): max||Co(X)-PXQ||=%.2e dim_S=%ld "
           "(=%ld) span_diff=%.2e\n", (long) n, (long) rkP, (long) rkQ, maxdiff,
           (long) dim_S, (long) (rkP * rkQ), spand);

    acb_clear(tip);
    acb_mat_clear(Diff); acb_mat_clear(PiR); acb_mat_clear(PiS);
    for (slong m = 0; m < dim_S; m++) acb_mat_clear(C[m]);
    if (C) flint_free(C);
    acb_mat_clear(diff); acb_mat_clear(XQ); acb_mat_clear(PXQ); acb_mat_clear(coX);
    acb_mat_clear(Co);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

static void test_t5(void)
{
    t5_oracle(2, 1, 1, 256);   /* P=diag(1,0) Q=diag(0,1), dim S = 1 */
    t5_oracle(3, 2, 1, 256);   /* P=diag(1,1,0) Q=diag(0,0,1), dim S = 2 */
}

/* ===================== T1: Co^2 = Co (exact idempotent) ================== */
/* On both the eta=0 M_n algebra and an eta>0 instance: ||Co^2-Co||_op ~ 0. */
static void t1_idempotent(const char *name, const aic_ecstar *A,
                          const acb_mat_t P, const acb_mat_t Q, slong prec)
{
    slong d = A->dim_A;
    acb_mat_t Co, Co2, Diff;
    acb_mat_init(Co, d, d);
    acb_mat_init(Co2, d, d);
    acb_mat_init(Diff, d, d);
    aic_corner_Co(Co, A, P, Q, prec);
    acb_mat_sqr(Co2, Co, prec);
    acb_mat_sub(Diff, Co2, Co, prec);
    double e = opnorm_d(Diff, prec);
    printf("T1 %-22s ||Co^2-Co||_op=%.3e\n", name, e);
    AIC_CHECK_MSG(e < 1e-9, "T1 %s: Co not idempotent (%.3e) -- theta snap missing?",
                  name, e);
    acb_mat_clear(Diff); acb_mat_clear(Co2); acb_mat_clear(Co);
}

static void test_t1(void)
{
    const slong P = 256;
    /* eta=0: A = M_3 from identity(3); P=diag(1,1,0), Q=diag(0,0,1). */
    {
        aic_ucp_kraus phi;
        make_identity(&phi, 3);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, P);
        acb_mat_t Pp, Qq;
        diag_proj(Pp, 3, 2);
        acb_mat_init(Qq, 3, 3); acb_mat_zero(Qq);
        acb_set_si(acb_mat_entry(Qq, 2, 2), 1);
        t1_idempotent("eta=0 M_3", &ae.A, Pp, Qq, P);
        acb_mat_clear(Qq); acb_mat_clear(Pp);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    /* eta>0: dep(4)+conj(dep(4)); P,Q the algebra's own delta-projectors. We use
     * exact diagonal projectors of the AMBIENT space (Hermitian; build_L/R assert
     * Hermiticity, not membership -- the construction is well-defined for any
     * Hermitian P,Q, and Co is still an EXACT idempotent after theta). */
    {
        aic_ucp_kraus phi;
        make_eta_family(&phi, 4, 0.06, P);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, P);
        acb_mat_t Pp, Qq;
        diag_proj(Pp, 4, 2);
        acb_mat_init(Qq, 4, 4); acb_mat_zero(Qq);
        for (slong i = 2; i < 4; i++) acb_set_si(acb_mat_entry(Qq, i, i), 1);
        t1_idempotent("eta>0 dep+conj(4)", &ae.A, Pp, Qq, P);
        acb_mat_clear(Qq); acb_mat_clear(Pp);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* Conjugate a Hermitian projector P by a fixed unitary V: P' = V^dag P V (still
 * an exact projector, but NON-diagonal -- so the L_P / R_Q asymmetry is genuinely
 * exercised, giving T2 teeth against an L<->R swap that diagonal P,Q would miss).
 * `Pout` init'd here. */
static void conj_proj(acb_mat_t Pout, const acb_mat_t P, const acb_mat_t V,
                      slong prec)
{
    slong n = acb_mat_nrows(P);
    acb_mat_init(Pout, n, n);
    acb_mat_t Vd, t;
    acb_mat_init(Vd, n, n);
    acb_mat_init(t, n, n);
    acb_mat_conjugate_transpose(Vd, V);
    acb_mat_mul(t, Vd, P, prec);        /* V^dag P */
    acb_mat_mul(Pout, t, V, prec);      /* V^dag P V */
    acb_mat_clear(t);
    acb_mat_clear(Vd);
}

/* ===================== T2: involution Co_{P,Q}(X)^dag = Co_{Q,P}(X^dag) === */
static void test_t2(void)
{
    const slong prec = 256;
    /* eta=0 M_4 (identity), NON-DIAGONAL projectors P = V^dag diag(1,1,0,0) V
     * (rank2), Q = V^dag |e2><e2| V (rank1), P != Q. The unitary rotation makes
     * P,Q non-diagonal so L_P and R_Q genuinely differ (a diagonal P,Q algebra is
     * too symmetric and would NOT catch an L<->R swap in build_R). */
    aic_ucp_kraus phi;
    make_identity(&phi, 4);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong n = 4, d = ae.A.dim_A;
    acb_mat_t Pdiag, Qdiag, V, P, Q;
    diag_proj(Pdiag, n, 2);
    diag_proj_at(Qdiag, n, 2);
    acb_mat_init(V, n, n);
    make_fixed_unitary(V, n, prec);
    conj_proj(P, Pdiag, V, prec);
    conj_proj(Q, Qdiag, V, prec);
    acb_mat_clear(V); acb_mat_clear(Qdiag); acb_mat_clear(Pdiag);

    acb_mat_t CoPQ, CoQP;
    acb_mat_init(CoPQ, d, d);
    acb_mat_init(CoQP, d, d);
    aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
    aic_corner_Co(CoQP, &ae.A, Q, P, prec);

    double maxdiff = 0.0, maxact = 0.0, maxcf = 0.0;
    acb_mat_t lhs, lhsd, Xd, rhs, diff, QXdP, cf, XdP;
    acb_mat_init(lhs, n, n);
    acb_mat_init(lhsd, n, n);
    acb_mat_init(Xd, n, n);
    acb_mat_init(rhs, n, n);
    acb_mat_init(diff, n, n);
    acb_mat_init(QXdP, n, n);
    acb_mat_init(cf, n, n);
    acb_mat_init(XdP, n, n);
    for (slong k = 0; k < d; k++) {
        aic_corner_apply(lhs, &ae.A, CoPQ, ae.A.B[k], prec);  /* Co_{P,Q}(B_k) */
        acb_mat_conjugate_transpose(lhsd, lhs);               /* Co_{P,Q}(B_k)^dag */
        acb_mat_conjugate_transpose(Xd, ae.A.B[k]);           /* B_k^dag */
        aic_corner_apply(rhs, &ae.A, CoQP, Xd, prec);         /* Co_{Q,P}(B_k^dag) */
        acb_mat_sub(diff, lhsd, rhs, prec);                   /* RELATION: == 0 */
        double dn = opnorm_d(diff, prec);
        if (dn > maxdiff) maxdiff = dn;
        /* non-vacuity: the maps are genuinely nonzero (else the relation is
         * trivially satisfied 0=0 and catches nothing). */
        double an = opnorm_d(rhs, prec);
        if (an > maxact) maxact = an;
        /* CLOSED-FORM cross-check (eta=0): Co_{Q,P}(X^dag) = Q X^dag P (since
         * Co_{Q,P}(Y)=QYP in M_n). This makes T2 catch a wrong Co (e.g. a P<->Q
         * swap in the build) -- the relation comparison alone is structurally
         * robust (a symmetric perturbation cancels on both sides), the closed form
         * is not. */
        acb_mat_mul(XdP, Xd, P, prec);          /* X^dag P */
        acb_mat_mul(QXdP, Q, XdP, prec);        /* Q X^dag P */
        acb_mat_sub(cf, rhs, QXdP, prec);
        double cn = opnorm_d(cf, prec);
        if (cn > maxcf) maxcf = cn;
    }
    printf("T2 involution eta=0 M_4: max ||Co_PQ(X)^dag - Co_QP(X^dag)||=%.3e "
           "(closed-form ||Co_QP(X^dag)-Q X^dag P||=%.2e, max act=%.3f)\n",
           maxdiff, maxcf, maxact);
    AIC_CHECK_MSG(maxdiff < 1e-9,
                  "T2: involution Co_PQ(X)^dag != Co_QP(X^dag) (%.3e)", maxdiff);
    AIC_CHECK_MSG(maxcf < 1e-9,
                  "T2: Co_QP(X^dag) != Q X^dag P (closed form, %.3e)", maxcf);
    AIC_CHECK_MSG(maxact > 0.1,
                  "T2: Co maps all basis elements to ~0 (relation vacuously 0=0)");

    acb_mat_clear(XdP); acb_mat_clear(cf); acb_mat_clear(QXdP);
    acb_mat_clear(diff); acb_mat_clear(rhs); acb_mat_clear(Xd);
    acb_mat_clear(lhsd); acb_mat_clear(lhs);
    acb_mat_clear(CoQP); acb_mat_clear(CoPQ);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T3: ||L_P R_Q - Co|| <= c(delta+eps) ============== */
static void test_t3(void)
{
    const slong prec = 256;
    printf("T3 ||L_P R_Q - Co|| and ||R_Q L_P - Co|| vs eta:\n");
    double cmax = 0.0;
    const double ts[] = {0.0, 0.02, 0.06};
    for (int it = 0; it < 3; it++) {
        aic_ucp_kraus phi;
        make_eta_family(&phi, 4, ts[it], prec);
        double eta = eta_proxy(&phi, prec);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        slong d = ae.A.dim_A;
        acb_mat_t P, Q;
        diag_proj(P, 4, 2);
        acb_mat_init(Q, 4, 4); acb_mat_zero(Q);
        for (slong i = 2; i < 4; i++) acb_set_si(acb_mat_entry(Q, i, i), 1);

        acb_mat_t Lp, Rq, LR, RL, Co, dLR, dRL;
        acb_mat_init(Lp, d, d); acb_mat_init(Rq, d, d);
        acb_mat_init(LR, d, d); acb_mat_init(RL, d, d);
        acb_mat_init(Co, d, d); acb_mat_init(dLR, d, d); acb_mat_init(dRL, d, d);
        aic_corner_build_L(Lp, &ae.A, P, prec);
        aic_corner_build_R(Rq, &ae.A, Q, prec);
        acb_mat_mul(LR, Lp, Rq, prec);
        acb_mat_mul(RL, Rq, Lp, prec);
        aic_corner_Co(Co, &ae.A, P, Q, prec);
        acb_mat_sub(dLR, LR, Co, prec);
        acb_mat_sub(dRL, RL, Co, prec);
        double eLR = opnorm_d(dLR, prec), eRL = opnorm_d(dRL, prec);
        double cLR = eta > 1e-14 ? eLR / eta : 0.0;
        double cRL = eta > 1e-14 ? eRL / eta : 0.0;
        if (cLR > cmax) cmax = cLR;
        if (cRL > cmax) cmax = cRL;
        printf("  t=%.2f eta=%.4e: ||LR-Co||=%.4e (c=%.3f) ||RL-Co||=%.4e (c=%.3f)\n",
               ts[it], eta, eLR, cLR, eRL, cRL);
        if (ts[it] == 0.0) {
            AIC_CHECK_MSG(eLR < 1e-9 && eRL < 1e-9,
                          "T3: at eta=0, L_P R_Q and R_Q L_P must equal Co (%.3e,%.3e)",
                          eLR, eRL);
        } else {
            /* O(delta+eps): the ratio to eta is a bounded constant (report it). */
            AIC_CHECK_MSG(cLR < 5.0 && cRL < 5.0,
                          "T3: ||LR-Co||/eta=%.3f or ||RL-Co||/eta=%.3f exceeds 5 "
                          "(t=%.2f) -- not O(delta+eps)?", cLR, cRL, ts[it]);
            AIC_CHECK_MSG(eLR > 1e-7 || eRL > 1e-7,
                          "T3: both LR-Co and RL-Co ~0 at eta>0 (t=%.2f) -- "
                          "expected O(eta) gap", ts[it]);
        }
        acb_mat_clear(dRL); acb_mat_clear(dLR); acb_mat_clear(Co);
        acb_mat_clear(RL); acb_mat_clear(LR); acb_mat_clear(Rq); acb_mat_clear(Lp);
        acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("  -> measured c (||LR-Co||/eta, ||RL-Co||/eta) max = %.3f\n", cmax);
}

/* ===================== T4: almost-containment (.tex:1074) ================ */
static void test_t4(void)
{
    const slong prec = 256;
    /* (i) X in S_{P,Q} (= Img Co) is EXACT-fixed (machine-zero): eta=0 M_3. */
    {
        aic_ucp_kraus phi;
        make_identity(&phi, 3);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        acb_mat_t P, Q;
        diag_proj(P, 3, 2);
        diag_proj_at(Q, 3, 2);
        acb_mat_t Co;
        acb_mat_init(Co, 9, 9);
        aic_corner_Co(Co, &ae.A, P, Q, prec);
        acb_mat_t *C;
        slong dim_S;
        aic_corner_extract(&C, &dim_S, Co, &ae.A, prec);
        double maxfix = 0.0;
        acb_mat_t coC, diff;
        acb_mat_init(coC, 3, 3);
        acb_mat_init(diff, 3, 3);
        for (slong m = 0; m < dim_S; m++) {
            aic_corner_apply(coC, &ae.A, Co, C[m], prec);
            acb_mat_sub(diff, coC, C[m], prec);
            double dn = opnorm_d(diff, prec);
            if (dn > maxfix) maxfix = dn;
        }
        printf("T4(i) X in S_{P,Q} fixed: max ||Co(X)-X||=%.3e (dim_S=%ld)\n",
               maxfix, (long) dim_S);
        AIC_CHECK_MSG(maxfix < 1e-9,
                      "T4(i): X in Img Co not exact-fixed (%.3e)", maxfix);
        acb_mat_clear(diff); acb_mat_clear(coC);
        for (slong m = 0; m < dim_S; m++) acb_mat_clear(C[m]);
        if (C) flint_free(C);
        acb_mat_clear(Co);
        acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    /* (ii) genuine almost-containment P1 subset P: in eta=0 M_4, P1=diag(1,0,0,0)
     * rank1 sub-projector of P=diag(1,1,0,0) rank2 (P1 P = P1 exactly), Q1=Q=
     * diag(0,0,1,0). X in S_{P1,Q1} => ||Co_{P,Q}(X)-X|| ~ 0 (eta=0 => the
     * O(delta+eps) bound is 0; the exact-C* corner S_{P1,Q1} subset S_{P,Q}). */
    {
        aic_ucp_kraus phi;
        make_identity(&phi, 4);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        acb_mat_t P, Q, P1;
        diag_proj(P, 4, 2);          /* rank 2 */
        diag_proj_at(Q, 4, 2);       /* |e2><e2| */
        diag_proj_at(P1, 4, 0);      /* |e0><e0|, P1 P = P1 */
        acb_mat_t CoP1Q1, CoPQ;
        slong d = ae.A.dim_A;
        acb_mat_init(CoP1Q1, d, d);
        acb_mat_init(CoPQ, d, d);
        aic_corner_Co(CoP1Q1, &ae.A, P1, Q, prec);
        aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
        acb_mat_t *C;
        slong dim_S1;
        aic_corner_extract(&C, &dim_S1, CoP1Q1, &ae.A, prec);
        double maxcont = 0.0;
        acb_mat_t coC, diff;
        acb_mat_init(coC, 4, 4);
        acb_mat_init(diff, 4, 4);
        for (slong m = 0; m < dim_S1; m++) {
            aic_corner_apply(coC, &ae.A, CoPQ, C[m], prec);  /* Co_{P,Q}(X) */
            acb_mat_sub(diff, coC, C[m], prec);
            double dn = opnorm_d(diff, prec);
            if (dn > maxcont) maxcont = dn;
        }
        printf("T4(ii) S_{P1,Q1} subset S_{P,Q}: max ||Co_PQ(X)-X||=%.3e "
               "(dim S_{P1,Q1}=%ld)\n", maxcont, (long) dim_S1);
        AIC_CHECK_MSG(dim_S1 == 1, "T4(ii): dim S_{P1,Q1}=%ld != 1", (long) dim_S1);
        AIC_CHECK_MSG(maxcont < 1e-9,
                      "T4(ii): S_{P1,Q1} not (almost-)contained in S_{P,Q} (%.3e)",
                      maxcont);
        acb_mat_clear(diff); acb_mat_clear(coC);
        for (slong m = 0; m < dim_S1; m++) acb_mat_clear(C[m]);
        if (C) flint_free(C);
        acb_mat_clear(CoPQ); acb_mat_clear(CoP1Q1);
        acb_mat_clear(P1); acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ===================== T6: degenerate dims (.tex:1066) =================== */
static void test_t6(void)
{
    const slong prec = 256;
    slong n = 3;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;

    /* P ~ 0: dim S_P = 0. */
    {
        acb_mat_t Pzero, Co;
        diag_proj(Pzero, n, 0);     /* all-zero Hermitian projector */
        acb_mat_init(Co, d, d);
        aic_corner_Co_P1(Co, &ae.A, Pzero, prec);   /* Co_{P,1}; S_P uses Co_{P,P}
                                                     * but for P=0, L_P=0 so same */
        acb_mat_t *C;
        slong dim_S;
        aic_corner_extract(&C, &dim_S, Co, &ae.A, prec);
        printf("T6 P~0: dim S_P=%ld (expect 0)\n", (long) dim_S);
        AIC_CHECK_MSG(dim_S == 0, "T6: dim S_P=%ld != 0 for P=0", (long) dim_S);
        if (C) { for (slong m = 0; m < dim_S; m++) acb_mat_clear(C[m]); flint_free(C); }
        acb_mat_clear(Co); acb_mat_clear(Pzero);
    }
    /* P ~ I: dim S_P = dim_A. Co_P = Co_{I,I} = theta(2 L_I R_I - 1); with P=I,
     * L_I R_I = id, so Co = id, S_P = A. */
    {
        acb_mat_t Imat, Co;
        acb_mat_init(Imat, n, n);
        acb_mat_one(Imat);
        acb_mat_init(Co, d, d);
        aic_corner_Co(Co, &ae.A, Imat, Imat, prec);
        acb_mat_t *C;
        slong dim_S;
        aic_corner_extract(&C, &dim_S, Co, &ae.A, prec);
        printf("T6 P~I: dim S_P=%ld (expect dim_A=%ld)\n", (long) dim_S, (long) d);
        AIC_CHECK_MSG(dim_S == d, "T6: dim S_P=%ld != dim_A=%ld for P=I",
                      (long) dim_S, (long) d);
        for (slong m = 0; m < dim_S; m++) acb_mat_clear(C[m]);
        if (C) flint_free(C);
        acb_mat_clear(Co); acb_mat_clear(Imat);
    }
    /* rank-1 P: dim S_P = 1 (one-dimensional delta-projection). Co_P = Co_{P,P}. */
    {
        acb_mat_t P1, Co;
        diag_proj_at(P1, n, 1);     /* |e1><e1| */
        acb_mat_init(Co, d, d);
        aic_corner_Co(Co, &ae.A, P1, P1, prec);
        acb_mat_t *C;
        slong dim_S;
        aic_corner_extract(&C, &dim_S, Co, &ae.A, prec);
        printf("T6 rank-1 P: dim S_P=%ld (expect 1)\n", (long) dim_S);
        AIC_CHECK_MSG(dim_S == 1, "T6: dim S_P=%ld != 1 for rank-1 P", (long) dim_S);
        for (slong m = 0; m < dim_S; m++) acb_mat_clear(C[m]);
        if (C) flint_free(C);
        acb_mat_clear(Co); acb_mat_clear(P1);
    }
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ===================== T7: OBLIQUE oracle (gives the suite teeth) ========= *
 * The T1-T6 fixtures are all DEGENERATE for two correctness questions:
 *   (1) build_L's star vs the plain product. T1(eta=0)/T2/T4/T5/T6 use the
 *       identity channel (Phi_tilde = id => star == plain product), and
 *       T3/T1(eta>0) use the dep+conj family whose Phi_tilde is ORTHOGONAL,
 *       on which the COORDINATE matrix L_P is machine-identical for star vs
 *       plain (the difference lives in A^perp, invisible to <B_i,.>_F). So a
 *       build_L that drops the star (uses acb_mat_mul) is INVISIBLE to T1-T6.
 *   (2) extract's LEFT vs RIGHT singular vectors. Every Co in T1-T6 is an
 *       ORTHOGONAL projector (sigma in {0,1}, left == right singular vectors),
 *       so swapping LEFT for RIGHT (range of Co^dag = co-range) is INVISIBLE.
 * Both blindnesses are closed by ONE oblique fixture: compress_idemp(4,2),
 * whose Phi_tilde is genuinely OBLIQUE (S_tilde sigma_max = sqrt(3) ~ 1.732)
 * and EXACTLY idempotent (a clean oracle, no eta-tolerance fuzz). With
 *   P = projector onto (e0+e1)/sqrt(2) (rank 1, IN M=span(e0,e1)),  Q = |e0><e0|
 * the compression Co_{P,Q} is itself OBLIQUE (measured sigma_max = 2/sqrt(3)
 * ~ 1.1547 > 1), so LEFT != RIGHT singular vectors. dim S_{P,Q} = 1.
 *
 * (a) build_L INDEPENDENT-ORACLE (kills bug 1). L_expected[i,j] = <B_i, Yj>_F
 *     with Yj = Phi_tilde(P B_j) computed via aic_assoc_superop_apply (routes
 *     through Phi_tilde EXACTLY like the star should). ||L_act - L_exp|| ~ 0
 *     (machine-zero, compress_idemp is exactly idempotent). NON-VACUITY GUARD:
 *     ||L_act - L_plain|| (L_plain = <B_i, P B_j>_F, NO Phi_tilde) is O(1)
 *     (~0.667 here), so the oracle check WOULD be RED if build_L dropped the
 *     star. MUTATION-PROVEN: changing aic_corner_build_L's aic_ecstar_star to
 *     plain acb_mat_mul(PBj,P,B_j) made (a) RED at ||L_act-L_exp|| = 6.67e-1
 *     (= the non-vacuity gap), restored.
 * (b) OBLIQUE Co fixed-point (kills bug 2). sigma_max(Co) = 2/sqrt(3) > 1
 *     confirms Co is oblique (else left/right is vacuous). Each extracted
 *     C_m in S_{P,Q} = Img Co is EXACT-fixed: ||Co(C_m) - C_m|| ~ 0 (machine-
 *     zero) -- TRUE iff the LEFT singular vectors (range of Co) were used; the
 *     RIGHT/co-range vectors are NOT fixed by the oblique Co. MUTATION-PROVEN:
 *     changing aic_corner_extract to RIGHT singular vectors (aic_latd_svd(sv,
 *     NULL, Vt, ...) with C_m = sum_k conj(Vt[m*d+k]) B_k) made (b)'s fixed-
 *     point residual jump to 4.08e-1 (RED), restored. */
static void test_t7(void)
{
    const slong prec = 256, n = 4;
    aic_ucp_kraus phi;
    make_compress_idemp(&phi, 4, 2);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;

    /* P = |v><v|, v = (e0+e1)/sqrt2 (rank-1 Hermitian projector in M); Q=|e0><e0|. */
    acb_mat_t P, Q;
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    for (slong i = 0; i < 2; i++)
        for (slong j = 0; j < 2; j++) acb_set_d(acb_mat_entry(P, i, j), 0.5);
    diag_proj_at(Q, n, 0);

    /* (a) independent build_L oracle + non-vacuity gap. */
    acb_mat_t Lact, Lexp, Lplain, PBj, Yj, Dla, Dlp;
    acb_mat_init(Lact, d, d); acb_mat_init(Lexp, d, d); acb_mat_init(Lplain, d, d);
    acb_mat_init(PBj, n, n); acb_mat_init(Yj, n, n);
    acb_mat_init(Dla, d, d); acb_mat_init(Dlp, d, d);
    aic_corner_build_L(Lact, &ae.A, P, prec);
    for (slong j = 0; j < d; j++) {
        acb_mat_mul(PBj, P, ae.A.B[j], prec);                /* P . B_j (plain) */
        aic_assoc_superop_apply(Yj, ae.S_tilde, PBj, prec);  /* Phi_tilde(P B_j) */
        for (slong i = 0; i < d; i++) {
            frob_ip(acb_mat_entry(Lexp, i, j), ae.A.B[i], Yj, prec);
            frob_ip(acb_mat_entry(Lplain, i, j), ae.A.B[i], PBj, prec);
        }
    }
    acb_mat_sub(Dla, Lact, Lexp, prec);
    acb_mat_sub(Dlp, Lact, Lplain, prec);
    double l_oracle = opnorm_d(Dla, prec), l_vac = opnorm_d(Dlp, prec);
    printf("T7 build_L oracle: ||L_act-L_exp||=%.3e  ||L_act-L_plain||=%.4f "
           "(non-vacuity gap)\n", l_oracle, l_vac);
    AIC_CHECK_MSG(l_oracle < 1e-9,
                  "T7(a): build_L != <B_i, Phi_tilde(P B_j)>_F (%.3e) -- star wrong?",
                  l_oracle);
    AIC_CHECK_MSG(l_vac > 0.1,
                  "T7(a): non-vacuity gap ||L_act-L_plain||=%.3e <= 0.1 -- fixture "
                  "not exercising the star (Phi_tilde ~ id on A?)", l_vac);

    /* (b) oblique Co + LEFT-singular-vector fixed-point. */
    acb_mat_t Co;
    acb_mat_init(Co, d, d);
    aic_corner_Co(Co, &ae.A, P, Q, prec);
    double co_smax = opnorm_d(Co, prec);    /* sigma_max(Co) = ||Co||_op */
    printf("T7 Co sigma_max=%.6f (oblique iff > 1)\n", co_smax);
    AIC_CHECK_MSG(co_smax > 1.0 + 1e-6,
                  "T7(b): Co sigma_max=%.6f <= 1 -- Co is orthogonal, left/right "
                  "test VACUOUS; pick a more oblique P,Q", co_smax);
    acb_mat_t *C;
    slong dim_S;
    aic_corner_extract(&C, &dim_S, Co, &ae.A, prec);
    double maxfix = 0.0;
    acb_mat_t coC, diff;
    acb_mat_init(coC, n, n);
    acb_mat_init(diff, n, n);
    for (slong m = 0; m < dim_S; m++) {
        aic_corner_apply(coC, &ae.A, Co, C[m], prec);
        acb_mat_sub(diff, coC, C[m], prec);
        double dn = opnorm_d(diff, prec);
        if (dn > maxfix) maxfix = dn;
    }
    printf("T7 oblique Co fixed-point: dim_S=%ld max ||Co(C_m)-C_m||=%.3e\n",
           (long) dim_S, maxfix);
    AIC_CHECK_MSG(dim_S == 1, "T7(b): dim S_{P,Q}=%ld != 1", (long) dim_S);
    AIC_CHECK_MSG(maxfix < 1e-9,
                  "T7(b): C_m in Img Co not exact-fixed (%.3e) -- RIGHT (co-range) "
                  "singular vectors used instead of LEFT?", maxfix);

    acb_mat_clear(diff); acb_mat_clear(coC);
    for (slong m = 0; m < dim_S; m++) acb_mat_clear(C[m]);
    if (C) flint_free(C);
    acb_mat_clear(Co);
    acb_mat_clear(Dlp); acb_mat_clear(Dla); acb_mat_clear(Yj); acb_mat_clear(PBj);
    acb_mat_clear(Lplain); acb_mat_clear(Lexp); acb_mat_clear(Lact);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ============ Increment 2a: T8 compressed product (.tex:1077-1082) ======== *
 * X . Y = Co_{P,R}(X * Y) : S_{P,Q} x S_{Q,R} -> S_{P,R} (eq compr_prod). The
 * STAR is load-bearing (CLAUDE.md callout). The suite's teeth (mirrors T7):
 *
 * (a) DEFINITION/INDEPENDENT-ORACLE on the OBLIQUE compress_idemp(4,2)
 *     (S_tilde sigma_max = sqrt(3), exactly idempotent => clean oracle). For X,Y
 *     in the relevant corners, aic_corner_cdot must equal the INDEPENDENT route
 *     Co_{P,R}(Phi_tilde(X.Y)) (X.Y plain product, Phi_tilde via
 *     aic_assoc_superop_apply -- the same independent Phi_tilde used in T7).
 *     NON-VACUITY GUARD: the PLAIN-product compression Co_{P,R}(X.Y) (NO
 *     Phi_tilde) differs from cdot by O(1) on this oblique algebra, so the oracle
 *     WOULD be RED if cdot dropped the star. MUTATION-PROVEN: replacing
 *     aic_ecstar_star with acb_mat_mul in aic_corner_cdot makes (a) RED at the
 *     non-vacuity gap (see test_t8 result line), restored.
 * (b) BOUND ||X . Y - X * Y|| <= c(delta+eps) ||X|| ||Y|| (.tex:1081): on the
 *     eta>0 dep+conj(4) family report c; EXACT (machine-zero) at eta=0. The
 *     eta>0 arm is ALSO the COMPRESSION-STEP teeth: at eta>0 the star X*Y =
 *     Phi_tilde(XY) is NOT exactly in S_{P,R}, so cdot = Co_{P,R}(X*Y) DIFFERS
 *     from X*Y by an O(eta) gap. (At eta=0, and for ALL valid corner X,Y, X*Y is
 *     already in S_{P,R} so Co_{P,R} is a no-op on it -- a THEOREM, not a test
 *     gap; this is exactly WHY cdot is "close to the ambient product", .tex:1081.
 *     So the Co_{P,R} step in cdot is invisible at eta=0 but VISIBLE at eta>0:
 *     the X.Y == X*Y check below has teeth precisely there.) MUTATION-PROVEN:
 *     making cdot SKIP the Co_{P,R} apply (return the bare star X*Y) makes the
 *     eta>0 arm RED ("X.Y == X*Y at eta>0 -- expected O(eta) gap"), restored.
 * (c) UNIT (.tex:1082): Ptilde = Co_P(P) is the unit of (S_P,.): ||Ptilde . X -
 *     X||, ||X . Ptilde - X|| <= c(delta+eps) ||X|| for X in S_P. The unit law
 *     needs P a genuine delta-PROJECTION IN A (|v01><v01| is NOT in A -- residual
 *     0.333 -- so its O(delta+eps) bound is vacuous); tested on identity(3)
 *     (P=diag2, exact) and the oblique compress(4,2) with P=|e1><e1| (in A). Both
 *     machine-zero. */

/* a single rank-1 Hermitian projector |v><v| onto the unit vector with entries
 * vr[0..n-1] (real). `Pout` init'd here. */
static void rank1_proj(acb_mat_t Pout, slong n, const double *vr)
{
    acb_mat_init(Pout, n, n);
    acb_mat_zero(Pout);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
            acb_set_d(acb_mat_entry(Pout, i, j), vr[i] * vr[j]);
}

/* T8(a) + the eta=0 arm of (b) on the oblique exactly-idempotent
 * compress_idemp(4,2). The eta>0 arm of (b) is test_t8_bound; (c) is test_t8_unit. */
static void test_t8_oblique(void)
{
    const slong prec = 256, n = 4;
    aic_ucp_kraus phi;
    make_compress_idemp(&phi, 4, 2);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;

    /* P=Q=R=|v><v|, v=(e0+e1)/sqrt2 (rank-1 in M). This is the combo for which
     * the plain product X.Y leaves A substantially (||Phi_tilde(XY)-XY|| = 0.167),
     * so the star-vs-plain gap SURVIVES the Co_{P,R} compression (non-vacuity gap
     * 0.0556, measured). Other rank-1 choices (e.g. v01,e0,e1) wash the star out
     * (XY already near A, gap ~0): the symmetric P=Q=R=|v><v| with v straddling
     * the M/Mperp split is what exercises it (CLAUDE.md: verify obliqueness is
     * genuinely exercised, not washed out by a symmetric combination). */
    double v01[4] = {0.70710678118654752, 0.70710678118654752, 0.0, 0.0};
    acb_mat_t P, Q, R;
    rank1_proj(P, n, v01);
    rank1_proj(Q, n, v01);
    rank1_proj(R, n, v01);

    /* corner bases of S_{P,Q} and S_{Q,R}; pick a representative X, Y. */
    acb_mat_t CoPQ, CoQR, CoPR;
    acb_mat_init(CoPQ, d, d); acb_mat_init(CoQR, d, d); acb_mat_init(CoPR, d, d);
    aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
    aic_corner_Co(CoQR, &ae.A, Q, R, prec);
    aic_corner_Co(CoPR, &ae.A, P, R, prec);
    acb_mat_t *CX, *CY;
    slong dXY, dYR;
    aic_corner_extract(&CX, &dXY, CoPQ, &ae.A, prec);
    aic_corner_extract(&CY, &dYR, CoQR, &ae.A, prec);
    AIC_CHECK_MSG(dXY >= 1 && dYR >= 1,
                  "T8: empty corner (dim S_{P,Q}=%ld, S_{Q,R}=%ld) -- pick P,Q,R "
                  "with nonzero corners", (long) dXY, (long) dYR);

    acb_mat_t X, Y, prod, Zexp, ZexpCo, Zplain, ZplainCo, dorc, dvac;
    acb_mat_init(X, n, n); acb_mat_init(Y, n, n);
    acb_mat_set(X, CX[0]); acb_mat_set(Y, CY[0]);
    acb_mat_init(prod, n, n);
    acb_mat_init(Zexp, n, n); acb_mat_init(ZexpCo, n, n);
    acb_mat_init(Zplain, n, n); acb_mat_init(ZplainCo, n, n);
    acb_mat_init(dorc, n, n); acb_mat_init(dvac, n, n);

    aic_corner_cdot(prod, &ae.A, CoPR, X, Y, prec);          /* production X . Y */
    /* INDEPENDENT: Co_{P,R}(Phi_tilde(X.Y)). */
    acb_mat_t XYplain;
    acb_mat_init(XYplain, n, n);
    acb_mat_mul(XYplain, X, Y, prec);                        /* X . Y (plain)    */
    aic_assoc_superop_apply(Zexp, ae.S_tilde, XYplain, prec);/* Phi_tilde(X.Y)   */
    aic_corner_apply(ZexpCo, &ae.A, CoPR, Zexp, prec);       /* Co_{P,R}(...)    */
    /* NON-VACUITY: Co_{P,R}(X.Y) with NO Phi_tilde. */
    aic_corner_apply(ZplainCo, &ae.A, CoPR, XYplain, prec);
    acb_mat_sub(dorc, prod, ZexpCo, prec);
    acb_mat_sub(dvac, prod, ZplainCo, prec);
    double oracle = mid_opnorm_d(dorc, prec), vac = mid_opnorm_d(dvac, prec);
    printf("T8(a) oblique cdot oracle: ||cdot - Co(Phi_tilde(XY))||=%.3e  "
           "||cdot - Co(plain XY)||=%.4f (non-vacuity gap)\n", oracle, vac);
    AIC_CHECK_MSG(oracle < 1e-9,
                  "T8(a): cdot != Co_{P,R}(Phi_tilde(X*Y)) (%.3e) -- star wrong?",
                  oracle);
    AIC_CHECK_MSG(vac > 0.03,
                  "T8(a): non-vacuity gap ||cdot - Co(plain XY)||=%.3e <= 0.03 -- "
                  "fixture not exercising the star in cdot", vac);

    /* (b) bound ||X.Y - X*Y|| on this exactly-idempotent fixture: machine-zero. */
    acb_mat_t star, db;
    acb_mat_init(star, n, n); acb_mat_init(db, n, n);
    aic_ecstar_star(star, &ae.A, X, Y, prec);     /* X * Y (ambient star product) */
    acb_mat_sub(db, prod, star, prec);
    double bnd0 = mid_opnorm_d(db, prec);
    printf("T8(b) eta=0 ||X.Y - X*Y||=%.3e (machine-zero)\n", bnd0);
    AIC_CHECK_MSG(bnd0 < 1e-9, "T8(b): at eta=0 X.Y != X*Y (%.3e)", bnd0);

    acb_mat_clear(db); acb_mat_clear(star);
    acb_mat_clear(XYplain);
    acb_mat_clear(dvac); acb_mat_clear(dorc);
    acb_mat_clear(ZplainCo); acb_mat_clear(Zplain);
    acb_mat_clear(ZexpCo); acb_mat_clear(Zexp);
    acb_mat_clear(prod); acb_mat_clear(Y); acb_mat_clear(X);
    for (slong m = 0; m < dYR; m++) acb_mat_clear(CY[m]);
    if (CY) flint_free(CY);
    for (slong m = 0; m < dXY; m++) acb_mat_clear(CX[m]);
    if (CX) flint_free(CX);
    acb_mat_clear(CoPR); acb_mat_clear(CoQR); acb_mat_clear(CoPQ);
    acb_mat_clear(R); acb_mat_clear(Q); acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T8(c) UNIT (.tex:1082): Ptilde = Co_P(P) is the unit of (S_P,.). The unit law
 * ||Ptilde . X - X|| <= O(delta+eps)||X|| holds only when P is a genuine
 * delta-PROJECTION IN A (a star-idempotent close to A); the proj_residual probe
 * shows |v01><v01| is NOT in A for compress_idemp (residual 0.333, ||P*P-P||=0.5),
 * so its O(delta+eps) bound is vacuous (delta = O(1)) and it is a wrong unit
 * fixture. We test the unit law where P IS a genuine projection in A:
 *   - identity(n): A=M_n, star=plain product, every ambient projector P is
 *     star-idempotent => Ptilde = P EXACTLY and the unit law is machine-exact;
 *   - the OBLIQUE compress_idemp(4,2) with P=|e1><e1| (proj_residual in A = 0,
 *     verified) so the oblique star path is also exercised on a clean unit. */
static void t8_unit_on(const char *name, aic_assoc_ecstar *ae, const acb_mat_t P,
                       slong prec)
{
    slong n = ae->A.n, d = ae->A.dim_A;
    acb_mat_t CoPP, Ptil, lu, ru, dlu, dru;
    acb_mat_init(CoPP, d, d);
    aic_corner_Co(CoPP, &ae->A, P, P, prec);
    acb_mat_t *CP;
    slong dP;
    aic_corner_extract(&CP, &dP, CoPP, &ae->A, prec);
    AIC_CHECK_MSG(dP >= 1, "T8(c) %s: dim S_P=%ld < 1", name, (long) dP);
    acb_mat_init(Ptil, n, n);
    aic_corner_Ptilde(Ptil, &ae->A, P, prec);
    acb_mat_init(lu, n, n); acb_mat_init(ru, n, n);
    acb_mat_init(dlu, n, n); acb_mat_init(dru, n, n);
    double lmax = 0.0, rmax = 0.0;
    for (slong m = 0; m < dP; m++) {
        aic_corner_cdot(lu, &ae->A, CoPP, Ptil, CP[m], prec);   /* Ptilde . X */
        aic_corner_cdot(ru, &ae->A, CoPP, CP[m], Ptil, prec);   /* X . Ptilde */
        acb_mat_sub(dlu, lu, CP[m], prec);
        acb_mat_sub(dru, ru, CP[m], prec);
        double le = mid_opnorm_d(dlu, prec), re = mid_opnorm_d(dru, prec);
        if (le > lmax) lmax = le;
        if (re > rmax) rmax = re;
    }
    printf("T8(c) unit %-18s dim S_P=%ld max||Ptilde.X-X||=%.3e "
           "max||X.Ptilde-X||=%.3e\n", name, (long) dP, lmax, rmax);
    AIC_CHECK_MSG(lmax < 1e-9 && rmax < 1e-9,
                  "T8(c) %s: Ptilde not the unit of (S_P,.) (%.3e, %.3e)", name,
                  lmax, rmax);
    acb_mat_clear(dru); acb_mat_clear(dlu); acb_mat_clear(ru); acb_mat_clear(lu);
    acb_mat_clear(Ptil);
    for (slong m = 0; m < dP; m++) acb_mat_clear(CP[m]);
    if (CP) flint_free(CP);
    acb_mat_clear(CoPP);
}

static void test_t8_unit(void)
{
    const slong prec = 256;
    /* identity(3): exact oracle. P = diag(1,1,0) (rank-2 genuine projector). */
    {
        aic_ucp_kraus phi;
        make_identity(&phi, 3);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        acb_mat_t P;
        diag_proj(P, 3, 2);
        t8_unit_on("id(3) P=diag2", &ae, P, prec);
        acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    /* oblique compress_idemp(4,2): P=|e1><e1| (genuine member of A, residual 0). */
    {
        aic_ucp_kraus phi;
        make_compress_idemp(&phi, 4, 2);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        acb_mat_t P;
        diag_proj_at(P, 4, 1);
        t8_unit_on("compress(4,2) e1", &ae, P, prec);
        acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* T8(b) eta>0: ||X.Y - X*Y|| / eta is a bounded constant c on dep+conj(4). */
static void test_t8_bound(void)
{
    const slong prec = 256;
    printf("T8(b) ||X.Y - X*Y|| vs eta (dep+conj(4)):\n");
    double cmax = 0.0;
    const double ts[] = {0.02, 0.06};
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phi;
        make_eta_family(&phi, 4, ts[it], prec);
        double eta = eta_proxy(&phi, prec);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        slong d = ae.A.dim_A, n = 4;
        /* P=Q=R=diag(1,1,0,0) (rank 2). X,Y in S_P; use basis elements. */
        acb_mat_t P;
        diag_proj(P, n, 2);
        acb_mat_t CoPP;
        acb_mat_init(CoPP, d, d);
        aic_corner_Co(CoPP, &ae.A, P, P, prec);
        acb_mat_t *CP;
        slong dP;
        aic_corner_extract(&CP, &dP, CoPP, &ae.A, prec);
        AIC_CHECK_MSG(dP >= 2, "T8(b) eta>0: dim S_P=%ld < 2", (long) dP);
        acb_mat_t prod, star, db;
        acb_mat_init(prod, n, n); acb_mat_init(star, n, n); acb_mat_init(db, n, n);
        double worst = 0.0;
        for (slong a = 0; a < dP; a++)
            for (slong b = 0; b < dP; b++) {
                aic_corner_cdot(prod, &ae.A, CoPP, CP[a], CP[b], prec);
                aic_ecstar_star(star, &ae.A, CP[a], CP[b], prec);
                acb_mat_sub(db, prod, star, prec);
                double e = mid_opnorm_d(db, prec);
                if (e > worst) worst = e;
            }
        double c = eta > 1e-14 ? worst / eta : 0.0;
        if (c > cmax) cmax = c;
        printf("  t=%.2f eta=%.4e: max||X.Y - X*Y||=%.4e (c=%.3f)\n",
               ts[it], eta, worst, c);
        AIC_CHECK_MSG(c < 5.0,
                      "T8(b): ||X.Y-X*Y||/eta=%.3f exceeds 5 (t=%.2f) -- not "
                      "O(delta+eps)?", c, ts[it]);
        AIC_CHECK_MSG(worst > 1e-7,
                      "T8(b): X.Y == X*Y at eta>0 (t=%.2f) -- expected O(eta) gap",
                      ts[it]);
        acb_mat_clear(db); acb_mat_clear(star); acb_mat_clear(prod);
        for (slong m = 0; m < dP; m++) acb_mat_clear(CP[m]);
        if (CP) flint_free(CP);
        acb_mat_clear(CoPP); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("  -> measured c (||X.Y-X*Y||/eta) max = %.3f\n", cmax);
}

/* ===================== T9: lem_alpha block bijection (.tex:1086-1119) ===== *
 * alpha = sum_{jk} Co_{P,Q}|_{S_{Pj,Qk}} : (+)_{jk} S_{Pj,Qk} -> S_{P,Q} is a
 * linear bijection (tex:1096) with ||alpha|| <= pq + O(pq(d+e)), ||alpha^{-1}|| <=
 * 1 + O(pq(d+e)) (tex:1097-1099). Tests:
 *   (i)  eta=0 EXACT oracle (tex:1124): exact orthogonal sub-projectors in M_n,
 *        dim S_{P,Q} = sum dim S_{Pj,Qk} (the dim-count assert inside alpha is
 *        itself this oracle; we also re-assert it), round-trip alpha alpha^{-1} ~
 *        I and alpha^{-1} alpha ~ I (machine-zero), ||alpha|| <= pq, ||alpha^{-1}||
 *        ~ 1, and ||gamma|| ~ 0 reported.
 *   (ii) eta>0 (dep+conj(4)): the bounds ||alpha|| <= pq + O(pq(d+e)), ||alpha^{-1}||
 *        <= 1 + O(pq(d+e)) hold; report the constants and ||gamma||.
 *   (iii) FAIL-LOUD teeth: OVERLAPPING projectors (||P1 P2|| not small) violate the
 *        lem_alpha hypothesis pq(d+e)<const => ||gamma|| >= 1 => the assert fires
 *        (forked child, SIGABRT; mirrors test_idemp test_entry_guard).
 * MUTATION (the .tex TYPO, mutation-proven): building beta_{jk} as Co_{Pj,Qj}
 * (the tex:1109 typo, Q_j not Q_k) on the eta=0 id(4) p=q=2 fixture with distinct
 * Q_1=|e2><e2|, Q_2=|e3><e3| makes beta wrong => beta alpha != I and ||gamma||
 * jumps from 0 to 1.0000 (CONFIRMED), tripping the ||gamma||<1 assert (SIGABRT).
 * Restored. Documented in src/aic_corner_product.c. */

/* Build alpha on the identity(n) algebra (A=M_n, star=plain) with the given
 * sub-projector supports: P_j = |e_{Psup[j]}><e_{Psup[j]}|, Q_k similarly.
 * Asserts the dim count, round-trips, and the norm bounds. The supports must be
 * DISJOINT (so P_j P_l = delta_jl P_l exactly) for the eta=0 oracle. */
static void t9_identity(const char *name, slong n, const slong *Psup, slong p,
                        const slong *Qsup, slong q, slong prec)
{
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);

    acb_mat_t *Pp = flint_malloc((size_t) p * sizeof(acb_mat_t));
    acb_mat_t *Qp = flint_malloc((size_t) q * sizeof(acb_mat_t));
    acb_mat_t P, Q;
    acb_mat_init(P, n, n); acb_mat_zero(P);
    acb_mat_init(Q, n, n); acb_mat_zero(Q);
    for (slong j = 0; j < p; j++) {
        diag_proj_at(Pp[j], n, Psup[j]);
        acb_add(acb_mat_entry(P, Psup[j], Psup[j]),
                acb_mat_entry(P, Psup[j], Psup[j]),
                acb_mat_entry(Pp[j], Psup[j], Psup[j]), prec);
    }
    for (slong k = 0; k < q; k++) {
        diag_proj_at(Qp[k], n, Qsup[k]);
        acb_add(acb_mat_entry(Q, Qsup[k], Qsup[k]),
                acb_mat_entry(Q, Qsup[k], Qsup[k]),
                acb_mat_entry(Qp[k], Qsup[k], Qsup[k]), prec);
    }

    slong N, dPQ;
    aic_corner_alpha_dims(&N, &dPQ, &ae.A, (const acb_mat_t *) Pp, p, P,
                          (const acb_mat_t *) Qp, q, Q, prec);
    AIC_CHECK_MSG(N == dPQ,
                  "T9 %s: dim count N=%ld != dim S_{P,Q}=%ld", name, (long) N,
                  (long) dPQ);
    /* exact M_n value: dim S_{P,Q} = rank(P) rank(Q) = p*q (each part rank 1). */
    AIC_CHECK_MSG(dPQ == p * q,
                  "T9 %s: dim S_{P,Q}=%ld != p*q=%ld (M_n exact)", name,
                  (long) dPQ, (long) (p * q));

    acb_mat_t alpha, alpha_inv, prod1, prod2, one1, one2;
    acb_mat_init(alpha, dPQ, N);
    acb_mat_init(alpha_inv, N, dPQ);
    arb_t na, nai, ng, ngrad;
    arb_init(na); arb_init(nai); arb_init(ng); arb_init(ngrad);
    slong N2, dPQ2;
    aic_corner_alpha(alpha, alpha_inv, na, nai, ng, ngrad, &N2, &dPQ2, &ae.A,
                     (const acb_mat_t *) Pp, p, P, (const acb_mat_t *) Qp, q, Q,
                     prec);

    acb_mat_init(prod1, dPQ, dPQ);
    acb_mat_init(prod2, N, N);
    acb_mat_init(one1, dPQ, dPQ);
    acb_mat_init(one2, N, N);
    acb_mat_mul(prod1, alpha, alpha_inv, prec);   /* alpha alpha^{-1} = I_dPQ */
    acb_mat_mul(prod2, alpha_inv, alpha, prec);   /* alpha^{-1} alpha = I_N   */
    acb_mat_one(one1); acb_mat_one(one2);
    acb_mat_sub(prod1, prod1, one1, prec);
    acb_mat_sub(prod2, prod2, one2, prec);
    double rt1 = mid_opnorm_d(prod1, prec), rt2 = mid_opnorm_d(prod2, prec);
    printf("T9 %-16s p=%ld q=%ld dim=%ld ||a||=%.4f ||a^-1||=%.4f ||g||=%.3e "
           "(rad=%.2e) rt(aa^-1)=%.2e rt(a^-1a)=%.2e\n", name, (long) p, (long) q,
           (long) dPQ, dd(na), dd(nai), dd(ng), dd(ngrad), rt1, rt2);
    /* The round-trips rt1, rt2 confirm acb_mat_solve SUCCEEDED and alpha is
     * SURJECTIVE (so the certified inverse exists): they are NOT a witness that
     * alpha's/beta's ENTRIES are numerically correct. Because alpha_inv =
     * (beta alpha)^{-1} beta, alpha^{-1} alpha = (beta alpha)^{-1} (beta alpha) =
     * I_N is an ALGEBRAIC IDENTITY for ANY alpha, beta (a beta->1.5 beta or an
     * alpha-entry mutation leaves it I_N), and alpha alpha^{-1} = I_dPQ likewise
     * holds whenever alpha has full row rank dPQ. The genuine entry-level oracle
     * is t9_alpha_entry_oracle (alpha^dag alpha = I_N + the per-column D_l
     * reconstruction), computed WITHOUT going through alpha_inv. */
    AIC_CHECK_MSG(rt1 < 1e-9 && rt2 < 1e-9,
                  "T9 %s: alpha solve failed / alpha not surjective (round-trip "
                  "%.3e, %.3e)", name, rt1, rt2);
    /* eta=0 bounds: ||alpha|| <= pq, ||alpha^{-1}|| ~ 1 (exact orthogonal). */
    AIC_CHECK_MSG(dd(na) <= (double) (p * q) + 1e-6,
                  "T9 %s: ||alpha||=%.4f > pq=%ld", name, dd(na), (long) (p * q));
    AIC_CHECK_MSG(dd(nai) <= 1.0 + 1e-6,
                  "T9 %s: ||alpha^{-1}||=%.4f > 1 (eta=0)", name, dd(nai));
    AIC_CHECK_MSG(dd(ng) < 1e-6,
                  "T9 %s: ||gamma||=%.3e not ~0 at eta=0", name, dd(ng));

    arb_clear(ngrad); arb_clear(ng); arb_clear(nai); arb_clear(na);
    acb_mat_clear(one2); acb_mat_clear(one1);
    acb_mat_clear(prod2); acb_mat_clear(prod1);
    acb_mat_clear(alpha_inv); acb_mat_clear(alpha);
    acb_mat_clear(Q); acb_mat_clear(P);
    for (slong k = 0; k < q; k++) acb_mat_clear(Qp[k]);
    for (slong j = 0; j < p; j++) acb_mat_clear(Pp[j]);
    flint_free(Qp); flint_free(Pp);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T9(ii) eta>0: bounds hold on dep+conj(4); report constants. The dep+conj(4)
 * algebra A is diagonal-dominant, so the only NONEMPTY corners S_{Pj,Qk} are the
 * diagonal ones (off-diagonal P{0,1} x Q{2,3} corners are 0-dim, measured). We
 * take P = Q = diag(1,1,0,0) with parts P_1=Q_1=|e0><e0|, P_2=Q_2=|e1><e1|
 * (disjoint rank-1 delta-projections in A): N = dim S_{P,Q} = 2 (the two diagonal
 * blocks). p=q=2 still exercises the four-block alpha (two blocks empty). */
static void t9_eta_pos(slong prec)
{
    aic_ucp_kraus phi;
    make_eta_family(&phi, 4, 0.04, prec);
    double eta = eta_proxy(&phi, prec);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong n = 4, p = 2, q = 2;
    acb_mat_t Pp[2], Qp[2], P, Q;
    diag_proj_at(Pp[0], n, 0); diag_proj_at(Pp[1], n, 1);
    diag_proj_at(Qp[0], n, 0); diag_proj_at(Qp[1], n, 1);
    diag_proj(P, n, 2);                         /* P = diag(1,1,0,0) */
    diag_proj(Q, n, 2);                         /* Q = diag(1,1,0,0) = P */

    slong N, dPQ;
    aic_corner_alpha_dims(&N, &dPQ, &ae.A, (const acb_mat_t *) Pp, p, P,
                          (const acb_mat_t *) Qp, q, Q, prec);
    AIC_CHECK_MSG(N == dPQ, "T9 eta>0: N=%ld != dPQ=%ld", (long) N, (long) dPQ);
    acb_mat_t alpha, alpha_inv;
    acb_mat_init(alpha, dPQ, N);
    acb_mat_init(alpha_inv, N, dPQ);
    arb_t na, nai, ng, ngrad;
    arb_init(na); arb_init(nai); arb_init(ng); arb_init(ngrad);
    slong N2, dPQ2;
    aic_corner_alpha(alpha, alpha_inv, na, nai, ng, ngrad, &N2, &dPQ2, &ae.A,
                     (const acb_mat_t *) Pp, p, P, (const acb_mat_t *) Qp, q, Q,
                     prec);
    /* The two surviving (diagonal) corners give ||alpha|| ~ ||alpha^{-1}|| ~ 1
     * (isometric at eta=0), so the genuine O(delta+eps) deviation is from 1; the
     * paper's pq-bound is a (loose-here) upper bound, asserted below. The reported
     * constants are |dev|/eta. */
    double ca = (dd(na) - 1.0) / eta;       /* (||a|| - 1)/eta   */
    double cai = (dd(nai) - 1.0) / eta;     /* (||a^-1|| - 1)/eta */
    double cg = dd(ng) / eta;               /* ||gamma||/eta      */
    printf("T9 eta>0 dep+conj(4): eta=%.4e ||a||=%.4f (dev/eta=%.2f) "
           "||a^-1||=%.4f (dev/eta=%.2f) ||g||=%.3e (/eta=%.2f) rad=%.2e [pq=%ld]\n",
           eta, dd(na), ca, dd(nai), cai, dd(ng), cg, dd(ngrad), (long) (p * q));
    AIC_CHECK_MSG(dd(na) <= (double) (p * q) + 5.0 * eta + 1e-9,
                  "T9 eta>0: ||alpha||=%.4f exceeds pq + 5 eta", dd(na));
    AIC_CHECK_MSG(dd(nai) <= 1.0 + 5.0 * eta + 1e-9,
                  "T9 eta>0: ||alpha^{-1}||=%.4f exceeds 1 + 5 eta", dd(nai));
    AIC_CHECK_MSG(dd(ng) < 1.0,
                  "T9 eta>0: ||gamma||=%.4f >= 1 (contraction lost)", dd(ng));
    arb_clear(ngrad); arb_clear(ng); arb_clear(nai); arb_clear(na);
    acb_mat_clear(alpha_inv); acb_mat_clear(alpha);
    acb_mat_clear(Q); acb_mat_clear(P);
    acb_mat_clear(Qp[1]); acb_mat_clear(Qp[0]);
    acb_mat_clear(Pp[1]); acb_mat_clear(Pp[0]);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T9(iii) FAIL-LOUD teeth: OVERLAPPING P_1,P_2 (P_1=|e0><e0|, P_2=|v><v|,
 * v=(e0+e1)/sqrt2, so ||P_1 P_2|| = 1/sqrt2, NOT small) violate the lem_alpha
 * hypothesis ||P_j P_l - delta_jl P_l|| <= delta. The construction must REJECT
 * this ill-posed input via SIGABRT. For GROSS overlap the EARLIEST hypothesis
 * check that fires is the prop_P basin inside aic_corner_Co: P = P_1 + P_2 with a
 * large overlap is far from a delta-projection, so the symmetric combination M =
 * 1/2(L_P R_P + R_P L_P) has ||M^2 - M|| >= 1/4 and aic_prop_P aborts.
 *
 * Finding 6 (corrected). prop_P gates only GROSS overlap; it is NOT the
 * load-bearing invertibility guard. The review found overlaps in [0.05, 0.20]
 * PASS prop_P and are CORRECTLY handled (not rejected) by the certified ||gamma||<1
 * guard inside aic_corner_alpha -- that guard is the load-bearing invertibility
 * gate. The lemma's named per-projector hypotheses (||P_j P_l - delta P_l|| <=
 * delta, etc., tex:1091) are SUFFICIENT CONDITIONS for ||gamma||<1 (tex:1114), not
 * separately asserted in code. The ||gamma||<1 assert is exercised on the VALID
 * path (||gamma|| reported in every T9 pass: 0 at eta=0, 7.3e-6 at eta=0.028) and
 * fires on the rare borderline input that slips past prop_P. Either way an
 * out-of-hypothesis input is fail-loud (Rule 4) -- that is the tooth. This test
 * uses ||P1 P2|| = 0.707 (gross), so prop_P fires here. Forked child (mirrors
 * test_idemp). */
static void test_t9_failloud(void)
{
    const slong prec = 256, n = 4;
    pid_t pid = fork();
    if (pid == 0) {
        aic_ucp_kraus phi;
        make_identity(&phi, n);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        acb_mat_t Pp[2], Qp[1], P, Q;
        diag_proj_at(Pp[0], n, 0);                    /* |e0><e0| */
        double v01[4] = {0.70710678118654752, 0.70710678118654752, 0, 0};
        rank1_proj(Pp[1], n, v01);                    /* |v><v|, overlaps P_0 */
        diag_proj_at(Qp[0], n, 2);
        /* P = P_0 + P_1 is NOT a clean delta-projection (overlap); build it anyway
         * (the construction's fail-loud is what we are testing). */
        acb_mat_init(P, n, n);
        acb_mat_add(P, Pp[0], Pp[1], prec);
        acb_mat_init(Q, n, n);
        acb_set_si(acb_mat_entry(Q, 2, 2), 1);
        slong N, dPQ;
        /* aic_corner_alpha_dims itself builds Co_{P,P}, so the prop_P basin assert
         * fires HERE (the earliest hypothesis check) on this overlapping P. */
        aic_corner_alpha_dims(&N, &dPQ, &ae.A, (const acb_mat_t *) Pp, 2, P,
                              (const acb_mat_t *) Qp, 1, Q, prec);
        acb_mat_t alpha;
        acb_mat_init(alpha, dPQ, N);
        arb_t ng;
        arb_init(ng);
        slong N2, dPQ2;
        aic_corner_alpha(alpha, NULL, NULL, NULL, ng, NULL, &N2, &dPQ2, &ae.A,
                         (const acb_mat_t *) Pp, 2, P, (const acb_mat_t *) Qp, 1,
                         Q, prec);    /* (or here: ||gamma||>=1 / dim mismatch) */
        _exit(0);                     /* reached only if NOT aborted */
    }
    int status = 0;
    waitpid(pid, &status, 0);
    AIC_CHECK_MSG(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT,
                  "T9(iii): overlapping projectors did not abort via SIGABRT "
                  "(WIFSIGNALED=%d, WTERMSIG=%d)", WIFSIGNALED(status),
                  WIFSIGNALED(status) ? WTERMSIG(status) : -1);
    printf("T9(iii) fail-loud: overlapping P_1,P_2 (||P1 P2||=0.707, GROSS) "
           "rejected via SIGABRT (prop_P basin gates gross overlap; ||gamma||<1 is "
           "the load-bearing invertibility guard for borderline inputs)\n");
}

/* T9(iv) INDEPENDENT eta=0 ENTRY oracle (FIX 1, bead aic — Increment-2a hostile
 * review). The t9_identity round-trips alpha alpha^{-1} ~ I cannot fail on alpha's
 * ENTRIES: with alpha_inv = (beta alpha)^{-1} beta, alpha^{-1} alpha = I_N is an
 * algebraic identity for ANY alpha, beta (a beta->1.5 beta mutation stays green; the
 * review proved only the ||gamma|| check caught it). This oracle tests alpha's
 * entries DIRECTLY, with NO path through alpha_inv / acb_mat_solve:
 *
 *   At eta=0 on the identity(4) algebra (A = M_4, star = plain product) with
 *   DISJOINT rank-1 parts P_j=|e_j><e_j| (j=0,1), Q_k=|e_{2+k}><e_{2+k}| (k=0,1),
 *   almost-containment is EXACT: each block S_{Pj,Qk} = span{|e_j><e_{2+k}|} (a
 *   matrix unit, dim 1, all 4 blocks populated -- the off-diagonal blocks are
 *   present), and Co_{P,Q}(C^{jk}_m) = P C^{jk}_m Q = C^{jk}_m EXACTLY. So
 *   alpha is the rectangular CHANGE-OF-BASIS ISOMETRY: column (b,m) of alpha is
 *   <D_l, C^{jk}_m>_F (re-coordinating the block basis element in {D_l}). Both
 *   {C^{jk}_m} (across blocks) and {D_l} are Frobenius-orthonormal bases of the
 *   SAME space S_{P,Q} (N = dPQ = 4), so alpha is UNITARY. We check:
 *     (a) ORTHONORMAL COLUMNS: ||alpha^dag alpha - I_N|| < tol. A column-scale or
 *         transpose bug in the alpha ASSEMBLY breaks this, independently of
 *         alpha_inv (the round-trip is blind to it).
 *     (b) RECONSTRUCTION: for each block column (b,m), sum_l alpha[l,off+m] D_l ==
 *         C^{jk}_m to machine-zero -- an entry-level witness that the block basis
 *         element is preserved, recomputed via an INDEPENDENT path ({D_l},
 *         {C^{jk}_m} re-extracted here, NOT the solve).
 * MUTATION (forked children, confirmed RED). At eta=0 alpha is a 4x4 PERMUTATION
 *   matrix (each block basis element = a matrix unit |e_j><e_{2+k}| maps to a single
 *   {D_l}), so two complementary mutations isolate the two checks:
 *     col_scale=1.5 on alpha column 0 -> (a) ||alpha^dag alpha - I|| = 1.5^2-1 =
 *       1.250e0 RED, (b) reconstruction residual 5.000e-1 RED (SIGABRT). Both fire.
 *     transpose_mut=1 swaps alpha COLUMNS 0,1 -> (a) is BLIND (permutation stays
 *       orthonormal, ||a^dag a - I||=0) but (b) RED: block 0 rebuilds C^{01} != C^{00},
 *       reconstruction residual 1.414e0 (SIGABRT). This proves (b) has teeth against
 *       a column permutation, the kind of bug (a) cannot see.
 *   With col_scale=1.0, transpose_mut=0 (no mutation) the oracle is GREEN (both
 *   residuals ~1e-74). See t9_entry_check_aborts. */
static void t9_entry_oracle(double col_scale, int transpose_mut, slong prec)
{
    const slong n = 4, p = 2, q = 2;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;

    /* disjoint rank-1 parts; P=sum Pj, Q=sum Qk. */
    acb_mat_t Pp[2], Qp[2], P, Q;
    diag_proj_at(Pp[0], n, 0); diag_proj_at(Pp[1], n, 1);
    diag_proj_at(Qp[0], n, 2); diag_proj_at(Qp[1], n, 3);
    diag_proj(P, n, 2);                          /* diag(1,1,0,0) */
    acb_mat_init(Q, n, n); acb_mat_zero(Q);
    acb_set_si(acb_mat_entry(Q, 2, 2), 1);
    acb_set_si(acb_mat_entry(Q, 3, 3), 1);       /* diag(0,0,1,1) */

    slong N, dPQ;
    aic_corner_alpha_dims(&N, &dPQ, &ae.A, (const acb_mat_t *) Pp, p, P,
                          (const acb_mat_t *) Qp, q, Q, prec);
    AIC_CHECK_MSG(N == dPQ && N == 4,
                  "T9(iv): expected N=dPQ=4, got N=%ld dPQ=%ld", (long) N,
                  (long) dPQ);

    acb_mat_t alpha;
    acb_mat_init(alpha, dPQ, N);
    arb_t na, nai, ng, ngrad;
    arb_init(na); arb_init(nai); arb_init(ng); arb_init(ngrad);
    slong N2, dPQ2;
    aic_corner_alpha(alpha, NULL, na, nai, ng, ngrad, &N2, &dPQ2, &ae.A,
                     (const acb_mat_t *) Pp, p, P, (const acb_mat_t *) Qp, q, Q,
                     prec);

    /* MUTATIONS (applied to the test's COPY of alpha; production is untouched). */
    if (col_scale != 1.0) {
        arb_t s;
        arb_init(s);
        arb_set_d(s, col_scale);
        for (slong l = 0; l < dPQ; l++)              /* scale alpha column 0 */
            acb_mul_arb(acb_mat_entry(alpha, l, 0), acb_mat_entry(alpha, l, 0),
                        s, prec);
        arb_clear(s);
    }
    if (transpose_mut) {
        /* swap alpha columns 0 and 1 (the block-0 and block-1 coordinate columns).
         * alpha here is a permutation matrix, so a COLUMN swap is invisible to (a)
         * (||a^dag a - I|| stays 0) but corrupts (b): block 0's reconstruction now
         * uses column 1's coords and rebuilds C^{01} != C^{00}. This is the
         * mutation that exercises the RECONSTRUCTION check specifically. */
        acb_t tmp;
        acb_init(tmp);
        for (slong l = 0; l < dPQ; l++) {
            acb_set(tmp, acb_mat_entry(alpha, l, 0));
            acb_set(acb_mat_entry(alpha, l, 0), acb_mat_entry(alpha, l, 1));
            acb_set(acb_mat_entry(alpha, l, 1), tmp);
        }
        acb_clear(tmp);
    }

    /* INDEPENDENT bases: D_l of S_{P,Q}, and C^{jk}_m per block (re-extracted, NOT
     * via alpha_inv). Mirrors aic_corner_build_blocks but kept local for
     * independence from the production builder. */
    acb_mat_t CoPQ;
    acb_mat_init(CoPQ, d, d);
    aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
    acb_mat_t *Dl;
    slong dD;
    aic_corner_extract(&Dl, &dD, CoPQ, &ae.A, prec);
    AIC_CHECK_MSG(dD == dPQ, "T9(iv): dim {D_l}=%ld != dPQ=%ld", (long) dD,
                  (long) dPQ);

    /* (a) orthonormal columns: ||alpha^dag alpha - I_N||_op. */
    acb_mat_t ad, gram, eye, gdiff;
    acb_mat_init(ad, N, dPQ);
    acb_mat_init(gram, N, N);
    acb_mat_init(eye, N, N);
    acb_mat_init(gdiff, N, N);
    acb_mat_conjugate_transpose(ad, alpha);
    acb_mat_mul(gram, ad, alpha, prec);          /* alpha^dag alpha */
    acb_mat_one(eye);
    acb_mat_sub(gdiff, gram, eye, prec);
    double orth = mid_opnorm_d(gdiff, prec);

    /* (b) reconstruction: for each block column (b,m), R = sum_l alpha[l,off+m] D_l
     * must equal C^{jk}_m. */
    double recon = 0.0;
    slong off = 0;
    for (slong j = 0; j < p; j++)
        for (slong k = 0; k < q; k++) {
            acb_mat_t CoBlk;
            acb_mat_init(CoBlk, d, d);
            aic_corner_Co(CoBlk, &ae.A, Pp[j], Qp[k], prec);
            acb_mat_t *Cm;
            slong dB;
            aic_corner_extract(&Cm, &dB, CoBlk, &ae.A, prec);
            for (slong m = 0; m < dB; m++) {
                acb_mat_t R, dr;
                acb_mat_init(R, n, n);
                acb_mat_init(dr, n, n);
                acb_mat_zero(R);
                for (slong l = 0; l < dPQ; l++) {  /* R = sum_l alpha[l,off+m] D_l */
                    acb_mat_t sc;
                    acb_mat_init(sc, n, n);
                    acb_mat_scalar_mul_acb(sc, Dl[l],
                                           acb_mat_entry(alpha, l, off + m), prec);
                    acb_mat_add(R, R, sc, prec);
                    acb_mat_clear(sc);
                }
                acb_mat_sub(dr, R, Cm[m], prec);   /* == 0 (block elt preserved) */
                double e = mid_opnorm_d(dr, prec);
                if (e > recon) recon = e;
                acb_mat_clear(dr); acb_mat_clear(R);
            }
            for (slong m = 0; m < dB; m++) acb_mat_clear(Cm[m]);
            if (Cm) flint_free(Cm);
            acb_mat_clear(CoBlk);
            off += dB;
        }

    printf("T9(iv) entry oracle [scale=%.2f tr=%d]: ||a^dag a - I||=%.3e "
           "max||sum_l a[l,*]D_l - C^{jk}||=%.3e\n", col_scale, transpose_mut,
           orth, recon);
    AIC_CHECK_MSG(orth < 1e-9,
                  "T9(iv): alpha columns not orthonormal (||a^dag a - I||=%.3e) -- "
                  "alpha-entry scale/transpose bug", orth);
    AIC_CHECK_MSG(recon < 1e-9,
                  "T9(iv): block basis element NOT reconstructed (max %.3e) -- "
                  "alpha entries wrong (independent of the solve)", recon);

    acb_mat_clear(gdiff); acb_mat_clear(eye); acb_mat_clear(gram); acb_mat_clear(ad);
    for (slong l = 0; l < dD; l++) acb_mat_clear(Dl[l]);
    if (Dl) flint_free(Dl);
    acb_mat_clear(CoPQ);
    arb_clear(ngrad); arb_clear(ng); arb_clear(nai); arb_clear(na);
    acb_mat_clear(alpha);
    acb_mat_clear(Q); acb_mat_clear(P);
    acb_mat_clear(Qp[1]); acb_mat_clear(Qp[0]);
    acb_mat_clear(Pp[1]); acb_mat_clear(Pp[0]);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T9(iv) MUTATION-PROOF: fork children that run t9_entry_oracle with a mutation
 * and confirm SIGABRT (one of (a)/(b) goes RED). Mirrors test_t9_failloud. Also
 * re-confirms the beta->1.5 beta mutation (which the round-trip was BLIND to) is
 * now caught at the ENTRY level: a beta scale changes alpha (its columns are
 * <D_l, Co_{P,Q}(C)>_F, computed from beta's blocks... ) -- but beta does NOT feed
 * alpha's columns, only the solve. So the relevant entry mutation is on alpha
 * itself (col_scale / transpose). The beta->1.5 beta soundness is instead covered
 * by the ||gamma|| guard (FIX 2 certified) + the round-trip's surjectivity role;
 * the entry oracle's job is to catch alpha-entry corruption the round-trip misses. */
static void t9_entry_check_aborts(double col_scale, int transpose_mut,
                                  const char *what)
{
    pid_t pid = fork();
    if (pid == 0) {
        t9_entry_oracle(col_scale, transpose_mut, 256);
        _exit(0);                /* reached only if NOT aborted */
    }
    int status = 0;
    waitpid(pid, &status, 0);
    AIC_CHECK_MSG(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT,
                  "T9(iv) mutation %s did not abort via SIGABRT (WIFSIGNALED=%d, "
                  "WTERMSIG=%d)", what, WIFSIGNALED(status),
                  WIFSIGNALED(status) ? WTERMSIG(status) : -1);
    printf("T9(iv) mutation %s rejected via SIGABRT (entry oracle has teeth)\n",
           what);
}

static void test_t9(void)
{
    const slong prec = 256;
    /* (i) eta=0 exact oracle. p=q=2: P=diag(1,1,0,0), Q=diag(0,0,1,1). */
    slong Psup2[2] = {0, 1}, Qsup2[2] = {2, 3};
    t9_identity("id(4) p=q=2", 4, Psup2, 2, Qsup2, 2, prec);
    /* p=2, q=1: Q a single rank-1 block. */
    slong Qsup1[1] = {2};
    t9_identity("id(4) p=2,q=1", 4, Psup2, 2, Qsup1, 1, prec);
    /* (ii) eta>0 bounds. */
    t9_eta_pos(prec);
    /* (iii) fail-loud teeth. */
    test_t9_failloud();
    /* (iv) INDEPENDENT eta=0 entry oracle (FIX 1) + its mutation-proof. */
    t9_entry_oracle(1.0, 0, prec);                 /* no mutation: GREEN */
    t9_entry_check_aborts(1.5, 0, "col0 scale 1.5 (breaks orthonormality (a))");
    t9_entry_check_aborts(1.0, 1, "swap cols 0,1 (breaks reconstruction (b))");
}

/* ===================== Increment 2b: §7 one-dimensional-Q machinery ======== *
 * T10 lem_PQ_Hilb inner product (.tex:1123-1132): Hermitian symmetry, linearity,
 *     |<X,X>-||X||^2| <= c(delta+eps)||X||^2 (report c), eta=0 GNS oracle, and the
 *     OBLIQUE star-vs-plain non-vacuity gap (the star is load-bearing).
 * T11 Ha-map (.tex:1146-1160): Ha_dag EXACT identity (with non-vacuity + a closed-
 *     form eta=0 cross-check so the identity has TEETH), ||Ha(Z)(X)-Z.X||_Euc bound
 *     (report c), and Ha^Q_{P,P} the O(delta+eps)-HOMOMORPHISM (the downstream-
 *     critical property; report c) with ||Ha(Ptilde)-1|| small.
 * T12 lem_PQR (.tex:1162-1177): | ||X.Y||-||X||||Y|| | <= c(delta+eps) (report c),
 *     machine-zero at eta=0.
 * T13 lem_1d_proj + equivalence (.tex:1179-1187): dim S_{P,Q} <= 1 for 1d P,Q;
 *     reflexive, symmetric, transitive.
 *
 * MUTATION-PROOFS (manual source perturbations, RED values documented; mirrors the
 * T8(a) star mutation discipline since these are NUMERICAL RED, not SIGABRT):
 *   - ip_1d star-drop (aic_ecstar_star -> acb_mat_mul in aic_corner_ip_1d): T10's
 *     oblique oracle goes RED at the non-vacuity gap (see test_t10 result line).
 *   - Ha cdot wrong-Co (CoQR -> CoPQ in the Y^dag.Z step of aic_corner_ha): T11's
 *     eta=0 Ha=Z. oracle goes RED (CONFIRMED, see test_t11 comment).
 *   - Ha Gram with the FROBENIUS ip instead of lem_PQ_Hilb ip: documented blind/RED
 *     verdict in test_t11 (the two ip's differ only at eta>0).
 *   All three confirmed below and restored. */

/* apply the Ha matrix (d_PQ x d_RQ in {C^{PQ}},{C^{RQ}} coords) to the m-th basis
 * element of S_{R,Q}: out = sum_l Ha[l,m] C^{PQ}_l. (The m-th input is C^{RQ}_m,
 * whose coord vector is the unit vector e_m, so the output is column m of Ha
 * re-expanded in the S_{P,Q} basis.) `out` init'd n x n. */
static void ha_apply_col(acb_mat_t out, const acb_mat_t Ha, acb_mat_t *CPQ,
                         slong dPQ, slong m, slong n, slong prec)
{
    acb_mat_zero(out);
    acb_t t;
    acb_init(t);
    for (slong l = 0; l < dPQ; l++)
        for (slong a = 0; a < n; a++)
            for (slong b = 0; b < n; b++) {
                acb_mul(t, acb_mat_entry(Ha, l, m), acb_mat_entry(CPQ[l], a, b), prec);
                acb_add(acb_mat_entry(out, a, b), acb_mat_entry(out, a, b), t, prec);
            }
    acb_clear(t);
}

/* T10 lem_PQ_Hilb inner product. */
static void test_t10(void)
{
    const slong prec = 256;
    /* (i) eta=0 GNS oracle on identity(n): A=M_n, star=plain, 1d Q=|e0><e0|,
     * P=diag(1..1) rank r => S_{P,Q} = {rows<r, col 0}, basis C_l=|e_l><e0|, so
     * <C_l,C_m>=delta_lm and <X,X>=||X||^2 EXACTLY (GNS). */
    {
        slong n = 4, r = 3;
        aic_ucp_kraus phi;
        make_identity(&phi, n);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        slong d = ae.A.dim_A;
        acb_mat_t P, Q;
        diag_proj(P, n, r);
        diag_proj_at(Q, n, 0);
        AIC_CHECK_MSG(aic_corner_dim_S(&ae.A, Q, Q, prec) == 1,
                      "T10: Q=|e0><e0| not 1d in M_%ld", (long) n);
        acb_mat_t CoPQ, CoQ, Qtilde;
        acb_mat_init(CoPQ, d, d); acb_mat_init(CoQ, d, d); acb_mat_init(Qtilde, n, n);
        aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
        aic_corner_Co(CoQ, &ae.A, Q, Q, prec);
        aic_corner_Ptilde(Qtilde, &ae.A, Q, prec);
        acb_mat_t *C;
        slong dPQ;
        aic_corner_extract(&C, &dPQ, CoPQ, &ae.A, prec);
        AIC_CHECK_MSG(dPQ == r, "T10: dim S_{P,Q}=%ld != r=%ld", (long) dPQ, (long) r);

        double herm = 0.0, gns = 0.0, lin = 0.0;
        acb_t ip, ipc, ipxx, ipyy, ipsum, ipsumdir, tmp;
        acb_init(ip); acb_init(ipc); acb_init(ipxx); acb_init(ipyy);
        acb_init(ipsum); acb_init(ipsumdir); acb_init(tmp);
        for (slong i = 0; i < dPQ; i++)
            for (slong j = 0; j < dPQ; j++) {
                aic_corner_ip_1d(ip, &ae.A, Qtilde, CoQ, C[j], C[i], prec); /* <C_i,C_j> */
                aic_corner_ip_1d(ipc, &ae.A, Qtilde, CoQ, C[i], C[j], prec);/* <C_j,C_i> */
                acb_conj(ipc, ipc);   /* conj<C_j,C_i> should == <C_i,C_j> */
                acb_sub(tmp, ip, ipc, prec);
                double h = abs_acb(tmp);
                if (h > herm) herm = h;
                /* GNS: <C_i,C_j> == delta_ij (orthonormal C_l). */
                if (i == j) acb_sub_si(ip, ip, 1, prec);
                double g = abs_acb(ip);
                if (g > gns) gns = g;
            }
        /* linearity in 2nd arg: <C0, C1 + 2 C2> == <C0,C1> + 2<C0,C2>. */
        if (dPQ >= 3) {
            acb_mat_t comb;
            acb_mat_init(comb, n, n);
            acb_mat_set(comb, C[1]);
            for (slong a = 0; a < n; a++)
                for (slong b = 0; b < n; b++) {
                    acb_mul_si(tmp, acb_mat_entry(C[2], a, b), 2, prec);
                    acb_add(acb_mat_entry(comb, a, b), acb_mat_entry(comb, a, b), tmp, prec);
                }
            aic_corner_ip_1d(ipsumdir, &ae.A, Qtilde, CoQ, comb, C[0], prec);
            aic_corner_ip_1d(ipxx, &ae.A, Qtilde, CoQ, C[1], C[0], prec);
            aic_corner_ip_1d(ipyy, &ae.A, Qtilde, CoQ, C[2], C[0], prec);
            acb_mul_si(ipyy, ipyy, 2, prec);
            acb_add(ipsum, ipxx, ipyy, prec);
            acb_sub(tmp, ipsum, ipsumdir, prec);
            lin = abs_acb(tmp);
            acb_mat_clear(comb);
        }
        /* CONJUGATION ORDER (FIX 2, bead aic). ip_1d computes <Y,X> = scalar
         * Co_Q(Y^dag * X): conjugate-linear in Y (the FIRST API slot), linear in X
         * (the SECOND). The herm/lin/gns checks above are all BLIND to the slot of
         * the conjugation -- swapping the impl to Co_Q(X^dag * Y) keeps them green
         * (every fixture uses real-symmetric corner elements where Y^dag X = X^dag Y
         * up to the real GNS value). We pin the order with a GENUINELY COMPLEX Y:
         *   <i*C0, C0> = conj(i) <C0,C0> = -i <C0,C0>   (conjugate-linear in Y),
         *   <C0, i*C0> =      i  <C0,C0> = +i <C0,C0>   (linear in X).
         * At eta=0 GNS <C0,C0> = 1 (real), so the two land at -i and +i; the X^dag*Y
         * swap flips both signs. MUTATION-PROVEN (source perturbation in
         * aic_corner_ip_1d: dagger X instead of Y and star Xd*Y, i.e. compute
         * Co_Q(X^dag*Y)): conj_err jumps from 0 to 2.000 (|-i - (+i)| = 2), RED,
         * while herm/gns/lin stay green (blind to the swap). Restored. */
        double conj_err = 0.0;
        {
            acb_mat_t C0i;
            acb_mat_init(C0i, n, n);
            acb_t iu, w;
            acb_init(iu);
            acb_init(w);
            acb_onei(iu);                              /* i */
            acb_mat_scalar_mul_acb(C0i, C[0], iu, prec); /* i*C0 */
            acb_t ip00, lhs, rhs, df;
            acb_init(ip00);
            acb_init(lhs);
            acb_init(rhs);
            acb_init(df);
            aic_corner_ip_1d(ip00, &ae.A, Qtilde, CoQ, C[0], C[0], prec); /* <C0,C0> */
            /* (1) conjugate-linear in Y (first slot): <i*C0, C0> == -i <C0,C0>. */
            aic_corner_ip_1d(lhs, &ae.A, Qtilde, CoQ, C[0], C0i, prec);   /* Y=i*C0,X=C0 */
            acb_neg(w, iu);                            /* -i */
            acb_mul(rhs, w, ip00, prec);               /* -i <C0,C0> */
            acb_sub(df, lhs, rhs, prec);
            double e1 = abs_acb(df);
            if (e1 > conj_err) conj_err = e1;
            /* (2) linear in X (second slot): <C0, i*C0> == +i <C0,C0>. */
            aic_corner_ip_1d(lhs, &ae.A, Qtilde, CoQ, C0i, C[0], prec);   /* X=i*C0,Y=C0 */
            acb_mul(rhs, iu, ip00, prec);              /* +i <C0,C0> */
            acb_sub(df, lhs, rhs, prec);
            double e2 = abs_acb(df);
            if (e2 > conj_err) conj_err = e2;
            acb_clear(df);
            acb_clear(rhs);
            acb_clear(lhs);
            acb_clear(ip00);
            acb_clear(w);
            acb_clear(iu);
            acb_mat_clear(C0i);
        }
        printf("T10(i) eta=0 GNS oracle (M_%ld, r=%ld): herm=%.2e |<C_i,C_j>-d_ij|"
               "=%.2e lin=%.2e conj_order=%.2e\n", (long) n, (long) r, herm, gns,
               lin, conj_err);
        AIC_CHECK_MSG(herm < 1e-9, "T10: <Y,X> not Hermitian-symmetric (%.3e)", herm);
        AIC_CHECK_MSG(gns < 1e-9, "T10: eta=0 not GNS (<C_i,C_j> != delta_ij, %.3e)", gns);
        AIC_CHECK_MSG(lin < 1e-9, "T10: <.,.> not linear in 2nd arg (%.3e)", lin);
        AIC_CHECK_MSG(conj_err < 1e-9,
                      "T10: ip_1d conjugation order wrong (conj-linear-in-Y, "
                      "linear-in-X violated, %.3e) -- Co_Q(X^dag*Y) instead of "
                      "Co_Q(Y^dag*X)?", conj_err);

        acb_clear(tmp); acb_clear(ipsumdir); acb_clear(ipsum);
        acb_clear(ipyy); acb_clear(ipxx); acb_clear(ipc); acb_clear(ip);
        for (slong m = 0; m < dPQ; m++) acb_mat_clear(C[m]);
        if (C) flint_free(C);
        acb_mat_clear(Qtilde); acb_mat_clear(CoQ); acb_mat_clear(CoPQ);
        acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }

    /* (ii) OBLIQUE non-vacuity (compress_idemp(4,2), exactly idempotent => clean).
     * The fixture must use GENUINE delta-projections IN A (the lem_PQ_Hilb
     * hypothesis); a diagonal projector that is NOT in A (e.g. |e0><e0|, proj
     * residual 0.667, ||Q*Q-Q|| = 1.0) gives a VACUOUS O(1) bound -- the same
     * trap as T8(c). The in-A 1d delta-projections of compress_idemp(4,2) are the
     * algebra elements B + B[0,0](1-Pi_M) for B a projector in B(M); we take
     *   P = |e1><e1|              (in A: B=|e1><e1|, B[0,0]=0; proj residual 0),
     *   Q = diag(1,0,1,1)         (in A: B=|e0><e0|, B[0,0]=1, the Mperp scalar;
     *                              a genuine 1d delta-projection, dim S_Q = 1).
     * Q touches the e2,e3 indices where Phi_tilde's OBLIQUE e2,e3 -> e0 coupling
     * lives, so the plain product Y^dag X leaves A and Co_Q(star) != Co_Q(plain):
     * measured star-vs-plain gap 0.6667 (the star is LOAD-BEARING in the inner
     * product). At eta=0 with genuine in-A projections |<X,X>-||X||^2| = 0 exactly.
     * MUTATION-PROVEN: replacing aic_ecstar_star with acb_mat_mul in
     * aic_corner_ip_1d makes ip_1d == the plain surrogate, so the gap is exactly
     * the RED amount of the (vac > 0.05) oracle. */
    {
        slong n = 4;
        aic_ucp_kraus phi;
        make_compress_idemp(&phi, 4, 2);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        slong d = ae.A.dim_A;
        acb_mat_t P, Q;
        diag_proj_at(P, n, 1);             /* P = |e1><e1|, in A */
        acb_mat_init(Q, n, n); acb_mat_zero(Q);
        acb_set_si(acb_mat_entry(Q, 0, 0), 1);
        for (slong i = 2; i < 4; i++) acb_set_si(acb_mat_entry(Q, i, i), 1); /* diag(1,0,1,1) */
        AIC_CHECK_MSG(aic_corner_dim_S(&ae.A, Q, Q, prec) == 1,
                      "T10(ii): Q=diag(1,0,1,1) not 1d in compress(4,2)");
        acb_mat_t CoPQ, CoQ, Qtilde;
        acb_mat_init(CoPQ, d, d); acb_mat_init(CoQ, d, d); acb_mat_init(Qtilde, n, n);
        aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
        aic_corner_Co(CoQ, &ae.A, Q, Q, prec);
        aic_corner_Ptilde(Qtilde, &ae.A, Q, prec);
        acb_mat_t *C;
        slong dPQ;
        aic_corner_extract(&C, &dPQ, CoPQ, &ae.A, prec);
        AIC_CHECK_MSG(dPQ >= 1, "T10(ii): empty corner");

        /* ip via the production ip_1d (star) vs a PLAIN-product surrogate:
         * ip_plain(Y,X) = <Qtilde, Co_Q(Y^dag X plain)>_F / <Qtilde,Qtilde>_F. */
        acb_t ipstar, ipplain, den, num, t;
        acb_init(ipstar); acb_init(ipplain); acb_init(den); acb_init(num); acb_init(t);
        frob_ip(den, Qtilde, Qtilde, prec);
        double vac = 0.0, x2gap = 0.0;
        for (slong m = 0; m < dPQ; m++) {
            aic_corner_ip_1d(ipstar, &ae.A, Qtilde, CoQ, C[m], C[m], prec);
            /* plain surrogate */
            acb_mat_t Yd, plain, W;
            acb_mat_init(Yd, n, n); acb_mat_init(plain, n, n); acb_mat_init(W, n, n);
            acb_mat_conjugate_transpose(Yd, C[m]);
            acb_mat_mul(plain, Yd, C[m], prec);             /* NO star */
            aic_corner_apply(W, &ae.A, CoQ, plain, prec);
            frob_ip(num, Qtilde, W, prec);
            acb_div(ipplain, num, den, prec);
            acb_sub(t, ipstar, ipplain, prec);
            double v = abs_acb(t);
            if (v > vac) vac = v;
            /* also report |<X,X> - ||X||^2| (the .tex:1130 closeness, eta=0 => ~0). */
            arb_t nx; arb_init(nx);
            aic_mat_opnorm(nx, C[m], prec);
            double nn = dd(nx) * dd(nx);
            double xxre = arf_get_d(arb_midref(acb_realref(ipstar)), ARF_RND_NEAR);
            double g = fabs(xxre - nn);
            if (g > x2gap) x2gap = g;
            arb_clear(nx);
            acb_mat_clear(W); acb_mat_clear(plain); acb_mat_clear(Yd);
        }
        printf("T10(ii) oblique non-vacuity (compress(4,2)): ||ip_star-ip_plain||"
               "=%.4f (star load-bearing) |<X,X>-||X||^2|=%.2e (eta=0)\n", vac, x2gap);
        AIC_CHECK_MSG(vac > 0.05,
                      "T10(ii): star-vs-plain ip gap %.3e <= 0.05 -- star not "
                      "exercised (fixture washes it out)", vac);
        AIC_CHECK_MSG(x2gap < 1e-9,
                      "T10(ii): |<X,X>-||X||^2|=%.3e not ~0 at eta=0", x2gap);

        acb_clear(t); acb_clear(num); acb_clear(den); acb_clear(ipplain); acb_clear(ipstar);
        for (slong m = 0; m < dPQ; m++) acb_mat_clear(C[m]);
        if (C) flint_free(C);
        acb_mat_clear(Qtilde); acb_mat_clear(CoQ); acb_mat_clear(CoPQ);
        acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }

    /* (iii) eta>0 closeness constant c: |<X,X>-||X||^2| <= c(eta) on dep+conj(4). */
    {
        slong n = 4;
        const double ts[] = {0.02, 0.06};
        double cmax = 0.0;
        for (int it = 0; it < 2; it++) {
            aic_ucp_kraus phi;
            make_eta_family(&phi, n, ts[it], prec);
            double eta = eta_proxy(&phi, prec);
            aic_assoc_ecstar ae;
            aic_assoc_ecstar_from_phi(&ae, &phi, prec);
            slong d = ae.A.dim_A;
            acb_mat_t P, Q;
            diag_proj(P, n, 2);
            diag_proj_at(Q, n, 0);
            if (aic_corner_dim_S(&ae.A, Q, Q, prec) != 1) {
                acb_mat_clear(Q); acb_mat_clear(P);
                aic_assoc_ecstar_clear(&ae); aic_ucp_kraus_clear(&phi);
                continue;
            }
            acb_mat_t CoPQ, CoQ, Qtilde;
            acb_mat_init(CoPQ, d, d); acb_mat_init(CoQ, d, d); acb_mat_init(Qtilde, n, n);
            aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
            aic_corner_Co(CoQ, &ae.A, Q, Q, prec);
            aic_corner_Ptilde(Qtilde, &ae.A, Q, prec);
            acb_mat_t *C;
            slong dPQ;
            aic_corner_extract(&C, &dPQ, CoPQ, &ae.A, prec);
            double worst = 0.0;
            acb_t xx; acb_init(xx);
            for (slong m = 0; m < dPQ; m++) {
                aic_corner_ip_1d(xx, &ae.A, Qtilde, CoQ, C[m], C[m], prec);
                arb_t nx; arb_init(nx);
                aic_mat_opnorm(nx, C[m], prec);
                double nn = dd(nx) * dd(nx);
                double xxre = arf_get_d(arb_midref(acb_realref(xx)), ARF_RND_NEAR);
                double g = fabs(xxre - nn);
                if (g > worst) worst = g;
                arb_clear(nx);
            }
            acb_clear(xx);
            double c = eta > 1e-14 ? worst / eta : 0.0;
            if (c > cmax) cmax = c;
            printf("T10(iii) t=%.2f eta=%.4e: max|<X,X>-||X||^2|=%.4e (c=%.3f)\n",
                   ts[it], eta, worst, c);
            AIC_CHECK_MSG(c < 20.0, "T10(iii): closeness c=%.3f exceeds 20", c);
            for (slong m = 0; m < dPQ; m++) acb_mat_clear(C[m]);
            if (C) flint_free(C);
            acb_mat_clear(Qtilde); acb_mat_clear(CoQ); acb_mat_clear(CoPQ);
            acb_mat_clear(Q); acb_mat_clear(P);
            aic_assoc_ecstar_clear(&ae);
            aic_ucp_kraus_clear(&phi);
        }
        printf("T10(iii) -> measured closeness c max = %.3f\n", cmax);
    }
}

/* T11 Ha-map: build Ha^Q_{P,R}(Z), test Ha_dag exact, the ||Ha(Z)(X)-Z.X|| bound,
 * and the Ha^Q_{P,P} homomorphism property.
 *
 * MUTATION-PROOFS (manual, RED documented):
 *  - wrong-Co (CoQR -> CoPQ in the Y^dag.Z cdot of aic_corner_ha): on the eta=0
 *    M_n homomorphism oracle below, the Ha(Z)(X)=Z.X residual jumps from ~0 to O(1)
 *    (CONFIRMED 5.00e-1 with P=R=diag(1,1,0,0),Q=e0), restored.
 *  - Frobenius-Gram (G = I, since {CPQ_l} are Frobenius-orthonormal, instead of the
 *    lem_PQ_Hilb Gram): at eta=0 the two Grams COINCIDE (GNS) so this is BLIND on
 *    every eta=0 oracle (CONFIRMED: with G=I forced, T11(a/b/c) all stay machine-
 *    zero). At eta>0 the lem_PQ_Hilb Gram deviates from I by O(eta) -- but on the
 *    dep+conj(4) fixture dim S_{P,Q}=1 (1d Q gives a 1d corner here), so G is 1x1
 *    and the G=I mutation only RESCALES a scalar: the Ha OUTPUT moves by < the
 *    c<50 threshold (Ha-vs-cdot c 0.000 -> 0.004, homdef/eta -> 0.004), so the
 *    mutation is effectively BLIND on the Ha output even at eta>0 (DOCUMENTED, not
 *    papered over). The substantive distinction IS witnessed: T11(eta) asserts the
 *    lem_PQ_Hilb Gram ||G-I|| is NONZERO (~1.5e-4) and O(eta), i.e. the two metrics
 *    provably DIFFER. Full multi-dim-corner output teeth await an eta>0 fixture with
 *    dim S_{P,Q}>1 for 1d Q (filed as a follow-up bead).
 *  - Ha_dag teeth: the identity is structurally satisfiable by Ha==0; the non-
 *    vacuity guard (||Ha||>0.1) + the eta=0 closed form Ha(Z)(X)=Z.X give it teeth. */
static void test_t11(void)
{
    const slong prec = 256;

    /* (a) eta=0 oracle on identity(4): A=M_n, star=plain. P=R=diag(1,1,0,0) (so
     * Ha^Q_{P,P} is the homomorphism case), 1d Q=|e0><e0|. Closed form: at eta=0
     * the algebra is an exact C* algebra, so Ha(Z)(X) = Z.X EXACTLY (the O(d+e)
     * bound is 0). Test column-by-column: ||Ha(Z)(C^{RQ}_m) - Z.C^{RQ}_m|| ~ 0,
     * for several Z in S_P. Also the homomorphism Ha(Z.W) = Ha(Z) Ha(W) EXACT, and
     * ||Ha(Ptilde) - I|| ~ 0. */
    slong n = 4;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;
    acb_mat_t P, Q;
    diag_proj(P, n, 2);          /* P=R=diag(1,1,0,0) */
    diag_proj_at(Q, n, 0);
    AIC_CHECK_MSG(aic_corner_dim_S(&ae.A, Q, Q, prec) == 1, "T11: Q not 1d");

    /* corner bases: S_{P,Q}=S_{R,Q} (P==R), S_P. */
    acb_mat_t CoPQ, CoPP;
    acb_mat_init(CoPQ, d, d); acb_mat_init(CoPP, d, d);
    aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
    aic_corner_Co(CoPP, &ae.A, P, P, prec);
    acb_mat_t *CPQ, *CP;
    slong dPQ, dP;
    aic_corner_extract(&CPQ, &dPQ, CoPQ, &ae.A, prec);
    aic_corner_extract(&CP, &dP, CoPP, &ae.A, prec);
    AIC_CHECK_MSG(dPQ >= 1 && dP >= 1, "T11: empty corner (dPQ=%ld dP=%ld)",
                  (long) dPQ, (long) dP);

    /* Ha(Z)(X) = Z.X oracle + non-vacuity (||Ha|| > 0.1). */
    double ha_vs_cdot = 0.0, ha_norm = 0.0;
    acb_mat_t Ha, hax, zx, dhz;
    acb_mat_init(Ha, dPQ, dPQ);     /* P==R so d_RQ == d_PQ */
    acb_mat_init(hax, n, n); acb_mat_init(zx, n, n); acb_mat_init(dhz, n, n);
    for (slong zi = 0; zi < dP; zi++) {
        aic_corner_ha(Ha, &ae.A, CP[zi], P, P, Q, prec);
        double hn = opnorm_d(Ha, prec);
        if (hn > ha_norm) ha_norm = hn;
        for (slong m = 0; m < dPQ; m++) {
            ha_apply_col(hax, Ha, CPQ, dPQ, m, n, prec);          /* Ha(Z)(C_m) */
            aic_corner_cdot(zx, &ae.A, CoPQ, CP[zi], CPQ[m], prec); /* Z . C_m */
            acb_mat_sub(dhz, hax, zx, prec);
            double e = mid_opnorm_d(dhz, prec);
            if (e > ha_vs_cdot) ha_vs_cdot = e;
        }
    }
    printf("T11(a) eta=0 Ha(Z)(X)=Z.X oracle (id(4)): max||Ha(Z)(X)-Z.X||=%.3e "
           "||Ha||_op=%.3f (non-vacuity)\n", ha_vs_cdot, ha_norm);
    AIC_CHECK_MSG(ha_vs_cdot < 1e-9,
                  "T11(a): Ha(Z)(X) != Z.X at eta=0 (%.3e) -- wrong Co in cdot chain?",
                  ha_vs_cdot);
    AIC_CHECK_MSG(ha_norm > 0.1,
                  "T11(a): ||Ha||_op=%.3f <= 0.1 -- Ha structurally ~0, identity "
                  "tests VACUOUS", ha_norm);

    /* Ha_dag EXACT (.tex:1153): Ha^Q_{R,P}(Z^dag) = Ha^Q_{P,R}(Z)^dag. With P==R
     * here Z^dag is in S_{R,P}=S_{P,P}=S_P, so both sides are dPQ x dPQ. Compare
     * the matrices (the Ha-map adjoint is w.r.t. the lem_PQ_Hilb inner product, so
     * the matrix adjoint is G^{-1} Ha^dag G; at eta=0 G=I exactly so it is the
     * plain conjugate-transpose). NON-VACUITY: an exact identity can be a "test that
     * can't fail" if it is M==M structurally. Here it is NOT: the Ha(Z) matrices are
     * genuinely NON-Hermitian (asym = max||Ha(Z)-Ha(Z)^H|| measured 1.0 on this M_4
     * fixture, dim S_{P,Q}=2), and Ha_dag relates Ha(Z^dag) to Ha(Z)^dag for DISTINCT
     * Z, Z^dag -- a real constraint. We assert asym > 0.1 so the identity has teeth;
     * the wrong-Co mutation (which corrupts Ha) also breaks this (CONFIRMED via
     * T11(a)). */
    double hadag = 0.0, asym = 0.0;
    acb_mat_t HaZ, Zd, HaZd, HaZdag, dd_, HaH;
    acb_mat_init(HaZ, dPQ, dPQ); acb_mat_init(Zd, n, n);
    acb_mat_init(HaZd, dPQ, dPQ); acb_mat_init(HaZdag, dPQ, dPQ);
    acb_mat_init(dd_, dPQ, dPQ); acb_mat_init(HaH, dPQ, dPQ);
    for (slong zi = 0; zi < dP; zi++) {
        aic_corner_ha(HaZ, &ae.A, CP[zi], P, P, Q, prec);    /* Ha^Q_{P,P}(Z)  */
        acb_mat_conjugate_transpose(Zd, CP[zi]);             /* Z^dag in S_P    */
        aic_corner_ha(HaZd, &ae.A, Zd, P, P, Q, prec);       /* Ha^Q_{P,P}(Z^dag) */
        acb_mat_conjugate_transpose(HaZdag, HaZ);            /* Ha(Z)^dag (G=I) */
        acb_mat_sub(dd_, HaZd, HaZdag, prec);
        double e = opnorm_d(dd_, prec);
        if (e > hadag) hadag = e;
        /* non-vacuity: is Ha(Z) genuinely non-Hermitian? */
        acb_mat_sub(HaH, HaZ, HaZdag, prec);
        double a = opnorm_d(HaH, prec);
        if (a > asym) asym = a;
    }
    printf("T11(b) Ha_dag exact (eta=0): max||Ha(Z^dag) - Ha(Z)^dag||=%.3e "
           "(non-vacuity: max||Ha(Z)-Ha(Z)^H||=%.3f)\n", hadag, asym);
    AIC_CHECK_MSG(hadag < 1e-9, "T11(b): Ha_dag identity violated (%.3e)", hadag);
    AIC_CHECK_MSG(asym > 0.1,
                  "T11(b): Ha matrices Hermitian (asym=%.3f) -- Ha_dag is M==M, "
                  "VACUOUS; need a fixture with non-Hermitian Ha", asym);

    /* Homomorphism (.tex:1157): Ha(Z.W) = Ha(Z) Ha(W) EXACT at eta=0; and
     * ||Ha(Ptilde) - I|| ~ 0. The matrices multiply (B(S_{P,Q}) composition). */
    double hom = 0.0;
    acb_mat_t HaW, HaZW, ZW, prodHa, dh;
    acb_mat_init(HaW, dPQ, dPQ); acb_mat_init(HaZW, dPQ, dPQ);
    acb_mat_init(ZW, n, n); acb_mat_init(prodHa, dPQ, dPQ); acb_mat_init(dh, dPQ, dPQ);
    for (slong zi = 0; zi < dP; zi++)
        for (slong wi = 0; wi < dP; wi++) {
            aic_corner_ha(HaZ, &ae.A, CP[zi], P, P, Q, prec);
            aic_corner_ha(HaW, &ae.A, CP[wi], P, P, Q, prec);
            aic_corner_cdot(ZW, &ae.A, CoPP, CP[zi], CP[wi], prec);  /* Z.W in S_P */
            aic_corner_ha(HaZW, &ae.A, ZW, P, P, Q, prec);           /* Ha(Z.W)   */
            acb_mat_mul(prodHa, HaZ, HaW, prec);                     /* Ha(Z)Ha(W)*/
            acb_mat_sub(dh, HaZW, prodHa, prec);
            double e = opnorm_d(dh, prec);
            if (e > hom) hom = e;
        }
    /* ||Ha(Ptilde) - I|| (.tex:1157): Ptilde the unit of S_P. */
    acb_mat_t Ptil, HaPt, eye, dpt;
    acb_mat_init(Ptil, n, n); acb_mat_init(HaPt, dPQ, dPQ);
    acb_mat_init(eye, dPQ, dPQ); acb_mat_init(dpt, dPQ, dPQ);
    aic_corner_Ptilde(Ptil, &ae.A, P, prec);
    aic_corner_ha(HaPt, &ae.A, Ptil, P, P, Q, prec);
    acb_mat_one(eye);
    acb_mat_sub(dpt, HaPt, eye, prec);
    double hapt = opnorm_d(dpt, prec);
    printf("T11(c) eta=0 homomorphism (id(4)): max||Ha(Z.W)-Ha(Z)Ha(W)||=%.3e "
           "||Ha(Ptilde)-I||=%.3e\n", hom, hapt);
    AIC_CHECK_MSG(hom < 1e-9, "T11(c): Ha^Q_{P,P} not a homomorphism at eta=0 (%.3e)", hom);
    AIC_CHECK_MSG(hapt < 1e-9, "T11(c): ||Ha(Ptilde)-I||=%.3e not ~0 at eta=0", hapt);

    acb_mat_clear(dpt); acb_mat_clear(eye); acb_mat_clear(HaPt); acb_mat_clear(Ptil);
    acb_mat_clear(dh); acb_mat_clear(prodHa); acb_mat_clear(ZW);
    acb_mat_clear(HaZW); acb_mat_clear(HaW);
    acb_mat_clear(HaH); acb_mat_clear(dd_); acb_mat_clear(HaZdag); acb_mat_clear(HaZd);
    acb_mat_clear(Zd); acb_mat_clear(HaZ);
    acb_mat_clear(dhz); acb_mat_clear(zx); acb_mat_clear(hax); acb_mat_clear(Ha);
    for (slong m = 0; m < dP; m++) acb_mat_clear(CP[m]);
    if (CP) flint_free(CP);
    for (slong m = 0; m < dPQ; m++) acb_mat_clear(CPQ[m]);
    if (CPQ) flint_free(CPQ);
    acb_mat_clear(CoPP); acb_mat_clear(CoPQ);
    acb_mat_clear(Q); acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T11(d) P != R arm (the GENERAL Ha^Q_{P,R}, .tex:1146). The homomorphism oracle
 * test_t11 uses P==R, where S_{P,Q}==S_{R,Q} and CoPQ==CoRQ, so the Z.X-step Co
 * confusion (CoPQ vs CoRQ) is INVISIBLE (a fixture degeneracy: that wrong-Co
 * mutation stays GREEN under P==R). With DISTINCT P=diag(1,1,0,0), R=diag(0,0,1,1),
 * Q=|e0><e0| on identity(4): S_{P,R} (Z's, dim 4), S_{R,Q} (X inputs, dim 2),
 * S_{P,Q} (Ha rows, dim 2), S_Q (dim 1) are all DISTINCT corners, so every Co in
 * the cdot chain is genuinely different. At eta=0 the C* algebra is exact so
 * Ha^Q_{P,R}(Z)(X) = Z.X EXACTLY. MUTATION-PROVEN (both RED at 5.00e-1 here,
 * BLIND under test_t11's P==R fixture where CoPQ==CoRQ and S_Q's corner collapses):
 *   (A) Z.X step CoPQ -> CoRQ in aic_corner_ha,
 *   (B) scalar term1 step CoQ -> CoQR in aic_corner_ha.
 * Both restored. (The CoQR step is already caught by test_t11's P==R arm at 5.00e-1.) */
static void test_t11_pneqr(void)
{
    const slong prec = 256, n = 4;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;
    acb_mat_t P, R, Q;
    diag_proj(P, n, 2);                                   /* diag(1,1,0,0) */
    acb_mat_init(R, n, n); acb_mat_zero(R);
    for (slong i = 2; i < 4; i++) acb_set_si(acb_mat_entry(R, i, i), 1); /* diag(0,0,1,1) */
    diag_proj_at(Q, n, 0);
    AIC_CHECK_MSG(aic_corner_dim_S(&ae.A, Q, Q, prec) == 1, "T11(d): Q not 1d");

    acb_mat_t CoPR, CoRQ, CoPQ;
    acb_mat_init(CoPR, d, d); acb_mat_init(CoRQ, d, d); acb_mat_init(CoPQ, d, d);
    aic_corner_Co(CoPR, &ae.A, P, R, prec);
    aic_corner_Co(CoRQ, &ae.A, R, Q, prec);
    aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
    acb_mat_t *CPR, *CRQ, *CPQ;
    slong dPR, dRQ, dPQ;
    aic_corner_extract(&CPR, &dPR, CoPR, &ae.A, prec);
    aic_corner_extract(&CRQ, &dRQ, CoRQ, &ae.A, prec);
    aic_corner_extract(&CPQ, &dPQ, CoPQ, &ae.A, prec);
    AIC_CHECK_MSG(dPR >= 1 && dRQ >= 1 && dPQ >= 1,
                  "T11(d): empty corner (dPR=%ld dRQ=%ld dPQ=%ld)",
                  (long) dPR, (long) dRQ, (long) dPQ);

    double worst = 0.0, hanorm = 0.0;
    acb_mat_t Ha, hax, zx, dh;
    acb_mat_init(Ha, dPQ, dRQ);
    acb_mat_init(hax, n, n); acb_mat_init(zx, n, n); acb_mat_init(dh, n, n);
    for (slong zi = 0; zi < dPR; zi++) {
        aic_corner_ha(Ha, &ae.A, CPR[zi], P, R, Q, prec);      /* Ha^Q_{P,R}(Z) */
        double hn = opnorm_d(Ha, prec);
        if (hn > hanorm) hanorm = hn;
        for (slong m = 0; m < dRQ; m++) {
            ha_apply_col(hax, Ha, CPQ, dPQ, m, n, prec);        /* Ha(Z)(C^{RQ}_m) */
            aic_corner_cdot(zx, &ae.A, CoPQ, CPR[zi], CRQ[m], prec); /* Z . X */
            acb_mat_sub(dh, hax, zx, prec);
            double e = mid_opnorm_d(dh, prec);
            if (e > worst) worst = e;
        }
    }
    printf("T11(d) P!=R eta=0 Ha^Q_{P,R}(Z)(X)=Z.X (id(4), dPR=%ld dRQ=%ld dPQ=%ld): "
           "max||Ha(Z)(X)-Z.X||=%.3e ||Ha||_op=%.3f\n", (long) dPR, (long) dRQ,
           (long) dPQ, worst, hanorm);
    AIC_CHECK_MSG(worst < 1e-9,
                  "T11(d): Ha^Q_{P,R}(Z)(X) != Z.X at eta=0 with P!=R (%.3e) -- "
                  "wrong Co in the cdot chain (Z.X step CoPQ vs CoRQ)?", worst);
    AIC_CHECK_MSG(hanorm > 0.1, "T11(d): ||Ha||=%.3f <= 0.1 VACUOUS", hanorm);

    acb_mat_clear(dh); acb_mat_clear(zx); acb_mat_clear(hax); acb_mat_clear(Ha);
    for (slong m = 0; m < dPQ; m++) acb_mat_clear(CPQ[m]);
    if (CPQ) flint_free(CPQ);
    for (slong m = 0; m < dRQ; m++) acb_mat_clear(CRQ[m]);
    if (CRQ) flint_free(CRQ);
    for (slong m = 0; m < dPR; m++) acb_mat_clear(CPR[m]);
    if (CPR) flint_free(CPR);
    acb_mat_clear(CoPQ); acb_mat_clear(CoRQ); acb_mat_clear(CoPR);
    acb_mat_clear(Q); acb_mat_clear(R); acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T11(disjoint) — the term2 Co-index teeth (FIX 3, bead aic). In aic_corner_ha the
 * scalar term2 Y^dag.(Z.X) is compressed with Co_Q (a S_Q target). All other Ha
 * fixtures (test_t11, test_t11_pneqr) have range(Q) SUBSET range(P): there
 * Co_{P,Q} = Co_Q on Q-supported elements, so mutating that term2 Co_Q -> Co_{P,Q}
 * is BLIND (stays green). This DISJOINT-SUPPORT fixture breaks that degeneracy:
 *   identity(4) (A = M_4, star = plain), P = diag(0,0,1,1), Q = |e0><e0|,
 * so range(Q) = span(e0) is DISJOINT from range(P) = span(e2,e3) -- Co_{P,Q} !=
 * Co_Q. dim S_Q = 1, dim S_{P,Q} = 2 (verified). At eta=0 the C* algebra is exact
 * so Ha^Q_{P,P}(Z)(X) = Z.X EXACTLY. MUTATION-PROVEN (source perturbation in
 * aic_corner_ha: the term2 cdot Co_Q -> CoPQ): the Ha(Z)(X)=Z.X residual jumps from
 * 0 to 5.000e-1 (RED) on THIS fixture, whereas it stays 0 under test_t11's
 * range(Q) subset range(P) fixture. Restored. */
static void test_t11_disjoint(void)
{
    const slong prec = 256, n = 4;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;

    acb_mat_t P, Q;
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    acb_set_si(acb_mat_entry(P, 2, 2), 1);
    acb_set_si(acb_mat_entry(P, 3, 3), 1);   /* P = diag(0,0,1,1) */
    diag_proj_at(Q, n, 0);                    /* Q = |e0><e0|, range disjoint from P */
    AIC_CHECK_MSG(aic_corner_dim_S(&ae.A, Q, Q, prec) == 1, "T11(disjoint): Q not 1d");
    slong dSPQ = aic_corner_dim_S(&ae.A, P, Q, prec);
    AIC_CHECK_MSG(dSPQ == 2,
                  "T11(disjoint): dim S_{P,Q}=%ld != 2 (range(Q)=span(e0) x "
                  "range(P)=span(e2,e3))", (long) dSPQ);

    acb_mat_t CoPQ, CoPP;
    acb_mat_init(CoPQ, d, d);
    acb_mat_init(CoPP, d, d);
    aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
    aic_corner_Co(CoPP, &ae.A, P, P, prec);
    acb_mat_t *CPQ, *CP;
    slong dPQ, dP;
    aic_corner_extract(&CPQ, &dPQ, CoPQ, &ae.A, prec);
    aic_corner_extract(&CP, &dP, CoPP, &ae.A, prec);
    AIC_CHECK_MSG(dPQ == 2 && dP >= 1, "T11(disjoint): dPQ=%ld dP=%ld", (long) dPQ,
                  (long) dP);

    double worst = 0.0, hanorm = 0.0;
    acb_mat_t Ha, hax, zx, dh;
    acb_mat_init(Ha, dPQ, dPQ);   /* R==P so d_RQ == d_PQ */
    acb_mat_init(hax, n, n);
    acb_mat_init(zx, n, n);
    acb_mat_init(dh, n, n);
    for (slong zi = 0; zi < dP; zi++) {
        aic_corner_ha(Ha, &ae.A, CP[zi], P, P, Q, prec);
        double hn = opnorm_d(Ha, prec);
        if (hn > hanorm) hanorm = hn;
        for (slong m = 0; m < dPQ; m++) {
            ha_apply_col(hax, Ha, CPQ, dPQ, m, n, prec);
            aic_corner_cdot(zx, &ae.A, CoPQ, CP[zi], CPQ[m], prec);
            acb_mat_sub(dh, hax, zx, prec);
            double e = mid_opnorm_d(dh, prec);
            if (e > worst) worst = e;
        }
    }
    printf("T11(disjoint) term2-Co teeth (id(4), P=diag(0,0,1,1) Q=|e0><e0|, "
           "range disjoint, dPQ=2): max||Ha(Z)(X)-Z.X||=%.3e ||Ha||=%.3f\n",
           worst, hanorm);
    AIC_CHECK_MSG(worst < 1e-9,
                  "T11(disjoint): Ha(Z)(X) != Z.X at eta=0 (%.3e) -- term2 Co_Q "
                  "wrong (e.g. CoPQ instead of CoQ), exposed by disjoint supports",
                  worst);
    AIC_CHECK_MSG(hanorm > 0.1, "T11(disjoint): ||Ha||=%.3f <= 0.1 VACUOUS", hanorm);

    acb_mat_clear(dh);
    acb_mat_clear(zx);
    acb_mat_clear(hax);
    acb_mat_clear(Ha);
    for (slong m = 0; m < dP; m++) acb_mat_clear(CP[m]);
    if (CP) flint_free(CP);
    for (slong m = 0; m < dPQ; m++) acb_mat_clear(CPQ[m]);
    if (CPQ) flint_free(CPQ);
    acb_mat_clear(CoPP);
    acb_mat_clear(CoPQ);
    acb_mat_clear(Q);
    acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T11 eta>0 arms: ||Ha(Z)(X)-Z.X||/eta and the homomorphism defect / eta bounded
 * (the downstream-critical O(delta+eps)-homomorphism, .tex:1157). Also the arm
 * where the FROBENIUS-Gram mutation goes RED (its c grows). */
static void test_t11_eta(void)
{
    const slong prec = 256, n = 4;
    const double ts[] = {0.02, 0.06};
    double cmax_prod = 0.0, cmax_hom = 0.0;
    for (int it = 0; it < 2; it++) {
        aic_ucp_kraus phi;
        make_eta_family(&phi, n, ts[it], prec);
        double eta = eta_proxy(&phi, prec);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        slong d = ae.A.dim_A;
        acb_mat_t P, Q;
        diag_proj(P, n, 2);
        diag_proj_at(Q, n, 0);
        if (aic_corner_dim_S(&ae.A, Q, Q, prec) != 1) {
            acb_mat_clear(Q); acb_mat_clear(P);
            aic_assoc_ecstar_clear(&ae); aic_ucp_kraus_clear(&phi);
            continue;
        }
        acb_mat_t CoPQ, CoPP, CoQ, Qtilde;
        acb_mat_init(CoPQ, d, d); acb_mat_init(CoPP, d, d);
        acb_mat_init(CoQ, d, d); acb_mat_init(Qtilde, n, n);
        aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
        aic_corner_Co(CoPP, &ae.A, P, P, prec);
        aic_corner_Co(CoQ, &ae.A, Q, Q, prec);
        aic_corner_Ptilde(Qtilde, &ae.A, Q, prec);
        acb_mat_t *CPQ, *CP;
        slong dPQ, dP;
        aic_corner_extract(&CPQ, &dPQ, CoPQ, &ae.A, prec);
        aic_corner_extract(&CP, &dP, CoPP, &ae.A, prec);

        /* GRAM-DEVIATION teeth for the lem_PQ_Hilb-vs-Frobenius distinction
         * (.tex:1130). The extracted {CPQ_l} are Frobenius-orthonormal, so the
         * FROBENIUS Gram is EXACTLY I; the lem_PQ_Hilb Gram G_{lm} = <C_l|C_m>
         * deviates from I by O(eta). We measure ||G - I||_op and assert it is (i)
         * NONZERO at eta>0 (so the two metrics genuinely DIFFER -- the substantive
         * content the G=I mutation erases) and (ii) <= c*eta (the .tex:1130 bound).
         * This is what gives the lem_PQ_Hilb-ip-vs-Frobenius-ip choice teeth: the
         * G=I mutation is BLIND on the Ha OUTPUT here only because dim S_{P,Q} = 1
         * for this 1d-Q fixture (a scalar rescaling), but the metrics provably
         * differ. The multi-dim-corner Ha output teeth await an eta>0 fixture with
         * dim S_{P,Q} > 1 for 1d Q (follow-up bead). */
        double gram_dev;
        {
            acb_mat_t G, GmI;
            acb_mat_init(G, dPQ, dPQ); acb_mat_init(GmI, dPQ, dPQ);
            for (slong l = 0; l < dPQ; l++)
                for (slong m = 0; m < dPQ; m++)
                    aic_corner_ip_1d(acb_mat_entry(G, l, m), &ae.A, Qtilde, CoQ,
                                     CPQ[m], CPQ[l], prec);
            acb_mat_one(GmI);
            acb_mat_sub(GmI, G, GmI, prec);
            gram_dev = mid_opnorm_d(GmI, prec);
            acb_mat_clear(GmI); acb_mat_clear(G);
        }

        acb_mat_t Ha, HaW, HaZW, ZW, prodHa, dh, hax, zx, dhz;
        acb_mat_init(Ha, dPQ, dPQ); acb_mat_init(HaW, dPQ, dPQ);
        acb_mat_init(HaZW, dPQ, dPQ); acb_mat_init(ZW, n, n);
        acb_mat_init(prodHa, dPQ, dPQ); acb_mat_init(dh, dPQ, dPQ);
        acb_mat_init(hax, n, n); acb_mat_init(zx, n, n); acb_mat_init(dhz, n, n);
        double wprod = 0.0, whom = 0.0;
        for (slong zi = 0; zi < dP; zi++) {
            aic_corner_ha(Ha, &ae.A, CP[zi], P, P, Q, prec);
            /* product closeness */
            arb_t nz; arb_init(nz);
            aic_mat_opnorm(nz, CP[zi], prec);
            double znorm = dd(nz);
            arb_clear(nz);
            for (slong m = 0; m < dPQ; m++) {
                ha_apply_col(hax, Ha, CPQ, dPQ, m, n, prec);
                aic_corner_cdot(zx, &ae.A, CoPQ, CP[zi], CPQ[m], prec);
                acb_mat_sub(dhz, hax, zx, prec);
                arb_t nxm; arb_init(nxm);
                aic_mat_opnorm(nxm, CPQ[m], prec);
                double e = mid_opnorm_d(dhz, prec) / (znorm * dd(nxm) + 1e-300);
                arb_clear(nxm);
                if (e > wprod) wprod = e;
            }
            /* homomorphism closeness */
            for (slong wi = 0; wi < dP; wi++) {
                aic_corner_ha(HaW, &ae.A, CP[wi], P, P, Q, prec);
                aic_corner_cdot(ZW, &ae.A, CoPP, CP[zi], CP[wi], prec);
                aic_corner_ha(HaZW, &ae.A, ZW, P, P, Q, prec);
                acb_mat_mul(prodHa, Ha, HaW, prec);   /* Ha holds Ha(Z) for this zi */
                acb_mat_sub(dh, HaZW, prodHa, prec);
                double e = opnorm_d(dh, prec);
                if (e > whom) whom = e;
            }
        }
        double cp = eta > 1e-14 ? wprod / eta : 0.0;
        double ch = eta > 1e-14 ? whom / eta : 0.0;
        if (cp > cmax_prod) cmax_prod = cp;
        if (ch > cmax_hom) cmax_hom = ch;
        double gc = eta > 1e-14 ? gram_dev / eta : 0.0;
        printf("T11(eta) t=%.2f eta=%.4e: ||Ha(Z)(X)-Z.X||/(eta||Z||||X||)=%.3f "
               "homdef/eta=%.3f ||G-I||=%.3e (Hilb!=Frob; /eta=%.3f)\n",
               ts[it], eta, cp, ch, gram_dev, gc);
        AIC_CHECK_MSG(cp < 50.0, "T11(eta): Ha-vs-cdot c=%.3f exceeds 50", cp);
        AIC_CHECK_MSG(ch < 50.0, "T11(eta): homomorphism c=%.3f exceeds 50", ch);
        AIC_CHECK_MSG(gram_dev > 1e-9,
                      "T11(eta): lem_PQ_Hilb Gram ||G-I||=%.3e ~ 0 -- the Hilbert ip "
                      "coincides with the Frobenius ip, so the ip choice is VACUOUS "
                      "(fixture has no eta-deviation)", gram_dev);
        AIC_CHECK_MSG(gc < 50.0,
                      "T11(eta): ||G-I||/eta=%.3f exceeds 50 -- not O(delta+eps)", gc);
        acb_mat_clear(dhz); acb_mat_clear(zx); acb_mat_clear(hax);
        acb_mat_clear(dh); acb_mat_clear(prodHa); acb_mat_clear(ZW);
        acb_mat_clear(HaZW); acb_mat_clear(HaW); acb_mat_clear(Ha);
        for (slong m = 0; m < dP; m++) acb_mat_clear(CP[m]);
        if (CP) flint_free(CP);
        for (slong m = 0; m < dPQ; m++) acb_mat_clear(CPQ[m]);
        if (CPQ) flint_free(CPQ);
        acb_mat_clear(Qtilde); acb_mat_clear(CoQ);
        acb_mat_clear(CoPP); acb_mat_clear(CoPQ);
        acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
    printf("T11(eta) -> Ha-vs-cdot c max=%.3f, homomorphism c max=%.3f\n",
           cmax_prod, cmax_hom);
}

/* T11(mixconj) — FIX 1 regression fixture + the lem_PQ_Hilb-vs-Frobenius Gram
 * teeth (bead aic). The eta>0 fixtures in test_t11_eta all have dim S_{P,Q} = 1
 * for the 1d-Q corner (the dep+conj(4) algebra is diagonal-dominant), so the
 * Ha Gram is 1x1 and BOTH the blocker (the ||G-I|| guard never sees a near-zero
 * OFF-diagonal) and the Frobenius-vs-lem_PQ_Hilb distinction are degenerate. The
 * mixconj(5,3,0.02) fixture (a convex mix of compress_idemp(5,3) with its unitary
 * conjugate) is the first eta>0 instance with a genuine dim S_{P,Q} = 2 corner:
 *   P = projector onto span(e1,e2)  (in-A residual ~ 9.2e-3),
 *   Q = |e1><e1|                    (in-A residual ~ 4.8e-3, dim S_Q = 1),
 * giving eta ~ 1.30e-2, dim S_{P,Q} = 2, lem_PQ_Hilb Gram off-diagonal G[0,1] ~
 * 1.7e-6, ||G-I|| ~ 1.3e-5 (trivially in-basin). This is the regime lem_extension
 * needs (Ha^Q_{P,P} a homomorphism on a >1-dim corner at eta>0).
 *
 * THE BLOCKER. Before FIX 1, aic_corner_ha's ||G-I|| guard called aic_mat_opnorm,
 * whose Gram (G-I)^dag(G-I) off-diagonal relative-Hermiticity check FALSE-FAILS on
 * the ~1.7e-6 off-diagonal (radius exceeds the absolute floor, bead aic-2yo) and
 * SIGABRTs. (Confirmed: the routine SIGABRTed on this exact fixture before the
 * opnorm swap.) After FIX 1 (certified mid+radius upper bound) it RUNS.
 *
 * THE FROBENIUS-GRAM MUTATION (the lem_PQ_Hilb Gram's first teeth). The {C^{PQ}_l}
 * are Frobenius-orthonormal, so the FROBENIUS Gram is EXACTLY I; the lem_PQ_Hilb
 * Gram deviates by O(eta). On a 1d corner G=I only rescales a scalar (BLIND on the
 * Ha output, documented in test_t11_eta); on THIS 2x2 corner the off-diagonal
 * G[0,1] ~ 1.7e-6 is genuine, so forcing G=I in aic_corner_ha (the source mutation
 * "acb_mat_one(G)" replacing the aic_corner_ip_1d Gram build) MOVES the Ha output.
 * MUTATION-PROVEN (manual source perturbation, RED documented): the max
 * ||Ha(Z)(X)-Z.X|| residual jumps from 2.276e-6 (production) to 1.177e-5 (G=I
 * forced), tripping the < 1.0e-5 absolute teeth assert below (RED). Restored. The
 * c<... and gram_dev>0 asserts mirror test_t11_eta. */
static void test_t11_mixconj(void)
{
    const slong prec = 256, n = 5;
    aic_ucp_kraus phi;
    make_mixconj(&phi, 5, 3, 0.02, prec);
    double eta = eta_proxy(&phi, prec);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);
    slong d = ae.A.dim_A;

    /* P = projector onto span(e1,e2), Q = |e1><e1| (1d). */
    acb_mat_t P, Q;
    acb_mat_init(P, n, n);
    acb_mat_zero(P);
    acb_set_si(acb_mat_entry(P, 1, 1), 1);
    acb_set_si(acb_mat_entry(P, 2, 2), 1);
    diag_proj_at(Q, n, 1);
    AIC_CHECK_MSG(aic_corner_dim_S(&ae.A, Q, Q, prec) == 1,
                  "T11(mixconj): Q=|e1><e1| not 1d in mixconj(5,3)");
    slong dSPQ = aic_corner_dim_S(&ae.A, P, Q, prec);
    /* (i) the corner is genuinely 2-dimensional (the regime the blocker hit). */
    AIC_CHECK_MSG(dSPQ == 2,
                  "T11(mixconj): dim S_{P,Q}=%ld != 2 -- fixture lost its 2d corner",
                  (long) dSPQ);

    acb_mat_t CoPQ, CoPP, CoQ, Qtilde;
    acb_mat_init(CoPQ, d, d);
    acb_mat_init(CoPP, d, d);
    acb_mat_init(CoQ, d, d);
    acb_mat_init(Qtilde, n, n);
    aic_corner_Co(CoPQ, &ae.A, P, Q, prec);
    aic_corner_Co(CoPP, &ae.A, P, P, prec);
    aic_corner_Co(CoQ, &ae.A, Q, Q, prec);
    aic_corner_Ptilde(Qtilde, &ae.A, Q, prec);
    acb_mat_t *CPQ, *CP;
    slong dPQ, dP;
    aic_corner_extract(&CPQ, &dPQ, CoPQ, &ae.A, prec);
    aic_corner_extract(&CP, &dP, CoPP, &ae.A, prec);
    AIC_CHECK_MSG(dPQ == 2 && dP >= 2,
                  "T11(mixconj): extracted dPQ=%ld dP=%ld (expect dPQ=2)",
                  (long) dPQ, (long) dP);

    /* lem_PQ_Hilb Gram ||G-I|| (the metric that the Frobenius G=I mutation erases):
     * NONZERO and O(eta) here (the two metrics provably differ on a 2d corner). */
    double gram_dev;
    {
        acb_mat_t G, GmI;
        acb_mat_init(G, dPQ, dPQ);
        acb_mat_init(GmI, dPQ, dPQ);
        for (slong l = 0; l < dPQ; l++)
            for (slong m = 0; m < dPQ; m++)
                aic_corner_ip_1d(acb_mat_entry(G, l, m), &ae.A, Qtilde, CoQ,
                                 CPQ[m], CPQ[l], prec);
        acb_mat_one(GmI);
        acb_mat_sub(GmI, G, GmI, prec);
        gram_dev = mid_opnorm_d(GmI, prec);
        acb_mat_clear(GmI);
        acb_mat_clear(G);
    }

    /* (ii) Ha(Z)(X) ~ Z.X to O(eta): RUNS (no abort -- the FIX 1 regression) and
     * the residual is small. Also the ABSOLUTE-teeth value for the Frobenius-Gram
     * mutation: production max||Ha(Z)(X)-Z.X|| = 2.276e-6; the G=I mutation moves it
     * to 1.177e-5 (> 1.0e-5 -> RED). (iii) Ha^Q_{P,P} the O(eta)-homomorphism. */
    acb_mat_t Ha, hax, zx, dhz;
    acb_mat_init(Ha, dPQ, dPQ);   /* P==R so d_RQ == d_PQ */
    acb_mat_init(hax, n, n);
    acb_mat_init(zx, n, n);
    acb_mat_init(dhz, n, n);
    double ha_abs = 0.0, cprod = 0.0, ha_norm = 0.0;
    for (slong zi = 0; zi < dP; zi++) {
        aic_corner_ha(Ha, &ae.A, CP[zi], P, P, Q, prec);   /* RUNS (FIX 1) */
        double hn = mid_opnorm_d(Ha, prec);
        if (hn > ha_norm) ha_norm = hn;
        arb_t nz;
        arb_init(nz);
        aic_mat_opnorm(nz, CP[zi], prec);   /* CP exact (radius 0), no false-fail */
        double znorm = dd(nz);
        arb_clear(nz);
        for (slong m = 0; m < dPQ; m++) {
            ha_apply_col(hax, Ha, CPQ, dPQ, m, n, prec);
            aic_corner_cdot(zx, &ae.A, CoPQ, CP[zi], CPQ[m], prec);
            acb_mat_sub(dhz, hax, zx, prec);
            double e = mid_opnorm_d(dhz, prec);
            if (e > ha_abs) ha_abs = e;
            arb_t nxm;
            arb_init(nxm);
            aic_mat_opnorm(nxm, CPQ[m], prec);   /* CPQ exact, no false-fail */
            double cc = e / (eta * znorm * dd(nxm) + 1e-300);
            arb_clear(nxm);
            if (cc > cprod) cprod = cc;
        }
    }

    /* homomorphism Ha(Z.W) ~ Ha(Z)Ha(W) to O(eta) on the 2x2 corner. */
    acb_mat_t HaZ, HaW, HaZW, ZW, prodHa, dh;
    acb_mat_init(HaZ, dPQ, dPQ);
    acb_mat_init(HaW, dPQ, dPQ);
    acb_mat_init(HaZW, dPQ, dPQ);
    acb_mat_init(ZW, n, n);
    acb_mat_init(prodHa, dPQ, dPQ);
    acb_mat_init(dh, dPQ, dPQ);
    double whom = 0.0;
    for (slong zi = 0; zi < dP; zi++)
        for (slong wi = 0; wi < dP; wi++) {
            aic_corner_ha(HaZ, &ae.A, CP[zi], P, P, Q, prec);
            aic_corner_ha(HaW, &ae.A, CP[wi], P, P, Q, prec);
            aic_corner_cdot(ZW, &ae.A, CoPP, CP[zi], CP[wi], prec);
            aic_corner_ha(HaZW, &ae.A, ZW, P, P, Q, prec);
            acb_mat_mul(prodHa, HaZ, HaW, prec);
            acb_mat_sub(dh, HaZW, prodHa, prec);
            double e = mid_opnorm_d(dh, prec);
            if (e > whom) whom = e;
        }
    double chom = whom / (eta + 1e-300);

    printf("T11(mixconj) FIX1 regression (mixconj(5,3,0.02), dim S_{P,Q}=2): "
           "eta=%.4e ||Ha||=%.3f\n  Ha-vs-cdot abs=%.3e (c=%.3f) homdef=%.3e "
           "(c=%.3f) ||G-I||_Hilb=%.3e (Frob G=I -> 0; /eta=%.3f)\n",
           eta, ha_norm, ha_abs, cprod, whom, chom, gram_dev,
           gram_dev / (eta + 1e-300));
    /* RUNS (no SIGABRT) + 2d corner -- the FIX 1 regression. */
    AIC_CHECK_MSG(ha_norm > 0.1,
                  "T11(mixconj): ||Ha||=%.3f <= 0.1 -- Ha structurally ~0, VACUOUS",
                  ha_norm);
    /* O(eta) bounds (loose constants; report the measured c). */
    AIC_CHECK_MSG(cprod < 50.0,
                  "T11(mixconj): Ha-vs-cdot c=%.3f exceeds 50 -- not O(delta+eps)",
                  cprod);
    AIC_CHECK_MSG(chom < 50.0,
                  "T11(mixconj): homomorphism c=%.3f exceeds 50 -- not O(delta+eps)",
                  chom);
    /* ABSOLUTE Ha-vs-cdot teeth: production 2.276e-6 < 1.0e-5; the Frobenius-Gram
     * (G=I) source mutation moves it to 1.177e-5 > 1.0e-5 (RED, mutation-proven). */
    AIC_CHECK_MSG(ha_abs < 1.0e-5,
                  "T11(mixconj): max||Ha(Z)(X)-Z.X||=%.3e >= 1.0e-5 -- Ha output "
                  "corrupted (e.g. Frobenius G=I instead of the lem_PQ_Hilb Gram)",
                  ha_abs);
    /* the lem_PQ_Hilb Gram genuinely DIFFERS from Frobenius (G != I) and is O(eta):
     * the substantive content the G=I mutation erases. */
    AIC_CHECK_MSG(gram_dev > 1e-9,
                  "T11(mixconj): lem_PQ_Hilb Gram ||G-I||=%.3e ~ 0 -- coincides with "
                  "the Frobenius ip, ip choice VACUOUS", gram_dev);
    AIC_CHECK_MSG(gram_dev < 50.0 * eta,
                  "T11(mixconj): ||G-I||=%.3e exceeds 50 eta -- not O(delta+eps)",
                  gram_dev);

    acb_mat_clear(dh);
    acb_mat_clear(prodHa);
    acb_mat_clear(ZW);
    acb_mat_clear(HaZW);
    acb_mat_clear(HaW);
    acb_mat_clear(HaZ);
    acb_mat_clear(dhz);
    acb_mat_clear(zx);
    acb_mat_clear(hax);
    acb_mat_clear(Ha);
    for (slong m = 0; m < dP; m++) acb_mat_clear(CP[m]);
    if (CP) flint_free(CP);
    for (slong m = 0; m < dPQ; m++) acb_mat_clear(CPQ[m]);
    if (CPQ) flint_free(CPQ);
    acb_mat_clear(Qtilde);
    acb_mat_clear(CoQ);
    acb_mat_clear(CoPP);
    acb_mat_clear(CoPQ);
    acb_mat_clear(Q);
    acb_mat_clear(P);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* T12 lem_PQR (.tex:1162-1177): | ||X.Y||-||X||||Y|| | bound, 1d Q. */
static void test_t12(void)
{
    const slong prec = 256, n = 4;
    /* eta=0 on identity(4): exact, ||X.Y|| = ||X|| ||Y||. P=diag(1,1,0,0),
     * Q=|e0><e0| (1d), R=diag(0,0,1,1). */
    {
        aic_ucp_kraus phi;
        make_identity(&phi, n);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, prec);
        acb_mat_t P, Q, R;
        diag_proj(P, n, 2);
        diag_proj_at(Q, n, 0);
        acb_mat_init(R, n, n); acb_mat_zero(R);
        for (slong i = 2; i < 4; i++) acb_set_si(acb_mat_entry(R, i, i), 1);
        arb_t def; arb_init(def);
        aic_corner_pqr_defect(def, &ae.A, P, Q, R, prec);
        printf("T12 eta=0 lem_PQR (id(4)): max| ||X.Y||-||X||||Y|| |=%.3e\n", dd(def));
        AIC_CHECK_MSG(arb_le_dbl(def, 1e-9),
                      "T12: ||X.Y|| != ||X||||Y|| at eta=0 (%.3e)", dd(def));
        arb_clear(def);
        acb_mat_clear(R); acb_mat_clear(Q); acb_mat_clear(P);
        aic_assoc_ecstar_clear(&ae); aic_ucp_kraus_clear(&phi);
    }
    /* eta>0 dep+conj(4): | ||X.Y||-||X||||Y|| | / eta bounded. */
    {
        const double ts[] = {0.02, 0.06};
        double cmax = 0.0;
        for (int it = 0; it < 2; it++) {
            aic_ucp_kraus phi;
            make_eta_family(&phi, n, ts[it], prec);
            double eta = eta_proxy(&phi, prec);
            aic_assoc_ecstar ae;
            aic_assoc_ecstar_from_phi(&ae, &phi, prec);
            acb_mat_t P, Q, R;
            diag_proj(P, n, 2);
            diag_proj_at(Q, n, 0);
            acb_mat_init(R, n, n); acb_mat_zero(R);
            for (slong i = 2; i < 4; i++) acb_set_si(acb_mat_entry(R, i, i), 1);
            if (aic_corner_dim_S(&ae.A, Q, Q, prec) != 1) {
                acb_mat_clear(R); acb_mat_clear(Q); acb_mat_clear(P);
                aic_assoc_ecstar_clear(&ae); aic_ucp_kraus_clear(&phi);
                continue;
            }
            arb_t def; arb_init(def);
            aic_corner_pqr_defect(def, &ae.A, P, Q, R, prec);
            double c = eta > 1e-14 ? dd(def) / eta : 0.0;
            if (c > cmax) cmax = c;
            printf("T12 t=%.2f eta=%.4e: max| ||X.Y||-||X||||Y|| |=%.4e (c=%.3f)\n",
                   ts[it], eta, dd(def), c);
            AIC_CHECK_MSG(c < 50.0, "T12: lem_PQR c=%.3f exceeds 50", c);
            arb_clear(def);
            acb_mat_clear(R); acb_mat_clear(Q); acb_mat_clear(P);
            aic_assoc_ecstar_clear(&ae); aic_ucp_kraus_clear(&phi);
        }
        printf("T12 -> measured lem_PQR c max = %.3f\n", cmax);
    }
}

/* T13 lem_1d_proj + equivalence (.tex:1179-1187) on identity(4) (A=M_4). 1d
 * projectors are rank-1; P~Q iff dim S_{P,Q}=1. Diagonal rank-1's are pairwise
 * equivalent (dim S_{P,Q}=1 for any rank-1 P,Q in M_n), so the relation is the
 * full single-class on 1d projectors -- reflexive, symmetric, transitive. */
static void test_t13(void)
{
    const slong prec = 256, n = 4;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, prec);

    acb_mat_t e0, e1, e2;
    diag_proj_at(e0, n, 0);
    diag_proj_at(e1, n, 1);
    diag_proj_at(e2, n, 2);
    /* a NON-diagonal rank-1 (conjugated) to exercise off-diagonal corners. */
    acb_mat_t V, Pd, Pv;
    acb_mat_init(V, n, n);
    make_fixed_unitary(V, n, prec);
    diag_proj_at(Pd, n, 0);
    conj_proj(Pv, Pd, V, prec);

    /* lem_1d_proj: dim S_{P,Q} <= 1 for all 1d P,Q. */
    slong d00 = aic_corner_dim_S(&ae.A, e0, e0, prec);
    slong d01 = aic_corner_dim_S(&ae.A, e0, e1, prec);
    slong d12 = aic_corner_dim_S(&ae.A, e1, e2, prec);
    slong d0v = aic_corner_dim_S(&ae.A, e0, Pv, prec);
    printf("T13 lem_1d_proj dims: S_{e0,e0}=%ld S_{e0,e1}=%ld S_{e1,e2}=%ld "
           "S_{e0,Pv}=%ld (all <= 1)\n", (long) d00, (long) d01, (long) d12, (long) d0v);
    AIC_CHECK_MSG(d00 <= 1 && d01 <= 1 && d12 <= 1 && d0v <= 1,
                  "T13: dim S_{P,Q} > 1 for 1d P,Q (lem_1d_proj)");

    /* equivalence: reflexive, symmetric, transitive. */
    int r00 = aic_corner_equiv_1d(&ae.A, e0, e0, prec);          /* reflexive */
    int e01 = aic_corner_equiv_1d(&ae.A, e0, e1, prec);
    int e10 = aic_corner_equiv_1d(&ae.A, e1, e0, prec);          /* symmetric */
    int e12 = aic_corner_equiv_1d(&ae.A, e1, e2, prec);
    int e02 = aic_corner_equiv_1d(&ae.A, e0, e2, prec);          /* transitive */
    printf("T13 equiv: e0~e0=%d e0~e1=%d e1~e0=%d e1~e2=%d e0~e2=%d\n",
           r00, e01, e10, e12, e02);
    AIC_CHECK_MSG(r00 == 1, "T13: equivalence not reflexive");
    AIC_CHECK_MSG(e01 == e10, "T13: equivalence not symmetric (e0~e1=%d e1~e0=%d)",
                  e01, e10);
    AIC_CHECK_MSG(!(e01 && e12) || e02,
                  "T13: equivalence not transitive (e0~e1, e1~e2, but not e0~e2)");
    /* all diagonal rank-1's ARE equivalent in M_n => the transitive premise is
     * non-vacuously satisfied (e01 && e12 both true), so e02 is genuinely tested. */
    AIC_CHECK_MSG(e01 == 1 && e12 == 1 && e02 == 1,
                  "T13: rank-1's in M_n not all equivalent (premise vacuous?)");

    acb_mat_clear(Pv); acb_mat_clear(Pd); acb_mat_clear(V);
    acb_mat_clear(e2); acb_mat_clear(e1); acb_mat_clear(e0);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    test_t5();   /* headline eta=0 oracle first */
    test_t1();
    test_t2();
    test_t3();
    test_t4();
    test_t6();
    test_t7();   /* oblique oracle: gives build_L/extract teeth */
    test_t8_oblique();  /* Increment 2a: compressed product oracle (oblique star) */
    test_t8_unit();     /* Increment 2a: Ptilde unit of (S_P,.) */
    test_t8_bound();    /* Increment 2a: ||X.Y - X*Y|| O(eta) bound */
    test_t9();          /* Increment 2a: lem_alpha block bijection */
    test_t10();         /* Increment 2b: lem_PQ_Hilb inner product */
    test_t11();         /* Increment 2b: Ha-map (eta=0 oracle + Ha_dag + homom.) */
    test_t11_pneqr();   /* Increment 2b: Ha-map P!=R (full Co-index bookkeeping) */
    test_t11_disjoint();/* FIX 3: term2 Co-index teeth (range(Q) disjoint from range(P)) */
    test_t11_eta();     /* Increment 2b: Ha-map eta>0 bounds */
    test_t11_mixconj(); /* FIX 1: dim S_{P,Q}=2 eta>0 corner (blocker) + Frob-Gram teeth */
    test_t12();         /* Increment 2b: lem_PQR norm multiplicativity */
    test_t13();         /* Increment 2b: lem_1d_proj + equivalence */
    aic_test_report("test_corner");
    printf("OK test_corner\n");
    return 0;
}
