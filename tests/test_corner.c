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
#include <stdio.h>

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

int main(void)
{
    test_t5();   /* headline eta=0 oracle first */
    test_t1();
    test_t2();
    test_t3();
    test_t4();
    test_t6();
    test_t7();   /* oblique oracle: gives build_L/extract teeth */
    aic_test_report("test_corner");
    printf("OK test_corner\n");
    return 0;
}
