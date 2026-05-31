/* test_opspace_choi.c — cross-checks for the O2.1 adjoint Choi assemblers (bead
 * aic-pjr): aic_opspace_choi_vstar / aic_opspace_choi_vinvstar
 * (docs/research/opspace_o2_design.md §0.5 PINNED CONVENTION; the GOLDEN RULE Choi
 * J(f)[s*out+i, t*out+j] = f(E_st)[i,j] for the cb-norm SDP). See
 * include/aic_opspace.h and src/aic_opspace_choi.c for the contract.
 *
 * The fixture is the eta=0 exact-*-isomorphism v: B -> A built EXACTLY as
 * tests/test_opspace.c T1 does (make_block_cond_exp -> aic_assoc_ecstar_from_phi ->
 * aic_cstar_build, iso_def ~ 0). On this fixture v is a unital *-iso, so v* and
 * (v^{-1})* are well-understood and the checks have clean closed forms.
 *
 * THE CHECKS (Rule 5/6 — every check is an invariant with teeth, not "runs without
 * error"; mutation-proven RED at the bottom of each block comment):
 *   T1 Hermiticity. v, v^{-1} are HP => J(v*), J((v^{-1})*) are HERMITIAN.
 *      max_{ij} |J[i,j] - conj(J[j,i])| ~ 1e-12. MUTATION: dropping the conj in the
 *      assembler (place vE[i][a,b] not conj) breaks Hermiticity on the complex
 *      fixture -> RED.
 *   T2 independent recompute (catches an index swap a<->r / b<->c). For a handful
 *      of (a,b,r,c) recompute J(v*)[a*nB+r, b*nB+c] by the INDEPENDENT trace-pairing
 *      route v*(E_ab)[r,c] = conj( v(E_{rc}^{(n_B)})[a,b] ), evaluating v(E_{rc}) via
 *      the production aic_dhom_v_apply (no index logic shared with the assembler).
 *      Agreement ~1e-12. MUTATION: swapping the assembler slot to
 *      Jvs[a*nB + c, b*nB + r] -> RED.
 *   T3 eta=0 trace. Tr(J(f)) = Tr(f(I_in)) (GOLDEN RULE diagonal: only s==t, i==j
 *      survive). The EXPECTED value is computed INDEPENDENTLY in-test (sum of
 *      f(E_ss) traces, f via the assembler formula re-derived from vE / vinvB), then
 *      the assembler's Tr(J) is asserted to match it AND the eta=0 closed forms
 *      Tr(J(v*)) = n_B, Tr(J((v^{-1})*)) = N (v a unital *-iso => v*(I_N)=I_{n_B},
 *      (v^{-1})*(I_{n_B}) = I_N, traces n_B and N).
 *   T4 cross-check J(v*) vs J(v). Build the FORWARD Choi J(v)[p*N+a, q*N+b] =
 *      v(E_pq)[a,b] (design §1.2, an independent construction) and verify the
 *      relation J(v*)[a*nB+p, b*nB+q] = conj( J(v)[p*N+a, q*N+b] ) implied by
 *      v*(E_ab)[p,q] = conj(v(E_pq)[a,b]) (the HS-adjoint). Two constructions agree
 *      ~1e-12. MUTATION: same index-swap as T2 makes the relation fail -> RED.
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
#include "aic_cstar.h"
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_opspace.h"
#include "aic_test.h"
#include "aic_ucp.h"
#include "test_idemp.h"

#define CPREC 256

static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

/* Suppress -Wunused for the test_idemp.h builders this TU does not use. */
__attribute__((unused))
static void suppress_unused(void)
{
    (void) make_trace_replace;
    (void) make_identity;
    (void) make_compress_idemp;
    (void) make_dephasing;
    (void) make_mixconj;
    (void) make_noiseless_subsystem;
    (void) make_conjugated;
    (void) make_fixed_unitary;
}

/* max_{ij} |M[i,j] - conj(M[j,i])| as a double (Hermiticity residual). */
static double herm_residual(const acb_mat_t M)
{
    slong n = acb_mat_nrows(M);
    double worst = 0.0;
    acb_t t;
    arb_t a;
    acb_init(t);
    arb_init(a);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            acb_conj(t, acb_mat_entry(M, j, i));
            acb_sub(t, acb_mat_entry(M, i, j), t, CPREC);
            acb_abs(a, t, CPREC);
            double d = dd(a);
            if (d > worst) worst = d;
        }
    arb_clear(a);
    acb_clear(t);
    return worst;
}

/* Build the FORWARD Choi J(v) of v: B -> A (design §1.2):
 *   J(v)[p*N + a, q*N + b] = v(E_pq^{(n_B)})[a,b],  p,q in [0,n_B) MAJOR, a,b in
 *   [0,N) MINOR. v(E_pq) = vE[i] when (p,q) = b_unit_pos(i), else 0. Built
 *   INDEPENDENTLY of the assembler via aic_dhom_v_apply on the ambient unit E_pq.
 * `Jv` caller-init'd (n_B*N) x (n_B*N). */
static void build_forward_choi(acb_mat_t Jv, const aic_dhom_v *v)
{
    slong N = v->n, nB = v->B->n_B;
    acb_mat_zero(Jv);
    acb_mat_t Epq, vEpq;
    acb_mat_init(Epq, nB, nB);
    acb_mat_init(vEpq, N, N);
    for (slong p = 0; p < nB; p++)
        for (slong q = 0; q < nB; q++) {
            acb_mat_zero(Epq);
            acb_one(acb_mat_entry(Epq, p, q));
            aic_dhom_v_apply(vEpq, v, Epq, CPREC); /* v(E_pq) (0 off-block) */
            for (slong a = 0; a < N; a++)
                for (slong b = 0; b < N; b++)
                    acb_set(acb_mat_entry(Jv, p * N + a, q * N + b),
                            acb_mat_entry(vEpq, a, b));
        }
    acb_mat_clear(vEpq);
    acb_mat_clear(Epq);
}

/* The fixture: an eta=0 exact-*-iso v (block_cond_exp(d,m)). Caller frees via
 * fixture_clear. */
typedef struct {
    aic_assoc_ecstar ae;
    aic_dhom_B B;
    aic_dhom_v v;
    arb_t iso;
} choi_fixture;

static void fixture_build(choi_fixture *fx, aic_ucp_kraus *phi, const char *name)
{
    aic_assoc_ecstar_from_phi(&fx->ae, phi, CPREC);
    arb_init(fx->iso);
    aic_cstar_build(&fx->B, &fx->v, fx->iso, &fx->ae.A, 0.0, CPREC);
    AIC_CHECK_MSG(dd(fx->iso) < 1e-9, "%s: eta=0 iso_def=%.3e not ~0", name,
                  dd(fx->iso));
    AIC_CHECK_MSG(fx->v.A->dim_A == fx->B.dim_B,
                  "%s: v not bijective dim_A=%ld dim_B=%ld", name,
                  (long) fx->v.A->dim_A, (long) fx->B.dim_B);
}

static void fixture_clear(choi_fixture *fx)
{
    arb_clear(fx->iso);
    aic_dhom_v_clear(&fx->v);
    aic_dhom_B_clear(&fx->B);
    aic_assoc_ecstar_clear(&fx->ae);
}

/* ===================== the four checks on one fixture ====================== */
static void run_one(const char *name, aic_ucp_kraus *phi)
{
    choi_fixture fx;
    fixture_build(&fx, phi, name);
    aic_dhom_v *v = &fx.v;
    slong N = v->n, nB = fx.B.n_B;
    printf("  [%s] N=%ld n_B=%ld dim_B=%ld\n", name, (long) N, (long) nB,
           (long) fx.B.dim_B);

    /* assemble J(v*) (N*n_B square) and J((v^{-1})*) (n_B*N square). */
    acb_mat_t Jvs, Jvis;
    acb_mat_init(Jvs, N * nB, N * nB);
    acb_mat_init(Jvis, nB * N, nB * N);
    aic_opspace_choi_vstar(Jvs, v, CPREC);
    aic_opspace_choi_vinvstar(Jvis, v, CPREC);

    /* ---- T1 Hermiticity ---- */
    double h_vs = herm_residual(Jvs);
    double h_vis = herm_residual(Jvis);
    printf("    T1 Hermiticity: |J(v*) - J(v*)^dag|=%.3e |J((v^-1)*) - ...|=%.3e\n",
           h_vs, h_vis);
    AIC_CHECK_MSG(h_vs < 1e-12, "%s T1: J(v*) not Hermitian (res %.3e) — a dropped "
                  "conj in the assembler", name, h_vs);
    AIC_CHECK_MSG(h_vis < 1e-12, "%s T1: J((v^-1)*) not Hermitian (res %.3e)", name,
                  h_vis);

    /* ---- T2 independent recompute of J(v*) entries ----
     * J(v*)[a*nB+r, b*nB+c] = v*(E_ab^{(N)})[r,c] = conj( v(E_{rc}^{(n_B)})[a,b] ),
     * v(E_rc) via aic_dhom_v_apply (production path, no shared index logic). Sweep a
     * handful of (a,b) and ALL (r,c). */
    {
        acb_mat_t Erc, vErc;
        acb_mat_init(Erc, nB, nB);
        acb_mat_init(vErc, N, N);
        acb_t expect, diff;
        arb_t mag;
        acb_init(expect);
        acb_init(diff);
        arb_init(mag);
        double worst = 0.0;
        slong checked = 0;
        for (slong r = 0; r < nB; r++)
            for (slong c = 0; c < nB; c++) {
                acb_mat_zero(Erc);
                acb_one(acb_mat_entry(Erc, r, c));
                aic_dhom_v_apply(vErc, v, Erc, CPREC); /* v(E_rc) (0 off-block) */
                for (slong a = 0; a < N; a++)
                    for (slong b = 0; b < N; b++) {
                        acb_conj(expect, acb_mat_entry(vErc, a, b));
                        acb_sub(diff, acb_mat_entry(Jvs, a * nB + r, b * nB + c),
                                expect, CPREC);
                        acb_abs(mag, diff, CPREC);
                        double dv = dd(mag);
                        if (dv > worst) worst = dv;
                        checked++;
                    }
            }
        printf("    T2 indep recompute J(v*): max diff=%.3e over %ld entries\n",
               worst, (long) checked);
        AIC_CHECK_MSG(worst < 1e-12, "%s T2: J(v*) disagrees with the independent "
                      "trace-pairing recompute by %.3e — an index swap (a<->r / "
                      "b<->c) in the assembler", name, worst);
        arb_clear(mag);
        acb_clear(diff);
        acb_clear(expect);
        acb_mat_clear(vErc);
        acb_mat_clear(Erc);
    }

    /* ---- T2b independent recompute of J((v^-1)*) entries ----
     * J((v^-1)*)[a*N+r, b*N+s] = (v^-1)*(E_ab^{(n_B)})[r,s] = conj( v^{-1}(E_rs^{(N)}
     * )[a,b] ) (double-adjoint trace pairing). v^{-1}(E_rs^{(N)}) = sum_k <B_k,E_rs>_F
     * vinvB[k] = sum_k conj(B_k[r,s]) vinvB[k], assembled here in a DIFFERENT loop
     * order (output-major) from the assembler (k-major). This is the entrywise tooth
     * for the vinvstar path — T3's trace alone is BLIND to an output-index transpose
     * (r,s)<->(s,r) in B_k (caught here; see report). Sweep all (r,s), a handful (a,b). */
    {
        acb_mat_t *vinvB = NULL;
        slong dimA = 0, nBchk = 0;
        aic_opspace_build_vinv(&vinvB, &dimA, &nBchk, v, CPREC);
        acb_mat_t vinv, term;
        acb_mat_init(vinv, nB, nB);
        acb_mat_init(term, nB, nB);
        acb_t ck, expect, diff;
        arb_t mag;
        acb_init(ck);
        acb_init(expect);
        acb_init(diff);
        arb_init(mag);
        double worst = 0.0;
        slong checked = 0;
        for (slong r = 0; r < N; r++)
            for (slong s = 0; s < N; s++) {
                /* v^{-1}(E_rs^{(N)}) = sum_k conj(B_k[r,s]) vinvB[k] */
                acb_mat_zero(vinv);
                for (slong k = 0; k < dimA; k++) {
                    acb_conj(ck, acb_mat_entry(v->A->B[k], r, s));
                    acb_mat_scalar_mul_acb(term, vinvB[k], ck, CPREC);
                    acb_mat_add(vinv, vinv, term, CPREC);
                }
                for (slong a = 0; a < nB; a++)
                    for (slong b = 0; b < nB; b++) {
                        acb_conj(expect, acb_mat_entry(vinv, a, b));
                        acb_sub(diff, acb_mat_entry(Jvis, a * N + r, b * N + s),
                                expect, CPREC);
                        acb_abs(mag, diff, CPREC);
                        double dv = dd(mag);
                        if (dv > worst) worst = dv;
                        checked++;
                    }
            }
        printf("    T2b indep recompute J((v^-1)*): max diff=%.3e over %ld entries\n",
               worst, (long) checked);
        AIC_CHECK_MSG(worst < 1e-12, "%s T2b: J((v^-1)*) disagrees with the "
                      "independent recompute by %.3e — an index swap / transpose in "
                      "the vinvstar assembler", name, worst);
        arb_clear(mag);
        acb_clear(diff);
        acb_clear(expect);
        acb_clear(ck);
        acb_mat_clear(term);
        acb_mat_clear(vinv);
        aic_opspace_vinv_clear(vinvB, dimA);
    }

    /* ---- T3 trace ---- *
     * GOLDEN-RULE diagonal: Tr(J(f)) = Tr(f(I_in)). The DOUBLE-adjoint identity
     * Tr(f*(I_in)) = Tr(f(I_out)) (from <I_out, f*(I_in)> = <f(I_out), I_in>) gives
     * the INDEPENDENT closed forms (NOT n_B / N — those assume I_A = I_N, which holds
     * for block_cond_exp but NOT for noiseless_subsystem, whose I_A is a rank-6
     * carrier in M_6 with n_B=3; the wrong hardcode mis-fires there, see report):
     *   Tr(J(v*))       = Tr(v(I_{n_B}))   (apply v to B's ambient identity),
     *   Tr(J((v^-1)*))  = Tr(v^{-1}(I_N))  (apply v^{-1} to A's ambient identity).
     * Both computed via the production v-apply / vinvB coords, INDEPENDENT of the
     * assembler's diagonal-summation. On the eta=0 *-iso Tr(v(I_B)) = Tr(I_A). */
    {
        acb_t tr_assemb;
        acb_init(tr_assemb);

        /* assembler Tr(J(v*)) */
        acb_zero(tr_assemb);
        for (slong i = 0; i < N * nB; i++)
            acb_add(tr_assemb, tr_assemb, acb_mat_entry(Jvs, i, i), CPREC);
        double tr_vs = arf_get_d(arb_midref(acb_realref(tr_assemb)), ARF_RND_NEAR);
        double tr_vs_im = arf_get_d(arb_midref(acb_imagref(tr_assemb)), ARF_RND_NEAR);

        /* independent Tr(v(I_{n_B})) via aic_dhom_v_apply on B's ambient identity */
        acb_mat_t Inb, vInb;
        acb_mat_init(Inb, nB, nB);
        acb_mat_init(vInb, N, N);
        acb_mat_one(Inb);
        aic_dhom_v_apply(vInb, v, Inb, CPREC);
        acb_t tr_indep;
        acb_init(tr_indep);
        acb_zero(tr_indep);
        for (slong r = 0; r < N; r++)
            acb_add(tr_indep, tr_indep, acb_mat_entry(vInb, r, r), CPREC);
        double tri_re = arf_get_d(arb_midref(acb_realref(tr_indep)), ARF_RND_NEAR);
        printf("    T3 Tr(J(v*)): assembler=%.9f%+.1ei  Tr(v(I_B))=%.9f\n",
               tr_vs, tr_vs_im, tri_re);
        AIC_CHECK_MSG(fabs(tr_vs - tri_re) < 1e-9 && fabs(tr_vs_im) < 1e-9,
                      "%s T3: assembler Tr(J(v*))=%.9f%+.1ei != Tr(v(I_B))=%.9f "
                      "(GOLDEN-RULE diagonal / adjoint-trace identity broken)",
                      name, tr_vs, tr_vs_im, tri_re);
        acb_mat_clear(vInb);
        acb_mat_clear(Inb);

        /* assembler Tr(J((v^-1)*)) */
        acb_zero(tr_assemb);
        for (slong i = 0; i < nB * N; i++)
            acb_add(tr_assemb, tr_assemb, acb_mat_entry(Jvis, i, i), CPREC);
        double tr_vis = arf_get_d(arb_midref(acb_realref(tr_assemb)), ARF_RND_NEAR);
        double tr_vis_im = arf_get_d(arb_midref(acb_imagref(tr_assemb)), ARF_RND_NEAR);

        /* independent Tr(v^{-1}(I_N)) = Tr( sum_k conj(Tr B_k) vinvB[k] )
         *   = sum_k conj(Tr B_k) Tr(vinvB[k]) (A-coords of I_N are <B_k,I_N>=conj Tr B_k). */
        acb_mat_t *vinvB = NULL;
        slong dimA = 0, nBchk = 0;
        aic_opspace_build_vinv(&vinvB, &dimA, &nBchk, v, CPREC);
        acb_zero(tr_indep);
        acb_t trBk, trvk, ck;
        acb_init(trBk);
        acb_init(trvk);
        acb_init(ck);
        for (slong k = 0; k < dimA; k++) {
            acb_zero(trBk);
            for (slong r = 0; r < N; r++)
                acb_add(trBk, trBk, acb_mat_entry(v->A->B[k], r, r), CPREC);
            acb_zero(trvk);
            for (slong a = 0; a < nB; a++)
                acb_add(trvk, trvk, acb_mat_entry(vinvB[k], a, a), CPREC);
            acb_conj(ck, trBk);          /* <B_k, I_N>_F = conj(Tr B_k)  */
            acb_mul(ck, ck, trvk, CPREC);
            acb_add(tr_indep, tr_indep, ck, CPREC);
        }
        double trvi_re = arf_get_d(arb_midref(acb_realref(tr_indep)), ARF_RND_NEAR);
        printf("    T3 Tr(J((v^-1)*)): assembler=%.9f%+.1ei  Tr(v^{-1}(I_N))=%.9f\n",
               tr_vis, tr_vis_im, trvi_re);
        AIC_CHECK_MSG(fabs(tr_vis - trvi_re) < 1e-9 && fabs(tr_vis_im) < 1e-9,
                      "%s T3: assembler Tr(J((v^-1)*))=%.9f%+.1ei != Tr(v^{-1}(I_N))="
                      "%.9f", name, tr_vis, tr_vis_im, trvi_re);
        /* eta=0 *-iso: v^{-1}(I_A)=I_B, but I_N != I_A in general; the value equals
         * Tr(I_B)=n_B ONLY when I_A = I_N (the unital-into-M_N case, block_cond_exp).
         * Assert that sub-case explicitly when it holds (Tr(v(I_B)) ~ N <=> I_A=I_N). */
        if (fabs(tri_re - (double) N) < 1e-9) {
            AIC_CHECK_MSG(fabs(tr_vis - (double) nB) < 1e-9,
                          "%s T3: I_A=I_N case but Tr(J((v^-1)*))=%.9f != n_B=%ld",
                          name, tr_vis, (long) nB);
            printf("      (I_A = I_N: Tr(J(v*))=N=%ld, Tr(J((v^-1)*))=n_B=%ld)\n",
                   (long) N, (long) nB);
        }
        acb_clear(ck);
        acb_clear(trvk);
        acb_clear(trBk);
        aic_opspace_vinv_clear(vinvB, dimA);
        acb_clear(tr_indep);
        acb_clear(tr_assemb);
    }

    /* ---- T4 cross-check J(v*) vs J(v) ----
     * J(v*)[a*nB+p, b*nB+q] = v*(E_ab)[p,q] = conj(v(E_pq)[a,b]) = conj(J(v)[p*N+a,
     * q*N+b]). Build J(v) independently and verify the relation entrywise. */
    {
        acb_mat_t Jv;
        acb_mat_init(Jv, nB * N, nB * N);
        build_forward_choi(Jv, v);
        acb_t expect, diff;
        arb_t mag;
        acb_init(expect);
        acb_init(diff);
        arb_init(mag);
        double worst = 0.0;
        for (slong a = 0; a < N; a++)
            for (slong b = 0; b < N; b++)
                for (slong p = 0; p < nB; p++)
                    for (slong q = 0; q < nB; q++) {
                        acb_conj(expect, acb_mat_entry(Jv, p * N + a, q * N + b));
                        acb_sub(diff, acb_mat_entry(Jvs, a * nB + p, b * nB + q),
                                expect, CPREC);
                        acb_abs(mag, diff, CPREC);
                        double dv = dd(mag);
                        if (dv > worst) worst = dv;
                    }
        printf("    T4 J(v*) vs conj(J(v)) (HS-adjoint relation): max diff=%.3e\n",
               worst);
        AIC_CHECK_MSG(worst < 1e-12, "%s T4: J(v*) != conj(J(v)) by %.3e — the "
                      "adjoint Choi and the forward Choi disagree (assembler index "
                      "or conj error)", name, worst);
        arb_clear(mag);
        acb_clear(diff);
        acb_clear(expect);
        acb_mat_clear(Jv);
    }

    acb_mat_clear(Jvis);
    acb_mat_clear(Jvs);
    fixture_clear(&fx);
}

int main(void)
{
    printf("=== test_opspace_choi: adjoint Choi assemblers (O2.1, aic-pjr) ===\n");
    aic_ucp_kraus phi;
    make_block_cond_exp(&phi, 4, 2);      /* M_2 (+) M_2: n_B=4, dim_B=8 */
    run_one("block_cond_exp(4,2)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_block_cond_exp(&phi, 6, 2);      /* M_2 (+) M_4: n_B=6, dim_B=20 */
    run_one("block_cond_exp(6,2)", &phi);
    aic_ucp_kraus_clear(&phi);
    make_noiseless_subsystem(&phi, 3, 2); /* M_3 carrier: n_B=3, dim_B=9 */
    run_one("noiseless_subsystem(3,2)", &phi);
    aic_ucp_kraus_clear(&phi);

    /* COMPLEX eta=0 fixture (the conj/index teeth). make_conjugated wraps an exactly-
     * idempotent base by a FIXED COMPLEX unitary V = exp(iH): Phi' = Ad_{V^dag} o Phi
     * o Ad_V is STILL exactly idempotent (Ad_V Ad_{V^dag}=id) but its image algebra is
     * NOT real, so vE / A->B carry genuine imaginary parts — without which the T1/T4
     * conj-drop mutation is INVISIBLE (the real fixtures above have conj == identity).
     * This is the fixture that gives the conj mutation teeth (see report). */
    {
        aic_ucp_kraus base, conj;
        make_block_cond_exp(&base, 4, 2);    /* M_2 (+) M_2, exactly idempotent */
        make_conjugated(&conj, &base, CPREC);
        run_one("conj(block_cond_exp(4,2)) [complex eta=0]", &conj);
        aic_ucp_kraus_clear(&conj);
        aic_ucp_kraus_clear(&base);
    }
    aic_test_report("test_opspace_choi");
    return 0;
}
