/* aic_cstar_extension.c — lem_extension, the inductive Stage-2 growth step of the
 * th_main master loop (bead aic-097, module "cstar_build", Increment 4). See
 * include/aic_cstar.h (Increment 4 banner) for the full why, the six-step route,
 * and the ownership/lifetime contract; design contract docs/research/
 * lem_extension_spec.md (the signature-verified spec). lem_extension is the only
 * consumer of lem_merging's two_block=0 LIVE-off-diagonal shape (FINDINGS §C9).
 *
 * approximate_algebras.tex:1378-1411 — Lemma lem_extension (the proof skeleton):
 *   ".tex:1383: dim S_{P,Q} = n.  .tex:1388: h_{jk} = Ha^Q_{P_j,P_k} with P_1=P,
 *    P_2=Q.  .tex:1391: mu_{11} = lem_approx(h_{11} v), an isomorphism M_n ->
 *    B(S_{P,Q}).  .tex:1404: mu_{11}(A) = U_1 A U_1^dag, U_1 : C^n -> S_{P,Q}
 *    unitary; U_2 : c |-> c Qtilde.  .tex:1405-1410: gamma_{11}(A11)=v(A11),
 *    gamma_{12}(A12)=U_1(A12), gamma_{21}(A21)=(U_1(A21^dag))^dag,
 *    gamma_{22}(A22)=A22 Qtilde.  .tex:1411: v_+ = lem_merging(gamma_{jk})."
 *
 * CONSTRUCTIVE ROUTE (CLAUDE.md Law 3). The proof asserts mu_{11} exists and is
 * inner (U_1 A U_1^dag) by representation theory; in finite dim lem_approx IS the
 * Newton iteration (aic_dhom_approx on the codomain M_n = aic_cstar_matrix_algebra,
 * a GENUINE C* algebra so eps_target=0), and U_1 is recovered constructively from
 * mu_{11}: since mu_{11}(E_l0) = (U_1 e_l)(U_1 e_0)^dag = u_l u_0^dag (E_00 a rank-1
 * projector), u_0 is the top LEFT singular vector of mu_{11}(E_00) = u_0 u_0^dag, and
 * u_l = mu_{11}(E_l0) u_0 / ||u_0||^2. The U_1 columns are S_{P,Q}-COORDINATE vectors
 * (the {C_l} corner basis is Frobenius-orthonormal, so the coordinate norm ~= the
 * lem_PQ_Hilb Hilbert norm up to O(delta+eps), FINDINGS §F-I4-2/§F-I4-4); the physical
 * operator is U1_op[l] = sum_m U1[m,l] C_m.
 *
 * THE Ha-MAP ARGUMENT ORDER (the highest-risk convention, spec §2). aic_corner_ha
 * (Ha, A, Z, P, R, Q) computes Ha^Q_{P,R}(Z) : S_{R,Q} -> S_{P,Q} (d_PQ x d_RQ). For
 * h_{11} = Ha^Q_{P,P} the call is (P, P, Q): output d_PQ x d_PQ (n x n). Feeding the
 * wrong pair (e.g. (Q,P,Q) = Ha^Q_{Q,P}, a 1 x n map) trips the aic_corner_ha
 * dimension assert (the Ha is pre-init'd n x n) — the spec §10.A / test tooth (a).
 *
 * THE STAR IS LOAD-BEARING (FINDINGS §C2/§C8). The output mult_def goes through A's
 * STAR (aic_dhom_defect_sweep); on in-A delta-projections the star-vs-plain gap is a
 * MAGNITUDE gap c = mult_def/eta (the §C8 tooth, in the test).
 *
 * LOC budget (Rule 10): three static helpers (build_h11v, extract_U1,
 * build_gamma_arrays) + the assembly. The A_cod codomain (M_n) is BORROWED by h11v
 * and must outlive its aic_dhom_v_clear (spec §6.3).
 */
#include <assert.h>
#include <complex.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_corner.h"
#include "aic_cstar.h"
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_latd.h"

/* Step 2+3 codomain map h_{11}v : M_n -> B(S_{P,Q}) = M_n (.tex:1390-1391). For
 * each matrix unit E_{lm} of M_n, h11v(E_lm) = Ha^Q_{P,P}(v(E_lm)) (n x n, the
 * d_PQ x d_PQ corner-basis matrix). v->vE[l*n+m] = v(E_lm) in S_P; A_parent carries
 * the corner machinery; A_cod = M_n is the genuine-C* codomain. The (P,P,Q) order
 * is load-bearing (file docstring). */
static void build_h11v(aic_dhom_v *h11v, const aic_dhom_v *v, const acb_mat_t P,
                       const acb_mat_t Q, const aic_ecstar *A_parent, slong n,
                       slong prec)
{
    for (slong l = 0; l < n; l++)
        for (slong m = 0; m < n; m++)
            aic_corner_ha(h11v->vE[l * n + m], A_parent, v->vE[l * n + m],
                          P, P, Q, prec);
}

/* Step 4: U_1 from mu_{11} (.tex:1404). Fills U1_ops[0..n-1] (each n x n, the
 * physical S_{P,Q} operator U_1(e_l) = sum_m U1[m,l] C_PQ[m]). C_PQ[] = the dim_S=n
 * Frobenius-orthonormal corner basis of S_{P,Q}. ASSERTS ||u0||^2 > 0.5 (the rank-1
 * SVD of mu_{11}(E_00) is nondegenerate, Rule 4). */
static void extract_U1(acb_mat_t *U1_ops, const aic_dhom_v *mu11,
                       acb_mat_t *const C_PQ, slong n, slong prec)
{
    /* u0 = top LEFT singular vector of mu_{11}(E_00) = mu11->vE[0] (= u0 u0^dag). */
    double _Complex *m00 = flint_malloc((size_t) (n * n) * sizeof(double _Complex));
    double _Complex *Ul  = flint_malloc((size_t) (n * n) * sizeof(double _Complex));
    double *sv = flint_malloc((size_t) n * sizeof(double));
    aic_latd_from_acb_mat(m00, mu11->vE[0]);
    aic_latd_svd(sv, Ul, NULL, m00, n, n);

    acb_mat_t u0, col, U1coord;     /* coordinate column vectors (length n) */
    acb_mat_init(u0, n, 1);
    acb_mat_init(col, n, 1);
    acb_mat_init(U1coord, n, n);    /* U1[m,l] = m-th coord of U_1(e_l) */
    double nrm2 = 0.0;
    for (slong i = 0; i < n; i++) {
        acb_set_d_d(acb_mat_entry(u0, i, 0), creal(Ul[i * n + 0]),
                    cimag(Ul[i * n + 0]));
        nrm2 += creal(Ul[i * n + 0]) * creal(Ul[i * n + 0]) +
                cimag(Ul[i * n + 0]) * cimag(Ul[i * n + 0]);
    }
    assert(nrm2 > 0.5 &&
           "extract_U1: ||u0||^2 <= 0.5 (mu_{11}(E_00) SVD degenerate)");

    /* U1_col_l = mu_{11}(E_l0) u0 / ||u0||^2, the coordinate column l of U_1. */
    acb_t inv; acb_init(inv); acb_set_d(inv, 1.0 / nrm2);
    for (slong l = 0; l < n; l++) {
        acb_mat_mul(col, mu11->vE[l * n + 0], u0, prec);  /* mu_{11}(E_l0) u0 */
        for (slong m = 0; m < n; m++) {
            acb_mul(acb_mat_entry(U1coord, m, l), acb_mat_entry(col, m, 0), inv,
                    prec);
        }
    }

    /* Physical operator U1_op[l] = sum_m U1[m,l] C_PQ[m] (.tex:1404; §F-I4-4). */
    acb_mat_t term;
    acb_mat_init(term, acb_mat_nrows(C_PQ[0]), acb_mat_ncols(C_PQ[0]));
    for (slong l = 0; l < n; l++) {
        acb_mat_zero(U1_ops[l]);
        for (slong m = 0; m < n; m++) {
            acb_mat_scalar_mul_acb(term, C_PQ[m], acb_mat_entry(U1coord, m, l),
                                   prec);
            acb_mat_add(U1_ops[l], U1_ops[l], term, prec);
        }
    }

    acb_mat_clear(term);
    acb_clear(inv);
    acb_mat_clear(U1coord);
    acb_mat_clear(col);
    acb_mat_clear(u0);
    flint_free(sv);
    flint_free(Ul);
    flint_free(m00);
}

/* Step 5: the four gamma_{jk} arrays (.tex:1405-1410). n1 = n (P block), n2 = 1 (Q
 * block). g11[a*n+b] = v(E_ab) = v->vE[a*n+b]; g12[l] = U1_op[l]; g21[l] =
 * g12[l]^dag; g22[0] = Qtilde. All n x n DEEP COPIES (the lem_merging contract). */
static void build_gamma_arrays(acb_mat_t *g11, acb_mat_t *g12, acb_mat_t *g21,
                               acb_mat_t *g22, const aic_dhom_v *v,
                               const acb_mat_t *U1_ops, const acb_mat_t Qtilde,
                               slong n)
{
    for (slong i = 0; i < n * n; i++) acb_mat_set(g11[i], v->vE[i]);
    for (slong l = 0; l < n; l++) {
        acb_mat_set(g12[l], U1_ops[l]);
        acb_mat_conjugate_transpose(g21[l], U1_ops[l]);
    }
    acb_mat_set(g22[0], Qtilde);
}

void aic_cstar_lem_extension(aic_dhom_B *B_out, aic_dhom_v *v_out,
                             arb_t mult_def, arb_t sigma_min,
                             const aic_dhom_v *v,
                             const acb_mat_t P, const acb_mat_t Q,
                             const aic_ecstar *A_parent, slong prec)
{
    assert(B_out != NULL && v_out != NULL && v != NULL && A_parent != NULL);
    assert(v->A == A_parent && "lem_extension: v->A must be A_parent");
    assert(v->B->num_blocks == 1 && "lem_extension: v->B must be a single block");
    slong n = v->B->d[0];
    slong nn = A_parent->n;
    assert(acb_mat_nrows(P) == nn && acb_mat_nrows(Q) == nn &&
           "lem_extension: P, Q must be A_parent->n x A_parent->n");

    /* Step 1 (.tex:1383): dim S_Q == 1 (the 1d-Q hypothesis, also the Ha-map's
     * precondition) and dim S_{P,Q} == n (the lem_add_dim CONCLUSION). The proof
     * derives the conclusion from the per-j "dim S_{P_j,Q} == 1 for each j" (P_j =
     * v(E_jj)); we assert the CONCLUSION directly, which is the operative guard.
     * The per-j sub-query is INTENTIONALLY NOT run: with P_j = v->vE[jj] = Ptilde_P
     * (an OBLIQUE element of S_P at eta>0), aic_corner_dim_S builds Co_{P_j,Q} whose
     * sgn-basin dispatch (aic_funcalc_int_def_X2 -> aic_mat_opnorm) hits the aic-qgs
     * Gram-Hermiticity FALSE-FAIL on the near-identity (2M-1)^2-I (radii ~1e-3,
     * genuinely IN-basin; SIGABRTs). The cleaner equivalent dim_S(P,Q) (P a genuine
     * Hermitian projection) does NOT trip it. See paper/FINDINGS.md §C5/§C10. */
    assert(aic_corner_dim_S(A_parent, Q, Q, prec) == 1 &&
           "lem_extension: dim S_Q != 1 (Q not one-dimensional)");
    slong dPQ = aic_corner_dim_S(A_parent, P, Q, prec);
    assert(dPQ == n && "lem_extension: dim S_{P,Q} != n (lem_add_dim conclusion; "
                       "either S_{P,Q}=0 hypothesis violated or P not rank-n)");

    /* Steps 2-3: A_cod = genuine M_n; h11v = Ha^Q_{P,P} v; mu_{11} = lem_approx. */
    aic_ecstar A_cod;
    aic_cstar_matrix_algebra(&A_cod, n, prec);
    aic_dhom_B B_dom;
    slong dims[1] = {n};
    aic_dhom_B_init(&B_dom, dims, 1);
    aic_dhom_v h11v;
    aic_dhom_v_init(&h11v, &B_dom, &A_cod);
    build_h11v(&h11v, v, P, Q, A_parent, n, prec);
    aic_dhom_approx(&h11v, /*eps_target=*/0.0, /*tol_abs=*/1e-12,
                    /*unit_tol=*/1e-10, /*max_steps=*/30, prec, NULL);

    /* Step 4: corner basis {C_PQ} of S_{P,Q}; U1_ops from mu_{11}. */
    acb_mat_t CoPQ;
    acb_mat_init(CoPQ, A_parent->dim_A, A_parent->dim_A);
    aic_corner_Co(CoPQ, A_parent, P, Q, prec);
    acb_mat_t *C_PQ;
    slong dim_S;
    aic_corner_extract(&C_PQ, &dim_S, CoPQ, A_parent, prec);
    assert(dim_S == n && "lem_extension: extracted dim S_{P,Q} != n");
    acb_mat_t *U1_ops = flint_malloc((size_t) n * sizeof(acb_mat_t));
    for (slong l = 0; l < n; l++) acb_mat_init(U1_ops[l], nn, nn);
    extract_U1(U1_ops, &h11v, C_PQ, n, prec);

    /* Step 5: Qtilde = Co_Q(Q); the four gamma arrays (n1=n, n2=1). */
    acb_mat_t Qtilde;
    acb_mat_init(Qtilde, nn, nn);
    aic_corner_Ptilde(Qtilde, A_parent, Q, prec);
    acb_mat_t *g11 = flint_malloc((size_t) (n * n) * sizeof(acb_mat_t));
    acb_mat_t *g12 = flint_malloc((size_t) n * sizeof(acb_mat_t));
    acb_mat_t *g21 = flint_malloc((size_t) n * sizeof(acb_mat_t));
    acb_mat_t g22[1];
    for (slong i = 0; i < n * n; i++) acb_mat_init(g11[i], nn, nn);
    for (slong l = 0; l < n; l++) { acb_mat_init(g12[l], nn, nn);
                                    acb_mat_init(g21[l], nn, nn); }
    acb_mat_init(g22[0], nn, nn);
    build_gamma_arrays(g11, g12, g21, g22, v, (const acb_mat_t *) U1_ops, Qtilde, n);

    /* Step 6 (.tex:1411): v_+ = lem_merging(two_block=0, n1=n, n2=1). The lem_merging
     * output guards (mult_def, sigma_min) certify the merged inclusion; the per-block
     * merging conditions follow automatically (F-I4-1). */
    aic_cstar_lem_merging(B_out, v_out, mult_def, sigma_min, /*merge_cond_max=*/NULL,
                          /*two_block=*/0, /*n1=*/n, /*n2=*/1,
                          (const acb_mat_t *) g11, (const acb_mat_t *) g12,
                          (const acb_mat_t *) g21, (const acb_mat_t *) g22,
                          A_parent, P, Q, prec);

    /* Cleanup. A_cod is borrowed by h11v; clear h11v first. */
    acb_mat_clear(g22[0]);
    for (slong l = 0; l < n; l++) { acb_mat_clear(g21[l]); acb_mat_clear(g12[l]); }
    for (slong i = 0; i < n * n; i++) acb_mat_clear(g11[i]);
    flint_free(g21); flint_free(g12); flint_free(g11);
    acb_mat_clear(Qtilde);
    for (slong l = 0; l < n; l++) acb_mat_clear(U1_ops[l]);
    flint_free(U1_ops);
    for (slong m = 0; m < dim_S; m++) acb_mat_clear(C_PQ[m]);
    if (C_PQ) flint_free(C_PQ);
    acb_mat_clear(CoPQ);
    aic_dhom_v_clear(&h11v);
    aic_dhom_B_clear(&B_dom);
    aic_ecstar_clear(&A_cod);
}
