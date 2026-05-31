/* aic_opspace_choi.c — the ADJOINT Choi assemblers for the opspace O2 certified
 * cb-norm bound (bead aic-pjr). Builds the Convention-A Choi matrices of v* and
 * (v^{-1})* in the arb (acb_mat) path, which the diamond-norm SDP certifier (a
 * later task) feeds to the Watrous SDP-restoration machinery. See
 * include/aic_opspace.h for the contract.
 *
 * WHY THE ADJOINT, AND WHY DIRECTLY. We certify the cb-SPECTRAL norm of the
 * th_main_ext isomorphism v: B -> A (.tex:1538-1540) via the adjoint duality
 *   ||v||_cb = ||v*||_⋄,   ||v^{-1}||_cb = ||(v^{-1})*||_⋄
 * (the cb-spectral / cb-trace bridge |||u|||_∞ = |||u*|||_1, Watrous 2009;
 * docs/research/opspace_o2_design.md §0.5 PINNED CONVENTION). The diamond norm is
 * the Watrous SDP the project already solves. We build the adjoint Choi DIRECTLY in
 * the INPUT-MAJOR Convention A — NOT by transposing J(v): the full transpose of J(v)
 * keeps the [out,in] block layout, i.e. it is v* in the OUTPUT-major convention and
 * would need swapped SDP dims; building v*'s Choi directly in input-major form is the
 * clean object the rect SDP expects (§0.5, supersedes §1.4's "J(v*) = J(v)^T").
 *
 * THE GOLDEN RULE (design §0.5). For f: M_in -> M_out the Convention-A Choi is
 *   J(f)[s*out + i, t*out + j] = f(E_st)[i,j]   (INPUT s,t MAJOR/stride out;
 *                                                OUTPUT i,j MINOR).
 * The SDP is then diamond_*_rect(J, d_maj = in, d_min = out): raw optval = ||f||_⋄
 * EXACTLY, normalization factor 1, NO 2/n. (The pin was set empirically by
 * tools/probe_o2_pin2.jl against an asymmetric NON-CP closed-form truth + the eta=0
 * complete-isometry oracle ||v||_cb = ||v^{-1}||_cb = 1.)
 *
 * THE HS-ADJOINT FORMULAS (the load-bearing index/conjugation content — a swapped
 * index or a dropped conj silently corrupts the SDP input, FINDINGS-class).
 *   v*: M_N -> M_{n_B}.  v(X) = sum_i x_i vE[i] is linear in the B matrix-unit
 *     coords x_i. The HS-adjoint w.r.t. <A,B>_HS = Tr(A^dag B): for Y in M_N,
 *       <E_i, v*(Y)>_HS = <v(E_i), Y>_HS = Tr(vE[i]^dag Y),
 *     and since {E_i} is the (Frobenius-ON) matrix-unit basis of B, the i-th coord of
 *     v*(Y) is Tr(vE[i]^dag Y). For Y = E_ab^{(N)} (single 1 at (a,b)),
 *       Tr(vE[i]^dag E_ab) = (vE[i]^dag)[b,a] = conj(vE[i][a,b]),
 *     hence  v*(E_ab^{(N)}) = sum_i conj(vE[i][a,b]) E_i.   ((a,b) INPUT, E_i OUTPUT.)
 *   (v^{-1})*: M_{n_B} -> M_N. Same derivation with v^{-1}(B_k) = vinvB[k] (each
 *     n_B x n_B in B) for A's Frobenius-ON basis {B_k}: v^{-1}(X) = sum_k <B_k,X>_F B
 *     -wait, v^{-1} maps A -> B and is linear in the A-coords. Its HS-adjoint on the
 *     B-ambient matrix unit E_ab^{(n_B)} is
 *       (v^{-1})*(E_ab^{(n_B)}) = sum_k conj(vinvB[k][a,b]) B_k.  ((a,b) INPUT in
 *     M_{n_B}, B_k OUTPUT in M_N.)
 *
 * The two routines are the SAME shape (zero the Choi; loop the linear map's defining
 * operators; place conj(coord) at the input-major / output-minor slot), differing in
 * the (in,out) dims and which operator family carries the output.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>

#include "aic_opspace.h"
#include "aic_opspace_map_internal.h"

/* J(v*) for v*: M_N -> M_{n_B}, d_maj = in = N, d_min = out = n_B.
 *   v*(E_ab^{(N)}) = sum_i conj(vE[i][a,b]) E_i,  E_i = single 1 at
 *   aic_opspace_b_unit_pos(B,i) = (r,c).
 * GOLDEN RULE slot: Jvs[a*n_B + r, b*n_B + c] = conj(vE[i][a,b]) (a,b major / (r,c)
 * minor). Each matrix unit E_i hits ONE (r,c), so no accumulation is needed (distinct
 * i -> distinct (r,c)); we set rather than add. */
void aic_opspace_choi_vstar(acb_mat_t Jvs, const aic_dhom_v *v, slong prec)
{
    const aic_dhom_B *B = v->B;
    slong N = v->n, nB = B->n_B, dimB = B->dim_B, M = N * nB;
    (void) prec; /* pure copy + conjugation; no rounding-prec arithmetic */
    assert(acb_mat_nrows(Jvs) == M && acb_mat_ncols(Jvs) == M
           && "choi_vstar: Jvs must be (N*n_B) x (N*n_B)");
    acb_mat_zero(Jvs);
    for (slong i = 0; i < dimB; i++) {
        slong r, c;
        aic_opspace_b_unit_pos(B, i, &r, &c);
        for (slong a = 0; a < N; a++)
            for (slong b = 0; b < N; b++)
                acb_conj(acb_mat_entry(Jvs, a * nB + r, b * nB + c),
                         acb_mat_entry(v->vE[i], a, b));
    }
}

/* J((v^{-1})*) for (v^{-1})*: M_{n_B} -> M_N, d_maj = in = n_B, d_min = out = N.
 *   (v^{-1})*(E_ab^{(n_B)}) = sum_k conj(vinvB[k][a,b]) B_k, vinvB[k] = v^{-1}(B_k)
 *   (each n_B x n_B), B_k = v->A->B[k] (each N x N).
 * GOLDEN RULE slot: Jvis[a*N + r, b*N + s] += conj(vinvB[k][a,b]) * B_k[r,s], summed
 * over k (a,b in [0,n_B) major; (r,s) in [0,N) minor). Here B_k is a DENSE N x N
 * output operator, so each (a,b) accumulates over k AND over the (r,s) support of B_k
 * — hence add, not set. */
void aic_opspace_choi_vinvstar(acb_mat_t Jvis, const aic_dhom_v *v, slong prec)
{
    const aic_ecstar *A = v->A;
    slong N = v->n, nB = v->B->n_B, M = nB * N;
    assert(A->dim_A == v->B->dim_B
           && "choi_vinvstar: v not bijective (dim_A != dim_B)");
    assert(acb_mat_nrows(Jvis) == M && acb_mat_ncols(Jvis) == M
           && "choi_vinvstar: Jvis must be (n_B*N) x (n_B*N)");

    acb_mat_t *vinvB = NULL;
    slong dimA = 0, nB_chk = 0;
    aic_opspace_build_vinv(&vinvB, &dimA, &nB_chk, v, prec);
    assert(nB_chk == nB && "choi_vinvstar: build_vinv n_B mismatch");

    acb_mat_zero(Jvis);
    acb_t coef, t;
    acb_init(coef);
    acb_init(t);
    for (slong k = 0; k < dimA; k++) {
        /* B_k = A->B[k] (N x N output operator); vinvB[k] = v^{-1}(B_k) (n_B x n_B
         * input-coord carrier). */
        for (slong a = 0; a < nB; a++)
            for (slong b = 0; b < nB; b++) {
                acb_conj(coef, acb_mat_entry(vinvB[k], a, b));
                if (acb_is_zero(coef)) continue;
                for (slong r = 0; r < N; r++)
                    for (slong s = 0; s < N; s++) {
                        acb_mul(t, coef, acb_mat_entry(A->B[k], r, s), prec);
                        acb_add(acb_mat_entry(Jvis, a * N + r, b * N + s),
                                acb_mat_entry(Jvis, a * N + r, b * N + s), t, prec);
                    }
            }
    }
    acb_clear(t);
    acb_clear(coef);
    aic_opspace_vinv_clear(vinvB, dimA);
}
