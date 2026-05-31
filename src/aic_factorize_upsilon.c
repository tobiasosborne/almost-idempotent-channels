/* aic_factorize_upsilon.c — the UCP DECODE map Upsilon, PART 1: the per-block
 * Choi data W_j, R_j, C_j, xi_j and the lem_RC certification (bead aic-tff,
 * module "factorize", Increment F3). EXTENDS aic_factorize.c (F1) +
 * aic_factorize_delta.c (F2). See include/aic_factorize.h (the F3 banner) and
 * docs/research/factorize_f3_synthesis.md (the LOCKED decisions D1-D6).
 *
 * approximate_algebras.tex:2831-2857 — th_factorization Step 5 (verbatim):
 *   Choi_Delta (.tex:2834):
 *     Delta(X) = sum_{j=1}^m W_j^dag (X_j (x) 1_{E_j}) W_j,
 *     W_j: H -> L_j (x) E_j,  sum_j W_j^dag W_j = 1_H.
 *   lem_RC (.tex:2840-2847):
 *     R_j = sum_s p_{js} (U_{js}^dag (x) 1_{E_j}) W_j W_j^dag (U_{js} (x) 1_{E_j})
 *     has the form 1_{L_j} (x) C_j for some C_j in B(E_j). Furthermore,
 *     1 - O(eta) <= ||C_j|| <= 1.
 *
 * THE LOCKED CONVENTIONS (the bug-class that bit cbnorm / opspace O2 — tensor /
 * partial-trace DIRECTION; pinned by the F3 standalone probes test_ptrace_pin /
 * test_centrality_pin in tests/test_factorize.c BEFORE this code was written):
 *
 *   D1 (the #1 risk). R_j in B(L_j (x) E_j) with L_j the LEFT/MAJOR factor (forced
 *     by the U_{js} (x) 1_{E_j} Kroneckers, left-major per aic_mat.h). lem_RC(i)
 *     R_j = 1_{L_j} (x) C_j, so C_j = (1/d_{L_j}) Tr_{L_j} R_j =
 *     (1/d_j) aic_mat_partial_trace_left(R_j, a=d_j, b=e_j) — trace the LEFT factor,
 *     keep the e_j x e_j RIGHT (E_j). (On R_j = 1_{L_j} (x) C_j, Tr_left =
 *     Tr(1_{L_j}) C_j = d_j C_j, so the 1/d_j recovers C_j exactly.)
 *
 *   D2. W_j is built L_j-MAJOR DIRECTLY: the Convention-A Choi of the j-th CP term
 *     Delta_j: M_{d_j} -> B(H) of the UCP Delta (F2) is extracted to Kraus
 *     {D_{j,c}: H -> L_j} (each d_j x N), then stacked
 *       W_j[a*e_j + c, p] = D_{j,c}[a, p]   (a in L_j, c in E_j, p in H).
 *     Do NOT route through aic_ucp_kraus_to_stinespring: its column-stack packs the
 *     ancilla (E_j) index LEFT, reversing the needed L_j-LEFT ordering (D2 warning).
 *     GAUGE: W_j is unique only up to 1_{L_j} (x) u_j on E_j; we ASSERT only the
 *     gauge-invariants (sigma_max(C_j), ||R_j - 1(x)C_j||, the final Upsilon').
 *
 *   D3. The U_{js} are the PER-BLOCK d_j x d_j Paulis S_ab = X^a Z^b (d_j^2 of
 *     them, weight d_j^{-2}; aic_dhom_pauli) on L_j = C^{d_j} — NOT F2's whole-B
 *     joint design (which at m>=2 mixes blocks so R_j != 1(x)C_j). They coincide
 *     only at m=1.
 *
 * VERIFICATION THAT THE STACK REPRODUCES Choi_Delta. With the L_j-major stack,
 *   [W_j^dag (X (x) 1_{E_j}) W_j]_{p,q}
 *     = sum_{a,b,c} conj(D_{j,c}[a,p]) X[a,b] D_{j,c}[b,q]
 *     = sum_c (D_{j,c}^dag X D_{j,c})[p,q] = Delta_j(X)[p,q],
 * exactly the Heisenberg action Delta_j(X) = sum_c D_{j,c}^dag X D_{j,c} that the
 * Convention-A Choi extraction (dim_K=d_j, dim_H=N) produces.
 */
#include <assert.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_dhom.h"
#include "aic_factorize.h"
#include "aic_mat.h"
#include "aic_ucp.h"

/* The PSD-CONE EXTRACTION TOLERANCE for the almost-idempotent W_j Choi
 * (paper/FINDINGS.md §C14). The UCP Delta of th_factorization is CP only to
 * O(eta^2): its per-block Choi Cj_choi has a small negative eigenvalue of order
 * O(eta^2) (measured worst -2.5e-6 at eta=2.3e-2, scaling as ~0.005*eta^2),
 * because Delta'(.tex:2788) is manifestly CP only when Dtilde=v is an EXACT
 * *-hom, and v is an O(eta)-iso (the star != ordinary product). We admit that
 * cone defect (project onto the PSD cone) but ABORT a genuine non-CP input:
 *
 *   AIC_FACTORIZE_CONE_NEG_TOL = 1e-3 — the per-eigenvalue abort floor passed to
 *     aic_ucp_choi_to_kraus_latd_tol. ABSOLUTE (not relative): chosen because the
 *     measured O(eta^2) defect ~1e-6 (at eta~1e-2) sits FAR below 1e-3, while a
 *     genuine O(1) non-CP eigenvalue, OR an unexpectedly-O(eta) negative one
 *     (eta~1e-2 >> 1e-3, signalling the construction is MORE broken than the
 *     O(eta)-CP theory predicts), trips the abort. "Measure then pin": 1e-3 is
 *     ~400x above the worst measured defect, ~10x below the smallest eta scale.
 */
#define AIC_FACTORIZE_CONE_NEG_TOL 1e-3

void aic_factorize_upsilon_Wj(acb_mat_t Wj, slong *e_j_out, double *clipped_out,
                              const aic_factorize *F, slong j, slong prec)
{
    assert(F != NULL && e_j_out != NULL);
    assert(F->delta_ready &&
           "upsilon_Wj: aic_factorize_delta_build (F2) must run first");
    const aic_dhom_B *B = F->v->B;
    assert(j >= 0 && j < B->num_blocks && "upsilon_Wj: block index out of range");
    slong N = F->N, dj = B->d[j];

    /* Convention-A Choi of Delta_j: M_{d_j} -> B(H) (the j-th Choi_Delta term),
     * Cj_choi[a*N+p, b*N+q] = Delta(E^{(j)}_ab)[p,q] (the F2 block-Choi LAYOUT, but
     * with the UCP Delta, not Delta'). dim_K = d_j (L_j, major), dim_H = N. */
    acb_mat_t Cj_choi, Eab, DE;
    acb_mat_init(Cj_choi, dj * N, dj * N);
    acb_mat_init(Eab, B->n_B, B->n_B);
    acb_mat_init(DE, N, N);
    acb_mat_zero(Cj_choi);
    for (slong a = 0; a < dj; a++)
        for (slong b = 0; b < dj; b++) {
            slong idx = aic_dhom_B_index(B, j, a, b);
            aic_dhom_B_matunit(Eab, B, idx);            /* E^{(j)}_ab in B */
            aic_factorize_delta(DE, F, Eab, prec);      /* Delta(E^{(j)}_ab) (UCP) */
            for (slong p = 0; p < N; p++)
                for (slong q = 0; q < N; q++)
                    acb_set(acb_mat_entry(Cj_choi, a * N + p, b * N + q),
                            acb_mat_entry(DE, p, q));
        }

    /* Extract Kraus {D_{j,c}: H -> L_j} (each d_j x N; e_j = rank). Delta_j's Choi
     * is degenerate -> the double path (aic-w4o.1 precedent). TOLERANCE-AWARE: the
     * UCP Delta is CP only to O(eta^2), so Cj_choi has a small O(eta^2) negative
     * eigenvalue (paper/FINDINGS.md §C14); the _tol extraction projects it onto the
     * PSD cone (drops the cone defect) and returns the clipped mass, but STILL
     * aborts on a genuine <= -AIC_FACTORIZE_CONE_NEG_TOL eigenvalue (Rule 4). */
    aic_ucp_kraus phi_j;
    double clipped = 0.0;
    aic_ucp_choi_to_kraus_latd_tol(&phi_j, Cj_choi, dj, N,
                                   AIC_FACTORIZE_CONE_NEG_TOL, &clipped);
    if (clipped_out) *clipped_out = clipped;
    slong e_j = phi_j.r;
    assert(e_j >= 1 && "upsilon_Wj: empty Stinespring (Delta_j is zero?)");

    /* W_j[a*e_j + c, p] = D_{j,c}[a, p] (L_j-major; D2). */
    acb_mat_init(Wj, dj * e_j, N);
    acb_mat_zero(Wj);
    for (slong c = 0; c < e_j; c++)
        for (slong a = 0; a < dj; a++)
            for (slong p = 0; p < N; p++)
                acb_set(acb_mat_entry(Wj, a * e_j + c, p),
                        acb_mat_entry(phi_j.K[c], a, p));

    *e_j_out = e_j;
    aic_ucp_kraus_clear(&phi_j);
    acb_mat_clear(DE); acb_mat_clear(Eab); acb_mat_clear(Cj_choi);
}

void aic_factorize_upsilon_Rj(acb_mat_t Rj, const aic_factorize *F, slong j,
                              const acb_mat_t Wj, slong e_j, slong prec)
{
    assert(F != NULL);
    const aic_dhom_B *B = F->v->B;
    slong N = F->N, dj = B->d[j];
    slong dim = dj * e_j;                       /* dim(L_j (x) E_j) */
    assert(acb_mat_nrows(Rj) == dim && acb_mat_ncols(Rj) == dim &&
           "upsilon_Rj: Rj must be (d_j*e_j) x (d_j*e_j)");
    assert(acb_mat_nrows(Wj) == dim && acb_mat_ncols(Wj) == N &&
           "upsilon_Rj: Wj must be (d_j*e_j) x N");

    /* WWdag = W_j W_j^dag in B(L_j (x) E_j). */
    acb_mat_t Wd, WWd, IE, U, Ud, UI, UdI, tmp, term;
    acb_mat_init(Wd, N, dim);
    acb_mat_init(WWd, dim, dim);
    acb_mat_conjugate_transpose(Wd, Wj);
    acb_mat_mul(WWd, Wj, Wd, prec);

    acb_mat_init(IE, e_j, e_j);
    acb_mat_one(IE);
    acb_mat_init(U, dj, dj);
    acb_mat_init(Ud, dj, dj);
    acb_mat_init(UI, dim, dim);                 /* U_{js} (x) 1_{E_j}      */
    acb_mat_init(UdI, dim, dim);                /* U_{js}^dag (x) 1_{E_j}  */
    acb_mat_init(tmp, dim, dim);
    acb_mat_init(term, dim, dim);

    /* p_{js} = 1/d_j^2 (uniform over the d_j^2 per-block Paulis; D3). */
    arb_t ps;
    arb_init(ps);
    arb_set_si(ps, 1);
    arb_div_si(ps, ps, dj * dj, prec);

    acb_mat_zero(Rj);
    for (slong sa = 0; sa < dj; sa++)
        for (slong sb = 0; sb < dj; sb++) {
            aic_dhom_pauli(U, dj, sa, sb, prec);            /* U_{js} = S_ab on L_j */
            acb_mat_conjugate_transpose(Ud, U);
            aic_mat_kronecker(UI, U, IE, prec);             /* U (x) 1_{E_j}     */
            aic_mat_kronecker(UdI, Ud, IE, prec);           /* U^dag (x) 1_{E_j} */
            /* term = (U^dag (x) 1) WWdag (U (x) 1) */
            acb_mat_mul(tmp, UdI, WWd, prec);
            acb_mat_mul(term, tmp, UI, prec);
            acb_mat_scalar_mul_arb(term, term, ps, prec);
            acb_mat_add(Rj, Rj, term, prec);
        }

    arb_clear(ps);
    acb_mat_clear(term); acb_mat_clear(tmp); acb_mat_clear(UdI); acb_mat_clear(UI);
    acb_mat_clear(Ud); acb_mat_clear(U); acb_mat_clear(IE);
    acb_mat_clear(WWd); acb_mat_clear(Wd);
}

void aic_factorize_upsilon_Cj(acb_mat_t Cj, const acb_mat_t Rj, slong d_j,
                              slong e_j, slong prec)
{
    assert(acb_mat_nrows(Cj) == e_j && acb_mat_ncols(Cj) == e_j &&
           "upsilon_Cj: Cj must be e_j x e_j (E_j)");
    /* C_j = (1/d_{L_j}) Tr_{L_j} R_j (D1): trace the LEFT factor L_j (a=d_j),
     * keep the e_j x e_j RIGHT (E_j). */
    aic_mat_partial_trace_left(Cj, Rj, d_j, e_j, prec);
    arb_t inv;
    arb_init(inv);
    arb_set_si(inv, 1);
    arb_div_si(inv, inv, d_j, prec);
    acb_mat_scalar_mul_arb(Cj, Cj, inv, prec);
    arb_clear(inv);
}
