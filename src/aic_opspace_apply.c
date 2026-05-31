/* aic_opspace_apply.c — the BLOCK-AMPLIATION primitives of the opspace module
 * (bead aic-zwo, increment O1): the ampliated map 1_{M_n} (x) f, its HS adjoint
 * 1_{M_n} (x) f^dag, and the project-onto-M_n(x)U operator. Split from the HOPM
 * kernel (aic_opspace_ampliate.c) to keep each file under the ~200 LOC core limit
 * (CLAUDE.md Rule 10). See include/aic_opspace.h for the contract and FINDINGS
 * §C12 (operator-norm not Frobenius); the map abstraction (the `opmap`) is in
 * src/aic_opspace_hopm.h.
 *
 * THESE THREE FUNCTIONS ARE THE MATHEMATICAL CONTENT OF th_main_ext
 * (approximate_algebras.tex:1447-1561, §10): the block-diagonal-INDEPENDENT
 * ampliation structure. A hostile review (2026-05-31, FINDINGS §C12) crippled
 * aic_opmap_apply_fn to process ONLY the (0,0) block and ALL 72 opspace checks
 * stayed GREEN — the real fixtures give a level-INDEPENDENT a_n ~ 1.0005, so a
 * structural bug that preserves the stretch RATIO was invisible. The fix is the
 * structural-ampliation cross-check tooth (test_opspace.c test_struct_ampliation),
 * which pins aic_opmap_apply_fn against an INDEPENDENT explicit Kronecker block
 * assembly at n>=2 with genuinely nonzero off-diagonal blocks: that cripple now
 * makes the tooth RED. Hence apply_fn is PUBLIC (aic_opmap_apply_fn), so the tooth
 * drives the SAME code path the HOPM kernel uses.
 *
 * THE AMPLIATION (def:opspace .tex:1453-1464). For a level-1 linear map
 * f: U -> V (the `opmap`: ON bases {U[a]} of U, {W[b]} of V, coordinate matrix
 * F[b,a] = <W[b], f(U[a])>_F), 1_{M_n} (x) f maps X in M_n (x) U (an (n*sU) x
 * (n*sU) block matrix, block (k,l) = X_{kl} in U) to the (n*sV) x (n*sV) matrix
 * with block (k,l) = f(X_{kl}). EVERY (k,l) block is mapped, NOT just (0,0). The
 * HS adjoint 1_{M_n} (x) f^dag maps blockwise by f^dag.
 */
#include "aic_opspace_hopm.h"

/* <P, Q>_F = sum conj(P[p]) Q[p] over s*s entries (row-major). */
static cx frob_ip(const cx *P, const cx *Q, slong ss)
{
    cx s = 0.0;
    for (slong p = 0; p < ss; p++) s += conj(P[p]) * Q[p];
    return s;
}

/* Apply f_n = 1_{M_n} (x) f: X (n*sU square) -> Y (n*sV square). cU = n*sU,
 * cV = n*sV. For EACH block (k,l): extract X_{kl} (sU x sU), its coords
 * <U[a],X_{kl}>, build f(X_{kl}) (sV x sV) into block (k,l) of Y. Scratch xb
 * (sU*sU), yb (sV*sV) caller-provided. (The (k,l) loop covers ALL n^2 blocks —
 * FINDINGS §C12 structural landmine: a (0,0)-only cripple here is caught by the
 * structural-ampliation tooth.) */
void aic_opmap_apply_fn(cx *Y, const opmap *m, const cx *X, slong n, cx *xb,
                        cx *yb)
{
    slong sU = m->sU, sV = m->sV, cU = n * sU, cV = n * sV;
    for (slong p = 0; p < cV * cV; p++) Y[p] = 0.0;
    for (slong k = 0; k < n; k++)
        for (slong l = 0; l < n; l++) {
            for (slong i = 0; i < sU; i++)
                for (slong j = 0; j < sU; j++)
                    xb[i * sU + j] = X[(k * sU + i) * cU + (l * sU + j)];
            for (slong p = 0; p < sV * sV; p++) yb[p] = 0.0;
            for (slong b = 0; b < m->dV; b++) {
                cx coef = 0.0;
                for (slong a = 0; a < m->dU; a++)
                    coef += m->F[b * m->dU + a] * frob_ip(m->U[a], xb, sU * sU);
                for (slong p = 0; p < sV * sV; p++) yb[p] += coef * m->W[b][p];
            }
            for (slong i = 0; i < sV; i++)
                for (slong j = 0; j < sV; j++)
                    Y[(k * sV + i) * cV + (l * sV + j)] = yb[i * sV + j];
        }
}

/* Apply f_n^dag = 1_{M_n} (x) f^dag: P (n*sV square) -> G (n*sU square), block
 * (k,l) = f^dag(P_{kl}), which lands in U. Scratch pb (sV*sV), gb (sU*sU). */
void aic_opmap_apply_fn_dag(cx *G, const opmap *m, const cx *P, slong n, cx *pb,
                            cx *gb)
{
    slong sU = m->sU, sV = m->sV, cU = n * sU, cV = n * sV;
    for (slong p = 0; p < cU * cU; p++) G[p] = 0.0;
    for (slong k = 0; k < n; k++)
        for (slong l = 0; l < n; l++) {
            for (slong i = 0; i < sV; i++)
                for (slong j = 0; j < sV; j++)
                    pb[i * sV + j] = P[(k * sV + i) * cV + (l * sV + j)];
            for (slong p = 0; p < sU * sU; p++) gb[p] = 0.0;
            for (slong a = 0; a < m->dU; a++) {
                cx coef = 0.0;
                for (slong b = 0; b < m->dV; b++)
                    coef += conj(m->F[b * m->dU + a]) * frob_ip(m->W[b], pb, sV * sV);
                for (slong p = 0; p < sU * sU; p++) gb[p] += coef * m->U[a][p];
            }
            for (slong i = 0; i < sU; i++)
                for (slong j = 0; j < sU; j++)
                    G[(k * sU + i) * cU + (l * sU + j)] = gb[i * sU + j];
        }
}

/* Project the ambient polar factor Pol (cU x cU) onto M_n (x) U IN PLACE: block
 * (k,l) <- sum_a <U[a], Pol_{kl}>_F U[a]. Scratch pb,gb (sU*sU). */
void aic_opmap_project_amp(cx *Pol, const opmap *m, slong n, cx *pb, cx *gb)
{
    slong sU = m->sU, cU = n * sU;
    for (slong k = 0; k < n; k++)
        for (slong l = 0; l < n; l++) {
            for (slong i = 0; i < sU; i++)
                for (slong j = 0; j < sU; j++)
                    pb[i * sU + j] = Pol[(k * sU + i) * cU + (l * sU + j)];
            for (slong p = 0; p < sU * sU; p++) gb[p] = 0.0;
            for (slong a = 0; a < m->dU; a++) {
                cx c = frob_ip(m->U[a], pb, sU * sU);
                for (slong p = 0; p < sU * sU; p++) gb[p] += c * m->U[a][p];
            }
            for (slong i = 0; i < sU; i++)
                for (slong j = 0; j < sU; j++)
                    Pol[(k * sU + i) * cU + (l * sU + j)] = gb[i * sU + j];
        }
}
