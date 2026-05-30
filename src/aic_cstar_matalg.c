/* aic_cstar_matalg.c — the genuine M_d matrix-algebra ecstar (bead aic-097,
 * module "cstar_build", Increment 1, piece B). Wraps M_d = B(C^d) as a GENUINE
 * finite-dimensional C* algebra in the aic_ecstar verification data model. See
 * include/aic_cstar.h; design contract docs/research/cstar_build_design.md §4.4
 * Step 2 + §G2 (the B(S_{P,Q}) ~= M_n codomain of lem_extension).
 *
 * approximate_algebras.tex:1371-1374 — the standard basis of M_n:
 *   "The algebra Ma{n} has a standard basis {E_{jk}: j,k=1,...,n} with the
 *    relations E_{jk}E_{lm} = delta_{kl} E_{jm}, E_{jk}^dag = E_{kj},
 *    I = sum_j E_{jj}."
 * Each matrix unit E_{lm} (a single 1 at row l, column m) is Frobenius-normalized:
 * <E_{lm},E_{lm}>_F = 1, and distinct units are Frobenius-orthogonal — so the d^2
 * matrix units are an orthonormal operator basis, exactly what aic_ecstar wants.
 *
 * THE STAR IS THE GENUINE PRODUCT. For a genuine C* algebra the Choi-Effros
 * regularizer Phi is the IDENTITY map, so
 *   X * Y = Phi(XY) = XY = plain matrix multiplication,
 * the honest product of M_d. We realize Phi = id through the star_phi seam with a
 * thunk that just COPIES its argument (aic_ecstar_star forms XY then calls the
 * thunk on it). Every eps-C* axiom defect (associativity, submultiplicativity, the
 * C* identity, the unit laws, the involution) is then machine-zero — the cleanest
 * genuine-C* oracle in the project (tests/test_cstar.c T4). phi stays NULL (no
 * Kraus map needed) and star_ctx is NULL (the thunk needs no data), mirroring
 * aic_assoc_ecstar's borrowed-star convention; aic_ecstar_clear frees only the
 * basis.
 */
#include <assert.h>

#include <flint/acb_mat.h>

#include "aic_cstar.h"
#include "aic_ecstar.h"

/* The genuine-product star: Phi = identity, so out = XY (the product already
 * formed by aic_ecstar_star). ctx/prec unused (id map, exact copy). */
static void matalg_star(acb_mat_t out, const acb_mat_t XY, void *ctx, slong prec)
{
    (void) ctx;
    (void) prec;
    acb_mat_set(out, XY);
}

void aic_cstar_matrix_algebra(aic_ecstar *out, slong d, slong prec)
{
    (void) prec;
    assert(out != NULL);
    assert(d >= 1 && "aic_cstar_matrix_algebra: need d >= 1");

    out->n = d;
    out->dim_A = d * d;
    out->phi = NULL;
    out->star_phi = matalg_star;
    out->star_ctx = NULL;

    /* The d^2 matrix units E_{lm}, row-major order index l*d + m. Each is a single
     * 1 at (l,m): Frobenius-orthonormal, ||E_lm||_F = 1. */
    out->B = (acb_mat_t *) flint_malloc((size_t) (d * d) * sizeof(acb_mat_t));
    assert(out->B != NULL);
    for (slong l = 0; l < d; l++)
        for (slong m = 0; m < d; m++) {
            acb_mat_t *E = &out->B[l * d + m];
            acb_mat_init(*E, d, d);
            acb_mat_zero(*E);
            acb_set_si(acb_mat_entry(*E, l, m), 1);
        }
}
