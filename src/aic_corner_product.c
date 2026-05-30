/* aic_corner_product.c — Increment 2a of the corner module: the COMPRESSED
 * PRODUCT (compr_prod, .tex:1077-1082) and its unit Ptilde (.tex:1082) of section
 * 7. Builds on Increment 1 (aic_corner_Co, aic_corner_apply); see
 * include/aic_corner.h for the API contracts. The BLOCK BIJECTION alpha
 * (lem_alpha) lives in the sibling src/aic_corner_alpha.c (split for the ~200 LOC
 * limit, Rule 10).
 *
 * approximate_algebras.tex:1078-1082 — the compressed product (verbatim):
 *   "(X,Y) |-> X . Y : S_{P,Q} x S_{Q,R} -> S_{P,R},   X . Y = Co_{P,R}(XY).
 *    It is close to the ambient product in A, namely, ||X . Y - XY|| <=
 *    O(delta+eps) ||X|| ||Y||. ... the compressed product turns the Banach space
 *    S_P ... into an O(delta+eps)-C* algebra with unit Ptilde = Co_P(P)."
 *
 * THE STAR, NOT THE PLAIN PRODUCT (CLAUDE.md callout, LOAD-BEARING). "XY" in the
 * paper's compr_prod is the product in A, i.e. the Choi-Effros STAR X*Y =
 * Phi_tilde(XY), aic_ecstar_star (NOT acb_mat_mul). So X . Y = Co_{P,R}(X * Y):
 * compute Z = X * Y (n x n, via the star) then apply the d x d Co_{P,R} to Z
 * (aic_corner_apply). Dropping the star (plain Z = X.Y) is RED on the oblique
 * compress_idemp(4,2) fixture (tests/test_corner.c T8(a), non-vacuity gap 0.0556).
 *
 * WHY Co_{P,R} IS A NO-OP ON A VALID X*Y (and where its teeth are). For X in
 * S_{P,Q}, Y in S_{Q,R}, the star X*Y = Phi_tilde(XY) lands (approximately) in
 * S_{P,R} already, so Co_{P,R}(X*Y) = X*Y up to O(delta+eps); at eta=0 it is EXACT
 * (measured comp_gap = 0 on the oblique fixture). This is the paper's "close to
 * the ambient product" (.tex:1081). So the Co_{P,R} step is invisible at eta=0 and
 * VISIBLE only at eta>0 (where X*Y leaves S_{P,R} by O(eta)); test_corner T8(b)'s
 * eta>0 arm is the compression-step teeth (skipping the Co apply makes it RED,
 * "X.Y == X*Y at eta>0"). The cdot-specific content is thus: (i) it uses the STAR
 * (T8(a)), (ii) it then COMPRESSES (T8(b) eta>0); the Co_{P,R} map itself is the
 * Increment-1 aic_corner_apply, already tested (T5/T7).
 */
#include <assert.h>

#include <flint/acb_mat.h>

#include "aic_corner.h"
#include "aic_ecstar.h"

void aic_corner_cdot(acb_mat_t out, const aic_ecstar *A, const acb_mat_t CoPR,
                     const acb_mat_t X, const acb_mat_t Y, slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(out != X && out != Y && "cdot: out must not alias X or Y");
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n);
    assert(acb_mat_nrows(CoPR) == d && acb_mat_ncols(CoPR) == d &&
           "cdot: CoPR must be dim_A x dim_A");

    /* Z = X * Y (the Choi-Effros STAR, .tex:1078 product-in-A; NOT acb_mat_mul). */
    acb_mat_t Z;
    acb_mat_init(Z, n, n);
    aic_ecstar_star(Z, A, X, Y, prec);
    aic_corner_apply(out, A, CoPR, Z, prec);   /* X . Y = Co_{P,R}(X * Y) */
    acb_mat_clear(Z);
}

void aic_corner_Ptilde(acb_mat_t out, const aic_ecstar *A, const acb_mat_t P,
                       slong prec)
{
    slong n = A->n, d = A->dim_A;
    assert(out != P && "Ptilde: out must not alias P");
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n);
    acb_mat_t CoPP;
    acb_mat_init(CoPP, d, d);
    aic_corner_Co(CoPP, A, P, P, prec);        /* Co_{P,P} */
    aic_corner_apply(out, A, CoPP, P, prec);   /* Ptilde = Co_P(P) (.tex:1082) */
    acb_mat_clear(CoPP);
}
