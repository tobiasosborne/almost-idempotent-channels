/* aic_funcalc_proj.c — theta(X) and the prop_P projection-extractor (beads
 * aic-4bg, "T-funcalc").
 *
 * approximate_algebras.tex:527-528 — theta and its use:
 *     Ptilde = theta(2P - I),   where   theta(X) = (1/2)(I + sgn(X)).
 * approximate_algebras.tex:524-533 — Proposition prop_P:
 *     Let P satisfy ||P^2 - P|| <= delta < 1/4, and Ptilde = theta(2P - I).
 *     Then Ptilde^2 = Ptilde, Ptilde commutes with P, and
 *     ||Ptilde - P|| <= ||2P - I|| O(delta).
 *   proof (tex:532): set X = 2P - I; then ||X^2 - I|| = 4||P^2 - P|| <= 4 delta < 1,
 *   so sgn(X) is in its validity domain and theta(X) maps it to an exact
 *   idempotent.
 *
 * Constructive route (CLAUDE.md Law 3). The paper's proof IS the algorithm:
 *   1. X = 2P - I.
 *   2. certify rho(P^2 - P) < 1/4  (equivalently rho(X^2 - I) = 4 rho(P^2-P) < 1).
 *   3. S = sgn(X)  (eig-free; aic_sgn auto-dispatches: Newton-Schulz in the
 *      op-norm basin, else the globally-convergent Newton, bead aic-8hz — handles
 *      the degenerate spectrum of a near-projector that the eig path cannot,
 *      file docstring of aic_funcalc_sgn.c, bead aic-w4o.1).
 *   4. Ptilde = (I + S)/2.
 *
 * THE SPECTRAL RELAXATION (bead aic-8hz). prop_P's hypothesis is
 * ||P^2-P|| <= delta < 1/4 (tex:525). For a Banach-algebra ELEMENT the norm is
 * the binding quantity, but for the finite-dim SUPEROPERATOR S_Phi the regularity
 * that actually governs the functional calculus is the SPECTRAL radius
 * rho(P^2-P) < 1/4 (sgn/theta need spec((2P-I)^2) off (-inf,0], tex:546-550, i.e.
 * 2P-I with no purely-imaginary eigenvalue, i.e. rho(I-(2P-I)^2) = 4 rho(P^2-P)<1).
 * The op-NORM guard ||P^2-P||_op < 1/4 is STRICTER than this and rejects valid
 * non-normal inputs (where ||.||_op >> rho), so we certify the SPECTRAL condition
 * via the eig-free Gelfand bound on M = I - (2P-I)^2 = 4(P^2-P): accept the first
 * k with ||M^k||_F^{1/k} certainly < 1 (rho(M)<1 <=> rho(P^2-P)<1/4). This is a
 * RELAXATION — it accepts strictly more P (the full non-normal eta<1/4 regime) —
 * and fails loud at the true rho(P^2-P) >= 1/4 boundary (Rule 4). The sign itself
 * is then handled by aic_sgn's auto-dispatch (global Newton out of the op-basin).
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_mat.h"
#include "aic/aic_funcalc.h"
#include "aic_funcalc_internal.h"

void aic_theta(acb_mat_t out, const acb_mat_t X, slong prec)
{
    slong n = acb_mat_nrows(X);
    assert(acb_mat_ncols(X) == n && "aic_theta: X must be square");
    assert(out != X && "aic_theta: out must not alias X");

    /* S = sgn(X) (asserts ||X^2 - I|| < 1 internally), then (I + S)/2. */
    acb_mat_t S, I;
    acb_mat_init(S, n, n);
    acb_mat_init(I, n, n);
    aic_sgn(S, X, prec);
    acb_mat_one(I);
    acb_mat_add(out, I, S, prec);                /* I + sgn(X) */
    acb_mat_scalar_mul_2exp_si(out, out, -1);    /* * 1/2 */

    acb_mat_clear(I);
    acb_mat_clear(S);
}

void aic_prop_P(acb_mat_t out, const acb_mat_t P, slong prec)
{
    slong n = acb_mat_nrows(P);
    assert(acb_mat_ncols(P) == n && "aic_prop_P: P must be square");
    assert(out != P && "aic_prop_P: out must not alias P");

    /* Certify the SPECTRAL prop_P hypothesis rho(P^2 - P) < 1/4 (tex:525),
     * EIG-FREE via the Gelfand bound on M = I - (2P-I)^2 = 4(P^2 - P): a certified
     * rho(M) < 1 is exactly rho(P^2 - P) < 1/4. This RELAXES the previous op-norm
     * guard ||P^2-P||_op < 1/4 to the spectral one (bead aic-8hz), reaching the
     * full non-normal regime; fail loud at the true boundary (Rule 4). */
    acb_mat_t P2, X, M;
    acb_mat_init(P2, n, n);
    acb_mat_init(X, n, n);
    acb_mat_init(M, n, n);

    /* X = 2P - I; M = I - X^2 = 4(P^2 - P). */
    acb_mat_scalar_mul_2exp_si(X, P, 1);         /* 2P */
    acb_mat_one(P2);                             /* reuse P2 as I */
    acb_mat_sub(X, X, P2, prec);                 /* X = 2P - I */
    acb_mat_sqr(M, X, prec);                     /* X^2 */
    acb_mat_sub(M, P2, M, prec);                 /* M = I - X^2 = 4(P^2 - P) */

    arb_t rho;
    arb_init(rho);
    slong kused = 0;
    if (!aic_funcalc_int_gelfand_rho(rho, &kused, M, 32, prec)) {
        fprintf(stderr, "aic_prop_P: cannot certify rho(P^2 - P) < 1/4 "
                "(prop_P hypothesis; tex:525). Certified Gelfand bound on "
                "rho(I-(2P-I)^2)=4 rho(P^2-P) did not fall < 1 by k=32 — P is "
                "at/over the eta=1/4 spectral boundary, or raise prec.\n");
        abort();
    }

    /* Ptilde = theta(2P - I); aic_sgn auto-dispatches (global Newton out of the
     * op-norm basin, in-basin Newton-Schulz). */
    aic_theta(out, X, prec);

    arb_clear(rho);
    acb_mat_clear(M);
    acb_mat_clear(X);
    acb_mat_clear(P2);
}
