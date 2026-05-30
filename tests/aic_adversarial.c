/* aic_adversarial.c — adversarial instance generators, part 1 of 2 (bead aic-dbo.2,
 * Increment 1). Structural / spectral families: gen1 jordan_t13, gen3
 * degenerate_proj, gen4 near_degen_herm, gen5 graded_diag. The norm/boundary
 * families (gen2, gen6, gen7) live in aic_adversarial_nla.c. See aic_adversarial.h
 * for the corpus rationale and the per-generator knob -> property contracts.
 *
 * Determinism: every entry is set from an exact closed form (integers, exact dyadic
 * Hadamard entries, or an arb power r^k computed at the caller's prec). No RNG here;
 * the only fixed-seed pattern is in aic_adversarial_nla.c (gen2).
 *
 * Measured anchors (prec=256 unless noted), asserted by tests/test_adversarial.c:
 *   gen1 t=1e-1: rho=0.4642 ; t=1e-9: rho=1.000e-3  (rho = t^{1/3})
 *   gen3 dim=4,rank=2: ||P^2-P||=0 exactly (dyadic 1/sqrt(4)), tr(P)=2,
 *        ||(2P-I)^2-I||=0 exactly; dim=8,rank=3: defect ~1e-75 (1/sqrt(8) irrational
 *        => prec-floor ball, still machine-floor; tol 1e-60)
 *   gen4 n=4 gap=1e-1: min-gap=0.1 ; gap=1e-9: min-gap=1e-9
 *   gen5 n=6 range=10: kappa=10 ; range=1e8: kappa=1e8 (Gram-free certify, see test)
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"

/* gen1 — Jordan t^{1/3} (nla 1a). tex:540-544:
 *   X(t) = [[0,1,0],[0,0,1],[t,0,0]],  char poly lambda^3 - t,
 *   spec = {t^{1/3}, omega t^{1/3}, omega^2 t^{1/3}},  rho(X) = t^{1/3}.
 * A size-t perturbation in entry (2,0) bursts the triple-zero eigenvalue out to
 * radius t^{1/3} — the non-Lipschitz 1/3-power sensitivity the paper warns of. */
void aic_adv_jordan_t13(acb_mat_t out, double t, slong prec)
{
    (void) prec;
    assert(acb_mat_nrows(out) == 3 && acb_mat_ncols(out) == 3 &&
           "aic_adv_jordan_t13: out must be 3x3");
    acb_mat_zero(out);
    acb_set_si(acb_mat_entry(out, 0, 1), 1);
    acb_set_si(acb_mat_entry(out, 1, 2), 1);
    acb_set_d(acb_mat_entry(out, 2, 0), t); /* the evil knob */
}

/* Sylvester-Hadamard column q_k (k-th column of the normalised H_dim), written as
 * the real entries +-1/sqrt(dim) into column `col` of M. dim must be a power of 2.
 * H[i][k] = (-1)^popcount(i & k); normalised by 1/sqrt(dim). 1/sqrt(dim) is exact
 * dyadic only when dim is a perfect square power of 2 (4, 16, 64, ...); for the
 * other powers it is irrational, so we set it as a certified arb sqrt — the entries
 * are then balls but P = Q Q^dag is still an EXACT projector to machine floor
 * because Q's columns are certified orthonormal. */
static void hadamard_col(acb_mat_t M, slong col, slong dim, slong k, slong prec)
{
    arb_t s;
    arb_init(s);
    arb_set_si(s, dim);
    arb_rsqrt(s, s, prec); /* 1/sqrt(dim) */
    for (slong i = 0; i < dim; i++) {
        unsigned long bits = (unsigned long) (i & k);
        int sign = (__builtin_popcountl(bits) & 1) ? -1 : 1;
        if (sign < 0)
            arb_neg(acb_realref(acb_mat_entry(M, i, col)), s);
        else
            arb_set(acb_realref(acb_mat_entry(M, i, col)), s);
        arb_zero(acb_imagref(acb_mat_entry(M, i, col)));
    }
    arb_clear(s);
}

/* gen3 — Exact-degeneracy projector (nla 4b, lethal #1). P = Q_r Q_r^dag where Q_r
 * holds the first `rank` Hadamard columns (orthonormal). Spectrum {1 (mult rank),
 * 0 (mult dim-rank)} — repeated, so FLINT's simple-spectrum eig ABORTS; the eig-free
 * sgn path is the only route. Non-diagonal (Hadamard basis) so it is a genuine
 * stress, not a trivial diag(1,..,0). dim must be a power of 2 (asserted). */
void aic_adv_degenerate_proj(acb_mat_t out, slong dim, slong rank, slong prec)
{
    assert(acb_mat_nrows(out) == dim && acb_mat_ncols(out) == dim &&
           "aic_adv_degenerate_proj: out must be dim x dim");
    assert(rank > 0 && rank < dim && "aic_adv_degenerate_proj: need 0 < rank < dim");
    assert((dim & (dim - 1)) == 0 && "aic_adv_degenerate_proj: dim must be 2^m");

    acb_mat_t Q, Qd;
    acb_mat_init(Q, dim, rank);
    acb_mat_init(Qd, rank, dim);
    for (slong c = 0; c < rank; c++)
        hadamard_col(Q, c, dim, c, prec); /* column c uses Hadamard index c */
    acb_mat_conjugate_transpose(Qd, Q);
    acb_mat_mul(out, Q, Qd, prec); /* P = Q Q^dag, rank-`rank` projector */
    acb_mat_clear(Qd);
    acb_mat_clear(Q);
}

/* gen4 — Near-degenerate Hermitian (nla 4a). H = diag(1-gap/2, 1+gap/2, 3, 4, ..., n)
 * (real diagonal => Hermitian). The two eigenvalues 1-gap/2 and 1+gap/2 are exactly
 * `gap` apart and are the minimally-separated pair (the rest are integers 3..n, gaps
 * >= 1). As gap -> 0 the simple-spectrum certifier acb_mat_eig_simple can no longer
 * isolate them and aic_mat_eig_hermitian aborts (Rule 4). */
void aic_adv_near_degen_herm(acb_mat_t out, slong n, double gap, slong prec)
{
    (void) prec;
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_adv_near_degen_herm: out must be n x n");
    assert(n >= 2 && "aic_adv_near_degen_herm: need n >= 2");
    acb_mat_zero(out);
    /* near-degenerate pair around 1 */
    arb_set_d(acb_realref(acb_mat_entry(out, 0, 0)), 1.0 - gap / 2.0);
    arb_set_d(acb_realref(acb_mat_entry(out, 1, 1)), 1.0 + gap / 2.0);
    /* well-separated remainder 3, 4, ..., n (skip 2 so no entry is near 1) */
    for (slong i = 2; i < n; i++)
        acb_set_si(acb_mat_entry(out, i, i), i + 1);
}

/* gen5 — Graded diagonal (nla 5c). D = diag(r^0, r^1, ..., r^{n-1}) with
 * r = range^{1/(n-1)}, so sigma_max/sigma_min = r^{n-1} = range exactly. r^k is
 * computed as a certified arb power (arb_pow) so the dynamic range is exact to the
 * working precision, not contaminated by double rounding of r^k for large k. */
void aic_adv_graded_diag(acb_mat_t out, slong n, double range, slong prec)
{
    assert(acb_mat_nrows(out) == n && acb_mat_ncols(out) == n &&
           "aic_adv_graded_diag: out must be n x n");
    assert(n >= 2 && "aic_adv_graded_diag: need n >= 2");
    assert(range >= 1.0 && "aic_adv_graded_diag: range must be >= 1");
    acb_mat_zero(out);
    arb_t r, rk, exph;
    arb_init(r);
    arb_init(rk);
    arb_init(exph);
    arb_set_d(r, range);
    arb_set_si(exph, n - 1);
    arb_inv(exph, exph, prec);
    arb_pow(r, r, exph, prec); /* r = range^{1/(n-1)} */
    arb_one(rk);
    for (slong i = 0; i < n; i++) {
        arb_set(acb_realref(acb_mat_entry(out, i, i)), rk); /* D[i][i] = r^i */
        arb_mul(rk, rk, r, prec);
    }
    arb_clear(exph);
    arb_clear(rk);
    arb_clear(r);
}
