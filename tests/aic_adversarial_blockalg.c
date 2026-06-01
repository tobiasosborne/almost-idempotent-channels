/* aic_adversarial_blockalg.c — adversarial CHANNEL generator targeting the
 * DIMENSION-BLOWUP regime eps ~ c/n (bead aic-cxo, the HIGH "3D" item). A UCP
 * self-map on B(C^N), N = k*d, whose associated eps-C* algebra is
 *     A = (+)_{j=1}^{k} M_d,    dim_A = k * d^2,
 * eta-idempotent with eta TUNABLE by a knob t (t=0 => EXACTLY idempotent). This
 * is the "dimension-dependent constant blowup" family catalogued at
 * docs/adversarial/domain.md:416-449 (family 3D). Sibling of the other channel
 * families (aic_adversarial_domain.c cb-op-gap/depol/unital-defect,
 * aic_adversarial_carrier.c near-degenerate carrier); kept in its own file
 * (Rule 10). See aic_adversarial.h for the corpus rationale (CLAUDE.md
 * adversarial-first directive, docs/adversarial/README.md).
 *
 * WHY THIS FAMILY (domain.md:426-431, tex:484). th_main (.tex:460) promises an
 * O(eps)-isomorphism with a UNIVERSAL, DIMENSION-INDEPENDENT constant c_0. The
 * paper explicitly WARNS (tex:484) that the NAIVE cohomological route — fixing
 * the multiplication via a Haar-measure diagonal D = int dU (U^dag (x) U) —
 * fails when eps ~ c/n, because a naive diagonal construction has error bounded
 * BELOW by something proportional to n = dim A. The correct route uses the
 * EXPLICIT generalized-Pauli (Heisenberg-Weyl) diagonal of the GENUINE C*
 * algebra B (tex:1249-1254), which has ||D|| = 1 regardless of n. This generator
 * is the witness whose dim_A = k*d^2 grows with k: feeding the family at a SWEEP
 * over k exposes whether cstar_build / errreduce / dhom keep c_0 flat (correct,
 * the explicit-Pauli-diagonal route) or let it grow ~n (the naive route the
 * paper warns against). The eta=0 instance is the cleanest: the recovered block
 * structure (+)_j M_d is EXACTLY pinnable against aic_idemp_decompose, so the
 * self-test asserts the block COUNT tracks k and dim_B == k*d^2.
 *
 * THE CONSTRUCTION (mix a k-block conditional expectation with its DFT conjugate).
 *   1. BLOCK PROJECTORS on C^N (N = k*d):
 *        P_j = sum_{a=0}^{d-1} |jd+a><jd+a|,   j = 0 .. k-1.
 *      The P_j partition the identity: sum_j P_j = I_N, P_j P_l = delta_{jl} P_j
 *      (disjoint contiguous index ranges [jd, jd+d)). Each P_j is a rank-d
 *      diagonal Hermitian projector.
 *   2. BASE = the k-block conditional expectation
 *        Phi0(X) = sum_j P_j X P_j.
 *      In the OBSERVABLE convention Phi(X) = sum_a K_a^dag X K_a (aic_ucp.h) this
 *      has Kraus K_j = P_j (P_j Hermitian, P_j^dag = P_j). Img Phi0 = the
 *      block-diagonal matrices = (+)_j M_d, dim = k*d^2. Phi0 is EXACTLY
 *      idempotent (Phi0^2 = Phi0, since P_j X P_j is already P_l-block-diagonal,
 *      so re-applying Phi0 is the identity on it) and UNITAL
 *      (sum_j K_j^dag K_j = sum_j P_j = I_N). This is the k-block generalization
 *      of test_idemp.h:make_block_cond_exp (which is the k=2 case).
 *   3. A FIXED, DETERMINISTIC, NON-BLOCK-DIAGONAL unitary V on C^N — the DFT:
 *        V[a,b] = exp(2 pi i a b / N) / sqrt(N),   a,b = 0 .. N-1.
 *      V is unitary (the discrete Fourier transform; V^dag V = I_N exactly in the
 *      ideal, certified to the arb ball here). It is NOT block-diagonal w.r.t. the
 *      k-block partition for k >= 2 (the DFT couples ALL coordinates), so the
 *      conjugate's image V^dag (Img Phi0) V genuinely DIFFERS from (+)_j M_d — the
 *      property that makes the mix have eta > 0 for t > 0 (a block-diagonal V
 *      would commute with the partition and give eta = 0 at all t). Built in arb
 *      (certified), the SAME for a given N (reproducible, no RNG).
 *   4. CONJUGATE Kraus K'_j = V^dag P_j V (the channel Ad_{V^dag} . Phi0 . Ad_V,
 *      mirroring test_idemp.h:make_conjugated). MIX: the output channel has the
 *      Kraus UNION
 *        { sqrt(1-t) * P_j : j } U { sqrt(t) * V^dag P_j V : j },
 *      Kraus rank r = 2k. So
 *        out(X) = (1-t) sum_j P_j X P_j  +  t sum_j (V^dag P_j V)^dag X (V^dag P_j V)
 *               = (1-t) Phi0(X)  +  t (Ad_{V^dag} Phi0 Ad_V)(X).
 *      UNITALITY (for ALL t in [0,1]):
 *        sum_a K_a^dag K_a = (1-t) sum_j P_j + t sum_j V^dag P_j V
 *                          = (1-t) I_N + t V^dag (sum_j P_j) V
 *                          = (1-t) I_N + t V^dag I_N V = (1-t) I_N + t I_N = I_N.
 *      So out is UCP and unital for every t. (This is the make_mixconj pattern of
 *      test_idemp.h generalized from a single compression block to k equal blocks;
 *      built SELF-CONTAINED here — no #include of the test-only statics, so the
 *      corpus is reusable without test-only headers.)
 *
 * t=0 REDUCTION (the eta=0 oracle, CLAUDE.md cross-check ladder #3). At t=0 the
 * second Kraus group has weight 0, so out = Phi0 = the k-block conditional
 * expectation, EXACTLY idempotent (Phi0^2 = Phi0), defect 0. Its associated
 * algebra is (+)_j M_d with dim_A = k*d^2 and exactly k blocks each of dim d —
 * the structural property the self-test pins (T1 oracle + T3 block structure).
 *
 * KNOB t in [0,1) tunes eta: at t > 0 the DFT conjugate's image is misaligned
 * with the block partition, so Phi0 and Ad_{V^dag} Phi0 Ad_V do NOT share an
 * idempotent, and ||out^2 - out||_cb > 0 grows with t (the corpus-can't-go-
 * toothless monotone, T2). t < 1 keeps the base block weight positive.
 *
 * USE k >= 2 FOR THE DIMENSION-BLOWUP PROPERTY. k = 1 is a legal unital channel
 * but DEGENERATE for this family's purpose: with a single block P_0 = I_N the
 * conjugate is V^dag P_0 V = V^dag I_N V = I_N, identical to the base, so eta = 0
 * for ALL t (no t-dependence, no blowup). The named adversarial property — dim_A
 * = k*d^2 growing while c_0 stays flat (tex:484) — needs k >= 2 (and a sweep over
 * k). k = 1 is allowed (the assert is k >= 1) only so the degenerate case is a
 * well-defined single M_d, not a crash.
 *
 * Asserts (Rule 4, fail loud): k >= 1, d >= 1, k*d >= 2, 0 <= t < 1. The output
 * `out` is aic_ucp_kraus_init'd HERE (dim_K = dim_H = N, r = 2k; caller clears
 * with aic_ucp_kraus_clear), matching the aic_adv_chan_* convention.
 *
 * Determinism: closed form (P_j exact integers, V the arb DFT, sqrt of doubles).
 * No RNG; the same (k,d,t,prec) yields a bit-identical channel.
 *
 * Measured anchors (prec=256), asserted by tests/test_adversarial.c:
 *   (k=2,d=2) eta bracket lower bounds (eig-free ||Phi^2-Phi||_cb, lo = ||J||_F/N):
 *     t=0 -> 0 (exact idempotent, bracket [0,0]); t=0.02 -> 0.012002;
 *     t=0.08 -> 0.045071 (strictly increasing, mild < lethal).
 *   recovered block structure at t=0: (k=2,d=2) num_blocks=2, d=[2,2], dim_B=8,
 *     iso_def 0; (k=3,d=2) num_blocks=3, d=[2,2,2], dim_B=12, iso_def 0. The
 *     block COUNT tracks k (dim_A = k*d^2 grows with k — the eps~c/n regime).
 */
#include <assert.h>
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_adversarial.h"
#include "aic/aic_ucp.h"

/* Build the DFT unitary V on C^N: V[a,b] = exp(2 pi i a b / N) / sqrt(N). Unitary
 * (certified to the arb ball at `prec`); NON-block-diagonal w.r.t. the k-block
 * partition for k >= 2. `V` must be init'd N x N. */
static void blockalg_dft(acb_mat_t V, slong N, slong prec)
{
    arb_t two_pi, invsqrtN, ang, c, s;
    arb_init(two_pi);
    arb_init(invsqrtN);
    arb_init(ang);
    arb_init(c);
    arb_init(s);
    arb_const_pi(two_pi, prec);
    arb_mul_2exp_si(two_pi, two_pi, 1);     /* 2 pi */
    arb_set_si(invsqrtN, N);
    arb_rsqrt(invsqrtN, invsqrtN, prec);    /* 1 / sqrt(N) */
    for (slong a = 0; a < N; a++)
        for (slong b = 0; b < N; b++) {
            /* angle = 2 pi (a*b mod N) / N — reduce a*b mod N first (exact). */
            slong m = (a * b) % N;
            arb_mul_si(ang, two_pi, m, prec);
            arb_div_si(ang, ang, N, prec);
            arb_sin_cos(s, c, ang, prec);   /* exp(i ang) = c + i s */
            arb_mul(c, c, invsqrtN, prec);
            arb_mul(s, s, invsqrtN, prec);
            acb_set_arb_arb(acb_mat_entry(V, a, b), c, s);
        }
    arb_clear(s);
    arb_clear(c);
    arb_clear(ang);
    arb_clear(invsqrtN);
    arb_clear(two_pi);
}

/* fam3D — dimension-blowup block algebra (domain.md:416-449, tex:484/1249). A
 * UCP self-map on B(C^N), N=k*d, whose associated eps-C* algebra is (+)_j M_d
 * (dim_A = k*d^2), eta-idempotent with eta tunable by t (t=0 => exact). Mix of the
 * k-block conditional expectation Phi0(X)=sum_j P_j X P_j (Kraus P_j) with its DFT
 * conjugate Ad_{V^dag} Phi0 Ad_V (Kraus V^dag P_j V): Kraus union
 *   { sqrt(1-t) P_j } U { sqrt(t) V^dag P_j V },  rank 2k, unital for all t.
 * KNOB t in [0,1). Asserts k>=1, d>=1, k*d>=2, 0<=t<1 (Rule 4, fail loud). `out`
 * is aic_ucp_kraus_init'd HERE (dim_K=dim_H=N, r=2k; caller clears it). See the
 * file docstring for the full Kraus derivation, unitality proof, and t=0 oracle. */
void aic_adv_chan_blockalg(aic_ucp_kraus *out, slong k, slong d, double t,
                           slong prec)
{
    assert(k >= 1 && "aic_adv_chan_blockalg: need k >= 1");
    assert(d >= 1 && "aic_adv_chan_blockalg: need d >= 1");
    assert(k * d >= 2 && "aic_adv_chan_blockalg: need N = k*d >= 2");
    assert(t >= 0.0 && t < 1.0 && "aic_adv_chan_blockalg: need 0 <= t < 1");

    slong N = k * d;
    /* dim_K = dim_H = N (self-map on B(C^N)); Kraus rank r = 2k (P_j and conj). */
    aic_ucp_kraus_init(out, N, N, 2 * k);

    /* V = DFT, V^dag. */
    acb_mat_t V, Vd, P, tmp;
    acb_mat_init(V, N, N);
    acb_mat_init(Vd, N, N);
    acb_mat_init(P, N, N);     /* scratch: P_j (rebuilt per block) */
    acb_mat_init(tmp, N, N);   /* scratch: V^dag P_j */
    blockalg_dft(V, N, prec);
    acb_mat_conjugate_transpose(Vd, V);

    arb_t w0, w1;
    arb_init(w0);
    arb_init(w1);
    arb_set_d(w0, sqrt(1.0 - t)); /* base block weight sqrt(1-t) */
    arb_set_d(w1, sqrt(t));       /* conj block weight sqrt(t)   */

    for (slong j = 0; j < k; j++) {
        /* P_j = sum_{a=0}^{d-1} |jd+a><jd+a| (diagonal projector onto block j). */
        acb_mat_zero(P);
        for (slong a = 0; a < d; a++)
            acb_set_si(acb_mat_entry(P, j * d + a, j * d + a), 1);

        /* K_j = sqrt(1-t) P_j. */
        acb_mat_scalar_mul_arb(out->K[j], P, w0, prec);

        /* K_{k+j} = sqrt(t) V^dag P_j V. */
        acb_mat_mul(tmp, Vd, P, prec);            /* V^dag P_j   */
        acb_mat_mul(out->K[k + j], tmp, V, prec); /* V^dag P_j V */
        acb_mat_scalar_mul_arb(out->K[k + j], out->K[k + j], w1, prec);
    }

    arb_clear(w1);
    arb_clear(w0);
    acb_mat_clear(tmp);
    acb_mat_clear(P);
    acb_mat_clear(Vd);
    acb_mat_clear(V);
}
