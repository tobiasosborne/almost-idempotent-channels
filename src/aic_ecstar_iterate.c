/* aic_ecstar_iterate.c — the HOPM inner restart loop + block updates with the
 * SUBSPACE-POLAR ACCEPT-GUARD (bead aic-knm, Cycle 2). The fast double-path
 * search engine; the certified driver is aic_ecstar_search.c. Split for Rule 10.
 *
 * Method 1 of the web leg (scale-invariant alternating maximization). For each
 * block (X, then Y, then Z) the Frobenius gradient C = sum_k c_k B_k of the
 * objective restricted to A is formed (c_k = <u, h(B_k, ...) v>, h the relevant
 * multilinear defect map); the candidate is X' = Pi_A(polar(C)) re-normalized to
 * the op-norm unit sphere, ACCEPTED only if the scale-invariant objective
 * strictly improves (else the old block is kept). This keeps every accepted
 * iterate exactly in A (valid witness) and makes the approximate-polar step
 * monotone. See include/aic_ecstar.h for the full rationale.
 *
 * Deterministic init: a fixed-counter LCG seeded from the restart index (no
 * wall-clock RNG, project rule). Half the restarts warm-start from the
 * Frobenius-relaxation top-singular direction; half from a random op in A.
 */
#include <complex.h>
#include <math.h>

#include "aic_ecstar_iterate.h"

/* xoshiro-ish LCG on a 64-bit state; returns a double in [-1,1). Deterministic. */
static double lcg_pm1(uint64_t *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    uint64_t x = *s >> 11; /* top 53 bits */
    return ((double) x / (double) (1ULL << 53)) * 2.0 - 1.0;
}

/* Random op in A: x_k = N(0,1)-ish via LCG, X = sum x_k B_k, normalized by op-norm.
 * Returns 0 on a numerically-zero draw (caller retries / falls back). */
static int rand_in_A(cx *X, const aic_ehk *h, uint64_t *seed, cx *scratch_cf)
{
    slong n = h->n, d = h->d;
    for (slong p = 0; p < n * n; p++) X[p] = 0.0;
    for (slong k = 0; k < d; k++) {
        cx c = lcg_pm1(seed) + lcg_pm1(seed) * I;
        scratch_cf[k] = c;
        for (slong p = 0; p < n * n; p++) X[p] += c * h->B[k][p];
    }
    double nrm = aic_ehk_opnorm(X, n);
    if (nrm < 1e-12) return 0;
    for (slong p = 0; p < n * n; p++) X[p] /= nrm;
    return 1;
}

/* Random unit vector in C^n via LCG. */
static void rand_unit_vec(cx *v, slong n, uint64_t *seed)
{
    double s = 0.0;
    for (slong k = 0; k < n; k++) {
        v[k] = lcg_pm1(seed) + lcg_pm1(seed) * I;
        s += creal(v[k] * conj(v[k]));
    }
    s = sqrt(s);
    if (s < 1e-300) { v[0] = 1.0; return; }
    for (slong k = 0; k < n; k++) v[k] /= s;
}

/* Build the A-gradient C = sum_k <u, term(B_k) v> B_k and write the candidate
 * block X' = Pi_A(polar(C))/||.||_op into cand. `term` fills w->g (n*n) with the
 * defect-map value when the active block is set to B_k. Returns 1 if cand is a
 * valid op-norm-1 matrix in A, else 0 (degenerate; keep the old block). The
 * gradient sign is +1 to MAXIMIZE, -1 to MINIMIZE (cstar). */
int aic_ehk_block_candidate(cx *cand, aic_ehk_work *w, const aic_ehk *h,
                            const cx *u, const cx *v, double sign,
                            void (*term)(cx *out, aic_ehk_work *w,
                                         const aic_ehk *h, slong k))
{
    slong n = h->n, d = h->d;
    /* C = sum_k c_k B_k, c_k = sign * <u, term(B_k) v> */
    for (slong p = 0; p < n * n; p++) w->C[p] = 0.0;
    for (slong k = 0; k < d; k++) {
        term(w->g, w, h, k);             /* w->g = defect map with block <- B_k */
        /* <u, g v> = u^dag g v */
        for (slong i = 0; i < n; i++) {
            cx s = 0.0;
            for (slong j = 0; j < n; j++) s += w->g[i * n + j] * v[j];
            w->vec[i] = s;
        }
        cx ip = 0.0;
        for (slong i = 0; i < n; i++) ip += conj(u[i]) * w->vec[i];
        cx ck = sign * ip;
        for (slong p = 0; p < n * n; p++) w->C[p] += ck * h->B[k][p];
    }
    /* polar(C) then project back onto A then normalize */
    aic_ehk_polar(w->pol, w->C, n, w->U, w->Vt, w->sv);
    aic_ehk_project_A(cand, w->cf, h, w->pol);
    double nrm = aic_ehk_opnorm(cand, n);
    if (nrm < 1e-12) return 0;
    for (slong p = 0; p < n * n; p++) cand[p] /= nrm;
    return 1;
}

/* Initialize block X for restart `ri`: warm-start (ri even) from the top-singular
 * direction of a random A-op's polar (a Frobenius-relaxation surrogate), or a
 * plain random op in A (ri odd). Falls back to the first basis op if degenerate. */
void aic_ehk_init_block(cx *X, const aic_ehk *h, uint64_t *seed, int warm,
                        aic_ehk_work *w)
{
    slong n = h->n;
    if (!rand_in_A(X, h, seed, w->cf)) {
        for (slong p = 0; p < n * n; p++) X[p] = h->B[0][p];
        double nb = aic_ehk_opnorm(X, n);
        if (nb > 1e-12) for (slong p = 0; p < n * n; p++) X[p] /= nb;
        return;
    }
    if (warm) {
        /* one polar+project pass tightens the random draw toward an extreme point
         * of the op-norm ball in A (cheap Frobenius-relaxation warm start). */
        aic_ehk_polar(w->pol, X, n, w->U, w->Vt, w->sv);
        aic_ehk_project_A(w->C, w->cf, h, w->pol);
        double nrm = aic_ehk_opnorm(w->C, n);
        if (nrm > 1e-12)
            for (slong p = 0; p < n * n; p++) X[p] = w->C[p] / nrm;
    }
}

/* Random unit vectors u,v for restart init. */
void aic_ehk_init_uv(cx *u, cx *v, slong n, uint64_t *seed)
{
    rand_unit_vec(u, n, seed);
    rand_unit_vec(v, n, seed);
}

/* Deterministic per-restart seed from the restart index (no wall clock). */
uint64_t aic_ehk_seed(int restart_index)
{
    return 0x9E3779B97F4A7C15ULL ^
           ((uint64_t) (restart_index + 1) * 0xD1B54A32D192ED03ULL);
}
