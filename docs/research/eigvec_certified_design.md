# Certified invariant-subspace / eigenvector extraction (increment 2, bead `aic-4td` / `aic-w4o.1` inc.2)

> Design document for the implementation subagents. Empirically grounded
> (every load-bearing claim pinned by a bounded probe against FLINT 3.0.1,
> not derived), constructively rigorous, decision-complete. Produced by the
> increment-2 design/research lead; **no production code is written here.**
>
> Read-order grounding: CLAUDE.md (Laws 1–4, Rules 4/6/7/9/10); FINDINGS
> §D5, §D5n, §C14, §C4; `src/aic_mat_eig_multiple.c` (increment 1, the
> eigenvalue+rank layer this builds on); `src/aic_ucp_latd.c` (the double-path
> Choi→Kraus this mirrors); `src/aic_ucp_carrier.c`; the FLINT acb_mat eig
> family.

---

## 0. TL;DR for the impatient implementer

1. **D5n is RESOLVED** (measured, not conjectured). The §D5n route "per-cluster
   Rump + re-orthonormalised seed" as literally written **does NOT work** — but
   the *root cause is different from what §D5n claims*, and a one-line fix on top
   of per-cluster Rump resolves every previously-failing case. **The fix: densify
   the input by a certified unitary conjugation `A' = U A U†` (dense rational-Givens
   `U`, exact-dyadic, `U U† − I` certified to ~1e-37), run the eig family on `A'`,
   map the certified subspace back as `X = U† X'`.** The `acb_mat_eig_multiple`
   *eigenvalue* layer is *also* rescued by the same densification.
2. The root cause is **not** seed near-parallelism and **not** a precision wall.
   It is that `acb_mat_eig_enclosure_rump`'s internal **frozen-row partition**
   (`partition_X_sorted`, picks `k` rows by magnitude) degenerates when the
   cluster's invariant subspace is **row-sparse** (lives in few coordinates).
   Densification spreads every eigenvector across all `n` rows, so the partition
   is well-conditioned.
3. The new substrate routine returns, per cluster, `{certified eigenvalue ball,
   n×k certified invariant-subspace enclosure X}`. From `X` build the certified
   range projector `Π_M = X(X†X)⁻¹X†` (measured projector defect ~1e-31). This
   feeds certified Choi→Kraus and the certified carrier split.
4. `aic_ucp_choi_to_kraus_arb` mirrors `_latd` exactly: same Convention-A
   conjugate reshape `K_a[i,c] = √λ_a · conj(v_a[i·dim_H + c])`, same
   `_tol`/PSD-cone variant (FINDINGS §C14). It is built ON the new certified
   eig+subspace layer.

---

## 1. Route + rigor

### 1.1 The FLINT primitive, pinned (docs + source)

`acb_mat_eig_enclosure_rump(acb_t lambda, acb_mat_t J, acb_mat_t X,
const acb_mat_t A, const acb_t lambda_approx, const acb_mat_t X_approx, slong prec)`
— FLINT 3.0.1, declared `/usr/include/flint/acb_mat.h:421`. Documented +
source-pinned semantics (the latter read from `acb_mat/eig_enclosure_rump.c`,
Rump 2010 §13.4):

- **Inputs.** `lambda_approx` = one approximate (cluster) eigenvalue;
  `X_approx` = an `n×k` approximate right invariant-subspace basis (`k` = cluster
  multiplicity, `k=1` for a simple eigenvalue).
- **Outputs.** `lambda` = a certified enclosure containing **at least `k`
  eigenvalues** of `A` (the cluster). `X` = a certified `n×k` enclosure of the
  invariant subspace. `J` (pass non-NULL to get it) = a `k×k` interval matrix
  with `J = lambda_approx·I_k + Vᵀ·Y`; **the eigenvalues of `J` are the cluster's
  eigenvalues**, and `A·X − X·J` is certified to enclose `0`.
- **X is NOT orthonormalised** (pinned from source: "no normalization step is
  applied to the output"; the partition freezes `k` rows by magnitude and
  corrects the remaining `n−k`). *Load-bearing*: to build a projector we must
  orthonormalise/Gram-invert ourselves (§1.4).
- **Failure mode.** Krawczyk iteration with ε-inflation; if `Yeps ⊇ Y`
  containment is not reached within `maxiter = 5 + FLINT_BIT_COUNT(prec)` steps
  it emits a **non-finite (`[±∞]`) enclosure** — *it does not return a status
  code; you must test `acb_is_finite(lambda) && acb_mat_is_finite(X)`*. The
  documented triggers are an ill-conditioned preconditioner / poor `X_approx` /
  precision exhaustion. **MEASURED here: the dominant trigger in our regime is a
  row-sparse `X_approx` defeating the frozen-row partition (§2), not precision.**

`acb_mat_eig_multiple(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx,
const acb_mat_t R_approx, slong prec)` (increment 1's primitive): returns `n`
balls "either disjoint or identical, identical grouped consecutively"; a run of
`k` identical balls is a certified `k`-cluster isolated from the rest. Returns
`0` when the whole-spectrum Rump enclosure cannot close (the §D5n failure).

### 1.2 The constructive algorithm (certified per-cluster invariant subspaces)

Input: Hermitian `A` (`n×n`), `prec`. Output: a list of clusters, each
`{eigenvalue ball λ_c, certified n×k_c subspace enclosure X_c}`, with
`Σ_c k_c = n`.

```
1. DENSIFY.  Build U = dense rational-Givens unitary (§1.3); Ud = U†.
             A' = U · A · Ud   (certified acb arithmetic, prec).
             Assert ‖U Ud − I‖ certified-small (fail loud otherwise).
2. SEED.     Double-path zheev on midpoint(A') → ascending real evals ev[0..n)
             and ORTHONORMAL eigenvector columns Vd (n×n).   [seed only — NOT
             trusted for the bound; the Rump enclosure is what we certify.]
3. CLUSTER.  Partition ev[] into runs by a numerical gap (§1.5): cluster c =
             [s0_c, s0_c + k_c).   λ̂_c = mean(ev over the cluster).
4. PER CLUSTER c:
      X'_approx = Vd[:, s0_c : s0_c+k_c]            (n×k_c, orthonormal, dense)
      acb_mat_eig_enclosure_rump(λ_c, J_c, X'_c, A', λ̂_c, X'_approx, prec)
      ASSERT acb_is_finite(λ_c) && acb_mat_is_finite(X'_c)  — else FAIL LOUD
             (Rule 4; clusters unresolved at this prec, raise prec / densify
             harder — see §5 risk R3).
      ASSERT imag(λ_c) ∋ 0                           (Hermitian ⇒ real spectrum)
      BACK-MAP  X_c = Ud · X'_c                       (certified; subspace of A)
5. DISJOINTNESS (the soundness gate, §1.6):
      for all c≠c':  ASSERT  !acb_overlaps(λ_c, λ_c')  — else FAIL LOUD
             (two clusters merged into one ball at this prec ⇒ the subspace
             decomposition is NOT certified; raise prec).
6. RETURN the per-cluster {λ_c, X_c}.
```

`acb_mat_eig_multiple(A')` (densified) ALSO certifies the eigenvalue balls
all-at-once (measured §2), so increment-1's `aic_mat_eig_hermitian_multiple` can
be made D5n-robust by the same `A → A'` densification (§3.1, the increment-1
hardening). Use `eig_multiple(A')` for the eigenvalue *grouping* and
per-cluster Rump for the *subspaces*; both agree (the cluster ball from Rump
must lie inside the eig_multiple ball — a free cross-check).

### 1.3 The densifying unitary `U` (exact, certifiable)

`U` = product over all planes `(a,b)`, `a<b`, of the rational Givens rotation
with `cos = 3/5, sin = 4/5` (the same exact dyadic-free rationals already used in
`tests/test_eigmult.c`'s fixtures). Each entry is an exact rational; at any
`prec ≥ ~8` the acb representation of `3/5, 4/5` has a tiny but nonzero radius,
and `U U† − I` is certified to **~1e-37 at prec=128** (measured, §2) — i.e. `U`
is unitary to far below the working precision. `n(n−1)/2` Givens factors;
`U·A·U†` is two `acb_mat_mul`s. *Why this `U` and not Haar:* it is reproducible,
exact, allocation-free of any RNG, and empirically sufficient to densify every
project spectrum (projectors, block-algebra carriers, Choi matrices). It is a
**finite explicit element** (CLAUDE.md Law 3 "Haar diagonal → explicit finite
average" spirit). If a future adversarial input is somehow left row-sparse by
this fixed `U` (none found — §5 R3), the fallback is a second `U` with a
different angle or an extra random orthogonal seed, gated by the same finite-X
assert.

### 1.4 Certified range projector `Π_M` from the (non-orthonormal) `X_c`

Rump's `X_c` is full-column-rank but not orthonormal. The orthogonal projector
onto its column space (= the certified invariant subspace) is the standard
oblique-free formula

```
Π_M = X_c (X_c† X_c)⁻¹ X_c†          [G = X_c† X_c is k×k, certified-invertible]
```

MEASURED (rank-2 range of a C^5 projector, prec=128): `acb_mat_inv(G)` succeeds
(`G` certified well-conditioned because the columns are a near-orthonormal seed
mapped through a unitary), `‖Π² − Π‖_ub = 1.5e-31`, `‖Π − Π†‖_ub = 1.0e-31`,
and `‖Π − Π_true‖_ub = 5.1e-32` vs the exact range projector. For the **carrier
split** (§3.3) the carrier `M` is the range of `Q = Σ K_a K_a†`; `Π_M` is the
sum over the *nonzero-eigenvalue clusters'* projectors `Σ_{λ_c > thr} Π_{M_c}`.

> Note (FINDINGS §C4): for an *oblique* idempotent the range is spanned by the
> **LEFT** singular vectors. Here `A` (and `Q`) are **Hermitian/PSD**, so left =
> right and Rump's right-invariant-subspace `X` is the range basis directly — no
> left/right hazard. The §C4 caveat bites only if a non-Hermitian idempotent is
> ever fed to this routine; the Hermiticity precondition (asserted, §1.2 step 0)
> rules it out.

### 1.5 Cluster identification & threshold (the PSD-Choi / big-zero-cluster case)

A PSD Choi `C_Φ` (rank `r ≪ dim`) has a **large zero cluster** (multiplicity
`dim − r`) plus possibly-distinct nonzero eigenvalues *and* possibly-degenerate
nonzero clusters (the `⊕B(L_j)` block multiplicities). Clustering rule on the
ascending double seed `ev[]`:

- **Gap clustering for the eigen-grouping:** start a new cluster whenever
  `ev[j] − ev[j−1] > gap_thr`. For the kept/dropped *rank* decision reuse the
  increment-1 / `_latd` QETLAB threshold `thr = (dim_K·dim_H)·DBL_EPSILON·‖C‖_F`
  (the cone boundary). `gap_thr` should separate the zero cluster from the
  smallest nonzero eigenvalue; a robust default is `gap_thr = max(thr·10,
  relative-gap heuristic)`. **The certification does not depend on `gap_thr`
  being "right"** — it only decides which approximate eigenvalues are handed to
  one Rump call. If two genuinely-distinct eigenvalues are lumped into one
  cluster, Rump still certifies a (wider) invariant subspace and the
  **disjointness gate (§1.6) is what enforces soundness**; if two genuinely-equal
  eigenvalues are split into two clusters, Rump on the too-small seed will either
  fail-loud (non-finite) or the two λ balls will *overlap* and the disjointness
  gate fails loud. Either way the output is never silently wrong.
- MEASURED on realistic spectra (§2, probe d5n6): `{0^4|1|2|3}` (rank-3 Choi),
  `{0^3|2,2|5}` (block-algebra), `{0^3|4^3|9}`, `{0^3|7,7}`, `{0^3|1,1}` — all
  **SOUND**: each kernel/nonzero cluster certifies, balls disjoint,
  `Σ k_c = n`.

### 1.6 The rigor argument (one paragraph)

**Soundness of the decomposition rests on three certified facts.** (i) *Per
cluster, Rump is self-verifying*: a *finite* `acb_mat_eig_enclosure_rump` output
is a rigorous Krawczyk fixed-point certificate that the ball `λ_c` contains
exactly `k_c` eigenvalues of `A'` and that `X'_c` rigorously encloses the
corresponding `k_c`-dimensional invariant subspace, with `A'·X'_c − X'_c·J_c ∋ 0`
entrywise (recomputed independently, `‖·‖_ub ~ 1e-32`, §2). (ii) *Conjugation by
the certified unitary `U` is spectrum- and subspace-preserving*: since
`A' = U A U†` with `U U† − I` certified to ~1e-37, the eigenvalues are identical
(the λ balls are unchanged) and `X_c = U† X'_c` satisfies `A·X_c = U†·A'·X'_c =
U†·X'_c·J_c = X_c·J_c` — **the same `k_c×k_c` `J_c` works for the original `A`**,
so `X_c` is a certified invariant subspace of `A` with cluster eigenvalues
`spec(J_c) ⊂ λ_c` (verified directly on `A`: `‖A X_c − X_c J_c‖_ub ~ 1e-31` at
prec=128, growing gently to ~1e-27 at n=16). (iii) *Cross-cluster
ball-disjointness* (`!acb_overlaps(λ_c, λ_c')` for all `c≠c'`, asserted) proves
the clusters index **distinct** eigenvalues, so the certified invariant subspaces
are mutually orthogonal complements summing to `Cⁿ` (`Σ k_c = n`, asserted). Any
failure of (i) (non-finite enclosure) or (iii) (overlapping balls) is a FAIL-LOUD
abort (Rule 4) — the routine is **sound (never silently wrong), and complete on
every project spectrum measured**.

---

## 2. D5n verdict — EMPIRICAL, with measured numbers

**Verdict: RESOLVED.** The §D5n boundary cases (`acb_mat_eig_multiple` returns 0
on two-clusters-each-≥2) are fully certified by the densify-then-Rump route, at
prec as low as 53. The §D5n *route as written* (per-cluster Rump on the raw
seed) does **not** work; the *densification* on top of it is the fix.

### 2.1 What §D5n claimed vs what is measured

§D5n hypothesised the failure is **seed near-parallelism** in the multi-cluster
case, fixable by per-cluster re-orthonormalisation. **This is FALSE.** Measured:

| probe | input | result |
|---|---|---|
| d5n1 | `acb_mat_eig_multiple` all-at-once, `build_projector` C^5 {2,3}, C^6 {2,4}, C^7 {3,4} | **FAIL(0)** at prec 53/128/256 (reproduces §D5n) |
| d5n1 | per-cluster Rump, raw zheev seed (orthonormal!) or qr seed, same inputs | **FAIL**: the SMALLER cluster yields `[±∞]` even with a clean isolated orthonormal seed, at prec 53/128/256 |
| d5n2 | single mult-2 cluster + simple rest `{0,0|1,2,3}` (n=5) | **FAIL** — so it is NOT multi-cluster |
| d5n3-A/B | single mult-2 cluster, axis-aligned AND disjoint-Givens, n=3..6 | **FAIL** in every case |
| d5n3-C | same spectrum, **dense Givens-chain** conjugation | **CERTIFIES** in every case |
| d5n3-D | whole-space cluster `k=n` (all evals equal) | **CERTIFIES** |
| d5n3-E | `k=n−1` cluster | **FAIL** |

**Root cause (pinned).** Rump's `partition_X_sorted` selects `k` "frozen" rows of
`X_approx` by magnitude; when the cluster's invariant subspace is **row-sparse**
(supported on few coordinates — the axis-aligned or disjoint-Givens projector
leaves whole rows ≈0 across all `k` columns), no well-conditioned `k`-row frozen
set exists and the Krawczyk preconditioner is singular → `[±∞]`. A `k=n` cluster
(every row used) or a dense subspace (every row nonzero) has a good partition.
This is **not** precision (fails at prec=256), **not** seed parallelism (fails
with orthonormal zheev seed), **not** inter-cluster gap (d5n2-H2: gap 1→20 does
not help).

### 2.2 The fix, measured SOUND end-to-end

Densify `A' = U A U†` (dense rational-Givens `U`), eig on `A'`, back-map
`X = U† X'`:

| probe | input | `eig_multiple(A')` | per-cluster subspace | disjoint | `Σk=n` | residual on **original A** | verdict |
|---|---|---|---|---|---|---|---|
| d5n4/5/7 | C^5 {3,2} | — | ✓ | ✓ | ✓ | `‖AX−XJ‖_ub = 1.7e-31` | **SOUND** |
| d5n5/7 | C^6 {4,2} | — | ✓ | ✓ | ✓ | `4.8e-32` | **SOUND** |
| d5n5/7 | C^7 {4,3} | — | ✓ | ✓ | ✓ | `9.3e-32` | **SOUND** |
| d5n7 | C^12 {7,5} | — | ✓ | ✓ | ✓ | `5.8e-29` | **SOUND** |
| d5n7 | C^16 {9,7} | — | ✓ | ✓ | ✓ | `1.4e-27` | **SOUND** |
| d5n6 | `{0^4|1|2|3}` rank-3 Choi | **OK** | ✓ | ✓ | ✓ | — | **SOUND** |
| d5n6 | `{0^3|2,2|5}` block-alg | **OK** | ✓ | ✓ | ✓ | — | **SOUND** |
| d5n6 | `{0^3|4^3|9}` | **OK** | ✓ | ✓ | ✓ | — | **SOUND** |
| d5n6 | `{0^3|7,7}` two-clust-no-anchor | **OK** | ✓ | ✓ | ✓ | — | **SOUND** |
| d5n6 | `{0^3|1,1}` (the §D5n killer) | **OK** | ✓ | ✓ | ✓ | — | **SOUND** |

λ-ball radii ~1e-32 at prec=128 (~1e-15 at prec=53); `U U† − I` certified to
1.8e-37 (n=5) … 4.6e-37 (n=7). Timing of the whole densify+eig chain: 1 ms
(n=5) … 28 ms (n=16) — negligible. The certified projector `Π_M` from the
non-orthonormal `X`: `‖Π²−Π‖_ub = 1.5e-31`, `‖Π−Π†‖_ub = 1.0e-31`,
`‖Π−Π_true‖_ub = 5.1e-32` (d5n-pim).

### 2.3 The exact recipe (decision)

**CHOSEN ROUTE: densify-then-per-cluster-Rump (§1.2).** No case remains a
documented limitation in the project's regime; the §D5n fail-loud teeth become a
*defence in depth* (R3) rather than the expected path. The increment-1 fail-loud
on `eig_multiple(A)`-returns-0 is *removed in favour of* densify-and-retry
(§3.1); only a genuinely unresolvable (overlapping-balls / non-finite at this
prec) input still aborts.

---

## 3. API contract (for the implementation subagents)

Conventions: every routine takes explicit `prec`; outputs are caller-initialised
unless the docstring says the routine `init`s them; arb/acb storage owned by the
caller; fail-loud (Rule 4) on any unresolved certification. `≤200 LOC/file`
(Rule 10) ⇒ the file split below.

### 3.1 New substrate: certified invariant subspaces — `src/aic_mat_eigvec.c` (NEW)

Header additions to `include/aic/aic_mat.h`. A small result struct mirrors the
per-cluster output:

```c
/* One certified eigen-cluster of a Hermitian matrix: a real eigenvalue ball and
 * an n×k certified invariant-subspace enclosure (k = multiplicity). X is NOT
 * orthonormal (FLINT Rump does not normalise); use aic_mat_cluster_projector to
 * get the orthogonal range projector. */
typedef struct {
    acb_t lambda;     /* certified eigenvalue ball (imag ∋ 0 asserted) */
    acb_mat_t X;      /* n×k certified invariant-subspace enclosure (A X ⊆ X) */
    slong k;          /* multiplicity (cluster size) */
} aic_mat_eigcluster;

/* Certified eigenvalue + INVARIANT-SUBSPACE decomposition of a Hermitian H,
 * degeneracy-robust via dense-unitary densification + per-cluster Rump
 * (acb_mat_eig_enclosure_rump). Resolves the §D5n two-clusters-each-≥2 wall
 * (measured). gap_thr clusters the seed eigenvalues (§1.5); pass <=0 for the
 * default. OUTPUT: *clusters is allocated HERE to *n_clusters entries (each
 * lambda/X init'd); free with aic_mat_eigcluster_free. ASSERTS (Rule 4):
 * Hermiticity; every per-cluster Rump enclosure finite (else "clusters
 * unresolved at prec"); cross-cluster lambda balls disjoint (else
 * "clusters overlap at prec"); Σ k == n. Sorted ascending by eigenvalue. */
void aic_mat_eig_hermitian_subspaces(aic_mat_eigcluster **clusters,
                                     slong *n_clusters, const acb_mat_t H,
                                     double gap_thr, slong prec);

void aic_mat_eigcluster_free(aic_mat_eigcluster *clusters, slong n_clusters);

/* Orthogonal projector Π = X (X†X)⁻¹ X† onto a certified invariant subspace
 * (the cluster's range). out must be init'd n×n. ASSERTS X†X certified-invertible
 * (acb_mat_inv != 0) — fail loud if the subspace enclosure is rank-deficient. */
void aic_mat_cluster_projector(acb_mat_t out, const aic_mat_eigcluster *cl,
                               slong prec);
```

Internal helper (static in `aic_mat_eigvec.c`, or shared via
`aic_mat_internal.h` if reused): `aic_mat_dense_unitary(acb_mat_t U, slong n,
slong prec)` building the rational-Givens densifier (§1.3).

**Increment-1 hardening (same file or `aic_mat_eig_multiple.c`):**
`aic_mat_eig_hermitian_multiple` (eigenvalue balls only) should gain the
densify-retry: on `acb_mat_eig_multiple(A)==0`, retry on `A' = U A U†` (the
spectrum is conjugation-invariant, so the balls are returned as-is) before
aborting. This retires the §D5n fail-loud for the *eigenvalue/rank* layer too
(carrier rank on real `⊕B(L_j)` inputs no longer aborts). Keep the abort only
for genuinely-unresolvable-at-prec inputs. **Mutation-test the retry** (force
the raw path to 0, confirm the densified path certifies).

> File-size note: `aic_mat_eig_multiple.c` is 151 LOC; adding the densify-retry
> (~25 LOC + the shared `dense_unitary` ~20 LOC) keeps it < 200. The subspace
> struct + the two new public routines (~120 LOC) go in the NEW
> `src/aic_mat_eigvec.c`. Put `aic_mat_dense_unitary` in `aic_mat_internal.h` so
> both files share it.

### 3.2 Certified Choi→Kraus — `src/aic_ucp_kraus_arb.c` (NEW)

Mirror `src/aic_ucp_latd.c` EXACTLY (the Convention-A conjugate reshape and the
`_tol` PSD-cone variant), but on the certified arb path via §3.1.

```c
/* Certified Choi→Kraus (arb counterpart to aic_ucp_choi_to_kraus_latd). PSD
 * eigendecomposition of C (Convention A) via aic_mat_eig_hermitian_subspaces;
 * for each cluster with λ_c CERTIFIED above keep_thr, the orthonormalised range
 * basis columns v give Kraus ops K_a[i,c] = sqrt(λ_a)·conj(v_a[i·dim_H + c])
 * (the conjugate reshape, file-docstring Convention A). keep_thr = an arb ball
 * (dim_K·dim_H)·2^-52·‖C‖_F matching the _latd gate. ASSERTS (Rule 4): C
 * Hermitian; every kept λ_c certified-above and every dropped certified-below
 * keep_thr (a STRADDLING cluster ⇒ rank unresolved ⇒ abort, like
 * aic_mat_certified_rank); a cluster λ certified < −neg_tol ⇒ non-CP abort.
 * phi is OUTPUT (aic_ucp_kraus_init'd HERE; caller clears). */
void aic_ucp_choi_to_kraus_arb(aic_ucp_kraus *phi, const acb_mat_t C,
                               slong dim_K, slong dim_H, slong prec);

/* TOLERANCE-AWARE PSD-cone variant (FINDINGS §C14), arb counterpart to
 * aic_ucp_choi_to_kraus_latd_tol. Three certified regimes per cluster:
 *   λ_c certified > keep_thr            → KEEP (sqrt-scaled conj-reshape Kraus);
 *   λ_c certified ∈ (−neg_tol,keep_thr] → DROP (cone-defect/noise); if its mid
 *                                         < 0 add k_c·|mid| to *clipped_neg_out;
 *   λ_c certified ≤ −neg_tol            → FAIL LOUD (genuine non-CP).
 * A cluster ball STRADDLING keep_thr or −neg_tol is UNRESOLVED → abort (raise
 * prec). neg_tol ≥ 0; strict aic_ucp_choi_to_kraus_arb delegates with
 * neg_tol = keep_thr. clipped_neg_out (nullable) = certified cone-defect mass. */
void aic_ucp_choi_to_kraus_arb_tol(aic_ucp_kraus *phi, const acb_mat_t C,
                                   slong dim_K, slong dim_H, double neg_tol,
                                   double *clipped_neg_out, slong prec);
```

> **Subtlety the implementer must respect (FINDINGS §C14 + the `_latd`
> conjugate-reshape Convention A).** The kept eigenvector columns must be an
> ORTHONORMAL basis of the cluster's range (so `kraus_to_choi(extracted) == C`).
> Rump's `X` is **not** orthonormal — orthonormalise per cluster before the
> reshape (a certified QR / `aic_mat_cluster_projector`-derived basis, or
> Gram-Schmidt with certified `acb` arithmetic). For a *simple* nonzero
> eigenvalue `k=1`, normalise the single certified eigenvector. **Mirror the
> ascending-sweep threshold logic of `_latd` so the two paths produce the same
> Kraus rank.** `neg_tol` semantics are byte-for-byte the `_latd_tol` contract.

### 3.3 Certified carrier subspace split — extend `src/aic_ucp_carrier.c`

```c
/* Certified RANGE projector Π_M of the carrier Q = Σ K_a K_a† (lem_carrier,
 * .tex:1724): the orthogonal projector onto the support (range) of Q, i.e. the
 * sum of the cluster projectors over Q's eigenvalue clusters certified ABOVE
 * thr = dim_K·2^-52·‖Q‖_F. Built via aic_mat_eig_hermitian_subspaces +
 * aic_mat_cluster_projector. out init'd dim_K×dim_K. The certified counterpart
 * to the double-path rank (aic_ucp_carrier_rank_latd) — gives the projector,
 * not just the dimension. ASSERTS (Rule 4) no Q-eigenvalue cluster straddles
 * thr (range unresolved at prec). Cross-checks: trace(Π_M) == carrier rank;
 * Π_M Q == Q (Q supported on M); for any X with Ker X ⊇ M the cor_carrier
 * annihilate-defect is 0. */
void aic_ucp_carrier_projector(acb_mat_t out, const acb_mat_t Q, slong dim_K,
                               slong prec);
```

`aic_ucp_carrier_rank` (increment 1) is now the `trace(Π_M)` of this — keep both;
the rank delegates to `aic_mat_certified_rank` (cheaper, no subspaces) and the
projector to the subspace layer; a test asserts `round(Tr Π_M) == rank`.

### 3.4 File split summary (Rule 10)

| file | status | contents | est. LOC |
|---|---|---|---|
| `src/aic_mat_eigvec.c` | NEW | `aic_mat_eig_hermitian_subspaces`, `aic_mat_eigcluster_free`, `aic_mat_cluster_projector` | ~150 |
| `src/aic_mat_eig_multiple.c` | EXTEND | densify-retry in `aic_mat_eig_hermitian_multiple` | 151 → ~180 |
| `include/aic_mat_internal.h` | EXTEND | `aic_mat_dense_unitary` (shared densifier) | +~6 |
| `src/aic_ucp_kraus_arb.c` | NEW | `aic_ucp_choi_to_kraus_arb`, `_arb_tol` | ~160 |
| `src/aic_ucp_carrier.c` | EXTEND | `aic_ucp_carrier_projector` | 151 → ~190 |
| `include/aic/aic_mat.h`, `include/aic/aic_ucp.h` | EXTEND | the prototypes + struct above | +~40 each |
| `tests/test_eigvec.c` | NEW | §4 cross-checks | ~250 (test file, exempt) |

---

## 4. Test plan (cross-checks BEFORE code — Rules 6/7)

Write these RED first; mutation-prove each tooth (perturb impl → RED → restore).
Reuse the `tests/test_eigmult.c` fork+SIGALRM watchdog for the fail-loud teeth.

**S1 — invariant-subspace residual (rung 2/3, the core certificate).** For
`build_projector` C^5{2,3}, C^6{2,4}, C^7{3,4} (the §D5n cases) and a few
`build_conj_diag` block spectra: run `aic_mat_eig_hermitian_subspaces`, assert
every cluster's `‖H X_c − X_c J_c‖_ub` (recompute `J_c` from the returned λ, or
expose it) is `< 1e-20` at prec=128, all λ balls disjoint, `Σ k_c = n`. Print the
worst residual (expect ~1e-31). *This is the case test_eigmult.c T5 asserts FAILS
LOUD — S1 asserts the SAME inputs now CERTIFY; update T5 accordingly (it is no
longer the boundary).*

**S2 — η=0 / exact-projection oracle (rung 3).** A rank-r orthogonal projector
`P`: the certified range projector `Π_M` from the λ=1 cluster equals `P` itself
(`‖Π_M − P‖_op < 1e-25`); `trace(Π_M) == r`; the kernel cluster projector
`= I − P`. For the depolarizing Choi `(1/d)I_{d²}` (single cluster): one cluster,
`k = d²`, `Π_M = I`.

**S3 — double-vs-arb@53 Choi→Kraus (rung 2).** For a complex asymmetric channel
(reuse `test_ucp.c`'s Convention-A oracle fixture) and a multi-block
`mixconj_blocks`-style channel: run `aic_ucp_choi_to_kraus_arb` at prec=53 and
`aic_ucp_choi_to_kraus_latd`; **compare AS CHANNELS** (rebuild Choi from each,
`‖C_arb − C_latd‖_op < 1e-12`) — NOT operator-by-operator (Kraus gauge freedom,
header warning). Recovered Kraus rank must match.

**S4 — rebuild-Choi round-trip (rung 3).** `kraus_to_choi(choi_to_kraus_arb(C))
== C` within the certified ball: `‖C_rebuilt − C‖_op` enclosure contains 0 (it is
exact up to the kept/dropped tail). Use an exact-rational PSD `C` so the ball is
tight.

**S5 — PSD-cone `_tol` variant (FINDINGS §C14).** Feed a Choi with a planted
O(η²) negative cluster (mid `−2.5e-6`, the measured factorize defect): with
`neg_tol` admitting it, `_arb_tol` drops it, accumulates `clipped_neg_out`
(certified mass), and the kept Kraus rebuild a CP map. With a planted O(1)
negative (`−0.3`): `_arb_tol` FAILS LOUD ("not CP"). Mutation: flip the
`λ ≤ −neg_tol` comparison → the O(1) case must stop aborting (RED).

**S6 — carrier cross-checks.** Build `Q` = rank-r projector; assert
`aic_ucp_carrier_projector(Q)` has `trace == aic_ucp_carrier_rank(Q) ==
aic_ucp_carrier_rank_latd(Q) == r`, `Π_M Q = Q`, and the cor_carrier annihilate
defect (`aic_ucp_carrier_annihilate_defect` with `X = I − Π_M`) is `0`.

**S7 — fail-loud teeth (Rule 4, fork watchdog).**
  (a) *Straddle:* a Choi whose smallest kept eigenvalue cluster sits ON `keep_thr`
      → `choi_to_kraus_arb` aborts ("straddle"/"unresolved").
  (b) *Cross-cluster overlap:* a Hermitian with two eigenvalues separated by far
      below `2^-prec` at low prec=24 → `aic_mat_eig_hermitian_subspaces` aborts
      ("clusters overlap"). Mutation: drop the disjointness assert → child
      exits 0 (RED).
  (c) *Non-finite enclosure:* if any densified-Rump still returns `[±∞]` (none
      found in probes, but the guard must exist) → abort ("clusters unresolved").
  (d) *Non-CP:* §C14 O(1)-negative Choi → strict `choi_to_kraus_arb` aborts.

---

## 5. Risk register / open items (beads)

- **R1 — (resolved here) D5n.** Resolved by densification (§2). Bead `aic-4td`
  deliverable (3). **Update FINDINGS §D5n** from "OPEN limitation, fail-loud" to
  "RESOLVED by dense-unitary preconditioning (§2 numbers)"; update the
  `aic_mat_eig_multiple.c` docstring (the "raising prec may not help / SEED-
  conditioning limit" prose is now superseded — the cause was the Rump frozen-row
  partition on row-sparse subspaces, fixed by densification, not seed
  re-orthonormalisation). Update `tests/test_eigmult.c` T5 (the C^5{2,3} case now
  CERTIFIES; keep a genuine unresolvable-at-prec fixture as the fail-loud tooth —
  e.g. S7(b) sub-`2^-prec` gap at prec=24).
- **R2 — orthonormalisation of Rump X for the Kraus reshape.** Rump `X` is not
  orthonormal; the Choi→Kraus reshape and `kraus_to_choi` round-trip REQUIRE an
  orthonormal cluster basis or the Choi will not rebuild. Spec'd (§3.2) but the
  implementer must verify the certified QR/Gram-Schmidt does not blow up the ball
  for a near-defective `X†X`. Probe it during impl (cross-check S4). *Low risk*
  (measured `X†X` well-conditioned, `acb_mat_inv` succeeded), but it is the one
  numerical step not yet end-to-end-probed for the Kraus path specifically.
- **R3 — densifier adequacy (defence in depth).** The fixed rational-Givens `U`
  densified every probed spectrum. A pathological input left row-sparse by this
  exact `U` is not impossible in principle; the routine FAILS LOUD (non-finite
  Rump enclosure) rather than returning garbage, so it is sound. *Mitigation if
  ever hit:* retry with a second angle or a random orthogonal `U` (gated by the
  same finite-X assert). File as a follow-up bead only if a real input triggers
  it. Keep the fail-loud tooth (S7c).
- **R4 — prec floor for the disjointness gate.** Two genuinely-distinct but
  O(2^-prec)-close eigenvalues will overlap and abort; raising prec resolves a
  true gap (unlike the old D5n seed limit which prec could not fix). This is
  correct fail-loud behaviour, not a wall — but the carrier `keep_thr`/`gap_thr`
  must be chosen so the **zero cluster vs smallest nonzero** gap is resolvable at
  the working prec (§1.5). Document the gap→prec relation in `ALGORITHM.md`.
- **R5 — deferred hostile review of increment 1** (bead `aic-4td` deliverable
  (1)): scrutinise increment-1's `arb_gt/arb_lt` straddle rigor, the T1
  `eigset_close` predicate, the threshold-ball radius honesty, and run the
  certified carrier rank on a REAL `⊕B(L_j)` block-algebra Choi (the §D5n
  boundary) — which the densify-retry (§3.1) should now make pass. This review is
  IN SCOPE for the increment and must run before the increment-2 code lands
  (Rule 9, Core change).
- **R6 — `ALGORITHM.md`** has no eig/subspace section yet; add the densify-Rump
  narrative + the §1.6 rigor argument cited to `.tex:1724` (lem_carrier) and the
  FLINT Rump semantics when the code lands (Rule 11).

### Provenance / citations

- FLINT semantics: `/usr/include/flint/acb_mat.h:417-432`; FLINT 3.0.1
  `acb_mat/eig_enclosure_rump.c` (Rump 2010 §13.4), `flintlib.org/doc/acb_mat.html`.
- Convention A conjugate reshape + `_tol` PSD-cone: `src/aic_ucp_latd.c:1-149`,
  `include/aic/aic_ucp.h:23-51,205-251`, FINDINGS §C14.
- Carrier `Q = Σ K_a K_a†` = range of `Tr_F(VV†)`: `.tex:1724` (lem_carrier),
  `src/aic_ucp_carrier.c:56-73`.
- Left-vs-right range hazard (N/A for Hermitian): FINDINGS §C4.
- Increment-1 eigenvalue/rank layer + §D5n: `src/aic_mat_eig_multiple.c`,
  FINDINGS §D5/§D5n, `tests/test_eigmult.c`.
- All numeric claims in §2: bounded probes run 2026-06-01 against FLINT 3.0.1
  (libflint18t64 3.0.1-3.1build1), prec ∈ {53,128,256}, `build_projector`
  rational-Givens fixtures + diagonal block spectra; probes deleted after
  measurement per probe-hygiene.
```
