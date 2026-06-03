# Adversarial / correctness stress suite — index

**Internal developer notes — not user documentation.** This is the design catalog
of the adversarial test corpus (the families of "evil" inputs the suite is built
around); the implementing test code lives in [`../../tests/`](../../tests/). For
the user documentation, start at [`../src/index.md`](../src/index.md).

> Directive (2026-05-28): the happy path is the least useful. Every benchmark
> and correctness test is built around punishing, adversarial instances.
> Correct behavior on an evil input = the certified bound provably holds, OR the
> routine fails loud (asserts/refuses) — never a plausible-but-wrong answer.

Bead track: epic `aic-dbo`. Feeds the shared benchmark corpus (`aic-f9u.1`).

## Catalogs
- [`nla.md`](nla.md) — general numerical-linear-algebra evil (8 families):
  Jordan/defective (incl. the paper's `t^{1/3}` example, tex:540), non-normal /
  pseudospectra (Davies, Grcar), rank-deficient / near-singular, degenerate &
  clustered spectra, ill-conditioned / graded (Hilbert, Kahan), gallery torture
  (Frank, companion, Toeplitz), near-boundary `‖X²−I‖→1`, arb ball-explosion.
- [`domain.md`](domain.md) — operator-algebra / channel evil (6 families):
  UCP/Choi (extremal, cb-vs-op gap, near-degenerate carrier), `η→1/4⁻`,
  ε-C\* (concentrated associativity defect, near-`C*`-violation, near-trivial
  projections, the `ε~c/n` blowup regime), projection/corner, δ-hom/main-loop,
  factorization composite bound. Plus the `η=0` oracle and a sweep matrix.

### Increment 2, tranche 1 (landed 2026-06-01, bead aic-dbo.2)

The first two channel/UCP-map generators are in `tests/aic_adversarial_domain.c`.
Self-test count: 65 (up from 33 after Increment 1). Both mutation-proven.

- **1B `chan_cb_op_gap`** (tex:366-388): measure-prepare cb-vs-op gap. Closed-form
  defect `||Phi^2-Phi||_cb = eta*sqrt(1-eta)`. Measured anchor (prec=256, d=2):
  eta=0.20 gives 0.178885, certified bracket [0.126491, 0.505964]. eta=0 reduces to
  exact dephasing (bracket [0,0]).
- **2A `chan_depol_boundary`** (tex:516-525): depolarizing, rho(Phi^2-Phi) = p(1-p),
  hitting 1/4 at p=1/2. Lower bound lo = p(1-p)*sqrt(3)/2 at d=2 (ratio lo/[p(1-p)]
  constant = 0.866025 across all p). Measured anchor (prec=256, d=2): p=0.5 gives
  bracket [0.216506, 0.866025], the exact basin edge. p in {0,1} reduce to exact
  idempotents (bracket [0, ~3e-16]).

### Increment 2, tranche 2-partial (landed 2026-06-01, bead aic-cxo)

Two more channel generators; self-test count 84.

- **1D `chan_unital_defect`** (tex:432/672): single Hermitian Kraus
  `K_0 = diag(sqrt(1+delta_u), sqrt(1-delta_u), 1,..)`, so `Phi(I) = I + delta_u*E`
  (E traceless, ||E||=1) and the certified unital defect equals `delta_u` exactly
  (measured == delta_u to ~1e-12 at delta_u in {1e-3,0.1,0.5}; delta_u=0 unital).
- **1C `chan_carrier_dropout`** (tex:1724/1731, new file `aic_adversarial_carrier.c`):
  single Hermitian Kraus `K_0 = diag(1,..,1,sqrt(gap))`, carrier `Q = sum K K^dag =
  diag(1,..,gap)`. Certified carrier rank = d for gap >> thr (gap in {0.5,1e-3,1e-6});
  smallest carrier eig == gap; gap=0 oracle certifies rank d-1. NOTE (hostile review):
  a DIAGONAL carrier only drives the *certify* path (point-ball eigenvalues, never
  straddles); the straddle/fail-loud hard case needs a densified carrier (covered at
  the eigvec level by test_eigvec S6b) — a channel-level densified generator is the
  follow-up bead aic-v5f.

Remaining inc-2 families (tranche 3): 3D eps~c/dim universality canary (unblocks
aic-dbo.3); a DENSIFIED near-degenerate carrier `U diag(1,..,gap) U^dag` (drives the
straddle/fail-loud path + is convention-sensitive — closes the 1C diagonal-carrier
gap); non-commutative calibrated eta->1/4; family 2B; the make_mixconj
corpus-unification.

## Most lethal instances (merged shortlist)
1. **Exact-degeneracy projector** (nla 4b) — eigenvalues {0,1} highly repeated;
   FLINT's certified *simple*-spectrum eig cannot isolate it. Proves the eig-free
   `sgn` (Newton–Schulz) is the only viable route; any regression to eig fails.
2. **`η → 1/4⁻` regularization boundary** (domain 2A) — `θ(2Φ−1)` series
   near-divergent; doubles cannot detect the blow-up. Attacks funcalc +
   contraction + assoc_ecsa at once.
3. **`ε ~ c/dim` blowup regime** (domain 3D) — the canary for the paper's
   *universality* claim: if any constant (esp. `c_0` from `cor_improvement`)
   grows with dimension, only a **dimension sweep** reveals it. Distinguishes the
   incremental construction from the naive Haar route (tex:484).
4. **`t^{1/3}` non-normal** at `t=1e-6` (nla 1a, tex:540) — the paper's own
   fragility warning; tests the `sgn` domain guard at the edge.
5. **Sub-ULP `‖X²−I‖ = 1` straddle** (nla 7a) — the guard must conservatively
   abort when the ball straddles the domain boundary; doubles are unconditionally
   wrong here.
6. **cb-vs-op ampliation-sensitive idempotence defect** (domain 1B) — `‖Φ²−Φ‖`
   tiny but `‖·‖_cb` large; defeats naive operator-norm idempotence estimates.
7. **Near-trivial projections `P ≈ ±I`** (domain 3C) — Route A spectral split has
   no usable gap; `lem_nontriv_projection` must still deliver or refuse.
8. **Worst-case composite `O(η)`** across the full factorization (domain 6A) —
   forces the unextracted composite constant into the open.

## Cross-cutting findings (drive the bead track)
- **`sgn`'s `‖X²−I‖<1` domain guard is the fragility chokepoint** — attacked by
  all NLA families; it depends on a certified `opnorm` that itself degrades on
  non-normal/ill-conditioned inputs. Harden and stress it first.
- **The dimension sweep is the universality canary.** A constant `∝ dim A` (e.g.
  from accidentally using the Haar diagonal instead of the explicit B-diagonal,
  tex:1249–1254) is invisible on any single instance. `aic-dbo.3`.
- **`factorize` and `cstar_build` are the most at-risk modules** (end-of-pipeline,
  outline proof, multiple near-singular inversions, composite constant).
- **The `η=0` exact-idempotent oracle is the zero-defect anchor** the evil sweeps
  perturb away from — every sweep should reduce to it continuously as the knob → 0.
