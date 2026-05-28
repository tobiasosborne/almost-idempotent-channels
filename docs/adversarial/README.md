# Adversarial / correctness stress suite — index

> Directive (2026-05-28): the happy path is the least useful. Every benchmark
> and correctness test is built around punishing, adversarial instances.
> Correct behavior on an evil input = the certified bound provably holds, OR the
> routine fails loud (asserts/refuses) — never a plausible-but-wrong answer.

Bead track: epic `aic-dbo`. Feeds the shared benchmark corpus (`aic-w4o.5`).

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
