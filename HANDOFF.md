# HANDOFF.md — almost-idempotent-channels

Five-minute orientation for a fresh agent. Last updated 2026-05-28.

## What this project is

Implement Alexei Kitaev, *Almost-idempotent quantum channels and approximate
$C^*$-algebras* ([arXiv:2405.02434](https://arxiv.org/abs/2405.02434)) as
**constructive finite-dimensional algorithms**. Architecture: tight C with two
numeric paths — **FLINT/arb** (certified arbitrary precision, the ground truth)
and **LAPACK double** (fast, uncertified) — later wrapped by a Julia `ccall`
package. The paper proves things non-constructively (Lefschetz–Hopf, Haar
diagonals, holomorphic calculus); our job is constructive finite-dim algorithms
whose outputs meet the paper's `O(ε)`/`O(η)` bound, which the arb path certifies.

## Read order (do this first)

1. `CLAUDE.md` — the contract: Laws 1–4 (ground-truth-by-citation; reuse FLINT/
   LAPACK; constructivize-don't-transcribe; **audition-don't-presume + the
   Pareto rule**), domain hallucination callouts, the paper structure map.
2. `MODULE_PLAN.md` — the module DAG, result→module map, build order.
3. `paper/shards/shard-*.md` — per-section analysis (statement+proof+
   constructivization). `paper/src/approximate_algebras.tex` is canonical ground
   truth (cite its line numbers in code).
4. `docs/adversarial/{README,nla,domain}.md` — the punishing-instance catalog
   (testing is adversarial-first; the happy path is a liability).

## How to build & test

```bash
make test     # 6 test binaries, all green, zero warnings (strict flags)
make bench    # ns/op microbenchmarks (double-vs-arb head-to-heads)
make clean
```
Deps: gcc/C11, FLINT 3.0.1 (arb/acb bundled), LAPACK+LAPACKE+BLAS (installed
2026-05-28). `make` links `-lflint -llapacke -llapack -lblas -lm`.

## What's built (all green, committed)

- **E0-infra epic — DONE.** `Makefile` (dual path), `include/aic.h`, the audition
  harness `tests/aic_test.h` (fail-loud `AIC_CHECK*`, `aic_acb_close`,
  `AIC_CHECK_DOUBLE_CLOSE`/`EIGSET_CLOSE` for double-vs-arb@53) and
  `bench/aic_bench.h` (`AIC_BENCH` → `BENCH … ns/op=…`).
- **`mat`** (`aic_mat_*`, arb substrate): Frobenius/operator norm, Hermitian eig
  (certified, **simple-spectrum only**), singular values, Kronecker, partial
  trace (left-major convention, matches tex:263–288). 14 cross-checks.
- **`funcalc`** (`aic_funcalc_*`): `sgn` (Newton–Schulz default + Denman–Beavers,
  **eig-free** so degeneracy-proof), `abs`, `theta`, `prop_P` (tex:524–533),
  `xpow`. 10 cross-checks incl. degenerate projector.
- **`contraction`** (`aic_contraction_*`): `lem_invfun` solver (tex:564–592),
  Picard + Anderson, callback-based, fail-loud guards. 30 cross-checks incl. the
  tex:594–599 sgn-system vs `aic_sgn` to 1e-25.
- **`latd`** (`aic_latd_*`, LAPACK double path): `LAPACKE_zheev` (Hermitian eig,
  **handles degenerate spectra**), `LAPACKE_zgesvd` (opnorm/SVD). 45 adversarial
  cross-checks; double-vs-arb@53 restored. 30×–370× faster than arb.

## Key decisions & findings (the non-obvious stuff — don't relearn the hard way)

- **Two complementary paths, neither dominates (Pareto).** Double(LAPACK): fast,
  uncertified, **handles degeneracy**. Arb(FLINT): certified `O(ε)` balls,
  **simple-spectrum only**. Keep both, dispatch per regime.
- **FLINT's certified eig requires a SIMPLE spectrum** — it aborts on repeated
  eigenvalues. But this project is *dominated by degenerate spectra* (projections
  have eigenvalues 0/1 with multiplicity; target algebras are `⊕ B(L_j)`). So:
  - certified degenerate eig is an open P1 (`aic-w4o.1`, needs
    `acb_mat_eig_multiple_rump`/cluster enclosures);
  - meanwhile LAPACK `zheev` is the fast uncertified route for degenerate inputs;
  - **functional calculus is eig-free** (Newton–Schulz), so `sgn`/`abs`/`theta`/
    projection-extract sidestep the whole problem — prefer those routes.
- **Audition methodology (Law 4 + Pareto).** No algorithm is presumed fit; every
  choice is red-green TDD auditioned on correctness AND performance, candidates
  kept alive. When none dominates, add perf dimensions (wall time, iterations,
  memory, accuracy@prec, robustness) and keep the Pareto frontier.
- **Testing is adversarial-first.** Correct behavior on an evil input = the bound
  provably holds OR the routine fails loud — never plausible-but-wrong. The
  `mat` eig bench's well-separated spectrum is the *anti-pattern* (flatters the
  easy case). See `docs/adversarial/`.
- **The dimension-sweep is the universality canary** (`aic-dbo.3`): the paper's
  `O(ε)` constant is universal (dim-independent). A constant `∝ dim A` (e.g. from
  the naive Haar diagonal instead of the explicit B-diagonal, tex:484 vs
  tex:1249–1254) is invisible on any single instance — only a dim sweep catches
  it. This is the single highest-value correctness test.
- **η=0 oracle** is the zero-defect anchor: an exactly idempotent Φ must drive
  the whole pipeline to the exact Choi–Effros / `th_idemp_structure` structure
  with zero defect. Every evil sweep should reduce to it continuously as knob→0.
- **`factorize` and `cstar_build` are the most at-risk modules** (end-of-pipeline,
  outline proof, near-singular inversions, the composite constant).

## What's next (ready work — `bd ready`)

Issue tracker is **beads**, prefix `aic` (persistent across sessions; JSONL at
`.beads/issues.jsonl` is committed). 45 issues, 9 closed.

- **`T-ucp` (`aic-c7n`, P1)** — the next module; completes E1-foundations. UCP map
  representations (Choi/Stinespring/Kraus) + conversions + cb-norm + carrier
  (shard G). Kraus extraction = eig of the (degenerate-prone) Choi matrix → use
  the `latd` double path / certify later.
- **`aic-w4o.1` (P1)** — certified degenerate Hermitian eig (arb).
- **`aic-dbo.2` (P1)** — implement the adversarial instance generators (the
  evil-matrix corpus; unify with `aic-f9u.1`).
- Then **E2-eps** (`aic-vs9`): `ecstar`, `unitfix`, `unitary`, `projection`
  (blocked by `aic-w4o.1`), `corner`; **E3-mainthm**; **E4-headline**;
  **E5-julia**; **E6-research** (`aic-1bc` c_0, `aic-1sk` factorization closure,
  etc.); **E7 = `aic-dbo`** (adversarial suite).
- Housekeeping: `aic-w4o.4` (split oversized test files: test_mat/funcalc/
  contraction/latd all > 200 LOC); `aic-w4o.2` (full arb SVD); `aic-w4o.3`
  (opnorm power-iteration audition); `aic-f9u.1` (shared corpus); `aic-f9u.2`
  (multi-dim Pareto bench reporting).

## Orchestration mode (how this was run, to continue it)

The user runs this serially: each module is delegated to a **subagent** (Opus
for coding, Sonnet for research/summarization), the orchestrator **verifies the
build independently** (don't trust the report — re-run `make test`), files beads
for issues that arise, commits atomically per closed bead with `.tex`-line
provenance, and reports at epic milestones. **No parallel Julia.** `TaskCreate`
is fine for ephemeral in-session tracking; beads is the only persistent tracker.
**No GitHub CI; no remote configured** (commit locally, don't push).

## Conventions / gotchas

- Cite `approximate_algebras.tex:<line>` + the verbatim equation in every math
  routine. `≤200 LOC` per source file (several test files violate this →
  `aic-w4o.4`). Arb path never links LAPACK. Fail-loud on precondition violation
  / LAPACK `info!=0` / lost precision. `acb_mat_eig_simple` needs a simple
  spectrum; `LAPACKE_zheev` does not. Conversions in `latd` are row-major.
