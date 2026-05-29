# HANDOFF.md — almost-idempotent-channels

Orientation for a fresh agent. Last updated **2026-05-29**, after the
**`ecstar` session** (the ε-$C^*$ algebra data model + axiom-defect estimators,
bead `aic-knm`, built in two committed increments). Current state: `master`
clean, **14 test binaries green, zero warnings**, 17 beads closed of 57. Prior
session built the η-defect cb-norm stack (`aic-m24` CLOSED). The η=0
vertical-slice milestone (`aic-9kk`) was achieved earlier. Read §"Channel-module
conventions" before touching `ucp`/`idemp`/`cbnorm` or building `assoc_ecsa`;
those conventions are load-bearing and prior hostile reviews learned them hard.

**⚠️ RESUME HERE — `aic-knm` (`ecstar`) is CLOSED; next is `aic-92f`
(`assoc_ecsa`), the headline step.** The session ended mid-orchestration when an
internet outage made subagent spawns unsafe (a half-spawned agent had already
left the tree dirty twice). `aic-knm` is complete: Cycles 1+2 committed & green,
Cycle 2's hostile review done inline (bead `aic-b7c`, CLOSED — all 7 surfaces
PASS, no blockers, the canary mutation-proven to fire on a d^{3/2} injection).
The next step toward the paper's headline is **`aic-92f` (`assoc_ecsa`)** —
Φ̃=θ(2Φ−1) on the superoperator `2Φ−1` (reuse `funcalc` θ), A=ImgΦ̃, the
approximate Choi–Effros star, `th_almost_idemp` (extended O(η)-$C^*$), proofs at
`.tex:2208/2239/2639`. It is unblocked (all four deps closed) and able to
MEASURE/CERTIFY η (`eta_idempotence` / `aic_cbnorm_certify`). assoc_ecsa is a
large module — START IT WITH A STABLE CONNECTION (research leg first: θ on the
non-Hermitian superoperator `2Φ−1`, its CP-ness / Hermicity, the Φ_assoc1/2
identities). It will also finally exercise the deferred `aic-3qq` (Π_A teeth, an
η>0 A) and `aic-4c7` (real O(η) canary). The `ecstar` estimators (`aic_ecstar_*`)
are the verification tooling assoc_ecsa feeds.

**`ecstar` (this session, `aic-knm` Cycles 1+2, committed `124ed68`/`4c623bf`).**
A = subspace of $M_n$ as an orthonormal operator basis {B_k} (reshaped from
`aic_idemp_decomp.Delta`) + a borrowed UCP map Φ; star X⋆Y=Φ(XY) (`aic_ucp_apply`,
NOT plain XY); norm = operator norm; unit I=1_n; all inherited from B(H)
(verified vs `.tex:2186-2215`). Two audition rungs for the defect estimators
(assoc/submult/cstar/involution/unit), per the Pareto discipline:
- **Cycle 1 — certified-arb BASIS-SWEEP** (`src/aic_ecstar{,_assoc,_defect,_involution}.c`):
  max over basis tuples, the cheap EXACT zero-detector for the η=0 oracle (machine-
  zero on 7 idemp channels incl. a complex-conjugated ASYMMETRIC one) + always-zero
  invariants (involution/unit are structural for HP unital Φ). Hostile review caught
  the signature "test that can't fail" twice (involution no-op; transposed reshape) —
  both fixed (shared generic-apply core; asymmetric channel) and re-mutation-proven.
- **Cycle 2 — faithful HOPM worst-case** (`src/aic_ecstar_{hopm,iterate,search,setup,certify}.c`):
  scale-invariant alternating max over the OPERATOR-norm ball (closed-form u,v power
  step + Π_A(polar(C)) block step with a monotone-ascent accept-guard), double-path
  search → arb-certified LOWER bound on the witness. The **universality canary**
  (HOPM defect/t flat across dim_A=8→45, vs the basis-sweep which drifts) is the
  payoff: HOPM avoids the d^{3/2} Frobenius trap the basis-sweep pays.
- Findings (load-bearing): submult is structurally ≤0 for ANY UCP star; the η=0
  oracle is structurally BLIND to whether Φ is in the star (only the perturbation
  teeth pin the Choi-Effros product); the Π_A accept-guard's CORRECTNESS half is
  untestable until a η>0 (non-polar-closed) A exists → bead `aic-3qq` (needs
  `aic-92f`). Cycle 3 (certified SDP UPPER bound, Watrous bilinear + SOS trilinear)
  = bead `aic-0at`. Research: `docs/research/ecstar_{paper,web}_leg.md`.

**cb-norm / η-defect stack (this session, `aic-m24` CLOSED + `aic-cne`).** The
central quantity η = ‖Φ²−Φ‖_cb now has a full certified pipeline:
- **eig-free bracket** `[‖J‖_F/n, 2‖J‖_F]` (`src/aic_cbnorm.c`) — always-valid,
  no solver, no eig (dodges `aic-w4o.1`); the fallback near η=0.
- **value entry point** `eta_idempotence(kraus)` (`julia/AlmostIdempotentChannels.jl`,
  ccalls the `libaic.so` shim `aic_ucp_choi_diff_d` then solves the Watrous SDP via
  Convex+MOSEK) — downstream modules MEASURE η through this.
- **tight certified-arb ball** `aic_cbnorm_certify` (`src/aic_cbnorm_certify*.c`) —
  rigorous `[lo,hi]`, widths ~1e-12 (~10¹²× tighter than eig-free), from two MOSEK
  feasible points (MAX primal + MIN dual solved separately in Julia, `src/sdp.jl`)
  via arb feasibility-restoration (convex-combination LOWER toward the Slater
  center; eigenvalue-shift UPPER; one-sided `herm_max_eig` PSD certs; dispatch +
  fail-loud). Hostile-reviewed (no blockers); all load-bearing choices
  mutation-proven (`tests/test_certify.c` + `tests/test_certify_teeth.c`).
- **independent oracle** `tools/diamond_oracle.wls` (wolframscript / Mathematica
  `SemidefiniteOptimization`, complex-native) confirms the MOSEK golden master to
  ≤1.1e-7. KEY: Mathematica SDP is **machine-precision only** (no `WorkingPrecision`)
  — certified high-precision η must come from the FLINT/arb route, not Mathematica.
- Design contract: `docs/cbnorm_tight_certifier.md` (TIB-grounded: Jansson SIAM
  2008, Watrous 1207.5726/TQI Ch.3, QETLAB). Follow-up `aic-ssu` (Julia end-to-end
  `certify(kraus)→(lo,hi)` wrapper). **Partial-trace direction in the dual is the
  load-bearing convention** (`aic_mat_partial_trace_right` / Convex sys 2, the
  MINOR/input factor) — pinned empirically by the asymmetric paper-example anchor;
  the design doc's first draft had it wrong.

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
make test     # 13 test binaries (mat/funcalc/contraction/harness/latd/smoke/
              #   ucp/idemp/ucp_d24/cbnorm/shim/certify/certify_teeth), all green
make lib      # build/libaic.so (shim symbols for the Julia ccall package)
              # zero warnings throughout (strict flags)
make bench    # ns/op microbenchmarks (double-vs-arb head-to-heads)
make fixtures # OPTIONAL: regenerate tests/fixtures_d24.inc.h via Julia+MOSEK
              #   (the committed .inc.h means `make test` does NOT need Julia)
make clean
```
Deps: gcc/C11, FLINT 3.0.1 (arb/acb bundled), LAPACK+LAPACKE+BLAS. `make` links
`-lflint -llapacke -llapack -lblas -lm`. For `make fixtures` only: Julia
(`~/.juliaup`) with the committed `julia/env` (Convex 0.16 + Mosek 11.2 +
MosekTools), MOSEK 11.1 at `~/mosek` + license `~/mosek/mosek.lic`. **No Python
anywhere** (user rule — Julia or C only; consumer Python bindings come later).

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
- **`ucp`** (`aic_ucp_*`, `aic-c7n`): UCP maps Φ:B(K)→B(H) as Kraus ops (Heisenberg
  `Φ(X)=ΣK_a†XK_a`); Kraus↔Choi (Convention A, `C=ΣE_ij⊗Φ(E_ij)`, K-major)↔
  Stinespring; CP cert via degeneracy-robust `herm_max_eig(−C)`; unital (two ways);
  carrier `Q=ΣK_aK_a†` + cor_carrier + PhiX_M; cb-norm CLOSED FORM `‖Φ(1)‖` (=1 UCP).
  373 checks. A hostile review caught a Choi-conjugation convention bug (now pinned
  by a `C[i,j]==Φ(E_ij)` oracle test). `‖Φ²−Φ‖_cb` is NOT here → `aic-d24` (SDP).
- **`idemp`** (`aic_idemp_*`, `aic-wuh`): `th_idemp_structure` for EXACT idempotent Φ
  (the η=0 oracle, tex:2055–2091). Double path extracts M=range(Q), A=ImgΦ (column-
  image SVD); arb path certifies the 5 relations (ΓC_MΔ=1_A, ΔΓC_M=Φ, w *-hom,
  block-diag, Γ-CP) machine-zero across 6 channels incl. M⊊H. 76 checks.
  **Milestone `aic-9kk` (η=0 vertical slice) ACHIEVED.**
- **`cbnorm`/η-defect** (`aic_ucp_compose`+`aic_ucp_choi_diff`, `aic-d24` increment 1):
  `‖Φ²−Φ‖_cb` (η). Λ=Φ²−Φ is NOT CP → needs the Watrous diamond-norm SDP. Route B:
  C core supplies `Choi(Φ²−Φ)` (both paths, cross-checked); the **SDP is in
  Julia+MOSEK** (`tools/gen_fixtures_d24.jl`, committed `julia/env`), NO Python.
  Load-bearing: `‖Λ‖_⋄ = (2/n)·SDP_optval` (Convention-A Choi has trace n). Golden
  master triple-verified (hand-derivation + Julia regen + independent CLARABEL)
  across n=2,3,4 + the paper's `η√(1−η)` example; 17 fixtures (incl. a complex one
  guarding Choi-conjugation, + a rank-deficient Λ); `test_ucp_d24` n=71. Certified-arb
  cb-ball + reusable value entry point → `aic-m24`.

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
- **The recurring failure mode is a test that *cannot fail*, not wrong math.**
  All three channel-module hostile reviews caught one: a Choi matrix that was the
  conjugate of the cited convention but passed every transpose-invariant check
  (`ucp`); a "Γ is CP" certificate that never read the stored Γ (`idemp`); a
  golden-master corpus of only real Kraus that could not catch a conjugation bug
  (`cbnorm`). **Always run a hostile review of a new core module, and demand each
  test be mutation-proven** (perturb the impl → RED → restore). The math is usually
  right; the test's *teeth* are what's missing.

## Channel-module conventions (ucp / idemp / cbnorm) — LOAD-BEARING, don't relearn

These pin the data model the headline path is built on. Verify against the `.tex`,
but start from here.

- **Heisenberg picture (observables), pinned.** A UCP map is `Φ: B(K)→B(H)`,
  `Φ(X)=Σ_a K_a† X K_a`, with `K_a: H→K` stored as a `dim_K × dim_H` `acb_mat`
  (`aic_ucp_kraus`). Unital ⇔ `Σ_a K_a†K_a = 1_H` (tex:1571). The state/Schrödinger
  dual is `T=Φ*`, `T(ρ)=Σ_a K_a ρ K_a†` — getting the dual backwards silently
  transposes everything (CLAUDE.md callout).
- **Choi Convention A (the bug that bit).** `C_Φ = Σ_ij E_ij ⊗ Φ(E_ij)`, K-factor
  MAJOR/LEFT: `C[i*dim_H+a, j*dim_H+b] = Φ(E_ij)[a,b] = Σ_x conj(K_x[i,a])·K_x[j,b]`
  — **the conjugation is on the FIRST (i,a) factor.** The first ucp impl conjugated
  the second factor (building `conj(C_Φ)=C_Φᵀ`); it passed unital/CP/round-trip
  (all transpose-invariant) and was caught only by the oracle test asserting a Choi
  block equals `aic_ucp_apply(Φ, E_ij)`. That oracle + a complex-Kraus fixture now
  guard it. Unital check via `aic_mat_partial_trace_left` (traces out the LEFT/K
  factor → keeps H).
- **Carrier `M = range(Q)`, `Q = Σ_a K_a K_a† = Φ*(1)`** (lem_carrier tex:1724).
  `M ⊊ H ⇔ Φ is NOT trace-preserving` (i.e. `Q≠I`). Do NOT assume `M=H` (a web
  source wrongly claimed "carrier=H for UCP" — it confused `Φ(I)=I` with the
  carrier). The certified rank of `Q` is gap-dependent → `aic-w4o.1`.
- **`idemp` (th_idemp_structure, the η=0 oracle).** `A = ImgΦ` is built by the
  **column-image SVD** (apply `Φ` to all `E_ij`, SVD the stacked outputs) — NOT the
  `n²×n²` superoperator matrix (that route risks a vec-convention bug; the
  column-image route reuses the tested `aic_ucp_apply`). `Γ` is built via the
  **`Λ=ΔΓ` factorization** (`Λ(Y)=Φ(J_M Y J_M†)`, `Γ=Δ†Λ` since `Δ` has orthonormal
  columns), NOT a `w`-pseudoinverse — this `Γ` is UCP by construction (tex:2088).
  `A`'s basis has a **unitary gauge freedom** → cross-check `A` as a SUBSPACE
  (`Π_A=ΔΔ†`), never operator-by-operator.
- **cb-norm / η-defect (Route B, no Python).** `η = ‖Φ²−Φ‖_cb`; `Λ=Φ²−Φ` is
  Hermiticity-preserving but **NOT CP**, so the CP closed form `‖Φ‖_cb=‖Φ(1)‖`
  does NOT apply — a **Watrous diamond-norm SDP** is required (tex:352). Split: the
  **C core** computes `Choi(Φ²−Φ)` (`aic_ucp_compose` gives `Φ²` Kraus `{K_bK_a}`;
  `aic_ucp_choi_diff`); the **SDP lives in Julia+MOSEK** (`tools/gen_fixtures_d24.jl`).
  **Load-bearing normalization: `‖Λ‖_⋄ = (2/n)·SDP_optval`** because Convention-A
  Choi has trace `n` (the Watrous SDP is trace-1 calibrated) — verified across
  n=2,3,4 (the dimension canary pins it). Analytic anchors for tests:
  `‖Dep_d−id_d‖_⋄ = 2(1−1/d²)`; `Φ_t=(1−t)id+t·Dep ⇒ ‖Φ_t²−Φ_t‖_⋄ = t(1−t)·2(1−1/d²)`;
  the paper's own example (tex:367) `‖Φ²−Φ‖_cb = η√(1−η)`.
- **The certified-arb degenerate-eig wall (`aic-w4o.1`) is the recurring deferral.**
  FLINT's certified eig (`acb_mat_eig_simple`) needs a SIMPLE spectrum; the project's
  spectra are degenerate (projections 0/1, `⊕B(L_j)`). So the *certified* halves of
  Kraus-extraction (`ucp`), `M`/`A`-subspace extraction (`idemp`), and `projection`
  all use the LAPACK double path now and defer certification to `aic-w4o.1`
  (tooling: `acb_mat_eig_multiple_rump`, `acb_mat_eig_global_enclosure`, or an
  eig-free Cholesky route). The *checks* that avoid full eig (PSD via
  `herm_max_eig(−C)`, relation defects) are already certified.

## What's next (ready work — `bd ready`)

Issue tracker is **beads**, prefix `aic` (persistent across sessions; JSONL at
`.beads/issues.jsonl` is committed). `bd ready` for the live list. **This session
BUILT `aic-knm` (`ecstar`) Cycles 1+2** (data model + basis-sweep estimators +
HOPM worst-case search, committed) and filed `aic-b7c`/`aic-3qq`/`aic-cr6`/
`aic-0at`/`aic-4c7`. Prior session CLOSED `aic-m24`/`aic-cne`; earlier
`aic-c7n`/`aic-wuh`/`aic-9kk`/`aic-d24`.

- **`aic-b7c` (P1) — RUN FIRST.** The deferred Cycle-2 hostile review (interrupted
  by the outage); its bead description is the full review brief. Blocks `aic-knm`
  closure. Fix any findings, then `bd close aic-knm`.
- **`aic-knm` (P2, in_progress)** — `ecstar`: Cycles 1+2 committed & green; remaining
  to close = `aic-b7c` (review) passes. Cycle 3 (certified SDP UPPER bound) is the
  separate `aic-0at`. See the playbook below.
- **`aic-92f`** (`assoc_ecsa`, almost-idempotent ε-C* via Φ̃=θ(2Φ−1)) — **the next
  headline step**, unblocked by `aic-knm`'s deliverables. Builds A=ImgΦ̃ (θ on the
  superoperator 2Φ−1, reuse `funcalc`), feeds the `ecstar` estimators, and finally
  exercises the Π_A teeth (`aic-3qq`) + the real O(η) canary (`aic-4c7`). η
  measurable/certifiable now (Julia+MOSEK + the arb certifier).
- **`aic-ssu` (P2, NEW)** — Julia end-to-end `certify(kraus)→(lo,hi)` wrapper:
  solve MAX primal + MIN dual (already in `src/sdp.jl`), ccall `aic_cbnorm_certify_d`,
  return the rigorous bracket; test vs the wolframscript oracle on live channels.
- **`aic-w4o.1` (P1)** — certified degenerate Hermitian eig (arb). The main
  certified-path debt: gates the certified extraction deferred in `ucp` (Kraus),
  `idemp` (M, A subspaces), and `projection`. Tooling: `acb_mat_eig_multiple_rump`
  (cluster enclosures), `acb_mat_eig_global_enclosure`; audition vs an eig-free
  Cholesky route (bead notes).
- **`aic-m24` (P2, NEW)** — certified-arb cb-ball (KKT/rational-certificate; no arb
  SDP solver — research escalation) + reusable diamond-norm value entry point (E5).
- Rest of **E2-eps** (`aic-vs9`): `unitfix`, `unitary`, `projection` (blocked by
  `aic-w4o.1`), `corner`. Then **E3-mainthm**, **E4-headline**, **E5-julia**,
  **E6-research** (`aic-1bc` c_0, `aic-1sk` factorization closure).
- **`aic-dbo.2` (P1)** — adversarial instance generators (evil-matrix corpus; unify
  with `aic-f9u.1`).
- **NEW follow-ups from this session:** `aic-ynu` (P2, Artin–Wedderburn block
  decomposition + prop_Gamma explicit form + Δ/Γ Kraus reps); `aic-kyj` (P3, idemp
  test teeth: single-split block-diag channel + genuine double-vs-arb subspace
  oracle when `aic-w4o.1` lands).
- Housekeeping: `aic-w4o.4` (split oversized test files — `test_mat/funcalc/
  contraction/latd/idemp` and `test_idemp.h` (208 LOC) all > 200 LOC); `aic-w4o.2`
  (full arb SVD); `aic-w4o.3` (opnorm power-iteration audition); `aic-f9u.1`
  (shared corpus); `aic-f9u.2` (multi-dim Pareto bench reporting).

## Next-step playbook — `aic-knm` (`ecstar`), the chosen next module

Goal: the **ε-$C^*$ algebra data model + axiom-defect estimators** (shard A; defs
at `tex:403–492`, the approximate axioms at `tex:407–439`). This is the last
blocker for `aic-92f` (`assoc_ecsa`), so it is the direct route to the headline.

- **Research first** (the discipline): one paper-leg (Sonnet) on shard A + `tex:403–492`
  — the exact ε-$C^*$ axioms (associativity, `‖XY‖≤(1+ε)‖X‖‖Y‖`, the $C^*$ identity,
  unit laws — ALL only up to ε), the δ-homomorphism defs (`tex:443–455`; keep ε vs δ
  distinct), and what an "axiom-defect estimator" must compute; plus one web-leg
  (best-in-class: estimating the associativity/multiplicativity defect of a finite-dim
  near-algebra — operator norm of the associator trilinear map; how to take the sup
  over the unit ball — basis sweep vs power iteration). Decide the route, then implement.
- **Data model.** An ε-$C^*$ algebra here is a subspace `A ⊆ M_N` with the
  **Choi–Effros star** `X⋆Y = Φ(XY)` (NOT plain `XY`; the product leaves `ImgΦ`),
  involution `X↦X†`, norm inherited from `M_N`. For the exact case `A = ImgΦ`
  (reuse `aic_idemp`); for the almost-idempotent case (built by `assoc_ecsa`)
  `A = ImgΦ̃`, `Φ̃ = θ(2Φ−1)` (reuse `funcalc` θ on the **superoperator** `2Φ−1`).
  Store `A` as its orthonormal operator basis (the `aic_idemp` `A_basis` shape).
- **Defect estimators** (return certified arb balls): associativity
  `‖(X⋆Y)⋆Z − X⋆(Y⋆Z)‖` over the unit ball; submultiplicativity slack
  `sup ‖X⋆Y‖−(1+ε)‖X‖‖Y‖`; $C^*$-identity `‖X⋆X†‖` vs `‖X‖²`; unit-law defect
  (the exact unit only exists after `prop_unit`, `tex:672` — `unitfix`, a sibling
  bead). Each star multiplication goes through `aic_ucp_apply`.
- **η=0 oracle** (mandatory): for an EXACT idempotent Φ, `A=ImgΦ` is a genuine
  $C^*$ algebra → all defects must be **machine-zero**. Reuse the six `test_idemp`
  channels. Then an ε-sweep (perturb toward almost-idempotent) must show every
  defect growing `= O(ε)`, and the **dimension sweep** must show the constant does
  NOT grow with `dim A` (the universality canary, `aic-dbo.3`).
- **Reuse** `aic_ucp` (the star), `aic_mat` (`opnorm`/`herm_max_eig`), `aic_idemp`
  (the exact-case `A`). η is measurable via the `aic-d24` Julia+MOSEK golden master.
- `bd show aic-knm` for the bead; then `bd update aic-knm --claim`.

## Orchestration mode (how this was run, to continue it)

The user runs this as a serial orchestration: the orchestrator stays in the loop,
delegates each step to a **subagent** (Opus for coding/hostile-review, Sonnet for
research/summarization), and **independently verifies** (don't trust a subagent
report — re-run `make clean && make test`, read the diff, check zero warnings).
The per-module pipeline used this session, which works:

1. **Research (parallel OK)** — for every core module, a paper-constructivization
   leg (shard + `.tex`) **and** a web-leg surveying best-in-class published
   algorithms/impls. *Web-research EVERY algorithm* (user directive) — cb-norm,
   fixed-point structure, etc. **Research may run in parallel; Julia may NOT.**
2. **Decide the route** (orchestrator), recording the audition slate (Law 4).
3. **Implement** (Opus subagent), tests-first TDD, both number paths, ≤200 LOC/file.
4. **Independent build-verify** (orchestrator).
5. **Hostile review** (Opus subagent) — ALWAYS for a core module; it has caught a
   real defect every single time (see the "tests that can't fail" note above).
6. **Fix pass** (Opus), mutation-proving the new teeth.
7. **Commit atomically** per closed bead with `.tex`-line provenance; file beads
   for every deferral as it arises.

Standing directives (also in the harness memory dir): **No Python in the codebase**
— Julia or C only (consumer Python bindings come later). **Golden masters** are
encouraged as an independent cross-check rung — MOSEK/Gurobi (local) for SDPs,
Mathematica/`wolframscript` (needs the user's VPN) for exact fun-calc, driven from
Julia or C. **No parallel Julia** (precompile/env fragility — one Julia process at
a time). `TaskCreate` is fine for ephemeral in-session tracking; **beads is the
only persistent tracker** (`bd ready`/`bd show`; export to `.beads/issues.jsonl`
before committing — the pre-commit hook also syncs it). **No GitHub CI; no remote
configured** (commit locally; do NOT push — the auto-block's push step is N/A).

## Conventions / gotchas

- Cite `approximate_algebras.tex:<line>` + the verbatim equation in every math
  routine (Law 1). `≤200 LOC` per **source** file (several test files exceed it →
  `aic-w4o.4`). Arb-path `.c` files never *call* LAPACK (the link line bundles both
  for tests; the rule is about which library each routine calls). Fail-loud on
  precondition violation / LAPACK `info!=0` / lost precision / an arb ball that
  straddles a decision boundary (e.g. PSD: `herm_max_eig(−C)` ball straddling 0).
- `acb_mat_eig_simple` needs a SIMPLE spectrum; `LAPACKE_zheev` does not. Vec/Choi
  is **row-major** (`vec(X)[a*n+b]=X[a,b]`; Choi K-major); `latd` conversions are
  row-major. The Choi conjugation is on the FIRST factor (see Channel-module
  conventions — this is the bug that bit).
- `tests/fixtures_d24.inc.h` is **generated but committed** (so `make test` needs
  no Julia); regenerate with `make fixtures` after touching the SDP/corpus, and the
  Makefile makes `tests/*.inc.h` a test prerequisite so a stale fixture forces a
  rebuild. Julia work is **serial** (no parallel Julia). **No Python.**
- Before committing, `bd export > .beads/issues.jsonl` so the bead state lands with
  the commit (the pre-commit hook in `core.hooksPath=.beads/hooks` also syncs it).
  Close beads with a detailed `--reason`; force-close (`--force`) only when a
  deliverable ships but a DAG dep is a deferred remainder (tracked elsewhere) — as
  `aic-wuh` was (its certified-extraction half is on `aic-w4o.1`).
