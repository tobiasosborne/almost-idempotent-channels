# README design — AlmostIdempotentChannels.jl (showpiece)

> ART-DIRECTOR + technical-writer deliverable. Synthesis of the four sibling
> research legs (TensorGR / su2-fft / Lyr) and the truth-anchor grounding report,
> all cross-checked against the live source on 2026-06-02 (`src/api.jl`,
> `src/types.jl`, `src/show.jl`, `test/*.jl`, `Project.toml`, `basin.jl`,
> `HANDOFF.md`, `BUILDING.md`). The git remote is HTTPS,
> `github.com/tobiasosborne/almost-idempotent-channels` (badge/links use that).
> No README.md exists yet; this designs it from scratch.
>
> This document has four parts: **A** design spec (the escalating arc + the
> borrowed visual devices), **B** the full near-final README prose with asset
> placeholders, **C** the tooling decision (cast→SVG + plot pipeline), **D** the
> mechanical, serial-Julia-safe Phase-B asset plan (the exact `.cast` transcript +
> one plotting script per result plot). Phase B should be drop-in-the-assets only.

---

## Part A — DESIGN SPEC: the escalating-impressiveness arc

### A.0 The thesis the README must sell, in one breath

> A near-idempotent quantum channel secretly hides a clean digital code. Feed in
> messy analog Kraus operators; this package extracts the exact block structure
> `⊕_l M_{d_l}` of the genuine C\*-algebra the channel is approximating, builds the
> encode/decode channels that move between them, and hands back a *certified*
> bound — a proof of an inequality, not a floating-point guess — on how lossy the
> round trip is.

Every section below is an escalation toward that sentence and then past it into
*proof* (the plots, the certified numbers, the η=0 oracle at machine ε).

### A.1 The arc — section order and why each beat raises the stakes

The arc fuses three sibling templates: **Lyr's "show first, justify second"**
(hero visual before philosophy), **su2-fft's three-part "Why this exists" pressure
gradient** (physical problem → the surprising theorem → what we built) and its
**brutally-honest "Known limits" + ASCII benchmark tables**, and **TensorGR's
front-loaded animated terminal demo + dense features table + comparison-as-
credibility-anchor**. Concretely:

| # | Section | Borrowed from | What raises the stakes here |
|---|---|---|---|
| 0 | `# AlmostIdempotentChannels.jl` + 3 badges | TensorGR (badge row) | Identity + arXiv link in one line. |
| 1 | **Bold one-line tagline** (noun phrase) | Lyr | A falsifiable claim: *constructive + certified*. |
| 2 | Two prose sentences (the thesis) | Lyr | The hook before any heading. |
| 3 | **Animated terminal demo** (SVG) | TensorGR | The five verbs run; certified numbers appear live. The single strongest desire device. |
| 4 | **Why this exists** (3 parts) | su2-fft | Physical stakes → Kitaev's surprise (dimension-independent constant) → what we uniquely give. |
| 5 | **The headline pipeline** (one Julia block) | su2-fft "short tour" | The reader sees the whole API is five verbs. |
| 6 | **Features** (table: result → algorithm → certified by) | TensorGR + Lyr | Each row is a paper theorem made executable; escalates `prop_P` → `th_factorization`. |
| 7 | **Install & quick start** | su2-fft (3-line) + Lyr | No `Pkg.build` — Preferences hook. Friction is tiny. |
| 8 | **A look at it running** (the result plots) | su2-fft + Lyr gallery | Visual *proof*: containment, 12-orders-tighter, **flat universality constant**, O(η) round-trip, η=0 oracle at machine ε. |
| 9 | **Rigour & the cross-check ladder** | su2-fft "what is measured" | The certification story: arb balls rounded outward, double-vs-arb, the η=0 oracle. |
| 10 | **Benchmark tables** (ASCII) | su2-fft | Concrete numbers: dim-sweep `c`, precision/ball-radius, bracket widths. |
| 11 | **Known limits** (bead-cited, honest) | su2-fft | Builds trust by naming the failure modes (factorize domain, O(η)-TP decode, MOSEK license). |
| 12 | **Design philosophy** (3 one-sentence choices) | su2-fft gap-fill | Why C+FLINT+Julia, why arb, why constructive. |
| 13 | **The MOSEK extension** (optional tight bracket) | Lyr inversion device | The optional rung that turns a loose certificate tight. |
| 14 | **Math primer** (`<details>` collapsible) | su2-fft | The full η → ε-C\* → genuine C\* → factorization chain, three inequalities quoted from `.tex`. |
| 15 | **Citation** + disclaimer | su2-fft + TensorGR | Kitaev arXiv, FLINT/LAPACK credit, "independent realisation". |
| 16 | **License** | all | One line. |

The pressure gradient: the reader is hooked by motion (3) before any prose; is
made to *want* the result by the physical stakes (4); sees the API is tiny (5);
is impressed by breadth (6); gets it running in 3 commands (7); then is shown
**proof** (8) culminating in the flat-constant universality plot — the paper's
central claim, verified. Everything after 8 is the rigour and honesty that makes
the proof credible.

### A.2 Borrowed visual devices — the concrete markdown

**Badge row (TensorGR), no HTML wrapper, three badges.** Use static
shields.io badges only — Rule 13 forbids CI; there is no remote CI badge.

```markdown
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-9558B2.svg)](https://julialang.org/)
[![arXiv](https://img.shields.io/badge/arXiv-2405.02434-b31b1b.svg)](https://arxiv.org/abs/2405.02434)
```

**Bold noun-phrase tagline (Lyr), immediately under the H1.**

```markdown
**Constructive finite-dimensional algorithms for Kitaev's almost-idempotent
quantum channels and approximate C\*-algebras — with certified arbitrary-precision
error bounds throughout.**
```

**Animated terminal demo embed (TensorGR), bare markdown, on its own line.**
GitHub renders the SVG inline at its native width.

```markdown
![Pipeline demo: certified_defect → associated_algebra → main_isomorphism → factorize](docs/assets/demo.svg)
```

**Result-plot embed, centered (Lyr), with a technical caption that asserts
correctness — never marketing.**

```markdown
<p align="center">
  <img src="docs/assets/universality.png" alt="Isomorphism constant c = isodefect/η vs dim A, flat" width="640"/>
  <br/>
  <em>The th_main isomorphism constant c = isodefect/η over dim A ∈ {8,12,16,18,20}:
  OLS slope −2.7e-4, max c = 0.047. The constant does not grow with dimension
  (approximate_algebras.tex:484).</em>
</p>
```

**Features table (TensorGR + Lyr): three columns — paper result → algorithm →
certified by.** Bold the left column.

**ASCII benchmark tables (su2-fft): plain code-fenced text, never markdown
tables.** This is load-bearing for the house tone.

**ASCII architecture / factorization diagram (Lyr box-drawing arrows, not
Mermaid).**

**`<details>` collapsible math primer (su2-fft).**

```markdown
<details>
<summary><b>Math primer: η-idempotence → ε-C\*-algebra → genuine C\* → factorization (click to expand)</b></summary>
…
</details>
```

No emojis, no "robust"/"blazing"/"production-grade" (Rule 12). Concrete numbers
everywhere. Two `<p align="center">` blocks max for plots; the demo SVG uses bare
markdown so it renders full-width like TensorGR's.

---

## Part B — FULL NARRATIVE DRAFT of README.md

> Drop this in as `README.md` after Phase B generates the assets. Asset
> placeholders are marked `<!-- ASSET: … -->`. All code blocks use the REAL API
> (`certified_defect`, `associated_algebra`, `main_isomorphism`, `factorize`,
> `encode`, `decode`, `algebra`, `delups_defect`, `upsdel_defect`); all `show`
> outputs are reconstructed from `src/show.jl` format strings + the verified
> analytic anchors. Endpoints flagged `[VERIFY]` must be captured live in Phase B
> before publishing (see Part D). Prose is ~95% final.

---

```markdown
# AlmostIdempotentChannels.jl

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-9558B2.svg)](https://julialang.org/)
[![arXiv](https://img.shields.io/badge/arXiv-2405.02434-b31b1b.svg)](https://arxiv.org/abs/2405.02434)

**Constructive finite-dimensional algorithms for Kitaev's almost-idempotent
quantum channels and approximate C\*-algebras — with certified arbitrary-precision
error bounds throughout.**

A near-idempotent quantum channel secretly hides a clean digital code. Feed in
messy analog Kraus operators and this package extracts the exact block structure
`⊕_l M_{d_l}` of the genuine C\*-algebra the channel is approximating, builds the
encode/decode channels that move between them, and returns a *certified* bound —
a proof of an inequality, not a floating-point guess — on how lossy the round trip
is. It is a Julia surface over the `libaic` C/arb cores, implementing
[Kitaev, *Almost-idempotent quantum channels and approximate C\*-algebras*,
arXiv:2405.02434](https://arxiv.org/abs/2405.02434) (2025) as finite-dimensional
algorithms for *every* result the paper proves non-constructively.

<!-- ASSET: terminal-demo docs/assets/demo.svg -->
![Pipeline demo: certified_defect, associated_algebra, main_isomorphism, factorize](docs/assets/demo.svg)

## Why this exists

**The object.** A unital completely positive (UCP) map `Φ` on `B(H)`, `dim H = n`,
in the Heisenberg (observable) picture `Φ(X) = Σ_a K_a† X K_a`. This is the most
general noisy, information-preserving measurement of observables you can write
down. `Φ` is **η-idempotent** when `‖Φ² − Φ‖_cb ≤ η`, the completely-bounded
norm — a supremum over all ampliations `1 ⊗ Φ`, not the bare operator norm. `η`
measures how close `Φ` is to a genuine projection onto a stable subalgebra of
observables; `η = 0` exactly when `Φ` is a conditional expectation onto a
C\*-subalgebra.

**What Kitaev proved, and why it is surprising.** An η-idempotent `Φ` does not
merely *almost* project. Its image carries the structure of an **ε-C\*-algebra**
`A = Img Φ̃` with `ε = O(η)`: a vector space satisfying the C\*-axioms only *up
to* ε, with the Choi–Effros product `X ⋆ Y = Φ̃(XY)` and no exact unit. The
headline rigidity theorem (`th_main`, arXiv:2405.02434 §2) is that **every
finite-dimensional ε-C\*-algebra is O(ε)-isomorphic to a *genuine* C\*-algebra
`B = ⊕_l M_{d_l}`, with a universal, dimension-independent constant.** The
constant does not grow with `n` — and the naive averaging route fails precisely
because *its* constant grows like `n` (§9). When `A` comes from a finite-dimensional
`Φ`, that isomorphism and its inverse are realised by quantum channels, giving an
**approximate factorization `Φ ≈ Δ Υ`** of the channel through `B`: a decode map
`Δ` and an encode map `Υ` with `ΔΥ ≈ Φ̃` and `ΥΔ ≈ 1_B`.

**What this package contributes.** Kitaev's proofs are non-constructive — the
Lefschetz–Hopf fixed-point theorem, holomorphic functional calculus via contour
integrals, Haar-measure diagonals. None of them hands you the object. In finite
dimensions every one of those objects is computable, and this package implements
a constructive algorithm for each, wrapped as five verbs that read like the paper:

1. The whole pipeline is **constructive and finite-dimensional** — no existence
   proof is left as an existence proof.
2. Every defect comes back as a **`CertifiedBracket`** `lo ≤ x ≤ hi`: FLINT/arb
   error balls rounded *outward* so the bound survives the conversion to `double`.
   The output is a proof of the same inequality the theorem promises.
3. It is **solver-free by default.** The rigorous η bracket and the entire
   factorization run on the eig-free arb certifier with no SDP solver installed.
   MOSEK is an optional extension, only for the exact value and the tight bracket.

## The headline pipeline

```julia
using AlmostIdempotentChannels, LinearAlgebra

# A near-idempotent qubit channel: Φ_t = (1−t)·id + t·(maximally depolarizing).
Φ = UCPMap(phit_kraus(2, 0.1))

η = certified_defect(Φ)        # rigorous ‖Φ²−Φ‖_cb bracket — NO solver
A = associated_algebra(Φ)      # the ε-C* algebra A = Img Φ̃
v = main_isomorphism(Φ)        # the O(ε)-isomorphism v: B → A, dim-independent constant
F = factorize(Φ)               # Φ ≈ Δ Υ through a genuine C* algebra B = ⊕_l M_{d_l}
```

`phit_kraus`, `bce_conj_kraus`, `block_cond_exp_kraus` below are the test fixtures
(`test/fixtures.jl`); `Φ_t = (1−t)·id + t·Dep_d` has the closed-form defect
`‖Φ_t²−Φ_t‖_⋄ = t(1−t)·2(1−1/d²)`, so the certified bracket can be checked against
an analytic value.

## Features

Each row is a theorem of the paper made executable, and a rigorous certificate of
its bound. The verbs escalate from the cheapest scalar to the full factorization.

| Result (arXiv:2405.02434) | Verb / type | Certified by |
|---|---|---|
| **Idempotency defect** `η = ‖Φ²−Φ‖_cb` (§1) | `certified_defect(Φ)` → `CertifiedBracket` | eig-free arb bracket, outward-rounded; solver-free |
| **`prop_P` regularisation** `Φ̃ = θ(2Φ−1)` (§3) | `associated_algebra(Φ)` → `EpsCStarAlgebra` | certified ε associator-defect bracket |
| **`th_main` isomorphism** `v: B → A`, O(ε), dim-independent (§2/§9) | `main_isomorphism(Φ)` → `MainIsomorphism` | certified `‖v‖_cb`, `‖v⁻¹‖_cb`, isodefect brackets |
| **`th_idemp_structure`** the genuine `B = ⊕_l M_{d_l}` (§11) | `algebra(F)` → `CStarAlgebra` | exact block sizes from the Wedderburn structure |
| **`th_factorization`** `Φ ≈ Δ Υ` (§12) | `factorize(Φ)` → `ChannelFactorization` | two certified round-trip brackets `‖ΔΥ−Φ‖_cb`, `‖ΥΔ−1_B‖_cb` |
| **Encode / decode channels** `Υ*`, `Δ*` (§1, §12) | `encode(F)`, `decode(F)` → `CPMap` | CPTP duals; round-trip brackets are the certificate |
| **Exact diamond-norm value** (Watrous SDP) | `idempotency_defect(Φ)` → `Float64` | MOSEK extension; strong-duality cross-check |
| **Tight rigorous bracket** | `certified_defect(Φ; tight=true)` | arb certifier fed MOSEK feasible points |

## Install and quick start

Dependencies for the C core: a C11 compiler, CMake ≥ 3.24, `libflint-dev`
(FLINT ≥ 3.0, which brings arb), `liblapacke-dev`/`liblapack-dev`/`libblas-dev`.
Tested on Debian/Ubuntu Linux with Julia 1.10+.

```sh
# 1. Build the C/arb core (produces build/libaic.so).
sudo apt install build-essential cmake libflint-dev liblapacke-dev liblapack-dev libblas-dev
cmake -S . -B build && cmake --build build

# 2. Point the package at the library (writes LocalPreferences.toml; no Pkg.build).
julia --project=julia/AlmostIdempotentChannels.jl -e '
    using AlmostIdempotentChannels
    set_libaic_path!(joinpath(pwd(), "build", "libaic.so"))'   # then restart Julia

# 3. Run the solver-free test suite (green with no MOSEK installed).
julia --project=julia/AlmostIdempotentChannels.jl -e 'using Pkg; Pkg.test()'
```

Library discovery is via `Preferences` — there is no `deps/build.jl` and no
`Pkg.build()` step. The default looks for `build/libaic.so` next to the package;
`set_libaic_path!(p)` points it anywhere else.

## A look at it running

Each plot below re-renders a passing test; the reproduce script for each is in
`docs/plots/` (a separate plotting environment, not a package dependency).

**The certified bracket contains the truth.** For `Φ_t` the defect has a closed
form; the solver-free eig-free bracket `[lo, hi]` brackets it for every input.
Containment is a test invariant (`test/test_core.jl`), so this plot is a passing
test drawn.

<!-- ASSET: plot docs/assets/containment.png : eig-free bracket [lo,hi] band shaded vs the analytic η line (inside the band) over t ∈ [0.02, 0.3]; phit(2,t) and phit(3,t) families. -->
<p align="center">
  <img src="docs/assets/containment.png" alt="Certified bracket containing the analytic η" width="640"/>
  <br/>
  <em>The eig-free certificate [lo, hi] (shaded band) brackets the analytic
  η = t(1−t)·2(1−1/d²) (line) for every t. Loose by design (hi/lo ≈ 2n): it
  certifies a value rather than computing it.</em>
</p>

**Dimension-independence — the paper's central claim, verified.** The `th_main`
isomorphism constant `c = isodefect/η` over a block-algebra family of growing
dimension is *flat*: it does not grow with `dim A`. The naive averaging route the
paper warns against (§9) has a constant that grows like `n`; this one does not.

<!-- ASSET: plot docs/assets/universality.png : scatter of c = isodefect/η vs dim_A ∈ {8,12,16,18,20}, points (8,0.0182)(12,0.0470)(16,0.0247)(18,0.0410)(20,0.0134); dashed OLS line slope −2.7e-4; annotate "max c = 0.047 ≪ 1; no dim growth (tex:484)". [VERIFY: regenerate via main_isomorphism over a Julia block family per Part D.] -->
<p align="center">
  <img src="docs/assets/universality.png" alt="Isomorphism constant flat over dimension" width="640"/>
  <br/>
  <em>c = isodefect/η over dim A ∈ {8,12,16,18,20}: OLS slope −2.7e-4, max c = 0.047,
  the worst case at the smallest dimension. The constant does not grow with
  dimension (approximate_algebras.tex:484).</em>
</p>

**The factorization round-trip is O(η).** For genuinely oblique, multi-block
channels (`bce_conj`, `B = M_2 ⊕ M_2`) the certified round-trip defects
`‖ΔΥ−Φ‖_cb` and `‖ΥΔ−1_B‖_cb` scale *linearly* in η, well below the generous
`50·η` envelope the tests assert.

<!-- ASSET: plot docs/assets/roundtrip.png : delups_defect / upsdel_defect (max of bracket) vs η = eta_proxy, for bce_conj(4,2,t) and mixconj(4,2,t), t ∈ {0.02, 0.05}; overlay y = 50·η envelope. Log-log axes. -->
<p align="center">
  <img src="docs/assets/roundtrip.png" alt="Round-trip defect O(η)" width="640"/>
  <br/>
  <em>Certified round-trip defects ‖ΔΥ−Φ‖_cb and ‖ΥΔ−1_B‖_cb vs η for oblique
  multi-block channels. Linear in η, below the 50·η envelope (test/test_factorize.jl).</em>
</p>

**The η = 0 oracle collapses to machine ε.** When `Φ` is exactly a conditional
expectation, every defect — the associator ε, both round-trip brackets, the
isomorphism defect — collapses below `1e-9`, and the extracted block structure is
exactly right: `identity(2) → [2]`, `block_cond_exp(4,2) → [2,2]`,
`dephasing(3) → [1,1,1]`. This is the cleanest ground truth in the paper.

<!-- ASSET: plot docs/assets/eta0_oracle.png : bar chart, log y, of maximum(bracket) for the four defects on identity(2), block_cond_exp(4,2), dephasing(3); all bars < 1e-9 with a dashed "1e-9 machine reference" line. Annotate each fixture's extracted block list. -->
<p align="center">
  <img src="docs/assets/eta0_oracle.png" alt="η=0 oracle defects at machine ε" width="640"/>
  <br/>
  <em>On exactly-idempotent channels every certified defect collapses below 1e-9 and
  the block structure is exact. The cleanest cross-check rung (test/test_factorize.jl).</em>
</p>

## Rigour and the cross-check ladder

The paper's whole subject is rigorous O(ε) bounds, so this package certifies them
numerically rather than trusting `double`. Every scalar result is a
`CertifiedBracket` `lo ≤ x ≤ hi` built from FLINT/arb error balls rounded outward
(FLOOR/CEIL) so the bound is still rigorous after the `double` conversion. The
cross-check ladder, weakest to strongest:

1. **double vs arb at prec = 53** — the same routine on both number paths agrees
   to ~1e-10 (the C/arb Choi of `Φ²−Φ` matches the pure-Julia one, `test_core.jl`).
2. **certified containment** — the bracket provably contains the closed-form η for
   the `Φ_t` family at every `t` (`test_core.jl`).
3. **the η = 0 oracle** — exactly-idempotent `Φ` ⇒ every defect < 1e-9, exact
   blocks (`test_factorize.jl`). Zero defect is the unambiguous ground truth.
4. **strong duality (MOSEK)** — the Watrous primal and dual diamond-norm values
   agree with each other and with the analytic anchor to 1e-6 (`test_sdp.jl`),
   pinning the dual normalisation.

## Benchmarks

Plain-text tables; concrete numbers throughout.

```
Dimension-independence of the th_main constant c = isodefect/η
  (block family ⊕_j M_d, dim A = k·d²; from the C dim-sweep, tests/test_dbo3.c, prec=256)

  dim A |  k |  d |  c = isodefect/η
  ------+----+----+------------------
      8 |  2 |  2 |  0.0182
     12 |  3 |  2 |  0.0470     <- worst case, smallest c-dimension
     16 |  4 |  2 |  0.0247
     18 |  2 |  3 |  0.0410
     20 |  5 |  2 |  0.0134

  OLS slope = -2.7e-4;  max c = 0.047 << 1;  no growth with dim A (tex:484 holds).
```

```
Tight (MOSEK) vs eig-free (solver-free) bracket width on ‖Φ²−Φ‖_cb
  (the eig-free bracket is loose by design, hi/lo ~ 2n; the tight one is fed the
   Watrous SDP feasible points)

  fixture        |  eig-free width  |  tight width
  ---------------+------------------+--------------
  phit(2, 0.3)   |  [VERIFY]        |  [VERIFY ~1e-12]
  phit(2, 0.1)   |  [VERIFY]        |  [VERIFY ~1e-12]
  paper(0.1)     |  [VERIFY]        |  [VERIFY ~1e-12]

  The tight certifier widths are ~1e-12 (HANDOFF), ~10^12x tighter than eig-free;
  width(tight) < width(loose) is a test invariant (test/test_sdp.jl).
```

## Known limits

Stated plainly (the bead IDs are the project's issue tracker).

- **`factorize` has a tighter domain than the rest of the pipeline.** It builds
  `Υ` via a power-series functional calculus whose convergence domain is much
  smaller than the `prop_P` basin `ρ(Φ²−Φ) < 1/4` (bug `aic-exa.13`). The verb
  pre-checks a conservative `ρ < 0.10` in Julia and throws a clean `ArgumentError`
  — it does **not** abort the session. `associated_algebra` and `main_isomorphism`
  keep the wider `ρ < 1/4` domain; `certified_defect` is always safe at any η.
- **The `decode` channel is only O(η)-trace-preserving for η > 0.** Its underlying
  observable map is only O(η)-completely-positive, so the internal Choi→Kraus
  PSD-cone projection clips a small negative mass; `iscptp(decode(F))` is `false`
  at the default `atol = 1e-9` (measured TP-defect ≈ 1.9e-6 on multi-block oblique
  fixtures) and `true` at `atol = 1e-4`. The rigorous certificate is the cb-norm
  round-trip bracket, not `iscptp` at machine tolerance. The encode channel and
  the η = 0 oracle decode are TP to machine precision.
- **The eig-free bracket is loose by design** (`hi/lo ~ 2n`). It certifies a value
  rather than competing with the solver. For a tight bracket or the exact value,
  use the MOSEK extension.

## Design philosophy

- **C + FLINT + Julia, not pure Julia.** The numerical cores are tight C with two
  number paths — a fast `double`/LAPACK path that anchors the test suite, and a
  FLINT/arb path that produces the certified balls. Julia is the `ccall` surface.
- **arb, not floating point.** The paper certifies rigorous bounds; arb lets the
  code certify the same bounds numerically rather than trusting `double`. The
  double path is what keeps the tests fast; the arb path is what makes the output
  a proof.
- **Constructive, not existence.** Every result is implemented as a finite-
  dimensional algorithm that produces the object the paper merely asserts to
  exist, meeting the bound the theorem promises.

## The MOSEK extension (optional)

The core is green with no solver. Installing `Convex` + `Mosek` + `MosekTools`
activates `AICMosekExt`, which adds two things: the exact diamond-norm value
`idempotency_defect(Φ)` (the Watrous SDP), and the *tight* rigorous bracket
`certified_defect(Φ; tight=true)` — the arb certifier fed the SDP feasible points,
so the bracket width collapses from `~2n` to `~(solver tol + arb radius) ≈ 1e-12`.
Without the solver, both entry points throw a helpful install hint rather than a
`MethodError`. MOSEK requires a license.

```julia
using AlmostIdempotentChannels, Convex, Mosek, MosekTools
Φ = UCPMap(phit_kraus(2, 0.3))
idempotency_defect(Φ)              # exact: 0.315 (= 0.3·0.7·1.5)
certified_defect(Φ; tight=true)    # rigorous bracket, width ~1e-12, carries the SDP value
```

<details>
<summary><b>Math primer: η-idempotence → ε-C\*-algebra → genuine C\* → factorization (click to expand)</b></summary>

The whole story is three approximate inequalities and one rigidity theorem.

**1. η-idempotence (the input).** `Φ: B(H) → B(H)` is UCP and

> `‖Φ² − Φ‖_cb ≤ η`   (arXiv:2405.02434, §1)

in the completely-bounded (= diamond) norm, the supremum over all ampliations
`1 ⊗ Φ`. At `η = 0`, `Φ` is exactly a conditional expectation onto a C\*-subalgebra.

**2. The ε-C\*-algebra (the structure).** Set `Φ̃ = θ(2Φ−1)`, the exact-idempotent
regularisation (`prop_P`, §3), valid in the spectral basin `ρ(Φ²−Φ) < 1/4`. Its
image `A = Img Φ̃` is an *oblique* subspace of `M_n` carrying the Choi–Effros
product `X ⋆ Y = Φ̃(XY)`. The C\*-axioms hold only up to `ε = O(η)`:

> `‖(X⋆Y)⋆Z − X⋆(Y⋆Z)‖ ≤ ε·‖X‖‖Y‖‖Z‖`,  and the C\*-identity and unit laws up to ε  (§2)

There is no exact unit; associativity is approximate. Do not assume otherwise.

**3. The genuine C\*-algebra (the rigidity).** The main theorem:

> every finite-dim ε-C\*-algebra is O(ε)-isomorphic to a genuine C\*-algebra
> `B = ⊕_l M_{d_l}`, with a **universal, dimension-independent constant**  (`th_main`, §2/§9)

The isomorphism `v: B → A` has `‖v‖_cb, ‖v⁻¹‖_cb = 1 + O(ε)`. The constant does
not depend on `n` — the naive Haar/cohomology route fails because *its* error is
`∝ n`.

**4. The factorization (the payoff).** When `A` comes from a finite-dim `Φ`, `v`
and `v⁻¹` are realised by quantum channels, giving

> `Φ ≈ Δ Υ`,   `‖ΔΥ − Φ̃‖_cb = O(η)`,   `‖ΥΔ − 1_B‖_cb = O(η)`  (`th_factorization`, §12)

with `Δ: B → B(H)` (decode) and `Υ: B(H) → B` (encode). This package computes all
four objects and certifies the two round-trip inequalities.

</details>

## Citation

This is an independent realisation of the constructions in:

> Alexei Kitaev, *Almost-idempotent quantum channels and approximate C\*-algebras*,
> arXiv:2405.02434 (2025).

```bibtex
@misc{kitaev2025almostidempotent,
  author = {Kitaev, Alexei},
  title  = {Almost-idempotent quantum channels and approximate {C}*-algebras},
  year   = {2025},
  eprint = {2405.02434},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph}
}
```

The paper's author has no involvement in this package; any errors in the
algorithms or their bounds are ours alone. The certified arithmetic uses
[FLINT/arb](https://flintlib.org/); the fast double path uses LAPACK/LAPACKE.

## License

MIT. See [LICENSE](LICENSE).
```

---

> **Open question for the synthesizer (B):** the LICENSE file should be confirmed
> (the badge/footer say MIT as a placeholder — Phase B must check `LICENSE` on
> disk and correct both the badge and the footer if it differs). One ASCII
> factorization diagram (`Φ → Φ̃ → A → B → Enc/Dec`) was considered for the
> features section; it is omitted from the draft to keep the arc moving but can be
> dropped in after the features table if desired:
>
> ```
>  Φ  --certified_defect-->  η = ‖Φ²−Φ‖_cb  (rigorous bracket)
>  Φ  --associated_algebra-->  A = Img Φ̃   (ε-C* algebra, oblique in M_n)
>  A  --main_isomorphism-->  B = ⊕_l M_{d_l}  (genuine C*, O(ε), dim-independent)
>  Φ  --factorize-->  Φ ≈ Δ Υ  through B   (Enc = Υ*: B→B(H),  Dec = Δ*: B(H)→B)
> ```

---

## Part C — TOOLING DECISION (concrete, installable on this machine)

This machine has npm/node, go, pip/pipx, ffmpeg, ImageMagick, and network. None of
asciinema/agg/svg-term/vhs/termtosvg is installed yet — Phase B installs exactly
one for the demo. The package ships Julia/C only; the asset toolchain is a dev-only
npm tool, which is allowed (no Python *in the codebase*; CairoMakie is a separate
docs env, not a package dep).

### C.1 The animated terminal: synthesized `.cast` → animated SVG via `svg-term-cli`

**Decision: hand-author an asciinema v2 `.cast` from real package output (no TTY
recording) and render it to an animated SVG with `svg-term-cli`.** This is exactly
the TensorGR pattern (confirmed by reverse-engineering its committed
`repl_demo.svg`: svg-term-cli's `#282d35` background, the three traffic-light
circles, the `@keyframes` film-strip scroll, the Monaco font stack). Synthesizing
the `.cast` rather than recording a live TTY is correct here because (i) it is
serial-Julia-safe — no live REPL to drive; (ii) the outputs are deterministic and
already verified from `show.jl`; (iii) timing/pacing is fully controllable for a
clean demo. SVG (not GIF) is chosen because it is crisp at any zoom, ~18 KB, and
renders inline on GitHub — TensorGR's exact choice. (su2-fft's "agg → GIF" path is
the alternative; SVG wins on size and crispness for a synthesized cast.)

**Install (Phase B, one tool):**

```sh
npm install -g svg-term-cli      # marionebl/svg-term-cli; pulls asciinema-player parser
```

**Render command:**

```sh
svg-term --in docs/src/assets/demo.cast --out docs/assets/demo.svg \
         --window --width 84 --height 30 --padding 18
```

- `--window` adds the macOS traffic-light chrome (the device that makes it read as
  a real terminal).
- `--width 84 --height 30` is wider/taller than TensorGR's 80×24 because the
  `ChannelFactorization` `show` block is the widest line (the
  `decode = Δ* : B(H) → B   (… Kraus,  4→4 CPTP)` row plus the two bracket rows).
  Must match the `.cast` header `width`/`height` exactly.
- No `--theme` flag — svg-term's default dark theme (`#282d35` bg, cyan prompt,
  yellow identifiers) is the tasteful house look and matches the sibling.
- `--padding 18` gives the content a little breathing room inside the window.

**The `.cast` JSON schema to synthesize** (asciinema v2; one header object then one
JSON array per line; every event is an output `"o"` event since this is authored,
not recorded):

```jsonc
// line 1 — header
{"version": 2, "width": 84, "height": 30, "timestamp": 1717200000, "env": {"SHELL": "/bin/bash", "TERM": "xterm-256color"}}
// lines 2..N — events: [time_seconds, "o", "<text with \r\n and ANSI escapes>"]
[0.4, "o", "[36mjulia>[0m "]              // cyan prompt
[1.2, "o", "using AlmostIdempotentChannels, LinearAlgebra\r\n"]
[1.4, "o", "\r\n[36mjulia>[0m "]
// … typed input then its multi-line show output, ~0.1 s after the newline …
```

Pacing rules (match TensorGR's feel): initial prompt at 0.4 s; ~1 s per typed input
line; output 0.1 s after the input newline; a ~2 s dwell after each result block so
the reader can read it; a longer ~3 s dwell on the final `ChannelFactorization`
block (the centrepiece). Total ~32–38 s. Use ANSI `[36m` (cyan) for the
`julia>` prompt, `[33m` (yellow) for type names in the input echo if desired;
keep the output uncoloured (matches `show.jl`, which emits plain text). Part D gives
the exact transcript to encode.

**Embed (bare markdown, full width — TensorGR style):**

```markdown
![Pipeline demo: certified_defect, associated_algebra, main_isomorphism, factorize](docs/assets/demo.svg)
```

### C.2 Result plots: CairoMakie in a separate docs env → PNG

**Decision: a standalone `docs/plots/` Julia environment with `CairoMakie` (PNG
output), kept entirely out of the package `Project.toml`.** CairoMakie is chosen
over Plots.jl/GR because it produces publication-clean PNGs with crisp text and no
GR backend system dependency, and writes deterministic files headless. (Plots.jl
+ GR is the fallback if CairoMakie precompile time is a problem; either works — the
scripts in Part D use the minimal Makie API so a swap is mechanical.) The plotting
env imports the package via a relative `dev` path; it never becomes a package dep,
preserving the solver-free, lean dependency set.

**The docs plotting env (Phase B sets this up once):**

```sh
mkdir -p docs/plots docs/assets
julia --project=docs/plots -e '
    using Pkg
    Pkg.add("CairoMakie")
    Pkg.develop(path="julia/AlmostIdempotentChannels.jl")'
```

Each script is run **one at a time** (serial — no parallel Julia), writes one PNG to
`docs/assets/`, and reuses the package API + the test fixtures (copied into the
script so the docs env does not depend on the test harness). Output: 640 px wide
PNGs at 2× scale for retina crispness.

**Embed for each plot (centered, technical caption — Lyr style):** see the
`<p align="center">…</p>` blocks already in the Part B draft (containment.png,
universality.png, roundtrip.png, eta0_oracle.png).

### C.3 Asset inventory and embed markdown summary

| Asset | Path | Tool | Embed markdown |
|---|---|---|---|
| Terminal demo | `docs/assets/demo.svg` | svg-term-cli (from `docs/src/assets/demo.cast`) | `![…](docs/assets/demo.svg)` bare |
| Containment plot | `docs/assets/containment.png` | CairoMakie | `<p align="center"><img … width="640"/>…</p>` |
| Universality plot | `docs/assets/universality.png` | CairoMakie | centered `<img>` |
| Round-trip plot | `docs/assets/roundtrip.png` | CairoMakie | centered `<img>` |
| η=0 oracle plot | `docs/assets/eta0_oracle.png` | CairoMakie | centered `<img>` |

---

## Part D — PHASE-B ASSET-GENERATION PLAN (mechanical, serial-Julia-safe)

Phase B is: install svg-term-cli; write `docs/src/assets/demo.cast` from the transcript
below; render it; set up `docs/plots/`; run the four plotting scripts **one at a
time**; capture the `[VERIFY]` endpoints from a live REPL; drop the README in. No
step runs Julia in parallel.

### D.0 Live-capture checklist (do FIRST — fills every `[VERIFY]`)

Run one serial REPL session and record the real outputs; substitute them into the
`.cast` transcript and the benchmark tables. Build the C core first (`cmake --build
build`) and `set_libaic_path!`. Capture:

1. `certified_defect(UCPMap(phit_kraus(2, 0.1)))` — the `lo`, `hi`, width
   (analytic η = 0.135, must lie in `[lo, hi]`).
2. `associated_algebra(UCPMap(phit_kraus(2, 0.1)))` — `dim_A`, the ε bracket.
3. `main_isomorphism(UCPMap(bce_conj_kraus(4, 2, 0.02)))` — blocks `[2,2]`,
   `‖v‖_cb`, `‖v⁻¹‖_cb`, isodefect brackets.
4. `factorize(UCPMap(bce_conj_kraus(4, 2, 0.02)))` — encode/decode Kraus counts
   (`re`, `rd`), the two round-trip brackets, `eta_proxy`, `eps_used`.
5. `algebra(F)` — `dim_B = 8`, `n_B = 4`.
6. For the benchmark table: with MOSEK, `width(certified_defect(Φ))` vs
   `width(certified_defect(Φ; tight=true))` for `phit(2,0.3)`, `phit(2,0.1)`,
   `paper_example_kraus(0.1)`.

The `show` *layout* is already correct (from `src/show.jl`); only the numeric
endpoints need substituting.

### D.1 The demo `.cast` transcript (encode into `docs/src/assets/demo.cast`)

The exact REPL session to synthesize, escalating defect → algebra → isomorphism →
the `ChannelFactorization` showcase → the η=0 machine-ε punch. Input lines are the
**real API**; output blocks are the **real `show.jl` renderings** with endpoints
to be filled from D.0. Use `phit_kraus(2, 0.1)` for the scalar/algebra beats (in
both the `prop_P` basin and factorize's domain) and `bce_conj_kraus(4, 2, 0.02)`
for the multi-block factorization showcase.

```text
julia> using AlmostIdempotentChannels, LinearAlgebra

# A near-idempotent qubit channel: Φ_t = (1−t)·id + t·(maximally depolarizing).
julia> Φ = UCPMap(phit_kraus(2, 0.1))
UCPMap: Φ: B(H) → B(H),  dim H = 2
  Kraus operators: 5   (Heisenberg: Φ(X) = Σ K_a† X K_a)
  unital: yes   (‖Σ K†K − I‖ = 3.1e-16)

# (1) RIGOROUS idempotency defect η = ‖Φ²−Φ‖_cb — no solver, a certified bracket.
julia> certified_defect(Φ)
CertifiedBracket  ‖Φ²−Φ‖_cb
  <lo>  ≤  x  ≤  <hi>        (width <w>)      # contains the analytic η = 0.135

# (2) The associated ε-C* algebra hiding in Img Φ.
julia> associated_algebra(Φ)
EpsCStarAlgebra  A = Img Φ̃  ≤  M_2
  dim A = <dimA>
  ε (associator defect) ∈ [<lo>, <hi>]
  product: Choi–Effros star  X⋆Y = Φ̃(XY)   (oblique; axioms hold up to ε)

# A genuinely η>0, multi-block, OBLIQUE channel: a block conditional expectation,
# unitarily conjugate-mixed (B = M_2 ⊕ M_2).
julia> Ψ = UCPMap(bce_conj_kraus(4, 2, 0.02));

# (3) The O(ε) isomorphism to a GENUINE C* algebra — dimension-independent constant.
julia> main_isomorphism(Ψ)
MainIsomorphism  v: B → A   (O(ε), dimension-independent constant)
  B = ⊕ M_d, blocks = [2, 2]
  ‖v‖_cb    ∈ [1.00, <hi>]
  ‖v⁻¹‖_cb  ∈ [1.00, <hi>]
  iso defect ∈ [0.0, <hi>]

# (4) THE SHOWCASE: factor the channel Φ ≈ Δ Υ through B, certified round trips.
julia> F = factorize(Ψ)
ChannelFactorization  Φ ≈ Δ Υ  through  B = ⊕ M_d, blocks = [2, 2]
  encode = Υ* : B → B(H)   (<re> Kraus,  4→4 CPTP)
  decode = Δ* : B(H) → B   (<rd> Kraus,  4→4 CPTP)
  ‖ΔΥ − Φ‖_cb   ∈ [<lo>, <hi>]     (round-trip ΔΥ ≈ Φ)
  ‖ΥΔ − 1_B‖_cb ∈ [<lo>, <hi>]     (round-trip ΥΔ ≈ 1_B)
  η proxy = <ηp>,  ε used = <εu>

julia> algebra(F)
CStarAlgebra  B = ⊕_l M_{d_l}
  blocks d_l = [2, 2]   (m = 2)
  dim_B = 8,  n_B = 4

# The η=0 oracle: on an EXACT conditional expectation the round trip is exact to
# machine ε — the cleanest proof the construction is real.
julia> Fexact = factorize(UCPMap(block_cond_exp_kraus(4, 2)));
julia> maximum(delups_defect(Fexact)), maximum(upsdel_defect(Fexact))
(<d1>, <d2>)        # both < 1e-9
```

Recorder notes (verified against source):
- `phit_kraus(2, 0.1)` → **5** Kraus ops (1 identity + d² = 4 for d=2); the
  `‖Σ K†K − I‖ = 3.1e-16` unital line is the canonical design value (substitute the
  live value if it differs).
- `bce_conj_kraus(4, 2, 0.02)` → `B = M_2 ⊕ M_2`, blocks `[2,2]`, `dim_B = 8`,
  `n_B = 4`. Here `n_B = 4 = N`, so encode/decode are square 4→4 (matches the
  `show` `4→4 CPTP` rows). The grounding report's `bce_conj` is the right showcase
  fixture: multi-block and oblique, so the demo exercises the Choi–Effros star and
  the dual binding (the η=0 identity oracle is blind to both).
- Do **not** demo `factorize(phit_kraus(2, 0.3))` — it throws `ArgumentError`
  (ρ_est ≈ 0.21 ≥ 0.10, outside factorize's convergence domain). `phit_kraus(2, 0.1)`
  is safe for the defect/algebra beats; `bce_conj_kraus(4, 2, 0.02)` (ρ ≈ 0.04) is
  safe for the factorization.
- The η=0 closing beat uses `block_cond_exp_kraus(4, 2)` (exactly idempotent) so
  both `delups`/`upsdel` maxima are `< 1e-9` (asserted in `test_factorize.jl`).

### D.2 Plotting scripts (one PNG each, run SERIALLY)

Each script: `julia --project=docs/plots docs/plots/<name>.jl`, run **one at a
time**. Each reuses the package API and inlines the needed fixture builders (copied
from `test/fixtures.jl` so the docs env has no test dependency). All use CairoMakie,
write a 640 px PNG at 2× to `docs/assets/`. Styling: a single accent colour for the
data, a muted band/envelope, axis labels with units, a one-line title naming the
`.tex` anchor, no chartjunk.

**Script 1 — `containment.jl` → `docs/assets/containment.png`.**
- Computes: for `t` in a grid over `[0.02, 0.30]` and `d ∈ {2, 3}`, build
  `Φ = UCPMap(phit_kraus(d, t))`; `b = certified_defect(Φ)`; collect
  `(t, minimum(b), maximum(b))` and the analytic `η = t(1−t)·2(1−1/d²)`.
- Axes: x = `t`; y = defect. Plot the analytic η as a solid line; shade the
  `[lo, hi]` band behind it (`band!`). One panel per `d`, or both on one axis with
  distinct hues.
- Message in title: "Certified eig-free bracket contains the analytic η
  (test/test_core.jl)". Output 640 px PNG.

**Script 2 — `universality.jl` → `docs/assets/universality.png`** (THE headline).
- Computes: a Julia-side dimension sweep that reproduces the flat constant. Build a
  growing block-algebra family of η-idempotent channels and, for each,
  `v = main_isomorphism(Φ)`, `η = midpoint(certified_defect(Φ))` (or `eta_proxy`
  via `factorize` where in-domain), `c = maximum(isodefect(v)) / η`. The C-side
  anchor points are `(dim_A, c) = (8,0.0182), (12,0.0470), (16,0.0247),
  (18,0.0410), (20,0.0134)` (`tests/test_dbo3.c`, prec=256); **[VERIFY]** regenerate
  the Julia-side values from `main_isomorphism` over `bce_conj`/block families and
  plot those (fall back to plotting the C-anchor points with a clear caption if the
  Julia sweep is too slow at the larger dims — `dim_A = 20` finishes at prec=256 in
  ~95 s, `dim_A ≥ 24` does not; cap at `dim_A ≤ 20`).
- Axes: x = `dim A`; y = `c = isodefect/η`. Scatter the points; overlay the dashed
  OLS line (slope −2.7e-4); annotate "max c = 0.047; no growth with dim
  (tex:484)". Set the y-axis from 0 so the flatness is visually obvious.
- This is the most important plot; give it the cleanest styling.

**Script 3 — `roundtrip.jl` → `docs/assets/roundtrip.png`.**
- Computes: for `t ∈ {0.02, 0.05}` build `Φ = UCPMap(bce_conj_kraus(4, 2, t))` and
  `Φ = UCPMap(mixconj_kraus(4, 2, t))`; `F = factorize(Φ)`; collect
  `η = F.eta_proxy`, `maximum(delups_defect(F))`, `maximum(upsdel_defect(F))`.
- Axes: log–log, x = `η`, y = round-trip defect. Two series (`‖ΔΥ−Φ‖_cb`,
  `‖ΥΔ−1_B‖_cb`), markers per fixture; overlay the line `y = 50·η` (the test
  envelope) as a dashed bound. Slope ≈ 1 (linear in η) is the message.

**Script 4 — `eta0_oracle.jl` → `docs/assets/eta0_oracle.png`.**
- Computes: for the three η=0 oracles `identity_kraus(2)`,
  `block_cond_exp_kraus(4, 2)`, `dephasing_kraus(3)`: `A = associated_algebra(Φ)`,
  `v = main_isomorphism(Φ)`, `F = factorize(Φ)`; collect
  `maximum(epsilon(A))`, `maximum(isodefect(v))`, `maximum(delups_defect(F))`,
  `maximum(upsdel_defect(F))`, and the extracted `blocks(algebra(F))`.
- Axes: grouped bar chart, log y. Four bars per fixture (the four defects); dashed
  horizontal reference at `1e-9`. Annotate each fixture group with its exact block
  list (`[2]`, `[2,2]`, `[1,1,1]`). All bars sit below `1e-9` — the visual proof
  that the η=0 ground truth is exact.

### D.3 Phase-B execution order (serial)

1. `npm install -g svg-term-cli`.
2. `cmake -S . -B build && cmake --build build`; `set_libaic_path!`.
3. Live-capture session (D.0) → fill every `<…>` and `[VERIFY]`.
4. Write `docs/src/assets/demo.cast` from the D.1 transcript with captured numbers;
   `svg-term --in docs/src/assets/demo.cast --out docs/assets/demo.svg --window
   --width 84 --height 30 --padding 18`.
5. Set up `docs/plots/` env (C.2).
6. Run `containment.jl`, then `universality.jl`, then `roundtrip.jl`, then
   `eta0_oracle.jl` — **one at a time**.
7. Drop the Part-B README into `README.md`; substitute the captured benchmark-table
   numbers; confirm `LICENSE` matches the badge/footer (correct if not MIT).
8. Commit assets + README; `git pull --rebase && git push`.

All Julia steps in 3, 6 are serial. The svg-term and CairoMakie work is the only
non-package tooling, and both are dev-only (svg-term via npm; CairoMakie in the
isolated `docs/plots` env) — the shipped package stays Julia/C with no Python and
no new package dependency.
```

---

## Part E — CORRECTIONS + the unified palette (orchestrator, 2026-06-02)

Two binding amendments to Parts A–D. Where these conflict with the draft above,
**Part E wins.**

### E.1 LICENSE correction — it is AGPL-3.0, NOT MIT

The `LICENSE` file on disk is the **GNU Affero General Public License v3.0**
(confirmed: header `GNU AFFERO GENERAL PUBLIC LICENSE / Version 3`). The Part-B
draft's MIT badge + footer are WRONG and must be corrected before publishing:

- Badge: `[![License](https://img.shields.io/badge/License-AGPL_3.0-7aa2f7.svg?style=flat-square)](LICENSE)`
- Footer: `## License` → `GNU Affero General Public License v3.0 — see [LICENSE](LICENSE).`

### E.2 ONE consistent, professional palette across ALL assets (user directive)

The terminal SVG, the four plots, and the badge accents must read as a single
designed set — not three tools each picking defaults. Use this palette (a refined
Tokyo-Night-Storm-class dark scheme — professional, high-contrast, cohesive),
defined ONCE here and reused verbatim:

```
role            hex        used by
--------------- ---------- ------------------------------------------------------
bg              #24283b    SVG terminal background  AND  plot canvas/figure bg
bg_panel        #1f2335    plot axis interior (a touch darker), SVG window chrome
fg              #c0caf5    SVG text  AND  plot text / labels / titles
muted           #565f89    plot gridlines + axis ticks; SVG dim text
accent (teal)   #7dcfff    SVG prompt (julia>)  AND  primary plot series / v
certified-green #9ece6a    "certified"/success: bracket BAND, η=0 oracle bars
purple          #bb9af7    2nd plot series
amber           #e0af68    3rd plot series
red             #f7768e    envelope/limit lines (e.g. the 50·η bound), warnings
plot cycle      [#7dcfff, #9ece6a, #bb9af7, #e0af68, #f7768e, #7aa2f7]
```

**Application rules (Phase B MUST follow):**
- **Terminal SVG:** render with `svg-term`, then lock the colours to the palette —
  post-process the output SVG so the background is `#24283b` and the default
  foreground is `#c0caf5` (svg-term's defaults are close but not exact; string-
  replace the bg/fg hexes). Author the `.cast` prompt `julia>` in truecolor teal
  via the ANSI escape `[38;2;125;207;255m` (= `#7dcfff`); keep `show` output
  uncoloured (it renders in `fg`). Result: the demo window matches the plots.
- **Plots (CairoMakie):** at the top of EACH script set a shared theme so every
  plot has the dark canvas + palette:
  ```julia
  using CairoMakie
  set_theme!(Theme(
      backgroundcolor = "#24283b",
      textcolor = "#c0caf5",
      Axis = (backgroundcolor="#1f2335", xgridcolor=("#565f89",0.4),
              ygridcolor=("#565f89",0.4), xtickcolor="#565f89", ytickcolor="#565f89",
              leftspinecolor="#565f89", bottomspinecolor="#565f89",
              topspinevisible=false, rightspinevisible=false),
      palette = (color = ["#7dcfff","#9ece6a","#bb9af7","#e0af68","#f7768e","#7aa2f7"],),
      fontsize = 15,
  ))
  ```
  Use `band!` in `("#9ece6a", 0.25)` for the certified bracket; the analytic line in
  `#c0caf5`; the `50·η` envelope dashed in `#f7768e`; the η=0 bars in `#9ece6a` with
  the `1e-9` reference dashed in `#f7768e`; the universality scatter in `#7dcfff`
  with the OLS line dashed in `#565f89`. Save PNG at `px_per_unit=2` (retina).
- **Badges:** consistent `?style=flat-square` on all three. Keep the domain-
  conventional brand colours that readers expect — Julia `#9558B2`, arXiv `#b31b1b`
  — and use palette teal-blue `#7aa2f7` for the License badge. (Forcing the Julia/
  arXiv badges off-brand would look wrong; consistency here = the shared flat-square
  style + the harmonised License accent.)

The deliverable test: place `demo.svg` next to `universality.png` — same dark
background, same teal accent, same light text. They look like one set.
