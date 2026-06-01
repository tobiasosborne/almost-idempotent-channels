# CLAUDE.md — almost-idempotent-channels

> This file is the authoritative instruction set for this repository. Read it
> top to bottom at the start of every session and after any context
> compression. The rules drift out of working memory faster than you think;
> that is why they are numbered. It is paired in spirit with the house rules of
> `../su2-fft` (the C+FLINT+Julia sibling), `../BennettVM.jl` (the Laws), and
> `../arithmetic-quantum-mechanics` (sharded-source discipline).

## What this is

An implementation — as **constructive, finite-dimensional algorithms** — of
the results in:

> Alexei Kitaev, *Almost-idempotent quantum channels and approximate
> $C^*$-algebras*, [`arXiv:2405.02434`](https://arxiv.org/abs/2405.02434)
> (v1, dated February 3, 2025).

The paper takes a unital completely positive (UCP) map $\Phi$ on $\Bo(\calH)$
that is **$\eta$-idempotent** ($\|\Phi^2-\Phi\|_\cb\le\eta$), builds an
associated **$\eps$-$C^*$ algebra** of almost-invariant observables
($\eps=O(\eta)$), proves that any finite-dimensional $\eps$-$C^*$ algebra is
**$O(\eps)$-isomorphic to a genuine $C^*$ algebra** $\calB$ (universal
constant, dimension-independent), and — when $\calA$ comes from a finite-dim
$\eta$-idempotent $\Phi$ — realizes the isomorphism and its inverse by UCP
maps. That yields an **approximate factorization of the channel $\Phi^*$** into
a decoding channel (producing a state on $\calB$) and an encoding channel.

The deliverable is a package:

- **tight C** for the numerical cores (functional calculus, contraction-map
  iterations, polar decomposition, the incremental $\calB$ construction,
  Choi/Stinespring, channel factorization);
- **FLINT / arb** (`acb_mat`, `arb`) for arbitrary precision with **certified
  error balls** — the paper's whole subject is rigorous $O(\eps)$ bounds, so we
  certify them numerically rather than trust `double`;
- a **Julia package** that `ccall`s the C cores (the `../su2-fft` pattern).

## THE MANDATE — constructive algorithms for every result

This is the project's reason to exist; everything else serves it.

**Many proofs in the paper are non-constructive and phrased for possibly
infinite-dimensional spaces. But every *theorem* applies to finite-dimensional
matrices, and in finite dimensions the objects those proofs merely assert to
exist are computable. Our goal is a constructive algorithm for *every* result.**

Concretely, for each result you must separate two things:

1. **The paper's proof technique** — and whether it is constructive. The paper
   leans on non-constructive tools: the Lefschetz–Hopf fixed-point theorem and
   H-space homotopy (existence of a nontrivial projection, §6); holomorphic
   functional calculus via contour integrals (§3); Haar-measure diagonals /
   Hochschild cohomology (error reduction, §2/§9); existence statements with no
   procedure attached.

2. **The constructive finite-dim algorithm we implement** — which need not
   follow the proof. Examples of the move:
   - *Nontrivial projection* (proved topologically via Lefschetz–Hopf) →
     in finite dim, construct one directly: spectral split of a Hermitian
     element, or `θ(2P−I)` applied to a near-idempotent, or minimize
     $\|P^2-P\|$ over Hermitian $P$ bounded away from $\pm I$.
   - *Holomorphic functional calculus contour integral* → eigendecomposition,
     or a Newton/Denman–Beavers iteration for $\sgn(X)$, or a rational
     (Zolotarev) approximation with a certified arb error.
   - *Haar-measure diagonal* $D=\int dU\,(U^\dag\otimes U)$ → a finite average
     over the regular representation / an explicit fixed element computable in
     finite dim.
   - *Implicit/inverse function theorem* (§3, `lem_invfun`) → it is **already
     constructive**: the Banach contraction iteration $x_n=g_y(x_{n-1})$ is the
     algorithm; implement it as such, with the paper's $c<1$ giving the
     convergence rate.

   **The paper's $O(\eps)$ bound is the algorithm's correctness specification.**
   The constructive routine must produce an output meeting the bound the
   theorem promises; that bound is the assertion the cross-check verifies.

If a result resists constructivization, that is a finding to escalate (a stop
condition), not a corner to paper over. File it; do not fake it.

## Read order

For any task, in order. Refuse to add code or math content before steps 1–2.

0. `HANDOFF.md` — five-minute orientation: current state, what's built, what's
   next, and the key findings/decisions. Read first when resuming work.
1. This file (`CLAUDE.md`).
2. `paper/src/approximate_algebras.tex` — the verbatim arXiv source. **Canonical
   ground truth.** When notes, memory, or a shard disagree with it, the
   `.tex` wins.
3. The relevant shard under `paper/shards/` — the per-region map (statement +
   its proof + dependencies + constructivization analysis + implementability).
4. `MODULE_PLAN.md` — the DAG of C/Julia modules and which result each realizes.
5. `ALGORITHM.md` (once it exists) — the canonical math+code narrative for the
   module you are touching, every formula cited to a `.tex` line.
6. `paper/FINDINGS.md` — the **living log of paper issues**: `.tex` typos/formula
   errors, non-constructive steps + our constructive routes, load-bearing
   subtleties (the "tests that can't fail" our hostile reviews keep catching), and
   open escalations. Skim before touching a new region of the paper so you inherit
   the known issues; **append to it** whenever research/implementation/review finds
   a new one, and cite `paper/FINDINGS.md §Xn` from the source comment where it bites.

If a claim in code is not anchored to a `paper/src/approximate_algebras.tex`
line number, it is undocumented; fix that before changing the code. When the code
deviates from a `.tex` formula (a typo, a non-constructive route, an oblique-image
subtlety), that deviation **must** have a `paper/FINDINGS.md` entry it cites.

## Architecture

Mirrors `../su2-fft`. Three layers, cross-validated against each other:

```
Julia (AlmostIdempotentChannels.jl)   <- ccall surface, high-level API, tests
        │  ccall
        ▼
C cores (src/*.c, include/*.h)         <- algorithms; two number paths:
        ├── double path  (BLAS/LAPACK)     fast, the prec=53 anchor
        └── arb path     (FLINT acb/arb)   certified error balls
```

Every nontrivial routine exists (eventually) in **both** number paths, and the
arb path at `prec=53` must agree with the double path — that is a primary
cross-check. The arb path is what *certifies* the paper's $O(\eps)$ bound; the
double path is what makes the test suite fast.

## The Laws (unconditional)

**Law 1 — Ground truth before code.** Every implementation decision cites the
`.tex` line it implements, with a verbatim copy of the equation/statement in a
source comment (the `../su2-fft` house style):

```c
/* approximate_algebras.tex:528 — Proposition prop_P:
 *   Ptilde = theta(2P - I),  theta(X) = (I + sgn(X))/2;
 *   then Ptilde^2 = Ptilde, ||Ptilde - P|| <= ||2P-I|| * O(delta),  delta<1/4.
 */
```

No paraphrasing the construction from memory. The paper is on disk — open it.

**Law 2 — Reuse before reinvention.** Do not reinvent numerical linear algebra.
Use **FLINT/arb** (`acb_mat_*`, matrix sign, eigenproblem helpers) and
**LAPACK** (double path: `zheev`, `zgesvd`, `zgees`) deliberately. The user
picked FLINT for certified arbitrary precision; do not propose alternatives or
"modernize" the arb path until you have shipped a measured improvement. Wrap,
don't reimplement.

**Law 3 — Constructivize, don't transcribe.** (See THE MANDATE.) The artifact
is a finite-dimensional algorithm meeting the paper's bound, not a literal
rendering of a non-constructive existence proof. Always record, per result, the
paper's proof technique *and* the constructive route you implemented.

**Law 4 — Audition, don't presume.** No algorithm is presumed fit for any
result. The project is exploratory: every algorithm choice (matrix-sign route,
projection route, polar method, …) is a **candidate auditioned via red-green TDD
against correctness AND performance benchmarks**, then selected on the evidence.
The shards and `MODULE_PLAN.md` list candidate routes (e.g. projection Route
A/B/C); treat them as the audition slate, not a ranked decision, even where a
subagent wrote "recommended." Keep multiple candidates alive until a benchmark
retires one. Performance is a first-class gate, not an afterthought (the
`su2-fft` `PROFILING.md` discipline).

When **no candidate dominates** on the current metric, do **not** pick
arbitrarily — **add performance dimensions** (wall time, iteration count, peak
memory, accuracy / certified-ball tightness at fixed precision, robustness /
convergence basin on adversarial inputs) and look for **Pareto domination**. If
none dominates across the regimes of interest, keep the entire **Pareto
frontier** in the code (named, dispatchable) and select per regime, recording
the frontier and the regime boundaries. Only declare a single default when one
candidate Pareto-dominates everywhere of interest.

## The Rules (numbered, non-negotiable)

0. **Laws 1–3 apply.** Always.

1. **Cite everything.** Every math routine cites `(a)` the `.tex` line and
   `(b)` the verbatim equation/statement. This is how a routine survives a
   rewrite.

2. **Skepticism.** Be skeptical of your own assumptions, subagent output,
   previous-session claims, the existing tests, and the paper itself (LaTeX
   typos happen; so do off-by-one constants). Verify against the `.tex`.
   LLM summaries of operator-algebra proofs are especially untrustworthy —
   re-read the source.

3. **All bugs are deep.** No bandaids. A constant that makes one test pass at
   $n=4$ but produces garbage at $n=20$, or a sign routine that works for
   well-separated spectra and NaNs near a spectral gap of size $\eps$, is a
   future incident with a long fuse. Find the root cause.

4. **Fail fast, fail loud.** `assert()` invariants; never silently return. If a
   precision is insufficient, a spectral gap is too small to resolve, a
   contraction constant $c\ge1$, or an input violates an $\eps$-axiom bound —
   abort with a clear message at the call site. Corrupted output is worse than
   a crash. In the arb path, a ball that has lost all precision is a loud
   failure, not a silent widening.

5. **"Runs without errors" is NOT a passing test.** Every test asserts an
   invariant against a known value or bound: the paper's $O(\eps)$ bound holds;
   two independent routines agree to tolerance; an exact special case
   reproduces a closed form; a round-trip (encode∘decode, or the
   $O(\eps)$-isomorphism and its inverse) returns to identity within bound.

6. **Cross-checks > unit tests.** The strongest tests here:
   - **double vs arb@prec=53** for the same routine (the `su2-fft` anchor);
   - **bound certification:** the arb path proves the output satisfies the
     theorem's $\le O(\eps)$ inequality for random inputs at the boundary of the
     hypotheses;
   - **exact-idempotent oracle:** when $\eta=0$ the constructions must reduce to
     the exact Choi–Effros / `th_idemp_structure` decomposition — a genuine
     $C^*$ algebra, $0$ defect. This is the cleanest ground truth in the paper.
   Unit tests catch typos; cross-checks catch algorithmic errors. Prefer the
   latter, and add the cross-check **before** the code.

7. **TDD — two valid shapes.**
   - *Spec-from-scratch:* RED → GREEN → refactor. Write the failing assertion
     (usually a bound or a cross-check) first.
   - *Port-and-verify:* for a routine realizing a paper result, implement the
     constructive algorithm, capture its invariants as tests, **mutation-prove**
     the tests catch regressions (perturb the impl, confirm RED, restore), and
     cross-validate against an independent oracle.
   The discipline is "a test has caught a real regression," not "the test was
   committed first."

8. **Get feedback fast.** Run the relevant test after every non-trivial change;
   check every ~50 LOC, don't code blind for 500. A quick `julia --project -e`
   probe or a single C test binary is fine.

9. **Tiered workflow + reviewer-gating.** Scale effort to change size:
   - *Trivial* (≤5 LOC; typo/comment): direct fix, reviewer-exempt.
   - *Small* (≤30 LOC; one function using existing abstractions): write the
     test, write the code, one reviewer subagent before declaring done.
   - *Core* (new algorithm/path/result, >30 LOC, math routine): research the
     `.tex` + shard first; decide the constructive route; write the cross-check
     before the code; hostile reviewer subagent always. For contested
     constructivization choices, run 2–3 research subagents independently first.

10. **LOC limit ~200 per source file.** When a module approaches it, split.
    Inherited from the sibling projects' Feynfeld convention.

11. **Literate programming.** Each source file opens with a docstring that
    explains *why* it is shaped as it is: which `.tex` result it realizes, the
    paper's proof technique vs the constructive route taken, which defensive
    checks guard which degeneracy, the precision argument. A fresh reader should
    read the file top-to-bottom like a chapter. When you touch the algorithm,
    update `ALGORITHM.md`; the docs are the contract.

12. **No emojis, no marketing.** Code and docs read like SQLite / TigerBeetle
    docs. Concrete numbers always: "max-diff 3e-16 at $n=12$, $\eps=10^{-6}$",
    not "very accurate". No "robust", "blazing fast", "production-grade".

13. **No GitHub CI.** Quality gates run locally. The user has rejected automated
    remote runs across all their projects; do not create `.github/workflows/`,
    do not propose "add CI" tasks.

14. **Persistent tracking = beads (initialized 2026-05-28, prefix `aic`).** `bd`
    is the only *persistent* cross-session tracker; the full module plan is
    registered as issues — `bd ready` (entry point: `T-build`), `bd show <id>`,
    `bd dep tree`. `TaskCreate`/`TaskUpdate` are allowed for *ephemeral
    in-session* state (user directive 2026-05-28) — **this supersedes the blanket
    "do NOT use TaskCreate" in the auto-generated Beads block at the bottom of
    this file.** No markdown TODO lists. The harness auto-memory (the external
    `memory/` directory) is a separate mechanism and is unaffected by the
    auto-block's "no MEMORY.md files" rule. Never `bd init --force`.

15. **Don't replace FLINT.** See Law 2. The arb path is the point of the project.

16. **Re-read this file** at session start, after `/clear`, after any
    compaction.

## Domain hallucination-risk callouts

Mistakes that look right but aren't. When you catch yourself about to do one,
stop and re-check the cited `.tex` line.

- **The $\eps$-$C^*$ axioms hold only *up to* $\eps$.** Associativity,
  $\|XY\|\le(1+\eps)\|X\|\|Y\|$, the $C^*$ identity, and the unit laws are all
  approximate (`.tex:407–439`). Do not assume exact associativity or an exact
  unit in any routine operating on $\calA$. The exact unit is only available
  after the $O(\eps)$-change of `prop_unit` (`.tex:672`).

- **The product on $\calA=\Img\Phi$ is the Choi–Effros product
  $X\star Y=\Phi(XY)$, not $XY$** (`.tex:342`). Plain operator multiplication
  leaves $\Img\Phi$. Every multiplication in the associated algebra goes through
  $\Phi$.

- **cb-norm ≠ operator norm.** $\eta$-idempotence is $\|\Phi^2-\Phi\|_\cb\le\eta$
  (`.tex:347`), a supremum over $1_{\Ma{n}}\otimes\Phi$ for all $n$. Estimating
  it requires the Choi-matrix / ampliation, not just $\|\Phi^2-\Phi\|$. Tensor
  extensions (§10) exist precisely because cb-norm sees all ampliations.

- **Spectra are perturbation-sensitive; do not trust naive eigenvalues.** The
  paper's own example (`.tex:540–544`): a $3\times3$ matrix with a $t$ entry has
  eigenvalues $\sim t^{1/3}$. For non-normal $X$, $\spec(X)$ is fragile — this is
  *why* the paper avoids holomorphic calculus and *why* we want certified arb
  balls. Prefer constructions on Hermitian/normal elements; certify everything
  else.

- **$\sgn(X)$, $|X|$, $\theta(X)$ need $\|X^2-I\|<1$ (resp. $\|X^2-x_0^2I\|<x_0^2$)**
  (`.tex:516`, `.tex:520`). Outside that ball the functional calculus is
  undefined; assert the precondition.

- **$\delta$-homomorphism vs $\eps$-algebra: two different small parameters.**
  $\eps$ measures how far $\calA$ is from a $C^*$ algebra; $\delta$ measures how
  far a map $v$ is from multiplicative (`.tex:443–455`). The main proof juggles
  several ($\eps,\delta,\delta_0,\eta$); keep them distinct in code and never
  silently identify them.

- **The isomorphism is $O(\eps)$, not exact.** `th_main` (`.tex:460`) gives an
  $O(\eps)$-isomorphism, with a *universal, dimension-independent* constant.
  Tests must assert the bound, not exact equality, and must check the constant
  does **not** grow with $n$ (the naive Haar/cohomology route fails exactly
  because its error $\propto n$ — `.tex:484`).

- **UCP (observables) vs CPTP (states) are dual; mind the direction.** $\Phi$
  acts on $\Bo(\calH)$; the channel $T=\Phi^*$ acts on states. `Enc`/`Dec` are
  duals of $\Gamma\Co_\calM$ / $\Delta$ (`.tex:334`). Getting the dual backwards
  silently transposes everything.

- **$\Img\Phi=\Ker(1-\Phi)$ only when $\Phi$ is *exactly* idempotent**
  (`.tex:344`). For $\eta$-idempotent $\Phi$ this is itself approximate — the
  associated algebra construction (§12.1) is the careful version; don't shortcut
  it with the exact identity.

## Arbitrary precision — why, and the cross-check ladder

The paper certifies rigorous bounds; arb lets us do the same numerically.
Precision is load-bearing where: (i) $\eps$ is tiny and the bound involves
differences of $O(1)$ quantities ($(XY)Z-X(YZ)$, $\Phi^2-\Phi$); (ii) a routine
inverts a near-singular operator ($(\La_{|X|}+\Ra_{|X|})^{-1}$, `.tex:615`);
(iii) a spectral gap is comparable to $\eps$. The cross-check ladder, weakest to
strongest:

1. double path internal sanity (norms, Hermiticity).
2. double vs arb@prec=53 agreement (~1e-12).
3. exact-idempotent oracle: $\eta=0 \Rightarrow$ exact $C^*$ structure, zero
   defect.
4. arb bound certification: the output ball provably satisfies the theorem's
   $\le O(\eps)$ inequality at the hypothesis boundary.

## Paper structure map (canonical line numbers)

`paper/src/approximate_algebras.tex`, 2928 lines. Two proofs are **deferred far
from their statements** — keep them connected:

| § | Title (label) | Lines | Key results | Proof location |
|---|---|---|---|---|
| 1 | Motivation; encode/decode (`sec_motivation`) | 251–402 | `prop_hom_structure`(p.273); `th_idemp_structure`(stmt) | **proof @2055 (§11.3)** |
| 2 | Main theorem + strategy (`sec_main_thm`) | 403–492 | Defs of $\eps$-$C^*$, $\delta$-hom; `th_main`(stmt) | **proof @1414 (§9)** |
| 3 | Analytic tools | 493–631 | `prop_P`(sgn/θ projection); `lem_invfun`(contraction); `rem_X2` | inline |
| 4 | Basic properties of $\eps$-Banach algebras | 632–689 | `prop_unit`(exact-unit fix) | inline |
| 5 | Approximate unitary group (`sec_unitaries`) | 690–914 | `lem_U_delta`; `lem_gV`; `prop_polar` | inline |
| 6 | Nontrivial projection (`sec_projection`) | 915–1051 | `lem_nontriv_projection`(Lefschetz–Hopf); `prop_H-group` | inline (@934,@1023) |
| 7 | Subspaces from projections (`sec_subspaces`) | 1052–1189 | `lem_alpha`,`lem_PQ_Hilb`,`lem_PQR`,`lem_1d_proj` | inline |
| 8 | $\delta$-homomorphisms | 1190–1320 | `prop_delta_hominc`; `lem_approx`; `cor_improvement`(error reduction) | inline (@1256) |
| 9 | Proof of main theorem (`sec_proof_main`) | 1321–1446 | `lem_merging`,`cor_merge_sum`,`lem_add_dim`,`lem_extension`; **proof of `th_main`** | inline (@1414) |
| 10 | Tensor extensions (`sec_tens_ext`) | 1447–1561 | `def:opspace`; `prop_inc_ext`; `lem_approx_ext`; `th_main_ext` | inline |
| 11 | Details on UCP maps (`sec_channels`) | 1562–2161 | Choi/Stinespring `prop_KLHG`; `lem_carrier`,`cor_carrier`; `lem_idemp`; **proof of `th_idemp_structure`**(@2055); `prop_Gamma` | inline |
| 12 | Almost-idempotent UCP maps (`sec_almost_idemp`) | 2162–2900 | `th_almost_idemp`(assoc $\eps$-$C^*$); `th_factorization`(headline); `lem_RC` | @2208/@2239/@2639 |

## File map (planned; create as needed)

```
CLAUDE.md                      # this file
MODULE_PLAN.md                 # module DAG, result→module map, build order
ALGORITHM.md                   # canonical math+code narrative (grows with code)
paper/
  source.tar.gz                # verbatim arXiv e-print (do not edit)
  src/approximate_algebras.tex # canonical ground truth
  src/*.sty, src/*.pdf         # style + figures
  shards/shard-*.md            # per-region maps (statement+proof+constructivization)
include/aic.h                  # public C API + storage-layout contract
src/                           # C cores, double + arb paths, ≤200 LOC each
tests/                         # C cross-checks (double vs arb, bound certification)
julia/AlmostIdempotentChannels.jl/   # Julia package, ccall wrappers, runtests.jl
```

## Build & test

The C core builds with **CMake** (bead aic-95g.1); the hand-written `Makefile`
is now a thin convenience wrapper that forwards to it. See `BUILDING.md` for the
full recipe (dependencies, install, the `find_package(AIC)`/pkg-config consumer
paths, and the Debian flint.pc / LAPACKE gotchas).

```bash
# C cores (CMake is canonical; `make` forwards to it)
cmake -S . -B build && cmake --build build   # or: make
ctest --test-dir build -L fast               # sub-second laptop gate (or: make test-fast)
ctest --test-dir build                        # full suite, parallel (or: make test)
ctest --test-dir build -R test_X              # one test
cmake --install build --prefix <prefix>       # set the prefix at CONFIGURE time too

# Julia package
julia --project=julia/AlmostIdempotentChannels.jl -e 'using Pkg; Pkg.test()'
```

A single-compile OBJECT library yields `libaic.so` + `libaic.a` from ONE compile
of each `src/*.c`. Every test carries a CTest `TIMEOUT` (fail-loud, no hangs);
the `fast`/`slow` labels split a sub-second iteration gate from the arb/MOSEK-
heavy drivers (`-L fast` is the laptop default; full `ctest` is the gate).
Dependencies: C11 compiler, CMake ≥ 3.24, `libflint-dev` (FLINT ≥ 3.0, brings
arb), `liblapacke-dev`/LAPACK/BLAS. Strict per-target flags following the
sibling: `-Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes -O2 -g -std=c11`
(keep `-Wshadow`; it has caught real bugs in the sibling). NO `-Werror`.

## Probe/sweep hygiene — don't hang the session (2026-05-31 incident)

A throwaway hostile-review probe hung an entire session and the in-flight
review's state was lost on exit. The avoidable mistakes, recorded so no future
agent repeats them:

- **Scratch/probe `.c` files SHOULD NOT match `tests/test_*.c`.** Both the CMake
  `tests/CMakeLists.txt` (`file(GLOB ... test_*.c)`) and the old Makefile glob
  auto-discover that pattern into a registered test, so a stray
  `tests/test_zprobe.c` still gets built and run by `ctest`/`make test`. Since
  the CMake migration every test carries a `TIMEOUT 300`, so a non-converging
  probe now fails LOUD after 5 min instead of hanging the gate forever (the
  aic-xo0 mitigation) — but still put throwaway code outside `tests/` (a `/tmp`
  build, or a name that does not start with `test_` that you delete
  immediately), never as `tests/test_zprobe.c`.
- **Stay inside the spectral `ρ(Φ²−Φ) < 1/4` regime — and know what actually
  leaves it.** `Φ̃ = θ(2Φ−1)` needs the SPECTRAL basin `ρ(Φ²−Φ) < 1/4`
  (`.tex:520-525`; since bead aic-8hz the guard is the eig-free Gelfand `ρ`, NOT
  the stricter op-norm `‖Φ²−Φ‖`). The pipeline now **fails loud** out-of-basin —
  it aborts in ~0.2 s at `aic_prop_P` and at the `aic_assoc_regularize` entry
  guard (bead aic-xo0); it does NOT hang (the old hang was the pre-aic-8hz
  op-norm route). **Correction (FINDINGS §C15):** `make_mixconj*` CANNOT leave
  the spectral basin at ANY `t` — it mixes two spectrally-idempotent maps, so
  `ρ(Φ²−Φ) < 1/4` for all `t` (measured `ρ≈0.165` even at `t=0.45`); the earlier
  "`t=0.45` is out-of-regime" claim was the stale op-norm proxy and is FALSE. A
  genuine out-of-regime fixture needs a non-idempotent spectrum, e.g. unitary
  conjugation by a reflection `Φ(X)=UXU†`, `U=diag(1,−1)` (`ρ(Φ²−Φ)=2`).
- **Bound every exploratory sweep.** Wrap each heavy pipeline build in a
  per-point `timeout`, or run it backgrounded, so one non-converging point
  cannot hang forever. Note `pkill -f "make test"` self-matches the shell
  (exit 144) — use `pkill -x make`.

## Session completion ("landing the plane")

This repo is **under git + beads** (initialized 2026-05-28). A **git remote now
exists** (`origin` → `git@github.com:tobiasosborne/almost-idempotent-channels.git`,
`master` tracks `origin/master`; user directive 2026-05-30: "push should be done
regularly and should work"). So **push after every commit** — `git pull --rebase &&
git push` — and certainly at session end (adopt the sibling push protocol: work is
not "done" until pushed). Before committing, `bd export > .beads/issues.jsonl` so
the bead state lands in the commit (the pre-commit hook also syncs it); the JSONL
in git IS the bead persistence — **no `bd dolt` remote is configured**, so `bd dolt
push` is a no-op/N/A (do not invent a dolt remote unless asked). At session/commit
boundaries: file/close beads for work surfaced, run the quality gates, commit
atomically with `.tex`-line provenance, push, and leave the working tree clean
(the untracked `.claude/docs`/`.claude/tools` are session-start hook artifacts —
leave them untracked, do not commit them).

## Stop conditions (escalate to the user)

- A result you cannot constructivize in finite dimensions (the central risk —
  see THE MANDATE). File it with the specific obstruction.
- A bound that grows with dimension where the paper claims it is universal — you
  have almost certainly taken the naive route the paper warns against
  (`.tex:484`).
- A `.tex` statement that appears to contradict another, or whose constant looks
  like a typo — flag it; do not silently "fix" the math.
- A precision wall: the arb path cannot certify a bound at any feasible `prec`,
  suggesting the construction is numerically unstable as written.
- A subagent returning a non-actionable "looks fine" with no per-result signoff
  — re-request with a checklist; if still vague, escalate.

When escalating, attach: the result label, the `.tex` line, the specific
obstruction, and the reproducible command.


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:7510c1e2 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
