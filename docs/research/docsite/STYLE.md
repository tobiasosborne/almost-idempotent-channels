# Documentation style contract — AlmostIdempotentChannels.jl docs site

Every page-writing subagent MUST follow this file exactly. It exists so ~24
pages written by different authors read as one voice and all build clean under
Documenter with `doctest = true`. Deviations are bugs.

The four source briefs live next to this file in `docs/research/docsite/`:
`brief-math.md`, `brief-api.md`, `brief-architecture.md`, `brief-examples.md`,
`brief-structure.md`. Draw facts and snippets from them; do not invent.

---

## 1. Voice (CLAUDE.md Rule 12 — non-negotiable)

- Read like SQLite / TigerBeetle docs: precise, calm, concrete.
- **No emojis. Anywhere.**
- **No marketing words:** forbidden — "powerful", "robust", "blazing", "seamless",
  "production-grade", "state of the art", "cutting-edge", "simply", "just",
  "effortless", "elegant" (as praise), "unique". Do not write "X is the best/only…".
- **Concrete numbers always**, never adjectives for performance/precision. Write
  "219 tests pass in ~2m09s with no solver", not "fast and well-tested". Write
  "the bracket is loose by design, hi/lo ~ 2n", not "the bracket is a bit loose".
- Acknowledge tradeoffs in the SQLite idiom: "solver-free by default; loose
  bracket by design (hi/lo ~ 2n); tight only with MOSEK."
- Prefer the active voice and short sentences. One idea per paragraph.

## 2. Diátaxis mode discipline (the one rule that matters most)

Each page is exactly ONE mode. Never mix.

- **Tutorial** (learning): a guided run that produces a visible result. May NAME a
  function and say "this returns a rigorous bracket". MUST NOT explain *why* it is
  rigorous — that is one sentence of Explanation; link to it instead.
- **How-to** (task): numbered steps to accomplish one real goal. No teaching, no
  theory. Opens with a one-line assumption ("This guide assumes you have a
  `UCPMap` in scope; see [Quick start](@ref)").
- **Explanation** (understanding): the why and how-it-fits. No install steps, no
  API signatures, no step-by-step. May contain formulas (with paper anchors).
- **Reference** (lookup): signatures, fields, constraints. No narrative.

If a sentence explains *why* something is true, it belongs in Explanation. Move it
there and cross-link.

## 3. Every page opens with a one-sentence contract

First line after the H1 title states what the reader will be able to DO (tutorial/
how-to) or UNDERSTAND (explanation) or LOOK UP (reference) after the page. Example:
"After this guide you can certify a channel's idempotency defect without any SDP
solver." Explanation pages additionally get a "you are here" sentence locating the
page: "This page explains the mathematics behind the five-verb pipeline; to run it,
see the [pipeline tutorial](@ref)."

## 4. Math: PLAIN KaTeX only (Documenter default). No custom macros.

The briefs use the PAPER's LaTeX macros (`\calH`, `\Bo`, `\cb`, `\eps`, `\wt`,
`\Ma`, `\Img`, `\Ker`, `\Tr`, `\sgn`, `\Co`, `\CC`, `\dag`, `\Diamond`, `\spec`,
`\Enc`, `\Dec`). KaTeX does NOT know these. You MUST translate them to plain KaTeX,
matching the existing `docs/src/index.md` convention:

| paper macro | write instead |
|---|---|
| `\calH \calA \calB \calM \calF` | `\mathcal{H} \mathcal{A} \mathcal{B} \mathcal{M} \mathcal{F}` |
| `\Bo(\calH)` | `\mathcal{B}(\mathcal{H})` (or plain `B(H)` in light inline math) |
| `\cb` (subscript) | `_{\mathrm{cb}}` e.g. `\|\Phi^2-\Phi\|_{\mathrm{cb}}` |
| `\Diamond` | `\diamond` |
| `\eps` | `\varepsilon` |
| `\wt\Phi` | `\tilde\Phi` |
| `\Ma{n}` | `M_n` |
| `\Img \Ker \Tr \sgn \spec` | `\operatorname{Img}` `\operatorname{Ker}` `\operatorname{Tr}` `\operatorname{sgn}` `\operatorname{spec}` |
| `\Co_\calM` | `\operatorname{Co}_{\mathcal{M}}` |
| `\CC` | `\mathbb{C}` |
| `\dag` | `\dagger` |
| `\Enc \Dec` | `\mathrm{Enc}` `\mathrm{Dec}` |
| `\star` | `\star` (this one is standard) |
| `\bigoplus_l M_{d_l}` | unchanged (standard) |

Display math uses ` ```math ` fenced blocks (Documenter) or `$$...$$`. Inline math
uses double backtick ``` ``...`` ``` the way `index.md` does, OR single `$...$`.
Match `index.md`: it uses ``` ``\Phi`` ``` inline and ` ```math ` for display.
When in doubt, keep formulas short and inline.

Every nontrivial equation in an Explanation page carries a paper anchor in prose:
"(`approximate_algebras.tex:347`)" or "(arXiv:2405.02434, §3, `prop_P`)". This is
CLAUDE.md Law 1 applied to docs.

## 5. Code examples: executed, never un-run

Three kinds of code block, in priority order:

- **`jldoctest`** — for STRUCTURALLY STABLE output only: dims, block lists like
  `[2, 2]`, type names, boolean flags, small exact integers/strings. Output is
  checked exactly at build time, so it must be reproducible. Use named modules
  (`jldoctest label`) when state must persist across blocks; use `setup = :( ... )`
  to seed definitions. Follow the existing `tutorial.md` convention (`tut`, `tut1`).
- **`@example`** — for ANYTHING float-heavy: a `CertifiedBracket` printout,
  `width`/`midpoint` values, `epsilon(A)`, any defect bracket. Output is embedded
  live (robust to last-digit wobble). Use a named block (`@example pagename`) so
  later blocks on the page see earlier state.
- **Plain ```` ```julia ```` fence (NON-executing)** — ONLY for shell commands and
  for MOSEK examples (which cannot run in the default solver-free build). Mark MOSEK
  blocks clearly as illustrative inside a `!!! tip "With MOSEK installed"`.

ALL the snippets you need are pre-verified in `brief-examples.md` §2, with the
exact expected output and a `[jldoctest]` / `[@example]` tag per block. USE THOSE.
Do not write a code example that is not derivable from `brief-examples.md` or the
existing `docs/src/{index,tutorial}.md`. If you are unsure whether output is stable,
use `@example` (it cannot fail the build on a float).

### The build-breaking pitfalls (brief-examples.md §4) — obey all ten

1. **Fixtures are NOT exported.** `phit_kraus`, `block_cond_exp_kraus`,
   `dephasing_kraus`, `identity_kraus`, `fixed_unitary`, `bce_conj_kraus`,
   `mixconj_kraus`, etc. live in the TEST suite. Every example that uses one must
   INLINE the full constructor body (copy from `brief-examples.md` §2). The
   `@setup` block is the place to hide this boilerplate.
2. **factorize domain:** never call `factorize` on a fixture with ρ ≥ 0.10.
   `phit_kraus(2, 0.3)` THROWS. Safe for factorize: any η=0 oracle,
   `phit_kraus(2, 0.1)`, `bce_conj_kraus(4, 2, 0.02)`, `mixconj_kraus(4, 2, 0.02)`.
3. **`iscptp(decode(F))` is `false` at the default `atol=1e-9`** for any η>0
   multi-block fixture. Show `iscptp(decode(F); atol = 1e-4)` and name the design.
4. **MOSEK examples cannot execute** in the default build — non-executing fences only.
5. **CertifiedBracket floats are not last-digit stable** — `@example`, never `jldoctest`.
6. The `∈` containment of the analytic η needs a 1e-7 prose margin (the test uses
   slack); say "within ~1e-7" rather than asserting bare `analytic in η` as exact.
7. `value(certified_defect(Φ))` is `nothing` on the solver-free path.
8. `blocks(v)` (from `main_isomorphism`) and `blocks(algebra(F))` (from `factorize`)
   are the same data via different objects — don't let setup variable names collide.
9. For `bce_conj(4,2,t)` encode/decode are square (4×4). For the rectangular case
   (n_B=2, N=4) use `mixconj(4,2,0.02)`.
10. `associated_algebra`/`main_isomorphism` ALSO throw on ρ ≥ 1/4 (the prop_P guard);
    only `certified_defect` is safe at any ρ.

## 6. Admonitions: three purposes only (no decoration)

- `!!! warning` — the two hard constraints: factorize's tighter domain; decode's
  `iscptp` tolerance. 
- `!!! note` — the one design choice that surprises: the eig-free bracket is loose
  by design (hi/lo ~ 2n).
- `!!! tip` — "If you have MOSEK installed…".

Never use an admonition to emphasise ordinary important prose (that is `**bold**`),
and never `!!! danger`.

## 7. The three Known Limits — canonical wording (quote consistently)

1. **factorize domain.** "`factorize` builds Υ via a power-series functional
   calculus whose convergence domain is much smaller than the prop_P basin
   ρ(Φ²−Φ) < 1/4. The verb pre-checks a conservative ρ < 0.10 in Julia and throws a
   clean `ArgumentError` — it does not abort the session. `associated_algebra` and
   `main_isomorphism` keep the wider ρ < 1/4 domain; `certified_defect` is always
   safe at any η." (bug `aic-exa.13`)
2. **decode is O(η)-TP.** "The `decode` channel is only O(η)-trace-preserving for
   η > 0 (an internal Choi→Kraus PSD-cone clip; measured clipped mass ≈ 3.7e-6 on
   multi-block oblique fixtures). `iscptp(decode(F))` is `false` at the default
   `atol = 1e-9` and `true` at `atol = 1e-4`. The rigorous certificate is the
   cb-norm round-trip bracket, not `iscptp` at machine tolerance. The encode channel
   and the η = 0 oracle decode are TP to machine precision."
3. **eig-free bracket loose by design.** "The eig-free bracket is loose by design
   (hi/lo ~ 2n); it certifies a value rather than computing it. For a tight bracket
   or the exact value, use the MOSEK extension."

## 8. Cross-links (bidirectional)

Use Documenter `@ref`. Link targets are page section headers or docstring names.
- How-to pages end with: "For why this works, see [the relevant Explanation](@ref).
  For the full signature, see the [API reference](@ref)."
- Explanation pages end with: "To apply this, see [the relevant How-to](@ref)."
- The η=0 oracle is a recurring motif: every page that produces a certified value
  should note "on an exactly-idempotent channel this collapses to 0 — see
  [The η = 0 oracle](@ref)."
- Use `[`certified_defect`](@ref)` to link a function to its docstring.

## 9. The η = 0 oracle motif

The cleanest ground truth in the project (exactly-idempotent Φ ⇒ every defect at
machine ε; measured maxima 4.4e-75 and 3.9e-75 at prec=128). Treat it as the
calibration device the reader returns to. Name it, link it, reuse it.

## 10. Do not duplicate the README

The README is the GitHub front door. The site Home (`index.md`) orients and routes;
it does not reproduce the README's benchmark tables verbatim. Pull the three-sentence
story arc and the five-verb table; link to the README for the install detail is
NOT needed (we have a real install page). Avoid copying numbers into two places that
would then drift — cite the test/source instead.

## 11. Filenames and the pages tree

Write to the EXACT path you are assigned. The `make.jl` pages tree references these
paths; a missing file fails the build. Section headers you create become `@ref`
targets, so give them stable, descriptive names.
