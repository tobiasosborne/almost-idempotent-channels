# Documentation site brief — AlmostIdempotentChannels.jl

Research conducted 2026-06-03.
Sources: diataxis.fr; Documenter.jl manual; JuMP.jl docs; Makie.jl docs;
SQLite documentation; project files (index.md, tutorial.md, api.md, make.jl, README.md).

---

## A. Diataxis cheat-sheet for this package

The four modes differ by user state and user need, not by topic or length.
The rule that separates excellent docs from average ones: **never let two modes
appear on the same page**. A page that teaches AND explains AND lists API
simultaneously confuses all three readers. (Source: diataxis.fr/tutorials/,
diataxis.fr/how-to-guides/, diataxis.fr/explanation/, diataxis.fr/reference/)

| Mode | User state | User need | BELONGS here for this package | DOES NOT belong here |
|---|---|---|---|---|
| **Tutorial** | Learning, no prior experience with the package | Guided experience that produces a visible result | The five-verb pipeline executed top to bottom; the η=0 oracle as a triumphant finish; copy-paste fixtures that just work; live `@example` output the reader can compare against their own run | Math derivations ("why prop_P works"); design decisions; API option tables; alternative fixtures; the word "note that…" followed by a theorem |
| **How-to guide** | Competent user, specific goal | Steps to accomplish a real task | "Certify the defect without MOSEK"; "Extract block sizes from a factorization"; "Interpret a CertifiedBracket"; "Check whether your channel is in factorize's domain"; "Install MOSEK and use the tight bracket"; "Build your own Kraus representation and check unitality" | Teaching what a CertifiedBracket is; explaining why arb is used; the full pipeline from scratch; reference tables |
| **Explanation** | Engaged, curious, away from keyboard | Understanding — the why, the how the pieces fit | The four-stage mathematical story (η-idempotence → ε-C*-algebra → th_main isomorphism → factorization); certified arithmetic and the cross-check ladder; the constructive mandate (why non-constructive proofs are replaced); dimension-independence and why the naive route fails; the C+FLINT+Julia architecture; the Choi-Effros product vs plain multiplication; cb-norm vs operator norm; known limits and their root causes | Any executable code steps; API signatures; install instructions |
| **Reference** | Working, knows what they want | Lookup — exact signature, semantics, constraints | Every exported name with its full docstring; type fields; accepted argument ranges; what each return type's fields mean; the three known limits as flat, scannable statements; the notation table with paper anchors | Narrative text; motivation; any sentence beginning "In order to understand…" |

**The one rule to enforce rigorously:** if a sentence explains *why* something
is true, it belongs in Explanation, not Tutorial or How-to. Tutorials may
say "this returns a rigorous bracket" but not "this is rigorous because arb
rounds error balls outward." The latter is one sentence of Explanation —
pull it out, link to it.

---

## B. Proposed site map / navigation

This is the concrete `pages = [...]` argument for `docs/make.jl`.
Annotations: (M) = must-have for a usable site; (N) = nice-to-have, add
in a second pass; (merge?) = note where two pages could combine if the
first pass feels too sparse.

```julia
pages = [
    "Home" => "index.md",                        # (M)

    "Getting started" => [                        # (M) — Diátaxis: Tutorial-adjacent
        "Installation" => "getting_started/install.md",   # (M)
        "Quick start"  => "getting_started/quickstart.md", # (M)
    ],

    "Tutorials" => [                              # (M) — pure Tutorial mode
        "The five-verb pipeline"         => "tutorials/pipeline.md",       # (M)
        "The η = 0 oracle"               => "tutorials/eta0_oracle.md",    # (M)
        "Multi-block channels"           => "tutorials/multiblock.md",     # (N) merge↑ if short
        "Tight brackets with MOSEK"      => "tutorials/mosek_tight.md",    # (N)
    ],

    "How-to guides" => [                          # (M) — pure How-to mode
        "Certify the defect (solver-free)"        => "howto/certify_defect.md",     # (M)
        "Interpret a CertifiedBracket"            => "howto/read_bracket.md",       # (M)
        "Extract block structure"                 => "howto/block_structure.md",    # (M)
        "Factorize a channel"                     => "howto/factorize.md",          # (M)
        "Check factorize's domain"                => "howto/domain_check.md",       # (M)
        "Check trace-preservation of decode"      => "howto/iscptp.md",             # (M)
        "Install and use MOSEK"                   => "howto/mosek_install.md",      # (M)
        "Build a UCPMap from Kraus operators"     => "howto/build_ucpmap.md",       # (N)
    ],

    "Explanation" => [                            # (M) — pure Explanation mode
        "The mathematics: four-stage story"       => "explanation/math_story.md",   # (M)
        "Certified arithmetic and the cross-check ladder" => "explanation/certified_arithmetic.md", # (M)
        "Dimension-independence"                  => "explanation/dim_independence.md", # (M)
        "Architecture: C + FLINT + Julia"         => "explanation/architecture.md", # (M)
        "The constructive mandate"                => "explanation/constructive.md",  # (N)
        "Design decisions and known limits"       => "explanation/design_limits.md", # (M)
    ],

    "Reference" => [                              # (M) — pure Reference mode
        "API" => "reference/api.md",              # (M) — the current api.md, promoted
        "Notation glossary"   => "reference/notation.md",   # (M)
        "Bibliography"        => "reference/bibliography.md", # (N) merge into notation if short
    ],
]
```

### Page-by-page purpose (one line each)

**Home** (`index.md`) — existing file, already strong; acts as orientation
hub and learning-path signpost; the four result-plots live here as the
"proof by picture" introduction. Retain. Do not duplicate the README
verbatim — pull from it, link back to it for install.

**Getting started / Installation** (`getting_started/install.md`) — exact
shell commands: apt dependencies, cmake build, `set_libaic_path!`, green
test run. No theory. Numbers: "~2m09s, 219 passes."

**Getting started / Quick start** (`getting_started/quickstart.md`) — a
single executed `@example` block: construct `phit_kraus(2, 0.1)`, run
`certified_defect`, show the bracket, done. Fits on one screen. The reader
does one thing and sees a real output. Links to Tutorial for next steps.

**Tutorial: The five-verb pipeline** (`tutorials/pipeline.md`) — the
existing `tutorial.md`, restructured so it has no explanation paragraphs.
Steps 1–5 of the existing file are already close to correct Tutorial mode;
strip or move the "why prop_P works" aside to Explanation.

**Tutorial: The η = 0 oracle** (`tutorials/eta0_oracle.md`) — could stay
as the last section of pipeline.md (merge?) OR become its own short page.
Standalone is better: it is the "ground truth" motif used repeatedly and
deserves its own URL to link back to from How-to guides and Explanation.

**Tutorial: Multi-block channels** (`tutorials/multiblock.md`) — walk the
oblique `bce_conj_kraus(4,2,0.02)` fixture from channel construction
through `factorize`, showing `algebra(F)` with `[2,2]` blocks. (N): merge
into pipeline.md if the Tutorial page becomes too thin.

**Tutorial: Tight brackets with MOSEK** (`tutorials/mosek_tight.md`) — (N)
for users with MOSEK installed. Not a blocker for the (M) set.

**How-to: Certify the defect (solver-free)** — 5-step recipe; no theory.
Cross-links to Explanation/certified_arithmetic.

**How-to: Interpret a CertifiedBracket** — what `in`, `width`, `midpoint`,
`lo`, `hi` mean operationally. The one factual statement ("the bracket
width ~ 2n by design") is a fact, not an explanation — it belongs here.
The *reason* the width is 2n lives in Explanation.

**How-to: Extract block structure** — how to get `blocks(algebra(F))`,
what the integers mean, how to use them downstream. Four steps, no prose.

**How-to: Factorize a channel** — the single `factorize(Φ)` call with
domain pre-check, the four accessors (`encode`, `decode`, `delups_defect`,
`upsdel_defect`), and what to do when factorize throws `ArgumentError`.

**How-to: Check factorize's domain** — operational: compute
`spectral_radius(Φ^2 - Φ)` (or `rho(Φ)`) and compare to the `< 0.10`
threshold. Links to Explanation/design_limits for the reason.

**How-to: Check trace-preservation of decode** — `iscptp(decode(F); atol)`
with the correct tolerance, and the sentence: "the rigorous certificate is
the round-trip bracket, not iscptp."

**How-to: Install and use MOSEK** — exactly the install steps, the two
functions that unlock, the expect-output. No theory.

**Explanation: The mathematics — four-stage story** — the arc from
η-idempotence → ε-C*-algebra (Choi-Effros product, the approximate axioms)
→ th_main rigidity → th_factorization, with paper section/line anchors.
This is the page the quantum-information reader wants before touching code.
No code blocks except formula illustrations.

**Explanation: Certified arithmetic and the cross-check ladder** — WHY arb,
HOW the bracket is constructed (outward rounding, the Gelfand-spectral-radius
guard), the four-rung ladder, and when each rung fires. The one page that
makes a numerical-analysis reader trust the output.

**Explanation: Dimension-independence** — the central claim of the paper
(tex:484), why it is surprising, why the naive Haar-averaging route fails,
and the measured evidence from the dim-sweep (OLS slope -2.7e-4). The
universality.png plot lives here.

**Explanation: Architecture: C + FLINT + Julia** — the three-layer stack,
why C+FLINT rather than pure Julia, the double/arb dual-path, the ccall
surface. For users who want to know what they are depending on.

**Explanation: The constructive mandate** — (N) why non-constructive proofs
(Lefschetz-Hopf, contour integrals, Haar measure) are replaced by
constructive algorithms, and what that means for correctness and reproducibility.
Worth having for the operator-algebra audience who read the paper.

**Explanation: Design decisions and known limits** — the three stated limits
(factorize domain, O(η)-TP for decode, eig-free bracket width), their root
causes, and any design decisions that deviate from naive transcription. The
"Quirks" page in SQLite idiom — lead with constraints, not capabilities.

**Reference: API** — the current `api.md`, promoted to `reference/api.md`.
Group by: types (`UCPMap`, `CPMap`, `CertifiedBracket`, `EpsCStarAlgebra`,
`CStarAlgebra`, `MainIsomorphism`, `ChannelFactorization`), then verbs
(`certified_defect`, `associated_algebra`, `main_isomorphism`, `factorize`,
`encode`, `decode`), then utilities. Explicit `@docs` blocks rather than
`@autodocs` — the order should be curated, not alphabetical. Keep the
`@index` for the flat lookup list.

**Reference: Notation glossary** — one table: symbol, meaning, paper anchor
(tex line), Julia type. Symbols: Φ, η, ε, A, B, ⊕, M_d, v, Υ, Δ, ‖·‖_cb,
CertifiedBracket, ρ. This table is the single most valuable cross-reference
device for a mixed-audience package.

**Reference: Bibliography** — (N) the Kitaev arXiv entry + 4-5 items
(Choi 1975, Effros-Kuan 1987, Watrous SDP, FLINT, arb). Could merge into
notation.md.

---

## C. Standout checklist

These are checkable, binary conditions. Each references the source of the
recommendation.

1. **Every page opens with a one-sentence statement of what the reader will
   be able to DO or UNDERSTAND after reading it.** (Diataxis rule: state the
   contract. Source: diataxis.fr/tutorials/) Tutorial pages: "do"; Explanation
   pages: "understand"; How-to pages: "accomplish X"; Reference: "look up X."

2. **A "start here" learning path on the Home page.** Two or three named
   routes: "I am a quantum-information researcher" → Explanation/math_story
   then Tutorial/pipeline; "I want to use the package now" → Quick start;
   "I want to understand the certified arithmetic" → Explanation/certified.
   (Source: JuMP.jl docs structure; Makie.jl progressive pathway.)

3. **No un-run code anywhere.** Every non-trivial code block is either a
   `jldoctest` (output checked exactly) or `@example` (output embedded live).
   `@setup` blocks hide boilerplate. Static code fences (```` ```julia ````)
   only for pseudocode or illustrative fragments that cannot be executed (e.g.
   shell commands). (Source: Documenter.jl manual; existing make.jl discipline.)

4. **A notation glossary with paper line-anchors.** One scannable table:
   symbol → meaning → `approximate_algebras.tex:NNN`. Without this, a reader
   who opens the docs after reading the arXiv preprint has no map between
   the two. (Source: domain-hallucination-risk section of CLAUDE.md; the
   paper's own notation density.)

5. **The three known limits as a named, findable section.** Not buried in
   tutorial footnotes — a flat page (Explanation/design_limits) and a short
   admonition on the Home page. Each limit states: what it is, the measured
   number ("atol=1e-4 not 1e-9"; "hi/lo~2n"; "ρ<0.10"), and a cross-link to
   the explanation. Model: SQLite "Quirks of SQLite" as a top-level item.
   (Source: sqlite.org architecture; CLAUDE.md Rule 12 "concrete numbers always".)

6. **Admonitions used for exactly three purposes — not decoration.** Use
   `!!! warning` for the two hard constraints (factorize domain; decode
   iscptp tolerance). Use `!!! note` for the one design choice that surprises
   (eig-free bracket loose by design). Use `!!! tip` for "if you have MOSEK."
   Never use `!!! danger` or `!!! note` for general elaboration — that belongs
   in prose or Explanation. (Source: JuMP.jl tutorial admonition audit;
   Documenter.jl syntax docs.)

7. **Explicit `@docs` blocks in curated order, not `@autodocs`.** The
   Reference page should group and order exports the way a reader thinks
   about them (types → verbs → utilities), not alphabetically. Keep `@index`
   for the flat lookup list. (Source: Documenter.jl guide; JuMP.jl API
   reference structure.)

8. **The four result-plots placed in Explanation, not buried in the home
   page only.** `containment.png` → Explanation/certified_arithmetic.
   `universality.png` → Explanation/dim_independence. `roundtrip.png` →
   Explanation/design_limits or Tutorial/pipeline. `eta0_oracle.png` →
   Tutorial/eta0_oracle. Each plot appears once, with a caption containing
   the concrete numbers from the README (OLS slope, measured max, etc.).
   The Home page may show thumbnails but links to the full plot + context.

9. **Cross-links are bidirectional: How-to ↔ Explanation ↔ Reference.**
   Each How-to page ends with: "For the reason this works, see Explanation/X.
   For the full signature, see Reference/api." Each Explanation page ends with:
   "To apply this, see How-to/Y." Each Reference docstring's "See also" links
   to the relevant Tutorial section. (Source: Diataxis rule on navigation;
   Documenter.jl `@ref` syntax.)

10. **The η = 0 oracle appears in Tutorial, Explanation, and as a cross-link
    in every How-to that produces a certified value.** It is the recurring
    "ground truth" motif: the only point where every defect is unambiguously
    correct. Every How-to that produces a bracket should say "on an exactly
    idempotent channel this bracket contains 0 — see Tutorial/eta0_oracle."
    This motif is also the discipline test for any future feature: if it does
    not give zero on η=0 input, something is wrong.

11. **Concrete measured numbers, never "fast", "accurate", "tight".** Every
    claim about performance or precision must be anchored to a specific test
    and a specific number: "219 tests in ~2m09s on a 2023 laptop"; "hi/lo~2n
    at prec=256"; "max primal-dual gap 1.23e-11 (MOSEK, test_sdp.jl)". The
    README already does this; the site must match the same discipline.
    (Source: CLAUDE.md Rule 12; SQLite documentation style.)

12. **A "How to trust the output" section = the cross-check ladder.** The
    four rungs (double vs arb; certified containment; η=0 oracle; strong
    duality) form the core of Explanation/certified_arithmetic. This is the
    single page that answers a numerical analyst's first question. It should
    be discoverable from the Home page "start here" block without three clicks.

13. **jldoctest blocks for stable structural output only; @example for floats.**
    Structural facts (block sizes, type names, dims, boolean flags) use
    `jldoctest` with exact-match output. Float-heavy results use `@example`
    so last-digit wobble never fails the build. `@setup` hides repeated fixture
    boilerplate. This discipline already exists in tutorial.md — enforce it
    uniformly across all new pages. (Source: existing make.jl comment; Documenter.jl doctest guide.)

14. **The Home page is NOT the Tutorial and NOT the README.** Home orients:
    one-paragraph statement of what the package does; the five-verb table;
    the "start here" routes; the Known limits admonition; the citation. The
    Tutorial teaches. The README is for GitHub. Home should link to both but
    duplicate neither. Currently index.md is close to this but contains
    Tutorial-mode content (the full pipeline @example blocks) that should
    migrate to Quick start or Tutorial/pipeline. (Source: Diataxis Home ≠
    Tutorial rule; README.md audit showing full pipeline already covered.)

15. **Named doctest modules share state within a Tutorial page, unnamed ones
    do not.** The existing convention (`tut`, `tut1`, `tut2`, `tut3`) is
    correct. New pages must follow: one `@example <pagename>` module per
    page for executed examples; separate `jldoctest <label>; setup = :(…)`
    for doctests needing prior definitions. Never use unnamed jldoctests where
    state is needed across blocks. (Source: Documenter.jl doctest guide
    on named vs unnamed modules.)

16. **A "You are here" sentence at the top of every Explanation page.**
    Because Explanation is read away from the product, the reader can land
    on it from a search engine with no context. One sentence locating the
    page in the Diataxis structure: "This page explains the mathematics
    behind the five-verb pipeline; for how to run it, see Tutorial/pipeline."
    (Source: Diataxis explanation mode rules; SQLite one-sentence openers.)

---

## D. Engagement devices appropriate for this package

These are structural decisions, not gimmicks. Each has a specific placement.

### The four result-plots

Each plot is a passing test drawn. The caption must contain the concrete
number from the test assertion. Placement:

- `containment.png` — Explanation/certified_arithmetic. Caption: "The eig-free
  certificate [lo, hi] brackets the analytic η = t(1−t)·2(1−1/d²) for all t;
  containment is a test invariant (test_core.jl)."
- `universality.png` — Explanation/dim_independence. Caption: "c = isodefect/η
  over dim A ∈ {8,12,16,18,20}; OLS slope -2.7e-4, max c = 0.047 (test_dbo3.c,
  prec=256). Tex:484 holds."
- `roundtrip.png` — Tutorial/pipeline (Step 4 section). Caption: "Certified
  round-trip defects ‖ΔΥ−Φ‖_cb and ‖ΥΔ−1_B‖_cb vs η; linear in η, below
  the 50·η test envelope (test_factorize.jl)."
- `eta0_oracle.png` — Tutorial/eta0_oracle. Caption: "On exactly-idempotent
  channels every certified defect collapses below 1e-9; measured maxima are
  4.4e-75 and 3.9e-75 (test_factorize.jl)."
- `demo.svg` — Home page, as the hero visual. No caption needed on Home;
  the five-verb table is the caption. The SVG is the "proof by animation"
  the README already uses.

### The "story arc" on the Home page

The Home page should open with a three-sentence narrative arc before any
code or tables: (1) the problem — a near-idempotent channel hides a clean
structure; (2) what Kitaev proved — rigidity, dimension-independent constant;
(3) what this package delivers — every non-constructive proof replaced by
a finite algorithm, every result certified by an arb error ball. The arc
is already in index.md/README.md — the task is to make it the first thing
the reader encounters, before the five-verb table, not after.

### The η = 0 oracle as recurring motif

Every tutorial section that produces a certified quantity ends with:
"On a channel that is exactly a conditional expectation (η = 0), this
quantity is 0. Run Tutorial/eta0_oracle to see all five collapse simultaneously."
This creates a recurring reference point that a reader can always return to
as a sanity check on their own channels.

### Worked end-to-end examples

The Tutorial pages are already nearly correct. The additional engagement
comes from:
- Showing the output of every `@example` block next to the code, so the
  reader sees the concrete `CertifiedBracket{...}` type repr rather than
  having to run it.
- Putting the `universality.png` plot in the Explanation with the dim-sweep
  table from the README (the plain-text benchmark table), so the claim is
  supported by both a figure and a number table.
- Including the `bce_conj_kraus(4,2,0.02)` oblique multi-block fixture in
  the Tutorial (it is already in tutorial.md Step 4) with a comment explaining
  WHY that fixture is more interesting than `phit_kraus`: it has a genuinely
  oblique image that is not a *-subalgebra on the nose — that is the generic
  case.

### Navigation signposts

Three explicit "you are here" routes on the Home page (see checklist item 2),
plus a "Where to go next" footer on every Tutorial page (already in tutorial.md
— standardize it). Every How-to page opens with a one-line assumption statement:
"This guide assumes you have run the Quick start and have a `UCPMap` in scope."

---

## E. Anti-patterns to avoid for this specific package

Listed by failure mode, with the specific Rule or finding it violates.

**1. Mixing Tutorial and Explanation on the same page.**
The existing `tutorial.md` already has one violation: the paragraph beginning
"associated_algebra builds A = Img Φ̃ with Φ̃ = θ(2Φ−1), the exact-idempotent
regularisation (prop_P, §3). A is an oblique subspace…" is Explanation, not
Tutorial. A Tutorial may name the function; it may not explain why the
construction is shaped as it is. Move the prop_P explanation to
Explanation/math_story. (Diataxis rule; identified during tutorial.md audit.)

**2. Duplicating the README.**
The index.md currently reproduces most of the README, including the full
pipeline `@example` blocks, the five-verb table, the Known limits section,
and the citation. On the site, the Home page and README serve different
audiences (site visitors vs GitHub README readers). The site Home should
orient and route; the README can stay as-is on GitHub. Duplication creates
maintenance debt: a number that changes (e.g. "219 passes") must be changed
in two places. Solution: the site Home links to the README for the detailed
pipeline demonstration and pulls only the three-sentence story arc and the
demo.svg. (CLAUDE.md Rule 12: concrete numbers must not drift.)

**3. Un-run code blocks.**
Any ```` ```julia ```` fence in a Tutorial or How-to page is a liability.
There are currently no such fences in tutorial.md (it is all `@example` or
`jldoctest`). New pages must maintain this discipline. If a code pattern
cannot be executed in the docs environment (e.g. it requires MOSEK), put it
in a `!!! tip "With MOSEK installed"` admonition with a static fence, and
say explicitly that the block is illustrative.

**4. Marketing language.**
Forbidden: "powerful", "robust", "blazing", "seamless", "production-grade",
"state of the art." Also forbidden: "Note that AlmostIdempotentChannels
is unique in providing…" CLAUDE.md Rule 12 is the authority. The SQLite
model: "Small. Fast. Reliable. Choose any three." — acknowledges tradeoffs
explicitly. For this package: "solver-free by default; loose bracket by design
(hi/lo ~ 2n); tight only with MOSEK."

**5. Math without a paper anchor.**
Any equation appearing in the docs (Explanation pages especially) must cite
the tex line: "(approximate_algebras.tex:347)" or "(arXiv:2405.02434, §3,
prop_P)". Without anchors, the docs accumulate unverifiable claims. This
is the same rule as CLAUDE.md Law 1 applied to documentation.

**6. Admonitions as emphasis.**
Do not use `!!! note` to highlight a paragraph that is just important prose.
That is `**bold**`. Admonitions are for three specific purposes (see checklist
item 6) — constraints, surprises, extensions. Overused admonitions are noise
and readers stop reading them.

**7. Alphabetical API ordering without grouping.**
`@autodocs` dumps exports alphabetically: `CertifiedBracket`, `CPMap`,
`CStarAlgebra`, `EpsCStarAlgebra`… A reader who wants to understand the
verbs must hunt. Use explicit `@docs` blocks grouped as: (a) the channel types,
(b) the five headline verbs, (c) the result types and their accessors,
(d) utilities. The `@index` block can remain alphabetical for lookup.

**8. A glossary that is prose, not a table.**
The notation is dense (Φ, η, ε, A, B, v, Υ, Δ, ‖·‖_cb, ρ, ⊕, M_d, Choi-
Effros product). A prose glossary requires the reader to read it linearly.
A two-column table (symbol | meaning + paper anchor) lets a reader scan in
two seconds. Build the table; not the essay.

**9. Treating `iscptp` return value as the rigorous certificate.**
The existing Known limits text already handles this correctly. The anti-pattern
to avoid in new How-to pages: saying "check with `iscptp(decode(F))`" without
immediately saying "at atol=1e-4, not the default 1e-9, and the rigorous
certificate is the round-trip bracket." This distinction is load-bearing
for anyone building on top of the package.

**10. `deploydocs` or any `.github/workflows/` file.**
CLAUDE.md Rule 13 is unconditional. The site must build with a single local
command:
```
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```
Do not add a deploy step even as a commented-out template; someone will
uncomment it.

---

## Implementation order

Suggested sequence so the site is usable at each step, not half-done:

1. (M, first) Create `docs/src/reference/` and move `api.md` there with
   curated `@docs` groups; update `make.jl` pages tree.
2. (M) Create `docs/src/reference/notation.md` — the notation table.
3. (M) Create `docs/src/getting_started/install.md` and `quickstart.md` —
   pull install steps from README; single-block quickstart.
4. (M) Refactor `docs/src/index.md` — strip Tutorial-mode content, add
   learning-path routes, keep demo.svg + story arc.
5. (M) Refactor `docs/src/tutorial.md` → split into `tutorials/pipeline.md`
   and `tutorials/eta0_oracle.md`; strip Explanation content into stubs.
6. (M) Create `docs/src/explanation/math_story.md` — the four-stage story
   with paper anchors; receive the content stripped from tutorial.md.
7. (M) Create `docs/src/explanation/certified_arithmetic.md` — the cross-check
   ladder with containment.png.
8. (M) Create `docs/src/explanation/dim_independence.md` — with universality.png
   and the dim-sweep number table.
9. (M) Create `docs/src/explanation/design_limits.md` — the three known limits.
10. (M) Create the seven How-to stubs with correct one-sentence openers and
    cross-links; flesh out in order of user frequency (certify_defect, read_bracket,
    factorize, domain_check, iscptp, block_structure, mosek_install).
11. (N) `explanation/architecture.md`, `explanation/constructive.md`,
    `tutorials/multiblock.md`, `tutorials/mosek_tight.md`, `howto/build_ucpmap.md`,
    `reference/bibliography.md`.

---

## Sources

- https://diataxis.fr/ — four-mode framework
- https://diataxis.fr/tutorials/ — Tutorial rules
- https://diataxis.fr/how-to-guides/ — How-to rules
- https://diataxis.fr/explanation/ — Explanation rules
- https://diataxis.fr/reference/ — Reference rules
- https://documenter.juliadocs.org/stable/man/guide/ — @ref, pages=, @contents
- https://documenter.juliadocs.org/stable/man/syntax/ — @example, @repl, @setup, @raw
- https://documenter.juliadocs.org/stable/man/doctests/ — jldoctest, named modules, filters
- https://jump.dev/JuMP.jl/stable/ — navigation structure, admonition usage, progressive scaffolding
- https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/ — tutorial structure
- https://docs.makie.org/stable/ — Diataxis-aligned navigation; backend-choice upfront; scaffolded pathway
- https://www.sqlite.org/docs.html — prose style; "Quirks" as top-level item; concrete numbers;
  "lead with constraints"; no marketing language
- Project files: docs/src/index.md, docs/src/tutorial.md, docs/src/api.md,
  docs/make.jl, README.md, CLAUDE.md
