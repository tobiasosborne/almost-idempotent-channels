# aic-dbo.3 — Dimension-sweep universality regression: design

> Status: DESIGN (no code). Bead `aic-dbo.3` — "Dimension-sweep regression:
> certified constants must NOT grow with dim (universality canary)". The
> highest-value correctness test surfaced by the research: it directly tests the
> paper's central universality claim, `approximate_algebras.tex:484` /
> `:461` — the implicit constant in `th_main`'s `O(eps)` is
> **dimension-independent**.
>
> This doc is a buildable blueprint. The implementing agent should read
> `CLAUDE.md` (Rules 6/7, the `tex:484` stop condition), the two existing
> canaries (`tests/test_cstar_build.c:353-511` T2b, `tests/test_errreduce.c:694-775`
> T4) and the API headers cited below before writing `tests/test_dbo3.c`.

---

## 0. The ground-truth claim under test

`approximate_algebras.tex:484` (verbatim, the failure mode we guard against):

> "For finite-dimensional `C*` algebras, a diagonal can be obtained as
> `D = ∫ dU (U† ⊗ U)` ... Unfortunately, naive constructions of the Haar
> measure (or just the diagonal) in the eps-associative setting have error
> bounds **proportional to `n = dim A`**. So the outlined procedure of fixing
> the multiplication works only if `eps < c n^{-1}`."

`approximate_algebras.tex:461` (the positive claim): the `th_main`
`O(eps)`-isomorphism has a **universal, dimension-independent constant**. The
paper avoids the `tex:484` trap by the *incremental* construction whose
`cor_improvement` error reduction (`tex:1317`) uses the **explicit
generalized-Pauli `B`-diagonal** (`tex:1249-1254`, `‖D‖ = 1` regardless of `n`),
NOT the eps-associative Haar diagonal.

The canary's job: build `A` over a **growing `dim_A`**, measure the realized
`O(eps)` constant `c` at each point, and assert `c` does **not** grow with
`dim_A`. A `c ∝ dim_A` (or even `∝ √dim_A`) is the `tex:484` failure mode and a
**stop condition** (escalate `aic-1bc`), not a test bug — see §7.

---

## 1. What aic-dbo.3 adds over the existing T2b canary

`tests/test_cstar_build.c:353-511` (T2b) already measures `c = iso_def/eta` and
runs the dual AND-gate predicate (§4) over a sweep. **But its fixture is
`make_mixconj(d, m, t)`** — a *single* compressed `M_m` block conjugated by a
fixed unitary (`test_idemp.h`; the docstring at `test_cstar_build.c:159-169`
confirms the base is a SINGLE `M_m` block, so every `eta>0` mixconj is **one
equivalence class**, `A ≅ M_m`). The sweep there varies the **ambient** `n`
(`d = 4..10`) at fixed block algebra `M_2` / `M_3`. That sweeps *fixture
geometry / compression overhead*, not `dim_A`: `dim_A = m²` is held at 4 or 9
across each sub-sweep. The T2b banner is explicit that the residual variation it
sees is "fixture-GEOMETRY SPREAD across different ambient n, NOT
dimension-GROWTH" (`test_cstar_build.c:355-365`).

`tests/test_errreduce.c:694-775` (T4) *does* grow `dim_A` (single-block
`B = M_2..M_5`, `dim_B = 4..25`, and a block-count sweep `m=1,2,3` of `M_2`) and
asserts the `c_0`-ratio bound — but it measures `c_0` **directly via
`aic_errreduce`** on a hand-`build_inclusion`'d `v`, i.e. it tests the
`cor_improvement` engine in isolation, NOT the full `aic_cstar_build` master
loop. And its A-family is again `make_mixconj`/`make_mixconj_bce` (compressed
blocks), with the block sizes encoding the algebra only obliquely.

**aic-dbo.3's new coverage, precisely:**

1. **The literal `tex:484` dimension-blowup family.** It uses the NEW
   equal-block generator `aic_adv_chan_blockalg(out, k, d, t, prec)` whose
   associated algebra is `A = ⊕_{j=1}^{k} M_d`, `dim_A = k·d²` — exactly the
   `domain.md:416-449` Family 3D construction `A = ⊕_{j=1}^{k} M_{d_j}` with
   `d_1=…=d_k=d`. This is the family the paper's warning is *about*; T2b's
   single-block `M_m` is not.

2. **Two independent dimension axes.** `dim_A = k·d²` grows two genuinely
   different ways:
   - **fix `d`, grow `k`** (more equal blocks): stresses Stage-3 **class
     merging** (`cstar_build` Stage 3, `tex:1443` — `k` equivalence classes
     merged via `cor_merge_sum`; the per-class block-diagonal `B`-diagonal is the
     embedded-sum diagonal whose `‖D‖_proj = m` worst-case bound is the
     `test_errreduce.c:750-752` T4(B) concern). A `c` that grew with **block
     count** would surface here.
   - **fix `k`, grow `d`** (bigger blocks): stresses Stage-2 **inductive matrix
     growth** (`cstar_build` Stage 2, `tex:1430-1441` — `M_{r-1} → M_r` via
     `lem_extension`, `d-1` extension steps per block) and the generalized-Pauli
     diagonal of `M_d` (`tex:1251`, the `d²` terms each `1/d²`). A `c` that grew
     with **block size** would surface here.

   T2b grows neither axis (it grows ambient `n` at fixed `M_m`); T4 grows them
   but only through `errreduce`, not the master loop. aic-dbo.3 grows **both**,
   through the **full `aic_cstar_build`** (Stage 1 skeleton + Stage 2 growth +
   Stage 3 merge), which is the faithful end-to-end `tex:484` test: every stage's
   error budget must stay dimension-flat for `c` to stay flat.

3. **End-to-end constant, not engine constant.** The headline metric is the
   `th_main` isomorphism defect `iso_def` from `aic_cstar_build`
   (`aic_cstar.h:489`), normalized to a constant. That folds in every stage's
   `cor_improvement` call; a dimension-leak in *any* stage (a Haar diagonal
   sneaking in, a `√d` Frobenius-to-operator inflation, a per-merge `m`-factor)
   shows up in the final `iso_def`. T4's per-call `c_0` cannot see a leak that
   only appears when stages compose.

**One-line summary for the implementer:** T2b = "constant flat as fixture
geometry varies at fixed `M_m`"; T4 = "`cor_improvement`'s `c_0` flat as `dim_B`
grows, engine in isolation"; **aic-dbo.3 = "the *end-to-end* `th_main` constant
flat as `dim_A` grows two ways on the literal `⊕ M_d` blowup family."**

---

## 2. The dimension sweep (compute-bounded for the laptop arb path)

### 2.1 The compute wall (read first)

`aic_cstar_build` runs on the parent eps-`C*` algebra `A` built by
`aic_assoc_ecstar_from_phi` (`aic_assoc.h:177`), which first forms the
**`N² × N²` superoperator** `S_Φ` (`N = k·d` the ambient Hilbert dimension),
regularizes it via `theta(2S-1)` (a Newton-Schulz `sgn` on an `N²×N²` non-normal
matrix), and SVD-extracts a `dim_A`-dim basis. Every Stage-1/2/3 step then does
corner compressions, `dhom` defect sweeps over `dim_A²` matrix-unit pairs, and
`errreduce` Newton iterations — all in arb at `CPREC`. **Cost is dominated by the
`N²×N²` arb linear algebra, growing like `(N²)³ = N⁶` per dense op.** This is why
the HANDOFF caps laptop `dim_A` at ~64 (`domain.md:448` lists `256` in the
catalogue, but that is out of laptop reach in arb).

The binding constraint is therefore **`N = k·d`, not `dim_A` directly.** At
`N=16` the superoperator is `256 × 256`; an arb SVD / Newton-Schulz at
`CPREC=256` on a `256×256` non-normal matrix, iterated across the master loop, is
the edge of the 600 s slow `TIMEOUT`. **Keep `N ≤ 16`** for every sweep point and
estimate-then-measure (§7).

### 2.2 Proposed schedule

Hold `t` fixed across the whole sweep (one `t` per axis; see §3.3) so `eta` — and
hence `eps` — is roughly constant and the **ratio `c` is the clean signal**.

**Axis A — fix `d=2`, grow `k` (block-count axis, Stage-3 stress):**

| point | `k` | `d` | `N=k·d` | `dim_A=k·d²` | superop size |
|------:|----:|----:|--------:|-------------:|-------------:|
| A0    | 2   | 2   | 4       | 8            | 16×16        |
| A1    | 4   | 2   | 8       | 16           | 64×64        |
| A2    | 8   | 2   | 16      | 32           | 256×256      |

`dim_A ∈ {8, 16, 32}` over a 4× range, `N ≤ 16`.

**Axis B — fix `k=1`, grow `d` (block-size axis, Stage-2 stress):**

| point | `k` | `d` | `N=k·d` | `dim_A=k·d²` | superop size |
|------:|----:|----:|--------:|-------------:|-------------:|
| B0    | 1   | 2   | 2       | 4            | 4×4          |
| B1    | 1   | 3   | 3       | 9            | 9×9          |
| B2    | 1   | 4   | 4       | 16           | 16×16        |
| B3    | 1   | 5   | 5       | 25           | 25×25        |

`dim_A ∈ {4, 9, 16, 25}` over a 6.25× range, `N ≤ 5` (cheap — the `√d` /
generalized-Pauli `d²`-term inflation lives entirely here, so push `d` as far as
the budget allows; `d=5` matches T4(A)'s widest point).

**Optional bridge point (if budget allows after measuring §7):** `(k=4, d=4)`
gives `N=16`, `dim_A=64` — the catalogue's mid target and the only point reaching
`dim_A=64`. It is a `256×256` superop AND a 4-class × 4-step master loop, the
single most expensive instance. **Treat it as opt-in:** add it behind a measured
wall-time check; if a solo run of `(8,2)` (A2) already approaches ~300 s, **do
NOT add `(4,4)`** and LOG the cap explicitly (see §2.4).

### 2.3 Recommendation: run BOTH axes

Run **both** Axis A and Axis B. They are independent dimension axes and catch
different leaks (block-count vs block-size); two clean axes is materially stronger
evidence of universality than one. They also share no expensive point except
possibly the optional bridge, so the combined cost is ≈ `cost(A) + cost(B)` and B
is nearly free (`N ≤ 5`). Run each axis through the **same** dual-gate predicate
(§4) as an independent sub-sweep, exactly as T2b runs `m=2` and `m=3` separately
(`test_cstar_build.c:440-466`).

### 2.4 Precision and caps

- **`CPREC = 256`.** Match T2b/T1 (`test_cstar_build.c:69 #define CPREC 256`).
  The T2b banner records that `prec 256 == prec 512 byte-for-byte` on this
  family's `c` (`test_cstar_build.c:357-358`), so 256 is sufficient and stable; do
  NOT drop below — the `iso_def`/`eps` ratio at small `eps` is a difference of
  `O(1)` quantities (CLAUDE.md "Arbitrary precision — why", case (i)).
- **No silent caps (CLAUDE.md probe hygiene).** If the budget forces dropping A2
  (`N=16`) or the `(4,4)` bridge, the test MUST `printf` the cap and the reason
  ("dim_A=64 dropped: (8,2) solo wall-time X s, (4,4) projected > TIMEOUT"), and a
  follow-up bead filed. A dropped point is a logged decision, never an
  unremarked omission.
- **Per-axis budget gate.** Wrap each sweep point's build so a non-converging
  point fails loud under the 600 s `TIMEOUT` rather than hanging (it will: every
  arb routine here is bounded, and `aic_cstar_build` asserts Stage-1 termination
  in `≤ dim_A+1` iterations, `aic_cstar.h:477`). No `sleep`, no unbounded loop.

---

## 3. Which constant(s) to measure

### 3.1 Headline metric — `c_thmain` (the `th_main` `O(eps)` constant)

At each sweep point, build `A` and run the master loop, then form the ratio:

```
aic_adv_chan_blockalg(&phi, k, d, t, CPREC);          // A = ⊕_j M_d
eta = eta_proxy(&phi, CPREC);                          // ‖S_Φ²−S_Φ‖_op midpoint route
aic_assoc_ecstar_from_phi(&ae, &phi, CPREC);           // A as aic_ecstar
aic_ecstar_defect_assoc(eps_ball, &ae.A, CPREC);       // eps = ax_assoc defect
eps = dd(eps_ball);
aic_cstar_build(&B, &v, iso, &ae.A, eps, CPREC);       // th_main isomorphism
c_thmain = dd(iso) / eps;                              // the O(eps) constant
```

This is the metric the bead names first ("the O(eps)/O(eta) constants in
`th_main`"). `iso_def` is the master loop's certified isomorphism defect
(`aic_cstar.h:470-472`, "the max of the four `aic_errreduce` inclusion-defects of
the FINAL `v_out`"); `iso_def/eps` is the realized `O(eps)` constant. A true
universal isomorphism has `c_thmain` bounded and flat across the sweep.

**Normalize by `eps`, not `eta`** here (differs from T2b, which normalizes by
`eta`). Rationale: `aic_cstar_build`'s `iso_def ≤ AIC_ERRREDUCE_C0_CERT·eps'`
(`aic_cstar.h:472`), so `iso_def/eps` is the directly-comparable `O(eps)`
constant the build itself certifies against. Use `eps = aic_ecstar_defect_assoc`
(the associativity defect, `aic_ecstar.h:136`), the same `eps` passed INTO the
build. (T2b normalizes by `eta` because it predates this convention and because on
its single-class fixture `eps_assoc ≈ 700×` below `eta`,
`test_cstar_build.c:164-168`; for the `⊕ M_d` family measure and report BOTH
`eta` and `eps` per point in the `printf` so the implementer can pin the cleaner
normalizer empirically — but gate on `iso_def/eps`.)

**Why the ratio is the right signal (and is `eps`-scale-free):** even if the
absolute `eps` drifts with `dim_A` (it may — `ax_assoc` is a basis sweep whose
magnitude can track `dim_A`), `c = iso_def/eps` divides that drift out. `c` stays
flat **iff** the *construction* is dimension-universal. A construction that leaked
a `dim_A` factor into `iso_def` would push `c` up even as `eps` rises — exactly
what we want to catch.

### 3.2 Optional secondary metric — `c_0` (`cor_improvement`)

The bead also names "`c_0` from `cor_improvement`." Two routes:

- **(a) via `aic_cstar_build` internals.** `aic_cstar_build` calls
  `aic_errreduce` at every stage but does **not** expose the per-stage `c_0` —
  only the final `iso_def` (`aic_cstar.h:470`). So `c_0` is *not* cleanly
  extractable from the public master-loop output. `iso_def` already folds the
  worst stage's `cor_improvement` defect in, so **`c_thmain` is the headline and
  subsumes `c_0` for the end-to-end test.**

- **(b) via the standalone `errreduce` path.** Measuring `c_0` directly means
  hand-building a `v : B → A` on the `⊕ M_d` family and calling `aic_errreduce`
  with `c0_out` (`aic_errreduce.h:117`), exactly the `sb_errreduce` /
  `mb_errreduce` helpers of `test_errreduce.c:275-319`. **This is already covered
  by `test_errreduce.c` T4** (`test_errreduce.c:694-775`), which sweeps `dim_B`
  and asserts the `c_0`-ratio < 1.6 with a `√dim`-injection mutation tooth
  (`test_errreduce.c:705-718`).

**Recommendation:** make `c_thmain` (route 3.1) the **headline and sole gated
metric** of `tests/test_dbo3.c`. Do NOT duplicate the `c_0` engine sweep — note
in the test docstring that the `cor_improvement` `c_0` axis is covered by
`test_errreduce.c` T4, and that aic-dbo.3's distinct contribution is the
**end-to-end `th_main` constant on the literal `⊕ M_d` blowup family** (§1). If
the implementer wants belt-and-braces `c_0` coverage on the new family
specifically, add it as a *reportable* (non-gated) `printf` via the
`mb_errreduce`-style helper, but the gate is `c_thmain`.

### 3.3 Holding `t` (hence `eta`, hence `eps`) fixed

Pick one `t` per axis and hold it across all points of that axis so the ratio `c`
is the clean signal. The generator is declared `t=0` EXACTLY idempotent, `eta`
tunable by `t`; pick a `t` giving `eta ~ 1e-2` (well inside the `ρ(Φ²−Φ) < 1/4`
basin, CLAUDE.md spectral-basin callout) and **non-degenerate** (`eta > 0`,
`eps > 0`, `iso_def > 0` — the non-vacuity check, §5). Start from the T2b values
(`t=0.03` for the `d=2` family, `t=0.02` for larger blocks,
`test_cstar_build.c:399-402`) and adjust if `eta` drifts across the sweep:
**`eta` staying roughly constant as `dim_A` grows is itself worth a `printf`** —
if `eta` (or `eps`) blows up with `dim_A` at fixed `t`, that is a finding about
the generator (flag it; it may mean `t` must be re-tuned per point, in which case
hold `eta` constant by tuning `t`, not `t` constant). The ratio `c` divides `eps`
out either way, but a constant `eps` keeps the mutation tooth (§5) interpretable.

---

## 4. The pass predicate (reuse the T2b dual AND-gate)

Reuse the T2b predicate model verbatim (`test_cstar_build.c:386-466`). For each
sub-sweep (Axis A, Axis B independently), with the per-point ratios
`c[i]` and ambient dimension proxy `dimv[i]` (use `dim_A` as the x-axis here, not
ambient `n` — `dim_A` IS the quantity `tex:484` warns scales the error):

**(i) ABSOLUTE boundedness:** `c_abs_max < C_ABS`. A `c` growing with `dim_A`
eventually exceeds any fixed ceiling.

**(ii) NO GENUINE upward trend with `dim_A`:**
`NOT( halves_ratio(c) > KAPPA AND slope(c vs dim_A) > SLOPE_MIN )`. A violation
requires BOTH the upper/lower-half mean ratio AND the OLS slope to fire — robust
to a lone hard-geometry outlier (which moves the endpoint-sensitive ratio but not
the slope), with teeth against a genuine `c = O(dim_A)` law (which moves both).

**Starting constants (reuse T2b's, `test_cstar_build.c:388-397`):**

```
const double C_ABS     = 5.0;   // absolute boundedness ceiling
const double KAPPA      = 1.25;  // halves-ratio gate (endpoint-sensitive)
const double SLOPE_MIN  = 0.28;  // OLS-slope gate
```

**Re-tuning note (mandatory).** These were pinned to `make_mixconj`'s measured
`c` distribution (peak 3.27, `test_cstar_build.c:359-362`). The `⊕ M_d` family's
`c` distribution is **unmeasured** and may differ — possibly *flatter* (the family
is the clean blowup instance, not a hard-compression geometry) or possibly
*peakier* (multi-class Stage-3 merges are exercised for the first time at `eta>0`,
the I5 coverage path of `test_cstar_build.c:161-168`). The implementer **measures
the real `c` over the sweep first**, then pins:
- `C_ABS` to ~1.5× the measured peak (T2b uses 1.53× margin),
- `SLOPE_MIN` to sit between the measured geometry-noise slope and the
  mutation-tooth's faithful-`O(dim)` slope (`cflat/dim_min`, §5),
- `KAPPA` left at 1.25 unless a lone outlier on this family inflates it.

Record the measured `c` table and the pinned constants in a banner comment (T2b
style) and in `paper/FINDINGS.md` (a new §, cited from the test). The slope's
x-axis is `dim_A`, so `SLOPE_MIN` scales as `cflat / dim_A_min`; for Axis A
(`dim_A_min = 8`) a faithful `c·(dim_A/8)` law has per-step slope
`cflat/8 ≈ 0.12·cflat`; for Axis B (`dim_A_min = 4`) it is `cflat/4`. **Tune
`SLOPE_MIN` per axis** (the two axes have different `dim_A_min`), or normalize the
slope by `dim_A_min` so one constant serves both — recommend the latter
(divide the OLS slope by `1/dim_A_min`, i.e. regress `c` against
`dim_A/dim_A_min`, giving a dimensionless per-doubling slope comparable across
axes and to T2b).

**Helpers — extract to a shared header.** `trend_slope` and
`trend_halves_ratio` are currently `static` in `test_cstar_build.c:124-155`. Both
T2b and aic-dbo.3 (and arguably a future ext canary) need them. **Recommendation:
extract them into a new shared `tests/aic_trend.h`** (header-only, `static
inline`), include it from both `test_cstar_build.c` and `test_dbo3.c`, and delete
the `static` copies from `test_cstar_build.c`. This is a small refactor (~35 LOC
moved) that avoids a copy-paste fork of two load-bearing statistical helpers.
*Fallback* if the implementer wants to avoid touching the green T2b test: copy
the two functions into `test_dbo3.c` and add a one-line comment cross-referencing
`test_cstar_build.c:124-155` as the source of truth. Prefer the extraction;
duplicated statistical helpers drift.

---

## 5. The mutation tooth (Rule 7 — load-bearing)

Mirror T2b's two-model mutation (`test_cstar_build.c:468-511`). After the real
sweep, re-run the predicate on **injected** data — NO injection touches the build
path; the tooth proves the gate *would* fire on the `tex:484` failure mode.

**Model A (literal `tex:484` scaling):** `c_injA[i] = c[i]·(dim_A[i]/dim_A_min)`.
Amplifies the real (possibly noisy) `c` by the dimension ratio. Assert the
**ABSOLUTE arm fires**: `max(c_injA) ≥ C_ABS`. (On the `⊕ M_d` family, if the
real `c` is flat and monotone-ish this may *also* trip the trend arm; that is
fine — assert at least the absolute arm, the guaranteed catch.)

**Model B (faithful `c = O(dim_A)` law):**
`c_injB[i] = cflat·(dim_A[i]/dim_A_min)`, `cflat = mean(c)`. A clean monotone
dimension-proportional law — exactly the `tex:484` constant-scaling-with-`dim`
failure. Assert the **TREND arm fires**: `halves_ratio(c_injB) > KAPPA AND
slope(c_injB) > SLOPE_MIN`. This is the arm that exists *specifically* to catch a
constant that scales with dimension.

Together: the canary catches the literal injection (absolute arm) AND a genuine
`c ∝ dim_A` law (trend arm). Run the mutation on **both** axes' data (or at least
the axis with the wider `dim_A` range — Axis B, 6.25× — gives the trend arm the
most leverage). Mirror `test_cstar_build.c:483-511` structurally.

**Second-order NON-VACUITY check (Rule 5 — "runs without errors" is not a test).**
Before trusting any green, assert the test is *non-degenerate*: at the **smallest**
sweep point, `eps > 0` AND `iso_def > 0` (hence `c > 0`). A degenerate instance
(`eta` accidentally 0, or `aic_cstar_build` returning a trivial `iso_def=0`) would
make the dual gate pass vacuously. Add:

```
AIC_CHECK_MSG(eps_min > 1e-12 && c_min > 0.0,
    "dbo3: smallest-dim point is degenerate (eps=%.3e c=%.3e) — the canary "
    "would pass VACUOUSLY; the generator t must give eta>0, eps>0, iso_def>0",
    eps_min, c_min);
```

This guards against the recurring "test that can't fail" failure mode CLAUDE.md
calls out. (The generator's own self-test in `test_adversarial.c` should already
assert `eta(t=0)=0` and `eta(t)>0`; this is the *consumer-side* non-vacuity guard
for the canary specifically.)

---

## 6. File / LOC plan

### 6.1 New file: `tests/test_dbo3.c` (recommended over extending `test_adversarial.c`)

**Justification.** `test_adversarial.c` is the **corpus self-test** — it asserts
each generator's *knob → property* contract (`test_adversarial.c:114-921`), is
already in the slow set, and its concern is "does the generator produce the
claimed adversarial matrix/channel." aic-dbo.3 is a different concern: "does
`aic_cstar_build` keep its constant dimension-flat **across** the corpus's 3D
family." That is a *downstream-routine* regression, not a generator self-test, and
it pulls in a much heavier include set (`aic_cstar.h`, `aic_assoc.h`,
`aic_ecstar.h`, `aic_dhom.h`, `aic_errreduce.h`) that `test_adversarial.c` does
not currently carry. Keeping it separate keeps each test's concern and link
surface clean, and matches the existing split (T2b lives in `test_cstar_build.c`,
T4 in `test_errreduce.c`, each next to the routine it canaries — aic-dbo.3 is the
cross-family canary, so it gets its own file).

### 6.2 Structure of `tests/test_dbo3.c` (≤ 200 LOC, CLAUDE.md Rule 10)

```
/* docstring: realizes aic-dbo.3; the tex:484/461 universality canary on the
 * ⊕_j M_d blowup family (domain.md:416-449); how it differs from T2b/T4 (§1);
 * the dual-gate predicate; the mutation tooth; cites FINDINGS §<new>. */

#include <math.h>, <stdio.h>
#include <flint/acb_mat.h>, <flint/arb.h>
#include "aic_adversarial.h"      // aic_adv_chan_blockalg
#include "aic/aic_assoc.h"        // aic_assoc_ecstar_from_phi
#include "aic/aic_cstar.h"        // aic_cstar_build
#include "aic/aic_ecstar.h"       // aic_ecstar_defect_assoc
#include "aic/aic_dhom.h"         // aic_dhom_B / aic_dhom_v (build outputs)
#include "aic/aic_mat.h"          // aic_mat_opnorm (eta_proxy)
#include "aic_test.h"             // AIC_CHECK_MSG
#include "aic_trend.h"            // trend_slope, trend_halves_ratio (extracted §4)

#define CPREC 256

static double dd(const arb_t x);                    // arb midpoint -> double
static double eta_proxy(const aic_ucp_kraus*, slong);// ‖S²−S‖_op midpoint (copy
                                                     // test_cstar_build.c:84-105)

/* measure c_thmain at one (k,d,t): build A, run master loop, return iso_def/eps;
 * write eta, eps, dim_A out via pointers for the printf + non-vacuity. */
static double measure_c(slong k, slong d, double t, double *eta_out,
                        double *eps_out, slong *dimA_out);

/* run one axis sub-sweep: fill nv[]/cv[], apply dual gate, run mutation tooth. */
static void run_axis(const char *name, const slong *kk, const slong *dd_,
                     const double *tt, int npts);

static void test_axis_A_blockcount(void);  // d=2, k in {2,4,8}
static void test_axis_B_blocksize(void);   // k=1, d in {2,3,4,5}

int main(void) {
    test_axis_B_blocksize();   // cheap (N<=5) — run first as smoke
    test_axis_A_blockcount();  // heavier (N up to 16)
    printf("test_dbo3: PASS\n");
    return 0;
}
```

LOC budget: the two `eta_proxy`/`dd` helpers ~25 (or share via `aic_trend.h`
companion), `measure_c` ~30, `run_axis` (gate + mutation) ~70, two axis wrappers
~15 each, `main` ~10, docstring ~30. Comfortably under 200. If it crowds 200,
split the mutation tooth into `aic_trend.h` as a reusable
`trend_mutation_check(...)` (it is identical logic to T2b's, so extracting it
deduplicates T2b too).

### 6.3 CMake wiring (`tests/CMakeLists.txt`)

`test_*.c` auto-globs into a per-test executable (`tests/CMakeLists.txt:49-61`),
so `tests/test_dbo3.c` is picked up automatically — **no glob edit needed**. The
**only** required change: add `test_dbo3` to the `AIC_SLOW_TESTS` list
(`tests/CMakeLists.txt:41-47`), so it gets `TIMEOUT 600` + `LABEL slow` instead of
the fast `TIMEOUT 120` (`tests/CMakeLists.txt:70-73`). It belongs in the slow set:
it runs `aic_cstar_build` in arb at `N` up to 16, far over the sub-second fast
budget. Add it next to the other master-loop drivers:

```cmake
set(AIC_SLOW_TESTS
    test_cstar_build test_cstar_extension test_cstar_merging test_cstar_merge
    test_cstar test_certify test_certify_teeth test_cbnorm test_opspace
    test_opspace_o2 test_opspace_choi test_factorize test_assoc test_assoc2
    test_adversarial test_ucp_d24 test_dbo3            # <-- add here
    test_errreduce test_dhom test_projection test_corner test_ecstar)
```

If `aic_trend.h` is added (§4): it is a header (`aic_*.h`), NOT a `test_*.c` /
`aic_*.c`, so it is NOT globbed into a driver or the support lib — it is just
`#include`d. No CMake change for the header beyond the existing
`-Itests` (`tests/CMakeLists.txt:57`).

---

## 7. Risks / open questions for the implementer

1. **The compute wall at `N ≥ 16` (the binding risk).** A2 (`k=8,d=2`, `N=16`,
   `256×256` superop, 8-class Stage-3 merge) and the optional `(4,4)` bridge are
   the cost peaks. **Estimate first:** measure a solo run of `(8,2)` before
   committing the sweep; if it approaches ~300 s, drop A2 to `(6,2)` (`N=12`,
   `dim_A=24`) or cap Axis A at `(4,2)` and lean on Axis B for the wide range —
   and LOG the cap (§2.4). The `(4,4)` bridge (`dim_A=64`) is opt-in only.

2. **The `⊕ M_d` family's `c` may NOT be as flat as `make_mixconj`'s — and that
   could be a REAL finding.** This is the central open question. T2b's `c` is
   bounded in `[0.25, 3.27]` and flat *on its single-class family*. The
   equal-block family exercises, for the first time at `eta>0`, the **multi-class
   Stage-3 merge** (the I5 path, `test_cstar_build.c:161-168` notes this is the
   *only* fixture hitting it) and the **embedded-sum diagonal** whose worst-case
   `‖D‖_proj = m` block-count factor is exactly the T4(B) concern
   (`test_errreduce.c:750-752`). **If the measured `c` genuinely climbs with `k`
   (block count) or `d` (block size), the canary SHOULD go RED — that is a
   `tex:484` stop condition (escalate `aic-1bc`), NOT a test bug.** The
   implementer must be ready to distinguish: (a) a lone hard-geometry spike (dual
   gate's robustness absorbs it — fine), vs (b) a systematic monotone climb (both
   gate arms fire — escalate). The mutation tooth (§5) is what proves the gate can
   tell them apart. **Flag in the test banner and FINDINGS that a RED here is a
   candidate real finding, to be triaged before being treated as a regression.**

3. **Precision sufficiency at larger `N`.** `CPREC=256` is verified stable on the
   single-class family (`test_cstar_build.c:357`). On the multi-class `⊕ M_d`
   family at `N=16` the Newton-Schulz `sgn` on a larger non-normal superop, and
   the deeper master-loop recursion, *could* need more precision (the arb ball on
   `iso_def` losing all precision is a loud failure — CLAUDE.md Rule 4 / arb
   ladder). **If `iso_def`'s ball is not tight (`dd` midpoint unreliable), bump
   `CPREC` and re-measure;** if no feasible `CPREC` certifies, that is the
   precision-wall stop condition (CLAUDE.md stop conditions). Report the
   `arb_rel_accuracy_bits(iso)` at the worst point in a `printf` so precision loss
   is visible.

4. **Does `eta` stay controllable across the sweep?** At fixed `t`, `eta` may
   drift with `dim_A` (more blocks / bigger blocks → different
   `‖S_Φ²−S_Φ‖_cb`). The ratio `c` divides `eps` out so the gate is robust to
   this, but: (a) the basin guard `ρ(Φ²−Φ) < 1/4` must hold at *every* point or
   `aic_assoc_regularize` aborts (`aic_assoc.h:103-108`) — verify `eta` stays
   small enough; (b) if `eta` drifts a lot, the mutation tooth's interpretation
   blurs. **Recommendation:** `printf` `eta` and `eps` at every point; if `eta`
   drifts > ~2× across an axis, re-tune `t` per point to *hold `eta` constant*
   (tune `t`, not hold `t`), which keeps `c` the cleanest possible signal. The
   generator's `t=0 ⇒ eta=0` contract guarantees a continuous knob to do this.

5. **`aic_adv_chan_blockalg` is being built concurrently — depend on the FIXED
   signature only.** `void aic_adv_chan_blockalg(aic_ucp_kraus *out, slong k,
   slong d, double t, slong prec)`, declared in `tests/aic_adversarial.h`,
   `t=0` EXACTLY idempotent, `A = ⊕_{j=1}^k M_d`, `dim_A = k·d²`. Do not assume
   anything about its Kraus structure beyond: it returns a self-map on `B(C^{k·d})`
   whose `aic_assoc_ecstar_from_phi` yields `dim_A = k·d²` with `k` equivalence
   classes of `M_d`. **Add a build-side sanity assert** at the first sweep point:
   `B.num_blocks == k` and every `B.d[l] == d` and `B.dim_B == k·d²` (the
   `aic_dhom_B` fields, `aic_dhom.h:104-110`) — this both documents the contract
   and catches a generator that silently produces the wrong algebra (which would
   make the canary measure the wrong thing). If the generator's self-test in
   `test_adversarial.c` already asserts this, keep the consumer-side assert anyway
   (cheap, and it is the canary's own correctness precondition).

---

## Appendix — API surface the implementer relies on (verified against headers)

| Call | Header:line | Note |
|---|---|---|
| `aic_adv_chan_blockalg(out,k,d,t,prec)` | `tests/aic_adversarial.h` (concurrent) | `A=⊕_{1}^{k} M_d`, `dim_A=k·d²`, `t=0` exact-idemp; `out` init'd by generator, caller `aic_ucp_kraus_clear`s |
| `aic_assoc_ecstar_from_phi(&ae,&phi,prec)` | `aic_assoc.h:177` | builds `S_Φ` (`N²×N²`), regularizes, extracts basis; `aic_assoc_ecstar_clear` |
| `aic_ecstar_defect_assoc(eps_ball,&ae.A,prec)` | `aic_ecstar.h:136` | the `eps` (associativity defect), `arb_t` out |
| `aic_cstar_build(&B,&v,iso,&ae.A,eps,prec)` | `aic_cstar.h:489` | master loop; `iso` = isomorphism defect, `eps` is `double` |
| `aic_dhom_B` fields `num_blocks,d[],dim_B,n_B` | `aic_dhom.h:104-110` | block structure of built `B` (contract assert) |
| `aic_errreduce(&vt,c0,&defs,&v,eps,tol,steps,prec)` | `aic_errreduce.h:117` | only if doing the optional §3.2(b) `c_0` route (covered by T4) |
| `trend_slope`, `trend_halves_ratio` | `test_cstar_build.c:124-155` | extract to `tests/aic_trend.h` (§4) |
| `eta_proxy`, `dd` | `test_cstar_build.c:71,84-105` | copy/share; `eta_proxy` = `‖S²−S‖_op` midpoint route |
| dual-gate predicate + mutation tooth | `test_cstar_build.c:386-511` | the model to reuse (§4, §5) |
| `C_ABS=5.0, KAPPA=1.25, SLOPE_MIN=0.28` | `test_cstar_build.c:388-390` | starting constants; re-pin per measured `⊕ M_d` `c` (§4) |
| slow-set registration | `tests/CMakeLists.txt:41-47` | add `test_dbo3` |
