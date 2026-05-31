# Applications & API-friction report — driving AIC as a research instrument

> **What this is.** An experience report, not a design doc. It records the
> friction encountered when pointing the *existing green pipeline* at concrete
> quantum-information / many-body research questions, plus a prioritized feature
> wishlist that would turn "the algorithms are built" into "a researcher can run
> a study." Grounded in one end-to-end spike (certified cb-mixing) and a survey
> of ten candidate applications.
>
> The wishlist items are written so they can be filed as beads (prefix `aic`);
> none is research-hard except W7 (noted). House style: concrete numbers, cite
> the module / `.tex` line, no marketing.

## 1. The validating spike (so the frictions are grounded, not hypothetical)

**Idea: certified cb-mixing via almost-idempotence of channel powers.** For a UCP
self-map `Phi` with a unique fixed-point algebra, `Phi^n -> Phi^inf` (a
conditional expectation, hence idempotent), so the idempotence defect
`eta(Phi^n) = ||Phi^{2n} - Phi^n||_cb` decays to 0 at a rate set by the
subdominant eigenvalue: `eta(Phi^n) ~ |lambda_2|^n`. Reading the certified
bracket on `eta(Phi^n)` is therefore a *rigorous* cb/diamond mixing-rate
estimate (stronger than a single-input trace-norm bound).

Tested on qubit/qudit Pauli channels (where `lambda_2` is known analytically and
the self-dual superoperator is Hermitian, so the project's degeneracy-robust eig
confirms it). Modules exercised: `aic_ucp_kraus`, `aic_ucp_compose`,
`aic_cbnorm_eigfree_ball`, `aic_latd_eig_hermitian`.

Result (eig-free certified bracket, pure C, no MOSEK):

| channel | analytic \|lambda_2\| | fitted decay rate exp(slope) | rel. err |
|---|---|---|---|
| dephasing p=0.35 | 0.300 | 0.3040 | 1.3% |
| dephasing p=0.20 | 0.600 | 0.6275 | 4.6% |
| dephasing p=0.10 | 0.800 | 0.877 | 9.6% (pre-asymptotic at n<=14) |
| generic Pauli (0.55,0.20,0.10,0.15) | 0.500 | plateau `mid/\|l2\|^n` -> 1.26 | — |

Superoperator-eig cross-check (dephasing p=0.20, Hermitian): eigenvalues
`{0.600, 0.600, 1.000, 1.000}` exactly — `|lambda_2| = 0.6` confirmed
independently. The `mid/|lambda_2|^n` ratio plateaus to a clean constant
(~1.768 for dephasing, ~1.26 for the biased Pauli), confirming
`eta(Phi^n) = C |lambda_2|^n`.

**Key methodological finding.** The eig-free bracket has a *fixed* width factor
`2*dim` (measured `hi/lo = 4.00` exactly for the qubit). That factor is
**irrelevant for the rate** — both endpoints scale as `|lambda_2|^n`, so the
slope (hence `lambda_2`) is recovered regardless — but it makes the bracket loose
for an **absolute** `eta` value. Conclusion: *mixing-rate* studies need no SDP;
*absolute* cb-norm values need the Julia+MOSEK rung.

## 2. Frictions encountered (ranked by severity)

### F1 — Kraus blow-up: no channel-power, no Kraus-rank compression (HIGH)
`aic_ucp_compose` multiplies Kraus counts (`r_1 * r_2`), and
`aic_cbnorm_eigfree_ball` composes *again* internally (`Psi o Psi`). So naive
`eta(Phi^n)` carries `r^{2n}` Kraus operators even though the channel's Choi
rank is `<= dim^2`. The first run was **OOM-killed at n=14** (qubit, r=2 ->
`2^28` Kraus 2x2 matrices).

*Workaround that worked:* a `kraus -> Choi -> kraus_latd` round-trip
(`aic_ucp_kraus_to_choi` + `aic_ucp_choi_to_kraus_latd`) after each power caps
the rank at `dim^2`. **Caveat:** `choi_to_kraus_latd` is the LAPACK double path,
so this silently drops the arb certification of the Kraus rep — acceptable for
rate estimation, but a *fully* certified study would need a certified
compression (gated by the degenerate-eig wall, bead `aic-w4o.1`).

There is no `aic_ucp_power` and no `aic_ucp_compress` primitive.

### F2 — `eta` is hardwired to `||Phi^2 - Phi||`; no general `||Phi - Psi||_cb` (MEDIUM)
The public `cbnorm` API computes the idempotence defect of a *single* map. To get
`eta(Phi^n) = ||Phi^{2n} - Phi^n||` I had to build `Psi = Phi^n` and call the
defect routine on it (which redundantly recomputes `Psi^2`). A general certified
diamond/cb *distance* `||Phi - Psi||_cb` — the instrument behind application #10
(certified channel discrimination, distance-to-nearest-unitary, comparing two
tomographic noise models) — does not have an entry point. The pieces exist:
`aic_ucp_choi_diff` + `aic_cbnorm_eigfree_ball_choi` are both already exposed, so
the C rung is glue, not new math.

### F3 — Tight/absolute certified bracket needs Julia+MOSEK, absent web-side (MEDIUM)
Pure C gives only the eig-free `[||J||_F/dim, 2||J||_F]` bracket. The tight
bracket (`aic_cbnorm_certify`) requires two MOSEK feasible points produced in
Julia (`src/sdp.jl`). Julia + MOSEK are not in the web container (and the project
forbids parallel Julia). So an *absolute* certified cb value is not a pure-C,
pure-web capability today. Rates are unaffected (F1-spike shows this).

### F4 — No non-Hermitian eigenvalue routine (`geev`) (MEDIUM)
The `lambda_2` cross-check only worked because Pauli channels are self-dual
(Hermitian superoperator -> `aic_latd_eig_hermitian`). A generic channel's
superoperator / an MPS transfer operator is non-normal; the API has only
Hermitian eig (`zheev`) and SVD (`zgesvd`), no `zgeev`. The "compare decay to the
subdominant eigenvalue / spectral gap" half of several applications (#2 off the
self-dual class, #5 MPS) is therefore not turnkey — one must reach around the
project to LAPACKE directly.

### F5 — Experiment build ergonomics (MEDIUM)
`make lib` builds `libaic.so` from **all** of `src/*.c` (`SRC := $(wildcard
src/*.c)`; Makefile:27,104-105), so it **does** export the full internal API —
every non-static symbol, including `aic_ucp_compose` and
`aic_cbnorm_eigfree_ball`, is in the `.so`. The real friction is that an
experiment built like the test targets (Makefile:55 compiles `$< $(SRC)` from
scratch) recompiles all ~70 `src/*.c` every iteration (~30-40 s) instead of
linking the prebuilt `.so`. What is missing is (a) a documented one-line recipe
to link an experiment against `build/libaic.so`, and (b) an umbrella header
aggregating the ~20 public `include/aic_*.h` (today each must be `#include`d
individually). Both would cut the edit-run loop to seconds — no new library
needed.

### F6 — arb<->double bridge has no helper (LOW)
To feed an `acb_mat` superoperator into the double-path Hermitian eig I hand-rolled
the extraction (`acb_get_real`/`acb_get_imag` + `arf_get_d` per entry). A pair of
`acb_mat -> double _Complex*` / `double _Complex* -> acb_mat` (midpoint)
converters would remove a recurring papercut at the arb/double seam.

### F7 — No physical-channel constructor library (MEDIUM)
The codebase represents channels as raw Kraus arrays; I hand-built every channel
(dephasing, Pauli) by setting `acb` entries. The `tests/aic_adversarial*` corpus
is *evil matrices* for NLA stress, not *physical channels*. Most applications
(#1 noise/DFS, #3 twirls, #6 dephasing, #7 conserved quantities, #8 codes) start
by instantiating a standard channel family with a physical knob; today that is
boilerplate per study.

### F8 — Scale ceiling (expectation-setting, not a defect)
Everything is dense `dim^2 x dim^2` superoperator / Choi arithmetic. Comfortable
range is `dim <= 8-10` (a few qubits, small qudits, MPS bond dimension <= ~8).
"Many-body" here means short chains / few-body — not the thermodynamic limit.
Worth stating up front so applications are scoped accordingly.

## 3. Feature wishlist (prioritized; mapped to frictions)

| ID | Feature | Kills | Effort | Notes |
|---|---|---|---|---|
| **W1** | `aic_ucp_power(out, phi, n)` + `aic_ucp_compress(phi)` (canonical Choi-rank reduction), both number paths | F1 | S–M | Double path: reuse `choi_to_kraus_latd`. Certified arb compression gated by `aic-w4o.1` (degenerate eig) — ship double first, certified later. |
| **W2** | `aic_cbnorm_distance(lo, hi, phi, psi, prec)` — general certified cb/diamond bracket | F2 | S | C rung is glue over the already-exposed `aic_ucp_choi_diff` + `aic_cbnorm_eigfree_ball_choi`. Directly unlocks application #10. |
| **W3** | Physical-channel constructors: `dephasing`, `depolarizing`, `pauli`, unital `mix_unitaries`, `group_twirl(rep)`, `cond_expectation(subalgebra)` | F7 | M | Mirror `make_*` test-fixture style but exported. Feeds #1,#3,#6,#7,#8. |
| **W4** | `aic_latd_eig_general` (LAPACKE `zgeev`): full complex spectrum + spectral-gap helper for a superoperator/transfer op | F4 | S | Double path only (non-normal -> arb cert is `aic-w4o.1`-class). Feeds #2 (general channels), #5 (MPS). |
| **W5** | Documented one-line recipe to link experiments against the existing `build/libaic.so` (already exports the full API) + an umbrella `aic_all.h` header | F5 | S | NOT a new library — `make lib` already links all `src/*.c`. Pure docs/header plumbing; cuts the experiment loop from ~40 s to ~1 s. |
| **W6** | `acb_mat`<->`double _Complex*` (midpoint) bridge converters | F6 | S | Removes the manual per-entry extraction at the arb/double seam. |
| **W7** | MOSEK-free certified *tight* cb bracket (pure-arb SDP feasibility), or a documented fallback | F3 | L (research) | Adjacent to beads `aic-m24`/`aic-0at`. Lower priority — rates work without it; only *absolute* values need it. |
| **W8** | Committed `examples/` (not gitignored) with the cb-mixing study as a worked, regression-guarded example | F5/F7 | S | Turns throwaway scratch into a runnable reference + a smoke test of the public API. |

## 4. Which features unblock which of the ten surveyed applications

The ten candidate applications (full list in the session notes): (1) DFS/noiseless-
subsystem robustness, (2) certified cb-mixing rates, (3) exact-vs-approximate
group twirls, (4) RG/coarse-graining fixed-point diagnostic, (5) MPS transfer-op
structure, (6) quantum-to-classical / pointer-basis sharpness, (7) approximate
conserved quantities, (8) AQEC code-algebra quality, (9) empirical sharpness of
`eps = O(eta)`, (10) certified diamond-norm instrument.

| Application | Runnable today? | Blocking wishlist items |
|---|---|---|
| #2 cb-mixing (self-dual) | yes (with F1 workaround) | W1 (clean), W4 (general channels) |
| #10 diamond-norm instrument | almost — needs glue | **W2** |
| #9 eps-vs-eta sharpness | yes (estimators exist) | W3 (channel families) |
| #1 DFS robustness | yes | W3 |
| #3 group twirls | yes | W3 (twirl constructor) |
| #6 pointer-basis sharpness | yes | W3 |
| #7 approximate conserved qty | yes | W3 |
| #4 RG fixed-point | yes | W3, W4 |
| #5 MPS transfer op | partial | W4 (non-normal spectrum), unital-gauge step |
| #8 AQEC (structure half) | yes; recovery channel half | W3; `factorize` (in progress) for Enc/Dec |

**Lowest-hanging given the current surface:** #2, #9, #10 — they reuse only
green modules; #10 needs just W2.

## 5. Relationship to existing beads

Adjacent open beads (do not duplicate): `aic-w4o.2` (full arb SVD),
`aic-0at` (certified SDP upper bound), `aic-m24` / `aic-ssu` (certified cb-ball
value entry point + Julia wrapper), `aic-w4o.1` (certified degenerate eig — the
common gate for certified extraction/compression). **New** relative to those:
W1 (channel power/compress), W2 (general cb distance), W3 (channel constructors),
W4 (`geev`), W5/W6 (experiment ergonomics), W8 (examples dir). Recommend filing
W1–W6 + W8 as beads; W7 folds into the `aic-m24`/`aic-0at` research line.
