# HANDOFF.md — almost-idempotent-channels

Orientation for a fresh agent. Last updated **2026-05-30**, after an orchestrated
session that built **almost the entire `th_main` pipeline** — the constructive
proof of the main theorem (`tex:460`: any finite-dim ε-C* algebra is
`O(ε)`-isomorphic to a genuine C* algebra, with a *dimension-independent*
constant). Four modules landed, each through the full research → implement →
verify → **hostile-review → fix** → commit pipeline (the hostile review caught a
**real blocker every time** — see the per-module blocks below):
- **`corner`** (§7, bead `aic-czm`, CLOSED) — `Co_{P,Q}` compressions, `S_{P,Q}`,
  compressed product, `lem_alpha`, `lem_PQ_Hilb`, the **Ha-map**, `lem_PQR`,
  `lem_1d_proj`. 3 increments, 134 checks.
- **`projection`** (§6, bead `aic-mqf`, CLOSED) — the **constructivization crux**
  (`lem_nontriv_projection`): ambient spectral split of a Hermitian element of A,
  projected into A via `Φ̃`. 58 checks.
- **`dhom`** (§8, bead `aic-c1n`, CLOSED) — `prop_delta_hominc` + `lem_approx`
  (the multiplicativity-defect Newton iteration with the explicit Pauli diagonal).
  89 checks.
- **`errreduce`** (§8/9, bead `aic-t81`, CLOSED) — `cor_improvement`; the **c₀
  universality canary** lives here and is FLAT (dimension-independent — the
  headline `tex:484` correctness check passes). 51 checks.

Also created **`paper/FINDINGS.md`** — a living log of `.tex` typos / formula
errors / non-constructive steps / load-bearing subtleties / open escalations,
wired into the CLAUDE.md read order (step 6). **Read it.** Two confirmed `.tex`
findings this session: §A1 (`tex:1109` `lem_alpha` β subscript `Q_j`→`Q_k`) and
§A2 (`tex:1254` direct-sum diagonal is non-central for the finite Pauli design;
the correct one is the embedded sum `Σ_l D_l`).

Current state: `master` clean (pushed to `origin`), **22 test binaries green,
zero warnings**, 28 beads closed of 62. `make test` ~5 min (the n=16 canary in
`test_assoc2` U3, bead `aic-erz`; per-binary `make build/test_X && ./build/test_X`
stays quick). Read §"Channel-module conventions" before touching
`ucp`/`idemp`/`cbnorm`/`ecstar`/`assoc_ecsa`, and `paper/FINDINGS.md` §C (the
"tests that can't fail" the reviews keep catching: the **star vs plain product**
blindness, the **oblique-`Φ̃`-not-`Π_A`** projection, **rem_X2** = no funcalc
inside A, the **`aic_corner_gamma_opnorm_ub`** workaround for the `aic-qgs`
opnorm Gram false-fail).

---

## ⚠️ RESUME HERE — `cstar_build` (§9 master loop, bead `aic-097`) = the proof of `th_main`

This is the **only remaining `th_main` step** and the culmination: it assembles
the four modules above into the constructive `O(ε)`-isomorphism `v: B → A`. The
proof structure is in **shard E** (`paper/shards/shard-E-delta-homs-main-proof.md`,
"The master loop") and `tex:1414-1444`; read that first, plus `MODULE_PLAN.md`
L5. The 3-stage loop:
1. **Stage 1 (commutative skeleton):** repeatedly apply `projection` (aic-mqf) to
   split A; force the projections to be 1-dimensional (`lem_1d_proj`, in `corner`);
   classify into equivalence classes (`lem_PQR` transitivity, in `corner`); merge
   via `cor_merge_sum` + reset error via `errreduce` (aic-t81).
2. **Stage 2 (extend each class to a full matrix block):** the inductive
   `lem_extension` (`tex:1378`) — uses the **Ha-map** `Ha^Q_{P,P}` (already built
   in `corner`, the homomorphism), `lem_approx` (dhom), and `lem_merging`; then
   `errreduce` after each step.
3. **Stage 3 (merge blocks):** `cor_merge_sum` over the equivalence classes +
   `errreduce`, giving `B = ⊕_C M_{|C|}` and the `c₀ε`-isomorphism `v`.

**What `cstar_build` must still BUILD (not yet implemented):** the §9 assembly
lemmas `lem_merging` (`tex:1325`), `cor_merge_sum` (`tex:1352`), `lem_add_dim`
(`tex:1363`), `lem_extension` (`tex:1378`) — all on top of the existing
`corner` (`S_{P,Q}`/compressed product/`lem_alpha`/Ha-map), `projection`, `dhom`,
`errreduce`. These are **constructive** (shard E marks them so) but `lem_extension`
is the substantive one (Ha-map linearization + `lem_approx` inside + `lem_merging`
assembly). It is a LARGE module — budget a focused session; research shard E +
`tex:1321-1446` first, decide the loop's data model (how to thread `B`'s growing
block structure + the running `v`), then build the assembly lemmas, then the loop.

**Reuse map (all CLOSED + green):** `aic_projection_nontrivial` (the nontrivial
projection), `aic_corner_*` (Co/S/cdot/alpha/ip_1d/**ha**/dim_S/equiv_1d),
`aic_dhom_approx`+`aic_dhom_prop_bounds` (lem_approx/prop_delta_hominc),
`aic_errreduce` (the error reset), `aic_dhom_v_sigma_min` (the true inclusion
lower bound). The η=0 oracle: an exact-idempotent Φ → A a genuine C* algebra →
`cstar_build` must recover `B≅A` with `v` an **exact** isomorphism (zero defect).
The universality canary (c₀ not growing with dim A) is the headline test — the
`errreduce`/`dhom` canaries already pass; `cstar_build` must preserve it
end-to-end.

**Known open escalations feeding cstar_build** (see `paper/FINDINGS.md` §D):
`aic-3qv` (does an Ω(1) gap always exist for `projection`? — per-instance
certified now, universal guarantee open); `aic-1bc` (the analytic `c₀` from
`cor_improvement` — only the *measured* c₀ exists today); `aic-w4o.1` (certified
degenerate eig — corner `dim S`, projection gap *enclosure*, errreduce σ_min are
double-path now). None block starting cstar_build (the constructive double-path
route works; certification defers).

**th_main_ext (§10, `opspace`, bead `aic-zwo`) and the headline `factorize`
(`th_factorization`, `aic-tff`) come AFTER cstar_build** — the user scoped this
work to plain `th_main`, not the tensor extension.

---

**`aic-8hz` (this session, CLOSED, committed).** funcalc now has a globally-
convergent non-normal matrix-`sgn` reaching the full SPECTRAL `ρ(I−X²)<1` regime,
eig-free (no `aic-w4o.1` dependency). Route: Newton `½(Y+Y⁻¹)` (Higham *Functions
of Matrices* Thm 5.6, global+quadratic) gated by a **Gelfand precondition**
`‖M^k‖_F^{1/k}<1` on `M=I−X²` (`k=1` = the old op/Frobenius basin; non-normal needs
a few powers; `k_max=32`, fail loud), + an a-posteriori `‖Y²−I‖`/`‖XY−YX‖`
certificate (fires at prec=53). `aic_sgn` auto-dispatches (in-basin Newton–Schulz
byte-for-byte, else global Newton); **`aic_prop_P` relaxed to the spectral
`ρ(P²−P)<1/4`**, so `theta`/`prop_P`/`aic_assoc_regularize` now reach the full
NON-NORMAL `η<1/4` regime (`assoc_ecsa` near-boundary/large-n robustness). New
`src/aic_funcalc_global.c` (200 LOC), `tests/test_funcalc_global.c` (n=12). Teeth:
2×2 closed-form oracle (Mathematica-confirmed 2 ways to 58 digits) + `test_assoc2`
**U6** assoc payoff (oblique channel, `‖S²−S‖_op=0.42≥¼` OLD aborts, `ρ=0.21` cert
`4ρ<1` at k=6 → `S̃²=S̃` to 8.9e-70) with a **unitality** sign witness `Φ̃(1)=1`
(`‖S̃·vec(I)−vec(I)‖=4e-70` vs `√n=2` for −sgn), mutation-proven. Hostile review:
no blockers; 3 non-blocking findings all fixed inline. Deferred: scaled-Newton
(Kenney–Laub) Pareto candidate → bead **`aic-68c`** (P3; only if iteration/k counts
climb near `η→¼`). Docs: `ALGORITHM.md` §"Module funcalc — global sgn",
`docs/research/sgn_global_research.md`.

**`aic-dbo.2` Increment 1 (this session, committed; bead stays in_progress for
Inc2).** The adversarial-instance corpus: `tests/aic_adversarial.{h,c}` +
`aic_adversarial_nla.c` + driver `tests/test_adversarial.c` (n=32 self-tests, all
mutation-proven). 7 representative NLA generators spanning the distinct attack
modes + the lethal shortlist, each with a tunable knob, deterministic construction
(gen2 = fixed-LCG seed), and a certified-arb self-test asserting the knob produces
the claimed property: gen1 `jordan_t13` (1a, tex:540; ρ=t^{1/3}), gen2
`nonnormal_tri` (2c; the aic-8hz op-vs-ρ gap regime, `‖[X,X†]‖_F` 0.07→172 as
c 0.1→10), gen3 `degenerate_proj` (4b, lethal #1; exact {0,1} projector the
simple-eig path can't isolate), gen4 `near_degen_herm` (4a; gap knob), gen5
`graded_diag` (5c; κ knob to 1e8), gen6 `boundary_x2I` (7a, lethal #5; `‖X²−I‖`
straddles 1), gen7 `propP_delta` (7b; `‖P²−P‖=δ` at the ¼ edge). The Makefile now
links `tests/aic_*.c` (TEST_HELPER_SRC) into every test/bench. **The corpus already
paid off**: gen5 surfaced a real substrate bug → bead **`aic-2yo`** (FIXED + closed
this session): `aic_mat_int_assert_hermitian`'s ABSOLUTE tol false-failed on a
graded/ill-conditioned Gram (the flagged asymmetry was the ball RADIUS, which grows
with magnitude), so `aic_mat_singular_values`/`opnorm` ABORTED for condition number
≳ 1e2 (precision-independent). Fix: a RELATIVE + absolute-floor tol
`tol·(1+|H_ij|+|H_ji|)`, refactored into a non-aborting predicate
`aic_mat_int_is_hermitian` so the teeth (genuine non-Hermitian still rejected) are
testable; `test_mat` test7 (n=29) mutation-proven, full `make test` green.
Increment 2 (still open on `aic-dbo.2`): the DOMAIN/channel families
(`docs/adversarial/domain.md` — several gated by unbuilt modules) + the remaining
NLA families; unblocks the universality canary `aic-dbo.3` once the `ε~c/dim` dim-
sweep instances land. `test_adversarial.c` (451 LOC) joins the oversized-test-file
pile `aic-w4o.4`.

**⚠️ RESUME HERE — `aic-92f` (`assoc_ecsa`, `th_almost_idemp`) is CLOSED; the
paper's `th_almost_idemp` is now VERIFIED end-to-end.** This session built
Φ̃=θ(2Φ−1) (regularization to an exact-idempotent superoperator), A=ImgΦ̃, the
Choi–Effros star X⋆Y=Φ̃(XY), and the full extended-O(η)-$C^*$ axiom verification —
including the **dimension-independence universality canary** (flat, so the paper's
universal-constant claim `.tex:2193` HOLDS; no stop-condition). Two Opus hostile
reviews each caught a "test that cannot fail" (untested foundational
`aic_mat_herm_max_eig` near-zero fallback; toothless canary threshold; unasserted
C* defect); all fixed + mutation-proven. Closed `aic-3qq` (Π_A teeth on the η>0
non-polar-closed A, `test_assoc2` U4) along the way.

**The paper's FINAL headline is `aic-tff` (`factorize`, `th_factorization`
`.tex:2730`)** — the approximate factorization Φ*≈Enc∘Dec through a genuine $C^*$
algebra B. It is **NOT directly unblocked**: th_factorization Step 3 invokes
`th_main_ext`, i.e. the entire incremental $C^*$-algebra construction (E2-eps:
`unitfix`/`unitary`/`projection`/`corner`; E3-mainthm: `dhom`/`errreduce`/
`cstar_build`; §10 `opspace`), and its proof is an OUTLINE with open escalations
(shard H §"Open implementation questions"; MODULE_PLAN escalations 1,4,5). So the
real next constructive frontier is the **th_main pipeline**. `bd ready` P1: the
recurring `aic-w4o.1` (certified DEGENERATE Hermitian eig — gates the certified
halves of `projection`/Kraus/subspace extraction, deferred all session) and
`aic-dbo.2` (adversarial instance generators). `aic-8hz` (globally-convergent
non-normal `sgn`) is now **CLOSED** (see the top block); `aic-erz` (cut the canary
runtime) stays open, and `aic-4c7` stays open (the U2/U3 O(η) sweep used an op-norm
η-PROXY; closing it needs the certified cb-η normalization to rule out a proxy↔cb-η
dimension factor).

**`assoc_ecsa` (this session, `aic-92f` Increments 1+2, commits `96f381e` + this
one).** Realizes `th_almost_idemp` (`.tex:2162-2237`) as a finite-dim algorithm.
- **Φ̃=θ(2Φ−1)** is `prop_P` applied to Φ at the **n²×n² superoperator** level:
  build S_Φ column-by-column via the tested `aic_ucp_apply` on E_pq (row-major
  vec; oracle-tested + Kronecker cross-check `S_Φ=Σ K_a†⊗K_aᵀ`), then
  `aic_assoc_regularize`=`aic_prop_P(S_Φ)` (reuse `funcalc`; default sgn=
  Newton-Schulz, eig-free → sound on the NON-NORMAL superoperator). Φ̃ is NOT CP
  (`.tex:363`). `include/aic_assoc.h`, `src/aic_assoc_superop.c`,
  `src/aic_assoc_regularize.c`; `tests/test_assoc.c` n=93.
- **A=ImgΦ̃** = range of the OBLIQUE idempotent S̃, extracted by thin SVD
  (`dim_A=round(trace)` cross-checked vs SVD-gap@0.5; top-dim_A left singular
  vectors → Frobenius-ON basis {B_k}). Nonzero σ(S̃)≥1, =1 iff Φ̃ HS-self-adjoint
  (self-dual Φ), >1 genuinely oblique (`compress_idemp`, σ=√3). `aic_assoc_extract.c`.
- **Star + verification:** the `ecstar` data model was generalized BACKWARD-
  COMPATIBLY (`star_phi`/`star_ctx` apply-fn seam; Kraus path byte-identical →
  `test_ecstar` n=109 unchanged) so X⋆Y=Φ̃(XY) feeds the existing axiom-defect
  estimators. `aic_assoc_algebra.c`; `tests/test_assoc2.c` n=139: U1 η=0 oracle
  (machine-zero, A matches `aic_idemp`'s ImgΦ), U2 real η>0 (assoc+cstar O(η),
  vanishing, teeth), U3 universality canary (flat over dim_A=4..16, threshold 2.0,
  mutation-proven vs linear+√(dim) injection), U4 Π_A on non-polar-closed A, U5
  oblique σ_max.
- **Load-bearing decisions/findings:** (i) funcalc's sgn asserts the FROBENIUS
  basin 4‖S²−S‖_F<1 — dimension-dependent, holds at η=0/small η, TRIPS near η→1/4
  / large n → bead `aic-8hz`. (ii) A is built from the MIDPOINT of S̃ (funcalc-
  iteration balls ~5e-75 trip `aic_mat_opnorm`'s Hermiticity tol when squared into
  the Gram); rigor lives in I1's η=0 + non-normal A1-vs-A2 cross-checks (same class
  as the double-extraction / `aic-w4o.1` deferral). (iii) A=ImgΦ̃ is EXACTLY
  idempotent, so its associator vanishes for self-dual/compatible families —
  genuine non-associativity needs two INCOMPATIBLE algebras mixed (the η>0 family
  is dep(d)⊕conj(dep(d)); obliqueness ≠ non-associativity, distinct properties).
  (iv) A foundational fix landed in `src/aic_mat_spectral.c`: `aic_mat_herm_max_eig`
  falls back to the rigorous eig-free [−‖H‖_F,‖H‖_F] on a non-finite certifier
  radius (near-zero degenerate Gram), now with a direct mutation-proof test in
  `test_mat.c` n=23. Research: `docs/research/assoc_ecsa_research.md`.

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
before committing — the pre-commit hook also syncs it). **No GitHub CI.** A git
**remote now exists** (`origin`, GitHub; user directive 2026-05-30) — **push after
every commit** (`git pull --rebase && git push`). The `.beads/issues.jsonl` in git
is the bead persistence; no `bd dolt` remote is configured (`bd dolt push` is N/A).

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
