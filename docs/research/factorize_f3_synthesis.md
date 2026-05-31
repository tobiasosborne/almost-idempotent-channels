# factorize F3 — synthesis & locked decisions (the implementer's blueprint)

> **Status:** DECISION RECORD. Synthesised by the orchestrator from THREE competing
> design specs run in parallel at different thinking levels (the user's "launch
> parallel designers in competition" directive), then verified against the actual
> API. The three specs:
> - `factorize_f3_spec.md`   — Designer A (Opus, deep blueprint).
> - `factorize_f3_spec_B.md` — Designer B (Opus, adversarial-empirical / pin-everything).
> - `factorize_f3_spec_C.md` — Designer C (Sonnet, independent literature-grounded).
>
> Realizes `th_factorization` Step 5 (the UCP decode map Υ via lem_RC),
> `approximate_algebras.tex:2840-2899`. F3 follows F1 (Δ̃,Υ̃) + F2 (UCP Δ), both green.

## 0. Why a competition: the risky conventions

F3 carries the bug-class that bit `cbnorm` and `opspace` O2: **tensor-factor /
partial-trace DIRECTION**. Three independent designers converged on the same
verdicts; where the `.tex` is genuinely ambiguous, the decision is "pin it
empirically with a shape-distinguishing oracle BEFORE writing dependent code"
(the O2 `§0.5` discipline, FINDINGS §C12.O2-PIN). The convergence + the API
verification below make these LOCKED, not guesses.

## 1. LOCKED decisions (3-way consensus, orchestrator-verified)

### D1 — Partial-trace direction for `C_j` (THE #1 risk). VERIFIED.
`R_j ∈ B(L_j⊗E_j)`, **L_j the LEFT/MAJOR factor** (forced by the `X_j⊗1_{E_j}`
and `U_{js}⊗1_{E_j}` Kroneckers, left-major per `aic_mat.h`). lem_RC(i):
`R_j = 1_{L_j}⊗C_j`, so trace out **L_j (left)**:
```
C_j = (1/d_{L_j}) · aic_mat_partial_trace_left(R_j, a = d_{L_j}, b = d_{E_j})
```
Verified against `include/aic_mat.h`: `partial_trace_left` traces factor 1
(LEFT/C^a), leaving a `b×b` operator on the RIGHT (E_j); `(Tr_1 M)[j,l] =
Σ_i M[i·b+j, i·b+l]`. On `R_j = 1_{L_j}⊗C_j`, `Tr_1 = Tr(1_{L_j})·C_j = d_{L_j}·C_j`,
so the `1/d_{L_j}` recovers `C_j` exactly. **MANDATORY empirical pin (all three
designers):** build the asymmetric fixture `d_{L_j}=3 ≠ d_{E_j}=2`, `C_true=diag(1,0.5)`;
the correct `partial_trace_left` returns a 2×2 = `C_true`, the WRONG
`partial_trace_right` returns a 3×3 = `1.5·I_3` — wrong SHAPE and value, caught
instantly. Build this probe (and bake it as a regression) BEFORE any L_j code.

### D2 — W_j extraction from the UCP Δ. Per-block, L_j-major, built DIRECTLY.
`Δ(X) = Σ_j W_j†(X_j⊗1_{E_j})W_j`, `W_j: H → L_j⊗E_j`, `Σ_j W_j†W_j = 1_H`
(`.tex:2831-2837`). Route: for each block j, build the Convention-A Choi of the
per-block CP map `Δ_j: M_{d_j}→B(H)` (reuse F2's `aic_factorize_delta_block_choi`
LAYOUT with the **UCP Δ** swapped in for Δ'), extract Kraus `{D_{j,c}}` via
`aic_ucp_choi_to_kraus_latd` (Δ's Choi is degenerate → double path), then stack
**L_j-major**: `W_j[a·e_j + c, p] = D_{j,c}[a,p]` (a∈L_j, c∈E_j, p∈H), with
`e_j = E_j` = the per-block Stinespring rank.
**CONVENTION WARNING (Designer C, independently A):** do NOT obtain W_j from
`aic_ucp_kraus_to_stinespring` — its column-stack packs the ANCILLA (Kraus/F)
index LEFT (`V[a·dim_K+i, j]=K_a[i,j]`), which would put E_j in the LEFT slot,
reversing the needed L_j-LEFT ordering. Build W_j L_j-major directly.
**Gauge:** W_j is unique only up to `1_{L_j}⊗u_j` on E_j; NEVER assert C_j/W_j
entrywise — assert only the gauge-invariants `σ_max(C_j)`, `‖R_j−1⊗C_j‖`, and the
final `Υ'_j`. The η=0 oracle pins the rest.

### D3 — The per-block design for R_j (Designer B's P4, a sharp catch). PER-BLOCK, NOT the F2 join.
lem_RC's `U_{js}` are the **per-block** d_j×d_j Paulis on `L_j=C^{d_j}` (d_j² of
them, weight `d_j^{-2}`) — NOT F2's whole-B joint `aic_factorize_design_unitary`
(∏_l d_l² terms). They coincide only at m=1; at m≥2 the joint sum mixes blocks
and `R_j ≠ 1⊗C_j`. F3 enumerates per-block Paulis via `aic_dhom_pauli(out, d_j, a, b)`
directly. Also (B.1): in `Δ(U_{js})` (used in L_j and the lower-bound), `U_{js}` is
the **block-j-embedded** element of B (= U_{js} in block j, ZERO elsewhere), NOT
identity-padded — so `Δ(U_{js}) = W_j†(U_{js}⊗1_{E_j})W_j` (only the j-th Choi term).

### D4 — ξ_j = top RIGHT singular vector of C_j (all three).
So `‖C_j ξ_j‖ = σ_max(C_j)`. The LEFT vector gives the wrong norm for non-normal
C_j. ASSERT `σ_max(C_j) ≥ 1−O(η)` (lem_RC(ii), Rule 4; abort if the ball straddles
the threshold). SVD via the double-path `aic_latd_svd` / `aic_mat_singular_values`.

### D5 — The F-ancilla ordering in L_j (Designer C's "highest-risk", all three touch).
`L_j = Σ_s p_{js}(Δ(U_{js}†)⊗1_F) V W_j†(U_{js}⊗ξ_j)`, `V` = Stinespring of the
ORIGINAL Φ (`aic_ucp_kraus_to_stinespring`), `F=C^r`. Because `V` packs the ancilla
F LEFT (`V[a·dim_K+i,j]=K_a[i,j]`), the `⊗1_F` factors MUST be written **F-LEFT**:
`1_F⊗Δ(U_{js}†)` and the `Φ(X)⊗1_F` in Υ'_j as `1_F⊗Φ(X)`, to match V's layout.
η=0 with r=1 is BLIND to this (F is 1-dim); the η>0 test with **r>1** must exercise
it (compare both orderings → only the F-LEFT one reconstructs `Υ'_j Δ ≈ 1_B`).
Pin empirically with an r>1 fixture.

### D6 — Υ' Choi-PSD verdict: §C5 midpoint+Weyl (the F2/O2 pattern, do not rediscover).
The arb-assembled Υ' Choi (deep matmul chain) WILL trip `aic_mat_herm_max_eig`'s
strict relative-Hermiticity guard on near-zero radius-heavy entries (FINDINGS §C5,
exactly as F2's Δ' Choi did). Use the rigorous verdict F2 established (test helper
`dprime_is_cp`): `Cmid=(mid(C)+mid(C)†)/2` (exactly Hermitian → passes the guard),
`R=‖C−Cmid‖_F` (Weyl), `mineig_lb = −maxeig(−Cmid) − R`, CP iff `mineig_lb ≥ −tol`.

## 2. The 1-design CENTRALITY teeth (orchestrator's F2-handoff requirement — all three)

F2's tests (Δ' CP, `‖Δ'−Δ̃‖=O(η)`) hold for ANY unitary set with `Σp_s=1` — they
are BLIND to the genuine 1-design **centrality** `diag_j2` (`.tex:2776`,
`Σ_s p_{js} X U_{js}†⊗U_{js} = Σ_s p_{js} U_{js}†⊗U_{js} X` ∀X). **lem_RC(i) is the
first place centrality is load-bearing** (`R_j = 1_{L_j}⊗C_j` requires it via
Schur/commutant). F3 MUST carry, per block:
- **(T-cent-direct)** `‖Σ_s p_{js} X₀ U_{js}†⊗U_{js} − Σ_s p_{js} U_{js}†⊗U_{js} X₀‖ < 1e-9`
  on a NON-scalar X₀; mutation (perturb a weight / drop a Pauli) → O(1) RED.
- **(T-cent-struct)** `‖R_j − 1_{L_j}⊗C_j‖` ≈ machine (η=0) / O(η) (η>0). This single
  check certifies lem_RC(i), the partial-trace DIRECTION (D1), AND centrality
  simultaneously — it RED-fires if the design isn't central or the trace is wrong.

## 3. η=0 oracle (the cleanest cross-check, Rule 6)

For exact-idempotent Φ, Υ must reduce to the `th_idemp_structure` decode
`Υ = Γ∘C_M` (`aic_idemp_decompose`'s `Gamma`, `C_M`; prop_Gamma `.tex:2106`). Assert
(gauge-invariant, as MAPS): `ΥΔ = 1_B` and `ΔΥ = Φ̃` (machine-zero), and lem_RC's
`R_j = 1_{L_j}⊗C_j` with `‖C_j‖ = 1` EXACT. A wrong partial-trace direction or a
broken design is caught HERE already (η=0), before any O(η) subtlety.

## 4. Test plan & file split

- **Source (≤200 LOC/file):** `src/aic_factorize_upsilon.c` (per-block design enumerator,
  R_j, C_j via D1 partial trace, ξ_j via D4 SVD, the lem_RC asserts) +
  `src/aic_factorize_upsilon2.c` (W_j extraction D2, L_j D5, Υ'_j, Υ unitalize D6).
  Header additions to `include/aic_factorize.h`.
- **Tests (extend `tests/test_factorize.c`):** T5 η=0 oracle (§3) + T6 η>0
  (mixconj 4,2 / 5,2: Υ UCP via §C5 midpoint+Weyl, `‖Υ−Υ̃‖=O(η)`, `σ_max(C_j)≥1−O(η)`,
  centrality teeth §2, the D1 partial-trace pin, the D5 r>1 F-ordering pin). Plus the
  standalone asymmetric `d_Lj≠d_Ej` partial-trace unit probe (D1) and the diag_j2 probe.
- **Mutation teeth (Rule 7, mutation-prove each):** wrong partial-trace direction
  (shape RED via D1 asymmetric fixture); F-ordering `Δ(U)⊗1_F` instead of `1_F⊗Δ(U)`
  (D5, RED at r>1); per-block→joint design swap (D3, RED at m≥2); broken centrality
  (§2, RED); ξ_j left-vector instead of right (D4, wrong σ).

## 5. Ranked risks

1. **Partial-trace direction (D1)** — silent + shape-blind UNLESS `d_Lj≠d_Ej`.
   Mitigation: build the asymmetric probe first; the η=0 `‖R_j−1⊗C_j‖` tooth + the
   shape mismatch catch it.
2. **F-ancilla ordering (D5)** — η=0/r=1 blind. Mitigation: r>1 η>0 fixture comparing
   both orderings against `Υ'Δ≈1_B`.
3. **W_j gauge / per-block extraction (D2)** — assert only gauge-invariants.
4. **§C5 Choi-PSD false-fail (D6)** — reuse the F2 midpoint+Weyl verdict.
5. **σ_max(C_j) certification near 1−O(η)** — fail-loud on straddle (Rule 4).
6. **Composite O(η) constant for `‖Υ−Υ̃‖`** — DEFERRED (measure per-instance +
   dimension-independence canary; analytic constant chained after `aic-1bc`,
   FINDINGS §D4). NOT a wall.

## 6. Orchestration note (race-condition caveat from the parallel run)

Designer A read F2's `src/aic_factorize_delta.c` while the concurrent F2 hostile
reviewer had a mutation applied (the "s=0-only" 1-design probe), and reported it as
a "live artifact." That was a transient read of a mid-mutation file, NOT a real
defect. **Before F3 implementation begins, the orchestrator independently verifies
F2's source is clean** (`delta_prime`'s loop is the full `s∈[0,nterms)`, no probe
artifacts, 36 checks green). F3 builds on the committed F2.
