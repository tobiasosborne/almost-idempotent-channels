# factorize F3 — DESIGN SPEC B (adversarial-empirical lens): the UCP decode map Υ via lem_RC

**Bead:** aic-tff (factorize F3; `th_factorization` Step 5, the UCP decode Υ; `lem_RC`).
**Date:** 2026-05-31
**Ground truth:** `paper/src/approximate_algebras.tex` (cited by line throughout).
**Status:** DESIGN (read-only on `.c`/`.h`; no build/Julia in this pass).
**Lens (B):** trust NOTHING; pin EVERY index/tensor/partial-trace/gauge convention
empirically against an INDEPENDENT closed-form oracle BEFORE committing. Mirrors the
opspace O2 convention-pinning discipline (`docs/research/opspace_o2_design.md §0.5`,
`paper/FINDINGS.md §C12.O2-PIN`): the design's "derived" convention was WRONG twice;
only the independent-oracle probe is trustworthy.

**Reads (cited):** `.tex:2755-2899` (tilde_Del/PhiDelta1-3 :2804-2829, Choi_Delta :2831-2837,
lem_RC :2840-2857, ξ_j/L_j/Υ' :2859-2870, the Υ'Δ≈1_B calc :2871-2893, Υ unitalize :2896-2899;
diag_j2 :2775-2778; D_j/joint-U_s :2771-2784; Δ' :2786-2801);
`paper/shards/shard-H-almost-idemp-factorization.md` (lem_RC L193-216, open Q #5/#6/#7);
`docs/research/factorize_d4_research.md §3` (lem_RC IN DETAIL), `§1 Step 5`, `§4`;
`docs/research/opspace_o2_design.md §0.5`; `paper/FINDINGS.md §C12.O2-PIN, §C5, §A2, §C2, §C3, §C12.O2-SUBTLETY`;
`include/aic_mat.h:43-60,126-144` (partial-trace + Kronecker convention);
`include/aic_ucp.h:6-50,89-93` (Choi Conv A, Stinespring V layout);
`src/aic_ucp_carrier.c:41-54` (the verbatim V column-stack);
`include/aic_idemp.h:80-93` (η=0 oracle Δ/Γ/W_j);
`include/aic_funcalc.h:88-105` (xpow -1/2 + basin);
`include/aic_factorize.h` (F1/F2 APIs to reuse) + `src/aic_factorize_delta.c:16-90,166-177`;
`include/aic_dhom.h:103-122` (B blocks, vE); `include/aic_opspace.h:182-198` (build_vinv).

---

## A. CONVENTION PIN TABLE (the headline — do this FIRST, it is the GO/NO-GO gate)

Every row is UNKNOWN until pinned by an independent closed-form oracle on a SMALL
ASYMMETRIC fixture where the two candidate readings give VISIBLY DIFFERENT,
hand-computable answers. **Do not write the F3 production code until every row reads
GREEN on its probe.** This is exactly the opspace O2 lesson (`§C12.O2-PIN`): the design's
*derived* partial-trace direction and normalization were both wrong; the probe corrected
them. F3 has the SAME landmines (the partial-trace direction bit cbnorm AND opspace O2,
`opspace_o2_design.md §2.5`), plus three new ones (W_j gauge, the per-block U_{js} indexing,
the ξ_j singular-vector side).

| # | Convention (where it bites) | Candidate readings | Distinguishing oracle (asymmetric fixture) | Expected value per reading |
|---|---|---|---|---|
| **P1 (HIGHEST RISK)** | **Partial-trace DIRECTION for `C_j = (1/d_{L_j}) Tr_{L_j}(R_j)`** (`.tex:2849`, shard L207). `R_j ∈ B(L_j ⊗ E_j)`; trace out the `L_j` factor, keep `E_j`. | (a) `L_j` is the LEFT/MAJOR factor of `R_j` → trace the LEFT → `aic_mat_partial_trace_left(C_j, R_j, d_Lj, d_Ej)` → keeps RIGHT = E_j (size d_Ej). (b) `L_j` is the RIGHT/MINOR factor → `aic_mat_partial_trace_right` → keeps LEFT. | Build `R_j = 1_{L_j} ⊗ C_j^true` DIRECTLY with an **asymmetric** hand-picked `C_j^true` (a non-scalar `d_Ej × d_Ej`, e.g. `diag(1, 0.5)` with `d_Ej=2`) and `d_Lj ≠ d_Ej` (e.g. `d_Lj=3`), using `aic_mat_kronecker(R_j, I_{d_Lj}, C_j^true)` (LEFT-major, `aic_mat.h:126-128`). Then trace and compare to `C_j^true`. | (a) `aic_mat_partial_trace_left(out, R_j, a=d_Lj, b=d_Ej)` keeps the RIGHT (E_j) factor → `out = Tr(I_{d_Lj})·C_j^true = d_Lj·C_j^true` (`aic_mat.h:59` `Tr_1(A⊗B)=Tr(A)B`), so `(1/d_Lj)·out = C_j^true` EXACTLY. **GREEN.** (b) `aic_mat_partial_trace_right(out, R_j, a=d_Lj, b=d_Ej)` keeps the LEFT (L_j) factor → `out = Tr(C_j^true)·I_{d_Lj} = 1.5·I_3` — a 3×3 scalar, WRONG shape (d_Lj×d_Lj not d_Ej×d_Ej) and WRONG value. **RED, visibly.** |
| **P2** | **The `R_j` Kronecker factor ORDER** — is `R_j` stored as `L_j ⊗ E_j` (L major) or `E_j ⊗ L_j`? This is set by `W_j: H → L_j ⊗ E_j` (`.tex:2837`) and the `1_{E_j}` placement in `U_{js}^† ⊗ 1_{E_j}` (`.tex:2843`). | (a) L_j MAJOR / E_j MINOR (matches `aic_mat.h:43-60` LEFT-major Kronecker AND `aic_ucp.h:28` Choi K-major). (b) E_j MAJOR / L_j MINOR. | `U_{js}^† ⊗ 1_{E_j}` must, when conjugating `W_j W_j^†`, act ONLY on the `L_j` factor and leave `E_j` alone (that is the WHOLE mechanism of `diag_j2` centrality, `.tex:2849`). Build `W_j W_j^†` block-diagonally-asymmetric and check `(U^†⊗1) M (U⊗1)` rotates the intended factor: pick `U = X^1 Z^0` (cyclic shift on L_j) and `M = P ⊗ Q` with `P,Q` distinct → `(U^†⊗1)(P⊗Q)(U⊗1) = (U^†PU)⊗Q`. | (a) result `= (U^†PU) ⊗ Q` with Q untouched → CORRECT (E_j inert). (b) the kron in `aic_mat_kronecker(out, Ujs, I_Ej)` would put `U` on the MAJOR factor; if the convention claimed E_j-major, the `⊗1_{E_j}` would mis-target → Q gets rotated, P inert → visibly different matrix. **The probe is `aic_mat_kronecker(KU, Ujs_dag, I_Ej)` then check it commutes with `1_{L_j}⊗(any C)`.** |
| **P3** | **`W_j` extraction GAUGE from Δ's Choi/Stinespring** (`Choi_Delta` `.tex:2831-2837`: `Δ(X)=Σ_j W_j^†(X_j⊗1_{E_j})W_j`, `W_j: H→L_j⊗E_j`). Δ is the UCP map F2 produced (`aic_factorize_delta`). | (a) `W_j` = the block of the Stinespring isometry `V_Δ` of Δ that lands in `L_j⊗E_j` (`aic_ucp_kraus_to_stinespring`, the `V[a·dim_K+i,j]` column-stack, **F-MAJOR/K-MINOR**, `aic_ucp.h:89-93`, `aic_ucp_carrier.c:52`). (b) `W_j` from a Choi-→-Kraus eigendecomposition of `C_Δ`. (c) any other isometry with the SAME `W_j^†(·⊗1)W_j` action (gauge freedom: `W_j → (1_{L_j}⊗u_j)W_j` for unitary `u_j` on E_j leaves Δ invariant but CHANGES `R_j` by `u_j`-conjugation of `C_j`). | **The gauge is NOT free for `C_j` — but `‖C_jξ_j‖` and the final `Υ'_j` ARE gauge-invariant** (proof: §C below). So the oracle is NOT "recover a specific W_j" but "the end product `Υ'_j(Δ(Y)) ≈ Y_j` is gauge-independent". Probe: build two W_j gauges (multiply by random unitary `u_j` on E_j), recompute `R_j → C_j → ξ_j → L_j → Υ'_j`, assert `Υ'_j` agrees to O(η). | Both gauges give the SAME `Υ'_j(Δ(Y))` (to O(η)); but `C_j` itself differs by `u_j C_j^{(b)} u_j^†`-type conjugation — so **never test C_j entrywise across gauges; test only the gauge-invariant outputs** (`σ_max(C_j)`, `‖R_j−1⊗C_j‖`, `Υ'`). |
| **P4** | **The per-block `U_{js}` indexing in `R_j` / `L_j`** — `.tex:2843,2862` sum over `s` with the PER-BLOCK design `p_{js}, U_{js}` on `B(L_j)` (`d_j²` terms, `d_j×d_j` unitaries), NOT the JOINT `U_s = U_{1,s_1}⊕…⊕U_{m,s_m}` that F2's `aic_factorize_design_unitary` builds (`.tex:2782-2783`, `src/aic_factorize_delta.c:61-90`). | (a) `U_{js}` = the genuine per-block `d_j×d_j` Pauli `S_ab` (from `aic_dhom_pauli(S, d_j, a, b)`), `s∈[0,d_j²)`, summed with weight `d_j^{-2}`. (b) `U_{js}` = the joint `n_B×n_B` `U_s` (F2's `design_unitary`) restricted/projected to block j. | At `m=1` (single block) (a) and (b) COINCIDE — direction-blind, like the symmetric-channel O2 case (`opspace_o2_design.md §2.5`). At **m≥2** they differ: the joint `U_s` has `∏_l d_l²` terms; the per-block has `d_j²`. Probe on `B = M_2 ⊕ M_3` (m=2): count terms (4 vs 36 for block 0) and check `R_j = 1_{L_j}⊗C_j` (P1's structural check) — only the PER-BLOCK sum yields a clean `1⊗C_j` for block j; the joint sum mixes blocks. | (a) `R_0` (block 0, d_0=2) is a clean `1_2⊗C_0` with `d_0²=4` terms. (b) using joint `U_s` (36 terms) `R_0` is NOT of the form `1⊗C` (centrality fails per-block) → `‖R_0−1⊗C_0‖` is O(1), **RED**. So **F3 must build the per-block design itself (a) — reuse `aic_dhom_pauli` per block, NOT F2's joint `design_unitary`.** |
| **P5** | **ξ_j singular-vector SIDE** — `ξ_j∈E_j` unit vector, `‖C_jξ_j‖≈1` (`.tex:2859`, shard #5). | (a) `ξ_j` = top RIGHT singular vector of `C_j` (so `‖C_jξ_j‖=σ_max(C_j)`). (b) top LEFT singular vector (gives `‖C_j^†ξ_j‖=σ_max`, NOT `‖C_jξ_j‖`). | For NON-normal `C_j`, left and right singular vectors differ; `‖C_jξ‖` for the LEFT vector ≠ σ_max in general. Build asymmetric NON-normal `C_j` (e.g. `[[1,0.3],[0,0.9]]`), compute both, evaluate `‖C_jξ‖`. | (a) RIGHT vector: `‖C_jξ_j‖ = σ_max(C_j)` EXACTLY (the defining property). (b) LEFT vector: `‖C_jξ‖ < σ_max` (measurably; here σ_max≈1.06 but `‖C_j·(left vec)‖` differs). **The `.tex` wants `‖C_jξ_j‖≈1` so (a) is required.** Note: `C_j` is generically NOT Hermitian (`R_j=Σ p U^† WW^† U` is PSD-conjugate-summed, but `C_j` extracted via partial trace of a PSD `R_j` IS PSD/Hermitian — see caveat below), so left=right when C_j IS Hermitian; the probe CONFIRMS Hermiticity first, then (a) is the safe universal choice. |

### A.0 The SINGLE highest-risk convention + its GO/NO-GO probe

**P1 — the partial-trace direction for `C_j = (1/d_{L_j}) Tr_{L_j}(R_j)`.** This is the
exact landmine that bit BOTH the original cbnorm work (`HANDOFF:340`) and opspace O2
(`opspace_o2_design.md §2.5`, `FINDINGS §C12.O2-PIN`): a wrong partial-trace direction is
SILENT (returns a plausible matrix) and only an asymmetric, hand-computed oracle catches
it. Here the trap has an extra tooth: getting it wrong on a `d_Lj = d_Ej` block is
**shape-blind** (both directions return a square of the same size), so the fixture MUST
have `d_Lj ≠ d_Ej`.

> **GO/NO-GO probe `probe_f3_ptrace` (run before ANY F3 code).**
> Pick `d_Lj = 3`, `d_Ej = 2`. Set `C_true = diag(1.0, 0.5)` (2×2, asymmetric, trace 1.5).
> Assemble `R = aic_mat_kronecker(I_3, C_true)` (6×6, LEFT-major per `aic_mat.h:126-128`,
> so `R[i·2+r, k·2+s] = I_3[i,k]·C_true[r,s]` — the `L_j ⊗ E_j = I_3 ⊗ C_true` layout the
> paper's `R_j ∈ B(L_j⊗E_j)` demands, P2-pinned).
> - Reading (a): `aic_mat_partial_trace_left(out, R, a=3, b=2)` returns a **2×2**;
>   `(1/3)·out` MUST equal `C_true = diag(1, 0.5)` to ~1e-15. **GO** ⟺ this holds.
> - Reading (b): `aic_mat_partial_trace_right(out, R, a=3, b=2)` returns a **3×3** =
>   `1.5·I_3`. Different shape, different value. If F3 code path produces this, **NO-GO** —
>   the direction is inverted; STOP and re-derive before proceeding.
>
> **Verdict from the header convention (`aic_mat.h:55-60`, to be CONFIRMED by the probe,
> NOT trusted):** "`Tr_1` keeps the SECOND (RIGHT) factor" and "`Tr_1(A⊗B)=Tr(A)B`". With
> `R_j = 1_{L_j}⊗C_j` (L_j LEFT, E_j RIGHT — P2), tracing `L_j` = tracing the LEFT = `Tr_1`
> = **`aic_mat_partial_trace_left(C_j, R_j, d_Lj, d_Ej)`**, giving `Tr(1_{L_j})·C_j =
> d_Lj·C_j`, so `C_j = (1/d_Lj)·partial_trace_left(R_j)`. The header NAMES it "left" because
> it traces the LEFT factor; the paper CALLS it "Tr_{L_j}" because L_j is what's traced —
> **same operation, and L_j sits in the LEFT slot, so the names agree, but THIS COINCIDENCE
> IS EXACTLY THE TRAP** (in O2 the input factor sat in the MAJOR slot and the naive reuse of
> `partial_trace_right` was wrong). The probe with `d_Lj≠d_Ej` removes all doubt.

---

## B. THE CONSTRUCTION, STEP BY STEP (mapped to existing APIs, `.tex`-cited)

Inputs available when F3 runs (all built/certified upstream): the `aic_factorize F`
handle (F1+F2) carrying `F->v` (the iso `v:B→A`, `aic_dhom_v`, `include/aic_dhom.h:117-122`),
`F->phi` (the ORIGINAL UCP `Φ`, `aic_ucp_kraus`), `F->Aec` (A = Img Φ̃ + S_tilde),
`F->vinvB` (`v^{-1}(B_k)`), and the UCP encode `Δ` via `aic_factorize_delta` (F2). Block
data: `B = ⊕_{j=1}^m M_{d_j}`, `d_j = F->v->B->d[j]`, `n_B = F->v->B->n_B`, `N = F->N`
(`include/aic_factorize.h:98-108`).

The decode is built per block `j`, then assembled as `Υ' = (Υ'_1,…,Υ'_m)`.

| step | object | formula (`.tex`) | shape | API / route |
|---|---|---|---|---|
| **S0** | `Δ` as a UCP `aic_ucp_kraus` `B(H)→B(H)`, or its Kraus list | F2's `Δ(X)=Δ'(I)^{-1/2}Δ'(X)Δ'(I)^{-1/2}` (`.tex:2799`) | `Δ: M_{n_B}→M_N` (domain B's rep, codomain H) | F2 gives `aic_factorize_delta` as a function pointer-able map. For Choi/Stinespring we need Δ's KRAUS form: assemble `C_Δ = Σ_{p,q} E_{pq}^{(n_B)}⊗Δ(E_{pq})` (Convention A, `aic_ucp.h:26`) then `aic_ucp_choi_to_kraus_latd` (double eig — Δ's Choi is generically degenerate, `aic_ucp.h:46-50`). **NEW helper** `aic_factorize_delta_kraus`. |
| **S1** | `W_j: H → L_j⊗E_j` — the Choi/Stinespring block of Δ landing in block j | `Choi_Delta` `.tex:2831-2837`, `Δ(X)=Σ_j W_j^†(X_j⊗1_{E_j})W_j`, `Σ_j W_j^†W_j=1_H` | `W_j` is `(d_j·d_{E_j}) × N` | From the Stinespring `V_Δ` (`aic_ucp_kraus_to_stinespring`, `V: H→K⊗F`, here `K=M_{n_B}`-rep). `W_j` = the rows of `V_Δ` whose `K`-index lies in block j's `L_j` range × the ancilla `F` (which IS `E_j`). **GAUGE: P3 — pin via the gauge-invariance oracle, not entrywise.** `d_{E_j} = r_Δ` (Δ's Kraus count / ancilla dim). |
| **S2** | `U_{js}`, `p_{js}` — PER-BLOCK design on `L_j` | `D_j=Σ_s p_{js}U_{js}^†⊗U_{js}`, `diag_j2` `.tex:2775-2778` | `U_{js}` is `d_j×d_j`; `s∈[0,d_j²)`; `p_{js}=d_j^{-2}` | **P4: build per-block** with `aic_dhom_pauli(U, d_j, a, b)` (`include/aic_dhom.h`), `s↔(a,b)=(s/d_j, s%d_j)`. **NOT** F2's joint `aic_factorize_design_unitary`. |
| **S3** | `R_j = Σ_s p_{js}(U_{js}^†⊗1_{E_j}) W_jW_j^† (U_{js}⊗1_{E_j})` | `R_def` `.tex:2843` | `(d_j·d_{E_j})²` | `M = W_jW_j^†` (`d_j·d_{E_j}` square, PSD). Per `s`: `KU = aic_mat_kronecker(U_{js}^†, I_{d_{E_j}})` (P2: U on L_j MAJOR), `term = KU·M·KU^†`, accumulate `p_{js}·term`. |
| **S4** | `C_j = (1/d_{L_j}) Tr_{L_j}(R_j)`, and structural check `R_j ≈ 1_{L_j}⊗C_j` | `.tex:2846` (i); shard L207 | `C_j` is `d_{E_j}×d_{E_j}` | **P1: `aic_mat_partial_trace_left(Cj_raw, R_j, d_Lj=d_j, d_Ej)`, then `C_j = (1/d_j)·Cj_raw`.** Then assert `‖R_j − aic_mat_kronecker(I_{d_j}, C_j)‖_op = O(η)` (the centrality teeth, §C). |
| **S5** | `ξ_j∈E_j`, `‖C_jξ_j‖ = σ_max(C_j) ≥ 1−O(η)` | `.tex:2859`; lem_RC(ii) `.tex:2846` | `ξ_j` is `d_{E_j}×1` | **P5: top RIGHT singular vector of `C_j`.** Double-path SVD (`C_j` may be degenerate); certify `σ_max(C_j)≥1−O(η)` in arb (`aic_mat_singular_values` if simple, else double + arb opnorm of `C_jξ_j`). **Rule 4: abort if the ball straddles `1−O(η)`** (a lem_RC precondition failure — shard #5, `factorize_d4_research.md §3`). |
| **S6** | `V` — Stinespring of the ORIGINAL `Φ`: `Φ(X)=V^†(X⊗1_F)V` | `.tex:2859` | `V: H→H⊗F`, `(N·r_Φ)×N` | `aic_ucp_kraus_to_stinespring(V, F->phi)` (`aic_ucp.h:89-93`). `F = C^{r_Φ}`, `r_Φ = F->phi->r`. **V layout F-MAJOR**: `V[a·N+i, j]=K_a[i,j]` (`aic_ucp_carrier.c:52`). |
| **S7** | `L_j = Σ_s p_{js} (Δ(U_{js}^†)⊗1_F) V W_j^† (U_{js}⊗ξ_j)` | `.tex:2862-2863` | `L_j: L_j → H⊗F`, `(N·r_Φ)×d_j` | Per `s`: (i) `Δ(U_{js}^†)` — but `U_{js}` is `d_j×d_j` on block j; embed into B's `n_B×n_B` rep (block j, identity elsewhere? **NO** — see §B.1 note) then `aic_factorize_delta`. (ii) `(U_{js}⊗ξ_j)`: `d_j×d_j` ⊗ `d_{E_j}×1` → `(d_j·d_{E_j})×d_j` (P2 order). (iii) products. Accumulate `p_{js}·(…)`. |
| **S8** | `Υ'_j(X) = L_j^† (Φ(X)⊗1_F) L_j` (manifestly CP) | `.tex:2866-2869` | `Υ'_j: B(H)→B(L_j)`, output `d_j×d_j` | `Φ(X)=aic_ucp_apply(F->phi, X)` (`aic_ucp.h:80-83`); `Φ(X)⊗1_F` = `aic_mat_kronecker(PhiX, I_{r_Φ})`; then `L_j^†·(…)·L_j`. CP by construction (a `L^†(·)L` congruence of the CP `X↦Φ(X)⊗1`). |
| **S9** | `Υ'(X) = (Υ'_1(X),…,Υ'_m(X))` as a block-diagonal element of B | `.tex:2867` | `n_B×n_B` block-diag | place each `Υ'_j(X)` (`d_j×d_j`) in block j. |
| **S10** | `Υ(X) = Υ'(1_H)^{-1/2} Υ'(X) Υ'(1_H)^{-1/2}` (UCP) | `.tex:2896-2899` | `n_B×n_B` | `Υ'(1_H)` is `n_B×n_B` block-diag; assert `‖Υ'(1_H)−1_B‖_op<1` (Rule 4, the xpow basin, `aic_funcalc.h:88-96`, shard #7), then `aic_funcalc_xpow(…, -0.5, 1.0)` per block (block-diagonal so the inverse-sqrt is per-block). Congruence preserves CP; unital by construction. |

### B.1 The `Δ(U_{js}^†)` subtlety (an indexing trap inside S7)

`U_{js}` lives on `L_j` (`d_j×d_j`), but `Δ` takes an element of `B = ⊕_l M_{d_l}` (the
`n_B×n_B` block-diagonal rep). In `.tex:2862`, `Δ(U_{js}^†)` means **Δ applied to the
B-element that is `U_{js}^†` in block j and ZERO in all other blocks** — i.e. `X_j = U_{js}^†`,
`X_{l≠j}=0` in the `Δ(X)=Σ_j W_j^†(X_j⊗1)W_j` decomposition. This is the `(X=(X_1,…,X_m)`
with a single nonzero `X_j)` reading the paper uses explicitly in the lem_RC proof
(`.tex:2849` "with a single nonzero `X_j`"). So embed `U_{js}^†` as `E_i`-coords supported on
block j only, NOT as an identity-padded block-diagonal unitary. **This is a SEPARATE pin
from P4** (P4 is which Pauli set; this is the zero-padding of the OTHER blocks). Probe:
`Δ(embed_block_j(U))` must equal `W_j^†(U⊗1_{E_j})W_j` (the single-block Choi_Delta term)
to O(η) — if you pad with identity instead of zero, you pick up `Σ_{l≠j} W_l^†(1⊗1)W_l`
spurious terms and the equality fails by O(1).

### B.2 File / API additions

All new code in `factorize`; reuses closed primitives. **NEW** functions:
`aic_factorize_delta_kraus` (S0, Δ→Kraus via Choi+latd), `aic_factorize_Wj` (S1),
`aic_factorize_Rj_Cj` (S3-S4, with the P1 partial trace), `aic_factorize_xi` (S5),
`aic_factorize_Lj` (S7), `aic_factorize_upsilon_prime` (S8-S9),
`aic_factorize_upsilon_build` + `aic_factorize_upsilon` (S10). The handle `aic_factorize`
gains `upsI_invsqrt` (Υ'(1_H)^{-1/2}, owned) and per-block cached `L_j`/`ξ_j`, mirroring F2's
`deltaI_invsqrt`/`delta_ready` (`include/aic_factorize.h:106-108`).

---

## C. THE 1-DESIGN CENTRALITY TEETH (the load-bearing structural check)

**Why this is where centrality becomes load-bearing.** F2's tests (the CP teeth, per-block
Choi PSD) **cannot distinguish a genuine 1-design from any unitary set with `Σp=1`** — a
broken "design" (e.g. a single Pauli, or a set that is `Σp=1` but NOT central) can still
give a CP `Δ'` (any `Σ p_s Φ(Z_s^†Z_s) ≥ 0` is CP regardless of centrality). The centrality
property `diag_j2` (`.tex:2775-2778`, `Σ_s p_{js} X U_{js}^†⊗U_{js} = Σ_s p_{js}
U_{js}^†⊗U_{js} X`) is FIRST load-bearing at **lem_RC(i)**: it is the SOLE reason
`R_j = 1_{L_j}⊗C_j` (`.tex:2849`: "Due to property `diag_j2`, `R_j` commutes with `X⊗1_{E_j}`
… hence `R_j=1_{L_j}⊗C_j`"). A non-central design gives an `R_j` that is NOT of the form
`1⊗C`, and the partial-trace extraction `C_j=(1/d)Tr_{L_j}(R_j)` then silently discards the
`L_j`-dependence — a "test that cannot fail" unless we check the structure directly.

**Three teeth (all required):**

1. **Direct `diag_j2` check (the design itself).** For the per-block design `{p_{js},U_{js}}`
   (S2), pick a random non-scalar `X∈M_{d_j}` and assert
   `‖Σ_s p_{js}(X U_{js}^†⊗U_{js}) − Σ_s p_{js}(U_{js}^†⊗U_{js} X)‖_op ≤ tol` (`.tex:2776`).
   This is a property of the DESIGN alone (no W_j), so it isolates "is the Pauli set actually
   a 1-design". The generalized-Pauli set satisfies it EXACTLY (`§A2`: the single-block Pauli
   sum reproduces the Haar second moment; centrality is the first-moment statement and holds
   per block). **Tooth value: machine-zero for the genuine design.**

2. **Structural `‖R_j − 1_{L_j}⊗C_j‖ ≈ 0` (lem_RC(i) on the REAL W_j).** After S3-S4,
   assert `‖R_j − aic_mat_kronecker(I_{d_j}, C_j)‖_op ≤ O(η)` (at η=0: EXACTLY 0, since
   `R_j=1⊗C_j` is exact when Δ is exactly the η=0 isometry-block — `factorize_d4_research.md
   §3` η=0 oracle). This is the CONSTRUCTIVE verification of lem_RC(i) (shard L207
   "cross-check `‖R_j−1⊗C_j‖=O(η)` certifies (i)"). It catches BOTH a non-central design AND
   a wrong W_j gauge that breaks the tensor structure.

3. **The mutation that breaks centrality (Rule 7 non-vacuity).** Replace the per-block design
   with a **single Pauli** `{U_{j,0}=I}` (or a `Σp=1` but non-central pair, e.g. `{0.5·I,
   0.5·X}` which is NOT a 1-design). Then: (mutation a) tooth #1 jumps to O(1)
   (`‖diag_j2 defect‖ ≈ 0.5–1`, the `§A2` measured scale for non-central sums); (mutation b)
   tooth #2 `‖R_j−1⊗C_j‖` jumps to O(1) — `R_j` is no longer `1⊗C`. **Both must go RED**, or
   the centrality tooth is inert. (Mirror `§A2`'s measured `‖XD−DX‖≈0.54` for the broken
   joint design.)

**Gauge-invariance note (ties P3 to C).** `C_j` is gauge-DEPENDENT (`W_j→(1⊗u_j)W_j` conjugates
`C_j` by `u_j`), but `σ_max(C_j)` and the centrality structure `‖R_j−1⊗C_j‖` are
gauge-INVARIANT, and so is the final `Υ'_j` (because `‖C_jξ_j‖` and the `L_j`-built congruence
absorb `u_j`). So tooth #2 and the σ_max bound (S5) are safe gauge-invariant assertions; do
NOT assert `C_j` entrywise (P3 reading).

---

## D. TEST PLAN

### D.1 η=0 oracle (the cleanest ground truth, CLAUDE.md cross-check ladder rung 3)

At η=0, `Φ` is exactly idempotent; `Δ = ι∘v` and `Υ = v^{-1}∘Φ̃` reduce to the EXACT
`th_idemp_structure` decode `Υ = Γ∘C_M` (`include/aic_idemp.h:9-12`, `Gamma_C_Delta`
`.tex:328`), and the chain collapses (`factorize_d4_research.md §3, §6`):

- **O1 (Υ = Γ∘C_M reduction).** Build Υ via the F3 pipeline on an exact-idempotent fixture
  (`block_cond_exp`, `noiseless_subsystem` — the same η=0 fixtures O2 used, `FINDINGS §C12.O2`)
  and assert it equals the `aic_idemp_decompose` `Γ∘C_M` decode to ~1e-12 (compare on a
  matrix-unit basis via `aic_mat_opnorm` of the difference; the `n²×n²` superop is never
  formed).
- **O2 (`ΥΔ = 1_B`, `ΔΥ = Φ` as maps).** Assert `‖Υ(Δ(E_i)) − E_i‖_op ≤ 1e-12` over B's
  matrix units (`.tex:2739`) and `‖Δ(Υ(E_{pq})) − Φ(E_{pq})‖_op ≤ 1e-12` over H's
  (`.tex:2733`). EXACT at η=0 (defects → 0, Rule 6).
- **O3 (`‖C_j‖ = 1` EXACT).** At η=0, `R_j = (1/d_{L_j}) 1_{L_j}⊗Tr_{L_j}(W_jW_j^†)`, `C_j`
  is a state on E_j, `‖C_j‖ = σ_max(C_j) = 1` to ~1e-15 (shard L216, `factorize_d4_research.md
  §3`). Certify in arb.
- **O4 (centrality tooth #1, #2 read machine-zero).** §C teeth at η=0.

### D.2 η>0 (mixconj family, m≥2 for the P4/centrality teeth to have teeth)

- **U1 (Υ is UCP — the §C5 midpoint+Weyl Choi-PSD route).** Assemble `C_Υ = Σ_{pq}
  E_{pq}⊗Υ(E_{pq})` (Convention A) and certify PSD via the **midpoint-collapse + Weyl-inflated
  eig** of `FINDINGS §C5 / §C12.O2-SUBTLETY`: `C_Υ` is assembled in arb over the cstar_build
  `v`, so it carries an ACCUMULATED radius (~1e-71 at prec=256) that EXCEEDS the
  `herm_max_eig` Hermiticity tol (`2^{-(prec-8)}`); feed the MIDPOINT `(C+C^†)/2` to
  `aic_mat_herm_max_eig(-C_mid)` (or the PSD-defect helper) so the genuinely-Hermitian Choi is
  not REJECTED (this is the EXACT wall O2 hit, `FINDINGS §C12.O2-SUBTLETY` — reuse the fix, do
  not rediscover it). Unitality: `‖Υ(1_H) − 1_B‖_op ≤ 1e-12` (by the S10 congruence).
- **U2 (`‖Υ − Υ̃‖_cb = O(η)`).** Υ̃ = `v^{-1}∘Φ̃` (F1's `aic_factorize_upsilon_tilde`).
  `‖Υ−Υ̃‖_cb` via the Choi-difference + diamond SDP (`aic_ucp_choi_diff` + `aic_cbnorm_*`,
  the `aic-d24`/O2 precedent). Assert `≤ C·η` with `C` bounded + dimension-independent
  (`FINDINGS §D2` robust canary: abs-max + no-upward-trend over an n-sweep, NOT a fragile
  ratio). This is the F3 piece of the composite-constant story (`factorize_d4_research.md §5`).
- **U3 (`σ_max(C_j) ≥ 1−O(η)`).** lem_RC(ii); certify in arb (S5). Mutation: shrink one Pauli
  weight (break `Σp=1`) → `σ_max(C_j)` drops below `1−O(η)` → Rule 4 abort (RED).
- **U4 (centrality teeth #1/#2/#3, §C).** On a `B=M_2⊕M_3` mixconj fixture (m=2). #3's two
  mutations must go RED.
- **U5 (the P1 partial-trace pin, as a permanent regression).** Bake `probe_f3_ptrace` (A.0)
  into the test suite as a standalone unit assertion (it needs no channel — pure kron/trace),
  AND mutation-prove the F3 path: swap `partial_trace_left`→`partial_trace_right` in S4 and
  confirm `‖R_j−1⊗C_j‖` and `σ_max(C_j)` go RED (the silent-direction tooth).

### D.3 Mutation teeth summary (Rule 7 — each must flip RED)

| mutation | where | expected RED signal |
|---|---|---|
| `partial_trace_left`→`_right` in S4 | P1 | wrong-shape/value `C_j`; `‖R_j−1⊗C_j‖` O(1); η=0 `‖C_j‖≠1` |
| single Pauli / non-central design | C tooth #3 | `diag_j2` defect O(1); `‖R_j−1⊗C_j‖` O(1) |
| joint `U_s` (F2) instead of per-block | P4 | block-mixing → `‖R_0−1⊗C_0‖` O(1) at m≥2 |
| LEFT singular vector for ξ_j | P5 | `‖C_jξ_j‖ < σ_max`; η=0 `‖C_jξ_j‖<1` |
| identity-pad instead of zero-pad `Δ(U_{js}^†)` | B.1 | `Δ(embed)≠W_j^†(U⊗1)W_j` O(1); ΥΔ≠1_B |
| skip S10 unitalization (use Υ') | S10 | `Υ(1_H)≠1_B`; not UCP (Tr_K C_Υ ≠ 1) |
| `Φ̃`(superop, not CP) for `Φ` in S8 | star §C2 | Υ'_j not CP (Choi not PSD) |

---

## E. RANKED RISKS + MITIGATIONS + FILE/LOC SPLIT

### Ranked risks

1. **(HIGHEST) P1 partial-trace direction silently wrong.** The exact bug that bit cbnorm
   and O2 twice. *Mitigation:* the A.0 GO/NO-GO probe with `d_Lj≠d_Ej`, run BEFORE any F3
   code and baked in as U5; the header convention (`aic_mat.h:55-60`) is a HYPOTHESIS to
   confirm, not trust (`§C12.O2-PIN` ethos).
2. **W_j extraction gauge + degeneracy (P3 + S0/S1).** Δ's Choi is generically degenerate, so
   `W_j` comes off the DOUBLE-path Choi→Kraus (`aic_ucp.h:46-50`), not the arb eig. The gauge
   is non-unique. *Mitigation:* never assert `C_j`/`W_j` entrywise; assert only the
   gauge-invariant `σ_max(C_j)`, `‖R_j−1⊗C_j‖`, and `Υ'_j` (P3 oracle). Certify the FINAL
   bounds in arb, accept the double-path W_j as a (gauge-arbitrary) intermediate.
3. **Centrality blindness (C).** F2's CP teeth can't see a broken design. *Mitigation:* the
   three §C teeth, m≥2 fixture, the two centrality-breaking mutations.
4. **The §C5/§C12.O2-SUBTLETY Choi-PSD arb-radius wall (U1).** `C_Υ` assembled in arb over `v`
   carries a radius exceeding the Hermiticity tol → spurious REJECT. *Mitigation:* reuse the
   PINNED midpoint-collapse+Weyl-inflation fix (`FINDINGS §C12.O2-SUBTLETY`); do NOT feed the
   raw-radius block to `herm_max_eig`.
5. **xpow basin for the two unitalizations (S10, F2's `deltaI`).** `aic_funcalc_xpow(-1/2)`
   needs `‖Υ'(1_H)−1_B‖_op<1` (`aic_funcalc.h:90`). At larger η this can fail → an explicit
   `η`-threshold (shard #7). *Mitigation:* assert the basin (Rule 4) before the power; if it
   fails, that is a genuine `η`-too-large stop condition, not a bug to paper over.
6. **σ_max(C_j) certification straddles `1−O(η)` (S5).** A precision-wall risk. *Mitigation:*
   Rule 4 abort with a clear message; escalate as a lem_RC precondition failure if it persists
   at feasible prec (CLAUDE.md stop conditions).
7. **Composite constant C unbounded / dimension-growing (U2).** The `.tex:484` failure mode.
   *Mitigation:* the `FINDINGS §D2` robust canary (abs-max + no-upward-trend); the analytic C
   is deferred (`factorize_d4_research.md §5`, posture (i)).

### File / LOC split (≤200 LOC/file, Rule 10)

| file | ~LOC | content |
|---|---|---|
| `src/aic_factorize_ups.c` (NEW) | ~180 | S0 `delta_kraus`, S1 `Wj`, S3-S4 `Rj_Cj` (with P1), S5 `xi`. The lem_RC core. |
| `src/aic_factorize_ups2.c` (NEW) | ~170 | S6 V, S7 `Lj`, S8-S9 `upsilon_prime`, S10 `upsilon_build`/`upsilon`. |
| `include/aic_factorize.h` (EDIT) | ~40 | F3 decls + handle fields (`upsI_invsqrt`, cached `L_j`/`ξ_j`, `ups_ready`); docstring (paper technique vs constructive route, the per-block-vs-joint P4 note, the P1 pin). |
| `tests/test_factorize_ups.c` (NEW) | ~190 | D.1 O1-O4 (η=0), D.2 U1-U5 (η>0), D.3 mutation teeth; the standalone `probe_f3_ptrace` regression. |

Split `ups.c`/`ups2.c` further if either exceeds 200 once written.

---

## Appendix — citation index

- partial-trace + Kronecker convention (P1/P2): `include/aic_mat.h:43-60,126-144`;
  asserted by `tests/test_mat.c:259-277` (`Tr_2(A⊗B)=Tr(B)A`, `Tr_1(A⊗B)=Tr(A)B`).
- Choi Conv A (K major), Stinespring V column-stack (F major) (P3/S1/S6):
  `include/aic_ucp.h:23-50,89-93`; `src/aic_ucp_carrier.c:41-54`.
- per-block Pauli design vs joint (P4): `.tex:2771-2784`; `src/aic_factorize_delta.c:16-90`;
  `aic_dhom_pauli` `include/aic_dhom.h`; `§A2` (non-central joint sum measured 0.54).
- lem_RC (R_j, C_j, ξ_j): `.tex:2840-2864`; shard L193-216; `factorize_d4_research.md §3`.
- diag_j2 centrality (C): `.tex:2775-2778`.
- Δ Choi/W_j (Choi_Delta): `.tex:2831-2837`; PhiDelta1/2/3 `.tex:2808-2816`; Delta_norm `.tex:2806`.
- xpow -1/2 + basin (S10): `include/aic_funcalc.h:88-105`.
- η=0 oracle Γ∘C_M (D.1): `include/aic_idemp.h:9-12,80-93`; `Gamma_C_Delta` `.tex:328`.
- §C5 Gram false-fail / §C12.O2-SUBTLETY midpoint+Weyl (U1): `paper/FINDINGS.md:142-171,606-621`.
- O2 convention-pin methodology: `docs/research/opspace_o2_design.md §0.5`; `FINDINGS §C12.O2-PIN`.
- F1/F2 reuse + handle pattern: `include/aic_factorize.h:90-182`; `src/aic_factorize_delta.c`.
- build_vinv (Υ̃): `include/aic_opspace.h:182-198`.
