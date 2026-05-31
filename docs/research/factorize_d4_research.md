# `factorize` (D4) — design / feasibility assessment of the headline th_factorization

> **Status:** DESIGN / ASSESSMENT (read-only). Bead `aic-tff` (module), gated by
> `aic-zwo` (opspace / th_main_ext, in flight) and `aic-1sk` (this research / the
> D4 closure). Realizes `th_factorization` (`approximate_algebras.tex:2730-2899`),
> `lem_RC` (`.tex:2840`). Shard: `paper/shards/shard-H-almost-idemp-factorization.md`.
> FINDINGS: §D4 (the open escalation this doc addresses), §C12/§D3 (the opspace
> interface), §C2 (the star), §C11 (master-loop coverage).
>
> **The central question (THE MANDATE, CLAUDE.md).** th_factorization's *proof*
> is an OUTLINE — the CP-ization (Steps 4-5) is prose, the composite `O(η)`
> constant is unstated. Is that a HARD WALL (genuinely non-constructive, a stop
> condition to escalate) or BUILDABLE-MODULO (constructive in finite dim despite
> the prose proof)? **Verdict, up front: BUILDABLE-MODULO.** Every object in
> Steps 4-5 is an explicit finite-dim matrix expression; the only genuinely open
> item is the *certified composite constant* (escalation 5), which is a
> certification gap, not a constructivity wall. Details below, per step.

---

## 0. Executive landscape: what is already built vs what `factorize` must add

`factorize` is the **last** module on the L7 row of the DAG (MODULE_PLAN.md). By
the time it runs, the hard machinery exists and is green:

| Pipeline ingredient | Module / symbol | State |
|---|---|---|
| `Φ̃ = θ(2Φ−1)` (superoperator) | `aic_assoc_regularize` (`src/aic_assoc_regularize.c`) | built |
| `A = Img Φ̃` + Choi-Effros star | `aic_assoc_ecstar_from_phi` (`include/aic_assoc.h`) | built |
| `(B, v)` from th_main (the iso) | `aic_cstar_build` (`include/aic_cstar.h`) | built (η>0, m≥3, multi-class — FINDINGS §C11) |
| th_main_ext cb-certificate of `v` | `aic_opspace_certify` (`include/aic_opspace.h`) | O1 built; O2 (SDP upper bd) a filed bead |
| Pauli / Heisenberg-Weyl 1-design | `aic_dhom_pauli`, `aic_dhom_diag_build` (`include/aic_dhom.h`) | built (the §12 1-design IS this object) |
| η=0 oracle (Δ=incl, Υ=Γ∘Co_M) | `aic_idemp_decompose` (`include/aic_idemp.h`) | built |
| Choi / Stinespring / compose of UCP maps | `aic_ucp_*` (`include/aic_ucp.h`) | built |

What `factorize` must **add** (none of it a new numerical primitive):

1. The non-UCP `Δ̃ = ι∘v` and `Υ̃ = v⁻¹∘Φ̃` (pure linear-algebra glue over `v`).
2. The CP-ization `Δ'` (Step 4) — a finite sum of `Φ(Δ̃(·)Δ̃(·))` over the 1-design.
3. `lem_RC`'s `R_j → C_j` and the leading-singular-vector `ξ_j` (Step 5).
4. The Kraus-like `L_j` and `Υ'_j(X) = L_j†(Φ(X)⊗1_F)L_j` (Step 5).
5. The two unitalizations `(·)^{-1/2}(·)(·)^{-1/2}` (certified inverse-sqrt).
6. The end-to-end verification `‖ΔΥ−Φ‖_cb` and `‖ΥΔ−1_B‖_cb` (the spec).

The deliverable is glue + verification, **modulo** the composite constant
(§5) and modulo opspace O2 supplying a *certified* (not just lower-bounded) cb
isometry on `v` (§4).

---

## 1. THE 7-STEP HEADLINE PIPELINE (mapped to `.tex` + modules)

Input: UCP `Φ: B(H)→B(H)`, `‖Φ²−Φ‖_cb ≤ η`, `η < 1/4`, `dim H = N` finite.

### Step 1 — Regularize. `Φ̃ = θ(2Φ−1)` (`.tex:2743`, `tilde_Phi` `.tex:2171-2174`)
`Φ̃ = ½(1 + (2Φ−1)(1−4(Φ−Φ²))^{−1/2})`, exact idempotent in the cb-algebra,
`‖Φ̃−Φ‖_cb = O(η)`, `Φ̃(I)=I`, `Φ̃(X†)=Φ̃(X)†`.
- **Constructive now?** YES. `aic_assoc_regularize` builds the `N²×N²`
  superoperator `S_Φ` and applies `aic_prop_P` (eig-free Newton-Schulz θ). The
  prop_P basin `‖S_Φ²−S_Φ‖_op < 1/4` is asserted (Rule 4). FINDINGS §C11
  (`aic-1vp`) resolved the wide-radius sgn wall — completes for the mixconj
  family up to at least n=10.
- **Module:** `assoc_ecsa` (Increment 1). **`.tex`:** 2171-2174, prop_P @524.

### Step 2 — `A = Img Φ̃`, star `Z★W = Φ̃(ZW)` (`.tex:2743`, `.tex:2184-2189`)
By `th_almost_idemp` (`.tex:2192`), `(A,★,‖·‖,†,I)` is an extended `O(η)`-C*
algebra (associativity from `Phi_assoc1/2`, C*-id from `PhiXdX`).
- **Constructive now?** YES. `aic_assoc_ecstar_from_phi` extracts the
  Frobenius-orthonormal range basis `{B_k}` of `S_tilde` (top-`dim_A` left
  singular vectors) and routes the star through the superop thunk. **The star is
  load-bearing (FINDINGS §C2):** every product in `A` is `Φ̃(XY)`, never plain.
- **Module:** `assoc_ecsa` (Increment 2). **`.tex`:** 2184-2192.

### Step 3 — `(B, v)` via th_main_ext (`.tex:2748`, th_main_ext `.tex:1538`)
`B = ⊕_{j=1}^m B(L_j)` genuine C*; `v: B→A` an **extended** `O(η)`-isomorphism
(every ampliation `1_{M_n}⊗v` is a δ-iso, same `δ=O(η)`, constant independent
of `n` and of `dim A`).
- **Constructive now?** YES for `v` itself (`aic_cstar_build`, the §9 master
  loop). The §10 "extended" upgrade does **not** recompute `v` (`.tex:1542,1557`)
  — it is a post-hoc certification (`aic_opspace_certify`) that the SAME `v` is a
  uniform cb-isometry. **Caveat:** O1 gives an operator-norm *lower* bound on the
  ampliated stretch (HOPM); the *certified* `‖v‖_cb ≤ 1+O(η)` upper bound is O2
  (Watrous cb-norm SDP, a filed bead). For `factorize`'s `tilde_DelUps` norm
  bounds (`‖Δ̃‖_cb, ‖Υ̃‖_cb ≤ 1+O(η)`) the O2 upper bound is the right input — see §4.
- **Module:** `cstar_build` (the iso) + `opspace` (the cb-certificate).
  **`.tex`:** 2748, 1538-1557.

### Step 4 — CP-ize `Δ̃` into a UCP `Δ` (`.tex:2771-2801`)
`Δ̃ = ι∘v: B→A→B(H)` (not UCP). Per block `j`, the unitary 1-design
`D_j = Σ_s p_{js} U_{js}†⊗U_{js}` (`diag_j2` centrality, `.tex:2776`) gives the
manifestly-CP averaged map
```
Δ'(X) = Σ_s p_s Φ( Δ̃(X U_s†) Δ̃(U_s) )        (.tex:2788)
```
(CP because positive `X=Y†Y` maps to `Σ_s p_s Φ_n(Δ̃_n(Y†(I⊗U_s†)) Δ̃_n((I⊗U_s)Y)) ≥ 0`,
`.tex:2791-2795`); `‖Δ'−Δ̃‖_cb = O(η)` by `tilde_Del2` (`.tex:2797`); then
unitalize `Δ(X) = Δ'(I_B)^{−1/2} Δ'(X) Δ'(I_B)^{−1/2}` (`.tex:2799`), UCP,
`‖Δ−Δ̃‖_cb = O(η)`.
- **PROSE-ONLY in the `.tex`** (this is half of the D4 gap), but every piece is
  an explicit finite matrix expression. See §2.
- **Module:** `factorize` (new), reusing `aic_dhom_diag_build` (the 1-design),
  `aic_ucp_apply` (`Φ`), `aic_funcalc_pow` (inverse sqrt). **`.tex`:** 2771-2801.

### Step 5 — Decode `Υ̃` into a UCP `Υ` (`.tex:2840-2899`)
`Υ̃ = v⁻¹∘Φ̃` (not UCP). Build via Choi `Φ(X)=V†(X⊗1_F)V` (`.tex:2859`), Choi of
`Δ` (`Choi_Delta` `.tex:2831`, `Δ(X)=Σ_j W_j†(X_j⊗1_{E_j})W_j`), `lem_RC`'s `C_j`
(`.tex:2840`), `ξ_j` = unit vector with `‖C_j ξ_j‖≈1`, then
```
L_j = Σ_s p_{js} (Δ(U_{js}†)⊗1_F) V W_j† (U_{js}⊗ξ_j)    (.tex:2862)
Υ'_j(X) = L_j† (Φ(X)⊗1_F) L_j  (manifestly CP)            (.tex:2869)
Υ(X) = Υ'(1_H)^{−1/2} Υ'(X) Υ'(1_H)^{−1/2}  (UCP)         (.tex:2897)
```
- **PROSE-ONLY in the `.tex`** (the other half of D4). All finite matrix ops;
  `lem_RC` is itself constructive (§3). See §2.
- **Module:** `factorize` (new), reusing `aic_ucp_kraus_to_stinespring` (`V`,
  `W_j`), `lem_RC` helper, `aic_funcalc_pow`. **`.tex`:** 2840-2899.

### Step 6 — Verify (the spec — CLAUDE.md Rule 5/6)
Assert the theorem's bounds:
```
‖ΔΥ − Φ‖_cb ≤ C·η            (DelUps, .tex:2733)
‖Υ_n(Δ_n(X)Δ_n(Y)) − XY‖ ≤ C·η‖X‖‖Y‖   (UpsDel2, .tex:2736)
  ⇒ ‖ΥΔ − 1_B‖_cb ≤ C·η      (.tex:2739)
```
- **Module:** `factorize` tests. `‖·‖_cb` of `ΔΥ−Φ` and `ΥΔ−1_B` is the
  diamond/cb-norm of a difference of (CP, here non-CP-difference) maps — the
  Choi-difference + diamond SDP (`aic_cbnorm_*`, bead `aic-d24` precedent). The
  η=0 oracle gives EXACT equality (defects → 0). **`.tex`:** 2733-2739.

### Step 7 — Dual channels (`.tex:2159`, Enc `.tex:280`, Dec `.tex:286`)
`Dec = Δ*` (states on B(H) → states on B), `Enc = Υ*` (states on B → states on
B(H)). `Dec·Enc = (ΥΔ)* ≈ 1_{B*}`, `Enc·Dec = (ΔΥ)* ≈ Φ*`. The approximate
factorization `Φ* ≈ Enc·Dec` up to `O(η)` in cb/diamond norm.
- **Constructive now?** YES — the dual of a UCP map is the transpose of its
  superoperator (the Schrödinger/Heisenberg flip, CLAUDE.md "UCP vs CPTP"). No
  new construction; `factorize` returns `Δ`, `Υ`, `B` and the duals are read off.
  **`.tex`:** 2159, 280-287.

---

## 2. THE D4 ASSESSMENT — per-step HARD-WALL vs BUILDABLE-MODULO

**Where the prose is.** The labelled proof block ends at `.tex:2769` with "To
complete the proof, we will approximate `Δ̃` and `Υ̃` by some UCP maps...". The
actual constructions live at `.tex:2771-2899` as *continuous prose, not flagged
as proof steps* (shard H gap #1). So Steps 1-3 are theorem-backed; **Steps 4-5
are the prose region** and are the locus of D4.

What is "missing" decomposes into three distinct items; I judge each.

### (D4-a) The unitary 1-design CP-ization (Step 4) — **BUILDABLE-MODULO (trivial)**
The 1-design is **finite, explicit, and already built.** `.tex:2780` says "(See
(Pauli_diag) for an explicit example)", and `Pauli_diag` (`.tex:1250-1254`) is
the generalized-Pauli / Heisenberg-Weyl design
`D_l = Σ_{j,k} d_l^{−2} S_{jk}†⊗S_{jk}`, `S_{jk}=X^jZ^k`. This is **exactly** the
object `aic_dhom_diag_build` already produces (`include/aic_dhom.h`, the
"dimension-independence crux"), stored as `(weight, U)` pairs `(p_{js}, U_{js})`
— literally the `Σ_s p_{js} U_{js}†⊗U_{js}` of `.tex:2771`.

- **Reuse check (shard H open-question #4):** the `aic_dhom_diag` list is per-block
  `d_l²` terms `(w = d_l^{−2}, U = Ŝ^{(l)}_{jk})`. For Step 4 we need, **per block
  separately**, the within-block design `D_j` on `B(L_j)`: `p_{js} = d_j^{−2}`,
  `U_{js} = S_{jk}` (the genuine `d_j×d_j` Pauli, **not** the embedded partial
  isometry `Ŝ`). The needed per-block primitive is `aic_dhom_pauli(out, d, j, k)`
  (already exposed) summed with weight `d^{−2}` — `d²` terms per block. So the
  1-design is reusable directly; one only must take the **genuine per-block**
  `S_{jk}` (from `aic_dhom_pauli`), not the cross-term-free embedded sum `D` that
  `aic_dhom_diag_build` assembles for the lem_approx error reduction (FINDINGS §A2
  — that embedded form is for a *different* purpose, the diagonal-of-the-whole-B
  in lem_approx; Step 4 wants the per-block design `D_j`).
- **CP-ness:** `Δ'(Y†Y) = Σ_s p_s Φ_n(Δ̃_n(Y†(I⊗U_s†))Δ̃_n((I⊗U_s)Y)) ≥ 0` because
  `Φ` is CP and the summand is `Φ` of `Z†Z` (with `Z = Δ̃_n((I⊗U_s)Y)`). This is a
  one-line verification, not an existence claim.
- **Verdict:** trivially constructive. The only "missing" thing in the `.tex` is
  the *closure statement* that `Δ'` is CP and `O(η)`-close — both are explicit
  and checkable. **NOT a wall.**

### (D4-b) `lem_RC` decode `C_j`, `ξ_j` (Step 5) — **BUILDABLE-MODULO (explicit)**
`lem_RC` (`.tex:2840`) is a *constructive lemma*, not an existence claim: see §3.
`R_j = Σ_s p_{js}(U_{js}†⊗1)W_jW_j†(U_{js}⊗1)` is an explicit finite sum;
`R_j = 1_{L_j}⊗C_j` and `C_j = d_{L_j}^{−1} Tr_{L_j}(R_j)` (partial trace over the
`L_j` factor — shard H notes this); `ξ_j` = leading right singular vector of
`C_j` (shard H #5, the SVD recipe), giving `‖C_jξ_j‖ = σ_max(C_j) ≥ 1−O(η)` by
`lem_RC`(ii). All finite, all matrix ops.
- **Verdict:** constructive. **NOT a wall.** (Certification of the `≥1−O(η)`
  lower bound needs arb — routine.)

### (D4-c) The composite `O(η)` constant — **the genuinely open item (escalation 5)**
This is NOT a constructivity wall: the *algorithm* runs and produces `Δ,Υ,B`. It
is a **certification** gap. The headline `O(η)` is the composition of:
`O(η)` from `Φ̃≈Φ` (prop_P / assoc) × `O(η)` from `th_main_ext` (the iso `v`,
itself a deferred chain through cor_improvement's `c_0`, FINDINGS §D2, still
analytically OPEN, bead `aic-1bc`) × `O(η)` from the CP-ization (Steps 4-5). The
paper writes every link as `O(η)` and never multiplies them out. See §5.
- **Verdict:** BUILDABLE-MODULO, but the *certified* end-to-end `C` is the real
  open work. The honest posture (matching opspace O1/O2): **measure** `C` per
  instance (`‖ΔΥ−Φ‖_cb / η`), assert it is bounded and dimension-independent (the
  `.tex:484` canary, FINDINGS §D2 robust form), and defer the *analytic* `C` to a
  research bead — exactly as `cstar_build` did with `c_0`.

### Summary table

| D4 item | `.tex` | Finite-dim constructive? | Verdict |
|---|---|---|---|
| (a) 1-design CP-ization of `Δ̃` | 2771-2801 | yes, design already built (`aic_dhom_pauli`) | BUILDABLE-MODULO (trivial) |
| (b) `lem_RC` / `C_j` / `ξ_j` decode | 2840-2864 | yes, partial trace + SVD | BUILDABLE-MODULO |
| (c) composite `O(η)` constant | 2733/2895 | algorithm runs; constant unstated | BUILDABLE-MODULO; certified `C` is the open work (escalation 5) |
| (closure) link constructions → DelUps/UpsDel2 | 2895 chain | the chain `Υ'≈Υ'Φ≈Υ'Φ̃=Υ'Δ̃Υ̃≈Υ'ΔΥ̃≈Υ̃` is explicit | reconstructed below |

**D4 is NOT a hard wall.** No object in Steps 4-5 is non-constructive in finite
dim; the proof-closure chain (`.tex:2895`) is a sequence of `O(η)`-triangle
inequalities each backed by an already-stated bound (`PhiDelta1`, `PhiDelta3`,
`lem_RC`, `tilde_DelUps`). The one genuinely open thing is the *value* of the
composite constant, which the project already handles per-instance + canary
elsewhere (`c_0`, opspace O1).

---

## 3. `lem_RC` (`.tex:2840`) IN DETAIL

**Statement (`.tex:2840-2846`).** Given the Choi rep of `Δ`,
`Δ(X)=Σ_j W_j†(X_j⊗1_{E_j})W_j` with `W_j: H→L_j⊗E_j`, `Σ_j W_j†W_j = 1_H`, and
the per-block 1-design `p_{js},U_{js}` on `B(L_j)`, define
```
R_j = Σ_s p_{js}(U_{js}†⊗1_{E_j}) W_j W_j† (U_{js}⊗1_{E_j})   ∈ B(L_j⊗E_j).   (R_def)
```
Then **(i)** `R_j = 1_{L_j}⊗C_j` for some `C_j ∈ B(E_j)`, and **(ii)**
`1−O(η) ≤ ‖C_j‖ ≤ 1`.

**Proof technique (`.tex:2848-2857`).**
- (i) By the 1-design centrality `diag_j2` (`.tex:2776`), `R_j` commutes with
  `X⊗1_{E_j}` for all `X∈B(L_j)`. The commutant of `B(L_j)⊗1` in `B(L_j⊗E_j)` is
  `1_{L_j}⊗B(E_j)` (Schur / double-commutant in finite dim), so `R_j = 1⊗C_j`.
- (ii) Upper: `‖C_j‖ = ‖R_j‖ ≤ Σ_s p_{js}‖W_j‖² ≤ 1` (each `U_{js}` unitary,
  `‖W_j‖≤1`). Lower: `Φ(W_j†R_jW_j) = Σ_s p_{js}Φ(Δ(U_{js}†)Δ(U_{js}))
  = Σ_s p_{js}Δ(U_{js}†U_{js}) + O(η) = Δ(1_{L_j}) + O(η)` (using `PhiDelta2`
  `.tex:2810`), and `‖Δ(1_{L_j})‖ ≥ 1−O(η)` (using `Delta_norm` `.tex:2806`);
  since `‖Φ(W_j†R_jW_j)‖ ≤ ‖C_j‖`, get `‖C_j‖ ≥ 1−O(η)`.

**Is the decode `Υ` it produces finite-dim-constructive? YES.**
- **`R_j`** is an explicit finite sum (`d_j²` Pauli terms × the Choi pieces of
  `Δ`). `W_j` is read from the Choi/Stinespring rep of the *already-built UCP* `Δ`
  (`aic_ucp_kraus_to_stinespring` on `Δ`'s Kraus form; `W_j` is the block of `V_Δ`
  landing in `L_j⊗E_j`). All available once `Δ` (Step 4) exists.
- **`C_j`** = `d_{L_j}^{−1} Tr_{L_j}(R_j)` (the commutant decomposition is computed
  by partial trace over the `L_j` factor — shard H `lem_RC` "Constructive in
  finite dim?" note; cross-check `‖R_j − 1⊗C_j‖ = O(η)` certifies (i)).
- **`ξ_j`** = top right singular vector of `C_j` (shard H #5); `‖C_jξ_j‖ = σ_max(C_j)`,
  and `lem_RC`(ii) certifies `σ_max(C_j) ≥ 1−O(η)`, so the normalization in `L_j`
  is well-conditioned. Arb certifies the lower bound (Rule 4: abort if `σ_max`
  ball straddles `1−O(η)` at feasible prec — a `lem_RC` precondition failure).
- The rest (`L_j`, `Υ'_j(X)=L_j†(Φ(X)⊗1_F)L_j`) are explicit triple products.

**η=0 oracle for `lem_RC`** (shard H): `R_j = d_{L_j}^{−1} 1_{L_j}⊗Tr_{L_j}(W_jW_j†)`,
`C_j` = a state on `E_j`, `‖C_j‖ = 1` exactly. (This dovetails with the
`th_idemp_structure` decode `Γ_j(X)=Tr_{E_j}(W_jXW_j†(1⊗γ_j))`, prop_Gamma
`.tex:2106-2113` — the `γ_j` density matrices are the η=0 limit of `C_j/Tr`.)

---

## 4. THE INTERFACE th_main_ext / opspace MUST EXPOSE TO §12 (shard H #4)

This is the **actionable payoff** of doing this assessment in parallel with
opspace (`aic-zwo`, in flight). `factorize` Step 5 needs, from `(B,v)`:

**(a) `v` at `n=1`** — the per-matrix-unit images `vE[i] = v(E_i)` of the
forward map `v: B→A ⊆ M_N`. **Already exposed:** `aic_dhom_v.vE`
(`include/aic_dhom.h`), the output `v_out` of `aic_cstar_build`. `Δ̃ = ι∘v` is
then `Δ̃(X) = Σ_i x_i vE[i]` (the ambient `N×N` matrix), and `Δ̃`'s Choi/Kraus
rep is obtained by writing `Δ̃` as a linear map `B→M_N` and Choi-izing. **No
change needed.**

**(b) the certified uniform cb-bound `a = 1−O(η)`** — `aic_opspace_result.a_cb`
and `cb_forward`. **Caveat that BITES `factorize` (flag NOW):** O1 reports `a_cb`
as an *upper* bound on the true lower-isometry `a_n` (HOPM lower-bounds
`‖v⁻¹_n‖_op`), and `cb_forward` as a *lower* bound on `‖v‖_cb`. The
`tilde_DelUps` norm bounds `‖Δ̃‖_cb, ‖Υ̃‖_cb ≤ 1+O(η)` (`.tex:2752`) are **UPPER**
bounds on `‖v‖_cb` and `‖v⁻¹‖_cb`. So **`factorize`'s `tilde_DelUps` certificate
needs opspace O2's SDP UPPER bound, not O1's HOPM lower bound.** O1 alone can
certify the η=0 complete-isometry oracle and the universality canary, but the
`1+O(η)` *upper* bound on `Δ̃,Υ̃` is gated on O2 (the Watrous cb-norm SDP, the
filed bead). This is consistent with FINDINGS §C12/§D3 and should be stated in
the `factorize` bead's dependency note.

**(c) the Pauli diagonal `D`** — `aic_dhom_diag_build` (`include/aic_dhom.h`).
But see §2 (D4-a): Step 4 wants the **genuine per-block** `D_j` on `B(L_j)`
(weights `d_j^{−2}`, unitaries the `d_j×d_j` `S_{jk}` from `aic_dhom_pauli`), NOT
the cross-term-free embedded sum `D = Σ_l D_l` that `aic_dhom_diag_build`
assembles for lem_approx (FINDINGS §A2). The block dims `d_j` come from `B`
(`aic_dhom_B.d`, the `aic_cstar_build` output). So the needed primitive
(`aic_dhom_pauli`) is exposed; `factorize` just iterates it per block.

### What opspace's output (O1/O2) should ADD NOW so §12 plugs in cleanly

This is the concrete, forward-looking recommendation (do it while opspace is
open, cheap now, expensive to retrofit):

1. **Expose `n_B = Σ_l d_l` and the block dims on the result** (or guarantee they
   are read off `v->B`). `factorize` needs the block structure of `B` and the
   Smith level `n_B` to size the per-block 1-designs and the `W_j` Choi blocks.
   `aic_opspace_result.N_max` is the forward Smith level `N`; ADD the inverse
   Smith level `n_B` explicitly (currently implicit). Low cost.
2. **Make O2's certified upper bound `‖v‖_cb`, `‖v⁻¹‖_cb` first-class result
   fields** (not just an internal SDP value), because `factorize`'s `tilde_DelUps`
   reads them directly. Name them e.g. `cb_forward_ub`, `cb_inverse_ub` so the
   O1 lower bounds (`cb_forward`) and the O2 upper bounds coexist (the
   "HOPM ≤ SDP" cross-check, FINDINGS §D3, and the actual certificate for §12).
3. **Expose `v⁻¹` (the inverse coordinate map) as a reusable object.** opspace
   O1 already inverts the coordinate matrix `M_1[k,i]=⟨B_k,vE[i]⟩_F` to build
   `v⁻¹` for the inverse stretch (`aic_opspace_inverse_stretch`). `Υ̃ = v⁻¹∘Φ̃`
   needs exactly this `v⁻¹` (as `dim_A` operators `v⁻¹(B_k) ∈ B`). Currently it
   is internal to `aic_opspace_inverse_stretch`; lifting it to a small public
   builder (`aic_opspace_build_vinv` returning the `dim_A` B-side images) saves
   `factorize` from re-deriving it and guarantees the SAME inverse is certified
   and used. **This is the single highest-value add.**
4. **Carry the measured `delta = iso_def` through.** `factorize`'s composite
   constant (§5) multiplies `iso_def` (the `v` defect) by the assoc `eps` and the
   CP-ization constant; opspace already takes `delta` as input — surface it on the
   result so `factorize` composes from one struct.

---

## 5. THE COMPOSITE `O(η)` CONSTANT (escalation 5)

**Can the per-lemma constants be composed into a certified
`‖ΔΥ−Φ‖_cb ≤ C·η`? In principle yes; the analytic `C` is unstated and is the
open work.** Sketch of the chain (each link an `O(η)` triangle step):

```
Φ  --(prop_P/assoc, .tex:2179)-->  Φ̃            ‖Φ̃−Φ‖_cb ≤ c1·η
Φ̃ = Δ̃ Υ̃        (exact, tilde_DelUps .tex:2751)
Δ̃ --(Step 4, .tex:2797,2801)-->    Δ            ‖Δ−Δ̃‖_cb ≤ c2·η   (uses tilde_Del2)
Υ̃ --(Step 5, .tex:2895,2899)-->    Υ            ‖Υ−Υ̃‖_cb ≤ c3·η   (uses PhiDelta1/3, lem_RC)
```
Then, with `Δ = Δ̃ + O(c2 η)`, `Υ = Υ̃ + O(c3 η)`, `Δ̃Υ̃ = Φ̃ = Φ + O(c1 η)`:
```
‖ΔΥ − Φ‖_cb ≤ ‖ΔΥ − Δ̃Υ̃‖_cb + ‖Δ̃Υ̃ − Φ‖_cb
            ≤ ‖Δ‖_cb‖Υ−Υ̃‖_cb + ‖Δ−Δ̃‖_cb‖Υ̃‖_cb + ‖Φ̃−Φ‖_cb
            ≤ (1+O(η))c3 η + c2 η(1+O(η)) + c1 η  =  C·η,
   C = c1 + c2 + c3 + O(η).
```
and similarly `‖ΥΔ−1_B‖_cb ≤ C'·η` from `UpsDel2`'s chain
`Υ'≈Υ'Φ≈Υ'Φ̃=Υ'Δ̃Υ̃≈Υ'ΔΥ̃≈Υ̃` (`.tex:2895`).

**Where the constant is unstated (the gaps):**
- `c1` (prop_P): from the `θ` Taylor series remainder for `η<1/4`. Bounded
  per-instance by the assoc layer (`aic_assoc_regularize` measures `‖Φ̃−Φ‖`); the
  *analytic* constant is unstated in the paper.
- The `v`-isomorphism constant `c0` inside `c2, c3` — th_main_ext's `O(η)` is
  itself a deferred chain through cor_improvement's `c_0`, **analytically OPEN**
  (FINDINGS §D2, bead `aic-1bc`; the master loop fixes `c_0` to the MEASURED
  first-call value). So `c2, c3` inherit this open constant.
- `c2, c3` (the CP-ization): the `O(η)` in `‖Δ'−Δ̃‖_cb` (from `tilde_Del2`) and
  in the `Υ'` chain are stated only as `O(η)` in the prose; never multiplied out.
- The two unitalizations `(Δ'(I_B))^{−1/2}`, `(Υ'(1_H))^{−1/2}`: invertibility
  needs `‖Δ'(I_B)−1_H‖, ‖Υ'(1_H)−1_B‖ = O(η) < 1` (shard H #7), giving an
  explicit `η`-threshold `η < 1/C` — derivable, not in the paper.

**Recommended posture (matches `c_0` and opspace O1).** `factorize` **measures**
`C = ‖ΔΥ−Φ‖_cb / η` per instance (the diamond-SDP on the Choi difference),
asserts it is bounded (e.g. `< C_max`, a generous fail-loud gate) and
**dimension-independent** (the FINDINGS §D2 robust canary: bounded abs-max +
no-upward-trend over an `n`-sweep, NOT a fragile within-family ratio). The
*analytic* `C` is a follow-up research bead (the analogue of `aic-1bc` for the
factorization chain). This is honest: the algorithm is certified to meet `O(η)`
*per instance with a measured, non-growing constant*, which is the project's
standard for the other deferred constants.

---

## 6. RECOMMENDATION + escalation

**Is `factorize` reachable after opspace O2? YES — BUILDABLE-MODULO.** D4 is not
a hard wall. Concretely:

- After **opspace O1** (in flight): the η=0 oracle path is fully reachable
  (`Δ = ι∘v`, `Υ = v⁻¹∘Φ̃` reduce to the exact `th_idemp_structure` `Δ, Γ∘Co_M`,
  `ΔΥ=Φ` and `ΥΔ=1_B` EXACTLY — the cleanest cross-check, Rule 6). The η>0 path
  *runs* and *measures* the composite constant.
- After **opspace O2** (the SDP upper bound, filed bead): the η>0 `tilde_DelUps`
  norm bounds (`‖Δ̃‖_cb,‖Υ̃‖_cb ≤ 1+O(η)`) are *certified*, and the end-to-end
  `‖ΔΥ−Φ‖_cb` certificate is trustworthy (the diamond SDP gives a certified
  upper bound on the cb-difference).

So `factorize` should be **gated on opspace O2** for the certified η>0 headline,
but the η=0 oracle and the η>0 *measured* path are reachable after O1.

### Increment skeleton (BUILDABLE-MODULO)

`factorize` is glue + verification; suggest 4 increments (each ≤ ~200 LOC, the
file-split rule), reusing closed modules, each with its η=0 oracle + an η>0
mixconj fixture + the star/non-vacuity teeth:

- **F1 — `Δ̃`, `Υ̃` (non-UCP) + `tilde_DelUps` certificate.** Build `Δ̃ = ι∘v`
  (read `vE` from `aic_cstar_build`) and `Υ̃ = v⁻¹∘Φ̃` (read `v⁻¹` from the
  opspace add #3). Certify `Δ̃Υ̃ = Φ̃` and `Υ̃Δ̃ = 1_B` (exact by construction),
  and the `‖Δ̃‖_cb,‖Υ̃‖_cb ≤ 1+O(η)` bounds from opspace O2. **`.tex`:** 2749-2767.
- **F2 — `Δ` (Step 4).** Per-block `D_j` (genuine Pauli via `aic_dhom_pauli`),
  `Δ'(X) = Σ_s p_s Φ(Δ̃(XU_s†)Δ̃(U_s))`, CP-check via Choi PSD
  (`aic_ucp_is_cp_choi`), unitalize with `aic_funcalc_pow` `(·)^{−1/2}`. Assert
  `‖Δ'(I_B)−1_H‖ = O(η) < 1` (invertibility, shard H #7). Certify
  `‖Δ−Δ̃‖_cb = O(η)` and the `PhiDelta1/2/3` properties. **`.tex`:** 2771-2829.
- **F3 — `lem_RC` + `Υ` (Step 5).** `R_j → C_j` (partial trace) → `ξ_j` (top
  SVD), assert `σ_max(C_j) ≥ 1−O(η)` (Rule 4). Build `L_j`, `Υ'_j`, unitalize.
  Certify `‖Υ−Υ̃‖_cb = O(η)`. **`.tex`:** 2840-2899.
- **F4 — end-to-end DelUps/UpsDel2 + duals + the composite-constant canary.**
  `‖ΔΥ−Φ‖_cb`, `‖ΥΔ−1_B‖_cb` via the diamond SDP on Choi differences; the
  measured `C = ‖ΔΥ−Φ‖_cb/η` bounded + dimension-independent (FINDINGS §D2
  robust canary); `Enc=Υ*`, `Dec=Δ*` read off. **`.tex`:** 2733-2739, 2159.

### Escalation (stop-condition format — the ONE genuinely open item)

> **Result:** `th_factorization` composite `O(η)` constant `C` (escalation 5).
> **`.tex` line:** `.tex:2733` (DelUps), `.tex:2895` (the closure chain); the
> per-link `O(η)` at `.tex:2179` (prop_P), `.tex:2797` (Δ'), `.tex:2895` (Υ').
> **Obstruction:** the paper writes every link as `O(η)` and never composes the
> constants; moreover one factor (th_main_ext's iso constant, via cor_improvement
> `c_0`) is itself analytically OPEN (FINDINGS §D2, bead `aic-1bc`). So a
> *certified analytic* `C` is not extractable from the `.tex` as written.
> **What a decision/repro needs:** EITHER (i) accept the project-standard posture
> — measure `C` per instance + assert bounded & dimension-independent (the
> `cstar_build`/opspace precedent), filing the analytic `C` as a research bead
> chained after `aic-1bc`; OR (ii) escalate to the user a request to derive the
> three constants (`c1` from the prop_P remainder, `c2,c3` from the CP-ization
> `O(η)` prose, the unitalization `η`-threshold) analytically. **Recommendation:
> posture (i)** — it is constructive, certified per-instance, and matches how the
> project already treats every other deferred universal constant. The CP-ization
> *constructions* themselves (D4-a, D4-b) are NOT a wall and need no escalation.

**Net:** D4 is BUILDABLE-MODULO, not a hard wall. The CP-ization (Steps 4-5) is
fully constructive in finite dim from objects already built; the only open item
is the *value* of the composite constant, which the project handles
per-instance + canary. The actionable opspace add (do it now): lift `v⁻¹` to a
public builder, make O2's certified cb *upper* bounds first-class result fields,
and expose the inverse Smith level `n_B`.
