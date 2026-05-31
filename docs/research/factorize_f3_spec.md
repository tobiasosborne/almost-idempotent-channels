# `factorize` F3 — DESIGN / DERIVATION SPEC: the UCP decode map Υ via lem_RC

> **Status:** DESIGN / DERIVATION (read-only on code; this doc is the only
> artifact). Bead `aic-tff` (module `factorize`), increment **F3** = Step 5 of
> `th_factorization` (`approximate_algebras.tex:2730`): the UCP **decode** map
> `Υ: B(H) → B`, built from `lem_RC` (`.tex:2840`). Gated on F2 (the UCP encode
> `Δ`, both `aic_factorize_delta*.c`) and on the closed modules `ucp`, `mat`,
> `dhom`, `funcalc`, `assoc`, `cstar_build`, `opspace`, `idemp`.
>
> **Read order honored:** `CLAUDE.md`; `.tex:2803-2899` (the verbatim Step-5
> construction); `shard-H-almost-idemp-factorization.md` (lem_RC); D4 research
> §3 (lem_RC in detail); `FINDINGS.md` §A2, §C2/C3/C5, §C12/D3, §D4.
>
> **One-line verdict (inherited from D4, FINDINGS §D4):** F3 is
> BUILDABLE-MODULO. Every object in Step 5 is an explicit finite-dim matrix
> expression. The single open item is the *value* of the composite `O(η)`
> constant (deferred to F4 / `aic-1bc`), not a constructivity wall. No new
> stop-condition is expected.

---

## 0. PRELIMINARY — a working-tree caveat the F3 implementer must clear first

F2's sources `src/aic_factorize_delta.c`, `src/aic_factorize_delta2.c` and the
F2 deltas to `include/aic_factorize.h` / `src/aic_factorize.c` are **untracked /
modified** in the working tree (not yet committed). More importantly,
`aic_factorize_delta_prime` (`src/aic_factorize_delta.c:148`) currently contains
a **live mutation-test artifact**:

```c
for (slong s = 0; s < 1; s++) {   /* HUNT1-MUT-B: s=0 only, weight unchanged */
```

i.e. the 1-design average runs **only the `s=0` term** (the identity Pauli),
which is NOT the full design `Σ_s p_s (...)`. F3 reads `Δ` heavily (Step 5 uses
`Δ(U_{js}†)`, `Δ(I_B)` via Choi, and the 1-design). **Before building F3,
restore the loop to `s < aic_factorize_design_nterms(F)`** (the mutation was
left mid-proof). Treat this as an F2 close-out item; F3's lem_RC(i) teeth
(`‖R_j − 1⊗C_j‖ ≈ 0`, §B below) will themselves *fail loud* on the truncated
design, because lem_RC(i) is exactly where the FULL design's centrality
(`diag_j2`) is load-bearing — so a not-restored mutation surfaces in F3 as a
RED centrality tooth, not as silent garbage. Note this in the bead.

---

## 1. Where F3 sits in Step 5 and what it inherits

`th_factorization` Step 5 (`.tex:2803-2899`) builds, in order:

1. The Choi rep of the UCP `Δ` (Step 4's output): `Δ(X)=Σ_j W_j†(X_j⊗1_{E_j})W_j`
   (`Choi_Delta`, `.tex:2831-2837`).
2. `lem_RC` (`.tex:2840`): `R_j = Σ_s p_{js}(U_{js}†⊗1)W_jW_j†(U_{js}⊗1)`,
   `R_j = 1_{L_j}⊗C_j`, `1−O(η) ≤ ‖C_j‖ ≤ 1`.
3. `ξ_j` (`.tex:2859`): unit vector with `‖C_jξ_j‖ ≈ 1`.
4. `L_j` (`.tex:2860-2864`): `L_j = Σ_s p_{js}(Δ(U_{js}†)⊗1_F) V W_j†(U_{js}⊗ξ_j)`.
5. `Υ'_j(X) = L_j†(Φ(X)⊗1_F)L_j` (`.tex:2866-2870`), manifestly CP, `Υ'=(Υ'_1,…,Υ'_m)`.
6. Unitalize `Υ(X)=Υ'(1_H)^{−1/2}Υ'(X)Υ'(1_H)^{−1/2}` (`.tex:2896-2899`).

Inherited from F1/F2 (all built):
- `F->v` (B→A, `aic_cstar_build`), `F->Aec` (A=Img Φ̃ + `S_tilde`),
  `F->phi` (the **original** UCP Φ), `F->vinvB[k]` (= v⁻¹(B_k)).
- `F->delta_ready`, `F->deltaI_invsqrt = Δ'(I_B)^{−1/2}`; the UCP `Δ` itself via
  `aic_factorize_delta(out, F, X, prec)`.
- The per-block 1-design: `aic_factorize_design_unitary(out,F,s,prec)` (the JOINT
  `U_s`) and `aic_factorize_design_nterms`. **CAUTION:** the joint `U_s` is the
  whole-B block-diagonal join (`.tex:2782`); lem_RC needs the **per-block**
  Paulis `U_{js}` on `L_j = C^{d_j}` — see §A.2. Get those directly from
  `aic_dhom_pauli(out, d_j, a, b, prec)` (`include/aic_dhom.h:166`).

The properties of `Δ` that F3 *uses* (not re-derives) are `Delta_norm`
(`.tex:2806`), `PhiDelta1` (`.tex:2809`), `PhiDelta2` (`.tex:2811`), `PhiDelta3`
(`.tex:2814`) — they back lem_RC(ii) and the `Υ'Δ≈1_B` calc.

---

## 2. THE PARTIAL-TRACE DIRECTION — the #1 risk, pinned here FRONT AND CENTER

This is the FINDINGS bug-class (§C12/D3, §C12.O2-PIN): a partial-trace /
Choi-convention direction is the recurring "test that can't fail" landmine. F3
has TWO direction choices — the `L_j ⊗ E_j` factor ordering inside `R_j`, and
the `K ⊗ H` factor ordering inside the Choi/Stinespring of `Δ`. Both are pinned
against the project conventions below; the cross-check `‖R_j − 1_{L_j}⊗C_j‖` and
an asymmetric fixture are the *teeth* that catch a wrong direction.

### 2.1 The conventions, verbatim from the headers

**`aic_mat.h` tensor convention (`include/aic_mat.h:43-60`, pinned, load-bearing).**
A matrix on `C^a ⊗ C^b` is stored `(a·b)×(a·b)` with the **LEFT factor MAJOR**
(row-major Kronecker, NumPy/textbook `kron`):

```
row I = i·b + j   (i ∈ [0,a) = LEFT, j ∈ [0,b) = RIGHT)
M[I,J] = <i,j| M |k,l>
Tr_2 M  ("trace out RIGHT/C^b factor, keep LEFT C^a") = aic_mat_partial_trace_right
Tr_1 M  ("trace out LEFT/C^a factor,  keep RIGHT C^b") = aic_mat_partial_trace_left
```

with `Tr_2(A⊗B)=Tr(B)·A` and `Tr_1(A⊗B)=Tr(A)·B` (tested in `test_mat.c`).

**`aic_ucp.h` Choi Convention A (`include/aic_ucp.h:23-42`).** For a UCP map
`Φ: B(K)→B(H)` (Kraus `K_a: H→K`, Heisenberg `Φ(X)=Σ_a K_a† X K_a`), the Choi is
`C_Φ = Σ_{ij} E_{ij} ⊗ Φ(E_{ij}) ∈ B(K)⊗B(H)`, **K factor LEFT/MAJOR, H factor
RIGHT**, size `(dim_K·dim_H)²`. The Stinespring `V: H → K⊗F` is the column-stack
`V[a·dim_K + i, j] = K_a[i,j]` (`include/aic_ucp.h:89-93`), `V†V = 1_H`.

### 2.2 The `R_j` direction VERDICT (the load-bearing line)

`R_j ∈ B(L_j ⊗ E_j)` (`.tex:2842`). In `Δ`'s Choi `Δ(X)=Σ_j W_j†(X_j⊗1_{E_j})W_j`,
`W_j: H → L_j ⊗ E_j` (`.tex:2837`). The block algebra acts on the **L_j factor**
(`X_j ⊗ 1_{E_j}`: `X_j` on `L_j`, identity on `E_j`), so in the `L_j ⊗ E_j`
storage **L_j is the LEFT/MAJOR factor, E_j is the RIGHT/MINOR factor** — to match
the `X_j ⊗ 1` Kronecker, `L_j` MUST be left (so that `X_j⊗1_{E_j}` built with
`aic_mat_kronecker(., X_j, I_{e_j})` is the correct embedding).

lem_RC(i): `R_j = 1_{L_j} ⊗ C_j`, `C_j ∈ B(E_j)`. To extract `C_j` we trace out
`L_j` (the LEFT factor) and KEEP `E_j` (the RIGHT factor). The shard's recipe is
`C_j = (1/dim L_j) Tr_{L_j}(R_j)` (shard H lem_RC "Constructive in finite dim?",
D4 §3). Tracing out the LEFT factor and keeping the RIGHT is exactly:

> **VERDICT: use `aic_mat_partial_trace_left(Cj, Rj, d_Lj, e_j, prec)`**
> (`include/aic_mat.h:143`, `Tr_1`, traces out LEFT `C^a`, keeps RIGHT `C^b`),
> with `a = d_Lj = dim L_j = d_j` and `b = e_j = dim E_j`, then **scale by
> `1/d_Lj`**: `C_j = (1/d_j) · Tr_1(R_j)`.

Sanity in the convention: if `R_j = 1_{L_j} ⊗ C_j` exactly (the η=0 oracle),
then `Tr_1(1_{L_j} ⊗ C_j) = Tr(1_{L_j})·C_j = d_j·C_j`, so `(1/d_j)Tr_1(R_j) =
C_j` exactly — confirms both the direction AND the `1/d_j` normalization.

**The trap if reversed:** `aic_mat_partial_trace_right` (`Tr_2`, keeps LEFT) would
return `(1/e_j)Tr_2(R_j)`. On `R_j = 1_{L_j}⊗C_j` that gives `Tr(C_j)·1_{L_j}/e_j`
— a `d_j×d_j` operator (WRONG shape: should be `e_j×e_j`), so a shape assert
catches the *gross* reversal immediately. But the subtle reversal hides if one
ALSO swaps the factor ordering (storing `E_j` left, `L_j` right); then the shapes
match and only an asymmetric fixture exposes it. Hence §2.3.

### 2.3 The empirical PIN (an asymmetric fixture; FINDINGS §C12.O2-PIN ethos)

Conventions are PINNED not derived (FINDINGS §C12.O2-PIN: "neither the design's
nor the LLM-research-leg's derivation was trustworthy; only the
independent-oracle measurement was"). Build a fixture where the WRONG direction
gives a VISIBLY wrong `C_j`:

- Take `d_Lj = 2`, `e_j = 3` (so the two factors have **different dims** — a
  symmetric `d=e` fixture is BLIND to a factor swap). Construct
  `R_test = 1_{L_j} ⊗ C0` with an **asymmetric, non-scalar** `C0 = diag(1, 0.5,
  0.2)` (Hermitian, distinct eigenvalues, `‖C0‖=1`). Build `R_test` via
  `aic_mat_kronecker(R_test, I_2, C0)`.
- `C_correct = (1/2) aic_mat_partial_trace_left(R_test, 2, 3)` must equal `C0`
  to machine precision (3×3). **`C_wrong = (1/3) aic_mat_partial_trace_right(R_test,
  2, 3)`** is `2×2` (shape assert RED if naively assigned to a 3×3 `Cj`) and,
  forced into shape, equals `Tr(C0)/3 · 1_2 = (1.7/3)·1_2 = 0.567·1_2` — a
  visibly wrong scalar, distinct eigenstructure. This is the pin.
- Add the *factor-swap* pin: build `R_swap = aic_mat_kronecker(C0, I_2)` (E left,
  L right) and confirm `partial_trace_left(R_swap, 3, 2)` gives `Tr(C0)/?·1_2`
  garbage while `partial_trace_right(R_swap, 3, 2)` recovers `C0`. This documents
  that the ORDERING (`L_j` left) is the load-bearing choice, not just the trace
  routine name.

Record the verdict + the fixture numbers in the F3 source docstring and cite
`include/aic_mat.h:43-60` and this §2.

---

## 3. THE CONSTRUCTION, STEP BY STEP — formulas, shapes, existing APIs

Throughout: `N = dim H = F->N`; `m = F->v->B->num_blocks`; `d_j = F->v->B->d[j]`
(= `dim L_j`); `n_B = F->n_B`; `dim_A = F->dim_A`. `Φ` = `F->phi` (original UCP,
`dim_K = dim_H = N`, `r` Kraus ops). The Choi of `Δ` introduces ancillas
`E_j`; the Choi of `Φ` introduces `F = C^r`.

### 3.1 The Choi / Stinespring rep of the UCP `Δ` (`.tex:2831-2837`) — THE HARD PART

`Δ: B = ⊕_j M_{d_j} → B(H)`, UCP (Step 4). We need
`Δ(X) = Σ_j W_j†(X_j⊗1_{E_j})W_j`, `W_j: H → L_j⊗E_j`, `Σ_j W_j†W_j = 1_H`.

**The structure (prop_hom_structure / prop_KLHG context).** `Δ` is UCP from a
*direct-sum* C* algebra `B`. The Choi/Stinespring of such a map block-decomposes:
the input matrix-unit `X_j` lives in `M_{d_j}` and the `W_j` is the Stinespring
isometry of `Δ` *restricted to block `j`*. The ancilla `E_j` collects Δ's
multiplicity/leakage for block `j`; `dim E_j = e_j` is the Stinespring rank of
`Δ|_{M_{d_j}}`.

**Two routes evaluated:**

- **Route (a) — full Stinespring of `Δ`, then identify blocks.** Build `Δ`'s
  whole Choi `C_Δ ∈ B(B)⊗B(H)` (Convention A, source `B`, target `H`); extract a
  Kraus set via `aic_ucp_choi_to_kraus_latd`; stack to `V_Δ: H → B-rep ⊗ F_Δ`.
  Then split `V_Δ` into per-block pieces landing in `L_j ⊗ E_j`. **Problem:** the
  source of `Δ` is the *direct-sum* `B` (block-diagonal `n_B`-rep), not a single
  `M_d`. The Convention-A Choi of a map out of `B` mixes the blocks in the
  `B`-factor index in a way that requires re-block-diagonalizing `V_Δ`; the
  ancilla `F_Δ` is shared across blocks and must be split into the `E_j`. This
  is doable but the block identification is fiddly and gauge-sensitive.

- **Route (b) — PER-BLOCK Stinespring (RECOMMENDED).** For each `j`, `Δ`
  restricted to block `j` is a CP map `Δ_j: M_{d_j} → B(H)`,
  `Δ_j(Y) = Δ(ι_j(Y))` where `ι_j: M_{d_j} ↪ B` embeds `Y` as the `j`-th block
  (zeros elsewhere). Build `Δ_j`'s Choi `C_{Δ_j} ∈ B(L_j)⊗B(H)` directly:
  ```
  C_{Δ_j}[a·N + p, b·N + q] = Δ(E^{(j)}_{ab})[p,q],   a,b ∈ [0,d_j), p,q ∈ [0,N)
  ```
  (Convention A: `L_j` factor LEFT/MAJOR, `H` factor RIGHT — IDENTICAL to F2's
  `aic_factorize_delta_block_choi`, `src/aic_factorize_delta.c:166`, except using
  the UCP `Δ` (`aic_factorize_delta`) not `Δ'`). Extract Kraus
  `{D_{j,c}: H → L_j}` via `aic_ucp_choi_to_kraus_latd` (`c ∈ [0,e_j)`,
  `e_j` = recovered rank). Then **`W_j` is the column-stack of those Kraus**:
  ```
  W_j: H → L_j ⊗ E_j,    W_j[c·d_j + a, p] = D_{j,c}[a, p],   E_j = C^{e_j},
  ```
  i.e. `W_j = aic_ucp_kraus_to_stinespring(V, Δ_j_kraus)` where the resulting
  `V` has shape `(d_j · e_j) × N`. **Verify** `W_j†(X_j⊗1_{E_j})W_j = Δ(ι_j(X_j))`
  on matrix units (the Choi-rep reconstruction tooth).

**RECOMMENDATION: Route (b).** Reasons:
- It reuses F2's `aic_factorize_delta_block_choi` pattern verbatim (same
  Convention-A `(d_j·N)` block-Choi layout), swapping `Δ'`→`Δ`. No new
  block-identification logic.
- The per-block ancilla `E_j` falls out naturally as the Stinespring rank of
  `Δ_j` (= `aic_ucp_choi_to_kraus_latd`'s recovered `r`), with `dim E_j = e_j`.
  Route (a) must *split* a shared `F_Δ`.
- The unitality `Σ_j W_j†W_j = 1_H` is `Σ_j Δ_j(1_{L_j}) = Δ(I_B) = 1_H` (Δ
  UCP) — a single check, not a gauge-fragile reassembly.

**Stinespring factor ordering for `W_j`.** `aic_ucp_kraus_to_stinespring` packs
`V[c·dim_K + a, p] = K_c[a,p]` with `dim_K = d_j` (the codomain of `Δ_j`'s Kraus
`D_{j,c}: H → L_j`). So the row index `c·d_j + a` has the ANCILLA `c` (∈ `E_j`)
MAJOR and `L_j` index `a` MINOR — i.e. **the column-stack convention puts `E_j`
LEFT, `L_j` RIGHT** in `L_j⊗E_j`... but `.tex` writes `W_j: H → L_j ⊗ E_j` with
`L_j` first. **This is a real ordering subtlety — pin it (§3.2).**

### 3.2 The W_j factor-order subtlety (a second direction pin)

`aic_ucp_kraus_to_stinespring` produces `V` with the ANCILLA index major
(`include/aic_ucp.h:89-93`: `V[a*dim_K + i, j]`, where there `a` = ancilla, `i` =
`K`-index). For `Δ_j`'s Kraus `D_{j,c}: H → L_j`, `dim_K = d_j`, so the stacked
`V_j[c·d_j + a, p]` has ancilla `c` MAJOR, `L_j`-index `a` MINOR. In the textbook
`kron` convention (`aic_mat.h`, LEFT major), that storage represents
`E_j ⊗ L_j` (ancilla left). The paper writes `L_j ⊗ E_j`.

**Two equivalent fixes; pick (i):**
- **(i) Build `W_j` directly with the `.tex` ordering (`L_j` left).** Do NOT use
  `aic_ucp_kraus_to_stinespring`; instead stack with `L_j` MAJOR:
  `W_j[a·e_j + c, p] = D_{j,c}[a, p]` (a ∈ [0,d_j), c ∈ [0,e_j)). Then
  `R_j ∈ B(L_j⊗E_j)` is in the `L_j`-major convention and §2.2 (trace out LEFT)
  applies verbatim. **RECOMMENDED** — it makes the factor ordering match the
  paper and §2.2, so the partial-trace direction analysis is airtight.
- (ii) Use `aic_ucp_kraus_to_stinespring` (ancilla-major = `E_j⊗L_j`) and then
  trace out the RIGHT factor (`L_j`) in `R_j` via `aic_mat_partial_trace_right`.
  This is internally consistent but INVERTS the §2.2 verdict, so it is
  error-prone — REJECTED to keep a single direction story.

Pin in the source: `W_j` is `L_j`-major (`W_j[a·e_j+c, p]`), `R_j` is
`L_j ⊗ E_j` with `L_j` LEFT, `C_j = (1/d_j)·partial_trace_left(R_j)`. One
ordering, one trace direction.

### 3.3 The Stinespring of the ORIGINAL `Φ` (`.tex:2859`)

`Φ(X) = V†(X⊗1_F)V`, `V: H → H⊗F`, `F = C^r`. Build with
`aic_ucp_kraus_to_stinespring(V, F->phi, prec)` → `V` shape `(N·r) × N`. Here
`dim_K = dim_H = N` (Φ self-map), so `V[a·N + i, j] = K_a[i,j]`, ancilla `a` ∈ F
MAJOR, `H`-index `i` MINOR — i.e. `F ⊗ H` in the storage. The paper writes
`X ⊗ 1_F` (X on H left, identity on F right) inside `Υ'_j` (§3.6). **So for `Φ`
the storage has F left, H right.** This is the OPPOSITE ordering from §3.2's
`W_j` choice, and that is FINE — they are independent tensor spaces. The point is
to be consistent *within each* expression. We handle this in §3.5/§3.6 by
building the `⊗1_F` embeddings explicitly with `aic_mat_kronecker` in the SAME
ordering `V` uses (F left). Concretely: anything multiplied by `V` on the
H-side must be embedded as `1_F ⊗ (H-operator)` to match `V`'s `(F major,
H minor)` row layout. **Flag this explicitly in the source; it is the third
direction pin.**

> **Cleaner alternative (RECOMMENDED to avoid the F-ordering bookkeeping):**
> do NOT form the `⊗1_F` Kroneckers at all. Use the Kraus form of `Φ` directly:
> `Φ(X)⊗1_F` appears only inside `Υ'_j(X) = L_j†(Φ(X)⊗1_F)L_j` and inside `L_j`
> via `V W_j†(...)`. Both can be rewritten in Kraus/Stinespring-free form (§3.5,
> §3.6) so that the only tensor object built explicitly is `R_j ∈ B(L_j⊗E_j)`
> (§3.4), whose ordering is pinned in §3.2. This removes two of the three
> direction pins. See §3.5/§3.6 for the rewritten contractions.

### 3.4 `R_j` and lem_RC (`.tex:2842-2846`, D4 §3)

Per block `j`, the per-block 1-design has `d_j²` Paulis: `U_{js} = S_{ab} = X^a
Z^b` (`a,b ∈ [0,d_j)`), `p_{js} = d_j^{-2}` (uniform). Get each via
`aic_dhom_pauli(S, d_j, a, b, prec)` (NOT the joint `aic_factorize_design_unitary`,
which gives the whole-B block-diagonal join — §A.2).

```
R_j = Σ_{a,b} d_j^{-2} (S_{ab}† ⊗ 1_{E_j}) W_j W_j† (S_{ab} ⊗ 1_{E_j})   ∈ B(L_j⊗E_j)
```
Shapes: `W_j W_j†` is `(d_j·e_j) × (d_j·e_j)`; `S_{ab} ⊗ 1_{E_j}` is the same,
built via `aic_mat_kronecker(SkE, S_ab, I_{e_j})` (S on `L_j` LEFT, identity on
`E_j` RIGHT — matches §3.2's `L_j`-major `W_j`). Accumulate the `d_j²`-term sum
with weight `d_j^{-2}`.

Then (§2.2):
```
C_j = (1/d_j) · aic_mat_partial_trace_left(R_j, d_j, e_j, prec)   ∈ B(E_j), e_j×e_j.
```

**lem_RC(ii) precondition (Rule 4, fail loud).** Assert `σ_max(C_j) ≥ 1 − O(η)`.
Route: `σ_max(C_j)² = λ_max(C_j† C_j)` via `aic_mat_herm_max_eig` on the (Hermitian
PSD) Gram `C_j† C_j` (degeneracy-robust). If the certified ball straddles the
`1 − c·η` threshold at feasible prec, ABORT (lem_RC precondition failure — a
straddling ball is a loud failure, not a silent pass; cf. §C5 / Rule 4). Use the
`arb_lt`/`arb_gt` loud-on-uncertainty form, as F2's `delta_build` does
(`src/aic_factorize_delta2.c:52`).

### 3.5 `ξ_j` (`.tex:2859`, shard H #5, D4 §3)

`ξ_j ∈ E_j` is a unit vector with `‖C_jξ_j‖ ≈ 1`; the recipe is the **leading
RIGHT singular vector of `C_j`** (then `‖C_jξ_j‖ = σ_max(C_j)`, which lem_RC(ii)
certifies `≥ 1−O(η)`).

**SVD API choice (a degeneracy subtlety, §C4/C5/idemp lesson).** The arb
`aic_mat_singular_values` (`include/aic_mat.h:124`) requires a SIMPLE Gram
spectrum and ABORTS on a repeated singular value — and `C_j` at η=0 is a state
on `E_j` (a density matrix, `‖C_j‖=1`) whose other singular values may be
degenerate. So for the *vector*, use the **double-path full SVD**
`aic_latd_svd(svals, NULL, Vt, C_j_arr, e_j, e_j)` (`include/aic_latd.h:103`),
which HANDLES degeneracy. The leading right singular vector is the **conjugate of
the first row of `Vt`** (`Vt(0,:) = v_0†`, so `ξ_j[c] = conj(Vt[0·e_j + c])`),
then re-embed into an `e_j`-vector and renormalize (it is already unit). Convert
`C_j` to the double array with `aic_latd_from_acb_mat`.

This mirrors the project's settled pattern: degenerate eig/SVD → LAPACK double
path (uncertified vector), arb path → CERTIFY the scalar bound (here
`σ_max ≥ 1−O(η)` via §3.4). Cross-check: `‖C_j ξ_j‖` (computed in arb) must equal
`σ_max(C_j)` (the certified ball) to ~1e-12 — ties the double-path vector back to
the certified bound (catches a wrong `Vt` row / a left-vs-right swap, the §C4
trap).

### 3.6 `L_j` and `Υ'_j` (`.tex:2860-2870`)

```
L_j = Σ_s p_{js} (Δ(U_{js}†) ⊗ 1_F) V W_j† (U_{js} ⊗ ξ_j)     L_j: L_j → H⊗F     (.tex:2862)
Υ'_j(X) = L_j† (Φ(X) ⊗ 1_F) L_j                                Υ'_j: B(H)→B(L_j)  (.tex:2869)
```

Shapes (per term, with `s = (a,b)`, `U_{js} = S_{ab}`):
- `U_{js} ⊗ ξ_j`: `L_j ⊗ E_j ← L_j`, i.e. a `(d_j·e_j) × d_j` matrix mapping
  `L_j → L_j⊗E_j` (a column-per-`L_j`-basis with the `ξ_j` ancilla attached).
  Build via `aic_mat_kronecker(SxiE, S_ab, ξ_j)` where `ξ_j` is the `e_j × 1`
  column (LEFT = `L_j`, RIGHT = `E_j`, matching §3.2).
- `W_j†`: `H ← L_j⊗E_j`, `N × (d_j·e_j)`. So `W_j†(U_{js}⊗ξ_j)`: `N × d_j`.
- `V`: `H⊗F ← H`, `(N·r) × N` (from §3.3). So `V W_j†(U_{js}⊗ξ_j)`: `(N·r) × d_j`.
- `Δ(U_{js}†) ⊗ 1_F`: `H⊗F ← H⊗F`, `(N·r) × (N·r)`. `Δ(U_{js}†)` is `N × N` via
  `aic_factorize_delta(out, F, U_{js}†, prec)`. **`U_{js}†` is the per-block Pauli
  embedded into the `j`-th block of an `n_B`-rep** (zeros elsewhere) — build
  `S_{ab}†` (`d_j×d_j`), place it in block `j` of an `n_B×n_B` matrix (offset
  `r_off = Σ_{l<j} d_l`), then `Δ(·)`. Embed `Δ(U_{js}†) ⊗ 1_F` to match `V`'s
  ordering (§3.3, `F` MAJOR, `H` MINOR) → use `aic_mat_kronecker(., I_r,
  Δ(U_{js}†))` (`1_F ⊗ Δ(...)`, F left). Result: `(N·r) × d_j`.
- Accumulate over the `d_j²` terms with weight `p_{js} = d_j^{-2}` → `L_j`,
  `(N·r) × d_j`.

`Υ'_j(X)`: `Φ(X) ⊗ 1_F` is `(N·r)×(N·r)`; build `Φ(X)` via `aic_ucp_apply(.,
F->phi, X, prec)` (`N×N`), embed `1_F ⊗ Φ(X)` (F left, §3.3). Then
`Υ'_j(X) = L_j† (1_F⊗Φ(X)) L_j`, a `d_j × d_j` operator on `B(L_j)` — manifestly
CP (a `S†(·)S`-congruence of the CP map `X↦Φ(X)⊗1_F`).

**The F-ordering elimination (the RECOMMENDED cleaner path, §3.3 alt).** Both
`L_j` and `Υ'_j` only ever pair `V` with `(Δ(U)⊗1_F)` on the left and `Φ(X)⊗1_F`
in the middle. Using `Φ(X) = V_Φ†(X⊗1_F)V_Φ` we have, with `V = V_Φ`:
`L_j† (Φ(X)⊗1_F) L_j`. Rather than tracking F-orderings, build everything from
`V` (the `(N·r)×N` Stinespring) and the `N×N` operators `Δ(U)`, `Φ(X)`, using
`aic_mat_kronecker(., I_r, ·)` UNIFORMLY for every `⊗1_F` (F left = `V`'s
ancilla-major layout). The ONLY non-`1_F` Kronecker is `R_j`/`U⊗ξ_j` in
`L_j⊗E_j` (L left, §3.2). Two orderings total, each pinned, each tested by a
reconstruction tooth (§4). Keep `r` small (Φ's Kraus rank); the `(N·r)`
dimension is the cost driver — see §E perf note.

`Υ' = (Υ'_1,…,Υ'_m): B(H) → B`. Assemble the `B`-rep output as the
block-diagonal `n_B×n_B` matrix with `Υ'_j(X)` in block `j`.

### 3.7 Unitalize `Υ` (`.tex:2896-2899`)

```
Υ'(1_H) ∈ B   (n_B×n_B block-diagonal, block j = Υ'_j(1_H)),
Υ(X) = Υ'(1_H)^{−1/2} Υ'(X) Υ'(1_H)^{−1/2}.
```
Mirror F2's `delta_build`/`delta` exactly (`src/aic_factorize_delta2.c`):
- ASSERT `‖Υ'(1_H) − 1_B‖_op < 1` (shard H #7, the inverse-sqrt basin for
  `aic_funcalc_xpow(., -0.5, 1.0)`) BEFORE the power — `arb_lt(e, one)`,
  loud-on-uncertainty. `Υ'(1_H) = Υ'(1_H rep on B) = (Υ'_j(1_H))_j`, and
  `1_H = aic_dhom_B_unit`? NO — `Υ'(1_H)` takes the **ambient** `1_N` (X = `1_N`
  on `B(H)`), output in `B`. Build `Υ'(1_N)` (the block-diagonal assembly with
  `X = I_N`), subtract `I_B = aic_dhom_B_unit`, op-norm, assert `< 1`.
- `Υ'(1_H)^{−1/2}` via `aic_funcalc_xpow(M, -0.5, 1.0, prec)` on the `n_B×n_B`
  `Υ'(1_N)`. Store as `F->upsI_invsqrt` (a new OWNED field, n_B×n_B).
- `Υ(X) = upsI_invsqrt · Υ'(X) · upsI_invsqrt` (n_B×n_B congruence) — UCP
  (unital by construction, CP because congruence preserves CP).

---

## 4. THE CENTRALITY TEETH (the orchestrator's specific requirement)

F2's tests verify `Δ'` CP + `O(η)`-close — which hold for ANY unitary set with
`Σp_s = 1`; they do NOT exercise the genuine **1-design centrality** `diag_j2`
(`.tex:2776`): `Σ_s p_{js} X U_{js}† ⊗ U_{js} = Σ_s p_{js} U_{js}† ⊗ U_{js} X`
for all `X ∈ B(L_j)`. But lem_RC(i) (`R_j = 1_{L_j}⊗C_j`) is **exactly** where
centrality becomes load-bearing (the proof, `.tex:2848-2849`: "Due to the
property (diag_j2) of the diagonal, `R_j` commutes with `X⊗1_{E_j}`"). So F3
MUST carry a centrality tooth. Two complementary checks:

### 4.1 DIRECT centrality check (machine-precision; the structural tooth)

For each block `j`, form both sides of `diag_j2` for a non-scalar Hermitian test
`X_0 ∈ M_{d_j}` (e.g. `diag(1, 2, …, d_j)`):
```
LHS = Σ_{a,b} d_j^{-2} (X_0 S_{ab}†) ⊗ S_{ab}
RHS = Σ_{a,b} d_j^{-2} (S_{ab}† ) ⊗ (S_{ab} X_0)
assert ‖LHS − RHS‖_op < 1e-9.
```
(The generalized-Pauli sum `Σ_{ab} d^{-2} S_{ab}†⊗S_{ab}` is the SINGLE-block
Haar diagonal and IS central — FINDINGS §A2 confirms the per-block design is
correct; the §A2 defect was only for the cross-sector JOINT direct-sum, which F3
does NOT use, because lem_RC works block-by-block.)

**MUTATION (Rule 7):** (a) perturb a weight (`p_{j,0} *= 1.3`, renormalize) →
the sum loses the first moment cancellation → `‖LHS−RHS‖ = O(1)` RED; (b) drop
one Pauli (skip `(a,b)=(0,1)`) → likewise RED. Confirm both, then restore.

### 4.2 STRUCTURAL centrality teeth via `R_j` (the load-bearing one)

```
assert ‖R_j − 1_{L_j} ⊗ C_j‖_op = O(η)   (η=0: machine-zero, EXACT 1⊗C_j).
```
Build `1_{L_j} ⊗ C_j` via `aic_mat_kronecker(reconstr, I_{d_j}, C_j)` (L left,
§3.2) and op-norm the difference. This single check certifies BOTH lem_RC(i) AND
the trace direction (§2.2) AND that the design is genuinely central — it FAILS if
the design is not central (the `s=0`-only HUNT1 mutation, §0; a dropped Pauli; a
perturbed weight) because then `R_j` is NOT in the commutant of `B(L_j)⊗1` and
has a nonzero off-`1⊗C_j` part. **MUTATION:** drop a Pauli from `R_j`'s sum →
`‖R_j − 1⊗C_j‖ = O(1)` RED.

Both teeth go in `test_factorize.c` (T6, η>0) and the η=0 oracle (T5).

---

## 5. THE TEST PLAN (Rule 5/6: every test asserts a bound or cross-check)

### T5 — η=0 oracle (the cleanest ground truth; FINDINGS ladder rung 3)

For an EXACTLY idempotent Φ (e.g. `block_cond_exp`, `noiseless_subsystem`),
Φ̃ = Φ, A = Img Φ is a genuine C*, and the construction must reduce to the
`th_idemp_structure` decode (`aic_idemp_decompose`). Assert:
- **`Υ Δ = 1_B` and `Δ Υ = Φ` EXACTLY as MAPS** (gauge-invariant — compare the
  superoperators / apply to matrix units, NOT basis-by-basis; cf. `aic_idemp`
  gauge note). The η=0 limit of `Υ` is `Γ ∘ C_M` (the `th_idemp_structure`
  decode, `aic_idemp_decomp.Gamma`, `.w`/`C_M`); `Δ` is the inclusion. So
  `Υ Δ − 1_B` and `Δ Υ − Φ` are machine-zero (defects → 0). Use the F4 cb-norm /
  apply-to-matrix-units route; for F3-local testing, `max_{ab} ‖Υ(Δ(E_{ab})) −
  E_{ab}‖_op` and `max_{ij} ‖Δ(Υ(E_{ij})) − Φ(E_{ij})‖_op`.
- **lem_RC at η=0:** `R_j = 1_{L_j} ⊗ C_j` with `‖C_j‖ = 1` EXACTLY
  (`σ_max(C_j) = 1`, the certified ball contains 1 with radius ~2^-prec). At
  η=0, `C_j = (1/d_j)Tr_{L_j}(W_jW_j†)` is a state on `E_j` (shard H η=0 oracle;
  D4 §3 dovetail with prop_Gamma's `γ_j`).
- `Υ'(1_H) = 1_B` exactly (so the unitalization is the identity congruence).

The η=0 oracle is BLIND to the star (§C2) and to the left/right SVD distinction
(§C4) — so T5 alone is NOT sufficient; T6 (η>0) is mandatory.

### T6 — η>0 (the genuine teeth; `make_mixconj(4,2)` and `(5,2)`)

- **`Υ` UCP:** unital (`Υ(I_N) = I_B`, by construction) AND CP. CP via the
  **per-block Choi PSD** (`Υ` lands in `B = ⊕_j M_{d_j}`, so `Υ` CP iff each
  block-component `Υ_j: B(H)→M_{d_j}` is CP iff its Convention-A Choi
  `C_j^{Υ} = Σ_{pq} E_{pq}⊗Υ_j(E_{pq}) ∈ B(H)⊗M_{d_j}` is PSD). Verdict via
  `aic_ucp_is_cp_choi` — **EXPECT the §C5 false-fail** and apply the
  midpoint+Weyl verdict F2 already established (§D risk 3 below).
- **`‖Υ − Υ̃‖_cb = O(η)`** (`.tex:2899`). Measure the cb-norm of the difference
  (F4's diamond-SDP machinery; for F3-local, the op-norm of the Choi difference
  as a coarse gate). `Υ̃ = aic_factorize_upsilon_tilde`. Assert
  `‖Υ−Υ̃‖_cb / η ≤ C_max` (a generous bounded gate, the §D4 posture).
- **`σ_max(C_j) ≥ 1 − O(η)`** (lem_RC(ii)) — the certified arb assertion (§3.4).
- **The centrality teeth (§4.1 + §4.2)** — both, with their mutations.
- **The partial-trace-direction PIN (§2.3)** — the asymmetric `2⊗3` fixture, in
  `test_mat.c`-style or inline in `test_factorize.c`.
- **`Υ' Δ ≈ 1_B`** at n=1 (`.tex:2871-2893`): `max_k ‖Υ'_j(Δ(B_k-block)) −
  (B_k)_j‖_op ≤ O(η)` — the explicit `Υ'Δ≈1_B` calc, a direct check of the §3.6
  contraction wired correctly.

### Mutation teeth (Rule 7, the "test has caught a real regression" discipline)

- Trace direction: swap `partial_trace_left`→`right` in `C_j` → T6 §4.2
  (`‖R_j−1⊗C_j‖`) RED + shape assert.
- `ξ_j` left↔right singular vector: take the leading LEFT singular vector
  (first column of `U`) instead of the right → `‖C_jξ_j‖ ≠ σ_max` → the
  cross-check (§3.5) RED. (The §C4 trap: blind on Hermitian `C_j`, so use a
  NON-Hermitian `C_j` fixture — at η>0 `C_j` is generically non-Hermitian.)
- Design centrality: §4.1/§4.2 mutations (drop Pauli / perturb weight).
- W_j ordering: swap `W_j` to ancilla-major (§3.2 route ii) without flipping the
  trace → `‖R_j−1⊗C_j‖` RED.
- Unitalization: skip the `^{−1/2}` congruence → `Υ` not unital → `Υ(I_N) ≠ I_B` RED.

---

## 6. RISKS / OPEN ITEMS (ranked, with mitigation + η=0 constraint)

**Rank 1 — partial-trace / factor-order direction (pin empirically FIRST).**
The recurring FINDINGS bug-class (§C12/D3, §C12.O2-PIN). THREE direction choices:
`R_j`'s `L_j⊗E_j` order (§2.2: L left, trace LEFT), `W_j`'s stacking order
(§3.2: build `L_j`-major directly), and `Φ`'s `⊗1_F` order (§3.3: F left = V's
ancilla-major; or eliminate via the §3.3-alt). *Mitigation:* the §2.3 asymmetric
`2⊗3` fixture + the §4.2 `‖R_j−1⊗C_j‖` tooth; build the §2.3 pin BEFORE any L_j
code. *η=0 constraint:* the oracle gives EXACT `R_j = 1⊗C_j`, `‖C_j‖=1`, so the
tooth is machine-zero — a wrong direction is caught at η=0 already (no need to
wait for η>0).

**Rank 2 — the `W_j` extraction from `Δ` (gauge/ordering, §3.1/§3.2).** Route (b)
(per-block Stinespring) recommended; the Kraus extraction `aic_ucp_choi_to_kraus_latd`
is gauge-free up to a unitary on `E_j` — but `R_j = Σ p (S†⊗1)W_jW_j†(S⊗1)` and
`C_j` depend ONLY on `W_jW_j†` (gauge-invariant: a unitary `u` on `E_j` gives
`W_j → (1⊗u)W_j`, so `W_jW_j† → (1⊗u)W_jW_j†(1⊗u)†`, and `C_j → u C_j u†`, same
`σ_max`, same `‖C_jξ_j‖`). So the gauge does NOT affect any tested quantity.
*Mitigation:* verify `W_j†(X_j⊗1)W_j = Δ(ι_j(X_j))` (reconstruction tooth) and
`Σ_j W_j†W_j = 1_H`. *η=0 constraint:* `W_jW_j†` is the range projection of
`W_j`; the oracle ties it to `th_idemp_structure`'s carrier blocks (D4 §3
prop_Gamma dovetail).

**Rank 3 — the §C5 false-fail on the `Υ'` / `Υ` Choi (it WILL bite).** F2's
lesson (FINDINGS §C5, §C12.O2 SUBTLETY): a deep arb matmul chain (here `L_j`
involves `V`, `W_j†`, `Δ(U)` products) accumulates a radius that trips
`aic_mat_opnorm`/`is_hermitian`/`is_cp_choi`'s Hermiticity predicate even on a
genuinely-Hermitian/PSD Choi. *Mitigation:* the established midpoint+Weyl verdict
— for the CP check, collapse the per-block Choi `C_j^{Υ}` to its MIDPOINT and
symmetrize `(C+C†)/2` before `aic_ucp_is_cp_choi`, with the Weyl `R`-inflation
(the `aic_mat_gram` fix is now in the substrate, FINDINGS §C5 RESOLVED — so plain
`aic_mat_opnorm` no longer false-fails; but `aic_ucp_is_cp_choi` runs
`herm_max_eig(-C)` which asserts Hermiticity at the tight tol, so the
midpoint-symmetrize-before-eig step from §C12.O2 SUBTLETY is still needed for the
arb-assembled Choi). *η=0 constraint:* at η=0 the Choi is exactly PSD with a
clean `{0,σ}` spectrum, low risk; the false-fail is an η>0 / deep-chain
phenomenon.

**Rank 4 — `σ_max(C_j)` certification near `1−O(η)` (§3.4).** The lower bound is
`1−O(η)`; if η is not tiny, the certified ball may straddle the threshold.
*Mitigation:* Rule 4 abort with a clear message ("lem_RC(ii) precondition: σ_max
ball straddles 1−c·η at prec P"); bump prec; for the test, choose mixconj
deviation small enough (0.02–0.03, as F2 used) that `σ_max ≳ 0.97`. *η=0
constraint:* `σ_max = 1` exactly, ball radius ~2^-prec — clears trivially.

**Rank 5 — the composite `O(η)` constant for `‖Υ−Υ̃‖` (defer to F4 / `aic-1bc`).**
NOT a constructivity wall (FINDINGS §D4, D4 §5): the algorithm runs and produces
`Υ`. The certified analytic `C` composes prop_P `c1`, the iso `c0`
(analytically OPEN, `aic-1bc`), and the CP-ization `c3` — none multiplied out in
the paper. *Mitigation/posture:* MEASURE `C = ‖Υ−Υ̃‖_cb / η` per instance, assert
bounded + dimension-independent (the §D2 robust canary: abs-max + no-upward-trend
over an n-sweep, NOT the fragile within-family ratio). File the analytic `C` as a
research bead chained after `aic-1bc`. *η=0 constraint:* `Υ = Υ̃` exactly, `C·0 =
0` — vacuously satisfied.

**No genuine wall.** Per FINDINGS §D4 (BUILDABLE-MODULO) expect NONE. If the
`W_j` extraction or the §C5 verdict turns out to need a certified degenerate SVD
(blocked on `aic-w4o.1`), that is a *precision-path* deferral (use the double
path for the vector + arb for the scalar bound, the settled pattern), NOT a
stop-condition.

---

## 7. FILE / LOC PLAN (Rule 10: ≤200 LOC/file)

F3 EXTENDS the `factorize` module. Proposed split:

- **`src/aic_factorize_upsilon.c`** (lem_RC core, ~180 LOC): the per-block
  Stinespring `W_j` (Route b, §3.1/§3.2), `R_j` (§3.4), `C_j` via
  `partial_trace_left` (§2.2), the `σ_max(C_j) ≥ 1−O(η)` assert (§3.4), `ξ_j` via
  `aic_latd_svd` (§3.5). Helpers: `aic_factorize_Wj` (build `W_j`, `(d_j·e_j)×N`,
  + recovered `e_j`), `aic_factorize_Rj_Cj` (build `R_j`, return `C_j` +
  `sigma_max` ball), `aic_factorize_xi_j` (leading right singular vector).
- **`src/aic_factorize_upsilon2.c`** (`L_j`, `Υ'_j`, unitalize; ~180 LOC): `V`
  from `aic_ucp_kraus_to_stinespring` (§3.3), `L_j` (§3.6), `Υ'_j(X)` and the
  block-diagonal assembly `Υ'(X)` (§3.6), `aic_factorize_upsilon_build`
  (precompute `Υ'(1_N)^{−1/2}` → `F->upsI_invsqrt`, with the `‖Υ'(1_H)−1_B‖<1`
  assert, §3.7), `aic_factorize_upsilon(out, F, W, prec)` (the UCP `Υ`, §3.7).
- **`include/aic_factorize.h`** additions: a new OWNED field `acb_mat_t
  upsI_invsqrt;` + `int upsilon_ready;` on `aic_factorize`; declarations for
  `aic_factorize_upsilon_build`, `aic_factorize_upsilon`, and the per-block
  helpers (`Wj`, `Rj_Cj`, `xi_j`, `Lj`, `upsilon_prime_block`,
  `upsilon_block_choi`) with the F3 banner (paper technique vs constructive
  route, §3.1 Route (b), the three direction pins). Update `aic_factorize_clear`
  to free `upsI_invsqrt` when `upsilon_ready`.
- **`tests/test_factorize.c`** additions: **T5** (η=0 oracle, §5) and **T6**
  (η>0 mixconj 4,2 / 5,2, §5) with all the §4 centrality teeth, the §2.3
  direction pin, and the §5 mutation teeth.

Cite, in every routine: the `.tex` line, the verbatim formula, and
`docs/research/factorize_f3_spec.md §N`; cite `FINDINGS.md §A2/§C5/§D4` where they
bite. Keep `include/aic_mat.h:43-60` (the tensor convention) cited at `R_j`/`C_j`.

---

## Appendix A — API mapping table

| Construction (`.tex`) | Formula | Existing API | Shape |
|---|---|---|---|
| `Δ`'s per-block Choi (§3.1 b) | `C_{Δj}[a·N+p,b·N+q]=Δ(E^{(j)}_{ab})[p,q]` | `aic_factorize_delta` + manual fill (cf. `delta_block_choi`) | `(d_j·N)²` |
| Kraus of `Δ_j` | PSD eig of `C_{Δj}` | `aic_ucp_choi_to_kraus_latd` | `{D_{j,c}: N→d_j}`, `c∈[0,e_j)` |
| `W_j` (§3.2, L-major) | `W_j[a·e_j+c,p]=D_{j,c}[a,p]` | manual stack (NOT `kraus_to_stinespring`) | `(d_j·e_j)×N` |
| per-block Pauli `U_{js}` | `S_{ab}=X^aZ^b` | `aic_dhom_pauli(S, d_j, a, b)` | `d_j×d_j` |
| `S_{ab}⊗1_{E_j}` | (L left) | `aic_mat_kronecker(., S_ab, I_{e_j})` | `(d_j·e_j)²` |
| `R_j` (§3.4) | `Σ_{ab} d_j^{-2}(S†⊗1)W_jW_j†(S⊗1)` | accumulate | `(d_j·e_j)²` |
| `C_j` (§2.2) | `(1/d_j)Tr_{L_j}(R_j)` | `aic_mat_partial_trace_left(., R_j, d_j, e_j)` ×`1/d_j` | `e_j×e_j` |
| `σ_max(C_j)` (§3.4) | `√λ_max(C_j†C_j)` | `aic_mat_herm_max_eig` on `C_j†C_j` | scalar ball |
| `ξ_j` (§3.5) | leading right s.v. of `C_j` | `aic_latd_svd(.,NULL,Vt,C_j,e_j,e_j)`, `ξ=conj(Vt[0,:])` | `e_j×1` |
| `V` (Φ Stinespring, §3.3) | `V[a·N+i,j]=K_a[i,j]` | `aic_ucp_kraus_to_stinespring(V, F->phi)` | `(N·r)×N` |
| `Δ(U_{js}†)` | UCP `Δ` of embedded `S_{ab}†` | `aic_factorize_delta(., F, Uemb, prec)` | `N×N` |
| `⊗1_F` (F left, §3.3) | `1_F ⊗ (·)` | `aic_mat_kronecker(., I_r, ·)` | `(N·r)²` |
| `L_j` (§3.6) | `Σ_s p_{js}(Δ(U†)⊗1_F)V W_j†(U⊗ξ_j)` | accumulate | `(N·r)×d_j` |
| `Φ(X)` | `Σ_a K_a†XK_a` | `aic_ucp_apply(., F->phi, X)` | `N×N` |
| `Υ'_j(X)` (§3.6) | `L_j†(Φ(X)⊗1_F)L_j` | matmul | `d_j×d_j` |
| `Υ'(1_H)^{−1/2}` (§3.7) | inverse-sqrt | `aic_funcalc_xpow(., -0.5, 1.0)` | `n_B×n_B` |
| `Υ(X)` (§3.7) | `M^{−1/2}Υ'(X)M^{−1/2}` | matmul | `n_B×n_B` |

## Appendix B — the per-block design vs the joint design (FINDINGS §A2)

lem_RC works **block by block** (`R_j` involves only `U_{js}` on `L_j`). So F3
uses the **single-block** generalized-Pauli design `{S_{ab} = X^aZ^b : a,b ∈
[0,d_j)}`, `p = d_j^{-2}`, via `aic_dhom_pauli(S, d_j, a, b)` — `d_j²` terms. This
single-block design IS central (`diag_j2`) and IS the Haar second moment (FINDINGS
§A2: "the per-block Pauli sum *does* correctly reproduce the single-block Haar
second moment — the single-block case is fine; only m≥2 is broken"). The §A2
non-centrality defect afflicts ONLY the cross-sector JOINT direct-sum
`U_{1s_1}⊕…⊕U_{ms_m}` (F2's `aic_factorize_design_unitary`), which F3 does NOT
use. Do NOT call `aic_factorize_design_unitary` or `aic_dhom_diag_build` in F3.
