# opspace O2 — certified Watrous-SDP cb-norm UPPER bound for `v: B → A`

**Bead:** aic-pjr (opspace O2; depends on aic-zwo = O1, reuses aic-m24 / aic-0at infra)
**Date:** 2026-05-31
**Ground truth:** `paper/src/approximate_algebras.tex` (cited by line throughout)
**Status:** DESIGN (read-only; no code/build/Julia run in this pass)
**Reads:** `docs/research/opspace_web_leg.md`, `opspace_paper_leg.md`; `paper/FINDINGS.md §C12 §D3`;
`include/aic_cbnorm.h`, `src/aic_cbnorm_certify.c`, `src/aic_cbnorm_certify_restore.c`,
`src/aic_cbnorm_internal.h`; `julia/AlmostIdempotentChannels.jl/src/sdp.jl`;
`tools/gen_fixtures_{d24,certify}.jl`; `include/aic_dhom.h`, `include/aic_ucp.h`,
`src/aic_ucp_core.c`/`aic_ucp.h`; `include/aic_opspace.h`, `src/aic_opspace_{map,cert}.c`.

---

## 0. Executive summary

O1 (aic-zwo) computes a HOPM **lower** bound on the ampliated operator-norm
max-stretch `‖v_n‖_op` (and, via `1/‖v_n^{-1}‖_op`, an **upper** bound on the
lower-isometry `a_n`). O2 must produce the matching **certified upper** bound on
`‖v‖_cb` and `‖v^{-1}‖_cb` so that `th_main_ext`'s claim
`‖v‖_cb ≤ 1+O(ε)`, `‖v^{-1}‖_cb ≤ 1+O(ε)` (`.tex:1538-1540`) becomes a rigorous
arb statement. The route (web leg §2; FINDINGS §D3) is the standard one:

> `‖v‖_cb = ‖v*‖_⋄` (cb-spectral norm = diamond norm of the adjoint; QETLAB `CBNorm`
> = `DiamondNorm` on the adjoint), and the diamond-norm SDP is exactly the one the
> project already solves in `src/sdp.jl` (Watrous 2012, arXiv:1207.5726).

So O2 = **feed `J(v*)` into the existing diamond-norm SDP, restore the feasible
points to exact feasibility in arb (the `aic_cbnorm_certify` machinery), read off
the rigorous ball.** The deltas vs the existing `‖Φ²−Φ‖_cb` pipeline:

1. **Input Choi** is `J(v*) = J(v)^T` (transpose), not `Choi(Φ²−Φ)`. `J(v)` is
   `N·n_B × N·n_B`, generically **rectangular block structure** (input space M_N,
   output space M_{n_B}) — NOT a square self-map. This is the single biggest
   structural change: `src/sdp.jl` and the C certifier are hard-wired to a square
   `n²×n²` self-map and must be generalized to two dims `(d_in, d_out)`.
2. **Normalization.** Derived in §2.4: the diamond norm of a map presented by a
   **Convention-A Choi with trace = d_in** is `‖·‖_⋄ = (2/d_in)·SDP_optval` for the
   MAX-primal, `optval/2` for the MIN-dual — *the same `(2/·)` family*, but with
   `d_in` = the **input** ambient dim of the map being normed (for `v*` that is `N`;
   for `(v^{-1})*` that is `n_B`). The existing self-map uses `d_in = n`. The eta=0
   oracle (`‖v‖_cb = 1` exactly) is the anchor that pins it.
3. **cb-SPECTRAL vs cb-TRACE.** `src/sdp.jl` computes the **diamond / cb-trace**
   norm `‖·‖_⋄ = |||·|||_1`. We want the **cb-spectral** norm `‖v‖_cb = |||v|||_∞`.
   The bridge is the adjoint duality `|||v|||_∞ = |||v*|||_1` (Watrous 2009 eq.(1)):
   compute the diamond norm **of the adjoint map `v*`**, i.e. feed `J(v*)`, not `J(v)`.
4. **Partial-trace direction is the load-bearing convention** (the bug that bit the
   original cbnorm work, HANDOFF:340). The dual objective `‖Tr_2(Y0)‖_∞+‖Tr_2(Y1)‖_∞`
   traces the **MINOR/input factor**. For `v*` the input factor is the **A-ambient
   (N)** factor; getting it wrong silently doubles `hi`. Pin it with the eta=0
   oracle AND an asymmetric eta>0 channel (the paper-example analogue).

**Cross-check (the load-bearing one, aic-0at):** `HOPM(O1) ≤ SDP(O2)` brackets the
true `‖v_n‖_op`; combined with Smith they bracket `‖v‖_cb`. Plus the eta=0 oracle
`‖v‖_cb = ‖v^{-1}‖_cb = 1` exactly (lo=hi=1 to arb radius).

**Smith is MOOT for O2** (§6.3). Unlike O1 — which builds an *explicit* ampliated
matrix and therefore needs Smith to know *which* level N to build — the Watrous SDP
on `J(v)` intrinsically captures `sup_k ‖1_{M_k}⊗v‖` with NO explicit ampliation
(the Choi matrix already encodes all k). The SDP returns `‖v‖_cb` directly; Smith's
`N_max` truncation is not used in O2's computation. Smith remains the *theorem* that
licenses O1's single-level HOPM and that says the two routes target the same number.

**Biggest risk:** `v` is **Hermiticity-preserving but NOT completely positive**, so
`J(v)` is Hermitian but **indefinite** (not PSD). The Watrous diamond-norm SDP is
valid for *any* linear (HP) map — but the existing pipeline was exercised only on
the HP-but-indefinite `Λ=Φ²−Φ`, which is the same regime, so this is *safe*; the
genuine new hazard is the **rectangular (non-self-map) generalization** of the SDP
and the partial-trace direction under it.

**Doc path:** `docs/research/opspace_o2_design.md`.

---

## 1. The Choi matrices `J(v)`, `J(v^{-1})`, `J(v*)`

### 1.1 The maps and their ambient dims

From `include/aic_dhom.h:103-122`:
- `B = ⊕_l M_{d_l}(C)`, block-diagonal in `M_{n_B}`, `n_B = Σ_l d_l`, `dim_B = Σ_l d_l²`.
  Matrix-unit basis `{E_i}`, `i = 0..dim_B-1`, each `E_i` an `n_B × n_B` block-diagonal
  matrix with a single `1` (linear index `μ(l,a,b)=offset_l + a·d_l + b`,
  `aic_dhom_B_index`).
- `A ⊆ M_N`, `N = v->n` (ambient). `dim_A = dim_B` in the bijective case
  (`aic_opspace_cert.c:45`).
- `v: B → A` is stored by `vE[i] = v(E_i)` (each `N × N`), and `v(X) = Σ_i x_i vE[i]`
  is **linear in the B matrix-unit coordinates** `x_i` (`aic_dhom_v_apply`,
  `aic_dhom.h:30-31, 180-184`).

So `v: M_{n_B} ⊇ B → A ⊆ M_N`: **input ambient `n_B`, output ambient `N`**.

### 1.2 `J(v)` in the project's Choi Convention A

The project's Convention A (`include/aic_ucp.h:24-42`, verbatim from the header for
a UCP map `Φ: B(K) → B(H)`):

```
C_Phi[i*dim_H + a, j*dim_H + b] = sum_x conj(K_x[i,a]) * K_x[j,b]
                                = Phi(E_{ij})[a,b],
```

i.e. `C_Φ = Σ_{ij} E_{ij} ⊗ Φ(E_{ij})` — **input/source index `(i,j)` in the LEFT/MAJOR
position, output/target index `(a,b)` in the RIGHT/MINOR position** (`aic_ucp.h:28-29`;
generator `choi_A` in `gen_fixtures_d24.jl:37-43`: `C[(i-1)n+a,(j-1)n+b] += conj(K[i,a])K[j,b]`,
1-based). For a *general linear* map this convention is just

```
   J(Φ)[ source_row · (out dim) + out_row ,  source_col · (out dim) + out_col ]
       =  Φ(E_{source_row, source_col})[out_row, out_col].
```

Specialize to `v` (source = B-ambient `n_B`, out = A-ambient `N`):

> **`J(v)[ p·N + a ,  q·N + b ]  =  v(E_{pq})[a,b]`**, with `p,q ∈ [0,n_B)` the
> B-ambient (source/MAJOR) indices and `a,b ∈ [0,N)` the A-ambient (output/MINOR)
> indices. Equivalently `J(v) = Σ_{p,q} E_{pq}^{(n_B)} ⊗ v(E_{pq})`.
>
> **Dimension:** `(n_B · N) × (n_B · N)`.

Here `E_{pq}` ranges over **all** `n_B²` matrix units of the *ambient* `M_{n_B}`, but
`v` is only defined on the `dim_B`-dimensional subspace `B`. The off-block units (those
`E_{pq}` with `(p,q)` not inside any block `l`) are **NOT in the domain of `v`**. Two
correct readings, both valid for the cb-norm-into-`M_N` we want:

- **(Recommended) Define `v` to be ZERO on off-block units.** `v(E_{pq}) = 0` whenever
  `(p,q)` is off-block. This makes `J(v)` a genuine `M_{n_B} → M_N` Choi and is exactly
  what Smith's `u: E → M_N` framing measures (`E = B ⊆ M_{n_B}` as an operator space, the
  off-block directions contribute 0 stretch). The cb-norm of this extended-by-zero map
  equals `‖v‖_cb` as a map on the operator space `B` (the zero block adds no stretch).
- The alternative ("restrict the SDP to the `B` subspace") needs a projected SDP and is
  not worth the complexity; the extend-by-zero `J(v)` is the clean object.

**C assembly from `aic_dhom_v`** (the ~50-LOC helper, §5). `vE[i] = v(E_i)` with `E_i`
the block-diagonal matrix unit at linear index `i = μ(l,a,b)`. The `(p,q)` ambient
position of `E_i` is `b_unit_pos` in `aic_opspace_map.c:30-44` (reuse it verbatim):
`(row,col) = (boff + loc/d_l, boff + loc%d_l)`. So:

```
J(v) = 0   (n_B*N square, zeroed)
for i in 0..dim_B-1:
    (p,q) = b_unit_pos(B, i)        # ambient source position of E_i
    for a in 0..N-1, b in 0..N-1:
        J(v)[p*N + a, q*N + b] = vE[i][a,b]
```

`J(v)` is Hermitian iff `v` is Hermiticity-preserving: `v` HP ⟺ `J(v) = J(v)^†`
(web leg §2.3; the involution-symmetry `v(E_{ba})^† = v(E_{ab})` is already asserted by
`aic_dhom_approx`, `aic_dhom.h:73`). HP holds for the approximate `*`-isomorphism, so
`J(v)` is Hermitian — but **indefinite** (`v` is not CP), exactly as `Λ=Φ²−Φ` is.

### 1.3 `J(v^{-1})`

`v^{-1}: A → B` (input ambient `N`, output ambient `n_B`). O1 already builds `v^{-1}`
as a coordinate map: invert the `dim_A × dim_B` coordinate matrix
`M_1[k,i] = ⟨B_k, vE[i]⟩_F` (`aic_opspace_map.c:46-69`, `build_M1`) via
`acb_mat_solve` (`aic_opspace_map.c:117-130`, fail-loud on singular). With
`M_1^{-1}` in hand, the images `v^{-1}(B_k)` (where `{B_k}` is A's Frobenius-orthonormal
basis, `aic_ecstar.h:91-92`) are recovered as B-elements:

```
v^{-1}(B_k) = Σ_i (M_1^{-1})[i,k] · E_i        (an element of B, n_B × n_B)
```

But the Choi matrix wants `v^{-1}` on the **A-ambient matrix units** `E_{pq}^{(N)}`,
not on the orthonormal basis `{B_k}`. Since `v^{-1}` is linear and `{B_k}` spans `A`,
write each ambient unit `E_{pq}^{(N)}` in the `{B_k}` basis *projected into A* (the
component orthogonal to A is annihilated — same extend-by-zero reading as §1.2):

```
v^{-1}(E_{pq}^{(N)}) = Σ_k ⟨B_k, E_{pq}^{(N)}⟩_F · v^{-1}(B_k)
                     = Σ_k conj(B_k[p,q]) · v^{-1}(B_k)
```

then

```
J(v^{-1})[ p·n_B + r ,  q·n_B + s ]  =  v^{-1}(E_{pq}^{(N)})[r,s],
```

`p,q ∈ [0,N)` (A-ambient MAJOR), `r,s ∈ [0,n_B)` (B-ambient MINOR), dimension
`(N · n_B) × (N · n_B)`. In practice it is cleaner to assemble `J(v^{-1})` the same
way as `J(v)` but with the roles swapped, reusing the already-built `M_1^{-1}` and the
`{B_k}` matrices — see §5 for the C helper signature (one assembler, parameterized by
"which map").

### 1.4 `J(v*) = J(v)^T`

The adjoint `v*: A → B` (note: same input/output ambient dims as `v^{-1}` — both
`M_N → M_{n_B}`, but `v* ≠ v^{-1}` in general). Its Choi is the **transpose** of
`J(v)` in the input↔output basis (web leg §2.1, lines 209-216):

> `J(v*) = J(v)^T` (transpose, **not** conjugate-transpose).

Derivation: `v*(Y) = Σ_{pq} Tr(v(E_{pq}) Y) E_{pq}` (the trace-pairing adjoint), so
`v*(E_{ab}^{(N)})[p,q] = Tr(v(E_{pq}) E_{ab}^{(N)}) = v(E_{pq})[b,a]`. Reading off the
Convention-A Choi of `v*` (now source = A-ambient `N`, out = B-ambient `n_B`):
`J(v*)[a·n_B + p, b·n_B + q] = v*(E_{ab})[p,q] = v(E_{pq})[b,a]`. Compare
`J(v)[p·N+a, q·N+b] = v(E_{pq})[a,b]`. The index map
`(a·n_B+p, b·n_B+q) ↦ (p·N+a, q·N+b)` with `[b,a]↔[a,b]` is exactly the global matrix
transpose of `J(v)`. Hence **`J(v*) = J(v)^T`** as `(N·n_B)×(N·n_B)` matrices. ✓

Since `v` is HP, `J(v)=J(v)^†`, so `J(v*) = J(v)^T = conj(J(v))` (web leg line 235);
either route gives the same matrix. In code, build `J(v)` then take a plain transpose
(`acb_mat_transpose`) — do **not** confuse with conjugate-transpose.

**For `‖v‖_cb`:** feed `J(v*) = J(v)^T` into the diamond-norm SDP with input dim `= N`
(the ambient of `v*`'s domain `A`). **For `‖v^{-1}‖_cb`:** feed
`J((v^{-1})*) = J(v^{-1})^T` with input dim `= n_B` (the ambient of `(v^{-1})*`'s
domain `B`).

---

## 2. The SDP (Watrous diamond-norm) and the precise deltas vs `src/sdp.jl`

### 2.1 The existing SDP (square self-map)

`src/sdp.jl` solves three forms of the Watrous 2012 program (arXiv:1207.5726) for a
**self-map** `Λ: M_n → M_n` with `J = Choi_A(Λ)`, `n²×n²` (`sdp.jl:55-111`):

```
MAX primal (X,P,Q):  max Re tr(J^† X)
                     s.t. [[P, X],[X^†, Q]] ⪰ 0,  P+Q = I_{n²}
                     ->  eta = (2/n) * optval

MIN dual  (Y0,Y1):   min  ‖Tr_2(Y0)‖_∞ + ‖Tr_2(Y1)‖_∞
                     s.t. [[Y0, -J],[-J^†, Y1]] ⪰ 0,  Y0,Y1 ⪰ 0
                     ->  eta = optval / 2
                     Tr_2 = partialtrace(., 2, [n,n])  (the MINOR/INPUT factor)
```

The arb certifier `aic_cbnorm_certify` (`aic_cbnorm_certify.c`,
`aic_cbnorm_certify_restore.c`) restores each feasible point to exact feasibility and
reads off the rigorous bracket. **All matrices are assumed `n²×n²`** (asserted
`aic_cbnorm_certify.c:44-51`) and the input/output ampliation factors are **both `n`**
(`partial_trace_right(T0, Y0, n, n, prec)`, `aic_cbnorm_certify_restore.c:146`).

### 2.2 Delta 1 — input is `J(v*)`, and it is a NON-SELF-MAP

For `‖v‖_cb` the SDP input is `J(v*)`, the Choi of the map `v*: M_N → M_{n_B}`:
- **input ambient** `d_in = N` (domain of `v*` = A-ambient),
- **output ambient** `d_out = n_B` (codomain of `v*` = B-ambient),
- Choi size `(d_out · d_in) × (d_out · d_in) = (n_B·N) × (n_B·N)`.

The Watrous SDP for a general `Φ: M_{d_in} → M_{d_out}` (arXiv:1207.5726 Thm 6;
web leg §2.1 lines 177-216) is:

```
MAX primal:  max Re⟨J, X⟩
             s.t. [[1_{d_out}⊗ρ_0,  X],[X^†, 1_{d_out}⊗ρ_1]] ⪰ 0
                  ρ_0, ρ_1 density ops on C^{d_in}     (Tr ρ = 1)
                  X ∈ L(M_{d_out} ⊗ M_{d_in})         # the (d_out·d_in)-square block var

MIN dual:    min  (1/2)(‖Tr_{d_in}(Y0)‖_∞ + ‖Tr_{d_in}(Y1)‖_∞)
             s.t. [[Y0, -J],[-J^†, Y1]] ⪰ 0,  Y0,Y1 ⪰ 0 on M_{d_out}⊗M_{d_in}
                  Tr_{d_in} = partial trace over the INPUT (d_in) factor.
```

So the generalization is **two dim parameters** `(d_in, d_out)` where the self-map had
`d_in = d_out = n`. The MAX-primal constraint `P+Q=I_{n²}` becomes the
`1_{d_out}⊗ρ_0 / 1_{d_out}⊗ρ_1` form (Watrous Thm 6 primal); the project's self-map
form `[[P,X],[X^†,Q]], P+Q=I` is the **special case** `d_out=d_in=n` with
`P = 1_n⊗ρ_0`, `Q = 1_n⊗ρ_1` (the `P+Q=I` constraint then encodes
`ρ_0+ρ_1` summing appropriately — see Watrous 2012 §3.2 for the equivalence). For O2 we
implement the general two-dim primal/dual directly.

> **Implementation note.** The existing `src/sdp.jl` functions take `(J, n)`. Add
> generalized variants `(J, d_in, d_out)` (or refactor the existing to call the general
> one with `d_in=d_out=n`). The block variables/marginals are over `M_{d_out}⊗M_{d_in}`;
> the partial trace is over the **`d_in`** factor.

### 2.3 Delta 2 — cb-SPECTRAL vs cb-TRACE (the adjoint duality, why `J(v*)` not `J(v)`)

`src/sdp.jl` computes the **diamond / cb-TRACE** norm `‖Φ‖_⋄ = |||Φ|||_1`
(`sdp.jl:1-2`). The quantity O2 wants is the **cb-SPECTRAL** norm
`‖v‖_cb = |||v|||_∞`. The bridge (web leg §2.1 lines 159-168; Watrous 2009 eq.(1)):

```
|||v|||_∞ = |||v*|||_1       (cb-spectral norm of v = diamond norm of the adjoint v*)
```

Hence: **compute the diamond norm of `v*` by feeding `J(v*)=J(v)^T`** into the (general)
SDP. The SDP machinery itself is unchanged — it always computes `|||·|||_1` of whatever
map it is handed; we hand it `v*` to obtain `|||v*|||_1 = ‖v‖_cb`. (This is precisely
what QETLAB's `CBNorm` does: call `DiamondNorm` on the adjoint — web leg ref 7,
HANDOFF mention.) Symmetrically `‖v^{-1}‖_cb = |||(v^{-1})*|||_1`, feed
`J((v^{-1})*) = J(v^{-1})^T`.

### 2.4 Delta 3 — the normalization factor (DERIVED)

The project's pinned normalization (`sdp.jl:12-16`, `gen_fixtures_d24.jl:10-16`) for a
**self-map** with Convention-A Choi of trace `n` is `‖Λ‖_⋄ = (2/n)·MAX_optval`
(equivalently `MIN_optval/2`). Derive the analogue for `v*`.

The factor arises because a **Convention-A Choi of a map `Φ: M_{d_in}→M_{d_out}` has
trace `= d_in`** when `Φ` is trace-preserving-like (for UCP `Φ`, `Tr J = Tr Φ(I_{d_in})·?`
— concretely for the self-map UCP `Tr C_Φ = n`). The Watrous SDP as written
(`maximize real(tr(J'X))` with `ρ` density / `P+Q=I`) is calibrated so that for a
**trace-1-normalized** Choi the optval equals `|||Φ|||_1/2`. The Convention-A Choi is
scaled by the trace, hence the `(2/trace)` rescale. The trace of a Convention-A Choi is
governed by the **input** dimension (the `Σ_p E_{pp}⊗Φ(E_{pp})` diagonal sum runs over
the `d_in` source units). Therefore:

> **`‖v*‖_⋄ = (2/N)·MAX_optval`** (input dim of `v*` is `N`), so `‖v‖_cb = (2/N)·MAX_optval`.
> **`‖(v^{-1})*‖_⋄ = (2/n_B)·MAX_optval`**, so `‖v^{-1}‖_cb = (2/n_B)·MAX_optval`.
> The MIN-dual gives `optval/2` (no dim factor — the `(1/2)` is absorbed into the
> objective as in `sdp.jl:110`), but its raw objective `‖Tr_{d_in}(Y0)‖+‖Tr_{d_in}(Y1)‖`
> must itself be consistent with the same trace scaling; see the empirical pin below.

**This must NOT be trusted as derived — it must be PINNED empirically on the eta=0
oracle**, exactly as the self-map factor was (`sdp.jl:26-30` warns the design-doc table
was wrong on the direction; only the eta=0 + asymmetric anchors are trustworthy). The
**clean anchor**: at η=0, `v` is an exact `*`-isomorphism (a complete isometry), so
`‖v‖_cb = ‖v^{-1}‖_cb = 1` **exactly**. Solve the SDP on `J(v*)` for an exact-idempotent
fixture; require `(2/N)·MAX_optval = 1` and `MIN_optval/2 = 1` to ~1e-9. If a different
constant (e.g. `2/n_B`, or `2/(N·something)`) is needed to hit 1.0, **that** is the
normalization — record the measured value, do not assume `2/N`. (Candidate subtlety:
whether the trace is `N` vs `n_B` vs `dim_A` for the extend-by-zero `J(v*)`; the eta=0
solve discriminates them since they differ unless `N=n_B`.)

### 2.5 Delta 4 — the partial-trace DIRECTION (the load-bearing bug)

This is the convention that bit the original cbnorm work (HANDOFF:340-341; the long
analysis in `sdp.jl:18-30` and `gen_fixtures_certify.jl:118-127, 180-188`). For the
self-map the MIN-dual traces **sys 2 = the MINOR/INPUT factor** (`a` in
`J=C[(i-1)n+a,(j-1)n+b]`); tracing sys 1 (the output factor) gives the WRONG answer
(ratio 4.0 vs 2.0 on the asymmetric paper example).

For `v*: M_N → M_{n_B}` with `J(v*)` size `(n_B·N)`:
- Convention A puts the **source/MAJOR index** of `v*` (the A-ambient `N`-index) LEFT
  and the **output/MINOR index** (the B-ambient `n_B`-index) RIGHT. Re-examine §1.4:
  `J(v*)[a·n_B+p, b·n_B+q]` — the source index `a ∈ [0,N)` is MAJOR, the output index
  `p ∈ [0,n_B)` is MINOR. So in `J(v*)` the **LEFT/MAJOR factor is the INPUT (N)** and
  the **RIGHT/MINOR factor is the OUTPUT (n_B)**.
- The Watrous dual traces the **INPUT** factor. For the self-map the input was the MINOR
  (right) factor — hence `partial_trace_right`. **But for `J(v*)` the input is the
  MAJOR (LEFT) factor.** So naively reusing `partial_trace_right` would trace the
  OUTPUT — the WRONG direction.

> **CRITICAL.** The partial-trace direction depends on which factor the Choi convention
> places the INPUT in. For the existing self-map `Λ`, input = minor = right ⇒
> `partial_trace_right`. For `J(v*)`, depending on how §1.4's transpose lands the
> indices, the input may be the LEFT/major factor ⇒ `partial_trace_left`, with marginal
> dims `(a=N, b=n_B)` rather than `(a=b=n)`.

**Pin it, do not derive it.** Two anchors:
1. **eta=0 oracle** (`‖v‖_cb=1`): necessary but possibly direction-blind if `J(v*)` is
   symmetric there (the self-map symmetric channels were direction-blind, `sdp.jl:27`).
2. **Asymmetric eta>0 fixture** (the analogue of `paper_example`): require the dual `hi`
   to MATCH the primal `lo` (and the HOPM lower bound) to ~1e-6; the WRONG direction
   doubles `hi` (`gen_fixtures_certify.jl:90,186`). This is the discriminating
   poison-pill — the O2 fixture corpus MUST include an asymmetric eta>0 channel whose
   `v` produces an asymmetric `J(v*)` marginal, or the direction test has no teeth.

**The clean way to remove all ambiguity:** in the generalized SDP, write the partial
trace explicitly as "trace the INPUT factor of dims `(d_in)`", parameterized so the
factor positions are explicit from `(d_in, d_out)` and the Convention-A index layout —
then the eta=0 + asymmetric anchors confirm it. Mirror `gen_fixtures_certify.jl`'s
`ptrace2` (minor) / `ptrace1` (major) helpers and the `asym = ‖Tr_R−Tr_L‖_F` poison
guard (`gen_fixtures_certify.jl:184-188`).

---

## 3. The arb certification reuse (`aic_cbnorm_certify`)

The arb feasibility-restoration certifier (`aic_cbnorm_certify.c` +
`aic_cbnorm_certify_restore.c`) is **dimension-agnostic in its math** and reusable with
**one parameterization change**: the hard-coded `N = n*n` and the
`partial_trace_right(., n, n, .)` must become the general `(d_in, d_out)` block
structure.

### 3.1 What is reused UNCHANGED

- **LOWER recipe** `aic_cbnorm_int_lower` (`restore.c:73-125`): symmetrize the equality
  (`E=P+Q−I; P'=P−E/2, Q'=Q−E/2`), `delta = max(0,−λ_min(block'))` via one-sided
  `herm_max_eig(−block')` (`psd_defect`, `restore.c:50-71`), `t=1/(1+2δ)`,
  `lo = t·(norm)·Re tr(J^† X)`. The **block assembly `block2`** (`restore.c:36-46`) and
  the **convex-combination toward the Slater center** are dimension-agnostic — the only
  change is the `(2/n)` factor → `(2/d_in)` and the block size `2(d_out·d_in)`.
- **UPPER recipe** `aic_cbnorm_int_upper` (`restore.c:127-165`): shift
  `eps=max(0,−λ_min(block_D))`, `(Y0+εI, Y1+εI)` feasible,
  `hi=(1/2)(λ_max(Tr_in Y0)+λ_max(Tr_in Y1))+ε·d_in`. **The `Tr_2` call and the
  `eps*n` shift term must become `Tr_{input}` with the correct direction (§2.5) and
  `eps*d_in`** (the shift contributes `ε·d_in` to the marginal because
  `Tr_{d_in}(εI_{d_out·d_in}) = ε·d_in·I_{d_out}`).
- **Dispatch + fail-loud** (`certify.c:39-114`): the eta=0 regime
  (`‖J‖_F < 1e-9·scale`) → eig-free `~[0,0]`; restoration-off / bracket-inverted →
  eig-free fallback; even that straddles → `abort()`. All reusable; the eta=0 threshold
  `1e-9*n` becomes `1e-9*sqrt(d_in·d_out)` (or just keep a scale proportional to the
  Choi size). **The eig-free fallback** `aic_cbnorm_eigfree_ball_choi` uses
  `[‖J‖_F/n, 2‖J‖_F]` (`aic_cbnorm.h:36-44`) — generalize the `/n` to `/sqrt(d_in·d_out)`
  or the appropriate scale; this fallback is loose by design so the exact constant is
  not load-bearing, but it must remain a *valid* outer bracket.

### 3.2 What is NEW

- A **generalized certifier entry** `aic_cbnorm_certify_rect(lo, hi, J, X, P, Q, Y0, Y1,
  d_in, d_out, norm_factor, prec)` (or refactor `aic_cbnorm_certify` to take
  `(d_in, d_out, norm_factor)` and keep the self-map wrapper passing `(n, n, 2.0/n)`).
- The **partial-trace direction** must be wired to trace the INPUT factor (§2.5); pass
  a flag or pick `partial_trace_left`/`partial_trace_right` per the pinned convention.

### 3.3 The new C glue (~50 LOC) + Julia wiring

- **C:** `aic_opspace_choi_v.c` — `void aic_opspace_choi_v(acb_mat_t Jv, const
  aic_dhom_v *v, slong prec)` assembling `J(v)` (§1.2) from `vE[i]` + `b_unit_pos`, and
  `void aic_opspace_choi_vinv(acb_mat_t Jvi, const aic_dhom_v *v, slong prec)` assembling
  `J(v^{-1})` (§1.3) reusing `build_M1` + `acb_mat_solve` (factor it out of
  `aic_opspace_map.c` or duplicate the ~20 lines). Then `J(v*) = transpose(J(v))` and
  `J((v^{-1})*) = transpose(J(v^{-1}))` (plain `acb_mat_transpose`). Each ~25 LOC ⇒ one
  `<200` LOC file.
- **Julia:** extend `src/sdp.jl` with `diamond_norm_watrous_primal_rect(J, d_in, d_out)`
  and `diamond_norm_dual_rect(J, d_in, d_out)` (exposing the feasible points like the
  existing `*_primal` / `_dual`). A committed-fixture generator
  `tools/gen_fixtures_opspace_o2.jl` builds `v` from a channel (via the existing
  cstar_build pipeline OR a directly-constructed `v`), forms `J(v*)`, solves both SDPs,
  asserts the eta=0 / asymmetric anchors, and emits `tests/fixtures_opspace_o2.inc.h`
  with the feasible points (so `make test` stays Julia-free — the
  `gen_fixtures_certify.jl` pattern, lines 154-276).

---

## 4. Cross-checks

### 4.1 The HOPM(O1) ≤ SDP(O2) bracket (load-bearing, aic-0at)

This is the headline cross-check (web leg §4; `aic_opspace.h:55-72`; bead aic-pjr
description). The two routes are honest about their directions:

- **O1 HOPM** gives a **LOWER** bound on `‖v_N‖_op` (and `‖v_{n_B}^{-1}‖_op`):
  `cb_forward = HOPM(‖v_N‖_op) ≤ ‖v_N‖_op` (`aic_opspace.h:55-58, 103-104`).
- **O2 SDP** gives a certified **UPPER** bound on `‖v‖_cb = sup_n ‖v_n‖_op` (Smith ⇒
  `= ‖v_N‖_op` for the forward map).

Therefore, by Smith (forward `N_max = N`, `‖v‖_cb = ‖v_N‖_op`):

```
HOPM_lower(‖v_N‖_op)   ≤   ‖v_N‖_op = ‖v‖_cb   ≤   SDP_upper(‖v‖_cb).
        (O1)                  (truth, Smith)            (O2)
```

The **assertion**: `O1.cb_forward ≤ O2.hi + slack` (with the arb radii + HOPM-gate
slack). For the inverse map: `1/O1.a_cb` lower-bounds `‖v_{n_B}^{-1}‖_op`, and
O2's `‖v^{-1}‖_cb` upper-bounds it; assert `1/O1.a_cb ≤ O2_inv.hi + slack`. When the
HOPM has converged to the true max-stretch and the SDP is tight, the bracket is narrow
— that mutual confirmation is the strongest test (the eta>0 regime).

### 4.2 The eta=0 oracle (`‖v‖_cb = ‖v^{-1}‖_cb = 1` exact)

At η=0, `v` is an **exact `*`-isomorphism**, hence a **complete isometry**:
`‖v_n‖_op = ‖v_n^{-1}‖_op = 1` for all n (`aic_opspace.h:60-65`; FINDINGS §C12 (a),
measured to 1e-12 in O1). So:

- `J(v*)` is a complete-isometry Choi; the SDP MUST return `‖v‖_cb = 1` exactly
  (`lo` and `hi` both enclose 1 to arb radius).
- This is the **normalization anchor** (§2.4) AND the cleanest ground truth
  (CLAUDE.md cross-check ladder rung 3). Use the `block_cond_exp` / `noiseless_subsystem`
  exact-idempotent fixtures O1 already uses (FINDINGS §C12 line 502-503).
- The eta=0 dispatch in `aic_cbnorm_certify` (`certify.c:53-68`) currently routes
  `‖J‖_F<thr → ~[0,0]`. But here `J(v*)` is NOT small at η=0 (`v` is a nonzero isometry);
  `‖v‖_cb=1`, not 0. **So the eta=0 regime branch is the WRONG analogy** — for O2 the
  "trivial" value is `1`, not `0`, and `J(v*)` has `O(1)` Frobenius norm. The certifier
  must NOT short-circuit to `[0,0]` here; the eta=0 oracle goes through the **full SDP +
  restoration** path and must return `[1,1]`. (Flag: the `‖J‖_F<1e-9·scale` dispatch is
  a `‖Λ‖`-defect heuristic specific to the idempotence-defect use; for the cb-norm-of-`v`
  use it should be removed or its threshold reasoned separately. Escalate if the
  restoration is unstable on the isometric `J(v*)`.)

### 4.3 Other rungs

- **double vs arb@53.** O1 already cross-checks the HOPM gate double-vs-arb to ~1e-15
  (FINDINGS §C12 line 507); O2's SDP value is double (MOSEK) restored in arb, the same
  ladder rung 2 the `aic_cbnorm_certify` tests use.
- **universality canary.** `‖v‖_cb` and `‖v^{-1}‖_cb` must NOT grow with `dim A` nor `n_B`
  (the `.tex:484` failure mode). Sweep dim and assert flat (abs-max + no-upward-trend, the
  robust form FINDINGS §D2 2026-05-31). O1's `a_cb_flat` is the lower-bound analogue.

---

## 5. Increment plan + file split (≤200 LOC/file)

| step | file | LOC | content |
|------|------|-----|---------|
| O2.1 | `src/aic_opspace_choi.c` (NEW) | ~80 | `aic_opspace_choi_v`, `aic_opspace_choi_vinv` (§1.2-1.3), and `aic_opspace_choi_adjoint` = transpose. Reuses `b_unit_pos`, `build_M1`, `acb_mat_solve` (factor the latter two out of `aic_opspace_map.c` into a shared static or a small `aic_opspace_coord.c`). |
| O2.2 | `include/aic_opspace.h` (EDIT) | ~30 | declare the three Choi builders + the O2 result fields (`cb_upper`, `cbinv_upper` arb_t). |
| O2.3 | `julia/.../src/sdp.jl` (EDIT) | ~60 | `diamond_norm_watrous_primal_rect(J,d_in,d_out)`, `diamond_norm_dual_rect(J,d_in,d_out)` (general two-dim, exposing feasible points). Keep the self-map wrappers. |
| O2.4 | `src/aic_cbnorm_certify_rect.c` (NEW) or EDIT `certify.c`/`certify_restore.c` | ~120 | generalized `(d_in,d_out,norm_factor,trace_dir)` certifier; the LOWER/UPPER recipes parameterized (§3.1-3.2). Keep the self-map wrapper `aic_cbnorm_certify`. Split if >200. |
| O2.5 | `tools/gen_fixtures_opspace_o2.jl` (NEW) | ~150 | build `v` for the corpus, form `J(v*)`/`J((v^{-1})*)`, solve both SDPs, assert eta=0 + asymmetric anchors (poison guard), emit `tests/fixtures_opspace_o2.inc.h`. NO PYTHON; the `gen_fixtures_certify.jl` pattern. |
| O2.6 | `src/aic_opspace_o2.c` (NEW) | ~120 | `aic_opspace_certify_cb(result, v, J_feasible_points..., prec)` — assemble `J(v*)`, call the rect certifier, fill `cb_upper`; same for inverse → `cbinv_upper`; assert the HOPM≤SDP bracket (§4.1) + eta=0 oracle. |
| O2.7 | `tests/test_opspace_o2.c` (NEW) | ~150 | (a) eta=0 oracle `lo,hi ∋ 1` to 1e-9 on `block_cond_exp`; (b) eta>0 channel: `‖v‖_cb,‖v^{-1}‖_cb ∈ [1, 1+O(η)]`, asymmetric-`J(v*)` direction tooth; (c) HOPM(O1) ≤ SDP(O2) bracket; mutation-prove each. |

**Build-order note.** O2.5 (the committed fixtures) is generated once with Julia+MOSEK;
thereafter `make test` reads `fixtures_opspace_o2.inc.h` and runs only C (the
`gen_fixtures_certify.jl` discipline). Regenerate with `make fixtures` after touching the
SDP/corpus.

---

## 6. Risks, obstructions, escalations

### 6.1 `v` is HP but NOT CP — does `J(v)` indefiniteness break `src/sdp.jl`?

**No.** The Watrous diamond-norm SDP (`sdp.jl`) is valid for **any linear map**, and the
project already runs it on `Λ=Φ²−Φ`, which is HP-but-**indefinite** (`aic_ucp.h:170-172`:
"the result is Hermitian … but in general NOT PSD … which is exactly why the diamond-norm
SDP is needed"). `J(v*)` is the same regime (Hermitian, indefinite). The MAX/MIN
programs constrain *block* PSD-ness of `[[P,X],[X†,Q]]` and `[[Y0,−J],[−J†,Y1]]`, never
`J ⪰ 0`. So CP is **not** assumed anywhere in `sdp.jl` — safe.

The genuine new hazard is the **rectangular generalization** (§2.2), not CP. The
existing certifier hard-codes `N=n²` and `(n,n)` partial-trace dims; the rect refactor
is where bugs hide.

### 6.2 No-parallel-Julia constraint

This DESIGN pass is read-only (no build, no Julia — a build would race the concurrent
hostile review). The **committed-fixture pattern** (O2.5) is exactly what makes the
eventual `make test` Julia-free, so the steady-state test does not need Julia at all; the
one-time fixture generation runs Julia serially (`gen_fixtures_d24.jl:8` "Run serially
(single Julia process)"). No conflict introduced by the design.

### 6.3 Is Smith truncation MOOT for O2? — YES (the key subtlety)

**Yes, Smith truncation is MOOT for the O2 computation.** The Watrous cb-norm SDP on
`J(v)` computes `‖v‖_cb = sup_k ‖1_{M_k}⊗v‖_op` **without any explicit ampliation**: the
Choi matrix `J(v)` already encodes the complete-bounded structure of `v`, and the SDP's
primal/dual variables live on the FIXED space `M_{d_out}⊗M_{d_in}` (size `n_B·N`),
**independent of `k`**. There is no "loop over ampliation levels" — the SDP returns the
all-`k` supremum intrinsically (this is the entire point of Watrous's SDP characterization
of cb-norms). Contrast O1: the HOPM builds the *explicit* matrix `1_{M_n}⊗v` and must be
told *which* `n` to build — there Smith (forward `N_max=N`, inverse `n_B`) is **load-
bearing** (it says a single finite level suffices and which one).

**Implication for the O2 design:** O2 does **not** sweep ampliation levels, does **not**
build `1_{M_n}⊗v`, and does **not** consult `N_max`. It assembles ONE Choi `J(v*)` (size
`n_B·N`) and solves ONE SDP. Smith's theorem still matters as the *theorem that the O1
HOPM lower bound and the O2 SDP upper bound are bracketing the SAME number* `‖v‖_cb`
(forward `=‖v_N‖_op`, inverse `=‖v_{n_B}^{-1}‖_op`) — i.e. Smith licenses the §4.1
cross-check, but is absent from O2's actual computation.

### 6.4 The eta=0 dispatch mismatch (flag — §4.2)

The existing `aic_cbnorm_certify` short-circuits `‖J‖_F<thr → [0,0]` because for the
idempotence-defect use the trivial value is `0`. For O2 the trivial (η=0) value is
`‖v‖_cb = 1` with `‖J(v*)‖_F = O(1)`. The rect certifier must **not** inherit that
short-circuit; the η=0 case must traverse the full restoration and return `[1,1]`. If the
restoration is numerically unstable on the near-isometric `J(v*)` (the Slater center for
an isometry may be ill-conditioned), **escalate** — that is a precision-wall stop
condition (CLAUDE.md stop conditions).

### 6.5 Normalization is empirically pinned, not derived (§2.4)

`sdp.jl:26-30, 32` are emphatic: the design-doc normalization table was WRONG; only the
analytic anchors + the asymmetric poison-pill are trustworthy. The `(2/N)` factor in §2.4
is a **derivation to be confirmed**, not an assumption. The eta=0 oracle (`‖v‖_cb=1`)
discriminates `2/N` vs `2/n_B` vs `2/dim_A` (they differ unless `N=n_B`). **Pin it on the
oracle; record the measured constant.** If none of the candidate constants hits 1.0
exactly, the `J(v*)` assembly or the transpose convention (§1.4) is wrong — re-derive
before trusting any eta>0 value.

### 6.6 Partial-trace direction (§2.5) — the load-bearing convention

Restating as a standalone risk because it is THE bug that bit the original work
(HANDOFF:340). The Convention-A index layout places `v*`'s **input (N) factor in the
MAJOR/LEFT** position (§1.4) — *opposite* to the self-map `Λ` (input minor/right). So the
rect dual likely needs `partial_trace_left` with marginal dims `(N, n_B)`, NOT
`partial_trace_right(.,n,n)`. **Pin with an asymmetric eta>0 fixture** (the WRONG
direction doubles `hi`); the corpus MUST contain one with an asymmetric `J(v*)` marginal
(the `asym=‖Tr_R−Tr_L‖_F` poison guard of `gen_fixtures_certify.jl:184-188`), else the
direction test is inert.

---

## Appendix — citation index

- Choi Convention A: `include/aic_ucp.h:24-42`; generator `tools/gen_fixtures_d24.jl:37-43`.
- `v` data model / `vE[i]=v(E_i)` / `b_unit_pos`: `include/aic_dhom.h:103-122, 180-184`;
  `src/aic_opspace_map.c:30-44, 46-69`.
- `M_1`, `M_1^{-1}` (coordinate map + inverse): `src/aic_opspace_map.c:46-134`.
- Existing diamond-norm SDP (3 forms) + normalization: `julia/.../src/sdp.jl:38-111`,
  esp. the `(2/n)` factor `sdp.jl:12-16` and the partial-trace-direction pinning
  `sdp.jl:18-30`.
- arb feasibility-restoration certifier: `src/aic_cbnorm_certify.c` (dispatch),
  `src/aic_cbnorm_certify_restore.c` (LOWER `int_lower`, UPPER `int_upper`),
  `include/aic_cbnorm.h:53-74`.
- Committed-fixture pattern: `tools/gen_fixtures_certify.jl` (feasible points + poison
  guards), `tools/gen_fixtures_d24.jl` (golden master).
- O1 contract + honesty (HOPM lower, SDP=O2): `include/aic_opspace.h:55-83`;
  pipeline `src/aic_opspace_cert.c`.
- Smith truncation, adjoint duality, Watrous SDP refs: `docs/research/opspace_web_leg.md`
  §1-2 (Pisier Prop 1.12; Watrous 2009 arXiv:0901.4709 eq.(1); Watrous 2012 arXiv:1207.5726
  Thm 6; QETLAB `CBNorm`).
- `th_main_ext` statement / `a_n` / `def:opspace`: `.tex:1453-1464, 1487-1505, 1538-1540,
  1542, 1557`; `docs/research/opspace_paper_leg.md` §1, §4.
- FINDINGS §C12 (operator-norm not Frobenius), §D3 (Smith): `paper/FINDINGS.md:488-617`.
