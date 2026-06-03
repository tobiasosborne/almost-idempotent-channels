# AlmostIdempotentChannels.jl — API Inventory and Docstring Audit

Generated from a read-only audit of all source files under
`julia/AlmostIdempotentChannels.jl/src/` and `ext/AICMosekExt.jl`.
Authoritative export list: `src/AlmostIdempotentChannels.jl` lines 84–107.

---

## 1. Public surface map

**Summary counts:** 44 exported names total — 7 value-types, 19 accessors/type-methods,
8 pipeline/defect verbs, 3 low-level ccall-surface names, 5 MOSEK-gated names,
2 library-discovery names.

### 1.1 Defect verbs

| Name | Full signature | Return type | Purpose | Paper result |
|------|---------------|-------------|---------|--------------|
| `certified_defect` | `certified_defect(Φ::UCPMap; prec::Int=106, tight::Bool=false)` | `CertifiedBracket` | Rigorous two-sided bracket on η = ‖Φ²−Φ‖_cb. `tight=false`: solver-free eig-free arb certifier (loose, ~2n width). `tight=true`: MOSEK-tight (extension-only). Always safe — no basin pre-check needed. | tex:347; th_factorization |
| `idempotency_defect` | `idempotency_defect(Φ::UCPMap; prec::Int=106)` | `Float64` | Exact η value via Watrous diamond-norm SDP. Extension-only (MOSEK); throws helpful error otherwise. | tex:347-354 |
| `idempotency_defect` | `idempotency_defect(kraus::Vector{Matrix{ComplexF64}}; prec::Int=106)` | `Float64` | Back-compat raw-Kraus overload (pre-UCPMap surface; exercises the same SDP path). | tex:347-354 |
| `eta_idempotence` | `eta_idempotence(kraus::Vector{Matrix{ComplexF64}})` | `Float64` | Back-compat alias of `idempotency_defect`. Extension-only. | tex:347-354 |
| `choi_diff` | `choi_diff(kraus::Vector{Matrix{ComplexF64}}, n::Int; prec::Int=106)` | `Matrix{ComplexF64}` | Choi(Φ²−Φ), Convention-A, n²×n², via arb core `aic_ucp_choi_diff_d`. | tex:347-354 |
| `eta_eigfree` | `eta_eigfree(kraus::Vector{Matrix{ComplexF64}}; prec::Int=106)` | `Tuple{Float64,Float64}` | Solver-free two-sided bracket (lo, hi) on η via `aic_cbnorm_eigfree_d`. No SDP. | tex:347-354 |

### 1.2 Pipeline verbs

| Name | Full signature | Return type | Purpose | Paper result |
|------|---------------|-------------|---------|--------------|
| `associated_algebra` | `associated_algebra(Φ::UCPMap; prec::Int=256)` | `EpsCStarAlgebra` | The associated ε-C* algebra A = Img Φ̃. Basin-guarded (throws `ArgumentError` if ρ(Φ²−Φ) ≥ 0.24). | tex:2184; th_almost_idemp |
| `main_isomorphism` | `main_isomorphism(Φ::UCPMap; prec::Int=256)` | `MainIsomorphism` | O(ε)-isomorphism v: B → A, with universal (dim-independent) constant. Basin-guarded. | tex:460; th_main |
| `factorize` | `factorize(Φ::UCPMap; prec::Int=256)` | `ChannelFactorization` | Approximate channel factorization Φ ≈ Δ Υ through a genuine C* algebra B. Tighter domain guard (_FACTORIZE_RHO_MAX = 0.10) than basin. Extends `LinearAlgebra.factorize`. | tex:2730-2899; th_factorization |

### 1.3 Per-type channel accessors

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `encode` | `encode(F::ChannelFactorization)` | `CPMap` | Enc = Υ* : B → B(H) (CPTP, dim_in=n_B, dim_out=N). |
| `decode` | `decode(F::ChannelFactorization)` | `CPMap` | Dec = Δ* : B(H) → B (CPTP, dim_in=N, dim_out=n_B). |
| `choi` | `choi(Φ::UCPMap)` | `Matrix{ComplexF64}` | Convention-A Choi of Φ itself (not Φ²−Φ). Pure Julia. |
| `iscptp` | `iscptp(c::CPMap; atol::Real=1e-9)` | `Bool` | Whether c is trace-preserving to atol. |

### 1.4 The seven value-types

| Name | Constructor summary |
|------|-------------------|
| `CertifiedBracket{T<:AbstractFloat}` | `CertifiedBracket(lo, hi; value=nothing, label="")` |
| `UCPMap` | `UCPMap(kraus; check=true, atol=1e-9)` or `UCPMap(; choi, n, check=true, atol=1e-9)` |
| `CPMap` | `CPMap(kraus, dim_in, dim_out; check=true, atol=1e-9)` |
| `EpsCStarAlgebra` | Internal struct; built only by `associated_algebra` |
| `CStarAlgebra` | `CStarAlgebra(blocks::AbstractVector{<:Integer})` |
| `MainIsomorphism` | Internal struct; built only by `main_isomorphism` |
| `ChannelFactorization` | Internal struct; built only by `factorize` |

### 1.5 Per-type accessors

**UCPMap accessors:**

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `n` | `n(Φ::UCPMap)` | `Int` | Hilbert space dimension (each Kraus op is n×n). |
| `nkraus` | `nkraus(Φ::UCPMap)` | `Int` | Number r of Kraus operators. |
| `kraus` | `kraus(Φ::UCPMap)` | `Vector{Matrix{ComplexF64}}` | The Kraus operators {K_a} (Heisenberg picture). |
| `isunital` | `isunital(Φ::UCPMap; atol::Real=1e-9)` | `Bool` | Whether ‖Σ K_a†K_a − I‖ ≤ atol. |
| `Base.adjoint` | `adjoint(Φ::UCPMap)` (as `Φ'`) | `UCPMap` | Dual map with Kraus {K_a†}; built with `check=false` (trace-preserving). **Extended, not separately exported.** |

**CertifiedBracket accessors (Base extensions plus package names):**

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `width` | `width(b::CertifiedBracket)` | `AbstractFloat` | hi − lo (bracket looseness). |
| `midpoint` | `midpoint(b::CertifiedBracket)` | `AbstractFloat` | (lo + hi) / 2. |
| `value` | `value(b::CertifiedBracket)` | `Union{AbstractFloat,Nothing}` | Point estimate (SDP value) or nothing. |
| `Base.in` | `in(x::Real, b::CertifiedBracket)` | `Bool` | lo ≤ x ≤ hi. **Base extension, not separately exported.** |
| `Base.minimum` | `minimum(b::CertifiedBracket)` | `AbstractFloat` | b.lo. **Base extension.** |
| `Base.maximum` | `maximum(b::CertifiedBracket)` | `AbstractFloat` | b.hi. **Base extension.** |

**CStarAlgebra accessors:**

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `blocks` | `blocks(B::CStarAlgebra)` | `Vector{Int}` | Block sizes d_l of B = ⊕_l M_{d_l}. |
| `dim_B` | `dim_B(B::CStarAlgebra)` | `Int` | Vector-space dimension Σ d_l². |
| `n_B` | `n_B(B::CStarAlgebra)` | `Int` | Block-diagonal representative size Σ d_l. |

**EpsCStarAlgebra accessors:**

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `dim` | `dim(A::EpsCStarAlgebra)` | `Int` | dim A (the ε-C* algebra dimension). |
| `ambient` | `ambient(A::EpsCStarAlgebra)` | `Int` | Ambient n (A ≤ M_n). |
| `epsilon` | `epsilon(A::EpsCStarAlgebra)` | `CertifiedBracket` | Certified associator-defect ε bracket. |

**MainIsomorphism accessors:**

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `blocks` | `blocks(v::MainIsomorphism)` | `Vector{Int}` | Block sizes d_l of target algebra B. |
| `cbnorm_forward` | `cbnorm_forward(v::MainIsomorphism)` | `CertifiedBracket` | Certified bracket on ‖v‖_cb. |
| `cbnorm_inverse` | `cbnorm_inverse(v::MainIsomorphism)` | `CertifiedBracket` | Certified bracket on ‖v⁻¹‖_cb. |
| `isodefect` | `isodefect(v::MainIsomorphism)` | `CertifiedBracket` | Certified isomorphism-defect bracket. |

**ChannelFactorization accessors:**

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `algebra` | `algebra(F::ChannelFactorization)` | `CStarAlgebra` | The genuine C* algebra B = ⊕_l M_{d_l} the factorization passes through. |
| `delups_defect` | `delups_defect(F::ChannelFactorization)` | `CertifiedBracket` | Certified bracket on ‖ΔΥ − Φ‖_cb (tex:2733). |
| `upsdel_defect` | `upsdel_defect(F::ChannelFactorization)` | `CertifiedBracket` | Certified bracket on ‖ΥΔ − 1_B‖_cb (tex:2739). |

### 1.6 Library-discovery names

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `libaic_path` | `libaic_path()` | `String` | Absolute path to the loaded `libaic.so` (Preferences value or default). |
| `set_libaic_path!` | `set_libaic_path!(p)` | `Nothing` | Persist a new libaic.so path as a LocalPreferences.toml entry. Restart required. |

### 1.7 MOSEK-gated names (exported from core as helpful-error stubs)

| Name | Full signature | Return type | Purpose |
|------|---------------|-------------|---------|
| `diamond_norm_watrous` | `diamond_norm_watrous(J::Matrix{ComplexF64}, n::Int)` | `(Float64, String)` | Watrous 2012 MAX-primal SDP value + status. Extension-only. |
| `diamond_norm_watrous_primal` | `diamond_norm_watrous_primal(J::Matrix{ComplexF64}, n::Int)` | `(Float64, Matrix, Matrix, Matrix, String)` | MAX-primal SDP exposing feasible point (X,P,Q) + eta + status. Extension-only. |
| `diamond_norm_dual` | `diamond_norm_dual(J::Matrix{ComplexF64}, n::Int)` | `(Float64, Matrix, Matrix, String)` | MIN-dual SDP exposing (Y0,Y1) for upper-bound certifier + eta_dual + status. Extension-only. |

Note: `idempotency_defect` and `certified_defect(tight=true)` are also solver-gated (see section 3).

### 1.8 Internal names (NOT exported — noted for documentation completeness)

These live in `ccall.jl` / `ccall_factorize.jl` and are underscore-prefixed module-private
wrappers. The docsite should NOT expose them. They are:
`_assoc_summary`, `_main_iso_summary`, `_factorize_artifacts`,
`_cbnorm_certify`, `_rho_defect_estimate`, `_assert_in_basin`,
`_assert_factorize_domain`, `_flatten_kraus`, `_unflatten_kraus`,
`_flatten_NN`, `_diamond_value_impl`, `_tight_bracket_impl`.

---

## 2. The seven value-types in detail

### 2.1 `CertifiedBracket{T<:AbstractFloat}`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `lo` | `T` | Rigorous lower bound (arb ball rounded FLOOR to double). |
| `hi` | `T` | Rigorous upper bound (arb ball rounded CEIL to double). |
| `value` | `Union{T,Nothing}` | Optional point estimate (e.g. MOSEK SDP value); `nothing` for solver-free brackets. |
| `label` | `String` | Display name of the quantity ("‖Φ²−Φ‖_cb", "ε_assoc", "‖v‖_cb", …). |

**Constructor signatures:**

```julia
CertifiedBracket(lo::T, hi::T; value=nothing, label::AbstractString="") where {T<:AbstractFloat}
CertifiedBracket(lo::Real, hi::Real; value=nothing, label::AbstractString="")  # promotes to float
```

**Validation:** lo ≤ hi (else `ArgumentError`). When value is given, lo-1e-9 ≤ value ≤ hi+1e-9
(1e-9 slack absorbs benign double-rounding of an SDP point estimate against an outward-rounded bracket).

**Accessors (exported):** `width`, `midpoint`, `value`.

**Base extensions (not separately exported):** `Base.in`, `Base.minimum`, `Base.maximum`.

**Show output shapes:**

Compact (one-line, used in arrays/interpolation):
```
[0.0840, 0.2700]              # no value
[0.0840, 0.1350, 0.2700]     # with value
```

Multi-line (REPL `text/plain`):
```
CertifiedBracket  ‖Φ²−Φ‖_cb
  0.0840  ≤  x  ≤  0.2700        (width 0.1860)
  value = 0.1350
```

---

### 2.2 `UCPMap`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `kraus` | `Vector{Matrix{ComplexF64}}` | r operators {K_a}, each n×n (Heisenberg picture: Φ(X) = Σ K_a† X K_a). |
| `n` | `Int` | Hilbert space dimension. |
| `r` | `Int` | Number of Kraus operators. |

**Constructor signatures:**

```julia
UCPMap(kraus::AbstractVector{<:AbstractMatrix}; check::Bool=true, atol::Real=1e-9)
UCPMap(; choi::AbstractMatrix, n::Integer, check::Bool=true, atol::Real=1e-9)
```

**Validation (Kraus constructor):**
- Empty list: `ArgumentError` "UCPMap: empty Kraus list (need ≥ 1 operator)".
- Non-square K_a: `DimensionMismatch`.
- Ragged K_a (differing sizes): `DimensionMismatch`.
- Non-unital beyond atol (when check=true): `ArgumentError` naming ‖Σ K†K − I‖.

**Choi constructor** reconstructs Kraus via PSD eigendecomposition of the
Convention-A Choi matrix `C[(i-1)n+a, (j-1)n+b] = Σ_x conj(K_x[i,a]) K_x[j,b]`.
Tolerance-aware (drops near-zero/negative eigenvalues). Then validates unitality
identically to the Kraus path (with the same `check`/`atol` kwargs).

**Accessors (exported):** `n`, `nkraus`, `kraus`, `isunital`.

**Base extension (not separately exported):** `Base.adjoint` (as `Φ'`).
The dual is built with `check=false` (its certificate is trace-preservation,
not Heisenberg unitality).

**Show output shapes:**

Compact:
```
UCPMap(n=2, r=4)
```

Multi-line:
```
UCPMap: Φ: B(H) → B(H),  dim H = 2
  Kraus operators: 4   (Heisenberg: Φ(X) = Σ K_a† X K_a)
  unital: yes   (‖Σ K†K − I‖ = 2.78e-17)
```
If non-unital (e.g. a Dec/Enc dual built with `check=false`):
```
  unital: no    (‖Σ K†K − I‖ = 1.94e-06;  trace-preserving dual?)
```

---

### 2.3 `CPMap`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `kraus` | `Vector{Matrix{ComplexF64}}` | Kraus ops K_a: each dim_out × dim_in (Schrödinger: ρ ↦ Σ K_a ρ K_a†). |
| `dim_in` | `Int` | Input space dimension. |
| `dim_out` | `Int` | Output space dimension. |

**Constructor signature:**
```julia
CPMap(kraus::AbstractVector{<:AbstractMatrix}, dim_in::Integer, dim_out::Integer;
      check::Bool=true, atol::Real=1e-9)
```

**Validation:**
- Empty list: `ArgumentError`.
- dim_in < 1 or dim_out < 1: `ArgumentError`.
- Wrong-shaped K_a (not dim_out × dim_in): `DimensionMismatch`.
- Not trace-preserving beyond atol (when check=true): `ArgumentError` naming ‖Σ K_a†K_a − I_{dim_in}‖.

**Key distinction from UCPMap:** Rectangular (dim_in ≠ dim_out in general).
Not named `Channel` (would collide with `Base.Channel`).

**Dual-direction binding for factorize:**
- `encode(F)` = Enc = Υ* : B → B(H), `dim_in = n_B`, `dim_out = N` (CPTP, n_B→N).
- `decode(F)` = Dec = Δ* : B(H) → B, `dim_in = N`, `dim_out = n_B` (CPTP, N→n_B).
These are NOT crossed. The dims are asserted in tests.

**Show output shapes:**

Compact:
```
CPMap(2 → 4, 8 Kraus)
```

Multi-line (three cases based on measured TP defect):
```
CPMap: CPTP  C^2 → C^4   (Schrödinger: ρ ↦ Σ K_a ρ K_a†)
  Kraus operators: 8   (each 4×2)
  trace-preserving: yes        (‖Σ K_a† K_a − I_2‖ = 3.14e-17)
```
For an η>0 `decode` with measured defect in (1e-9, 1e-4]:
```
  trace-preserving: approx O(η) (‖Σ K_a† K_a − I_4‖ = 1.90e-06;  the cb-norm round-trip bracket is the certificate)
```
For clearly non-TP:
```
  trace-preserving: no         (‖Σ K_a† K_a − I_4‖ = 0.4300)
```

---

### 2.4 `EpsCStarAlgebra`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `n` | `Int` | Ambient dimension (A is an oblique subspace of M_n). |
| `dim_A` | `Int` | Dimension of the ε-C* algebra. |
| `basis` | `Vector{Matrix{ComplexF64}}` | Frobenius-orthonormal snapshot {B_k} (n×n each). **Oblique**: the algebra product is the Choi-Effros star X⋆Y = Φ̃(XY), NOT plain XY. Do not multiply basis elements directly in Julia. |
| `eps` | `CertifiedBracket{Float64}` | Certified associator-defect ε bracket (how far A is from a genuine C* algebra). |
| `source` | `UCPMap` | The originating channel. |

**Constructor:** Internal; built only by `associated_algebra`. No public constructor.

**Accessors (exported):** `dim`, `ambient`, `epsilon`.

**Show output shapes:**

Compact:
```
EpsCStarAlgebra(dim A=4 ≤ M_2)
```
Multi-line:
```
EpsCStarAlgebra  A = Img Φ̃  ≤  M_2
  dim A = 4
  ε (associator defect) ∈ [3.12e-04, 8.40e-02]
  product: Choi–Effros star  X⋆Y = Φ̃(XY)   (oblique; axioms hold up to ε)
```

---

### 2.5 `CStarAlgebra`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `blocks` | `Vector{Int}` | Block sizes d_l of B = ⊕_l M_{d_l}. |
| `dim_B` | `Int` | Vector-space dimension Σ d_l² (computed by constructor). |
| `n_B` | `Int` | Block-diagonal representative size Σ d_l (computed by constructor). |

**Constructor signature:**
```julia
CStarAlgebra(blocks::AbstractVector{<:Integer})
```
Validates: non-empty (else `ArgumentError`); all d_l ≥ 1 (else `ArgumentError`).

**Accessors (exported):** `blocks`, `dim_B`, `n_B`.

**Show output shapes:**

Compact:
```
⊕ M_d [1, 1]
```
Multi-line:
```
CStarAlgebra  B = ⊕_l M_{d_l}
  blocks d_l = [1, 1]   (m = 2)
  dim_B = 2,  n_B = 2
```

---

### 2.6 `MainIsomorphism`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `source` | `UCPMap` | The originating channel. |
| `B` | `CStarAlgebra` | The genuine finite-dim C* algebra B = ⊕_l M_{d_l}. |
| `cbnorm_fwd` | `CertifiedBracket{Float64}` | Certified bracket on ‖v‖_cb. |
| `cbnorm_inv` | `CertifiedBracket{Float64}` | Certified bracket on ‖v⁻¹‖_cb (≈ 1 + O(ε)). |
| `isodefect` | `CertifiedBracket{Float64}` | Certified isomorphism-defect bracket. |

**Constructor:** Internal; built only by `main_isomorphism`. No public constructor.

**Accessors (exported):** `blocks` (delegates to B), `cbnorm_forward`, `cbnorm_inverse`, `isodefect`.

**Show output shapes:**

Compact:
```
MainIsomorphism(B = ⊕ M_d [1, 1])
```
Multi-line:
```
MainIsomorphism  v: B → A   (O(ε), dimension-independent constant)
  B = ⊕ M_d, blocks = [1, 1]
  ‖v‖_cb    ∈ [1.00, 1.04]
  ‖v⁻¹‖_cb  ∈ [1.00, 1.04]
  iso defect ∈ [2.10e-12, 4.20e-12]
```

---

### 2.7 `ChannelFactorization`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `source` | `UCPMap` | The originating channel. |
| `B` | `CStarAlgebra` | The genuine C* algebra B the factorization passes through. |
| `encode` | `CPMap` | Enc = Υ* : B → B(H), `dim_in=n_B`, `dim_out=N`. |
| `decode` | `CPMap` | Dec = Δ* : B(H) → B, `dim_in=N`, `dim_out=n_B`. |
| `delups` | `CertifiedBracket{Float64}` | ‖ΔΥ − Φ‖_cb (tex:2733). |
| `upsdel` | `CertifiedBracket{Float64}` | ‖ΥΔ − 1_B‖_cb (tex:2739). |
| `eta_proxy` | `Float64` | Sentinel from C pipeline (≈ 0 at η = 0). |
| `eps_used` | `Float64` | The ε used by the C pipeline (read back from shim). |

**Constructor:** Internal; built only by `factorize`. No public constructor.

**Accessors (exported):** `algebra`, `delups_defect`, `upsdel_defect`.
Channels accessed via `encode(F)` / `decode(F)` (exported functions, not field names).

**Dual-direction binding (LOAD-BEARING — Appendix B1, FINDINGS §C13):**
`encode(F)` = Enc = Υ* (B → B(H), n_B→N); `decode(F)` = Dec = Δ* (B(H) → B, N→n_B).
The physical names ("encode" = "compress into B", "decode" = "expand from B")
bind to the OBSERVABLE duals (Υ*/Δ*), NOT to the observable maps Υ/Δ.
"Encode" is NEVER placed next to "Dec=Δ*" in display or docs.

**Show output shapes:**

Compact:
```
ChannelFactorization(Φ ≈ Δ Υ, B = ⊕ M_d [1, 1])
```
Multi-line (Kraus counts right-justified to align columns):
```
ChannelFactorization  Φ ≈ Δ Υ  through  B = ⊕ M_d, blocks = [1, 1]
  encode = Υ* : B → B(H)   (4 Kraus,  2→4 CPTP)
  decode = Δ* : B(H) → B   (4 Kraus,  4→2 CPTP)
  ‖ΔΥ − Φ‖_cb   ∈ [0.00, 4.20e-12]     (round-trip ΔΥ ≈ Φ)
  ‖ΥΔ − 1_B‖_cb ∈ [0.00, 4.20e-12]     (round-trip ΥΔ ≈ 1_B)
  η proxy = 0.0,  ε used = 0.0
```

---

## 3. Solver-gated behaviour

### 3.1 Names requiring the MOSEK extension

| Name | Gating mechanism |
|------|-----------------|
| `idempotency_defect(Φ::UCPMap; ...)` | Routes to `_diamond_value_impl`; core stub throws. |
| `idempotency_defect(kraus::Vector{Matrix{ComplexF64}}; ...)` | Same. |
| `eta_idempotence(kraus)` | Alias of `idempotency_defect`; same gating. |
| `certified_defect(Φ; tight=true)` | Routes to `_tight_bracket_impl`; core stub throws. |
| `diamond_norm_watrous(J, n)` | Direct stub in `sdp_stubs.jl`. |
| `diamond_norm_watrous_primal(J, n)` | Direct stub in `sdp_stubs.jl`. |
| `diamond_norm_dual(J, n)` | Direct stub in `sdp_stubs.jl`. |

**Always safe (NO solver needed):**
`certified_defect(Φ)` (default `tight=false`), `choi_diff`, `eta_eigfree`,
`associated_algebra`, `main_isomorphism`, `factorize`.

### 3.2 The helpful-error mechanism

**Pattern (design doc §3.3, "register-by-method-override"):**

The core package defines catch-all fallbacks in `sdp_stubs.jl`:
```julia
_diamond_value_impl(::Any, ::Int) = error("idempotency_defect " * _MOSEK_HINT)
_tight_bracket_impl(::Any, ::Int) = error("certified_defect(Φ; tight=true) " * _MOSEK_HINT)
diamond_norm_watrous(::Any, ::Any) = error("diamond_norm_watrous " * _MOSEK_HINT)
diamond_norm_watrous_primal(::Any, ::Any) = error("diamond_norm_watrous_primal " * _MOSEK_HINT)
diamond_norm_dual(::Any, ::Any) = error("diamond_norm_dual " * _MOSEK_HINT)
```

The constant `_MOSEK_HINT` is:
```
requires the MOSEK extension. Install + load the solver:
    using Pkg; Pkg.add(["Convex","Mosek","MosekTools"])
    using Convex, Mosek, MosekTools
(Mosek needs a license: https://www.mosek.com/products/academic-licenses/)
```

Additionally, `__init__` registers a `Base.Experimental.register_error_hint` callback
for `MethodError` that appends `_MOSEK_HINT` whenever a `MethodError` is thrown on
any of the gated functions — belt-and-suspenders, since the stubs above throw `error`
(not `MethodError`).

When `AICMosekExt` is loaded, the extension adds concrete methods on the specific types
(`UCPMap`, `Matrix{ComplexF64}`, `Int`) which are strictly MORE SPECIFIC than the
catch-all stubs and win via Julia's normal method table, without overwriting the stubs.

### 3.3 MOSEK extension's additional internal names (not exported)

The extension adds two functions used only internally for the O2 (rectangular) diamond norm
in tooling/probing contexts (not part of the documented public API):
- `diamond_dual_rect(J, d_maj, d_min, tr_sys)` — rectangular MIN-dual SDP.
- `diamond_primal_rect(J, d_maj, d_min, rho_on)` — rectangular MAX-primal SDP.

These are NOT exported and should NOT appear in the API reference.

### 3.4 MOSEK size bound and graceful degradation

The `_tight_bracket_impl` path has a hard bound `_TIGHT_DUAL_NMAX = 5`:
- For `n > 5`: the MIN-dual SDP stalls/times out (aic-jhe / bug D6). The extension
  falls back to the eig-free bracket silently (with a `@warn`), never hangs.
- For `n ≤ 5` with a non-converged SDP (status not OPTIMAL/SLOW_PROGRESS): same
  eig-free fallback with `@warn`.
- The diamond VALUE (`idempotency_defect`) uses only the MAX-primal SDP (no dual),
  so it does not hit the `_TIGHT_DUAL_NMAX` bound.

---

## 4. Domain / safety preconditions

### 4.1 `certified_defect`

No domain restriction. Always safe. No basin pre-check needed. Returns a `CertifiedBracket`.

### 4.2 `associated_algebra`

**Pre-check location:** `basin.jl:_assert_in_basin`.

**Method:** Power-iteration estimate of ρ(Φ²−Φ) (80 iterations, deterministic Hermitian seed,
no n²×n² superoperator formed). Threshold constant: `_BASIN_RHO_MAX = 0.24` (just below 1/4).

**Exception type and message on failure:**
```
ArgumentError: "Φ is out of the prop_P basin: estimated ρ(Φ²−Φ) ≈ <ρ> ≥ 1/4 " *
    "(approximate_algebras.tex:520-525, FINDINGS §C15). The exact-idempotent " *
    "regularisation Φ̃ = θ(2Φ−1) is undefined there, and the C pipeline would abort. " *
    "Use certified_defect(Φ) (always safe, eig-free) to inspect the idempotency defect instead."
```

### 4.3 `main_isomorphism`

**Pre-check:** Same `_assert_in_basin` / `_BASIN_RHO_MAX = 0.24` as `associated_algebra`.
Identical exception message.

### 4.4 `factorize`

**Pre-check location:** `basin.jl:_assert_factorize_domain`.

**Method:** Same power-iteration estimate of ρ(Φ²−Φ). Threshold constant:
`_FACTORIZE_RHO_MAX = 0.10` (tighter than the prop_P basin; determined by empirical
sweep showing C-abort boundary at ρ_est ≈ 0.1275 across d∈{2,3,4}).

**Exception type and message on failure:**
```
ArgumentError: "Φ is outside factorize's convergence domain: estimated ρ(Φ²−Φ) ≈ " *
    "<ρ> ≥ 0.10. factorize builds Υ via a " *
    "power-series functional calculus (aic_funcalc_xpow inside upsilon_build) " *
    "whose convergence domain is TIGHTER than the prop_P basin ρ < 1/4 (bug " *
    "aic-exa.13): out there the C pipeline ABORTS the session. ..."
```

**Best-effort caveat:** The guard is best-effort; a severely out-of-domain channel of an
untested shape (where the estimate under-reads) could still abort. The error message
names this explicitly.

### 4.5 Domain summary

| Verb | Domain | Guard threshold | Exception on failure |
|------|--------|----------------|---------------------|
| `certified_defect` | Unrestricted | — | — |
| `choi_diff` / `eta_eigfree` | Unrestricted | — | — |
| `associated_algebra` | ρ(Φ²−Φ) < 1/4 | `_BASIN_RHO_MAX = 0.24` | `ArgumentError` |
| `main_isomorphism` | ρ(Φ²−Φ) < 1/4 | `_BASIN_RHO_MAX = 0.24` | `ArgumentError` |
| `factorize` | small-η only | `_FACTORIZE_RHO_MAX = 0.10` | `ArgumentError` |

---

## 5. Docstring audit

Priority: **(A) MISSING** > **(B) thin/no example** > **(C) missing tex citation** > **(D) inconsistency**.

### Priority A — MISSING docstrings

**A1. `EpsCStarAlgebra` constructor (public fields and how to obtain it)**

The struct docstring (types_algebra.jl:14-26) documents the type, but the
constructor is internal and there is no docstring explaining that the type is
ONLY obtainable via `associated_algebra`. The field `basis` is mentioned but
the critical warning "do not multiply basis elements directly in Julia" is only
in an inline comment, not in the docstring. Users who encounter the type in REPL
output need to understand they cannot construct it directly.

**What to add:** A note in the struct docstring: "Constructed only via
`associated_algebra(Φ)`; no public constructor. The `basis` field is a snapshot
for display; the algebra's operations (Choi-Effros star, tex:342) go through the
C pipeline, not through direct Julia multiplication."

---

**A2. `MainIsomorphism` constructor (same issue as A1)**

The struct docstring (types_algebra.jl:108-115) describes the type but does not
state it is only obtainable via `main_isomorphism`. The `source` field is listed
without explanation of what it is used for (it exists to allow downstream
routines to re-derive the isomorphism from the same Φ).

**What to add:** "Constructed only via `main_isomorphism(Φ)`; no public
constructor. The `source` field allows downstream arb certifiers to rebuild the
isomorphism from the same Φ with a consistent `eps_used`."

---

**A3. `ChannelFactorization` constructor (same issue as A1/A2)**

Same pattern: struct docstring exists but no statement that it is only
obtainable via `factorize`.

**What to add:** "Constructed only via `factorize(Φ)`; no public constructor."

---

**A4. `iscptp` — missing a tex or ALGORITHM.md citation**

The `iscptp` docstring (types_channel.jl:63-79) is thorough on the O(η)-TP
caveat (FINDINGS §C14) but carries no tex citation for the CPTP condition
(tex:334 Enc/Dec as UCP duals, or the general CPTP definition). The docstring
is correct but would benefit from at least one citation.

---

**A5. `isunital` — missing a tex citation**

The `isunital` docstring (types_core.jl:183-186) is functionally complete but
carries no tex citation (Heisenberg unitality Σ K_a†K_a = I appears at tex:8
in aic_ucp.h; the paper's definition of UCP self-map is tex:347-351). Add at
minimum a reference to the paper's UCP self-map definition.

---

**A6. `choi` — note about Convention-A placement**

The `choi` docstring (api.jl:188-195) correctly notes "not of Φ²−Φ; the latter
is the low-level `choi_diff`" but does not clarify WHY the convention matters
for downstream use (the Watrous SDP normalization factor `2/n` depends on the
Choi having trace n for a UCP map). A one-sentence note would close this.

---

**A7. `libaic_path` and `set_libaic_path!` — missing example**

Both have correct docstrings but no usage example. `set_libaic_path!` especially
needs a one-liner showing the typical out-of-tree build case:
```julia
AlmostIdempotentChannels.set_libaic_path!("/path/to/build/libaic.so")
# restart Julia to apply
```

---

### Priority B — thin or missing examples

**B1. `certified_defect` — no example showing both `tight=false` and `tight=true`**

The docstring explains both paths clearly but provides no concrete usage example.
A user inspecting the REPL help needs to see:
```julia
Φ = UCPMap(phit_kraus(2, 0.05))
b = certified_defect(Φ)          # -> CertifiedBracket, [lo, hi] on eta
value(b)                          # -> nothing (solver-free)
b2 = certified_defect(Φ; tight=true)  # (requires MOSEK extension)
value(b2)                         # -> the SDP point estimate
```

---

**B2. `UCPMap` — no example of the Choi constructor**

The Kraus-based constructor is the common path, but the Choi constructor
`UCPMap(; choi=C, n=2)` is undiscoverable without an example. It is the
round-trip from `choi(Φ)` back to a `UCPMap` and deserves one concrete example.

---

**B3. `CStarAlgebra` — no example showing derived fields**

The constructor docstring is minimal. No example showing that
`CStarAlgebra([2,1])` gives `dim_B=5`, `n_B=3`. Trivial to add; useful for
users writing tests.

---

**B4. `factorize` — the domain caveat is thorough but has no positive example**

The docstring is the longest in the codebase and focuses heavily on failure modes
(out-of-domain) while omitting a single happy-path usage example. Users need:
```julia
F = factorize(Φ)
encode(F)         # -> CPMap (Υ*)
decode(F)         # -> CPMap (Δ*)
delups_defect(F)  # -> CertifiedBracket on ‖ΔΥ−Φ‖_cb
upsdel_defect(F)  # -> CertifiedBracket on ‖ΥΔ−1_B‖_cb
```

---

**B5. `associated_algebra` and `main_isomorphism` — no examples**

Both have correct mathematical docstrings but no usage example. In particular,
the user cannot see without an example that `associated_algebra` gives them
`dim(A)`, `epsilon(A)`, and a `basis` snapshot.

---

**B6. `eta_eigfree` — the "loose by design" warning buries the use case**

The docstring says "the bracket is loose by design (ratio hi/lo ~ 2n); it
certifies the MOSEK value rather than competing with it." This is accurate but
could mislead a user into thinking the bracket is useless. A concrete example
showing a typical bracket width at n=2 (e.g. hi/lo ≈ 4) would calibrate
expectations.

---

### Priority C — missing tex citations

**C1. `associated_algebra` — cites tex:2184 (th_almost_idemp) but not the ε axioms**

The docstring says "the axioms hold only UP TO ε (no exact unit)" but doesn't
cite the ε-C* axioms (tex:407-439) or the exact-unit proposition (tex:672,
prop_unit). A reader wanting to understand the precise sense of "ε-C*" needs
these.

---

**C2. `main_isomorphism` — cites tex:460 (th_main) but not the FINDINGS item for the universal constant**

The docstring says "UNIVERSAL, dimension-INDEPENDENT constant (FINDINGS §D2 —
the constant must NOT grow with n)" but does not cite the tex line that states
the constant is universal. The exact statement is at tex:460-461. The FINDINGS
cross-reference is good; add the tex citation.

---

**C3. `factorize` — tex:2730-2899 range is correct but should cite the two specific theorems**

The relevant results are `th_factorization` (tex:2730) and the two round-trip
bounds (tex:2733 and tex:2739). The docstring cites the range but only uses the
name `th_factorization`. The round-trip bound lines should be cited individually
since they are the docstring's primary correctness claim.

---

**C4. `delups_defect` and `upsdel_defect` — cite tex lines**

Both accessor docstrings say "tex:2733" and "tex:2739" respectively — they
already have citations. No action needed.

---

**C5. `diamond_norm_watrous` stub — correct but thin**

The stub docstring says "Watrous 2012 MAX-primal diamond-norm SDP for the map
with Convention-A Choi `J`." It does not cite arXiv:1207.5726. The extension
comment block does cite it, but the exported docstring does not. The docstring
is what the user sees via `?diamond_norm_watrous`; add the citation.

---

### Priority D — inconsistencies / signature drift / stale claims

**D1. `idempotency_defect` — two docstrings for the same exported name; the stub version is stale**

The stub in `sdp_stubs.jl:28-38` has a docstring describing the raw-Kraus
signature `(kraus::Vector{Matrix{ComplexF64}})`. The method in `api.jl:57-64`
has a docstring describing the UCPMap signature `(Φ::UCPMap)`. Both appear
under `?idempotency_defect`. The stub docstring no longer accurately describes
the primary (UCPMap) use case; a user reading it first will be confused. The
stub docstring should be updated or its note redirected to the UCPMap form.

**What to fix:** Either (a) remove the docstring from the stub and let the
api.jl docstring be the primary, or (b) retitle the stub docstring as "raw-Kraus
back-compat overload" and cross-reference the UCPMap method.

---

**D2. `eta_idempotence` docstring — says "Back-compat alias" but doesn't say of what signature**

The docstring (AlmostIdempotentChannels.jl:194-202) says "Back-compat alias of
`idempotency_defect`." It does not say which overload it aliases (the raw-Kraus
`Vector{Matrix{ComplexF64}}` one, NOT the UCPMap one). The signature shown in
the docstring `(kraus::Vector{Matrix{ComplexF64}})` is correct, but the "alias"
statement is incomplete. A user might try `eta_idempotence(Φ::UCPMap)` and get
a method error.

**What to fix:** "Back-compat alias of `idempotency_defect(kraus::Vector{Matrix{ComplexF64}})`
— takes the raw Kraus list, not a `UCPMap`. Use `idempotency_defect(Φ::UCPMap)` for the
idiomatic path."

---

**D3. `choi_diff` — takes `kraus::Vector{Matrix{ComplexF64}}` not a `UCPMap`**

The function signature `choi_diff(kraus::Vector{Matrix{ComplexF64}}, n::Int; prec=106)`
is a lower-level entry point. The docstring is correct. But its exported position
in the API is inconsistent with the pattern that all high-level verbs take a `UCPMap`.
There is no `choi_diff(Φ::UCPMap)` overload. This is a design choice (choi_diff is
a low-level building block), but the docstring should note: "For the high-level path,
use `choi(Φ::UCPMap)` for the Choi of Φ itself; `choi_diff` is for the difference Φ²−Φ
and is called internally by the arb certifiers."

---

**D4. `blocks(v::MainIsomorphism)` — docstring says "C* block sizes d_l of the target algebra B" but doesn't name the return type**

The docstring (types_algebra.jl:124-128) is functionally correct but does not
state the return type (`Vector{Int}`). The accessor for `CStarAlgebra` does:
"The C* block sizes d_l of B = ⊕_l M_{d_l}." The MainIsomorphism version is
shorter. Add the return type annotation.

---

**D5. `width`, `midpoint` — technically correct but no example of when width is large**

Not a strict inconsistency, but the only context in which these accessors are
meaningful is when the bracket is loose (the eig-free path, where `hi/lo ~ 2n`).
Without an example, a user might use `width` to compare bracket quality without
understanding that the ratio `hi/lo` (not `width`) is the normalized looseness
measure. The docstrings could add: "For the solver-free eig-free bracket,
`width` scales as ~2n·η; use `midpoint/width` or compare to `value` for
normalized looseness."

---

## 6. Suggested @docs grouping for the API reference page

Instead of a flat `@autodocs` block (which would dump all 44 exported names in
alphabetical order), organize into the following logical sections. Each section
becomes one `@docs` block in the Documenter page.

### Group 1 — Computing the defect (no solver needed)

```
@docs
certified_defect
eta_eigfree
choi_diff
choi
```

Narrative: start with the solver-free defect tools; the user's first interaction
with the package is always `certified_defect(Φ)`.

---

### Group 2 — The pipeline verbs (the Kitaev construction)

```
@docs
associated_algebra
main_isomorphism
factorize
```

Narrative: the three headline verbs in paper order (§12 → th_almost_idemp → th_main → th_factorization).

---

### Group 3 — Reading a factorization

```
@docs
encode
decode
iscptp
delups_defect
upsdel_defect
algebra
```

Narrative: the accessors a user needs immediately after `factorize(Φ)`.

---

### Group 4 — The input type: UCPMap

```
@docs
UCPMap
n
nkraus
kraus
isunital
```

Note: `Base.adjoint` (as `Φ'`) is a Base extension — include a manual note
rather than a `@docs` entry, since `adjoint` is not in the export list.

---

### Group 5 — Result types

```
@docs
CertifiedBracket
width
midpoint
value
EpsCStarAlgebra
dim
ambient
epsilon
CStarAlgebra
blocks
dim_B
n_B
MainIsomorphism
cbnorm_forward
cbnorm_inverse
isodefect
ChannelFactorization
CPMap
```

Note: `blocks(v::MainIsomorphism)` and `blocks(B::CStarAlgebra)` share a name;
both should appear in the same section.

---

### Group 6 — Solver-gated: exact value (MOSEK extension)

```
@docs
idempotency_defect
eta_idempotence
diamond_norm_watrous
diamond_norm_watrous_primal
diamond_norm_dual
```

This group should have a clear admonition block: "The names in this section
require the AICMosekExt package extension. Without it, calling them throws a
helpful install hint. The solver-free alternative is `certified_defect(Φ)`,
which is always available."

---

### Group 7 — Library discovery

```
@docs
libaic_path
set_libaic_path!
```

This group belongs at the bottom ("installation troubleshooting") and should
cross-reference `BUILDING.md` for the cmake build recipe.

---

### Group 8 — Internal ccall wrappers (NOT in the API reference)

`ccall.jl` / `ccall_factorize.jl` functions (`_assoc_summary`, `_main_iso_summary`,
`_factorize_artifacts`, `_cbnorm_certify`) are intentionally internal and should
NOT appear in the docsite. They are underscore-prefixed and not exported.

---

*End of audit. File path: `docs/research/docsite/brief-api.md`.*
