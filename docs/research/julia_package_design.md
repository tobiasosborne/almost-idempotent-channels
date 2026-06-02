# AlmostIdempotentChannels.jl — master design doc (aic-exa.2 [D])

> Lead designer: Opus. Gates all of aic-exa.3 [C], .4 [J1], .5 [J2], .6 [J3],
> .7 [J4], .8 [J5], .9 [T], .10 [X].
>
> Read order: CLAUDE.md → paper/FINDINGS.md §C → docs/research/julia_package_research_{su2fft,bestpractice}.md → THIS DOC.
>
> Ground truth read first-hand for this doc: the current package
> (`Project.toml`, `src/AlmostIdempotentChannels.jl`, `src/sdp.jl`,
> `test/runtests.jl`, `deps/build.jl`); the headers `aic.h`, `aic_ucp_shim.h`,
> `aic_factorize_shim.h`, `aic_opspace_shim.h`, `aic_factorize.h`, `aic_ecstar.h`,
> `aic_cstar.h`, `aic_idemp.h`, `aic_assoc.h`, `aic_dhom.h`, `aic_ucp.h`; the two
> research legs; FINDINGS §A1–A2, §C1–C21, §D1–D2.

---

## 0. The one-paragraph design

`AlmostIdempotentChannels.jl` exposes the constructive Kitaev pipeline as a
small set of value-types and verbs that read like the paper. The headline verb
chain is

```julia
Φ  = UCPMap(kraus)              # a UCP self-map Φ: B(H)→B(H)
η  = certified_defect(Φ)        # rigorous bracket on ‖Φ²−Φ‖_cb, NO solver (default)
A  = associated_algebra(Φ)      # the ε-C* algebra A = Img Φ̃  (ε = O(η))
v  = main_isomorphism(Φ)        # O(ε)-isomorphism v: B → A, dim-independent constant
F  = factorize(Φ)               # Factorization(B, Δ, Υ) + certified round-trip defects
```

Every defect comes back as a `CertifiedBracket` (`lo ≤ value ≤ hi`, the C arb
balls rounded outward) — that is the rigorous, solver-free spine of the package.
The **exact cb-norm value** (`idempotency_defect`, the Watrous diamond-norm SDP)
needs MOSEK and lives **only** in a package extension. The core package
precompiles, loads, and tests green with **no solver installed**.

The ABI strategy is: **reuse `aic_factorize_choi_shim_d` as the single pipeline
rebuild**, and add a small number of **flat-double shim entry points** (no FLINT
type, no opaque handle) that emit the artifacts the Julia types want — Kraus of
Δ/Υ/Dec/Enc, B's block sizes, the algebra summary, the v / v⁻¹ cb-defect
brackets. Opaque handles are rejected (rationale §4.1): the artifacts are
read-once value objects; a flat-array snapshot is simpler, finalizer-free, and
matches the existing six `_d` shims byte-for-byte.

---

## 1. THE PUBLIC API SURFACE

Naming principle (CLAUDE.md Rule 12, "read like the paper"): **nouns are the
math objects** (`UCPMap`, `EpsCStarAlgebra`, `CStarAlgebra`, `Factorization`,
`MainIsomorphism`, `CertifiedBracket`); **verbs are the constructions**
(`associated_algebra`, `main_isomorphism`, `factorize`, `encode`, `decode`,
`certified_defect`, `idempotency_defect`). Singular nouns (one object), verbs
present-tense. `Φ` is the universal first argument.

### 1.1 Exports

```julia
export UCPMap, EpsCStarAlgebra, CStarAlgebra, Factorization, MainIsomorphism,
       CertifiedBracket
export certified_defect, idempotency_defect,
       associated_algebra, main_isomorphism, factorize,
       encode, decode, choi, kraus, isunital, iscptp
export libaic_path, set_libaic_path!
# kept for back-compat / power users (documented as low-level):
export choi_diff, eta_eigfree, eta_idempotence
```

### 1.2 `UCPMap` — the channel (the universal input)

```julia
"""
    UCPMap(kraus::Vector{<:AbstractMatrix}; check=true, atol=1e-9)
    UCPMap(; choi::AbstractMatrix, n::Int, check=true)

A unital completely-positive (UCP) self-map Φ: B(H) → B(H), dim H = n, in the
HEISENBERG picture  Φ(X) = Σ_a K_a† X K_a  (observables;
aic_ucp.h:8). `kraus` are the {K_a}, each n×n. The dual channel (states) is
Φ* (Schrödinger). Unitality (Σ_a K_a† K_a = I, asserted up to `atol` when
`check`) is the defining UCP-self-map hypothesis the whole pipeline needs.

The Choi constructor reconstructs {K_a} from a Convention-A Choi
(C[i*n+a, j*n+b] = Σ_x conj(K_x[i,a]) K_x[j,b]; conjugation on the FIRST factor,
aic_ucp.h:36) by a tolerance-aware PSD eigen-extraction.
"""
struct UCPMap ... end          # fields §2.1

n(Φ::UCPMap)::Int              # dim H
nkraus(Φ::UCPMap)::Int         # number of Kraus ops r
kraus(Φ::UCPMap)::Vector{Matrix{ComplexF64}}
choi(Φ::UCPMap; prec=106)::Matrix{ComplexF64}   # Convention-A Choi of Φ (not Φ²−Φ)
isunital(Φ::UCPMap; atol=1e-9)::Bool
adjoint(Φ::UCPMap)::UCPMap     # Base.adjoint: the Schrödinger dual Φ* (Kraus K_a†)
```

Example:

```julia
julia> using AlmostIdempotentChannels, LinearAlgebra
julia> Φ = UCPMap([Matrix{ComplexF64}(I,2,2)]);   # identity channel
julia> n(Φ), nkraus(Φ)
(2, 1)
```

Validation: empty kraus → `ArgumentError`; ragged/non-square → `DimensionMismatch`;
non-unital beyond `atol` → `ArgumentError` naming the measured `‖Σ K†K − I‖`.

### 1.3 The defect API — TWO entry points, one default

```julia
"""
    certified_defect(Φ::UCPMap; prec=106) -> CertifiedBracket

A RIGOROUS two-sided bracket on the idempotency defect η = ‖Φ²−Φ‖_cb
(approximate_algebras.tex:347). Solver-free: the eig-free arb certifier
(aic_cbnorm_eigfree_d) rounded outward, so lo ≤ η ≤ hi holds to the arb working
precision. THE DEFAULT — works with no MOSEK. The bracket is loose by design
(hi/lo ~ 2n); it certifies a value rather than computing it. `value` field is
`nothing` (no point estimate without the SDP). See FINDINGS — cb-norm ≠ op-norm.
"""
certified_defect(Φ::UCPMap; prec=106)::CertifiedBracket

"""
    idempotency_defect(Φ::UCPMap; prec=106) -> Float64        # MOSEK extension only

The EXACT value η = ‖Φ²−Φ‖_cb = ‖Φ²−Φ‖_⋄ via the Watrous diamond-norm SDP.
Requires Convex + Mosek + MosekTools (a package extension). Without them this
throws a helpful error telling the user to `using Mosek, MosekTools, Convex`.
"""
idempotency_defect(Φ::UCPMap; prec=106)::Float64
```

Rationale for the split (the prompt's hard constraint): `certified_defect` is the
core, solver-free, rigorous answer and is the package's headline; `idempotency_defect`
is the opt-in exact value. We do NOT make the bare value the default — a value
with no certificate and a mandatory heavyweight solver is the wrong default for a
package whose subject is rigorous bounds.

Example:

```julia
julia> b = certified_defect(UCPMap([Matrix{ComplexF64}(I,2,2)]));
julia> 0.0 in b
true
```

### 1.4 `associated_algebra` → `EpsCStarAlgebra`

```julia
"""
    associated_algebra(Φ::UCPMap; prec=256) -> EpsCStarAlgebra

The associated ε-C* algebra A = Img Φ̃, Φ̃ = θ(2Φ−1) the exact-idempotent
regularization (th_almost_idemp, approximate_algebras.tex:2184-2195; built by
aic_assoc_ecstar_from_phi). A is an OBLIQUE subspace of M_n (FINDINGS §C3); its
product is the Choi–Effros star X⋆Y = Φ̃(XY), NOT plain XY (FINDINGS §C2); the
axioms hold only UP TO ε (no exact unit). Carries dim A, the ambient n, and a
certified ε = associator defect (aic_ecstar_defect_assoc). FAILS LOUD if Φ is
out of the prop_P basin ρ(Φ²−Φ) ≥ 1/4 (aic_assoc.h; the error names ρ).
"""
associated_algebra(Φ::UCPMap; prec=256)::EpsCStarAlgebra

dim(A::EpsCStarAlgebra)::Int        # dim_A
ambient(A::EpsCStarAlgebra)::Int    # n  (A ≤ M_n)
epsilon(A::EpsCStarAlgebra)::CertifiedBracket   # the associator-defect ε bracket
```

We deliberately do NOT expose a `star(A, X, Y)` Julia method in 0.2: the star is
oblique + only-ε-associative and exposing it invites the §C2 "plain XY" bug at
the user boundary. The algebra object is a *summary + handle to the C-built
data*, used internally by `main_isomorphism`/`factorize`. (Revisit in 0.3 with a
guarded `⋆` operator if users ask.)

### 1.5 `main_isomorphism` → `MainIsomorphism`

```julia
"""
    main_isomorphism(Φ::UCPMap; prec=256) -> MainIsomorphism

The O(ε)-isomorphism v: B → A of th_main (approximate_algebras.tex:460), B a
GENUINE finite-dim C* algebra B = ⊕_l M_{d_l}, with a UNIVERSAL, dimension-
INDEPENDENT constant (aic_cstar_build). Carries B's block sizes, and the
certified cb-norm brackets ‖v‖_cb, ‖v⁻¹‖_cb (opspace O2; ≈ 1+O(ε)) plus the
isomorphism defect bracket. The dim-independence of the constant is the headline
guarantee — assert it in tests, never let it grow with n (FINDINGS §D2).
"""
main_isomorphism(Φ::UCPMap; prec=256)::MainIsomorphism

blocks(v::MainIsomorphism)::Vector{Int}   # d_l, the C* block sizes of B
cbnorm_forward(v::MainIsomorphism)::CertifiedBracket    # ‖v‖_cb
cbnorm_inverse(v::MainIsomorphism)::CertifiedBracket    # ‖v⁻¹‖_cb
isodefect(v::MainIsomorphism)::CertifiedBracket
```

### 1.6 `factorize` → `Factorization` (the headline)

```julia
"""
    factorize(Φ::UCPMap; prec=256) -> Factorization

The approximate channel factorization Φ ≈ Δ Υ through a genuine C* algebra B
(th_factorization, approximate_algebras.tex:2730-2899). Δ: B → B(H) (encode,
UCP) and Υ: B(H) → B (decode, UCP) satisfy the EXACT-by-construction targets
ΔΥ = Φ̃, ΥΔ = 1_B, certified by two cb-norm round-trip brackets:
  ‖ΔΥ − Φ‖_cb   (DelUps, .tex:2733)   — encode∘decode returns to Φ within O(η)
  ‖ΥΔ − 1_B‖_cb (UpsDel, .tex:2739)   — decode∘encode returns to the identity on B
Reuses aic_factorize_choi_shim_d, which already rebuilds assoc→cstar→factorize.
The dual channels Dec = Δ*, Enc = Υ* (.tex:2159) are the state-space (CPTP) maps
(`encode`/`decode` below).
"""
factorize(Φ::UCPMap; prec=256)::Factorization

algebra(F::Factorization)::CStarAlgebra      # B = ⊕ M_{d_l}
encode(F::Factorization)::UCPMap             # Dec = Δ* : states on B(H) ← states on B  (CPTP)
decode(F::Factorization)::UCPMap             # Enc = Υ* : states on B(H) → states on B  (CPTP)
delups_defect(F::Factorization)::CertifiedBracket   # ‖ΔΥ − Φ‖_cb bracket
upsdel_defect(F::Factorization)::CertifiedBracket   # ‖ΥΔ − 1_B‖_cb bracket
```

**Naming the dual direction (FINDINGS / CLAUDE.md UCP-vs-CPTP callout).** Δ/Υ are
the Heisenberg (observable) maps; the *channels* are their duals. The paper sets
`Dec = Δ*`, `Enc = Υ*` (.tex:2159, aic_factorize.h:306-316). To avoid the
"transpose everything" trap we name the EXPORTED verbs by their physical role on
states: `encode(F)` returns the channel that takes a B-state into B(H) (= `Dec` =
`Δ*`, dims dim_K=n_B→dim_H=N… wait: aic_factorize_dec_kraus gives dim_K=n_B,
dim_H=N) and `decode(F)` the channel B(H)→B (= `Enc` = `Υ*`). **Implementation
contract for [J4]:** wire `encode → aic_factorize_dec_kraus`, `decode →
aic_factorize_enc_kraus`, and assert the dims in a test so a swap is caught
(round-trip `decode∘encode ≈ 1_B`). This direction binding is load-bearing; the
[C] shim returns both Kraus sets tagged with their dims so Julia cannot guess
wrong.

Example (the headline tour):

```julia
julia> Φ = UCPMap(phit_kraus(2, 0.1));     # near-idempotent
julia> certified_defect(Φ)                 # rigorous η bracket, no solver
CertifiedBracket: 0.084 ≤ ‖Φ²−Φ‖_cb ≤ 0.27
julia> F = factorize(Φ);
julia> algebra(F)
CStarAlgebra: ⊕ M_d, blocks = [1, 1]  (dim_B = 2, n_B = 2)
julia> delups_defect(F)                    # encode∘decode ≈ Φ
CertifiedBracket: 0.0 ≤ ‖ΔΥ − Φ‖_cb ≤ 0.21
```

### 1.7 `CertifiedBracket` — the rigorous-result type

```julia
"""
    CertifiedBracket{T<:AbstractFloat}

A RIGOROUS two-sided bracket lo ≤ x ≤ hi on a scalar (the arb balls rounded
outward, FLOOR/CEIL — aic_ucp_shim.h:46-52). Optional `value` (a point estimate,
e.g. the MOSEK SDP value) when one is available; `nothing` for solver-free
brackets. `label` names the quantity for show ("‖Φ²−Φ‖_cb").
"""
struct CertifiedBracket{T} ... end          # fields §2.5

Base.in(x::Real, b::CertifiedBracket)        # lo ≤ x ≤ hi
Base.minimum(b)/Base.maximum(b)              # lo / hi
width(b::CertifiedBracket)                   # hi − lo
midpoint(b::CertifiedBracket)                # (lo+hi)/2
value(b::CertifiedBracket)                   # the point estimate or `nothing`
```

---

## 2. THE TYPE SYSTEM

All value types are **immutable** `struct`s (no `mutable struct`, no finalizers —
§4.1 rejects opaque handles). `ComplexF64`/`Float64` storage; matrices are dense
`Matrix{ComplexF64}` in the Convention-A / K-major flat layouts the shims already
use. Parametrize only `CertifiedBracket` (on the float type) — the rest are
concrete.

### 2.1 `UCPMap`

```julia
struct UCPMap
    kraus::Vector{Matrix{ComplexF64}}   # r operators, each n×n (K_a: H→H, Heisenberg)
    n::Int                              # dim H
    r::Int                              # number of Kraus ops
end
```

Constructors normalize ragged input to `Vector{Matrix{ComplexF64}}`, validate,
and (when `check`) assert unitality. `adjoint(Φ)` returns `UCPMap([K' for K in Φ.kraus])`.

### 2.2 `EpsCStarAlgebra`

```julia
struct EpsCStarAlgebra
    n::Int                       # ambient dim, A ≤ M_n
    dim_A::Int                   # dim A
    basis::Vector{Matrix{ComplexF64}}  # Frobenius-orthonormal {B_k}, n×n (from [C] shim)
    eps::CertifiedBracket{Float64}     # associator-defect ε bracket (label "ε_assoc")
    source::UCPMap               # the Φ it came from (for re-derivation / show)
end
```

The `basis` is a snapshot of A's `B_k` (aic_ecstar.B). We store it so the object
is self-contained and showable, but the algebra's *operations* are not re-done in
Julia — `main_isomorphism`/`factorize` rebuild from `source` through the C
pipeline (which is cheap and avoids re-marshalling the oblique star).

### 2.3 `CStarAlgebra`

```julia
struct CStarAlgebra
    blocks::Vector{Int}   # d_l, the genuine C* block sizes (B = ⊕_l M_{d_l})
    dim_B::Int            # Σ d_l²
    n_B::Int              # Σ d_l   (block-diagonal representative size)
end
```

Pure structural summary of `aic_dhom_B` — no matrices needed.

### 2.4 `MainIsomorphism` and `Factorization`

```julia
struct MainIsomorphism
    source::UCPMap
    B::CStarAlgebra
    cbnorm_fwd::CertifiedBracket{Float64}   # ‖v‖_cb
    cbnorm_inv::CertifiedBracket{Float64}   # ‖v⁻¹‖_cb
    isodefect::CertifiedBracket{Float64}
end

struct Factorization
    source::UCPMap
    B::CStarAlgebra
    Δ::Vector{Matrix{ComplexF64}}   # encode Kraus  (Dec = Δ*: dim_K=n_B, dim_H=N)  aic_factorize.h:306
    Υ::Vector{Matrix{ComplexF64}}   # decode Kraus  (Enc = Υ*: dim_K=N,   dim_H=n_B)
    delups::CertifiedBracket{Float64}   # ‖ΔΥ − Φ‖_cb
    upsdel::CertifiedBracket{Float64}   # ‖ΥΔ − 1_B‖_cb
    eta_proxy::Float64                  # out_eta[0] from the shim (sanity, ≈0 at η=0)
    eps_used::Float64                   # out_eta[1] (the ε actually fed to cstar_build)
end
```

`encode(F)`/`decode(F)` wrap `F.Δ`/`F.Υ` (with the dual dims) into `UCPMap`-like
channels — but note these are CPTP state maps, not unital self-maps, so they are
returned as a `UCPMap` with `check=false` (unitality is the dual's trace-preservation
and is checked by the round-trip test, not the constructor). **[J4] decision:**
return them as `UCPMap` (Kraus container) and document that for Dec/Enc the
"unital" predicate is replaced by trace-preservation; the round-trip defect is
the real certificate.

### 2.5 `CertifiedBracket`

```julia
struct CertifiedBracket{T<:AbstractFloat}
    lo::T
    hi::T
    value::Union{T,Nothing}     # point estimate (MOSEK) or nothing (solver-free)
    label::String               # "‖Φ²−Φ‖_cb", "ε_assoc", "‖ΔΥ−Φ‖_cb", ...
    function CertifiedBracket(lo::T, hi::T; value=nothing, label="") where {T<:AbstractFloat}
        lo ≤ hi || throw(ArgumentError("CertifiedBracket: lo=$lo > hi=$hi"))
        value === nothing || lo - 1e-9 ≤ value ≤ hi + 1e-9 ||
            throw(ArgumentError("CertifiedBracket: value=$value outside [$lo,$hi]"))
        new{T}(lo, hi, value, label)
    end
end
```

We do **not** depend on IntervalArithmetic.jl (§7) — this 6-line type is enough,
and an IEEE-1788 interval would imply rounded arithmetic semantics we do not
provide (the C side already did the rigorous rounding).

### 2.6 `Base.show` renderings

Compact `show` (one line, used in arrays / interpolation) and multi-line
`show(::MIME"text/plain")` (the REPL display). Sketches of actual output:

`UCPMap`, compact: `UCPMap(n=2, r=4)`. Multi-line:

```
UCPMap: Φ: B(H) → B(H),  dim H = 2
  Kraus operators: 4   (Heisenberg: Φ(X) = Σ K_a† X K_a)
  unital: yes   (‖Σ K†K − I‖ = 3.1e-16)
```

`CertifiedBracket`, compact: `[0.084, 0.27]` (or with value `[0.084, 0.135, 0.27]`).
Multi-line (label drives the line):

```
CertifiedBracket  ‖Φ²−Φ‖_cb
  0.0840  ≤  η  ≤  0.2700        (width 0.186)
```

with a `value =` line appended when present:

```
  value = 0.1350   (Watrous SDP)
```

`EpsCStarAlgebra`, multi-line:

```
EpsCStarAlgebra  A = Img Φ̃  ≤  M_2
  dim A = 2
  ε (associator defect) ∈ [2.1e-12, 4.4e-12]
  product: Choi–Effros star  X⋆Y = Φ̃(XY)   (oblique; axioms hold up to ε)
```

`CStarAlgebra`, compact `⊕ M_d [1,1]`; multi-line:

```
CStarAlgebra  B = ⊕_l M_{d_l}
  blocks d_l = [1, 1]   (m = 2)
  dim_B = 2,  n_B = 2
```

`MainIsomorphism`, multi-line:

```
MainIsomorphism  v: B → A   (O(ε), dimension-independent constant)
  B = ⊕ M_d, blocks = [1, 1]
  ‖v‖_cb    ∈ [1.00, 1.04]
  ‖v⁻¹‖_cb  ∈ [1.00, 1.05]
  iso defect ∈ [0.0, 0.02]
```

`Factorization`, multi-line (the showcase):

```
Factorization  Φ ≈ Δ Υ   through  B = ⊕ M_d, blocks = [1, 1]
  encode  Δ: B → B(H)   (12 Kraus,  Dec = Δ*)
  decode  Υ: B(H) → B   ( 8 Kraus,  Enc = Υ*)
  ‖ΔΥ − Φ‖_cb   ∈ [0.0,  0.21]     (encode∘decode ≈ Φ)
  ‖ΥΔ − 1_B‖_cb ∈ [0.0,  0.18]     (decode∘encode ≈ 1_B)
  η proxy = 0.031,  ε used = 0.031
```

`show` must use only ASCII-safe-on-failure Unicode that the project already uses
in comments (Φ, Υ, Δ, ⊕, ≤, ∈) — these render in the project's existing docs.

---

## 3. THE EXTENSION ARCHITECTURE (MOSEK optional)

### 3.1 Project.toml block (exact)

```toml
name = "AlmostIdempotentChannels"
uuid = "6bf14b85-acf1-40b3-a7a0-aef96f56c734"
authors = ["Tobias Osborne <tobias.j.osborne@gmail.com>"]
version = "0.2.0"

[deps]
Libdl         = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Preferences   = "21216c6a-2e73-6563-6e65-726566657250"

[weakdeps]
Convex     = "f65535da-76fb-5f13-bab9-19810c17039a"
Mosek      = "6405355b-0ac2-5fba-af84-adbd65488c0e"
MosekTools = "1ec41992-ff65-5c91-ac43-2df89e9693a4"

[extensions]
AICMosekExt = ["Convex", "Mosek", "MosekTools"]

[compat]
julia       = "1.10"
Preferences = "1"
Convex      = "0.16"
Mosek       = "11"
MosekTools  = "0.15, 0.16"

[extras]
Test   = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Aqua   = "4c88cf16-eb10-579e-8560-4a9242c79595"
# the weakdeps are ALSO listed in extras so Pkg.test installs them and the
# extension activates during the (optional) SDP test pass:
Convex     = "f65535da-76fb-5f13-bab9-19810c17039a"
Mosek      = "6405355b-0ac2-5fba-af84-adbd65488c0e"
MosekTools = "1ec41992-ff65-5c91-ac43-2df89e9693a4"

[targets]
test = ["Test", "Random", "Aqua", "Convex", "Mosek", "MosekTools"]
```

`julia = "1.10"` (not 1.9): the project already pins 1.10 and Arblib uses 1.10;
parallel GC is irrelevant to a serial-Julia project but 1.10 is the safe floor.
JET is intentionally NOT a hard test dep in 0.2 (run opportunistically; the
research leg flags it as flaky as a gate — §8.2).

### 3.2 ext/ file layout

```
ext/AICMosekExt.jl     # the whole Watrous SDP (current src/sdp.jl moves here verbatim)
```

The extension extends stubs declared in the core (next paragraph). It owns the
five `diamond_*` functions AND the solver-gated public verbs
`idempotency_defect` and the tight-bracket path.

### 3.3 Dispatch + graceful degradation

Core declares the solver-gated symbols as **stubs that throw a helpful error**:

```julia
# src/sdp_stubs.jl  (loaded by the core)
const _MOSEK_HINT = """
    requires the MOSEK extension. Install + load the solver:
        using Pkg; Pkg.add(["Convex","Mosek","MosekTools"])
        using Convex, Mosek, MosekTools
    (Mosek needs a license: https://www.mosek.com/products/academic-licenses/)
"""
idempotency_defect(Φ::UCPMap; kw...) = error("idempotency_defect " * _MOSEK_HINT)
diamond_norm_watrous(::Matrix{ComplexF64}, ::Int) = error("diamond_norm_watrous " * _MOSEK_HINT)
# ... the other four diamond_* stubs ...
```

The extension defines `AlmostIdempotentChannels.idempotency_defect(Φ; prec)` etc.,
which shadow the stubs once `Convex+Mosek+MosekTools` are all loaded. This is the
register-by-method-override pattern (research §1.3); when the extension is loaded
the more specific method (still the same signature, but defined later in the
extension module) wins via Julia's normal world-age method table — for a clean
override we instead have the stub call an internal `_idempotency_defect_impl(Φ,
prec)` that the extension *adds a method to* (the extension owns the only real
method; the stub's fallback fires only when none exists). Concretely:

```julia
# core:
idempotency_defect(Φ::UCPMap; prec=106) = _diamond_value_impl(Φ, prec)
_diamond_value_impl(::UCPMap, ::Int) = error("idempotency_defect " * _MOSEK_HINT)
# ext:
AlmostIdempotentChannels._diamond_value_impl(Φ::UCPMap, prec::Int) = ... # real SDP
```

This keeps the public name + docstring in the core (so it is discoverable and
Documenter-visible), with the implementation gated. Also register a
`Base.Experimental.register_error_hint` in the core `__init__` so a `MethodError`
on the `diamond_*` stubs prints the same install hint (belt and suspenders).

### 3.4 What is core vs extension

| Capability | Core (no solver) | Extension (MOSEK) |
|---|---|---|
| `certified_defect` (rigorous bracket) | ✅ default | — |
| `associated_algebra`, `main_isomorphism`, `factorize` | ✅ (full C pipeline) | — |
| `delups_defect`/`upsdel_defect`/`isodefect` brackets | ✅ (eig-free arb) | — |
| `idempotency_defect` (exact η value) | ❌ throws hint | ✅ Watrous SDP |
| tight `certified_defect(...; tight=true)` | ❌ throws hint | ✅ aic_cbnorm_certify_d + MOSEK feasible points |
| `diamond_norm_watrous*`, `diamond_*_rect` | ❌ throws hint | ✅ |

The entire headline pipeline is solver-free. Only the cb-norm *value* and the
*tight* bracket need MOSEK.

---

## 4. THE ABI STRATEGY (gates aic-exa.3 [C])

### 4.1 Flat-double arrays, NOT opaque handles — for every artifact

Decision: **all new artifacts cross the boundary as flat `double*` re/im arrays +
`int` scalars**, identical in spirit to the six existing `_d` shims. No opaque
`Ptr{Cvoid}` handle + finalizer.

Rationale:
- The artifacts are **read-once value snapshots** (Kraus lists, block sizes,
  bracket endpoints). Julia copies them into immutable structs and never calls
  back into C with them. There is no long-lived C object to manage, so a handle
  buys nothing and costs a finalizer + thread-safety reasoning (research §4.3).
- It matches the **existing contract** exactly (`aic_*_shim_d`: int dims +
  `Ptr{Cdouble}` re/im, caller-owns, FLINT types never cross — aic_ucp_shim.h,
  aic_factorize_shim.h). One marshalling idiom for the whole package.
- Sizes are bounded and known up front (n, r, n_B all ≤ n; safe upper bound
  `n^4` for Choi-sized buffers as the factorize shim already documents), so the
  caller can pre-allocate. **The one place this is awkward** — the number of
  Dec/Enc Kraus operators is data-dependent — is handled by a **two-call
  protocol** (query sizes, then fill), §4.3.
- The serial-Julia constraint and the absence of any mutate-in-place C object
  remove the only motivations for a handle.

(An opaque handle would only pay off if we later expose *interactive* algebra
operations — applying the star, composing maps — in C with state kept alive
across calls. That is out of 0.2 scope; revisit if 0.3 adds an interactive
algebra surface.)

### 4.2 Reuse `aic_factorize_choi_shim_d`

The existing `aic_factorize_choi_shim_d` already rebuilds the WHOLE pipeline
(`kraus_from_flat → aic_assoc_ecstar_from_phi → aic_cstar_build →
aic_factorize_build/_delta_build/_upsilon_build`) and returns only the two SDP-
feeding Choi. The new shims **call the same rebuild** and emit the artifacts the
Julia types need instead of (or in addition to) the Choi. Each new shim is "the
same pipeline build + a different extraction tail," keeping them ≤200 LOC (Rule
10) and sharing the in-basin/eps-sentinel behavior already documented.

### 4.3 The new C shim entry points (the precise contract for aic-exa.3)

All follow the house contract: `int` return (0 = OK; fail-loud `assert` on bad
shape per Rule 4), caller owns all buffers, FLINT types never cross, Convention-A
/ K-major layouts preserved, cite `.tex` + `FINDINGS §Cn`.

**(C1) `aic_cbnorm_eta_bracket_d`** — certified idempotency bracket (rename/alias
of the existing `aic_cbnorm_eigfree_d`; already present). NO new code needed;
`certified_defect` wraps the existing symbol. Listed here for completeness.

```c
int aic_cbnorm_eigfree_d(double *lo, double *hi,
                         const double *kraus_re, const double *kraus_im,
                         int n, int r, int prec);   /* EXISTS — reuse */
```

**(C2) `aic_assoc_summary_d`** — associated-algebra summary + basis snapshot.

```c
/* Build A = Img Φ̃ from Φ's flat Kraus (aic_assoc_ecstar_from_phi), report:
 *   *out_n     = n (ambient), *out_dimA = dim A;
 *   eps_assoc[0] = lo, eps_assoc[1] = hi  : certified ε = aic_ecstar_defect_assoc
 *                bracket (arf lower/upper rounded FLOOR/CEIL, as the eigfree shim).
 *   basis_re/basis_im : the dim_A Frobenius-orthonormal B_k, each n×n, laid out
 *                [k*n*n + i*n + j]  (K-major like the Kraus layout). Caller pre-
 *                allocates length n^4 (dim_A ≤ n²; safe upper bound (n²)*(n*n)=n^4).
 * Wraps: aic_assoc_ecstar_from_phi + aic_ecstar_defect_assoc. FAILS LOUD out of
 * the prop_P basin (the assoc_regularize assert) — the caller pre-checks via the
 * eta bracket and surfaces a helpful Julia error. Cites .tex:2184, FINDINGS §C2/§C3.
 * Returns 0; asserts n≥1, r≥1, prec≥2, non-null. */
int aic_assoc_summary_d(int *out_n, int *out_dimA,
                        double *eps_assoc,          /* length 2 */
                        double *basis_re, double *basis_im,  /* length n^4 each */
                        const double *kraus_re, const double *kraus_im,
                        int n, int r, double eps, int prec);
```

**(C3) `aic_main_iso_summary_d`** — th_main isomorphism summary.

```c
/* Rebuild v: B → A (aic_assoc_ecstar_from_phi → aic_cstar_build) and the opspace
 * cb-norm brackets, report:
 *   *out_nB, *out_dimB ; blocks[0..m-1] = B's d_l, *out_m = num_blocks (caller
 *      pre-allocates blocks length n, since m ≤ n_B ≤ n);
 *   cbfwd[2], cbinv[2], isodef[2]  : certified [lo,hi] brackets for ‖v‖_cb,
 *      ‖v⁻¹‖_cb, and the iso defect (eig-free arb, FLOOR/CEIL rounding).
 *   out_eps[0]=eps_used (the eps sentinel, as the factorize shim).
 * Wraps: aic_cstar_build (iso_def) + the opspace O2 eig-free cb brackets
 * (structurally ‖v‖_cb / ‖v⁻¹‖_cb; aic_opspace_* eig-free path, NOT the SDP).
 * NOTE — the opspace SDP VALUE needs MOSEK; here we emit only the eig-free
 * BRACKETS so the core is solver-free. FINDINGS §C12 (use the OPERATOR-norm
 * inclusion inf, not the Frobenius σ_min ampliation). Cites .tex:460, :484, §D2.
 * Returns 0; asserts shapes, non-null. eps sentinel identical to the factorize shim. */
int aic_main_iso_summary_d(int *out_nB, int *out_dimB, int *out_m,
                           int *blocks,                 /* length n */
                           double *cbfwd, double *cbinv, double *isodef,  /* len 2 each */
                           double *out_eps,             /* length 1 */
                           const double *kraus_re, const double *kraus_im,
                           int n, int r, double eps, int prec);
```

**(C4) `aic_factorize_artifacts_sizes_d`** — query Dec/Enc Kraus counts + B shape
(the first half of the two-call protocol; data-dependent sizes).

```c
/* Rebuild the factorize pipeline (reusing the aic_factorize_choi_shim_d build
 * path) and report ONLY the SIZES needed to allocate the artifact buffers:
 *   *out_N, *out_nB, *out_dimB, *out_m, blocks[0..m-1],
 *   *out_rDec, *out_rEnc  : the number of Dec (= Δ*) and Enc (= Υ*) Kraus ops
 *      (aic_factorize_dec_kraus / _enc_kraus return aic_ucp_kraus with a data-
 *      dependent r; the caller needs these to size the Kraus buffers in C5).
 * No heavy SDP. Cites .tex:2730, :2159. Returns 0; asserts shapes, non-null.
 * eps sentinel identical to aic_factorize_choi_shim_d. */
int aic_factorize_artifacts_sizes_d(int *out_N, int *out_nB, int *out_dimB,
                                    int *out_m, int *blocks,    /* blocks length n */
                                    int *out_rDec, int *out_rEnc,
                                    const double *kraus_re, const double *kraus_im,
                                    int n, int r, double eps, int prec);
```

**(C5) `aic_factorize_artifacts_d`** — the full factorization artifacts (the
second half; caller pre-allocates from C4's sizes).

```c
/* Rebuild the factorize pipeline and emit the full artifacts:
 *   Dec (= Δ*) Kraus: dec_re/dec_im, rDec ops, each dim_K=n_B × dim_H=N, laid out
 *      [a*(nB*N) + i*N + j]  (K-major; dim_K major rows). (aic_factorize.h:306-316.)
 *   Enc (= Υ*) Kraus: enc_re/enc_im, rEnc ops, each dim_K=N × dim_H=n_B, laid out
 *      [a*(N*nB) + i*nB + j].
 *   delups[2], upsdel[2] : certified eig-free brackets on ‖ΔΥ−Φ‖_cb (.tex:2733)
 *      and ‖ΥΔ−1_B‖_cb (.tex:2739) — aic_factorize_eigfree_delups/_upsdel, arf
 *      lower/upper FLOOR/CEIL. (The TIGHT cb VALUE is the MOSEK extension's job,
 *      via aic_factorize_choi_shim_d + the SDP — kept solver-free here.)
 *   out_eta[2] : [eta_proxy, eps_used] (identical to aic_factorize_choi_shim_d).
 * The caller MUST have called C4 first and sized dec_*/enc_* to rDec*(nB*N) and
 * rEnc*(N*nB) respectively (the shim re-asserts the counts match what it rebuilds
 * — fail loud if the caller mis-sized). Wraps the SAME build as the choi shim plus
 * aic_factorize_dec_kraus / _enc_kraus / _eigfree_*. FINDINGS §C13 (F-ancilla
 * ordering — Dec/Enc dims are load-bearing, the round-trip test catches a swap),
 * §C14 (Δ′ only O(η)-CP; the Choi→Kraus PSD-cone projection is internal). Cites
 * .tex:2730-2899, :2159. Returns 0; asserts shapes, non-null, the rebuilt rDec/rEnc
 * match the caller's. eps sentinel identical to the choi shim. */
int aic_factorize_artifacts_d(double *dec_re, double *dec_im, int rDec,
                              double *enc_re, double *enc_im, int rEnc,
                              double *delups, double *upsdel,  /* len 2 each */
                              double *out_eta,                  /* len 2 */
                              const double *kraus_re, const double *kraus_im,
                              int n, int r, double eps, int prec);
```

**Memory-ownership contract (all C2–C5).** Caller allocates every buffer (Julia
`Vector{Float64}`); the shim only writes into them; nothing is heap-returned, so
nothing is freed across the boundary. Output `int*` scalars are caller-owned. The
internal `aic_ucp_kraus` Dec/Enc are built and `aic_ucp_kraus_clear`'d inside the
shim (never returned). This is exactly the existing `_d`-shim ownership model.

**The eps sentinel** (reused verbatim from aic_factorize_shim.h:78-89): `eps==0.0`
→ η=0 oracle path; `eps==-2.0` → multi-block η>0 (eps_used = eta_proxy, the §C11
caveat); `eps<0` (≠−2) → derive `aic_ecstar_defect_assoc`; `eps>0` → verbatim.
Julia always passes the matching sentinel and reads back `eps_used` to keep the
certifier reproducible.

### 4.4 C-test obligations for aic-exa.3 (mutation-proof, Rule 7)

- **η=0 oracle:** identity / a genuine exact-idempotent Φ → every defect bracket
  `[lo,hi]` straddles 0 with `hi` at machine-ε, B-block sizes match
  `aic_idemp_wedderburn` (Rule 6 canonical cross-check).
- **η>0 oblique fixture** (`compress_idemp(4,2)` or a `dep+conj(dep)` mix, §C2):
  every star-touching path is exercised; mutate the Dec/Enc dim swap → round-trip
  `decode∘encode` defect jumps O(1) → RED (catches §C13 ordering).
- **Bracket containment:** the new eig-free `delups`/`upsdel` brackets contain the
  MOSEK tight value computed by the existing `aic_factorize_choi_shim_d` + SDP
  fixture (ties the solver-free bracket to the value).
- **Dim-independence canary** (FINDINGS §D2): the iso defect / cb brackets do not
  grow with n across n=2,3,4 fixtures.

---

## 5. FILE / MODULE LAYOUT of the Julia package

```
julia/AlmostIdempotentChannels.jl/
  Project.toml                    # §3.1 (weakdeps + extensions + Preferences)
  src/
    AlmostIdempotentChannels.jl   # module: includes below; __init__ (dlopen+RTLD); Preferences libpath
    libaic.jl                     # LIBAIC_PATH discovery + set_libaic_path! + Ref{Ptr} handles + dlsym
    ccall.jl                      # the @ccall/ccall wrappers over C1–C5 + the 6 existing _d shims; _flatten_kraus, _unflatten_*
    types.jl                      # UCPMap, EpsCStarAlgebra, CStarAlgebra, MainIsomorphism, Factorization, CertifiedBracket
    show.jl                       # Base.show (compact + MIME"text/plain") for every type
    api.jl                        # certified_defect, associated_algebra, main_isomorphism, factorize, encode/decode, choi/kraus
    sdp_stubs.jl                  # idempotency_defect + diamond_* stubs (helpful-error fallbacks)
  ext/
    AICMosekExt.jl                # the Watrous SDP (current src/sdp.jl verbatim) + _diamond_value_impl + tight bracket
  test/
    runtests.jl                   # entry: core testset always; SDP testset try/loaded-or-skip
    test_core.jl                  # T0–T4, T6–T8 (no solver)
    test_factorize.jl             # T7 factorize round-trip (η=0 oracle + η>0 oblique)
    test_sdp.jl                   # T5 strong duality (MOSEK; skipped if absent)
    test_aqua.jl                  # Aqua.test_all
  deps/
    build.jl                      # KEPT for dev: cmake --build (with make fallback), writes nothing if Preferences set
  docs/                           # aic-exa.10 (Documenter): make.jl, src/index.md, src/api.md, src/tutorial.md
  DESIGN.md                       # the marshalling/layout contract (this doc's §2/§4 distilled)
```

`src/sdp.jl` is **deleted** (moves to `ext/`); `deps/deps.jl` is **deleted**
(Preferences replaces it). Each `src/*.jl` stays well under the 200-LOC house
convention by splitting on concern (types vs show vs api vs ccall).

---

## 6. BUILD / DEV STORY

### 6.1 libaic discovery (Preferences + graceful make/cmake fallback)

`src/libaic.jl`:

```julia
using Preferences
const _DEFAULT_LIBAIC = let
    repo = abspath(joinpath(@__DIR__, "..", "..", ".."))   # pkg is 3 levels deep
    joinpath(repo, "build", "libaic.so")
end
const LIBAIC_PATH = @load_preference("libaic_path", _DEFAULT_LIBAIC)

set_libaic_path!(p) = (@set_preferences!("libaic_path" => string(p));
    @info "libaic_path set; restart Julia to apply")
libaic_path() = LIBAIC_PATH
```

`__init__` (keep the load-bearing libgmp preload + RTLD flags verbatim from the
current module — research §1c, su2-fft bead su2fft-e5z):

```julia
function __init__()
    if !isfile(LIBAIC_PATH)
        error("libaic.so not found at \"$LIBAIC_PATH\".\n" *
              "Build it:  cmake -S . -B build && cmake --build build   (in the repo root)\n" *
              "or set the path:  AlmostIdempotentChannels.set_libaic_path!(\"/abs/path/libaic.so\")")
    end
    try Libdl.dlopen("/lib/x86_64-linux-gnu/libgmp.so.10", Libdl.RTLD_NOW|Libdl.RTLD_GLOBAL)
    catch; try Libdl.dlopen("libgmp.so.10", Libdl.RTLD_NOW|Libdl.RTLD_GLOBAL); catch; end; end
    flags = Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL | Libdl.RTLD_DEEPBIND
    _LIBAIC_HANDLE[] = Libdl.dlopen(LIBAIC_PATH, flags)
    for (sym, ref) in _ALL_SYMS           # one Ref{Ptr{Cvoid}} per C entry point
        ref[] = Libdl.dlsym(_LIBAIC_HANDLE[], sym)
    end
    Base.Experimental.register_error_hint(_mosek_hint, MethodError)
end
```

ccall keeps the **dlsym-handle 4-arg form** for the existing six symbols and the
new C2–C5 (research §1e, §2g: `@ccall` cannot take a `Ptr{Cvoid}`; the dlsym
handle is required to combine the RTLD flags with type-safe ccall). Wrap every
buffer pass in `GC.@preserve` (research §4.2) for documentation + future-proofing.

`deps/build.jl` is kept for the dev convenience path but made non-authoritative:
it tries `cmake --build build`, falls back to `make lib`, and **only writes a
Preferences entry if the default path does not resolve** — Preferences is the
source of truth, `deps.jl` is gone. A user who has run cmake once needs no
`Pkg.build` at all (the default path resolves).

### 6.2 Yggdrasil JLL (aic-95g.2) — a documented LATER milestone, NOT a blocker

The end state is `AIC_jll` (BinaryBuilder/Yggdrasil): then `[deps]` gains
`AIC_jll`, `__init__` drops the dlopen+RTLD workaround (the JLL ships a matching
libgmp), and ccall becomes `@ccall libaic.fn(...)::Ret`. This is **out of 0.2
scope**; documented in `DESIGN.md` / `PACKAGING.md` as the migration path
(research §3.1, §2a). Do not block 0.2 on it.

---

## 7. CERTIFIED-RESULTS ERGONOMICS

`CertifiedBracket` (§2.5) is the rigorous-result type. Decision points:

- **`value::Union{T,Nothing}`** field: a bracket is solver-free by default
  (`value === nothing`); the MOSEK extension fills `value` with the SDP point
  estimate. `show` renders the `value` line only when present. This lets the same
  type carry both "rigorous bracket only" (core) and "bracket + exact value"
  (extension) without two types.
- **No IntervalArithmetic.jl dependency.** Recommended AGAINST (research §6.1):
  (i) it is a large dependency tree for a 6-line need; (ii) IEEE-1788 intervals
  imply *interval arithmetic semantics* (rounded ops, set operations) we do not
  want — the rigor was already established in the C/arb layer, and re-wrapping in
  an arithmetic interval would invite users to do `b1 + b2` expecting rigorous
  propagation we are not providing. A plain immutable bracket with `in`,
  `width`, `midpoint`, `value` is the honest surface. If a user wants interval
  arithmetic they can `import IntervalArithmetic; interval(b.lo, b.hi)`.
- `Base.in(x, b)` is the primary assertion idiom in tests (`@test η in b`).
- The bracket is constructed ONLY from the C arb FLOOR/CEIL output — Julia never
  recomputes endpoints (no double-rounding that could break rigor).

---

## 8. VERSION / RELEASE PLAN (semver)

- **0.2.0 — this epic.** The first coherent release: weakdeps/extension MOSEK,
  Preferences libpath, the type system + show, the C2–C5 ABI, the full headline
  API (`certified_defect`, `associated_algebra`, `main_isomorphism`, `factorize`,
  `encode`/`decode`), core-green-without-solver, Aqua gate, Documenter site +
  doctests. Breaking vs 0.1 (return types changed: `eta_eigfree` →
  `CertifiedBracket`; MOSEK no longer a hard dep) — bump minor (pre-1.0 minor is
  the breaking axis). Each child bead lands an internally-coherent slice (the
  package stays loadable + the core suite green at every step; §9).
- **0.3.0 (later).** Optional: a guarded `⋆` star operator on `EpsCStarAlgebra`;
  an interactive algebra surface (would justify opaque handles, §4.1); tensor-
  extension API (§10 of the paper).
- **0.4 / 1.0 (later).** `AIC_jll` via Yggdrasil (drops dlopen workaround);
  registration in General; the Python consumer bindings (out of this project's
  Julia/C-only scope, per the no-Python directive — consumers come later).

`0.1.0` was never registered; 0.2.0 is the first intended-public version.

---

## 9. ORDERED IMPLEMENTATION CHECKLIST (bead-by-bead contracts)

Each step keeps the package loadable and the core suite green. Serial Julia
throughout (one process); never run a [C] cmake build while a Julia agent holds
`build/`.

### aic-exa.3 [C] — extend the ABI (no Julia; cmake+ctest only)
1. Add C2 `aic_assoc_summary_d`, C3 `aic_main_iso_summary_d`, C4
   `aic_factorize_artifacts_sizes_d`, C5 `aic_factorize_artifacts_d` (§4.3). Each
   reuses the existing pipeline build (§4.2); each ≤200 LOC; each cites `.tex` +
   `FINDINGS §Cn`. Put them in new `src/aic_*_artifacts_shim.c` files, declared in
   matching `include/aic/aic_*_artifacts_shim.h`.
2. Header docstrings: layout (K-major Kraus, Convention-A Choi), eps sentinel,
   FLOOR/CEIL rigor, ownership (§4.3).
3. C tests (§4.4): η=0 oracle (machine-ε defects, B-blocks vs Wedderburn), η>0
   oblique fixture (Dec/Enc dim-swap mutation → RED), bracket-contains-SDP-value,
   dim-independence canary. CTest TIMEOUT + fast/slow label.
4. Hostile review MANDATORY (>30 LOC, new ABI).
   **Output contract for downstream:** the five symbol names + exact signatures in
   §4.3 are frozen; .6 [J3] depends on them verbatim.

### aic-exa.4 [J1] — Project.toml + Preferences + ext skeleton
1. Rewrite `Project.toml` to §3.1 (weakdeps/extensions, Preferences, julia 1.10).
2. Add `src/libaic.jl` (§6.1): Preferences libpath, `set_libaic_path!`, the
   `Ref{Ptr{Cvoid}}` handle table, `__init__` (keep RTLD/libgmp verbatim).
3. Add `ext/AICMosekExt.jl` skeleton (empty module that imports the package; real
   methods land in .8). Add `src/sdp_stubs.jl` with the helpful-error stubs (§3.3).
4. Delete `src/sdp.jl` content into `ext/` is .8's job; here just move it so the
   core does not `using Convex/Mosek`. Delete `deps/deps.jl`.
5. **Gate:** `julia --project -e 'using AlmostIdempotentChannels'` precompiles +
   loads with NO solver present (one process).

### aic-exa.5 [J2] — type system + Base.show
1. `src/types.jl`: the six structs (§2), immutable, validated constructors
   (ArgumentError/DimensionMismatch with measured values, research §1h).
2. `src/show.jl`: compact + multi-line `show` for each (§2.6 sketches are the
   spec — match them).
3. Pure-Julia; no ccall yet (constructors that need C data take the raw fields;
   `UCPMap`'s Kraus/Choi constructors can be pure-Julia — Choi→Kraus via a
   Julia PSD eig is fine for validation, the C path is for the pipeline).
4. **Gate:** serial julia smoke — construct + `show` each type.

### aic-exa.6 [J3] — @ccall wrapper layer over the new ABI
1. `src/ccall.jl`: one wrapper per C1–C5 + reuse the existing
   `choi_diff`/`eta_eigfree` wrappers; `_flatten_kraus` (exists) +
   `_unflatten_kraus` / `_unflatten_choi` helpers; two-call protocol for C4→C5.
2. Dlsym-handle 4-arg `ccall`, `GC.@preserve`, fail-loud on rc≠0 (research §1e,h).
3. **Gate:** round-trip each wrapper against the C test values (the η=0 oracle
   numbers from the [C] tests) in a serial julia run. Depends on [C] landed; not
   while a [C] build runs.

### aic-exa.7 [J4] — high-level idiomatic API (the joy layer)
1. `src/api.jl`: `certified_defect` (wraps eigfree), `associated_algebra` (C2),
   `main_isomorphism` (C3), `factorize` (C4→C5), `encode`/`decode` (wrap Δ/Υ
   Kraus with the DUAL dims — §1.6 direction binding; assert dims), `choi`,
   `kraus`, `isunital`, `iscptp`.
2. Helpful errors: out-of-basin Φ → catch the C fail-loud, re-raise naming
   ρ(Φ²−Φ) and pointing at `certified_defect`; non-unital → ArgumentError.
3. Sensible prec defaults (106 for cb-norm brackets; 256 for the pipeline,
   matching the C tests' CPREC).
4. **Gate:** serial julia smoke on live channels (GADC, phi_t, paper example) —
   the §1.6 headline tour runs end-to-end.

### aic-exa.8 [J5] — MOSEK extension
1. Move `sdp.jl` content into `ext/AICMosekExt.jl` verbatim (the five `diamond_*`).
2. Add `AlmostIdempotentChannels._diamond_value_impl(Φ, prec)` (the real
   `idempotency_defect`) and the tight-bracket path
   (`aic_cbnorm_certify_d` fed MOSEK feasible points — `diamond_norm_watrous_primal`
   + `diamond_norm_dual`). Mind the MOSEK-hang risk (aic-jhe): tightened tols
   (already in sdp.jl), bound the solves, fail loud on non-OPTIMAL.
3. **Gate:** load the ext, run one anchor (`idempotency_defect(phit_kraus(2,0.3))
   ≈ 0.315`); core still green without it.

### aic-exa.9 [T] — test restructure + Aqua/JET
1. Split `test/runtests.jl`: always-run core testset (`test_core.jl`,
   `test_factorize.jl`); SDP testset (`test_sdp.jl`, the current T5) guarded by a
   try/`using Mosek` + `@test_skip` (research §1.4).
2. Core green with NO MOSEK: T1 ccall round-trip; T2 anchors via `certified_defect`
   containment (not the SDP value — that moved to the gated suite); T3 bracket;
   T7 factorize round-trip (η=0 oracle + η>0 oblique, `delups`/`upsdel` straddle 0
   / are O(η)); T6 `UCPMap` validation; T8 `CertifiedBracket` containment.
3. `test_aqua.jl`: `Aqua.test_all` (ambiguities, stale_deps, deps_compat,
   piracies; undocumented_names off until docstrings settle). JET opportunistic,
   not a gate (§8.2).
4. Mutation-proof ≥1 core assertion (perturb → RED → restore).
5. **Gate:** `Pkg.test()` green WITHOUT a solver (one process).

### aic-exa.10 [X] — docs + DX polish
1. Documenter site: `index.md`, getting-started, `api.md` (`@autodocs`),
   `tutorial.md` telling the headline story (define near-idempotent Φ → measure
   `certified_defect` → `factorize` → round-trip).
2. README 60-second tour; every export has a runnable `jldoctest`.
3. Polish `show` aesthetics + error messages; concrete numbers, no marketing
   (Rule 12).
4. **Gate:** doctests pass (serial julia).

---

## 10. HONORING THE FINDINGS TRAPS (cross-cutting)

- **Choi–Effros star, not plain XY (§C2):** the star is never exposed as a naive
  Julia op in 0.2; all multiplication stays inside the C pipeline. The η>0
  oblique test fixture is mandatory wherever a star-touching artifact is emitted
  (the η=0 identity oracle is structurally blind to star bugs).
- **A is oblique (§C3):** `EpsCStarAlgebra.basis` is a snapshot only; projection
  into A is Φ̃ (done in C), never a Julia Frobenius projector.
- **Axioms hold only up to ε (callout):** no `EpsCStarAlgebra` method implies an
  exact unit or exact associativity; `epsilon(A)` surfaces the defect as a
  bracket so the approximate-ness is visible, not hidden.
- **cb-norm ≠ op-norm ≠ Frobenius (callout, §C12):** every `*_defect` /
  `cbnorm_*` is a cb-norm quantity; the value needs the SDP (extension), the
  bracket is the eig-free arb. The main-iso cb brackets use the OPERATOR-norm
  inclusion inf, not the Frobenius σ_min ampliation (§C12, "test that cannot
  fail").
- **Choi Convention-A in marshalling (callout):** the K-major Kraus layout
  `[a*n*n + i*n + j]` and Convention-A Choi `[p*N+q]` (conjugation on the FIRST
  factor) are preserved by reusing the existing `_flatten_kraus` and the shim
  layouts verbatim; the round-trip test (`choi_diff` vs pure-Julia) guards it.
- **UCP vs CPTP direction (callout, §C13):** `encode`/`decode` bind to
  `aic_factorize_dec_kraus`/`_enc_kraus` with asserted dual dims; the round-trip
  test catches a swap.
- **Universal, dim-independent constant (§D2):** `main_isomorphism`/`factorize`
  brackets must not grow with n — the dim-independence canary is a [C] test and a
  Julia test across n=2,3,4.
- **Fail loud out of basin (Rule 4):** `associated_algebra`/`factorize` surface
  the C `prop_P` basin abort as a Julia error naming ρ(Φ²−Φ) and pointing the
  user at `certified_defect`; no silent NaN.

---

## Appendix — answers to the prompt's nine required specs (index)

1. Public API surface → §1 (every export, signature, return type, docstring
   intent, example).
2. Type system → §2 (struct fields, immutability, parametrization, show sketches).
3. Extension architecture → §3 (exact Project.toml, ext layout, dispatch,
   graceful degradation).
4. ABI strategy → §4 (flat-double vs handle decision + rationale; the five new C
   shim signatures, what each wraps/returns, ownership).
5. File/module layout → §5.
6. Build/dev story → §6 (Preferences + make/cmake fallback; JLL as later
   milestone).
7. Certified-results ergonomics → §7 (`CertifiedBracket`; no IntervalArithmetic,
   justified).
8. Version/release plan → §8 (semver; 0.2 contents).
9. Ordered implementation checklist mapped to beads → §9 (.3 [C] … .10 [X]).

---

## Appendix B — orchestrator review addendum (ground-truthed 2026-06-02)

Reviewed against `include/aic/aic_factorize.h` and `.tex:2159` first-hand. The design is
**adopted** with these CORRECTIONS/confirmations. Where this addendum conflicts with §1–§10,
**this addendum wins** (it is ground-truthed).

**B1 — CORRECTED encode/decode ↔ Dec/Enc binding (the §C13 dual-direction trap; §1.6 + §2.4
prose were INVERTED).** Ground truth (`aic_factorize.h:306-316`):
- Δ ("Delta") = the UCP *encode* OBSERVABLE map B→B(H); Υ ("Upsilon") = the UCP *decode*
  OBSERVABLE map B(H)→B. `ΔΥ≈Φ̃`, `ΥΔ≈1_B`.
- The STATE channels are the duals: `Dec = Δ*` (channel B(H)→B, "produces a state on B",
  Kraus n_B×N, `aic_factorize_dec_kraus`); `Enc = Υ*` (channel B→B(H), Kraus N×n_B,
  `aic_factorize_enc_kraus`).
- Therefore the EXPORTED Julia verbs bind (note the names line up):
  - **`decode(F)` = Dec = Δ\* = `aic_factorize_dec_kraus`** (B(H)→B; dim_K=n_B, dim_H=N).
  - **`encode(F)` = Enc = Υ\* = `aic_factorize_enc_kraus`** (B→B(H); dim_K=N, dim_H=n_B).
  This matches CLAUDE.md ("decoding channel producing a state on B; encoding channel").
  [C] (C5) must tag each Kraus set with (dim_K,dim_H); [J4] asserts the dims and a
  `decode∘encode ≈ 1_B` round-trip so a swap goes RED (the design's safeguard is correct).
  In `struct Factorization`, store the Kraus by their CHANNEL role and dims, do not conflate
  Δ (observable) with Dec (channel).

**B2 — solver-free brackets CONFIRMED present (no new C math, just wrapping):**
- `aic_factorize_eigfree_delups` / `_upsdel` (arb lo,hi) already give the rigorous eig-free
  `‖ΔΥ−Φ‖_cb` / `‖ΥΔ−1_B‖_cb` brackets — C5 just rounds them FLOOR/CEIL to doubles.
- The general `aic_cbnorm_eigfree_ball_choi(J, n)` gives an eig-free two-sided bracket for ANY
  Choi J. So C3's `‖v‖_cb` / `‖v⁻¹‖_cb` brackets are obtained solver-free by forming the
  Convention-A Choi J(v) / J(v⁻¹) (opspace adjoint Choi, built directly per the O2 GOLDEN RULE,
  NOT by transpose) and calling that core. If no opspace eig-free entry point exists yet, [C]
  forms J(v) and calls `aic_cbnorm_eigfree_ball_choi` directly. Escalate ONLY if a solver-free
  rigorous UPPER bound on ‖v‖_cb proves unavailable (not expected).

**B3 — verify the weakdeps-as-test-deps pattern on Julia 1.10 ([J1]/[T]).** §3.1 lists
Convex/Mosek/MosekTools in BOTH `[weakdeps]` and `[extras]`+`[targets].test` so `Pkg.test`
activates the extension. This is the documented pattern but has historically been finicky.
[J1] must EMPIRICALLY confirm `Pkg.resolve()`/`Pkg.test()` succeeds on 1.10 with this layout;
if it does not, fall back to a dedicated `test/Project.toml` test environment that depends on
the weakdeps. Do not assume — verify in a serial julia run.

**B4 — minor:** `adjoint(Φ)`, `encode(F)`, `decode(F)` build `UCPMap`-shaped Kraus containers
whose duals are TRACE-PRESERVING, not unital — construct them with `check=false` (the
round-trip defect, not the constructor, is their certificate). Strip the stray "wait:" prose
from §1.6 when [X] polishes docs.

**B6 — LinearAlgebra export clashes (caught in [J2] smoke; the design's §1.1 exports collide).** The
exported `Factorization` collides with `LinearAlgebra.Factorization` and the verb `factorize` collides with
`LinearAlgebra.factorize` — a user doing `using AlmostIdempotentChannels, LinearAlgebra` (the common case,
since they build ComplexF64 matrices) gets ambiguous bindings. Resolution (authoritative, supersedes §1.1/
§1.6/§2.4 naming):
- Rename the type `Factorization` → **`ChannelFactorization`** (exported; clash-free). Used everywhere the
  body said `Factorization`.
- The verb stays `factorize` but is provided by **extending `LinearAlgebra.factorize`** (`import
  LinearAlgebra: factorize`; add `factorize(Φ::UCPMap; …)::ChannelFactorization`; `export factorize`). Because
  we re-export LinearAlgebra's *same* binding, `using both` is NOT ambiguous (ambiguity only arises between
  two *different* bindings of one name). Not piracy — the method's first arg is our `UCPMap` (Aqua-clean).
- [J2] does the type rename + audits `names(AlmostIdempotentChannels) ∩ names(LinearAlgebra)` for any other
  different-binding clash. [J4] does the `factorize` extension. [X] reflects the names in docs/doctests.

**B5 — adopted as-is:** flat-double ABI (no handles, §4.1); reuse of the factorize pipeline
build (§4.2); the C2–C5 signatures (§4.3) MODULO B1/B2; `CertifiedBracket` without
IntervalArithmetic (§7); the file layout (§5); Preferences libpath + dlopen/RTLD verbatim
(§6.1); MOSEK extension via `_diamond_value_impl` + error hints (§3); the bead-by-bead
checklist (§9). Serial-Julia + no-build-during-julia is the global execution invariant.
