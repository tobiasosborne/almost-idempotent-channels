# types.jl — the six immutable value-types of the package (design doc §2).
#
# These read like the paper: the nouns are the math objects (UCPMap,
# EpsCStarAlgebra, CStarAlgebra, MainIsomorphism, ChannelFactorization) and the
# rigorous scalar result (CertifiedBracket). ALL are immutable `struct`s with concrete
# ComplexF64/Float64 storage; only CertifiedBracket is parametrised (on the float
# type). No mutable state, no finalizers — the artifacts are read-once value
# snapshots (design §4.1: flat-double ABI, no opaque handles).
#
# THIS BEAD (J2) is PURE JULIA: no ccall. The constructors that the C pipeline
# populates (EpsCStarAlgebra/MainIsomorphism/ChannelFactorization) take their fields
# directly so bead J4 (api.jl) can call them with C-derived data; here we only
# DEFINE them + sane constructors. UCPMap's Kraus/Choi constructors ARE pure
# Julia (LinearAlgebra only — for construction/validation; the *certified* Choi
# path is C/J4).
#
# Accessors live next to their type. Name clashes with Base (`dim`, `value`,
# `width`, `midpoint`, `epsilon`) are NOT imported from Base — they are defined
# fresh in the module namespace and exported (avoiding type piracy). `Base.in`,
# `Base.minimum`, `Base.maximum`, `Base.adjoint` ARE legitimate Base extensions
# (we add methods for our own types only).

# ============================================================================
# CertifiedBracket — the rigorous two-sided result type (design §2.5, §1.7)
# ============================================================================

"""
    CertifiedBracket{T<:AbstractFloat}

A RIGOROUS two-sided bracket `lo ≤ x ≤ hi` on a scalar (the C arb balls rounded
outward, FLOOR/CEIL — aic_ucp_shim.h:46-52, so the bound survives the double
conversion). Optional `value` is a point estimate (e.g. the Watrous SDP value)
when one is available; `nothing` for the solver-free brackets that are the
package default. `label` names the quantity for `show` ("‖Φ²−Φ‖_cb", "ε_assoc",
…).

    CertifiedBracket(lo, hi; value=nothing, label="")

The inner constructor validates `lo ≤ hi` (else `ArgumentError`) and, when a
`value` is given, that it lies in `[lo-1e-9, hi+1e-9]` (a 1e-9 slack absorbs the
benign double-rounding of an SDP point estimate against an outward-rounded
bracket).
"""
struct CertifiedBracket{T<:AbstractFloat}
    lo::T
    hi::T
    value::Union{T,Nothing}
    label::String
    function CertifiedBracket(lo::T, hi::T; value=nothing, label::AbstractString="") where {T<:AbstractFloat}
        lo ≤ hi || throw(ArgumentError("CertifiedBracket: lo=$lo > hi=$hi (a bracket must have lo ≤ hi)"))
        v = value === nothing ? nothing : convert(T, value)
        v === nothing || (lo - 1e-9 ≤ v ≤ hi + 1e-9) ||
            throw(ArgumentError("CertifiedBracket: value=$v outside [$lo, $hi]"))
        new{T}(lo, hi, v, String(label))
    end
end

# Promote a mixed (lo, hi) pair to a common AbstractFloat type so
# CertifiedBracket(0, 1.0) and CertifiedBracket(0.0, 1) both work.
function CertifiedBracket(lo::Real, hi::Real; value=nothing, label::AbstractString="")
    T = float(promote_type(typeof(lo), typeof(hi)))
    vv = value === nothing ? nothing : convert(T, value)
    return CertifiedBracket(convert(T, lo), convert(T, hi); value=vv, label=label)
end

Base.in(x::Real, b::CertifiedBracket) = b.lo ≤ x ≤ b.hi
Base.minimum(b::CertifiedBracket) = b.lo
Base.maximum(b::CertifiedBracket) = b.hi

"""
    width(b::CertifiedBracket) -> AbstractFloat

The bracket width `hi − lo` — how loose the certificate is.
"""
width(b::CertifiedBracket) = b.hi - b.lo

"""
    midpoint(b::CertifiedBracket) -> AbstractFloat

The bracket midpoint `(lo + hi) / 2`.
"""
midpoint(b::CertifiedBracket) = (b.lo + b.hi) / 2

"""
    value(b::CertifiedBracket) -> Union{AbstractFloat,Nothing}

The point estimate (e.g. the Watrous SDP value) if one was supplied, else
`nothing` (the solver-free brackets carry no point estimate).
"""
value(b::CertifiedBracket) = b.value

# ============================================================================
# UCPMap — the channel, the universal input (design §2.1, §1.2)
# ============================================================================

"""
    UCPMap(kraus::Vector{<:AbstractMatrix}; check=true, atol=1e-9)
    UCPMap(; choi::AbstractMatrix, n::Int, check=true, atol=1e-9)

A unital completely-positive (UCP) self-map Φ: B(H) → B(H), dim H = n, in the
HEISENBERG picture Φ(X) = Σ_a K_a† X K_a (observables; aic_ucp.h:8). `kraus` are
the {K_a}, each n×n. Unitality (Σ_a K_a† K_a = I, asserted up to `atol` when
`check`) is the defining UCP-self-map hypothesis the whole pipeline needs.

Validation: empty kraus → `ArgumentError`; ragged / non-square → `DimensionMismatch`;
non-unital beyond `atol` (when `check`) → `ArgumentError` naming the measured
`‖Σ K†K − I‖`.

The Choi constructor reconstructs {K_a} from a Convention-A Choi
`C[(i-1)n+a, (j-1)n+b] = Σ_x conj(K_x[i,a]) K_x[j,b]` (conjugation on the FIRST
factor, aic_ucp.h:36) by a tolerance-aware PSD Hermitian eigendecomposition
(LinearAlgebra only — this is for construction/validation; the certified path is
C/J4).
"""
struct UCPMap
    kraus::Vector{Matrix{ComplexF64}}   # r operators, each n×n (K_a: H→H, Heisenberg)
    n::Int                              # dim H
    r::Int                              # number of Kraus ops
end

function UCPMap(kraus::AbstractVector{<:AbstractMatrix}; check::Bool=true, atol::Real=1e-9)
    isempty(kraus) && throw(ArgumentError("UCPMap: empty Kraus list (need ≥ 1 operator)"))
    nn = size(kraus[1], 1)
    for (a, K) in enumerate(kraus)
        size(K, 1) == size(K, 2) ||
            throw(DimensionMismatch("UCPMap: Kraus op $a is $(size(K)), not square"))
        size(K, 1) == nn ||
            throw(DimensionMismatch("UCPMap: Kraus op $a is $(size(K)), expected ($nn, $nn) to match op 1"))
    end
    K64 = Matrix{ComplexF64}[Matrix{ComplexF64}(K) for K in kraus]
    if check
        S = sum(K' * K for K in K64)
        defect = opnorm(S - I)
        defect ≤ atol || throw(ArgumentError(
            "UCPMap: not unital — ‖Σ K†K − I‖ = $defect > atol = $atol " *
            "(Heisenberg unitality Σ K_a† K_a = I is required; pass check=false " *
            "to skip, e.g. for the trace-preserving Dec/Enc duals)"))
    end
    return UCPMap(K64, nn, length(K64))
end

function UCPMap(; choi::AbstractMatrix, n::Integer, check::Bool=true, atol::Real=1e-9)
    n ≥ 1 || throw(ArgumentError("UCPMap(choi=…): n must be ≥ 1, got $n"))
    size(choi, 1) == size(choi, 2) ||
        throw(DimensionMismatch("UCPMap(choi=…): Choi must be square, got $(size(choi))"))
    size(choi, 1) == n * n ||
        throw(DimensionMismatch("UCPMap(choi=…): Choi is $(size(choi)), expected ($(n*n), $(n*n)) for n=$n"))
    C = Matrix{ComplexF64}(choi)
    # Convention-A Choi C[(i-1)n+a, (j-1)n+b] = Σ_x conj(K_x[i,a]) K_x[j,b].
    # Hermitise (the input may carry rounding asymmetry), then PSD-eigendecompose.
    H = Hermitian((C + C') / 2)
    F = eigen(H)                         # ascending eigenvalues, λ_k with vectors v_k
    λmax = isempty(F.values) ? 0.0 : maximum(F.values)
    tol = max(atol, 1e-12 * max(λmax, 1.0))
    krs = Matrix{ComplexF64}[]
    for k in eachindex(F.values)
        λ = F.values[k]
        λ > tol || continue              # drop the (near-)zero / negative PSD tail
        v = F.vectors[:, k]
        # v is a vectorised Kraus op in the SAME index order as the Choi: the
        # composite index p = (i-1)n + a (1-based) carries row i, col a. So
        # K_x[i,a] = sqrt(λ) * v[(i-1)n + a].  (Convention-A: conj on the first
        # factor (i,a) means the column space of C, i.e. the eigenvectors, are
        # the un-conjugated K_x in (i,a) order.)
        s = sqrt(λ)
        K = Matrix{ComplexF64}(undef, n, n)
        for i in 1:n, a in 1:n
            K[i, a] = s * v[(i - 1) * n + a]
        end
        push!(krs, K)
    end
    isempty(krs) && (krs = Matrix{ComplexF64}[zeros(ComplexF64, n, n)])  # the zero map
    return UCPMap(krs; check=check, atol=atol)
end

"""
    n(Φ::UCPMap) -> Int

The dimension of the Hilbert space H (Φ acts on B(H), each Kraus op is n×n).
"""
n(Φ::UCPMap) = Φ.n

"""
    nkraus(Φ::UCPMap) -> Int

The number r of Kraus operators {K_a}.
"""
nkraus(Φ::UCPMap) = Φ.r

"""
    kraus(Φ::UCPMap) -> Vector{Matrix{ComplexF64}}

The Kraus operators {K_a} of Φ (Heisenberg picture, Φ(X) = Σ_a K_a† X K_a).
"""
kraus(Φ::UCPMap) = Φ.kraus

"""
    isunital(Φ::UCPMap; atol=1e-9) -> Bool

Whether Φ is unital to `atol`, i.e. `‖Σ_a K_a† K_a − I‖ ≤ atol`. Heisenberg
unitality is the UCP-self-map hypothesis; the Dec/Enc duals (built with
`check=false`) are trace-preserving instead and report `false` here.
"""
function isunital(Φ::UCPMap; atol::Real=1e-9)
    S = sum(K' * K for K in Φ.kraus)
    return opnorm(S - I) ≤ atol
end

"""
    adjoint(Φ::UCPMap) -> UCPMap

The dual map with Kraus {K_a†} (Appendix B4). The dual of a UCP Heisenberg map
is the trace-preserving Schrödinger channel, so it is built with `check=false`
(its certificate is trace-preservation / the round-trip defect, not Heisenberg
unitality). Usable as `Φ'`.
"""
Base.adjoint(Φ::UCPMap) = UCPMap(Matrix{ComplexF64}[Matrix(K') for K in Φ.kraus]; check=false)

# ============================================================================
# EpsCStarAlgebra — the associated ε-C* algebra A = Img Φ̃ (design §2.2, §1.4)
# ============================================================================

"""
    EpsCStarAlgebra(n, dim_A, basis, eps, source)

The associated ε-C* algebra A = Img Φ̃ (Φ̃ = θ(2Φ−1) the exact-idempotent
regularisation; th_almost_idemp, approximate_algebras.tex:2184). A is an OBLIQUE
subspace of M_n (FINDINGS §C3); its product is the Choi–Effros star X⋆Y = Φ̃(XY),
NOT plain XY (FINDINGS §C2); the axioms hold only UP TO ε (no exact unit).

Fields: `n` ambient dim (A ≤ M_n), `dim_A`, `basis` the Frobenius-orthonormal
{B_k} snapshot (n×n each), `eps` the certified associator-defect bracket, and
`source` the UCPMap it came from. The basis is a snapshot for showability; the
algebra's operations are not re-done in Julia (main_isomorphism/factorize rebuild
through the C pipeline).
"""
struct EpsCStarAlgebra
    n::Int
    dim_A::Int
    basis::Vector{Matrix{ComplexF64}}
    eps::CertifiedBracket{Float64}
    source::UCPMap
end

"""
    dim(A::EpsCStarAlgebra) -> Int

The dimension dim A of the ε-C* algebra.
"""
dim(A::EpsCStarAlgebra) = A.dim_A

"""
    ambient(A::EpsCStarAlgebra) -> Int

The ambient dimension n (A is an oblique subspace of M_n).
"""
ambient(A::EpsCStarAlgebra) = A.n

"""
    epsilon(A::EpsCStarAlgebra) -> CertifiedBracket

The certified associator-defect ε bracket (how far A is from a genuine C*
algebra; the axioms hold only up to ε).
"""
epsilon(A::EpsCStarAlgebra) = A.eps

# ============================================================================
# CStarAlgebra — the genuine target C* algebra B = ⊕_l M_{d_l} (design §2.3, §1.5)
# ============================================================================

"""
    CStarAlgebra(blocks::AbstractVector{<:Integer})

A genuine finite-dimensional C* algebra B = ⊕_l M_{d_l}, given by its block sizes
`d_l`. The derived `dim_B = Σ d_l²` (vector-space dimension) and `n_B = Σ d_l`
(block-diagonal representative size) are computed by the constructor. A pure
structural summary of `aic_dhom_B` — no matrices needed.
"""
struct CStarAlgebra
    blocks::Vector{Int}   # d_l, the genuine C* block sizes
    dim_B::Int            # Σ d_l²
    n_B::Int              # Σ d_l
end

function CStarAlgebra(blocks::AbstractVector{<:Integer})
    isempty(blocks) && throw(ArgumentError("CStarAlgebra: empty block list (need ≥ 1 block)"))
    all(d -> d ≥ 1, blocks) || throw(ArgumentError("CStarAlgebra: every block size d_l must be ≥ 1, got $blocks"))
    b = Int[Int(d) for d in blocks]
    return CStarAlgebra(b, sum(d -> d^2, b), sum(b))
end

"""
    blocks(B::CStarAlgebra) -> Vector{Int}

The C* block sizes d_l of B = ⊕_l M_{d_l}.
"""
blocks(B::CStarAlgebra) = B.blocks

"""
    dim_B(B::CStarAlgebra) -> Int

The vector-space dimension dim_B = Σ d_l².
"""
dim_B(B::CStarAlgebra) = B.dim_B

"""
    n_B(B::CStarAlgebra) -> Int

The block-diagonal representative size n_B = Σ d_l.
"""
n_B(B::CStarAlgebra) = B.n_B

# ============================================================================
# MainIsomorphism — the O(ε)-isomorphism v: B → A (design §2.4, §1.5)
# ============================================================================

"""
    MainIsomorphism(source, B, cbnorm_fwd, cbnorm_inv, isodefect)

The O(ε)-isomorphism v: B → A of th_main (approximate_algebras.tex:460), with a
UNIVERSAL, dimension-INDEPENDENT constant. `B` is the genuine C* algebra, and the
three `CertifiedBracket` fields are the certified cb-norm brackets ‖v‖_cb,
‖v⁻¹‖_cb (≈ 1 + O(ε)) and the isomorphism defect. `source` is the originating
UCPMap.
"""
struct MainIsomorphism
    source::UCPMap
    B::CStarAlgebra
    cbnorm_fwd::CertifiedBracket{Float64}
    cbnorm_inv::CertifiedBracket{Float64}
    isodefect::CertifiedBracket{Float64}
end

"""
    blocks(v::MainIsomorphism) -> Vector{Int}

The C* block sizes d_l of the target algebra B.
"""
blocks(v::MainIsomorphism) = v.B.blocks

"""
    cbnorm_forward(v::MainIsomorphism) -> CertifiedBracket

The certified bracket on ‖v‖_cb.
"""
cbnorm_forward(v::MainIsomorphism) = v.cbnorm_fwd

"""
    cbnorm_inverse(v::MainIsomorphism) -> CertifiedBracket

The certified bracket on ‖v⁻¹‖_cb.
"""
cbnorm_inverse(v::MainIsomorphism) = v.cbnorm_inv

"""
    isodefect(v::MainIsomorphism) -> CertifiedBracket

The certified isomorphism-defect bracket.
"""
isodefect(v::MainIsomorphism) = v.isodefect

# ============================================================================
# CPMap — a RECTANGULAR CPTP map (Appendix B7/B10); UCPMap is square-only
# ============================================================================

"""
    CPMap(kraus::Vector{<:AbstractMatrix}, dim_in::Int, dim_out::Int; check=true, atol=1e-9)

A completely-positive trace-preserving (CPTP) quantum channel between possibly
DIFFERENT spaces, acting on states in the SCHRÖDINGER picture

    ρ ↦ Σ_a K_a ρ K_a†          (K_a is `dim_out` × `dim_in`).

This is the type of the encode/decode maps of a [`factorize`](@ref): they map
between the algebra B (size `n_B`) and B(H) (size `N`), so their Kraus are
RECTANGULAR and they do NOT fit the square self-map `UCPMap` (Appendix B7).

The type is named `CPMap` (completely-positive map; mirrors `UCPMap` — an
observable unital-CP self-map ↔ a state CP channel), NOT `Channel`: a type
named `Channel` would collide with `Base.Channel` (the concurrency primitive),
making a bare `Channel` ambiguous under `using` (Appendix B10).

Trace preservation `Σ_a K_a† K_a = I_{dim_in}` is asserted up to `atol` when
`check` (the Schrödinger dual of Heisenberg unitality). Empty/ragged/mis-shaped
Kraus throw `ArgumentError`/`DimensionMismatch` naming the offending operator.
"""
struct CPMap
    kraus::Vector{Matrix{ComplexF64}}   # K_a: dim_out × dim_in (Schrödinger)
    dim_in::Int
    dim_out::Int
end

function CPMap(kraus::AbstractVector{<:AbstractMatrix}, dim_in::Integer, dim_out::Integer;
                 check::Bool=true, atol::Real=1e-9)
    isempty(kraus) && throw(ArgumentError("CPMap: empty Kraus list (need ≥ 1 operator)"))
    (dim_in ≥ 1 && dim_out ≥ 1) ||
        throw(ArgumentError("CPMap: dim_in=$dim_in, dim_out=$dim_out must both be ≥ 1"))
    for (a, K) in enumerate(kraus)
        size(K) == (dim_out, dim_in) || throw(DimensionMismatch(
            "CPMap: Kraus op $a is $(size(K)), expected ($dim_out, $dim_in) " *
            "(K_a maps the input space dim_in into the output space dim_out, Schrödinger)"))
    end
    K64 = Matrix{ComplexF64}[Matrix{ComplexF64}(K) for K in kraus]
    if check
        S = sum(K' * K for K in K64)             # Σ K_a† K_a, must be I_{dim_in}
        defect = opnorm(S - I)
        defect ≤ atol || throw(ArgumentError(
            "CPMap: not trace-preserving — ‖Σ K_a† K_a − I_$(dim_in)‖ = $defect > atol = $atol " *
            "(CPTP requires Σ K_a† K_a = I on the INPUT space; pass check=false to skip)"))
    end
    return CPMap(K64, Int(dim_in), Int(dim_out))
end

"""
    iscptp(c::CPMap; atol=1e-9) -> Bool

Whether `c` is trace-preserving to `atol`, i.e. `‖Σ_a K_a† K_a − I_{dim_in}‖ ≤ atol`
(the defining CPTP condition on the INPUT space).

CAVEAT for a [`factorize`](@ref) `decode` channel (Dec = Δ*) with η>0: its
trace-preservation is only O(η)-exact (the §C14 PSD-cone clip — measured TP-defect
≈ 1e-6 on multi-block oblique fixtures), so `iscptp` is `false` at the default
`atol = 1e-9` and you must loosen `atol` (~1e-5) to confirm. The rigorous
round-trip certificate is the cb-norm bracket [`upsdel_defect`](@ref), not this
machine-tolerance predicate. The `encode` channel and the η=0 oracle decode are
TP to machine precision.
"""
function iscptp(c::CPMap; atol::Real=1e-9)
    S = sum(K' * K for K in c.kraus)
    return opnorm(S - I) ≤ atol
end

# ============================================================================
# ChannelFactorization — the headline Φ ≈ Δ Υ factorization (design §2.4, §1.6, B1)
# ============================================================================

"""
    ChannelFactorization(source, B, encode, decode, delups, upsdel, eta_proxy, eps_used)

The approximate channel factorization Φ ≈ Δ Υ through a genuine C* algebra B
(th_factorization, approximate_algebras.tex:2730-2899). The observable maps
Δ: B → B(H) (encode) and Υ: B(H) → B (decode) satisfy the by-construction targets
ΔΥ ≈ Φ̃, ΥΔ ≈ 1_B, certified by two cb-norm round-trip brackets
`delups = ‖ΔΥ − Φ‖_cb` (.tex:2733) and `upsdel = ‖ΥΔ − 1_B‖_cb` (.tex:2739).
`eta_proxy`/`eps_used` are the shim sentinels (≈ 0 at η = 0).

The STORED objects are the state-space CHANNELS — the duals (Appendix B1, B7):
the `encode` CPMap = Enc = Υ* : B → B(H) (CPTP, n_B→N) and the `decode` CPMap
= Dec = Δ* : B(H) → B (CPTP, N→n_B). These are RECTANGULAR CPTP maps between
DIFFERENT spaces, so they are `CPMap`s, NOT square `UCPMap`s. The accessors
[`encode`](@ref)/[`decode`](@ref) expose them.

The name is `ChannelFactorization`, NOT `Factorization`: the latter clashes with
the abstract type `LinearAlgebra.Factorization`, which `using` re-exports — a user
doing `using AlmostIdempotentChannels, LinearAlgebra` would otherwise hit an
ambiguous binding. The verb (bead J4) shares LinearAlgebra's binding by EXTENDING
`LinearAlgebra.factorize` to return a `ChannelFactorization`, so it never clashes.
"""
struct ChannelFactorization
    source::UCPMap
    B::CStarAlgebra
    encode::CPMap   # Enc = Υ* : B → B(H)   (CPTP, dim_in=n_B, dim_out=N)
    decode::CPMap   # Dec = Δ* : B(H) → B   (CPTP, dim_in=N,   dim_out=n_B)
    delups::CertifiedBracket{Float64}   # ‖ΔΥ − Φ‖_cb
    upsdel::CertifiedBracket{Float64}   # ‖ΥΔ − 1_B‖_cb
    eta_proxy::Float64
    eps_used::Float64
end

"""
    algebra(F::ChannelFactorization) -> CStarAlgebra

The genuine C* algebra B = ⊕_l M_{d_l} the factorization passes through.
"""
algebra(F::ChannelFactorization) = F.B

"""
    delups_defect(F::ChannelFactorization) -> CertifiedBracket

The certified bracket on ‖ΔΥ − Φ‖_cb (encode∘decode returns to Φ within O(η)).
"""
delups_defect(F::ChannelFactorization) = F.delups

"""
    upsdel_defect(F::ChannelFactorization) -> CertifiedBracket

The certified bracket on ‖ΥΔ − 1_B‖_cb (decode∘encode returns to 1_B within O(η)).
"""
upsdel_defect(F::ChannelFactorization) = F.upsdel
