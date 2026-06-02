# types_channel.jl — the two state-space channel value-types (design doc §2): the
# rectangular CPTP channel `CPMap` and the headline `ChannelFactorization`
# (Φ ≈ Δ Υ). Split out of the original `types.jl` (house Rule 10, ≤200 LOC).
# `ChannelFactorization` references `CStarAlgebra` (in types_algebra.jl), `CPMap`
# (defined here), and `CertifiedBracket` (in types_core.jl) as field types, so
# this file is included AFTER types_algebra.jl. Same conventions as the sibling
# types files: immutable structs, concrete storage, accessors next to their type,
# no ccall (bead J2).

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
