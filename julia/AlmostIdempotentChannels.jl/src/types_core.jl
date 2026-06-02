# types_core.jl — the two foundational value-types (design doc §2): the rigorous
# scalar `CertifiedBracket` and the channel input `UCPMap`. Split out of the
# original `types.jl` (house Rule 10, ≤200 LOC); types_algebra.jl/types_channel.jl
# reference these as field types, so this is included FIRST. PURE JULIA, no ccall.

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
