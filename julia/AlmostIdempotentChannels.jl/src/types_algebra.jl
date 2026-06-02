# types_algebra.jl — the three algebra-side value-types (design doc §2): the
# associated ε-C* algebra `EpsCStarAlgebra`, the genuine target `CStarAlgebra`,
# and the O(ε)-isomorphism `MainIsomorphism` between them. Split out of the
# original `types.jl` (house Rule 10, ≤200 LOC). These reference
# `CertifiedBracket` and `UCPMap` (as field types), so this file is included
# AFTER types_core.jl. Same conventions as types_core.jl: immutable structs,
# concrete storage, accessors next to their type, no ccall (bead J2).

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
