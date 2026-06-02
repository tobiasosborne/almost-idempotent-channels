# sdp_stubs.jl — solver-gated public surface as helpful-error STUBS.
#
# The EXACT cb-norm value (the Watrous diamond-norm SDP) needs Convex + Mosek +
# MosekTools and lives ONLY in the package extension AICMosekExt. The core
# declares the public names here so they are discoverable and Documenter-visible,
# with the implementation gated behind an internal dispatch hook the extension
# adds a method to (design doc §3.3).
#
# Pattern (register-by-method-override): the PUBLIC verb idempotency_defect calls
# the internal _diamond_value_impl; the core defines _diamond_value_impl's
# helpful-error fallback; the extension adds the REAL method, which is more
# specific (concrete arg types) and wins via Julia's normal method table. The
# diamond_* names are likewise stubs (the extension shadows them). __init__ also
# registers a MethodError hint (libaic.jl) as belt-and-suspenders.
#
# SIGNATURES: the type system (UCPMap) is a LATER bead (J2). For 0.2-J1 we stub
# on the RAW signatures the current surface uses — kraus as
# Vector{Matrix{ComplexF64}} (idempotency_defect / _diamond_value_impl), and the
# Convention-A Choi Matrix{ComplexF64} + dim n for the diamond_* SDP entry points
# (matching ext/AICMosekExt.jl's eventual src/sdp.jl methods). When UCPMap lands,
# J5 re-points idempotency_defect at it; the hint mechanism is unchanged.

const _MOSEK_HINT = """requires the MOSEK extension. Install + load the solver:
    using Pkg; Pkg.add(["Convex","Mosek","MosekTools"])
    using Convex, Mosek, MosekTools
(Mosek needs a license: https://www.mosek.com/products/academic-licenses/)"""

"""
    idempotency_defect(kraus::Vector{Matrix{ComplexF64}}; prec=106) -> Float64

The EXACT idempotency defect η = ‖Φ²−Φ‖_cb = ‖Φ²−Φ‖_⋄ of the UCP self-map Φ
given by `kraus` (Heisenberg picture), via the Watrous diamond-norm SDP
(approximate_algebras.tex:347-354). Requires Convex + Mosek + MosekTools (the
AICMosekExt package extension); WITHOUT them this throws a helpful error pointing
at the install step. The SOLVER-FREE rigorous bracket is `eta_eigfree` (the
package default for the defect). See design doc §1.3, §3.3.
"""
idempotency_defect(kraus::Vector{Matrix{ComplexF64}}; prec::Int = 106)::Float64 =
    _diamond_value_impl(kraus, prec)

# Internal dispatch hook. The extension adds the real (SDP) method; this fallback
# fires only when the extension is not loaded.
_diamond_value_impl(::Any, ::Int) = error("idempotency_defect " * _MOSEK_HINT)

# ---- the Watrous diamond-norm SDP entry points (extension-only) ----
# Stubbed on the Convention-A Choi + dim signature src/sdp.jl uses; the extension
# shadows each with the real Convex.jl + MOSEK method (bead J5, aic-exa.8).

"""
    diamond_norm_watrous(J::Matrix{ComplexF64}, n::Int)

Watrous 2012 MAX-primal diamond-norm SDP for the map with Convention-A Choi `J`.
Extension-only (Convex + Mosek + MosekTools); throws a helpful error otherwise.
"""
diamond_norm_watrous(::Matrix{ComplexF64}, ::Int) =
    error("diamond_norm_watrous " * _MOSEK_HINT)

"""
    diamond_norm_watrous_primal(J::Matrix{ComplexF64}, n::Int)

MAX-primal SDP exposing the feasible point (X,P,Q) for the arb lower-bound
certifier. Extension-only; throws a helpful error otherwise.
"""
diamond_norm_watrous_primal(::Matrix{ComplexF64}, ::Int) =
    error("diamond_norm_watrous_primal " * _MOSEK_HINT)

"""
    diamond_norm_dual(J::Matrix{ComplexF64}, n::Int)

QETLAB/Watrous MIN-dual SDP exposing (Y0,Y1) for the arb upper-bound certifier.
Extension-only; throws a helpful error otherwise.
"""
diamond_norm_dual(::Matrix{ComplexF64}, ::Int) =
    error("diamond_norm_dual " * _MOSEK_HINT)
