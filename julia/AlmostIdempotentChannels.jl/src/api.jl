# api.jl — the HIGH-LEVEL idiomatic API (bead aic-exa.7 [J4]): the layer humans
# and agents actually touch. The verbs read like the paper:
#
#     Φ = UCPMap(kraus)
#     η = certified_defect(Φ)        # rigorous ‖Φ²−Φ‖_cb bracket, NO solver (default)
#     A = associated_algebra(Φ)      # the ε-C* algebra A = Img Φ̃
#     v = main_isomorphism(Φ)        # O(ε)-isomorphism v: B → A, dim-independent
#     F = factorize(Φ)               # Φ ≈ Δ Υ through a genuine C* algebra B
#
# Built on the J3 internal wrappers (_assoc_summary, _main_iso_summary,
# _factorize_artifacts) and the J2 value-types. Design: julia_package_design.md
# §1 (the surface), §9 [J4] (the contract), Appendices B1/B6/B7/B8.
#
# THREE LOAD-BEARING CORRECTNESS TRAPS guarded across this file + basin.jl:
#   B1  the encode/decode ↔ Enc/Dec dual-direction binding (§C13): encode = Enc =
#       Υ* (B→B(H), N×n_B, from the `enc` buffer); decode = Dec = Δ* (B(H)→B,
#       n_B×N, from the `dec` buffer). NOT crossed. The dims are asserted in tests.
#   B7  encode/decode are RECTANGULAR CPTP maps between DIFFERENT spaces, so they
#       are a `CPMap` (dim_in≠dim_out), NOT a square `UCPMap` (the B10 type name:
#       `Channel` clashes with `Base.Channel`).
#   B8  the C pipeline HARD-ASSERTS (SIGABRT, kills Julia) when ρ(Φ²−Φ) ≥ 1/4 (the
#       prop_P basin, §C15). Julia CANNOT catch a C abort, so associated_algebra/
#       main_isomorphism/factorize PRE-CHECK ρ in Julia first (see basin.jl,
#       included before this file) and throw an ArgumentError out of basin.

# ============================================================================
# The defect verbs
# ============================================================================

"""
    certified_defect(Φ::UCPMap; prec=106) -> CertifiedBracket

A RIGOROUS two-sided bracket on the idempotency defect η = ‖Φ²−Φ‖_cb
(approximate_algebras.tex:347). Solver-free — the eig-free arb certifier
(`aic_cbnorm_eigfree_d`) rounded outward, so `lo ≤ η ≤ hi` holds to the arb
working precision. THE DEFAULT: works with no MOSEK and NEVER aborts (no basin
pre-check needed). The bracket is loose by design (`hi/lo ~ 2n`); it certifies a
value rather than computing it, so `value` is `nothing`. The exact η VALUE is
[`idempotency_defect`](@ref) (the Watrous SDP, MOSEK extension).
"""
function certified_defect(Φ::UCPMap; prec::Int=106)::CertifiedBracket
    lo, hi = eta_eigfree(Φ.kraus; prec=prec)
    return CertifiedBracket(lo, hi; label="‖Φ²−Φ‖_cb")
end

"""
    idempotency_defect(Φ::UCPMap; prec=106) -> Float64

The EXACT idempotency defect η = ‖Φ²−Φ‖_cb = ‖Φ²−Φ‖_⋄ via the Watrous
diamond-norm SDP (approximate_algebras.tex:347). Requires the MOSEK extension
(Convex + Mosek + MosekTools); WITHOUT it this throws a helpful install hint
(NOT a `MethodError`). The solver-free rigorous alternative is
[`certified_defect`](@ref) (the package default).
"""
idempotency_defect(Φ::UCPMap; prec::Int=106)::Float64 = _diamond_value_impl(Φ, prec)

# NOTE (review N-3): the helpful-error fallback for `_diamond_value_impl` is the
# single `(::Any, ::Int)` method in sdp_stubs.jl — it already covers `UCPMap`.
# We deliberately do NOT define a `(::UCPMap, ::Int)` fallback here: the J5 MOSEK
# extension adds exactly that concrete method, so a fallback on the identical
# signature would be REDEFINED/overwritten by the extension (a "method
# overwritten" coupling). With only the `(::Any, ::Int)` stub in the core, the
# extension's `(::UCPMap, ::Int)` is strictly more specific and wins cleanly with
# no overwrite (design §3.3 register-by-method-override).

# ============================================================================
# The pipeline verbs (all basin-guarded, Appendix B8)
# ============================================================================

"""
    associated_algebra(Φ::UCPMap; prec=256) -> EpsCStarAlgebra

The associated ε-C* algebra A = Img Φ̃, Φ̃ = θ(2Φ−1) the exact-idempotent
regularisation (th_almost_idemp, approximate_algebras.tex:2184). A is an OBLIQUE
subspace of M_n (FINDINGS §C3); its product is the Choi–Effros star X⋆Y = Φ̃(XY),
NOT plain XY (FINDINGS §C2); the axioms hold only UP TO ε (no exact unit).
Carries `dim A`, the ambient `n`, and a certified ε = associator-defect bracket.
Pre-checks the prop_P basin in Julia (Appendix B8) and throws an `ArgumentError`
naming ρ(Φ²−Φ) if Φ is out of basin, so the C pipeline cannot abort the session.
"""
function associated_algebra(Φ::UCPMap; prec::Int=256)::EpsCStarAlgebra
    _assert_in_basin(Φ)
    s = _assoc_summary(Φ.kraus, Φ.n, Φ.r; eps=-2.0, prec=prec)
    return EpsCStarAlgebra(s.n, s.dim_A, s.basis, s.eps, Φ)
end

"""
    main_isomorphism(Φ::UCPMap; prec=256) -> MainIsomorphism

The O(ε)-isomorphism v: B → A of th_main (approximate_algebras.tex:460), with B a
GENUINE finite-dim C* algebra B = ⊕_l M_{d_l} and a UNIVERSAL, dimension-
INDEPENDENT constant (FINDINGS §D2 — the constant must NOT grow with n). Carries
B's block sizes and the certified cb-norm brackets ‖v‖_cb, ‖v⁻¹‖_cb (≈ 1 + O(ε))
plus the isomorphism-defect bracket. Pre-checks the prop_P basin (Appendix B8).
"""
function main_isomorphism(Φ::UCPMap; prec::Int=256)::MainIsomorphism
    _assert_in_basin(Φ)
    s = _main_iso_summary(Φ.kraus, Φ.n, Φ.r; eps=-2.0, prec=prec)
    return MainIsomorphism(Φ, CStarAlgebra(s.blocks), s.cbfwd, s.cbinv, s.isodef)
end

"""
    factorize(Φ::UCPMap; prec=256) -> ChannelFactorization

The approximate channel factorization Φ ≈ Δ Υ through a genuine C* algebra B
(th_factorization, approximate_algebras.tex:2730-2899), certified by two cb-norm
round-trip brackets ‖ΔΥ − Φ‖_cb (`delups_defect`) and ‖ΥΔ − 1_B‖_cb
(`upsdel_defect`). The state-space CHANNELS are the duals (Appendix B1):

    encode(F) = Enc = Υ*  : B → B(H)   (dim_in=n_B, dim_out=N)
    decode(F) = Dec = Δ*  : B(H) → B   (dim_in=N,   dim_out=n_B)

This method EXTENDS `LinearAlgebra.factorize` (Appendix B6) — `using
AlmostIdempotentChannels, LinearAlgebra` is unambiguous because we share that one
binding (the first argument `Φ::UCPMap` is ours, so it is not type piracy).

DOMAIN (TIGHTER than the prop_P basin). factorize builds Υ via a power-series
functional calculus (`aic_funcalc_xpow` inside `upsilon_build`) whose convergence
domain is MUCH tighter than the prop_P basin ρ(Φ²−Φ) < 1/4 used by
`associated_algebra` / `main_isomorphism` (bug aic-exa.13). So this verb pre-checks
a SEPARATE, conservative threshold (`_FACTORIZE_RHO_MAX = 0.10`, Appendix B8 /
basin.jl) and throws an `ArgumentError` — naming the estimated ρ — when Φ is out of
domain, instead of letting the C pipeline SIGABRT the session. The factorization is
meaningful only in the small-η (O(η)) regime anyway; for larger η use
[`certified_defect`](@ref) (always safe), [`associated_algebra`](@ref), or
[`main_isomorphism`](@ref) (the wider ρ < 1/4 domain). The guard is BEST-EFFORT: a
severely out-of-domain channel of an untested shape could still abort (the abort
boundary is dimension-dependent and the ρ estimate can under-read — bug aic-exa.13).
"""
function LinearAlgebra.factorize(Φ::UCPMap; prec::Int=256)::ChannelFactorization
    _assert_factorize_domain(Φ)
    a = _factorize_artifacts(Φ.kraus, Φ.n, Φ.r; eps=-2.0, prec=prec)
    B = CStarAlgebra(a.blocks)
    # B1: encode = Enc = Υ* from the `enc` buffer (N×n_B, B→B(H)); decode = Dec =
    # Δ* from the `dec` buffer (n_B×N, B(H)→B). NOT crossed. CPTP duals are
    # trace-preserving, so build with check=false (the round-trip brackets, not the
    # constructor, are their certificate — Appendix B4).
    enc = CPMap(a.enc, B.n_B, a.N; check=false)   # B → B(H)
    dec = CPMap(a.dec, a.N, B.n_B; check=false)   # B(H) → B
    return ChannelFactorization(Φ, B, enc, dec, a.delups, a.upsdel, a.eta_proxy, a.eps_used)
end

# ============================================================================
# ChannelFactorization channel accessors (Appendix B1/B7)
# ============================================================================

"""
    encode(F::ChannelFactorization) -> CPMap

The encode CHANNEL Enc = Υ* : B → B(H) (CPTP, `dim_in = n_B`, `dim_out = N`).
Takes a state on the genuine C* algebra B into a state on B(H) (Appendix B1).
"""
encode(F::ChannelFactorization) = F.encode

"""
    decode(F::ChannelFactorization) -> CPMap

The decode CHANNEL Dec = Δ* : B(H) → B (CPTP, `dim_in = N`, `dim_out = n_B`).
Takes a state on B(H) into a state on the genuine C* algebra B (Appendix B1).

TRACE-PRESERVATION IS ONLY O(η)-EXACT for η>0 inputs (FINDINGS §C14). The
underlying observable map Δ′ is only O(η)-completely-positive, so the internal
Choi→Kraus PSD-cone projection (which clips the small negative cone mass) leaves
`Σ K† K = I` satisfied only to ≈ O(η·machine) — empirically ~1e-6 on multi-block
oblique fixtures, NOT 1e-9. Therefore `iscptp(decode(F))` returns `false` at the
default `atol = 1e-9` for a realistic η>0 factorization; pass a loosened tolerance
(`iscptp(decode(F); atol = 1e-5)`) to confirm it is TP to the artifact's accuracy.
The RIGOROUS certificate of the round-trip is [`upsdel_defect`](@ref) /
[`delups_defect`](@ref) — the certified cb-norm brackets — not `iscptp` at machine
tolerance (Appendix B4). The η=0 oracle decode IS exactly TP (defect 0).
"""
decode(F::ChannelFactorization) = F.decode

# ============================================================================
# choi(Φ) — the Convention-A Choi of Φ ITSELF (not Φ²−Φ; choi_diff is that)
# ============================================================================

"""
    choi(Φ::UCPMap) -> Matrix{ComplexF64}

The Convention-A Choi matrix of Φ ITSELF (n²×n²),
`C[(i-1)n+a, (j-1)n+b] = Σ_x conj(K_x[i,a]) K_x[j,b]` — conjugation on the FIRST
factor (aic_ucp.h:36). NOTE: this is the Choi of Φ, not of Φ²−Φ; the latter is
the low-level [`choi_diff`](@ref). Pure Julia (LinearAlgebra only).
"""
function choi(Φ::UCPMap)::Matrix{ComplexF64}
    nn = Φ.n
    C = zeros(ComplexF64, nn * nn, nn * nn)
    for K in Φ.kraus, i in 1:nn, a in 1:nn, j in 1:nn, b in 1:nn
        C[(i - 1) * nn + a, (j - 1) * nn + b] += conj(K[i, a]) * K[j, b]
    end
    return C
end
