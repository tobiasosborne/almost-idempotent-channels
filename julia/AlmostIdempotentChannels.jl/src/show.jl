# show.jl — compact `show(io, x)` (one line, used in arrays / interpolation) and
# multi-line `show(io, ::MIME"text/plain", x)` (the REPL display) for every value
# type. The renderings MATCH the design doc §2.6 sketches; the ChannelFactorization
# and CertifiedBracket renderings are the showcase.
#
# Unicode is restricted to the glyphs the project already uses in comments and
# docs (Φ Υ Δ ⊕ ≤ ∈ ε); no emoji, no marketing (Rule 12). Numbers are concrete.

# ----- number formatting helpers -----

# Format a scalar for a bracket endpoint: small/large magnitudes in scientific
# notation, mid-range with a few decimals, exact 0 as "0.0". Keeps the §2.6 look
# (0.0840, 2.1e-12, 0.27).
function _fmtnum(x::Real)
    x == 0 && return "0.0"
    ax = abs(x)
    if ax < 1e-3 || ax ≥ 1e5
        return Printf.@sprintf("%.2e", x)
    elseif ax < 1
        return Printf.@sprintf("%.4f", x)            # 0.0840, 0.2700
    elseif ax < 1e3
        return Printf.@sprintf("%.2f", x)            # 1.00, 1.04, 13.00 → near-1 brackets keep 2 dp
    else
        return Printf.@sprintf("%.4g", x)
    end
end

# ----- CertifiedBracket -----

# compact: "[0.084, 0.27]"  or with a value "[0.084, 0.135, 0.27]"
function Base.show(io::IO, b::CertifiedBracket)
    if b.value === nothing
        print(io, "[", _fmtnum(b.lo), ", ", _fmtnum(b.hi), "]")
    else
        print(io, "[", _fmtnum(b.lo), ", ", _fmtnum(b.value), ", ", _fmtnum(b.hi), "]")
    end
end

function Base.show(io::IO, ::MIME"text/plain", b::CertifiedBracket)
    head = isempty(b.label) ? "CertifiedBracket" : "CertifiedBracket  $(b.label)"
    println(io, head)
    print(io, "  ", _fmtnum(b.lo), "  ≤  x  ≤  ", _fmtnum(b.hi),
          "        (width ", _fmtnum(width(b)), ")")
    if b.value !== nothing
        print(io, "\n  value = ", _fmtnum(b.value))
    end
end

# ----- UCPMap -----

# compact: "UCPMap(n=2, r=4)"
Base.show(io::IO, Φ::UCPMap) = print(io, "UCPMap(n=", Φ.n, ", r=", Φ.r, ")")

function Base.show(io::IO, ::MIME"text/plain", Φ::UCPMap)
    println(io, "UCPMap: Φ: B(H) → B(H),  dim H = ", Φ.n)
    println(io, "  Kraus operators: ", Φ.r, "   (Heisenberg: Φ(X) = Σ K_a† X K_a)")
    S = sum(K' * K for K in Φ.kraus)
    defect = opnorm(S - I)
    if defect ≤ 1e-9
        print(io, "  unital: yes   (‖Σ K†K − I‖ = ", _fmtnum(defect), ")")
    else
        print(io, "  unital: no    (‖Σ K†K − I‖ = ", _fmtnum(defect),
              "; trace-preserving dual?)")
    end
end

# ----- EpsCStarAlgebra -----

# compact: "EpsCStarAlgebra(dim A=2 ≤ M_2)"
Base.show(io::IO, A::EpsCStarAlgebra) =
    print(io, "EpsCStarAlgebra(dim A=", A.dim_A, " ≤ M_", A.n, ")")

function Base.show(io::IO, ::MIME"text/plain", A::EpsCStarAlgebra)
    println(io, "EpsCStarAlgebra  A = Img Φ̃  ≤  M_", A.n)
    println(io, "  dim A = ", A.dim_A)
    println(io, "  ε (associator defect) ∈ ", sprint(show, A.eps))
    print(io, "  product: Choi–Effros star  X⋆Y = Φ̃(XY)   (oblique; axioms hold up to ε)")
end

# ----- CStarAlgebra -----

# compact: "⊕ M_d [1, 1]"
Base.show(io::IO, B::CStarAlgebra) = print(io, "⊕ M_d ", B.blocks)

function Base.show(io::IO, ::MIME"text/plain", B::CStarAlgebra)
    println(io, "CStarAlgebra  B = ⊕_l M_{d_l}")
    println(io, "  blocks d_l = ", B.blocks, "   (m = ", length(B.blocks), ")")
    print(io, "  dim_B = ", B.dim_B, ",  n_B = ", B.n_B)
end

# ----- MainIsomorphism -----

# compact: "MainIsomorphism(B = ⊕ M_d [1, 1])"
Base.show(io::IO, v::MainIsomorphism) =
    print(io, "MainIsomorphism(B = ⊕ M_d ", v.B.blocks, ")")

function Base.show(io::IO, ::MIME"text/plain", v::MainIsomorphism)
    println(io, "MainIsomorphism  v: B → A   (O(ε), dimension-independent constant)")
    println(io, "  B = ⊕ M_d, blocks = ", v.B.blocks)
    println(io, "  ‖v‖_cb    ∈ ", sprint(show, v.cbnorm_fwd))
    println(io, "  ‖v⁻¹‖_cb  ∈ ", sprint(show, v.cbnorm_inv))
    print(io,   "  iso defect ∈ ", sprint(show, v.isodefect))
end

# ----- CPMap (rectangular CPTP, Appendix B7/B10) -----

# compact: "CPMap(2 → 4, 8 Kraus)"
Base.show(io::IO, c::CPMap) =
    print(io, "CPMap(", c.dim_in, " → ", c.dim_out, ", ", length(c.kraus), " Kraus)")

function Base.show(io::IO, ::MIME"text/plain", c::CPMap)
    println(io, "CPMap: CPTP  C^", c.dim_in, " → C^", c.dim_out,
            "   (Schrödinger: ρ ↦ Σ K_a ρ K_a†)")
    println(io, "  Kraus operators: ", length(c.kraus),
            "   (each ", c.dim_out, "×", c.dim_in, ")")
    defect = opnorm(sum(K' * K for K in c.kraus) - I)
    # An η>0 factorize `decode` channel is only O(η)-TP (the §C14 PSD-cone clip),
    # so a measured defect of ~1e-6 is EXPECTED, not a bug — say "approx (O(η))"
    # rather than a bald "no" so the show is honest about the certificate level.
    if defect ≤ 1e-9
        print(io, "  trace-preserving: yes        (‖Σ K_a† K_a − I_", c.dim_in, "‖ = ",
              _fmtnum(defect), ")")
    elseif defect ≤ 1e-4
        print(io, "  trace-preserving: approx O(η) (‖Σ K_a† K_a − I_", c.dim_in, "‖ = ",
              _fmtnum(defect), ";  the cb-norm round-trip bracket is the certificate)")
    else
        print(io, "  trace-preserving: no         (‖Σ K_a† K_a − I_", c.dim_in, "‖ = ",
              _fmtnum(defect), ")")
    end
end

# ----- ChannelFactorization (the showcase) -----
#
# UNAMBIGUOUS by construction (Appendix B7): the header carries the OBSERVABLE
# identity Φ ≈ Δ Υ; the body lists the two CHANNELS by their physical role and
# dual binding (encode = Υ* : B → B(H); decode = Δ* : B(H) → B). The word
# "encode" is NEVER placed next to "Dec=Δ*" (the §C13 conflation trap).

# compact: "ChannelFactorization(Φ ≈ Δ Υ, B = ⊕ M_d [1, 1])"
Base.show(io::IO, F::ChannelFactorization) =
    print(io, "ChannelFactorization(Φ ≈ Δ Υ, B = ⊕ M_d ", F.B.blocks, ")")

function Base.show(io::IO, ::MIME"text/plain", F::ChannelFactorization)
    enc, dec = F.encode, F.decode
    println(io, "ChannelFactorization  Φ ≈ Δ Υ  through  B = ⊕ M_d, blocks = ", F.B.blocks)
    # right-justify the two Kraus counts so the columns line up (§2.6).
    re, rd = length(enc.kraus), length(dec.kraus)
    w = max(ndigits(re), ndigits(rd))
    println(io, "  encode = Υ* : B → B(H)   (", lpad(re, w), " Kraus,  ",
            enc.dim_in, "→", enc.dim_out, " CPTP)")
    println(io, "  decode = Δ* : B(H) → B   (", lpad(rd, w), " Kraus,  ",
            dec.dim_in, "→", dec.dim_out, " CPTP)")
    println(io, "  ‖ΔΥ − Φ‖_cb   ∈ ", sprint(show, F.delups), "     (round-trip ΔΥ ≈ Φ)")
    println(io, "  ‖ΥΔ − 1_B‖_cb ∈ ", sprint(show, F.upsdel), "     (round-trip ΥΔ ≈ 1_B)")
    print(io,   "  η proxy = ", _fmtnum(F.eta_proxy), ",  ε used = ", _fmtnum(F.eps_used))
end
