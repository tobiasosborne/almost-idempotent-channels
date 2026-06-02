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

# ----- ChannelFactorization (the showcase) -----

# compact: "ChannelFactorization(Φ ≈ Δ Υ, B = ⊕ M_d [1, 1])"
Base.show(io::IO, F::ChannelFactorization) =
    print(io, "ChannelFactorization(Φ ≈ Δ Υ, B = ⊕ M_d ", F.B.blocks, ")")

function Base.show(io::IO, ::MIME"text/plain", F::ChannelFactorization)
    println(io, "ChannelFactorization  Φ ≈ Δ Υ   through  B = ⊕ M_d, blocks = ", F.B.blocks)
    # right-justify the two Kraus counts so the columns line up (§2.6).
    rΔ, rΥ = length(F.Δ), length(F.Υ)
    w = max(ndigits(rΔ), ndigits(rΥ))
    println(io, "  encode  Δ: B → B(H)   (", lpad(rΔ, w), " Kraus,  Dec = Δ*)")
    println(io, "  decode  Υ: B(H) → B   (", lpad(rΥ, w), " Kraus,  Enc = Υ*)")
    println(io, "  ‖ΔΥ − Φ‖_cb   ∈ ", sprint(show, F.delups), "     (encode∘decode ≈ Φ)")
    println(io, "  ‖ΥΔ − 1_B‖_cb ∈ ", sprint(show, F.upsdel), "     (decode∘encode ≈ 1_B)")
    print(io,   "  η proxy = ", _fmtnum(F.eta_proxy), ",  ε used = ", _fmtnum(F.eps_used))
end
