# fixtures.jl — shared channel constructors for the test suite (bead aic-exa.9 [T]).
#
# These build the UCP self-maps the cross-checks exercise, as Kraus lists in the
# HEISENBERG picture Φ(X) = Σ_a K_a† X K_a. They are pure Julia (LinearAlgebra
# only), deterministic (no RNG), and reused across test_core.jl / test_factorize.jl
# / test_sdp.jl. NO PYTHON. The η=0 oracles (identity, block_cond_exp, dephasing)
# are EXACTLY idempotent; the η>0 obliques (mixconj, bce_conj) are the §C2/§C13
# fixtures the η=0 oracle is structurally BLIND to (the design's mandatory
# oblique fixtures, julia_package_design.md §C2 / Appendix B1/B4).
#
# The oblique fixtures mirror the C test_idemp.h builders byte-for-byte in
# construction (make_compress_idemp + make_fixed_unitary + make_conjugated +
# make_mixconj; make_block_cond_exp), so the Julia cross-checks line up with the
# C-test oracle values in tests/test_shim_assoc_iso.c / test_shim_factorize_artifacts.c.

using LinearAlgebra

# ---- basis vector ----
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)

# ============================================================================
# η = 0 EXACT-IDEMPOTENT oracles (the cleanest ground truth, CLAUDE.md Rule 6)
# ============================================================================

# Identity Φ = id on C^d (test_idemp.h make_identity). dim_A = d², B-blocks = [d].
identity_kraus(d) = [Matrix{ComplexF64}(I, d, d)]

# Block conditional expectation on C^d, first m vs the rest (test_idemp.h
# make_block_cond_exp): K_0 = diag(1..1,0..0) (m ones), K_1 = the complement.
# EXACTLY idempotent; multi-block B = M_m ⊕ M_{d-m} (blocks [m, d-m]).
function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d)
    K1 = zeros(ComplexF64, d, d)
    for i in 1:m
        K0[i, i] = 1
    end
    for i in (m + 1):d
        K1[i, i] = 1
    end
    return Matrix{ComplexF64}[K0, K1]
end

# Dephasing Φ(X) = diag(X) on C^d (test_idemp.h make_dephasing): K_i = |e_i><e_i|.
# EXACTLY idempotent; commutative B = ⊕ M_1 (blocks [1,…,1], d ones).
function dephasing_kraus(d)
    K = Matrix{ComplexF64}[]
    for i in 1:d
        Ki = zeros(ComplexF64, d, d)
        Ki[i, i] = 1
        push!(K, Ki)
    end
    return K
end

# ============================================================================
# Depolarizing / Φ_t family + the paper example (analytic-η anchors)
# ============================================================================

depolarizing_kraus(d) = [(e(d, i) * e(d, j)') / sqrt(d) for i in 1:d for j in 1:d]

# Φ_t = (1-t) id + t Dep_d.  ‖Φ_t²−Φ_t‖_⋄ = t(1-t)·2(1-1/d²) (closed form).
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

# Paper example (.tex:367-378): Φ(X) = P0 Tr(g0 X) + P1 Tr(g1 X) on C².
# ‖Φ²−Φ‖_cb = η√(1-η).
function paper_example_kraus(etap)
    g0 = ComplexF64[1-etap sqrt(etap * (1 - etap)); sqrt(etap * (1 - etap)) etap]
    g1 = ComplexF64[0 0; 0 1]
    P0 = ComplexF64[1 0; 0 0]
    P1 = ComplexF64[0 0; 0 1]
    function prep(P, g)
        j = findfirst(c -> norm(P[:, c]) > 0, 1:2)
        u = P[:, j] / norm(P[:, j])
        ev = eigen(Hermitian(g))
        Ks = Matrix{ComplexF64}[]
        for k in 1:2
            ev.values[k] <= 1e-14 && continue
            push!(Ks, sqrt(ev.values[k]) * (ev.vectors[:, k] * u'))
        end
        return Ks
    end
    return vcat(prep(P0, g0), prep(P1, g1))
end

# Complex asymmetric unital qubit channel (the T5 discriminating anchor).
complex_qubit_kraus() = [ComplexF64[0.6 0; 0.8im 0], ComplexF64[0 0.8im; 0 0.6]]

# GADC Kraus (live, non-fixture; η=0 at γ=0,1, interior > 0).
function gadc_kraus(gamma::Float64, p::Float64)
    sp, sq = sqrt(p), sqrt(1 - p)
    sg, sg1 = sqrt(gamma), sqrt(1 - gamma)
    A0 = ComplexF64[sp 0; 0 sp*sg1]
    A1 = ComplexF64[0 sp*sg; 0 0]
    A2 = ComplexF64[sq*sg1 0; 0 sq]
    A3 = ComplexF64[0 0; sq*sg 0]
    return Matrix{ComplexF64}[A0, A1, A2, A3]
end

# ============================================================================
# η > 0 OBLIQUE fixtures (the §C2/§C13 fixtures; mirror C test_idemp.h)
# ============================================================================

# make_compress_idemp(d,m): exactly-idempotent, proper carrier M = span(e0..e_{m-1}).
# K_0 = Π_M, K_{1+b} = |e_0><e_{m+b}| (test_idemp.h:123).
function compress_idemp_kraus(d, m)
    K = Matrix{ComplexF64}[zeros(ComplexF64, d, d)]
    for i in 1:m
        K[1][i, i] = 1
    end
    for b in 0:(d - m - 1)
        Kb = zeros(ComplexF64, d, d)
        Kb[1, m + b + 1] = 1
        push!(K, Kb)
    end
    return K
end

# make_fixed_unitary(n): V = exp(iH) for a fixed Hermitian H (real super-diagonal
# 0.4, complex 0.3i second-super-diagonal, diagonal phases 0.1·k) — non-real,
# basis-mixing (test_idemp.h:139). Used to conjugate an idempotent into an OBLIQUE
# image (the §C2/§C13 tooth the η=0 identity oracle is blind to).
function fixed_unitary(nn)
    H = zeros(ComplexF64, nn, nn)
    for p in 1:(nn - 1)
        H[p, p + 1] = 0.4
        H[p + 1, p] = 0.4
    end
    for p in 1:(nn - 2)
        H[p, p + 2] = 0.3im
        H[p + 2, p] = -0.3im
    end
    for p in 1:nn
        H[p, p] = 0.1 * p
    end
    return exp(im * Matrix(Hermitian(H)))
end

# Convex mix of an exactly-idempotent `base` (on C^d) with its unitary conjugate
# {V† K_a V}: Φ_t = (1-t) Φ_base + t Ad_{V†}∘Φ_base∘Ad_V, via the Kraus union
# {√(1-t) K_a} ∪ {√t V† K_a V} (test_idemp.h make_mixconj / make_conjugated). A
# GENUINELY η>0, oblique, non-commutative almost-idempotent channel.
function _mix_conj_of(base::Vector{Matrix{ComplexF64}}, d, t)
    V = fixed_unitary(d)
    Vd = V'
    K = Matrix{ComplexF64}[]
    for Ka in base
        push!(K, sqrt(1 - t) * Ka)
    end
    for Ka in base
        push!(K, sqrt(t) * (Vd * Ka * V))
    end
    return K
end

# mixconj(d,m,t): the C make_mixconj — compress_idemp(d,m) conjugate-mixed. A
# SINGLE-BLOCK oblique η>0 fixture (B = M_m). In factorize's small-η domain for
# small t (ρ(Φ²−Φ) ≈ 0.04 at t=0.02). The PSD-cone clip is empty for the single
# block, so its decode is TP to machine precision.
mixconj_kraus(d, m, t) = _mix_conj_of(compress_idemp_kraus(d, m), d, t)

# bce_conj(d,m,t): block_cond_exp(d,m) conjugate-mixed — a MULTI-BLOCK (B = M_m ⊕
# M_{d-m}) oblique η>0 fixture. The genuine multi-block oblique structure triggers
# the §C14 Choi→Kraus PSD-cone clip in upsilon_build, so the decode channel is only
# O(η)-trace-preserving (measured TP-defect ≈ 1.9e-6 at t=0.02, ≈ 2.2e-5 at t=0.05):
# iscptp(decode; atol=1e-9) is FALSE, iscptp(decode; atol=1e-4) is TRUE. This is the
# Appendix-B4 fixture: the cb-norm round-trip bracket (upsdel/delups), NOT iscptp at
# machine tolerance, is the rigorous certificate. In factorize's domain for small t.
bce_conj_kraus(d, m, t) = _mix_conj_of(block_cond_exp_kraus(d, m), d, t)

# Trace-preservation defect ‖Σ_a K_a† K_a − I‖ of a CPMap's Kraus (for the B4
# O(η)-TP probe; CPMap has no public `kraus` accessor — read the field).
_tp_defect(c::CPMap) = opnorm(sum(K' * K for K in c.kraus) - I)
