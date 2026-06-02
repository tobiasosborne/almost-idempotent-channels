# test_core.jl — the SOLVER-FREE core suite (bead aic-exa.9 [T]). Green with NO
# MOSEK: the whole point of the optional-extension design (julia_package_design.md
# §3). Every assertion is solver-free; it touches the C/arb cores via ccall and the
# pure-Julia value-types, never the Watrous SDP. Cross-check ladder (CLAUDE.md
# Rule 6): the ccall round-trip (C/arb vs pure-Julia Choi), the certified eig-free
# bracket CONTAINMENT of the analytic η (with the §B9 slack), the η=0 oracle
# straddle, the main_isomorphism block cross-check against the C-test oracle, and
# the basin guards throwing cleanly IN-PROCESS (no SIGABRT).
#
# Assumes `using AlmostIdempotentChannels, Test, LinearAlgebra` + fixtures.jl
# already included by runtests.jl.

# ---- pure-Julia Choi (Convention A) + Heisenberg self-composition ----
function choi_A(Kraus::Vector{Matrix{ComplexF64}}, nn::Int)
    C = zeros(ComplexF64, nn * nn, nn * nn)
    for K in Kraus, i in 1:nn, a in 1:nn, j in 1:nn, b in 1:nn
        C[(i - 1) * nn + a, (j - 1) * nn + b] += conj(K[i, a]) * K[j, b]
    end
    return C
end
compose_self(Kr::Vector{Matrix{ComplexF64}}) = [Kb * Ka for Ka in Kr for Kb in Kr]
lambda_choi(Kr::Vector{Matrix{ComplexF64}}, nn::Int) =
    choi_A(compose_self(Kr), nn) - choi_A(Kr, nn)

dep_id_norm(d) = 2 * (1 - 1 / d^2)   # ‖Dep_d − id_d‖_⋄

@testset "T1 ccall round-trip: choi_diff (C/arb) vs pure-Julia Choi(Φ²−Φ)" begin
    # Cross-checks the end-to-end ccall marshalling: the C/arb shim's Convention-A
    # Choi of Λ = Φ²−Φ must equal the pure-Julia one to ~1e-10 (the marshalling is
    # byte-exact; the residual is the arb→double rounding).
    #
    # The first four fixtures are REAL-valued (max|imag(K)|=0), so for them the
    # conj() in the Convention-A Choi is a no-op — they pin only the K-major index
    # layout / transpose / scale, NOT the conjugation convention. complex_qubit is
    # the GENUINELY COMPLEX fixture (max|imag(K)|=0.8, max|imag(Choi(Φ²−Φ))|≈0.61):
    # it is the ONLY T1 case that catches a conj-on-the-wrong-factor marshalling bug
    # in the choi_diff shim. (NIT 1, aic-exa.9 hostile review.) It is η>0 but
    # choi_diff only forms Choi(Φ²−Φ) — no basin guard — so any UCP map is admissible.
    worst = 0.0
    cmplx_seen = 0.0
    for (name, kr, nn) in (("identity_2", identity_kraus(2), 2),
                           ("phi_t_2_0p3", phit_kraus(2, 0.3), 2),
                           ("paper_0p1", paper_example_kraus(0.1), 2),
                           ("dephasing_3", dephasing_kraus(3), 3),
                           ("complex_qubit", complex_qubit_kraus(), 2))
        Jc = choi_diff(kr, nn)        # C/arb via ccall
        Jj = lambda_choi(kr, nn)      # pure Julia
        d = maximum(abs.(Jc .- Jj))
        worst = max(worst, d)
        @test d <= 1e-10
        if name == "complex_qubit"
            # guard the guard: confirm this fixture actually exercises conj (a
            # real-valued Jj would silently make conj a no-op and re-blind T1).
            cmplx_seen = maximum(abs.(imag.(Jj)))
            @test cmplx_seen > 1e-2
        end
    end
    @info "T1 round-trip worst max-entry diff (C/arb ccall vs pure Julia)" worst complex_qubit_max_imag_choi = cmplx_seen
end

@testset "T6 UCPMap construction + validation" begin
    Φ = UCPMap(identity_kraus(2))
    @test n(Φ) == 2
    @test nkraus(Φ) == 1
    @test isunital(Φ)
    @test length(kraus(Φ)) == 1
    # non-unital → ArgumentError naming the measured defect
    nonunital = Matrix{ComplexF64}[ComplexF64[2 0; 0 1]]   # Σ K†K = diag(4,1) ≠ I
    err = try
        UCPMap(nonunital)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("not unital", err.msg)
    @test occursin("Σ K†K", err.msg) || occursin("K†K", err.msg)
    # ragged / non-square → DimensionMismatch
    @test_throws DimensionMismatch UCPMap(Matrix{ComplexF64}[ComplexF64[1 0; 0 1],
                                                              ComplexF64[1 0 0; 0 1 0]])
    @test_throws DimensionMismatch UCPMap(Matrix{ComplexF64}[ComplexF64[1 0 0; 0 1 0]])
    # empty → ArgumentError
    @test_throws ArgumentError UCPMap(Matrix{ComplexF64}[])
    # Choi round-trip constructor: build a UCPMap from the Convention-A Choi of a
    # known channel, recover a UCP self-map that is unital + same dim.
    Φb = block_cond_exp_kraus(4, 2)
    C = choi_A(Φb, 4)
    Φr = UCPMap(; choi=C, n=4)
    @test n(Φr) == 4
    @test isunital(Φr; atol=1e-9)
end

@testset "T8 CertifiedBracket: in / width / midpoint / value + validating ctor" begin
    b = CertifiedBracket(0.1, 0.3; label="x")
    @test 0.2 in b
    @test !(0.5 in b)
    @test minimum(b) == 0.1
    @test maximum(b) == 0.3
    @test width(b) ≈ 0.2
    @test midpoint(b) ≈ 0.2
    @test value(b) === nothing
    bv = CertifiedBracket(0.1, 0.3; value=0.2, label="x")
    @test value(bv) == 0.2
    # validating ctor throws: lo > hi
    @test_throws ArgumentError CertifiedBracket(0.5, 0.1)
    # validating ctor throws: value outside [lo-1e-9, hi+1e-9]
    @test_throws ArgumentError CertifiedBracket(0.1, 0.3; value=0.9)
end

@testset "T2/T3 certified_defect CONTAINMENT of the analytic η (§B9 slack)" begin
    # The eig-free bracket [lo,hi] must CONTAIN the analytic η. We assert with the
    # §B9 / aic-exa.12 N-2 SLACK: minimum(b)-tol ≤ η ≤ maximum(b)+tol — NOT strict
    # Base.in — because the arb FLOOR-rounded lo equals the analytic η to ~15 digits
    # and can sit a hair ABOVE it at near-η lower bounds. NO solver: the value comes
    # from the closed form, the bracket from aic_cbnorm_eigfree_d.
    tol = 1e-7
    margins = NamedTuple[]
    cases = (("phi_t_2_0p3", phit_kraus(2, 0.3), 0.3 * 0.7 * dep_id_norm(2)),
             ("phi_t_2_0p1", phit_kraus(2, 0.1), 0.1 * 0.9 * dep_id_norm(2)),
             ("phi_t_3_0p05", phit_kraus(3, 0.05), 0.05 * 0.95 * dep_id_norm(3)),
             ("paper_0p1", paper_example_kraus(0.1), 0.1 * sqrt(0.9)))
    for (name, kr, η) in cases
        b = certified_defect(UCPMap(kr))
        @test value(b) === nothing                       # solver-free: no point estimate
        @test minimum(b) - tol <= η <= maximum(b) + tol   # §B9 slack containment
        push!(margins, (name=name, lo=minimum(b), η=η, hi=maximum(b)))
    end
    @info "T2/T3 certified_defect brackets (lo ≤ η ≤ hi, §B9 slack)" margins
end

@testset "associated_algebra: dim_A + ε bracket straddles 0 at η=0" begin
    A = associated_algebra(UCPMap(identity_kraus(2)))
    @test ambient(A) == 2
    @test dim(A) == 4                       # identity on C²: dim_A = n² = 4
    @test minimum(epsilon(A)) <= 0.0 <= maximum(epsilon(A))   # straddle 0
    @test maximum(epsilon(A)) < 1e-9        # η=0 → ε ~ machine
    Abce = associated_algebra(UCPMap(block_cond_exp_kraus(4, 2)))
    @test minimum(epsilon(Abce)) <= 0.0 <= maximum(epsilon(Abce))
    @test maximum(epsilon(Abce)) < 1e-9
end

@testset "main_isomorphism: C-oracle block lists + cbfwd/cbinv contain 1" begin
    # Cross-check the B-block sizes against the test_shim_assoc_iso.c oracle:
    #   identity(2) → [2], block_cond_exp(4,2) → [2,2], dephasing(3) → [1,1,1].
    for (name, kr, want) in (("identity(2)", identity_kraus(2), [2]),
                             ("block_cond_exp(4,2)", block_cond_exp_kraus(4, 2), [2, 2]),
                             ("dephasing(3)", dephasing_kraus(3), [1, 1, 1]))
        v = main_isomorphism(UCPMap(kr))
        @test sort(blocks(v)) == sort(want)
        # complete-isometry (η=0): ‖v‖_cb = ‖v⁻¹‖_cb = 1 lies in the brackets.
        @test 1.0 in cbnorm_forward(v)
        @test 1.0 in cbnorm_inverse(v)
        @test maximum(isodefect(v)) < 1e-9   # η=0 iso defect ~ machine
    end
end

@testset "choi / kraus / isunital" begin
    Φ = UCPMap(block_cond_exp_kraus(4, 2))
    C = choi(Φ)
    @test size(C) == (16, 16)
    @test C ≈ C' atol = 1e-12                # Choi of a CP map is Hermitian
    @test isunital(Φ)
    @test length(kraus(Φ)) == nkraus(Φ)
end

@testset "BASIN GUARDS throw cleanly IN-PROCESS (no SIGABRT)" begin
    # factorize(phit(2,0.3)): in the prop_P basin but OUT of factorize's tighter
    # xpow domain (ρ_est ≈ 0.21 ≥ 0.10) → ArgumentError, NOT a C abort that kills
    # the session. The test running past this line IS the no-SIGABRT proof.
    err1 = try
        factorize(UCPMap(phit_kraus(2, 0.3)))
        nothing
    catch e
        e
    end
    @test err1 isa ArgumentError
    @test occursin("convergence domain", err1.msg)
    # factorize(UCPMap([diag(1,-1)])): reflection conjugation, ρ(Φ²−Φ)=2, far OUT
    # of the prop_P basin → ArgumentError (the prop_P / xpow guard), in-process.
    refl = UCPMap([ComplexF64[1 0; 0 -1]])
    err2 = try
        factorize(refl)
        nothing
    catch e
        e
    end
    @test err2 isa ArgumentError
    # associated_algebra on the reflection hits the wider prop_P guard.
    err3 = try
        associated_algebra(refl)
        nothing
    catch e
        e
    end
    @test err3 isa ArgumentError
    @test occursin("prop_P basin", err3.msg)
    # certified_defect stays SAFE on both (eig-free, never aborts) — the always-safe
    # escape hatch the guard messages point at.
    @test certified_defect(UCPMap(phit_kraus(2, 0.3))) isa CertifiedBracket
    bref = certified_defect(refl)
    @test bref isa CertifiedBracket
    @test minimum(bref) > 1.0               # ρ=2 reflection: η is large, bracket > 1
    @test isunital(refl)                    # diag(1,-1) is unitary ⇒ unital
end
