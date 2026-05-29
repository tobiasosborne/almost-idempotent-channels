# runtests.jl — tests for AlmostIdempotentChannels.jl (bead aic-m24, increment 2b).
#
# Cross-check ladder (CLAUDE.md Rule 6) for the Julia entry point:
#   T1 ccall round-trip:  choi_diff(kraus) (C/arb via ccall) == pure-Julia
#                         Choi(Phi^2-Phi) to ~1e-10 — checks the marshalling.
#   T2 analytic anchors:  eta_idempotence == the EXACT closed forms (id->0;
#                         Phi_t -> t(1-t)*2(1-1/d^2); paper -> eta*sqrt(1-eta)).
#   T3 certified bracket: lo <= eta_idempotence <= hi (eig-free C bracket brackets
#                         the MOSEK value across the ccall boundary).
#   T4 GADC (live, non-fixture): eta~0 at gamma=0,1 (idempotent endpoints),
#                         interior > 0.
#
# Channel constructors + pure-Julia Choi are copied from tools/gen_fixtures_d24.jl;
# GADC Kraus from tools/smoke_gadc_cbnorm.jl. NO PYTHON. Single Julia process.

using AlmostIdempotentChannels
using Test
using LinearAlgebra

# ---- channel constructors (verbatim from tools/gen_fixtures_d24.jl) ----
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
identity_kraus(d) = [Matrix{ComplexF64}(I, d, d)]
depolarizing_kraus(d) = [ (e(d, i) * e(d, j)') / sqrt(d) for i in 1:d for j in 1:d ]
function phit_kraus(d, t)  # Phi_t = (1-t) id + t Dep_d
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end
# Paper example (.tex:367-378): Phi(X)=P0 Tr(g0 X)+P1 Tr(g1 X) on C^2.
function paper_example_kraus(etap)
    g0 = ComplexF64[1-etap sqrt(etap*(1-etap)); sqrt(etap*(1-etap)) etap]
    g1 = ComplexF64[0 0; 0 1]
    P0 = ComplexF64[1 0; 0 0]; P1 = ComplexF64[0 0; 0 1]
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

# ---- pure-Julia Choi (Convention A) + Heisenberg self-composition (verbatim) ----
function choi_A(Kraus::Vector{Matrix{ComplexF64}}, n::Int)
    C = zeros(ComplexF64, n * n, n * n)
    for K in Kraus, i in 1:n, a in 1:n, j in 1:n, b in 1:n
        C[(i - 1) * n + a, (j - 1) * n + b] += conj(K[i, a]) * K[j, b]
    end
    return C
end
compose_self(Kr::Vector{Matrix{ComplexF64}}) =
    [Kb * Ka for Ka in Kr for Kb in Kr]
lambda_choi(Kr::Vector{Matrix{ComplexF64}}, n::Int) =
    choi_A(compose_self(Kr), n) - choi_A(Kr, n)

# ---- GADC Kraus (verbatim from tools/smoke_gadc_cbnorm.jl) ----
function gadc_kraus(gamma::Float64, p::Float64)
    sp, sq = sqrt(p), sqrt(1 - p)
    sg, sg1 = sqrt(gamma), sqrt(1 - gamma)
    A0 = ComplexF64[sp 0; 0 sp*sg1]
    A1 = ComplexF64[0 sp*sg; 0 0]
    A2 = ComplexF64[sq*sg1 0; 0 sq]
    A3 = ComplexF64[0 0; sq*sg 0]
    return Matrix{ComplexF64}[A0, A1, A2, A3]
end

@testset "AlmostIdempotentChannels" begin

    @testset "T1 ccall round-trip: choi_diff vs pure-Julia Choi(Phi^2-Phi)" begin
        # Cross-checks the end-to-end ccall marshalling: the C/arb shim's
        # Convention-A Choi of Lambda=Phi^2-Phi must equal the pure-Julia one.
        worst = 0.0
        for (name, kr, n) in (("identity_2", identity_kraus(2), 2),
                              ("phi_t_2_0p3", phit_kraus(2, 0.3), 2))
            Jc = choi_diff(kr, n)             # C/arb via ccall
            Jj = lambda_choi(kr, n)           # pure Julia
            d = maximum(abs.(Jc .- Jj))
            worst = max(worst, d)
            @test d <= 1e-10
        end
        @info "T1 round-trip worst max-entry diff (C/arb ccall vs pure Julia)" worst
    end

    @testset "T2 eta_idempotence matches exact analytic anchors" begin
        dep_id_norm(d) = 2 * (1 - 1 / d^2)   # ||Dep_d - id_d||_diamond
        # identity -> 0
        @test eta_idempotence(identity_kraus(2)) ≈ 0.0 atol = 1e-8
        # Phi_t = (1-t) id + t Dep_d  =>  ||Phi_t^2-Phi_t||_diamond = t(1-t)*2(1-1/d^2)
        @test eta_idempotence(phit_kraus(2, 0.3)) ≈ 0.3 * 0.7 * dep_id_norm(2) atol = 1e-6  # = 0.315
        @test eta_idempotence(phit_kraus(2, 0.1)) ≈ 0.1 * 0.9 * dep_id_norm(2) atol = 1e-6  # = 0.135
        # paper example (.tex:378): ||Phi^2-Phi||_cb = eta*sqrt(1-eta)
        @test eta_idempotence(paper_example_kraus(0.1)) ≈ 0.1 * sqrt(0.9) atol = 1e-6
        @info "T2 anchors" id = eta_idempotence(identity_kraus(2)) phi_t_0p3 = eta_idempotence(phit_kraus(2, 0.3)) phi_t_0p1 = eta_idempotence(phit_kraus(2, 0.1)) paper_0p1 = eta_idempotence(paper_example_kraus(0.1))
    end

    @testset "T3 certified eig-free bracket brackets the SDP value" begin
        # The certified C bracket [lo,hi] (eig-free, no solver) must contain the
        # MOSEK diamond-norm value — ties the certified C bracket to the SDP across
        # the ccall boundary. Slack 1e-7 for the SDP solver tolerance.
        margins = NamedTuple[]
        for (name, kr) in (("phi_t_2_0p3", phit_kraus(2, 0.3)),
                           ("phi_t_2_0p1", phit_kraus(2, 0.1)),
                           ("paper_0p1", paper_example_kraus(0.1)),
                           ("dep_2", depolarizing_kraus(2)),
                           ("phi_t_3_0p2", phit_kraus(3, 0.2)))
            lo, hi = eta_eigfree(kr)
            eta = eta_idempotence(kr)
            @test lo - 1e-7 <= eta <= hi + 1e-7
            push!(margins, (name = name, lo = lo, eta = eta, hi = hi,
                            below = eta - lo, above = hi - eta))
        end
        @info "T3 bracket margins (lo <= eta <= hi)" margins
    end

    @testset "T4 GADC (live channel): eta=0 at gamma=0,1; interior > 0" begin
        # eta=0 anchors: gamma=0 -> identity; gamma=1 -> reset-to-fixed-state
        # (T^2=T). Interior 0<gamma<1 is not idempotent -> eta > 0.
        p = 0.5
        eta0 = eta_idempotence(gadc_kraus(0.0, p))
        eta1 = eta_idempotence(gadc_kraus(1.0, p))
        etamid = eta_idempotence(gadc_kraus(0.5, p))
        @test eta0 ≈ 0.0 atol = 1e-7
        @test eta1 ≈ 0.0 atol = 1e-7
        @test etamid > 1e-3
        @info "T4 GADC (p=0.5)" eta_gamma0 = eta0 eta_gamma1 = eta1 eta_gamma0p5 = etamid
    end

end
