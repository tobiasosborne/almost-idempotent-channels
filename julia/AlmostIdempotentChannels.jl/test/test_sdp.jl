# test_sdp.jl — the MOSEK-ONLY exact-value suite (bead aic-exa.9 [T], design §9
# step 1). Runs only under HAVE_MOSEK (runtests.jl includes it conditionally), so a
# user with no solver still gets a fully GREEN package (the core + factorize + aqua
# suites). Three rungs of the cross-check ladder that need the SDP:
#   - T5 strong duality: the MAX-primal value == the MIN-dual value == the exact
#     closed-form η for every analytic anchor (pins the QETLAB dual normalization;
#     the ASYMMETRIC paper example is the discriminating anchor).
#   - idempotency_defect matches the analytic anchors (the public exact-value verb).
#   - certified_defect(Φ; tight=true) CONTAINS the SDP value and is TIGHTER than the
#     solver-free eig-free bracket (ties the rigorous bracket to the value).
#
# Assumes `using AlmostIdempotentChannels, Test, LinearAlgebra`, fixtures.jl, AND
# `using Convex, Mosek, MosekTools` (so the AICMosekExt extension is active).

# pure-Julia Choi (Convention A) for the T5 raw-Choi anchors.
function _choi_A(Kraus::Vector{Matrix{ComplexF64}}, nn::Int)
    C = zeros(ComplexF64, nn * nn, nn * nn)
    for K in Kraus, i in 1:nn, a in 1:nn, j in 1:nn, b in 1:nn
        C[(i - 1) * nn + a, (j - 1) * nn + b] += conj(K[i, a]) * K[j, b]
    end
    return C
end
_compose_self(Kr::Vector{Matrix{ComplexF64}}) = [Kb * Ka for Ka in Kr for Kb in Kr]
_lambda_choi(Kr::Vector{Matrix{ComplexF64}}, nn::Int) =
    _choi_A(_compose_self(Kr), nn) - _choi_A(Kr, nn)

_dep_id_norm(d) = 2 * (1 - 1 / d^2)

@testset "idempotency_defect matches the exact analytic anchors" begin
    @test idempotency_defect(identity_kraus(2)) ≈ 0.0 atol = 1e-7
    @test idempotency_defect(phit_kraus(2, 0.3)) ≈ 0.3 * 0.7 * _dep_id_norm(2) atol = 1e-6
    @test idempotency_defect(phit_kraus(2, 0.1)) ≈ 0.1 * 0.9 * _dep_id_norm(2) atol = 1e-6
    @test idempotency_defect(paper_example_kraus(0.1)) ≈ 0.1 * sqrt(0.9) atol = 1e-6
    # GADC live anchors: η=0 at γ=0,1; interior > 0.
    @test idempotency_defect(gadc_kraus(0.0, 0.5)) ≈ 0.0 atol = 1e-7
    @test idempotency_defect(gadc_kraus(1.0, 0.5)) ≈ 0.0 atol = 1e-7
    @test idempotency_defect(gadc_kraus(0.5, 0.5)) > 1e-3
end

@testset "certified_defect(tight=true) contains the value + tighter than eig-free" begin
    # The MOSEK-tight rigorous bracket must (i) contain the exact SDP value and
    # (ii) be strictly TIGHTER than the solver-free eig-free bracket (width hi-lo).
    for (name, kr) in (("phi_t_2_0p3", phit_kraus(2, 0.3)),
                       ("phi_t_2_0p1", phit_kraus(2, 0.1)),
                       ("paper_0p1", paper_example_kraus(0.1)))
        Φ = UCPMap(kr)
        loose = certified_defect(Φ)                 # eig-free, solver-free
        tight = certified_defect(Φ; tight=true)     # MOSEK feasible points
        v = idempotency_defect(Φ)
        @test minimum(tight) - 1e-7 <= v <= maximum(tight) + 1e-7
        @test width(tight) < width(loose)
        @test value(tight) !== nothing               # carries the SDP point estimate
        @info "tight vs loose bracket" name v width_tight=width(tight) width_loose=width(loose)
    end
end

@testset "T5 strong duality: η_primal == η_dual == η_ref (pins the dual norm)" begin
    cases = Tuple{String,Matrix{ComplexF64},Int,Float64}[
        ("id_2",          _lambda_choi(identity_kraus(2), 2),                                  2, 0.0),
        ("dep2-id2",      _choi_A(depolarizing_kraus(2), 2) - _choi_A(identity_kraus(2), 2),  2, _dep_id_norm(2)),
        ("dep3-id3",      _choi_A(depolarizing_kraus(3), 3) - _choi_A(identity_kraus(3), 3),  3, _dep_id_norm(3)),
        ("dep4-id4",      _choi_A(depolarizing_kraus(4), 4) - _choi_A(identity_kraus(4), 4),  4, _dep_id_norm(4)),
        ("phi_t2_0p3",    _lambda_choi(phit_kraus(2, 0.3), 2),                                 2, 0.3 * 0.7 * _dep_id_norm(2)),
        ("phi_t3_0p05",   _lambda_choi(phit_kraus(3, 0.05), 3),                                3, 0.05 * 0.95 * _dep_id_norm(3)),
        ("phi_t4_0p2",    _lambda_choi(phit_kraus(4, 0.2), 4),                                 4, 0.2 * 0.8 * _dep_id_norm(4)),
        ("paper_0p1",     _lambda_choi(paper_example_kraus(0.1), 2),                           2, 0.1 * sqrt(0.9)),
        ("paper_0p01",    _lambda_choi(paper_example_kraus(0.01), 2),                          2, 0.01 * sqrt(0.99)),
        ("complex_qubit", _lambda_choi(complex_qubit_kraus(), 2),                              2, 1.28),
    ]
    rows = NamedTuple[]
    maxgap = 0.0
    for (name, J, nn, eref) in cases
        ep, _X, _P, _Q, ps = diamond_norm_watrous_primal(J, nn)
        ed, _Y0, _Y1, ds = diamond_norm_dual(J, nn)
        gap = abs(ep - ed)
        maxgap = max(maxgap, gap)
        @test ps in ("OPTIMAL", "SLOW_PROGRESS")
        @test ds in ("OPTIMAL", "SLOW_PROGRESS")
        @test isapprox(ep, eref; atol=1e-6)
        @test isapprox(ed, eref; atol=1e-6)    # pins the dual normalization
        @test isapprox(ep, ed; atol=1e-6)
        push!(rows, (name=name, eta_primal=ep, eta_dual=ed, eta_ref=eref, gap=gap))
    end
    @info "T5 strong-duality table (η_primal, η_dual, η_ref, gap)" maxgap
    println("\nT5 strong-duality / dual-normalization table:")
    println(rpad("channel", 16), rpad("eta_primal", 16), rpad("eta_dual", 16),
            rpad("eta_ref", 16), "gap")
    for r in rows
        println(rpad(r.name, 16),
                rpad(string(round(r.eta_primal; sigdigits=9)), 16),
                rpad(string(round(r.eta_dual; sigdigits=9)), 16),
                rpad(string(round(r.eta_ref; sigdigits=9)), 16),
                string(round(r.gap; sigdigits=3)))
    end
    println("max primal-dual gap = ", round(maxgap; sigdigits=3), "\n")
end
