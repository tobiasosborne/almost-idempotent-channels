# runtests.jl — the top-level test entry point (bead aic-exa.9 [T], design §5/§9).
#
# THE OPTIONAL-MOSEK DESIGN. The CORE + factorize + Aqua testsets ALWAYS run and
# are SOLVER-FREE — a user with no solver gets a green, useful package (the headline
# pipeline certified_defect / associated_algebra / main_isomorphism / factorize is
# entirely solver-free; only the exact cb-norm VALUE needs MOSEK). The SDP testset
# (the exact-value cross-checks) runs ONLY `if HAVE_MOSEK`, with @test_skip markers
# otherwise so the count stays honest.
#
# HAVE_MOSEK: an extension activates only when its trigger packages are LOADED in
# the session (not merely present in the env), so we `using` them here. They are
# test-target deps (Project.toml [targets].test); if the load fails (no solver
# installed / no license) the SDP testset is SKIPPED. The override env var
# AIC_FORCE_NO_MOSEK=1 forces HAVE_MOSEK=false for the GATE-CORE no-solver gate.
#
# Files: fixtures.jl (shared channel constructors), test_core.jl (solver-free core),
# test_factorize.jl (solver-free round-trip), test_sdp.jl (MOSEK-only), test_aqua.jl.

using AlmostIdempotentChannels
using Test
using LinearAlgebra
using Random
using Aqua

include("fixtures.jl")

const _FORCE_NO_MOSEK = get(ENV, "AIC_FORCE_NO_MOSEK", "0") == "1"
const HAVE_MOSEK = if _FORCE_NO_MOSEK
    @info "AIC_FORCE_NO_MOSEK=1 — forcing HAVE_MOSEK=false (GATE-CORE no-solver gate)"
    false
else
    try
        @eval using Convex, Mosek, MosekTools
        Base.get_extension(AlmostIdempotentChannels, :AICMosekExt) !== nothing
    catch err
        @warn "MOSEK extension unavailable; SDP suite SKIPPED" err
        false
    end
end

@testset "AlmostIdempotentChannels" begin
    @testset "core (solver-free)" begin
        include("test_core.jl")
    end

    @testset "factorize (solver-free)" begin
        include("test_factorize.jl")
    end

    @testset "Aqua" begin
        include("test_aqua.jl")
    end

    @testset "SDP (MOSEK)" begin
        if HAVE_MOSEK
            include("test_sdp.jl")
        else
            @info "MOSEK absent — SDP suite skipped (exact-value cross-checks)"
            @test_skip "idempotency_defect analytic anchors (MOSEK)"
            @test_skip "certified_defect(tight=true) contains value (MOSEK)"
            @test_skip "T5 strong duality (MOSEK)"
        end
    end
end
