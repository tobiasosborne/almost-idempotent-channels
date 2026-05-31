# probe_o2_diag2.jl — DECISIVE convention oracle for opspace O2 (bead aic-pjr).
# diag1 showed the rectangular SDP value depends on the DECLARED (d_maj,d_min)
# and that the *-isometry (||.||_⋄=1) is too symmetric to disambiguate. This
# probe uses an ASYMMETRIC CP map with an INDEPENDENT CLOSED-FORM diamond norm,
# so no convention coincidence can survive:
#
#   Psi: M_3 -> M_2,  Psi(Y) = A' Y A,  A a random 3x2 (NOT an isometry).
#   Psi is CP, so  ||Psi||_⋄ = max_{rho in D(M_3)} Tr[Psi(rho)]
#                            = max Tr[rho A A'] = lambda_max(A A') = sigma_max(A)^2.
#
# We sweep every feed/declared-dims/placement/direction combo and flag those whose
# RAW optval equals sigma_max(A)^2 (ratio 1) — that is the correct convention +
# (no) normalization. Then confirm it also gives 1.0 on the *-iso oracle.
#   julia --project=julia/env tools/probe_o2_diag2.jl    (serial; NO PYTHON)

using Convex, MosekTools, LinearAlgebra, Printf, Random
include(joinpath(@__DIR__, "..", "julia", "AlmostIdempotentChannels.jl",
                 "src", "sdp.jl"))

function choi_convA(fmap, a::Int, b::Int)      # J[s*b+i,t*b+j]=f(E_st)[i,j]
    J = zeros(ComplexF64, a * b, a * b)
    for s in 1:a, t in 1:a
        Est = zeros(ComplexF64, a, a); Est[s, t] = 1.0
        F = fmap(Est)
        for i in 1:b, j in 1:b
            J[(s - 1) * b + i, (t - 1) * b + j] = F[i, j]
        end
    end
    return J
end

Random.seed!(20260531)
println("="^78)
println("O2 DECISIVE oracle — asymmetric CP map Psi(Y)=A'YA, M_3->M_2")
println("="^78)

# A : 3x2 random (NOT isometry). Psi: M_3 -> M_2, input dim 3, output dim 2.
A = randn(ComplexF64, 3, 2)
psi(Y) = A' * Y * A                            # Y in M_3 -> M_2
truth = maximum(real.(eigvals(Hermitian(A * A')))) # sigma_max(A)^2 = ||Psi||_⋄
@printf("  A is 3x2 (rank %d); sigma_max(A)^2 = ||Psi||_⋄ = %.8f  (INDEPENDENT truth)\n",
        rank(A), truth)

# J(Psi) Convention A: source=input=3 (major), target=output=2 (minor). 6x6.
Jpsi = choi_convA(psi, 3, 2)
in_dim, out_dim = 3, 2

# Sweep: feed in {J, transpose(J)} x declared dims in {(3,2),(2,3)} x rho_on.
feeds = [("J", Jpsi), ("transpose(J)", permutedims(Jpsi))]
dimset = [(3, 2), (2, 3)]
println("\n  PRIMAL raw optvals (flag = matches truth $(round(truth,digits=6)) to 1e-5):")
@printf("    %-14s %-10s %-8s %-14s %s\n", "feed", "dims", "rho_on", "raw", "ratio raw/truth")
prim_hits = Tuple{String,Tuple{Int,Int},Symbol,Float64}[]
for (flbl, M) in feeds, dm in dimset, ro in (:major, :minor)
    try
        r = diamond_primal_rect(Matrix{ComplexF64}(M), dm[1], dm[2], ro)
        rat = r[1] / truth
        flag = abs(r[1] - truth) <= 1e-5 ? "  <== TRUTH" : ""
        abs(r[1] - truth) <= 1e-5 && push!(prim_hits, (flbl, dm, ro, r[1]))
        @printf("    %-14s %-10s %-8s %-14.8f %.6f%s\n", flbl, string(dm), String(ro), r[1], rat, flag)
    catch e
        @printf("    %-14s %-10s %-8s FAILED %s\n", flbl, string(dm), String(ro), sprint(showerror, e))
    end
end

println("\n  DUAL raw optvals (flag = matches truth to 1e-5):")
@printf("    %-14s %-10s %-8s %-14s %s\n", "feed", "dims", "tr_sys", "raw", "ratio raw/truth")
dual_hits = Tuple{String,Tuple{Int,Int},Int,Float64}[]
for (flbl, M) in feeds, dm in dimset, ts in (1, 2)
    try
        r = diamond_dual_rect(Matrix{ComplexF64}(M), dm[1], dm[2], ts)
        rat = r[1] / truth
        flag = abs(r[1] - truth) <= 1e-5 ? "  <== TRUTH" : ""
        abs(r[1] - truth) <= 1e-5 && push!(dual_hits, (flbl, dm, ts, r[1]))
        @printf("    %-14s %-10s %-8d %-14.8f %.6f%s\n", flbl, string(dm), ts, r[1], rat, flag)
    catch e
        @printf("    %-14s %-10s %-8d FAILED %s\n", flbl, string(dm), ts, sprint(showerror, e))
    end
end

println("\n  PRIMAL hits (combos == truth):")
for h in prim_hits; @printf("    feed=%s dims=%s rho_on=%s\n", h[1], string(h[2]), String(h[3])); end
println("  DUAL hits (combos == truth):")
for h in dual_hits; @printf("    feed=%s dims=%s tr_sys=%d\n", h[1], string(h[2]), h[3]); end

# ---- confirm the found convention also gives the *-iso truth 1.0 ----
# Use the matching primal+dual convention on a NON-square *-iso (emb 2->3).
println("\n  CONFIRM on *-iso embedding (2->3), expect ||v*||_⋄ = 1:")
W0 = ComplexF64[1 0; 0 1; 0 0]
vemb(X) = W0 * X * W0'                          # M_2 -> M_3
Jv = choi_convA(vemb, 2, 3)                      # source n_B=2 (major), target N=3 (minor)
Jvs = permutedims(Jv)                            # J(v*): source N=3 major, target n_B=2 minor
# For v*: input=N=3, output=n_B=2. Mirror Psi (input=3,output=2): feed J(v*) the
# SAME way the winning Psi combo fed J(Psi). Report all four combos for clarity.
for (flbl, M) in [("J(v*)", Jvs), ("transpose=J(v)", permutedims(Jvs))], dm in dimset
    try
        r = diamond_primal_rect(Matrix{ComplexF64}(M), dm[1], dm[2], :major)
        @printf("    primal feed=%-14s dims=%-8s -> %.8f%s\n", flbl, string(dm), r[1],
                abs(r[1]-1.0)<1e-5 ? "  <== 1.0" : "")
    catch e
        @printf("    primal feed=%-14s dims=%-8s FAILED %s\n", flbl, string(dm), sprint(showerror,e))
    end
end
println("="^78)
println("READ: the (feed, dims, rho_on/tr_sys) that hits TRUTH on the ASYMMETRIC CP")
println("map AND gives 1.0 on the *-iso is the correct O2 convention (no coincidence).")
println("Truth = sigma_max(A)^2 = ", round(truth, digits=8))
