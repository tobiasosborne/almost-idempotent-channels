# probe_o2_pin2.jl — DEFINITIVE O2 convention pin (bead aic-pjr), using the
# GOLDEN RULE established by the CP oracle (probe_o2_diag2.jl):
#
#   To compute ||f||_⋄ for f: M_in -> M_out, build the Convention-A Choi
#   J = choi_convA(f, in, out)  (INPUT-major, OUTPUT-minor), feed to
#   diamond_*_rect(J, in, out, ...). Raw optval = ||f||_⋄ EXACTLY (factor = 1).
#
# Verified by the asymmetric CP map Psi(Y)=A'YA (in=3,out=2): raw = sigma_max(A)^2
# (independent truth) with dims (in,out)=(3,2). This probe (a) re-confirms the
# rule on the *-iso for ||v||_cb and ||v^-1||_cb (build the ADJOINT's Choi
# directly, no transpose), and (b) PINS the dual partial-trace direction tr_sys
# on a genuinely asymmetric NON-CP HP map with the CORRECT dims, via strong
# duality (the tr_sys whose dual MIN matches the primal MAX is correct).
#   julia --project=julia/env tools/probe_o2_pin2.jl     (serial; NO PYTHON)

using Convex, MosekTools, LinearAlgebra, Printf, Random
include(joinpath(@__DIR__, "..", "julia", "AlmostIdempotentChannels.jl",
                 "src", "sdp.jl"))

choi_convA(fmap, a::Int, b::Int) = begin       # J[s*b+i,t*b+j]=f(E_st)[i,j]
    J = zeros(ComplexF64, a * b, a * b)
    for s in 1:a, t in 1:a
        Est = zeros(ComplexF64, a, a); Est[s, t] = 1.0
        F = fmap(Est)
        for i in 1:b, j in 1:b
            J[(s - 1) * b + i, (t - 1) * b + j] = F[i, j]
        end
    end
    J
end

Random.seed!(20260531)
println("="^78); println("O2 DEFINITIVE PIN (golden rule: choi_convA(f,in,out), dims (in,out), factor 1)"); println("="^78)

# ---- (a) *-iso, n_B=2 -> N=3 ----
N, n_B = 3, 2
W = Matrix(qr(randn(ComplexF64, N, N)).Q)[:, 1:n_B]      # 3x2 isometry
@assert norm(W' * W - I) < 1e-12
# ||v||_cb = ||v*||_⋄,  v*: M_N -> M_{n_B},  v*(Y) = W' Y W  (in=N, out=n_B).
vstar(Y) = W' * Y * W
Jvstar = choi_convA(vstar, N, n_B)                        # in=N=3 major, out=n_B=2 minor
# ||v^-1||_cb = ||(v^-1)*||_⋄, (v^-1)*: M_{n_B} -> M_N, (v^-1)*(X)=W X W' (in=n_B,out=N).
vinvstar(X) = W * X * W'
Jvinvstar = choi_convA(vinvstar, n_B, N)                  # in=n_B=2 major, out=N=3 minor
println("\n(a) *-iso: ||v||_cb and ||v^-1||_cb must both be 1.0 (factor 1, dims=(in,out)):")
for (lbl, J, din, dout) in [("||v||_cb   (v*)", Jvstar, N, n_B),
                            ("||v^-1||_cb ((v^-1)*)", Jvinvstar, n_B, N)]
    rp = diamond_primal_rect(Matrix{ComplexF64}(J), din, dout, :major)
    rd1 = diamond_dual_rect(Matrix{ComplexF64}(J), din, dout, 1)
    rd2 = diamond_dual_rect(Matrix{ComplexF64}(J), din, dout, 2)
    @printf("  %-22s dims=(%d,%d): primal=%.8f  dual1=%.8f  dual2=%.8f\n",
            lbl, din, dout, rp[1], rd1[1], rd2[1])
end

# ---- (b) CP oracle re-confirm (independent truth) ----
A = randn(ComplexF64, 3, 2); psi(Y) = A' * Y * A
truthCP = maximum(real.(eigvals(Hermitian(A * A'))))
Jcp = choi_convA(psi, 3, 2)
rp = diamond_primal_rect(Jcp, 3, 2, :major)
rd1 = diamond_dual_rect(Jcp, 3, 2, 1); rd2 = diamond_dual_rect(Jcp, 3, 2, 2)
println("\n(b) CP oracle Psi (in=3,out=2): truth sigma_max(A)^2 = $(round(truthCP,digits=8))")
@printf("  primal=%.8f  dual1=%.8f  dual2=%.8f  (all == truth -> factor 1 confirmed)\n",
        rp[1], rd1[1], rd2[1])

# ---- (c) DIRECTION pin: asymmetric NON-CP HP map, correct dims (3,2) ----
# Jg = Jcp + Hermitian perturbation (keeps it a valid HP-map Choi, breaks CP & symmetry).
Hh = randn(ComplexF64, 6, 6); Jg = Jcp + 0.5 * (Hh + Hh') / 2
@assert norm(Jg - Jg') < 1e-12
println("\n(c) asymmetric NON-CP HP map (dims (3,2)): correct tr_sys = the one matching primal")
pm = diamond_primal_rect(Matrix{ComplexF64}(Jg), 3, 2, :major)
pmi = diamond_primal_rect(Matrix{ComplexF64}(Jg), 3, 2, :minor)
d1 = diamond_dual_rect(Matrix{ComplexF64}(Jg), 3, 2, 1)
d2 = diamond_dual_rect(Matrix{ComplexF64}(Jg), 3, 2, 2)
truth = pm[1]
@printf("  primal rho_on=:major = %.8f  (TRUTH, no partial trace)\n", pm[1])
@printf("  primal rho_on=:minor = %.8f  (should match :major if placement irrelevant)\n", pmi[1])
@printf("  dual tr_sys=1 = %.8f   gap to truth = %.3e\n", d1[1], abs(d1[1]-truth))
@printf("  dual tr_sys=2 = %.8f   gap to truth = %.3e\n", d2[1], abs(d2[1]-truth))
cdir = abs(d1[1]-truth) <= abs(d2[1]-truth) ? 1 : 2
teeth = abs(d1[1]-d2[1])
@printf("  => CORRECT dual tr_sys = %d ;  teeth |dual1-dual2| = %.6f\n", cdir, teeth)
if teeth < 1e-3
    println("  *** WARNING: teeth tiny — direction inert on this fixture; strengthen perturbation.")
end

# tr_sys -> traced factor.  dims=[in=d_maj, out=d_min]: sys 1 = input(major), sys 2 = output(minor).
println("\n" * "="^78)
println("PINNED O2 CONVENTION:")
println("  recipe: build J = choi_convA(adjoint-map, in, out) [INPUT-major]; feed diamond_*_rect(J, in, out)")
println("  normalization factor : 1  (raw optval = ||.||_⋄ directly)")
@printf("  dual tr_sys          : %d  (traces the %s factor)\n", cdir, cdir==1 ? "INPUT/major" : "OUTPUT/minor")
println("  primal density on    : :major (= INPUT factor)  [verify :major==:minor above]")
println("  ||v||_cb  = ||v*||_⋄        : feed choi_convA(v*,   N,   n_B), dims (N,   n_B)")
println("  ||v^-1||_cb = ||(v^-1)*||_⋄ : feed choi_convA((v^-1)*, n_B, N), dims (n_B, N)")
println("="^78)
