# smoke_gadc_cbnorm.jl — quick smoke test: eta-idempotence defect
# ||Phi^2 - Phi||_cb of the generalised amplitude damping channel (GADC).
#
# Context (see HANDOFF.md / ALGORITHM.md module cb-norm, bead aic-d24):
#   The plain cb-norm of a channel is trivially 1; the quantity this project
#   computes is the eta-idempotence defect ||Phi^2-Phi||_cb = ||Phi^2-Phi||_diamond
#   (approximate_algebras.tex:347-354). This script reuses the SAME Watrous 2012
#   SDP (Convex.jl + MOSEK, complex-native) and the SAME Convention-A Choi /
#   Heisenberg-composition conventions as tools/gen_fixtures_d24.jl, so it is a
#   genuine exercise of the cb-norm machinery on a real, non-fixture channel.
#
#   NORMALIZATION (load-bearing, same as the generator): a Convention-A Choi has
#   trace n for a UCP map, but the Watrous SDP is trace-1 calibrated, so
#   ||Lambda||_diamond = (2/n) * SDP_optval. (n=2 here.)
#
# NO PYTHON (user directive). Run serially (single Julia process):
#   julia --project=julia/env tools/smoke_gadc_cbnorm.jl
#
# GADC (qubit, parameters gamma in [0,1] damping, p in [0,1] ground population).
# STATE-picture Kraus T(rho)=sum A_i rho A_i^dag:
#   A0 = sqrt(p)   [[1,0],[0,sqrt(1-g)]]    A1 = sqrt(p)   [[0,sqrt(g)],[0,0]]
#   A2 = sqrt(1-p) [[sqrt(1-g),0],[0,1]]    A3 = sqrt(1-p) [[0,0],[sqrt(g),0]]
# T is trace-preserving (sum A_i^dag A_i = 1). The OBSERVABLE map this project
# uses is the Heisenberg dual Phi = T^*, Phi(X) = sum_i A_i^dag X A_i, i.e. the
# aic_ucp Kraus K_a are exactly the A_i (Phi(X)=sum_a K_a^dag X K_a). Phi is UCP
# and UNITAL. cb-norm is adjoint-invariant, so ||Phi^2-Phi||_cb = ||T^2-T||_cb.
#
# Two exact eta=0 anchors (the eta=0 oracle, CLAUDE.md cross-check ladder rung 3):
#   gamma=0  -> A0=sqrt(p) 1, A2=sqrt(1-p) 1  -> Phi = identity      -> eta = 0
#   gamma=1  -> T(rho) = Tr(rho) * diag(p,1-p) (reset to fixed state) -> T^2 = T
#                                                                      -> eta = 0
# Interior 0<gamma<1 is NOT idempotent -> eta > 0 (rises from 0, returns to 0).

using Convex, MosekTools, LinearAlgebra, Printf

# ---- Watrous diamond-norm SDP (verbatim from gen_fixtures_d24.jl) ----
function diamond_norm_watrous(J::Matrix{ComplexF64}, n::Int)
    n2 = n * n
    X = ComplexVariable(n2, n2)
    P = HermitianSemidefinite(n2)
    Q = HermitianSemidefinite(n2)
    blk = [P X; X' Q]
    prob = maximize(real(tr(J' * X)),
                    [isposdef(blk),
                     P + Q == Matrix{ComplexF64}(I, n2, n2)])
    solve!(prob, MosekTools.Mosek.Optimizer; silent = true)
    return (2.0 / n) * prob.optval, string(prob.status)
end

# ---- Convention-A Choi and Heisenberg self-composition (verbatim) ----
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
function eta_of(Kr::Vector{Matrix{ComplexF64}}, n::Int)
    v, st = diamond_norm_watrous(lambda_choi(Kr, n), n)
    st == "OPTIMAL" || error("SDP not OPTIMAL (status=$st)")
    return v
end

# ---- GADC Kraus (aic_ucp Heisenberg-dual layout K_a = A_i) ----
function gadc_kraus(gamma::Float64, p::Float64)
    sp, sq = sqrt(p), sqrt(1 - p)
    sg, sg1 = sqrt(gamma), sqrt(1 - gamma)
    A0 = ComplexF64[sp 0; 0 sp*sg1]
    A1 = ComplexF64[0 sp*sg; 0 0]
    A2 = ComplexF64[sq*sg1 0; 0 sq]
    A3 = ComplexF64[0 0; sq*sg 0]
    return Matrix{ComplexF64}[A0, A1, A2, A3]
end

# unital sanity: sum_a K_a^dag K_a == 1_2 (Phi UCP unital, .tex:1571)
function unital_defect(Kr::Vector{Matrix{ComplexF64}})
    S = sum(K' * K for K in Kr)
    return opnorm(S - Matrix{ComplexF64}(I, 2, 2))
end

# ---- run: gamma sweep at two p values, plus the two eta=0 anchors ----
println("GADC  eta = ||Phi^2 - Phi||_cb  (Watrous SDP + MOSEK, n=2)\n")
println(rpad("gamma", 8), rpad("p", 6), rpad("eta", 16), "unital_defect")

anchor_tol  = 1e-7    # gamma in {0,1} must hit machine/SDP zero
interior_tol = 1e-3   # interior must be clearly nonzero

for p in (0.5, 0.25)
    etas = Float64[]
    for gamma in 0.0:0.1:1.0
        Kr = gadc_kraus(gamma, p)
        ud = unital_defect(Kr)
        ud < 1e-12 || error("GADC not unital at gamma=$gamma p=$p (defect=$ud)")
        eta = eta_of(Kr, 2)
        push!(etas, eta)
        @printf("%-8.1f%-6.2f%-16.9e%.2e\n", gamma, p, eta, ud)
    end
    # eta=0 anchors (gamma=0 identity, gamma=1 reset-to-fixed-state idempotent)
    abs(etas[1])   <= anchor_tol || error("ANCHOR FAIL p=$p gamma=0: eta=$(etas[1]) (expected ~0, identity)")
    abs(etas[end]) <= anchor_tol || error("ANCHOR FAIL p=$p gamma=1: eta=$(etas[end]) (expected ~0, idempotent reset)")
    # interior must be nonzero (GADC is not idempotent for 0<gamma<1)
    maximum(etas) >= interior_tol || error("ANCHOR FAIL p=$p: interior eta never exceeds $interior_tol")
    gpk = (0.0:0.1:1.0)[argmax(etas)]
    @printf("  -> p=%.2f: eta(0)=%.2e (id), eta(1)=%.2e (idemp), peak eta=%.6e at gamma=%.1f\n\n",
            p, etas[1], etas[end], maximum(etas), gpk)
end

println("smoke_gadc_cbnorm: OK (eta=0 anchors at gamma=0,1; interior eta>0 confirmed)")
