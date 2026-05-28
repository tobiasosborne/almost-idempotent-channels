# gen_fixtures_d24.jl — golden-master generator for the eta-idempotence defect
# ||Phi^2-Phi||_cb = ||Phi^2-Phi||_diamond (bead aic-d24, approximate_algebras.tex
# :347-354). Emits tests/fixtures_d24.inc.h (committed, generated): for each test
# channel, the Kraus operators (matching the C aic_ucp layout), the Convention-A
# Choi of Lambda=Phi^2-Phi (so the C test can cross-check its OWN Choi substrate),
# and the reference diamond-norm eta from the Watrous SDP solved with MOSEK.
#
# NO PYTHON (user directive). Run serially (single Julia process):
#   julia --project=julia/env tools/gen_fixtures_d24.jl
#
# The SDP is the Watrous 2012 program (arXiv:1207.5726), Convex.jl + MOSEK,
# complex-native. NORMALIZATION (load-bearing): Convention-A Choi has trace=n for
# a UCP map; the Watrous SDP is trace-1 calibrated, so ||Lambda||_diamond =
# (2/n)*SDP_optval. Before writing ANY fixture this script ASSERTS the analytic
# anchors (id->0; Dep_d-id_d -> 2(1-1/d^2); Phi_t -> t(1-t)*1.5; paper example
# .tex:378 -> eta*sqrt(1-eta)) — it fails loud (error()) if the SDP disagrees,
# so a broken solver/normalization can never silently poison the fixtures.
using Convex, MosekTools, LinearAlgebra
using Printf

# ---- the Watrous diamond-norm SDP (VERIFIED in assert_anchors() below) ----
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

# ---- Choi (Convention A) and Heisenberg self-composition ----
# C[i*n+a, j*n+b] = sum_x conj(Kx[i,a]) Kx[j,b]  (1-based here).
function choi_A(Kraus::Vector{Matrix{ComplexF64}}, n::Int)
    C = zeros(ComplexF64, n * n, n * n)
    for K in Kraus, i in 1:n, a in 1:n, j in 1:n, b in 1:n
        C[(i - 1) * n + a, (j - 1) * n + b] += conj(K[i, a]) * K[j, b]
    end
    return C
end
# Kraus of Phi o Phi (Heisenberg): {K_b K_a}.
compose_self(Kr::Vector{Matrix{ComplexF64}}) =
    [Kb * Ka for Ka in Kr for Kb in Kr]

# Lambda=Phi^2-Phi Choi and its diamond norm.
function lambda_choi(Kr::Vector{Matrix{ComplexF64}}, n::Int)
    return choi_A(compose_self(Kr), n) - choi_A(Kr, n)
end
function eta_of(Kr::Vector{Matrix{ComplexF64}}, n::Int)
    v, st = diamond_norm_watrous(lambda_choi(Kr, n), n)
    st == "OPTIMAL" || error("SDP not OPTIMAL (status=$st)")
    return v
end

# ---- channel constructors (return Vector{Matrix{ComplexF64}} Kraus ops) ----
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
identity_kraus(d) = [Matrix{ComplexF64}(I, d, d)]
dephasing_kraus(d) = [e(d, a) * e(d, a)' for a in 1:d]
function depolarizing_kraus(d)
    [ (e(d, i) * e(d, j)') / sqrt(d) for i in 1:d for j in 1:d ]
end
function block_cond_exp_kraus()  # C^4 = C^2 (+) C^2, project to block-diagonal
    P0 = ComplexF64[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]
    P1 = ComplexF64[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]
    [P0, P1]
end
function phit_kraus(d, t)  # Phi_t = (1-t) id + t Dep_d
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end
# Paper example (.tex:367-378): Phi(X)=P0 Tr(g0 X)+P1 Tr(g1 X) on C^2.
# Kraus for X -> P Tr(g X) with P=u u^dag rank-1: K_k = sqrt(lam_k) v_k u^dag,
# g = sum_k lam_k v_k v_k^dag. Then sum_k K_k^dag X K_k = Tr(g X) P (probe-verified).
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
# A genuinely rank-deficient (low-rank) Lambda fixture: a UNITAL channel on C^3
# that acts nontrivially only on the upper 2-block and as the identity on e2.
# Phi = block-diag(Phi_t on span(e0,e1),  id on span(e2)) with Phi_t=(1-t)id+t Dep_2
# embedded in the 2x2 block; an extra |2><2| Kraus carries the identity on e2.
# Because Phi acts trivially on e2, Lambda=Phi^2-Phi vanishes on that subspace, so
# Choi(Lambda) is genuinely RANK-DEFICIENT (real zero eigenvalues): at t=0.2 the
# 9x9 Lambda-Choi has rank 4 (eigvals {-0.24, 0(x5), 0.08(x3)}). Unital
# (sum K_a^dag K_a = 1_3 to 1e-16). This replaces the earlier mislabeled
# phit_kraus(3,0.05) (which had a FULL-rank Lambda-Choi).
function block_lowrank3_kraus(t)
    embed(M) = (Z = zeros(ComplexF64, 3, 3); Z[1:2, 1:2] = M; Z)
    Ks = Matrix{ComplexF64}[embed(sqrt(1 - t) * Matrix{ComplexF64}(I, 2, 2))]
    for i in 1:2, j in 1:2
        push!(Ks, embed(sqrt(t) * (e(2, i) * e(2, j)') / sqrt(2)))
    end
    push!(Ks, e(3, 3) * e(3, 3)')  # identity on the e2 block (|2><2|)
    return Ks
end
# Genuinely COMPLEX, asymmetric, UNITAL qubit channel (matches the C test's
# make_complex_channel): K_0=[[0.6,0],[0.8i,0]], K_1=[[0,0.8i],[0,0.6]].
# sum_a K_a^dag K_a = 1_2 (unital). Phi is NON-idempotent and — crucially — its
# Lambda=Phi^2-Phi Choi has genuinely NONZERO IMAGINARY off-diagonal entries
# (max|Im| ~ 0.61), so the C-vs-Julia Choi cross-check now exercises the
# conjugation in aic_ucp_kraus_to_choi: dropping the conj changes the Lambda-Choi
# by ~1.28 in operator norm, far above the 1e-11 fixture tolerance. max|Im(Kraus)|
# = 0.8 > 0.1. eta from MOSEK (~1.28); validated by the SDP + the Choi cross-check.
function complex_qubit_kraus()
    K0 = ComplexF64[0.6 0; 0.8im 0]
    K1 = ComplexF64[0 0.8im; 0 0.6]
    return [K0, K1]
end

# ---- fixture corpus ----
struct Fixture
    name::String
    n::Int
    Kr::Vector{Matrix{ComplexF64}}
    eta_is_zero::Bool
end
fixtures = Fixture[
    Fixture("id_1",            1, identity_kraus(1),     true),
    Fixture("id_2",            2, identity_kraus(2),     true),
    Fixture("id_3",            3, identity_kraus(3),     true),
    Fixture("dep_2",           2, depolarizing_kraus(2), true),
    Fixture("dep_3",           3, depolarizing_kraus(3), true),
    Fixture("dephasing_M2",    2, dephasing_kraus(2),    true),
    Fixture("block_cond_exp4", 4, block_cond_exp_kraus(),true),
    Fixture("phi_t_0p3",       2, phit_kraus(2, 0.3),    false),
    Fixture("phi_t_0p1",       2, phit_kraus(2, 0.1),    false),
    Fixture("phi_t_0p01",      2, phit_kraus(2, 0.01),   false),
    Fixture("phi_t_0p001",     2, phit_kraus(2, 0.001),  false),
    Fixture("paper_ex_0p1",    2, paper_example_kraus(0.1),   false),
    Fixture("paper_ex_0p01",   2, paper_example_kraus(0.01),  false),
    Fixture("phi_t3_0p05",     3, phit_kraus(3, 0.05),        false),  # n=3 canary
    Fixture("phi_t4_0p2",      4, phit_kraus(4, 0.2),    false),       # n=4 canary
    Fixture("block_lowrank3",  3, block_lowrank3_kraus(0.2),  false),  # rank-deficient Lambda-Choi
    Fixture("complex_qubit",   2, complex_qubit_kraus(),      false),  # complex Kraus: guards Choi conjugation
]

# ---- ASSERT the analytic anchors BEFORE emitting anything (fail loud) ----
# Poison-pill anchors covering EVERY dimension-varying fixture (d=2,3,4). The
# general formulae the SDP must reproduce (else the fixtures are silently
# poisoned): ||Dep_d - id_d||_diamond = 2(1 - 1/d^2), and since
# Phi_t = (1-t)id + t Dep_d gives Phi_t^2 - Phi_t = t(1-t)(Dep_d - id_d), also
# ||Phi_t^2 - Phi_t||_diamond = t(1-t) * 2(1 - 1/d^2). Both d=3 (16/9) and d=4
# (15/8 = 1.875) are now guarded, not just d=2.
function assert_anchors()
    tol = 1e-6
    chk(name, got, want) = abs(got - want) <= tol ||
        error("ANCHOR FAIL $name: got $got want $want")
    dep_id_norm(d) = 2 * (1 - 1 / d^2)            # ||Dep_d - id_d||_diamond
    eta_dep_diff(d) = begin                       # SDP of Dep_d - id_d
        C = choi_A(depolarizing_kraus(d), d) - choi_A(identity_kraus(d), d)
        v, _ = diamond_norm_watrous(C, d); v
    end
    chk("id_2->0", eta_of(identity_kraus(2), 2), 0.0)
    # Dep_d - id_d anchor: 1.5 (d=2), 16/9 (d=3), 15/8 (d=4).
    for d in (2, 3, 4)
        chk("Dep$d-id$d=2(1-1/$d^2)", eta_dep_diff(d), dep_id_norm(d))
    end
    # Phi_t poison-pill across all three dimensions that appear in the corpus.
    for (d, ts) in ((2, (0.3, 0.1, 0.01, 0.001)), (3, (0.05,)), (4, (0.2,)))
        for t in ts
            chk("phi_t d=$d t=$t", eta_of(phit_kraus(d, t), d),
                t * (1 - t) * dep_id_norm(d))
        end
    end
    for ep in (0.1, 0.01)
        chk("paper t=$ep", eta_of(paper_example_kraus(ep), 2), ep * sqrt(1 - ep))
    end
    # The complex fixture must have genuinely complex Kraus ops, else the C-vs-
    # Julia Choi cross-check cannot guard the conjugation in aic_ucp_kraus_to_choi.
    maxim = maximum(maximum(abs.(imag.(K))) for K in complex_qubit_kraus())
    maxim > 0.1 || error("ANCHOR FAIL complex_qubit: max|Im(Kraus)|=$maxim <= 0.1")
    println("complex_qubit max|Im(Kraus)| = ", maxim, " (> 0.1, guards Choi conj)")
    println("ANCHORS OK")
end

# ---- emit the C header ----
function emit(io, fixtures)
    println(io, "/* Auto-generated by tools/gen_fixtures_d24.jl - DO NOT EDIT.")
    println(io, " * Golden master for bead aic-d24 (||Phi^2-Phi||_cb, tex:347-354):")
    println(io, " * Kraus ops (aic_ucp layout), Convention-A Choi of Lambda=Phi^2-Phi,")
    println(io, " * and the diamond-norm reference eta from the Watrous SDP + MOSEK.")
    println(io, " * The C test cross-checks its OWN Choi(Lambda) against choi_re/im;")
    println(io, " * the eta value is the Julia golden master (validated vs analytic")
    println(io, " * anchors at generation time). NO PYTHON. */")
    println(io, "#ifndef AIC_FIXTURES_D24_INC_H")
    println(io, "#define AIC_FIXTURES_D24_INC_H\n")
    println(io, "typedef struct {")
    println(io, "    const char *name;")
    println(io, "    long n;            /* dim_K == dim_H */")
    println(io, "    long r;            /* number of Kraus ops */")
    println(io, "    const double *kraus_re;  /* flat [op*n*n + i*n + j] */")
    println(io, "    const double *kraus_im;")
    println(io, "    const double *choi_re;   /* Choi(Phi^2-Phi), flat [p*N + q], N=n*n */")
    println(io, "    const double *choi_im;")
    println(io, "    double eta_ref;          /* ||Phi^2-Phi||_diamond (Watrous SDP) */")
    println(io, "    int eta_is_zero;         /* 1 for idempotent (eta=0) oracles */")
    println(io, "} aic_d24_fixture;\n")

    # solve the SDP ONCE per fixture (it is the slow step); reuse the value.
    etas = [eta_of(f.Kr, f.n) for f in fixtures]

    for (idx, f) in enumerate(fixtures)
        n = f.n; N = n * n; r = length(f.Kr)
        eta = etas[idx]
        L = lambda_choi(f.Kr, n)
        # Kraus flat arrays
        kr = Float64[]; ki = Float64[]
        for K in f.Kr, i in 1:n, j in 1:n
            push!(kr, real(K[i, j])); push!(ki, imag(K[i, j]))
        end
        cr = Float64[]; ci = Float64[]
        for p in 1:N, q in 1:N
            push!(cr, real(L[p, q])); push!(ci, imag(L[p, q]))
        end
        wr(io, name, v) = begin
            print(io, "static const double $(name)[] = {")
            print(io, join((@sprintf("%.17e", x) for x in v), ", "))
            println(io, "};")
        end
        wr(io, "aic_d24_kr_$idx", kr)
        wr(io, "aic_d24_ki_$idx", ki)
        wr(io, "aic_d24_cr_$idx", cr)
        wr(io, "aic_d24_ci_$idx", ci)
        @printf("  %-16s n=%d r=%d  eta=%.6e  (%s)\n",
                f.name, n, r, eta, f.eta_is_zero ? "oracle" : "nonzero")
        println(io)
    end

    println(io, "static const aic_d24_fixture aic_d24_fixtures[] = {")
    for (idx, f) in enumerate(fixtures)
        n = f.n; r = length(f.Kr)
        eta = etas[idx]
        row = @sprintf("    {\"%s\", %d, %d, aic_d24_kr_%d, aic_d24_ki_%d, ",
                       f.name, n, r, idx, idx)
        row *= @sprintf("aic_d24_cr_%d, aic_d24_ci_%d, %.17e, %d},",
                        idx, idx, eta, f.eta_is_zero ? 1 : 0)
        println(io, row)
    end
    println(io, "};")
    println(io, "#define AIC_D24_NFIX ((int)(sizeof(aic_d24_fixtures)/sizeof(aic_d24_fixtures[0])))\n")
    println(io, "#endif /* AIC_FIXTURES_D24_INC_H */")
end

assert_anchors()
out = joinpath(@__DIR__, "..", "tests", "fixtures_d24.inc.h")
open(out, "w") do io
    emit(io, fixtures)
end
println("wrote ", normpath(out), "  (", length(fixtures), " fixtures)")
