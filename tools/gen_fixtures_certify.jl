# gen_fixtures_certify.jl — golden-master feasible-point generator for the TIGHT
# certified-arb cb-ball on eta = ||Phi^2-Phi||_cb (bead aic-m24, step 3a;
# design doc docs/research/cbnorm_tight_certifier.md). Emits tests/fixtures_certify.inc.h
# (committed, generated): for a small representative+adversarial subset of test
# channels, the SDP FEASIBLE POINTS that the arb certifier (step 3b) certifies by
# weak duality to bracket eta — so the C certifier test stays Julia-free.
#
# Two solves per channel on J = Choi(Phi^2-Phi) (Convention A):
#   MAX primal (X,P,Q): rigorous LOWER bound  eta >= (2/n) Re tr(J^dag X).
#   MIN dual  (Y0,Y1) : rigorous UPPER bound  eta <= (1/2)(||Tr_2(Y0)|| + ||Tr_2(Y1)||).
# DUAL normalization DETERMINED EMPIRICALLY (NOT the design-doc table): partial
# trace over sys 2 (the MINOR/input factor a in J = C[(i-1)n+a,(j-1)n+b]) and the
# constant 1/2. Pinned by the ASYMMETRIC paper example (sys 1 -> ratio 4.0 WRONG;
# sys 2 -> ratio 2.0 == eta_ref). See src/sdp.jl for the full pinning analysis.
#
# Before writing ANY fixture this script ASSERTS, for every channel, that BOTH
# emitted feasible points reproduce eta_ref to ~1e-6 — a bad solve cannot poison
# the substrate. Mirrors tools/gen_fixtures_d24.jl emit style. NO PYTHON.
#   julia --project=julia/env tools/gen_fixtures_certify.jl
using Convex, MosekTools, LinearAlgebra, Printf

# ---- tightened-MOSEK optimizer (verbatim from src/sdp.jl) ----
const MOSEK_TOL = (
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS"   => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS"   => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_MU_RED"  => 1.0e-14,
)
function mosek_tight()
    o = MosekTools.Mosek.Optimizer()
    for (k, v) in MOSEK_TOL
        MosekTools.MOI.set(o, MosekTools.MOI.RawOptimizerAttribute(k), v)
    end
    return o
end

# ---- Watrous MAX primal, EXPOSING (X,P,Q); eta = (2/n) optval ----
function watrous_primal(J::Matrix{ComplexF64}, n::Int)
    n2 = n * n
    X = ComplexVariable(n2, n2)
    P = HermitianSemidefinite(n2)
    Q = HermitianSemidefinite(n2)
    prob = maximize(real(tr(J' * X)),
                    [isposdef([P X; X' Q]), P + Q == Matrix{ComplexF64}(I, n2, n2)])
    solve!(prob, mosek_tight; silent = true)
    return (2.0 / n) * prob.optval, Matrix{ComplexF64}(evaluate(X)),
           Matrix{ComplexF64}(evaluate(P)), Matrix{ComplexF64}(evaluate(Q)),
           string(prob.status)
end

# ---- QETLAB/Watrous MIN dual, EXPOSING (Y0,Y1); eta = optval/2 ----
function watrous_dual(J::Matrix{ComplexF64}, n::Int)
    n2 = n * n
    Y0 = HermitianSemidefinite(n2)
    Y1 = HermitianSemidefinite(n2)
    obj = opnorm(partialtrace(Y0, 2, [n, n])) + opnorm(partialtrace(Y1, 2, [n, n]))
    prob = minimize(obj, [isposdef([Y0 -J; -J' Y1])])
    solve!(prob, mosek_tight; silent = true)
    return prob.optval / 2.0, Matrix{ComplexF64}(evaluate(Y0)),
           Matrix{ComplexF64}(evaluate(Y1)), string(prob.status)
end

# ---- Choi (Convention A), Heisenberg self-composition (verbatim, gen_fixtures_d24) ----
function choi_A(Kraus::Vector{Matrix{ComplexF64}}, n::Int)
    C = zeros(ComplexF64, n * n, n * n)
    for K in Kraus, i in 1:n, a in 1:n, j in 1:n, b in 1:n
        C[(i - 1) * n + a, (j - 1) * n + b] += conj(K[i, a]) * K[j, b]
    end
    return C
end
compose_self(Kr::Vector{Matrix{ComplexF64}}) = [Kb * Ka for Ka in Kr for Kb in Kr]
lambda_choi(Kr::Vector{Matrix{ComplexF64}}, n::Int) =
    choi_A(compose_self(Kr), n) - choi_A(Kr, n)

# ---- channel constructors (verbatim from gen_fixtures_d24.jl) ----
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
identity_kraus(d) = [Matrix{ComplexF64}(I, d, d)]
function phit_kraus(d, t)  # Phi_t = (1-t) id + t Dep_d
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end
complex_qubit_kraus() = [ComplexF64[0.6 0; 0.8im 0], ComplexF64[0 0.8im; 0 0.6]]
# Paper example (.tex:367-378): Phi(X) = P0 Tr(g0 X) + P1 Tr(g1 X) on C^2 (verbatim
# from tools/gen_fixtures_d24.jl / runtests.jl). ASYMMETRIC channel: its MIN-dual
# optimum has marginals with ||Tr_R(Y)-Tr_L(Y)||_F ~ 0.134 (>> the symmetric corpus'
# ~1e-13), so the partial-trace DIRECTION (Risk #5) acquires teeth — swapping
# partial_trace_right->left doubles hi (Tr_L gives hi = 2*eta), breaking tightness.
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

# ---- partial trace over sys 2 (minor/input factor) for a plain N x N matrix ----
# Tr_2(M)[i,j] = sum_a M[(i-1)n+a, (j-1)n+a]. Used only to re-verify the dual bound.
function ptrace2(M::Matrix{ComplexF64}, n::Int)
    R = zeros(ComplexF64, n, n)
    for i in 1:n, j in 1:n, a in 1:n
        R[i, j] += M[(i - 1) * n + a, (j - 1) * n + a]
    end
    return R
end
# ---- partial trace over sys 1 (major/output factor) — the WRONG direction; used
# ONLY to MEASURE the marginal asymmetry ||Tr_R - Tr_L||_F that gives the direction
# its teeth. Tr_1(M)[a,b] = sum_i M[(i-1)n+a, (i-1)n+b].
function ptrace1(M::Matrix{ComplexF64}, n::Int)
    R = zeros(ComplexF64, n, n)
    for a in 1:n, b in 1:n, i in 1:n
        R[a, b] += M[(i - 1) * n + a, (i - 1) * n + b]
    end
    return R
end
maxeig_herm(M) = maximum(real.(eigvals(Hermitian((M + M') / 2))))

# ---- fixture corpus: representative + adversarial subset (design doc validation) ----
# phi_t_0p3    n=2 nonzero baseline
# complex_qubit n=2 conjugation guard (genuinely complex J^dag, Tr_2)
# phi_t4_0p2   n=4 dimension canary (largest n -> 16x16 blocks)
# id_2         eta=0 oracle (LOWER bound must be ~0; dispatch boundary)
# paper_0p1    n=2 ASYMMETRIC-marginal fixture (closes GAP 1): the MIN-dual optimum
#              has ||Tr_R(Y)-Tr_L(Y)||_F ~ 0.134, giving the partial-trace DIRECTION
#              its teeth (Risk #5); swapping right->left makes hi = 2*eta.
struct CFix
    name::String
    n::Int
    Kr::Vector{Matrix{ComplexF64}}
    eta_ref::Float64
    eta_is_zero::Bool
end
dep_id_norm(d) = 2 * (1 - 1 / d^2)
fixtures = CFix[
    CFix("phi_t_0p3",     2, phit_kraus(2, 0.3),    0.3 * 0.7 * dep_id_norm(2), false),
    CFix("complex_qubit", 2, complex_qubit_kraus(), 1.28,                       false),
    CFix("phi_t4_0p2",    4, phit_kraus(4, 0.2),    0.2 * 0.8 * dep_id_norm(4), false),
    CFix("paper_0p1",     2, paper_example_kraus(0.1), 0.1 * sqrt(0.9),         false),
    CFix("id_2",          2, identity_kraus(2),     0.0,                        true),
]

# ---- solve + ASSERT both feasible points reproduce eta_ref BEFORE emitting ----
struct Solved
    f::CFix
    J::Matrix{ComplexF64}
    X::Matrix{ComplexF64}; P::Matrix{ComplexF64}; Q::Matrix{ComplexF64}
    Y0::Matrix{ComplexF64}; Y1::Matrix{ComplexF64}
    lo::Float64; hi::Float64   # bounds RECOMPUTED from the emitted matrices
end
function solve_and_check(f::CFix)
    J = lambda_choi(f.Kr, f.n); n = f.n
    ep, X, P, Q, ps = watrous_primal(J, n)
    ed, Y0, Y1, ds = watrous_dual(J, n)
    ps in ("OPTIMAL", "SLOW_PROGRESS") || error("$(f.name): primal status $ps")
    ds in ("OPTIMAL", "SLOW_PROGRESS") || error("$(f.name): dual status $ds")
    # Recompute the weak-duality bounds straight from the EMITTED matrices (this is
    # exactly what the arb certifier 3b does), then assert they bracket eta_ref.
    lo = (2.0 / n) * real(tr(J' * X))
    hi = (maxeig_herm(ptrace2(Y0, n)) + maxeig_herm(ptrace2(Y1, n))) / 2.0
    tol = 1e-6
    abs(lo - f.eta_ref) <= tol ||
        error("FIXTURE POISON $(f.name): LOWER from emitted X = $lo, eta_ref = $(f.eta_ref)")
    abs(hi - f.eta_ref) <= tol ||
        error("FIXTURE POISON $(f.name): UPPER from emitted Y = $hi, eta_ref = $(f.eta_ref)")
    # the SOLVED optvals must also agree (catches a normalization drift)
    abs(ep - f.eta_ref) <= tol || error("FIXTURE POISON $(f.name): eta_primal=$ep ref=$(f.eta_ref)")
    abs(ed - f.eta_ref) <= tol || error("FIXTURE POISON $(f.name): eta_dual=$ed ref=$(f.eta_ref)")
    # MARGINAL ASYMMETRY of the MIN-dual optimum (the partial-trace DIRECTION teeth,
    # Risk #5 / GAP 1): ||Tr_R(Y)-Tr_L(Y)||_F. ~1e-13 on a symmetric channel (the
    # swap is then SILENT); paper_0p1 must be >> that or the direction test has no
    # teeth. We POISON-GUARD any fixture whose name advertises asymmetry.
    asym = norm(ptrace2(Y0, n) - ptrace1(Y0, n)) + norm(ptrace2(Y1, n) - ptrace1(Y1, n))
    if startswith(f.name, "paper")
        asym > 0.05 || error("ASYMMETRY POISON $(f.name): marginal asym=$asym <= 0.05 " *
                             "(the direction-swap teeth would be inert)")
    end
    @printf("  %-14s n=%d  eta_ref=%.6e  LOWER(X)=%.6e  UPPER(Y)=%.6e  gap=%.2e  asym=%.4e\n",
            f.name, n, f.eta_ref, lo, hi, abs(ep - ed), asym)
    return Solved(f, J, X, P, Q, Y0, Y1, lo, hi)
end

# ---- emit the C header ----
function flat_rowmajor(M::Matrix{ComplexF64})  # [p*N + q], N = size
    N = size(M, 1)
    re = Float64[]; im = Float64[]
    for p in 1:N, q in 1:N
        push!(re, real(M[p, q])); push!(im, imag(M[p, q]))
    end
    return re, im
end
function wr(io, name, v)
    print(io, "static const double $(name)[] = {")
    print(io, join((@sprintf("%.17e", x) for x in v), ", "))
    println(io, "};")
end

function emit(io, solved::Vector{Solved})
    println(io, "/* Auto-generated by tools/gen_fixtures_certify.jl - DO NOT EDIT.")
    println(io, " * Golden-master SDP FEASIBLE POINTS for the tight certified-arb")
    println(io, " * cb-ball on eta = ||Phi^2-Phi||_cb (bead aic-m24 step 3a; design doc")
    println(io, " * docs/research/cbnorm_tight_certifier.md). For each channel: Kraus ops, eta_ref,")
    println(io, " * and two SDP feasible points the arb certifier (3b) certifies by weak")
    println(io, " * duality to bracket eta:")
    println(io, " *   MAX primal (X,P,Q): LOWER eta >= (2/n) Re tr(J^dag X).")
    println(io, " *   MIN dual  (Y0,Y1) : UPPER eta <= (1/2)(||Tr_2(Y0)|| + ||Tr_2(Y1)||),")
    println(io, " *                       Tr_2 = partial trace over sys 2 (input factor a).")
    println(io, " * All matrices are N x N (N = n*n), flat ROW-MAJOR [p*N + q]. The")
    println(io, " * generator ASSERTS both points reproduce eta_ref (~1e-6) before writing,")
    println(io, " * so a bad solve cannot poison the substrate. NO PYTHON. */")
    println(io, "#ifndef AIC_FIXTURES_CERTIFY_INC_H")
    println(io, "#define AIC_FIXTURES_CERTIFY_INC_H\n")
    println(io, "typedef struct {")
    println(io, "    const char *name;")
    println(io, "    long n;             /* dim_K == dim_H */")
    println(io, "    long r;             /* number of Kraus ops */")
    println(io, "    const double *kraus_re;  /* flat [op*n*n + i*n + j] */")
    println(io, "    const double *kraus_im;")
    println(io, "    /* MAX-primal feasible point (LOWER bound), each N x N, [p*N + q], N=n*n */")
    println(io, "    const double *X_re;  const double *X_im;")
    println(io, "    const double *P_re;  const double *P_im;")
    println(io, "    const double *Q_re;  const double *Q_im;")
    println(io, "    /* MIN-dual feasible point (UPPER bound), each N x N, [p*N + q] */")
    println(io, "    const double *Y0_re; const double *Y0_im;")
    println(io, "    const double *Y1_re; const double *Y1_im;")
    println(io, "    double eta_ref;      /* ||Phi^2-Phi||_diamond (Watrous SDP) */")
    println(io, "    int eta_is_zero;     /* 1 for idempotent (eta=0) oracles */")
    println(io, "} aic_certify_fixture;\n")

    for (idx, s) in enumerate(solved)
        f = s.f; n = f.n; r = length(f.Kr)
        kr = Float64[]; ki = Float64[]
        for K in f.Kr, i in 1:n, j in 1:n
            push!(kr, real(K[i, j])); push!(ki, imag(K[i, j]))
        end
        wr(io, "aic_cf_kr_$idx", kr); wr(io, "aic_cf_ki_$idx", ki)
        for (tag, M) in (("X", s.X), ("P", s.P), ("Q", s.Q), ("Y0", s.Y0), ("Y1", s.Y1))
            re, im = flat_rowmajor(M)
            wr(io, "aic_cf_$(tag)r_$idx", re); wr(io, "aic_cf_$(tag)i_$idx", im)
        end
        println(io)
    end

    println(io, "static const aic_certify_fixture aic_certify_fixtures[] = {")
    for (idx, s) in enumerate(solved)
        f = s.f; n = f.n; r = length(f.Kr)
        row  = @sprintf("    {\"%s\", %d, %d, aic_cf_kr_%d, aic_cf_ki_%d, ", f.name, n, r, idx, idx)
        row *= @sprintf("aic_cf_Xr_%d, aic_cf_Xi_%d, aic_cf_Pr_%d, aic_cf_Pi_%d, ", idx, idx, idx, idx)
        row *= @sprintf("aic_cf_Qr_%d, aic_cf_Qi_%d, ", idx, idx)
        row *= @sprintf("aic_cf_Y0r_%d, aic_cf_Y0i_%d, aic_cf_Y1r_%d, aic_cf_Y1i_%d, ", idx, idx, idx, idx)
        row *= @sprintf("%.17e, %d},", f.eta_ref, f.eta_is_zero ? 1 : 0)
        println(io, row)
    end
    println(io, "};")
    println(io, "#define AIC_CERTIFY_NFIX ((int)(sizeof(aic_certify_fixtures)/sizeof(aic_certify_fixtures[0])))\n")
    println(io, "#endif /* AIC_FIXTURES_CERTIFY_INC_H */")
end

println("solving + asserting feasible points reproduce eta_ref (fail-loud poison guard):")
solved = [solve_and_check(f) for f in fixtures]
out = joinpath(@__DIR__, "..", "tests", "fixtures_certify.inc.h")
open(out, "w") do io
    emit(io, solved)
end
println("wrote ", normpath(out), "  (", length(solved), " channels)")
