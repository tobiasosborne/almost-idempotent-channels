# gen_fixtures_opspace_o2.jl — committed-fixture generator for the opspace O2
# certified cb-norm bound on the th_main_ext isomorphism v: B -> A (bead aic-pjr;
# approximate_algebras.tex:1538-1540; design docs/research/opspace_o2_design.md).
# Emits tests/fixtures_opspace_o2.inc.h (committed, generated): for a small corpus
# of UCP fixture channels, the Watrous diamond-norm SDP FEASIBLE POINTS (Y0,Y1) +
# reference values that the arb certifier (next task) restores to exact
# feasibility — so the eventual `make test` stays Julia-free.
#
# THE PINNED CONVENTION (design §0.5 GOLDEN RULE — SUPERSEDES the 2/N / 2/n_B
# hypotheses; pinned empirically by tools/probe_o2_pin2.jl against an asymmetric
# NON-CP closed-form truth + the eta=0 complete-isometry oracle):
#   ||v||_cb   = ||v*||_⋄        : feed J(v*)      to diamond_*_rect(d_maj=N,   d_min=n_B)
#   ||v^-1||_cb = ||(v^-1)*||_⋄  : feed J((v^-1)*) to diamond_*_rect(d_maj=n_B, d_min=N)
# with FACTOR 1 (raw optval = cb-norm EXACTLY), dual tr_sys = 2 (trace the MINOR
# factor, size d_min), primal rho_on = :major (density on the MAJOR/input factor).
# J(v*), J((v^-1)*) are built in C (aic_opspace_choi_v{star,invstar}) and handed
# over by the flat-double shim aic_opspace_choi_shim_d (build/libaic.so).
#
# POISON GUARDS before emitting ANY fixture (fail-loud, mirrors
# gen_fixtures_certify.jl): recompute the UPPER bound from the emitted (Y0,Y1) and
# assert it reproduces eta_ref to ~1e-6; assert primal ~= dual (strong duality);
# for eta=0 oracles assert eta_ref == 1.0; for the eta>0 fixture assert the
# WRONG-direction dual differs from eta_ref (the direction tooth has teeth).
#
# NO PYTHON. Determinism. Run serially (single Julia process), AFTER `make lib`:
#   make lib && julia --project=julia/env tools/gen_fixtures_opspace_o2.jl
using Convex, MosekTools, LinearAlgebra, Printf, Libdl

include(joinpath(@__DIR__, "..", "julia", "AlmostIdempotentChannels.jl",
                 "src", "sdp.jl"))   # diamond_dual_rect / diamond_primal_rect (PINNED)

const LIBAIC = abspath(joinpath(@__DIR__, "..", "build", "libaic.so"))
const CPREC = 256          # matches test_opspace.c CPREC

# ---- dlopen libaic (mirror the package __init__: system libgmp first, DEEPBIND) ----
function load_libaic()
    isfile(LIBAIC) || error("build/libaic.so not found at $LIBAIC — run `make lib` first")
    try
        Libdl.dlopen("/lib/x86_64-linux-gnu/libgmp.so.10",
                     Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL)
    catch
        try
            Libdl.dlopen("libgmp.so.10", Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL)
        catch
        end
    end
    h = Libdl.dlopen(LIBAIC, Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL | Libdl.RTLD_DEEPBIND)
    return Libdl.dlsym(h, :aic_opspace_choi_shim_d)
end

# ---- ccall the shim: rebuild v from flat Kraus, get J(v*), J((v^-1)*), N, n_B ----
# Returns (Jvs, Jvis, N, nB, dimB, iso) with both Choi as (N*n_B) x (N*n_B)
# ComplexF64 matrices, unflattened from the C row-major [p*sz+q] arrays.
function shim_choi(sym, kraus::Vector{Matrix{ComplexF64}}, n::Int, r::Int,
                   eps::Float64, prec::Int)
    # Flatten Kraus: kraus[a*n*n + i*n + j] = K_a[i+1,j+1] (Convention A, K-major).
    kr = Vector{Float64}(undef, r * n * n); ki = similar(kr)
    for a in 1:r
        K = kraus[a]; size(K) == (n, n) || error("Kraus $a size $(size(K)) != ($n,$n)")
        base = (a - 1) * n * n
        for i in 1:n, j in 1:n
            idx = base + (i - 1) * n + j
            kr[idx] = real(K[i, j]); ki[idx] = imag(K[i, j])
        end
    end
    n4 = n * n * n * n
    jvs_re = Vector{Float64}(undef, n4); jvs_im = Vector{Float64}(undef, n4)
    jvis_re = Vector{Float64}(undef, n4); jvis_im = Vector{Float64}(undef, n4)
    oN = Ref{Cint}(0); onB = Ref{Cint}(0); odimB = Ref{Cint}(0)
    oiso = Vector{Float64}(undef, 2)        # [0]=iso_def, [1]=eps actually used
    rc = ccall(sym, Cint,
               (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
                Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cdouble, Cint),
               jvs_re, jvs_im, jvis_re, jvis_im, oN, onB, odimB, oiso,
               kr, ki, n, r, eps, prec)
    rc == 0 || error("aic_opspace_choi_shim_d returned $rc")
    N = Int(oN[]); nB = Int(onB[]); dimB = Int(odimB[]); sz = N * nB
    unflat(re, im) = begin
        M = Matrix{ComplexF64}(undef, sz, sz)
        for p in 1:sz, q in 1:sz
            idx = (p - 1) * sz + q          # C [p*sz+q] (0-based) -> 1-based +1
            M[p, q] = complex(re[idx], im[idx])
        end
        M
    end
    return unflat(jvs_re, jvs_im), unflat(jvis_re, jvis_im), N, nB, dimB,
           oiso[1], oiso[2]                 # iso_def, eps_used
end

# ---- partial traces of a (d_maj*d_min)-square matrix, INPUT-major index layout ----
# M index = (s-1)*d_min + i  (s = MAJOR in [1,d_maj], i = MINOR in [1,d_min]).
# ptrace2 = trace the MINOR (sys 2, size d_min) -> d_maj x d_maj (the PINNED dual dir).
function ptrace2(M::Matrix{ComplexF64}, d_maj::Int, d_min::Int)
    R = zeros(ComplexF64, d_maj, d_maj)
    for s in 1:d_maj, t in 1:d_maj, i in 1:d_min
        R[s, t] += M[(s - 1) * d_min + i, (t - 1) * d_min + i]
    end
    return R
end
# ptrace1 = trace the MAJOR (sys 1, size d_maj) -> d_min x d_min (the WRONG direction;
# used only to MEASURE direction teeth via the wrong-direction dual value).
function ptrace1(M::Matrix{ComplexF64}, d_maj::Int, d_min::Int)
    R = zeros(ComplexF64, d_min, d_min)
    for i in 1:d_min, j in 1:d_min, s in 1:d_maj
        R[i, j] += M[(s - 1) * d_min + i, (s - 1) * d_min + j]
    end
    return R
end
maxeig_herm(M) = maximum(real.(eigvals(Hermitian((M + M') / 2))))

# ---- channel constructors (verbatim Kraus of the C test_idemp.h builders) ----
# block_cond_exp(d,m): Phi(X)=P1 X P1 + P2 X P2, P1=proj[0,m), P2=proj[m,d). eta=0.
function block_cond_exp_kraus(d, m)
    K1 = zeros(ComplexF64, d, d); K2 = zeros(ComplexF64, d, d)
    for i in 1:m;     K1[i, i] = 1.0; end
    for i in (m+1):d; K2[i, i] = 1.0; end
    return [K1, K2]
end
# noiseless_subsystem(dL,dE): H=C^dL (x) C^dE, K_{jk}=(1_L (x) |e_j><e_k|)/sqrt(dE),
# H row a*dE+j (left-major). eta=0. N=dL*dE ambient; image ~ B(C^dL).
function noiseless_subsystem_kraus(dL, dE)
    n = dL * dE; inv = 1.0 / sqrt(dE)
    Kr = Matrix{ComplexF64}[]
    for j in 0:(dE-1), k in 0:(dE-1)
        K = zeros(ComplexF64, n, n)
        for a in 0:(dL-1); K[a*dE + j + 1, a*dE + k + 1] = inv; end
        push!(Kr, K)
    end
    return Kr
end
# mixconj(d,m,t): convex mix of compress_idemp(d,m) and its conjugate by a fixed
# complex unitary V=exp(iH) — a genuinely eta>0, oblique, non-commutative channel.
# Kraus union {sqrt(1-t) K_a} U {sqrt(t) V' K_a V}. MUST match test_idemp.h exactly.
function compress_idemp_kraus(d, m)
    r = 1 + (d - m)
    Kr = [zeros(ComplexF64, d, d) for _ in 1:r]
    for i in 1:m; Kr[1][i, i] = 1.0; end
    for b in 0:(d-m-1); Kr[2 + b][1, m + b + 1] = 1.0; end
    return Kr
end
function fixed_unitary(n)   # V = exp(iH), H the test_idemp.h make_fixed_unitary
    H = zeros(ComplexF64, n, n)
    for p in 1:(n-1); H[p, p+1] = 0.4; H[p+1, p] = 0.4; end
    for p in 1:(n-2); H[p, p+2] = 0.3im; H[p+2, p] = -0.3im; end
    for p in 1:n;     H[p, p] = 0.1 * p; end
    return exp(im * Matrix(H))
end
function mixconj_kraus(d, m, t)
    base = compress_idemp_kraus(d, m)
    V = fixed_unitary(d); Vd = V'
    conj = [Vd * K * V for K in base]
    return vcat([sqrt(1 - t) * K for K in base], [sqrt(t) * K for K in conj])
end

# ---- the corpus (start small; all must build) ----
struct Chan
    name::String
    kraus::Vector{Matrix{ComplexF64}}
    n::Int
    eps::Float64
    eta_is_zero::Bool
end
corpus = Chan[
    Chan("block_cond_exp_4_2",     block_cond_exp_kraus(4, 2),         4,  0.0, true),
    Chan("noiseless_subsystem_3_2", noiseless_subsystem_kraus(3, 2),   6,  0.0, true),
    Chan("mixconj_6_2_0p03",       mixconj_kraus(6, 2, 0.03),          6, -1.0, false),
]
# eps convention (matches the shim sentinel, aic_opspace_shim.h):
#   eta=0 oracles  -> eps = 0.0 (exact-idempotent path).
#   mixconj (eta>0) -> eps = -1.0 SENTINEL: the shim DERIVES eps =
#     aic_ecstar_defect_assoc(A) internally, EXACTLY as test_opspace.c T2 does, so
#     the certifier test (rebuilding v with the same derived eps) gets the IDENTICAL
#     v / Choi. The derived value comes back as eps_used and is what we emit.
# If a channel fails to build (e.g. an out-of-basin eta), it is DROPPED with a loud
# note (no fake fixture). The mixconj(6,2,0.03) / (6,3,0.02) builds are confirmed by
# the cstar_build / test_opspace T-series.

# ---- solve both SDPs for one (Choi, dims) direction, with poison guards ----
struct Dir
    tag::String              # "fwd" or "inv"
    J::Matrix{ComplexF64}
    d_maj::Int; d_min::Int
    Y0::Matrix{ComplexF64}; Y1::Matrix{ComplexF64}
    eta_ref::Float64
    teeth::Float64           # |wrong-dir dual - eta_ref| (direction tooth margin)
end
function solve_dir(tag, J, d_maj, d_min, eta_is_zero)
    # PINNED: dual tr_sys=2 (minor), primal rho_on=:major, factor 1.
    ed, Y0, Y1, ds = diamond_dual_rect(J, d_maj, d_min, 2)
    ds in ("OPTIMAL", "SLOW_PROGRESS") || error("$tag dual status $ds")
    ep, _, _, _, ps = diamond_primal_rect(J, d_maj, d_min, :major)
    ps in ("OPTIMAL", "SLOW_PROGRESS") || error("$tag primal status $ps")
    eta_ref = ed
    # strong-duality / consistency guard
    abs(ep - ed) <= 1e-6 ||
        error("FIXTURE POISON $tag: primal=$ep dual=$ed gap=$(abs(ep-ed)) > 1e-6")
    # recompute UPPER from the EMITTED (Y0,Y1) — exactly what the arb certifier does:
    #   hi = (1/2)(lambda_max(Tr2(Y0)) + lambda_max(Tr2(Y1))),  Tr2 = trace MINOR.
    hi = (maxeig_herm(ptrace2(Y0, d_maj, d_min)) +
          maxeig_herm(ptrace2(Y1, d_maj, d_min))) / 2.0
    abs(hi - eta_ref) <= 1e-6 ||
        error("FIXTURE POISON $tag: UPPER from emitted Y = $hi, eta_ref = $eta_ref")
    # eta=0 oracle: a *-iso is a COMPLETE ISOMETRY, ||v||_cb = ||v^-1||_cb = 1 exactly.
    if eta_is_zero
        abs(eta_ref - 1.0) <= 1e-6 ||
            error("FIXTURE POISON $tag: eta=0 oracle eta_ref=$eta_ref != 1.0")
    end
    # DIRECTION TOOTH: the WRONG-direction dual (tr_sys=1, traces MAJOR) should differ
    # from eta_ref so the certifier's tr_sys=2 choice has teeth. Measure the margin.
    edw, _, _, dsw = diamond_dual_rect(J, d_maj, d_min, 1)
    teeth = (dsw in ("OPTIMAL", "SLOW_PROGRESS")) ? abs(edw - eta_ref) : NaN
    return Dir(tag, J, d_maj, d_min, Y0, Y1, eta_ref, teeth)
end

# ---- emit helpers (mirror gen_fixtures_certify.jl style) ----
function flat_rowmajor(M::Matrix{ComplexF64})   # [p*sz + q], sz = size
    sz = size(M, 1)
    re = Float64[]; im = Float64[]
    for p in 1:sz, q in 1:sz
        push!(re, real(M[p, q])); push!(im, imag(M[p, q]))
    end
    return re, im
end
function wr(io, name, v)
    print(io, "static const double $(name)[] = {")
    print(io, join((@sprintf("%.17e", x) for x in v), ", "))
    println(io, "};")
end

struct Solved
    c::Chan
    N::Int; nB::Int; dimB::Int; iso::Float64; eps_used::Float64
    fwd::Dir; inv::Dir
end
function solve_channel(sym, c::Chan)
    n = c.n; r = length(c.kraus)
    Jvs, Jvis, N, nB, dimB, iso, eps_used = shim_choi(sym, c.kraus, n, r, c.eps, CPREC)
    # FORWARD ||v||_cb = ||v*||_⋄: J(v*), dims (d_maj=N, d_min=n_B).
    fwd = solve_dir("$(c.name).fwd", Jvs, N, nB, c.eta_is_zero)
    # INVERSE ||v^-1||_cb = ||(v^-1)*||_⋄: J((v^-1)*), dims (d_maj=n_B, d_min=N).
    inv = solve_dir("$(c.name).inv", Jvis, nB, N, c.eta_is_zero)
    @printf("  %-22s N=%d n_B=%d dim_B=%d eps=%.3e iso=%.2e | fwd eta=%.6f teeth=%.4e | inv eta=%.6f teeth=%.4e\n",
            c.name, N, nB, dimB, eps_used, iso, fwd.eta_ref, fwd.teeth, inv.eta_ref, inv.teeth)
    return Solved(c, N, nB, dimB, iso, eps_used, fwd, inv)
end

function emit(io, solved::Vector{Solved})
    println(io, "/* Auto-generated by tools/gen_fixtures_opspace_o2.jl - DO NOT EDIT.")
    println(io, " * Committed SDP feasible points + reference cb-norms for the opspace O2")
    println(io, " * certified bound on the th_main_ext isomorphism v: B -> A (bead aic-pjr;")
    println(io, " * approximate_algebras.tex:1538-1540). For each channel: Kraus ops (to")
    println(io, " * rebuild v in C), eps_used, the ambient dims (N, n_B, dim_B), and PER")
    println(io, " * DIRECTION (fwd ||v||_cb, inv ||v^-1||_cb) the Convention-A adjoint Choi J,")
    println(io, " * a Watrous MIN-dual feasible point (Y0,Y1), and eta_ref = the certified")
    println(io, " * target cb-norm.")
    println(io, " *")
    println(io, " * THE PINNED CONVENTION (design opspace_o2_design.md §0.5 — factor 1, dual")
    println(io, " * tr_sys = 2 = the MINOR factor): the UPPER bound is recomputed in C from")
    println(io, " *   hi = (1/2)(lambda_max(Tr_2(Y0)) + lambda_max(Tr_2(Y1))),")
    println(io, " * where Tr_2 traces the MINOR factor of dims [d_maj, d_min]:")
    println(io, " *   fwd (J(v*)):      d_maj = N,   d_min = n_B   (v*: M_N -> M_{n_B}).")
    println(io, " *   inv (J((v^-1)*)): d_maj = n_B, d_min = N     ((v^-1)*: M_{n_B} -> M_N).")
    println(io, " * Each J / Y0 / Y1 is (d_maj*d_min) x (d_maj*d_min), flat ROW-MAJOR [p*sz+q].")
    println(io, " * The generator ASSERTS, before writing: primal ~= dual (strong duality),")
    println(io, " * the UPPER recomputed from (Y0,Y1) reproduces eta_ref (~1e-6), eta=0 oracles")
    println(io, " * give eta_ref == 1.0, and the wrong-direction dual differs (teeth). NO PYTHON. */")
    println(io, "#ifndef AIC_FIXTURES_OPSPACE_O2_INC_H")
    println(io, "#define AIC_FIXTURES_OPSPACE_O2_INC_H\n")
    println(io, "typedef struct {")
    println(io, "    const char *name;")
    println(io, "    long n;              /* UCP self-map ambient dim (dim_K == dim_H)        */")
    println(io, "    long r;              /* number of Kraus ops                              */")
    println(io, "    const double *kraus_re;  /* flat [op*n*n + i*n + j] (Convention A)       */")
    println(io, "    const double *kraus_im;")
    println(io, "    double eps_used;     /* eps fed to aic_cstar_build                       */")
    println(io, "    long N;              /* A ambient dim   (= v.n)                          */")
    println(io, "    long n_B;            /* B ambient dim   (= v.B->n_B)                     */")
    println(io, "    long dim_B;          /* dim of A == dim_B (= v.B->dim_B)                 */")
    println(io, "    /* FORWARD ||v||_cb = ||v*||_⋄: J(v*), dims d_maj=N, d_min=n_B; sz=N*n_B */")
    println(io, "    long fwd_dmaj; long fwd_dmin;")
    println(io, "    const double *fwd_J_re;  const double *fwd_J_im;   /* sz x sz, [p*sz+q]  */")
    println(io, "    const double *fwd_Y0_re; const double *fwd_Y0_im;")
    println(io, "    const double *fwd_Y1_re; const double *fwd_Y1_im;")
    println(io, "    double fwd_eta_ref;")
    println(io, "    /* INVERSE ||v^-1||_cb = ||(v^-1)*||_⋄: J((v^-1)*), d_maj=n_B, d_min=N    */")
    println(io, "    long inv_dmaj; long inv_dmin;")
    println(io, "    const double *inv_J_re;  const double *inv_J_im;")
    println(io, "    const double *inv_Y0_re; const double *inv_Y0_im;")
    println(io, "    const double *inv_Y1_re; const double *inv_Y1_im;")
    println(io, "    double inv_eta_ref;")
    println(io, "    int eta_is_zero;     /* 1 for exact-idempotent (complete-isometry) oracles */")
    println(io, "} aic_opspace_o2_fixture;\n")

    for (idx, s) in enumerate(solved)
        c = s.c; n = c.n; r = length(c.kraus)
        kr = Float64[]; ki = Float64[]
        for K in c.kraus, i in 1:n, j in 1:n
            push!(kr, real(K[i, j])); push!(ki, imag(K[i, j]))
        end
        wr(io, "aic_o2_kr_$idx", kr); wr(io, "aic_o2_ki_$idx", ki)
        for (pre, d) in (("fwd", s.fwd), ("inv", s.inv))
            for (tag, M) in (("J", d.J), ("Y0", d.Y0), ("Y1", d.Y1))
                re, im = flat_rowmajor(M)
                wr(io, "aic_o2_$(pre)_$(tag)r_$idx", re)
                wr(io, "aic_o2_$(pre)_$(tag)i_$idx", im)
            end
        end
        println(io)
    end

    println(io, "static const aic_opspace_o2_fixture aic_opspace_o2_fixtures[] = {")
    for (idx, s) in enumerate(solved)
        c = s.c; n = c.n; r = length(c.kraus)
        row  = @sprintf("    {\"%s\", %d, %d, aic_o2_kr_%d, aic_o2_ki_%d, %.17e, ",
                        c.name, n, r, idx, idx, s.eps_used)
        row *= @sprintf("%d, %d, %d, ", s.N, s.nB, s.dimB)
        row *= @sprintf("%d, %d, aic_o2_fwd_Jr_%d, aic_o2_fwd_Ji_%d, ",
                        s.fwd.d_maj, s.fwd.d_min, idx, idx)
        row *= @sprintf("aic_o2_fwd_Y0r_%d, aic_o2_fwd_Y0i_%d, aic_o2_fwd_Y1r_%d, aic_o2_fwd_Y1i_%d, %.17e, ",
                        idx, idx, idx, idx, s.fwd.eta_ref)
        row *= @sprintf("%d, %d, aic_o2_inv_Jr_%d, aic_o2_inv_Ji_%d, ",
                        s.inv.d_maj, s.inv.d_min, idx, idx)
        row *= @sprintf("aic_o2_inv_Y0r_%d, aic_o2_inv_Y0i_%d, aic_o2_inv_Y1r_%d, aic_o2_inv_Y1i_%d, %.17e, ",
                        idx, idx, idx, idx, s.inv.eta_ref)
        row *= @sprintf("%d},", c.eta_is_zero ? 1 : 0)
        println(io, row)
    end
    println(io, "};")
    println(io, "#define AIC_OPSPACE_O2_NFIX ((int)(sizeof(aic_opspace_o2_fixtures)/sizeof(aic_opspace_o2_fixtures[0])))\n")
    println(io, "#endif /* AIC_FIXTURES_OPSPACE_O2_INC_H */")
end

# ---- main ----
println("loading build/libaic.so + solving the PINNED O2 SDP (factor 1, dual tr_sys=2):")
sym = load_libaic()
solved = Solved[]
for c in corpus
    try
        push!(solved, solve_channel(sym, c))
    catch err
        @warn "channel $(c.name) failed to build/solve — DROPPED (no fake fixture)" exception=err
    end
end
isempty(solved) && error("no channels survived — refusing to emit an empty fixture")

# ---- DIRECTION-TOOTH guard (the load-bearing partial-trace-direction check) ----
# The certifier's tr_sys=2 (MINOR) choice must be discriminating on at least ONE
# committed fixture, else swapping it could go unnoticed. The wrong-direction dual
# (tr_sys=1) differs from eta_ref by `teeth`; we require some direction on some
# fixture to have teeth > 1e-3. The eta>0 mixconj fixture is the intended carrier
# (its J(v*) / J((v^-1)*) marginals are asymmetric); NOTE loudly if it is inert.
const TEETH_MIN = 1e-3
all_teeth = vcat([[s.fwd.teeth, s.inv.teeth] for s in solved]...)
max_teeth = maximum(filter(!isnan, all_teeth); init = 0.0)
for s in solved
    if !s.c.eta_is_zero
        t = max(s.fwd.teeth, s.inv.teeth)
        if t <= TEETH_MIN
            @warn "DIRECTION TOOTH INERT on eta>0 fixture $(s.c.name): max teeth " *
                  "= $t <= $TEETH_MIN — the tr_sys-direction swap would be SILENT " *
                  "on it. The corpus tooth must live elsewhere."
        else
            @printf("  direction tooth on eta>0 fixture %s: teeth=%.6f (fwd=%.4e inv=%.4e)\n",
                    s.c.name, t, s.fwd.teeth, s.inv.teeth)
        end
    end
end
max_teeth > TEETH_MIN ||
    error("DIRECTION TOOTH POISON: NO committed fixture has teeth > $TEETH_MIN " *
          "(max over corpus = $max_teeth). The certifier's tr_sys=2 direction choice " *
          "would be untested — add an asymmetric eta>0 fixture.")
@printf("direction-tooth guard PASSED: max teeth over corpus = %.6f (> %.1e)\n",
        max_teeth, TEETH_MIN)

out = joinpath(@__DIR__, "..", "tests", "fixtures_opspace_o2.inc.h")
open(out, "w") do io
    emit(io, solved)
end
println("wrote ", normpath(out), "  (", length(solved), " channels)")
