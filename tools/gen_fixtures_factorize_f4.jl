# ============================================================================
# STATUS (2026-05-31): F4.2 DEFERRED TO v0.2 (bead aic-bag). This generator is
# SCAFFOLDING, not a working committed-fixture producer yet. The PRIMAL self-map
# diamond SDP converges (OPTIMAL = the exact cb-norm), but the DUAL stalls
# (SLOW_PROGRESS) at n>=6 in Convex.jl, so the rigorous two-sided certificate is
# blocked at the sweep dims and the strong-duality poison guard rightly refuses to
# emit (no fake fixture). See paper/FINDINGS.md D6. The headline th_factorization
# is CLOSED on F4.1 (aic-tff). v0.2 routes: reformulate the dual / move to direct
# JuMP+MOSEK / primal-only canary. The eager-flush logging, F4_ONLY filter,
# relaxed-tolerance _mosek_tight override, and GC below are the v0.2 starting point.
# ============================================================================
# gen_fixtures_factorize_f4.jl — committed-fixture generator for the factorize F4.2
# TIGHT cb-norm certification of the th_factorization headline bounds (bead aic-tff,
# increment F4.2; approximate_algebras.tex:2730-2739; design
# docs/research/factorize_f4_design.md §C.3/§C.4). A CLOSE MIRROR of
# tools/gen_fixtures_opspace_o2.jl — same dlopen/ccall/marshalling/poison-guard/emit
# structure; the differences from O2 are flagged inline.
#
# WHAT IT EMITS: tests/fixtures_factorize_f4.inc.h (committed, generated). For each
# corpus channel the shim aic_factorize_choi_shim_d rebuilds the FULL factorize
# pipeline and returns TWO square Convention-A Choi DIFFERENCE matrices:
#   J_DelUps = Choi(Delta Upsilon) - Choi(Phi)   : N^2 x N^2   (DelUps, .tex:2733)
#   J_UpsDel = Choi(Upsilon Delta - 1_B)         : n_B^2 x n_B^2 (UpsDel2, .tex:2739)
# Both are SQUARE HP self-maps (on M_N and M_{n_B} respectively), so we use the
# SQUARE self-map Watrous SDP (the (2/n) normalization, NOT the rect factor-1 one
# the O2 generator uses). For EACH square J we solve BOTH:
#   eta_primal, X,P,Q = diamond_norm_watrous_primal(J,n)  -> MAX primal LOWER feasible
#   eta_dual,   Y0,Y1 = diamond_norm_dual(J,n)            -> MIN dual  UPPER feasible
# and emit both feasible points + eta_ref = eta_dual (the certified upper). The arb
# certifier (the C consumer, test_factorize_f4.c) restores these to exact feasibility
# so `make test` stays Julia-free.
#
# CRITICAL DIFFERENCE FROM O2 (eta=0 oracle): O2's eta=0 oracle is a COMPLETE
# ISOMETRY, so ||v||_cb = 1 exactly. HERE the eta=0 oracle is the FACTORIZATION
# DEFECT Delta Upsilon - Phi (resp. Upsilon Delta - 1_B), which is ZERO at eta=0
# (ladder rung 3, design §B.1) — so the oracle guard asserts eta_ref < 1e-6, NOT
# eta_ref == 1.0. (Loud note below.)
#
# POISON GUARDS before emitting ANY fixture (fail-loud; mirror
# gen_fixtures_opspace_o2.jl:195-213 + gen_fixtures_certify.jl square-case helpers):
# strong duality (primal ~= dual); recompute the UPPER from the emitted (Y0,Y1) and
# the LOWER from the emitted (X,P,Q), both must reproduce eta_ref to ~1e-6; for eta=0
# oracles assert eta_ref < 1e-6. A channel that fails ANY guard (or fails to
# build/solve) is DROPPED with a loud note — never a fake fixture.
#
# NO PYTHON. Determinism. Run serially (single Julia process), AFTER `make lib`:
#   make lib && julia --project=julia/env tools/gen_fixtures_factorize_f4.jl
using Convex, MosekTools, LinearAlgebra, Printf, Libdl

include(joinpath(@__DIR__, "..", "julia", "AlmostIdempotentChannels.jl",
                 "src", "sdp.jl"))   # diamond_norm_watrous_primal / diamond_norm_dual

# --- F4.2 SDP solver tuning (generator-LOCAL override; sdp.jl is NOT edited) ---
# Julia resolves the `_mosek_tight` global at CALL time, so redefining it here makes
# the diamond_norm_* functions (also in Main, via include) pick up this RELAXED
# factory. WHY: sdp.jl's 1e-12 feas/gap + 1e-14 MU_RED make Convex.jl/MOSEK GRIND
# (SLOW_PROGRESS) ~100s and ~20GB at n=6 and OOM at n=7 for these self-map diamond
# SDPs — the 1e-14 central-path target is below MOSEK's numerical floor so it stalls.
# Relaxing MU_RED to 1e-10 (a REACHABLE target) lets it converge to OPTIMAL fast at
# bounded memory. Tolerances 1e-9 keep STRONG DUALITY clean (the primal==dual poison
# guard, 1e-6) — do NOT use a MAX_ITERATIONS cap: it truncated the DUAL minimization
# short (primal 0.0130 vs truncated dual 0.0237, an 83% gap the guard rightly caught).
# The canary itself needs only ~1e-4 (design §F.3 T9), so 1e-9 is ample. The committed
# O2 / cbnorm certifiers (separate runs) keep sdp.jl's tight tolerances unchanged.
function _mosek_tight()
    o = MosekTools.Mosek.Optimizer()
    for (k, v) in ("MSK_DPAR_INTPNT_CO_TOL_DFEAS"   => 1.0e-9,
                   "MSK_DPAR_INTPNT_CO_TOL_PFEAS"   => 1.0e-9,
                   "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1.0e-9,
                   "MSK_DPAR_INTPNT_CO_TOL_MU_RED"  => 1.0e-10)
        MosekTools.MOI.set(o, MosekTools.MOI.RawOptimizerAttribute(k), v)
    end
    return o
end

const LIBAIC = abspath(joinpath(@__DIR__, "..", "build", "libaic.so"))
const CPREC = 256          # matches test_factorize CPREC

# Eager-flush logger (CLAUDE.md: verbose eager flush + fail fast/loud). Julia
# block-buffers stdout to a redirected (non-TTY) file, so a crash/hang shows NO
# output unless we flush after every progress line. Use logf() for ALL progress.
logf(msg) = (println(msg); flush(stdout))

# ---- dlopen libaic (mirror the package __init__: system libgmp first, DEEPBIND) ----
# (Verbatim from gen_fixtures_opspace_o2.jl:36-49, only the dlsym symbol differs.)
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
    return Libdl.dlsym(h, :aic_factorize_choi_shim_d)
end

# ---- ccall the shim: rebuild the factorize pipeline, get J_DelUps, J_UpsDel ----
# Returns (Jdu, Jud, N, nB, dimB, eta, eps_used) with Jdu an (N*N) x (N*N) and Jud
# an (n_B*n_B) x (n_B*n_B) ComplexF64 matrix, unflattened from the C row-major
# [p*sz+q] arrays (sz_du = N*N, sz_ud = n_B*n_B). The ccall tuple/marshalling MIRRORS
# gen_fixtures_opspace_o2.jl:54-89 EXACTLY (same buffer lengths n^4, same Kraus
# layout, same Ref{Cint}/length-2 out buffer; only the two output Choi have DIFFERENT
# sizes — sz_du=N^2 vs sz_ud=n_B^2 — so they are unflattened with their own sz).
function shim_choi(sym, kraus::Vector{Matrix{ComplexF64}}, n::Int, r::Int,
                   eps::Float64, prec::Int)
    # Flatten Kraus: kraus_re[(a-1)*n*n + (i-1)*n + j] = K_a[i,j] (Convention A, K-major).
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
    jdu_re = Vector{Float64}(undef, n4); jdu_im = Vector{Float64}(undef, n4)
    jud_re = Vector{Float64}(undef, n4); jud_im = Vector{Float64}(undef, n4)
    oN = Ref{Cint}(0); onB = Ref{Cint}(0); odimB = Ref{Cint}(0)
    oeta = Vector{Float64}(undef, 2)        # [0]=eta_proxy, [1]=eps actually used
    rc = ccall(sym, Cint,
               (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
                Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cdouble, Cint),
               jdu_re, jdu_im, jud_re, jud_im, oN, onB, odimB, oeta,
               kr, ki, n, r, eps, prec)
    rc == 0 || error("aic_factorize_choi_shim_d returned $rc")
    N = Int(oN[]); nB = Int(onB[]); dimB = Int(odimB[])
    sz_du = N * N; sz_ud = nB * nB
    unflat(re, im, sz) = begin
        M = Matrix{ComplexF64}(undef, sz, sz)
        for p in 1:sz, q in 1:sz
            idx = (p - 1) * sz + q          # C [p*sz+q] (0-based) -> 1-based +1
            M[p, q] = complex(re[idx], im[idx])
        end
        M
    end
    return unflat(jdu_re, jdu_im, sz_du), unflat(jud_re, jud_im, sz_ud),
           N, nB, dimB, oeta[1], oeta[2]     # eta_proxy, eps_used
end

# ---- SQUARE self-map partial trace + helpers (verbatim from gen_fixtures_certify.jl:
#      109-128). Tr_2(M)[i,j] = sum_a M[(i-1)n+a, (j-1)n+a] traces the MINOR factor
#      (sys 2 of [n,n]), matching the SDP's partialtrace(Y,2,[n,n]). For the square
#      self-map BOTH factors are size n. ----
function ptrace2(M::Matrix{ComplexF64}, n::Int)
    R = zeros(ComplexF64, n, n)
    for i in 1:n, j in 1:n, a in 1:n
        R[i, j] += M[(i - 1) * n + a, (j - 1) * n + a]
    end
    return R
end
maxeig_herm(M) = maximum(real.(eigvals(Hermitian((M + M') / 2))))

# ---- channel constructors (Kraus built IN JULIA; reuse the o2 generator's
#      verbatim-from-C builders). block_cond_exp / fixed_unitary are identical to the
#      o2 generator's; mixconj_blocks is the block-base analogue of o2's mixconj
#      (compress-base) — it MIRRORS test_cstar_build.c:170-189 make_mixconj_blocks. ----

# block_cond_exp(d,m): Phi(X)=P1 X P1 + P2 X P2, P1=proj[0,m), P2=proj[m,d). eta=0,
# >=2 genuine matrix blocks (dim_A = m^2+(d-m)^2). Verbatim from
# gen_fixtures_opspace_o2.jl:114-119 == C make_block_cond_exp (test_idemp.h:30-35).
function block_cond_exp_kraus(d, m)
    K1 = zeros(ComplexF64, d, d); K2 = zeros(ComplexF64, d, d)
    for i in 1:m;     K1[i, i] = 1.0; end
    for i in (m+1):d; K2[i, i] = 1.0; end
    return [K1, K2]
end

# fixed_unitary(n) = exp(iH), H the test_idemp.h make_fixed_unitary (off-diagonal
# real 0.4 super-diagonal + 0.3i second-super-diagonal, phases 0.1*p on the
# diagonal). Verbatim from gen_fixtures_opspace_o2.jl:142-148 == C make_fixed_unitary
# (test_idemp.h:106-133). The 1-based loops here reproduce the C 0-based (p+1)
# diagonal phase exactly.
function fixed_unitary(n)
    H = zeros(ComplexF64, n, n)
    for p in 1:(n-1); H[p, p+1] = 0.4; H[p+1, p] = 0.4; end
    for p in 1:(n-2); H[p, p+2] = 0.3im; H[p+2, p] = -0.3im; end
    for p in 1:n;     H[p, p] = 0.1 * p; end
    return exp(im * Matrix(H))
end

# mixconj_blocks(d,m,t): convex mix of block_cond_exp(d,m) [TWO genuine matrix
# blocks => >=2 equivalence classes at eta>0] and its conjugate by the fixed unitary
# V=exp(iH). Kraus union {sqrt(1-t) K_a} U {sqrt(t) V' K_a V}, r=4. MIRRORS
# test_cstar_build.c:170-189 make_mixconj_blocks EXACTLY (base = make_block_cond_exp,
# conj = make_conjugated via K'_a = V^dag K_a V; out->K = [sqrt(1-t)*base ;
# sqrt(t)*conj]). This is the o2 mixconj_kraus recipe with the BLOCK base swapped in
# for the compress base (the §C13(c) multi-block eta>0 fixture). eta ~ 0.02 (<1/4).
function mixconj_blocks_kraus(d, m, t)
    base = block_cond_exp_kraus(d, m)
    V = fixed_unitary(d); Vd = V'
    conj = [Vd * K * V for K in base]           # K'_a = V^dag K_a V  (make_conjugated)
    return vcat([sqrt(1 - t) * K for K in base], [sqrt(t) * K for K in conj])
end

# ---- the corpus ----
# eps convention (the shim sentinel, aic_factorize_shim.h:78-89):
#   eta=0 ORACLES        -> eps = 0.0  (exact-idempotent path; eps_used = 0.0).
#   eta>0 multi-block     -> eps = -2.0 SENTINEL: the shim sets eps_used = eta_proxy
#     (the §C11/§C13(c) caveat — pass eta, NOT the ~700x smaller eps_assoc, else the
#     Stage-1 errreduce C0 gate fires and the build aborts). The derived eps_used
#     comes back in out_eta[1] and is what we emit.
struct Chan
    name::String
    kraus::Vector{Matrix{ComplexF64}}
    n::Int
    eps::Float64
    eta_is_zero::Bool
    is_canary::Bool
end
corpus = Chan[
    # (A) eta=0 ORACLES (multi-block, >=2 blocks, eta_ref must be ~0 for BOTH J).
    Chan("block_cond_exp_4_2", block_cond_exp_kraus(4, 2), 4, 0.0, true,  false),
    Chan("block_cond_exp_6_3", block_cond_exp_kraus(6, 3), 6, 0.0, true,  false),
    # (B) eta>0 CANARY SWEEP (multi-block, INCREASING dim so dim_B grows; >=4 distinct
    #     dim_B values for the halves-ratio to be robust, §D2). All at t=0.03.
    Chan("mixconj_blocks_4_2_0p03", mixconj_blocks_kraus(4, 2, 0.03), 4, -2.0, false, true),
    Chan("mixconj_blocks_5_2_0p03", mixconj_blocks_kraus(5, 2, 0.03), 5, -2.0, false, true),
    Chan("mixconj_blocks_6_2_0p03", mixconj_blocks_kraus(6, 2, 0.03), 6, -2.0, false, true),
    Chan("mixconj_blocks_6_3_0p03", mixconj_blocks_kraus(6, 3, 0.03), 6, -2.0, false, true),
    Chan("mixconj_blocks_7_3_0p03", mixconj_blocks_kraus(7, 3, 0.03), 7, -2.0, false, true),
]

# ---- solve the SQUARE self-map SDP for one Choi (tag "du" or "ud"), poison guards ----
# Mirrors gen_fixtures_opspace_o2.jl solve_dir (195-213), adapted to the SQUARE
# self-map programs (diamond_norm_watrous_primal / diamond_norm_dual, (2/n) norm)
# instead of the rect ones, and emitting (X,P,Q) AND (Y0,Y1) (the self-map SDP
# carries both bounds). n = the side of the M_n the map acts on (N for du, n_B for ud).
struct Quant
    tag::String                 # "du" or "ud"
    n::Int                      # side: M_n self-map (J is n^2 x n^2)
    J::Matrix{ComplexF64}
    X::Matrix{ComplexF64}; P::Matrix{ComplexF64}; Q::Matrix{ComplexF64}
    Y0::Matrix{ComplexF64}; Y1::Matrix{ComplexF64}
    eta_ref::Float64
end
function solve_quant(tag, J, n, eta_is_zero)
    size(J, 1) == n * n && size(J, 2) == n * n ||
        error("$tag: J shape $(size(J)) != (n^2, n^2) = ($(n*n), $(n*n))")
    # MAX primal (LOWER feasible (X,P,Q)) and MIN dual (UPPER feasible (Y0,Y1)).
    local eta_primal, X, P, Q, ps, eta_dual, Y0, Y1, ds
    logf("      [$tag] primal SDP n=$n J=$(size(J,1))^2 solving...")
    tp = @elapsed ((eta_primal, X, P, Q, ps) = diamond_norm_watrous_primal(J, n))
    logf("      [$tag] primal done $(round(tp,digits=2))s eta=$(round(eta_primal,sigdigits=6)) status=$ps")
    ps in ("OPTIMAL", "SLOW_PROGRESS") || error("$tag primal status $ps")
    logf("      [$tag] dual SDP solving...")
    td = @elapsed ((eta_dual, Y0, Y1, ds) = diamond_norm_dual(J, n))
    logf("      [$tag] dual done $(round(td,digits=2))s eta=$(round(eta_dual,sigdigits=6)) status=$ds")
    ds in ("OPTIMAL", "SLOW_PROGRESS") || error("$tag dual status $ds")
    eta_ref = eta_dual                      # the certified UPPER is the reference
    # GUARD 1 — strong duality: primal == dual.
    abs(eta_primal - eta_dual) <= 1e-6 ||
        error("FIXTURE POISON $tag: primal=$eta_primal dual=$eta_dual " *
              "gap=$(abs(eta_primal-eta_dual)) > 1e-6")
    # GUARD 2 — recompute UPPER from emitted (Y0,Y1) (what the arb certifier does):
    #   hi = (1/2)(lambda_max(Tr_2(Y0)) + lambda_max(Tr_2(Y1))),  Tr_2 = trace MINOR.
    #   (diamond_norm_dual returns optval/2 with optval = ||Tr2 Y0|| + ||Tr2 Y1||.)
    hi = (maxeig_herm(ptrace2(Y0, n)) + maxeig_herm(ptrace2(Y1, n))) / 2.0
    abs(hi - eta_ref) <= 1e-6 ||
        error("FIXTURE POISON $tag: UPPER from emitted Y = $hi, eta_ref = $eta_ref")
    # GUARD 3 — recompute LOWER from emitted (X,P,Q):  lo = (2/n) Re tr(J' X).
    lo = (2.0 / n) * real(LinearAlgebra.tr(J' * X))
    abs(lo - eta_ref) <= 1e-6 ||
        error("FIXTURE POISON $tag: LOWER from emitted (X,P,Q) = $lo, eta_ref = $eta_ref")
    # GUARD 4 — eta=0 ORACLE (CRITICAL DIFFERENCE FROM O2): the factorization defect
    #   Delta Upsilon - Phi (resp. Upsilon Delta - 1_B) is ZERO at eta=0, so eta_ref
    #   ~ 0 here — NOT 1.0 like O2's complete-isometry oracle (||v||_cb = 1). The O2
    #   oracle measures a *-iso (a complete isometry); THIS oracle measures a DEFECT
    #   that vanishes. See design §B.1 / §C.3.
    if eta_is_zero
        eta_ref < 1e-6 ||
            error("FIXTURE POISON $tag: eta=0 oracle eta_ref=$eta_ref >= 1e-6 " *
                  "(the factorization defect must VANISH at eta=0)")
    end
    return Quant(tag, n, J, X, P, Q, Y0, Y1, eta_ref)
end

# ---- emit helpers (mirror gen_fixtures_opspace_o2.jl:217-229) ----
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
    N::Int; nB::Int; dimB::Int; eta::Float64; eps_used::Float64
    du::Quant; ud::Quant
end
function solve_channel(sym, c::Chan)
    n = c.n; r = length(c.kraus)
    logf("[$(c.name)] ccall shim (n=$n r=$r eps=$(c.eps) prec=$CPREC) rebuilding pipeline...")
    local Jdu, Jud, N, nB, dimB, eta, eps_used
    tc = @elapsed ((Jdu, Jud, N, nB, dimB, eta, eps_used) =
                   shim_choi(sym, c.kraus, n, r, c.eps, CPREC))
    logf("[$(c.name)] shim done $(round(tc,digits=2))s: N=$N n_B=$nB dim_B=$dimB " *
         "eta=$(round(eta,sigdigits=6)) ||Jdu||_F=$(round(norm(Jdu),sigdigits=4)) " *
         "||Jud||_F=$(round(norm(Jud),sigdigits=4))")
    # DelUps: square self-map on M_N (J_DelUps is N^2 x N^2).
    du = solve_quant("$(c.name).du", Jdu, N, c.eta_is_zero)
    GC.gc()                                   # free the Convex problem before the next SDP
    # UpsDel: square self-map on M_{n_B} (J_UpsDel is n_B^2 x n_B^2).
    ud = solve_quant("$(c.name).ud", Jud, nB, c.eta_is_zero)
    GC.gc()
    @printf("  %-26s N=%d n_B=%d dim_B=%d eta=%.6e eps_used=%.6e | du_eta_ref=%.6e ud_eta_ref=%.6e\n",
            c.name, N, nB, dimB, eta, eps_used, du.eta_ref, ud.eta_ref)
    return Solved(c, N, nB, dimB, eta, eps_used, du, ud)
end

function emit(io, solved::Vector{Solved})
    println(io, "/* Auto-generated by tools/gen_fixtures_factorize_f4.jl - DO NOT EDIT.")
    println(io, " * Committed TIGHT-SDP feasible points + reference cb-norms for the")
    println(io, " * factorize F4.2 certification of the th_factorization headline bounds")
    println(io, " * (bead aic-tff; approximate_algebras.tex:2730-2739). For each channel:")
    println(io, " * Kraus ops (to rebuild the WHOLE factorize pipeline in C), eps_used, the")
    println(io, " * dims (N, n_B, dim_B), eta, and per QUANTITY (du = DelUps ||Delta Upsilon")
    println(io, " * - Phi||_cb on M_N, ud = UpsDel ||Upsilon Delta - 1_B||_cb on M_{n_B}) the")
    println(io, " * square Convention-A Choi DIFFERENCE J, a Watrous MAX-primal feasible point")
    println(io, " * (X,P,Q) (LOWER), a MIN-dual feasible point (Y0,Y1) (UPPER), and eta_ref =")
    println(io, " * the certified target cb-norm (= the dual optval).")
    println(io, " *")
    println(io, " * THE CONVENTION (inherited; design factorize_f4_design.md §E, §C.1): both")
    println(io, " * J_DelUps and J_UpsDel are SQUARE HP self-maps (on M_N, M_{n_B}), so the")
    println(io, " * SQUARE self-map (2/n) normalization applies (NOT the O2 rect factor-1):")
    println(io, " *   LOWER  lo = (2/n) Re tr(J' X)                 (from the MAX primal).")
    println(io, " *   UPPER  hi = (1/2)(lambda_max(Tr_2(Y0)) + lambda_max(Tr_2(Y1)))")
    println(io, " *             Tr_2 traces the MINOR factor of [n,n]   (from the MIN dual).")
    println(io, " *   du: n = N (J is N^2 x N^2);  ud: n = n_B (J is n_B^2 x n_B^2).")
    println(io, " *")
    println(io, " * eta_is_zero ORACLES (block_cond_exp): the factorization defect VANISHES at")
    println(io, " * eta=0 (Delta Upsilon = Phi, Upsilon Delta = 1_B), so eta_ref ~ 0 — UNLIKE")
    println(io, " * the O2 complete-isometry oracle (||v||_cb = 1). is_canary fixtures")
    println(io, " * (mixconj_blocks) form the INCREASING-dim sweep for the dimension-")
    println(io, " * independence canary (C = eta_ref/eta vs dim_B; design §C.4).")
    println(io, " *")
    println(io, " * The generator ASSERTS, before writing: primal ~= dual (strong duality), the")
    println(io, " * UPPER from (Y0,Y1) and LOWER from (X,P,Q) both reproduce eta_ref (~1e-6),")
    println(io, " * and eta=0 oracles give eta_ref < 1e-6. NO PYTHON.")
    println(io, " *")
    println(io, " * Regenerate:  make lib && julia --project=julia/env tools/gen_fixtures_factorize_f4.jl */")
    println(io, "#ifndef AIC_FIXTURES_FACTORIZE_F4_INC_H")
    println(io, "#define AIC_FIXTURES_FACTORIZE_F4_INC_H\n")
    println(io, "typedef struct {")
    println(io, "    const char *name;")
    println(io, "    long n; long r;")
    println(io, "    const double *kraus_re; const double *kraus_im;   /* flat [a*n*n+i*n+j] */")
    println(io, "    double eps_used;                                   /* positive; rebuild with this */")
    println(io, "    long N; long n_B; long dim_B;")
    println(io, "    double eta;                                         /* eta_proxy */")
    println(io, "    int eta_is_zero;                                    /* 1 => oracle */")
    println(io, "    int is_canary;                                      /* 1 => in the dim sweep */")
    println(io, "    /* DelUps: square self-map on M_N, sz=N*N */")
    println(io, "    long du_sz;")
    println(io, "    const double *du_J_re;  const double *du_J_im;")
    println(io, "    const double *du_X_re;  const double *du_X_im;")
    println(io, "    const double *du_P_re;  const double *du_P_im;")
    println(io, "    const double *du_Q_re;  const double *du_Q_im;")
    println(io, "    const double *du_Y0_re; const double *du_Y0_im;")
    println(io, "    const double *du_Y1_re; const double *du_Y1_im;")
    println(io, "    double du_eta_ref;")
    println(io, "    /* UpsDel: square self-map on M_{n_B}, sz=n_B*n_B */")
    println(io, "    long ud_sz;")
    println(io, "    const double *ud_J_re;  const double *ud_J_im;")
    println(io, "    const double *ud_X_re;  const double *ud_X_im;")
    println(io, "    const double *ud_P_re;  const double *ud_P_im;")
    println(io, "    const double *ud_Q_re;  const double *ud_Q_im;")
    println(io, "    const double *ud_Y0_re; const double *ud_Y0_im;")
    println(io, "    const double *ud_Y1_re; const double *ud_Y1_im;")
    println(io, "    double ud_eta_ref;")
    println(io, "} aic_factorize_f4_fixture;\n")

    for (idx, s) in enumerate(solved)
        c = s.c; n = c.n
        kr = Float64[]; ki = Float64[]
        for K in c.kraus, i in 1:n, j in 1:n
            push!(kr, real(K[i, j])); push!(ki, imag(K[i, j]))
        end
        wr(io, "aic_f4_kr_$idx", kr); wr(io, "aic_f4_ki_$idx", ki)
        for (pre, qd) in (("du", s.du), ("ud", s.ud))
            for (tag, M) in (("J", qd.J), ("X", qd.X), ("P", qd.P), ("Q", qd.Q),
                             ("Y0", qd.Y0), ("Y1", qd.Y1))
                re, im = flat_rowmajor(M)
                wr(io, "aic_f4_$(pre)_$(tag)r_$idx", re)
                wr(io, "aic_f4_$(pre)_$(tag)i_$idx", im)
            end
        end
        println(io)
    end

    println(io, "static const aic_factorize_f4_fixture aic_f4_fixtures[] = {")
    for (idx, s) in enumerate(solved)
        c = s.c; n = c.n; r = length(c.kraus)
        row  = @sprintf("    {\"%s\", %d, %d, aic_f4_kr_%d, aic_f4_ki_%d, %.17e, ",
                        c.name, n, r, idx, idx, s.eps_used)
        row *= @sprintf("%d, %d, %d, %.17e, %d, %d, ",
                        s.N, s.nB, s.dimB, s.eta,
                        c.eta_is_zero ? 1 : 0, c.is_canary ? 1 : 0)
        # DelUps block (sz = N*N).
        row *= @sprintf("%d, aic_f4_du_Jr_%d, aic_f4_du_Ji_%d, ", s.du.n * s.du.n, idx, idx)
        row *= @sprintf("aic_f4_du_Xr_%d, aic_f4_du_Xi_%d, aic_f4_du_Pr_%d, aic_f4_du_Pi_%d, ",
                        idx, idx, idx, idx)
        row *= @sprintf("aic_f4_du_Qr_%d, aic_f4_du_Qi_%d, aic_f4_du_Y0r_%d, aic_f4_du_Y0i_%d, ",
                        idx, idx, idx, idx)
        row *= @sprintf("aic_f4_du_Y1r_%d, aic_f4_du_Y1i_%d, %.17e, ", idx, idx, s.du.eta_ref)
        # UpsDel block (sz = n_B*n_B).
        row *= @sprintf("%d, aic_f4_ud_Jr_%d, aic_f4_ud_Ji_%d, ", s.ud.n * s.ud.n, idx, idx)
        row *= @sprintf("aic_f4_ud_Xr_%d, aic_f4_ud_Xi_%d, aic_f4_ud_Pr_%d, aic_f4_ud_Pi_%d, ",
                        idx, idx, idx, idx)
        row *= @sprintf("aic_f4_ud_Qr_%d, aic_f4_ud_Qi_%d, aic_f4_ud_Y0r_%d, aic_f4_ud_Y0i_%d, ",
                        idx, idx, idx, idx)
        row *= @sprintf("aic_f4_ud_Y1r_%d, aic_f4_ud_Y1i_%d, %.17e},", idx, idx, s.ud.eta_ref)
        println(io, row)
    end
    println(io, "};")
    println(io, "#define AIC_F4_NFIX ((int)(sizeof(aic_f4_fixtures)/sizeof(aic_f4_fixtures[0])))\n")
    println(io, "#endif /* AIC_FIXTURES_FACTORIZE_F4_INC_H */")
end

# ---- main ----
println("loading build/libaic.so + solving the SQUARE self-map factorize SDP ((2/n) norm):")
println("  NOTE: eta=0 oracle gives eta_ref ~ 0 (factorization defect VANISHES) — NOT 1.0")
println("        like the O2 complete-isometry oracle. Guard 4 asserts eta_ref < 1e-6.")
flush(stdout)
sym = load_libaic()
# Diagnostic / incremental filter: F4_ONLY=<substr> runs only matching channels.
only = get(ENV, "F4_ONLY", "")
run_corpus = isempty(only) ? corpus : filter(c -> occursin(only, c.name), corpus)
logf("processing $(length(run_corpus))/$(length(corpus)) channels" *
     (isempty(only) ? "" : " (F4_ONLY=\"$only\")"))
solved = Solved[]
for c in run_corpus
    try
        push!(solved, solve_channel(sym, c))
    catch err
        @warn "channel $(c.name) failed to build/solve — DROPPED (no fake fixture)" exception=err
        flush(stderr)
    end
end
isempty(solved) && error("no channels survived — refusing to emit an empty fixture")

# ---- canary-sweep sanity: report the (dim_B, C=du_eta_ref/eta) trend the C consumer
#      will assert (design §C.4). Generator does NOT gate on it (the C canary does);
#      it just SURFACES the sweep so a poisoned shim/SDP is visible at generation. ----
canary = [s for s in solved if s.c.is_canary]
distinct_dimB = unique(sort([s.dimB for s in canary]))
if !isempty(canary)
    println("canary sweep (C = du_eta_ref / eta vs dim_B; design §C.4):")
    for s in sort(canary; by = x -> x.dimB)
        C = s.eta > 1e-12 ? s.du.eta_ref / s.eta : NaN
        @printf("  %-26s dim_B=%-3d eta=%.4e du_eta_ref=%.4e  C=%.4f\n",
                s.c.name, s.dimB, s.eta, s.du.eta_ref, C)
    end
    length(distinct_dimB) >= 4 ||
        @warn "canary has only $(length(distinct_dimB)) distinct dim_B values " *
              "(< 4) — the halves-ratio (§D2) may be fragile; consider widening the sweep."
end

out = joinpath(@__DIR__, "..", "tests", "fixtures_factorize_f4.inc.h")
open(out, "w") do io
    emit(io, solved)
end
@printf("wrote %s  (%d channels: %d oracle, %d canary; %d distinct canary dim_B)\n",
        normpath(out), length(solved),
        count(s -> s.c.eta_is_zero, solved), length(canary), length(distinct_dimB))
