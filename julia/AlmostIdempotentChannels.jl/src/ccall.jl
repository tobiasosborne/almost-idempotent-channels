# ccall.jl — low-level @ccall wrappers over the NEW C ABI shims C2/C3 (bead
# aic-exa.6 [J3]). C4/C5 (the two-call factorize-artifacts protocol) live in the
# sibling ccall_factorize.jl to keep each file under the 200-LOC house limit.
#
# These wrappers are INTERNAL (underscore-prefixed, NOT exported): bead J4 (api.jl)
# builds the public verbs (associated_algebra, main_isomorphism, factorize) on top.
# Each uses the dlsym-handle 4-arg ccall form (research §1e: @ccall cannot take a
# raw Ptr{Cvoid}; the dlsym handle in libaic.jl is what combines the load-bearing
# RTLD flags with a type-safe ccall). GC.@preserve guards every buffer; fail loud
# (Rule 4) on rc != 0 (naming the C function + rc) and on inconsistent dims.
#
# MARSHALLING (must match the shims byte-for-byte; layout bugs are this project's
# worst historical defect class):
#   * Kraus IN: K-major [a*n*n + i*n + j] via _flatten_kraus (main module).
#   * basis OUT (C2): dim_A ops, each n x n, K-major [k*n*n + i*n + j] -> _unflatten_kraus.
#   * brackets arrive as length-2 [lo,hi] doubles, already FLOOR/CEIL-rigorous from
#     C; wrapped in CertifiedBracket (src/types.jl) with the project Unicode label.
#   * out_eps/out_eta carry eps_used (read back, NOT the sentinel).
#
# Design: docs/research/julia_package_design.md §4.3 (C2/C3 signatures), §6.1
# (the dlsym-handle ccall form), §9 [J3] (the step contract); headers
# include/aic/aic_assoc_shim.h, aic_iso_shim.h. The eps SENTINEL (§4.3): eps==0.0
# eta=0 oracle; eps==-2.0 multi-block eta>0; eps<0 (!=-2) derive; eps>0 verbatim.
# Default eps=-2.0 (multiblock) here; J4 chooses the sentinel per call.

# ----- unflatten helpers (flat C K-major arrays -> Julia ComplexF64 matrices) -----

# Inverse of _flatten_kraus for a set of `count` matrices each `row` x `col`, laid
# out K-major [a*(row*col) + i*col + j] (0-based C; 1-based Julia). Returns a
# Vector{Matrix{ComplexF64}}. Used for the C2 n x n basis (row=col=n) and the C5
# dec/enc Kraus (rectangular).
function _unflatten_kraus(re::Vector{Float64}, im::Vector{Float64},
                          count::Int, row::Int, col::Int)
    ops = Vector{Matrix{ComplexF64}}(undef, count)
    rc = row * col
    for a in 1:count
        M = Matrix{ComplexF64}(undef, row, col)
        base = (a - 1) * rc
        for i in 1:row, j in 1:col
            idx = base + (i - 1) * col + j   # 1-based; = (a-1)*rc + (i-1)*col + (j-1) + 1
            M[i, j] = complex(re[idx], im[idx])
        end
        ops[a] = M
    end
    return ops
end

# ----- C2: aic_assoc_summary_d -----

# _assoc_summary(kraus, n, r; eps, prec) -> NamedTuple
#   (n::Int, dim_A::Int, eps::CertifiedBracket, basis::Vector{Matrix{ComplexF64}}).
# Wraps aic_assoc_summary_d (aic_assoc_shim.h:48): builds A = Img Phi_tilde (the
# eps-C* algebra) from Phi's flat Kraus and reports the ambient n, dim_A, the
# certified associator-defect eps bracket [lo,hi] (FLOOR/CEIL rigorous), and the
# dim_A Frobenius-orthonormal basis {B_k} (n x n, K-major; a SNAPSHOT — its product
# is the Choi-Effros star, NOT plain XY, so do NOT multiply these in Julia).
function _assoc_summary(kraus::Vector{Matrix{ComplexF64}}, n::Int, r::Int;
                        eps::Float64 = -2.0, prec::Int = 256)
    n >= 1 || throw(ArgumentError("_assoc_summary: n must be >= 1, got $n"))
    r >= 1 || throw(ArgumentError("_assoc_summary: need >= 1 Kraus op, got $r"))
    prec >= 2 || throw(ArgumentError("_assoc_summary: prec must be >= 2, got $prec"))
    kr, ki = _flatten_kraus(kraus, n, r)
    out_n = Ref{Cint}(0)
    out_dimA = Ref{Cint}(0)
    eps_assoc = Vector{Cdouble}(undef, 2)
    n4 = n * n * n * n
    basis_re = Vector{Cdouble}(undef, n4)
    basis_im = Vector{Cdouble}(undef, n4)
    rc = GC.@preserve kr ki eps_assoc basis_re basis_im ccall(
        _SYM_ASSOC_SUMMARY_D[], Cint,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
         Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cdouble, Cint),
        out_n, out_dimA, eps_assoc, basis_re, basis_im, kr, ki, n, r, eps, prec)
    rc == 0 || error("aic_assoc_summary_d returned $rc")
    nn = Int(out_n[])
    dim_A = Int(out_dimA[])
    nn == n || error("aic_assoc_summary_d: out_n=$nn != input n=$n")
    (dim_A >= 1 && dim_A <= n * n) ||
        error("aic_assoc_summary_d: dim_A=$dim_A out of range [1, $(n*n)]")
    (isfinite(eps_assoc[1]) && isfinite(eps_assoc[2])) ||
        error("aic_assoc_summary_d: non-finite eps_assoc bracket " *
              "[$(eps_assoc[1]), $(eps_assoc[2])] (lost-precision arb ball?)")
    eps_bracket = CertifiedBracket(eps_assoc[1], eps_assoc[2]; label = "ε_assoc")
    basis = _unflatten_kraus(basis_re, basis_im, dim_A, n, n)
    return (n = nn, dim_A = dim_A, eps = eps_bracket, basis = basis)
end

# ----- C3: aic_main_iso_summary_d -----

# _main_iso_summary(kraus, n, r; eps, prec) -> NamedTuple
#   (n_B::Int, dim_B::Int, blocks::Vector{Int}, cbfwd::CertifiedBracket,
#    cbinv::CertifiedBracket, isodef::CertifiedBracket, eps_used::Float64).
# Wraps aic_main_iso_summary_d (aic_iso_shim.h:55): rebuilds the O(eps)-isomorphism
# v: B -> A (B = (+)_l M_{d_l}, a genuine C* algebra) and reports n_B = Sum d_l,
# dim_B = Sum d_l^2, the block sizes d_l, the certified cb-norm brackets ||v||_cb
# (cbfwd) and ||v^-1||_cb (cbinv) (eig-free, solver-free, OPERATOR-faithful diamond
# bracket per Appendix B2 / FINDINGS §C12), the isomorphism-defect bracket isodef,
# and eps_used (read back so a downstream certifier rebuilds the IDENTICAL v).
function _main_iso_summary(kraus::Vector{Matrix{ComplexF64}}, n::Int, r::Int;
                           eps::Float64 = -2.0, prec::Int = 256)
    n >= 1 || throw(ArgumentError("_main_iso_summary: n must be >= 1, got $n"))
    r >= 1 || throw(ArgumentError("_main_iso_summary: need >= 1 Kraus op, got $r"))
    prec >= 2 || throw(ArgumentError("_main_iso_summary: prec must be >= 2, got $prec"))
    kr, ki = _flatten_kraus(kraus, n, r)
    out_nB = Ref{Cint}(0)
    out_dimB = Ref{Cint}(0)
    out_m = Ref{Cint}(0)
    blocks = Vector{Cint}(undef, n)    # caller pre-allocates length n (m <= n_B <= n)
    cbfwd = Vector{Cdouble}(undef, 2)
    cbinv = Vector{Cdouble}(undef, 2)
    isodef = Vector{Cdouble}(undef, 2)
    out_eps = Vector{Cdouble}(undef, 1)
    rc = GC.@preserve kr ki blocks cbfwd cbinv isodef out_eps ccall(
        _SYM_MAIN_ISO_SUMMARY_D[], Cint,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
         Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint,
         Cdouble, Cint),
        out_nB, out_dimB, out_m, blocks, cbfwd, cbinv, isodef, out_eps,
        kr, ki, n, r, eps, prec)
    rc == 0 || error("aic_main_iso_summary_d returned $rc")
    nB = Int(out_nB[])
    dimB = Int(out_dimB[])
    m = Int(out_m[])
    (m >= 1 && m <= nB && nB <= n) ||
        error("aic_main_iso_summary_d: inconsistent dims m=$m n_B=$nB n=$n")
    isfinite(out_eps[1]) ||
        error("aic_main_iso_summary_d: non-finite eps_used=$(out_eps[1])")
    blk = Int[Int(blocks[l]) for l in 1:m]
    return (n_B = nB, dim_B = dimB, blocks = blk,
            cbfwd = CertifiedBracket(cbfwd[1], cbfwd[2]; label = "‖v‖_cb"),
            cbinv = CertifiedBracket(cbinv[1], cbinv[2]; label = "‖v⁻¹‖_cb"),
            isodef = CertifiedBracket(isodef[1], isodef[2]; label = "iso_def"),
            eps_used = out_eps[1])
end
