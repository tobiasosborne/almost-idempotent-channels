"""
    AlmostIdempotentChannels

Julia bindings for the libaic C/arb cores — the reusable eta-idempotence value
entry point of Kitaev, *Almost-idempotent quantum channels and approximate
C*-algebras* (arXiv:2405.02434). The headline quantity is the idempotency defect

    eta = ||Phi^2 - Phi||_cb = ||Phi^2 - Phi||_diamond     (approximate_algebras.tex:347-354)

of a UCP self-map Phi: B(H) -> B(H) given by Kraus operators K_a (Heisenberg
picture, Phi(X) = sum_a K_a^dag X K_a). Downstream modules (assoc_ecsa,
factorize) consume `eta_idempotence(kraus) -> Float64`.

Architecture (Route B, bead aic-d24/aic-m24): the certified Choi matrix
J = Choi(Phi^2 - Phi) is built in the C/arb core via the flat-double ccall shim
(`aic_ucp_choi_diff_d` in build/libaic.so), then the Watrous 2012 diamond-norm
SDP is solved in Julia + MOSEK (src/sdp.jl). A second shim
(`aic_cbnorm_eigfree_d`) returns a RIGOROUS eig-free bracket [lo, hi] on eta with
no solver, used to cross-check the SDP value.

Mirrors the ../su2-fft/julia ccall-C package pattern (dlopen + dlsym). NO PYTHON.
Bead: aic-m24, increment 2b.
"""
module AlmostIdempotentChannels

using Libdl
using Convex, Mosek, MosekTools, LinearAlgebra

# Load LIBAIC (absolute path to build/libaic.so, set by deps/build.jl).
const DEPS_FILE = joinpath(@__DIR__, "..", "deps", "deps.jl")
if !isfile(DEPS_FILE)
    error("AlmostIdempotentChannels: deps.jl not found at $DEPS_FILE. " *
          "Run `Pkg.build(\"AlmostIdempotentChannels\")` first.")
end
include(DEPS_FILE)

# The Watrous diamond-norm SDP (verbatim from tools/gen_fixtures_d24.jl).
include("sdp.jl")

# Module-private handle to the dlopen'd libaic.so and dlsym'd shim symbols.
# RTLD flags copied from ../su2-fft/julia/src/SU2FFT.jl: RTLD_DEEPBIND so
# libaic's transitive deps (libflint, libgmp) resolve against the system search
# order, not Julia's pre-loaded symbol table.
const _LIBAIC_HANDLE       = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_CHOI_DIFF_D     = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_CBNORM_EIGFREE_D = Ref{Ptr{Cvoid}}(C_NULL)

function __init__()
    # Pre-load the SYSTEM libgmp.so.10 with RTLD_GLOBAL before libaic/libflint
    # are mapped. This works around the Julia bundled-libgmp ABI mismatch: on
    # this system RTLD_DEEPBIND on libaic alone was insufficient because libflint
    # resolves __gmpn_* symbols against the already-loaded Julia-bundled libgmp;
    # pre-loading the system libgmp keeps its definitions in the global symbol
    # table first (the exact issue ../su2-fft handles, bead su2fft-e5z).
    try
        Libdl.dlopen("/lib/x86_64-linux-gnu/libgmp.so.10",
                     Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL)
    catch
        try
            Libdl.dlopen("libgmp.so.10", Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL)
        catch
            # silently ignore on systems without a discoverable system libgmp
        end
    end
    flags = Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL | Libdl.RTLD_DEEPBIND
    _LIBAIC_HANDLE[]        = Libdl.dlopen(LIBAIC, flags)
    _SYM_CHOI_DIFF_D[]      = Libdl.dlsym(_LIBAIC_HANDLE[], :aic_ucp_choi_diff_d)
    _SYM_CBNORM_EIGFREE_D[] = Libdl.dlsym(_LIBAIC_HANDLE[], :aic_cbnorm_eigfree_d)
end

export choi_diff, eta_eigfree, eta_idempotence, diamond_norm_watrous,
       diamond_norm_watrous_primal, diamond_norm_dual, libaic_path

# ----- marshalling helpers (Julia ComplexF64 matrices <-> flat C arrays) -----

# Flatten r Kraus matrices (each n x n, K_a: H->K) into the C row-major layout
# kraus[a*n*n + i*n + j] = K_a[i+1, j+1] (aic_ucp_shim.h:18-19; 0-based C, 1-based
# Julia). Returns (kraus_re, kraus_im) as flat Vector{Float64} of length r*n*n.
function _flatten_kraus(kraus::Vector{Matrix{ComplexF64}}, n::Int, r::Int)
    kr = Vector{Float64}(undef, r * n * n)
    ki = Vector{Float64}(undef, r * n * n)
    for a in 1:r
        K = kraus[a]
        size(K) == (n, n) || throw(ArgumentError(
            "Kraus op $a has size $(size(K)), expected ($n, $n)"))
        base = (a - 1) * n * n
        for i in 1:n, j in 1:n
            idx = base + (i - 1) * n + j   # 1-based; = (a-1)nn + (i-1)n + (j-1) + 1
            kr[idx] = real(K[i, j])
            ki[idx] = imag(K[i, j])
        end
    end
    return kr, ki
end

# ----- ccall wrappers -----

"""
    choi_diff(kraus::Vector{Matrix{ComplexF64}}, n::Int; prec::Int=106) -> Matrix{ComplexF64}

Build J = Choi(Phi^2 - Phi) (Convention-A, n^2 x n^2) for the UCP self-map Phi
given by `kraus`, via the certified arb core (`aic_ucp_choi_diff_d`). `prec` is
the arb working precision (default 106 for headroom; 53 cross-checks the double
path). Convention-A: C[i*n+a, j*n+b] = sum_x conj(Kx[i,a]) Kx[j,b] — the
conjugation is on the FIRST (i,a) factor (aic_ucp_shim.h:20-22).
"""
function choi_diff(kraus::Vector{Matrix{ComplexF64}}, n::Int; prec::Int = 106)::Matrix{ComplexF64}
    n >= 1 || throw(ArgumentError("n must be >= 1, got $n"))
    r = length(kraus)
    r >= 1 || throw(ArgumentError("need >= 1 Kraus op, got $r"))
    prec >= 2 || throw(ArgumentError("prec must be >= 2, got $prec"))
    kr, ki = _flatten_kraus(kraus, n, r)
    N = n * n
    choi_re = Vector{Float64}(undef, N * N)
    choi_im = Vector{Float64}(undef, N * N)
    rc = ccall(_SYM_CHOI_DIFF_D[], Cint,
               (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                Cint, Cint, Cint),
               choi_re, choi_im, kr, ki, n, r, prec)
    rc == 0 || error("aic_ucp_choi_diff_d returned $rc")
    # Unflatten C row-major choi[p*N + q] (0-based) into a Julia (N x N) matrix
    # J[p+1, q+1]. C index p*N + q (0-based) -> flat[p*N + q + 1] (1-based).
    J = Matrix{ComplexF64}(undef, N, N)
    for p in 1:N, q in 1:N
        idx = (p - 1) * N + q
        J[p, q] = complex(choi_re[idx], choi_im[idx])
    end
    return J
end

"""
    eta_eigfree(kraus::Vector{Matrix{ComplexF64}}; prec::Int=106) -> (lo::Float64, hi::Float64)

Certified eig-free two-sided bracket on eta = ||Phi^2-Phi||_cb, via
`aic_cbnorm_eigfree_d`. Returns RIGOROUS bounds with `lo <= eta <= hi` (the C
shim rounds the arb balls outward with FLOOR/CEIL to keep the bracket rigorous
after the double conversion). No SDP, no eigendecomposition. The bracket is loose
by design (ratio hi/lo ~ 2n); it certifies the MOSEK value rather than competing
with it. See ALGORITHM.md module cbnorm.
"""
function eta_eigfree(kraus::Vector{Matrix{ComplexF64}}; prec::Int = 106)::Tuple{Float64, Float64}
    n = size(kraus[1], 1)
    r = length(kraus)
    r >= 1 || throw(ArgumentError("need >= 1 Kraus op, got $r"))
    prec >= 2 || throw(ArgumentError("prec must be >= 2, got $prec"))
    kr, ki = _flatten_kraus(kraus, n, r)
    lo = Ref{Cdouble}(0.0)
    hi = Ref{Cdouble}(0.0)
    rc = ccall(_SYM_CBNORM_EIGFREE_D[], Cint,
               (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                Cint, Cint, Cint),
               lo, hi, kr, ki, n, r, prec)
    rc == 0 || error("aic_cbnorm_eigfree_d returned $rc")
    return (lo[], hi[])
end

"""
    eta_idempotence(kraus::Vector{Matrix{ComplexF64}}) -> Float64

The eta-idempotence defect eta = ||Phi^2 - Phi||_cb = ||Phi^2 - Phi||_diamond
(approximate_algebras.tex:347-354) of the UCP self-map Phi given by `kraus`
(Heisenberg picture, Phi(X) = sum_a K_a^dag X K_a). Builds J = Choi(Phi^2-Phi)
in the certified arb core (`choi_diff`), then solves the Watrous 2012 diamond-norm
SDP in Julia + MOSEK (`diamond_norm_watrous`). This is the reusable value entry
point downstream modules (assoc_ecsa, factorize) call.

Fails loud (assert) if the SDP status is not OPTIMAL.
"""
function eta_idempotence(kraus::Vector{Matrix{ComplexF64}})::Float64
    n = size(kraus[1], 1)
    J = choi_diff(kraus, n)
    val, status = diamond_norm_watrous(J, n)
    status == "OPTIMAL" || error("Watrous SDP not OPTIMAL (status=$status)")
    return val
end

# ----- diagnostics -----

"""
    libaic_path() -> String

Absolute path to the dynamically-loaded build/libaic.so, as recorded by
deps/build.jl. For debugging linkage.
"""
libaic_path() = LIBAIC

end # module AlmostIdempotentChannels
