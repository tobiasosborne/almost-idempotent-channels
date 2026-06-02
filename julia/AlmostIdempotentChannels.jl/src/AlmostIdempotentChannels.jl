"""
    AlmostIdempotentChannels

Julia bindings for the libaic C/arb cores — the constructive Kitaev pipeline of
Kitaev, *Almost-idempotent quantum channels and approximate C*-algebras*
(arXiv:2405.02434). The headline quantity is the idempotency defect

    eta = ||Phi^2 - Phi||_cb = ||Phi^2 - Phi||_diamond     (approximate_algebras.tex:347-354)

of a UCP self-map Phi: B(H) -> B(H) given by Kraus operators K_a (Heisenberg
picture, Phi(X) = sum_a K_a^dag X K_a).

SOLVER-FREE BY DEFAULT. The core package precompiles, loads, and works with NO
solver installed:
  * `choi_diff`   — J = Choi(Phi^2 - Phi), Convention-A, via aic_ucp_choi_diff_d.
  * `eta_eigfree` — a RIGOROUS eig-free two-sided bracket [lo, hi] on eta, via
                    aic_cbnorm_eigfree_d. No SDP, no eigendecomposition.
The EXACT cb-norm VALUE (`eta_idempotence` / `idempotency_defect`, the Watrous
diamond-norm SDP) needs Convex + Mosek + MosekTools and lives ONLY in the
package extension AICMosekExt; without it the value entry points throw a helpful
install hint (src/sdp_stubs.jl) rather than a MethodError or a missing-solver
crash.

Architecture: the native-library lifecycle (Preferences libpath discovery, the
dlopen'd handle, the dlsym'd symbol table, __init__ with the load-bearing libgmp
preload + RTLD_DEEPBIND) lives in src/libaic.jl. Mirrors the ../su2-fft ccall-C
package pattern. NO PYTHON. Design: docs/research/julia_package_design.md §3, §6.
"""
module AlmostIdempotentChannels

using LinearAlgebra
using Printf

# Native-library lifecycle: Preferences libpath, dlopen handle, dlsym'd symbol
# table, __init__ (libgmp preload + RTLD_DEEPBIND), the MethodError hint.
# Exports nothing itself; defines libaic_path, set_libaic_path!, the _SYM_* Refs.
include("libaic.jl")

# Solver-gated public surface (idempotency_defect + diamond_*) as helpful-error
# stubs; the AICMosekExt extension adds the real methods.
include("sdp_stubs.jl")

# The six immutable value-types (UCPMap, EpsCStarAlgebra, CStarAlgebra,
# MainIsomorphism, ChannelFactorization, CertifiedBracket) + their accessors (bead J2).
# Pure Julia; no ccall (the C-data-driven verbs are bead J4).
include("types.jl")

# Compact + multi-line Base.show for every type (bead J2; design §2.6).
include("show.jl")

# Low-level @ccall wrappers over the C2-C5 ABI shims (bead J3 aic-exa.6). INTERNAL
# (underscore-prefixed, NOT exported); bead J4 (api.jl) builds the public verbs on
# top. ccall.jl = C2/C3 (_assoc_summary, _main_iso_summary) + the _unflatten_kraus
# helper; ccall_factorize.jl = C4/C5 (_factorize_artifacts, the two-call protocol).
# Included AFTER types.jl/show.jl: they construct CertifiedBracket and reuse
# _flatten_kraus (defined in the main module body below — resolved at call time).
include("ccall.jl")
include("ccall_factorize.jl")

export choi_diff, eta_eigfree, eta_idempotence, idempotency_defect,
       diamond_norm_watrous, diamond_norm_watrous_primal, diamond_norm_dual,
       libaic_path, set_libaic_path!

# ----- the value-types + their accessors (bead J2) -----
# NOTE (incremental, per the bead): we export only the types and the accessors
# J2 OWNS. The high-level verbs that need C data (certified_defect,
# associated_algebra, main_isomorphism, factorize, encode, decode, choi) are bead
# J4 and are NOT exported here (no method yet).
export UCPMap, EpsCStarAlgebra, CStarAlgebra, MainIsomorphism, ChannelFactorization,
       CertifiedBracket
# UCPMap accessors:
export n, nkraus, kraus, isunital            # Base.adjoint is extended, not exported
# CertifiedBracket accessors (Base.in/minimum/maximum are Base extensions):
export width, midpoint, value
# CStarAlgebra accessors:
export blocks, dim_B, n_B
# EpsCStarAlgebra accessors:
export dim, ambient, epsilon
# MainIsomorphism accessors:
export cbnorm_forward, cbnorm_inverse, isodefect
# ChannelFactorization accessors (encode/decode are bead J4):
export algebra, delups_defect, upsdel_defect

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

# ----- ccall wrappers (read the dlsym'd symbol Refs from libaic.jl) -----

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
    rc = GC.@preserve kr ki choi_re choi_im ccall(_SYM_CHOI_DIFF_D[], Cint,
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
with it. SOLVER-FREE — the package default for the idempotency defect. See
ALGORITHM.md module cbnorm.
"""
function eta_eigfree(kraus::Vector{Matrix{ComplexF64}}; prec::Int = 106)::Tuple{Float64, Float64}
    n = size(kraus[1], 1)
    r = length(kraus)
    r >= 1 || throw(ArgumentError("need >= 1 Kraus op, got $r"))
    prec >= 2 || throw(ArgumentError("prec must be >= 2, got $prec"))
    kr, ki = _flatten_kraus(kraus, n, r)
    lo = Ref{Cdouble}(0.0)
    hi = Ref{Cdouble}(0.0)
    rc = GC.@preserve kr ki ccall(_SYM_CBNORM_EIGFREE_D[], Cint,
               (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                Cint, Cint, Cint),
               lo, hi, kr, ki, n, r, prec)
    rc == 0 || error("aic_cbnorm_eigfree_d returned $rc")
    return (lo[], hi[])
end

"""
    eta_idempotence(kraus::Vector{Matrix{ComplexF64}}) -> Float64

The EXACT eta-idempotence defect eta = ||Phi^2 - Phi||_cb = ||Phi^2 - Phi||_diamond
(approximate_algebras.tex:347-354) of the UCP self-map Phi given by `kraus`
(Heisenberg picture). The exact VALUE needs the Watrous diamond-norm SDP and so
requires the MOSEK extension (Convex + Mosek + MosekTools); without it this
throws a helpful install hint. The SOLVER-FREE rigorous bracket is `eta_eigfree`.
Back-compat alias of `idempotency_defect`.
"""
eta_idempotence(kraus::Vector{Matrix{ComplexF64}}) = idempotency_defect(kraus)

end # module AlmostIdempotentChannels
