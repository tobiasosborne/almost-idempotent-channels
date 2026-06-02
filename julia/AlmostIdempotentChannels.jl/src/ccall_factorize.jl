# ccall_factorize.jl — low-level @ccall wrappers over the C4/C5 factorize-artifacts
# shims (bead aic-exa.6 [J3]). Split from ccall.jl to keep each file under the
# 200-LOC house limit (Rule 10). INTERNAL (underscore-prefixed, NOT exported);
# bead J4 (api.jl) builds the public `factorize`/`encode`/`decode` verbs on top.
#
# THE TWO-CALL PROTOCOL (the Dec/Enc Kraus counts are DATA-DEPENDENT, design §4.3):
#   1. aic_factorize_artifacts_sizes_d (C4): rebuild the pipeline, report SIZES only
#      (N, n_B, dim_B, num_blocks, the block sizes, AND rDec/rEnc — the Dec/Enc Kraus
#      counts the Choi->Kraus PSD-eigen extraction produces data-dependently).
#   2. allocate dec_*/enc_* from rDec/rEnc, then aic_factorize_artifacts_d (C5) to
#      fill (C5 re-asserts the rebuilt rDec/rEnc match what the caller passed).
#
# THE B1 CHANNEL-DIRECTION BINDING (Appendix B1, the §C13 dual-direction trap —
# LOAD-BEARING). dec_* comes from aic_factorize_dec_kraus (Dec = Delta*, channel
# B(H)->B, n_B x N, layout [a*(nB*N) + i*N + j]); enc_* from aic_factorize_enc_kraus
# (Enc = Upsilon*, channel B->B(H), N x n_B, layout [a*(N*nB) + i*nB + j]). This
# wrapper returns dec from the dec_* args and enc from the enc_* args — NOT crossed.
# A cross flips the dims (dec would be N x n_B): the gate asserts the dims to catch it.
#
# Headers: include/aic/aic_factorize_artifacts_shim.h. Design §4.3 (C4/C5), §9 [J3].
# The eps SENTINEL is identical to C2/C3; out_eta[2] = [eta_proxy, eps_used].

# _factorize_artifacts(kraus, n, r; eps, prec) -> NamedTuple
#   (N::Int, n_B::Int, dim_B::Int, blocks::Vector{Int},
#    dec::Vector{Matrix{ComplexF64}} (each n_B x N),
#    enc::Vector{Matrix{ComplexF64}} (each N x n_B),
#    delups::CertifiedBracket, upsdel::CertifiedBracket,
#    eta_proxy::Float64, eps_used::Float64).
function _factorize_artifacts(kraus::Vector{Matrix{ComplexF64}}, n::Int, r::Int;
                              eps::Float64 = -2.0, prec::Int = 256)
    n >= 1 || throw(ArgumentError("_factorize_artifacts: n must be >= 1, got $n"))
    r >= 1 || throw(ArgumentError("_factorize_artifacts: need >= 1 Kraus op, got $r"))
    prec >= 2 || throw(ArgumentError("_factorize_artifacts: prec must be >= 2, got $prec"))
    kr, ki = _flatten_kraus(kraus, n, r)

    # ----- C4: query the data-dependent sizes -----
    out_N = Ref{Cint}(0)
    out_nB = Ref{Cint}(0)
    out_dimB = Ref{Cint}(0)
    out_m = Ref{Cint}(0)
    blocks = Vector{Cint}(undef, n)   # caller pre-allocates length n (m <= n_B <= n)
    out_rDec = Ref{Cint}(0)
    out_rEnc = Ref{Cint}(0)
    rc = GC.@preserve kr ki blocks ccall(
        _SYM_FACTORIZE_ARTIFACTS_SIZES_D[], Cint,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
         Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cdouble, Cint),
        out_N, out_nB, out_dimB, out_m, blocks, out_rDec, out_rEnc,
        kr, ki, n, r, eps, prec)
    rc == 0 || error("aic_factorize_artifacts_sizes_d returned $rc")
    N = Int(out_N[])
    nB = Int(out_nB[])
    dimB = Int(out_dimB[])
    m = Int(out_m[])
    rDec = Int(out_rDec[])
    rEnc = Int(out_rEnc[])
    (m >= 1 && m <= nB && nB <= n) ||
        error("aic_factorize_artifacts_sizes_d: inconsistent dims m=$m n_B=$nB n=$n")
    (N >= 1 && rDec >= 1 && rEnc >= 1) ||
        error("aic_factorize_artifacts_sizes_d: bad sizes N=$N rDec=$rDec rEnc=$rEnc")
    blk = Int[Int(blocks[l]) for l in 1:m]

    # ----- C5: allocate from the sizes, then fill -----
    # dec: rDec ops each n_B x N, layout [a*(nB*N) + i*N + j].
    # enc: rEnc ops each N x n_B, layout [a*(N*nB) + i*nB + j].
    dec_re = Vector{Cdouble}(undef, rDec * nB * N)
    dec_im = Vector{Cdouble}(undef, rDec * nB * N)
    enc_re = Vector{Cdouble}(undef, rEnc * N * nB)
    enc_im = Vector{Cdouble}(undef, rEnc * N * nB)
    delups = Vector{Cdouble}(undef, 2)
    upsdel = Vector{Cdouble}(undef, 2)
    out_eta = Vector{Cdouble}(undef, 2)
    rc = GC.@preserve kr ki dec_re dec_im enc_re enc_im delups upsdel out_eta ccall(
        _SYM_FACTORIZE_ARTIFACTS_D[], Cint,
        (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint,
         Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
         Cint, Cint, Cdouble, Cint),
        dec_re, dec_im, rDec, enc_re, enc_im, rEnc, delups, upsdel, out_eta,
        kr, ki, n, r, eps, prec)
    rc == 0 || error("aic_factorize_artifacts_d returned $rc")

    (isfinite(out_eta[1]) && isfinite(out_eta[2])) ||
        error("aic_factorize_artifacts_d: non-finite out_eta " *
              "[eta_proxy=$(out_eta[1]), eps_used=$(out_eta[2])]")
    # B1: dec is n_B x N (from dec_*), enc is N x n_B (from enc_*) — NOT crossed.
    dec = _unflatten_kraus(dec_re, dec_im, rDec, nB, N)
    enc = _unflatten_kraus(enc_re, enc_im, rEnc, N, nB)
    return (N = N, n_B = nB, dim_B = dimB, blocks = blk, dec = dec, enc = enc,
            delups = CertifiedBracket(delups[1], delups[2]; label = "‖ΔΥ−Φ‖_cb"),
            upsdel = CertifiedBracket(upsdel[1], upsdel[2]; label = "‖ΥΔ−1_B‖_cb"),
            eta_proxy = out_eta[1], eps_used = out_eta[2])
end
