# basin.jl — the prop_P basin pre-check (bead aic-exa.7 [J4], Appendix B8).
#
# WHY THIS EXISTS. associated_algebra / main_isomorphism / factorize call C cores
# that HARD-ASSERT (SIGABRT) when ρ(Φ²−Φ) ≥ 1/4 — the prop_P basin where the
# exact-idempotent regularisation Φ̃ = θ(2Φ−1) is defined (.tex:520-525, FINDINGS
# §C15). A C abort() KILLS the user's Julia session, and Julia CANNOT catch it.
# So those verbs PRE-CHECK ρ IN JULIA and throw a clean ArgumentError first.
#
# THE ESTIMATE. ρ(Φ²−Φ) is the spectral radius of the defect superoperator
# M = Φ²−Φ. We estimate it by power iteration that APPLIES M to a Hermitian iterate
# — Φ(X) = Σ_a K_a† X K_a, M(X) = Φ(Φ(X)) − Φ(X) — and NEVER forms the n²×n²
# superoperator (the vec/Kron trap CLAUDE.md warns of). The seed is deterministic
# (no RNG) so the guard is reproducible; each iterate is symmetrised because M maps
# Hermitian→Hermitian and the prop_P radius lives on the Hermitian subspace.

# Apply the Heisenberg map Φ(X) = Σ_a K_a† X K_a.
_apply_phi(kraus::Vector{Matrix{ComplexF64}}, X::Matrix{ComplexF64}) =
    sum(K' * X * K for K in kraus)

# Power-iteration estimate of ρ(Φ²−Φ). ρ ≈ ‖M(X)‖_F / ‖X‖_F at the fixed point.
function _rho_defect_estimate(Φ::UCPMap; iters::Int=80)
    nn = Φ.n
    K = Φ.kraus
    # deterministic Hermitian seed (fixed pattern, no RNG): X = H + H'.
    H = ComplexF64[(i + 2j) + (i - j) * im for i in 1:nn, j in 1:nn]
    X = H + H'
    X ./= norm(X)
    ρ = 0.0
    for _ in 1:iters
        Y = _apply_phi(K, _apply_phi(K, X)) - _apply_phi(K, X)   # M(X) = Φ²(X) − Φ(X)
        Y = (Y + Y') / 2                                          # stay Hermitian
        ny = norm(Y)
        ny == 0 && return 0.0
        ρ = ny / norm(X)
        X = Y / ny
    end
    return ρ
end

# The prop_P basin is ρ(Φ²−Φ) < 1/4. Refuse with a CONSERVATIVE margin (err toward
# refusing) so a slight under-read of the estimate still keeps the C SIGABRT away.
# certified_defect is the always-safe escape hatch named in the message. Severely
# out-of-basin input may still abort if the estimate badly under-reads — best-effort.
const _BASIN_RHO_MAX = 0.24    # < 1/4

function _assert_in_basin(Φ::UCPMap)
    ρ = _rho_defect_estimate(Φ)
    ρ < _BASIN_RHO_MAX || throw(ArgumentError(
        "Φ is out of the prop_P basin: estimated ρ(Φ²−Φ) ≈ $(round(ρ; sigdigits=4)) ≥ 1/4 " *
        "(approximate_algebras.tex:520-525, FINDINGS §C15). The exact-idempotent " *
        "regularisation Φ̃ = θ(2Φ−1) is undefined there, and the C pipeline would abort. " *
        "Use certified_defect(Φ) (always safe, eig-free) to inspect the idempotency defect instead."))
    return ρ
end

# factorize's safe domain is MUCH TIGHTER than the prop_P basin (bug aic-exa.13).
# It calls aic_factorize_upsilon_build → aic_funcalc_xpow, whose power-series
# functional calculus aborts (SIGABRT — series tail not below tol in 200 terms,
# q at the convergence boundary) well INSIDE the prop_P basin, for in-basin but
# larger-η channels. _assert_in_basin (ρ < 0.24) does NOT keep that abort away.
#
# A SEPARATE-PROCESS sweep of factorize(phit(d,t)) (d∈{2,3,4}, t∈[0.02,0.30],
# 2026-06-02) located the raw C-abort boundary by rho_est = _rho_defect_estimate:
#   d=2: last OK rho_est=0.1275 (t=0.15); SIGABRT from rho_est=0.16 (t=0.20).
#   d=3,4: last OK rho_est=0.1056 (t=0.12); SIGABRT from rho_est=0.1275 (t=0.15).
# So the SMALLEST observed aborting rho_est across the grid is 0.1275. We set the
# factorize guard at 0.10 — below that abort boundary with a ~0.028 margin (and
# below the largest still-OK rho_est 0.1056, i.e. we err toward refusing). This is
# BEST-EFFORT: the boundary is dimension-dependent and the estimate can under-read,
# so a severely out-of-domain channel of an untested SHAPE could still abort.
# certified_defect is the always-safe escape hatch; associated_algebra /
# main_isomorphism keep the wider prop_P guard (_BASIN_RHO_MAX, they do NOT hit
# the xpow abort).
const _FACTORIZE_RHO_MAX = 0.10

function _assert_factorize_domain(Φ::UCPMap)
    ρ = _rho_defect_estimate(Φ)
    ρ < _FACTORIZE_RHO_MAX || throw(ArgumentError(
        "Φ is outside factorize's convergence domain: estimated ρ(Φ²−Φ) ≈ " *
        "$(round(ρ; sigdigits=4)) ≥ $(_FACTORIZE_RHO_MAX). factorize builds Υ via a " *
        "power-series functional calculus (aic_funcalc_xpow inside upsilon_build) " *
        "whose convergence domain is TIGHTER than the prop_P basin ρ < 1/4 (bug " *
        "aic-exa.13): out there the C pipeline ABORTS the session. The channel " *
        "factorization is meaningful only in the small-η (O(η)) regime anyway. " *
        "For larger η, use certified_defect(Φ) (always safe, eig-free) to inspect " *
        "the idempotency defect, or associated_algebra(Φ) / main_isomorphism(Φ), " *
        "which have the wider prop_P domain (ρ < 1/4)."))
    return ρ
end
