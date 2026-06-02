# AICMosekExt.jl — the MOSEK package extension (design doc §3.2/§3.3/§3.4).
#
# Activated automatically when Convex + Mosek + MosekTools are ALL loaded
# alongside AlmostIdempotentChannels. It OWNS the only things the solver-free
# core cannot do (the core-vs-ext capability table, §3.4):
#   * the EXACT cb-norm VALUE  — the Watrous 2012 diamond-norm SDP — exposed via
#     the internal hook AlmostIdempotentChannels._diamond_value_impl (so
#     idempotency_defect(Φ) resolves), and the five diamond_* functions (moved
#     here VERBATIM from the former src/sdp.jl);
#   * the TIGHT certified bracket — certified_defect(Φ; tight=true) — via the
#     internal hook AlmostIdempotentChannels._tight_bracket_impl: it computes the
#     SDP feasible points (MAX-primal X,P,Q; MIN-dual Y0,Y1) and feeds them to the
#     CORE solver-free marshalling wrapper AlmostIdempotentChannels._cbnorm_certify
#     (over aic_cbnorm_certify_d), so the rigorous bracket collapses to the
#     MOSEK-tight ball instead of the eig-free ~2n-wide one.
#
# The core stays loadable + green with NO solver: the stubs (src/sdp_stubs.jl)
# fire until this extension is loaded; this extension only ADDS methods to the
# core's existing dispatch hooks (it does NOT make the core depend on MOSEK).
#
# MOSEK-HANG GUARD (aic-jhe / bug D6/aic-bag): the cb-distance MIN dual can STALL
# for n>=6. Two guards: (1) a MOSEK SOLVER TIME LIMIT (MSK_DPAR_OPTIMIZER_MAX_TIME)
# in _mosek_tight so a stall FAILS LOUD instead of hanging; (2) the tight-bracket
# path BOUNDS n — for n>=6 it falls back to the eig-free bracket (documented note),
# never invoking the dual SDP. The diamond VALUE uses only the MAX primal (no dual).

module AICMosekExt

using AlmostIdempotentChannels
using Convex, Mosek, MosekTools, LinearAlgebra

import AlmostIdempotentChannels: choi_diff, n, kraus, eta_eigfree,
    _diamond_value_impl, _tight_bracket_impl, _cbnorm_certify, CertifiedBracket,
    UCPMap, diamond_norm_watrous, diamond_norm_watrous_primal, diamond_norm_dual

# =============================================================================
# The Watrous 2012 diamond-norm SDP (arXiv:1207.5726), Convex.jl + MOSEK,
# complex-native. This is the eta = ||Phi^2-Phi||_cb = ||Lambda||_diamond engine,
# in three forms (bead aic-m24, increment 3a; design doc docs/cbnorm_tight_certifier.md):
#
#   1. diamond_norm_watrous(J, n)        — MAX primal, value only (the original).
#   2. diamond_norm_watrous_primal(J, n) — MAX primal, EXPOSES the feasible point
#                                          (X, P, Q) for the arb LOWER-bound certifier.
#   3. diamond_norm_dual(J, n)           — QETLAB/Watrous MIN dual, solved as its
#                                          OWN SDP so (Y0, Y1) are PRIMAL variables
#                                          for the arb UPPER-bound certifier.
#
# NORMALIZATION (load-bearing, see HANDOFF.md / ALGORITHM.md module cb-norm):
# a Convention-A Choi has trace n for a UCP map, but the Watrous SDP is trace-1
# calibrated, so ||Lambda||_diamond = (2/n) * MAX_optval. This (2/n) factor is the
# same one the fixture generator and the GADC smoke test use; it is what ties
# eta_idempotence here to the committed golden master eta_ref values.
#
# DUAL normalization DETERMINED EMPIRICALLY (NOT from the design-doc table, which
# was wrong on the partial-trace direction). The MIN dual objective is
#   ||Tr_2(Y0)||_inf + ||Tr_2(Y1)||_inf   (operatornorm = largest eig for HermPSD)
# and  eta = optval / 2.  Pinned against EVERY analytic anchor (id->0; Dep_d-id_d
# at d=2,3,4; phi_t; paper) and the gen_fixtures golden master to ~1e-11. The
# ASYMMETRIC paper example is the discriminating poison-pill: with the partial
# trace over factor 1 (LEFT) its ratio is 4.0 (WRONG); over factor 2 (RIGHT) it
# is 2.0 (CORRECT, == eta_ref). Convention-A J = C[(i-1)n+a,(j-1)n+b] puts the
# K/output index i in the MAJOR (left) position and the H/input index a in the
# MINOR (right) position; partialtrace(.,2,[n,n]) traces the MINOR (input) factor,
# leaving an n x n Hermitian PSD on the output space. DIRECTION = sys 2, CONST = 1/2.
# (The design doc claimed "LEFT (K) factor" + eta=optval/2; the constant is right,
# the direction is sys 2 (RIGHT), as the paper anchor proves.)
#
# Do NOT re-derive the conventions — they are matched to the committed fixtures.
# NO PYTHON (user directive).

# Tightened MOSEK interior-point tolerances. Both solves are small (n^2 x n^2),
# so the cost is negligible and the tighter dual-feasibility cert helps the arb
# step. The primal-dual gap achieved over the anchor set is ~1e-11 (max).
# MSK_DPAR_OPTIMIZER_MAX_TIME is a HARD wall-clock cap so a stalled dual SDP
# (aic-jhe / D6: cb-distance dual can hang for n>=6) FAILS LOUD (MOSEK returns a
# non-OPTIMAL status, the callers assert on it) instead of hanging the session.
const _MOSEK_TOL = (
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS"   => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS"   => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_MU_RED"  => 1.0e-14,
    "MSK_DPAR_OPTIMIZER_MAX_TIME"    => 60.0,
)

# Optimizer factory with the tightened tolerances applied. Returned as a thunk so
# Convex's solve! constructs a fresh optimizer per solve.
function _mosek_tight()
    o = MosekTools.Mosek.Optimizer()
    for (k, v) in _MOSEK_TOL
        MosekTools.MOI.set(o, MosekTools.MOI.RawOptimizerAttribute(k), v)
    end
    return o
end

# ---- 1. Watrous MAX primal, value only (VERIFIED by the fixture generator) ----
function diamond_norm_watrous(J::Matrix{ComplexF64}, n::Int)
    n2 = n * n
    X = ComplexVariable(n2, n2)
    P = HermitianSemidefinite(n2)
    Q = HermitianSemidefinite(n2)
    blk = [P X; X' Q]
    prob = maximize(real(tr(J' * X)),
                    [isposdef(blk),
                     P + Q == Matrix{ComplexF64}(I, n2, n2)])
    solve!(prob, _mosek_tight; silent = true)
    return (2.0 / n) * prob.optval, string(prob.status)
end

# ---- 2. Watrous MAX primal, EXPOSING the feasible point (X, P, Q) ----
# A primal-feasible (X, P, Q) yields the rigorous LOWER bound (arb step 3b)
#   eta >= (2/n) * Re tr(J^dag X).
# Returns (eta, X_val, P_val, Q_val, status). The matrices are evaluate(...)'d to
# dense ComplexF64; eta == (2/n) * prob.optval (== diamond_norm_watrous value).
function diamond_norm_watrous_primal(J::Matrix{ComplexF64}, n::Int)
    n2 = n * n
    X = ComplexVariable(n2, n2)
    P = HermitianSemidefinite(n2)
    Q = HermitianSemidefinite(n2)
    blk = [P X; X' Q]
    prob = maximize(real(tr(J' * X)),
                    [isposdef(blk),
                     P + Q == Matrix{ComplexF64}(I, n2, n2)])
    solve!(prob, _mosek_tight; silent = true)
    Xv = Matrix{ComplexF64}(evaluate(X))
    Pv = Matrix{ComplexF64}(evaluate(P))
    Qv = Matrix{ComplexF64}(evaluate(Q))
    return (2.0 / n) * prob.optval, Xv, Pv, Qv, string(prob.status)
end

# ---- 3. QETLAB/Watrous MIN dual, solved as its OWN SDP ----
# min  ||Tr_2(Y0)||_inf + ||Tr_2(Y1)||_inf   s.t.  [[Y0, -J],[-J^dag, Y1]] >= 0,
#      Y0, Y1 >= 0.   eta = optval / 2.   A dual-feasible (Y0, Y1) yields the
# rigorous UPPER bound (arb step 3b)  eta <= (1/2)(||Tr_2(Y0)|| + ||Tr_2(Y1)||).
# Y0, Y1 are PRIMAL variables of THIS SDP (accessible via evaluate) — we do NOT
# extract MOSEK dual variables (fragile, MOI Bridges bug; design doc risk 1).
# Partial trace over sys 2 (the MINOR/input factor) — empirically pinned above.
# Returns (eta_dual, Y0_val, Y1_val, status).
function diamond_norm_dual(J::Matrix{ComplexF64}, n::Int)
    n2 = n * n
    Y0 = HermitianSemidefinite(n2)
    Y1 = HermitianSemidefinite(n2)
    blk = [Y0 -J; -J' Y1]
    # opnorm(M) = largest singular value = largest eig for a Hermitian PSD marginal.
    obj = opnorm(partialtrace(Y0, 2, [n, n])) +
          opnorm(partialtrace(Y1, 2, [n, n]))
    prob = minimize(obj, [isposdef(blk)])
    solve!(prob, _mosek_tight; silent = true)
    Y0v = Matrix{ComplexF64}(evaluate(Y0))
    Y1v = Matrix{ComplexF64}(evaluate(Y1))
    return prob.optval / 2.0, Y0v, Y1v, string(prob.status)
end

# =============================================================================
# opspace O2 — RECTANGULAR Watrous diamond-norm SDP (bead aic-pjr; design doc
# docs/research/opspace_o2_design.md §2.2). Generalizes the square self-map
# programs above to a map  f: M_{d_in} -> M_{d_out}  whose Convention-A Choi
# J = J(f) is (d_in * d_out) x (d_in * d_out) with the SOURCE index in the
# MAJOR/LEFT kron factor (size d_in) and the TARGET index in the MINOR/RIGHT
# factor (size d_out) — i.e. J(f) = sum_{s,t} E_{st}^{(d_in)} (x) f(E_{st}).
#
# These compute ||f||_diamond = |||f|||_1 (cb-trace norm). For O2 we feed
# J(v*) = J(v)^T to obtain |||v*|||_1 = ||v||_cb (cb-spectral norm, adjoint
# duality, design doc §2.3).
#
# NEITHER the normalization factor NOR the dual partial-trace direction is
# baked in here: both functions return the RAW prob.optval. They were PINNED
# EMPIRICALLY by tools/probe_o2_pin2.jl (+ probe_o2_diag2.jl) against an
# INDEPENDENT closed-form truth (asymmetric CP map ||Psi||_⋄ = sigma_max(A)^2)
# and a complete-isometry oracle (||v||_cb = 1 exactly). PINNED RESULT (corrects
# the design doc's 2/N AND the research leg's 2/n_B):
#   GOLDEN RULE: to get ||f||_⋄ for f: M_in -> M_out, build J = choi_convA(f,in,out)
#   (INPUT-major) and call diamond_*_rect(J, d_maj=in, d_min=out, ...).
#     * normalization factor = 1  (raw optval = ||f||_⋄ EXACTLY).
#     * dual tr_sys = 2  (trace the MINOR/OUTPUT factor, size d_min).
#     * primal rho_on = :major  (density on the INPUT factor; :major != :minor on
#       an asymmetric map, so placement is load-bearing).
#   ||v||_cb   = ||v*||_⋄        : J = choi_convA(v*,    in=N,   out=n_B), dims (N,   n_B).
#   ||v^-1||_cb = ||(v^-1)*||_⋄  : J = choi_convA((v^-1)*, in=n_B, out=N),  dims (n_B, N).
# Do NOT re-derive — match these pinned values.
#
# The functions are parameterized by GENERIC factor dims (d_maj = major/left,
# d_min = minor/right) rather than (d_in, d_out) so the probe can sweep BOTH
# the partial-trace subsystem (tr_sys in {1,2}) and the primal density
# placement (rho_on in {:major,:minor}) without baking in which factor is the
# input. For J(v*) (design doc §1.4) the source/input (A-ambient, size N) is
# MAJOR and the target/output (B-ambient, size n_B) is MINOR, so the caller
# passes d_maj = N, d_min = n_B; the probe then decides tr_sys / rho_on.
# =============================================================================

# ---- Rectangular Watrous MIN dual (UPPER bound), parameterized direction ----
# J is (d_maj*d_min) x (d_maj*d_min). Y0,Y1 are HermitianSemidefinite on that
# space; the block [[Y0,-J],[-J',Y1]] >= 0 is the dual feasibility constraint.
# The objective traces the `tr_sys` factor of the structure dims=[d_maj,d_min]:
#   obj = 1/2 ( ||Tr_{tr_sys}(Y0)||_inf + ||Tr_{tr_sys}(Y1)||_inf ).
# (The leading 1/2 is in the objective, as in diamond_norm_dual above, so the
# probe compares prob.optval directly to the primal.) Returns the RAW
# prob.optval (the probe applies the pinned normalization factor on top).
# Returns (raw_optval, Y0_val, Y1_val, status).
function diamond_dual_rect(J::Matrix{ComplexF64}, d_maj::Int, d_min::Int, tr_sys::Int)
    d = d_maj * d_min
    Y0 = HermitianSemidefinite(d)
    Y1 = HermitianSemidefinite(d)
    blk = [Y0 -J; -J' Y1]
    # opnorm(M) = largest singular value = largest eig for a Hermitian PSD marginal.
    obj = 0.5 * (opnorm(partialtrace(Y0, tr_sys, [d_maj, d_min])) +
                 opnorm(partialtrace(Y1, tr_sys, [d_maj, d_min])))
    prob = minimize(obj, [isposdef(blk)])
    solve!(prob, _mosek_tight; silent = true)
    Y0v = Matrix{ComplexF64}(evaluate(Y0))
    Y1v = Matrix{ComplexF64}(evaluate(Y1))
    return prob.optval, Y0v, Y1v, string(prob.status)
end

# ---- Rectangular Watrous MAX primal (the TRUE value), density form ----
# QETLAB-style Watrous 2012 Thm 6 primal. Two density operators rho0,rho1 on
# the factor named by `rho_on` (:major -> size d_maj, :minor -> size d_min);
# the OTHER factor carries the identity. P,Q are the assembled block diagonals
#   P = kron(rho0, I)  / kron(I, rho0)   (rho_on = :major / :minor)
#   Q = kron(rho1, I)  / kron(I, rho1)
# with the LEFT kron factor MAJOR (matching Convention A's major = leftmost).
# X is the full (d_maj*d_min)-square block variable; [[P,X],[X',Q]] >= 0.
#   maximize Re tr(J' X)   s.t.  Tr rho0 = Tr rho1 = 1.
# The primal value is the TRUE diamond norm (direction-independent — no partial
# trace), so the probe uses it as the strong-duality reference for the asym
# fixture. Returns the RAW prob.optval (probe applies the pinned factor).
# Returns (raw_optval, X_val, P_val, Q_val, status); P,Q evaluated to dense
# ComplexF64 (the assembled kron(rho,I) blocks, for the arb certifier later).
function diamond_primal_rect(J::Matrix{ComplexF64}, d_maj::Int, d_min::Int, rho_on::Symbol)
    d = d_maj * d_min
    rho_dim = (rho_on == :major ? d_maj : d_min)
    rho0 = HermitianSemidefinite(rho_dim)
    rho1 = HermitianSemidefinite(rho_dim)
    Imaj = Matrix{ComplexF64}(I, d_maj, d_maj)
    Imin = Matrix{ComplexF64}(I, d_min, d_min)
    if rho_on == :major
        P = kron(rho0, Imin)
        Q = kron(rho1, Imin)
    elseif rho_on == :minor
        P = kron(Imaj, rho0)
        Q = kron(Imaj, rho1)
    else
        error("diamond_primal_rect: rho_on must be :major or :minor, got $rho_on")
    end
    X = ComplexVariable(d, d)
    blk = [P X; X' Q]
    prob = maximize(real(tr(J' * X)),
                    [isposdef(blk), tr(rho0) == 1.0, tr(rho1) == 1.0])
    solve!(prob, _mosek_tight; silent = true)
    Xv = Matrix{ComplexF64}(evaluate(X))
    Pv = Matrix{ComplexF64}(evaluate(P))
    Qv = Matrix{ComplexF64}(evaluate(Q))
    return prob.optval, Xv, Pv, Qv, string(prob.status)
end

# =============================================================================
# The core dispatch hooks the extension fills (design §3.3/§3.4)
# =============================================================================

# The EXACT cb-norm VALUE eta = ||Phi^2-Phi||_cb = ||Lambda||_diamond. The real
# method for AlmostIdempotentChannels._diamond_value_impl — more specific than the
# core's (::Any,::Int) fallback, so it wins when this extension is loaded. This is
# what idempotency_defect(Φ) calls. Uses only the MAX primal (no dual SDP, so no
# n-bound needed). FAIL LOUD if the solver did not reach OPTIMAL/SLOW_PROGRESS.
function _diamond_value_impl(Φ::UCPMap, prec::Int)
    return _diamond_value_impl(kraus(Φ), n(Φ), prec)
end

# Back-compat raw-Kraus path: idempotency_defect(kraus::Vector{Matrix{ComplexF64}})
# / eta_idempotence(kraus) route here (sdp_stubs.jl:38, the pre-UCPMap surface the
# committed test/runtests.jl T2-T4 still exercise). More specific than the core's
# (::Any,::Int) fallback. The dim n is inferred from the first Kraus op.
function _diamond_value_impl(kraus_ops::Vector{Matrix{ComplexF64}}, prec::Int)
    return _diamond_value_impl(kraus_ops, size(kraus_ops[1], 1), prec)
end

# The shared worker: build J = Choi(Φ²−Φ) and run the MAX-primal Watrous SDP.
# FAIL LOUD if the solver did not reach OPTIMAL/SLOW_PROGRESS (Rule 4).
function _diamond_value_impl(kraus_ops::Vector{Matrix{ComplexF64}}, nn::Int, prec::Int)
    J = choi_diff(kraus_ops, nn; prec = prec)
    eta, status = diamond_norm_watrous(J, nn)
    status in ("OPTIMAL", "SLOW_PROGRESS") || error(
        "idempotency_defect: Watrous diamond-norm SDP did not converge " *
        "(MOSEK status=$status; n=$nn). The cb-norm SDP may have stalled " *
        "(aic-jhe / MSK_DPAR_OPTIMIZER_MAX_TIME). Use certified_defect (solver-free, " *
        "always safe) for a rigorous bracket instead.")
    return eta
end

# The largest ambient dim n for which the tight-bracket path runs the Watrous
# SDPs. The cb-distance MIN dual STALLS for n>=6 (aic-jhe / bug D6/aic-bag), and
# on this hardware even the MAX primal (block (2n^2)x(2n^2)) hits the 60s
# MSK_DPAR_OPTIMIZER_MAX_TIME cap at n=6 (status TIME_LIMIT; measured: n=5 primal
# OPTIMAL in ~41s, n=6 primal TIME_LIMIT at 66s). So n=5 is the largest both SDPs
# reliably clear. ABOVE this we do NOT touch the solver at all (no stall, no
# wasted ~60s) and return the rigorous eig-free bracket with a documented note.
const _TIGHT_DUAL_NMAX = 5

# The MOSEK-tight rigorous bracket on eta = ||Phi^2-Phi||_cb (certified_defect(Φ;
# tight=true), design §3.4). For n<=_TIGHT_DUAL_NMAX: compute the SDP feasible
# points (MAX-primal X,P,Q; MIN-dual Y0,Y1) and feed them + J to the CORE
# solver-free certifier _cbnorm_certify (over aic_cbnorm_certify_d); value = the
# SDP point estimate, and the bracket collapses to the MOSEK-tight ball. For n
# above the bound (the SDP stalls / times out, aic-jhe), DEGRADE GRACEFULLY to the
# rigorous eig-free bracket (still lo<=eta<=hi, just loose) WITHOUT a point
# estimate, and @warn the user — never hang, never return an uncertified value.
# A non-converged SDP at n<=bound is the same class of failure: we degrade to
# eig-free rather than feed a possibly-infeasible point to the certifier (which
# would break rigor), again with a @warn.
function _tight_bracket_impl(Φ::UCPMap, prec::Int)
    nn = n(Φ)
    if nn > _TIGHT_DUAL_NMAX
        @warn "certified_defect(Φ; tight=true): n=$nn exceeds the tight-bracket " *
              "SDP bound (_TIGHT_DUAL_NMAX=$(_TIGHT_DUAL_NMAX)); the Watrous SDP " *
              "stalls/times out for large n (aic-jhe). Returning the rigorous " *
              "eig-free bracket (loose, no point estimate) instead."
        lo, hi = eta_eigfree(kraus(Φ); prec = prec)
        return CertifiedBracket(lo, hi; label = "‖Φ²−Φ‖_cb")
    end
    J = choi_diff(kraus(Φ), nn; prec = prec)
    eta_p, Xv, Pv, Qv, ps = diamond_norm_watrous_primal(J, nn)
    _eta_d, Y0v, Y1v, ds = diamond_norm_dual(J, nn)   # dual is the certifying point
    if !(ps in ("OPTIMAL", "SLOW_PROGRESS")) || !(ds in ("OPTIMAL", "SLOW_PROGRESS"))
        @warn "certified_defect(Φ; tight=true): a Watrous SDP did not converge " *
              "(MAX-primal=$ps, MIN-dual=$ds; n=$nn) — possibly the aic-jhe stall " *
              "(MSK_DPAR_OPTIMIZER_MAX_TIME cap). Returning the rigorous eig-free " *
              "bracket (loose, no point estimate) instead."
        lo, hi = eta_eigfree(kraus(Φ); prec = prec)
        return CertifiedBracket(lo, hi; label = "‖Φ²−Φ‖_cb")
    end
    # The point estimate is the MAX primal (the TRUE value; the dual is the
    # certifying feasible point). They agree to the solver gap (~1e-11).
    return _cbnorm_certify(J, Xv, Pv, Qv, Y0v, Y1v, nn;
                           prec = prec, value = eta_p, label = "‖Φ²−Φ‖_cb")
end

end # module AICMosekExt
