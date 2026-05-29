# sdp.jl — the Watrous 2012 diamond-norm SDP (arXiv:1207.5726), Convex.jl + MOSEK,
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
const _MOSEK_TOL = (
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS"   => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS"   => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1.0e-12,
    "MSK_DPAR_INTPNT_CO_TOL_MU_RED"  => 1.0e-14,
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
