# sdp.jl — the Watrous 2012 diamond-norm SDP (arXiv:1207.5726), copied VERBATIM
# from tools/gen_fixtures_d24.jl (bead aic-d24). Convex.jl + MOSEK, complex-native.
#
# NORMALIZATION (load-bearing, see HANDOFF.md / ALGORITHM.md module cb-norm):
# a Convention-A Choi has trace n for a UCP map, but the Watrous SDP is trace-1
# calibrated, so ||Lambda||_diamond = (2/n) * SDP_optval. This (2/n) factor is the
# same one the fixture generator and the GADC smoke test use; it is what ties
# eta_idempotence here to the committed golden master eta_ref values.
#
# Do NOT re-derive — this is a verbatim copy so the SDP value matches the fixtures.
# NO PYTHON (user directive).

# ---- the Watrous diamond-norm SDP (VERIFIED by the fixture generator anchors) ----
function diamond_norm_watrous(J::Matrix{ComplexF64}, n::Int)
    n2 = n * n
    X = ComplexVariable(n2, n2)
    P = HermitianSemidefinite(n2)
    Q = HermitianSemidefinite(n2)
    blk = [P X; X' Q]
    prob = maximize(real(tr(J' * X)),
                    [isposdef(blk),
                     P + Q == Matrix{ComplexF64}(I, n2, n2)])
    solve!(prob, MosekTools.Mosek.Optimizer; silent = true)
    return (2.0 / n) * prob.optval, string(prob.status)
end
