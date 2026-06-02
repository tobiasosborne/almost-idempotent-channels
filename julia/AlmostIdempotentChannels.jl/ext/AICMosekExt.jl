# AICMosekExt.jl — the MOSEK package extension (design doc §3.2/§3.3).
#
# Activated automatically when Convex + Mosek + MosekTools are ALL loaded
# alongside AlmostIdempotentChannels. It OWNS the exact cb-norm value: the
# Watrous diamond-norm SDP (currently src/sdp.jl), the real method for the
# internal dispatch hook AlmostIdempotentChannels._diamond_value_impl (so
# idempotency_defect resolves), and the tight-bracket path.
#
# SKELETON (bead J1, aic-exa.4): this is a valid, loadable extension that imports
# the package and the three weakdeps so the weakdeps-as-test-deps wiring is
# verified end to end. The real methods — the five diamond_* functions and the
# _diamond_value_impl SDP — land VERBATIM (from src/sdp.jl) in bead aic-exa.8
# [J5]. Until then loading this extension simply confirms the extension graph
# resolves; the core stubs (src/sdp_stubs.jl) still fire because no real method
# is defined yet.

module AICMosekExt

using AlmostIdempotentChannels
using Convex, Mosek, MosekTools, LinearAlgebra

# Real methods land in bead aic-exa.8 [J5]:
#   AlmostIdempotentChannels._diamond_value_impl(kraus, prec)  — the Watrous SDP
#   AlmostIdempotentChannels.diamond_norm_watrous / _primal / _dual
#   the tight-bracket path (aic_cbnorm_certify_d fed MOSEK feasible points)

end # module AICMosekExt
