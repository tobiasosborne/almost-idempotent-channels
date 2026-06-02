# test_aqua.jl — the Aqua.jl quality gate (bead aic-exa.9 [T], design §9 step 3).
# Aqua catches the structural defects unit tests miss: method ambiguities, undefined
# exports, stale / under-constrained deps, and TYPE PIRACY. The load-bearing piracy
# check here is the `factorize` LinearAlgebra extension (Appendix B6): we EXTEND
# `LinearAlgebra.factorize` with a method whose first argument is our own `UCPMap`,
# so it is NOT piracy (the method is "owned" by our type). Aqua confirms this.
#
# undocumented_names is OFF: the docstrings are finalized in bead aic-exa.10 [X]
# (Documenter + doctests), so gating on them now would be premature. TODO(aic-exa.10
# [X]): re-enable `undocumented_names=true` once the docs land and every export has
# a finalized docstring.
#
# Assumes `using AlmostIdempotentChannels, Test` + `using Aqua`.

@testset "Aqua quality gate" begin
    Aqua.test_all(AlmostIdempotentChannels;
                  ambiguities=true,
                  unbound_args=true,
                  undefined_exports=true,
                  stale_deps=true,
                  deps_compat=true,
                  piracies=true,
                  project_extras=true,
                  undocumented_names=false)   # TODO(aic-exa.10 [X]): enable when docs land
end
