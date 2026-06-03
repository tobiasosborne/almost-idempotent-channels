# Certify the idempotency defect (solver-free)

After this guide you can obtain a rigorous two-sided bracket on
``\eta = \|\Phi^2 - \Phi\|_{\mathrm{cb}}`` without any SDP solver.

**Assumes:** `AlmostIdempotentChannels` is installed and you have a
`UCPMap` in scope; see [Quick start](@ref).

## Steps

**1. Build a `UCPMap` from Kraus operators.**

Inline the constructor body (fixtures are not exported — see step 2 in
[Build a UCPMap from Kraus operators](@ref)):

```@example certify
using AlmostIdempotentChannels, LinearAlgebra

e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

Φ = UCPMap(phit_kraus(2, 0.1))
```

**2. Call `certified_defect`.**

```@example certify
η = certified_defect(Φ)
```

The default `prec = 106` (bits) is sufficient for all standard fixtures.
Pass a higher value for tighter arb balls when the bracket width matters
more than speed:

```julia
η_precise = certified_defect(Φ; prec = 256)
```

`certified_defect` has no domain restriction — it is safe to call at any
``\eta``, including channels that would throw from [`factorize`](@ref).

**3. Read the bracket.**

```@example certify
minimum(η), maximum(η)
```

```@example certify
width(η), midpoint(η)
```

**4. Test containment of a known value.**

The analytic defect for ``\Phi_t`` on ``\mathbb{C}^d`` is
``t(1-t) \cdot 2(1-1/d^2)``. For ``t = 0.1``, ``d = 2``:

```jldoctest certify_contain; setup = :(using AlmostIdempotentChannels, LinearAlgebra; e(d,i)=(v=zeros(ComplexF64,d);v[i]=1;v); function phit_kraus(d,t); K=Matrix{ComplexF64}[sqrt(1-t)*Matrix{ComplexF64}(I,d,d)]; for i in 1:d, j in 1:d; push!(K,sqrt(t)*(e(d,i)*e(d,j)')/sqrt(d)); end; return K; end; Φ=UCPMap(phit_kraus(2,0.1)); η=certified_defect(Φ))
julia> analytic = 0.1 * 0.9 * 2 * (1 - 1/2^2)   # 0.135
0.135

julia> minimum(η) - 1e-7 ≤ analytic ≤ maximum(η) + 1e-7   # the bracket contains it
true
```

The bracket is rigorous, so it contains the true defect; the `1e-7` margin
absorbs the benign double-rounding of the analytic value against the
outward-rounded endpoints.

**5. Note that `value` is `nothing` on the solver-free path.**

```jldoctest certify_value; setup = :(using AlmostIdempotentChannels, LinearAlgebra; e(d,i)=(v=zeros(ComplexF64,d);v[i]=1;v); function phit_kraus(d,t); K=Matrix{ComplexF64}[sqrt(1-t)*Matrix{ComplexF64}(I,d,d)]; for i in 1:d, j in 1:d; push!(K,sqrt(t)*(e(d,i)*e(d,j)')/sqrt(d)); end; return K; end; Φ=UCPMap(phit_kraus(2,0.1)); η=certified_defect(Φ))
julia> value(η)   # nothing — only the MOSEK tight bracket carries a point estimate
```

(`value(η)` returns `nothing`; Documenter suppresses `nothing` output by default.)

!!! note "Loose bracket by design"
    The eig-free bracket is loose by design (hi/lo ~ 2n); it certifies a
    rigorous interval rather than computing the exact value. For ``n = 2``
    the solver-free bracket width is measured at ~2.34e-01.  For a tight
    bracket (width ~1.60e-13) use `certified_defect(Φ; tight = true)` with
    the MOSEK extension. See [Install and use the MOSEK extension](@ref).

On an exactly-idempotent channel ``\eta = 0`` and the bracket collapses to a
neighbourhood of 0 at machine precision — see [The η = 0 oracle](@ref).

## See also

- [Interpret a CertifiedBracket](@ref) — what `lo`/`hi`/`width`/`midpoint`/`value`/`∈` mean
- [Certified arithmetic](@ref) — why the bracket is a proof
- [The five-verb pipeline](@ref) — full walkthrough including algebra extraction and factorization
