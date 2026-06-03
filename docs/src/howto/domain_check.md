# Stay in factorize's domain

After this guide you can determine whether a `UCPMap` is in [`factorize`](@ref)'s
convergence domain and choose the right escape hatch when it is not.

**Assumes:** `AlmostIdempotentChannels` is installed; see [Quick start](@ref).

!!! warning "Conservative threshold"
    `factorize` pre-checks ``\rho(\Phi^2-\Phi) < 0.10`` and throws a clean
    `ArgumentError` when this fails. `associated_algebra` and `main_isomorphism`
    keep the wider ``\rho < 1/4`` basin; `certified_defect` is always safe at any
    ``\eta``. (bug `aic-exa.13`)

## Steps

**1. Know the threshold.**

`factorize` uses a conservative threshold of ``\rho(\Phi^2-\Phi) < 0.10``
(the tighter functional-calculus convergence domain, not the wider `prop_P`
basin ``\rho < 1/4``).

**2. Know which fixtures are in and out.**

| Channel | ``\rho(\Phi^2-\Phi)`` | In factorize domain? |
|---|---|---|
| `phit(2, 0.1)` | ~0.09 | yes |
| `phit(2, 0.3)` | ~0.21 | no |
| `bce_conj(4, 2, 0.02)` | ~0.04 | yes |
| `bce_conj(4, 2, 0.05)` | (small) | yes |
| Any ``\eta = 0`` oracle | 0 | yes |

**3. See the `ArgumentError` for an out-of-domain channel.**

```@example domain_check
using AlmostIdempotentChannels, LinearAlgebra

e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

# phit(2, 0.3): ρ(Φ²−Φ) ≈ 0.21 — outside factorize's threshold of 0.10.
try
    factorize(UCPMap(phit_kraus(2, 0.3)))
catch e
    println(typeof(e))
    println(occursin("convergence domain", e.msg))
end
```

**4. Use `certified_defect` as the escape hatch.**

`certified_defect` has no domain restriction and is always safe.  Even with
``t = 0.3`` (outside factorize's domain):

```jldoctest domain_safe; setup = :(using AlmostIdempotentChannels, LinearAlgebra; e(d,i)=(v=zeros(ComplexF64,d);v[i]=1;v); function phit_kraus(d,t); K=Matrix{ComplexF64}[sqrt(1-t)*Matrix{ComplexF64}(I,d,d)]; for i in 1:d, j in 1:d; push!(K,sqrt(t)*(e(d,i)*e(d,j)')/sqrt(d)); end; return K; end)
julia> b = certified_defect(UCPMap(phit_kraus(2, 0.3)));

julia> b isa CertifiedBracket
true
```

The analytic defect ``0.315`` is contained in the bracket:

```jldoctest domain_contain; setup = :(using AlmostIdempotentChannels, LinearAlgebra; e(d,i)=(v=zeros(ComplexF64,d);v[i]=1;v); function phit_kraus(d,t); K=Matrix{ComplexF64}[sqrt(1-t)*Matrix{ComplexF64}(I,d,d)]; for i in 1:d, j in 1:d; push!(K,sqrt(t)*(e(d,i)*e(d,j)')/sqrt(d)); end; return K; end)
julia> analytic = 0.3 * 0.7 * 2 * (1 - 1/2^2)   # 0.315
0.315

julia> b = certified_defect(UCPMap(phit_kraus(2, 0.3)));

julia> minimum(b) - 1e-7 ≤ analytic ≤ maximum(b) + 1e-7   # the bracket contains it
true
```

**5. Use `associated_algebra` or `main_isomorphism` for the wider domain.**

If you need the ``\varepsilon``-C\*-algebra structure but not the channel
factorization, both verbs accept any channel with ``\rho(\Phi^2-\Phi) < 1/4``
(threshold `_BASIN_RHO_MAX = 0.24`):

```julia
# phit(2, 0.3) is out of factorize's domain but inside the wider prop_P basin.
A = associated_algebra(UCPMap(phit_kraus(2, 0.3)))   # succeeds
v = main_isomorphism(UCPMap(phit_kraus(2, 0.3)))     # succeeds
```

Channels with ``\rho \ge 1/4`` (e.g.\ a reflection ``\Phi(X) = UXU^\dagger``
with ``U = \mathrm{diag}(1, -1)``, where ``\rho = 2``) throw from all three
pipeline verbs; only `certified_defect` is safe there.

On an exactly-idempotent channel ``\rho = 0``; every verb is in domain — see
[The η = 0 oracle](@ref).

## See also

- [Design decisions and known limits](@ref) — the full rationale for the tighter threshold
- [Factorize a channel](@ref) — what to do once you are in domain
- [Certify the idempotency defect (solver-free)](@ref) — the always-safe alternative
