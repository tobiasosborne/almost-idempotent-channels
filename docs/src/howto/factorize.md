# Factorize a channel

After this guide you can call [`factorize`](@ref) on an almost-idempotent
channel and use the four accessors to read the encode/decode channels and their
certified round-trip defects.

**Assumes:** your `UCPMap` is in factorize's domain; see
[Stay in factorize's domain](@ref).

!!! warning "factorize has a tighter domain"
    `factorize` builds ``\Upsilon`` via a power-series functional calculus whose
    convergence domain is much smaller than the `prop_P` basin
    ``\rho(\Phi^2-\Phi) < 1/4``. The verb pre-checks a conservative
    ``\rho < 0.10`` in Julia and throws a clean `ArgumentError` — it does not
    abort the session. `associated_algebra` and `main_isomorphism` keep the wider
    ``\rho < 1/4`` domain; `certified_defect` is always safe at any ``\eta``.
    (bug `aic-exa.13`)

## Steps

**1. Pre-check that your channel is in domain.**

Verify ``\rho(\Phi^2-\Phi) < 0.10`` before calling `factorize`.  The cleanest
way is to attempt `certified_defect` first and inspect the bracket midpoint, or
to call `factorize` in a try/catch and redirect to `certified_defect` on
`ArgumentError`.  See [Stay in factorize's domain](@ref) for the full procedure.

**2. Call `factorize`.**

```@example factorize_how
using AlmostIdempotentChannels, LinearAlgebra

function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end
function fixed_unitary(nn)
    H = zeros(ComplexF64, nn, nn)
    for p in 1:(nn - 1); H[p, p + 1] = 0.4; H[p + 1, p] = 0.4; end
    for p in 1:(nn - 2); H[p, p + 2] = 0.3im; H[p + 2, p] = -0.3im; end
    for p in 1:nn; H[p, p] = 0.1 * p; end
    return exp(im * Matrix(Hermitian(H)))
end
function bce_conj_kraus(d, m, t)
    base = block_cond_exp_kraus(d, m); V = fixed_unitary(d); Vd = V'
    K = Matrix{ComplexF64}[]
    for Ka in base; push!(K, sqrt(1 - t) * Ka); end
    for Ka in base; push!(K, sqrt(t) * (Vd * Ka * V)); end
    return K
end

Ψ = UCPMap(bce_conj_kraus(4, 2, 0.02))
F = factorize(Ψ)
```

**3. Read the encode and decode channels.**

`encode(F)` is ``\Upsilon^* : B \to \mathcal{B}(\mathcal{H})``, a CPTP map
from ``B`` into the ambient space.
`decode(F)` is ``\Delta^* : \mathcal{B}(\mathcal{H}) \to B``, a CPTP map
from the ambient space into ``B``.

```@example factorize_how
encode(F)
```

```@example factorize_how
decode(F)
```

**4. Read the two round-trip defect brackets.**

``\|\Delta\Upsilon - \Phi\|_{\mathrm{cb}}`` — how well `decode ∘ encode`
recovers ``\Phi``:

```@example factorize_how
delups_defect(F)
```

``\|\Upsilon\Delta - 1_B\|_{\mathrm{cb}}`` — how well `encode ∘ decode` is the
identity on ``B``:

```@example factorize_how
upsdel_defect(F)
```

Both are ``O(\eta)`` for ``\eta > 0`` channels; the test bound is
`maximum(defect) < 50η`.

**5. Read the block structure of the target algebra.**

```jldoctest fac_blocks; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; Ψ=UCPMap(bce_conj_kraus(4,2,0.02)); F=factorize(Ψ))
julia> blocks(algebra(F))
2-element Vector{Int64}:
 2
 2
```

**6. Handle an `ArgumentError`.**

If `factorize` throws, your channel is outside the tighter convergence domain.
`certified_defect` is always safe and returns a bracket with no domain
restriction:

```julia
try
    F = factorize(Φ)
catch err
    if err isa ArgumentError && occursin("convergence domain", err.msg)
        # Out of factorize's domain. Use the solver-free bracket instead.
        η = certified_defect(Φ)
        @warn "factorize out of domain; solver-free η bracket: $η"
    else
        rethrow(err)
    end
end
```

For the ``C^*``-algebra structure without the channel factorization, call
`associated_algebra(Φ)` or `main_isomorphism(Φ)` — both use the wider
``\rho < 1/4`` basin.

On an exactly-idempotent channel both defect brackets sit below 1e-9; see
[The η = 0 oracle](@ref).

## See also

- [Stay in factorize's domain](@ref) — how to pre-check the domain condition
- [Check trace-preservation of a channel](@ref) — why `iscptp(decode(F))` can be `false`
- [The mathematics](@ref) — the paper result `th_factorization` this implements
- [API reference](@ref) — full signatures for `factorize`, `encode`, `decode`, `delups_defect`, `upsdel_defect`
