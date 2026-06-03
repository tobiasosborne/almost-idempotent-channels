# Build a UCPMap from Kraus operators

After this guide you can construct a [`UCPMap`](../reference/api.md) from a list of Kraus
operators and read its basic accessors.

**Assumes:** `AlmostIdempotentChannels` and `LinearAlgebra` are loaded.

## Steps

**1. Know the Heisenberg-picture convention.**

A `UCPMap` represents a unital completely positive (UCP) self-map
``\Phi : \mathcal{B}(\mathcal{H}) \to \mathcal{B}(\mathcal{H})`` in the
Heisenberg picture:

```math
\Phi(X) = \sum_{a=1}^r K_a^\dagger X K_a
```

Each Kraus operator ``K_a`` is an ``n \times n`` complex matrix.  Unitality
means ``\sum_a K_a^\dagger K_a = I``.

**2. Assemble a `Vector{Matrix{ComplexF64}}`.**

The constructor accepts any `AbstractVector` of `AbstractMatrix` and promotes to
`ComplexF64`:

```@example build_ucpmap
using AlmostIdempotentChannels, LinearAlgebra

# Block conditional expectation on C^4: first 2 vs last 2.
# K0 = diag(1,1,0,0),  K1 = diag(0,0,1,1).
function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end

kraus = block_cond_exp_kraus(4, 2)
Φ = UCPMap(kraus)
```

**3. Read the structural accessors.**

```jldoctest build_ucpmap_acc; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; Φ=UCPMap(block_cond_exp_kraus(4,2)))
julia> n(Φ)       # ambient Hilbert space dimension
4

julia> nkraus(Φ)  # number of Kraus operators
2

julia> isunital(Φ)   # ‖Σ K_a†K_a − I‖ ≤ 1e-9
true
```

`kraus(Φ)` returns the stored `Vector{Matrix{ComplexF64}}` (the same objects
passed in, not copies).

**4. Read the constructor's unitality check.**

By default (`check = true`) the constructor verifies ``\|\sum_a K_a^\dagger K_a
- I\| \le \mathrm{atol}`` (default `atol = 1e-9`) and throws `ArgumentError`
naming the norm if it fails.  To suppress the check (e.g.\ for a dual map
known to be trace-preserving rather than unital) pass `check = false`:

```julia
Φ_unchecked = UCPMap(kraus; check = false)
```

**5. Use a second constructor: from the Choi matrix.**

If you have a Convention-A Choi matrix ``C`` (size ``n^2 \times n^2``) you can
round-trip through the Choi path:

```julia
C = choi(Φ)            # Convention-A Choi of Φ
Φ2 = UCPMap(; choi = C, n = n(Φ))
```

This reconstructs Kraus operators via PSD eigendecomposition and then validates
unitality identically to the Kraus path.

**6. Work through a second example: a convex mix with a conjugated channel.**

The oblique fixtures in the test suite follow this pattern (inlined here — they
are NOT exported):

```@example build_ucpmap
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

Ψ = UCPMap(bce_conj_kraus(4, 2, 0.02))   # η>0 oblique multi-block
```

```jldoctest build_ucpmap_psi; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; Ψ=UCPMap(bce_conj_kraus(4,2,0.02)))
julia> n(Ψ), nkraus(Ψ)
(4, 4)

julia> isunital(Ψ)
true
```

On an exactly-idempotent channel every downstream certified defect is below
1e-9 — see [The η = 0 oracle](../tutorials/eta0_oracle.md).

## See also

- [API reference](../reference/api.md) — full `UCPMap` signature including the Choi constructor and `Base.adjoint`
- [The mathematics](../explanation/math_story.md) — the Heisenberg vs Schrödinger picture and the UCP condition
- [The five-verb pipeline](../tutorials/pipeline.md) — end-to-end walkthrough using `UCPMap` as the entry point
