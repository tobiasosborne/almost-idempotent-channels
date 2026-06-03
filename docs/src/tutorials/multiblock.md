# Multi-block channels

After this page you will have factorized a genuinely oblique, multi-block channel
and read off the block structure ``B = M_2 \oplus M_2``.

## The fixture

`bce_conj_kraus(4, 2, 0.02)` is a convex mix of `block_cond_exp_kraus(4, 2)`
with its unitary conjugate. It is genuinely ``\eta > 0`` and oblique — its image
is not a ``*``-subalgebra of ``M_4`` on the nose; that is the generic case
(see [The mathematics](@ref) for why obliqueness is the interesting regime rather
than the special case). The three constructors it needs are inlined:

```@example multiblock
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
nothing # hide
```

## Block structure

Build the channel and factorize it:

```jldoctest mb_blocks; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; Ψ = UCPMap(bce_conj_kraus(4,2,0.02)); F = factorize(Ψ))
julia> blocks(algebra(F))
2-element Vector{Int64}:
 2
 2
```

The genuine target algebra is ``B = M_2 \oplus M_2``, recovered correctly from
the oblique channel.

## Encode and decode dimensions

For `bce_conj(4, 2, 0.02)` the algebra dimension is ``n_B = 2 + 2 = 4`` and the
ambient dimension is ``N = 4``, so both maps are square:

```jldoctest mb_dims; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; Ψ = UCPMap(bce_conj_kraus(4,2,0.02)); F = factorize(Ψ); en = encode(F); de = decode(F))
julia> (en.dim_in, en.dim_out)
(4, 4)

julia> (de.dim_in, de.dim_out)
(4, 4)
```

## Round-trip defects

The certified round-trip defects are ``O(\eta)``, below the ``50\eta`` test
envelope:

```jldoctest mb_roundtrip; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; Ψ = UCPMap(bce_conj_kraus(4,2,0.02)); F = factorize(Ψ); η = F.eta_proxy)
julia> maximum(delups_defect(F)) < 50 * η
true

julia> maximum(upsdel_defect(F)) < 50 * η
true
```

## Trace-preservation of decode

```jldoctest mb_cptp; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; Ψ = UCPMap(bce_conj_kraus(4,2,0.02)); F = factorize(Ψ))
julia> iscptp(encode(F))
true

julia> iscptp(decode(F))
false

julia> iscptp(decode(F); atol = 1e-4)
true
```

!!! warning "decode is O(η)-trace-preserving for η > 0"
    The `decode` channel is only ``O(\eta)``-trace-preserving for ``\eta > 0``
    (an internal Choi→Kraus PSD-cone clip; measured clipped mass ``\approx
    3.7\times10^{-6}`` on multi-block oblique fixtures). `iscptp(decode(F))` is
    `false` at the default `atol = 1e-9` and `true` at `atol = 1e-4`. The
    rigorous certificate is the cb-norm round-trip bracket, not `iscptp` at
    machine tolerance. The encode channel and the ``\eta = 0`` oracle decode are
    TP to machine precision.

On the ``\eta = 0`` oracle every defect is below `1e-9` and decode is TP to
machine precision — see [The η = 0 oracle](@ref).

## Where to go next

- [Design decisions and known limits](@ref) — root causes for the three known
  limits.
- [The mathematics](@ref) — the four-stage mathematical story with paper anchors.
