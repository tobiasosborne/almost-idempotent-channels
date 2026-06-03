# Extract the C\*-algebra block structure

After this guide you can read off the block dimensions ``d_l`` of the genuine
C\*-algebra ``B = \bigoplus_l M_{d_l}`` that an almost-idempotent channel
factors through.

**Assumes:** you have a `UCPMap` in domain; see [Stay in factorize's domain](@ref).

## Steps

**1. Choose your entry point.**

Two verbs expose the block list.  Use whichever fits your workflow:

- `factorize(Φ)` → `blocks(algebra(F))` — the full pipeline; also gives the
  encode/decode channels and round-trip defect brackets.
- `main_isomorphism(Φ)` → `blocks(v)` — the isomorphism stage only; skips the
  channel factorization.

Both return the same `Vector{Int}` of block sizes.

**2. Factorize path — multi-block oblique channel.**

```jldoctest bs_fac; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end)
julia> Ψ = UCPMap(bce_conj_kraus(4, 2, 0.02));

julia> F = factorize(Ψ);

julia> blocks(algebra(F))
2-element Vector{Int64}:
 2
 2
```

**3. Isomorphism path — same data, different object.**

`main_isomorphism` uses the wider ``\rho < 1/4`` basin (not the tighter
factorize threshold) and gives the blocks without the encode/decode channels:

```jldoctest bs_iso; setup = :(using AlmostIdempotentChannels, LinearAlgebra; e(d,i)=(v=zeros(ComplexF64,d);v[i]=1;v); function phit_kraus(d,t); K=Matrix{ComplexF64}[sqrt(1-t)*Matrix{ComplexF64}(I,d,d)]; for i in 1:d, j in 1:d; push!(K,sqrt(t)*(e(d,i)*e(d,j)')/sqrt(d)); end; return K; end)
julia> Φ = UCPMap(phit_kraus(2, 0.1));

julia> v = main_isomorphism(Φ);

julia> blocks(v)
1-element Vector{Int64}:
 2
```

Note: `blocks(v)` and `blocks(algebra(F))` are the same data via different
objects.  Do not let setup variable names collide if you build both in one
session.

**4. Interpret the integers.**

Each ``d_l`` in the list is the side length of a full matrix block.
``B = \bigoplus_l M_{d_l}`` has vector-space dimension ``\sum_l d_l^2``
(``\mathtt{dim\_B}``) and block-diagonal representative size ``\sum_l d_l``
(``\mathtt{n\_B}``).  For the ``[2, 2]`` example: ``\dim_B = 8``,
``n_B = 4``.

**5. Use the ``\eta = 0`` oracle for ground-truth verification.**

When ``\Phi`` is exactly idempotent the block structure is exact and every
certified defect is below 1e-9 (measured 4.4e-75 and 3.9e-75 at prec=128):

```jldoctest bs_oracle; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end)
julia> F0 = factorize(UCPMap(block_cond_exp_kraus(4, 2)));

julia> blocks(algebra(F0))
2-element Vector{Int64}:
 2
 2

julia> maximum(delups_defect(F0)) < 1e-9
true

julia> maximum(upsdel_defect(F0)) < 1e-9
true
```

See [The η = 0 oracle](@ref) for more on this calibration device.

## See also

- [The mathematics](@ref) — why ``B = \bigoplus_l M_{d_l}`` is a genuine C\*-algebra
- [Multi-block channels](@ref) — when and why multiple blocks appear
- [The five-verb pipeline](@ref) — end-to-end worked example including block extraction
