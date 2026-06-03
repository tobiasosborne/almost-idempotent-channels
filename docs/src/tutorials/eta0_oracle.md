# The η = 0 oracle

After this page you will have watched every certified defect collapse to machine
``\varepsilon`` on exactly-idempotent channels — the cleanest ground truth in the
paper.

This is the calibration device for the package. Every other page that produces a
certified quantity should reduce to 0 here. If it does not, something is wrong.

## The three standard oracles

When ``\Phi`` is exactly a conditional expectation (``\eta = 0``), the whole
construction collapses to the genuine Choi–Effros structure: every certified defect
drops to machine ``\varepsilon`` and the extracted block structure is exact.

The three constructors below live in `test/fixtures.jl` and are inlined here:

```jldoctest oracle_blocks
julia> using AlmostIdempotentChannels, LinearAlgebra

julia> identity_kraus(d) = [Matrix{ComplexF64}(I, d, d)];

julia> function block_cond_exp_kraus(d, m)
           K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
           for i in 1:m; K0[i, i] = 1; end
           for i in (m + 1):d; K1[i, i] = 1; end
           return Matrix{ComplexF64}[K0, K1]
       end;

julia> function dephasing_kraus(d)
           K = Matrix{ComplexF64}[]
           for i in 1:d
               Ki = zeros(ComplexF64, d, d); Ki[i, i] = 1; push!(K, Ki)
           end
           return K
       end;

julia> F1 = factorize(UCPMap(identity_kraus(2)));

julia> blocks(algebra(F1))
1-element Vector{Int64}:
 2

julia> F2 = factorize(UCPMap(block_cond_exp_kraus(4, 2)));

julia> blocks(algebra(F2))
2-element Vector{Int64}:
 2
 2

julia> F3 = factorize(UCPMap(dephasing_kraus(3)));

julia> blocks(algebra(F3))
3-element Vector{Int64}:
 1
 1
 1
```

The block structures are exact:
- `identity(2)` → ``B = M_2``
- `block_cond_exp(4, 2)` → ``B = M_2 \oplus M_2``
- `dephasing(3)` → ``B = M_1 \oplus M_1 \oplus M_1`` (commutative)

## Defect collapse

On the ``\eta = 0`` oracle every certified defect is well below `1e-9`. For
`block_cond_exp(4, 2)`:

```jldoctest oracle_defects; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d, m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; F2 = factorize(UCPMap(block_cond_exp_kraus(4, 2))))
julia> maximum(delups_defect(F2)) < 1e-9
true

julia> maximum(upsdel_defect(F2)) < 1e-9
true

julia> iscptp(decode(F2))
true

julia> iscptp(encode(F2))
true
```

At `prec=128` the measured maxima are 4.4e-75 (`delups_defect`) and 3.9e-75
(`upsdel_defect`) — well below any floating-point floor (`README.md`
cross-check ladder). These numbers are prose, not jldoctest literals: the test
assertions only bound them at `< 1e-9`.

On the ``\eta = 0`` oracle the `decode` channel is TP to machine precision
(no PSD-cone clip), so `iscptp(decode(F2))` is `true` at the default
`atol = 1e-9`. Contrast with the ``\eta > 0`` multi-block case, where decode is
only ``O(\eta)``-TP — see [Multi-block channels](multiblock.md) and the Known limits in
[Design decisions and known limits](../explanation/design_limits.md).

![On exactly-idempotent channels every certified defect collapses below 1e-9; measured maxima are 4.4e-75 and 3.9e-75 (test_factorize.jl)](../assets/eta0_oracle.png)

## The oracle as a calibration motif

Every page in this documentation that produces a certified quantity notes that the
same quantity is 0 on an ``\eta = 0`` oracle. This page is the proof-by-picture
behind that note. If you add a new channel fixture and any defect does not collapse
below `1e-9` when ``\eta = 0``, treat it as a regression.

For the cross-check ladder — where the ``\eta = 0`` oracle is rung 3 — see
[Certified arithmetic](../explanation/certified_arithmetic.md).

## Where to go next

- [The five-verb pipeline](pipeline.md) — the full chain on a genuinely ``\eta > 0``
  channel.
- [Certified arithmetic](../explanation/certified_arithmetic.md) — the four-rung cross-check ladder where this
  oracle is rung 3.
