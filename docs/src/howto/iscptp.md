# Check trace-preservation of a channel

After this guide you can use [`iscptp`](../reference/api.md) correctly, understand when it
returns `false` by design, and know what the rigorous certificate actually is.

**Assumes:** you have a `CPMap` — typically `encode(F)` or `decode(F)` from a
[`ChannelFactorization`](../reference/api.md) — in scope.

**Warning — iscptp at machine tolerance is not the rigorous certificate.** For ``\eta > 0`` multi-block channels, `iscptp(decode(F))` returns `false`
at the default `atol = 1e-9` by design. The rigorous certificate for the
round-trip quality is the cb-norm defect brackets `delups_defect` /
`upsdel_defect`, not `iscptp` at machine tolerance.

## Steps

**1. Know what `iscptp` checks.**

`iscptp(c; atol = 1e-9)` returns `true` when
``\|\sum_a K_a^\dagger K_a - I\|_{\mathrm{op}} \le \mathrm{atol}``, i.e.\ the
Schrödinger-picture trace-preservation condition.  It is a convenience check,
not a certified bound.

**2. Build a factorization and check each channel.**

```jldoctest iscptp_how; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; F0=factorize(UCPMap(block_cond_exp_kraus(4,2))); F=factorize(UCPMap(bce_conj_kraus(4,2,0.02))))
julia> iscptp(decode(F0))   # η=0 oracle: exact TP
true

julia> iscptp(encode(F0))   # encode always TP
true

julia> iscptp(encode(F))    # encode always TP (even η>0)
true

julia> iscptp(decode(F))    # η>0 multi-block: false at default atol=1e-9
false

julia> iscptp(decode(F); atol = 1e-4)   # true at loosened tolerance
true
```

**3. Understand why `decode` is only ``O(\eta)``-trace-preserving.**

The `decode` channel is only ``O(\eta)``-trace-preserving for ``\eta > 0`` (an
internal Choi→Kraus PSD-cone clip; measured clipped mass ~3.7e-6 on multi-block
oblique fixtures). `iscptp(decode(F))` is `false` at the default `atol = 1e-9`
and `true` at `atol = 1e-4`. The rigorous certificate is the cb-norm round-trip
bracket, not `iscptp` at machine tolerance. The encode channel and the
``\eta = 0`` oracle decode are TP to machine precision.

Concretely: for `bce_conj(4, 2, 0.02)` the measured TP defect
``\|\sum K_a^\dagger K_a - I\|`` of `decode` is ~1.9e-6.  For
`bce_conj(4, 2, 0.05)` it is ~2.2e-5.  Both sit in `(1e-9, 1e-4)`, which is
why the `atol = 1e-4` loosened check passes while the default fails.

**4. Use the round-trip brackets as the rigorous certificate.**

```@example iscptp_brackets
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

F = factorize(UCPMap(bce_conj_kraus(4, 2, 0.02)))
delups_defect(F), upsdel_defect(F)
```

Both brackets certify the round-trip accuracy rigorously; they do not depend on
a TP tolerance.

On an exactly-idempotent channel decode is exactly TP — see
[The η = 0 oracle](../tutorials/eta0_oracle.md).

## See also

- [Multi-block channels](../tutorials/multiblock.md) — why the PSD-cone clip happens and when it matters
- [Design decisions and known limits](../explanation/design_limits.md) — the full O(η)-TP design note
- [API reference](../reference/api.md) — full signature of `iscptp`
