# The five-verb pipeline

After this page you will have run the whole chain — certify the defect, extract
the ``\varepsilon``-C\*-algebra, factor the channel through a genuine C\*-algebra
— on a real near-idempotent channel.

Every code block on this page is executed at documentation-build time. Float-heavy
outputs use `@example` blocks (embedded live, robust to last-digit wobble);
structural outputs — block sizes, dims, booleans — use `jldoctest` blocks whose
output is checked exactly.

## Fixtures

The constructors below live in the test suite (`test/fixtures.jl`), not in the
package API. They are inlined here so the page is self-contained:

```@example tut
using AlmostIdempotentChannels, LinearAlgebra

e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)

# Φ_t = (1−t)·id + t·(maximally depolarizing) on C^d.
# Closed-form defect: ‖Φ_t²−Φ_t‖_cb = t(1−t)·2(1−1/d²).
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

# Block conditional expectation on C^d: first m coords vs the rest. η = 0.
function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end

# A fixed basis-mixing unitary V = exp(iH).
function fixed_unitary(nn)
    H = zeros(ComplexF64, nn, nn)
    for p in 1:(nn - 1); H[p, p + 1] = 0.4; H[p + 1, p] = 0.4; end
    for p in 1:(nn - 2); H[p, p + 2] = 0.3im; H[p + 2, p] = -0.3im; end
    for p in 1:nn; H[p, p] = 0.1 * p; end
    return exp(im * Matrix(Hermitian(H)))
end

# bce_conj(d,m,t): block_cond_exp(d,m) mixed with its unitary conjugate.
# Genuinely η>0, oblique, multi-block; in factorize's domain for small t.
function bce_conj_kraus(d, m, t)
    base = block_cond_exp_kraus(d, m); V = fixed_unitary(d); Vd = V'
    K = Matrix{ComplexF64}[]
    for Ka in base; push!(K, sqrt(1 - t) * Ka); end
    for Ka in base; push!(K, sqrt(t) * (Vd * Ka * V)); end
    return K
end
nothing # hide
```

## Step 1 — build the channel and certify the defect

[`UCPMap`](@ref) wraps the Kraus operators and checks Heisenberg unitality
``\sum_a K_a^\dagger K_a = I``:

```@example tut
Φ = UCPMap(phit_kraus(2, 0.1))
```

The channel size and Kraus count are stable structural facts:

```jldoctest tut1
julia> using AlmostIdempotentChannels, LinearAlgebra

julia> e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v);

julia> function phit_kraus(d, t)
           K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
           for i in 1:d, j in 1:d
               push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
           end
           return K
       end;

julia> Φ = UCPMap(phit_kraus(2, 0.1));

julia> n(Φ), nkraus(Φ)
(2, 5)
```

[`certified_defect`](@ref) returns a rigorous ``\|\Phi^2-\Phi\|_{\mathrm{cb}}``
bracket — no SDP solver required:

```@example tut
η = certified_defect(Φ)
```

Test containment and read the bracket geometry:

```@example tut
(0.0 in η, width(η), midpoint(η))
```

!!! note
    The eig-free bracket is loose by design (`hi/lo ~ 2n`). It certifies a value
    rather than computing it. For a tight bracket or the exact value, use the MOSEK
    extension — see [Tight brackets with MOSEK](@ref).

On an exactly-idempotent channel (``\eta = 0``), the bracket collapses to machine
``\varepsilon``. See [The η = 0 oracle](@ref) for a worked demonstration.

## Step 2 — the associated ``\varepsilon``-C\*-algebra

[`associated_algebra`](@ref) builds the ``\varepsilon``-C\*-algebra ``A``; for
what that means, see [The mathematics](@ref):

```@example tut
A = associated_algebra(Φ)
```

```@example tut
dim(A), ambient(A), epsilon(A)
```

## Step 3 — the ``O(\varepsilon)``-isomorphism to a genuine C\*-algebra

[`main_isomorphism`](@ref) produces ``v: B \to A``, where
``B = \bigoplus_l M_{d_l}`` is a genuine finite-dimensional C\*-algebra:

```@example tut
v = main_isomorphism(Φ)
```

The block sizes of ``B`` are a structural fact (for this qubit channel,
``B = M_2``):

```jldoctest tut1
julia> v = main_isomorphism(Φ);

julia> blocks(v)
1-element Vector{Int64}:
 2
```

## Step 4 — factorize the channel

[`factorize`](@ref) realises the isomorphism and its inverse as quantum channels,
giving ``\Phi \approx \Delta\Upsilon`` through ``B``.

!!! warning "factorize has a tighter domain"
    `factorize` builds ``\Upsilon`` via a power-series functional calculus whose
    convergence domain is much smaller than the `prop_P` basin
    ``\rho(\Phi^2-\Phi) < 1/4``. The verb pre-checks a conservative
    ``\rho < 0.10`` and throws a clean `ArgumentError` for out-of-domain inputs —
    it does **not** abort the session. [`associated_algebra`](@ref) and
    [`main_isomorphism`](@ref) keep the wider ``\rho < 1/4`` domain;
    [`certified_defect`](@ref) is always safe at any ``\eta``. (bug `aic-exa.13`)

Use a genuinely oblique, multi-block fixture so the full machinery is exercised
(``B = M_2 \oplus M_2``):

```@example tut
Ψ = UCPMap(bce_conj_kraus(4, 2, 0.02))
F = factorize(Ψ)
```

The encode and decode channels are the duals ([`CPMap`](@ref), not square
`UCPMap`s):

```@example tut
encode(F), decode(F)
```

The two certified round-trip brackets:

```@example tut
delups_defect(F), upsdel_defect(F)
```

The genuine target algebra and its block structure:

```jldoctest tut2; setup = :(using AlmostIdempotentChannels, LinearAlgebra; e(d,i)=(v=zeros(ComplexF64,d);v[i]=1;v); function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; function fixed_unitary(nn); H=zeros(ComplexF64,nn,nn); for p in 1:(nn-1);H[p,p+1]=0.4;H[p+1,p]=0.4;end; for p in 1:(nn-2);H[p,p+2]=0.3im;H[p+2,p]=-0.3im;end; for p in 1:nn;H[p,p]=0.1*p;end; exp(im*Matrix(Hermitian(H))); end; function bce_conj_kraus(d,m,t); base=block_cond_exp_kraus(d,m);V=fixed_unitary(d);Vd=V'; K=Matrix{ComplexF64}[]; for Ka in base;push!(K,sqrt(1-t)*Ka);end; for Ka in base;push!(K,sqrt(t)*(Vd*Ka*V));end; K; end; Ψ = UCPMap(bce_conj_kraus(4,2,0.02)); F = factorize(Ψ))
julia> algebra(F)
CStarAlgebra  B = ⊕_l M_{d_l}
  blocks d_l = [2, 2]   (m = 2)
  dim_B = 8,  n_B = 4

julia> blocks(algebra(F))
2-element Vector{Int64}:
 2
 2
```

The certified round-trip defects are linear in ``\eta`` and lie below the
``50\eta`` test envelope (`test_factorize.jl`):

![Certified round-trip defects vs η; linear in η, below the 50·η envelope (test_factorize.jl)](../assets/roundtrip.png)

For ``\eta > 0`` the `decode` channel is only ``O(\eta)``-trace-preserving (an
internal Choi→Kraus PSD-cone clip). [`iscptp`](@ref) is `false` at the default
`atol = 1e-9` and `true` once loosened. The rigorous certificate is the cb-norm
round-trip bracket above, not `iscptp` at machine tolerance. The `encode` channel
is TP to machine precision:

```@example tut
(iscptp(encode(F)),
 iscptp(decode(F)),               # false at atol=1e-9
 iscptp(decode(F); atol = 1e-4))  # true once loosened
```

## Where to go next

- [The η = 0 oracle](@ref) — watch every certified defect collapse to machine
  ``\varepsilon`` on exactly-idempotent channels.
- [Multi-block channels](@ref) — a dedicated walkthrough of the oblique
  ``B = M_2 \oplus M_2`` case.
- [The mathematics](@ref) — the four-stage mathematical story with paper anchors.
- [API reference](@ref) — every exported docstring.
