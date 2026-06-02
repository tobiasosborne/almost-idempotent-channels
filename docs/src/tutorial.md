# Tutorial

A runnable walkthrough of the headline story: build a near-idempotent channel,
certify its defect, extract the associated ``\varepsilon``-C\*-algebra, factorize
the channel through a genuine C\*-algebra, and finish with the ``\eta = 0``
oracle — the cleanest ground truth in the paper.

Every code block on this page is executed at documentation-build time. The
float-heavy ones are `@example` blocks (their output is embedded live and is
robust to last-digit float wobble); the structural ones — block sizes, dims, type
reprs — are `jldoctest` blocks whose output is checked exactly.

## Fixtures

The package's test suite ships several channel constructors in `test/fixtures.jl`.
Those are part of the *test harness*, not the package API, so here we inline the
two or three small ones the tutorial uses (Heisenberg-picture Kraus,
``\Phi(X) = \sum_a K_a^\dagger X K_a``):

```@example tut
using AlmostIdempotentChannels, LinearAlgebra

e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)

# Φ_t = (1−t)·id + t·(maximally depolarizing) on C^d. In factorize's domain for
# small t; the defect has the closed form ‖Φ_t²−Φ_t‖_⋄ = t(1−t)·2(1−1/d²).
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

# Block conditional expectation on C^d: first m vs the rest. EXACTLY idempotent
# (η = 0); a multi-block C* algebra B = M_m ⊕ M_{d−m}.
function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end

# A fixed basis-mixing unitary V = exp(iH), used to conjugate an idempotent into
# an OBLIQUE image (so the image of Φ̃ is not a *-subalgebra of M_n on the nose).
function fixed_unitary(nn)
    H = zeros(ComplexF64, nn, nn)
    for p in 1:(nn - 1); H[p, p + 1] = 0.4; H[p + 1, p] = 0.4; end
    for p in 1:(nn - 2); H[p, p + 2] = 0.3im; H[p + 2, p] = -0.3im; end
    for p in 1:nn; H[p, p] = 0.1 * p; end
    return exp(im * Matrix(Hermitian(H)))
end

# bce_conj(d,m,t): block_cond_exp(d,m) convex-mixed with its unitary conjugate.
# A GENUINELY η>0, oblique, multi-block almost-idempotent channel; in factorize's
# small-η domain for small t.
function bce_conj_kraus(d, m, t)
    base = block_cond_exp_kraus(d, m); V = fixed_unitary(d); Vd = V'
    K = Matrix{ComplexF64}[]
    for Ka in base; push!(K, sqrt(1 - t) * Ka); end
    for Ka in base; push!(K, sqrt(t) * (Vd * Ka * V)); end
    return K
end
nothing # hide
```

## Step 1 — build the channel and certify its defect

A [`UCPMap`](@ref) wraps the Kraus operators and (by default) checks Heisenberg
unitality ``\sum_a K_a^\dagger K_a = I``:

```@example tut
Φ = UCPMap(phit_kraus(2, 0.1))
```

Its size and Kraus count are stable structural facts — here is the exact-match
doctest form (one identity Kraus + ``d^2 = 4`` depolarizing Kraus = 5):

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

The idempotency defect ``\eta = \|\Phi^2 - \Phi\|_{cb}`` comes back as a rigorous
[`CertifiedBracket`](@ref) — solver-free, no SDP:

```@example tut
η = certified_defect(Φ)
```

The bracket is the package's rigorous spine. Test containment with `∈`, and read
its loose-by-design width:

```@example tut
(0.0 in η, width(η), midpoint(η))
```

## Step 2 — the associated ``\varepsilon``-C\*-algebra

[`associated_algebra`](@ref) builds ``A = \operatorname{Img}\tilde\Phi`` with
``\tilde\Phi = \theta(2\Phi - 1)``, the exact-idempotent regularisation
(`prop_P`, §3). ``A`` is an *oblique* subspace of ``M_n`` carrying the
Choi–Effros product ``X \star Y = \tilde\Phi(XY)``; the C\*-axioms hold only up
to a certified ``\varepsilon``:

```@example tut
A = associated_algebra(Φ)
```

```@example tut
dim(A), ambient(A), epsilon(A)
```

## Step 3 — the ``O(\varepsilon)``-isomorphism to a genuine C\*-algebra

[`main_isomorphism`](@ref) produces ``v: B \to A``, where
``B = \bigoplus_l M_{d_l}`` is a *genuine* finite-dimensional C\*-algebra. The
certified ``\|v\|_{cb}``, ``\|v^{-1}\|_{cb}`` brackets are ``\approx 1 + O(\varepsilon)``,
and the constant is dimension-independent — the paper's central rigidity claim:

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

[`factorize`](@ref) is the headline. It realises ``v`` and ``v^{-1}`` as quantum
channels, giving ``\Phi \approx \Delta\Upsilon`` through ``B``, with two certified
round-trip brackets ``\|\Delta\Upsilon - \Phi\|_{cb}`` and
``\|\Upsilon\Delta - 1_B\|_{cb}``.

!!! warning "factorize has a tighter domain"
    `factorize` builds ``\Upsilon`` via a power-series functional calculus whose
    convergence domain is much smaller than the `prop_P` basin
    ``\rho(\Phi^2-\Phi) < 1/4``. The verb pre-checks a conservative
    ``\rho < 0.10`` and throws a clean `ArgumentError` for out-of-domain inputs —
    it does **not** abort the session. Use the in-domain fixtures `phit_kraus(2, 0.1)`,
    `bce_conj_kraus(4, 2, 0.02)`, or the ``\eta = 0`` oracles below. A
    larger-``\eta`` channel like `phit_kraus(2, 0.3)` will throw.

Use a genuinely oblique, multi-block fixture so the full machinery is exercised
(``B = M_2 \oplus M_2``):

```@example tut
Ψ = UCPMap(bce_conj_kraus(4, 2, 0.02))
F = factorize(Ψ)
```

The two encode/decode channels are the duals (rectangular CPTP maps, hence
[`CPMap`](@ref), not square `UCPMap`s):

```@example tut
encode(F), decode(F)
```

```@example tut
delups_defect(F), upsdel_defect(F)
```

The genuine target algebra and its block structure are exact:

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

For ``\eta > 0`` the `decode` channel is only ``O(\eta)``-trace-preserving (an
internal Choi→Kraus PSD-cone clip), so [`iscptp`](@ref) is `false` at the default
`atol = 1e-9` and `true` once loosened — the rigorous certificate is the cb-norm
round-trip bracket above, not `iscptp` at machine tolerance. The `encode` channel
is TP to machine precision:

```@example tut
(iscptp(encode(F)),
 iscptp(decode(F)),               # false at atol=1e-9
 iscptp(decode(F); atol = 1e-4))  # true once loosened
```

## Step 5 — the ``\eta = 0`` oracle

When ``\Phi`` is *exactly* a conditional expectation, the whole construction
collapses to the genuine Choi–Effros structure: every certified defect drops to
machine ``\varepsilon`` and the extracted block structure is exactly right. This
is the cleanest cross-check in the paper.

```@example tut
Φ0 = UCPMap(block_cond_exp_kraus(4, 2))   # exactly idempotent: η = 0
F0 = factorize(Φ0)
```

Both round-trip brackets sit at machine ``\varepsilon`` (well below `1e-9`), and
the ``\eta = 0`` decode is exactly trace-preserving:

```@example tut
(maximum(delups_defect(F0)) < 1e-9,
 maximum(upsdel_defect(F0)) < 1e-9,
 iscptp(decode(F0)))
```

The block structure of ``B`` is exact — ``M_2 \oplus M_2`` from a 4-dimensional
conditional expectation onto two 2-blocks:

```jldoctest tut3; setup = :(using AlmostIdempotentChannels, LinearAlgebra; function block_cond_exp_kraus(d,m); K0=zeros(ComplexF64,d,d);K1=zeros(ComplexF64,d,d); for i in 1:m;K0[i,i]=1;end; for i in (m+1):d;K1[i,i]=1;end; Matrix{ComplexF64}[K0,K1]; end; Φ0 = UCPMap(block_cond_exp_kraus(4,2)); F0 = factorize(Φ0))
julia> blocks(algebra(F0))
2-element Vector{Int64}:
 2
 2
```

## Where to go next

- The [API reference](@ref) renders every exported docstring.
- The optional MOSEK extension adds the exact diamond-norm value
  [`idempotency_defect`](@ref) and a tight bracket
  (`certified_defect(Φ; tight = true)`); without a solver both throw a helpful
  install hint rather than a `MethodError`.
