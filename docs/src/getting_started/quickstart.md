# Quick start

In five minutes you will build a near-idempotent qubit channel and certify its
idempotency defect with a rigorous, solver-free bracket.

## Build the channel

``\Phi_t = (1-t)\cdot\mathrm{id} + t\cdot\mathrm{Dep}_d`` is a convex mix of
the identity and the maximally depolarizing map on ``\mathbb{C}^d``. It has a
closed-form defect ``\|\Phi_t^2 - \Phi_t\|_{\mathrm{cb}} = t(1-t)\cdot 2(1-1/d^2)``,
so the certified bracket can be checked against an analytic value.

The fixtures `phit_kraus` and the helpers it uses live in the test suite, not
in the package API. Inline them here (pitfall P1 in the style guide):

```@example quickstart
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

Check the structural facts before doing anything numerical:

```jldoctest quickstart_struct
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
n(Φ), nkraus(Φ)   # 1 identity + 2²=4 depolarizing Kraus operators
# output
(2, 5)
```

`n(Φ) = 2` is the ambient Hilbert-space dimension; `nkraus(Φ) = 5` is the
number of Kraus operators.

## Certify the defect

```@example quickstart
η = certified_defect(Φ)
```

`certified_defect` calls the C/arb core, constructs a rigorous bracket from
FLINT/arb error balls rounded outward, and returns a [`CertifiedBracket`](@ref).
No SDP solver is needed.

## Read the bracket

The bracket is the proof: a certified inequality ``lo \le \|\Phi^2-\Phi\|_{\mathrm{cb}} \le hi``,
not a floating-point guess.

The analytic defect for ``t = 0.1``, ``d = 2`` is
``0.1 \cdot 0.9 \cdot 2 \cdot (1 - 1/4) = 0.135``. The bracket contains it,
to within a 1e-7 margin (the test invariant in `test/test_core.jl`):

```@example quickstart
analytic_eta = 0.1 * 0.9 * 2 * (1 - 1/2^2)   # 0.135

minimum(η) ≤ analytic_eta ≤ maximum(η)         # bracket contains the analytic value
```

The bracket geometry:

```@example quickstart
width(η), midpoint(η)
```

!!! note "The eig-free bracket is loose by design"
    The solver-free bracket width satisfies `hi/lo ~ 2n`; it certifies a value
    rather than computing it. For a tight bracket (width ~1e-13 for this fixture)
    or the exact diamond-norm value, use the MOSEK extension.

`value(η)` returns `nothing` on the solver-free path — there is no point
estimate, only the interval. On an exactly-idempotent channel (``\eta = 0``)
the bracket collapses to near machine ``\varepsilon`` — see [The η = 0 oracle](@ref).

## Where to go next

- [The five-verb pipeline](@ref) — `certified_defect`, `associated_algebra`,
  `main_isomorphism`, `factorize`, `encode`/`decode` in sequence on the same
  channel.
- [Interpret a CertifiedBracket](@ref) — what `lo`, `hi`, `width`, `midpoint`,
  `value`, and `∈` mean, and when `value` is non-nothing.
- [The mathematics](@ref) — the paper results each verb implements, with
  `approximate_algebras.tex` anchors.
