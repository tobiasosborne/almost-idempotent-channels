# AlmostIdempotentChannels.jl

**Constructive finite-dimensional algorithms for Kitaev's almost-idempotent
quantum channels and approximate C\*-algebras — with certified arbitrary-precision
error bounds throughout.**

A near-idempotent quantum channel secretly hides a clean digital code. Feed in
messy analog Kraus operators and this package extracts the exact block structure
``\bigoplus_l M_{d_l}`` of the genuine C\*-algebra the channel is approximating,
builds the encode/decode channels that move between them, and returns a
*certified* bound — a proof of an inequality, not a floating-point guess — on how
lossy the round trip is. It is a Julia surface over the `libaic` C/arb cores,
implementing
[Kitaev, *Almost-idempotent quantum channels and approximate C\*-algebras*,
arXiv:2405.02434](https://arxiv.org/abs/2405.02434) (2025) as finite-dimensional
algorithms for *every* result the paper proves non-constructively.

## Why this exists

**The object.** A unital completely positive (UCP) map ``\Phi`` on ``B(H)``,
``\dim H = n``, in the Heisenberg (observable) picture
``\Phi(X) = \sum_a K_a^\dagger X K_a``. This is the most general noisy,
information-preserving measurement of observables you can write down. ``\Phi`` is
**``\eta``-idempotent** when ``\|\Phi^2 - \Phi\|_{cb} \le \eta``, the
completely-bounded norm — a supremum over all ampliations ``1 \otimes \Phi``, not
the bare operator norm. ``\eta`` measures how close ``\Phi`` is to a genuine
projection onto a stable subalgebra of observables; ``\eta = 0`` exactly when
``\Phi`` is a conditional expectation onto a C\*-subalgebra.

**What Kitaev proved, and why it is surprising.** An ``\eta``-idempotent
``\Phi`` does not merely *almost* project. Its image carries the structure of an
**``\varepsilon``-C\*-algebra** ``A = \operatorname{Img}\tilde\Phi`` with
``\varepsilon = O(\eta)``: a vector space satisfying the C\*-axioms only *up to*
``\varepsilon``, with the Choi–Effros product ``X \star Y = \tilde\Phi(XY)`` and
no exact unit. The headline rigidity theorem (`th_main`, arXiv:2405.02434 §2) is
that **every finite-dimensional ``\varepsilon``-C\*-algebra is
``O(\varepsilon)``-isomorphic to a *genuine* C\*-algebra
``B = \bigoplus_l M_{d_l}``, with a universal, dimension-independent constant.**
The constant does not grow with ``n`` — and the naive averaging route fails
precisely because *its* constant grows like ``n`` (§9). When ``A`` comes from a
finite-dimensional ``\Phi``, that isomorphism and its inverse are realised by
quantum channels, giving an **approximate factorization
``\Phi \approx \Delta\Upsilon``** of the channel through ``B``: a decode map
``\Delta`` and an encode map ``\Upsilon`` with
``\Delta\Upsilon \approx \tilde\Phi`` and ``\Upsilon\Delta \approx 1_B``.

**What this package contributes.** Kitaev's proofs are non-constructive — the
Lefschetz–Hopf fixed-point theorem, holomorphic functional calculus via contour
integrals, Haar-measure diagonals. None of them hands you the object. In finite
dimensions every one of those objects is computable, and this package implements
a constructive algorithm for each, wrapped as five verbs that read like the paper.

Every defect comes back as a [`CertifiedBracket`](@ref) ``lo \le x \le hi``:
FLINT/arb error balls rounded *outward* so the bound survives the conversion to
`double`. The output is a proof of the same inequality the theorem promises. The
package is **solver-free by default** — the rigorous ``\eta`` bracket and the
entire factorization run on the eig-free arb certifier with no SDP solver
installed.

## The headline pipeline

The whole verb chain, executed live at documentation-build time (so this output
can never rot):

```@example pipeline
using AlmostIdempotentChannels, LinearAlgebra

# An in-domain near-idempotent qubit channel:
# Φ_t = (1−t)·id + t·(maximally depolarizing), t = 0.1.
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

```@example pipeline
η = certified_defect(Φ)        # rigorous ‖Φ²−Φ‖_cb bracket — NO solver
```

```@example pipeline
A = associated_algebra(Φ)      # the ε-C* algebra A = Img Φ̃
```

```@example pipeline
v = main_isomorphism(Φ)        # the O(ε)-isomorphism v: B → A, dim-independent constant
```

`certified_defect` is always safe at any ``\eta``;
[`factorize`](@ref) has a tighter, small-``\eta`` domain (see the
[Tutorial](@ref) and the [Known limits](#Known-limits) below).

## How to read the results

Every scalar comes back as a [`CertifiedBracket`](@ref) `lo ≤ x ≤ hi` you can
test with `∈`:

```@example pipeline
0.0 in certified_defect(Φ)   # the defect is non-negative; is 0 in the bracket?
```

```@example pipeline
width(η), midpoint(η)        # how loose the certificate is, and its centre
```

The bracket is loose by design (`hi/lo ~ 2n`) — it *certifies* a value rather
than competing with the solver. The optional MOSEK extension adds the exact
diamond-norm value ([`idempotency_defect`](@ref)) and a tight bracket
(`certified_defect(Φ; tight=true)`); without a solver both throw a helpful
install hint rather than a `MethodError`.

## The five verbs

Each row is a theorem of the paper made executable, and a rigorous certificate of
its bound.

| Result (arXiv:2405.02434) | Verb / type | Certified by |
|---|---|---|
| **Idempotency defect** ``\eta = \|\Phi^2-\Phi\|_{cb}`` (§1) | [`certified_defect`](@ref) → [`CertifiedBracket`](@ref) | eig-free arb bracket, outward-rounded; solver-free |
| **`prop_P` regularisation** ``\tilde\Phi = \theta(2\Phi-1)`` (§3) | [`associated_algebra`](@ref) → [`EpsCStarAlgebra`](@ref) | certified ``\varepsilon`` associator-defect bracket |
| **`th_main` isomorphism** ``v: B \to A``, ``O(\varepsilon)``, dim-independent (§2/§9) | [`main_isomorphism`](@ref) → [`MainIsomorphism`](@ref) | certified ``\|v\|_{cb}``, ``\|v^{-1}\|_{cb}``, isodefect brackets |
| **`th_idemp_structure`** the genuine ``B = \bigoplus_l M_{d_l}`` (§11) | [`algebra`](@ref) → [`CStarAlgebra`](@ref) | exact block sizes from the Wedderburn structure |
| **`th_factorization`** ``\Phi \approx \Delta\Upsilon`` (§12) | [`factorize`](@ref) → [`ChannelFactorization`](@ref) | two certified round-trip brackets |
| **Encode / decode channels** ``\Upsilon^*``, ``\Delta^*`` (§1, §12) | [`encode`](@ref), [`decode`](@ref) → [`CPMap`](@ref) | CPTP duals; round-trip brackets are the certificate |

The [Tutorial](@ref) walks the whole chain — including the ``\eta = 0`` oracle,
the cleanest ground truth in the paper. The [API reference](@ref) renders every
exported docstring.

## Known limits

Stated plainly:

- **`factorize` has a tighter domain than the rest of the pipeline.** It builds
  ``\Upsilon`` via a power-series functional calculus whose convergence domain is
  much smaller than the `prop_P` basin ``\rho(\Phi^2-\Phi) < 1/4``. The verb
  pre-checks a conservative ``\rho < 0.10`` in Julia and throws a clean
  `ArgumentError` — it does **not** abort the session.
  [`associated_algebra`](@ref) and [`main_isomorphism`](@ref) keep the wider
  ``\rho < 1/4`` domain; [`certified_defect`](@ref) is always safe at any
  ``\eta``.
- **The `decode` channel is only ``O(\eta)``-trace-preserving for ``\eta > 0``.**
  Its internal Choi→Kraus PSD-cone projection clips a small negative mass
  (``\approx 3.7\times10^{-6}`` on multi-block oblique fixtures); so
  `iscptp(decode(F))` is `false` at the default `atol = 1e-9` and `true` at
  `atol = 1e-4`. The rigorous certificate is the cb-norm round-trip bracket,
  not `iscptp` at machine tolerance. The encode channel and the ``\eta = 0``
  oracle decode are TP to machine precision.
- **The eig-free bracket is loose by design** (`hi/lo ~ 2n`). For a tight bracket
  or the exact value, use the MOSEK extension.

## Citation

This is an independent realisation of the constructions in:

> Alexei Kitaev, *Almost-idempotent quantum channels and approximate C\*-algebras*,
> arXiv:2405.02434 (2025).

The paper's author has no involvement in this package; any errors in the
algorithms or their bounds are ours alone. The certified arithmetic uses
[FLINT/arb](https://flintlib.org/); the fast double path uses LAPACK/LAPACKE.
The package is licensed under the GNU Affero General Public License v3.0.
