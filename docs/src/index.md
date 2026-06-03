# AlmostIdempotentChannels.jl

**Constructive finite-dimensional algorithms for Kitaev's almost-idempotent
quantum channels and approximate ``C^*``-algebras — with certified
arbitrary-precision error bounds throughout.**

A near-idempotent quantum channel secretly hides a clean digital code. Feed in
messy analog Kraus operators and this package extracts the exact block structure
``\bigoplus_l M_{d_l}`` of the genuine ``C^*``-algebra the channel is approximating,
builds the encode/decode channels that move between them, and returns a *certified*
bound — a proof of an inequality, not a floating-point guess — on how lossy the
round trip is. It is a Julia surface over the `libaic` C/arb cores, implementing
[Kitaev, *Almost-idempotent quantum channels and approximate ``C^*``-algebras*,
arXiv:2405.02434](https://arxiv.org/abs/2405.02434) (2025) as finite-dimensional
algorithms for *every* result the paper proves non-constructively.

```@raw html
<p align="center">
  <img src="assets/demo.svg" alt="The pipeline: certified_defect, associated_algebra, main_isomorphism, factorize" width="760"/>
</p>
```

## The story in three sentences

A unital completely positive map ``\Phi`` that is *almost* idempotent
(``\|\Phi^2-\Phi\|_{\mathrm{cb}} \le \eta``) does not merely almost-project — its
image carries the structure of an **``\varepsilon``-``C^*``-algebra**, an algebra
obeying the ``C^*`` axioms only up to ``\varepsilon = O(\eta)``. Kitaev's rigidity
theorem says every finite-dimensional ``\varepsilon``-``C^*``-algebra is
``O(\varepsilon)``-isomorphic to a *genuine* ``C^*``-algebra
``B = \bigoplus_l M_{d_l}``, with a **universal, dimension-independent constant** —
and when the algebra comes from a finite-dimensional ``\Phi``, that isomorphism is
realised by quantum channels, factoring ``\Phi \approx \Delta\Upsilon`` through
``B``. This package turns each of those non-constructive proofs into a finite
algorithm and certifies every bound it produces with a FLINT/arb error ball.

## The five verbs

Each row is a theorem of the paper made executable, and a rigorous certificate of
its bound. The verbs read like the paper — nouns are the math objects, verbs are
the constructions.

| Result (arXiv:2405.02434) | Verb / type | Certified by |
|---|---|---|
| **Idempotency defect** ``\eta = \|\Phi^2-\Phi\|_{\mathrm{cb}}`` (§1) | [`certified_defect`](@ref) → [`CertifiedBracket`](@ref) | eig-free arb bracket, outward-rounded; solver-free |
| **`prop_P` regularisation** ``\tilde\Phi = \theta(2\Phi-1)`` (§3) | [`associated_algebra`](@ref) → [`EpsCStarAlgebra`](@ref) | certified ``\varepsilon`` associator-defect bracket |
| **`th_main` isomorphism** ``v: B \to A``, ``O(\varepsilon)``, dim-independent (§2/§9) | [`main_isomorphism`](@ref) → [`MainIsomorphism`](@ref) | certified ``\|v\|_{\mathrm{cb}}``, ``\|v^{-1}\|_{\mathrm{cb}}`` brackets |
| **`th_idemp_structure`** the genuine ``B = \bigoplus_l M_{d_l}`` (§11) | [`algebra`](@ref) → [`CStarAlgebra`](@ref) | exact block sizes from the Wedderburn structure |
| **`th_factorization`** ``\Phi \approx \Delta\Upsilon`` (§12) | [`factorize`](@ref) → [`ChannelFactorization`](@ref) | two certified round-trip brackets |
| **Encode / decode channels** ``\Upsilon^*``, ``\Delta^*`` (§1, §12) | [`encode`](@ref), [`decode`](@ref) → [`CPMap`](@ref) | CPTP duals; round-trip brackets are the certificate |

The package is **solver-free by default**: the rigorous ``\eta`` bracket and the
entire factorization run on the eig-free arb certifier with no SDP solver
installed. The optional MOSEK extension adds the exact diamond-norm value
([`idempotency_defect`](@ref)) and a tight bracket.

## Start here

Three routes through the documentation, by what you want:

- **"I want to use the package now."** → [Installation](@ref) → [Quick start](@ref)
  → [The five-verb pipeline](@ref). Five minutes to a certified bound.
- **"I read the paper and want the constructive story."** → [The mathematics](@ref)
  → [The constructive mandate](@ref) → [Dimension-independence](@ref). The
  non-constructive proofs, made finite.
- **"I care whether the numbers are trustworthy."** →
  [Certified arithmetic](@ref) (the four-rung cross-check ladder) →
  [The η = 0 oracle](@ref) (the cleanest ground truth in the paper).

The [How-to guides](@ref "Certify the idempotency defect (solver-free)") are task-oriented
recipes; the [API reference](@ref) renders every exported docstring; the
[Notation glossary](@ref) maps every symbol to its paper line.

## Reading the results

Every scalar comes back as a [`CertifiedBracket`](@ref) ``\mathrm{lo} \le x \le
\mathrm{hi}`` — a proof of that inequality, the FLINT/arb error ball rounded
*outward* so the bound survives the conversion to `double`. It is loose by design
(``\mathrm{hi}/\mathrm{lo} \sim 2n``): it certifies a value rather than competing
with a solver. See [Interpret a CertifiedBracket](@ref) and
[Certified arithmetic](@ref).

## Known limits

Stated plainly; the reasons live in [Design decisions and known limits](@ref).

!!! warning "Two hard constraints"
    - **`factorize` has a tighter domain than the rest of the pipeline.** It
      pre-checks a conservative ``\rho(\Phi^2-\Phi) < 0.10`` and throws a clean
      `ArgumentError` for out-of-domain inputs — it does not abort the session.
      [`associated_algebra`](@ref) and [`main_isomorphism`](@ref) keep the wider
      ``\rho < 1/4`` domain; [`certified_defect`](@ref) is always safe at any ``\eta``.
    - **The `decode` channel is only ``O(\eta)``-trace-preserving for ``\eta > 0``.**
      `iscptp(decode(F))` is `false` at the default `atol = 1e-9` and `true` at
      `atol = 1e-4`. The rigorous certificate is the cb-norm round-trip bracket,
      not `iscptp` at machine tolerance.

!!! note "The eig-free bracket is loose by design"
    The default bracket has ``\mathrm{hi}/\mathrm{lo} \sim 2n``. For a tight bracket
    (``\sim 10^{-13}``) or the exact value, use the [MOSEK extension](@ref "Install and use the MOSEK extension").

## Citation

This is an independent realisation of the constructions in:

> Alexei Kitaev, *Almost-idempotent quantum channels and approximate ``C^*``-algebras*,
> arXiv:2405.02434 (2025).

The paper's author has no involvement in this package; any errors in the algorithms
or their bounds are ours alone. The certified arithmetic uses
[FLINT/arb](https://flintlib.org/); the fast double path uses LAPACK/LAPACKE. The
package is licensed under the GNU Affero General Public License v3.0. See the
[Bibliography](@ref) for the full reference list.
