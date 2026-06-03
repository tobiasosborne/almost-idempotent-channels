# AlmostIdempotentChannels.jl

**Constructive finite-dimensional algorithms for Kitaev's almost-idempotent
quantum channels and approximate C\*-algebras — with certified
arbitrary-precision error bounds throughout.**

A near-idempotent quantum channel secretly hides a clean digital code. Feed in
messy analog Kraus operators and this package extracts the exact block structure
``\bigoplus_l M_{d_l}`` of the genuine ``C^*``-algebra the channel is approximating,
builds the encode/decode channels that move between them, and returns a *certified*
bound — a proof of an inequality, not a floating-point guess — on how lossy the
round trip is. It is a Julia surface over the `libaic` C/arb cores, implementing
Kitaev, _Almost-idempotent quantum channels and approximate C\*-algebras_
([arXiv:2405.02434](https://arxiv.org/abs/2405.02434), 2025) as finite-dimensional
algorithms for *every* result the paper proves non-constructively.

![The pipeline: certified_defect, associated_algebra, main_isomorphism, factorize](assets/demo.svg)

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
| **Idempotency defect** ``\eta = \|\Phi^2-\Phi\|_{\mathrm{cb}}`` (§1) | [`certified_defect`](reference/api.md) → [`CertifiedBracket`](reference/api.md) | eig-free arb bracket, outward-rounded; solver-free |
| **`prop_P` regularisation** ``\tilde\Phi = \theta(2\Phi-1)`` (§3) | [`associated_algebra`](reference/api.md) → [`EpsCStarAlgebra`](reference/api.md) | certified ``\varepsilon`` associator-defect bracket |
| **`th_main` isomorphism** ``v: B \to A``, ``O(\varepsilon)``, dim-independent (§2/§9) | [`main_isomorphism`](reference/api.md) → [`MainIsomorphism`](reference/api.md) | certified ``\|v\|_{\mathrm{cb}}``, ``\|v^{-1}\|_{\mathrm{cb}}`` brackets |
| **`th_idemp_structure`** the genuine ``B = \bigoplus_l M_{d_l}`` (§11) | [`algebra`](reference/api.md) → [`CStarAlgebra`](reference/api.md) | exact block sizes from the Wedderburn structure |
| **`th_factorization`** ``\Phi \approx \Delta\Upsilon`` (§12) | [`factorize`](reference/api.md) → [`ChannelFactorization`](reference/api.md) | two certified round-trip brackets |
| **Encode / decode channels** ``\Upsilon^*``, ``\Delta^*`` (§1, §12) | [`encode`](reference/api.md), [`decode`](reference/api.md) → [`CPMap`](reference/api.md) | CPTP duals; round-trip brackets are the certificate |

The package is **solver-free by default**: the rigorous ``\eta`` bracket and the
entire factorization run on the eig-free arb certifier with no SDP solver
installed. The optional MOSEK extension adds the exact diamond-norm value
([`idempotency_defect`](reference/api.md)) and a tight bracket.

## Start here

Three routes through the documentation, by what you want:

- **"I want to use the package now."** → [Installation](getting_started/install.md) → [Quick start](getting_started/quickstart.md)
  → [The five-verb pipeline](tutorials/pipeline.md). Five minutes to a certified bound.
- **"I read the paper and want the constructive story."** → [The mathematics](explanation/math_story.md)
  → [The constructive mandate](explanation/constructive.md) → [Dimension-independence](explanation/dim_independence.md). The
  non-constructive proofs, made finite.
- **"I care whether the numbers are trustworthy."** →
  [Certified arithmetic](explanation/certified_arithmetic.md) (the four-rung cross-check ladder) →
  [The η = 0 oracle](tutorials/eta0_oracle.md) (the cleanest ground truth in the paper).

The [How-to guides](howto/certify_defect.md) are task-oriented
recipes; the [API reference](reference/api.md) renders every exported docstring; the
[Notation glossary](reference/notation.md) maps every symbol to its paper line.

## Reading the results

Every scalar comes back as a [`CertifiedBracket`](reference/api.md) ``\mathrm{lo} \le x \le
\mathrm{hi}`` — a proof of that inequality, the FLINT/arb error ball rounded
*outward* so the bound survives the conversion to `double`. It is loose by design
(``\mathrm{hi}/\mathrm{lo} \sim 2n``): it certifies a value rather than competing
with a solver. See [Interpret a CertifiedBracket](howto/read_bracket.md) and
[Certified arithmetic](explanation/certified_arithmetic.md).

## Known limits

Stated plainly; the reasons live in [Design decisions and known limits](explanation/design_limits.md).

**Warning — Two hard constraints.**

- **`factorize` has a tighter domain than the rest of the pipeline.** It
  pre-checks a conservative ``\rho(\Phi^2-\Phi) < 0.10`` and throws a clean
  `ArgumentError` for out-of-domain inputs — it does not abort the session.
  [`associated_algebra`](reference/api.md) and [`main_isomorphism`](reference/api.md) keep the wider
  ``\rho < 1/4`` domain; [`certified_defect`](reference/api.md) is always safe at any ``\eta``.
- **The `decode` channel is only ``O(\eta)``-trace-preserving for ``\eta > 0``.**
  `iscptp(decode(F))` is `false` at the default `atol = 1e-9` and `true` at
  `atol = 1e-4`. The rigorous certificate is the cb-norm round-trip bracket,
  not `iscptp` at machine tolerance.

**Note — The eig-free bracket is loose by design.** The default bracket has ``\mathrm{hi}/\mathrm{lo} \sim 2n``. For a tight bracket
(``\sim 10^{-13}``) or the exact value, use the [MOSEK extension](howto/mosek_install.md).

## Citation

This is an independent realisation of the constructions in:

> Alexei Kitaev, _Almost-idempotent quantum channels and approximate C\*-algebras_,
> arXiv:2405.02434 (2025).

The paper's author has no involvement in this package; any errors in the algorithms
or their bounds are ours alone. The certified arithmetic uses
[FLINT/arb](https://flintlib.org/); the fast double path uses LAPACK/LAPACKE. The
package is licensed under the GNU Affero General Public License v3.0. See the
[Bibliography](reference/bibliography.md) for the full reference list.
