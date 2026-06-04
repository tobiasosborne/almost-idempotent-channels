# The mathematics

This page explains *what* the package computes and *why*, following

> Alexei Kitaev, _Almost-idempotent quantum channels and approximate
> C\*-algebras_, [arXiv:2405.02434](https://arxiv.org/abs/2405.02434) (2025).

After reading it you will understand the four-stage story the five verbs walk:
a near-idempotent channel, regularised into an exact idempotent, whose image is
an *approximate* ``C^*``-algebra, which is rigidly close to a *genuine* one,
through which the channel approximately factors. This is the understanding-page;
to *run* the pipeline see [The five-verb pipeline](../tutorials/pipeline.md), and to look up a symbol
see the [Notation glossary](../reference/notation.md).

Every nontrivial claim is anchored to a line of the canonical source
`paper/src/approximate_algebras.tex`, cited as `tex:NNN`. Subtleties the
project's own reviews keep catching are cited to `paper/FINDINGS.md`. The package
implements **constructive finite-dimensional algorithms**: the paper's proofs are
often non-explicit — a Lefschetz–Hopf fixed-point count with no witness, or
objects written infinite-dimensionally (a contour-integral functional calculus, a
Haar-measure diagonal) — but in finite dimensions the objects those proofs merely
assert to exist are computable. For each stage we separate **the
paper's proof technique** from **the constructive route the package takes** — and
the paper's ``O(\varepsilon)`` bound is the algorithm's correctness specification.

## The four stages at a glance

| Stage | Object | Key result |
|---|---|---|
| 0 | a UCP map ``\Phi`` and its idempotency defect ``\eta = \|\Phi^2-\Phi\|_{\mathrm{cb}}`` | the completely-bounded norm (`tex:347`) |
| 1 | the regularisation ``\tilde\Phi = \theta(2\Phi-1)`` and its image ``A`` | `prop_P` (§3, `tex:524`) |
| 2 | ``A`` is an ``\varepsilon``-``C^*``-algebra, ``\varepsilon = O(\eta)`` | `th_almost_idemp` (§12, `tex:2192`) |
| 3 | ``A`` is ``O(\varepsilon)``-isomorphic to a *genuine* ``B = \bigoplus_l M_{d_l}`` | `th_main` (§2, `tex:460`) |
| 4 | the channel factors: ``\Phi \approx \Delta\Upsilon`` through ``B`` | `th_factorization` (§12, `tex:2730`) |

## Stage 0 — the object and its defect

The package works in the **Heisenberg (observable) picture**: a *unital
completely positive* (UCP) map

```math
\Phi\colon \mathcal{B}(\mathcal{H}) \to \mathcal{B}(\mathcal{H}),
\qquad \Phi(X) = \sum_a K_a^\dagger X K_a, \qquad \sum_a K_a K_a^\dagger = 1,
```

with ``\dim\mathcal{H} = n``. This is the most general noisy,
information-preserving transformation of observables. States live in the dual; a
**quantum channel** is a CPTP map of density matrices, dual to a UCP map of
observables (`tex:277`). The two pictures are dual, and the direction matters —
the encode/decode *channels* of Stage 4 are the duals of the observable-side maps
(`tex:334`); getting the dual backwards silently transposes everything.

The single quantity that governs the whole theory is the **completely bounded
norm** (`tex:347`):

```math
\|\Lambda\|_{\mathrm{cb}} = \sup_{m}\ \sup_{X \neq 0}
\frac{\|(1_{M_m}\otimes\Lambda)(X)\|}{\|X\|}.
```

It is a supremum over **all ampliations** ``1_{M_m}\otimes\Lambda``, not the bare
operator norm; its dual on functionals is the diamond norm ``\|\cdot\|_\diamond``
(`tex:352`). A map ``\Phi`` is

```math
\Phi \text{ is } \eta\text{-idempotent} \iff \|\Phi^2 - \Phi\|_{\mathrm{cb}} \le \eta
\qquad (\texttt{tex:354}).
```

**Why the cb-norm and not the operator norm.** The hypothesis must survive
*extending the physical system with new, non-interacting parts* — exactly an
ampliation ``1_{M_m}\otimes\Phi``, which the cb-norm sees and the operator norm
does not (`tex:346`). This is what lets every statement in the paper hold
*uniformly in ``m``* (`tex:398`), and it is why the tensor extensions of §10 exist
at all. The practical consequence: the spectral radius ``\rho(\Phi^2-\Phi)``, the
operator norm, and ``\eta = \|\Phi^2-\Phi\|_{\mathrm{cb}}`` are **three different
numbers** that coincide only when the defect superoperator is normal
(`FINDINGS §C17`). The package estimates ``\eta`` from the Choi matrix, never from
``\rho`` alone — see [Certified arithmetic](certified_arithmetic.md).

At ``\eta = 0``, ``\Phi`` is *exactly* a conditional expectation onto a
``C^*``-subalgebra. That is the cleanest ground truth in the whole subject, and
the package's strongest cross-check — the [The η = 0 oracle](../tutorials/eta0_oracle.md).

## Stage 1 — regularisation: from near-idempotent to exactly idempotent

To extract an algebra you need a genuinely *idempotent* map. The paper applies
functional calculus to ``\Phi`` itself (`tex:356`, `tex:2171`):

```math
\tilde\Phi = \theta(2\Phi - 1) = \tfrac12\bigl(1 + \operatorname{sgn}(2\Phi-1)\bigr),
\qquad \theta(x) = \tfrac12(1 + \operatorname{sgn} x).
```

The defining series in ``4(\Phi-\Phi^2)`` **converges when ``\eta < 1/4``**
(`tex:2176`), and then ``\tilde\Phi`` is *exactly* idempotent with
``\tilde\Phi^2 = \tilde\Phi``, ``\|\tilde\Phi - \Phi\|_{\mathrm{cb}} \le O(\eta)``,
``\tilde\Phi(1) = 1`` and ``\tilde\Phi(X^\dagger) = \tilde\Phi(X)^\dagger``
(`tex:2177`). One caveat carried through the rest of the theory: ``\tilde\Phi`` is
**not guaranteed completely positive** (`tex:363`, with an explicit
2-dimensional counterexample at `tex:365`) — which is why the *channel*
factorization of Stage 4 is a separate, more delicate construction than the
algebra.

The engine is **Proposition `prop_P`** (`tex:524`): if a Hermitian ``P`` satisfies
``\|P^2 - P\| \le \delta < 1/4`` then ``\tilde P = \theta(2P-I)`` is an exact
projection, commutes with ``P``, and ``\|\tilde P - P\| \le \|2P-I\|\cdot O(\delta)``.
The functional calculus has a **precondition**: ``\operatorname{sgn}(X)``,
``|X|``, ``\theta(X)`` are defined only for ``\|X^2 - I\| < 1`` (`tex:516`,
`tex:520`); the package asserts it and fails loud outside it.

The algebra is then the image, with the **Choi–Effros product**:

```math
A = \operatorname{Img}\tilde\Phi = \operatorname{Ker}(1 - \tilde\Phi),
\qquad X \star Y = \tilde\Phi(XY)\quad(\texttt{tex:341},\ \texttt{tex:2188}).
```

**Warning — The product is ``\tilde\Phi(XY)``, not ``XY``.** The multiplication on ``A`` is the Choi–Effros star ``X \star Y =
\tilde\Phi(XY)``, **not** the plain operator product, which would leave
``\operatorname{Img}\tilde\Phi``. The ``\eta = 0`` identity-channel oracle is
*structurally blind* to this (there ``\tilde\Phi = \mathrm{id}``, so
``\star`` is plain multiplication); every claim about the star must be checked
on a genuinely oblique ``\eta > 0`` fixture (`FINDINGS §C2`).

**Paper's technique vs the constructive route.** The paper frames
``\operatorname{sgn}`` through the resolvent contour integral but explicitly
*declines* to rely on it, because the spectrum of a non-normal operator is
perturbation-fragile — its own ``3\times3`` example has eigenvalues moving like
``t^{1/3}`` (`tex:535`–`tex:544`). The package instead computes ``\tilde\Phi`` by
a certified matrix-sign iteration (Newton–Schulz / Denman–Beavers) gated by an
a-posteriori certificate, with the basin guarded by the **eig-free Gelfand
spectral radius** ``\rho(\Phi^2-\Phi) < 1/4`` rather than the stricter operator
norm (`FINDINGS §C15`). Two facts about ``A`` are load-bearing downstream:

- ``\operatorname{Img}\Phi = \operatorname{Ker}(1-\Phi)`` holds only when ``\Phi``
  is *exactly* idempotent (`tex:344`); for the ``\eta``-idempotent input the
  identity uses the *regularised* ``\tilde\Phi`` — that is the entire reason for
  Stage 1.
- ``A = \operatorname{Img}\tilde\Phi`` is an **oblique** image: ``\tilde\Phi`` is
  Hermiticity-preserving but not Hilbert–Schmidt self-adjoint, so the
  Frobenius-orthogonal projector onto ``A`` does *not* respect the star structure.
  Project into ``A`` with ``\tilde\Phi`` itself (`FINDINGS §C3`).

## Stage 2 — the approximate algebra (``\varepsilon = O(\eta)``)

An **``\varepsilon``-``C^*``-algebra** satisfies the ``C^*`` axioms only *up to*
``\varepsilon`` (`tex:407`–`tex:439`):

```math
\|XY\| \le (1+\varepsilon)\|X\|\,\|Y\|, \qquad
\|(XY)Z - X(YZ)\| \le \varepsilon\,\|X\|\,\|Y\|\,\|Z\|, \qquad
\|X^\dagger X\| \ge (1-\varepsilon)\|X\|^2,
```

and the unit laws ``\|XI - X\| \le \varepsilon\|X\|`` hold only approximately —
there is **no exact unit** by default; an exact unit becomes available only after
the ``O(\varepsilon)``-change of `prop_unit` (`tex:441`, `tex:672`).

**Note — The axioms hold only up to ``\varepsilon``.** Associativity, submultiplicativity, the ``C^*`` identity, and the unit laws
are all approximate. No routine on ``A`` may assume exact associativity or an
exact unit. For a compressed subalgebra the unit is not ``1_n`` but a
projection ``\tilde P`` (`FINDINGS §C7`).

**Theorem `th_almost_idemp`** (`tex:2192`) is that ``A = \operatorname{Img}
\tilde\Phi``, with the inherited norm and involution and the star product, is an
*extended* ``O(\eta)``-``C^*``-algebra — "extended" meaning the same holds
uniformly for ``M_m \otimes A`` (`tex:2209`). This is exactly where the cb-norm
hypothesis of Stage 0 pays off, and it is why **``\varepsilon = O(\eta)``**
(`tex:398`).

A *separate* small parameter governs maps. A **``\delta``-homomorphism**
``v\colon A' \to A''`` satisfies ``\|v(I)-I\| \le \delta`` and ``\|v(XY) -
v(X)v(Y)\| \le \delta\|X\|\|Y\|`` (`tex:443`); a ``\delta``-isomorphism is a
bijective ``\delta``-inclusion (`tex:455`). Keep ``\varepsilon`` (how far ``A`` is
from a ``C^*``-algebra) and ``\delta`` (how far a map is from multiplicative)
distinct — the main proof juggles several small parameters and never identifies
them.

The package materialises ``A`` as a Frobenius-orthonormal basis of
``\operatorname{Img}\tilde\Phi`` and certifies ``\varepsilon`` directly as a
**certified associator-defect bracket**: a sweep of ``\|(X\star Y)\star Z -
X\star(Y\star Z)\|`` over the basis, using the star throughout. That bracket,
divided by ``\eta``, is the executable form of ``\varepsilon = O(\eta)``.

## Stage 3 — the rigidity theorem ``\texttt{th\_main}``

> **`th_main`** (`tex:460`). Any finite-dimensional ``\varepsilon``-``C^*``-algebra
> ``A`` is ``O(\varepsilon)``-isomorphic to some genuine ``C^*``-algebra ``B``.
> The implicit constant **does not depend on ``A`` or its dimension.**

``B`` is an honest finite-dimensional ``C^*``-algebra, ``B = \bigoplus_l M_{d_l}``
(Wedderburn structure, `tex:257`), and the map ``v\colon B \to A`` is a
``\delta``-isomorphism with ``\delta = O(\varepsilon)`` and a **universal**
constant. This is the surprising part: a *quantitative, dimension-free* rigidity.
An algebra that obeys the ``C^*`` axioms only to within ``\varepsilon`` is within
``O(\varepsilon)`` of an honest one, with a constant that does not blow up as the
matrices grow.

**Why the naive route fails — and why this matters for the tests.** The
cohomological idea (`tex:464`–`tex:484`) is to trivialise the associativity defect
as a coboundary using a *diagonal* ``D = \int dU\,(U^\dagger\otimes U)`` (Haar
average). But:

> "naive constructions of the Haar measure (or just the diagonal) in the
> ``\varepsilon``-associative setting have error bounds proportional to
> ``n = \dim A``. So the outlined procedure works only if ``\varepsilon < c\,n^{-1}``."
> (`tex:484`)

That is dimension-*dependent* — the opposite of the theorem. So any honest test of
`th_main` must check that the constant **does not grow with ``n``**; a test that
merely confirms an ``O(\varepsilon)``-isomorphism *exists* for one fixed ``n`` is a
test that cannot fail. This is the subject of [Dimension-independence](dim_independence.md), and
a constant that grows with dimension would be a genuine stop condition, not a bug
to patch.

**Paper's technique vs the constructive route.** The proof (deferred to §9,
`tex:1414`) builds ``B`` incrementally: find a maximal commutative subalgebra whose
projections are one-dimensional; the first split ``\mathbb{C} \to \mathbb{C}\oplus
\mathbb{C}`` needs a **nontrivial projection** in ``A``, whose *existence*
(`lem_nontriv_projection`, `tex:931`) is proved by a pure-existence argument — the
**Lefschetz–Hopf fixed-point theorem** and H-space homotopy (§6) — that exhibits no
witness. Each non-explicit ingredient is replaced by a finite-dimensional algorithm
meeting the same bound (see [The constructive mandate](constructive.md)):

- **Nontrivial projection** ← an *ambient spectral split*: take a non-scalar
  Hermitian ``H \in A`` with the largest interior eigenvalue gap, form
  ``X = \operatorname{sgn}(s(H - tI))`` in ``M_n``, set ``P = \tfrac12(I+X)``, and
  project into ``A`` with ``\tilde\Phi``. This avoids the unitary-group machinery
  and stays in ``M_n`` where the calculus is valid (`FINDINGS §B1`). The gap
  audition must be **unit-aware** (`FINDINGS §C16`).
- **Haar diagonal** ← an explicit Heisenberg–Weyl closed form
  ``D = \sum_{jk} d^{-2} S_{jk}^\dagger \otimes S_{jk}`` with generalized Paulis —
  exact and eig-free (`FINDINGS §B2`).
- **The inverse-function lemma `lem_invfun`** is *already* constructive: the Banach
  contraction iteration ``x_{n} = g_y(x_{n-1})`` is the algorithm (`tex:564`).

The error-control engine is `cor_improvement` (`tex:1317`): a ``\delta``-inclusion
of an exact ``C^*``-algebra into ``A`` is improvable to a ``\delta_0``-inclusion
with ``\delta_0 = O(\varepsilon)`` *independent of ``\delta``* — a Newton-type
contraction down to the ``O(\varepsilon)`` floor.

## Stage 4 — the exact structure and the channel factorization

### The exact case ``\eta = 0`` (``\texttt{th\_idemp\_structure}``)

> **`th_idemp_structure`** (statement `tex:318`, proof `tex:2055`). For an
> *idempotent* UCP map ``\Phi`` there are a subspace ``\mathcal{M} \subseteq
> \mathcal{H}``, a genuine ``C^*``-algebra ``A``, and UCP maps ``\Delta`` and
> ``\Gamma`` with ``\Gamma\operatorname{Co}_{\mathcal{M}}\Delta = 1_A`` and
> ``\Delta\Gamma\operatorname{Co}_{\mathcal{M}} = \Phi`` (`tex:328`).

This is the ``\eta = 0`` oracle: the whole construction collapses to the exact
Choi–Effros decomposition — a genuine ``C^*``-algebra with **zero defect**.
Combined with `prop_hom_structure` (`tex:259`), the image is identified as
``\bigoplus_j \mathcal{B}(\mathcal{L}_j)`` with explicit isometries, and the
conditional-expectation density matrices are uniquely extractable (`FINDINGS §C20`).
It is the cleanest cross-check in the paper and the discipline test for every
feature: if a construction does not give zero on ``\eta = 0`` input, something is
wrong (see [The η = 0 oracle](../tutorials/eta0_oracle.md)).

### The almost-idempotent case (``\texttt{th\_factorization}``)

> **`th_factorization`** (`tex:2730`). For ``\Phi`` with ``\|\Phi^2-\Phi\|_{\mathrm{cb}}
> \le \eta`` there are a finite-dimensional ``C^*``-algebra ``B`` and UCP maps
> ``\Delta\colon B \to \mathcal{B}(\mathcal{H})``, ``\Upsilon\colon
> \mathcal{B}(\mathcal{H}) \to B`` with
> ``\|\Delta\Upsilon - \Phi\|_{\mathrm{cb}} \le O(\eta)`` (`tex:2734`) and
> ``\|\Upsilon\Delta - 1_B\|_{\mathrm{cb}} \le O(\eta)`` (`tex:2739`).

So ``\Phi`` approximately factors *through* the genuine algebra ``B``:
``\Upsilon`` **encodes** ``\mathcal{B}(\mathcal{H})`` into ``B``, ``\Delta``
**decodes** ``B`` back, with ``\Delta\Upsilon \approx \tilde\Phi`` and
``\Upsilon\Delta \approx 1_B``. The duals ``\mathrm{Dec} = \Delta^*`` and
``\mathrm{Enc} = \Upsilon^*`` are the approximate decode/encode **channels** of
states. The package computes all four objects and certifies the two round-trip
inequalities as brackets.

The paper gives this as an *outline* (`tex:2742`, `FINDINGS §D4`): every object —
the non-UCP ``\tilde\Delta = v``, ``\tilde\Upsilon = v^{-1}\tilde\Phi``, then their
CP-ization via per-block unitary 1-designs and `lem_RC` — is an explicit
finite-dimensional matrix expression, so the package's gap is *closure*, not
*constructivity*. Two findings bite here, and both surface as documented limits:

- **``\Delta'`` is only ``O(\eta)``-completely-positive**, not exactly CP, because
  ``v`` is only an extended ``O(\eta)``-isomorphism (the star is not the ordinary
  product). The package projects onto the CP cone with a tolerance-aware
  Choi``\to``Kraus extraction that clips the ``O(\eta^2)`` defect but still aborts
  on a genuine ``O(\eta)`` or ``O(1)`` negative (`FINDINGS §C14`). The consequence
  is that the `decode` channel is only ``O(\eta)``-trace-preserving for
  ``\eta > 0`` — the rigorous certificate is the cb-norm round-trip bracket, not
  `iscptp` at machine tolerance.
- **The round-trip target is the conditional expectation ``P_B``, not the full
  identity ``1_B``.** ``\Delta`` reads only ``B``'s block-diagonal coordinates, so
  ``\Upsilon\Delta = P_B`` (zero on off-block units); comparing the round trip to
  the full ``1_B`` would read ``1`` on every off-block unit even at ``\eta = 0``
  (`FINDINGS §C14`).

**Warning — ``factorize`` has a tighter domain than the rest of the pipeline.** ``\Upsilon`` is built by a power-series functional calculus whose convergence
domain is much smaller than the prop_P basin ``\rho(\Phi^2-\Phi) < 1/4``. The
`factorize` verb pre-checks a conservative ``\rho < 0.10`` and throws a clean
`ArgumentError` rather than aborting (bug `aic-exa.13`). `associated_algebra`
and `main_isomorphism` keep the wider ``\rho < 1/4`` domain;
`certified_defect` is always safe at any ``\eta``. See [Stay in factorize's
domain](../reference/api.md) and [Design decisions and known limits](design_limits.md).

## What remains open

Two project-level questions are honestly open, certified *per instance* rather than
proven a priori (consistent with the paper, not a proof of it):

1. Whether a dimension-independent spectral gap ``\Omega(1)`` always exists for
   ``\dim A > 1`` — the role Lefschetz–Hopf plays in the paper (`FINDINGS §D1`).
2. The explicit value of the universal constant ``c_0`` from `cor_improvement`
   (`FINDINGS §D2`).

The package fails loud on a degenerate spectrum and measures the constant bounded
with no dimensional trend; it does not claim the a-priori guarantee.

## Where to go next

- Run it: [The five-verb pipeline](../tutorials/pipeline.md) and [The η = 0 oracle](../tutorials/eta0_oracle.md).
- Trust it: [Certified arithmetic](certified_arithmetic.md) (the cross-check ladder).
- The central claim: [Dimension-independence](dim_independence.md).
- Look up a symbol: the [Notation glossary](../reference/notation.md).
- The constructivisation in detail: [The constructive mandate](constructive.md).
