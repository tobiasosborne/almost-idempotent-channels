# The constructive mandate

After this page you will understand what "constructive" means in this package,
which non-constructive steps of the paper it replaces, and why that replacement
is what makes the output checkable.

This page explains why the package replaces the paper's non-constructive proofs
with finite-dimensional algorithms; for the mathematics those algorithms compute,
see [The mathematics](math_story.md). For the claim that the resulting constant does not
grow with the matrix size, see [Dimension-independence](dim_independence.md).

## The mandate

The package exists to answer one question: where the paper proves that an object
*exists* without telling you how to find it, can you find it anyway? In finite
dimensions the answer is yes, and that is the whole point (`CLAUDE.md:39-79`,
verbatim):

> **Many proofs in the paper are non-constructive and phrased for possibly
> infinite-dimensional spaces. But every *theorem* applies to finite-dimensional
> matrices, and in finite dimensions the objects those proofs merely assert to
> exist are computable. Our goal is a constructive algorithm for *every* result.**

So for each result two things are kept strictly separate: the paper's proof
*technique*, and the finite-dimensional algorithm the package *implements*. The
algorithm need not follow the proof — it only needs to land on the same object
and meet the same bound (`CLAUDE.md:148-151`).

## The bound is the specification

The correctness specification is not "matches the proof step by step." It is the
inequality the theorem promises (`CLAUDE.md:74-76`, verbatim):

> **The paper's ``O(\varepsilon)`` bound is the algorithm's correctness specification.**
> The constructive routine must produce an output meeting the bound the theorem
> promises; that bound is the assertion the cross-check verifies.

This is why a constructive replacement is *safer* than a transcription, not just
more usable. Every routine produces a number, and that number is checked against
the theorem's own ``\le O(\eta)`` inequality on adversarial inputs at the
boundary of the hypotheses. A proof you cannot run cannot catch a sign error in
its own constants; a bound you can certify can. The cleanest instance is the
``\eta = 0`` case, where every defect must collapse to machine precision — see
[The η = 0 oracle](../tutorials/eta0_oracle.md).

## The four canonical moves

The paper's non-constructive tools fall into four families. For each, the package
takes a finite-dimensional route that produces the same object
(`CLAUDE.md:58-72`).

| Paper's non-constructive tool | Constructive finite-dim route | Anchor |
|---|---|---|
| **Nontrivial projection** via the Lefschetz–Hopf fixed-point theorem and H-space homotopy (§6): build the approximate unitary group as a manifold and count fixed points of inversion | Spectral split of a Hermitian element of the ambient ``M_n``: pick a non-scalar Hermitian ``H`` with the largest interior eigenvalue gap, form ``X = \operatorname{sgn}(s(H - tI))``, set ``P_{\mathrm{amb}} = \tfrac12(I + X)``, then project into ``\mathcal{A}`` with ``\tilde\Phi`` | `FINDINGS §B1`, `tex:931` |
| **Holomorphic functional calculus** via a contour integral ``f(X)=\tfrac1{2\pi i}\oint_C (zI-X)^{-1}f(z)\,dz`` (§3) | Eigendecomposition; or a Newton / Denman–Beavers sign iteration for ``\operatorname{sgn}(X)``; or a Zolotarev rational approximation, each with a certified arb error | `tex:535`, `tex:548` |
| **Haar-measure diagonal** ``D = \int dU\,(U^\dagger \otimes U)`` (§2/§9) | An explicit finite average over the regular representation: for a matrix block ``M_d`` the closed form ``D = \sum_{jk} d^{-2} S_{jk}^\dagger \otimes S_{jk}`` over the Heisenberg–Weyl operators ``S_{jk} = X^j Z^k`` (generalized Paulis) | `FINDINGS §B2`, `tex:484` |
| **Implicit / inverse function theorem** ``\texttt{lem\_invfun}`` (§3): asserts a solution exists | *Already* constructive: the Banach contraction iteration ``x_n = g_y(x_{n-1})`` **is** the algorithm, with the paper's ``c < 1`` giving the convergence rate ``\|x_n - x_{n-1}\| < r\,c^{n-1}`` | `tex:564`, `tex:591` |

Two of these deserve a sentence of why the naive route is wrong, not just slow.

The paper itself declines the contour-integral calculus, and says why
(`approximate_algebras.tex:535`, verbatim):

> "We have not found use of these methods due to their fragility in the
> approximate setting and the lack of a ``*``-representation, which could
> alleviate it."

The reason is that the spectrum of a non-normal operator is
perturbation-fragile: the paper's own ``3 \times 3`` example has eigenvalues
moving like ``t^{1/3}`` in a ``t``-sized entry (`tex:540`). So the package does
its functional calculus on Hermitian or normal elements in the exact ambient
``M_n``, where the calculus is valid, and certifies everything else with arb
error balls.

The Haar diagonal needs equal care. The literal Cartesian-product transcription
of the per-block formula to a direct sum is **non-central** — the finite Pauli
sum has a nonzero first moment, leaving spurious cross-sector terms. The correct
finite central diagonal is the cross-term-free embedded sum ``D = \sum_l \sum_{jk}
d_l^{-2} (\hat S^{(l)}_{jk})^\dagger \otimes \hat S^{(l)}_{jk}`` with each ``\hat
S^{(l)}`` supported on block ``l`` and zero elsewhere. This is a flagged ``.tex``
discrepancy: the formula at `tex:1254` is wrong for two or more blocks
(`FINDINGS §A2`, `FINDINGS §B2`).

## Two questions that remain open

The mandate is honest about where it certifies per instance rather than proves a
priori. Two project-level questions are **open**, and the package treats both the
same way: it certifies them for the instance in front of it (fail-loud on
degeneracy, measured-bounded constant) and does not claim a general proof.

**The ``\Omega(1)`` spectral gap.** The nontrivial-projection construction needs
a Hermitian element with an interior eigenvalue gap bounded away from zero. The
package certifies that gap *per instance* and aborts fail-loud on a degenerate
spectrum. Whether a dimension-independent gap *always* exists for ``\dim
\mathcal{A} > 1`` is exactly what the paper invokes Lefschetz–Hopf to guarantee,
and an explicit finite-dimensional a-priori guarantee is **open** (`FINDINGS
§D1`). Empirically a genuine ``\ge 2``-dimensional ``\varepsilon``-``C^*`` algebra always
delivers a projection at an ``O(1)`` gap, so the small-gap regime does not arise
at the algebra level — but this is measured, not proven.

**The universal constant ``c_0``.** The rigidity theorem's isomorphism is
``O(\varepsilon)`` with a constant that is universal and dimension-independent
(`approximate_algebras.tex:460`). The package measures ``c = \text{isodefect}/\eta``
per instance and checks it is bounded with no upward trend across a dimension
sweep — measured ``c \le 0.047`` with OLS slope ``-2.7\mathrm{e}{-4}`` over
``\dim \mathcal{A} \in \{8, 12, 16, 18, 20\}``. The *analytic* value of the
universal constant ``c_0`` from the error-reduction lemma `cor_improvement`
remains **open** (`FINDINGS §D2`). The naive averaging route is rejected
precisely because *its* constant grows ``\propto n`` (`tex:484`); a bound that
grows with dimension is therefore a stop condition, not a constant to patch — see
[Dimension-independence](dim_independence.md).

When a result genuinely resists constructivization, that is a finding to escalate
with the specific obstruction, not a corner to paper over (`CLAUDE.md:78-79`).

## See also

- [The mathematics](math_story.md) — the four-stage story the constructive routes realise
- [Dimension-independence](dim_independence.md) — why the constant must not grow with ``n``,
  and how the package checks it
- [The η = 0 oracle](../tutorials/eta0_oracle.md) — the exact case where every constructive
  output collapses to machine precision
