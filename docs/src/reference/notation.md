# Notation glossary

Look up any symbol used in the docs and the paper, with its meaning, paper line, and Julia type or function.

For the full mathematical story behind these symbols, see [The mathematics](../explanation/math_story.md).

| Symbol | Meaning | Paper anchor | Julia type / function |
|--------|---------|-------------|----------------------|
| ``\Phi`` | UCP map ``\mathcal{B}(\mathcal{H}) \to \mathcal{B}(\mathcal{H})`` on observables (Heisenberg picture), ``\Phi(X) = \sum_a K_a^\dagger X K_a``, ``\Phi(1) = 1`` | `approximate_algebras.tex:316`, `tex:1568` | `UCPMap` |
| ``\eta`` | Idempotency defect: ``\|\Phi^2 - \Phi\|_{\mathrm{cb}} \le \eta``; require ``\eta < 1/4`` | `approximate_algebras.tex:354` | `certified_defect`, `CertifiedBracket` |
| ``\tilde\Phi`` | Regularised idempotent ``\theta(2\Phi - 1)``; exactly idempotent, ``\|\tilde\Phi - \Phi\|_{\mathrm{cb}} \le O(\eta)`` | `approximate_algebras.tex:356`, `tex:2171` | `associated_algebra` (builds ``\mathcal{A} = \operatorname{Img}\tilde\Phi``) |
| ``\|\cdot\|_{\mathrm{cb}}`` | Completely bounded norm ``\sup_n \sup_{X \neq 0} \|(1_{M_n} \otimes \Lambda)X\| / \|X\|``; dual = diamond norm ``\|\cdot\|_\diamond`` | `approximate_algebras.tex:347` | `certified_defect`, `eta_eigfree`, `choi_diff` |
| ``\rho`` | Gelfand spectral radius (eig-free); basin guard ``\rho(\Phi^2 - \Phi) < 1/4`` | `approximate_algebras.tex:520`, FINDINGS §C15 | `associated_algebra` (pre-checks ``\rho < 0.24``); `factorize` (pre-checks ``\rho < 0.10``) |
| ``\operatorname{sgn}``, ``\theta`` | ``\operatorname{sgn}(X) = X(X^2)^{-1/2}``; ``\theta(X) = \tfrac12(I + \operatorname{sgn}\,X)``; require ``\|X^2 - I\| < 1`` | `approximate_algebras.tex:514`, `tex:528` | — (C core `aic_prop_P`) |
| ``\star`` | Choi–Effros product ``X \star Y = \tilde\Phi(XY)`` (not plain ``XY``) on ``\mathcal{A}`` | `approximate_algebras.tex:341`, `tex:2188` | `EpsCStarAlgebra` (`basis` elements; do not multiply directly in Julia) |
| ``\mathcal{A}`` | The ``\varepsilon``-C* algebra ``\operatorname{Img}\tilde\Phi = \operatorname{Ker}(1 - \tilde\Phi)`` (oblique image) | `approximate_algebras.tex:391`, `tex:2185` | `EpsCStarAlgebra`, `associated_algebra` |
| ``\varepsilon`` | How far ``\mathcal{A}`` is from a C* algebra; ``\varepsilon``-C* axioms hold up to ``\varepsilon``; ``\varepsilon = O(\eta)`` | `approximate_algebras.tex:407`, `tex:398` | `epsilon(A::EpsCStarAlgebra)` → `CertifiedBracket` |
| ``\delta`` | How far a map ``v`` is from multiplicative (``\delta``-homomorphism / -inclusion / -isomorphism) | `approximate_algebras.tex:443`, `tex:451`, `tex:455` | — |
| ``\delta_0`` | Fixed reduced inclusion defect, ``\delta_0 = O(\varepsilon)``, independent of ``\delta`` (error reduction) | `approximate_algebras.tex:490`, `tex:1317` | — |
| ``\mathcal{B}`` | The genuine finite-dim C* algebra ``\bigoplus_l M_{d_l}`` that ``\mathcal{A}`` is ``O(\varepsilon)``-isomorphic to | `approximate_algebras.tex:257`, `tex:461` | `CStarAlgebra` |
| ``\bigoplus_l M_{d_l}`` | Block-diagonal matrix algebra with block sizes ``d_l``; ``\dim \mathcal{B} = \sum_l d_l^2``, ``n_B = \sum_l d_l`` | `approximate_algebras.tex:257`, `tex:2771` | `CStarAlgebra`; `blocks`, `dim_B`, `n_B` |
| ``v`` | The ``O(\varepsilon)``-isomorphism ``\mathcal{B} \to \mathcal{A}`` from `th_main`; universal (dimension-independent) constant | `approximate_algebras.tex:400`, `tex:460` | `MainIsomorphism`, `main_isomorphism` |
| ``\Delta`` | Decode (observable side) ``\mathcal{B} \to \mathcal{B}(\mathcal{H})``, UCP; dual ``\mathrm{Dec} = \Delta^*`` | `approximate_algebras.tex:319`, `tex:2730` | `decode(F::ChannelFactorization)` → `CPMap` |
| ``\Upsilon`` | Encode (observable side) ``\mathcal{B}(\mathcal{H}) \to \mathcal{B}``, UCP; dual ``\mathrm{Enc} = \Upsilon^*`` | `approximate_algebras.tex:400`, `tex:2730`, `tex:2897` | `encode(F::ChannelFactorization)` → `CPMap` |
| ``\Gamma`` | UCP ``\mathcal{B}(\mathcal{M}) \to \mathcal{A}`` (conditional expectation) in the exact ``\eta = 0`` case | `approximate_algebras.tex:319`, `tex:2106` | — (exact-oracle path; ``\eta = 0`` only) |
| ``\operatorname{Co}_{\mathcal{M}}`` | Compression ``X \mapsto J_{\mathcal{M}}^\dagger X J_{\mathcal{M}} : \mathcal{B}(\mathcal{H}) \to \mathcal{B}(\mathcal{M})`` | `approximate_algebras.tex:325` | — |
