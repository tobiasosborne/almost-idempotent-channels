# Bibliography

Primary sources and the libraries the package builds on.

---

**Kitaev, Alexei.** *Almost-idempotent quantum channels and approximate C\*-algebras.* arXiv:2405.02434, 2025.

This is the paper the package implements. Every algorithm is a constructive finite-dimensional realization of a result in this paper; every source comment cites a line number from its LaTeX source.

```bibtex
@misc{kitaev2025almostidempotent,
  author        = {Kitaev, Alexei},
  title         = {Almost-idempotent quantum channels and approximate {C}*-algebras},
  year          = {2025},
  eprint        = {2405.02434},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph}
}
```

---

**Watrous, John.** *Simpler semidefinite programs for completely bounded norms.* arXiv:1207.5726; published in *Chicago Journal of Theoretical Computer Science*, 2013.

Source of the MAX-primal and MIN-dual semidefinite programs for the diamond norm (``\|\cdot\|_\diamond = \|\cdot\|_{\mathrm{cb}}``) used by `diamond_norm_watrous`, `diamond_norm_watrous_primal`, `diamond_norm_dual`, and `idempotency_defect` in the MOSEK extension.

---

**Choi, Man-Duen.** *Completely positive linear maps on complex matrices.* Linear Algebra and its Applications, 10(3):285–290, 1975.

Source of the Choi matrix characterization of completely positive maps. The Convention-A Choi matrix used throughout this package is the standard object from this paper; it encodes CP as positivity and unitality as a partial-trace constraint.

---

**The Choi–Effros product** ``X \star Y = \tilde\Phi(XY)`` on ``\operatorname{Img}\tilde\Phi`` (see arXiv:2405.02434 §1 and references therein).

The algebra structure on the image of an idempotent UCP map was identified by Choi and Effros; this package implements the approximate version for ``\eta``-idempotent maps. The `EpsCStarAlgebra` type carries this product.

---

**The FLINT project / arb.** *Fast Library for Number Theory*, including the `acb_mat` and `arb` arbitrary-precision ball arithmetic. [https://flintlib.org/](https://flintlib.org/)

The arb path of the C core uses FLINT's certified ball arithmetic (`acb_mat_*`, `arb_*`) for all nontrivial functional-calculus and norm-certification steps. Certified outward-rounded brackets on ``\eta``, ``\varepsilon``, ``\|v\|_{\mathrm{cb}}``, and the factorization round-trip defects are computed in this library.

---

**LAPACK / LAPACKE.** Anderson et al., *LAPACK Users' Guide*, SIAM, 1999. [https://www.netlib.org/lapack/](https://www.netlib.org/lapack/)

The fast double path in the C core uses LAPACK routines (`zheev`, `zgesvd`, `zgees`) for eigendecomposition, singular value decomposition, and Schur form. The double path serves as the speed anchor and the prec=53 cross-check target for the arb path.

---

**License.** This package is distributed under the GNU Affero General Public License v3.0 (AGPL-3.0); see the `LICENSE` file at the repository root.

This is an independent realisation of the algorithms in arXiv:2405.02434. The paper's author has no involvement in this package; any errors in the algorithms or their bounds are ours.
