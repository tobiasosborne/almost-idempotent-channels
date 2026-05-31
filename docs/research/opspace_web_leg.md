# Web Research: cb-norm and cb-inclusion for the opspace module

**Bead:** aic-zwo (th_main_ext / opspace)
**Date:** 2026-05-31
**Scope:** Best-in-class literature survey for computing the completely bounded
(cb) norm and cb-inclusion of the approximate isomorphism `v: B -> A`, where
`B = ⊕_j M_{d_j}` (a genuine C* algebra) and `A ⊆ M_N` (a subspace of N×N
matrices). `v` is stored by `vE[i] = v(E_i)` (N×N matrices), is
Hermiticity-preserving, and is NOT completely positive.

---

## Executive Summary (≤8 lines)

**Smith's theorem (1983) resolves D3 — YES, with caveats.** For `v: B -> M_N`,
the cb-spectral norm `||v||_cb = sup_n ||1_{M_n} ⊗ v||_op` is attained at
ampliation level N (the codomain/ambient dimension) by Pisier Proposition 1.12
(citing Smith [Sm2]) and Watrous TQI Theorem 3.46 + duality. This converts D3's
conjecture into a theorem for the operator-norm cb-inclusion: `a_n = ||v^{-1}||_cb`
is attained at level n_B = sum_j d_j (the "size" of B, the domain). **The
cb-norm SDP is directly reusable from the existing pipeline — YES.** The
project already solves the Watrous (2009/2012) diamond-norm SDP, which equals
the cb-spectral norm SDP for Hermiticity-preserving maps via the adjoint
duality `||Phi||_cb = ||Phi*||_diamond` (Watrous TQI Section 3.3.4). Computing
`||v||_cb` requires feeding `J(v)` (the Choi matrix of `v`) as the input to
the existing `src/sdp.jl` MAX-primal solver, with one sign change: the
cb-spectral norm SDP maximizes over the unit spectral-norm ball (not the
trace-norm ball). **The recommended route:** Smith-truncate at N_max = N
(forward) and N_max = n_B (inverse) + Watrous cb-spectral-norm SDP (Watrous
2012 Theorem 6 / TQI Section 3.3.4) via the adjoint, certified by the existing
`aic_cbnorm_certify` arb pipeline with minimal adaptation.

---

## 1. Smith's Theorem (the Truncation Theorem — resolves D3?)

### 1.1 The exact statement

The result is due to R.R. Smith (1983). It appears in the literature in two
equivalent forms:

**Primary citation:**
R.R. Smith, "Completely bounded maps between C*-algebras," *Journal of the
London Mathematical Society* (2) **27** (1983) 157–166.
DOI: 10.1112/jlms/s2-27.1.157

**Textbook statement (Pisier, Proposition 1.12):**
G. Pisier, *Introduction to Operator Space Theory*, Cambridge University Press,
2003. ISBN 978-0-521-81165-1. Proposition 1.12 (p. 26):

> Consider E ⊂ B(H) and u: E -> M_N = B(C^N, C^N). Then:
>   ||u||_cb = ||u_N||_{M_N(E) -> M_N(M_N)},
> where u_N = id_{M_N} ⊗ u (the N-fold ampliation).

The proof in Pisier uses: for any n×n matrix (a_{ij}) in M_n(E) and unit
vectors x_1,...,x_n in C^N, there exist scalar matrices b = (b_{jk}) with
||b||_op ≤ 1 and vectors x'_1,...,x'_N in C^N with ||x'_k|| ≤ 1 such that
x_j = sum_k b_{jk} x'_k. Applying this to the n-ampliation reduces any
n-level test matrix to the N-level one.

**Textbook statement (Watrous TQI, Theorem 3.46 + duality):**
J. Watrous, *The Theory of Quantum Information*, Cambridge University Press,
2018. Section 3.3.2 and Section 3.3.4.

- Theorem 3.46 states: for Phi in T(X, Y) and any space Z,
  `||Phi ⊗ 1_{L(Z)}||_1 ≤ |||Phi|||_1`, with equality when `dim(Z) ≥ dim(X)`.
  That is, the cb-trace norm (diamond norm) is attained at ampliation level
  dim(X) = domain dimension, and is not increased by larger ampliation.
  (Corollary 3.47 states the equality `||Phi ⊗ 1_{L(Z)}||_1 = |||Phi|||_1`.)

- By the adjoint duality `|||Phi|||_1 = |||Phi*|||_inf` (Watrous 2009, equation
  (1); Watrous 2012, equation (1)), the cb-spectral norm `|||Phi|||_inf =
  ||Phi||_cb` is attained at `dim(Y)` = codomain dimension (applying
  Theorem 3.46 to `Phi*`).

- **Smith's result in TQI** (p. 200, Theorem 3.62): Watrous's TQI book notes
  "Theorem 3.62 is due to Smith (1983)" — this is the spectral-norm
  characterization `|||Phi|||_1 = inf_(A0,A1) ||A0|| ||A1||` over Stinespring
  pairs. Smith's 1983 paper's main contribution (the one cited as Proposition
  1.12 in Pisier) is specifically the codomain-truncation: for maps INTO M_N,
  the cb-norm is determined by the N-fold ampliation alone.

### 1.2 Exact hypotheses

**Domain:** E is any operator space (a subspace of B(H)). In particular, it
does NOT need to be a C*-algebra or operator system — it applies to operator
spaces. For our v: B -> A ⊆ M_N, B is a C*-algebra (a subspace of M_{n_B}),
which is an operator space; A ⊆ M_N is an operator space. Both sides satisfy
the hypothesis.

**Codomain:** The codomain must be a concrete matrix algebra M_N (or equivalently
a subspace of M_N for the spectral-norm version, since the ambient ampliation
handles the rest via isometric embedding). The hypothesis is: the codomain is
M_N for some specific N (or embedded in M_N).

**No complete positivity required.** Smith's theorem applies to any linear map,
not just CP or HP maps. This is essential since v is not CP.

### 1.3 Application to v: B -> A ⊆ M_N and v^{-1}: A -> B

**Forward map v: B -> A ⊆ M_N.**

Since A ⊆ M_N, the codomain of v (viewed as mapping into M_N) has ambient
dimension N. By Smith (Pisier Proposition 1.12):
```
||v||_cb = ||1_{M_N} ⊗ v||_op     (attained at ampliation level N)
```
This is the operator-norm of a map from M_N ⊗ B ⊆ M_{N*n_B} to M_N ⊗ M_N =
M_{N^2}. The sup over all n is attained at n = N, a finite level.

**Inverse map v^{-1}: A -> B ⊆ M_{n_B}.**

Here B = ⊕_j M_{d_j} is faithfully represented as n_B × n_B block-diagonal
matrices (n_B = sum_j d_j). The codomain of v^{-1} is M_{n_B} (the ambient).
By Smith:
```
||v^{-1}||_cb = ||1_{M_{n_B}} ⊗ v^{-1}||_op    (attained at ampliation level n_B)
```
The cb-inclusion lower bound is:
```
a_cb = 1 / ||v^{-1}||_cb = 1 / (||1_{M_{n_B}} ⊗ v^{-1}||_op)
```
This is a rigorous finite-N certificate.

**D3 verdict (corrected, definitive):** D3's conjecture "n ≤ dim A suffices" is
now a THEOREM (in fact a stronger one): the truncation is N_max = N (ambient
dimension of A) for the forward map, and N_max = n_B = sum_j d_j for the
inverse map. Both are finite. The "for all n" in th_main_ext is genuinely
handled by Smith's theorem, with no further conjecture.

**Subtlety — subspace vs full M_k codomain.** Smith's theorem (Pisier Prop 1.12)
is stated for maps INTO M_N. For v: B -> A where A ⊊ M_N, the cb-norm of v
(mapping into M_N, ignoring that the image lands in A) equals the N-fold
ampliation. The cb-norm measuring the image in the operator-space structure
of A (not M_N) would be different. For the th_main_ext cb-isomorphism bound,
v is measured as a map from B into M_N (the ambient), which is the correct
quantity: the paper's prop_inc_ext measures `||(1_{M_n} ⊗ v)(X)||_op` where
the op-norm is the ambient M_{nN} norm.

---

## 2. The cb-Norm as an SDP

### 2.1 The Watrous SDP for the cb-spectral norm

**Primary citation:**
J. Watrous, "Semidefinite programs for completely bounded norms," *Theory of
Computing* **5** (2009) 217–238. arXiv:0901.4709.

J. Watrous, "Simpler semidefinite programs for completely bounded norms,"
*Chicago Journal of Theoretical Computer Science* **2013** (2013) Article 8.
arXiv:1207.5726.

J. Watrous, *The Theory of Quantum Information*, Cambridge University Press,
2018. Section 3.3.4 (SDP for the diamond norm / cb-trace norm).

**Key structural fact (Watrous 2009, equation around line 207–215; Watrous
2012, equations 1 and following):**
```
|||Phi|||_inf  (cb-spectral norm) = sup_k ||1_{M_k} ⊗ Phi||_op
|||Phi|||_1    (cb-trace norm, diamond norm) = sup_k ||1_{M_k} ⊗ Phi||_1
|||Phi|||_1 = |||Phi*|||_inf      (adjoint duality)
```
The cb-spectral norm is `||Phi||_cb` in standard notation; the diamond norm is
`||Phi||_diamond = |||Phi|||_1`. They are related by:
```
||Phi||_cb = ||Phi*||_diamond    for every linear map Phi.
```

**The Watrous 2012 SDP for the cb-TRACE norm (diamond norm) of a general
(not necessarily CP) linear map Phi: M_n -> M_m (Watrous 2012, Theorem 6,
Section 3.2):**

Given Phi with Choi matrix J(Phi) ∈ L(M_m ⊗ M_n) (using the convention
J(Phi) = sum_{ij} Phi(E_ij) ⊗ E_ij, dimension nm × nm):

**Primal problem** (maximizes |||Phi|||_1):
```
maximize:    Re⟨J(Phi), X⟩        (= (1/2) Re Tr(J(Phi)† X))
subject to:  [1_m ⊗ rho_0    X   ]  >= 0
             [X†          1_m ⊗ rho_1]
             rho_0, rho_1 ∈ D(C^n)    (density operators on the n-dim input space)
             X ∈ L(M_m ⊗ M_n)
```
(The objective is (1/2)(⟨J(Phi), X⟩ + ⟨J(Phi)*, X*⟩) = Re⟨J(Phi), X⟩.)

**Dual problem** (minimizes |||Phi|||_1):
```
minimize:    (1/2)||Tr_n(Y_0)||_inf + (1/2)||Tr_n(Y_1)||_inf
subject to:  [Y_0     -J(Phi) ]  >= 0
             [-J(Phi)†   Y_1  ]
             Y_0, Y_1 ∈ Pos(M_m ⊗ M_n)
```
where Tr_n is the partial trace over the right (n-dimensional input) factor,
and ||·||_inf is the operator norm (largest singular value).

Strong duality holds for all Phi; both primal and dual optima are attained
(Watrous 2012, Section 3.1 "Strong duality"). The optimal value is
`|||Phi|||_1 / 2` (i.e., the SDP optimum is `|||Phi|||_1`; the normalization
depends on the exact form of the objective). See note on project normalization
below.

**The SDP for the cb-SPECTRAL norm ||Phi||_cb:**
By adjoint duality `||Phi||_cb = |||Phi|||_inf = |||Phi*|||_1`, one computes:
```
||Phi||_cb = SDP_optval(Phi*)
```
where SDP_optval(·) is the Watrous diamond-norm SDP applied to the adjoint
map Phi*. The Choi matrix of Phi* is J(Phi)† (conjugate-transpose of the
Choi matrix of Phi), up to the transposition convention. Concretely:
```
J(Phi*) = J(Phi)^T    (transpose, not conjugate-transpose)
```
in the convention where J(Phi) = sum_{ij} Phi(E_ij) ⊗ E_ij. (The adjoint
Phi*: M_m -> M_n satisfies Phi*(Y) = sum_{ij} Tr(Phi(E_ij) Y) E_ij;
the Choi matrix is the transpose of J(Phi) in the input-output basis.)

### 2.2 Is this the same SDP family as the project's existing diamond-norm SDP?

YES. The project (`src/sdp.jl`, `include/aic_cbnorm.h`, bead aic-m24) already
solves EXACTLY this Watrous diamond-norm SDP for the Hermiticity-preserving
map Lambda = Phi^2 - Phi. The SDP structure is identical:
- The primal and dual programs are the same (same constraint structure).
- The normalization is `||Lambda||_diamond = (2/n) * SDP_optval` for the
  project's Convention-A Choi (where Choi has trace n), pinned empirically.
- The `aic_cbnorm_certify` arb restoration pipeline applies without change
  to any Choi matrix input.

**To compute ||v||_cb** (the cb-spectral norm of v: B -> A ⊆ M_N):
1. Form J(v) = sum_i v(E_i) ⊗ E_i^T ∈ M_{N*n_B} (the Choi matrix of v,
   where {E_i} is a basis for B ⊆ M_{n_B} — the matrix units). This is
   just the stack of the vE[i] images.
2. Compute J(v*) = J(v)^T (the Choi matrix of the adjoint map v*).
   Alternatively: since v is HP (Hermiticity-preserving), J(v) = J(v)†
   as a block matrix, so J(v*) = J(v)^T = conj(J(v)).
3. Feed J(v*) into `src/sdp.jl` as the input Choi matrix (the same solver
   the project uses for Lambda = Phi^2 - Phi), get SDP_optval.
4. Apply the normalization: `||v||_cb = (2/n_B) * SDP_optval` (where n_B
   is the domain dimension — the input space of v*). Verify the normalization
   on a known oracle (v = exact isomorphism on M_d gives ||v||_cb = 1).

**To compute ||v^{-1}||_cb** (the cb-spectral norm of v^{-1}: A -> B ⊆ M_{n_B}):
Same procedure, with v replaced by v^{-1}: form J(v^{-1}), compute its
adjoint's Choi matrix, feed to the SDP. Note v^{-1} as a linear map on
coordinates can be computed from the coordinate matrix M (the dim_A × dim_B
matrix, `aic_dhom_v_sigma_min`'s M_1) by inverting M.

**What changes relative to the existing `aic_cbnorm.h` pipeline:**
- Input Choi matrix: J(v*) instead of J(Phi^2 - Phi). The matrix size is
  N*n_B × N*n_B (potentially larger than the n^2 × n^2 size for a self-map).
- Normalization: the (2/n) factor uses n_B (domain of v*) rather than n
  (the Phi self-map size). Must verify empirically on the eta=0 oracle.
- The `aic_cbnorm_certify` arb restoration pipeline (convex-combination
  primal + eigenvalue-shift dual) works unchanged, since it operates on
  generic PSD Choi matrices.
- Smith's theorem gives the finite truncation: no need to iterate over
  ampliation levels beyond N (for forward) and n_B (for inverse). This
  removes the "for all n" verification burden entirely.

### 2.3 The Choi matrix of v given vE[i] = v(E_i)

The Choi matrix of v: B -> M_N, with {E_i} = standard matrix units of B
(viewed as a subspace of M_{n_B}) and vE[i] = v(E_i) ∈ M_N:
```
J(v) = sum_i vE[i] ⊗ E_i^T  ∈ M_{N * n_B}
```
In block form: the (k,l) block of J(v) (each of size n_B × n_B) is:
```
J(v)[k,l block] = sum_{i: E_i[k,l] = 1} vE[i] * (delta_{i = k*n_B + l})
```
More directly, J(v) is the nm × nm matrix (m = n_B) where the (p,q) entry
of the (a,b) block is v(E_{ab})[p,q]. This is a direct generalization of the
project's existing `aic_ucp_choi_diff` construction (which builds J(Lambda)
for a self-map Lambda: M_n -> M_n).

For a self-adjoint (HP) map v, J(v) is Hermitian:
```
v HP  ⟺  J(v) = J(v)†.
```
Since v is the approximate *-isomorphism (Hermiticity-preserving), J(v) is
Hermitian, and the project's existing Choi conventions apply directly.

---

## 3. The cb-Inclusion / Min-Stretch

### 3.1 Standard approach: invert v and take ||v^{-1}||_cb

The lower isometry constant at ampliation level n:
```
a_n = inf_{X ≠ 0} ||(1_{M_n} ⊗ v)(X)||_op / ||X||_op
    = 1 / ||(1_{M_n} ⊗ v)^{-1}||_op
    = 1 / ||(1_{M_n} ⊗ v^{-1})||_op    (when v is bijective)
```
The cb-inclusion constant is:
```
a_cb = 1 / ||v^{-1}||_cb    (the cb-spectral norm of v^{-1})
```
By Smith's theorem, `||v^{-1}||_cb` is attained at ampliation level n_B
(the codomain dimension of v^{-1} is B ⊆ M_{n_B}):
```
||v^{-1}||_cb = ||1_{M_{n_B}} ⊗ v^{-1}||_op
```
This is a single operator-norm computation at a FIXED finite ampliation level
n_B, not a sup over all n.

### 3.2 Direct SDP for the min-stretch

There is no standard separate SDP for the min-stretch; the standard route is:
1. Compute v^{-1} as a linear map (invert the coordinate matrix M_1 of v,
   or solve the linear system v(X) = Y for each basis element of A).
2. Compute `||v^{-1}||_cb` via the Watrous cb-spectral-norm SDP (Section 2).
3. Return `a_cb = 1 / ||v^{-1}||_cb`.

**Frobenius lower bound (already in the project, NOT operator-norm a_cb):**
The project's `aic_dhom_v_sigma_min` computes `sigma_min(M_1)` (the Frobenius
inclusion inf, FINDINGS §C6). This is a sound collapse-detector but is NOT
the operator-norm a_cb. The relationship is:
```
sigma_min(M_1) / sqrt(dim_A)  ≤  a_cb (operator norm)  ≤  sigma_min(M_1)
```
(norm equivalence, dimension-dependent). The Frobenius bound therefore gives
a valid LOWER bound on a_cb, but with a sqrt(dim_A) gap. For the th_main_ext
cb-isomorphism certification, the operator-norm a_cb (via SDP) is needed.

**Confirmed finding from FINDINGS §C12 (opspace_design.md correction):** The
Frobenius sigma_min route produces a vacuous ampliation test because
sigma_min(I_{n^2} ⊗ M_1) = sigma_min(M_1) trivially. The operator-norm route
(HOPM or SDP via v^{-1}) is the faithful implementation.

### 3.3 Confirmation via Smith

Since `||v^{-1}||_cb = ||1_{M_{n_B}} ⊗ v^{-1}||_op` by Smith (attained at n_B),
the min-stretch `a_cb` is the smallest singular value of the linear map
`1_{M_{n_B}} ⊗ v^{-1}: M_{n_B} ⊗ A -> M_{n_B} ⊗ B`. In coordinate matrix
form, this is the operator-norm min-stretch of the n_B-fold ampliated
coordinate matrix, which has a concrete O(N^2 * n_B^2) matrix representation
and can be measured by `aic_mat_opnorm` / LAPACK singular values on that
explicit matrix.

---

## 4. HOPM Alternative (Fallback if SDP Doesn't Fit)

The higher-order power method (HOPM) / alternating-maximization approach
(Johnston, Kribs, Paulsen 2009, arXiv:0711.3636; QETLAB CBNorm documentation)
computes the cb-norm via alternating optimization over the ampliated unit ball.
For the cb-spectral norm `||v||_cb`, HOPM iterates:
```
u_{k+1} = argmax_{||u||_op ≤ 1} ||(1_n ⊗ v)(u)||_op
```
The ecstar module (bead aic-knm) already implements a HOPM-style defect search
(`aic_ecstar_hopm_*`) for the associativity defect of the ε-C* algebra. The
same pattern would apply here. Key limitations:

- HOPM gives only a LOWER bound on ||v||_cb (a witness, not a proof of the
  maximum). A certified UPPER bound requires the SDP dual or Smith's theorem.
- HOPM may not converge to the global maximum (non-convex); the SDP gives the
  global optimum directly.
- HOPM is useful as a fast initial estimate or double-path check before solving
  the SDP.

For the project's purposes: HOPM is the double-path surrogate (fast, uncertified
lower bound), and the Watrous SDP + arb-certify is the arb-path certified value.
This matches the existing architecture (ecstar Cycle 2 HOPM + deferred Cycle 3
SDP, bead aic-0at). The `aic_ecstar_hopm_*` kernel can be adapted for the
ampliated operator-norm optimization, reusing the `aic_mat_opnorm` / LAPACK
substrate.

---

## 5. Recommendation

Given the existing project infrastructure (Julia+MOSEK SDP + `aic_cbnorm_certify`
arb-restoration pipeline, already certifying the Watrous diamond-norm SDP for
`||Phi^2 - Phi||_cb`), the recommended route for a certified cb-isomorphism
bound for v: B -> A is:

**Route (i): Smith-truncate + Watrous cb-spectral-norm SDP on J(v*) and J((v^{-1})*), certified via the existing arb pipeline.**

Concretely:
1. **Compute J(v)** from the vE[i] images (a straightforward block-matrix
   construction, analogous to `aic_ucp_choi_diff`).
2. **Compute J(v^{-1})** by inverting the coordinate matrix M_1 (already
   available from `aic_dhom_v_sigma_min`) and assembling the inverse map's
   Choi matrix.
3. **Compute ||v||_cb = ||v*||_diamond** by feeding J(v*) = J(v)^T into
   `src/sdp.jl` (the existing Watrous diamond-norm MAX-primal/MIN-dual
   solver). Smith's theorem guarantees this n_B-fold ampliation is exact.
4. **Compute ||v^{-1}||_cb = ||( v^{-1} )*||_diamond** analogously.
5. **Arb-certify both** via `aic_cbnorm_certify` with the SDP feasible points.
   The normalization `(2/n_domain) * SDP_optval` must be re-verified on the
   eta=0 oracle (where ||v||_cb = 1 exactly).
6. **Report a_cb = 1 / ||v^{-1}||_cb** as the certified cb-inclusion bound.

**What changes from the existing pipeline:**
- The input Choi matrix size is N * n_B × N * n_B (not n^2 × n^2 for a self-map).
- For `v` forward: the adjoint map v*: M_N -> B is NOT a self-map, so the
  solver dimension parameters need updating (n_input = N, n_output = n_B for v*).
- The `aic_cbnorm_certify` arb-restoration logic is dimension-agnostic and
  reusable without change.
- A new helper `aic_opspace_choi_v(v, J_out, prec)` builds J(v) from the
  `aic_dhom_v` struct (a ~50 LOC addition).

**Trade-offs vs alternatives:**
- Route (i) SDP: rigorous, reuses existing solver + certifier, O(1/SDP-gap)
  runtime, certified to any precision by arb. Requires MOSEK (already present).
- Route (ii) HOPM: fast lower bound, no MOSEK needed, matches `aic_ecstar`
  pattern. Gives only LOWER bound; Smith gives the upper bound = the HOPM value
  at the optimal iterate (when HOPM converges). Use as double-path check.
- Route (iii) Direct induction verification: verify `a_{2n} >= a_n/2` for
  n=1,...,n_B in the operator norm. Requires computing operator-norm singular
  values of the ampliated map at each level. More expensive (n_B evaluations
  vs 1 SDP), and gives a sequence of bounds rather than a single certified value.
  Useful as a cross-check on the prop_inc_ext inductive structure.

**The recommended implementation order:**
O1 (new): `aic_opspace_choi_v` — assemble J(v) from `aic_dhom_v`.
O2 (new): `aic_opspace_cbnorm` — call `src/sdp.jl` via Julia ccall shim with
  J(v*), return SDP feasible points; call `aic_cbnorm_certify` for arb ball.
O3 (new): `aic_opspace_cbinclusion` — same for v^{-1}: compute M_1^{-1},
  assemble J((v^{-1})*), certify.
Tests: eta=0 oracle (||v||_cb = 1, a_cb = 1), eta>0 sweep (||v||_cb ~= 1+O(eta),
  a_cb ~= 1-O(eta)), universality canary (bounds flat over dim_A).

---

## Key Citations (full list)

1. **Smith (1983) — cb-norm attained at codomain dimension:**
   R.R. Smith, "Completely bounded maps between C*-algebras," *J. London Math.
   Soc.* (2) **27** (1983) 157–166. DOI: 10.1112/jlms/s2-27.1.157.
   [Primary source; not freely available online.]

2. **Pisier (2003) — Proposition 1.12 (Smith's theorem for operator spaces):**
   G. Pisier, *Introduction to Operator Space Theory*, Cambridge University
   Press, 2003. Sample available at:
   https://assets.cambridge.org/97805218/11651/sample/9780521811651ws.pdf
   Proposition 1.12, p. 26: `||u||_cb = ||u_N||_{M_N(E)->M_N(M_N)}` for
   u: E -> M_N. Proof uses scalar matrix factorization of unit ball vectors.

3. **Watrous (2009) — SDP for cb-norm (original):**
   J. Watrous, "Semidefinite programs for completely bounded norms," *Theory of
   Computing* **5** (2009) 217–238. arXiv:0901.4709.
   https://arxiv.org/abs/0901.4709
   Key facts: cb-spectral norm = |||Phi|||_inf = ||Phi ⊗ 1_{L(Y)}||_op (attained
   at codomain dim(Y)); diamond norm = |||Phi|||_1 = ||Phi ⊗ 1_{L(X)}||_1
   (attained at domain dim(X)); adjoint duality |||Phi|||_1 = |||Phi*|||_inf.

4. **Watrous (2012) — Simpler SDP for cb-norm:**
   J. Watrous, "Simpler semidefinite programs for completely bounded norms,"
   *Chicago J. Theoretical Computer Science* **2013** Article 8. arXiv:1207.5726.
   https://arxiv.org/abs/1207.5726
   Theorem 6 (Section 3.1, p. 10): the cb-trace norm SDP in terms of J(Phi)
   directly (no Stinespring factorization needed):
     Primal: max Re⟨J(Phi), X⟩ s.t. [[1_Y⊗rho_0, X],[X†,1_Y⊗rho_1]] >= 0
     Dual:   min (1/2)(||Tr_X(Y_0)||_inf + ||Tr_X(Y_1)||_inf)
             s.t. [[Y_0,-J(Phi)],[-J(Phi)†,Y_1]] >= 0, Y_0,Y_1 >= 0.
   This is EXACTLY the SDP the project already solves (src/sdp.jl, aic-m24).

5. **Watrous TQI (2018) — Theorem 3.46 (truncation) and Theorem 3.62 (Smith):**
   J. Watrous, *The Theory of Quantum Information*, Cambridge University Press,
   2018. Section 3.3.2 and 3.3.4.
   https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf (draft)
   - Theorem 3.46: cb-trace norm attained at dim(X) (domain).
   - Corollary 3.47: equality holds for any dim(Z) >= dim(X).
   - Theorem 3.62 (Smith): cb-trace norm = inf_{Stinespring} ||A0|| ||A1||.
   - p. 200: "Theorem 3.62 is due to Smith (1983)."
   - Section 3.3.4: cb-spectral norm = |||Phi*|||_1 = ||Phi*||_diamond
     (the completely bounded norm of Phi equals the diamond norm of the adjoint).

6. **Paulsen (2002) — standard CB maps reference:**
   V. Paulsen, *Completely Bounded Maps and Operator Algebras*, Cambridge
   University Press, 2002. ISBN 978-0-521-81694-6.
   https://www.cambridge.org/core/books/completely-bounded-maps-and-operator-algebras
   Contains the full theory of cb-norms and the finite-dimensional truncation.
   (Specific theorem number for Smith's result is cited differently in different
   editions; the statement is Pisier Prop 1.12 = Smith 1983.)

7. **QETLAB CBNorm:**
   https://qetlab.com/CBNorm
   Implements the cb-spectral norm via the identity
   `||Phi||_cb = ||Phi†||_diamond` (adjoint duality), calling the DiamondNorm
   SDP on the adjoint map. This confirms that the standard computational
   approach to ||v||_cb is exactly the project's existing diamond-norm SDP
   applied to v* (the adjoint of v).

8. **Johnston, Kribs, Paulsen (2009) — HOPM algorithm for cb-norm:**
   N. Johnston, D.W. Kribs, V.I. Paulsen, "Computing stabilized norms for
   quantum operations via the theory of completely bounded maps," *Quantum
   Inf. Comput.* **9** (2009) 16–35. arXiv:0711.3636.
   https://arxiv.org/abs/0711.3636
   Presents the alternating-maximization algorithm (HOPM) for the cb-norm,
   with the code available at njohnston.ca. HOPM gives a lower bound; the SDP
   gives a certified upper bound.

---

## Self-correction record

**§3.2 vs opspace_design.md §3.2:** The prior design doc (§3.2, the "ORCHESTRATOR
CORRECTION" section at the end of `docs/research/opspace_design.md`) correctly
identified that the Frobenius sigma_min route is vacuous for the ampliation test.
This web leg confirms: the correct route is the operator-norm cb-inclusion via
the Watrous SDP on J(v^{-1}*), NOT sigma_min of I_{n^2} ⊗ M_1.

**Smith's theorem vs D3 conjecture:** D3 (FINDINGS §D3) was classified
"OPEN (bead aic-2jd), conjecture n ≤ dim A suffices." This web leg proves
the conjecture is a THEOREM (Smith 1983, Pisier Prop 1.12): the truncation is
N_max = N (ambient codomain dim) for the forward map and N_max = n_B for the
inverse map. The conjecture's form "n ≤ dim A" was a conservative estimate;
the correct bound is N (for forward) and n_B (for inverse), both finite.
The D3 bead (aic-2jd) should be updated to RESOLVED-THEOREM status.
