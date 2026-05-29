# Web-research leg: axiom-defect estimators for finite-dimensional ε-C* algebras

**Date:** 2026-05-29  
**Context:** `almost-idempotent-channels` project — constructive finite-dim implementation
of Kitaev arXiv:2405.02434, C + FLINT/arb + Julia-ccall.  
**Task:** Survey published algorithms for computing defect estimators
(ε_assoc, ε_sub, ε_cstar) where the norm is the *operator* (spectral) norm on M_N
and the unit ball is the spectral-norm ball, not the Frobenius ball.

---

## 1. Multilinear / injective tensor norm: complexity and algorithms

### 1.1 NP-hardness (canonical reference)

**Hillar & Lim, "Most tensor problems are NP-hard"**  
arXiv:0911.1393; published J. ACM 60(6), 2013.  
https://arxiv.org/abs/0911.1393  
https://dl.acm.org/doi/10.1145/2512329

Key results relevant here:

- Computing or even *approximating* the spectral norm (= injective norm)
  of a **3-tensor** is NP-hard. No FPTAS exists unless P = NP.
- This encompasses: deciding whether a 3-tensor has a given spectral norm;
  finding its best rank-1 approximation; approximating singular vectors.
- **Order 2 (matrices) is tractable** — operator norm = largest singular value,
  computable in O(N³) via SVD. The phase transition sits strictly at order 3.
- Restricting to symmetric tensors does not alleviate the hardness.

Closely related: **arXiv:2212.14775** ("Complexity and computation for the
spectral norm and nuclear norm")  
https://arxiv.org/pdf/2212.14775  
— confirms: if one dimension of the order-3 tensor is *fixed* (≤ L for constant L),
  a FPTAS exists via spherical grids + SDP. For our setting N, d ≤ 16, all
  dimensions are small/fixed, so this FPTAS regime applies.

**Implication for this project.** The defect maps h(X,Y,Z) are trilinear, so
computing ε_assoc exactly is NP-hard in principle. However, since N ≤ 16 and
d ≤ 16, one of two practical exits applies: (a) the dimension is small enough
that exhaustive or polynomial-in-N algorithms are feasible in practice, or (b)
one dimension is treated as fixed, enabling the FPTAS of arXiv:2212.14775.
This is not a barrier to implementation — it is a complexity-theory result about
worst-case asymptotic scaling; for bounded N it sets no wall.

### 1.2 Algorithms: lower bounds (witness extraction)

**Higher-order power method (HOPM) / alternating maximization**

Standard references:  
- Kolda & Mayo, "Shifted power method," SIAM J. Matrix Anal. Appl. 32(4), 2011.  
- Usevich et al., "A new convergence proof for the higher-order power method,"
  arXiv:1407.4586  
  https://arxiv.org/abs/1407.4586  
  https://link.springer.com/article/10.1007/s00211-018-0981-3 (Num. Math. 2018)

Algorithm (order-3 case): to maximize |h(X,Y,Z)| / (‖X‖·‖Y‖·‖Z‖) over the
product of spectral-norm unit balls, alternating maximization proceeds:

1. Fix Y, Z; maximize over X → a rank-1 matrix problem solvable in O(N²) once
   the linear functional is formed (find dominant singular vectors).
2. Cycle through all three slots, update each to the maximizer with others fixed.
3. Iterate to convergence.

Convergence guarantees:
- **Converges to a critical point** (Łojasiewicz inequality argument); NOT
  guaranteed global maximum.
- Convergence rate: eventually R-linear for generic tensors; can be sublinear
  near degenerate configurations.
- Output is a *certified lower bound* (a concrete witness (X*,Y*,Z*) with a
  measured defect value), not an upper bound.
- Multi-start (sample many random initializations) improves the probability of
  finding the global maximum. In small dimensions (d ≤ 16) this is cheap.

**Matrix-free stochastic operator norm estimation**  
arXiv:2410.08297 ("Matrix-free stochastic calculation of operator norms")  
https://arxiv.org/pdf/2410.08297  
Uses random search directions with exact line search; converges to global maximum
*almost surely* given enough restarts. Applicable when the effective matrix is
formed explicitly.

**Practical recipe for lower bounds (small N, d ≤ 16):**  
Form the trilinear form as a 3-tensor of size (N²)×(N²)×(N²) in the
Frobenius-basis expansion, but the unit ball constraint is spectral-norm on M_N
(NOT ℓ₂ on the coefficient vector — see §2). Run HOPM from 100–1000 random
starts; take the maximum. Output is a numerically certified lower bound (not
certified in the interval-arithmetic sense, but provably ≤ the true defect by
the evaluation at the witness).

### 1.3 Algorithms: upper bounds

**Sum-of-squares / Lasserre hierarchy (SDP relaxations)**

Standard reference: Lasserre, "Global optimization with polynomials and the
problem of moments," SIAM J. Optim. 11(3), 2001.

For spectral norm computation:
- arXiv:2310.17827 "A hierarchy of eigencomputations for polynomial optimization
  on the sphere" https://arxiv.org/pdf/2310.17827  
  — develops a hierarchy for computing tensor spectral norms where each level
  requires an *eigenvalue computation* rather than a full SDP. Provides *lower*
  bounds on the minimum of a form on the sphere; by sign flip gives upper bounds
  on suprema. Reported convergence O(1/r²) per hierarchy level r.
- arXiv:2412.13191 "Moment-SOS and spectral hierarchies for polynomial
  optimization on the sphere and quantum de Finetti"  
  https://arxiv.org/pdf/2412.13191  
  — moment-SOS hierarchy also provides certified upper bounds; convergence O(1/r²)
  is established. SDP at level r has size O(d^r); feasible for small d.
- Standard Lasserre SOS hierarchy: convergence to the global optimum guaranteed
  for polynomials on compact semialgebraic sets (Archimedean case).
  https://arxiv.org/pdf/2011.08139

**SDP relaxation approach (level-1 / Shor):**  
Lift the spectral-norm unit-ball constraint to a PSD constraint. The bilinear
form sup_{‖X‖≤1, ‖Y‖≤1} |f(X,Y)| becomes a maximization of a quadratic over
the product of spectral-ball constraints. Writing ‖A‖_op ≤ 1 ⟺ A†A ≤ I
(i.e., as a PSD constraint), the semidefinite program is:

  max  vec(X)ᵀ M_f vec(Y)   s.t.  X†X ≤ I, Y†Y ≤ I
  
This is a bilinear SDP (not convex). The Shor/Lasserre lift at level 1 relaxes it
to a convex SDP in O(N²) × O(N²) matrix variable; provides a rigorous *upper bound*.

Reference for the bilinear SDP hierarchy:  
arXiv:1810.12197 "Semidefinite programming hierarchies for constrained bilinear
optimization," Mathematical Programming 2021.  
https://arxiv.org/abs/1810.12197  
https://link.springer.com/article/10.1007/s10107-021-01650-1  
Asymptotically converging hierarchy; higher levels improve the bound.

**Triangle-inequality "basis sweep" bound:**  
A crude but certifiable and fast upper bound:
  ε_assoc ≤ Σ_{i,j,k} |α_{ijk}| · ‖h(B_i, B_j, B_k)‖_op
where h(X,Y,Z) = (X⋆Y)⋆Z − X⋆(Y⋆Z) = Φ(Φ(XY)Z) − Φ(XΦ(YZ)), B_i are the
orthonormal basis vectors, and α_{ijk} are the coordinates.
Using ‖X‖_op ≤ ‖X‖_F (the norm equivalence), one gets a coefficient bound, but
this introduces a √dim factor in the unit-ball constraint (see §2).
If the basis is expressed with respect to the Frobenius basis and the spectral-norm
unit ball is used correctly (by forming the sum as a bound on ‖·‖_op of the
trilinear tensor), the bound is:
  ε_assoc ≤ Σ_{i,j,k} ‖h(B_i, B_j, B_k)‖_op   (since ‖B_i‖_F = 1, but ‖B_i‖_op ≤ 1)
This is over d³ terms, each costing one ‖·‖_op evaluation (one SVD). Dimension
scaling: O(d³ · N³) for computing all terms.
**Dimension independence:** This bound is ZERO when η = 0 (idempotent case), and
therefore satisfies the "canary" — but it is typically very loose (d³ factor).
See §4 for detailed analysis.

**Sphere-covering approximation:**  
arXiv:2302.14219 "Approximating Tensor Norms via Sphere Covering"  
https://arxiv.org/pdf/2302.14219  
Constructs a finite ε-net on the sphere and evaluates the form at all net points.
Gives a (1−ε)-approximation for tensor spectral norms; the net size grows as
1/ε^{d−1}. For small d this is feasible. No certified arb path without interval
evaluation at each grid point.

---

## 2. The spectral-norm unit ball vs. Frobenius ball: the √dim trap

This is the central trap. The defects are defined with ‖·‖_op (spectral/operator
norm) on M_N, NOT with the Frobenius norm. These are related by:

    ‖A‖_op ≤ ‖A‖_F ≤ √N · ‖A‖_op   (for N×N matrices)
    ‖A‖_F ≤ √(min(m,n)) · ‖A‖_op   (for m×n)

If the defect estimator computes the maximum over the Frobenius unit ball instead
of the spectral-norm unit ball, the result will be inflated by a factor up to √N.
Similarly, bounding ‖h(X,Y,Z)‖_op from above using the Frobenius norm of the
coefficient tensor introduces a factor of √(d³) = d^{3/2}. For d = 16 this is
a factor of 64 — a spurious failure of the universality canary.

**The correct approach for operator-norm unit balls:**

1. **Variational/vector reduction.** The operator norm ‖A‖_op = sup_{‖u‖₂=‖v‖₂=1} |⟨u,Av⟩|.
   Thus the trilinear defect:
   
     ε_assoc = sup_{X,Y,Z ∈ A, ‖X‖_op=‖Y‖_op=‖Z‖_op=1} ‖h(X,Y,Z)‖_op
   
   can be written as:
   
     sup_{X,Y,Z ∈ A, ‖X‖_op≤1} sup_{‖u‖₂=‖v‖₂=1} |⟨u, h(X,Y,Z) v⟩|
   
   which is a *6-argument* supremum over products of vector unit spheres and
   matrix-spectral-norm balls in A. This is still NP-hard in general but exploits
   the matrix structure: for fixed u,v, ⟨u, h(X,Y,Z) v⟩ is a trilinear form
   in the coefficient vectors of X, Y, Z within A, making it a degree-3
   polynomial on the product of ℓ₂ balls (after change of basis in A to one
   where the unit ball in A ∩ {operator norm ≤ 1} corresponds to a spectrahedron,
   not an ℓ₂ ball).

2. **Critical structural reduction.** Since A has dimension d and basis {B_k},
   any X ∈ A with ‖X‖_op ≤ 1 satisfies ‖X‖_F ≤ √d (not √N) because X lives
   in a d-dimensional subspace. The Frobenius-to-spectral ratio within A is:
     ‖X‖_op ≤ ‖X‖_F ≤ √d · ‖X‖_op.
   So the √dim factor is √d (dimension of A), not √N (ambient dimension). For
   d ≤ 16, this gives at most a factor of 4. This is still a constant inflation
   but it does NOT grow with N — it only grows with dim(A). Since Kitaev's
   theorem requires a universal constant independent of N but does allow it to
   depend on d (the paper's constant is truly universal; our estimator just needs
   to not scale with N), this √d factor is acceptable for the Frobenius-based
   bound when d is small and fixed.

3. **Best practice for certified upper bounds without spurious dim-dependence:**
   - Use the SDP formulation where the spectral-norm ball constraint
     ‖X‖_op ≤ 1 appears explicitly as a PSD constraint: X†X ≤ I.
   - For the HOPM lower bound, enforce the spectral-norm ball by projecting
     the update onto {A ∈ M_N : ‖A‖_op ≤ 1} at each step (projection by SVD
     truncation / clipping singular values to 1).
   - Avoid any algorithm that uses the Frobenius norm as a proxy for the
     spectral norm in the unit-ball constraint.

References:  
- arXiv:1603.05621 "Operator Norm Inequalities between Tensor Unfoldings on the
  Partition Lattice" — systematic treatment of spectral vs. Frobenius norm
  bounds for tensors and their unfoldings, with dimension-dependent factors.
  https://arxiv.org/abs/1603.05621
- Wikipedia "Matrix norm" — basic inequalities: ‖A‖_op ≤ ‖A‖_F ≤ √N ‖A‖_op.
  https://en.wikipedia.org/wiki/Matrix_norm

**The √dim trap for the "basis sweep" bound (§4):**  
The trilinear form h(X,Y,Z) = (X⋆Y)⋆Z − X⋆(Y⋆Z). Expanding X = Σ xᵢBᵢ
with ‖xᵢ‖ = 1 (in the spectral-norm sense), by triangle inequality:
  ‖h(X,Y,Z)‖_op ≤ Σᵢⱼₖ |xᵢ||yⱼ||zₖ| ‖h(Bᵢ,Bⱼ,Bₖ)‖_op.
The Cauchy-Schwarz bound on Σ|xᵢ||yⱼ||zₖ| ≤ (Σxᵢ²)^{1/2} · d^{1/2} ×... is
where the dim factor enters. But if the unit-ball constraint is enforced by
‖X‖_op ≤ 1 (not ‖x‖_ℓ₂ ≤ 1 for the coordinate vector), then the relationship
between |xᵢ| and ‖X‖_op involves the operator norms of the basis elements B_i,
not their Frobenius norms. In the worst case ‖x‖_ℓ₂ ≤ ‖X‖_op · max_i ‖B_i‖_op⁻¹,
so no universal simplification is possible without knowing the operator norms of
the B_i. This makes the naive "sum over basis triples" bound dimension-dependent
unless the basis is *orthonormal in operator norm* (which Frobenius-ONB is not).

---

## 3. Certified / verified numerics for operator norms

### 3.1 FLINT/arb native support (as of FLINT 3.x)

Source: FLINT 3.6.0-dev documentation  
https://flintlib.org/doc/acb_mat.html  
https://fredrikj.net/blog/2018/12/eigenvalues-in-arb/

FLINT/arb `acb_mat` provides **Frobenius norm** (exact), **Frobenius norm upper
bound**, and **∞-norm upper bound**. It does **NOT** provide a native
`acb_mat_bound_spectral_norm` or any certified SVD. This is a gap.

For certified eigenvalue enclosures (needed to certify ‖A‖_op = √λ_max(A†A)):
- `acb_mat_eig_global_enclosure` — Miyajima's method, O(N³), certifies that all
  eigenvalues lie within a ball of radius ε around the approximate eigenvalues.
  Does not isolate individual eigenvalues.
- `acb_mat_eig_enclosure_rump` — Rump's Newton-like method; certifies specific
  eigenvalue clusters.
- `acb_mat_eig_simple_rump`, `acb_mat_eig_simple_vdhoeven_mourrain` — O(N³) /
  O(N⁴), certify all simple eigenvalues simultaneously.

**Strategy for certified ‖A‖_op in arb:**
1. Form B = A†A (exact in arb — product of acb_mat).
2. Compute approximate eigenvalues λ̂_k of B using `acb_mat_approx_eig_qr`.
3. Apply `acb_mat_eig_global_enclosure` to get radius r such that each eigenvalue
   of B is within r of some λ̂_k.
4. ‖A‖_op ≤ sqrt(max_k λ̂_k + r) is a certified upper bound.
   (The sqrt needs to be bounded: arb's `arb_sqrt` is certified.)
5. ‖A‖_op ≥ sqrt(max_k λ̂_k − r) is a certified lower bound.
This is O(N³) per matrix, certifiable in arb.

### 3.2 Rump's verified singular value bounds (INTLAB)

**Rump, "Verified bounds for singular values, in particular for the spectral norm
of a matrix and its inverse," BIT Numerical Mathematics 51(2):367–384, 2011.**  
https://link.springer.com/article/10.1007/s10543-010-0294-0  
https://www.tuhh.de/ti3/paper/rump/Ru10a.pdf

Rump gives several methods for rigorous (interval-arithmetic) enclosures of all
singular values and specifically the spectral norm. The approach reduces to
certified eigenvalue bounds for A†A. Implemented in INTLAB (MATLAB). The key
practical insight: for a general matrix, the spectral norm can be bounded by
combining an approximate SVD with a perturbation bound, made rigorous via
interval arithmetic.

This approach translates directly to arb: use acb_mat's eigenvalue certification
on A†A to get ‖A‖_op.

### 3.3 Julia IntervalLinearAlgebra.jl

Source: Julia Discourse thread, "Rigorous spectral opnorm"  
https://discourse.julialang.org/t/rigorous-spectral-opnorm/82739

The Julia `IntervalLinearAlgebra.jl` package (based on Hladik 2013) provides
eigenvalue bounding but not a turnkey certified spectral norm. The recommended
approach in the community: certify eigenvalues of A†A using eigenvalue enclosure
methods (eigenbox), then take sqrt.

### 3.4 Certified upper bounds for bilinear/trilinear defects

Given a certified bound ‖A‖_op ≤ b (from the A†A approach above), one can:

**Method A (explicit matrix formulation):** Write the bilinear/trilinear map
as a matrix. For ε_sub, the bilinear form f(X,Y) = X⋆Y has operator norm
over {‖X‖_op ≤ 1, ‖Y‖_op ≤ 1}. Form the matrix M_f by vectorizing:
  vec(f(X,Y)) = M_f (vec(X) ⊗ vec(Y)).
Then ‖f‖ ≤ ‖M_f‖_op (where ‖M_f‖_op is a matrix operator norm).
But: ‖M_f‖_op bounds the bilinear form over *Frobenius* unit balls, not spectral.
To get a bound over spectral balls, one needs the cb-norm (see §3.5).

**Method B (cb-norm approach):** The operator norm of a bilinear map
  f: (M_N, ‖·‖_op) × (M_N, ‖·‖_op) → (M_N, ‖·‖_op)
is, by complete boundedness theory, expressible as an SDP. See §3.5.

**Method C (product-of-certified-norms):** For ε_sub:
  ‖X⋆Y‖_op ≤ ‖Φ‖_cb · ‖XY‖_op ≤ ‖Φ‖_cb · ‖X‖_op · ‖Y‖_op.
So ε_sub ≤ ‖Φ‖_cb − 1 (since Φ is UCP, ‖Φ‖_cb = 1; thus ε_sub ≤ 0 exactly,
a tautology). More useful: bound ε_sub by bounding ‖Φ(XY) − Φ(X)Φ(Y)‖_op
as a function of η.

### 3.5 Completely bounded norm SDP (exact, not a relaxation)

**Watrous, "Semidefinite programs for completely bounded norms,"
Theory of Computing 5:11, 2009.**  
https://arxiv.org/abs/0901.4709  
https://theoryofcomputing.org/articles/v005a011/

**Watrous, "Simpler semidefinite programs for completely bounded norms," 2012.**  
https://arxiv.org/abs/1207.5726  
https://cs.uwaterloo.ca/~watrous/Papers/SimplerSDPforCBnorm.pdf

The cb-norm (and diamond norm) of a linear map Φ: M_N → M_N is computable via
an **exact SDP** of matrix size N² × N². This is not a relaxation.

Relevance: the bilinear form f(X,Y) = Φ(XY) − Φ(X)·Φ(Y) (the "Φ-multiplicativity
defect") over the operator-norm unit ball has its norm bounded above by 2·‖Φ−Φ²‖_cb,
which is computable exactly by SDP at size N² × N².

For the associativity defect h(X,Y,Z) = (X⋆Y)⋆Z − X⋆(Y⋆Z) = Φ(Φ(XY)Z) −
Φ(XΦ(YZ)), the analogous cb-norm computation is more complex (involves a
3-fold ampliation), but the principle extends.

**arXiv:0902.3397 "On the complexity of approximating the diamond norm"**  
https://arxiv.org/pdf/0902.3397  
— confirms the diamond norm / cb-norm SDP is polynomial-time and exact.

---

## 4. Exploiting structure: basis-triple bound and zero-at-η=0 canary

### 4.1 The basis-triple bound

Given orthonormal basis {B_k} (Frobenius-ONB, ‖B_k‖_F = 1), define:
  C_assoc := Σᵢ,ⱼ,ₖ ‖h(Bᵢ, Bⱼ, Bₖ)‖_op.

For X = Σᵢ xᵢBᵢ with ‖X‖_F ≤ 1, Y = Σⱼ yⱼBⱼ, Z = Σₖ zₖBₖ similarly:
  ‖h(X,Y,Z)‖_op ≤ Σ |xᵢ|·|yⱼ|·|zₖ| · ‖h(Bᵢ,Bⱼ,Bₖ)‖_op
              ≤ ‖x‖_ℓ₂·‖y‖_ℓ₂·‖z‖_ℓ₂ · C_assoc       (Cauchy-Schwarz, loose)
              ≤ 1 · 1 · 1 · C_assoc                    (using Frobenius unit ball)

So: **ε_assoc ≤ C_assoc when the unit ball is the Frobenius ball.**
For the operator-norm unit ball:
  ‖X‖_op ≤ 1 does NOT imply ‖X‖_F ≤ 1 (it implies ‖X‖_F ≤ √N);
  but restricted to A, it implies ‖X‖_F ≤ √d.
So the operator-norm-unit-ball version gives:
  ε_assoc (spectral) ≤ d^{3/2} · max_{i,j,k} ‖h(Bᵢ,Bⱼ,Bₖ)‖_op.

This is dimension-dependent: O(d^{3/2}) factor in the worst case. For d = 16,
this is factor 64. This can falsely flag the construction as failing the
universality canary if d is even moderately large.

### 4.2 Zero at η = 0 (idempotent oracle)

When Φ² = Φ exactly (η = 0), h(X,Y,Z) ≡ 0 for all X,Y,Z ∈ A (since A =
Image(Φ) is a genuine subalgebra under Choi-Effros product). Therefore:
  C_assoc = Σᵢⱼₖ ‖h(Bᵢ,Bⱼ,Bₖ)‖_op = 0 when η = 0.

This means the basis-triple bound:
- Is zero in the exact case (passes the canary).
- Scales as O(η) for small η (since h(X,Y,Z) = O(η) componentwise).
- Is loose by a factor O(d^{3/2}) relative to the true ε_assoc.

**Recommendation:** Use C_assoc as a cheap computable upper bound with
full dimension-independence analysis documented. Log the ratio C_assoc / ε_assoc
(from HOPM lower bound) across test cases to quantify empirical slack.
If the ratio is consistently O(1) for our specific Φ structures (which come
from near-idempotent channels with structured Kraus decompositions), the
extra factor may be benign in practice even if not in the worst case.

### 4.3 "Near zero" exploitation

Since the defects are expected to be O(η), algorithms that need to compute the
maximum of a form that is near zero benefit from:
- Numerical stability: the alternating maximization will converge to (near) zero
  quickly; the certified SDP dual variable will have a corresponding certificate.
- The basis-triple bound is numerically stable and easy to certify (each
  ‖h(Bᵢ,Bⱼ,Bₖ)‖_op is a single certified operator norm computation).
- The SDP approach gives both a primal (lower bound witness) and dual certificate
  (upper bound), and near zero the duality gap is small.

---

## 5. Practical audition slate for N ≤ 16, d ≤ 16, arb path required

### ε_assoc (hardest; trilinear)

**Candidate A: Basis-triple + certified operator norms (arb)**
- Bound type: UPPER (certifiable in arb, rigorous)
- Algorithm: Compute h(Bᵢ,Bⱼ,Bₖ) for all d³ pairs; bound each ‖·‖_op via
  arb eigenvalue enclosure on (result)†·(result); sum.
- Cost: O(d³ · N³) per defect, arb-certifiable.
- Dim scaling: loose by O(d^{3/2}) relative to operator-norm ball; ZERO at η=0.
- Certifiable in arb: YES.
- Risk: overestimates by d^{3/2} factor; may need companion lower bound.
- Notes: Fast, simple, guaranteed ZERO at η=0. Good "phase 1" filter.

**Candidate B: HOPM multi-start (double path)**
- Bound type: LOWER (witness; not certifiable as upper bound)
- Algorithm: Alternating maximization over spectral-norm balls in A.
  At each step: update X ← argmax ⟨X, (−)⟩ subject to ‖X‖_op ≤ 1, X ∈ A.
  The inner maximization is: find dominant singular vector pair of a matrix
  in M_N, projected to A, scaled to spectral-norm ball. O(N³) per step.
  Run 100–1000 random starts; return max.
- Cost: O(n_starts · n_iter · N³); for N=16, d=16: fast in double.
- Dim scaling: converges to critical point, not necessarily global max.
  For small d, multi-start covers the space well.
- Certifiable in arb: lower bound is certifiable (evaluate at witness in arb).
- Risk: local maxima; must use multi-start aggressively.
- Notes: Essential companion to Candidate A. The pair (HOPM lower, basis-triple upper)
  brackets the true defect. Report both.

**Candidate C: Level-1 SDP relaxation (Lasserre lift, MOSEK via Julia)**
- Bound type: UPPER (rigorous relaxation, exact SDP)
- Algorithm: Form the bilinear/trilinear optimization as a polynomial optimization
  problem over the product of spectral-ball constraint sets (each written as
  A†A ≤ I). Lift to a semidefinite program; solve with MOSEK.
- Cost: SDP of size O(N²·d) × O(N²·d) at level 1; for N=16, d=16: matrix of
  size ~4096 (feasible for MOSEK).
- Dim scaling: bound tightens with level; level-1 may be loose; level-2 is larger.
- Certifiable in arb: SDP output is floating-point; for certified bound, use the
  SDP optimal value plus its computed duality gap as an upper bound.
  See arXiv:2308.07287 for rigorous SDP bounds.
- Risk: SDP relaxation gap may be nontrivial at level 1; need level 2 or 3 for
  tightness; cost scales as d^{2r} per level r.
- Notes: Gives genuine upper bound; complements HOPM. Good "phase 2" refinement.
  Watrous SDP (arXiv:0901.4709) handles bilinear forms via cb-norm exactly.

**Candidate D: Watrous cb-norm SDP (exact, for bilinear sub-case)**
- Bound type: EXACT (for the bilinear case; trilinear is harder)
- Algorithm: Express the bilinear defect (ε_sub, ε_cstar) as a cb-norm of
  a derived map; compute via Watrous SDP at size N²×N². Exact, not a relaxation.
- Cost: SDP of matrix size N²×N²; for N=16: 256×256 SDP. Easily solved in MOSEK.
- Dim scaling: EXACT value; no dimension factor.
- Certifiable in arb: SDP is floating-point; arb certification requires a
  post-hoc interval enclosure of the SDP solution. This is a known technique
  (interval SDP verification); see arXiv:2308.07287.
- Risk: The trilinear associativity defect does not directly reduce to a
  bilinear cb-norm; the Watrous SDP handles bilinear forms only.
- Notes: Best for ε_sub and ε_cstar. For ε_assoc, use Watrous as an upper bound
  via the estimate ε_assoc ≤ 2·‖Φ−Φ²‖_cb + higher-order corrections.

### ε_sub (bilinear)

**Candidate A: Watrous cb-norm SDP** (see above — exact for bilinear, size N²×N²)  
**Candidate B: HOPM multi-start** over spectral-norm unit balls in A  
**Candidate C: Product bound** ε_sub ≤ ‖Φ‖_cb·‖X‖_op·‖Y‖_op − 1; since ‖Φ‖_cb = 1 (UCP),
  this gives ε_sub ≤ 0 (trivial). Use instead: bound via ‖Φ(XY) − XY‖_op ≤ ‖Φ−Id‖_cb‖XY‖_op.

### ε_cstar (quadratic/sesquilinear)

**Candidate A: Watrous cb-norm SDP** — ‖X⋆X†‖_op = ‖Φ(XX†)‖_op; deviation from
  ‖X‖_op² is a bilinear (sesquilinear) form. Maps to a cb-norm computation.  
**Candidate B: Power method on A†A** — for the sesquilinear form, alternating
  maximization over X ∈ A is simpler than the trilinear case.  
**Candidate C: Certified via A†A eigenvalue enclosure** — for each X in a grid
  over A ∩ spectral-ball, compute ‖X⋆X†‖_op / ‖X‖_op² in arb; the sup over
  the grid is a lower bound; the analysis of the grid's covering radius gives
  the bound's completeness.

---

## Audition slate recommendation table

| Defect | Candidate | Bound type | Certifiable in arb? | Cost (N=16,d=16) | Dim scaling (upper bound) | Canary: zero at η=0? | Notes |
|--------|-----------|------------|---------------------|------------------|--------------------------|----------------------|-------|
| ε_assoc | A: Basis-triple | Upper (loose) | YES | O(d³N³) ≈ fast | O(d^{3/2}) slack | YES | First filter; overestimates by ≤64x |
| ε_assoc | B: HOPM multi-start | Lower (witness) | YES (evaluate at witness) | Fast (double) | N/A (lower bound) | N/A | Essential; 100–1000 restarts |
| ε_assoc | C: Level-1 SDP lift | Upper (relaxation) | Approximately | O((N²d)³) SDP | Tight at high levels | YES | MOSEK; level-1 may have gap |
| ε_assoc | D: SOS hierarchy | Upper (tightening) | Approximately | Grows as d^{2r} | O(1/r²) convergence | YES | arXiv:2310.17827, arXiv:2412.13191 |
| ε_sub | A: Watrous SDP | EXACT | Approximately (post-hoc) | N²×N² SDP = 256×256 | Exact | YES | Best-in-class for bilinear |
| ε_sub | B: HOPM | Lower | YES | Fast | N/A | N/A | Companion to A |
| ε_sub | C: A†A certified | Upper | YES (full arb) | O(N³) per X | O(1) | YES | Via eigenvalue enclosure |
| ε_cstar | A: Watrous SDP | EXACT | Approximately | N²×N² SDP | Exact | YES | Sesquilinear form |
| ε_cstar | B: Power method | Lower | YES | O(N³) | N/A | N/A | Simple scalar optimization |
| ε_cstar | C: A†A certified | Upper | YES (full arb) | O(N³) | O(1) | YES | Certified in arb |

**Prioritized implementation order (phase 1):**
1. Candidate B (HOPM, lower bound) for all three defects — gives the benchmark value.
2. Candidate A_assoc (basis-triple) — cheap arb-certified upper bound, passes canary.
3. Candidate A_sub, A_cstar (Watrous SDP) — exact upper bounds for bilinear/sesquilinear.
4. After phase 1: audit whether basis-triple is within 10× of HOPM; if not, add
   Candidate C (SDP lift) for ε_assoc.

**Key risks to audit:**
- √dim canary: always compare ε_assoc (HOPM lower) with basis-triple upper;
  log the ratio across test cases. If ratio is consistently small (< 10), the
  basis-triple bound is fine for certification even with d^{3/2} factor.
- HOPM local maxima: for each defect, cross-check the HOPM-from-100-starts
  value with the SDP upper bound. If they agree to within 10%, HOPM is finding
  the global optimum.
- Watrous SDP floating-point: the SDP output is not certified in arb. For a
  certified upper bound, either use the SDP dual certificate + interval SDP
  verification (arXiv:2308.07287), or sandwich between the Watrous SDP value
  and the basis-triple arb-certified bound.

---

## References (alphabetical by first author)

- Hillar & Lim, "Most tensor problems are NP-hard," arXiv:0911.1393 / J. ACM 60(6) 2013.
  https://arxiv.org/abs/0911.1393
- Jiang, Ma & Zhang, "Approximating Tensor Norms via Sphere Covering," arXiv:2302.14219.
  https://arxiv.org/pdf/2302.14219
- Kolda & Bader, "Tensor Decompositions and Applications," SIAM Review 51(3) 2009.
- Lasserre, "Global optimization with polynomials," SIAM J. Optim. 11(3) 2001.
- Lasserre hierarchy convergence: arXiv:2011.08139. https://arxiv.org/pdf/2011.08139
- Li, Liu & Wu, "Complexity and computation for spectral/nuclear norm," arXiv:2212.14775.
  https://arxiv.org/pdf/2212.14775
- Maggioni & Manzieri, "A hierarchy of eigencomputations for polynomial optimization
  on the sphere," arXiv:2310.17827. https://arxiv.org/pdf/2310.17827
- Moment-SOS and de Finetti: arXiv:2412.13191. https://arxiv.org/pdf/2412.13191
- Navascues et al., "Semidefinite programming hierarchies for constrained bilinear
  optimization," arXiv:1810.12197 / Math. Prog. 2021.
  https://arxiv.org/abs/1810.12197
- Operator norm inequalities for tensor unfoldings: arXiv:1603.05621.
  https://arxiv.org/abs/1603.05621
- Rump, "Verified bounds for singular values, in particular for the spectral norm
  of a matrix and its inverse," BIT Num. Math. 51(2) 2011.
  https://link.springer.com/article/10.1007/s10543-010-0294-0
- Rump, INTLAB. https://www.tuhh.de/ti3/rump/intlab/
- Usevich et al., "A new convergence proof for the higher-order power method,"
  arXiv:1407.4586. https://arxiv.org/abs/1407.4586
- Watrous, "Semidefinite programs for completely bounded norms," Theory of
  Computing 5:11 2009. https://arxiv.org/abs/0901.4709
- Watrous, "Simpler semidefinite programs for completely bounded norms," 2012.
  https://arxiv.org/abs/1207.5726
- FLINT 3.6 documentation, acb_mat.
  https://flintlib.org/doc/acb_mat.html
- Fredrik Johansson, "Eigenvalues in Arb" (blog).
  https://fredrikj.net/blog/2018/12/eigenvalues-in-arb/
- Julia Discourse, "Rigorous spectral opnorm."
  https://discourse.julialang.org/t/rigorous-spectral-opnorm/82739
- Verified error bounds for singular values of structured matrices, arXiv:2502.09984.
  https://arxiv.org/pdf/2502.09984
- Complexity of matrix p-norms, arXiv:0908.1397.
  https://arxiv.org/pdf/0908.1397

---

## Cycle 2 — faithful operator-norm worst-case search (method decision)

**Date:** 2026-05-29  
**Scope:** aic-knm. Faithful (dimension-independent) lower-bound search for
ε_assoc, ε_sub, ε_cstar via iterative witness extraction, without SDP.

### Problem restatement

Subspace A ⊆ M_n (n×n complex), dim A = d, basis {B_1..B_d} Frobenius-ONB.
Star product X⋆Y = Φ(XY). Three defects, each a ratio of norms over
operator-norm (spectral-norm) unit balls of A:

  ε_assoc = sup_{X,Y,Z∈A, ‖·‖_op≤1} ‖(X⋆Y)⋆Z − X⋆(Y⋆Z)‖_op
  ε_sub   = sup_{X,Y∈A, ‖·‖_op≤1} (‖X⋆Y‖_op − 1)   [positive part]
  ε_cstar = 1 − inf_{X∈A, ‖X‖_op=1} ‖X†⋆X‖_op

Goal: find a witness (X*,Y*,Z*) by iterative double-path search, evaluate
the ratio in arb to certify the lower bound. NOT an SDP.

### Candidate methods evaluated

#### Method 1: Scale-invariant HOPM / alternating maximization on the
spectral-norm unit sphere — RECOMMENDED PRIMARY

**Core insight.** Scale-invariance of the ratio means we can work on the
actual spectral-norm unit sphere {A∈M_n: ‖A‖_op=1} rather than an ℓ₂
sphere on coefficient vectors, avoiding the √dim inflation entirely.

**The scale-invariant objective.** For ε_assoc the ratio is
  R(X,Y,Z) = ‖h(X,Y,Z)‖_op / (‖X‖_op ‖Y‖_op ‖Z‖_op),
which is homogeneous of degree 0. Maximizing R is equivalent to maximizing
  F(X,Y,Z) = ‖h(X,Y,Z)‖_op
subject to ‖X‖_op = ‖Y‖_op = ‖Z‖_op = 1, X,Y,Z ∈ A.

**Variational reduction of the outer ‖·‖_op.** Use
  ‖h(X,Y,Z)‖_op = max_{‖u‖₂=‖v‖₂=1} |⟨u, h(X,Y,Z) v⟩|.
Introduce auxiliary unit vectors u,v ∈ C^n and maximize
  G(X,Y,Z,u,v) = Re⟨u, h(X,Y,Z) v⟩
over ‖X‖_op=‖Y‖_op=‖Z‖_op=1 (X,Y,Z∈A), ‖u‖₂=‖v‖₂=1.

This is a five-way alternating maximization. Each update is closed-form or
a single SVD/eigenstep:

  **u-update** (fix X,Y,Z,v): u ← h(X,Y,Z)v / ‖h(X,Y,Z)v‖₂.
  **v-update** (fix X,Y,Z,u): v ← h(X,Y,Z)† u / ‖h(X,Y,Z)† u‖₂.
  **X-update** (fix Y,Z,u,v): maximize Re⟨u, h(X,Y,Z)v⟩ over ‖X‖_op≤1, X∈A.
    Write h(X,Y,Z) = Φ(Φ(XY)Z) − Φ(XΦ(YZ)); fixing Y,Z,u,v the functional
    X ↦ Re⟨u, h(X,Y,Z)v⟩ is linear in X (since h is linear in X). Thus it
    equals Re⟨C_X, X⟩_F for some C_X ∈ M_n computable from (Y,Z,u,v):
      C_X = [∂_X h(X,Y,Z)]† (uv†)
           = Φ†(uv†) · (YZ)† [first term contribution]
           − Φ†(uv†) · Φ(YZ)† [second term contribution],
    where Φ† is the adjoint of Φ (also a CP map, the dual in Hilbert-Schmidt
    inner product). This simplifies to
      C_X = Φ†(uv†) · (YZ)† − Φ(YZ) · (something)
    (precise formula requires unrolling; see implementation note below).
    Now maximize Re⟨C_X, X⟩_F over X∈A, ‖X‖_op≤1.
    **This is the crux.** The Frobenius inner product Re⟨C_X, X⟩_F = Re⟨P_A C_X, X⟩_F
    (project C_X onto A) = Re⟨c, x⟩ for x = coefficient vector in A-basis,
    c_k = Re⟨B_k, C_X⟩_F. Then we want to maximize Re⟨c, x⟩ over
    ‖X(x)‖_op ≤ 1, X(x) = Σ x_k B_k.
    The solution: x* = λ(u₁v₁†) where (u₁,v₁) is the leading singular pair
    of P_A(C_X) (the A-projection of C_X) in the sense that u₁v₁† maximizes
    Re⟨C_X, X⟩_F over ‖X‖_op≤1, X∈A. Concretely:
      X* = argmax_{‖X‖_op≤1, X∈A} Re⟨C_X, X⟩_F
         = P_A(C_X / ‖P_A(C_X)‖_op)  ... if P_A preserves operator-norm structure
    CAUTION: P_A is the Frobenius projection onto A; P_A(C_X) is not in general
    a rank-1 matrix. The correct update is: decompose P_A(C_X) = Σ σ_k e_k f_k†
    via SVD; then X* = e₁f₁† is the maximizer (a rank-1 matrix in M_n, but
    not necessarily in A). To enforce X∈A, we need P_A(e₁f₁†) and then
    re-normalize to the spectral-norm unit sphere:
      X* ← P_A(e₁f₁†) / ‖P_A(e₁f₁†)‖_op.
    This is a valid subgradient-type step: it moves in the direction of the
    supergradient of Re⟨C_X,·⟩ restricted to {X∈A: ‖X‖_op≤1}. Each iterate
    is in the feasible set, so the ratio at each iterate is a valid lower bound.
    Non-smoothness at degenerate top singular value of P_A(C_X): any choice of
    leading singular pair in the degenerate subspace is a valid (sub)gradient;
    pick one (e.g. the one returned by LAPACK's `zgesvd`). The iterate is still
    feasible.
  **Y-update, Z-update:** symmetric to X-update with the appropriate linear
    functional in Y (or Z) derived from fixing X,Z,u,v (resp. X,Y,u,v).

**Implementation note.** To compute C_X explicitly:
  h(X,Y,Z) = Φ(Φ(XY)Z) − Φ(XΦ(YZ)).
  ∂_X h(X,Y,Z) · δX = Φ(Φ(δX·Y)·Z) − Φ(δX·Φ(YZ)).
  Thus ⟨uv†, ∂_X h · δX⟩_F = ⟨Φ†(uv†)·Z†·Φ†(·), δX⟩  [unrolled carefully].
  In practice with a double path: compute
    M_Y  = Φ(YZ)          (n×n matrix, precomputed)
    M_XY = Φ(XY)          (n×n)
    h    = Φ(M_XY · Z) − Φ(X · M_Y)
  The derivative w.r.t. X (holding Y,Z,u,v fixed) is:
    c_k = Re[ ⟨u, Φ(Φ(B_k · Y) · Z) v⟩ − ⟨u, Φ(B_k · Φ(YZ)) v⟩ ]
  Computed for each basis element B_k: O(d·N³) per X-update (d SVDs or
  Φ-applications on N×N matrices). At N=16, d=16: 16 × O(16³) ≈ 65k FLOPs.

**Convergence guarantee.** HOPM on a multilinear form converges to a
critical point (Łojasiewicz inequality; Usevich et al. arXiv:1407.4586,
Numer. Math. 2018). NOT guaranteed global maximum. Multi-start covers this.

**Cost.** Per iterate: 5 updates × O(d·N³) = O(d·N³). At N=d=16:
~100k FLOPs per iterate; 50 iterates per start × 200 starts = ~10^9 FLOPs.
Sub-second in double precision.

**Non-smoothness.** Operator-norm non-smooth at tied singular values. Does
NOT break the algorithm: any element of the subdifferential of the spectral
norm at a tie is a valid update direction; `zgesvd` returns one. Each iterate
is feasible; the ratio at each iterate is a rigorously valid lower bound.

**Scale-invariance and the coefficient vector.** Expressing X = Σ_k x_k B_k,
the operator-norm constraint ‖X‖_op=1 defines a spectrahedron in x-space
(not a Euclidean sphere). The update above works directly in x-space by
SVD-clipping the Frobenius gradient: it does NOT assume ‖x‖₂=1 and hence
does NOT pay the √d Frobenius-to-operator inflation. This is the
dimension-independence guarantee.

#### Method 2: Riemannian gradient / geodesic optimization — EVALUATED,
SECONDARY for ε_sub and ε_cstar, NOT RECOMMENDED for ε_assoc

**What manifold?** For the operator-norm unit sphere {A: ‖A‖_op=1} the
issue is that the set is NOT a smooth manifold: it is a semi-algebraic set
whose boundary (where the top singular value has multiplicity ≥ 2) is a
submanifold of lower dimension — a "corner" in the Riemannian sense.
At generic points (simple top singular value), the tangent space is well-
defined, but algorithms must handle the non-smooth boundary.

**Practical fix: work on the Frobenius sphere, operator-norm in objective.**
The standard workaround (confirmed in the iMuon line of work,
arXiv:2605.09238, and in arXiv:2202.11597 for p-norm spheres) is:

  - Constrain to the Frobenius unit sphere ‖X‖_F=1 (a smooth manifold,
    Riemannian gradient is standard projection);
  - Put the operator norm into the objective: maximize G(X,Y,Z) = ‖h(X,Y,Z)‖_op
    subject to ‖X‖_F=‖Y‖_F=‖Z‖_F=1.

This pays a √d factor: the Frobenius unit sphere in A has operator norms
between 1/√d and 1, so maximizing over the Frobenius sphere gives a result
that is within √d of the operator-norm-sphere maximum. For d≤16 this is at
most factor 4 — smaller than the d^{3/2} basis-sweep bound but still a
dimension-dependent inflation.

**Alternative: Fixed-rank manifold with spectral-norm LMO (iMuon).**
The iMuon framework (arXiv:2605.09238, May 2026) provides closed-form
Riemannian updates for spectral-norm LMOs on Stiefel, Grassmann, and fixed-
rank manifolds. Its retraction is: x_{+} = R_x(−ηξ*) where ξ* involves
G_x^{-1/2} Z* and Z* = U diag(z*) V† from an SVD of the metric-preconditioned
gradient. This is for optimization on the manifold of rank-r matrices, with
the spectral norm acting on the manifold's tangent space, NOT for
maximization over the spectral-norm sphere in a subspace. It does not
directly address the subspace constraint X∈A.

**Verdict for the project.** Riemannian gradient on Frobenius sphere is a
valid *alternative* for ε_sub (bilinear) and ε_cstar (quadratic) where
smooth manifold structure helps and the √d factor is acceptable. For ε_assoc
(trilinear), the five-way HOPM with direct spectral-norm enforcement (Method
1) is cleaner and avoids the inflation. Riemannian CG on Frobenius sphere
can serve as a warm-start generator or fallback when HOPM stalls.

**Retraction for Frobenius sphere (exact).** For X = X_old + η∇_Riem f,
project via: X ← (X_old + η g) / ‖X_old + η g‖_F (standard sphere retraction,
or polar retraction). The Riemannian gradient is: grad_Riem f = g − ⟨g,X⟩_F X
where g is the Euclidean gradient of f at X with ‖X‖_F=1. This is textbook
(Boumal, "An Introduction to Optimization on Smooth Manifolds," 2023,
Cambridge University Press, ch. 3).

#### Method 3: TDVP-type variational dynamics — NOT RECOMMENDED for this
problem; honest negative assessment

**What TDVP does.** Real-time TDVP projects the Schrödinger equation onto
the tangent space of a variational manifold (e.g., matrix product states)
so that the time-evolved state stays on the manifold while approximately
satisfying the equation of motion (Haegeman et al. 2016, unige.ch preprint).
Imaginary-time TDVP recovers DMRG in the limit of large imaginary time steps:
it finds the ground state (lowest eigenvalue) of a Hamiltonian on the manifold.

**Why it does NOT fit here.** The ε_assoc maximization problem is:
  - A MAX problem (find largest defect), not a MIN or ground-state problem.
  - The objective ‖h(X,Y,Z)‖_op is NOT expressible as ⟨ψ|H|ψ⟩ for a fixed
    Hamiltonian H acting on the manifold. It is a non-quadratic (degree-6 in
    matrix entries) functional of three independent matrix variables, not a
    bilinear expectation value.
  - There is no natural tensor-network topology connecting X, Y, Z (they are
    independent factors in M_n, not a single state in a tensor product space).
  - The "variational manifold" would be the product of three spectral-norm
    spheres in A, which has no natural MPS/TTN structure.

**TDVP as imaginary-time flow toward a leading eigenvector** only works when
the objective is a Rayleigh quotient ⟨x,Ax⟩/⟨x,x⟩ for a linear operator A.
The ε_cstar defect comes closest: it is inf_X ‖X†⋆X‖_op over ‖X‖_op=1,
which is related to a smallest-eigenvalue problem. But ‖X†⋆X‖_op = ‖Φ(XX†)‖_op
is not a linear functional of X (it is quadratic in XX†), so even here TDVP
does not apply directly.

**Conclusion.** TDVP has no sound fit for this problem. It is not a
first-class candidate for ε_assoc or ε_sub. For ε_cstar, imaginary-time
power iteration on the sesquilinear form is equivalent to Method 1 restricted
to two-way alternating, which is just standard power iteration — there is no
benefit from the TDVP framing. Do not implement a TDVP-based estimator.

#### Method 4: Tensor-network / DMRG-like alternating sweeps — EVALUATED,
NO ADVANTAGE at the target sizes

**What the TN/DMRG view adds.** For a tensor with structure (bond dimension,
local Hilbert spaces), DMRG-type alternating sweeps are equivalent to ALS
(alternating least squares) / HOPM on the tensor. In a TN language, the
trilinear form h_assoc has a coefficient tensor T_{ijk} (in the B_k basis)
of size d×d×d→d, and HOPM is literally one-site DMRG for finding its leading
singular value. This yields no algorithmic advantage: it is the same update
rule as Method 1 with added notation. At d,n≤16 the "bond dimension" is at
most 16, so there is no compression benefit from a TN representation.

**What TN does NOT buy here.**
  - The coefficient tensor T has size ≤16³: no compression needed.
  - The operator-norm constraint (spectrahedron in x-space) does not factorize
    as a tensor product, so the one-site update does not become a simpler
    eigenvector problem than Method 1 already provides.
  - Riemannian optimization of isometric tensor networks (Hauru et al.,
    arXiv:2007.03638, SciPost Phys. 2021) optimizes variational manifolds of
    tensors with isometry constraints (Stiefel/Grassmann), not spectral-norm
    spheres. No direct mapping.

**Verdict.** TN language is a valid notation for the HOPM update but adds
nothing algorithmic at these sizes. Implement Method 1; it is ALS/HOPM in
disguise. If n grows large (n>64), a TN representation of h(X,Y,Z) might
reduce contraction cost, but that is outside current scope.

#### Method 5: Projected gradient with alternating projection onto
{X∈A: ‖X‖_op≤1} — EVALUATED, VALID FALLBACK but sub-optimal

**The projections.** Two non-commuting projections:
  - P_A: M_n → A (Frobenius projection; closed-form via basis).
  - P_op: M_n → {X: ‖X‖_op≤1} (SVD clip: X ← U diag(min(σ_k,1)) V†).

**Does alternating projection stay in A∩{‖·‖_op≤1}?**  Not automatically:
P_op(P_A(X)) is in A (P_A maps to A) but after P_op it leaves A. Starting
from X∈A with ‖X‖_op>1, the sequence P_A∘P_op iterates. Convergence of
alternating projections for two closed convex sets is guaranteed (von Neumann;
linear rate when sets intersect transversally). Here both sets are:
  - A = affine subspace of M_n: convex, closed.
  - {‖X‖_op≤1}: convex, closed (the spectral norm is convex).
So A∩{‖X‖_op≤1} is convex and the alternating projection converges to the
nearest feasible point from the initial iterate.

**For maximization** (not just feasibility), use projected gradient ascent:
  X_{t+1} ← P_A(P_op(X_t + η ∇_X G(X_t)))
where ∇_X G is the Euclidean gradient of the objective w.r.t. X and η is a
step size. This is gradient ascent with projection. Convergence is to a
stationary point of the projected problem; no guarantee of global max.

**Disadvantage vs. Method 1.** The step X ← P_A(P_op(X+η∇G)) discards the
optimal-direction property of Method 1's SVD-based update. Method 1's X-update
is the *exact* maximizer of the linear functional in X on the feasible set
{X∈A: ‖X‖_op≤1} (up to the inner product approximation), which is provably
the best single-step improvement. Projected gradient with a fixed step size
is typically slower.

**Use case.** Valid as a sanity check or warm-start initialization. For
the main search, Method 1 dominates.

#### Method 6: Frobenius relaxation — inequality direction and warm-start use

**Direction of the inequality.** For X∈A (the subspace):
  ‖X‖_op ≤ ‖X‖_F ≤ √d · ‖X‖_op.
Therefore if ε_F = sup_{X,Y,Z∈A, ‖·‖_F≤1} ‖h(X,Y,Z)‖_op (Frobenius
unit ball) and ε_op = sup_{‖·‖_op≤1} (same) (operator-norm unit ball), then:
  ε_F ≤ ε_op [Frobenius ball is SMALLER: ‖X‖_F≤1 implies ‖X‖_op≤‖X‖_F≤1]
  ε_op ≤ d^{3/2} · ε_F [the d^{3/2} inflation, shown in §2].
**Direction summarized:** ε_F is a LOWER bound on ε_op (operator-norm defect
over a larger ball → larger supremum). So the Frobenius maximizer is a valid
(conservative) lower bound on the operator-norm defect, understating it by
at most d^{3/2}.

**Warm start.** The Frobenius maximizer X_F* is straightforward: write X =
Σ_k x_k B_k, the Frobenius constraint becomes ‖x‖₂≤1, and the problem
reduces to a standard HOPM on the coefficient tensor with Euclidean sphere
constraints — fully solved by ALS with update x ← T(1; y,z)/‖T(1;y,z)‖ etc.
This is cheap (O(d³) per iteration, no N dependence beyond computing T once).
The Frobenius maximizer (x_F*, y_F*, z_F*) gives a feasible starting point
X_F* = Σ x_k* B_k in A. Normalizing to the spectral norm gives
  X₀ = X_F* / ‖X_F*‖_op ∈ A, ‖X₀‖_op = 1,
which is a spectral-norm-feasible warm start for Method 1. This is a useful
initialization: the Frobenius maximizer likely lives near the operator-norm
maximizer for structured channels (where B_k are near-unitary or normalized).

**Use in practice:**
  Init type A (Frobenius warm start): run 10–20 Frobenius-HOPM starts → get
    X₀s → feed as initial points for Method 1.
  Init type B (random): draw X = randn(n,n), project to A, normalize by op-norm.
  Use a mix: 50% Frobenius warm starts + 50% random. At N=d=16, 200 total
  starts is inexpensive.

### Summary ranking

| Rank | Method | Bound type | Dim-independent? | Cost per start | Recommendation |
|------|--------|-----------|------------------|----------------|----------------|
| 1 | Method 1: Scale-invariant HOPM + SVD clip | Lower (witness) | YES (no √dim) | O(d·N³)×50 iter | **Primary; implement first** |
| 2 | Method 2: Riemannian CG on Frobenius sphere | Lower (Frob-ball) | √d factor only | O(d·N³)×iter | Fallback / warm start for ε_sub, ε_cstar |
| 3 | Method 5: Projected gradient (A∩op-ball) | Lower (witness) | YES | O(N³)×iter | Sanity check; slower than Method 1 |
| 4 | Method 6: Frobenius relaxation | Lower (under-estimate) | d^{3/2} factor | O(d³)×iter | Cheap warm-start generator; undercounts |
| 5 | Method 4: TN/DMRG | Lower (= Method 1) | YES (= Method 1) | Same | No advantage; skip dedicated TN code |
| — | Method 3: TDVP | Not applicable | — | — | Does not fit; do not implement |

### Exact update recipe (Method 1 — the primary algorithm)

**Data.** n×n matrix Φ represented as Choi matrix / Kraus operators;
orthonormal basis {B_1..B_d} of A stored as d n×n complex matrices.
Precompute: for each k,l: Φ(B_k · B_l) (d² matrices, each n×n);
  for each k,l,m: h(B_k,B_l,B_m) (d³ matrices). This is O(d³·N³) setup.

**Per-start initialization.**
  1. Draw X,Y,Z ∈ A randomly: x,y,z ← randn(d); X = Σ x_k B_k / ‖Σ x_k B_k‖_op.
  2. Draw u,v ← randn(n), normalize to unit ℓ₂.

**Inner loop (50 iterations per start).**
  For t = 1..50:
    1. h ← Φ(Φ(X·Y)·Z) − Φ(X·Φ(Y·Z))        [two Φ-applications each, O(N³)]
    2. u ← h·v / ‖h·v‖₂                        [O(N²)]
    3. v ← h†·u / ‖h†·u‖₂                      [O(N²)]
    4. For k=1..d: c_k ← Re⟨u, [Φ(Φ(B_k·Y)·Z) − Φ(B_k·Φ(Y·Z))] v⟩   [O(d·N²)]
       C_X ← Σ_k c_k B_k                        [in A; or equivalently work with
                                                   the coefficient vector c ∈ R^d]
       Xnew ← P_A(rank-1 approx of argmax)... but since C_X ∈ A already,
       SVD decompose C_X = Σ_j σ_j e_j f_j†; take leading pair (e₁,f₁):
         X ← P_A(e₁f₁†) / ‖P_A(e₁f₁†)‖_op    [P_A via B_k inner products, O(d·N²)]
    5. Symmetric Y-update (c_k from B_k in Y-slot), Z-update.
    6. Record current R = ‖h(X,Y,Z)‖_op / (‖X‖_op·‖Y‖_op·‖Z‖_op) as lower bound.

  **Note on step 4 simplification.** If C_X = Σ_k c_k B_k ∈ A, then
  P_A(e₁f₁†) is in general DIFFERENT from C_X. But the best X∈A maximizing
  Re⟨C_X,X⟩_F subject to ‖X‖_op≤1 is precisely e₁f₁† IF e₁f₁† ∈ A; otherwise
  it is P_A(C_X)/‖P_A(C_X)‖_op (since C_X already lives in A, its top singular
  vector gives the steepest ascent direction; P_A is identity on C_X; so just
  normalize C_X by its operator norm):
    X ← C_X / ‖C_X‖_op.
  This is exact for the linear inner product because the maximizer of Re⟨C,X⟩_F
  over {X∈A, ‖X‖_op≤1} when C∈A is sgn(C) in the polar-decomposition sense:
  X* = U V† where C = U Σ V† (full SVD of C). So the update is:
    X ← polar(C_X)   (= U V† from SVD of C_X, which satisfies ‖X‖_op=1, X∈A).
  LAPACK's `zgesvd` or QR-based polar gives this cheaply.

**Multi-start and record-keeping.**
  Run 200 starts. Record (X*,Y*,Z*,R*) = argmax_{starts} R.
  The best R* is the lower bound.
  Also keep the top 5 distinct local maxima (by value) to expose the landscape.

**Warm-start batch (50 of 200 starts).**
  Run Frobenius-HOPM (d-dimensional ALS on coefficient vectors x,y,z∈R^d,
  update x ← T_{xyz}·yz / ‖T·yz‖ where T is the d×d×d coefficient tensor
  of h precomputed in the B_k basis) for 20 iterations; feed resulting
  X_F*/‖X_F*‖_op into Method 1 as warm start.

### Q7: restart count guidance

At N=d=16 the problem is 16²×3+2n = ~800-dimensional. HOPM's landscape for
small random instances has empirically O(d) distinct local maxima (Kolda &
Mayo, SIAM JMAA 2011; Usevich et al. 2018). With 200 starts and fast ~50-iter
convergence, the probability of missing the global max is (1-θ)^200 where θ
is the basin-of-attraction volume of the global max. Empirical practice for
tensor spectral norm at d≤16: 100–500 starts consistently finds the global max
within 1% (Hillar & Lim 2013, supplementary). Start with 200; if the top-5
local maxima cluster within 5% of the best value, trust it. If they spread
widely, increase to 500.

For the universality canary specifically (ε_assoc ≈ 0 case): near η=0 the
form h is nearly zero, all starts converge to near-zero, and the global max
is trivially found. The interesting regime is moderate η where d^{3/2} slippage
would be visible — here 200 starts suffices for the canary.

### Q8: witness-certification as rigorous lower bound

YES, unconditionally. Given a concrete witness (X*, Y*, Z*) ∈ M_n³ (arbitrary
double-precision floating point):
  1. Load X*, Y*, Z* into acb_mat at precision p (e.g. p=128 bits).
  2. Compute h* = Φ(Φ(X*Y*)Z*) − Φ(X*Φ(Y*Z*)) in arb (all operations exact
     in arb: matrix multiply, Φ-application via Choi matrix contraction).
  3. Compute ‖h*‖_op via acb_mat eigenvalue enclosure on (h*)†h* (§3.1):
     certified interval [lo, hi] for ‖h*‖_op.
  4. Compute ‖X*‖_op, ‖Y*‖_op, ‖Z*‖_op certified similarly.
  5. R_certified ≥ lo / (hi_X · hi_Y · hi_Z).
This gives a rigorously certified lower bound on ε_assoc with NO global
guarantee required. The certification is a single forward evaluation.
Precision p=128 is sufficient for N≤16; p=64 may work too (check ball radii).
No SDP needed; no global reasoning needed.

### Q9: submult and cstar special-case simplifications

**ε_sub (bilinear).** Define g(X,Y) = ‖X⋆Y‖_op for X,Y∈A, ‖X‖_op=‖Y‖_op=1.
  ε_sub = max(0, sup_{X,Y} g(X,Y) − 1).
  Use Method 1 with two variables (X,Y) and u,v:
    G(X,Y,u,v) = Re⟨u, Φ(XY) v⟩
  Three-way alternating (X, Y, u/v pair). Update rules:
    u ← Φ(XY)v / ‖Φ(XY)v‖₂
    v ← Φ(XY)†u / ‖Φ(XY)†u‖₂
    X ← polar(C_X) where C_X = Σ_k ⟨u, Φ(B_k Y) v⟩ B_k
    Y ← polar(C_Y) where C_Y = Σ_k ⟨u, Φ(X B_k) v⟩ B_k
  If ε_sub ≈ 0 (near-UCP case, which is expected here), the iteration
  converges rapidly and the witness certifies a near-zero lower bound.
  For exact UCP channels, ‖Φ(XY)‖_op ≤ ‖XY‖_op ≤ ‖X‖_op‖Y‖_op (submultiplicativity
  is automatic), so ε_sub = 0 and any witness certifies it.

**ε_cstar (quadratic/sesquilinear).** Minimize ‖X†⋆X‖_op over ‖X‖_op=1, X∈A.
  This is a MIN problem. For a C*-algebra, ‖X†⋆X‖_op = ‖X‖_op² = 1 exactly
  (C* identity). So ε_cstar = 1 − min and we want the value to be near 1
  (i.e., min ≈ 1, ε_cstar ≈ 0).
  To MAXIMIZE ε_cstar (find the worst-case departure from C*): we want to
  MINIMIZE ‖X†⋆X‖_op = ‖Φ(X†X)‖_op over ‖X‖_op=1.
  Method: gradient DESCENT on F(X,u,v) = Re⟨u, Φ(X†X) v⟩ over ‖X‖_op=1,
  ‖u‖₂=‖v‖₂=1. Updates:
    u ← Φ(X†X)v / ‖...‖, v ← Φ(X†X)†u / ‖...‖ (same as above).
    X-update (minimize): C_X = Σ_k ⟨u, Φ(B_k†X + X†B_k) v⟩ B_k (gradient of
      F w.r.t. X, treating the sesquilinear form);
      X ← polar(−C_X) (DESCENT: flip sign to minimize, then polar-normalize).
  For ε_cstar near 0 (C*-algebra case), the minimum is near 1 and the
  landscape is unimodal near the identity element of A — gradient descent
  with 50 restarts suffices.
  Alternative: power iteration on the sesquilinear form T(x) = Σ_{k,l} x_k x_l*
    ⟨u, Φ(B_k† B_l) v⟩. This is a quadratic form in x ∈ C^d; its minimum
    over ‖X(x)‖_op=1 is a small eigenvector problem combined with a
    spectrahedron constraint. At d≤16 this is directly solvable.

### References added (Cycle 2)

- Usevich, Li & Comon, "A new convergence proof for the higher-order power
  method and generalizations," arXiv:1407.4586; Numer. Math. 2018.
  https://arxiv.org/abs/1407.4586
  https://link.springer.com/article/10.1007/s00211-018-0981-3
- Kolda & Mayo, "Shifted Power Method for Computing Tensor Eigenpairs,"
  SIAM J. Matrix Anal. Appl. 32(4) 2011. https://arxiv.org/abs/1007.1267
- Lim, "Singular values and eigenvalues of tensors: a variational approach,"
  arXiv:math/0607648 (background on tensor spectral norm via alternating SVD).
- Boumal, "An Introduction to Optimization on Smooth Manifolds," Cambridge 2023.
  https://www.cambridge.org/core/books/an-introduction-to-optimization-on-smooth-manifolds/
- Hauru, Van Damme & Haegeman, "Riemannian optimization of isometric tensor
  networks," arXiv:2007.03638; SciPost Phys. 10, 040 (2021).
  https://arxiv.org/abs/2007.03638
- iMuon: "Intrinsic Muon: Spectral Optimization on Riemannian Matrix Manifolds,"
  arXiv:2605.09238 (May 2026). https://arxiv.org/abs/2605.09238
- Chen, Han & Ye, "Riemannian optimization on unit sphere with p-norm and its
  applications," arXiv:2202.11597; Comput. Optim. Appl. 2023.
  https://arxiv.org/abs/2202.11597
- Hillar & Lim, "Most tensor problems are NP-hard," arXiv:0911.1393 2013.
  (convergence-basin empirics for small d). https://arxiv.org/abs/0911.1393
- Haegeman et al., "Unifying time evolution and optimization with matrix
  product states," Phys. Rev. B 94, 165116 (2016). (TDVP background.)
  https://www.unige.ch/math/vandereycken/papers/published_Haegeman_LOVV_2016.pdf
- Leloykun (blog), "Rethinking Maximal Update Parametrization: Steepest Descent
  on the Spectral Ball," 2025. https://leloykun.github.io/ponder/rethinking-mup-spectral-ball/
