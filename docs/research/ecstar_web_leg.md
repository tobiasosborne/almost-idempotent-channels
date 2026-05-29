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
