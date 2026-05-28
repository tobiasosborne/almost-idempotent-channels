# Shard G — UCP map details & idempotent-structure proof (L1562–2161)

## Scope

§11 "Details on unital completely positive maps" (sec_channels, L1562–2161).

- §11.1 (L1566–1687): Choi representation (finite-dim), Stinespring representation (general), canonical Stinespring, prop_KLHG.
- §11.2 (L1689–1742): Schwarz inequality for UCP maps, cb-norm bound, carrier definition, lem_carrier, cor_carrier, equation PhiX_M.
- §11.3 (L1744–2161): idempotent UCP map structure: lem_idemp, deferred proof of th_idemp_structure (L318/L2055), prop_Gamma, dual Enc/Dec forms.

---

## Core data structures (how to represent a UCP map in code)

### Choi representation (L1568–1617, eq. Choi, L1569)

For finite-dimensional H, K, any UCP map Phi: B(K) -> B(H) is represented by an isometry V: H -> K ⊗ F:

```
Phi(X) = V† (X ⊗ 1_F) V,   V†V = 1_H,   V: H -> K ⊗ F.
```

The paper calls (F, V) the "Choi representation". The auxiliary space F is constructed canonically as F = K* ⊗ H (L1575). In matrix terms, if dim K = k and dim H = h, then dim F ≤ kh, and V is a kh × h isometry matrix (column-isometry). The paper notes this is equivalent to Choi's original matrix construction [Cho75] but written abstractly.

**Kraus form** (not named "Kraus" in the paper): choosing an ONB {f_a} of F and writing V as block column V = col(K_a) where K_a: H -> K, the Choi representation becomes:
```
Phi(X) = Σ_a  K_a† X K_a,   Σ_a K_a K_a† = 1_H   (unitality).
```
This follows directly from expanding V†(X ⊗ 1_F)V in the F basis; it is standard but the paper does not separately state it. The K_a are the "Kraus operators" in quantum information language.

**Choi matrix** (classical Choi [Cho75], implicit at L1575): fix an ONB {e_i} of K; define
```
C_Phi = Σ_{i,j} E_{ij} ⊗ Phi(E_{ij})   ∈ B(K) ⊗ B(H),
```
where E_{ij} = |e_i><e_j|. Then C_Phi ≥ 0 iff Phi is CP; Tr_K(C_Phi) = dim(K) * 1_H iff Phi is unital (after normalization). The paper uses V directly rather than C_Phi as its primary object; the two are related by C_Phi = (1_K ⊗ V)(Omega)(1_K ⊗ V†) where Omega is the maximally entangled state on K ⊗ K*.

**Canonical construction** (L1575): V can always be taken with F = K* ⊗ H, dim F = k*h. Concretely, index F by pairs (i, alpha) with i in {1,...,k} and alpha in {1,...,h}; then V|psi> = Σ_i |e_i> ⊗ |phi_{i,psi}> where {phi_{i,alpha}} encodes Phi. This is not the minimal Kraus rank; minimizing rank(F) is equivalent to finding a minimal Kraus decomposition (by eigendecomposing C_Phi and keeping nonzero eigenvectors).

### Stinespring representation (L1621–1633, eq. Stinespring_def, L1622)

For a general C* algebra A (possibly infinite-dim):
```
Phi(X) = V† u(X) V,
```
where u: A -> B(G) is a *-homomorphism and V: H -> G is an isometry. The canonical Stinespring uses G = completion of (A ⊗ H)/null with inner product <B⊗eta | A⊗xi> = eta† Phi(B†A) xi (L1626–1628), and u(C)(A⊗xi) = CA⊗xi, V xi = I⊗xi (L1632).

In finite dimensions A = B(K), the canonical Stinespring coincides (up to the choice of F) with the Choi representation at L1568.

### cb-norm via Choi representation (L1717–1719)

For any UCP map Phi, ||Phi||_cb ≤ 1. Proof: the Schwarz inequality Phi(X†)Phi(X) ≤ Phi(X†X) (L1692–1714) implies ||Phi(X)|| ≤ ||X||; replacing Phi by 1_Mn ⊗ Phi gives the cb-norm bound. The cb-norm equals 1 for any UCP map (not merely ≤ 1; equality is achieved by X = 1). Computing ||Phi||_cb in finite dimensions is a semidefinite program [Watrous, Chapter 3] cited at L352.

### Conversions between representations

| From | To | Method |
|------|----|--------|
| Choi matrix C_Phi (kh × kh PSD matrix) | Kraus operators K_a | Eigendecompose C_Phi; nonzero eigenvectors (reshaped) give K_a |
| Kraus operators K_a | Isometry V | V = col(K_a): H -> K ⊗ F, F = C^r, r = number of Kraus ops |
| Isometry V | Choi matrix C_Phi | C_Phi = (1_K ⊗ V)(max-entangled state)(1_K ⊗ V†) |
| Stinespring (G, u, V) with A = B(K) | Isometry V' for Choi rep | u decomposes by Schur; compress to minimal invariant subspace |

**Checking UCP axioms from Choi matrix:**
- CP: C_Phi ≥ 0 (PSD).
- Unital: Tr_K(C_Phi) = 1_H (partial trace over the K factor = identity on H).

---

## The idempotent-structure algorithm (th_idemp_structure, proof L2055)

**Input**: idempotent UCP map Phi: B(H) -> B(H) in finite dim, given as isometry V: H -> H ⊗ F (Choi rep).

**Statement** (th_idemp_structure, L318–332): produce subspace M ⊆ H, C* algebra A, UCP maps Delta: A -> B(H) and Gamma: B(M) -> A such that:
```
Gamma C_M Delta = 1_A,    Delta Gamma C_M = Phi,
```
where C_M: B(H) -> B(M) is the compression X |-> J_M† X J_M. Furthermore, w = C_M Delta: A -> B(M) is a *-homomorphism, and Im(Delta) is block-diagonal w.r.t. H = M ⊕ M⊥.

### Algorithm steps

**Step 1 — Compute the carrier M** (lem_carrier, L1724).

M is the support of Phi*(rho_0) for any full-rank rho_0 on H. Equivalently, M is the smallest subspace such that M ⊗ F ⊇ Im(V).

*Matrix algorithm*: form the matrix V (kh × h after choosing dim H = h = dim K for the self-map case). Compute Im(V); then M = range of the projection Pi_M where Pi_M ⊗ 1_F is the projection onto Im(V) in H ⊗ F restricted to the H factor. Concretely: form the h × h matrix Q = (partial trace over F) of V V†, i.e., Q_{alpha,beta} = Σ_a V_{a,alpha} V*_{a,beta} for a ranging over the F basis. Then M = Im(Q) = column space of Q. Since Phi is unital, V†V = 1_H, so Tr(Q) = dim(F)*? — actually Q = Σ_a K_a K_a† = 1_H by unitality; thus M = all of H only if Im(V) already spans all of H ⊗ F, but M is generally a proper subspace.

Correction: lem_carrier says M is the smallest M' ⊆ K (here K = H) such that M' ⊗ F ⊇ Im(V). Compute Im(V) ⊆ H ⊗ F; then M = image under Tr_F (projection onto H) of Im(V), i.e., M = span{(Pi_H ⊗ bra_f) v : v ∈ Im(V), f ∈ F}. In matrix terms: form the basis matrix B of Im(V) (a subspace of C^{dim(H)*dim(F)}); reshape each basis vector as a dim(H) × dim(F) matrix; M = column span of all these matrices (union over basis vectors).

**Step 2 — Identify A as the image of Phi** (L2056).

A = Im(Phi) = Ker(1 - Phi) as a vector subspace of B(H). Basis: solve (1 - Phi)(X) = 0 for X ∈ B(H). Since Phi is a linear map on the dim(H)^2 dimensional space of matrices, form the matrix of the superoperator Phi (in terms of vec-ing), compute eigenspace for eigenvalue 1. This gives a basis {A_1, ..., A_d} for A.

**Step 3 — Define Delta and verify it is injective** (L2056, L2074).

Delta: A -> B(H) is just inclusion: Delta(A) = A for A ∈ A ⊆ B(H). Verify injectivity: Delta is injective by definition.

**Step 4 — Define Gamma: B(M) -> A** (L2061–2065).

Define Lambda: B(M) -> B(H) by Lambda(Y) = Phi(J_M Y J_M†). Then Lambda factors through A (since Im(Phi) = A), so Lambda = Delta ∘ Gamma for a unique Gamma: B(M) -> A. Since Delta is injective, Gamma is well-defined by: for each Y ∈ B(M), Gamma(Y) is the unique element of A such that Delta(Gamma(Y)) = Phi(J_M Y J_M†).

*Matrix algorithm*: form the superoperator matrix for Y |-> Phi(J_M Y J_M†): B(M) -> B(H). Its range lies in A. Express in terms of the basis of A obtained in Step 2 to get the matrix of Gamma.

**Step 5 — Verify the factorization identities** (L2072–2081).

Check Phi(X) = Lambda(C_M(X)) = Phi(Pi_M X Pi_M) (this is eq. PhiX_M, L1740, proved at L1742).
Check Gamma C_M Delta = 1_A: since Delta(A) = Phi(Delta(A)) = Delta(Gamma C_M Delta(A)) and Delta is injective, Gamma C_M Delta = 1_A (L2074).

**Step 6 — Identify w = C_M Delta: A -> B(M) as a *-homomorphism** (L2084–2088).

Set Psi = C_M ∘ Lambda = w ∘ Gamma. Show Psi is idempotent (L2084). Then Im(w) = Im(Psi) = Ker(1 - Psi) ⊆ B(M), and lem_idemp (L1916, Step 7 below) shows w carries the Choi-Effros product to operator product. Hence Im(w) is a C* subalgebra of B(M), and w is a *-homomorphism (L2088).

**Step 7 — Use lem_idemp to establish multiplication** (lem_idemp, L1916).

lem_idemp says: for X,Y ∈ A = Im(Phi), Phi(XY)|_M = (XY)|_M = X|_M Y|_M. This shows w(A ∘_{CE} B) = w(A) w(B) where ∘_{CE} is the Choi-Effros product X ∘_{CE} Y = Phi(XY) on A. Thus w intertwines the C* algebra structure on A (Choi-Effros) with ordinary matrix multiplication on Im(w) ⊆ B(M).

**Step 8 — Identify the C* algebra structure of A by prop_hom_structure** (L2093–2096).

Since w: A -> B(M) is a C* algebra inclusion and dim M < ∞, apply prop_hom_structure (L259, Shard A): there exist Hilbert spaces L_j, E_j and isometries W_j: M -> L_j ⊗ E_j (with W unitary M -> ⊕_j L_j ⊗ E_j) such that
```
w(A_1,...,A_m) = Σ_j W_j† (A_j ⊗ 1_{E_j}) W_j.
```
This uniquely determines A ≅ ⊕_j B(L_j) as a C* algebra. The W_j are computed by simultaneous block-diagonalization of Im(w).

**Step 9 — Determine Delta on M⊥ by eq. Delta_structure** (L2098–2103, eq. Delta_structure).

For A ∈ A, Delta(A)|_M = w(A) (known). Delta(A)|_{M⊥} = Sigma(A) for a UCP map Sigma: A -> B(M⊥). The full form is:
```
Delta(A) = Σ_j J_M W_j† (A_j ⊗ 1_{E_j}) W_j J_M†  +  J_{M⊥} Sigma(A) J_{M⊥}†.
```
Sigma is an arbitrary UCP map A -> B(M⊥); it is constrained by requiring Phi = Delta Gamma C_M.

**Step 10 — Determine Gamma by prop_Gamma** (L2106–2113).

prop_Gamma gives the unique form of any Gamma: B(M) -> A satisfying Gamma w = 1_A:
```
Gamma_j(X) = Tr_{E_j}( W_j X W_j† (1_{L_j} ⊗ gamma_j) )
```
for some density matrices gamma_j on E_j.

**Step 11 — Recover dual Enc/Dec** (L2159).

Enc = (Gamma C_M)* has form eq. Enc (L280); Dec = Delta* has form eq. Dec (L286). Verify Dec Enc = 1_{A*} and Enc Dec = Phi*.

---

## Results

### prop_KLHG — Proposition "Kitaev-Lashkevich-Hastings-Gottesman" (statement L1635–1651, proof L1654–1687)

- **Statement (faithful):** Let J: L -> K be an isometry, w: B(L) -> B(H) a *-homomorphism, and (G, u, V) the canonical Stinespring of Phi: X |-> w(J† X J). Then:
  - Diagram 1 (L1638–1643): u(X) V = V w(J† X J) for all X ∈ B(K), i.e., the canonical representation factors as Phi = (Y |-> V† Y V) ∘ u.
  - Diagram 2 (L1644–1650, specific to canonical Stinespring): u(J Y J†) = V w(Y) V† for all Y ∈ B(L).
  The second diagram says the canonical Stinespring representation has an additional "reverse" compatibility: the left multiplication by J Y J† in the Stinespring space is intertwined by V with w(Y) on H.

- **Proof:** L1654–1687. First diagram: immediate from Stinespring form Phi(X) = V† u(X) V and hypothesis. Second diagram: show s'(Y) = u(JYJ†) equals s''(Y) = V w(Y) V† by checking inner products against spanning vectors B⊗eta. Compute <B⊗eta | u(JYJ†)(A⊗xi)> = eta† Phi(B† J Y J† A) xi = eta† w(J†B†J) w(Y) w(J†AJ) xi (L1665–1673). Compute <B⊗eta | V Z V†(A⊗xi)> = eta† Phi(B†) Z Phi(A) xi for any Z by rank-1 linearity argument (L1675–1685). Substitute Z = w(Y) and equate.

- **Constructive in finite dim?:** Yes. Given explicit matrices for J, w (on generators), and V (from the canonical Stinespring construction), both diagrams can be verified numerically. The canonical Stinespring is constructed as follows: form G = B(K) ⊗ H (dim = dim(K)^2 * dim(H)); define the Gram matrix M_{(A,xi),(B,eta)} = eta† Phi(B†A) xi; eigendecompose to find the null space; quotient; this is finite-dimensional and fully algorithmic.

- **Depends on:** Stinespring theorem (L1621), canonical construction (L1625–1633).
- **Used by:** th_idemp_structure proof context (structure of Delta, Gamma), lem_idemp proof (via eq. eq_idemp).
- **Arbitrary precision needed?:** No for verification; yes if numerically constructing the canonical Stinespring (Gram matrix eigendecomposition requires arb for certified null-space identification).
- **Cross-check / oracle:** Verify diagram 1: compute V† u(X) V and w(J†XJ) as explicit matrices; compare. Verify diagram 2: compute u(JYJ†) and V w(Y) V† as matrices; compare.
- **Implementation notes / risks:** The canonical Stinespring in finite dim is just a Gram matrix eigendecomposition; LAPACK dsyev/zheev suffices. Risk: near-zero eigenvalues of Gram matrix require certified thresholding (arb interval arithmetic) to identify the null subspace exactly.

---

### lem_carrier — Lemma "Carrier as image support" (statement L1724–1725, proof L1727–1729)

- **Statement (faithful):** Let Phi: B(K) -> B(H) be a UCP map with Choi rep (F, V). The carrier of Phi (= support of Phi*(rho_0) for any full-rank rho_0 on H) is the smallest subspace M ⊆ K such that M ⊗ F ⊇ Im(V). Equivalently, M' ⊇ carrier iff M' ⊗ F ⊇ Im(V).

- **Proof:** L1727–1729. Let Pi_{M'} be the projection onto M'. M' ⊇ carrier iff Tr((1-Pi_{M'}) Phi*(rho_0)) = 0 iff Phi(1 - Pi_{M'}) = 0 (since rho_0 has full rank and Phi(1-Pi_{M'}) ≥ 0) iff Z†Z = 0 for Z = ((1-Pi_{M'}) ⊗ 1_F)V iff Im(V) ⊆ M' ⊗ F. The carrier is the smallest such M', i.e., the projection onto the H-marginal of Im(V).

- **Constructive in finite dim?:** Yes. Algorithm: form V as a kh × h matrix (K-index first, F-index second); reshape each column of V into a k × h matrix; the carrier M is the row-span of these matrices viewed as maps F -> K, i.e., M = Im(Pi_K) where Pi_K is the partial trace projection. Concretely: form the k × k matrix Q = V_reshaped * V_reshaped†  (summing over F and H indices) = Σ_{a,alpha} v_{ia,alpha} v*_{ja,alpha}; M = column space of Q. Note Q = Phi(1_K): since Phi(1_K) = V†(1_K ⊗ 1_F)V = 1_H ... wait, that's not right. Actually carrier = support of Phi*(rho_0); the correct matrix is: write V as a (dim K * dim F) × dim H matrix; form P = V V† (acting on K ⊗ F); then M = Im( Tr_F(P) ) where Tr_F is partial trace over F. Compute eigendecomposition of Tr_F(V V†) to get M.

- **Depends on:** Choi representation (L1568), definition of carrier (L1722).
- **Used by:** cor_carrier (L1731), lem_idemp (L1916, L1926), th_idemp_structure proof (L2056).
- **Arbitrary precision needed?:** Yes — finding M requires identifying zero vs. nonzero eigenvalues of Tr_F(VV†). Use arb interval arithmetic to certify eigenvalue clustering.
- **Cross-check / oracle:** Verify: for any rho_0 full-rank on H, compute rho_1 = Phi*(rho_0) by rho_1 = Tr_H(V (rho_0)^T V† ⊗ 1_K) (Choi duality); support of rho_1 should equal M.
- **Implementation notes / risks:** Partial trace is a standard BLAS operation (contraction over tensor indices). The key risk is correct identification of near-zero eigenvalues; arb needed for certified rank.

---

### cor_carrier — Corollary "Carrier annihilation" (statement L1731–1732, proof L1735–1737)

- **Statement (faithful):** With M = carrier of Phi (Choi rep (F,V)): X ∈ B(K) annihilates all vectors in M (i.e., Ker X ⊇ M) iff (X ⊗ 1_F) V = 0.

- **Proof:** L1735–1737. "Ker X ⊇ M" means M ⊆ Ker X, and "(Ker X) ⊗ F ⊇ Im(V)" is the condition from lem_carrier. These are equivalent: Ker X ⊇ M iff (Ker X) ⊗ F ⊇ M ⊗ F ⊇ Im(V) (using lem_carrier for the last step in reverse). One-line proof by substituting M' = Ker X into lem_carrier.

- **Constructive in finite dim?:** Yes. Numerically: form (X ⊗ 1_F) V as a matrix and check norm = 0 (or below threshold). This is a single matrix multiply.
- **Depends on:** lem_carrier (L1724).
- **Used by:** lem_idemp proof (L1926, L1997, L1999), th_idemp_structure proof (L2074).
- **Arbitrary precision needed?:** Same as lem_carrier — need certified zero detection for (X ⊗ 1_F) V.
- **Cross-check / oracle:** Given X in image of Phi (hence X|_M should be nonzero for nontrivial X), confirm (X ⊗ 1_F) V ≠ 0.
- **Implementation notes / risks:** Straightforward BLAS matrix multiply. Numerical tolerance must match that used for lem_carrier eigendecomposition.

---

### lem_idemp — Lemma "Idempotent structure on carrier" (statement L1916–1922, proof L1925–2053)

- **Statement (faithful):** Let H be finite-dim, Phi: B(H) -> B(H) idempotent UCP, A = Im(Phi), M = carrier. Then:
  1. Operators in A do not mix M and M⊥: for all X ∈ A, (1 - Pi_M) X Pi_M = 0.
  2. The Choi-Effros product restricted to M equals the operator product: for X,Y ∈ A,
     ```
     Phi(XY)|_M = (XY)|_M = X|_M Y|_M,
     ```
     where X|_M = J_M† X J_M.

- **Proof:** L1925–2053. Long proof using tensor-network (string diagram) calculus.
  - Part 1 (L1926–1997): set C = ((1-Pi_M) X ⊗ 1_F) V. Compute C†C using idempotency Phi^2 = Phi (written as two applications of V) and eq. eq_idemp (L1777, derived from Phi^2 = Phi), plus ((1-Pi_M) ⊗ 1_F)V = 0 (from lem_carrier, since Im(V) ⊆ M ⊗ F). All four terms in C†C cancel (L1994), so C†C = 0, hence C = 0, hence (1-Pi_M)X Pi_M = 0 by cor_carrier.
  - Part 2 (L1999–2053): use C_M = w, eq. eq_idemp again, and X = Phi(X), Y = Phi(Y), Phi^2 = Phi to show (Phi(XY) ⊗ 1_F)V = (XY ⊗ 1_F)V. Conclude Phi(XY)|_M = (XY)|_M by cor_carrier.

  Key intermediate identity eq_idemp (L1777):
  ```
  (X ⊗ 1_F) V V† (Y ⊗ 1_F) V = V V†(X ⊗ 1_F)(Y ⊗ 1_F) V ... [diagrammatic]
  ```
  derived from Phi^2 = Phi. This is an algebraic identity on V,V† that holds iff Phi is idempotent; it is the core engine of the entire proof.

- **Constructive in finite dim?:** Yes — every step is a matrix identity checkable numerically. The proof is non-constructive in style (showing C†C = 0 implies C = 0) but that implication is constructive.
- **Depends on:** cor_carrier (L1731), eq_idemp (L1777, derived from Phi^2=Phi), Choi representation (L1568).
- **Used by:** th_idemp_structure proof (L2055, L2088–2090), identifying A as C* algebra.
- **Arbitrary precision needed?:** Only if checking identities numerically — arb needed to certify C†C = 0 to machine precision.
- **Cross-check / oracle:** Given A = Im(Phi) explicitly, verify for random X,Y ∈ A: compute Phi(XY) and XY; restrict both to M (by conjugation with J_M); check agreement. Also verify block-diagonal structure of each element of A by checking off-diagonal block (1-Pi_M) X Pi_M = 0.
- **Implementation notes / risks:** eq_idemp (L1777) is the key structural identity — implement and test it independently first. The string diagram proof is long but reduces to repeated matrix multiplications using V, V†, Pi_M.

---

### th_idemp_structure — Theorem "Structure of idempotent UCP maps" (statement L318–332, proof L2055–2091)

See "The idempotent-structure algorithm" section above for the full algorithm. Below: result-entry format.

- **Statement (faithful):** (See L318–332 and the algorithm section above.) Given idempotent UCP Phi: B(H) -> B(H) (finite-dim), there exist M ⊆ H, C* algebra A ≅ ⊕_j B(L_j), and UCP maps Delta: A -> B(H), Gamma: B(M) -> A with:
  - Gamma C_M Delta = 1_A
  - Delta Gamma C_M = Phi
  - w = C_M Delta: A -> B(M) is a *-homomorphism
  - Im(Delta) is block-diagonal w.r.t. H = M ⊕ M⊥

- **Proof method (L2055–2091):** Constructive. Steps:
  1. M = carrier(Phi) (lem_carrier).
  2. A = Im(Phi) as vector space.
  3. Delta = inclusion A ↪ B(H).
  4. Lambda(Y) = Phi(J_M Y J_M†); factored as Delta Gamma.
  5. Phi = Lambda C_M (by PhiX_M, L1740).
  6. Gamma C_M Delta = 1_A from idempotency + injectivity of Delta (L2074).
  7. w = C_M Delta; set Psi = C_M Lambda = w Gamma; show Psi idempotent.
  8. Im(w) is C* subalgebra by lem_idemp Part 2; block-diagonality by lem_idemp Part 1 (L2088–2090).

- **Constructive in finite dim?:** Yes — all steps are explicit matrix computations. See algorithm section.
- **Depends on:** lem_carrier (L1724), cor_carrier (L1731), lem_idemp (L1916), PhiX_M (L1740), prop_hom_structure (L259).
- **Used by:** prop_Gamma (L2106), th_almost_idemp (L2192), th_factorization (L2730). The Enc/Dec structure at L2159 follows directly.
- **Arbitrary precision needed?:** Yes — Steps 1 and 8 require certified eigenvalue clustering (arb).
- **Cross-check / oracle:** Verify Delta Gamma C_M = Phi by computing the superoperator matrix of Delta Gamma C_M and comparing to Phi. Verify w is a *-homomorphism by checking w(AB) = w(A)w(B) for a basis of A.
- **Implementation notes / risks:** The main risk is that A is identified as a subspace of B(H) without an explicit basis; forming the superoperator for Phi and extracting its eigenspace-1 must be numerically stable. Use arb for certified null-space extraction.

---

### prop_Gamma — Proposition "Form of Gamma" (statement L2106–2113, proof L2117–2157)

- **Statement (faithful):** Let w: A -> B(M) be as in eq. Aw (L2094–2096): w(A_1,...,A_m) = Σ_j W_j†(A_j ⊗ 1_{E_j})W_j, where W = (W_1,...,W_m): M -> ⊕_j L_j ⊗ E_j is unitary. Any UCP map Gamma: B(M) -> A = ⊕_j B(L_j) satisfying Gamma w = 1_A has the form:
  ```
  Gamma_j(X) = Tr_{E_j}( W_j X W_j† (1_{L_j} ⊗ gamma_j) )
  ```
  for some density matrices gamma_j on E_j. (This is the general form of a conditional expectation in finite dimensions, per L2115.)

- **Proof:** L2117–2157. Write each Gamma_j in Choi form: Gamma_j(X) = L_j†(X ⊗ 1_{F_j})L_j with L_j: L_j -> M ⊗ F_j isometry. Substitute into Gamma w = 1_A; the condition becomes: for M_{kj} = (W_k ⊗ 1_{F_j}) L_j, one has Tr(M_{kj}†(A_k ⊗ 1_{E_k ⊗ F_j})M_{kj} B_j) = Tr(A_j B_j) for all A,B ∈ A. This forces R_{kj}†R_{kj} = delta_{kj} 1_{L_j} Tr_{L_j} (L2147), where R_{kj} is a "partial isometry" map L_j ⊗ L_k* -> E_k ⊗ F_j. Solution: R_{jj} = xi_j Tr_{L_j} for some xi_j ∈ E_j ⊗ F_j (unit vector), R_{kj} = 0 for k ≠ j. Hence M_{kj} = delta_{kj}(1_{L_j} ⊗ xi_j), so L_j = (W_j† ⊗ 1_{F_j})(1_{L_j} ⊗ xi_j). Plugging into Gamma_j formula and taking partial trace gives gamma_j = Tr_{F_j}(xi_j xi_j†).

- **Constructive in finite dim?:** Yes. Given W_j (from prop_hom_structure applied to w), the allowed Gamma_j are parametrized by density matrices gamma_j on E_j. To find a valid Gamma, choose any gamma_j (e.g., gamma_j = 1_{E_j}/dim(E_j)); then compute Gamma_j explicitly.
- **Depends on:** prop_hom_structure (L259), Choi representation of Gamma_j.
- **Used by:** th_idemp_structure (dual maps Enc, Dec, L2159), th_almost_idemp context.
- **Arbitrary precision needed?:** No — given W_j and gamma_j exactly (or to working precision), Gamma_j is an explicit trace-and-multiply formula.
- **Cross-check / oracle:** Verify Gamma w = 1_A explicitly: for each basis element A_j ∈ B(L_j) ⊆ A, compute Gamma(w(A)) and check = A. This is a finite set of matrix identity checks.
- **Implementation notes / risks:** The partial trace Tr_{E_j}(W_j X W_j† (1_{L_j} ⊗ gamma_j)) is a standard quantum channel (CPTP map); implement as a Kraus operator sum. The gamma_j are free parameters encoding the "ancilla reset state" — any density matrix works. In application to almost-idempotent maps, gamma_j will be determined (approximately) by the actual state, so this freedom is important.

---

## Key definitions & notation introduced here

| Symbol | Meaning | First appearance |
|--------|---------|-----------------|
| (F, V) | Choi representation of Phi: B(K)->B(H); V: H->K⊗F isometry | L1568 (eq. Choi) |
| Stinespring_def | General Stinespring: Phi(X) = V† u(X) V, u *-hom | L1622 |
| Carrier of Phi | Support of Phi*(rho_0) for full-rank rho_0; = smallest M s.t. M⊗F ⊇ Im(V) | L1722 |
| C_M | Compression map: B(H) -> B(M), X |-> J_M† X J_M | L2059, eq. th_idemp_structure |
| Lambda | B(M) -> B(H), Y |-> Phi(J_M Y J_M†) | L2061 |
| Gamma | B(M) -> A, right inverse of C_M∘Delta up to composition | L2057, eq. Gamma_C_Delta |
| Delta | A -> B(H), inclusion of image | L2056 |
| w = C_M Delta | A -> B(M), *-homomorphism | L2084, eq. Aw |
| Sigma | A -> B(M⊥), UCP map for the M⊥ part of Delta | L2098, eq. Delta_structure |
| gamma_j | Density matrix on E_j parametrizing Gamma | L2113, eq. Gamma |
| Phi†(X)Phi(X) ≤ Phi(X†X) | Schwarz inequality for UCP maps | L1692, eq. PhiXdX |
| ||Phi||_cb ≤ 1 | cb-norm bound for UCP | L1718 |
| X ★ Y = Phi(XY) | Choi-Effros product on Im(Phi) | L340 (Shard A), used in lem_idemp |
| eq_idemp | Key identity from Phi^2=Phi (tensor diagram) | L1777 |
| PhiX_M | Phi(X) = Phi(Pi_M X Pi_M), carrier localization | L1740 |

---

## Dependencies OUT / INTO

**Into this shard (used here, proved elsewhere):**
- prop_hom_structure (L259, Shard A): structure of *-homomorphisms; used at L2093.
- th_idemp_structure statement (L318, Shard A): the theorem proved here.
- Enc (L280), Dec (L286): forms recovered at L2159.

**Out of this shard (proved here, used later):**
- lem_carrier (L1724): used by lem_idemp (L1926), th_idemp_structure proof, th_almost_idemp (§12).
- cor_carrier (L1731): used by lem_idemp proof (L1926, L1997), th_idemp_structure proof.
- lem_idemp (L1916): used by th_idemp_structure proof (L2088–2090).
- th_idemp_structure proof (L2055): closes the deferred proof from L318; result used by th_almost_idemp (L2192), th_factorization (L2730).
- prop_Gamma (L2106): used by th_almost_idemp context, Enc/Dec at L2159.
- Choi representation (F,V) (L1568): used throughout §12 (almost-idempotent maps).
- prop_KLHG (L1635): used in §12 via Stinespring factorization arguments.

---

## Open implementation questions / escalations

1. **Minimal Choi representation vs. canonical**: the canonical Choi uses F = K* ⊗ H (dim = k*h, potentially large). For practical computation, minimize: eigendecompose C_Phi and keep only eigenvectors with nonzero eigenvalue. This gives a Kraus rank-r representation with r ≤ k*h. What is r for generic Phi? Needs investigation for the almost-idempotent case.

2. **Carrier computation robustness**: lem_carrier requires identifying Im(V) ⊆ H ⊗ F. Near-idempotent Phi will have a carrier that is only approximately well-defined (eigenvalues of Tr_F(VV†) nearly 0 or 1 but not exactly). The threshold for calling an eigenvalue "zero" must be certified using arb and tied to the almost-idempotency parameter eta. Need to decide: certify carrier exactly (requires interval arithmetic + gap condition), or track carrier approximation error as a separate parameter.

3. **Identifying A = Im(Phi) numerically**: forming the superoperator of Phi as a (dim H)^2 × (dim H)^2 complex matrix and finding its eigenspace-1 is a standard eigenvalue problem, but near-idempotent Phi will have eigenvalues near 1 rather than exactly 1. Need a gap condition. Escalate to §12 (th_almost_idemp) which is specifically about this.

4. **eq_idemp (L1777) as a standalone test**: the identity (Phi^2 = Phi <=> eq_idemp holds) should be implemented as an explicit matrix-level test and verified first before attempting the full th_idemp_structure algorithm.

5. **Sigma (the M⊥ part of Delta, L2098)**: Sigma is an arbitrary UCP map A -> B(M⊥). In the exact idempotent case, it is determined by Phi itself (via Delta(A)|_{M⊥} = Phi projected to M⊥). But in the almost-idempotent case, the choice of Sigma will affect the quality of the approximate decomposition. This requires clarification in §12.

6. **prop_KLHG applicability**: prop_KLHG (L1635) applies to canonical Stinespring of maps of the form X |-> w(J†XJ). In §12, Phi is only approximately of this form. Does prop_KLHG extend approximately? Not addressed in §11; escalate to §12.
