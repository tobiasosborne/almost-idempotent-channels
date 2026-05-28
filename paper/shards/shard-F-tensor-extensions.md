# Shard F — Tensor extensions (L1447–1561)

## Scope

§10 "Tensor extensions" (sec_tens_ext), L1447–1561. Upgrades the main theorem (th_main, L460) from a plain O(eps)-isomorphism to an *extended* O(eps)-isomorphism: the map v: B → A and all its amplifications 1_{M_n} ⊗ v: M_n ⊗ B → M_n ⊗ A are simultaneously O(eps)-isomorphisms, with a constant independent of n. This is exactly the cb-norm / operator-space control that §12 (UCP maps) needs: a UCP map Φ is completely positive precisely because 1_{M_n} ⊗ Φ is positive for every n, so bounding all ampliations is not optional.

---

## Results

### def:opspace — Definition "Operator space" (L1453–1464)

- **Statement (faithful):** A complex vector space L is an *operator space* if each M_n ⊗ L (n = 1,2,...) carries a norm ‖·‖_n satisfying two Ruan axioms:
  - (R1) ‖AXB‖_n ≤ ‖A‖ ‖X‖_k ‖B‖   for A ∈ M_{n,k}, B ∈ M_{k,n}, X ∈ M_k ⊗ L.
  - (R2) ‖diag(X,Y)‖_{k+n} = max{‖X‖_k, ‖Y‖_n}   for X ∈ M_k ⊗ L, Y ∈ M_n ⊗ L.
  The n=1 case recovers the norm on L. L is *self-adjoint* if there is a conjugate-linear involution † preserving all ‖·‖_n.

  Consequences (L1467–1475): The inclusion M_{n'} ⊗ L ↪ M_n ⊗ L (n' ≤ n, pad with zeros) is isometric. The norm ‖·‖_{n',n''} on M_{n',n''} ⊗ L is defined by embedding into M_n ⊗ L for any n ≥ max{n',n''}. The map X ↦ I_n ⊗ X is isometric: ‖I_n ⊗ X‖_n = ‖X‖.

- **Proof:** Definition only; no proof. The Ruan representation theorem (Rua88, Paulsen Thm 13.4) says conversely that every operator space embeds isometrically into B(H) for some H, and in the self-adjoint case the involution is represented by Hermitian conjugation.

- **Constructive in finite dim?:** Yes. In finite dim, L ≅ C^d; each M_n ⊗ L ≅ M_{nd} (after fixing a basis). All norms ‖·‖_n are operator norms on matrix algebras. The Ruan axioms are checkable as linear-algebraic conditions. No infinite-dimensional machinery needed.

- **Depends on:** Paulsen [Paulsen], Ruan [Rua88] (background only; not algorithmically needed).
- **Used by:** Definition L1477 (extended eps-C* algebra), prop_inc_ext (L1483), lem_approx_ext (L1508), th_main_ext (L1538).

- **Arbitrary precision needed?:** No, the definition is structural. But norm computation for ‖·‖_n requires computing operator norms of block matrices — standard LAPACK (largest singular value), certified with arb if needed.

- **Cross-check / oracle:** Verify Ruan axioms (R1), (R2) on explicit examples: B(H) with ‖·‖_n = operator norm of n×n matrices over B(H).

- **Implementation notes / risks:** In practice, L is always a subspace of some M_d; the norms ‖·‖_n are inherited as largest singular values of nd×nd matrices. No extra data structure needed beyond the ambient matrix algebra.

---

### def:extended_eps_cstar — Definition "Extended eps-C* algebra / extended delta-homomorphism" (L1477–1479)

- **Statement (faithful):** An *extended eps-C* algebra* is a complete self-adjoint operator space A with multiplication and unit such that each M_n ⊗ A is an eps-C* algebra (i.e., all the eps-approximate algebra axioms hold at each ampliation level). An *extended delta-homomorphism* v: A' → A'' is a linear map such that for every n, the ampliation 1_{M_n} ⊗ v: M_n ⊗ A' → M_n ⊗ A'' is a delta-homomorphism. An extended delta-inclusion / extended delta-isomorphism are defined via prop_inc_ext (L1483) and bijectivity.

- **Proof:** Definition only.

- **Constructive in finite dim?:** Yes. In finite dim, A = C^d; checking "each M_n ⊗ A is an eps-C* algebra" reduces to checking the eps-C* conditions on nd × nd matrices for every n. But "every n" is infinite — see cb-norm truncation question below.

- **Depends on:** def:opspace (L1453), eps-C* algebra definition (L435–L455, Shard A).
- **Used by:** prop_inc_ext (L1483), lem_approx_ext (L1508), th_main_ext (L1538).

- **Arbitrary precision needed?:** No (structural), but the truncation question is real.

- **Cross-check / oracle:** Any exact C* algebra is trivially an extended 0-C* algebra (all ampliations are exact C* algebras).

- **Implementation notes / risks:** The "for each n" quantifier is the central implementation risk. See Open Questions below.

---

### prop_inc_ext — Proposition "Extended inclusion from lower bound" (statement L1483–1484, proof L1486–1506)

- **Statement (faithful):** Let v: A' → A'' be an extended delta-homomorphism of eps-C* algebras. If ‖v(X)‖ ≥ η‖X‖ for some η > 2δ and all X ∈ A', then for all n and all X ∈ M_n ⊗ A':
  ```
  ‖(1_{M_n} ⊗ v)(X)‖_n ≥ a ‖X‖_n,   a = 1 − O(δ + eps),
  ```
  where a is independent of n.

- **Proof:** L1486–1506. Method: define a_n = inf{‖(1_{M_n} ⊗ v)(X)‖_n / ‖X‖_n : X ≠ 0}. The sequence a_n is non-increasing (L1491). Key step: show a_{2n} ≥ a_n / 2 by decomposing a 2n×2n block matrix X into four n×n blocks, using the operator space axioms to bound ‖X‖_{2n} ≤ 2 ‖X‖_{n,max} (where ‖X‖_{n,max} = max_{p,q} ‖X_{pq}‖_n), and using a_n on each block. Then apply prop_delta_hominc (L1194) at level n=1: a_1 ≥ η > 2δ implies a_1 ≥ 1 − δ' for δ' = O(δ+eps). By induction on n (using a_{2n} ≥ a_n/2 and the smallness of δ,eps), a_n ≥ 1 − δ' for all n.

  Adapts prop_delta_hominc (L1194–1196): the scalar (n=1) case is proved there; prop_inc_ext lifts it to all n uniformly by the doubling argument.

- **Constructive in finite dim?:** Partial. The *proof* uses induction over all n; for implementation, the bound a_n ≥ 1 − δ' is established for all n by the inductive structure, but numerically one only needs finite n. See truncation question.

- **Depends on:** prop_delta_hominc (L1194), def:opspace (L1453), def:extended_eps_cstar (L1477).
- **Used by:** th_main_ext (L1538, via cor_improvement adapted to extended inclusions, L1557).

- **Arbitrary precision needed?:** No (the constant a = 1 − O(δ+eps) is explicit; computing a_n for specific n is a singular-value computation).

- **Cross-check / oracle:** At n=1 reduces to prop_delta_hominc. Check: for an exact *-homomorphism (δ=eps=0) that is injective, a_n = 1 for all n (isometry of *-homomorphisms).

- **Implementation notes / risks:**
  - The doubling argument (a_{2n} ≥ a_n/2) is the key: it shows the sequence cannot decay faster than geometrically, and combined with the n=1 lower bound, pins a_n ≥ 1−δ' uniformly.
  - For a constructive check of the extended isomorphism property, compute a_n for n = 1, 2, 4, ..., N_max using the block decomposition; see truncation question for N_max.
  - LAPACK primitive: largest and smallest singular values of (nd × nd) matrices.

---

### lem_approx_ext — Lemma "Approximate extended homomorphism" (statement L1508–1510, proof L1512–1536)

- **Statement (faithful):** For any extended delta-homomorphism v from a finite-dimensional C* algebra B to an eps-C* algebra A, there is an extended O(eps)-homomorphism ṽ: B → A such that the ampliations v_n = 1_{M_n} ⊗ v and ṽ_n = 1_{M_n} ⊗ ṽ satisfy:
  ```
  ‖ṽ_n(X) − v_n(X)‖ ≤ O(δ) ‖X‖   for all X ∈ M_n ⊗ B, all n.
  ```

- **Proof:** L1512–1536. Adapts lem_approx (L1224–1316) to the ampliated setting. The core extension:
  1. Define multiplicativity defect of v_n: g_n(X,Y) = v_n(XY) − v_n(X)v_n(Y) for X,Y ∈ M_n ⊗ B.
  2. Use the diagonal D = Σ_j A_j ⊗ B_j ∈ B ⊗ B (computed from B, not from A; standard Pauli diagonal, L1249) to define the correcting map:
     ```
     w_n'(X) = Σ_j v_n(I_n ⊗ A_j) g_n(I_n ⊗ B_j, X).
     ```
  3. The key generalization (L1521–1534): extend the identity Σ_j v(XA_j) g(B_j,Y) = Σ_j v(A_j) g(B_j X, Y) from X,Y ∈ B to X,Y ∈ M_n ⊗ B. This is done by passing to matrix elements [·]_{kp} ∈ B and applying the scalar identity at each matrix entry, using the diagonal property Σ_j Z A_j ⊗ B_j = Σ_j A_j ⊗ B_j Z for scalar Z. The rest of the proof is identical to lem_approx (L1256ff).
  4. The output ṽ = v + w satisfies ‖ṽ_n(X) − v_n(X)‖ ≤ O(δ)‖X‖ and ṽ_n is an O(eps)-homomorphism for each n.

  Note: the correction w_n' depends on n but is derived from a single diagonal D ∈ B ⊗ B and the single map v. Thus ṽ (defined at level n=1) automatically generates ṽ_n = 1_{M_n} ⊗ ṽ with the required properties.

- **Constructive in finite dim?:** Yes. The diagonal D is explicit (Pauli/generalized Pauli, L1249–1254). The correction formula is a finite sum of matrix products. Computing ṽ requires:
  1. Compute D = Σ_j A_j ⊗ B_j for B (closed form, L1249).
  2. Compute g = multiplicativity defect of v on a basis of B.
  3. Compute the linear correction w via the diagonal formula (lem_approx, L1273ff adapted).
  All operations are matrix multiplications + norm estimates.

- **Depends on:** lem_approx (L1224), diagonal construction (L1245–1254), def:extended_eps_cstar (L1477).
- **Used by:** th_main_ext (L1538, L1557): cor_improvement must be adapted to extended inclusions using lem_approx_ext and prop_inc_ext.

- **Arbitrary precision needed?:** No. The correction is a finite algebraic computation over C. Arb interval arithmetic for certifying ‖w‖ ≤ O(δ)‖X‖.

- **Cross-check / oracle:** At n=1 reduces exactly to lem_approx. Check that ṽ_n = 1_{M_n} ⊗ ṽ for all n (the extended structure is automatic once ṽ is defined at n=1).

- **Implementation notes / risks:**
  - The multiplicativity defect g_n at level n is fully determined by g at level n=1: g_n(I_n ⊗ X, I_n ⊗ Y) = I_n ⊗ g(X,Y). The w_n' formula only uses I_n ⊗ A_j, so it reduces to a block-diagonal action. No new data is needed beyond v and D.
  - FLINT/arb: compute the diagonal D once (closed form for M_d: d^{-2} Σ_{j,k} S_{jk}^† ⊗ S_{jk}). The correction w is a single matrix (n=1 computation); its ampliation is automatic.

---

### th_main_ext — Theorem "Extended main theorem" (statement L1538–1540, proof L1542–1558)

- **Statement (faithful):** For any finite-dimensional extended eps-C* algebra A, there exist a C* algebra B and an extended O(eps)-isomorphism v: B → A. (The implicit constant in O(eps) is independent of A and its dimension.)

- **Proof:** L1542–1558. Adapts the proof of th_main (L1414–1444) to the extended setting. Structure:

  1. **New objects constructed for the original A** (L1542): nontrivial delta-projections, compression maps, subspaces — all built from A at the scalar level, exactly as in th_main. The tensor extension to I_n ⊗ P is free: if P^† = P and ‖P^2 − P‖ ≤ δ in A, then (I_n ⊗ P)^† = I_n ⊗ P and ‖(I_n ⊗ P)^2 − (I_n ⊗ P)‖ = ‖I_n ⊗ (P^2 − P)‖ = ‖P^2 − P‖ ≤ δ in M_n ⊗ A.

  2. **Compression maps** (L1544): Co_{P,Q}: A → A extends to 1_{M_n} ⊗ Co_{P,Q} = Co_{I_n⊗P, I_n⊗Q} with the same error bounds, because I_n ⊗ P is a delta-projection in the eps-C* algebra 1_{M_n} ⊗ A. Similarly, S_{P,Q} extends to M_n ⊗ S_{P,Q} = S_{I_n⊗P, I_n⊗Q}.

  3. **One-dimensional projections and Hilbert spaces** (L1546–1555): lem_PQ_Hilb (L1123) gives a Hermitian inner product on S_{P,Q} when Q is one-dimensional. The extended space C^n ⊗ S_{P,Q} = M_{n,1} ⊗ S_{P,Q} satisfies the same inner product formula (L1547–1548):
     ```
     Y^† · X = ⟨Y,X⟩ Q̃,   where ⟨Y,X⟩ = Σ_l ⟨[Y]_{l1}, [X]_{l1}⟩.
     ```
     The norm equivalence ‖X‖_Euc = √⟨X,X⟩ ≈ ‖X‖ (up to 1±O(eps+δ)) still holds (L1552–1553). The map Ha^Q_{P,R} (L1555) extends to 1_{M_n} ⊗ Ha^Q_{P,R}: M_n ⊗ S_{P,R} → B(C^n ⊗ S_{R,Q}, C^n ⊗ S_{P,Q}).

  4. **Error reduction** (L1557): cor_improvement (L1317) is adapted to extended inclusions using lem_approx_ext and prop_inc_ext; the improvement step works at all ampliation levels simultaneously.

  5. **Section §9 arguments** (L1557): "require only trivial modifications, namely, one should use the norms ‖·‖_n in certain places." All structural lemmas (lem_merging, lem_add_dim, lem_extension, cor_merge_sum) carry through because their proofs are norm inequalities that hold level-by-level.

  Connection to th_main (L460): th_main_ext strictly strengthens th_main — its output v is not merely a delta-isomorphism but an *extended* one. The same B works; the proof constructs B identically and shows the existing v has cb-norm control.

  Connection to §9 proof (L1414): The §9 proof is reprised step by step; the only change is that every use of ‖·‖ is replaced by ‖·‖_n and the ampliation 1_{M_n} ⊗ (·) is tracked throughout.

- **Constructive in finite dim?:** Yes, with the truncation caveat. The algorithm is: run the th_main algorithm (which is already constructive per Shard E / §9) and verify that the output v satisfies the extended conditions. The construction of v does not change; only the correctness certificate must check all n. See truncation question.

- **Depends on:** th_main (L460) and its §9 proof (L1414), lem_approx_ext (L1508), prop_inc_ext (L1483), cor_improvement (L1317), lem_PQ_Hilb (L1123), lem_merging (L1325), lem_add_dim (L1363), lem_extension (L1378), cor_merge_sum (L1352).
- **Used by:** §12 (UCP maps, sec_channels, L1562ff) — the extended isomorphism is what allows the UCP realization to be completely positive, not merely positive.

- **Arbitrary precision needed?:** Yes, for certifying the O(eps) bound uniformly over n. The bound itself is independent of n, so a single certified computation at a sufficiently large N_max suffices (see truncation question).

- **Cross-check / oracle:** At n=1 reduces to th_main. For an exact C* algebra (eps=0), the extended 0-isomorphism is just a *-isomorphism, which is automatically completely isometric.

- **Implementation notes / risks:**
  - The theorem says B and v exist; constructively, B = ⊕_j M_{d_j} (determined by th_main) and v is the map constructed in §9. No new object is created.
  - The extended conditions are verified post-hoc by checking ‖(1_{M_n} ⊗ v)(X)‖_n / ‖X‖_n ≥ a and the extended eps-C* conditions on M_n ⊗ A for n = 1, ..., N_max.
  - Coupling to §12 (UCP maps): A UCP map Φ: B(K) → B(H) corresponds (via Choi-Stinespring, L1568) to an isometry V: H → K ⊗ F. Complete positivity of Φ is equivalent to positivity of 1_{M_n} ⊗ Φ for all n. The extended isomorphism from th_main_ext is exactly what is needed to transport this structure through the approximate algebra.

---

## Key definitions & notation introduced here (operator space, matrix norms)

| Symbol | Meaning | First appearance |
|--------|---------|-----------------|
| ‖·‖_n | norm on M_n ⊗ L (the n-th matrix norm of L) | L1454 |
| (R1), (R2) | Ruan axioms for operator spaces | L1456, L1460 |
| ‖X‖_{n',n''} | norm on M_{n',n''} ⊗ L, defined by embedding into M_n ⊗ L | L1473 |
| ‖X‖_{n,max} | max_{p,q} ‖X_{pq}‖_n for a block matrix X | L1496 |
| extended eps-C* algebra | operator space A with each M_n ⊗ A an eps-C* algebra | L1477 |
| extended delta-homomorphism | v s.t. 1_{M_n} ⊗ v is a delta-homomorphism for all n | L1478 |
| extended delta-inclusion / -isomorphism | as above + lower bound / bijectivity | L1481 |
| v_n, ṽ_n | shorthand for 1_{M_n} ⊗ v, 1_{M_n} ⊗ ṽ | L1509 |
| g_n | multiplicativity defect of v_n: g_n(X,Y) = v_n(XY)−v_n(X)v_n(Y) | L1515 |
| a_n | lower isometry constant: inf ‖v_n(X)‖/‖X‖ | L1489 |
| I_n | identity in M_n (abbrev. for I_{M_n}, introduced L1451) | L1451 |
| 0_{n,k} | zero n×k matrix | L1451 |

**cb-norm context:** The completely bounded norm ‖v‖_{cb} = sup_n ‖1_{M_n} ⊗ v‖. The theorem asserts ‖1_{M_n} ⊗ v‖ and ‖(1_{M_n} ⊗ v)^{-1}‖ are both 1+O(eps) uniformly in n, i.e., ‖v‖_{cb} ≤ 1+O(eps) and ‖v^{-1}‖_{cb} ≤ 1+O(eps). This is strictly stronger than the n=1 bound of th_main.

---

## Dependencies OUT / INTO

**INTO this section (used here):**
- th_main / proof (L460 / L1414): blueprint for th_main_ext proof.
- prop_delta_hominc (L1194): scalar lower-bound propagation; used inside prop_inc_ext.
- lem_approx (L1224): scalar version; generalized by lem_approx_ext.
- cor_improvement (L1317): error reduction; adapted to extended inclusions at L1557.
- lem_PQ_Hilb (L1123): Hermitian inner product on S_{P,Q}; extended at L1546.
- lem_merging (L1325), cor_merge_sum (L1352), lem_add_dim (L1363), lem_extension (L1378): all used in §9 proof; carried over with trivial modification.
- Diagonal D for B (L1245–1254): explicit Pauli diagonal, used in lem_approx_ext.
- Ruan [Rua88], Paulsen §13: background on operator spaces.

**OUT from this section (used later):**
- th_main_ext (L1538): used by §12 (sec_channels, L1562) for UCP map realization — the extended isomorphism provides the complete positivity control.
- def:extended_eps_cstar (L1477): structural definition referenced throughout §12.
- prop_inc_ext (L1483): the "extended inclusion" characterization needed whenever §12 invokes lower-bound arguments on UCP maps.

---

## Open implementation questions / escalations

1. **cb-norm truncation — which finite n suffices?**
   The theorem asserts bounds for *all* n. In practice, for a d-dimensional algebra A, the structure is determined at n ≤ d (any map M_n ⊗ B → M_n ⊗ A factors through the structure at n ≤ d by a matrix-rank argument). More precisely: if B = ⊕_j M_{d_j} with Σ d_j = D, then any X ∈ M_n ⊗ B has rank at most n·D as an operator, and the operator space structure of a D-dimensional space L stabilizes at n = D (by the completely bounded Hahn-Banach theorem / Wittstock). **Implementation question:** What is the explicit N_max (in terms of dim A) such that verifying ‖(1_{M_n} ⊗ v)(X)‖_n ≥ a‖X‖_n for n ≤ N_max certifies the uniform bound? The paper does not state this explicitly; it must be derived from the Ruan representation theorem. **Escalation:** This needs a separate analysis before implementation. Candidate answer: N_max = dim(A) suffices for the Ruan representation, but a rigorous finite-n certificate is needed. Flag for escalation to the formal verification layer.

2. **Extended eps-C* check:** Verifying that A is an extended eps-C* algebra (all M_n ⊗ A are eps-C* algebras) requires checking the eps-C* conditions at each n. Same truncation question applies: what finite n suffices? For exact C* algebras (eps=0), n=1 suffices (since C* norms are unique). For eps > 0, unclear — needs analysis.

3. **Computing the matrix norms ‖·‖_n:** For a concrete A ⊆ M_d, the norm ‖X‖_n for X ∈ M_n ⊗ A is the operator norm of X viewed as an nd × nd matrix. This is the largest singular value — computable in O((nd)^3) via LAPACK/arb. Cost grows as n^3 d^3, which is manageable for small n but expensive if N_max is large.

4. **Coupling to §12 UCP maps:** The extended isomorphism from th_main_ext feeds into the Choi-Stinespring construction (L1568). The implementation must verify that the isometry V in the Choi representation (L1572) is certified to be approximately isometric, with error O(eps), for the completely positive approximation. This requires the cb-norm bound from th_main_ext as input. The interface between this shard and shard for §12 must pass: (a) the map v at n=1, (b) the certified uniform bound a = 1 − O(eps), and (c) the diagonal D (for the correction formula).

5. **Inductive doubling argument in prop_inc_ext:** The induction a_{2n} ≥ a_n/2 combined with a_1 ≥ 1−δ' requires (1−δ')/2 > 2δ, i.e., δ < (1−δ')/4, which holds for small δ,eps. In an implementation with certified δ and eps, this condition must be checked explicitly before invoking the bound. If δ is too large, the induction breaks and the result does not apply — the algorithm must flag this as a failure mode.
