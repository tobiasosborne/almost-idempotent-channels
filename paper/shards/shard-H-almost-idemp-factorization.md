# Shard H — Almost-idempotent UCP maps & factorization (L2162–2900)

## Scope

§12 "Almost-idempotent UCP maps" (L2162–2900), two subsections:
- §12.1 "The associated eps-C* algebra" (L2164–2725): th_almost_idemp and its multi-part proof.
- §12.2 "Approximate factorization through a C* algebra" (L2726–2900): th_factorization, supporting constructions, lem_RC.

This is the headline result of the paper: an eta-idempotent UCP map Phi admits an approximate factorization through a genuine C* algebra B, with all maps UCP and accuracy O(eta).

---

## End-to-end pipeline (the headline algorithm)

Input: UCP map Phi: B(H) -> B(H) with ||Phi^2 - Phi||_cb <= eta, eta < 1/4. H finite-dimensional.

**Step 1. Regularize Phi to an exact idempotent (in the CB algebra).**

Compute Phi_tilde = theta(2*Phi - 1) = (1/2)(1 + sgn(2*Phi - 1)) via the Taylor expansion
  Phi_tilde = (1/2)(1 + (2*Phi - 1)(1 - 4(Phi - Phi^2))^{-1/2})
converging for eta < 1/4. This uses prop_P (L524) applied to X = 2*Phi - 1.
Properties: Phi_tilde^2 = Phi_tilde, ||Phi_tilde - Phi||_cb = O(eta), Phi_tilde(I) = I, Phi_tilde(X^dag) = Phi_tilde(X)^dag.

**Step 2. Define the eps-C* algebra A = Img(Phi_tilde).**

A = Ker(1 - Phi_tilde) subset B(H); closed subspace, contains I, invariant under *-involution.
Equip A with:
- norm, unit, involution inherited from B(H),
- multiplication X star Y = Phi_tilde(XY) for X, Y in A  [Choi-Effros product, L2188-2189].

By th_almost_idemp (L2192): (A, star, ||.||, dag, I) is an extended O(eta)-C* algebra.
Key step powering this: the associativity equations (Phi_assoc1)/(Phi_assoc2), proved using Choi or Stinespring calculus (see th_almost_idemp proof parts below).

**Step 3. Find a genuine C* algebra B and an extended O(eta)-isomorphism v: B -> A.**

Apply th_main_ext (L1538) to A (a finite-dimensional extended O(eta)-C* algebra).
th_main_ext invokes the full incremental construction from th_main (L460): nontrivial projection search (prop_H-group, lem_nontriv_projection), merging/extension (lem_merging, lem_extension), error reduction (cor_improvement), all generalized to tensor extensions.
Output: B = direct sum of full matrix algebras B(L_j), and v: B -> A with ||v||_cb, ||v^{-1}||_cb = 1 + O(eta) and v approximately preserving the star-product structure in the CB sense.

**Step 4. Construct the approximate encoding map Delta: B -> B(H).**

Set Delta_tilde = inclusion_of_A_in_B(H) after v: B -> A -> B(H). Not UCP.
CP-ize: for each component j, use the unitary 1-design D_j = sum_s p_{js} U_{js}^dag tensor U_{js} on B(L_j) (L2771-2783) to define
  Delta'(X) = sum_s p_s Phi(Delta_tilde(X U_s^dag) Delta_tilde(U_s))
which is manifestly CP, commutes with *, and satisfies ||Delta' - Delta_tilde||_cb = O(eta) (uses tilde_Del2, L2761).
Unitalize: Delta(X) = (Delta'(I_B))^{-1/2} Delta'(X) (Delta'(I_B))^{-1/2}.
Delta is UCP with ||Delta - Delta_tilde||_cb = O(eta).

**Step 5. Construct the approximate decoding map Upsilon: B(H) -> B.**

Use Choi representation of Phi: Phi(X) = V^dag (X tensor 1_F) V.
Use Choi representation of Delta: Delta(X) = sum_j W_j^dag (X_j tensor 1_{E_j}) W_j, where W_j: H -> L_j tensor E_j, sum_j W_j^dag W_j = I_H.
Define auxiliary operators L_j: L_j -> H tensor F via
  L_j = sum_s p_{js} (Delta(U_{js}^dag) tensor 1_F) V W_j^dag (U_{js} tensor xi_j)
where xi_j in E_j is a unit vector with ||C_j xi_j|| approx 1 (xi_j from lem_RC, L2840).
Define CP map Upsilon'_j(X) = L_j^dag (Phi(X) tensor 1_F) L_j; unitalize similarly:
  Upsilon(X) = (Upsilon'(I_H))^{-1/2} Upsilon'(X) (Upsilon'(I_H))^{-1/2}.
Upsilon is UCP with ||Upsilon - Upsilon_tilde||_cb = O(eta).

**Step 6. Verify the factorization.**

Output: UCP maps Delta: B -> B(H) and Upsilon: B(H) -> B satisfying:
  (a) ||Delta Upsilon - Phi||_cb = O(eta)   [eq. DelUps, L2733]
  (b) ||Upsilon_n(Delta_n(X) Delta_n(Y)) - XY|| <= O(eta)||X||||Y||  [eq. UpsDel2, L2736]
      In particular, ||Upsilon Delta - 1_B||_cb = O(eta).

**Step 7 (dual picture). Decode/encode channels on states.**

Phi^* is a quantum channel (CPTP map) on density matrices on H.
Dec = Delta^*: states on B(H) -> states on B.
Enc = Upsilon^*: states on B -> states on B(H).
Then: Dec Enc = (Upsilon Delta)^* approx 1_{B^*} and Enc Dec = (Delta Upsilon)^* approx Phi^*.
This is the approximate factorization Phi^* = Enc Dec up to O(eta) in cb norm.

---

## Results

### th_almost_idemp — Theorem 12.1 "Associated eps-C* algebra" (statement L2192–2193)

- **Statement (faithful):**
  Let H be an arbitrary nonzero Hilbert space, Phi: B(H) -> B(H) a UCP map with ||Phi^2 - Phi||_cb <= eta. Let Phi_tilde = theta(2*Phi - 1) (exact idempotent, ||Phi_tilde - Phi||_cb = O(eta)), A = Img(Phi_tilde), and X star Y = Phi_tilde(XY) for X, Y in A.
  Then (A, star, ||.||, dag, I) is an extended O(eta)-C* algebra; specifically:
  - ||X star I - X|| = ||I star X - X|| = 0 (unit exact).
  - (X star Y)^dag = Y^dag star X^dag (involution exact).
  - ||X star Y|| <= (1 + O(eta))||X||||Y|| (approximate submultiplicativity).
  - ||(X star Y) star Z - X star (Y star Z)|| <= O(eta)||X||||Y||||Z|| (approximate associativity).
  - ||X^dag star X|| >= (1 - O(eta))||X||^2 (approximate C*-identity lower bound).
  The constant in O(eta) is universal (no dependence on dim H or A).

- **Proof — Part 1: Using eqs (Phi_assoc1) and (Phi_assoc2) (L2208–2237).**
  Assumes (Phi_assoc1) and (Phi_assoc2) hold for all eta-idempotent UCP maps and their tensor extensions 1_{M_n} tensor Phi.
  - Unit and involution axioms: trivial from Phi_tilde(I) = I and *-preservation.
  - Approximate submultiplicativity: ||X star Y|| = ||Phi_tilde(XY)|| <= ||Phi(XY)|| + O(eta)||XY|| <= (1 + O(eta))||X||||Y|| using ||Phi_tilde - Phi||_cb = O(eta) and ||Phi||_cb = 1.
  - Approximate associativity: replace Phi by Phi_tilde in (Phi_assoc1) and (Phi_assoc2) (error O(eta)), subtract to get ||(Phi_tilde(Phi_tilde(X)Phi_tilde(Y))Phi_tilde(Z)) - (Phi_tilde(X Phi_tilde(Phi_tilde(Y)Phi_tilde(Z))))|| <= O(eta)||X||||Y||||Z||. For X, Y, Z in A, Phi_tilde acts as identity, giving associativity defect O(eta).
  - Approximate C*-identity: from inequality Phi(X^dag)Phi(X) <= Phi(X^dag X) [PhiXdX, L1692], replacing Phi with Phi_tilde: ||X^dag star X|| >= (1 - O(eta))||X||^2.
  The same argument applies to 1_{M_n} tensor Phi_tilde, establishing the extended version.

- **Proof — Part 2: (Phi_assoc1) and (Phi_assoc2) in the finite-dimensional case (L2239–2585).**
  Uses Choi representation Phi(X) = V^dag(X tensor 1_F)V with Pi = VV^dag projection, 1 - Pi is the "leakage" projector.
  Key intermediate result (L2273–2296): ||(1-Pi) u_1(Phi(X)) V^(two copies)|| = O(sqrt(eta))||X||. Shown by computing the square of its norm = Phi^2(Phi(X^dag)Phi(X)) - Phi(Phi^2(X^dag)Phi^2(X)) = O(eta)||X||^2 (using ||Phi^2 - Phi||_cb <= eta).
  Key tensor diagram W at L2385–2422: a diagram with two (1-Pi) insertions around Y is O(eta)||X||||Y||||Z||. Expanding (1-Pi) = 1 - VV^dag gives a telescoping sum of four terms (L2425–2582) that collapses to Phi(Phi(Phi(X)Phi(Y))Phi(Z)) - Phi(Phi(X)Phi(Y)Phi(Z)) + O(eta). Combined with W = O(eta), this gives (Phi_assoc1). Equation (Phi_assoc2) obtained similarly.

- **Proof — Part 3: Remark (L2628–2637).**
  Discusses how to extend the Choi-based argument to infinite dimensions using the Stinespring "stack": inductive construction H_0 = H, with Stinespring representations (H_n, u_n, V_n) such that Phi(X) = V_{n+1}^dag u_{n+1}(X) V_{n+1} and u_n(V_n^dag X V_n) = V_{n+1}^dag u_{n+1}(X) V_{n+1}. Eqs (uVdV, L2601) and (uVVd, L2606); multi-step versions (uVdV-multi, uVVd-multi, L2617–2625). Introduces H_infty as completion of directed union.

- **Proof — Part 4: General proof of (Phi_assoc1) and (Phi_assoc2) (L2639–2724).**
  Replaces the finite-dim (1-Pi)-insertion by R_k(X) = u_{2+k,1}((I_1 - V_1 V_1^dag) u_1(Phi(X))) V_{2+k} V_{1+k} (L2642-2644). Shows ||R_k(X)||^2 = O(eta)||X||^2 via computation using (uVdV-multi) (L2646-2660). Analogue of diagram W at L2670-2672: W = V_1^dag R_1(X^dag)^dag u_{3,0}(Phi(Y)) V_3 R_0(Z) = O(eta)||X||||Y||||Z||. Expanded form of W (L2675–2723) reduces to Phi(Phi(Phi(X)Phi(Y))Phi(Z)) - Phi(Phi(X)Phi(Y)Phi(Z)) + O(eta)||X||||Y||||Z||, giving (Phi_assoc1). Uses only (uVdV) and (uVdV-multi), not (uVVd).

- **Constructive in finite dim?:** Yes. All steps are explicit matrix computations:
  - Phi as Choi matrix: O(n^4) complex numbers.
  - Compute Phi^2 - Phi to verify eta bound.
  - Compute Phi_tilde via truncated power series in (Phi - Phi^2); series truncated at precision O(eta^K) for desired accuracy.
  - A = Img(Phi_tilde) computed as range of Phi_tilde (column space); O(n^2) dimensional.
  - Verify eps-C* axioms numerically. Eps = O(eta) with explicit constants from the proof.

- **Depends on:** prop_P (L524, for the theta-function idempotentization), PhiXdX (L1692, for C*-identity bound), prop_KLHG (L1635, for the Stinespring stack in the general proof).

- **Used by:** th_factorization (L2730, A is the intermediate algebra), all of §12.2.

- **Arbitrary precision needed?:** Yes. The Taylor series for (1 - 4(Phi - Phi^2))^{-1/2} must be computed to sufficient depth; arb/FLINT ball arithmetic needed to certify convergence and the O(eta) bound constants.

- **Cross-check / oracle:**
  - eta = 0: Phi exactly idempotent. Phi_tilde = Phi exactly. A = Img(Phi) is a genuine C* algebra (Choi-Effros theorem). Axiom bounds all become 0.
  - Verify: ||(X star Y) star Z - X star (Y star Z)|| computed directly on test elements should be <= C*eta*||X||||Y||||Z|| for an explicit C.

---

### th_factorization — Theorem 12.2 "Approximate factorization" (statement L2730–2739)

- **Statement (faithful):**
  Let H be a nonzero finite-dimensional Hilbert space, Phi: B(H) -> B(H) UCP with ||Phi^2 - Phi||_cb <= eta. Then there exist a finite-dimensional C* algebra B and UCP maps
    Delta: B -> B(H),  Upsilon: B(H) -> B
  such that:
    ||Delta Upsilon - Phi||_cb <= O(eta),                                          [DelUps, L2733]
    ||Upsilon_n(Delta_n(X) Delta_n(Y)) - XY|| <= O(eta)||X||||Y||   for X,Y in M_n tensor B.    [UpsDel2, L2736]
  In particular [L2739]: ||Upsilon Delta - 1_B||_cb <= O(eta).

- **Proof: "Discussion and an outline of the proof" (L2742–2769). THIS PROOF IS AN OUTLINE. Gaps flagged below.**

  Outline structure:
  1. (L2743) Form Phi_tilde and A = Img(Phi_tilde) as extended O(eta)-C* algebra (th_almost_idemp).
  2. (L2748) Apply th_main_ext to get B and extended O(eta)-isomorphism v: B -> A. This is asserted but th_main_ext's proof (L1542) is itself a substantial construction; all constants must be made explicit.
  3. (L2749) Define Delta_tilde = (v followed by inclusion A -> B(H)) and Upsilon_tilde = (Phi_tilde with target A, followed by v^{-1}). State properties: tilde_DelUps (L2750-2752), tilde_Del1 (L2756), tilde_Del2 (L2759-2761).
  4. (L2768) Assert: "To complete the proof, we will approximate Delta_tilde and Upsilon_tilde by UCP maps Delta and Upsilon in the completely bounded norm." The actual constructions are given in the subsequent text (L2771-2899) but are NOT part of the labelled proof block.

  Construction of Delta (L2771–2801):
  - B = direct_sum_j B(L_j). Define unitary 1-designs D_j via Pauli-type construction (diag_j01, diag_j2, L2772–2778).
  - Delta'(X) = sum_s p_s Phi(Delta_tilde(X U_s^dag) Delta_tilde(U_s)); manifestly CP (L2791-2795).
  - ||Delta' - Delta_tilde||_cb = O(eta) from tilde_Del2.
  - Delta(X) = (Delta'(I_B))^{-1/2} Delta'(X) (Delta'(I_B))^{-1/2} is UCP with ||Delta - Delta_tilde||_cb = O(eta).
  - Key properties of Delta: Delta_norm (L2806), PhiDelta1 (L2808), PhiDelta2 (L2810-2812), PhiDelta3 (L2813-2815).
  - Choi representation of Delta: L2831–2837.

  Construction of Upsilon (L2859–2899):
  - Requires lem_RC (L2840) to get C_j in B(E_j) with 1 - O(eta) <= ||C_j|| <= 1.
  - Choose xi_j in E_j unit vector with ||C_j xi_j|| near 1.
  - L_j = sum_s p_{js} (Delta(U_{js}^dag) tensor 1_F) V W_j^dag (U_{js} tensor xi_j).
  - Upsilon'_j(X) = L_j^dag (Phi(X) tensor 1_F) L_j; manifestly CP.
  - Show ||Upsilon' Delta - 1_B||_cb = O(eta) by explicit calculation (L2871–2893), using PhiDelta1, PhiDelta3, lem_RC.
  - Upsilon_tilde approx: show Upsilon' approx Upsilon' Phi approx Upsilon' Phi_tilde = Upsilon' Delta_tilde Upsilon_tilde approx Upsilon' Delta Upsilon_tilde approx Upsilon_tilde with O(eta) accuracy.
  - Unitalize: Upsilon(X) = (Upsilon'(I_H))^{-1/2} Upsilon'(X) (Upsilon'(I_H))^{-1/2} is UCP.

  **GAPS AND OUTLINE ISSUES:**
  1. The "proof" block L2742–2769 ends with "we will approximate..." without executing the approximation. The constructions in L2771–2899 are presented as continuous prose, not flagged as proof steps. There is no explicit closure statement proving (DelUps) and (UpsDel2) from those constructions.
  2. Constant tracking: O(eta) constants from th_main_ext (which itself invokes th_main via a deferred proof chain through prop_H-group, multiple lemmas) are not made explicit. The headline O(eta) in th_factorization is the composition of O(eta) from th_almost_idemp, O(eta) from th_main_ext, and O(eta) from the CP-ization steps. The composite constant is unspecified.
  3. The choice of xi_j ("a unit vector such that ||C_j xi_j|| - 1 = O(eta)") is asserted to exist via lem_RC but no algorithm is given. For implementation: find the leading singular vector of C_j.
  4. The unitary 1-design construction is referenced to equation (Pauli_diag) elsewhere in the paper; the explicit construction must be located and implemented.
  5. (UpsDel2) is deduced from the chain Upsilon approx Upsilon_tilde and the properties of v^{-1}, but the tensorized version (for all n) is only partially verified. The proof says "The same relations hold for u_n and Delta_n" (L2747) but does not recheck UpsDel2 explicitly.
  6. The entire pipeline depends on th_main_ext (L1538) being fully constructive. That theorem's proof (L1542–1557) says "we adapt the proof of the main theorem" with several assertions left to "straightforward modifications". The main theorem proof (L1414) is itself a long incremental construction; all constants must be tracked.

- **Constructive in finite dim?:** Partial. The explicit constructions (Phi_tilde, A, Delta', Upsilon') are all finite-dimensional linear-algebraic operations. The gap is in th_main_ext: the incremental C* algebra construction requires: finding nontrivial projections (non-constructive in general via Lefschetz-Hopf; for finite dim may use SVD/eigenvalue methods instead), merging/extension steps, error reduction (cohomological). These need constructive replacements.

- **Depends on:** th_almost_idemp (L2192), th_main_ext (L1538), prop_P (L524), lem_RC (L2840), PhiDelta1/2/3 properties of Delta.

- **Used by:** This is the headline result; no further results depend on it in the paper.

- **Arbitrary precision needed?:** Yes. The CP-ization step uses (Delta'(I_B))^{-1/2} and (Upsilon'(I_H))^{-1/2}, both require certified inversion. All O(eta) bounds must be tracked with explicit constants for a verified implementation.

- **Cross-check / oracle:**
  - eta = 0: Phi exactly idempotent. th_idemp_structure (L318) gives exact B (= A as a C* algebra), Delta = the inclusion Delta from th_idemp_structure, Upsilon = Gamma_Co_M from th_idemp_structure, with Delta Upsilon = Phi exactly and Upsilon Delta = 1_B exactly.
  - Verify: ||Delta Upsilon - Phi||_cb <= C*eta and ||Upsilon Delta - 1_B||_cb <= C*eta with explicit C.
  - Verify: Phi^* = Enc Dec up to O(eta) in diamond/cb norm on states.

- **Implementation notes / risks:**
  - The proof outline does not label which lines establish each of (DelUps) and (UpsDel2); implementer must reconstruct the proof closure from prose.
  - The map v: B -> A from th_main_ext has no polynomial-time algorithm stated; the incremental construction of B is a sequence of pivoting steps whose count depends on dim(A).
  - The unitary 1-design step (L2771) requires an explicit example for each B(L_j); for M_n, the Heisenberg-Weyl (Pauli) group works [referenced in paper as eq. (Pauli_diag), location TBD in paper].
  - PhiDelta3 (L2813) requires the three-fold product estimate; its proof invokes (Phi_assoc1) applied to 1_{M_n} tensor Phi, which requires the tensor extension of th_almost_idemp — check that the O(eta) constants survive tensorization.

---

### lem_RC — Lemma 12.3 "Structure of R_j" (statement L2840–2846, proof L2848–2857)

- **Statement (faithful):**
  Given Choi representation Delta(X) = sum_j W_j^dag (X_j tensor 1_{E_j}) W_j with W_j: H -> L_j tensor E_j and sum_j W_j^dag W_j = I_H, and the unitary 1-design p_{js}, U_{js} on B(L_j), define
    R_j = sum_s p_{js} (U_{js}^dag tensor 1_{E_j}) W_j W_j^dag (U_{js} tensor 1_{E_j}) in B(L_j tensor E_j).
  Then:
  (i) R_j = 1_{L_j} tensor C_j for some C_j in B(E_j).
  (ii) 1 - O(eta) <= ||C_j|| <= 1.

- **Proof (L2848–2857):**
  (i) R_j commutes with X tensor 1_{E_j} for all X in B(L_j) because of the 1-design property (diag_j2). By Schur's lemma / commutant characterization in finite dimensions, R_j = 1_{L_j} tensor C_j.
  (ii) Upper bound: ||C_j|| = ||R_j|| <= sum_s p_{js} ||(U_{js}^dag tensor 1_{E_j}) W_j W_j^dag (U_{js} tensor 1_{E_j})|| <= ||W_j||^2 <= 1 (since W_j is a partial isometry with W_j^dag W_j = partial of I_H).
  Lower bound: Compute Phi(W_j^dag R_j W_j) = sum_s p_{js} Phi(Delta(U_{js}^dag) Delta(U_{js})) = sum_s p_{js} Delta(U_{js}^dag U_{js}) + O(eta) = Delta(1_{L_j}) + O(eta) using (PhiDelta2). Then ||Phi(W_j^dag R_j W_j)|| >= ||Delta(1_{L_j})|| - O(eta) >= 1 - O(eta) by (Delta_norm). Since ||Phi(W_j^dag R_j W_j)|| <= ||C_j||, we get ||C_j|| >= 1 - O(eta).

- **Constructive in finite dim?:** Yes. R_j is an explicit sum of conjugates of W_j W_j^dag; the commutant decomposition R_j = 1 tensor C_j is found by tracing over L_j: C_j = (1/dim L_j) Tr_{L_j}(R_j). Bounds verified numerically.

- **Depends on:** Choi representation of Delta (L2831-2837), 1-design property (diag_j2, L2776-2778), PhiDelta2 (L2810), Delta_norm (L2806).

- **Used by:** th_factorization construction of Upsilon (L2859).

- **Arbitrary precision needed?:** Yes, for certified lower bound on ||C_j||; needed to certify invertibility of C_j^dag C_j (used at L2891) and to bound ||xi_j|| after normalization.

- **Cross-check / oracle:**
  - eta = 0: R_j = (1/dim L_j) 1_{L_j} tensor Tr_{L_j}(W_j W_j^dag). If Phi is exactly idempotent with carrier M, then sum_j W_j^dag W_j = I_M, and W_j W_j^dag is the range projection of W_j in L_j tensor E_j; C_j is the partial trace, and ||C_j|| = 1.

---

## Key definitions & notation introduced here

- **eta-idempotent UCP map (L2166-2169):** Phi: B(H) -> B(H) UCP with ||Phi^2 - Phi||_cb <= eta. The cb-norm ||Phi^2 - Phi||_cb = sup_n ||(1_{M_n} tensor (Phi^2 - Phi))(X)||/||X||.

- **Phi_tilde (L2171-2175, eq. tilde_Phi):** Phi_tilde = theta(2*Phi - 1) = (1/2)(1 + (2*Phi-1)(1 - 4(Phi - Phi^2))^{-1/2}). Exact idempotent in the CB-algebra, O(eta)-close to Phi. NOT guaranteed UCP (may not be completely positive).

- **A = Img(Phi_tilde) (L2184-2186):** The associated approximate algebra. Closed subspace of B(H), contains I, invariant under dag.

- **Star product (Choi-Effros product) on A (L2188-2189, eq. Choi-Effros):** X star Y = Phi_tilde(XY) for X, Y in A. Reduces to the standard Choi-Effros product X star Y = Phi(XY) when Phi is exactly idempotent. This is the key construction: Phi_tilde supplies the retraction onto A, and using it to define multiplication gives approximate associativity from the approximate idempotency.

- **(Phi_assoc1) (L2198-2200):**
  Phi(Phi(Phi(X)Phi(Y))Phi(Z)) = Phi(Phi(X)Phi(Y)Phi(Z)) + O(eta)||X||||Y||||Z||.
  Meaning: the image of Phi is approximately closed under product in B(H), and the "nested" vs "flat" triple products agree to O(eta). This is what makes the star product approximately associative.

- **(Phi_assoc2) (L2202-2204):**
  Phi(Phi(X)Phi(Phi(Y)Phi(Z))) = Phi(Phi(X)Phi(Y)Phi(Z)) + O(eta)||X||||Y||||Z||.
  Right-associative version of the same.

- **Stinespring stack (L2587-2626):** Inductive sequence H_0 = H, H_1, H_2, ... with isometries V_n: H_{n-1} -> H_n and *-homomorphisms u_n. Used in the general (infinite-dim) proof of (Phi_assoc1)/(Phi_assoc2). In finite dim, this stack collapses to the Choi representation with tensor powers of F.

- **Delta_tilde, Upsilon_tilde (L2749):** Non-UCP approximations to Delta and Upsilon, satisfying Delta_tilde Upsilon_tilde = Phi_tilde, Upsilon_tilde Delta_tilde = 1_B, ||Delta_tilde||_cb = 1 + O(eta), ||Upsilon_tilde||_cb = 1 + O(eta).

- **Unitary 1-design D_j (L2771-2784):** A probability distribution p_{js} over unitaries U_{js} on L_j such that sum_s p_{js} X U_{js}^dag tensor U_{js} = sum_s p_{js} U_{js}^dag tensor U_{js} X for all X in B(L_j). Used to CP-ize Delta_tilde. Explicit example: Heisenberg-Weyl (Pauli) group, referenced in paper as eq. (Pauli_diag).

- **R_j (L2843-2845, eq. R_def):** R_j = sum_s p_{js} (U_{js}^dag tensor 1_{E_j}) W_j W_j^dag (U_{js} tensor 1_{E_j}). Equals 1_{L_j} tensor C_j. Key object linking the Choi structure of Delta to the 1-design averaging.

- **L_j (L2862-2863):** L_j = sum_s p_{js} (Delta(U_{js}^dag) tensor 1_F) V W_j^dag (U_{js} tensor xi_j). The Kraus-like operator defining Upsilon'_j.

---

## Dependencies OUT / INTO

**Depends on (incoming):**
- prop_P (L524): theta-function / sign-function calculus for approximate projections; gives Phi_tilde.
- PhiXdX (L1692): Kadison inequality Phi(X^dag)Phi(X) <= Phi(X^dag X); gives C*-identity lower bound.
- th_main_ext (L1538): existence of B and v: B -> A as extended O(eta)-isomorphism; the single most expensive dependency.
- prop_KLHG (L1635): commutativity of Stinespring diagrams; used in general proof of (Phi_assoc1)/(Phi_assoc2).
- lem_idemp (L1916): exact idempotent case: Choi-Effros product on Img(Phi) coincides with operator product restricted to carrier M; provides the eta=0 oracle.
- th_idemp_structure (L318): complete structure of exact idempotent UCP maps; provides the eta=0 baseline for the full pipeline.
- prop_Gamma (L2106): form of UCP maps Gamma satisfying Gamma w = 1_A (used in §11, provides context for Upsilon).
- prop_hom_structure (L259): structure of *-homomorphisms into B(M); gives the W_j decomposition of Delta.

**Produces (outgoing):**
- th_almost_idemp -> th_factorization (A is the intermediate algebra).
- th_factorization: the end-to-end result. Nothing further in the paper depends on it; it is the paper's headline theorem.
- lem_RC -> construction of Upsilon in th_factorization.

---

## Open implementation questions / escalations

1. **th_main_ext constructivity (ESCALATE).** The headline pipeline stalls at Step 3 if th_main_ext cannot be made constructive. The proof (L1542) says "we adapt the proof of th_main". The proof of th_main uses: (a) existence of nontrivial projections via a Lefschetz-Hopf topological argument (non-constructive over general fields); (b) incremental construction with merging/extension/error-reduction steps. For finite-dimensional matrices, the nontrivial projection may be found constructively (e.g., compute spectral gap, find approx projection via eigenvalue decomposition), but this needs verification against the paper's bounds.

2. **Constant tracking.** The composite O(eta) in th_factorization is a product of constants from: (a) prop_P, (b) (Phi_assoc1)/(Phi_assoc2) proofs, (c) th_main_ext, (d) Delta' construction, (e) Upsilon' construction. None of these individual constants are given explicitly. For certified arb/FLINT arithmetic, all must be extracted and composed. This requires re-reading all dependent proofs.

3. **Proof closure of th_factorization.** The labelled proof block ends at L2769 without proving (DelUps) and (UpsDel2). The subsequent constructions (L2771-2899) must be read as the continuation of the proof; implementer must identify exactly which lines establish each numbered equation.

4. **Unitary 1-design implementation.** The explicit 1-design for B(L_j) with dim L_j = d requires the Heisenberg-Weyl group (d^2 elements, each a Pauli-like unitary). Location of (Pauli_diag) in the paper must be identified; the construction needs to be extracted for each block j.

5. **Choice of xi_j.** lem_RC guarantees ||C_j|| >= 1 - O(eta) but does not specify xi_j algorithmically. For implementation: compute SVD of C_j, take xi_j = leading right singular vector. Verify that the resulting ||C_j xi_j||^2 = largest singular value of C_j^dag C_j is indeed >= 1 - O(eta).

6. **UpsDel2 for tensor extensions.** The paper proves Upsilon' Delta approx 1_B at n=1 (L2871-2893) and asserts "The same relations hold for u_n and Delta_n" (L2747). The tensorized UpsDel2 must be verified explicitly; the argument requires that all steps (PhiDelta1, PhiDelta3, lem_RC) tensorize correctly, which they do since they all involve cb norms, but this should be made rigorous.

7. **Invertibility of Delta'(I_B) and Upsilon'(I_H).** Both unitalization steps require these operators to be invertible. The paper argues ||Delta'(I_B) - I_H||, ||Upsilon'(I_H) - I_B|| = O(eta), so invertibility holds for small enough eta. An explicit eta threshold must be derived (likely eta < 1/C for some explicit C).
