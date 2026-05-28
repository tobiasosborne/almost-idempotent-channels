# Shard A — Motivation & main results (L251–492)

## Scope

Covers §1 "Motivation: How to recognize an encoded subsystem?" (L251–401) and §2 "Main theorem on eps-C* algebras and the proof strategy" (L403–491). §1 establishes the exact structure of *-homomorphisms and idempotent UCP maps as motivation, defines the Choi-Effros product, and introduces eps-idempotent UCP maps. §2 defines eps-Banach/*eps-Banach/eps-C* algebras and delta-homomorphisms/inclusions/isomorphisms, states the main theorem (th_main), and gives an informal proof-strategy overview (cohomology roadmap, why naive Haar fails, and the incremental construction outline). All deferred proofs are in §9 (th_main) and §11.3 (th_idemp_structure).

---

## Results

### prop_hom_structure — Proposition "*-homomorphism structure" (statement L259–271, proof L273–275 inline)

- **Statement (faithful):** Let L_1,...,L_m and M be finite-dimensional Hilbert spaces, A = ⊕_j B(L_j), and w: A → B(M) a *-homomorphism. Then there exist auxiliary Hilbert spaces E_j and a unitary W = (W_1,...,W_m): M → ⊕_j (L_j ⊗ E_j) such that
  ```
  w(A_1,...,A_m) = Σ_j W_j† (A_j ⊗ 1_{E_j}) W_j   for all A ∈ A.
  ```
  The homomorphism w is injective iff E_j ≠ 0 for all j.

- **Proof:** L273–275, inline. Method: M decomposes into irreducible representations of A (complement of any invariant subspace is invariant, so full decomposition into irreps follows). Each B(L_j) appears n_j times; set E_j = C^{n_j}. The Wigner-type identification then gives the stated form. Entirely standard representation theory.

- **Constructive in finite dim?:** Yes — already an algorithm. Decompose M into irreps of A by block-diagonalizing the commutant of the image of w (via simultaneous block-diagonalization). Count multiplicities n_j; set E_j = C^{n_j}. The isometries W_j are the columns of the diagonalizing unitary restricted to each block.
  - Inputs: explicit matrices for generators of A (n×n complex), explicit matrix for w on generators.
  - Outputs: integers n_j (dim E_j), isometries W_j (explicit matrices).
  - Key primitive: joint Schur decomposition / simultaneous block-diagonalization of the image of w.

- **Depends on:** Standard finite-dim representation theory of semisimple algebras; no paper-internal dependencies.

- **Used by:** th_idemp_structure (L2093), prop_Gamma (L2107), th_factorization (§12.2). Also the blueprint for the approximate version throughout.

- **Arbitrary precision needed?:** Yes, in the sense that near-degenerate multiplicities (n_j counts) require certified eigenvalue clustering. With arb/FLINT, use interval arithmetic to certify that the commutant's eigenvalues cluster into distinct groups; each group's cardinality gives n_j exactly (an integer).

- **Cross-check / oracle:** (1) Verify W is unitary: ||W†W - I|| < tol. (2) Check reconstruction: ||w(A) - Σ_j W_j†(A_j ⊗ 1)W_j|| < tol for a test set of generators. (3) eta=0 reduction: if A = B(C^n) (single block), output is one W_j = W unitary with E_j = C^{n_1} where n_1 = dim M / n. (4) Two independent implementations: one via commutant diagonalization, one via Wedderburn decomposition.

- **Implementation notes / risks:** Degeneracy risk: if two blocks L_j and L_k have the same dimension (n_j = n_k), standard eigensolvers may mix them. Must use certified clustering (e.g., arb's `acb_mat_eig`) or structure-aware algorithms. LAPACK primitive: `zheev`/`zhegvd` for commutant diagonalization, then `zungtr` for basis. In the arb layer, use `acb_mat_eig_simple` with error bound checks.

---

### Example "Physical example" (L305–312) — UMTC/topological context, not a theorem

Not a theorem or definition with algorithmic content. Describes a spin system on a disk with two holes whose low-energy subspace decomposes as M ≅ ⊕_{a,b,c} V^{ab}_c ⊗ E^c_{ab}, where the fusion spaces V^{ab}_c = C(c, a⊗b) come from a UMTC C. The logical algebra A = ⊕_{a,b,c} B(V^{ab}_c). The idempotent channel T = T_- T_+ is constructed via methods of [SKK19]. This is motivation only; no bounds or constructive content for the implementation.

---

### th_idemp_structure — Theorem "Idempotent UCP map structure" (statement L318–332, proof DEFERRED to §11.3 L2055–2091)

- **Statement (faithful):** Let H be a finite-dimensional Hilbert space and Φ: B(H) → B(H) an idempotent UCP map. Then there exist a subspace M ⊆ H (with inclusion J_M: M → H), a C*-algebra A, and UCP maps Δ: A → B(H) and Γ: B(M) → A such that the diagram
  ```
         B(H)
        ↑     ↓ C_M
  A --Δ--   B(M) --Γ--> A
  ```
  with C_M: X ↦ J_M† X J_M satisfies:
  ```
  Γ C_M Δ = 1_A,    Δ Γ C_M = Φ.
  ```
  Moreover, w = C_M Δ: A → B(M) is a C*-algebra homomorphism, and operators in Im(Δ) are block-diagonal with respect to H = M ⊕ M⊥.

- **Proof:** DEFERRED to §11.3, L2055–2091. Method sketch:
  1. Define M = carrier of Φ (range of the projection determined by Φ), A = Im(Φ) as a vector space, Δ = inclusion map.
  2. Define Λ: X ↦ Φ(J_M X J_M†): B(M) → B(H); factor Λ = ΔΓ via A = Im(Φ).
  3. Use idempotency of Φ to derive Γ C_M Δ = 1_A (L2074).
  4. Show Ψ = C_M Λ = wΓ is idempotent, Im(w) = Im(Ψ) is a C*-subalgebra of B(M) closed under product, unit, and involution (uses lem_idemp L1916).
  5. Identify A with Im(w), making w a C*-algebra inclusion and Γ a UCP map.
  6. Block-diagonality of Im(Δ) is the first part of lem_idemp (L1916).
  Connection to prop_hom_structure: once Im(w) is a C*-subalgebra of B(M), prop_hom_structure describes its explicit decomposition into blocks (L2093–2104), yielding the concrete forms of Δ (L2099–2104) and Γ (prop_Gamma L2106).

- **Constructive in finite dim?:** Yes, but proof is non-constructive in infinite dimensions (uses abstract idempotent UCP theory). In finite dim, the construction is:
  1. Compute eigenvectors of the cb-norm superoperator Φ with eigenvalue 1 to find M and A = Im(Φ). In finite dim, M = range of the projection Π_M = lim_{k→∞} Φ^k (converges in finite dim since Φ is idempotent; directly compute M = ker(I - Φ) via null space of (I - Φ) as a matrix on vec(B(H))).
  2. Form w = C_M Δ as explicit matrices.
  3. Diagonalize Im(w) via prop_hom_structure to get L_j, E_j, W_j.
  4. Construct Γ via prop_Gamma: Γ_j(X) = Tr_{E_j}(W_j X W_j† (1_{L_j} ⊗ γ_j)) for chosen density matrices γ_j.
  - Inputs: Φ as an n²×n² matrix (superoperator on B(H)).
  - Outputs: M (subspace, encoded as projector), A (basis), Δ, Γ as explicit matrices.

- **Depends on:** lem_idemp (L1916), prop_hom_structure (L259), prop_Gamma (L2106).

- **Used by:** prop_Gamma (L2106), the approximate factorization result (§12.2 L2730).

- **Arbitrary precision needed?:** Yes, to certify that M = ker(I - Φ) is the correct subspace (eigenvalue = 1 versus < 1 must be certified). Interval arithmetic on eigenvalues of the (I-Φ) superoperator matrix.

- **Cross-check / oracle:** (1) Check Γ C_M Δ = 1_A and Δ Γ C_M = Φ explicitly. (2) eta=0 case: Φ itself idempotent, so output should satisfy exact equalities. (3) Randomized test: pick random density matrices ρ, check Dec(Enc(ρ)) = ρ where Enc = (Γ C_M)*, Dec = Δ*.

- **Implementation notes / risks:** lem_idemp (L1916) must be implemented first (shard covering §11.3). The carrier M is defined via lem_carrier (L1724). Risk: if Φ is only approximately idempotent, M is not well-defined; this theorem is strictly for exact Φ. The approximate version is the job of th_main + §12.

---

### def:eps-algebras — Definition block "eps-Banach / *eps-Banach / eps-C* algebra" (L407–440)

Not a theorem; a definition block. Key axioms (all for eps ≥ 0 sufficiently small):

**eps-Banach algebra** (L407–414): Banach space A with bilinear multiplication satisfying:
- Sub-multiplicativity: ||XY|| ≤ (1+eps)||X|||Y||  [eq ax_prodnorm, L410–411]
- Approximate associativity: ||(XY)Z - X(YZ)|| ≤ eps ||X|||Y|||Z||  [eq ax_assoc, L412–413]

**eps'-commutativity** (L415–419): ||XY - YX|| ≤ eps'||X|||Y||  [eq ax_comm, L417–418]

***eps-Banach algebra** (L420–424): complex eps-Banach algebra with conjugate-linear involution X ↦ X† satisfying ||X†|| = ||X|| and (XY)† = Y†X†  [eq ax_*, L422–423] (both exact).

**eps-C* algebra** (L425–429): *eps-Banach algebra additionally satisfying:
- Approximate C* identity: ||X†X|| ≥ (1-eps)||X||²  [eq ax_C*, L427–428]
- Note: upper bound ||X†X|| ≤ (1+eps)||X||² follows from ax_prodnorm + ax_*.

**Unit conditions** (L430–439): default approximate unit [eq ax_eps_unit, L432–434]: ||XI - X|| ≤ eps||X||, ||IX - X|| ≤ eps||X||, |||I|| - 1| ≤ eps; or exact unit [eq ax_exact_unit, L435–437]: XI = X, IX = X, ||I|| = 1 (achievable by O(eps)-change of unit and multiplication via prop_unit L672).

---

### def:delta-morphisms — Definition block "delta-homomorphism / delta-inclusion / delta-isomorphism" (L443–456)

Not a theorem; a definition block.

**delta-homomorphism** (L443–450): bounded linear v: A' → A'' (between eps'-Banach and eps''-Banach algebras) with:
- ||v(I) - I|| ≤ delta  [eq hom_unit, L446–447]
- ||v(XY) - v(X)v(Y)|| ≤ delta ||X|||Y||  [eq hom_mult, L448–449]
Non-unital variant: only hom_mult. In *-algebra setting: also v(X†) = v(X)†.

**delta-inclusion** (L451–454): delta-homomorphism with (1-delta)||X|| ≤ ||v(X)|| ≤ (1+delta)||X||.

**delta-isomorphism** (L455): bijective delta-inclusion. Note (L458): inverse of a delta-isomorphism is a (delta + O(delta²))-isomorphism; each O() is a concrete function independent of additional data.

---

### th_main — Theorem "Main theorem on eps-C* algebras" (statement L460–462, proof DEFERRED to §9 L1414–1444)

- **Statement (faithful):** Any finite-dimensional eps-C* algebra A is O(eps)-isomorphic to some C* algebra B. The implicit constant in O(eps) does not depend on A or its dimensionality.

- **Proof:** DEFERRED to §9, L1414–1444. Method sketch (3 stages):
  1. **Stage 1 (L1415–1426): Find maximal commutative C*-subalgebra.** Construct a c_0*eps-inclusion v_comm: B_comm → A where B_comm is a commutative C*-algebra of maximum dimensionality. Show that the images P_j = v_comm(Π_j) of the projection basis are all one-dimensional (if any P_j spans a space of dim > 1, use lem_nontriv_projection L931 to split it, then cor_merge_sum L1352 to merge, then cor_improvement L1317 to reduce error, contradicting maximality).
  2. **Stage 2 (L1430–1441): Incremental matrix-algebra inclusions.** Group the indices j into equivalence classes C where dim(S_{P_j,P_k}) = 1 iff j,k in same class (else 0, by lem_add_dim L1363). For each class C = {1,...,s}, inductively build c_0*eps'-isomorphisms v_r: M_r → A_r (r = 1,...,s) using lem_extension L1378 to add one dimension at a time, then cor_improvement L1317 to reduce error at each step.
  3. **Stage 3 (L1443): Merge across classes.** Successively merge the isomorphisms v_C: M_{|C|} → S_{P_C} for distinct classes C using cor_merge_sum L1352, followed by cor_improvement L1317 each time.
  Final output: a c_0*eps-isomorphism v: B → A where B = ⊕_C M_{|C|}.

- **Constructive in finite dim?:** Yes, but the proof is non-constructive in the key lemmas it depends on (lem_nontriv_projection L931 uses Lefschetz-Hopf theorem; prop_H-group L971 uses H-space homotopy theory). Those lemmas must be replaced by constructive algorithms in finite dim. Constructive route for th_main itself, given constructive subroutines:
  - Inputs: A as an n×n matrix algebra with explicit multiplication table (or basis + structure constants), eps.
  - Stage 1: enumerate projections in A (self-adjoint elements with ||P²-P|| small); find a maximal commutative set; reduce via cor_improvement.
  - Stage 2: for each equivalence class, run the iterative extension loop (lem_extension is constructive in finite dim — explicit matrix construction); apply error reduction after each step.
  - Stage 3: run merging (cor_merge_sum is constructive in finite dim); apply error reduction.
  - Output: explicit isomorphism v (as a linear map, i.e. a matrix), identification of B = ⊕_j M_{n_j}.
  The main non-constructive bottleneck to resolve: lem_nontriv_projection (escalate to shard covering §6 L931).

- **Depends on:** def:eps-algebras (L407), def:delta-morphisms (L443), lem_nontriv_projection (L931), prop_H-group (L971), cor_improvement (L1317), lem_merging (L1325), cor_merge_sum (L1352), lem_add_dim (L1363), lem_extension (L1378).

- **Used by:** th_main_ext (L1538, tensor extension), th_factorization (L2730, approximate factorization of Phi), th_almost_idemp (L2192, approximate idempotent result).

- **Arbitrary precision needed?:** Yes, for two reasons: (1) checking the O(eps) bound constants in cor_improvement requires certified norm computations (cb-norm via SDP); (2) the equivalence class partition in Stage 2 requires certifying that dim(S_{P_j,P_k}) is exactly 0 or 1 (integer decision, requires certified eigenvalue gaps).

- **Cross-check / oracle:** (1) Verify v is a delta-isomorphism: check hom_unit, hom_mult, and near-isometry conditions explicitly for a basis. (2) eta=0 reduction: if A is already a C* algebra (eps=0), the output B should be isomorphic to A with delta=0 (identity map). (3) Low-dimensional sanity check: for A = M_2(C) with eps-perturbed multiplication, reconstruct the standard M_2(C). (4) Two independent implementations: one following the paper's 3-stage algorithm, one using a global optimization (minimize ||v(XY)-v(X)v(Y)|| over linear maps v).

- **Implementation notes / risks:** The constant c_0 from cor_improvement (L1317) must be made explicit (not just existential). The paper says "does not depend on A or its dimensionality" but the actual constant must be extracted from cor_improvement's proof. Risk: the merging step (cor_merge_sum) compounds errors; track the exact error accumulation across Stage 3 iterations. LAPACK/FLINT primitives: SVD for near-isometry checks (`zgesvd`), SDP solver for cb-norm (external: SDPA or SCS), arb for certified norms.

---

### Informal proof-strategy overview (L463–491) — roadmap for implementation

This is not a theorem but the central roadmap. Key points:

1. **Cohomological approach and why it fails naively (L464–484):**
   The associativity defect h(X,Y,Z) = (XY)Z - X(YZ) satisfies the 3-cocycle equation (L470–471). If h were a coboundary — h(X,Y,Z) = X g(Y,Z) - g(XY,Z) + g(X,YZ) - g(X,Y)Z + O(eps²)||X|||Y|||Z|| — then the modified multiplication X·Y = XY + g(X,Y) would be O(eps²)-associative, and iteration would give exact associativity.
   For an exact C* algebra, a **diagonal** D = Σ_j A_j ⊗ B_j ∈ A ⊗ A exists satisfying Σ_j A_j ⊗ B_j X = Σ_j X A_j ⊗ B_j and Σ_j A_j B_j = I. Given D, any cocycle h has coboundary g(X,Y) = Σ_j A_j h(B_j,X,Y). For finite-dim C* algebras the diagonal is D = ∫ dU (U† ⊗ U) (Haar measure on unitary group).
   **Why naive Haar fails:** constructing the diagonal in the eps-associative setting gives error bounds proportional to n = dim A. So the procedure works only if eps < c/n — dimension-dependent, not universal.

2. **Incremental construction overview (L486–490):**
   - Build v: B → A incrementally, maintaining v as an O(eps)-inclusion at each step.
   - **Nontrivial projection (L486, lem_nontriv_projection §6):** To start (B_0 = C → B_1 = C⊕C), need a P ∈ A with P† = P, ||P²-P|| ≤ delta, and P not close to 0 or I. Strategy: construct approximate unitary group U ⊆ A as a C^1 manifold (§5, implicit function theorem). Fixed points of the inversion map sigma: U ↦ U^{-1} on U correspond to U = 2P-I with U² = I. Use **Lefschetz-Hopf theorem** + H-space properties (§6) to count fixed points of sigma; shows more than 2 fixed points exist when dim A > 1, giving a nontrivial P.
   - **Merging and extension (L488):** B is built via a partitioned index set with basis E_{jk} satisfying E_{jk}E_{lm} = delta_{kl}E_{jm}. Process: start with all singletons (B commutative). Each "splitting" step: replace j with {j',j''} and E_{jj} with E_{j'j'} + E_{j''j''}. When splitting impossible, add off-diagonal E_{jk}. Structured via cor_merge_sum (merging) and lem_extension (extension).
   - **Error reduction (L490, cor_improvement §8):** After each extension/merge step, error delta may exceed delta_0 = O(eps). cor_improvement: if a delta-inclusion of C* algebra B into eps-C* algebra A exists (delta below a threshold), there also exists a delta_0-inclusion with delta_0 = O(eps) independent of delta. Proved by cohomological method exploiting B's diagonal. Inspired by [Chr80, Joh88]; simpler here due to finite dimensionality.

---

## Key definitions & notation introduced here

- **C_M** (compression map, L325): X ↦ J_M† X J_M: B(H) → B(M). Used throughout.
- **Choi-Effros product** (L341–343): X ★ Y = Phi(XY) on Im(Phi) for idempotent Phi. Key object showing Im(Phi) is a C* algebra.
- **Completely bounded norm ||·||_cb** (L347–352): sup_n sup_{X≠0} ||(1_{M_n}⊗Λ)(X)||/||X||. Efficiently computable by SDP in finite dim. Dual: diamond norm ||·||_◇.
- **eta-idempotent UCP map** (L354): ||Phi²-Phi||_cb ≤ eta.
- **theta-idempotent approximant** (L355–363): Phi_tilde = theta(2Phi-1) where theta = indicator of positive reals. Well-defined if eta < 1/4; idempotent; ||Phi_tilde - Phi||_cb ≤ O(eta). NOT guaranteed to be completely positive.
- **eps-Banach algebra** (L407–414): axioms ax_prodnorm (L410) and ax_assoc (L412).
- ***eps-Banach algebra** (L420–424): adds exact involution conditions ax_* (L422).
- **eps-C* algebra** (L425–429): adds approximate C* identity ax_C* (L427).
- **Approximate/exact unit** (L430–439): ax_eps_unit (L432) vs ax_exact_unit (L435).
- **delta-homomorphism** (L443–450): hom_unit (L446) + hom_mult (L448).
- **delta-inclusion** (L451–454): delta-homomorphism that is also (1±delta)-isometric.
- **delta-isomorphism** (L455): bijective delta-inclusion.
- **Associativity defect** (L464–467): h(X,Y,Z) = (XY)Z - X(YZ), satisfies 3-cocycle equation (L470–471).
- **Diagonal** (L480): D = Σ_j A_j ⊗ B_j ∈ A⊗A with bimodule and normalization conditions; trivializes Hochschild cohomology.

---

## Dependencies OUT / Dependencies INTO this region

**Dependencies OUT** (results in this region that later shards use):
- prop_hom_structure (L259): used by th_idemp_structure proof (L2093), prop_Gamma (L2107), §12.
- th_idemp_structure (L318): used by prop_Gamma (L2106), §12.2 factorization.
- def:eps-algebras (L407): used by every subsequent section.
- def:delta-morphisms (L443): used by every subsequent section.
- th_main (L460): used by th_main_ext (L1538), th_almost_idemp (L2192), th_factorization (L2730).
- Informal proof strategy (L463–491): roadmap for §5 (unitaries), §6 (nontrivial projection), §8 (error reduction), §9 (proof of th_main).

**Dependencies INTO this region** (things from outside L251–492 that are cited here):
- prop_P (L524): cited at L359 to justify properties of Phi_tilde. Must be available when implementing the theta(2Phi-1) construction.
- prop_unit (L672): cited at L441 for exact-unit normalization.
- lem_nontriv_projection (L931): cited at L486 as key lemma for finding nontrivial projections.
- prop_H-group (L971): cited implicitly at L486 (H-space homotopy for fixed-point count).
- cor_improvement (L1317): cited at L490 as error reduction subroutine.
- cor_merge_sum (L1352): cited at L488 as merging procedure.
- lem_extension (L1378): cited at L488 as extension procedure.
- lem_idemp (L1916): used in the deferred proof of th_idemp_structure (L2084–2090).
- lem_carrier (L1724): used to define M in the deferred proof of th_idemp_structure.
- [ChEf77] Choi-Effros theorem: cited at L344 for exact Choi-Effros product.
- [Joh88] Johnson theorem: cited at L490 as precedent for error reduction.
- [Chr80] Christensen theorem: cited at L490 as precedent for approximate subalgebra.

---

## Open implementation questions / escalations

1. **Explicit constant in O(eps) for th_main (L461):** The paper says the constant is universal but does not state it explicitly in the theorem statement. The constant is inherited from c_0 in cor_improvement (L1317). **Escalate:** shard covering §8 (cor_improvement, L1317) must extract the explicit numerical constant c_0.

2. **Constructive nontrivial projection (lem_nontriv_projection, L931):** The proof uses the Lefschetz-Hopf fixed-point theorem on a topological space (U ⊆ A as C^1 manifold) and H-space theory. This is non-constructive. **Escalate:** shard covering §6 (L931–969) must propose a constructive finite-dim alternative (e.g., exhaustive search over approximate self-adjoint elements, or spectral methods on the approximate unitary group manifold).

3. **cb-norm computation (L347–352):** The paper references SDP computation [Watrous Chapter 3]. The implementation needs a certified SDP solver compatible with arb interval arithmetic, or a post-hoc certification step. **Escalate to implementation layer:** determine whether SDPA/SCS output can be certified with arb, or whether a direct matrix-norm formula suffices for the finite cases arising here.

4. **theta(2Phi-1) is not CP (L363):** The paper notes Phi_tilde = theta(2Phi-1) is well-defined and idempotent but not guaranteed completely positive. The eps-C* algebra built from Phi_tilde is the input to th_main, not Phi_tilde itself as a channel. **No escalation needed** for the algebra construction, but the approximate factorization result (§12) requires UCP maps — handled by th_factorization (L2730).

5. **Open problem cited (L390):** Whether all eta-idempotent UCP maps can be approximated by idempotent UCP maps with accuracy O(sqrt(eta)) or better (dimension-independently) is stated as open. **No implementation dependency** but flag for future work.

6. **Diagonal construction for eps-associative algebras (L484):** The naive Haar measure construction fails with error O(n*eps). The paper's solution (incremental construction) avoids this. However, the error reduction step (cor_improvement) still uses the diagonal of the exact C* algebra B — this is well-defined since B is exact. **Confirm in shard covering §8** that the diagonal of B is explicitly constructible (it is: B = ⊕_j M_{n_j}, so D = Σ_j Σ_{k,l} E_{kl}^{(j)} ⊗ E_{lk}^{(j)} / n_j in each block).

7. **Matrix norm vs operator norm in finite dim:** The paper uses ||·|| for the operator (spectral) norm on B(H) and ||·||_cb for the cb-norm. In finite dim these are related but distinct. The implementation must track which norm is used in each bound. In particular, ax_C* (L427) and ax_assoc (L412) use the operator norm; cb-norm bounds on Phi feed into eps-C* bounds via Section §11.3. **Confirm:** shard covering §11.3 should state how ||Phi²-Phi||_cb translates to the eps in def:eps-algebras.
