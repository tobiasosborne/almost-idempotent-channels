# Shard E — delta-homomorphisms & main-theorem proof (L1190–1446)

## Scope

Lines 1190–1446 of `approximate_algebras.tex`. Covers §8 "Properties and approximation
of delta-homomorphisms" and §9 "Proof of the main theorem" (label `sec_proof_main`).
These two sections contain the central constructive algorithm of the whole paper: given
a finite-dimensional eps-C*-algebra A, build a genuine C*-algebra B and an O(eps)-
isomorphism v: B -> A by an incremental loop that is stabilised at each step by
error-reduction via the diagonal/cohomology trick.

~30-line context reads: §7 (L1080–1188) defines S_{P,Q}, the compressed product, the
Ha-map (L1146–1160), and one-dimensional projection lemmas (lem_PQR L1162, lem_1d_proj
L1179, equivalence of 1d-projections L1187). These are hard dependencies for §8–9.
th_main is stated at L460; its deferred proof begins at L1414.

---

## The master loop (th_main proof)

The proof at L1414 constructs a c_0*eps-isomorphism v: B -> A in three explicit stages.
c_0 is the explicit constant from cor_improvement (L1317). The loop is:

**Pre-loop setup.**
Fix eps_max, delta_max, c_0 from cor_improvement. The target is the eps-C*-algebra A
with eps < eps_max. Set eps' = O(eps) to be a uniform upper bound valid inside all
compressed algebras S_P encountered below.

**Stage 1 — build the maximal commutative skeleton.**

1. [Find commutative B_comm.] Among all c_0*eps-inclusions of commutative C*-algebras
   into A, pick one v_comm: B_comm -> A of maximum dimension.
   Power: lem_nontriv_projection (L931); cor_improvement (L1317).
   Error maintained: v_comm is a c_0*eps-inclusion by construction.

2. [Force 1d projections.] Let {Pi_1,...,Pi_m} be the projection basis of B_comm;
   set P_j = v_comm(Pi_j). Claim: each P_j is one-dimensional (dim S_{P_j} = 1).
   Proof by contradiction (L1417–1426): if some P_m has dim S_{P_m} > 1, apply
   lem_nontriv_projection to S_{P_m} to split P_m into P', P''. Merge the resulting
   O(eps)-inclusions via cor_merge_sum (L1352) to get an O(eps)-inclusion of a
   strictly larger commutative algebra; then cor_improvement yields a c_0*eps-inclusion,
   contradicting maximality. Hence all P_j are one-dimensional O(eps)-projections.

3. [Classify equivalence classes.] Partition {1,...,m} into equivalence classes by:
   j ~ k  iff  dim S_{P_j,P_k} = 1  (L1428, using lem_1d_proj L1179 and lem_PQR L1162
   for transitivity). Denote class C, projection P_C = sum_{j in C} P_j, compressed
   algebra A_C = S_{P_C}. Each P_C is a c_0*eps-projection and A_C is an eps'-C*-algebra
   (L1428).

**Stage 2 — extend each block to a full matrix algebra.**

For each equivalence class C = {1,...,s} (treated independently, L1430):

4. [Base case r=1.] The one-dimensional projection P_1 gives A_1 = S_{P_1} ~ CC.
   Set v_1: M_1 -> A_1 to be the obvious c_0*eps'-isomorphism (scalar multiples of
   the unit tilde{P_1}).

5. [Inductive step r-1 -> r, for r = 2,...,s.] (L1435–1441)
   a. Compress v_{r-1}: M_{r-1} -> A_{r-1} into A_r = S_{P_{[1,r]}} via the
      compression map Co_{P_{[1,r]}}, obtaining v'_{r-1}: M_{r-1} -> A'_{r-1}
      and compressed projections P'_{[1,r-1]}, P'_r.
   b. Check dim S_{P'_{[1,r-1]}, P'_r} != 0 (guaranteed because j ~ r for j in C,
      i.e. dim S_{P_j, P_r} = 1 for each j, and lem_add_dim (L1363) gives
      dim S_{P,Q} = sum_{j,k} dim S_{P_j,Q_k} = r-1 != 0 after compression).
   c. Apply lem_extension (L1378): extend the O(eps)-isomorphism v'_{r-1}: M_{r-1} ->
      A'_{r-1} to an O(eps)-isomorphism v^+_{r-1}: M_r -> A_r.
      Power: lem_approx (L1224) inside the proof of lem_extension approximates the
      O(eps)-homomorphism h_11 v by a genuine homomorphism mu_11; lem_merging (L1325)
      assembles the block maps gamma_{jk}.
   d. Apply cor_improvement (L1317): replace v^+_{r-1} with a c_0*eps'-isomorphism v_r.
      O(eps') error maintained after this step.

   After s steps, have c_0*eps'-isomorphism v_C: M_s -> A_C.

**Stage 3 — merge blocks.**

6. [Merge over equivalence classes.] (L1443)
   The off-diagonal spaces S_{P_C, P_D} = 0 for C != D (lem_add_dim).
   Successively merge the isomorphisms v_C via cor_merge_sum (L1352), then apply
   cor_improvement after each merge to keep error at c_0*eps'.
   Final result: a c_0*eps-isomorphism v: B -> A where B = direct_sum_C M_{|C|}.

**O(eps) invariant.** The constant c_0 is fixed once at the start (from cor_improvement).
Every use of cor_merge_sum introduces at most O(eps') error; cor_improvement immediately
reduces this back to c_0*eps'. Since eps' = O(eps) and c_0 is universal (independent of
dim A), the final isomorphism satisfies ||v(X) - ...|| <= c_0*eps * ||X|| with an
explicit, dimension-independent constant. This is the claim of th_main (L460).

---

## Results

### prop_delta_hominc — Proposition "Norm bounds for delta-homomorphisms" (statement L1194–1195)

- **Statement (faithful):** Let v: A' -> A'' be a non-unital delta-homomorphism of
  eps-C*-algebras. Then (i) ||v|| <= 1 + O(delta + eps); (ii) if ||v(X)|| >= eta*||X||
  for some eta > 2*delta and all X, then ||v(X)|| >= (1 - O(delta+eps))*||X|| for all X;
  (iii) if ||v(I_{A'}) - I_{A''}|| < some absolute positive constant, then
  ||v(I_{A'}) - I_{A''}|| <= O(delta + eps).

- **Proof:** L1198–1222. Technique: bootstrap fixed-point inequalities.
  From the delta-homomorphism condition ||v(X*X) - v(X)v(X)*|| <= delta*||X||^2 and the
  C*-property, bound the inf a and sup b of ||v(X)||/||X|| by quadratic inequalities in
  a and b (L1210–1215). The unit bound uses lem_invfun (L564) applied to x^2 - x near
  I (L1217–1221).
  Technique label: bootstrap / implicit-function.

- **Constructive in finite dim?:** Yes. The bounds are explicit polynomial inequalities
  in (a, b, delta, eps). No oracle needed; just evaluate ||v||, check the hypotheses.

- **Depends on:** lem_invfun (L564), basic C*-property and eps-C*-algebra definition.
- **Used by:** lem_approx (L1224, cited at L1260, L1311), cor_improvement (L1317),
  lem_merging (L1348), lem_extension (L1391).

- **Arbitrary precision needed?:** No. All bounds are closed-form in delta, eps.

- **Cross-check / oracle:** Set delta=0, eps=0: recover ||v|| = 1 and isometry, i.e.,
  exact *-homomorphism is an isometry on the unit ball.

- **Implementation notes / risks:** The three conclusions require three separate
  preconditions; implement as three distinct predicates. The "certain positive constant"
  for conclusion (iii) is not named explicitly; tie it to the implicit function theorem
  constant from lem_invfun (L564/579), which itself should be certified.

---

### lem_approx — Lemma "Approximation of delta-homomorphisms" (statement L1224–1226)

- **Statement (faithful):** For any delta-homomorphism v from a finite-dimensional C*-
  algebra B to an eps-C*-algebra A, there exists an O(eps)-homomorphism tilde_v: B -> A
  such that ||tilde_v(X) - v(X)|| <= O(delta)*||X|| for all X in B.

- **Proof:** L1256–1314. Technique: **multiplicativity defect** + **diagonal** +
  Newton iteration.
  (1) Define multiplicativity defect g = G_v: B x B -> A by g(X,Y) = v(XY) - v(X)v(Y);
      ||g(X,Y)|| <= delta*||X||*||Y|| (L1258–1259).
  (2) Show g satisfies an approximate 2-cocycle equation up to O(eps)*||X||*||Y||*||Z||
      error (L1261–1266), using associativity of B, eps-associativity of A, and
      ||v|| <= const from prop_delta_hominc.
  (3) For correction w with ||w|| <= O(delta), compute G_{v+w}(X,Y) = g(X,Y) -
      F_w(X,Y) + O(delta^2)*||X||*||Y|| where F_w(X,Y) = v(X)w(Y) - w(XY) + w(X)v(Y)
      (L1268–1274).
  (4) Choose w'(X) = sum_j v(A_j) g(B_j, X) using the diagonal D = sum_j A_j ⊗ B_j
      of B (finite-dimensional formula L1246–1247). Use the diagonal identity
      sum_j v(XA_j)g(B_j,Y) = sum_j v(A_j)g(B_jX,Y) (from XD=DX) and the 2-cocycle
      equation to show F_{w'}(X,Y) = g(X,Y) + O(delta^2 + eps)*||X||*||Y|| (L1282–1303).
  (5) Symmetrise: set v^(1) = v + (1/2)(w' + w'') where w''(X) = w'(X*)* to restore
      involution-compatibility (L1305–1311). Then G_{v^(1)}(X,Y) <= O(delta^2+eps)*...
  (6) Iterate as Newton's method: delta_s ~ delta^{2^s} until delta_s <= O(eps);
      stop after finitely many steps. If eps=0, take limit (L1312–1313).

- **Constructive in finite dim?:** YES, fully constructive. The diagonal D for a
  finite-dim C*-algebra B = ⊕_l M_{d_l} is given explicitly by the generalised Pauli
  formula (L1250–1254):
    D = ⊕_l sum_{j,k=0}^{d_l-1} d_l^{-2} S_{jk}^* ⊗ S_{jk},
  with S_{jk} the shift-clock (Heisenberg-Weyl) operators. For a direct sum, the
  component diagonals are tensored (L1254). The correction w'(X) involves evaluating v
  on a fixed finite list of unitaries {A_j} and computing g(B_j, X) = v(B_j X) - v(B_j)v(X).
  Number of Newton steps S = O(log(delta/eps)) (with exact arithmetic) or S = O(1)
  with arb/FLINT interval arithmetic since eps > 0.

- **Depends on:** prop_delta_hominc (L1194, for ||v|| bound); diagonal of finite-dim
  C*-algebra (L1245–1254); Johnson's theorem [Joh88] cited for eps=0 limit.
- **Used by:** cor_improvement (L1317); lem_extension (L1382, cited at L1391).

- **Arbitrary precision needed?:** Yes for the Newton iteration: need to certify that
  delta_s <= O(eps) to stop. With arb ball arithmetic, each step doubles the "good bits"
  minus a log(n) overhead (but n = dim B is fixed, not the ambient dim).

- **Cross-check / oracle:** eps=0: tilde_v is an exact homomorphism of C*-algebras
  (Johnson 1988 Thm 3.1). For the Pauli diagonal: verify XD = DX and pi(D) = I
  symbolically on the Pauli basis.

- **Implementation notes / risks:**
  - Pauli diagonal has d_l^2 terms per block; for d_l large this is expensive but
    exact. Alternative: Caratheodory representation (L1245) with fewer terms is
    less canonical but reduces cost; need a separate convex-hull algorithm.
  - The correction w' is a linear map on B; represent as a dim(B) x dim(B) matrix
    acting on vec(X). Newton steps are matrix updates.
  - The O() in "O(delta^2 + eps)" hides a concrete constant that must be tracked for
    the termination guarantee. Risk: without explicit constants, the iteration may not
    certifiably terminate.
  - FLINT/arb: represent elements of A as complex matrices; products and norms via arb.

---

### cor_improvement — Corollary "Error reduction" (L1317–1319)

- **Statement (faithful):** There exist explicit positive constants eps_max, delta_max,
  c_0 such that: for all eps < eps_max, if a finite-dimensional C*-algebra B is
  delta_max-included into an eps-C*-algebra A, then there also exists a c_0*eps-inclusion
  of B into A. If the original inclusion is bijective (i.e. a delta_max-isomorphism),
  the new inclusion is also bijective.

- **Proof:** L1317 (one sentence). Immediate combination of lem_approx (L1224) and
  prop_delta_hominc (L1194). lem_approx gives tilde_v with multiplicativity defect
  O(eps); prop_delta_hominc (conclusion iii) then shows tilde_v is an O(eps)-inclusion.
  Bijectivity: lem_approx preserves the lower norm bound from the original inclusion
  (prop_delta_hominc conclusion ii), so tilde_v remains injective.
  Technique label: **cohomology-diagonal** (the diagonal of B used inside lem_approx
  is the same object that trivialises Hochschild cohomology; L480–482 gives the
  conceptual explanation; the finite-dim version is constructive).

- **Constructive in finite dim?:** YES. The computation is identical to lem_approx:
  build D from Pauli operators of B, run Newton iteration to reduce delta to O(eps),
  certify the norm bounds with arb. The constants eps_max, delta_max, c_0 are in
  principle extractable from the proof of lem_approx; the paper says they are explicit
  but does not state them numerically—this is the main open question (see below).

- **Depends on:** lem_approx (L1224); prop_delta_hominc (L1194).
- **Used by:** proof of th_main at L1426, L1441, L1443 (every place in Stages 1–3
  where error is reset).

- **Arbitrary precision needed?:** Yes. c_0 must be a certified upper bound; use arb
  interval throughout the Newton loop of lem_approx.

- **Cross-check / oracle:** Apply to the identity map v = id: B -> B (which is a
  0-inclusion). Output should still be a 0-inclusion. c_0*eps should be independent of
  dim B; verify numerically for small examples.

- **Implementation notes / risks:**
  - c_0 is the pivotal constant for the entire algorithm; a loose estimate makes the
    algorithm vacuous for large eps. Extracting c_0 explicitly from the Newton-step
    analysis of lem_approx is essential prior work.
  - Risk: the paper notes (L484) that "naive constructions of the diagonal in the
    eps-associative setting have error bounds proportional to n = dim A." For
    cor_improvement this is not a problem because B is a GENUINE C*-algebra and its
    Pauli diagonal is exact. Document this clearly in the implementation.

---

### lem_merging — Lemma "Merging block maps" (statement L1325–1346, proof L1348–1350)

- **Statement (faithful):** Let Pi_1, Pi_2 be complementary projections in a C*-algebra
  B (Pi_1 + Pi_2 = I), and P_1, P_2 delta-projections in an eps-C*-algebra A with
  ||P_1 + P_2 - I|| <= delta. Given linear maps gamma_{jk}: S_{Pi_j,Pi_k} -> S_{P_j,P_k}
  (j,k in {1,2}) satisfying:
    (merging0) gamma_{kj}(X*) = gamma_{jk}(X)*,
    (merging1) ||gamma_{jl}(XY) - gamma_{jk}(X) . gamma_{kl}(Y)|| <= delta*||X||*||Y||,
    (merging2) ||gamma_{jj}(Pi_j) - P_j|| <= delta,
    (merging3) (1-delta)||X|| <= ||gamma_{jk}(X)|| <= (1+delta)||X||,
  the combined map gamma: (X_{jk}) |-> sum_{j,k} gamma_{jk}(X_{jk}): B -> A is an
  O(delta+eps)-inclusion. Bijectivity of all gamma_{jk} implies bijectivity of gamma.

- **Proof:** L1348–1350. Uses canonical bijection mu: B -> ⊕_{j,k} S_{Pi_j,Pi_k}
  (maximum norm), and the map alpha from lem_alpha (L1086). Conditions (merging0–2)
  give an O(delta)-homomorphism; (merging3) + lem_alpha give norm equivalence; then
  prop_delta_hominc upgrades to O(delta+eps)-inclusion.
  Technique label: block-matrix assembly / lem_alpha norm estimate.

- **Constructive in finite dim?:** Yes. All maps are explicitly specified; verification
  of conditions (merging0–3) is a finite computation. The assembly gamma is just
  summing the block maps.

- **Depends on:** lem_alpha (L1086); prop_delta_hominc (L1194); def of S_{P,Q} and
  compressed product (L1080).
- **Used by:** cor_merge_sum (L1352); lem_extension (L1382, L1411).

- **Arbitrary precision needed?:** Only to certify the O() bound in O(delta+eps); same
  level as prop_delta_hominc.

- **Cross-check / oracle:** Take B = M_2, P_1 = E_{11}, P_2 = E_{22} in exact C*-algebra.
  Then delta=0=eps; gamma should be an exact isomorphism.

- **Implementation notes / risks:** The "dot" in (merging1) is the compressed product
  (L1080), not the ambient product in A. Must implement Co_{P,R}(XY) correctly.

---

### cor_merge_sum — Corollary "Merging direct sums" (L1352–1358)

- **Statement (faithful):** Let B_1, B_2 be C*-algebras, P_1, P_2 delta-projections in
  an eps-C*-algebra A with ||P_1+P_2-I|| <= delta, and v_j: B_j -> S_{P_j} delta-
  inclusions (j=1,2). The combined map v: (X_1,X_2) |-> v_1(X_1)+v_2(X_2): B_1 ⊕ B_2
  -> A is an O(delta+eps)-inclusion. If v_1, v_2 are bijective and S_{P_1,P_2}=0, then
  v is bijective.

- **Proof:** L1352–1358 (one paragraph). Immediate from lem_merging: set gamma_{jj} = v_j
  and gamma_{12}=gamma_{21}=0. Bijectivity of gamma requires S_{P_1,P_2}=0 to ensure
  the off-diagonal blocks vanish.
  Technique label: direct application of lem_merging.

- **Constructive in finite dim?:** Yes. Just concatenate the two delta-inclusion maps.

- **Depends on:** lem_merging (L1325).
- **Used by:** proof of th_main Stage 1 (L1426) and Stage 3 (L1443).

- **Arbitrary precision needed?:** No beyond what lem_merging requires.

- **Cross-check / oracle:** v should recover exact direct-sum isomorphism when delta=eps=0.

- **Implementation notes / risks:** The condition S_{P_1,P_2}=0 must be certified
  computationally (dim S_{P_1,P_2} = 0). This requires checking that ||Co_{P_1,P_2}(X)||
  is small for all X in a basis; may need a threshold and arb certification.

---

### lem_add_dim — Lemma "Dimension decomposition of S_{P,Q}" (statement L1363–1365, proof L1367–1369)

- **Statement (faithful):** Let B and C be commutative C*-algebras with projection bases
  {Pi_1,...,Pi_p} and {Sigma_1,...,Sigma_p}, and let v: B -> A, w: C -> A be non-unital
  delta-inclusions into an eps-C*-algebra A. Set P_j = v(Pi_j), Q_k = w(Sigma_k),
  P = v(I_B), Q = w(I_C). Then there is a bounded linear bijection
    S_{P,Q}  ~=  ⊕_{j,k} S_{P_j,Q_k},
  in particular dim S_{P,Q} = sum_{j,k} dim S_{P_j,Q_k}.

- **Proof:** L1367–1369. Inductive application of lem_alpha (L1086) to partial sums
  P_{[1,j]} and Q_{[1,k]}; near-orthogonality of the P_j (resp. Q_k) from the
  delta-inclusion conditions.
  Technique label: iterated lem_alpha / dimension counting.

- **Constructive in finite dim?:** Yes. The bijection is the map alpha from lem_alpha
  applied iteratively; it is explicit. dim S_{P,Q} is computable as rank of the
  compression Co_{P,Q}.

- **Depends on:** lem_alpha (L1086); non-unital delta-inclusion definition (L451).
- **Used by:** lem_extension (L1382, L1383) to show dim S_{P,Q} = n; proof of th_main
  Stage 2 (L1435) to ensure the extension hypothesis is met.

- **Arbitrary precision needed?:** The bijection norm bound involves O(pq*(delta+eps));
  since p,q <= 2 throughout actual use, constants are harmless.

- **Cross-check / oracle:** For exact orthogonal projections P_j, Q_k in M_n: recover
  dim S_{P,Q} = sum_{j,k} dim S_{P_j,Q_k} exactly.

- **Implementation notes / risks:** The "bounded linear bijection" is only O(1)-bi-
  Lipschitz (not isometric); the implementation must track the Lipschitz constants for
  downstream norm estimates.

---

### lem_extension — Lemma "Extension" (statement L1378–1380, proof L1382–1412)

- **Statement (faithful):** Let A be an eps-C*-algebra, P, Q delta-projections with
  ||P+Q-I|| <= delta. Suppose v: M_n -> S_P is a delta-isomorphism, dim S_Q = 1, and
  S_{P,Q} != 0. Then v extends to an O(delta+eps)-isomorphism v_+: M_{n+1} -> A.

- **Proof:** L1382–1412. Multi-step construction:
  (1) Show dim S_{P,Q} = n. The images P_j = v(E_{jj}) are equivalent 1d O(delta+eps)-
      projections; each dim S_{P_j,Q} in {0,1}; S_{P,Q} != 0 forces all = 1; lem_add_dim
      gives dim S_{P,Q} = n (L1383).
  (2) Define Ha-maps h_{jk} = Ha^Q_{P_j,P_k} (L1389). h_{11}: S_P -> B(S_{P,Q}) is an
      O(delta+eps)-homomorphism; h_{12}, h_{21}, h_{22} are O(delta+eps)-close to identity.
  (3) Apply lem_approx to approximate h_{11}*v: M_n -> B(S_{P,Q}) by an exact homomorphism
      mu_{11} with O(delta+eps) accuracy. Since M_n is simple, mu_{11} is an isomorphism
      (L1391). Therefore h_{11} (and by construction each h_{jk}) is an O(delta+eps)-
      isomorphism.
  (4) Write mu_{11}(A) = U_1 A U_1* for unitary U_1: C^n -> S_{P,Q} (L1404).
  (5) Define gamma_{11}=v, gamma_{12}(A)=U_1(A), gamma_{21}(A)=(U_1(A*))*, gamma_{22}(A)=A*tilde_Q.
      Check that h_{jk}*gamma_{jk} is O(delta+eps)-close to mu_{jk} (L1411), hence
      (merging0–3) hold with delta' = O(delta+eps).
  (6) Apply lem_merging (L1325) to get v_+ (L1411).
  Technique label: Ha-map linearisation + lem_approx inside + lem_merging assembly.

- **Constructive in finite dim?:** YES, fully constructive.
  - dim S_{P,Q} is computed as a rank (step 1).
  - h_{11} is the explicit map Z |-> Ha^Q_{P,P}(Z) which is computable from inner
    products on S_{P,Q} (def at L1146–1149).
  - lem_approx is constructive (Newton + Pauli diagonal).
  - U_1 is found by computing the unitary in the polar decomposition of the isomorphism
    mu_{11} written in a basis of S_{P,Q}; this is a standard SVD (LAPACK: zgesvd).
  - The assembly of gamma_{jk} and v_+ is explicit matrix arithmetic.

- **Depends on:** lem_add_dim (L1363); lem_PQ_Hilb (L1123); Ha-map construction
  (L1146–1160); lem_approx (L1224); lem_merging (L1325); lem_1d_proj (L1179);
  lem_PQR (L1162).
- **Used by:** proof of th_main Stage 2, step (c) (L1441).

- **Arbitrary precision needed?:** Yes at the lem_approx step inside; same caveat as
  lem_approx.

- **Cross-check / oracle:** Start with v: M_1 -> S_P (scalar), Q 1d, S_{P,Q} != 0.
  Output should be v_+: M_2 -> A. If A = M_2 exactly, v_+ should be an exact isomorphism.

- **Implementation notes / risks:**
  - The Ha-map requires Euclidean inner products on S_{P,Q} (lem_PQ_Hilb); these
    themselves have O(delta+eps) error in the norm equivalence. Track carefully.
  - LAPACK zgesvd for the polar decomposition step (U_1). Condition number of mu_{11}
    is 1 + O(delta+eps) so this is well-conditioned.
  - The step "h_{11} is O(delta+eps)-close to mu_{11}*v^{-1}" requires v^{-1}; since v
    is a delta-isomorphism, v^{-1} is a (delta+O(delta^2))-isomorphism (L458); use arb
    to certify the inverse.

---

### Proof of th_main — (statement L460–461, proof L1414–1444)

- **Statement (L460–461):** Any finite-dimensional eps-C*-algebra A is O(eps)-isomorphic
  to some C*-algebra B. The implicit constant in O(eps) is independent of A and dim A.

- **Proof structure:** See "The master loop" section above. Three explicit stages
  (L1415–1443). The key tools at each stage are:
    Stage 1: lem_nontriv_projection (L931), cor_merge_sum (L1352), cor_improvement (L1317).
    Stage 2: lem_extension (L1378), cor_improvement (L1317).
    Stage 3: cor_merge_sum (L1352), cor_improvement (L1317), lem_add_dim (L1363).
  Technique label: incremental C*-algebra construction / Lefschetz + cohomology-diagonal.

- **Constructive in finite dim?:** Partially. All three stages are constructive EXCEPT:
  - lem_nontriv_projection (Stage 1) relies on the Lefschetz-Hopf theorem via a
    topological argument (L934–968). The projection P itself is O(eps)-computable once
    we know it exists (by searching fixed points of sigma on the approximate unitary group
    calU). But the existence proof is non-constructive. In finite dim with eps
    small enough, one can find P by exhaustive search or eigenvalue methods on the
    approximate unitaries — this is an ESCALATION (see Open questions below).

- **Depends on:** All results in this shard; also lem_nontriv_projection (L931),
  prop_H-group (L971), lem_alpha (L1086), lem_PQ_Hilb (L1123), S_{P,Q} machinery.
- **Used by:** th_idemp_structure (L318, proof at L2055); th_main_ext (L1538).

- **Arbitrary precision needed?:** Yes throughout; see individual lemmas.

- **Cross-check / oracle:** B must be a direct sum of matrix algebras (th_main gives no
  more). Certify: (i) v is a c_0*eps-isomorphism by checking all four conditions of
  def delta-inclusion (L447–455); (ii) B = ⊕ M_{|C|} has the correct block structure.

---

## Key definitions & notation introduced here

**Multiplicativity defect** (L1256–1259): For a linear map u: B -> A,
  G_u(X,Y) = u(XY) - u(X)u(Y): B x B -> A.
A map is a delta-homomorphism iff ||G_u(X,Y)|| <= delta*||X||*||Y||.

**Diagonal D of a finite-dim C*-algebra** (L1239–1254):
  D = sum_j A_j ⊗ B_j in B ⊗_proj B satisfying XD = DX for all X and pi(D) = I_B
  where pi(sum A_j ⊗ B_j) = sum A_j B_j.
  For a genuine finite-dim C*-algebra, D = integral dU (U* ⊗ U) over Haar measure;
  by Caratheodory, this is a finite convex combination D = sum_j p_j U_j* ⊗ U_j with
  sum p_j = 1, p_j >= 0, U_j unitary, and sum p_j ||U_j*|| ||U_j|| = 1 (L1246–1247).
  For M_d, explicit formula via generalised Pauli operators S_{jk} (L1250–1253).
  The unique norm-1 diagonal; not to be confused with "diagonal subalgebra".

**F_w** (L1273–1274): Linearised change in multiplicativity defect,
  F_w(X,Y) = v(X)w(Y) - w(XY) + w(X)v(Y).
  This is the Hochschild coboundary of w viewed as a 1-cochain.

**Approximate 2-cocycle equation** (L1261–1266): The multiplicativity defect g of a
delta-homomorphism satisfies v(X)g(Y,Z) - g(XY,Z) + g(X,YZ) - g(X,Y)v(Z) = O(eps).

**Compressed product (dot product)** (L1080): X . Y = Co_{P,R}(XY) for X in S_{P,Q},
Y in S_{Q,R}. Used in lem_merging condition (merging1).

**Ha-map** (L1146–1149): For one-dim delta-projection Q, and delta-projections P, R,
  Ha^Q_{P,R}: S_{P,R} -> B(S_{R,Q}, S_{P,Q})
defined by the symmetrised inner-product equation (L1148). Satisfies Ha_{R,P}(Z*)=Ha_{P,R}(Z)*
(L1152) and is O(delta+eps)-close to left-multiplication (L1155).

**Equivalence of 1d projections** (L1187): P ~ Q (both 1d delta-projections) iff
dim S_{P,Q} = 1. Equivalence relation by lem_PQR (L1162).

**Projection basis** (L1361): For commutative C*-algebra B, a set {Pi_1,...,Pi_n} with
Pi_j* = Pi_j, Pi_j Pi_k = delta_{jk} Pi_k, sum Pi_j = I.

**Compression map** (L1421): Co_P(X) = P X + X P - P X P (approximately); more precisely
Co_{P,Q}(X) defined at L1080. Used to "fit" maps into compressed subalgebras.

---

## Dependencies OUT / INTO

**Dependencies in (needed by this shard):**
- lem_alpha (L1086): bijection S_P -> ⊕ S_{P_j,Q_k}, norm bounds.
- lem_PQ_Hilb (L1123): Hilbert space structure on S_{P,Q} when Q is 1d.
- Ha-map definition and properties (L1146–1160).
- lem_PQR (L1162): near-isometry of compressed product through 1d projection.
- lem_1d_proj (L1179): dim S_{P,Q} <= 1 when both P, Q are 1d.
- lem_nontriv_projection (L931): existence of nontrivial O(eps)-projection (used in
  Stage 1, non-constructive via Lefschetz).
- lem_invfun (L564): implicit function theorem used in prop_delta_hominc.

**Dependencies out (this shard used by):**
- th_idemp_structure (L318, proof L2055): uses th_main as a black box.
- th_main_ext (L1538): tensor extension of th_main result; uses the v: B->A from th_main.
- prop_inc_ext (L1483): uses cor_improvement to extend inclusions.
- lem_approx_ext (L1508): uses lem_approx for operator-space maps.

---

## Open implementation questions / escalations

1. **ESCALATION: constructive nontrivial projection.** lem_nontriv_projection (L931)
   is proved via Lefschetz-Hopf, a purely existential argument. In Stage 1 we need an
   actual matrix P. Proposed constructive substitute: enumerate approximate unitaries U
   in calU satisfying ||U* - U|| <= O(eps), compute P = (2I + U + U*)/4, keep those
   that are nontrivial. Requires: (a) a constructive description of calU (already built
   in §5–6); (b) either exhaustive search (exponential) or eigenvalue methods on the
   approximate unitary group manifold. This is a non-trivial algorithmic step; requires
   separate shard or explicit algorithm design.

2. **Explicit constants eps_max, delta_max, c_0.** The proof of cor_improvement says
   "there exist" but gives no numerical values. These constants must be extracted from
   the Newton convergence analysis of lem_approx (quadratic convergence kicks in once
   delta < sqrt(eps); the concrete threshold depends on the Lipschitz constant of F_w,
   which depends on ||v|| <= 1 + O(delta+eps)). This extraction is essential before
   any implementation can certify correctness.

3. **Pauli diagonal vs. Caratheodory representation.** For d_l large, the Pauli diagonal
   has d_l^2 terms per block; Caratheodory gives at most dim(B)^2 + 1 terms total. For
   implementation, the Pauli formula (L1250–1253) is preferred (explicit, no convex-hull
   computation), but costs O(d_l^4) flops per Newton step. For the matrix algebras M_{|C|}
   encountered in Stage 2, |C| = s <= dim A, so d_l = s can be large. Profile this.

4. **Norm of v^{-1} in lem_extension.** The proof uses v^{-1}: S_P -> M_n, which is a
   (delta+O(delta^2))-isomorphism (L458). Tracking the exact inverse norm through the
   compression step (Co_{P_{[1,r]}} applied to v_{r-1}) requires careful bookkeeping;
   each compression introduces an O(eps') additive error. Across s steps this could
   accumulate. Verify that the O() in O(delta+eps) after lem_extension absorbs this.

5. **S_{P_C, P_D} = 0 certification in Stage 3.** The proof asserts this by lem_add_dim
   and the definition of equivalence classes. Computationally: dim S_{P_C, P_D} = sum_{j
   in C, k in D} dim S_{P_j, P_k} = 0 because j !~ k. Each dim S_{P_j,P_k} is a rank
   of Co_{P_j,P_k} on a finite-dim space; certify this rank is 0 up to a numerical
   threshold (arb lower bound on singular values).

6. **Order of merges in Stage 3.** The proof says "successively merge the isomorphisms
   v_C for different C." The order matters for intermediate error bounds; each merge
   introduces O(eps') error then cor_improvement brings it back to c_0*eps'. Confirm
   that after at most m-1 merges (m = number of equivalence classes), the total error
   is still O(eps) with the same c_0. Likely yes since each reduction uses the same c_0,
   but verify that c_0*c_0 = c_0*O(1) is absorbed by the O() in O(eps).
