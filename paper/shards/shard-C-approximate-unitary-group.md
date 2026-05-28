# Shard C — Approximate unitary group (L690–914)

## Scope

§5 of Kitaev arXiv:2405.02434.  Input: an eps-C*-algebra A with exact unit I (exact unit assumed throughout; prop_unit from §4 / L672 provides this if needed).  Output: the set U of unitaries as a C^1 manifold, a tubular neighbourhood parametrised by a polar decomposition, and an approximate group structure (multiplication mu, inversion sigma) on U.  These objects are consumed directly by §6 (lem_nontriv_projection, prop_H-group) to extract nontrivial projections.

---

## Results

### lem_U_delta — Lemma "Approximate unitaries: norm and La_X invertibility" (statement L699–705)

- **Statement (faithful):**  Let X in U_bar_delta (i.e. ||X†X - I|| <= 2*delta and X has a right inverse).  Then:
  1. ||X|| <= 1 + O(delta + eps).
  2. La_X (left-multiplication-by-X operator on A) is invertible with ||La_X^{-1}|| <= 1 + O(delta + eps).
  3. ||La_X^{-1} - La_{X†}|| <= O(delta + eps).  (La_{X†} acts as an approximate inverse of La_X.)

- **Proof:** L707–726.  Purely norm/operator-algebra estimates; no topology.
  - ||X|| bound: C*-like inequality ||X†X|| >= (1-eps)||X||^2 plus ||X†X|| <= 1+2*delta.
  - La_X invertibility: construct explicit right inverse La_Y (La_X La_Y)^{-1} using right inverse Y of X; construct left inverse (La_{X†} La_X)^{-1} La_{X†} using ||La_{X†} La_X - 1|| <= 2*delta + O(eps) < 1.  Bounds follow by Neumann series.
  - **Constructive**: every step is an explicit matrix/operator computation.

- **Constructive in finite dim?:** yes-already.  In finite dim, La_X is the matrix of left-multiplication (an n^2 x n^2 matrix if A is n-dim); inversion is exact linear algebra.  Bounds give certified error on the inverse.

- **Depends on / Used by:**
  - Depends on: definition of eps-C*-algebra (exact unit, C*-like norm axiom), §4 operator-algebra notation (La_X, Ra_X).
  - Used by: lem_gV (L758) constructs local manifold charts using La_V^{-1}; prop_polar (L809) repeatedly invokes the bound ||La_V^{-1}|| = 1 + O(delta+eps); the inversion map sigma = u(X†) defined at L859 requires La_{U†} ~ La_U^{-1} from bound (3) above; §6 uses sigma via L939–943.

- **Arbitrary precision needed?:** Yes.  The bound O(delta+eps) is a certified numerical bound; verifying X in U_bar_delta requires computing ||X†X - I|| to guaranteed precision.  arb interval arithmetic around LAPACK/BLAS routines suffices for finite-dim matrices.

- **Cross-check / oracle:**  For exact unitaries (eps=delta=0), La_X^{-1} = La_{X^{-1}} = La_{X†} exactly.  Test: random near-unitary matrix, compare computed La_X^{-1} against La_{X†}, verify ||difference|| = O(delta+eps).

- **Implementation notes / risks:**
  - La_X is an n^2 x n^2 matrix (vectorised left-multiplication map); explicit inversion via LU is O(n^6).  For matrix algebras A = M_n(C), La_X(Z) = X*Z, so La_X^{-1}(Z) = X^{-1}*Z — reduces to one n x n linear solve, O(n^3).
  - The "right inverse" hypothesis is trivially satisfied in finite dim if La_X has nonzero determinant; no separate check needed beyond verifying injectivity.
  - FLINT/arb: compute X†X - I with arb, verify norm bound, then compute La_X^{-1} via arb matrix inversion with rigorous enclosure.

---

### lem_gV — Lemma "Local normal form: Hermitian correction g_V" (statement L758–763)

- **Statement (faithful):**  Let V in U_bar_delta.  Parametrise elements of A near V by local coordinates A = phi_V(X) = La_V^{-1}(X-V), split as A = A^parallel (anti-Hermitian) + A^perp (Hermitian).  For each anti-Hermitian A^parallel with ||A^parallel|| < 2*delta, there exists a unique Hermitian A^perp = g_V(A^parallel) solving f_V(A^parallel + A^perp) = 0, where f(X) = (1/2)(X†X - I).  The function g_V is C^1 and satisfies:
  - g_V(A^parallel) = -f(V) + O(eps*delta + delta^2),
  - partial g_V / partial A^parallel = O(delta + eps).

- **Proof:** L765–793.  Applies lem_invfun (the quantitative implicit function theorem, L564) to the map A^perp -> f_V(A^parallel + A^perp).  The key estimate (L768–771) shows:
  f_V(A) = (1/2)(V†V - I) + A^perp + O(eps*delta + delta^2),
  so partial f_V / partial A^perp = 1 + O(delta+eps), making lem_invfun applicable.  The derivative bound follows from the implicit function derivative formula (L787–791).
  - **Constructive**: lem_invfun is a Banach fixed-point / contraction-mapping argument.  In finite dim this is a Newton iteration converging quadratically.

- **Constructive in finite dim?:** yes-already.  g_V(A^parallel) is computed by Newton iteration on the Hermitian equation A^perp = A^perp - (partial f_V / partial A^perp)^{-1} f_V(A^parallel + A^perp), starting from -f(V).  One or two Newton steps suffice at eps,delta << 1 due to the O(delta+eps) derivative bound.

- **Depends on / Used by:**
  - Depends on: lem_invfun (L564), lem_U_delta (L699) for the bound on La_V^{-1}, definition of local coordinate charts phi_V (L728–756).
  - Used by: prop_polar (L809) — the entire polar decomposition proof pivots on lem_gV to parametrise unitaries near a given element; the manifold tangent space / Maurer-Cartan form (L795–807) uses g_V to compute the chart derivative; the inversion map sigma and its derivative (L889–892) implicitly use the chart structure.

- **Arbitrary precision needed?:** Yes.  Newton iteration in arb; the convergence radius depends on ||f(V)|| < delta, which must be certified.  The bound O(eps*delta + delta^2) on g_V controls the "thickness" of the manifold in the normal direction.

- **Cross-check / oracle:**  For exact V (V†V = I), g_V(A^parallel) = 0 exactly for all A^parallel, since f(V) = 0 and the manifold is tangent to i*H at I.  Test: for nearly-unitary V, verify ||g_V(0) + f(V)|| = O(eps*delta + delta^2).

- **Implementation notes / risks:**
  - The splitting A = A^parallel + A^perp is just A^parallel = (A - A†)/2, A^perp = (A + A†)/2; cheap.
  - The Newton solve for g_V uses the (n x n real-dimensional) Hermitian subspace; in finite dim the Jacobian partial f_V / partial A^perp is (dim H) x (dim H) and need not be square in coordinates unless A is a matrix algebra.  For A = M_n(C), dim A = n^2 (real 2n^2), dim H = n^2 (real), dim iH = n^2 (real skew-Hermitian).
  - Risk: the proof's O(.) constants are not tracked explicitly; an implementation must carry explicit constants through the lem_invfun contraction bound (c < 1 requires explicit threshold).

---

### prop_polar — Proposition "Approximate polar decomposition" (statement L809–811)

- **Statement (faithful):**  Define the polar map p: U x Ba^H_delta(I) -> A by p(U, H) = U*H (where H is Hermitian, Ba^H_delta(I) = {H Hermitian: ||H - I|| < delta}).  Then p is a diffeomorphism onto a set S with:
  U_{delta - O(eps*delta + delta^2)} ⊆ S ⊆ U_{delta + O(eps*delta + delta^2)},
  where U_delta = {X: ||X†X - I|| < 2*delta, X has right inverse}.

  Equivalently: every X in U_delta has a unique factorisation X = u(X) * h(X) with u(X) in U (exact unitary), h(X) Hermitian near I, and:
  - ||X - u(X)|| <= delta + O(eps*delta + delta^2).
  - Derivative of (u(X), h(X)) w.r.t. X is I + O(eps + delta) (L852–855).

- **Proof:** L813–843.  Uses the variant map (U, Q) -> U(I+Q) for Q Hermitian near 0.  By lem_gV, unitaries near V are parametrised by A^parallel; the map p_V: (A^parallel, Q) -> (B^parallel, B^perp) (local coordinate expression) has derivative = I + O(eps+delta) (L836–840).  Applies lem_invfun to conclude: (i) p_V is a local diffeomorphism; (ii) f(V) < delta_- implies ||Q|| < delta; (iii) ||Q|| < delta implies f(V) < delta + O(eps*delta + delta^2).  Constructive: contraction mapping / Newton iteration throughout.

- **Constructive in finite dim?:** yes-already.  The unitary factor u(X) is the standard polar factor: u(X) = X (X†X)^{-1/2}, computable via SVD (X = W Sigma V†, then u(X) = W V†) or by Newton iteration for (X†X)^{-1/2}.  The paper's approach via local coordinates provides certified bounds; SVD is the practical algorithm.

- **Depends on / Used by:**
  - Depends on: lem_gV (L758), lem_invfun (L564), lem_U_delta (L699), local coordinate charts (L728–756).
  - Used by:
    - L857–859: defines the group operations mu(U,V) = u(UV), sigma(U) = u(U†).  sigma is the key object for §6.
    - L861–868: certifies that UV in U_bar_{delta_0} and U† in U_bar_{delta_0}, so u is applicable.
    - L845–855: gives derivative bounds on u(X) and h(X) that feed into the approximate group law derivatives (L880–892).
    - §6 / lem_nontriv_projection (L931): the inversion map sigma and its derivative (sigma'(U) ~ -La_U Ra_U^{-1}) are used to identify fixed points of sigma as approximate projections.

- **Arbitrary precision needed?:** Yes.  SVD in arb (FLINT has arb_mat_approx_eig; may need custom Newton-Schulz iteration for certified (X†X)^{-1/2}).  The bound ||X - u(X)|| <= delta + O(eps*delta + delta^2) is the certified output spec.

- **Cross-check / oracle:**  Standard polar decomposition: for any invertible matrix X, u(X) = XV^{-1} where V = (X†X)^{1/2}.  Test: compute polar factor via SVD, verify u(X)†u(X) = I to machine precision, verify ||X - u(X)|| matches the delta bound.

- **Implementation notes / risks:**
  - SVD path: LAPACK dgesvd / zgesvd, then u(X) = W * V†.  Gives machine-precision u(X) directly without needing the local-coordinate Newton iteration.  Use arb_mat wrappers for certified bounds.
  - Newton-Schulz iteration for X(X†X)^{-1/2}: quadratic convergence, numerically stable, no square root of eigenvalues needed.  Iterate: U_{k+1} = (3 U_k - U_k U_k† U_k) / 2.  Converges when ||U_0† U_0 - I|| < 1.
  - The paper's variant map (U, Q) -> U(I+Q) uses H = I+Q near identity, so this is the "H close to I" regime.  Full polar decomposition (arbitrary Hermitian factor) follows by rescaling.
  - Degeneracy risk: if X is close to singular (||X†X - I|| near 1), the polar factor is ill-conditioned.  The proposition assumes delta < delta_max (small constant); we must enforce this precondition.
  - The group operations mu and sigma are defined L857–859.  In an implementation: mu(U,V) = polar_unitary(U*V); sigma(U) = polar_unitary(U†).  Both are O(n^3) per call.

---

### Implicit group structure (L857–893, no separate label)

- **What is constructed (L857–859):**  On U, define:
  - e = I (unit),
  - mu(U, V) = u(UV) (projected product),
  - sigma(U) = u(U†) (approximate inversion).
  These satisfy approximate group laws to O(eps) (L870–878).  Derivatives of mu and sigma in local coordinates (L880–892):
  - partial mu / partial (phi_U(X)) = La_V^{-1} Ra_V + O(eps)  [left derivative],
  - partial mu / partial (phi_V(Y)) = 1 + O(eps)  [right derivative],
  - partial sigma / partial (phi_U(X)) = -La_U Ra_U^{-1} + O(eps).

- **H-space structure (L895–912):**  The section ends by noting U is an associative H-space with left inversion map sigma.  The homotopy for associativity is constructed by projecting straight paths onto U via u (L906).  This topology is used in §6 — specifically only the existence of the left homotopy sigma(x) * x ~ e is needed (L912), not full associativity.

- **Constructive note:**  The H-space/homotopy language (L895–912) is non-constructive topology.  It is used in §6's proof of lem_nontriv_projection via Lefschetz-Hopf.  The constructive content already extracted above (sigma, mu, their derivatives) is fully algorithmic; the homotopy argument in §6 is the non-constructive step.

---

## Key definitions & notation introduced here

| Symbol | Definition | Line |
|--------|-----------|------|
| U_bar_delta | {X in A : ||X†X - I|| <= 2*delta, X has right inverse} | L693–696 |
| U_delta | same with strict inequality | L696 |
| U (calU) | exact unitaries = U_bar_0 | L692 |
| H (calH) | Hermitian elements of A | L691 |
| iH | anti-Hermitian elements of A | L691 |
| phi_V(X) | La_V^{-1}(X - V), local coordinate chart at V | L729–731 |
| phi_V^parallel(X), phi_V^perp(X) | anti-Hermitian / Hermitian parts of phi_V(X) | L743–747 |
| f(X) | (1/2)(X†X - I), "unitarity defect" map f: A -> H | L750 |
| f_V(A) | f(phi_V^{-1}(A)) = f(V(I+A)), f in local coordinates | L753–756 |
| g_V(A^parallel) | unique Hermitian solution to f_V(A^parallel + A^perp) = 0 | L759 |
| omega_U | Maurer-Cartan form: Z -> (La_U^{-1} Z)^parallel, omega_U: T_U(U) -> iH | L802–807 |
| u(X) | unitary factor in polar decomp X = u(X) h(X) | L845 |
| h(X) | Hermitian factor in polar decomp, h(X) near I | L845 |
| mu(U,V) | u(UV), approximate group multiplication on U | L857–859 |
| sigma(U) | u(U†), approximate inversion map on U | L857–859 |
| e | I, unit of U | L857–859 |
| T_U(U) | tangent space of manifold U at point U | L795 |

---

## Dependencies OUT / INTO

**Dependencies IN (what this section uses from earlier):**

- `lem_invfun` (L564–592): quantitative implicit / inverse function theorem.  Core engine for both lem_gV and prop_polar.  Every "unique solution" and "diffeomorphism" claim reduces to lem_invfun.
- `prop_unit` (L672–687): ensures A has exact unit (used in L692 "exact unit assumed").
- `rem_X2` (L628): explains why exact X^2 = I solutions may not exist, motivating the polar approach.
- eps-C*-algebra axioms: C*-like norm axiom ||(X†X)|| >= (1-eps)||X||^2 (used in lem_U_delta proof L708), approximate associativity bounds on La_X La_Y - La_{XY} (L655–660).
- Left/right multiplication operators La_X, Ra_X and their norm bounds (L632–670).
- Frechet derivative bounds partial_X(X^2-X) = La_X + Ra_X - 1 (L663–670).

**Dependencies OUT (what §6 and later sections consume from here):**

- `sigma(U) = u(U†)`: the inversion map.  §6 (lem_nontriv_projection, L931–943) characterises fixed points of sigma as O(eps)-projections via P = (1/4)(2I + U + U†) when sigma(U) = U.  The derivative partial sigma / partial phi_U = -La_U Ra_U^{-1} + O(eps) at L892 is needed to show ±I are isolated fixed points.
- `mu(U,V) = u(UV)`: used in prop_H-group (L971) for the H-space structure argument.
- Manifold structure (charts phi_V, tangent space T_U, Maurer-Cartan form omega_U): used in §6's orientability argument (L950) and the quotient manifold construction (L945).
- The bound ||X - u(X)|| <= delta + O(eps*delta + delta^2) from prop_polar: used wherever approximate unitaries must be retracted to exact unitaries.
- U_delta tubular neighbourhood: §6 works inside U_{delta_max} throughout.

---

## Open implementation questions / escalations

1. **Explicit O(.) constants.**  The paper uses O(eps), O(delta+eps), O(eps*delta+delta^2) throughout without naming the constants.  An arb implementation must fix these constants to certify correctness.  Needs: a pass through each proof extracting explicit numeric prefactors (e.g. the "2" in ||X†X - I|| <= 2*delta, the factor in lem_invfun's c < 1 condition).

2. **Newton-Schulz vs. SVD for u(X).**  SVD is cleaner and gives exact polar factor; Newton-Schulz is iterative and more amenable to certified arb bounds.  Decision: use SVD for practical speed, add a certified Newton refinement step for arb output.  Verify both give ||u(X)†u(X) - I|| < machine_eps.

3. **La_X as n^2 x n^2 matrix vs. n x n linear solve.**  For A = M_n(C), La_X(Z) = XZ, so La_X^{-1}(Z) = X^{-1}Z.  The shard should explicitly state: for matrix algebras reduce all La_X^{-1} calls to n x n linear solves, avoiding the n^2 x n^2 representation.  General eps-C*-algebras require the full n^2 x n^2 treatment.

4. **delta_max / eps_max thresholds.**  Proposition prop_polar requires eps, delta < some absolute constant delta_max (L697).  The paper does not state delta_max explicitly.  Must either: (a) extract the condition from lem_invfun (c < 1 in the contraction, tracing through all inequalities), or (b) treat delta_max as a parameter verified numerically.  Escalate: tracked in implementation as a named constant to be calibrated.

5. **H-space / Lefschetz-Hopf in §6.**  The constructive content of §5 (sigma, mu, polar decomp) is clean.  But §6's use of the Lefschetz-Hopf fixed-point theorem (L957–960) to guarantee existence of a nontrivial fixed point of sigma is non-constructive.  The constructive route for §6 is: enumerate fixed points of sigma by numerical root-finding (Newton's method on sigma(U) - U = 0), then certify each with arb.  This is flagged as an escalation for shard D (§6).

6. **Right-inverse hypothesis in U_bar_delta.**  In finite dim, "X has a right inverse" iff det(X) ≠ 0.  The implementation must check this (or bound away from zero det) before asserting X in U_bar_delta.

7. **Approximate group laws: bound tracking.**  The group laws (L870–878) hold to O(eps).  An algorithm computing in U will accumulate O(eps) errors per group operation.  Need error budget: how many mu/sigma calls before accumulated error exceeds delta_max?  This is application-dependent but must be documented.
