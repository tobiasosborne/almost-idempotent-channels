# Shard D — Nontrivial projection & subspaces (L915–1189)

## Scope

§6 "The existence of a nontrivial projection" (L915–969) and §7 "Subspaces associated with
projections" (L1052–1188).  The section boundary between §6 and §7 is at L1052; the prop_H-group
proof (L1023–1050) straddles it but belongs to §6 logically.

---

## Results

### lem_nontriv_projection — Lemma "nontrivial O(ε)-projection exists" (statement L931–932)

- **Statement (faithful):**
  Any ε-C* algebra A with 1 < dim A < ∞ has a nontrivial O(ε)-projection P.
  "Nontrivial" means both P and I−P are nonvanishing, i.e. ‖P‖ = 1+O(ε) and ‖I−P‖ = 1+O(ε)
  (equivalently, ‖P‖ is bounded away from 0 and from 1 by Ω(1) once ε is small).
  The bound δ = O(ε) satisfies ‖P²−P‖ ≤ δ and P† = P (L918–929).

- **Proof:** L934–969.  Technique: **Lefschetz–Hopf fixed-point theorem + H-space homotopy
  (purely non-constructive topological)**.
  Steps:
  1. WLOG assume A has exact unit (prop_unit).
  2. Search for fixed points of the inversion map σ: U → U† on the approximate unitary group
     𝒰 (as a C¹ manifold, prop_polar).  A fixed point U satisfies σ(U) = U, i.e. U†=U, so
     P = ¼(2I + U + U†) = ½(I+U) is an O(ε)-projection (L935–943).
  3. Pass to the quotient manifold Ǔ = 𝒰_e/U(1) (L944–945) to kill the trivial fixed points
     ±I (which are isolated, with derivative of σ equal to −I+O(ε) in local coords, L943, L951).
  4. Apply Lefschetz–Hopf: if Ǔ were the only fixed point of σ̌, the sum of fixed-point
     indices would equal 1.  But prop_H-group shows Λ(σ̌) = Σ_k dim H^k(Ǔ;ℝ) ≥ 2
     (because Ǔ is a non-zero-dimensional compact connected manifold).  Contradiction ⟹
     another fixed point exists.
  5. That non-trivial fixed point Ǔ ≠ ě lifts to some U ∈ 𝒰_e with σ(U) = e^{iφ}U;
     then P = ½(I + e^{−iφ/2}U) is nontrivial.

- **Constructive in finite dim?:  NO — the proof is purely existential.  Three concrete
  constructive routes follow (see also §"Open implementation questions").**

- **Depends on:** prop_P (L524), prop_unit (L672), prop_polar (L809), lem_U_delta (L699),
  lem_gV (L758), rem_X2 (L628), prop_H-group (L971).
- **Used by:** th_main (L1414, structural decomposition needs one nontrivial projection to
  start the induction on dim A).

- **Arbitrary precision needed?:**  Yes.  The spectral gap of X = 2P−I near ±1 (eigenvalues
  close to ±1) controls whether the lift P is truly nontrivial vs. near-zero or near-I.
  Interval arithmetic (arb/acb) required to certify the gap.

- **Cross-check / oracle:**  For a concrete matrix algebra M_n(ℂ) the answer is obvious
  (any rank-k projector for 0 < k < n).  For a perturbed algebra the constructive route
  should reproduce a projector close to a known one.

- **Implementation notes / risks:**
  - LAPACK: `zheev` on a Hermitian element X ∈ A to get eigenvalues λ_i.
  - Danger: eigenvalue very close to 0 or 1 after the spectral split — need arb/acb gap
    certification.
  - The U(1)-quotient and Lefschetz index calculation are not needed in the implementation;
    only the resulting object P matters.

---

### prop_H-group — Proposition "H-space cohomology trace" (statement L971–973)

- **Statement (faithful):**
  Let M be a connected CW complex with dim H*(M;ℝ) < ∞.  If M is an H-space with a left
  inversion map σ, then Tr σ^{*k} = (−1)^k dim H^k(M;ℝ) for every k.

- **Proof:** L1023–1050.  Technique: **Hopf structure theorem for connected commutative
  associative bialgebras + augmentation filtration (algebraic topology, non-constructive)**.
  Key steps:
  - The finite-dimensional cohomology algebra is an exterior algebra on odd-degree generators
    x₁,…,x_m (Hopf's theorem, L1014).
  - The left-inversion axiom forces σ*(x_l) ≡ −x_l mod (A⁺)², and hence
    σ*(x₁^{p₁}⋯x_m^{p_m}) ≡ (−1)^p x₁^{p₁}⋯x_m^{p_m} in the associated graded (L1043–1048).
  - Trace formula follows from this triangular form (L1049).

- **Constructive in finite dim?:**  Not needed for implementation; this result is only used
  inside the proof of lem_nontriv_projection.  **Escalate** — purely topological; no
  constructive finite-dim algorithm required or meaningful.

- **Depends on:** Hopf's theorem (Hatcher 3C.4), Künneth theorem; the H-space structure of
  Ǔ = 𝒰_e/U(1) established in §6.
- **Used by:** lem_nontriv_projection (L968, to obtain Λ(σ̌) ≥ 2).

- **Arbitrary precision needed?:**  N/A (not implemented directly).

- **Cross-check / oracle:**  N/A.

- **Implementation notes / risks:**  None; pure math used only in existence proof.

---

### lem_alpha — Lemma "block-sum compression bijection" (statement L1086–1102)

- **Statement (faithful):**
  Let P₁,…,P_p and Q₁,…,Q_q be Hermitian elements of an ε-C* algebra A with
  ‖P_j P_l − δ_{jl} P_l‖ ≤ δ,  ‖Q_k Q_m − δ_{km} Q_m‖ ≤ δ.
  Suppose P = Σ_j P_j and Q = Σ_k Q_k are δ-projections with ‖P_j P − P_j‖ ≤ δ and
  ‖Q_k Q − Q_k‖ ≤ δ, and pq(δ+ε) < c₀ (some absolute constant).
  Define α_{jk}: 𝒮_{P_j,Q_k} → 𝒮_{P,Q} by α_{jk}(X) = Co_{P,Q}(X), and
  α = Σ_{j,k} α_{jk}: ⊕_{j,k} 𝒮_{P_j,Q_k} → 𝒮_{P,Q} (⊕ with max-norm).
  Then α is a linear bijection satisfying
    ‖α‖ ≤ pq + O(pq(δ+ε)),   ‖α⁻¹‖ ≤ 1 + O(pq(δ+ε)).

- **Proof:** L1104–1119.  Technique: **explicit / constructive**.
  Constructs β_{jk}(X) = Co_{P_j,Q_k}(X) as an approximate left inverse of α and shows
  βα = 1 + γ with ‖γ‖ ≤ pq·O(δ+ε) < 1, giving invertibility via Neumann series.  Same
  for αβ.  Bound on ‖α⁻¹‖ follows from ‖β‖ ≤ 1 + O(δ+ε).

- **Constructive in finite dim?:**  Yes.  α and β are explicit maps given by compression
  operators Co_{P,Q}.  The Neumann series for (βα)^{−1} is finite in finite dimension.
  Direct route: compute α and β as matrices (acting on the finite-dim spaces 𝒮_{P_j,Q_k});
  verify ‖γ‖ < 1 numerically; invert.

- **Depends on:** Co_{P,Q} definition (L1054–1064), almost-containment inequality
  ‖Co_{P,Q}(X)−X‖ ≤ O(δ+ε)‖X‖ for X ∈ 𝒮_{P₁,Q₁} (L1068–1076).
- **Used by:** prop_delta_hominc (L1194), lem_approx (L1224), proof of th_main (L1414) —
  specifically the block-decomposition step.

- **Arbitrary precision needed?:**  Moderate.  The constant c₀ (pq(δ+ε) < c₀) must be
  certified; for p,q ≤ 2 (the only case used, L1084) this is a numerically mild condition.
  Interval arithmetic advisable for ‖γ‖ < 1 check.

- **Cross-check / oracle:**  For exact projections (δ=0, ε=0) α is an isometric bijection
  (classical corner decomposition of a C*-algebra).

- **Implementation notes / risks:**
  - Build the matrix of Co_{P,Q} as an explicit linear map on A (n²×n² matrix if A ⊆ M_n).
  - LAPACK `dgesvd` / `zgesvd` to get singular values; check ‖γ‖ < 1.
  - Risk: if δ+ε is not small relative to (pq)^{−1}, inversion blows up.

---

### lem_PQ_Hilb — Lemma "one-dim-Q corner is a Hilbert space" (statement L1123–1132)

- **Statement (faithful):**
  Let P, Q be δ-projections in an ε-C* algebra, with Q one-dimensional (dim 𝒮_Q = 1).
  Then 𝒮_{P,Q} is a Hilbert space with inner product ⟨Y,X⟩ defined by
    Y†·X = ⟨Y,X⟩ Q̃   (X,Y ∈ 𝒮_{P,Q}),
  where Y†·X = Co_Q(Y†X) and Q̃ = Co_Q(Q).
  Furthermore |⟨X,X⟩ − ‖X‖²| ≤ O(δ+ε)‖X‖²  (L1129–1131).

- **Proof:** L1135–1144.  Technique: **explicit / constructive**.
  - Linearity and conjugate-symmetry are immediate.
  - Positivity: suppose ⟨X,X⟩ < 0 for some X.  The function f(t) = ⟨X+tQ',X+tQ'⟩ is
    continuous, negative at t=0, positive for large t; by IVT it vanishes at some t>0,
    forcing X = −tQ' and ⟨X,X⟩ > 0, contradiction.  (The IVT step is existential but only
    needed for the proof that the inner product is positive — a property that can be
    independently verified numerically.)
  - Completeness follows from (L1129–1131) since the inner-product norm and the ambient norm
    are equivalent up to 1±O(δ+ε).

- **Constructive in finite dim?:**  Yes.  Given a finite-dim A, Q one-dimensional:
  1. Compute Q̃ = Co_Q(Q) (apply prop_P to Q).
  2. For any X ∈ 𝒮_{P,Q} compute Co_Q(X†X); since dim 𝒮_Q = 1 this is a scalar multiple
     of Q̃; extract the scalar as ⟨X,X⟩.
  The Hilbert space structure is concretely computable.

- **Depends on:** Co_{P,Q} definition (L1054–1064), almost-containment (L1074–1076),
  prop_P (L524).
- **Used by:** lem_PQR (L1162), lem_1d_proj (L1179), Ha^Q_{P,R} construction (L1146–1160),
  and downstream in proof of th_main.

- **Arbitrary precision needed?:**  Mild.  The bound O(δ+ε) in ‖⟨X,X⟩−‖X‖²‖ ≤ O(δ+ε)‖X‖²
  must be tracked; arb sufficient.

- **Cross-check / oracle:**  For exact unit projections (δ=ε=0) this recovers the standard
  identification 𝒮_{P,Q} ≅ Pℋ (as a Hilbert space) via the GNS-type construction.

- **Implementation notes / risks:**
  - Since dim A < ∞ and dim 𝒮_Q = 1, the Co_Q map is a rank-1 projector; the scalar
    extraction is numerically stable.
  - Risk: if Q is only approximately one-dimensional (dim 𝒮_Q not exactly 1), must first
    verify dim 𝒮_Q = 1 via rank of Co_Q.

---

### lem_PQR — Lemma "compressed product norm with 1-dim middle" (statement L1162–1168)

- **Statement (faithful):**
  Let P, Q, R be δ-projections in an ε-C* algebra, Q one-dimensional.  Then for all
  X ∈ 𝒮_{P,Q} and Y ∈ 𝒮_{Q,R}:
    |‖X·Y‖ − ‖X‖‖Y‖| ≤ O(δ+ε)‖X‖‖Y‖,
  where X·Y = Co_{P,R}(XY) is the compressed product (L1080).

- **Proof:** L1170–1177.  Technique: **explicit computation**.
  Chain of approximations:
    ‖X·Y‖ ≈ ‖XY‖ ≈ √‖(XY)†(XY)‖.
  Then (XY)†(XY) = Y†(X†X)Y + O(ε)‖X‖²‖Y‖²,
  and X†X ≈ ⟨X,X⟩Q̃ ≈ ‖X‖²Q̃ (by lem_PQ_Hilb, L1128–1131),
  giving (XY)†(XY) ≈ ‖X‖² Y†Y and hence ‖X·Y‖ ≈ ‖X‖‖Y‖.

- **Constructive in finite dim?:**  Yes.  All steps are explicit norm computations.

- **Depends on:** lem_PQ_Hilb (L1123), compressed product definition (L1080),
  ε-C* algebra axiom (L1172: ‖Y†((X†X)Y) − Y†·((X†X)·Y)‖ ≤ ε‖X†X‖‖Y‖²).
- **Used by:** lem_1d_proj (L1183), transitivity of equivalence of 1-dim projections (L1187).

- **Arbitrary precision needed?:**  No — constant-factor bounds; standard floating point
  adequate with O(δ+ε) tracking.

- **Cross-check / oracle:**  For exact projections, ‖XY‖ = ‖X‖‖Y‖ is the standard
  isometric property of the "multiplication by a minimal projection" map in a C*-algebra.

- **Implementation notes / risks:**
  - Pure norm computation; no eigensolver needed.
  - Risk: Q is only approximately one-dimensional.

---

### lem_1d_proj — Lemma "dim S_{P,Q} ≤ 1 for two 1-dim projections" (statement L1179–1181)

- **Statement (faithful):**
  If P and Q are both one-dimensional δ-projections in an ε-C* algebra, then dim 𝒮_{P,Q} ≤ 1.

- **Proof:** L1183–1185.  Technique: **proof by contradiction + lem_PQR**.
  If dim 𝒮_{P,Q} ≥ 2, pick two orthonormal (w.r.t. ⟨·,·⟩ from lem_PQ_Hilb) elements
  X, Y ∈ 𝒮_{P,Q} with Y†·X = 0 and ‖X‖, ‖Y‖ ≈ 1.  Apply lem_PQR to the triple (Q,P,Q):
  since P is one-dimensional, ‖X†·Y‖ ≈ ‖X‖‖Y‖ ≈ 1, contradicting Y†·X = 0.

- **Constructive in finite dim?:**  Yes (as a verification/oracle).  Compute dim 𝒮_{P,Q}
  via rank of Co_{P,Q}; check ≤ 1.  No existential steps needed.

- **Depends on:** lem_PQ_Hilb (L1123), lem_PQR (L1162).
- **Used by:** Equivalence relation on 1-dim projections (L1187); downstream in th_main
  for the "matrix units" / Artin–Wedderburn reconstruction.

- **Arbitrary precision needed?:**  Mild.  Rank determination of Co_{P,Q} may require
  gap certification; arb advisable.

- **Cross-check / oracle:**  In M_n(ℂ) two rank-1 projections span a space of dimension ≤ 1
  iff they are scalar multiples — straightforward to verify.

- **Implementation notes / risks:**
  - LAPACK `zheev` or SVD on the matrix of Co_{P,Q}; count singular values above threshold.
  - Threshold for "near-zero" singular value requires ε, δ-dependent tolerance.

---

## Key definitions & notation introduced here

| Symbol / term | Introduced at | Meaning |
|---|---|---|
| δ-projection | L917–919 | Hermitian P with ‖P²−P‖ ≤ δ |
| nonvanishing δ-projection | L928 | δ-projection with ‖P‖ = 1±O(δ+ε) |
| nontrivial δ-projection | L929 | δ-projection with both P and I−P nonvanishing |
| 𝒰, 𝒰_δ, Ū_δ | L692–696 | (approximate) unitary group of A |
| inversion map σ | L939–941 | σ(U) = U† on 𝒰; fixed pts satisfy U†=U |
| Ǔ = 𝒰_e/U(1) | L944 | quotient manifold used in topological argument |
| Co_{P,Q} | L1055–1057 | compression map θ(L_P R_Q + R_Q L_P − 1): exact idempotent |
| 𝒮_{P,Q} | L1064 | image of Co_{P,Q}; the "corner" subspace PAQ |
| compressed product X·Y | L1078–1081 | Co_{P,R}(XY) for X ∈ 𝒮_{P,Q}, Y ∈ 𝒮_{Q,R} |
| one-dimensional δ-projection | L1066 | P with dim 𝒮_P = 1 |
| ⟨·,·⟩ (Hilbert ip on 𝒮_{P,Q}) | L1125–1127 | Y†·X = ⟨Y,X⟩Q̃ when Q is 1-dim |
| Ha^Q_{P,R} | L1146–1150 | symmetrized "right-action" map 𝒮_{P,R}→Bo(𝒮_{R,Q},𝒮_{P,Q}) |
| equivalence of 1-dim projections | L1187 | P ~ Q iff dim 𝒮_{P,Q} = 1 |

---

## Dependencies OUT / INTO

**INTO this shard (used here, proved elsewhere):**
- prop_P (L524): θ(2P−I) produces an exact idempotent from a near-idempotent; used to
  define Co_{P,Q} (L1057) and the projection-snapping step.
- prop_unit (L672): WLOG exact unit for lem_nontriv_projection.
- prop_polar (L809): 𝒰 is a C¹ manifold; used in lem_nontriv_projection proof.
- lem_U_delta (L699): norm bounds on approximate unitaries.
- lem_gV (L758): local parametrization of 𝒰 near V.
- rem_X2 (L628): motivation for why exact projections are hard; justifies working with
  δ-projections throughout.

**OUT of this shard (proved here, used elsewhere):**
- lem_nontriv_projection → th_main (L1414): the induction on dim A starts here.
- Co_{P,Q}, 𝒮_{P,Q}, compressed product → lem_alpha (L1086), lem_approx (L1224),
  lem_merging (L1325), lem_add_dim (L1363), prop_delta_hominc (L1194).
- lem_alpha → prop_delta_hominc (L1194), lem_approx (L1224).
- lem_PQ_Hilb, lem_PQR, lem_1d_proj → equivalence of 1-dim projections (L1187), Ha^Q
  construction (L1146), proof of th_main.

---

## Open implementation questions / escalations

### Constructivization of lem_nontriv_projection (THE KEY PROBLEM)

The paper's proof certifies existence of a nontrivial O(ε)-projection via a topological
fixed-point argument on Ǔ.  For finite-dimensional matrix algebras, three concrete
constructive routes are available:

---

**Route A — Hermitian eigensolver spectral split (recommended)**

Given ε-C* algebra A (concretely: a finite-dim operator space A ⊆ M_N(ℂ) with approximate
structure).

1. Find any Hermitian element H ∈ A with H ≠ O(ε)·I and H ≠ I − O(ε)·I.
   (Such an H exists because dim A > 1: pick any non-central self-adjoint element, or
   form H = (1/2)(A + A†) for a non-scalar A ∈ A.)
2. Compute eigenvalues of H (as a matrix in M_N): λ₁ ≤ … ≤ λ_N.
3. Find a spectral gap: locate t such that no eigenvalue lies in (t−g, t+g) for some gap
   g > 0.  Take P_spec = 1_{(−∞,t]}(H) (projection onto eigenvalues ≤ t).
4. Snap to A: set X = 2P_spec−I and apply prop_P to map U = sgn(X) ∈ 𝒰; then P = ½(I+U).
   By the argument of L935–943, P is an O(ε + 1/g)-projection in A.
5. Check nontriviality: verify ‖P‖ ≈ 1 and ‖I−P‖ ≈ 1.

Trade-offs:
- Pro: Entirely explicit; O(N³) via LAPACK `zheev`.
- Pro: Works for any A with a Hermitian element having a spectral gap.
- Con: Requires a spectral gap g = Ω(1); if all eigenvalues of all Hermitian elements in A
  cluster near 0 or near 1, this fails — but that would contradict dim A > 1 for an
  ε-C* algebra (a deeper structural argument shows such a gap must exist).
- Con: Finding H ∈ A with a useful gap may require scanning; for dim A small (≤ ~100) this
  is cheap.

---

**Route B — Apply θ to a near-involution (paper's own implicit route)**

The paper (L935–943) already contains the key idea: any U ∈ 𝒰 with ‖U†−U‖ ≤ δ produces a
δ-projection P = ¼(2I + U + U†).  So:

1. Search 𝒰 for a Hermitian element: minimize ‖U†−U‖ over U ∈ 𝒰 subject to U ≠ ±I+O(ε).
2. Equivalently: in the coordinate chart φ_I on 𝒰 (tangent space = iℋ, L788–799), find a
   near-real point; project to ℋ ∩ 𝒰.
3. Once such U is found, P = ½(I+U) is an O(ε)-projection.

Concrete algorithm: parametrize 𝒰 locally via the exponential map exp(iH) for H ∈ ℋ.
The inversion map σ acts by H → −H in this chart.  Fixed points of σ are H = 0 (mod 2π),
but at H = π·v for a unit vector v they give non-trivial P = ½(I + exp(iπv)).  The
problem reduces to finding H ∈ ℋ with ‖H‖ ≈ π and exp(iH) ∈ 𝒰.

Trade-offs:
- Pro: Very direct; equivalent to Route A in the spectral domain.
- Con: Optimization over 𝒰 in a non-convex space is needed; may need gradient descent.
- Con: Correctness certification (is U truly far from ±I?) requires arb interval arithmetic.

---

**Route C — Minimize ‖P²−P‖ over Hermitian P ∈ A with ‖P‖ and ‖I−P‖ bounded away from 0**

Pose the constrained optimization:
  minimize ‖P²−P‖   subject to P† = P, P ∈ A, ‖P‖ ≥ c₁, ‖I−P‖ ≥ c₁
for some c₁ = Θ(1).

1. Project the feasible set to a parametrized family: P = Σ_j a_j B_j where {B_j} is an
   orthonormal basis of the Hermitian subspace of A.
2. Use a semidefinite or eigenvalue-constrained optimizer (or coordinate descent) to find a
   local minimizer.
3. Check ‖P²−P‖ ≤ O(ε).

Trade-offs:
- Pro: Does not require knowing a spectral gap in advance; finds the "best" approximate
  projection automatically.
- Con: Non-convex optimization; convergence not guaranteed in general.
- Con: The existence of a feasible point is exactly what the paper proves non-constructively;
  so this route needs a warm-start from Route A or B, or must be run with random restarts.
- Recommended use: as a refinement step after Route A or B has found a candidate P.

---

**Recommended strategy:**  Use Route A as the primary algorithm (Hermitian eigensolver, O(N³),
certifiable via interval arithmetic on the spectral gap), optionally followed by Route C as a
polishing step.  Route B is equivalent to Route A when A = M_n(ℂ) but more expensive in the
general approximate-algebra setting.

---

**Spectral gap certification (arb/acb):**

For all routes, the critical numerical risk is that the "gap" g between eigenvalue clusters
near ±1 is too small to certify that P is truly nontrivial.  In arb/acb:
- Compute eigenvalues of H as balls [λ_i ± r_i]; verify that for some t, all balls lie
  strictly below t or strictly above t.
- The paper's bound ε is the a priori gap lower bound (from dim A > 1); use this as a
  floor for g.
- If no certified gap is found, escalate: the input algebra may be ε-too-large or not
  genuinely an ε-C* algebra.

---

**Question for th_main implementors:**  The proof of th_main uses lem_nontriv_projection
as an entry point (split A by P, recurse on the two corners 𝒮_P and 𝒮_{I−P}).  The
constructive Route A therefore requires the implemented inverse of prop_P (i.e., the map
from a spectral projector to an element of A) to be tracked carefully through the
induction.  Flag as a dependency when implementing L1414.
