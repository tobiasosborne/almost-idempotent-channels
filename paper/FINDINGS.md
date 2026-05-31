# FINDINGS.md — a living log of paper issues, constructivizations, and load-bearing subtleties

> **Purpose.** As we implement Kitaev, *Almost-idempotent quantum channels and
> approximate $C^*$-algebras* ([arXiv:2405.02434](https://arxiv.org/abs/2405.02434))
> as constructive finite-dim algorithms, we keep hitting (a) **typos / formula
> errors** in the `.tex`, (b) **non-constructive steps** that need a constructive
> replacement, (c) **load-bearing subtleties** that are easy to get wrong (and that
> our hostile reviews keep catching as "tests that can't fail"), and (d) **open
> questions** the paper leaves implicit. This file collects them all in one place
> so they are not rediscovered the hard way. **Append to it as you go** — every
> module's research/implementation/review should add its findings here, and cite
> this file from the source comment where the finding bites.
>
> Canonical ground truth is still `paper/src/approximate_algebras.tex` (cite line
> numbers). This file records *deviations from* and *clarifications of* it.
>
> **Status legend:** `CONFIRMED` (verified, e.g. by numerical measurement or a
> proof) · `FLAGGED` (suspected, not yet hard-verified) · `RESOLVED` (handled in
> code, route documented) · `OPEN` (escalation, unresolved).

---

## A. Typos / formula errors in the `.tex`

Per CLAUDE.md stop conditions, we **flag** these and document the correct reading;
we do **not** silently "fix" the math. Where the correct reading is unambiguous we
implement it (citing the justification); where it is a genuine error we escalate.

### A1. `tex:1109` — `lem_alpha`'s β subscript: `Co_{P_j,Q_j}` should be `Co_{P_j,Q_k}`
- **Status:** CONFIRMED (clear typo; correct reading unambiguous). Surfaced: `corner` Increment 2a (bead aic-czm).
- **What the `.tex` says:** `β_{jk}: X ↦ Co_{P_j,Q_j}(X) : S_{P,Q} → S_{P_j,Q_k}`.
- **Why it's wrong:** the stated codomain is `S_{P_j,Q_k}`, and the proof at `tex:1114`
  needs `β_{jk}α_{lm}(X_{lm}) = δ_{jl}δ_{km}X_{lm}+…` — the Kronecker `δ_{km}` forces
  the second index to be `Q_k`.
- **Correct:** `β_{jk} = Co_{P_j,Q_k}`. Implemented as `Q_k` (`src/aic_corner_alpha.c`),
  with the discrepancy in the source comment; the `Q_j` build is mutation-RED (‖γ‖
  jumps 0→1).

### A2. `tex:1254` — the direct-sum diagonal formula is non-central for the finite Pauli 1-design
- **Status:** CONFIRMED (numerically measured). Surfaced: `dhom` (bead aic-c1n), hostile review.
- **What the `.tex` says:** the diagonal of `⊕_l M_{d_l}` is the **Cartesian product**
  over per-block Pauli tuples `j=(j_1,…,j_m)` with the **joint** block-diagonal
  unitary `U_{j_1…j_m}=U_{1j_1}⊕…⊕U_{mj_m}` and weight `∏_l p_{lj_l}`.
- **Why it's wrong (for the *finite* Pauli design):** a diagonal must satisfy
  `XD=DX` (`tex:1241`). The Haar diagonal `D=∫dU U†⊗U` over `∏_l U(d_l)` has its
  cross-sector (`a≠b`) terms vanish because the unitary **first moment**
  `∫_{U(d)} U dU = 0`. But the finite generalized-Pauli sum has **nonzero first
  moment** `Σ_{jk} S_{jk}` (a rank-1 matrix, *not* 0), so the joint Cartesian
  product leaves **spurious cross-sector terms** → **non-central**. Measured on the
  unmodified literal transcription: `‖XD−DX‖_op ≈ 0.54` for `M₂⊕M₂`, `0.61` for
  `M₂⊕M₃`. This silently degrades `lem_approx` from quadratic to linear convergence.
  (The per-block Pauli sum *does* correctly reproduce the single-block Haar second
  moment — the single-block case is fine; only `m≥2` is broken.)
- **Correct finite central diagonal:** the **cross-term-free embedded sum**
  `D = Σ_l D_l = Σ_l Σ_{jk} d_l⁻² (Ŝ^{(l)}_{jk})† ⊗ Ŝ^{(l)}_{jk}`, where
  `Ŝ^{(l)}_{jk}` is `S^{(l)}_{jk}` on block `l` and **zero elsewhere** (a partial
  isometry — *not* the joint unitary). This is exactly the Haar diagonal (cross
  terms absent by construction). `nterms = Σ_l d_l²` (a **sum**, not the product
  `∏_l d_l²`). Implemented in `src/aic_dhom_diag.c` (the `.tex:1254 DISCREPANCY`
  box documents it in full).
- **Consequence to watch:** in this cross-term-free representation `‖D‖_proj = m`
  (the **block count**), not 1 (the norm-1 property `tex:1247` belongs to the unique
  *optimal* representation; the cross-sector terms that would compress it to 1 are
  exactly the non-central ones we removed). So the `w'` bound becomes `O(mδ)` — a
  number-of-blocks dependence. Empirically (`test_dhom` T5 block-count canary
  m=1,2,3) the th_main constant does **not** grow with `m`, but this is tracked
  against the universal-constant claim (`tex:484/460`; see D2).

---

## B. Non-constructive steps and our constructive routes (Law 3)

The paper leans on non-constructive tools; in finite dimensions we replace them.
Each entry: the paper's technique → our constructive route.

### B1. `lem_nontriv_projection` (`tex:931`) — Lefschetz–Hopf → ambient spectral split
- **Status:** RESOLVED (route in `src/aic_projection*.c`, bead aic-mqf). The
  *existence* guarantee (Ω(1) gap) remains OPEN — see D1.
- **Paper:** existence of a nontrivial projection via a Lefschetz–Hopf fixed-point
  argument on the approximate-unitary quotient manifold (`tex:944-969`) — purely
  existential.
- **Ours:** reduction (`tex:935`) `P=½(I+X)` for a Hermitian near-involution `X∈A`.
  Get `X` as the **ambient** sign of a gap-shifted Hermitian element of A:
  pick non-scalar Hermitian `H∈A` (largest interior eigenvalue gap), form
  `X=sgn(s(H−tI))` in `M_n` (eig-free funcalc), `P_amb=½(I+X)`, then project into A
  via `Φ̃` (see C3). Avoids the §5 unitary group entirely.
- **Why ambient, not in-A:** rem_X2 (see C1).

### B2. The diagonal `D=∫dU U†⊗U` (`tex:1245`, Haar) → explicit Pauli closed form
- **Status:** RESOLVED (`src/aic_dhom_diag.c`). But the direct-sum form needed the
  A2 correction.
- **Paper:** `D` defined as a Haar integral (non-constructive), then asserted to be
  a finite convex combination by Carathéodory.
- **Ours:** the explicit generalized-Pauli (Heisenberg–Weyl) form `tex:1250`,
  `D=Σ_jk d⁻²S_jk†⊗S_jk`, `S_jk=X^jZ^k`. Exact, eig-free. (Direct sum: see A2.)

---

## C. Load-bearing subtleties / hallucination risks (confirmed in practice)

These are the things that *look* right but aren't, and that the hostile reviews
keep catching as "tests that cannot fail." Several restate CLAUDE.md callouts but
with the concrete evidence from where they bit.

### C1. rem_X2 (`tex:628`) — no functional calculus *inside* the ε-C* algebra
- **Status:** CONFIRMED load-bearing (projection crux, bead aic-mqf; 3 independent
  research legs + web survey).
- The matrix-sign / `prop_P` functional calculus does **not** generalize to the
  ε-C* algebra A (the star is only ε-associative; `Y⋆Y=I` may have no solution near
  a target). Do the calculus **ambient** (in `M_n`, an exact C* algebra), then
  project the result into A. No published work supports sign-iteration convergence
  under ε-associativity. (My first instinct — star-Newton-Schulz inside A — was
  wrong; legs 2+3 caught it.)

### C2. The product on A is the **Choi–Effros star** `X⋆Y=Φ̃(XY)`, not plain `XY`
- **Status:** CONFIRMED; the #1 test-blindness in *every* module (`tex:341`, `:2189`).
- Plain `acb_mat_mul` for an "A-product" leaves `Img Φ̃` and is a bug. **The η=0
  identity-channel oracle is structurally BLIND to this** (there `Φ̃=id` so star≡plain).
  So every module that uses the star MUST test on an **oblique / genuine-η>0**
  fixture (`make_compress_idemp`, `make_mixconj`, `make_eta_family`) with a
  **non-vacuity guard** (mutate star→plain → confirm RED). Caught in corner
  (Inc 1, both star and left/right-singular-vector blind), corner cdot (Inc 2a),
  dhom (T4). Use `acb_mat_mul` only for genuine ambient/coordinate products.

### C3. `A=Img Φ̃` is an **oblique** image — project into A with `Φ̃`, not the Frobenius `Π_A`
- **Status:** CONFIRMED (projection module, bead aic-mqf — corrected the research spec).
- `Φ̃` is Hermicity-preserving but **not** HS-self-adjoint, so the
  Frobenius-orthogonal projector `Π_A(W)=Σ_k⟨B_k,W⟩_F B_k` does **not** respect A's
  star structure — it leaves the star defect `‖P⋆P−P‖=O(1)` (~0.5, constant in η).
  The correct projection onto `A=Img Φ̃` is **`Φ̃` itself**, available through the
  public star API as `W⋆I = Φ̃(W)`, which gives the `O(η)` defect.

### C4. Oblique idempotents: nonzero singular values are `≥1`; extract the range via LEFT singular vectors
- **Status:** CONFIRMED (assoc, corner; web leg Antezana–Corach 2009).
- An oblique (non-Hermitian) idempotent `E` has `σ_i≥1` on its range (`=1` iff
  orthogonal), `≈0` on its kernel — so `round(trace)` + an SVD-gap-at-0.5 count both
  give the rank, and the **top-`rank` LEFT singular vectors** span `Img E` (the RIGHT
  ones span the co-range `Img E†` — a different space; reversing them is a bug
  invisible to orthogonal-projector fixtures). The η=0 identity/orthogonal fixtures
  are blind to the left/right distinction — needs an oblique fixture.

### C5. The `aic_mat_opnorm` Gram-path Hermiticity false-fail (implementation, recurring)
- **Status:** RESOLVED (bead **aic-qgs**; root-cause fix in `src/aic_mat_norms.c`
  `aic_mat_gram`; regression `tests/test_mat.c` test8, mutation-proven). The
  `aic_corner_gamma_opnorm_ub` mid+radius workaround can now be **progressively
  retired** (NOT in the aic-qgs change — a separate cleanup).
- `aic_mat_opnorm` / `aic_mat_singular_values` form the Gram `M†M` and route through
  the relative-Hermiticity predicate (`aic_mat_int_is_hermitian`), which **false-failed
  and SIGABRTed** when a deep matmul chain that CANCELS feeds in a matrix whose Gram
  entries are SMALL in magnitude (`~1e-6`) yet carry a LARGE matmul-**accumulated**
  arb radius (`~1e-72`): the magnitude-relative floor `tol·(1+|G_ij|+|G_ji|)` is then
  `~tol` (since `|G|≪1`), `~1e3×` smaller than the radius. The predicate tests the
  radius of the DIFFERENCE ball `|G_ij−conj(G_ji)|`, so it fired even though `G=M†M`
  is genuinely Hermitian. Hit by corner (`lem_alpha` γ, Ha-map `G−I`, `lem_PQR`),
  projection (basin assert), dhom, and the oblique S_P-wrapper corner (§C10).
- **The diagnosed fix (a) — symmetrize `G←(G+G†)/2` — is INSUFFICIENT** (the bead's
  diagnosis missed this): the off-diagonal midpoints of `M†M` are ALREADY exact
  conjugates, so the residual the predicate flags is purely the difference-ball
  RADIUS, which no midpoint manipulation removes (arb subtraction adds the radii;
  verified empirically). **The root-cause fix (rigorous, the bead's option (b) in
  spirit):** `aic_mat_gram` splits the certified Gram into an EXACTLY-Hermitian
  midpoint matrix `Gmid` (zero radius → the predicate passes for the right reason)
  plus a rigorous Weyl perturbation bound `R ≥ ‖G_true−Gmid‖_op ≤ ‖G_true−Gmid‖_F`;
  the certified eig runs on `Gmid` and the eigenvalues are inflated by `R` (Weyl's
  inequality for Hermitian matrices) before the `sqrt`. This is the substrate version
  of the `aic_corner_gamma_opnorm_ub` mid+radius idea. It changes NO value on tight
  inputs (`R~2^-prec` there, below the existing global-enclosure radius — the
  test_mat exact special cases, precision ladder, singular values, and aic-2yo graded
  Gram all stay green with the SAME numbers) and does NOT touch the
  `herm_max_eig`/`is_hermitian` guard, intact for its DIRECT callers (CP-cert
  `herm_max_eig(-C)`, the aic-2yo teeth).

### C6. δ-inclusion lower bound: the basis sweep is BASIS-BLIND — use σ_min of the coordinate matrix
- **Status:** RESOLVED (route in `src/aic_dhom_sigmin.c`; guard switch in
  `src/aic_errreduce.c`). Surfaced: `errreduce` (bead aic-t81), hostile review (F1
  BLOCKER, a soundness hole).
- **The trap.** The δ-inclusion hypothesis (`tex:451-453`) is
  `(1−δ)‖X‖ ≤ ‖v(X)‖ ≤ (1+δ)‖X‖` over the **OPERATOR norm**, i.e. a bound on the
  unit-ball infimum `inf_{X≠0} ‖v(X)‖/‖X‖`. The natural cheap surrogate is the
  basis sweep `min_i ‖v(E_i)‖` (dhom's `aic_dhom_prop_bounds` `norm_lb`). That is
  **NOT** the inclusion infimum: a `v` bounded below on every basis element can
  still **collapse a general combination**. Witness (review): `B = C⊕C → A = M₂`,
  `v(E₀)=diag(1,0)`, `v(E₁)=|u⟩⟨u|` at angle 0.1. Each `‖v(Eᵢ)‖_op = 1` (so the
  basis sweep reads `a = 1.0`, **passes** the `≥0.5` guard), but
  `‖v(E₀−E₁)‖_op = |sin 0.1| = 0.0998` while `‖E₀−E₁‖ = 1`, so the true inclusion
  inf is `≤0.0998` — the hypothesis is violated and `aic_errreduce` silently
  certified it as a 0-inclusion (`c₀=0`, no abort). A "test that can't fail."
- **The sound route — σ_min of the coordinate matrix.** `v(X) = Σᵢ xᵢ v(Eᵢ)` is
  linear. With A's Frobenius-orthonormal basis `{B_k}` and B's Frobenius-orthonormal
  matrix units `{Eᵢ}`, assemble `M` (`dim_A × dim_B`), `M[k,i] = ⟨B_k, v(Eᵢ)⟩_F`.
  Then `‖v(X)‖_F = ‖Mx‖₂` and `‖X‖_F = ‖x‖₂`, so
  `σ_min(M) = inf_{X≠0} ‖v(X)‖_F/‖X‖_F` — the **exact unit-ball inclusion infimum
  in the Frobenius/coordinate norm**. It SEES all combinations (`σ_min = 0` iff `v`
  collapses a direction), so it is a **sound collapse detector**.
- **Frobenius vs operator (the documented caveat).** `σ_min(M)` is the
  Frobenius/coordinate inf, not the exact operator-norm inclusion inf the `.tex`
  states; they differ by `≤√dim` factors (norm equivalence). It is therefore a true
  Frobenius unit-ball lower bound and a correct REJECTER of non-inclusions, but not
  the precise operator-norm `a`. The faithful operator-norm worst-case (HOPM, like
  `aic_ecstar`) is a later cycle (bead aic-0at). The σ_min GATE uses the **double
  midpoint** (`aic_latd_singular_values` on `mid(M)`, uncertified — certified
  enclosure defers to aic-w4o.1/.2), a coarse fail-loud gate adequate for the 0.5
  threshold (cf. the projection-nontriviality gate).
- **Where it bites.** Switched BOTH `aic_errreduce` guards (input δ-inclusion check
  AND the certification `lower_gap = max(0, 1−a)`) and the `aic_errreduce_is_bijective`
  injectivity test from the basis sweep to `aic_dhom_v_sigma_min`. `test_errreduce`
  T6 is the witness fixture + abort test; mutation-proven (revert to the basis sweep
  → the collapse is no longer caught → RED).

### C7. `aic_ecstar_defect_unit` tests the AMBIENT `1_n` unit — wrong for a compressed S_P (unit is `Ptilde`)
- **Status:** CONFIRMED (cstar_build I1, bead aic-097, S_P wrapper genuine-C* oracle).
- The unit-law estimator `aic_ecstar_defect_unit` (`include/aic_ecstar.h:183-189`) is
  HARDCODED to the ambient unit `I = 1_n`: it checks `‖1_n − Π_A(1_n)‖`, `‖B_k⋆I − B_k‖`,
  `‖I⋆B_k − B_k‖`. That is correct only when the algebra's unit IS `1_n` (`A = Img Φ̃`,
  whose unit is the inherited `1_n`, `tex:2186-2187`; and the genuine `M_d`,
  `aic_cstar_matrix_algebra`). It is **WRONG for a compressed subalgebra `S_P`**: the
  unit of `(S_P, ⋅)` is `Ptilde = Co_P(P)` (`tex:1082`), NOT `1_n` — and `1_n ∉ S_P`
  in general. On the η=0 identity-channel `A=M_3`, `P=diag(1,1,0)` → `S_P ≅ M_2`
  wrapper, the four unit-INDEPENDENT defects (assoc, submult, C*, involution) all read
  machine-zero (it IS a genuine C* algebra), but `aic_ecstar_defect_unit` reads **1.0**
  because `Ptilde ≠ 1_3`. This is the same load-bearing distinction CLAUDE.md flags
  ("the exact unit is only available after the O(ε)-change of `prop_unit`").
- **The route.** Test the genuine-C* unit law for `S_P` against `Ptilde` directly
  (`‖Ptilde⋆C_m − C_m‖`, `‖C_m⋆Ptilde − C_m‖` over the corner basis), NOT via
  `aic_ecstar_defect_unit`. `tests/test_cstar.c` (T1) uses `max_axiom_defect_no_unit`
  + a `ptilde_unit_defect` helper; both read machine-zero at η=0 (measured 0.0).
  A unit-parametrized `aic_ecstar_defect_unit(..., unit)` would let any subalgebra
  reuse the estimator; deferred (not needed for I1, the helper suffices). Note the same
  estimator also hits the `aic_mat_opnorm` Gram false-fail (§C5) on the OBLIQUE wrapper,
  so the η>0 S_P star-defect (T2) is measured via the midpoint-opnorm route too.

### C8. Merged-`v` star teeth: on IN-A delta-projection fixtures the η-shrink DIRECTION is a weak tooth — use the `c = defect/η` MAGNITUDE bound (a §C2 corollary)
- **Status:** CONFIRMED (cstar_build I2, bead aic-097, `cor_merge_sum` B2 star teeth).
- When the merged-`v` multiplicativity defect (`aic_dhom_defect_sweep`, which uses A's
  STAR) is exercised on a fixture whose `P_j` are GENUINE in-A delta-projections
  (residual `≈0`), mutating the star → plain product does **NOT** make the defect blow
  up to O(1): `Ptilde_j` stays near-genuine off the oblique direction, so the plain-
  product merged defect ALSO trends down with η. Measured (mixconj(5,3),
  `P_1=span(e1),P_2=span(e2)`): star `c=defect/η ≈ 0.017`, plain `c ≈ 0.43` — both
  shrink in absolute terms, so a "defect shrinks as η→0" assertion passes for BOTH (a
  weak/half-blind tooth). The SHARP discriminant is the `c`-ratio **magnitude** (25×
  gap: a `c < 0.2` bound is RED for the plain mutation). 
- **How to apply.** Any test asserting a merged/assembled-`v` defect is O(η) with the
  STAR must gate on the `c = defect/η` magnitude (a small constant), NOT merely on the
  η→0 shrink direction. This recurs in I3 (`lem_merging`) and I5 (the loop) — both
  certify merged-`v` defects. (The complementary route: test on a genuinely-oblique
  `P` not in A, where plain-vs-star is O(1) even in direction — but then the unit
  defect itself is O(1), so the two concerns split across two `P`'s, as in I1 T2.)
- **I3 confirmation.** `lem_merging` C2 (`two_block` shape, mixconj(5,3)) measured
  star `c ≈ 0.017 → 0.005`; the plain-product mutation (overwrite `ae.A.star_phi`)
  gives `c ≈ 0.244 → 0.230` (`> 0.2`) — RED, mutation-proven by hand (`/tmp/teeth_c2`).
  (Magnitudes smaller than I2's 0.43 because `lem_merging`'s two-block defect sweep
  has only the single cross-pair, like `cor_merge_sum`; the `c < 0.2` bound still
  cleanly separates star from plain.)

### C9. `lem_merging`'s `B` has TWO shapes — single matrix block (lem_extension) vs two-block direct sum (cor_merge_sum); a single block with zeroed off-diagonal is an INVALID input
- **Status:** CONFIRMED (cstar_build I3, bead aic-097, `lem_merging`). The increment
  prompt's "single `M_{n1+n2}` block + `gamma_12=gamma_21=0` reduces to cor_merge_sum"
  reading is **mathematically incorrect**; corrected here.
- `lem_merging` (`tex:1325`) only requires `B` to be "a C* algebra" with two
  COMPLEMENTARY projections `Π_1+Π_2=I`. Two shapes arise in the master loop:
  - **single block** (`two_block=0`): `B = M_{n1+n2}`, `Π_1` = proj onto the first
    `n1` coords, `Π_2` onto the last `n2`. The off-diagonal sectors `S_{Π_1,Π_2}`,
    `S_{Π_2,Π_1}` are the **LIVE** rectangular blocks of `M_{n1+n2}` (the units `E_lm`
    with `l,m` across the `n1` boundary EXIST), so `gamma_12, gamma_21` must be
    NONZERO. This is what **lem_extension** (I4, `tex:1378`, `v_+ : M_{n+1} → A`)
    needs (`gamma_12 = U_1`, `gamma_21 = U_1(·)^†`).
  - **two block** (`two_block=1`): `B = M_{n1} ⊕ M_{n2}`, `Π_j = I` of block `j`. The
    off-diagonal sectors are **EMPTY** (block-diagonal `B` has no cross-block units),
    so `gamma_12 = gamma_21 = 0` (zero-dim domain). This is what **cor_merge_sum**
    (`tex:1352`) uses; the merged `v` then equals `aic_cstar_merge_sum`'s concat
    EXACTLY (I3 C2 verified to machine precision: `mult_def`/`sigma_min`/`vE` all match
    `cor_merge_sum` at η ∈ {3.7e-2, 1.3e-2}).
- **The trap (a "test that can't fail").** A single `M_{n1+n2}` block with
  `gamma_12=gamma_21=0` is NOT a valid `lem_merging` input: the LIVE pair
  `(E_{0,n1}, E_{n1,0})` has `E·E = E_{00}` in `M_d` but `gamma(E)=0`, so merging1 is
  violated and `mult_def = ‖Ptilde_1‖ ≈ 1.0` — **O(1)**, NOT the `cor_merge_sum`
  cross-defect `‖Ptilde_1 ⋆ Ptilde_2‖ ≈ 1.9e-5`. Measured directly (`/tmp/probe_c2`):
  `‖Ptilde_1⋆Ptilde_2‖ = 1.90e-5` vs `‖Ptilde_1‖ = 1.00`. Conflating the two shapes
  silently mis-assembles `B` and is exactly the off-by-`n1` a hostile review hunts.
- **The route.** `aic_cstar_lem_merging` takes an explicit `int two_block` selector;
  `B_out` is built as `M_{n1+n2}` (single) or `M_{n1}⊕M_{n2}` (two) and the routing +
  `merge_cond` sweep branch on it. I4 will pass `two_block=0`; Stage-3 of the loop
  (which is `cor_merge_sum`) is already served by `aic_cstar_merge_sum` directly, so
  `two_block=1` is mainly the cross-check seam.

### C10. The OBLIQUE `S_P`-wrapper corner path hits the aic-qgs Gram false-fail in the `sgn`-basin opnorm (a new §C5 manifestation); blocks the genuinely-oblique `lem_extension` end-to-end on a compressed parent
- **Status:** RESOLVED by the §C5 / aic-qgs `aic_mat_gram` fix (exact-Hermitian
  midpoint Gram + Weyl `R` inflation). The blocker below is GONE: `aic_corner_dim_S`
  / `aic_corner_Co` now COMPLETE on the oblique `S_P` wrapper of `make_mixconj(5,3,
  0.06)` over the 2-dim corner `span(e1,e2)` — verified by the aic-qgs RED→GREEN
  probe (pre-fix: SIGABRT in the opnorm Gram-Hermiticity path; post-fix: `dim
  S_{P,Q}=1, dim S_P=1, dim S_Q=1`). **Follow-up:** the genuinely-oblique
  `S_P`-WRAPPER-as-parent end-to-end `lem_extension` (deferred from I4 to aic-qgs)
  can now be implemented as a `test_cstar_extension` T3 leg, and the
  `aic_corner_gamma_opnorm_ub` workaround retired — both separate cleanups, NOT in
  the aic-qgs substrate change. — Original report (cstar_build I4, bead aic-097,
  `lem_extension`): a manifestation of the §C5 / aic-qgs
  `aic_mat_opnorm` Gram-Hermiticity false-fail,
  surfaced because I4 is the FIRST customer to run the corner machinery (`aic_corner
  _Co` / `dim_S` / `ha`) on a **genuinely-oblique (η>0) parent** — either an `S_P`
  wrapper (`aic_cstar_subalg`) used as `A_parent`, or a raw `ae.A` with an OBLIQUE
  projection (`Ptilde_P = v(E_jj)`) as the first `dim_S` argument. (`test_cstar.c` T3
  only exercised the η=0 identity-channel wrapper, whose Gram is exactly identity.)
- **Where it bites.** `aic_corner_Co` → `aic_prop_P` → `aic_theta` → `aic_sgn` →
  `aic_funcalc_int_in_op_basin` → `aic_funcalc_int_def_X2` → `aic_mat_opnorm` on
  `(2M-1)^2 - I` (`M = ½(L_P R_Q + R_Q L_P)`). For the oblique corner this matrix is a
  near-identity with tiny matmul-ACCUMULATED arb radii on its off-diagonals; the
  relative-Hermiticity check in `aic_mat_opnorm`'s Gram path FALSE-FAILS and SIGABRTs.
  **Measured (mixconj(5,3), wrapper over span(e1,e2)):** `‖(2M-1)^2 - I‖_ub ≈ 1.08e-3
  (Co_{P,P}), 4.46e-3 (Co_{Q,Q}), 3.12e-3 (Co_{P,Q})` — all DEEPLY in-basin (< 1); the
  abort is purely the Hermiticity false-fail, NOT a genuine out-of-basin input.
- **The route taken (no infra bandaid; Rule 3 / stop conditions).** `aic_prop_P`
  itself uses the EIG-free Gelfand spectral route (no opnorm), so the only offender is
  the `sgn`-basin DISPATCH opnorm. The proper fix is aic-qgs (make `in_op_basin` use a
  certified UB / Frobenius bound like `aic_corner_gamma_opnorm_ub`); NOT touched in I4.
  In `aic_cstar_lem_extension` the Step-1 per-`j` guard `dim S_{P_j,Q}==1` (with
  `P_j = v(E_jj)`, oblique) is REPLACED by the equivalent checks `dim S_Q==1` and the
  CONCLUSION `dim S_{P,Q}==n` (both queried with GENUINE Hermitian projections `P`,`Q`,
  which do NOT trip the false-fail) — the conclusion is the operative Rule-4 guard; the
  per-`j` is the proof's derivation path, redundant as a guard. The Ha-map loop
  (`ha(v(E_lm), P, P, Q)`) and `aic_corner_Ptilde` accept the oblique `v(E_lm)` as the
  Z OPERAND without aborting (only the first-PROJECTION-arg `dim_S`/`Co` path trips),
  so `lem_extension` RUNS end-to-end on the raw oblique `ae.A` (test T3b).
- **Consequence for I4 testing.** The genuinely-oblique `S_P`-WRAPPER-as-parent
  end-to-end `lem_extension` is DEFERRED to aic-qgs. T3 instead splits the two concerns
  (the §C7/§C8 split pattern): T3a (η=0 wrapper-as-parent, proves the compressed-parent
  path with `‖P+Q-unit‖=0` exactly) + T3b (genuinely-oblique end-to-end on raw `ae.A`,
  the §C8 `c=mult_def/η` star MAGNITUDE tooth, `c_star≈0.017 ≪ 0.2 ≪ c_plain≈0.50`,
  a ~29× gap). T3b is a CORNER-LOCAL extension (`P+Q = span(e1,e2) ≠ 1_5`, printed
  honestly); the valid non-vacuous content is the merging1 multiplicativity star tooth.

### C11. I5 master loop (`aic_cstar_build`, proof of th_main) — load-bearing subtleties, the promoted §G findings, and the open coverage gaps
- **Status:** CONFIRMED (cstar_build I5, bead aic-097). This entry is the canonical
  home for the design-doc §G-series findings that became load-bearing in shipped I5
  code (the code cites `FINDINGS §C11`; the detailed derivations remain in
  `docs/research/cstar_build_design.md` §6 G1/G5 and `cstar_masterloop_spec.md`).
- **Nontriviality-on-wrapper (promoted design §G1).** `aic_projection_nontrivial`'s
  internal nontriviality assert uses the AMBIENT `‖1_n − P‖` (vacuous on an S_P
  wrapper, whose unit is `Ptilde_m ≠ 1_n`). The Stage-1 greedy loop therefore
  verifies nontriviality ITSELF after the projector call: `‖Ptilde_m − P'‖_op ≥ 0.15`
  via the §C5 midpoint opnorm (route (b): no projection-module change). A split
  returning the wrapper unit (`P' = Ptilde_m`) trips this RED (mutation-proven).
- **`aic_cstar_errreduce_unit` generalizes §C7 to a non-`1_n` unit (load-bearing at
  I5).** Stock `aic_errreduce` hardcodes `unit_def = ‖v(I_B) − 1_n‖` (FINDINGS §C7),
  which is ~1 and spuriously trips `AIC_ERRREDUCE_C0_CERT` for the Stage-2 maps into
  an S_P wrapper (unit `Ptilde`) and the Stage-3 intermediate merges (unit = the
  running `P_total ≠ 1_n`). `errreduce_unit` certifies against the SUPPLIED unit,
  reusing every `aic_errreduce` primitive. CAVEAT (review finding, follow-up bead):
  it relaxes `aic_dhom_approx`'s `unit_tol` to 2.0, which also loosens that routine's
  involution-symmetry assert (a genuine involution break in `[eps, 2.0]` is not caught
  HERE; the Newton contraction + max_steps guards DO still catch divergence).
- **c_0 is the MEASURED first-call constant (§D2, promoted design §G5).** Fixed from
  the first `aic_errreduce`; later calls gated `c0 ≤ 3·nominal + 5` (coarse last-resort;
  the precise universality canary is the test-level `c_hi/c_lo ≤ 2.5` ratio, the
  `.tex:484` check). The analytic c_0 stays open (aic-1bc).
- **OPEN coverage gaps (η=0-only at present; follow-up beads):** every η>0 fixture
  that stays in the funcalc basin is a SINGLE equivalence class, so (a) the Stage-3
  MULTI-CLASS merge and (b) the `errreduce_unit` Stage-3 running-`P_total≠1_n` branch
  run only at η=0 (machine-zero defects). A η>0 fixture with ≥2 distinct classes that
  stays in-basin is needed to exercise them — blocked by the next item.
- **The `aic_sgn` convergence wall — RESOLVED for the radius-floor case (bead aic-1vp,
  2026-05-31).** Diagnosed: NOT a basin gap and NOT precision. The Stage-2 oblique-
  wrapper `lem_extension` corner matrices were genuinely IN-BASIN (`‖I−X²‖_op ~ 0.1–0.6
  ≪ 1`, spec cleanly near {0,1}, gap ≈ 0.998) but carried a WIDE inherited arb radius
  (~2.7e-70 from the upstream corner star-product matmul chain). The Newton-Schulz /
  global-Newton iteration reached the involution by k≈6–7 (`‖Y²−I‖_mid ~ 1e-72`), but
  the inherited radius INFLATES ~2.75×/step through `acb_mat_sqr/mul`; once it dominates
  the iterate's deviation from the true sign (k≈7), the ball-arithmetic MIDPOINT itself
  drifts off the involution, so the midpoint step never falls below `tol = 2^-(prec-8)
  ≈ 2.2e-75` and the 100-iter cap fired. Higher prec made it WORSE (lower tol, further
  below the radius floor ~1e-72). **Fix (the house midpoint-iteration + a-posteriori-
  certificate pattern):** `aic_sgn_newton_schulz`/`_denman_beavers`/`_newton_global` now
  iterate on `mid(X)` (radius stripped once at entry, `aic_funcalc_int_to_midpoint`),
  then GATE the result by a shared a-posteriori certificate `aic_funcalc_int_certify_sign`
  (`‖Y²−I‖_F` AND `‖YX−XY‖_F < tol`, `tol = max(2^−(prec/2), in_rad)·8√n` — a SANITY
  BACKSTOP, prec-floor-dominated in practice (~2.9e-39 at prec=256 ≫ the ~1e-70 input
  radius); the `‖YX−XY‖` arm is computed on the radius-carrying X. The actual SOUNDNESS
  is NOT this tol magnitude — it is the away-from-0 basin/Gelfand precondition (on X,
  before midpointing: `ρ(I−X²)<1` keeps spec off the imaginary axis so `sgn(mid X)=
  sgn(X_true)`) + the `Y₀=mid(X)` seeding (Higham Thm 5.6, correct inertia; the cert
  alone cannot tell +sgn from −sgn). Hostile-review-corrected: earlier comments
  oversold this as "from the input radius" — the prec-floor dominates.) A genuinely
  out-of-basin X still FAILS LOUD (basin guard / Gelfand precondition / certificate). On a zero-radius tight
  input the strip is a no-op → BYTE-FOR-BYTE the prior behavior (existing in-basin
  cross-checks unchanged; T-global-3's `acb_mat_equal` dispatch check still green).
  **CONFIRMED:** `aic_cstar_build` on `make_mixconj(6,2,0.03)` now COMPLETES
  (num_blocks=1, d=[2], dim_B=4, iso_def=5.35e-2, iso_def/eta=2.57 — O(eta), v bijective,
  σ_min=0.974); the n=6,m=2 wall is GONE (also n=4,5 m=2 all complete). Regression +
  mutation-proof in `tests/test_funcalc.c` (`test_wide_radius`, `test_out_of_basin_failloud`).
  **WAS OPEN, now RESOLVED (the `m≥3` frontier, 2026-05-31, this entry's next item):**
  `make_mixconj(6,3,0.02)` got past the sgn wall but hit a DIFFERENT abort in
  `aic_dhom_approx` — the `lem_extension` `lem_approx` call. Resolved below.

- **The `m≥3` `lem_extension` frontier — RESOLVED (2026-05-31, bead aic-097 I4).** Once
  `m≥3` (a single equivalence class `M_m`) or any `dim S_{P,Q} ≥ 2` corner is reached,
  the Stage-2 `aic_cstar_lem_extension` `aic_dhom_approx(&h11v, eps_target=0, unit_tol=
  1e-10, …)` call aborted with `involution-symmetry 1.51e-1 > tol 1e-10 at step 1` (on
  `make_mixconj(6,3,0.02)`). Diagnosis (classification C, proven sound end-to-end): TWO
  independent caller-parameter bugs in the SAME call, neither in `aic_dhom_approx` itself.
  - **(i) The `{C_l}` corner basis is FROBENIUS-orthonormal, NOT operator-ON — the old
    spec §F-I4-2 claim is WRONG for `dim S_{P,Q} ≥ 2`.** `docs/research/lem_extension_spec.md
    §F-I4-2` claimed `‖C_l‖_op ≈ ‖C_l‖_F` and `<C_l,C_m>_Euc ≈ δ_lm + O(δ+ε)` (the
    coordinate norm ≈ the lem_PQ_Hilb Hilbert norm, `.tex:1146`). MEASURED on
    `make_mixconj(6,3)`: `‖C_0‖_op = 0.925` (NOT ≈ `‖C_0‖_F = 1`), and the lem_PQ_Hilb
    Euclidean Gram `G` differs from `I` by `‖G−I‖ ≈ 0.14` — an `O(1)` difference, not
    `O(δ+ε)`. The `O(δ+ε)`-equivalence the spec invoked holds only for `dim S_{P,Q} = 1`
    (where `S_{Q,Q} = ℂ Q̃`); for `dim S_{P,Q} ≥ 2` the off-diagonal corner products make
    `G ≠ I` by `O(1)`. (The `U_1` SVD extraction is still sound — the `‖u_0‖² > 0.5`
    guard tolerates the `O(1)` twist — but the per-coordinate-norm claim is false.)
  - **(ii) The Ha involution is G-TWISTED for `dim S_{P,Q} ≥ 2`.** `h11v = Ha^Q_{P,P}∘v`
    maps into `B(S_{P,Q}) ≅ M_n` whose lem_PQ_Hilb inner product has Gram `G ≠ I`. The
    map's TRUE Hilbert-space adjoint is `G^{-1}(·)^H G` (exact to 1e-16). But
    `aic_dhom_approx`'s involution-symmetry assert measures the PLAIN matrix-conjugate-
    transpose defect `‖v(E_ba)^H − v(E_ab)‖ ≈ 0.15` — an `O(1)` STRUCTURAL ARTIFACT of
    the G-twist, NOT a `*`-linearity break. It does NOT shrink with `δ₀` (it is `O(1)`,
    not `O(δ+ε)`), so an `8·δ₀`-style ceiling would also trip; the symmetrized Newton
    update cannot change it (frozen, harmless). FIX: pass `unit_tol = 2.0` ONLY at this
    caller (the `errreduce`/`errreduce_unit` callers already pass `1.0`/`2.0` for the
    SAME Ha-codomain reason). The downstream `aic_errreduce_is_bijective` certificate +
    the Newton-contraction guard still catch real divergence. Do NOT relax the assert in
    `aic_dhom_approx` (it is meaningful for the `eta=0` / genuine-codomain callers).
  - **(iii) `eps_target=0` mis-target (a separate real bug).** The call passed
    `eps_target=0.0`, treating the codomain `M_n` as exact-C*. But `h11v = Ha∘v` is only
    an `O(δ+ε)`-HOMOMORPHISM (`.tex:1157`: `Ha^Q_{P,P}` is an `O(δ+ε)`-hom), so its
    multiplicativity defect floors at `O(δ+ε) ≈ O(η)`, NEVER 0 → with `eps_target=0` the
    Newton iteration churns to `max_steps` and the cap-not-reached assert fires. FIX:
    `eps_target = 2·aic_ecstar_defect_assoc(A_parent)` (the parent's `ε = O(η)`), so
    `lem_approx` terminates cleanly at the `O(η)` floor. (`aic_ecstar_defect_assoc` is a
    `d³` star sweep, cheap for the small per-class `n`; works on the wrapper parent.)
  **CONFIRMED end-to-end:** with both fixed, `aic_cstar_build` on `make_mixconj(6,3,0.02)`
  COMPLETES (num_blocks=1, d=[3], dim_B=9, iso_def=5.58e-3, iso_def/η=0.43 — O(η), v
  bijective, σ_min=0.997); the `id≈0.15` involution defect stays frozen (harmless), the
  multiplicativity defect drops to its O(η²+ε) floor. T2b/T3 (`tests/test_cstar_build.c`)
  now sweep `m=3` (n=6,7); per-family c-ratios 1.76 (m=2)/2.05 (m=3), abs-max c=1.79 < 5
  (the `.tex:484` dimension-independence canary, extended to m≥3). Fix is in
  `src/aic_cstar_extension.c` (the `aic_dhom_approx` call ~line 195; no signature change
  — the `eps` is read from the parent in-place).
  **DEFERRED DEEPER FIX (bead aic-5aq):** make `aic_corner_extract` return an
  OPERATOR-orthonormal basis (so `G = I` and the plain conjugate-transpose IS the true
  involution, removing the G-twist entirely and letting this caller use a tight
  `unit_tol`). Not needed for correctness (the bijectivity certificate guards), but it
  would restore the involution invariant's meaning at this caller and let §F-I4-2's
  norm-equivalence claim hold by construction.

- **MULTI-CLASS Stage-3 merge at η>0 is now COVERED (2026-05-31, was an OPEN gap above).**
  The old "OPEN coverage gaps" item noted every in-basin η>0 fixture was a SINGLE class.
  `make_mixconj_blocks` (a `make_block_cond_exp(d,m)` base → `M_m ⊕ M_{d−m}`, conjugate-
  mixed; local to `tests/test_cstar_build.c`, T4) is a genuinely oblique η>0 channel with
  **2 equivalence classes** (confirmed: 4 one-dim projections → 2 classes at η>0), so the
  Stage-3 `cor_merge_sum` + the `errreduce_unit` running-`P_total≠1_n` branch RUN at η>0.
  T4 asserts: build completes, num_blocks=2 (sizes [2,2]/[2,3]/[2,2]), iso_def/η bounded
  (~0.004, O(η)), v bijective. CAVEAT (measured): this near-block-diagonal fixture has
  associativity defect `eps_assoc ≈ 2.2e-5 ≪ η ≈ 1.6e-2` (~700×), so T4 passes `η` as the
  build's `eps` (a faithful O(η) scale); the assoc defect would make the Stage-1
  `errreduce` C0 gate fire (`10·eps < the true O(η) inclusion defect`). No star tooth on
  this fixture (nearly block-diagonal → PLAIN-product c ≈ 0.28, does NOT fire the >20
  magnitude tooth; that discriminant stays on the single-block mixconj fixtures, T3).

### C12. `th_main_ext` (§10 opspace) cb-defect MUST be the OPERATOR-norm inclusion inf — the Frobenius coordinate σ_min ampliation is a "test that cannot fail"
- **Status:** CONFIRMED (opspace design review, bead aic-zwo; pre-implementation —
  flagged so the next session does NOT build a vacuous module). Cross-ref §C6, D3.
- **The trap.** th_main_ext's content is that the OPERATOR-norm cb-inclusion
  `a_n = inf ‖(1_{M_n}⊗v)(X)‖_op/‖X‖_op` stays `≥ 1−δ'` uniformly in n (the non-trivial
  prop_inc_ext induction `a_{2n}≥a_n/2`, `.tex:1493-1503`). The opspace design draft
  (`docs/research/opspace_design.md` §3.2) proposed to COMPUTE `a_n` as `σ_min` of the
  COORDINATE matrix `I_{n²}⊗M_1` (the Frobenius inclusion inf, the same quantity
  `aic_dhom_v_sigma_min` uses, §C6). But `σ_min(I_{n²}⊗M_1)=σ_min(M_1)` TRIVIALLY for
  ANY linear v (the singular values of `I_{n²}⊗M_1` are those of `M_1` repeated n²×),
  INDEPENDENT of v's quality. So "the cb-check reduces to the n=1 check" is vacuous: it
  is dimension-independent by pure linear algebra and tests NOTHING about the
  operator-norm ampliation. A Frobenius-σ_min "universality canary" for th_main_ext
  cannot fail (CLAUDE.md Rule 5).
- **The route.** Measure the OPERATOR-norm `a_n` (the genuine cb-inclusion): either the
  operator-norm worst-case (HOPM) over the ampliated unit ball `M_n⊗B` (the faithful
  worst-case deferred to **aic-0at** — th_main_ext is the first place it is LOAD-BEARING,
  since the Frobenius version is vacuous here), OR directly verify the prop_inc_ext
  doubling step `a_{2n}≥a_n/2` + uniformity in the operator norm for n=1..N. The η=0
  oracle (a_n=1 exact) + the universality canary (a_n dimension-independent in dim A AND
  the level n) must be on the OPERATOR-norm `a_n`. Recorded as the binding ORCHESTRATOR
  CORRECTION at the end of `opspace_design.md`. D3 stays BUILDABLE (the induction is
  rigorous, no hard stop) but the implementation is NON-trivial (operator-norm, depends
  on aic-0at or the induction-step verification — NOT the trivial `sigma_min` reuse).

---

## D. Open questions / escalations (unresolved)

### D1. Does a dimension-independent spectral gap (`Ω(1)`) always exist for `dim A>1`?
- **Status:** OPEN (bead **aic-3qv**). The projection construction certifies the gap
  **per-instance** and fail-loud-aborts on a degenerate spectrum; the *universal*
  a-priori guarantee is exactly what the paper needs Lefschetz for and is not proven
  constructively. The defect is `O(ε + ε/g)`, `=O(ε)` iff `g=Ω(1)`. Empirically the
  constant is dimension-independent (projection canary to d=9), but no proof.

### D2. The universal constant `c_0` (`cor_improvement`, `tex:1317`) is unstated
- **Status:** OPEN (the ANALYTIC `c_0` defers to **aic-1bc** research). The
  errreduce module (**aic-t81**, `src/aic_errreduce.c`, `cor_improvement`) is BUILT and
  returns the MEASURED `c_0` per instance (= max(inclusion-defect of `ṽ`)/ε), not the
  analytic constant. The paper says "there exist constants `ε_max, δ_max, c_0`"
  without numerical values; the analytic extraction must come from the `lem_approx`
  Newton analysis (the `δ_{s+1}≤C(δ_s²+ε)` constant and the `prop_delta_hominc`
  bounds). The `‖D‖_proj=m` consequence of A2 feeds into this (the `w'` bound is
  `O(mδ)`).
- **MEASURED (errreduce, `test_errreduce` T3/T4, 2026-05-30):** the empirical `c_0` is
  `≈ 2.0–2.7` and does NOT grow with dimension. Block-dim sweep B=M₂/M₃/M₄/M₅ (dim_B
  4/9/16/25, a 6.25× range): `c_0 = 2.714 / 2.218 / 2.069 / 1.962` (DECREASING with
  dim), ratio `c0_max/c0_min = 1.384 < 1.6` (T4(A) threshold tightened in F2 to catch
  even sublinear/sqrt-dim growth — a `c_0 *= √dim_B` injection drives the ratio to
  1.807 > 1.6, RED; mutation-proven 2026-05-30, independently of T3). If T4(A)
  RED-fires the `tex:484` failure mode is back. Block-count sweep m=1,2,3 of M₂:
  `c_0 = 2.714 / 1.238 / 1.564` — does NOT grow with m, so the `‖D‖_proj=m → O(mδ)` w'
  bound does NOT manifest as growth in the achieved `c_0` (the per-step quadratic
  Newton contraction reaches the O(ε) floor regardless of m). These match the dhom
  layer's `C=defect/ε` exactly (errreduce just wraps `aic_dhom_approx`). The
  lem_approx termination floor is set to `AIC_ERRREDUCE_EPS_FLOOR=4`×ε (NOT bare ε:
  the defect cannot beat the algebra's intrinsic O(ε) obstruction, so driving toward
  bare ε stalls and bounces — fail-loud in `aic_dhom_approx`'s contraction guard). The
  certification ceiling `AIC_ERRREDUCE_C0_CERT=10` (tightened 50→10 in F3) is a
  generous fail-loud guard for "is `ṽ` an O(ε)-inclusion", NOT the analytic `c_0`
  (the worst achieved max-defect over the corpus is `≈ 2.71`×ε at T4(A) M₂, a ~3.7×
  margin under 10×ε; the machine-floor cases T1/T2 clear it with ≥30× margin).

### D3. cb-norm truncation `N` (shard F, `tex:1447-1561`)
- **Status:** OPEN (bead **aic-2jd**). "for all n" in the cb-norm must be truncated;
  conjecture `n≤dim A` suffices. Needs proof or a certified-N procedure before
  `opspace`/`th_main_ext`. (Out of scope for plain th_main.)

### D4. `th_factorization` is an outline (`tex:2742`, shard H)
- **Status:** OPEN (bead **aic-1sk**). The labelled proof block ends without executing
  the CP-ization (Steps 4–5); the constructions in `tex:2771-2899` are prose. The
  composite `O(η)` constant is unspecified. Reconstruct the closure before
  `factorize` is trustworthy.

### D5. The certified degenerate-eigenvalue wall (`aic-w4o.1`)
- **Status:** OPEN (the recurring deferral). FLINT's certified eig needs a *simple*
  spectrum; the project's spectra are degenerate (projections 0/1, `⊕B(L_j)`). So
  certified **rank/subspace extraction** (corner `dim S`, idemp/assoc subspaces,
  projection gap *enclosure*) uses the LAPACK double path now and defers
  certification to aic-w4o.1 (Rump cluster enclosures / eig-free Cholesky). The
  *defects* that avoid full eig are already arb-certified.

---

## How to use / extend this file

- When a module's research or hostile review finds a `.tex` typo, a non-constructive
  step, a subtle pitfall, or an open question: **add an entry here** (right
  category, with `.tex` line + status + where it surfaced), and reference
  `paper/FINDINGS.md §Xn` from the source comment where it bites.
- When an OPEN item is resolved, update its status and point to the resolving
  commit/module (don't delete it — the history is the value).
- This file is part of the **read order** (CLAUDE.md / AGENTS.md): skim it before
  touching a new region of the paper, so you inherit the known issues.
