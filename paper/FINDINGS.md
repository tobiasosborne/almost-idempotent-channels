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

- **CORRECTION (2026-05-31): the T2b "per-family c-ratio ≤ 2.5" metric above was FLAWED
  — measures fixture-GEOMETRY SPREAD, not dim-GROWTH; replaced by a robust bounded +
  no-trend canary.** The claim five paragraphs up ("per-family c-ratios 1.76 (m=2)/2.05
  (m=3), abs-max c=1.79 < 5") was computed from a TWO-POINT within-family ratio
  `c_hi/c_lo` over a tiny fixture set (m=2: n=4,5; m=3: n=6,7). That metric measures the
  SPREAD of `c=iso_def/η` across different ambient `n`, NOT whether `c` GROWS with `n`
  (the genuine `.tex:484` / §D2 failure mode). On the committed HEAD it read RED
  DETERMINISTICALLY: the m=3 ratio is 6.9035 > 2.5 (a hard-geometry outlier, n=7 c=1.87,
  sitting over a favorable n=6 c=0.47 — pure placement, no dim content). An extended,
  PRECISION-STABLE sweep (orchestrator diagnostic; prec 256 == prec 512 byte-for-byte)
  established the truth: `c` is BOUNDED in `[0.25, 3.27]` and does **NOT** grow with `n`.
  m=2 (B=M₂, t=0.03) over n=4..10: `c = 1.01, 1.79, 2.56, 3.27, 0.39, 0.49, 0.61`. m=3
  (B=M₃, t=0.02) over n=6..10: `c = 0.47, 1.87, 0.25, 0.74, 0.79`. BOTH families PEAK at
  the n=7 Heisenberg-Weyl/compression geometry (a hard corner) then CRASH at n=8 — the
  LARGEST-`n` points have the SMALLEST `c`. All isos are HEALTHY (bijective,
  σ_min ∈ [0.963, 0.999]; the low-`c` points have the HIGHEST σ_min — low `c` = a GOOD
  iso with small defect, never a degenerate collapse). th_main is SOUND (bounded,
  dimension-independent constant); the variation is fixture-geometry noise.
  **T2b now uses a ROBUST canary** (`tests/test_cstar_build.c`, FINDINGS §D2): the
  EXPANDED sweep (m=2 n=4..10, m=3 n=6..10) fed to (i) absolute boundedness abs-max
  `c < 5.0` (measured 3.27, 1.53× margin) AND (ii) a no-upward-trend halves-ratio
  `mean(c|hi-half n)/mean(c|lo-half n) ≤ 1.25` (measured m=2 = 0.28, m=3 = 0.65; the
  least-squares slope β is reported too: -0.21 / -0.05, NEGATIVE). The halves aggregate
  dilutes the n=7 spike; a genuine c=O(n) law survives the aggregation (halves-ratio
  → 1.80 for m=2, mutation-proven: a `c_flat·(n/n_min)` injection trips the trend arm,
  and the literal `c·(n/n_min)` trips the absolute arm at abs-max 5.73 > 5.0). The old
  geometry-fragile within-family ratio is GONE.

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
- **Status:** CONFIRMED + RESOLVED for O1 (opspace operator-norm ampliation machinery
  IMPLEMENTED + fix-pass hardened, bead aic-zwo; `include/aic_opspace.h`,
  `src/aic_opspace_ampliate.c` (HOPM kernel) + `src/aic_opspace_apply.c` (the
  `1_{M_n}⊗f` block-ampliation primitives, Rule 10 split) + `src/aic_opspace_map.c`
  (opmap builders) + `src/aic_opspace_entry.c` (public entry points + the I1/I2
  factorize adds, Rule 10 split) + `src/aic_opspace_cert.c`; `tests/test_opspace.c`,
  89 checks). The vacuous Frobenius-σ_min route was NOT built. The header
  PRECISION-POSTURE note was corrected (finding 4): the HOPM kernel is PURE DOUBLE
  (`out` = arb_set_d, a zero-radius wrap), so T-cross tests the prec of the M_1 / M_1⁻¹
  ASSEMBLY (~1e-15), a coarse gate consistent with `aic_dhom_v_sigma_min` — NOT a
  genuine arb-vs-double ALGORITHM cross-check. Cross-ref §C6, D3.
- **O1 IMPLEMENTED (the faithful operator-norm route, 2026-05-31).** `aic_opspace_*`
  measures the OPERATOR-norm ampliated max-stretch via a scale-invariant HOPM over the
  op-norm unit ball of `M_n⊗B` (forward `‖v_n‖_op`) / `M_n⊗A` (inverse `‖v_n^{-1}‖_op`,
  the eps-C* subspace → polar-then-PROJECT accept guard, the ecstar pattern). a_n =
  1/‖v_n^{-1}‖_op. Smith truncation: forward N_max = N = v->n, inverse n_B = Σ_l d_l
  (D3). **HONESTY:** HOPM is a LOWER bound on the op-norm stretch (good witness, maybe
  suboptimal), so O1 certifies the η=0 complete isometry, the prop_inc_ext doubling,
  and the universality canary — NOT `‖v‖_cb ≤ 1+O(eps)` (the SDP UPPER bound is the
  SEPARATE O2 increment, a bead to be filed). MEASURED: η=0 oracle (block_cond_exp,
  noiseless_subsystem) → `‖v_n‖_op == ‖v_n^{-1}‖_op == 1` to 1e-12 for n=1,2,4,...,n_B
  (a COMPLETE ISOMETRY, exactly §C12 (a)); η>0 mixconj → a_n = 0.97–0.998 (O(η)-close
  to 1), a_cb_flat (a_n across the ampliation level n) = 1.00004–1.0005 (genuinely
  level-INDEPENDENT, the prop_inc_ext/Smith content); doubling a_{2n} ≥ a_n/2 holds.
  double vs arb@53 agree to ~1e-15.
- **§C12 NON-VACUITY tooth (measured, sharp).** On mixconj(6,3,0.02), scaling one
  vE[0] by 1.6 INFLATES the operator-norm forward stretch `‖v_n‖_op` from 1.000006 to
  1.600141 (Δ ≈ +0.60 at n=1 AND n=2 — the op-norm HOPM catches the inflated witness
  direction), while the Frobenius `σ_min(M_1)` is UNCHANGED (1e-14 = machine noise) —
  AND `σ_min(I_{n²}⊗M_1) = σ_min(M_1)` is level-independent for ANY v. So a
  Frobenius-σ_min ampliation canary is provably BLIND to this operator-norm effect (the
  exact "test that cannot fail"), and the op-norm route is faithful. This is the
  mutation-proven discriminant (`test_opspace.c` test_c12_nonvacuity).
- **STRUCTURAL-AMPLIATION "test that cannot fail" CAUGHT + CLOSED (O1 hostile review,
  2026-05-31, fix-pass aic-zwo).** The first opspace test suite (72 checks) had a
  §C12-class hole: a hostile review crippled `apply_fn` (the `1_{M_n}⊗v` block map,
  `src/aic_opspace_apply.c`) to process ONLY the (0,0) block — i.e. computed the
  ampliation AS IF it zeroed every off-diagonal block — and ALL 72 checks stayed
  GREEN. The entire mathematical content of th_main_ext (the block ampliation
  structure) was tested by nothing with teeth, because the real fixtures give a
  level-INDEPENDENT `a_n ≈ 1.0005`, so a structural bug that preserves the stretch
  RATIO is invisible to the stretch/flatness checks. CLOSED by the
  **structural-ampliation cross-check tooth** (`test_opspace.c` test_struct_ampliation,
  + the §C12 n=2 arm tied to it): it pins the production ampliation
  `aic_opspace_ampliate_forward` (driving the SAME `aic_opmap_apply_fn` the HOPM uses)
  against an INDEPENDENT explicit Kronecker block assembly at n=2,3 with genuinely
  nonzero off-diagonal blocks (reference block (k,l) = `aic_dhom_v_apply(v, X_{kl})`, a
  path with no opmap coordinate matrix). Measured agreement: max|ampliate − block_ref|
  = 0 (η=0) / 1.1e-15 (η>0), off-diagonal reference content ≈ 1.0–1.4 (genuinely
  nonzero). **Mutation-proven:** the (0,0)-block-only cripple (MUTATION-D) now makes
  the structural tooth RED (max-diff 1.165 ≫ 1e-9) AND the §C12 n=2 arm RED (level-2
  struct-diff 1.808, the dropped inflated off-diagonal `1.6·vE[0]` block); reverted.
  Other O1 fix-pass items: the §C12 n=2 stretch arm threshold tightened 0.3 → 0.5
  toward the genuine +0.60 move; the AXIS-D confounded KAPPA "halves-ratio(hi/lo dim)"
  arm DROPPED (the fixture array mixes dim_A=4/9 and the d=7 N-spike confounds a
  position-split — a genuine c∝dim_A law gives halves-ratio 0.206, does NOT trip
  KAPPA=1.6; the rigorous dim-A TREND sweep lives in cstar_build T2b on the SAME v,
  opspace inherits it), keeping the mutation-proven abs-max `c < C_ABS` boundedness arm.
  Check count 72 → 89.
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
- **§C12.O2-PIN — the O2 SDP convention was pinned EMPIRICALLY, correcting BOTH the
  design doc and the research leg (bead aic-pjr, 2026-05-31).** `‖v‖_cb = ‖v*‖_⋄` is
  computed by feeding the adjoint's Convention-A Choi into the Watrous diamond-norm SDP.
  TWO conventions were uncertain and the prior analyses DISAGREED: the design doc
  (`opspace_o2_design.md` §2.4) said normalization `2/N`; the Sonnet research leg derived
  `2/n_B`. **Both were WRONG.** A convention-sweep probe (`tools/probe_o2_pin2.jl` +
  `probe_o2_diag2.jl`) pinned it against an INDEPENDENT closed-form truth — an asymmetric
  CP map `Ψ(Y)=A†YA: M_3→M_2` with `‖Ψ‖_⋄ = σ_max(A)²` — plus a complete-isometry oracle
  (`‖v‖_cb=1` exactly). **PINNED GOLDEN RULE (design §0.5):** to get `‖f‖_⋄` for
  `f: M_in→M_out`, build `J = choi_convA(f, in, out)` (INPUT-major,
  `J[s·out+i,t·out+j]=f(E_st)[i,j]`), feed the SDP with `(d_maj=in, d_min=out)` →
  **raw optval = `‖f‖_⋄` EXACTLY, normalization FACTOR = 1** (no `2/n` at all — the `2/n`
  was an artifact of the self-map's `P+Q=I` primal form; the rectangular density-form
  primal needs no factor). Dual traces the **MINOR/OUTPUT** factor (`tr_sys=2`); primal
  density on `:major` (= input; `:major`≠`:minor` on an asymmetric map, so placement is
  load-bearing). Build the adjoint's Choi DIRECTLY (`v*(E_ab)=Σ_i conj(vE[i][a,b])E_i`),
  NOT `transpose(J(v))` (the full transpose keeps the wrong `[n_B,N]` block layout — it is
  `v*` in the OUTPUT-major convention, needing dims `(n_B,N)`, which mis-grouped the
  factors and gave the W-dependent garbage 1.76/2.0 that first surfaced the bug). The
  lesson restates the project ethos (design §6.5): normalization/direction are PINNED not
  derived — neither the design's nor the LLM-research-leg's derivation was trustworthy;
  only the independent-oracle measurement was. (Standing rule: Sonnet for survey,
  Opus/empirics for derivation.)
- **§C12.O2 — the certified cb-norm UPPER bound is now built (bead aic-pjr, O2.4/6/7,
  2026-05-31).** `aic_cbnorm_certify_rect_upper` (`src/aic_cbnorm_certify_rect.c`) is the
  RECTANGULAR generalization of the self-map dual restorer: it certifies `‖f‖_⋄` for
  `f: M_in→M_out` from a committed Watrous MIN-dual point, via the PINNED convention
  (design §0.5: normalization FACTOR 1, dual `tr_sys=2` = the MINOR/output factor
  `partial_trace_right(.,d_maj,d_min)`, shift `eps·d_min`). `aic_opspace_certify_cb_upper`
  (`src/aic_opspace_o2.c`) assembles `J(v*)`/`J((v⁻¹)*)`, certifies, and asserts the
  HOPM(O1)≤SDP(O2) bracket. MEASURED (`tests/test_opspace_o2.c`): η=0 oracles
  (block_cond_exp 4×4, noiseless_subsystem 6×3) give `hi=[1,1]` fwd+inv; mixconj(6,2,0.03)
  fwd `‖v‖_cb=1.0019683734` (HOPM 1.001431 ≤ this), inv `‖v⁻¹‖_cb=1.5353598357`
  (HOPM 1/a_cb=1.018942 ≤ this). The §C12 non-vacuity is SHARP at O2: the cb
  `‖v⁻¹‖_cb=1.535` vs the vacuous Frobenius `1/σ_min(M_1)=1.027` (gap 0.51). Restoration
  PSD-defect `eps`: 0 (block_cond_exp), ≤9.9e-13 (noiseless), ≤8.9e-11 (mixconj) — NO
  precision wall (design §6.4 cleared).
- **§C12.O2 SUBTLETY (load-bearing, the arb-radius vs Hermiticity-tol wall).** The rect
  certifier's `J(v*)`/`J((v⁻¹)*)` are ASSEMBLED in arb over a `cstar_build` `v`, so their
  entries carry an ACCUMULATED radius (~1e-71 at prec=256, ≫ a single-rounding ulp). The
  dual block `[[Y0,−J],[−J†,Y1]]` is then fed to `aic_mat_herm_max_eig`, which asserts
  Hermiticity RIGOROUSLY at the prec-tight tol `2^-(prec-8)` (~1e-75 at prec=256, §C5/
  aic-2yo). The off-diagonal `−J`/`−J†` are INDEPENDENT balls whose asymmetry ball ~ J's
  radius EXCEEDS that tol → the genuinely-Hermitian block is REJECTED (the self-map
  `int_upper` never hit this: its committed Choi is a zero-radius double). FIX (in
  `aic_cbnorm_certify_rect_upper`): collapse the block to its MIDPOINT then symmetrize
  `(blk+blk†)/2` before the defect eig — the midpoint differs from the true block by ≤ the
  radius (1e-71) in op-norm, far below the 1e-4 tightness tol, so `eps` stays a rigorous
  defect bound for the (midpoint) feasible point — exactly the self-map's zero-radius-
  committed-J posture. The committed self-map path is UNAFFECTED (`test_certify` 34 checks
  stay green). Note: this is why a feasible-point seed assembled in arb is NOT directly
  interchangeable with a committed double seed for the rigorous eig — the radius must be
  collapsed first.

### C13. `factorize` F3 (lem_RC / UCP decode Υ, `tex:2840-2899`) — the F-ancilla ordering deviation, the structurally-vacuous ξ_j tooth, and the m≥2∧η>0 coverage debt
- **Status:** CONFIRMED (factorize F3, bead aic-tff; `src/aic_factorize_upsilon{,2,3}.c`,
  `tests/test_factorize.c` T5/T6 + the P-trace/P-cent/D5 pins). Hostile-reviewed SHIP.

- **(a) THE F-ANCILLA ORDERING — code deviates from the `.tex` Kronecker side (Law-1 deviation).**
  `.tex:2862` writes `L_j = Σ_s p_{js}(Δ(U_{js}†)⊗1_F) V W_j†(U_{js}⊗ξ_j)` and `.tex:2869`
  writes `Υ'_j(X)=L_j†(Φ(X)⊗1_F)L_j` — the ancilla `F` factor on the RIGHT. **The code
  writes those `⊗1_F` factors F-LEFT: `1_F⊗Δ(U_{js}†)` and `1_F⊗Φ(X)`.** This is REQUIRED,
  not a typo-fix: `aic_ucp_kraus_to_stinespring` builds `V` ancilla-MAJOR
  (`V[a·dim_K+i,j]=K_a[i,j]`, `aic_ucp.h:14,90`), i.e. `Φ(X)=V†(1_F⊗X)V`, so EVERY `⊗1_F`
  must be F-LEFT to match V's layout. The paper's `(·)⊗1_F` is convention-relative to its
  own `V: H→H⊗F` (H-left); our V is F-left, so the sides swap. **η=0 / r=1 is BLIND** (F is
  1-dim, `1_F⊗M == M⊗1_F`); the **r>1 D5 pin** (`aic_factorize_upsilon_d5_pin`, exercised at
  T6 r=6/8) is the discriminant: F-LEFT gives `‖Υ'Δ−1_B‖ ≈ 4.4e-2 = O(η)`, the WRONG F-RIGHT
  gives `≈ 7.6e-1 = O(1)` (~17× separation). A production mutation to F-RIGHT aborts
  fail-loud in the Υ' unitalization basin. The source docstrings (upsilon2.c:10-18,
  upsilon3.c:8-24) cite THIS entry.

- **(b) The ξ_j right-vs-left singular-vector tooth is VACUOUS IN PRINCIPLE (an unreachable
  "test that cannot fail").** `.tex:2859` picks ξ_j with `‖C_jξ_j‖≈1`; the code uses the top
  RIGHT singular vector (so `‖C_jξ_j‖=σ_max(C_j)`). But `C_j = (1/d_{L_j})Tr_{L_j}(R_j)` with
  `R_j = Σ_s p_s(U_s†⊗1)W_jW_j†(U_s⊗1)` is a positively-weighted sum of congruences of the
  PSD `W_jW_j†`, hence **structurally Hermitian PSD for ANY W_j** (measured `‖C_j−C_j†‖_F ~
  1e-86`, `|⟨u₀,v₀⟩|=1` exactly on every fixture). So left ≡ right and a right→left mutation
  changes NOTHING (verified: 70/70 stay green). The right-vector choice is paper-correct but
  a non-normal C_j is NOT constructible from this exact R_j — the tooth is correctly dormant,
  not a coverage debt to repay.

- **(c) The `‖R_j−1_{L_j}⊗C_j‖` centrality tooth is EXACT (not O(η)) but VACUOUS at η=0;
  the m≥2 ∧ η>0 regime is unexercised (coverage debt → F4).** lem_RC(i) `R_j=1_{L_j}⊗C_j`
  is EXACT regardless of η (the full Heisenberg-Weyl Pauli set is an EXACT unitary 1-design,
  so the twirl projects exactly onto the `{U⊗1}`-commutant); the O(η) of lem_RC lives entirely
  in `‖C_j‖`=σ_max via lem_RC(ii) (`σ_max ≈ 1−5.4e-3 ~ 1−O(η)`). The check HAS teeth — a
  drop-Pauli / U→I mutation drives `‖R_j−1⊗C_j‖` to `5.1e-1` (O(1) RED) — but ONLY at η>0:
  at η=0 the exact-idempotent W_j already yields a central `W_jW_j†`, so the twirl is not
  load-bearing and the residual is 0 for ANY design (verified: U→I gives residual 0 at m=2
  η=0 vs 5.1e-1 at m=1 η>0). Since the η>0 fixtures (`make_mixconj`) are all SINGLE-block
  (m=1), the **per-block-Pauli-vs-whole-B-join distinction (D3) at m≥2 ∧ η>0 is untested**.
  **F4 must add a multi-block η>0 fixture** (`make_mixconj_blocks` in `tests/test_cstar_build.c`
  is ready; honor its §C11 caveat — pass `eta`, NOT the ~700×-smaller `eps_assoc`, as the
  build's eps scale). This is the F2→F3→F4 coverage debt (the F2 review flagged the same).

---

### C14. `factorize` F2/F3/F4 — Δ′ is only **O(η)-completely-positive**, not exactly CP (the multi-block η>0 Kraus extraction needs PSD-cone projection); and the §A.2 `1_B`=`P_B` conditional-expectation fix
- **Status:** CONFIRMED (factorize F4, bead aic-tff; surfaced by the FIRST m≥2 ∧ η>0 test, the §C13(c) coverage debt). Root cause verified against `.tex:2786-2796`. The constructive route (PSD-cone projection in the Choi→Kraus extraction) is BUILDABLE-MODULO — **not a wall**.

- **(a) THE PAPER PRESENTS Δ′ AS EXACTLY CP, BUT IT IS CP ONLY TO O(η).** The 1-design CP-ization
  `Δ′(X)=Σ_s p_s Φ(Δ̃(XU_s†)Δ̃(U_s))` (`.tex:2786-2789`) is argued completely positive at
  `.tex:2791-2796` by writing, for `X=Y†Y`, `Δ′_n(Y†Y)=Σ_s p_s Φ_n(Δ̃_n(Y†(I⊗U_s†))·Δ̃_n((I⊗U_s)Y))≥0`
  (manifestly PSD: each term is `Φ(Z_s†Z_s)` with `Z_s=Δ̃((I⊗U_s)Y)`, `Φ` CP). **That manifest form
  requires the EXACT term-equality** `Δ̃(Y†YU_s†)Δ̃(U_s)=Δ̃(Y†U_s†)Δ̃(U_sY)`, which holds iff `Δ̃` is
  an EXACT unital homomorphism (then both sides `=h(Y†Y)` since `h(U_s†)h(U_s)=h(I_B)=1`). But
  `Δ̃=v` is only an *extended O(η)-isomorphism* (`v(XY)=v(X)⋆v(Y)+O(η)`, the Choi–Effros star
  `≠` the ordinary product — §C2). So the equality holds only to **O(η)**, and **Δ′ (hence the
  unitalized UCP `Δ=Δ′(I)^{-1/2}Δ′(·)Δ′(I)^{-1/2}`) is completely positive only to O(η).** At
  m≥2 ∧ η>0 the per-block Choi of `Δ` carries a small NEGATIVE eigenvalue of order **O(η²)**
  (MEASURED `−2.5e-6` at η=2.3e-2, `−8.7e-8` at η=4e-3 → ratio 28.7 ≈ η²-ratio 33; both arb
  `herm_max_eig` and LAPACK agree, ≫ machine noise). **SINGLE-block Δ′ is exactly CP** (T4 minEig
  `+9.4e-6`): for one block the per-block Pauli twirl makes `Δ̃` effectively multiplicative on the
  relevant combination, so the O(η) hom-defect cancels — which is WHY this never surfaced before F4.

- **The abort it caused.** `aic_ucp_choi_to_kraus_latd` (`src/aic_ucp_latd.c:79-89`) fail-loud-aborts
  on any eigenvalue `< −thr` (`thr ~ 1e-9`). The W_j Stinespring extraction
  (`aic_factorize_upsilon.c`, `.tex:2831-2838` Choi_Delta) and the F4.1 dual read-off
  (`aic_factorize_dual.c`, Dec=Δ*/Enc=Υ* via the full-UCP Choi) both hit the O(η²) negative at
  multi-block η>0 and aborted — blocking the §C13(c) debt AND F4.2's dim-independence canary
  (which sweeps `make_mixconj_blocks` at η>0).

- **CONSTRUCTIVE ROUTE (not a bandaid, Rule 3; not a wall).** A tolerance-aware PSD-cone-projecting
  Choi→Kraus extraction `aic_ucp_choi_to_kraus_latd_tol(.., neg_tol, *clipped_neg_out)`: clip
  eigenvalues in `(−neg_tol, keep_thr]` to 0 (the cone defect), accumulate the clipped negative
  mass, but **STILL fail-loud abort on `λ ≤ −neg_tol`** (a genuine non-CP — an O(1) or O(η)-scale
  negative). The strict `aic_ucp_choi_to_kraus_latd` is preserved BYTE-FOR-BYTE by delegating to
  `_tol` with `neg_tol=thr` (existing callers unchanged). The factorize sites pass a `neg_tol` that
  ADMITS the O(η²) defect but ABORTS an O(1) one, and `aic_factorize_upsilon_build` ASSERTS the
  accumulated clipped mass is O(η²)-small (fail loud otherwise — the genuine-bug guard). The W_j
  become the Kraus of the NEAREST genuinely-CP map (`Δ` projected onto the CP cone), within O(η²)
  of `Δ`. **Why this is sound, not a fudge:** th_factorization only requires `‖Δ−Δ̃‖_cb≤O(η)`
  (`.tex:2801`), so even an O(η) cone-projection stays within spec — and the measured defect is the
  far-smaller O(η²). lem_RC's R_j/C_j and Υ′ then proceed; their O(η) bounds absorb the O(η²)
  perturbation.

- **(b) THE §A.2 AMBIENT-`M_{n_B}` CHOI OF `ΥΔ−1_B`: the `−1_B` term is the CONDITIONAL EXPECTATION
  `P_B`, not the full identity.** `factorize_f4_design.md` §A.2/§B.1 says "subtract `1_B(E_pq)=E_pq`
  for every ambient unit `E_pq`." That is WRONG for off-block-diagonal `E_pq`: `Δ` reads only the
  block-diagonal coordinates of its `M_{n_B}` input (so `ΥΔ(E_pq)=0` for off-block `E_pq`), so
  subtracting `E_pq` leaves `−E_pq` on every off-block column → `‖J_UpsDel‖_F=√(#off-block units)
  =2.83` even at η=0 (a "test that cannot pass", §C12 class). The correct `1_B` extension to ambient
  `M_{n_B}` is the **conditional expectation `P_B` onto the block-diagonal** (`P_B(E_pq)=E_pq` for
  in-block, `=0` for off-block), making `ΥΔ−P_B` zero on off-block input (both terms drop it) and
  `=ΥΔ−id_B` on in-block input. Then the ambient cb-norm `=` the in-`B` cb-norm (the off-block
  columns of `J_UpsDel` are zero — this IS why route (i) ≡ route (ii); §D probe `offblk=0.000`).
  Mutation-proven (revert `P_B`→full-`I` → `‖J_UpsDel‖_F=2.83` RED). Source: `aic_factorize_verify.c`.

### C15. `make_mixconj` CANNOT produce an out-of-regime (`η ≥ 1/4`) channel — the "t≈0.45 is out-of-regime" claim is FALSE for the SPECTRAL basin; the real fail-loud fixture is unitary conjugation by a reflection
- **Status:** CONFIRMED + RESOLVED (bead aic-xo0, the fail-loud entry guard). Corrects a load-bearing
  premise in the bead, in CLAUDE.md "Probe/sweep hygiene", and in the aic-xo0 NOTES.
- **The trap.** The bead + CLAUDE.md cite `make_mixconj(5,3,t≈0.45)` as an OUT-OF-REGIME channel that
  hangs the pipeline, via the OP-NORM eta-proxy (`~6.5 t`, so t=0.45 ⇒ proxy ~2.9 ≫ 1/4). But the
  proxy is `‖S_Φ²−S_Φ‖_op`, NOT the spectral radius `rho(S_Φ²−S_Φ)` that `aic_prop_P` actually
  certifies post-aic-8hz. `make_mixconj` is a CONVEX MIX of an exactly-idempotent compression and its
  UNITARY conjugate — BOTH summands are spectrally idempotent — so the mix stays SPECTRALLY
  near-idempotent at EVERY `t`. **Measured** `rho(4(S²−S))` (Gelfand ub) over `t ∈ {0, .02, .06, .45,
  .7, .85, .95, 1.0}`: `{0, .095, .27, .66, .56, .62, .23, ~0}` — ALL `< 1` (`rho(S²−S) < 1/4`), so
  `make_mixconj` is IN regime for the whole `t`-range (it is exactly idempotent at both t=0 and t=1).
  Feeding any `make_mixconj` to `aic_assoc_regularize` SUCCEEDS; it never trips the basin guard. So
  any "evil η-sweep" over `make_mixconj` `t` is exercising IN-regime inputs (the original session
  "hang" was the OLD op-norm guard + radius inflation on a large-`n` non-normal `S`, a path the
  aic-8hz spectral relaxation removed — the current Gelfand guard at `k_max=32` aborts in ~0.2 s, it
  does NOT hang).
- **The genuine out-of-regime fixture.** The `*`-automorphism `Φ(X)=U X U^†` with `U=diag(1,−1)` (a
  Pauli-Z reflection): superop `S = U ⊗ conj(U)` has eigenvalues `u_i conj(u_j) ∈ {+1,−1}`, and the
  `λ=−1` eigenvalues give `(S²−S)` the eigenvalue `(−1)²−(−1)=2`, so `rho(S²−S)=2 ≥ 1/4` and
  `rho(4(S²−S))=8` (measured Gelfand ub ~8.09, NOT certified `<1`). It is unital, CP, trace-
  preserving — a VALID UCP map, just not almost-idempotent. This is the fail-loud fixture in
  `tests/test_xo0_failloud.c`.
- **Where the guard lives.** Added at the `aic_assoc_regularize` PUBLIC entry
  (`src/aic_assoc_regularize.c`, `aic_assoc_int_assert_eta_basin`): certifies `rho(4(S²−S)) < 1` via
  the EXISTING eig-free Gelfand certifier `aic_funcalc_int_gelfand_rho` (NOT a cheap op-norm gate,
  which would FALSE-reject the §C/U6 spectral-relaxation regime where `‖S²−S‖_op ≥ 1/4` but
  `rho < 1/4`, aic-8hz), and aborts with a CHANNEL-LEVEL message naming `aic_assoc_regularize`, the
  `η < 1/4` hypothesis (`.tex:2168/2176`), and the prop_P basin (`.tex:524-525`). The deeper
  `aic_prop_P` Gelfand guard ALREADY aborted fast, but with a generic "P²−P / prop_P" message that did
  not attribute the failure to the input channel; the entry guard is the API-boundary attribution.
  The cstar_build master loop has NO unbounded loop (Stage-1 capped at `dim_A+1`, errreduce at
  `max_steps`, the only `while` is a bounded union-find); the funcalc sgn iterations cap at 100/200.

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
- **The universality canary's ROBUST form (2026-05-31, `test_cstar_build` T2b):** the
  dimension-independence check on the th_main constant `c=iso_def/η` is **bounded
  abs-max + no-upward-trend**, NOT the geometry-fragile within-family `c_hi/c_lo` ratio
  (which measures fixture-geometry SPREAD across ambient `n`, not dim-GROWTH — it read
  RED 6.90 on a hard-geometry n=7 outlier with no `.tex:484` content; see §C11's
  2026-05-31 correction). The robust canary: over an EXPANDED ambient-n sweep (m=2 B=M₂
  n=4..10, m=3 B=M₃ n=6..10), assert (i) abs-max `c < 5.0` AND (ii) the halves-ratio
  `mean(c|hi-half n)/mean(c|lo-half n) ≤ 1.25` (the halves aggregate dilutes a single
  hard-geometry spike; a genuine c=O(n) law survives it). Measured: abs-max 3.27,
  halves-ratio 0.28 (m=2) / 0.65 (m=3), slope -0.21 / -0.05 (NEGATIVE). Mutation-proven:
  a `c_flat·(n/n_min)` (genuine c=O(n)) injection drives the halves-ratio to 1.80 > 1.25
  (trend arm RED), and the literal `c·(n/n_min)` injection drives abs-max to 5.73 > 5.0
  (absolute arm RED). The errreduce-layer canary (T4(A) above) keeps the dim-sweep ratio
  form because its `c_0` is a smooth, monotone-in-dim quantity there (no geometry spike).

### D3. cb-norm truncation `N` (shard F, `tex:1447-1561`)
- **Status:** RESOLVED (2026-05-31, bead **aic-2jd** / **aic-zwo**); O1 IMPLEMENTED.
  `src/aic_opspace_*` truncates the operator-norm ampliated stretch at the Smith
  levels: forward `N_max = N = v->n`, inverse `n_B = Σ_l d_l` (the doubling sweep
  `1,2,4,...,n_B`). The op-norm route is the FAITHFUL one (§C12); the Frobenius
  σ_min reuse §3.2 of `opspace_design.md` proposed was NOT built (vacuous). The
  conjecture `n≤dim A` is superseded by a THEOREM. Two independent research legs
  (`docs/research/opspace_paper_leg.md`, `opspace_web_leg.md`) converged on **Smith's
  lemma** (R.R. Smith, "Completely bounded maps between C*-algebras," *J. London Math.
  Soc.* (2) **27** (1983) 157–166; textbook: Pisier, *Introduction to Operator Space
  Theory* (2003) **Proposition 1.12**, p. 26; Watrous, *TQI* (2018) **Thm 3.46 + Cor
  3.47** via adjoint duality): for ANY linear map `u: E → M_N` from an operator space
  `E` into the matrix algebra `M_N` (no CP / C*-algebra hypothesis),
  `‖u‖_cb = ‖1_{M_N}⊗u‖_op` — the cb-norm is ATTAINED at ampliation level `N` (the
  codomain ambient dimension), and is not increased by any larger ampliation.
- **Application to our `v: B → A`.** `A ⊆ M_N` (ambient `N = v.n`), so
  `‖v‖_cb = ‖1_{M_N}⊗v‖_op` — finite truncation `N_max = N`. The inverse
  `v⁻¹: A → B ⊆ M_{n_B}` (`n_B = Σ_l d_l`) has `‖v⁻¹‖_cb = ‖1_{M_{n_B}}⊗v⁻¹‖_op`,
  truncation `N_max = n_B`. The cb-inclusion lower bound `a_n = inf ‖v_n(X)‖_op/‖X‖_op`
  satisfies `a_n = 1/‖v_n⁻¹‖_op` (v bijective), so its all-`n` value `= 1/‖v⁻¹‖_cb`
  is attained at `n = n_B`. **No "for all n" verification problem remains** — and this
  is INDEPENDENT of the prop_inc_ext induction (`.tex:1486-1505`), which already proves
  the uniform bound analytically; Smith gives the exact finite certification level.
- **Caveat (carry into implementation).** The norms above are the OPERATOR norm
  (§C12), NOT the Frobenius coordinate σ_min — the design-doc `σ_min(I_{n²}⊗M_1)`
  route is vacuous. `‖v_N‖_op` is a genuine operator-norm worst-case (max-stretch of a
  linear map), computed either by the Watrous cb-norm SDP (certified UPPER bound,
  reusing `src/sdp.jl` + `aic_cbnorm_certify` via `‖v‖_cb = ‖v*‖_⋄`; the O2 increment)
  or a HOPM operator-norm worst-case (LOWER bound; the O1 structural core). The
  CERTIFIED cb-bound needs the SDP upper bound; HOPM is the cross-check lower bound
  (HOPM ≤ SDP, the aic-0at agreement check). See `opspace_web_leg.md` §1–§2.

### D4. `th_factorization` is an outline (`tex:2742`, shard H)
- **Status:** ASSESSED — **BUILDABLE-MODULO, not a hard wall** (2026-05-31, bead
  **aic-1sk** research; `docs/research/factorize_d4_research.md`). The labelled proof
  block ends without executing the CP-ization (Steps 4–5); the constructions in
  `tex:2771-2899` are prose. But every object in Steps 4–5 is an explicit finite-dim
  matrix expression, so the prose is a *closure* gap, not a constructivity wall.
  Per-step verdict (research §2): (D4-a) the unitary 1-design CP-ization of `Δ̃`
  (`tex:2771-2801`) is BUILDABLE-MODULO (trivial) — the per-block Heisenberg-Weyl
  design is already built (`aic_dhom_pauli`, the genuine per-block `S_{jk}`, NOT the
  embedded sum); (D4-b) `lem_RC`'s `C_j = d_{L_j}^{−1} Tr_{L_j}(R_j)` + `ξ_j` (top SVD)
  (`tex:2840-2864`) is BUILDABLE-MODULO (partial trace + SVD, all finite matrix ops,
  `lem_RC` is itself constructive); the closure chain (`tex:2895`) is a sequence of
  `O(η)`-triangle steps each backed by an already-stated bound. **The ONLY open item is
  the composite `O(η)` CONSTANT `C` (research §5, escalation 5)** — a CERTIFICATION gap,
  NOT a constructivity wall: the algorithm runs and produces `Δ,Υ,B`; `C = ‖ΔΥ−Φ‖_cb/η`
  is composed of the prop_P `c1`, the iso `c0` (itself analytically OPEN via
  cor_improvement `c_0`, §D2/bead `aic-1bc`), and the CP-ization `c2,c3`, none
  multiplied out in the paper. **Recommended posture (project standard):** measure `C`
  per instance + assert it is bounded & dimension-independent (the §D2 robust canary),
  handled per-instance + canary like `c_0` — file the analytic `C` as a research bead
  chained after `aic-1bc`. So `factorize` is reachable after opspace O1 (η=0 oracle +
  η>0 measured path) / O2 (the certified η>0 `tilde_DelUps` upper bound). The
  increment skeleton (F1–F4) is in the research doc §6.
- **opspace interface adds landed (2026-05-31, fix-pass aic-zwo).** The research §4
  highest-value adds are now exposed so `factorize` plugs in cleanly: (I1) `v⁻¹` is a
  PUBLIC builder `aic_opspace_build_vinv` (`include/aic_opspace.h`) returning `v⁻¹(B_k)
  ∈ B` as `dim_A` operators — the SAME M_1⁻¹ the inverse stretch certifies, so
  `Υ̃ = v⁻¹∘Φ̃` uses the certified-and-used inverse (round-trip test `v⁻¹(v(E_i)) = E_i`
  to 7.8e-17, `test_opspace` test_vinv_roundtrip); (I2) the inverse Smith level
  `n_B = Σ_l d_l` is a first-class `aic_opspace_result.n_B` field. Still gated on O2 for
  the certified `‖Δ̃‖_cb,‖Υ̃‖_cb ≤ 1+O(η)` UPPER bound (research §4(b); O1's `cb_forward`
  is a LOWER bound).

### D5. The certified degenerate-eigenvalue wall (`aic-w4o.1`)
- **Status:** OPEN (the recurring deferral). FLINT's certified eig needs a *simple*
  spectrum; the project's spectra are degenerate (projections 0/1, `⊕B(L_j)`). So
  certified **rank/subspace extraction** (corner `dim S`, idemp/assoc subspaces,
  projection gap *enclosure*) uses the LAPACK double path now and defers
  certification to aic-w4o.1 (Rump cluster enclosures / eig-free Cholesky). The
  *defects* that avoid full eig are already arb-certified.

### D6. `factorize` F4.2 — the diamond-norm DUAL SDP stalls (SLOW_PROGRESS) at n≥6 in Convex.jl
- **Status:** OPEN, DEFERRED to v0.2 (bead `aic-bag`). Surfaced 2026-05-31 building the
  F4.2 dimension-independence canary; the headline (`th_factorization`, `aic-tff`) is
  CLOSED on F4.1 (the constructive content; see below).
- The F4.2 canary measures `C = ‖ΔΥ−Φ‖_cb / η` over a dim sweep via the SQUARE self-map
  Watrous diamond-norm SDP on `J_DelUps`/`J_UpsDel` (`tools/gen_fixtures_factorize_f4.jl`,
  `src/aic_factorize_shim.c`). The **PRIMAL** (`diamond_norm_watrous_primal`,
  `maximize Re tr(J†X)`) converges **OPTIMAL** and its optimum **IS** the diamond/cb-norm
  (Watrous) — e.g. `mixconj_blocks(6,2,0.03)`: primal `‖ΔΥ−Φ‖_cb = 0.012955` OPTIMAL,
  reproduced across runs. The **DUAL** (`diamond_norm_dual`, `min` over
  `opnorm(partialtrace(·))` epigraphs) **STALLS at SLOW_PROGRESS** for n≥6, returning a
  loose `0.0237` it cannot drive down to the primal optimum — independent of MOSEK
  tolerance: tight `1e-12`/`1e-14` → SLOW_PROGRESS + ~100s + ~20 GB at n=6 and **OOM at
  47 GB for n=7**; relaxed `1e-9` (no iter cap) → still SLOW_PROGRESS `0.0237`. The
  strong-duality poison guard CORRECTLY refuses to emit (no fake fixture).
- **Root cause:** the QETLAB `opnorm(partialtrace(·))` dual epigraph is ill-conditioned
  for Convex.jl's DCP→MOI translation at these dims; the `1e-14` MU_RED target is below
  MOSEK's numerical floor (the stall). NOT a math error — the primal value is the
  faithful cb-norm; only the redundant tight-upper *certificate* won't converge.
- **Decision (user, 2026-05-31):** land the headline on **F4.1** (committed, green:
  certified `Δ,Υ,B`; η=0-exact oracle; per-instance rigorous `O(η)` eig-free bound; the
  four dual channels). Defer the **rigorous two-sided diamond-norm certification** + the
  faithful canary to **v0.2** (`aic-bag`). The headline's constructive content is closed;
  the canary is the universality-certification refinement (the project already certifies
  dimension-independence at th_main `cstar_build` T2b and th_main_ext `opspace`).
- **v0.2 routes (in `aic-bag`):** (a) reformulate the dual without the `opnorm` epigraph;
  (b) move the self-map diamond SDP from Convex.jl to **direct JuMP + MOSEK** (lower memory,
  controllable conic form) — recommended; (c) primal-only canary (primal optimum IS the
  cb-norm) + eig-free T8 as the rigorous per-instance upper. Scaffolding committed:
  `src/aic_factorize_shim.{c,h}` (green), `tools/gen_fixtures_factorize_f4.jl` (eager-flush,
  `F4_ONLY` filter, relaxed-tol override, GC).

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
