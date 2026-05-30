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
- **Status:** CONFIRMED recurring; worked around; root-cause OPEN (bead **aic-qgs**;
  neighbour of aic-2yo).
- `aic_mat_opnorm` forms the Gram `M†M` and routes through a relative-Hermiticity
  check that **false-fails** when an off-diagonal Gram entry has midpoint ~0 but a
  matmul-**accumulated** arb radius exceeding the absolute floor — even though `M`
  is tight and genuinely Hermitian. Hit by corner (`lem_alpha` γ, Ha-map `G−I`,
  `lem_PQR`), projection (basin assert), dhom (everywhere). **Workaround:** use the
  certified mid+radius upper bound `aic_corner_gamma_opnorm_ub` (`‖mid(M)‖_op +
  ‖rad(M)‖_F`) for op-norms of near-zero-off-diagonal / star-defect matrices.
  Until aic-qgs is fixed, new modules must use that helper.

---

## D. Open questions / escalations (unresolved)

### D1. Does a dimension-independent spectral gap (`Ω(1)`) always exist for `dim A>1`?
- **Status:** OPEN (bead **aic-3qv**). The projection construction certifies the gap
  **per-instance** and fail-loud-aborts on a degenerate spectrum; the *universal*
  a-priori guarantee is exactly what the paper needs Lefschetz for and is not proven
  constructively. The defect is `O(ε + ε/g)`, `=O(ε)` iff `g=Ω(1)`. Empirically the
  constant is dimension-independent (projection canary to d=9), but no proof.

### D2. The universal constant `c_0` (`cor_improvement`, `tex:1317`) is unstated
- **Status:** OPEN (beads **aic-t81** errreduce / **aic-1bc** research). The paper says
  "there exist constants `ε_max, δ_max, c_0`" without numerical values; they must be
  extracted from the `lem_approx` Newton analysis (the `δ_{s+1}≤C(δ_s²+ε)` constant
  and the `prop_delta_hominc` bounds). The `‖D‖_proj=m` consequence of A2 feeds into
  this (the `w'` bound is `O(mδ)`).

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
