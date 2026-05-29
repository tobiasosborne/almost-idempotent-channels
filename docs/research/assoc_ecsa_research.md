# Research capture: assoc_ecsa Increment 1 — Φ̃ = θ(2Φ−1) regularization

**Date:** 2026-05-29
**Context:** `almost-idempotent-channels` — constructive finite-dim implementation
of Kitaev arXiv:2405.02434, C + FLINT/arb + Julia-ccall.
**Task:** STEP 1 of `th_almost_idemp` (`.tex:2162-2237`): regularize an
η-idempotent UCP map `Φ` to an EXACTLY idempotent superoperator `Φ̃ = θ(2Φ−1)`,
the matrix `prop_P` (`.tex:524-533`) applied to `Φ` at the superoperator level.

---

## Route decision (the audition slate)

The regularization is `Φ̃ = θ(2Φ−1) = ½(1 + sgn(2Φ−1)) = ½(1 + (2Φ−1)(1−4(Φ−Φ²))^{−1/2})`
(`.tex:2171-2174`). In finite dim a UCP self-map on `M_n` is an `n²×n²` matrix
`S_Φ`, so this is matrix functional calculus on the (NON-NORMAL) matrix `2S_Φ−1`.

Candidate routes for the matrix `sgn`/`θ` (all eig-free; the project's `funcalc`
audition already chose Newton-Schulz as the default `aic_sgn`):

- **A1 — `aic_prop_P`(S_Φ)** = `θ(2S−1)` via `aic_sgn` (Newton-Schulz iteration
  `Y_{k+1} = ½ Y_k(3I − Y_k²)`). DEFAULT. Inverse-free; converges for
  `‖I−X²‖_op < 1`; no normality assumption (matrix multiply/add only).
- **A2 — binomial series** via `aic_funcalc_xpow`: `(1−4(Φ−Φ²))^{−1/2}` as the
  Taylor series `Σ C(−½,k)(M−I)^k`, `M = (2S−1)²`. Independent of A1 (Newton vs
  Taylor). USED AS THE NON-NORMAL CROSS-CHECK (test T6), not the default.
- **A3 — Denman-Beavers** `Y_{k+1} = ½(Y_k + Y_k^{−1})` (one inverse/step). Already
  in `funcalc`; not separately exercised here (A1 is the bench winner; A2 is the
  independent oracle).
- **B1 — eigendecomposition of `2S−1`** then scalar `θ` on eigenvalues. REJECTED:
  `2S−1` is non-normal (fragile spectrum, `.tex:540-544`), and the certified arb
  eig path needs a SIMPLE spectrum (the η=0 superoperator is maximally degenerate).
  This is exactly the route `funcalc` rejected for the Hermitian case.

**Decision:** A1 is the default build (`aic_assoc_regularize`); A2 is the
independent algorithmic cross-check that ALSO fills funcalc's never-before-tested
non-normal `sgn` gap (T6). A1 and A2 must agree to ~1e-25 at prec=256 on a
genuinely non-normal in-basin channel.

---

## Key citations (research legs)

### Matrix sign / Newton-Schulz convergence

- **Higham, *Functions of Matrices: Theory and Computation* (SIAM 2008), Ch. 5.**
  Thm 5.6: the (full) Newton sign iteration `X_{k+1} = ½(X_k + X_k^{−1})` converges
  globally and quadratically to `sgn(A)` for any `A` with no pure-imaginary
  eigenvalue. The Newton-Schulz variant `X_{k+1} = ½ X_k(3I − X_k²)` is
  inverse-free but converges ONLY in the basin `‖I − A²‖ < 1` (§5.3). For
  `A = 2Φ−1` the basin is `‖(2Φ−1)² − I‖ = 4‖Φ−Φ²‖ < 1`, i.e. `η < ¼` — exactly
  the paper's hypothesis (`.tex:2176`). No normality is needed: both iterations
  are defined by polynomials/rational functions in `A`, so they commute with `A`
  and act on its Jordan form; non-normality only affects the transient, not the
  fixed point. This justifies using Newton-Schulz on the non-normal `S_Φ`.

- **Newton-Schulz Frobenius basin.** The Frobenius-norm sufficient condition
  `‖I − A²‖_F < 1` is the standard conservative basin (Higham §5.3, and the
  inverse-free-iteration literature). `funcalc`'s `aic_sgn` asserts the
  operator-norm form `‖I − A²‖_op < 1` (tighter, still rigorous). assoc keeps test
  channels with `4‖S²−S‖_op < 1` and asserts it live (T6); out-of-basin is bead
  aic-8hz.

### Superoperator / vec conventions

- **Wood, Biamonte, Cory, "Tensor networks and graphical calculus for open quantum
  systems," arXiv:1111.6950** (Quantum Inf. Comput. 2015). Canonical reference for
  the row-vectorization / column-vectorization conventions and the resulting
  superoperator (natural representation) `S = Σ_a K̄_a ⊗ K_a` (column-vec,
  Schrödinger) vs the row-vec form. We use the HEISENBERG action
  `Φ(X) = Σ_a K_a^† X K_a` with ROW-MAJOR vec `vec_r(X)[i·n+j] = X[i,j]` (matching
  the project's `aic_idemp`/`aic_ecstar`), which gives `S_Φ = Σ_a K_a^† ⊗ K_a^T`
  (left factor `K_a^†`, right factor `K_a^T`, left-major Kronecker). Derived by
  index AND oracle-tested two ways (T1 column-image vs T2 Kronecker), per the
  project hard rule (the vec convention is the #1 historical bug; `idemp`
  deliberately avoided forming the superoperator).

### Why a Hermitian-only purification route does NOT apply

- **McWeeny purification** (`P_{k+1} = 3P_k² − 2P_k³`) and the broader
  density-matrix-purification family drive a near-idempotent toward an exact
  projector. They are designed for HERMITIAN (symmetric) `P` (electronic-structure
  density matrices) and rely on the real eigenvalue ordering on `[0,1]`. The
  superoperator `S_Φ` is NON-Hermitian / non-normal, so McWeeny is INAPPLICABLE as
  a regularization route here — the correct general tool is `sgn`/`θ` (Higham
  Ch. 5), which `prop_P` (`.tex:524-533`) uses. Recorded so a future agent does not
  reach for purification by analogy with the "idempotent" keyword.

---

## Findings surfaced during implementation

1. **A2 sign is load-bearing.** `M = 1 − 4(Φ−Φ²)`, so `D = 4(S − S²)`, NOT
   `4(S² − S)`. The η=0 case (`D=0`) is blind to the sign; the η>0 NON-NORMAL T6
   channel pins it (a flipped sign gives `sgn² ≠ I` ⇒ A1/A2 disagree by `O(0.1)`).
   This is precisely why a non-normal eta>0 cross-check was mandated.

2. **`aic_mat_herm_max_eig` near-zero fragility (root-cause fixed).** For a
   near-zero Hermitian `H` (the Gram of `S²−S` for an EXACT idempotent, entries
   `~1e-31`), `acb_mat_eig_global_enclosure` returns a NON-FINITE radius ⇒
   `aic_mat_opnorm` = `nan±inf` ⇒ `prop_P`'s basin guard spuriously aborts on the
   cleanest ground-truth input. Fixed in `src/aic_mat_spectral.c`: when the
   certified radius is non-finite, fall back to the rigorous eig-free bound
   `|λ| ≤ ‖H‖_F`. No existing test regressed.

3. **Inexact "exact" idempotents.** `trace_replace`, `noiseless_subsystem` use
   `1/√k` from a `double`, so `Φ` is idempotent only to `~2.7e-16` (not machine-
   zero at prec=256); the η=0 oracle (T3) therefore uses a defect-aware tolerance
   `100·ubound(‖S²−S‖_F) + 1e-22` = the `prop_P` bound with the input's true `δ`.
   Integer-Kraus channels (`block_cond_exp`, `identity`, `dephasing`,
   `compress_idemp`) ARE exact (`δ ~ 1e-30`).
