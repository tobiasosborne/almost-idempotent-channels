# Matrix-free / low-rank / certified Φ̃ — deep-research findings (aic-tg8)

Deep-research sweep (2026-06-04, 101 agents, 19 primary sources, 25 claims
adversarially verified 3-vote → 21 confirmed / 4 killed). Question: can the
O(d⁶) dense-superoperator computations be made matrix-free / low-rank scaling
with Kraus rank r, **while preserving the paper's certified O(ε) bound** (FLINT
arb/acb)? Bead `aic-tg8` (epic `aic-g0z`). See PROFILING.md §9.

## Bottom line

The dense O(d⁶) work **can** be made matrix-free → **O(#matvecs · r·d³)** (each
Φ-application is O(r·d³)), **but no off-the-shelf method delivers a certified
O(ε) bound for free.** Every fast matrix-free method surveyed (FEAST,
Krylov-Ritz sign, Zolotarev/QDWH) ships only a *residual-threshold or sgn²=I
consistency heuristic* — not a validated/interval enclosure. So:

- **Fast (uncertified) path**: any of them serves directly → feeds the fast-path
  epic `aic-dsj`.
- **Certified path**: requires a **separate rigorous a-posteriori enclosure
  layer** (verified-eigensolver / Davis–Kahan / resolvent-gap, in ball
  arithmetic) *layered on the fast output*. This is exactly the project's
  existing "fast guess at low prec + cheap rigorous a-posteriori certify"
  idiom — that two-layer split **is** the architecture. The enclosure layer is
  the open research problem (below).

## The three matrix-free primitives for Φ̃ = (I + sgn(2S_Φ − 1))/2

**#1 — Rational Arnoldi for sgn(2S_Φ−1)·v** (the only one with an *a-priori*
non-normal bound).
- *Cost*: one shifted linear solve per iteration; #matvecs·O(r·d³), no dense O(N³).
- *Certification*: the **Crouzeix numerical-range a-priori bound**
  `‖f(A)b − f_m^RA‖ ≤ 2C‖b‖ min_{r_m} ‖f − r_m‖_Σ`, universal `C ≤ 11.08`
  (`C=1` Hermitian), valid for general non-normal A with f analytic on a compact
  Σ ⊇ W(A). **The only A-dependent term is a scalar best-rational-approximation
  factor.** *But*: Güttel states it is **explicitly NOT sharp for highly
  non-normal A** ("the gap can be quite large even for symmetric matrices") — so
  it likely cannot certify a tight O(ε) target alone; back it with an a-posteriori
  enclosure. Source: [Güttel, *Rational Krylov approximation of matrix functions*,
  Cor. 3.4 / Eq. 9](https://ora.ox.ac.uk/objects/uuid:cd6a7402-1d50-4055-ac28-9b53eace04e1).

**#2 — Multishift / nested Krylov-Ritz matrix sign with deflation** (the
demonstrated matrix-free non-normal sign workhorse).
- `sgn(z)=z/√(z²)=sgn(Re z)` (cut on the negative real axis); projected Ritz sign
  `sgn(H_k)` via quadratic Roberts–Higham Newton (7–10 iters). Rational sign
  approximations → **multishift**: `sgn(A)v` = Σ shifted resolvents solved in one
  shift-invariant Krylov subspace.
- *Cost driver*: inner `sgn(H_k)` is **O(k³)**; for accuracy k grows large.
  Published fix: **nested** projection onto an ℓ≪k subspace → O(ℓ³)+O(kℓ).
- *Certification*: ships only the **sgn²=I heuristic** (`ε=|x̃−x|/2|x|`, validated
  empirically only) — NOT rigorous. Demonstrated at scale (overlap Dirac, lattices
  to 10⁴). Sources: [Bloch et al. 0912.4457](https://arxiv.org/pdf/0912.4457),
  [Bloch et al. 0910.2927](https://arxiv.org/pdf/0910.2927).

**#3 — FEAST / contour-integral spectral projector** (matrix-free invariant
subspace).
- Applies `P=(1/2πi)∮(μI−A)⁻¹dμ` to a random block (randomized range finding);
  per iteration = `n_e` shifted solves with `m0` RHS, never forming P. Non-normal
  extension: separate **left + right** projectors (≈2× solves).
- *Convergence*: **linear**, rate `|ρ_a(λ_{m0+1})/ρ_a(λ_i)|`; **requires m0 ≥ m**
  (eigenvalue count in the contour), recommend `m0=2m`. The near-1 cluster of S_Φ
  is exactly the slow regime — cluster size sets the cost.
- *Certification*: only a relative-residual threshold + a **probabilistic**
  eigenvalue-count (exponential miscount decay under a uniform model) — NOT a
  validated enclosure. Sources:
  [Kestyn/Polizzi/Tang 1506.04463](https://arxiv.org/pdf/1506.04463),
  [Ye/Xia/Chan, SIMAX 2017](https://www.math.purdue.edu/~xiaj/work/feasteig.pdf).

## Dominant failure mode (all three)

Spectral fragility of the **non-normal** S_Φ with eigenvalues clustered near 0/1.
The shifted operator `2S_Φ−1` has spectrum near ±1, but eigenvalues drifting
toward the **sign discontinuity** (S_Φ eigenvalues near ½ → imaginary axis) force
a low-order polynomial to make a steep −1→+1 jump → impractically large Krylov
dimension. Non-normal ⇒ **long recurrences** (no short recurrence, Faber–Manteuffel)
⇒ restarts. Standard cure: **deflation** of the near-discontinuity eigenvectors
(compute their sign contribution exactly). This is the CLAUDE.md spectral-fragility
callout, sharpened.

## Traps caught by the adversarial verification (do NOT step on these)

1. **HSS rank ≠ Kraus rank.** The Ye/Xia "O(r·n²) all-eigenpairs" speedup uses
   `r` = HSS / off-diagonal numerical rank, **wholly unrelated to the Choi/Kraus
   rank**. It does NOT transfer to S_Φ unless S_Φ is *separately* shown to be
   HSS-rank-structured (open question). A real trap.
2. **Zolo-pd (QDWH) is Hermitian-only for the speedup.** The 2-iteration rate-17
   Zolotarev algorithm computes the **polar** decomposition; `sign = unitary polar
   factor` holds **only for Hermitian A**, so the headline speedup does NOT give
   `sign(S_Φ)` for the non-normal superoperator. (The projector *formula*
   ½(I±sgn) is general; only the fast Zolo route is Hermitian-restricted.)
3. **A claimed "rigorous a-priori+a-posteriori rational-Krylov bound"
   (arXiv:2311.02701) was REFUTED** (0-3) — do not rely on it.
4. **Szehr's resolvent bound carries a dimension-dependent prefactor**
   (`cot(π/4n)` classical, `n²` quantum) — *conflicts with the project's
   dimension-INDEPENDENT universal-constant requirement* (`th_main`, tex:484). It
   certifies perturbation sensitivity, not a projector, and its dimension scaling
   must be reconciled. [Szehr, J. Spectral Theory 2014](https://arxiv.org/pdf/1305.7208).

## Coverage gaps (honest — what the sweep did NOT settle)

- **The CRUX is unsolved by any single source.** No turnkey validated/interval
  enclosure (Rump/Miyajima/Behnke, INTLAB verifyeig) for a Krylov/FEAST result on
  a **non-normal** operator was found as a *verified claim*. The verified-eigensolver
  sources were fetched (`21M1451440`, `S0024379513007969`, `s11075-017-0332-y`,
  `2008.04140`, `120894683`, `1307.7219`) but produced no surviving top-25 claim —
  applicability/cost on S_Φ is **unestablished**. This is open question #1.
- **Candidates D (low-rank Choi / Kraus-rank fixed-point-algebra) and E
  (randomized low-rank) produced NO surviving verified claims.** Sources were
  fetched (`1212.3839`, `1803.00109` for fixed-point/DFS; `2211.04676`,
  `1904.06277` for randomized) but didn't survive the verification budget. **The
  report explicitly flags candidate D as possibly "the more natural route to
  Kraus-rank scaling" than spectral-projector methods** — because the spectral
  route (A/B/C) is matrix-free but its #matvecs can be large (clustered spectrum)
  and is NOT obviously *Kraus-rank* scaling, whereas the rank-r Choi /
  Choi–Effros / Lindblad fixed-point structure (which the project ALREADY computes
  at η=0 via Wedderburn) might give genuine r-scaling. **This is the highest-value
  follow-up.**

## Open questions (carried forward)

1. Cost & tightness of a **rigorous a-posteriori enclosure** (verified-eigensolver
   / Davis–Kahan / resolvent-gap, in arb/acb) layered on a fast Krylov/FEAST
   invariant-subspace guess for the non-normal d²×d² S_Φ — can it be done in the
   target O(r·d³) extra cost and meet a *dimension-independent* O(ε)? **(The crux.)**
2. Are there **Kraus-rank-scaling fixed-point-algebra / DFS / attractor-subspace**
   algorithms exploiting the rank-r Choi directly (η=0 Choi–Effros / Lindblad
   structure)? Candidate D — likely the more natural r-scaling route; warrants its
   own focused sweep.
3. Is S_Φ (or 2S_Φ−1) actually HSS / off-diagonal-rank-structured in any basis, so
   the Ye/Xia contour-eigensolver complexity could transfer?
4. Does the near-fixed-point spectral gap `min|1−λ_i|` of an η-idempotent Φ stay
   bounded below by a function of **η alone** (dimension-independent)? Needed for
   any resolvent-gap certificate.

## Strategic implications for the beads

- **`aic-xsv` (Krylov Φ̃) is viable but must be a TWO-LAYER design**: (i) fast
  matrix-free sign — rational Arnoldi (#1) or multishift Krylov-Ritz+deflation (#2);
  (ii) a SEPARATE rigorous a-posteriori enclosure (open question #1). The
  enclosure, not the fast layer, is the bottleneck for the certified path.
- The fast layer alone **directly serves `aic-dsj`** (the uncertified fast path).
- **`aic-d7a` (Choi-rank exploitation) is likely the more promising Kraus-rank
  route** and is under-researched — a focused follow-up deep-research on candidate
  D (fixed-point-algebra / DFS / Lindblad-attractor algorithms + their
  certification) is the recommended next research step.
- Do NOT pursue Zolo-pd for sign(S_Φ) (Hermitian-only) or assume the Ye/Xia
  r-scaling (HSS rank, not Kraus rank).

### Sources (all primary unless noted)
Güttel rational Krylov (ORA Oxford); Bloch et al. [0912.4457], [0910.2927];
Kestyn/Polizzi/Tang [1506.04463]; Ye/Xia/Chan (Purdue, SIMAX 2017);
Nakatsukasa/Freund [10.1137/140990334]; Szehr [1305.7208]. Verified-eigensolver +
candidate-D/E sources listed in the gaps section (fetched, not yet verified).
