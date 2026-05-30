# Research + decision: globally-convergent non-normal matrix-sign (bead aic-8hz)

Two Sonnet research legs (web), 2026-05-30. Goal: relax funcalc's `sgn` precondition
from the operator-norm Newton-Schulz basin `‖I−X²‖_op<1` to the SPECTRAL
no-purely-imaginary-eigenvalue condition, so the full `η<1/4` regime of the
NON-NORMAL superoperator `S = 2S_Φ−1` is reachable (assoc_ecsa, aic-92f I1).

## Findings (legs 1 + 2)

1. **Global-convergence theorem (Higham, *Functions of Matrices*, SIAM 2008,
   Ch.5, Thm 5.6).** The Newton sign iteration `Y_{k+1}=½(Y_k+Y_k^{-1})`, `Y_0=X`
   converges globally and quadratically to `sgn(X)` for ANY `X` with no
   purely-imaginary eigenvalue. Hypothesis is purely SPECTRAL (not a norm bound).
   Newton-Schulz `½Y(3I−Y²)` is only LOCALLY convergent (basin `‖I−X²‖<1`).
   This is exactly our `aic_sgn_denman_beavers` iteration — but it is currently
   gated by the op-norm basin assert, artificially.
2. **Our regime is favorable.** Eigenvalues `λ` of `X=2S_Φ−1` satisfy
   `λ²−1 = 4μ(μ−1)`, `μ∈spec(S_Φ)`, so `ρ(I−X²)=4ρ(S_Φ²−S_Φ)≤4η<1` and the
   spectrum sits in the disk `|z−1|<1 ⊂ {Re>0}` — already near `±1`, so UNSCALED
   Newton converges in a handful of steps. Scaling (Kenney–Laub spectral-ratio
   `μ_k=(‖Y_k^{-1}‖/‖Y_k‖)^{1/2}`, Byers determinantal `|det Y_k|^{-1/n}`) only
   earns its keep near `η→1/4`; KEPT as a Pareto candidate, NOT built yet (Law 4).
3. **Zolotarev / QDWH (Nakatsukasa–Freund, SIAM Rev. 2016) are INAPPLICABLE** —
   they require Hermitian/normal input (polar = sign only for Hermitian). The
   non-normal superoperator rules them out.
4. **Non-normal stability (Bai–Demmel, SIMAX 1998).** Newton is asymptotically
   backward-stable; non-normality inflates the sign condition number `~1/sep`,
   `sep≥2√(1−4η)` here (diverges as `η→1/4`, as expected). Mitigation in arb:
   raise prec; the certified ball radius reports if more is needed.

## Decision (orchestrator)

**Precondition certified RIGOROUSLY and EIG-FREE via the Gelfand bound.** The
spectral condition `ρ(I−X²)<1` is certified by `‖(I−X²)^k‖_F^{1/k} < 1` for
increasing `k` (Gelfand: `‖M^k‖^{1/k}→ρ(M)`, each power a rigorous upper bound;
`k=1` is exactly the current op-norm/Frobenius check). This converges DOWN to the
true spectral radius, reaches the full `ρ<1` regime, needs no eigensolver, and so
**sidesteps the open degenerate-eig debt `aic-w4o.1`**. Cap `k` and fail loud.

**A-posteriori certificate.** After convergence: certified `‖Y²−I‖` (near an
involution) + `‖XY−YX‖` (commutes with X) + per-step non-singular-inverse guard.
With the certified precondition and `Y_0=X` basin invariance (Higham Thm 5.6),
these pin `Y=sgn(X)` (rules out `−sgn`, wrong-inertia involutions). Fail loud on
any straddling ball (Rule 4).

**Iteration = unscaled Newton `½(Y+Y^{-1})`** (the DB core), reused. Scaling
deferred to a future audition if a benchmark shows too many iterations near 1/4.

**Wiring:** `aic_sgn` auto-dispatches — fast Newton-Schulz when the op-basin is
certified, else the global Newton variant (strict improvement: previously
aborted). `aic_prop_P` relaxed to the Gelfand spectral precondition so
`aic_theta`/`prop_P`/`aic_assoc_regularize` reach the full non-normal `η<1/4`.

**Teeth (exact oracle).** Synthetic non-normal `X=[[a,c],[0,−b]]`, `a,b≈1`,
large `c`: `ρ(I−X²)=max(|1−a²|,|1−b²|)` small but `‖I−X²‖_op≈|c(a−b)|` ≥1, and
`sgn(X)=[[1, 2c/(a+b)],[0,−1]]` EXACTLY (2×2 poly: `sgn=(2X−(a−b)I)/(a+b)`). Old
op-basin `sgn` ABORTS; global Newton returns the closed form to ~1e-70. Plus the
η=0 oracle (`sgn(2P−I)=2P−I` machine-zero) and in-basin NS-vs-Newton agreement.

## References
- Higham, *Functions of Matrices*, SIAM 2008, Ch.5 (Thm 5.6, uniqueness via
  spectral half-planes, gap-Lipschitz bound).
- Kenney & Laub, SIMAX 13(3) 1992 (scaling); Byers, LAA 85 1987 (det scaling).
- Bai & Demmel, SIMAX 19(1) 1998 (sign for invariant subspaces; non-normal cond).
- Nakatsukasa & Freund, SIAM Rev. 58(3) 2016 (Zolotarev — Hermitian only).
- Rump, *Acta Numerica* 19 (2010), "Verification methods" (a-posteriori certs).
- FLINT/arb `acb_mat` eigenvalue enclosures (Johansson 2018) — noted, not needed
  given the Gelfand route.
