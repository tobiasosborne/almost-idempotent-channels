# Tight certified-arb cb-ball for `eta = ||Phi^2-Phi||_cb` (aic-m24 increment 3)

Implementation contract for the TIGHT two-sided certified ball `[eta_lo, eta_hi]`
on the eta-idempotence defect. TIB-grounded against the actual literature
(2026-05-29); supersedes the earlier from-memory reconstruction in the research
notes. The always-valid eig-free bracket `[||J||_F/n, 2||J||_F]`
(`src/aic_cbnorm.c`, increment 1) is the fallback this dispatches to near `eta=0`.

## Sources (grounded, fetched via TIB/web)
- **Jansson, Chaykin, Keil**, "Rigorous error bounds for the optimal value in
  semidefinite programming," SIAM J. Numer. Anal. 46(1):180-200, 2008,
  DOI 10.1137/050622870. Thm 3.2 (rigorous lower bound from approximate dual),
  Thm 4.1 (rigorous upper bound from approximate primal).
- **Watrous**, "Simpler semidefinite programs for completely bounded norms,"
  arXiv:1207.5726 (2012), Sec. 3 — the Choi-based primal/dual diamond-norm SDP.
- **Watrous**, *The Theory of Quantum Information*, Sec. 3.3 (diamond / cb-trace norm).
- **QETLAB** `DiamondNorm.m` — the dual (minimize) form used in practice.

## The architecture: solve in double (MOSEK), certify in arb
Two INDEPENDENT feasible points, each certified in arb (weak duality brackets eta):

- Solve the **Watrous MAX primal** (= the `gen_fixtures_d24.jl` / `sdp.jl` SDP):
  `max Re tr(J^dag X)` s.t. `[[P,X],[X^dag,Q]]>=0, P+Q=I_{n^2}`.
  A primal-feasible `(P,X,Q)` gives a **rigorous LOWER bound**
  `eta >= (2/n) * Re tr(J^dag X)`  (the `(2/n)` is the Convention-A trace-n calibration).

- Solve the **QETLAB/Watrous MIN dual** as its OWN SDP (so `Y0,Y1` are PRIMAL
  variables — NO MOSEK dual-variable extraction, which is fragile, MOI Bridges
  bug):  `min (1/2)(||Tr_Y(Y0)||_inf + ||Tr_Y(Y1)||_inf)` s.t.
  `[[Y0,-J],[-J^dag,Y1]]>=0, Y0,Y1>=0`.
  A dual-feasible `(Y0,Y1)` gives a **rigorous UPPER bound**
  `eta <= (1/2)(||Tr_Y(Y0)||_inf + ||Tr_Y(Y1)||_inf)`  (no `2/n` here — QETLAB
  convention already absorbs it; `eta = optval_dual/2`). **`Tr_Y` traces the
  RIGHT/MINOR (input H) factor, KEEPING the LEFT (output K) factor ->
  `aic_mat_partial_trace_right`** (= Convex `partialtrace(.,2,[n,n])`), yielding an
  `n x n` Hermitian PSD. EMPIRICALLY PINNED in 3a: the asymmetric paper anchor gives
  ratio 4.0 for the LEFT/sys-1 trace (WRONG) vs 2.0 == eta_ref for sys-2 (RIGHT).
  (An earlier draft of this doc said LEFT/`partial_trace_left` — that was wrong.)

Strong duality holds (Slater: `P=Q=I/2, X=0` strictly feasible), so the two
bounds meet at the true `eta`; the gap is the MOSEK duality gap (~solver tol).

### Normalization table (load-bearing)
| Form | input `J` | optval | eta |
|---|---|---|---|
| Watrous (trace-1 Choi) | `J/n` | `eta` | `=optval` |
| gen_fixtures MAX (Convention-A, trace n) | `J` | `(n/2)eta` | `eta=(2/n)optval` |
| QETLAB MIN dual (trace n) | `J` | `2eta` | `eta=optval/2` |
| eig-free lower / upper | `J` | — | `||J||_F/n <= eta <= 2||J||_F` |

## Jansson correction (rigor when the slack is not certified PSD) — Thm 3.2
MOSEK's solution is only approximately feasible. If the dual slack
`D = [[Y0,-J],[-J^dag,Y1]]` is NOT certified PSD (its `herm_max_eig(-D)` ball
straddles 0), the bound is corrected, not abandoned. Jansson Thm 3.2 (mapped to
our LOWER bound): with `d = ` certified lower bound on `lambda_min(D)` (negative
when infeasible) and the a-priori primal bound `R = 2*||J||_inf * n^2` (Watrous's
trace bound on the feasible set),
```
eta_lo = (2/n) * Re tr(J^dag X)  -  (2/n) * n^2 * |d|^- * R      ( |d|^- = max(0,-d) )
```
Symmetric mu-shift for the UPPER bound: if `Y0 -> Y0 + mu*I_{n^2}` restores PSD,
then `Tr_Y(Y0+mu I_{n^2}) = Tr_Y(Y0) + mu*n*I_n`, so
```
eta_hi = (1/2)(||Tr_Y(Y0)||_inf + ||Tr_Y(Y1)||_inf) + mu*n.
```
(The earlier claim that the dual shift leaves the objective unchanged is WRONG —
the partial trace picks up `mu*n`.)

## arb certifier — exact op sequence (reuse existing primitives, NO full eig)
Inputs (double arrays from the two Julia solves) loaded to `acb_mat` at `prec`;
`J` from `aic_ucp_choi_diff` (arb). Every PSD test is a ONE-SIDED
`aic_mat_herm_max_eig(-M)<=tol` (degeneracy-robust global enclosure) — no full
degenerate eigendecomposition, so the `aic-w4o.1` wall is dodged.

LOWER: `Jdag=adjoint(J)`; `M=Jdag*X` (`acb_mat_mul`); `re=Re(acb_mat_trace(M))`;
`lo_raw=(2/n)*re`; certify `[[P,X],[X^dag,Q]]>=0` and `||P+Q-I||_F<tol`; apply the
Jansson correction if the block isn't certified PSD.

UPPER: `P0=partial_trace_right(Y0)`, `P1=partial_trace_right(Y1)` (n x n Herm PSD;
traces the MINOR/input factor = Convex sys 2, pinned by the paper anchor in 3a);
`s0=herm_max_eig(P0)`, `s1=herm_max_eig(P1)` (op-norm = max eig for PSD);
`hi_raw=(s0+s1)/2`; certify `[[Y0,-J],[-J^dag,Y1]]>=0`, `Y0,Y1>=0`; apply the
`+mu*n` correction if a shift was needed.

Return `[arb_lower(lo), arb_upper(hi)]`.

## Dispatch / precision-wall stop condition (Rule 4)
- **Tight ball** when both feasibilities certify (or correct cleanly) AND
  `(hi-lo)/lo < 0.5`: return `[lo,hi]`.
- **Fall back to eig-free** `[||J||_F/n, 2||J||_F]` when `||J||_F < ~1e-9*n`
  (the `eta_is_zero` regime) OR a feasibility ball straddles 0 at the current
  prec and a `mu`-shift correction would exceed ~10% of `eta`.
- **Fail loud** when even at `prec>=256` the tight certifier inverts (`lo>hi`)
  and the eig-free `||J||_F` ball itself straddles 0 — escalate (eta value, n,
  the `||J||_F` ball at prec=256, repro command).
- **eta->0 floor caveat:** through the double-Kraus ccall path, `||J||_F` floors
  at `~n^2 * DBL_EPSILON` (~1e-13..1e-16) regardless of arb prec — the floor is
  set by the DOUBLE input, not arb. For certified tiny eta, feed the arb Choi
  (`acb_mat`) directly instead of the double round-trip.

## Validation plan
- 17-fixture golden master: `arb_lower(lo) <= eta_ref <= arb_upper(hi)` (slack
  1e-7 for MOSEK tol); lower/upper each tight to ~MOSEK tol; `hi/lo < 2` for
  non-zero fixtures.
- wolframscript oracle (`tools/diamond_oracle.wls`, machine-prec, agrees to
  ~1e-7): the arb ball must contain the oracle value.
- eta=0 oracle (7 fixtures): dispatch to eig-free, `hi < 1e-9`.
- Universality canary (aic-dbo.3): `hi/lo` must NOT grow with n across n=2,3,4.
- Adversarial: `phi_t_0p001` (tiny eta, prec 53 vs 256), `complex_qubit`
  (conjugation guard for `J^dag` and `Tr_Y`), `block_lowrank3` (degenerate
  `Y0,Y1` -> the PSD ball must resolve 0 vs negative), `phi_t4_0p2` (largest n).

## Risk register (carry these; escalate where noted)
1. **MOSEK dual extraction fragility** -> mitigated by solving the dual form as a
   second SDP (Y0,Y1 = its primal vars). Escalate if its PSD cert fails consistently.
2. **Degenerate-slack PSD cert at optimality** (`lambda_min(D)~0`): increase prec,
   else `mu`-shift with the `+mu*n` correction, else eig-free fallback. Escalate
   if `mu > eta/10`.
3. **eta->0 double-Kraus floor** (~1e-16): dispatch to eig-free below the floor;
   feed arb Choi directly for certified tiny eta.
4. **MOSEK tolerance** too loose for the slack PSD cert: tighten
   `MSK_DPAR_INTPNT_CO_TOL_DFEAS` to 1e-12; the Jansson correction term keeps the
   bound rigorous regardless. Escalate if the correction exceeds 50% of eta.
5. **Partial-trace direction** — `Tr_Y` over the RIGHT/MINOR (input H) factor
   (`aic_mat_partial_trace_right`), EMPIRICALLY PINNED in 3a (asymmetric paper
   anchor: sys-1 ratio 4.0 WRONG vs sys-2 ratio 2.0 == eta_ref). Mutation-test
   left-vs-right on an ASYMMETRIC channel (`complex_qubit`/`paper`) — a wrong
   direction must turn the test RED (a symmetric channel cannot catch it).

## Build order (connection-aware; "no parallel Julia")
- **3a (Julia):** add the QETLAB MIN dual solve to `AlmostIdempotentChannels.jl`
  (returns `(Y0,Y1)` + exposes the primal `(X,P,Q)`); tighten MOSEK tol; assert
  primal optval == dual optval (strong duality) and both match `eta_ref`.
- **3b (C shim + arb certifier):** `aic_cbnorm_certify(...)` taking the two
  feasible points as double arrays -> the arb op sequence above -> `[lo,hi]`,
  with the correction + dispatch + fail-loud; extend the shim ABI.
- **3c (validation):** the plan above, mutation-proven, vs golden master +
  wolframscript oracle + adversarial corpus + universality canary.
