# ALGORITHM.md — canonical math+code narrative

Per CLAUDE.md Rule 11, this is the math+code narrative for each implemented
module: which `paper/src/approximate_algebras.tex` result it realizes, the
paper's proof technique vs the constructive route taken, the defensive checks,
and the precision argument. Every formula is cited to a `.tex` line. Sections
are appended as modules land.

---

## Module `ucp` — UCP-map representations and the carrier (bead aic-c7n)

Realizes §11.1–11.2: the Choi/Stinespring representations (`.tex:1566-1633`,
`prop_KLHG` `.tex:1635`), the cb-norm bound (`.tex:1717-1719`), and the carrier
results `lem_carrier` (`.tex:1724`), `cor_carrier` (`.tex:1731`), `PhiX_M`
(`.tex:1740`).

Files: `include/aic_ucp.h`, `src/aic_ucp_core.c` (arb path, 193 LOC),
`src/aic_ucp_carrier.c` (arb path, 123 LOC), `src/aic_ucp_latd.c` (double path,
122 LOC), tests `tests/test_ucp.c` (368 checks).

### Data model — Kraus operators, Heisenberg picture (load-bearing)

A UCP map `Phi : B(K) -> B(H)` is stored as `r` Kraus operators `K_a : H -> K`
(each `dim_K x dim_H`). The action is the OBSERVABLE (Heisenberg) form

```
Phi(X) = sum_a K_a^dag X K_a          (X in B(K)),
```

which is the `.tex:1570` isometry form `Phi(X) = V^dag (X (x) 1_F) V` (with
`V^dag V = 1_H`, `.tex:1571`) written in an F-basis with V the column-stack
`V[a*dim_K + i, j] = K_a[i,j]`, `V : H -> K (x) F`, `F = C^r`. Unitality is
`Phi(1_K) = sum_a K_a^dag K_a = 1_H` (`.tex:1571`).

This is the OBSERVABLE map, NOT the Schrodinger/state dual
`T = Phi^* (X) = sum_a K_a X K_a^dag`. The unital test
`sum_a K_a^dag K_a = 1_H` pins the direction (CLAUDE.md "UCP vs CPTP" callout)
and the tests assert it both from Kraus and from the Choi partial trace, then
cross-check the two agree — a wrong partial-trace DIRECTION makes them disagree.
That cross-check is transpose-invariant, though, so it does NOT see the per-entry
conjugation convention; the conjugation/index convention is validated separately
by the Choi-convention oracle test (Choi block `(i,j)` == `aic_ucp_apply(E_ij)`,
with a complex asymmetric channel; mutation-proven against the conj-factor flip).

### Choi matrix — Convention A (paper / Watrous / QETLAB / Qiskit)

With ONB `{e_i}` of K and `E_{ij} = |e_i><e_j|`,

```
C_Phi = sum_{i,j} E_{ij} (x) Phi(E_{ij})  in B(K) (x) B(H),   [shard G; .tex:1575]
```

K factor LEFT/MAJOR, H factor RIGHT (matching `aic_mat.h`'s left-major,
row-major Kronecker). So `C_Phi` is `(dim_K*dim_H) x (dim_K*dim_H)`. Hand-deriving
from `Phi(X) = sum_x K_x^dag X K_x`: since
`(K_x^dag E_{ij} K_x)[a,b] = sum_{p,q} conj(K_x[p,a]) E_{ij}[p,q] K_x[q,b]` and
`E_{ij}[p,q] = delta_{pi} delta_{qj}`,

```
Phi(E_{ij})[a,b] = sum_x conj(K_x[i,a]) * K_x[j,b]    (conj on the FIRST factor),
```

so

```
C_Phi[i*dim_H + a, j*dim_H + b] = sum_x conj(K_x[i,a]) * K_x[j,b].
```

- CP iff `C_Phi >= 0` (PSD).
- Unital iff `Tr_K(C_Phi) = 1_H` — tracing out the LEFT (K) factor leaves the
  `dim_H` identity, via `aic_mat_partial_trace_left` (traces the LEFT, keeps the
  RIGHT). `aic_ucp_kraus_to_choi` builds it in one matrix loop (no eigenproblem;
  exact up to rounding).

### Carrier — `lem_carrier` (`.tex:1724`)

Paper technique: the carrier M is the support of `Phi^*(rho_0)` for a full-rank
state `rho_0`; `lem_carrier` characterizes it as the smallest `M` with
`M (x) F >= Im V`. The proof is an existence/equivalence statement.

Constructive finite-dim route: M is the K-marginal of `Im V`, i.e. the support
(range) of `Q = Tr_F(V V^dag)`. In Kraus coords
`(V V^dag)[(i,a),(j,b)] = (K_a K_b^dag)[i,j]`, so the F-partial-trace (`a=b`)
gives

```
Q = sum_a K_a K_a^dag          (dim_K x dim_K, Hermitian PSD).     [aic_ucp_carrier_Q]
```

The carrier is `range(Q)`. SKEPTICISM NOTE (CLAUDE.md Rule 2): the shard's first
formula `Q = sum_a K_a^dag K_a` equals `1_H` by unitality — the WRONG marginal;
the shard self-corrects, and `Q = sum_a K_a K_a^dag` is the right one. This is
mutation-proven: substituting `K^dag K` makes the eta=0 oracle's carrier-rank
check fail (rank 4 instead of the carrier dim 2).

The RANGE/RANK of Q (the carrier dimension) needs a degenerate-spectrum
eigensolver (Q has eigenvalues 0 with multiplicity `dim_K - dim M`). The arb
certified simple-eig path aborts on repeated eigenvalues (bead aic-w4o.1), so
the rank is computed on the LAPACK double path
(`aic_ucp_carrier_rank_latd`, `LAPACKE_zheev`, threshold
`dim_K * DBL_EPSILON * ||Q||_F`). The CERTIFIED carrier rank is gap-dependent
and is a documented module TODO (blocked on aic-w4o.1).

### `cor_carrier` (`.tex:1731`) and `PhiX_M` (`.tex:1740`)

`cor_carrier`: X annihilates the carrier iff `(X (x) 1_F) V = 0`. In Kraus
coords `(X (x) 1_F) V` is the column-stack of `(X K_a)`, so the defect is
`||stack_a(X K_a)||_op` (0 iff `Ker X >= M`). Tested both ways: `X = 1 - Pi_M`
gives defect ~0, `X = Pi_M` gives defect bounded away from 0.

`PhiX_M`: `Phi(X) = Phi(Pi_M X Pi_M)` for ANY UCP map (proof `.tex:1742` inserts
`Pi_M (x) 1_F` on both sides of V using `(Pi_M (x) 1_F)V = V`). The routine
returns `||Phi(X) - Phi(Pi_M X Pi_M)||_op` (0 up to rounding).

### CP certificate and cb-norm

CP (`C_Phi >= 0`): `aic_ucp_is_cp_choi` encloses the largest eigenvalue of
`-C_Phi` with `aic_mat_herm_max_eig` (degeneracy-robust global enclosure, NOT
the simple-eig path the degenerate Choi spectrum would abort). Verdict:
`+1` if the ball is certainly `<= tol` (PSD certified), `0` if certainly
`> tol` (not PSD), ABORT if it straddles tol (gap unresolved at this prec —
Rule 4, fail loud, never guess). Mutation context: a Choi with
`lambda_min = -1e-6` (far above prec=53 noise) is certified NOT PSD.

cb-norm: closed form `||Phi||_cb = ||Phi(1_K)||_op` (Paulsen Prop 3.6;
`.tex:1717-1719` gives `<= 1`, with `= 1` for UCP per `.tex:359`). This is the
ONLY cb-norm in scope. The idempotency defect `||Phi^2-Phi||_cb` needs the
Watrous diamond-norm SDP (bead aic-d24, out of scope, blocked infra); no
plausible-but-wrong Choi-opnorm surrogate is shipped.

### Choi -> Kraus extraction (double path)

`C_Phi = sum_a lambda_a v_a v_a^dag` (PSD eigendecomp), i.e.
`C[i*dh+a, j*dh+b] = sum_a lambda_a v_a[i*dh+a] conj(v_a[j*dh+b])`. Matching this
against Convention A `C[i*dh+a, j*dh+b] = sum_x conj(K_x[i,a]) K_x[j,b]`
term-by-term forces `conj(K_a[i,c]) = sqrt(lambda_a) v_a[i*dh+c]`, so

```
K_a[i,c] = sqrt(lambda_a) * conj(v_a[i*dim_H + c])    (CONJUGATE reshape).
```

The CONJUGATE is load-bearing and must MATCH the `kraus_to_choi` conjugation
(H2 lockstep): flip either one and the complex round-trip `C != C'`. (Real
channels have real eigenvectors, so the conjugation is a no-op there — only the
complex asymmetric channel exercises it.) The reshaped eigenvectors (row index
`i*dim_H + c -> (i,c)`) ARE the Kraus ops up to the `sqrt(eigenvalue)` scale.
Keep `lambda_a > (dim_K*dim_H)*DBL_EPSILON*||C||_F` (QETLAB-style); drop the
numerically-zero tail. The recovered set is a valid Kraus rep of the SAME channel
only up to the unitary gauge freedom of Kraus decompositions — compared AS
CHANNELS (rebuild Choi), never operator-by-operator. Uses `LAPACKE_zheev`
(degenerate Choi spectra; bead aic-w4o.1). Mutations proven RED: missing `sqrt`,
transposed reshape, flipped reshape conjugation (complex round-trip breaks).

### Precision argument

The unital/CP/carrier quantities are differences of O(1) operators, so `prec=53`
resolves them to ~1e-15; the tests cross-check double-vs-arb@53 (Choi entries,
carrier-Q top eigenvalue, cbnorm) within ~1e-12, the round-trip Choi within
~1e-10 (double-path extraction in the middle). Load-bearing precision (tiny
`eta`, near-singular operators) arises in the §12 almost-idempotent modules, not
here.

### Cross-check ladder coverage (CLAUDE.md Rule 6)

1. Internal sanity: unitality two ways agree; V isometry; CP certified PSD.
2. Round-trip: Kraus->Choi->Kraus->Choi reproduces C_Phi (rank-deficient and
   maximally-degenerate inputs).
3. Double-vs-arb@53: Choi entries, carrier-Q top eigenvalue, cbnorm.
4. eta=0 exact-idempotent oracle (compression-idempotent Phi): Phi^2=Phi,
   PhiX_M exact, carrier spectrum {0,1} (rank = carrier dim), CP+unital
   certified, cbnorm=1.
5. Closed forms: identity channel -> Choi = unnormalized maximally-entangled
   projector (`C[i*d+a,j*d+b] = [i==a][j==b]`), cbnorm=1; state channel
   `Phi(X)=Tr(rho X)` (dim_H=1) -> carrier = supp(rho).
6. Adversarial: rank-deficient Choi; depolarizing (all-equal Choi eigenvalues);
   proper-subspace carrier (dim_K=4, M=2); non-self-map K!=H (B(C^2)->B(C^3));
   near-non-CP rejected loudly.

### Open items / TODO (for the orchestrator to bead)

- Certified arb carrier RANK / range: blocked on aic-w4o.1 (degenerate certified
  eig). Only the uncertified double-path rank is shipped.
- Certified arb Choi->Kraus extraction: same blocker (Choi spectrum degenerate).
  A Cholesky-based eig-free extraction (arb-friendly) was NOT implemented this
  session — deferred as an audition follow-up (Law 4).
- `||Phi^2-Phi||_cb` diamond-norm SDP: bead aic-d24 — FIRST INCREMENT shipped
  (the Choi substrate + Julia golden master; see the `cb-norm` section below).

---

## Module `cb-norm` / eta-idempotence defect (bead aic-d24)

Computes the central quantity of the paper, the eta-idempotence defect

> `.tex:354`, verbatim: "A UCP map Phi: B(H)->B(H) is called eta-idempotent if
> ||Phi^2-Phi||_cb <= eta."

`Lambda = Phi^2 - Phi` is Hermiticity-preserving but NOT completely positive, so
the CP closed form `||Phi||_cb = ||Phi(1)||_op` (`.tex:1717-1719`, used by
`aic_ucp_cbnorm_cp`) does NOT apply. The cb-norm equals the diamond norm and is
"efficiently computable using semidefinite programming" (`.tex:352`); for a
Hermiticity-preserving map the ampliation truncation `N = dim(input)` is rigorous
(Watrous 2012), so a single SDP at the channel's own dimension suffices.

### Architecture (Route B) — C computes the Choi, Julia+MOSEK solves the SDP

No arb-native SDP solver exists, so the work splits:

- **C core** (`src/aic_ucp_compose.c`, 104 LOC; decls in `include/aic_ucp.h`)
  supplies the SDP its input — the Choi matrix of `Lambda` — and nothing else.
  Two routines, both arb/acb (explicit `prec`):
  - `aic_ucp_compose(out, phi, psi)` — Kraus of `Phi o Psi`. Heisenberg
    composition (derived in the header, verbatim there):
    `(Phi o Psi)(X) = sum_a K_a^dag (sum_b L_b^dag X L_b) K_a
                    = sum_{a,b} (L_b K_a)^dag X (L_b K_a)`,
    so `Kraus(Phi o Psi) = { L_b K_a } = { psi.K[b] @ phi.K[a] }`. Requires
    `phi->dim_K == psi->dim_H` (asserted, fail loud); result shape
    `(psi->dim_K, phi->dim_H)`, `r = phi->r * psi->r`. `Phi^2 = compose(Phi,Phi)`.
  - `aic_ucp_choi_diff(C, phi1, phi2)` — `C = Choi(phi1) - Choi(phi2)`,
    Convention A. Choi is linear in the map, so `Choi(Phi^2 - Phi) =
    Choi(Phi^2) - Choi(Phi)`; pass `phi1 = compose(Phi,Phi)`, `phi2 = Phi`. The
    result is Hermitian but generically indefinite (why the SDP is needed).
- **Julia golden master** (`tools/gen_fixtures_d24.jl`, env `julia/env/`):
  Convex.jl 0.16.6 + MosekTools 0.15.10 + Mosek 11.2.0 (NO PYTHON). The Watrous
  2012 SDP (complex-native): vars `X` (n^2 x n^2 complex), `P,Q` HermPSD;
  `maximize Re tr(J' X)` s.t. `[[P,X],[X',Q]] >= 0`, `P + Q = I_{n^2}`.
  NORMALIZATION (load-bearing): Convention-A Choi has trace `n` for a UCP map but
  the SDP is trace-1 calibrated, so `||Lambda||_diamond = (2/n) * SDP_optval`.

### Verification — analytic anchors and the C-vs-Julia Choi cross-check

The generator ASSERTS the analytic anchors before emitting any fixture (fails
loud otherwise), so a broken solver/normalization cannot silently poison the
golden master. Anchors (all reproduced to <1e-6 at generation):

- `id_d`: `Lambda = 0` => `eta = 0`.
- `Dep_d - id_d`: `2(1 - 1/d^2)` (1.5 at d=2; 16/9 at d=3; 15/8=1.875 at d=4 —
  the dimension canary confirming the `2/n` normalization, now asserted at ALL
  three dimensions that appear in the corpus, not just d=2).
- `Phi_t = (1-t) id + t Dep`: `Phi_t^2 - Phi_t = t(1-t)(Dep - id)`, so
  `eta = t(1-t) * 2(1-1/d^2)` (0.315 / 0.135 / 0.01485 / 0.0014985 at d=2,
  t=0.3/0.1/0.01/0.001; 0.084444 at d=3,t=0.05; 0.30 at d=4,t=0.2). The
  poison-pill `assert_anchors()` checks this general formula at d=2,3,4 so a
  broken solver/normalization cannot silently poison any dimension-varying
  fixture.
- **The paper's own example** (`.tex:367-378`): the C^2 channel
  `Phi(X)=P0 Tr(g0 X)+P1 Tr(g1 X)` has `||Phi^2-Phi||_cb = eta*sqrt(1-eta)`
  (verbatim `.tex:378`); reproduced as 0.09486833 at eta=0.1, 0.009949873 at 0.01.

`tests/fixtures_d24.inc.h` (generated, committed) carries, per fixture: `n, r`,
the Kraus ops (aic_ucp layout), the diamond-norm `eta_ref`, and the reference
`Choi(Phi^2-Phi)`. `tests/test_ucp_d24.c` rebuilds each channel, computes `Phi^2`
via `aic_ucp_compose` and `Choi(Lambda)` via `aic_ucp_choi_diff`, and asserts it
matches the Julia reference — a genuine **C-vs-Julia cross-check of the Choi
substrate** (max-diff 2.7e-15 over the 17-fixture corpus; worst at the n=4
`phi_t4_0p2`), plus the eta=0 oracles carry `eta_ref < 1e-9`. Regenerate with
`make fixtures` (serial Julia; not a `make test` prerequisite — the .inc.h is
committed). Two of the 17 fixtures are dedicated hardening cases:
- `complex_qubit` (`K_0=[[0.6,0],[0.8i,0]]`, `K_1=[[0,0.8i],[0,0.6]]`, unital,
  eta≈1.28): the FIRST fixture with genuinely complex Kraus ops
  (max|Im(Kraus)|=0.8) and a Lambda-Choi with nonzero imaginary off-diagonals
  (max|Im|≈0.61), so the C-vs-Julia cross-check now GUARDS the conjugation in
  `aic_ucp_kraus_to_choi`. Teeth-proven: dropping the `conj` makes the C Choi
  diverge from the Julia reference by 1.28 in operator norm (≫ the 1e-11 fixture
  tolerance), turning the assertion RED; restored to 5.9e-16.
- `block_lowrank3` (`Phi_t` on the upper 2-block of C^3 ⊕ identity on e2, unital,
  eta=0.16): a genuinely RANK-DEFICIENT Lambda-Choi (9×9, rank 4, with real zero
  eigenvalues), replacing the earlier mislabeled `partial_dephase3` whose Lambda
  was full rank.

### Cross-check ladder coverage (Rule 6) and mutation proofs

- compose correctness: `apply(compose(Phi,Psi),X) == Phi(Psi(X))` (the definition,
  an independent path), incl. non-self-map chaining `B(C^3)->B(C^2)->B(C^3)`.
- compose preserves unitality; `Phi^2` Kraus count `r^2`, shape correct.
- eta=0 Choi oracle: idempotent `Phi` => `Choi(Phi^2-Phi) == 0` to machine zero,
  on identity / Dep_2 / Dep_3 / dephasing_M2 / block_cond_exp.
- double-vs-arb@53 on the composed Kraus ops and on the Choi-diff entries.
- MUTATION PROOF (Rule 7): the wrong product order `K_a L_b` (instead of
  `L_b K_a`) breaks the compose cross-check on the non-commuting
  complex-channel/dephasing pair — verified by a source mutation that turned the
  test RED, then restored.

### Deferred (beaded / escalation)

- The **certified arb cb-ball**: no arb-native SDP solver exists; the credible
  route is solve-in-double then a rational-certificate / verified-SDP
  post-processing (KKT + duality-gap in FLINT/arb). NOT implemented; escalation
  note — the certified ladder rung (4) for `eta` is unreachable without it.
- The **C-callable `aic_ucp_idemp_defect_cb`** that invokes the SDP via `ccall`:
  belongs in the Julia package (E5, bead aic-obc). For now the diamond-norm VALUE
  is produced only by the Julia generator / golden master.

### cbnorm (eig-free ball) — the always-valid fallback rung (bead aic-m24, incr. 1)

An EIG-FREE, certified, two-sided bracket on the eta-idempotence defect
`eta = ||Phi^2-Phi||_cb = ||Phi^2-Phi||_diamond` (`.tex:347-354`), using ONLY the
certified Frobenius norm of `J = Choi(Lambda)`, `Lambda = Phi^2 - Phi`
(Convention-A, `n^2 x n^2`, Hermitian). No SDP solver, no eigendecomposition (so
it dodges the degenerate-eig wall, bead aic-w4o.1) and it NEVER aborts (the only
assert is a fail-loud shape check, Rule 4). This is the cheapest rung of the
certified-bound ladder and the fallback the future tight certifier dispatches to
near `eta=0`. Files: `include/aic_cbnorm.h`, `src/aic_cbnorm.c` (108 LOC, arb
path only), `tests/test_cbnorm.c` (125 checks).

**The math (derivation).** The Watrous 2012 SDP (cb-norm section above;
`tools/gen_fixtures_d24.jl`) gives `eta = (2/n)*optval`, where
`optval = max Re tr(J^dag X)` s.t. `[[P,X],[X^dag,Q]] >= 0`, `P+Q = I_{n^2}`,
`P,Q >= 0`.

- *Upper bound (eig-free).* `P+Q=I_{n^2}` with `P,Q>=0` forces `0<=P,Q<=I`
  (`||P||_op,||Q||_op<=1`) and `tr(Q) <= tr(P+Q) = n^2`. The block-PSD condition
  Schur-factors `X = P^{1/2} K Q^{1/2}` with `||K||_op<=1`, so
  `||X||_F <= ||P^{1/2}||_op ||K||_op ||Q^{1/2}||_F <= 1*1*sqrt(tr Q) <= n`.
  Then `optval <= ||J||_F ||X||_F <= n ||J||_F`, giving `eta <= 2 ||J||_F`.
- *Lower bound (eig-free).* `J = (id_n (x) Lambda)(omega)`, `omega=|Omega><Omega|`,
  `|Omega>=sum_i |i>|i>` (unnormalized). For the normalized maximally-entangled
  `rho = omega/n`, `(id (x) Lambda)(rho) = J/n`, so (diamond-norm sup over input
  states; truncation `N=n` rigorous for a Herm-preserving map, Watrous 2012, bead
  aic-d24) `eta >= ||(id (x) Lambda)(rho)||_1 = ||J||_1/n >= ||J||_2/n = ||J||_F/n`
  (Schatten `||.||_1 >= ||.||_2 = ||.||_F`).

So the certified bracket is `[ ||J||_F/n , 2*||J||_F ]`. Ratio `hi/lo = 2n` —
valid but loose by design; the tight ball is a later increment and is NOT
attempted here.

**API.** `aic_cbnorm_eigfree_ball_choi(lo,hi,J,n,prec)` (from the `n^2 x n^2`
Choi) and `aic_cbnorm_eigfree_ball(lo,hi,phi,prec)` (from a UCP self-map: builds
`Phi^2` via `aic_ucp_compose`, `J` via `aic_ucp_choi_diff`, then defers).
Implementation: `fro = ||J||_F` via `aic_mat_frobenius_norm` (certified ball);
`lo = fro/n` (`arb_div_ui`), `hi = 2*fro` (`arb_mul_2exp_si` by +1, exact).
Returning the arb BALLS suffices for rigor: the true `||J||_F` lies in the `fro`
ball, so by monotonicity `arb_lower(lo) <= eta` and `arb_upper(hi) >= eta`.

**Verification — two assertion families over all 17 d24 fixtures.**
- (A) CORRECTNESS / validity: build `J` from each fixture's golden Choi, assert
  the bracket CONTAINS `eta_ref` (`lo <= eta_ref <= hi`, slack `1e-7` for the
  SDP solver tolerance). A wrong `n`-power or a sign error makes some fixture's
  bracket miss `eta_ref`. Tightest upper margin `hi - eta_ref = 1.96e-3`, on the
  smallest nonzero `eta` (`phi_t_0p001`, `eta_ref=1.499e-3`, `2*fF=3.46e-3`).
- (B) TEETH / exactness: INDEPENDENTLY recompute `fF = ||J||_F` by a direct
  double sum over the stored Choi arrays, then assert `lo == fF/n` and
  `hi == 2*fF` to `~1e-12`. This pins BOTH constants. Measured worst-case
  midpoint errors: `|lo - fF/n| = 1.11e-16`, `|hi - 2fF| = 8.88e-16` over the
  corpus.
- (C) eta=0 oracle: for the seven `eta_is_zero` fixtures `J~0`, so the bracket
  collapses to `~[0,0]`: `arb_upper(hi) < 1e-9` (all `hi = 0` exactly).
- (D) precision self-consistency: `prec=53` vs `prec=256` on `phi_t_0p3` and
  `phi_t4_0p2` both bracket `fF/n` and `2*fF`.
- (E) wrapper path: `aic_cbnorm_eigfree_ball` (rebuilds `Phi`, runs
  compose+choi_diff) agrees with the choi-level path within `1e-9` on
  `phi_t_0p3 / complex_qubit / phi_t4_0p2 / block_lowrank3` — cross-checks the
  whole `Phi -> Phi^2 -> J -> ball` pipeline.

**Mutation proofs (Rule 7).** (1) `2*fro -> 1*fro` (drop the factor 2): caught
by family (B) at `phi_t_0p3` (`hi == 2*fF = 0.727461` failed; the mutated
`hi = fF = 0.363731`). Family (A) does NOT catch this mutation — the bracket is
loose enough that `fF >= eta_ref` still holds on every fixture, so the
constant-pinning family (B) is exactly what gives this teeth. (2) `fro/n -> fro`
(drop the `1/n`): caught by family (A) at `phi_t_0p3`: the mutated
`lo = fF = 0.363731` exceeds `eta_ref = 0.315`, violating the rigorous lower-bound
containment `lo <= eta`. Both restored. No teeth gaps: each constant is pinned by
at least one family on this corpus.

**Precision argument.** `||J||_F` is a sum of squares of `O(1)` Choi entries, so
`prec=53` resolves `fF` to `~1e-15` (family-D 53-vs-256 agreement confirms); no
near-singular inversion or spectral gap arises (this rung is deliberately
eig/SDP-free). The arb balls are the rigorous object; the bound holds by ball
monotonicity at any `prec`.

**Deferred (now BUILT).** The TIGHT certified ball (ratio `-> 1`) is the next
subsection (`cbnorm tight certifier`, aic-m24 incr. 3). It takes the
solve-in-double (MOSEK) feasible points and restores them to exact feasibility in
arb (the "verified-SDP" route). This eig-free bracket is the always-valid
fallback the tight certifier dispatches to near `eta=0` and when the MOSEK points
cannot be restored.

### cbnorm tight certifier — MOSEK feasible points + arb restoration (bead aic-m24, incr. 3b)

The TIGHT certified two-sided ball `[lo,hi]` on `eta = ||Phi^2-Phi||_cb`
(`.tex:347-354`), `arb_lower(lo) <= eta <= arb_upper(hi)`, gap `~` MOSEK duality
gap `+` arb radius — orders of magnitude tighter than the eig-free `2n`-wide
bracket. Files: `include/aic_cbnorm.h` (`aic_cbnorm_certify` decl),
`src/aic_cbnorm_certify.c` (114 LOC, the dispatcher),
`src/aic_cbnorm_certify_restore.c` (165 LOC, the two restoration recipes),
`src/aic_cbnorm_internal.h` (the split's internal header, also the surface the
GAP-2 tests use to read the restoration defects),
`tests/test_certify.c` (34 checks), and the flat-double `aic_cbnorm_certify_d`
shim in `src/aic_ucp_shim.c` (exported in `libaic.so` for the 3c Julia ccall).
Design doc: `docs/cbnorm_tight_certifier.md`. NO Julia run here — the C/arb
certifier consumes the committed golden-master feasible points (step 3a).

**Strategy.** Two INDEPENDENT MOSEK feasible points bracket `eta` by weak
duality: the Watrous MAX-primal `(X,P,Q)` -> LOWER, the QETLAB/Watrous MIN-dual
`(Y0,Y1)` -> UPPER. Each MOSEK point is only APPROXIMATELY feasible, so each is
RESTORED to EXACT feasibility in arb before its objective is read off; only then
is the bound rigorous. Every PSD test is the ONE-SIDED `aic_mat_herm_max_eig(-M)`
global enclosure (degeneracy-robust; NO full eig, dodging the aic-w4o.1 wall).

**LOWER recipe (convex-combination restoration toward the Slater center
`(0, I/2, I/2)`).** `max Re tr(J^dag X)` s.t. `[[P,X],[X^dag,Q]] >= 0`,
`P+Q = I_{n^2}`; `eta = (2/n)*optval` (Convention-A trace-`n` calibration).
1. symmetrize the equality: `E = P+Q-I`; `P' = P-E/2`, `Q' = Q-E/2` (so `P'+Q'=I`
   exactly in arb).
2. `delta =` rigorous UPPER bound on `max(0, -lambda_min(block'))`,
   `block' = [[P',X],[X^dag,Q']]`, via `herm_max_eig(-block')` (arb_upper, clamp 0).
3. `t = 1/(1 + 2*delta)`; the convex point `(tX, tP'+(1-t)I/2, tQ'+(1-t)I/2)` is
   FEASIBLE (`lambda_min(t*block' + (1-t)I/2) >= -t*delta + (1-t)/2 >= 0`;
   equality holds; block PSD => diagonal blocks PSD).
4. `lo = t*(2/n)*Re tr(J^dag X)` (`<= eta`).

**UPPER recipe (eigenvalue-shift restoration).** `min (1/2)(||Tr_2 Y0||_inf +
||Tr_2 Y1||_inf)` s.t. `[[Y0,-J],[-J^dag,Y1]] >= 0`, `Y0,Y1 >= 0`;
`eta = optval/2` (QETLAB convention).
1. `eps =` rigorous UPPER bound on `max(0, -lambda_min(block_D))`,
   `block_D = [[Y0,-J],[-J^dag,Y1]]`.
2. `(Y0+eps I, Y1+eps I)` is dual-feasible. `Tr_2 = aic_mat_partial_trace_right`
   (RIGHT/MINOR input factor, PINNED in 3a). The marginal is Hermitian PSD, so
   `||.||_inf = lambda_max = herm_max_eig`.
3. `Tr_2(Y0+eps I_{n^2}) = Tr_2(Y0) + eps*n*I_n`, so
   `hi = (1/2)(lambda_max(Tr_2 Y0) + lambda_max(Tr_2 Y1)) + eps*n` (`>= eta`).

**Dispatch / fail-loud (Rule 4).** `||J||_F < ~1e-9*n` (eta=0 regime) ->
`aic_cbnorm_eigfree_ball_choi` (~`[0,0]`). `delta` or `eps > 0.1` (MOSEK badly
off / precision wall) OR `arb_lower(lo) > arb_upper(hi)` (inverted) -> eig-free
bracket as the rigorous fallback; if even THAT straddles, ABORT with
`eta/n/prec` (the restoration recipe is wrong; orchestrator re-derives).

**Measured results (per fixture, prec=128).** Containment holds for all 5
fixtures; widths are tiny:
- `phi_t_0p3`  (n=2): `[0.31499999999, 0.31500000000]`, width `5.8e-13`, eta_ref `0.315`.
- `complex_qubit` (n=2): `[1.28000000000, 1.28000000000]`, width `6.9e-14`, eta_ref `1.28`.
- `phi_t4_0p2` (n=4): `[0.29999999994, 0.30000000001]`, width `6.0e-12`, eta_ref `0.300`.
- `paper_0p1` (n=2, ASYMMETRIC): `[0.094868329805, 0.094868329805]`, width `4.5e-14`,
  eta_ref `0.1*sqrt(0.9) = 0.0948683`. THE direction-teeth fixture (see below).
- `id_2` (eta=0): `[0,0]` exactly (eig-free dispatch), `hi < 1e-9`.

On the natural MOSEK feasible points the solution is so accurate that `delta = 0`
and `eps <= 7e-16` for the whole corpus: the restoration corrections (`t ~ 1`,
`eps*n ~ 1e-16`) are NEGLIGIBLE, well below the `1e-7` containment slack. The
corrections are EXERCISED by synthetic perturbation in test family (G), not by
the natural corpus (see the GAP-2 closure below). **Universality canary:** the
relative width `(hi-lo)/eta_ref` is `1.8e-12` (n=2) vs `2.0e-11` (n=4) — it does
NOT grow with `n` (the bound's dimension-independence is preserved, `.tex:484`).
**double-vs-arb:** prec=64 and prec=256 give the SAME width (`5.98e-12`) on
`phi_t4_0p2` — the floor is the DOUBLE MOSEK feasible point, NOT the arb prec, so
higher prec does not tighten further (consistent with the design-doc
double-Kraus-floor caveat).

**Mutation-proof results (Rule 7; RESTORED after each).**
- *(3a) drop `(2/n)` in LOWER* -> RED, caught by `phi_t_0p3` (tightness: width
  jumps to `0.55` as `lo` becomes `2*Re tr ~ 2*eta`).
- *(3b) use `optval` not `optval/2` in UPPER* (drop the `/2`) -> RED, caught by
  `phi_t_0p3` (tightness: width `0.315` as `hi` becomes `2*eta`).

**GAP 1 CLOSED — partial-trace DIRECTION teeth (asymmetric fixture).** The earlier
4-channel corpus had NUMERICALLY SYMMETRIC dual marginals
(`||Tr_R(Y0) - Tr_L(Y0)||_F ~ 1e-13` for all four), so swapping
`partial_trace_right -> partial_trace_left` was SILENT — the load-bearing
direction (design-doc Risk #5) had no teeth. The `paper_0p1` fixture (the paper
example `Phi(X)=P0 Tr(g0 X)+P1 Tr(g1 X)`, `.tex:367-378`, `etap=0.1`) has a
GENUINELY ASYMMETRIC MIN-dual optimum: `||Tr_R(Y0)-Tr_L(Y0)||_F = 1.342e-01`
(vs the symmetric corpus' `1e-13`..`1e-16`), measured in test family (F) and
asserted there (`paper_0p1` asym `>= 0.05`; the other fixtures asym `< 1e-9`).
With this fixture present:
  - *swap `partial_trace_right -> partial_trace_left` in `aic_cbnorm_int_upper`*
    -> **RED**, caught by `paper_0p1` tightness (B): `Tr_L` gives `hi = 2*eta`, so
    width jumps to `9.487e-02 ~ eta_ref` (`!< 1e-4`). The symmetric fixtures stay
    tight, so `paper_0p1` is the channel that catches the direction. Restored.
  The generator (`tools/gen_fixtures_certify.jl`) re-asserts the asymmetry (`>0.05`)
  on any `paper*` fixture before emitting, so a future regen that loses it fails
  loud rather than silently disarming the teeth.

**GAP 2 CLOSED — restoration-correction teeth (synthetic perturbation).** MOSEK is
too accurate to hand us an infeasible point, so test family (G) MANUFACTURES one
in-test (calling the internal `aic_cbnorm_int_lower`/`_upper` to read `delta`/`eps`):
  - *LOWER:* scale `X -> 2X`. The primal block goes non-PSD by a known amount:
    `delta = 0.5` (so `t = 1/(1+2*0.5) = 0.5`). The raw `(2/n)Re tr(J^dag 2X) =
    0.63 = 2*eta` would EXCEED eta; the `t`-scaling pulls `lo` back to `0.315 =
    eta_ref` EXACTLY. Asserts `delta > 0`, `lo < raw` (t<1), `lo <= eta_ref`.
    *Mutation: set `t=1` (drop the t-scaling)* -> **RED** (`lo=0.63 !< raw`, then
    `lo > eta_ref`). Restored.
  - *UPPER:* shift `Y0,Y1 -> Y0-sI, Y1-sI` with `s=0.05`. `block_D` goes non-PSD:
    `eps = 0.05`. The marginal `lambda_max` drops by `s*n=0.1` to `0.215`; the
    `+eps*n` shift restores `hi = 0.315 = eta_ref` EXACTLY. Asserts `eps > 0`,
    `hi >= eta_ref`. *Mutation: drop the `+eps*n` shift* -> **RED**
    (`hi=0.215 < eta_ref=0.315`). Restored.
  Both corrections are now PROVEN: the defect goes `>0`, the correction activates,
  rigor is preserved, and dropping each correction breaks the rigor assertion.

**Precision argument.** The bounds are differences/sums of `O(1)` quantities
(`Re tr(J^dag X)`, marginal `lambda_max`) plus the tiny restoration corrections;
no near-singular inversion arises (the restoration AVOIDS the `(L+R)^{-1}` route).
`prec=53..128` resolves everything to `~1e-13`; the floor is the double MOSEK
input, not arb. The arb balls are the rigorous object; the bound holds by ball
monotonicity at any `prec`.

### ucp shim / libaic.so — flat-double ccall surface (bead aic-m24, incr. 2a)

The C substrate for the reusable `eta_idempotence(kraus)` entry point. A
FLAT-double ABI over the arb cores so the (later, increment 2b) Julia package can
`ccall` them WITHOUT any FLINT type crossing the language boundary. Files:
`include/aic_ucp_shim.h`, `src/aic_ucp_shim.c` (206 LOC, arb path internally; the
incr.-3b `aic_cbnorm_certify_d` was added here), `tests/test_shim.c` (170 checks),
plus the `lib` target in the `Makefile`.

**Why a shim (the ccall contract).** Julia `ccall` passes `int` / `double` /
`Ptr{Cdouble}` cleanly, but NOT FLINT structs — `acb_mat_t` is an opaque array
type, `slong` is FLINT's machine-word int. So the public boundary uses `int`
dimensions and `Ptr{Cdouble}` arrays only; the shim converts to/from the internal
`aic_ucp_kraus` / `acb_mat` AT the boundary, runs the SAME cores the C tests
exercise (`aic_ucp_compose` + `aic_ucp_choi_diff`; `aic_cbnorm_eigfree_ball`),
and writes the result back into caller-owned flat double arrays. Pure
marshalling — NO new math.

**ABI / layout (matches `fixtures_d24.inc.h` and Convention-A, `aic_ucp.h:23-42`).**
- `kraus_re[a*n*n + i*n + j]`, `kraus_im[...]` — `r` ops, each `n x n` (`K_a`).
- `choi_re[p*N + q]`, `choi_im[...]` — `N = n*n`, row-major.
- `int aic_ucp_choi_diff_d(choi_re, choi_im, kraus_re, kraus_im, n, r, prec)`
  builds `Phi`, `Phi^2 = compose(Phi,Phi)`, `J = choi_diff(Phi^2, Phi) =
  Choi(Phi^2 - Phi)` (Convention-A, `n^2 x n^2`), extracts `J` to the double
  arrays (`arf_get_d` on the real/imag midpoints). Returns 0; fail-loud asserts
  on bad shape (`n>=1`, `r>=1`, `prec>=2`, non-null arrays).
- `int aic_cbnorm_eigfree_d(lo, hi, kraus_re, kraus_im, n, r, prec)` wraps
  `aic_cbnorm_eigfree_ball` and converts the arb balls to a RIGOROUS double
  bracket: `*lo = FLOOR(arb_lower(lo_ball))`, `*hi = CEIL(arb_upper(hi_ball))`
  (`arb_get_lbound_arf`/`arb_get_ubound_arf` -> `arf_get_d` with
  `ARF_RND_FLOOR`/`ARF_RND_CEIL`) — the only rounding directions that keep
  `*lo <= eta <= *hi` after the double conversion. Julia gets a certified bracket
  via `ccall` with NO SDP solver.
- `int aic_cbnorm_certify_d(lo, hi, J_*, X_*, P_*, Q_*, Y0_*, Y1_*, n, prec)`
  (incr. 3b) wraps `aic_cbnorm_certify`: loads the flat `N x N` `[p*N+q]` (re,im)
  arrays for `J` and the two MOSEK feasible points into `acb_mat`, runs the TIGHT
  certifier, and converts the balls with the SAME rigorous FLOOR/CEIL rounding.
  Julia (3c) gets the MOSEK-tight certified bracket via `ccall`. Exported in
  `libaic.so` (verified: `nm -D` shows `aic_cbnorm_certify_d`).

**`libaic.so`.** `make lib` -> `build/libaic.so` via a dedicated recipe that
adds `-fPIC -shared` LOCALLY (NOT to the global `CFLAGS`, to keep the test/bench
builds untouched). `$(SRC)` globs all `src/*.c`, so the shim is auto-included.
Verified: `nm -D build/libaic.so | grep aic_ucp_choi_diff_d` shows
`T aic_ucp_choi_diff_d` (and `T aic_cbnorm_eigfree_d`).

**Round-trip cross-check (`test_shim.c`, all 17 d24 fixtures).** For each
fixture the shim's double-out Choi is compared against BOTH:
- (RT-golden) the fixture's stored golden `choi_re/choi_im` (the Julia/MOSEK
  Convention-A Choi) — worst max abs entry diff `5.551e-16` over the corpus;
- (RT-arb) the pure-arb path (build Kraus as `acb`, `compose` + `choi_diff`,
  extract to double) at `prec=53` — worst diff `0.000e+00` (the shim runs the
  identical arb cores, so the marshalling adds no error). This is the
  double-vs-arb@53 ladder rung 2 isolating the shim's own marshalling.

Plus the BONUS `aic_cbnorm_eigfree_d`: its `(lo,hi)` bracket contains `eta_ref`
(slack `1e-7` for the SDP tol) on every fixture, and its rounded doubles match
the C-level `aic_cbnorm_eigfree_ball` endpoints to `1e-12`.

**Mutation proof (Rule 7).** Negate the imaginary part on extraction
(`choi_im = -arf_get_d(...)`, i.e. DROP the Convention-A conjugation sign):
caught by the `complex_qubit` fixture (the only fixture with nonzero imaginary
Choi entries) — `shim Choi vs golden max-diff = 1.229e+00 > 1e-12`, abort. The
real-only fixtures do NOT catch it (their `choi_im` is all zero, so the sign flip
is a no-op), confirming `complex_qubit` is the load-bearing conjugation guard
(the documented Choi-conjugation bug class, `aic_ucp.h:104-107`). Restored to
GREEN (170 checks).

**Precision argument.** The shim adds no numerical error of its own beyond the
double round-trip at the boundary: it runs the certified arb cores at `prec`
(default callers pass 106 for headroom; the test uses 53 for the rung-2 anchor),
then `arf_get_d` rounds the `O(1)` Choi midpoints to double (`~1e-16`). The
cbnorm bracket stays rigorous after the double conversion by the FLOOR/CEIL
rounding above.

### AlmostIdempotentChannels.jl (Julia entry point) — bead aic-m24, incr. 2b

The Julia package `julia/AlmostIdempotentChannels.jl/` exposes the reusable
eta-idempotence VALUE entry point `eta_idempotence(kraus) -> Float64`
= `||Phi^2-Phi||_diamond` (`.tex:347-354`) that downstream modules
(`assoc_ecsa`, `factorize`) consume. It `ccall`s the flat-double shim in
`build/libaic.so` (`aic_ucp_choi_diff_d`, `aic_cbnorm_eigfree_d`) for the
certified arb `Choi(Phi^2-Phi)` and the eig-free bracket, then solves the Watrous
diamond-norm SDP in Julia + MOSEK. The certified TIGHT cb-ball is a later
increment (the SDP value is the golden-master rung here). NO PYTHON. Mirrors the
`../su2-fft/julia` ccall-C package pattern exactly.

**Layout.** `Project.toml` (deps Convex/Libdl/LinearAlgebra/Mosek/MosekTools;
PINNED `Convex = "0.16"`, `Mosek = "11"`, `MosekTools = "0.15, 0.16"` to match
`julia/env` so the SDP value equals the committed fixtures); `Manifest.toml`
(seeded from `julia/env/Manifest.toml`, resolved — Convex 0.16.6, Mosek 11.2.0,
MosekTools 0.15.10, identical to the fixture generator); `deps/build.jl`
(`make -C <repo_root> lib`, then writes `deps/deps.jl` with the abspath
`LIBAIC`; the repo root is THREE levels up — the package sits two levels deep,
one more `..` than su2-fft); `src/sdp.jl` (the Watrous SDP, verbatim from
`tools/gen_fixtures_d24.jl`, with the load-bearing `(2/n)` normalization);
`src/AlmostIdempotentChannels.jl` (the module); `test/runtests.jl`.

**dlopen / ABI.** `__init__` PRE-LOADS the system `libgmp.so.10`
(`RTLD_NOW|RTLD_GLOBAL`, full path then bare-name fallback, in try/catch) BEFORE
`dlopen(LIBAIC, RTLD_NOW|RTLD_GLOBAL|RTLD_DEEPBIND)` and `dlsym` of the two shim
symbols into `Ref{Ptr{Cvoid}}`. This is the su2-fft fix for the Julia
bundled-libgmp ABI clash with FLINT: `RTLD_DEEPBIND` on libaic alone is
insufficient because libflint resolves `__gmpn_*` against the already-loaded
Julia-bundled libgmp; pre-loading the system libgmp keeps its definitions in the
global symbol table first. On this system (Julia 1.12.5, FLINT 18, libgmp.so.10
at `/lib/x86_64-linux-gnu/`) the flags worked on the FIRST try — no further ABI
work needed.

**Marshalling.** `choi_diff(kraus, n; prec=106)` flattens the `r` Kraus matrices
into the C row-major layout `kraus[a*n*n + i*n + j]` (1-based Julia -> 0-based C
explicitly, NOT a column-major reinterpret), `ccall`s the shim, and unflattens
`choi[p*N + q]` (`N=n*n`) into a Julia `N x N` `ComplexF64`. `eta_eigfree(kraus)`
returns the rigorous `(lo, hi)` bracket. `eta_idempotence(kraus)` =
`diamond_norm_watrous(choi_diff(kraus, n), n)`, asserting SDP status `OPTIMAL`.

**Cross-check ladder (`test/runtests.jl`, 14/14 pass, ~52 s).**
- (T1) ccall ROUND-TRIP: `choi_diff(kraus, n)` (C/arb via ccall) vs a pure-Julia
  `Choi(Phi^2-Phi)` (`choi_A` + `compose_self`, copied from the fixture
  generator) for `identity(2)` and `phi_t(2,0.3)` — worst max-entry diff
  `8.33e-17`. Exercises the marshalling end-to-end.
- (T2) ANALYTIC anchors: `eta_idempotence` vs the EXACT closed forms —
  `id(2) -> 0` (got `0.0`); `phi_t(2,0.3) -> 0.315` (got `0.31499999928`);
  `phi_t(2,0.1) -> 0.135` (got `0.13499998256`); paper example `0.1 ->
  0.1*sqrt(0.9) = 0.0948683` (got `0.09486832704`). All within `1e-6`.
- (T3) the certified eig-free bracket BRACKETS the SDP value across the ccall
  boundary: `lo <= eta_idempotence <= hi` for `phi_t(2,0.3)`, `phi_t(2,0.1)`,
  `paper_ex(0.1)`, `dep_2`, `phi_t(3,0.2)`. Tightest below-margin
  `eta - lo = 0.0278` (paper_ex_0p1), tightest above-margin `hi - eta = 0.173`
  (paper_ex_0p1) — loose by design (ratio `~2n`), it certifies the MOSEK value
  rather than competing with it.
- (T4) GADC (live, non-fixture, from `tools/smoke_gadc_cbnorm.jl`):
  `eta_idempotence` at `gamma=0` (`2.9e-28`) and `gamma=1` (`7.3e-29`) ~0
  (the two idempotent endpoints), interior `gamma=0.5 -> 0.332 > 1e-3`.

---

## Module `idemp_structure` — structure of exactly-idempotent UCP maps (bead aic-wuh)

Realizes `th_idemp_structure` (`.tex:318` statement, `.tex:2055` proof), the
eta=0 oracle (milestone aic-9kk) — the cleanest ground truth in the project. The
proof IS the algorithm. Supporting: `lem_idemp` (`.tex:1916`), `lem_carrier`
(`.tex:1724`), `PhiX_M` (`.tex:1740`), `prop_Gamma` (`.tex:2106`, deferred).

Files: `include/aic_idemp.h`, `src/aic_idemp_internal.h`, plus six cores
(carrier/image on the double path, the rest arb): `aic_idemp_carrier.c` (118),
`aic_idemp_image.c` (108), `aic_idemp_maps.c` (the orchestrator),
`aic_idemp_build.c` (steps 5-7), `aic_idemp_verify.c`
(GCD/DGC), `aic_idemp_verify2.c` (w-hom/block-diag), `aic_idemp_cp.c`
(Gamma-unital + Psi-Choi + the w-Gamma=C_M-Lambda tie-check). Tests:
`tests/test_idemp.c` (channels in `tests/test_idemp.h`).

### The theorem and the construction (the proof is the algorithm, `.tex:2055-2091`)

Given an idempotent UCP self-map `Phi : B(H) -> B(H)` (`dim H = n`):

1. **Carrier `M = range(Q)`**, `Q = sum_a K_a K_a^dag` (`lem_carrier .tex:1724`;
   reuses `aic_ucp_carrier_Q`). `J_M` = eigenvectors of `Q` with eigenvalue
   `> thr = n*eps_mach*||Q||_F` (`n x dim_M` isometry), `J_Mperp` = the rest,
   `Pi_M = J_M J_M^dag`. **`M = H` iff `Phi` is trace-preserving** (`Q = 1`); `M`
   can be a PROPER subspace when `Phi` is unital but not TP (the `compress_idemp`
   oracle, `dim_M = 2 < n = 4`). Double path (`LAPACKE_zheev`): `Q`'s spectrum is
   degenerate. Fails loud on a non-PSD `Q` or an eigenvalue straddling `thr`.

2. **`A = Img Phi`** via the COLUMN-IMAGE SVD (not the `n^2 x n^2` superoperator
   — that risks a vec-convention bug): apply the tested `aic_ucp_apply` to each
   matrix unit `E_ij`, stack `vec(Phi(E_ij))` (row-major, `vec(Y)[a*n+b]=Y[a,b]`)
   as columns of `P` (`n^2 x n^2`), SVD, keep the left singular vectors with
   `sigma > thr`. `dim_A` = numerical rank; `B_k` = reshape of `u_k`. `Phi(B_k) =
   B_k` (asserted — fixed points). Double path (`LAPACKE_zgesvd`; degenerate
   singular spectrum).

3. **`Delta : A -> B(H)`** = inclusion = `[vec B_1|...|vec B_d]` (`n^2 x dim_A`),
   ORTHONORMAL columns (`Delta^dag Delta = 1_A`), so `Im Delta = Img Phi`.

4. **`C_M : X |-> J_M^dag X J_M`** (compression).

5. **`Lambda : B(M) -> B(H)`**, `Lambda(Y) = Phi(J_M Y J_M^dag)` (`.tex:2061`); CP
   + unital. Stored as `Lambda_mat` (`n^2 x dim_M^2`).

6. **`Gamma : B(M) -> A`** via `Lambda = Delta Gamma` (`.tex:2065`). Since
   `Delta^dag Delta = 1_A`, `Gamma_mat = Delta^dag Lambda_mat` (NO pseudoinverse;
   the paper guarantees this `Gamma` is UCP, `.tex:2088`).

7. **`w = C_M Delta`** (`.tex:2084`), `w_mat` column `k = vec(J_M^dag B_k J_M)`; a
   `*`-homomorphism (`lem_idemp`).

Maps are stored as their matrix in the row-major vec convention; a map
`B(C^p)->B(C^q)` is a `(q*q)x(p*p)` matrix.

### Constructivization note (Law 3)

The paper's proof is already constructive in finite dimensions: `M` is the
support of `Phi^*(rho_0)` (an existence statement, but `= range(Q)` computably),
and `A`, `Delta`, `Gamma`, `w` are the explicit maps above. The only
non-elementary step is the rank-revealing eig/SVD, handled by the degeneracy-
robust LAPACK double path (certified degenerate eig/SVD is blocked on aic-w4o.1;
the arb path here CERTIFIES the relations given the double-path subspaces, it does
not re-extract).

### Verification — the eta=0 oracle (`.tex:2080-2090`)

Certified arb defects, each `== 0` for an exact idempotent; maps are applied to
matrix units / basis coords and compared via `aic_mat_opnorm` (the `n^2 x n^2`
superoperator is NEVER formed — convention-safe). Coordinate maps use that
`Delta` has orthonormal columns: the `A`-coordinates of `Z in A` are
`Delta^dag vec(Z)`.

- `||Gamma C_M Delta - 1_A||_op` (`.tex:2081`),
- `max_ij ||Delta Gamma C_M(E_ij) - Phi(E_ij)||_op` (`.tex:2072`, `= Phi`),
- `max_jk ||w(B_j star B_k) - w(B_j)w(B_k)||_op`, `star = Phi(B_j B_k)` the
  Choi-Effros product (`lem_idemp`: `w` is a `*`-homomorphism),
- `max_k (||(1-Pi_M)B_k Pi_M|| + ||Pi_M B_k(1-Pi_M)||)` (block-diagonal, `.tex:2090`),
- `||Delta Gamma(1_M) - 1_H||_op` (Gamma unital).
- **Gamma CP — a PAIR of certificates** (`.tex:2084`, `.tex:2088`). The earlier
  single Choi check was VACUOUS for the stored `Gamma`: `Psi = C_M Lambda` is a
  composition of manifestly-CP maps, so its Choi is PSD by construction from
  `Lambda`, `J_M` alone — it never reads `Gamma` (zeroing `Gamma` left it green).
  We now certify two things:
  (i) the `C_M o Lambda` FACTORIZATION is CP — the abstract algebra `A`'s
      Choi-Effros product is CP — via the Choi matrix of `Psi : B(M)->B(M)` PSD
      (`aic_idemp_psi_choi` + `aic_ucp_is_cp_choi`), and
  (ii) the STORED `Gamma` is tied to that CP object:
      `||w Gamma - C_M Lambda||_op = 0` (`aic_idemp_defect_wG_eq_CML`, reads BOTH
      `w` and `Gamma`). Since `w` is an injective `*`-monomorphism (completely
      isometric, `.tex:2088`) and `Psi = w Gamma` is CP, `Gamma` is CP. A corrupt
      `Gamma` (scaled by 2, or a column zeroed) makes defect (ii) `O(1)` — RED.
      Working through `Psi` on `B(M)` is gauge-invariant.

Measured at `prec=53`, all five defects are `<= 5.96e-16` across all six oracle
channels (most exactly 0): block_cond_exp(5,2) `(dim_M,dim_A)=(5,13)`,
trace_replace(3) `(3,1)`, noiseless_subsystem(2,2) `(4,4)`, identity(3) `(3,9)`,
dephasing(4) `(4,4)`, compress_idemp(4,2) `(2,4)` [M proper]. Dim counts all match
the closed-form structure.

### Gauge freedom (load-bearing for cross-checks)

The `A`-basis `{B_k}` is unique only up to a unitary on `C^{dim_A}`. Double-vs-arb
/ cross-run comparisons compare the SUBSPACE via `Pi_A = Delta Delta^dag`, never
basis vectors. The five relations are gauge-invariant and asserted directly.

### Mutation proofs (Rule 7)

Three source mutations each turned a relation RED (then restored): (i) zeroing
`Lambda` -> GCD fails; (ii) scaling `Gamma` by 2 -> GCD fails; (iii) taking `J_M`
from the dropped (small-eigenvalue) `Q` block -> GCD fails on the M-proper case
(when `M = H` the block choice is immaterial — only the proper-carrier instance
catches a wrong carrier). The in-suite teeth test confirms a non-idempotent
`Phi` (`0.5 X + 0.5 diag(X)`) aborts at the entry guard.

### Cross-check ladder coverage (Rule 6)

1. Internal sanity: `J_M^dag J_M = 1_M`, `Delta^dag Delta = 1_A`, `Pi_M^2 = Pi_M`,
   `Phi(B_k) = B_k`.
2. acb@53 vs acb@256: `Pi_A` and the relation defects agree (the extraction is
   double-path regardless of prec, so this is the precision-self-consistency rung).
3. eta=0 oracle: the five relations `== 0` to `~1e-15`.
4. dim counts vs the known closed-form structure for each oracle channel.

### Precision

The relations are differences of `O(1)` quantities, so `prec=53` resolves them to
`~1e-15`. No near-singular inversion arises here (the `eta>0` machinery's
`(L+R)^{-1}` and the `O(eps)` cancellations belong to the §12 modules). The entry
guard rejects any `Phi` with `max_ij ||Phi^2(E_ij)-Phi(E_ij)||_op > 1e-9`.

### Deferred (beaded)

- Artin-Wedderburn block decomposition of `A`, `prop_Gamma`'s explicit `Tr_{E_j}`
  conditional-expectation form, and Kraus reps of `Delta`/`Gamma` -> bead aic-ynu.
- Certified degenerate eig/SVD for the extraction -> bead aic-w4o.1.
- `eta > 0` (almost-idempotent, via `Phi_tilde`) is `assoc_ecsa`'s job (the entry
  guard rejects it here).

## Module `ecstar` — ε-C* algebra data model + axiom-defect estimators (bead aic-knm)

Files: `include/aic_ecstar.h`, `src/aic_ecstar.c` (data model + star +
projection residual), `src/aic_ecstar_assoc.c` (the trilinear associator),
`src/aic_ecstar_defect.c` (submult, C*, involution, unit), shared helper
`src/aic_ecstar_internal.h`; tests `tests/test_ecstar.c`.

### Data model

An ε-C* algebra is a subspace `A ⊆ M_n` (`dim A = d`) given by a Frobenius-
ORTHONORMAL operator basis `{B_1,…,B_d}` (`<B_j,B_k>_F = δ_jk`, each `n×n`),
together with a unital CP map `Φ : B(C^n) → B(C^n)` (an `aic_ucp_kraus`,
`dim_K == dim_H == n`). The multiplication is the **Choi–Effros star**

```
X ⋆ Y = Φ(X·Y)            (approximate_algebras.tex:341-342, :2187-2189)
```

with `X·Y` the ORDINARY matrix product and `Φ` applied via `aic_ucp_apply`. The
involution is the matrix adjoint `X ↦ X†`; the unit is `I = 1_n` (the ambient
identity, inherited, .tex:2186-2187); the norm is the OPERATOR norm
`aic_mat_opnorm` (inherited from `B(H)`, .tex:2192). The cb-norm appears ONLY in
the definition of `η` (.tex:347) and is NOT used inside the estimators.

`aic_ecstar_from_idemp` fills `B_k[i,j] = d->Delta[i*n+j, k]` (reshape column `k`
ROW-MAJOR, the project's vec convention), the exact `η=0` case `A = Img Φ`; it
ASSERTS the columns are Frobenius-orthonormal on entry (the inherited
`Delta† Delta = 1_A`). `phi` is BORROWED (caller-owned, must outlive the struct).

### The estimators (basis sweeps; certified arb balls)

Each is a sweep over the stored basis using `aic_mat_opnorm`; each is the eta=0
EXACT-ZERO detector / a basis-sweep LOWER bound on the true sup-over-unit-ball
`ε`. The faithful worst-case search (HOPM) and the certified SDP upper bound are
later cycles (bead aic-0at).

| estimator | axiom (.tex) | formula (basis sweep) | η=0 reading | kind |
|---|---|---|---|---|
| `assoc` | ax_assoc :412-413 | `max_ijk ‖(B_i⋆B_j)⋆B_k − B_i⋆(B_j⋆B_k)‖_op` | ~0 | EXACT zero-detector (associator is trilinear ⟹ basis triples ⟺ identically) |
| `submult` | ax_prodnorm :410-411 | `max_jk max(0, ‖B_j⋆B_k‖_op − ‖B_j‖_op‖B_k‖_op)` | ~0 | LOWER bound (clamped ≥0) |
| `cstar` | ax_C* :427-428 | `max_k max(0, 1 − ‖B_k†⋆B_k‖_op / ‖B_k‖_op²)` | ~0 | LOWER bound (clamped ≥0); fails loud if `‖B_k‖_op < 1/(2√n)` |
| `involution` | ax_* :422-423 | `max_jk ‖(B_j⋆B_k)† − B_k†⋆B_j†‖_op` | 0 | ALWAYS-zero invariant for any Hermicity-preserving Φ (.tex:2208/:2211) |
| `unit` | ax_eps_unit :432-434 | `max( ‖1_n−Π_A(1_n)‖_op, max_k‖B_k⋆I−B_k‖_op, max_k‖I⋆B_k−B_k‖_op )` | 0 | ALWAYS-zero invariant for unital Φ, `A=Img Φ` (.tex:2211 `X⋆I=X=I⋆X`) |

A structural identity, recorded as a FINDING: **a UCP-defined star is always
submultiplicative** (`‖Φ(XY)‖ ≤ ‖Φ‖_cb‖XY‖ ≤ ‖X‖‖Y‖`, `‖Φ‖_cb=1`), so
`ax_prodnorm` holds with `ε ≤ 0` for ANY UCP `Φ` — the `submult` defect cannot be
made positive by mixing UCP maps; its teeth require a NON-contractive map.

### Cross-check ladder coverage (measured, prec=53)

Rung 3 (η=0 oracle) — all five defects per channel, run via
`aic_idemp_decompose → aic_ecstar_from_idemp`:

```
block_cond_exp(5,2)     assoc 0        submult 0  cstar 0  invol 0  unit 0
block_cond_exp(3,1)     assoc 0        submult 0  cstar 0  invol 0  unit 0
trace_replace(3)        assoc 0        submult 0  cstar 0  invol 0  unit 5.08e-16
noiseless_subsystem(2,2) assoc 0       submult 0  cstar 0  invol 0  unit 3.01e-16
identity(3)             assoc 0        submult 0  cstar 0  invol 0  unit 0
dephasing(4)            assoc 0        submult 0  cstar 0  invol 0  unit 0
compress_idemp(4,2)     assoc 0        submult 0  cstar 0  invol 0  unit 1.67e-16
```

The measured oracle ceiling is `≤ 5.08e-16` (unit, trace_replace(3)); assoc /
submult / cstar read EXACTLY 0 on every channel. Gate: involution + unit at
1e-12, assoc/submult/cstar at 1e-9 (the latter route Φ through the LAPACK-
extracted Δ-basis; in practice they stay machine-zero, well inside both gates).

Rung 2 (arb@53 vs arb@256): max defect diff `0.000e+00` on block_cond_exp(5,2).

Gauge invariance (a real plane-rotation chain on the `d` basis vectors,
`B'_k = Σ_l U[l,k] B_l`): max defect diff `1.923e-16` — confirms no estimator
secretly depends on the basis choice.

### Mutation / teeth map (Rule 7 — mandatory; the project's #1 bug is a test that cannot fail)

Each estimator is shown to MOVE under a deliberate non-C* perturbation:

- **`Φ_t = (1−t)Φ + t·Dep`** (UCP + Hermicity-preserving mix; Dep = trace-replace),
  basis from the exact `A=Img Φ`, star uses `Φ_t`. `B_k` is no longer a `Φ_t`
  fixed point, so:
  ```
  t=1e-3:  assoc 2.50e-4  cstar 7.50e-4  unit 1.00e-3   (submult 0, invol 0)
  t=1e-2:  assoc 2.48e-3  cstar 7.50e-3  unit 1.00e-2   (≈10× each)
  ```
  `involution` stays 0 (Hermicity preserved); `submult` stays 0 (UCP ⟹
  submultiplicative — the structural finding above).
- **`submult`** teeth need a NON-contractive map `Φ' = (1+δ)Φ` (scale each Kraus
  op by `√(1+δ)`): `submult` reads `0 → 1.00e-3 (δ=1e-3) → 1.00e-2 (δ=1e-2)`, 10×.
- **basis off A**: push `B_0` of the dephasing algebra by `t·E_{01}` (off the
  diagonal subalgebra): `unit` term (b) `‖Φ(B_0)−B_0‖` reads
  `0 → 1.00e-3 → 1.00e-2`, 10×.
- **projection residual / unit term (a)**: a 1-d algebra `span{E_00}` (Φ=id) does
  NOT contain `1_2`, so `aic_ecstar_proj_residual(1_2) = 1.000` and the unit
  defect = `1.000` — proves term (a) is non-vacuous.
- **`involution`**: `aic_ucp_kraus` can only encode Hermicity-preserving maps
  (`Φ(X)=ΣK†XK` is HP for any Kraus set), so the production estimator is a
  STRUCTURAL always-zero. The mutation-proof recomputes the same defect with a
  synthetic NON-HP map (doubling the imaginary part of one operand): the metric
  jumps `0.00 (HP) → 5.00e-1 (non-HP)`, confirming the adjoint bookkeeping is real.
- **plain-XY star mutant** (drop `Φ` from `aic_ecstar_star`): the η=0 oracle STILL
  PASSES (the exact-idempotent basis is closed/associative under plain `XY` too,
  lem_idemp), but the `Φ_t` teeth go RED (nothing routes through `Φ`). This is
  exactly why Part B (teeth) is mandatory: the oracle alone does not pin the
  Choi–Effros product.

### Precision

The defects are differences of `O(1)` operator-norm quantities, so `prec=53`
resolves them to `~1e-16` (oracle ceiling 5.08e-16) and arb@53 vs arb@256 agree
to `0`. No near-singular inversion arises in the basis sweep (the `(Λ+R)^{-1}` and
`O(ε)` cancellations belong to the §12 modules). `cstar` fails loud if any
`‖B_k‖_op` drops below `1/(2√n)` (a Frobenius-unit op has `‖·‖_op ≥ 1/√n`).

### Cycle 2 — HOPM faithful worst-case search (the second audition candidate)

Files: `src/aic_ecstar_hopm.{c,h}` (double-path matrix kernel: Φ-apply, star,
associator, opnorm, u/v power-step, project-to-A, polar), `src/aic_ecstar_iterate.{c,h}`
(restart init + the block-update with the accept-guard), `src/aic_ecstar_search.{c,h}`
(the multi-restart engine + `objval`), `src/aic_ecstar_setup.{c,h}` (double snapshot
of the `aic_ecstar` + scratch + the six defect-map term thunks),
`src/aic_ecstar_certify.c` (public API + arb certification). All ≤200 LOC.

Public API (`include/aic_ecstar.h`):
`aic_ecstar_defect_{assoc,submult,cstar}_hopm(arb_t lo, A, n_restarts, n_iters, prec)`
return a RIGOROUS LOWER BOUND `lo` on the true sup over the operator-norm unit
ball of `A`; `aic_ecstar_defect_assoc_hopm_witness` additionally exposes the
witness (for the soundness test).

| defect | objective (sup over op-norm-1 in A) | .tex |
|---|---|---|
| `assoc`   | `‖(X⋆Y)⋆Z − X⋆(Y⋆Z)‖_op` | :412-413 |
| `submult` | `‖X⋆Y‖_op − 1` (expected ~0: UCP star is submultiplicative) | :410-411 |
| `cstar`   | `1 − min_{‖X‖_op=1} ‖X†⋆X‖_op` (a MINIMIZATION inside) | :427-428 |

**The method (web leg "Cycle 2 method decision", Method 1).** Scale-invariant
alternating maximization. `X = Σ_k x_k B_k` so every iterate is automatically in
`A` (valid witness); `‖h‖_op = max_{‖u‖=‖v‖=1}|⟨u,hv⟩|` introduces auxiliary unit
vectors. Per inner iteration: (1) `u ← hv/‖hv‖`, `v ← h†u/‖h†u‖` (closed form);
(2) per block, form the A-gradient `C = Σ_k ⟨u, h(B_k,·,·) v⟩ B_k` (in A), then the
candidate `X' = Π_A(polar(C))` re-normalized to op-norm 1, **accepted only if the
scale-invariant objective strictly improves**. `cstar` minimizes (`sign=−1`,
accept on decrease); `submult` drops the Z block. Deterministic restarts (LCG
seeded from the restart index — no wall-clock RNG); half warm-start from a
polar-projected random A-op, half plain random in A.

**The load-bearing SUBSPACE-POLAR subtlety + accept-guard.** `polar(C)=UV†` is the
exact maximizer of `Re⟨C,X⟩` over the FULL `M_n` op-norm ball, but a *genuine* C\*
algebra is polar/von-Neumann closed while our `A` is only an *ε*-C\* algebra
(η>0), so `polar(C) ∉ A` in general. Using it directly would produce a witness
OUTSIDE A → an INVALID lower bound. Hence `Π_A(polar(C))` + the monotone-ascent
accept-guard: this keeps every accepted iterate exactly in A AND makes the
approximate-polar step robust. (Confirmed: the soundness test recomputes
`proj_residual(witness) < 1e-30`.)

**arb-vs-double split + certification.** The search runs entirely in the fast
double path (LAPACKE SVD via `latd`, C99 `double _Complex` Φ-apply on ball
midpoints) — decisions use midpoints. The FINAL witness (explicit, in A) is
re-evaluated ONCE in the certified arb path (`aic_ecstar_star` + `aic_mat_opnorm`),
and `lo` is `‖num‖_op.lo / (Π ‖den‖_op.hi)` (arb rounding directions chosen so the
quotient provably under-estimates the value at the witness ≤ the sup). So `lo` is
a rigorous lower bound regardless of how the witness was found (web leg Q8).

**Timing (prec=53, this box).** Search cost ≈ `O(n_restarts · n_iters · nblk · d · r · n³)`
(the `d·r` factor is the per-block gradient looping over the `d` basis × `r` Kraus).
Convergence near η=0 is fast: 10 restarts × 20 iters already reaches the converged
value on bce(4,2) (0.35s). The full `test_ecstar` (Cycle 1 + Cycle 2, 109 checks)
runs in ~12s on this box; the HOPM calls themselves are 0.1-1s each (the d=9
dim_A=45 dephasing HOPM is 0.9s). NOTE: the slow part of the canary is the
Cycle-1 BASIS sweep, not HOPM — `aic_ecstar_defect_assoc` is `O(dim_A³)` arb star
products (91125 at dim_A=45 ≈ 190s), so the canary skips it at the large-dim
point. The dominant HOPM
cost is the per-`k` gradient recomputing the full associator; a `Φ†`-adjoint
gradient (no `k` loop) is the obvious future speedup (deferred — it is exactly the
kind of adjoint bookkeeping where a sign/transpose bug hides, Rule 3).

### The audition verdict (Law 4 — both candidates kept, Pareto)

| | basis sweep (Cycle 1) | HOPM (Cycle 2) |
|---|---|---|
| what | max over basis triples/pairs/singletons | scale-invariant alternating max over the op-norm unit ball of A |
| bound | exact-zero detector / loose LOWER bound | faithful LOWER bound (tighter) |
| cost | `O(d³ n³)` one pass, no iteration | `O(restarts·iters·d·r·n³)` |
| dim-independence | can DRIFT with dim (`√d`/`d^{3/2}`, web §2) | flat — works on the op-norm sphere, no Frobenius inflation |
| certified | yes (arb) | yes (arb witness re-eval) |

Neither dominates: the basis sweep is the **cheap zero-detector** (η=0 oracle in
one pass, no search); HOPM is the **faithful, dimension-independent worst-case
lower bound** but iterative/slower. Both are kept and dispatchable. Measured
audition payoff and canary (prec=53):

```
AUDITION (bce(4,2) + trace_replace mix):
  t=1e-3   assoc basis_sweep=2.498e-4   HOPM_lo=1.449e-3   (HOPM 5.80x BS)
  t=1e-2   assoc basis_sweep=2.481e-3   HOPM_lo=1.443e-2   (HOPM 5.81x BS)

CANARY (A) dimension sweep, dephasing-mix, assoc/t, fixed t=1e-2:
  bce(4,2) dim_A= 8   HOPM/t=1.000
  bce(6,3) dim_A=18   HOPM/t=1.319
  bce(9,3) dim_A=45   HOPM/t=1.302     -> spread factor 1.32 (FLAT, dim-independent)

CANARY (B) contrast, trace_replace-mix, assoc/t, fixed t=1e-2:
  bce(4,2) dim_A= 8   basis_sweep/t=0.2481   HOPM/t=1.4030
  bce(5,2) dim_A=13   basis_sweep/t=0.1984   HOPM/t=1.3331
  -> basis_sweep/t DRIFTS DOWN 20% with dim; HOPM/t stays flat (~1.3-1.4)
```

The HOPM/t ratio staying within a factor 1.32 across `dim_A ∈ {8,18,45}` is the
universality canary passing (HANDOFF "the single highest-value test"); the
trace_replace contrast shows the basis sweep's ratio drifting where HOPM does not.
(The contrast's second point is `bce(5,2)`, dim_A=13, NOT d=9: the basis-sweep
assoc is `O(dim_A³)` arb star products — 91125 at dim_A=45, ~190s — so the
canary skips it at the large-dim point, where only HOPM/t flatness is asserted.
The `O(dim_A³)` arb basis sweep is the practical ceiling on the basis-sweep rung,
a finding in its own right.)

η=0 oracle (HOPM, prec=53): assoc_lo `= 0` exactly; submult_lo / cstar_lo in
`[-1.4e-14, -7e-16]` (slightly negative = the genuine C\* algebra has zero/negative
defect) on the 7+1 channels — all under the `1e-7` gate.

Mutation-proof (Rule 7) — the weakest search (1 restart, 0 iters) finds STRICTLY
less than the full search on a perturbed instance, so a "return 0 / never iterate"
stub is RED:
```
SEARCHES assoc  weakest(1,0)=4.698e-3  full=1.443e-2   (also full >= basis-sweep witness 2.481e-3)
SEARCHES submult weakest(1,0)=-9.476e-3 full=-6.286e-4
SEARCHES cstar  weakest(1,0)=-5.1e-15  full=6.102e-3
```
Soundness: `lo = ratio@witness = 1.442605e-2`, witness `proj_residual < 1e-30`.

### Deferred (beaded)

- The certified SDP UPPER bound (Watrous cb-norm SDP for the bilinear
  submult/cstar; Lasserre/SOS lift for the trilinear assoc) -> bead aic-0at.
  Cycle 1 (basis sweep) + Cycle 2 (HOPM) bracket the sup from below; the SDP is
  the matching upper bracket.
- A `Φ†`-adjoint gradient for the block update (eliminates the per-`k` loop -> the
  main speedup for larger `n`,`d`). Filed as a follow-up bead.

## Module `assoc_ecsa` (Increment 1) — regularization Φ̃ = θ(2Φ−1) (bead aic-92f)

STEP 1 of `th_almost_idemp` (`approximate_algebras.tex:2162-2237`): regularize an
η-idempotent UCP map `Φ` to an EXACTLY idempotent superoperator `Φ̃`. Increment 1
builds the superoperator matrix `S_Φ`, the regularization `Φ̃ = θ(2Φ−1)`, and the
apply helper; `A = Img Φ̃` and the Choi-Effros star are Increment 2.

### Paper technique vs constructive route

The regularization is (`.tex:2171-2174`, eq `tilde_Phi`):

```
Φ̃ = θ(2Φ − 1) = (1/2)(1 + sgn(2Φ − 1))
              = (1/2)(1 + (2Φ − 1)(1 − 4(Φ − Φ²))^{−1/2}).
```

The paper frames this as the functional calculus of the **cb-norm Banach algebra**
of completely bounded maps on `B(H)` (`.tex:354-359`) — a Taylor expansion in
`4(Φ − Φ²)` converging for `η < 1/4`. This is Proposition `prop_P`
(`.tex:524-533`: `‖P²−P‖ ≤ δ < 1/4 ⇒ P̃ = θ(2P−I)` is an exact idempotent) applied
to `P = Φ`. We do **not** transcribe the infinite-dim cb-algebra calculus. In
finite dimensions a UCP self-map on `M_n` IS an `n²×n²` matrix `S_Φ`, so
`θ(2 S_Φ − 1)` is an ORDINARY matrix-functional-calculus call, served by the
eig-free certified `funcalc` engine: `aic_prop_P` → `aic_theta` → `aic_sgn`
(Newton-Schulz). The paper's bound `‖Φ̃−Φ‖_cb ≤ O(η)` (`.tex:2179`) is the spec;
T7 asserts a computable proxy of it (linear-in-defect growth, vanishing at η=0,
dimension-independent constant).

### Superoperator + the row-major oracle

vec convention (matches `aic_idemp.h`, `aic_ecstar.h`): `vec_r(X)[i·n+j] = X[i,j]`;
the superoperator `S` of `T : M_n → M_n` has `S[a·n+b, p·n+q] = T(E_{pq})[a,b]`, so
column `(p·n+q)` of `S` is `vec_r(T(E_{pq}))`.

- PRIMARY build `aic_assoc_superop_from_ucp` (lowest vec-bug risk): column-by-column,
  reshaping `aic_ucp_apply(Φ, E_{pq})` row-major. Reuses the TESTED apply + reshape
  (idemp/ecstar precedent).
- SECONDARY build `aic_assoc_superop_kron` (the cross-check): the index-derived
  `S_Φ = Σ_a K_a^† ⊗ K_a^T` (left-major Kronecker, left `K_a^†`, right `K_a^T`).
  Derivation: `Φ(X)[a,b] = Σ_a Σ_{p,q} (K_a^†)[a,p] X[p,q] (K_a^T)[q,b]`, so
  `S[a·n+b, p·n+q] = Σ_a (K_a^†)[a,p] (K_a^T)[b,q] = Σ_a (K_a^† ⊗ K_a^T)[a·n+b,p·n+q]`.

The vec convention is the project's #1 historical bug (idemp deliberately AVOIDED
forming the superoperator). So `S_Φ` is oracle-tested two independent ways:
- **T1 oracle:** `S_Φ vec_r(X) == vec_r(aic_ucp_apply(Φ,X))` for 3 fixed `X` (complex
  non-Hermitian, shifted, Hermitian) on ≥2 channels incl. a complex/asymmetric one
  (conjugated dephasing). Machine-zero. Mutation-proven: a **transposed reshape**
  (store `Φ(E_{pq})[b,a]`) → RED at entry `(0,1)`.
- **T2 Kronecker:** column-image `S_Φ == Σ_a K_a^† ⊗ K_a^T`, machine-zero.
  Mutation-proven: **swapped Kronecker factors** (`K_a^T ⊗ K_a^†`) → RED at `(0,1)`.

### The η=0 reduction (headline oracle, T3)

If `Φ² = Φ` exactly then `Φ − Φ² = 0 ⇒ (1−4(Φ−Φ²))^{−1/2} = I`, so
`Φ̃ = (1/2)(1 + (2Φ−1)) = Φ`; equivalently `(2S_Φ−1)² = I` exactly ⇒
`sgn(2S_Φ−1) = 2S_Φ−1`. **T3:** for each exact-idempotent test channel
`S_Φ̃ == S_Φ`. Subtlety (a real finding): some "exact" idempotents carry
IRRATIONAL Kraus entries (`trace_replace`, `noiseless_subsystem` use `1/√k` set
from a `double`), so `Φ` is idempotent only to its own residual
`δ = ‖S²−S‖ ≈ 2.7e-16`; integer-Kraus channels (`block_cond_exp`, `identity`,
`dephasing`, `compress_idemp`) are exact (`δ ≈ 1e-30`). The honest T3 bound is
`prop_P`'s `‖S̃−S‖ ≤ ‖2S−I‖·O(δ)` with the input's TRUE `δ` (a defect-aware tol
= `100·ubound(‖S²−S‖_F) + 1e-22`, NOT a constant) — the channel reduces to its OWN
idempotence precision. Measured: `0.00e+00` (integer-Kraus), `2.69e-16`
(`trace_replace(3)`), `3.55e-16` (`noiseless_subsystem(2,2)`), `1.01e-75`
(conjugated dephasing) across 7 channels.

### Property proofs (T4, T5)

`.tex:2177-2182` — `Φ̃²=Φ̃`, `Φ̃(1)=1`, `Φ̃(X^†)=Φ̃(X)^†`.
- **T4 idempotency, EXACT (machine-floor `1e-22`):** `sgn²=I` exactly ⇒ `θ²=θ`
  exactly, INDEPENDENT of the input defect. Measured `‖S̃²−S̃‖_F ≤ 2e-74`.
  Mutation-proven: `aic_assoc_regularize` returning `S_Φ` unmodified → RED (caught
  on both `trace_replace` T4 and the dedicated non-idempotent T4 teeth, where the
  no-op leaves `‖S̃²−S̃‖_F = ‖S²−S‖ ≫ tol`).
- **T5 unit (defect-aware):** `Φ̃(I)=I` because `I` is a `Φ`-fixed point (up to δ).
- **T5 Hermicity, EXACT (machine-floor):** `θ` of an HP superoperator is HP, so
  `Φ̃(X^†)=Φ̃(X)^†` machine-zero regardless of δ.

### The non-normality cross-check (T6) — funcalc's first non-normal customer

`S_Φ` is NON-NORMAL (`Φ` is Hermicity-preserving but NOT HS-self-adjoint), and
`aic_sgn` (Newton-Schulz) was auditioned only on Hermitian/near-normal inputs.
T6 fills the gap AND is an independent algorithmic cross-check:
- **Route A1** = `aic_prop_P(S_Φ)` (Newton-Schulz `sgn`).
- **Route A2** = the explicit binomial series (`aic_funcalc_xpow`):
  `D = 4(Φ−Φ²) = 4(S−S²)`; `M = I − D = (2S−I)²`; `M^{−1/2} = xpow(M,−1/2,1)`;
  `S̃ = (1/2)(I + (2S−I) M^{−1/2})`. **Sign is load-bearing:** `M = 1−4(Φ−Φ²)`, so
  `D = 4(S−S²)` not `4(S²−S)`; the η=0 case (`D=0`) is blind to the sign, the η>0
  non-normal T6 channel PINS it (a flipped sign gives `sgn²≠I` ⇒ A1/A2 disagree by
  `O(0.1)` — this was caught during bring-up).
- TWO structurally-different non-normal in-basin channels (each runs the same
  `T6_nonnormal` helper, which asserts non-normality + in-basin live):
  - `compress_idemp(3,1)` mixed `0.05` with `dephasing(3)`: `‖[S,S^†]‖_F = 3.126`,
    `4‖S²−S‖_op = 0.329`, A1-vs-A2 `1.92e-74`.
  - `compress_idemp(4,1)` mixed `0.05` with `block_cond_exp(4,2)` (different base
    dim + a different, two-block mix family ⇒ a different non-normal structure):
    `‖[S,S^†]‖_F = 4.421`, `4‖S²−S‖_op = 0.380`, A1-vs-A2 `3.00e-74`.
- Newton-Schulz is CORRECT on non-normal input (it uses only matrix multiply/add,
  no normality assumption). Both channels mutation-proven: flipping the A2 sign
  (`D=S²−S`) makes A1≠A2 by `O(0.1)` ⇒ RED on each channel independently.

### The O(η) bound (T7) — a PRELIMINARY check, not the universality claim

`.tex:2179` gives `‖Φ̃−Φ‖_cb ≤ O(η)`, but **prop_P (`.tex:524-533`) only proves
`‖Φ̃−Φ‖ ≤ ‖2Φ−1‖·O(δ)` at the superoperator level**, and `‖2Φ−1‖` need NOT be
dimension-independent. So `‖Φ̃−Φ‖` dimension-independence is **NOT** the paper's
universality claim — that claim is about the **algebra's axiom-defect constant**
being dimension-independent, which is **Increment 2's canary** (`aic-3qq`/
`aic-4c7`), not this. T7 is a preliminary sanity check using the computable proxy
`η_proxy = ‖S²−S‖_op`, `dist = ‖S̃−S‖_op` (no cb-norm SDP for the trend), in two
honest parts:

- **T7a (linear growth + vanishing).** Smooth family
  `Φ_t = (1−t)·dephasing(3) + t·trace_replace(3)`, `t = 0.02..0.18` at ONE dim.
  `dist/η_proxy ∈ [1.02, 1.22]`, `dist` strictly increasing with `t` (asserted)
  and `→ 0` as `t → 0`. This is the genuinely-meaningful O(η) behavior, on a family
  where the op-norm proxy is smooth.
- **T7b (dimension-boundedness canary, strengthened).** The earlier
  `dephasing⊕trace` family was DIMENSION-TRIVIAL (byte-IDENTICAL ratios at every
  `d`; the hostile review showed it would pass even with a dim-dependent bug). T7b
  instead uses `compress_idemp(d,1) ⊕ dephasing(d)`, `d = 3,4,6,8`, `t = 0.02,0.04`
  (all in-basin), whose non-normality `‖[S,S^†]‖_F` GROWS with `d` (`3.3 → 10.2`),
  so `dist/η_proxy` is genuinely MEASURED (`∈ [0.52, 1.04]`; e.g. the `d=8,t=.02`
  ratio `0.52` ≠ the `d=3,t=.02` ratio `1.02`), not constant-by-construction. The
  ratio stays BOUNDED across dims (`≤ 5` gate, met at `~1`), i.e. does not grow
  with dim. Two INTEGRITY guards keep T7b from itself becoming toothless, BOTH
  mutation-proven (swap to `dephasing⊕trace` ⇒ RED): (i) `‖[S,S^†]‖_F > 1e-3` at
  every dim AND grows ≥2× from `d=3` to `d=8` (a NORMAL family has `‖[S,S^†]‖_F=0`
  and is rejected); (ii) the `d=3` and `d=8` ratios DIFFER by `> 1e-3` (a constant-
  by-construction family is rejected). The ratio bound and T7a monotonicity were
  also mutation-proven (tighten the gate / reverse the sweep ⇒ RED).

### Fail-loud basin (T9) + its dimension dependence

`θ`/`sgn` are well-defined when `2S−1` has no purely-imaginary eigenvalue
(`.tex:516,520-521,546-550`). The ORIGINAL `aic_prop_P` asserted the OPERATOR-norm
basin `‖S²−S‖_op < 1/4` (`.tex:525`), which TRIPS near `η → 1/4` or for large-`n` /
NON-NORMAL channels even though the SPECTRAL condition `ρ(S²−S) < 1/4` still holds.
**Bead aic-8hz (LANDED)** relaxed this: `aic_prop_P` now certifies the spectral
condition EIG-FREE via the Gelfand bound `‖M^k‖_F^{1/k} < 1` on `M = I−(2S−1)² =
4(S²−S)` (accept the first `k` whose certified ball is `< 1`; `k=1` is the old
Frobenius check), and `aic_sgn` auto-dispatches to a globally-convergent Newton
iteration `½(Y+Y⁻¹)` (Higham *Functions of Matrices* Thm 5.6) out of the op-norm
basin. The payoff is `test_assoc2` **U6**: a genuinely oblique channel
(`compress_idemp(4,1)×dephasing(4)`, `t=0.30`) with `‖S²−S‖_op = 0.420 ≥ 1/4` (OLD
`prop_P` aborts) but `ρ(S²−S) = 0.210 < 1/4` (Gelfand certifies `4ρ<1` at `k=6`,
`ρ_ub=0.975`); the relaxed `aic_assoc_regularize` succeeds, `‖S̃²−S̃‖_op = 8.9e-70`.
A sign witness pins the correct branch via unitality `Φ̃(1)=1` (`.tex:2179`):
`‖S̃·vec(I)−vec(I)‖ = 4.0e-70`, vs `√n = 2.0` for the `−sgn` complement (mutation-
proven). See `## Module funcalc — global sgn` below and `docs/research/sgn_global_research.md`.

### Precision ladder (T8) + the near-zero opnorm fix

Following funcalc's arb-only precedent (no double `sgn`), the prec=53-vs-256 rung
substitutes for double-vs-arb. **T8:** `Φ̃@53` vs `Φ̃@256` agree to `2.26e-14`
(within `1e-12` gate).

A near-zero-input fragility in `aic_mat_herm_max_eig` was found and ROOT-CAUSE
fixed here (Rule 3, not a bandaid). For a near-zero Hermitian `H` — e.g. the Gram
of `S²−S` for an EXACT idempotent (entries `~1e-31`, the η=0 oracle) —
`acb_mat_eig_global_enclosure` can return a NON-FINITE radius `ε`, which `sqrt`'d
in `aic_mat_opnorm` gives `nan±inf` and spuriously trips `prop_P`'s basin guard on
the cleanest ground-truth input. Fix: when `ε` is non-finite, fall back to the
rigorous eig-free bound `|λ| ≤ ‖H‖_op ≤ ‖H‖_F` ⇒ `λ_max ∈ [−‖H‖_F, ‖H‖_F]` (the
heuristic QR midpoint is not trusted when the certifier failed). Sound, always
finite, tight where used (near-zero `H` ⇒ near-zero `‖H‖_F`). This fallback now has
a DIRECT, mutation-proven correctness test in `tests/test_mat.c` (test 6): the
trigger (found empirically) is the Gram `H = D^†D`, `D = S²−S` of the exact
idempotent `noiseless_subsystem(2,2)` superoperator (entries `~1e-31`); plain
near-zero diagonals do NOT fire the branch (FLINT returns a finite `ε` for them),
so the degenerate Gram structure is what overflows. The test asserts (1) the
fallback FIRED (mid=0, radius `== ‖H‖_F = 6.287e-32`, fails loud if a future FLINT
stops triggering it), (2) the ball is SOUND (finite, radius `≥ ‖H‖_F`), (3) the
ball CONTAINS the true `λ_max = 3.144e-32` (LAPACK double-path oracle), which is
exactly `0.5·‖H‖_F`. Mutation-proven: `arb_zero(out)` (point-ball `[0,0]`) and
`0.25·‖H‖_F` (too-small radius) BOTH fail to contain `λ_max` ⇒ RED. The
inflated-radius opnorm
(funcalc iterations grow `S̃` radii to `~5e-75` at prec=256, exceeding the
Hermiticity tol `2^{−(prec−8)} ≈ 2.2e-75` when squared into the Gram) is handled in
the TREND measurements (T7, defect proxy) by `midpoint_opnorm` (strip radii first;
rigor lives in T3/T6, not the trend).

### Files + numbers

- `include/aic_assoc.h`, `src/aic_assoc_superop.c` (116 LOC), `src/aic_assoc_regularize.c`
  (61 LOC), `tests/test_assoc.c` (93 checks). The `aic_mat_herm_max_eig` fallback
  test lives in `tests/test_mat.c` (test 6; `test_mat` is now 23 checks).
- T3 η=0 oracle: max `‖S̃−S_Φ‖_F = 3.55e-16` across 7 channels at prec=256 (machine-
  zero for integer-Kraus). T4: `‖S̃²−S̃‖_F ≤ 2e-74`. T6 (two non-normal channels)
  A1-vs-A2 `1.92e-74` (`‖[S,S^†]‖_F=3.126`, `4‖S²−S‖_op=0.329`) and `3.00e-74`
  (`‖[S,S^†]‖_F=4.421`, `4‖S²−S‖_op=0.380`). T7a `dist/η ∈ [1.02,1.22]` (monotone);
  T7b `dist/η ∈ [0.52,1.04]` across `d=3,4,6,8` on a dim-nontrivial family
  (`‖[S,S^†]‖_F 3.3→10.2`). T8: `2.26e-14`.

### Deferred (beaded)

- Out-of-basin globally-convergent `sgn` (η near 1/4 / large n) → bead **aic-8hz**
  CLOSED 2026-05-30 (Gelfand spectral precondition + global Newton; see the T9
  section above and the dedicated funcalc-global module section).

## Module `assoc_ecsa` (Increment 2) — A = Img Φ̃, the Choi–Effros star, th_almost_idemp (bead aic-92f)

Increment 2 completes `th_almost_idemp` (`.tex:2184-2237`): build `A = Img Φ̃ =
Ker(1−Φ̃)` (`.tex:2185-2186`), the approximate Choi–Effros product `X⋆Y = Φ̃(XY)`
(`.tex:2189`), and CERTIFY that `(A, ⋆, ‖·‖, †, I)` is an extended `O(η)`-`C*`
algebra by running the `ecstar` axiom-defect estimators on the genuinely `η>0` `A`.

### A-extraction via oblique-projector SVD (Route-B1)

`Φ̃` is the `n²×n²` superoperator `S̃` (Increment 1). `A = range(S̃)` is extracted
(`src/aic_assoc_extract.c`) as a Frobenius-orthonormal operator basis `{B_k}`:

- `dim_A = round(Re Tr S̃)`. `S̃²=S̃` ⇒ eigenvalues `∈ {0,1}` ⇒ `Tr S̃ = #(unit
  eigenvalues) = rank = dim A` (an integer). CROSS-CHECKED against the thin-SVD gap:
  the count of singular values above `0.5` must equal `round(trace)` or **fail loud**
  (Rule 4). The gap genuineness is also asserted (`σ[dim_A−1] > 0.5 > σ[dim_A]`).
- `B_k = reshape_row-major(top-k LEFT singular vector of S̃)`. Left singular vectors
  are orthonormal in `C^{n²}` = Frobenius-orthonormal as operators (asserted,
  `‖⟨B_j,B_k⟩_F − δ_jk‖ ≤ 1e-9`).

**The projector-SVD subtlety (load-bearing).** `S̃` is idempotent (`S̃²=S̃`) but in
general NOT Hermitian (`Φ̃` is HP but not always HS-self-adjoint, Increment 1's
non-normality note). For an idempotent the nonzero SINGULAR values are `≥ 1`, with
EQUALITY (all `= 1`, an ORTHOGONAL projector, SVD spectrum exactly `{1,…,1,0,…}`)
when `Φ̃` IS HS-self-adjoint — the self-dual channels: **dephasing**, AND (a FINDING,
test U5) the η>0 **`dep(d)+conj(dep(d))`** family, since a convex mix of HS-self-adjoint
maps stays HS-self-adjoint, so its `S̃` is ORTHOGONAL with `σ_max = 1`. The sigmas are
STRICTLY `> 1` only when `Φ̃` is genuinely OBLIQUE (non-HS-self-adjoint), e.g. the
compression channel **`compress_idemp(4,2)`**, whose `S̃` inflates `σ_max = √3 ≈ 1.732`
(measured, U5). The earlier note that `dep+conj` is "oblique" was WRONG — its `σ_max`
is exactly 1; that family is the U2/U3/U4 *non-associativity* witness, which is a
SEPARATE property from obliqueness. Either way `rank = trace` still holds, `range(S̃)`
is correctly spanned by the top-`dim_A` LEFT singular vectors, and the `0.5` gap
separates nonzero (`≥ 1`) from zero (`~0`) cleanly in both cases. Extraction is the
fast DOUBLE path (LAPACK `zgesvd` via `aic_latd_svd`; certified degenerate eig/SVD is
the deferred `aic-w4o.1` wall — extraction is double, the DEFECT checks are
arb-certified), the basis returned as `acb_mat` at `prec`.

### The star + the ecstar generalization (the apply-fn seam)

`Φ̃` is a superoperator, NOT a UCP/Kraus map (`.tex:363` — not even CP), but the
`ecstar` star was hardwired to a Kraus `Φ` via `aic_ucp_apply`. The star is
generalized through the apply-fn seam `aic_ecstar.h` already shipped, MINIMALLY and
BACKWARD-COMPATIBLY: two optional fields `star_phi`/`star_ctx`. `aic_ecstar_star`:
if `star_phi != NULL`, compute `XY` then `star_phi(out, XY, star_ctx)`; else the
existing `aic_ucp_apply(A->phi, XY)`. `aic_ecstar_init` sets `star_phi=star_ctx=NULL`
so the Kraus path is byte-identical (test_ecstar n=109, **unchanged**, is the
regression guard). `aic_ecstar_defect_involution` routes through `star_phi` when set
(for the superop case `Φ̃` is HP, so the residual stays ~0, but it is MEASURED through
the generic core, not assumed). The HOPM double-path engine (`aic_ehk`) gains a
superop field `S` (built once from `star_phi` applied to the n² matrix units in
`aic_ecstar_setup.c`, fully decoupling the engine from the assoc module); `aic_ehk_star`
dispatches on `S != NULL`.

The assoc owner `aic_assoc_ecstar` (`src/aic_assoc_algebra.c`) holds the OWNED `S̃` +
ctx + the embedded `aic_ecstar` whose `star_phi` = a superop-apply thunk over `S̃`
(via the tested `aic_assoc_superop_apply`), `phi = NULL`. `aic_ecstar` only BORROWS
`star_phi/star_ctx` (mirroring its borrowed-`phi` rule), so the owner frees them in
`aic_assoc_ecstar_clear`.

**The midpoint-`S̃` decision (load-bearing, NOT a bandaid).** The regularized `S̃`
carries Newton–Schulz error balls (~5e-75 radius at prec=256). The defect estimators
apply `S̃` to a product and feed the result to `aic_mat_opnorm`, which forms the Gram
`M†M`; squaring a ~5e-75 radius exceeds `aic_mat_opnorm`'s strict Hermiticity tol
`2^-(prec-8) ≈ 2.2e-75` (the asymmetry is an `acb_mat_mul` artifact, not real
non-Hermiticity). So `A` is built from the MIDPOINT of `S̃` (a zero-radius idempotent
superop): the estimators then return certified balls of a well-defined exact-idempotent
superop, exactly as `B_k` is a zero-radius double-path extraction. The certified rigor
of the regularization itself is Increment 1's T3/T6 cross-checks (a separate ladder
rung); tightening to a radius-aware opnorm is the same deferral as Increment 1's
inflated-radius wall.

### The η=0 consistency oracle (U1) — the headline cross-check

For each of 8 exact-idempotent channels, `Φ̃=Φ`, so `A = Img Φ̃` must MATCH `ecstar`'s
`A = Img Φ` (`aic_ecstar_from_idemp`): same `dim_A`; same SUBSPACE (the gauge-invariant
projectors `Π_A = Σ_k vec(B_k)vec(B_k)†` agree — `A` has unitary gauge freedom, so
compared as subspaces, never operator-by-operator); ALL FIVE axiom defects (basis-sweep
+ assoc/cstar HOPM) machine-zero. **Measured:** subspace diff `0` (≤ 8.5e-16 for the
irrational-Kraus `noiseless_subsystem`); all five defects `≤ 1.69e-32` (machine-zero;
the `noiseless` cstar/unit floor `~3-6e-16` is the double-path SVD floor). This ties
Increment 2 to the trusted η=0 oracle.

**A genuine finding.** `A = Img Φ̃` with `X⋆Y = Φ̃(XY)` is an EXACTLY idempotent
structure, so its associator is exactly zero for ANY family whose `Φ̃` image happens to
be a real algebra. Most natural perturbation families (`compress+dephasing`,
`block_cond_exp+trace_replace`, `id+dephasing`) "snap" to a real `C*` algebra under
`θ(2Φ−1)` ⇒ associator ~0 even at η>0. A genuinely-nonzero associator needs two
DIFFERENT (incompatible) algebras mixed: `dep(d) ⊕ conj(dep(d))` (a diagonal algebra
and a unitarily-rotated diagonal) is the chosen `η>0` family — `Φ̃`'s image is then a
genuine NON-algebra and the associator is nonzero `O(η)`.

### The O(η) results (U2) + the universality canary (U3)

U2 (`dep(4)+conj(dep(4))`, t-sweep): the five defects are `≤ C·η`, vanish as `t→0`,
NONZERO for `t>0` (genuinely approximate — teeth). Both the associator AND the
`C*`-identity (`.tex:2233-2236`) are now ASSERTED at η>0 (FIX A): **measured** assoc
HOPM `3.7e-5 → 1.5e-3` (`t=0.01→0.06`, `assoc/η ∈ [0.005,0.037]`, bounded, monotone
up); cstar HOPM `6.0e-3 → 3.4e-2` — a GENUINE strong tooth (`Φ̃` is HP but NOT
positive, `.tex:363`, so `‖X†⋆X‖` drops `O(η)` below `‖X‖²` for the worst-case X),
with `cstar/η ∈ [0.75,0.82]` (bounded `≤ 1.5`, monotone up, vanishing at `t→0`, SAME
η proxy as assoc — apples-to-apples); submult slack stays ~0 (`≤ 2e-3`, the star is
submultiplicative, `.tex:2216-2219`); involution exactly 0 (`Φ̃` HP, `.tex:2182`);
unit `~1e-16`. The cstar assertion is mutation-proven: `aic_ecstar_defect_cstar_hopm`
→ `return 0` fires the `c > 1e-6` tooth RED.

U3 THE UNIVERSALITY CANARY (`aic-dbo.3`): at fixed `t=0.05` across a WIDE sweep
`dim_A = 4,6,8,10,12,16` (FIX B — widened from `4,6,8,10` so a dim-dependent constant
is amplified across a factor-4 span), HOPM `assoc/η` does NOT grow with `dim_A`.
**Measured `assoc/η ∈ [0.033,0.047]`, real spread factor 1.42** — flat,
dimension-independent (the naive Haar route fails `∝ dim`, `.tex:484`); per-dim η is
`~constant` (`3.52e-2 → 3.76e-2`, ~7% over the sweep — so the ratio's flatness is the
DEFECT's, apples-to-apples). The threshold is **2.0** (a modest margin above the real
1.42), mutation-proven to CATCH both pathologies: injecting a LINEAR `∝ dim/4` factor
drives the spread to **4.67** (RED); a milder `∝ sqrt(dim/4)` factor drives it to
**2.33** (RED). The upper end `dim_A=16` is bounded by the Newton–Schulz regularization
cost (~68s at a 256×256 superop, prec=128 — the dominant cost, NOT the HOPM search), so
the intermediate `dim_A=14` is dropped; six points still trace a genuine sweep at
prec=128 (lower prec destabilizes the HOPM worst-case search and INFLATES the spread, a
false trip — so prec=128 is kept). INTEGRITY guard: the STRUCTURAL witness (which
search-count noise cannot fake) is that `dim_A` genuinely SPANS a range (`dA_max >
dA_min`); a dim-trivial family collapses it to `[4,4]` ⇒ RED.

U5 THE PROJECTOR-SVD SPECTRUM (FIX C): confirms the two extraction regimes. For an
idempotent `S̃` the nonzero singular values are `≥ 1`: `σ_max = 1` (orthogonal) for
HS-self-adjoint `Φ̃` — **dephasing(4)** AND **`dep(4)+conj(dep(4))`** both measure
`1.000000` (a FINDING: the η>0 family is a convex mix of HS-self-adjoint maps, hence
itself HS-self-adjoint, hence an ORTHOGONAL projector — it is NOT oblique, contrary to
an earlier docstring; it is the *non-associativity* witness, a separate property); and
`σ_max > 1` (genuinely oblique) for non-HS-self-adjoint `Φ̃` — **`compress_idemp(4,2)`**
measures `1.732051 = √3`, exercising the strictly->1 extraction path.

### U4 — the Π_A accept-guard on a non-polar-closed A (`aic-3qq`)

The η>0 `A` is genuinely NOT polar/von-Neumann closed (a real `C*` algebra is; ours is
only ε-close), so the HOPM block step's `Π_A(polar(C))` + monotone-accept CORRECTNESS
half is finally exercised (structurally untestable while only the polar-closed η=0 `A`
existed). **Measured:** the assoc-HOPM witness is IN `A` (`proj_residual < 5.9e-16`,
certified arb) and `lo = ratio@witness = 1.55e-3`; AND `‖polar(wX) − Π_A polar(wX)‖ =
1.86e-2 > 0`, proving `A` is non-polar-closed so `Π_A` is doing real work (without it
the iterate would leave `A`).

### Mutation proofs (Rule 7 — the genuinely new teeth)

- (a) BOTTOM-`dim_A` left singular vectors in `aic_assoc_extract_range` (the wrong
  subspace) → U1 subspace match RED (`‖dΠ‖=1.0`). Restored.
- (b) `star_phi` applying `S̃` to `X` instead of `XY` (broken Choi–Effros product) →
  U1 η=0 oracle RED (unit defect `1.6`). Restored.
- (c) dim-TRIVIAL canary family (`make_eta_family` forced to `d=4`) → U3 integrity
  guard RED (`dim_A range [4,4]`; the `dim_A`-span guard is the load-bearing teeth —
  ratio-noise alone would NOT catch it). Restored.
- (FIX A) `aic_ecstar_defect_cstar_hopm` → `arb_zero(lo); return;` → U2 `c > 1e-6`
  cstar tooth RED (`cstar ~0 at t=0.01`). Estimator restored byte-identical.
- (FIX B-i) inject a LINEAR `∝ dim/4` factor into the U3 assoc/η ratio → spread `4.67
  ≥ 2.0` → RED. Restored.
- (FIX B-ii) inject a `∝ sqrt(dim/4)` factor → spread `2.33 ≥ 2.0` → RED (the milder
  pathology is caught by the factor-4 dim span + threshold 2.0). Restored.

### Files + numbers

- `include/aic_assoc.h`, `src/aic_assoc_extract.c`, `src/aic_assoc_algebra.c` (113
  LOC); `tests/test_assoc2.c` (n=139). The ecstar generalization touched
  `include/aic_ecstar.h`, `src/aic_ecstar.c` (178), `src/aic_ecstar_involution.c`,
  `src/aic_ecstar_setup.c`, `src/aic_ecstar_hopm.{c,h}` (all `≤ 200`; the Kraus path
  byte-identical — `test_ecstar` n=109 unchanged).
- U1 max axiom defect `1.69e-32`, subspace diff `≤ 8.5e-16` across 8 channels. U2
  `assoc/η ∈ [0.005,0.037]`, `cstar/η ∈ [0.75,0.82]` over t (both asserted). U3
  `assoc/η ∈ [0.033,0.047]`, real spread 1.42 over `dim_A=4..16`, threshold 2.0 (the
  universality canary; linear injection → 4.67, sqrt injection → 2.33). U4 witness
  `proj_residual < 5.9e-16`, `‖polar(wX)−Π_A polar(wX)‖ = 1.86e-2`. U5 `σ_max(S̃)`:
  dephasing `1.0`, dep+conj `1.0` (both orthogonal), compress_idemp `1.732` (oblique).

### Deferred (beaded)

- Certified cb-norm `η` (not the op-norm proxy) for the U2/U3 ratios → the
  Julia+MOSEK `eta_idempotence` / `aic_cbnorm_certify` route (`aic-m24`/`aic-ssu`).
- Radius-aware opnorm so the defect estimators consume the full certified `S̃` ball
  (not the midpoint) → same deferral as the inflated-radius wall (`aic-w4o.1`).
- The downstream channel factorization (`th_factorization`) → `aic-tff`.

## Module `funcalc` — globally-convergent non-normal `sgn` (bead aic-8hz)

The original `funcalc` `sgn` candidates (`aic_sgn_newton_schulz`,
`aic_sgn_denman_beavers`) both assert the OPERATOR-norm basin `‖X²−I‖_op < 1`
(`.tex:520-521`), the Newton–Schulz local-convergence region. For the
`assoc_ecsa` superoperator `X = 2S_Φ − 1` (NON-NORMAL: `Φ` is HP-preserving but
not HS-self-adjoint), the eigenvalues satisfy `λ²−1 = 4μ(μ−1)`, `μ ∈ spec(S_Φ)`,
so `ρ(I−X²) = 4ρ(S_Φ²−S_Φ) ≤ 4η < 1` — the SPECTRAL radius is below 1, but for a
non-normal `X` the op-norm `‖I−X²‖_op` can exceed 1, tripping the basin guard near
`η → 1/4` / large `n` even though `sgn(X)` is well-defined.

**Constructive route (Law 3).** Higham, *Functions of Matrices* (SIAM 2008),
Ch. 5 Thm 5.6: the Newton sign iteration `Y₀=X, Y_{k+1}=½(Y_k+Y_k⁻¹)` converges
GLOBALLY + quadratically to `sgn(X)` for any `X` with no purely-imaginary
eigenvalue. (Zolotarev/QDWH are Hermitian-only — inapplicable. Scaled Newton,
Kenney–Laub 1992, is a deferred Pareto candidate, `aic-8hz` notes / a follow-up
bead, not needed since our spectrum starts near `±1`.)

**Precondition, eig-free + rigorous (the Gelfand certifier).** The implemented
basin is `ρ(I−X²) < 1` — a SUFFICIENT sub-case of Higham's hypothesis (NOT the
whole: `diag(2,−3)` is Higham-valid yet aborts here), exactly the disk our
`X=2S_Φ−1` inputs occupy. `aic_funcalc_int_gelfand_rho` certifies `ρ(M)<1` for
`M=I−X²` via the Gelfand bound: `‖M^k‖_F^{1/k}` is a rigorous upper bound on `ρ(M)`
for every `k` and decreases to it, so the first `k≤k_max=32` with a certified
`‖M^k‖_F^{1/k} < 1` (whole arb ball `< 1`) proves it. `k=1` is the old Frobenius
check; non-normal `M` needs a few powers past the transient (`k=2` for the 2×2
teeth, `k=6` for the U6 superoperator). No eigensolver ⇒ sidesteps the
degenerate-eig debt `aic-w4o.1`. Fail loud if no `k` certifies (at/over `ρ=1`).

**A-posteriori certificate.** After convergence, `aic_sgn_newton_global` certifies
`‖Y²−I‖ < tol` (near an involution) AND `‖XY−YX‖ < tol` (commutes with `X`);
combined with the certified precondition and `Y₀=X` basin invariance these pin
`Y=sgn(X)` (rule out `−sgn` and wrong-inertia involutions). A straddling ball or a
singular per-step inverse aborts (Rule 4). These checks genuinely fire: at prec=53
on `[[1.05,c],[0,−0.95]]` they reject via `‖XY−YX‖` at `c=1e9` and `‖Y²−I‖` at
`c=1e12` (hostile review).

**Dispatch + prop_P.** `aic_sgn` auto-dispatches: in-basin (`‖X²−I‖_op<1` certified
by the non-aborting probe `aic_funcalc_int_in_op_basin`) → Newton–Schulz, BYTE-FOR-
BYTE the prior behavior (existing tests unchanged); else → `aic_sgn_newton_global`.
`aic_prop_P` likewise certifies the spectral `ρ(P²−P)<1/4` via Gelfand on
`M=4(P²−P)`, so `aic_theta`/`prop_P`/`aic_assoc_regularize` reach the full
non-normal `η<1/4` regime.

**Tests (`tests/test_funcalc_global.c`, n=12).** T-global-1 (teeth, exact closed-
form oracle `sgn([[1.05,20],[0,−0.95]])=[[1,20],[0,−1]]`, INDEPENDENTLY confirmed
by Mathematica `X(X²)^{-1/2}` + eigendecomposition to ~58 digits): op-basin probe
FALSE (`‖I−X²‖_op=2.005`), Gelfand `ρ_ub=0.132` at `k=2`, global Newton 7 iters,
max-diff vs oracle `< 1e-60`. T-global-2 (η=0 oracle `sgn(2P−I)=2P−I`). T-global-3
(in-basin global/DB/Newton–Schulz agree to `1e-25`; `aic_sgn` dispatch byte-equal).
T-global-4 (Gelfand rejects `ρ≥1`). The −sgn / wrong-inertia distinction is
mutation-proven: a sign flip turns T-global-1/2/3 RED. `test_assoc2` U6 is the
assoc payoff with the unitality sign witness (see the Increment-1 T9 section).

**Files.** `src/aic_funcalc_global.c` (200 LOC), `src/aic_funcalc_domain.c`
(shared `aic_funcalc_int_step_norm`), `src/aic_funcalc_sgn.c` (dispatch),
`src/aic_funcalc_proj.c` (prop_P spectral relax), `include/aic_funcalc.h`,
`src/aic_funcalc_internal.h`. Research: `docs/research/sgn_global_research.md`.

## Module `corner` (Increment 1) — compression maps `Co_{P,Q}` and corner subspaces `S_{P,Q}` (bead aic-czm)

Realizes §7 "Subspaces associated with projections" (`.tex:1052–1076`), the
substrate the th_main master loop (`cstar_build`) and §8 (`dhom`) are built on.
This is the first step of the **th_main critical path**
(`corner → projection → dhom → errreduce → cstar_build`); the §5 approximate
unitary group and §4 `unitfix` are auxiliary (the constructive Route-A projection
bypasses the §6 Lefschetz machinery, so it does not need them).

### Paper technique vs constructive route (Law 3)
The paper's `Co_{P,Q}` is an operator on the (possibly infinite-dim) Banach
algebra A, snapped to an exact idempotent by the holomorphic functional calculus
of `prop_P`. In finite dim A is a d-dim space (d = dim A) with a
Frobenius-orthonormal operator basis `{B_k}`, so an operator on A **is** a d×d
matrix. We build the left/right star-multiplication superoperators `L_P`, `R_Q`
as d×d coordinate matrices, form `2M−1 = L_P R_Q + R_Q L_P − I` (with
`M = ½(L_P R_Q + R_Q L_P)`), and call `aic_prop_P` (`θ(2M−1)`, the eig-free
Newton–Schulz route) to get the exact idempotent `Co` (d×d). The paper's
`O(δ+ε)` bound is the spec; the tests assert it.

### The construction (`.tex:1054–1066`)
`Co_{P,Q} = θ(L_P R_Q + R_Q L_P − 1)` where `(L_P)_{ij} = ⟨B_i, P⋆B_j⟩_F`,
`(R_Q)_{ij} = ⟨B_i, B_j⋆Q⟩_F`. **The star, not the plain product** (load-bearing
domain callout): `P⋆B_j = aic_ecstar_star(A,P,B_j) = Φ̃(P·B_j)` — every "product
in A" routes through the Choi–Effros star; `acb_mat_mul` appears only for the
Frobenius inner products and the d×d coordinate-matrix products `L_P R_Q`,
`R_Q L_P`. The identity `L_P R_Q + R_Q L_P − 1 = 2M − I` is exactly the `prop_P`
substitution, so the basin assert `‖M²−M‖ < 1/4` (= `‖(2M−1)²−I‖/4`) guards it
(fail-loud, Rule 4). One-sided variants `Co_{P,1} = θ(2L_P−1) = aic_prop_P(L_P)`
(R₁=id ⇒ symmetric combo collapses to `L_P`) and `Co_{1,Q} = aic_prop_P(R_Q)`;
`Co_P := Co_{P,P}`. Relations (`.tex:1061–1064`): `Co²=Co` (exact after θ),
`Co_{P,Q}(X)† = Co_{Q,P}(X†)`, `‖L_P R_Q − Co‖, ‖R_Q L_P − Co‖ ≤ O(δ+ε)`.

### `S_{P,Q} = Img Co = Ker(1−Co)` extraction (`.tex:1064`)
`dim S = round(Re Tr Co)`, **cross-checked** against the count of singular values
> 0.5 of `Co` (fail loud on mismatch, Rule 4). The basis `C_m = Σ_k u_m[k] B_k`
from the top-`dim S` **LEFT** singular vectors of `Co` (the operator→coord bridge
`x_k = ⟨B_k,X⟩_F`). `Co` is in general an **oblique** idempotent (`L_P,R_Q` need
not commute / be self-adjoint in coords), so its nonzero singular values are ≥ 1
and the **left** singular vectors span `Img Co` (range), not the right ones
(co-range = range of `Co†`) — the same subtlety as `assoc_extract`. Extraction is
the fast double path (LAPACK `zgesvd` via `aic_latd_svd`; the certified
degenerate-SVD wall `aic-w4o.1` deferred); relation defects are arb-certified.

### Cross-checks (Rule 6) and the hostile-review blocker
`tests/test_corner.c`, 33 checks at prec=256: T1 `Co²=Co` (η=0 machine-zero,
η>0 `3.3e-74`); T2 involution (relation + **closed-form** `Co_{Q,P}(X†)=QX†P`
companion + non-vacuity guard — the bare relation is structurally robust, the
historical "test that can't fail" trap); T3 `‖L_P R_Q − Co‖/η ≈ 0.007`
(measured O(δ+ε) constant, non-vacuity asserted); T4 almost-containment
(`.tex:1074`, exact-fixed + `P₁⊊P`); T5 **η=0 oracle** (identity channel,
A=M_n: `Co(X)=PXQ` machine-zero, `dim S = rank(P)·rank(Q)`, subspace match via
the gauge-free projector `Π_S`); T6 degenerate dims (0 / dim_A / 1, via the
fail-loud extract).

The hostile review (mandatory, found a blocker as on every core module) proved
T1–T6 were **blind** to two load-bearing facts: (1) star-vs-plain product —
every fixture used either the identity channel (star≡plain) or the dep+conj
family (orthogonal Φ̃, where the d×d `L_P` is machine-identical for star vs
plain), so a `build_L` that dropped the star kept all 28 checks green
byte-for-byte; (2) left-vs-right singular vectors — every `Co` was orthogonal
(σ∈{0,1}, left==right). **Fix:** T7, an **oblique** fixture `compress_idemp(4,2)`
(Φ̃ oblique, σ_max=√3, exactly idempotent → clean oracle), `P` = rank-1 projector
onto `(e₀+e₁)/√2`, `Q=|e₀⟩⟨e₀|`, giving `Co` σ_max = 2/√3 > 1. (a) independent
`build_L` oracle `L_exp[i,j] = ⟨B_i, Φ̃(P B_j)⟩_F` via `aic_assoc_superop_apply`
(machine-zero) with non-vacuity gap `‖L_act − L_plain‖ = 0.667`; (b) oblique-`Co`
fixed-point `Co(C_m)=C_m` (machine-zero). Both mutation-proven RED: star→plain
makes (a) RED at 0.667; left→right makes (b) RED at 0.408. **Surprise:** the
reviewer's suggested diagonal `P,Q` gave an *orthogonal* `Co` (the symmetric
combination washes out obliqueness) — a rotated rank-1 `P` was needed.

### Files + deferred
`include/aic_corner.h`, `src/aic_corner_compress.c` (`build_L`/`build_R`/`Co`/
one-sided), `src/aic_corner_extract.c` (coord↔operator apply + SVD extraction),
`src/aic_corner_internal.h` (shared `aic_corner_frob_ip`), `tests/test_corner.c`
(779 LOC — joins the oversized-test split list `aic-w4o.4`). Certified `dim S`
rank defers to `aic-w4o.1` (double-path SVD now, like `idemp`/`assoc`). `.tex`
typo escalation (Increment 2, lem_alpha): `tex:1109` `β_{jk}=Co_{P_j,Q_j}` should
read `Co_{P_j,Q_k}` (proof `tex:1114` `δ_{jl}δ_{km}` requires `Q_k`).

## Module `corner` (Increment 2a) — compressed product + `lem_alpha` block bijection (bead aic-czm)

Adds the §7 product machinery on top of Increment 1.

### Compressed product (`.tex:1077–1082`)
`X·Y = Co_{P,R}(X⋆Y)` for `X∈S_{P,Q}`, `Y∈S_{Q,R}` — the star `X⋆Y=Φ̃(XY)`
(`aic_ecstar_star`) then the d×d `Co_{P,R}` via `aic_corner_apply`
(`aic_corner_cdot`). `(S_P,·)` is an O(δ+ε)-C* algebra with unit
`P̃=Co_P(P)` (`aic_corner_Ptilde`). **Finding (load-bearing for the test design):
`Co_{P,R}` is a NO-OP on a valid `X⋆Y` at η=0** — the star already lands
(approximately) in `S_{P,R}`, so `Co_{P,R}(X⋆Y)=X⋆Y` exactly at η=0 (the paper's
"close to the ambient product", `.tex:1081`); the compression step is visible
only at η>0. So the two cdot teeth are distinct: T8(a) catches a dropped *star*
(non-vacuity gap 0.0556 on the oblique `compress_idemp(4,2)` fixture, mutation-RED),
T8(b)'s η>0 arm catches a dropped *compression* (`X·Y==X⋆Y` at η>0, RED).
Measured `‖X·Y−X⋆Y‖/η ≈ 0.006`.

### `lem_alpha` block bijection (`.tex:1086–1119`)
`α=Σ α_{jk}`, `α_{jk}=Co_{P,Q}|_{S_{Pj,Qk}} : ⊕_{jk} S_{Pj,Qk}→S_{P,Q}`;
`β_{jk}=Co_{Pj,Qk}`; `βα=1+γ`, the contraction `‖γ‖<1` makes `βα` invertible and
`α⁻¹=(βα)⁻¹β` (certified `acb_mat_solve`). Built as explicit matrices in the
extracted corner coordinates: `α[l,⟨jk⟩m]=⟨D_l,Co_{P,Q}(C^{jk}_m)⟩_F`,
`β[⟨jk⟩m,l]=⟨C^{jk}_m,Co_{Pj,Qk}(D_l)⟩_F` (block index `b=j·q+k`, running column
offset). Used only for `p,q≤2` (`.tex:1084`). **`.tex` typo (documented, Law 1):**
`tex:1109` prints `β_{jk}=Co_{P_j,Q_j}`; the codomain `S_{Pj,Qk}` and the proof
`tex:1114` (`δ_{jl}δ_{km}`) require `Co_{P_j,Q_k}` — implemented `Q_k`, the `Q_j`
build is mutation-RED (`‖γ‖` jumps 0→1). The dim-count bijection oracle
`N=Σ dim S_{Pj,Qk}=dim S_{P,Q}` is a fail-loud assert (`.tex:1124`).

### Two soundness fixes from the hostile review (no blockers)
(1) **The `α⁻¹·α≈I` round-trip is a tautology** (`α⁻¹=(βα)⁻¹β` ⇒ `α⁻¹α=I` for
any α,β — a β→1.5β mutation kept it green). Replaced by an **independent η=0 entry
oracle**: α at η=0 is the exact change-of-basis isometry, so `‖α†α−I‖<tol`
(orthonormal columns) AND `Σ_l α[l,·]D_l == C^{jk}_m` (reconstruction, no `α⁻¹`
path) — mutation-RED under a column-scale AND a column-swap (the swap is invisible
to `α†α` but caught by reconstruction). (2) **The `‖γ‖<1` gate is now a CERTIFIED
upper bound** (`aic_corner_gamma_opnorm_ub`): `‖M‖_op ≤ ‖mid(M)‖_op + ‖rad(M)‖_F`,
the radius term via rigorous `mag_t` arithmetic — the decision uses the upper end
< 1 (honors the arb ladder), sidestepping the `aic_mat_opnorm` Gram false-fail on
near-zero off-diagonals (bead aic-2yo). Measured radius contribution ~5e-72 at
the tested η (negligible); self-mutation-proven (mid=0.999/rad=0.01 → 1.009 fails).

### Files (Increment 2a)
`src/aic_corner_product.c` (cdot + Ptilde, 67 LOC), `src/aic_corner_alpha.c`
(α assembly + γ + solve, 183 LOC), `src/aic_corner_alpha_build.c`
(build_blocks/free_blocks/alpha_dims + the certified-norm helper, 160 LOC),
header + tests extended (73 checks total). `lem_alpha` hypotheses are not
separately asserted — they are sufficient conditions for `‖γ‖<1`, which is the
load-bearing guard (`tex:1114`); the `prop_P` basin assert in `aic_corner_Co`
gates gross overlap upstream.

## Module `corner` (Increment 2b) — the 1d-projection Hilbert/Ha cluster (bead aic-czm)

Completes §7: the one-dimensional-`Q` machinery (`.tex:1123–1189`). This is what
the downstream `lem_extension` (in `cstar_build`) leans on — `Ha^Q_{P,P}` as an
O(δ+ε)-homomorphism is the linearization that extends an isomorphism by one
matrix dimension.

### `lem_PQ_Hilb` inner product (`.tex:1123–1132`, `aic_corner_ip_1d`)
For a **one-dimensional** δ-projection `Q` (assert `dim S_Q=1`, fail-loud),
`S_{P,Q}` is a Hilbert space with `⟨Y,X⟩` defined by `Co_Q(Y†⋆X)=⟨Y,X⟩Q̃`,
`Q̃=Co_Q(Q)`. Since `S_Q=ℂQ̃` is 1-dim, the scalar is `⟨Q̃, Co_Q(Y†⋆X)⟩_F /
⟨Q̃,Q̃⟩_F`. Conjugate-linear in `Y`, linear in `X`; `|⟨X,X⟩−‖X‖²|≤O(δ+ε)‖X‖²`
(measured c≈0.007). The star `Y†⋆X` is load-bearing (oblique non-vacuity gap
0.667, mutation-RED).

### The Ha-map (`.tex:1146–1160`, `aic_corner_ha`) — the intricate one
`Ha^Q_{P,R}(Z): S_{R,Q}→S_{P,Q}`, defined by the implicit equation
`(Y†·Z)·X + Y†·(Z·X) = 2⟨Y, Ha(Z)(X)⟩Q̃` for all `Y∈S_{P,Q}`. Solved as a **Gram
system**: with `{C_l}` the basis of `S_{P,Q}`, `G_{lm}=⟨C_l,C_m⟩` the
**lem_PQ_Hilb inner product** (NOT Frobenius — `G=I+O(δ+ε)`, assert `‖G−I‖<1`),
and `b_l=½[scalar((C_l†·Z)·X)+scalar(C_l†·(Z·X))]`, then `Ha(Z)(X)=Σ_l(G⁻¹b)_l C_l`
(certified `acb_mat_solve`). **The Co-index bookkeeping is the highest
convention-risk in the module** — each compressed product `·` targets a specific
corner: `Y†·Z`→Co_{Q,R}, `(Y†·Z)·X`→Co_Q, `Z·X`→Co_{P,Q}, `Y†·(Z·X)`→Co_Q.
Properties: `Ha_dag` exact (`Ha^Q_{R,P}(Z†)=Ha^Q_{P,R}(Z)†`), `Ha(Z)(X)≈Z·X`,
and `Ha^Q_{P,P}` an O(δ+ε)-homomorphism (`Ha(Z·W)≈Ha(Z)·Ha(W)`, measured c≈0.001).

### `lem_PQR` / `lem_1d_proj` / equivalence (`.tex:1162–1187`)
For 1d `Q`: `|‖X·Y‖−‖X‖‖Y‖|≤O(δ+ε)‖X‖‖Y‖` (`aic_corner_pqr_defect`, norm test,
exact at η=0). For 1d `P` AND `Q`: `dim S_{P,Q}≤1` (`lem_1d_proj`). Equivalence
`P~Q ⇔ dim S_{P,Q}=1` (`aic_corner_equiv_1d`), reflexive/symmetric/transitive
(transitivity via `lem_PQR`).

### The hostile-review BLOCKER (the consequential one) + fixes
The review found `aic_corner_ha` **SIGABRTs** on a valid in-hypothesis
`dim S_{P,Q}=2, η>0` corner — exactly the regime `lem_extension` needs — because
its `‖G−I‖<1` assert called the certified `aic_mat_opnorm`, whose Gram path hits
the **aic-2yo Hermiticity false-fail** on the near-zero off-diagonal entries (arb
radius exceeds the absolute floor). The Ha math was correct; the *guard*
false-aborted. **Fix:** use the certified mid+radius `aic_corner_gamma_opnorm_ub`
(the alpha module's same-issue workaround) — applied to `aic_corner_ha` AND, a
**second instance the fix caught**, `aic_corner_pqr_defect` (the cdot output
carries a ~1e-70 accumulated radius that trips the same check despite `‖X·Y‖≈1`).
The review also constructed the missing regression cell — `make_mixconj` (convex
mix of `compress_idemp(5,3)` with its unitary conjugate, η≈0.013, genuine
δ-projections, `dim S_{P,Q}=2`, off-diagonal Gram `‖G−I‖≈1.3e-5`): the routine
now runs (before SIGABRT, after exit 0), and this fixture gives the
**lem_PQ_Hilb-vs-Frobenius Gram distinction its first teeth** (replacing `G` with
the Frobenius `I` moves the Ha output 2.3e-6→1.18e-5, mutation-RED). Two further
can't-fail test-gaps closed: the ip_1d **conjugation order** (all fixtures used
real symmetric corner elements → blind; a complex `⟨i·C0,C0⟩=−i` test now pins
conjugate-linearity-in-Y, swap RED at 2.0) and the **term2 Co-index** (blind
because every fixture had `range(Q)⊆range(P)`; a disjoint-support fixture
`P=diag(0,0,1,1),Q=|e0⟩⟨e0|` makes `Co_Q→Co_{P,Q}` RED at 0.5).

### Files (Increment 2b)
`src/aic_corner_hilbert.c` (ip_1d + pqr_defect + dim_S + equiv_1d, 203 LOC),
`src/aic_corner_ha.c` (the Ha-map, 221 LOC — marginally over the ~200 soft limit,
in line with `extract.c` 205 precedent; literate docstring, single cohesive
routine; noted on `aic-w4o.4`). `test_corner.c` now 134 checks total. The corner
module (§7) is COMPLETE; next on the th_main critical path is `projection` (§6,
bead aic-mqf, Route-A Hermitian-eigensolve nontrivial-projection finder).

## Module `projection` — constructive nontrivial O(ε)-projection (bead aic-mqf, the constructivization crux)

Realizes `lem_nontriv_projection` (`.tex:931`): any ε-C* algebra A with
`1 < dim A < ∞` has a nontrivial O(ε)-projection. The paper proves *existence*
non-constructively (Lefschetz–Hopf on the approximate-unitary quotient manifold,
`.tex:944–969`); we replace that with a direct finite-dim construction (Law 3).
The §6 entry point of the th_main master loop.

### The route (decided by 3 independent research legs, Rule 9)
**Reduction (`.tex:935–939`):** a nontrivial projection is `P=½(I+X)` for a
Hermitian near-involution `X∈A` (`X⋆X≈I`, `X≠±I`). **The key correction the
independent legs converged on:** by **rem_X2 (`.tex:628`)** the functional
calculus does *not* generalize to the ε-C* algebra — so do the sign calculus
**AMBIENT** (in M_n, an exact C* algebra), never inside A. (Leg 1 first proposed
in-A star-Newton-Schulz; legs 2+3 + the web survey rejected it — no literature
supports sign-iteration under ε-associativity, and ambient Hermitian spectra are
*stable*, the fragility being non-normal-only.)

**Route A-ambient (`src/aic_projection.c` + `src/aic_projection_find.c`):**
1. **Pick H** (`aic_projection_pick_H`): scan the Hermitianized basis
   `H_k=½(B_k+B_k†)`, choose the one with the **largest interior eigenvalue gap**
   (tie-break by spectral spread; spread via certified `aic_mat_herm_max_eig`).
   Fail-loud if all `H_k` near-scalar (a stop condition: dim A>1 but no non-scalar
   Hermitian element).
2. **Gap** (`aic_projection_gap`): eigenVALUES of H via LAPACK zheev (double path —
   Hermitian ⟹ stable spectrum; only the values, not eigenvectors, so no
   `aic-w4o.1` dependency); largest interior gap `[λ_m,λ_{m+1}]`, threshold
   `t=½(λ_m+λ_{m+1})`. Fail-loud only if **no positive interior gap** exists (the
   genuine aic-3qv degenerate-spectrum stop condition).
3. **Ambient projector** (`ambient_projector`): `Y=s(H−tI)`, `s=1/max(t−λ_min,λ_max−t)`;
   assert the eig-free sgn basin `‖Y²−I‖<1` (certified via
   `aic_corner_gamma_opnorm_ub` to dodge the aic-qgs Gram false-fail; holds for any
   `g>0` since `‖Y²−I‖=1−(s·g/2)²`); `X=aic_sgn(Y)` (AMBIENT, eig-free);
   `P_amb=½(I+X)` — an exact ambient idempotent onto `λ>t`.
4. **Project into A** (`project_into_A`) — **the load-bearing discovery
   (corrected the research spec):** `A=ImgΦ̃` is an *oblique* image (Φ̃ is HP but
   not HS-self-adjoint), so the Frobenius-orthogonal projector `Π_A` does NOT
   respect the star structure — it leaves `‖P⋆P−P‖=O(1)` (~0.5, constant in η).
   The correct projection is **Φ̃ itself**, available through the public star API
   as `P = P_amb⋆I = Φ̃(P_amb)`, giving the O(η) defect the lemma promises.
5. **Certify (arb):** `δ=‖P⋆P−P‖_op` (star, via `aic_corner_gamma_opnorm_ub`);
   nontriviality `‖P‖_op,‖I−P‖_op` (fail-loud if either ≤0.3); membership
   `aic_ecstar_proj_residual(P)`.

rem_X2-safe (sgn only on ambient Y), eig-free for the projector (eigenVALUES only,
for gap-finding), no §5 unitary group (no aic-q2x dependency). **Audition (Law 4):**
Route A-ambient primary; Route B (σ-fixed-point, needs aic-q2x) and Route C
(optimize `‖P⋆P−P‖`, needs warm start) are Pareto alternatives, deferred.

### The Ω(1)-gap question (aic-3qv) — honestly bounded
The defect is `δ=O(ε+ε/g)`, = O(ε) iff the gap `g=Ω(1)`. The construction
certifies the gap **per-instance** and fail-loud-aborts only on a genuinely
degenerate spectrum (no positive interior gap). The universal a-priori guarantee
("dim A>1 ⟹ some H has an Ω(1)-relative gap") is what the paper needs Lefschetz
for and is NOT proven constructively — escalation aic-3qv. Empirically the
constant is dimension-independent (canary below).

### Cross-checks (Rule 6/7) — `tests/test_projection.c`, 58 checks
- T1 **η=0 oracle**: M_2/M_3 (identity), block-d4 → `δ` machine-zero, `‖P‖=‖I−P‖=1`,
  P∈A. (Blind to the oblique fix, since Φ̃=id at η=0 — only T2 catches it.)
- T2 **η>0** (mixconj): `δ/η` small and flat; non-vacuity rewritten to be
  fixture-robust — (A) the most-oblique `P_amb` has its star defect reduced **459×**
  by Φ̃, (B) `P=Φ̃(P_amb)` for the chosen H (mutation: skip Φ̃ → RED at 3.7e-3).
- T3 **universality canary** (aic-dbo.3, the highest-value test, tex:484): η>0 sweep
  over 5 dims (dim_A∈{4,9,16}), `(d,C=δ/η)` = (4,5.6e-2),(5,3.6e-4),…,(8,3.0e-3);
  **slope assert** `C(d_max)/C(d_min)<1.5` (mutation: inject `C=0.05·d` → slope 2.0
  RED). C does not grow with d (it shrinks) — the universal-constant claim holds.
- T4 **fail-loud** (message-grepped, fork+stderr-capture): dim-1 abort
  (`"dim A = 1 <= 1"`) and the **aic-3qv no-interior-gap abort**
  (`"NO positive interior spectral gap"`), each mutation-proven RED on guard deletion.
- T5 double-vs-arb@53 agreement; T6 asymmetric-spectrum fixture ({−5,0,1,2},{−2,−1,3})
  asserts `m==count(λ≥t)` and `rank(P_amb)==m` (mutation: off-by-one in m → RED).

Hostile review (no blockers): core verified rem_X2-safe, m/threshold correct on
asymmetric spectra, dimension-independent to d=9, Φ̃-vs-Frobenius protected by T2.
Four SHOULD-FIX applied: the `g_min` floor was `0.05·spread` (spurious abort on a
well-conditioned uniform spectrum at n≈22, in-range) → now a positive-gap floor
`1e-9·max(1,spread)` + gap-ranked H-pick (uniform spectrum succeeds to n=24); T4
cannot-fail (dim-1 aborts via any guard) → message-grep + the aic-3qv path tested;
the canary print-only → slope assert; + the asymmetric/oblique coverage.

### Files
`include/aic_projection.h`, `src/aic_projection.c` (orchestrator + ambient_projector
+ project_into_A + star_defect + nontriviality, 207 LOC), `src/aic_projection_find.c`
(H-pick + gap, 195 LOC), `src/aic_projection_internal.h`, `tests/test_projection.c`.
Certified gap ENCLOSURE (Rump) defers to `aic-w4o.1`; gap-finding uses double-path
zheev now, the defect is arb-certified. Next th_main step: `dhom` (§8, bead aic-c1n).

## Module `dhom` — δ-homomorphism approximation (bead aic-c1n, §8)

Realizes `prop_delta_hominc` (`.tex:1194`) + `lem_approx` (`.tex:1224`): any
δ-homomorphism `v` from a genuine finite-dim C* algebra B into an ε-C* algebra A
is `O(δ)`-close to an `O(ε)`-homomorphism `ṽ`. The §8 tool the master loop's
`lem_extension`/`cor_improvement` lean on. **Fully constructive, eig-free**
(pure products + norms).

### Data model
B = ⊕_l M_{d_l} (block dims + matrix-unit basis `E^{(l)}_{ab}`, block matmul,
block adjoint, unit ⊕I_{d_l}). `v: B→A` stored as `dim_B` n×n matrices
`vE[i]=v(E_i)`; `v(X)=Σ_i x_i vE[i]`. **The multiplicativity defect**
`G_v(X,Y)=v(XY)−v(X)⋆v(Y)` uses *two* products: `XY` is B's block matmul, but
`v(X)⋆v(Y)` is A's **Choi–Effros star** (`aic_ecstar_star`, load-bearing — never
plain matmul; the η=0 identity-channel oracle is blind to this, only the η>0 T4
catches it, the corner lesson).

### The generalized-Pauli diagonal (dimension-independence crux, tex:484)
`D=Σ_jk d⁻² S_jk†⊗S_jk`, `S_jk=X^jZ^k` (shift·clock Heisenberg–Weyl). The
**exact** diagonal (central `XD=DX`, `π(D)=I`, `‖D‖_proj=1`) is what makes the
th_main constant dimension-independent (a naive ε-averaged diagonal carries O(n)
error). **A real `.tex` finding (see `paper/FINDINGS.md`):** the paper's `tex:1254`
direct-sum formula (Cartesian product of joint block-diagonal unitaries) is
**non-central for the finite Pauli 1-design** (`Σ_jk S_jk≠0` unlike Haar
`∫U dU=0`, leaving spurious cross-sector terms; measured `‖XD−DX‖=0.54` for
M₂⊕M₂). The correct finite central diagonal is the cross-term-free **embedded
sum** `D=Σ_l D_l` (per-block Pauli embedded as a partial isometry in sector (l,l)),
`nterms=Σ_l d_l²` — this is the Haar diagonal `∫dU U†⊗U`. Consequence: `‖D‖_proj=m`
(block count), not 1, so the `w'` bound is `O(mδ)`; the canary (below) checks this
m-dependence does not surface in the constant.

### lem_approx Newton iteration
`w'(X)=Σ_t v⁽ˢ⁾(A_t)⋆g(B_t,X)` (star; uses the *current* v⁽ˢ⁾), `w''(X)=w'(X†)†`,
`v⁽ˢ⁺¹⁾=v⁽ˢ⁾+½(w'+w'')` — defect `O(δ)→O(δ²+ε)` per step (the F_w cancellation
`tex:1282-1304` relies on `XD=DX`), iterate to `≤O(ε)`, fail-loud on the cap.
`prop_delta_hominc` = 3 predicates (‖v‖ upper, lower-bound surrogate,
unit-preservation). Certified op-norms via `aic_corner_gamma_opnorm_ub` (aic-qgs).

### Cross-checks — `tests/test_dhom.c`, 89 checks
T1 Pauli identities (unitarity, π(D)=I, **centrality incl. the direct-sum case**
M₂⊕M₂/M₂⊕M₃, ‖D‖_proj; mutation: clock-phase k=0 → centrality RED 0.48; the
old joint Cartesian build → direct-sum centrality RED 0.71). T2 η=0 oracle
(`v(X)=UXU†` exact hom → G_v=0, 0 steps). T3 Johnson (δ>0,ε=0 → exact hom in 4
steps, ‖ṽ−v‖=1.02δ). T4 ε>0 (mixconj, C=defect/ε=2.2; **star-mutation RED**:
defect floors at 1.11, both star call sites independently caught). T5
universality canary (block-dim sweep M₂/M₃/M₄ ratio 1.31; **block-count sweep**
m=1,2,3: C={2.71,1.24,1.56} — does NOT grow with m; Pauli-mutation → C grows RED).
T6 prop bounds. Hostile review found the direct-sum non-centrality BLOCKER (T1
never tested multi-block diagonals — a "test that cannot fail"); fixed + the
direct-sum centrality test added.

### Files
`include/aic_dhom.h`, `src/aic_dhom_B.c` (B + matrix units), `src/aic_dhom_diag.c`
(Pauli diagonal + the embedded-sum fix + the `.tex:1254` discrepancy box),
`src/aic_dhom_defect.c` (v + G_v + prop_delta_hominc), `src/aic_dhom_approx.c`
(Newton), `src/aic_dhom_internal.h`. No eig anywhere. `c_0` (cor_improvement) is
unstated in the paper → the next module `errreduce` (aic-t81) / research aic-1bc.
Next th_main step: `errreduce` (§8/9, aic-t81) → `cstar_build` (§9 master loop,
aic-097).

## Module `errreduce` — error reduction / `cor_improvement` (bead aic-t81, §8/9)

Realizes `cor_improvement` (`.tex:1317`): a δ-inclusion of a genuine C* algebra B
into an ε-C* algebra A is upgraded to a `c₀ε`-inclusion (bijective preserved).
The error-reset the th_main master loop applies after every merge/extend. A
**thin orchestration of `dhom`**: `ṽ=lem_approx(v)` (an O(ε)-homomorphism with
`‖ṽ−v‖≤O(δ)`), then `prop_delta_hominc` applied to `ṽ` (whose own defect is O(ε))
certifies `ṽ` is an O(ε)-inclusion. The MEASURED `c₀=max-defect/ε` is returned;
the analytic `c₀` defers to aic-1bc (FINDINGS §D2).

### The c₀ universality canary (the headline th_main correctness test, tex:484)
`errreduce` is where `c₀` is measured, so its canary (T4) is *the* check that the
th_main constant is **dimension-independent**. Block-dim sweep B=M₂/M₃/M₄/M₅
(dim_B 4/9/16/25): `c₀=2.71/2.22/2.07/1.96` — **decreases** with dim; ratio
`c₀_max/c₀_min=1.38`, asserted `<1.6` over the 6.25× dim range (a √dim-growth
injection trips it at 1.807 — mutation-proven). Block-count sweep m=1,2,3 flat
too. So the explicit-Pauli-diagonal route (dhom) delivers the universal constant
the paper claims; no escalation.

### The hostile-review BLOCKER + fix (the soundness hole)
The δ-inclusion guards (input "is v a δ-inclusion?" + the output `lower_gap`
certification + `is_bijective` injectivity) originally used the **basis-sweep**
`min_i‖v(E_i)‖` — **blind to a map that collapses a non-basis direction**
(`v(E_0)=diag(1,0), v(E_1)=`rank-1 at angle 0.1: basis-sweep `a=1.0` passes, but
`‖v(E_0−E_1)‖=0.10` violates the hypothesis → silently certified `c₀=0`). Fix:
`aic_dhom_v_sigma_min` (`src/aic_dhom_sigmin.c`) — assemble the coordinate matrix
`M[k,i]=⟨B_k,vE[i]⟩_F` (both bases Frobenius-orthonormal) so
`σ_min(M)=inf_{X≠0}‖v(X)‖_F/‖X‖_F` (true unit-ball lower bound in the coordinate
norm, via double-path SVD `aic_latd_singular_values`); it **sees all combinations**
and catches collapses (the collapsing-v fixture now aborts, σ_min=0.0998<0.5;
genuine inclusions pass at σ_min≈1). The exact *operator-norm* inclusion inf
(HOPM, like ecstar) and the certified σ_min enclosure defer to aic-w4o.1/aic-w4o.2;
the coordinate-norm σ_min is a sound collapse-detector now.

### Files
`include/aic_errreduce.h`, `src/aic_errreduce.c` (the orchestration + the 4
inclusion-defects + the O(ε)-inclusion certification + measured c₀ + bijectivity),
`src/aic_dhom_sigmin.c` (the σ_min collapse-detector, shared with dhom),
`tests/test_errreduce.c` (51 checks: T1 η=0 oracle, T2 Johnson, T3 ε>0, T4 the
c₀ canary, T5 bijectivity, T6 the collapse-abort). `C0_CERT=10` (output-inclusion
ceiling), `EPS_FLOOR=4` (lem_approx termination floor — the defect cannot beat
~c₀ε). No eig in errreduce; σ_min uses double-path SVD.

**th_main status:** `corner`✓ `projection`✓ `dhom`✓ `errreduce`✓ — only the §9
master loop **`cstar_build`** (aic-097) remains; it assembles all four into the
proof of `th_main`. See HANDOFF.md "RESUME HERE".
