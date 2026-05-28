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
