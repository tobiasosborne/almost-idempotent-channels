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
- `||Phi^2-Phi||_cb` diamond-norm SDP: bead aic-d24 (out of scope).
