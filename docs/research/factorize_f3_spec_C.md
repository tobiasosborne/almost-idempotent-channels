# F3 Design Spec — Designer C (independent re-derivation + literature conventions)

> Competing designer C for bead `aic-tff`, Increment F3: the UCP decode map
> Upsilon via lem_RC (Lemma 12.3, `.tex:2840-2899`).
> Read-only on `.c`/`.h`; this doc is a clean independent derivation from
> `.tex` and standard QI references, serving as a cross-check on conventions.
>
> Standard references used throughout:
> - Watrous, *Theory of Quantum Information* (2018) — "TQI"
> - Paulsen, *Completely Bounded Maps and Operator Algebras* (2002)
> - Hayden-Leung-Winter / Dankert et al. unitary 1-design / twirl literature
>
> All `.tex` line citations are to `paper/src/approximate_algebras.tex`.

---

## A. Construction step-by-step

### A.0. Prerequisite objects from prior increments

F3 requires the following, all produced by F1+F2 and the upstream modules:

- `phi`: the original UCP map Phi: B(H) -> B(H), dim_H = N, Kraus operators
  `K_a: H -> H` (`aic_ucp_kraus`, `r` operators, each `N x N`). Available in
  `F->phi`.
- `V`: the Stinespring isometry of Phi, `V: H -> H tensor F`, assembled via
  `aic_ucp_kraus_to_stinespring`. Shape: `(N*r) x N` with left-major indexing
  `V[a*N + i, j] = K_a[i,j]` (the `aic_ucp.h` convention). Unitality:
  `V^dag V = 1_H`.
- `Delta` (F2's UCP encode): B -> B(H), available via `aic_factorize_delta`.
  Its Kraus rep / Choi is the input to lem_RC.
- Per-block 1-design `{(p_{js}, U_{js})}`: the Heisenberg-Weyl (Pauli) design
  on each `B(L_j)` (dimension `d_j x d_j`), `d_j^2` terms, weights
  `p_{js} = d_j^{-2}`. Built from `aic_dhom_pauli`. Same design used in F2.

The genuine C* algebra is `B = (+)_{j=1}^m B(L_j)` (block dims `d_j`, block
count `m`), accessed via `F->v->B`.

---

### A.1. Choi representation of Delta and the Stinespring V

`.tex:2831-2837` (eq. `Choi_Delta`):

```
Delta(X) = sum_{j=1}^m  W_j^dag (X_j tensor 1_{E_j}) W_j       (Choi_Delta)

where  W_j : H -> L_j tensor E_j,   sum_j W_j^dag W_j = 1_H.
```

Here X = (X_1, ..., X_m) in B = (+)_j B(L_j); X_j in B(L_j); 1_{E_j} is the
identity on an ancilla space E_j. The W_j are partial Kraus operators of Delta,
and the full Stinespring isometry of Delta is the column-stack of the W_j.

**Matrix shape and extraction.** Delta has Kraus operators
`{G_{j,ab}}_{j,a,b}` where `(a,b)` ranges over `B(L_j)` matrix-unit pairs; the
W_j are assembled as sub-blocks of Delta's Stinespring isometry `V_Delta`:
`V_Delta: H -> (+)_j (L_j tensor E_j)`. In the left-major (`aic_ucp.h`)
convention, `V_Delta` is `(sum_j d_j * r_j) x N` and the j-th block of rows
corresponds to `W_j: H -> L_j tensor E_j`.

**Practical extraction.** Build Delta's Choi matrix `C_Delta = sum_{a,b} E_{ab}
tensor Delta(E_{ab}^{(j)})` (Convention A), run `aic_ucp_choi_to_kraus_latd` to
recover the Kraus operators, then group by which block j of B each Kraus operator
maps into; reshape each group as `W_j` (`d_j * dim_E_j x N`). Alternatively, a
direct column-extraction from `V_Delta` suffices. The exact splitting is a
gauge choice; what matters is that `sum_j W_j^dag W_j = 1_H` and
`W_j W_j^dag` lives in `B(L_j tensor E_j)`.

---

### A.2. R_j and the lem_RC centrality argument

`.tex:2840-2857` (lem_RC, eq. `R_def`):

```
R_j = sum_s p_{js} (U_{js}^dag tensor 1_{E_j}) W_j W_j^dag (U_{js} tensor 1_{E_j})

     in B(L_j tensor E_j).
```

**Shape.** `R_j` is `(d_j * dim_E_j) x (d_j * dim_E_j)`, Hermitian PSD
(a weighted average of PSD operators).

**Why R_j = 1_{L_j} tensor C_j (the Schur / commutant argument, lem_RC(i)).**

The 1-design centrality property `diag_j2` (`.tex:2776`) is:

```
sum_s p_{js} X U_{js}^dag tensor U_{js}  =  sum_s p_{js} U_{js}^dag tensor U_{js} X
    for all X in B(L_j).
```

This is exactly the statement that the operator `D_j = sum_s p_{js} U_{js}^dag
tensor U_{js}` commutes with every `X tensor 1` in `B(L_j) tensor 1`. In
Schur / representation-theory terms: the twirl

```
T_j(A) = sum_s p_{js} (U_{js}^dag tensor 1) A (U_{js} tensor 1)
```

is a projection onto the commutant of `{X tensor 1 : X in B(L_j)}` inside
`B(L_j tensor E_j)`. By Schur's lemma (finite-dimensional version), since
`B(L_j)` acts irreducibly on `L_j`, its commutant in `B(L_j tensor E_j)` is

```
commutant of {X tensor 1_{E_j} : X in B(L_j)}  =  1_{L_j} tensor B(E_j).
```

This is a standard result in quantum information: see TQI, Theorem 2.22
(Schur's lemma in the finite-dimensional setting); also Paulsen Chapter 1 for
the algebraic version. The twirl T_j onto this commutant gives
`T_j(W_j W_j^dag) = 1_{L_j} tensor C_j` for some `C_j in B(E_j)`.

**Proof that R_j is the twirl of W_j W_j^dag.** Directly from the `R_def`
formula: R_j = sum_s p_{js} (U_{js}^dag tensor 1)(W_j W_j^dag)(U_{js} tensor 1).
This is the twirl T_j applied to the PSD operator W_j W_j^dag. So R_j = T_j(W_j
W_j^dag) = 1_{L_j} tensor C_j. QED (lem_RC(i)).

**Extraction of C_j in finite dim.** Given R_j = 1_{L_j} tensor C_j, the operator
C_j in B(E_j) is recovered by a partial trace over the L_j factor:

```
C_j = (1 / d_{L_j}) Tr_{L_j}(R_j)
```

where d_{L_j} = dim L_j = d_j. This follows from the identity
`Tr_{L_j}(1_{L_j} tensor C_j) = d_{L_j} C_j`.

**Partial-trace direction — the critical convention question; full analysis in
section C below.**

**lem_RC(ii): 1 - O(eta) <= ||C_j|| <= 1.** Upper: ||C_j|| = ||R_j|| <= 1
since ||W_j|| <= 1 (W_j is part of a partial isometry with sum_j W_j^dag W_j =
1_H) and U_{js} are unitaries. Lower: the paper's proof at `.tex:2849-2857`
uses PhiDelta2 (`.tex:2810`) and Delta_norm (`.tex:2806`) — the argument is
reproduced in shard H's lem_RC proof note and needs no re-derivation here.

---

### A.3. xi_j — leading right singular vector of C_j

`.tex:2859`:

```
xi_j in E_j  a unit vector such that  | ||C_j xi_j|| - 1 | <= O(eta).
```

Since ||C_j|| >= 1 - O(eta) (lem_RC(ii)) and ||C_j|| = sigma_max(C_j), the
unit vector that maximizes ||C_j v|| is the leading right singular vector of
C_j. So:

```
xi_j = v_1(C_j)     (leading right singular vector of C_j)
||C_j xi_j|| = sigma_max(C_j) >= 1 - O(eta)
```

**Algorithm:** run SVD of C_j (dimension `dim_E_j x dim_E_j`; uncertified LAPACK
double path via `aic_mat_singular_values` or the full `LAPACKE_zgesvd`); take
the leading right singular vector. The `aic_mat_singular_values` currently
requires DISTINCT singular values (see `aic_mat.h:120`); for the singular vector
itself, use LAPACK's `zgesvd` (double path) since the degenerate case (multiple
singular values at sigma_max) is allowable — any xi_j in the leading singular
subspace achieves ||C_j xi_j|| = sigma_max(C_j). Use the double path; certify the
lower bound `sigma_max(C_j) >= 1 - O(eta)` via the arb path.

---

### A.4. L_j — the Kraus-like operator

`.tex:2860-2863`:

```
L_j : L_j -> H tensor F,

L_j = sum_s p_{js} (Delta(U_{js}^dag) tensor 1_F) V W_j^dag (U_{js} tensor xi_j)
```

where V: H -> H tensor F is the Stinespring isometry of Phi.

**Matrix shape analysis.**

- `U_{js} tensor xi_j: L_j -> L_j tensor E_j` is the map `u |-> U_{js} u tensor xi_j`,
  a `(d_j * dim_E_j) x d_j` matrix with left-major indexing
  `(U_{js} tensor xi_j)[a*dim_E_j + b, c] = U_{js}[a,c] * xi_j[b]`.
- `W_j^dag: L_j tensor E_j -> H` is `dim_E_j*d_j x N` (transposed from `W_j`),
  so `W_j^dag (U_{js} tensor xi_j): L_j -> H` is `N x d_j`.
- `V: H -> H tensor F` is `(N*r) x N`.
  `V W_j^dag (U_{js} tensor xi_j): L_j -> H tensor F` is `(N*r) x d_j`.
- `Delta(U_{js}^dag): B(H) -> B(H)` is an N x N matrix (the image of U_{js}^dag
  under the UCP map Delta); here we apply it as a LEFT multiplier tensored with
  1_F via `(Delta(U_{js}^dag) tensor 1_F): H tensor F -> H tensor F`,
  a `(N*r) x (N*r)` matrix.
  `(Delta(U_{js}^dag) tensor 1_F) V W_j^dag (U_{js} tensor xi_j): L_j -> H tensor F`
  is `(N*r) x d_j`.

So `L_j` is `(N*r) x d_j`. It maps L_j (dimension d_j) into H tensor F
(dimension N*r), as stated in `.tex:2861`.

**Convention note.** `Delta(U_{js}^dag)` is applied as an N x N matrix acting
on the H factor of H tensor F, hence the tensor extension
`Delta(U_{js}^dag) tensor 1_F = kron(Delta(U_{js}^dag), 1_F)` in the
LEFT-major convention: rows and columns indexed as `(i*r + a, j*r + b)` for
i,j in [0,N) and a,b in [0,r). This matches `aic_mat_kronecker` with A = Delta(U)
and B = 1_F.

---

### A.5. Upsilon'_j — the manifestly CP component map

`.tex:2865-2869`:

```
Upsilon'_j : B(H) -> B(L_j),

Upsilon'_j(X) = L_j^dag (Phi(X) tensor 1_F) L_j
```

**CP structure.** This is manifestly CP: write `L_j^dag (Phi(X) tensor 1_F) L_j
= L_j^dag (sum_a K_a^dag X K_a tensor 1_F) L_j`. Expanding
`Phi(X) tensor 1_F = sum_a (K_a tensor 1_F)^dag (X tensor 1_F) (K_a tensor 1_F)`
(using `(K_a tensor 1_F)^dag = K_a^dag tensor 1_F` — the complex conjugate of K_a
acts on H), the map becomes a composition of the CP map `X |-> (X tensor 1_F)`
ampliated and sandwiched by L_j. So Upsilon'_j is a sum of maps
`X |-> (L_j^dag K_a^dag tensor 1)(X tensor 1)(K_a tensor 1 L_j)`... but the
simpler direct argument: `X |-> L_j^dag (sum_a K_a^dag X K_a tensor 1_F) L_j`
is a sum of PSD-valued maps `X |-> (K_a tensor 1_F L_j)^dag X (K_a tensor 1_F
L_j)` — all positive, so CP.

**Practical computation.** `Phi(X) tensor 1_F` is `(N*r) x (N*r)`. Assemble the
Stinespring V once; compute Phi(X) via `aic_ucp_apply`; then form
`kron(Phi(X), 1_F)` (size `(N*r) x (N*r)`, left-major) and sandwich by L_j.
Alternatively, exploit linearity: the Choi matrix of Upsilon'_j is
`C_j = sum_{a,b in [0,N)} E_{ab} tensor Upsilon'_j(E_{ab})` (Convention A,
left-major, source space H of dimension N). This is a `(N * d_j) x (N * d_j)`
PSD matrix since Upsilon'_j is CP; use `aic_ucp_is_cp_choi` for the CP
certificate (as in FINDINGS §C5 for the midpoint+Weyl pattern).

---

### A.6. The full Upsilon': B(H) -> B, assembly and unitalization

`.tex:2865-2899`:

```
Upsilon' = (Upsilon'_1, ..., Upsilon'_m)  :  B(H) -> B = (+)_j B(L_j)
```

where each component Upsilon'_j is as in A.5.

**Unitalization (`.tex:2896-2898`):**

```
Upsilon(X) = (Upsilon'(1_H))^{-1/2} Upsilon'(X) (Upsilon'(1_H))^{-1/2}
```

where `Upsilon'(1_H) = (Upsilon'_1(1_H), ..., Upsilon'_m(1_H))` and each
component `Upsilon'_j(1_H) = L_j^dag (Phi(1_H) tensor 1_F) L_j = L_j^dag
(1_H tensor 1_F) L_j = L_j^dag L_j` (since Phi is unital: Phi(1_H) = 1_H).

So `Upsilon'_j(1_H) = L_j^dag L_j`, a `d_j x d_j` PSD matrix. The paper's
bound `||Upsilon'(1_H) - 1_B||_cb = O(eta)` (implied by `||Upsilon' - tilde_Upsilon||_cb
= O(eta)` and `tilde_Upsilon(1_H) = 1_B`) ensures invertibility for small eta.

**Invertibility assertion (Rule 4, shard H #7).** Assert
`||Upsilon'_j(1_H) - 1_{L_j}||_op < 1` for each j before computing the
inverse-sqrt; abort if not (FINDINGS §D4, escalation 5 for the eta threshold).
The inverse-sqrt of `Upsilon'_j(1_H)` uses `aic_funcalc_xpow` (the
`(-1/2)`-power of a PSD matrix close to identity, `.tex:2897`).

---

## B. lem_RC(i) — Schur/commutant proof intuition and the centrality teeth

### B.1. WHY the 1-design forces R_j = 1_{L_j} tensor C_j

The key is that the 1-design average `T_j(A) = sum_s p_{js} (U_{js}^dag tensor 1)
A (U_{js} tensor 1)` is a PROJECTOR onto the commutant of the algebra `{X tensor 1
: X in B(L_j)}` in `B(L_j tensor E_j)`.

The precise quantum-information version (see e.g. Hayden-Leung-Winter,
"Aspects of Generic Entanglement", CMP 2006, or the 1-design characterization
in Dankert et al., PRA 2009): a set of unitaries `{U_s}` with weights `{p_s}` is
a unitary 1-design on U(d) iff

```
sum_s p_s U_s^dag X U_s = Tr(X)/d * 1_d   for all X in B(L).
```

The centrality property `diag_j2` (`.tex:2776`) is the two-sided version:
`sum_s p_s X U_s^dag tensor U_s = sum_s p_s U_s^dag tensor U_s X` for all X in
`B(L_j)`. This is equivalent to saying the operator `D_j = sum_s p_s U_s^dag
tensor U_s` commutes with every `X tensor 1_{E_j}`.

For the twirl `T_j(A) = sum_s p_s (U_s^dag tensor 1) A (U_s tensor 1)`:
- If X in B(L_j), `T_j(X tensor A') = (sum_s p_s U_s^dag X U_s) tensor A'`.
  By the 1-design property, `sum_s p_s U_s^dag X U_s = (Tr X / d_j) 1_{L_j}`.
  So `T_j(X tensor A') = (Tr X / d_j) 1_{L_j} tensor A'`.
- The image of T_j lands in `{1_{L_j} tensor A' : A' in B(E_j)}` — the L_j-trivial
  subalgebra.
- Since T_j is a conditional expectation onto the commutant, and `1 tensor B(E_j)`
  is exactly the commutant of `B(L_j) tensor 1` by Schur (irreducible action of
  B(L_j) on L_j), T_j maps any operator in B(L_j tensor E_j) into `1 tensor B(E_j)`.

Applied to A = W_j W_j^dag: T_j(W_j W_j^dag) = 1_{L_j} tensor C_j with
`C_j = d_j^{-1} Tr_{L_j}(W_j W_j^dag)` ...which is NOT the direct formula; the
formula from the R_j definition gives `C_j = d_j^{-1} Tr_{L_j}(R_j)`. (They
agree: in the eta=0 case, R_j = T_j(W_j W_j^dag) = 1 tensor C_j, so C_j =
d_j^{-1} Tr_{L_j}(R_j).) For general eta, R_j itself involves U_{js} from the
SAME 1-design used in Delta (a coupling between lem_RC and Delta's 1-design) —
see A.2 above.

### B.2. The centrality teeth F3 must carry (what F2 cannot test)

F2's Delta construction uses the 1-design only through the property `sum_s p_s
Delta'(X) >= 0` (CP by averaging). It does NOT verify the commutation identity
`diag_j2`.

**F3 must test lem_RC(i) has TEETH — three required checks:**

**(T1) Direct centrality (diag_j2) check.** For each j and a random `X in B(L_j)`:

```
|| R_j (X tensor 1_{E_j}) - (X tensor 1_{E_j}) R_j ||_op  <=  O(eta)
```

This is `diag_j2` applied to `W_j W_j^dag`-pieces. For the exact 1-design this
is ZERO (the twirl lands in the commutant exactly). For eta > 0 it is O(eta).

**(T2) Direct product-form check.**

```
|| R_j - 1_{L_j} tensor C_j ||_op  <=  O(eta)
```

where `C_j = d_j^{-1} Tr_{L_j}(R_j)`. This is the constructive version of
lem_RC(i). At eta = 0 this is zero. At eta > 0 it should be O(eta); assert
`ratio = || R_j - 1 tensor C_j || / eta < C_max` for a generous constant.

**(T3) Breaking mutation — tooth for (T2).** Replace `C_j` by
`C_j + epsilon * random_PD` (a perturbation that breaks the product form).
Verify that (T2) is now RED (||R_j - 1 tensor C_j_perturbed|| >> O(eta)).
This confirms (T2) is not vacuously satisfied.

---

## C. Convention cross-checks (partial-trace direction, L_j tensor E_j ordering, Kronecker)

### C.1. What the standard convention says (Watrous)

TQI (Watrous 2018), Chapter 1, uses the following tensor-product convention:
for spaces X = C^m, Y = C^n with orthonormal bases {e_i} resp. {f_j}, the
tensor product X tensor Y has basis {e_i tensor f_j} ordered lexicographically
(i first, j second). A vector u in X tensor Y has index `u[i*n + j]`. The partial
trace over Y: `(Tr_Y A)[i,k] = sum_j A[i*n+j, k*n+j]` — this traces the SECOND
(Y, right) factor, keeping the first (X, left) factor. (TQI, Definition 1.3.6.)

Watrous's Choi-Jamiolkowski isomorphism (TQI, Definition 2.2.1) uses the
convention: for a map F: L(X) -> L(Y), the Choi matrix is
`J(F) = sum_{i,j} E_{ij} tensor F(E_{ij})` where E_{ij} in L(X) are the
standard matrix units. So the LEFT/MAJOR factor is the SOURCE space (X) and the
RIGHT/MINOR factor is the TARGET space (Y).

### C.2. The project's API convention (aic_mat.h / aic_ucp.h)

From `aic_mat.h` docstring:
- LEFT factor MAJOR: `row = i*b + j` with i in [0,a), j in [0,b).
- `aic_mat_partial_trace_right`: traces OUT the RIGHT (b) factor, keeps LEFT (a).
- `aic_mat_partial_trace_left`: traces OUT the LEFT (a) factor, keeps RIGHT (b).

From `aic_ucp.h` docstring:
- Choi matrix Convention A: `C_Phi = sum_{i,j} E_{ij} tensor Phi(E_{ij})` in
  `B(K) tensor B(H)`, with the K (SOURCE/domain) factor LEFT/MAJOR and the H
  (TARGET/codomain) factor RIGHT.
- Unitality check: `Tr_K(C_Phi) = 1_H` — trace out the LEFT (K) factor, keep
  the RIGHT (H) factor.
- The unitality check uses `aic_mat_partial_trace_left` (traces the LEFT/K factor).

**CONCLUSION: the project's Convention A matches Watrous.** Source factor is
left/major, target is right. Partial-trace-left for unitality. This is correct.

### C.3. Partial-trace direction for C_j = d_j^{-1} Tr_{L_j}(R_j)

**The space layout of R_j.** R_j lives in `B(L_j tensor E_j)`. With the
left-major convention: rows and columns indexed as `(i*dim_E_j + a, k*dim_E_j
+ b)` with `i,k in [0,d_j)` (L_j indices) and `a,b in [0,dim_E_j)` (E_j
indices). The product form `R_j = 1_{L_j} tensor C_j` means:
`R_j[i*dim_E_j + a, k*dim_E_j + b] = delta_{ik} C_j[a,b]`.

**What trace gives C_j?** We want `C_j in B(E_j)`:
`C_j[a,b] = (1/d_j) sum_i R_j[i*dim_E_j + a, i*dim_E_j + b]`.
This is `(1/d_j) Tr_{L_j}(R_j)` where the trace is over the LEFT (L_j) factor
of the left-major representation. So we need **`aic_mat_partial_trace_left`**
with `a = d_j` (dim_L) and `b = dim_E_j`, giving a `dim_E_j x dim_E_j` output.

**Formula:** `C_j = (1/d_j) aic_mat_partial_trace_left(R_j, d_j, dim_E_j)`.

**Watrous cross-check.** In Watrous's left-major convention (which the project
follows), R_j's L_j factor is LEFT. Tracing out LEFT = tracing L_j = giving
C_j in B(E_j). This is consistent: `Tr_{L_j}(1_{L_j} tensor C_j) = d_j C_j`,
so `C_j = (1/d_j) Tr_{L_j}(R_j)`. CONFIRMED: use `aic_mat_partial_trace_left`.

**FLAG — the direction is easy to get backwards.** The docstring for
`aic_mat_partial_trace_left` says "traces out the LEFT (C^a) factor, keeps C^b."
For R_j: LEFT = L_j (the factor being traced, dimension d_j), RIGHT = E_j
(the factor to keep, dimension dim_E_j). Pass `a = d_j`, `b = dim_E_j`.
Passing `(a, b) = (dim_E_j, d_j)` would trace E_j and give an d_j x d_j object
— NOT C_j. This is a swapped-factor bug that is INVISIBLE to the eta=0 oracle
when d_j = dim_E_j (equal-size blocks). **A non-square (d_j != dim_E_j) fixture
is needed to catch this (F3 must have such a fixture).**

### C.4. The L_j tensor E_j ordering in W_j

`.tex:2837`: `W_j: H -> L_j tensor E_j`. In the left-major convention, the
target of W_j is the `(d_j * dim_E_j)`-dimensional space with L_j LEFT and E_j
RIGHT. The row index of `V_Delta` for the j-th block is `i*dim_E_j + a` with
`i in [0,d_j)` (L_j) and `a in [0,dim_E_j)` (E_j). This must match the layout
assumed in the `U_{js} tensor 1_{E_j}` application:
`(U_{js} tensor 1_{E_j})[i*dim_E_j + a, k*dim_E_j + b] = U_{js}[i,k] delta_{ab}`,
which is correct under left-major with L_j LEFT.

**FLAG (potential convention mismatch).** If W_j is extracted from Delta's Kraus
ops in a different block order (E_j LEFT, L_j RIGHT), then `(U_{js}^dag tensor
1_{E_j})` and the partial-trace-left direction are both wrong. The safest
construction: build W_j EXPLICITLY as `W_j = column_stack of {K_{j,a,b}}` where
each Kraus op `K_{j,a,b}: H -> L_j` is right-padded to `H -> L_j tensor E_j` by
the E_j basis vector `e_b`, so `(W_j v)[i*dim_E_j + b] = K_{j,i,b}[row] * v`.
Then L_j is LEFT by construction.

The `aic_ucp.h` Stinespring convention (`V[a*dim_K + i, j] = K_a[i,j]`) puts the
Kraus INDEX (a) first and the target space (dim_K) second. For Delta's Stinespring
`V_Delta: H -> (+)_j (L_j tensor E_j)`, the row index in the j-th block must be
`i_Lj * dim_E_j + i_Ej` (L_j left, E_j right). The `aic_ucp.h` stacking
convention would put the Kraus index a first and i_K = i_{L_j} second, giving
row `a*d_j + i_Lj` — this is E_j LEFT if `a` corresponds to the E_j index and
i_Lj to L_j. **DIVERGENCE FLAG:** the `aic_ucp.h` Stinespring stacking puts the
ancilla F index first (`a * dim_K + i`, i.e. F LEFT), whereas for the L_j tensor
E_j decomposition of Delta, we need L_j LEFT. If E_j = F (same ancilla), the
two conventions conflict: F LEFT in V, but L_j LEFT needed for W_j.

**Resolution.** The W_j are NOT the Stinespring V restricted to one block of the
target space with a reordering — they are DEFINED by the Choi decomposition of
Delta (`.tex:2834`). The right construction is to work directly from Delta's
Kraus operators grouped by block. Each Kraus operator `G_s` of Delta maps
`H -> B` (as an ambient N x n_B matrix); project onto the j-th block to get
the contributions to W_j. Concretely: for each Kraus op `G_s` of Delta (N x N,
selecting the j-th diagonal block gives a `d_j x N` matrix),
define `W_j` as the column-stack of these `d_j x N` blocks over all s — giving
a `(d_j * r_j) x N = (dim_L tensor dim_E_j) x N` matrix with dim_E_j = r_j
(number of Kraus ops contributing to block j). This places L_j (dim d_j)
LEFT and E_j (dim r_j) RIGHT, consistent with left-major.

**Conclusion:** build W_j by blocking Delta's Kraus operators per target block j,
stacking as L_j-LEFT, E_j-RIGHT. Partial trace uses `aic_mat_partial_trace_left(R_j, d_j, dim_E_j)`.

### C.5. The `(Delta(U_{js}^dag) tensor 1_F)` in L_j

V is Phi's Stinespring: `V: H -> H tensor F`, shape `(N*r) x N`, with H LEFT
and F RIGHT in V's row index (`a*N + i` for Kraus index a, H index i). So V
stores F LEFT (`a`) and H RIGHT (`i`) — WAIT: from `aic_ucp.h`:
`V[a*dim_K + i, j] = K_a[i, j]`, so row index `a*dim_K + i`. With dim_K = N
(self-map), row = `a*N + i` with `a` the ANCILLA index (F) and `i` the
domain-space index (H). So F is LEFT (Kronecker outer factor), H is RIGHT (inner
factor) in V's row space. This means V: H -> F tensor H (in the Kronecker
ordering), NOT H -> H tensor F as the paper states!

**DIVERGENCE FLAG.** The paper writes `V: H -> H tensor F` (`.tex:2859`), placing
H LEFT and F RIGHT. But `aic_ucp.h`'s Stinespring packs F (Kraus index) LEFT and
H (domain) RIGHT: `V[a*N+i, j] = K_a[i,j]` gives V: H -> F tensor H in the
left-major convention. This is a convention flip on V's output space.

For the product `V W_j^dag (U_{js} tensor xi_j)` in L_j's formula:
- `(U_{js} tensor xi_j): L_j -> L_j tensor E_j` (L_j LEFT, E_j RIGHT) maps
  `u |-> U_{js} u tensor xi_j`, shape `(d_j * dim_E_j) x d_j`.
- `W_j^dag: L_j tensor E_j -> H` is `N x (d_j * dim_E_j)` (H = target).
- `V: H -> F tensor H` (in project convention) is `(N*r) x N`.
- So `V W_j^dag (U_{js} tensor xi_j): L_j -> F tensor H`, shape `(N*r) x d_j`.

For `(Delta(U_{js}^dag) tensor 1_F)` to act on this, we need it to act on the
F tensor H space: `(Delta(U) tensor 1_F)` in the project convention should be
`kron(1_F, Delta(U))` — the 1_F acts on the LEFT F factor and Delta(U) acts on
the H (RIGHT) factor. But in left-major Kronecker, `kron(A, B)` puts A LEFT
and B RIGHT, so the correct matrix is `kron(1_F, Delta(U_{js}^dag))` (shape
`(N*r) x (N*r)`).

**THIS IS A BUG-PRONE SPOT.** The paper's `.tex:2862` writes
`(Delta(U_{js}^dag) tensor 1_F)` but in the project's left-major convention with
F LEFT in V, the correct matrix product is `kron(1_F, Delta(U_{js}^dag))`.
Writing `kron(Delta(U), 1_F)` would put Delta(U) on the F factor (wrong).

**Recommended probe.** At eta = 0: Phi = id, V = 1_H (r=1, single Kraus = 1_H,
Stinespring = column-vector of 1_H = 1_N), F = 1-dim. Then `kron(1_F, Delta(U))
= Delta(U)` (a scalar times 1, trivial). The non-trivial test is eta > 0 with
r > 1. Build both `kron(1_F, Delta(U))` and `kron(Delta(U), 1_F)`, verify that
only the correct one gives `||Upsilon Delta - 1_B|| = O(eta)`.

**Summary of C section.** The identified potential divergences:

| Convention point | Paper | Watrous TQI | Project API | Verdict |
|---|---|---|---|---|
| Choi convention (source LEFT) | E_{ij} tensor Phi(E_{ij}) | same | Convention A, source LEFT | AGREE |
| Partial trace for unitality | Tr_K(C) = 1_H | Tr_Y (right) | `partial_trace_left` (traces left/K) | AGREE (API correctly traces source = left) |
| C_j extraction | Tr_{L_j}(R_j)/d_j | trace 1st factor | `partial_trace_left(R_j, d_j, dim_E_j)` | AGREE, but verify L_j is LEFT in R_j layout |
| V: H -> H tensor F | H LEFT, F RIGHT | H LEFT, F RIGHT | V packs F LEFT (`a*N+i`), H RIGHT | DIVERGENCE: project has F LEFT |
| `(Delta(U) tensor 1_F)` in L_j | Delta LEFT, 1_F RIGHT | same | needs `kron(1_F, Delta(U))` given F LEFT | DIVERGENCE / BUG-PRONE |

---

## D. Test plan

### D.1. eta = 0 oracle (exact idempotent)

When eta = 0, Phi is exactly idempotent. By th_idemp_structure (`.tex:318`), the
exact factorization is Phi = Delta Gamma C_M with Delta = inclusion, Gamma = the
conditional expectation, and Upsilon = Gamma C_M.

F3's Upsilon must reduce to `Gamma C_M` in this limit (shard H cross-check):

**(O1) Upsilon(X) = Gamma(C_M(X)) + O(machine_eps).** Compare F3's Upsilon output
against `aic_idemp_decomp.Gamma` applied to the compression `C_M(X) = J_M^dag X J_M`.
Gate: `max_ij ||Upsilon(E_ij) - Gamma(C_M(E_ij))|| < 1e-10`.

**(O2) Upsilon Delta = 1_B exactly.** `|| (Upsilon Delta)(E_ij) - E_ij ||_op < 1e-10`
for all basis elements E_ij of B. This is a HARD requirement; at eta = 0 there
is no O(eta) slack.

**(O3) Delta Upsilon = Phi exactly.** `|| (Delta Upsilon)(E_ij) - Phi(E_ij) ||_op < 1e-10`.

**(O4) ||C_j|| = 1 exactly.** For eta = 0, `sigma_max(C_j) = 1` to machine precision.

**(O5) R_j = 1 tensor C_j exactly.** `|| R_j - kron(1_{L_j}, C_j) ||_op < 1e-10`.

Fixtures: `make_block_cond_exp(d, m)` (noiseless subsystem, the standard oracle
from `test_opspace.c`). Also `aic_idemp_decompose` from `include/aic_idemp.h`.

### D.2. eta > 0 (mixconj fixture)

Use `make_mixconj(n, m, t)` (the standard near-identity-idempotent fixture):

**(H1) ||Upsilon Delta - 1_B||_cb = O(eta).** Measure ratio `r1 = ||Upsilon Delta -
1_B||_cb / eta`; assert `r1 < C_H1 = 20` (generous fail-loud gate). Certify via
the Choi difference + `aic_ucp_is_cp_choi`-style estimate.

**(H2) ||Delta Upsilon - Phi||_cb = O(eta).** Measure ratio `r2 < C_H2 = 20`.

**(H3) sigma_max(C_j) >= 1 - C_RC * eta** for an explicit `C_RC`. Assert this per
block j (Rule 4 abort if violated — a lem_RC precondition failure).

**(H4) Upsilon is UCP.** (a) `||sum_j L_j^dag L_j - 1_H||_op < tol` (unitality of
Upsilon'_j(1_H) before normalization). (b) CP via per-block Choi PSD check
(FINDINGS §C5 midpoint+Weyl pattern as in F2).

**(H5) ||Upsilon - tilde_Upsilon||_cb = O(eta).** Measure `||Upsilon - tilde_Upsilon||_cb / eta < C_H5`.
tilde_Upsilon is from `aic_factorize_upsilon_tilde`.

### D.3. Centrality teeth (lem_RC)

**(C1) R_j commutes with X tensor 1_{E_j}** (see B.2 T1 above). At eta = 0:
exact to machine precision. At eta > 0: O(eta). Use random Hermitian X in B(L_j).

**(C2) ||R_j - 1_{L_j} tensor C_j||_op / eta < C_C2 at eta > 0.** (B.2 T2.)

**(C3) Breaking mutation (B.2 T3).** Replace C_j with `C_j + eps_mut * rand_PD`;
assert (C2) fires RED (defect grows by O(eps_mut) which >> O(eta) for eps_mut chosen
appropriately).

**(C4) Non-square fixture (partial trace direction).** At least one fixture with
`d_j != dim_E_j` to catch the swapped-factor bug (C.3).

**(C5) Stinespring convention probe (C.5 divergence).** Directly verify that
`kron(1_F, Delta(U))` is the correct operator in L_j (not `kron(Delta(U), 1_F)`)
by checking `||Upsilon Delta - 1_B||_op` with both orderings; the incorrect
ordering gives O(1) defect.

### D.4. Star tooth (FINDINGS §C2 / §C8)

The star-blindness finding applies: the eta = 0 oracle is blind to whether the
Choi-Effros star is used correctly. F3 uses `aic_factorize_delta_prime` (Delta'
via Phi, the CP-ized version) which internally calls `aic_ucp_apply` — this is
the PLAIN matrix product inside Phi (the observable form), NOT the star. The star
appears in verifying `Delta'(X) approx tilde_Delta(X) approx v(X)`. Confirm:
mutate Delta' -> plain ambient product (remove the Phi wrapper) and confirm
||Upsilon Delta - 1_B|| goes from O(eta) to O(1).

---

## E. Ranked risks and file/LOC split

### E.1. Ranked risks

| Rank | Risk | Source | Mitigation |
|---|---|---|---|
| 1 | `kron(Delta(U), 1_F)` vs `kron(1_F, Delta(U))` in L_j (C.5) | V's left-major F-LEFT vs paper's H-LEFT convention | explicit convention probe (D.3 C5); the eta=0 oracle is BLIND here if r=1 |
| 2 | Partial-trace direction for C_j (C.3/C.4) | `aic_mat_partial_trace_left` vs `_right`, easy to swap | non-square fixture d_j != dim_E_j (D.3 C4) |
| 3 | W_j block ordering (L_j LEFT vs E_j LEFT, C.4) | Stinespring stacking convention vs `.tex:2837` | build W_j from Delta's Kraus blocks directly, NOT from V_Delta's stacking |
| 4 | `Upsilon'_j(1_H) = L_j^dag L_j` invertibility | eta threshold for shard H #7 | assert `||L_j^dag L_j - 1_{L_j}||_op < 1` per j; abort if not (Rule 4) |
| 5 | sigma_max(C_j) < 1 - O(eta) abort | lem_RC(ii) precondition; fail loud | Rule 4 abort with clear message citing `.tex:2846` |
| 6 | Composite O(eta) constant untracked | FINDINGS §D4 escalation 5 | measure ratio ||Upsilon Delta - 1_B|| / eta per instance; assert bounded + no-trend canary (§D2 robust form) |
| 7 | Degenerate singular values in C_j | `aic_mat_singular_values` requires distinct values | use LAPACK `zgesvd` (double path) for xi_j extraction; C_j SVD is small (dim_E_j x dim_E_j) |
| 8 | `aic_mat_opnorm` Gram false-fail (FINDINGS §C5) | accumulated arb radius in L_j^dag L_j | apply §C5 midpoint+Weyl pattern on `Upsilon'_j(1_H)` before CP check |

### E.2. File/LOC split (each <= 200 LOC, CLAUDE.md Rule 10)

Suggested split for F3 (extends `aic_factorize.h`/`aic_factorize_delta.c`):

**`src/aic_factorize_upsilon.c`** (~180 LOC):
- `aic_factorize_lemRC_R_j(R_j, F, j, prec)` — compute R_j for block j.
- `aic_factorize_lemRC_C_j(C_j, R_j, d_j, dim_Ej, prec)` — partial-trace-left
  extraction of C_j.
- `aic_factorize_lemRC_xi_j(xi_j, C_j, prec)` — leading right singular vector
  (LAPACK double path).
- `aic_factorize_L_j(L_j, F, V, W_j, U_js, xi_j, j, s, prec)` — single term
  of the L_j sum.
- Internal helpers for the sum over s.

**`src/aic_factorize_upsilon_build.c`** (~160 LOC):
- `aic_factorize_upsilon_prime_j(out, F, W_j, L_j, X, j, prec)` — Upsilon'_j(X).
- `aic_factorize_upsilon_build(F, prec)` — precompute Upsilon'(1_H)^{-1/2} per
  block; store in new `upsilon_invsqrt[j]` fields.
- `aic_factorize_upsilon(out, F, X, prec)` — the final UCP map.
- Assertions: ||C_j|| lower bound, ||Upsilon'(1_H) - 1_B|| < 1.

**`tests/test_factorize_upsilon.c`** (~190 LOC):
- T1: eta=0 oracle (O1-O5 above).
- T2: eta>0 mixconj fixture (H1-H5 above).
- T3: centrality teeth (C1-C5 above).
- T4: convention probe for kron order (C5).
- T5: dimension-independence canary for composite constant (§D2 robust form).

**Header additions to `include/aic_factorize.h`** (~30 LOC):
- `aic_factorize_upsilon_build(F, prec)`.
- `aic_factorize_upsilon(out, F, X, prec)`.
- New struct fields for `upsilon_invsqrt[j]` (per-block) and `upsilon_ready`.

---

## Summary: 15-line verdict for the calling agent

**Partial-trace direction (C_j).** Use `aic_mat_partial_trace_left(R_j, d_j, dim_Ej)`
— this traces the L_j (left) factor, leaving C_j in B(E_j). This is consistent
with Watrous's convention AND the project's left-major layout since L_j is the
left factor of `L_j tensor E_j` in R_j. CONFIDENT. The non-square fixture (d_j
!= dim_E_j) is the mandatory mutation catch.

**W_j route.** Build W_j by grouping Delta's Kraus operators by target block j,
stacking with L_j LEFT and E_j RIGHT. Do NOT derive W_j from the Stinespring
`V_Delta` column-stack without verifying the L_j/E_j ordering matches (the
`aic_ucp.h` ancilla-LEFT convention may conflict). CONFIDENT on the partial-trace
direction; EMPIRICAL PROBE REQUIRED on the W_j blocking.

**Centrality teeth design.** Three checks needed: (1) R_j commutes with X tensor 1;
(2) ||R_j - 1 tensor C_j|| = O(eta); (3) breaking mutation (perturb C_j, confirm
RED). F2 cannot test lem_RC(i) — it does not use R_j at all. F3 is the first
module that exercises the Schur/commutant identity.

**Convention divergences from standard QI (the bug-prone spots).**
(a) `V`'s Stinespring convention in `aic_ucp.h` packs F (Kraus index) LEFT and
H (domain) RIGHT — the paper writes `V: H -> H tensor F` (H LEFT, F RIGHT). This
is a pure convention difference but it means `(Delta(U) tensor 1_F)` in L_j's
formula must be `kron(1_F, Delta(U))` in the project's left-major with F-LEFT
— NOT `kron(Delta(U), 1_F)`. This is a HIGH-RISK point: it does not show up at
eta = 0 with r = 1 (trivial F), and the eta > 0 oracle test must explicitly
compare both orderings. (b) The partial-trace direction for unitality in `aic_ucp.h`
is correctly `partial_trace_left` (tracing the source/K factor); the same
pattern applies for C_j. These two are the only identified divergences; both are
indexing/ordering bugs rather than mathematical errors.

**Confidence summary.** High confidence on: Schur/commutant structure (B.1),
partial-trace formula for C_j (C.3), xi_j via leading right singular vector, Upsilon'_j
manifestly CP, unitalization form. Medium confidence (needs empirical probe) on:
the kron ordering in L_j's formula (C.5 — the single highest risk), the W_j
block-extraction layout (C.4). The eta=0 oracle is partially blind to both of
these (F-space trivializes at r=1; equal d_j=dim_E_j hides the partial-trace
direction). Design the test suite so these are NOT blind spots.
