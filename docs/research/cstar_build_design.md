# Design Document: `cstar_build` — Constructive proof of `th_main`

**Module bead:** `aic-097`
**Tex citation:** `approximate_algebras.tex:1321–1444` (§9)
**Depends on:** `corner` (aic-czm, CLOSED), `projection` (aic-mqf, CLOSED),
`dhom` (aic-c1n, CLOSED), `errreduce` (aic-t81, CLOSED)
**Date:** 2026-05-30

---

## 0. Reading-order note

This document verifies and refines a proposed design against the `.tex`. Cite
format: `.tex:LINE`. Every API reference is to an existing, tested header.
"Closed" headers = green test suite, no sorries.

---

## 1. Verdict on the S_P-as-aic_ecstar wrapper

### 1.1 What the proposed wrapper says

The loop's inner objects are compressed subalgebras S_{P_m} that must be
treated as eps'-C* algebras in their own right (so that `projection` and
`lem_extension` can be called recursively on them). The proposal wraps each
S_P as a **derived** `aic_ecstar` by:

- `n` = same as parent A (all operators still n x n ambient),
- basis `{C_m}` = corner-extracted Frobenius-orthonormal basis of S_P
  (the top-dim_S left singular vectors of Co_P, via `aic_corner_extract`),
- `dim_A` = dim S_P,
- star: `star_phi` = a thunk computing `Co_P(Phi_tilde_parent(XY))` for
  arguments X, Y in S_P,
- unit: `Ptilde = Co_P(P) = aic_corner_Ptilde(A_parent, P)`.

### 1.2 Verification against the .tex

**The star formula.** The compressed product `.tex:1080` gives
```
  X . Y = Co_{P,R}(X * Y) = Co_{P,P}(Phi_tilde(XY))   for X, Y in S_P
```
(using Co_{P,P} = Co_P and * = the parent star = Phi_tilde). The compressed
product is also what `aic_corner_cdot` computes, via its `CoPR` argument.
So the proposed `star_phi` thunk is literally:

```c
  out = aic_corner_cdot(A_parent, Co_P, X, Y)
      = Co_P applied (via aic_corner_apply) to (X * Y)
      = Co_P(Phi_tilde(XY))     (the compressed product, .tex:1080)
```

This is **correct**: it is the Choi-Effros star of the compressed algebra
S_P, which is an O(delta+eps)-C* algebra by `.tex:1082` ("the compressed
product turns S_P into an O(delta+eps)-C* algebra with unit Ptilde").

**The unit.** `.tex:1082` gives the unit as `Ptilde = Co_P(P)`, exactly
what `aic_corner_Ptilde` computes. This is the unit of (S_P, .), **not**
the ambient `1_n`. This distinction is load-bearing for several downstream
steps; see §2.2.

**VERDICT: The wrapper's data model is correct.** The n stays the same
(all matrices are n x n), the basis is the extracted corner basis, and the
star is the compressed product. The `aic_ecstar` struct supports this through
the `star_phi`/`star_ctx` seam (the same path used by `aic_assoc_ecstar`
for non-Kraus stars).

### 1.3 Concretely: what to store in the wrapper struct

The S_P wrapper requires the following owned data:
- the d x d matrix `Co_P` (= `aic_corner_Co(A_parent, P, P)`), which serves
  as the idempotent that defines the compression,
- a `ctx` struct containing: a pointer to `A_parent` and the `Co_P` matrix,
- `star_phi` = `sp_star_thunk` (a static function computing
  `aic_corner_apply(A_parent, Co_P, tmp)` where `tmp = aic_ecstar_star(A_parent, X, Y)`).

The `n`, `dim_A`, `B[]` fields are filled by:
```
  aic_corner_extract(&C_out, &dim_S, Co_P, A_parent, prec)
```
giving the n x n `C_m` operators. No new n^2 x n^2 superoperator is needed.

---

## 2. The projection-integration claim

### 2.1 The claim

The projection module's `project_into_A` uses `W * 1_n` (hardcoded to the
ambient identity, `aic_projection.c` step 4: `P = Phi_tilde(P_amb) = P_amb * I`
where I = 1_n). The claim is: with the S_P wrapper, this automatically gives
the correct projection into S_P, because the wrapper's star computes

```
  P_amb * I  (star in S_P wrapper)
  = Co_P(Phi_tilde(P_amb · 1_n))   (compressed product with 1_n)
  = Co_P(Phi_tilde(P_amb))         (since Phi_tilde(X · 1_n) = Phi_tilde(X))
  = Co_P(P_from_parent)
```

which is the oblique projection of `P_amb` into S_P via Phi_tilde followed
by Co_P. The claim is that this gives the **correct oblique projection** into
S_P (not the Frobenius projector Pi_A, which would be wrong per FINDINGS §C3).

### 2.2 Verification — step 4 (projecting into A)

**The parent's projector.** In `aic_projection.c`, step 4 computes
`P = W * I_A` where I_A is the unit of A and * is the star of A. With
the S_P wrapper, the unit field of the ecstar struct is `Ptilde = Co_P(P)`,
NOT `1_n`. So `W * I_{S_P_wrapper}` computes:

```
  W * Ptilde  (in the S_P wrapper's star)
  = Co_P(Phi_tilde(W · Ptilde))
  = Co_P(Phi_tilde(W) · Ptilde) + O(eps)      (by .tex:1082: S_P is eps-C*)
  ≈ Co_P(Co_P(Phi_tilde(W)))
  = Co_P(Phi_tilde(W))
```

because Co_P^2 = Co_P (Co_P is an exact idempotent after aic_corner_Co
applies prop_P). This gives the oblique projection of W into S_P, with
O(delta+eps) error from the eps-C* approximations.

**CRITICAL CORRECTION:** The projection source code computes `P_amb * I`
where `I` is the unit passed to `aic_projection_nontrivial`. Looking at the
header:

```c
void aic_projection_nontrivial(acb_mat_t P, arb_t delta, arb_t Pnorm,
                               arb_t ImPnorm, const aic_ecstar *A, ...)
```

Step 4 in `aic_projection.c` uses `aic_ecstar_star(A, P_amb, I_n)` where
`I_n` is the n x n identity built inside the function (not the wrapper's unit
Ptilde). Let us check whether the claim that this gives the right projection
into S_P is actually sound.

From `.tex:2211` (via the assoc module design): `X * I_A = Phi_tilde(X)` when
the unit is the ambient `1_n` and the star is `Phi_tilde`. For the S_P
wrapper:

```
  aic_ecstar_star(S_P_wrapper, P_amb, 1_n)
  = Co_P(Phi_tilde(P_amb · 1_n))
  = Co_P(Phi_tilde(P_amb))
```

which IS the oblique projection of P_amb into S_P. **The claim is sound.**
The key point: the S_P wrapper's star_phi is `Co_P(Phi_tilde(·))`, so
multiplying by the ambient `1_n` (not the wrapper's Ptilde) gives
`Co_P(Phi_tilde(P_amb · 1_n)) = Co_P(Phi_tilde(P_amb))`, which is exactly
the oblique projection into S_P. This is the correct route (FINDINGS §C3).

**However: there is a subtlety with `I_n` vs `Ptilde`.** The projection
code builds the ambient `I_n` (the n x n identity). It does NOT use the
ecstar's unit field. This is load-bearing: if the code were changed to use
the ecstar's stored unit (which is `Ptilde` for the S_P wrapper), it would
give `Co_P(Phi_tilde(P_amb · Ptilde)) ≈ Co_P(Co_P(Phi_tilde(P_amb))) =
Co_P(Phi_tilde(P_amb))`, the same result up to O(eps). Either route works
for the output; the ambient-1_n route is simpler (no approximation). The
current source code (ambient 1_n) requires NO changes for the S_P wrapper.

### 2.3 Other projection steps with the S_P wrapper

**Step 1: pick_H over basis `{C_m}`.** The S_P wrapper's basis is the
extracted `{C_m}` (Frobenius-orthonormal n x n operators in S_P). The
Hermitian elements `H_k = (C_k + C_k^dag)/2` are computed and have
their ambient eigenvalues found by LAPACK zheev. This is unaffected by
whether the basis spans S_P or A — the ambient matrices are just n x n.
**Sound, no change needed.**

**Step 2: spectral gap finding, step 3: ambient sgn.** These operate on
the Hermitian matrix H in M_n (the ambient exact C* algebra), entirely
via LAPACK and `aic_sgn`. FINDING §C1 (rem_X2: no funcalc inside A)
is automatically satisfied because the funcalc is in M_n. **Sound.**

**Step 5: certify.** The star-defect `||P * P - P||_op` uses
`aic_ecstar_star(S_P_wrapper, P, P)` which is `Co_P(Phi_tilde(P^2))`.
The membership test `aic_ecstar_proj_residual(S_P_wrapper, P)` computes the
Frobenius residual against the basis `{C_m}` of S_P. Both are well-defined.
**Sound.**

**Nontriviality bounds: `||P||_op` and `||I - P||_op`.** The nontriviality
in the S_P sub-algebra means `||P||_op` and `||Ptilde - P||_op` are bounded
away from 0, where Ptilde is the unit of S_P. The current code checks
`||I_n - P||_op` (ambient). For the S_P wrapper this is correct at the
nontriviality level: if P is nontrivial in S_P, then Co_P(P) ≠ 0 and
Co_P(P) ≠ Ptilde (both nonvanishing), and the ambient norms approximate
these up to O(delta+eps). There is a minor RISK that the nontriviality
check `||1_n - P||_op >= 0.3` could pass trivially for a "nontrivial in S_P
but trivial in A" projection. This must be guarded: the test should check
`||Ptilde - P||_op >= 0.3` (S_P unit minus P, not ambient 1_n minus P).
**CORRECTION NEEDED:** the projection code's nontriviality check uses the
ambient 1_n; for S_P wrappers, `||Ptilde - P||_op` is the correct comparand.
This requires either modifying `aic_projection_nontrivial` to accept a unit
argument, or the cstar_build loop checking nontriviality itself after calling
the projector.

**VERDICT: The projection integration claim is MOSTLY SOUND**, with one
exception: the nontriviality check needs to use the S_P unit Ptilde, not
the ambient 1_n. All other steps are compatible with the S_P wrapper.

---

## 3. The aic_dhom_v direct-construction claim

### 3.1 The proposal

`cor_merge_sum` concatenates two delta-inclusions `v_1 : B_1 -> S_{P_1}`
and `v_2 : B_2 -> S_{P_2}` by building a combined `aic_dhom_v` over
`B = B_1 (+) B_2` (a new `aic_dhom_B` with the concatenated block-dim
list) and `vE[]` = the concatenated basis images.

### 3.2 Verification against the headers

**`aic_dhom_B`**: block dims `d[0..m-1]`, offsets, `dim_B`, `n_B`. For
`B = B_1 (+) B_2` the block list is the concatenation of the two block lists.
`aic_dhom_B_init(B, dims, m)` accepts any dims array; concatenation is just
building the merged dims array. **Sound.**

**`aic_dhom_v`**: stores `B` (borrowed), `A` (borrowed), `n = A->n`, and
`vE[0..dim_B-1]` (one n x n operator per matrix unit of B). For the merged
`B_1 (+) B_2`, the merged `vE[]` is exactly `[v_1.vE[0..dim_B1-1],
v_2.vE[0..dim_B2-1]]`. `aic_dhom_v_init` allocates `dim_B` operators.
After `init`, the caller writes into `vE[i]` directly. This is legal:
`vE` is a flint_malloc'd array of `acb_mat_t`, each n x n, writable.
**Sound — direct field assignment is the intended construction.**

**CORRECTION (minor):** The merged `A` must be the PARENT A (not S_{P_1}
or S_{P_2}), since both `v_1` and `v_2` map into a subspace of the parent.
The `aic_dhom_v` struct stores a single `A`; for a cor_merge_sum result that
lives in the parent, `A` should be the parent. If `v_1` and `v_2` have
different `A` fields (pointing to their respective S_P wrappers), the merged
v must point to the parent A. This is an API-level bookkeeping issue, not
a data-structure gap.

---

## 4. Assembly lemmas — constructive algorithms as existing-API call sequences

### 4.1 lem_add_dim (.tex:1363–1369)

**Statement (faithful).** Let v: B -> A and w: C -> A be non-unital
delta-inclusions from commutative C* algebras with projection bases
{Pi_j} and {Sigma_k}. Set P_j = v(Pi_j), Q_k = w(Sigma_k), P = v(I_B),
Q = w(I_C). Then dim S_{P,Q} = sum_{jk} dim S_{P_j,Q_k}.

**Constructive algorithm.** This is purely a DIMENSION QUERY; no new map is
returned by this lemma (the proof constructs a bijection internally). The
algorithm is:

1. For each (j,k): call `aic_corner_dim_S(A, P_j, Q_k, prec)` — returns
   `round(Re Tr Co_{P_j,Q_k})`, cross-checked against SVD-gap count.
2. Sum: `N = sum_{jk} dim S_{P_j,Q_k}`.
3. Call `aic_corner_dim_S(A, P, Q, prec)` — returns `d_PQ`.
4. ASSERT `N == d_PQ` (fail loud, the lem_add_dim conclusion).

Or use `aic_corner_alpha_dims(A, P_parts, p, P, Q_parts, q, Q, N_out, dPQ_out, prec)`,
which does steps 1–4 internally and returns both N and d_PQ. **Already available.**

**Stage-3 use (S_{P_C, P_D} = 0).** For distinct equivalence classes C, D:
dim S_{P_C, P_D} = sum_{j in C, k in D} dim S_{P_j, P_k}. Since j !~ k,
each dim S_{P_j, P_k} = 0 by definition. So:

```
  aic_corner_dim_S(A, P_C, P_D, prec) -> ASSERT 0
```

is the certification. This uses lem_add_dim implicitly: the sum over (j,k)
is 0, and lem_add_dim says the total dim equals the sum.

**eta=0 oracle.** For an exact idempotent Phi, dim S_{P,Q} = rank of
Co_{P,Q}. In exact C* algebra M_n with orthogonal P, Q: dim S_{P,Q} =
rank(P) * rank(Q). For 1-dim P_j and 1-dim Q_k: dim S_{P_j,Q_k} in {0,1},
and the sum equals dim S_{P,Q} exactly (the bijection is isometric).

**Teeth (what must turn RED).** A mutation replacing `aic_corner_dim_S(A,
P_j, Q_k)` with `aic_corner_dim_S(A, P, Q)` (ignoring the sum) would either
over- or under-count. The assertion `N == d_PQ` fires it RED.

### 4.2 cor_merge_sum (.tex:1352–1358)

**Statement.** Let P_1, P_2 be delta-projections with ||P_1+P_2-I|| <= delta,
and v_j: B_j -> S_{P_j} be delta-inclusions. Then
v: (X_1, X_2) |-> v_1(X_1) + v_2(X_2) : B_1 (+) B_2 -> A
is an O(delta+eps)-inclusion.

**Constructive algorithm** (direct construction — no new API function call
for the "merging" itself; it is set up, then verified):

1. Build `B_12 : aic_dhom_B` with block dims = concat(B_1.d, B_2.d).
2. `aic_dhom_v_init(&v12, &B_12, A_parent, prec)` (allocates `dim_B12` slots).
3. Copy `v1.vE[0..dim_B1-1]` into `v12.vE[0..dim_B1-1]` (acb_mat_set).
4. Copy `v2.vE[0..dim_B2-1]` into `v12.vE[dim_B1..dim_B12-1]`.
5. Certify O(delta+eps)-inclusion: call `aic_dhom_defect_sweep` (mult defect)
   and `aic_dhom_v_sigma_min` (inclusion lower bound).
6. If bijective: ASSERT dim S_{P_1,P_2} == 0 (aic_corner_dim_S) and
   dim_B1 + dim_B2 == dim A (bijectivity condition).

**The merging conditions (.tex:1338–1344) for cor_merge_sum.** The lemma's
lem_merging is invoked with gamma_{11}=v_1, gamma_{22}=v_2, gamma_{12}=0,
gamma_{21}=0. The conditions (merging0–3) reduce to:
- (merging0): v_j(X^dag) = v_j(X)^dag — the star-compatibility of v_j.
- (merging1): ||v_j(XY) - v_j(X) . v_j(Y)|| <= delta||X||||Y|| — the
  delta-homomorphism condition of v_j (already stored as its sweep defect).
- (merging2): ||v_j(Pi_j) - P_j|| <= delta — the unit condition of v_j.
- (merging3): (1-delta)||X|| <= ||v_j(X)|| <= (1+delta)||X|| — the
  inclusion bounds, certified by sigma_min + norm_ub.
All of these are ALREADY STORED in v1, v2 (the delta-inclusion was
verified by aic_errreduce). No additional computation beyond copying the vE[].

**eta=0 oracle.** For exact C* structures, v = v_1 (+) v_2 is an exact
isomorphism (zero defect). `aic_dhom_defect_sweep` returns machine-zero.

**Teeth (star-vs-plain blindness, FINDINGS §C2).** The multiplicativity
defect `G_v(E_i, E_j) = v(E_i E_j) - v(E_i) * v(E_j)` uses the STAR in A.
Mutating `*` -> plain matmul in `aic_dhom_Gv` leaves G_v ~0 at eta=0 but
diverges at eta>0. Test fixture must use an oblique (eta>0) channel; the
`make_compress_idemp` family in existing tests covers this.

### 4.3 lem_merging (.tex:1325–1350)

**Statement.** Given gamma_{jk}: S_{Pi_j,Pi_k} -> S_{P_j,P_k} satisfying
(merging0–3), the combined map gamma: B -> A is an O(delta+eps)-inclusion.

**Constructive algorithm.** lem_merging is the MOST COMPLEX of the three
pre-extension lemmas. Its proof uses the canonical bijection mu and lem_alpha;
we must explicitly compute gamma and verify the inclusion.

The combined map gamma is:
```
  gamma(X) = sum_{jk} gamma_{jk}(X_{jk})
```
where X_{jk} = Pi_j X Pi_k are the block components of X in B (exactly
computable since B is a genuine C* algebra with exact Pi_j).

**Constructive steps:**

1. Extract the block decomposition: for each matrix unit E_{mu(l,a,b)} of B,
   determine which block (j,k) it falls in (j corresponds to Pi_1 sector,
   k to Pi_2 sector). In the 2x2 case (Pi_1, Pi_2 complementary):
   - E_i in S_{Pi_1,Pi_1}: block (1,1), so gamma maps it via gamma_{11}.
   - E_i in S_{Pi_2,Pi_2}: block (2,2), via gamma_{22}.
   - E_i in S_{Pi_1,Pi_2} or S_{Pi_2,Pi_1}: off-diagonal.

2. Build the combined `aic_dhom_v`: for each matrix unit E_i of B, compute
   `gamma(E_i) = gamma_{jk}(E_i)` by applying the relevant gamma_{jk}.
   This gives a flat array `vE[i] = gamma(E_i)` — exactly the `aic_dhom_v`
   representation.

3. Certify: `aic_dhom_defect_sweep` checks the multiplicativity defect
   `G_gamma(E_i, E_j)`. The proof of lem_merging says conditions (merging0–2)
   give an O(delta)-homomorphism; (merging3) + lem_alpha give norm equivalence;
   then prop_delta_hominc upgrades to O(delta+eps)-inclusion. This certification
   is done by:
   - `aic_dhom_defect_sweep` -> defect <= O(delta);
   - `aic_dhom_v_sigma_min` -> sigma_min >= 1 - O(delta) (from (merging3)).

**What "gamma_{jk}" are explicitly.** In lem_extension (the only use in
cstar_build), the gamma_{jk} are given by .tex:1406–1409:
```
  gamma_{11}(A_{11}) = v(A_{11}),       v : M_n -> S_P
  gamma_{12}(A_{12}) = U_1(A_{12}),     U_1 : C^n -> S_{P,Q} (a linear map)
  gamma_{21}(A_{21}) = (U_1(A_{21}^dag))^dag
  gamma_{22}(A_{22}) = A_{22} * Qtilde, Qtilde in S_Q (scalar mult)
```
These are EXPLICIT and computable from v and U_1. Building the merged
`aic_dhom_v` from these requires applying each gamma_{jk} to each
basis element of M_{n+1}.

**API status.** No existing function computes the merged gamma from four
sub-maps. This is NEW CODE required for `cstar_build`. However it is
pure array assembly (apply each gamma_{jk} to basis elements) with no
new FLINT/LAPACK primitives.

**eta=0 oracle.** For exact P_j in M_n: gamma is an exact isomorphism,
defect = 0.

**Teeth.** Mutating the compressed product `.` -> plain `*` in
`aic_corner_cdot` (used inside `aic_dhom_defect_sweep`) is RED at eta>0.
Mutating gamma_{12} -> 0 (suppressing the off-diagonal block) is RED
(sigma_min -> 0 for any X_{12} != 0).

### 4.4 lem_extension (.tex:1378–1412)

This is the SUBSTANTIVE lemma. The proof has 6 steps; we map each to API calls.

**Setup.** Inputs: A (eps-C* algebra), P and Q (delta-projections, ||P+Q-I|| <=
delta), v: M_n -> S_P (delta-isomorphism), dim S_Q = 1, S_{P,Q} != 0.
Output: v_+: M_{n+1} -> A (O(delta+eps)-isomorphism).

#### Step 1: dim S_{P,Q} = n

**(.tex:1383.)** The images P_j = v(E_{jj}) are equivalent 1d O(delta+eps)-
projections. Each dim S_{P_j, Q} in {0, 1}; S_{P,Q} != 0 forces each = 1;
lem_add_dim gives dim S_{P,Q} = n * 1 = n.

**API calls:**
```c
  // For j = 0..n-1:
  slong dj = aic_corner_dim_S(A, P_j, Q, prec);
  // ASSERT dj == 1 (fail loud if any is 0 -- the S_{P,Q}!=0 hypothesis)
  // The sum is n * 1 = n.
  slong dPQ = aic_corner_dim_S(A, P, Q, prec);
  // ASSERT dPQ == n (lem_add_dim conclusion)
```

Or use `aic_corner_alpha_dims` with P_parts = {P_j} and Q_parts = {Q}
(p=n, q=1) to get the sum and the total simultaneously.

**New code needed:** the loop over j to check each dim S_{P_j,Q} = 1.
No new API function is required.

#### Step 2: define the Ha-maps h_{jk}

**(.tex:1388–1390.)** Set P_1 = P, P_2 = Q:
```
  h_{jk} = Ha^Q_{P_j, P_k} : S_{P_j,P_k} -> B(S_{P_k,Q}, S_{P_j,Q})
  h_{11} = Ha^Q_{P,P} : S_P -> B(S_{P,Q})  (an O(delta+eps)-homomorphism)
  h_{12} = Ha^Q_{P,Q} : S_{P,Q} -> S_{P,Q} (O(delta+eps)-close to identity)
  h_{21} = Ha^Q_{Q,P} : S_{Q,P} -> S_{Q,P} (O(delta+eps)-close to identity)
  h_{22} = Ha^Q_{Q,Q} : S_Q -> B(S_Q)      (O(delta+eps)-close to identity)
```

**API call:**
```c
  // h_{11} on Z in S_P: aic_corner_ha(Ha11, A, Z, P, P, Q, prec)
  // Ha11 is a d_PQ x d_PQ matrix (operator on S_{P,Q}).
```

But h_{11} is a MAP from S_P to B(S_{P,Q}), i.e., for EACH Z in S_P it
returns a d_PQ x d_PQ matrix. For the lem_approx step, we need the composition
h_{11} o v : M_n -> B(S_{P,Q}).

To apply lem_approx to h_{11}v, we must represent h_{11}v as an
`aic_dhom_v` from M_n to the eps-C* algebra B(S_{P,Q}).

**THE CODOMAIN eps-C* algebra B(S_{P,Q}).** This is the space of n x n
bounded linear maps from S_{P,Q} to itself. As a finite-dim C* algebra of
dimension n^2 (dim S_{P,Q} = n), it is ISOMORPHIC to M_n (this is what the
proof uses: both M_n and B(S_{P,Q}) are n^2-dimensional, so any injective
homomorphism from M_n is an isomorphism). In the algorithm, B(S_{P,Q}) is
represented as M_{d_PQ} (d_PQ = n matrices in the extracted basis {C_l}
of S_{P,Q}). Its product is PLAIN matrix multiplication of the d_PQ x d_PQ
coordinate matrices. Its norm is induced by the Euclidean norm on S_{P,Q}
(the lem_PQ_Hilb inner product, which differs from the operator norm by a
1+O(delta+eps) factor, .tex:1146).

**Representing B(S_{P,Q}) as an aic_ecstar.** This is required to call
`aic_dhom_approx` (which expects codomain `aic_ecstar *A`). Options:

Option A: Build B(S_{P,Q}) as a GENUINE C* algebra aic_ecstar (using
aic_dhom_B for the internal product, not aic_ecstar_star). Since B(S_{P,Q})
is a genuine matrix algebra, its star is exact: Phi = id (the identity
superoperator), so the aic_ecstar's `star_phi` is `NULL` and `phi` is a
Kraus map with a single Kraus operator I_{d_PQ} (the identity on C^{d_PQ}).
The basis is the standard matrix-unit basis E_{lm} (d_PQ^2 elements, each
a d_PQ x d_PQ matrix with a single 1 at (l,m)).

Option B: Do NOT use aic_ecstar for B(S_{P,Q}) at all; instead, represent
h_{11}v directly as a flat linear map from M_n (coordinate basis {E_{lm}})
to d_PQ x d_PQ matrices, and run a custom Newton iteration that treats
the codomain as plain matrix multiplication.

**Recommended route: Option A.** Build B(S_{P,Q}) as a genuine aic_ecstar
with identity Kraus (exact star = plain matmul on d_PQ x d_PQ matrices).
Then `aic_dhom_approx` works unchanged.

**MISSING API:** There is no existing function that wraps a plain d x d
matrix algebra as an `aic_ecstar`. This requires NEW CODE:
`aic_ecstar_init_matrix_algebra(out, d, prec)` — builds an aic_ecstar
whose n = d, basis = the d^2 matrix units E_{lm} (Frobenius-orthonormal:
||E_{lm}||_F = 1), and Kraus phi = a single Kraus operator I_d (so the
star is X*Y = I_d (XY) I_d = XY). This is a SMALL new function (~20 LOC).

#### Step 3: approximate h_{11}v by mu_{11} via lem_approx

**(.tex:1391.)** h_{11}v: M_n -> B(S_{P,Q}) is an O(delta+eps)-homomorphism.
lem_approx (already implemented as aic_dhom_approx) gives an exact
homomorphism mu_{11} with O(delta+eps) accuracy.

**API calls:**
```c
  // 1. Build the h11v representation as aic_dhom_v:
  //    - B_domain: aic_dhom_B with block dims [n] (single M_n block).
  //    - A_codomain: the B(S_{P,Q}) aic_ecstar (genuine, d_PQ x d_PQ).
  //    - vE[mu(0,l,m)] = h_{11}(v(E_{lm}))  for l,m = 0..n-1.
  //    Each h_{11}(v(E_{lm})) is a d_PQ x d_PQ matrix:
  //    call aic_corner_ha(Ha, A_parent, v_E_lm, P, P, Q, prec) to get the
  //    d_PQ x d_PQ matrix of Ha^Q_{P,P}(v(E_{lm})) in the {C_l} basis.
  //
  // 2. Run lem_approx:
  //    aic_dhom_approx(&h11v_dhom, eps_target, tol_abs, unit_tol, max_steps, prec, &st)
  //    Result: mu11_dhom, an exact hom M_n -> B(S_{P,Q}).
```

**Note on the aic_corner_ha interface.** The current `aic_corner_ha` signature:
```c
void aic_corner_ha(acb_mat_t Ha, const aic_ecstar *A, const acb_mat_t Z,
                   const acb_mat_t P, const acb_mat_t R, const acb_mat_t Q, prec)
```
returns the d_PQ x d_RQ matrix of Ha^Q_{P,R}(Z). For h_{11}(Z) = Ha^Q_{P,P}(Z),
call with P=P, R=P, Q=Q; the result is d_PQ x d_PQ. This is ALREADY the
right interface; no change needed.

#### Step 4: U_1 from the polar/SVD of mu_{11}

**(.tex:1404.)** mu_{11}: M_n -> B(S_{P,Q}) is an isomorphism of finite-dim
C* algebras. Every such isomorphism has the form mu_{11}(A) = U_1 A U_1^dag
for a unitary U_1: C^n -> S_{P,Q} (.tex:1404). Constructively: mu_{11} is
represented by its values on the E_{lm} basis; the map E_{lm} |->
U_1 e_l e_m^* U_1^dag says that column l of U_1 is the first column of
mu_{11}(E_{l0}) (normalized). More precisely:

```
  mu_{11}(E_{lm}) = U_1 E_{lm} U_1^dag = (U_1 col_l)(col_m^*)
```

So the columns of U_1 are the n x 1 column vectors u_l such that
mu_{11}(E_{l0}) = u_l * u_0^dag. To extract U_1:

1. From mu_{11}.vE[mu(0,0,0)] (= mu_{11}(E_{00}) in B(S_{P,Q})),
   extract the first column: this gives u_0 (up to phase).
2. From mu_{11}.vE[mu(0,l,0)], extract the l-th column: gives u_l = c_l u_0
   for some scalar c_l (from the E_{l0} = E_{l0} E_{00} factorization).
3. Normalize to get U_1 as a d_PQ x n matrix.

Alternatively: mu_{11}(E_{00}) is a d_PQ x d_PQ rank-1 matrix (since E_{00}
is a rank-1 projector). Its SVD gives U_1 col_0 as the left singular vector.
Then mu_{11}(E_{l0}) = u_l u_0^dag gives u_l = mu_{11}(E_{l0}) u_0 / ||u_0||^2.

**API:** `aic_latd_svd` (LAPACK zgesvd) is available for the SVD.

**U_1 as a linear map from C^n to S_{P,Q}.** U_1 is a d_PQ x n matrix
in the {C_l} coordinate basis of S_{P,Q}. As an n x n operator in M_n,
U_1(e_l) is the operator `sum_m U_1[m,l] C_m`. The map gamma_{12}(A_{12}) =
U_1(A_{12}) for a column vector A_{12} in C^n is then an n x n operator
in S_{P,Q} (a linear combination of the C_m). Concretely:
```
  gamma_{12}(A_{12}) = sum_{l=0..n-1} A_{12}[l] * (sum_m U_1[m,l] C_m)
                     = sum_{l,m} U_1[m,l] A_{12}[l] C_m.
```

#### Step 5: define the gamma_{jk} maps (.tex:1406–1409)

The four maps are:
```
  gamma_{11}(A_{11}) = v(A_{11}),            A_{11} in M_n (n x n)
  gamma_{12}(A_{12}) = U_1(A_{12}),          A_{12} in M_{n,1} = C^n (column)
  gamma_{21}(A_{21}) = (U_1(A_{21}^dag))^dag, A_{21} in M_{1,n} (row)
  gamma_{22}(A_{22}) = A_{22} * Qtilde,       A_{22} in M_1 = C (scalar)
```

Applied to the matrix-unit basis of M_{n+1}:

- E_{l,m} for l,m in [1,n]: block (1,1), gamma_{11}(E_{lm}) = v_E[mu(l,m)].
- E_{l,0} for l in [1,n]:   block (1,2), gamma_{12}(e_l) = U_1 col_l
  (the l-th column of U_1 in the {C_m} representation).
- E_{0,m} for m in [1,n]:   block (2,1), gamma_{21}(e_m^dag) = (U_1 col_m)^dag.
- E_{0,0}: block (2,2), gamma_{22}(1) = Qtilde = aic_corner_Ptilde(A, Q).

**Building the merged aic_dhom_v for v_+:**

1. `B_plus: aic_dhom_B` with block dims [n+1] (single M_{n+1} block).
2. `aic_dhom_v_init(&vplus, &B_plus, A_parent)`.
3. For each matrix unit E_{lm} of M_{n+1} (l,m in [0,n]):
   - Compute `vplus.vE[mu(0,l,m)]` from the appropriate gamma_{jk}.
4. Call `aic_dhom_defect_sweep` (the merging conditions imply O(delta+eps)
   defect) and `aic_dhom_v_sigma_min` to verify inclusion lower bound.
5. Optionally call `aic_errreduce` to reduce to c_0*eps inclusion.

#### Step 6: verify merging conditions and assemble v_+ via lem_merging

**(.tex:1411.)** The proof argues that ||h_{jk} gamma_{jk} - mu_{jk}|| <=
O(delta+eps) implies conditions (merging0–3) hold with delta' = O(delta+eps).
The constructive computation verifies these conditions numerically:

- (merging3): check ||gamma_{jk}(E)|| in [(1-delta), (1+delta)] ||E|| for each
  basis element — pure norm computation.
- (merging1): compute G_gamma(E_i, E_j) for all matrix-unit pairs — uses
  `aic_dhom_defect_sweep` on v_+.
- (merging0): check ||v_+(E_{ji})^dag - v_+(E_{ij}^dag)|| — involution compatibility.
- (merging2): check ||v_+(E_{00}) - Qtilde|| and ||v_+(I_{M_n block}) - Ptilde||.

**Highest-risk convention in lem_extension: the Ha-map index bookkeeping.**
The aic_corner.h header explicitly flags this:
```
  THE Co-INDEX BOOKKEEPING (the module's highest convention-risk):
    Y^dag . Z    : S_{Q,P} x S_{P,R} -> S_{Q,R}   via Co_{Q,R}
    (Y^dag.Z).X  : S_{Q,R} x S_{R,Q} -> S_Q       via Co_Q
    Z . X        : S_{P,R} x S_{R,Q} -> S_{P,Q}   via Co_{P,Q}
    Y^dag.(Z.X)  : S_{Q,P} x S_{P,Q} -> S_Q       via Co_Q
```
For lem_extension, we need h_{jk} = Ha^Q_{P_j, P_k} with P_1=P, P_2=Q:
- h_{11} = Ha^Q_{P,P}: call `aic_corner_ha(Ha, A, Z, P, P, Q, prec)`.
- h_{12} = Ha^Q_{P,Q}: call `aic_corner_ha(Ha, A, Z, P, Q, Q, prec)`.
- h_{21} = Ha^Q_{Q,P}: call `aic_corner_ha(Ha, A, Z, Q, P, Q, prec)`.
- h_{22} = Ha^Q_{Q,Q}: call `aic_corner_ha(Ha, A, Z, Q, Q, Q, prec)`.
The index convention (P,R,Q) in `aic_corner_ha` maps to (P_j, P_k, Q) in
lem_extension notation, where j=row, k=column. VERIFY this against the
actual aic_corner_ha signature before writing the code.

**Dimension of the B(S_{P,Q}) codomain.** d_PQ = n at this point (Step 1).
So B(S_{P,Q}) is represented as M_n (d_PQ x d_PQ = n x n matrices in the
{C_l} basis of S_{P,Q}). The lem_approx call is from M_n (genuine C*) to
M_n (genuine C*, exact star), so the Newton iteration on B(S_{P,Q}) involves
the Pauli diagonal of a single M_n block — this has n^2 terms and cost O(n^4)
per Newton step. For n up to ~10 this is acceptable.

---

## 5. The master loop / proof of th_main (.tex:1414–1444)

### 5.1 Stage 1 — Maximal commutative skeleton

**The "maximum dimensionality" phrasing (.tex:1417).** The proof asserts
existence of a maximum-dimension c_0*eps-inclusion of a commutative C* algebra
into A. The constructive substitute (GREEDY LOOP) is:

1. Start with v_comm: M_1 -> A (the trivial inclusion of the unit).
   Actually, start with the nontrivial projection finder:

```
WHILE any P_j has dim S_{P_j} > 1:
  find m such that dim S_{P_m} > 1
  build the S_{P_m} wrapper (aic_ecstar with compressed star)
  P' = aic_projection_nontrivial(S_{P_m}_wrapper)    [nontrivial projection in S_{P_m}]
  P'' = Ptilde_m - P'                                [the complement]
  merge: cor_merge_sum({P_1,...,P_{m-1},P',P'',P_{m+1},...}) + cor_improvement
  update the projection list
```

The CONTRADICTION ARGUMENT of .tex:1417–1426 proves that this loop
TERMINATES: at termination, all P_j have dim S_{P_j} = 1. The loop runs
at most dim A - m steps (each step increases the number of projections
by 1, and we start with m projections and can have at most dim A).

**Verifying the contradiction argument is constructive (.tex:1419–1426).**
The argument uses `lem_nontriv_projection` (available as
`aic_projection_nontrivial`) and `cor_merge_sum` + `cor_improvement`.
Both are implemented and tested. The contradiction drives the algorithm
(maximum dimensionality = loop termination). The ALGORITHM is the loop;
the contradiction is the termination proof.

**Error maintained.** After each cor_improvement step, v_comm is a c_0*eps-
inclusion. The eps' = O(eps) bound (.tex:1428) accounts for the compressed
algebras S_P being O(eps')-C* with eps' >= c_0*eps.

### 5.2 Stage 2 — Extend each class to a full matrix block

**Inductive step (.tex:1435–1441).** For each equivalence class C = {1,...,s}:

Base case r=1: v_1: M_1 -> A_1 = S_{P_1}. The obvious inclusion maps the
scalar 1 to Ptilde_1 = Co_{P_1}(P_1). No lem_extension needed.

Inductive step r -> r+1:
1. Form A_r = S_{P_{[1,r]}} as an aic_ecstar wrapper (S_P wrapper described
   in §1).
2. The compression (.tex:1437):
```
  v_{r-1}': X |-> Co_{P_{[1,r]}}(v_{r-1}(X))
```
   This gives v_{r-1}': M_{r-1} -> A_{r-1}' where A_{r-1}' = S_{P_{[1,r-1]}}
   seen inside S_{P_{[1,r]}}.

   **API:** For each v_{r-1}.vE[i], compute
   `aic_corner_apply(A_r, Co_{P_{[1,r]}}, v_{r-1}.vE[i], prec)`.
   This gives a new `aic_dhom_v` v_{r-1}' over the S_{P_{[1,r]}} wrapper.

3. Projections inside A_r: P_{[1,r-1]}' = Co_{P_{[1,r]}}(P_{[1,r-1]}),
   P_r' = Co_{P_{[1,r]}}(P_r). Use `aic_corner_apply`.

4. Verify S_{P_{[1,r-1]}', P_r'} != 0 (lem_add_dim: dim = r-1 != 0 since
   j ~ r for all j in C).

5. Call lem_extension (the 6-step procedure of §4.4):
   v_{r-1}^+ : M_r -> A_r = S_{P_{[1,r]}}.

6. Call `aic_errreduce` on v_{r-1}^+ to produce v_r: c_0*eps'-isomorphism
   M_r -> A_r.

**Error bookkeeping.** The error after each step:
- After lem_extension: O(eps') from the O(delta+eps) of v_{r-1}' (.tex:1441).
- After cor_improvement: c_0 * eps' (.tex:1441 "Corollary cor_improvement to
  replace v_{r-1}^+ with a c_0*eps'-isomorphism v_r").
- eps' = O(eps), set at the start of Stage 2 (.tex:1428: "eps' >= c_0*eps,
  and moreover, S_P for any c_0*eps-projection P is an eps'-C* algebra").

The c_0 is FIXED ONCE from `aic_errreduce` (the measured constant from the
first call; the paper says "there exist constants" — we use the measured c_0
throughout). As long as eps' = O(eps) and c_0 is bounded, the final
isomorphism satisfies the paper's bound. The universality canary (c_0 not
growing with dim) was verified for the errreduce module (FINDINGS D2).

### 5.3 Stage 3 — Merge blocks

**(.tex:1443.)** After Stage 2 we have c_0*eps'-isomorphisms v_C: M_{|C|} -> S_{P_C}
for all equivalence classes C. Stage 3 merges them:

```
FOR each pair of classes C, D (C != D):
  ASSERT aic_corner_dim_S(A, P_C, P_D) == 0   [lem_add_dim, by equivalence classes]
  v_merged = cor_merge_sum(v_C, v_D)          [the concat construction of §4.2]
  v_merged = aic_errreduce(v_merged)           [back to c_0*eps']
```

After m-1 merges (m = number of classes), the result is v: B -> A where
B = (+)_C M_{|C|}, a c_0*eps'-isomorphism.

**Error accumulation (.tex:1443).** Each merge introduces O(eps') error;
cor_improvement immediately reduces to c_0*eps'. Since c_0 is fixed and
each reduction applies the same cor_improvement, the final error is c_0*eps'
with the SAME c_0 throughout. The O() in the paper's "O(eps)" absorbs
the c_0 factor (.tex:1415: "construct a c_0*eps-isomorphism").

**The final B.** B = (+)_C M_{|C|}. Its block dims are the equivalence class
sizes. The `aic_dhom_B` for the final v has `num_blocks = m` (the number of
classes) and `d[C] = |C|` for each class C.

---

## 6. Open questions / escalations

### G1. The S_P unit in nontriviality checks (CORRECTION NEEDED)

As identified in §2.3: `aic_projection_nontrivial` checks `||1_n - P||_op >= 0.3`
(ambient). For the S_P wrapper, the correct check is `||Ptilde - P||_op >= 0.3`
where Ptilde = aic_corner_Ptilde(A_parent, P). The code must either:
(a) accept a "unit" argument to `aic_projection_nontrivial`, or
(b) perform the nontriviality check in cstar_build after calling the projector.

Route (b) is preferred (keeps the projection module clean; the cstar_build
loop already computes Ptilde from aic_corner_Ptilde).

### G2. The B(S_{P,Q}) aic_ecstar wrapper (MISSING API — small)

A new function `aic_ecstar_init_matrix_algebra(out, d, prec)` is needed to
represent a d x d matrix algebra as an `aic_ecstar` with exact star (plain
matmul, represented as a single-Kraus-operator I_d UCP map). ~20 LOC.

**Alternatively** (and preferable if simpler): represent B(S_{P,Q}) as an
`aic_dhom_B` with a single M_{d_PQ} block, and pass it directly to a modified
lem_approx that accepts an `aic_dhom_B` for the codomain rather than an
`aic_ecstar`. This avoids creating the wrapper but requires extending
`aic_dhom_approx`'s interface. Either route needs ~20–40 LOC of new code.

### G3. No polar decomposition helper is exposed

`aic_latd_svd` computes the full SVD (LAPACK zgesvd) but does not return
the polar factor U = U_left U_right^dag directly. For Step 4 of lem_extension,
the polar factor of mu_{11} is needed. This requires extracting U_1 from
the SVD as described in §4.4 (Step 4). ~10 LOC of new code after the SVD call.

### G4. Compressing v_{r-1} into A_r: the S_{P_{[1,r]}} wrapper's Co

The compression map `.tex:1437` applies Co_{P_{[1,r]}} to each vE[i]. This
requires `aic_corner_apply(A_parent, Co_{P_{[1,r]}}, vE[i], prec)` for each i.
The Co matrix must be prebuilt: `aic_corner_Co(A_parent, P_{[1,r]}, P_{[1,r]}, prec)`.
No missing API, but the cstar_build loop must maintain the running Co_{P_{[1,r]}}
across inductive steps (rebuild it each time r increments).

### G5. The analytic c_0 is still unmeasured (FINDINGS D2)

`aic_errreduce` returns the MEASURED c_0 per instance. For the master loop,
c_0 must be fixed ONCE at the start (.tex:1415: "Let c_0 be the constant from
cor_improvement"). The implementation uses the first `aic_errreduce` call's
measured c_0 as the nominal constant, and verifies at each subsequent call
that the measured c_0 does not exceed it by more than a small factor.
This is a heuristic; the analytic c_0 (bead aic-1bc) would make it rigorous.

### G6. The Ha-map evaluation cost in lem_extension

`aic_corner_ha` is called for each basis element E_{lm} of M_n to build the
h_{11}v `aic_dhom_v`. For n^2 basis elements, each `aic_corner_ha` call
solves a d_PQ x d_PQ Gram system. Total cost: O(n^2 * d_PQ^3) = O(n^5)
(since d_PQ = n). For n up to ~6 this is fast; for n > 10 this could be
slow. Profile and flag if the Stage 2 inductive loop is bottlenecked here.

### G7. lem_approx termination in the B(S_{P,Q}) setting

`aic_dhom_approx` terminates when the basis-sweep defect <= max(tol_abs, eps_target).
For the B(S_{P,Q}) codomain (a genuine M_n), eps_target should be set to
the current eps' (the eps of the S_{P_{[1,r]}} wrapper's star, which is
O(delta+eps) from the v_{r-1}' compression). The termination criterion must
use the *codomain* algebra's eps (the eps of B(S_{P,Q}) as an eps-C*
algebra — which is 0, since B(S_{P,Q}) is a GENUINE C* algebra). So
eps_target = 0 in principle, but in finite precision the iteration stops
when the defect is < tol_abs (machine floor). This is the same as the
eta=0 oracle path (exact hom approximation). **Sound; no new code.**

### G8. dim S_{P_j, Q} = 1 assertion at the start of Stage 2

The proof (.tex:1435) does not explicitly verify that after Stage 1, each
dim S_{P_j, P_r} = 1 for j in [1,r-1] and r the current "next" index.
This follows from the equivalence class definition (j ~ r iff dim S_{P_j,P_r}=1).
But the equivalence-class partition was computed from `aic_corner_equiv_1d`
on the Stage-1 projections. In the implementation:

- Compute the equivalence classes at the END of Stage 1 using
  `aic_corner_equiv_1d(A, P_j, P_k, prec)` for all pairs (j,k).
- At the START of Stage 2 induction step r: ASSERT
  `aic_corner_dim_S(A, P_j, P_r) == 1` for all j in [1,r-1].

This assertion is a fail-loud guard (Rule 4), and verifies the hypothesis
of lem_extension (S_{P_{[1,r-1]}, P_r} != 0, from lem_add_dim).

---

## 7. Increment plan

### Proposed order (I1–I5) and critique

**I1: S_P-as-aic_ecstar wrapper + B(S_{P,Q}) init**

Files: `src/aic_cstar_sp_wrap.c` (the S_P wrapper constructor/destructor,
~60 LOC) and `src/aic_cstar_bsq.c` (B(S_{P,Q}) as genuine aic_ecstar, ~30 LOC).

Tests: eta=0 oracle (exact idempotent, S_P wrapper star = exact compressed
product = 0 defect); oblique channel (star via compressed product gives
nonzero defect ~O(eta)); projection works on S_P wrapper (nontrivial proj
found in S_P). Tier: Core (hostile review required). Dependency: none (reuses
existing corner + ecstar).

**CORRECTION:** Also implement the nontriviality check fix (§G1) here:
either add a `Ptilde` argument to `aic_projection_nontrivial` or check it
in the wrapper test.

**I2: lem_add_dim + cor_merge_sum**

Files: `src/aic_cstar_merge.c` (functions `aic_cstar_lem_add_dim` and
`aic_cstar_merge_sum`, ~80 LOC total). These are the SIMPLER of the
assembly lemmas.

Tests: eta=0 oracle (exact P, Q orthogonal in M_n: dim S_{P,Q} = sum dim
S_{P_j,Q_k} exactly); oblique channel (verify dim count + inclusion defect
O(eta)); S_{P_C,P_D}=0 assertion on distinct classes.

Tier: Core. Dependency: I1 (for the S_P wrapper tests).

**I3: lem_merging**

Files: `src/aic_cstar_merging.c` (~80 LOC). Takes four gamma_{jk} maps (as
arrays of n x n operators indexed by matrix units) and assembles the merged
aic_dhom_v.

Tests: eta=0 oracle (M_2 decomposed as two M_1 blocks: gamma_{11}=v, gamma_{12}=0,
etc. -> exact isomorphism); oblique fixture (merging conditions hold with
O(eta) defect); mutation: gamma_{12} -> 0 when S_{Pi_1,Pi_2} != 0 -> RED
(sigma_min drops).

Tier: Core. Dependency: I2 (for cor_merge_sum verification pattern).

**I4: lem_extension**

Files: `src/aic_cstar_extension.c` (~120 LOC). The substantive lemma.
Calls: aic_corner_ha (h_{11}), aic_dhom_approx (lem_approx), aic_latd_svd
(U_1 extraction), lem_merging (I3).

Tests: Base case n=1: v: M_1 -> S_P (scalar), Q 1d, S_{P,Q} != 0 -> v_+:
M_2 -> A (exact M_2 if A=M_2, eta=0). Inductive case n=2: M_2 extended
to M_3. Ha-map bookkeeping cross-check: h_{11} is the O(delta+eps)-homomorphism;
mu_{11} has defect < tol_abs on the exact C* codomain.

Critical teeth: ha-map index permutation (swap P and R in aic_corner_ha) ->
RED (h_{11} fails to be a homomorphism). Star-vs-plain in the h_{11}v
construction -> RED at eta>0. U_1 phase ambiguity (two distinct global phases
of U_1 give the same mu_{11} but different gamma_{12}/gamma_{21}) -> GREEN
(the merging conditions depend on ||h_{jk} gamma_{jk} - mu_{jk}||, not on
U_1 directly; both phases pass).

Tier: Core (hostile review required — highest-risk increment). Dependency: I3.

**I5: Master loop (proof of th_main)**

Files: `src/aic_cstar_build.c` (~150 LOC) + `include/aic_cstar.h` (~80 LOC).
The main loop: Stage 1 greedy projection, Stage 2 inductive extension, Stage 3
merge. Input: `aic_ecstar *A`. Output: `aic_dhom_B *B` and `aic_dhom_v *v`
with v being a c_0*eps-isomorphism.

Tests:
- eta=0 oracle: exact idempotent channel -> B = (+)_C M_{dim C} isomorphic
  to A with zero defect and the exact th_idemp_structure block sizes.
- Universality canary: sweep dim A from 4 to 16, verify c_0 does not grow
  (the headline test, inheriting from errreduce's T4 canary).
- Oblique fixture: make_compress_idemp / make_mixconj at several eta values.
- B structure: verify B is a direct sum of matrix algebras (num_blocks x
  d[l]^2 sum = dim A as a consistency check).

Critical teeth: Star-vs-plain in any product inside the loop -> RED at eta>0.
Projection of S_P using Frobenius Pi_A instead of Co_P(Phi_tilde(.)) ->
RED (the star defect of the found P stays O(1) instead of O(eta), FINDINGS §C3).

Tier: Core (hostile review required). Dependency: I4.

### Dependency critique

The proposed I1-I2-I3-I4-I5 order is correct: each increment is a prerequisite
for the next. I2 (cor_merge_sum) is logically independent of I3 (lem_merging)
in the .tex, but I3 is needed by I4, and I4 consumes lem_merging via the
gamma_{jk} assembly. The B(S_{P,Q}) init can be done in I1 (same session as
the S_P wrapper) since they are related infrastructure.

**I1 is the critical path gater.** If the S_P wrapper's star gives the wrong
oblique projection (FINDINGS §C3 risk), it will appear as a non-O(eta) star
defect on the S_P wrapper itself, and every subsequent increment will fail.
Verify this with explicit arb certification before proceeding to I2.

**I4 is the highest-risk increment** (see §8 below).

---

## 8. Summary

**(a) Verdict on the S_P-wrapper data model + projection integration.**

The S_P wrapper is CORRECT in principle: the star `Co_P(Phi_tilde(XY))` is
the compressed product of S_P (.tex:1080), the basis is the extracted corner
basis, and the unit is Ptilde = Co_P(P). The projection integration is MOSTLY
SOUND: the ambient-1_n projection step gives the correct oblique projection
into S_P (via the star_phi seam). ONE CORRECTION is required: the
nontriviality check `||1_n - P||_op >= 0.3` should be `||Ptilde - P||_op >= 0.3`
for the S_P wrapper. This must be handled in cstar_build (not in the
projection module, to avoid API changes).

**(b) Highest-risk part of lem_extension.**

The highest-risk part is **Step 2–3**: representing h_{11}v as an `aic_dhom_v`
with the B(S_{P,Q}) genuine C* algebra as codomain, and verifying the Ha-map
index bookkeeping. The `aic_corner_ha` signature has three projection arguments
(P, R, Q) that must map to (P_1=P, P_2=Q, Q=Q) correctly for each of the four
h_{jk} maps. The FINDINGS §C2 (star-vs-plain blindness) and the aic_corner.h
"highest convention-risk" warning both point here. The mutation `aic_corner_ha(A,
Z, P, R, Q)` with P and R transposed (giving h_{11} -> h_{11}^T) must be a
RED test.

**(c) Missing existing-API capabilities.**

1. `aic_ecstar_init_matrix_algebra(out, d, prec)` — wrap M_d as a genuine
   aic_ecstar with exact star (identity Kraus). Required by I4. ~20 LOC.
2. A polar-factor extractor after `aic_latd_svd` — needed to get U_1 in
   Step 4 of lem_extension. ~10 LOC, NEW inline in `aic_cstar_extension.c`.
3. The nontriviality-with-Ptilde check — not an aic_projection.h change, but
   a cstar_build-level assert. ~5 LOC.

No other MISSING API capabilities were identified. All FLINT/LAPACK primitives
needed (acb_mat_set, aic_corner_Co, aic_corner_apply, aic_corner_extract,
aic_corner_ha, aic_corner_dim_S, aic_corner_equiv_1d, aic_dhom_approx,
aic_errreduce, aic_latd_svd, aic_projection_nontrivial) are CLOSED + green.

**(d) Recommended increment plan.**

I1 (S_P wrapper + B(S_{P,Q}) init, Core tier), I2 (lem_add_dim + cor_merge_sum,
Core), I3 (lem_merging, Core), I4 (lem_extension, Core, hostile review mandatory),
I5 (master loop, Core, hostile review mandatory). The proposed I1–I5 ordering
is correct. No reordering is beneficial. Each increment should include the
eta=0 oracle test AND an oblique eta>0 fixture with star-vs-plain mutation
proof, following the FINDINGS §C2 discipline.

---

*End of design document.*
