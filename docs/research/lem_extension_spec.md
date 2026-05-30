# Implementation spec: `aic_cstar_extension` — Increment 4 (`lem_extension`)

**Module bead:** `aic-097` (cstar_build), Increment 4
**Tex citation:** `approximate_algebras.tex:1378-1412`
**Depends on:** I1 (`aic_cstar_subalg_build`, `aic_cstar_matrix_algebra`),
               I2 (`aic_cstar_lem_add_dim`, `aic_cstar_merge_sum`),
               I3 (`aic_cstar_lem_merging`)
**Date:** 2026-05-30

---

## 0. Reading and verification record

Files read verbatim:
- `.tex:1052-1189` (§7: S_{P,Q} definition through lem_1d_proj + Ha-map definition
  at `.tex:1146-1160`)
- `.tex:1325-1412` (§9 opening: lem_merging, cor_merge_sum, lem_add_dim,
  lem_extension statement+proof)
- `include/aic_corner.h` (full; `aic_corner_ha` signature and Co-index bookkeeping)
- `include/aic_dhom.h` (full; `aic_dhom_approx` parameters)
- `include/aic_cstar.h` (full; `aic_cstar_matrix_algebra`, `aic_cstar_lem_merging`
  including routing convention)
- `include/aic_latd.h` (full; `aic_latd_svd` convention)
- `include/aic_ecstar.h` (full; struct layout)
- `docs/research/cstar_build_design.md` §4.4 and §6 (the proposed 6-step plan and
  open questions G1-G8)
- `paper/FINDINGS.md` (full)
- `tests/test_cstar_merging.c` (full; I3 fixture and convention patterns)

---

## 1. Faithfulness check: design §4.4 vs `.tex:1378-1412`

### 1.1 The `.tex` proof structure

The proof of `lem_extension` (`.tex:1382-1411`) proceeds as follows. Quoting the
key moves with line numbers:

**`.tex:1383`** — dimensionality of S_{P,Q}:
> "Since v: M_n -> S_P is a delta-isomorphism, the elements P_j = v(E_{jj}) are
> one-dimensional delta-projections in S_P, which are equivalent to each other.
> Hence P_1,...,P_n are equivalent one-dimensional O(delta+eps)-projections in A,
> implying that either dim S_{P_j,Q} = 0 for all j or dim S_{P_j,Q} = 1 for all j.
> Using Lemma lem_add_dim, we conclude that dim S_{P,Q} is either 0 or n.
> The first possibility contradicts the hypothesis; thus dim S_{P,Q} = n."

**`.tex:1386-1389`** — definition of h_{jk}:
```
  P_1 = P,  P_2 = Q,
  h_{jk} = Ha^Q_{P_j,P_k} : S_{P_j,P_k} -> B(S_{P_k,Q}, S_{P_j,Q}).
```
> "h_{11}: S_P -> B(S_{P,Q}) is an O(delta+eps)-homomorphism, and h_{12}, h_{22},
> h_{21} are O(delta+eps)-close to identity maps."

**`.tex:1391`** — lem_approx applied to h_{11}v:
> "The composition h_{11}v: M_n -> B(S_{P,Q}) is an O(delta+eps)-homomorphism of
> C* algebras. By Lemma lem_approx, it can be approximated by a homomorphism
> mu_{11}: M_n -> B(S_{P,Q}) with O(delta+eps) accuracy. Since M_n has no
> nontrivial ideals, mu_{11} is injective. It is actually an isomorphism because
> both M_n and B(S_{P,Q}) are n^2-dimensional."

**`.tex:1404`** — U_1 extraction:
> "Being an isomorphism of finite-dimensional C* algebras, mu_{11} has the form
> mu_{11}(A) = U_1 A U_1^dag for some unitary map U_1: C^n -> S_{P,Q}."

Also introduces U_2: c |-> c Qtilde from C to S_{Q,Q}, and the four mu_{jk}(A) =
U_j A U_k^dag defining an isomorphism M_{n+1} -> B(S_{P,Q} (+) S_{Q,Q}).

**`.tex:1405-1410`** — the four gamma_{jk} (verbatim):
```
  gamma_{11}(A_{11}) = v(A_{11}),
  gamma_{12}(A_{12}) = U_1(A_{12}),
  gamma_{21}(A_{21}) = (U_1(A_{21}^dag))^dag,
  gamma_{22}(A_{22}) = A_{22} Qtilde,
```
where A_{11} is n x n, A_{12} in M_{n,1} = C^n (column vector), A_{21} in
M_{1,n} = (C^n)^* (row vector), A_{22} is a complex number.

**`.tex:1411`** — lem_merging applied:
> "It follows that ||h_{jk} gamma_{jk} - mu_{jk}|| <= O(delta+eps). Thus
> equations (merging0)-(merging3) from the hypothesis of Lemma lem_merging (with
> delta replaced by delta' = O(delta+eps)) follow from (merging0h)-(merging3h).
> Applying the lemma, we obtain the required O(delta+eps)-isomorphism v_+: M_{n+1}
> -> A, namely, v_+((A_{11} A_{12} / A_{21} A_{22})) = sum_{jk} gamma_{jk}(A_{jk})."

### 1.2 Faithfulness verdict

The design §4.4's 6-step decomposition (Steps 1-6) maps faithfully to the proof.

**Steps 1-2 are FAITHFUL.** Step 1 (dim S_{P,Q} = n, `.tex:1383`) and Step 2 (the
h_{jk} definitions, `.tex:1386-1389`) match exactly.

**Step 3 is FAITHFUL.** The h_{11}v composition and lem_approx call match
`.tex:1391`. The design correctly identifies B(S_{P,Q}) as M_n (a genuine C*
algebra) and that mu_{11} is an isomorphism (both sides are n^2-dimensional).

**Step 4 is FAITHFUL.** The U_1 extraction from mu_{11} matches `.tex:1404`.
The design's recipe (SVD of mu_{11}(E_{00}) for col_0, then mu_{11}(E_{l0}) u_0 /
||u_0||^2 for col_l) correctly realizes `mu_{11}(A) = U_1 A U_1^dag`. See §4
below for the exact math.

**Step 5 is FAITHFUL.** The four gamma_{jk} at `.tex:1405-1410` match the design's
§4.4 Step 5 formulas. The design correctly notes gamma_{12}(A_{12}) = U_1(A_{12})
for A_{12} a column vector.

**Step 6 is FAITHFUL.** The design correctly routes through `aic_cstar_lem_merging`
with `two_block=0` and n1=n, n2=1 to assemble v_+.

**ONE DISCREPANCY FOUND** (FINDING F-I4-1 below): the design's Step 6 discusses
verifying the merging conditions as a separate step, but the proof's logic is that
the conditions *follow automatically* from ||h_{jk} gamma_{jk} - mu_{jk}|| <=
O(delta+eps) (`.tex:1411`, which uses the h_{jk}^{-1} bounds at `.tex:1392-1401`).
The constructive content is: the defect_sweep of the merged v will report O(delta+
eps) (the merging1 violation goes through the h-map invertibility bound), and
sigma_min reports >= 1 - O(delta+eps) (from merging3 applied to the gamma_{jk}
norms). No separate per-block condition sweep is required before calling
lem_merging; the output guards (mult_def, sigma_min, merge_cond_max) from
`aic_cstar_lem_merging` already REPORT the joint result. This is a simplification,
not an error in the design.

---

## 2. The Ha-map argument order (HIGHEST-RISK item)

### 2.1 The `aic_corner_ha` signature (verbatim from header)

```c
void aic_corner_ha(acb_mat_t Ha, const aic_ecstar *A, const acb_mat_t Z,
                   const acb_mat_t P, const acb_mat_t R, const acb_mat_t Q,
                   slong prec);
```

From the header docstring:
> "The Ha-map Ha^Q_{P,R}(Z) (.tex:1146-1160). [...] Returned as a d_PQ x d_RQ
> matrix `Ha` in the EXTRACTED corner bases {C^{PQ}_l} of S_{P,Q} (rows) and
> {C^{RQ}_m} of S_{R,Q} (columns)."
> "Ha must be init'd d_PQ x d_RQ (d_PQ = dim S_{P,Q}, d_RQ = dim S_{R,Q})."
> "Z is n x n in S_{P,R}."

So the function computes Ha^Q_{P,R}(Z), a linear map S_{R,Q} -> S_{P,Q},
returned as a d_PQ x d_RQ matrix. The argument order `(A, Z, P, R, Q)` maps to
the notation `Ha^Q_{P,R}(Z)`: output is (d_{PQ} x d_{RQ}), operand is Z in
S_{P,R}, lower-left index is P, lower-right index is R.

### 2.2 The four h_{jk} calls

From `.tex:1388`, with P_1 = P and P_2 = Q:

```
h_{jk} = Ha^Q_{P_j, P_k}  : S_{P_j, P_k} -> B(S_{P_k, Q}, S_{P_j, Q})
```

So h_{jk}(Z) is a LINEAR MAP from S_{P_k, Q} to S_{P_j, Q}, returned as a
d_{P_j,Q} x d_{P_k,Q} matrix. Mapping to the `aic_corner_ha(Ha, A, Z, P, R, Q)`
call: P = P_j (row index), R = P_k (column index), the fixed Q = Q.

**The four calls, each for a given Z in S_{P_j, P_k}:**

| Map         | `.tex` notation       | `aic_corner_ha` arguments (P, R, Q) | Output shape    |
|-------------|----------------------|--------------------------------------|-----------------|
| h_{11}(Z)   | Ha^Q_{P,P}(Z)        | `(P, P, Q)`                          | d_PQ x d_PQ     |
| h_{12}(Z)   | Ha^Q_{P,Q}(Z)        | `(P, Q, Q)`                          | d_PQ x d_QQ     |
| h_{21}(Z)   | Ha^Q_{Q,P}(Z)        | `(Q, P, Q)`                          | d_QQ x d_PQ     |
| h_{22}(Z)   | Ha^Q_{Q,Q}(Z)        | `(Q, Q, Q)`                          | d_QQ x d_QQ     |

Where d_PQ = dim S_{P,Q} = n (established in Step 1) and d_QQ = dim S_{Q,Q} = 1
(since Q is a one-dimensional delta-projection by hypothesis). So:
- h_{11}(Z) is an n x n matrix (the homomorphism S_P -> B(S_{P,Q}))
- h_{12}(Z) is an n x 1 column (S_{P,Q} -> S_{P,Q} as a scalar, "~identity")
- h_{21}(Z) is a 1 x n row (S_{Q,P} -> S_{Q,P} as a scalar, "~identity")
- h_{22}(Z) is a 1 x 1 scalar (~identity)

### 2.3 Confirming vs the design doc

The design doc (§4.4 Step 2) states:
```
h_{11} = Ha^Q_{P,P}: call aic_corner_ha(Ha, A, Z, P, P, Q, prec)
h_{12} = Ha^Q_{P,Q}: call aic_corner_ha(Ha, A, Z, P, Q, Q, prec)
h_{21} = Ha^Q_{Q,P}: call aic_corner_ha(Ha, A, Z, Q, P, Q, prec)
h_{22} = Ha^Q_{Q,Q}: call aic_corner_ha(Ha, A, Z, Q, Q, Q, prec)
```

This is **CONFIRMED CORRECT**. The design's `(P, R, Q)` = `(P_j, P_k, Q)` mapping
is exactly what the header documents.

### 2.4 Mutation specification

The REQUIRED mutation that MUST turn the I4 test RED:

**Mutation (a): swap P and R in the h_{11} call.**
Replace `aic_corner_ha(Ha11, A, Z, P, P, Q, prec)` with
`aic_corner_ha(Ha11, A, Z, Q, P, Q, prec)`.

This computes Ha^Q_{Q,P}(Z) instead of Ha^Q_{P,P}(Z). The output is now a
d_QQ x d_PQ = 1 x n matrix (wrong shape: should be n x n). Even if shape-matched
by accident, the homomorphism property of h_{11} fails: Ha^Q_{P,P} is an
O(delta+eps)-homomorphism S_P -> B(S_{P,Q}) (`.tex:1160`), while Ha^Q_{Q,P}
maps S_{Q,P} -> B(S_{P,Q}, S_{Q,Q}) (a different map entirely).

**Predicted symptom:**
- ASSERT in `aic_corner_ha` fires on Ha dimension mismatch (Ha init'd n x n
  but the correct output is 1 x n), **OR**
- if the caller passes a 1 x n Ha: the h_{11}v `aic_dhom_v` is assembled with
  1 x n matrices (wrong), and the subsequent lem_approx iteration either
  (a) asserts a dimension mismatch on the B(S_{P,Q}) codomain's n x n basis,
  or (b) returns a mu_{11} with mult_def >> O(delta+eps) and sigma_min << 1.
- In the test: the Step-3 lem_approx defect sweep over the assembled h11v
  dhom_v fails to converge or the initial defect reads >> 0.1 (O(1)), RED.

---

## 3. mu_{11} construction: `lem_approx` on h_{11}v

### 3.1 Domain and codomain

- **Domain B:** `aic_dhom_B` with block dims `[n]`, a single M_n block.
  Build via `aic_dhom_B_init(&B_domain, &n, 1)`.
- **Codomain A_cod:** `aic_cstar_matrix_algebra` with d = n (the M_n genuine
  C* algebra B(S_{P,Q}) in the {C_l} coordinate basis of S_{P,Q}).
  Star = plain matrix multiplication on n x n matrices. Every eps-C* axiom defect
  is machine-zero.
- **The map h11v:** `aic_dhom_v_init(&h11v, &B_domain, &A_cod)`, then for each
  matrix unit E_{lm} of M_n (linear index mu(0,l,m) = l*n + m):
  ```
  aic_corner_ha(h11v.vE[l*n + m], A_parent, v->vE[l*n + m], P, P, Q, prec)
  ```
  where `v->vE[l*n+m]` is the image v(E_{lm}) in S_P (an n x n operator in the
  parent A). The result Ha is an n x n matrix in A_cod's coordinate basis
  (columns = coords of Ha^Q_{P,P}(v(E_{lm}))(C_j^{RQ}) in the {C_l^{PQ}} basis
  of S_{P,Q}). Since dim S_{P,Q} = n = d_PQ = d_RQ, each Ha is n x n. CORRECT.

### 3.2 `aic_dhom_approx` parameters

The full signature (from `aic_dhom.h`):
```c
void aic_dhom_approx(aic_dhom_v *v, double eps_target, double tol_abs,
                     double unit_tol, slong max_steps, slong prec,
                     aic_dhom_approx_stats *st);
```

Parameter choices for the lem_extension call:
- `eps_target` = **0.0** (the codomain A_cod is a GENUINE C* algebra — machine-zero
  star-defect — so the O(eps) floor is 0; the iteration should reach tol_abs).
  Compare design §G7: "eps_target = 0 in principle [...] iteration stops when
  defect < tol_abs (machine floor)."
- `tol_abs` = **1e-12** at `prec=256` (consistent with the established project
  pattern: `test_cstar_merging.c` uses the `< 1e-12` threshold for eta=0 oracle).
- `unit_tol` = **1e-10** (per-step assert on unit preservation; generous relative
  to tol_abs; any value in (tol_abs, 1e-8) is safe).
- `max_steps` = **30** (the Newton iteration is quadratic; from delta_0 ~ O(eta)
  with eta ~ 0.1, we need log2(log2(0.1/1e-12)) ~ 6 steps; 30 is a safe cap).

Note: h_{11}v is an O(delta+eps)-homomorphism. The defect after lem_approx is
O(delta+eps); on the eta=0 oracle path delta ~ 0 and the iteration reaches machine
zero. The `AIC_ERRREDUCE_EPS_FLOOR = 4*eps` rule (`.tex:1302` comment in dhom
module) does NOT apply here: that floor applies when the codomain is itself an
O(eps)-C* algebra; here A_cod is genuine C*, so the floor is literally 0.

---

## 4. U_1 extraction (polar of mu_{11})

### 4.1 Mathematical basis (`.tex:1404`)

`mu_{11}(A) = U_1 A U_1^dag` for a UNITARY U_1: C^n -> S_{P,Q}.

In the {C_l^{PQ}} coordinate basis of S_{P,Q} (which has n elements since
dim S_{P,Q} = n), U_1 is represented as an n x n unitary matrix.

From the formula mu_{11}(E_{lm}) = U_1 E_{lm} U_1^dag = (U_1 col_l)(col_m U_1^dag):
- E_{00} is a rank-1 projector. mu_{11}(E_{00}) = u_0 u_0^dag (rank-1).
- E_{l0} = E_{l0} E_{00}, so mu_{11}(E_{l0}) = u_l u_0^dag.
- Hence u_l = mu_{11}(E_{l0}) u_0 / <u_0, u_0> = mu_{11}(E_{l0}) u_0 / ||u_0||^2.

### 4.2 Extraction algorithm (~10 LOC)

```
Step (i): mu11_00 = mu11.vE[0]   (the n x n matrix mu_{11}(E_{00}))
Step (ii): compute SVD of mu11_00 via aic_latd_svd
           svals[min(n,n)=n], U_left (n x n), Vt (n x n)
           mu11_00 = U_left * diag(svals) * Vt
           Note: svals[0] ~ 1 (the only nonzero singular value; mu_{11}(E_{00}) is
           rank-1 since E_{00} is a rank-1 projector and mu_{11} is an isomorphism)
           u0 = U_left[:, 0]  (the first LEFT singular vector, a length-n column)

Step (iii): for l = 0..n-1:
           mu11_l0 = mu11.vE[l * n + 0]     (n x n matrix mu_{11}(E_{l0}))
           norm_u0_sq = sum_m |u0[m]|^2 (~1; ASSERT > 0.5)
           U1_col_l = mu11_l0 * u0 / norm_u0_sq
           (result: a length-n column vector, the l-th column of U_1)

Step (iv): Assemble U1 as an n x n matrix: U1[:, l] = U1_col_l.
```

U_1 is a d_PQ x n = n x n unitary matrix IN THE COORDINATE BASIS of S_{P,Q}
(NOT a physical n x n operator in M_n(C); it is an n x n matrix whose columns
are the S_{P,Q}-coordinates of the images U_1(e_l)).

**Converting U_1 to a physical operator map.** The element `U_1(e_l)` as a
physical n x n operator in M_n is:
```
  U1_op[l] = sum_{m=0}^{n-1} U1[m, l] * C_m   (C_m = m-th corner basis operator)
```
where `{C_m}` are the extracted Frobenius-orthonormal operators of S_{P,Q}
(obtained from `aic_corner_extract` on the Co_{P,Q} matrix).

### 4.3 Convention: U_left from LAPACK `aic_latd_svd`

The `aic_latd_svd` convention (header):
> "A(i,j) = sum_k U(i,k) * svals[k] * Vt(k,j). [...] U: row-major m*m, the left
> singular vectors as COLUMNS (jobu='A')."

So `U[i*n + k] = U_{ik}`, and column k = `{U[i*n + k] : i=0..n-1}`. Column 0 is
`u0[i] = U[i*n + 0]`. This is the first LEFT singular vector. The LAPACK
convention is CONFIRMED (left singular vectors as columns of U).

**Note on C_l basis.** The `{C_m}` operators of S_{P,Q} are obtained from
`aic_corner_extract` on `aic_corner_Co(A, P, Q, prec)`. The `aic_corner_ha`
function uses these INTERNALLY (it "builds all Co's and bases internally from
P,R,Q" per the header). The Ha output is already in the {C^{PQ}_l} basis; the
coordinate representation is consistent.

### 4.4 Dimension note

This constructs a d_PQ x n = n x n isometry U_1 (in the S_{P,Q} coordinate basis).
In the coordinate basis it is unitary (n x n) because dim S_{P,Q} = n = dim C^n.
The extraction is well-defined and square.

---

## 5. Block-index convention for M_{n+1} and the four gamma_{jk} arrays

### 5.1 The chosen convention

**CONVENTION**: n1 = n (the M_n block, associated to P / P_1 = P),
               n2 = 1 (the M_1 block, associated to Q / P_2 = Q).

In `aic_cstar_lem_merging` with `two_block=0`, `n1=n`, `n2=1`:
- The global matrix index (l, m) of M_{n+1} = M_{n1+n2} is mapped to a block by:
  `l < n1 AND m < n1`  -> block (1,1), local (a,b) = (l, m)
  `l < n1 AND m >= n1` -> block (1,2), local (a,b) = (l, m-n1) = (l, m-n)
  `l >= n1 AND m < n1` -> block (2,1), local (a,b) = (l-n1, m) = (l-n, m)
  `l >= n1 AND m >= n1`-> block (2,2), local (a,b) = (l-n, m-n) = (0, 0)
- The linear index of v_out->vE[i] = i = l*(n1+n2) + m = l*(n+1) + m.

**The Q/M_1 block is at GLOBAL INDICES l=n..n, m=n..n (last row and column).**
The P/M_n block occupies indices l=0..n-1, m=0..n-1.

This matches the `.tex` embedding M_n ⊂ M_{n+1} (`.tex:1371`: "M_n has standard
basis {E_{jk}: j,k=1,...,n}" extended by the (n+1)-th index for Q). In 0-indexed
code: P occupies indices 0..n-1, Q occupies index n.

**Reconciling with the design doc's "E_{00} for the Q/M_1 block"**: The design
§4.4 Step 5 uses E_{00} for the Q block ("E_{0,0}: block (2,2), gamma_{22}(1) =
Qtilde"). This is the design's *within-block* local index (local index (0,0) for
the scalar 1 in M_1), NOT the global index. The global position of the Q block is
(n, n). The design is internally consistent; it names local block-indices. The
physical `aic_cstar_lem_merging` routing uses GLOBAL indices, and the Q block sits
at global (n, n). CONFIRMED CONSISTENT.

### 5.2 The four gamma_{jk} arrays (exact formulas)

Following the `aic_cstar_lem_merging` header contract:
```
gamma11 : n1*n1 = n*n operators, gamma11[a*n1 + b] = gamma_{11}(E^{(11)}_{ab})
gamma12 : n1*n2 = n*1 operators, gamma12[a*n2 + 0] = gamma_{12}(E^{(12)}_{a,0})
gamma21 : n2*n1 = 1*n operators, gamma21[0*n1 + b] = gamma_{21}(E^{(21)}_{0,b})
gamma22 : n2*n2 = 1*1 operators, gamma22[0*n2 + 0] = gamma_{22}(E^{(22)}_{0,0})
```

From `.tex:1405-1410`:

**gamma11[a*n + b] = gamma_{11}(E_{ab}) = v(E_{ab}) = v->vE[a*n + b]**
(a copy of the delta-isomorphism v's a*n+b image; deep copy via acb_mat_set)

**gamma12[a] = gamma_{12}(e_a) = U_1(e_a) = U_1 col_a as a physical operator**
= sum_{m=0}^{n-1} U1[m, a] * C_m,
where C_m are the physical n x n operators of S_{P,Q}'s corner basis. This is
an n x n operator IN A (the physical image of the a-th standard basis vector
of C^n under U_1, landing in S_{P,Q}).

**gamma21[b] = gamma_{21}(e_b^dag) = (U_1(e_b))^dag = (gamma12[b])^dag**
(acb_mat_conjugate_transpose of gamma12[b])

Note: `.tex:1407-1408` gives gamma_{21}(A_{21}) = (U_1(A_{21}^dag))^dag for
A_{21} a row vector. For A_{21} = e_b^dag (a row with 1 at position b),
A_{21}^dag = e_b (column vector), so gamma_{21}(e_b^dag) = (U_1(e_b))^dag =
(gamma12[b])^dag. Index b runs 0..n-1, stored as gamma21[0*n + b] = gamma21[b].

**gamma22[0] = gamma_{22}(1) = 1 * Qtilde = aic_corner_Ptilde(A_parent, Q)**
(Qtilde = Co_Q(Q) = aic_corner_Ptilde; computed once and stored)

### 5.3 Index correctness summary

The lem_merging routing table (from `aic_cstar.h`):
```
  l <  n, m <  n : block (1,1), (a,b)=(l,m),   -> gamma11[l*n + m]   ✓
  l <  n, m == n : block (1,2), (a,b)=(l, 0),   -> gamma12[l*1 + 0]  ✓
  l == n, m <  n : block (2,1), (a,b)=(0, m),   -> gamma21[0*n + m]  ✓
  l == n, m == n : block (2,2), (a,b)=(0, 0),   -> gamma22[0*1 + 0]  ✓
```
(using n1=n, n2=1 so m >= n1 means m >= n means m == n, and l >= n1 means l == n)

**The global linear index for v_out->vE at row l, col m is l*(n+1)+m** (since
aic_dhom_B_index(B, 0, l, m) = offset[0] + l*(n+1) + m = l*(n+1) + m for a
single block M_{n+1}).

---

## 6. New code inventory (I4)

### 6.1 New functions

**File: `src/aic_cstar_extension.c`** (target ~120 LOC, Rule 10)

```
aic_cstar_lem_extension(B_out, v_out, mult_def, sigma_min,
                        v, P, Q, A_parent, prec)
```
The top-level function. Inputs: `v` (the delta-isomorphism M_n -> S_P as
`aic_dhom_v`), `P`, `Q` (n x n Hermitian delta-projections), `A_parent` (the
parent eps-C* algebra), `prec`. Outputs: `B_out` (M_{n+1} as `aic_dhom_B`),
`v_out` (the assembled v_+ as `aic_dhom_v`), `mult_def` and `sigma_min` (the
certified output guards from `aic_cstar_lem_merging`).

Internal helpers (static, same file or aic_cstar_extension_internal.h if LOC
budget requires splitting):

- **`build_h11v`** (~25 LOC): given `v`, `P`, `Q`, `A_parent`, the M_n codomain
  `A_cod`, prec, fills a `aic_dhom_v h11v` (dim_B=n^2, codomain A_cod). Calls
  `aic_corner_ha(h11v.vE[l*n+m], A_parent, v->vE[l*n+m], P, P, Q, prec)` for
  each (l,m). [NOTE: the `Ha` matrix is init'd n x n BEFORE the call; `aic_corner_ha`
  asserts the output is d_PQ x d_RQ = n x n — MUST be pre-allocated correctly.]

- **`extract_U1`** (~20 LOC): given `mu11` (`aic_dhom_v`, the lem_approx output),
  the corner-basis operators `C_PQ[]` (dim S_{P,Q} = n operators), fills a
  contiguous array of n `acb_mat_t` (each n x n), `U1_ops[l]` = the physical
  n x n operator representing U_1(e_l). Uses `aic_latd_svd` for the rank-1 SVD
  of mu11.vE[0], then the column-extraction formula from §4.2. Returns these n
  physical operators.

- **`build_gamma_arrays`** (~30 LOC): given `v`, `U1_ops[0..n-1]`, `Qtilde`,
  allocates and fills the four flat arrays gamma11[n*n], gamma12[n], gamma21[n],
  gamma22[1] of `acb_mat_t` as described in §5.2.

No other new functions are needed. The main assembly logic calls (in order):
1. `aic_corner_dim_S` (loop for P_j = v->vE[j*n+j], assert each == 1)
2. `aic_corner_dim_S` on (P,Q), assert == n
3. `aic_cstar_matrix_algebra` for A_cod (M_n)
4. `build_h11v`
5. `aic_dhom_approx` on h11v (lem_approx)
6. `aic_corner_extract` for C_PQ basis
7. `extract_U1`
8. `aic_corner_Ptilde` for Qtilde
9. `build_gamma_arrays`
10. `aic_cstar_lem_merging` with `two_block=0`, n1=n, n2=1

### 6.2 Reused (CLOSED, green)

All reused and tested. No changes needed:
- `aic_corner_ha`, `aic_corner_dim_S`, `aic_corner_extract`, `aic_corner_Ptilde`
  (include/aic_corner.h — CLOSED)
- `aic_dhom_approx`, `aic_dhom_v_init`, `aic_dhom_v_clear`, `aic_dhom_B_init`,
  `aic_dhom_B_index` (include/aic_dhom.h — CLOSED)
- `aic_cstar_matrix_algebra`, `aic_cstar_lem_merging`
  (include/aic_cstar.h, I1+I3 — CLOSED)
- `aic_latd_svd` (include/aic_latd.h — CLOSED)

### 6.3 Missing FLINT/LAPACK primitives

None. Design §8(c) claimed "only 3 small new pieces" — confirmed:
1. `build_h11v` (the Ha-map loop, ~25 LOC)
2. `extract_U1` (the SVD + column-extraction, ~20 LOC)
3. `build_gamma_arrays` (the array assembly, ~30 LOC)

No additional FLINT or LAPACK calls beyond those already green.

**Note on A_cod lifetime.** The `A_cod = aic_cstar_matrix_algebra(n)` ecstar is
used as the codomain of h11v and is passed to `aic_dhom_v_init`. The `aic_dhom_v`
BORROWS its A pointer. A_cod must be kept alive until `aic_dhom_v_clear(&h11v)` and
`aic_dhom_v_clear(&mu11)` have been called. It is local to `aic_cstar_lem_extension`
and is cleared before the function returns.

---

## 7. Test design (`tests/test_cstar_extension.c`)

All tests must satisfy Rule 5 (every test asserts an invariant) and Rule 6
(cross-checks > unit tests). Following the house style of `test_cstar_merging.c`.

### 7.1 T1: eta=0 oracle, base case n=1

**Setup:** A = genuine M_2 (`aic_cstar_matrix_algebra(2)`), n=1 (extending M_1 to
M_2). P = diag(1,0), Q = diag(0,1) (exact projections with ||P+Q-I||=0). v: M_1
-> S_P is the inclusion v(1) = Ptilde = P (on the genuine C* algebra S_P = span(P)
= span(E_{00}), Ptilde = Co_P(P) = P exactly). S_{P,Q} ≅ span(E_{01}) (1-dim).

**Expected outputs:**
- dim S_{P,Q} = 1 = n (Step 1 assertion passes)
- dim S_{P_0,Q} = dim S_{P,Q} = 1 (the single j=0 check)
- mu_{11} has defect ~machine-zero (h_{11}v is an exact hom since A is genuine C*)
- U_1 is a 1x1 unitary (a phase), U_1_op[0] = phase * C_0 where C_0 spans S_{P,Q}
- v_+ (= v_+: M_2 -> M_2) has mult_def ~ machine-zero and sigma_min ~ 1.0
- v_+(E_{lm}) == E_{lm} for l,m in {0,1} (v_+ is the IDENTITY inclusion into M_2,
  up to the global phase of U_1 which cancels in gamma_{12}/gamma_{21})

**Key assertion:**
```c
AIC_CHECK(mult_def < 1e-12, "T1: eta=0 mult_def not ~0");
AIC_CHECK(fabs(sigma_min - 1.0) < 1e-9, "T1: eta=0 sigma_min != 1");
/* v_+(E_lm) == E_lm: */
for (l=0; l<2; l++) for (m=0; m<2; m++) {
    AIC_CHECK(opnorm(v_out.vE[l*2+m] - E_{lm}) < 1e-11, ...);
}
```

**Note on phase.** The U_1 phase is determined by the left singular vector of
mu_{11}(E_{00}). On the exact A=M_2 oracle the phase is 1 (the SVD of E_{00}
has u_0 = e_0 with zero phase). If not, the CANCELLATION in
gamma_{12}(e_0) = U_1_op[0] and gamma_{21}(e_0^dag) = U_1_op[0]^dag means v_+
is still an exact isomorphism regardless of phase. Test the MERGED v_+, not U_1
directly, to avoid phase sensitivity.

### 7.2 T2: eta=0 oracle, inductive case n=2 (M_2 -> M_3)

**Setup:** A = M_3, P = diag(1,1,0), Q = diag(0,0,1). v: M_2 -> S_P is the
natural inclusion v(E_{lm}) = E_{lm} (the upper-left 2x2 block). The S_P wrapper
has dim_S = 4 (= 2^2). S_{P,Q} ≅ span(E_{02}, E_{12}), dim = 2 = n.

**Expected outputs:**
- dim S_{P,Q} = 2 = n (Step 1)
- dim S_{P_j, Q} = 1 for j=0,1 (P_0 = E_{00}, P_1 = E_{11}, Q = E_{22})
- v_+: M_3 -> M_3 has mult_def ~ machine-zero, sigma_min ~ 1
- v_+(E_{lm}) == E_{lm} for l,m in {0,1,2}

**Key assertion:** same pattern as T1 but for 3x3 matrices.

**dim_B of v_out:** B_out should be M_3 (single block, dim_B = 9, num_blocks=1, d[0]=3).

### 7.3 T3: oblique eta>0 fixture with c=mult_def/eta magnitude tooth (FINDINGS §C8)

**Setup:** A = `make_mixconj(5,3,t,prec)` for t in {0.06, 0.02}, eta computed as
||S_Phi^2 - S_Phi||_op (the proxy used in I2/I3 tests, see `test_cstar_merging.c`
line 54). Find 1d projections P, Q in A with P+Q close to I and S_{P,Q} != 0
(e.g. P = Ptilde for a nontrivial projection found by `aic_projection_nontrivial`
in A, Q = another from the Stage-1 loop).

**Constructing a valid (P, Q, v) input for lem_extension:** use a genuine M_2
sub-structure in `make_mixconj(5,3,t)`. On this 5-dim oblique algebra, find two
1d delta-projections P, Q with ||P+Q-I_A|| < delta and S_{P,Q} != 0 (confirmed by
aic_corner_dim_S = 1). Then v: M_1 -> S_P is v(1) = Ptilde_P. This gives n=1,
and lem_extension produces v_+: M_2 -> A.

**Expected outputs:**
- mult_def = O(eta): assert c = mult_def/eta < C_max for some C_max (~10; oblique
  fixture, the constant grows relative to the exact case but is O(1)).
- sigma_min > 0.5.
- The merge_cond_max from aic_cstar_lem_merging is O(eta): mc/eta < C_max.

**The §C8 tooth (STAR MUTATION):** After building v_+ with the correct star in
A_parent, record mult_def_star. Then temporarily mutate `A_parent.star_phi` to
point to a thunk that does PLAIN matrix multiplication (bypassing Phi_tilde). Re-run
the defect sweep on the assembled v_out. Assert:
- c_star = mult_def_star / eta < 0.2 (correct star)
- c_plain = mult_def_plain / eta > 0.2 (plain product; the §C8 gap of ~25x)
This is the DEFERRED §C8 tooth that the I3 review explicitly deferred to I4.

**Implementation note.** The mutation target is the STAR in the codomain A_parent,
which is used by `aic_dhom_defect_sweep(v_out, prec)` to compute
`v(E_i * E_j) - v(E_i) * v(E_j)`. Mutating A_parent->star_phi to the plain-product
thunk (the same thunk `aic_cstar_matrix_algebra` uses) exercises the C2 loading
path. After the mutation, restore the original star_phi.

### 7.4 Mutation teeth (must fire RED)

All four mutations must be CONFIRMED RED by hand (Rule 7 "port-and-verify"
mutation-prove discipline) before the test file is committed.

**(a) Ha-map P<->R swap in h_{11}.** Change `aic_corner_ha(Ha, A, Z, P, P, Q)`
to `aic_corner_ha(Ha, A, Z, Q, P, Q)`.
PREDICTED: `Ha` is init'd n x n but Ha^Q_{Q,P}(Z) has shape d_QQ x d_PQ = 1 x n.
The `aic_corner_ha` assert on the Ha dimensions fires (Rule 4), aborting with a
clear message. RED via assert abort.
If the code checks the asserts at compile time only (and runs with NDEBUG), the
1 x n matrix is stored in the n x n Ha, producing garbage in the upper n-1 rows.
The h11v defect sweep then reads >> 0.1 on the T1 oracle (E_{00}*E_{00} = E_{00}
but h11v(E_{00}) * h11v(E_{00}) != h11v(E_{00}) in the garbled map). RED.

**(b) star -> plain in the defect sweep.** As described in T3 §C8. The c=mult_def/eta
tooth: plain gives c > 0.2 while correct star gives c < 0.2. RED.

**(c) gamma_12 -> 0.** Replace all gamma12[l] with the zero n x n matrix.
The off-diagonal units E_{l,n} have v_+(E_{l,n}) = 0. But E_{l,n} * E_{n,l}^dag =
E_{l,l} (a nonzero rank-1 projector), while v_+(E_{l,l}) = v(E_{ll}) != 0.
So mult_def = ||v_+(E_{l,n} * E_{n,l}) - v_+(E_{l,n}) * v_+(E_{n,l})|| =
||v(E_{ll})|| ~ 1. Also sigma_min drops to 0 (the entire off-diagonal direction
collapses). PREDICTED RED on both mult_def > 0.5 and sigma_min < 0.5.

**(d) U_1 not normalized (double the first singular vector norm).** Divide
u0 by 2 before the column extraction (so U1_col_l = mu11.vE[l*n+0] * (u0/2) /
||(u0/2)||^2 = mu11.vE[l*n+0] * 2*u0 / ||u0||^2 = 2 * U1_col_l_correct).
This gives gamma12 images scaled by 2, while gamma21 images are also scaled by 2
(via acb_mat_conjugate_transpose of gamma12). Then:
- merging3 violation: ||gamma12[l]|| ~ 2 instead of ~1 (out of the [1-delta, 1+delta]
  band).
- mult_def stays ~machine-zero on the T1 oracle (the scaling cancels in products;
  gamma_12(e_l) * gamma_21(e_l^dag) = 2*X * 2*X^dag = 4 X X^dag, but the S_Q unit
  gamma_22(1) = Qtilde has norm ~1, so merging1 is violated for the cross block).
Actually: the merged v_+ has v_+(E_{l,n}) = 2 * gamma12[l] and
v_+(E_{l,n}^dag * E_{l,n}) = v_+(E_{nn}) = Qtilde ~ 1, while
v_+(E_{l,n}) * v_+(E_{l,n}^dag) ~ 4 C_l C_l^dag ~ 4 Ptilde_l.
So mult_def ~ ||Qtilde - 4 Ptilde_l|| ~ O(1) for n >= 2 (the 4 vs 1 mismatch).
PREDICTED RED on mult_def > 0.5 for n >= 2 (T2 inductive case).

---

## 8. FINDINGS-worthy items

### F-I4-1: The merging-condition verification is IMPLICIT in lem_merging's output guards (design correctness nuance, NOT a bug)

The design §4.4 Step 6 describes separately verifying the merging conditions
(merging0-3) before calling lem_merging. In the `.tex` proof (`.tex:1411`), these
conditions FOLLOW AUTOMATICALLY from ||h_{jk} gamma_{jk} - mu_{jk}|| <= O(delta+eps)
and the h_{jk}^{-1} bounds. The constructive realization is: the output guards
`mult_def` and `sigma_min` of `aic_cstar_lem_merging` ALREADY certify the merged
v_+ satisfies the O(delta+eps)-inclusion property. A separate per-block
condition sweep (merging0-3 individually) is REDUNDANT with these output guards
(the `merge_cond_max` output additionally reports the per-block conditions for
diagnostic purposes but is not the primary certification). The implementation
should gate on `mult_def` and `sigma_min`, not try to re-derive the merging
conditions from the gamma_{jk} inputs. **This is a simplification** that makes
the code ~20 LOC shorter than a literal reading of the design's Step 6.

### F-I4-2: lem_PQ_Hilb inner product vs operator norm: the `{C_l}` basis coordinate representation gives FROBENIUS, not EUCLIDEAN, norms

The paper (`.tex:1146`) equips S_{P,Q} with the Hilbert-space norm from lem_PQ_Hilb:
||X||_Euc = sqrt(<X,X>) where <Y,X> Qtilde = Y^dag . X. The header notes
"||X||_Euc = ||X|| up to 1 +/- O(delta+eps)". The `aic_corner_ha` output is in
the `{C_l}` corner-basis coordinate representation, where the {C_l} are
Frobenius-orthonormal. Since <C_l, C_m>_F = delta_{lm} and
<C_l, C_m>_Euc ~ delta_{lm} + O(delta+eps) (lem_PQ_Hilb, `.tex:1130`), the
coordinate inner product on C^n (standard dot product of n-vectors) approximates
the Euclidean inner product up to O(delta+eps). For the U_1 extraction this means
the column orthogonality ||u_l||^2 ~ ||U1_op[l]||_op^2 up to O(delta+eps), so the
SVD-based extraction is sound. Fail-loud: assert ||u_0||^2 > 0.5 before dividing.

### F-I4-3: `aic_corner_ha` builds its Co/basis INTERNALLY — no preallocated basis needed

The header states `aic_corner_ha` "Builds all Co's and bases internally from P,R,Q."
This means the `{C_l^{PQ}}` and `{C_m^{RQ}}` bases are rebuilt fresh on every call.
For the `build_h11v` loop (n^2 calls for n^2 matrix units), this rebuilds the same
Co_{P,Q} and basis n^2 times. At n=6 (d_PQ=6), each `aic_corner_ha` call solves a
6x6 Gram system; the loop cost is n^4 solves total (~1296 at n=6). This is the §G6
concern from the design doc. For n <= 6 it is fast; for n > 8 the loop may become
a bottleneck (profile in I4 and flag if > 100ms per call). No code change is needed
for I4 but the profiling gate should be in the test report.

**Optimization note (not for I4, deferred):** extract Co_{P,Q} and {C_l^{PQ}}
ONCE before the loop and pass them to a variant `aic_corner_ha_with_basis` that
accepts preallocated bases. This is a future cycle.

### F-I4-4: U_1 is n x n UNITARY in the coordinate basis (not n x n ambient isometry)

U_1 maps C^n -> S_{P,Q}. In the {C_l} coordinate basis of S_{P,Q}, the image of
e_l has coordinates U1[:, l] (a column in C^n). Since the {C_l} are
Frobenius-orthonormal and the Hilbert norm ~= Frobenius norm (§F-I4-2), U_1 in
the coordinate basis is UNITARY (U1^dag U1 = I_n in the coordinate sense). The
physical operators U1_op[l] = sum_m U1[m,l] C_m are Frobenius-orthonormal as
physical operators. This must be asserted (||<U1_op[l], U1_op[m]>_F - delta_{lm}||
< 0.01 for the T1/T2 eta=0 oracle; at eta>0 up to O(delta+eps)).

### F-I4-5: The `aic_cstar_lem_extension` prototype (preliminary)

```c
/* aic_cstar_lem_extension (.tex:1378-1411): extend the delta-isomorphism
 *   v : M_n -> S_P  (passed as aic_dhom_v)
 * to an O(delta+eps)-isomorphism
 *   v_+ : M_{n+1} -> A
 * given dim S_Q = 1 and S_{P,Q} != 0.
 *
 * INPUTS:
 *   v       : BORROWED delta-isomorphism M_n -> S_P (as aic_dhom_v; v->A must
 *             point to A_parent; v->B must be a single M_n block).
 *   P, Q    : n x n Hermitian delta-projections in A_parent with ||P+Q-I||<=delta.
 *   A_parent: the parent eps-C* algebra (BORROWED).
 *
 * OUTPUTS:
 *   B_out, v_out : the assembled M_{n+1} and v_+ : M_{n+1} -> A_parent.
 *                  v_out->B points at B_out (keep alive); v_out->A = A_parent.
 *                  Free with aic_dhom_B_clear + aic_dhom_v_clear.
 *   mult_def     : OUTPUT arb ball, aic_dhom_defect_sweep(v_out) w.r.t. A_parent.
 *                  NULL to skip.
 *   sigma_min    : OUTPUT arb ball, aic_dhom_v_sigma_min(v_out). NULL to skip.
 *
 * ASSERTS (fail loud, Rule 4):
 *   v->B->num_blocks == 1 and v->B->d[0] == n (single M_n block)
 *   v->A == A_parent
 *   dim S_{P_j,Q} == 1 for each j (aic_corner_dim_S, the S_{P,Q}!=0 hypothesis)
 *   dim S_{P,Q} == n (lem_add_dim conclusion)
 *   lem_approx converges (aic_dhom_approx does not hit max_steps)
 *   ||u_0||^2 > 0.5 (the rank-1 SVD of mu_{11}(E_{00}) is nondegenerate)
 */
void aic_cstar_lem_extension(aic_dhom_B *B_out, aic_dhom_v *v_out,
                             arb_t mult_def, arb_t sigma_min,
                             const aic_dhom_v *v,
                             const acb_mat_t P, const acb_mat_t Q,
                             const aic_ecstar *A_parent, slong prec);
```

This should be declared in `include/aic_cstar.h` in an `==== Increment 4 ====`
block, consistent with the existing I2/I3 blocks.

---

## 9. Summary: resolved items and deferred findings

### Resolved

1. **Ha-map argument order** (the highest-risk item): CONFIRMED. The design's
   `aic_corner_ha(Ha, A, Z, P, R, Q)` with `(P, R, Q) = (P_j, P_k, Q_fixed)` is
   correct. For the four h_{jk}: h_{11} uses (P,P,Q), h_{12} uses (P,Q,Q),
   h_{21} uses (Q,P,Q), h_{22} uses (Q,Q,Q). The mutation (swap P↔R in h_{11})
   will trip the Ha dimension assert (n×n init'd, 1×n required).

2. **M_{n+1} block-index convention**: n1=n (P block, indices 0..n-1), n2=1 (Q
   block, index n). The Q block is at global (n,n). gamma12[l] corresponds to
   E_{l,n} (global). The design's "E_{00} for Q block" refers to LOCAL block-index
   (0,0) within M_1, which maps to global (n,n). CONSISTENT.

3. **mu_{11} via lem_approx**: eps_target=0.0, tol_abs=1e-12, unit_tol=1e-10,
   max_steps=30. The codomain A_cod is `aic_cstar_matrix_algebra(n)` (genuine M_n).
   No AIC_ERRREDUCE_EPS_FLOOR applies.

4. **U_1 extraction**: `aic_latd_svd` on mu11.vE[0] gives u_0 = U_left[:,0].
   Then U1_col_l = mu11.vE[l*n+0] * u_0 / ||u_0||^2 (coordinate vector in C^n).
   Physical operator: U1_op[l] = sum_m U1[m,l] C_m_PQ. Both n x n.

5. **gamma array layout**: confirmed against `aic_cstar_lem_merging` header. See
   §5.2 for exact formulas.

6. **New code**: 3 static helpers + 1 exported function, ~120 LOC. No missing
   FLINT/LAPACK primitives.

7. **The §C8 deferred tooth (c=defect/eta magnitude)** is I4's responsibility,
   realized in T3.

### New findings logged

- F-I4-1: merging conditions are implicitly certified by lem_merging output guards.
- F-I4-2: lem_PQ_Hilb inner product ~= Frobenius norm up to O(delta+eps) in coords.
- F-I4-3: aic_corner_ha rebuilds Co/basis on every call — n^4 total; profile.
- F-I4-4: U_1 unitary in coordinate basis; physical ops Frobenius-orthonormal.
- F-I4-5: Preliminary `aic_cstar_lem_extension` prototype.

---

## 10. Orchestrator verification addendum (2026-05-30)

The orchestrator independently re-verified §2 (Ha-map order) and §5 (block/gamma
conventions) against `include/aic_corner.h` and `include/aic_cstar.h` and confirms
them. Four refinements/cautions for the implementer, BINDING on top of the spec:

**A. The h_{11} mutation label is wrong (substance is fine).** §2.4 / §7.4(a) call
the tooth "swap P↔R in the h_{11} call." h_{11} is `aic_corner_ha(...,P,P,Q)` — the
two swapped args are BOTH `P`, so a literal P↔R swap is a NO-OP and would NOT be RED.
The agent's *written result* `(Q,P,Q)` is a genuinely different mutation (computes
Ha^Q_{Q,P}, output 1×n into an n×n slot → trips the `aic_corner_ha` dimension assert
→ RED). Use THAT as the tooth, and label it correctly: "feed h_{11} the wrong
projection pair (Q,P,Q) instead of (P,P,Q)". ALSO add a subtler, dimension-preserving
tooth that exercises the index convention without a gross shape error: e.g. swap the
roles of h_{12} and h_{21} conventions, or transpose the assembled h11v image — these
must drive the lem_approx defect O(1) (RED) on the T2 oracle.

**B. T3 oblique fixture — ‖P+Q−unit_A‖ ≤ δ is a HARD hypothesis (do not produce a
vacuous tooth).** `lem_extension` requires P (rank n) and Q (rank 1) to be
COMPLEMENTARY: P+Q ≈ the unit of the algebra they live in (`.tex:1378`,
‖P+Q−I‖≤δ). Two 1-dim projections in the rank-5 `make_mixconj(5,3)` algebra
(unit = 1_5) do NOT satisfy this, so the naive §7.3 setup would feed an INVALID
input and the c=mult_def/η tooth would be meaningless (a "test that can't fail").
The implementer MUST construct a fixture where P+Q genuinely approximates the
relevant unit, e.g.:
  - build an S_P wrapper (via `aic_cstar_subalg_build`) of a 2-dimensional corner of
    an oblique channel, whose unit is Ptilde, with P, Q complementary 1-d projections
    summing to Ptilde (P+Q ≈ Ptilde); v: M_1 → S_P, then v_+ : M_2 → S_P-context; OR
  - choose an oblique channel whose A=ImgΦ̃ has a genuine 2-dim block and set P,Q to
    its two 1-d sub-projections.
The test MUST PRINT ‖P+Q−unit‖ and ASSERT it is ≤ O(η) (fail loud) BEFORE trusting
mult_def. If you cannot construct a valid oblique n=1 input cheaply, escalate to the
orchestrator rather than ship a vacuous tooth — but a 2×2 oblique block is
constructible (the assoc/idemp test channels have block structure).

**C. Verify `aic_corner_ha` on the genuine matrix algebra early.** The η=0 oracle
T1/T2 use A = `aic_cstar_matrix_algebra` (genuine M_d, star = plain product). The
corner machinery (Co_{P,Q}, Ha-map) must reduce to the exact compression on this
genuine algebra. Smoke-test `aic_corner_ha` on A=M_2 with P=diag(1,0),Q=diag(0,1)
BEFORE wiring the full pipeline, so a corner-vs-genuine incompatibility surfaces
immediately rather than as a confusing downstream defect.

**D. Profiling (F-I4-3).** The n^4 `aic_corner_ha` rebuild cost is acceptable for I4;
do NOT optimize. Just report per-call wall time at the largest n you test (flag in the
increment report if > ~100 ms/call). The `aic_corner_ha_with_basis` optimization is a
deferred future cycle, not I4.

---

*End of spec.*
