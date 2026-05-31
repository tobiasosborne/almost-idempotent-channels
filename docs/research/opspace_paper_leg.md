# Paper leg: `opspace` — mathematical content of §10 (`th_main_ext`)

**Ground truth:** `paper/src/approximate_algebras.tex` (cited by line number throughout)
**Related files:** `paper/FINDINGS.md §C12 §D3`; `paper/shards/shard-F-tensor-extensions.md`
**Date:** 2026-05-31
**Purpose:** Pin down the exact mathematical content of §10 before implementation,
confirm which cb-norm route is faithful, verify the vacuity finding of §C12.

---

## Executive summary (5 lines)

(a) **The cb-inclusion `a_n` IS the operator-norm inclusion inf**, not the
Frobenius coordinate σ_min. `.tex:1489` defines `a_n = inf{‖(1_{M_n}⊗v)(X)‖_n /
‖X‖_n : X≠0}` where `‖·‖_n` is the operator norm on M_{nN} (largest singular
value), as established by `def:opspace` `.tex:1454–1462` (Ruan axioms R1, R2 for
the concrete operator-space structure). FINDINGS §C12 is correct and binding.

(b) **Two operator-norm quantities to measure:** the forward norm `‖1_{M_n}⊗v‖_op`
(max-stretch, i.e. `σ_max` of the ampliated map, certifies `‖v‖_cb ≤ 1+O(ε)`) and
the lower isometry constant `a_n = inf{‖(1_{M_n}⊗v)(X)‖_n / ‖X‖_n}` (min-stretch,
certifies `‖v^{-1}‖_cb ≤ 1+O(ε)`). When `v_n` is bijective, `a_n = 1/‖v_n^{-1}‖_op`
exactly, and the two quantities are each other's reciprocals.

(c) **v is reused unchanged from `aic_cstar_build`.** `.tex:1542,1557` confirm: B and v
come from the §9 master loop unchanged; `opspace` adds only a cb-norm certification
pass over the existing v. No new construction.

(d) **Precondition guards:** `a_1 ≥ η > 2δ` (`.tex:1504–1505`) is the induction
precondition; the induction also requires `(1−δ')/2 > 2δ` (`δ,ε` sufficiently small,
`.tex:1505`). Both must be checked explicitly with fail-loud aborts.

(e) **The truncation-N question is OPEN.** The paper gives no explicit finite N_max;
the induction proves the all-n bound analytically, and the finite-N certificate
depends on an operator-space theorem (Smith's lemma / Wittstock) that may make it
a THEOREM rather than a conjecture. This is the key open question for a web research
leg (FINDINGS §D3, bead aic-2jd).

---

## 1. The norm in prop_inc_ext and `def:opspace`

### 1.1 What `def:opspace` says (`.tex:1453–1464`)

The verbatim definition (`.tex:1453–1464`):

> A complex vector space L is called an *operator space* if each space
> M_n ⊗ L (n = 1,2,...) is equipped with a norm ‖·‖_n satisfying the
> following axioms:
>
>   (R1)  ‖AXB‖_n ≤ ‖A‖ ‖X‖_k ‖B‖    (A∈M_{n,k}, B∈M_{k,n}, X∈M_k⊗L)
>
>   (R2)  ‖diag(X,Y)‖_{k+n} = max{‖X‖_k, ‖Y‖_n}   (X∈M_k⊗L, Y∈M_n⊗L)

The norm on L itself is the n=1 case. (R1) is the natural submultiplicativity
of operator norms under left/right multiplication by matrices; (R2) says the
block-diagonal norm is the max of the blocks.

For the concrete case A ⊆ M_N (a subspace of an ambient matrix algebra),
`.tex:1465` (immediately after the definition) states: "In particular, all
linear subspaces of C* algebras satisfy axioms (R1), (R2)." And `.tex:1475`:
"the map X ↦ I_n ⊗ X is an isometric inclusion of L into M_n ⊗ L." The
inherited norm on M_n ⊗ A ⊆ M_{nN}(C) is the **operator norm** (largest
singular value) on nN×nN matrices. There is no freedom to choose a different
norm; the operator norm is the unique norm on a C*-algebra (and the induced
operator norm is the only one consistent with the C* identity).

**Conclusion.** `‖·‖_n` in `def:opspace` and throughout §10 is the **operator
norm** of an nN×nN matrix (nN = n × (ambient dim of H)). It is NOT the
Frobenius norm, not the coordinate norm, not σ_min of a Frobenius coordinate
matrix. This confirms FINDINGS §C12 and refutes the `opspace_design.md §3.2`
route.

### 1.2 The `a_n` definition (`.tex:1487–1491`)

Verbatim from the proof of prop_inc_ext (`.tex:1487–1491`):

```
a_n = inf { ‖(1_{M_n} ⊗ v)(X)‖_n / ‖X‖_n  :  X ∈ M_n ⊗ A', X ≠ 0 }

     (a_1 ≥ a_2 ≥ a_3 ≥ ... ≥ 0)
```

Both the numerator and denominator use `‖·‖_n`, the operator norm on M_{nN}.
This is the operator-norm lower isometry constant of the ampliated map
`1_{M_n}⊗v : M_n⊗A' → M_n⊗A''`, i.e. the smallest singular value of that map
when viewed as a linear operator between operator spaces. It is **not** the
Frobenius-coordinate singular value used in `aic_dhom_v_sigma_min` (§C6).

The vacuity of the Frobenius route (FINDINGS §C12): if one computes `a_n` as
σ_min of the ampliated coordinate matrix `I_{n²}⊗M_1`, one gets
`σ_min(I_{n²}⊗M_1) = σ_min(M_1)` for all n by elementary linear algebra
(the singular values of `I_{n²}⊗M_1` are those of `M_1` repeated n² times).
This is true for ANY linear v regardless of its quality — it tests nothing about
the operator-norm ampliation. It is the exact "test that cannot fail" pattern
(CLAUDE.md Rule 5).

---

## 2. The doubling step `a_{2n} ≥ a_n/2` (`.tex:1493–1503`)

### 2.1 The two ingredients

The proof uses two operator-space facts:

**Fact 1 (upper bound on `‖X‖_{2n}`, `.tex:1494–1498`):** For any
X = [[X_{11}, X_{12}], [X_{21}, X_{22}]] with X_{pq} ∈ M_n⊗A',

```
‖X‖_{2n} ≤ 2 ‖X‖_{n,max},    where ‖X‖_{n,max} = max_{p,q} ‖X_{pq}‖_n.
```

Derivation (`.tex:1498`): decompose

```
X = [[X_{11}, 0], [0, X_{22}]]  +  [[0, I_n], [I_n, 0]] [[X_{21}, 0], [0, X_{12}]]
```

The first term has norm `max{‖X_{11}‖_n, ‖X_{22}‖_n} ≤ ‖X‖_{n,max}` by
axiom (R2). The second term: the swap matrix `[[0,I_n],[I_n,0]]` has operator
norm 1 (it is a unitary permutation), and `[[X_{21},0],[0,X_{12}]]` has norm
`‖X‖_{n,max}` by (R2). By (R1) applied to the product, the second term has
norm ≤ 1 · `‖X‖_{n,max}` · 1 = `‖X‖_{n,max}`. Triangle inequality gives
`‖X‖_{2n} ≤ 2‖X‖_{n,max}`.

**Fact 2 (lower bound on `‖(1_{M_{2n}}⊗v)(X)‖_{2n}`, `.tex:1499–1501`):**

```
‖(1_{M_{2n}}⊗v)(X)‖_{2n} ≥ ‖(1_{M_{2n}}⊗v)(X)‖_{n,max} ≥ a_n ‖X‖_{n,max}
```

The first inequality: `‖Y‖_{2n} ≥ ‖Y_{pq}‖_n` for any block Y (this is the
reversal of (R2): the block-max norm is a lower bound). Formally, by (R1)
applied to `Y_{pq} = e_p^† Y e_q` (where e_p are the standard basis column
vectors), `‖Y_{pq}‖_n ≤ ‖e_p^†‖ ‖Y‖_{2n} ‖e_q‖ = ‖Y‖_{2n}`, so
`‖Y‖_{n,max} ≤ ‖Y‖_{2n}`. The second: apply the definition of `a_n` to each
block `(1_{M_n}⊗v)(X_{pq})`.

### 2.2 Combining to get `a_{2n} ≥ a_n/2`

From Facts 1 and 2:

```
‖(1_{M_{2n}}⊗v)(X)‖_{2n} ≥ a_n ‖X‖_{n,max} ≥ (a_n/2) ‖X‖_{2n}
```

Since this holds for all X≠0, `a_{2n} ≥ a_n/2`. (`.tex:1503`)

### 2.3 Is 1/2 tight?

The constant 1/2 is tight in the following sense: it comes from the exact
decomposition at `.tex:1498`, which uses exactly 2 terms. The factor 2 in
`‖X‖_{2n} ≤ 2‖X‖_{n,max}` and the factor 2 in `a_{2n} ≥ a_n/2` are the same
constant from the same decomposition. For the X that achieves `‖X‖_{2n} =
2‖X‖_{n,max}` (a rank-1 off-diagonal block X with equal-norm diagonal and
off-diagonal parts), the inequality is sharp. So 1/2 is tight as a bound on the
doubling step; it is not merely an artifact of a crude estimate.

---

## 3. prop_delta_hominc: the exact `δ'` and the induction

### 3.1 The prop_delta_hominc bound (`.tex:1194–1215`)

prop_delta_hominc (`.tex:1194–1196`) states: if `‖v(X)‖ ≥ η‖X‖` for some
η > 2δ, then `‖v(X)‖ ≥ a‖X‖` for `a = 1 − O(δ+ε)`. The proof derives (`.tex:1210–1215`)
the inequality

```
a ≥ sqrt( ((1−ε)a − δ) / (1+ε) )
```

and shows this implies `a ≥ 1 − O(δ+ε)`. The paper states the bound as
`1 − O(δ+ε)` with an unspecified analytic constant (the section preamble at
`.tex:1192` warns: "Only in Corollary cor_improvement, which is directly used
in the proof of the main theorem, the constants are mentioned explicitly"). So
`δ' = O(δ+ε)` but the exact rational function of δ, ε is implicit.

**Extracting the analytic form:** from the fixed-point inequality
`a ≥ sqrt(((1−ε)a−δ)/(1+ε))`, squaring and rearranging:
`a²(1+ε) ≥ (1−ε)a − δ` → `a²(1+ε) − (1−ε)a + δ ≥ 0`. At `a = 1 − δ'`
with `δ' = O(δ+ε)` small, the dominant terms give `δ' ≈ (2ε+δ)/(2−ε) ≈
δ + 2ε` to leading order. The analytic constant is not made explicit by the
paper at this point; FINDINGS §D2 tracks the open question of the exact `c_0`
(bead aic-1bc).

The key qualitative fact is that `δ'` depends ONLY on `δ` and `ε`, not on n,
not on the dimension of A, and not on any properties of v beyond that it is a
δ-homomorphism with a ≥ η > 2δ.

### 3.2 How the induction works (`.tex:1504–1505`)

Verbatim (`.tex:1504–1505`):

> By Proposition prop_delta_hominc, if a_n > 2δ, then a_n ≥ 1−δ' for some
> δ' = O(δ+ε) that does not depend on n. We have a_1 ≥ η > 2δ, and hence,
> a_1 ≥ 1−δ'. It follows by induction that a_n ≥ 1−δ' for all n, provided
> (1−δ')/2 > 2δ. The last inequality holds under the standard assumption that
> δ and ε are sufficiently small.

The induction structure:
- **Base:** `a_1 ≥ η > 2δ`. Apply prop_delta_hominc at n=1 (a δ-homomorphism
  at the scalar level): gives `a_1 ≥ 1−δ'`.
- **Step:** Suppose `a_n ≥ 1−δ'`. Then `a_{2n} ≥ a_n/2 ≥ (1−δ')/2 > 2δ`
  (by the precondition `(1−δ')/2 > 2δ`). Now `1_{M_{2n}}⊗v` is itself a
  δ-homomorphism (v is an extended δ-homomorphism by hypothesis, `.tex:1478`),
  and `a_{2n} > 2δ`. Apply prop_delta_hominc at ampliation level 2n: gives
  `a_{2n} ≥ 1−δ'`.
- **Coverage:** powers of 2 get n = 1,2,4,8,.... For non-powers-of-2, use
  that a_n is non-increasing (`a_1 ≥ a_2 ≥ ...`, `.tex:1491`) combined with
  `a_{2n} ≥ 1−δ'` for all powers of 2: for any n, find the smallest power of 2
  k ≥ n, so `a_n ≥ a_k ≥ 1−δ'`.

Wait — the non-increasing claim needs care. The sequence is stated as
`a_1 ≥ a_2 ≥ a_3 ≥ ...` at `.tex:1491`. Given `a_{2n} ≥ 1−δ'` for all
powers of 2, and non-increasing, any n satisfies `a_n ≥ a_{2^k} ≥ 1−δ'`
where 2^k ≥ n (the smallest such power). This covers all n. The induction is
complete.

### 3.3 The precondition `(1−δ')/2 > 2δ`

This is the key smallness condition (`.tex:1505`). Rearranging:
`1 − δ' > 4δ`, i.e. `δ' < 1 − 4δ`, i.e. `O(δ+ε) < 1 − 4δ`, which holds
for sufficiently small δ, ε. The paper gives no explicit threshold; it is
governed by the analytic form of δ' from prop_delta_hominc (see §3.1).
The implementer must CHECK this condition with a fail-loud abort: if the
measured `a_1` and `iso_def` (= δ from cstar_build) satisfy `(a_1−δ')/2 ≤ 2δ`,
the induction breaks and the module must abort (Rule 4, stop condition).

In practice: `aic_cstar_build` already asserts `iso_def = c_0·ε` with `c_0 ≈
2–3` (FINDINGS §D2), and `a_1 = σ_min ≈ 1 − O(ε)` (e.g. `σ_min = 0.974–0.997`
in the test corpus). With ε small these conditions hold with large margins.

---

## 4. The cb-isomorphism goal and the two quantities to measure

### 4.1 What th_main_ext asserts (`.tex:1538–1540`)

Verbatim:

> For any finite-dimensional extended ε-C* algebra A, there exist a C* algebra
> B and an extended O(ε)-isomorphism v: B → A. (The implicit constant in O(ε)
> does not depend on A or its dimensionality.)

"Extended O(ε)-isomorphism" means (per `def:extended_eps_cstar`, `.tex:1477–1481`,
and the shard-F cb-norm context note):

1. `v: B → A` is a δ-isomorphism at n=1 (bijectivity + both n=1 bounds).
2. For every n: `1_{M_n}⊗v` is a δ-isomorphism with the SAME δ = O(ε),
   constant independent of n.

This is exactly `‖v‖_cb ≤ 1+O(ε)` AND `‖v^{-1}‖_cb ≤ 1+O(ε)` uniformly in n.

### 4.2 The two operator-norm quantities

**Forward max-stretch:**
```
‖1_{M_n}⊗v‖_op  =  sup_{X≠0} ‖(1_{M_n}⊗v)(X)‖_n / ‖X‖_n
                 =  σ_max of the ampliated map (largest singular value)
```
This certifies `‖v‖_cb ≤ 1 + O(ε)`. From prop_delta_hominc (`.tex:1194–1196`):
`‖v‖ ≤ 1 + O(δ+ε)` at n=1, and by extension to all n for an extended δ-hom
(the same bound applies level-by-level since `1_{M_n}⊗v` is itself a
δ-homomorphism).

**Lower isometry constant (min-stretch):**
```
a_n  =  inf_{X≠0} ‖(1_{M_n}⊗v)(X)‖_n / ‖X‖_n
      =  smallest singular value of 1_{M_n}⊗v as a linear map on M_n⊗B
```
This certifies `‖v^{-1}‖_cb ≤ 1/(1−O(ε)) ≤ 1+O(ε)`. By prop_inc_ext,
`a_n ≥ 1−δ'` uniformly in n.

### 4.3 The identity `a_n = 1/‖(1_{M_n}⊗v)^{-1}‖_op`

This identity holds exactly when `v_n = 1_{M_n}⊗v` is bijective. Proof:
by definition, `a_n = inf_{X≠0} ‖v_n(X)‖/‖X‖`. When v_n is bijective,
for any Y = v_n(X) one has X = v_n^{-1}(Y), so
`‖v_n(X)‖/‖X‖ = ‖Y‖/‖v_n^{-1}(Y)‖`, and
`a_n = inf_{Y≠0} ‖Y‖/‖v_n^{-1}(Y)‖ = 1/sup_{Y≠0} ‖v_n^{-1}(Y)‖/‖Y‖
= 1/‖v_n^{-1}‖_op`.

**Conditions for bijectivity of v_n:** `v_n = 1_{M_n}⊗v` is bijective iff v
itself is bijective (since tensoring with I_{M_n} over a field is exact —
M_n⊗(·) is an exact functor on finite-dimensional vector spaces). So v_n is
bijective for all n simultaneously, contingent on v being bijective. The paper
asserts v is an extended δ-isomorphism (hence bijective), and `aic_cstar_build`
already certifies bijectivity via `σ_min > 0` in the n=1 Frobenius check (a
genuine bijectivity certificate for n=1; for n>1 the same exactness argument
applies).

---

## 5. Does the construction of v change? (`.tex:1542, 1557`)

### 5.1 `.tex:1542` verbatim

> We adapt the proof of the main theorem. All new objects, such as nontrivial
> δ-projections, are constructed for the original algebra A, but we need to show
> that their tensor extensions have similar properties. In the case of
> δ-projections, that is straightforward: if P ∈ A satisfies P† = P and
> ‖P²−P‖ ≤ δ, then so does the element I_n⊗P ∈ M_n⊗A because the map P ↦ I_n⊗P
> commutes with the involution and the multiplication as well as preserving the
> norm.

**v is NOT recomputed.** The same B and the same v from §9 (`aic_cstar_build`)
are used. The proof of th_main_ext is a CERTIFICATION argument: it shows that
the existing §9 construction satisfies the extended conditions. No new object
is created for the cb-norm control.

### 5.2 `.tex:1557` verbatim

> Corollary cor_improvement (error reduction) should be adapted to extended
> inclusions using Lemma lem_approx_ext and Proposition prop_inc_ext. The
> arguments in Section sec_proof_main require only trivial modifications, namely,
> one should use the norms ‖·‖_n in certain places.

"Only trivial modifications" and "one should use the norms ‖·‖_n": this
confirms that the §9 argument is reprised with the ampliated norms but no new
objects. The existing v from `aic_cstar_build` is the output; `opspace` is a
post-hoc certification module.

### 5.3 lem_approx_ext adds no new computation at n=1 (`.tex:1508–1536`)

lem_approx_ext says: the correction ṽ is defined at n=1 and the ampliation
`ṽ_n = 1_{M_n}⊗ṽ` is automatic. The proof (`.tex:1512–1535`) reduces the n>1
correction to the n=1 scalar identity by working matrix-element-wise
(`[X]_{kp} ∈ B`). Specifically, the correction formula at level n:

```
w_n'(X) = Σ_j v_n(I_n⊗A_j) g_n(I_n⊗B_j, X)
```

(`.tex:1519`) uses only `v(A_j)` and `g(B_j, [X]_{kp})` — both level-1
quantities. No additional data or computation at level n is needed; `ṽ_n =
1_{M_n}⊗ṽ` by construction. The `aic_cstar_build` pipeline (which calls
`aic_dhom_approx` = lem_approx) already produces ṽ; the extended version is
the same ṽ viewed as generating `ṽ_n` for free.

---

## 6. The truncation question (D3)

### 6.1 Does §10 state any finite N?

No. The paper gives no explicit finite N_max. The prop_inc_ext proof covers
"for all n" by induction. The theorem statement (`.tex:1538`) quantifies over
all n implicitly. The paper cites the Ruan representation theorem
(`.tex:1465`, `[Rua88]`, `[Paulsen, Thm 13.4]`) as background, but does not
use it to derive a finite truncation bound.

### 6.2 The shard-F candidate: N_max = dim A

Shard-F open question 1 argues: any d-dimensional operator space embeds
isometrically into B(C^d) (Ruan 1988), so the operator-space structure is
determined by d = dim(A). For the cb-norm of v: B → A with dim B = dim A
(the bijective case), a matrix-rank argument suggests the worst-case ampliation
is at most level dim A. This gives N_max = dim A as a candidate.

### 6.3 Smith's lemma: would it make this a THEOREM?

**Smith's lemma** (classical operator-space theory) states: for a bounded linear
map φ: A → M_k (where A is an operator space and M_k is a matrix algebra),
`‖φ‖_cb = ‖φ_k‖_op = ‖(1_{M_k}⊗φ)‖_op`. In other words, for a map INTO an
N×N matrix algebra, the cb-norm is attained at ampliation level N. If this
applies to the present setting (v: B → A, with A itself embedded in M_N), it
would give:

```
‖v‖_cb = sup_n ‖1_{M_n}⊗v‖_op  =  ‖1_{M_N}⊗v‖_op   (where N = dim H, ambient)
```

making the truncation problem a THEOREM (N_max = dim H, the ambient Hilbert
space dimension, not merely the algebra dimension). Whether Smith's lemma applies
directly here depends on: (a) whether the target A is being viewed as M_N (the
ambient) or as A ⊆ M_N (a subspace), and (b) the precise formulation of Smith's
lemma for maps between operator spaces vs. maps into a matrix algebra.

**Flag for the web leg.** This is the key question: does Smith's lemma (or a
generalization — see Paulsen's book, chapter 7/8, or Effros–Ruan's operator
space theory) give an explicit N_max as a function of dim B and dim A? If yes,
the D3 "conjecture" becomes a theorem and the implementation has a certified
finite N. If no, a more careful analysis of the Ruan representation (`.tex:1465`)
is needed. **Do NOT resolve this here; flag it as the primary open question
for a web research leg (FINDINGS §D3, bead aic-2jd).**

### 6.4 The induction is still rigorous without the truncation

The key point: the prop_inc_ext induction already proves the all-n bound
unconditionally. The truncation question affects only the finite-N CERTIFICATION
(how many levels must be checked numerically to guarantee the bound for all n).
If Smith's lemma applies at N_max = dim A, the certification is exact. If not,
the induction proof still holds and the algorithm can check a conservative
N_max = dim A as a guard.

---

## Summary of findings

| Question | Answer | .tex cite |
|----------|--------|-----------|
| Norm in a_n | Operator norm (largest singular value of nN×nN matrix) | 1453–1464, 1489 |
| a_n = Frobenius σ_min? | NO. Vacuous for any linear v (§C12). | FINDINGS §C12 |
| Two quantities to measure | ‖1_{M_n}⊗v‖_op (forward) and a_n (lower) | 1489, 1200–1202 |
| a_n = 1/‖v_n^{-1}‖_op? | Yes, exactly, when v_n bijective (all n since v bijective) | — |
| a_{2n} ≥ a_n/2 tight? | Yes; 1/2 is sharp from the block decomposition | 1493–1503 |
| δ' explicit? | No; O(δ+ε) with unspecified constant (see §D2) | 1215, 1505 |
| Precondition for induction | a_1 > 2δ AND (1−δ')/2 > 2δ (δ,ε small) | 1505 |
| v changes for th_main_ext? | No; same B and v from aic_cstar_build | 1542, 1557 |
| lem_approx_ext new at n=1? | No new computation; ṽ_n = 1_{M_n}⊗ṽ automatic | 1508–1536 |
| Finite N_max in paper? | None stated; D3 open | 1538–1540 |
| Smith's lemma applies? | Possibly → THEOREM, not conjecture; web leg needed | §D3 |

---

## Implementation consequence (for the orchestrator)

The `opspace` module is a **certification module over the existing v from
`aic_cstar_build`**. Its primary job is to measure the OPERATOR-norm quantities:

1. **`a_n = inf ‖(1_{M_n}⊗v)(X)‖_op / ‖X‖_op`**: the genuine cb-lower bound.
   For a finite-dim bijective v, this is the smallest singular value of `v_n`
   viewed as a linear map on `M_n⊗B` with the operator norms. Computing it
   requires the operator-norm HOPM (bead aic-0at) or equivalent, NOT the
   Frobenius σ_min reuse from `aic_dhom_v_sigma_min` (which is vacuous for n>1).

2. **`‖1_{M_n}⊗v‖_op`**: the largest singular value of the same map. This is
   bounded by `1+O(δ+ε)` from prop_delta_hominc and is straightforwardly
   computed as σ_max of the nN×nN ampliated operator.

The `σ_min(I_{n²}⊗M_1) = σ_min(M_1)` route in `opspace_design.md §3.2` is
VACUOUS and must not be implemented. The ORCHESTRATOR CORRECTION at the end of
`opspace_design.md` is correct and binding.
