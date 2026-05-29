# ecstar paper leg — ε-C* axioms, defect estimators, and the η=0 oracle

Research subagent analysis for bead `aic-knm` (module `ecstar`).
All line numbers refer to `paper/src/approximate_algebras.tex`.
All quotes are verbatim from that file.

---

## 1. The exact ε-C* axioms (tex 407–439)

The definition block runs from L407 to L440, with the unit conditions at L430–439.

### ax_prodnorm (L410–411)

Label: `ax_prodnorm`

Verbatim:
```
\|XY\| \le (1+\eps)\ts\|X\|\ts\|Y\|\qquad (X,Y\in \calA)
```

- **Norm:** operator/spectral norm on M_N (inherited from Bo(H)). Not cb-norm. The
  cb-norm only appears in defining η-idempotence of Φ (L347–352); the ε-C* axioms
  themselves use the ordinary operator norm inherited from Bo(H) (confirmed by the
  proof of th_almost_idemp at L2213–2218, where `‖X★Y‖` is an operator norm and
  `‖Φ‖_cb ≤ 1` is used to bound it — two distinct norms).
- **Defect side:** the quantity `‖XY‖/(‖X‖‖Y‖) - 1` should be ≤ ε. The defect is
  the amount by which submultiplicativity is violated upward.

### ax_assoc (L412–413)

Label: `ax_assoc`

Verbatim:
```
\|(XY)Z-X(YZ)\| \le \eps\ts\|X\|\ts\|Y\|\ts\|Z\|\qquad (X,Y,Z\in \calA)
```

where `XY` denotes the Choi–Effros star product `X★Y = Φ̃(XY)` when A = ImgΦ̃.

- **Norm:** operator norm on M_N.
- **Defect side:** symmetric; both `(X★Y)★Z` and `X★(Y★Z)` may deviate from each
  other. The defect is the supremal ratio `‖(X★Y)★Z − X★(Y★Z)‖ / (‖X‖‖Y‖‖Z‖)`.

### ax_comm (L417–418)

Label: `ax_comm`

Verbatim:
```
\|XY-YX\|\le \eps'\ts\|X\|\ts\|Y\|\qquad (X,Y\in \calA)
```

- **Norm:** operator norm. Note ε' is a separate parameter from ε (the algebra may
  be ε-C* without being ε'-commutative for any small ε'). The paper uses ε' at
  L415: "Such an algebra is called ε'-commutative if...". Commutativity is an
  additional qualifier, not part of the core ε-C* axioms.
- **Status for our data model:** Optional estimator. Not part of ax_prodnorm /
  ax_assoc / ax_* / ax_C* / ax_eps_unit.

### ax_* (L422–423)

Label: `ax_*`

Verbatim:
```
\|X^\dag\|=\|X\|,\qquad (XY)^\dag=Y^\dag X^\dag\qquad (X,Y\in \calA)
```

- **Norm:** operator norm on M_N.
- **Both conditions are EXACT equalities**, not approximate. The text at L420 says
  "satisfying the equations" (not inequalities). This is confirmed by th_almost_idemp
  proof at L2211: "it is also evident that ... `(X★Y)†=Y†★X†` for all X,Y∈A" (exact).

Derivation that `(X★Y)† = Y†★X†` is exact for the Choi–Effros star when Φ̃ is
Hermicity-preserving:

```
(X★Y)† = (Φ̃(XY))† = Φ̃((XY)†)   [Φ̃(Z†)=Φ̃(Z)†, stated at tex L2181]
        = Φ̃(Y†X†) = Y†★X†.
```

This uses `Φ̃(X†)=Φ̃(X)†` (L2181, exact) and `(XY)†=Y†X†` in Bo(H) (exact, matrix
adjoint). So both ax_* conditions are exact for the Choi–Effros product, with no
ε correction — **the defect for ax_* is identically zero** when Φ̃ is
Hermicity-preserving, which holds (L2181).

Also `‖X†‖ = ‖X‖` is exact in Bo(H) (the operator norm is unitarily invariant).

### ax_C* (L427–428)

Label: `ax_C*`

Verbatim:
```
\|X^{\dag}X\|\ge (1-\eps)\ts\|X\|^{2}\qquad (X\in \calA)
```

- **Norm:** operator norm on M_N.
- **One-sided lower bound only.** The text at L430 immediately notes: "A bound from
  the other side, `‖X†X‖ ≤ (1+ε)‖X‖²`, follows from ax_prodnorm and ax_*." So the
  upper bound is free and the implementation does not need to estimate it separately;
  it follows from ax_prodnorm.
- **Defect direction:** the ε for which the lower bound holds is
  `ε = sup_{X≠0} (1 − ‖X†★X‖ / ‖X‖²)`. This is the "C*-identity slack". Note
  the star products: `X†★X = Φ̃(X†X)`.

### ax_eps_unit (L432–434)

Label: `ax_eps_unit` (default); `ax_exact_unit` (L435–437, if specified)

Verbatim (approximate unit, default):
```
\|XI-X\|\le\eps\ts\|X\|,\qquad
\|IX-X\|\le\eps\ts\|X\|,\qquad
\bigl|\|I\|-1\bigr|\le\eps
```

Verbatim (exact unit):
```
XI=X,\qquad IX=X,\qquad \|I\|=1
```

- **Norm:** operator norm.
- In our setting A = ImgΦ̃, the unit element is I = 1_H (the identity operator on H),
  since Φ̃(1)=1 (stated at L2180: `Φ̃(1)=1`) and 1 ∈ ImgΦ̃. For X ∈ A,
  `X★I = Φ̃(XI) = Φ̃(X) = X` (exact, since X∈A means Φ̃(X)=X). Similarly
  `I★X = X` (exact). This is confirmed at L2211: "it is also evident that
  `X★I = X = I★X`...for all X∈A" (exact).
- **So ax_eps_unit holds with ε=0 for A=ImgΦ̃.** The unit IS exact (in the
  ax_exact_unit sense) in this algebra. However, for a general ε-C* algebra the
  approximate unit only becomes exact after prop_unit (see §4 below).

---

## 2. Defect functionals

For each axiom, the defect functional an estimator must compute (working with the
stored basis {B_k}, k=1..dim_A, and the Choi–Effros star product via Φ̃):

### assoc_defect

```
assoc_defect(A) = sup_{X,Y,Z ∈ A \ {0}} ‖(X★Y)★Z − X★(Y★Z)‖ / (‖X‖‖Y‖‖Z‖)
```

In finite dimension this supremum is attained at basis vectors; in practice, evaluate
over the operator-basis set and all triples. The paper (L412) says this must be ≤ ε.

### prodnorm_defect (submultiplicativity slack)

```
prodnorm_defect(A) = sup_{X,Y ∈ A \ {0}} (‖X★Y‖ / (‖X‖‖Y‖)) − 1
```

The paper's ax_prodnorm (L410) says `‖XY‖ ≤ (1+ε)‖X‖‖Y‖`, so this supremum
equals the tightest ε for ax_prodnorm. If negative (star-product is strictly
submultiplicative), defect is 0.

### cstar_defect (C*-identity slack)

```
cstar_defect(A) = sup_{X ∈ A \ {0}} (1 − ‖X†★X‖ / ‖X‖²)
```

Per ax_C* (L427–428): must be ≤ ε. Note `X†★X = Φ̃(X†X)`. A negative value here
means the lower bound holds strictly; the defect is 0 in that case.

### involution_defect

For ax_* (L422–423), both conditions are exact for the Choi–Effros star when Φ̃ is
Hermicity-preserving. Nevertheless, for a GENERAL ε-C* algebra (not necessarily
arising from a UCP map), these could be approximate. The estimator is:

```
invnorm_defect(A) = sup_{X ∈ A \ {0}} | ‖X†‖ − ‖X‖ | / ‖X‖
```

```
invprod_defect(A) = sup_{X,Y ∈ A \ {0}} ‖(X★Y)† − Y†★X†‖ / (‖X‖‖Y‖)
```

Both must be exactly 0 for A=ImgΦ̃; in the η=0 oracle they must read 0.0.

### unit_defect

```
unit_right_defect(A) = sup_{X ∈ A \ {0}} ‖X★I − X‖ / ‖X‖
unit_left_defect(A)  = sup_{X ∈ A \ {0}} ‖I★X − X‖ / ‖X‖
unit_norm_defect(A)  = | ‖I‖ − 1 |
```

For A=ImgΦ̃ all three are exactly 0 (see §1 analysis above and L2211). For a
general ε-C* algebra, these are the ax_eps_unit defects.

---

## 3. ε vs δ vs η — distinct small parameters

### η (eta)

Introduced at L354:

```
A UCP map Φ: Bo(H) → Bo(H) is called η-idempotent if ‖Φ² − Φ‖_cb ≤ η.
```

The norm is the **cb-norm** (L347–352):
```
‖Λ‖_cb = sup_n sup_{X≠0} ‖(1_{M_n}⊗Λ)(X)‖ / ‖X‖,   X ∈ M_n ⊗ A'.
```

η measures the idempotence defect of the UCP map Φ. This is the PRIMARY input to
the system.

### ε (epsilon)

Introduced at L407–440 (the ε-C* axiom block). ε measures "how far A is from a
C* algebra" — it is the maximum of all the axiom defects listed in §1 above.

### δ (delta)

Introduced at L443–455 (the δ-homomorphism block). δ measures "how far a linear
map v between algebras is from being multiplicative (and unital)":

```
‖v(I) − I‖ ≤ δ                              [hom_unit, L446]
‖v(XY) − v(X)v(Y)‖ ≤ δ ‖X‖‖Y‖             [hom_mult, L448]
```

A δ-inclusion additionally satisfies `(1−δ)‖X‖ ≤ ‖v(X)‖ ≤ (1+δ)‖X‖` (L453).

### Relationship: η → ε for A=ImgΦ̃ (th_almost_idemp, L2192–2237)

Theorem statement at L2192–2193:

```
The space A with the norm, involution, and unit inherited from Bo(H) and the
multiplication (X,Y) ↦ X★Y is an extended O(η)-C* algebra.
```

"Extended" means not just A but also all M_n⊗A satisfy the ε-C* axioms for the
same ε = O(η), uniformly in n (definition at L1477–1478). This is stronger than
a plain ε-C* algebra — it is what one gets precisely because η-idempotence is
stated in the cb-norm.

The specific O(η) claims from the proof (L2208–2237):

- **ax_prodnorm** (L2215–2218): `‖X★Y‖ ≤ (1+O(η))‖X‖‖Y‖`
  — uses `‖Φ̃−Φ‖_cb ≤ O(η)` and `‖Φ‖_cb ≤ 1`
- **ax_assoc** (L2229–2230): `‖(X★Y)★Z − X★(Y★Z)‖ ≤ O(η)‖X‖‖Y‖‖Z‖`
  — uses equations (Phi_assoc1) and (Phi_assoc2), each O(η), via cb-norm bound
- **ax_C*** (L2234–2235): `‖X†★X‖ ≥ (1−O(η))‖X‖²`
  — uses the UCP inequality Φ(X†)Φ(X) ≤ Φ(X†X) [eq PhiXdX, L1692]
- **ax_*, unit:** exact (L2211), with 0 defect

The paper does not state the explicit constant in O(η). The O() notation is an
unspecified universal constant (L458: "each instance of big-O ... stands for a
concrete function, not depending on any additional data"). The constant must be
extracted from the proof steps above; it depends on the constant in
`‖Φ̃−Φ‖_cb ≤ c₀η` and on the constants in equations (Phi_assoc1/2).

---

## 4. The exact-unit fix (prop_unit, L672–687)

Proposition statement (L672–677):

```
Let A be an ε-Banach algebra with unit I. Then there exist a new unit J∈A and
a new multiplication denoted by a dot that make A into an O(ε)-Banach algebra
with exact unit while being O(ε)-close to the original unit and multiplication:

‖J−I‖ ≤ O(ε),  ‖X·Y − XY‖ ≤ O(ε)‖X‖‖Y‖  (X,Y∈A).

This transformation respects the involution if one is present, i.e. J†=J and
(X·Y)† = Y†·X†.
```

The construction (L682–686): find J with J²=J (existence follows from the implicit
function theorem, lem_invfun, applied to f: X↦X²−X near I, using d_X2X at L668–669);
then define `X·Y = Ra_J^{-1}(X) La_J^{-1}(Y)`.

**Implication for the defect estimator:**

1. Before applying prop_unit, the unit_defect may be up to O(ε) (the ax_eps_unit
   defects). The estimator should report the actual defect values.
2. After prop_unit, the unit_defect for ax_exact_unit is zero by construction.
3. **For our algebra A=ImgΦ̃ the unit is ALREADY exact** (L2211, as shown above),
   so prop_unit is not needed — but the bead `aic-rt1` (unitfix) should confirm this.
4. The exact unit only guarantees ax_exact_unit; the multiplication changes by
   O(ε) which means the OTHER defects (assoc_defect, cstar_defect, prodnorm_defect)
   change by at most O(ε) — still O(ε) overall.

---

## 5. The η=0 oracle: exact idempotent Φ and ALL defects = 0

For an EXACTLY idempotent UCP map Φ (η=0), the subspace A = ImgΦ = Ker(1−Φ) with
the Choi–Effros star product `X★Y = Φ(XY)` is a genuine C* algebra. This is the
Choi–Effros theorem, cited at L344:

```
Choi and Effros showed that for any idempotent UCP map Φ: Bo(H) → Bo(H) and any
underlying Hilbert space H (not only a finite-dimensional one), the subspace
A = ImgΦ = Ker(1−Φ) ⊆ Bo(H) equipped with this product, together with the norm
and the involution X↦X† inherited from Bo(H), satisfies all axioms of a C*
algebra [Theorem 3.1, ChEf77].
```

This means for η=0, all axiom defects must equal EXACTLY zero (not just small).
Let us derive each:

### ax_* (exact): zero by derivation above.

### ax_eps_unit (exact, ε=0): Φ(1)=1 for unital Φ, so 1∈A. For X∈A: `X★I = Φ(XI) = Φ(X) = X` exactly (since X∈ImgΦ = Ker(1−Φ) means Φ(X)=X). Unit defect = 0.

### ax_assoc (exact, ε=0): For X,Y,Z∈A (so Φ(X)=X, Φ(Y)=Y, Φ(Z)=Z):

```
(X★Y)★Z = Φ(Φ(XY) · Z)
```

Since X★Y = Φ(XY) ∈ A (A is closed under ★, so Φ(XY)∈A), we have
Φ(Φ(XY)) = Φ(XY) (idempotent). Therefore:

```
(X★Y)★Z = Φ(Φ(XY) · Z)  = Φ(XY · Z)   [using lem_idemp L1916–1922]
```

The key fact used here is lem_idemp (L1916–1922): for idempotent Φ and X,Y∈A,
`Φ(XY)|_M = (XY)|_M` — i.e. the Choi–Effros product restricted to the carrier M
equals the ordinary matrix product. More precisely, lem_idemp states (L1918–1919):

```
Φ(XY)|_M = (XY)|_M = X|_M Y|_M    (X,Y∈A)
```

So for X∈A, Φ(XY) can be replaced by the ordinary product when subsequently
wrapped in Φ: `Φ(Φ(XY)·Z) = Φ(XY·Z)` because Φ(XY)=XY as elements of A (their
restriction to M agrees, and A is determined by its action on M via th_idemp_structure).

Therefore:
```
(X★Y)★Z = Φ(XYZ) = X★(Y★Z)
```

exactly. The assoc_defect = 0 when η=0.

### ax_C* (exact, ε=0): From the UCP inequality Φ(X†)Φ(X) ≤ Φ(X†X) (L1692) applied
with X∈A (so Φ(X)=X):
```
‖Φ(X†X)‖ = ‖X†★X‖ ≥ ‖Φ(X†)Φ(X)‖ = ‖X†X‖ = ‖X‖²
```

so `‖X†★X‖ ≥ ‖X‖²`. The C*-identity slack = 0. (Upper bound: from ax_prodnorm
with ε=0, `‖X†★X‖ ≤ ‖X†‖‖X‖ = ‖X‖²`, so equality holds.) cstar_defect = 0.

### ax_prodnorm (exact, ε=0): From `‖Φ‖_cb ≤ 1` (L1718), for X,Y∈A:
```
‖X★Y‖ = ‖Φ(XY)‖ ≤ ‖XY‖ ≤ ‖X‖‖Y‖
```

so `prodnorm_defect = sup (‖X★Y‖/(‖X‖‖Y‖)) − 1 ≤ 0`, i.e. defect = 0.

**Summary for η=0 oracle:** All six estimators (assoc_defect, prodnorm_defect,
cstar_defect, invnorm_defect, invprod_defect, unit_defect × 3) must read ≤
machine-epsilon (exactly 0 in exact arithmetic). The cross-check passes iff all
read ≤ O(double-precision eps) ≈ 1e-15.

---

## 6. Universality canary: why the constant must not grow with dim A

The paper explicitly warns at L484:

```
For finite-dimensional C* algebras, a diagonal can be obtained as
D = ∫ dU (U†⊗U), where the integral is taken with respect to the Haar measure
on the unitary group. Unfortunately, naive constructions of the Haar measure (or
just the diagonal) in the ε-associative setting have error bounds proportional to
n = dim A. So the outlined procedure of fixing the multiplication works only if
ε < cn^{-1} for some constant c.
```

The paper's whole point is that the incremental construction of th_main (L460–462)
circumvents this: "the implicit constant in O(ε) does not depend on A or its
dimensionality." Likewise, in th_almost_idemp (L2192–2193): the O(η) constant is
universal.

**Implication for the defect estimator tests:**

1. Fix η (by taking η-idempotent Φ at a fixed η value, e.g. η=0.01).
2. Vary dim A (e.g. dim A = 4, 9, 16, 25 by choosing Φ on M_N for N=2,3,4,5
   with a rank-N idempotent perturbation).
3. For each dim A, compute assoc_defect / η and prodnorm_defect / η.
4. The ratios must REMAIN BOUNDED as dim A grows. If any ratio grows with dim A,
   the estimator has taken the naive Haar-diagonal route (or the construction is
   wrong).

This is the "universality canary": the normalized defect/η ratio is the test
statistic; it must be O(1) uniformly in dim A.

---

## Implementation implications for ecstar

### Data model

The `ecstar` struct holds:
- The Frobenius-orthonormal basis {B_k}, k=1..d, B_k ∈ M_N.
- The UCP map Φ̃ (as the superoperator, N²×N² matrix) that defines the star product.
- The underlying N.

All defect estimators work with the operator norm (spectral norm, largest singular
value) on M_N. The cb-norm is NOT used inside the estimator — that was needed to
define η. The estimators consume A with its inherited operator norm.

### Defect estimator signatures (prose)

1. **assoc_defect(A, Phi_tilde)** → double
   Computes sup over basis triples of `‖(B_j★B_k)★B_l − B_j★(B_k★B_l)‖ / (‖B_j‖‖B_k‖‖B_l‖)`.
   Norm used: operator/spectral norm (LAPACK `dgesvd` or `zgesvd`).
   η=0 oracle: must return ≤ 1e-14.

2. **prodnorm_defect(A, Phi_tilde)** → double
   Computes `max_{j,k} ‖B_j★B_k‖ / (‖B_j‖‖B_k‖) − 1`, lower-bounded at 0.
   Norm: operator norm.
   η=0 oracle: must return ≤ 1e-14.

3. **cstar_defect(A, Phi_tilde)** → double
   Computes `sup_{X in basis} (1 − ‖X†★X‖ / ‖X‖²)`, lower-bounded at 0.
   Note: X†★X = Φ̃(X†X); X† is the M_N adjoint of X.
   Norm: operator norm.
   η=0 oracle: must return ≤ 1e-14.

4. **invnorm_defect(A)** → double
   Computes `max_k | ‖B_k†‖ − ‖B_k‖ | / ‖B_k‖`.
   Norm: operator norm. Should be 0 by construction (M_N adjoint preserves op-norm).
   η=0 oracle: must return 0.0 (exact).

5. **invprod_defect(A, Phi_tilde)** → double
   Computes `max_{j,k} ‖(B_j★B_k)† − B_k†★B_j†‖ / (‖B_j‖‖B_k‖)`.
   Norm: operator norm. Should be 0 by ax_* exactness.
   η=0 oracle: must return ≤ 1e-14.

6. **unit_right_defect(A, Phi_tilde, I)** → double
   Computes `max_k ‖B_k★I − B_k‖ / ‖B_k‖`.
   Norm: operator norm. I = 1_N (identity matrix).
   η=0 oracle: must return 0.0 (exact).

7. **unit_left_defect(A, Phi_tilde, I)** → double
   Computes `max_k ‖I★B_k − B_k‖ / ‖B_k‖`.
   η=0 oracle: must return 0.0.

8. **unit_norm_defect(A, I)** → double
   Computes `|‖I‖ − 1|` where I = 1_N.
   Should be 0.0 (operator norm of identity = 1 exactly).

### Which norm in which estimator

All eight estimators use the **operator (spectral) norm** on M_N, computed as the
largest singular value. This is the norm the paper uses throughout §2's ε-C*
definition. The cb-norm is NOT used here; it appears only in the definition of η and
in the certification that A is an extended O(η)-C* algebra (the uniformity in n
from th_almost_idemp's "extended" qualifier). The extended property should be checked
separately by a cb-norm test on Φ̃ directly.

### η=0 oracle summary

When Φ is EXACTLY idempotent (η=0), th_idemp_structure (L318–344) and the
Choi–Effros theorem (L344, [ChEf77]) guarantee A=ImgΦ is a genuine C* algebra.
All estimators above must return values ≤ O(machine epsilon):

| estimator        | expected η=0 value | reason                              |
|------------------|-------------------|-------------------------------------|
| assoc_defect     | 0                 | lem_idemp (L1916): Φ(Φ(XY)Z)=Φ(XYZ) |
| prodnorm_defect  | 0                 | ‖Φ‖_cb ≤ 1 implies ‖X★Y‖≤‖X‖‖Y‖    |
| cstar_defect     | 0                 | PhiXdX + Φ(X†X)=X†★X for X∈A        |
| invnorm_defect   | 0                 | M_N op-norm is unitarily invariant   |
| invprod_defect   | 0                 | Φ̃ Hermicity-preserving + (XY)†=Y†X† |
| unit_right/left  | 0                 | Φ(X)=X for X∈A, Φ(1)=1              |
| unit_norm        | 0                 | ‖1_N‖ = 1 exactly                    |

Any nonzero reading in the η=0 oracle is a bug in the estimator, not a
numerical-precision issue (for reasonable N ≤ 20).

### Critical distinction for implementation

The product in all estimators is `X★Y = Φ̃(XY)` — the Choi–Effros star, where `XY`
is the ordinary matrix product in M_N and Φ̃ is applied afterward. Using `XY`
directly (without Φ̃) would exit ImgΦ̃ and give wrong answers. This is the
hallucination-risk at CLAUDE.md §"Domain hallucination-risk callouts" item 2.
