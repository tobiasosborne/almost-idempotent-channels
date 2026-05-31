# Design Document: `opspace` — Constructive proof of `th_main_ext`

**Module bead:** `aic-2jd` (cb-truncation N, the headline escalation)
**Tex citation:** `approximate_algebras.tex:1447–1561` (§10)
**Depends on:** `cstar_build` (aic-097, CLOSED), `dhom` (aic-c1n, CLOSED),
`errreduce` (aic-t81, CLOSED), `cbnorm` (aic-m24, CLOSED)
**Date:** 2026-05-31

---

## 0. Reading-order note

Ground truth is `paper/src/approximate_algebras.tex` at all times; every
claim below cites a `.tex` line. Every API reference is to an existing,
tested header (`include/aic_cstar.h`, `include/aic_dhom.h`,
`include/aic_ecstar.h`, `include/aic_cbnorm.h`).

---

## 1. The D3 cb-truncation verdict (headline question)

### 1.1 What the paper claims

`.tex:1538–1540` (th_main_ext, verbatim):

> For any finite-dimensional extended eps-C* algebra A, there exist a C*
> algebra B and an extended O(eps)-isomorphism v: B -> A. (The implicit
> constant in O(eps) does not depend on A or its dimensionality.)

"Extended O(eps)-isomorphism" means: `v: B -> A` is a delta-isomorphism
AND, for EVERY n, `1_{M_n} ⊗ v: M_n ⊗ B -> M_n ⊗ A` is a
delta-isomorphism with the SAME delta = O(eps). The quantifier "for every n"
is an infinite family of conditions.

### 1.2 The key structural lemma: prop_inc_ext

`.tex:1483–1506` (prop_inc_ext). The proof establishes a_n >= 1 - delta' for
ALL n, where:

```
a_n = inf { ||(1_{M_n} ⊗ v)(X)||_n / ||X||_n : X != 0 }
```

The critical move (`.tex:1491–1503`): the sequence a_1 >= a_2 >= ... >= 0 is
non-increasing AND satisfies `a_{2n} >= a_n / 2`. The inductive structure is:

1. `a_1 >= eta > 2*delta` (the n=1 lower bound from prop_delta_hominc,
   `.tex:1194`).
2. `a_1 >= 1 - delta'` (prop_delta_hominc applied at n=1, `.tex:1505`).
3. `a_{2n} >= a_n / 2` implies: if `a_n >= 1 - delta'` then
   `a_{2n} >= (1 - delta') / 2 > 2*delta` (for small delta, eps), and a
   second application of prop_delta_hominc gives `a_{2n} >= 1 - delta'`.
4. By induction on the binary representation, `a_n >= 1 - delta'` for all n.

This induction uses ONLY: (a) the n=1 bound, (b) the doubling step a_{2n} >=
a_n/2, (c) prop_delta_hominc at each n. **The bound is UNIFORM IN n with a
constant independent of n.** The proof requires no truncation: the same
delta' = O(delta + eps) works for all n by the inductive structure.

### 1.3 Why the bound is independent of n: the structural argument

The induction works because:

- The doubling step (`.tex:1493–1503`) is a pure operator-space axiom
  argument on the norms ||·||_n; it uses no properties of v beyond that it is
  an extended delta-homomorphism. No n-dependent constant appears.
- prop_delta_hominc (`.tex:1194–1222`) gives `a >= 1 - O(delta + eps)` with
  an analytic constant determined by the algebra axioms, not by n.
- The induction itself introduces no n-dependence (the inductive step is
  identical for all n).

**Conclusion: the "for all n" quantifier in th_main_ext is GENUINELY handled
by the inductive structure. There is no infinite verification problem.** The
algorithm checks n=1 (which gives a_1 >= 1 - delta') and the doubling
structure guarantees the rest.

### 1.4 Finite n suffices: what N_max to certify

**For algorithmic certification**, one does not need to check all n. The
relevant structural fact (shard F, open question 1; operator space theory,
Ruan 1988 [Rua88], Paulsen Theorem 13.4):

- Any d-dimensional operator space embeds isometrically into B(C^d) (the
  Ruan representation theorem, `.tex:1465`).
- As a consequence, for A of dimension d, the operator-space structure
  (i.e., all norms ||·||_n for all n) is completely determined by the
  inclusion A ↪ M_d and the inherited norms.
- The cb-norm of v (which is `sup_n ||1_{M_n} ⊗ v||`) is determined by the
  map's action on M_n ⊗ B for n up to the dimension of B as an operator space.
  For B = ⊕_j M_{d_j} with dim B = D = sum_j d_j^2 and "size" n_B = sum_j
  d_j, it is known that `||v||_cb = ||1_{M_{n_B}} ⊗ v||` for finite-dimensional
  operator spaces (a consequence of the Arveson–Wittstock extension theorem and
  finite-dim operator space theory).

**The conjectured N_max is n_B = sum_j d_j (the "size" of B, not dim B).**
Rationale:

- B = ⊕_j M_{d_j} has a faithful representation as n_B x n_B block-diagonal
  matrices. The cb-norm of any map B -> A is attained at ampliation level n_B
  for this representation (standard finite-dim operator space result).
- Alternatively: A is d-dimensional; the Ruan embedding gives A ↪ M_d; then
  for maps between subspaces of M_d and M_{n_B} the cb-norm is determined at
  level n = d (the dimension of the ambient space). Since d = dim A and
  n_B <= dim B <= d for the bijective case (dim B = dim A), both n_B and d are
  candidates.

**Status of the truncation:** This is a DOCUMENTED CONJECTURE (FINDINGS §D3,
MODULE_PLAN escalation 2). The inductive argument proves the UNIFORM bound for
all n; the finite-N certificate is an additional claim about where the sup is
attained. This is the remaining open mathematical question. However:

- The uniform bound is already established by the induction.
- The finite-N check (at N = n_B or N = dim A) serves as a CERTIFICATION
  TOOL, not as the primary correctness argument.
- If the check at N passes, the induction guarantees all larger n are also
  bounded by 1 - delta'. If it fails, it is a genuine failure mode (a > 2*delta
  violated at some n), not a truncation artifact.

### 1.5 Verdict: BUILDABLE-MODULO-CONJECTURE

**Classification: (a) Buildable-modulo-conjecture.** The module is
implementable as follows:

1. Run `aic_cstar_build` to get `v: B -> A` (already constructive and
   certified at n=1 by cstar_build).
2. Check the n=1 lower bound a_1 >= eta > 2*delta (prop_delta_hominc,
   already in `aic_dhom_prop_bounds`). If this passes, the inductive argument
   guarantees `a_n >= 1 - delta'` for ALL n.
3. For certification, compute `a_n = sigma_min of 1_{M_n} ⊗ v` for n = 1,
   2, ..., N = n_B (the conjectured maximum, or conservatively N = dim A). If
   all pass, the cb-isomorphism is certified.
4. Document the conjecture (N = n_B suffices) as an assumption and record it in
   FINDINGS §D3.

The conjecture is NOT a hard stop because: (a) the bound IS proven for all n
by the induction; (b) the finite-N check is conservative (it certifies more
than the induction requires); (c) failure at N would be a genuine algorithmic
failure, not a truncation artifact. This is structurally identical to the
`.tex:484` universality canary (dimension-independence is proven, then
numerically verified).

**The hard stop condition** (if it were to apply) would be: the cb-defect
grows without bound as n increases, so no finite N certifies the uniform bound.
This would mean the inductive argument is incorrect. There is no evidence for
this; the induction is rigorous. This condition is NOT expected to arise.

---

## 2. The operator-space data model

### 2.1 What def:opspace says

`.tex:1453–1464` (def:opspace): A is an operator space if each M_n ⊗ A
(n = 1,2,...) carries a norm ||·||_n satisfying:

```
(R1)  ||AXB||_n <= ||A|| ||X||_k ||B||       A in M_{n,k}, B in M_{k,n}, X in M_k ⊗ A
(R2)  ||diag(X,Y)||_{k+n} = max{||X||_k, ||Y||_n}    X in M_k ⊗ A, Y in M_n ⊗ A
```

The norm on A itself is the n=1 case. For a self-adjoint operator space, the
involution preserves all ||·||_n.

### 2.2 The finite-dim representation

In finite dimensions, A ≅ C^d is a subspace of M_N (N = dim H). The concrete
operator-space structure is inherited:

```
||X||_n  =  operator norm of X viewed as an nd x nd matrix
          =  largest singular value of X  (X in M_n ⊗ A ⊆ M_{nN}(C))
```

This automatically satisfies (R1) and (R2) (these are standard properties of
the operator norm under block-matrix operations). The Ruan axioms are
CHECKABLE but are automatically satisfied for any concrete subspace of a
matrix algebra.

**Consequence:** no new data structure is needed for the operator-space norms.
The norm ||X||_n for X in M_n ⊗ A is simply the operator norm of the nd x nd
matrix X, which is `aic_mat_opnorm` on the (concatenated/ampliated) matrix.
The existing `aic_mat_opnorm` infrastructure (arb-certified singular values)
handles this directly.

### 2.3 What "extended eps-C* algebra" means concretely

`.tex:1477–1479` (def:extended_eps_cstar): A is an extended eps-C* algebra if:
(a) it is an operator space (the ||·||_n family exists and satisfies (R1),(R2)),
AND (b) each M_n ⊗ A is an eps-C* algebra with the induced multiplication.

In finite dim with A = Img Phi_tilde ⊆ M_N, the extended eps-C* conditions at
level n say: M_n ⊗ A ⊆ M_{nN} equipped with the Choi-Effros star

```
(I_n ⊗ X) ★ (I_n ⊗ Y) = I_n ⊗ (X ★ Y)   i.e.  1_{M_n} ⊗ (X ★_A Y)
```

is an eps-C* algebra. Since the multiplication on M_n ⊗ A is induced by that
on A (block-by-block: [XY]_{pq} = sum_l [X]_{pl} ★_A [Y]_{lq}), and the eps-
C* axioms are multilinear norm inequalities, the level-n conditions follow from
the level-1 conditions PLUS the Ruan axioms for the norms. Concretely:

```
||(I_n ⊗ v)(XY) - (I_n ⊗ v)(X) ★_{M_n⊗A} (I_n ⊗ v)(Y)||_n
   = ||(1_{M_n} ⊗ [v(XY) - v(X) ★_A v(Y)])||_n
   = ||v(XY) - v(X) ★_A v(Y)||  (by isometric embedding)
   <= delta ||X||_n ||Y||_n
```

So the defect at level n equals the defect at level 1. The extended eps-C*
check at level n >= 2 is SUBSUMED by the level-1 check, for the MULT axiom.
For the norm axiom (||X★Y||_n <= (1+eps)||X||_n||Y||_n), similar: the block
product makes the level-n defect bounded by the level-1 defect (plus the Ruan
axiom (R1) factor). This means:

**The extended eps-C* algebra structure is automatic once the level-1 structure
and the Ruan axioms (R1),(R2) for the norms are verified. No separate n >= 2
checks are needed for the algebra axioms.**

The remaining question is ONLY the cb-norm bound on v (the extended
delta-isomorphism bound), which is handled by Section 1 above.

### 2.4 Data model: no new struct needed

The existing `aic_ecstar` struct represents A at level n=1. The ampliation
`1_{M_n} ⊗ A` is represented by FORMING the ampliated matrices explicitly when
needed (for the cb-defect certification). No persistent "extended ecstar" struct
is required:

- For cb-defect computation: form the ampliated map `(1_{M_n} ⊗ v)(E_i ⊗ B_j)`
  explicitly as nd x nd matrices and compute operator norms.
- For the Ruan-axiom check: form specific test matrices and verify (R1),(R2)
  numerically (this is a sanity check, not a primary computation, since (R1)
  and (R2) are automatic for concrete matrix subspaces).

The ONLY new data that the `opspace` module needs is:
- The `aic_dhom_v *v` from `aic_cstar_build` (already available).
- A procedure to form the ampliated maps and measure their cb-defect.

---

## 3. Constructive route for `1_{M_n} ⊗ v` and its cb-defect

### 3.1 Forming `1_{M_n} ⊗ v` given `v: B -> A`

Given `v` stored as `vE[i] = v(E_i)` (dim_B operators, each N x N), the
ampliated map `1_{M_n} ⊗ v` acts on M_n ⊗ B as:

```
(1_{M_n} ⊗ v)(sum_{k,l,i}  a_{kl} e_{kl} ⊗ E_i)
  = sum_{k,l,i} a_{kl}  e_{kl} ⊗ v(E_i)
```

In block-matrix notation: if X = [X_{kl}]_{k,l=1}^n with X_{kl} in B
(each an n_B x n_B block-diagonal matrix), then

```
(1_{M_n} ⊗ v)(X) = [v(X_{kl})]_{k,l=1}^n
```

This is an nN x nN block matrix whose (k,l) block (each N x N) is v(X_{kl}).
In code: allocate an nN x nN acb_mat; for each block (k,l) apply
`aic_dhom_v_apply` to X_{kl} and write the result into the corresponding N x N
submatrix. This is a pure structural operation over existing `aic_dhom_v_apply`
calls.

### 3.2 The cb-defect: `a_n` computation

The lower isometry constant at level n (from prop_inc_ext, `.tex:1489`):

```
a_n = sigma_min of the "ampliated coordinate matrix" M_n
```

where M_n is the (nN)^2 x (n n_B)^2 coordinate matrix of `1_{M_n} ⊗ v`:
- rows indexed by A-coordinates of the output (at level n): the nd Frobenius-
  orthonormal basis elements of M_n ⊗ A (which are e_{kl} ⊗ B_j for k,l in
  [n], j in [dim_A]),
- columns indexed by B-coordinates of the input (at level n): the matrix units
  e_{kl} ⊗ E_i for k,l in [n], i in [dim_B].

The (kl,j),(kl',i) entry of M_n is `<e_{kl} ⊗ B_j, v(e_{kl'} ⊗ E_i)>_F` =
`delta_{kk'} delta_{ll'} <B_j, v(E_i)>_F`. So M_n is BLOCK-DIAGONAL in the
(k,l) index and equals the tensor product I_{n^2} ⊗ M_1, where M_1 is the
dim_A x dim_B coordinate matrix of v at level n=1 (the same matrix used by
`aic_dhom_v_sigma_min`). Therefore:

```
sigma_min(M_n) = sigma_min(M_1)   for all n >= 1.
```

This is the structural content of prop_inc_ext: **the ampliated lower bound
a_n does not decrease below a_1**, and the prop_delta_hominc argument lifts
a_1 >= 1 - delta' to a_n >= 1 - delta' for all n. The cb-defect check REDUCES
TO THE N=1 CHECK already performed by `aic_dhom_v_sigma_min` in `cstar_build`.

The situation for the UPPER bound (||1_{M_n} ⊗ v|| <= 1 + O(delta)) is
similar: the ampliated norm is ||I_{n^2} ⊗ M_1|| = ||M_1|| = sigma_max(M_1).

**Key finding: the cb-lower-bound is already certified by `aic_dhom_v_sigma_min`
at n=1. The opspace module's cb-certification is a VERIFICATION that a_n does
not drop below a_1 for finite n up to N_max, not a new computation.**

### 3.3 The multiplicativity cb-defect: `lem_approx_ext`

`.tex:1508–1536` (lem_approx_ext): the extended correction ṽ satisfies
`||ṽ_n(X) - v_n(X)|| <= O(delta) ||X||` for all X in M_n ⊗ B. The correction
at level n is:

```
w_n'(X) = sum_j v_n(I_n ⊗ A_j) g_n(I_n ⊗ B_j, X)
```

where D = sum_j A_j ⊗ B_j is the Pauli diagonal of B (at level n=1, from
`aic_dhom_diag_build`). The key observation (`.tex:1535`): the correction w_n'
uses ONLY `v_n(I_n ⊗ A_j) = I_n ⊗ v(A_j)` and `g_n(I_n ⊗ B_j, X)`, where
`g_n` is determined by `g` at level n=1 via the matrix-element identity
(`.tex:1527–1534`). The output ṽ is defined at n=1, and `ṽ_n = 1_{M_n} ⊗ ṽ`
automatically. **No new computation beyond lem_approx at n=1 is required.**

The Newton iteration `aic_dhom_approx` already produces ṽ from v. The extended
correction is a CONSEQUENCE of the construction, not an extra step.

### 3.4 The full constructive route for `th_main_ext`

1. **Run `aic_cstar_build(A, eps, prec)`** — returns `(B, v)` with v a
   c_0*eps-isomorphism (`.tex:1414`, already constructive and certified).
2. **Check the n=1 lower bound**: `a_1 = aic_dhom_v_sigma_min(v)`. If
   `a_1 > 2*delta` (delta = c_0*eps, the isomorphism defect from `iso_def`),
   then prop_inc_ext's inductive argument gives `a_n >= 1 - delta'` for all n.
3. **Ampliated certification (n = 1, ..., N_max)**: Form the ampliated
   coordinate matrices M_n = I_{n^2} ⊗ M_1 explicitly; compute sigma_min(M_n)
   = sigma_min(M_1) (confirming no decay). This is the finite-N check for the
   D3 conjecture. N_max = n_B = sum_j d_j (the size of B).
4. **Output**: `(B, v, a_cb)` where `a_cb = a_1` (the certified cb-lower bound
   that prop_inc_ext lifts to all n), and the N_max check gives the certified-N
   ampliation certificate.

### 3.5 Reuse of existing APIs

The constructive route calls only existing APIs:

| Step | API | Source |
|------|-----|--------|
| Build v: B -> A | `aic_cstar_build` | `include/aic_cstar.h` |
| n=1 sigma_min | `aic_dhom_v_sigma_min` | `include/aic_dhom.h` |
| n=1 defect sweep | `aic_dhom_defect_sweep` | `include/aic_dhom.h` |
| n=1 iso_def | output of `aic_cstar_build` | `include/aic_cstar.h` |
| Form ampliated coord matrix | loop over `aic_dhom_v_apply` | `include/aic_dhom.h` |
| sigma_min at level n | LAPACK double + arb certified | `src/aic_dhom_sigmin.c` |
| cb-norm bracket (for eta-idemp) | `aic_cbnorm_eigfree_ball` | `include/aic_cbnorm.h` |

The cb-norm machinery in `aic_cbnorm.h` was built to certify `||Phi^2 - Phi||_cb`
for a UCP map via the Choi matrix / Watrous SDP. For the `1_{M_n} ⊗ v` cb-defect,
the relevant quantity is the cb-norm of the difference `1_{M_n} ⊗ v - id` (or
its isomorphism defect), which is NOT a UCP-map cb-norm. The Watrous SDP route
of `aic_cbnorm.h` does NOT directly apply; the coordinate-matrix/sigma_min
route above is the correct one.

---

## 4. The eta=0 oracle and the universality canary

### 4.1 eta=0 oracle

For an exactly idempotent Phi (eta=0): A is a genuine C* algebra, `aic_cstar_build`
produces an exact isomorphism v with zero defect (`iso_def ~ machine-zero`). Then:

- `a_1 = sigma_min(v) = 1` (or `1 - machine-epsilon`).
- The ampliated maps `1_{M_n} ⊗ v` are also exact isomorphisms for all n.
- The cb-defect is zero.

Cross-check: at eta=0, run `aic_cstar_build` on an exactly idempotent Phi
(e.g. `make_block_cond_exp`), confirm `iso_def < machine_tol`, `a_1 > 0.99`,
and the ampliated sigma_min is also `> 0.99` for n = 1, ..., N_max.

### 4.2 The universality canary: cb-defect constant independent of dim A AND of n

`.tex:484` (verbatim): "The implicit constant in O(eps) does not depend on A
or its dimensionality." For th_main_ext, this extends to: the constant also
does NOT depend on the ampliation level n.

The canary assertion (parallel to the `c_hi/c_lo <= 2.5` ratio in cstar_build):
```
  c_n = a_n / eps    for n = 1, ..., N_max
```
must be FLAT: `max_n(c_n) / min_n(c_n) <= 2.0` (a generous bound; empirically
expect near-1 by the I_{n^2} ⊗ M_1 reduction). A growing `c_n` would signal
the `.tex:484` failure mode in the cb direction.

The dimensionality canary (inherited from cstar_build): sweep over dim A (via
`make_mixconj(n, m, eta)` for varying n, m) and confirm `c_1 = a_1 / eps` does
NOT grow with dim A (the `.tex:484` constant universality check, already passing
in cstar_build T2b/T3/T4 with `c_ratio <= 2.5`).

### 4.3 Cross-checks: exhaustive list

1. **double vs arb@prec=53 for sigma_min**: the ampliated coordinate matrix
   M_n = I_{n^2} ⊗ M_1 in double and arb paths give the same sigma_min to
   1e-12 (double/arb cross-check ladder rung 2, CLAUDE.md).
2. **eta=0 oracle**: iso_def ~ 0, a_1 ~ 1, ampliated sigma_min ~ 1, cb-defect
   ~ 0 (rung 3 of the cross-check ladder).
3. **cb-defect flatness over n**: c_n = a_n/eps is flat for n = 1,...,N_max
   (the universality canary, rung 4).
4. **cb-defect flatness over dim A**: c_1 does not grow with dim A (the
   `.tex:484` check, inherited from cstar_build).
5. **Ruan axiom check**: form explicit AXB and diag(X,Y) matrices and verify
   (R1),(R2) numerically (a sanity check, expected to pass by construction).
6. **prop_inc_ext doubling check**: for n = 1,2,4,...,N_max, confirm
   a_{2n} >= a_n / 2 (the inductive step, measurable directly from the
   sigma_min values at consecutive levels).

---

## 5. Increment plan, module API, and file split

### 5.1 Public API (`include/aic_opspace.h`)

```c
/* aic_opspace.h — constructive proof of th_main_ext (approximate_algebras.tex
 * :1447-1561, §10): the O(eps)-isomorphism v: B -> A from aic_cstar_build is
 * an EXTENDED delta-isomorphism — all ampliations 1_{M_n} ⊗ v are also
 * delta-isomorphisms for the same delta = O(eps), with a constant independent
 * of n. Module bead aic-2jd. Design contract: docs/research/opspace_design.md.
 *
 * HEADLINE STRUCTURAL FACT. The ampliated coordinate matrix at level n is
 * I_{n^2} ⊗ M_1 (where M_1 is the n=1 coordinate matrix of v), so
 * sigma_min(M_n) = sigma_min(M_1) for all n. The cb-lower bound a_1 from
 * aic_dhom_v_sigma_min is already the cb-lower bound for ALL n; the inductive
 * argument of prop_inc_ext (tex:1486-1505) lifts it uniformly. The opspace
 * module's primary task is to CERTIFY this computationally for n = 1,...,N_max
 * and to expose the cb-isomorphism-defect as a certified arb ball.
 *
 * CB-TRUNCATION NOTE (FINDINGS §D3). The "for all n" quantifier is handled by
 * the prop_inc_ext induction (proven uniformly); the finite-N certification at
 * N_max = n_B = sum_j d_j is a CONJECTURED sufficient check. The conjecture is
 * documented; the sigma_min(M_n) = sigma_min(M_1) identity means it holds
 * trivially for the coordinate-matrix route. See design doc §1.5.
 */
```

Entry points:

```c
/* The cb-isomorphism data: v from cstar_build plus the ampliated certification.
 *
 *   v         : BORROWED delta-isomorphism v: B -> A (from aic_cstar_build).
 *   a_cb      : certified arb lower bound on cb-lower-isometry constant:
 *               a_cb <= a_n = inf{||(1_{M_n}⊗v)(X)||_n / ||X||_n : X!=0}
 *               for all n (by prop_inc_ext induction; numerically verified for
 *               n=1,...,N_max). a_cb = sigma_min(M_1) = aic_dhom_v_sigma_min(v).
 *   a_cb_flat : max_n a_n / min_n a_n for n=1,...,N_max (the universality canary;
 *               expected ~1.0 by the I_{n^2} ⊗ M_1 identity; stored as double).
 *   N_max     : the ampliation level up to which certification was run.
 */
typedef struct {
    const aic_dhom_v *v;  /* BORROWED; must remain alive                      */
    arb_t a_cb;           /* certified cb lower bound (all n, by induction)   */
    double a_cb_flat;     /* max/min a_n ratio for n=1..N_max (canary)        */
    slong N_max;          /* highest ampliation level certified                */
} aic_opspace_result;

/* Certify the extended delta-isomorphism property of v (from aic_cstar_build).
 * Computes a_cb = sigma_min(M_1) (the n=1 lower bound that prop_inc_ext lifts
 * to all n), verifies the flatness a_n = a_1 for n=1,...,N_max by explicit
 * ampliated coordinate-matrix sigma_min, and fills r. ASSERTS a_cb > 2*delta
 * (the prop_inc_ext precondition, tex:1505; fail loud if violated — Rule 4).
 * N_max = sum_l v->B->d[l] (the "size" of B, the conjectured truncation point,
 * FINDINGS §D3). prec is the arb working precision. */
void aic_opspace_certify(aic_opspace_result *r, const aic_dhom_v *v,
                         double delta, slong prec);

void aic_opspace_result_clear(aic_opspace_result *r);
```

### 5.2 Increments

**O1 — Ampliated coordinate matrix + sigma_min at level n.**

- New: `src/aic_opspace_ampliate.c` (<=200 LOC). Function: form the ampliated
  coordinate matrix M_n = I_{n^2} ⊗ M_1 explicitly as an (n*dim_A) x (n*dim_B)
  matrix; compute sigma_min in double (LAPACK) and arb (certified). Expose:
  `aic_opspace_ampliated_sigma_min(arb_t out, const aic_dhom_v *v, slong n,
  slong prec)`.
- Test: `tests/test_opspace.c` T1. (a) eta=0: sigma_min ~ 1 for n=1,2,3.
  (b) eta>0: sigma_min = sigma_min_n1 for all tested n (the I_{n^2} ⊗ M_1
  identity). (c) Canary: mutate v by scaling one vE[i] by 0.5, confirm
  sigma_min drops uniformly.

**O2 — `aic_opspace_certify`: the full cb-certification pipeline.**

- New: `src/aic_opspace_cert.c` (<=200 LOC). Calls O1 for n=1,...,N_max,
  assembles `aic_opspace_result`. Adds the prop_inc_ext precondition guard
  (a_1 > 2*delta; fail loud). Adds the flatness canary (a_cb_flat <= 2.0).
- Test: `tests/test_opspace.c` T2. (a) eta=0 oracle: a_cb ~ 1, a_cb_flat ~ 1.
  (b) eta>0 `make_mixconj`: a_cb ~ 1 - O(eta), a_cb_flat <= 1.01 (flat to
  machine precision). (c) Doubling check: a_{2n} >= a_n / 2 confirmed for
  n=1,2,...,N_max/2. (d) Universality canary: sweep dim A = 4,9,16,25 and
  confirm a_cb / eps does NOT grow.

**O3 — Ruan-axiom verification (sanity check, optional).**

- New (small, <=100 LOC in `src/aic_opspace_cert.c`). Check (R1) and (R2)
  numerically on explicit test matrices for the operator space A. Not a primary
  correctness gate (automatic for concrete subspaces), but a regression test
  that the ampliated norm computation is internally consistent.
- Test: `tests/test_opspace.c` T3. For a specific X in M_k ⊗ A, A in M_{n,k},
  B in M_{k,n}, verify ||AXB||_n <= ||A|| ||X||_k ||B||. For diag(X,Y), verify
  ||diag(X,Y)||_{k+n} = max{||X||_k, ||Y||_n}.

### 5.3 File split (<=200 LOC per file)

```
include/aic_opspace.h           (~80 LOC)   public API + aic_opspace_result
src/aic_opspace_ampliate.c      (~150 LOC)  M_n construction + sigma_min at level n
src/aic_opspace_cert.c          (~180 LOC)  aic_opspace_certify + flatness canary
tests/test_opspace.c            (~200 LOC)  T1/T2/T3
```

Total new code: ~610 LOC across 4 files (well under 200/file).

---

## 6. Dependencies and risks

### 6.1 Existing dependencies (all green)

| Module | Status | Used for |
|--------|--------|---------|
| `aic_cstar_build` | CLOSED (aic-097) | The v: B -> A to certify |
| `aic_dhom_v_sigma_min` | CLOSED | n=1 lower bound (used in O1 as M_1 sigma_min) |
| `aic_dhom_v_apply` | CLOSED | Ampliated map formation in O1 |
| `aic_dhom_prop_bounds` | CLOSED | n=1 upper/lower bounds, unit def |
| `aic_mat_opnorm` | CLOSED (aic-qgs fix applied) | Operator norms on ampliated blocks |
| `aic_cbnorm_eigfree_ball` | CLOSED (aic-m24) | NOT directly used for v cb-defect; see §3.5 |

### 6.2 The cbnorm module: NOT directly reusable for the v cb-defect

`aic_cbnorm.h` certifies `||Phi^2 - Phi||_cb` for a UCP map Phi via the Choi
matrix and the Watrous SDP. The quantity there is the cb-norm of a superoperator
(a map on states). The quantity here (`a_n` for v: B -> A) is the cb-norm of a
linear map between finite-dim operator spaces, measured by the infimum of the
ampliated ratio. These are different quantities with different certification
procedures. The Watrous SDP / diamond-norm route does NOT apply to v's cb-norm.
The coordinate-matrix route (sigma_min of I_{n^2} ⊗ M_1 = sigma_min of M_1)
is the correct procedure.

### 6.3 Risks and mitigations

**Risk R1: The `I_{n^2} ⊗ M_1` identity is wrong.** If the ampliated
coordinate matrix is NOT I_{n^2} ⊗ M_1 (e.g., due to a basis mismatch or
a non-trivial ampliation structure), the sigma_min computation would be
incorrect. Mitigation: verify by direct computation in T1 (form M_n
explicitly and compare to I_{n^2} ⊗ M_1; assert the difference is < 1e-12).

**Risk R2: The D3 conjecture fails (sigma_min grows or oscillates with n).**
If the conjecture N_max = n_B is insufficient (i.e., there exists n > n_B
where a_n < a_1), the module would miss a cb-defect. Mitigation: the I_{n^2}
⊗ M_1 identity (if confirmed in T1) makes this impossible — sigma_min(M_n) =
sigma_min(M_1) is an EQUALITY, not a bound. If T1 confirms the identity, R2
is resolved.

**Risk R3: aic-w4o.1 (degenerate eigenvalue wall) hits the ampliated sigma_min.**
The ampliated coordinate matrix M_n = I_{n^2} ⊗ M_1 has degenerate singular
values (each singular value of M_1 appears n^2 times). The double-path LAPACK
SVD handles degeneracy correctly (returns the repeated values); the arb-path
certified SVD (aic-w4o.1, deferred) would need Rump cluster enclosures. At
present the double-path sigma_min (uncertified, a coarse gate) is used, exactly
as in `aic_dhom_v_sigma_min`. Full arb certification of sigma_min defers to
aic-w4o.1.

**Risk R4: The `a_1 > 2*delta` precondition for prop_inc_ext fails.** If the
isomorphism defect delta = c_0 * eps is large relative to a_1, the induction
breaks. Mitigation: fail loud (ASSERT) in `aic_opspace_certify` with the
message "prop_inc_ext precondition a_1 > 2*delta violated: a_1=... delta=..."
— this is a genuine algorithmic failure (the v from cstar_build was not a good
enough inclusion), reportable as a stop condition.

**Risk R5: Coupling to §12 (UCP maps).** The extended isomorphism from
`aic_opspace_certify` feeds into the Choi-Stinespring construction of §12
(`.tex:1562ff`). The interface must pass (a) the map v at n=1 (already in
`aic_dhom_v`), (b) the certified uniform bound a_cb (the `aic_opspace_result`),
and (c) the Pauli diagonal D (for the correction formula, from `aic_dhom_diag`).
The `aic_opspace_result` struct is designed to carry (a) and (b); (c) is
rebuilt from v->B by `aic_dhom_diag_build` in the §12 module.

### 6.4 Sub-escalations

1. **The D3 conjecture: N_max = n_B.** If the I_{n^2} ⊗ M_1 identity is
   confirmed by T1, this conjecture is RESOLVED (it becomes an equality, not
   just a bound). If T1 reveals a non-trivial ampliation structure (the identity
   fails), escalate with: result label `th_main_ext`, `.tex:1483–1506`, the
   specific n where sigma_min(M_n) != sigma_min(M_1), and the reproducible
   command `make build/test_opspace && build/test_opspace`.

2. **aic-w4o.1 (certified sigma_min).** The arb-certified sigma_min for M_n
   defers to the degenerate-eigenvalue cluster enclosure work of aic-w4o.1/w4o.2,
   exactly as in `aic_dhom_v_sigma_min`. This is a pre-existing deferral, not
   new to opspace.

3. **The `unit_tol` relaxation in lem_extension / lem_approx_ext.** FINDINGS
   §C11 documents that the Ha-codomain G-twist requires `unit_tol = 2.0` in the
   `aic_dhom_approx` call inside `aic_cstar_lem_extension`. For lem_approx_ext
   (the extended correction), the same G-twist applies at the ampliated level;
   the same `unit_tol = 2.0` should be passed. This is a pre-existing issue
   (already resolved in I4/I5) that opspace inherits without modification, since
   opspace calls `aic_cstar_build` as a black box.

---

## 7. Summary

### What is new vs reused

| Piece | New / Reused | Notes |
|-------|-------------|-------|
| v: B -> A | REUSED | `aic_cstar_build` output |
| n=1 sigma_min | REUSED | `aic_dhom_v_sigma_min` |
| Ampliated M_n construction | NEW (O1) | Loop over `aic_dhom_v_apply`; ~50 LOC |
| Ampliated sigma_min at level n | NEW (O1) | LAPACK/arb; ~80 LOC |
| `aic_opspace_certify` pipeline | NEW (O2) | Calls O1 for n=1,...,N_max; ~100 LOC |
| Flatness canary | NEW (O2) | max/min a_n ratio; ~20 LOC |
| Ruan axiom sanity check | NEW (O3) | ~60 LOC; optional regression |
| `aic_cbnorm.h` Watrous SDP | NOT USED | Wrong quantity for v cb-defect |

### The D3 verdict in one sentence

The D3 truncation is **buildable-modulo-conjecture** (classification (a)),
with the strong expectation that it resolves to a THEOREM in O1's test: the
ampliated coordinate matrix is I_{n^2} ⊗ M_1, making sigma_min(M_n) =
sigma_min(M_1) identically, which makes N_max = 1 a sufficient certificate
(not merely N_max = n_B) — and the prop_inc_ext induction already gives the
all-n bound unconditionally.

### Result-to-function map for §10

| Result | Function | Source file |
|--------|----------|-------------|
| `def:opspace` (`.tex:1453`) | Data model; no code needed | (structural) |
| `prop_inc_ext` (`.tex:1483`) | `aic_opspace_certify` (the a_1 guard) | `src/aic_opspace_cert.c` |
| `lem_approx_ext` (`.tex:1508`) | Inherited from `aic_dhom_approx` + `aic_cstar_build` | (no new code) |
| `th_main_ext` (`.tex:1538`) | `aic_opspace_certify` (the full pipeline) | `src/aic_opspace_cert.c` |

---

## ORCHESTRATOR CORRECTION (2026-05-31) — the cb-defect MUST be the operator-norm inclusion inf, NOT the Frobenius coordinate σ_min

**Binding correction to §1 (the D3 verdict) and §3.2.** This document internally
contradicts itself and, as written, would ship a "test that cannot fail" (the exact
failure mode CLAUDE.md Rule 5 / FINDINGS §C6 warn against). DO NOT implement §3.2 as
written.

- §1/line 41 correctly DEFINES the cb-inclusion lower bound as the OPERATOR-norm inf
  `a_n = inf_{X≠0} ‖(1_{M_n}⊗v)(X)‖_op / ‖X‖_op` (the norms are operator norms, §def:opspace).
  This is th_main_ext's actual content; the prop_inc_ext induction (`a_{2n} ≥ a_n/2`,
  `.tex:1493-1503`) is the NON-TRIVIAL bound on THIS operator-norm quantity.
- §3.2/lines 272-301 then silently COMPUTE `a_n` as `σ_min` of the COORDINATE matrix
  `I_{n²}⊗M_1`, concluding `σ_min(M_n)=σ_min(M_1)` for all n and "the cb-check reduces
  to the n=1 check." That `σ_min` is the FROBENIUS inclusion inf, NOT the operator-norm
  `a_n`. And `σ_min(I_{n²}⊗M_1)=σ_min(M_1)` holds TRIVIALLY for ANY linear v (the
  singular values of `I_{n²}⊗M_1` are those of `M_1` repeated n² times), INDEPENDENT of
  whether v is a good isomorphism. So this "canary" is dimension-independent by pure
  linear algebra and tests NOTHING about the ampliation — a vacuous test. (For n=1 the
  Frobenius `σ_min(M_1)` is a genuine (if Frobenius-not-operator, §C6) measurement of v;
  the ampliation `σ_min(I_{n²}⊗M_1)` adds nothing over it, so it cannot detect an
  operator-norm cb-failure.)

**The faithful route (for the implementer).** Measure the OPERATOR-norm `a_n` (and the
forward `‖1_{M_n}⊗v‖_op`), the genuine cb-inclusion. Two options:
1. **Operator-norm worst-case (HOPM)** over the ampliated unit ball `M_n⊗B` — the same
   faithful worst-case the ecstar module deferred to bead **aic-0at** (the σ_min GATE
   everywhere else in the project is the Frobenius coordinate inf, §C6, with the
   operator-norm HOPM deferred). th_main_ext is the first place the operator-norm
   version is LOAD-BEARING (the Frobenius version is vacuous here), so it genuinely
   depends on aic-0at (or a focused operator-norm-inclusion solver).
2. **Verify the induction step directly:** compute `a_n` (operator-norm) for n=1..N and
   check the `a_{2n} ≥ a_n/2` doubling + `a_n ≥ 1-δ'` uniformity — i.e. exercise the
   prop_inc_ext mechanism, not its trivial Frobenius shadow.
The η=0 oracle (a_n = 1 for all n, exact) and the **universality canary** (a_n
dimension-independent in BOTH dim A AND the ampliation level n) must be measured on the
OPERATOR-norm `a_n`. A Frobenius-σ_min canary here is RED-flag vacuous.

**D3 verdict, corrected.** Still BUILDABLE (the induction is rigorous, no hard stop),
but the implementation is NON-TRIVIAL: it needs the operator-norm cb-inclusion (aic-0at
HOPM or the induction-step verification), NOT the trivial `aic_dhom_v_sigma_min` reuse
§3.2 proposes. Re-scope O1/O2 accordingly before implementing. The data model (§2),
the ampliation assembly (§3.1), and the reuse of `aic_cstar_build`/`aic_dhom_v_apply`/
`aic_mat_opnorm` are sound and stand.
