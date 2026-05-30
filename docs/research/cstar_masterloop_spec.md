# Spec: I5 — Master Loop (`aic_cstar_build`), proof of `th_main`

**Module bead:** `aic-097`
**Tex citation:** `approximate_algebras.tex:1414-1444` (proof), `:460-491` (statement),
`:1317-1319` (cor_improvement / c_0), `:1428` (eps' + equivalence classes)
**Depends on:** I1-I4 (all CLOSED + green), aic-qgs substrate (CLOSED: unblocks oblique corner path)
**Implements:** `aic_cstar_build` in `src/aic_cstar_build.c` (entry point)
**Date:** 2026-05-31

---

## 0. Read-order note

This spec is authoritative for I5. All `.tex` citations are to
`paper/src/approximate_algebras.tex`. Every design §-reference is to
`docs/research/cstar_build_design.md`. Every finding reference is to
`paper/FINDINGS.md §Xn`.

---

## 1. Faithfulness verdict: does the proposed design match the `.tex:1414-1444` proof?

**VERDICT: YES, with the clarifications below.**

The design's three-stage decomposition is a faithful constructivization of the
paper's proof. Key matches:

- Stage 1 (greedy loop) = `.tex:1417-1426` (the maximum-dimensionality
  contradiction argument).
- Stage 2 (inductive extension per class) = `.tex:1430-1441` (the M_{r-1} -> M_r
  induction with compression, lem_extension, and cor_improvement).
- Stage 3 (merge classes) = `.tex:1443` (cor_merge_sum over classes + cor_improvement).

**Key moves in the proof, quoted verbatim for alignment:**

`.tex:1417-1426` — termination argument: "Let v_comm : B_comm -> A be a
c_0*eps-inclusion of a commutative C* algebra of maximum dimensionality...
Without loss of generality, we may assume that the last projection fails to be
one-dimensional, i.e. dim S_{P_m} > 1... Merging v_comm^(1) and v_comm^(2) using
Corollary cor_merge_sum, we obtain an O(eps)-inclusion v_comm^+ : B_comm^+ -> A...
But this contradicts the maximum dimensionality assumption."

`.tex:1428` — setup for Stages 2-3: "Let eps' = O(eps) denote a suitable upper
bound such that eps' >= c_0*eps, and moreover, S_P for any c_0*eps-projection P
is an eps'-C* algebra. The indices j of the one-dimensional projections P_j are
subdivided into equivalence classes such that dim S_{P_j,P_k} is 1 if j and k
are in the same class and 0 otherwise."

`.tex:1435-1441` — Stage 2 induction: "To obtain v_r from v_{r-1}, we first fit
v_{r-1}, P_{[1,r-1]}, and P_r into A_r by applying the compression map
Co_{P_{[1,r]}}... Then we extend the O(eps)-isomorphism v_{r-1}' to v_{r-1}^+ using
Lemma lem_extension. Finally, we use Corollary cor_improvement to replace v_{r-1}^+
with a c_0*eps'-isomorphism v_r."

`.tex:1443` — Stage 3: "Note that by Lemma lem_add_dim, S_{P_C,P_D}=0 if the
classes C and D are distinct. This allows us to successively merge the isomorphisms
v_C for different C."

**One design clarification (Stage 1 initialization):** The design (§5.1) says
"start with v_comm: M_1 -> A (the trivial inclusion of the unit)." The paper says
"let v_comm be a c_0*eps-inclusion of a commutative C* algebra of maximum
dimensionality" (existence by the greedy-termination argument itself). The constructive
substitute — starting from a single nontrivial projection P found by
`aic_projection_nontrivial(A)` and initializing {P_1=P, P_2=Ptilde-P} — is
faithful: the initial c_0*eps-inclusion has B_comm = C (+) C (dimension 2) and
v_comm(Pi_j) = P_j (the two projections). See Stage 1 spec below for the exact
initialization.

**No discrepancies requiring a stop condition were found.**

---

## 2. Stage 1 — greedy commutative skeleton loop

### 2.1 Initialization

**The first call is to `aic_projection_nontrivial(A)`** (not the "trivial M_1 inclusion").

The paper's "maximum dimensionality" starting point is the endpoint of the loop;
constructively we initialize with the FIRST nontrivial split:

```
IF dim A == 1: FAIL LOUD (no nontrivial projection; th_main is vacuous, B = C)
// This also catches the already-scalar case.

// Find a nontrivial projection in A itself:
P_1 = aic_projection_nontrivial(A, ...)     // .tex:931; returns delta-proj, O(delta) defect
P_2 = Ptilde_A - P_1                        // Ptilde_A = unit of A (its aic_ecstar unit field)

// Apply cor_improvement to the initial 2-projection inclusion:
v_comm = aic_errreduce({P_1, P_2} as a v : C(+)C -> A)
// This normalizes the projections to a c_0*eps-inclusion.

// Initial list: {P_1, P_2}  (m = 2 projections, indices 1..2)
```

The initial v_comm maps Pi_1 -> P_1, Pi_2 -> P_2 where Pi_1, Pi_2 are the
standard basis projections of C (+) C = C^2. Its vE[] are {P_1, P_2} (two
operators; aic_dhom_B has num_blocks=2, d[0]=d[1]=1, dim_B=2).

After `aic_errreduce`, the two projections satisfy the c_0*eps-inclusion
bounds: mult_def <= c_0*eps, sigma_min >= 1-c_0*eps, unit_def <= c_0*eps.

**Degenerate case dim A = 1:** B = C (a single 1-dim block), v : C -> A,
v(1) = Ptilde_A (the unit of A). This is an EXACT isomorphism (zero defect,
since A = C * Ptilde_A). No projection splitting needed. The loop body is
skipped; proceed directly to Stage 2 (which is vacuous for a 1-dim class) and
Stage 3 (a single class of size 1 needs no merge).

### 2.2 Loop body

```
WHILE any j has dim S_{P_j} > 1:
    find m such that dim S_{P_m} > 1    // aic_corner_dim_S(A, P_m, P_m, prec) > 1
    build wrapper: aic_cstar_subalg_build(&sp, A, P_m, prec)
    // sp.A is the S_{P_m} eps'-C* algebra
    
    // Find nontrivial projection in S_{P_m}:
    aic_projection_nontrivial(&P_prime, &delta_p, &Pnorm, &ImPnorm, &sp.A, NULL, prec)
    // P_prime is in S_{P_m}, an O(eps')-projection
    
    // Nontriviality CHECK in S_{P_m} (FINDINGS §G1 / §C7 / design §2.3):
    // aic_projection_nontrivial's internal guard uses ||1_n - P||_op (ambient),
    // which is VACUOUS for an S_{P_m} projection that satisfies ||1_n - P||_op ~ 1
    // trivially. The CORRECT check is in S_{P_m}'s unit Ptilde_m:
    ASSERT arb_ge_d(ImPnorm, 0.15)     // the header's own >= 0.3 bound applies
                                        // to the S_P unit, not the ambient 1_n;
                                        // ImPnorm = ||I_A - P||_op reported by the
                                        // projector uses the ambient 1_n BUT
                                        // Ptilde_m is stored in sp.Ptilde, so we
                                        // additionally check:
    // Compute ||sp.Ptilde - P_prime||_op (mid+radius UB, aic_corner_gamma_opnorm_ub)
    // ASSERT that UB >= 0.15 (so P_prime != Ptilde_m in S_{P_m})
    // (This is the §G1 route (b): cstar_build-level nontriviality verify.)
    
    P_double_prime = sp.Ptilde - P_prime  // Ptilde_m - P' (n x n operator)
    // Note: Ptilde_m = aic_corner_Ptilde(A, P_m) stored in sp.Ptilde
    
    // Build two delta-inclusions:
    // v_1: C (+)_{j<m} M_1 -> S_{P_{[1,m-1]}} via the existing projections
    // v_2: C (+) C -> S_{P_m} via P_prime, P_double_prime
    //
    // Concretely: build the new flat projection list
    // P_new = {P_1, ..., P_{m-1}, P_prime, P_double_prime, P_{m+1}, ..., P_k}
    // (replacing P_m with the two new projections)
    
    // Merge via cor_merge_sum + errreduce:
    // Stage 1 uses cor_merge_sum to concatenate:
    //   v_1 = current v_comm restricted to indices != m (P_[1,m-1] + P_[m+1,k] part)
    //   v_2 = trivial C(+)C inclusion v_2(Pi') = P_prime, v_2(Pi'') = P_double_prime
    // The merged B_comm = B_1 (+) (C (+) C), with P_{[1,m-1]} + P_prime + P_double_prime
    // + P_{[m+1,k]} as the combined projection list.
    
    aic_cstar_merge_sum(&B_new, &v_new, NULL, NULL, &v_part1, &v_2, A, prec)
    // (v_part1 = current v_comm with P_m's slot removed and P_[m+1..k] re-indexed)
    
    // Apply cor_improvement to restore c_0*eps error:
    aic_errreduce(&v_comm, &c0_measured, NULL, &v_new, eps, tol_abs, max_steps, prec)
    
    aic_cstar_subalg_clear(&sp)
    // Update projection list from v_comm.vE[] (diagonal elements)
```

**Implementation note on "v_part1":** The current v_comm has projections
P_1, ..., P_{m-1}, P_m, P_{m+1}, ..., P_k. When P_m is split, we build two sub-maps:
(a) the restriction of v_comm to the first m-1 + last k-m blocks (a new aic_dhom_v
over a B with those blocks, copying vE[mu(l,0,0)] for l != m-1 in 0-indexed),
(b) a fresh 2-block v_2 with vE[0] = P_prime, vE[1] = P_double_prime.
Then aic_cstar_merge_sum concatenates (a) and (b).

### 2.3 Termination guard and loop bound

The contradiction argument (.tex:1417-1426) proves the greedy loop terminates:
each split increases the number of projections by 1 (replacing 1 projection with 2).
Starting from 2 projections (or 1 for dim A = 1), we reach at most dim A projections
(since after dim A - 1 splits we have dim A 1-dimensional projections and cannot
split further). Therefore:

```c
#define STAGE1_MAX_ITERS (A->dim_A + 1)
slong iter = 0;
while (any_dim_S_gt_1) {
    ASSERT(++iter <= STAGE1_MAX_ITERS,
           "Stage 1 loop exceeded dim_A + 1 iterations — contradiction: "
           "th_main Stage 1 termination violated (bead aic-3qv stop condition)");
    ...
}
```

This is a fail-loud guard (Rule 4). Exceeding `STAGE1_MAX_ITERS` iterations is
a stop condition: it cannot happen if cor_improvement and lem_nontriv_projection
behave as specified by the paper (.tex:461: dimension-independent constant).

### 2.4 dim S_{P_j} query

```c
slong d_j = aic_corner_dim_S(A, P_j, P_j, prec);
```

With the aic-qgs substrate fix, this now works on genuinely oblique P_j (e.g.
P_j = v_comm(Pi_j) after errreduce, which is in A but may be a non-Hermitian
oblique projection; FINDINGS §C10 RESOLVED). The `aic_corner_Co` path no longer
SIGABRTs on near-zero Gram off-diagonals.

**The projections P_j stored in the v_comm.vE[] are the diagonal elements**
`v_comm.vE[mu(j,0,0)]` for each 1x1 block j (single-block C with d[l]=1). The
Stage-1 loop iterates over these in-place.

---

## 3. Equivalence classes after Stage 1

After the Stage 1 loop, all projections P_1, ..., P_m satisfy dim S_{P_j} = 1
for all j. Now compute equivalence classes:

**Algorithm (union-find on {1, ..., m}):**

```c
slong parent[m];
for (slong j = 0; j < m; j++) parent[j] = j;   // initial: each j its own class

// For all pairs (j,k) with j < k:
for (slong j = 0; j < m; j++) {
    for (slong k = j + 1; k < m; k++) {
        int equiv = aic_corner_equiv_1d(A, P_j, P_k, prec);
        // = 1 iff dim S_{P_j,P_k} == 1  (.tex:1428)
        if (equiv) {
            // union j and k in parent[]
            union_find_merge(parent, j, k);
        }
    }
}
// Classes: C_0, C_1, ..., C_{num_classes-1}, each a subset of {0,...,m-1}
// (0-indexed for C; the .tex uses 1-indexed {1,...,m})
// Class representatives: sort each class so that C = {r, ...} with r the
// find(root) index; order projections within class by index for reproducibility.

// For each class C: P_C = sum_{j in C} P_j (n x n, Hermitian to O(eps))
```

**Off-diagonal certification (Stage-3 precondition):**

After partitioning, for all pairs of distinct classes C, D:

```c
ASSERT aic_cstar_off_diag_zero(A, P_C, P_D, prec)
    // = aic_corner_dim_S(A, P_C, P_D) == 0
    // By lem_add_dim (.tex:1443): dim S_{P_C,P_D} = sum_{j in C, k in D} dim S_{P_j,P_k} = 0
    // since j not ~ k for all j in C, k in D (distinct classes).
    // FAIL LOUD if this assert fires — the equivalence-class partition is inconsistent.
```

---

## 4. Stage 2 — inductive extension per equivalence class

For each class C = {r_0, r_1, ..., r_{s-1}} (sorted by index, 0-based), s = |C|:

### 4.1 Base case (r = 1, i.e. s = 1 or the first induction step)

For C with s = 1 (a class containing one projection P_{r_0}):
- v_1 : M_1 -> S_{P_{r_0}} is trivial: v_1(1) = Ptilde_{r_0} = Co_{P_{r_0}}(P_{r_0})
  (the unit of S_{P_{r_0}}).
- Build: aic_dhom_B B_1 with num_blocks=1, d[0]=1; v_1.vE[0] = Ptilde_{r_0}.
- Apply aic_errreduce to normalize (at eps' not eps; see §4.4 below).

For C with s >= 2, the base case is r=1:
- A_1 = S_{P_{r_0}}, wrapper `sp1`.
- v_1 : M_1 -> S_{P_{r_0}}, v_1(1) = Ptilde_{r_0} = sp1.Ptilde.
- B_1: num_blocks=1, d[0]=1, dim_B=1. vE[0] = sp1.Ptilde.

### 4.2 Inductive step r-1 -> r (r = 2, ..., s for class C = {0,...,s-1})

Let P_r denote P_{r_{r-1}} (0-indexed; the r-th element of the sorted class).
Let P_{[1,r-1]} = sum of the first r-1 class projections (P_{r_0} + ... + P_{r_{r-2}}).
Let P_{[1,r]} = P_{[1,r-1]} + P_r.

**Step 0: assert S_{P_r} is 1-dim and class-consistent.**

```c
ASSERT aic_corner_dim_S(A, P_r, P_r, prec) == 1      // P_r is 1-dimensional
ASSERT aic_corner_equiv_1d(A, P_{r_0}, P_r, prec) == 1  // j ~ r for base j=r_0
// The design §G8 assert: dim S_{P_{[1,r-1]}, P_r} = r-1 (by lem_add_dim).
// Check:
slong d_off = aic_corner_dim_S(A, P_{[1,r-1]}, P_r, prec);
ASSERT(d_off == r - 1,
       "Stage 2 induction: dim S_{P_{[1,r-1]},P_r} != r-1 at step r=%ld", r);
// This is the S_{P,Q}!=0 hypothesis of lem_extension (.tex:1379).
```

**Step 1: build A_r = S_{P_{[1,r]}} wrapper.**

```c
aic_cstar_subalgebra sp_r;
aic_cstar_subalg_build(&sp_r, A, P_sum_r, prec);
// sp_r.A is the S_{P_{[1,r]}} eps'-C* algebra wrapper.
// sp_r.ctx->Co_P = Co_{P_{[1,r]}, P_{[1,r]}}  (d x d coordinate idempotent).
```

**Step 2: compress v_{r-1} into A_r (.tex:1437).**

```c
// v_{r-1}' : M_{r-1} -> A_{r-1}' where A_{r-1}' = S_{P_{[1,r-1]}} inside S_{P_{[1,r]}}
aic_dhom_v vr_prime;
aic_dhom_v_init(&vr_prime, v_{r-1}.B, &sp_r.A);
for (slong i = 0; i < v_{r-1}.B->dim_B; i++) {
    aic_corner_apply(vr_prime.vE[i], &sp_r.A, sp_r.ctx->Co_P, v_{r-1}.vE[i], prec);
    // = Co_{P_{[1,r]}}(v_{r-1}(E_i)) (.tex:1437)
}
// Compress the projections:
aic_corner_apply(P_sum_r_minus1_prime, &sp_r.A, sp_r.ctx->Co_P, P_sum_r_minus1, prec);
aic_corner_apply(Pr_prime,             &sp_r.A, sp_r.ctx->Co_P, P_r,            prec);
// P_{[1,r-1]}' = Co_{P_{[1,r]}}(P_{[1,r-1]}) (.tex:1438)
// P_r'         = Co_{P_{[1,r]}}(P_r)         (.tex:1439)
```

**Step 3: call lem_extension on the A_r wrapper.**

```c
// The A_r wrapper is now unblocked: aic-qgs resolved §C10 (FINDINGS §C10 RESOLVED).
// The S_P-wrapper-as-parent path for lem_extension runs end-to-end for genuinely
// oblique eta>0 (note: the compressed P_r' is an oblique element; the dim_S(P,Q)=n
// CONCLUSION guard, not the per-j guard, is used — FINDINGS §C10).
aic_dhom_B Br;
aic_dhom_v vr_plus;
arb_t mult_def_r, sigma_min_r;
// ...init mult_def_r, sigma_min_r...
aic_cstar_lem_extension(&Br, &vr_plus,
                        mult_def_r, sigma_min_r,
                        &vr_prime, P_sum_r_minus1_prime, Pr_prime,
                        &sp_r.A, prec);
// v_{r-1}^+ : M_r -> S_{P_{[1,r]}} (O(eps')-isomorphism)
```

**Step 4: apply cor_improvement to get v_r.**

```c
aic_dhom_v vr;
arb_t c0_r;
aic_errreduce(&vr, &c0_r, NULL, &vr_plus, eps_prime, tol_abs, max_steps, prec);
// v_r : M_r -> S_{P_{[1,r]}} (c_0*eps'-isomorphism, .tex:1441)
// c0_r is the measured c_0 for this step (GATE against nominal c_0; see §6).
```

**Step 5: store v_r as v_{r-1} for the next iteration; clear intermediates.**

```c
aic_cstar_subalg_clear(&sp_r);
aic_dhom_v_clear(&vr_prime);
aic_dhom_B_clear(&Br_plus);
aic_dhom_v_clear(&vr_plus);
// Move vr -> v_{r-1} for next step (or accumulate v_C if r == s).
```

**At r = s:** v_C = v_s is a c_0*eps'-isomorphism M_{|C|} -> S_{P_C}.
The associated B_C has num_blocks=1, d[0]=|C|.

### 4.3 Classes of size 1

For a class C = {j} with |C| = 1: v_C : M_1 -> S_{P_j}. The construction is
trivial (no lem_extension needed): v_C(1) = Ptilde_j. After aic_errreduce,
v_C is a c_0*eps'-isomorphism (measured defect ~ machine-zero since A_1 = S_{P_j}
is a genuine 1-dim C* algebra = C * Ptilde_j, and v_C is the exact inclusion).
B_C: num_blocks=1, d[0]=1.

### 4.4 Error bookkeeping: eps' management

The parameter eps' is defined at `.tex:1428`:
"Let eps' = O(eps) denote a suitable upper bound such that eps' >= c_0*eps, and
moreover, S_P for any c_0*eps-projection P is an eps'-C* algebra."

Constructively, eps' is NOT computed analytically. It is set as:

```c
double eps_prime = AIC_STAGE2_EPS_PRIME_FACTOR * eps_A;
// AIC_STAGE2_EPS_PRIME_FACTOR = 4.0 (a dimension-independent constant empirically
// satisfying the "S_P is eps'-C* for any c_0*eps projection P" requirement;
// the actual defect of the S_P wrapper's star is O(eta) with small constant;
// 4.0 gives a 50% margin over the empirical ~2.7 * eps upper bound from errreduce).
```

where `eps_A` is the associativity defect of the parent A (the `eps` argument to
`aic_cstar_build`). The aic_errreduce calls in Stage 2 pass `eps_prime` (not `eps_A`)
as the floor, since they operate on sub-algebras A_r which are eps'-C* algebras.

The nominal c_0 is fixed by the FIRST aic_errreduce call (see §6).

---

## 5. Stage 3 — merge equivalence classes

After Stage 2: num_classes = #{C_0, C_1, ..., C_{q-1}} (q classes).
Have v_{C_i} : M_{|C_i|} -> S_{P_{C_i}} for each i.

**Iterative merge:**

```c
aic_dhom_B B_merged;
aic_dhom_v v_merged;
// Initialize with the first class:
B_merged = B_{C_0};  v_merged = v_{C_0};

FOR i = 1, ..., q-1:
    // Certify off-diagonal zero (.tex:1443, lem_add_dim corollary):
    ASSERT aic_cstar_off_diag_zero(A, P_{C_0 through C_{i-1} sum}, P_{C_i}, prec)
    // (Checked against each new class before merging; the sum is the running P_total.)
    
    aic_dhom_B B_new;
    aic_dhom_v v_new;
    aic_cstar_merge_sum(&B_new, &v_new, NULL, NULL,
                        &v_merged, &v_{C_i}, A, prec);
    // merged v: (X, Y) |-> v_merged(X) + v_{C_i}(Y) : B_merged (+) M_{|C_i|} -> A
    
    // Apply cor_improvement to reduce to c_0*eps':
    aic_dhom_v v_next;
    arb_t c0_i;
    aic_errreduce(&v_next, &c0_i, NULL, &v_new, eps_prime, tol_abs, max_steps, prec);
    // GATE c0_i against nominal c_0 (see §6)
    
    aic_dhom_B_clear(&B_merged);  aic_dhom_v_clear(&v_merged);
    aic_dhom_B_clear(&B_new);     aic_dhom_v_clear(&v_new);
    B_merged = B_new_after_errreduce;  v_merged = v_next;
```

**After q-1 merges:** v_final = v : B -> A, a c_0*eps'-isomorphism where
B = (+)_{C} M_{|C|} (num_blocks = q, d[i] = |C_i|).

**ASSERT bijectivity (cor_improvement last sentence, `.tex:1318`):**

```c
arb_t a_bij;
int bij = aic_errreduce_is_bijective(&a_bij, &v_final, prec);
ASSERT(bij, "Stage 3: v_final is not bijective (sigma_min ~ 0 or dim_B != dim_A); "
            "stop condition — consult FINDINGS §D1/§D2");
// bij iff sigma_min(coord matrix of v) > 0 AND dim_B == dim_A.
```

---

## 6. c_0 management (FINDINGS §G5/§D2)

**The c_0 contract:**

`.tex:1415`: "Let c_0 be the constant from Corollary cor_improvement."
The paper treats c_0 as a fixed, universal constant. The implementation uses the
MEASURED c_0 from the first aic_errreduce call, per FINDINGS §D2.

```c
arb_t c0_nominal;   // fixed after first aic_errreduce call
double c0_nominal_d = 0.0;

// First aic_errreduce call (Stage 1, after the initial 2-projection split):
aic_errreduce(&v_comm, &c0_nominal, NULL, &v_raw, eps, tol_abs, max_steps, prec);
c0_nominal_d = arf_get_d(arb_midref(c0_nominal), ARF_RND_NEAR);

// Every subsequent aic_errreduce call:
arb_t c0_this;
aic_errreduce(&v_out, &c0_this, NULL, &v_in, eps_prime, tol_abs, max_steps, prec);
double c0_this_d = arf_get_d(arb_midref(c0_this), ARF_RND_NEAR);
// GATE: fail loud if this step's c_0 exceeds the nominal by more than 3x:
ASSERT(c0_this_d <= 3.0 * c0_nominal_d + 5.0,   // 5.0 for the tol_abs regime
       "c_0 drift: step c_0 = %.3f, nominal c_0 = %.3f (>3x); "
       "stop condition — FINDINGS §D2 universality canary (bead aic-1bc)", ...);
```

The `3.0` factor gives tolerance for variation across sub-algebras (which may
have slightly different eps' values). The `5.0` addend handles the eta~0 regime
where c_0_this_d is near 0 and the ratio is noisy.

**The universality canary** (the headline test, see §8.2): the measured c_0 must
not grow with dim A. The cstar_build loop inherits the errreduce T4 canary; the
I5 test extends it to sweep the FULL master-loop output.

**FINDINGS §D2 note:** The analytic c_0 is numerically ~2.0-2.7 (FINDINGS §D2
measured corpus). AIC_ERRREDUCE_C0_CERT = 10.0 is the generous fail-loud ceiling.
The master-loop c_0 gate uses 3*nominal to catch drift while not over-constraining
the slightly different eps' regime.

---

## 7. Public API

### 7.1 Entry-point signature

Add to `include/aic_cstar.h` (after the Increment 4 block):

```c
/* ============================ Increment 5 =============================== *
 * aic_cstar_build (.tex:1414-1444): the MASTER LOOP — constructive proof of
 * th_main (.tex:460): every finite-dim eps-C* algebra A is O(eps)-isomorphic to
 * a genuine C* algebra B = (+)_C M_{|C|}. The isomorphism v : B -> A is a
 * c_0*eps'-isomorphism (c_0 from cor_improvement, .tex:1317; eps' = O(eps)).
 * The constant c_0 is DIMENSION-INDEPENDENT (.tex:461; see FINDINGS §D2).
 *
 * THREE-STAGE LOOP (.tex:1414-1443):
 *   Stage 1: greedy commutative skeleton (aic_projection_nontrivial + cor_merge_sum
 *            + cor_improvement, iterated <= dim_A times; terminates by the
 *            maximum-dimensionality contradiction .tex:1417-1426).
 *   Stage 2: inductive extension per equivalence class (aic_cstar_lem_extension
 *            applied inductively M_{r-1} -> M_r, .tex:1435-1441; aic-qgs unblocked
 *            the oblique-wrapper-as-parent path, FINDINGS §C10 RESOLVED).
 *   Stage 3: merge classes via cor_merge_sum + cor_improvement (.tex:1443).
 *
 * INPUTS:
 *   A   : BORROWED eps-C* algebra (the aic_ecstar from aic_assoc_ecstar_from_phi
 *         or a genuine M_d via aic_cstar_matrix_algebra for testing). ASSERT A->n > 0.
 *   eps : the eps of A (e.g. aic_ecstar_defect_assoc for a genuine eps-C* algebra).
 *         This is the O(eps) denominator for the isomorphism certificate and c_0.
 *   prec: arb working precision.
 *
 * OUTPUTS:
 *   B_out    : OUTPUT OWNED genuine C* algebra B = (+)_C M_{|C|} (aic_dhom_B).
 *              Allocated here. num_blocks = number of equivalence classes;
 *              d[l] = |C_l|; dim_B = sum_l d[l]^2 = dim_A (bijective case).
 *              Free with aic_dhom_B_clear(B_out).
 *   v_out    : OUTPUT OWNED c_0*eps'-isomorphism v : B -> A (aic_dhom_v).
 *              Allocated here. v_out->B points at B_out (keep B_out alive for
 *              v_out's lifetime); v_out->A = A. Free with aic_dhom_v_clear(v_out).
 *   iso_def  : OUTPUT certified arb ball: the isomorphism defect of v, defined as
 *              max(mult_def, lower_gap, norm_excess, unit_def) of the final v_out
 *              (= the max of the four aic_errreduce_defects). A true c_0*eps'-
 *              isomorphism has iso_def <= c_0*eps' = AIC_ERRREDUCE_C0_CERT * max(eps', tol).
 *              Caller-initialised arb_t; pass NULL to skip.
 *
 * ASSERTS (fail loud, Rule 4):
 *   A->n >= 1 and A->dim_A >= 1;
 *   Stage 1 loop terminates in <= dim_A + 1 iterations (contradiction, .tex:1417);
 *   Each aic_projection_nontrivial call succeeds (fail loud if spectrum degenerate:
 *     the aic-3qv escalation, FINDINGS §D1);
 *   Each aic_errreduce call certifies a c_0*eps'-inclusion (AIC_ERRREDUCE_C0_CERT gate);
 *   c_0 across all steps stays within 3x the nominal (universality canary gate);
 *   Stage 3 off-diagonal zero: aic_cstar_off_diag_zero for all distinct class pairs;
 *   Final v is bijective: aic_errreduce_is_bijective (dim_B == dim_A AND sigma_min > 0).
 *
 * On the eta=0 oracle (A from an exact idempotent Phi): B matches the block structure
 * of th_idemp_structure (aic_idemp_decompose) exactly (num_blocks, d[l] in bijection),
 * iso_def ~ machine-zero, v bijective. (The canonical cross-check, Rule 6.)
 * Free: aic_dhom_B_clear(B_out) then aic_dhom_v_clear(v_out).
 */
void aic_cstar_build(aic_dhom_B *B_out, aic_dhom_v *v_out, arb_t iso_def,
                     const aic_ecstar *A, double eps, slong prec);
```

### 7.2 Ownership and lifetime contract

- `B_out` is OWNED by the caller. Free with `aic_dhom_B_clear(B_out)`.
- `v_out` is OWNED by the caller. `v_out->B` points at `B_out` (BORROWED);
  keep `B_out` alive for `v_out`'s lifetime.
- `v_out->A` points at `A` (BORROWED); keep `A` alive for `v_out`'s lifetime.
- `iso_def` is an arb_t that the caller init'd and owns; pass NULL to skip.

---

## 8. Test design (Rules 6/7 — write the spec before the code)

### 8.1 T1 — eta=0 oracle (the cleanest ground truth)

**Fixture:** An exact idempotent channel whose image algebra has multiple blocks.

Concrete choices (from `test_idemp.h`):
- `make_block_cond_exp(phi, d=4, m=2)`: Img Phi = B(C^2) (+) B(C^2),
  dim_A = 4 + 4 = 8, two blocks of size 2.
- `make_block_cond_exp(phi, d=6, m=2)`: dim_A = 4 + 16 = 20, blocks of size 2 and 4.
- `make_noiseless_subsystem(phi, dL=3, dE=2)`: Img Phi = B(C^3), dim_A = 9,
  single block of size 3.
- `make_block_cond_exp(phi, d=6, m=3)`: dim_A = 9 + 9 = 18, two blocks of size 3.

**Procedure:**

```c
// (a) Run aic_cstar_build:
aic_assoc_ecstar_from_phi(&ae, phi, prec);
aic_cstar_build(&B_main, &v_main, &iso_def, &ae.A, 0.0, prec);

// (b) Run aic_idemp_decompose independently:
aic_idemp_decompose(&decomp, phi, prec);
// decomp gives: dim_A, dim_M, and the structure of ImgPhi.

// For make_block_cond_exp(d,m): th_idemp_structure gives
//   blocks: one of dim m^2 and one of dim (d-m)^2.
// For make_noiseless_subsystem(dL,dE): one block of dim dL^2.

// (c) ASSERTIONS:
// iso_def ~ machine-zero (< 1e-10 at prec=128):
ASSERT arb_lt_d(iso_def, 1e-10)

// B structure matches th_idemp_structure block sizes:
// Map B_main.num_blocks and B_main.d[] against the known block structure.
// For make_block_cond_exp(4,2): ASSERT B_main.num_blocks == 2,
//   {B_main.d[0], B_main.d[1]} is {2,2} (in some order).
// For make_noiseless_subsystem(3,2): ASSERT B_main.num_blocks == 1, B_main.d[0] == 3.

// v is bijective:
arb_t a_bij; arb_init(a_bij);
ASSERT aic_errreduce_is_bijective(&a_bij, &v_main, prec)
ASSERT arb_gt_d(a_bij, 0.5)   // sigma_min > 0.5
arb_clear(a_bij);
```

**Why this is the cleanest ground truth:** At eta=0, A is a GENUINE C* algebra.
The master loop must produce B isomorphic to A with zero defect. If B has the
wrong block structure (e.g. a single M_8 block instead of M_2 (+) M_2), the
isomorphism defect is O(1). The `aic_idemp_decompose` oracle is completely
INDEPENDENT of `aic_cstar_build` (different algorithm, different code path).

### 8.2 T2 — universality canary (the headline test, `.tex:484` failure mode)

**Purpose:** Verify that the measured c_0 = iso_def / eps does NOT grow with dim A.
A c_0 growing like dim A or sqrt(dim A) is the `.tex:484` failure mode (the
naive Haar/cohomology route) — a stop condition (FINDINGS §D2).

**Fixture:** A family of exact idempotent channels of increasing dimension, each
passed through `aic_assoc_ecstar_from_phi` at a small t to produce a genuine eta>0
algebra (so eps > 0 and the c_0 ratio is well-defined).

```c
// Sweep: make_mixconj(phi, d, m=1, t=0.06, prec) for d = 3, 4, 5, 6, 7, 8
// (dim A ~ d^2 to ~2d, varying from 9 to 64).
// Alt: make_block_cond_exp(phi, d, m=1) perturbed via t-mix for dim_A = 1 + (d-1)^2,
// varying from 4 to 25 for d=2..6.

double c0_arr[num_dims];
for (slong di = 0; di < num_dims; di++) {
    // build phi, build ae.A (eps_di = aic_ecstar_defect_assoc)
    aic_cstar_build(&B, &v, &iso_def, &ae.A, eps_di, prec);
    c0_arr[di] = arf_get_d(arb_midref(iso_def), ARF_RND_NEAR) / eps_di;
    // cleanup
}
double c0_max = max(c0_arr), c0_min = min(c0_arr);
// ASSERT: ratio c0_max / c0_min <= 1.6 (the errreduce T4 threshold, FINDINGS §D2)
ASSERT(c0_max / (c0_min + 1e-10) <= 1.6,
       "universality canary: c_0 ratio %.3f > 1.6; tex:484 failure mode — "
       "stop condition (FINDINGS §D2)", c0_max / c0_min);
```

**Mutation tooth:** Inject `c_0 *= sqrt(dim_B)` into one aic_errreduce call in
the loop. The canary ratio becomes sqrt(dim_max/dim_min) ~ 2.6 for the dim
range, exceeding 1.6 -> RED. (Per FINDINGS §D2: the errreduce T4 tooth already
mutation-proves this for the errreduce layer; I5 repeats it at the loop level.)

### 8.3 T3 — oblique eta>0 fixture

**Fixture:** `make_mixconj(phi, d=5, m=3, t, prec)` at t in {0.01, 0.03, 0.06}.

```c
// eta_proxy ~ 6.5 * t (from test_idemp.h make_mixconj comment)
// At t=0.06: eta ~ 0.039

aic_cstar_build(&B, &v, &iso_def, &ae.A, eps, prec);

// (a) iso_def is certified O(eta):
double c = arf_get_d(arb_midref(iso_def), ARF_RND_NEAR) / eta;
ASSERT(c < 20.0,    // generous bound; true c_0 ~ 2-5 from errreduce corpus
       "iso_def / eta = %.3f exceeds bound 20.0 (star path failed?)", c);

// (b) Star-vs-plain MAGNITUDE tooth (FINDINGS §C8):
// Mutate A.star_phi -> NULL (plain matmul) inside aic_ecstar_star.
// The iso_def with PLAIN product gives c_plain = iso_def_plain / eta.
// ASSERT c_plain > 0.2 (5x gap from c_star < 0.2).
// (The mutation is CONFIRMED RED by hand per §C8 discipline; not left in the suite.)
```

### 8.4 T4 — mutation teeth

These MUST be confirmed RED (by hand if the automated suite is kept clean):

**(a) star -> plain anywhere in the loop:** In any aic_dhom_defect_sweep call
inside the loop (via the vE star-multiplication path), mutate `aic_ecstar_star`
to `acb_mat_mul`. At eta>0, the iso_def/eta ratio exceeds 0.2 -> RED (FINDINGS §C8).

**(b) Wrong B block count:** Force B to have one fewer block (merge two classes
that should be distinct). The iso_def jumps to O(1) because the isomorphism
collapses a direction (sigma_min -> 0) -> bijectivity assert fires -> RED.

**(c) Skip errreduce in Stage 1 / 2 / 3:** Remove any one aic_errreduce call.
The c_0 for SUBSEQUENT errreduce calls grows beyond 3*nominal -> the c_0 drift
gate fires -> RED. (Or: the intermediate defect is too large for the next
lem_extension call, which ASSERTs a delta_max bound on its input.)

**(d) Frobenius Pi_A projection instead of Co_P(Phi_tilde):** The S_P wrapper
uses `Phi_tilde(X) = X * 1_n` (the oblique projection). If the compression were
replaced with the Frobenius projector (plain Pi_A), the star defect of the
projected P would be O(1) instead of O(eta) (FINDINGS §C3). The aic_projection_nontrivial
step 4 would produce a P with ||P*P-P||_op ~ O(1) -> the projection's own
ASSERT(delta < 0.1) fires -> RED.

**(e) Wrong equivalence class partition:** Force j !~ k to be classed together
(set parent[j] = parent[k] even when dim S_{P_j,P_k} = 0). Stage-3's
off-diagonal-zero assert fires (aic_cstar_off_diag_zero returns 0) -> RED.

---

## 9. New code inventory + file split

### 9.1 File plan

The master loop is ~150 LOC for the loop logic plus ~80 LOC for helper functions
(union-find, eps' bookkeeping, projection list management). Total I5 new code:
~230 LOC across new source + header additions. To respect the 200 LOC / file
limit (Rule 10):

```
src/aic_cstar_build.c   (~170 LOC)
    - aic_cstar_build()   (the public API entry point)
    - stage1_greedy()     (Stage 1 loop, ~60 LOC)
    - stage2_extend()     (Stage 2 per-class induction, ~50 LOC)
    - stage3_merge()      (Stage 3 class merge, ~30 LOC)
    - static helpers:
        - find_dim_gt1()       (scan P_j list for dim S_{P_j} > 1, ~8 LOC)
        - compute_classes()    (union-find partition, ~25 LOC)
    Total: ~173 LOC (within limit)

include/aic_cstar.h     (add the Increment 5 block: API declaration + docstring,
                         ~80 LOC added; total file ~500 LOC, Rule 10 exception:
                         header-only additions are exempt from the 200 LOC limit
                         per the "split source files" convention)

tests/test_cstar_build.c   (~180 LOC: T1 eta=0 oracle, T2 universality canary,
                             T3 oblique, T4 mutation teeth comments)
```

**No new source file split beyond aic_cstar_build.c is needed** (the three
stages together are ~150 LOC, fitting with helper functions in one file).

### 9.2 Reused (all CLOSED)

Every primitive called by I5 is from a CLOSED, green module:

| Primitive | Module | Used for |
|---|---|---|
| `aic_projection_nontrivial` | aic-mqf (CLOSED) | Stage 1 split |
| `aic_cstar_subalg_build/clear` | aic-097 I1 (CLOSED) | S_P wrapper |
| `aic_corner_dim_S` | aic-czm (CLOSED) | dim S_{P_j} query, off-diag |
| `aic_corner_equiv_1d` | aic-czm (CLOSED) | equivalence class computation |
| `aic_corner_Ptilde` | aic-czm (CLOSED) | Ptilde_j for base case |
| `aic_corner_apply` | aic-czm (CLOSED) | Stage 2 compression (.tex:1437) |
| `aic_cstar_merge_sum` | aic-097 I2 (CLOSED) | Stage 1 and Stage 3 merge |
| `aic_cstar_off_diag_zero` | aic-097 I2 (CLOSED) | Stage 3 off-diag cert |
| `aic_cstar_lem_extension` | aic-097 I4 (CLOSED) | Stage 2 induction |
| `aic_errreduce` | aic-t81 (CLOSED) | Error reduction at all stages |
| `aic_errreduce_is_bijective` | aic-t81 (CLOSED) | Final bijectivity check |
| `aic_dhom_v_init/clear` | aic-c1n (CLOSED) | v management |
| `aic_dhom_B_init/clear` | aic-c1n (CLOSED) | B management |
| `aic_dhom_B_index` | aic-c1n (CLOSED) | Matrix-unit indexing |

### 9.3 Genuinely new code in aic_cstar_build.c

1. **`compute_classes()` (~25 LOC):** A simple union-find on the m projection indices.
   Standard: path-compression + union by rank. Allocates `parent[m]`, `rank[m]`,
   computes class partition via `aic_corner_equiv_1d` for all pairs. Returns a
   `class_id[m]` array and `num_classes`. No new FLINT/LAPACK primitive.

2. **`stage2_extend()` (~50 LOC):** The r=1..s inductive loop per class. Allocates
   the running `P_sum[r]` operator (n x n) updated by adding `P_{r_l}` at each step,
   the Co_{P_{[1,r]}} compression (via the wrapper's ctx->Co_P), and v_prime (the
   compressed v_{r-1}). Calls `aic_cstar_lem_extension` and `aic_errreduce`. No new
   primitive; pure loop + lifetime management.

3. **Projection list management (~15 LOC):** Maintain a `flint_malloc`'d array of
   `acb_mat_t` (one per Stage-1 projection) that is updated when a projection P_m is
   split into P' and P''. Resize via `flint_realloc` (safe since each element is an
   initialized `acb_mat_t` moved by value). The "split" operation is: remove index m,
   insert P' at m and P'' at m+1 (shift the tail).

4. **`aic_cstar_build()` (~30 LOC):** Driver: calls stage1, compute_classes, stage2
   per class, stage3, bijectivity assert, fills iso_def.

---

## 10. Open escalations carried forward (do NOT fake)

### §D1: the Omega(1) spectral gap is per-instance-certified, not universally proven

`aic_projection_nontrivial` fail-loud-aborts (the aic-3qv stop condition) if a
Stage-1 sub-algebra has a degenerate spectrum (all eigenvalues coincident) where
dim > 1. The paper's Lefschetz-Hopf argument proves such a gap must exist for
any eps-C* algebra of dim > 1 (.tex:931: "any eps-C* algebra A such that
1 < dim A < infinity has a nontrivial O(eps)-projection"). The constructive route
CERTIFIES the gap per-instance and aborts loudly if it fails.

**If the abort fires on a genuinely >1-dimensional sub-algebra, it is a stop
condition:** attach the (A, P_m, sp.A) data and escalate with bead aic-3qv.
It cannot be silently skipped. This is a fundamental theoretical gap that
could, in principle, arise if the paper's Lefschetz argument fails
constructively for some finite-dim instance.

### §G5/§D2: analytic c_0 is unmeasured

The master loop uses the MEASURED c_0 from the first aic_errreduce call. The
analytic c_0 (bead aic-1bc) would make the universality claim rigorous. The
measured c_0 ~2-3 (FINDINGS §D2) is the empirical witness; the I5 universality
canary test (T2) is the machine-testable surrogate.

---

## 11. Summary: what the implementer does, step by step

1. Add the Increment 5 block to `include/aic_cstar.h` (signature + docstring
   above), keeping the I1-I4 blocks intact.

2. Create `src/aic_cstar_build.c` (~170 LOC):
   - Static `compute_classes(A, P_arr, m, prec)` -> `class_id[]`, `num_classes`.
   - Static `stage1_greedy(A, eps, prec)` -> returns the final `{P_j}` list and
     the corresponding `v_comm` (aic_dhom_v). Guards the loop count, calls
     `aic_projection_nontrivial` + S_P nontriviality check + `aic_cstar_merge_sum`
     + `aic_errreduce`.
   - Static `stage2_extend(A, class_C, len_C, v_comm_proj_subset, eps_prime, prec)`
     -> returns `v_C` (aic_dhom_v) and `B_C` (aic_dhom_B) for class C.
   - Static `stage3_merge(A, v_C_arr, B_C_arr, num_classes, eps_prime, prec)` ->
     returns final `v` and `B`.
   - Public `aic_cstar_build(B_out, v_out, iso_def, A, eps, prec)`.

3. Create `tests/test_cstar_build.c` (~180 LOC): T1 (eta=0 oracle with
   `make_block_cond_exp(4,2)` and `make_noiseless_subsystem(3,2)`), T2
   (universality canary, mixconj sweep d=3..7), T3 (oblique eta>0 with
   `make_mixconj(5,3,0.06)`), plus comments for the T4 mutation teeth.

4. Run `make test` (all existing tests stay green). Then run
   `build/test_cstar_build` (the new tests).

---

*End of spec.*
