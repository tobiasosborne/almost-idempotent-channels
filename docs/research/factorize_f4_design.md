# `factorize` F4 — design / feasibility spec: end-to-end verification, dual channels, dimension-independence canary

> **Status:** DESIGN / FEASIBILITY (read-only on `.c`/`.h`; this is the only new
> file). Bead `aic-tff` (module `factorize`), increment **F4**, the FINAL one.
> Realizes the headline `th_factorization` end-to-end bounds (`approximate_algebras.tex:2730-2739`),
> the u-isomorphism corollary (`.tex:2742-2747`), and the dual channels
> `Enc=Υ*`, `Dec=Δ*` (`.tex:2159`, `Enc` `.tex:280`, `Dec` `.tex:286`).
> Builds on F1/F2/F3 (committed, green, 70 checks in `test_factorize`).
> Prerequisites read: `factorize_d4_research.md` §5/§6, FINDINGS §C12/§D2/§D4/§C13(c)/§C12.O2-PIN,
> `aic_cbnorm.h`, `aic_ucp.h`, `aic_factorize.h`, `aic_dhom.h`, `aic_opspace.h`,
> `tools/gen_fixtures_opspace_o2.jl`, `tests/test_opspace_o2.c`.
>
> **Verdict, up front.** The orchestrator's proposed **F4.1 (no Julia) / F4.2
> (Julia+MOSEK)** split is **CORRECT and adopted** with one refinement below.
> F4.1 **closes the constructive headline** (produces certified `Δ, Υ, B`, the
> η=0-exact oracle, a per-instance certified `O(η)` UPPER bound, and the four
> dual channels). F4.2 is the **universality certification** — the only thing
> that can faithfully test dimension-independence, because the eig-free bracket's
> looseness is `O(N)` (dim-dependent) and would FAKE/MASK a `c=O(N)` law (the
> §C12 trap). **No NEW convention needs an empirical pin** — every convention F4
> touches is already pinned in `aic_cbnorm`/O2 (§E). The `ΥΔ−1_B`-on-`⊕M_{d_l}`
> cb-norm needs a definition CHOICE (not a competition): I recommend the
> **square-`M_{n_B}` Convention-A Choi with `Δ`/`Υ` viewed as maps on the ambient
> `M_{n_B}` block-diagonal representative**, reconciled with `UpsDel2` by verifying
> the cleaner `‖ΥΔ−1_B‖_cb` equivalent (§A.2, §C.canary).

---

## 0. What F1–F3 give F4 (the inputs, all built + green)

From `aic_factorize.h` (the `aic_factorize` handle `F`) and the closed modules:

| Object | How to obtain | Shape / type |
|---|---|---|
| `Δ̃ = ι∘v` (non-UCP) | `aic_factorize_delta_tilde(out, F, X, prec)` | `B → B(H)`, X: `n_B×n_B`, out: `N×N` |
| `Υ̃ = v⁻¹∘Φ̃` (non-UCP) | `aic_factorize_upsilon_tilde(out, F, W, prec)` | `B(H) → B`, W: `N×N`, out: `n_B×n_B` |
| `Δ` (UCP encode) | `aic_factorize_delta(out, F, X, prec)` after `aic_factorize_delta_build(F, prec)` | `B → B(H)` |
| `Υ` (UCP decode) | `aic_factorize_upsilon(out, F, W, prec)` after `aic_factorize_upsilon_build(F, tol, prec)` | `B(H) → B` |
| `B = ⊕_l M_{d_l}` | `F->v->B` (`aic_dhom_B`: `num_blocks=m`, `d[l]`, `n_B`, `dim_B`) | block-diag rep size `n_B`, dim `dim_B = Σ d_l²` |
| `Φ` (original UCP) | `F->phi` (`aic_ucp_kraus`, `dim_K==dim_H==N`, `r` Kraus) | `B(H) → B(H)` |
| matrix units of B | `aic_dhom_B_matunit(out, B, i)` (linear `i ∈ [0,dim_B)`), `aic_dhom_B_unit(out, B)` = `I_B` | `n_B×n_B` block-diag |
| `I_B = ⊕_l I_{d_l}` | `aic_dhom_B_unit` | `n_B×n_B` |

The key precision-relevant facts (CLAUDE.md callouts; FINDINGS §C2/§C3): `Φ̃(W)`
is the **oblique** image of `W` on `A` (`aic_assoc_superop_apply` over
`Aec->S_tilde`, NOT the Frobenius projection); the product on `A` is the
Choi–Effros star, never plain matmul; and `Δ`, `Υ` are UCP (`Δ(I_B)=1_H`,
`Υ(1_H)=I_B`, both CP by construction). All four maps are exact-by-construction
at η=0 and `O(η)`-defective at η>0.

`F->phi`, `F->v`, `F->Aec` are borrowed and must outlive `F`. The convention for
`Φ`: Kraus stored Heisenberg-style `Φ(X)=Σ_a K_a† X K_a`, Choi is Convention A
(`aic_ucp.h:23-43`), `aic_ucp_kraus_to_choi`/`aic_ucp_compose`/`aic_ucp_choi_diff`
are available.

---

## A. The composed maps `ΔΥ`, `ΥΔ`, and their Choi matrices

### A.1 `ΔΥ : B(H) → B(H)` and `J_DelUps = Choi(ΔΥ) − Choi(Φ)`

`ΔΥ` is the Heisenberg-picture self-map on `B(H) = M_N`:

```
ΔΥ(W) = Δ(Υ(W))        (W ∈ M_N)
```

In code: `aic_factorize_upsilon(U, F, W, prec)` (→ `n_B×n_B`), then
`aic_factorize_delta(out, F, U, prec)` (→ `N×N`). Both already require their
respective `_build` calls; F4's verify routine asserts `F->delta_ready` and
`F->upsilon_ready` up front (Rule 4).

**Convention-A Choi of `ΔΥ`** (`aic_ucp.h:23-29`, the paper/Watrous/QETLAB
convention F2 already uses in `aic_factorize_delta_block_choi`):

```
Choi(ΔΥ)[p*N + a, q*N + b] = ΔΥ(E_pq)[a,b]       (p,q ∈ [0,N) INPUT/major,
                                                  a,b ∈ [0,N) OUTPUT/minor)
```

i.e. `Choi(ΔΥ) = Σ_{p,q} E_pq ⊗ ΔΥ(E_pq)`, an `N²×N²` matrix. Build it by looping
`E_pq` over the `N²` ambient matrix units of `M_N`, applying `ΔΥ`, and writing the
`N×N` output block into position `(p,q)` (the F2 block-Choi pattern, generalized
to the full `M_N` input basis). This is a HP map's Choi, Hermitian but indefinite.

**`Choi(Φ)`** is `aic_ucp_kraus_to_choi(Cphi, F->phi, prec)` — `N²×N²`, Convention A.

```
J_DelUps = Choi(ΔΥ) − Choi(Φ)            (N²×N², Hermitian, NOT PSD)
```

Because `ΔΥ` is a difference-of-no, but `ΔΥ` itself is UCP (`Δ`,`Υ` UCP ⇒ `ΔΥ`
UCP) — so `Choi(ΔΥ)` IS PSD. `J_DelUps` is the difference of two PSD Convention-A
Choi matrices (both unital), hence Hermitian and indefinite, exactly the shape
`aic_ucp_choi_diff` produces for `Φ²−Φ`. `‖ΔΥ−Φ‖_cb` is the diamond/cb-norm of
this HP map (DelUps, `.tex:2733`).

> **Note (reuse):** `ΔΥ` is *not* an `aic_ucp_kraus` object (it is `Δ∘Υ` applied
> functionally), so we cannot route through `aic_ucp_compose` to get its Choi. We
> build `Choi(ΔΥ)` directly by the `E_pq`-loop. (One *could* extract Kraus of `Δ`
> and `Υ` and compose, but the functional `E_pq`-loop is simpler, exact, and is
> what F2's block-Choi already does — reuse the pattern, not a new primitive.)

### A.2 `ΥΔ : B → B`, `ΥΔ − 1_B`, and THE SUBTLETY (ranked risk #1)

`ΥΔ(X) = Υ(Δ(X))` is a map on the direct-sum algebra `B = ⊕_l M_{d_l}`, NOT on a
full `M_N`. In code: `aic_factorize_delta(D, F, X, prec)` (X: `n_B×n_B`, D:
`N×N`), then `aic_factorize_upsilon(out, F, D, prec)` (out: `n_B×n_B`). `1_B` is
`aic_dhom_B_unit` and `ΥΔ(I_B) = I_B` exactly (both maps unital).

**The subtlety.** "What is `‖ΥΔ − 1_B‖_cb`?" The map's domain/codomain is
`B = ⊕_l M_{d_l}`, an operator space embedded block-diagonally in `M_{n_B}`. There
are three candidate cb-norm definitions:

(i) **Ambient-`M_{n_B}` Choi** (RECOMMENDED). View `ΥΔ` as a linear map
`M_{n_B} → M_{n_B}` by feeding it the block-diagonal representative. Build the
Convention-A Choi over the **full `M_{n_B}` matrix-unit basis** (`n_B²` units),
`J_UpsDel = Choi(ΥΔ − 1_B) = Σ_{p,q} E_pq^{(n_B)} ⊗ (ΥΔ−1_B)(E_pq)`,
an `n_B²×n_B²` Hermitian matrix. `‖ΥΔ − 1_B‖_cb` = the diamond norm of this map.

(ii) **Block-restricted Choi** over only the `dim_B = Σ d_l²` matrix units of `B`
(`aic_dhom_B_matunit`, the in-algebra basis), a `dim_B`-indexed object. This is
the cb-norm of `ΥΔ−1_B` as a map on the *operator space* `B` proper (Smith's
lemma still applies; the cb-norm is attained at the codomain ambient level).

(iii) **Per-block** cb-norms `max_l ‖(ΥΔ−1_B)|_{M_{d_l}}‖_cb`, which is WRONG —
it ignores the cross-block ampliation structure that the direct-sum cb-norm sees
(it would be a "test that cannot fail" in the §C12 sense if `ΥΔ` leaked between
blocks).

**Recommendation: route (i) — the ambient `M_{n_B}` Convention-A Choi.** Reasons:

- **It is what the existing toolbox certifies.** `aic_cbnorm_certify` (square
  case) and `aic_cbnorm_eigfree_ball_choi` both consume an `n²×n²` Convention-A
  Choi of a map on `M_n`; here `n = n_B`. Route (i) plugs directly in with `n =
  n_B`, no new primitive. The diamond-SDP convention (2/n, Tr_2=`partial_trace_right`)
  is the SQUARE self-map convention already pinned (§E) — and `ΥΔ−1_B` is a square
  HP self-map on `M_{n_B}`, exactly like `Φ²−Φ`, so the `aic_cbnorm_certify`
  (2/n) convention applies verbatim (NOT the rect `factor-1` convention, which is
  for the asymmetric `v*`/`(v⁻¹)*`).
- **It is a valid upper bound on the operator-space cb-norm (ii).** For a
  block-diagonal-valued map, the ambient cb-norm equals the operator-space cb-norm:
  feeding a non-block-diagonal `M_{n_B}` input to `ΥΔ`(via `Δ` ∘ `Υ`) — both `Δ`
  and `Υ` are defined on the full ambient (`Υ` takes any `N×N`; `Δ` takes the
  `n_B×n_B` rep and only reads the matrix-unit coordinates `aic_dhom_v_apply`
  uses). The off-block-diagonal entries of the `M_{n_B}` input are dropped by `Δ`
  (it reads only the block-diagonal coordinates), so the ambient Choi's
  off-block-diagonal input columns contribute zero and the ambient cb-norm equals
  the in-`B` cb-norm. This is the same reason the η=0 oracle gives `[0,0]` for both
  (i) and (ii). **Caveat to verify empirically (a probe, see §D):** that the
  ambient and block-restricted Choi give the SAME diamond norm at η>0 — if `Δ`
  does NOT zero the off-block-diagonal input (e.g. if a future `Δ` reads them),
  routes (i) and (ii) diverge and (ii) is the correct definition. The η=0 oracle
  is BLIND to this (`ΥΔ=1_B` exactly), so this needs an η>0 probe before F4.2
  fixes the convention. **This is the single empirical probe F4 needs** (it is a
  definition-correctness check, not a new convention pin).

**Reconcile with `UpsDel2` (`.tex:2736`).** The theorem states the bilinear bound
`‖Υ_n(Δ_n(X)Δ_n(Y)) − XY‖ ≤ O(η)‖X‖‖Y‖` for `X,Y ∈ M_n⊗B`, and the parenthetical
`.tex:2739` says it IMPLIES `‖ΥΔ−1_B‖_cb ≤ O(η)`. F4 has a CHOICE: verify `UpsDel2`
**directly** (a bilinear bound over `B`-matrix-unit pairs, ampliated) or verify the
equivalent **`‖ΥΔ−1_B‖_cb`** (`.tex:2739`).

**Recommend: verify `‖ΥΔ−1_B‖_cb` (the cleaner route).** Rationale:

- `UpsDel2` is a *bilinear* (two-argument) bound requiring an ampliated double
  loop over `M_n⊗B` matrix-unit pairs — `O(n² dim_B²)` map applications, and the
  ampliation `n` is itself a Smith-truncation question. `‖ΥΔ−1_B‖_cb` is a single
  linear-map cb-norm with a known finite truncation (Smith: attained at the
  codomain ambient level `n_B`, FINDINGS §D3) and a single committed SDP point.
- The paper's own framing (`.tex:2739`) treats `‖ΥΔ−1_B‖_cb` as the *headline*
  consequence; `UpsDel2` is the intermediate. The theorem's contract is the
  cb-norm bound.
- `‖ΥΔ−1_B‖_cb` is what the dual `Dec∘Enc ≈ 1_{B*}` (§B.duals) certifies, so it
  is the load-bearing quantity for the factorization story.

F4 may ADD a `UpsDel2` spot-check (a few `X,Y` matrix-unit pairs at `n=1`, the
un-ampliated version `‖Υ(Δ(X)Δ(Y)) − XY‖ ≤ O(η)`) as a cheap independent
cross-check on the same `Δ,Υ` (it reuses F3's existing `aic_dhom_B_mul` for `XY`
and `aic_ucp_apply`-free plain-matmul for `Δ(X)Δ(Y)`), but the **certified
headline is `‖ΥΔ−1_B‖_cb`**. (Note: `Δ(X)Δ(Y)` here is the ORDINARY `M_N` product,
per the F2 docstring `.tex:2788` — not the star; `Υ` then maps back to `B`.)

---

## B. F4.1 (no Julia) — the constructive headline

F4.1 is the part that runs in `make test` with no Julia/MOSEK. It produces the
certified-per-instance headline and the four dual channels.

### B.1 η=0 ORACLE (the cleanest cross-check, Rule 6 / ladder rung 3)

At η=0, `Φ̃ = Φ`, `Δ̃ = ι∘v` reduces to the exact `th_idemp_structure` inclusion,
and `Υ̃ = v⁻¹∘Φ` to the exact `Γ∘Co_M` decode (F1 T1 already verifies the tilde
maps; F2 T3 and F3 T5 verify `Δ`,`Υ` reduce to the tilde maps). Therefore as MAPS:

```
ΔΥ = Φ          ⟹  J_DelUps = Choi(ΔΥ) − Choi(Φ) = 0   (to machine)
ΥΔ = 1_B        ⟹  J_UpsDel = Choi(ΥΔ − 1_B)      = 0   (to machine)
```

**Assertions (T7):** on the η=0 oracle fixtures (`make_block_cond_exp(4,2)`,
`make_block_cond_exp(6,3)` — the multi-block η=0 cases; and a `noiseless_subsystem`
if convenient):

- `‖J_DelUps‖_F < 1e-10` (Frobenius norm of the `N²×N²` difference).
- `‖J_UpsDel‖_F < 1e-10` (the `n_B²×n_B²` difference).
- `aic_cbnorm_eigfree_ball_choi(lo, hi, J_DelUps, N, prec)` returns `hi < 1e-9`
  (the bracket's `2‖J‖_F` collapses to ~0; the eig-free routine gives `[0,0]` at
  `J=0`). Same for `J_UpsDel` with `n = n_B`.
- (optional) the un-ampliated `UpsDel2` spot-check `‖Υ(Δ(X)Δ(Y)) − XY‖ < 1e-10`
  over a few B-matrix-unit pairs.

The η=0 oracle is the strongest single check: it pins the WHOLE pipeline (all
three CP-izations + the two unitalizations + the v-inversion) to exact reduction.

### B.2 PER-INSTANCE certified UPPER bound (T8, η>0, eig-free)

For an η>0 fixture, build `J_DelUps` and `J_UpsDel` as above and call the
always-valid eig-free bracket (no SDP, no Julia, never aborts):

```
aic_cbnorm_eigfree_ball_choi(lo_du, hi_du, J_DelUps, N,   prec);   // ‖ΔΥ−Φ‖_cb ≤ hi_du = 2‖J_DelUps‖_F
aic_cbnorm_eigfree_ball_choi(lo_ud, hi_ud, J_UpsDel, n_B, prec);   // ‖ΥΔ−1_B‖_cb ≤ hi_ud = 2‖J_UpsDel‖_F
```

These are RIGOROUS upper bounds (`aic_cbnorm.h:39`). At a FIXED instance with a
MEASURED `η` (`eta_proxy(phi, prec)`, the `test_factorize` helper), assert:

- `hi_du ≤ C_max · η` and `hi_ud ≤ C_max · η` for a generous fail-loud
  `C_max` (e.g. `40` — the eig-free `2N` looseness means `C_max` is inflated by
  `~2N`; the assertion is a coarse boundedness gate, NOT a tight constant). This
  certifies the `O(η)` SCALING per instance (the algorithm meets the theorem's
  bound shape at this instance).

**State clearly (the honest posture, FINDINGS §D4 / §C12 trap):** the eig-free
bracket `hi = 2‖J‖_F` is per-instance-VALID but **NOT dimension-faithful** — its
`2N`-wide looseness is `O(N)`, so `hi/η` GROWS with `N` even when the true
`‖ΔΥ−Φ‖_cb/η` is bounded. Therefore T8 can assert per-instance `O(η)` but **cannot
test dimension-independence** — that is the §C12 trap (a dim-dependent proxy fakes
universality) and is deferred to F4.2's tight SDP. The T8 `C_max` gate must be
read as "the algorithm produces an `O(η)` output at this instance," nothing more.

### B.3 The four DUAL channels (T7/T8 read-off; `.tex:2159, 280-287`)

The dual of a UCP (Heisenberg/observable) map is a CPTP (Schrödinger/state) map
— the superoperator adjoint (CLAUDE.md "UCP vs CPTP"). The paper's `.tex:2159`:

```
Dec = Δ*   : B(H)* → B*        (states on B(H) → states on B; the DECODING channel)
Enc = Υ*   : B*    → B(H)*     (states on B → states on B(H); the ENCODING channel)
Dec·Enc = (ΥΔ)* ≈ 1_{B*}       (.tex:2159 with Δ,Υ in place of Γ∘Co_M, Δ)
Enc·Dec = (ΔΥ)* ≈ Φ*           (the approximate factorization Φ* ≈ Enc∘Dec)
```

> **Mind the direction (CLAUDE.md callout, ranked risk note).** `Δ: B→B(H)` is
> the observable-encode; its dual `Dec = Δ*: B(H)*→B*` is the state-DECODE. `Υ:
> B(H)→B` is the observable-decode; its dual `Enc = Υ*: B*→B(H)*` is the
> state-ENCODE. So the dual SWAPS the words encode/decode relative to the
> observable maps — `Dec` is the dual of the *encode* `Δ`. Getting this backwards
> silently transposes everything. The paper's `.tex:2159` is the authority:
> `Dec=Δ*`, `Enc=(Γ Co_M)*` there; our `Υ` plays the role of `Γ Co_M`, so
> `Enc=Υ*`.

**How to read off the duals (constructive, no new construction).** The HS-adjoint
of a linear map is the conjugate-transpose of its superoperator matrix; for a UCP
map in Kraus form `Φ(X)=Σ_a K_a† X K_a` the dual is the Kraus-adjoint
`Φ*(ρ)=Σ_a K_a ρ K_a†` (the Schrödinger Kraus, the SAME `K_a` with the dagger
moved). So:

- If we have Kraus `{D_a}` of `Δ` (extract via `aic_ucp_choi_to_kraus_latd` on
  `Choi(Δ)`), then `Dec = Δ*` has Schrödinger action `Dec(ρ) = Σ_a D_a ρ D_a†`.
  Equivalently and more cheaply: the dual is the **superoperator transpose** of
  the Convention-A Choi (Choi-matrix of `Φ*` is the SWAP-conjugate of `Choi(Φ)`),
  so F4 need not extract Kraus — it can expose the duals as the transposed
  superoperators and verify via the Choi.
- For F4's PURPOSES the duals do not need a new struct: the verification quantities
  are `‖Dec·Enc − 1_{B*}‖_⋄` and `‖Enc·Dec − Φ*‖_⋄`, and by the duality
  `‖Λ*‖_⋄ = ‖Λ‖_cb` (Watrous; the SAME duality the O2 rect certifier uses,
  `aic_cbnorm.h:89`), these EQUAL `‖ΥΔ−1_B‖_cb` and `‖ΔΥ−Φ‖_cb` — **already
  computed in B.1/B.2**. So the dual-channel assertions are:

```
‖Dec·Enc − 1_{B*}‖_⋄ = ‖ΥΔ − 1_B‖_cb ≤ hi_ud      (B.2 bound)
‖Enc·Dec − Φ*‖_⋄     = ‖ΔΥ − Φ‖_cb   ≤ hi_du       (B.2 bound)
```

and the **η=0 exactness** `Dec·Enc = 1_{B*}`, `Enc·Dec = Φ*` follows from the
η=0 oracle (`ΥΔ=1_B`, `ΔΥ=Φ`) by the same duality (`J=0 ⟹ J*=0`).

**Recommendation for the dual API.** F4 should expose a small read-off that
returns the dual channels as Kraus (`Dec` from `Choi(Δ)` via
`aic_ucp_choi_to_kraus_latd`; `Enc` from `Choi(Υ)`) so a downstream caller can
APPLY the channels to density matrices (the deliverable is "an approximate
factorization of `Φ*`"). The VERIFICATION, though, goes through the cb-norm of the
observable maps (no separate diamond-SDP on the duals needed — the duality makes
them identical, and computing them twice would be wasteful). Add a cheap CPTP
sanity check on the read-off Kraus: `Dec` trace-preserving ⟺ `Σ_a D_a† D_a = 1_H`
= `Δ` unital (already F2-asserted), and likewise `Enc` TP ⟺ `Υ` unital. So the
dual TP-ness is FREE from `Δ`,`Υ` unitality — assert it as a counted check for
visibility (mutation-proof: a non-unital `Δ` would make `Dec` non-TP).

---

## C. F4.2 (Julia+MOSEK) — the FAITHFUL canary

F4.2 is the universality-certification refinement. It is the ONLY part that can
faithfully test dimension-independence (the eig-free bound's `O(N)` looseness
cannot — §B.2, §C12). It mirrors the O2 fixture dance exactly.

### C.1 The tight `‖ΔΥ−Φ‖_cb` via `aic_cbnorm_certify` (square N case)

`ΔΥ−Φ` is a **square HP self-map on `M_N`**, so it uses the SQUARE self-map
certifier `aic_cbnorm_certify` (NOT the rect one) with the **(2/N) normalization,
Tr_2 = `partial_trace_right`** convention already pinned in `aic_cbnorm`
(`aic_cbnorm.h:53-74`). It needs BOTH a MAX-primal `(X,P,Q)` and a MIN-dual
`(Y0,Y1)` feasible point for `J_DelUps` (the self-map SDP has `P+Q=I` so it
carries both bounds; cf. `aic_cbnorm.h:58-61`). The committed point gives a
rigorous TIGHT two-sided ball `lo ≤ ‖ΔΥ−Φ‖_cb ≤ hi`.

For `‖ΥΔ−1_B‖_cb` (also a square self-map, on `M_{n_B}`): same `aic_cbnorm_certify`
with `J_UpsDel` and `n = n_B`.

### C.2 The shim (C builds J, exposes to Julia) — mirror `aic_opspace_shim`

Add `src/aic_factorize_shim.c` + `include/aic_factorize_shim.h`, a flat-double
ABI mirroring `aic_opspace_shim.h` (no FLINT type crosses the boundary):

```
int aic_factorize_choi_shim_d(
    double *jdu_re, double *jdu_im,    /* J_DelUps = Choi(ΔΥ)−Choi(Φ), N²×N², [p*sz+q] sz=N² */
    double *jud_re, double *jud_im,    /* J_UpsDel = Choi(ΥΔ−1_B),    n_B²×n_B², [p*sz+q] sz=n_B² */
    int *out_N, int *out_nB, int *out_dimB,
    double *out_eta,                   /* LENGTH-2: [0]=eta_proxy, [1]=eps_used (sentinel like O2) */
    const double *kraus_re, const double *kraus_im,
    int n, int r, double eps, int prec);
```

Internally it runs the EXACT C cores the C test exercises:
`kraus_from_flat → aic_assoc_ecstar_from_phi → aic_cstar_build → aic_factorize_build
→ aic_factorize_delta_build → aic_factorize_upsilon_build`, then assembles
`J_DelUps` (the `E_pq` loop, §A.1) and `J_UpsDel` (§A.2 route (i)), writes them
flat row-major. The `eps` sentinel convention is IDENTICAL to O2's
(`aic_opspace_shim.h:69-76`): `eps=0` for η=0 oracles, `eps<0` derives
`aic_ecstar_defect_assoc` internally so the C certifier test (rebuilding the same
pipeline) gets the IDENTICAL `J`. **CAVEAT (§C13(c) / §C11):** for the multi-block
η>0 fixture pass `eta`, NOT the ~700× smaller `eps_assoc`, as the build's eps scale
(the shim sentinel handles this if we add an "eps = eta" mode, or pass `eta`
explicitly as a positive `eps`).

> **Buffer sizing.** `J_DelUps` is `N²×N²` (length `N⁴`); `J_UpsDel` is `n_B²×n_B²`
> (length `n_B⁴`). Since `n_B ≤ N`, allocate both at `N⁴`. (The shim asserts the
> shapes, Rule 4.)

### C.3 The Julia generator — mirror `tools/gen_fixtures_opspace_o2.jl`

Add `tools/gen_fixtures_factorize_f4.jl`:

- `dlopen build/libaic.so`, `ccall aic_factorize_choi_shim_d` to get `J_DelUps`,
  `J_UpsDel`, `N`, `n_B`, `eta`, `eps_used` for each corpus channel.
- For each square `J` (both `J_DelUps` at size `N²` and `J_UpsDel` at size `n_B²`):
  solve the Watrous **self-map** diamond-SDP via `src/sdp.jl` (the `diamond_primal`
  / `diamond_dual` SELF-MAP form, with the (2/n) normalization — NOT the rect
  `diamond_*_rect`). Emit BOTH the MAX-primal `(X,P,Q)` and MIN-dual `(Y0,Y1)`
  feasible points + `eta_ref` = the SDP optval = the tight cb-norm.
- **POISON GUARDS before emitting** (mirror `gen_fixtures_opspace_o2.jl:196-213`):
  recompute the bound from the emitted `(X,P,Q,Y0,Y1)` and assert it reproduces
  `eta_ref` to ~1e-6; assert primal ≈ dual (strong duality); for η=0 oracles assert
  `eta_ref < 1e-6` (the η=0 cb-norm is ZERO here — unlike O2 where the complete-
  isometry oracle gives 1, the factorization defect `ΔΥ−Φ` is ZERO at η=0; NOTE
  this difference loudly in the generator). DROP any channel that fails to
  build/solve (no fake fixture).
- Emit `tests/fixtures_factorize_f4.inc.h`: per channel, the Kraus (to rebuild the
  whole pipeline in C), `eps_used`, dims `N, n_B, dim_B`, `eta`, and per quantity
  (`du` = DelUps, `ud` = UpsDel) the Convention-A Choi `J`, the MAX-primal
  `(X,P,Q)`, the MIN-dual `(Y0,Y1)`, and `eta_ref`. (Structure mirrors the O2
  fixture struct.)

The C consumer (`tests/test_factorize.c` T9, or a new `test_factorize_f4.c` TU)
rebuilds the pipeline from the committed Kraus + `eps_used`, ASSERTS the rebuilt
`J` reproduces the committed one to ~1e-9 (the O2 `jdiff` sanity, so the committed
feasible point is feasible for the assembled `J`), then calls `aic_cbnorm_certify`
and asserts the tight ball matches `eta_ref`. **`make test` stays Julia-free.**

### C.4 The dimension-independence canary (FINDINGS §D2 ROBUST form)

The canary measures `C = ‖ΔΥ−Φ‖_cb / η` (the tight SDP cb-norm, NOT the eig-free
bound) over a SWEEP of multi-block η>0 fixtures at INCREASING dimension, and
asserts:

- (i) **bounded abs-max:** `C < C_max` over the whole sweep (a generous fail-loud
  gate; the magnitude may be larger than naive — see risk #4).
- (ii) **no upward trend:** the halves-ratio `mean(C | hi-half dim) / mean(C |
  lo-half dim) ≤ 1.25` (the §D2 robust form: the halves aggregate dilutes a
  single hard-geometry spike, but a genuine `c=O(dim)` law survives the
  aggregation). Report the least-squares slope `β` too (should be ≤ 0).

**The sweep fixtures: `make_mixconj_blocks(d, m, t)`** (`test_cstar_build.c:170`,
the multi-block η>0 fixture, §C13(c)). It yields `≥2` equivalence classes at η>0
(the ONLY fixture exercising the multi-class merge + the per-block-Pauli-vs-whole-B
distinction at `m≥2 ∧ η>0`, the F2→F3→F4 coverage debt). The sweep increases `d`
(hence `N`, `n_B`, `dim_B`) at a few points: e.g. `(4,2,t)`, `(5,2,t)`, `(6,2,t)`,
`(6,3,t)`, `(7,3,t)` — chosen so `n_B`/`dim_B` grows while staying in the build
basin. **§C11 eps-scale caveat:** pass `eta` (≈1.6e-2), NOT the ~700× smaller
`eps_assoc` (≈2e-5), as the build's eps argument — else the Stage-1 errreduce C0
gate fires (`10·eps < the true O(η) inclusion defect`). The shim must honor this.

**MUTATION-PROOF (Rule 7, the §D2 discipline):** inject a `C_flat·(dim/dim_min)`
(genuine `c=O(dim)`) factor into the measured `C` and confirm the trend arm goes
RED (halves-ratio → ~1.80 > 1.25); inject the literal `C·(dim/dim_min)` and confirm
the absolute arm goes RED (`abs-max > C_max`). This is the SAME mutation the
`cstar_build` T2b canary uses (FINDINGS §D2, mutation-proven).

**WHY the eig-free bound CANNOT do this canary (the §C12 trap — state explicitly).**
`aic_cbnorm_eigfree_ball_choi` returns `hi = 2‖J‖_F`. The ratio `‖J‖_F / ‖J‖_⋄`
for an `N²×N²` Choi can be as large as `O(N)` (the Frobenius norm sums `N²`
singular values; the diamond norm sees the worst ampliated stretch). So
`hi/η = 2‖J‖_F/η` carries a SPURIOUS `O(N)` factor that GROWS with dimension even
when the TRUE `‖ΔΥ−Φ‖_cb/η` is bounded. Using the eig-free bound in the canary
would either (a) FAKE a `c=O(N)` law (false positive: the canary RED-fires on a
sound, dimension-independent algorithm — the §C12 "the Frobenius proxy fakes
non-universality") or (b) MASK a genuine `c=O(N)` law under an even-larger looseness
(false negative). Either way it is a "test that cannot fail meaningfully" — the
canary REQUIRES the faithful tight SDP cb-norm. This is the load-bearing reason
F4.2 exists.

---

## D. THE ONE EMPIRICAL PROBE F4 NEEDS (before fixing the §A.2 convention)

Route (i) (ambient `M_{n_B}` Choi) for `‖ΥΔ−1_B‖_cb` rests on the claim that `Δ`
zeroes the off-block-diagonal entries of its `M_{n_B}` input (so ambient = in-B
cb-norm). `aic_dhom_v_apply` reads only the matrix-unit coordinates of the
block-diagonal `X`, so this SHOULD hold — but the η=0 oracle is BLIND to it
(`ΥΔ=1_B` regardless). **The probe:** at an η>0 fixture, build `J_UpsDel` BOTH ways
— route (i) over the full `n_B²` ambient matrix units, and route (ii) over the
`dim_B` in-B matrix units — and compare the eig-free `hi` (or, in F4.2, the tight
SDP `eta_ref`). If they AGREE, route (i) is sound and is the recommended definition
(it plugs into `aic_cbnorm_certify` directly). If they DIFFER, route (ii) is the
correct cb-norm and the certifier must consume the `dim_B`-basis Choi (which is
NOT square `n_B²` — it would need the rect machinery or a basis-aware variant). I
expect agreement (the off-block-diagonal columns of `Choi(ΥΔ−1_B)` should be zero
because `Δ` drops them), but it MUST be measured because it is exactly the kind of
ambient-vs-operator-space subtlety §C12 warns is silently vacuous. This is a
**definition-correctness probe, NOT a convention pin** — the conventions (2/n,
Tr_2 direction) are inherited; only WHICH Choi (full `n_B²` vs `dim_B`-basis) is
open, and the probe resolves it.

---

## E. NEW convention risk? — NO. Everything inherited / safe.

Per the orchestrator's §E ask, I checked every convention F4 touches against
what is ALREADY pinned in `aic_cbnorm`/O2 (FINDINGS §C12.O2-PIN, the empirical pin):

| Convention | F4 use | Pinned where | Status |
|---|---|---|---|
| Convention-A Choi (input-major, conj on first factor) | `Choi(ΔΥ)`, `Choi(ΥΔ−1_B)`, `Choi(Φ)` | `aic_ucp.h:23-43`; F2 `aic_factorize_delta_block_choi` already uses it | INHERITED |
| Square self-map diamond-SDP **(2/n)** normalization, Tr_2 = `partial_trace_right` | `aic_cbnorm_certify(J_DelUps, N)`, `(J_UpsDel, n_B)` | `aic_cbnorm.h:53-74` (the self-map `P+Q=I` form) | INHERITED — `ΔΥ−Φ`, `ΥΔ−1_B` ARE square HP self-maps (exactly like `Φ²−Φ`) so the (2/n) self-map convention applies, NOT the rect factor-1 one |
| eig-free bracket `[‖J‖_F/n, 2‖J‖_F]` | T8 per-instance upper | `aic_cbnorm.h:36-44` | INHERITED |
| Kraus Heisenberg convention, dual = Kraus-adjoint | duals `Dec=Δ*`, `Enc=Υ*` | `aic_ucp.h:18-21` (the unital test pins direction) | INHERITED |
| diamond-of-adjoint duality `‖Λ*‖_⋄ = ‖Λ‖_cb` | dual-channel verification = observable cb-norm | `aic_cbnorm.h:89` (O2 rect uses it) | INHERITED |
| shim flat-double ABI, `eps` sentinel | F4.2 shim | `aic_opspace_shim.h` / `aic_ucp_shim.h` | INHERITED |
| Julia fixture dance (POISON guards, committed `.inc.h`) | F4.2 generator | `gen_fixtures_opspace_o2.jl` | INHERITED |

**The one thing that is NOT a pre-pinned convention** is the §A.2 definition choice
(ambient `M_{n_B}` Choi vs `dim_B`-basis Choi for the direct-sum `ΥΔ−1_B`). I
classify this as a **definition choice with a one-shot empirical probe** (§D), NOT
a 3-way design competition like F3's `⊗1_F` ordering got. Reason: it is not an
ambiguous SIGN/direction with two plausible answers and an O(1)-vs-O(η)
discriminant; it is a "do these two Choi give the same diamond norm" check whose
answer I can predict (yes, because `Δ` drops off-block-diagonal input) and which
the η>0 probe confirms in one run. If the probe SURPRISES (they differ), THEN it
escalates to a design decision — but I do not anticipate that, and pre-committing
to a competition would be over-engineering. **Recommended: route (i) + the §D
probe as a committed teeth-bearing check** (assert ambient ≈ in-B at η>0).

> **One subtlety to FLAG (not a new convention, but a sharp edge):** the §C13(c)
> distinction — the per-block Pauli twirl (`R_j`) vs the whole-B 1-design join — is
> the F3 `lem_RC` internal; F4's `J_UpsDel` does not re-touch it, but the
> multi-block η>0 fixture (`make_mixconj_blocks`) is the FIRST place `m≥2 ∧ η>0`
> exercises the full `ΥΔ` round-trip. So the F4.2 canary fixture corpus is ALSO
> the coverage-debt repayment for §C13(c). Two birds.

---

## F. Risks (ranked) + file/LOC plan + test plan

### F.1 Ranked risks + mitigations (each constrained by the η=0 oracle)

**Risk #1 — the `ΥΔ−1_B`-on-direct-sum-`B` cb-norm definition (§A.2).**
The map lives on `B = ⊕M_{d_l}`, not `M_N`; getting the cb-norm definition wrong
(per-block route (iii), or ambient-vs-in-B mismatch) gives a "test that cannot
fail" (§C12). *Mitigation:* route (i) (ambient `M_{n_B}` Choi → `aic_cbnorm_certify`
with `n=n_B`) + the §D empirical probe (assert ambient ≈ in-B at η>0, a counted
teeth-bearing check). *η=0 constraint:* the oracle gives `J_UpsDel=0` for BOTH
routes, so the oracle ALONE cannot pin route (i) vs (ii) — the §D η>0 probe is
REQUIRED (the oracle is blind here, exactly the kind of gap §C12 catches).

**Risk #2 — the MOSEK fixture corpus + the shim (F4.2).**
The shim must rebuild the WHOLE pipeline (assoc → cstar → factorize delta+upsilon)
and assemble TWO square Choi; a bug there produces a `J` the committed feasible
point is not feasible for, silently weakening the certificate. *Mitigation:* the
O2 `jdiff < 1e-9` sanity (rebuilt `J` reproduces committed `J`) in the C consumer;
the generator's POISON guards (strong duality, recompute bound from emitted point);
DROP-not-fake on build failure. *η=0 constraint:* the η=0 oracle fixtures
(`eps=0` path) give `eta_ref ≈ 0` exactly — a generator guard `eta_ref < 1e-6`
catches a shim that mis-assembles `J` (it would not be ~0).

**Risk #3 — the dim-canary fixture (multi-block η>0, §C11 eps-scale).**
`make_mixconj_blocks` has `eps_assoc ≈ 700×` below `η`; passing `eps_assoc` as the
build's eps fires the errreduce C0 gate (build aborts) — and silently using the
wrong scale would make `C = cb/η` meaningless. *Mitigation:* pass `eta` as the
build's eps (the §C11 caveat, honored in the shim's eps mode); the generator
DROPS+notes any fixture that fails to build. The dim sweep must use ≥4 points to
make the halves-ratio robust (§D2). *η=0 constraint:* none directly — this is an
η>0-only fixture; the η=0 oracle uses the SINGLE/multi-block `block_cond_exp`
fixtures.

**Risk #4 — the composite constant magnitude (NOT smallness).**
FINDINGS notes `‖v⁻¹‖_cb` can be `1.535` for oblique `A` (O2 mixconj), and the
composite `C = c1+c2+c3+O(η)` (research §5) inherits this — so `C` may be LARGER
than a naive reader expects (the canary must confirm dimension-INDEPENDENCE, NOT
smallness). *Mitigation:* set `C_max` generously (e.g. 40 in T8's eig-free gate;
in F4.2's tight canary, base `C_max` on the MEASURED first-sweep value with a
~1.5× margin, the `cstar_build` T2b posture — measured `c ∈ [0.25, 3.27]` there,
so `C_max ~ 5` for the tight quantity, but F4's COMPOSITE may run higher; PIN it
from the measured sweep, do not pre-guess). The canary's load-bearing assertion is
the no-upward-trend halves-ratio, NOT the absolute magnitude. *η=0 constraint:* the
oracle gives `C=0/0` — undefined; skip the canary at η=0 (it is an η>0 sweep).

**Risk #5 (minor) — `ΔΥ` Choi is `N²×N²` (memory/time).**
For `N` up to ~10, `J_DelUps` is `100×100` and the self-map SDP is `10⁴` variables
— feasible (O2's mixconj `N=6` Choi is `36×36`, solved fine). The canary sweep top
end (`d=7` → `N=7`) gives a `49×49` Choi, the self-map SDP `2401`-dim — within
MOSEK reach. *Mitigation:* cap the sweep at the proven build basin (`d ≤ 7-ish`);
the eig-free T8 (no SDP) runs at any `N`.

### F.2 Source split (≤200 LOC/file, Rule 10)

- **`src/aic_factorize_verify.c`** (NEW): the composed-map Choi builders +
  eig-free bound + dual read-off. Public functions:
  - `aic_factorize_choi_delups(acb_mat_t J, const aic_factorize *F, slong prec)`
    — `J_DelUps = Choi(ΔΥ) − Choi(Φ)`, `N²×N²` (the `E_pq` loop + `aic_ucp_kraus_to_choi`).
  - `aic_factorize_choi_upsdel(acb_mat_t J, const aic_factorize *F, slong prec)`
    — `J_UpsDel = Choi(ΥΔ − 1_B)`, `n_B²×n_B²` (route (i), the ambient `M_{n_B}`
    `E_pq` loop minus the `1_B` block-diagonal).
  - `aic_factorize_eigfree_delups(arb_t lo, arb_t hi, const aic_factorize *F, slong prec)`
    and `_upsdel` — wrap `aic_cbnorm_eigfree_ball_choi`.
  - `aic_factorize_dec_kraus`, `aic_factorize_enc_kraus` — read off `Dec=Δ*`,
    `Enc=Υ*` as Kraus (via `aic_ucp_choi_to_kraus_latd` on `Choi(Δ)`/`Choi(Υ)`).
    (Choi of `Δ`,`Υ` are built by the same `E_pq` pattern.)
  This file may push 200 LOC; if so, split the dual read-off into
  `src/aic_factorize_dual.c`.
- **`include/aic_factorize.h`** (EXTEND): the F4 declarations above + an
  `aic_factorize` field note (no new owned state needed — the Choi/duals are
  built on demand into caller buffers).
- **`include/aic_factorize_shim.h` + `src/aic_factorize_shim.c`** (NEW, F4.2):
  the flat-double `aic_factorize_choi_shim_d` (§C.2). ~150 LOC, mirrors
  `aic_opspace_shim.c`.
- **`tools/gen_fixtures_factorize_f4.jl`** (NEW, F4.2): mirrors
  `gen_fixtures_opspace_o2.jl`. ~250 LOC (tools exempt from Rule 10).
- **`tests/fixtures_factorize_f4.inc.h`** (NEW, GENERATED): committed feasible
  points.

### F.3 Test plan

In `tests/test_factorize.c` (extend; T7/T8 are Julia-free) and a new
`tests/test_factorize_f4.c` (T9, consumes the committed MOSEK fixtures — mirrors
`test_opspace_o2.c`):

- **T7 — η=0 ORACLE (Julia-free).** On `block_cond_exp(4,2)`, `block_cond_exp(6,3)`
  (multi-block η=0): `‖J_DelUps‖_F < 1e-10`, `‖J_UpsDel‖_F < 1e-10`, eig-free
  `hi_du < 1e-9`, `hi_ud < 1e-9`; the un-ampliated `UpsDel2` spot-check
  `< 1e-10`; dual TP-ness (`Dec`,`Enc` trace-preserving = `Δ`,`Υ` unital).
  MUTATION: drop the `Choi(Φ)` subtraction → `J_DelUps = Choi(ΔΥ)` (a UCP Choi,
  `‖·‖_F = O(1)`) → RED.
- **T8 — η>0 per-instance eig-free (Julia-free).** On `make_mixconj` (single-block)
  AND `make_mixconj_blocks` (multi-block, §C13(c) coverage): `hi_du ≤ C_max·η`,
  `hi_ud ≤ C_max·η` (`C_max` generous, e.g. 40); print `hi_du/η`, `hi_ud/η`.
  The §D ambient-vs-in-B probe: `‖J_UpsDel^{(i)} − J_UpsDel^{(ii)}‖_F < 1e-12`
  (route definitions agree) — a counted teeth-bearing check.
  MUTATION: route (iii) per-block max → would miss cross-block leakage (verify by
  injecting an artificial cross-block term).
- **T9 — MOSEK-certified tight bound + canary (consumes committed fixtures).**
  Rebuild the pipeline from committed Kraus + `eps_used`; assert rebuilt `J`
  reproduces committed (`< 1e-9`); `aic_cbnorm_certify(J_DelUps, N)` gives a tight
  ball matching `eta_ref` (`hi ≥ eta_ref − 1e-7`, `hi − eta_ref ≤ 1e-4`); same for
  `J_UpsDel`. The CANARY: over the `make_mixconj_blocks` sweep, `C = eta_ref/η`
  bounded (`< C_max_tight`) + halves-ratio `≤ 1.25` + slope `≤ 0`. MUTATIONS
  (§D2): `C·(dim/dim_min)` → absolute arm RED; `C_flat·(dim/dim_min)` → trend arm
  RED.

### F.4 Does F4.1 alone CLOSE the constructive headline? — YES.

F4.1 (T7 + T8 + the duals, all Julia-free) produces:
- The certified `Δ: B→B(H)` (UCP), `Υ: B(H)→B` (UCP), and `B = ⊕M_{d_l}` (genuine
  C*) — the headline objects `th_factorization` asserts EXIST.
- The η=0-EXACT oracle (`ΔΥ=Φ`, `ΥΔ=1_B` to machine) — ladder rung 3.
- A per-instance RIGOROUS certified `O(η)` upper bound (`‖ΔΥ−Φ‖_cb ≤ 2‖J‖_F`,
  `‖ΥΔ−1_B‖_cb ≤ 2‖J‖_F`), valid at every instance, with `hi/η` a finite
  per-instance constant.
- The four dual channels (`Enc=Υ*`, `Dec=Δ*`, `Dec·Enc≈1_{B*}`, `Enc·Dec≈Φ*`),
  the deliverable "approximate factorization of `Φ*`."

That IS the constructive content of `th_factorization` — an algorithm that, for any
η-idempotent `Φ`, builds the factorization and certifies it meets the `O(η)` bound
PER INSTANCE. **F4.2 is the universality CERTIFICATION refinement** (the
dimension-INDEPENDENT, faithful `C`), which the project standard (FINDINGS §D2/§D4,
the `c_0`/opspace O1 precedent) treats as a per-instance-measured + canary-asserted
quantity with the ANALYTIC constant deferred to a research bead (chained after
`aic-1bc`). So: F4.1 closes the headline; F4.2 closes the universality canary; the
analytic composite `C` remains the one filed-research escalation (research §5,
already an open bead posture), NOT a wall.

---

## G. Summary of recommendations

1. **Adopt the F4.1/F4.2 split** as proposed. F4.1 = Julia-free (η=0 oracle +
   eig-free per-instance upper + duals), CLOSES the constructive headline. F4.2 =
   Julia+MOSEK (tight `aic_cbnorm_certify` + the dim-independence canary), the
   universality certification.
2. **`ΥΔ−1_B` cb-norm = route (i):** the ambient `M_{n_B}` Convention-A Choi,
   `n_B²×n_B²`, fed to the SQUARE self-map `aic_cbnorm_certify` with `n=n_B` and
   the (2/n) convention. Verify `‖ΥΔ−1_B‖_cb` (cleaner) as the headline; add an
   un-ampliated `UpsDel2` spot-check as a cross-check.
3. **One empirical probe (§D):** at η>0, confirm route (i) (ambient) ≡ route (ii)
   (in-B basis) Choi give the same diamond norm (a definition-correctness check the
   η=0 oracle is blind to). NOT a convention pin, NOT a 3-way competition.
4. **No NEW convention needs pinning** — every convention is inherited from
   `aic_cbnorm`/O2 (§E). The (2/n) self-map convention (not the rect factor-1)
   applies, because `ΔΥ−Φ` and `ΥΔ−1_B` are square HP self-maps.
5. **The canary (§C.4):** `C = ‖ΔΥ−Φ‖_cb/η` (TIGHT SDP, never eig-free) over a
   `make_mixconj_blocks` dim sweep; bounded abs-max + no-upward-trend halves-ratio
   (§D2 robust form), mutation-proven against a `c=O(dim)` injection. The eig-free
   bound CANNOT do this (its `O(N)` looseness fakes/masks a dim-dependent `C` — the
   §C12 trap).
