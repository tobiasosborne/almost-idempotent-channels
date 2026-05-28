# MODULE_PLAN.md — almost-idempotent-channels

Module plan for implementing arXiv:2405.02434 (Kitaev) as constructive
finite-dimensional algorithms. Read `CLAUDE.md` first (the mandate, the Laws,
the rules); read the relevant `paper/shards/shard-*.md` before writing a module.
This file is the DAG and the result→module map; it is durable and is updated
when the structure changes.

> **Registered in beads** (2026-05-28, prefix `aic`): every module/research item
> below is a tracked issue with the build-order DAG encoded as blocking
> dependencies. `bd ready` (entry point `T-build`), `bd dep tree <id>`.
>
> **Audition, not decision (Law 4).** Algorithm names in the tables/diagram below
> (Denman–Beavers, Newton, SVD-polar, projection Route A, …) are **candidates to
> audition** via red-green TDD against correctness *and* performance benchmarks —
> not chosen implementations. Keep candidates alive until a benchmark retires
> one. Nothing here is presumed fit.

## Goal recap

Given a finite-dimensional UCP map $\Phi$ with $\|\Phi^2-\Phi\|_\cb\le\eta$,
produce: (1) the associated $\eps$-$C^*$ algebra $\calA$ ($\eps=O(\eta)$);
(2) a genuine $C^*$ algebra $\calB$ and an $O(\eps)$-isomorphism $v:\calB\to\calA$
with a **universal, dimension-independent** constant; (3) UCP maps realizing $v$
and $v^{-1}$, giving the approximate factorization of the channel $\Phi^*$ into
decode∘encode. Tight C + FLINT/arb (certified) cores; Julia `ccall` package.

**The mandate (see CLAUDE.md):** the paper's proofs are largely
non-constructive / infinite-dim; we implement constructive finite-dim
algorithms whose outputs meet the paper's $O(\eps)$ bound. The bound is the
spec; the arb path certifies it.

## Shard index (the per-result analysis lives here)

| Shard | Region | Lines | Drives modules |
|---|---|---|---|
| A | §1 Motivation, §2 Main thm + strategy | 251–492 | (problem statement; defs) `ecstar` |
| B | §3 Analytic tools, §4 ε-Banach basics | 493–689 | `funcalc`, `contraction`, `unitfix` |
| C | §5 Approximate unitary group | 690–914 | `unitary` |
| D | §6 Nontrivial projection, §7 Subspaces | 915–1189 | `projection`, `corner` |
| E | §8 δ-homs, §9 Main-theorem proof | 1190–1446 | `dhom`, `errreduce`, `cstar_build` |
| F | §10 Tensor extensions | 1447–1561 | `opspace` |
| G | §11 UCP map details, idemp-structure proof | 1562–2161 | `ucp`, `idemp_structure` |
| H | §12 Almost-idempotent maps, factorization | 2162–2900 | `assoc_ecsa`, `factorize` |

## Result → module map

| Result (`.tex` line) | Module | Constructive? |
|---|---|---|
| fc Taylor / `|X|` / `sgn(X)` (514–522) | `funcalc` | yes (Denman–Beavers / Newton) |
| `prop_P` projection-extract θ(2P−I) (524) | `funcalc` | yes |
| `lem_invfun` contraction (564) | `contraction` | **yes — the iteration is the algorithm** |
| `prop_unit` exact unit (672) | `unitfix` | yes |
| ε-/δ- definitions (407, 443) | `ecstar` (types) | n/a (data model) |
| `lem_U_delta`,`lem_gV`,`prop_polar` (699,758,809) | `unitary` | yes (SVD/Newton polar + contraction charts) |
| `lem_nontriv_projection` (931) | `projection` | **proof non-constructive → Route A** |
| `prop_H-group` (971) | — (proof-only) | not implemented (existence-proof internal) |
| `lem_alpha`,`lem_PQ_Hilb`,`lem_PQR`,`lem_1d_proj` (1086…1179) | `corner` | yes (explicit compressions/ranks) |
| `prop_delta_hominc`,`lem_approx` (1194,1224) | `dhom` | yes (Newton on mult. defect) |
| `cor_improvement` error reduction (1317) | `errreduce` | **yes — explicit B-diagonal, no n-blowup** |
| `lem_merging`,`cor_merge_sum`,`lem_add_dim`,`lem_extension` (1325…1378) | `cstar_build` | yes |
| proof of `th_main` (1414) | `cstar_build` (master loop) | yes (modulo `projection`) |
| `def:opspace`, `prop_inc_ext`,`lem_approx_ext`,`th_main_ext` (1453…1538) | `opspace` | yes (modulo cb-truncation N) |
| Choi/Stinespring/Kraus, `prop_KLHG` (1568…1635) | `ucp` | yes (matrix constructions) |
| `lem_carrier`,`cor_carrier` (1724,1731) | `ucp` | yes (Im V / partial trace) |
| `lem_idemp`, proof of `th_idemp_structure` (1916,2055), `prop_Gamma` (2106) | `idemp_structure` | yes |
| `th_almost_idemp` + Phi_assoc1/2 (2192) | `assoc_ecsa` | yes |
| `th_factorization` (2730), `lem_RC` (2840) | `factorize` | **partial — proof is an outline (see escalations)** |

## Module DAG (build order = topological)

Each module ships **two number paths** (double/LAPACK, arb/`acb`) wherever it
does real arithmetic, cross-checked double-vs-arb@prec=53. `≤200 LOC` per file
(split as needed). Arrows = "depends on".

```
L0  mat ──────────────────────────────────────────────────────────────────┐
      complex matrix; ‖·‖_op (SVD), ‖·‖_F; Hermitian eig (zheev/acb);      │
      SVD; Kronecker ⊗; partial trace; reshape. double + acb_mat.          │
        │             │                │                                    │
L1  funcalc      contraction          ucp                                   │
      sgn,|X|,θ,    Banach fixed-pt    Choi/Stinespring/Kraus;              │
      x^α; precond  solver (lem_invfun) conversions; Choi-PSD/unital       │
      ‖X²−I‖<1      x_n=g_y(x_{n-1})   checks; cb-norm; carrier             │
        │   │            │              (lem_carrier/cor_carrier)           │
        │   └──────┐     │                 │                                │
L2  ecstar      unitary  │            idemp_structure                       │
      A⊆M_N + ⋆-prod;  polar u(X),    th_idemp_structure algo:             │
      axiom-defect     μ,σ, charts     Φ idemp → M,A,Δ,Γ;                  │
      estimators;      (prop_polar)    lem_idemp, prop_Gamma               │
      unitfix(prop_unit)  │                                                 │
        │   │            │                                                  │
L3  projection ◄─────────┘         corner                                   │
      Route A: Herm eig +            Co_{P,Q}, S_{P,Q}, compressed          │
      gap + θ-snap +                 product; lem_alpha (block bijection);  │
      arb gap cert                   lem_PQ_Hilb/PQR/1d_proj; equiv classes │
      (lem_nontriv_projection)         │                                    │
        │                              │                                    │
L4      └──────────┬───────────────────┘                                    │
              dhom            errreduce                                      │
            mult. defect;     cor_improvement via explicit                  │
            lem_approx        Heisenberg–Weyl diagonal of B                 │
            (Newton)          (constructive; dimension-independent)         │
                  │                │                                        │
L5            cstar_build  ◄────────┘   ◄── projection, corner              │
                THE MASTER LOOP (proof of th_main):                         │
                B = partitioned index set + matrix units E_jk;              │
                Stage1 commutative skeleton → Stage2 extend blocks          │
                (lem_extension) → Stage3 merge (cor_merge_sum);             │
                cor_improvement after each step keeps δ_0 = c_0·ε           │
                  │                                                          │
L6            opspace  ◄────────────────────────────────────────────────────┘
                ampliations 1_{M_n}⊗(·); cb-norm; prop_inc_ext,
                lem_approx_ext; cstar_build_ext = th_main_ext
                  │
L7  assoc_ecsa ───┴──── factorize                          ◄── ucp, idemp_structure
      Φ̃=θ(2Φ−1);          7-step headline pipeline (shard H):
      A=Img Φ̃; ⋆;          regularize → A → (th_main_ext) B,v →
      Phi_assoc1/2;         CP-ize Δ (unitary 1-design) → decode Υ
      th_almost_idemp       (lem_RC) → unitalize → verify ‖ΔΥ−Φ‖_cb=O(η);
                            dual: Enc=Υ*, Dec=Δ*  (th_factorization)
```

## C ABI surface (sketch — `include/aic.h`)

Opaque handles; explicit precision argument on arb entry points; fail-fast
(`assert` + error return). Storage-layout contract documented in the header (the
`su2-fft` discipline).

```c
/* representations */
aic_ucp  *aic_ucp_from_kraus(const acb_mat_t *K, int r, int dimH);
aic_ucp  *aic_ucp_from_choi (const acb_mat_t C, int dimK, int dimH);
double    aic_ucp_idemp_defect_cb(const aic_ucp *Phi);      /* ‖Φ²−Φ‖_cb */

/* core constructions */
aic_ecstar *aic_assoc_ecstar(const aic_ucp *Phi, slong prec);   /* §12.1 */
aic_iso    *aic_main_iso    (const aic_ecstar *A, slong prec);  /* th_main(_ext): B, v */
aic_fact   *aic_factorize   (const aic_ucp *Phi, slong prec);   /* §12.2: Δ,Υ,B */

/* certified bounds (arb): returns an upper ball on the defect */
void aic_fact_channel_defect_cb(arb_t out, const aic_fact *F);  /* ‖ΔΥ−Φ‖_cb */
```

## Julia layer (`julia/AlmostIdempotentChannels.jl`)

Mirrors `../su2-fft/julia/SU2FFT.jl`: `ccall` wrappers + accessors, local
`libaic.so` build via `deps/build.jl`, gold-standard cross-check in
`test/runtests.jl`. High-level API: `UCPMap` (from Kraus/Choi), `idemp_defect`,
`associated_ecstar`, `main_isomorphism`, `factorize` → `(B, Δ, Υ)`, plus
`channel_defect` returning the certified `O(η)` ball.

## Arb-precision-critical modules (where `double` will lie to you)

1. `funcalc` — `sgn`/`|X|` near a singular `|X|` (small spectral gap); the
   `(L_{|X|}+R_{|X|})^{-1}` inverse (`.tex:615`).
2. `projection` — the spectral gap that certifies non-triviality is `~ε`
   (shard D); arb interval on eigenvalues is the certificate.
3. `errreduce` — cancellation producing the `O(ε²)` residual that drives the
   Newton/cohomology iteration.
4. `opspace` — cb-norm = sup over ampliations; balls keep the truncation honest.
5. `factorize` — composing many `O(η)` bounds; only arb gives a trustworthy
   end-to-end certificate.

## Cross-check strategy

Per-module ladder (CLAUDE.md §"cross-check ladder"): internal sanity → double
vs arb@prec=53 → exact special case → arb bound certification. Plus two global
oracles:

- **η=0 oracle (the cleanest ground truth):** an exactly idempotent Φ must drive
  the whole pipeline to the *exact* Choi–Effros / `th_idemp_structure`
  decomposition — a genuine `C^*` algebra with **zero** associativity/C*/unit
  defect and `v` an exact isomorphism. Every module on the path asserts its
  defect → 0 as η → 0.
- **round-trip:** `Dec∘Enc = (ΥΔ)^* ≈ 1_{B^*}` and `Enc∘Dec = (ΔΥ)^* ≈ Φ^*`,
  both certified `O(η)` in cb-norm.

## Open escalations (carry these; do not paper over)

1. **`lem_nontriv_projection` constructivization (resolved → verify).** Route A
   (Hermitian eigensolve + spectral gap + `θ(2P−I)` snap + arb gap cert, shard
   D). Open: a structural argument that an `Ω(1)` gap always exists for
   `dim A>1`; until proven, certify per-instance and escalate on failure.
2. **cb-norm truncation N (shard F).** "for all n" must be truncated; conjecture
   `n ≤ dim A` suffices. Needs proof or a certified-N procedure before `opspace`
   can claim a cb-bound.
3. **Universal constant `c_0` (shards A, E).** `cor_improvement`'s constant fixes
   the `O(ε)` constant in `th_main`. Extract it explicitly; tests must confirm
   the constant does **not** grow with `dim A` (the failure mode of the naive
   Haar route, `.tex:484`).
4. **`th_factorization` is an outline (shard H).** Steps 4–5 (CP-ization via a
   unitary 1-design; decode `Υ` via `lem_RC`) are prose, not a closed proof; the
   closure linking the constructions to the stated bounds (`DelUps`, `UpsDel2`)
   must be reconstructed before `factorize` is trustworthy.
5. **Composite `O(η)` constant (shard H).** Re-derive and compose the per-lemma
   constants along the full chain for a certified end-to-end
   `‖ΔΥ−Φ‖_cb ≤ C·η`.

## Suggested first vertical slice

Don't build bottom-up to completion before any integration. Build the **η=0
spine** first: `mat → funcalc → ucp → idemp_structure`, validated by the η=0
oracle (exact idempotent Φ → exact M,A,Δ,Γ; reconstruct Φ and compare). That
exercises the data structures and the cleanest theorem end-to-end, and gives a
fast regression base before the harder `projection`/`cstar_build` work.
