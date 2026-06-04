# PROFILING.md — almost-idempotent-channels

> The performance counterpart to `ALGORITHM.md`. Records, with concrete numbers,
> where the time goes, what each optimization lever is *measured* to buy, and the
> Law-4 audition slate per lever. House style: concrete numbers only, no
> marketing; measured vs estimated vs research-backed are labeled distinctly;
> external claims cite a source. Companion to bead epic **`aic-xv4`**.
>
> Campaign started 2026-06-04 (user directive: "go all out on performance").
> Methodology: 5 bench primitives + targeted FLINT/aic probes + a full serial
> `ctest` hot-path ranking + 4 deep-research reports (2 codebase sweeps, 2
> literature). Deep attribution (callgrind/DHAT/cachegrind) per the playbook in §7.

## 0. Hardware / build under test

- **CPU**: Intel Xeon E5-2699 v3 (Haswell), 18 physical cores / 36 threads,
  2.30 GHz, **AVX2 + FMA** (no AVX-512), L1d 32 KiB/core, L2 256 KiB/core,
  **L3 45 MiB shared**. Under **WSL2** (Microsoft kernel — `perf` hardware
  counters unavailable; valgrind is the deep profiler).
- **Numerical stack**: FLINT 3.0.1 (arb/acb inside it), **reference netlib BLAS +
  LAPACK** (the Debian `libblas3`/`liblapack3` default), LAPACKE.
- **Build**: CMake `RelWithDebInfo` = `-O2 -g`, **`-DNDEBUG` stripped** (342
  asserts stay live — a correctness choice, Rule 4), no `-march=native`, no LTO,
  single-threaded (`flint_set_num_threads` never called → default 1).

## 1. The two-worlds framing (the central finding)

There are two disjoint performance regimes with **different levers**. Conflating
them is the main trap.

| | **Arb path** (certified) | **Double path** (LAPACK anchor) |
|---|---|---|
| Backs | the slow drivers: cstar_build, opspace, errreduce, factorize, assoc | the fast gate, Julia default `certified_defect`, dim-sweep canaries |
| Arithmetic | pure FLINT ball arithmetic (GMP/MPFR limbs) — **never touches BLAS** | reference netlib BLAS/LAPACK |
| Dominant cost | certified Hermitian eig + matrix-sign on an **n²×n² superoperator** | small `zheev`/`zgesvd` (n≤32) |
| Lever | task-parallelism / **precision** / **arb algorithm** (NOT BLAS, NOT threads) | backend choice (but see §5: OpenBLAS is *not* a clear win here) |

**Consequence**: OpenBLAS — the obvious "reference BLAS is slow" lever — does
**not** touch the minute-scale drivers (they are arb-path, BLAS-independent), and
the audition (§5) shows it is not even a clear win on the double path at these
sizes. The real arb-path levers are task-level parallelism and precision.

## 2. Hot-path ranking (measured — full serial `ctest`, 2026-06-04)

Total suite wall time **1379 s** (23 min, serial). Top consumers:

| rank | test (driver exercised) | wall (s) | % of suite |
|---|---|---|---|
| 1 | `test_cstar_build` (`aic_cstar_build`, proof of th_main) | **327** | 24% |
| 2 | `test_opspace` (`aic_opspace_*` cb-norm via HOPM ampliations) | **139** | 10% |
| 3 | `test_errreduce` (`aic_errreduce`, cor_improvement) | **105** | 8% |
| 4 | `test_factorize` (`aic_factorize` pipeline) | **96** | 7% |
| 5 | `test_adversarial` | 52 | 4% |
| 6 | `test_shim_factorize_artifacts` | 35 | 3% |
| 7 | `test_dhom` | 25 | 2% |
| 8 | `test_projection` | 24 | 2% |
| 9 | `test_shim_assoc_iso` | 23 | 2% |

**Top 4 = 667 s = 48% of the suite.** These are the targets. (`test_opspace` at
#2 was under-emphasized in prior handoffs — the HOPM cb-norm ampliation is a
first-class cost.)

> NOTE (Rule 2): these per-test wall times are the *test harness* cost, a proxy
> for the library hot paths. True per-function call counts inside the recursive
> drivers need callgrind (§7) — static call-tree counting undercounts the
> th_main master-loop recursion onto S_P wrapper subalgebras.

> KNOWN RED: `test_opspace_o2` (T2 mixconj) aborts in the full suite — committed
> SDP dual infeasible for the assembled J (diff 1.728). Bead **`aic-nsb`**.
> Pre-existing (not a perf change); the fast gate stays green.

## 3. Measured primitive costs

### 3a. Bench baselines (prec=53, the cheap rung)

| op | n=4 | n=8 | n=16 | n=32 |
|---|---|---|---|---|
| certified Hermitian eig (arb) | 0.23 ms | 1.43 ms | 8.6 ms | 63 ms |
| operator norm (→ herm_max_eig) | 0.22 ms | 1.43 ms | 9.1 ms | 61 ms |
| sgn Newton–Schulz (arb) | 0.45 ms | 2.55 ms | 16.7 ms | — |
| acb 4×4 mul | 3.6 µs | — | — | — |

### 3b. Double vs arb gap (certified eig) — and it *widens* with n

| n | double eig | arb@53 eig | ratio |
|---|---|---|---|
| 4 | 4.1 µs | 234 µs | 58× |
| 16 | 70 µs | 8.4 ms | 120× |
| 32 | 424 µs | 63 ms | **149×** |

### 3c. Precision is a real cost knob (opnorm n=16, 1 thread)

| prec | 53 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|
| us/op | 8 367 | 12 100 | 27 430 | 51 300 | 125 100 |

**~3.3× from 256→53; ~2× per halving.** Adaptive precision (compute low, certify
the ball, bump only on straddle) is a 2–3× lever where applicable. Sources
([Johansson, *Using ball arithmetic*](https://flintlib.org/doc/using.html)):
the guess-and-verify precision-doubling loop keyed on `acb_rel_accuracy_bits` is
the canonical arb idiom.

### 3d. FLINT intra-op threading — capped at ≤1.9× at EVERY scale

prec=256, `flint_set_num_threads(t)`:

| op | 1 thr | best | speedup | saturates |
|---|---|---|---|---|
| certified eig / opnorm n=32 | 210 ms | 190 ms | **1.1×** | — (does not parallelize) |
| sgn n=32 | 606 ms | 472 ms | 1.4× | 4–8 |
| acb_mat_mul n=32 | 3.0 ms | 1.6 ms | 1.9× | 4–8 |
| acb_mat_mul n=128 | 34.8 ms | 20.1 ms | 1.73× | 16 |
| acb_mat_mul n=256 (superop scale) | 145 ms | 80 ms | **1.81×** | 16 |

**Verdict (corroborated by [Johansson 2019, arXiv:1901.04289 §VI](https://arxiv.org/abs/1901.04289)
and [Parallelizing Arb](https://fredrikj.net/blog/2022/05/parallelizing-arb-for-fun-and-profit/)):**
FLINT's matrix code is single-threaded; multithreaded mat-mul is "typically not
useful" at these sizes; threading wins only appear at ≥10⁴ digits. **Stop chasing
`flint_set_num_threads`.** Parallelize across *independent* problems at the aic
level (§4) with the inner libs pinned to 1 thread.

### 3e. Compiler flags on the glue — 0–7% (noise)

`-O3 -march=native -mtune=native -funroll-loops -flto` vs `-O2`, re-running the
arb probe: mul/eig/sgn/precision all within ±7% (noise). **Confirmed: the time is
inside precompiled libflint; recompiling the aic glue cannot touch it.** A local
FLINT/GMP rebuild with `-march=native` is a separate (heavier) candidate — and
per [FLINT building docs](https://flintlib.org/doc/building.html) FLINT does
*compile-time* ISA selection and `--enable-avx2` only accelerates the small-prime
FFT (a large-integer path, NOT our n≤32/prec≤256 critical path), so expect little.

### 3f. FLINT allocation churn — 3.5–7% per op (minor)

Microbench: `acb_mat_init`+`acb_mat_mul`+`acb_mat_clear` vs reusing the matrix
(prec=256): churn overhead n=8 **7.2%**, n=16 **3.5%**, n=32 **4.9%**; bare
init/clear is 0.5–10 µs. A real but secondary lever — hoist temporaries / use
`_acb_vec_init` + `TMP_INIT` in hot loops ([arXiv:1901.04289 §V.D](https://arxiv.org/abs/1901.04289)).
DHAT (§7) will give per-site authoritative numbers.

### 3g. OpenBLAS audition (double path) — NOT a clear win (clean run)

`OPENBLAS_NUM_THREADS=1 LD_PRELOAD=libopenblas.so.0` (deterministic config per
R4) vs reference netlib, `zheev`:

| n | reference | OpenBLAS-1thr | result |
|---|---|---|---|
| 4 | 22.2 µs | 4.6 µs | OpenBLAS **4.8× faster** |
| 8 | 18.8 µs | 15.6 µs | 1.2× faster |
| 16 | 83.3 µs | 69.9 µs | 1.19× faster |
| 32 | 387 µs | 508 µs | OpenBLAS **1.3× SLOWER** |

**Counterintuitive but explained** ([BLASFEO arXiv:1704.02457](https://arxiv.org/abs/1704.02457),
[OpenBLAS #1614](https://github.com/xianyi/OpenBLAS/issues/1614)): at small n,
`zheev` is tridiagonalization-bound (Level-2/memory-bound), *not* GEMM-bound, so
OpenBLAS's GEMM advantage barely applies and its LAPACK loses to netlib at n=32.
The big n=4 win is real but on already-microsecond absolute times. **OpenBLAS does
NOT Pareto-dominate** for this workload; the double path is also not a bottleneck
(it's 58–149× faster than arb). Decision: keep reference as the deterministic
default (it also feeds the double-vs-arb cross-check, which needs run-to-run
determinism — multithreaded OpenBLAS has no reproducibility guarantee,
[OpenBLAS #5267](https://github.com/OpenMathLib/OpenBLAS/issues/5267)). Use
OpenBLAS-1thr *selectively* only if an n≤8 double hot spot is found.

### 3h. Callgrind: where the instructions actually go (capstone)

One `aic_sgn` + one `aic_mat_eig_hermitian` + one `aic_mat_opnorm` (n=16, prec=256)
under callgrind — self-cost by function:

| function | lib | Ir % |
|---|---|---|
| `_arb_dot_addmul_generic` | flint | 17.1% |
| `acb_approx_dot` | flint | 16.3% |
| `__gmpn_mul_basecase` | gmp | 15.8% |
| `__gmpn_rshift` | gmp | 9.9% |
| `acb_dot` | flint | 9.6% |
| `_arf_set_round_mpn` | flint | 5.6% |
| `__gmpn_mul_n` / `add_n` / `sub_n` | gmp | 8.0% |
| `_arf_add_mpn` | flint | 2.3% |

- **~43% is the arb fused dot product** (`acb_dot`/`acb_approx_dot`/
  `_arb_dot_addmul_generic`) — matmul, eig, and sgn all bottom out here.
- **~34% is raw GMP limb arithmetic** (mantissa mul/shift/add at prec=256 ≈ 4 limbs).
- **~8% is arf rounding.**
- **No aic-glue function appears in the top 11** — it is *all* libflint/libgmp.

**Implication**: there is no micro-optimizable hotspot in our code. The only ways
to move the needle are (i) **fewer** dot products (fewer iterations → scaled-NS,
L2.3), (ii) **smaller** ones (lower precision → L3), (iii) **parallel** ones
(task-level → L1). The 34% GMP share is the one place a `-march=native` GMP
rebuild (Haswell `mulx`/`adcx`/`adox` bignum kernels) *could* help — but Debian's
libgmp uses runtime CPU dispatch and is likely already using them; audition before
investing (§3e caveat). DHAT cross-check: 52,661 allocations / 4.5 MB for those 3
ops — many but cheap blocks (the §3f 3.5–7% time cost), not a primary lever.

## 4. Where the wins actually are (lever ranking)

Ranked by (measured or research-backed) leverage × feasibility. Estimated wins
from the codebase sweep (R1/R2) are labeled **[est]** and must be measured before
trusting.

### L1 — Task-level parallelism across independent work (the biggest arb lever)

The bottleneck primitive (certified eig) won't thread internally (§3d), **but it
is invoked per-block / per-ampliation / per-sweep-point independently**, on an
otherwise idle 18-core box. Pattern ([scikit-learn parallelism](https://scikit-learn.org/stable/computing/parallelism.html),
[FLINT threading](https://flintlib.org/doc/threading.html)): `#pragma omp
parallel for schedule(dynamic)` over tasks, **inner libs pinned to 1 thread**
(`flint_set_num_threads(1)`, `openblas_set_num_threads(1)`, `OMP_MAX_ACTIVE_LEVELS=1`)
to avoid P×B×F oversubscription; per-thread `flint_cleanup()`, one
`flint_cleanup_master()` at exit. FLINT is thread-safe in pthreads builds with
thread-local caches. There is **no CPU batched Hermitian eig** — task-parallel
loop is the route. Identified independent loops (R1/R2):

| site | independent unit | iters | win [est] |
|---|---|---|---|
| `aic_opspace_cert.c:78` ampliation levels | level n∈{1,2,4,…,n_B} | ~log n_B (4–6) | 4–6× the HOPM stretch [est] |
| `aic_opspace_ampliate.c` HOPM restarts | independent restarts | nrest 4–8 | 3–4× [est] |
| `aic_factorize_upsilon3.c:139` block loop | per-block defect | num_blocks 2–10 | 3–6× [est] |
| `aic_cstar_stages.c` Stage-2 classes | per equivalence class | num_classes 2–10 | 4–8× [est] |
| `aic_cstar_build.c:150` `find_dim_gt1` scan | per-projection dim_S | O(dim_A) | parallel scan [est] |
| `aic_assoc_superop.c:48` superop build | per-column Φ(E_pq) | n² columns | n²-way [est] |

### L2 — Arb-algorithm auditions (faster certified primitives)

From [Johansson 2019](https://arxiv.org/abs/1901.04289) (measured tables) — these
are the highest-ceiling arb-path changes; each is a Law-4 audition:

1. **certified eig: vdHM O(n³) default, Rump O(n⁴) only as cluster fallback.**
   `acb_mat_eig_simple` (vdHM) is 1.25–1.6× over bare QR; `acb_mat_eig_multiple`
   /Rump is 2.3×(n=10)→6×(n=100) slower. *Status in code*: `aic_mat_eig_hermitian`
   already uses vdHM `acb_mat_eig_simple`; `aic_mat_eig_hermitian_multiple` uses
   Rump for the **degenerate block-algebra spectra** (where vdHM returns 0). Verify
   `acb_mat_eig_multiple` already does simple-first (it claims to); if a hot path
   forces Rump where a cluster split isn't needed, that's a 2–6× win.
2. **fast-guess + certify-once.** Run iterations (sgn Newton–Schulz, solves) in
   `acb_mat_approx_*` (drops radius, ~2× faster at p≤128, avoids per-iteration
   radius blow-up), then ONE ball residual to certify. The sgn loop currently runs
   fully in ball arithmetic. **[est ~2× on sgn + tighter balls].**
3. **scaled Newton–Schulz for sgn** (Chen–Chow, [PDF](https://faculty.cc.gatech.edu/~echow/pubs/chen-chow-2014.pdf)):
   `X_{k+1}=½α_k X_k(3I−α_k²X_k²)`, α from a scalar recurrence (one initial
   smallest-eigenvalue estimate, no per-step eigensolve) → **~2× fewer iterations**
   near a small gap. **MUST use the §4 stabilized α** (raw optimal is numerically
   unstable — a Rule-3 trap). This is bead **`aic-68c`** (Kenney–Laub), now
   confirmed as a real lever.
4. **ρ/opnorm need only a spectral bound, not full eig.** `aic_mat_herm_max_eig`
   already uses `acb_mat_eig_global_enclosure` (Miyajima, one O(n³) pass) — good —
   but still pays for a full `acb_mat_approx_eig_qr`. For *just* λ_max a verified
   power iteration + one O(n²) ball residual (Weyl/Kato–Temple) is cheaper. Audition
   for the `ρ(Φ²−Φ)<1/4` Gelfand guard and opnorm.
5. **certified solves: preconditioned (Hansen–Smith), never ball-LU** (ball-LU
   loses O(n) digits). `acb_mat_solve` already auto-selects; verify the factorize
   `M1inv` path uses it.

### L3 — Adaptive precision (2–3×, §3c)

Wrap drivers in guess-and-verify precision-doubling keyed on
`acb_rel_accuracy_bits` with a hard prec cap that fails loud (Rule 4). Candidate
sites (R1): `aic_errreduce` (prec=256 fixed, called O(dim_A)×), `aic_opspace_cert`
HOPM, the assoc Gelfand/sgn on the superoperator. **[est 25–45%]** — measure.

### L4 — FLINT intra-op threads (1.1–1.9×, §3d) — free but weak; mind oversubscription with L1.

### L5 — Alloc churn (3.5–7%, §3f), compiler flags (0–7%, §3e), local FLINT rebuild — minor.

## 5. Build / backend recommendations

- **BLAS**: keep **reference netlib** as the deterministic default (§3g; feeds the
  cross-check oracle). OpenBLAS-1thr only for a measured n≤8 double hot spot. MKL
  only if a contractual reproducibility guarantee is needed (`MKL_CBWR=AVX2,STRICT`).
- **NEVER `-ffast-math`/`-Ofast` on any TU in the arb/ball path** — it enables
  `-fno-rounding-math`, letting the compiler fold FP under round-to-nearest and
  silently invalidate the certified enclosure ([Simon Byrne](https://simonbyrne.github.io/notes/fastmath/),
  [CGAL #3180](https://github.com/CGAL/cgal/issues/3180)). Hard rule; set flags
  per-target. (Currently moot — build is `-O2` — but a guardrail for any future
  "speed up the glue" change.)
- **Glue**: `-O2 -march=native -flto` is low-risk but low-yield (§3e). PGO is
  optional, low end of 5–20% (hot math is precompiled).
- **Threading**: pin inner libs to 1 thread; parallelize at the task level (L1).

## 6. Audition slate (Law 4 — keep candidates alive until a benchmark retires one)

| lever | candidate A | candidate B | metric to decide |
|---|---|---|---|
| sgn | Newton–Schulz (current) | scaled-NS stabilized / DWH-Halley | iters AND certified-ball width (Pareto) |
| certified eig | vdHM-first/Rump-fallback | Rump always | wall × n, on degenerate block spectra |
| λ_max / ρ | global_enclosure (current) | verified power iteration O(n²) | wall, ball tightness on the guard |
| iteration certification | full ball arithmetic | approx + certify-once | wall × ball width |
| precision | fixed 256 | adaptive doubling | wall × cert success rate |
| parallelism | serial (current) | OpenMP task-level | scaling vs core count, determinism |

## 7. Profiling playbook (the deep-attribution workflow)

`perf` is unavailable under WSL2 (no hardware PMU). Valgrind is the instrument
(pure user-space, deterministic, sees into libflint with `-g`). Run on a
*reduced* corpus (10–50× slowdown) and trust relative attribution.

1. **Leak-clean baseline** first (`flint_cleanup_master()` at exit) so alloc
   profiles aren't teardown noise: `valgrind --leak-check=full ./driver`.
2. **Callgrind** — where the time goes, *inside* libflint:
   `valgrind --tool=callgrind --dump-instr=yes --collect-jumps=yes --cache-sim=yes ./driver`
   then `callgrind_annotate --inclusive=yes --tree=both`. (KCachegrind for the graph.)
3. **DHAT / heaptrack** — allocation-site churn (the §3f story at per-site
   resolution): `valgrind --tool=dhat ./driver`.
4. **Cachegrind** — confirm the n²×n² superoperator working set vs the 45 MiB L3.
5. **Massif** — peak heap over a dim-sweep (the n²×n² superop is the memory
   hotspot, ~5 MB per acb 256×256 at prec=256).

(A first callgrind+DHAT pass on a 1-sgn/1-eig/1-opnorm driver is in
`build-perf/cg_top.txt` / `dhat.out` — see §8 once it lands.)

## 8. Status / open beads

- Epic **`aic-xv4`** (this campaign). Children to file: per-lever beads (L1
  OpenMP task-parallel opspace/factorize/cstar; L2 sgn scaled-NS audition =
  `aic-68c`; L2 eig vdHM/Rump verify; L2 approx+certify-once; L3 adaptive
  precision; the shared bench corpus `aic-f9u.1` + Pareto metrics `aic-f9u.2`).
- **`aic-nsb`** — the `test_opspace_o2` full-suite failure (P1, surfaced here).
- Existing relevant: `aic-68c` (Kenney–Laub sgn), `aic-f9u.1/.2` (bench corpus +
  Pareto reporter), `aic-erz`/`aic-rcm` (cut canary wall cost).

> The headline: this is an **arb-path, certified-numerics** performance problem,
> not a BLAS-or-threads problem. The wins are task-level parallelism over
> independent work (L1), faster certified primitives (L2), and adaptive precision
> (L3) — in that order. The "obvious" levers (OpenBLAS, `flint_set_num_threads`,
> `-march=native`) are measured to be marginal or counterproductive here.
