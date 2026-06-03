# Architecture, certified arithmetic, and the constructive mandate

A citeable brief for the documentation site. Three explanation pages in one file:
the three-layer architecture (A), the certified-arithmetic story (B–C), and the
constructive mandate (D), plus a practical build/consume section (E). Every claim
is anchored to a source file and, where the math is load-bearing, to a
`paper/src/approximate_algebras.tex` line. Numbers are measured values from the
test corpus, not estimates.

Ground truth for the math is the arXiv source
[`arXiv:2405.02434`](https://arxiv.org/abs/2405.02434), Alexei Kitaev,
*Almost-idempotent quantum channels and approximate C\*-algebras* (v1).

---

## A. Three-layer architecture

The package mirrors the `../su2-fft` sibling: three layers, cross-validated
against each other (`CLAUDE.md:108-124`).

```
                 ┌──────────────────────────────────────────────────────┐
                 │  AlmostIdempotentChannels.jl                          │
                 │    high-level API, CertifiedBracket result type,      │
                 │    the Watrous diamond-norm SDP (MOSEK extension),    │
                 │    runtests.jl cross-checks                           │
                 └───────────────────────┬──────────────────────────────┘
                                         │  @ccall, flat-double ABI only
                                         │  (int dims + Ptr{Cdouble} arrays;
                                         │   NO acb_mat_t / slong cross the line)
                                         ▼
                 ┌──────────────────────────────────────────────────────┐
                 │  flat-double shims  (src/aic_*_shim.c, C2–C5)         │
                 │    pure marshalling: double arrays <-> acb_mat        │
                 └───────────────────────┬──────────────────────────────┘
                                         │  internal acb_mat / aic_ucp_kraus
                                         ▼
   ┌─────────────────────────────────────────────────────────────────────────┐
   │  C cores  (src/*.c, include/aic/*.h)   — two number paths                │
   │                                                                          │
   │    double path  (BLAS / LAPACK: zheev, zgesvd, zgees / LAPACKE)          │
   │        fast; the prec=53 ANCHOR; degenerate-spectrum eigensolvers        │
   │                                                                          │
   │    arb path     (FLINT acb_mat / arb, explicit prec)                     │
   │        CERTIFIED error balls; proves the paper's O(ε) bounds             │
   └─────────────────────────────────────────────────────────────────────────┘
```

### Why both number paths

`CLAUDE.md:121-124`, verbatim:

> Every nontrivial routine exists (eventually) in **both** number paths, and the
> arb path at `prec=53` must agree with the double path — that is a primary
> cross-check. The arb path is what *certifies* the paper's $O(\eps)$ bound; the
> double path is what makes the test suite fast.

The double path is the fast `prec=53` anchor and carries the degeneracy-robust
LAPACK eigensolvers (the certified simple-eig arb path aborts on the generically
degenerate Choi spectrum — `aic_ucp.h:44-51`). The arb path is the point of the
project: it returns rigorous error balls so the package can certify Kitaev's
$O(\eps)$ inequalities rather than trust `double` (`CLAUDE.md:34-36`, Law 15
`CLAUDE.md:266`).

### The ABI: flat-double `@ccall` shims (C2–C5)

Julia's `ccall` passes `int` / `double` / `Ptr{Cdouble}` cleanly but **not**
FLINT structs — `acb_mat_t` is an opaque array type and `slong` is FLINT's
machine-word int. So the public language boundary uses `int` dimensions and
`Ptr{Cdouble}` arrays only; the shim converts to/from the internal
`aic_ucp_kraus` / `acb_mat` at the boundary, runs the *same* arb cores the C
tests exercise, and writes the result back into caller-owned flat double arrays
(`aic_ucp_shim.h:1-23`, `ALGORITHM.md:601-639`). It is pure marshalling — no new
math.

The shim symbols dlsym'd into Julia (`libaic.jl:67-84`):
`aic_ucp_choi_diff_d`, `aic_cbnorm_eigfree_d`, `aic_cbnorm_certify_d`,
`aic_assoc_summary_d`, `aic_main_iso_summary_d`,
`aic_factorize_artifacts_sizes_d`, `aic_factorize_artifacts_d`.

### Storage-layout contract — Kraus, Heisenberg picture

This is load-bearing and pinned in the header (`aic_ucp.h:6-21`). A UCP map
$\Phi : B(K) \to B(H)$ is stored as $r$ Kraus operators $K_a : H \to K$, each a
`dim_K × dim_H` matrix. The action is the **observable (Heisenberg)** form

```
Φ(X) = Σ_a K_a† X K_a            (X in B(K)),         [.tex:1570]
```

with unitality $\Phi(1_K) = \sum_a K_a^\dagger K_a = 1_H$ (`.tex:1571`). This is
*not* the Schrödinger/state dual $T = \Phi^*(X) = \sum_a K_a X K_a^\dagger$ —
getting the dual backwards silently transposes everything (`aic_ucp.h:18-21`,
`CLAUDE.md:315-318`). The flat ABI layout is row-major:
`kraus_re[a*n*n + i*n + j]`, `r` ops each `n × n` (`aic_ucp_shim.h:17-22`).

### Convention-A Choi (paper / Watrous / QETLAB / Qiskit)

With ONB $\{e_i\}$ of $K$ and $E_{ij} = |e_i\rangle\langle e_j|$,
$C_\Phi = \sum_{ij} E_{ij} \otimes \Phi(E_{ij})$ in $B(K) \otimes B(H)$, $K$
factor LEFT/MAJOR (`.tex:1575`, shard G). Hand-deriving $\Phi(E_{ij})[a,b]$ from
the Heisenberg action gives the pinned per-entry formula (`aic_ucp.h:30-39`,
verbatim):

```
C_Phi[i*dim_H + a, j*dim_H + b] = sum_x conj(K_x[i,a]) * K_x[j,b].
```

The conjugation is on the **first** $(i,a)$ factor. This conjugation is
load-bearing: it is exercised only by the complex asymmetric channel
(`complex_qubit`), and dropping it makes the C Choi diverge from the Julia
reference by `1.28` in operator norm — far above the `1e-11` fixture tolerance
(`ALGORITHM.md:264-266`, `aic_ucp.h:104-107`). $\Phi$ is CP iff $C_\Phi \succeq 0$;
$\Phi$ is unital iff $\mathrm{Tr}_K(C_\Phi) = 1_H$ (trace out the LEFT/$K$
factor).

### Native-library lifecycle

The Julia package finds `libaic.so` via **Preferences**: the default is
`build/libaic.so` relative to the repo root, overridable by
`set_libaic_path!(p)` which writes a `LocalPreferences.toml` entry read at
precompile time (`libaic.jl:7-12,27-58`). On load, `__init__`
(`libaic.jl:104-135`):

1. **pre-loads the system `libgmp.so.10`** (`RTLD_NOW | RTLD_GLOBAL`, full path
   then bare-name fallback, in try/catch) **before** libaic/libflint are mapped;
2. `dlopen(LIBAIC_PATH, RTLD_NOW | RTLD_GLOBAL | RTLD_DEEPBIND)`;
3. `dlsym`s each shim symbol into a module-private `Ref{Ptr{Cvoid}}`.

The libgmp preload is the `../su2-fft` fix (bead su2fft-e5z) for the Julia
bundled-libgmp ABI clash with FLINT: `RTLD_DEEPBIND` on libaic alone is
insufficient because libflint resolves `__gmpn_*` against the already-loaded
Julia-bundled libgmp; pre-loading the system libgmp keeps its definitions first
in the global symbol table (`libaic.jl:111-127`, `ALGORITHM.md:699-708`).

---

## B. Certified arithmetic — why arb, not float

### What a `CertifiedBracket` is

Every scalar result is a `CertifiedBracket` `lo ≤ x ≤ hi` — a **proof of an
inequality, not a float guess** (`types_core.jl:11-37`, design §2.5). The C arb
error balls are rounded **outward** (FLOOR on the lower endpoint, CEIL on the
upper) so the bound is still rigorous after the `double` conversion:

```
*lo = floor(arb_lower(lo_ball))   (arb_get_lbound_arf -> arf_get_d, ARF_RND_FLOOR)
*hi = ceil (arb_upper(hi_ball))   (arb_get_ubound_arf -> arf_get_d, ARF_RND_CEIL)
```

These are the only two rounding directions that keep `*lo ≤ x ≤ *hi` after the
double conversion (`aic_ucp_shim.h:46-52`, `ALGORITHM.md:627-633`). The Julia
inner constructor validates `lo ≤ hi`, and when an SDP point estimate `value` is
attached, that it lies in `[lo − 1e-9, hi + 1e-9]` (the slack absorbs the benign
double-rounding of an SDP estimate against an outward-rounded ball,
`types_core.jl:22-36`).

### Where precision is load-bearing

`double` is fine almost everywhere, but the bound itself fails silently in three
regimes (`CLAUDE.md:327-331`):

1. $\eps$ is tiny and the bound is a **difference of $O(1)$ quantities** — the
   associator $(XY)Z - X(YZ)$, or $\Phi^2 - \Phi$. Catastrophic cancellation
   eats the answer in `double`; the arb ball tracks the true uncertainty.
2. A routine **inverts a near-singular operator** —
   $(\La_{|X|} + \Ra_{|X|})^{-1}$ at `.tex:615`.
3. A **spectral gap comparable to $\eps$** — the verdict (rank, PSD, basin) would
   be a guess in `double`; the arb path fails loud instead (Rule 4,
   `CLAUDE.md:191-196`).

This is also why the constructions prefer Hermitian/normal elements and certify
everything else: for non-normal $X$ the spectrum is fragile (the paper's own
$3\times3$ example has eigenvalues $\sim t^{1/3}$, `.tex:540-544`,
`CLAUDE.md:292-297`).

### The eig-free Watrous bracket — rigorous by an inequality chain, loose by design

The central quantity is the $\eta$-idempotence defect (`.tex:354`, verbatim):

> "A UCP map $\Phi: B(H)\to B(H)$ is called $\eta$-idempotent if
> $\|\Phi^2-\Phi\|_\cb \le \eta$."

$\Lambda = \Phi^2 - \Phi$ is Hermiticity-preserving but **not** CP, so the CP
closed form $\|\Phi\|_\cb = \|\Phi(1)\|_\op$ does not apply; the cb-norm equals
the diamond norm and is "efficiently computable using semidefinite programming"
(`.tex:347-352`). For a Hermiticity-preserving map the ampliation truncation
$N = \dim(\text{input}) = n$ is rigorous (Watrous 2012,
`aic_cbnorm.h:13-16`).

The cheapest, **always-valid, never-aborting** rung uses *only* the certified
Frobenius norm of $J = \mathrm{Choi}(\Lambda)$ (the $n^2 \times n^2$
Convention-A Choi) — no SDP solver, no eigendecomposition, so it dodges the
degenerate-eig wall (`aic_cbnorm.h:17-22`,
`aic_cbnorm_eigfree_ball_choi`). The exact inequality chain
(`ALGORITHM.md:307-326`):

- **Upper.** From $P+Q = I_{n^2}$ with $P,Q \succeq 0$, Schur-factoring the
  block-PSD constraint $X = P^{1/2} K Q^{1/2}$ with $\|K\|_\op \le 1$ gives
  $\|X\|_F \le \sqrt{\mathrm{tr}\,Q} \le n$, so
  $\mathrm{optval} \le \|J\|_F\,\|X\|_F \le n\|J\|_F$, hence
  $\eta \le 2\|J\|_F$.
- **Lower.** $(\mathrm{id}\otimes\Lambda)(\rho) = J/n$ for the normalized
  maximally-entangled $\rho$, and
  $\eta \ge \|J/n\|_1 \ge \|J\|_2/n = \|J\|_F/n$ (Schatten
  $\|\cdot\|_1 \ge \|\cdot\|_2$).

So $\eta \in [\,\|J\|_F/n,\ 2\|J\|_F\,]$, ratio $\mathrm{hi}/\mathrm{lo} = 2n$ —
rigorous but loose by $\sim 2n$ **by design**. It is built to *certify* (contain
the truth, never abort), not to *compete* on tightness; it is the fallback the
tight certifier dispatches to near $\eta = 0$ (`aic_cbnorm.h:17-22`,
`ALGORITHM.md:324-326`). For a general rectangular Hermiticity-preserving map
$f: M_{d_{in}} \to M_{d_{out}}$ the same derivation generalizes to
$\|f\|_\diamond \in [\,\|J\|_F/d_{in},\ d_{out}\|J\|_F\,]$ — the rectangular core
needed for the $\|v\|_\cb / \|v^{-1}\|_\cb$ isomorphism brackets, since the
self-map $P+Q=I_{n^2}$ primal does not hold for a rectangular map
(`aic_cbnorm.h:53-72`, FINDINGS §C12, §C22(a)).

### The tight rung — arb fed MOSEK feasible points

The tight certified ball collapses the $2n$ width to roughly (MOSEK duality gap +
arb radius). Two **independent** MOSEK feasible points bracket $\eta$ by weak
duality (`aic_cbnorm.h:92-113`, `ALGORITHM.md:481-530`):

- MAX-primal $(X,P,Q)$ → rigorous **lower** bound $\eta \ge (2/n)\,\mathrm{Re}\,\mathrm{tr}(J^\dagger X)$;
- MIN-dual $(Y_0,Y_1)$ → rigorous **upper** bound
  $\eta \le \tfrac{1}{2}(\|\mathrm{Tr}_2 Y_0\|_\infty + \|\mathrm{Tr}_2 Y_1\|_\infty)$.

Each MOSEK point is only *approximately* feasible, so each is **restored to exact
feasibility in arb** — a convex combination toward the Slater center for the
primal, an eigenvalue shift for the dual — before its objective is read off; only
then is the bound rigorous. Every PSD test is the one-sided
`herm_max_eig(−M)` global enclosure (degeneracy-robust, no full eig, dodging the
aic-w4o.1 wall). Dispatch falls back to the eig-free bracket near $\eta = 0$ or
when the MOSEK points cannot be restored; it fails loud only if even the fallback
straddles (Rule 4).

Measured tight widths at `prec=128` are `~1e-13`/`1e-14`, about $10^{12}\times$
tighter than the eig-free bracket (`README.md:216-228`):

```
fixture        |  eig-free width  |  tight width
---------------+------------------+--------------
phit(2, 0.3)   |  5.46e-01        |  5.76e-13
phit(2, 0.1)   |  2.34e-01        |  1.60e-13
paper(0.1)     |  2.01e-01        |  4.60e-14
```

The floor is the double MOSEK feasible point, not the arb prec: `prec=64` and
`prec=256` give the same width on `phi_t4_0p2` (`ALGORITHM.md:548-551`).

---

## C. The cross-check ladder — how to trust the output

Cross-checks beat unit tests here: unit tests catch typos, cross-checks catch
algorithmic errors, and the cross-check is added **before** the code (Rule 6,
`CLAUDE.md:204-213`). The ladder, weakest to strongest (`README.md:178-195`,
`CLAUDE.md:334-339`):

**Rung 1 — double vs arb at prec = 53.** The same routine on both number paths
agrees to `~1e-10`: the C/arb Choi of $\Phi^2-\Phi$ matches the pure-Julia one
(`test_core.jl`). *Catches:* a marshalling or convention slip in either path that
the other path does not share — e.g. a transposed reshape or a dropped
conjugation. This is the `su2-fft` anchor.

**Rung 2 — certified containment of the closed form.** The bracket provably
contains the closed-form $\eta$ for the $\Phi_t = (1-t)\,\mathrm{id} + t\,\mathrm{Dep}$
family at every $t$, with $\eta = t(1-t)\cdot 2(1-1/d^2)$ (`test_core.jl`;
generator anchors `ALGORITHM.md:240-249`). *Catches:* a wrong dimension power or
a sign error in the bound that still "runs" — the bracket then misses
$\eta_\mathrm{ref}$. The $2(1-1/d^2)$ dimension canary (1.5 at $d=2$, 16/9 at
$d=3$, 15/8 at $d=4$) pins the $2/n$ normalization at every dimension in the
corpus.

**Rung 3 — the $\eta = 0$ exact-idempotent oracle.** When $\Phi$ is exactly a
conditional expectation, every defect collapses below `1e-9` and the extracted
block structure is exactly right: `identity(2) → [2]`,
`block_cond_exp(4,2) → [2,2]`, `dephasing(3) → [1,1,1]` (`README.md:165-169`,
`test_factorize.jl`). This is the cleanest ground truth in the paper — zero
defect is unambiguous. Measured maxima: `4.4e-75` for $\|\Delta\Upsilon - \Phi\|_\cb$
and `3.9e-75` for $\|\Upsilon\Delta - P_B\|_\cb$ (`README.md:190-192`).
*Catches:* the whole class of "test that cannot fail" — a construction that looks
right but never reduces to the exact Choi–Effros decomposition when $\eta=0$.

**Rung 4 — MOSEK strong duality.** The Watrous primal and dual diamond-norm
values agree with each other and with the analytic anchor to `1e-6`; measured max
primal–dual gap `1.23e-11` (`README.md:193-195`, `test_sdp.jl`). *Catches:* a
wrong SDP normalization — the gap closes only when the $2/n$ (primal) and $/2$
(dual) calibrations are both right; it pins the dual normalization the eig-free
rung cannot see.

The C side runs an identical ladder per module (`ALGORITHM.md:164-178`,
`272-283`), adding the `prec=53` vs `prec=256` precision self-consistency check
and source-level **mutation proofs** (perturb the impl, confirm RED, restore) so
each constant in a bound is pinned by at least one family (Rule 7,
`CLAUDE.md:215-223`).

---

## D. The constructive mandate

The project's reason to exist (`CLAUDE.md:39-79`, THE MANDATE, verbatim):

> **Many proofs in the paper are non-constructive and phrased for possibly
> infinite-dimensional spaces. But every *theorem* applies to finite-dimensional
> matrices, and in finite dimensions the objects those proofs merely assert to
> exist are computable. Our goal is a constructive algorithm for *every* result.**

For each result two things are kept separate: the paper's proof *technique* (and
whether it is constructive) and the finite-dim algorithm we *implement* — which
need not follow the proof (Law 3, `CLAUDE.md:148-151`). The canonical moves
(`CLAUDE.md:58-72`):

| Paper's non-constructive tool | Constructive finite-dim route |
|---|---|
| **Nontrivial projection** via Lefschetz–Hopf fixed-point + H-space homotopy (§6) | Spectral split of a Hermitian element; or $\theta(2P-I)$ on a near-idempotent ($\theta(X)=(I+\mathrm{sgn}\,X)/2$); or minimize $\|P^2-P\|$ over Hermitian $P$ bounded off $\pm I$ (FINDINGS §B1) |
| **Holomorphic functional calculus** via contour integrals (§3) | Eigendecomposition; or a Newton / Denman–Beavers iteration for $\mathrm{sgn}(X)$; or a Zolotarev rational approximation with a certified arb error |
| **Haar-measure diagonal** $D=\int dU\,(U^\dagger\otimes U)$ (§2/§9) | An explicit finite average over the regular representation — for the Pauli case a closed form (FINDINGS §B2) |
| **Implicit / inverse function theorem** (§3, `lem_invfun`) | *Already* constructive: the Banach contraction iteration $x_n = g_y(x_{n-1})$ **is** the algorithm, with the paper's $c<1$ giving the convergence rate |

The correctness specification is the paper's bound itself (`CLAUDE.md:74-76`,
verbatim):

> **The paper's $O(\eps)$ bound is the algorithm's correctness specification.**
> The constructive routine must produce an output meeting the bound the theorem
> promises; that bound is the assertion the cross-check verifies.

The isomorphism of the main theorem is $O(\eps)$, never exact, and the constant
is **universal and dimension-independent** (`th_main`, `.tex:460`). Tests assert
the bound and that the constant does *not* grow with $n$ — the naive
Haar/cohomology route is rejected precisely because its error grows $\propto n$
(`.tex:484`, `CLAUDE.md:309-313`). The measured constant $c = \text{isodefect}/\eta$
stays $\le 0.047 \ll 1$ with OLS slope $-2.7\mathrm{e}{-4}$ across $\dim A \in
\{8,12,16,18,20\}$ (`README.md:201-213`).

If a result resists constructivization, that is a **stop condition** to escalate
with the specific obstruction, not a corner to paper over (`CLAUDE.md:78-79`,
`454-469`).

---

## E. Build and consume

### Dependencies (`BUILDING.md:8-22`)

| Component | Minimum | Debian/Ubuntu package |
|---|---|---|
| C11 compiler | — | `gcc` / `clang` |
| CMake | 3.24 | `cmake` |
| FLINT (incl. arb/acb) | 3.0 | `libflint-dev` |
| LAPACKE (C interface) | — | `liblapacke-dev` |
| LAPACK (Fortran) | — | `liblapack-dev` |
| BLAS | — | `libblas-dev` |

```sh
sudo apt install build-essential cmake libflint-dev liblapacke-dev liblapack-dev libblas-dev
```

arb is *inside* FLINT 3.x — there is no separate `libarb`.

### Build and test (`BUILDING.md:24-73`, `CLAUDE.md:385-400`)

```sh
cmake -S . -B build            # CMake is canonical; `make` forwards to it
cmake --build build            # -> build/libaic.so (SONAME libaic.so.0) + build/libaic.a

ctest --test-dir build -L fast # sub-second laptop gate (10 drivers, ~0.5 s)
ctest --test-dir build         # full suite (fast + slow: arb / MOSEK-fixture / large-n)
ctest --test-dir build -R test_funcalc   # one test by name
```

The ~85 `src/*.c` compile **exactly once** into an `OBJECT` library (`aic_objs`)
from which both the shared and static libraries are assembled — the sources are
never compiled twice. Every test carries a `TIMEOUT` (fail-loud, no hangs); the
`fast`/`slow` labels split the inner-loop gate from the heavy drivers. `assert()`
stays live in **every** build type — the build strips `-DNDEBUG` from the
optimized configs, since the library uses `assert()` for its 342 fail-loud Rule-4
invariants (`BUILDING.md:149-153`).

### Install and consume (`BUILDING.md:93-138`)

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=<prefix>   # set the prefix at CONFIGURE time
cmake --build build
cmake --install build
```

`CMAKE_INSTALL_PREFIX` must be set at configure time: the generated `aic.pc` and
`AICConfig.cmake` bake the prefix at configure, so a late `--install --prefix`
can leave them pointing at the wrong location.

```cmake
find_package(AIC CONFIG REQUIRED)
target_link_libraries(app PRIVATE AIC::aic)          # or AIC::aic_static
```
```sh
cc app.c $(pkg-config --cflags --libs aic) -o app            # shared
cc app.c $(pkg-config --cflags --libs --static aic) -o app   # static: --static is required
```

Two Debian gotchas, both handled by the build (`BUILDING.md:141-148`): Debian's
`flint.pc` omits `-lflint` (`Libs:` is only `-lgmp -lmpfr`), so `aic.pc` lists
`-lflint` itself in `Libs.private` and the exported CMake targets bake the
absolute `libflint.so` path; and `find_package(LAPACK)` provides Fortran LAPACK
only, not LAPACKE — LAPACKE is found separately via `pkg_check_modules`.

### Julia: finding `libaic.so`, and the solver stance

The Julia package finds the native library through **Preferences**: the default
is `build/libaic.so` relative to the repo root, so building once with CMake needs
no `Pkg.build`; to point elsewhere call
`AlmostIdempotentChannels.set_libaic_path!("/abs/path/libaic.so")` (writes
`LocalPreferences.toml`; restart Julia to apply) (`libaic.jl:7-12,40-58`).

```sh
julia --project=julia/AlmostIdempotentChannels.jl -e 'using Pkg; Pkg.test()'
```

The package is **solver-free by default**: the eig-free certified bracket and the
flat-double ccall surface need no SDP solver. Only the Watrous diamond-norm SDP
(`idempotency_defect`, the tight bracket, the `diamond_*` entry points) requires
MOSEK, and that lives **only** in the optional package extension `AICMosekExt`
(MOSEK is a `[weakdeps]`). Without the extension loaded those functions throw a
helpful error pointing at the install hint, registered as a
`Base.Experimental` error hint; with it loaded the extension's concrete methods
shadow the stubs (`sdp_stubs.jl:1-50`, `Project.toml:9-26`,
`libaic.jl:86-100`). NO PYTHON.

---

## Notes flagged for verification

Items where the source docs are internally inconsistent, use a notation the
math contradicts, or could read as stale; recorded so they can be checked before
publication.

1. **`P_B` vs `1_B` in the rung-3 round-trip label.** `README.md:192` writes the
   $\eta=0$ oracle maximum as `‖ΥΔ−1_B‖_cb = 3.9e-75`, but FINDINGS §C22(b)
   (`paper/FINDINGS.md:930`) and §C14(b) establish that the correct round-trip
   target is the **conditional expectation `P_B`** (block-diagonal), *not* the
   full identity `1_B` — $\Upsilon\Delta = P_B$, and comparing to the full
   identity is "the exact §C14(b) test that cannot pass" (reads `1.0` on every
   off-block unit even at $\eta=0$). The README's `1_B` notation is therefore at
   best loose; the measured `3.9e-75` is consistent only with the `P_B` target.
   I used `P_B` in section C rung 3. Worth aligning the README label.

2. **`make_mixconj` "out-of-regime" claim.** `CLAUDE.md:421-432` already carries
   the FINDINGS §C15 correction that `make_mixconj*` *cannot* leave the spectral
   basin $\rho(\Phi^2-\Phi) < 1/4$ at any $t$ (measured $\rho\approx0.165$ even at
   $t=0.45$), retracting an earlier "$t=0.45$ is out-of-regime" claim. The
   CLAUDE.md text is self-correcting and current; flagging only so a doc-site
   author does not resurrect the stale claim from older notes. (Separately,
   FINDINGS §C17 notes the *cb-norm* — a different quantity — can reach the $1/4$
   boundary.)

3. **cb-norm "out of scope" comment in `aic_ucp.h` reads as stale.**
   `aic_ucp.h:121-126` says the idempotency defect $\|\Phi^2-\Phi\|_\cb$ "needs
   the diamond-norm SDP (bead aic-d24, out of scope)", but that defect is in fact
   shipped — the eig-free bracket (`aic_cbnorm.h`), the tight MOSEK-fed certifier
   (`aic_cbnorm_certify`), and the Julia `eta_idempotence`/`idempotency_defect`
   all compute it. The "out of scope" phrasing is scoped to the original `ucp`
   module increment and is superseded by `aic_cbnorm.*`; harmless but potentially
   confusing to a reader landing in `aic_ucp.h` first.

4. **Eig-free ratio stated two ways, both correct.** The bracket ratio is given
   as `2n` (`aic_cbnorm.h:21`, `ALGORITHM.md:324`) and as `hi/lo ~ 2n`
   (`README.md:218`). These agree ($\mathrm{hi}/\mathrm{lo} = 2\|J\|_F /
   (\|J\|_F/n) = 2n$); noting it only so the "~2n" in the benchmark table is not
   read as a separate looser claim.
