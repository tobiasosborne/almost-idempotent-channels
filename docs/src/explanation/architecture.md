# Architecture: C + FLINT + Julia

After this page you will understand the three-layer software stack the package
is built on, why it is C + FLINT + Julia rather than pure Julia, and what you
depend on when you load it.

This page explains the software stack the package is built on; to install and run
it, see [Installation](../getting_started/install.md). It is a companion to [Certified arithmetic](certified_arithmetic.md),
which covers the bracket mathematics and the cross-check ladder in depth — here we
summarise *why* there are two number paths and link there for the ladder.

## The three layers

The package mirrors its `su2-fft` sibling: three layers, each cross-validated
against the next (`CLAUDE.md:108-124`).

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

**Top — the Julia surface.** This is the API you call: `UCPMap`, the five-verb
pipeline, and the `CertifiedBracket` result type. It also holds the Watrous
diamond-norm SDP, which lives in an optional extension (below). It calls down
through one narrow boundary and never touches a FLINT type directly.

**Middle — the flat-double shims (C2–C5).** A thin marshalling layer in
`src/aic_*_shim.c`. It converts between flat `double` arrays and the internal
`acb_mat` / Kraus structures, runs the *same* arb cores the C test suite
exercises, and writes the result back into caller-owned `double` arrays
(`aic_ucp_shim.h:1-23`). It contains no new mathematics.

**Bottom — the C cores, two number paths.** The roughly 85 `src/*.c` files
implement the algorithms (functional calculus, the projection construction,
Choi/Stinespring, the channel factorization). Each routine has a double path on
BLAS/LAPACK and an arb path on FLINT's `acb_mat`/`arb` with explicit precision.

## Why both number paths

The two paths are not redundant; they have different jobs (`CLAUDE.md:121-124`,
verbatim):

> Every nontrivial routine exists (eventually) in **both** number paths, and the
> arb path at `prec=53` must agree with the double path — that is a primary
> cross-check. The arb path is what *certifies* the paper's ``O(\varepsilon)`` bound; the
> double path is what makes the test suite fast.

The **double path** is the fast `prec=53` anchor. It carries the
degeneracy-robust LAPACK eigensolvers, because the generic Choi spectrum is
degenerate and a certified simple-eigenvalue arb routine would abort on it
(`aic_ucp.h:44-51`). It makes the 219-test solver-free suite run in about 2m09s.

The **arb path** is the reason the project exists. It returns rigorous error
balls, so the package can certify Kitaev's ``O(\eta)`` inequalities rather than
trust `double` (`CLAUDE.md:34-36`). `double` is adequate almost everywhere, but
the *bound itself* fails silently in three regimes — a difference of ``O(1)``
quantities like the associator ``(XY)Z - X(YZ)`` or ``\Phi^2-\Phi`` where
cancellation eats the answer; the inversion of a near-singular operator; and a
spectral gap comparable to ``\eta`` where a `double` verdict on rank or positivity
would be a guess (`CLAUDE.md:327-331`). In those regimes the arb ball tracks the
true uncertainty, and an enclosure that has lost all precision is a loud failure,
not a silent widening. The two paths agreeing at `prec=53` is the cheapest rung
of the [Certified arithmetic](certified_arithmetic.md) cross-check ladder; the arb path's outward-rounded ball is the certificate at the top of it.

## Why C + FLINT, not pure Julia

The certified arithmetic is the constraint that fixes the language choice. The
rigorous error balls come from FLINT's `arb`/`acb` types — interval-style
arbitrary precision with provable enclosures — which is a C library. The cores
are written in C against that library directly (`CLAUDE.md` Law 2, Law 15:
"Don't replace FLINT"). The double path uses the same BLAS/LAPACK that any
numerical stack does. Julia sits on top as the high-level surface and the
cross-check harness, exactly the `su2-fft` arrangement. There is no Python
anywhere in the stack.

## The ABI: flat-double `@ccall` shims

Julia's `@ccall` passes `int`, `double`, and `Ptr{Cdouble}` cleanly, but **not**
FLINT structs: `acb_mat_t` is an opaque array type and `slong` is FLINT's
machine-word integer. So the language boundary carries `int` dimensions and
`Ptr{Cdouble}` arrays only — no `acb_mat_t` and no `slong` ever cross it. The
shim does the conversion at the boundary (`aic_ucp_shim.h:1-23`,
`ALGORITHM.md:601-639`). The seven shim entry points dlsym'd into Julia are
`aic_ucp_choi_diff_d`, `aic_cbnorm_eigfree_d`, `aic_cbnorm_certify_d`,
`aic_assoc_summary_d`, `aic_main_iso_summary_d`,
`aic_factorize_artifacts_sizes_d`, and `aic_factorize_artifacts_d`
(`libaic.jl:67-84`). The flat array layout is row-major:
`kraus_re[a*n*n + i*n + j]` for `r` operators each `n × n`
(`aic_ucp_shim.h:17-22`).

## The Kraus storage contract — Heisenberg picture

The storage layout is load-bearing and pinned in the header (`aic_ucp.h:6-21`). A
UCP map ``\Phi`` is stored as ``r`` Kraus operators ``K_a``, and its action is the
**observable (Heisenberg)** form (`approximate_algebras.tex:1570`):

```math
\Phi(X) = \sum_a K_a^\dagger X K_a,
```

with unitality ``\sum_a K_a^\dagger K_a = 1`` (`approximate_algebras.tex:1571`).
This is **not** the Schrödinger/state dual ``T = \Phi^*(X) = \sum_a K_a X
K_a^\dagger`` (`aic_ucp.h:18-21`). Getting the dual backwards silently transposes
everything, so the picture is fixed once, at the storage layer, and every routine
above inherits it. The dual *channel* on states is recovered explicitly where it
is needed (the `encode`/`decode` maps), never by accident.

## The Convention-A Choi matrix

Many routines work through the **Choi matrix**, in the convention shared by the
paper, Watrous, QETLAB, and Qiskit (`approximate_algebras.tex:1575`): with an
orthonormal basis ``\{e_i\}`` and ``E_{ij} = |e_i\rangle\langle e_j|``,

```math
C_\Phi = \sum_{ij} E_{ij} \otimes \Phi(E_{ij}),
```

with the input factor on the **left/major** index. Deriving ``\Phi(E_{ij})``
from the Heisenberg action gives the per-entry formula pinned in the header
(`aic_ucp.h:30-39`):

```
C[i*n + a, j*n + b] = sum_x conj(K_x[i,a]) * K_x[j,b].
```

The conjugation sits on the **first** ``(i,a)`` factor. This is not cosmetic: it
is exercised only by the complex asymmetric `complex_qubit` channel, and dropping
it makes the C Choi diverge from the Julia reference by `1.28` in operator norm,
far above the `1e-11` fixture tolerance (`ALGORITHM.md:264-266`,
`aic_ucp.h:104-107`). With this convention, ``\Phi`` is completely positive iff
``C_\Phi \succeq 0`` and unital iff ``\operatorname{Tr}_K(C_\Phi) = 1`` (tracing
out the left factor). On an exactly-idempotent ``\Phi`` the defect Choi
``\operatorname{Choi}(\Phi^2-\Phi)`` collapses to zero — see
[The η = 0 oracle](../tutorials/eta0_oracle.md).

## The native-library lifecycle

The Julia package loads `libaic.so` through three mechanisms, in order.

**Discovery via Preferences.** The default native-library path is
`build/libaic.so` relative to the repo root, so a single CMake build needs no
`Pkg.build`. To point elsewhere,
`AlmostIdempotentChannels.set_libaic_path!("/abs/path/libaic.so")` writes a
`LocalPreferences.toml` entry read at precompile time (`libaic.jl:7-12,27-58`).

**The libgmp preload.** On load, `__init__` first pre-loads the system
`libgmp.so.10` with `RTLD_NOW | RTLD_GLOBAL` (full path, then bare-name fallback,
in a try/catch) **before** libaic and libflint are mapped (`libaic.jl:104-135`).
This is the `su2-fft` fix for a clash between Julia's bundled libgmp and the one
FLINT was built against: libflint resolves its `__gmpn_*` symbols against
whichever libgmp is already loaded, so pre-loading the system copy keeps its
definitions first in the global symbol table (`libaic.jl:111-127`,
`ALGORITHM.md:699-708`).

**The dlopen and symbol binding.** Then `dlopen(LIBAIC_PATH, RTLD_NOW |
RTLD_GLOBAL | RTLD_DEEPBIND)` maps the library, and each shim symbol is `dlsym`'d
into a module-private `Ref{Ptr{Cvoid}}`. `RTLD_DEEPBIND` on libaic alone is
*insufficient* for the libgmp clash — that is why the preload step exists — but
it is still used so libaic prefers its own dependency symbols.

## Solver-free by default

The package needs **no SDP solver** for its default operation: the eig-free
certified bracket and the flat-double ccall surface stand alone. Only the Watrous
diamond-norm SDP — the exact `idempotency_defect`, the tight bracket, the
`diamond_*` entry points — requires MOSEK, and that lives **only** in the
optional package extension `AICMosekExt` (MOSEK is a `[weakdeps]`,
`Project.toml:9-26`). Without the extension loaded, those functions throw a
helpful error registered through `Base.Experimental.register_error_hint`; with
`using Convex, Mosek, MosekTools` loaded, the extension's concrete methods shadow
the stubs (`sdp_stubs.jl:1-50`, `libaic.jl:86-100`). The dual-path rationale for
*why* the solver-free bracket is loose but rigorous — and how MOSEK tightens it —
is in [Certified arithmetic](certified_arithmetic.md).

## See also

- [Certified arithmetic](certified_arithmetic.md) — the `CertifiedBracket` inequality, outward
  rounding, the eig-free Watrous bracket, and the cross-check ladder
- [The mathematics](math_story.md) — what the five verbs compute and why
- [The constructive mandate](constructive.md) — why the cores are finite-dim algorithms
- [Installation](../getting_started/install.md) — building `libaic` and pointing Julia at it
