# Julia package research — su2-fft idioms for AlmostIdempotentChannels.jl
# aic-exa.1 (research leg)

> Read order: CLAUDE.md -> this doc -> then docs/research/julia_package_design.md (aic-exa.2).
> Companion: docs/research/packaging_julia.md (earlier packaging research leg).

## What was studied

- `/home/tobias/Projects/su2-fft/julia/Project.toml`
- `/home/tobias/Projects/su2-fft/julia/src/SU2FFT.jl`
- `/home/tobias/Projects/su2-fft/julia/deps/build.jl` + `deps/deps.jl`
- `/home/tobias/Projects/su2-fft/julia/test/runtests.jl`
- `/home/tobias/Projects/su2-fft/julia/DESIGN.md`
- `/home/tobias/Projects/su2-fft/include/su2.h` (C ABI shape)
- The existing AIC package at
  `julia/AlmostIdempotentChannels.jl/` (the gap to close).

---

## 1. Idioms worth COPYING (with file:line)

### 1a. Project.toml: minimal deps, extras-only for test-only packages

`su2-fft/julia/Project.toml` has zero non-stdlib `[deps]` — only `Libdl`
(stdlib). `Test`, `LinearAlgebra`, `Random` live entirely under `[extras]` +
`[targets]`. No `[weakdeps]`, no `[extensions]` — because the package has no
optional solver.

```toml
[deps]
Libdl = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[targets]
test = ["Test", "LinearAlgebra", "Random"]
```

For AIC: `Mosek`/`MosekTools`/`Convex` must move to `[weakdeps]` + a package
extension, matching this pattern. `LinearAlgebra` stays as a real dep (used in
the core type system + marshalling), or moves to extras if only needed in tests.

### 1b. build.jl: shells to make, writes absolute path to deps.jl

`su2-fft/julia/deps/build.jl` (all 25 lines):

```julia
const PROJECT_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const SO_PATH      = joinpath(PROJECT_ROOT, "build", "libsu2.so")
run(`make -C $PROJECT_ROOT $MAKE_TARGET`)
isfile(SO_PATH) || error("Build failed: $SO_PATH not found after make")
open(DEPS_FILE, "w") do io
    println(io, "const LIBSU2 = ", repr(SO_PATH))
end
```

AIC's `build.jl` correctly copies this with one extra `..` level (the package
sits three levels deep, not two). The pattern is correct and worth keeping for
the active-development phase. For a registrable release, replace with a JLL
(see packaging_julia.md §"The JLL pipeline").

### 1c. __init__: RTLD_GLOBAL libgmp preload + RTLD_DEEPBIND on the lib

`su2-fft/julia/src/SU2FFT.jl:51-80`:

```julia
function __init__()
    try
        Libdl.dlopen("/lib/x86_64-linux-gnu/libgmp.so.10",
                     Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL)
    catch; end
    flags = Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL | Libdl.RTLD_DEEPBIND
    _LIBSU2_HANDLE[] = Libdl.dlopen(LIBSU2, flags)
    _SYM_FFT[]       = Libdl.dlsym(_LIBSU2_HANDLE[], :su2_fft)
    ...
end
```

This is LOAD-BEARING on this machine: Julia 1.12 bundles a libgmp that conflicts
with the system libgmp that libflint resolves `__gmpn_*` against. AIC already
copies this pattern (with the fallback `Libdl.dlopen("libgmp.so.10", ...)`).
Both must do the RTLD_GLOBAL preload or libflint silently uses the wrong libgmp.

### 1d. Const Ref sentinel pattern for dlsym handles

```julia
const _LIBSU2_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_FFT       = Ref{Ptr{Cvoid}}(C_NULL)
```

Populated in `__init__`. Allows type-stable ccall via `_SYM_FFT[]` without
re-dlsym on every call. Every new ABI entry point added to AIC must follow this
pattern: declare a `Ref{Ptr{Cvoid}}` constant, dlsym in `__init__`.

### 1e. ccall via dlsym handle, NOT string/symbol literal

SU2FFT uses `ccall(_SYM_FFT[], ...)` (the dlsym handle), NOT
`ccall((:su2_fft, LIBSU2), ...)`. This matters because the string/symbol literal
form re-looks up the symbol on every call; the handle is resolved once in
`__init__`. The AIC package already copies this. Keep it.

### 1f. Flat arrays for C interop, no structs crossing the boundary

`su2.h` exports only functions taking `int N`, `double _Complex *f`,
`double _Complex *fhat` — no opaque handles, no FLINT types. The Julia side
allocates `Vector{ComplexF64}` and passes it via `Ptr{ComplexF64}`. This is the
same design as the AIC shims (flat `double *re` / `double *im` split arrays, int
dims). The shim design is correct and should be extended, not replaced.

### 1g. Return-type annotation on every public function

Every exported SU2FFT function carries an explicit return type:

```julia
function total_coeffs(N::Integer)::Int
function wigner_d(l::Integer, n::Integer, m::Integer, theta::Real)::ComplexF64
function fft(f::Array{ComplexF64,3})::Vector{ComplexF64}
```

This keeps ccall return types visible, enables inference, and documents the ABI
contract. AIC must do the same on all public functions.

### 1h. ArgumentError (not @assert) for user-facing validation

```julia
N >= 2 || throw(ArgumentError("N must be >= 2, got $N"))
size(f) == (N, N, N) || throw(ArgumentError("samples must be (N,N,N), got $(size(f))"))
```

`@assert` is stripped by `--check-bounds=no`; `ArgumentError` is always raised.
Use `ArgumentError` for invalid user inputs, `error()` for C return-code != 0.
The AIC `choi_diff` and `eta_eigfree` already do this correctly.

### 1i. Docstring discipline: describe the C convention, not just the signature

SU2FFT docstrings name the C function being wrapped, state the layout convention
(e.g. `f[j2+1, k+1, j1+1] = sample at Euler-grid (j1, k, j2)`), cite the
relevant note/design-doc section, and mention known limitations (e.g. the
Riemann quadrature floor). AIC must do the same: cite the `.tex` line, name the
convention (Convention-A Choi, K-major Kraus), state what the C return code
means.

### 1j. Tests: cross-check structure (independent oracle vs ccall path)

`runtests.jl` for SU2FFT has:
- exact-special-case tests (zero input, delta_l00 -> constant 1)
- cross-check test: `fft` vs `ft_direct` to 1e-10 on random input (mirrors the C
  test `tests/test_fft.c::test_fft_matches_direct_random` verbatim)
- layout/accessor invariants (`fhat_at`, `fhat_block`)
- input validation (ArgumentError paths)
- `@info` reporting concrete numbers (max diff, leakage floor) on every testset

AIC's existing tests (T1-T5) follow the same structure. The pattern is correct;
extend it for each new ABI entry point added in aic-exa.3 onward.

### 1k. DESIGN.md: memory-layout contract documented in one place

`su2-fft/julia/DESIGN.md` pins: the C row-major layout, the Julia column-major
index permutation that avoids a copy, the coefficient flat layout, the ccall
signatures, and the test strategy. AIC should maintain an equivalent
`docs/julia_api_design.md` (or expand ALGORITHM.md) that pins these for the new
ABI entries.

### 1l. libXXX_path() diagnostic function

```julia
libsu2_path() = LIBSU2
```

Trivial, but useful when debugging linkage. AIC already has `libaic_path()`.
Keep it in all future versions.

---

## 2. What NOT to copy (dated/fragile patterns)

### 2a. deps/build.jl is NOT registrable

`deps/build.jl` (the `Pkg.build()` hook) has been deprecated for registrable
packages since at least Julia 1.6. For the active-development phase it is fine.
Do NOT submit to the General registry until this is replaced by a JLL
(packaging_julia.md §"The JLL pipeline"). The JLL approach also eliminates the
`RTLD_DEEPBIND` + system-libgmp workaround, since BinaryBuilder ships a
consistent libgmp that matches the FLINT_jll.

### 2b. RTLD_GLOBAL on the library itself (fragile)

SU2FFT passes `RTLD_GLOBAL` on the library handle itself (`flags = RTLD_NOW |
RTLD_GLOBAL | RTLD_DEEPBIND`). This puts all libsu2 symbols into the global
symbol table, which can cause silent symbol collisions if two packages do it.
The correct approach for production is to use RTLD_LOCAL on the library and
rely on the JLL's pre-loaded dependency tree. For development (local build), the
current approach works but is platform-fragile. Accept it for now; revisit when
moving to JLL.

### 2c. No package extension / weakdeps mechanism

SU2FFT has no optional solver, so there is no extension. The AIC equivalent
MUST use `[weakdeps]` + `[extensions]` for the MOSEK/Convex/MosekTools path.
Do not copy the "hard deps on everything" pattern from the CURRENT AIC
`Project.toml`. The new AIC `Project.toml` should look more like SU2FFT's
(core deps only in `[deps]`) plus a `[weakdeps]` block for the SDP solver.

### 2d. No Preferences.jl for libpath

SU2FFT hardcodes the absolute path in `deps/deps.jl`. This is fine for
development but locks the path to the build machine. Once a JLL is added,
use `Preferences.jl` to let users override the libpath:

```julia
const LIBAIC = @load_preference("libaic", AlmostIdempotentChannels_jll.libaic)
```

This is the HDF5.jl pattern (cited in packaging_julia.md). Skip for now; add
when the JLL lands.

### 2e. No Aqua.jl / JET.jl quality checks

Neither SU2FFT nor the current AIC package runs Aqua or JET. Both should.
Aqua catches ambiguities, undefined exports, missing compat bounds; JET catches
type-inference gaps. Add them as test extras when the type system is settled
(aic-exa.5 onward).

### 2f. No Base.show, no type system

SU2FFT has no custom types — everything is a `Vector{ComplexF64}` or a plain
scalar. For AIC's richer domain (UCP channels, certified brackets, factorization
results), plain vectors are inadequate. The new AIC package must add:
- `struct UCP` (wraps Kraus operators, validates UCP conditions)
- `struct CertifiedInterval` (the [lo, hi] bracket from the eig-free shim)
- `struct FactorizationResult` (B, Delta, Upsilon Choi matrices, n, n_B, eta proxy)
- `Base.show` for each

There is no pattern in SU2FFT to copy for this; design from scratch in aic-exa.5.

### 2g. ccall over dlsym handles vs @ccall macro

SU2FFT uses `ccall(_SYM[][], RetType, (ArgTypes...), args...)` — the pre-Julia-1.9
4-argument form via a dlsym handle. For new code, the `@ccall` macro (Julia 1.9+,
compat already set to 1.10) is cleaner:

```julia
@ccall _LIBAIC_HANDLE[].aic_fn(arg1::Cint, arg2::Ptr{Cdouble})::Cint
```

However `@ccall` requires a module/function path, not a `Ptr{Cvoid}` directly, so
the dlsym-handle pattern cannot be expressed as `@ccall` without an intermediate
`ccall`. The 4-arg form via dlsym handles is the only way to combine (a) dlopen +
dlsym (required for the RTLD flags) with (b) type-safe ccall. Keep the existing
pattern; do NOT switch to `@ccall` for the dlsym path.

If/when moving to a JLL, `@ccall libaic.fn(...)::Ret` becomes available (no
dlopen needed) and is preferred.

---

## 3. Concrete checklist for aic-exa.2 (master design doc)

### libpath discovery

- [ ] Keep `deps/build.jl` + `deps/deps.jl` for development phase (mirrors su2-fft).
- [ ] Document the three-level path offset in a comment (already in AIC build.jl).
- [ ] Add Preferences.jl hook as prep for JLL migration (no functional change now).

### ccall style

- [ ] All 6 existing shim entry points: keep `Ref{Ptr{Cvoid}}` + dlsym handle pattern.
- [ ] New shim entry points from aic-exa.3: same — one const Ref per symbol, dlsym in `__init__`.
- [ ] Return-type annotation (`::ReturnType`) on every public function.
- [ ] `ArgumentError` for input validation; `error()` for C rc != 0.

### Type system (aic-exa.5)

- [ ] `struct UCP` — wraps `Vector{Matrix{ComplexF64}}` Kraus, validates sum K_a^dag K_a = I.
- [ ] `struct CertifiedInterval{T<:Real}` — `lo::T`, `hi::T`, `method::Symbol`.
- [ ] `struct FactorizationResult` — `B_choi`, `delta_choi`, `upsilon_choi`, `n`, `n_B`, `dim_B`, `eta_proxy`, `eps_used`.
- [ ] `Base.show` for each: one-line "UCP(n=4, r=6)", "CertifiedInterval([0.1234, 0.1237])",
  "FactorizationResult(n=4 -> n_B=4, dim_B=2, eta_proxy=0.031)".
- [ ] No `mutable struct` for value types; `mutable struct` + `finalizer` only for C heap handles.

### Project.toml (aic-exa.4)

- [ ] Move `Mosek`, `MosekTools`, `Convex` from `[deps]` to `[weakdeps]`.
- [ ] Add `[extensions]` block: `MosekExt = ["Mosek", "MosekTools", "Convex"]`.
- [ ] Extension lives in `ext/MosekExt.jl`; re-exports `diamond_norm_watrous` etc.
- [ ] Core package must test GREEN with no MOSEK present.
- [ ] `LinearAlgebra` stays in `[deps]` (used in marshalling + type system).
- [ ] `Libdl` stays in `[deps]` (load path).

### Tests (aic-exa.9)

- [ ] T0 (new): `libaic_path()` returns a readable file.
- [ ] T1: choi_diff C/arb vs pure-Julia (already exists).
- [ ] T2: analytic anchors (already exists).
- [ ] T3: eig-free bracket (already exists — requires no MOSEK).
- [ ] T4: GADC (already exists — requires no MOSEK).
- [ ] T5: strong duality table (MOSEK, move to extension test).
- [ ] T6 (new): `UCP` constructor validation (non-UCP Kraus throws).
- [ ] T7 (new): `factorize` round-trip for the eta=0 oracle.
- [ ] T8 (new): `CertifiedInterval` containment: `lo <= midpoint <= hi`.
- [ ] Add `using Aqua; Aqua.test_all(AlmostIdempotentChannels)` (extras = Aqua).

### Build system

- [ ] `deps/build.jl` shells to `cmake --build build` (not `make lib`) once aic-exa.3 CMake target is set.
- [ ] Or keep `make lib` if Makefile target is reliable — verify both produce the same `libaic.so`.

### Docs (aic-exa.10)

- [ ] Module docstring: cite tex:347 for eta definition, tex:460 for th_main, tex:2730 for th_factorization.
- [ ] Each public function: cite tex line, state Convention-A layout, state what rc=0 means.
- [ ] `DESIGN.md` (or `docs/julia_api_design.md`): pin the flat-array layout, Kraus convention, Choi convention for every shim.

---

## 4. ABI gap: what the new shims must emit (aic-exa.3)

The 6 existing shims only return Choi matrices (for the SDP generator). The
headline pipeline artifacts are NOT yet accessible from Julia. The new ABI
extension must emit, per `aic_factorize_choi_shim_d`, the actual factorization
artifacts:

- `n`, `n_B`, `dim_B` (integers) — already returned via `out_N`, `out_nB`,
  `out_dimB` output pointers.
- `B_choi` (n_B^2 x n_B^2) — the Choi of the C* algebra B (UCP id_B via Choi).
- `Delta_choi` (n^2 x n_B^2) — encoding map Delta: B(H_B) -> B(H).
- `Upsilon_choi` (n_B^2 x n^2) — decoding map Upsilon: B(H) -> B(H_B).
- `eta_proxy`, `eps_used` — already in `out_eta[0]`, `out_eta[1]`.

The Julia side assembles these into a `FactorizationResult` struct. The ccall
signature for the extended shim (to be designed in aic-exa.3) should follow the
same pattern: flat `double *` re/im buffers, int output scalars, caller
pre-allocates to the safe upper bound `n^4`.

---

## 5. Key differences: su2-fft vs AIC (the harder problem)

| Dimension | su2-fft | AIC (new) |
|---|---|---|
| C ABI surface | 15 flat `double _Complex *` functions | 6 shims (flat double), need extension |
| Type system needed? | No: just arrays + scalars | Yes: UCP, CertifiedInterval, FactorizationResult |
| Optional solver | None | MOSEK (weakdep + extension) |
| Math objects | Real-valued arrays, no certified bounds | Certified [lo, hi] brackets, eps-C* algebra artifacts |
| FLINT/arb types in C | Yes (internal only, not crossing boundary) | Yes (internal only, not crossing boundary) |
| Build | deps/build.jl -> make -> libsu2.so | Same (deps/build.jl -> cmake/make -> libaic.so) |
| libgmp workaround | Yes (RTLD_GLOBAL preload) | Yes (already copied) |
| Aqua/JET | No | Add in aic-exa.9 |

The su2-fft pattern is a FOUNDATION, not a complete blueprint. The main additions
for AIC are the type system (aic-exa.5) and the MOSEK package extension (aic-exa.8).
Everything else is a straight port.
