# Julia package best-practice survey — AlmostIdempotentChannels.jl refactor
# 2026-06-02

> Scope: current (2024-2026) best practice for a Julia package that wraps a
> C library with heavy native dependencies (FLINT/arb + LAPACK), with an
> optional heavyweight solver (MOSEK/Watrous SDP) behind a package extension.
> Primary sources: Pkg.jl docs, BinaryBuilder.jl docs, Arblib.jl, Nemo.jl,
> Preferences.jl, Aqua.jl, JET.jl, Documenter.jl.

---

## 1. Package extensions + [weakdeps] (Julia >= 1.9): making MOSEK optional

### 1.1 Canonical Project.toml layout

The entire Watrous SDP (Convex + Mosek + MosekTools) moves to a package
extension. The main package keeps only what is needed without a solver.

```toml
name = "AlmostIdempotentChannels"
uuid = "6bf14b85-acf1-40b3-a7a0-aef96f56c734"
version = "0.2.0"

[deps]
Libdl      = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Preferences = "21216c6a-2e73-6563-6e65-726566657250"

# ---------- optional solver (Watrous diamond-norm SDP) ----------
[weakdeps]
Convex      = "f65535da-76fb-5f13-bab9-19810c17039a"
Mosek       = "6405355b-0ac2-5fba-af84-adbd65488c0e"
MosekTools  = "1ec41992-ff65-5c91-ac43-2df89e9693a4"

[extensions]
AICMosekExt = ["Convex", "Mosek", "MosekTools"]

[compat]
julia      = "1.9"
Preferences = "1"
# compat entries for weakdeps still belong here:
Convex      = "0.16"
Mosek       = "11"
MosekTools  = "0.15, 0.16"

# ---------- test-only ----------
[extras]
Test          = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
Random        = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Aqua          = "4c88cf16-eb10-579e-8560-4a9242c79595"
JET           = "c3a54625-cd67-489e-a8e7-0a5a0ff4e31c"

[targets]
test = ["Test", "Random", "Aqua", "JET"]
```

Key rules (Pkg.jl docs, pkgdocs.julialang.org/v1/creating-packages/):
- `[weakdeps]` lists optional packages — they are NOT installed automatically.
  Version constraints still belong in `[compat]`.
- `[extensions]` maps the extension module name to a list (or single string)
  of triggering weakdeps. All listed weakdeps must be loaded for the extension
  to activate.
- Extension tests: add the triggering weakdeps to `[extras]` + `[targets]`
  so `Pkg.test()` installs them and the extension activates during the test run.
  (For Julia >= 1.12 with workspace support, use `test/Project.toml` instead.)
- Backward compat (pre-1.9): if you must support Julia < 1.9, also add the
  weakdeps to `[extras]` and use `Requires.jl` with
  `!isdefined(Base, :get_extension)` guard. For this package Julia >= 1.9 is
  the minimum, so Requires.jl is not needed.

### 1.2 Extension file skeleton

`ext/AICMosekExt.jl`:

```julia
module AICMosekExt

# All three triggering packages must be loaded for this module to be active.
import AlmostIdempotentChannels
using Convex, Mosek, MosekTools, LinearAlgebra

# Extend the stub defined in the main package (see §1.3).
function AlmostIdempotentChannels.diamond_norm_watrous(
        J::Matrix{ComplexF64}, n::Int)
    # ... (current sdp.jl implementation moves here verbatim) ...
end

function AlmostIdempotentChannels.diamond_norm_watrous_primal(
        J::Matrix{ComplexF64}, n::Int)
    # ...
end

function AlmostIdempotentChannels.diamond_norm_dual(
        J::Matrix{ComplexF64}, n::Int)
    # ...
end

function AlmostIdempotentChannels.diamond_dual_rect(
        J::Matrix{ComplexF64}, d_maj::Int, d_min::Int, tr_sys::Int)
    # ...
end

function AlmostIdempotentChannels.diamond_primal_rect(
        J::Matrix{ComplexF64}, d_maj::Int, d_min::Int, rho_on::Symbol)
    # ...
end

end # module AICMosekExt
```

Naming convention: `<Package><Trigger>Ext` or `<Package><Feature>Ext`.
Module name must match the filename exactly.

### 1.3 Stubs in the main package with informative error hints

In `src/AlmostIdempotentChannels.jl`, declare the functions with a fallback
that fires a human-readable error when the extension is not loaded:

```julia
# Stubs — implementations live in ext/AICMosekExt.jl.
# error_hint approach from Base.Experimental (discourse.julialang.org/t/104156):

function diamond_norm_watrous(J::Matrix{ComplexF64}, n::Int)
    error("""
    diamond_norm_watrous requires Convex, Mosek, and MosekTools.
    Load them first:
        using Convex, Mosek, MosekTools
    """)
end
# Repeat for diamond_norm_watrous_primal, diamond_norm_dual, etc.
```

Alternative: use `Base.Experimental.register_error_hint` on `MethodError`
to attach the message to the Julia error system rather than defining a
catch-all method. That keeps the stub namespace clean but is more elaborate.
The explicit-error-method pattern is simpler and sufficient here.

### 1.4 Testing extensions

In `test/runtests.jl`:

```julia
# Tests T1-T4 run with no solver (exercise only the C/arb path).
# T5 (strong duality / dual normalization) requires the solver extension.
# Structure: try to load the solver; skip gracefully if not available.

using Test, AlmostIdempotentChannels

@testset "Core (no solver)" begin
    include("test_core.jl")   # T1-T4
end

@testset "Watrous SDP (requires MOSEK extension)" begin
    try
        using Convex, Mosek, MosekTools
        include("test_sdp.jl")  # T5 (current content)
    catch e
        @warn "Skipping SDP tests: $(e.msg)"
        @test_skip "SDP tests skipped (Mosek not available)"
    end
end
```

In CI / `Pkg.test()` with the solver extension available, add Convex/Mosek/
MosekTools to `[extras]`+`[targets]` (already shown above) so the extension
activates automatically.

---

## 2. Preferences.jl for libaic.so path override

### 2.1 Current problem

`deps/build.jl` shells out to `make` and writes an absolute path to
`deps/deps.jl`, which is then `include`d at package load. This pattern
(from 2016-era Pkg2) is brittle: the path is baked at `Pkg.build` time, the
.so path is never in version control, and there is no way to override it
without re-running build.

### 2.2 Idiomatic 2024 replacement

Use Preferences.jl (`@load_preference`) at package load time. Preferences are
stored in `LocalPreferences.toml` (next to the active `Project.toml`), which
is `.gitignore`d. The preference is read at precompilation time — a change
requires restarting Julia to take effect.

In `src/AlmostIdempotentChannels.jl`:

```julia
using Preferences

# Compile-time preference: changing requires Julia restart.
# Default: the build/libaic.so path relative to the package root,
#          resolved at precompile time.
const _DEFAULT_LIBAIC = let
    repo = abspath(joinpath(@__DIR__, "..", "..", ".."))
    joinpath(repo, "build", "libaic.so")
end

const LIBAIC_PATH = @load_preference("libaic_path", _DEFAULT_LIBAIC)

function __init__()
    isfile(LIBAIC_PATH) || error(
        "libaic.so not found at \"$LIBAIC_PATH\".\n" *
        "Either run `cmake --build build` in the repo root, " *
        "or set the path with AlmostIdempotentChannels.set_libaic_path!(path).")
    # ... dlopen as before ...
end
```

Expose a convenience setter:

```julia
"""
    set_libaic_path!(path::AbstractString)

Override the path to libaic.so. Writes to LocalPreferences.toml in the
active project. Restart Julia for the change to take effect.

Example:
    AlmostIdempotentChannels.set_libaic_path!("/opt/aic/lib/libaic.so")
"""
function set_libaic_path!(path::AbstractString)
    @set_preferences!("libaic_path" => string(path))
    @info "libaic_path set to \"$path\". Restart Julia for this to take effect."
end
```

`LocalPreferences.toml` format (user-editable directly):

```toml
[AlmostIdempotentChannels]
libaic_path = "/home/tobias/Projects/almost-idempotent-channels/build/libaic.so"
```

Sources: Preferences.jl README (github.com/JuliaPackaging/Preferences.jl);
JuMP custom binaries guide (jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/).

### 2.3 Removing deps/build.jl entirely

Once Preferences.jl is in place, `deps/build.jl` and `deps/deps.jl` can be
deleted. The cmake build (`cmake --build build`) is already documented in
`BUILDING.md`; users run it once, then the path is either the discoverable
default or explicitly set via `set_libaic_path!`. This eliminates the
`Pkg.build` callback and the stale-path class of bugs.

---

## 3. JLL / BinaryBuilder + Overrides.toml for local dev

### 3.1 Where a JLL fits in the long run

BinaryBuilder.jl + Yggdrasil submission is the correct end state for public
release: a `AIC_jll` package contains the pre-built `libaic.so` for each
platform, installed by Pkg like any other binary. Then the package's `[deps]`
gains `AIC_jll` and the dlopen call becomes `using AIC_jll; libaic` (the JLL
exports the path constant).

### 3.2 Minimal local story (pre-Yggdrasil)

Three mechanisms, from lightest to heaviest:

**A. Preferences.jl path override (recommended for now, see §2).**
No JLL, no BinaryBuilder. User builds with cmake, sets path once.
Zero friction for the primary developer.

**B. dev'd JLL override directory.**
Stub a minimal JLL package locally (`Pkg.dev`), create
`~/.julia/dev/AIC_jll/override/lib/libaic.so` pointing at the cmake output
(symlink or copy). The JLL's `dev_jll()` function sets this up.
Slightly more infrastructure; good if you want to test the JLL loading path.

**C. Artifacts Overrides.toml.**
`~/.julia/artifacts/Overrides.toml` replaces an entire artifact tree by UUID:

```toml
# UUID from AIC_jll's Project.toml (when it exists)
[<artifact-uuid>]
AIC = "/home/tobias/Projects/almost-idempotent-channels/build"
```

The path must replicate the artifact directory structure (`lib/libaic.so`,
`include/`, etc.). Heavy but robust for testing cross-platform packaging.

Source: BinaryBuilder.jl JLL docs (docs.binarybuilder.org/stable/jll/).

**Decision for this project:**
Stick with option A (Preferences.jl) until a public release is in scope.
Document option C in a `PACKAGING.md` note for when BinaryBuilder submission
happens.

---

## 4. Library loading: @ccall, dlopen, GC.@preserve, finalizers

### 4.1 @ccall vs ccall

Julia 1.5+ provides the `@ccall` macro as cleaner syntax. For dlsym-obtained
function pointers (the current pattern, load-bearing because of the RTLD
flags), `ccall(ptr, RetType, (ArgTypes...,), args...)` is still required —
`@ccall` does not accept a `Ptr{Cvoid}` directly. Keep `ccall` for the dlsym
path; use `@ccall` for any new direct-by-name calls.

```julia
# Current pattern (correct for dlsym-obtained pointers — keep as-is):
rc = ccall(_SYM_CHOI_DIFF_D[], Cint, (Ptr{Cdouble}, ...), args...)

# For a future direct-by-name call (if the libpath is stable enough):
@ccall libaic.aic_factorize_choi_shim_d(...)::Cint
```

### 4.2 GC.@preserve

Use `GC.@preserve` whenever you pass a Julia-managed buffer (a `Vector` or
`Array`) as a `Ptr` to C and the C function might allocate or trigger GC:

```julia
GC.@preserve kr ki choi_re choi_im begin
    rc = ccall(_SYM_CHOI_DIFF_D[], Cint,
               (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                Cint, Cint, Cint),
               choi_re, choi_im, kr, ki, n, r, prec)
end
```

The current code passes `Vector{Float64}` buffers directly to `ccall`; Julia
pins them for the duration of the ccall automatically (the `unsafe_convert`
chain), so GC.@preserve is not strictly needed for synchronous ccalls that do
not call back into Julia. It is still good practice for documentation clarity
and future-proofing if the C function is ever called from a threaded context.

### 4.3 Opaque handle structs + finalizers

For any future C objects that must be explicitly freed (e.g., `aic_plan_t*`
once a pipeline plan struct is exposed in the ABI), the standard pattern is a
mutable struct + finalizer:

```julia
mutable struct AICPlan
    ptr::Ptr{Cvoid}   # opaque C handle; C_NULL after free
    function AICPlan(ptr::Ptr{Cvoid})
        ptr == C_NULL && throw(ArgumentError("null AICPlan handle"))
        obj = new(ptr)
        finalizer(obj) do p
            if p.ptr != C_NULL
                ccall((:aic_plan_free, LIBAIC_PATH), Cvoid, (Ptr{Cvoid},), p.ptr)
                p.ptr = C_NULL
            end
        end
        return obj
    end
end

# Deterministic cleanup (preferred over relying on GC):
Base.close(p::AICPlan) = (p.ptr != C_NULL &&
    ccall((:aic_plan_free, LIBAIC_PATH), Cvoid, (Ptr{Cvoid},), p.ptr);
    p.ptr = C_NULL)
```

Rules (from discourse.julialang.org/t/86233 and Vulkan.jl wrapper docs):
- Finalizers ONLY work on mutable structs.
- Finalizers may run on any thread (Julia 1.9+ parallel GC). If the C free
  function is not thread-safe, guard with a lock or use `close()` as the
  primary cleanup path and treat the finalizer as a safety net only.
- Do NOT call back into Julia from a finalizer (will deadlock).
- The deterministic `close()` / `do`-block pattern is always safer than
  relying on the GC for external resources.

### 4.4 Thread safety and FLINT

FLINT is thread-safe per-context but shares global state for some caches.
The existing RTLD_DEEPBIND + system libgmp RTLD_GLOBAL preload is the
correct workaround for the Julia bundled-libgmp ABI mismatch (documented in
the su2-fft bead su2fft-e5z pattern). Keep it.

For the Julia GC: the current single-threaded serial-Julia-process constraint
(CLAUDE.md) means parallel-GC finalizer races are not a current concern.
If Julia threading is ever enabled, add `gc_safe = true` annotation to the
heavy C calls that do not call back into Julia:

```julia
# Future threaded annotation (not needed now):
rc = ccall(_SYM_CHOI_DIFF_D[], Cint, ...; gc_safe = true)
```

---

## 5. Arblib.jl pattern as a reference implementation

Arblib.jl v1.7.0 (github.com/kalmarek/Arblib.jl) is the closest public
analogue to this project — it wraps FLINT/arb with certified arb balls,
exposes a low-level ccall layer and a high-level Julia interface. Key decisions
worth copying:

- **`[deps]`**: `FLINT_jll` (stdlib), `LinearAlgebra`, `ScopedValues`,
  `SpecialFunctions`. No `deps/build.jl`.
- **`[weakdeps]`**: `ForwardDiff` only (for AD support). Extension:
  `ArblibForwardDiffExt`.
- **Minimum Julia**: 1.10 (safe baseline; Julia 1.9 is the extension minimum,
  but 1.10 ships parallel GC which is relevant).
- **Test deps**: Aqua, Documenter, ForwardDiff, Test.
- The JLL (`FLINT_jll`) is the idiomatic binary dependency — no `build.jl`.
  Until `AIC_jll` exists, Preferences.jl fills the same role.
- Low-level wrappers: one file per C module (e.g., `src/arb.jl` wraps arb_*),
  with auto-generated ccall wrappers from parsing C headers (Clang.jl pattern).
  For AIC the six shim entry points are few enough to write manually.

---

## 6. Rich types and idiomatic API

### 6.1 A certified interval type

The C eig-free bracket returns `(lo::Float64, hi::Float64)`. Returning a named
struct instead of a bare tuple is more self-documenting and supports `Base.show`.

Recommendation: define a minimal local `CertifiedBound` type. Do NOT add
`IntervalArithmetic.jl` as a hard dependency — it is IEEE 1788-2015 compliant
but large and pulls in many dependencies. The tiny local type is sufficient.

```julia
"""
    CertifiedBound{T<:AbstractFloat}

A rigorous two-sided bracket [lo, hi] on a scalar quantity, produced by the
eig-free C certifier. The C shim rounds outward (FLOOR/CEIL) so lo <= true_value
<= hi is guaranteed within the arb working precision.
"""
struct CertifiedBound{T<:AbstractFloat}
    lo::T
    hi::T
    function CertifiedBound(lo::T, hi::T) where T<:AbstractFloat
        lo <= hi || throw(ArgumentError("CertifiedBound: lo=$lo > hi=$hi"))
        new{T}(lo, hi)
    end
end

Base.in(x::Real, b::CertifiedBound) = b.lo <= x <= b.hi
Base.width(b::CertifiedBound) = b.hi - b.lo
midpoint(b::CertifiedBound) = (b.lo + b.hi) / 2

function Base.show(io::IO, b::CertifiedBound{T}) where T
    print(io, "CertifiedBound{$T}(", b.lo, " ≤ x ≤ ", b.hi, ")")
end
function Base.show(io::IO, ::MIME"text/plain", b::CertifiedBound{T}) where T
    print(io, "CertifiedBound [$( b.lo), $( b.hi)]  (width=$(width(b)))")
end
```

`eta_eigfree` then returns `CertifiedBound{Float64}` instead of a bare tuple.

### 6.2 Keyword defaults and error messages

All public functions should:
1. Accept keyword arguments with sensible defaults (`prec::Int = 106`).
2. Validate inputs early and raise `ArgumentError` (not generic `error`).
3. Include the function name and the violating value in the message.

```julia
function eta_eigfree(kraus::Vector{Matrix{ComplexF64}};
                     prec::Int = 106)::CertifiedBound{Float64}
    isempty(kraus) && throw(ArgumentError("eta_eigfree: need >= 1 Kraus operator"))
    prec >= 2      || throw(ArgumentError("eta_eigfree: prec=$prec, must be >= 2"))
    n = size(kraus[1], 1)
    all(K -> size(K) == (n, n), kraus) ||
        throw(DimensionMismatch("eta_eigfree: all Kraus operators must be $n × $n"))
    # ...
end
```

### 6.3 Base.show for the module-level handle

The current `libaic_path()` diagnostic function is good. Add `Base.show` for
any future opaque-handle structs (§4.3).

---

## 7. Documenter.jl, doctests, @autodocs

### 7.1 Setup

`docs/Project.toml` (separate environment):

```toml
[deps]
AlmostIdempotentChannels = "6bf14b85-acf1-40b3-a7a0-aef96f56c734"
Documenter = "e30172f5-a6a5-5a46-863b-614d45cd2de4"

[compat]
Documenter = "1"
```

`docs/make.jl`:

```julia
using Documenter, AlmostIdempotentChannels

DocMeta.setdocmeta!(AlmostIdempotentChannels, :DocTestSetup,
    :(using AlmostIdempotentChannels); recursive = true)

makedocs(
    modules   = [AlmostIdempotentChannels],
    sitename  = "AlmostIdempotentChannels.jl",
    pages     = ["Home" => "index.md", "API" => "api.md"],
    doctest   = true,
)
```

`docs/src/api.md`:

```markdown
# API Reference

```@autodocs
Modules = [AlmostIdempotentChannels]
```
```

### 7.2 Doctest format

REPL-style is preferred for short examples (verified by Documenter):

```julia
"""
    eta_eigfree(kraus; prec=106) -> CertifiedBound{Float64}

...

# Example
```jldoctest
julia> using AlmostIdempotentChannels, LinearAlgebra
julia> K = [Matrix{ComplexF64}(I,2,2)];
julia> b = eta_eigfree(K);
julia> 0.0 in b
true
```
"""
```

To fix failing doctests after API changes: `doctest(AlmostIdempotentChannels, fix=true)`.
Always inspect the diff before committing.

### 7.3 DocumenterVitepress

DocumenterVitepress (github.com/LuxDL/DocumenterVitepress.jl) produces a
Vitepress-based site instead of the Documenter default. It is production-quality
for new packages (e.g., SciML ecosystem). For this project the default
Documenter HTML backend is sufficient and has zero extra dependencies;
defer Vitepress unless the docs site needs richer interactivity.

---

## 8. Quality gates: Aqua.jl and JET.jl

### 8.1 Aqua.jl

Aqua.test_all runs 9 checks (juliatesting.github.io/Aqua.jl/stable/test_all/):
1. `test_ambiguities` — conflicting method definitions
2. `test_unbound_args` — unbound type parameters
3. `test_undefined_exports` — exports without definitions
4. `test_project_extras` — Project.toml extras consistency
5. `test_stale_deps` — unused [deps] entries
6. `test_deps_compat` — all [deps] have [compat] entries
7. `test_piracies` — type piracy (extending foreign types with foreign args)
8. `test_persistent_tasks` — tasks that block precompilation
9. `test_undocumented_names` — public names without docstrings (default off)

Minimum snippet for `test/runtests.jl`:

```julia
using Aqua

@testset "Aqua" begin
    Aqua.test_all(AlmostIdempotentChannels;
        ambiguities = true,
        stale_deps  = true,
        deps_compat = true,
        piracies    = true,
        undocumented_names = false,   # enable once all public API is docstring'd
    )
end
```

### 8.2 JET.jl

JET uses the compiler's type inference to find type errors and instabilities.
The 2024 `report_package` API requires loading the target module first:

```julia
using JET

@testset "JET" begin
    # report_package finds all method definitions and checks them.
    # ignored_modules suppresses false-positives from Base/stdlib.
    rep = JET.report_package(AlmostIdempotentChannels;
        ignored_modules = (Base, LinearAlgebra))
    @test length(JET.get_reports(rep)) == 0
end
```

Note: JET is most useful after the codebase is type-stable. Run it
opportunistically rather than as a hard gate until the public API settles.
Use `@report_opt` interactively in the REPL to find specific hot spots.

### 8.3 Test organization recommendation

```
test/
  runtests.jl        # entry: runs all testsets, skips SDP if no MOSEK
  test_core.jl       # T1-T4 (no solver)
  test_sdp.jl        # T5 (Watrous SDP, requires extension)
  test_aqua.jl       # Aqua + JET
```

---

## 9. Summary of changes required

| Item | Current state | Target state |
|------|---------------|--------------|
| MOSEK hard-dep | `[deps]` Convex+Mosek+MosekTools | `[weakdeps]` + ext/AICMosekExt.jl |
| Library loading | `deps/build.jl` + `deps/deps.jl` | `Preferences.jl` + `@load_preference` |
| T1-T4 run without MOSEK | No (all tests import Convex/Mosek) | Yes |
| T5 (SDP) | Always required | Optional, skipped gracefully |
| `eta_eigfree` return type | `Tuple{Float64,Float64}` | `CertifiedBound{Float64}` |
| Quality gates | None | Aqua.test_all + JET.report_package |
| Doctests | None | jldoctest blocks on all public functions |
| Minimum Julia | 1.10 | 1.9 (extensions) or keep 1.10 (parallel GC) |

---

## 10. Concrete files to create/modify

```
julia/AlmostIdempotentChannels.jl/
  Project.toml                 # add [weakdeps]+[extensions]; drop Convex/Mosek/MosekTools from [deps]
  src/AlmostIdempotentChannels.jl  # Preferences libpath; CertifiedBound; SDP stubs
  ext/AICMosekExt.jl           # new: SDP implementations (current sdp.jl content)
  src/sdp.jl                   # delete (content moves to ext/)
  deps/build.jl                # delete
  deps/deps.jl                 # delete (generated; gitignored anyway)
  test/runtests.jl             # split into test_core + test_sdp; add Aqua/JET
  test/test_core.jl            # new: T1-T4
  test/test_sdp.jl             # new: T5 (current test_sdp content)
  test/test_aqua.jl            # new: Aqua + JET
```

---

## Sources

- [Pkg.jl Creating Packages — Package Extensions](https://pkgdocs.julialang.org/v1/creating-packages/)
- [Package Extensions blog post (GLCS)](https://blog.glcs.io/package-extensions)
- [Extensions/weakdeps design discussion (Discourse)](https://discourse.julialang.org/t/extensions-weakdeps-loading-and-package-design/104156)
- [Testing a package extension (Discourse)](https://discourse.julialang.org/t/testing-a-package-extension/107846)
- [BinaryBuilder.jl — JLL packages](https://docs.binarybuilder.org/stable/jll/)
- [JuMP custom solver binaries (Preferences + Overrides)](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
- [Preferences.jl README](https://github.com/JuliaPackaging/Preferences.jl/blob/master/README.md)
- [Arblib.jl](https://github.com/kalmarek/Arblib.jl)
- [Nemo.jl developer introduction](https://nemocas.github.io/Nemo.jl/dev/developer/introduction/)
- [Aqua.jl test_all reference](https://juliatesting.github.io/Aqua.jl/stable/test_all/)
- [JET.jl tutorial](https://aviatesk.github.io/JET.jl/dev/tutorial/)
- [Documenter.jl doctests](https://documenter.juliadocs.org/stable/man/doctests/)
- [GC / finalizer C wrapper discussion (Discourse)](https://discourse.julialang.org/t/best-way-to-have-gc-manage-freeing-c-allocated-storage/86233)
- [Julia C interface docs](https://docs.julialang.org/en/v1/base/c/)
- [Vulkan.jl wrapper types](https://juliagpu.github.io/Vulkan.jl/dev/reference/wrapper_types/)
