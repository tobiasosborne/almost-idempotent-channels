# Packaging research — Python bindings + polyglot release engineering

> Research leg 3 of 3 (2026-05-31). Sources cited inline. Companion docs:
> `packaging_c_core.md` (the CMake linchpin), `packaging_julia.md`. Bead:
> "Professional polyglot packaging" epic. Part A = Python bindings; Part B =
> monorepo + release engineering across C + Julia + Python.

---

## PART A — Python bindings for the pure-C ABI

### Binding technology: **cffi (API mode, out-of-line)**

| Tool | Verdict |
|---|---|
| nanobind | reject — C++17 only, no pure-C story |
| pybind11 | reject — C++11, forces C++ wrapper glue around a C header |
| Cython | viable; generates `.c` from `.pyx`, fast, but adds a language + regen step |
| **cffi (API, out-of-line)** | **chosen** — needs only a C compiler; paste C decls from `aic.h`; the generated `.c` can be pre-committed so end users build with any C compiler; faster than ABI mode; used by cryptography/PyNaCl |
| ctypes | probes only — per-call libffi overhead, brittle struct marshalling |

Workflow: `python/src/aic/_build.py` does `ffi.cdef(<decls from include/aic/aic.h>)`
+ `ffi.set_source(...)`; the build runs `ffi.compile()` (a CMake custom command)
→ `_aic.<abi>.so`; `aic/__init__.py` wraps it into ergonomic objects. The C header
stays the single source of truth — no `.pyx`, no C++. If numpy buffer integration
is later needed, add a thin Cython layer *on top* of the cffi layer.

### Build backend: **scikit-build-core**

Reuses the **same CMake** that builds the C core (no second build definition,
unlike meson-python which would need a Meson rewrite). Was an early supporter of
free-threaded 3.13t.

```toml
# python/pyproject.toml
[build-system]
requires = ["scikit-build-core>=0.10", "cffi>=1.17"]
build-backend = "scikit_build_core.build"

[project]
name = "almost-idempotent-channels"
version = "0.1.0"
requires-python = ">=3.11"

[tool.scikit-build]
cmake.build-type = "Release"
wheel.packages = ["python/src/aic"]
```

The C core's `libaic` is an `IMPORTED`/`add_subdirectory` CMake target; one
`cmake --build` builds both it and the cffi extension. No double build.

### THE HARD PART: shipping FLINT/arb + LAPACK in a wheel

manylinux images are too old for FLINT 3.x, so **build GMP → MPFR → FLINT from
source in cibuildwheel's `before-all`** into a local prefix, then let
`auditwheel`/`delocate`/`delvewheel` bundle them into the wheel and rewrite
RPATHs. This is exactly the **python-flint** recipe — copy its
`bin/build_dependencies_unix.sh` (downloads source tarballs, `./configure
--enable-shared --prefix=$(pwd)/.local`, `make`, `make install`; sets
`PKG_CONFIG_PATH`). python-flint's bundled wheel (MPFR+FLINT) is ~9 MB.

```toml
[tool.cibuildwheel.linux]
before-all = "bash bin/build_deps_linux.sh"
environment = { PKG_CONFIG_PATH = "$(pwd)/.local/lib/pkgconfig", LD_LIBRARY_PATH = "$(pwd)/.local/lib" }
repair-wheel-command = "auditwheel repair -w {dest_dir} {wheel}"
[tool.cibuildwheel.macos]
before-all = "bash bin/build_deps_macos.sh"
repair-wheel-command = "delocate-wheel --require-archs {delocate_archs} -w {dest_dir} {wheel}"
```

**LAPACK in the wheel:** either bundle a reference/OpenBLAS built in `before-all`
(scipy's approach), or — pragmatic for a lib that uses LAPACK only for its double
path — link against numpy's bundled BLAS at runtime (declare `numpy` a dep).
**Do NOT depend on python-flint's wheel `.so` at runtime** — there is no standard
for cross-wheel native `.so` deps (Python Discourse #23913); auditwheel rejects
unresolved external `.so` refs. Bundle your own GMP/MPFR/FLINT.

Platform matrix (cibuildwheel 3.x): manylinux_2_28 x86_64/aarch64 (+ musllinux),
macOS arm64/x86_64, Windows x86_64 (delvewheel; build deps via MSYS2/MinGW64).
A **conda-forge** feedstock (`meta.yaml` depending on conda-forge `flint`+`lapack`)
is a lower-effort parallel path for conda users.

### Hygiene

- Type stubs: `stubgen` draft → hand-annotate the public API `.pyi`; ship it +
  `py.typed` (PEP 561) in the wheel.
- `uv` for dev (`uv sync`/`uv run pytest`; delegates to scikit-build-core);
  workspaces group multiple Python packages under one lockfile.
- pytest in `python/tests/` asserting **invariants** (η=0 oracle, double-vs-arb),
  not "runs without error". Cross-validate against the C suite.
- Publish via `twine`/`uv publish` on a version tag.

---

## PART B — Polyglot monorepo + release engineering

### Layout (one repo)

```
/  (git tag = version authority)
  VERSION                # single source of truth: "0.1.0"
  CMakeLists.txt         # root: finds GMP/MPFR/FLINT/LAPACK; exports AIC::aic
  include/aic/  src/  tests/        # C core (per packaging_c_core.md)
  julia/AlmostIdempotentChannels.jl/
  python/
    pyproject.toml       # scikit-build-core + cffi
    src/aic/{__init__.py,_build.py,py.typed,_aic.pyi}
    tests/
    bin/build_deps_{linux,macos,windows}.sh
  paper/  docs/  CLAUDE.md ...
```

The root CMake exports `AIC::aic`. Julia `ccall`s the JLL/`Overrides.toml` lib
(no CMake at Julia level). Python's scikit-build-core consumes `AIC::aic` (via
`add_subdirectory` or `find_package(AIC)`). One CMake tree, three consumers — the
Nemo.jl + FLINT_jll + python-flint pattern; Arrow/Z3 "C core as imported target".

### Single-source versioning

Git tag is canonical. CMake `file(READ VERSION …)` before `project()`. Python:
scikit-build-core `metadata.version.provider` (regex on `VERSION`, or
`setuptools_scm` from the git tag). Julia `Project.toml` `version` is the Julia
registry's source of truth — bump it from the tag at release. One edit per
release (Arrow's unified-version model).

### Coordinated release flow

```
1. Bump VERSION (or just tag).
2. C tests:      cmake --build build && ctest --test-dir build -j$(nproc)
3. Julia tests:  julia --project=julia/... -e 'using Pkg; Pkg.test()'
4. Python tests: cd python && uv run pytest
5. git commit -am "chore: release v0.1.0"; git tag v0.1.0; git push --tags
6. cibuildwheel (on tag) → wheels → twine/uv publish
7. Julia: @JuliaRegistrator register
8. bd export > .beads/issues.jsonl && git push   (session-close protocol)
```

Even with local quality gates, the wheel-build+upload can be a tag-triggered
CI job (no test suite remotely — just build + publish).

### Tests

C cross-checks `tests/test_*.c` (CTest); Julia `test/runtests.jl`; Python
`python/tests/`; plus an optional `tests/cross_lang/` that runs the same input
through the C binary, Julia, and Python and asserts agreement (the polyglot
regression).

### Exemplars

- **python-flint** — the wheel-bundling recipe (build FLINT/GMP/MPFR from source
  in `before-all`, bundle via auditwheel). Copy directly.
- **Apache Arrow** — unified release version across languages; C++ core as a
  CMake imported target consumed by all bindings; split languages to separate
  repos only once they diverged (this project is small + cohesive → stay one
  repo).
- **Z3** — one CMakeLists; bindings as subdirectories; only the C API (`z3.h`) is
  exposed to external consumers (mirror: only expose `aic.h`).

---

## Combined decisions

| Question | Answer |
|---|---|
| Python binding | cffi, API mode, out-of-line |
| Python build backend | scikit-build-core (reuses the C-core CMake) |
| FLINT/GMP/MPFR in wheel | build from source in cibuildwheel `before-all`; bundle via auditwheel/delocate (python-flint recipe) |
| LAPACK in wheel | bundle OpenBLAS, or link numpy's BLAS at runtime |
| Cross-wheel native dep | NO — never depend on python-flint's `.so`; bundle your own |
| Monorepo | one repo; root CMake exports `AIC::aic`; `/include /src /tests /julia /python` |
| Versioning | git tag canonical; `VERSION` read by CMake + scikit-build-core; Julia `Project.toml` bumped at release |
| Wheels | cibuildwheel 3.x, manylinux_2_28 |

## Sources

- cffi modes: https://cffi.readthedocs.io/en/stable/overview.html
- nanobind (why not): https://nanobind.readthedocs.io/en/latest/benchmark.html
- scikit-build-core: https://scikit-build-core.readthedocs.io/en/stable/
- SciPy-2024 scikit-build-core: https://proceedings.scipy.org/articles/FMKR8387
- python-flint: https://github.com/flintlib/python-flint
- python-flint build: https://python-flint.readthedocs.io/en/latest/build.html
- cibuildwheel: https://cibuildwheel.pypa.io/en/stable/
- cross-wheel native deps (no standard): https://discuss.python.org/t/native-dependencies-in-other-wheels-how-i-do-it-but-maybe-we-can-standardize-something/23913
- single-source version: https://packaging.python.org/en/latest/discussions/single-source-version/
- Apache Arrow: https://github.com/apache/arrow
- Z3: https://github.com/Z3Prover/z3
- compiled packaging guide: https://learn.scientific-python.org/development/guides/packaging-compiled/
- uv workspaces: https://docs.astral.sh/uv/concepts/projects/workspaces/
- PEP 561 (py.typed): https://peps.python.org/pep-0561/
