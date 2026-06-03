# Building almost-idempotent-channels

The C core is built with **CMake** (`CMakeLists.txt` is the canonical build).
The top-level `Makefile` is a thin convenience wrapper that forwards to
`cmake`/`ctest` so the old `make` / `make test` commands still work; it is not
the build recipe.

## Dependencies

| Component | Minimum | Why | Debian/Ubuntu package |
|---|---|---|---|
| C11 compiler | — | `-std=c11`, no GNU extensions | `gcc` / `clang` |
| CMake | 3.24 | the build | `cmake` |
| FLINT (incl. arb/acb) | 3.0 | certified arbitrary-precision path; arb is *inside* FLINT 3.x, there is no separate libarb | `libflint-dev` |
| LAPACKE | — | C interface to LAPACK used by the fast double path (`LAPACKE_zheev`, `zgesvd`) | `liblapacke-dev` |
| LAPACK | — | Fortran LAPACK (not provided by LAPACKE) | `liblapack-dev` |
| BLAS | — | dense linear algebra under LAPACK | `libblas-dev` |
| Ninja | optional | faster generator; NOT used by default (see Gotchas) | `ninja-build` |

```sh
sudo apt install build-essential cmake libflint-dev liblapacke-dev liblapack-dev libblas-dev
```

## Build

```sh
cmake -S . -B build
cmake --build build
```

or, via the wrapper:

```sh
make            # configure + build
```

Both produce `build/libaic.so` (versioned, SONAME `libaic.so.0`) and
`build/libaic.a`. The ~85 `src/*.c` are compiled **exactly once** into an
`OBJECT` library (`aic_objs`) from which both the shared and static libraries are
assembled — the laptop never compiles the sources twice. The version string
reported by `aic_version()` is single-sourced from `project(AIC VERSION ...)`:
CMake injects `-DAIC_VERSION_STRING="aic <version>"` onto the object library, so
the package, the SONAME, and the symbol can never disagree.

## Test

Tests are run with **CTest**. Every test has a `TIMEOUT` (fail-loud — a
non-converging point is killed, never hangs the gate), and is labelled `fast` or
`slow`.

```sh
ctest --test-dir build              # full suite (fast + slow), parallel
ctest --test-dir build -L fast      # the sub-second laptop gate (fast only)
```

via the wrapper:

```sh
make test           # full suite, parallel
make test-fast      # the sub-second gate (fast only)
```

The `fast` label is the snappy laptop gate (≈0.5 s wall on this box, 10 drivers).
The `slow` label covers the arb / MOSEK-fixture / large-`n` and measured-heavy
drivers; run the full suite (`make test`) before declaring work done, but use
`make test-fast` for the inner loop.

To build/run a single test:

```sh
ctest --test-dir build -R test_funcalc            # run one by name
cmake --build build --target test_funcalc         # just build it
```

## Benchmarks

The performance benches (`bench/bench_*.c`) are **OFF by default** so the normal
build/test path stays bench-free. They build only with `-DAIC_BUILD_BENCH=ON` and
are registered under the CTest `bench` label:

```sh
make bench                                         # configure ON, build, run
# or, by hand:
cmake -S . -B build -DAIC_BUILD_BENCH=ON
cmake --build build
ctest --test-dir build -L bench --output-on-failure
```

A bench that exits non-zero fails (a measurement that crashed is not a result).
Some benches are minutes-long; they are not part of any default gate.

## Julia package and documentation

The Julia package `AlmostIdempotentChannels.jl` calls into `libaic.so` via
`@ccall`, so build the C core first (above), then:

```sh
# Point the package at the library you just built (the default already looks
# for build/libaic.so next to the package, so this is only needed elsewhere).
julia --project=julia/AlmostIdempotentChannels.jl -e '
    using AlmostIdempotentChannels
    set_libaic_path!(joinpath(pwd(), "build", "libaic.so"))'   # then restart Julia

# Run the test suite (green with NO solver installed; ~2m09s, 219 passes).
julia --project=julia/AlmostIdempotentChannels.jl -e 'using Pkg; Pkg.test()'
```

Library discovery is via `Preferences` (no `deps/build.jl`, no `Pkg.build()`): the
default path is `build/libaic.so` next to the package, and `set_libaic_path!(p)`
overrides it (writes `LocalPreferences.toml`). The optional MOSEK extension
(`Convex` + `Mosek` + `MosekTools`) adds the exact diamond-norm value and a tight
bracket; without it the package is fully functional and the suite stays green.

To build the documentation site (executes every example against `libaic`):

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
# -> docs/build/index.html
```

See [`docs/README.md`](docs/README.md) for the documentation layout.

## Install

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build build
cmake --install build
```

**`CMAKE_INSTALL_PREFIX` must be set at CONFIGURE time**, not only at install
time: the generated `aic.pc` and `AICConfig.cmake` bake the prefix at configure,
so a late `--install --prefix` can leave the `.pc`/config pointing at the wrong
location. Set it on the `cmake -S . -B build` line as above. The install lays
down: `lib/libaic.so*` + `lib/libaic.a`, `include/aic/*.h` (internal `*_internal.h`
headers excluded), `lib/cmake/AIC/` (CMake package config + exported targets), and
`lib/pkgconfig/aic.pc`.

## Consuming the installed library

### (a) CMake — `find_package(AIC CONFIG)`

```cmake
find_package(AIC CONFIG REQUIRED)
add_executable(app main.c)
target_link_libraries(app PRIVATE AIC::aic)         # or AIC::aic_static
```

```c
#include <aic/aic.h>
```

`AIC::aic` carries the include dirs and the full dependency link line
(flint + lapacke + lapack + blas + m) as usage requirements, including the
absolute path to `libflint.so` — so the consumer needs no manual `-lflint`.
Point CMake at your install prefix with `-DCMAKE_PREFIX_PATH=<prefix>` if it is
not a system default.

### (b) pkg-config

```sh
cc app.c $(pkg-config --cflags --libs aic) -o app          # shared link
cc app.c $(pkg-config --cflags --libs --static aic) -o app # static link
```

For a **static** link you must pass `--static`: that pulls in `Libs.private`,
which spells out `-lflint -llapacke -llapack -lblas -lm` explicitly. `aic.pc`
carries `-lflint` itself because Debian's `flint.pc` omits it (see Gotchas) — so a
plain `Requires: flint` would not put `libflint` on a static link line.

## Gotchas

- **Debian `flint.pc` omits `-lflint`.** Its `Libs:` is only `-lgmp -lmpfr`
  (even with `--static`). The build resolves `libflint` via `find_library` and
  bakes the absolute path into the exported CMake targets; `aic.pc` lists
  `-lflint` explicitly in `Libs.private`. Handled — you do not need to add it.
- **`find_package(LAPACK)` does NOT provide LAPACKE.** It gives Fortran LAPACK
  only; the C-interface LAPACKE is found separately via `pkg_check_modules`. The
  final link needs both. Handled.
- **`assert()` stays live in every build type.** The library uses `assert()` for
  its fail-loud Rule-4 invariants (342 sites) and some tests assert directly, so
  the build strips `-DNDEBUG` from all optimized configs
  (Release/RelWithDebInfo/MinSizeRel). Do not re-add `-DNDEBUG` globally — it
  would silently compile out every guard.
- **Ninja "manifest dirty" / clock-skew failures.** On a box with clock skew,
  the Ninja generator can refuse to build ("manifest dirty"). The default
  **Unix Makefiles** generator (what `cmake -S . -B build` and `make` use here)
  only warns and proceeds. If you opted into Ninja and hit this, either
  `rm -rf build` and reconfigure, or drop back to the default generator.
