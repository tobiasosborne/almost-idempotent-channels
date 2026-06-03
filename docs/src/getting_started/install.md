# Installation

After this page you will have a built `libaic.so`, the Julia package pointed at
it, and the solver-free test suite green (219 tests, no MOSEK required).

## Step 1 — C-core dependencies

The C/arb core requires a C11 compiler, CMake >= 3.24, FLINT >= 3.0 (which
brings arb), and LAPACK/LAPACKE. Tested on Debian/Ubuntu Linux with Julia 1.10+.

```sh
sudo apt install build-essential cmake libflint-dev liblapacke-dev liblapack-dev libblas-dev
```

On other distributions, install the equivalent packages that provide `cmake`,
`libflint` (FLINT >= 3.0 / arb), and a LAPACK implementation with its LAPACKE
C interface.

## Step 2 — Build the C/arb core

From the root of the repository:

```sh
cmake -S . -B build && cmake --build build
```

This produces `build/libaic.so` (and `build/libaic.a`). The build takes under
a minute on a laptop; the CMake output lists each compiled source file.

## Step 3 — Point the package at the library

Library discovery uses `Preferences` — there is no `deps/build.jl` and no
`Pkg.build()` step.

The default looks for `build/libaic.so` relative to the repository root, which
is where Step 2 put it. If you are working from that directory, no extra
configuration is needed.

If you built to a different location, or installed `libaic` system-wide, tell
the package where to find it:

```julia
using AlmostIdempotentChannels
set_libaic_path!(joinpath(pwd(), "build", "libaic.so"))
# then restart Julia
```

`set_libaic_path!` writes `LocalPreferences.toml` next to the package's
`Project.toml`. Restart Julia after calling it so the new preference takes
effect.

## Step 4 — Run the solver-free test suite

```sh
julia --project=julia/AlmostIdempotentChannels.jl -e 'using Pkg; Pkg.test()'
```

A passing run ends with:

```
Test Summary:                      | Pass  Total  Time
AlmostIdempotentChannels           |  219    219  ~2m09s
```

All 219 tests pass with no solver installed.

## Optional: the MOSEK extension

A MOSEK license unlocks two additional verbs: `idempotency_defect` (the exact
diamond-norm value via the Watrous SDP) and `certified_defect(Φ; tight=true)`
(a tight bracket fed the SDP feasible points, ~10^12x narrower than the
solver-free bracket). See [Install and use the MOSEK extension](../howto/mosek_install.md) for the install steps
and [Tight brackets with MOSEK](../tutorials/mosek_tight.md) for usage examples.

The solver-free pipeline — `certified_defect`, `associated_algebra`,
`main_isomorphism`, and `factorize` — works without a solver at any step.

---

**Next:** [Quick start](quickstart.md) for a first certified run, or [The five-verb pipeline](../tutorials/pipeline.md)
for the full chain.
