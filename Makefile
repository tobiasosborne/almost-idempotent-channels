# Makefile — thin convenience wrapper over the CMake build (bead aic-95g.1).
#
# CMake (CMakeLists.txt) is the CANONICAL build; this Makefile only forwards to
# cmake/ctest so the documented `make` / `make test` muscle-memory still works.
# It is NOT the build recipe any more — see CMakeLists.txt and BUILDING.md.
#
# Old -> new command mapping:
#   old `make` / `make all`     -> configure + build the libraries (libaic.so + .a)
#   old `make test`             -> ctest, full suite, parallel, fail-loud TIMEOUT
#   new `make test-fast`        -> ctest -L fast (sub-second laptop gate)
#   old `make build/test_X`     -> `ctest --test-dir build -R test_X`
#                                  or `cmake --build build --target test_X`
#   old `make bench`            -> build+run benches (behind -DAIC_BUILD_BENCH=ON)
#   old `make lib`              -> `make` (CMake always builds libaic.so AND .a)
#   old `make fixtures`         -> unchanged (Julia + MOSEK golden-master generator)
#   old `make clean`            -> remove the build dir
#
# Generator: the cmake DEFAULT (Unix Makefiles) — we deliberately do NOT force
# Ninja. This box has shown a clock-skew "manifest dirty" Ninja failure; the
# Make generator is robust against it (BUILDING.md "Gotchas").

BUILD ?= build
JOBS  ?= $(shell nproc 2>/dev/null || echo 2)

.PHONY: all configure test test-fast bench fixtures install clean

configure:
	cmake -S . -B $(BUILD)

all: configure
	cmake --build $(BUILD) -j $(JOBS)

# Full suite: every test, parallel, print the log of any failure. Each test has
# a CTest TIMEOUT (fail-loud, no hangs).
test: all
	ctest --test-dir $(BUILD) --output-on-failure -j $(JOBS)

# Sub-second laptop gate: only the `fast`-labelled tests.
test-fast: all
	ctest --test-dir $(BUILD) -L fast --output-on-failure -j $(JOBS)

# Benches are OFF in the default configure; reconfigure with the option ON,
# build, then run the `bench`-labelled CTest entries.
bench:
	cmake -S . -B $(BUILD) -DAIC_BUILD_BENCH=ON
	cmake --build $(BUILD) -j $(JOBS)
	ctest --test-dir $(BUILD) -L bench --output-on-failure

# Install to PREFIX (or the configure-time CMAKE_INSTALL_PREFIX if unset). NOTE:
# for aic.pc / AICConfig to bake the right prefix, set it at CONFIGURE time:
#   cmake -S . -B build -DCMAKE_INSTALL_PREFIX=<prefix>   (see BUILDING.md).
install: all
	cmake --install $(BUILD) $(if $(PREFIX),--prefix $(PREFIX),)

# Regenerate the committed golden-master fixtures (bead aic-d24). Unchanged from
# the old Makefile: serial Julia + MOSEK generator, run by hand when the corpus
# changes. Needs julia/env instantiated and a MOSEK license.
fixtures:
	julia --project=julia/env tools/gen_fixtures_d24.jl

clean:
	rm -rf $(BUILD)
