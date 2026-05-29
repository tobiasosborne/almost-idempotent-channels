# Makefile — almost-idempotent-channels C cores (beads aic-aox, "T-build").
#
# Layer-0 scaffolding: compile the src cores, auto-discover and run every
# tests/test_*.c, and reserve a bench target. Strict flags per CLAUDE.md
# §"Build & test" (mirrors ../su2-fft); -Wshadow is kept deliberately.
#
# Link policy (beads aic-5ty, updated aic-w4o.5): -lflint for the CERTIFIED
# arb/acb path (FLINT 3.0.1 bundles arb/acb, so <flint/acb_mat.h> + -lflint is
# the foundation) PLUS -llapacke -llapack -lblas for the FAST, UNCERTIFIED
# double path (module latd: include/aic_latd.h, src/aic_latd*.c). LAPACKE is now
# present and links/runs on this box (verified). The LAPACK libs are linked into
# EVERY target — harmless for targets that do not call them — so the
# auto-discovered tests/benches need no per-target link tweaks. -lm is added for
# the C99 complex/real math (cabs/fabs/conj) the double-path tests use.
#
# Harness include paths (beads aic-73v): tests use tests/aic_test.h, benchmarks
# use bench/aic_bench.h; both are header-only. -Itests/-Ibench are added per
# target so a test cannot accidentally see the bench header and vice versa.
# Benchmarks need clock_gettime: under -std=c11 that requires _POSIX_C_SOURCE,
# added only on bench compiles.

CC ?= cc
CFLAGS = -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes -O2 -g -std=c11 -Iinclude
LIBS = -lflint -llapacke -llapack -lblas -lm

BUILD := build
SRC := $(wildcard src/*.c)
TEST_SRC := $(wildcard tests/test_*.c)
TEST_BIN := $(patsubst tests/%.c,$(BUILD)/%,$(TEST_SRC))
BENCH_SRC := $(wildcard bench/bench_*.c)
BENCH_BIN := $(patsubst bench/%.c,$(BUILD)/%,$(BENCH_SRC))
# Generated/committed fixture includes (e.g. tests/fixtures_d24.inc.h, the
# aic-d24 golden master from tools/gen_fixtures_d24.jl). Made a prerequisite of
# every test binary so regenerating a fixture forces the dependent test to
# rebuild — otherwise a stale binary would test against old golden values.
TEST_INC := $(wildcard tests/*.inc.h)

.PHONY: all test bench clean fixtures lib

all: $(TEST_BIN)

$(BUILD):
	mkdir -p $(BUILD)

# Each test binary links the test driver, all src cores, and FLINT. -Itests
# exposes the header-only aic_test.h and the generated tests/*.inc.h fixtures
# (which are prerequisites so a regenerated fixture rebuilds its test).
$(BUILD)/%: tests/%.c $(SRC) $(TEST_INC) | $(BUILD)
	$(CC) $(CFLAGS) -Itests $< $(SRC) $(LIBS) -o $@

# Each bench binary links the bench driver, all src cores, and FLINT. -Ibench
# exposes aic_bench.h; _POSIX_C_SOURCE enables clock_gettime under -std=c11.
$(BUILD)/%: bench/%.c $(SRC) | $(BUILD)
	$(CC) $(CFLAGS) -Ibench -D_POSIX_C_SOURCE=200809L $< $(SRC) $(LIBS) -o $@

# Build then run every test; report per-test pass/fail and fail the target if
# any binary exits non-zero (fail-loud, CLAUDE.md Rule 4).
test: $(TEST_BIN)
	@fail=0; \
	for t in $(TEST_BIN); do \
		echo "=== running $$t ==="; \
		if ./$$t; then \
			echo "PASS: $$t"; \
		else \
			echo "FAIL: $$t (exit $$?)"; \
			fail=1; \
		fi; \
	done; \
	if [ $$fail -ne 0 ]; then echo "make test: FAILED"; exit 1; fi; \
	echo "make test: all tests passed"

# Build then run every bench/bench_*.c; each prints its own greppable BENCH
# line(s) (beads aic-73v). A benchmark that exits non-zero fails the target
# (fail-loud, CLAUDE.md Rule 4) — a measurement that crashed is not a result.
bench: $(BENCH_BIN)
	@if [ -z "$(BENCH_BIN)" ]; then echo "no benchmarks found"; exit 0; fi; \
	fail=0; \
	for b in $(BENCH_BIN); do \
		echo "=== running $$b ==="; \
		if ! ./$$b; then echo "FAIL: $$b (exit $$?)"; fail=1; fi; \
	done; \
	if [ $$fail -ne 0 ]; then echo "make bench: FAILED"; exit 1; fi; \
	echo "make bench: all benchmarks ran"

# Regenerate the committed golden-master fixtures (bead aic-d24) via the serial
# Julia + MOSEK generator. OPTIONAL and NOT a prerequisite of `make test` (the
# generated .inc.h is committed); run it by hand when the corpus changes. Needs
# the julia/env environment instantiated and a MOSEK license.
fixtures:
	julia --project=julia/env tools/gen_fixtures_d24.jl

# Shared library for the Julia ccall surface (bead aic-m24, increment 2a). The
# flat-double shim (src/aic_ucp_shim.c) is the ccall entry point; $(SRC) already
# globs all src/*.c so it is auto-included. -fPIC is kept LOCAL to this recipe
# (NOT added to the global CFLAGS) to minimize blast radius on the test/bench
# builds. The Julia package (increment 2b) dlopen's build/libaic.so.
lib: $(BUILD)/libaic.so
$(BUILD)/libaic.so: $(SRC) | $(BUILD)
	$(CC) $(CFLAGS) -fPIC -shared $(SRC) $(LIBS) -o $@

clean:
	rm -rf $(BUILD)
