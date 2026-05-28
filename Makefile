# Makefile — almost-idempotent-channels C cores (beads aic-aox, "T-build").
#
# Layer-0 scaffolding: compile the src cores, auto-discover and run every
# tests/test_*.c, and reserve a bench target. Strict flags per CLAUDE.md
# §"Build & test" (mirrors ../su2-fft); -Wshadow is kept deliberately.
#
# Link policy (beads aic-5ty): ONLY -lflint. FLINT 3.0.1 bundles arb/acb here,
# so <flint/acb_mat.h> + -lflint is the whole foundation. LAPACK has no dev
# symlink on this box; the double/fast path is a deferred audition. Do NOT add
# -llapack/-lblas/-llapacke — it breaks the build.
#
# Harness include paths (beads aic-73v): tests use tests/aic_test.h, benchmarks
# use bench/aic_bench.h; both are header-only. -Itests/-Ibench are added per
# target so a test cannot accidentally see the bench header and vice versa.
# Benchmarks need clock_gettime: under -std=c11 that requires _POSIX_C_SOURCE,
# added only on bench compiles.

CC ?= cc
CFLAGS = -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes -O2 -g -std=c11 -Iinclude
LIBS = -lflint

BUILD := build
SRC := $(wildcard src/*.c)
TEST_SRC := $(wildcard tests/test_*.c)
TEST_BIN := $(patsubst tests/%.c,$(BUILD)/%,$(TEST_SRC))
BENCH_SRC := $(wildcard bench/bench_*.c)
BENCH_BIN := $(patsubst bench/%.c,$(BUILD)/%,$(BENCH_SRC))

.PHONY: all test bench clean

all: $(TEST_BIN)

$(BUILD):
	mkdir -p $(BUILD)

# Each test binary links the test driver, all src cores, and FLINT. -Itests
# exposes the header-only aic_test.h.
$(BUILD)/%: tests/%.c $(SRC) | $(BUILD)
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

clean:
	rm -rf $(BUILD)
