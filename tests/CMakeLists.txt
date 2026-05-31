# tests/CMakeLists.txt — CTest wiring for the AIC C core (bead aic-95g.1).
#
# Replaces the Makefile's serial `make test` glob-and-exec loop. Each
# tests/test_*.c becomes one executable + one add_test; CTest runs them in
# parallel and — crucially — never auto-executes a stray probe (the hang-probe
# class, bead aic-xo0). Every test gets a TIMEOUT so a non-converging point
# fails loud (CLAUDE.md Rule 4) instead of hanging the gate forever.
#
# Tests link the STATIC library (libaic.a), not the shared one: this keeps FULL
# internal-symbol access once the later increment turns on -fvisibility=hidden
# (a hidden symbol is unreachable across a .so boundary but always reachable
# when statically archived into the test binary). The static lib already
# carries the flint+lapacke+lapack+blas+m link line and the include dirs as
# PUBLIC usage requirements, so tests inherit them.

set(AIC_WARN_FLAGS -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes)

# --- Shared non-driver test helpers -----------------------------------------
# tests/aic_*.c are NOT test drivers (no main()) — they are the shared
# adversarial-instance corpus (aic_adversarial.c, aic_adversarial_nla.c). The
# Makefile compiles them into every test binary; here we archive them once into
# a static support lib that every test links. The glob deliberately matches
# `aic_*.c` (helpers), never `test_*.c` (drivers).
file(GLOB AIC_TEST_HELPERS CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/aic_*.c)
add_library(aic_testsupport STATIC ${AIC_TEST_HELPERS})
target_compile_options(aic_testsupport PRIVATE ${AIC_WARN_FLAGS})
# -Itests so the helpers see aic_test.h / aic_adversarial.h; aic_static carries
# the project include dirs + the dep link line.
target_include_directories(aic_testsupport PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(aic_testsupport PUBLIC aic_static)

# --- Per-test executables + add_test ----------------------------------------
# `slow` set: arb / MOSEK-fixture / large-n drivers, PLUS the measured-heavy
# drivers whose wall time would make `ctest -L fast` a poor laptop gate. The
# fast set is everything else — all sub-second on this box (measured via CTest
# CostData: the 10 fast drivers run in <1 s each; the five reclassified here
# cost test_errreduce 115 s, test_dhom 30 s, test_projection/test_corner 15 s,
# test_ecstar 10 s). So `ctest -L fast` is the snappy gate; `ctest` (no -L) is
# the full suite. Refinable by the orchestrator.
set(AIC_SLOW_TESTS
    test_cstar_build test_cstar_extension test_cstar_merging test_cstar_merge
    test_cstar test_certify test_certify_teeth test_cbnorm test_opspace
    test_opspace_o2 test_opspace_choi test_factorize test_assoc test_assoc2
    test_adversarial test_ucp_d24
    # measured-heavy (not MOSEK/arb-fixture, just long): keep them out of -L fast
    test_errreduce test_dhom test_projection test_corner test_ecstar)

file(GLOB AIC_TEST_SRCS CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/test_*.c)

foreach(_src ${AIC_TEST_SRCS})
  get_filename_component(_name ${_src} NAME_WE)
  add_executable(${_name} ${_src})
  target_compile_options(${_name} PRIVATE ${AIC_WARN_FLAGS})
  # -Itests so each test sees aic_test.h and the committed *.inc.h fixtures.
  target_include_directories(${_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  # Link the STATIC lib (internal symbols) + the helper corpus.
  target_link_libraries(${_name} PRIVATE aic_static aic_testsupport)

  add_test(NAME ${_name} COMMAND ${_name})

  # TIMEOUT on EVERY test — fail-loud kill of the hang class (bead aic-xo0).
  # slow gets 600 s, not 300: the heaviest arb drivers run ~200 s SOLO and the
  # 300 s cap produced FALSE timeouts under -j contention (hostile review
  # aic-95g.1, test_assoc2 timed out at 300 s under -j6 but passes in 203 s
  # solo). The timeout exists to catch HANGS (infinite), so a generous 600 s
  # still fails loud without flaking under load. fast stays tight at 120 s (all
  # those drivers are sub-second).
  if(${_name} IN_LIST AIC_SLOW_TESTS)
    set_tests_properties(${_name} PROPERTIES TIMEOUT 600 LABELS "slow")
  else()
    set_tests_properties(${_name} PROPERTIES TIMEOUT 120 LABELS "fast")
  endif()
endforeach()
