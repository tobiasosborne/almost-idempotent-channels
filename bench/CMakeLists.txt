# bench/CMakeLists.txt — CTest wiring for the AIC performance benches
# (bead aic-95g.1, Inc4). Built ONLY when -DAIC_BUILD_BENCH=ON (the parent
# CMakeLists guards the add_subdirectory), so the normal laptop build/test path
# never pays for them.
#
# Each bench/bench_*.c is one executable + one add_test, labelled `bench` so
# `ctest -L bench` runs the audition slate (CLAUDE.md Law 4). Mirrors
# tests/CMakeLists.txt, but simpler: the benches use no shared aic_*.c corpus
# (there is none under bench/) and need no fixture .inc.h, so there is no
# support library — each bench links the STATIC libaic.a directly. aic_static
# already carries the flint+lapacke+lapack+blas+m link line and the project
# include dirs as PUBLIC usage requirements, so the benches inherit them.
#
# Two bench-specific compile requirements vs the tests:
#   -Ibench                      so each bench sees the header-only aic_bench.h
#                                (the AIC_BENCH timing macro). Kept LOCAL to the
#                                bench targets so a test can never see it.
#   -D_POSIX_C_SOURCE=200809L    aic_bench.h calls clock_gettime(CLOCK_MONOTONIC),
#                                which -std=c11 (no GNU extensions) does not
#                                expose without the POSIX feature-test macro.

set(AIC_WARN_FLAGS -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes)

file(GLOB AIC_BENCH_SRCS CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/bench_*.c)

foreach(_src ${AIC_BENCH_SRCS})
  get_filename_component(_name ${_src} NAME_WE)
  add_executable(${_name} ${_src})
  target_compile_options(${_name} PRIVATE ${AIC_WARN_FLAGS})
  target_compile_definitions(${_name} PRIVATE _POSIX_C_SOURCE=200809L)
  target_include_directories(${_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${_name} PRIVATE aic_static)

  add_test(NAME ${_name} COMMAND ${_name})
  # A bench that exits non-zero fails (fail-loud). TIMEOUT guards a runaway loop
  # the same way the tests are guarded (bead aic-xo0).
  set_tests_properties(${_name} PROPERTIES TIMEOUT 300 LABELS "bench")
endforeach()
