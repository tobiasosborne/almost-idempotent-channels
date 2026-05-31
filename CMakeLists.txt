# CMakeLists.txt — almost-idempotent-channels C core (bead aic-95g.1).
#
# FIRST CMake build for this repository. It realizes the library + test wiring
# of the polyglot-packaging design in docs/research/packaging_c_core.md (the
# "linchpin" CMake the Julia JLL and Python wheel will both consume, and the
# velocity fix that retires the hang-prone serial `make test`). This increment
# delivers ONLY the library targets and CTest wiring — NO install()/export()/
# pkg-config (a later increment), NO -fvisibility=hidden and NO AIC_EXPORT
# annotation yet (also later). The hand-written Makefile remains the ground
# truth for the build/link recipe and is left untouched this increment.
#
# House style (CLAUDE.md §"The Rules" 11/12): literate, no marketing; every
# non-obvious choice is commented to its reason. The two load-bearing choices
# here are (a) the LAPACKE link — the C-interface that find_package(LAPACK) does
# NOT provide — and (b) the single-compile OBJECT library so the laptop does not
# compile the ~85 sources twice (once for .so, once for .a).

cmake_minimum_required(VERSION 3.24)
project(AIC VERSION 0.1.0 LANGUAGES C)

# Strict C11, no GNU extensions — mirrors the Makefile's -std=c11. The library
# compiles clean under this today; keep it that way.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Default to RelWithDebInfo for single-config generators (Ninja/Make): that is
# -O2 -g, matching the Makefile's optimization+debug posture.
get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT CMAKE_BUILD_TYPE AND NOT _is_multi_config)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Build type (RelWithDebInfo => -O2 -g, like the Makefile)" FORCE)
endif()

# --- Dependencies -----------------------------------------------------------
# FLINT (>=3.0) carries the certified arb/acb path; arb is INSIDE FLINT 3.x, so
# there is no separate libarb. LAPACKE is the C interface to LAPACK that the
# double path (module latd, src/aic_latd*.c — LAPACKE_zheev/zgesvd) calls.
find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT   REQUIRED IMPORTED_TARGET flint>=3.0)
pkg_check_modules(LAPACKE REQUIRED IMPORTED_TARGET lapacke)
find_package(LAPACK REQUIRED)   # LAPACK::LAPACK — Fortran LAPACK (not LAPACKE)
find_package(BLAS   REQUIRED)   # BLAS::BLAS

# GOTCHA 1 (the LAPACKE link): `pkg-config --libs lapacke` on this box yields
# only `-llapacke` — it does NOT pull in -llapack/-lblas. find_package(LAPACK)
# gives the Fortran LAPACK but NOT the C-interface LAPACKE. So the final link
# needs ALL of: flint + lapacke + lapack + blas + m (mirrors the Makefile's
# `-lflint -llapacke -llapack -lblas -lm`).
#
# GOTCHA 2 (the missing -lflint): on Debian/Ubuntu the installed flint.pc has
#   Libs: -L${libdir} -lgmp -lmpfr
# — it OMITS -lflint itself (the .pc lists only FLINT's own deps). So
# PkgConfig::FLINT contributes the include dir + -lgmp -lmpfr but NOT libflint.
# We resolve the actual shared object by hand and link it explicitly. The
# Makefile sidesteps this by hardcoding -lflint.
find_library(FLINT_LIBRARY NAMES flint REQUIRED)
message(STATUS "AIC: using libflint at ${FLINT_LIBRARY}")

# Common usage requirements shared by the shared+static libs and the object lib.
# m is the C99 complex/real math (cabs/fabs/conj) the double path uses.
set(AIC_DEP_LIBS
    PkgConfig::FLINT
    ${FLINT_LIBRARY}
    PkgConfig::LAPACKE
    LAPACK::LAPACK
    BLAS::BLAS
    m)

# Strict warnings on the library compile (CLAUDE.md §"Build & test"). NO -Werror.
set(AIC_WARN_FLAGS -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes)

# --- Library (single-compile via OBJECT library) ----------------------------
# The OBJECT library compiles each src/*.c exactly once; both libaic.so and
# libaic.a are assembled from $<TARGET_OBJECTS:aic_objs>. This is the hard
# requirement: the laptop must not compile ~85 files twice.
#
# GLOB matches the project's existing Makefile convention ($(wildcard src/*.c));
# CONFIGURE_DEPENDS re-globs on rebuild so a new src/ file is picked up.
file(GLOB AIC_SRCS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)

add_library(aic_objs OBJECT ${AIC_SRCS})
set_target_properties(aic_objs PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(aic_objs PRIVATE ${AIC_WARN_FLAGS})
# The sources need the dep include dirs (flint/, lapacke) AND the public+
# generated include dirs to compile. PUBLIC here so they propagate to anything
# that reuses aic_objs' interface, but aic_objs is NOT itself put in an export
# set (it never leaves the build tree).
target_include_directories(aic_objs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)
target_link_libraries(aic_objs PUBLIC ${AIC_DEP_LIBS})

# Shared libaic.so and static libaic.a, both assembled from the SAME objects.
add_library(aic        SHARED $<TARGET_OBJECTS:aic_objs>)
add_library(aic_static STATIC $<TARGET_OBJECTS:aic_objs>)
set_target_properties(aic_static PROPERTIES OUTPUT_NAME aic)  # -> libaic.a

# SONAME: libaic.so.0.1.0 with SONAME libaic.so.0 (SOVERSION = major).
set_target_properties(aic PROPERTIES
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

# PUBLIC usage requirements on BOTH real libraries: include dirs (build-interface
# only this increment — no INSTALL_INTERFACE until the install() increment) and
# the full dep link line. A consumer linking aic/aic_static thus inherits the
# headers and flint+lapacke+lapack+blas+m automatically.
foreach(_lib aic aic_static)
  target_include_directories(${_lib} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)
  target_link_libraries(${_lib} PUBLIC ${AIC_DEP_LIBS})
endforeach()

# --- Export header INFRA (generated now, not yet used) ----------------------
# Generate aic_export.h so the macro (AIC_EXPORT/AIC_NO_EXPORT) exists for the
# later visibility increment. We do NOT set CMAKE_C_VISIBILITY_PRESET this
# increment, so default visibility => the library exports all symbols (like the
# current `make lib`), and tests + future Julia ccall keep full symbol access.
include(GenerateExportHeader)
generate_export_header(aic
    BASE_NAME AIC
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/include/aic/aic_export.h)

# --- Tests / CTest ----------------------------------------------------------
include(CTest)                          # defines BUILD_TESTING (default ON)
option(AIC_BUILD_TESTS "Build the AIC test suite" ON)
if(BUILD_TESTING AND AIC_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
