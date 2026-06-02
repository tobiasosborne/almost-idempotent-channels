# libaic.jl — libaic.so discovery (Preferences), the dlopen'd handle, the
# dlsym'd symbol table for every C entry point, and __init__.
#
# This file OWNS the native-library lifecycle (moved out of the main module so
# the module body stays solver-free and declarative). Design doc §6.1.
#
# LIBPATH DISCOVERY (Preferences, replacing the old deps/deps.jl):
#   default = build/libaic.so relative to the repo root. This package lives
#   THREE directory levels below the repo root (julia/AlmostIdempotentChannels.jl/src),
#   so the repo root is joinpath(@__DIR__, "..", "..", ".."). A user who has run
#   cmake once needs no Pkg.build; to point elsewhere call set_libaic_path!(p)
#   (writes a LocalPreferences.toml entry; restart Julia to apply).
#
# SYMBOL TABLE: one Ref{Ptr{Cvoid}} per C entry point. __init__ dlsym's each into
# its Ref; the ccall wrappers (in the main module / later beads) read the Ref and
# use the dlsym-handle 4-arg ccall form (research §1e: @ccall cannot take a raw
# Ptr{Cvoid}, and the dlsym handle is what lets us combine the load-bearing RTLD
# flags with a type-safe ccall). Entry points wrapped:
#   EXISTING (J0): aic_ucp_choi_diff_d, aic_cbnorm_eigfree_d, aic_cbnorm_certify_d
#   NEW (frozen ABI, aic-exa.3 commit cf7ebc8; ccalls land in bead J3 aic-exa.6):
#     aic_assoc_summary_d, aic_main_iso_summary_d,
#     aic_factorize_artifacts_sizes_d, aic_factorize_artifacts_d

using Libdl
using Preferences

# ----- libpath discovery (Preferences) -----

const _DEFAULT_LIBAIC = let
    repo = abspath(joinpath(@__DIR__, "..", "..", ".."))   # pkg is 3 levels deep
    joinpath(repo, "build", "libaic.so")
end

# @load_preference is read at PRECOMPILE time; the value is baked into the
# precompiled image. set_libaic_path! invalidates the cache so the next load
# (after restart) picks up the new path.
const LIBAIC_PATH = @load_preference("libaic_path", _DEFAULT_LIBAIC)

"""
    set_libaic_path!(p) -> Nothing

Persist an absolute path to `libaic.so` as a package Preference
(LocalPreferences.toml). Use when the default `build/libaic.so` (relative to the
repo root) does not resolve — e.g. an out-of-tree build directory or an installed
prefix. Restart Julia to apply (the path is baked into the precompiled image).
"""
function set_libaic_path!(p)
    @set_preferences!("libaic_path" => string(p))
    @info "libaic_path set to $(string(p)); restart Julia to apply."
    return nothing
end

"""
    libaic_path() -> String

Absolute path to the dynamically-loaded `libaic.so` (the Preferences value, or
the default `build/libaic.so` relative to the repo root). For debugging linkage.
"""
libaic_path() = LIBAIC_PATH

# ----- the dlopen'd handle + dlsym'd symbol table -----

const _LIBAIC_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)

# Symbol-name => Ref{Ptr{Cvoid}} for every C entry point. __init__ dlsym's each.
# The Refs are module-private; the ccall wrappers index this table by name.
const _SYM_CHOI_DIFF_D               = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_CBNORM_EIGFREE_D          = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_CBNORM_CERTIFY_D          = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_ASSOC_SUMMARY_D           = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_MAIN_ISO_SUMMARY_D        = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_FACTORIZE_ARTIFACTS_SIZES_D = Ref{Ptr{Cvoid}}(C_NULL)
const _SYM_FACTORIZE_ARTIFACTS_D     = Ref{Ptr{Cvoid}}(C_NULL)

# (symbol :name, the Ref to fill) — iterated by __init__.
const _ALL_SYMS = (
    (:aic_ucp_choi_diff_d,             _SYM_CHOI_DIFF_D),
    (:aic_cbnorm_eigfree_d,            _SYM_CBNORM_EIGFREE_D),
    (:aic_cbnorm_certify_d,            _SYM_CBNORM_CERTIFY_D),
    (:aic_assoc_summary_d,             _SYM_ASSOC_SUMMARY_D),
    (:aic_main_iso_summary_d,          _SYM_MAIN_ISO_SUMMARY_D),
    (:aic_factorize_artifacts_sizes_d, _SYM_FACTORIZE_ARTIFACTS_SIZES_D),
    (:aic_factorize_artifacts_d,       _SYM_FACTORIZE_ARTIFACTS_D),
)

# ----- the MOSEK error hint (belt-and-suspenders MethodError hint) -----

# Registered in __init__ as a Base.Experimental.register_error_hint so a
# MethodError on the solver-gated stubs (sdp_stubs.jl) prints the install hint.
# The string body is the const _MOSEK_HINT defined in sdp_stubs.jl (included
# after this file), so reference it lazily inside the hint closure.
function _mosek_error_hint(io, exc, argtypes, kwargs)
    fn = exc.f
    gated = (idempotency_defect, _diamond_value_impl,
             diamond_norm_watrous, diamond_norm_watrous_primal, diamond_norm_dual)
    if any(g -> fn === g, gated)
        print(io, "\n", _MOSEK_HINT)
    end
    return nothing
end

# ----- native-library lifecycle -----

function __init__()
    if !isfile(LIBAIC_PATH)
        error("AlmostIdempotentChannels: libaic.so not found at \"$LIBAIC_PATH\".\n" *
              "Build it:  cmake -S . -B build && cmake --build build   (in the repo root)\n" *
              "or point at an existing build:  " *
              "AlmostIdempotentChannels.set_libaic_path!(\"/abs/path/libaic.so\")  (then restart Julia)")
    end
    # Pre-load the SYSTEM libgmp.so.10 with RTLD_GLOBAL before libaic/libflint
    # are mapped. This works around the Julia bundled-libgmp ABI mismatch: on
    # this system RTLD_DEEPBIND on libaic alone was insufficient because libflint
    # resolves __gmpn_* symbols against the already-loaded Julia-bundled libgmp;
    # pre-loading the system libgmp keeps its definitions in the global symbol
    # table first (the exact issue ../su2-fft handles, bead su2fft-e5z).
    # KEEP VERBATIM — load-bearing on this machine.
    try
        Libdl.dlopen("/lib/x86_64-linux-gnu/libgmp.so.10",
                     Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL)
    catch
        try
            Libdl.dlopen("libgmp.so.10", Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL)
        catch
            # silently ignore on systems without a discoverable system libgmp
        end
    end
    flags = Libdl.RTLD_NOW | Libdl.RTLD_GLOBAL | Libdl.RTLD_DEEPBIND
    _LIBAIC_HANDLE[] = Libdl.dlopen(LIBAIC_PATH, flags)
    for (sym, ref) in _ALL_SYMS
        ref[] = Libdl.dlsym(_LIBAIC_HANDLE[], sym)
    end
    Base.Experimental.register_error_hint(_mosek_error_hint, MethodError)
    return nothing
end
