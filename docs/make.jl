# docs/make.jl — the Documenter.jl build for AlmostIdempotentChannels.jl
# (bead aic-exa.10 [X], the final J-series piece).
#
# The site is built from an ISOLATED docs environment (docs/Project.toml:
# Documenter + a dev'd AlmostIdempotentChannels + LinearAlgebra), NOT from the
# package's own deps — exactly like docs/plots. Documenter is therefore not a
# package dependency.
#
# `doctest = true`: every `jldoctest` block is EXECUTED and its output checked.
# We keep the exact-match `jldoctest` blocks to STABLE STRUCTURAL outputs (block
# sizes, dims, type names) and use `@example` blocks (which execute and embed the
# real output, robust to float last-digits) for everything float-heavy. The
# `@example` blocks are the live, never-rotting versions of the README narrative.
#
# Rule 13: NO `deploydocs`, NO CI. The site must just BUILD locally
# (`julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`).

using Documenter
using AlmostIdempotentChannels

# Make the package's own docstrings doctestable (so any `jldoctest` in a docstring
# runs against the package namespace).
DocMeta.setdocmeta!(
    AlmostIdempotentChannels,
    :DocTestSetup,
    :(using AlmostIdempotentChannels, LinearAlgebra);
    recursive = true,
)

makedocs(
    sitename = "AlmostIdempotentChannels.jl",
    authors = "Tobias Osborne",
    modules = [AlmostIdempotentChannels],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = nothing,
        edit_link = nothing,
        repolink = "https://github.com/tobiasosborne/almost-idempotent-channels",
        sidebar_sitename = true,
        # We are not deploying (Rule 13), so the auto-extracted Project.toml
        # version is unavailable (the docs env is isolated). Pin the inventory
        # version explicitly to silence the two cosmetic HTMLWriter warnings.
        inventory_version = "0.2.0",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "API reference" => "api.md",
    ],
    doctest = true,
    # `linkcheck` would hit the network (Rule 13: local-only); keep it off.
    # `warnonly`: do NOT fail the build on a missing-docstring cross-reference or a
    # docstring that lives only in the package (the autodocs page sweeps the
    # exported ones). Doctest FAILURES still error (they are not in this list), so
    # the doctests are a real gate.
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

# No deploydocs — Rule 13 (no CI / no auto-deploy); the site builds locally only.
