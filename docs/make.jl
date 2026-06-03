# docs/make.jl — the Documenter.jl build for AlmostIdempotentChannels.jl
# (bead aic-exa.10 [X] established the site; bead aic-6k6 grew it into the full
# Diátaxis structure: Getting started / Tutorials / How-to / Explanation / Reference).
#
# The site is built from an ISOLATED docs environment (docs/Project.toml:
# Documenter + a dev'd AlmostIdempotentChannels + LinearAlgebra), NOT from the
# package's own deps — exactly like docs/plots. Documenter is therefore not a
# package dependency.
#
# `doctest = true`: every `jldoctest` block is EXECUTED and its output checked.
# We keep the exact-match `jldoctest` blocks to STABLE STRUCTURAL outputs (block
# sizes, dims, type names, booleans) and use `@example` blocks (which execute and
# embed the real output, robust to float last-digits) for everything float-heavy.
# The `@example` blocks are the live, never-rotting versions of the narrative.
#
# Math is PLAIN KaTeX (Documenter's default mathengine) — no custom macros — so
# every page renders without a macro-config dependency (see docs/research/docsite/
# STYLE.md §4 for the paper-macro → KaTeX translation the pages follow).
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
        collapselevel = 1,
        # We are not deploying (Rule 13), so the auto-extracted Project.toml
        # version is unavailable (the docs env is isolated). Pin the inventory
        # version explicitly to silence the two cosmetic HTMLWriter warnings.
        inventory_version = "0.2.0",
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => [
            "Installation" => "getting_started/install.md",
            "Quick start"  => "getting_started/quickstart.md",
        ],
        "Tutorials" => [
            "The five-verb pipeline"    => "tutorials/pipeline.md",
            "The η = 0 oracle"          => "tutorials/eta0_oracle.md",
            "Multi-block channels"      => "tutorials/multiblock.md",
            "Tight brackets with MOSEK" => "tutorials/mosek_tight.md",
        ],
        "How-to guides" => [
            "Certify the defect (solver-free)" => "howto/certify_defect.md",
            "Interpret a CertifiedBracket"     => "howto/read_bracket.md",
            "Extract block structure"          => "howto/block_structure.md",
            "Factorize a channel"              => "howto/factorize.md",
            "Stay in factorize's domain"       => "howto/domain_check.md",
            "Check trace-preservation"         => "howto/iscptp.md",
            "Build a UCPMap from Kraus"        => "howto/build_ucpmap.md",
            "Install and use MOSEK"            => "howto/mosek_install.md",
        ],
        "Explanation" => [
            "The mathematics"                   => "explanation/math_story.md",
            "Certified arithmetic"              => "explanation/certified_arithmetic.md",
            "Dimension-independence"            => "explanation/dim_independence.md",
            "Architecture: C + FLINT + Julia"   => "explanation/architecture.md",
            "The constructive mandate"          => "explanation/constructive.md",
            "Design decisions and known limits" => "explanation/design_limits.md",
        ],
        "Reference" => [
            "API"               => "reference/api.md",
            "Notation glossary" => "reference/notation.md",
            "Bibliography"      => "reference/bibliography.md",
        ],
    ],
    doctest = true,
    # `linkcheck` would hit the network (Rule 13: local-only); keep it off.
    # `warnonly`: do NOT fail the build on a missing-docstring cross-reference or a
    # docstring that lives only in the package. Doctest FAILURES still error (they
    # are not in this list), so the doctests remain a real gate.
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

# No deploydocs — Rule 13 (no CI / no auto-deploy); the site builds locally only.
