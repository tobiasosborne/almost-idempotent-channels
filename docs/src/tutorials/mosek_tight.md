# Tight brackets with MOSEK

After this page you will know how to get the exact diamond-norm value and a tight
rigorous bracket (width ``\sim10^{-13}``) instead of the solver-free bracket
(width ``\sim 2n``) — when MOSEK is installed.

The examples on this page cannot execute in the default solver-free build (pitfall
P4 in the style guide). The code blocks are illustrative; the outputs shown are
from passing tests (`test_sdp.jl`, `README.md`).

## What changes with MOSEK

The solver-free path certifies a value via the eig-free arb certifier. The bracket
width is ``\sim 2n`` by design — it proves containment, not equality. MOSEK adds
two things:

1. The **exact diamond-norm value** via the Watrous SDP (`idempotency_defect`).
2. A **tight rigorous bracket** fed by the SDP primal/dual feasible points
   (`certified_defect(Φ; tight=true)`), with width ``\sim10^{-13}``.

!!! tip "With MOSEK installed"
    Activate the `AICMosekExt` extension by loading `Convex`, `Mosek`, and
    `MosekTools` before `AlmostIdempotentChannels`:

    ```julia
    using AlmostIdempotentChannels, LinearAlgebra
    using Convex, Mosek, MosekTools   # activates AICMosekExt

    e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
    function phit_kraus(d, t)
        K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
        for i in 1:d, j in 1:d
            push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
        end
        return K
    end

    Φ = UCPMap(phit_kraus(2, 0.3))

    # Exact diamond-norm value (Watrous SDP):
    idempotency_defect(Φ)
    # output: 0.3149999999994273
    # = 0.3 · 0.7 · 2 · (1 − 1/4) = 0.315 to solver precision

    # Tight rigorous bracket:
    tight = certified_defect(Φ; tight = true)
    width(tight)
    # output: ~5.76e-13
    value(tight)
    # output: ≈ 0.315  (the SDP point estimate)
    ```

    Without MOSEK both functions throw a helpful install hint rather than a
    `MethodError`.

## Solver-free vs MOSEK bracket widths

The table below compares the eig-free (solver-free) bracket width with the tight
MOSEK bracket width for three fixtures. Numbers from `README.md` benchmark table:

| Fixture | Eig-free width | MOSEK tight width |
|---|---|---|
| `phit(2, 0.3)` | 5.46e-01 | 5.76e-13 |
| `phit(2, 0.1)` | 2.34e-01 | 1.60e-13 |
| `paper_example(0.1)` | 2.01e-01 | 4.60e-14 |

The solver-free bracket is loose by design (`hi/lo ~ 2n`). The MOSEK bracket is
``\sim10^{13}`` times narrower and carries a point estimate in `value(tight)`.

## Strong duality cross-check

The Watrous SDP for diamond norm is self-dual. The maximum primal-dual gap across
all test fixtures is 1.23e-11 (`test_sdp.jl`). This is the strong-duality
cross-check: if the gap is large, the SDP solve failed and the tight bracket
is unreliable.

!!! note
    `phit(2, 0.3)` is outside `factorize`'s domain (``\rho \approx 0.21 > 0.10``)
    but is safe for `certified_defect` at any ``\eta``. The MOSEK tight bracket
    and `idempotency_defect` work at any ``\eta``.

## Where to go next

- [Install and use the MOSEK extension](@ref) — step-by-step install instructions.
- [Certified arithmetic](@ref) — the four-rung cross-check ladder, where the MOSEK
  strong-duality check is rung 4.
