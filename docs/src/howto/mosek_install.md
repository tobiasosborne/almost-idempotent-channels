# Install and use the MOSEK extension

After this guide you can install the MOSEK extension, obtain a license, and use
the two entry points it unlocks: [`idempotency_defect`](../reference/api.md) (exact
diamond-norm value) and `certified_defect(Φ; tight = true)` (tight bracket).

**Assumes:** `AlmostIdempotentChannels` is already installed; see
[Quick start](../getting_started/quickstart.md).

**Tip — MOSEK is optional.** All pipeline verbs — `certified_defect` (default), `associated_algebra`,
`main_isomorphism`, `factorize` — work without MOSEK. Without the extension,
`idempotency_defect` and `certified_defect(Φ; tight = true)` throw a helpful
install hint, not a `MethodError`. Install MOSEK only when you need the exact
diamond-norm value or a tight bracket (width ~5e-13 vs ~2e-01 solver-free).

## Steps

**1. Add the three solver packages.**

```julia
]add Convex Mosek MosekTools
```

Or equivalently:

```julia
using Pkg
Pkg.add(["Convex", "Mosek", "MosekTools"])
```

**2. Obtain a MOSEK license.**

Academic licenses are free at <https://www.mosek.com/products/academic-licenses/>.
Place the downloaded `mosek.lic` file at `~/mosek/mosek.lic` (the default path
that the MOSEK Julia wrapper looks for).

**3. Load the extension.**

Loading all three packages in the same Julia session activates `AICMosekExt`
automatically:

```julia
using AlmostIdempotentChannels, LinearAlgebra
using Convex, Mosek, MosekTools
```

No explicit `using AICMosekExt` is needed; Julia's package-extension mechanism
handles activation.

**4. Call `idempotency_defect` for the exact diamond-norm value.**

The following illustrates the expected output (this block does not execute in
the default no-solver build):

```julia
# Requires: using Convex, Mosek, MosekTools
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

Φ = UCPMap(phit_kraus(2, 0.3))
idempotency_defect(Φ)
# 0.3149999999994273
# (analytic: 0.3 · 0.7 · 2 · (1 − 1/4) = 0.315; measured primal-dual gap ≤ 1.23e-11)
```

**5. Call `certified_defect(Φ; tight = true)` for a tight bracket.**

```julia
tight = certified_defect(Φ; tight = true)
width(tight)    # ~5.76e-13  (vs ~5.46e-01 solver-free)
value(tight)    # ≈ 0.315 (the SDP point estimate)
```

Measured tight bracket widths from the benchmark table:

| Channel | Tight width | Solver-free width |
|---|---|---|
| `phit(2, 0.1)` | ~1.60e-13 | ~2.34e-01 |
| `phit(2, 0.3)` | ~5.76e-13 | ~5.46e-01 |
| `paper_example(0.1)` | ~4.60e-14 | ~2.01e-01 |

**6. Note the size bound on the tight bracket.**

The tight-bracket path (`_tight_bracket_impl`) has a hard bound at ``n = 5``
(the MIN-dual SDP stalls for larger systems). For ``n > 5`` the tight-bracket
call falls back to the solver-free eig-free bracket with a `@warn`; it does not
hang. The exact diamond-norm value (`idempotency_defect`) uses only the MAX-primal
SDP and is not affected by this bound.

**7. Verify the extension is active.**

Without MOSEK, calling a gated function throws an error containing the install
hint.  With MOSEK loaded, the error is suppressed and the function returns its
result.  You can check:

```julia
# With MOSEK loaded this returns a Float64; without it throws with an install hint.
idempotency_defect(UCPMap(phit_kraus(2, 0.1)))
# 0.13499999999...
```

## See also

- [Tight brackets with MOSEK](../tutorials/mosek_tight.md) — the mathematics of the MIN-dual certifier
- [Certify the idempotency defect (solver-free)](certify_defect.md) — the always-safe alternative
- [API reference](../reference/api.md) — full signatures for `idempotency_defect` and `certified_defect`
