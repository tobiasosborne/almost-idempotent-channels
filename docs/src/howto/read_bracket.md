# Interpret a CertifiedBracket

After this guide you can read every field of a [`CertifiedBracket`](../reference/api.md) and
know what the bracket proves.

**Assumes:** you have a `CertifiedBracket` value in scope, either from
[`certified_defect`](../reference/api.md) or from one of the pipeline accessors.

## Steps

**1. Understand what the bracket asserts.**

A `CertifiedBracket(lo, hi)` is a proof that the true value ``x`` satisfies
``lo \le x \le hi``.  The bounds are outward-rounded double-precision images of
arb balls, so the inequality holds exactly in exact arithmetic; no tolerance is
needed on the containment side.

**2. Construct a bracket.**

```@example rb
using AlmostIdempotentChannels

b = CertifiedBracket(0.1, 0.3; label="demo")
```

**3. Read the geometric accessors.**

`minimum` and `maximum` are the proof bounds — exact for this hand-constructed
bracket:

```jldoctest rb_acc; setup = :(using AlmostIdempotentChannels; b = CertifiedBracket(0.1, 0.3; label="demo"))
julia> minimum(b)
0.1

julia> maximum(b)
0.3
```

`width` is the bracket looseness ``hi - lo`` and `midpoint` is the centre. Both
are ordinary floating-point results (so ``0.3 - 0.1`` lands at
``0.19999999999999998``, not a rounded ``0.2`` — the bracket is over `Float64`):

```@example rb
width(b), midpoint(b)
```

**4. Test containment with `∈`.**

`lo ≤ x ≤ hi` is spelled `x in b` in Julia:

```jldoctest rb_in; setup = :(using AlmostIdempotentChannels; b = CertifiedBracket(0.1, 0.3; label="demo"))
julia> 0.2 in b
true

julia> 0.5 in b
false
```

**5. Read `value` — the optional point estimate.**

`value` is `nothing` on the solver-free path.  Only `certified_defect(Φ; tight
= true)` (MOSEK required) sets a non-`nothing` value:

```jldoctest rb_val; setup = :(using AlmostIdempotentChannels; b = CertifiedBracket(0.1, 0.3; label="demo"))
julia> value(b)   # solver-free bracket carries no point estimate

julia> bv = CertifiedBracket(0.1, 0.3; value = 0.2, label="demo");

julia> value(bv)
0.2
```

(`value(b)` returns `nothing`; Documenter suppresses `nothing` output.)

**6. Know the width of the solver-free bracket.**

The eig-free bracket is loose by design (hi/lo ~ 2n); it certifies a value
rather than computing it. Measured widths: ~2.34e-01 for ``\Phi_t(2, 0.1)``
and ~5.46e-01 for ``\Phi_t(2, 0.3)``; compare to the tight MOSEK bracket
widths of ~1.60e-13 and ~5.76e-13. The reason the solver-free bracket is loose
lives in [Certified arithmetic](../explanation/certified_arithmetic.md).

**7. Handle the §B9 slack when testing against an analytic value.**

`Base.in` is exact (`lo ≤ x ≤ hi`). The arb FLOOR-rounded `lo` can sit a
hair above the analytic ``\eta`` due to outward rounding. When asserting
containment programmatically, add a small tolerance:

```julia
analytic_eta = 0.1 * 0.9 * 2 * (1 - 1/2^2)   # 0.135
tol = 1e-7
@assert minimum(η) - tol ≤ analytic_eta ≤ maximum(η) + tol
```

On an exactly-idempotent channel every bracket collapses to a neighbourhood of
0 at machine precision — see [The η = 0 oracle](../tutorials/eta0_oracle.md).

## See also

- [Certified arithmetic](../explanation/certified_arithmetic.md) — why the arb rounding guarantees the containment
- [API reference](../reference/api.md) — full signature of `CertifiedBracket`, `width`, `midpoint`, `value`
- [The five-verb pipeline](../tutorials/pipeline.md) — end-to-end walk through bracket creation and use
