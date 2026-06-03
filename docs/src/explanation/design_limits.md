# Design decisions and known limits

After this page you will understand the three known limits of the package, the
root cause of each, and the design choices that produced them.

This page leads with the constraints, in the spirit of SQLite's "Quirks" page:
the things that surprise a reader who expected exact behaviour. Each limit is
real, each has a measured number, and each is a consequence of a deliberate
choice rather than an unfinished corner. To see *why* the underlying mathematics
forces these tradeoffs, see [The mathematics](@ref); for the certificate
machinery, see [Certified arithmetic](@ref).

## Known limits

### `factorize` has a tighter domain than the rest of the pipeline

`factorize` builds the encode map ``\Upsilon`` via a power-series functional
calculus whose convergence domain is much smaller than the `prop_P` basin
``\rho(\Phi^2-\Phi) < 1/4``. The verb pre-checks a conservative ``\rho < 0.10``
in Julia and throws a clean `ArgumentError` — it does not abort the session.
`associated_algebra` and `main_isomorphism` keep the wider ``\rho < 1/4`` domain;
`certified_defect` is always safe at any ``\eta``.

**Measured boundary.** `phit_kraus(2, 0.3)` has ``\rho \approx 0.21`` and throws;
`phit_kraus(2, 0.1)` has ``\rho \approx 0.09`` and is inside; `mixconj_kraus(4,
2, 0.02)` has ``\rho \approx 0.04`` and is inside (`test_core.jl`, `fixtures.jl`).

**Root cause.** ``\Upsilon`` is assembled from a Taylor/power-series expansion
``\texttt{xpow}`` of the regularising functional calculus, and that series
converges on a strictly smaller disc than the spectral basin the projection
construction `prop_P` lives in. The basin guard ``\rho(\Phi^2-\Phi) < 1/4`` is the
right gate for the algebra-side verbs but is too generous for the series-side
encode map; the conservative ``\rho < 0.10`` pre-check is the documented fix (bug
`aic-exa.13`, `FINDINGS §D4`/`§C13`). The error message names the convergence
domain so the reader is pointed back to `certified_defect`, which has no domain
restriction.

!!! warning "`factorize` domain"
    Do not call `factorize` on a fixture with ``\rho(\Phi^2-\Phi) \ge 0.10``; it
    throws `ArgumentError`. Safe inputs: any ``\eta = 0`` oracle,
    `phit_kraus(2, 0.1)`, `bce_conj_kraus(4, 2, 0.02)`, `mixconj_kraus(4, 2,
    0.02)`. For the staying-in-domain recipe see
    [Stay in factorize's domain](@ref).

### `decode` is only ``O(\eta)``-trace-preserving for ``\eta > 0``

The `decode` channel is only ``O(\eta)``-trace-preserving when ``\eta > 0``.
`iscptp(decode(F))` is `false` at the default `atol = 1e-9` and `true` at `atol =
1e-4`. The rigorous certificate is the cb-norm round-trip bracket, not `iscptp`
at machine tolerance. The encode channel and the ``\eta = 0`` oracle decode are
trace-preserving to machine precision.

**Measured number.** The clipped mass on multi-block oblique fixtures is
``\approx 3.7\mathrm{e}{-6}`` (`README.md` Known limits); the trace-preservation
defect of `decode` on `bce_conj_kraus(4, 2, 0.02)` is ``1.9\mathrm{e}{-6}``, and
``2.2\mathrm{e}{-5}`` at ``t = 0.05`` (`fixtures.jl`).

**Root cause.** The decode map ``\Delta'`` is only ``O(\eta)``-completely-positive,
not exactly CP. The paper's manifest-PSD argument
(`approximate_algebras.tex:2786-2796`) needs ``\tilde\Delta = v`` to be an *exact*
unital homomorphism, but ``v`` is only an extended ``O(\eta)``-isomorphism — the
Choi–Effros star is not the ordinary product — so the per-block Choi carries a
small ``O(\eta^2)`` negative eigenvalue once there are two or more blocks and
``\eta > 0``. The package projects onto the CP cone with a tolerance-aware
Choi-to-Kraus extraction that clips the ``O(\eta^2)`` defect but still aborts on a
genuine ``O(1)`` or ``O(\eta)`` negative eigenvalue. The clip stays within the
``\|\Delta - \tilde\Delta\|_{\mathrm{cb}} \le O(\eta)`` specification, so the
round-trip bracket remains the honest certificate (`FINDINGS §C14`). The
``3.7\mathrm{e}{-6}`` is the mass the clip removes, not an error in the bound.

!!! warning "`decode` trace-preservation tolerance"
    On an ``\eta > 0`` multi-block fixture, `iscptp(decode(F))` returns `false`
    at the default `atol = 1e-9`. Call `iscptp(decode(F); atol = 1e-4)` and read
    the round-trip bracket as the certificate. For the full recipe see
    [Check trace-preservation of a channel](@ref).

### The eig-free bracket is loose by design

The eig-free bracket is loose by design (``\mathrm{hi}/\mathrm{lo} \sim 2n``); it
certifies a value rather than computing it. For a tight bracket or the exact
value, use the MOSEK extension.

**Measured widths.** At ``\mathrm{prec} = 128``, `phit(2, 0.3)` brackets to width
``5.46\mathrm{e}{-1}`` eig-free versus ``5.76\mathrm{e}{-13}`` tight; `phit(2,
0.1)` to ``2.34\mathrm{e}{-1}`` versus ``1.60\mathrm{e}{-13}``; `paper(0.1)` to
``2.01\mathrm{e}{-1}`` versus ``4.60\mathrm{e}{-14}`` (`README.md` benchmark
table). The tight rung is roughly ``10^{12}\times`` narrower.

**Root cause.** The cheapest, always-valid rung of the certifier uses only the
certified Frobenius norm of ``J = \operatorname{Choi}(\Phi^2-\Phi)``, with no SDP
solver and no eigendecomposition, so it never hits the degenerate-eigenvalue wall.
The Watrous-style inequality chain bounds the diamond norm both ways:

```math
\eta \in \left[\,\frac{\|J\|_F}{n},\ 2\|J\|_F\,\right],
\qquad \frac{\mathrm{hi}}{\mathrm{lo}} = 2n,
```

(`approximate_algebras.tex:347-352`, the upper bound from a Schur factorisation
of the block-PSD constraint, the lower from the maximally-entangled state). The
ratio ``2n`` is the price of certifying without solving: the bracket is built to
*contain* the truth and never abort, not to compete on tightness. The MOSEK rung
feeds two independent feasible points, each restored to exact feasibility in arb,
and collapses the width to roughly the MOSEK duality gap plus the arb radius. The
mathematics of both rungs is in [Certified arithmetic](@ref).

!!! note "Loose by design"
    The eig-free bracket width grows like ``2n`` and is meant to. It is the
    fallback the tight certifier dispatches to near ``\eta = 0``, and the default
    on the solver-free path. For tightness, load MOSEK.

## Design decisions

These are choices, not limits — places where the obvious thing would be wrong.

**Solver-free by default; MOSEK optional.** The default path needs no SDP solver:
the eig-free bracket and the flat-double ccall surface stand alone, and the
solver-free suite of 219 tests passes in about 2m09s. The Watrous diamond-norm
SDP — the exact value and the tight bracket — lives only in the optional
`AICMosekExt` extension. The tradeoff is stated plainly: solver-free by default,
loose bracket by design, tight only with MOSEK.

**The round-trip target is ``P_B``, not ``1_B``.** The factorization round-trip
``\Upsilon\Delta`` returns the block-diagonal conditional expectation ``P_B``,
not the full identity ``1_B``. ``\Delta`` reads only ``\mathcal{B}``'s
block-diagonal coordinates, so ``\Upsilon\Delta`` is zero on off-block units.
Comparing the round-trip to the full ``1_B`` reads ``1.0`` on every off-block
unit even at ``\eta = 0`` — a test that cannot pass — so the package sweeps only
in-block units and targets ``P_B`` (`FINDINGS §C14(b)`, `§C22(b)`). The ``\eta =
0`` oracle confirms this: the round-trip defect ``\|\Upsilon\Delta - P_B\|_{\mathrm{cb}}``
collapses to a measured maximum of ``3.9\mathrm{e}{-75}`` at
``\mathrm{prec} = 128`` — see [The η = 0 oracle](@ref).

**The product is the Choi–Effros star ``X \star Y = \tilde\Phi(XY)``, not ``XY``.**
Every multiplication inside the associated algebra ``\mathcal{A} =
\operatorname{Img}\tilde\Phi`` goes through ``\tilde\Phi``, because plain operator
multiplication leaves the image. This is the single most common test-blindness in
the project: the ``\eta = 0`` identity oracle is structurally blind to it (there
``\tilde\Phi = \mathrm{id}``, so ``\star`` is plain multiplication), so every star
claim is checked on a genuinely oblique ``\eta > 0`` fixture (`FINDINGS §C2`).

**The hypothesis is in cb-norm, not operator norm.**
``\eta``-idempotence is ``\|\Phi^2-\Phi\|_{\mathrm{cb}} \le \eta``, a supremum
over all ampliations ``1_{M_n} \otimes \Phi``, not the bare operator norm. The
spectral radius ``\rho(\Phi^2-\Phi)``, the operator norm, and ``\eta_{\mathrm{cb}}``
are three distinct numbers that coincide only when the defect superoperator is
normal; a construction can sit inside the spectral basin ``\rho < 1/4`` yet have
``\eta_{\mathrm{cb}}`` on the hypothesis boundary. Adversarial fixtures are
therefore calibrated on ``\eta_{\mathrm{cb}}``, never on ``\rho`` (`FINDINGS §C17`).

## The round-trip picture

The ``\eta``-dependence of the round-trip defect — including the ``P_B`` target
and the ``O(\eta)``-TP behaviour above — is plotted in the five-verb pipeline
tutorial (`docs/assets/roundtrip.png`); see [The five-verb pipeline](@ref) rather
than duplicating the figure here.

## See also

- [The mathematics](@ref) — why the star, the cb-norm, and the ``P_B`` target
  are forced by the construction
- [Certified arithmetic](@ref) — the bracket inequality chain and the tight
  MOSEK rung
- [Stay in factorize's domain](@ref) — keeping ``\rho`` inside the convergence
  disc
- [Check trace-preservation of a channel](@ref) — the `iscptp` tolerance recipe
- [The η = 0 oracle](@ref) — the exact case where every defect collapses
  to machine precision
