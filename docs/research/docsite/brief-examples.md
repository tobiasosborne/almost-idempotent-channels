# brief-examples.md — verified docsite example corpus

Sources read (read-only): `test/fixtures.jl`, `test/test_core.jl`,
`test/test_factorize.jl`, `test/test_sdp.jl`, `test/test_aqua.jl`,
`test/runtests.jl`, `docs/plots/{containment,universality,roundtrip,eta0_oracle}.jl`,
`docs/src/{index,tutorial}.md`, `README.md`.

All numbers below come from those files verbatim; none are invented.
When a float is marked "not stated; must be generated at build time" it means the
source files contain the assertion invariant but not the printed value.

---

## 1. Fixture catalog

One row per constructor in `test/fixtures.jl`.
"Oracle" = cleanest η=0 ground truth per CLAUDE.md Rule 6.
"In factorize domain" = ρ(Φ²−Φ) < 0.10 conservative threshold (the tighter
`factorize` gate, not the wider prop_P basin ρ < 1/4).

| Name | Signature | Channel (words + math) | η=0? | B-blocks (when determinable) | In factorize domain? | Closed-form defect |
|---|---|---|---|---|---|---|
| `identity_kraus` | `identity_kraus(d)` | Identity map Φ=id on C^d; single Kraus K=I_d; Φ(X)=X | YES (η=0 ORACLE) | [d] | YES (η=0) | 0 |
| `block_cond_exp_kraus` | `block_cond_exp_kraus(d, m)` | Block conditional expectation; K0=diag(1..1,0..0) (m ones), K1=complement; Φ(X)=K0†XK0 + K1†XK1 | YES (η=0 ORACLE) | [m, d-m] | YES (η=0) | 0 |
| `dephasing_kraus` | `dephasing_kraus(d)` | Dephasing / diagonal projection; K_i = \|e_i><e_i\|; Φ(X)=diag(X); commutative B | YES (η=0 ORACLE) | [1,1,...,1] (d ones) | YES (η=0) | 0 |
| `depolarizing_kraus` | `depolarizing_kraus(d)` | Maximally depolarizing; K_{ij}=(e_i e_j†)/√d; Φ(X)=Tr(X)I/d | YES (η=0; exact idempotent) | [1] (scalar multiples of I) | YES (η=0) | 0 |
| `phit_kraus` | `phit_kraus(d, t)` | Φ_t=(1-t)·id + t·Dep_d; convex mix of identity and depolarizing | NO (η>0 for 0<t<1) | [d] (single block, B=M_d) | Small t only: t=0.1 is IN (ρ≈0.09); t=0.3 is OUT (ρ≈0.21) | t(1-t)·2(1-1/d²) |
| `paper_example_kraus` | `paper_example_kraus(etap)` | Kitaev paper example (.tex:367–378); Φ(X)=P0·Tr(g0·X)+P1·Tr(g1·X) on C²; asymmetric qubit channel | NO (η>0 for 0<etap<1) | [2] (single 2-block) | Small etap only (not stated numerically) | etap·√(1-etap) |
| `complex_qubit_kraus` | `complex_qubit_kraus()` | Fixed complex asymmetric unital qubit channel; K0=[0.6,0;0.8i,0], K1=[0,0.8i;0,0.6]; genuinely complex Choi (max\|Im(Choi(Φ²−Φ))\|≈0.61) | NO (η>0) | not stated | Not stated | 1.28 (from T5 SDP, test_sdp.jl) |
| `gadc_kraus` | `gadc_kraus(gamma, p)` | Generalized amplitude damping channel; 4 Kraus operators A0..A3; η=0 at γ=0 and γ=1, interior η>0 | Boundary only (γ=0 or γ=1) | not stated | not stated | None (analytic form not given; η=0 at boundary only) |
| `compress_idemp_kraus` | `compress_idemp_kraus(d, m)` | Proper-carrier compression: K0=Π_M (projection onto first m dims), K_{1+b}=\|e_0><e_{m+b}\|; exactly idempotent | YES (η=0 ORACLE) | [m] (single block) | YES (η=0) | 0 |
| `fixed_unitary` | `fixed_unitary(nn)` | Fixed basis-mixing unitary V=exp(iH); H has real super-diagonal 0.4, complex 0.3i second super-diagonal, diagonal phases 0.1k; used as conjugator only, not a channel | n/a (not a channel) | n/a | n/a | n/a |
| `mixconj_kraus` | `mixconj_kraus(d, m, t)` | Convex mix of compress_idemp(d,m) with its V-conjugate; Φ_t=(1-t)Φ_base + t·Ad_{V†}∘Φ_base∘Ad_V; SINGLE-BLOCK oblique η>0 | NO (η>0 for t>0) | [m] (single block B=M_m) | t=0.02: IN (ρ≈0.04); larger t not stated | Not stated (must be generated at build time) |
| `bce_conj_kraus` | `bce_conj_kraus(d, m, t)` | Convex mix of block_cond_exp(d,m) with its V-conjugate; MULTI-BLOCK oblique η>0; triggers §C14 Choi→Kraus PSD-cone clip so decode is only O(η)-TP | NO (η>0 for t>0) | [m, d-m] (multi-block B=M_m⊕M_{d-m}) | t=0.02 and t=0.05: IN | Not stated (must be generated at build time) |

### Oracle vs oblique classification

**Clean η=0 oracles** (iscptp(decode)=true at 1e-9; every defect < 1e-9 after
factorize; block structure exactly as listed): `identity_kraus`, `block_cond_exp_kraus`,
`dephasing_kraus`, `depolarizing_kraus`, `compress_idemp_kraus`.

**Genuinely η>0 oblique** (mandatory for exercising Choi–Effros star and dual
direction; decode only O(η)-TP for multi-block): `phit_kraus` (small t),
`bce_conj_kraus` (small t), `mixconj_kraus` (small t), `gadc_kraus`
(interior γ), `paper_example_kraus`, `complex_qubit_kraus`.

---

## 2. Verified example snippets

Each snippet is self-contained (inlines all fixtures it needs; fixtures are NOT
exported). Expected output is described after each block. The tag [jldoctest] means
the output is structurally stable and safe for exact-match checking; [@example]
means the output contains floats and must use the Documenter @example block form.

### 2a. Quickstart — first certified bound (doc page: index.md / quickstart)

Tag: [@example] for the bracket display; [jldoctest] for the (n, nkraus) line.

```julia
using AlmostIdempotentChannels, LinearAlgebra

# Φ_t = (1−t)·id + t·(maximally depolarizing) on C^d.
# Closed-form defect: ‖Φ_t²−Φ_t‖_cb = t(1−t)·2(1−1/d²).
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

Φ = UCPMap(phit_kraus(2, 0.1))
n(Φ), nkraus(Φ)   # structural fact: (2, 5)  ← 1 identity + 2²=4 depolarizing Kraus
```

Expected [jldoctest] output: `(2, 5)`

```julia
η = certified_defect(Φ)   # rigorous ‖Φ²−Φ‖_cb bracket — NO solver required
```

Expected [@example]: prints a CertifiedBracket display. The bracket must CONTAIN
the analytic value `0.1·0.9·2·(1-1/4) = 0.135`. The exact printout is not stable
to last digits; use @example.

```julia
# The bracket is the proof: test containment and read its geometry.
analytic_eta = 0.1 * 0.9 * 2 * (1 - 1/2^2)   # = 0.135
analytic_eta in η   # true (the test invariant; §B9 slack: minimum(η)-1e-7 ≤ analytic_eta)
```

Expected [jldoctest]-safe structural check: `true`.

```julia
width(η), midpoint(η)
```

Expected [@example]: two floats. The eig-free bracket is loose by design (hi/lo ~ 2n);
the exact numbers are not stated in the sources and must be generated at build time.

### 2b. η=0 oracle deep-dive (doc page: tutorial.md Step 5 / standalone oracle page)

These three are ALREADY partially present in tutorial.md (the block_cond_exp(4,2)
oracle in Step 5). The snippets below extend to all three standard oracles and add
the defect assertions. All structural outputs are [jldoctest]-safe.

Tag: [jldoctest] for block lists and boolean assertions; [@example] for raw defect floats.

```julia
using AlmostIdempotentChannels, LinearAlgebra

# Three exactly-idempotent conditional expectations (η = 0 ground truth).
identity_kraus(d) = [Matrix{ComplexF64}(I, d, d)]

function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end

function dephasing_kraus(d)
    K = Matrix{ComplexF64}[]
    for i in 1:d
        Ki = zeros(ComplexF64, d, d); Ki[i, i] = 1; push!(K, Ki)
    end
    return K
end

# Oracle 1: identity on C^2 → B = M_2 (single 2-block)
F1 = factorize(UCPMap(identity_kraus(2)))
blocks(algebra(F1))   # [2]

# Oracle 2: block conditional expectation on C^4 → B = M_2 ⊕ M_2
F2 = factorize(UCPMap(block_cond_exp_kraus(4, 2)))
blocks(algebra(F2))   # [2, 2]

# Oracle 3: dephasing on C^3 → B = M_1 ⊕ M_1 ⊕ M_1 (commutative)
F3 = factorize(UCPMap(dephasing_kraus(3)))
blocks(algebra(F3))   # [1, 1, 1]
```

Expected [jldoctest] output for `blocks(algebra(F1))`: `1-element Vector{Int64}: 2`
Expected [jldoctest] output for `blocks(algebra(F2))`: `2-element Vector{Int64}: 2  2`
Expected [jldoctest] output for `blocks(algebra(F3))`: `3-element Vector{Int64}: 1  1  1`

```julia
# Every certified defect is below 1e-9 (the test invariant):
maximum(delups_defect(F2)) < 1e-9   # true; measured ≈ 4.4e-75 at prec=128
maximum(upsdel_defect(F2)) < 1e-9   # true; measured ≈ 3.9e-75 at prec=128
iscptp(decode(F2))                   # true at default atol=1e-9 (η=0 oracle: no PSD clip)
iscptp(encode(F2))                   # true
```

Expected [jldoctest] output: all `true`. The measured sub-floor values (4.4e-75,
3.9e-75) come from `README.md` "cross-check ladder" section, attributed to the
η=0 oracle test at prec=128; they are the values printed in the README but the
exact digits are not reproduced in any test assertion (tests only assert `< 1e-9`).
Use them for documentation prose, not as jldoctest numeric literals.

### 2c. Reading a CertifiedBracket (doc page: index.md or API reference)

Tag: [jldoctest] for structural/boolean outputs; [@example] for float values.

```julia
using AlmostIdempotentChannels

# CertifiedBracket: lo ≤ x ≤ hi.
b = CertifiedBracket(0.1, 0.3; label="demo")
minimum(b)    # 0.1
maximum(b)    # 0.3
width(b)      # 0.2
midpoint(b)   # 0.2
0.2 in b      # true
0.5 in b      # false
value(b)      # nothing (no point estimate; only MOSEK tight bracket carries value)
```

Expected [jldoctest] outputs: `0.1`, `0.3`, `0.2`, `0.2`, `true`, `false`, `nothing`.
All are exact; this is the purest [jldoctest] target in the corpus.

```julia
# With a point estimate (the MOSEK tight bracket sets value):
bv = CertifiedBracket(0.1, 0.3; value=0.2, label="demo")
value(bv)   # 0.2
```

Expected [jldoctest] output: `0.2`.

```julia
# Validating constructor: lo > hi throws ArgumentError.
CertifiedBracket(0.5, 0.1)   # throws ArgumentError
```

Expected: `ArgumentError`. Use @example or a try/catch jldoctest form.

### 2d. Algebra block-structure extraction (doc page: tutorial.md Steps 2–3 / API)

Tag: [jldoctest] for block lists, dims, ambient; [@example] for epsilon bracket.

```julia
using AlmostIdempotentChannels, LinearAlgebra

e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

Φ = UCPMap(phit_kraus(2, 0.1))

A = associated_algebra(Φ)
dim(A)      # 4  (identity on C²: dim_A = n² = 4 when η small, image fills M_n)
ambient(A)  # 2  (the ambient Hilbert space dimension n = 2)
# epsilon(A) is a CertifiedBracket: 0.0 ∈ epsilon(A) is FALSE for η>0 (not an exact
# conditional expectation). The bracket width and midpoint are float-heavy (@example).

v = main_isomorphism(Φ)
blocks(v)   # [2] — the single M_2 block of B for a qubit channel
```

Expected [jldoctest] output for `dim(A)`: `4`.
Expected [jldoctest] output for `ambient(A)`: `2`.
Expected [jldoctest] output for `blocks(v)`: `1-element Vector{Int64}: 2`.

Note: `dim(A) = 4` for `phit_kraus(2, 0.1)` is confirmed by T6 (`associated_algebra`
test) via `dim(A) == 4` assertion in `test_core.jl` for `identity_kraus(2)`;
for phit(2,0.1) the dim is not explicitly stated in a test assertion — it is
"not stated; must be verified at build time". Use the identity(2) oracle for
the dim=4 jldoctest instead:

```julia
A0 = associated_algebra(UCPMap(identity_kraus(2)))
dim(A0)      # 4 (confirmed by test_core.jl T6)
ambient(A0)  # 2
maximum(epsilon(A0)) < 1e-9   # true (η=0 → ε ~ machine)
```

Expected [jldoctest]: `4`, `2`, `true`.

### 2e. Multi-block oblique factorize (doc page: tutorial.md Step 4)

This is ALREADY covered in `tutorial.md` Steps 4 and 5. The snippet below is the
minimal self-contained form for a standalone example page or API docstring.
Tag: [jldoctest] for `blocks`, dims, boolean `iscptp`; [@example] for defect brackets.

```julia
using AlmostIdempotentChannels, LinearAlgebra

function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end
function fixed_unitary(nn)
    H = zeros(ComplexF64, nn, nn)
    for p in 1:(nn - 1); H[p, p + 1] = 0.4; H[p + 1, p] = 0.4; end
    for p in 1:(nn - 2); H[p, p + 2] = 0.3im; H[p + 2, p] = -0.3im; end
    for p in 1:nn; H[p, p] = 0.1 * p; end
    return exp(im * Matrix(Hermitian(H)))
end
function bce_conj_kraus(d, m, t)
    base = block_cond_exp_kraus(d, m); V = fixed_unitary(d); Vd = V'
    K = Matrix{ComplexF64}[]
    for Ka in base; push!(K, sqrt(1 - t) * Ka); end
    for Ka in base; push!(K, sqrt(t) * (Vd * Ka * V)); end
    return K
end

# bce_conj(4, 2, 0.02): multi-block oblique η>0 channel; B = M_2 ⊕ M_2.
Ψ = UCPMap(bce_conj_kraus(4, 2, 0.02))
F = factorize(Ψ)
blocks(algebra(F))   # [2, 2]
```

Expected [jldoctest]: `2-element Vector{Int64}: 2  2`.

```julia
# B1 dimension assertions (a Dec/Enc swap would make these fail):
en = encode(F); de = decode(F)
(en.dim_in, en.dim_out)   # (4, 4)  — n_B(B) = 2+2 = 4, N = 4
(de.dim_in, de.dim_out)   # (4, 4)
```

Note: for bce_conj(4,2,0.02) n_B=4 and N=4 so both CPMaps are square (dim_in=dim_out=4).
The swap detection relies on the ordering: (en.dim_in, en.dim_out) == (de.dim_out, de.dim_in)
(confirmed by test_factorize.jl). For a non-square example use mixconj(4,2,0.02)
where n_B=2, N=4, giving en=(2,4) and de=(4,2).

```julia
# Round-trip defects are O(η), not 0 and not O(1):
η = F.eta_proxy
maximum(delups_defect(F)) < 50 * η   # true (test bound, test_factorize.jl)
maximum(upsdel_defect(F)) < 50 * η   # true

# The decode is only O(η)-trace-preserving (§C14 PSD-cone clip):
iscptp(encode(F))                     # true (encode always TP)
iscptp(decode(F))                     # false at default atol=1e-9
iscptp(decode(F); atol = 1e-4)        # true at loosened tolerance
```

Expected [jldoctest]: all five booleans `true, true, true, false, true` in order.

### 2f. Out-of-domain ArgumentError demo (doc page: tutorial.md warning box / Known limits)

Tag: [@example] with try/catch, or a jldoctest with `# output` matching the error type.
The exact error message contains "convergence domain" (confirmed by test_core.jl).

```julia
using AlmostIdempotentChannels, LinearAlgebra

e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

# phit(2, 0.3): in the prop_P basin (ρ(Φ²−Φ) ≈ 0.21) but OUTSIDE factorize's
# tighter convergence domain (ρ < 0.10). Throws ArgumentError, not a C abort.
try
    factorize(UCPMap(phit_kraus(2, 0.3)))
catch e
    println(typeof(e))         # ArgumentError
    println(occursin("convergence domain", e.msg))   # true
end
```

Expected [@example] output:
```
ArgumentError
true
```

```julia
# certified_defect is ALWAYS safe at any η — the escape hatch the error message
# recommends. Even with t=0.3 (out of factorize domain):
b = certified_defect(UCPMap(phit_kraus(2, 0.3)))
b isa CertifiedBracket   # true
```

Expected [jldoctest]: `true`.

```julia
# The specific closed-form value (needs MOSEK extension):
# idempotency_defect(UCPMap(phit_kraus(2, 0.3))) == 0.3149999999994273
# = 0.3 * 0.7 * 2 * (1 - 1/4) = 0.315 exactly.
# Without MOSEK, the solver-free certified_defect bracket still CONTAINS 0.315.
analytic = 0.3 * 0.7 * 2 * (1 - 1/2^2)   # 0.315
analytic in certified_defect(UCPMap(phit_kraus(2, 0.3)))   # true (§B9 slack)
```

Expected [jldoctest]: `true`.

### 2g. Trace-preservation / iscptp atol behaviour (doc page: Known limits / tutorial)

This is ALREADY demonstrated in tutorial.md Step 4. Minimal standalone form:

Tag: [jldoctest] for all five booleans.

```julia
# η=0 oracle: both encode and decode are TP to machine precision.
# (see snippet 2b above)
iscptp(decode(F2))   # true (η=0: no PSD-cone clip)
iscptp(encode(F2))   # true

# η>0 multi-block oblique: decode is only O(η)-TP.
# (see snippet 2e above)
iscptp(encode(F))    # true
iscptp(decode(F))    # false  (atol=1e-9 default)
iscptp(decode(F); atol = 1e-4)   # true (O(η)-TP; measured TP defect 1.9e-6 at t=0.02)
```

The measured TP defect 1.9e-6 at t=0.02 and 2.2e-5 at t=0.05 come from
`fixtures.jl` comments. These are stated as measured values in the fixture
docstring; the tests only bound them as `1e-9 < tpd < 1e-4` (test_factorize.jl).

### 2h. MOSEK extension — guarded example (doc page: MOSEK section of index.md / API)

These examples CANNOT run in the default no-solver build. They must be either:
- Non-executing code blocks (` ```julia ` without `@example` / `jldoctest`), or
- `@example` blocks with a `# hide` fence that only runs when MOSEK is available.
Documenter does not support conditional block execution natively; use non-executing
code blocks and note that the output is from a passing test.

```julia
# Requires: using Convex, Mosek, MosekTools   (activates AICMosekExt)
using AlmostIdempotentChannels, LinearAlgebra
# ... (phit_kraus inline as above) ...

Φ = UCPMap(phit_kraus(2, 0.3))

# Exact diamond-norm value (Watrous SDP):
idempotency_defect(Φ)   # 0.3149999999994273  (= 0.3·0.7·1.5; from README.md)

# Tight rigorous bracket (arb fed SDP feasible points):
tight = certified_defect(Φ; tight=true)
width(tight)    # ~5.76e-13  (from README.md benchmark table; measured, not exact)
value(tight)    # carries the SDP point estimate (≈ 0.315)
```

Numbers from `README.md` benchmark table and MOSEK section:
- `idempotency_defect(phit(2,0.3))` = `0.3149999999994273` (README.md MOSEK section)
- tight bracket width for phit(2,0.3): `5.76e-13` (README.md benchmark table)
- tight bracket width for phit(2,0.1): `1.60e-13`
- tight bracket width for paper(0.1): `4.60e-14`
- measured max primal-dual gap (T5, test_sdp.jl): `1.23e-11` (README.md cross-check ladder)

---

## 3. Numbers table with provenance

Every concrete measured number worth quoting in the documentation.

| Number | What it measures | Source file | Conditions |
|---|---|---|---|
| `0.0` | Exact defect of identity_kraus, block_cond_exp_kraus, dephasing_kraus, depolarizing_kraus, compress_idemp_kraus | fixtures.jl comments | η=0 definition |
| `t(1-t)·2(1-1/d²)` | Closed-form ‖Φ_t²−Φ_t‖_cb = ‖Φ_t²−Φ_t‖_⋄ for Φ_t | fixtures.jl comment, test_core.jl, test_sdp.jl | phit_kraus(d,t) |
| `0.1·0.9·2·(1-1/4) = 0.135` | Analytic η for phit(2,0.1) | test_core.jl T2/T3 | d=2, t=0.1 |
| `0.3·0.7·2·(1-1/4) = 0.315` | Analytic η for phit(2,0.3) | test_core.jl T2/T3, test_sdp.jl | d=2, t=0.3 |
| `0.3149999999994273` | `idempotency_defect(phit(2,0.3))` from MOSEK SDP | README.md MOSEK section | MOSEK required |
| `etap·√(1-etap)` | Closed-form ‖Φ²−Φ‖_cb for paper_example_kraus(etap) | fixtures.jl comment, test_sdp.jl | paper_example |
| `0.1·√0.9 ≈ 0.09487` | Analytic η for paper_example_kraus(0.1) | test_core.jl T2/T3, test_sdp.jl | etap=0.1 |
| `1.28` | η for complex_qubit_kraus (from T5 SDP reference) | test_sdp.jl cases table | complex_qubit |
| `1e-10` | Max-entry diff threshold, C/arb ccall vs pure-Julia Choi(Φ²−Φ) | test_core.jl T1 | round-trip cross-check |
| `1e-7` | §B9 slack tolerance for certified_defect containment of analytic η | test_core.jl T2/T3 | solver-free eig-free |
| `1e-9` | Threshold: maximum(epsilon(A)) < 1e-9 for η=0 oracles | test_core.jl associated_algebra testset | η=0 |
| `1e-9` | Threshold: maximum(delups/upsdel) < 1e-9 for η=0 factorize oracle | test_factorize.jl | η=0 |
| `1e-9` | Default atol for iscptp | test_factorize.jl, tutorial.md | all |
| `1e-4` | Loosened atol for iscptp on η>0 multi-block decode | test_factorize.jl, tutorial.md, fixtures.jl | bce_conj η>0 |
| `50·η` | Generous O(η) envelope: `maximum(delups/upsdel) < 50·η` | test_factorize.jl | bce_conj, mixconj; t=0.02,0.05 |
| `1.9e-6` | Measured TP defect ‖ΣK†K−I‖ of decode for bce_conj(4,2,0.02) | fixtures.jl comment | multi-block oblique |
| `2.2e-5` | Measured TP defect ‖ΣK†K−I‖ of decode for bce_conj(4,2,0.05) | fixtures.jl comment | multi-block oblique |
| `3.7e-6` | "Measured clipped mass" on multi-block oblique decode | README.md Known limits | bce_conj |
| `4.4e-75` | Measured max delups_defect for η=0 oracle at prec=128 | README.md cross-check ladder | prec=128, block_cond_exp(4,2) |
| `3.9e-75` | Measured max upsdel_defect for η=0 oracle at prec=128 | README.md cross-check ladder | prec=128, block_cond_exp(4,2) |
| `1.23e-11` | Max primal-dual gap (T5 Watrous SDP strong duality) | README.md cross-check ladder | MOSEK; all T5 fixtures |
| `5.76e-13` | Tight bracket width for phit(2,0.3) | README.md benchmark table | MOSEK |
| `1.60e-13` | Tight bracket width for phit(2,0.1) | README.md benchmark table | MOSEK |
| `4.60e-14` | Tight bracket width for paper(0.1) | README.md benchmark table | MOSEK |
| `5.46e-01` | Eig-free (solver-free) bracket width for phit(2,0.3) | README.md benchmark table | no solver |
| `2.34e-01` | Eig-free bracket width for phit(2,0.1) | README.md benchmark table | no solver |
| `2.01e-01` | Eig-free bracket width for paper(0.1) | README.md benchmark table | no solver |
| `0.018, 0.047, 0.025, 0.041, 0.013` | c=isodefect/η at dim_A=8,12,16,18,20 | docs/plots/universality.jl (from test_dbo3, prec=256, 2026-06-02) | C dim-sweep |
| `0.047` | Max c=isodefect/η over dim_A sweep (worst case at dim_A=12) | docs/plots/universality.jl, README.md | C dim-sweep, prec=256 |
| `-2.707458e-4` | OLS slope of c vs dim_A (flat = dimension-independent) | docs/plots/universality.jl | C dim-sweep |
| `20.0` | Upper bound for c_ups and c_del in §D2 dim-independence canary | test_factorize.jl §D2 | Julia mixconj t=0.02, n=4,5,6 |
| `1.5` | Max-to-min ratio for c_ups and c_del across n=4,5,6 (flatness gate) | test_factorize.jl §D2 | mixconj t=0.02 |
| `0.10` | Conservative ρ threshold for factorize domain | docs/src/index.md, tutorial.md warning | factorize |
| `0.21` | ρ(Φ²−Φ) estimate for phit(2,0.3) (outside factorize domain) | test_core.jl comment | phit d=2 t=0.3 |
| `0.04` | ρ(Φ²−Φ) estimate for mixconj(4,2,0.02) (inside factorize domain) | fixtures.jl comment | mixconj t=0.02 |
| `2` | ρ(Φ²−Φ) for reflection conjugation Φ(X)=UXU†, U=diag(1,-1) | test_core.jl BASIN GUARDS comment | genuine out-of-basin |
| `1e-3` | Lower bound for η_proxy of oblique fixtures in factorize | test_factorize.jl | bce_conj, mixconj; all t |
| `2` (n_B) | n_B(B) for mixconj(4,2,0.02) → B=M_2 | test_factorize.jl B1 dims | single-block |
| `4` (N) | ambient N for mixconj(4,2,0.02) on C^4 | test_factorize.jl B1 dims | single-block |
| `~1e-12` | double vs arb@prec=53 agreement threshold (cross-check rung 2) | CLAUDE.md cross-check ladder | internal |
| `219` | Test count (PASS) for the solver-free suite | README.md install section | no MOSEK, ~2m09s |

---

## 4. Pitfalls for doc examples

These are the things most likely to break a documentation build or produce
misleading output. Each is actionable.

**P1. Fixtures are not exported — must inline.**
`phit_kraus`, `block_cond_exp_kraus`, `dephasing_kraus`, `identity_kraus`,
`fixed_unitary`, `bce_conj_kraus`, `mixconj_kraus`, `compress_idemp_kraus`,
`paper_example_kraus`, `complex_qubit_kraus`, `gadc_kraus` all live in
`test/fixtures.jl` which is not part of the package API. Every doc example that
uses one must inline the full constructor. The tutorial.md already does this
correctly; the trap is a copy-paste that omits the function body.

**P2. factorize domain — do not use large-η fixtures in factorize examples.**
`phit_kraus(2, 0.3)` has ρ ≈ 0.21, which is outside the factorize threshold of
ρ < 0.10. Using it in a factorize example will produce an `ArgumentError`, not
a result. Safe choices for factorize: any η=0 oracle, `phit_kraus(2, 0.1)`,
`bce_conj_kraus(4, 2, 0.02)`, `bce_conj_kraus(4, 2, 0.05)`,
`mixconj_kraus(4, 2, 0.02)`. The ρ threshold for `associated_algebra` and
`main_isomorphism` is the wider ρ < 1/4 basin; `certified_defect` has no domain
restriction.

**P3. iscptp(decode(F)) is false at the default atol=1e-9 for η>0 multi-block.**
Any example that calls `iscptp(decode(F))` on an η>0 multi-block fixture
(bce_conj) without specifying `atol` will print `false`. This is correct
behaviour (the §C14 PSD-cone clip) but looks like a broken result to a reader
expecting `true`. Either (a) use `atol = 1e-4` and note the O(η)-TP design, or
(b) use an η=0 oracle where decode IS TP to machine precision.

**P4. MOSEK examples cannot run in the default build — must be non-executing blocks.**
`idempotency_defect`, `certified_defect(Φ; tight=true)`, `diamond_norm_watrous_primal`,
`diamond_norm_dual` all require the `AICMosekExt` extension (activated by
`using Convex, Mosek, MosekTools`). Without a MOSEK license they throw a
`MethodError` or a helpful install hint, not a result. Mark all MOSEK examples as
plain ` ```julia ` (non-executing) or wrap in a HAVE_MOSEK guard. Do NOT use
`@example` or `jldoctest` for MOSEK examples unless Documenter is configured to
skip them conditionally.

**P5. CertifiedBracket floats are not stable to last digits — use @example not jldoctest.**
The bracket lo and hi values printed by `certified_defect`, `epsilon(A)`,
`isodefect`, `delups_defect`, `upsdel_defect`, `cbnorm_forward`, `cbnorm_inverse`
all depend on the arb working precision and the outward rounding; they are not
reproducible to the last ULP across platform and precision settings. Use `@example`
blocks (which embed the actual output at build time) for any snippet that prints
a bracket. Only structural outputs (block lists, dims, boolean assertions, type
names) are safe for exact-match `jldoctest`.

**P6. The §B9 slack means `η ∈ certified_defect(Φ)` can be false without a tolerance.**
The test uses `minimum(b) - tol <= η <= maximum(b) + tol` with `tol = 1e-7`, not
`η in b` (which calls `Base.in` with no slack). The arb FLOOR-rounded lo can sit
a hair above the analytic η. A doc example that writes `analytic_eta in b` without
the tolerance may print `false` for some floating-point alignments. Either use the
tolerance explicitly or state "within a 1e-7 margin" in prose.

**P7. The `value` field of the solver-free CertifiedBracket is `nothing`.**
`value(certified_defect(Φ))` returns `nothing` (solver-free path; no point
estimate). Only `certified_defect(Φ; tight=true)` (MOSEK required) sets a
non-nothing value. A jldoctest example that calls `value(b)` on a solver-free
bracket and expects a float will print `nothing` instead.

**P8. `blocks(v)` vs `blocks(algebra(F))` — different methods, same data.**
`v = main_isomorphism(Φ)` returns a `MainIsomorphism`; `blocks(v)` gives the block
list. `F = factorize(Φ)` returns a `ChannelFactorization`; `blocks(algebra(F))`
gives the same block list. Both are confirmed by the tests. In jldoctest setup
blocks that build both, make sure the variable names do not collide.

**P9. For bce_conj(4,2,t) both encode and decode are square (dim_in=dim_out=4) because n_B=N=4.**
The dim-swap detection test only catches a bug when en.dim_in == de.dim_out and
en.dim_out == de.dim_in AND these pairs differ. For bce_conj(4,2,t) n_B=4=N so
all dims are 4 and the assertion passes trivially for both orderings. To see the
rectangular case (where a swap would give different numbers) use mixconj(4,2,t):
n_B=2, N=4, so encode=(2,4) and decode=(4,2).

**P10. `associated_algebra` and `main_isomorphism` hit the prop_P guard for ρ ≥ 1/4.**
The reflection `UCPMap([ComplexF64[1 0; 0 -1]])` (ρ=2) throws `ArgumentError` with
"prop_P basin" in the message even for `associated_algebra`, not just `factorize`.
Only `certified_defect` is safe for any ρ. Do not use the reflection or other
genuinely out-of-basin channels as inputs to any verb except `certified_defect`.
