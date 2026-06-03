# API reference

Every exported name of `AlmostIdempotentChannels`, grouped by role, rendered from its docstring.

## The module

```@docs
AlmostIdempotentChannels
```

## Index

```@index
Modules = [AlmostIdempotentChannels]
```

## Channel types

```@docs
UCPMap
CPMap
```

## The headline verbs

In pipeline order: certify the defect, build the algebra, find the isomorphism, factorize.

```@docs
certified_defect
associated_algebra
main_isomorphism
factorize
encode
decode
```

## Result types

```@docs
CertifiedBracket
EpsCStarAlgebra
CStarAlgebra
MainIsomorphism
ChannelFactorization
```

## Accessors

### UCPMap

```@docs
n
nkraus
kraus
isunital
choi
iscptp
Base.adjoint(::UCPMap)
```

`Base.adjoint` (as `Φ'`) is a Base extension returning the dual `UCPMap` with Kraus `{K_a†}`; it is not separately exported.

### CertifiedBracket

```@docs
width
midpoint
value
```

`Base.minimum`, `Base.maximum`, and `Base.in` are extended but not separately exported.

### EpsCStarAlgebra

```@docs
dim
ambient
epsilon
```

### CStarAlgebra

```@docs
blocks
dim_B
n_B
```

### MainIsomorphism

```@docs
cbnorm_forward
cbnorm_inverse
isodefect
```

`blocks(v::MainIsomorphism)` delegates to the underlying `CStarAlgebra` and is listed under [CStarAlgebra](@ref) above.

### ChannelFactorization

```@docs
algebra
delups_defect
upsdel_defect
```

## Low-level and library discovery

```@docs
choi_diff
eta_eigfree
libaic_path
set_libaic_path!
```

## MOSEK extension (solver-gated)

These names require the `AICMosekExt` package extension. Without it they throw a helpful install hint. The solver-free alternative is [`certified_defect`](@ref) (always available, no solver needed). See [Install and use the MOSEK extension](@ref).

```@docs
idempotency_defect
eta_idempotence
diamond_norm_watrous
diamond_norm_watrous_primal
diamond_norm_dual
```
