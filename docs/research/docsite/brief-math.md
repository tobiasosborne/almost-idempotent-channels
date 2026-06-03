# Background: the mathematics behind `AlmostIdempotentChannels.jl`

This page is the mathematical brief for the package. It explains, in four
stages, what the package computes and *why* — following

> Alexei Kitaev, *Almost-idempotent quantum channels and approximate
> $C^*$-algebras*, [arXiv:2405.02434](https://arxiv.org/abs/2405.02434) (2025).

Every nontrivial claim is anchored to a line of the canonical source
`paper/src/approximate_algebras.tex`, cited as `tex:NNN`. Subtleties that the
project's own hostile reviews keep catching are cited to the living log
`paper/FINDINGS.md`, as `FINDINGS §Cn` (load-bearing) or `FINDINGS §Dn` (open
question). The `.tex` is ground truth: where a paraphrase here would risk an
error, the inequality is quoted rather than restated.

The package implements **constructive finite-dimensional algorithms**. The
paper's proofs are frequently non-constructive (Lefschetz–Hopf fixed points,
holomorphic functional calculus by contour integral, Haar-measure diagonals);
in finite dimensions the objects those proofs merely assert to exist are
*computable*, and the package computes them directly. For each stage we
therefore separate **(b) the paper's proof technique** from **(c) the
constructive route the package takes**. The paper's $O(\eps)$ bound is the
correctness specification: the constructive routine must produce an output that
meets the bound the theorem promises.

---

## The four-stage story

| Stage | Object | Key result | `tex:` |
|---|---|---|---|
| 0 | UCP map $\Phi:\Bo(\calH)\to\Bo(\calH)$, cb-norm, $\eta$-idempotence | Choi/Stinespring rep | `tex:354`, `tex:1568` |
| 1 | regularisation $\wt\Phi=\theta(2\Phi-1)$; $\calA=\Img\wt\Phi$, star product | `prop_P` (§3) | `tex:524`, `tex:2171` |
| 2 | the $\eps$-$C^*$ algebra $\calA$, $\eps=O(\eta)$ | def. $\eps$-$C^*$ / $\delta$-hom (§2) | `tex:407`, `tex:443` |
| 3 | rigidity: $\calA$ is $O(\eps)$-isomorphic to a genuine $\calB$ | `th_main` (§2; proof §9) | `tex:460` |
| 4 | the channel factorization $\Phi\approx\Delta\Upsilon$ | `th_idemp_structure` (§11), `th_factorization` (§12) | `tex:318`, `tex:2730` |

---

## Stage 0 — the object: UCP maps and the completely bounded norm

### (a) Statement

We work in the **Heisenberg (observable) picture**: a *unital completely
positive* (UCP) map

$$\Phi\colon \Bo(\calH)\to\Bo(\calH),\qquad \dim\calH = n,$$

acting on observables, not states. In finite dimensions every UCP map has a
**Choi representation** (`tex:1568`, eq. `Choi`):

$$\Phi(X) = V^\dag (X\otimes 1_\calF)\,V,\qquad V^\dag V = 1_\calH,\qquad
V\colon\calH\to\calH\otimes\calF,$$

with auxiliary space $\calF$ (canonically $\calF=\calH^*\otimes\calH$,
`tex:1575`). Choosing an orthonormal basis $\{f_a\}$ of $\calF$ and writing
$V$ as a block column $V=\mathrm{col}(K_a)$ gives the **Kraus form**

$$\Phi(X) = \sum_a K_a^\dag X K_a,\qquad \sum_a K_a K_a^\dag = 1_\calH,$$

where the second equation is unitality. (The paper writes $V$ directly and does
not separately name "Kraus"; the equivalence is standard — shard G, `tex:1568`.)

The **Choi matrix** is $C_\Phi=\sum_{ij}E_{ij}\otimes\Phi(E_{ij})$. Then
$C_\Phi\ge 0$ iff $\Phi$ is completely positive, and the partial trace
$\Tr_{\calH}(C_\Phi)=1_\calH$ encodes unitality (shard G; classical Choi
[Cho75], implicit at `tex:1575`).

The norm that governs everything is the **completely bounded norm** (`tex:347`,
eq. unlabeled):

$$\|\Lambda\|_\cb = \sup_n\;\sup_{X\neq 0}
\frac{\|(1_{\Ma{n}}\otimes\Lambda)(X)\|}{\|X\|},\qquad
X\in\Ma{n}\otimes\calA.$$

It is a supremum over **all ampliations** $1_{\Ma{n}}\otimes\Lambda$, not the
bare operator norm. In finite dimensions it is computable by semidefinite
programming ([Watrous, Ch. 3], cited `tex:352`); its dual on linear functionals
is the **diamond norm** $\|\cdot\|_\Diamond$ (`tex:352`). $\Phi$ is

$$\boxed{\text{$\eta$-idempotent}\iff \|\Phi^2-\Phi\|_\cb\le\eta}\qquad(\texttt{tex:354}).$$

For any UCP map $\|\Phi\|_\cb=1$ (Schwarz inequality, `tex:1692`–`tex:1719`;
shard G).

### Why the cb-norm, not the operator norm

The hypothesis is phrased in $\|\cdot\|_\cb$ on purpose. The point of
$\eta$-idempotence is that it must survive **extending the physical system with
new, non-interacting parts** (`tex:346`); that is exactly an ampliation
$1_{\Ma{n}}\otimes\Phi$, which the cb-norm sees and the operator norm does not.
This is what lets the paper prove an *extended* (tensor-stable) version of every
statement, uniformly in $n$ (`tex:398`). Tensor extensions (§10) exist
*precisely because* the cb-norm sees all ampliations.

### The dual (state) picture and the direction of duality

States live in $\calB^*$. **Quantum channels** are CPTP maps of density
matrices, dual to UCP maps of observables (`tex:277`). $\Phi$ is idempotent iff
the channel $T=\Phi^*$ is idempotent (`tex:316`). The encode/decode channels of
Stage 4, $\Enc=(\Gamma\Co_\calM)^*$ and $\Dec=\Delta^*$, are the **duals** of
the observable-side maps $\Gamma\Co_\calM$ and $\Delta$ (`tex:334`). The
direction matters: UCP acts on observables, CPTP on states; getting the dual
backwards silently transposes everything (`FINDINGS §C22`, callout in CLAUDE.md).

### (b) Paper's technique vs (c) constructive route

The Choi/Stinespring representations are *already* algorithmic in finite
dimensions. The package builds $V$ from the Kraus list ancilla-major
(`V[a\cdot n + i, j] = K_a[i,j]`, so internally $\Phi(X)=V^\dag(1_\calF\otimes
X)V$ — note the ancilla side is swapped relative to the paper's $H$-left layout;
this is a deliberate, documented convention, `FINDINGS §C13(a)`). The Choi
matrix is assembled directly; CP is checked by a certified smallest-eigenvalue
of $C_\Phi$, unitality by a certified partial trace.

The genuinely subtle quantity is $\eta=\|\Phi^2-\Phi\|_\cb$ itself. The default
package route is **solver-free**: an eig-free arb bracket
$[\,\|J\|_F/N,\;2\|J\|_F\,]$ (rounded outward) that *certifies* $\eta$ without
an SDP, where $J$ is the Choi of the Hermiticity-preserving map $\Phi^2-\Phi$
(`FINDINGS §C17`, §C22(a)). The bracket is loose by design ($\mathrm{hi}/\mathrm{lo}\sim 2n$);
the optional MOSEK extension supplies the exact diamond-norm value.

### Subtleties that bite here

- **cb-norm $\neq$ operator norm.** Estimating $\eta$ requires the Choi
  matrix / ampliation, not just $\|\Phi^2-\Phi\|_{op}$. The spectral radius
  $\rho(\Phi^2-\Phi)$, the operator norm, and $\eta_\cb$ are three different
  numbers, and a construction can sit inside the spectral basin
  ($\rho<1/4$, Stage 1) yet have $\eta_\cb$ on the hypothesis boundary
  ($\eta_\cb=1/4-\kappa$) — they coincide only when the defect superoperator
  is normal (`FINDINGS §C17`). Calibrate adversarial fixtures on $\eta_\cb$,
  never on $\rho$.

---

## Stage 1 — regularisation: $\wt\Phi=\theta(2\Phi-1)$ and the oblique image

### (a) Statement

To extract an algebra one needs a genuinely *idempotent* map. The paper applies
functional calculus to the completely-bounded operator $\Phi$ itself
(`tex:354`–`tex:363`, `tex:2171`):

$$\wt\Phi = \theta(2\Phi-1)
= \tfrac12\bigl(1+\sgn(2\Phi-1)\bigr)
= \tfrac12\Bigl(1+(2\Phi-1)\bigl(1-4(\Phi-\Phi^2)\bigr)^{-1/2}\Bigr),
\qquad \theta(x)=\begin{cases}1&x>0\\0&x<0\end{cases}.$$

The Taylor series in $4(\Phi-\Phi^2)$ **converges if $\eta<1/4$** (`tex:2176`).
Then (`tex:2177`–`tex:2182`):

$$\wt\Phi^2=\wt\Phi,\qquad \|\wt\Phi-\Phi\|_\cb\le O(\eta),\qquad
\wt\Phi(1)=1,\qquad \wt\Phi(X^\dag)=\wt\Phi(X)^\dag.$$

$\wt\Phi$ is **not guaranteed completely positive** (`tex:363`, with the
explicit 2-dimensional counterexample at `tex:365`–`tex:388`). This is why the
$\eps$-$C^*$ algebra is built from $\wt\Phi$ but the *channel* factorization
(Stage 4) is a separate, more delicate construction.

The engine is **Proposition `prop_P`** (`tex:524`, the canonical projection
construction):

> Let $P$ satisfy $\|P^2-P\|\le\delta<1/4$, and let
> $\wt P=\theta(2P-I)$, $\theta(X)=\tfrac12(I+\sgn(X))$. Then $\wt P^2=\wt P$,
> $\wt P$ commutes with $P$, and $\|\wt P-P\|\le\|2P-I\|\cdot O(\delta)$.

Its proof (`tex:532`) applies the $\sgn$ calculus to $X=2P-I$, which satisfies
$\|X^2-I\|\le 4\delta<1$. The functional calculus has a **precondition**: $|X|$,
$\sgn(X)$, $\theta(X)$ are defined for $\|X^2-x_0^2 I\|<x_0^2$ (`tex:514`–`tex:516`),
i.e. $\|X^2-I\|<1$ for the involution case (`tex:520`). The sign of a near-involution
satisfies $\|\sgn(X)-X\|\le\|X\|\,O(\|X^2-I\|)$ (`tex:520`).

The algebra is then the image (`tex:2183`–`tex:2190`):

$$\calA=\Img\wt\Phi=\Ker(1-\wt\Phi)\subseteq\Bo(\calH),\qquad
X\star Y=\wt\Phi(XY)\quad\text{(Choi–Effros product)}.$$

The product is $X\star Y=\wt\Phi(XY)$, **not** the plain operator product $XY$,
which would leave $\Img\wt\Phi$ (`tex:341`, `tex:2189`).

### (b) Paper's technique: holomorphic functional calculus

The paper frames $\sgn$, $|\cdot|$, $\theta$ through the resolvent / contour
integral $f(X)=\frac{1}{2\pi i}\oint_C (zI-X)^{-1}f(z)\,dz$ (`tex:548`), with the
spectral-mapping theorem. It explicitly *declines* to rely on this:

> "We have not found use of these methods due to their fragility in the
> approximate setting and the lack of a $*$-representation, which could
> alleviate it." (`tex:535`)

It gives the reason: the spectrum of a non-normal operator is
perturbation-fragile — the $3\times3$ example $X$ with a $t$-entry has
eigenvalues $\sim t^{1/3}$ (`tex:540`–`tex:544`). The convergent route the paper
*does* use is the Taylor/Newton-style series above, plus the
inverse-function-theorem `lem_invfun` (`tex:564`) — which is already
constructive: the Banach contraction iteration $x_n=g_y(x_{n-1})$ with
$\|x_n-x_{n-1}\|<rc^{n-1}$ *is* the algorithm (`tex:591`).

### (c) Constructive route in the package

$\wt\Phi$ is computed by a **certified matrix-sign iteration** on the
superoperator $S=2S_\Phi-1$ (Newton–Schulz / Denman–Beavers / global Newton),
gated by an a-posteriori certificate. The basin guard is the **eig-free
Gelfand spectral radius** $\rho(\Phi^2-\Phi)<1/4$ (i.e. $\rho(4(\Phi-\Phi^2))<1$,
the spectral basin), *not* the stricter operator norm $\|\Phi^2-\Phi\|_{op}$
(`FINDINGS §C15`; the relaxation from op-norm to Gelfand $\rho$ is bead aic-8hz).
Out of basin, the pipeline **fails loud** at the `aic_assoc_regularize` entry
and at `aic_prop_P` (`FINDINGS §C15`) — it does not hang.

A subtlety in the certified path: when the input matrix carries a wide inherited
arb radius from an upstream matmul chain, the ball midpoint can drift off the
involution. The fix iterates on $\mathrm{mid}(X)$ (radius stripped once at entry)
and gates by a shared a-posteriori certificate $\|Y^2-I\|_F$ and $\|YX-XY\|_F$;
soundness comes from the away-from-zero basin precondition plus correct seeding
($Y_0=\mathrm{mid}(X)$, Higham Thm 5.6), not from the tolerance magnitude
(`FINDINGS §C11`, the `aic_sgn` convergence wall).

### Subtleties that bite here

- **$\Img\Phi=\Ker(1-\Phi)$ only when $\Phi$ is *exactly* idempotent**
  (`tex:344`). For $\eta$-idempotent $\Phi$ the identity $\calA=\Img\wt\Phi
  =\Ker(1-\wt\Phi)$ uses the *regularised* $\wt\Phi$ (`tex:392`), which is exactly
  idempotent — that is the whole reason for Stage 1.
- **$\calA=\Img\wt\Phi$ is an *oblique* image.** $\wt\Phi$ is
  Hermiticity-preserving but **not** Hilbert–Schmidt self-adjoint, so the
  Frobenius-orthogonal projector onto $\calA$ does *not* respect the star
  structure (it leaves a defect $\|P\star P-P\|=O(1)$). Project into $\calA$
  with $\wt\Phi$ itself, available as $W\star I=\wt\Phi(W)$ (`FINDINGS §C3`).
- **An oblique idempotent's range is spanned by its *left* singular vectors.**
  Nonzero singular values are $\ge 1$; the top-rank left singular vectors span
  $\Img E$, the right ones span $\Img E^\dag$ — a different space (`FINDINGS §C4`).
- **No functional calculus *inside* $\calA$.** The $\sgn$/`prop_P` calculus does
  **not** generalize to the $\eps$-$C^*$ algebra (`rem_X2`, `tex:628`): even
  $Y^2=I$ may have no solution near a target when associativity holds only up to
  $\eps$. Do the calculus *ambient* in $\Ma{n}$ (an exact $C^*$ algebra), then
  project the result into $\calA$ with $\wt\Phi$ (`FINDINGS §C1`).
- **The star is $\wt\Phi(XY)$, not $XY$** — the single most common test-blindness
  in the whole project. The $\eta=0$ identity-channel oracle is *structurally
  blind* to it (there $\wt\Phi=\mathrm{id}$, so $\star\equiv$ plain); every star
  test needs a genuinely-oblique $\eta>0$ fixture with a non-vacuity mutation
  (`FINDINGS §C2`).

---

## Stage 2 — the $\eps$-$C^*$ algebra ($\eps=O(\eta)$)

### (a) Statement: the approximate axioms

An **$\eps$-Banach algebra** is a Banach space with bilinear multiplication
satisfying (`tex:407`–`tex:413`):

$$\|XY\|\le(1+\eps)\,\|X\|\,\|Y\|\quad(\text{ax\_prodnorm, }\texttt{tex:410}),
\qquad
\|(XY)Z-X(YZ)\|\le\eps\,\|X\|\,\|Y\|\,\|Z\|\quad(\text{ax\_assoc, }\texttt{tex:412}).$$

A **$*\eps$-Banach algebra** adds a conjugate-linear involution with the
**exact** properties $\|X^\dag\|=\|X\|$ and $(XY)^\dag=Y^\dag X^\dag$
(`tex:420`–`tex:423`). An **$\eps$-$C^*$ algebra** adds the approximate $C^*$
identity (`tex:425`–`tex:428`):

$$\|X^\dag X\|\ge(1-\eps)\,\|X\|^2\qquad(\text{ax\_C}^*,\ \texttt{tex:427}),$$

the upper bound $\|X^\dag X\|\le(1+\eps)\|X\|^2$ following from ax\_prodnorm and
ax\_$*$. The **unit** is, by default, only approximate (`tex:430`–`tex:434`,
ax\_eps\_unit):

$$\|XI-X\|\le\eps\|X\|,\quad \|IX-X\|\le\eps\|X\|,\quad \bigl|\|I\|-1\bigr|\le\eps,$$

with an *exact* unit (`tex:435`, ax\_exact\_unit) available **only after**
the $O(\eps)$-change of both unit and multiplication of `prop_unit` (`tex:441`,
`tex:672`).

**Theorem `th_almost_idemp`** (`tex:2192`): the space $\calA=\Img\wt\Phi$ with
the inherited norm/involution/unit and the star product is an *extended*
$O(\eta)$-$C^*$ algebra. "Extended" means the same holds uniformly for
$\Ma{n}\otimes\calA$ (`tex:2209`) — the cb-norm hypothesis from Stage 0 is
exactly what makes the ampliated statement go through. So **$\eps=O(\eta)$**
(`tex:398`).

A separate small parameter governs maps. A **$\delta$-homomorphism**
$v\colon\calA'\to\calA''$ (`tex:443`–`tex:449`) satisfies

$$\|v(I)-I\|\le\delta\quad(\text{hom\_unit}),\qquad
\|v(XY)-v(X)v(Y)\|\le\delta\,\|X\|\,\|Y\|\quad(\text{hom\_mult}).$$

A **$\delta$-inclusion** additionally has $(1-\delta)\|X\|\le\|v(X)\|\le
(1+\delta)\|X\|$ (`tex:451`); a **$\delta$-isomorphism** is a bijective
$\delta$-inclusion (`tex:455`). The inverse of a $\delta$-isomorphism is a
$(\delta+O(\delta^2))$-isomorphism (`tex:458`).

### (b) Paper's technique

`th_almost_idemp` (`tex:2208`) reduces the multiplicative axioms to two
identities about $\Phi$ (`tex:2197`–`tex:2205`):

$$\Phi\bigl(\Phi(\Phi(X)\Phi(Y))\,\Phi(Z)\bigr)
=\Phi\bigl(\Phi(X)\Phi(Y)\Phi(Z)\bigr)+O(\eta)\|X\|\|Y\|\|Z\|,$$

and the mirror identity. These are proved in the finite-dim case via the Choi
representation $\Pi=VV^\dag$ (`tex:2239`–`tex:2272`). The non-multiplicative
axioms (norm, involution, $C^*$ identity) are immediate from
$\|\wt\Phi-\Phi\|_\cb\le O(\eta)$ and the Schwarz inequality
$\Phi(X^\dag)\Phi(X)\le\Phi(X^\dag X)$ (`tex:2233`).

### (c) Constructive route in the package

The package materializes $\calA$ by computing a Frobenius-orthonormal basis
$\{B_k\}$ of $\Img\wt\Phi$ (left singular vectors of the oblique idempotent,
`FINDINGS §C4`), and certifies $\eps$ directly as a **certified associator-defect
bracket**: a $d^3$ sweep of $\|(X\star Y)\star Z-X\star(Y\star Z)\|$ over the
basis (with the star, never the plain product). The defect $/\eta$ is asserted
bounded — this is the executable form of $\eps=O(\eta)$.

### Subtleties that bite here

- **The axioms hold only *up to* $\eps$.** Associativity, submultiplicativity,
  the $C^*$ identity, and the unit laws are all approximate (`tex:407`–`tex:439`).
  Never assume exact associativity or an exact unit in any routine on $\calA$.
- **The exact unit is not $1_\calH$ for a compressed subalgebra.** The unit of a
  corner/compressed subalgebra $S_P$ is $\wt P=\Co_P(P)$ (`tex:1082`), not $1_n$;
  a unit-law estimator hardcoded to $1_n$ reads $1.0$ for a genuine $C^*$
  subalgebra whose unit is $\wt P\neq 1_n$ (`FINDINGS §C7`, `FINDINGS §C11`).
- **$\eps$ and $\delta$ are distinct.** $\eps$ measures how far $\calA$ is from a
  $C^*$ algebra; $\delta$ measures how far a map $v$ is from multiplicative
  (`tex:443`–`tex:455`). The main proof juggles several
  ($\eps,\delta,\delta_0,\eta$); never silently identify them (callout in
  CLAUDE.md).
- **The $\delta$-inclusion lower bound is over the operator norm.** It is
  $\inf_{X\neq0}\|v(X)\|/\|X\|$, *not* the basis sweep $\min_i\|v(E_i)\|$: a $v$
  bounded below on every basis element can still collapse a general combination
  (`FINDINGS §C6`). The sound surrogate is $\sigma_{\min}$ of the coordinate
  matrix (the exact Frobenius unit-ball inf); the faithful operator-norm
  worst-case is `FINDINGS §C12`.

---

## Stage 3 — the rigidity theorem `th_main`

### (a) Statement

> **`th_main`** (`tex:460`). Any finite-dimensional $\eps$-$C^*$ algebra $\calA$
> is $O(\eps)$-isomorphic to some $C^*$ algebra $\calB$. The implicit constant in
> $O(\eps)$ **does not depend on $\calA$ or its dimensionality.**

$\calB$ is a genuine finite-dimensional $C^*$ algebra,
$\calB=\bigoplus_l M_{d_l}$ (`tex:257`; Wedderburn structure). The map
$v\colon\calB\to\calA$ is a $\delta$-isomorphism for $\delta=O(\eps)$, and the
constant is **universal**. The extended version `th_main_ext` (`tex:1538`) gives
a single $v$ with all $1_{\Ma{n}}\otimes v$ being $\delta$-isomorphisms,
$\delta$ depending only on $\eps$ (`tex:400`).

**Why this is surprising.** It says a *quantitative, dimension-free* rigidity
holds: an algebra that satisfies the $C^*$ axioms only to within $\eps$ is
within $O(\eps)$ of an honest one, with a constant that does not blow up as the
matrices get large. The naive route fails exactly because *its* constant grows
with $n$.

### Why the naive averaging / Haar route fails

The cohomological idea (`tex:464`–`tex:484`): the associativity defect
$h(X,Y,Z)=(XY)Z-X(YZ)$ is a 3-cocycle (`tex:470`); if it were an $O(\eps)$-exact
coboundary $h=\partial g$, the modified product $X\cdot Y=XY+g(X,Y)$ would be
$O(\eps^2)$-associative, and iteration would give exact associativity (`tex:478`).
For an exact $C^*$ algebra the trivialising **diagonal** is
$D=\int dU\,(U^\dag\otimes U)$ (Haar, `tex:484`), giving $g(X,Y)=\sum_j A_j
h(B_j,X,Y)$ (`tex:482`). But:

> "naive constructions of the Haar measure (or just the diagonal) in the
> $\eps$-associative setting have error bounds proportional to $n=\dim\calA$. So
> the outlined procedure works only if $\eps<cn^{-1}$." (`tex:484`)

That is **dimension-dependent**, the opposite of the theorem. Any test of
`th_main` must therefore assert the constant does **not** grow with $n$ — a
test that merely confirms an $O(\eps)$-isomorphism *exists* for fixed $n$ is a
test that cannot fail (`FINDINGS §D2`, `FINDINGS §C11`).

### (b) Paper's technique (proof deferred to §9, `tex:1414`)

A three-stage incremental construction (`tex:486`–`tex:490`; shard A, shard E):

1. **Stage 1 (`tex:1415`).** Find a maximal commutative $C^*$-subalgebra; show
   its projection basis is all one-dimensional. The first step
   $\calB_0=\CC\to\calB_1=\CC\oplus\CC$ needs a **nontrivial projection** in
   $\calA$ — Hermitian, $\|P^2-P\|\le\delta$, not close to $\pm I$. Its existence
   (`lem_nontriv_projection`, `tex:931`) is proved **non-constructively**: build
   the approximate unitary group $\calU\subseteq\calA$ as a $C^1$ manifold (§5,
   implicit function theorem), note $U=2P-I$ is a fixed point of the inversion
   $\sigma\colon U\mapsto U^{-1}$, and count fixed points with the
   **Lefschetz–Hopf theorem** and H-space homotopy (§6, `tex:486`).
2. **Stage 2 (`tex:1430`).** Group indices into equivalence classes
   ($\dim S_{P_j,P_k}=1$ within a class, $0$ across, by `lem_add_dim`,
   `tex:1363`) and build each matrix block $M_r$ one dimension at a time with
   `lem_extension` (`tex:1378`), reducing error after each step.
3. **Stage 3 (`tex:1443`).** Merge the per-class isomorphisms with
   `cor_merge_sum` (`tex:1352`), again reducing error.

The error-control engine is **`cor_improvement`** (`tex:1317`): a $\delta$-inclusion
of a $C^*$ algebra $\calB$ into $\calA$ can be improved to a $\delta_0$-inclusion
with $\delta_0=O(\eps)$ *independent of $\delta$* — proved cohomologically using
$\calB$'s own diagonal (which exists, because $\calB$ is exact), in the style of
[Joh88], simpler here by finite-dimensionality (`tex:490`).

### (c) Constructive route in the package

Each non-constructive ingredient is replaced by a finite-dimensional algorithm
that meets the same bound:

- **Nontrivial projection** (Lefschetz–Hopf $\to$ ambient spectral split,
  `FINDINGS §B1`). Use the reduction $P=\tfrac12(I+X)$ for a Hermitian
  near-involution $X\in\calA$. Get $X$ as the **ambient** sign of a gap-shifted
  Hermitian element: pick a non-scalar Hermitian $H\in\calA$ with the largest
  interior eigenvalue gap, form $X=\sgn(s(H-tI))$ in $\Ma{n}$ (eig-free
  funcalc), $P_{\mathrm{amb}}=\tfrac12(I+X)$, then project into $\calA$ via
  $\wt\Phi$ (the oblique-image route, `FINDINGS §C3`). This avoids the §5 unitary
  group entirely and stays in $\Ma{n}$ where the calculus is valid (`FINDINGS §C1`).
  The gap audition must be **unit-aware** — on a wrapper $S_P$ a unit-blind
  selection collapses to the wrapper unit $\wt P_m$ (`FINDINGS §C16`).
- **Haar diagonal** (`tex:484` $\to$ explicit Heisenberg–Weyl closed form,
  `FINDINGS §B2`). For an exact $C^*$ block $M_d$, $D=\sum_{jk}d^{-2}S_{jk}^\dag
  \otimes S_{jk}$ with $S_{jk}=X^jZ^k$ generalized Paulis — exact, eig-free. For
  a direct sum the literal Cartesian-product transcription is **non-central**
  (it leaves spurious cross-sector terms because the finite Pauli sum has nonzero
  first moment); the correct finite central diagonal is the cross-term-free
  embedded sum $D=\sum_l\sum_{jk}d_l^{-2}(\hat S^{(l)}_{jk})^\dag\otimes
  \hat S^{(l)}_{jk}$ with $\hat S$ supported on block $l$ and zero elsewhere
  (`FINDINGS §A2` — a flagged `.tex` discrepancy at `tex:1254`; the `.tex` formula
  is wrong for $m\ge2$ blocks).
- **`lem_invfun` / contraction** — *already constructive* (`tex:564`): the Banach
  fixed-point iteration $x_n=g_y(x_{n-1})$ is the algorithm, with the paper's
  $c<1$ giving the rate (THE MANDATE, CLAUDE.md).
- The merging / add-dim / extension lemmas are explicit finite matrix
  constructions; error reduction (`cor_improvement`) is a Newton contraction
  $\delta_{s+1}\le C(\delta_s^2+\eps)$ to the $O(\eps)$ floor.

### Subtleties that bite here

- **The constant must not grow with $n$.** The package measures
  $c=\mathrm{iso\_def}/\eta$ per instance and checks it is **bounded and has no
  upward trend** across a dimension sweep (the robust canary: abs-max plus
  no-trend, *not* a within-family $c_{hi}/c_{lo}$ ratio, which measures
  fixture-geometry spread, not dimension growth and reads RED on geometry
  outliers). Empirically $c$ is bounded with no trend; the **analytic universal
  constant $c_0$ remains OPEN** (`FINDINGS §D2`, `FINDINGS §C11`).
- **The $\Omega(1)$ spectral gap is not proven a priori.** The projection
  construction certifies the gap *per instance* and aborts fail-loud on a
  degenerate spectrum; the *universal* a-priori guarantee is exactly what the
  paper needs Lefschetz–Hopf for, and is **OPEN** (`FINDINGS §D1`, bead aic-3qv).
  Empirically a genuine $\ge2$-dimensional $\eps$-$C^*$ algebra always *delivers*
  a projection at an $O(1)$ gap (a gauge artifact of orthonormalization), so the
  small-gap regime does not arise at the algebra level — but this is measured,
  not proven.
- **A growing constant is a stop condition.** A bound that grows with dimension
  where the paper claims it is universal means the naive route of `tex:484` has
  crept in; escalate, do not patch (CLAUDE.md stop conditions).

---

## Stage 4 — the exact structure (`th_idemp_structure`) and the channel factorization (`th_factorization`)

### (a) The exact case $\eta=0$: `th_idemp_structure` (the cleanest ground truth)

> **`th_idemp_structure`** (statement `tex:318`, proof `tex:2055`). For an
> *idempotent* UCP map $\Phi$ there exist a subspace $\calM\subseteq\calH$
> (inclusion $J_\calM$), a $C^*$ algebra $\calA$, and UCP maps
> $\Delta\colon\calA\to\Bo(\calH)$ and $\Gamma\colon\Bo(\calM)\to\calA$ with
> $\Co_\calM\colon X\mapsto J_\calM^\dag X J_\calM$, such that
> $$\Gamma\Co_\calM\Delta=1_\calA,\qquad \Delta\Gamma\Co_\calM=\Phi
> \qquad(\texttt{tex:328}).$$
> Moreover $w=\Co_\calM\Delta\colon\calA\to\Bo(\calM)$ is a $C^*$-algebra
> homomorphism, and operators in $\Img\Delta$ are block-diagonal w.r.t.
> $\calH=\calM\oplus\calM^\perp$ (`tex:331`).

This is the **$\eta=0$ oracle**: the constructions reduce to the exact
Choi–Effros decomposition, a genuine $C^*$ algebra with **zero defect** — the
cleanest ground truth in the paper and the strongest cross-check (CLAUDE.md
cross-check ladder, item 3). Combined with `prop_hom_structure` (`tex:259`),
$\Img w$ is identified as $\bigoplus_j\Bo(\calL_j)$ with explicit isometries
$W_j$ (`tex:2093`), $\Delta$ as `tex:2099`, and $\Gamma$ as the general
conditional expectation $\Gamma_j(X)=\Tr_{\calE_j}(W_j X W_j^\dag(1_{\calL_j}
\otimes\gamma_j))$ (`prop_Gamma`, `tex:2106`).

**Paper's technique** (`tex:2055`): $\calM$ is the **carrier** of $\Phi$ (support
of $\Phi^*(\rho_0)$ for full-rank $\rho_0$, `lem_carrier` `tex:1724`); $\calA=
\Img\Phi$; $\Lambda\colon X\mapsto\Phi(J_\calM X J_\calM^\dag)$ factors as
$\Delta\Gamma$; idempotency gives $\Gamma\Co_\calM\Delta=1_\calA$ (`tex:2074`);
$\Psi=\Co_\calM\Lambda=w\Gamma$ is idempotent and $\Img w$ is a $C^*$-subalgebra
(via `lem_idemp`, `tex:1916`).

**Constructive route** (shard G): $\calM=\Img(\Tr_\calF$ of $\Img V)$ from the
Choi isometry; $\calA=\Ker(1-\Phi)$ as the eigenvalue-1 eigenspace of the
superoperator $\Phi$ on $\mathrm{vec}(\Bo(\calH))$ (certified eigenvalue
clustering distinguishes $\lambda=1$ from $\lambda<1$); block sizes $d_l$ from
the Wedderburn / `prop_hom_structure` decomposition (`tex:259`); the conditional
expectation density matrices $\gamma_j$ are **uniquely extractable** from the
stored $\Gamma$ by a data-driven solve — maximally mixed only for the *uniform*
noiseless-subsystem oracle, non-uniform in general (`FINDINGS §C20`).

### (b) The almost-idempotent case: `th_factorization`

> **`th_factorization`** (`tex:2730`). For $\Phi$ with $\|\Phi^2-\Phi\|_\cb\le
> \eta$, there are a finite-dimensional $C^*$ algebra $\calB$ and UCP maps
> $\Delta\colon\calB\to\Bo(\calH)$, $\Upsilon\colon\Bo(\calH)\to\calB$ with
> $$\|\Delta\Upsilon-\Phi\|_\cb\le O(\eta)\quad(\texttt{tex:2734}),
> \qquad
> \bigl\|\Upsilon_n(\Delta_n(X)\Delta_n(Y))-XY\bigr\|\le O(\eta)\|X\|\|Y\|
> \quad(\texttt{tex:2736}),$$
> the second implying $\|\Upsilon\Delta-1_\calB\|_\cb\le O(\eta)$ (`tex:2739`).

So $\Phi\approx\Delta\Upsilon$ factors *through* the genuine algebra $\calB$:
$\Upsilon$ **encodes** $\Bo(\calH)$ into $\calB$, $\Delta$ **decodes** $\calB$
back into $\Bo(\calH)$, with $\Delta\Upsilon\approx\wt\Phi$ and
$\Upsilon\Delta\approx 1_\calB$. The dual maps $\Dec=\Delta^*$ and
$\Enc=\Upsilon^*$ are the approximate decode/encode **channels** of states
(`tex:400`, `tex:2159`).

**Paper's technique** (`tex:2742`, an *outline* — the labelled proof block stops
before the CP-ization). Start from the non-UCP maps $\wt\Delta$ ($v$ then the
inclusion $\calA\hookrightarrow\Bo(\calH)$) and $\wt\Upsilon$ ($\wt\Phi$ with
target $\calA$, then $v^{-1}$), which exactly satisfy
$\wt\Delta\wt\Upsilon=\wt\Phi$, $\wt\Upsilon\wt\Delta=1_\calB$,
$\|\wt\Delta\|_\cb,\|\wt\Upsilon\|_\cb\le1+O(\eta)$ (`tex:2750`). Then
**CP-ize** them in cb-norm: $\Delta'$ via a per-block unitary 1-design
(Heisenberg–Weyl twirl, `tex:2771`–`tex:2801`), $\Upsilon'$ via `lem_RC` with
$C_j=d_{L_j}^{-1}\Tr_{L_j}(R_j)$ and a top-singular-vector $\xi_j$
(`tex:2840`–`tex:2871`), unitalized as
$\Upsilon\colon X\mapsto(\Upsilon'(1))^{-1/2}\Upsilon'(X)(\Upsilon'(1))^{-1/2}$
(`tex:2897`).

**Constructive route** (shard H; `FINDINGS §C13`, §C14, §D4). Every object in the
outline is an explicit finite-dim matrix expression — partial traces, SVDs,
Heisenberg–Weyl twirls — so the gap is *closure*, not *constructivity*
(`FINDINGS §D4`, BUILDABLE-MODULO). Two findings are load-bearing:

- **The ancilla side swaps** ($1_\calF$ on the left, not the right) because the
  package builds $V$ ancilla-major, opposite to the paper's $H$-left layout. The
  F-left form gives $\|\Upsilon'\Delta-1_\calB\|=O(\eta)$; F-right gives $O(1)$
  (`FINDINGS §C13(a)`).
- **$\Delta'$ is only $O(\eta)$-completely-positive, not exactly CP.** The
  paper's manifest-PSD argument (`tex:2786`–`tex:2796`) needs $\wt\Delta=v$ to be
  an *exact* unital homomorphism; $v$ is only an extended $O(\eta)$-isomorphism
  (the star $\neq$ ordinary product), so the per-block Choi carries a small
  $O(\eta^2)$ negative eigenvalue at $m\ge2$, $\eta>0$. The package projects onto
  the CP cone with a tolerance-aware Choi$\to$Kraus extraction that clips the
  $O(\eta^2)$ defect but still aborts on a genuine $O(1)$ or $O(\eta)$ negative;
  this stays within the $\|\Delta-\wt\Delta\|_\cb\le O(\eta)$ spec
  (`FINDINGS §C14`). Consequently the `decode` channel is only
  $O(\eta)$-trace-preserving for $\eta>0$ — the rigorous certificate is the
  cb-norm round-trip bracket, not `iscptp` at machine tolerance (index.md "Known
  limits").

### Subtleties that bite here

- **The round-trip target is the conditional expectation $P_B$, not the full
  identity $1_\calB$.** $\Delta$ reads only $\calB$'s block-diagonal coordinates,
  so $\Upsilon\Delta=P_B$ (zero on off-block units). Comparing
  $\mathrm{Dec}\circ\mathrm{Enc}$ to the full $1_{M_{n_B}}$ reads $1.0$ on every
  off-block unit even at $\eta=0$ — a test that cannot pass; sweep only in-block
  units (`FINDINGS §C14(b)`, §C22(b)).
- **`factorize` has a tighter domain than the rest of the pipeline.** $\Upsilon$
  is built by a power-series functional calculus whose convergence domain is much
  smaller than the `prop_P` basin $\rho(\Phi^2-\Phi)<1/4$; the package pre-checks
  a conservative $\rho<0.10$ and throws a clean error rather than aborting
  (index.md). `associated_algebra` and `main_isomorphism` keep the wider
  $\rho<1/4$ domain; `certified_defect` is always safe at any $\eta$.
- **The composite $O(\eta)$ constant is OPEN.** $C=\|\Delta\Upsilon-\Phi\|_\cb/\eta$
  composes the `prop_P` constant, the iso constant $c_0$ (itself analytically open
  via `cor_improvement`, `FINDINGS §D2`), and the CP-ization constants, none
  multiplied out in the paper. Measured per-instance and asserted bounded +
  dimension-independent; the analytic value is a research bead chained after the
  $c_0$ work (`FINDINGS §D4`).

---

## Notation table

| Symbol | Meaning | `tex:` anchor |
|---|---|---|
| $\Phi$ | UCP map $\Bo(\calH)\to\Bo(\calH)$ on observables (Heisenberg picture), $\Phi(X)=\sum_a K_a^\dag X K_a$, $\Phi(1)=1$ | `tex:316`, `tex:1568` |
| $\eta$ | idempotency defect: $\|\Phi^2-\Phi\|_\cb\le\eta$; require $\eta<1/4$ | `tex:354` |
| $\wt\Phi$ | regularised idempotent $\theta(2\Phi-1)$; exactly idempotent, $\|\wt\Phi-\Phi\|_\cb\le O(\eta)$, not always CP | `tex:356`, `tex:2171` |
| $\|\cdot\|_\cb$ | completely bounded norm $\sup_n\sup_{X\neq0}\|(1_{\Ma{n}}\otimes\Lambda)X\|/\|X\|$; dual = diamond norm $\|\cdot\|_\Diamond$ | `tex:347` |
| $\rho$ | Gelfand spectral radius (eig-free); the package's basin guard $\rho(\Phi^2-\Phi)<1/4$ | `tex:520`, `FINDINGS §C15` |
| $\sgn,\theta$ | $\sgn(X)=X(X^2)^{-1/2}$; $\theta(X)=\tfrac12(I+\sgn X)$; need $\|X^2-I\|<1$ | `tex:514`, `tex:528` |
| $\star$ | Choi–Effros product $X\star Y=\wt\Phi(XY)$ (not plain $XY$) | `tex:341`, `tex:2188` |
| $\calA$ | the $\eps$-$C^*$ algebra $\Img\wt\Phi=\Ker(1-\wt\Phi)$, oblique image | `tex:391`, `tex:2185` |
| $\eps$ | how far $\calA$ is from a $C^*$ algebra; the $\eps$-$C^*$ axioms; $\eps=O(\eta)$ | `tex:407`, `tex:398` |
| $\delta$ | how far a map $v$ is from multiplicative ($\delta$-hom / -inclusion / -iso) | `tex:443`, `tex:451`, `tex:455` |
| $\delta_0$ | the fixed reduced inclusion defect, $\delta_0=O(\eps)$, independent of $\delta$ (error reduction) | `tex:490`, `tex:1317` |
| $\calB$ | the genuine finite-dim $C^*$ algebra $\bigoplus_l M_{d_l}$ that $\calA$ is $O(\eps)$-isomorphic to | `tex:257`, `tex:461` |
| $M_{d_l}$ | matrix block $\Bo(\calL_l)$, $\dim\calL_l=d_l$, in $\calB=\bigoplus_l M_{d_l}$ | `tex:257`, `tex:2771` |
| $v$ | the $O(\eps)$-isomorphism $\calB\to\calA$ from `th_main` | `tex:400`, `tex:460` |
| $\Co_\calM$ | compression $X\mapsto J_\calM^\dag X J_\calM:\Bo(\calH)\to\Bo(\calM)$ | `tex:325` |
| $\Delta$ | decode (observable side) $\calB\to\Bo(\calH)$, UCP; dual $\Dec=\Delta^*$ | `tex:319`, `tex:2730` |
| $\Upsilon$ | encode (observable side) $\Bo(\calH)\to\calB$, UCP; dual $\Enc=\Upsilon^*$ | `tex:400`, `tex:2730`, `tex:2897` |
| $\Gamma$ | UCP $\Bo(\calM)\to\calA$ (conditional expectation) in the exact case | `tex:319`, `tex:2106` |

---

## Subtleties a reader MUST internalise

1. **"$\Phi$ is $\eta$-idempotent" is a cb-norm statement, not an operator-norm
   one.** $\|\Phi^2-\Phi\|_\cb\le\eta$ ($\sup$ over all ampliations) — *not*
   $\|\Phi^2-\Phi\|_{op}\le\eta$. Estimating it needs the Choi matrix; the
   spectral radius $\rho$, the operator norm, and $\eta_\cb$ are three distinct
   numbers that coincide only when the defect is normal (`tex:347`,
   `FINDINGS §C17`).

2. **The product on $\calA$ is $X\star Y=\wt\Phi(XY)$, not $XY$.** Plain
   multiplication leaves $\Img\wt\Phi$. The $\eta=0$ identity oracle is blind to
   this ($\wt\Phi=\mathrm{id}$), so every star claim must be checked on an
   oblique $\eta>0$ fixture (`tex:341`, `FINDINGS §C2`).

3. **The $\eps$-$C^*$ axioms hold only up to $\eps$, and there is no exact unit
   by default.** Associativity, submultiplicativity, the $C^*$ identity, the unit
   laws are all approximate; the exact unit appears only after the $O(\eps)$-change
   of `prop_unit`, and a compressed subalgebra's unit is $\wt P$, not $1_\calH$
   (`tex:407`–`tex:441`, `FINDINGS §C7`).

4. **The isomorphism is $O(\eps)$, not exact — with a universal, dimension-free
   constant.** Tests assert the bound, *and* that the constant does not grow with
   $n$. The naive Haar/cohomology route is wrong precisely because its error
   $\propto n$ (`tex:484`); a constant that grows with dimension is a stop
   condition (`tex:460`, `FINDINGS §D2`).

5. **$\Img\Phi=\Ker(1-\Phi)$ only when $\Phi$ is *exactly* idempotent.** For the
   $\eta$-idempotent $\Phi$ this identity uses the *regularised* $\wt\Phi$
   ($\calA=\Img\wt\Phi=\Ker(1-\wt\Phi)$); do not shortcut the construction with
   the raw $\Phi$ (`tex:344`, `tex:392`).

6. **No functional calculus *inside* $\calA$ — do it ambient, then project.** The
   $\sgn$/`prop_P` calculus does not generalize to the $\eps$-associative star
   (`rem_X2`); compute in $\Ma{n}$ and project into $\calA$ with $\wt\Phi$, which
   is the *oblique* projector — the Frobenius-orthogonal projector leaves an
   $O(1)$ star defect (`tex:628`, `FINDINGS §C1`, §C3).

7. **Spectra are perturbation-fragile; the algorithms avoid naive eigenvalues.**
   For non-normal $X$, $\spec(X)$ moves like $t^{1/3}$ (`tex:540`); this is *why*
   the paper avoids holomorphic calculus and *why* the package prefers
   constructions on Hermitian/normal elements and certifies everything else with
   arb error balls (`tex:535`–`tex:544`).

8. **UCP (observables) and CPTP (states) are dual — mind the direction.** $\Phi$
   acts on $\Bo(\calH)$; the channel $T=\Phi^*$ acts on states; $\Enc,\Dec$ are
   duals of $\Gamma\Co_\calM$, $\Delta$. Getting the dual backwards silently
   transposes everything (`tex:277`, `tex:334`).

9. **The factorization round-trip returns $P_B$, not $1_\calB$, and `factorize`
   has a tighter domain.** $\Upsilon\Delta=P_B$ (the block-diagonal conditional
   expectation), zero on off-block units; `decode` is only $O(\eta)$-TP for
   $\eta>0$; and $\Upsilon$'s power-series domain is narrower than the
   $\rho<1/4$ basin (`tex:2739`, `FINDINGS §C14`, §D4).

10. **Two project-level questions remain OPEN.** Whether a dimension-independent
    spectral gap $\Omega(1)$ always exists for $\dim\calA>1$ (the role of
    Lefschetz–Hopf; `FINDINGS §D1`), and the explicit value of the universal
    constant $c_0$ from `cor_improvement` (`FINDINGS §D2`). The package certifies
    both *per instance* (fail-loud on degeneracy; measured-bounded constant), not
    a priori — consistent with the paper, not a proof of it.

---

*Ground truth: `paper/src/approximate_algebras.tex`. Living issue log:
`paper/FINDINGS.md`. Per-region maps: `paper/shards/shard-{A,B,C,D,E,F,G,H}-*.md`.*
