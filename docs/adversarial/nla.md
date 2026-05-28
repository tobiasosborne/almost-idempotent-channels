# Adversarial NLA Instance Catalog — almost-idempotent-channels

**Scope.** General numerical-linear-algebra adversaries for the AIC operation set.
The operator-algebra / UCP / cb-norm domain-specific half is a separate catalog;
do not duplicate it here. This file is pure NLA: non-normality, precision, spectral
pathology, cancellation, ball blow-up.

**AIC ops referenced throughout** (from `MODULE_PLAN.md` and the implemented `src/`):
- `opnorm` — largest singular value (Gram eigenvalue, `aic_mat_opnorm`)
- `frob` — Frobenius norm
- `herm_eig` — Hermitian eigendecomposition, SIMPLE-spectrum certified path
  (`aic_mat_eig_hermitian`; requires isolated eigenvalues, aborts on degeneracy)
- `svd` — singular values
- `kron` — Kronecker product
- `ptrace` — partial trace
- `sgn` — matrix sign, Newton-Schulz / Denman-Beavers, eig-free (`aic_sgn`)
- `abs` — `|X| = (X^2)^{1/2}` (tex:514)
- `theta` — `theta(X) = (I + sgn(X))/2` (tex:527)
- `prop_P` — projection extractor `theta(2P-I)` (tex:524)
- `pow` — `x^alpha` Taylor series (`aic_funcalc_pow`)
- `contraction` — Banach fixed-point solver (`aic_contraction_solve`, Picard/Anderson)

**Notation.** `||.||` = operator norm unless specified. "Ball blow-up" = arb radius
grows until `arb_lt(bound, threshold)` can no longer be certified, causing an abort
that is correct behavior but may need precision escalation. "Silent garbage" = what
a `double` routine would return without alerting the caller.

---

## Family 1 — Jordan / Defective Blocks

### 1a. The Paper's Canonical t^{1/3} Example (tex:540–544)

**Math.** The 3x3 companion-style matrix
```
X(t) = [[0, 1, 0],
         [0, 0, 1],
         [t, 0, 0]]
```
has `spec(X) = {t^{1/3}, e^{2pi i/3} t^{1/3}, e^{-2pi i/3} t^{1/3}}` (tex:542–543).
At `t = 0` the matrix is nilpotent with a single eigenvalue at 0 (multiplicity 3);
for small `t > 0` the three eigenvalues burst out to radius `t^{1/3}`.

**Why it is evil.** A perturbation of size `t` in one entry shifts eigenvalues by
`O(t^{1/3})` — a 1/3-power sensitivity that is _not_ Lipschitz. The eigenvector
matrix becomes ill-conditioned as `t -> 0` (the eigenvectors coalesce into the
single Jordan chain direction `[1, 0, 0]^T`). Any algorithm that passes through
eigenvalue locations inherits this fragility.

**AIC ops attacked.**
- `herm_eig` (not applicable directly — not Hermitian — but the analogous Hermitian
  near-degenerate case is Family 4). The _lesson_ is why `herm_eig` is used only
  after confirming SIMPLE spectrum.
- `sgn(X(t))`: for small `t`, `X(t)^2` has a near-zero spectrum. `||X^2 - I||`
  is close to 1 (since `||X^2|| ~ t^{2/3} << 1`, so `||X^2 - I|| ~ 1`). The
  domain guard `aic_funcalc_int_assert_domain` triggers: `X(t)` for small `t` is
  **outside the validity domain** `||X^2 - I|| < 1` and must be rejected. Correct
  behavior: loud abort. A `double` implementation might proceed and produce
  garbage.
- `contraction` (lem_invfun, tex:564–591): if `f(x) = X^2 - I` and `V = I`, the
  contraction constant `c` near `X(t)` involves `||d_X(X^2)||` which is large when
  eigenvectors are ill-conditioned. The Picard iteration stalls or diverges.

**Generator recipe (acb_mat).**
```
n = 3
acb_mat_zero(A)
acb_mat_entry(A, 0, 1) = 1
acb_mat_entry(A, 1, 2) = 1
acb_mat_entry(A, 2, 0) = t      <- evil knob: t in [1e-1, 1e-9]
```
Evil knob: `t`. Sweep `t = 10^k` for `k = -1, -3, -5, -7, -9`.

**Expected correct behavior.** `aic_funcalc_int_assert_domain` aborts for small `t`
(the domain guard fires). For `t ~ 1`, the eigenvalues are O(1) and the domain
check passes; sgn and abs should converge. The routine must never return
plausible-but-wrong output. The certified CORRECT behavior is: abort with a domain
violation message; after raising precision (which does not help here since the
geometry is the obstruction, not precision), the same abort fires.

**Why arb matters.** A `double` implementation computes `||X^2 - I|| ~ 0.9998`
and proceeds; the Newton-Schulz iteration nearly stalls and after 100 steps returns
a matrix that is not close to any involutory matrix. The arb ball for
`||X^2 - I||` will straddle 1 or exceed it, triggering the `arb_lt` guard that
refuses to proceed.

---

### 1b. Graded Jordan Blocks — n x n Nilpotent + Scalar

**Math.** For an n x n Jordan block `J_n(lambda)` with `lambda` small but nonzero:
```
J_n(lambda) = lambda * I_n + N_n,    N_n = superdiagonal of ones
```
The eigenvalue `lambda` has multiplicity `n`; for the sign function `||J^2 - I||`
needs to be `< 1`, requiring `|lambda^2 - 1| + O(n * |lambda|)` — but the
nilpotent part `N_n` inflates the power-series radius requirements. Specifically,
`||J_n(lambda)^2 - I|| = ||(2 lambda N + N^2 + (lambda^2-1)I)|| >= n-1` for
`lambda = 1` because `||N||_op = 1` but `||N^k||_op` decays slowly for
upper-triangular `N`.

**Why it is evil.**
- The eigenvector matrix of `J_n(lambda)` is the n x n Vandermonde-style matrix
  `V(lambda, e^{2pi i/n} lambda, ...)` with condition number `kappa(V)` growing
  like `n!` (known result for Jordan blocks). Any eig-based route is catastrophic.
- Newton-Schulz convergence is governed by `||I - X^2||`; for `X = J_n(1)`, that
  norm grows with `n` due to the off-diagonal nilpotent contributions, potentially
  violating the domain condition for all `n > n_crit(prec)`.
- `pow(X, alpha)` Taylor series: the binomial coefficients times `||N||^k` sum
  to `(1 + 1)^n = 2^n` in the worst case (Cayley-Hamilton is not exploited by the
  generic series), so the series ball radius explodes exponentially in `n`.

**AIC ops attacked.**
- `sgn`, `abs`, `theta`, `prop_P`: domain check; ball radius grows with `n`; at
  moderate `n` (~8–12) the arb ball precision needed to certify convergence
  exceeds practical limits.
- `herm_eig`: only Hermitian inputs are routed here, but near-degenerate Hermitian
  analogues (see Family 4) exhibit the same behavior.
- `pow`: Taylor series ball inflation — `O(2^n)` worst-case radius growth.

**Generator recipe.**
```
n x n, lambda in (0, 1)
acb_mat_zero(A)
for i in 0..n-1:  acb_mat_entry(A, i, i) = lambda
for i in 0..n-2:  acb_mat_entry(A, i, i+1) = 1
```
Evil knobs: `n` (3 to 16), `lambda` (1.0 down to 0.5). At `lambda = 1` this is
the Jordan block `J_n(1)`; `sgn(J_n(1))` should be `I` exactly, but the ball
radius from Newton-Schulz grows with `n`. For `lambda` slightly less than 1,
`||X^2 - I|| = ||(X-I)(X+I)|| ~ 2 |lambda - 1| * n` so the domain check requires
`|lambda - 1| < 1/(2n)`.

**Expected correct behavior.** For large `n` or `|lambda - 1|` too large: domain
guard fires (correct abort). For `lambda = 1` exactly and small `n`: `sgn(J_n(1))
= I`; Newton-Schulz converges in O(log prec) steps. For `n >= 8` at `prec = 64`:
ball radius likely explodes before convergence; the midpoint-convergence check
still fires but the certified ball radius exceeds tolerance — the solver returns
with a wide ball, which is the honest arb answer. A `double` routine returns `I`
silently with no indication that the iteration barely converged.

**Why arb matters.** The ball width is the diagnostic: doubles hide the radius
growth. The certified ball for `||sgn(X) - I||` at `n = 12, lambda = 1, prec = 64`
will be O(1), indicating the result is not trustworthy without precision escalation.

---

### 1c. Defective 2x2 — Minimal Jordan Block

**Math.** The simplest defective matrix:
```
J_2(eps) = [[eps, 1],
             [0, eps]]
```
with `eps` varying. At `eps = 0`: the only Jordan block with zero eigenvalue.
At `eps = 1`: `J_2(1)` has `sgn = I` (eigenvalue +1, degree 2 Jordan block).

**Why it is evil.** Minimal example catching the Jordan phenomenon with the
smallest matrix. For `eps` near 0: `||X^2 - I|| ~ 1 + 2*eps*eps + eps^4 ~ 1`, so
the domain guard is right at the boundary. The transition from "outside domain"
to "inside domain" happens at `eps ~ (1/2)^{1/2}` and is numerically treacherous
— the arb ball for `||X^2 - I||` straddles 1 for a range of `eps` values.

**AIC ops attacked.** `sgn`, `theta`, `prop_P` domain boundary test; the
`arb_lt(bound, one)` guard in `aic_funcalc_int_assert_domain` (see
`src/aic_funcalc_domain.c:60`).

**Generator recipe.**
```
A[0][0] = A[1][1] = eps;  A[0][1] = 1;  A[1][0] = 0
```
Evil knob: `eps` in `[0, 1]`. At `eps = 1/sqrt(2) - delta` for small `delta`:
`||X^2 - I||` just barely exceeds 1. At `eps = 1/sqrt(2) + delta`: just inside.

**Expected correct behavior.** Smooth transition: for `eps < 1/sqrt(2)`, abort
with domain violation; for `eps > 1/sqrt(2)`, proceed and converge. The arb guard
handles the boundary correctly (straddles => abort). A `double` routine would
proceed for all `eps > 0` and produce a sequence of plausible-looking iterates
that do not converge to any involutory matrix.

---

## Family 2 — Highly Non-Normal Matrices / Pseudospectra

### 2a. Davies-Type Pseudo-Bidiagonal (Large Departure from Normality)

**Math.** The n x n bidiagonal matrix
```
D_n(c) = diag(1, 2, ..., n) + c * N,    N = superdiagonal of ones
```
has eigenvalues `{1, 2, ..., n}` but departure from normality
`||D^* D - D D^*|| = O(c^2 n^2)`. The pseudospectrum `sigma_eps(D)` extends
far beyond the true spectrum for `eps << c n`. The eigenvector condition number
is `kappa(V) ~ (c*n)^n / n!` (stretched Vandermonde).

**Why it is evil.**
- `herm_eig` does not apply directly (not Hermitian), but the Hermitian variant
  `(D + D^*)/2` has the same eigenvalues at `{1,...,n}` but eigenvectors that are
  a Gram-Schmidt orthogonalization of the ill-conditioned columns of `V` — very
  sensitive to perturbation when `c*n` is large.
- `sgn(D_n(c))`: the eigenvalues `{1,...,n}` are all positive, so `sgn = I`
  exactly. But the Newton-Schulz iteration `Y_{k+1} = Y_k(3I - Y_k^2)/2` is a
  polynomial iteration; at each step the operator norm `||Y_k - I||` is
  amplified by the non-normality, causing the iteration to temporarily diverge
  before converging (the "hump" phenomenon). For `c*n` large, the hump can push
  `||Y_k^2 - I||` above 1 transiently, which the midpoint convergence check
  incorrectly reports as non-convergence if checked too early.
- `pow`: the Taylor series bound `||f(X) - f(x_0 I)|| <= sum |a_k| ||X - x_0 I||^k`
  uses the Banach-algebra operator norm, which for non-normal `X` is much larger
  than necessary; the certified bound is very loose.

**AIC ops attacked.** `sgn` (iteration hump and premature abort); `pow` (loose
certified bound); `opnorm` (pseudospectrum inflates Gram matrix eigenvalues).

**Generator recipe.**
```
n x n diagonal: A[i][i] = i+1  (i = 0..n-1)
upper-bidiagonal: A[i][i+1] = c
```
Evil knobs: `n` (4 to 16), `c` (0.1 to 10). For `c * n > 10`, the hump in
Newton-Schulz exceeds the `||X^2 - I|| < 1` domain precondition transiently.

**Expected correct behavior.** `sgn(D_n(c)) = I` for all `c` (eigenvalues are
positive real). The Newton-Schulz midpoint should eventually converge. The correct
behavior is convergence with a possibly large number of iterations; the 100-step
cap (`AIC_SGN_MAX_ITERS`) should not fire for moderate `n, c`. If it fires, it is
a false alarm caused by slow convergence due to non-normality, which is a bug:
the domain condition `||X^2 - I|| < 1` is satisfied (since `spec(D^2) = {1,4,...,n^2}`
are all >= 1 so `D^2 - I` has eigenvalues `{0,3,...,n^2-1}` — wait, this means
`||D_n(c)^2 - I||_op >= n^2 - 1` for n >= 2, so `D_n(c)` is **outside the domain**
unless the eigenvalues are near +-1. Correct behavior: domain check fires. The
evil aspect is that the eigenvectors make this obvious failure mode invisible to a
naive check.

**Why arb matters.** The arb ball for `||D_n(c)^2 - I||` will correctly be
`>= n^2 - 1 >> 1`, and the domain guard fires. A `double` norm computation might
underestimate the operator norm (by using an incomplete power iteration) and proceed
into a divergent Newton-Schulz.

---

### 2b. Grcar Matrix

**Math.** The n x n Grcar matrix `G_n`:
```
G_n[i][i] = 1 (diagonal), G_n[i][i+1] = G_n[i][i+2] = G_n[i][i+3] = 1 (upper),
G_n[i][i-1] = -1 (subdiagonal)
```
has eigenvalues distributed around a lens-shaped curve in the complex plane.
The pseudospectrum `sigma_eps(G_n)` for `eps ~ 10^{-10}` covers the entire lens
area, while the true eigenvalues are confined to the boundary. The departure from
normality `||G^* G - G G^*||` grows with `n`.

**Why it is evil.** For certifying `sgn(G_n)` (which would be applicable only if
`||G_n^2 - I|| < 1`, not generally true), the pseudospectrum means that
perturbations at the level of rounding error can move eigenvalues dramatically.
More directly: `opnorm(G_n)` is well-defined, but the SVD-based Gram matrix
`G_n^* G_n` has a condition number that grows with `n`, making the power-iteration
for the largest singular value slow to certify. The `aic_mat_opnorm` (which goes
through `aic_mat_herm_max_eig` on the Gram matrix) will have a wide arb ball for
the operator norm.

**AIC ops attacked.**
- `opnorm`: Gram matrix eigenvalue certification; wide balls.
- `herm_eig` (on `G^* G`): the condition number of the Gram matrix is
  `kappa(G_n)^2`; for large `n` the certified isolation of the largest eigenvalue
  from the second requires resolving a tiny gap.
- `svd`: closely related; the gap between largest and second-largest singular value
  determines how fast the power iteration certifies.

**Generator recipe.**
```
n x n (n = 8, 12, 16, 20)
A[i][i] = 1
A[i][j] = 1  for j in {i+1, i+2, i+3}, j < n
A[i][i-1] = -1  for i > 0
```
Evil knob: `n`. Severity is monotone in `n`.

**Expected correct behavior.** `opnorm` should return a certified upper bound.
For large `n`, the ball radius grows; at some `n`, the arb `global_enclosure` for
the Gram matrix's maximum eigenvalue has a radius exceeding the gap to the
second eigenvalue, and the `eig_simple` isolation fails. Correct behavior: fail
loud at the `acb_mat_eig_simple` step inside `aic_mat_int_eig_certified`
(`src/aic_mat_spectral.c:91`). A `double` routine returns a number silently.

---

### 2c. Random Triangular with Off-Diagonal Scaling (Departure from Normality)

**Math.** Upper triangular `T_n(alpha)` with diagonal `{i/n}` (i=1..n) and
off-diagonal entries `alpha * U[i,j]` where `U` is a fixed random upper triangular
matrix with entries in `[-1,1]`. The departure from normality is
`||T^* T - T T^*||_F = O(alpha^2 * n^2)`. The eigenvalues are exactly `{i/n}` and
are Lipschitz in the diagonal entries but exponentially sensitive in `alpha`.

**Why it is evil.** The condition number of the eigenvector matrix grows as
`kappa ~ exp(alpha * n)`. For `prop_P` applied to a near-idempotent `P` with this
structure (e.g., `P = T_n(alpha)` with diagonal near `{0,1}`), the Newton-Schulz
iteration for `sgn(2P - I)` is not eig-based and avoids the worst sensitivity,
but the certified ball radius in the arb path grows with `alpha * n` due to
accumulated rounding in the matrix multiplications.

**AIC ops attacked.** `prop_P`: the `delta = ||P^2 - P||_op` computation via
`aic_mat_opnorm` involves the Gram matrix of the defect, which has the same
eigenvector ill-conditioning. `sgn`: ball inflation in Newton-Schulz.

**Generator recipe.**
```
n x n upper triangular
A[i][i] = (i+1)/n  (eigenvalues near {0,...,1})
A[i][j] = alpha * fixed_rational[i][j]  for j > i
  (use fixed rationals to make the instance reproducible)
```
Evil knobs: `n` (4 to 16), `alpha` (0.01 to 5).

**Expected correct behavior.** For `alpha = 0`: exact diagonal, `sgn` = diagonal
sign, convergence in 1–2 steps. For large `alpha*n`: wide balls but correct
midpoints; at extreme `alpha`, the domain check may fire if off-diagonal entries
push `||X^2 - I||` above 1.

---

## Family 3 — Rank-Deficient and Near-Singular

### 3a. Near-Singular Input to sgn/|X| (Small Minimum Singular Value)

**Math.** The matrix
```
X(sigma) = diag(1, 1, ..., 1, sigma)
```
with `sigma -> 0`. Here `|X| = X` (diagonal, positive) and `sgn(X) = I`; but the
computation `(X^2)^{-1/2}` requires inverting the near-zero eigenvalue `sigma^2`,
and `(L_{|X|} + R_{|X|})^{-1}` (the operator in `d_X sgn(X)`, tex:623) has
smallest eigenvalue `2 * sigma -> 0`, making the derivative ill-conditioned.

**Why it is evil.**
- `abs` = `(X^2)^{1/2}`: the entry `sigma^2` near zero means the Taylor series for
  `x^{1/2}` around `x_0^2 = 1` has `||X^2 - I|| = |sigma^2 - 1| ~ 1`. Domain
  guard likely fires.
- Denman-Beavers iteration: at step 1, `Y_0^{-1}` involves inverting the near-zero
  entry `sigma`. `acb_mat_inv` for Y with a near-zero entry will fail to certify
  invertibility and abort (correct behavior: loud failure in
  `aic_sgn_denman_beavers`, line 167).
- Newton-Schulz: at step 1, `Y_1 = X(3I - X^2)/2`; the `sigma` entry gives
  `Y_1[n,n] = sigma(3 - sigma^2)/2 ~ 3*sigma/2`, still near zero. The iteration
  does not blow up (Newton-Schulz is inverse-free) but converges to `sgn` which is
  `I` and the `sigma` entry converges from `sigma -> 3sigma/2 -> ...` — this is
  actually convergence of the scalar Newton iteration for `sgn(sigma)` starting
  at `sigma`; for `sigma > 0` the iteration converges to 1, not 0. The concern is
  the arb ball width for the `sigma^{1/3}` intermediate iterates.
- `contraction` solver: if `V = L_{|X|} + R_{|X|}` and `|X|` has a near-zero
  eigenvalue, `||V^{-1}||` is large, making the contraction constant `c` close to
  1 or exceeding 1 — the `aic_contraction_assert_c` guard fires.

**AIC ops attacked.** `abs`, `sgn` (Denman-Beavers certifiability), `contraction`
(contraction constant blowup near `c -> 1`), `prop_P` (if the near-idempotent P
has near-zero singular values in `2P - I`).

**Generator recipe.**
```
n x n diagonal
A[i][i] = 1  for i < n-1
A[n-1][n-1] = sigma
```
Evil knob: `sigma` in `[1, 1e-9]`. Also graded: `A[i][i] = sigma^i` for
geometric spectrum.

**Expected correct behavior.** For `sigma` near 1: correct computation, tight
ball. For `sigma -> 0`: if `||X^2 - I||` exceeds 1 (happens at `sigma < sqrt(0) = 0`,
actually `||X^2 - I|| = 1 - sigma^2`, which is always `< 1` for `sigma > 0`), so
the domain guard does NOT fire. Newton-Schulz should converge (inverse-free, and
all entries go to `sgn(diagonal[i]) = {+1,...,+1}`). The arb ball for the
`sigma` entry converges slowly; at `sigma = 1e-9, prec = 64`, the ball for
`sgn(sigma)` entry is wide because the Newton-Schulz iterates for scalar
`x -> x(3-x^2)/2` starting at `1e-9` need many steps to converge (the scalar
Newton-Schulz for `sgn` on near-zero starting point converges but slowly:
`x_{k+1} ~ 3x_k/2` for small `x`). The 100-step cap may fire. This is a
legitimate precision limit; correct behavior: abort with the iteration cap message,
instructing the caller to use a different centering `x_0` for the power series.

**Why arb matters.** A `double` diagonal sign function returns `diag(1,...,1,1)` in
3 floating-point Newton-Schulz steps (the iteration from `sigma = 1e-9` converges
to 1 quickly once it escapes the near-zero regime). The arb ball is wide throughout
and correctly signals that the result is not certified.

---

### 3b. Near-Singular Gram Matrix (Attack on Operator Norm)

**Math.** `X(eps) = U diag(1, 1, ..., 1, eps) V^*` where `U, V` are fixed unitary
matrices (e.g., DFT matrices or products of Householder reflectors with rational
parameters). The Gram matrix `X^* X = V diag(1,...,1,eps^2) V^*` has smallest
eigenvalue `eps^2`.

**Why it is evil.** `aic_mat_opnorm` calls `aic_mat_herm_max_eig` on the Gram
matrix `X^* X`, which uses `acb_mat_eig_global_enclosure`. The global enclosure
radius `delta` from `acb_mat_eig_global_enclosure` depends on the resolvent norm
`||(zI - A)^{-1}||` along the contour; near the near-zero eigenvalue `eps^2`, this
is `O(eps^{-2})`. The certified ball on `lambda_max = 1` acquires a radius
`O(eps^{-2}) * (floating-point-error)`. For `eps = 1e-6` and `prec = 64`, this
can cause the ball on `opnorm` to be wider than 1, making any downstream comparison
unreliable.

**AIC ops attacked.** `opnorm` (certified ball inflation); `prop_P` (the `delta <
1/4` check in `aic_prop_P` uses `aic_mat_opnorm` on the defect).

**Generator recipe.**
```
Build U, V as products of k Householder reflectors with entries in rational
arithmetic (e.g., Hadamard factor for n=2,4,8,16):
A = U * diag(1,...,1,eps) * V^dagger
```
Evil knob: `eps` in `[1, 1e-8]`. For reproducibility, use
`U = DFT_n / sqrt(n)`, `V = DFT_n^* / sqrt(n)` (circulant structure).

**Expected correct behavior.** `opnorm` = 1 exactly; the certified ball radius
grows as `eps` shrinks. At small enough `eps`, the ball radius exceeds a user
threshold and the caller should raise precision. No incorrect answer should be
returned; the ball is always an overestimate of the true norm, so the computation
is correct but uninformatively wide.

---

### 3c. Near-Singular Inverse in Denman-Beavers

**Math.** Input `X = I + alpha * N` where `N` is strictly upper triangular with
`N^2 = 0` (2x2 block: `N = [[0,1],[0,0]]`). Then `X^2 = I + 2*alpha*N`,
`||X^2 - I|| = 2*alpha`, so for `alpha < 1/2` the domain check passes. The
Denman-Beavers step requires `Y_k^{-1}`; after step 1, `Y_1 = (I + alpha*N +
(I - alpha*N)/1) / 2` — wait, `Y_0^{-1} = I - alpha*N` (since `(I+alpha*N)^{-1} =
I - alpha*N` for `N^2=0`). So `Y_1 = (I + alpha*N + I - alpha*N)/2 = I`. DB
converges in exactly one step. This is the "easy" case; the evil variant is below.

**Math (evil variant).** `X = diag(eps, 1-eps)` for `eps -> 0`. Then `X^2 = diag(eps^2,
(1-eps)^2)`, `||X^2 - I|| = max(|eps^2-1|, |(1-eps)^2-1|) = 1 - eps^2 < 1`. DB
step: `Y_0^{-1} = diag(1/eps, 1/(1-eps))`; `Y_1 = (diag(eps, 1-eps) + diag(1/eps,
1/(1-eps)))/2`. The `(0,0)` entry is `(eps + 1/eps)/2 ~ 1/(2*eps)`, which grows
unboundedly. But `Y_1` is a legitimate intermediate iterate; the iteration
converges to `sgn(X) = diag(1, 1) = I`. The issue is that the intermediate
iterates have entries of size `O(1/eps)`, so the arb balls for those entries are
O(1/eps) wide, and the final ball for `sgn(X) - I` is wide even though the
midpoint is near 0.

**AIC ops attacked.** `sgn` via Denman-Beavers: large intermediate iterates blow
up the arb ball; Newton-Schulz is inverse-free and does not have this problem (this
is one reason Newton-Schulz was chosen as the default, per `aic_funcalc_sgn.c:43`).

**Generator recipe.**
```
2x2 diagonal: A[0][0] = eps, A[1][1] = 1-eps
```
Evil knob: `eps` in `[0.4, 1e-6]`.

**Expected correct behavior.** Newton-Schulz: converges without large intermediate
values, tight balls. Denman-Beavers: correct answer but wide ball for small `eps`;
the cross-check (Newton-Schulz vs DB, `tests/test_funcalc.c`) should flag large
disagreement in ball radii even when midpoints agree.

---

## Family 4 — Degenerate and Clustered Spectra

### 4a. Near-Degenerate Hermitian Spectrum (Attack on herm_eig)

**Math.** `H(gap) = diag(1, 1+gap, 2, 3, ..., n)` for `gap -> 0`. The pair
of nearly-equal eigenvalues `1` and `1 + gap` have a spectral gap `gap`.
`acb_mat_eig_simple` requires proving all eigenvalues are isolated; it must certify
that the balls around `1` and `1+gap` do not overlap. For small `gap`, this
requires `prec >= -log2(gap) + c` for some constant `c`.

**Why it is evil.** The direct attack on `aic_mat_eig_hermitian` and
`aic_mat_int_eig_certified`: `acb_mat_eig_simple` fails to isolate the spectrum
and aborts (correct behavior). But the issue is that `gap` might be estimated as
`O(machine_eps)` when it is actually `O(eps_problem)`, and the abort message does
not distinguish between "genuinely degenerate" (should fail loud) and "needs more
precision" (should escalate precision). The arb path needs to distinguish these.

**AIC ops attacked.**
- `herm_eig`: abort from `acb_mat_eig_simple` when gap is below precision limit.
  `src/aic_mat_spectral.c:91`: the `!ok` branch fires.
- `opnorm` (via `aic_mat_herm_max_eig`): uses `acb_mat_eig_global_enclosure` which
  does NOT need isolation; this path is robust to degeneracy. But the ball on
  `lambda_max` grows when two eigenvalues are nearly equal.
- `prop_P`: when `P` has nearly equal eigenvalues near 0 and near 1 (a spectral gap
  of `gap` in `{eigenvalues near 0, eigenvalues near 1}`), the stability of
  `theta(2P-I)` is controlled by the gap of `2P - I` near 0, which is `2*gap`.
  For `gap < delta/2`, the certified idempotent `Ptilde` has a ball width of
  `O(delta/gap)`.

**Generator recipe.**
```
n x n diagonal Hermitian
H[i][i] = i + 1  for i >= 1
H[0][0] = 1 + gap/2, H[1][1] = 1 - gap/2  (the near-degenerate pair)
```
Evil knob: `gap` in `[1, 1e-9]`. Conjugate by a random unitary (Householder
product with rational entries) to make it non-diagonal.

**Expected correct behavior.**
- `herm_eig` with `gap > 2^{-prec+c}`: succeeds, returns sorted eigenvalues with
  certified enclosures.
- `herm_eig` with `gap < 2^{-prec+c}`: `acb_mat_eig_simple` fails to isolate;
  `aic_mat_int_eig_certified` aborts with "clustered eigenvalues? raise prec".
  Correct behavior: fail loud. The caller must either raise precision or use the
  `eig_multiple` path (audition not yet implemented; see `aic_mat_spectral.c:13`).
- A `double` routine silently returns two eigenvalues numerically equal to the 7th
  decimal place with no indication of the isolation failure.

---

### 4b. Exact-Degeneracy Projector (Attack on the sgn Path)

**Math.** An orthogonal projector `P = diag(1, 1, 0, 0)` (or any rank-k projector
in an arbitrary basis). The spectrum is `{0, 1}` with multiplicities `k` and `n-k`.
`X = 2P - I = diag(1, 1, -1, -1)` has `X^2 = I` exactly (in exact arithmetic;
in arb, `||X^2 - I||_op` is certifiably 0 if computed from the exact integer
matrix).

**Why it is evil.** The `herm_eig` path `aic_mat_eig_hermitian` routes through
`acb_mat_eig_simple`, which REQUIRES simple (non-degenerate) eigenvalues and
ABORTS on a projector spectrum `{0 (mult k), 1 (mult n-k)}`. This is documented
in `aic_funcalc_sgn.c:14–16` as the primary reason Newton-Schulz was chosen:
the eig-based sgn is not auditionable on projector inputs. The adversarial
instance proves that the chosen algorithm (Newton-Schulz) correctly handles the
degenerate case while the eig-based alternative would abort.

**AIC ops attacked.**
- `herm_eig`: abort on degeneracy (this is the _expected_ abort, not a bug).
- `sgn` via Newton-Schulz: `Y_0 = 2P-I`, `Y_1 = (2P-I)(3I - (2P-I)^2)/2 =
  (2P-I)(3I - I)/2 = (2P-I)` — converges in 1 step (since `X^2 = I` already,
  Newton-Schulz is at the fixed point). Correct behavior: fast convergence.

**Generator recipe.**
```
n x n Hermitian projector, rank k
Method: P = Q[0:k, :].T @ Q[0:k, :]  where Q is an n x n unitary with
rational entries (e.g., Hadamard / DFT normalized).
```
Parameters: `n` in `{4, 6, 8}`, `k` in `{1, 2, ..., n-1}`.

**Expected correct behavior.** `sgn(2P-I) = 2P-I` exactly (or within arb
tolerance); `theta(2P-I) = P` exactly. `prop_P` with `delta = 0`: converges in
1 Newton-Schulz step; arb ball on `||sgn(X) - X||` certifiably 0.

---

### 4c. Nearly-Commuting Matrix Pairs (Attack on Associativity Defect)

**Math.** `A = diag(1, 1 + gap, 2)`, `B = A + gap * [[0,1,0],[0,0,1],[0,0,0]]`.
The matrices `A` and `B` nearly commute: `||AB - BA||_F = O(gap)`. For the
epsilon-C* defect computation `||(XY)Z - X(YZ)||`, we need
`||A(AB) - (AA)B|| = ||[A, A]B|| + ||A^2 B - A A B|| = 0 + ||A[A,B]||`. For
non-commuting `A, B`, this defect is `O(gap * ||A||^2)`, which equals the
precision floor when `gap = machine_eps`.

**AIC ops attacked.** `(XY)Z - X(YZ)` associativity defect computation: for
near-commuting pairs the defect is `O(gap)` and must be distinguished from the
`O(machine_eps)` floating-point noise. The arb ball on the defect includes both
the mathematical defect and the accumulated rounding; if these are comparable, the
certified defect is `O(gap + machine_eps)`, which is an honest answer but may be
too wide to satisfy the bound in the eps-C* axiom check.

**Generator recipe.**
```
n = 3, A = diag(1, 1+gap, 2)
B = A + gap * N  (N strictly upper triangular, entries 1)
X = A (x) B (x) A  (Kronecker products for cb-norm ampliation)
```
Evil knob: `gap` in `[1, 1e-10]`.

**Expected correct behavior.** For `gap >> machine_eps`: certified defect is
`O(gap)`, ball tightly bounding the mathematical value. For `gap ~ 1e-15 * ||A||^2`:
the mathematical defect is below floating-point noise; the arb ball at `prec = 53`
cannot resolve the defect from the rounding. Correct behavior: the certified ball
for the defect includes 0, meaning the bound `||defect||_op < eps` is satisfied
for any positive `eps`. A `double` routine reports a number with large relative
uncertainty.

---

## Family 5 — Ill-Conditioned, Graded, Wide-Dynamic-Range

### 5a. Hilbert Matrix

**Math.** `H[i][j] = 1/(i+j+1)` for `i,j = 0..n-1`. The condition number
`kappa(H_n) ~ exp(3.52 n)`. The smallest singular value `sigma_min ~ (1 + sqrt(2))^{-4n}`.

**Why it is evil.** Any operation involving `H_n^{-1}` suffers exponential
ill-conditioning. For AIC:
- `sgn(H_n)`: since `H_n` is Hermitian PSD, `sgn(H_n) = I` if all eigenvalues are
  positive (which they are). But `||H_n^2 - I||_op` is large (smallest eigenvalue
  of `H_n` is `~ kappa^{-1}` so `||H_n^2 - I|| ~ 1 - kappa^{-2} < 1` barely).
  Actually `||H_n^2 - I|| = ||H_n - I|| * ||H_n + I|| <= (1 + ||H_n||) * ||H_n - I||`.
  For `n >= 4`: `||H_n|| ~ 1.5` (from largest eigenvalue) but `||H_n - I||` may
  exceed 1 due to the eigenvalue distribution. Domain check fires for moderate `n`.
- `opnorm(H_n)`: the largest singular value of `H_n` is straightforward but the
  global enclosure on the Gram matrix has a wide ball due to the tiny gap between
  eigenvalues.

**AIC ops attacked.** `opnorm` (wide ball on Gram matrix max eig due to clustered
small eigenvalues); `sgn` (likely domain check fires for `n >= 5`); `pow`.

**Generator recipe.**
```
n x n (n = 4, 6, 8, 10)
A[i][j] = 1/(i+j+1)  (use exact rational acb entries: acb_set_fmpq)
```
Evil knob: `n`. The exact rational entries avoid any initial rounding, making
the test purely about algorithm behavior rather than input precision.

**Expected correct behavior.** `sgn(H_n)` should abort with domain violation for
`n >= 6` (to verify: compute `||H_n^2 - I||_op` symbolically and check). For
`n <= 4`, the operation should succeed. `opnorm` should return a certified ball;
the ball width is the test diagnostic.

---

### 5b. Kahan Matrix (Catastrophic Cancellation in Products)

**Math.** The n x n upper triangular Kahan matrix
```
K_n(c, s) with c^2 + s^2 = 1:
K[i][i] = s^i,  K[i][j] = -c * s^(i + j - 1)  for j > i
```
for `c = cos(theta)`, `s = sin(theta)`. It has `kappa_2(K_n) ~ (2c/s^{n-1})^n`
and is designed to defeat column-pivoting QR in LAPACK (the original Kahan pathology).
The product `K^T K` has diagonals that cancel catastrophically: the `(n,n)` entry
involves `s^{2(n-1)}`, which is `<< 1` for small `s`.

**Why it is evil.**
- `(XY)Z - X(YZ)` for `X = Y = Z = K_n(c, s)`: each triple product involves
  `s^{3k}` entries for various `k`; the difference of two nearly-equal products
  requires computing each to `3k * log2(1/s)` bits of precision to certify the
  difference is not just rounding noise. At `prec = 64` bits, for `n = 8` and
  `s = 0.9`, `k = 7`, we need `21 * log2(10/9) ~ 2.1` bits — OK. For `s = 0.01`,
  `k = 7`: `21 * log2(100) ~ 140` bits — arb at `prec = 64` loses all precision.
- `Phi^2 - Phi` computation: same cancellation pattern.

**AIC ops attacked.** Any computation involving `||AB - BA||` or `||A^2 - A||`
on Kahan matrices will suffer catastrophic cancellation. Specifically: `prop_P`
(the `delta = ||P^2 - P||_op` check where `P` is related to a Kahan-structured
near-projector).

**Generator recipe.**
```
n x n upper triangular Kahan matrix
c = cos(theta), s = sin(theta), c^2 + s^2 = 1
K[i][i] = s^i
K[i][j] = -c * s^(i+j-1)  for j > i
K[i][j] = 0  for j < i
```
Evil knobs: `n` (4 to 16), `theta` (0.01 to 1.5 radians). Small `theta`
(nearly orthogonal, `s ~ 1, c ~ theta`) is mild; large `theta` (near `pi/2`, `s ~ 0`)
is lethal.

**Expected correct behavior.** At `prec = 64`: for large `theta` (small `s`),
the arb balls for `(K^3)_{ij}` entries become O(1) (all precision lost). The
certified bound on `||(KK)K - K(KK)||_op` becomes `O(||K||^3)` — uninformative
but honest. Correct behavior: return the wide ball; a precision guard comparing
ball_radius / midpoint against a threshold should detect precision loss and request
escalation. A `double` routine returns a number near 0 silently (incorrect if the
true defect is 0 but the double lies about which direction the cancellation goes).

---

### 5c. Graded Diagonal (Wide Dynamic Range, Precision Consumption)

**Math.** `D(r) = diag(1, r, r^2, ..., r^{n-1})` for `r > 1`. The dynamic range
is `r^{n-1}`; the smallest eigenvalue is 1, the largest is `r^{n-1}`.

**Why it is evil.** For `r = 10` and `n = 20`, the diagonal entries span 19 orders
of magnitude. `acb_mat` operations on this matrix at `prec = 64` bits (53-bit
mantissa) lose all relative precision for entries smaller than `r^{n-1} / 2^{53}
= 10^{19} / 10^{16} = 1000` times machine precision — i.e., entries `r^k` for
`k < 19 - 53*log10(2) ~ -27` lose all bits. In practice at `prec = 64`, `r = 2`,
`n = 70`: same issue.

The direct attack: `sgn(D(r))` = `I` for all positive `r`. Newton-Schulz on
`D(r)` with `||D(r)^2 - I||_op = r^{2(n-1)} - 1` fails the domain check for any
`r > 1` and `n >= 2`. The correct recipe is `sgn(D(r)) = I` via centering at `x_0
= (r^{(n-1)/2})` — but the scalar centering does not generalize to this diagonal
matrix (the center must be a scalar multiple of `I`).

**AIC ops attacked.** `sgn`, `abs`, `pow` (domain check fires or Taylor series
diverges); `opnorm` (certified largest singular value `r^{n-1}` has a ball width
growing as `r^{n-1} * 2^{-prec}`).

**Generator recipe.**
```
n x n diagonal
D[i][i] = r^i  (use acb_set_d or exact arb_t power for certification)
```
Evil knobs: `n` (4 to 16), `r` (1.01 to 100).

**Expected correct behavior.** `sgn` domain check fires for `r > 1`, `n >= 2`.
`opnorm` returns a certified ball with width `O(r^{n-1} * 2^{-prec})`; the
relative ball width is always `O(2^{-prec})`, so the certified norm is accurate
regardless of `r` (this is a strength of arb). The test should verify that the
ball width divided by the midpoint stays bounded by `O(2^{-prec+k})` for small `k`.

---

## Family 6 — Classic Gallery Torture Matrices

### 6a. Companion Matrix with Clustered Roots

**Math.** The companion matrix of `p(x) = (x - lambda)^n - eps` for small `eps > 0`:
```
C_n(lambda, eps) = the standard companion matrix of
p(x) = x^n - n*lambda x^{n-1} + ... + ((-1)^n * lambda^n - eps)
```
Wait — more directly: companion matrix of `prod_{k=0}^{n-1}(x - (lambda + delta_k))
- eps` where `delta_k = eps^{1/n} * e^{2pi i k/n}` are the perturbation magnitudes.
The roots are `{lambda + eps^{1/n} e^{2pi i k/n}}`: a ring of `n` roots tightly
clustered near `lambda` with gap `|delta_k - delta_j| ~ eps^{1/n}`.

**Why it is evil.** Same t^{1/3} mechanism as Family 1a, generalized. The companion
matrix is highly non-normal (condition number `kappa(V)` grows as `n!`). The
certified eigenvalue isolation by `acb_mat_eig_simple` needs to separate roots at
mutual distance `O(eps^{1/n})`, requiring `prec >= n * log2(1/eps) + c`. At
`eps = 1e-9, n = 6`: `prec >= 6 * 30 + c = 180 + c` bits minimum.

**AIC ops attacked.** `herm_eig` (via the Gram matrix `C^* C` having tightly
clustered eigenvalues corresponding to the clustered roots); `opnorm` (wide ball);
`sgn` if the matrix is close to the sign domain.

**Generator recipe.**
```
n = 4, 5, 6
lambda = 1
Roots r_k = 1 + eps^{1/n} * exp(2*pi*i*k/n)  for k = 0..n-1
Companion matrix from these roots (construct polynomial coefficients via Vieta)
```
Evil knobs: `n` (3 to 7), `eps` (`1e-3` to `1e-12`).

**Expected correct behavior.** At `prec = 64`: `acb_mat_eig_simple` fails for
small `eps`. Correct behavior: abort from `aic_mat_int_eig_certified` with
"clustered eigenvalues? raise prec". The test should verify that at `prec = 200`
the eigenvalues are successfully isolated for `eps = 1e-9, n = 4`.

---

### 6b. Frank Matrix

**Math.** The n x n Frank matrix:
```
F_n[i][j] = min(i, j) + 1   (1-indexed)
```
or equivalently, `F_n[i][j] = n - max(i,j) + 1`. It has
eigenvalues that are positive and pair up: `lambda_k * lambda_{n+1-k} = 1`
approximately (actually exact for the "symmetric Frank matrix"). The condition
number is `kappa_2(F_n) ~ exp(2.7 n)`.

**Why it is evil.**
- `herm_eig` on `F_n` (Hermitian variant: `(F_n + F_n^T)/2`): the eigenvalues
  near 1 are clustered for large `n`; isolation fails without high precision.
- `opnorm(F_n)`: the Gram matrix `F^T F` has condition number `kappa(F)^2 ~
  exp(5.4 n)`. The global enclosure ball for `lambda_max` grows with this
  condition number.
- `sgn((F_n - lambda I) / mu)` for appropriate centering: if the eigenvalues of
  `F_n` span a range straddling 0 after centering, the sign function changes sign
  and is sensitive to eigenvalue location.

**Generator recipe.**
```
n x n (n = 6, 8, 10, 12)
A[i][j] = min(i, j) + 1  (0-indexed: min(i,j) + 1)
```
Evil knob: `n`.

**Expected correct behavior.** `opnorm` returns a certified value; ball width
diagnostic. `herm_eig` fails for `n >= 10` at `prec = 64` due to clustered
eigenvalues near 1 (verify by computing the minimum eigenvalue gap symbolically).

---

### 6c. Toeplitz Near-Singular

**Math.** The n x n Toeplitz matrix with entries `T[i][j] = c^{|i-j|}` for
`|c| < 1`. As `c -> 1`, this approaches the rank-1 all-ones matrix. The smallest
singular value `sigma_min(T) ~ (1 - |c|)`.

**Why it is evil.** `opnorm(T)`: the largest singular value is `n * (1-c^n)/(1-c)
~ n` for `c ~ 1 - 1/n`. `herm_eig` on `T^* T`: near-singular spectrum with
smallest eigenvalue `~ (1-c)^2 ~ 1/n^2`, tightly clustered near 0.

**Generator recipe.**
```
n x n Toeplitz with c[k] = c^k for k = 0..n-1
T[i][j] = c^{|i-j|}
```
Evil knobs: `n` (4 to 16), `c` (0.5 to `1 - 1/n`).

**Expected correct behavior.** `opnorm` correct with moderate ball width.
`herm_eig` on `T^* T`: fails to isolate eigenvalues near 0 for `c` near `1 - 1/n`
at `prec = 64`. Correct behavior: abort from eig isolation. For `pow(T, alpha)`:
if `T` has zero or near-zero eigenvalues, the domain check for `pow` at `x_0 =
lambda_max` fires when `||T - x_0 I||_op >= x_0` (the domain condition for the
fractional power series), which is always true for `c -> 1` since `T` has both
large and near-zero eigenvalues.

---

## Family 7 — Near-Boundary Inputs for Validity Domains

### 7a. Just-Inside / Just-Outside ||X^2 - I|| = 1

**Math.** Construct a sequence of inputs parametrized by `alpha` such that
`||X(alpha)^2 - I||_op = alpha`. The simplest: `X(alpha) = diag(sqrt(1+alpha),
1, ..., 1)` gives `||X^2 - I|| = alpha` exactly (the `(0,0)` entry of `X^2 - I`
is `(1+alpha) - 1 = alpha`, all others 0).

**Why it is evil.** Tests the `arb_lt(bound, one)` domain guard in
`aic_funcalc_int_assert_domain` (`src/aic_funcalc_domain.c:60`). At `alpha = 1 -
eps_arb` (just inside), the guard must pass. At `alpha = 1 + eps_arb` (just
outside), it must fail. At `alpha` in the arb ball straddling 1, the guard must
fail (conservative). The transition width is `O(2^{-prec})`.

**AIC ops attacked.** `sgn`, `abs`, `theta`, `prop_P`: all check the domain guard.
The boundary of the valid domain for `prop_P` is `delta = 1/4` (since `||X^2 -
I|| = 4 delta`, the prop_P domain boundary is `delta = 1/4 <=> ||X^2 - I|| = 1`).

**Generator recipe.**
```
Parametric: X(alpha) = diag(sqrt(1 + alpha - eps_shift), 1, ..., 1)
where eps_shift scans from -2^{-prec+2} to +2^{-prec+2}
```
This tests the guard at sub-ULP resolution. Use `acb_set_arb` with an explicit
arb ball for `sqrt(1+alpha)` to control the input uncertainty.

Evil knob: `eps_shift` in `[-2^{-52}, +2^{-52}]` for `prec = 64` (this is the
transition zone where floating-point cannot determine which side of the boundary
the input is on).

**Expected correct behavior.**
- `alpha` certified `< 1`: domain guard passes, `sgn` proceeds.
- `alpha` certified `> 1`: domain guard fails, abort.
- `alpha` in ball straddling 1 (`arb_lt` returns false): domain guard fails, abort.
  A `double` routine silently proceeds for `alpha < 1.0 + 1e-16` and can return
  the wrong answer for inputs in the transition zone.

---

### 7b. prop_P Delta Near 1/4

**Math.** `P(delta_P) = diag(delta_P, 0, ..., 0)` — a rank-1 near-projector with
`P^2 = diag(delta_P^2, ...) = diag(delta_P, 0, ...) = P` only if `delta_P = 0`
or `1`. For `delta_P = 1/2`: `||P^2 - P|| = ||diag(1/4, 0, ...)|| = 1/4` exactly.

The general construction: `P = diag(1/2 + t/2, 0, ..., 0)` gives `||P^2 - P|| =
|(1/2+t/2)^2 - (1/2+t/2)| = |(1/2+t/2)||(1/2+t/2) - 1|| = (1/2+t/2)(1/2-t/2) =
(1-t^2)/4`. At `t = 0`: `delta = 1/4` exactly (boundary). At `t -> 1`: `delta -> 0`.

**Why it is evil.** The `aic_prop_P` guard `arb_lt(delta, quarter)` (`src/
aic_funcalc_proj.c:75`) tests `delta < 1/4`. At `delta = 1/4` exactly (a rational
arb value), `arb_lt` returns false (since `1/4` is not `< 1/4`). Just inside:
passes. Just outside: fails. The transition zone is the boundary adversary.

**Generator recipe.**
```
P = diag(1/2 + t/2, 0, ..., 0)  for t in [0, 1]
t = 0 is the boundary (delta = 1/4 exactly)
t = 1e-6, 1e-9 are just inside
t = -1e-6 is just outside (delta > 1/4)
```
Use exact rational `acb_set_fmpq` to avoid input rounding contaminating the boundary.

**Expected correct behavior.** `t > 0` (delta < 1/4): `prop_P` succeeds. `t = 0`
(delta = 1/4): `arb_lt` returns false, `prop_P` aborts. `t < 0` (delta > 1/4):
abort. Verify the abort message cites tex:525.

---

### 7c. Contraction Constant c Near 1

**Math.** Construct a map `f` with Lipschitz constant `c = 1 - eps_c` for small
`eps_c`. The paper's bound `||x_n - x_{n-1}|| < r * c^{n-1}` (tex:591) means
convergence requires `n ~ 1/eps_c` iterations. At `c = 1 - 1e-4`, convergence
to `tol = 1e-10` requires `n ~ 4 * log(10) / 1e-4 ~ 92000` iterations — the
100-step cap fires.

**Why it is evil.** The `aic_contraction_picard` iteration with `max_iter = 100`
(or similar cap) will abort even though the mathematical iteration converges.
The evil is distinguishing "converging slowly" from "not converging at all."

**Generator recipe.**
```
f(x) = (1 - eps_c) * x  on acb_mat 1x1 matrices
V = I, x0 = 1, y = 0, r = 2, c = 1 - eps_c
```
Evil knob: `eps_c` in `[0.1, 1e-5]`.

**Expected correct behavior.** For `eps_c = 0.1` (c = 0.9): 100 iterations gives
`0.9^{100} ~ 2.7e-5` reduction; converges to `tol = 1e-10` in `~230` iterations —
exceeds 100-step cap. Correct behavior: abort. The test should verify the abort is
triggered and the error message names c and the tolerance. For `eps_c = 0.5`
(c = 0.5): converges in `33` iterations, passes. The Anderson accelerator
(`aic_contraction_anderson`) should handle near-1 `c` better; this is a
comparative benchmark between Picard and Anderson candidates.

---

## Family 8 — Arb-Precision Adversarial (Ball Explosion at Fixed Precision)

### 8a. Input Engineered to Saturate Ball Radius at prec=64

**Math.** Start with a matrix `X_0` with `||X_0^2 - I|| = 0.99` (just inside the
domain). After one Newton-Schulz step, `||Y_1^2 - I||` should decrease quadratically.
But if the off-diagonal entries of `X_0^2 - I` have large imaginary parts (from
accumulated rounding in `acb_mat_sqr`), the arb ball for `||X_0^2 - I||` can
overestimate the true norm.

Construct: `X_0` = a matrix whose entries are rationals with 30-bit numerators and
denominators, arranged so that `X_0^2 - I` has entries with 60-bit numerators
(requiring 60 bits to represent exactly), which at `prec = 64` is representable but
at `prec = 53` forces rounding, inflating the ball.

**Why it is evil.** The ball inflation may cause `arb_lt(||X^2-I||, 1)` to fail
even though the true value is `< 1`. The test verifies that at `prec = 64` the
guard passes while at `prec = 53` it (correctly) fails due to ball overlap.

**Generator recipe.**
```
X = matrix of rational entries with ~30-bit numerators/denominators
chosen so X^2 - I has entries requiring ~62 bits to represent exactly.
Concretely: let p = 2^31 - 1 (a prime); use X[i][j] = (a_ij / p) + i*(b_ij/p)*eps
for fixed integers a_ij, b_ij.
||X^2 - I||_op = 1 - 2^{-30} (just inside at infinite precision).
```
Evil knob: `prec` (sweep `prec = 53, 64, 128, 256`). At `prec = 53`:
guard fails (correct; ball straddles 1). At `prec = 64`: guard passes. At
`prec = 128`: guard passes with tight ball.

**Expected correct behavior.** At `prec = 53`: domain guard aborts (correct —
cannot certify). At `prec = 64`: domain guard passes, Newton-Schulz converges with
a ball of radius `O(2^{-64})`. At `prec = 128`: tighter ball. A `double` routine
would proceed at prec=53 with the wrong answer.

---

### 8b. Accumulating Rounding in Long Newton-Schulz Chains

**Math.** A matrix `X` with `||X^2 - I|| = 0.5` (moderately inside the domain)
but requiring `K = 20` Newton-Schulz steps to converge. Each step of Newton-Schulz
`Y_{k+1} = Y_k(3I - Y_k^2)/2` performs two `acb_mat_mul` operations, each adding
`O(n^2 * 2^{-prec})` to the ball radius (from `acb_mat_mul`'s Frobenius-bound
rounding error). After 20 steps, the accumulated radius is `O(20 * n^2 * 2^{-prec})`.

**Why it is evil.** For `n = 16, prec = 64`: accumulated radius `~ 20 * 256 *
2^{-64} ~ 1.1e-15` — negligible. For `n = 32, prec = 128`: `~ 20 * 1024 *
2^{-128} ~ 6e-37` — also fine. The evil case is `n = 64, prec = 64` with a slowly
converging input: `20 * 4096 * 2^{-64} ~ 4.4e-13` — the ball radius for the
final iterate is `~1e-12`, meaning only 12 digits of precision remain.

**Generator recipe.**
```
X = Q D Q^*  where D = diag(d_1, ..., d_n)
with |d_k^2 - 1| = 0.5 for all k (e.g., d_k = sqrt(1.5) for all k,
or d_k = sqrt(1.5) * exp(2*pi*i*k/n) with varying arguments)
Q = random unitary (product of Householder reflectors with rational entries)
```
Evil knobs: `n` (8 to 64), `prec` (64 to 512). The slow-convergence case is when
`||Y_k^2 - I||` decreases slowly due to the non-normal `Q`-conjugation.

**Expected correct behavior.** The midpoint of `sgn(X)` should be `Q sgn(D) Q^*`
to `O(2^{-prec})`. The ball radius should grow as `O(K * n^2 * 2^{-prec})` where
`K` is the iteration count. The test should assert `ball_radius <= C * K * n^2 *
2^{-prec}` for an explicit constant `C` (the arb multiplication constant). A
`double` routine returns a result whose error cannot be bounded.

---

### 8c. Precision Guard Trigger: Insufficient prec Forces Abort

**Math.** Construct a matrix `X` such that at `prec = 53` the computation of
`||X^2 - I||_op` via `aic_mat_opnorm` returns an arb ball `[0.9 +/- 0.2]` (straddles
1), while at `prec = 200` the ball is `[0.95 +/- 1e-50]` (clearly inside). The
matrix is designed so that double-precision Frobenius-bound inflation causes the
`arb_lt` check to fail at `prec = 53`.

**Why it is evil.** This directly tests the precision-guard semantics: the guard
must fire at `prec = 53` even though the mathematical answer is "inside the domain",
and must not fire at `prec = 200`. It also ensures the caller cannot "paper over"
the failure by catching the abort and retrying at higher precision — the test suite
must explicitly test the retry logic.

**Generator recipe.** Same as 8a but with a shallower precision requirement:
`X = I + 0.3 * A` where `A` is a random matrix with entries in `[-0.1, 0.1]`
(rational, 20-bit numerators). Then `||X^2 - I|| = ||2*0.3*A + 0.09*A^2||_op ~
0.6 * ||A||_op <= 0.6 * 0.1 * n` — for `n = 3`, `||A||_op ~ sqrt(n) * 0.1 ~ 0.17`,
so `||X^2 - I|| ~ 0.1`. The arb ball at `prec = 53`: `0.1 +/- 1e-14` — passes.
For a carefully chosen `A` near `||X^2 - I|| = 0.99`, the ball straddles at
`prec = 53`.

**Expected correct behavior.** The test asserts:
1. At `prec = 53`: `aic_funcalc_int_assert_domain` aborts.
2. At `prec = 200`: completes and returns a result with ball radius `< 2^{-190}`.
The `double` path at `prec = 53` is the fail-loud diagnostic; the arb path at
`prec = 200` is the certified result.

---

## Sweep Design

### Parameter Ranges

For systematic coverage, sweep the following knobs:

| Knob | Range | Step | Priority |
|------|-------|------|----------|
| Jordan/nilpotent size `n` | 2, 3, 4, 6, 8, 10, 12 | multiplicative | 1 |
| Eigenvalue near zero `sigma` | 1, 1e-2, 1e-4, 1e-6, 1e-8 | decade | 1 |
| Spectral gap `gap` | 1, 1e-1, 1e-3, 1e-6, 1e-9 | decade | 1 |
| Near-boundary `alpha = ||X^2-I||` | 0.5, 0.9, 0.99, 0.999, 1.0, 1.001 | near-1 | 2 |
| Non-normality `c * n` | 1, 3, 5, 10, 20 | linear | 2 |
| Hilbert/Frank size `n` | 4, 6, 8, 10, 12 | additive | 2 |
| Kahan angle `theta` | 0.05, 0.3, 0.7, 1.2 | — | 3 |
| Denman-Beavers off-diagonal `eps` | 0.4, 0.1, 0.01, 1e-4 | decade | 3 |
| Precision `prec` | 53, 64, 128, 256, 512 | doubling | 1 |
| t^{1/3} perturbation `t` | 1e-1, 1e-3, 1e-6, 1e-9 | decade | 1 |

From mild to lethal: the "evil knob" moves parameters toward their lethal
values. The transition from "passes" to "aborts" as the knob increases is itself
a test: verify the transition is monotone and that the abort message is informative.

### Top 12 Single Most Punishing Instances

Ranked by expected adversarial impact on the AIC implementation:

1. **4b Exact-Degeneracy Projector** — `P = diag(1,1,0,0)` (or larger). Directly
   exercises the fact that `aic_mat_eig_hermitian` must NOT be used on projector
   inputs; validates the Newton-Schulz eig-free path as the only viable route.
   Tests the correct algorithm selection documented in `aic_funcalc_sgn.c:14–16`.
   Any regression to an eig-based route fails here immediately.

2. **1a t^{1/3} with t = 1e-6** — The paper's own canonical fragility example
   (tex:540–544). Domain guard must fire. Tests the interface between `||X^2-I||`
   norm certification and the `arb_lt` guard. Direct test of tex:516 precondition.

3. **7a Boundary alpha = 1.0 +/- sub-ULP** — The most direct test of the
   `aic_funcalc_int_assert_domain` guard. Input precision set so that the ball
   straddles the boundary. At sub-ULP resolution, a `double` routine cannot
   distinguish in/out; arb must abort conservatively.

4. **4a Near-Degenerate Hermitian gap = 1e-9** — Direct attack on
   `acb_mat_eig_simple` inside `aic_mat_int_eig_certified`. The abort path
   ("clustered eigenvalues?") must fire and be legible. Also tests that the
   `acb_mat_eig_global_enclosure` path in `aic_mat_herm_max_eig` stays robust
   (it does not need isolation).

5. **7b prop_P delta = 1/4 exactly** — The boundary of `prop_P`'s validity.
   With exact rational inputs, `arb_lt(delta, 1/4)` must correctly return false.
   Catches off-by-one errors in the `<` vs `<=` guard logic.

6. **8a prec=53 ball inflation** — Validates the precision-ladder: the same
   mathematical input passes at `prec=200` and fails at `prec=53`. Tests the
   "raise prec" suggestion in the abort message. Critical for the double-vs-arb
   cross-check (the double path is incorrect here, arb is correct).

7. **3a Near-Singular sigma = 1e-9, n=4** — Newton-Schulz slow convergence for
   near-zero entry. Tests whether the 100-step cap fires (it should for very small
   `sigma`), verifying the cap is in the right place. Also tests that the
   Denman-Beavers path aborts loud on the `acb_mat_inv` failure.

8. **5b Kahan n=8, theta=1.2** — Catastrophic cancellation in `||P^2 - P||_op`
   computation. The arb ball for the defect becomes O(1); tests whether the
   precision guard in the opnorm computation detects precision loss and requests
   escalation. A `double` routine returns 0 silently (incorrect).

9. **4a herm_eig on random unitary conjugate of gap=1e-6 diagonal** — More
   realistic than the diagonal case (Family 4a, generator recipe step 2: conjugate
   by Householder). Tests the full `aic_mat_int_eig_certified` path including the
   heuristic `approx_eig_qr` seed. The heuristic may return approximate eigenvalues
   that are accidentally isolated; the certified step must still fail.

10. **1b Jordan block J_8(1)** — Newton-Schulz on an exact sign matrix (should
    converge in 1 step) but with large ball radius at step 0. Tests the midpoint
    convergence check vs. ball convergence check distinction (see `aic_funcalc_sgn.c:66–88`).
    Validates that the midpoint-based `step_norm` function correctly reports
    convergence even when the arb ball radius is large.

11. **7c Contraction c = 0.9** — Picard iteration hits the 100-step cap before
    convergence. Tests that the abort fires, that Anderson acceleration handles
    this case better (comparative), and that the error message is informative about
    the contraction constant.

12. **2a Davies-type D_n(c) for n=10, c=2** — Domain check fires (operator norm
    of `D^2 - I` is `n^2 - 1 = 99 >> 1`). Tests that large-opnorm inputs are
    rejected early rather than entering a divergent iteration. Also a regression
    test for operator norm computation on ill-conditioned Gram matrices.

---

## Prioritized AIC Operation Fragility Assessment

**Most fragile op: `sgn` (via `aic_funcalc_int_assert_domain`).**
The domain guard `||X^2 - I|| < 1` is the chokepoint: it is attacked by all of
Families 1, 2, 3, 5, 7. The guard is correct by design but depends on a certified
`opnorm` computation, which itself can fail (Family 2b, 4a); the guard is therefore
only as reliable as the Gram-matrix eigenvalue certification it rests on. For
non-normal or degenerate inputs the guard aborts conservatively, which is the
right behavior but means many valid inputs are rejected at low precision.

**Second most fragile: `herm_eig` / `aic_mat_int_eig_certified`.**
The simple-spectrum requirement is a hard wall: any Hermitian input with a spectral
gap below `O(2^{-prec})` aborts. The wall is correctly placed but the abort must be
distinguished from "truly degenerate" vs "needs more precision". The `eig_multiple`
fallback is not yet implemented (`aic_mat_spectral.c:13`).

**Third most fragile: `prop_P` delta guard.**
The `delta < 1/4` check in `aic_prop_P` (tex:525) is a derived condition; errors
in the `opnorm` of the defect `P^2 - P` propagate here. For ill-conditioned `P`
(Families 2c, 5a, 5b), the certified delta may be larger than the true delta,
causing a false abort.
