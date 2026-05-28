# Shard B — Analytic tools & eps-Banach basics (L493–689)

## Scope

§3 "Some analytic tools" (L493–630) and §4 "Basic properties of eps-Banach algebras" (L632–689). Covers: Taylor/holomorphic functional calculus setup and bounds; the workhorse operations |X| and sgn(X); projection extraction (prop_P); a warning on spectral fragility; the quantitative inverse function theorem (lem_invfun); an illustrative derivation of d/dX sgn(X); a critical caveat that prop_P does not lift to eps-Banach algebras (rem_X2); and exact-unit normalisation (prop_unit). Most of §5–9 depends on results in this shard.

---

## Results

### fc-Taylor — Unnumbered "Taylor/power-series functional calculus" (L503–511)

- **Statement (faithful):** Let B be a Banach algebra and f(x) = sum_{n>=0} a_n (x - x_0)^n with convergence radius >= r. Then for all X in B with ||X - x_0 I|| < r:
  - f(X) = a_0 I + sum_{n>=1} a_n (X - x_0 I)^n  (L506–507)
  - ||f(X) - f(x_0 I)|| <= sum_{n>=1} |a_n| ||X - x_0 I||^n  (L508–510)
  - f(X) commutes with X.
  - Special case: f(x) = x^alpha is well-defined if ||X - x_0 I|| < x_0 (x_0 > 0 real). (L511)
  - For inverses: X^{-1} = X_0^{-1}(I + (X - X_0)X_0^{-1})^{-1} exists if ||X - X_0|| < ||X_0^{-1}||^{-1}. (L511)

- **Proof:** Standard Banach-algebra power-series convergence; no proof given in text.

- **Constructive in finite dim?:** yes-already. Inputs: matrix X (n x n, fp entries), x_0, coefficients {a_n} up to truncation error. Operations: repeated matrix multiply and add; norm via operator norm (largest singular value). Truncate series when |a_n| ||X - x_0 I||^n < tolerance.

- **Depends on / Used by:** Used by fc-abs-sgn (directly), prop_P (via sgn), prop_unit (via lem_invfun), lem_gV (§5), prop_polar (§5). Basis for everything in §5–9 that manipulates matrix functions.

- **Arbitrary precision needed?:** Yes, if ||X - x_0 I|| is close to r (series near divergence boundary). In exact arithmetic the bound is tight; with floating point, catastrophic cancellation occurs when convergence is marginal. arb ball arithmetic gives certified radius of convergence checks.

- **Cross-check / oracle:** Compare against SciPy/LAPACK matrix power for small n; compare series truncation against eigendecomposition result for normal matrices.

- **Implementation notes / risks:**
  - FLINT/arb: `acb_mat_pow_series` or manual loop via `acb_mat_mul`.
  - Operator norm via `acb_mat_bound_inf_norm` or largest singular value (arb SVD).
  - Convergence check ||X - x_0 I|| < r must be certified (arb gives enclosures).
  - For non-normal X in finite dim, series may converge but eigenvalue-based alternative is faster; both routes should be compared.

---

### fc-abs-sgn — Unnumbered "|X|, sgn(X) via power series" (L513–522)

- **Statement (faithful):**
  - |X| = (X^2)^{1/2},  sgn(X) = X (X^2)^{-1/2},  valid when ||X^2 - x_0^2 I|| < x_0^2 (x_0 > 0). (L514–516)
  - sgn(X)^2 = I for all X in the domain.
  - If ||X^2 - I|| < 1 then ||sgn(X) - X|| <= ||X|| O(||X^2 - I||). (L519–521)
  - sgn(X) commutes with X and is the nearest exact solution of Y^2 = I to X (within the domain).

- **Proof:** Direct application of fc-Taylor with f(x) = x^{1/2} and f(x) = x^{-1/2} centred at x_0^2 I. The bound (L519–521) follows from ||Z^{1/2} - I|| <= 1 - (1 - ||Z - I||)^{1/2} applied to Z = X^2 / x_0^2. (This bound is also used later at L615.)

- **Constructive in finite dim?:** yes-already. Two routes:
  1. Power series: compute X^2, then apply sqrt and inverse-sqrt series centred at x_0^2 I.
  2. Newton/Denman-Beavers iteration for matrix square root: Z_{k+1} = (Z_k + Z_k^{-1})/2 converges quadratically to sgn(X) from Z_0 = X when X has no pure-imaginary eigenvalues. This is the preferred algorithm for large matrices.
  - Inputs: X (n x n matrix), x_0 (real scalar centring point), tolerance.
  - Outputs: |X| and sgn(X) as certified matrix enclosures.
  - Precondition check: certify ||X^2 - x_0^2 I|| < x_0^2 before calling.

- **Depends on / Used by:** Depends on fc-Taylor. Used by: prop_P (directly), d-sgn derivation (L594–626), rem_X2, prop_polar (§5), lem_nontriv_projection (§6), prop_H-group (§6). Critical workhorse throughout.

- **Arbitrary precision needed?:** Yes. Near ||X^2 - x_0^2 I|| ~ x_0^2 (boundary of domain) sgn(X) is nearly singular or nearly imaginary; the Newton iteration may stagnate. Also when X has near-zero or near-imaginary eigenvalues, |X| ~ 0 and sgn is ill-conditioned. arb certification essential to confirm domain membership.

- **Cross-check / oracle:** For normal X: eigendecompose, apply scalar sgn to eigenvalues, recombine. Match against Newton iteration result.

- **Implementation notes / risks:**
  - Denman-Beavers (DB) iteration: Z_{k+1} = (Z_k + Z_k^{-1})/2, W_{k+1} = (W_k + W_k^{-1})/2 with Z_0 = X, W_0 = I; Z converges to sgn(X). Each step: one `acb_mat_inv` + one `acb_mat_add` + scale.
  - DB requires ||X^2 - I|| < 1 or equivalently no purely imaginary eigenvalues. Monitor convergence via ||Z_{k+1} - Z_k||.
  - For |X|: once sgn(X) = S is known, |X| = S*X = X*S (both equal, since S commutes with X).
  - Risk: if X is nearly defective or has nearly-imaginary eigenvalues, convergence is slow. Switch to power-series if x_0 is known and domain condition is certified.
  - FLINT: `acb_mat_sqrtm` (if available) or manual DB loop.

---

### prop_P — Proposition "projection extraction from approximate projection" (statement L524–530, proof L532–533)

- **Statement (faithful):** Let P be an element of a Banach algebra with ||P^2 - P|| <= delta < 1/4. Define
  - P_tilde = theta(2P - I),  theta(X) = (I + sgn(X)) / 2.
  - Then: P_tilde^2 = P_tilde (exact idempotent), P_tilde commutes with P, and ||P_tilde - P|| <= ||2P - I|| O(delta).

- **Proof:** L532–533. One-line: set X = 2P - I. Then ||X^2 - I|| = ||(2P-I)^2 - I|| = ||4P^2 - 4P|| = 4||P^2 - P|| <= 4 delta < 1. Apply fc-abs-sgn bound (L519–521) to get ||sgn(X) - X|| <= ||X|| O(4 delta), then theta maps sgn(X) to an exact idempotent.

- **Constructive in finite dim?:** yes-already. Full algorithm:
  1. Compute X = 2P - I.
  2. Certify ||X^2 - I|| < 1 (fails if delta >= 1/4).
  3. Compute S = sgn(X) via Denman-Beavers or power series (fc-abs-sgn).
  4. Output P_tilde = (I + S) / 2.
  - Inputs: P (n x n matrix), delta (certified bound on ||P^2 - P||).
  - Outputs: P_tilde (exact idempotent), certified bound on ||P_tilde - P||.

- **Depends on / Used by:** Depends on fc-abs-sgn. Used by: lem_nontriv_projection (§6, L931), proof of th_main (needs nontrivial projections). Note: does NOT apply in eps-Banach setting — see rem_X2.

- **Arbitrary precision needed?:** Yes if delta is close to 1/4 (domain boundary). Near delta = 1/4 the sgn computation is near-singular. arb needed to certify the precondition and bound the output error.

- **Cross-check / oracle:** For normal P: eigendecompose, project eigenvalues to {0,1}, recombine. Check ||P_tilde - P|| <= ||2P-I|| * C * delta for explicit C.

- **Implementation notes / risks:**
  - The condition delta < 1/4 is strict; floating point may underestimate ||P^2 - P||. Use certified norm.
  - P_tilde is an exact algebraic idempotent in the Banach algebra (not merely approximate). In finite dim matrices, verify P_tilde^2 = P_tilde to machine precision after computation.
  - The commutativity P_tilde P = P P_tilde holds by construction (sgn(X) commutes with X = 2P - I, hence with P).

---

### spec-hfc — Unnumbered "Spectrum & holomorphic functional calculus — fragility warning" (L535–552)

- **Statement (faithful):**
  - Spectrum: spec_B(X) = {lambda in C : X - lambda I not invertible in B}. Nonempty closed subset of {|lambda| <= ||X||}. (L536–539)
  - Fragility example (L540–544): the 3x3 matrix X = [[0,1,0],[0,0,1],[t,0,0]] has spec(X) = {t^{1/3}, e^{2pi i/3} t^{1/3}, e^{-2pi i/3} t^{1/3}}. For small t, eigenvalues are O(t^{1/3}), extremely sensitive to perturbations of order t.
  - Holomorphic calculus: for f holomorphic on open S containing spec(X), f(X) = (1/2pi i) integral_C (zI - X)^{-1} f(z) dz, contour C encloses spec(X) in S. (L547–549)
  - Spectral mapping theorem: spec(f(X)) = f(spec(X)). (L550)
  - Taylor expansion around X_0 exists for ||X - X_0|| < delta where delta = min_{z in C} ||(zI - X_0)^{-1}||^{-1}, but this number is often hard to estimate. (L552)
  - **KEY WARNING (L535):** Kitaev explicitly states he has NOT used holomorphic functional calculus in the paper due to spectral fragility in the approximate setting and the absence of a *-representation.

- **Proof:** N/A (expository, references Kadison-Ringrose §3).

- **Constructive in finite dim?:** partial. The contour integral f(X) = (1/2pi i) int (zI-X)^{-1} f(z) dz is in principle computable via numerical quadrature, but the fragility warning is the central point: for non-normal matrices, eigenvalues can shift by O(t^{1/3}) under O(t) perturbations (the 3x3 example). This makes certified-precision holomorphic calculus expensive and unreliable as a general tool. The Taylor/Newton routes in fc-Taylor and fc-abs-sgn are preferred precisely because they avoid eigenvalue location.

- **Depends on / Used by:** This is a negative result / warning. It motivates the choice of power-series and Newton-iteration routes throughout. Relevant to: any step in §5–11 that could naively use eigendecomposition. The 3x3 example is the canonical counterexample for fragility arguments.

- **Arbitrary precision needed?:** Yes, strongly. The 3x3 fragility example is the motivating case for arb: a floating-point eigendecomposition of X with t = 1e-9 gives eigenvalues at O(1e-3) but cannot certify their exact location at O(eps^{1/3}) perturbation. arb ball arithmetic can certify eigenvalue enclosures only if gaps are certified.

- **Cross-check / oracle:** Verify the 3x3 eigenvalue formula symbolically: det(X - lambda I) = -lambda^3 + t = 0 => lambda = t^{1/3} * roots of unity. Check numerically for t = 1e-6.

- **Implementation notes / risks:**
  - Do not use LAPACK `dgeev`/`zgeev` for certified computations near non-normal matrices. Eigenvalues are not Lipschitz for non-normal matrices.
  - If holomorphic calculus is ever needed (e.g., for a function not accessible via Newton), use `acb_mat_eig` with certified enclosures and verify gap conditions.
  - The contour integral route requires: (a) certified eigenvalue enclosures, (b) verified gap between spec(X) and the contour, (c) certified quadrature (arb integration). All three are expensive.
  - For the algorithms in this project, stick to Taylor/Newton methods throughout §5–9.

---

### lem_invfun — Lemma "quantitative inverse function theorem" (statement L564–577, proof L579–592)

- **Statement (faithful):** Let V: A -> B be a linear isomorphism of Banach spaces. Let f: Ball_r(x_0) subset A -> B be C^1 such that ||V^{-1} d_x f(x) - 1|| <= c for all x in Ball_r(x_0), where 0 <= c < 1. Then:
  1. f is injective and (1-c)||x_1 - x_2|| <= ||V^{-1}(f(x_1) - f(x_2))|| <= (1+c)||x_1 - x_2||  for all x_1, x_2 in Ball_r(x_0). (L573–575)
  2. The image of f contains f(x_0) + V(Ball_{(1-c)r}(0)). (L576)

- **Proof:** L579–592. Banach contraction. For each target y, rewrite f(x) = y as x = g_y(x) where g_y(x) = x + V^{-1}(y - f(x)). Since ||d_x g_y|| <= c < 1, g_y is a contraction. Part 1 follows from ||a - b|| <= c||a|| where a = x_1 - x_2, b = V^{-1}(f(x_1) - f(x_2)). Part 2: for y with ||V^{-1}(y - f(x_0))|| < (1-c)r, show g_y maps Ball_r(x_0) to itself, then fixed-point iteration x_n = g_y(x_{n-1}) converges since ||x_n - x_{n-1}|| < r c^{n-1}. (L591)

- **Constructive in finite dim?:** yes-already — the proof IS the algorithm. The iteration x_n = g_y(x_{n-1}) = x_{n-1} + V^{-1}(y - f(x_{n-1})) is Newton-like and directly implementable.
  - Inputs: f (computable C^1 map on matrices), V (linear isomorphism, given by its matrix representation), x_0 (initial point), y (target), r (ball radius), c (Lipschitz constant bound, must certify c < 1).
  - Outputs: x* such that f(x*) = y, certified by ||x_n - x*|| <= r c^n / (1-c).
  - Operations: evaluate f, compute V^{-1} (matrix inverse of the operator V), iterate until ||x_n - x_{n-1}|| < tolerance.
  - Convergence: geometric, rate c < 1. If c << 1 (near-isometry), very fast.

- **Depends on / Used by:** Depends on: Banach space structure only (no algebra required). Used by: prop_unit (L682), lem_gV (§5, L765), prop_polar (§5), d-sgn derivation (L594–626 uses implicit function theorem which is derived from this). The iteration g_y is the core computational primitive for §5.

- **Arbitrary precision needed?:** Yes. The contraction constant c must be certified < 1. In floating point, c might appear < 1 but be >= 1 due to rounding. arb interval arithmetic certifies ||V^{-1} d_x f(x) - 1|| <= c_upper < 1. Also needed to certify the convergence radius (1-c)r.

- **Cross-check / oracle:** For f(x) = x^2 - x near x = I in M_n(C): V = identity, c = O(eps), x* = J (exact idempotent). Matches prop_unit. Verify convergence geometrically: ||x_n - x*|| ~ c^n.

- **Implementation notes / risks:**
  - V is typically chosen as the Frechet derivative df(x_0) itself, making this a standard Newton iteration. The lemma is more general: V can be an approximation to df(x_0) as long as the Lipschitz bound holds.
  - In finite dim, V^{-1} requires computing the inverse of an n^2 x n^2 linear operator (acting on Mat_n) or an n x n matrix inverse depending on the application. Be explicit about which Banach space A is.
  - For prop_unit: A = B = M_n(C) as Banach space, f(X) = X^2 - X, V = 1 (identity), g_y(x) = x + (y - (x^2 - x)) = x + y - x^2 + x.
  - Risk: if c is close to 1, convergence is slow; need many iterations. The paper's bound (1-c)r for the image radius degrades as c -> 1.

---

### d-sgn — Unnumbered "Derivative d/dX sgn(X) via implicit function theorem" (L594–626)

- **Statement (faithful):** Define f: B x B x B -> B x B by f(X, Y_1, Y_2) = (Y_1 Y_2 - I, X Y_2 - Y_1 X). The equation f(X, Y_1, Y_2) = 0 has solution Y_1 = Y_2 = sgn(X) when sgn(X) is well-defined. By the implicit function theorem, (Y_1, Y_2) = (sgn(X), sgn(X)) is a C^1 function of X when the partial Jacobian in (Y_1, Y_2) is invertible. That Jacobian is the block operator [[R_{Y_2}, L_{Y_1}], [-R_X, L_X]] where L_A, R_A are left/right multiplication operators. (L606) Setting Y_1 = Y_2 = sgn(X), all blocks commute and the inverse is [[L_X, -L_{sgn(X)}], [R_X, R_{sgn(X)}]] times (L_X R_{sgn(X)} + L_{sgn(X)} R_X)^{-1}. The denominator operator factors as L_{sgn(X)} (L_{|X|} + R_{|X|}) R_{sgn(X)}, which is invertible when ||X^2 - I|| < 1 (using ||Z^{1/2} - I|| <= 1 - (1 - ||Z - I||)^{1/2} at Z = X^2 implies |||X| - I|| < 1). (L615) Result:
  - d_X sgn(X) = (1 - L_{sgn(X)} R_{sgn(X)}) / (L_{|X|} + R_{|X|})  when ||X^2 - I|| < 1. (L623–625)

- **Proof:** L594–626 is a worked example of the implicit function theorem (derived from lem_invfun). No separate proof block; the derivation is self-contained in the text.

- **Constructive in finite dim?:** yes-but-proof-nonconstructive(route). The implicit function theorem route is nonconstructive, but the explicit formula at L623 is computable:
  - Inputs: X (n x n matrix), sgn(X) (from fc-abs-sgn), |X| (from fc-abs-sgn).
  - The operator (L_{|X|} + R_{|X|}) acts on M_n as the linear map Z -> |X| Z + Z |X| (Lyapunov-type operator). Its inverse can be computed via Sylvester/Lyapunov equation solver.
  - The formula d_X sgn(X) applied to a perturbation dX gives: d(sgn(X)) = (L_{|X|} + R_{|X|})^{-1} (dX - sgn(X) dX sgn(X)).
  - Operations: solve Lyapunov equation (|X| Z + Z |X| = RHS) for Z; multiply by (I - L_{sgn(X)} R_{sgn(X)}).

- **Depends on / Used by:** Depends on: lem_invfun (implicit function theorem), fc-abs-sgn. Used by: prop_polar (§5, gradient computation), lem_gV (§5). Important for any sensitivity analysis or Newton step involving sgn(X).

- **Arbitrary precision needed?:** Yes. The operator (L_{|X|} + R_{|X|}) is ill-conditioned when |X| has near-zero eigenvalues (i.e., X has near-zero singular values). arb needed to certify invertibility.

- **Cross-check / oracle:** For normal X = U diag(sigma_i) U*: d sgn(X) applied to dX = U dD U* (diagonal) gives d(sgn) = U diag(d sigma_i / |sigma_i|) U*. Check formula against scalar derivative sgn'(x) = 1/|x|.

- **Implementation notes / risks:**
  - Lyapunov equation A Z + Z A = C (with A = |X| symmetric PSD) has unique solution if spec(A) + spec(A) does not include 0, i.e., no zero eigenvalue of |X|. Use LAPACK `dtrsyl` or arb matrix equation solver.
  - Risk: |X| near-singular when X has near-zero singular values. The domain condition ||X^2 - I|| < 1 rules out ||X|| < (1 - 1)^{1/2} = 0 but not near-zero singular values of X.
  - This formula is needed only if Newton refinement of sgn is required; for the main algorithms, the direct DB iteration may suffice without explicitly computing the derivative.

---

### rem_X2 — Remark "prop_P does not extend to eps-Banach algebras" (L628–630)

- **Statement (faithful):** In an eps-Banach algebra, the system Y_1 Y_2 = I, X Y_2 = Y_1 X has a unique solution with Y_1, Y_2 close to X (when ||X^2 - I|| is small), but in general Y_1 != Y_2. Even the single equation Y^2 = I may have no solution close to X. Reason: if Y^2 = I and ||Y - X|| <= O(eps), then representing X = Y + Z with ||Z|| = O(eps) gives ||X X^2 - X^2 X|| <= O(eps^2). But such a commutativity bound need not hold in a general eps-Banach algebra. Therefore prop_P does NOT generalise to the eps-Banach or eps-C* setting; later sections must work with approximate projections throughout.

- **Proof:** L628–630. Short self-contained argument: the assumed closeness Y_1, Y_2 -> X forces XX^2 approx X^2 X up to O(eps^2), a constraint that eps-Banach algebras do not satisfy.

- **Constructive in finite dim?:** N/A (this is a negative/caveat result). The implication for implementation is: do not attempt to sharpen approximate projections via sgn in the eps-algebra setting; the algorithm would fail. All §6–9 code must carry projections as approximate (not exact).

- **Depends on / Used by:** Depends on: prop_P, fc-abs-sgn (uses the same infrastructure, identifies where it breaks down). Motivates: lem_nontriv_projection (§6), all of §6–9 approximate-projection machinery. Affects: any algorithm that tries to "clean up" approximate projections in an eps-algebra context.

- **Arbitrary precision needed?:** Not directly (negative result), but the argument shows that the O(eps^2) bound for ||XX^2 - X^2 X|| is the key discriminator; arb can certify whether a given matrix satisfies this or not.

- **Cross-check / oracle:** Construct an explicit eps-Banach algebra (e.g., matrix algebra with perturbed multiplication) and exhibit an X with ||X^2 - I|| = O(eps) but no Y near X with Y^2 = I exactly.

- **Implementation notes / risks:** CRITICAL DESIGN NOTE. Any code path that invokes prop_P must first verify it is operating in a genuine (exact) Banach algebra, not an eps-algebra. Add a runtime assertion or type-level distinction.

---

### La-Ra-bounds — Unnumbered "Left/right multiplication operators in eps-Banach algebras" (L634–670)

- **Statement (faithful):** For X in an eps-Banach algebra A:
  - ||L_X|| <= (1+eps)||X||,  ||R_X|| <= (1+eps)||X||.  (L644)
  - ||L_X|| >= (1 - O(eps))||X||,  ||R_X|| >= (1 - O(eps))||X||  (general). (L648–649)
  - ||L_X|| >= ||X||,  ||R_X|| >= ||X||  (exact unit). (L651–652)
  - ||L_X L_Y - L_{XY}|| <= eps||X||||Y||,  similarly for R. (L657–658)
  - ||L_X R_Y - R_Y L_X|| <= eps||X||||Y||. (L659)
  - Frechet derivative bound: ||d_X(X^2 - X) - 1|| <= O(delta + eps) for X in closed ball of radius delta around I. (L668–669, eq. d_X2X)

- **Proof:** L634–670. Direct from eps-Banach axioms (ax_prodnorm, ax_assoc). The derivative bound uses ||L_I - 1|| <= eps (since ||IY - Y|| <= eps||Y||) and the triangle inequality.

- **Constructive in finite dim?:** yes-already. These are explicit norm bounds. No iteration required. Given eps (certified from the algebra's definition), compute bounds from norms of X and Y.

- **Depends on / Used by:** Depends on eps-Banach axioms. Used directly by: prop_unit (needs d_X2X bound at L682), lem_U_delta (§5, bounds on L_X^{-1}), lem_gV (§5), prop_polar (§5). Fundamental to all of §4–9.

- **Arbitrary precision needed?:** Yes, in the sense that eps must be certified; if eps is only a floating-point estimate, the bounds (1 +/- O(eps)) may be wrong. arb tracks eps as a rigorous enclosure.

- **Cross-check / oracle:** For exact algebra (eps = 0), verify L_X L_Y = L_{XY} exactly. For constructed eps-algebra with known eps, verify bounds numerically.

- **Implementation notes / risks:**
  - In finite dim, L_X and R_X are n^2 x n^2 matrices acting on vec(A). Operator norms are largest singular values of these matrices (expensive for large n; use Kronecker-product structure).
  - The key quantity d_X2X: d_X(X^2 - X) = L_X + R_X - 1 (as operator on A). Its deviation from the identity is bounded by O(delta + eps). This is the precondition for lem_invfun in prop_unit.

---

### prop_unit — Proposition "exact unit via O(eps)-modification" (statement L672–679, proof L681–687)

- **Statement (faithful):** Let A be an eps-Banach algebra with unit I. Then there exist J in A and a new multiplication X . Y such that:
  - ||J - I|| <= O(eps),  J^2 = J (exact idempotent). (L675)
  - ||X . Y - XY|| <= O(eps)||X||||Y||  for all X, Y. (L675–676)
  - (A, ., J) is an O(eps)-Banach algebra with exact unit J.
  - If A has an involution, then J^dag = J and (X . Y)^dag = Y^dag . X^dag. (L678)
  - Explicit new multiplication: X . Y = R_J^{-1}(X) L_J^{-1}(Y). (L684)

- **Proof:** L681–687. Uses d_X2X bound (eq. d_X2X, L669) and lem_invfun: for delta > 2 eps, the map f: X -> X^2 - X satisfies ||V^{-1} d_x f(x) - 1|| <= O(delta + eps) < 1 in Ball_delta(I) with V = 1. By lem_invfun part 2, the image of Ball_delta(I) contains 0, so there exists J with J^2 = J and ||J - I|| <= O(eps). Then defines X . Y = R_J^{-1}(X) L_J^{-1}(Y) and verifies J is the unit.

- **Constructive in finite dim?:** yes-already. Full algorithm:
  1. Certify ||d_X(X^2 - X)|_{X=I} - 1|| <= c < 1 (use La-Ra-bounds; c = O(eps)).
  2. Run lem_invfun iteration: x_n = x_{n-1} + V^{-1}(0 - f(x_{n-1})) = x_{n-1} - (x_{n-1}^2 - x_{n-1}) starting from x_0 = I, targeting y = 0.
  3. Output J = fixed point (certified ||J - I|| <= O(eps), J^2 = J).
  4. Compute R_J^{-1} and L_J^{-1} as n x n matrix operations (since J is a near-identity idempotent, L_J and R_J are invertible near the identity).
  5. Define new multiplication table: (X . Y)_{ij} = [R_J^{-1}(X) L_J^{-1}(Y)]_{ij}.
  - Inputs: A (eps-Banach algebra, eps certified), I (unit element), tolerance.
  - Outputs: J (exact idempotent), new multiplication operator or multiplication table.

- **Depends on / Used by:** Depends on: lem_invfun, La-Ra-bounds (d_X2X). Used by: prop_polar (§5, assumes exact unit), lem_U_delta (§5), and effectively all of §5–11 since the exact-unit normalisation is applied at the start. The master th_main proof (L1414) implicitly relies on this.

- **Arbitrary precision needed?:** Yes. The iteration for J must converge to a point with J^2 = J certified to full precision. In floating point, J^2 - J will be O(machine_eps) but not exactly zero; arb gives a rigorous enclosure proving ||J^2 - J|| < bound.

- **Cross-check / oracle:** For an exact algebra (eps = 0): J = I (no change). For a constructed eps-algebra with known J_true: verify ||J - J_true|| <= O(eps). Test that X . J = X exactly (within precision) for several X.

- **Implementation notes / risks:**
  - The new multiplication X . Y = R_J^{-1}(X) L_J^{-1}(Y) changes the basis: it is equivalent to conjugating by the invertible element sqrt(J) (in exact C*-algebra case). In finite dim, store as a new n x n x n x n multiplication tensor or as a change-of-basis matrix.
  - R_J and L_J are invertible because J is a near-identity idempotent; their inverses are close to the identity. Compute via (I + small perturbation)^{-1} = I - small perturbation + O(eps^2).
  - Involution compatibility: verify (X . Y)^dag = Y^dag . X^dag by checking (R_J^{-1}(X) L_J^{-1}(Y))^dag = (L_J^{-1}(Y))^dag (R_J^{-1}(X))^dag. This requires J^dag = J, which holds if I^dag = I and J is in the self-adjoint component.
  - Risk: if the eps-algebra multiplication is stored as a table and the new multiplication requires R_J^{-1} and L_J^{-1}, the cost is O(n^3) per multiplication. For large n, precompute R_J^{-1} and L_J^{-1} as n x n matrices.

---

## Key definitions & notation introduced here

- Ball notation (L556–558): Ball^M_r(x_0) = {x in M : d(x_0, x) < r}, closed ball with bar. Used throughout §3–9.
- Frechet derivative notation (L560): d_x f(x) or df(x)/dx; linear map A -> B.
- L_X, R_X (L634–636): left/right multiplication operators on an eps-Banach algebra A; n^2 x n^2 matrices in finite dim.
- eq. d_X2X (L668–669): ||d_X(X^2 - X) - 1|| <= O(delta + eps) for X in Bar_delta(I). Tagged equation, referenced by prop_unit proof.
- |X| = (X^2)^{1/2}, sgn(X) = X(X^2)^{-1/2} (L514–516): functional-calculus definitions valid on domain ||X^2 - x_0^2 I|| < x_0^2.
- theta(X) = (I + sgn(X)) / 2 (L527–528): "spectral projection onto positive half" — maps elements with X^2 near I to exact idempotents.
- Iteration g_y(x) = x + V^{-1}(y - f(x)) (L582): the Banach contraction map for the inverse function theorem.
- eps-Banach algebra axioms (ax_prodnorm, ax_assoc): referenced at L642, L655 but defined in §2 (L≈200–260, not in this shard). Shard B uses them but does not define them.

---

## Dependencies OUT / INTO

### Depends on (from earlier sections):
- §2 definitions: eps-Banach algebra axioms ax_prodnorm, ax_assoc, exact unit definition. Involution axioms.
- §1.1 prop_hom_structure (L259/273): structural results on eps-C* algebras (background context, not directly invoked here).

### Used by (downstream sections):
- **§5 lem_U_delta (L699/707):** uses La-Ra-bounds directly.
- **§5 lem_gV (L758/765):** uses lem_invfun (the g_y iteration) and La-Ra-bounds.
- **§5 prop_polar (L809/813):** uses lem_invfun, fc-abs-sgn (|X|, sgn), La-Ra-bounds; requires exact unit (prop_unit applied first).
- **§6 lem_nontriv_projection (L931/934):** uses prop_P as the mechanism for extracting a projection from an approximate one (in the exact Banach case).
- **§6 prop_H-group (L971/1023):** uses the unitary group structure built on top of prop_unit, lem_invfun.
- **§7 lem_alpha (L1086/1104):** uses La-Ra-bounds and lem_invfun.
- **§8 prop_delta_hominc (L1194/1198), lem_approx (L1224/1256):** use lem_invfun, La-Ra-bounds.
- **§9 proof th_main (L1414):** indirectly depends on all of §3–4 via §5–8.
- **rem_X2 (L628):** specifically constrains §6–9 to work with approximate projections rather than exact ones.

---

## Open implementation questions / escalations

1. **Matrix sgn via DB vs. power series:** Denman-Beavers is faster (quadratic convergence) but requires no purely imaginary eigenvalues. Power series requires explicit x_0 and domain certification. Need a decision procedure: given X with certified enclosure, which algorithm to invoke? Recommend: check ||X^2 - I|| < 1 first; if yes, use DB. If x_0 is available and ||X^2 - x_0^2 I|| < x_0^2, use power series.

2. **Operator norm of L_X, R_X in finite dim:** Computing ||L_X|| = ||R_X|| as the spectral norm of an n^2 x n^2 matrix is O(n^6) naively (or O(n^3) via Kronecker-product SVD). For large n, use the bound ||L_X|| <= (1+eps)||X|| (L644) with certified ||X|| instead.

3. **Invertibility of L_J, R_J in prop_unit:** The proof asserts these are invertible but doesn't give an explicit inverse formula. In code: compute (L_J)^{-1} as the solution of J Z = W (i.e., Z = J^{-1} W in the algebra sense). Since J is a near-identity idempotent, L_J = 1 + (L_J - 1) with ||L_J - 1|| = O(eps), so (L_J)^{-1} = sum_{k>=0} (1 - L_J)^k converges for eps < 1. Certify convergence via arb.

4. **New multiplication table storage:** The new multiplication X . Y = R_J^{-1}(X) L_J^{-1}(Y) changes the algebra structure. In the eventual C code, represent this as a change-of-basis: define new basis e'_i = L_J^{-1}(e_i) and store the structure constants in the new basis. This avoids recomputing R_J^{-1} and L_J^{-1} repeatedly.

5. **Certified eps for eps-Banach algebras:** All bounds in §4 (La-Ra-bounds, d_X2X) require a certified value of eps (the approximate-associativity constant). In practice, eps is computed from the given multiplication table as max_{||X||=||Y||=||Z||=1} ||(XY)Z - X(YZ)|| / (||X||||Y||||Z||). This is a norm computation over the unit ball — in finite dim, it reduces to a max-singular-value problem on the associativity defect tensor. Implement as a certified optimisation in arb.

6. **Fragility threshold for holomorphic calculus:** The 3x3 example (spec-hfc) shows eigenvalue sensitivity of order t^{1/3} for perturbations of order t. If any future algorithm requires holomorphic calculus, a gap condition spec(X) separated from the contour by >= gap must be certified. This requires certified eigenvalue enclosures (arb `acb_mat_eig`), which are O(n^3) per matrix and expensive. Escalate to architectural decision if holomorphic calculus becomes necessary.

7. **Involution compatibility of prop_unit:** The proof of L686 asserts "all other requirements are also met" without explicit verification of involution compatibility. The claim (X . Y)^dag = Y^dag . X^dag follows from (R_J^{-1}(X))^dag = L_J^{-1}(X^dag) when J^dag = J. Verify this algebraically before implementation.
