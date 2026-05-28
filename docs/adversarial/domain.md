# Adversarial Domain Catalog — AIC Project

Operator-algebra / quantum-channel domain instances designed to punish the
AIC pipeline. General NLA torture matrices (non-normality, condition-number
sweeps, precision walls) are handled by a companion catalog; this file stays in
the domain lane: UCP maps, eps-C* algebras, projections, factorization.

Every instance has a tunable knob so the generator sweeps from mild to lethal.
"Expected correct behavior" means the certified arb path produces output
satisfying the paper's O(eps)/O(eta) bound (or fails loud); never silently
wrong.

Canonical line-number refs: `paper/src/approximate_algebras.tex` (2928 lines).
Short-hand: tex:NNN = line NNN; prop_P = tex:524; lem_invfun = tex:564;
cor_improvement = tex:1317; lem_nontriv_projection = tex:931;
th_almost_idemp = tex:2192; th_factorization = tex:2730.

---

## Family 1. UCP / Choi adversarial

### 1A. Extremal channel — near-minimal Kraus rank

**Description.** A UCP map whose Choi matrix C_Phi is PSD but has r-1 eigenvalues
near machine-zero for Kraus rank r << d^2. Constructed by taking a rank-1
channel (pure state, r=1) and blowing it up toward an extremal point of the
convex set of UCP maps by mixing in a tiny amount of the identity channel.
Evil knob: `rho_mix` in (0, eps_max); at 0 the map is pure/extremal.

Generator (acb_mat):
```
C = (1 - rho_mix) * |v><v| (rank-1 Choi matrix for a pure channel)
  + rho_mix * (I_{d^2} / d)    (maximally mixed channel Choi)
```
Normalize so that Tr_1(C) = I_H (unitality).

**Why evil.** The Choi eigendecomposition has r=1 eigenvalue near 1 and
d^2-1 eigenvalues near rho_mix/d. `lem_carrier` (tex:1724) locates the
carrier M as the column space of the partial trace of VV†; for rho_mix -> 0
the carrier collapses to a 1-dim subspace. Certified rank identification in
`lem_carrier` requires arb to certify that the small eigenvalues are genuinely
zero vs. nonzero. The idempotence defect eta = ||Phi^2 - Phi||_cb is non-zero
for rho_mix > 0 but small; the carrier rank toggles between 1 and d when
rho_mix crosses the arb-certification threshold.

**Which AIC ops attacked.** `ucp` (lem_carrier, cor_carrier, Choi->Kraus
conversion), `idemp_structure` (carrier-based decomposition), `assoc_ecsa`
(Step 1 of pipeline: Phi_tilde formation; near-zero eigenvalues in (1-4(Phi-
Phi^2))^{-1/2} Taylor series tex:2174), `factorize` (W_j extraction, lem_RC
lower bound on ||C_j||).

Predicted failure modes:
- `double` path: silently picks wrong carrier rank when rho_mix < 1e-12.
- `arb` path at low prec: Taylor series for Phi_tilde diverges slowly instead
  of refusing; the ball widens rather than failing.
- `lem_RC`: ||C_j|| lower bound (tex:2846) becomes 1 - O(eta) - O(rho_mix);
  if sum is negative the inversion of C_j^dag C_j blows up.

**Expected correct behavior.** For rho_mix >= eta (i.e., the mixing dominates
the idempotence defect), the carrier is correctly identified as full-dimensional
and the pipeline proceeds. For rho_mix << eta, the arb path should certify
rank-1 carrier and either succeed cleanly or fail-loud (insufficient gap),
never silently wrong.

**Why arb matters.** Certifying dim(carrier) = 1 vs. d requires an arb
eigenvalue gap certificate on Tr_F(VV†): the gap between the near-1 and
near-rho_mix/d eigenvalues is O(rho_mix/d) which can be below double precision
for d >= 8 and rho_mix < 1e-15.

Knob sweep: rho_mix in {1e-1, 1e-4, 1e-8, 1e-14, 0.0}.
Dimensions: d in {2, 4, 8, 16}.

---

### 1B. cb-norm vs operator-norm gap — ampliation-sensitive idempotence defect

**Description.** A UCP map Phi where ||Phi^2 - Phi|| (operator norm on B(H))
is tiny but ||Phi^2 - Phi||_cb is O(1/eps). Constructed by embedding a
near-idempotent error into off-diagonal blocks that only amplify under
1_{M_n} tensor Phi.

Generator:
```
Phi(X) = sum_{j=1}^{d} P_j Tr(gamma_j X)   (entanglement-breaking, Kraus rank d)
gamma_j chosen so that Phi^2 - Phi has operator norm O(eps) but
(1_{M_2} tensor (Phi^2 - Phi)) on off-diagonal 2x2 blocks achieves O(1).
```
Concretely: use the paper's own motivating example (tex:366-388) with eta =
eps_level:
```
P_0 = |0><0|,  P_1 = |1><1|
gamma_0 = ((1-eta, sqrt(eta(1-eta))), (sqrt(eta(1-eta)), eta))
gamma_1 = |1><1|
```
This gives ||Phi^2 - Phi||_cb = eta*sqrt(1-eta) (tex:378) and the difference
between cb-norm and operator norm is maximized at eta ~ 1/4.

**Why evil.** The only correct bound for eta-idempotence is the cb-norm
(tex:347, tex:354). Any routine that estimates ||Phi^2 - Phi|| in operator norm
and treats it as eta will underestimate the defect by up to a factor of n (the
ampliation level), missing the cb-norm by exactly the factor that the paper
warns about (tex:484). This directly attacks the `opspace` ampliation machinery.

**Which AIC ops attacked.** `ucp` (cb-norm computation via SDP or ampliation
sup; computing ||Phi^2 - Phi||_cb), `assoc_ecsa` (input hypothesis is cb-norm
eta; if naively estimated via operator norm the theorem's hypothesis is not
verified), `opspace` (1_{M_n} tensor construction, prop_inc_ext tex:1453).

Predicted failure modes:
- `ucp` double path: reports eta_op << eta_cb; downstream modules see
  "safe" eta and proceed; output violates the O(eta_cb) bound.
- `assoc_ecsa` with naive ||Phi^2 - Phi||: associativity defect at n=4
  ampliation is 4x larger than predicted; test fails.

**Expected correct behavior.** The `aic_ucp_idemp_defect_cb` function uses the
SDP (or an ampliation up to n=dim_H) and returns eta_cb; the pipeline is then
fed this value. The arb path certifies the ampliation loop over n = 1,...,d.

**Why arb matters.** The ampliation sup over n is truncated at some N; balls
from the arb path bound the truncation error. Without arb, the sup is computed
at a few n values and may systematically miss the maximum.

Knob sweep: eta in {0.20, 0.15, 0.10, 0.05, 1e-3}; d in {2, 4, 8}.

---

### 1C. Near-degenerate carrier — carrier rank drop

**Description.** A UCP map Phi where carrier dimension almost drops: the
subspace M has an eigenvalue of Tr_F(VV†) at distance `gap` from zero.
As gap -> 0, the carrier "nearly" loses a dimension.

Generator:
```
V: H -> H tensor F, built from
V = sqrt(1-gap) * (|0><0| tensor |f_1>) + sqrt(gap) * (|d-1><d-1| tensor |f_2>)
  + ... (other basis vectors at O(1) contribution)
```
The key: one singular value of Tr_F(VV†) equals `gap` and the rest are O(1).

**Why evil.** `lem_carrier` (tex:1724) must certify that M is exactly rank
(d-1) not rank d. In `idemp_structure`, the sub-block structure depends on
this rank exactly. `cor_carrier` (tex:1731) uses (X tensor 1_F)V = 0 as an
exact check; near-singular V means near-zero products that are not quite zero.

**Which AIC ops attacked.** `ucp` (lem_carrier, cor_carrier), `idemp_structure`
(block-diagonality of Im(Delta) relies on exact carrier; tex:1917),
`assoc_ecsa` (Phi_tilde image space A = Ker(1-Phi_tilde): if carrier is
borderline the eigenspace-1 of Phi_tilde fluctuates).

**Expected correct behavior.** For gap > eta (the idempotence defect), the
carrier dimension is certified correctly by arb. For gap < eta, the carrier is
not certifiable separately from the defect; fail loud with a carrier-gap
certificate failure, not a silent wrong-rank answer.

**Why arb matters.** gap = 1e-10 and d = 4 means the relevant singular value is
below double-precision noise; arb at prec = 200 can separate it from zero.

Knob sweep: gap in {0.5, 0.1, 1e-4, 1e-8, 1e-14}.

---

### 1D. Unital-but-barely

**Description.** Phi(I) = I + delta_u * E, where E is a traceless Hermitian
matrix with ||E|| = 1, and delta_u -> tolerance boundary.

Generator: take any CPTP map Phi_0 with Phi_0(I) = I exactly. Perturb:
Phi(X) = Phi_0(X) + delta_u * Tr(rho_0 X) * E for a fixed state rho_0 and
Hermitian E with Tr(E) = 0 and ||E|| = 1. Choose rho_0, E so the result is
still CP (possible for small delta_u).

**Why evil.** The unit axiom ax_eps_unit (tex:432) requires ||XI - X|| <= eps||X||
in A. If Phi(I) is only approximately I, the image A = Img(Phi_tilde) contains
an approximate unit with ||I_A - I_H|| = O(delta_u). The `prop_unit` fix
(tex:672) absorbs delta_u into eps; but if delta_u ~ eps_max the absorption
fails. The `assoc_ecsa` module depends on exact unitality of Phi_tilde (tex:2181:
Phi_tilde(1) = 1 exactly, because Phi_tilde = theta(2Phi-1) preserves I). The
perturbation breaks this.

**Which AIC ops attacked.** `ecstar` (unit axiom check), `unitfix` (prop_unit:
the exact-unit fix is an O(eps)-change; at delta_u ~ eps_max the change may
exceed the algebra's radius), `assoc_ecsa` (the proof of th_almost_idemp at
tex:2211 assumes exact unit).

**Expected correct behavior.** The unitality check in `ucp` detects delta_u
and propagates it as an additional eps-unit-violation term; the pipeline
either absorbs it correctly or refuses if delta_u > eps_max.

Knob sweep: delta_u in {0.0, 1e-6, 1e-3, 0.01, 0.10}.

---

## Family 2. eta-idempotent near the boundary

### 2A. eta -> 1/4^- — domain boundary of theta(2Phi-1)

**Description.** An eta-idempotent UCP map with eta = 1/4 - kappa for kappa
shrinking to 0.

The Taylor series for (1 - 4(Phi - Phi^2))^{-1/2} (tex:2174) converges iff
4eta < 1, i.e., eta < 1/4. At eta = 1/4 - kappa, the radius of convergence is
1 - 4kappa -> 0. The series requires O(1/kappa) terms to achieve fixed absolute
accuracy.

Generator: use the explicit example at tex:366-388 with eta = 0.24.
Alternatively: take a direct sum Phi = (1/2+eta/2) * Phi_idemp + (1/2-eta/2) *
(I - Phi_idemp) where Phi_idemp is exactly idempotent; tune to get ||Phi^2 -
Phi||_cb exactly equal to eta.

**Why evil.** The Phi_tilde computation (assoc_ecsa Step 1) uses the Taylor
series tex:2174. At kappa = 0.01, the series needs O(100) terms; each term
accumulates rounding error. The contraction solver lem_invfun (tex:564) has
c = 4eta -> 1; the iteration x_n = g_y(x_{n-1}) converges at rate (4eta)^n
which becomes arbitrarily slow. sgn(2Phi-1) is computed via Denman-Beavers or
Newton; the precondition ||X^2 - I|| < 1 (tex:520) becomes 4eta = 1 - 4kappa
which requires the precondition to be verified as strictly < 1 by arb.

**Which AIC ops attacked.** `funcalc` (sgn/abs/theta: precondition tex:516-520,
convergence of Taylor/Newton for ||X^2-I|| near 1), `contraction` (lem_invfun:
c = 4eta near 1, iteration rate -> 0), `assoc_ecsa` (Phi_tilde: series
truncation vs. convergence), `ecstar` (eps = O(eta) from tex:361 becomes O(1/4)
which exceeds eps_max for most downstream lemmas).

Predicted failure modes:
- `funcalc` double path: Newton for sgn converges but is unaware of near-unit
  ||X^2-I||; returns a "converged" value that is not actually sgn.
- `contraction` double path: declares convergence after 100 iterations with
  residual that looks small but is actually c^100 = (4eta)^100 away from the
  true fixed point.
- `assoc_ecsa`: Taylor series returns a value with no certified error ball;
  the arb path should detect that the ball has not converged.

**Expected correct behavior.** For kappa > 1e-6 at prec = 256: arb certifies
convergence and returns a tight ball. For kappa < eta_min(prec), the arb path
refuses with a "series not converged / precondition violated" error (fail loud),
never a plausible-but-wrong answer.

**Why arb matters.** The precondition ||X^2 - I|| < 1 is a strict inequality.
At kappa = 1e-8, double arithmetic will compute ||X^2 - I|| = 0.9999999904
(< 1 numerically) even if the true value exceeds 1; arb provides a certified
upper bound.

Knob sweep: kappa in {0.10, 0.01, 1e-3, 1e-4, 1e-8}. Pair with dim d in
{2, 4, 8} — larger d makes Taylor coefficients grow.

---

### 2B. Phi^2 - Phi concentrated on a degenerate subspace

**Description.** An eta-idempotent Phi where Phi^2 - Phi = eta * |v><v| tensor
rho_sub for some unit vector |v> and a density matrix rho_sub supported on a
near-degenerate subspace of B(H).

Generator:
```
Phi = Phi_idemp + eta * delta_channel
delta_channel(X) = |v><v| Tr(rho_sub X) - Phi_idemp(|v><v| Tr(rho_sub X))
```
(This is chosen to satisfy ||Phi^2 - Phi||_cb = eta while the defect is rank-1
in the superoperator sense.)

**Why evil.** When Phi^2 - Phi is concentrated on a near-degenerate subspace,
the Phi_tilde image A = Ker(1-Phi_tilde) is close to the image of Phi_idemp
but differs by a rank-1 correction. The associativity equations
(Phi_assoc1)/(Phi_assoc2) at tex:2198-2204 involve a cancellation between
two terms of size O(1) to produce a result of size O(eta); that cancellation
can catastrophically fail in double arithmetic when the subspace has a near-
zero gap. Specifically, the intermediate quantity W (tex:2385-2422) used in the
Phi_assoc1 proof involves (1-Pi) insertions where Pi = VV†; when Im(V) nearly
contains the defect subspace, (1-Pi) applied to the defect is near zero,
making the bound ||W|| = O(eta) hard to verify without cancellation-safe
arithmetic.

**Which AIC ops attacked.** `assoc_ecsa` (proof of th_almost_idemp relies on
tex:2296 bound ||R_0(X)|| = O(sqrt(eta)); this bound is tight when the defect
is concentrated), `ecstar` (axiom-defect estimators: the associativity defect
(XY)Z - X(YZ) may be atypically large for the specific X,Y,Z in the defect
subspace), `factorize` (the Phi_assoc equations power PhiDelta3 at tex:2813
which is needed for the Upsilon construction).

**Expected correct behavior.** The O(eta) bound holds universally; the arb
path certifies the bound holds even for this concentrated defect, by tracking
the cancellation with certified ball arithmetic.

Knob sweep: eta in {0.20, 0.05, 0.01, 1e-4}; gap_sub (eigenvalue gap of the
degenerate subspace) in {0.5, 0.1, 1e-3, 1e-8}.

---

## Family 3. eps-C* algebra adversarial

### 3A. Maximally concentrated associativity defect

**Description.** An eps-C* algebra A ⊆ M_N with associativity defect
h(X,Y,Z) = (XY)Z - X(YZ) that is NOT diffuse but concentrated: there exists
a single triple (X_0, Y_0, Z_0) with ||X_0|| = ||Y_0|| = ||Z_0|| = 1 and
||h(X_0, Y_0, Z_0)|| = eps (the maximum allowed by ax_assoc tex:412).

Generator: take B = M_2 and perturb the multiplication law:
```
mu_eps(X, Y) = XY + eps * f(X) tensor g(Y) * Z_0
```
where f, g are linear forms and Z_0 is a fixed matrix chosen to maximize
||h(X_0, Y_0, Z_0)||. The result satisfies the eps-C* axioms if eps is small
enough.

**Why evil.** The error-reduction step (cor_improvement, tex:1317) uses the
diagonal D = sum_j A_j tensor B_j of B (a genuine C* algebra) to construct a
new map w from an old delta-homomorphism v. The key calculation (tex:1282-1303)
uses the 2-cocycle equation for the multiplicativity defect g = G_v. If h is
maximally concentrated at one triple, the coboundary calculation
g(X,Y) = sum_j A_j h(B_j,X,Y) (tex:482) produces a w with the maximum
possible error constant c_0 (tex:1318), stressing the universal constant
claim. The cstar_build master loop must carry this constant faithfully through
all iteration steps without letting it grow with dim A.

**Which AIC ops attacked.** `errreduce` (cor_improvement: explicit diagonal
coboundary; the error constant c_0 must not depend on where h is concentrated),
`dhom` (lem_approx, tex:1224: Newton iteration's convergence rate depends on
the 2-cocycle magnitude), `cstar_build` (master loop: each cor_improvement
call uses the same c_0; if c_0 was inflated by the concentrated defect, errors
accumulate over the merging steps).

**Expected correct behavior.** The arb path certifies that the output
delta_0-inclusion satisfies delta_0 = c_0 * eps with a concrete c_0 that does
NOT depend on where in A the defect is concentrated, and does NOT depend on
dim A. This is the paper's key universality claim (tex:461, tex:484).

**Why arb matters.** The diagonal coboundary involves a sum over d^2 terms;
each term is O(delta/d^2) but summing to O(delta). Catastrophic cancellation
is possible in double arithmetic; arb balls track it cleanly.

Knob sweep: eps in {0.10, 0.05, 0.01, 1e-3}; dim A in {4, 9, 16, 25}
(dim = d^2 for M_d blocks). Test that c_0 is constant across the dim sweep.

---

### 3B. Near-violation of the C* axiom ||X†X|| >= (1-eps)||X||^2

**Description.** An eps-C* algebra where the C* axiom ax_C* (tex:428) is
tight: there exists X_0 with ||X_0†X_0|| = (1-eps)||X_0||^2 + kappa for
kappa -> 0.

Generator: take A = M_2 with slightly perturbed inner product:
```
||X||^2_A = ||X||^2_op + delta_C* * <X, X>_extra
```
where the extra term reduces the C* ratio for a specific direction X_0.
Equivalently: construct the multiplication directly so that X_0† star X_0
has spectral norm (1-eps + kappa)||X_0||^2.

**Why evil.** The approximate C* identity is used at tex:2235 in the proof of
th_almost_idemp (via PhiXdX tex:1692) and at tex:1204-1214 in prop_delta_hominc.
The lower bound propagates through lem_approx (tex:1224) as the key ingredient
for ||tilde_v||_lb (prop_delta_hominc: if a delta-hom v satisfies ||v(X)|| >=
eta||X|| for eta > 2*delta, then v is an (1-O(delta+eps))-isometry). Near-
violation of ax_C* makes this lower bound nearly vacuous.

**Which AIC ops attacked.** `ecstar` (axiom-defect estimators: must certify
ax_C* is satisfied within eps), `dhom` (prop_delta_hominc: isometry bound
collapses if ax_C* is near-violated), `cstar_build` (the isometry lower bound
at each extension step relies on ax_C*).

**Expected correct behavior.** For kappa > eps (i.e., safely within ax_C*
hypothesis), the pipeline proceeds. For kappa < 0 (genuine violation), the
ecstar module detects the violation and refuses (fail loud).

Knob sweep: kappa in {0.05, 0.01, 0, -0.001, -0.01}; eps in {0.01, 0.05}.

---

### 3C. Near-trivial projections — Route A spectral split failure

**Description.** An eps-C* algebra A where every Hermitian element H in A has
spectrum clustered near {0, 1}, with the spectral gap between the cluster near
t and 1-t bounded by gap_spec = eps + kappa for kappa -> 0.

Generator: take A = B(C^n) for small n, choose a basis of Hermitian elements
that all nearly commute with a fixed projection P_0 of rank k, so that each
H's eigenvalues are within eps/2 of the spectrum {0, ..., 0, 1, ..., 1}.
Perturb to reduce the gap below eps.

**Why evil.** Route A for lem_nontriv_projection (shard-D, the only constructive
route): Step 3 requires "locate t such that no eigenvalue lies in (t-g, t+g)
for some gap g > 0". If all eigenvalues of all Hermitian elements of A are
within eps of either 0 or 1, there is no certified gap g >> eps. The sgn-snap
(prop_P) requires gap g = Omega(1) via ||X^2 - I|| < 1 (tex:516); if g < eps,
prop_P's precondition fails. The arb path must either certify the gap or
report an escalation condition.

**Which AIC ops attacked.** `projection` (lem_nontriv_projection: Route A Step
3; arb gap certification), `funcalc` (prop_P: sgn/theta precondition check for
||X^2 - I|| < 1 tex:520), `cstar_build` (master loop start: the first
nontrivial projection is needed to begin the induction).

Predicted failure modes:
- `double` path: `zheev` returns eigenvalues with floating-point noise
  obscuring the gap; projection snaps to a near-trivial value.
- `arb` path: eigenvalue balls overlap; cannot certify nontriviality; correct
  behavior is to escalate per MODULE_PLAN escalation #1.

**Expected correct behavior.** When gap_spec > C * eps (C a safe constant),
arb certifies the gap and Route A succeeds. When gap_spec < eps, the module
raises the "no certified spectral gap" escalation condition rather than
returning a silently near-trivial P.

**Why arb matters.** The gap between eigenvalue balls is the core certificate;
arb's interval arithmetic is the only way to know the balls are disjoint.

Knob sweep: gap_spec / eps ratio in {10, 3, 1.5, 1.0, 0.5}; eps in {0.01, 0.05};
dim n in {4, 9, 16}.

---

### 3D. eps ~ c/n — dimension-dependent constant blowup regime

**Description.** An eps-C* algebra A with dim A = n and eps = c/n for a
constant c near 1. This is the regime where the naive Haar-diagonal construction
fails (tex:484) because the error grows proportional to n.

Generator: take A = ⊕_{j=1}^{k} M_{d_j} with d_1 = ... = d_k = d and n = k * d^2.
Choose the multiplication deformation so that each block has associativity
defect ~ eps = 1/(2n) per tex:484 note.

**Why evil.** The paper explicitly warns (tex:484) that the naive cohomological
approach (diagonal via Haar measure) fails when eps ~ c/n; only the incremental
construction with cor_improvement gives a dimension-independent constant c_0.
An instance with eps ~ c/n exposes whether the implementation's errreduce /
cstar_build is actually using the correct constructive diagonal (tex:1249-1254)
versus accidentally re-deriving something proportional to n.

**Which AIC ops attacked.** `errreduce` (cor_improvement: must use the explicit
B-diagonal, not Haar; the B-diagonal is the Pauli/Heisenberg-Weyl tensor
product tex:1249 which has ||D|| = 1 regardless of n), `cstar_build` (master
loop: if c_0 grows with dim at any stage, the loop's error budget blows up),
`dhom` (lem_approx: relies on the same diagonal).

**Expected correct behavior.** The delta_0-inclusion from cor_improvement
satisfies delta_0 = c_0 * eps with a constant c_0 that is the SAME for n = 4,
9, 16, 25. The test suite must assert this numerically: sweep n, measure c_0,
confirm it is constant (not growing). Failure mode: c_0 doubles when n doubles.

**Why arb matters.** The B-diagonal is a sum of d^2 terms each of size 1/d^2
(tex:1251); cancellations in computing ||D|| = 1 require ball arithmetic to
certify ||D|| is not > 1 + O(n*eps_machine).

Knob sweep: dim A in {4, 16, 64, 256}; eps = 0.5/dim_A (at the critical ratio);
also eps = 0.1/dim_A (safely below critical). Monitor c_0 as a function of dim.

---

## Family 4. Projection / corner adversarial

### 4A. delta-projection with delta -> 1/4 boundary (prop_P stress)

**Description.** A delta-projection P (Hermitian, ||P^2 - P|| = delta) with
delta = 1/4 - kappa for kappa -> 0.

prop_P (tex:524) requires delta < 1/4. The snap Ptilde = theta(2P-I) uses
sgn(2P-I) which requires ||(2P-I)^2 - I|| = 4||P^2 - P|| < 1 (tex:516). At
delta = 1/4 this bound is tight.

Generator: P = diag(1/2 + sqrt(1/4 - delta), 1/2 - sqrt(1/4 - delta), 0, ...,
0) in M_d. This is a Hermitian matrix with ||P^2 - P|| = delta exactly.

**Why evil.** The Taylor series for (1-4(P^2-P))^{-1/2} (tex:2174; the same
series used in Phi_tilde) diverges at delta = 1/4. The funcalc sgn routine
uses either Newton/Denman-Beavers or the direct Taylor; both require the
precondition to be strictly satisfied. Near delta = 1/4, the bound ||sgn(X) -
X|| = O(||X^2 - I||) (tex:520) becomes O(1), meaning the snap is useless for
certifying nontriviality.

**Which AIC ops attacked.** `funcalc` (prop_P: sgn/theta, tex:524-530;
precondition ||X^2-I|| < 1 is tight), `projection` (Route A: snapping
spectral projector to exact projection), `corner` (compression map Co_{P,Q}
= theta(L_P R_Q + R_Q L_P - 1) at tex:1057 uses the same theta function;
the operand L_P R_Q + R_Q L_P - 1 is nearly +-I when P is near-trivial).

Predicted failure modes:
- `funcalc` double path: Newton for sgn converges to a fixed point of X^2=I
  but the fixed point is not near the "correct" branch; silently returns wrong
  sgn.
- `arb` path: correct behavior is fail-loud when kappa < kappa_min(prec),
  not a ball that brackets 0 and 1 simultaneously.

**Expected correct behavior.** arb certifies for each kappa whether
4*delta < 1 is provably true (tight ball); if not, refuses rather than proceeding.
The output Ptilde satisfies ||Ptilde - P|| <= ||2P-I|| * O(delta) (tex:530)
with a certified ball on this bound.

Knob sweep: kappa in {0.10, 0.01, 1e-3, 1e-5, 0.0 (exact boundary)};
dim d in {2, 4, 8}.

---

### 4B. Two projections with near-rank-deficient corner S_{P,Q}

**Description.** Two nonvanishing delta-projections P, Q in an eps-C* algebra
where the corner subspace S_{P,Q} has dim 1 but with the one non-zero singular
value of Co_{P,Q} near zero (i.e., the corner exists but is nearly empty).

Generator in M_n: P = diag(1, 0, ..., 0) + eps_perturb * H_P,
Q = diag(0, 1, 0, ..., 0) + eps_perturb * H_Q for small perturbations H_P, H_Q.
The exact S_{P,Q} has dimension 1 (spanned by |0><1|) but the compression
map Co_{P,Q} = theta(L_P R_Q + R_Q L_P - 1) has an eigenvalue near 0 because
P and Q are "almost orthogonal" in the (1-P) subspace.

**Why evil.** lem_1d_proj (tex:1179): "If P and Q are both one-dimensional
delta-projections, then dim S_{P,Q} <= 1." The near-rank-deficiency makes the
inner product on S_{P,Q} (lem_PQ_Hilb, tex:1123) nearly degenerate: the inner
product <X,X> from tex:1126 satisfies |<X,X> - ||X||^2| <= O(delta+eps)||X||^2
which means <X,X> could be near zero while ||X|| is O(1). The compressed
product norm bound (lem_PQR, tex:1162) uses ||X.Y|| ≈ ||X||||Y||; near-rank-
deficiency makes this bound saturate from below rather than above.

**Which AIC ops attacked.** `corner` (Co_{P,Q} compression map; lem_1d_proj;
lem_PQ_Hilb inner product degeneration; lem_PQR product norm bound), `cstar_build`
(Stage 2 uses equivalence classes of 1-dim projections via lem_1d_proj;
near-rank-deficiency of S_{P,Q} makes the equivalence decision borderline).

**Expected correct behavior.** For dim S_{P,Q} = 1 with a certified non-zero
singular value (bounded away from 0 by some gap > O(delta+eps)), the inner
product is certified positive definite by arb. If the singular value is
below O(delta+eps), the equivalence decision escalates (fail loud).

Knob sweep: gap (singular value of Co_{P,Q}) in {0.5, 0.1, 0.01, delta+eps, 0.1*(delta+eps)};
delta in {1e-3, 0.05, 0.1}.

---

### 4C. Clustered spectrum of 2P-I — eig-free sgn must carry it

**Description.** A delta-projection P where the spectrum of 2P-I has clusters
near +1 and -1, with eigenvalues at distances gap_+ from +1 and gap_- from -1,
both small.

Generator: P = (I + H)/2 where H = diag(1-2*gap_+, 1-2*gap_+, -1+2*gap_-, ...,
-1+2*gap_-) in M_d, so ||P^2 - P|| = delta = O(gap_+, gap_-).

**Why evil.** For sgn(H) where H = 2P-I, the eigenvalues near +-1 require the
algorithm to resolve which "side" of 0 they are on. The holomorphic contour
approach (tex:546-550) places the contour between the clusters; near-clustering
means the contour must thread a narrow gap. The Newton / Denman-Beavers
iterations converge faster when eigenvalues are well-separated from 0; near-
unit eigenvalues slow convergence dramatically. The eig-free sgn route
(Denman-Beavers or Zolotarev rational approx) must handle eigenvalues near +-1
without a contour; this is exactly the case where arb's eigenvalue-ball
arithmetic is the correct certificate.

**Which AIC ops attacked.** `funcalc` (sgn near-unit spectrum; the Newton
convergence rate degrades as eigenvalues approach +-1; the Taylor series
(I + (X^2-I))^{-1/2} around X=I diverges for eigenvalues near -1), `projection`
(Route A: certifying the spectral gap between eigenvalue clusters near +1 and -1
of H = 2P-I so that prop_P's output Ptilde is nontrivial).

Predicted failure modes:
- `funcalc` double path: Newton for sgn "converges" to a matrix that is close
  to one of the +-I trivial fixed points rather than the correct sgn(H);
  non-normal perturbations make this hard to detect without arb.
- `arb` path: eigenvalue balls for near-unit eigenvalues are wide; may not
  certify which side of 0 they are on; correct behavior is fail loud.

**Expected correct behavior.** For gap_+ = gap_- = g: if g > C*delta (for a
concrete C derived from prop_P), arb certifies the sign assignment and returns
Ptilde with certified ||Ptilde - P|| <= O(delta). For g < delta, escalate.

Knob sweep: gap_+ = gap_- = g in {0.5, 0.1, 0.01, delta, 0.1*delta};
delta in {1e-3, 0.05}; dim d in {4, 8, 16}.

---

## Family 5. delta-homomorphism / main-loop adversarial

### 5A. delta-homomorphism near cor_improvement threshold

**Description.** A delta-homomorphism v: B -> A (B a finite-dim C* algebra,
A an eps-C* algebra) with delta = delta_max - kappa for kappa -> 0, where
delta_max is the threshold in cor_improvement (tex:1318).

Generator: take B = M_2, A = M_2 with a small eps-perturbation.
Construct v(X) = X + delta * R(X) where R is a random linear map with
||R|| = 1. Tune delta so that ||v(XY) - v(X)v(Y)|| <= delta ||X||||Y|| exactly,
and v is a delta-inclusion.

**Why evil.** cor_improvement (tex:1317) says: if there exists a delta_max-
inclusion of B into A, there is also a c_0*eps-inclusion. The explicit constant
c_0 and delta_max are not stated in the theorem; they must be extracted from
the proof of lem_approx (tex:1224-1314). The key convergence condition for the
Newton iteration (tex:1311-1313) is that the multiplicativity defect satisfies
||G_v||_op <= delta and the 2-cocycle equation reduces G_{v+w} to O(delta^2 +
eps). Near delta = delta_max, the first Newton step may not decrease the defect.

**Which AIC ops attacked.** `dhom` (lem_approx Newton iteration; convergence
when delta is near delta_max), `errreduce` (cor_improvement: the threshold
condition; extracting c_0 and delta_max explicitly), `cstar_build` (master loop:
error-reduction is called after each extension step; if delta is near threshold
at any step, the loop may not converge in the claimed number of iterations).

**Expected correct behavior.** For delta = delta_max - kappa: if kappa > 0 and
eps << kappa, the Newton iteration in lem_approx converges and produces a
c_0*eps-homomorphism. If delta >= delta_max, the module refuses (fails loud:
"delta exceeds threshold delta_max"). The value delta_max must be certified by
arb as part of the module's precondition check.

**Why arb matters.** The 2-cocycle residual after one Newton step is O(delta^2
+ eps); near delta_max, delta^2 ~ delta_max^2 and eps ~ delta_max are
comparable; the ball arithmetic tracks whether the residual is genuinely < delta.

Knob sweep: delta / delta_max in {0.5, 0.8, 0.95, 1.0, 1.05 (over threshold)};
eps in {0.001, 0.01}; dim B in {1, 4, 9}.

---

### 5B. High-multiplicity block structure — degenerate eigenvalue merging

**Description.** An eps-C* algebra A ≅ B(C^d)^{oplus k} (k identical copies of
M_d) with d = 2, k large (k in {8, 16, 32}). The incremental construction in
th_main must merge k isomorphic blocks via cor_merge_sum (tex:1352).

Generator: build A as a direct sum of k copies of M_d with a shared small eps
perturbation to the multiplication (so they are not exactly identical):
```
A = ⊕_{j=1}^{k} M_d + O(eps) cross-coupling terms
```
The cross-coupling terms make the equivalence class partition (dim S_{P_j,P_l}
= 0 or 1) computationally ambiguous.

**Why evil.** The cstar_build master loop (proof of th_main at tex:1414)
groups 1-dim projections into equivalence classes by checking dim S_{P_j, P_l}.
With k identical blocks and O(eps) cross-coupling, the off-diagonal S_{P_j,P_l}
have singular values exactly at O(eps). The merging step (lem_merging tex:1325)
assumes the maps gamma_{jk}: S_{Pi_j, Pi_k} -> S_{P_j, P_k} satisfy the
conditions tex:1328-1334 with a delta that is O(eps) after cor_improvement;
but merging k copies compounds the error by k, potentially overrunning delta_max.
The explicit bound in lem_merging at tex:1325-1351 uses pq(delta+eps) < c_0;
for p = q = k and large k this fails unless delta = O(eps/k^2).

**Which AIC ops attacked.** `corner` (S_{P_j,P_l} dimension certification;
lem_1d_proj for each pair j,l), `cstar_build` (equivalence class partition;
cor_merge_sum applied k times; error accumulation over k iterations of the
master loop), `errreduce` (cor_improvement must reduce delta to c_0*eps before
each merging step; with k large blocks, k successive applications of
cor_improvement are needed).

**Expected correct behavior.** The loop terminates in O(k^2) iterations with a
final isomorphism delta_0-close to the exact B = M_d^{oplus k}, with delta_0 =
c_0*eps independent of k (the universality claim tex:461). The test must assert
delta_0 is constant as k grows from 2 to 32.

Knob sweep: k in {2, 4, 8, 16, 32}; d in {2, 3}; eps in {0.01, 0.05};
cross-coupling strength in {0, 0.001*eps, 0.01*eps, eps}.

---

## Family 6. Factorization end-to-end adversarial

### 6A. Worst-case composite O(eta) bound across the full pipeline

**Description.** An eta-idempotent Phi engineered so that EVERY sub-bound
tex:2733 (DelUps), tex:2739 (UpsDel2) achieves its O(eta) bound simultaneously.
This maximizes the composite constant C in ||DeltaUpsilon - Phi||_cb <= C*eta.

Generator: compose the adversarial instances from families 2A, 3A, and 5A:
```
Phi: eta-idempotent with eta = 0.20 (near the 1/4 boundary)
A = assoc_ecsa(Phi): associativity defect concentrated at a single triple (Family 3A)
v: B -> A: delta-homomorphism near delta_max (Family 5A)
```
Each sub-step contributes its maximum O(eta) constant; the composite is their
product.

**Why evil.** The paper states th_factorization with O(eta) but does not give
the composite constant explicitly (shard H escalation #2). For the factorize
module to be trustworthy, ||DeltaUpsilon - Phi||_cb <= C*eta must be verified
with a concrete C. This instance maximizes C and forces the arb path to track
the product of all intermediate bounds through tex:2742-2769.

**Which AIC ops attacked.** `factorize` (end-to-end: all five steps of the
pipeline; the unitary 1-design step at tex:2786-2800 is sensitive to Delta_tilde
norm; the lem_RC lower bound at tex:2840 requires ||C_j|| >= 1 - O(eta) across
all blocks j), `assoc_ecsa` + `cstar_build` + `errreduce` (all feedin to
factorize), `ucp` (final cb-norm certification of ||DeltaUpsilon - Phi||_cb).

**Expected correct behavior.** The arb path certifies ||DeltaUpsilon - Phi||_cb
<= C*eta and ||UpsilonDelta - 1_B||_cb <= C*eta with a concrete C. The
double path should agree within a small factor. The test must confirm:
(a) Enc = Upsilon*, Dec = Delta* are valid quantum channels (CPTP maps, certified);
(b) Dec Enc = (UpsilonDelta)* satisfies ||Dec Enc - 1_{B*}||_cb <= C*eta;
(c) Enc Dec = (DeltaUpsilon)* satisfies ||Enc Dec - Phi*||_cb <= C*eta.

**Why arb matters.** The pipeline composes O(1) steps each introducing O(eta)
error; the arb path provides a single certified ball on the final error that
accounts for all intermediate rounding.

Knob sweep: eta in {0.20, 0.10, 0.05, 0.01, 1e-4}; dim H in {4, 9, 16, 25}.

---

### 6B. lem_RC near-failure — ||C_j|| near zero

**Description.** A Phi such that in the factorization pipeline, the operator
C_j from lem_RC (tex:2843-2846) satisfies ||C_j|| = 1 - O(eta) with the
constant tight: ||C_j|| = 1 - A*eta for A approaching 1/O(eta_threshold).

lem_RC proves 1 - O(eta) <= ||C_j|| <= 1; the lower bound relies on
||Delta(1_{L_j})||_cb >= 1 - O(eta) (Delta_norm, tex:2806) and PhiDelta2
(tex:2810). Near-failure of Delta_norm means ||Delta_tilde|| is barely within
(1-O(eta), 1+O(eta)).

Generator: choose Phi such that the isomorphism v: B -> A from th_main_ext has
v^{-1} with large cb-norm near the threshold (1+delta_max).

**Why evil.** The Upsilon construction at tex:2859-2897 requires choosing
xi_j with ||C_j xi_j|| near 1. If ||C_j|| = 1 - O(eta) is achieved with a
tight bound, the choice of xi_j = leading singular vector of C_j gives
||C_j xi_j|| = sigma_max(C_j) = ||C_j|| = 1 - O(eta). The final step
(tex:2888-2892) approximates Y_j tensor (xi_j† C_j† C_j xi_j) ≈ Y_j;
the error is (1 - ||C_j xi_j||^2) which is O(eta). If tight, Upsilon's
normalization step tex:2897 inverts (Upsilon'(I_H))^{-1/2} = (C_j† C_j xi_j
term)^{-1/2} with an operator near zero, causing precision loss.

**Which AIC ops attacked.** `factorize` (lem_RC bounds; Upsilon construction;
normalization inversion of (Upsilon'(I_H))^{-1/2}), `ucp` (final Upsilon is
UCP; normalization near singular makes the UCP check marginal).

**Expected correct behavior.** For ||C_j|| >= 1 - C0*eta with C0 an explicit
constant: arb certifies invertibility of Upsilon'(I_H) (certified spectral gap
> 0 for the inverse), and the Upsilon construction succeeds. For ||C_j|| below
the invertibility threshold, fail loud.

Knob sweep: A (constant in ||C_j|| = 1 - A*eta) in {0.1, 0.5, 1.0, 2.0, 5.0};
eta in {0.05, 0.10, 0.15}.

---

## eta=0 exact-idempotent oracle

The eta=0 / exact-idempotent oracle is the zero-defect anchor from which all
adversarial sweeps perturb away.

When Phi is exactly idempotent (||Phi^2 - Phi||_cb = 0):

1. Phi_tilde = Phi exactly (theta(2Phi-1) = Phi since 2Phi-1 has eigenvalues
   exactly +-1 as a cb-algebra element).
2. A = Im(Phi) is a genuine C* algebra (Choi-Effros theorem, tex:344).
3. All eps-C* axiom defects are 0: associativity defect = 0, C* identity is
   exact, unit is exact.
4. th_main applied to A (eps=0) must return an isomorphism v with delta_0 = 0:
   v is an exact *-isomorphism B -> A.
5. th_factorization gives Delta, Upsilon with ||DeltaUpsilon - Phi||_cb = 0
   and ||UpsilonDelta - 1_B||_cb = 0 exactly — the Choi-Effros / th_idemp_structure
   decomposition (tex:318).

Every module in the AIC pipeline must reduce to the exact case at eta=0.
Specifically:
- `funcalc` sgn/theta: exact on elements with exact +-1 spectrum.
- `projection`: exact spectral projector; Ptilde = P exactly.
- `corner`: Co_{P,Q} is an exact idempotent; S_{P,Q} is an exact corner.
- `cstar_build`: produces B = A (no error, no error-reduction needed).
- `factorize`: Delta = exact inclusion, Upsilon = exact left inverse.
- All defect estimators return 0.

The oracle test: generate a random exact idempotent Phi (e.g., Phi(X) = V†XV
for a random isometry V with dim H = 4, dim carrier = 2), run the full pipeline,
assert all defect norms are < 1e-50 in double and that the arb path certifies
them as zero (ball contains 0 and has radius < 1e-100 at prec=256).

This is the PRIMARY regression check: any perturbation of the oracle that
causes a non-zero output (beyond O(eta)) is a bug.

---

## Sweep Design

### Knobs and ranges

| Knob | Symbol | Range | Purpose |
|------|--------|-------|---------|
| Idempotence defect | eta | 0 (oracle), 1e-6, 1e-3, 0.05, 0.10, 0.20, 1/4-1e-8 | Main domain parameter |
| eps-C* defect | eps | 0, 1e-4, 1e-2, 0.05, 0.10 | eps-algebra axiom quality |
| delta-projection defect | delta | 0, 1e-4, 0.05, 0.10, 1/4-1e-6 | prop_P domain boundary |
| delta-homomorphism defect | delta_hom | 0, 0.5*delta_max, 0.9*delta_max, delta_max | cor_improvement threshold |
| Block multiplicity | k | 1, 2, 4, 8, 16, 32 | Merging stress; universality |
| Block dimension | d | 2, 3, 4, 8 | n = k*d^2; combined with k |
| Spectral gap | gap | 0 (degenerate), 1e-8, 1e-4, 0.01, 0.1, 0.5 | Route A certificability |
| Carrier rank drop gap | gap_c | 0 (rank drop), 1e-10, 1e-4, 0.1 | lem_carrier stress |
| 1/4-boundary proximity | kappa | 0 (exact boundary), 1e-8, 1e-4, 0.01, 0.10 | Domain edge for theta |
| Defect concentration | conc | 0 (diffuse), 0.5 (mixed), 1.0 (maximally concentrated) | Coherency of error structure |
| Precision | prec | 53 (double baseline), 128, 256, 512 | arb certification depth |

### Instance sweep matrix

Critical combinations to sweep systematically:
1. (eta, kappa): eta in {0, 0.05, 0.15, 0.24}, kappa in {0.01, 0.001, 1e-6} — tests the
   1/4 boundary across all downstream modules.
2. (k, eps): k in {1, 4, 16, 32}, eps = 0.5/k^2 (at dimension-dependent blowup) — tests
   universality of c_0.
3. (delta, gap_spec): delta in {0, 0.05, 0.1}, gap_spec / eps in {0.5, 1.0, 2.0, 10.0} —
   tests projection Route A certificability.
4. (eta, conc): eta in {0.01, 0.10}, conc in {0, 0.5, 1.0} — tests associativity defect
   concentration through factorize.
5. (Kraus_rank, gap_c): Kraus_rank in {1, 2, d}, gap_c in {0, 1e-10, eta} — carrier
   degeneration across ucp and idemp_structure.

---

## Prioritized Shortlist: 8–12 Single Most Punishing Domain Instances

Ranked by predicted lethality (most likely to expose a real bug):

1. **2A at kappa = 1e-6** (eta = 1/4 - 1e-6, d = 8): Direct attack on the
   domain boundary of the whole pipeline; funcalc, contraction, and assoc_ecsa
   must all handle near-divergent series. The arb-vs-double disagreement here
   is catastrophic in double, detectable only by arb.

2. **3D at dim A = 64, eps = 0.5/64** (dimension-dependent blowup regime):
   Directly tests the universality claim (tex:461 c_0 independent of dim).
   If errreduce uses the wrong diagonal, c_0 doubles per block. Sweeping dim
   from 4 to 256 confirms or refutes universality quantitatively.

3. **1B at eta = 0.20, d = 8** (cb-norm vs operator norm gap): Any routine
   that calls ||Phi^2 - Phi|| instead of ||Phi^2 - Phi||_cb silently
   underestimates eta by a factor of d. Exposing this in the assoc_ecsa and
   opspace modules is the highest-ROI test for cb-norm correctness.

4. **3C at gap_spec / eps = 1.0, d = 9** (no-certified-gap case):
   Tests the exact boundary where Route A for lem_nontriv_projection
   should escalate vs. proceed. Any off-by-one in the gap certificate
   produces a near-trivial projection that silently corrupts the
   cstar_build master loop.

5. **5A at delta / delta_max = 0.95, eps = 0.01, k = 4** (near-threshold
   delta-hom): Tests whether cor_improvement correctly identifies its own
   threshold and whether the Newton iteration in lem_approx converges near
   the threshold. The most direct stress on dhom and errreduce.

6. **6A at eta = 0.20, dim H = 16** (worst-case composite bound): Forces
   every sub-step of th_factorization to its O(eta) bound simultaneously.
   Reveals whether the composite C in ||DeltaUpsilon - Phi||_cb <= C*eta is
   genuinely constant or hides dimension-dependent constants.

7. **5B at k = 32, d = 2, eps = 0.05** (high-multiplicity merging): The
   hardest test for cstar_build's master loop; k^2 = 1024 equivalence-class
   pairs to check, each near the O(eps) threshold for S_{P_j, P_l}.

8. **4A at kappa = 1e-5, d = 8** (prop_P domain boundary): The funcalc
   module's most basic correctness test at the hardest single-module boundary.
   Denman-Beavers / Newton for sgn must refuse for kappa below prec-dependent
   threshold, not silently converge.

9. **1C at gap_c = 1e-8, d = 4** (near carrier rank drop): Most likely to
   reveal a silent wrong-rank answer in the carrier computation, which then
   corrupts everything downstream in idemp_structure and assoc_ecsa.

10. **2B at eta = 0.05, gap_sub = 1e-4** (concentrated defect): The most
    demanding test for the O(sqrt(eta)) cancellation bound in the Phi_assoc1
    proof (tex:2296); double arithmetic loses all significant digits in this
    cancellation; only arb tracks it correctly.

11. **4B at gap = delta + eps, d = 6** (marginal 1-dim corner): Tests lem_1d_proj
    and the equivalence-class partition in cstar_build at the exact equivalence
    boundary; a one-off error in the rank check corrupts the final B structure.

12. **eta=0 oracle at d = 8** (zero-defect anchor): Must pass perfectly before
    any other instance is meaningful; it is the cheapest oracle and the one
    most likely to catch basic implementation errors in the pipeline.

---

## Most at-risk AIC module

Based on the above analysis, `factorize` is the most at risk. It is the final
composition of every other module, has an outline (not a closed proof) at
tex:2742-2769 (shard H escalation #3), depends on an unextracted composite
constant (shard H escalation #2), and requires at least four separate certified
inversions (Delta'(I_B)^{-1/2}, Upsilon'(I_H)^{-1/2}, v^{-1} from th_main_ext,
C_j^{-1/2} from lem_RC). Any one of these near-singular inversions can silently
return a plausible-but-wrong result in double arithmetic while only the arb
path certifies correctness.

The second-most-at-risk module is `cstar_build`, because the universality of
c_0 (the constant in cor_improvement) is the paper's main claim and the most
likely source of a dimension-dependent bug in an implementation that
accidentally uses the Haar-measure diagonal instead of the explicit B-diagonal.
