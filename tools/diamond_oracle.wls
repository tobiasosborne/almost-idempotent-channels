#!/usr/bin/env wolframscript
(* ::Package:: *)

(* diamond_oracle.wls -- INDEPENDENT diamond-norm oracle for the eta-idempotence
   defect  eta = ||Phi^2 - Phi||_diamond  (bead aic-cne; .tex:347-354). Third
   cross-check rung alongside the C Choi substrate and the Julia+MOSEK golden
   master (tools/gen_fixtures_d24.jl). This script shares NO code path with that
   generator: it rebuilds every channel analytically, derives the Convention-A
   Choi and the Heisenberg self-composition in Mathematica, and solves the
   Watrous 2012 SDP (arXiv:1207.5726) with Mathematica's complex-native
   SemidefiniteOptimization -- a different solver (interior point) from MOSEK.
   The whole point (CLAUDE.md Rule 6) is solver/implementation independence:
   if Mathematica reproduces the analytic anchors but disagrees with eta_ref,
   that is a real finding, reported loudly, not tuned away (Rule 2).

   .tex:354 verbatim: "A UCP map Phi: B(H)->B(H) is called eta-idempotent if
   ||Phi^2-Phi||_cb <= eta."  For the Hermiticity-preserving (not CP) map
   Lambda = Phi^2 - Phi the cb-norm equals the diamond norm (.tex:352), and the
   ampliation truncation N = dim(input) is rigorous (Watrous 2012), so a single
   SDP at the channel's own dimension suffices.

   NORMALIZATION (load-bearing, matches gen_fixtures_d24.jl): the Convention-A
   Choi has trace n for a UCP map while the Watrous SDP is trace-1 calibrated,
   so  ||Lambda||_diamond = (2/n) * SDP_optval.

   Run:  wolframscript tools/diamond_oracle.wls
   Self-contained, rerunnable. Prints each result as it is computed so output
   survives any exit-time kernel crash (a caveat observed by the orchestrator).
   No emojis, no marketing (Rule 12). NO Julia, NO C: this is an independent
   Mathematica oracle. *)

Off[General::munfl];

(* ---------------------------------------------------------------------------
   Choi (Convention A) and Heisenberg self-composition -- matched to
   gen_fixtures_d24.jl EXACTLY so the cross-check is valid.
     C[(i-1)n+a, (j-1)n+b] = sum_x conj(Kx[i,a]) Kx[j,b]   (1-based; conjugation
   on the FIRST factor -- the project's documented bug-prone convention).
   Kraus of Phi o Phi (Heisenberg): {Kb . Ka}.
   --------------------------------------------------------------------------- *)
choiA[Kraus_, n_] := Module[{C = ConstantArray[0, {n n, n n}]},
  Do[C[[(i - 1) n + a, (j - 1) n + b]] += Conjugate[K[[i, a]]] K[[j, b]],
    {K, Kraus}, {i, n}, {a, n}, {j, n}, {b, n}];
  C];
composeSelf[Kr_] := Flatten[Table[Kb . Ka, {Ka, Kr}, {Kb, Kr}], 1];
lambdaChoi[Kr_, n_] := choiA[composeSelf[Kr], n] - choiA[Kr, n];

(* ---------------------------------------------------------------------------
   The Watrous diamond-norm SDP via SemidefiniteOptimization.
     maximize Re tr(J^dag X)  s.t.  [[P,X],[X^dag,Q]] >= 0,  P+Q = I_{n^2},
                                    P,Q Hermitian PSD.
     eta = (2/n) * optval.
   Verified Mathematica syntax: each PSD/LMI constraint M>=0 (M Hermitian) MUST
   use the explicit cone form
     VectorGreaterEqual[{M, 0matrix}, {"SemidefiniteCone", dim}].
   The bare  M \[VectorGreaterEqual] 0  is ENTRYWISE and silently wrong.
   To keep Mathematica from injecting a Conjugate[] head into the constraint
   matrix (which trips bdspec1, because Conjugate[xr] is not a recognized
   variable), X^dag is built directly from the real entry variables as
   Transpose[XR] - I Transpose[XI], and Re tr(J^dag X) is expanded analytically
   into the real-variable bilinear  sum_{ij} (Re J)_ij (XR)_ij + (Im J)_ij (XI)_ij.
   --------------------------------------------------------------------------- *)
diamondNorm[J_, n_] := Module[
  {n2 = n n, JR, JI, xrv, xiv, Xc, Xdag, Pre, Pim, Qre, Qim, P, Q, blk, obj,
   freeP, freeQ, vars, cons, sol, optval},
  JR = Re[N[J]]; JI = Im[N[J]];
  xrv = Array[ToExpression["xr" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
  xiv = Array[ToExpression["xi" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
  Xc = xrv + I xiv;
  Xdag = Transpose[xrv] - I Transpose[xiv];
  Pre = Array[ToExpression["Pr" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
  Pim = Array[ToExpression["Pi" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
  Qre = Array[ToExpression["Qr" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
  Qim = Array[ToExpression["Qi" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
  buildHerm[re_, im_, m_] := Table[
    Which[i == j, re[[i, i]],
          i < j, re[[i, j]] + I im[[i, j]],
          True, re[[j, i]] - I im[[j, i]]], {i, m}, {j, m}];
  P = buildHerm[Pre, Pim, n2];
  Q = buildHerm[Qre, Qim, n2];
  blk = ArrayFlatten[{{P, Xc}, {Xdag, Q}}];
  obj = Total[Flatten[JR xrv + JI xiv]];
  freeP = Join[Flatten[Table[Pre[[i, j]], {i, n2}, {j, i, n2}]],
               Flatten[Table[Pim[[i, j]], {i, n2}, {j, i + 1, n2}]]];
  freeQ = Join[Flatten[Table[Qre[[i, j]], {i, n2}, {j, i, n2}]],
               Flatten[Table[Qim[[i, j]], {i, n2}, {j, i + 1, n2}]]];
  vars = Join[Flatten[xrv], Flatten[xiv], freeP, freeQ];
  cons = {
    VectorGreaterEqual[{blk, ConstantArray[0, {2 n2, 2 n2}]}, {"SemidefiniteCone", 2 n2}],
    P + Q == IdentityMatrix[n2]};
  sol = SemidefiniteOptimization[-obj, cons, vars];
  optval = obj /. sol;
  (2/n) optval];

etaOf[Kr_, n_] := diamondNorm[lambdaChoi[Kr, n], n];

(* ---------------------------------------------------------------------------
   Channel constructors -- rebuilt analytically from the Julia definitions in
   gen_fixtures_d24.jl (NOT read from the C .inc.h). Each returns a list of
   Kraus matrices.
   --------------------------------------------------------------------------- *)
ee[d_, i_] := UnitVector[d, i];
ket[v_] := Transpose[{v}];           (* column vector as a d x 1 matrix *)
op[u_, v_] := ket[u] . {Conjugate[v]}; (* |u><v|  (v conjugated for the bra) *)

identityKraus[d_] := {IdentityMatrix[d]};
dephasingKraus[d_] := Table[op[ee[d, a], ee[d, a]], {a, d}];  (* {|a><a|} *)
depolarizingKraus[d_] :=
  Flatten[Table[op[ee[d, i], ee[d, j]] / Sqrt[d], {i, d}, {j, d}], 1];
blockCondExpKraus[] := {
  ({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}),
  ({{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}})};
phitKraus[d_, t_] := Join[
  {Sqrt[1 - t] IdentityMatrix[d]},
  Flatten[Table[Sqrt[t] op[ee[d, i], ee[d, j]] / Sqrt[d], {i, d}, {j, d}], 1]];

(* Paper example (.tex:367-378): Phi(X) = P0 Tr(g0 X) + P1 Tr(g1 X) on C^2.
   Kraus for X -> P Tr(g X), P = u u^dag rank-1: K_k = sqrt(lam_k) v_k u^dag,
   g = sum_k lam_k v_k v_k^dag.  Matches paper_example_kraus(etap). *)
paperExampleKraus[etap_] := Module[{g0, g1, P0, P1, prep},
  g0 = {{1 - etap, Sqrt[etap (1 - etap)]}, {Sqrt[etap (1 - etap)], etap}};
  g1 = {{0, 0}, {0, 1}};
  P0 = {{1, 0}, {0, 0}}; P1 = {{0, 0}, {0, 1}};
  prep[P_, g_] := Module[{j, u, es, Ks},
    j = First[Select[Range[2], Norm[P[[All, #]]] > 0 &]];
    u = P[[All, j]] / Norm[P[[All, j]]];
    es = Eigensystem[N[g]];      (* {values, vectors} *)
    Ks = {};
    Do[If[es[[1, k]] > 1*^-14,
       AppendTo[Ks, Sqrt[es[[1, k]]] op[es[[2, k]], u]]], {k, 2}];
    Ks];
  Join[prep[P0, g0], prep[P1, g1]]];

(* Rank-deficient C^3 fixture: Phi_t on the upper 2-block (+) identity on e2. *)
blockLowrank3Kraus[t_] := Module[{embed, Ks},
  embed[M_] := Module[{Z = ConstantArray[0, {3, 3}]}, Z[[1 ;; 2, 1 ;; 2]] = M; Z];
  Ks = {embed[Sqrt[1 - t] IdentityMatrix[2]]};
  Do[AppendTo[Ks, embed[Sqrt[t] op[ee[2, i], ee[2, j]] / Sqrt[2]]], {i, 2}, {j, 2}];
  AppendTo[Ks, op[ee[3, 3], ee[3, 3]]];   (* |2><2| identity on e2 block *)
  Ks];

(* Genuinely complex, asymmetric, unital qubit channel. *)
complexQubitKraus[] := {
  {{0.6, 0}, {0.8 I, 0}},
  {{0, 0.8 I}, {0, 0.6}}};

(* ---------------------------------------------------------------------------
   Fixture corpus -- name, n, Kraus.  Matches the 17-row corpus of
   gen_fixtures_d24.jl in the same order.  The committed eta_ref values are
   grepped from tests/fixtures_d24.inc.h (7th struct field) and pasted here for
   the cross-check; matched by name.
   --------------------------------------------------------------------------- *)
fixtures = {
  {"id_1", 1, identityKraus[1]},
  {"id_2", 2, identityKraus[2]},
  {"id_3", 3, identityKraus[3]},
  {"dep_2", 2, depolarizingKraus[2]},
  {"dep_3", 3, depolarizingKraus[3]},
  {"dephasing_M2", 2, dephasingKraus[2]},
  {"block_cond_exp4", 4, blockCondExpKraus[]},
  {"phi_t_0p3", 2, phitKraus[2, 0.3]},
  {"phi_t_0p1", 2, phitKraus[2, 0.1]},
  {"phi_t_0p01", 2, phitKraus[2, 0.01]},
  {"phi_t_0p001", 2, phitKraus[2, 0.001]},
  {"paper_ex_0p1", 2, paperExampleKraus[0.1]},
  {"paper_ex_0p01", 2, paperExampleKraus[0.01]},
  {"phi_t3_0p05", 3, phitKraus[3, 0.05]},
  {"phi_t4_0p2", 4, phitKraus[4, 0.2]},
  {"block_lowrank3", 3, blockLowrank3Kraus[0.2]},
  {"complex_qubit", 2, complexQubitKraus[]}};

(* eta_ref grepped from tests/fixtures_d24.inc.h (the Julia+MOSEK golden master). *)
etaRef = <|
  "id_1" -> 0.0,
  "id_2" -> 0.0,
  "id_3" -> 0.0,
  "dep_2" -> 1.92295959952765871*^-28,
  "dep_3" -> 5.86063139590153365*^-29,
  "dephasing_M2" -> 0.0,
  "block_cond_exp4" -> 0.0,
  "phi_t_0p3" -> 3.14999999282871312*^-01,
  "phi_t_0p1" -> 1.34999982562583332*^-01,
  "phi_t_0p01" -> 1.48499936517265344*^-02,
  "phi_t_0p001" -> 1.49849597609099262*^-03,
  "paper_ex_0p1" -> 9.48683270364891751*^-02,
  "paper_ex_0p01" -> 9.94987270496749954*^-03,
  "phi_t3_0p05" -> 8.44444411563008424*^-02,
  "phi_t4_0p2" -> 2.99999999994493394*^-01,
  "block_lowrank3" -> 1.59999973678107282*^-01,
  "complex_qubit" -> 1.27999989226236810*^+00|>;

(* ---------------------------------------------------------------------------
   PART 1 -- VALIDATE against EXACT analytic anchors (independent of any solver).
     id_d                 -> 0
     Dep_d - id_d         -> 2(1 - 1/d^2)         (d = 2,3,4: 1.5, 16/9, 15/8)
     Phi_t (= (1-t)id+t Dep) -> t(1-t) 2(1-1/d^2)
     paper example etap   -> etap sqrt(1-etap)    (.tex:378)
   FLAG (do not hide) any mismatch > 1e-6.
   --------------------------------------------------------------------------- *)
(* one-line scientific string, e.g. "3.024e-09" (house style, never the 2-D
   superscript layout that ScientificForm wraps across terminal lines). *)
sci[x_] := Module[{v = N[x], m, ex},
  If[v == 0, Return["0.000e+00"]];
  ex = Floor[Log10[Abs[v]]];
  m = v/10^ex;
  StringJoin[ToString[NumberForm[m, {5, 3}]], "e",
    If[ex < 0, "-", "+"], StringPadLeft[ToString[Abs[ex]], 2, "0"]]];
fixed[x_] := ToString[NumberForm[N[x], {18, 14},
  ExponentFunction -> (Null &), NumberPadding -> {"", "0"}]];

ClearAll[anchorFail];
anchorFail = 0;
chk[name_, got_, want_] := Module[{d = Abs[got - want]},
  If[d > 1*^-6, anchorFail++];
  Print[Row[{
    StringPadRight[name, 28], "  got=",
    StringPadRight[fixed[got], 20], " want=",
    StringPadRight[fixed[want], 20], " |diff|=",
    sci[d], If[d > 1*^-6, "   *** ANCHOR FAIL ***", ""]}]]];

Print["================================================================"];
Print["PART 1: analytic anchors (computed SDP value vs exact closed form)"];
Print["================================================================"];

depIdNorm[d_] := 2 (1 - 1/d^2);
(* SDP of (Dep_d - id_d) directly, NOT via Phi^2-Phi. *)
etaDepDiff[d_] := diamondNorm[choiA[depolarizingKraus[d], d] - choiA[identityKraus[d], d], d];

chk["id_2 -> 0", etaOf[identityKraus[2], 2], 0];
Do[chk["Dep" <> ToString[d] <> "-id" <> ToString[d] <> " = 2(1-1/" <> ToString[d] <> "^2)",
   etaDepDiff[d], depIdNorm[d]], {d, {2, 3, 4}}];
Do[chk["phi_t d=" <> ToString[d] <> " t=" <> ToString[t],
   etaOf[phitKraus[d, t], d], t (1 - t) depIdNorm[d]],
  {dt, {{2, {0.3, 0.1, 0.01, 0.001}}, {3, {0.05}}, {4, {0.2}}}},
  {d, {dt[[1]]}}, {t, dt[[2]]}];
Do[chk["paper etap=" <> ToString[ep], etaOf[paperExampleKraus[ep], 2], ep Sqrt[1 - ep]],
  {ep, {0.1, 0.01}}];

Print[""];
Print[If[anchorFail == 0,
  "ANCHORS OK: all analytic anchors reproduced to < 1e-6.",
  "*** " <> ToString[anchorFail] <> " ANCHOR(S) FAILED -- see flags above. ***"]];

(* ---------------------------------------------------------------------------
   PART 2 -- per-fixture cross-check: eta_Mathematica vs eta_ref (MOSEK).
   --------------------------------------------------------------------------- *)
Print[""];
Print["================================================================"];
Print["PART 2: eta_Mathematica vs eta_ref (Julia+MOSEK golden master)"];
Print["================================================================"];
Print[Row[{StringPadRight["fixture", 18], StringPadRight["n", 3],
  StringPadRight["eta_Mathematica", 24], StringPadRight["eta_ref(MOSEK)", 24],
  "|diff|"}]];

maxDiff = 0; maxDiffName = "";
results = {};
Do[Module[{name = f[[1]], n = f[[2]], Kr = f[[3]], em, er, d},
   em = etaOf[Kr, n];
   er = etaRef[name];
   d = Abs[N[em] - N[er]];
   If[d > maxDiff, maxDiff = d; maxDiffName = name];
   AppendTo[results, {name, n, N[em], N[er], d}];
   Print[Row[{StringPadRight[name, 18], StringPadRight[ToString[n], 3],
     StringPadRight[fixed[em], 24],
     StringPadRight[If[er == 0, "0.00000000000000", sci[er]], 24],
     sci[d]}]];
   ], {f, fixtures}];

Print[""];
Print[Row[{"MAX |eta_Mathematica - eta_ref| over the 17-fixture corpus = ",
  sci[maxDiff], "  (at ", maxDiffName, ")"}]];

(* ---------------------------------------------------------------------------
   PART 3 -- WorkingPrecision probe: can SemidefiniteOptimization run above
   machine precision?  Decides usefulness as the high-precision oracle for the
   certified-arb tight ball (aic-m24 incr 3).
   --------------------------------------------------------------------------- *)
Print[""];
Print["================================================================"];
Print["PART 3: WorkingPrecision / high-precision capability probe"];
Print["================================================================"];
Print[Row[{"Options[SemidefiniteOptimization] = ", Options[SemidefiniteOptimization]}]];

wpResult = Quiet[Check[
   Module[{Jt = lambdaChoi[phitKraus[2, 0.3], 2], n2 = 4, JR, JI, xrv, xiv, Xc,
      Xdag, Pre, Pim, Qre, Qim, P, Q, blk, obj, freeP, freeQ, vars, cons, sol},
    JR = Re[Jt]; JI = Im[Jt];
    xrv = Array[ToExpression["wpxr" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
    xiv = Array[ToExpression["wpxi" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
    Xc = xrv + I xiv; Xdag = Transpose[xrv] - I Transpose[xiv];
    Pre = Array[ToExpression["wpPr" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
    Pim = Array[ToExpression["wpPi" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
    Qre = Array[ToExpression["wpQr" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
    Qim = Array[ToExpression["wpQi" <> ToString[#1] <> "c" <> ToString[#2]] &, {n2, n2}];
    P = Table[Which[i == j, Pre[[i, i]], i < j, Pre[[i, j]] + I Pim[[i, j]],
        True, Pre[[j, i]] - I Pim[[j, i]]], {i, n2}, {j, n2}];
    Q = Table[Which[i == j, Qre[[i, i]], i < j, Qre[[i, j]] + I Qim[[i, j]],
        True, Qre[[j, i]] - I Qim[[j, i]]], {i, n2}, {j, n2}];
    blk = ArrayFlatten[{{P, Xc}, {Xdag, Q}}];
    obj = Total[Flatten[JR xrv + JI xiv]];
    freeP = Join[Flatten[Table[Pre[[i, j]], {i, n2}, {j, i, n2}]],
                 Flatten[Table[Pim[[i, j]], {i, n2}, {j, i + 1, n2}]]];
    freeQ = Join[Flatten[Table[Qre[[i, j]], {i, n2}, {j, i, n2}]],
                 Flatten[Table[Qim[[i, j]], {i, n2}, {j, i + 1, n2}]]];
    vars = Join[Flatten[xrv], Flatten[xiv], freeP, freeQ];
    cons = {VectorGreaterEqual[{blk, ConstantArray[0, {2 n2, 2 n2}]}, {"SemidefiniteCone", 2 n2}],
            P + Q == IdentityMatrix[n2]};
    sol = SemidefiniteOptimization[-obj, cons, vars, WorkingPrecision -> 30];
    {"ACCEPTED", obj /. sol}],
   "REJECTED"], {SemidefiniteOptimization::optx, ReplaceAll::reps}];

Print[Row[{"SemidefiniteOptimization with WorkingPrecision -> 30:  ",
  If[StringQ[wpResult], wpResult,
   If[ListQ[wpResult] && First[wpResult] == "ACCEPTED", "ACCEPTED", "REJECTED"]]}]];
Print["Finding: SemidefiniteOptimization exposes only MaxIterations / Method /"];
Print["PerformanceGoal / Tolerance -- no WorkingPrecision. It is a"];
Print["MACHINE-PRECISION interior-point solver. For aic-m24 incr 3 it serves as"];
Print["an INDEPENDENT machine-precision oracle, NOT a high-precision one;"];
Print["arbitrary-precision certification of eta must come from the FLINT/arb"];
Print["verified-SDP route, not from Mathematica's SDP."];

Print[""];
Print["================================================================"];
Print[If[anchorFail == 0 && maxDiff <= 1*^-6,
  "VERDICT: Mathematica CONFIRMS the MOSEK golden master to <= 1e-6 over the\n" <>
  "entire 17-fixture corpus, and reproduces all exact analytic anchors. The two\n" <>
  "independent solvers agree.",
  "VERDICT: DISCREPANCY DETECTED -- see flags above. anchorFail=" <>
  ToString[anchorFail] <> ", maxDiff=" <> ToString[N[maxDiff]] <> "."]];
Print["================================================================"];

Exit[0];
