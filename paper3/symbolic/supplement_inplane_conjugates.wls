#!/usr/bin/env wolframscript
(* supplement_inplane_conjugates.wls  --  Paper III, supplement (proposition on the certificates, frames)
   a = 0, equatorial 2+1, g_opt = -deta^2 + N(eta,r) h0, N = 1/(1 - A(eta)^2 f/Ehat^2).
   In-plane null Jacobi equation y'' + Ric(k,k) y = 0 (3D: Weyl = 0, screen 1-dim);
   null conjugate points are conformally invariant, so any frame may be used.
   Checks:
     1. frame ghat = g_opt/N: Ric(k,k) contains no eta-derivative of A;
     2. same for every frame e^{2w(r)} ghat;
     3. CONTROL (must fail): in the frame g_opt itself A' and A'' do appear;
     4. at constant A the ultrastatic frame gives exactly the optical Gauss curvature K_t of Paper II;
     5. ghat frame: Ric(k,k) <= 0 for all r > 2M and all directions  <=>  E_eff^2 >= 2 (Reduce);
     6. diagnostic grid: windows of E_eff^2 for reference frames w = half log of N(r) at reference energy Es. *)
$Failures = 0;
check[label_, cond_] := (Print[If[TrueQ[cond], "PASS  ", "FAIL  "], label]; If[! TrueQ[cond], $Failures++]);
x = {eta, r, ph}; n = 3;
ricci[g_] := Module[{gi, G},
  gi = Simplify[Inverse[g]];
  G = Table[Simplify[(1/2) Sum[gi[[a, d]] (D[g[[d, b]], x[[c]]] + D[g[[d, c]], x[[b]]] - D[g[[b, c]], x[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  Table[Simplify[Sum[D[G[[a, b, c]], x[[a]]] - D[G[[a, b, a]], x[[c]]] + Sum[G[[a, a, d]] G[[d, b, c]] - G[[a, c, d]] G[[d, b, a]], {d, n}], {a, n}]], {b, n}, {c, n}]];
nullk[g_] := Module[{gtt = -g[[1, 1]]}, {1, Cos[al] Sqrt[gtt/g[[2, 2]]], Sin[al] Sqrt[gtt/g[[3, 3]]]}];
hasAdot[ex_] := ! FreeQ[ex, Derivative[1][A]] || ! FreeQ[ex, Derivative[2][A]];
f = 1 - 2/r;  (* M = 1 *)
Nn = 1/(1 - A[eta]^2 f/Eh^2);
gopt = DiagonalMatrix[{-1, Nn/f^2, Nn r^2/f}];
ghat = DiagonalMatrix[{-1/Nn, 1/f^2, r^2/f}];
gw = Exp[2 w[r]] ghat;
Rhat = Simplify[nullk[ghat] . ricci[ghat] . nullk[ghat]];
check["1. ghat = g_opt/N: Ric(k,k) free of A', A''", ! hasAdot[Rhat]];
Rw = Simplify[nullk[gw] . ricci[gw] . nullk[gw]];
check["2. e^{2w(r)} ghat: Ric(k,k) free of A', A''", ! hasAdot[Rw]];
Ropt = Simplify[nullk[gopt] . ricci[gopt] . nullk[gopt]];
check["3. CONTROL: in g_opt itself A' or A'' appear", hasAdot[Ropt]];
NE = 1/(1 - f/Ee^2);
gU = DiagonalMatrix[{-1, NE/f^2, NE r^2/f}];
Kt = (16 + 2*(9*Ee^2 - 14)*r + (3*Ee^4 - 19*Ee^2 + 16)*r^2 - (2*Ee^2 - 3)*(Ee^2 - 1)*r^3)/(Ee^2 r^5 (2 + (Ee^2 - 1) r));
check["4. constant A, ultrastatic frame: Ric(k,k) = K_t (optical Gauss curvature of Paper II)",
  Simplify[nullk[gU] . ricci[gU] . nullk[gU] - Kt] === 0];
RhE = Simplify[Rhat /. A[eta] -> Eh/Ee];
tang = Simplify[RhE /. al -> Pi/2]; rad = Simplify[RhE /. al -> 0];
check["5a. ghat: Ric(k,k) is affine in sin^2(al) (extremes at radial/tangential)",
  Simplify[D[RhE /. {Cos[al] -> Sqrt[1 - s], Sin[al] -> Sqrt[s]} /. Cos[2 al] -> 1 - 2 s, {s, 2}]] === 0];
check["5b. ghat: radial direction Ric(k,k) < 0 for all r > 2M, Ee > 1",
  Reduce[Exists[r, r > 2 && Ee > 1 && rad >= 0], Reals] === False];
check["5c. ghat: tangential Ric(k,k) <= 0 for all r > 2M  <=>  Ee^2 >= 2",
  Reduce[ForAll[r, r > 2, tang <= 0] && Ee > 1, Ee, Reals] === (Ee >= Sqrt[2])];
(* 6. diagnostic: reference frames *)
Rref = Simplify[Rw /. A[eta] -> Eh/Ee /. w -> Function[rr, Log[1/(1 - (1 - 2/rr)/Es^2)]/2]];
F = Compile[{{e, _Real}, {es, _Real}, {rr, _Real}, {a, _Real}}, Evaluate[Rref /. {Ee -> e, Es -> es, r -> rr, al -> a}]];
rs = Join[Range[2.0005, 3, 0.005], Range[3, 30, 0.1], {40, 60, 100, 300, 1000, 10^4}];
as = Range[0, Pi, Pi/24.];
mx[e_, es_] := Max[Flatten[Table[F[e, es, rr, a], {rr, rs}, {a, as}]]];
check["6a. sanity: reference frame at Ee = Ee*, Ee^2 = 1.4 has Ric(k,k) > 0 somewhere", mx[Sqrt[1.4], Sqrt[1.4]] > 0];
Do[okE = Select[Range[1.01, 3.0, 0.01], mx[#, es] <= 0 &];
  Print["      diagnostic: Ee*^2 = ", es^2, " -> window E_eff^2 in ", If[okE === {}, "none", {Min[okE]^2, Max[okE]^2}], " (grid top 9)"],
  {es, {1.25, 1.3, Sqrt[2], 1.5, 2.}}];
Print[]; If[$Failures == 0, Print["ALL CHECKS PASSED"], Print[$Failures, " CHECK(S) FAILED"]];
Exit[If[$Failures == 0, 0, 1]];
