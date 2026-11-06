(* Paper III, symbolic checks in Wolfram Language: common functions.
   chk[name, expr, assumptions]: PASS if FullSimplify[expr, assumptions] == 0 (or True);
   chkFalse[name, expr, assumptions]: a check that MUST fail (PASS if expr does not vanish; FAIL if undecided);
   chkN[name, val, tol]: numerical check |val| < tol. *)
$fails = {}; $npass = 0;
simplifyToZero[e_, ass_] := Module[{s = Quiet@FullSimplify[e, ass]},
   If[s === 0 || s === True, s, Quiet@FullSimplify[Together[TrigToExp[s]], ass]]];
chk[name_, e_, ass_: True] := Module[{s = TimeConstrained[simplifyToZero[e, ass], 300, $Aborted]},
   If[s === 0 || s === True, $npass++; Print["  [PASS] ", name],
     AppendTo[$fails, name]; Print["  [FAIL] ", name, "   residual: ", Short[s, 3]]]];
chkFalse[name_, e_, ass_: True] := Module[{s = TimeConstrained[simplifyToZero[e, ass], 300, $Aborted]},
   Which[s === 0 || s === True, AppendTo[$fails, name]; Print["  [FAIL] (should have failed) ", name],
     s === $Aborted, AppendTo[$fails, name]; Print["  [FAIL] (undecided: simplification aborted) ", name],
     True, $npass++; Print["  [PASS] (fails as expected) ", name]]];
chkN[name_, v_, tol_] := If[Abs[v] < tol, $npass++; Print["  [PASS] ", name, "   (", ScientificForm[N[Abs[v]], 3], ")"],
   AppendTo[$fails, name]; Print["  [FAIL] ", name, "   value: ", N[v]]];
summary[] := (Print["\nSUMMARY: ", $npass, " PASS, ", Length[$fails], " FAIL"];
   If[$fails =!= {}, Print["FAILED: ", $fails]]);
(* geometry: Christoffel symbols and Riemann R^i_{jkl} = d_k G^i_{lj} - d_l G^i_{kj} + G^i_{km}G^m_{lj} - G^i_{lm}G^m_{kj}.
   NOTE: the summation indices have reserved names (ii, jj, kk, ll, mm, nn) so that they do not capture the symbols of
   the metrics (indices called a, b, c once captured the spin a). *)
christoffel[g_, x_] := Module[{gi = Simplify[Inverse[g]], n = Length[x]},
   Table[Simplify[1/2 Sum[gi[[ii, ll]] (D[g[[ll, jj]], x[[kk]]] + D[g[[ll, kk]], x[[jj]]] - D[g[[jj, kk]], x[[ll]]]), {ll, n}]],
     {ii, n}, {jj, n}, {kk, n}]];
riemann[g_, x_] := Module[{G = christoffel[g, x], n = Length[x]},
   Table[D[G[[ii, ll, jj]], x[[kk]]] - D[G[[ii, kk, jj]], x[[ll]]] +
      Sum[G[[ii, kk, mm]] G[[mm, ll, jj]] - G[[ii, ll, mm]] G[[mm, kk, jj]], {mm, n}], {ii, n}, {jj, n}, {kk, n}, {ll, n}]];
ricci[g_, x_] := Module[{Rm = riemann[g, x], n = Length[x]}, Table[Sum[Rm[[kk, ii, kk, jj]], {kk, n}], {ii, n}, {jj, n}]];
(* R(X,Y,Z,W) = g_{im} R^m_{jkl} X^i Y^j Z^k W^l *)
rq[g_, Rm_, X_, Y_, Z_, W_] := Module[{n = Length[g]}, Sum[g[[ii, mm]] Rm[[mm, jj, kk, ll]] X[[ii]] Y[[jj]] Z[[kk]] W[[ll]],
   {ii, n}, {mm, n}, {jj, n}, {kk, n}, {ll, n}]];
(* internal test of the machinery: Schwarzschild is Ricci-flat *)
autotest[] := Module[{t, rr, th, ph, MM, gS, Ric},
   gS = DiagonalMatrix[{-(1 - 2 MM/rr), 1/(1 - 2 MM/rr), rr^2, rr^2 Sin[th]^2}];
   Ric = Simplify[ricci[gS, {t, rr, th, ph}]];
   Ric === ConstantArray[0, {4, 4}]];
