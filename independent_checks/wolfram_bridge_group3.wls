#!/usr/bin/env wolframscript
(* Independent Wolfram check of the main certificates of group 3 (../certificates/group3_{I2,I1plus,I1minus}_main.json):
   gbar_inf = -(1 - f u(eta)) (deta - beta)^2 + h0 on the equator of conformal Kerr with a time-dependent factor, Ricci
   tensor computed from scratch, k = e0 + c e_r + s e_phi.
   (1) Ric(k,k) = R_stat + D1 u' + D2 [u'' + f u'^2/(2w)] exactly (no other jet monomials), with
       R_stat = wN/w, D1 = c (al0 + al1 s)/w, D2 = -c^2 be/w (formulas of the text, exact identities);
   (2) the published P0 + P1 sqrt(g) of each main certificate equals a positive factor times the inequality of the text
       (I2: (1/3) w(-N) - c^2 be; I1 plus/minus: v (5/4 w(-N) -/+ c(al0 + al1 s)), v = sqrt(w S)), and g = w S.
   The delegated sub-certificates C1, C2A, C2B are composed by the verifier (composition data checked there).
   Does not use the Python code of the project. M = 1.
   Usage: wolframscript -file wolfram_bridge_group3.wls [../certificates] *)
cdir = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], FileNameJoin[{DirectoryName[$InputFileName], "..", "certificates"}]];
pass = 0; fail = 0;
ok[name_, c_] := (If[TrueQ[c], pass++; Print["PASS ", name], fail++; Print["FAIL ", name]]);

Print["--- (1) curvature and jets"];
co = {eta, r, ph};
f = 1 - 2/r; Del = r^2 - 2 r + a^2; V2 = 1 - f uu[eta];
hrr = r^2/(f Del); hpp = Del/f^2; bp = -2 a/(r f);
g = {{-V2, 0, V2 bp}, {0, hrr, 0}, {V2 bp, 0, hpp - V2 bp^2}};
gi = Simplify[Inverse[g]];
Gam = Table[Simplify[1/2 Sum[gi[[i, l]] (D[g[[l, j]], co[[k]]] + D[g[[l, k]], co[[j]]] - D[g[[j, k]], co[[l]]]), {l, 3}]], {i, 3}, {j, 3}, {k, 3}];
Ric = Table[Simplify[Sum[D[Gam[[i, j, k]], co[[i]]] - D[Gam[[i, j, i]], co[[k]]] +
      Sum[Gam[[i, i, l]] Gam[[l, j, k]] - Gam[[i, k, l]] Gam[[l, j, i]], {l, 3}], {i, 3}]], {j, 3}, {k, 3}];
e0 = {1/Sqrt[V2], 0, 0}; er = {0, 1/Sqrt[hrr], 0}; ep = {bp, 0, 1}/Sqrt[hpp];
jet = {uu''[eta] -> u2, uu'[eta] -> u1, uu[eta] -> u0};
Kv = e0 + cc er + ss ep;
RicK = Simplify[(Kv.Ric.Kv) /. jet];
Rst = Simplify[RicK /. {u1 -> 0, u2 -> 0}];
c2 = Simplify[Coefficient[RicK, u2]]; c11 = Simplify[Coefficient[RicK, u1, 2]];
c1 = Simplify[Coefficient[RicK, u1, 1] /. u2 -> 0];
w0 = 1 - f u0;
ok["Ric(k,k) = R_stat + D1 u' + D3 u'^2 + D2 u'' exactly (no other jet monomials)", Simplify[RicK - (Rst + c1 u1 + c11 u1^2 + c2 u2)] === 0];
ok["grouping: D3 = D2 f/(2w)", Simplify[c11 - c2 f/(2 w0)] === 0];
(* formulas of the text: x = 2/r, z = a^2, y^2 = 1 - x, v = sqrt(w S) *)
x = 2/r; z = a^2; ww = 1 - u0 + u0 x; SS = x^2 z - 4 x + 4; y2 = 1 - x; vv = Sqrt[ww SS];
nU = (2*u0^2*x^5*z - 4*u0^2*x^4*z + 2*u0^2*x^3*z + 6*u0*x^4*z - 12*u0*x^3*z + 6*u0*x^2*z + 4*u0*x^2 - 8*u0*x + 4*u0 + 5*x^3*z - 6*x^2*z - 6*x^2 + 14*x - 8);
nV = (4*u0^3*x^6*z - 12*u0^3*x^5*z + 12*u0^3*x^4*z - 4*u0^3*x^3*z + 24*u0^2*x^5*z - 66*u0^2*x^4*z - 32*u0^2*x^4 + 60*u0^2*x^3*z + 128*u0^2*x^3 - 18*u0^2*x^2*z - 192*u0^2*x^2 + 128*u0^2*x - 32*u0^2 + 31*u0*x^4*z - 61*u0*x^3*z - 48*u0*x^3 + 30*u0*x^2*z + 144*u0*x^2 - 144*u0*x + 48*u0 + 10*x^3*z - 12*x^2*z - 12*x^2 + 28*x - 16);
LL = 8 u0 x^2 - 14 u0 x + 6 u0 + 5 x - 6;
TT = 3*z*u0*x^4 - 4*z*u0*x^3 + z*u0*x^2 + 4*z*x^3 - z*x^2 - 16*u0*x^3 + 44*u0*x^2 - 40*u0*x + 12*u0 - 20*x^2 + 32*x - 12;
wU = x^3 nU/(32 y2); wV0 = x^3 nV/(64 y2 ww); wJ = -a x^4 LL vv/(32 y2);
wN = wU cc^2 + wV0 ss^2 + wJ ss;
al0 = 3 z x^5 Sqrt[y2] ww/(8 vv); al1 = -a x^3 Sqrt[y2] TT/(8 SS ww); be = z x^4 y2/(2 SS);
as = {r > 2, 0 < a < 1, 0 < u0 < 1/2, -1 < ss < 1, -1 < cc < 1};
ok["R_stat = w N / w (exact)", FullSimplify[(Rst - wN/ww) /. cc^2 -> 1 - ss^2, as] === 0];
ok["D1 = c (al0 + al1 s)/w (exact)", FullSimplify[c1 - cc (al0 + al1 ss)/ww, as] === 0];
ok["D2 = -c^2 be/w (exact)", FullSimplify[(c2 + cc^2 be/ww) /. cc^2 -> 1 - ss^2, as] === 0];

Print["--- (2) published main certificates"];
Clear[x, z, ww, SS, nU, nV, LL, TT, wU, wV0, wJ, wN, al0, al1, be];
(* rationalising variables of the certificates: y = sqrt(1 - x), t = tan(alpha/2), s = 2t/(1+t^2), c = (1-t^2)/(1+t^2);
   v = sqrt(w S) kept as a symbol *)
x = 1 - y^2; z = a^2; ww = 1 - u + u x; SS = x^2 z - 4 x + 4;
nU = (2*u^2*x^5*z - 4*u^2*x^4*z + 2*u^2*x^3*z + 6*u*x^4*z - 12*u*x^3*z + 6*u*x^2*z + 4*u*x^2 - 8*u*x + 4*u + 5*x^3*z - 6*x^2*z - 6*x^2 + 14*x - 8);
nV = (4*u^3*x^6*z - 12*u^3*x^5*z + 12*u^3*x^4*z - 4*u^3*x^3*z + 24*u^2*x^5*z - 66*u^2*x^4*z - 32*u^2*x^4 + 60*u^2*x^3*z + 128*u^2*x^3 - 18*u^2*x^2*z - 192*u^2*x^2 + 128*u^2*x - 32*u^2 + 31*u*x^4*z - 61*u*x^3*z - 48*u*x^3 + 30*u*x^2*z + 144*u*x^2 - 144*u*x + 48*u + 10*x^3*z - 12*x^2*z - 12*x^2 + 28*x - 16);
LL = 8 u x^2 - 14 u x + 6 u + 5 x - 6;
TT = 3*z*u*x^4 - 4*z*u*x^3 + z*u*x^2 + 4*z*x^3 - z*x^2 - 16*u*x^3 + 44*u*x^2 - 40*u*x + 12*u - 20*x^2 + 32*x - 12;
s = 2 t/(1 + t^2); c = (1 - t^2)/(1 + t^2);
wU = x^3 nU/(32 y^2); wV0 = x^3 nV/(64 y^2 ww); wJ = -a x^4 LL v/(32 y^2);
mwN = -(wU c^2 + wV0 s^2 + wJ s);
cal = c 3 z x^5 y ww v/(8 ww SS) + (-c s a x^3 y TT/(8 SS ww));   (* c (al0 + al1 s), using 1/v = v/(w S) *)
be = z x^4 y^2/(2 SS);
ineq = <|"group3_I2_main" -> mwN/3 - c^2 be, "group3_I1plus_main" -> v (5/4 mwN - cal), "group3_I1minus_main" -> v (5/4 mwN + cal)|>;
split[e_] := Module[{q = Together[e], nn}, nn = Expand[Numerator[q]] /. v^2 -> ww SS; nn = Expand[nn] /. v^3 -> ww SS v;
   nn = Expand[nn] /. v^2 -> ww SS; {Coefficient[nn, v, 0]/Denominator[q], Coefficient[nn, v, 1]/Denominator[q]}];
box = {y -> Y, t -> 2 T - 1, u -> 2/5 Uv, a -> 33/50 A};
rho[k_] := k (1 - 2 T + 2 T^2)^2 Y^2 (2 Uv Y^2 - 5) (1089 A^2 (1 - Y^2)^2 + 10^4 Y^2)/(3125 (Y^2 - 1)^3);
factor = <|"group3_I2_main" -> rho[96], "group3_I1plus_main" -> rho[128], "group3_I1minus_main" -> rho[128]|>;
Do[Module[{cc = Import[FileNameJoin[{cdir, n <> ".json"}], "RawJSON"], P, A0, A1, r0, r1},
   P[k_] := Total[(ToExpression[#[[2]]] Times @@ ({Y, T, Uv, A}^#[[1]])) & /@ cc["polynomials"][k]];
   {A0, A1} = split[ineq[n]] /. box;
   r0 = Together[P["P0"]/A0]; r1 = Together[P["P1"]/A1];
   ok[n <> ": P0/A0 = P1/A1 (one common factor)", Together[r0 - r1] === 0];
   ok[n <> ": the common factor is " <> If[n === "group3_I2_main", "96", "128"] <> " (1-2T+2T^2)^2 Y^2 (2UY^2-5)(1089A^2(1-Y^2)^2+10^4Y^2)/(3125(Y^2-1)^3)",
     Together[r0 - factor[n]] === 0];
   ok[n <> ": g_json = w S in the box map", Expand[P["g"] - (ww SS /. box)] === 0]],
  {n, Keys[ineq]}];
ok["the common factors are positive on the open box: 2UY^2 - 5 < 0 and (Y^2 - 1)^3 < 0 for 0 <= U <= 1, 0 < Y < 1; the other factors are positive",
  Resolve[ForAll[{Uv, Y}, 0 <= Uv <= 1 && 0 < Y < 1, 2 Uv Y^2 - 5 < 0 && (Y^2 - 1)^3 < 0], Reals] &&
  Resolve[ForAll[{T, A, Y}, 0 <= T <= 1 && 0 <= A <= 1 && 0 < Y < 1, 1 - 2 T + 2 T^2 > 0 && 1089 A^2 (1 - Y^2)^2 + 10^4 Y^2 > 0], Reals]];
Print["SUMMARY: ", pass, " PASS, ", fail, " FAIL"];
