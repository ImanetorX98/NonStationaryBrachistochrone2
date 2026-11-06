#!/usr/bin/env wolframscript
(* Independent Wolfram check of the certificate of group 2 (../certificates/group2_Kinf.json): h = G/N on the equator of
   conformal Kerr with constant factor (eq:kerr-G), ZAMO frame of h, Ricci tensor computed from scratch. The numerators
   P_U, P_V, P_J^2 are DERIVED here from the curvature (and shown to be polynomials); P_U and P_J agree with the formulas
   of the text (proposition on K_inf); the six published polynomials PV, P3, P4, P5, mPU, P6 equal the combinations of
   P_U, P_V, P_J^2 in the declared box map; the conversion factors are positive; and each of the three sufficient sets
   implies max_s Ric(k,k) <= 0. Does not use the Python code of the project. M = 1.
   Usage: wolframscript -file wolfram_bridge_group2.wls [../certificates] *)
cdir = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], FileNameJoin[{DirectoryName[$InputFileName], "..", "certificates"}]];
pass = 0; fail = 0;
ok[name_, c_] := (If[TrueQ[c], pass++; Print["PASS ", name], fail++; Print["FAIL ", name]]);
co = {eta, r, ph};
f = 1 - 2/r; Del = r^2 - 2 r + a^2; Nn = 1/(1 - f u);
Gpp = r^2 + a^2 (1 + 2/r) + u Del/(1 - f u);
G = {{-f, 0, -2 a/r}, {0, Nn r^2/Del, 0}, {-2 a/r, 0, Gpp}};
h = G/Nn; hi = Simplify[Inverse[h]];
Gam = Table[Simplify[1/2 Sum[hi[[i, l]] (D[h[[l, j]], co[[k]]] + D[h[[l, k]], co[[j]]] - D[h[[j, k]], co[[l]]]), {l, 3}]], {i, 3}, {j, 3}, {k, 3}];
Ric = Table[Simplify[Sum[D[Gam[[i, j, k]], co[[i]]] - D[Gam[[i, j, i]], co[[k]]] +
      Sum[Gam[[i, i, l]] Gam[[l, j, k]] - Gam[[i, k, l]] Gam[[l, j, i]], {l, 3}], {i, 3}]], {j, 3}, {k, 3}];
om = -h[[1, 3]]/h[[3, 3]]; al2 = -(h[[1, 1]] - h[[1, 3]]^2/h[[3, 3]]);
e0 = {1, 0, om}/Sqrt[al2]; e1 = {0, 1, 0}/Sqrt[h[[2, 2]]]; e2 = {0, 0, 1}/Sqrt[h[[3, 3]]];
as = {r > 1 + Sqrt[1 - a^2], 0 < a < 1, 0 < u < 1};
ok["ZAMO frame of h orthonormal", Simplify[{e0.h.e0, e1.h.e1, e2.h.e2, e0.h.e2, e0.h.e1, e1.h.e2}, as] === {-1, 1, 1, 0, 0, 0}];
ok["Ric(e0, e1) = Ric(e1, e2) = 0", Simplify[{e0.Ric.e1, e1.Ric.e2}, as] === {0, 0}];
(* k = e0 + c e1 + s e2: Ric(k,k) = U c^2 + V s^2 + J s *)
Um = e0.Ric.e0 + e1.Ric.e1; Vm = e0.Ric.e0 + e2.Ric.e2; Jm = 2 e0.Ric.e2;
(* variables of the text: R = r_+, a^2 = R(2-R), x = R/r *)
wt = R (1 - u) + 2 u x; DD = (1 - x) (R (1 + x) - 2 x);
Hh = (r^4 + a^2 r^2 + 2 a^2 r + 4 a^2 u) x^4/R;
sub[e_] := Simplify[Together[(e /. a -> Sqrt[R (2 - R)]) /. r -> R/x], {0 < x < 1, 1 < R < 2, 0 < u < 1}];
PUm = sub[Um R^3 wt Hh/(x^3 DD)]; PVm = sub[-Vm R^4 wt^2 Hh/x^3]; PJ2m = sub[Jm^2 R^7 wt^3 Hh^2/(4 x^8 (2 - R) DD)];
ok["P_U = U R^3 w~ H/(x^3 D) is a polynomial in (x, R, u)", PolynomialQ[PUm, {x, R, u}]];
ok["P_V = -V R^4 w~^2 H/x^3 is a polynomial in (x, R, u)", PolynomialQ[PVm, {x, R, u}]];
ok["P_J^2 = J^2 R^7 w~^3 H^2/(4 x^8 (2-R) D) is a polynomial in (x, R, u)", PolynomialQ[PJ2m, {x, R, u}]];
(* formulas of the text *)
PU = (-3*R^3*u*x^2 - R^3*u + 3*R^3*x^2 - 16*R^2*u^2*x^3 + 21*R^2*u*x^3 + 6*R^2*u*x^2 - 6*R^2*x^2 + 28*R*u^2*x^4 + 32*R*u^2*x^3 - 42*R*u*x^3 - 56*u^2*x^4);
PJ = (3*R^4*u^2*x^2 - 3*R^4*u^2 - 6*R^4*u*x^2 + 6*R^4*u + 3*R^4*x^2 - 3*R^4 - 21*R^3*u^2*x^3 - 6*R^3*u^2*x^2 + 19*R^3*u^2*x + 21*R^3*u*x^3 + 12*R^3*u*x^2 - 19*R^3*u*x - 6*R^3*x^2 - 20*R^2*u^3*x^4 + 48*R^2*u^2*x^4 + 42*R^2*u^2*x^3 - 24*R^2*u^2*x^2 - 42*R^2*u*x^3 + 36*R*u^3*x^5 + 40*R*u^3*x^4 - 96*R*u^2*x^4 - 72*u^3*x^5);
ok["P_U equals the formula of the text", Expand[PUm - PU] === 0];
ok["P_J^2 equals the square of the formula of the text", Expand[PJ2m - PJ^2] === 0];
ok["conversion factors w~, D, H positive for (x, R, u) in (0,1) x (1,2) x [0,1)",
  Resolve[ForAll[{x, R, u}, 0 < x < 1 && 1 < R < 2 && 0 <= u < 1, wt > 0 && DD > 0], Reals] &&
  Resolve[ForAll[{rr, aa, uu}, rr > 1 && 0 < aa < 1 && 0 <= uu < 1, rr^4 + aa^2 rr^2 + 2 aa^2 rr + 4 aa^2 uu > 0], Reals]];
(* combinations: V^2 - J^2 = x^6 P3/(R^8 w~^4 H^2), -4U(U-V) - J^2 ~ P4, V - U ~ P5, U^2 - J^2 ~ P6 (positive factors) *)
P3 = Expand[PVm^2 - 4 R (2 - R) x^2 DD wt PJ2m];
P4 = Expand[-PUm (R wt DD PUm + PVm) - x^2 (2 - R) PJ2m];
P5 = Expand[-(PVm + R wt DD PUm)]; P6 = Expand[R wt DD PUm^2 - 4 x^2 (2 - R) PJ2m];
Uf = x^3 DD PUm/(R^3 wt Hs); Vf = -x^3 PVm/(R^4 wt^2 Hs); J2f = 4 x^8 (2 - R) DD PJ2m/(R^7 wt^3 Hs^2);
ratio[e_, p_] := Factor[Together[e/p]];
ok["V^2 - J^2 = x^6 P3/(R^8 w~^4 H^2)", Simplify[ratio[Vf^2 - J2f, P3] - x^6/(R^8 wt^4 Hs^2)] === 0];
ok["-4U(U-V) - J^2 = 4 x^6 D P4/(R^7 w~^3 H^2)", Simplify[ratio[-4 Uf (Uf - Vf) - J2f, P4] - 4 x^6 DD/(R^7 wt^3 Hs^2)] === 0];
ok["V - U = x^3 P5/(R^4 w~^2 H)", Simplify[ratio[Vf - Uf, P5] - x^3/(R^4 wt^2 Hs)] === 0];
ok["U^2 - J^2 = x^6 D P6/(R^7 w~^3 H^2)", Simplify[ratio[Uf^2 - J2f, P6] - x^6 DD/(R^7 wt^3 Hs^2)] === 0];
(* published polynomials in the box map x = X, R = 1 + Rv, u = (3/5) U *)
cc = Import[FileNameJoin[{cdir, "group2_Kinf.json"}], "RawJSON"];
box = {x -> X, R -> 1 + Rv, u -> 3/5 Uv};
att = <|"PV" -> PVm, "P3" -> P3, "P4" -> P4, "P5" -> P5, "mPU" -> -PUm, "P6" -> P6|>;
Do[Module[{p = Total[(ToExpression[#[[2]]] Times @@ ({X, Rv, Uv}^#[[1]])) & /@ cc["polynomials"][n]]},
   ok["group2_Kinf, " <> n <> ": JSON polynomial = derived one in the box map", Expand[p - (att[n] /. box)] === 0]],
  {n, cc["polynomial_order"]}];
(* the three sufficient sets, each with V <= 0 and V^2 >= J^2, imply max_s N <= 0, N = U(1-s^2) + V s^2 + J s *)
Nk[s_] := Ak (1 - s^2) + Bk s^2 + Jv s; base = Bk <= 0 && Bk^2 >= Jv^2;
ok["set 0 (value at the vertex): J^2 <= -4U(U-V) => max_s N <= 0", Resolve[ForAll[{Ak, Bk, Jv, s}, (base && Jv^2 <= -4 Ak (Ak - Bk) && -1 <= s <= 1) \[Implies] Nk[s] <= 0], Reals]];
ok["set 1: V >= U => max_s N <= 0", Resolve[ForAll[{Ak, Bk, Jv, s}, (base && Bk >= Ak && -1 <= s <= 1) \[Implies] Nk[s] <= 0], Reals]];
ok["set 2: U <= 0, J^2 <= U^2 => max_s N <= 0", Resolve[ForAll[{Ak, Bk, Jv, s}, (base && Ak <= 0 && Jv^2 <= Ak^2 && -1 <= s <= 1) \[Implies] Nk[s] <= 0], Reals]];
Print["SUMMARY: ", pass, " PASS, ", fail, " FAIL"];
