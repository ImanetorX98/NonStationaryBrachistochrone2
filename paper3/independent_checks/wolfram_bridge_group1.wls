#!/usr/bin/env wolframscript
(* Independent Wolfram check of the certificates of group 1 (../certificates/group1_*.json): the curvature of the
   equatorial 2+1 metric gbar_inf = -V^2 (deta - beta)^2 + h0 of conformal Kerr with constant factor (q = Ehat^2/A^2
   constant), computed from scratch (own Christoffel symbols and Ricci tensor, own orthonormal frame), gives the
   quantities U, V, J of the text; the four published polynomials are -nU, -nV, P3, P6 in the declared box map; their
   Bernstein coefficients on the box are non-negative (each certificate is a single leaf); and the four signs imply
   max_s Ric(k,k) <= 0. Does not use the Python code of the project. M = 1.
   Usage: wolframscript -file wolfram_bridge_group1.wls [../certificates] *)
cdir = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], FileNameJoin[{DirectoryName[$InputFileName], "..", "certificates"}]];
pass = 0; fail = 0;
ok[name_, c_] := (If[TrueQ[c], pass++; Print["PASS ", name], fail++; Print["FAIL ", name]]);
load[name_] := Import[FileNameJoin[{cdir, name <> ".json"}], "RawJSON"];
poly[cc_, k_, vars_] := Total[(ToExpression[#[[2]]] Times @@ (vars^#[[1]])) & /@ cc["polynomials"][k]];

co = {eta, r, ph};
f = 1 - 2/r; Del = r^2 - 2 r + a^2; V2 = 1 - f/q;
hrr = r^2/(f Del); hpp = Del/f^2; bp = -2 a/(r f);
g = {{-V2, 0, V2 bp}, {0, hrr, 0}, {V2 bp, 0, hpp - V2 bp^2}};
gi = Simplify[Inverse[g]];
Gam = Table[Simplify[1/2 Sum[gi[[i, l]] (D[g[[l, j]], co[[k]]] + D[g[[l, k]], co[[j]]] - D[g[[j, k]], co[[l]]]), {l, 3}]], {i, 3}, {j, 3}, {k, 3}];
Ric = Table[Simplify[Sum[D[Gam[[i, j, k]], co[[i]]] - D[Gam[[i, j, i]], co[[k]]] +
      Sum[Gam[[i, i, l]] Gam[[l, j, k]] - Gam[[i, k, l]] Gam[[l, j, i]], {l, 3}], {i, 3}]], {j, 3}, {k, 3}];
e0 = {1/Sqrt[V2], 0, 0}; er = {0, 1/Sqrt[hrr], 0}; ep = {bp, 0, 1}/Sqrt[hpp];
as = {r > 2, q > 1, 0 < a < 1};
ok["orthonormal frame e0, e_r, e_phi", Simplify[{e0.g.e0, er.g.er, ep.g.ep, e0.g.ep, e0.g.er, er.g.ep}, as] === {-1, 1, 1, 0, 0, 0}];
ok["Ric(e0, e_r) = Ric(e_r, e_phi) = 0", Simplify[{e0.Ric.er, er.Ric.ep}, as] === {0, 0}];
(* k = e0 + c e_r + s e_phi: Ric(k,k) = U c^2 + V s^2 + J s *)
Um = e0.Ric.e0 + er.Ric.er; Vm = e0.Ric.e0 + ep.Ric.ep; Jm = 2 e0.Ric.ep;
(* quantities of the text: x = 2/r, u = 1/q, z = a^2, w = 1 - u + u x *)
rules = {x -> 2/r, u -> 1/q, z -> a^2};
w = 1 - u + u x; S = x^2 z - 4 x + 4; L = 8 u x^2 - 14 u x + 6 u + 5 x - 6;
nU = (2*u^2*x^5*z - 4*u^2*x^4*z + 2*u^2*x^3*z + 6*u*x^4*z - 12*u*x^3*z + 6*u*x^2*z + 4*u*x^2 - 8*u*x + 4*u + 5*x^3*z - 6*x^2*z - 6*x^2 + 14*x - 8);
nV = (4*u^3*x^6*z - 12*u^3*x^5*z + 12*u^3*x^4*z - 4*u^3*x^3*z + 24*u^2*x^5*z - 66*u^2*x^4*z - 32*u^2*x^4 + 60*u^2*x^3*z + 128*u^2*x^3 - 18*u^2*x^2*z - 192*u^2*x^2 + 128*u^2*x - 32*u^2 + 31*u*x^4*z - 61*u*x^3*z - 48*u*x^3 + 30*u*x^2*z + 144*u*x^2 - 144*u*x + 48*u + 10*x^3*z - 12*x^2*z - 12*x^2 + 28*x - 16);
Ut = x^3 nU/(32 (1 - x) w); Vt = x^3 nV/(64 (1 - x) w^2); J2t = z x^8 S L^2/(1024 (1 - x)^2 w);
ok["U = x^3 nU/(32 (1-x) w)", Simplify[Um - (Ut /. rules), as] === 0];
ok["V = x^3 nV/(64 (1-x) w^2)", Simplify[Vm - (Vt /. rules), as] === 0];
ok["J^2 = z x^8 S L^2/(1024 (1-x)^2 w)", Simplify[Jm^2 - (J2t /. rules), as] === 0];
K = z x^2 S L^2; P3 = Expand[nV^2 - 4 K w^3]; P6 = Expand[4 nU^2 - K w];
ok["V^2 - J^2 = x^6 P3/(4096 (1-x)^2 w^4)", Simplify[(Vt^2 - J2t) - x^6 P3/(4096 (1 - x)^2 w^4)] === 0];
ok["4U^2 - J^2 = x^6 P6/(1024 (1-x)^2 w^2)", Simplify[(4 Ut^2 - J2t) - x^6 P6/(1024 (1 - x)^2 w^2)] === 0];
(* published polynomials in the box map x = X, u = U/2, z = (1089/2500) Z *)
box = {x -> X, u -> Uv/2, z -> 1089/2500 Z};
expected = {"group1_minus_nU" -> -nU, "group1_minus_nV" -> -nV, "group1_P3" -> P3, "group1_P6" -> P6};
bmat[n_] := Table[If[j <= i, Binomial[i, j]/Binomial[n, j], 0], {i, 0, n}, {j, 0, n}];
bern3[p_] := Module[{degs = Exponent[p, #] & /@ {X, Uv, Z}, A},
  A = CoefficientList[p, {X, Uv, Z}, degs + 1];
  Do[A = TensorContract[TensorProduct[bmat[degs[[k]]], A], {{2, k + 2}}]; A = Transpose[A, Join[Range[2, k], {1}, Range[k + 1, 3]]], {k, 3}];
  Flatten[A]];
Do[Module[{cc = load[n], p, b},
   p = poly[cc, "P", {X, Uv, Z}];
   ok[n <> ": JSON polynomial = expected one in the box map", Expand[p - ((n /. expected) /. box)] === 0];
   ok[n <> ": tree is a single leaf", cc["tree"] === {{"l"}}];
   b = bern3[p];
   ok[n <> ": all " <> ToString[Length[b]] <> " Bernstein coefficients on [0,1]^3 >= 0 (minimum " <> ToString[N[Min[b], 6]] <> ")", Min[b] >= 0]],
  {n, Keys[expected]}];
(* signs: -nU >= 0, -nV >= 0, P3 >= 0, P6 >= 0 give U <= 0, V <= 0, V^2 >= J^2, J^2 <= 4U^2 (positive factors), hence: *)
ok["deduction: U <= 0, V <= 0, V^2 >= J^2, J^2 <= 4U^2  =>  max over s in [-1,1] of U(1-s^2) + V s^2 + J s <= 0",
  Resolve[ForAll[{UU, VV, JJ, s}, (UU <= 0 && VV <= 0 && VV^2 >= JJ^2 && JJ^2 <= 4 UU^2 && -1 <= s <= 1) \[Implies] UU (1 - s^2) + VV s^2 + JJ s <= 0], Reals]];
ok["conversion factors positive on the box: 1 - x > 0, w > 0 for 0 < x < 1, 0 <= u <= 1/2",
  Resolve[ForAll[{x, u}, 0 < x < 1 && 0 <= u <= 1/2, 1 - x > 0 && 1 - u + u x > 0], Reals]];
Print["SUMMARY: ", pass, " PASS, ", fail, " FAIL"];
