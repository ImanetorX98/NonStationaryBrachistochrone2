(* Bridge between the curvature certificates of groups 4 and 5 (../certificates) and the geometry, from scratch in
   Wolfram Language, independently of the Python chain in ../certificate_generators.

   Group 4. Weighted representative g_p = e^{2 sigma} h, sigma = -(1/2) log(1 - R/r) (so e^{2 sigma} = r/(r - R)),
   h = G/N on the equator of conformal Kerr with a time-dependent factor u = u(eta), M = 1, R = r_+, a^2 = R(2 - R).
   Ric of g_p is computed directly from its Christoffel symbols (no conformal-change formula), contracted with
   k = E0/sqrt(lap) + c E1/sqrt(h_rr) + s E2/sqrt(h_pp) (ZAMO frame of h, null for h and for g_p), and split in jets:
   Ric(k,k) = R_sig + D1 u' + D2 u'' + D3 u'^2. The inequalities of the certificates are
     S: -R_sig >= 0,   I2: K2(-R_sig) - C c^2 >= 0 with C c^2 = -D2,   I1: K1^2 R_sig^2 - D1^2 >= 0.
   Each is written as A0 + A1 vh with vh = x a sqrt(Delta w), x = R/r, t = tan(alpha/2), and mapped to the box of its
   band. The check: P0_json/A0 = P1_json/A1 =: F (one common factor), g_json = vh^2, and F has constant sign on the box
   (its factors are listed and evaluated at the box corners and centre).
   Group 5. Null tangents X+- = d_eta + Om+- d_phi of G; d_r G(X+-, X+-) = alpha +- beta sqrt(N Delta), Q = alpha^2 -
   beta^2 N Delta. The check: -alpha = minus_alpha_json/den_alpha_json and Q = Q_json/den_Q_json exactly.
   Usage: wolframscript -file bridge_group4_group5.wls ../certificates *)
cdir = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], FileNameJoin[{DirectoryName[$InputFileName], "..", "certificates"}]];
pass = 0; fail = 0;
ok[label_, cond_] := (If[TrueQ[cond], pass++; Print["PASS ", label], fail++; Print["FAIL ", label]]);
(* progress bar, kept out of the log: written to $PROGRESS_FILE if set, otherwise to standard error *)
progFile = Environment["PROGRESS_FILE"]; progT0 = AbsoluteTime[]; progDone = 0; progTotal = 14;
progress[label_] := Module[{el, n, line},
   progDone++; el = Round[AbsoluteTime[] - progT0]; n = Floor[30 progDone/progTotal];
   line = "bridge [" <> StringRepeat["#", n] <> StringRepeat["-", 30 - n] <> "] " <> ToString[progDone] <> "/" <>
     ToString[progTotal] <> " (" <> label <> "), " <> ToString[el] <> " s elapsed, about " <>
     ToString[Round[el (progTotal - progDone)/progDone]] <> " s left";
   If[StringQ[progFile], Export[progFile, line <> "\n", "Text"], WriteString["stderr", "\r" <> line]]];
load[name_] := Import[FileNameJoin[{cdir, name <> ".json"}], "RawJSON"];
poly[cc_, k_, vars_] := Total[(ToExpression[#[[2]]] Times @@ (vars^#[[1]])) & /@ cc["polynomials"][k]];

(* ------------------------------------------------------------------ group 4: curvature from scratch *)
co = {eta, r, ph};
f = 1 - 2/r; Del = r^2 - 2 r + a^2;
w[x_] := 1 - f x;
hpp[x_] := r^2 + a^2 (1 + 2/r) + 4 a^2 x/r^2;
h[x_] := {{-f w[x], 0, -2 a w[x]/r}, {0, r^2/Del, 0}, {-2 a w[x]/r, 0, hpp[x]}};
gp = r/(r - R) h[uu[eta]];
gi = Together[Inverse[gp]];
Gam = Table[Together[1/2 Sum[gi[[i, l]] (D[gp[[l, j]], co[[k]]] + D[gp[[l, k]], co[[j]]] - D[gp[[j, k]], co[[l]]]), {l, 3}]], {i, 3}, {j, 3}, {k, 3}];
jet = {uu''[eta] -> u2, uu'[eta] -> u1, uu[eta] -> u0};
Ric = Table[Together[(Sum[D[Gam[[i, j, k]], co[[i]]] - D[Gam[[i, j, i]], co[[k]]] +
        Sum[Gam[[i, i, l]] Gam[[l, j, k]] - Gam[[i, k, l]] Gam[[l, j, i]], {l, 3}], {i, 3}]) /. jet], {j, 3}, {k, 3}];
Print["Ricci of the weighted representative computed"];
progress["Ricci computed"];
h0 = h[u0]; hee = h0[[1, 1]]; hep = h0[[1, 3]]; hrr = h0[[2, 2]]; hp0 = h0[[3, 3]];
lap = Together[-(hee - hep^2/hp0)];
ok["lap h_pp = Delta w (so 1/sqrt(lap h_pp) = 1/sqrt(Delta w))", Together[lap hp0 - Del w[u0]] === 0];
E0 = {1, 0, -hep/hp0}; E1 = {0, 1, 0}; E2 = {0, 0, 1};
kv = il E0 + cc ir E1 + ss ip E2;
ok["k is null for h (with c^2 + s^2 = 1)", Together[Expand[kv.h0.kv /. {il -> 1/Sqrt[lap], ir -> 1/Sqrt[hrr], ip -> 1/Sqrt[hp0]}] /. ss^2 -> 1 - cc^2] === 0];
RicK = Expand[kv.Ric.kv];
(* reduce monomials il^i ir^j ip^k: squares to 1/lap, 1/hrr, 1/h_pp; the odd part il ip to V = 1/sqrt(Delta w) *)
red[e_] := Module[{rules = CoefficientRules[Expand[e], {il, ir, ip}], out = 0},
   Do[Module[{ex = ru[[1]], cf = ru[[2]], odd},
      odd = Mod[ex, 2];
      If[!MemberQ[{{0, 0, 0}, {1, 0, 1}}, odd], Print["unexpected odd monomial ", ex]; fail++];
      out += cf lap^(-Quotient[ex[[1]], 2]) hrr^(-Quotient[ex[[2]], 2]) hp0^(-Quotient[ex[[3]], 2]) If[odd == {1, 0, 1}, VV, 1]],
     {ru, rules}];
   out];
Rsig = red[RicK /. {u1 -> 0, u2 -> 0}];
D1 = Coefficient[RicK, u1, 1] /. u2 -> 0;
D1sq = red[Expand[D1^2]];
D2 = red[Coefficient[RicK, u2, 1]];
D3 = red[Coefficient[RicK, u1, 2]];
ok["Ric(k,k) has no other jet monomials", Expand[RicK - ((RicK /. {u1 -> 0, u2 -> 0}) + D1 u1 + Coefficient[RicK, u1, 2] u1^2 + Coefficient[RicK, u2, 1] u2)] === 0];
ok["D2 is proportional to c^2 on c^2 + s^2 = 1 (C c^2 = -D2 >= 0 needs only the sign of C)",
   Together[(Expand[Numerator[Together[D2]]] /. ss^2 -> 1 - cc^2) /. cc -> 0] === 0 && FreeQ[Together[D2], VV]];
ok["D3 = D2 f/(2w): the u'^2 term groups with u'' as u'' + f u'^2/(2w)", Together[Expand[Numerator[Together[D3 - D2 f/(2 w[u0])]]] /. ss^2 -> 1 - cc^2] === 0];
(* to the variables of the certificates: r = R/x, a = sqrt(R(2-R)), u0 = u, c, s from t, V = x a vh/(vh^2) with vh^2 = gh *)
gh = Expand[R (2 - R) (1 - x) (R (1 + x) - 2 x) (R (1 - u) + 2 u x)];
ok["gh = x^2 a^2 Delta w (with r = R/x, a^2 = R(2-R))", Simplify[Together[(x^2 a^2 Del w[u0] /. {r -> R/x, u0 -> u}) /. a -> Sqrt[R (2 - R)]] - gh] === 0];
toBox[e_] := Module[{z},
   z = e /. {VV -> x a vh/gh, cc -> (1 - t^2)/(1 + t^2), ss -> 2 t/(1 + t^2)};
   z = z /. {r -> R/x, u0 -> u};
   z = Together[z /. a -> Sqrt[R (2 - R)]];
   z = Together[Expand[Numerator[z]] /. vh^2 -> gh /. vh^3 -> gh vh /. vh^4 -> gh^2]/Denominator[z];
   {Together[Coefficient[Numerator[z], vh, 0]/Denominator[z]], Together[Coefficient[Numerator[z], vh, 1]/Denominator[z]]}];
mR = toBox[-Rsig]; Cc2 = toBox[-D2]; DD = toBox[D1sq];
Print["group 4 quantities in (x, R, u, t, vh) ready"];
(* products of two A0 + A1 vh forms *)
mult[{p0_, p1_}, {q0_, q1_}] := {Together[p0 q0 + p1 q1 gh], Together[p0 q1 + p1 q0]};
ineqs[K1_, K2_] := <|"S" -> mR, "I2" -> (K2 mR - Cc2), "I1" -> (K1^2 mult[mR, mR] - DD)|>;
factorData = {};
bands = <|"B1" -> {8/5, 44/25, 3/2, 7/2}, "B2" -> {143/100, 8/5, 3, 9}, "B3" -> {131/100, 143/100, 5, 19}, "B4" -> {57/50, 131/100, 14, 110}|>;
Do[Module[{Rlo, Rhi, K1, K2, xb, box, ins},
   {Rlo, Rhi, K1, K2} = bands[b]; xb = 10 Rlo/(11 Rlo - 1);
   box = {x -> xb X, R -> Rlo + (Rhi - Rlo) Rv, u -> 2/5 U, t -> 2 T - 1};
   ins = ineqs[K1, K2];
   Do[Module[{cc4 = load["group4_" <> b <> "_" <> nm], P0, P1, g, A0, A1, F0, F1, Ff, pts, vals},
      P0 = poly[cc4, "P0", {X, Rv, U, T}]; P1 = poly[cc4, "P1", {X, Rv, U, T}]; g = poly[cc4, "g", {X, Rv, U, T}];
      {A0, A1} = ins[nm] /. box;
      F0 = Together[P0/A0]; F1 = Together[P1/A1];
      ok["group4_" <> b <> "_" <> nm <> ": P0/A0 = P1/A1 (one common factor F)", Together[F0 - F1] === 0];
      ok["group4_" <> b <> "_" <> nm <> ": g_json = vh^2 on the box", Expand[g - (gh /. box)] === 0];
      Ff = Factor[F0];
      Print["   F = ", Ff];
      AppendTo[factorData, <|"certificate" -> "group4_" <> b <> "_" <> nm,
         "factors" -> ({ToString[#[[1]], InputForm], #[[2]]} & /@ FactorList[Ff])|>];
      pts = Tuples[{0, 1/2, 1}, 4];
      vals = (Ff /. Thread[{X, Rv, U, T} -> #]) & /@ Select[pts, #[[1]] > 0 &];
      ok["group4_" <> b <> "_" <> nm <> ": F > 0 at the 54 sample points with X > 0 (sign check; the factors above are of known sign)",
         AllTrue[vals, # > 0 &]]; progress["group4_" <> b <> "_" <> nm]],
     {nm, {"S", "I2", "I1"}}]],
  {b, Keys[bands]}];

(* ------------------------------------------------------------------ group 5: radial confinement, from G *)
Clear[r, a, u];
f5 = 1 - 2/r; Del5 = r^2 - 2 r + a^2; W5 = r - u (r - 2); Hp5 = r^4 + a^2 r^2 + 2 a^2 r + 4 a^2 u; N5 = 1/(1 - f5 u);
Gee = -f5; Gep = -2 a/r; Gpp = Hp5/(r W5);
ok["G_pp = N h_pp", Together[Gpp - N5 (r^2 + a^2 (1 + 2/r) + 4 a^2 u/r^2)] === 0];
ok["G_ep^2 - G_ee G_pp = N Delta = r Delta/W", Together[Gep^2 - Gee Gpp - r Del5/W5] === 0 && Together[N5 Del5 - r Del5/W5] === 0];
Om = (-Gep + sg SQ)/Gpp;                      (* SQ = sqrt(N Delta), sg = +-1 *)
dG = D[Gee, r] + 2 D[Gep, r] Om + D[Gpp, r] Om^2;
dGx = Expand[Together[dG] /. SQ^2 -> r Del5/W5];
al = Together[Coefficient[Numerator[Together[dGx]], SQ, 0]/Denominator[Together[dGx]]] /. sg^2 -> 1;
be = Together[(Coefficient[Numerator[Together[dGx]], SQ, 1]/Denominator[Together[dGx]]) /. sg -> 1];
ok["d_r G(X+-, X+-) = alpha +- beta sqrt(N Delta)", Simplify[Together[(dGx /. sg^2 -> 1) - (al + sg be SQ)]] === 0];
Q5 = Together[al^2 - be^2 r Del5/W5];
box5 = {r -> R + Y (R - 1)/10, a -> Sqrt[R (2 - R)], R -> 57/50 + 31/50 Rv, u -> 2/5 U};
sub5[e_] := Together[(e /. r -> R + Y (R - 1)/10 /. a -> Sqrt[R (2 - R)]) /. {R -> 57/50 + 31/50 Rv, u -> 2/5 U}];
p5[nm_] := poly[load["group5_" <> nm], "P", {Y, Rv, U}];
ok["group5: -alpha = minus_alpha_json/den_alpha_json exactly", Together[sub5[-al] - p5["minus_alpha"]/p5["den_alpha"]] === 0];
ok["group5: Q = Q_json/den_Q_json exactly", Together[sub5[Q5] - p5["Q"]/p5["den_Q"]] === 0];
Print["group5: with den_alpha, minus_alpha, den_Q, Q > 0 (certified), -alpha > 0 and alpha^2 > beta^2 N Delta, hence alpha +- beta sqrt(N Delta) < 0."];
progress["group 5"];
Export[FileNameJoin[{DirectoryName[$InputFileName], "..", "results", "bridge_group4_factors.json"}], factorData, "RawJSON"];
Print["factor lists written to ../results/bridge_group4_factors.json (sign of the odd-power factors: check_factor_signs.py)"];
Print["SUMMARY: ", pass, " PASS, ", fail, " FAIL"];
