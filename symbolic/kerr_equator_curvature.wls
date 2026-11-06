#!/usr/bin/env wolframscript
(* Paper III, section on the Kerr equator, curvature: Ric_h(k,k) for h = G/N (constant and variable A), structure of
   U, V, J, asymptotics, independent numerical check of K_inf, grouping of the jets, conformal formula in 2+1. M = 1. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Kerr equator, curvature =="];
(* check of the curvature machinery with the spin a in the metric: Kerr is Ricci-flat *)
Module[{S = rr^2 + a^2 Cos[tt]^2, Dk = rr^2 - 2 rr + a^2, gK},
  gK = {{-(1 - 2 rr/S), 0, 0, -2 a rr Sin[tt]^2/S}, {0, S/Dk, 0, 0}, {0, 0, S, 0}, {-2 a rr Sin[tt]^2/S, 0, 0, (rr^2 + a^2 + 2 a^2 rr Sin[tt]^2/S) Sin[tt]^2}};
  chkN["curvature machinery: Kerr Ricci-flat (point a = 0.7, r = 3.3, theta = 1.1)", Max[Abs[N[ricci[gK, {t0, rr, tt, pp}] /. {a -> 7/10, rr -> 33/10, tt -> 11/10}]]], 10^-10]];
xs = {eta, r, ph};
f = 1 - 2/r; Dl = r^2 - 2 r + a^2;
hmet[uu_] := Module[{w = 1 - f uu, Gpp}, Gpp = r^2 + a^2 (1 + 2/r) + uu Dl/w;
   w {{-f, 0, -2 a/r}, {0, r^2/(w Dl), 0}, {-2 a/r, 0, Gpp}}];   (* h = G/N = (1 - f u) G *)
(* --- constant A --- *)
h = hmet[u];
Ric = Simplify[ricci[h, xs]];
(* ZAMO frame of h *)
om = -h[[1, 3]]/h[[3, 3]];
al2 = -(h[[1, 1]] - h[[1, 3]]^2/h[[3, 3]]);
e0 = {1, 0, om}/Sqrt[al2]; e1 = {0, 1, 0}/Sqrt[h[[2, 2]]]; e2 = {0, 0, 1}/Sqrt[h[[3, 3]]];
ass = 0 < a <= 1 && r > 1 + Sqrt[1 - a^2] && 0 <= u < 1;
chk["ZAMO: e0 unit timelike, e0 _|_ e2", {e0.h.e0 + 1, e0.h.e2, e1.h.e1 - 1, e2.h.e2 - 1} // Total[Abs[#]] &, ass];
R00 = e0.Ric.e0; R11 = e1.Ric.e1; R22 = e2.Ric.e2; R01 = e0.Ric.e1; R02 = e0.Ric.e2; R12 = e1.Ric.e2;
chk["Ric_h(k,k) = U c^2 + V s^2 + J s: Ric_01 = Ric_12 = 0 (constant A)", Abs[R01] + Abs[R12], ass];
U = R00 + R11; V = R00 + R22; J = 2 R02;
(* substitutions of the text: R = r_+, a^2 = R(2-R), x = R/r *)
subs = {a -> Sqrt[RR (2 - RR)], r -> RR/x};
assx = 1 < RR < 2 && 0 < x < 1 && 0 <= u < 1;
wt = RR (1 - u) + 2 u x; DD = (1 - x) (RR (1 + x) - 2 x); HH = (r^4 + a^2 r^2 + 2 a^2 r + 4 a^2 u) x^4/RR /. subs;
PU = Together[(U /. subs) RR^3 wt HH/(x^3 DD)];
PV = Together[-(V /. subs) RR^4 wt^2 HH/x^3];
PJ2 = Together[(J^2 /. subs) RR^7 wt^3 HH^2/(4 x^8 (2 - RR) DD)];
chk["structure: U R^3 w~ H/(x^3 D) is a polynomial in (x, R, u)", Boole[!PolynomialQ[Simplify[PU, assx], {x, RR, u}]], assx];
chk["structure: -V R^4 w~^2 H/x^3 is a polynomial in (x, R, u)", Boole[!PolynomialQ[Simplify[PV, assx], {x, RR, u}]], assx];
chk["structure: J^2 R^7 w~^3 H^2/(4 x^8 (2-R) D) is the square of a polynomial",
  Module[{pj = Simplify[PJ2, assx], sq}, sq = Simplify[Sqrt[Factor[pj]], assx]; Boole[!(PolynomialQ[sq, {x, RR, u}] || PolynomialQ[-sq, {x, RR, u}])]], assx];
chk["structure: w~, D, H positive on (x, R, u) in (0,1) x (1,2) x [0,1)",
  Reduce[Exists[{x, RR, u}, 0 < x < 1 && 1 < RR < 2 && 0 <= u < 1 && (wt <= 0 || DD <= 0 || HH <= 0)], Reals] === False, True];
(* asymptotics *)
zser[e_] := Normal[Series[e /. r -> 1/z, {z, 0, 3}]];
chk["asymptotics: U ~ -u z^3/(1-u)", Coefficient[zser[U], z, 3] + u/(1 - u) + 0 Coefficient[zser[U], z, 2], 0 < a < 1 && 0 <= u < 1];
chk["asymptotics: U, V without terms z^0, z^1, z^2", Total[Abs[Table[Coefficient[zser[X], z, j], {X, {U, V}}, {j, 0, 2}]], 2], 0 < a < 1 && 0 <= u < 1];
chk["asymptotics: V ~ z^3 (5u - 3)/(1-u)", Coefficient[zser[V], z, 3] - (5 u - 3)/(1 - u), 0 < a < 1 && 0 <= u < 1];
(* K_inf: independent numerical check (not certified) *)
fK = Compile[{{aa, _Real}, {rr, _Real}, {uu, _Real}, {ss, _Real}}, Evaluate[(U (1 - s^2) + V s^2 + J s) /. {a -> aa, r -> rr, u -> uu, s -> ss}]];
scan[umax_] := Module[{best = -Infinity, pt = {}, v, rng}, SeedRandom[7];
   Do[Module[{av = RandomReal[{0, 1}], uu = RandomReal[{0, umax}], ss = RandomReal[{-1, 1}], rr, rp},
      rp = 1 + Sqrt[1 - av^2]; rr = rp + 10^RandomReal[{-4, 2.5}];
      v = fK[av, rr, uu, ss] rr^3;
      If[v > best, best = v; pt = {av, rr, uu, ss}]], {200000}]; {best, pt}];
k35 = scan[3/5]; k23 = scan[2/3];
Print["   max r^3 Ric_h(k,k) over 2e5 random points, u <= 3/5: ", k35];
Print["   max r^3 Ric_h(k,k) over 2e5 random points, u <= 2/3: ", k23];
chk["prop:kerr-Kinf (numerical, not certified): Ric_h(k,k) <= 0 for u <= 3/5", k35[[1]] <= 10^-12, True];
chk["prop:kerr-Kinf (control): for u up to 2/3 positive values appear", k23[[1]] > 0, True];
(* --- variable A: jets --- *)
hv = hmet[uf[eta]];
Ricv = ricci[hv, xs];
omv = -hv[[1, 3]]/hv[[3, 3]]; al2v = -(hv[[1, 1]] - hv[[1, 3]]^2/hv[[3, 3]]);
k = {1, 0, omv}/Sqrt[al2v] + c {0, 1, 0}/Sqrt[hv[[2, 2]]] + s {0, 0, 1}/Sqrt[hv[[3, 3]]];
Rk = k.Ricv.k;
jet = {uf'[eta] -> u1, uf''[eta] -> u2, uf[eta] -> u};
Rkj = Rk /. jet;
Rstat = Rkj /. {u1 -> 0, u2 -> 0};
wv = 1 - f u; Cc = 2 a^2/(r^2 Dl wv);
rest = Rkj - Rstat + Cc c^2 (u2 + f u1^2/(2 wv));
chk["Ric_h jets: Ric_h = R_stat + D1 u' - C c^2[u'' + f u'^2/(2w)], no other terms in u'^2, u''",
  {Coefficient[Expand[rest], u2], Coefficient[Expand[rest], u1, 2]} // Total[Abs[#]] &, ass && c^2 + s^2 == 1];
D1 = Coefficient[Expand[Simplify[rest, ass]], u1, 1];
chk["Ric_h jets: D1 proportional to c (D1 = 0 for c = 0)", D1 /. c -> 0, ass];
chk["Ric_h jets: static part = U c^2 + V s^2 + J s", (Rstat - (U c^2 + V s^2 + J s)) /. s -> Sqrt[1 - c^2], ass && -1 < c < 1];
(* conformal formula in 2+1 for null k and sigma(r) *)
sig = sg[r];
hc = Exp[2 sig] h;
Ricc = ricci[hc, xs];
kn = e0 + c e1 + s e2;
Hr = Table[D[r, xs[[i]], xs[[j]]] - Sum[christoffel[h, xs][[l, i, j]] D[r, xs[[l]]], {l, 3}], {i, 3}, {j, 3}];
chk["conformal formula: Ric_{e^{2 sig} h}(k,k) = Ric_h(k,k) - sig' Hess_h r(k,k) + (sig'^2 - sig'')(k^r)^2 (null k)",
  (kn.Ricc.kn - kn.Ric.kn + sg'[r] kn.Hr.kn - (sg'[r]^2 - sg''[r]) kn[[2]]^2) /. s -> Sqrt[1 - c^2], ass && -1 < c < 1];
summary[];
