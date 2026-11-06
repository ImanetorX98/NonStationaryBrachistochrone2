#!/usr/bin/env wolframscript
(* Paper III, sections on persistence and on ingoing Vaidya (arrival branch; non-decreasing and decreasing mass). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Ingoing Vaidya, arrival branch =="];
(* costs (eq:vaidya-costs) and confinement of cor:vaidya *)
fs = 1 - 2 m0/r;
wE = E0^2 - fs;
chk["cor:vaidya: -w + E sqrt(w) = f sqrt(w)/(E + sqrt(w))", -wE + E0 Sqrt[wE] - fs Sqrt[wE]/(E0 + Sqrt[wE]), E0 > 1 && r > 2 m0 && m0 > 0];
chk["cor:vaidya: f sqrt(w)/(E + sqrt(w)) < f/2 (sqrt(w) < E)",
  Reduce[Exists[{ee, ff}, ee > 1 && 0 < ff < 1 && ff Sqrt[ee^2 - ff]/(ee + Sqrt[ee^2 - ff]) >= ff/2], Reals] === False, True];
Lc = 2 (rf - rho) + 4 mm Log[(rf - 2 mm)/(rho - 2 mm)];
chk["cor:vaidya: L(rho, r_f) = int_rho^{r_f} dr/((1 - 2m/r)/2)", Lc - Integrate[2/(1 - 2 mm/x), {x, rho, rf}, Assumptions -> rf > rho > 2 mm > 0], rf > rho > 2 mm > 0];
chk["cor:vaidya: L -> oo for rho -> 2m (L >= 4m log((r_f-2m)/eps) with rho = 2m + eps)",
  Limit[Lc /. rho -> 2 mm + eps, eps -> 0, Direction -> "FromAbove", Assumptions -> rf > 2 mm > 0] === DirectedInfinity[1] || 
  Limit[Lc /. {rho -> 2 mm + eps, mm -> 1, rf -> 6}, eps -> 0, Direction -> "FromAbove"] === Infinity, True];
(* eq:vaidya-G *)
q0 = E0^2;
wv = q0 - 1 + 2 m[v]/r;
Gm = {{-wv + wv^2/q0, wv/q0, 0, 0}, {wv/q0, 1/q0, 0, 0}, {0, 0, r^2, 0}, {0, 0, 0, r^2 Sin[th]^2}};
xs = {v, r, th, ph};
ass = q0 > 1 && r > 2 m[v] && m[v] > 0;
chk["eq:vaidya-G: G(1, s, zeta) = -w + (s+w)^2/E^2 + r^2 |zeta|^2", {1, s, 0, zp}.Gm.{1, s, 0, zp} /. th -> Pi/2 // (# - (-wv + (s + wv)^2/q0 + r^2 zp^2)) &, ass];
chk["eq:vaidya-G: det of the (v, r) block = -w/E^2", Det[Gm[[1 ;; 2, 1 ;; 2]]] + wv/q0, ass];
chk["eq:vaidya-G: G(d_v, d_v) = -w f/E^2", Gm[[1, 1]] + wv (1 - 2 m[v]/r)/q0, ass];
chk["eq:vaidya-G: G^{-1}(dv, dv) = -1/w", Inverse[Gm][[1, 1]] + 1/wv, ass];
(* constant mass: v = t + r_*, G = (w f/E^2)[-dt^2 + N h0], N = E^2/w *)
Module[{wc = q0 - 1 + 2 m0/r, fc = 1 - 2 m0/r, Gt},
  (* dv = dt + dr/f: G in (t, r) *)
  Gt = Module[{J = {{1, 1/fc}, {0, 1}}, Gvr = {{-wc + wc^2/q0, wc/q0}, {wc/q0, 1/q0}}}, Transpose[J].Gvr.J];
  chk["eq:vaidya-G: constant mass => (t, r) block = (w f/E^2) diag(-1, N/f^2), N = E^2/w",
    Simplify[Gt - (wc fc/q0) DiagonalMatrix[{-1, (q0/wc)/fc^2}], r > 2 m0 && m0 > 0 && q0 > 1] // Flatten // Total[Abs[#]] &, r > 2 m0 && m0 > 0 && q0 > 1];
  chk["eq:vaidya-G: constant mass => angular part r^2 = (w f/E^2) N r^2/f", r^2 - (wc fc/q0) (q0/wc) r^2/fc, r > 2 m0 && m0 > 0]];
(* thm:vaidya-arrival: x' <= f/2 = x/(2r) *)
chk["thm:vaidya-arrival: f/2 = (r - 2m)/(2r)", (1 - 2 m[v]/r)/2 - (r - 2 m[v])/(2 r)];
(* thm:vaidya-sturm: Riemann tensor of G with arbitrary m(v) *)
Rm = riemann[Gm, xs];
ke = kv; k = {kv, kr, 0, kp};
nul = Solve[k.Gm.k == 0 /. th -> Pi/2, kr];
eperp = {0, 0, 1/r, 0};
eparS = {c1, c2, 0, c3};
Pv = 3 (2 q0 - 3) (q0 - 1) r^2 + (33 q0 - 38) m[v] r + 40 m[v]^2;
target = (r^2 kp)^2/(r^4 (q0 r + 2 m[v] - r)^2) (m[v] Pv/r + m'[v] r (3 (q0 - 1) r + 4 m[v]));
Do[
  kk = k /. nul[[i]];
  epar = eparS /. Solve[{eparS.Gm.kk == 0, eparS.Gm.eparS == 1, c1 == 0} /. th -> Pi/2, {c1, c2, c3}][[1]];
  Kperp = rq[Gm, Rm, eperp, kk, eperp, kk] /. th -> Pi/2;
  Kpar = rq[Gm, Rm, epar, kk, epar, kk] /. th -> Pi/2;
  chk["thm:vaidya-sturm: K_perp - K_par = L^2/(r^4 (E^2 r + 2m - r)^2)[m P/r + m' r (3(q-1)r + 4m)] (branch " <> ToString[i] <> ")",
    Kperp - Kpar - target, ass && kv > 0 && kp > 0],
  {i, Length[nul]}];
chkFalse["control: without the m' term (must fail)", Kperp - Kpar - (target /. m'[v] -> 0), ass && kv > 0 && kp > 0];
(* s_+, eq:vaidya-alpha, angular cone *)
chk["s_+ = -w + sqrt(q w) = w f/(sqrt(q w) + w)", (-wv + Sqrt[q0 wv]) - wv (1 - 2 m[v]/r)/(Sqrt[q0 wv] + wv), ass];
chk["s_+ <= f/2 (w <= q)", Reduce[Exists[{qq, ww}, qq > 1 && qq - 1 < ww < qq && ww (qq - ww)/(Sqrt[qq ww] + ww) > (qq - ww)/2], Reals] === False, True];
chk["null circle at constant r: r^2 phi'^2 = w f/q", (wv - wv^2/q0) - wv (1 - 2 m[v]/r)/q0, ass];
chk["eq:vaidya-alpha: the direction is null", (((-wv + Sqrt[q0 wv] Cos[al]) + wv)^2/q0 + r^2 (Sqrt[wv] Sin[al]/r)^2 - wv), ass];
(* lem:vaidya-direction through the Hamiltonian flow of G with m(v) *)
Gi = Simplify[Inverse[Gm] /. th -> Pi/2];
Hm = 1/2 {pv, pr, 0, pp}.Gi.{pv, pr, 0, pp};
dt[F_] := D[F, v] D[Hm, pv] + D[F, r] D[Hm, pr] - D[F, pv] D[Hm, v] - D[F, pr] D[Hm, r];
vd = D[Hm, pv]; rd = D[Hm, pr]; phd = D[Hm, pp];
cosal = (rd/vd + wv)/Sqrt[q0 wv];
sinal = (phd/vd) r/Sqrt[wv];
Bw = cosal (2 wv - q0 + 1) - 1/2 Sqrt[q0/wv] (3 wv - q0 + 1);
nulP = Solve[Hm == 0, pv];
Do[chk["lem:vaidya-direction: d(cos al)/dv = -sin^2(al) B(w, cos al)/r, arbitrary m(v) (root " <> ToString[i] <> ")",
    (dt[cosal]/vd + sinal^2 Bw/r) /. nulP[[i]], ass && pp > 0], {i, Length[nulP]}];
Bf[ww_, c_] := c (2 ww - q0 + 1) - 1/2 Sqrt[q0/ww] (3 ww - q0 + 1);
chk["lem:vaidya-direction: B(q, 1) = 1/2", Bf[q0, 1] - 1/2, q0 > 1];
chk["lem:vaidya-direction: B(q-1, c) = (q-1)c - sqrt(q(q-1))", Bf[q0 - 1, c] - ((q0 - 1) c - Sqrt[q0 (q0 - 1)]), q0 > 1];
(* decreasing mass: boundary, Hessian, H_q *)
Xb = {1, 2 m'[v], 0, 0};
chk["lem:vaidya-dec-boundary: G(X, X) = 4 m'(1 + m'/q) on r = 2m", (Xb.Gm.Xb /. th -> Pi/2 /. r -> 2 m[v]) - 4 m'[v] (1 + m'[v]/q0), q0 > 1 && m[v] > 0];
kb = {1, 2 m'[v], 0, kpb};
kpsol = Solve[(kb.Gm.kb /. th -> Pi/2 /. r -> 2 m[v]) == 0, kpb];
chk["lem:vaidya-dec-boundary: (2m)^2 |k^Omega|^2 = -4 m'(1 + m'/q) (k^v = 1)", ((2 m[v])^2 kpb^2 + 4 m'[v] (1 + m'[v]/q0)) /. kpsol[[1]], q0 > 1 && m[v] > 0];
Gc = christoffel[Gm, xs];
phiF = r - 2 m[v];
hessF = Table[D[phiF, xs[[a]], xs[[b]]] - Sum[Gc[[c, a, b]] D[phiF, xs[[c]]], {c, 4}], {a, 4}, {b, 4}];
Fv = -2 m''[v] + (m'[v] + (4 + 5/q0) m'[v]^2 + 4 (1 + q0) m'[v]^3/q0^2)/m[v];
chk["lem:vaidya-dec-boundary: Hess_G(r - 2m)(k, k) = (k^v)^2 F(v) for null k tangent to the boundary",
  ((kb.hessF.kb /. th -> Pi/2 /. r -> 2 m[v]) /. kpsol[[1]]) - Fv, q0 > 1 && m[v] > 0 && -q0 < m'[v] < 0];
chkFalse["control: with 5/q -> 4/q in F (must fail)", ((kb.hessF.kb /. th -> Pi/2 /. r -> 2 m[v]) /. kpsol[[1]]) - (Fv /. 5/q0 -> 4/q0), q0 > 1 && m[v] > 0];
Hq[z_] := z (3 (2 q0 - 3) (q0 - 1) + (33 q0 - 38) z + 40 z^2)/(3 (q0 - 1) + 4 z);
chk["lem:vaidya-dec-sturm: m P/r + m' r (3(q-1)r + 4m) = r (3(q-1)r + 4m)[m' + H_q(m/r)]",
  (m[v] Pv/r + m'[v] r (3 (q0 - 1) r + 4 m[v])) - r (3 (q0 - 1) r + 4 m[v]) (m'[v] + Hq[m[v]/r]), ass];
chk["lem:vaidya-dec-sturm: H_q strictly increasing on z > 0 for q >= 3/2",
  Reduce[Exists[{qq, z}, qq >= 3/2 && z > 0 && (D[Hq[zz] /. q0 -> qq, zz] /. zz -> z) <= 0], Reals] === False, True];
chk["prop:vaidya-dec-count: s_- = -w - sqrt(q w) <= -(q-1) - sqrt(q(q-1))",
  Reduce[Exists[{qq, ww}, qq > 1 && ww > qq - 1 && -ww - Sqrt[qq ww] > -(qq - 1) - Sqrt[qq (qq - 1)]], Reals] === False, True];
(* admissible profile *)
sv = 1/(1 + Exp[(v - 5)/5]); mp = 1 + sv/5;
dd = -D[mp, v];
chk["profile: -m' = d = s(1-s)/25", dd - sv (1 - sv)/25];
chk["profile: d <= 1/100", MaxValue[dd, v] - 1/100 <= 0, True];
chk["profile: |m''| <= d/5", Reduce[Exists[v, Abs[D[mp, {v, 2}]] > dd/5], v, Reals] === False, True];
chk["profile: d <= e^{1 - v/5}/25", Reduce[Exists[v, dd > Exp[1 - v/5]/25], v, Reals] === False, True];
(* in s in (0,1): ds/dv = -s(1-s)/5, m' = -s(1-s)/25, m'' = (1-2s) s (1-s)/125, m = 1 + s/5: all rational *)
mpS = -sS (1 - sS)/25; mppS = (1 - 2 sS) sS (1 - sS)/125; mS = 1 + sS/5; dS = sS (1 - sS)/25;
chk["profile: m'' as a function of s equals the derivative in v", (D[mp, {v, 2}] - (mppS /. sS -> sv)), True];
Do[
  FS = -2 mppS + (mpS + (4 + 5/qv) mpS^2 + 4 (1 + qv) mpS^3/qv^2)/mS;
  chk["profile: F + d/60 <= 0 for s in (0,1), q = " <> ToString[qv, InputForm], Reduce[0 < sS < 1 && FS + dS/60 > 0, sS, Reals] === False, True];
  chk["profile (stronger than needed): F <= -(13/30) d, sharp constant (sup as s -> 1), q = " <> ToString[qv, InputForm], Reduce[0 < sS < 1 && FS + 13 dS/30 > 0, sS, Reals] === False, True];
  chkFalse["control: F <= -(1/2) d does not hold (must fail), q = " <> ToString[qv, InputForm], Reduce[0 < sS < 1 && FS + dS/2 > 0, sS, Reals] === False, True];
  chk["profile: -q < m' < 0 for s in (0,1), q = " <> ToString[qv, InputForm], Reduce[0 < sS < 1 && (mpS >= 0 || mpS <= -qv), sS, Reals] === False, True],
  {qv, {3/2, 2}}];
chk["profile: H_2(z) = z + z^2(24 + 40 z)/(3 + 4z)", (Hq[z] /. q0 -> 2) - (z + z^2 (24 + 40 z)/(3 + 4 z))];
chk["profile: H_{3/2}(z) = 23/3 z^2 + 28/3 z^3/(3/2 + 4z)", (Hq[z] /. q0 -> 3/2) - (23/3 z^2 + 28/3 z^3/(3/2 + 4 z))];
chk["profile: R e^{-v/5} decreasing for v >= 0 (R = 6 + v/2)", Reduce[Exists[v, v >= 0 && D[(6 + v/2) Exp[-v/5], v] >= 0], v, Reals] === False, True];
chk["profile: R^2 e^{-v/5} decreasing for v >= 0", Reduce[Exists[v, v >= 0 && D[(6 + v/2)^2 Exp[-v/5], v] >= 0], v, Reals] === False, True];
chk["profile: 6e/25 < 1 and 36e/25 < 23/3", 6 E/25 < 1 && 36 E/25 < 23/3, True];
chk["profile: Sturm condition checked directly: d < H_q(m/R) for q in {3/2, 2}, v in [0, 400]",
  And @@ Table[NMinValue[{(Hq[mp/(6 + v/2)] /. q0 -> qv) - dd, 0 <= v <= 400}, v] > 0, {qv, {3/2, 2}}], True];
chk["H_q(z) >= (2q-3) z for q >= 3/2, z > 0", Reduce[Exists[{qq, z}, qq >= 3/2 && z > 0 && (Hq[z] /. q0 -> qq) < (2 qq - 3) z], Reals] === False, True];
summary[];
