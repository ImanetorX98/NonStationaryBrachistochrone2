#!/usr/bin/env wolframscript
(* Paper III, section on the proper-time branch with a constant factor A. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Proper-time branch, constant A =="];
ass = M > 0 && r > 2 M && A > 0 && E0 > A;
f = 1 - 2 M/r; Nn = 1/(1 - A^2 f/E0^2);
chk["lem:tau-density: (A^2 f sqrt(N)/E)^2 h0 = (A^4 N/E^2)(dr^2 + f r^2 dOmega^2) (radial and angular)",
  {(A^2 f Sqrt[Nn]/E0)^2 (1/f^2) - A^4 Nn/E0^2, (A^2 f Sqrt[Nn]/E0)^2 (r^2/f) - A^4 Nn f r^2/E0^2} // Total[Abs[#]] &, ass];
chk["lem:tau-density: limits as r -> 2M (radial -> A^4/E^2, angular -> 0)",
  {Limit[A^4 Nn/E0^2, r -> 2 M] - A^4/E0^2, Limit[A^4 Nn f r^2/E0^2, r -> 2 M]} // Total[Abs[#]] &, ass];
(* dive: cost of the arc, conformal times *)
chk["thm:dive: cost of the arc = A^2 sqrt(N) sqrt(f) r pi/E", (A^2 f Sqrt[Nn]/E0) (r/Sqrt[f]) Pi - A^2 Sqrt[Nn] Sqrt[f] r Pi/E0, ass];
chk["conformal time of the arc ~ eps^(-1/2): sqrt(N) r/sqrt(f) * sqrt(eps) -> finite constant",
  Module[{l = Limit[(Sqrt[Nn] r/Sqrt[f] Pi /. r -> 2 M + eps) Sqrt[eps], eps -> 0, Direction -> "FromAbove"]}, Boole[!(0 < l < Infinity)] /. True -> 0], ass];
(* proper acceleration of the circle near 2M (constant A): |a| ~ E^2 f^{-3/2}/(4 A^3 M) *)
Module[{g, ut, up, ar, mag},
  (* circular rail: A^2 f u^t = E, A^2(-f ut^2 + r^2 up^2) = -1; radial acceleration a^r = Gamma^r_tt ut^2 + Gamma^r_pp up^2 *)
  ut = E0/(A^2 f); up = Sqrt[(A^2 f ut^2 - 1)/A^2]/r;
  ar = (M f/r^2) ut^2 - r f up^2;            (* Christoffel symbols of Schwarzschild, invariant for constant A *)
  mag = A ar/Sqrt[f];                         (* |a| = sqrt(g_rr) |a^r|, g_rr = A^2/f *)
  chk["acceleration of the circle: |a| f^{3/2} -> E^2/(4 A^3 M) as r -> 2M",
    Limit[mag (f)^(3/2), r -> 2 M, Direction -> "FromAbove"] - E0^2/(4 A^3 M), ass]];
(* completion *)
cc = A^2/E0;
chk["prop:completion: h_tau = dR^2 + rho^2 dOmega^2 with R' = c sqrt(N), rho = c r sqrt(N f)",
  {(cc Sqrt[Nn])^2 - A^4 Nn/E0^2, (cc r Sqrt[Nn f])^2 - A^4 Nn f r^2/E0^2} // Total[Abs[#]] &, ass];
chk["prop:completion: rho^2 ~ 2 M c R near 2M (R ~ c (r - 2M))",
  Limit[((cc r Sqrt[Nn f])^2)/(2 M cc (cc (r - 2 M))), r -> 2 M, Direction -> "FromAbove"] - 1, ass];
qq = E0^2/A^2;
drho = D[r Sqrt[Nn f], r]/Sqrt[Nn];
chk["thm:antipodal: d rho/dR = (r - M)/(r sqrt f) + M sqrt(f)/(r (q - f))", drho - ((r - M)/(r Sqrt[f]) + M Sqrt[f]/(r (qq - f))), ass];
chk["thm:antipodal: (r - M)^2 - r(r - 2M) = M^2", (r - M)^2 - r (r - 2 M) - M^2];
(* Paper II *)
bb = qq - 1; Vp = r^2 (r - 2 M)/(bb r + 2 M);
chk["rho^2 = A^2 V (potential of Paper II)", (cc r Sqrt[Nn f])^2 - A^2 Vp, ass];
chk["c sqrt(N)/rho = 1/sqrt(r(r-2M))", cc Sqrt[Nn]/(cc r Sqrt[Nn f]) - 1/Sqrt[r (r - 2 M)], ass];
Kf = 1/Sqrt[r (r - 2 M)]; Wf = Kf/D[Vp, r];
chk["w(A sqrt V) = 2 V W, W = K/V'", (cc Sqrt[Nn])/D[A Sqrt[Vp], r] - 2 Vp Wf, ass];
chk["half angle: L w/(rho sqrt(rho^2 - L^2)) d rho = sqrt(V_t) dr/(sqrt(r(r-2M)) sqrt(V - V_t)), L = A sqrt(V_t)",
  ((A Sqrt[Vt]) (cc Sqrt[Nn]) /(cc r Sqrt[Nn f] Sqrt[(cc r Sqrt[Nn f])^2 - A^2 Vt])) - Sqrt[Vt]/(Sqrt[r (r - 2 M)] Sqrt[Vp - Vt]), ass && 0 < Vt < Vp];
(* regular kernel (M = A = 1) *)
Vq[x_] := x^2 (x - 2)/((qs - 1) x + 2);
chk["prop:certified: rho^2(r) - rho^2(t) = (r - t) W/D", Vq[r] - Vq[t] - (r - t) ((qs - 1) r t (r + t - 2) + 2 (r^2 + r t + t^2 - 2 r - 2 t))/(((qs - 1) r + 2) ((qs - 1) t + 2))];
(* lem:dips: comparison integral *)
chk["lem:dips: 2 int_L^rho0 L/(rho sqrt(rho^2 - L^2)) = 2 arccos(L/rho0)", 2 Integrate[L/(x Sqrt[x^2 - L^2]), {x, L, r0}, Assumptions -> r0 > L > 0] - 2 ArcCos[L/r0], r0 > L > 0];
(* rmk spin *)
Module[{rs = MM + Sqrt[MM^2 - aa^2 Cos[th]^2], Sg, Dk},
  Sg = (rs + ep)^2 + aa^2 Cos[th]^2; Dk = (rs + ep)^2 - 2 MM (rs + ep) + aa^2;
  chk["rmk spin: Sigma - 2Mr = eps(2 r_s - 2M + eps) on r = r_s + eps", Sg - 2 MM (rs + ep) - ep (2 rs - 2 MM + ep), MM > aa > 0];
  chk["rmk spin: Delta >= a^2 sin^2 on r >= r_s", Reduce[Exists[{m1, a1, t1, e1}, m1 > a1 > 0 && 0 < t1 < Pi && e1 >= 0 &&
      (Dk /. {MM -> m1, aa -> a1, th -> t1, ep -> e1}) < a1^2 Sin[t1]^2], Reals] === False ||
    Simplify[Dk - aa^2 Sin[th]^2 - ((rs + ep - MM)^2 - (MM^2 - aa^2 Cos[th]^2)) - 0, MM > aa > 0] === 0, True];
  chk["rmk spin: r_s'^2/Delta <= a^2/(M^2 - a^2)  (with Delta >= a^2 sin^2)",
    Reduce[Exists[{m1, a1, t1}, m1 > a1 > 0 && 0 < t1 < Pi &&
      (D[rs, th]^2/(aa^2 Sin[th]^2) /. {MM -> m1, aa -> a1, th -> t1}) > a1^2/(m1^2 - a1^2)], Reals] === False, True]];
summary[];
