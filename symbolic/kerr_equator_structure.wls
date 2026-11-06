#!/usr/bin/env wolframscript
(* Paper III, section on the Kerr equator (sec:kerr-eq), algebraic part: extended metric G, structure, primitives, lift,
   L = 0, collar. M = 1. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Thakurta-Kerr on the equator, structure =="];
rp = 1 + Sqrt[1 - a^2];
ass = 0 < a <= 1 && r > rp && 0 <= u < 1;
f = 1 - 2/r; Dl = r^2 - 2 r + a^2; Nn = 1/(1 - f u);
Gpp = r^2 + a^2 (1 + 2/r) + u Dl/(1 - f u);
G = {{-f, 0, -2 a/r}, {0, Nn r^2/Dl, 0}, {-2 a/r, 0, Gpp}};   (* coordinates (eta, r, phi) *)
(* G = f g_opt for r > 2 *)
be = -2 a/(r f); h0rr = r^2/(f Dl); h0pp = Dl/f^2;
gopt = {{-1, 0, be}, {0, Nn h0rr, 0}, {be, 0, -be^2 + Nn h0pp}};
chk["eq:kerr-G: G = f g_opt (r > 2)", Total[Flatten[Abs[Simplify[G - f gopt, ass && r > 2]]]], ass && r > 2];
chk["eq:kerr-G: identity f(r^2 + a^2 + 2a^2/r) + 4a^2/r^2 = Delta", f (r^2 + a^2 + 2 a^2/r) + 4 a^2/r^2 - Dl];
chk["lem:kerr-G: 1 - f u > 0 for u in [0,1), r > r_+", Reduce[Exists[{av, rr, uu}, 0 < av <= 1 && rr > 1 + Sqrt[1 - av^2] && 0 <= uu < 1 && 1 - (1 - 2/rr) uu <= 0], Reals] === False, True];
chk["lem:kerr-G: det G = -N^2 r^2", Det[G] + Nn^2 r^2, ass];
chk["lem:kerr-G: det of the (eta, phi) block = -N Delta", Det[G[[{1, 3}, {1, 3}]]] + Nn Dl, ass];
Gi = Inverse[G];
chk["lem:kerr-G: G^{rr} = Delta/(N r^2)", Gi[[2, 2]] - Dl/(Nn r^2), ass];
chk["lem:kerr-G: G^{eta eta} = -(r^4 + a^2 r^2 + 2a^2 r + 4a^2 u)/(r^2 Delta)", Gi[[1, 1]] + (r^4 + a^2 r^2 + 2 a^2 r + 4 a^2 u)/(r^2 Dl), ass];
chk["lem:kerr-G: G_phiphi >= r^2", Reduce[Exists[{av, rr, uu}, 0 < av <= 1 && rr > 1 + Sqrt[1 - av^2] && 0 <= uu < 1 &&
     (rr^2 + av^2 (1 + 2/rr) + uu (rr^2 - 2 rr + av^2)/(1 - (1 - 2/rr) uu)) < rr^2], Reals] === False, True];
chk["lem:kerr-G: completing the square (tangential form)",
  ({de, dph}.G[[{1, 3}, {1, 3}]].{de, dph}) - (-(Nn Dl/Gpp) de^2 + Gpp (dph - 2 a/(r Gpp) de)^2), ass];
Om[s_] := (2 a/r + s Sqrt[Nn Dl])/Gpp;
chk["lem:kerr-G: X_pm = d_eta + Omega_pm d_phi null", Table[{1, 0, Om[s]}.G.{1, 0, Om[s]}, {s, {1, -1}}] // Total[Abs[#]] &, ass];
chk["lem:kerr-G: Delta <= r^2 (a^2 <= 2r)", Reduce[Exists[{av, rr}, 0 < av <= 1 && rr > 1 + Sqrt[1 - av^2] && rr^2 - 2 rr + av^2 > rr^2], Reals] === False, True];
chk["lem:kerr-G: |dr/deta| <= Delta/(r sqrt(G_phiphi)) from the completed square (null vector with optimal fixed dphi)",
  Module[{rd}, rd = Sqrt[(Nn Dl/Gpp)/(Nn r^2/Dl)]; rd - Dl/(r Sqrt[Gpp])], ass];
chk["lem:kerr-G: |Omega_pm| <= 2a/r^3 + sqrt(N Delta)/r^2 <= 2 + sqrt(N_max)",
  Reduce[Exists[{av, rr, uu}, 0 < av <= 1 && rr > 1 + Sqrt[1 - av^2] && 0 <= uu < 1 &&
     2 av/rr^3 + Sqrt[(rr^2 - 2 rr + av^2)/(1 - (1 - 2/rr) uu)]/rr^2 > 2 + Sqrt[1/(1 - uu)]], Reals] === False, True];
(* primitive F of r^2/Delta *)
rm = 1 - Sqrt[1 - a^2];
Fp = r + rp^2/(rp - rm) Log[r - rp] - rm^2/(rp - rm) Log[r - rm];
chk["F' = r^2/Delta (|a| < 1)", D[Fp, r] - r^2/Dl, 0 < a < 1 && r > rp];
chk["F' = r^2/Delta (|a| = 1)", D[r + 2 Log[r - 1] - 1/(r - 1), r] - (r^2/(r^2 - 2 r + 1)), r > 1];
chk["d_u G_phiphi = Delta/(1 - f u)^2", D[Gpp, u] - Dl/(1 - f u)^2, ass];
chk["omega_* = 2a/(r S_*) is the midpoint of (Omega_-, Omega_+)", (Om[1] + Om[-1])/2 - 2 a/(r Gpp), ass];
chk["for r > 2, 0 in (Omega_-, Omega_+) (G_eta eta < 0)", Simplify[Om[1] Om[-1] < 0, ass && r > 2], True];
(* null lift *)
chk["prop:kerr-class: eta' = (b + sqrt(b^2 + H c))/H solves G(eta' T + V, eta' T + V) = 0",
  Module[{e = (bb + Sqrt[bb^2 + HT cc])/HT}, -HT e^2 + 2 bb e + cc], HT > 0 && cc >= 0 && bb \[Element] Reals];
(* lem:kerr-Z *)
Module[{kph, kr2, sol},
  kph = (2 a/r)/Gpp;                  (* L = G(k, d_phi) = 0 with k^eta = 1: -2a/r + G_phiphi k^phi = 0 *)
  kr2 = Solve[{1, kr, kph}.G.{1, kr, kph} == 0, kr];
  chk["lem:kerr-Z: with L = 0, (dr/deta)^2 = Delta^2/(r^2 G_phiphi)", (kr^2 /. kr2[[1]]) - Dl^2/(r^2 Gpp), ass]];
(* bands and thresholds *)
chk["bands: a(R) = sqrt(R(2-R)); B1 starts at sqrt(264)/25 (R = 44/25)", Sqrt[44/25 (2 - 44/25)] - Sqrt[264]/25];
chk["bands: B4 ends at sqrt(2451)/50 (R = 57/50)", Sqrt[57/50 (2 - 57/50)] - Sqrt[2451]/50];
chk["bands: approximate values of a at the edges (0.650, 0.800, 0.903, 0.951, 0.9902)",
  Max[Abs[N[Sqrt[# (2 - #)] & /@ {44/25, 8/5, 143/100, 131/100, 57/50}] - {0.650, 0.800, 0.903, 0.951, 0.9902}]] < 0.0006, True];
chk["sqrt(2451)/50 = 0.990152...", Abs[N[Sqrt[2451]/50] - 0.990152] < 10^-6, True];
chk["rmk: -f/(2w) <= 1/2 on r > r_+", Reduce[Exists[{av, rr, uu}, 0 < av <= 1 && rr > 1 + Sqrt[1 - av^2] && 0 <= uu < 1 &&
     -(1 - 2/rr)/(2 (1 - (1 - 2/rr) uu)) > 1/2], Reals] === False, True];
(* radial acceleration at the minimum; boundary *)
Xh = {1, 0, Om[1]} /. r -> rp;
chk["Omega_+ = Omega_- at r_+ (Delta = 0)", (Om[1] - Om[-1]) /. r -> rp, 0 < a < 1];
Module[{dG = D[G, r], X = {1, 0, Om[1]}, XH, val},
  XH = {1, 0, 2 a/(rp (Gpp /. r -> rp))};
  val = (XH.(dG /. r -> rp).XH);
  chk["boundary: d_r G(X_H, X_H) = -(r_+ - 1)/(2(1 - f(r_+) u))", val + (rp - 1)/(2 (1 - (1 - 2/rp) u)), 0 < a < 1 && 0 <= u < 1]];
(* collar: a^*(q) *)
astar[q_] := Sqrt[4 q/(Sqrt[(4 q + 1)^2 + 3] + 4 q + 1)];
chk["collar: a^*(q)^2 is the positive root of 3a^4 u^2 + 8a^2 u + 32 a^2 - 16 with u = 1/q",
  (3 z^2 u^2 + 8 z u + 32 z - 16 /. {z -> astar[q]^2, u -> 1/q}), q > 1];
chk["collar: a^*(1) = 0.6234...", Abs[N[astar[1]] - 0.6234] < 0.00005, True];
chk["collar: a^*(q) -> 1/sqrt(2)", Limit[astar[q], q -> Infinity] - 1/Sqrt[2]];
chk["collar: a^* increasing in q > 1", Reduce[Exists[q, q > 1 && D[astar[qq], qq] /. qq -> q <= 0], Reals] === False ||
   Reduce[q > 1 && Together[D[astar[q]^2, q]] <= 0, q, Reals] === False, True];
chk["collar: 33/50 < a^*(5/2) ~ 0.672", 33/50 < astar[5/2] && Abs[N[astar[5/2]] - 0.672] < 0.0005, True];
(* numerator of the prograde tangential cost: sign of d_r (1/Omega_+) at r = 2 *)
Module[{costo = 1/Om[1], dc},
  dc = Simplify[D[costo, r] /. r -> 2, 0 < a <= 1 && 0 <= u < 1];
  Print["   d_r(1/Omega_+) at r = 2: ", Factor[dc]];
  chk["collar: sign of d_r(1/Omega_+)|_{r=2} = sign of 3a^4u^2 + 8a^2u + 32a^2 - 16 (positive ratio)",
    Reduce[Exists[{av, uu}, 0 < av <= 1 && 0 <= uu < 1 && (dc /. {a -> av, u -> uu})/(3 av^4 uu^2 + 8 av^2 uu + 32 av^2 - 16) <= 0], Reals] === False, True]];
chk["collar: the numerator is increasing in u", D[3 a^4 u^2 + 8 a^2 u + 32 a^2 - 16, u] > 0 // Simplify[#, 0 < a <= 1 && u >= 0] &, True];
summary[];
