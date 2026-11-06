#!/usr/bin/env wolframscript
(* Paper III, sections on the proper-time branch with a variable factor and on the proper-time branch of ingoing
   Vaidya (symbolic identities). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Proper-time branch with A(eta) =="];
ass = M > 0 && r > 2 M && A > 0 && E0 > A;
f = 1 - 2 M/r; Nn = 1/(1 - A^2 f/E0^2); om = A^2 f/E0;
chk["rail as a graph: |xdot|_h0 = N^{-1/2} (cone: d eta = sqrt(N)|dx|)", Sqrt[Nn] (1/Sqrt[Nn]) - 1, ass];
chk["lem:radial-red: the angular speed w = (sqrt f/r) sqrt(1/N - rdot^2/f^2) saturates the cone",
  (rd^2/f^2 + (r^2/f) ((Sqrt[f]/r) Sqrt[1/Nn - rd^2/f^2])^2) - 1/Nn, ass && rd^2 < f^2/Nn];
chk["prop:tau-rail: d_r omega = 2 M A^2/(E r^2) > 0", D[om, r] - 2 M A^2/(E0 r^2), ass];
(* PMP: H = q + v_s |p| - omega with p0 = -1, v_s = N^{-1/2} *)
Af = AA[eta]; omf = Af^2 f/E0; vs = Sqrt[1 - Af^2 f/E0^2]; Nf = 1/(1 - Af^2 f/E0^2);
pnorm = Sqrt[f^2 pr^2 + f L^2/r^2];
Hm = qc + vs pnorm - omf;
chk["prop:tau-pmp: |p| = sqrt(N)(omega - q) on H = 0", (Solve[Hm == 0, qc][[1]] /. Rule -> Subtract)[[1]] /. qc -> 0 // (pnorm - Sqrt[Nf] (omf - qc) /. (Solve[Hm == 0, qc][[1]])) &, ass && AA[eta] > 0 && E0 > AA[eta]];
chk["prop:tau-pmp: qdot = -d_eta H = (A A' f/E)[N(omega - q)/E + 2] on H = 0",
  ((-D[Hm, eta]) - (Af AA'[eta] f/E0) (Nf (omf - qc)/E0 + 2)) /. Solve[Hm == 0, qc][[1]], ass && AA[eta] > 0 && E0 > AA[eta]];
chk["prop:tau-pmp: cos psi = f p_r/|p|, sin psi = sqrt(f) L/(r |p|) (velocity = normalised d H/d p)",
  Module[{xr = D[Hm, pr], xp = D[Hm, L]}, {xr/(vs f) - f pr/pnorm, xp/(vs Sqrt[f]/r) - Sqrt[f] L/(r pnorm)}] // Total[Abs[#]] &, ass && AA[eta] > 0 && E0 > AA[eta]];
chk["prop:tau-pmp: Clairaut rho_eff sin psi = L, rho_eff = r sqrt(N/f)(omega - q)",
  (r Sqrt[Nf/f] (omf - qc)) (Sqrt[f] L/(r Sqrt[Nf] (omf - qc))) - L, ass && AA[eta] > 0 && E0 > AA[eta] && omf > qc];
chk["prop:tau-pmp: constant A => rho_eff = (A^2/E) r sqrt(N f)", (r Sqrt[Nn/f] om) - (A^2/E0) r Sqrt[Nn f], ass];
chk["prop:tau-pmp: cost of waiting: omega - q = |p|/sqrt(N)", (omf - qc) - Sqrt[Nf] (omf - qc)/Sqrt[Nf], ass];
chk["eq:Ra-barrier: f |dx|_h0 >= |dr| and density A^2 f sqrt(N)/E |dx| >= A^2 sqrt(N)/E |dr|", (A^2 f Sqrt[Nn]/E0) (1/f) - A^2 Sqrt[Nn]/E0, ass];
chk["eq:Ra-barrier: N increasing in A", Simplify[D[Nn, A] > 0, ass], True];
chk["thm:tau-P(3): -d_eps H = (A0 f (N0 + 1)/E) alpha a eps = 0 (q = 0, |p| = sqrt(N0) omega0)",
  Module[{Ae = A0 + ep al, H0, pn0},
    pn0 = Sqrt[1/(1 - A0^2 f/E0^2)] A0^2 f/E0;
    H0 = qc + Sqrt[1 - Ae^2 f/E0^2] pn0 - Ae^2 f/E0;
    (-D[H0, ep] /. ep -> 0) - (A0 f (1/(1 - A0^2 f/E0^2) + 1)/E0) al], M > 0 && r > 2 M && 0 < A0 < E0];
chk["thm:tau-P: N0 A0^2 f/E^2 = N0 - 1", (1/(1 - A0^2 f/E0^2)) A0^2 f/E0^2 - (1/(1 - A0^2 f/E0^2) - 1), M > 0 && r > 2 M && 0 < A0 < E0];
chk["thm:escape: |rdot| <= f (v = f/sqrt(N) = f sqrt(1 - A^2 f/E^2) <= f)",
  Reduce[Exists[{mm1, rr, aa1, ee}, mm1 > 0 && rr > 2 mm1 && aa1 > 0 && ee > aa1 && Sqrt[1 - aa1^2 (1 - 2 mm1/rr)/ee^2] > 1], Reals] === False, True];
chk["thm:escape: time to descend to 2M + eps >= r0 - r_eps + 2M log((r0-2M)/eps)",
  Integrate[1/(1 - 2 M/x), {x, 2 M + ep, r0}, Assumptions -> r0 > 2 M + ep && ep > 0 && M > 0] - (r0 - (2 M + ep) + 2 M Log[(r0 - 2 M)/ep]), r0 > 2 M + ep && ep > 0 && M > 0];
Print["== Proper-time branch of ingoing Vaidya =="];
assV = E0 > 1 && mm > 0 && r > 2 mm;
fv = 1 - 2 mm/r; wv = E0^2 - fv;
chk["R(r) = int ds/sqrt(E^2 - 1 + 2m0/s) equals c int sqrt(N) (A = 1, M = m0)",
  (1/E0) Sqrt[1/(1 - fv/E0^2)] - 1/Sqrt[E0^2 - 1 + 2 mm/r], assV];
(* Vaidya rail: g = -f dv^2 + 2 dv dr + r^2 dOmega^2 *)
Module[{vd, sol},
  (* E = f vdot - rdot, normalisation; s = rdot/vdot *)
  chk["lem:vaidya-rail: d tau/dv = (f - s)/E  (from E = f vdot - rdot)", ((fv - s)/E0) - 1/(E0/(fv - s)), assV && fv > s];
  chk["lem:vaidya-rail: r^2 |zeta|^2 = w - (s + w)^2/E^2 from the normalisation",
    (fv - 2 s - (fv - s)^2/E0^2) - (wv - (s + wv)^2/E0^2), assV]];
chk["lem:vaidya-rail: on K, f - s >= E(E - sqrt w) (s <= s_+ = -w + E sqrt w)", (fv - (-wv + E0 Sqrt[wv])) - E0 (E0 - Sqrt[wv]), assV];
chk["lem:vaidya-rail: d_r omega = 2m/(E r^2), omega = f/E", D[fv/E0, r] - 2 mm/(E0 r^2), assV];
chk["lem:vaidya-rail: radial cost per unit radius = 1/sqrt(w) at the extreme speeds",
  {(fv - (-wv + E0 Sqrt[wv]))/(E0 (-wv + E0 Sqrt[wv])) - 1/Sqrt[wv], (fv - (-wv - E0 Sqrt[wv]))/(E0 (wv + E0 Sqrt[wv])) - 1/Sqrt[wv]} // Total[Abs[#]] &, assV && fv > 0];
chkN["support function of K: h = -w p_r + sqrt(w) S, S = sqrt(E^2 p_r^2 + |p_Omega|^2/r^2) (numerical maximum)",
  Module[{sp, zp, sol, val},
    (* max of p_r s + p_Om zeta on (s + w)^2/E^2 + r^2 zeta^2 = w: parametrisation s = -w + E sqrt(w) cos t, zeta = sqrt(w) sin t/r *)
    val = MaxValue[{pr (-wv + E0 Sqrt[wv] Cos[t]) + po Sqrt[wv] Sin[t]/r, 0 <= t <= 2 Pi}, t] /. {pr -> 3/10, po -> 7/5, E0 -> Sqrt[2], mm -> 1, r -> 6};
    N[val - ((-wv pr + Sqrt[wv] Sqrt[E0^2 pr^2 + po^2/r^2]) /. {pr -> 3/10, po -> 7/5, E0 -> Sqrt[2], mm -> 1, r -> 6})]], 10^-12];
chk["thm:vaidya-P: launch: det d(p_r0, L)/d(psi, q_in) = r0 S0/(E b(psi))",
  Module[{w0 = wv /. r -> r0, b, S0, pr0, L0},
    b = Sqrt[w0] - (w0/E0) Cos[ps]; S0 = (om0 - qin)/b;
    pr0 = (S0/E0) Cos[ps]; L0 = S0 r0 Sin[ps];
    Det[{{D[pr0, ps], D[pr0, qin]}, {D[L0, ps], D[L0, qin]}}] - r0 S0/(E0 b)], E0 > 1 && mm > 0 && r0 > 2 mm];
chk["thm:vaidya-P: with h = omega, bracket -d_eps H/beta = p_r/r - (E^2 + w)/(E r w)",
  Module[{me = mm + ep be, wE, fE, Sx, hE, HE},
    wE = E0^2 - (1 - 2 me/r); fE = 1 - 2 me/r;
    Sx = Sqrt[E0^2 pr^2 + po^2/r^2];
    HE = qc + (-wE pr + Sqrt[wE] Sx) - fE/E0;
    ((-D[HE, ep] /. ep -> 0)/be /. po -> Sqrt[r^2 (((fv/E0 + wv pr)/Sqrt[wv])^2 - E0^2 pr^2)]) - (pr/r - (E0^2 + wv)/(E0 r wv))], assV && be != 0 && pr > 0];
chk["thm:vaidya-P: p_r >= 0 => p_r <= 1/E + 1/sqrt(w) (from omega = h >= (-w + E sqrt w) p_r)",
  ((fv/E0)/(-wv + E0 Sqrt[wv])) - (1/E0 + 1/Sqrt[wv]), assV && fv > 0];
chk["thm:vaidya-P: hence bracket <= (sqrt(w) - E)/(r w) < 0",
  ((1/E0 + 1/Sqrt[wv])/r - (E0^2 + wv)/(E0 r wv)) - (Sqrt[wv] - E0)/(r wv), assV];
(* escape price *)
chk["escape: c_- = E + sqrt(w) on the maximal infall", (fv - (-wv - E0 Sqrt[wv]))/E0 - (E0 + Sqrt[wv]), assV];
chk["escape: c_+ = E - sqrt(w) on the maximal rise", (fv - (-wv + E0 Sqrt[wv]))/E0 - (E0 - Sqrt[wv]), assV];
chk["escape: s_- <= -(E^2 - 1)", Reduce[Exists[{ee, ww}, ee > 1 && ee^2 - 1 < ww < ee^2 && -ww - ee Sqrt[ww] > -(ee^2 - 1)], Reals] === False, True];
chk["lem:vaidya-calib: d_r sqrt(w) = -m/(r^2 sqrt w)", D[Sqrt[wv], r] + mm/(r^2 Sqrt[wv]), assV];
chk["lem:vaidya-calib: d_r s_- = (2m/r^2)(1 + E/(2 sqrt w))", D[-wv - E0 Sqrt[wv], r] - (2 mm/r^2) (1 + E0/(2 Sqrt[wv])), assV];
chk["lem:vaidya-calib: contact: 2 sqrt(w) d_r v_c = xi/E with d_r v_c = xi/(s_+ - s_-)",
  2 Sqrt[wv] (xi/((-wv + E0 Sqrt[wv]) - (-wv - E0 Sqrt[wv]))) - xi/E0, assV];
chk["lem:vaidya-calib: final algebra d_r W - 1/E = (1/E) int (2m/rho^2) xi (pointwise integrands)",
  Module[{I1 = (2 mm/r^2) (1 + E0/(2 Sqrt[wv])) xi, I2 = mm xi/(r^2 Sqrt[wv])}, (1/E0) I1 - I2 - (1/E0) (2 mm/r^2) xi], assV];
chk["lem:vaidya-calib: W_v + s W_r + (f-s)/E = (s - s_-)(W_r - 1/E) given W_v + s_- W_r + (f - s_-)/E = 0",
  Module[{Wv = -sm Wr - (fv - sm)/E0}, (Wv + s Wr + (fv - s)/E0) - (s - sm) (Wr - 1/E0)], True];
summary[];
