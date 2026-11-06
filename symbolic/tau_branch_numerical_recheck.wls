#!/usr/bin/env wolframscript
(* Paper III, independent numerical recomputations (Wolfram Language, extended precision; NOT certified):
   table of delta_c (prop:certified), Maxwell angle (prop:maxwell-angle), numbers of the proper-time branch of Vaidya. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Numerical recomputations of the proper-time branch =="];
wp = 25;
Vq[q_, x_] := x^2 (x - 2)/((q - 1) x + 2);
kern[q_, r0_, t_, u_] := Module[{x = t + (r0 - t) u^2, W, Dd},
   W = (q - 1) x t (x + t - 2) + 2 (x^2 + x t + t^2 - 2 x - 2 t); Dd = ((q - 1) x + 2) ((q - 1) t + 2);
   {2 Sqrt[r0 - t]/(Sqrt[x (x - 2)] Sqrt[W/Dd]), 2 Sqrt[r0 - t] Vq[q, x]/(Sqrt[x (x - 2)] Sqrt[W/Dd])}];
Delta[q_?NumericQ, r0_?NumericQ, t_?NumericQ] := 2 Sqrt[Vq[q, t]] NIntegrate[kern[q, r0, t, u][[1]], {u, 0, 1}, WorkingPrecision -> wp, PrecisionGoal -> 14];
ell[q_?NumericQ, r0_?NumericQ, t_?NumericQ] := 2 NIntegrate[kern[q, r0, t, u][[2]], {u, 0, 1}, WorkingPrecision -> wp, PrecisionGoal -> 14];
R0[q_, r0_] := NIntegrate[Sqrt[x/((q - 1) x + 2)], {x, 2, r0}, WorkingPrecision -> wp, PrecisionGoal -> 16];
deltac[q_, r0_] := Module[{tm, ts, r2 = 2 R0[q, r0], grid, vals},
   grid = Table[2 + (r0 - 2) k/40, {k, 1, 39}];
   vals = Delta[q, r0, #] & /@ grid;
   tm = grid[[Ordering[vals, -1][[1]]]];
   ts = t /. FindRoot[ell[q, r0, t] - r2, {t, (tm + r0)/2, tm, r0 - 10^-6 (r0 - 2)}, WorkingPrecision -> wp, AccuracyGoal -> 14, PrecisionGoal -> 14];
   Delta[q, r0, ts]/Pi];
tab = {{6/5, {{0.478032, 0.478041}, {0.560637, 0.560647}, {0.628871, 0.628880}, {0.691007, 0.691015}}},
   {3/2, {{0.486257, 0.486267}, {0.575456, 0.575466}, {0.652406, 0.652416}, {0.725796, 0.725805}}},
   {2, {{0.493846, 0.493856}, {0.588406, 0.588417}, {0.671513, 0.671523}, {0.751130, 0.751141}}},
   {4, {{0.504191, 0.504201}, {0.605064, 0.605075}, {0.694412, 0.694424}, {0.778986, 0.778998}}},
   {10, {{0.509874, 0.509885}, {0.613780, 0.613792}, {0.705767, 0.705779}, {0.792024, 0.792036}}}};
Do[Module[{q = tab[[i, 1]], r0 = {4, 6, 10, 20}[[j]], enc = tab[[i, 2, j]], dc},
   dc = deltac[q, r0];
   chk["tab:deltac q = " <> ToString[q, InputForm] <> ", r0 = " <> ToString[r0] <> "M: delta_c/pi = " <> ToString[NumberForm[N[dc], 9]] <> " in " <> ToString[enc],
     enc[[1]] <= dc <= enc[[2]], True]], {i, 5}, {j, 4}];
(* Maxwell angle *)
Module[{q = 1001/1000, r0 = 1500, t1, t2, G, ds, r2},
  t1[d_?NumericQ] := t /. FindRoot[Delta[q, r0, t] - d, {t, 25, 21.33, 30}, WorkingPrecision -> 20];
  t2[d_?NumericQ] := t /. FindRoot[Delta[q, r0, t] - d, {t, 120, 105, 140}, WorkingPrecision -> 20];
  G[d_?NumericQ] := ell[q, r0, t2[d]] - ell[q, r0, t1[d]];
  ds = d /. FindRoot[G[d], {d, 0.67711 Pi}, WorkingPrecision -> 20, AccuracyGoal -> 12];
  r2 = 2 R0[q, r0];
  Print["   Maxwell: delta*/pi = ", N[ds/Pi, 12], ", periastra ", N[{t1[ds], t2[ds]}, 8], ", length ", N[ell[q, r0, t1[ds]], 10], " < 2R0 = ", N[r2, 10]];
  chk["prop:maxwell-angle: delta*/pi in [0.677106398, 0.677112765]", 0.677106398 <= ds/Pi <= 0.677112765, True];
  chk["prop:maxwell-angle: the two dips are shorter than the radial path", ell[q, r0, t1[ds]] < r2, True]];
(* Vaidya, m = 1 + tanh(v/5)/2, E^2 = 2, r0 = rf = 6 *)
Module[{E2 = 2, mf, wf, sp, sm, bT, b0s, Rinf, rho, vc, Wv, bsol, rsol, BV},
  mf[v_] := 1 + Tanh[v/5]/2; wf[v_, r_] := E2 - 1 + 2 mf[v]/r;
  sp[v_, r_] := -wf[v, r] + Sqrt[E2 wf[v, r]]; sm[v_, r_] := -wf[v, r] - Sqrt[E2 wf[v, r]];
  b0s = Table[bsol = NDSolveValue[{b'[v] == sp[v, b[v]], b[TT] == 6}, b, {v, 0, TT}, WorkingPrecision -> 20, PrecisionGoal -> 12, AccuracyGoal -> 12, MaxSteps -> 10^6]; {TT, bsol[0]}, {TT, {200, 400}}];
  Print["   b_T(0) for T = 200, 400: ", N[b0s, 10]];
  bsol = NDSolveValue[{b'[v] == sp[v, b[v]], b[400] == 6}, b, {v, 0, 400}, WorkingPrecision -> 20, PrecisionGoal -> 12, AccuracyGoal -> 12, MaxSteps -> 10^6];
  chk["Vaidya tau: b_inf(0) ~ 2.6014", Abs[bsol[0] - 2.6014] < 0.00006, True];
  rsol = NDSolveValue[{rr'[v] == sm[v, rr[v]], rr[0] == 6, WhenEvent[rr[v] - bsol[v] <= 0, "StopIntegration"]}, rr, {v, 0, 50}, WorkingPrecision -> 20, PrecisionGoal -> 12, AccuracyGoal -> 12];
  vc = rsol["Domain"][[1, 2]];
  (* the frontier b_T with T = 400 approximates b_inf only for v << T: integrate up to 100 and 150 and compare *)
  Wv = NIntegrate[Sqrt[E2] + Sqrt[wf[v, rsol[v]]], {v, 0, vc}, WorkingPrecision -> 15] + NIntegrate[Sqrt[E2] - Sqrt[wf[v, bsol[v]]], {v, vc, 100}, WorkingPrecision -> 15];
  Print["   tail of the frontier: int_100^150 = ", NIntegrate[Sqrt[E2] - Sqrt[wf[v, bsol[v]]], {v, 100, 150}, WorkingPrecision -> 15]];
  Rinf[r_] := NIntegrate[1/Sqrt[1 + 3/s], {s, 3, r}, WorkingPrecision -> 20];
  BV = Wv + Rinf[6];
  Print["   v_c = ", N[vc, 8], ", W(0,6) = ", N[Wv, 10], ", B_V = ", N[BV, 10], ", 2 R_inf(6) = ", N[2 Rinf[6], 10]];
  chk["Vaidya tau: B_V ~ 5.21295", Abs[BV - 5.21295] < 0.000006, True];
  chk["Vaidya tau: R_inf(r0) + R_inf(rf) ~ 4.62260", Abs[2 Rinf[6] - 4.62260] < 0.000006, True]];
summary[];
