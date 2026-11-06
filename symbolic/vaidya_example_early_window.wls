(* Paper III, constants of the Vaidya example (supplement S11): early window and the certificate for delta_c at the
   final mass m1 = 3/2 ([N] = high-precision numerical integrals, not certified; the certified bounds are in
   ../riemann_sums/). *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
(* units M=1, r0=4 (physical: m1=3/2, r0=6) *)
rho2[r_] := r^2 (r - 2)/(r + 2); Rp[r_] := Sqrt[r/(r + 2)];
ok["rho^2 = r^2 f/(E^2-f) with E^2=2, M=1", Simplify[rho2[r] - r^2 (1 - 2/r)/(2 - (1 - 2/r))] == 0];
(* physical m1 = 3/2: rho(5)^2/rho(7)^2 *)
rhoP2[r_] := r^2 (1 - 3/r)/(2 - (1 - 3/r));
ok["rho(5)^2/rho(7)^2 = 125/392 > 1/4 (m1=3/2)", rhoP2[5]/rhoP2[7] == 125/392 && 125/392 > 1/4];
ok["arccos(sqrt(125/392)) < pi/3", ArcCos[Sqrt[125/392]] < Pi/3];
ok["a^2 f = a(a-3) >= 10 on [5,7]", Resolve[ForAll[a, 5 <= a <= 7, a (a - 3) >= 10]]];
ok["R'/rho' = ... : R'/rho <= 1/(a sqrt f) (m1=3/2)", Simplify[(1/Sqrt[2 - 1 + 3/a])/Sqrt[rhoP2[a]] - 1/(a Sqrt[1 - 3/a]), a > 3] == 0];
ok["final margin 157/100 - 3h/10 - (39/25 + 4h/3) = 251/30000 at h=1/1000", (157/100 - 3/10000) - (39/25 + 4/3000) == 251/30000 && Pi > 157/50];
ok["|r'| < 2E^2: s+ and s- in [-2E^2, E^2/4]", Resolve[ForAll[{x, e}, e > 1 && 0 < x < e, -x^2 - e x > -2 e^2 && -x^2 + e x <= e^2/4]]];
(* delta_c: Delta(t) and l(t), 2R(4) [N] high precision *)
L[t_] := Sqrt[rho2[t]];
Del[t_?NumericQ] := 2 NIntegrate[L[t] Rp[r]/(Sqrt[rho2[r]] Sqrt[rho2[r] - rho2[t]]), {r, t, 4}, WorkingPrecision -> 40, MaxRecursion -> 40];
len[t_?NumericQ] := NIntegrate[Rp[r] Sqrt[rho2[r]]/Sqrt[rho2[r] - rho2[t]], {r, t, 4}, WorkingPrecision -> 40, MaxRecursion -> 40] 2;
Print["[N] Delta(29/10) = ", N[Del[29/10], 15], "   Delta(297/100) = ", N[Del[297/100], 15]];
Print["[N] l(297/100) = ", N[len[297/100], 15], "   2R(4) = ", N[2 NIntegrate[Rp[r], {r, 2, 4}, WorkingPrecision -> 40], 15]];
ok["[N] (real part) Delta(29/10) > 39/25 > Delta(297/100), l(297/100) > 2R(4)", Re[Del[29/10]] > 39/25 > Re[Del[297/100]] && Re[len[297/100]] > 2 NIntegrate[Rp[r], {r, 2, 4}, WorkingPrecision -> 40]];
