(* Paper III, constants of the Vaidya example (supplement S11; theorem on late accretion, D = 91, s = 11/5). M = 1, E^2 = 2. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
Fa[r_, a_] := Sqrt[r (r + a)] - a ArcSinh[Sqrt[r/a]];
R2R6 = 2 (Fa[6, 2] - Fa[2, 2]); Lam = Fa[6, 2] - Fa[11/5, 2] + Fa[6, 3] - Fa[11/5, 3];
Print["2R(6) = ", N[R2R6, 12], "   Lambda(11/5) = ", N[Lam, 12]];
ok["2R(6) > Lambda(11/5) > 5.94 > 5.839", N[R2R6 - Lam, 30] > 0 && N[Lam - 594/100, 30] > 0];
ok["11 sqrt2 * 5839/1000 < 91", 11 Sqrt[2] 5839/1000 < 91];
ok["5839/1000 < 91/(11 sqrt2)", 5839/1000 < 91/(11 Sqrt[2])];
ok["sextic: rho^2(r)-rho^2(t) = (r-t)((t+2)r^2+(t^2-4)r+2t^2-4t)/((r+2)(t+2))", Simplify[r^2 (r - 2)/(r + 2) - t^2 (t - 2)/(t + 2) - (r - t) ((t + 2) r^2 + (t^2 - 4) r + 2 t^2 - 4 t)/((r + 2) (t + 2))] == 0];
