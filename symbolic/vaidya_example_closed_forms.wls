(* Paper III, constants of the Vaidya example (supplement S11; theorem on late accretion, D = 91, s = 11/5). M = 1, E^2 = 2. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
F[r_, a_] := Sqrt[r (r + a)] - a ArcSinh[Sqrt[r/a]];
ok["d/dr [sqrt(r(r+a)) - a arcsinh(sqrt(r/a))] = sqrt(r/(r+a))", Simplify[D[F[r, a], r] - Sqrt[r/(r + a)], {r > 0, a > 0}] == 0];
Lam = (F[6, 2] - F[11/5, 2]) + (F[6, 3] - F[11/5, 3]);
Print["Lambda(11/5) in closed form = ", N[Lam, 15], "   (Riemann-sum enclosure: 5.94106, 5.94135)"];
Print["2R(4), M=1 = ", N[2 (F[4, 2] - F[2, 2]), 15], "   (Riemann-sum enclosure: 3.081722, 3.081750)"];
