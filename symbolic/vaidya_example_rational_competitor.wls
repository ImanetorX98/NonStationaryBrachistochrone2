(* Paper III, constants of the Vaidya example (supplement S11; theorem on late accretion, D = 91, s = 11/5). M = 1, E^2 = 2. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
A[r_] := 2 ArcCosh[Sqrt[r/2]]; F[r_] := Sqrt[r (r + 2)] - 2 ArcSinh[Sqrt[r/2]];
rho[r_] := r Sqrt[(r - 2)/(r + 2)]; L = 26258028/10^7;
t0 = r /. FindRoot[rho[r] == L, {r, 4}, WorkingPrecision -> 30];
n = 4; bsN = Table[t0 + (6 - t0) (1 - i/n)^2, {i, 0, n}]; mids = (Most[bsN] + Rest[bsN])/2;
tsN = ((L/rho[#])/Sqrt[1 - (L/rho[#])^2] & /@ mids);
c = (Pi/2)/(2 Total[tsN (A[Most[bsN]] - A[Rest[bsN]])]); tsN = c tsN;
Print["numerical data: b = ", N[bsN, 8], "  tan = ", N[tsN, 8]];
(* rationalisation: radii rounded to 1/100 (end point 6 exact), tangents rounded UP to 1/100 *)
bs = Join[{6}, Rationalize[Round[Rest[bsN], 1/100], 0]];
ang0 = 2 Total[tsN (A[Most[bs]] - A[Rest[bs]])];
ts = Ceiling[tsN (Pi/2)/ang0, 1/100];
Print["rational data: b = ", bs, "  tan = ", ts];
ang = 2 Total[ts (A[Most[bs]] - A[Rest[bs]])];
len = 2 Total[Sqrt[1 + ts^2] (F[Most[bs]] - F[Rest[bs]])];
Print["total angle = ", N[ang, 15], "  (pi/2 = ", N[Pi/2, 15], ")"];
Print["length = ", N[len, 15]];
ok["angle >= pi/2 (sufficient budget)", N[ang - Pi/2, 40] > 0];
ok["length < 91/(11 sqrt2) = 5.8497", N[len - 91/(11 Sqrt[2]), 40] < 0];
ok["length < Lambda(11/5)", N[len - ((F[6] - F[11/5]) + (Sqrt[6 (6 + 3)] - 3 ArcSinh[Sqrt[6/3]] - (Sqrt[11/5 (11/5 + 3)] - 3 ArcSinh[Sqrt[(11/5)/3]]))), 40] < 0];
ok["radii above 2 and decreasing", And @@ Thread[Most[bs] > Rest[bs]] && Last[bs] > 2];
Print["exact expression of the length: ", len];
