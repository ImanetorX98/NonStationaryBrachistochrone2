#!/usr/bin/env wolframscript
(* Paper III, examples with a time-dependent factor (supplement S9): exact checks of the constants of the pulse at a
   small separation and of the drop at a large separation. M = 1, Ehat = Sqrt[2], r0 = r_f = 6. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
(* coefficients of the costate bound *)
ok["e^(3/100) <= 100/97", Exp[3/100] <= 100/97];
ok["(1/4)e^(3/100) <= 25/97 < 1/3", (1/4) Exp[3/100] <= 25/97 < 1/3];
ok["uniform corrected coefficient: (1/4)(100/97) = 25/97 < 1/3", (1/4)(100/97) == 25/97 && 25/97 < 1/3];
ok["(d/4)e^(3d/25) <= 25 d/97 for 0<=d<=1/4", Resolve[ForAll[d, 0 <= d <= 1/4, d/4 Exp[3 d/25] <= 25 d/97]]];
ok["4(e^(300h)-1) <= 12/97 for 0<=h<=1e-4", 4 (Exp[3/100] - 1) <= 12/97];
(* the angular speed near r = 6 and the chord *)
f[r_] := 1 - 2/r; hp[r_] := r^2/f[r];
w[r_, A_] := Sqrt[1 - A^2 f[r]/2]/Sqrt[hp[r]];
rb = 1801/300; ru = 1799/300;
ok["w(rbar,1/2) > 1302/10^4", w[rb, 1/2] > 1302/10^4];
ok["w(rbar,1) > 1110/10^4", w[rb, 1] > 1110/10^4];
ok["Jensen upper bound < 61/500", Sqrt[1 - 19 f[ru]/64]/Sqrt[hp[ru]] < 61/500];
ok["w decreasing in r on [ru,rb]", Resolve[ForAll[{r, A}, ru <= r <= rb && 1/2 <= A <= 1, D[w[r, A]^2, r] < 0]]];
ok["w concave in A^2", Resolve[ForAll[{r, x}, ru <= r <= rb && 1/4 <= x <= 1, D[Sqrt[1 - x f[r]/2]/Sqrt[hp[r]], {x, 2}] <= 0]]];
ok["int alpha^2 = 19/32", Integrate[(1 - Sin[Pi q]^2/2)^2, {q, 0, 1}] == 19/32];
ok["chord combination is an equality: (13/24)(1302/10^4) + (11/24)(1110/10^4) = 607/5000", (13/24)(1302/10000) + (11/24)(1110/10000) == 607/5000];
ok["upper switch enclosure: (1/50)/(607/5000 - 27/200000) = 4000/24253", (1/50)/(607/5000 - 27/200000) == 4000/24253];
(* the cubic drop profile *)
aa = 1 - (3 z^2 - 2 z^3)/2;
ok["drop is non-increasing", Resolve[ForAll[z, Implies[0 <= z <= 1, D[aa, z] <= 0]], Reals]];
ok["drop joins C1 to both plateaus", {aa /. z -> 0, aa /. z -> 1, D[aa, z] /. z -> 0, D[aa, z] /. z -> 1} == {1, 1/2, 0, 0}];
(* margins *)
ok["terminal derivative margin: (21/100)(46/11)(1100/147) = 46/7", (21/100)(46/11)(1100/147) == 46/7];
ok["strict triangular margin: 9/100 - (21/500)(21/50) > 7/100", 9/100 - (21/500)(21/50) > 7/100];
