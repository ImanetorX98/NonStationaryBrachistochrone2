(* Paper III, proper-time branch: general bounds used for the variable factor and for Vaidya (warp drho/dR > 1, angular
   speed, monotonicity of the radial density in A). *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
f = 1 - 2 M/r; nN = 1/(1 - A^2 f/e^2);
rho = A^2/e r Sqrt[nN f]; Rp = A^2/e Sqrt[nN];
q = Simplify[D[rho, r]/Rp, {r > 2 M, M > 0, 0 < A < e}];
Print["drho/dR = ", q];
ok["drho/dR > 1 for r > 2M, 0 < A < E (rationalised form)", Resolve[ForAll[{r, M, A, e}, r > 2 M && M > 0 && 0 < A < e, (r - M)^2 > r^2 f && (q > 1)]]];
ok["max sqrt(f)/r = 1/(3 sqrt3 M)", Simplify[Maximize[{Sqrt[1 - 2/x]/x, x > 2}, x][[1]] == 1/(3 Sqrt[3])]];
ok["R_Ainf' <= (Ainf/A0)^2 R_A0' (N increasing in A)", Simplify[D[nN, A] > 0, {r > 2 M, M > 0, 0 < A < e}]];
