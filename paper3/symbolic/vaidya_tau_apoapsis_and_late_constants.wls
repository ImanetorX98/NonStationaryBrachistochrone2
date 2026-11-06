(* Paper III, proper-time branch of ingoing Vaidya: radial acceleration at a turning point on the constraint
   h = omega - q (the local apoapsis with rddot = -23/54 of the text), constants of the late-change construction
   (alpha = g/(4 pi + 8), r_s = 2M + M alpha^2/4), bound T0 < 91, and a polynomial competitor [N]. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
f = 1 - 2 m/r; w = e^2 - f; om = f/e; S = Sqrt[e^2 pr^2 + J^2/r^2]; h = -w pr + Sqrt[w] S;
as = {e > 1, m > 0, r > 2 m, J > 0, pr > 0};
pstar = Sqrt[w] J/(e r Sqrt[f]);
acc = (D[h, {pr, 2}] (-D[h - om, r])) /. pr -> pstar;
Jc = r (f - e q)/Sqrt[f w];
accq = Simplify[acc /. J -> Jc, {e > 1, m > 0, r > 2 m, 0 <= q < f/e}];
P = 8 m^3 + 2 m^2 (-6 + e^2 + 4 e q) r + 3 (e^2 - 2) m (e q - 1) r^2 - (e^2 - 1) (e q - 1) r^3;
ok["turning point: h_pp p_r' = (r-2m) P/(E^2 r^4 [r(1-Eq)-2m])", Simplify[accq - (r - 2 m) P/(e^2 r^4 (r (1 - e q) - 2 m)), {e > 1, m > 0, r > 2 m, 0 <= q < f/e}] == 0];
ok["turning point: P(0)/m^3 = a u^3 + 3(a+1)u^2 + 2(a+1)u", Simplify[(P /. q -> 0)/m^3 - (a u^3 + 3 (a + 1) u^2 + 2 (a + 1) u) /. {a -> e^2 - 1, u -> (r - 2 m)/m}] == 0];
ok["turning point: P(omega) = 2m(r-2m)[(E^2-1)r+2m]", Simplify[(P /. q -> f/e) - 2 m (r - 2 m) ((e^2 - 1) r + 2 m)] == 0];
ok["turning point: P affine in q", Simplify[D[P, {q, 2}]] == 0];
ok["turning point: d_q (r h_pp p_r') formula", Simplify[D[r accq, q] - 2 m (r - 2 m)^2 ((e^2 - 1) r + 2 m)/(e r^2 (r (1 - e q) - 2 m)^2), {e > 1, m > 0, r > 2 m, 0 <= q < f/e}] == 0];
val = r accq /. {e -> Sqrt[2], m -> 5/4, r -> 3, q -> ((1 - 5/(2*3))/Sqrt[2])/2};
ok["turning point: example r h_pp p_r' = 2/9 (m=5/4, r=3, q=omega/2)", Simplify[val] == 2/9];
ok["turning point: r'' = 2/9/3 - (3/2)/3 = -23/54", Simplify[val/3 - (3/2)/3] == -23/54];
(* late change: b(r)^2 = 1 - M^2/(r-M)^2 ; r_s = 2M + M al^2/4 ; b(r_s)^2 < al^2 *)
ok["late change: (r sqrt f/(r-M))^2 = 1 - M^2/(r-M)^2", Simplify[r^2 (1 - 2 M/r)/(r - M)^2 - (1 - M^2/(r - M)^2)] == 0];
bs2 = 1 - M^2/(M + M al^2/4)^2;
ok["late change: b(r_s)^2 <= al^2/2 + al^4/16 < al^2 for 0<al<1", Resolve[ForAll[al, 0 < al < 1, Simplify[bs2, M > 0] <= al^2/2 + al^4/16 < al^2]]];
(* late change: 4(pi al/2 + al) = g/2 with al = g/(4 pi + 8) *)
ok["late change: 4(pi al/2 + al) = g/2", Simplify[4 (Pi al/2 + al) - g/2 /. al -> g/(4 Pi + 8)] == 0];
(* T0 < 11 sqrt2 * 5.82 < 91 *)
ok["T0 bound: 11 sqrt2 582/100 < 91", 11 Sqrt[2] 582/100 < 91];
(* length of a polynomial competitor [N] *)
rr[x_] := 43/10 + x^2 + 7/10 x^4;
F[x_] := Sqrt[rr[x]/(rr[x] + 2) rr'[x]^2 + rr[x]^2 (rr[x] - 2)/(rr[x] + 2) Pi^2/16];
Print["polynomial competitor [N]: length = ", N[2 NIntegrate[F[x], {x, 0, 1}, WorkingPrecision -> 30], 15]];
