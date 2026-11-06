(* Paper III, arrival branch of ingoing Vaidya for 1 < q = E^2 < 3/2: comparison of the modes and absence of return
   (bound on the accretion rate below 3/2). *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
B[w_, c_] := c (2 w - q + 1) - 1/2 Sqrt[q/w] (3 w - q + 1);
(* 1. turning point: rdot = -w + sqrt(q w) c = 0  <=>  c = sqrt(w/q); rddot = -w' + (sqrt(q)/(2 sqrt w)) w' c - sqrt(q w) sin(al) al' *)
(* with w' = 2 mp/r - 2 m rdot/r^2 and al' = sin(al) B/r ; at rdot = 0: *)
c0 = Sqrt[w/q]; s2 = 1 - w/q;
rdd = -(2 mp/r) + Sqrt[q]/(2 Sqrt[w]) (2 mp/r) c0 - Sqrt[q w] s2 B[w, c0]/r;
Qw = 2 w^2 - w (5 q/2 - 1) + q (q - 1)/2;
ok["sqrt(q w) B(w, sqrt(w/q)) = Q(w)", Simplify[Sqrt[q w] B[w, c0] - Qw, {q > 1, w > 0}] == 0];
ok["r'' = -(1/r)[m' + (1 - w/q) Q(w)] at the turning points", Simplify[rdd + (mp + s2 Qw)/r, {q > 1, w > 0, r > 0}] == 0];
ok["Q(q - 1/3) = -1/9 (r = 3m)", Simplify[Qw /. w -> q - 1/3] == -1/9];
ok["Q(q - 1) = -(q - 1) (r -> infinity)", Simplify[(Qw /. w -> q - 1) + (q - 1)] == 0];
(* 2. threshold on the accretion: mu*(q) = inf_{w in (q-1, q-1/3]} (1 - w/q)(-Q(w)) *)
g = Simplify[(1 - w/q) (-Qw)];
mu = Minimize[{g /. q -> 6/5, 1/5 < w <= 13/15}, w];
Print["q = 6/5: mu* = ", mu, "  ~ ", N[mu[[1]]]];
Do[Print["   q = ", qq, " : mu* = ", N[NMinimize[{g /. q -> qq, qq - 1 < w <= qq - 1/3}, w][[1]], 8], "   1/(27q) = ", N[1/(27 qq), 8], "   (q-1)/q = ", N[(qq - 1)/qq, 8]], {qq, {101/100, 11/10, 6/5, 13/10, 7/5, 149/100}}];
ok["g >= mu_*(q) = Min[1/(27 q), (q - 1)/q] on r >= 3m, for 1 < q < 3/2 (Resolve; the minimum is (q - 1)/q for q < 28/27 and 1/(27 q) otherwise)", Resolve[ForAll[{q, w}, 1 < q < 3/2 && q - 1 < w <= q - 1/3, g >= Min[1/(27 q), (q - 1)/q]]]];
(* 3. comparison of the modes: same polynomial as for Thakurta-Schwarzschild with M -> m *)
P[r_, m_] := 3 (2 q - 3) (q - 1) r^2 + (33 q - 38) m r + 40 m^2;
ok["P(r; m) = m^2 P_q(r/m; 1) (scaling)", Simplify[P[r, m] - m^2 P[r/m, 1]] == 0];
