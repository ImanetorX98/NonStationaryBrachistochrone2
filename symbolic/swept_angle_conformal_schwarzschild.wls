#!/usr/bin/env wolframscript
(* Paper III, sections on conformal Schwarzschild and on the swept-angle theorem. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Conformal Schwarzschild =="];
assR = M > 0 && r > 2 M;
f = 1 - 2 M/r;
chk["lem:h0-complete: rho' = 1/f with rho = r - 3M + 2M ln((r-2M)/M)", D[r - 3 M + 2 M Log[(r - 2 M)/M], r] - 1/f, assR];
chk["lem:h0-complete: rho(3M) = 0", (r - 3 M + 2 M Log[(r - 2 M)/M]) /. r -> 3 M, M > 0];
chk["rmk:freezing: Z = 1 at the radius r_F = 2M/(1 - E^2/A^2)", (A^2 f/E0^2 - 1) /. r -> 2 M/(1 - E0^2/A^2), A > E0 > 0 && M > 0];
chk["rmk:freezing: maximal radial speed f/sqrt(N) from the cone (h0_rr = 1/f^2)",
  Module[{v}, v = f/Sqrt[NN]; -1 + NN v^2/f^2], assR && NN > 0];
Print["== Swept angle =="];
qq = q[eta]; V = Sqrt[1 - f/qq];
chk["eq:tk-alpha: the direction (f V cos al, sqrt(f) V sin al / r) is null for g_inf",
  -V^2 + (f V Cos[al])^2/f^2 + (r^2/f) (Sqrt[f] V Sin[al]/r)^2, assR && q[eta] > 1];
(* lem:tk-sturm: direct computation of the Riemann tensor of g_inf with arbitrary q(eta) *)
xs = {eta, r, th, ph};
g = DiagonalMatrix[{-(1 - f/q[eta]), 1/f^2, r^2/f, r^2 Sin[th]^2/f}];
Rm = riemann[g, xs];
(* null vector in the equatorial plane: k = (ke, kr, 0, kp), ke fixed by the null condition *)
ke = Sqrt[(kr^2/f^2 + r^2 kp^2/f)/(1 - f/q[eta])];
k = {ke, kr, 0, kp};
eperp = {0, 0, Sqrt[f]/r, 0};
(* e_par: unit, orthogonal to k and to e_perp; chosen in the span of d_eta and d_r plus a phi component *)
epar0 = {c1, c2, 0, c3};
sol = Solve[{epar0.g.k == 0, epar0.g.epar0 == 1, c1 == 0}, {c1, c2, c3}];
epar = epar0 /. sol[[1]];
eqv = th -> Pi/2;
Kperp = rq[g, Rm, eperp, k, eperp, k] /. eqv;
Kpar = rq[g, Rm, epar, k, epar, k] /. eqv;
Pq = 3 (2 q[eta] - 3) (q[eta] - 1) r^2 + (33 q[eta] - 38) M r + 40 M^2;
target = kp^2 M Pq/(r (2 M + (q[eta] - 1) r)^2);
chk["lem:tk-sturm: K_perp - K_par = (dphi/dlam)^2 M P_q/(r(2M+(q-1)r)^2), any q(eta)",
  Kperp - Kpar - target, assR && q[eta] > 1 && kp > 0 && kr \[Element] Reals];
chk["lem:tk-sturm: no derivative of q in the result", FreeQ[Simplify[(Kperp - Kpar) /. eqv, assR && q[eta] > 1], Derivative], True];
chkFalse["control: P_q with the opposite sign of the 40M^2 term (must fail)",
  Kperp - Kpar - kp^2 M (Pq - 80 M^2)/(r (2 M + (q[eta] - 1) r)^2), assR && q[eta] > 1 && kp > 0];
chk["lem:tk-sturm: q >= 3/2 => coefficients of P_q non-negative", 
  Reduce[ForAll[qv, qv >= 3/2, 3 (2 qv - 3) (qv - 1) >= 0 && 33 qv - 38 > 0], Reals], True];
(* lem:tk-direction *)
Bq = (q[eta] - 1) r^2 + (6 - 3 q[eta]) M r - 8 M^2;
Cq = -Bq/(r^2 (2 M + (q[eta] - 1) r));
(* planar null geodesic of g_inf: Hamiltonian with p_eta = -Eps (not conserved if q depends on eta) *)
Hm = 1/2 (-pe^2/(1 - f/q[eta]) + f^2 pr^2 + f pp^2/r^2);
(* cos al = (dr/deta)/(f V); dr/deta = (dH/dpr)/(dH/dpe) *)
dot[F_] := D[F, r] D[Hm, pr] + D[F, eta] D[Hm, pe] - D[F, pr] D[Hm, r] - D[F, pe] D[Hm, eta] + D[F, pp] 0;
etad = D[Hm, pe];
cosal = (D[Hm, pr]/etad)/(f Sqrt[1 - f/q[eta]]);
sinal = (D[Hm, pp]/etad) r/(Sqrt[f] Sqrt[1 - f/q[eta]]);
nullo = pe -> -Sqrt[(1 - f/q[eta]) (f^2 pr^2 + f pp^2/r^2)];
chk["lem:tk-direction: d(cos al)/deta = -sin^2(al) V C (that is, d al/deta = V sin al C), any q(eta)",
  ((dot[cosal]/etad) + sinal^2 Sqrt[1 - f/q[eta]] Cq) /. nullo, assR && q[eta] > 1 && pp > 0 && pr \[Element] Reals];
chk["lem:tk-direction: B(q,2M) = -2 q M^2", (Bq /. r -> 2 M) + 2 q[eta] M^2];
chk["lem:tk-direction: B(q,3M) = M^2", (Bq /. r -> 3 M) - M^2];
chk["lem:tk-direction: B_q = r(r - 3M)", D[Bq /. q[eta] -> qs, qs] - r (r - 3 M)];
fr = D[f, r];
chk["lem:tk-direction: C = f_r/2 - f f_r/(2(q-f)) - f/r", Cq - (fr/2 - f fr/(2 (q[eta] - f)) - f/r), assR && q[eta] > 1];
chk["lem:tk-direction: C <= 3M/r^2 - 1/r for q > 1, r > 2M",
  Reduce[Exists[{rr, qs}, rr > 2 && qs > 1 && (-((qs - 1) rr^2 + (6 - 3 qs) rr - 8)/(rr^2 (2 + (qs - 1) rr))) > 3/rr^2 - 1/rr], Reals] === False, True];
chk["lem:tk-direction: 3M/r^2 - 1/r <= -1/(2r) for r >= 6M", Reduce[Exists[rr, rr >= 6 && 3/rr^2 - 1/rr > -1/(2 rr)], Reals] === False, True];
chk["lem:no-return: C <= -(r - 3M)/r^2 (same inequality)", (3 M/r^2 - 1/r) + (r - 3 M)/r^2, assR];
chk["lem:no-return: V^2 > 2M/r for q > 1", Reduce[Exists[{rr, qs}, rr > 2 && qs > 1 && 1 - (1 - 2/rr)/qs <= 2/rr], Reals] === False, True];
chk["thm:angle: q >= 3/2 => V^2 >= 1/3 and N <= 3", Reduce[Exists[{rr, qs}, rr > 2 && qs >= 3/2 && 1 - (1 - 2/rr)/qs < 1/3], Reals] === False, True];
chk["thm:angle: f V <= 1 (s_+ <= 1)", Reduce[Exists[{rr, qs}, rr > 2 && qs > 1 && (1 - 2/rr)^2 (1 - (1 - 2/rr)/qs) > 1], Reals] === False, True];
(* lem:radius *)
P[qs_, rr_] := 3 (2 qs - 3) (qs - 1) rr^2 + (33 qs - 38) M rr + 40 M^2;
Rq[qs_] := M (33 qs - 38 + Sqrt[(33 qs - 38)^2 + 480 (3 - 2 qs) (qs - 1)])/(6 (3 - 2 qs) (qs - 1));
chk["lem:radius: P_q(R(q)) = 0", P[qs, Rq[qs]], M > 0 && 1 < qs < 3/2];
chk["lem:radius: R(1+) = 8M", Limit[Rq[qs], qs -> 1, Direction -> "FromAbove"] - 8 M, M > 0];
chk["lem:radius: R(q) (3/2 - q) -> 23M/6", Limit[Rq[qs] (3/2 - qs), qs -> 3/2, Direction -> "FromBelow"] - 23 M/6, M > 0];
dR = Simplify[D[Rq[s] /. M -> 1, s]];
chk["lem:radius: R' > 0 on (1, 3/2) (Reduce on the rationalised numerator)",
  Resolve[Exists[s, 1 < s < 3/2 && dR <= 0], Reals] === False, True];
chk["lem:radius: d_q P = 3r[(4q-5)r + 11M]", D[P[qs, r], qs] - 3 r ((4 qs - 5) r + 11 M)];
chk["lem:radius: d_q^2 P = 12 r^2", D[P[qs, r], {qs, 2}] - 12 r^2];
chk["lem:radius: P_q(11M/(5-4q)) = -M^2(86q^2 - 72q + 1)/(4q-5)^2", P[qs, 11 M/(5 - 4 qs)] + M^2 (86 qs^2 - 72 qs + 1)/(4 qs - 5)^2];
chk["lem:radius: 86q^2 - 72q + 1 > 0 on (1, 5/4)", Reduce[Exists[qs, 1 < qs < 5/4 && 86 qs^2 - 72 qs + 1 <= 0], Reals] === False, True];
chk["lem:radius: d_r P at the zero = -M sqrt(b^2 + 160 a)",
  (D[P[qs, rr], rr] /. rr -> Rq[qs]) + M Sqrt[(33 qs - 38)^2 + 160*3 (3 - 2 qs) (qs - 1)], M > 0 && 1 < qs < 3/2];
(* Clairaut and the Fermat plane (constant A) *)
rho = r/(Sqrt[f] Sqrt[1 - f/qs]);
chk["Clairaut: C = -f rho_r/rho (constant q)", (Cq /. q[eta] -> qs) + f D[rho, r]/rho, assR && qs > 1];
w = M/rho;
chk["eq:fermat-plane: w = x sqrt((1-2x)(q-1+2x))/sqrt(q), x = M/r",
  (w /. r -> M/x) - x Sqrt[(1 - 2 x) (qs - 1 + 2 x)]/Sqrt[qs], 0 < x < 1/2 && qs > 1 && M > 0];
Fx = (1/D[rho, r])/(f Sqrt[1 - f/qs]);
ser = Series[(Fx /. r -> M/x) , {x, 0, 2}] // Normal;
(* F = 1 + A1 w + A2 w^2 + ...: expanded in w by inverting w(x) *)
wx = x Sqrt[(1 - 2 x) (qs - 1 + 2 x)]/Sqrt[qs];
xw = InverseSeries[Series[wx, {x, 0, 3}], ww] // Normal;
Fw = Series[(Fx /. r -> M/x) /. x -> xw, {ww, 0, 2}] // Normal // Simplify[#, qs > 1] &;
chk["eq:fermat-plane: F = 1 + A1 w + ..., A1 = sqrt(q)(2q-3)/(q-1)^(3/2)",
  Coefficient[Fw, ww, 1] - Sqrt[qs] (2 qs - 3)/(qs - 1)^(3/2), qs > 1];
chk["eq:fermat-plane: F(0) = 1", Coefficient[Fw, ww, 0] - 1, qs > 1];
chk["rmk:sharp32: coefficient A2 = q(15q^2 - 42q + 35)/(2(q-1)^3)",
  Coefficient[Fw, ww, 2] - qs (15 qs^2 - 42 qs + 35)/(2 (qs - 1)^3), qs > 1];
chk["rmk:sharp32: q -> oo: A1 -> 2, A2 -> 15/2 (deflection 4M/b + 15 pi M^2/(4 b^2))",
  {Limit[Sqrt[qs] (2 qs - 3)/(qs - 1)^(3/2), qs -> Infinity] - 2, Limit[qs (15 qs^2 - 42 qs + 35)/(2 (qs - 1)^3), qs -> Infinity] - 15/2} // Total[Abs[#]] &];
chk["rmk:sharp32: second-order deflection: 2 int_1^oo (A2 w^2)/(t sqrt(t^2-1)) = pi A2/2 (w = 1/(b t))",
  2 Integrate[1/(t^3 Sqrt[t^2 - 1]), {t, 1, Infinity}] - Pi/2];
(* prop:sharp32 *)
chk["prop:sharp32: 2 int_1^{D/b} dt/(t sqrt(t^2-1)) = 2 arccos(b/D)",
  2 Integrate[1/(t Sqrt[t^2 - 1]), {t, 1, Dd/b}, Assumptions -> Dd > b > 0] - 2 ArcCos[b/Dd], Dd > b > 0];
chk["prop:sharp32: A1 term: 2 int (A1 M/(b t)) dt/(t sqrt(t^2-1)) = 2 A1 M sqrt(D^2-b^2)/(b D)",
  2 A1 M/b Integrate[1/(t^2 Sqrt[t^2 - 1]), {t, 1, Dd/b}, Assumptions -> Dd > b > 0] - 2 A1 M Sqrt[Dd^2 - b^2]/(b Dd), Dd > b > 0];
chk["prop:sharp32: int_1^oo dt/(t^3 sqrt(t^2-1)) = pi/4", Integrate[1/(t^3 Sqrt[t^2 - 1]), {t, 1, Infinity}] - Pi/4];
Phi0[b_] := 2 ArcCos[b/Dd] - 2 cc M Sqrt[Dd^2 - b^2]/(b Dd);
chk["prop:sharp32: Phi_b (explicit part) = -2/sqrt(D^2-b^2) + 2cM sqrt(D^2-b^2)/(b^2 D) - 2cM/(D sqrt(D^2-b^2))",
  D[Phi0[b], b] - (-2/Sqrt[Dd^2 - b^2] + 2 cc M Sqrt[Dd^2 - b^2]/(b^2 Dd) + 2 cc M/(Dd Sqrt[Dd^2 - b^2])), Dd > b > 0];
Print["   note: the term 2cM/(D sqrt(D^2-b^2)) is absorbed in -2F(D)/sqrt(D^2-b^2) with F(D) = 1 - cM/D + G(D)"];
chk["prop:sharp32: limit D Phi_b -> -2 + 2c/lambda^2 with b = lambda sqrt(M D)",
  Limit[Dd ((-2 (1 - cc M/Dd)/Sqrt[Dd^2 - b^2] + 2 cc M Sqrt[Dd^2 - b^2]/(b^2 Dd)) /. b -> lam Sqrt[M Dd]), Dd -> Infinity] - (-2 + 2 cc/lam^2), M > 0 && lam > 0 && cc > 0];
chk["prop:sharp32: q = 6/5 => A1 = -3 sqrt(6)", (Sqrt[qs] (2 qs - 3)/(qs - 1)^(3/2) /. qs -> 6/5) + 3 Sqrt[6]];
chk["prop:sharp32: q = 6/5, F(x) = sqrt(1-2x)(1+10x)/(1+12x-40x^2)",
  ((Fx /. r -> M/x) /. qs -> 6/5) - Sqrt[1 - 2 x] (1 + 10 x)/(1 + 12 x - 40 x^2), 0 < x < 1/10 && M > 0];
chk["prop:sharp32: rational parametrisation of sqrt(1-2x), sqrt(1+10x) (same x)",
  Module[{s1 = (p^2 - 12 p + 30)/(p^2 - 10 p + 30), s2 = (30 - p^2)/(p^2 - 10 p + 30)}, (1 - s1^2)/2 - (s2^2 - 1)/10]];
xp = Module[{s1 = (p^2 - 12 p + 30)/(p^2 - 10 p + 30)}, Together[(1 - s1^2)/2]];
chk["prop:sharp32: x(0) = 0, x(11/10) >= 1/10 (covers r >= 10M; reaches r ~ 9.71M)", (xp /. p -> 0) == 0 && (xp /. p -> 11/10) >= 1/10, True];
chk["prop:sharp32: x(p) strictly increasing on [0, 11/10] (monotone covering)",
  Module[{dx = Together[D[xp, p]]}, Reduce[0 <= p <= 11/10 && dx <= 0, p, Reals] === False], True];
chk["prop:sharp32: sqrt(1-2x) and sqrt(1+10x) positive on [0, 11/10]",
  Resolve[Exists[p, 0 <= p <= 11/10 && ((p^2 - 12 p + 30) <= 0 || (30 - p^2) <= 0 || (p^2 - 10 p + 30) <= 0)], Reals] === False, True];
chk["prop:sharp32: rho(8200M)^2 = 6 8200^4 M^2/(8198 8210) > (2 10^4 M)^2",
  {(rho^2 /. {qs -> 6/5, r -> 8200 M}) - 6 8200^4 M^2/(8198 8210), 6 8200^4/(8198 8210) > (2 10^4)^2 // Boole // (1 - #) &} // Total, M > 0];

(* prop:sharp32, q = 6/5: certificate K = 465 in the parameter p (exact rational arithmetic, Reduce) *)
s1p = (p^2 - 12 p + 30)/(p^2 - 10 p + 30); s2p = (30 - p^2)/(p^2 - 10 p + 30);
wp = Together[xp s1p s2p/Sqrt[6]];          (* w = M/rho = x sqrt(1-2x) sqrt(1+10x)/sqrt(6) *)
Fp = Together[s1p (1 + 10 xp)/(1 + 12 xp - 40 xp^2)];
A1v = -3 Sqrt[6];
Gp = Together[Fp - 1 - A1v wp];               (* G as a function of p; rho = M/w *)
r2G = Together[Gp/wp^2];                      (* rho^2 G/M^2 *)
dGdrho = Together[D[Gp, p]/D[1/wp, p]];       (* dG/d(rho/M) *)
r3dG = Together[dGdrho/wp^3];                 (* rho^3 G'/M^2 *)
chk["prop:sharp32 cert: rho^2 G >= 0 for p in [0, 11/10]", Reduce[0 <= p <= 11/10 && r2G < 0, p, Reals] === False, True];
chk["prop:sharp32 cert: 465 M^2 - rho^2 G >= 0", Reduce[0 <= p <= 11/10 && 465 - r2G < 0, p, Reals] === False, True];
chk["prop:sharp32 cert: -rho^3 G' >= 0", Reduce[0 <= p <= 11/10 && -r3dG < 0, p, Reals] === False, True];
chk["prop:sharp32 cert: 930 M^2 + rho^3 G' >= 0", Reduce[0 <= p <= 11/10 && 930 + r3dG < 0, p, Reals] === False, True];
chkFalse["control: with K = 400 the first bound fails (must fail)", Reduce[0 <= p <= 11/10 && 400 - r2G < 0, p, Reals] === False, True];
Print["   max rho^2 G/M^2 = ", N[MaxValue[{r2G, 0 <= p <= 11/10}, p]], ",  max(-rho^3 G'/M^2) = ", N[MaxValue[{-r3dG, 0 <= p <= 11/10}, p]]];
(* bounds on D Phi_b at b_- = 300 M sqrt(D/D0), b_+ = 2000 M sqrt(D/D0), D >= D0 = 2 10^4 M (M = 1) *)
cv = 3 Sqrt[6]; Kv = 465; D0 = 2 10^4;
lowM[Dv_] := Module[{b = 300 Sqrt[Dv/D0]}, Dv (-2 (1 - cv/Dv + Kv/Dv^2)/Sqrt[Dv^2 - b^2] + 2 cv Sqrt[Dv^2 - b^2]/(b^2 Dv)) - Dv Pi Kv/b^3];
uppP[Dv_] := Module[{b = 2000 Sqrt[Dv/D0]}, Dv (-2 (1 - cv/Dv - Kv/Dv^2)/Sqrt[Dv^2 - b^2] + 2 cv Sqrt[Dv^2 - b^2]/(b^2 Dv)) + Dv Pi Kv/b^3];
mL = NMinValue[{lowM[Dv], Dv >= D0}, Dv, WorkingPrecision -> 30]; mU = NMaxValue[{uppP[Dv], Dv >= D0}, Dv, WorkingPrecision -> 30];
Print["   min_D D Phi_b(b_-) >= ", N[mL, 8], ";  max_D D Phi_b(b_+) <= ", N[mU, 8], ";  limits D -> oo: ", N[{Limit[lowM[Dv], Dv -> Infinity], Limit[uppP[Dv], Dv -> Infinity]}, 8]];
chk["prop:sharp32: D Phi_b(b_-) > 0.18 for every D >= D0", mL > 18/100, True];
chk["prop:sharp32: D Phi_b(b_+) < -1.9 for every D >= D0", mU < -19/10, True];
chk["prop:sharp32: Phi < pi on [b_-, b_+] (b <= D/2): -sqrt(3) c/b + pi K/(2 b^2) < 0 for b >= 300",
  Reduce[bb >= 300 && -Sqrt[3] cv/bb + Pi Kv/(2 bb^2) >= 0, bb, Reals] === False, True];
summary[];
