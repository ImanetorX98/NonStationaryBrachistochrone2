(* Paper III, finite drops at a large separation (theorem on finite drops, supplement S9): proper-time branch of
   Thakurta-Schwarzschild, M = 1, Ehat = Sqrt[2]. Parts: (1) single periastron, (2) frozen problem at A = 1/2,
   (3) metric and connection on the domain D, (4) canonical field in the layer, (5) late drops. *)
ok[name_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], name];
f[r_] := 1 - 2/r;
nN[r_, A_] := 1/(1 - A^2 f[r]/2);
om[r_, A_] := A^2 f[r]/Sqrt[2];
aH[r_, A_] := f[r]^2/nN[r, A];
bH[r_, A_] := f[r]/(nN[r, A] r^2);
S[r_, pr_, J_, A_] := Sqrt[aH[r, A] pr^2 + bH[r, A] J^2];
H[r_, pr_, J_, A_] := -om[r, A] + S[r, pr, J, A];
Rp[r_, A_] := A^2/Sqrt[2] Sqrt[nN[r, A]];            (* R_A' *)
rho[r_, A_] := A^2/Sqrt[2] r Sqrt[nN[r, A] f[r]];
RR[r_?NumericQ, A_] := NIntegrate[Rp[s, A], {s, 2, r}, WorkingPrecision -> 40];

Print["--- (1) single periastron"];
(* clock costate: q' = -dH/deta = -A' dH/dA, with q = -H, i.e. S = om - q *)
qp = Simplify[-D[H[r, pr, J, A], A] /. S[r, pr, J, A] -> (om[r, A] - q), {r > 2, A > 0}];
qpClaim = A f[r]/Sqrt[2] (nN[r, A] (om[r, A] - q)/Sqrt[2] + 2);
assum = {r > 2, 0 < A <= 1, J > 0, Element[pr, Reals]};
(* clock costate: dS/dA = -A f N S/2, hence -dH/dA = A f [sqrt2 + N S/2] = (A f/sqrt2)[2 + N S/sqrt2],
   that is q' = (A A' f/sqrt2)[N(om-q)/sqrt2 + 2] on H = 0 *)
ok["dS/dA = -A f N S/2", Simplify[D[S[r, pr, J, A], A] + A f[r] nN[r, A] S[r, pr, J, A]/2, assum] == 0];
ok["-dH/dA = (A f/sqrt2)[2 + N S/sqrt2]", Simplify[-D[H[r, pr, J, A], A] - A f[r]/Sqrt[2] (2 + nN[r, A] S[r, pr, J, A]/Sqrt[2]), assum] == 0];
(* p_r' = -H_r at p_r = 0 equals om_r - gamma S, gamma = (log bH)_r/2 *)
ok["p_r' at p_r=0 = om_r - gamma S",
   Simplify[(-D[H[r, pr, J, A], r] /. pr -> 0) - (D[om[r, A], r] - D[Log[bH[r, A]], r]/2 S[r, 0, J, A]), {r > 2, A > 0, J > 0}] == 0];
(* rho_A = om/sqrt(bH): comparison of the squares, both positive *)
ok["rho^2 = om^2/bH", Simplify[rho[r, A]^2 - om[r, A]^2/bH[r, A], assum] == 0];
ok["[N] rho = om/sqrt(bH) at 200 points", Max[Table[Abs[N[rho[r, A] - om[r, A]/Sqrt[bH[r, A]], 30]], {r, 21/10, 12, 1/2}, {A, 1/10, 1, 1/10}]] < 10^-25];
drr = Simplify[D[rho[r, A], r]/Rp[r, A], {r > 2, 0 < A <= 1}];
ok["drho/dR formula", Simplify[drr - ((r - 1)/(r Sqrt[f[r]]) + A^2 Sqrt[f[r]]/(r (2 - A^2 f[r]))), {r > 2, 0 < A <= 1}] == 0];
ok["drho/dR > 1 for r>2, 0<A<=1", Resolve[ForAll[{r, A}, r > 2 && 0 < A <= 1, (r - 1)^2 > r^2 f[r]]]];
ok["om_r - gamma om > 0", Resolve[ForAll[{r, A}, r > 2 && 0 < A <= 1,
     Simplify[D[om[r, A], r] - D[Log[bH[r, A]], r]/2 om[r, A]] > 0]]];

Print["--- (2) frozen A=1/2 (Q=8)"];
A0 = 1/2;
(* closed form of R_1/2(r) = int_2^r R_1/2' : R_1/2'^2 = r/(4(7r+2)) *)
Rc[r_] := (7 (Sqrt[r (7 r + 2)] - 4 Sqrt[2]) + 2 Sqrt[7] ArcTanh[Sqrt[7] (4 Sqrt[2] - Sqrt[r (7 r + 2)])/(7 r + 16)])/98;
ok["closed form: Rc(2) = 0 and Rc' = R_1/2' for r > 2", Rc[2] == 0 && FullSimplify[D[Rc[r], r] - Rp[r, A0], r > 2] === 0];
R6 = Rc[6]; Print["R_1/2(6) = ", N[R6, 12], "  (NIntegrate: ", N[RR[6, A0], 12], ")"];
ok["2R(6) > 254/(125 sqrt2) > 359/250 (closed form)", 2 R6 > 254/(125 Sqrt[2]) > 359/250];
ok["rho(6) = 3/sqrt11", Simplify[rho[6, A0] - 3/Sqrt[11]] == 0];
ok["3pi/(2sqrt11) < 33/(7sqrt11) < 711/500", 3 Pi/(2 Sqrt[11]) < 33/(7 Sqrt[11]) < 711/500];
ok["radial density < 1/5 on D", Maximize[{Rp[r, A0], 599/100 <= r <= 601/100}, r][[1]] < 1/5];
ok["143/200 < R(5.99), R(6.01) < 81/100 (closed form)", 143/200 < Rc[599/100] && Rc[601/100] < 81/100];
ok["L^2 > 1", 2 (143/200)^2 - 2 (81/100)^2 (1/100) > 1 && Sin[1/100] < 1/100];
Kc = Simplify[-D[rho[r, A0], {r, 2}]/(rho[r, A0] Rp[r, A0]^2) + D[rho[r, A0], r] D[Rp[r, A0], r]/(rho[r, A0] Rp[r, A0]^3), r > 2];
ok["K formula Q=8", Simplify[Kc - 4 (7 r^3 + 14 r^2 + 68 r + 8)/(r^3 (r - 2)^2 (7 r + 2)), r > 2] == 0];
ok["0 < K <= 1/2 on r>=4", Resolve[ForAll[r, r >= 4, 0 < 4 (7 r^3 + 14 r^2 + 68 r + 8)/(r^3 (r - 2)^2 (7 r + 2)) <= 1/2]]];
w4 = Rp[4, A0]/D[rho[r, A0], r] /. r -> 4;
ok["w(4)^2 = 450/529 > (9/10)^2 (exact)", Simplify[w4^2] == 450/529 && 450/529 > 81/100]; Print["   w(4) = ", N[w4, 10]];
ok["rho(4)^2/rho(5.99)^2 < (3/5)^2 (exact)", Simplify[rho[4, A0]^2/rho[599/100, A0]^2] < 9/25];
ok["rho(4)^2 = 4/15", Simplify[rho[4, A0]^2 - 4/15] == 0];
ok["cos(9/10) > 3/5 and 81/50 > pi/2+1/100", Cos[9/10] > 3/5 && 81/50 > Pi/2 + 1/100];
ok["J > 1/5", (Pi/2 - 1/100) (4/15)/(3/2) > 1/5];
(* [N] shallow periastron and angle of the dip at t=4, frozen problem r0=rf=6 *)
halfang[t_?NumericQ, rf_?NumericQ] := Module[{Jt = rho[t, A0]},
   NIntegrate[Jt Rp[s, A0]/(rho[s, A0] Sqrt[rho[s, A0]^2 - Jt^2]), {s, t, rf}, WorkingPrecision -> 30, Method -> "DoubleExponential"]];
Print["   [N] Delta(4) (r0=rf=6) = ", N[2 halfang[4, 6], 10]];
tsh = t /. FindRoot[2 halfang[t, 6] == Pi/2, {t, 9/2, 4.3, 5.9}, WorkingPrecision -> 25];
Print["   [N] shallow periastron at delta=pi/2: ", N[tsh, 10]];

Print["--- (3) metric and connection on D"];
grr = Together[Rp[r, A0]^2]; gpp = Together[rho[r, A0]^2];
ok["g_rr = r/(4(7r+2)), g_phiphi = r^2(r-2)/(4(7r+2))", Simplify[{grr - r/(4 (7 r + 2)), gpp - r^2 (r - 2)/(4 (7 r + 2))}] === {0, 0}];
inD[c_] := Resolve[ForAll[r, 599/100 <= r <= 601/100, c], Reals];
ok["1/32 < g_rr < 1/25 on D (Resolve)", inD[1/32 < grr < 1/25]];
ok["g_phiphi < 1 on D (Resolve)", inD[gpp < 1]];
ok["|G^r_rr| < 1/100 on D (Resolve)", inD[Abs[D[grr, r]/(2 grr)] < 1/100]];
ok["|G^phi_rphi| < 1/4 on D (Resolve)", inD[Abs[D[gpp, r]/(2 gpp)] < 1/4]];
ok["|G^r_phiphi| < 8 on D (Resolve)", inD[Abs[D[gpp, r]/(2 grr)] < 8]];
(* |Hess_ij| <= ||Hess|| sqrt(g_ii g_jj), ||Hess|| < 16/13, with the bounds just proved *)
hb = 16/13;
Vrr = hb (1/25) + 1/100 (1/5);
Vrp = hb Sqrt[(1/25) 1] + (1/4) 1;
Vpp = hb 1 + 8 (1/5);
Print["   bounds V_rr,V_rphi,V_phiphi: ", {Vrr, Vrp, Vpp}, " = ", N[{Vrr, Vrp, Vpp}, 5]];
ok["||D^2 V||_inf < 7/2 (rational bounds above)", Max[Vrr + Vrp, Vrp + Vpp] < 7/2];
ok["Jacobi: 1-s^3/12 and 1-s^2/4 at s=3/2", 3/2 - (3/2)^3/12 >= 13/16 (3/2) && 1 - (3/2)^2/4 >= 7/16];

Print["--- (4) canonical field in the layer"];
Y = {r, ph, pr, J};
G = {D[H[r, pr, J, A], pr], D[H[r, pr, J, A], J], -D[H[r, pr, J, A], r], 0};
DG = D[G, {Y}];
rowsum = Total /@ Abs[DG];
best = Table[NMaximize[{rowsum[[i]], 5 <= r <= 7, -1 <= pr <= 1, 1/10 <= J <= 2, 1/2 <= A <= 1}, {r, pr, J, A}, Method -> "DifferentialEvolution", WorkingPrecision -> 20][[1]], {i, 1, 3}];
Print["   [N] max row sums of |DG|: ", N[best, 6]];
ok["[N] ||DG||_inf < 300 (NMaximize)", Max[best] < 300];
ok["[N] ||G||_inf < 1", NMaximize[{Max[Abs /@ G], 5 <= r <= 7, -1 <= pr <= 1, 1/10 <= J <= 2, 1/2 <= A <= 1}, {r, pr, J, A}, WorkingPrecision -> 20][[1]] < 1];
(* exact bounds on the box 5 <= r <= 7, |p_r| <= 1, 1/10 <= J <= 2, 1/2 <= A <= 1, with d = 1 - A^2 f/2,
   a = f^2 d, b = f d/r^2, S^2 = a p_r^2 + b J^2; derivatives at fixed A *)
Print["   exact bounds for the canonical field:"];
dd = 1 - A^2 f[r]/2; aa = f[r]^2 dd; bb = f[r] dd/r^2;
ok["structure: aH = f^2 d, bH = f d/r^2, H = -om + Sqrt[aH pr^2 + bH J^2]",
  Simplify[{aH[r, A] - aa, bH[r, A] - bb}] === {0, 0}];
box = 5 <= r <= 7 && 1/2 <= A <= 1;
onBox[c_] := Resolve[ForAll[{r, A}, box, c], Reals];
ok["3/5 <= f <= 5/7, 9/14 <= d <= 1", onBox[3/5 <= f[r] <= 5/7 && 9/14 <= dd <= 1]];
ok["a <= 25/49, 1/128 < b <= 1/25", onBox[aa <= 25/49 && 1/128 < bb <= 1/25]];
ok["1/120 < S < 1: S^2 >= b J^2 > (1/128)(1/100) > 1/120^2 and S^2 <= 25/49 + 4/25 < 1", (1/128) (1/100) > 1/14400 && 25/49 + 4/25 < 1];
ok["|f'| <= 2/25, |f''| <= 4/125, |d'| <= 1/25, |d''| <= 2/125", onBox[Abs[D[f[r], r]] <= 2/25 && Abs[D[f[r], {r, 2}]] <= 4/125 && Abs[D[dd, r]] <= 1/25 && Abs[D[dd, {r, 2}]] <= 2/125]];
ok["|(log a)'|, |(log b)'| < 3/5 and |(log a)''|, |(log b)''| < 1 (exact on the box)",
  onBox[Abs[D[Log[aa], r]] < 3/5 && Abs[D[Log[bb], r]] < 3/5 && Abs[D[Log[aa], {r, 2}]] < 1 && Abs[D[Log[bb], {r, 2}]] < 1]];
F1 = (2/25)/(3/5); D1 = (1/25)/(9/14); F2 = (4/125)/(3/5) + F1^2; D2 = (2/125)/(9/14) + D1^2;
ok["the same from the elementary bounds: 2F1 + D1, F1 + D1 + 2/5 < 3/5; 2F2 + D2, F2 + D2 + 2/25 < 1",
  2 F1 + D1 < 3/5 && F1 + D1 + 2/5 < 3/5 && 2 F2 + D2 < 1 && F2 + D2 + 2/25 < 1];
(* with the weights al = a p_r^2/S^2, be = b J^2/S^2 (al + be = 1):
   S_r/S = (al (log a)' + be (log b)')/2, S_rr/S = (al a''/a + be b''/b)/2 - (S_r/S)^2, c''/c = (log c)'' + (log c)'^2 *)
SS = Sqrt[aa pr^2 + bb J^2];
ok["S_r/S = (al (log a)' + be (log b)')/2", Simplify[D[SS, r]/SS - (aa pr^2 D[Log[aa], r] + bb J^2 D[Log[bb], r])/(2 SS^2), r > 2] === 0];
ok["S_rr/S = (al a''/a + be b''/b)/2 - (S_r/S)^2", Simplify[D[SS, {r, 2}]/SS - ((aa pr^2 D[aa, {r, 2}]/aa + bb J^2 D[bb, {r, 2}]/bb)/(2 SS^2) - (D[SS, r]/SS)^2), r > 2] === 0];
ok["|S_r| < 3/10, |S_rr| < (1 + 9/25)/2 + 9/100 < 1 - 4/125 (S < 1)", (3/5)/2 == 3/10 && (1 + 9/25)/2 + (3/10)^2 < 1 - 4/125];
ok["|om_r| <= 2/25, |om_rr| <= 4/125, so |H_r| < 1, |H_rr| < 1", onBox[Abs[D[om[r, A], r]] <= 2/25 && Abs[D[om[r, A], {r, 2}]] <= 4/125] && 2/25 + 3/10 < 1];
ok["|H_p| <= sqrt(c) and H_rp = H_p ((log c)' - S_r/S): |H_rp| <= (3/5 + 3/10) sqrt(25/49) < 1",
  Simplify[D[aa pr/SS, r] - (aa pr/SS) (D[Log[aa], r] - D[SS, r]/SS), r > 2] === 0 &&
  Simplify[D[bb J/SS, r] - (bb J/SS) (D[Log[bb], r] - D[SS, r]/SS), r > 2] === 0 && (3/5 + 3/10) (5/7) < 1];
Dm = DiagonalMatrix[{Sqrt[aa], Sqrt[bb]}]; xi = Dm . {pr, J}/SS;
ok["D^2_pp S = S^-1 D (I - xi xi^T) D", Simplify[D[SS, {{pr, J}, 2}] - Dm . (IdentityMatrix[2] - Outer[Times, xi, xi]) . Dm/SS, r > 2] === {{0, 0}, {0, 0}}];
ok["entries of D^2_pp S at most max(a,b)/S < 120 (25/49); row sums of DG: 1 + 2 (120)(25/49) < 124 < 300, and 3", 1 + 2 120 25/49 < 124 && 3 < 124];
ok["||G||_inf < 1: |H_pr| <= 5/7, |H_J| <= 1/5, |H_r| < 1", 5/7 < 1 && 1/5 < 1];
ok["4(e^(3/100)-1) <= 12/97 < 1", 4 (Exp[3/100] - 1) <= 12/97 < 1];
ok["|Omega dot| <= 1/(3 sqrt3) < 1/5", Maximize[{Sqrt[f[r]]/r, r > 2}, r][[1]] == 1/(3 Sqrt[3]) && 1/(3 Sqrt[3]) < 1/5];
ok["waiting cost < 711/500 + h/2 < 359/250 (h=1e-4)", om[6, 1] < 1/2 && 711/500 + 1/20000 < 359/250];

Print["--- (5) late drops"];
ok["R_1' = sqrt(r/(r+2)) in [1/sqrt2, sqrt3/2] on [2,6]", Simplify[Rp[r, 1] - Sqrt[r/(r + 2)], r > 2] == 0 && Sqrt[6/8] == Sqrt[3]/2];
ok["density R_1/2 <= R_1/4", Resolve[ForAll[r, r > 2, Rp[r, 1/2] <= Rp[r, 1]/4]]];
ok["sqrt2 > 5/4 and (5/4)2sqrt3 < 13/3", Sqrt[2] > 5/4 && 5 Sqrt[3]/2 < 13/3];
ok["(13/3)sqrt2 < 13/2", 13/3 Sqrt[2] < 13/2];
ok["4 sqrt6/3 < 4", 4 Sqrt[6]/3 < 4];
ok["52/(10^7-13/2) < 10^-5", 52/(10^7 - 13/2) < 1/100000];
re = 2 + 1/100000;
ok["om_1(r(d)) < 1/200000", om[re, 1] < 1/200000];
ok["rho_1^2 = r^2(r-2)/(r+2)", Simplify[rho[r, 1]^2 - r^2 (r - 2)/(r + 2), r > 2] == 0];
ok["rho_1(2+1e-5) < 1/300", rho[re, 1] < 1/300];
rs = 201/100;
ok["(r sqrt f/(r-1))^2 = 201/10201 < 1/49 at r_s", rs^2 f[rs]/(rs - 1)^2 == 201/10201 && 201/10201 < 1/49];
ok["rho_1/2(r_s) > 1/40", rho[rs, 1/2] > 1/40]; Print["   rho_1/2(2.01) = ", N[rho[rs, 1/2], 8]];
ok["rho_1(r_s) > 1/10, om_1(r_s) > 1/400", rho[rs, 1] > 1/10 && om[rs, 1] > 1/400];
ok["(1/10)(1-(1/200000)/(1/400)) > 1/40", (1/10) (1 - 400/200000) > 1/40];
ok["arcsin(2/15) < 7/50", ArcSin[2/15] < 7/50];
ok["2pi/7+14/25+1/500 < pi/2", 2 Pi/7 + 14/25 + 1/500 < Pi/2];
(* w_q <= w_1 for q in [0, om_1): identity d rho_q/dR_q - d rho_1/dR_1 = rho_1 q om'/(om^2 R_1'(1-q/om)); all factors > 0 *)
rq = rho[r, 1] (1 - q/om[r, 1]); Rq = Rp[r, 1] (1 - q/om[r, 1]);
diff = D[rq, r]/Rq - D[rho[r, 1], r]/Rp[r, 1];
claim = rho[r, 1] q D[om[r, 1], r]/(om[r, 1]^2 Rp[r, 1] (1 - q/om[r, 1]));
ok["identity of the comparison w_q <= w_1", Simplify[diff - claim, {r > 2, 0 <= q < om[r, 1]}] == 0];
ok["om_1' > 0", Simplify[D[om[r, 1], r] > 0, r > 2]];
ok["[N] difference >= 0 on a grid", Min[Table[N[diff /. q -> s om[r, 1], 30], {r, 201/100, 12, 1/4}, {s, 0, 99/100, 1/20}]] >= 0];
ok["R_1' = rho_1/(r sqrt f)", Simplify[Rp[r, 1] - rho[r, 1]/(r Sqrt[f[r]]), r > 2] == 0];
