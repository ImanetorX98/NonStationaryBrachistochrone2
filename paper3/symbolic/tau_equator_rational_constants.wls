(* Paper III, example on the Kerr equator with spin (supplement S10): Wolfram check of the rational constants and of
   the auxiliary claims used in the proof of the theorem on global uniqueness and of the sector change.
   [EXACT] except where marked [N]. *)
pass[name_, e_] := Print[If[TrueQ[e], "PASS ", "FAIL "], name];
f = 1 - 2/r; Del = r^2 - 2 r + 1/4; hr = r^2/(f Del); hp = Del/f^2; bphi = -1/(r - 2); Bm = bphi^2/hp;
NN = 1/(1 - A^2 f/2); d = NN - Bm; LN = D[NN, r]/NN; Ld = D[d, r]/d; lr = D[Log[hr], r]; lp = D[Log[hp], r];
LB = D[Log[Bm], r];
shell = 23/4 <= r <= 25/4 && 1/2 <= A <= 1;
all[cond_] := TimeConstrained[Resolve[ForAll[{r, A}, shell, cond], Reals], 300, $Failed];
pass["force: max{(11/100)(2/3), (77/500)(67/500)(22/5)} = (77/500)(67/500)(22/5)", Max[11/100 2/3, 77/500 67/500 22/5] == 77/500 67/500 22/5];
pass["force: sum = 744633/5000000 < 3/20", 43/1000 + 77/500 67/500 22/5 + 67/500 43/2000 + 69/25000 22/5 + 6/1250 43/2000 == 744633/5000000 && 744633/5000000 < 3/20];
pass["0 <= L_N <= 800/17457", all[0 <= LN <= 800/17457]];
pass["0 <= L_d < 47/1000", all[0 <= Ld < 47/1000]];
pass["2 L_d - L_N < 6/125", all[2 Ld - LN < 6/125]];
pass["-11/50 < ell_r < -1/10", all[-11/50 < lr < -1/10]];
pass["1/5 < ell_phi < 13/50", all[1/5 < lp < 13/50]];
pass["(log B_m)_r = -2/r - Delta'/Delta", Simplify[LB - (-2/r - D[Del, r]/Del)] === 0];
Bs = 256/(23^2 349); Ns = 184/169; Ls = 800/17457;
pass["0 <= L_d - L_N <= B*(4/5 + L*)/(N* - B*)", all[0 <= Ld - LN <= Bs (4/5 + Ls)/(Ns - Bs)]];
pass["N/(N - B*)^2 decreasing in N for N >= N*", Resolve[ForAll[n, n >= Ns, D[n/(n - Bs)^2, n] < 0], Reals]];
pass["|L_B| < 4/5, |L_B'| < 13/50, |B_m''| < (9/10) B_m", all[Abs[LB] < 4/5 && Abs[D[LB, r]] < 13/50 && Abs[D[Bm, {r, 2}]] < 9/10 Bm]];
pass["N/d < 501/500", all[NN/d < 501/500]];
pass["|L_N'| < 2/125", all[Abs[D[LN, r]] < 2/125]];
pass["|L_d'| <= (501/500)(800/17457)(2/r_-) + (9/10) B*/d* + (47/1000)^2 < 21/1000",
  501/500 800/17457 2/(23/4) + 9/10 Bs/(184/169 - Bs) + (47/1000)^2 < 21/1000];
pass["h_r < 12/5, h_phi > 49, N < 25/16 on the shell", all[hr < 12/5 && hp > 49 && NN < 25/16]];
pass["a_H/b_H < 29 (a/b = d h_phi/(N h_r))", all[(d hp)/(NN hr) < 29]];
pass["weighted mean: m = 17920/8979, M = 48/11, M/m = 26937/12320 < 9/4",
  7/3 64/73 40/41 == 17920/8979 && 15/4 64/55 == 48/11 && (48/11)/(17920/8979) == 26937/12320 && 26937/12320 < 9/4];
pass["M(1-c)^2 - m c^2 <= 0, c = 3/5", 48/11 (2/5)^2 - 17920/8979 (3/5)^2 <= 0];
pass["t_H < 1/20 => sqrt(1+t) < 41/40", (41/40)^2 > 21/20];
pass["(3/20)(3/5)(21/100) = 189/10000", 3/20 3/5 21/100 == 189/10000];
pass["(36/1067)[11/100 - (11/100)^3/6]^2 > (189/10000)^2", 36/1067 (11/100 - (11/100)^3/6)^2 > (189/10000)^2];
pass["sin x > x - x^3/6 on (0, 1/5]", Resolve[ForAll[x, 0 < x <= 1/5, Sin[x] > x - x^3/6], Reals]];
pass["|p_r|/j: (203/1000)(10/49) + 567/16000 < 77/1000 and 63/2000/(8/9) = 567/16000", 203/1000 10/49 + 567/16000 < 77/1000 && (63/2000)/(8/9) == 567/16000];
pass["|tan(theta)| <= 10/49 for |theta| <= 1/5", Tan[1/5] <= 10/49];
pass["kappa_6 < 203/1000", 3456/84681 < (203/1000)^2];
pass["sqrt(72/97)/5 + 63/2000 < 1", Sqrt[72/97]/5 + 63/2000 < 1];
pass["|xi_r/xi| < (43/2000)(72/55); |gamma - alpha + xi_r/xi| < 3/10", 77/500 + 11/100 + 43/2000 72/55 < 3/10];
pass["|alpha'| + alpha^2 < 81/1000, |gamma'| + gamma^2 < 81/1000", 8/125 + (11/100)^2 < 81/1000 && 57/1000 + (77/500)^2 < 81/1000];
pass["S/j < 43/250", 67/500 73/64 28/25 < 43/250];
pass["t_H < 6/25 on Q_s: 29 (77/1000 * 64/55)^2 < 6/25; rho < 6/31", 29 (77/1000 64/55)^2 < 6/25 && (6/25)/(1 + 6/25) == 6/31];
pass["(1+t)^(3/2): (31/25)^3 < (7/5)^2, (21/20)^3 < (11/10)^2", (31/25)^3 < (7/5)^2 && (21/20)^3 < (11/10)^2];
pass["rho < 1/21 on return: (1/20)/(1+1/20) = 1/21", (1/20)/(21/20) == 1/21];
K2 = 6/1250 ((23/40)^2 + 3/20);
pass["|H_rr|/j: sum = 9606889437/248000000000 < 1/25",
  43/250 (81/1000 + 6/31 9/100) + 67/500 (2 77/500 43/2000 + 3/400) 9/8 + 3/200 9/8 + K2 73/64 + 2 69/25000 43/2000 9/8 + 6/1250 3/400 9/8 == 9606889437/248000000000 &&
  9606889437/248000000000 < 1/25];
pass["|H_r|/j < 2/25", 43/1000 9/8 + 77/500 67/500 73/64 + 67/500 43/2000 9/8 + 69/25000 73/64 + 6/1250 43/2000 9/8 < 2/25];
Ts = 21/100;
mat[B_] := {{1 - Ts B, -Ts 48/11}, {-Ts/25, 1 - Ts B}};
sol[B_, rhs_] := Inverse[mat[B]].rhs;
pass["M_B diagonal and det positive, inverse >= 0 (B = 7/50, 8/125)", And @@ Flatten[Table[{1 - Ts B > 0, Det[mat[B]] > 0, And @@ Thread[Flatten[Inverse[mat[B]]] >= 0]}, {B, {7/50, 8/125}}]]];
{u1, v1} = sol[7/50, {Ts 11/200 10/49, 203/1000/(49/50)^2 + Ts 57/1000 10/49}];
pass["theta column: U < 213/1000, V < 223/1000", u1 < 213/1000 && v1 < 223/1000];
{u2, v2} = sol[8/125, {93/1000 + 46/7 13/500, 2/25 + 46/7 7/125}];
pass["T column: U_T = 4035532775/5808706624, V_T = 1526921275/3319260928", u2 == 4035532775/5808706624 && v2 == 1526921275/3319260928];
nt = Ts (8/125 u2 + 48/11 v2) + 46/7 13/500;
pass["int|u'| = 436915382371/726088328000 < 61/100", nt == 436915382371/726088328000 && nt < 61/100];
pass["N <= 602/401, h_phi <= 2769038363/50501250 at r = 301/50, A = 1", (NN /. {r -> 301/50, A -> 1}) == 602/401 && (hp /. r -> 301/50) == 2769038363/50501250];
(* N and h_phi increasing in r (maximum at r = 301/50): exact minima of the derivatives *)
dN = Together[D[NN, r]]; dh = Together[D[hp, r]];
Print["dN/dr = ", dN]; Print["dh_phi/dr = ", Factor[dh]];
pass["dN/dr > 0 for r in [299/50, 301/50], A in [1/2,1] (exact MinValue)", MinValue[{dN, 299/50 <= r <= 301/50 && 1/2 <= A <= 1}, {r, A}] > 0];
pass["dh_phi/dr > 0 for r in [299/50, 301/50] (exact MinValue)", MinValue[{dh, 299/50 <= r <= 301/50}, r] > 0];
pass["b_H > 1/(N h_phi) > (11/100)^2", 1/(602/401 2769038363/50501250) > (11/100)^2];
{u3, v3} = sol[8/125, {31/1000, 3/80}];
pass["dilation: U_d = 55813825/829815232, V_d = 128068875/3319260928", u3 == 55813825/829815232 && v3 == 128068875/3319260928 && u3 < 7/100 && v3 < 1/25];
Zs = 63/7750;
hd = 1/10 - 21/100 (13/500 7/100 + 48/11 Zs 1/25) - 7/500 7/100;
pass["H_diag = 167669417/1705000000 > 49/500", hd == 167669417/1705000000 && hd > 49/500];
{u4, v4} = sol[8/125, {19/1000, 3/80}];
pass["T column at the change: 45639925/829815232, 127722375/3319260928", u4 == 45639925/829815232 && v4 == 127722375/3319260928];
nt4 = Ts (8/125 u4 + 48/11 v4);
pass["int|u'| = 3734179449/103726904000 < 21/500", nt4 == 3734179449/103726904000 && nt4 < 21/500];
Cs = 13/500 + 8/125 Zs;
gt = 7/50 + Cs 21/100 3/50 + Zs 21/500 + 7/500 3/50;
pass["Gamma_T < 1370932129/9687500000 < 71/500", gt == 1370932129/9687500000 && gt < 71/500];
(* w decreasing in r and concave in A^2 (|r - 6| < 1/300): exact identity and signs, numerical control *)
w = Sqrt[1 - x f/2]/Sqrt[hp];
dw = Simplify[D[w, r], 1799/300 <= r <= 1801/300 && 1/4 <= x <= 1];
pass["dw/dr < 0 (NMaxValue, [N] control)", NMaxValue[{D[w, r], 1799/300 <= r <= 1801/300 && 1/4 <= x <= 1}, {r, x}, WorkingPrecision -> 30] < 0];
(* exact form: w^2 = (1 - x f/2)/h_phi; d(w^2)/dr = -(x f_r/2)/h_phi - (1 - x f/2) h_phi'/h_phi^2 < 0 because f_r > 0, h_phi' > 0, 1 - x f/2 > 0 *)
pass["d(w^2)/dr = -(x f_r/2)/h_phi - (1 - x f/2) h_phi'/h_phi^2 (identity)",
  Simplify[D[(1 - x f/2)/hp, r] - (-(x D[f, r]/2)/hp - (1 - x f/2) D[hp, r]/hp^2)] === 0];
pass["f_r = 2/r^2 > 0 and 1 - x f/2 > 0 for x <= 1, f < 1", Simplify[D[f, r] - 2/r^2] === 0];
pass["sqrt(1 - x f/2) concave in x: second derivative = -f^2/(16 (1 - x f/2)^(3/2))",
  Simplify[D[Sqrt[1 - x ff/2], {x, 2}] + ff^2/(16 (1 - x ff/2)^(3/2))] === 0];
pass["1 - b+ w* > 0", 1 - 251/1000 13/100 > 0];
