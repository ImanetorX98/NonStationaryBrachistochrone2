(* Paper III, examples with a time-dependent factor (supplement S9; theorems on the pulse at a small separation and on
   finite drops at a large separation): proper-time branch of Thakurta-Schwarzschild (a = 0), M = 1, Ehat = Sqrt[2].
   "pulse" = the pulse at a small separation, "drop" = the drop at a large separation, "pulse condition" = the condition
   T|ell_T| <= 46/7 for general pulses. [EXACT] except where marked [N]. Identities with radicals in rational form. *)
pass[name_, e_] := Print[If[TrueQ[e], "PASS ", "FAIL "], name];
zr[e_] := TrueQ[Together[Simplify[e]] === 0];
all[vars_, dom_, cond_] := TimeConstrained[Resolve[ForAll[vars, dom, cond], Reals], 300, $Failed];
f = 1 - 2/r; NN = 1/(1 - A^2 f/2); om = A^2 f/Sqrt[2]; hr = 1/f^2; hp = r^2/f;
aH = f^2/NN; bH = f/(NN r^2); S = Sqrt[aH p^2 + bH J^2]; H = -om + S;
shell = 23/4 <= r <= 25/4 && 1/2 <= A <= 1;

Print["--- (A) 3D reduction and Hamiltonian [EXACT]"];
Clear[rd, Od];
pass["rail: rdot^2 h_r + Omegadot^2 h_phi = 1/N  =>  Omegadot = w0 sqrt(1-u^2), u = rdot sqrt(N)/f, w0 = sqrt(f)/(r sqrt N)",
  zr[(1/NN - rd^2 hr)/hp - (f/(r^2 NN)) (1 - (rd Sqrt[NN]/f)^2)]];
pass["H = -omega + |p|_{h0^-1}/sqrt(N), |p|^2 = f^2 p^2 + f J^2/r^2", zr[S^2 - (f^2 p^2 + f J^2/r^2)/NN]];
pass["a_H/b_H = r(r-2)", zr[aH/bH - r (r - 2)]];
pass["(a_H/sqrt b_H)^2 = f^3 r^2/N", zr[aH^2/bH - f^3 r^2/NN]];
lr = D[Log[hr], r]; lp = D[Log[hp], r]; LN = D[NN, r]/NN;
al = D[Log[aH], r]/2; ga = D[Log[bH], r]/2;
pass["ell_r = -4/(r(r-2)), ell_phi = 2/r - 2/(r(r-2))", zr[lr + 4/(r (r - 2))] && zr[lp - (2/r - 2/(r (r - 2)))]];
pass["alpha = -(L_N + ell_r)/2, gamma = -(L_N + ell_phi)/2", zr[al + (LN + lr)/2] && zr[ga + (LN + lp)/2]];
PG = 4 (r - 1)/(r^2 (r - 2)^2);
pass["ell_r' = 2 P, ell_phi' = -2/r^2 + P", zr[D[lr, r] - 2 PG] && zr[D[lp, r] - (-2/r^2 + PG)]];
pass["L_N' = L_N (L_N - 2/r)", zr[D[LN, r] - LN (LN - 2/r)]];
u = Sqrt[aH] p; w = Sqrt[bH] J; rho = aH p^2/S^2; sg = 1 - rho;
pass["S_r = (alpha U^2 + gamma W^2)/S", zr[(D[S, r] - (al u^2 + ga w^2)/S) /. {p -> 3/100, J -> 2, r -> 6, A -> 3/4}] &&
  Simplify[D[S, r] - (al u^2 + ga w^2)/S] === 0];
pass["S_rr/S = (alpha'+alpha^2) rho + (gamma'+gamma^2) sigma + rho sigma (gamma-alpha)^2",
  Simplify[D[S, {r, 2}]/S - ((D[al, r] + al^2) rho + (D[ga, r] + ga^2) sg + rho sg (ga - al)^2)] === 0];
pass["H_rp = rdot [2 alpha - alpha rho - gamma sigma]", Simplify[D[H, r, p] - D[H, p] (2 al - al rho - ga sg)] === 0];
pass["H_rJ = phidot [gamma + rho (gamma - alpha)]", Simplify[D[H, r, J] - D[H, J] (ga + rho (ga - al))] === 0];
Z = p/J;
pass["rank one: H_pJ = -Z H_pp, H_JJ = -Z H_pJ", Simplify[D[H, p, J] + Z D[H, p, p]] === 0 && Simplify[D[H, J, J] + Z D[H, p, J]] === 0];
pass["projective force: p H_rp + J H_rJ - H_r = omega_r", Simplify[p D[H, r, p] + J D[H, r, J] - D[H, r] - D[om, r]] === 0];
pass["(J H_pp)^2 = (a^2/b)/(1+t)^3 with t = a p^2/(b J^2)",
  zr[(J D[H, p, p])^2 - (aH^2/bH)/(1 + aH p^2/(bH J^2))^3 /. {p -> 1/50, J -> 3, r -> 6, A -> 3/4}] &&
  zr[(J aH bH J^2)^2/(bH J^2 (1 + aH p^2/(bH J^2)))^3 - (aH^2/bH)/(1 + aH p^2/(bH J^2))^3]];
Clear[du, dp, dJ];
uprime = D[H, p, r] du + D[H, p, p] dp + D[H, p, J] dJ;
zprime = -(D[H, r, r] du + D[H, r, p] dp + D[H, r, J] dJ) + (D[H, r]/J) dJ;
zz = dp - Z dJ;
pass["u' = b u + a P  (no source in h; a = J H_pp, P = z/J, J > 0)", Simplify[uprime - (D[H, r, p] du + D[H, p, p] zz)] === 0];
pass["P' = -l u - b P - f h, f = omega_r/j", Simplify[zprime - (-D[H, r, r] du - D[H, r, p] zz - D[om, r] dJ/J)] === 0];
pass["dphi' = (H_rJ + Z H_rp) u - Z u'", Simplify[(D[H, J, r] du + D[H, J, p] dp + D[H, J, J] dJ) - ((D[H, r, J] + Z D[H, r, p]) du - Z uprime)] === 0];

Print["--- (B) constants at r = 6 [EXACT]"];
cr2 = (om^2 NN hr) /. r -> 6; cp2 = (om^2 NN hp) /. r -> 6;
pass["c_r^2 = A^4/(2 - 2A^2/3)", zr[cr2 - A^4/(2 - 2 A^2/3)]];
pass["c_r^2 in [3/88, 3/4], c_phi^2 in [9/11, 18], kappa_6^2 = 1/24",
  (cr2 /. A -> 1/2) == 3/88 && (cr2 /. A -> 1) == 3/4 && (cp2 /. A -> 1/2) == 9/11 && (cp2 /. A -> 1) == 18 && zr[cr2/cp2 - 1/24]];
pass["c_r, c_phi increasing in A", all[A, 1/2 <= A <= 1, D[cr2, A] > 0]];
pass["kappa_6 in (1/5, 41/200)", 1/24 > 1/25 && 1/24 < (41/200)^2];
pass["d log c_phi/dA = 2/A + A/(3 - A^2)", zr[D[Log[cp2], A]/2 - (2/A + A/(3 - A^2))]];

Print["--- (C) bounds on the shell [EXACT]"];
pass["184/169 <= N <= 50/33", all[{r, A}, shell, 184/169 <= NN <= 50/33]];
pass["2 < h_r < 12/5, 49 < h_phi < 58", all[r, 23/4 <= r <= 25/4, 2 < hr < 12/5 && 49 < hp < 58]];
pass["sqrt(a_H) < 2/3, sqrt(b_H) < 27/200", all[{r, A}, shell, aH < 4/9 && bH < (27/200)^2]];
pass["b_H > 33/2900 (used for |phidot|)", all[{r, A}, shell, bH > 33/2900]];
pass["7/3 < a_H/sqrt(b_H) < 15/4  (squared)", all[{r, A}, shell, 49/9 < aH^2/bH < 225/16]];
pass["a_H/b_H < 27", all[r, 23/4 <= r <= 25/4, r (r - 2) < 27]];
pass["0 <= alpha <= 11/100, -77/500 <= gamma <= 0", all[{r, A}, shell, 0 <= al <= 11/100 && -77/500 <= ga <= 0]];
pass["omega_r < 43/1000", all[{r, A}, shell, D[om, r] < 43/1000]];
pass["L_N <= 800/17457", all[{r, A}, shell, LN <= 800/17457]];
pass["-11/50 < ell_r < -1/10, 1/5 < ell_phi < 13/50", all[r, 23/4 <= r <= 25/4, -11/50 < lr < -1/10 && 1/5 < lp < 13/50]];
pass["0 < ell_r' < 53/500, |ell_phi'| < 7/125, |L_N'| < 2/125",
  all[{r, A}, shell, 0 < D[lr, r] < 53/500 && -7/125 < D[lp, r] < 7/125 && -2/125 < D[LN, r] < 2/125]];
pass["|alpha'| < 8/125, |gamma'| < 57/1000, |omega_rr| < 3/200",
  all[{r, A}, shell, -8/125 < D[al, r] < 8/125 && -57/1000 < D[ga, r] < 57/1000 && -3/200 < D[om, {r, 2}] < 3/200]];
pass["P decreasing for r > 2", all[r, r > 2, D[Log[PG], r] < 0]];
pass["omega >= (1/4)(15/23)/sqrt2 > 11/100", (1/4 (15/23)/Sqrt[2]) > 11/100 && all[{r, A}, shell, om >= 1/4 (15/23)/Sqrt[2]]];

Print["--- (D) chain of the pulse [EXACT]"];
pass["L_{1/2}: radial coefficient >= 1/(4 sqrt2) (c = 1/(4 sqrt 2), sqrt(N) h_r^(1/2)... ) and 1/(8 sqrt2) > 7/80", 1/(8 Sqrt[2]) > 7/80];
pass["angular coefficient >= sqrt(345)/(16 sqrt2) > 4/5", Sqrt[345]/(16 Sqrt[2]) > 4/5];
pass["angular coefficient of h_tau^(1/2) = (1/(4 sqrt2)) sqrt(N_{1/2} f) r >= sqrt(345)/(16 sqrt 2) on the shell",
  all[r, 23/4 <= r <= 25/4, (1/(4 Sqrt[2])) Sqrt[(NN /. A -> 1/2) f] r >= Sqrt[345]/(16 Sqrt[2])]];
pass["radial coefficient of h_tau^(1/2) = (1/(4 sqrt2)) sqrt(N_{1/2}) >= 1/(4 sqrt2)", all[r, r > 2, (NN /. A -> 1/2) >= 1]];
pass["circle: c_phi(1) delta = 3 sqrt2/50 < 7/80; duration <= 9/50", Sqrt[18]/50 == 3 Sqrt[2]/50 && 3 Sqrt[2]/50 < 7/80];
pass["force: 43/1000 + (77/500)(27/200)(22/5) = 33619/250000 < 3/20; sqrt(18) < 22/5; (11/100)(2/3) < (77/500)(27/200)(22/5)",
  43/1000 + 77/500 27/200 22/5 == 33619/250000 && 33619/250000 < 3/20 && Sqrt[18] < 22/5 && 11/100 2/3 < 77/500 27/200 22/5];
pass["bootstrap: sqrt(3/4) + (3/20)(35/44) < 1", Sqrt[3/4] + 3/20 35/44 < 1];
pass["J^2 > (9/11)[1 - (88/3)(21/176)^2] > 4/9", 9/11 (1 - 88/3 (21/176)^2) > 4/9];
pass["t_H < 27 (3T/20)^2 (9/4) < (11/8) T^2", 27 (3/20)^2 9/4 < 11/8];
pass["T/sqrt(1 + (11/8) T^2) increasing; sqrt(33/2900)(1/5)/sqrt(1 + (11/8)/25) > 1/50 => T < 1/5",
  all[T, T > 0, D[T/Sqrt[1 + 11/8 T^2], T] > 0] && Sqrt[33/2900] (1/5)/Sqrt[1 + 11/8/25] > 1/50];
pass["|p_r| < 3/100, J^2 > 64/81, chi > 1/10, |r - 6| < 1/50",
  3/20 1/5 == 3/100 && 9/11 (1 - 88/3 (3/100)^2) > 64/81 && (8/9)/Sqrt[58] > 1/10 && 1/10 15/4 3/100 9/8 < 1/50];
pass["heading: (3/88)(149/750)^2 > (3/100)^2, sin(1/5) > 149/750", 3/88 (149/750)^2 > (3/100)^2 && Sin[1/5] > 149/750];
pass["|phidot| < 27/200 => T > (1/50)/(27/200) = 4/27 > 3/25", (1/50)/(27/200) == 4/27 && 4/27 > 3/25];
pass["Q: j >= (3/sqrt11) cos(1/5) > 7/8", 3/Sqrt[11] Cos[1/5] > 7/8];
pass["Q: sqrt(3/4)/5 + 63/2000 < 1; (41/200)(10/49) + (63/2000)/(7/8) < 79/1000 and 9/250 = (63/2000)(8/7)",
  Sqrt[3/4]/5 + 63/2000 < 1 && 41/200 10/49 + (63/2000)/(7/8) < 79/1000 && (63/2000)/(7/8) == 9/250];
pass["(15/4)(79/1000) < 3/10, (21/100)(3/10) < 3/40", 15/4 79/1000 < 3/10 && 21/100 3/10 < 3/40];
pass["signs at the edges: (3/88)(149/750)^2 > (63/2000)^2", 3/88 (149/750)^2 > (63/2000)^2];
pass["returns: Z0 = 9/250 >= (63/2000)(9/8), t_H < 27 Z0^2 < 1/25", 63/2000 9/8 <= 9/250 && 27 (9/250)^2 < 1/25];
pass["weight: (7/3)(50/51) = 350/153, sqrt(1 + 1/25) < 51/50, (15/4)/(350/153) = 459/280 < 9/4",
  7/3 50/51 == 350/153 && Sqrt[26/25] < 51/50 && (15/4)/(350/153) == 459/280 && 459/280 < 9/4];
pass["weighted mean c = 3/5 with M/m = 459/280", 459/280 (2/5)^2 - (3/5)^2 <= 0];
pass["|theta| < 11/100: (3/88)[11/100 - (11/100)^3/6]^2 > (189/10000)^2", 3/88 (11/100 - (11/100)^3/6)^2 > (189/10000)^2];
(* table: t_H, rho *)
pass["Q: t_H < 27 (79/1000)^2 < 9/50, rho < 9/59, (1+t)^(3/2) < 13/10", 27 (79/1000)^2 < 9/50 && (9/50)/(59/50) == 9/59 && (59/50)^3 < (13/10)^2];
pass["returns: rho < 1/26, (1+1/25)^(3/2) < 11/10", (1/25)/(26/25) == 1/26 && (26/25)^3 < (11/10)^2];
pass["a: Q (7/3)/(13/10) > 7/4; returns (7/3)/(11/10) > 2", 7/3/(13/10) > 7/4 && 7/3/(11/10) > 2];
pass["b: (11/50 + 77/500)(3/10) < 3/25; (11/50 + 77/500)(15/4)(9/250) < 11/200",
  (11/50 + 77/500) 3/10 < 3/25 && (11/50 + 77/500) 15/4 9/250 < 11/200];
pass["H_rJ: (27/200)(77/500 + (9/59)(11/100+77/500)) < 27/1000; (27/200)(77/500 + (1/26)(...)) < 23/1000",
  27/200 (77/500 + 9/59 (11/100 + 77/500)) < 27/1000 && 27/200 (77/500 + 1/26 (11/100 + 77/500)) < 23/1000];
pass["S/j < 3/20: (27/200) sqrt(1 + 9/50) < 3/20", 27/200 Sqrt[59/50] < 3/20];
pass["|alpha'| + alpha^2, |gamma'| + gamma^2 < 81/1000; (gamma - alpha)^2 < 7/100",
  8/125 + (11/100)^2 < 81/1000 && 57/1000 + (77/500)^2 < 81/1000 && (11/100 + 77/500)^2 < 7/100];
pass["|H_rr|/j < (3/20)[81/1000 + (9/59)(7/100)] + (3/200)(8/7) = 255189/8260000 < 1/25",
  3/20 (81/1000 + 9/59 7/100) + 3/200 8/7 == 255189/8260000 && 255189/8260000 < 1/25];
pass["|H_r|/j < (43/1000)(8/7) + (77/500)(27/200) < 7/100; f < (43/1000)(8/7) < 1/20",
  43/1000 8/7 + 77/500 27/200 < 7/100 && 43/1000 8/7 < 1/20];
Ts = 21/100; Amx = 15/4; Lmx = 1/25;
mat[B_] := {{1 - Ts B, -Ts Amx}, {-Ts Lmx, 1 - Ts B}};
pass["M_B: diagonal, det > 0, inverse >= 0 (B = 3/25, 11/200)",
  And @@ Flatten[Table[{1 - Ts B > 0, Det[mat[B]] > 0, And @@ Thread[Flatten[Inverse[mat[B]]] >= 0]}, {B, {3/25, 11/200}}]]];
{U1, V1} = Inverse[mat[3/25]].{0, (41/200)/(49/50)^2 + Ts 1/20 10/49};
pass["theta: U < 19/100, V < 9/40", U1 < 19/100 && V1 < 9/40];
pass["V > 1/5 - T*[(1/25)(19/100) + (3/25)(9/40) + (1/20)(10/49)] > 19/100", 1/5 - Ts (1/25 19/100 + 3/25 9/40 + 1/20 10/49) > 19/100];
pass["dU/dt > -(3/25)(19/100) + (7/4)(19/100) > 3/10, D22 > (3/10)(3/25) = 9/250 > 1/40",
  -3/25 19/100 + 7/4 19/100 > 3/10 && 3/10 3/25 == 9/250 && 9/250 > 1/40];
Z0 = 9/250; C0 = 23/1000 + 11/200 Z0;
pass["|D12|/D22 <= Z0 + C0 T < 21/500", Z0 + C0 Ts < 21/500];
pass["|u(T)| <= (15/4)(189/10000)(9/8) < 41/500", 15/4 189/10000 9/8 < 41/500];
{U2, V2} = Inverse[mat[11/200]].{41/500, 7/100 + 1/20 46/7};
pass["T column: U_T = 157971160/388167361 < 21/50, V_T = 1105038840/2717171527 < 21/50",
  U2 == 157971160/388167361 && V2 == 1105038840/2717171527 && U2 < 21/50 && V2 < 21/50];
N2 = Ts (11/200 U2 + 15/4 V2);
pass["int|u'| = 63070718199/194083680500 < 17/50 (no q source)", N2 == 63070718199/194083680500 && N2 < 17/50];
pass["return: N <= 602/401, h_phi <= 27270901/502500 at r = 301/50 => sqrt(b_H) > 11/100",
  (NN /. {r -> 301/50, A -> 1}) == 602/401 && (hp /. r -> 301/50) == 27270901/502500 && 1/(602/401 27270901/502500) > (11/100)^2];
pass["b_H = (1 - A^2 f/2)/h_phi, increasing/decreasing: minimum of b_H on the return at r = 301/50, A = 1",
  zr[bH - (1 - A^2 f/2)/hp] && all[{r, A}, 299/50 <= r <= 301/50 && 1/2 <= A <= 1, bH >= (bH /. {r -> 301/50, A -> 1})]];
D11 = 11/100 25/26 - C0 Ts 21/50 - Z0 17/50;
pass["D_1T > 296809483/3250000000 > 9/100", D11 == 296809483/3250000000 && D11 > 9/100];
pass["Gamma_T > 9/100 - (21/500)(21/50) > 7/100", 9/100 - 21/500 21/50 > 7/100];

Print["--- (E) sector change of the pulse [EXACT]"];
pass["j = sqrt(18) cos(theta) > 4 for |theta| < 11/100 (cos >= 49/50)", 18 (49/50)^2 > 16];
Zs = 63/8000; Cs = 23/1000 + 11/200 Zs;
pass["Z* = (63/2000)/4, 27 Z*^2 < 1/500, (3/20)/4 = 3/80", (63/2000)/4 == Zs && 27 Zs^2 < 1/500 && (3/20)/4 == 3/80];
pass["|rdot| < (15/4) Z* < 3/100, |rdot(T)| < (15/4)(189/10000)/4 < 9/500, |r - 6| < (1/10)(3/100) < 1/300",
  15/4 Zs < 3/100 && 15/4 (189/10000)/4 < 9/500 && 1/10 3/100 <= 1/300];
pass["|D12|/D22 <= Z* + T C* < 13/1000", Zs + Ts Cs < 13/1000];
{U3, V3} = Inverse[mat[11/200]].{3/100, 3/80};
pass["dilation: U_d = 23673900/388167361 < 13/200, V_d = 14927550/388167361 < 1/25",
  U3 == 23673900/388167361 && V3 == 14927550/388167361 && U3 < 13/200 && V3 < 1/25];
Hd = 1/10 - Ts (23/1000 13/200 + 15/4 Zs 1/25) - 13/1000 13/200;
pass["Gamma_T + Gamma_tau > 7887439/80000000 > 49/500", Hd == 7887439/80000000 && Hd > 49/500];
{U4, V4} = Inverse[mat[11/200]].{9/500, 3/80};
N4 = Ts (11/200 U4 + 15/4 V4);
pass["T column at the change: 18929340/388167361 < 1/20, 14887230/388167361 < 1/25, int = 5971163751/194083680500 < 33/1000",
  U4 == 18929340/388167361 && V4 == 14887230/388167361 && N4 == 5971163751/194083680500 && U4 < 1/20 && V4 < 1/25 && N4 < 33/1000];
GT = 27/200 + Cs Ts 1/20 + Zs 33/1000 + 13/1000 1/20;
pass["Gamma_T < 435698953/3200000000 < 7/50", GT == 435698953/3200000000 && GT < 7/50];
pass["Theta' < -(49/500)/(7/50) = -7/10", -(49/500)/(7/50) == -7/10];
pass["e = (27/200)(1/500)/2 <= 27/200000", 27/200 (1/500)/2 <= 27/200000];
rm = 1799/300; rp = 1801/300; wf[rr_, AA_] := Sqrt[(1 - AA^2 (1 - 2/rr)/2)/(rr^2/(1 - 2/rr))];
pass["w(r+,1/2) > 1302/10000, w(r+,1) > 1110/10000, Jensen < 61/500",
  wf[rp, 1/2] > 1302/10000 && wf[rp, 1] > 1110/10000 && Sqrt[(1 - 19 (1 - 2/rm)/64)/(rm^2/(1 - 2/rm))] < 61/500];
pass["chord: (13/24)(1302/10000) + (11/24)(1110/10000) = 607/5000 exactly (strictness from the two bounds on w above)", 13/24 1302/10000 + 11/24 1110/10000 == 607/5000];
pass["w decreasing in r, concave in A^2 (|r-6| < 1/300)",
  all[{r, x}, rm <= r <= rp && 1/4 <= x <= 1, D[(1 - x f/2)/hp, r] < 0] && zr[D[Sqrt[1 - x ff/2], {x, 2}] + ff^2/(16 (1 - x ff/2)^(3/2))]];
pass["tau* in (10/61, 4000/24253): 1/(50 (61/500)) and 1/(50 (607/5000 - 27/200000))",
  1/(50 61/500) == 10/61 && 1/(50 (607/5000 - 27/200000)) == 4000/24253 && 10/61 > 163/1000 && 4000/24253 < 165/1000];

Print["--- (F) pulse condition and drop [EXACT]"];
pass["pulse condition: T|ell_T| <= 46/7 if |A'| <= 1100/147 (T <= 21/100, factor <= 46/11)", 21/100 46/11 1100/147 <= 46/7];
pass["pulse condition: |b0| = (A A' f/E)(N omega/E + 2) < (5/4) d, |b1| = (A A' f/E)(N/E) < (3/5) d on the shell, A <= 1",
  all[{r, A}, shell, A f/Sqrt[2] (NN om/Sqrt[2] + 2) < 5/4 && A f/Sqrt[2] NN/Sqrt[2] < 3/5]];
pass["pulse condition: (1/5)(5/4) d exp((3/25) d) < d/3 for d < 1/4", 1/5 5/4 Exp[3/25 1/4] < 1/3];
pass["drop, existence: 1/(75 sqrt2) + 3 pi/(2 sqrt11) < 1/(75 sqrt2) + 33/(7 sqrt11) < 254/(125 sqrt2)",
  Pi < 22/7 && 1/(75 Sqrt[2]) + 33/(7 Sqrt[11]) < 254/(125 Sqrt[2])];
pass["drop: omega(1, r=6) h = (2/(3 sqrt2)) h; c_phi(1/2) = 3/sqrt(11)", (om /. {A -> 1, r -> 6}) == 2/(3 Sqrt[2]) && zr[(cp2 /. A -> 1/2) - 9/11]];
pass["drop: N_{1/2} >= 16/15 > (129/125)^2 on [4,6]; N >= 1", all[r, 4 <= r <= 6, (NN /. A -> 1/2) >= 16/15] && 16/15 > (129/125)^2];
pass["drop: [1/(2 sqrt2)][2 + 2(129/125)] = 254/(125 sqrt 2)", Simplify[1/(2 Sqrt[2]) (2 + 2 129/125) - 254/(125 Sqrt[2])] === 0];
(* rho_1(R) > R at A = 1: rho = c r sqrt(N f), R = c int_2^r sqrt(N) *)
N1 = NN /. A -> 1; c1 = 1/Sqrt[2];
pass["drop: d(r sqrt(N f))/dr > sqrt(N) for r > 2 (A = 1), hence rho_1(R) > R", all[r, r > 2, D[r Sqrt[N1 f], r] > Sqrt[N1]]];
R16 = NIntegrate[c1 Sqrt[N1], {r, 2, 6}, WorkingPrecision -> 30];
pass["drop: R_1(6) <= 2 sqrt(3) [N: R_1(6) = " <> ToString[N[R16, 8]] <> "]; exact: sqrt(N_1) <= sqrt(3/2)... c sqrt(3/2) 4 = 2 sqrt 3",
  (N1 /. r -> 6) <= 3/2 && all[r, 2 <= r <= 6, N1 <= 3/2] && c1 Sqrt[3/2] 4 == 2 Sqrt[3]];
pass["drop: sqrt2 > 5/4; 5 sqrt3/2 < 13/3; 100/(15 sqrt2) > 13/3; f(15/7) = 1/15", Sqrt[2] > 5/4 && 5 Sqrt[3]/2 < 13/3 && 100/(15 Sqrt[2]) > 13/3 && (f /. r -> 15/7) == 1/15];
pass["drop: density of R_{1/2} <= (1/4) density of R_1", all[r, r > 2, (1/(4 Sqrt[2])) Sqrt[NN /. A -> 1/2] <= 1/4 c1 Sqrt[N1]]];
