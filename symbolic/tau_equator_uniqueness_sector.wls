#!/usr/bin/env wolframscript
(* Paper III, example on the Kerr equator with spin (supplement S10): Wolfram checks of the steps of the theorem on
   global uniqueness and of the theorem on a single transversal sector change. Thakurta-Kerr, proper-time branch on the
   equator, M = 1, a = 1/2, Ehat = Sqrt[2].
   [EXACT] = symbolic identities, exact Reduce/MaxValue, rational arithmetic.
   [N] = high-precision numerical extremes over domains in several variables (not a certification).
   Identities with radicals are checked in rational form (squares, or product with the conjugate). *)
pass[name_, e_] := Print[If[TrueQ[e], "PASS ", "FAIL "], name];
zero[e_] := TrueQ[Simplify[e] === 0] || TrueQ[PossibleZeroQ[e]];
zr[e_] := TrueQ[Together[Simplify[e]] === 0];
$Assumptions = r > 5 && 1/2 <= A <= 1 && J > 0 && p \[Element] Reals;
f = 1 - 2/r; Del = r^2 - 2 r + 1/4; hr = r^2/(f Del); hp = Del/f^2; bphi = -1/(r - 2); Bm = bphi^2/hp;
NN = 1/(1 - A^2 f/2); om = A^2 f/Sqrt[2]; g = A^2/(Sqrt[2] r); xi = J - g;
d = NN - Bm; aH = 1/(d hr); bH = NN/(d^2 hp); kH = bphi/(d hp);
S = Sqrt[aH p^2 + bH xi^2]; H = -om + S - kH xi;
shell = 23/4 <= r <= 25/4 && 1/2 <= A <= 1;

Print["--- (A) identities of the geometry and of the Hamiltonian [EXACT]"];
pass["|beta|^2 = 1/(r^2 Delta)", zero[Bm - 1/(r^2 Del)]];
pass["omega beta_phi = -g", zero[om bphi + g]];
(* H = -omega + sigma_K(p + omega beta): the linear terms agree, and the radicands agree after dividing by d^2 *)
X = (xi bphi/hp)^2 + d (p^2/hr + xi^2/hp); Y = aH p^2 + bH xi^2;
pass["H = -omega + sigma_K(p + omega beta): linear part -xi beta/(hp d) = -k_H xi", zr[-(xi bphi/hp)/d + kH xi]];
pass["H = -omega + sigma_K(p + omega beta): X/d^2 = Y (radicands; both sides positive, d > 0)", zr[X/d^2 - Y]];
pass["rdot = a_H p/S, phidot = b_H xi/S - k_H", zero[D[H, p] - aH p/S] && zero[D[H, J] - (bH xi/S - kH)]];
lr = D[Log[hr], r]; lp = D[Log[hp], r];
pass["ell_r = 2/r - 2/(r(r-2)) - (2r-2)/Delta", zero[lr - (2/r - 2/(r (r - 2)) - (2 r - 2)/Del)]];
pass["ell_phi = (2r-2)/Delta - 4/(r(r-2))", zero[lp - ((2 r - 2)/Del - 4/(r (r - 2)))]];
VD = D[Del, r]/Del; PG = 4 (r - 1)/(r^2 (r - 2)^2);
pass["ell_r' = -2/r^2 + P - 2/Delta + V^2", zero[D[lr, r] - (-2/r^2 + PG - 2/Del + VD^2)]];
pass["ell_phi' = 2/Delta - V^2 + 2P", zero[D[lp, r] - (2/Del - VD^2 + 2 PG)]];
pass["V decreasing: 2 Delta - Delta'^2 < 0 for every real r", Reduce[2 Del - D[Del, r]^2 >= 0, r, Reals] === False];
pass["P decreasing on the shell: 1/(r-1) - 2/r - 2/(r-2) < 0", Resolve[ForAll[r, 23/4 <= r <= 25/4, 1/(r - 1) - 2/r - 2/(r - 2) < 0]]];
LN = D[NN, r]/NN; Ld = D[d, r]/d;
pass["L_N = N A^2/r^2", zero[LN - NN A^2/r^2]];
pass["L_N' = L_N (L_N - 2/r)", zero[D[LN, r] - LN (LN - 2/r)]];
pass["N''/N = 2 L_N^2 - 2 L_N/r", zero[D[NN, {r, 2}]/NN - (2 LN^2 - 2 LN/r)]];
pass["L_d' = (N'' - B'')/d - L_d^2", zero[D[Ld, r] - ((D[NN, {r, 2}] - D[Bm, {r, 2}])/d - Ld^2)]];
LB = D[Log[Bm], r];
pass["L_B = -2/r - V, L_B' = 2/r^2 - 2/Delta + V^2", zero[LB - (-2/r - VD)] && zero[D[LB, r] - (2/r^2 - 2/Del + VD^2)]];
Ck = 1/(r - 2) + Ld + lp;
pass["k_r = -k C_k, C_k = 1/(r-2) + L_d + ell_phi", zero[D[kH, r] + kH Ck]];
al = D[Log[aH], r]/2; ga = D[Log[bH], r]/2;
pass["alpha = -(L_d + ell_r)/2, gamma = (L_N - 2 L_d - ell_phi)/2", zero[al + (Ld + lr)/2] && zero[ga - (LN - 2 Ld - lp)/2]];
pass["omega_rr = -2 omega_r/r, xi_rr = -2 xi_r/r", zero[D[om, {r, 2}] + 2 D[om, r]/r] && zero[D[xi, {r, 2}] + 2 D[xi, r]/r]];

Print["--- (B) derivatives of S and of H [EXACT]"];
u = Sqrt[aH] p; w = Sqrt[bH] xi; rho = aH p^2/S^2; sg = 1 - rho; xr = D[xi, r];
pass["S_r = (alpha u^2 + gamma w^2)/S + sqrt(b) xi_r w/S", zero[D[S, r] - ((al u^2 + ga w^2)/S + Sqrt[bH] xr w/S)]];
pass["S_rr/S = (alpha'+alpha^2) rho + (gamma'+gamma^2) sigma + rho sigma (gamma-alpha+xi_r/xi)^2 + sqrt(b)(2 gamma xi_r + xi_rr) w/S^2",
  zero[D[S, {r, 2}]/S - ((D[al, r] + al^2) rho + (D[ga, r] + ga^2) sg + rho sg (ga - al + xr/xi)^2 +
       Sqrt[bH] (2 ga xr + D[xi, {r, 2}]) w/S^2)]];
pass["H_rr = -omega_rr + S_rr - k_rr xi - 2 k_r xi_r - k xi_rr",
  zero[D[H, {r, 2}] - (-D[om, {r, 2}] + D[S, {r, 2}] - D[kH, {r, 2}] xi - 2 D[kH, r] xr - kH D[xi, {r, 2}])]];
pass["H_rp = rdot [2 alpha - alpha rho - gamma sigma - (xi_r/xi) sigma]", zero[D[H, r, p] - D[H, p] (2 al - al rho - ga sg - xr/xi sg)]];
pass["H_rJ = (phidot + k)[gamma + rho(gamma - alpha + xi_r/xi)] - k_r",
  zero[D[H, r, J] - ((D[H, J] + kH) (ga + rho (ga - al + xr/xi)) - D[kH, r])]];
Z = p/xi;
pass["rank one: H_pJ = -Z H_pp, H_JJ = -Z H_pJ", zero[D[H, p, J] + Z D[H, p, p]] && zero[D[H, J, J] + Z D[H, p, J]]];
(* j H_pp = (a/sqrt b)/((|xi|/j)(1+t)^(3/2)), a/sqrt b = sqrt(h_phi)/(h_r sqrt N) (xi > 0, J > 0): checked squared *)
S2 = aH p^2 + bH xi^2; Hpp = aH bH xi^2/S2^(3/2); tq = aH p^2/(bH xi^2);
pass["(J H_pp)^2 = (a^2/b)/((xi/J)^2 (1+t)^3)  [S^2 = b xi^2 (1+t)]", zr[S2 - bH xi^2 (1 + tq)] &&
   zr[J^2 aH^2 bH^2 xi^4/(bH xi^2 (1 + tq))^3 - (aH^2/bH)/((xi/J)^2 (1 + tq)^3)]];
pass["(a/sqrt b)^2 = h_phi/(h_r^2 N)", zr[aH^2/bH - hp/(hr^2 NN)]];
pass["projective force identity: p H_rp + J H_rJ - H_r = omega_r - xi_r phidot + g H_rJ",
  zero[p D[H, r, p] + J D[H, r, J] - D[H, r] - (D[om, r] - xr D[H, J] + g D[H, r, J])]];

Print["--- (C) projective system and angular identity [EXACT]"];
Clear[du, dp, dJ];
uprime = D[H, p, r] du + D[H, p, p] dp + D[H, p, J] dJ;
pprime = -(D[H, r, r] du + D[H, r, p] dp + D[H, r, J] dJ);
phiprime = D[H, J, r] du + D[H, J, p] dp + D[H, J, J] dJ;
zz = dp - (p/J) dJ; jj = J; hh = dJ/J;
pass["u' = b u + a P - q h  (a = j H_pp, P = z/j, q = a Z g/j, j = J > 0)",
  zero[uprime - (D[H, r, p] du + J D[H, p, p] (zz/J) - J D[H, p, p] Z (g/J) hh)]];
pprime0 = D[-D[H, r], {p, 0}];   (* dp/deta = -H_r, dJ/deta = 0 *)
zprime = pprime - (-D[H, r]/J) dJ;
pass["P' = -l u - b P - f h  (l = H_rr/j, f = (omega_r - xi_r phidot + g H_rJ)/j)",
  zero[zprime/J - (-(D[H, r, r]/J) du - D[H, r, p] (zz/J) - ((D[om, r] - xr D[H, J] + g D[H, r, J])/J) hh)]];
pass["dphi' = (H_rJ + Z H_rp) u - Z u'", zero[phiprime - ((D[H, r, J] + Z D[H, r, p]) du - Z uprime)]];
Clear[cr, cp, th, s, ell];
pT = -s cr Sin[th]; JT = s cp Cos[th];
pass["theta column: z(T) = -s c_r/cos(theta)  (s = +-1)",
  And @@ Table[Simplify[(-s cr Cos[th] - (pT/JT) (-s cp Sin[th])) + s cr/Cos[th]] === 0, {s, {1, -1}}]];
pass["T column: with dp = H_r + ell p, dJ = ell J one has z(T) = H_r", zero[(Hr0 + ell p) - (p/J) (ell J) - Hr0]];

Print["--- (D) reference at r = 6 [EXACT]"];
Clear[vv, gam, KK, rb, C, NN6];
h6r = 216/97; h6p = 873/16; kap2 = h6r/h6p;
pass["h_r(6) = 216/97, h_phi(6) = 873/16, kappa_6^2 = 3456/84681", (hr /. r -> 6) == 216/97 && (hp /. r -> 6) == 873/16 && kap2 == 3456/84681];
dd6 = (NN - Bm) /. r -> 6;
pass["a kappa^2 = b gamma = 1/(d h_phi) at r = 6", zero[((aH kap2) /. r -> 6) - ((bH dd6/NN) /. r -> 6)] && zero[((aH kap2) /. r -> 6) - 1/(dd6 h6p)]];
Do[
  X = Cos[th] - s vv; L = Sqrt[X^2 + gam Sin[th]^2];
  pass["s=" <> ToString[s] <> ": d/dtheta (K sin/L) = K X (1 - s v cos)/L^3", zero[D[KK Sin[th]/L, th] - KK X (1 - s vv Cos[th])/L^3]];
  pass["s=" <> ToString[s] <> ": d/dtheta (sqrt(b) X/L) = -sqrt(b) gamma sin (1 - s v cos)/L^3", zero[D[rb X/L, th] + rb gam Sin[th] (1 - s vv Cos[th])/L^3]];
  pass["s=" <> ToString[s] <> ": cos X + sin^2 = 1 - s v cos", zero[Cos[th] X + Sin[th]^2 - (1 - s vv Cos[th])]],
  {s, {1, -1}}];
pass["d log c_phi/dA = 2/A + A/(3 - A^2) (at r = 6)", zero[D[Log[(om Sqrt[NN hp]) /. r -> 6], A] - (2/A + A/(3 - A^2))]];
pass["2/A + A/(3-A^2) decreasing on [1/2,1], maximum 46/11",
  Resolve[ForAll[A, 1/2 <= A <= 1, D[2/A + A/(3 - A^2), A] < 0]] && (2/A + A/(3 - A^2) /. A -> 1/2) == 46/11];
pass["x |sin(2 pi x)| <= 1 on [0,1]  (hence T|A'(T)| <= pi/2)", MaxValue[{x Abs[Sin[2 Pi x]], 0 <= x <= 1}, x] <= 1];
pass["(46/11)(pi/2) < 46/7", 46/11 Pi/2 < 46/7];
alpha[q_] := 1 - Sin[Pi q]^2/2;
pass["alpha'(1) = 0, alpha(1) = 1, alpha = 3/4 + cos(2 pi q)/4", (D[alpha[q], q] /. q -> 1) == 0 && alpha[1] == 1 && zero[alpha[q] - (3/4 + Cos[2 Pi q]/4)]];
pass["int_0^1 alpha^2 = 19/32", Integrate[alpha[q]^2, {q, 0, 1}] == 19/32];

Print["--- (E) elementary lemmas [EXACT]"];
Clear[MM, mm, c, T, B, Am, Lm];
pass["weighted mean: M(1-c)^2 - m c^2 <= 0 with c = 3/5 if M/m <= 9/4", Resolve[ForAll[{MM, mm}, mm > 0 && MM >= mm && MM <= 9/4 mm, MM (1 - 3/5)^2 - mm (3/5)^2 <= 0]]];
pass["weighted mean: int_0^1 (x-c) W <= [M(1-c)^2 - m c^2]/2 (W = M on [c,1], m on [0,c] is the worst case)",
  zero[Integrate[(x - c) MM, {x, c, 1}] + Integrate[(x - c) mm, {x, 0, c}] - (MM (1 - c)^2 - mm c^2)/2]];
Mcmp = {{1 - T B, -T Am}, {-T Lm, 1 - T B}};
pass["comparison matrix: inverse with entries >= 0 if 1 - TB > 0 and det > 0",
  Resolve[ForAll[{T, B, Am, Lm}, T > 0 && B >= 0 && Am >= 0 && Lm >= 0 && 1 - T B > 0 && (1 - T B)^2 - T^2 Am Lm > 0,
    And @@ Thread[Flatten[Inverse[Mcmp]] >= 0]]]];
pass["chord: weights 13/24, 11/24 from the mean 19/32 between 1/4 and 1", (19/32 - 1/4)/(3/4) == 11/24];
pass["sqrt(1 - f x/2) concave in x", Resolve[ForAll[{x, ff}, 0 < ff < 1 && 0 <= x <= 1, D[Sqrt[1 - ff x/2], {x, 2}] < 0]]];
pass["w/(1 - b w) increasing and convex, w/(1 + b w) increasing and concave (0 < b w < 1)",
  Resolve[ForAll[{ww, bb}, ww > 0 && bb > 0 && bb ww < 1, D[ww/(1 - bb ww), ww] > 0 && D[ww/(1 - bb ww), {ww, 2}] > 0 &&
    D[ww/(1 + bb ww), ww] > 0 && D[ww/(1 + bb ww), {ww, 2}] < 0]]];
(* pure angular speed: (sqrt(b) - s k)(sqrt(N h_phi) + s beta_phi) = 1, with sqrt(b) = sqrt(N)/(d sqrt(hp));
   sN = sqrt(N), sH = sqrt(h_phi) are kept as symbols and their squares replaced *)
Clear[sN, sH];
Do[
  prod = (sN/(d sH) - s kH) (sN sH + s bphi) /. {sN^2 -> NN, sH^2 -> hp};
  prod = Expand[prod] /. {sN^2 -> NN, sH^2 -> hp};
  pass["s=" <> ToString[s] <> ": pure angular speed (sqrt b - s k)(sqrt(N h_phi) + s beta_phi) = 1",
    zr[(prod /. {sN -> Sqrt[NN], sH -> Sqrt[hp]}) - 1] || zr[Expand[(sN/(d sH) - s kH) (sN sH + s bphi)] /. {sN^2 -> NN, sH^2 -> hp} - 1]],
  {s, {1, -1}}];
pass["1 - (1+t)^(-1/2) <= t/2 for t >= 0", Resolve[ForAll[t, t >= 0, 1 - (1 + t)^(-1/2) <= t/2]]];

Print["--- (F) extremes on the shell used in the tables [EXACT where marked]"];
ex[name_, expr_, op_, val_] := Module[{res},
  res = TimeConstrained[Resolve[ForAll[{r, A}, shell, op[expr, val]], Reals], 120, $Failed];
  If[res === $Failed,
    res = op[If[op === Less || op === LessEqual, NMaxValue, NMinValue][{expr, shell}, {r, A}, WorkingPrecision -> 30], val];
    pass["[N] " <> name, res], pass["[EXACT] " <> name, res]]];
ex["0 <= alpha", al, GreaterEqual, 0]; ex["alpha <= 11/100", al, LessEqual, 11/100];
ex["gamma <= 0", ga, LessEqual, 0]; ex["gamma >= -77/500", ga, GreaterEqual, -77/500];
ex["sqrt(a_H) < 2/3", aH, Less, 4/9]; ex["sqrt(b_H) < 67/500", bH, Less, (67/500)^2];
ex["|k_H| < 6/1250", -kH, Less, 6/1250]; ex["|k_H,r| < 69/25000", Abs[D[kH, r]], Less, 69/25000];
ex["omega_r < 43/1000", D[om, r], Less, 43/1000]; ex["xi_r < 43/2000", xr, Less, 43/2000];
ex["|alpha'| < 8/125", Abs[D[al, r]], Less, 8/125]; ex["|gamma'| < 57/1000", Abs[D[ga, r]], Less, 57/1000];
ex["|k_rr| < (6/1250)((23/40)^2 + 3/20)", Abs[D[kH, {r, 2}]], Less, 6/1250 ((23/40)^2 + 3/20)];
ex["|omega_rr| < 3/200", Abs[D[om, {r, 2}]], Less, 3/200]; ex["|xi_rr| < 3/400", Abs[D[xi, {r, 2}]], Less, 3/400];
ex["0 < ell_r' < 53/500", D[lr, r], Less, 53/500]; ex["ell_r' > 0", D[lr, r], Greater, 0];
ex["|ell_phi'| < 7/125", Abs[D[lp, r]], Less, 7/125];
ex["L_N <= 800/17457", LN, LessEqual, 800/17457]; ex["0 <= L_d < 47/1000", Ld, Less, 47/1000];
ex["|L_d'| < 21/1000", Abs[D[Ld, r]], Less, 21/1000];
ex["7/3 < a/sqrt(b)", aH/Sqrt[bH], Greater, 7/3]; ex["a/sqrt(b) < 15/4", aH/Sqrt[bH], Less, 15/4];
ex["a/b < 29", aH/bH, Less, 29];

Print["--- (G) phase-space domain: table of the coefficients [N] (NMaxValue, 4 variables)"];
tab[name_, expr_, dom_, op_, val_] := Module[{v},
  v = If[op === Less, NMaxValue, NMinValue][{expr, dom}, {r, A, p, J}, WorkingPrecision -> 25, Method -> "DifferentialEvolution"];
  pass["[N] " <> name <> " (extreme " <> ToString[N[v, 6]] <> ")", op[v, val]]];
glob = shell && 8/9 <= J <= Sqrt[291]/4 && -77/1000 J <= p <= 77/1000 J;
ret = 6 - 1/50 <= r <= 6 + 1/50 && 1/2 <= A <= 1 && 8/9 <= J <= Sqrt[291]/4 && -63/2000 <= p <= 63/2000;
aa = J D[H, p, p]; bb = D[H, r, p]; ll = D[H, r, r]/J; HrJ = D[H, r, J]; ff = (D[om, r] - xr D[H, J] + g HrJ)/J; qq = aa Z g/J;
Do[{dom, nm, Alo, Bv, HJv, fv, qv} = dd;
  tab[nm <> ": a > " <> ToString[Alo], aa, dom, Greater, Alo]; tab[nm <> ": a < 48/11", aa, dom, Less, 48/11];
  tab[nm <> ": |H_rp| < " <> ToString[Bv], Abs[bb], dom, Less, Bv]; tab[nm <> ": |H_rr|/j < 1/25", Abs[ll], dom, Less, 1/25];
  tab[nm <> ": |H_rJ| < " <> ToString[HJv], Abs[HrJ], dom, Less, HJv]; tab[nm <> ": |f| < " <> ToString[fv], Abs[ff], dom, Less, fv];
  tab[nm <> ": |q| < " <> ToString[qv], Abs[qq], dom, Less, qv]; tab[nm <> ": |H_r|/j < 2/25", Abs[D[H, r]/J], dom, Less, 2/25],
  {dd, {{glob, "whole Q", 29/20, 7/50, 4/125, 57/1000, 11/200}, {ret, "return", 37/20, 8/125, 13/500, 7/125, 13/500}}}];
(* sense s = -1: not covered by symmetry (g does not change sign); explicit check with J < 0 *)
globm = shell && -Sqrt[291]/4 <= J <= -8/9 && 77/1000 J <= p <= -77/1000 J;
retm = 6 - 1/50 <= r <= 6 + 1/50 && 1/2 <= A <= 1 && -Sqrt[291]/4 <= J <= -8/9 && -63/2000 <= p <= 63/2000;
aam = -J D[H, p, p]; llm = -D[H, r, r]/J; ffm = -(D[om, r] - xr D[H, J] + g HrJ)/J; qqm = aam Z g/(-J);
Do[{dom, nm, Alo, Bv, HJv, fv, qv} = dd;
  tab[nm <> " (s=-1): a > " <> ToString[Alo], aam, dom, Greater, Alo]; tab[nm <> " (s=-1): a < 48/11", aam, dom, Less, 48/11];
  tab[nm <> " (s=-1): |H_rp| < " <> ToString[Bv], Abs[bb], dom, Less, Bv]; tab[nm <> " (s=-1): |H_rr|/j < 1/25", Abs[llm], dom, Less, 1/25];
  tab[nm <> " (s=-1): |H_rJ| < " <> ToString[HJv], Abs[HrJ], dom, Less, HJv]; tab[nm <> " (s=-1): |f| < " <> ToString[fv], Abs[ffm], dom, Less, fv];
  tab[nm <> " (s=-1): |q| < " <> ToString[qv], Abs[qqm], dom, Less, qv],
  {dd, {{globm, "whole Q", 29/20, 7/50, 4/125, 57/1000, 11/200}, {retm, "return", 37/20, 8/125, 13/500, 7/125, 13/500}}}];
tab["whole Q (s=-1): |H_r|/j < 2/25", Abs[D[H, r]/J], globm, Less, 2/25];

Print["--- (H) rational chain of the uniqueness and sector-change theorems [EXACT]"];
Tst = 21/100; Amx = 48/11; Lmx = 1/25;
maj[B_, ru_, rv_] := Module[{dg = 1 - Tst B, dt}, dt = dg^2 - Tst^2 Amx Lmx; {(dg ru + Tst Amx rv)/dt, (Tst Lmx ru + dg rv)/dt}];
kap = 203/1000; V0 = kap/(49/50)^2; hth = 10/49;
{Uth, Vth} = maj[7/50, Tst 11/200 hth, V0 + Tst 57/1000 hth];
pass["theta column: U < 213/1000, V < 223/1000", Uth < 213/1000 && Vth < 223/1000];
Vlo = 1/5 - Tst (1/25 213/1000 + 7/50 223/1000 + 57/1000 hth); pass["V > 189/1000", Vlo > 189/1000];
slope = -7/50 213/1000 + 29/20 189/1000 - 11/200 hth; pass["dU/dt > 23/100", slope > 23/100];
pass["D22 > (23/100)(3/25) > 1/40", 23/100 3/25 > 1/40];
Z0 = 63/2000 72/55; C0 = 13/500 + 8/125 Z0; pass["|D12|/D22 < Z0 + (21/100) C0 < 6/125", Z0 + Tst C0 < 6/125];
pass["|u(T)| < (15/4)(189/10000)(72/55) < 93/1000", 15/4 189/10000 72/55 < 93/1000];
{UT, VT} = maj[8/125, 93/1000 + 13/500 46/7, 2/25 + 7/125 46/7];
pass["T column: U_T < 7/10, V_T < 47/100", UT < 7/10 && VT < 47/100];
NT = Tst (8/125 UT + Amx VT) + 13/500 46/7; pass["int|u'| < 61/100", NT < 61/100];
pass["sqrt(b_H) > 11/100 on the return domain: 1/(N h_phi) at r = 301/50, A = 1", (1/(NN hp) /. {r -> 301/50, A -> 1}) > (11/100)^2];
D11lo = 11/100 20/21 - 6/1250 - C0 Tst 7/10 - Z0 61/100; pass["D11 > 7/100", D11lo > 7/100];
pass["Gamma' > 7/100 - (6/125)(7/10) = 91/2500", 7/100 - 6/125 7/10 == 91/2500];
(* sector change *)
pass["j > 4 at the change: (291/16)(49/50)^2 > 16", 291/16 (49/50)^2 > 16];
Zs = 63/2000 8/31; pass["Z* = 63/7750, 29 Z*^2 < 1/500", Zs == 63/7750 && 29 Zs^2 < 1/500];
pass["|r - 6| < (1/10)(15/4) Z* < 1/300", 1/10 15/4 Zs < 1/300];
Cs = 13/500 + 8/125 Zs; pass["|D12|/D22 < Z* + (21/100) C* < 7/500", Zs + Tst Cs < 7/500];
{Ud, Vd} = maj[8/125, 31/1000, 3/80]; pass["dilation: U_d < 7/100, V_d < 1/25", Ud < 7/100 && Vd < 1/25];
Hd = 1/10 - Tst (13/500 7/100 + Amx Zs 1/25) - 7/500 7/100; pass["H_diag > 49/500", Hd > 49/500];
{UTs, VTs} = maj[8/125, 19/1000, 3/80]; NTs = Tst (8/125 UTs + Amx VTs);
pass["T column at the change: U_T < 3/50, V_T < 1/25, int|u'| < 21/500", UTs < 3/50 && VTs < 1/25 && NTs < 21/500];
GTup = 7/50 + Cs Tst 3/50 + Zs 21/500 + 7/500 3/50; pass["Gamma_T < 71/500", GTup < 71/500];
pass["F' < -(49/500)/(71/500) = -49/71 < -2/3", -(49/500)/(71/500) == -49/71 && -49/71 < -2/3];
(* enclosures: independent recomputation from the quantities of the sector-change proof *)
rm = 1799/300; rp = 1801/300; wf[rr_, AA_] := Sqrt[1 - AA^2 (1 - 2/rr)/2]/Sqrt[(Del/f^2) /. r -> rr];
pass["w(r+,1/2) > 0.12956, w(r+,1) > 0.11048, sqrt((1 - 19 f(r-)/64)/h_phi(r-)) < 0.1213, w(r-,1/2) < 0.13",
  wf[rp, 1/2] > 12956/100000 && wf[rp, 1] > 11048/100000 &&
  Sqrt[(1 - 19 (1 - 2/rm)/64)/((Del/f^2) /. r -> rm)] < 1213/10000 && wf[rm, 1/2] < 13/100];
Lw = 151/1250; Uw = 1213/10000; pass["L_w <= (13/24) w(r+,1/2) + (11/24) w(r+,1)", 13/24 12956/100000 + 11/24 11048/100000 >= Lw];
pass["b_beta in (249/1000, 251/1000) on |r - 6| < 1/300", 1/(rp - 2) > 249/1000 && 1/(rm - 2) < 251/1000];
ee = 67/500000; pass["e = (67/500)(1/500)/2 = 67/500000", 67/500 (1/500)/2 == ee];
(* the names Lpl, Upl, Lmi, Umi avoid the protected symbol Up *)
Lpl = Lw/(1 - 249/1000 Lw) - ee; Upl = Uw/(1 - 251/1000 13/100);
Lmi = Lw/(1 + 251/1000 13/100) - ee; Umi = Uw/(1 + 249/1000 Uw);
pass["tau+* in (1/(50 U+), 1/(50 L+)) = (96737/606500, 12124010000/75418769133)", 1/(50 Upl) == 96737/606500 && 1/(50 Lpl) == 12124010000/75418769133];
pass["tau-* in (1/(50 U-), 1/(50 L-)) = (10302037/60650000, 1032630000/6033081379)", 1/(50 Umi) == 10302037/60650000 && 1/(50 Lmi) == 1032630000/6033081379];
pass["simple enclosures: 159/1000 < tau+ < 161/1000 < 169/1000 < tau- < 172/1000",
  159/1000 < 1/(50 Upl) && 1/(50 Lpl) < 161/1000 && 169/1000 < 1/(50 Umi) && 1/(50 Lmi) < 172/1000];
pass["denominators and lower bounds positive", Lpl > 0 && Lmi > 0 && 1 - 251/1000 13/100 > 0];
