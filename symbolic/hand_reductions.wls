#!/usr/bin/env wolframscript
(* Paper III: Wolfram confirmation of the reasons written by hand in the text in place of decisions by Reduce or
   NMinValue. Every check is an identity (Simplify to 0) or an inequality decided exactly; two checks must fail. *)
pass = 0; fail = 0;
ok[name_, test_] := If[TrueQ[test], pass++; Print["PASS  ", name], fail++; Print["FAIL  ", name]];
id[name_, expr_, ass_: True] := ok[name, Simplify[expr, ass] === 0];

(* 1. lem:radius: 86q^2 - 72q + 1, zeros (36 +- sqrt(1210))/86 < 1 *)
id["radius: zeros of 86q^2-72q+1", Expand[86 (q - (36 + Sqrt[1210])/86) (q - (36 - Sqrt[1210])/86) - (86 q^2 - 72 q + 1)]];
ok["radius: (36+sqrt(1210))/86 < 1", (36 + Sqrt[1210])/86 < 1];
ok["radius: positive on (1,5/4) (Reduce, control)", Reduce[1 < q < 5/4 && 86 q^2 - 72 q + 1 <= 0, q, Reals] === False];

(* 2. lem:tk-direction: C = f_r/2 - f f_r/(2(q-f)) - f/r; f_r/2 - f/r = 3M/r^2 - 1/r; middle term <= 0 *)
f = 1 - 2 M/r; fr = D[f, r];
id["tk-direction: f_r/2 - f/r = 3M/r^2 - 1/r", fr/2 - f/r - (3 M/r^2 - 1/r)];
ok["tk-direction: f f_r/(2(q-f)) >= 0 for r > 2M, q > 1", Reduce[Exists[{rr, qq, MM}, MM > 0 && rr > 2 MM && qq > 1 &&
     ((f fr/(2 (q - f))) /. {r -> rr, q -> qq, M -> MM}) < 0], Reals] === False];

(* 3. lem:no-return: V^2 = 1 - f/q > 1 - f = 2M/r for q > 1, 0 < f < 1 *)
id["no-return: 1 - f = 2M/r", 1 - f - 2 M/r];
ok["no-return: 1 - f/q > 1 - f for q > 1, 0 < f < 1", Reduce[Exists[{ff, qq}, 0 < ff < 1 && qq > 1 && 1 - ff/qq <= 1 - ff], Reals] === False];
(* thm:angle: N = 1/(1 - f/q) < 1/(1 - 1/q) <= 3 for q >= 3/2 *)
ok["angle: 1/(1-1/q) <= 3 for q >= 3/2", Reduce[Exists[qq, qq >= 3/2 && 1/(1 - 1/qq) > 3], Reals] === False];

(* 4. excursion: theta = 1 - 2M|a|/(r sqrt(Delta)), r sqrt(Delta) increasing; theta(2M) = 0 *)
Dl = r^2 - 2 M r + a^2;
ok["excursion: d/dr (r^2 Delta) > 0 for r > 2M, 0 < a <= M", Reduce[Exists[{rr, aa, MM}, MM > 0 && 0 < aa <= MM && rr > 2 MM &&
     (D[r^2 Dl, r] /. {r -> rr, a -> aa, M -> MM}) <= 0], Reals] === False];
id["excursion: theta(2M) = 0", (1 - 2 M a/(r Sqrt[Dl])) /. r -> 2 M, M > 0 && a > 0];

(* 5. lem:kerr-G (M = 1): G_phiphi - r^2 = a^2(1 + 2/r) + u Delta/(1 - f u) >= 0; Omega bound *)
ok["kerr-G: a^2(1+2/r) + u Delta/(1-fu) >= 0 (r > r_+, 0 <= u < 1)",
  Reduce[Exists[{aa, rr, uu}, 0 < aa <= 1 && rr > 1 + Sqrt[1 - aa^2] && 0 <= uu < 1 &&
     aa^2 (1 + 2/rr) + uu (rr^2 - 2 rr + aa^2)/(1 - (1 - 2/rr) uu) < 0], Reals] === False];
ok["kerr-G: 2a/r^3 <= 2 and sqrt(Delta)/r^2 <= 1 (r > r_+ >= 1, a <= 1)",
  Reduce[Exists[{aa, rr}, 0 < aa <= 1 && rr > 1 + Sqrt[1 - aa^2] && (2 aa/rr^3 > 2 || Sqrt[rr^2 - 2 rr + aa^2]/rr^2 > 1)], Reals] === False];

(* 6. lem:vaidya-dec-sturm: H_q = z(c0 + c1 z + 40 z^2)/(c2 + 4z); (N'D - N D') with positive coefficients *)
c0 = 3 (2 q - 3) (q - 1); c1 = 33 q - 38; c2 = 3 (q - 1);
num = z (c0 + c1 z + 40 z^2); den = c2 + 4 z;
id["H_q: N'D - ND' = c0 c2 + 2 c1 c2 z + (4c1 + 120 c2) z^2 + 320 z^3",
  Expand[D[num, z] den - num D[den, z] - (c0 c2 + 2 c1 c2 z + (4 c1 + 120 c2) z^2 + 320 z^3)]];
ok["H_q: c0 >= 0, c1 > 0, c2 > 0 for q >= 3/2", Reduce[Exists[qq, qq >= 3/2 && ((c0 /. q -> qq) < 0 || (c1 /. q -> qq) <= 0 || (c2 /. q -> qq) <= 0)], Reals] === False];
ok["H_q: control that must fail, c1 < 0 at q = 11/10 (FAIL expected)", (c1 /. q -> 11/10) > 0];

(* 7. Kerr equator: -f/(2w) <= 1/2, w = 1 - f u, f > -1 on r > r_+ >= 1 *)
ok["-f/(2w) <= 1/2: f > -1, 0 <= u < 1, w > 0", Reduce[Exists[{ff, uu}, ff > -1 && ff < 1 && 0 <= uu < 1 && 1 - ff uu > 0 &&
     -ff/(2 (1 - ff uu)) > 1/2], Reals] === False];
ok["f = 1 - 2/r > -1 for r > 1", Reduce[Exists[rr, rr > 1 && 1 - 2/rr <= -1], Reals] === False];

(* 8. example on the equator with spin, shell [23/4, 25/4], M = 1, a = 1/2: logarithmic derivatives *)
Ds = r^2 - 2 r + 1/4; fs = 1 - 2/r;
id["shell: h_r = r^2/(f Delta) = r^3/((r-2)Delta)", r^2/(fs Ds) - r^3/((r - 2) Ds), r > 2];
id["shell: d log h_r = 3/r - 1/(r-2) - 2(r-1)/Delta", D[Log[r^3/((r - 2) Ds)], r] - (3/r - 1/(r - 2) - 2 (r - 1)/Ds), r > 2];
id["shell: 3/r - 1/(r-2) - 2/(r-1) = (6-4r)/(r(r-1)(r-2))", 3/r - 1/(r - 2) - 2/(r - 1) - (6 - 4 r)/(r (r - 1) (r - 2))];
ok["shell: Delta < (r-1)^2", Simplify[(r - 1)^2 - Ds] === 3/4];
id["shell: d log h_phi = 2(r-1)/Delta - 4/(r(r-2))", D[Log[Ds/fs^2], r] - (2 (r - 1)/Ds - 4/(r (r - 2))), r > 2];
id["shell: 2/(r-1) - 4/(r(r-2)) = 2(r^2-4r+2)/(r(r-1)(r-2))", 2/(r - 1) - 4/(r (r - 2)) - 2 (r^2 - 4 r + 2)/(r (r - 1) (r - 2))];
ok["shell: 23/4 > 2 + sqrt(2)", 23/4 > 2 + Sqrt[2]];
ok["shell: h_r decreasing, h_phi increasing on [23/4,25/4] (direct check)",
  Reduce[Exists[rr, 23/4 <= rr <= 25/4 && ((D[r^3/((r - 2) Ds), r] /. r -> rr) >= 0 || (D[Ds/fs^2, r] /. r -> rr) <= 0)], Reals] === False];

(* 9. S6: chain for D Phi_b, q = 6/5, c = 3 sqrt(6), K = 465, D0 = 2 10^4 (M = 1) *)
c = 3 Sqrt[6]; K = 465; D0 = 2 10^4;
ok["S6: K M^2/D^2 <= c M/D for D >= D0 (D0 > K/c, M = 1)", D0 > K/c];
ok["S6: F_min = 1 - c/D0 - K/D0^2 > 0", 1 - c/D0 - K/D0^2 > 0];
epsm = 300^2/D0^2; epsp = 2000^2/D0^2;
lowB = -2/Sqrt[1 - epsm] + 2 c D0 Sqrt[1 - epsm]/300^2 - Pi K D0/300^3;
uppB = -2 (1 - c/D0 - K/D0^2) + 2 c D0/2000^2 + Pi K D0/2000^3;
Print["   lower bound = ", N[lowB, 12], ";  upper bound = ", N[uppB, 12]];
ok["S6: lower bound > 1832/10000 (printed margin)", lowB > 1832/10000];
ok["S6: upper bound < -1922/1000 (printed margin)", uppB < -1922/1000];
ok["S6: control that must fail, lower bound > 19/100 (FAIL expected)", lowB > 19/100];
(* monotonicity in t = D/D0 >= 1 of each term, with b = beta sqrt(t): b^2/D^2 = beta^2/(D0^2 t) *)
ok["S6: 1/sqrt(1 - beta^2/(D0^2 t)) decreasing in t", Reduce[Exists[{t, be}, t >= 1 && 0 < be <= 2000 &&
     (D[1/Sqrt[1 - be^2/(D0^2 tt)], tt] /. tt -> t) >= 0], Reals] === False];
ok["S6: 1/sqrt(t) <= 1 for t >= 1", Reduce[Exists[t, t >= 1 && 1/Sqrt[t] > 1], Reals] === False];
ok["S6: b_+ = 2000 sqrt(t) <= D/2 = 10^4 t for t >= 1", Reduce[Exists[t, t >= 1 && 2000 Sqrt[t] > 10^4 t], Reals] === False];
ok["S6: 300 > pi K/(2 sqrt(3) c)", 300 > Pi K/(2 Sqrt[3] c)];
ok["S6: rho(8200)^2 = 6 8200^4/(8198 8210) > D0^2", 6 8200^4/(8198 8210) > D0^2];

(* 10. S6, domain of the remainders: q = 6/5, rho^2/M^2 = 6R^4/((R-2)(R+10)) *)
qq = 6/5; fR = 1 - 2/R; rho2 = R^2/(fR (1 - fR/qq));
id["S6: rho^2 = 6R^4/((R-2)(R+10)) for q = 6/5 (rho = r/(sqrt(f) V))", rho2 - 6 R^4/((R - 2) (R + 10)), R > 2];
id["S6: rho(10M)^2 = 375, rho_* = 5 sqrt(15) M", (6 R^4/((R - 2) (R + 10)) /. R -> 10) - 375];
id["S6: sign of (rho^2)' = sign of R^2 + 12R - 40", Factor[D[6 R^4/((R - 2) (R + 10)), R]] - 12 R^3 (R^2 + 12 R - 40)/((R - 2)^2 (R + 10)^2)];
id["S6: positive zero of R^2 + 12R - 40 = 2 sqrt(19) - 6", Expand[(R^2 + 12 R - 40) /. R -> 2 Sqrt[19] - 6]];
ok["S6: rho(3M)^2 = 486/13 < 49", (6 R^4/((R - 2) (R + 10)) /. R -> 3) == 486/13 && 486/13 < 49];
ok["S6: 300 > 5 sqrt(15)", 300 > 5 Sqrt[15]];

Print[""]; Print["PASS ", pass, "  FAIL ", fail, "  (2 FAIL expected: the checks that must fail)"];
