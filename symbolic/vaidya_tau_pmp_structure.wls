(* Paper III, proper-time branch of ingoing Vaidya (sec:vaidya-tau): structure of the Pontryagin maximum principle.
   Notation of the text: w = E^2 - f, f = 1-2m/r, omega = f/E, h = -w p_r + sqrt(w) S, S = sqrt(E^2 p_r^2 + J^2/r^2),
   Hext = q + h - omega = 0. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
f = 1 - 2 m/r; w = e^2 - f; om = f/e; S = Sqrt[e^2 pr^2 + J^2/r^2]; h = -w pr + Sqrt[w] S;
as = {e > 1, m > 0, r > 2 m, J > 0, Element[pr, Reals]};
(* 1. dm(h - omega) with h = omega - q *)
dm = D[h - om, m];
ok["dm(h-om) = (2/r)[-p_r + S/(2 sqrt w) + 1/E]", Simplify[dm - (2/r) (-pr + S/(2 Sqrt[w]) + 1/e), as] == 0];
(* bound: for q>=0 (h <= omega) and p_r > 0, -p_r/2 + h/(2w) + 1/E >= (E - sqrt w)/(2 E w)... check of the key step *)
x =.; ok["1/E[1-(E+x)(2x-E)/(2x^2)] = (E-x)/(2 x^2)", Simplify[1/e (1 - (e + x) (2 x - e)/(2 x^2)) - (e - x)/(2 x^2)] == 0];
(* direct numerical check: dm(h-om) > 0 on h = omega - q, 0 <= q < omega *)
cnt = 0; bad = 0; SeedRandom[7];
Do[Module[{ee = RandomReal[{1.01, 3}], mm = 1, rr, JJ, qq, sol, prs},
   rr = RandomReal[{2.001, 40}]; JJ = RandomReal[{0, 20}]; qq = RandomReal[{0, 0.999}] (1 - 2/rr)/ee;
   prs = pr /. NSolve[(h /. {e -> ee, m -> mm, r -> rr, J -> JJ}) == (1 - 2/rr)/ee - qq, pr, Reals, WorkingPrecision -> 30];
   If[ListQ[prs], Do[cnt++; If[(dm /. {e -> ee, m -> mm, r -> rr, J -> JJ, pr -> p}) <= 0, bad++], {p, prs}]]], {4000}];
Print["   [N] samples on the surface h=omega-q, q in [0,omega): ", cnt, " ; dm(h-om)<=0: ", bad];
ok["[N] dm(h-om) > 0 if q >= 0", bad == 0 && cnt > 1000];
(* 2. radial speed at p_r = 0: still negative (translated body) *)
rdot = D[h, pr];
ok["rdot(p_r=0) = -w < 0", Simplify[(rdot /. pr -> 0) + w, as] == 0];
(* 3. p_r' = -d_r(h - omega) at p_r = 0 is positive *)
ok["p_r'(p_r=0) > 0", Simplify[(-D[h - om, r] /. pr -> 0) > 0, as]];
