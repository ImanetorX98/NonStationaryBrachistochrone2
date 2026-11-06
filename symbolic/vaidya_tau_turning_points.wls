(* Paper III, proper-time branch of ingoing Vaidya: radial acceleration at the turning points rdot = 0 along the PMP
   flow. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
f = 1 - 2 m/r; w = e^2 - f; om = f/e; S = Sqrt[e^2 pr^2 + J^2/r^2]; h = -w pr + Sqrt[w] S;
as = {e > 1, m > 0, r > 2 m, J > 0, pr > 0};
rd = D[h, pr];
pstar = Sqrt[w] J/(e r Sqrt[f]);
ok["rdot = 0 at p_r = p*", Simplify[(rd /. pr -> pstar), as] == 0];
hpp = Simplify[D[h, {pr, 2}] /. pr -> pstar, as];
prd = Simplify[-D[h - om, r] /. pr -> pstar, as];
dmr = Simplify[D[rd, m] /. pr -> pstar, as];
ok["d_m rdot = -1/r at rdot = 0", Simplify[dmr + 1/r, as] == 0];
Print["h_pp = ", FullSimplify[hpp, as]];
Print["p_r' = ", FullSimplify[prd, as]];
acc = FullSimplify[hpp prd, as];
Print["h_pp p_r' = ", acc];
ok["control: without the constraint h = omega - q, h_pp p_r' > 0 for all r > 2m, J > 0 is false (must fail)", Resolve[ForAll[{r, J}, r > 2 m && J > 0, acc > 0] /. {m -> 1, e -> Sqrt[2]}] === False];
(* rddot|_{rdot=0} = h_pp p_r' - m'/r : threshold of m' *)
Print["threshold m' < r h_pp p_r' : ", FullSimplify[r acc, as]];
(* with the constraint h = omega - q at the turning point: J = r (f - E q)/sqrt(f w) *)
Jc = r (f - e q)/Sqrt[f w];
ok["h(p*) = sqrt(f w) J/(E r)", Simplify[(h /. pr -> pstar) - Sqrt[f w] J/(e r), as] == 0];
accq = FullSimplify[(acc /. J -> Jc), {e > 1, m > 0, r > 2 m, 0 <= q < f/e}];
Print["h_pp p_r' on the constraint = ", accq];
ok["q=0: h_pp p_r' > 0 for every r > 2m (E^2=2, m=1)", Resolve[ForAll[r, r > 2, (accq /. {q -> 0, m -> 1, e -> Sqrt[2]}) > 0]]];
ok["q in [0, omega): h_pp p_r' > 0 (E^2=2, m=1)", Resolve[ForAll[{r, q}, r > 2 && 0 <= q < (1 - 2/r)/Sqrt[2], (accq /. {m -> 1, e -> Sqrt[2]}) > 0]]];
ok["q=0, generic E>1, m=1", Resolve[ForAll[{r, e}, r > 2 && e > 1, (accq /. {q -> 0, m -> 1}) > 0]]];
Print["threshold of m' for no apoapsis (q=0): m' < r h_pp p_r' = ", FullSimplify[r accq /. q -> 0, {e > 1, m > 0, r > 2 m}]];
