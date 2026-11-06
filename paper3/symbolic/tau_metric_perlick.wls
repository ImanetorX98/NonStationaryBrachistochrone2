(* Paper III, proper-time branch: the metric h_tau equals the metric h_C of Perlick 1991, Proposition 3.1:
   h_C = h/(e^{2(C-V)} - 1), with g = h - e^{2V}(psi + dt)^2.
   Thakurta-Kerr with constant A on the equator: g = A^2 g_Kerr; e^{2V} = A^2 f; e^C = Ehat; h = A^2 (f h0),
   h0 = r^2/(f Dl) dr^2 + Dl/f^2 dphi^2. *)
f = 1 - 2 M/r; Dl = r^2 - 2 M r + a^2;
h0rr = r^2/(f Dl); h0pp = Dl/f^2;
e2V = A^2 f; eC = E;
hCrr = A^2 f h0rr/(eC^2/e2V - 1); hCpp = A^2 f h0pp/(eC^2/e2V - 1);
htrr = A^4/(E^2 - A^2 f) f r^2/Dl; htpp = A^4/(E^2 - A^2 f) Dl;
Print[If[Simplify[{hCrr - htrr, hCpp - htpp}] == {0, 0}, "PASS h_tau = h_C of Perlick (Prop. 3.1)", "FAIL"]];
(* control that must fail: the metric of Prop. 3.2 (h~_C = e^{2C-4V} h_C) does not coincide *)
(* exact witness: at M = 1, a = 0, A = 1, E = 2, r = 4 the rr components differ *)
wit = (eC^2/e2V^2 hCrr - htrr) /. {M -> 1, a -> 0, A -> 1, E -> 2, r -> 4};
Print[If[NumericQ[wit] && wit != 0, "PASS negative control: h~_C (arrival-time branch) is different (at M=1, a=0, A=1, E=2, r=4 the rr components differ by " <> ToString[wit, InputForm] <> ")", "FAIL expected failure did not occur"]];
Print["limit at the boundary f -> 0 (a = 0): h_C,pp = ", Limit[hCpp /. a -> 0, r -> 2 M, Direction -> "FromAbove"]];
