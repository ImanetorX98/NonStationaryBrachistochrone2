(* Paper III, Thakurta-Kerr on the equator, proper-time branch: purely angular motion at radius r (the ring of the
   Maxwell rule with a ring). g = A^2 g_Kerr, rail: d eta - beta = sqrt(N)|dx|_h0, beta = b dphi, b = -2 M a/(r f),
   h0_pp = Dl/f^2, N = 1/(1 - A^2 f/E^2); d tau = (A^2 f/E)(d eta - beta). Prints the quantities; no PASS/FAIL. *)
f = 1 - 2 M/r; Dl = r^2 - 2 M r + a^2; NN = 1/(1 - A^2 f/E^2); b = -2 M a/(r f);
dtau = A^2 f/E Sqrt[NN] Sqrt[Dl]/f;            (* per radian, both senses *)
detaPro = Sqrt[NN] Sqrt[Dl]/f + b;               (* d eta per radian, prograde sense (dphi > 0), a > 0 *)
detaRetro = Sqrt[NN] Sqrt[Dl]/f - b;
asm = M > 0 && 0 < a < M && 0 < A < E;
Print["tau cost per radian: ", Simplify[dtau, asm && r > 2 M], "   limit r -> 2M: ", Limit[dtau, r -> 2 M, Direction -> "FromAbove", Assumptions -> asm]];
Print["prograde clock at the limit: ", Simplify[Limit[detaPro, r -> 2 M, Direction -> "FromAbove", Assumptions -> asm], asm]];
Print["retrograde clock at the limit: ", Limit[detaRetro, r -> 2 M, Direction -> "FromAbove", Assumptions -> asm]];
(* is the angular cost minimal at the ring? derivative in r positive for r > 2M *)
dd = Simplify[D[dtau^2, r], asm && r > 2 M];
Print["d/dr (angular cost)^2 = ", Factor[dd]];
Print["positive for r > 2M (a < M): ", Resolve[ForAll[{r}, r > 2, (dd /. {M -> 1, a -> 1/2, A -> 1, E -> 3/2}) > 0], Reals],
  "  (check at fixed parameters)"];
