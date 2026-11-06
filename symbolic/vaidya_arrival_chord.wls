(* Paper III, arrival branch of ingoing Vaidya below 3/2: the chord bound g_q(f) >= mu_*(q) used in the text. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
Qw = 2 w^2 - w (5 q/2 - 1) + q (q - 1)/2;
g = (1 - w/q) (-Qw) /. w -> q - f;
ok["g_q(f) = -2f^3/q + (3/2 + 1/q) f^2 - f/2", Simplify[g - (-2 f^3/q + (3/2 + 1/q) f^2 - f/2)] == 0];
ok["g(1/3) = 1/(27q), g(1) = (q-1)/q", Simplify[(g /. f -> 1/3) - 1/(27 q)] == 0 && Simplify[(g /. f -> 1) - (q - 1)/q] == 0];
chord = 3/2 (1 - f) 1/(27 q) + 3/2 (f - 1/3) (q - 1)/q + (f - 1/3) (1 - f) (12 f + 10 - 9 q)/(6 q);
ok["identity of the chord", Simplify[g - chord] == 0];
ok["12f + 10 - 9q > 1/2 for f >= 1/3, q < 3/2", Resolve[ForAll[{f, q}, f >= 1/3 && q < 3/2, 12 f + 10 - 9 q > 1/2]]];
ok["mu* = (q-1)/q per q <= 28/27, 1/(27q) oltre", Simplify[((q - 1)/q - 1/(27 q)) /. q -> 28/27] == 0];
ok["max_q mu* = 1/28", 1/(27 28/27) == 1/28];
