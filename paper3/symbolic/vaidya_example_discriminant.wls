(* Paper III, constants of the Vaidya example (supplement S11; theorem on late accretion, D = 91, s = 11/5). M = 1, E^2 = 2. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
B = (t + 2) r^2 + (t^2 - 4) r + 2 t^2 - 4 t; H = r (r - 2) (r + 2) (r - t) B;
ok["Disc_r H = 2^24 t^10 (t-2)^5 (t+2)^3 (t^2-8t-4)(t^2+2t-4)^2", Factor[Discriminant[H, r]] == Factor[2^24 t^10 (t - 2)^5 (t + 2)^3 (t^2 - 8 t - 4) (t^2 + 2 t - 4)^2]];
ok["no factor vanishes on [29/10, 297/100]", Resolve[ForAll[t, 29/10 <= t <= 297/100, t (t - 2) (t + 2) (t^2 - 8 t - 4) (t^2 + 2 t - 4) != 0]]];
ok["angular differential = t sqrt(t-2) (r+2)/sqrt(H)", Simplify[(Sqrt[t^2 (t - 2)/(t + 2)] Sqrt[r/(r + 2)]/(Sqrt[r^2 (r - 2)/(r + 2)] Sqrt[r^2 (r - 2)/(r + 2) - t^2 (t - 2)/(t + 2)])) - t Sqrt[t - 2] (r + 2)/Sqrt[H], r > t > 2] == 0];
ok["91^2 - 2 (11*5839/1000)^2 = 15135559/500000", 91^2 - 2 (11 5839/1000)^2 == 15135559/500000];
