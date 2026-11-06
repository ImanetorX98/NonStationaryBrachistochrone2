(* Paper III, proper-time branch of ingoing Vaidya: kinematics of the velocity body K. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
(* body K: (s+w)^2/E^2 + r^2 zeta^2 <= w ; w = E^2 - 1 + 2m/r *)
sp = -w + e Sqrt[w]; sm = -w - e Sqrt[w];
ok["angular speed r|zeta| <= sqrt(w)", Resolve[ForAll[{s, z, w, e}, w > 0 && e > 0 && (s + w)^2/e^2 + z^2 <= w, z^2 <= w]]];
ok["s+ decreases in w if sqrt(w) > E/2", Simplify[D[sp, w] < 0, {e > 0, w > e^2/4}]];
ok["s- decreases in w", Simplify[D[sm, w] < 0, {e > 0, w > 0}]];
ok["w > E^2/4 always if E^2 > 4/3 (w >= E^2 - 1)", Resolve[ForAll[e, e^2 > 4/3, e^2 - 1 > e^2/4]]];
ok["0 interior to K iff 0 < f < E^2 (w < E^2 and w > 0): (0+w)^2/E^2 < w", Resolve[ForAll[{w, e}, 0 < w < e^2 && e > 0, w^2/e^2 < w]]];
