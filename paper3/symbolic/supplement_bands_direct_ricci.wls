(* supplement_bands_direct_ricci.wls -- Paper III, supplement (proposition on the certificates, frames b = 2 and
   b = 5/3, bands q >= 5/3 and [11/7, 2]). Independent check of the certified signs: Ricci tensor computed directly
   for the frame N_b ghat with constant b. *)
x = {eta, r, ph}; n = 3;
ricci[g_] := Module[{gi, G},
  gi = Simplify[Inverse[g]];
  G = Table[Simplify[(1/2) Sum[gi[[a, d]] (D[g[[d, b]], x[[c]]] + D[g[[d, c]], x[[b]]] - D[g[[b, c]], x[[d]]]), {d, n}]], {a, n}, {b, n}, {c, n}];
  Table[Simplify[Sum[D[G[[a, b, c]], x[[a]]] - D[G[[a, b, a]], x[[c]]] + Sum[G[[a, a, d]] G[[d, b, c]] - G[[a, c, d]] G[[d, b, a]], {d, n}], {a, n}]], {b, n}, {c, n}]];
f = 1 - 2/r;
Do[
  Nb = 1/(1 - f/b);
  V2 = 1 - f/q[eta];
  g = DiagonalMatrix[{-Nb V2, Nb/f^2, Nb r^2/f}];
  R = ricci[g];
  hasQdot = ! FreeQ[R, Derivative[1][q]] || ! FreeQ[R, Derivative[2][q]];
  (* null k with k^eta = 1, direction alpha measured h0-orthonormally *)
  kr = Cos[al] Sqrt[V2] f; kp = Sin[al] Sqrt[V2] Sqrt[f]/r;
  k = {1, kr, kp};
  Rkk = Simplify[k . R . k /. q[eta] -> Q];
  Cr = Simplify[Rkk /. al -> 0]; Cp = Simplify[Rkk /. al -> Pi/2];
  aff = Simplify[Rkk - (Cr Cos[al]^2 + Cp Sin[al]^2)];
  Print["b = ", b, ": Ric free of q', q'' : ", ! hasQdot, ";  Ric(k,k) = Cr cos^2 + Cp sin^2 : ", aff === 0];
  rng = If[b == 2, Q >= 5/3, 11/7 <= Q <= 2];
  Print["   Cr < 0 on r>2, ", rng, " : ", Reduce[Exists[r, r > 2 && rng && Cr >= 0], Reals] === False];
  Print["   Cp < 0 on r>2, ", rng, " : ", Reduce[Exists[r, r > 2 && rng && Cp >= 0], Reals] === False];
  lo = If[b == 2, 1.6, 1.55]; hi = If[b == 2, 1.66, 2.05];
  Print["   control just outside: max Cp at Q=", lo, " over r: ", NMaxValue[{Cp /. Q -> lo, r > 2.001, r < 10^5}, r], ";  at Q=", hi, ": ", NMaxValue[{Cp /. Q -> hi, r > 2.001, r < 10^5}, r]],
  {b, {2, 5/3}}];
