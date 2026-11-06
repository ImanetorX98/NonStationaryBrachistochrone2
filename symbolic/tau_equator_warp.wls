(* Paper III, Thakurta-Kerr on the equator, proper-time branch: drho/dR = [r-M + M Dl/(r^2(Q-f))]/(r Sqrt[f]) > 1
   (used with the theorem on a single maximum and the Maxwell rule with a ring). *)
f = 1 - 2 M/r; Dl = r^2 - 2 M r + a^2; NN = Q/(Q - f); c = cc;
rho = c Sqrt[NN Dl]; Rr = c r Sqrt[NN f]/Sqrt[Dl];
asm = M > 0 && a > 0 && Q > 1 && r > 2 M && cc > 0;
cod = (r - M + M Dl/(r^2 (Q - f)))/(r Sqrt[f]);
Print[If[Simplify[D[rho, r]/Rr - cod, asm] == 0, "PASS", "FAIL"], " drho/dR = the stated formula"];
Print[If[Simplify[(r - M)^2 - r^2 f - M^2] == 0, "PASS", "FAIL"], " (r-M)^2 - r^2 f = M^2"];
Print[If[Resolve[ForAll[{r}, r > 2 M && M > 0, r - M > r Sqrt[1 - 2 M/r]], Reals], "PASS", "FAIL"], " r - M > r Sqrt[f] for r > 2M"];
(* negative control: the formula with the opposite sign differs from the true one by the residual
   2 M Dl/(r^{3/2} Sqrt[r-2M] (2M + (Q-1) r)); the residual and its strict positivity are checked; any non-Boolean
   outcome is a FAIL *)
res = Simplify[D[rho, r]/Rr - (r - M - M Dl/(r^2 (Q - f)))/(r Sqrt[f]), asm];
chk1 = Simplify[res - 2 M (r^2 - 2 M r + a^2)/(r^(3/2) Sqrt[r - 2 M] (2 M + (Q - 1) r)), asm];
chk2 = Simplify[res > 0, asm];
Print[If[chk1 === 0 && chk2 === True, "PASS", "FAIL"], " negative control: residual of the opposite sign = 2M Dl/(r^(3/2)Sqrt[r-2M](2M+(Q-1)r)) > 0"];
(* representation of Paper II with spin: rho^2 = A^2 V_a, R'/rho = K_a, with
   V_a = r Dl/((Q-1) r + 2M), K_a = Sqrt[r(r-2M)]/Dl; a = 0 gives back V, K of the text (thm:single-max) *)
Va = r Dl/((Q - 1) r + 2 M); Ka = Sqrt[r (r - 2 M)]/Dl;
Print[If[Simplify[{rho^2 - (cc^2 Q) Va, Rr/rho - Ka}, asm] === {0, 0}, "PASS", "FAIL"],
  " rho^2 = A^2 V_a (A^2 = cc^2 Q) and R'/rho = K_a"];
Print[If[Simplify[{Va /. a -> 0, Ka /. a -> 0} - {r^2 (r - 2 M)/((Q - 1) r + 2 M), 1/Sqrt[r (r - 2 M)]}, asm] === {0, 0}, "PASS", "FAIL"],
  " a = 0: V and K of the text"];
Print[If[Simplify[(Va /. r -> 2 M) - a^2/Q, asm] === 0, "PASS", "FAIL"], " V_a(2M) = a^2/Q (ring)"];
