(* supplement_frame_family.wls -- Paper III, supplement: continuous family of frame factors b.
   Frame gbar_b = -N_b V^2 deta^2 + N_b h0, V^2 = 1 - f/q, N_b = 1/(1 - f/b), M = 1.
   Formula (eq:ric-static): Ric(k,k) = Vt^2 K_{ht} + Vt Hess_{ht} Vt (e_perp, e_perp),
   Vt^2 = N_b V^2, ht = N_b h0, k^eta = 1.  Checks:
   1. general C_r(q,b,r), C_phi(q,b,r) from the formula agree with the direct Ricci for sample (q,b);
   2. at q = b, C_r = C_phi = K_t(b) (ultrastatic frame);
   3. asymptotics: lim r^3 C_phi and lim r^3 C_r as r -> oo, as functions of (q,b);
   4. limits at r -> 2M;
   5. for b > 3/2: at q = b both r^3 C are < 0 on [2, oo] (including the limits) -> by compactness
      there is eps(b) > 0 with Ric(k,k) < 0 for |q - b| < eps(b);
   6. numerical bands [q_-(b), q_+(b)] for sample b (diagnostic). *)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l]; If[! TrueQ[c], $Failures++]);
f = 1 - 2/r;
Nb = 1/(1 - f/b); V2 = 1 - f/q;
Vt = Sqrt[Nb V2];
(* 2D metric ht = E dr^2 + G dphi^2 *)
EE = Nb/f^2; GG = Nb r^2/f;
Kht = Simplify[-1/(2 Sqrt[EE GG]) D[D[GG, r]/Sqrt[EE GG], r]];
(* Hessian of Vt (function of r): components Hrr = Vt'' - Gamma^r_rr Vt', Hpp = - Gamma^r_pp Vt' *)
Grrr = D[EE, r]/(2 EE); Grpp = -D[GG, r]/(2 EE);
Hrr = D[Vt, {r, 2}] - Grrr D[Vt, r]; Hpp = -Grpp D[Vt, r];
(* e radial -> e_perp tangential: Hess(e_phi,e_phi) = Hpp/GG ; e tangential -> e_perp radial: Hrr/EE *)
Cr = Simplify[Vt^2 Kht + Vt Hpp/GG];
Cp = Simplify[Vt^2 Kht + Vt Hrr/EE];
(* 1. direct Ricci for a sample *)
x = {eta, rr, ph};
ricci[g_] := Module[{gi, G},
  gi = Simplify[Inverse[g]];
  G = Table[Simplify[(1/2) Sum[gi[[a, d]] (D[g[[d, bb]], x[[c]]] + D[g[[d, c]], x[[bb]]] - D[g[[bb, c]], x[[d]]]), {d, 3}]], {a, 3}, {bb, 3}, {c, 3}];
  Table[Simplify[Sum[D[G[[a, bb, c]], x[[a]]] - D[G[[a, bb, a]], x[[c]]] + Sum[G[[a, a, d]] G[[d, bb, c]] - G[[a, c, d]] G[[d, bb, a]], {d, 3}], {a, 3}]], {bb, 3}, {c, 3}]];
Do[
  {qs, bs} = s;
  fr = 1 - 2/rr; Nbr = 1/(1 - fr/bs); V2r = 1 - fr/qs;
  g = DiagonalMatrix[{-Nbr V2r, Nbr/fr^2, Nbr rr^2/fr}];
  R = ricci[g];
  kr = {1, Sqrt[V2r] fr, 0}; kp = {1, 0, Sqrt[V2r] Sqrt[fr]/rr};
  check["1. sample (q,b) = " <> ToString[s] <> ": formula = direct Ricci (radial and tangential)",
    Simplify[(kr . R . kr) - (Cr /. {q -> qs, b -> bs, r -> rr})] === 0 && Simplify[(kp . R . kp) - (Cp /. {q -> qs, b -> bs, r -> rr})] === 0],
  {s, {{2, 7/4}, {17/10, 2}}}];
Kt = (16 + 2*(9*q - 14)*r + (3*q^2 - 19*q + 16)*r^2 - (2*q - 3)*(q - 1)*r^3)/(q r^5 (2 + (q - 1) r));
check["2. q = b: C_r = C_phi = K_t", Simplify[(Cr /. b -> q) - Kt] === 0 && Simplify[(Cp /. b -> q) - Kt] === 0];
aR = Simplify[Limit[r^3 Cr, r -> Infinity], q > 1 && b > 1];
aP = Simplify[Limit[r^3 Cp, r -> Infinity], q > 1 && b > 1];
Print["      lim r^3 C_r   = ", Factor[aR]];
Print["      lim r^3 C_phi = ", Factor[aP]];
lR = Simplify[Limit[Cr, r -> 2, Direction -> "FromAbove"], q > 1 && b > 1];
lP = Simplify[Limit[Cp, r -> 2, Direction -> "FromAbove"], q > 1 && b > 1];
Print["      C_r(r->2M) = ", Factor[lR], ";   C_phi(r->2M) = ", Factor[lP]];
check["5a. q = b > 3/2: lim r^3 C < 0 at infinity", Reduce[Exists[q, q > 3/2 && ((aR /. b -> q) >= 0 || (aP /. b -> q) >= 0)], Reals] === False];
check["5b. q = b > 3/2: C < 0 at r = 2M", Reduce[Exists[q, q > 3/2 && ((lR /. b -> q) >= 0 || (lP /. b -> q) >= 0)], Reals] === False];
(* 6. numerical bands *)
Fr = Compile[{{qq, _Real}, {bb, _Real}, {r0, _Real}}, Evaluate[r^3 Cr /. {q -> qq, b -> bb, r -> r0}]];
Fp = Compile[{{qq, _Real}, {bb, _Real}, {r0, _Real}}, Evaluate[r^3 Cp /. {q -> qq, b -> bb, r -> r0}]];
rs = Join[Range[2.0001, 3, 0.002], Range[3, 50, 0.05], {100, 300, 1000, 10^4, 10^5}];
ok[qq_, bb_] := Max[Table[Max[Fr[qq, bb, r0], Fp[qq, bb, r0]], {r0, rs}]] < 0;
Do[
  qs = Select[Range[1.501, 4, 0.002], ok[#, bb] &];
  Print["      b = ", bb, ": band (grid) = ", If[qs === {}, "none", {Min[qs], Max[qs]}]],
  {bb, {1.52, 1.55, 1.6, 1.65, 5/3. , 1.75, 2., 2.5, 3.}}];
Print[]; If[$Failures == 0, Print["ALL CHECKS PASSED"], Print[$Failures, " CHECK(S) FAILED"]];
Exit[If[$Failures == 0, 0, 1]];
