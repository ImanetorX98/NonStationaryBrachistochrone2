(* supplement_frozen_slice.wls -- Paper III, supplement (frozen factor).
   Frozen A, a = 0: slice metric a_T = N h0 of the ultrastatic frame -deta^2 + N h0.
   1. its Gauss curvature equals the optical Gauss curvature K_t of Paper II;
   2. K_t < 0 for all r > 2M when q >= 3/2, and K_t > 0 somewhere far out when 1 < q < 3/2;
   3. completeness: radial distance  int sqrt(N)/f dr  diverges at r -> 2M and at infinity;
      angular radius sqrt(N) r / sqrt(f) bounded below (no escape in the angle);
   4. CONTROL (must fail): at q = 1.4 the curvature is positive at r = 100 M.  *)
$Failures = 0;
check[l_, c_] := (Print[If[TrueQ[c], "PASS  ", "FAIL  "], l]; If[! TrueQ[c], $Failures++]);
f = 1 - 2/r; Nq = 1/(1 - f/q);
a = {{Nq/f^2, 0}, {0, Nq r^2/f}};
(* Gauss curvature of a diagonal 2-metric E dr^2 + G dphi^2 depending on r only *)
EE = a[[1, 1]]; GG = a[[2, 2]];
K = Simplify[-1/(2 Sqrt[EE GG]) D[D[GG, r]/Sqrt[EE GG], r]];
Kt = (16 + 2*(9*q - 14)*r + (3*q^2 - 19*q + 16)*r^2 - (2*q - 3)*(q - 1)*r^3)/(q r^5 (2 + (q - 1) r));
check["1. K(N h0) = K_t (optical Gauss curvature of Paper II, M = 1)", Simplify[K - Kt] === 0];
check["2a. K_t < 0 for all r > 2 when q >= 3/2", Reduce[Exists[r, r > 2 && q >= 3/2 && Kt >= 0], Reals] === False];
check["2b. for 1 < q < 3/2, K_t > 0 somewhere", Reduce[ForAll[q, 1 < q < 3/2, Exists[r, r > 2, Kt > 0]], Reals] === True];
radial = Simplify[Sqrt[Nq]/f, r > 2 && q > 1];
lim2 = Simplify[Limit[(r - 2) radial, r -> 2, Direction -> "FromAbove"], q > 1];
check["3a. radial integrand ~ c/(r-2), c > 0: log divergence at r -> 2 (c = " <> ToString[lim2] <> ")", Simplify[lim2 > 0, q > 1]];
liminf = Simplify[Limit[radial, r -> Infinity], q > 1];
check["3b. radial integrand -> positive constant at infinity: linear divergence (" <> ToString[liminf] <> ")", Simplify[liminf > 0, q > 1]];
ang = Simplify[Sqrt[Nq] r/Sqrt[f], r > 2 && q > 1];
check["3c. angular radius bounded below by a positive constant (q = 3/2, r > 2)", NMinimize[{ang /. q -> 3/2, r > 2.0001}, r][[1]] > 1];
check["4. CONTROL: q = 1.4, r = 100: K_t > 0", (Kt /. {q -> 1.4, r -> 100}) > 0];
Print[]; If[$Failures == 0, Print["ALL CHECKS PASSED"], Print[$Failures, " CHECK(S) FAILED"]];
Exit[If[$Failures == 0, 0, 1]];
