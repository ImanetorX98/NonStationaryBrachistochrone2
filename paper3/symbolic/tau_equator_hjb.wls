(* Paper III, Thakurta-Kerr on the equator, proper-time branch with variable A: exact Hamilton-Jacobi-Bellman equation,
   covariant form, checks. h0-orthonormal frame in the plane. *)
pass[name_, e_] := Print[If[TrueQ[e], "PASS ", "FAIL "], name];
(* 1. support function of K = {v : Sqrt[n]|v| + b.v <= 1}, |b|^2 < n *)
sig[xi_, b_, n_] := (Sqrt[(xi.b)^2 + (n - b.b) xi.xi] - xi.b)/(n - b.b);
SeedRandom[7]; ok = True;
Do[
  n = RandomReal[{1, 3}]; b = RandomReal[{-1, 1}, 2] 0.95/Sqrt[2]; xi = RandomReal[{-2, 2}, 2];
  (* boundary of K: v = u/(Sqrt[n] + b.u), u unit vector; maximum over th with a grid and refinement *)
  h[th_?NumericQ] := xi.{Cos[th], Sin[th]}/(Sqrt[n] + b.{Cos[th], Sin[th]});
  th0 = First@MaximalBy[Range[0., 2 Pi, 2 Pi/2000], h];
  num = h[th /. Last@FindMaximum[h[th], {th, th0}, WorkingPrecision -> 40, PrecisionGoal -> 15, AccuracyGoal -> 15]];
  If[Abs[num - sig[xi, b, n]] > 10^-10, ok = False], {40}];
pass["sigma_K = [sqrt((xi.b)^2 + (n-|b|^2)|xi|^2) - xi.b]/(n-|b|^2)  (40 random cases, refined maximum on the boundary)", ok];
(* 2. equivalence: sigma_K(xi) = lam  <=>  |xi - lam b| = lam Sqrt[n]  (lam > 0) *)
Clear[x1, x2, b1, b2, n, lam];
xi = {x1, x2}; b = {b1, b2};
s = sig[xi, b, n];
pass["|xi - sigma b|^2 = n sigma^2", Simplify[(xi - s b).(xi - s b) - n s^2, n > b1^2 + b2^2] === 0];
(* 3. HJB <=> covariant eikonal. V remaining cost: V_eta + omega - sigma_K(omega beta - gradV) = 0
      <=> |gradV + V_eta beta| = Sqrt[n] (omega + V_eta)  <=>  gopt^{-1}(theta + dV, theta + dV) = 0, theta = omega (deta - beta) *)
Clear[om, Ve, g1, g2];
lhsH = sig[om b - {g1, g2}, b, n] - (om + Ve);
cov = -(om + Ve)^2 + ({g1, g2} + Ve b).({g1, g2} + Ve b)/n;
(* on the surface lam = om + Ve > 0: sigma(xi) = lam  <=>  |xi - lam b|^2 = n lam^2 ; xi - lam b = -(gradV + Ve b) *)
pass["xi - (om+Ve) b = -(gradV + Ve b) with xi = om b - gradV", Simplify[(om b - {g1, g2}) - (om + Ve) b + ({g1, g2} + Ve b)] === {0, 0}];
pass["n * cov = |gradV + Ve b|^2 - n (om+Ve)^2 (covariant form)", Simplify[n cov - ((({g1, g2} + Ve b).({g1, g2} + Ve b)) - n (om + Ve)^2)] === 0];
(* inverse of gopt in the coframe {e0 = deta - beta, dx}: gopt^{-1} = -E0 E0 + n^{-1} h0^{ij} E_i E_j, E0 = d_eta, E_i = d_i + beta_i d_eta *)
G = {{-1 + b.b/n, b1/n, b2/n}, {b1/n, 1/n, 0}, {b2/n, 0, 1/n}};
gopt = {{-1, b1, b2}, {b1, n - b1^2, -b1 b2}, {b2, -b1 b2, n - b2^2}};   (* -(deta - b.dx)^2 + n |dx|^2 *)
pass["inverse of gopt", Simplify[G.gopt] === IdentityMatrix[3]];
alpha = {om + Ve, -om b1 + g1, -om b2 + g2};   (* theta + dV, theta = om(deta - b.dx) *)
pass["gopt^{-1}(theta+dV, theta+dV) = cov", Simplify[alpha.G.alpha - cov] === 0];
(* 4. frozen case: V = static d, |grad d|_h0 = om Sqrt[n] (eikonal of h_tau = om^2 n h0) solves the HJB *)
pass["frozen: Ve = 0, |grad d| = om Sqrt[n]  =>  sigma(om b - grad d) = om  (exact in 64 directions)",
  And @@ Table[FullSimplify[(sig[om b - om Sqrt[n] {Cos[t], Sin[t]}, b, n] - om) /. {n -> 2, b1 -> 3/10, b2 -> -1/5, om -> 7/10}] == 0, {t, 0, 2 Pi - Pi/32, Pi/32}]];
(* 5. subsolution from a frozen metric with alpha <= A:  |grad d_alpha| = om_a Sqrt[n_a] <= om_A Sqrt[n_A],
      and |gradV + Ve b| <= |gradV| + Ve |b| <= Sqrt[n](om + Ve) if Ve >= 0 and |b| < 1 <= Sqrt[n] *)
Clear[A, al, EE, f];
dens[x_] := x^2/Sqrt[EE^2 - x^2 f];    (* om Sqrt[n] / f  = A^2/(E Sqrt[1 - A^2 f/EE^2]) /E ... density per unit of |.|_k *)
pass["tau density increasing in A: d/dA [A^2/Sqrt[EE^2-A^2 f]] > 0",
  Simplify[D[A^2/Sqrt[EE^2 - A^2 f], A] > 0, EE > 0 && A > 0 && 0 <= f <= 1 && EE^2 > A^2]];
(* 6. ratio of the densities *)
pass["density ratio = (a^2/A^2) Sqrt[(EE^2-A^2 f)/(EE^2-a^2 f)] (squares, both positive)",
  Simplify[((al^2/Sqrt[EE^2 - al^2 f])/(A^2/Sqrt[EE^2 - A^2 f]))^2 - ((al^2/A^2))^2 (EE^2 - A^2 f)/(EE^2 - al^2 f), EE > 0 && A > 0 && al > 0 && 0 <= f <= 1 && EE^2 > A^2 && EE^2 > al^2] === 0];
