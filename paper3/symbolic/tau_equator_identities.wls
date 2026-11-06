(* Paper III, Thakurta-Kerr on the equator, proper-time branch: exact identities.
   Equatorial data: f = 1-2M/r, Dl = r^2-2Mr+a^2, h0 = r^2/(f Dl) dr^2 + Dl/f^2 dphi^2, beta = bp dphi, bp = -2Ma/(r f),
   N = 1/(1-A^2 f/E^2), omega = A^2 f/E. Rail: deta = beta(dx) + Sqrt[N]|dx|_h0. Relaxed cost: omega (deta - beta(dx)). *)
f = 1 - 2 M/r; Dl = r^2 - 2 M r + a^2; NN = 1/(1 - A^2 f/E^2); bp = -2 M a/(r f); om = A^2 f/E;
h0rr = r^2/(f Dl); h0pp = Dl/f^2;
asm = M > 0 && a > 0 && a < M && A > 0 && E > A && r > 2 M;
pass[name_, expr_] := Print[If[TrueQ[expr], "PASS ", "FAIL "], name];

(* 1. Perlick's metric and the function rho *)
htrr = om^2 NN h0rr; htpp = om^2 NN h0pp;
pass["h_tau = A^4/(E^2-A^2 f)[f r^2/Dl dr^2 + Dl dphi^2]",
  Simplify[{htrr - A^4/(E^2 - A^2 f) f r^2/Dl, htpp - A^4/(E^2 - A^2 f) Dl}, asm] == {0, 0}];
rho = Sqrt[htpp];
pass["rho^2 increasing: d(N Dl)/dr = N' Dl + N Dl' with N' > 0, Dl' > 0",
  Simplify[D[NN, r] - NN^2 A^2 2 M/(E^2 r^2), asm] == 0 && Simplify[D[Dl, r] - 2 (r - M)] == 0];
Print["rho at the boundary r=2M: ", Simplify[(rho /. r -> 2 M), asm]];
(* expansions near the boundary: s = r - 2M *)
serR = Series[Sqrt[htrr] /. r -> 2 M + s, {s, 0, 1}, Assumptions -> s > 0 && asm];
Print["sqrt(h_tau,rr) near the boundary: ", Normal[serR] // Simplify];
serRho = Series[rho /. r -> 2 M + s, {s, 0, 1}];
Print["rho near the boundary: ", Simplify[Normal[serRho], asm]];

(* 2. cost per radian and clock per radian on the ring *)
pass["tau cost per radian = rho", Simplify[om Sqrt[NN] Sqrt[h0pp] - rho, asm] == 0];
kpro = Limit[(Sqrt[NN] Sqrt[Dl] + bp f)/f, r -> 2 M, Direction -> "FromAbove", Assumptions -> asm];
Print["prograde clock per radian at the boundary: ", Simplify[kpro, asm]];
pass["kappa = a + a A^2/(2E^2) + 2M^2/a", Simplify[kpro - (a + a A^2/(2 E^2) + 2 M^2/a), asm] == 0];
Print["retrograde clock per radian, coefficient of 1/s: ",
  Simplify[Limit[s (Sqrt[NN] Sqrt[Dl] - bp f)/f /. r -> 2 M + s, s -> 0, Direction -> "FromAbove", Assumptions -> asm], asm]];

(* 3. norm of beta and maximal radial speed *)
gpp = r^2 + a^2 + 2 M a^2/r;
bnorm2 = bp^2/h0pp;
pass["|beta|^2 = 4M^2a^2/(r^2 Dl)", Simplify[bnorm2 - 4 M^2 a^2/(r^2 Dl), asm] == 0];
pass["1-|beta|^2 = f g_phiphi/Dl", Simplify[1 - bnorm2 - f gpp/Dl, asm] == 0];
(* maximum of |X| on { -b Y + Sqrt[N] Sqrt[X^2+Y^2] <= 1 }: Y* = b/(N-b^2), Xmax = 1/Sqrt[N-b^2] *)
Clear[b, n, X, Y];
FF[Y_] := (1 + b Y)^2/n - Y^2;
Ys = Y /. First@Solve[D[FF[Y], Y] == 0, Y];
pass["Y* = b/(n-b^2) and Xmax^2 = 1/(n-b^2)", Simplify[{Ys - b/(n - b^2), FF[Ys] - 1/(n - b^2)}] == {0, 0}];
pass["constraint saturated at Y*: 1 + b Y* >= 0", Simplify[1 + b Ys - n/(n - b^2)] == 0];
rdotmax = Sqrt[f Dl]/r/Sqrt[NN - bnorm2];
Print["maximal |rdot| at the boundary: ", Simplify[Limit[rdotmax, r -> 2 M, Direction -> "FromAbove", Assumptions -> asm], asm]];
Print["purely radial |rdot| (Y=0) near the boundary, coefficient of Sqrt[s]: ",
  Simplify[Limit[(Sqrt[f Dl]/r/Sqrt[NN] /. r -> 2 M + s)/Sqrt[s], s -> 0, Direction -> "FromAbove", Assumptions -> asm], asm]];

(* 4. singular arc of the PMP: p~ = p + omega beta = 0 on a set of positive measure *)
Clear[rr, ph, et, Af];
omS = Af[et]^2 (1 - 2 M/rr)/E; bpS = -2 M a/(rr - 2 M);   (* beta_phi = -2Ma/(r f) = -2Ma/(r-2M) *)
(* equations: (r) d_r omega (1 - bp phidot) = omega bp' phidot ;  (phi) bp d_eta omega + rdot d_r(omega bp) = 0 *)
eqr = D[omS, rr] (1 - bpS phd) - omS D[bpS, rr] phd;
eqp = bpS D[omS, et] + rd D[omS bpS, rr];
sol = Solve[{eqr == 0, eqp == 0}, {phd, rd}];
Print["singular arc: ", Simplify[sol, M > 0 && a > 0 && rr > 2 M]];
pass["singular arc: phidot = 1/a, rdot = 2 r A'/A",
  Simplify[({phd, rd} /. First[sol]) - {1/a, 2 rr Af'[et]/Af[et]}, M > 0 && a > 0 && rr > 2 M] == {0, 0}];
(* admissibility: beta(v) + Sqrt[N] |v|_h0 <= 1 would require, already with rdot = 0 and N = 1, Sqrt[Dl] <= a *)
pass["with phidot = 1/a: bp phidot + Sqrt[Dl]/(f a) - 1 = r (Sqrt[Dl] - a)/((r-2M) a) > 0",
  Simplify[(bp/a + Sqrt[Dl]/(f a) - 1) - r (Sqrt[Dl] - a)/((r - 2 M) a), asm] == 0 &&
  Resolve[ForAll[{r}, r > 2 M && M > 0 && a > 0, r^2 - 2 M r + a^2 > a^2], Reals]];

(* 5. derivation of the equations (r),(phi) from the condition d/deta (p + omega beta) = 0 with the adjoint equation
      p_i' = d_i omega - d_i(omega beta_j) xdot^j  (valid where p~ = 0) *)
Clear[pr, pp];
xd = {rd, phd}; bet = {0, bpS}; xs = {rr, ph};
padj = Table[D[omS, xs[[i]]] - Sum[D[omS bet[[j]], xs[[i]]] xd[[j]], {j, 2}], {i, 2}];
dtot[g_] := D[g, et] + rd D[g, rr] + phd D[g, ph];
cond = Table[padj[[i]] + dtot[omS bet[[i]]], {i, 2}];
pass["component r = eqr", Simplify[cond[[1]] - eqr] == 0];
pass["component phi = eqp", Simplify[cond[[2]] - eqp] == 0];

(* 6. Clairaut constant: p_phi of the PMP = Clairaut constant of h_tau (constant A) *)
Clear[vr, vp];
v0 = Sqrt[h0rr vr^2 + h0pp vp^2];
pphi = om Sqrt[NN] h0pp vp/v0;
Jt = htpp vp/Sqrt[htrr vr^2 + htpp vp^2];
pass["p_phi = J (Clairaut of h_tau)", Simplify[pphi - Jt, asm && vp > 0] == 0];

(* 7. dynamics of the clock costate: dH/deta = -lambda d_eta Sqrt[N] |v| - d_eta omega (1 - beta(v)), on a rail *)
Clear[Ae];
NNe = 1/(1 - Ae[et]^2 f/E^2); ome = Ae[et]^2 f/E;
lhs = D[ome, et] + (ome - qq) D[NNe, et]/(2 NNe);
pass["q' = (1-beta(v)) (A A' f/E)[2 + N(omega - q)/E]",
  Simplify[lhs - Ae[et] Ae'[et] f/E (2 + NNe (ome - qq)/E), asm] == 0];

(* 8. waiting near the boundary: cost per unit clock omega = A^2 f/E -> 0 *)
Print["omega at the boundary: ", Limit[om, r -> 2 M]];

(* 9. set of limit velocities at the ring: with s = r - 2M, scale rdot ~ 1, phidot ~ 1 *)
Clear[u, w];
constr = (bp w + Sqrt[NN] Sqrt[h0rr u^2 + h0pp w^2] - 1) f;   (* <= 0, multiplied by f > 0 *)
(* expansion: the leading term vanishes for w >= 0; look at the form (bp w + Sqrt[N h0pp] w) for w > 0 and the rest *)
lim = Series[(constr /. r -> 2 M + s), {s, 0, 1}, Assumptions -> s > 0 && asm && w > 0];
Print["constraint multiplied by f, order s (w > 0): ", Simplify[Normal[lim], s > 0 && asm && w > 0]];
