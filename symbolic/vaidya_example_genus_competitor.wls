(* Paper III, constants of the Vaidya example (supplement S11; theorem on late accretion, D = 91, s = 11/5). M = 1, E^2 = 2. *)
ok[n_, c_] := Print[If[TrueQ[c], "PASS ", "FAIL "], n];
(* M = 1, E^2 = 2: rho^2 = r^2 (r-2)/(r+2), R' = sqrt(r/(r+2)) *)
rho2[r_] := r^2 (r - 2)/(r + 2); Rp[r_] := Sqrt[r/(r + 2)];
ok["R'/rho = 1/sqrt(r(r-2))", Simplify[Rp[r]/Sqrt[rho2[r]] - 1/Sqrt[r (r - 2)], r > 2] == 0];
(* 1. genus of the integrals of the dips: angular integrand L R'/(rho sqrt(rho^2-L^2)) *)
num = Factor[Together[(rho2[r] - rho2[t])]];
Print["rho^2(r) - rho^2(t) = ", num];
rad = Factor[r (r - 2) Numerator[num] (r + 2)];
Print["radicand (up to squares) = ", rad, "   degree in r: ", Exponent[rad, r]];
Print["   => degree 6 in r without generic double roots: hyperelliptic integral of genus 2, not elliptic"];
ok["discriminant non-zero for t = 3 (simple roots)", Discriminant[rad /. t -> 3, r] != 0];
(* 2. competitor with constant direction: phi' = tan(th) R'/rho, ds = sec(th) R' dr *)
A[r_] := 2 ArcCosh[Sqrt[r/2]];
ok["dA/dr = 1/sqrt(r(r-2))", Simplify[D[A[r], r] - 1/Sqrt[r (r - 2)], r > 2] == 0];
F[r_] := Sqrt[r (r + 2)] - 2 ArcSinh[Sqrt[r/2]];
ok["dF/dr = R'", Simplify[D[F[r], r] - Rp[r], r > 2] == 0];
(* dip with constant direction from 6 to r1 and back: angle 2 tan(th)(A(6)-A(r1)) = pi/2, length 2 sec(th)(F(6)-F(r1)) *)
len[r1_] := Module[{tt = (Pi/4)/(A[6] - A[r1])}, 2 Sqrt[1 + tt^2] (F[6] - F[r1])];
best = NMinimize[{len[r1], 2 < r1 < 6}, r1, WorkingPrecision -> 30];
Print["competitor with constant direction: minimal length = ", N[best[[1]], 12], " at r1 = ", N[r1 /. best[[2]], 10]];
Print["   thresholds: Lambda(11/5) = 5.9412053, 91/(11 sqrt2) = ", N[91/(11 Sqrt[2]), 10], ", C0 [N] = 5.8130366"];
