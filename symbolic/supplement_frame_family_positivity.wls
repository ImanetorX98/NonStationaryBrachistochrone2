(* supplement_frame_family_positivity.wls -- Paper III, supplement (proposition on the family of frames): exact sign
   over the whole family (b, q) of reference frames.
   Claim: for 3/2 < b and q_-(b) <= q <= q_+(b), q_- = (4b-3)/(2b-1), q_+ = (3-b)/(4-2b) (b<2) or oo (b>=2),
   C_r < 0 and C_phi < 0 for every r > 2 (M = 1). *)
f = 1 - 2/r; Nb = 1/(1 - f/b); V2 = 1 - f/q; Vt = Sqrt[Nb V2];
EE = Nb/f^2; GG = Nb r^2/f;
Kht = -1/(2 Sqrt[EE GG]) D[D[GG, r]/Sqrt[EE GG], r];
Grrr = D[EE, r]/(2 EE); Grpp = -D[GG, r]/(2 EE);
Hrr = D[Vt, {r, 2}] - Grrr D[Vt, r]; Hpp = -Grpp D[Vt, r];
Cr = Together[Simplify[Vt^2 Kht + Vt Hpp/GG, r > 2 && q > 1 && b > 1]];
Cp = Together[Simplify[Vt^2 Kht + Vt Hrr/EE, r > 2 && q > 1 && b > 1]];
nr = Numerator[Cr]; dr = Denominator[Cr]; np = Numerator[Cp]; dp = Denominator[Cp];
Print["denominators: ", Factor[dr], "  |  ", Factor[dp]];
Print["numerators: ", Factor[nr], "  |  ", Factor[np]];
band = (b > 3/2 && b < 2 && (4 b - 3)/(2 b - 1) <= q <= (3 - b)/(4 - 2 b)) || (b >= 2 && q >= (4 b - 3)/(2 b - 1));
t1 = AbsoluteTiming[Reduce[Exists[{r, q, b}, r > 2 && band && Cr >= 0], Reals]];
Print["C_r >= 0 somewhere in the family band? ", t1];
t2 = AbsoluteTiming[Reduce[Exists[{r, q, b}, r > 2 && band && Cp >= 0], Reals]];
Print["C_phi >= 0 somewhere in the family band? ", t2];
