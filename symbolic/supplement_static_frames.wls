#!/usr/bin/env wolframscript
(* Paper III, supplement (convex fronts): Ric in the static frames g_b, certificates, family of frames, bands. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Supplement: static frames =="];
(* Ric(k,k) in g_b = -N_b V^2 deta^2 + N_b h0 (plane), generic q(eta): check of lem:static-frame and of the formulas *)
xs = {eta, r, ph};
f = 1 - 2 M/r;
gb[b_] := Module[{Nb = If[b === Infinity, 1, 1/(1 - f/b)]}, DiagonalMatrix[{-Nb (1 - f/q[eta]), Nb/f^2, Nb r^2/f}]];
Ckk[b_] := Module[{g = gb[b], Ric, V2, e},
   Ric = ricci[g, xs];
   V2 = -g[[1, 1]];
   (* k = (1, sqrt(V2) e), e = cos al e_r + sin al e_phi (h_b-unit) *)
   e = {0, Cos[al]/Sqrt[g[[2, 2]]], Sin[al]/Sqrt[g[[3, 3]]]};
   k = {1, 0, 0} + Sqrt[V2] e;
   Simplify[k.Ric.k, M > 0 && r > 2 M && q[eta] > 1]];
ass = M > 0 && r > 2 M && q[eta] > 1;
Cinf = Ckk[Infinity];
chk["lem:static-frame: Ric_{g_inf}(k,k) without derivatives of q", Boole[!FreeQ[Cinf, Derivative]], ass];
Crinf = Cinf /. al -> 0; Cpinf = Cinf /. al -> Pi/2;
chk["Ric = C_r cos^2 + C_phi sin^2 (no mixed term)", Cinf - (Crinf Cos[al]^2 + Cpinf Sin[al]^2), ass];
qq = q[eta];
chk["prop:certificates(1): C_r = -M{Mq + (2q-1)(r-2M)}/(q r^4)", Crinf + M (M qq + (2 qq - 1) (r - 2 M))/(qq r^4), ass];
Pd = Module[{d = qq - 2, x = r/M - 2}, (2 d^2 + 2 d) x^3 + (9 d^2 + 24 d + 12) x^2 + (12 d^2 + 48 d + 48) x + 4 d^2 + 16 d + 16];
chk["prop:certificates(1): C_phi = -M P/{q r^5 [2M + (q-1) r]} with the given P/M^3", Cpinf + M (M^3 Pd)/(qq r^5 (2 M + (qq - 1) r)), ass];
C2 = Ckk[2]; R = r/M; x = R - 2;
Pr2 = x^4 + (12 qq - 4) x^3 + (72 qq - 36) x^2 + (112 qq - 32) x + 32 qq;
Pp2 = (3 qq - 5) (qq - 1) x^4 + (20 qq^2 - 37 qq + 14) x^3 + (48 qq^2 - 42 qq) x^2 + 48 qq^2 x + 16 qq^2;
chk["prop:certificates(2): b = 2, C_r = -P_r/[M^2 q R^5 (R+2)^2]", (C2 /. al -> 0) + Pr2/(M^2 qq R^5 (R + 2)^2), ass];
chk["prop:certificates(2): b = 2, C_phi = -P_phi/[M^2 q R^5 (R+2)((q-1)R+2)]", (C2 /. al -> Pi/2) + Pp2/(M^2 qq R^5 (R + 2) ((qq - 1) R + 2)), ass];
C53 = Ckk[5/3];
Qr = 4 (2 - qq) x^4 + (37 qq - 3) x^3 + (390 qq - 210) x^2 + (700 qq - 200) x + 200 qq;
Qp = (7 qq - 11) (qq - 1) x^4 + (47 qq^2 - 87 qq + 34) x^3 + (114 qq^2 - 102 qq) x^2 + 116 qq^2 x + 40 qq^2;
chk["prop:certificates(3): b = 5/3, C_r = -Q_r/[4 M^2 q R^5 (R+3)^2]", (C53 /. al -> 0) + Qr/(4 M^2 qq R^5 (R + 3)^2), ass];
chk["prop:certificates(3): b = 5/3, C_phi = -Q_phi/[2 M^2 q R^5 (R+3)((q-1)R+2)]", (C53 /. al -> Pi/2) + Qp/(2 M^2 qq R^5 (R + 3) ((qq - 1) R + 2)), ass];
(* case 4: b = q constant: Ric = K_t; numerator of -K_t *)
Module[{Cq, num, d, xx},
  Cq = Ckk[qs] /. q[eta] -> qs;
  chk["prop:certificates(4): b = q constant => Ric(k,k) isotropic (independent of alpha)", D[Cq, al], M > 0 && r > 2 M && qs > 1];
  num = Together[-(Cq /. al -> 0) /. {qs -> 3/2 + d, r -> M (2 + xx)}];
  Print["   -K_t = ", Factor[num]];
  chk["prop:certificates(4): numerator of -K_t proportional to the given polynomial (ratio constant in x, d and positive)",
    Module[{P4 = 2 d^2 xx^3 + 9 d^2 xx^2 + 12 d^2 xx + 4 d^2 + d xx^3 + 16 d xx^2 + 34 d xx + 12 d + 23/4 xx^2 + 24 xx + 9, ratio},
      ratio = Factor[Numerator[num]/P4];
      Print["   ratio numerator/P = ", ratio, "   denominator = ", Factor[Denominator[num]]];
      Boole[!FreeQ[ratio, xx] || !FreeQ[ratio, d]]], True]];
chk["prop:certificates: 20q^2 - 37q + 14 >= 0 for q >= 5/3", Reduce[Exists[qv, qv >= 5/3 && 20 qv^2 - 37 qv + 14 < 0], Reals] === False, True];
chk["prop:certificates: 47q^2 - 87q + 34 >= 0 for 11/7 <= q", Reduce[Exists[qv, qv >= 11/7 && 47 qv^2 - 87 qv + 34 < 0], Reals] === False, True];
chk["prop:certificates(3): 4(2-q) >= 0, 37q-3, 390q-210, 700q-200 > 0 on [11/7, 2]", Reduce[Exists[qv, 11/7 <= qv <= 2 && (2 - qv < 0 || 390 qv - 210 <= 0 || (7 qv - 11) (qv - 1) < 0)], Reals] === False, True];
chk["prop:certificates(2): (3q-5)(q-1) >= 0 and 48q^2 - 42q > 0 for q >= 5/3", Reduce[Exists[qv, qv >= 5/3 && ((3 qv - 5) (qv - 1) < 0 || 48 qv^2 - 42 qv <= 0)], Reals] === False, True];
(* rmk:frames-sharp *)
chk["rmk:frames-sharp: lim r^3 C_phi (b = oo) = -2M(q-2)/q", Limit[r^3 Cpinf, r -> Infinity] + 2 M (qq - 2)/qq, ass];
(* lem:hadamard *)
Module[{h0 = DiagonalMatrix[{1/f^2, r^2/f}], K},
  K = Simplify[(ricci[h0, {r, ph}][[1, 1]])/h0[[1, 1]], M > 0 && r > 2 M];
  chk["lem:hadamard: K_{h0} = M(3M - 2r)/r^4", K - M (3 M - 2 r)/r^4, M > 0 && r > 2 M]];
chk["lem:hadamard: r^2/f - 27M^2 = (r - 3M)^2 (r + 6M)/(r - 2M)", r^2/f - 27 M^2 - (r - 3 M)^2 (r + 6 M)/(r - 2 M), M > 0 && r > 2 M];
chk["lem:hadamard: sqrt(N_b)/f ~ 2M/(r - 2M) as r -> 2M", Limit[(Sqrt[1/(1 - f/bb)]/f) (r - 2 M)/(2 M), r -> 2 M, Direction -> "FromAbove"] - 1, M > 0 && bb > 1];
(* family of frames *)
Print["== Supplement: family of frames =="];
Module[{Cb, Prb, Ppb, PrF, PpF, xx},
  Cb = Simplify[Ckk[bb] /. {q[eta] -> qs, M -> 1}, bb > 1 && qs > 1 && r > 2];
  PrF = (bb - 1) (2 bb qs - bb - 4 qs + 3) xx^4 + (13 bb^2 qs - 6 bb^2 - 24 bb qs + 11 bb + 8 qs - 2) xx^3 + 6 bb (5 bb qs - 2 bb - 4 qs + 1) xx^2 + 4 bb^2 (7 qs - 2) xx + 8 bb^2 qs;
  PpF = (qs - 1) (2 bb qs - 4 bb - qs + 3) xx^4 + (13 bb qs^2 - 24 bb qs + 8 bb - 6 qs^2 + 11 qs - 2) xx^3 + 6 qs (5 bb qs - 4 bb - 2 qs + 1) xx^2 + 4 qs^2 (7 bb - 2) xx + 8 bb qs^2;
  chk["sec:family: C_r = -P_r/[q r^5 (2 + (b-1) r)^2] (M = 1, r = 2 + x)", ((Cb /. al -> 0) + PrF/(qs r^5 (2 + (bb - 1) r)^2)) /. xx -> r - 2, bb > 1 && qs > 1 && r > 2];
  chk["sec:family: C_phi = -P_phi/[q r^5 (2 + (b-1) r)(2 + (q-1) r)]", ((Cb /. al -> Pi/2) + PpF/(qs r^5 (2 + (bb - 1) r) (2 + (qs - 1) r))) /. xx -> r - 2, bb > 1 && qs > 1 && r > 2];
  chk["prop:family: coefficients of x^0..x^3 positive for b, q >= 3/2",
    And @@ Flatten[Table[Reduce[Exists[{bv, qv}, bv >= 3/2 && qv >= 3/2 && (Coefficient[P, xx, j] /. {bb -> bv, qs -> qv}) <= 0], Reals] === False, {P, {PrF, PpF}}, {j, 0, 3}]], True];
  chk["prop:family: radial coefficient of x^3 in (u, v) = 13u^2 v + 27/2 u^2 + 15uv + 31/2 u + 5/4 v + 23/8",
    Expand[Coefficient[PrF, xx, 3] /. {bb -> 3/2 + u, qs -> 3/2 + v}] - (13 u^2 v + 27/2 u^2 + 15 u v + 31/2 u + 5/4 v + 23/8)];
  chk["prop:family: radial x^4 = (b-1)((2b-4)q + 3 - b), tangential x^4 = (q-1)((2b-1)q + 3 - 4b)",
    {Coefficient[PrF, xx, 4] - (bb - 1) ((2 bb - 4) qs + 3 - bb), Coefficient[PpF, xx, 4] - (qs - 1) ((2 bb - 1) qs + 3 - 4 bb)} // Total[Abs[#]] &];
  chk["prop:family: at q = q_-(b), (2b-4)q + 3 - b = 3(b-1)(2b-3)/(2b-1)", ((2 bb - 4) qs + 3 - bb /. qs -> (4 bb - 3)/(2 bb - 1)) - 3 (bb - 1) (2 bb - 3)/(2 bb - 1)];
  chk["prop:family: lim r^3 C_r = -(3 - b - 4q + 2bq)/((b-1)q)", Limit[r^3 (Cb /. al -> 0), r -> Infinity] + (3 - bb - 4 qs + 2 bb qs)/((bb - 1) qs), bb > 3/2 && qs > 1];
  chk["prop:family: lim r^3 C_phi = -(3 - 4b - q + 2bq)/((b-1)q)", Limit[r^3 (Cb /. al -> Pi/2), r -> Infinity] + (3 - 4 bb - qs + 2 bb qs)/((bb - 1) qs), bb > 3/2 && qs > 1]];
qm[b_] := (4 b - 3)/(2 b - 1); qp[b_] := (3 - b)/(4 - 2 b);
chk["q_-'(b) = 2/(2b-1)^2, q_+'(b) = 2/(4-2b)^2", {D[qm[b], b] - 2/(2 b - 1)^2, D[qp[b], b] - 2/(4 - 2 b)^2} // Total[Abs[#]] &];
chk["best frame: q_-(b_*) = q_min with b_* = (3 - q_min)/(4 - 2 q_min)", qm[(3 - qn)/(4 - 2 qn)] - qn, 3/2 < qn < 5/3];
chk["q_+(b_*) = (9 - 5 q_min)/(10 - 6 q_min)", qp[(3 - qn)/(4 - 2 qn)] - (9 - 5 qn)/(10 - 6 qn), 3/2 < qn < 5/3];
chk["examples: q_min = 1.55 -> 1.786, 1.6 -> 2.5, 1.65 -> 7.5",
  Max[Abs[N[(9 - 5 # )/(10 - 6 #) & /@ {155/100, 16/10, 165/100}] - {1.786, 2.5, 7.5}]] < 0.0006, True];
summary[];
