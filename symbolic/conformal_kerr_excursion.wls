#!/usr/bin/env wolframscript
(* Paper III, section on conformal Kerr: cost of an excursion. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Conformal Kerr, barrier =="];
ass = M > 0 && r > 2 M && 0 < aa <= M && 0 <= x <= 1;
Dl = r^2 - 2 M r + aa^2; f = 1 - 2 M/r;
thF = 1 - 2 M aa/(r Sqrt[Dl]);
chk["eq:excursion: P0 = h0_rr on the equator = r^2/(f Delta), c = theta r/sqrt(f Delta)", thF Sqrt[r^2/(f Dl)] - thF r/Sqrt[f Dl], ass];
chk["eq:excursion: theta(2M) = 0", thF /. r -> 2 M, M > 0 && 0 < aa <= M];
chk["eq:excursion: theta -> 1 for r -> oo", Limit[thF, r -> Infinity] - 1, M > 0 && 0 < aa <= M];
dth = Together[D[1 - 2 av/(y Sqrt[y^2 - 2 y + av^2]), y]];
chk["eq:excursion: theta increasing on r > 2M (M = 1)", Reduce[y > 2 && 0 < av <= 1 && dth <= 0, {y, av}, Reals] === False, True];
chk["eq:excursion: C(2M+, r2) finite for a != 0 (numerical integral at a = 1/2, r2 = 6M)",
  Module[{v = NIntegrate[(thF r/Sqrt[f Dl]) /. {M -> 1, aa -> 1/2}, {r, 2, 6}]}, If[NumericQ[v] && v < 10, 0, 1]]];
(* lem:worst *)
S = r^2 + aa^2 - aa^2 x^2; k = 2 M aa r/Sqrt[Dl];
hrr = S/(Dl (1 - 2 M r/S));
bnorm = 2 M aa r x/(S Sqrt[Dl]);
chk["lem:worst: |beta|_h0 = k x/Sigma", bnorm - k x/S, ass];
Fw = (1 - bnorm) Sqrt[hrr];
G = (S - k x)/Sqrt[S - 2 M r];
chk["lem:worst: sqrt(Delta) F = G(x) = (Sigma - k x)/sqrt(Sigma - 2Mr)", Sqrt[Dl] Fw - G, ass && S > 2 M r];
chk["lem:worst: (Sigma - 2Mr)^{3/2} G'(x) = a^2 x(4Mr - Sigma) - k Delta", (S - 2 M r)^(3/2) D[G, x] - (aa^2 x (4 M r - S) - k Dl), ass && S > 2 M r];
chk["lem:worst: 4M^2 Delta - a^2(4M - r)^2 = s(8M^3 + 4M a^2) + s^2(4M^2 - a^2), r = 2M + s",
  (4 M^2 Dl - aa^2 (4 M - r)^2 /. r -> 2 M + s) - (s (8 M^3 + 4 M aa^2) + s^2 (4 M^2 - aa^2))];
chk["lem:worst: on the equator F = c(r)", (Fw /. x -> 1) - thF r/Sqrt[f Dl], ass];
chk["lem:worst: F non-increasing in x = sin(theta) (Reduce, M = 1)",
  Reduce[Exists[{rr, av, xx}, rr > 2 && 0 < av <= 1 && 0 <= xx <= 1 &&
     (aa^2 x (4 M r - S) - k Dl /. {M -> 1, aa -> av, r -> rr, x -> xx}) > 0], Reals] === False, True];
(* numerical margins 2C(2M+, 6M) *)
Scan[Function[{pr}, Module[{av = pr[[1]], mv = pr[[2]], v}, v = 2 NIntegrate[(thF r/Sqrt[f Dl]) /. {M -> 1, aa -> av}, {r, 2, 6}, WorkingPrecision -> 20];
   Print["   2C(2M+, 6M) at a = ", av, ": ", N[v, 6]];
   chk["margin 2C(2M+, 6M) = " <> ToString[mv] <> "M for a = " <> ToString[av, InputForm], Abs[v - mv] < 0.05, True]]],
  {{1/2, 16.7}, {9/10, 12.9}, {99/100, 12.4}}];
summary[];
