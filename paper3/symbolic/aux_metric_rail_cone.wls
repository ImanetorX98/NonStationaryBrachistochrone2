#!/usr/bin/env wolframscript
(* Paper III, section on the auxiliary metric (sec:aux): threading form of Kerr, auxiliary metric g_opt, rail = null cone. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "common.wl"}]];
Print["== Auxiliary metric: the rail as the null cone of g_opt =="];
ass = M > 0 && r > 2 M && 0 < th < Pi && a \[Element] Reals && Abs[a] < M && A > 0 && E0 > 0;
S = r^2 + a^2 Cos[th]^2; Dl = r^2 - 2 M r + a^2;
gtt = -(1 - 2 M r/S); gtp = -2 M a r Sin[th]^2/S; gpp = (r^2 + a^2 + 2 M a^2 r Sin[th]^2/S) Sin[th]^2;
f = 1 - 2 M r/S; bph = gtp/f;
hpp = gpp + gtp^2/f;
chk["eq:threading: g_Kerr = -f(dt - beta)^2 + h (components tt, t phi, phi phi)",
  {gtt + f, gtp - f bph, gpp - (-f bph^2 + hpp)} // Simplify // Total[#^2] &, ass];
chk["beta = -2 M a r sin^2/(Sigma f) dphi", bph + 2 M a r Sin[th]^2/(S f), ass];
chk["identity f g_phiphi + g_tphi^2 = Delta sin^2", f gpp + gtp^2 - Dl Sin[th]^2, ass];
chk["eq:h: h_phiphi = Delta sin^2/f", hpp - Dl Sin[th]^2/f, ass];
chkFalse["control: h_phiphi != g_phiphi (must fail)", hpp - gpp, ass];
(* eq:equatorial-data *)
eq = th -> Pi/2;
chk["h0_rr on the equator = r^2/(f Delta)", ((S/Dl)/f /. eq) - r^2/((1 - 2 M/r) Dl), ass];
chk["h0_phiphi on the equator = Delta/f^2", ((hpp/f) /. eq) - Dl/(1 - 2 M/r)^2, ass];
chk["beta on the equator = -2 M a/(r f)", (bph /. eq) + 2 M a/(r (1 - 2 M/r)), ass];
(* lem:gopt-basic *)
b2 = bph^2/(hpp/f);
chk["lem:gopt-basic(1): 1 - |beta|^2_h0 = f g_phiphi/(Delta sin^2)", 1 - b2 - f gpp/(Dl Sin[th]^2), ass];
chk["lem:gopt-basic(1): |beta|^2_h0 = 4 M^2 a^2 r^2 sin^2/(Sigma^2 Delta)", b2 - 4 M^2 a^2 r^2 Sin[th]^2/(S^2 Dl), ass];
Gb = {{-1, bph}, {bph, NN hpp/f - bph^2}};
chk["lem:gopt-basic(2): g_opt^{eta eta} = -1 + |beta|^2/N", Inverse[Gb][[1, 1]] - (-1 + b2/NN), ass && NN >= 1];
(* lem:gopt-basic(3): g_opt(v,v) and the future cone *)
chk["lem:gopt-basic(3): g_opt(v,v) = -(v^eta - beta(v))^2 + N |v|^2",
  ({ve, vp}.Gb.{ve, vp}) - (-(ve - bph vp)^2 + NN hpp/f vp^2), ass];
(* prop:rail-null *)
Z = A^2 f/E0^2; NNv = 1/(1 - Z);
w = E0/(A^2 f); hu = f w^2 - A^-2;
chk["prop:rail-null(1): h(u,u)/(f w^2) = 1 - Z = 1/N", hu/(f w^2) - 1/NNv, ass];
lam = E0/(A^2 f wv);
chk["prop:rail-null(2): g(lambda v, lambda v) = -1 with w_v^2 = N h0(v,v)",
  (A^2 lam^2 f (-wv^2 + wv^2/NNv)) + 1, ass && wv > 0];
chk["prop:rail-null: N = n^2 with n = E/sqrt(E^2 - A^2 f)", (E0/Sqrt[E0^2 - A^2 f])^2 - NNv, ass && E0^2 > A^2 f];
summary[];
