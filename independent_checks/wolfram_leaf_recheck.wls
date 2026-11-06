#!/usr/bin/env wolframscript
(* Independent Wolfram recheck of leaves of a certificate of type 'radical' (P0 + P1 Sqrt[g] >= 0). Reads only the JSON
   (does not use the Python code of the project): rebuilds the cells from the cuts of the tree, restricts P0, P1, g to
   each cell by an affine substitution, computes the tensor Bernstein coefficients (matrices C(i,j)/C(n,j)) and repeats
   the checks of a radical leaf: envelopes s0^2 <= g0, s1^2 <= g1, sigma > 0; chord and tangent; coefficients >= 0.
   Usage: wolframscript -file wolfram_leaf_recheck.wls CERTIFICATE.json N
   checks the first leaf, the deepest leaf and N random leaves (N = all: every leaf). Control: P0 - 1/10 on the first
   leaf (it must be rejected when the margin is smaller than 1/10; a rejection is not guaranteed otherwise, and the
   control does not enter the result). The last line is RESULT: PASS (exit status 0) only if the whole tree was read,
   the cells tile the unit box and every checked leaf passes; otherwise RESULT: FAIL (exit status 1). This is an
   additional check of certificates already accepted by ../certificates/verify_certificate.py: the lower bound of g on
   a cell is taken as max(min Bernstein coefficient, 0), which relies on g >= 0 as checked by the main verifier. *)
file = $ScriptCommandLine[[2]]; nsample = If[$ScriptCommandLine[[3]] === "all", Infinity, ToExpression[$ScriptCommandLine[[3]]]];
cert = Import[file, "RawJSON"];
vars = cert["variables"]; d = Length[vars]; X = Table[Symbol["v" <> ToString[k]], {k, d}];
rat[s_String] := ToExpression[s];
poly[name_] := Total[(rat[#[[2]]] Times @@ (X^#[[1]])) & /@ cert["polynomials"][name]];
P0 = poly["P0"]; P1 = poly["P1"]; g = poly["g"]; P1g = Expand[P1 g];
degs = Table[Max[Exponent[#, X[[k]]] & /@ {P0, P1, P1g, g}], {k, d}];
Print["certificate ", cert["name"], "  degrees ", degs];
(* walk of the tree: cells rebuilt from the cuts *)
tree = cert["tree"]; pos = 0; leaves = {}; dels = {};
walk[cell_] := Module[{tok, ax, a, b, m},
  pos++; tok = tree[[pos]];
  Which[tok[[1]] == "s", ax = tok[[2]] + 1; {a, b} = cell[[ax]]; m = (a + b)/2;
    walk[ReplacePart[cell, ax -> {a, m}]]; walk[ReplacePart[cell, ax -> {m, b}]],
   tok[[1]] == "l", AppendTo[leaves, {cell, tok}],
   tok[[1]] == "d", AppendTo[dels, {cell, tok}]]];
walk[ConstantArray[{0, 1}, d]];
vol[c_] := Times @@ (#[[2]] - #[[1]] & /@ c);
Print["nodes read ", pos, " of ", Length[tree], "; leaves ", Length[leaves], "; delegations ", Length[dels],
  "; total volume ", Total[vol /@ Join[leaves[[All, 1]], dels[[All, 1]]]]];
(* tensor Bernstein coefficients on the cell *)
bmat[n_] := Table[If[j <= i, Binomial[i, j]/Binomial[n, j], 0], {i, 0, n}, {j, 0, n}];
mats = bmat /@ degs;
bern[p_, cell_] := Module[{t = Table[Symbol["t" <> ToString[k]], {k, d}], q, A},
  q = Expand[p /. Thread[X -> (cell[[All, 1]] + (cell[[All, 2]] - cell[[All, 1]]) t)]];
  A = CoefficientList[q, t, degs + 1];
  Do[A = TensorContract[TensorProduct[mats[[k]], A], {{2, k + 2}}];
     A = Transpose[A, Join[Range[2, k], {1}, Range[k + 1, d]]], {k, d}];
  Flatten[A]];
check[{cell_, tok_}] := Module[{s0 = rat[tok[[2]]], s1 = rat[tok[[3]]], sg = rat[tok[[4]]], bg, g0, g1, Bc, Ac, out},
  bg = bern[g, cell]; g0 = Max[Min[bg], 0]; g1 = Max[bg];
  If[! (s0 >= 0 && s0^2 <= g0 && s1 >= 0 && s1^2 <= g1 && sg > 0), Return[{False, "envelope"}]];
  Bc = If[g1 > g0, (s1 - s0)/(g1 - g0), 0]; Ac = s0 - Bc g0;
  out = Min[bern[P0 + #[[1]] P1 + #[[2]] P1g, cell]] & /@ {{Ac, Bc}, {sg/2, 1/(2 sg)}};
  {AllTrue[out, # >= 0 &], N[out, 3]}];
SeedRandom[20261102];
deep = First[MaximalBy[leaves, -vol[#[[1]]] &]];
sample = If[nsample === Infinity, Join[{First[leaves], deep}, leaves], Join[{First[leaves], deep}, RandomSample[leaves, Min[nsample, Length[leaves]]]]];
(* progress bar, kept out of the log: written to $PROGRESS_FILE if set, otherwise to standard error (at most every 2 s) *)
progFile = Environment["PROGRESS_FILE"]; progLast = 0;
progress[k_, tot_, el_] := Module[{n = Floor[30 k/tot], line},
   If[k < tot && AbsoluteTime[] - progLast < 2, Return[Null]]; progLast = AbsoluteTime[];
   line = "leaves [" <> StringRepeat["#", n] <> StringRepeat["-", 30 - n] <> "] " <> ToString[k] <> "/" <> ToString[tot] <>
     ", " <> ToString[Round[el]] <> " s elapsed, about " <> ToString[Round[el (tot - k)/k]] <> " s left";
   If[StringQ[progFile], Export[progFile, line <> "\n", "Text"], WriteString["stderr", "\r" <> line <> If[k == tot, "\n", ""]]]];
t0 = AbsoluteTime[];
res = Table[With[{r = check[sample[[i]]]}, progress[i, Length[sample], AbsoluteTime[] - t0]; r], {i, Length[sample]}];
Do[Print[If[res[[i, 1]], "PASS ", "FAIL "], "cell ", InputForm[sample[[i, 1]]], "  min coefficient (chord, tangent) ", res[[i, 2]]],
  {i, Length[sample]}];
Print["leaves checked ", Length[sample], ", PASS ", Count[res[[All, 1]], True], ", time ", Round[AbsoluteTime[] - t0], " s"];
(* control: the same leaf with P0 - 1/10 *)
P0 = P0 - 1/10; Print["tampered P0 - 1/10 on the first leaf, still accepted: ", First[check[First[leaves]]]];
okTree = pos === Length[tree] && Total[vol /@ Join[leaves[[All, 1]], dels[[All, 1]]]] === 1;
okLeaves = AllTrue[res[[All, 1]], TrueQ];
Print["RESULT: ", If[okTree && okLeaves, "PASS", "FAIL"], " (whole tree read and unit volume: ", okTree,
  "; all checked leaves pass: ", okLeaves, ")"];
Exit[If[okTree && okLeaves, 0, 1]];
