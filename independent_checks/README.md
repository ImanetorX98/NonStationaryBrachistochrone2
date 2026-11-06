# Independent checks

Second implementations of some computer-assisted steps, of two kinds.

- **Python, written through OpenAI Codex**, while the main scripts were written through Claude Code (see the paragraph
  *Verification methods and use of AI tools* in the paper). Each check is standalone and imports no other file of
  this repository (apart from the optional progress bar of `../tools/`).
- **Wolfram Language, a second computer algebra system** (`wolfram_*.wls`): the bridge between the curvature
  certificates and the geometry, rebuilt from scratch (own Christoffel symbols, Ricci tensor and frames), and a
  leaf-by-leaf recheck of the certificate trees. They read only the published JSON files and do not use the Python
  chain of `../certificate_generators/`.

**What the Wolfram checks cover.** The bridge scripts, together with `../symbolic/bridge_group4_group5.wls`,
independently reconstruct the geometric identities for the five groups: the published polynomials are the curvature
expressions of the text, multiplied by explicit positive factors, in the declared box maps. Complete acceptance of the
certificates also requires the main verifier `../certificates/verify_certificate.py` (proof trees, radicands, delegated
certificates and their composition identities), the steps that the supplement proves separately (paragraph *What the
verifier does not check*: the implications from the sufficient sets, the covering of region C2 of group 3 by C2A and
C2B, continuity at the boundary points of the boxes), and for group 4 the sign argument of
`../symbolic/check_factor_signs.py`. The Wolfram leaf recheck is an additional check of the ordinary radical leaves of
certificates already accepted by the main verifier; it does not replace the validation of the whole tree, of the
radicand or of the delegations.

| File | Step of the paper | Method | Expected output |
|---|---|---|---|
| `sharp32_constants.py` | constants of the sharpness example below q = 3/2 (K = 465, the bounds 0.1832 and -1.922 for D Phi_b, rho(8200)^2 > D0^2) and the elementary reasons written in the main text in place of decision procedures | rational intervals, Machin's formula for pi, SymPy identities; three negative controls (a slightly stronger bound, or the opposite monotonicity, is shown to be false) | `../results/independent_sharp32_constants.txt`, 34 PASS |
| `equator_spin_example_constants.py` | rational constants of the example on the equator with spin (localization, uniqueness) | exact fractions and SymPy identities | `../results/independent_equator_spin_example_constants.txt`, 56 PASS |
| `vaidya_example_constants_exact.py` | constants of the Vaidya example: late accretion (D = 91), critical separation at mass 3/2, early window | exact elementary enclosures (rational square roots, atanh series with tail); monotone Riemann sums with 16384 intervals | `../results/independent_vaidya_example_constants_exact.txt` and `.json`, 38 PASS |
| `unimodality_leaf_archives.py` | archived leaves of the unimodality certificates | reads the JSON archives only: tiling, signs of the printed bounds, the criterion `[PC]*[CN]*` | 20 archives PASS |
| `wolfram_bridge_group1.wls` | curvature certificates, group 1 (constant factor, frame of gbar_inf) | Ricci from scratch; U, V, J of the text; the four JSON polynomials equal -nU, -nV, P3, P6 in the box map; Bernstein coefficients of the single leaves; the deduction of the sign of Ric(k,k) by `Resolve` | `../results/independent_wolfram_bridge_group1.txt`, 21 PASS |
| `wolfram_bridge_group2.wls` | curvature certificates, group 2 (h = G/N, ZAMO frame) | Ricci from scratch; P_U, P_V, P_J^2 derived from the curvature (polynomials) and equal to the formulas of the text; the six JSON polynomials equal their combinations; the three sufficient sets imply the sign by `Resolve` | `../results/independent_wolfram_bridge_group2.txt`, 21 PASS |
| `wolfram_bridge_group3.wls` | curvature certificates, group 3 (time-dependent factor, jets u', u'') | Ricci from scratch; jet structure and the formulas of R_stat, D1, D2; P0 + P1 sqrt(g) of each main certificate is an explicit positive factor times the inequality of the text, and g = w S | `../results/independent_wolfram_bridge_group3.txt`, 15 PASS |
| `wolfram_leaf_recheck.wls` | ordinary leaves of a certificate of type radical, already accepted by the main verifier | rebuilds the cells from the tree, tensor Bernstein coefficients, the chord and tangent checks of each leaf; the lower bound of the radicand g on a cell is taken as max(min Bernstein coefficient, 0), which relies on g >= 0 as checked by the main verifier; a tampered leaf (P0 - 1/10) is printed for comparison, its rejection is not guaranteed; the last line `RESULT: PASS` (exit status 0) requires the whole tree read, cells tiling the unit box and every checked leaf passing | `../results/independent_wolfram_leaf_recheck_*.txt`; with `N = all` every ordinary leaf plus two repeated ones (the first and the deepest), so 2231 checks for the 2229 leaves of `group3_I2_main` (about an hour) and 50 for the 48 leaves of `group4_B1_I1` (about 80 minutes), all PASS; the tampered leaf is rejected for `group3_I2_main` and still accepted for `group4_B1_I1` |

The main implementations of the same steps are `../riemann_sums/vaidya_example_constants.py` (a different method:
cellwise interval enclosure after the substitution r = t + (4 - t) sin^2 theta), the scripts in `../symbolic/`, and
`../interval_arithmetic/`.

Usage of the Wolfram files (from this folder; pass the output through a pipe, see `../symbolic/README.md`):

```
wolframscript -file wolfram_bridge_group1.wls ../certificates 2>&1 | cat
wolframscript -file wolfram_leaf_recheck.wls ../certificates/group3_I2_main.json 10 | cat    # 10 random leaves
wolframscript -file wolfram_leaf_recheck.wls ../certificates/group3_I2_main.json all | cat   # every leaf
```

The leaf recheck shows a progress bar on standard error (or writes it to the file named by the environment variable
`PROGRESS_FILE`), so only standard output goes through the pipe; the bar does not change the output.
