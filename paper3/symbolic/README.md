# Symbolic checks (Wolfram Language)

Exact identities, decisions by `Reduce`/`Resolve`/`MinValue`, rational arithmetic, and a few high-precision numerical
recomputations marked `[N]` (not certified), for the steps of the paper and of the supplementary material that were
done or confirmed with a computer algebra system. Each check prints `PASS` or `FAIL` with a label; most labels start
with the LaTeX label of the statement they support (for example `lem:radius`, `prop:sharp32`).

## Running

```
./run_symbolic.sh                 # all scripts; logs in ../results/symbolic_<name>.txt
./run_symbolic.sh hand_reductions # one or more scripts (name without .wls)
./run_symbolic.sh --check         # check the existing logs, without running Wolfram
./run_symbolic.sh --status        # from another terminal: the running script
```

The launcher accepts a script only if wolframscript exits with status 0, the number of `PASS` lines equals the one in
`expected_outcomes.tsv`, the `FAIL` lines are exactly the deliberate negative controls listed in
`expected_failures.txt` (two lines of `hand_reductions.wls`; they must be present), and the log ends with the expected
last line, so that an empty or truncated log is rejected. It exits with status 0 only if every script is accepted.

Requirements: `wolframscript` with Wolfram Engine or Mathematica (tested with Mathematica 13.3.1). The whole set takes
about ten minutes; `bridge_group4_group5.wls` alone (it reads `../certificates`) takes 5–25 minutes and shows a progress
bar on standard error.

**Capture the output through a pipe.** With `wolframscript -file x.wls > log`, printed lines were occasionally lost in
our runs (one script lost one or two lines in 3 of 3 runs), while `wolframscript -file x.wls 2>&1 | cat > log` and
`WolframKernel -script x.wls` were complete in 5 of 5 runs. The computations are not affected, only the capture of the
output; `run_symbolic.sh` uses the pipe and checks the completeness of each log as described above.

## Map: script → statements

| Script | Statements of the paper (main text unless marked S = supplement) | Checks |
|---|---|---|
| `common.wl` | shared functions (`chk`, `chkFalse`, `chkN`, Christoffel, Riemann, Ricci) | — |
| `aux_metric_rail_cone.wls` | auxiliary metric, rail = null cone: eq:threading, eq:h, eq:equatorial-data, lem:gopt-basic, prop:rail-null | 15 |
| `swept_angle_conformal_schwarzschild.wls` | conformal Schwarzschild and the swept angle: lem:h0-complete, rmk:freezing, eq:tk-alpha, lem:tk-sturm, lem:tk-direction, lem:no-return, thm:angle, lem:radius, eq:fermat-plane, rmk:sharp32, prop:sharp32 (certificate K = 465) | 56 |
| `vaidya_arrival_branch.wls` | persistence and the arrival branch of ingoing Vaidya: cor:vaidya, eq:vaidya-G, thm:vaidya-arrival, thm:vaidya-sturm, eq:vaidya-alpha, lem:vaidya-direction, lem:vaidya-dec-boundary, lem:vaidya-dec-sturm, prop:vaidya-dec-count, admissible profile | 49 |
| `vaidya_arrival_low_energy.wls`, `vaidya_arrival_chord.wls` | arrival branch of Vaidya below 3/2: modes, absence of return, the chord bound with 12f + 10 − 9q > 1/2 | 6 + 6 |
| `conformal_kerr_excursion.wls` | conformal Kerr, cost of an excursion: eq:excursion, lem:worst, margins 2C(2M+, 6M) | 14 |
| `kerr_equator_structure.wls` | Kerr equator, extended metric: eq:kerr-G, lem:kerr-G, lem:kerr-Z, prop:kerr-class, bands, collar a*(q) | 34 |
| `kerr_equator_curvature.wls` | Kerr equator, curvature: structure of U, V, J, asymptotics, prop:kerr-Kinf ([N] sample), jets, conformal formula in 2+1 | 16 |
| `tau_branch_constant_factor.wls` | proper-time branch, constant factor: lem:tau-density, thm:dive, prop:completion, thm:antipodal, prop:certified (kernel), lem:dips | 18 |
| `tau_branch_variable_factor_vaidya.wls` | proper-time branch with A(η) and for Vaidya: lem:radial-red, prop:tau-rail, prop:tau-pmp, eq:Ra-barrier, thm:tau-P, thm:escape, lem:vaidya-rail, thm:vaidya-P, lem:vaidya-calib | 34 |
| `tau_branch_numerical_recheck.wls` | [N] recomputation of tab:deltac (20 cells), prop:maxwell-angle, b_∞(0), B_V, 4.62260 | 25 |
| `tau_branch_general_bounds.wls` | dρ/dR > 1, max √f/r, monotonicity of R_A in A | 3 |
| `tau_metric_perlick.wls` | h_τ equals the metric of Perlick 1991, Prop. 3.1 | 2 |
| `tau_equator_angular_motion.wls`, `tau_equator_identities.wls`, `tau_equator_warp.wls`, `tau_equator_hjb.wls` | Thakurta–Kerr equator, proper-time branch: the ring, costs and clocks per radian, singular arc, Clairaut constant, costate dynamics, thm:single-max warp, HJB equation and its covariant form | 0 + 14 + 7 + 9 (the first prints values) |
| `tau_equator_uniqueness_sector.wls` | S10, example on the equator with spin: steps of thm:tau-eq-unique and thm:tau-eq-sector (identities, extremes on the shell, [N] table of coefficients, rational chain, enclosures of τ±*) | 139 |
| `tau_equator_rational_constants.wls` | S10: rational constants and auxiliary claims of the same proofs | 52 |
| `pulse_drop_examples.wls` | S9: pulse at a small separation and drop at a large separation (reduction, constants, chain, sector change τ* ∈ (10/61, 4000/24253)) | 104 |
| `pulse_drop_constants.wls` | S9: constants 25/97, 607/5000, 1302/10⁴, 1110/10⁴, 61/500, 19/32, the cubic drop, margins | 17 |
| `drop_large_separation.wls` | thm:tau-var-drop (S9): single periastron, frozen problem at A = 1/2 (closed form of R_1/2), metric and connection on the domain (exact bounds), canonical field in the layer (exact bounds on G and D_Y G), late drops | 69 |
| `vaidya_tau_pmp_structure.wls`, `vaidya_tau_turning_points.wls`, `vaidya_tau_kinematics.wls`, `vaidya_tau_apoapsis_and_late_constants.wls` | proper-time branch of Vaidya: structure of the PMP, turning points, velocity body, the apoapsis with r̈ = −23/54, α = g/(4π+8) | 5 + 7 + 5 + 11 |
| `vaidya_example_*.wls` (6 files) | S11, constants of the Vaidya example: closed forms Λ(11/5), 2R(4), 2R(6); genus 2 of the dip integrals; rational banded competitor (length < 91/(11√2)); 11√2 · 5839/1000 < 91; discriminant; early window and δ_c at mass 3/2 | 1 + 4 + 4 + 4 + 4 + 8 |
| `hand_reductions.wls` | reasons written by hand in the text in place of `Reduce`/`NMinValue` decisions (lem:radius, lem:tk-direction, lem:no-return, thm:angle, lem:kerr-G, lem:vaidya-dec-sturm, the shell of S10, S6) | 39 + 2 that must fail |
| `supplement_static_frames.wls` | S: lem:static-frame, prop:certificates (cases 1–4), lem:hadamard, rmk:frames-sharp, sec:family, prop:family, best frame | 30 |
| `supplement_inplane_conjugates.wls`, `supplement_bands_direct_ricci.wls`, `supplement_frozen_slice.wls`, `supplement_frame_family.wls`, `supplement_frame_family_positivity.wls` | S, prop:certificates and prop:family: no time derivative of A in the frames, direct Ricci of the frames b = 2 and 5/3, the frozen slice (optical Gauss curvature of Paper II), eq:ric-static, the exact sign over the whole family (b, q) | 8 + printed + 7 + 5 + printed |
| `bridge_group4_group5.wls` + `check_factor_signs.py` | S7, curvature certificates of groups 4 and 5: geometry → published polynomials, and the sign of the common factor F | 47 + 12 |

`[N]` marks checks that rest on grids, sampling, numerical integration or numerical optimisation (high-precision
extremes, integrals and grid bounds): they confirm the stated constants, they are not certificates. Checks labelled
"must fail" or "control" are negative controls: their expected outcome is stated in the label. The independent Wolfram
bridge for the certificates of groups 1–3 and a leaf-by-leaf recheck are in `../independent_checks/`.
