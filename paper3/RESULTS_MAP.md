# Results map

Every computer-assisted step of the paper and of its supplementary material, in the order of the text, with the script
that checks it and the log of our run. Statements are named by their LaTeX labels (those of the supplementary material
in the second table). Paths are relative to the folder `paper3/`; logs are in `results/`.

Kind of check:

- **E**: exact symbolic computation (identities, `Reduce`/`Resolve` decisions, rational arithmetic);
- **C**: exact certificate (Bernstein subdivision in rational arithmetic), checked by an independent verifier;
- **I**: outward-rounded interval arithmetic (`mpmath.iv`);
- **N**: floating-point or high-precision numerical computation, without an error bound: the diagnostics of the text,
  labelled there "not certified" or "floating-point", and, inside the symbolic scripts, the checks marked `[N]`
  (numerical integrals, extremes, grids and samples). They are not certificates.

Steps that the text proves by hand without a computer are not listed. Libraries imported by the listed scripts (for
example `interval_arithmetic/jets_nested.py`, `certificate_generators/search.py`) are described in the README of their
folder. Where a step is also checked by an independent implementation, the second one is listed in the last section.

## Main text

| Section | Statement | What is checked | Script | Log (`results/`) | Kind |
|---|---|---|---|---|---|
| The rail as the null cone | eq:threading, eq:h, eq:equatorial-data, lem:gopt-basic, prop:rail-null | the auxiliary metric and its null cone | `symbolic/aux_metric_rail_cone.wls` | `symbolic_aux_metric_rail_cone.txt` | E |
| | prop:rail-null, prop:extremals | rail = null cone in 3+1 Kerr; the Hamiltonian H_eta of Paper II is its support function (the null covector of the proof of prop:extremals) | `python_checks/auxiliary_metric_null_cone.py` | `python_auxiliary_metric_null_cone.txt` | E |
| Causal tools | global hyperbolicity of the auxiliary metric, rmk:freezing | reduction at a = 0, N >= 1, speed bound, radial completeness of h0, the freezing front; a first-arrival example | `python_checks/global_hyperbolicity.py` | `python_global_hyperbolicity.txt` | E (+ N example) |
| Conformal Schwarzschild; the swept angle | lem:h0-complete, rmk:freezing, eq:tk-alpha, lem:tk-sturm, lem:tk-direction, lem:no-return, thm:angle, lem:radius, eq:fermat-plane | screen comparison, direction equation, radius bound, no return beyond 3M | `symbolic/swept_angle_conformal_schwarzschild.wls` | `symbolic_swept_angle_conformal_schwarzschild.txt` | E |
| | lem:radius, lem:tk-direction, lem:no-return, thm:angle | the elementary reasons written in the text in place of `Reduce` decisions | `symbolic/hand_reductions.wls` | `symbolic_hand_reductions.txt` | E |
| | rmk:sharp32, prop:sharp32 | sharpness of q = 3/2: certificate K = 465 and its constants | `symbolic/swept_angle_conformal_schwarzschild.wls`; constants also in the supplement (below) | `symbolic_swept_angle_conformal_schwarzschild.txt` | E |
| | remark after prop:sharp32 | in-plane conjugate points before pi (Δφ = 0.9451π at r = 530.09M; Φ = 0.949π at 1280M; earlier searches at 20M) | `diagnostics/conjugate_case_q6_5.py`, `conjugate_symmetric_asymptotics.py`, `arrival_ray_count_below_three_halves.py` | `diagnostics_*.txt` | N |
| | Figure 2 | radius bound R(q_min) (closed form) and markers of numerical conjugate points | `figures/make_fig_radius.py`; markers from `diagnostics/conjugate_threshold_scan.py`, `conjugate_threshold_near_three_halves.py` (check: `conjugate_envelope_check.py`) | `figures/fig2_radius.pdf`, `diagnostics_conjugate_*.txt` | E (curve), N (markers) |
| | swept angle for variable factors (not quoted in the text) | root searches for four profiles of A (q from 2 down to 1.512) at r0 = 6M, 40 sampled angles δ in [0.02π, 0.98π]: exactly one same-sense ray at every sampled angle for r_T = 6M and 10M; for r_T = 4M two profile–angle pairs where the search found none (0.438π, constant A; 0.414π, rising A) | `diagnostics/arrival_ray_count.py` | `diagnostics_arrival_ray_count.txt` | N |
| Persistence; ingoing Vaidya, arrival branch | eq:vaidya-costs, cor:vaidya, eq:vaidya-G, thm:vaidya-arrival, thm:vaidya-sturm, eq:vaidya-alpha, lem:vaidya-direction, lem:vaidya-dec-boundary, lem:vaidya-dec-sturm, prop:vaidya-dec-count | costs, screen comparison, direction equation, decreasing mass | `symbolic/vaidya_arrival_branch.wls` | `symbolic_vaidya_arrival_branch.txt` | E |
| | lem:vaidya-dec-sturm | elementary reason in place of a decision procedure | `symbolic/hand_reductions.wls` | `symbolic_hand_reductions.txt` | E |
| | lem:vaidya-no-return, thm:vaidya-radius | below 3/2: modes, absence of return, the threshold mu_*(q) = min{1/(27q), (q-1)/q} with maximum 1/28, the chord bound with 12f + 10 - 9q > 1/2 | `symbolic/vaidya_arrival_low_energy.wls`, `symbolic/vaidya_arrival_chord.wls` | `symbolic_vaidya_arrival_low_energy.txt`, `symbolic_vaidya_arrival_chord.txt` | E |
| | thm:vaidya-antipodal-cut (illustration after the theorem) | root searches, Ehat² = 2, r0 = 6, 60 sampled lifted angles in [0.02π, 1.9π], constant mass, m = 1 + tanh(v/5)/2 and a jump from 1 to 2: exactly one same-sense ray at every sampled angle for r_T ∈ {6, 10}, as stated in the text; the extra runs at r_T = 4 (not quoted) miss rays at some angles | `diagnostics/vaidya_ray_count.py` (library `vaidya_arrival_cut.py`) | `diagnostics_vaidya_ray_count.txt` | N |
| Conformal Kerr: a barrier | eq:excursion, lem:worst, margins 2C(2M+, 6M) | cost of an excursion | `symbolic/conformal_kerr_excursion.wls` | `symbolic_conformal_kerr_excursion.txt` | E + N (numerical values of the margins) |
| | thm:kerr | \|beta\|/sqrt(Q0), lapse-shift form; the margins 16.7M, 12.9M, 12.4M quoted after the theorem; random test of the excursion inequality | `python_checks/kerr_barrier.py` | `python_kerr_barrier.txt` | E + N |
| | lem:worst | the polar minimum of the cost density is at the equator, so the 3D margins equal the equatorial ones | `python_checks/kerr_3d_excursion.py` | `python_kerr_3d_excursion.txt` | E + N |
| Conformal Kerr on the equator | eq:kerr-G, lem:kerr-G, lem:kerr-Z, prop:kerr-class | the extended metric, its structure, bands, the collar a*(q) | `symbolic/kerr_equator_structure.wls`; `symbolic/hand_reductions.wls` (lem:kerr-G) | `symbolic_kerr_equator_structure.txt` | E + N (decimal values of the band ends and of the collar) |
| | prop:kerr-Kinf | U, V, J of the representative h = G/N; asymptotics; a numerical sample | `symbolic/kerr_equator_curvature.wls` | `symbolic_kerr_equator_curvature.txt` | E (+ N sample) |
| | prop:kerr-Kinf | Ric_h(k,k) <= 0 for q >= 5/3, every \|a\| <= M: certificate of group 2 | `certificates/group2_Kinf.json`, verifier `certificates/verify_certificate.py` | `certificates_verify_all.txt` | C |
| | prop:kerr-Kvar (1) | R_stat < 0 (q >= 2, \|a\| <= 33/50): group 1 (the certificates give -n_U, -n_V, P_3, P_6 >= 0 on the box; the strict sign for 0 < X < 1 follows because the Bernstein coefficients of the slices of index 1 of -n_U and -n_V and of index 2 of P_3 are all positive, with minima 4, 2 and 56733/27500, so U < 0 and V < -\|J\|, and R_stat <= U(1-\|s\|)^2 < 0 for \|s\| < 1, R_stat <= V + \|J\| < 0 for \|s\| = 1); \|D_1\| <= (5/4)(-R_stat), \|D_2\| <= (1/3)(-R_stat) (q >= 5/2): group 3 | `certificates/group1_*.json`, `certificates/group3_*.json` | `certificates_verify_all.txt` | C |
| | prop:kerr-Kvar (2) | the three inequalities on the spin bands B1–B4: group 4 | `certificates/group4_*.json` | `certificates_verify_all.txt` | C |
| | radial confinement in case (2) of prop:kerr-Kvar | signs of the denominators, -alpha > 0, alpha^2 - beta^2 N Delta > 0: group 5 | `certificates/group5_*.json` | `certificates_verify_all.txt` | C |
| | the same certificates | geometry → published polynomials (identification up to positive factors) | `certificate_generators/bridge.py`; `symbolic/bridge_group4_group5.wls` with `symbolic/check_factor_signs.py` | `certificate_generators_bridge.txt`, `bridge_group4_group5.txt`, `bridge_group4_factor_signs.txt` | E |
| The proper-time branch | lem:tau-density, thm:dive, prop:completion, thm:antipodal, lem:dips, prop:certified (kernel) | density, the dive, completion, antipodal targets | `symbolic/tau_branch_constant_factor.wls` | `symbolic_tau_branch_constant_factor.txt` | E |
| | prop:certified, tab:deltac | for the 20 grid cases: the edges, unimodality of Delta (labels P, C, N and the criterion `[PC]*[CN]*`), the enclosure of delta_c/pi | `interval_arithmetic/edges.py`, `unimodality.py`, `critical_separation.py` | `interval_edges.txt`, `interval_unimodality.txt`, `interval_critical_separation.txt`; leaves in `interval_arithmetic/leaves_unimodality/` | I |
| | prop:maxwell-angle | interior Maxwell angle at q = 1.001, r0 = 1500M: delta*/pi in [0.677106398319, 0.677112764518] | `interval_arithmetic/maxwell_angle.py` | `interval_maxwell_angle.txt`; leaves in `interval_arithmetic/leaves_maxwell/` | I |
| | tab:deltac, prop:maxwell-angle | high-precision recomputation of the 20 cells, of delta* and of b_inf(0), B_V, 4.62260 | `symbolic/tau_branch_numerical_recheck.wls` | `symbolic_tau_branch_numerical_recheck.txt` | N |
| | thm:maxwell | q = 1.01, r0 = 10⁴M: two maxima of Delta, the equal-area rule | `diagnostics/maxwell_rule_scan.py`, `maxwell_bridge_scan.py`, `maxwell_transition_scan.py` | `diagnostics_maxwell_*.txt` | N |
| Time-dependent factor | the proper-time metric h_tau | equals the metric of Perlick 1991, Prop. 3.1 (and differs from that of Prop. 3.2) | `symbolic/tau_metric_perlick.wls` | `symbolic_tau_metric_perlick.txt` | E |
| | lem:radial-red, prop:tau-rail, prop:tau-pmp, eq:Ra-barrier, thm:tau-P, thm:escape | radial reduction, rails, extremals, the barrier, escape price, perturbative selection | `symbolic/tau_branch_variable_factor_vaidya.wls` | `symbolic_tau_branch_variable_factor_vaidya.txt` | E |
| | auxiliary bounds | dρ/dR > 1, max √f/r, monotonicity of R_A in A | `symbolic/tau_branch_general_bounds.wls` | `symbolic_tau_branch_general_bounds.txt` | E |
| | rmk:tau-var-num | shooting against direct optimisation (10⁻⁶–10⁻⁵), crossing near eta_d = 8.03, fold near eta_d = 10.35, B = 3.938 at eta_d = 15, the semi-Lagrangian programme | `diagnostics/tau_var_numbers.py` and the `diagnostics/tau_var_*.py` scripts it summarises | `diagnostics_tau_var_*.txt` | N |
| | thm:tau-var-pulse, cor:tau-var-rho, prop:tau-var-delta | the pulse at a small separation (details in S, below) | `symbolic/pulse_drop_examples.wls`, `symbolic/pulse_drop_constants.wls` | `symbolic_pulse_drop_examples.txt`, `symbolic_pulse_drop_constants.txt` | E (+ N marked) |
| | lem:tau-var-peri, thm:tau-var-drop, thm:tau-var-late | single periastron, the frozen problem at A = 1/2, the metric on the domain, the canonical field in the layer, late drops | `symbolic/drop_large_separation.wls` | `symbolic_drop_large_separation.txt` | E: the bounds used in the proof (closed form of R_1/2, metric and Christoffel symbols on D, the canonical field: \|\|G\|\| < 1, \|\|D_Y G\|\| < 124); N marked: confirmations (dip angle, shallow periastron, NMaximize) |
| | thm:tau-var-late | the constants alpha = g/(4π + 8) and r_s | `symbolic/vaidya_tau_apoapsis_and_late_constants.wls` | `symbolic_vaidya_tau_apoapsis_and_late_constants.txt` | E |
| Equator with spin | lem:tau-eq-ring, prop:tau-eq-price, lem:tau-eq-kin | the ring, costs and clocks per radian, singular arc, Clairaut constant, costate dynamics | `symbolic/tau_equator_angular_motion.wls` (prints the quantities), `symbolic/tau_equator_identities.wls` | `symbolic_tau_equator_angular_motion.txt`, `symbolic_tau_equator_identities.txt` | E |
| | thm:single-max, the warp of the equatorial problem | thm:single-max in the warped metric | `symbolic/tau_equator_warp.wls` | `symbolic_tau_equator_warp.txt` | E |
| | eq:tau-eq-hjb | HJB equation and its covariant form | `symbolic/tau_equator_hjb.wls` | `symbolic_tau_equator_hjb.txt` | E |
| | thm:tau-eq-maxwell, Figure 7 | delta_ring, Phi_max, delta_c; minimisers and lengths of the figure | `diagnostics/equator_ring_frozen.py`; `figures/make_fig_tau_ring.py` | `diagnostics_equator_ring_frozen.txt`, `figures/fig7_tau_ring.pdf` | N |
| | rmk:tau-eq-example | one accepted solution at each of twelve tau in [0.12, 0.20] in each sense | `diagnostics/tau_eq_free_arrival_shooting.py` | `diagnostics_tau_eq_free_arrival_shooting.txt` | N |
| | prop:tau-eq-loc, prop:tau-eq-late, lem:tau-eq-triang, thm:tau-eq-unique, cor:tau-eq-sector, thm:tau-eq-sector | the example: steps of the proofs (details in S, below) | `symbolic/tau_equator_uniqueness_sector.wls`, `symbolic/tau_equator_rational_constants.wls` | `symbolic_tau_equator_uniqueness_sector.txt`, `symbolic_tau_equator_rational_constants.txt` | E (+ N marked) |
| | thm:tau-eq-unique | the floating-point checks quoted after the theorem: D_{2θ}/T >= 0.52, Γ' >= 0.10 | `diagnostics/tau_eq_uniqueness_checks.py` (guides `tau_eq_arrival_shooting.py`, `tau_eq_injectivity_scan.py`) | `diagnostics_tau_eq_*.txt` | N |
| | thm:tau-eq-sector | the root-finding quoted after the theorem: tau₊* = 0.160195, tau₋* = 0.170226, T_s' - 1 = -1.10, -1.09 | `diagnostics/tau_eq_sector_change_checks.py` | `diagnostics_tau_eq_sector_change_checks.txt` | N |
| Proper-time branch of ingoing Vaidya | lem:vaidya-rail, thm:vaidya-P, lem:vaidya-calib | rail, cost and velocity body; selection near constant mass; calibration | `symbolic/tau_branch_variable_factor_vaidya.wls` | `symbolic_tau_branch_variable_factor_vaidya.txt` | E |
| | lem:vaidya-rail (velocity body), lem:vaidya-turn | structure of the extremals of the Vaidya proper-time problem, turning points, kinematics of the velocity body | `symbolic/vaidya_tau_pmp_structure.wls` (with a sampled diagnostic), `vaidya_tau_turning_points.wls`, `vaidya_tau_kinematics.wls` | `symbolic_vaidya_tau_*.txt` | E (+ N sample) |
| | lem:vaidya-turn, thm:vaidya-late | the local apoapsis with r̈ = -23/54 | `symbolic/vaidya_tau_apoapsis_and_late_constants.wls` | `symbolic_vaidya_tau_apoapsis_and_late_constants.txt` | E |
| | thm:vaidya-escape, thm:vaidya-exist-inf | b_inf(0) = 2.6014, lower bound 4.62260, B_V = 5.21295, crossing at 0.522π; existence below the price | `diagnostics/vaidya_tau_escape_price_bounds.py`, `vaidya_tau_verification_function.py`, `vaidya_tau_calibration_identity.py`, `vaidya_tau_existence_no_deadline.py` | `diagnostics_vaidya_tau_*.txt` | N |
| | rmk:vaidya-var-example (thm:vaidya-early with h_* = 1/1000, thm:vaidya-late with D = 91, thm:vaidya-contact) | delta_c < 39/25 at mass 3/2, the early window, late accretion with D = 91 | `riemann_sums/vaidya_example_constants.py` (I), `symbolic/vaidya_example_*.wls` (E, with N confirmations; see the supplementary material, constants of the Vaidya example) | `riemann_sums_vaidya_example_constants.txt`, `symbolic_vaidya_example_*.txt` | I + E |
| | rmk:vaidya-var-example (numerical part) | step of duration 1/2: crossing near d = 0.6, frozen minimiser for d >= 14, deep family above the price | `diagnostics/vaidya_step_continuation.py`, `vaidya_step_branches.py` | `diagnostics_vaidya_step_*.txt` | N |

## Supplementary material

| Section | Statement | What is checked | Script | Log (`results/`) | Kind |
|---|---|---|---|---|---|
| Static conformal frames | lem:static-frame, prop:certificates (cases 1–4), rmk:frames-sharp, lem:hadamard, sec:family, prop:family | the frames, their certificates, Hadamard slices, the best frame | `symbolic/supplement_static_frames.wls` | `symbolic_supplement_static_frames.txt` | E + N (illustrative values of the frames) |
| | prop:certificates | no time derivative of A in the frames; direct Ricci of the frames b = 2 and 5/3; the frozen slice (optical Gauss curvature of Paper II) | `symbolic/supplement_inplane_conjugates.wls`, `supplement_bands_direct_ricci.wls`, `supplement_frozen_slice.wls` | `symbolic_supplement_*.txt` | E |
| | eq:ric-static, prop:family | the exact sign over the whole family (b, q) | `symbolic/supplement_frame_family.wls`, `supplement_frame_family_positivity.wls` | `symbolic_supplement_frame_family*.txt` | E |
| | Figure (bands of the frames) | bands in (q_min, q_max) | `figures/make_fig_bands_supplement.py` | `figures/fig2_bands_supplement.pdf` | E |
| Constants in the sharpness example | prop:sharp32 | K = 465, the bounds 0.1832 and -1.922 for D Phi_b, rho(8200)^2 > D0^2, and the elementary reasons of the main text | `independent_checks/sharp32_constants.py` | `independent_sharp32_constants.txt` | E (rational intervals) |
| | the shell of S10, S6 | elementary reasons in place of decision procedures | `symbolic/hand_reductions.wls` | `symbolic_hand_reductions.txt` | E |
| Curvature certificates on the equator | tab:kerr-cert | the 33 certificates (see the main text, prop:kerr-Kinf and prop:kerr-Kvar) | `certificates/` | `certificates_verify_all.txt` | C |
| The certified critical separations and the Maxwell example | the certified grid; the interior Maxwell angle | as prop:certified and prop:maxwell-angle above | `interval_arithmetic/` | `interval_*.txt` | I |
| The examples with a time-dependent factor | the pulse at a small separation | reduction, constants, chain, sector change tau* in (10/61, 4000/24253); constants 25/97, 607/5000, 1302/10⁴, 1110/10⁴, 61/500, 19/32, the cubic drop, margins | `symbolic/pulse_drop_examples.wls`, `symbolic/pulse_drop_constants.wls` | `symbolic_pulse_drop_examples.txt`, `symbolic_pulse_drop_constants.txt` | E |
| | the drop at a large separation | as thm:tau-var-drop above | `symbolic/drop_large_separation.wls` | `symbolic_drop_large_separation.txt` | E: the bounds used in the proof (closed form of R_1/2, metric and Christoffel symbols on D, the canonical field: \|\|G\|\| < 1, \|\|D_Y G\|\| < 124); N marked: confirmations (dip angle, shallow periastron, NMaximize) |
| The example on the equator with spin | localization, uniqueness, the sector change | identities, extremes on the shell, the rational chain, enclosures of tau±* ([N] table of coefficients marked) | `symbolic/tau_equator_uniqueness_sector.wls`, `symbolic/tau_equator_rational_constants.wls` | `symbolic_tau_equator_uniqueness_sector.txt`, `symbolic_tau_equator_rational_constants.txt` | E (+ N marked) |
| Constants of the Vaidya example | late accretion (D = 91), critical separation at mass 3/2, early window | interval enclosures of Δ(29/10), Δ(297/100), ℓ(297/100), 2R(4), Λ(11/5) and of the competitors | `riemann_sums/vaidya_example_constants.py` | `riemann_sums_vaidya_example_constants.txt`, `vaidya_example_constants.json` | I |
| | | closed forms Λ(11/5), 2R(4), 2R(6); genus 2 of the dip integrals; the banded and rational competitors (length < 91/(11√2)); 11√2 · 5839/1000 < 91; the discriminant | `symbolic/vaidya_example_closed_forms.wls`, `vaidya_example_genus_competitor.wls`, `vaidya_example_discriminant.wls`, `vaidya_example_early_window.wls` | `symbolic_vaidya_example_*.txt` | E (+ N marked) |
| | | the same constants compared at 30–40 digits (`N[...]`): confirmations, the certified enclosures are those of the interval program and of the independent check | `symbolic/vaidya_example_numbers.wls`, `vaidya_example_rational_competitor.wls` | `symbolic_vaidya_example_numbers.txt`, `symbolic_vaidya_example_rational_competitor.txt` | N |
| | the angle bound with arccos √(125/392) | Δ(29/10) > 39/25 > Δ(297/100), ℓ(297/100) > 2R(4) | `symbolic/vaidya_example_early_window.wls` | `symbolic_vaidya_example_early_window.txt` | E + N (marked) |

## Independent second implementations

| Step | Main implementation | Independent check | Log (`results/`) |
|---|---|---|---|
| curvature certificates, groups 1–3 | `certificate_generators/` + `certificates/verify_certificate.py` | `independent_checks/wolfram_bridge_group1.wls`, `wolfram_bridge_group2.wls`, `wolfram_bridge_group3.wls` (geometry rebuilt in Wolfram Language) | `independent_wolfram_bridge_group{1,2,3}.txt` |
| curvature certificates, leaves | `certificates/verify_certificate.py` | `independent_checks/wolfram_leaf_recheck.wls` (every ordinary leaf of `group3_I2_main` and of `group4_B1_I1`) | `independent_wolfram_leaf_recheck_*.txt` |
| unimodality leaves | `interval_arithmetic/unimodality.py` | `independent_checks/unimodality_leaf_archives.py` (reads the archives only) | `independent_unimodality_leaf_archives.txt` |
| constants of the sharpness example | `symbolic/swept_angle_conformal_schwarzschild.wls` | `independent_checks/sharp32_constants.py` | `independent_sharp32_constants.txt` |
| rational constants of the example on the equator with spin | `symbolic/tau_equator_uniqueness_sector.wls`, `tau_equator_rational_constants.wls` | `independent_checks/equator_spin_example_constants.py` | `independent_equator_spin_example_constants.txt` |
| constants of the Vaidya example | `riemann_sums/vaidya_example_constants.py` (cellwise interval enclosure) | `independent_checks/vaidya_example_constants_exact.py` (exact elementary enclosures, monotone Riemann sums) | `independent_vaidya_example_constants_exact.txt` |

What the Wolfram checks cover, and what still rests on the main verifier and on the arguments of the supplement, is
stated in `independent_checks/README.md`.

## Checks of the verifier and of reproducibility

These are not proofs of individual inequalities: they test the verifier and show that the published files are
reproduced from the geometry.

| What is checked | Script | Log (`results/`) |
|---|---|---|
| 15 deliberate alterations of the certificates, all rejected by the verifier | `certificates/tamper_tests.py` | `certificates_tamper_tests.txt` |
| rebuilding every certificate polynomial and every proof tree from the geometry, identical to the published files | `certificate_generators/build_certificates.py --search` | `certificate_generators_build_search.txt` |
| the same without the tree search | `certificate_generators/build_certificates.py` | `certificate_generators_build_check.txt` |
| a second run of the unimodality and Maxwell programs reproduces their archives (byte for byte, apart from the hash of the generator) | `interval_arithmetic/` | `interval_unimodality.txt`, `interval_maxwell_angle.txt` |
