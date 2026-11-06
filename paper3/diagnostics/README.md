# Floating-point diagnostics

The scripts in this folder produce the numbers that the paper labels as numerical observations, "not certified" or
"floating-point". They are **not proofs**: no error bound is attached to them, the searches are finite, and the
optimisations are local. The certified statements of the paper rest on `../certificates/`, `../interval_arithmetic/`,
`../riemann_sums/` and `../symbolic/`.

Each script states in its docstring which passage of the paper it supports and what it computes. All logs are in
`../results/diagnostics_<job>.txt`.

## Running

```
./run_diagnostics.sh --list          # jobs and approximate durations
./run_diagnostics.sh                 # everything, in order (about 7 hours on one core)
./run_diagnostics.sh tau_var_numbers # one or more named jobs
PYTHON="nice -n 10 python3" ./run_diagnostics.sh ...
```

Requirements: Python 3 with numpy, scipy, sympy, mpmath (see `../requirements.txt`); set `PYTHON` if the default
`python3` lacks them. Standard error, which contains
only warnings of numpy/scipy in our runs, goes to `../results/diagnostics_<job>.stderr` when it is not empty.
`tau_var_numbers` reads the log of `tau_var_direct_optimisation`, so it runs after it.

## Map: passage of the paper → script

| Passage of the paper | Numbers | Script(s) |
|---|---|---|
| Remark after the proposition on sharpness at q = 3/2 | Δφ = 0.9451π at r = 530.09M (q = 6/5, r₀ = 10⁴M) | `conjugate_case_q6_5.py` |
| same remark | q = 1.2, r₀ = r_T = 1280M, Φ = 0.949π | `conjugate_symmetric_asymptotics.py` (line (3) at d = 1280) |
| same remark, "earlier searches at 20M" | counts of arrival rays below 3/2 | `arrival_ray_count_below_three_halves.py` |
| Figure 2 (markers) | thresholds of in-plane conjugate points | `conjugate_threshold_scan.py`, `conjugate_threshold_near_three_halves.py`; library `conjugate_points_static.py`; independent check `conjugate_envelope_check.py` |
| Theorem on the swept angle for variable factors (not quoted in the text) | four profiles of A, 40 sampled angles in [0.02π, 0.98π]: one same-sense ray at each for r_T = 6M, 10M; for r_T = 4M two profile–angle pairs without a ray found (0.438π, 0.414π) | `arrival_ray_count.py` |
| Theorem on the Maxwell rule | q = 1.01, r₀ = 10⁴M, two maxima of Δ | `maxwell_rule_scan.py`; scans `maxwell_bridge_scan.py`, `maxwell_transition_scan.py` |
| Theorem on the Maxwell rule with a ring, Figure 7 | δ_ring, Φ_max, δ_c | `equator_ring_frozen.py` |
| Theorem on the antipodal cut (ingoing Vaidya) | 60 sampled lifted angles in [0.02π, 1.9π]: one same-sense ray at each for r_T ∈ {6, 10}, as in the text; the extra runs at r_T = 4 miss rays at some angles | `vaidya_ray_count.py`; library `vaidya_arrival_cut.py` |
| Remark on numerical observations (variable factor) | agreement of shooting and direct optimisation to 10⁻⁶–10⁻⁵; crossing near η_d = 8.03; fold near η_d = 10.35; B = 3.938 at η_d = 15; semi-Lagrangian programme | `tau_var_numbers.py` (the four numbers); `tau_var_shooting.py`, `tau_var_direct_optimisation.py`, `tau_var_drop_time.py`, `tau_var_fold.py`, `tau_var_asymmetric_dive.py`, `tau_var_late_drop.py`, `tau_var_value_function.py`; saturation lemma: `tau_var_saturation.py` |
| Example on the Kerr equator with spin | one accepted solution at each of twelve τ in [0.12, 0.20] in each sense | `tau_eq_free_arrival_shooting.py` |
| Theorem on global uniqueness in the example | D_{2θ}/T ≥ 0.52, Γ′ ≥ 0.10 | `tau_eq_uniqueness_checks.py`; guides `tau_eq_arrival_shooting.py`, `tau_eq_injectivity_scan.py` |
| Theorem on a single transversal sector change | τ₊* = 0.160195, τ₋* = 0.170226, T_s′ − 1 = −1.10, −1.09 | `tau_eq_sector_change_checks.py` |
| Theorem on the escape price (ingoing Vaidya) | b_∞(0) = 2.6014, lower bound 4.62260, B_V = 5.21295, crossing at 0.522π | `vaidya_tau_escape_price_bounds.py`, `vaidya_tau_verification_function.py`, `vaidya_tau_calibration_identity.py`; existence below the price: `vaidya_tau_existence_no_deadline.py` |
| Remark on an example with three phases | crossing near d = 0.6, frozen minimiser for d ≥ 14, deep family above the price | `vaidya_step_continuation.py` (m₁ = 3/2, ε = 1/2), `vaidya_step_branches.py` |

## Notes on individual scripts

- `maxwell_rule_scan.py`: the equal-area rule holds by construction on the hull; the independent trapezoidal check of the
  areas is coarse (log-spaced grid of 160 points in the periastron) and is reported only as a sanity check; at
  q = 1.01 the bridge extends to periastron t = 1582.8 and the printed residual is −15.8.
- `tau_var_direct_optimisation.py` and `tau_var_late_drop.py` are local optimisations (SLSQP with a few starts); at
  η_d = 15 SLSQP does not converge and its iterates drift towards r = 2M, as stated in the paper.
- `tau_var_value_function.py`: the value of the discrete programme is a numerical estimate, neither a proved upper nor
  a proved lower bound of the continuous value (bilinear interpolation and Euler steps have errors of uncontrolled sign).
- `tau_var_numbers.py` reproduces the shooting/direct comparison with the same raw factor in both (normalising A to
  A(0) = 1 in one of them and not in the other changes N and gives a spurious difference of 3.5·10⁻⁴), the crossing
  with the radial-descent dive, the fold and the domains of the semi-Lagrangian programme.
- `tau_eq_injectivity_scan.py` writes `../results/tau_eq_injectivity_data.npy` (a pickled NumPy object array).
- The `vaidya_tau_*` scripts import one another; `vaidya_tau_existence_no_deadline.py` and
  `vaidya_tau_verification_function.py` can be imported without running their checks.
