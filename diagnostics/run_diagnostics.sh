#!/usr/bin/env bash
# Runs the floating-point diagnostics (not certificates) and writes each log to ../results/diagnostics_<name>.txt.
# Usage: ./run_diagnostics.sh              all jobs, in order (several hours on one core)
#        ./run_diagnostics.sh NAME ...     only the named jobs
#        ./run_diagnostics.sh --list       list the jobs with their approximate duration
#        ./run_diagnostics.sh --status     show the running job and its progress bar (from another terminal)
# Standard error goes to ../results/diagnostics_<name>.stderr, removed when empty. PYTHON overrides the interpreter.
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
OUT=../results
# name | script and arguments | duration measured on one core (Apple silicon laptop, 2026)
JOBS=(
"conjugate_points_static|conjugate_points_static.py|40 s"
"conjugate_threshold_scan|conjugate_threshold_scan.py|45 s"
"conjugate_threshold_near_three_halves|conjugate_threshold_near_three_halves.py|10 s"
"conjugate_envelope_check|conjugate_envelope_check.py|20 s"
"conjugate_symmetric_asymptotics|conjugate_symmetric_asymptotics.py|5 min"
"conjugate_case_q6_5|conjugate_case_q6_5.py|5 s"
"equator_ring_frozen|equator_ring_frozen.py|20 s"
"arrival_ray_count|arrival_ray_count.py|65 min"
"arrival_ray_count_below_three_halves|arrival_ray_count_below_three_halves.py|45 min"
"maxwell_rule_scan|maxwell_rule_scan.py|10 s"
"maxwell_bridge_scan|maxwell_bridge_scan.py|25 s"
"maxwell_transition_scan|maxwell_transition_scan.py|25 s"
"vaidya_arrival_cut|vaidya_arrival_cut.py|15 min"
"vaidya_ray_count|vaidya_ray_count.py|12 min"
"vaidya_ray_count_low|vaidya_ray_count.py low|45 s"
"tau_var_shooting|tau_var_shooting.py|1 min"
"tau_var_saturation|tau_var_saturation.py|15 s"
"tau_var_asymmetric_dive|tau_var_asymmetric_dive.py|1 s"
"tau_var_drop_time|tau_var_drop_time.py|10 s"
"tau_var_fold|tau_var_fold.py|10 s"
"tau_var_direct_optimisation|tau_var_direct_optimisation.py|47 min"
"tau_var_late_drop|tau_var_late_drop.py|25 min"
"tau_var_numbers|tau_var_numbers.py|3 min"
"tau_var_value_function_frozen_1|tau_var_value_function.py frozen 0.05,0.025,0.05,11|5 s"
"tau_var_value_function_frozen_2|tau_var_value_function.py frozen 0.025,0.0125,0.025,21|30 s"
"tau_var_value_function_frozen_3|tau_var_value_function.py frozen 0.0125,0.00625,0.0125,41|7 min"
"tau_var_value_function_drop4_1|tau_var_value_function.py 4 0.05,0.025,0.05,11 2.566213 2.39,12.9,88.0|10 s"
"tau_var_value_function_drop4_2|tau_var_value_function.py 4 0.025,0.0125,0.025,21 2.566213 2.39,12.9,88.0|90 s"
"tau_var_value_function_drop4_3|tau_var_value_function.py 4 0.0125,0.00625,0.0125,41 2.566213 2.39,12.9,88.0|22 min"
"tau_var_value_function_up4_1|tau_var_value_function.py up4 0.05,0.025,0.05,11 9.821112 2.35,11.46,93.0|10 s"
"tau_var_value_function_up4_2|tau_var_value_function.py up4 0.025,0.0125,0.025,21 9.821112 2.35,11.46,93.0|80 s"
"tau_var_value_function_up4_3|tau_var_value_function.py up4 0.0125,0.00625,0.0125,41 9.821112 2.35,11.46,93.0|20 min"
"tau_eq_free_arrival_shooting|tau_eq_free_arrival_shooting.py|39 min"
"tau_eq_arrival_shooting|tau_eq_arrival_shooting.py|32 min"
"tau_eq_injectivity_scan|tau_eq_injectivity_scan.py|83 min"
"tau_eq_uniqueness_checks|tau_eq_uniqueness_checks.py|16 min"
"tau_eq_sector_change_checks|tau_eq_sector_change_checks.py|20 s"
"vaidya_tau_escape_price_bounds|vaidya_tau_escape_price_bounds.py|1 s"
"vaidya_tau_existence_no_deadline|vaidya_tau_existence_no_deadline.py|30 s"
"vaidya_tau_verification_function|vaidya_tau_verification_function.py|1 min"
"vaidya_tau_calibration_identity|vaidya_tau_calibration_identity.py|20 s"
"vaidya_step_continuation|vaidya_step_continuation.py 1.5 0.5|30 s"
"vaidya_step_continuation_m1.25|vaidya_step_continuation.py 1.25 0.5|30 s"
"vaidya_step_branches|vaidya_step_branches.py 1.5 0.5|7 min"
)
STATUS=$OUT/.diagnostics_status; export PROGRESS_FILE=$OUT/.diagnostics_progress
if [ "${1:-}" = "--status" ]; then
  if [ -f "$STATUS" ]; then cat "$STATUS"; else echo "no run in progress"; exit 0; fi
  read -r t0 < "$OUT/.diagnostics_start" 2>/dev/null && echo "  job running for $(( ($(date +%s) - t0)/60 )) min"
  [ -f "$PROGRESS_FILE" ] && sed 's/^/  /' "$PROGRESS_FILE"
  exit 0
fi
if [ "${1:-}" = "--list" ]; then
  for j in "${JOBS[@]}"; do IFS='|' read -r n c t <<< "$j"; printf '%-42s %-8s %s\n' "$n" "$t" "$c"; done; exit 0
fi
selected=()
for j in "${JOBS[@]}"; do
  IFS='|' read -r n c t <<< "$j"
  if [ $# -gt 0 ]; then case " $* " in *" $n "*) ;; *) continue ;; esac; fi
  selected+=("$j")
done
fail=0; k=0
for j in "${selected[@]}"; do
  IFS='|' read -r n c t <<< "$j"; k=$((k + 1))
  start=$(date +%s); rm -f "$PROGRESS_FILE"; echo "$start" > "$OUT/.diagnostics_start"
  echo "job $k/${#selected[@]}: $n (about $t), started $(date '+%H:%M')" | tee "$STATUS"
  # shellcheck disable=SC2086
  $PY $c > "$OUT/diagnostics_$n.txt" 2> "$OUT/diagnostics_$n.stderr"; rc=$?
  [ -s "$OUT/diagnostics_$n.stderr" ] || rm -f "$OUT/diagnostics_$n.stderr"
  # warnings quote the absolute path of the script: keep only the path inside the repository
  [ -f "$OUT/diagnostics_$n.stderr" ] && sed -i.bak "s#$(cd .. && pwd)/##g" "$OUT/diagnostics_$n.stderr" && rm -f "$OUT/diagnostics_$n.stderr.bak"
  echo "$n: exit $rc, $(( $(date +%s) - start )) s"
  [ $rc -eq 0 ] || fail=1
done
rm -f "$STATUS" "$PROGRESS_FILE" "$OUT/.diagnostics_start"
exit $fail
