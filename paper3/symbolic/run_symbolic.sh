#!/usr/bin/env bash
# Runs the Wolfram Language checks and writes each log to ../results/symbolic_<name>.txt.
# Usage: ./run_symbolic.sh              all scripts (about 10 minutes; bridge_group4_group5 alone about 25)
#        ./run_symbolic.sh NAME ...     only the named scripts (without .wls)
#        ./run_symbolic.sh --check      check the existing logs against the expected outcomes, without running
#                                       (and without touching any file)
#        ./run_symbolic.sh --check-log NAME FILE   check one log file as if it were the log of NAME
#        ./run_symbolic.sh --status     show the running script (from another terminal)
# Requires wolframscript (Wolfram Engine or Mathematica). WOLFRAMSCRIPT overrides the command.
# Exit status 0 only if every script is accepted (see expected_outcomes.tsv for the rules).
#
# The output goes through a pipe on purpose: with `wolframscript -file x.wls > log`, printed lines were occasionally
# lost in our runs (3 of 3 runs of one script lost one or two lines), while `... | cat > log` and `WolframKernel
# -script` were complete in 5 of 5 runs. The computations are not affected, only the capture of the output.
set -u
cd "$(dirname "$0")"
WS=${WOLFRAMSCRIPT:-wolframscript}
OUT=../results
STATUS=$OUT/.symbolic_status

lines() { [ -n "$1" ] && printf '%s\n' "$1"; }
# normal form of a line: no leading blanks, timings {0.12, False} -> {T, False}, Wolfram temporaries d$1377 -> d$K
normline() { sed -E -e 's/^[[:space:]]+//' -e 's/\{[0-9.]+(`[0-9.]*)?, (True|False)\}/{T, \2}/g' -e 's/\$[0-9]+/$K/g'; }
logfile() { if [ "$1" = bridge_group4_group5 ]; then echo "$OUT/bridge_group4_group5.txt"; else echo "$OUT/symbolic_$1.txt"; fi; }

# check_log NAME FILE: prints "accepted" or the reasons for rejection; returns 0 only if accepted
check_log() {
  local n=$1 log=$2 row npass last got glast exp_fail got_fail bad=""
  row=$(grep -v '^#' expected_outcomes.tsv | awk -F'\t' -v n="$n" '$1 == n')
  if [ -z "$row" ]; then echo "no expected outcome for $n in expected_outcomes.tsv"; return 1; fi
  npass=$(printf '%s' "$row" | cut -f2); last=$(printf '%s' "$row" | cut -f3)
  if [ ! -s "$log" ]; then echo "log missing or empty"; return 1; fi
  got=$(grep -E '^[[:space:]]*(\[PASS\]|PASS)' "$log" | grep -c -v -E '^[[:space:]]*PASS [0-9]+ +FAIL [0-9]+')
  [ "$got" -eq "$npass" ] || bad="$bad; $got PASS lines, expected $npass"
  exp_fail=$(grep -v '^#' expected_failures.txt | awk -F'\t' -v n="$n" '$1 == n {print $2}' | sort)
  got_fail=$(grep -E '^[[:space:]]*(\[FAIL\]|FAIL)' "$log" | sed -E 's/^[[:space:]]+//' | sort)
  if [ "$got_fail" != "$exp_fail" ]; then
    bad="$bad; FAIL lines differ from expected_failures.txt:"
    bad="$bad$(diff <(lines "$exp_fail") <(lines "$got_fail") | grep -E '^[<>]' | sed -e 's/^</ missing:/' -e 's/^>/ unexpected:/' | tr '\n' ';')"
  fi
  glast=$(grep -v -E '^[[:space:]]*$|^[A-Za-z$]+::[A-Za-z]+:' "$log" | tail -1 | normline)
  [ "$glast" = "$last" ] || bad="$bad; last line differs from the expected one (log truncated?)"
  if [ -n "$bad" ]; then echo "REJECTED${bad}"; return 1; fi
  echo "accepted ($got PASS lines${exp_fail:+, expected FAIL lines present})"
}

case "${1:-}" in
  --status) [ -f "$STATUS" ] && cat "$STATUS" || echo "no run in progress"; exit 0 ;;
  --check-log) check_log "$2" "$3"; exit $? ;;
esac
check_only=0
if [ "${1:-}" = "--check" ]; then check_only=1; shift; fi
if [ $check_only -eq 0 ] && ! command -v "$WS" > /dev/null 2>&1; then
  echo "wolframscript not found ($WS). Install Wolfram Engine or Mathematica, or set WOLFRAMSCRIPT=/path/to/wolframscript."
  echo "Nothing was run; the logs in results/ are unchanged."; exit 1
fi
fail=0; k=0
if [ $# -gt 0 ]; then names=("$@")
else
  # the manifest expected_outcomes.tsv is the list of scripts; every listed script must exist and every script
  # in this folder must be listed
  names=(); while IFS= read -r n; do names+=("$n"); done < <(grep -v '^#' expected_outcomes.tsv | cut -f1)
  for n in "${names[@]}"; do [ -f "$n.wls" ] || { echo "$n: listed in expected_outcomes.tsv but $n.wls is missing"; fail=1; }; done
  for f in *.wls; do awk -F'\t' -v n="${f%.wls}" '$1 == n {found = 1} END {exit !found}' expected_outcomes.tsv || { echo "${f%.wls}.wls is not listed in expected_outcomes.tsv"; fail=1; }; done
fi
for n in "${names[@]}"; do
  k=$((k + 1)); log=$(logfile "$n")
  if [ $check_only -eq 1 ]; then
    if [ ! -f "$n.wls" ]; then echo "$n: script missing"; fail=1; continue; fi
    res=$(check_log "$n" "$log") || fail=1; echo "$n: $res"; continue
  fi
  start=$(date +%s)
  echo "script $k/${#names[@]}: $n, started $(date '+%H:%M')" | tee "$STATUS"
  if [ "$n" = bridge_group4_group5 ]; then
    PROGRESS_FILE=$OUT/.symbolic_progress $WS -file "$n.wls" ../certificates < /dev/null 2>&1 | cat > "$log"
  else
    $WS -file "$n.wls" < /dev/null 2>&1 | cat > "$log"
  fi
  st=("${PIPESTATUS[@]}")   # read at once: any later command overwrites PIPESTATUS
  res=$(check_log "$n" "$log"); ok=$?
  if [ "${st[0]}" -ne 0 ] || [ "${st[1]}" -ne 0 ]; then ok=1; res="REJECTED (non-zero exit status); log: $res"; fi
  echo "$n: wolframscript exit ${st[0]}, capture exit ${st[1]}, $(( $(date +%s) - start )) s, $res"
  [ $ok -eq 0 ] || fail=1
done
[ $check_only -eq 1 ] || rm -f "$STATUS" "$OUT/.symbolic_progress"
[ $fail -eq 0 ] && echo "ALL SCRIPTS ACCEPTED" || echo "SOME SCRIPTS REJECTED"
exit $fail
