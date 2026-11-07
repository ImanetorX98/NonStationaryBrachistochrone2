#!/usr/bin/env bash
# Reruns the checks of the repository and verifies each outcome. Logs are rewritten in results/.
#
# Usage: ./run_all.sh            the standard set, about one hour (certificates, bridges, symbolic checks, Python
#                                checks, independent checks, Riemann sums, edges and critical separations)
#        ./run_all.sh --full     also the long runs: proof-tree search, unimodality (about 2 h on 10 processes),
#                                Maxwell angle (about 75 min), full Wolfram leaf rechecks (1-4 h each), figures and
#                                the floating-point diagnostics (about 7 h)
#        ./run_all.sh --list     list the jobs without running them
#        ./run_all.sh --check    apply the acceptance rules to the existing logs, without running anything
#        ./run_all.sh [--full] --from K   start from job K of the list (to resume an interrupted run)
# Environment: PYTHON (default python3, with the packages of requirements.txt); WOLFRAMSCRIPT (default wolframscript).
# Without wolframscript the Wolfram jobs are reported as SKIPPED.
# Exit status: 0 if every job ran and passed, 1 if a job failed, 2 if none failed but some were skipped.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
WS=${WOLFRAMSCRIPT:-wolframscript}
R=$ROOT/results
mode=standard; from=1; nmode=0
while [ $# -gt 0 ]; do
  case $1 in --full) mode=full; nmode=$((nmode + 1)) ;; --list) mode=list; nmode=$((nmode + 1)) ;;
    --check) mode=check; nmode=$((nmode + 1)) ;; --from) from=${2-}; shift ;;
    *) echo "unknown option $1"; exit 1 ;; esac
  shift
done
if [ $nmode -gt 1 ]; then echo "--full, --list and --check are alternatives: give at most one of them"; exit 1; fi
case $from in ''|*[!0-9]*) echo "--from needs a positive integer"; exit 1 ;; esac

# job table: set | name | folder | log (in results/) | rule | command
# A rule is a list of conditions joined by " && "; besides the conditions, the command and the capture of its log must
# exit with status 0, the log must not be empty and no line may start with FAIL or [FAIL].
#   exit          nothing more (jobs that write several logs of their own)
#   count:RE:N    exactly N lines match the extended regular expression RE
#   last:TEXT     the last non-empty line equals TEXT, after removing timings such as " (723s)" or "; 4432s"
#   names:RE      the 33 certificates listed in certificates/SHA256SUMS must be distinct and be exactly the JSON files
#                 of certificates/; for each of them exactly one line matches RE, with %s
#                 replaced by its name (without .json)
#   grid:RE       for each of the 20 grid cases q in {1.2, 1.5, 2, 4, 10}, r0 in {4, 6, 10, 20} exactly one line
#                 "q= Q r0= R:" is followed, right after the colon and one blank, by RE
# With --check the exit status is not available, so only the content of the logs is checked.
JOBS=$(cat <<'EOF'
standard|verify certificates|certificates|certificates_verify_all.txt|names:^VERIFIED \(local\) %s\.json: && last:OUTCOME: all certificates verified|$PY verify_certificate.py *.json
standard|tampering tests|certificates|certificates_tamper_tests.txt|count:REJECTED \(expected\):15 && last:SUMMARY: 15/15 tamperings rejected|$PY tamper_tests.py
standard|certificate bridge (Python)|certificate_generators|certificate_generators_bridge.txt|count:^ *\[?PASS:12 && last:SUMMARY: 12/12 PASS|$PY bridge.py
standard|rebuild certificates|certificate_generators|certificate_generators_build_check.txt|names:^MATCH +%s: && last:SUMMARY: 33/33 certificates rebuilt identically|$PY build_certificates.py
standard|interval edges|interval_arithmetic|interval_edges.txt|grid:lower edge .* OK .*upper edge .* OK|$PY edges.py
standard|critical separations|interval_arithmetic|interval_critical_separation.txt|grid:CERTIFIED delta_c/pi|$PY critical_separation.py
standard|Riemann sums (Vaidya example)|riemann_sums|riemann_sums_vaidya_example_constants.txt|count:^ *\[?PASS:24 && last:SUMMARY: 24/24 PASS|$PY vaidya_example_constants.py
standard|auxiliary metric null cone|python_checks|python_auxiliary_metric_null_cone.txt|count:^ *\[?PASS:13 && last:all checks passed|$PY auxiliary_metric_null_cone.py
standard|global hyperbolicity|python_checks|python_global_hyperbolicity.txt|count:^ *\[?PASS:11 && last:all checks passed|$PY global_hyperbolicity.py
standard|Kerr barrier|python_checks|python_kerr_barrier.txt|count:^ *\[?PASS:9 && last:all checks passed|$PY kerr_barrier.py
standard|Kerr 3D excursion|python_checks|python_kerr_3d_excursion.txt|count:^ *\[?PASS:9 && last:all checks passed|$PY kerr_3d_excursion.py
standard|independent: sharpness constants|independent_checks|independent_sharp32_constants.txt|count:^ *\[?PASS:35 && last:PASS 34; FAIL 0 (negative controls included)|$PY sharp32_constants.py
standard|independent: equator example|independent_checks|independent_equator_spin_example_constants.txt|count:^ *\[?PASS:56 && last:No real-model derivative budget established by this script.|$PY equator_spin_example_constants.py
standard|independent: Vaidya constants|independent_checks|independent_vaidya_example_constants_exact.txt|count:^ *\[?PASS:38 && last:SUMMARY: 38/38 PASS|$PY vaidya_example_constants_exact.py
standard|independent: unimodality archives|independent_checks|independent_unimodality_leaf_archives.txt|count:^PASS leaves:20 && last:SUMMARY: 20/20 archives PASS|$PY unimodality_leaf_archives.py
wolfram|symbolic checks (40 scripts)|symbolic|-|exit|./run_symbolic.sh
standard|signs of the group 4 factors|symbolic|bridge_group4_factor_signs.txt|count:^PASS group4_:12 && last:SUMMARY: all 12 factors F are positive on their boxes|$PY check_factor_signs.py
wolfram|Wolfram bridge, group 1|independent_checks|independent_wolfram_bridge_group1.txt|count:^PASS:21 && last:SUMMARY: 21 PASS, 0 FAIL|$WS -file wolfram_bridge_group1.wls ../certificates
wolfram|Wolfram bridge, group 2|independent_checks|independent_wolfram_bridge_group2.txt|count:^PASS:21 && last:SUMMARY: 21 PASS, 0 FAIL|$WS -file wolfram_bridge_group2.wls ../certificates
wolfram|Wolfram bridge, group 3|independent_checks|independent_wolfram_bridge_group3.txt|count:^PASS:15 && last:SUMMARY: 15 PASS, 0 FAIL|$WS -file wolfram_bridge_group3.wls ../certificates
full|proof-tree search|certificate_generators|certificate_generators_build_search.txt|names:^MATCH +%s: .*tree identical && last:SUMMARY: 33/33 certificates rebuilt identically|$PY build_certificates.py --search
full|unimodality (20 cases)|interval_arithmetic|interval_unimodality.txt|grid:CERTIFIED unimodal|$PY unimodality.py
full|independent: unimodality archives (after the rerun)|independent_checks|independent_unimodality_leaf_archives.txt|count:^PASS leaves:20 && last:SUMMARY: 20/20 archives PASS|$PY unimodality_leaf_archives.py
full|Maxwell angle|interval_arithmetic|interval_maxwell_angle.txt|count:^C: middle cover OK:1 && last:RESULT: CERTIFIED interior Maxwell angle delta* in [0.677106398319, 0.677112764518] pi|$PY maxwell_angle.py
fullwolfram|Wolfram leaf recheck, group3_I2_main|independent_checks|independent_wolfram_leaf_recheck_group3_I2_main.txt|count:^PASS cell:2231 && last:RESULT: PASS (whole tree read and unit volume: True; all checked leaves pass: True)|$WS -file wolfram_leaf_recheck.wls ../certificates/group3_I2_main.json all
fullwolfram|Wolfram leaf recheck, group4_B1_I1|independent_checks|independent_wolfram_leaf_recheck_group4_B1_I1.txt|count:^PASS cell:50 && last:RESULT: PASS (whole tree read and unit volume: True; all checked leaves pass: True)|$WS -file wolfram_leaf_recheck.wls ../certificates/group4_B1_I1.json all
full|figure 2|figures|-|exit|$PY make_fig_radius.py
full|figure 7|figures|-|exit|$PY make_fig_tau_ring.py
full|figure of the supplement|figures|-|exit|$PY make_fig_bands_supplement.py
full|floating-point diagnostics|diagnostics|-|exit|PYTHON="$PY" ./run_diagnostics.sh
EOF
)

normlast() { sed -E -e 's/ \([0-9.]+ ?s\)$//' -e 's/[;,] [0-9.]+ ?s$//'; }
check_cond() {  # condition log -> 0 if satisfied
  local c=$1 log=$2 re n got last name q r0
  case $c in
    exit) return 0 ;;
    count:*) re=${c#count:}; n=${re##*:}; re=${re%:*}
             got=$(grep -c -E "$re" "$log"); [ "$got" -eq "$n" ] ;;
    last:*) last=$(grep -v '^[[:space:]]*$' "$log" | tail -1 | normlast); [ "$last" = "${c#last:}" ] ;;
    names:*) re=${c#names:}
             # 33 distinct certificate names, the same set as the JSON files present in certificates/
             listed=$(awk '$2 ~ /\.json$/ {print $2}' "$ROOT/certificates/SHA256SUMS" | sort)
             [ "$(printf '%s\n' "$listed" | sort -u | grep -c .)" -eq 33 ] || return 1
             [ "$listed" = "$(cd "$ROOT/certificates" && ls *.json | sort)" ] || return 1
             for name in $(printf '%s\n' "$listed" | sed 's/\.json$//'); do
               got=$(grep -c -E "$(printf "$re" "$name")" "$log"); [ "$got" -eq 1 ] || return 1
             done ;;
    grid:*) re=${c#grid:}
            for q in 1.2 1.5 2 4 10; do for r0 in 4 6 10 20; do
              got=$(grep -c -E "^q= *${q//./\\.} r0= *$r0\.0: $re" "$log"); [ "$got" -eq 1 ] || return 1
            done; done ;;
    *) return 1 ;;
  esac
}
check_rule() {  # rule log status -> 0 if accepted
  local rule=$1 log=$2 st=$3 rest c
  [ "$st" -eq 0 ] || return 1
  if [ "$log" != /dev/null ]; then
    [ -s "$log" ] || return 1
    grep -q -E '^[[:space:]]*\[?FAIL' "$log" && return 1
  fi
  rest=$rule
  while [ -n "$rest" ]; do
    c=${rest%% && *}; check_cond "$c" "$log" || return 1
    [ "$c" = "$rest" ] && break; rest=${rest#* && }
  done
  return 0
}

have_ws=0; command -v "$WS" > /dev/null 2>&1 && have_ws=1

# Before running anything: the Python must exist and import every package of requirements.txt.
if [ "$mode" = standard ] || [ "$mode" = full ]; then
  if ! command -v "$PY" > /dev/null 2>&1; then
    echo "Python not found: $PY. Set PYTHON=/path/to/python3 (see README.md, Setting up)."; exit 1
  fi
  missing=$("$PY" - "$ROOT/requirements.txt" 2>&1 <<'PY_CHECK'
import importlib.util, sys
names = []
for line in open(sys.argv[1]):
    line = line.split('#')[0].strip()
    if line:
        names.append(line.split('==')[0].split('>=')[0].strip())
print(' '.join(n for n in names if importlib.util.find_spec(n.lower()) is None))
PY_CHECK
) || { echo "Python at $(command -v "$PY") could not check the packages: $missing"; exit 1; }
  echo "Python:        $(command -v "$PY") ($("$PY" -c 'import sys; print(sys.version.split()[0])'))"
  if [ -n "$missing" ]; then
    echo "This Python lacks: $missing"
    echo "Install them with: $PY -m pip install -r requirements.txt"
    echo "or choose another Python with PYTHON=/path/to/python3 (see README.md, Setting up). Nothing was run."
    exit 1
  fi
  if [ $have_ws = 1 ]; then echo "wolframscript: $(command -v "$WS")"
  else echo "wolframscript: not found ($WS); the Wolfram jobs will be skipped (set WOLFRAMSCRIPT to use another one)"; fi
fi
selected=$(printf '%s\n' "$JOBS" | awk -F'|' -v m="$mode" '$1 == "standard" || $1 == "wolfram" || (m != "standard" && ($1 == "full" || $1 == "fullwolfram"))')
total=$(printf '%s\n' "$selected" | grep -c .)
if [ "$from" -lt 1 ] || [ "$from" -gt "$total" ]; then echo "--from must be between 1 and $total"; exit 1; fi
npass=0; nfail=0; nskip=0; k=0; failed=""; t_all=$SECONDS
while IFS='|' read -r set name dir log rule cmd; do
  k=$((k + 1))
  if [ "$mode" = list ]; then
    case $set in full|fullwolfram) tag="--full" ;; *) tag="" ;; esac
    printf '%2d. %-6s %-52s %s/  %s\n' "$k" "$tag" "$name" "$dir" "$cmd"; continue
  fi
  if [ "$k" -lt "$from" ]; then continue; fi
  if [ "$mode" = check ]; then
    if [ "$dir" = symbolic ] && [ "$log" = - ]; then
      if (cd "$ROOT/symbolic" && ./run_symbolic.sh --check > /dev/null); then res=PASS; npass=$((npass + 1)); else res="FAIL (see symbolic/run_symbolic.sh --check)"; nfail=$((nfail + 1)); failed="$failed\n  $name"; fi
      printf '%2d. %-55s %s\n' "$k" "$name" "$res"; continue
    fi
    if [ "$log" = - ]; then printf '%2d. %-55s (no single log; see its folder)\n' "$k" "$name"; nskip=$((nskip + 1)); continue; fi
    if check_rule "$rule" "$R/$log" 0; then res=PASS; npass=$((npass + 1)); else res="FAIL ($rule)"; nfail=$((nfail + 1)); failed="$failed\n  $name"; fi
    printf '%2d. %-55s %s\n' "$k" "$name" "$res"; continue
  fi
  if { [ "$set" = wolfram ] || [ "$set" = fullwolfram ]; } && [ $have_ws -eq 0 ]; then
    printf '[%d/%d] %-55s SKIPPED (no %s)\n' "$k" "$total" "$name" "$WS"; nskip=$((nskip + 1)); continue
  fi
  printf '[%d/%d] %-55s started %s\n' "$k" "$total" "$name" "$(date '+%H:%M')"
  t0=$SECONDS
  if [ "$log" = - ]; then
    (cd "$ROOT/$dir" && eval "$cmd" < /dev/null); st=$?; logpath=/dev/null
  else
    logpath=$R/$log
    (cd "$ROOT/$dir" && eval "$cmd" < /dev/null | cat > "$logpath"; s=("${PIPESTATUS[@]}"); [ "${s[0]}" -eq 0 ] && [ "${s[1]}" -eq 0 ]); st=$?
  fi
  if check_rule "$rule" "$logpath" "$st"; then
    res=PASS; npass=$((npass + 1))
  else
    res="FAIL (exit status $st; rule $rule)"; nfail=$((nfail + 1)); failed="$failed\n  $name"
  fi
  printf '       %-55s %s, %d s\n' "" "$res" $((SECONDS - t0))
done <<< "$selected"
[ "$mode" = list ] && exit 0
# every selected job must have been accounted for (a command that read standard input once swallowed the job list)
if [ $((npass + nfail + nskip)) -ne $((total - from + 1)) ]; then
  echo "SUMMARY: only $((npass + nfail + nskip)) of $((total - from + 1)) jobs were run"; exit 1
fi
echo "SUMMARY: $npass passed, $nfail failed, $nskip skipped, $(( (SECONDS - t_all) / 60 )) min"
[ $nfail -gt 0 ] && { printf "failed:$failed\n"; exit 1; }
[ $nskip -gt 0 ] && exit 2
exit 0
