#!/usr/bin/env bash
# Table 1 benchmark: compare MST / GraySynth / PhasePoly on paper circuits.
set -euo pipefail

QSYN="${QSYN_BIN:-/tmp/qsyn-build/qsyn}"
ROOT="${QSYN_ROOT:-/app/qsyn}"
cd "$ROOT"

if [[ ! -x "$QSYN" ]]; then
  echo "qsyn binary not found: $QSYN" >&2
  exit 1
fi

# Paper Table 1 circuit set (pyzx-decomposed Clifford+T .qc files).
CIRCUITS=(
  tof_3_pyzx.qc
  tof_4_pyzx.qc
  tof_5_pyzx.qc
  tof_10_pyzx.qc
  barenco_tof_3_pyzx.qc
  barenco_tof_4_pyzx.qc
  barenco_tof_5_pyzx.qc
  barenco_tof_10_pyzx.qc
  grover_5_pyzx.qc
  ham15-low_pyzx.qc
  ham15-med_pyzx.qc
  ham15-high_pyzx.qc
  adder_8_pyzx.qc
  Adder8_pyzx.qc
  Adder16_pyzx.qc
  Adder32_pyzx.qc
  Adder64_pyzx.qc
  vbe_adder_3_pyzx.qc
  rc_adder_6_pyzx.qc
  qcla_adder_10_pyzx.qc
  mod_red_21_pyzx.qc
  mod_mult_55_pyzx.qc
  mod_adder_1024_pyzx.qc
  mod5_4_pyzx.qc
  hwb6_pyzx.qc
  hwb8_pyzx.qc
)

parse_stat() {
  local label="$1"
  local total=$(( ${1:-0} + 0 )) # unused, keep shellcheck quiet
  local clifford two_qubit t_family depth
  clifford=$(grep -E '^Clifford' | awk '{print $3}')
  two_qubit=$(grep -E '^2-qubit' | awk '{print $3}')
  t_family=$(grep -E '^T-family' | awk '{print $3}')
  depth=$(grep -E '^Depth' | awk '{print $3}')
  local gates=$(( clifford + t_family ))
  echo "${label},${gates},${two_qubit},${t_family},${depth}"
}

run_strategy() {
  local circuit="$1"
  local strategy="$2"
  local out
  out=$("$QSYN" -v <<EOF 2>/dev/null
qcir read benchmark/qc/optimized/${circuit}
convert qc tableau
tableau opt phasepoly todd
convert tableau qcir -r ${strategy}
qcir print --stat
EOF
)
  local clifford two_qubit t_family depth gates
  clifford=$(echo "$out" | grep -E '^Clifford' | awk '{print $3}')
  two_qubit=$(echo "$out" | grep -E '^2-qubit' | awk '{print $3}')
  t_family=$(echo "$out" | grep -E '^T-family' | awk '{print $3}')
  depth=$(echo "$out" | grep -E '^Depth' | awk '{print $3}')
  gates=$(( clifford + t_family ))
  if [[ -z "$clifford" || -z "$two_qubit" ]]; then
    echo "ERR,0,0,0,0"
    echo "$out" >&2
    return 1
  fi
  echo "${gates},${two_qubit},${t_family},${depth}"
}

echo "circuit,mst_gates,mst_cnot,mst_t,gray_gates,gray_cnot,gray_t,pp_gates,pp_cnot,pp_t"

tot_mst_g=0 tot_mst_c=0
tot_gray_g=0 tot_gray_c=0
tot_pp_g=0 tot_pp_c=0
ok=0 fail=0

for c in "${CIRCUITS[@]}"; do
  path="benchmark/qc/optimized/${c}"
  if [[ ! -f "$path" ]]; then
    echo "${c},SKIP" >&2
  fi
  printf "%s" "$c"
  for strat in mst graysynth phasepoly; do
    if stats=$(run_strategy "$c" "$strat"); then
      IFS=',' read -r g cx t d <<< "$stats"
      printf ",%s,%s,%s" "$g" "$cx" "$t"
      if [[ "$strat" == "mst" ]]; then
        tot_mst_g=$(( tot_mst_g + g ))
        tot_mst_c=$(( tot_mst_c + cx ))
      elif [[ "$strat" == "graysynth" ]]; then
        tot_gray_g=$(( tot_gray_g + g ))
        tot_gray_c=$(( tot_gray_c + cx ))
      else
        tot_pp_g=$(( tot_pp_g + g ))
        tot_pp_c=$(( tot_pp_c + cx ))
        ok=$(( ok + 1 ))
      fi
    else
      printf ",ERR,ERR,ERR"
      fail=$(( fail + 1 ))
    fi
  done
  echo ""
done

echo ""
echo "=== TOTALS (${ok} circuits) ==="
echo "MST+Todd   gates=${tot_mst_g} cnot=${tot_mst_c}"
echo "Gray+Todd  gates=${tot_gray_g} cnot=${tot_gray_c}"
echo "PhasePoly  gates=${tot_pp_g} cnot=${tot_pp_c}"
if [[ $tot_mst_g -gt 0 ]]; then
  pct=$(( 100 * tot_pp_g / tot_mst_g ))
  echo "pp/mst gates ratio: ${pct}%"
fi
if [[ $tot_gray_g -gt 0 ]]; then
  pct=$(( 100 * tot_pp_g / tot_gray_g ))
  echo "pp/gray gates ratio: ${pct}%"
fi
if [[ $tot_mst_c -gt 0 ]]; then
  red=$(( 100 - 100 * tot_pp_c / tot_mst_c ))
  echo "CNOT reduction vs MST: ${red}%"
fi
