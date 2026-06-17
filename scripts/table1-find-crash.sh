#!/usr/bin/env bash
set -uo pipefail
QSYN=/tmp/qsyn-build/qsyn
cd /app/qsyn
for c in $(/tmp/qsyn-build/qsyn -q -c 'logger error; benchmark table1 --circuits x' 2>&1 | head -1); do true; done
CIRCUITS=(
  tof_3_pyzx.qc tof_4_pyzx.qc tof_5_pyzx.qc tof_10_pyzx.qc
  barenco_tof_3_pyzx.qc barenco_tof_4_pyzx.qc barenco_tof_5_pyzx.qc barenco_tof_10_pyzx.qc
  grover_5_pyzx.qc ham15-low_pyzx.qc ham15-med_pyzx.qc ham15-high_pyzx.qc
  adder_8_pyzx.qc Adder8_pyzx.qc Adder16_pyzx.qc Adder32_pyzx.qc Adder64_pyzx.qc
  vbe_adder_3_pyzx.qc rc_adder_6_pyzx.qc qcla_adder_10_pyzx.qc
  mod_red_21_pyzx.qc mod_mult_55_pyzx.qc mod_adder_1024_pyzx.qc mod5_4_pyzx.qc
  hwb6_pyzx.qc hwb8_pyzx.qc
)
for c in "${CIRCUITS[@]}"; do
  printf "%-25s " "$c"
  if timeout 300 "$QSYN" -q -c "logger error; benchmark table1 --circuits $c" >/tmp/r.txt 2>/tmp/e.txt; then
    pp=$(grep '^| ' /tmp/r.txt | grep -v circuit | grep -v TOTAL | grep -v Paper | awk -F'|' '{gsub(/ /,"",$6); print $6}')
    echo "OK pp=$pp"
  else
    ec=$?
    if [[ $ec -eq 124 ]]; then echo "TIMEOUT"
    elif grep -q Assertion /tmp/e.txt 2>/dev/null; then echo "ASSERT: $(head -1 /tmp/e.txt)"
    else echo "FAIL ec=$ec: $(head -1 /tmp/e.txt)"
    fi
  fi
done
