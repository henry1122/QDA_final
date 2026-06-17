#!/usr/bin/env bash
# Quick benchmark: pyzx .qc + OLSQ qasm after PhasePoly multiblock opt.
set -euo pipefail
QSYN="${QSYN_BIN:-/tmp/qsyn-build/qsyn}"
cd "${QSYN_ROOT:-/app/qsyn}"

bench() {
  local f="$1"
  local name
  name=$(basename "$f")
  name=${name%.*}
  local before after b_cx a_cx
  before=$($QSYN -q -c "logger error; qcir read $f; qcir print --stat" 2>/dev/null)
  after=$($QSYN -q -c "logger error; qcir read $f; qcir optimize phasepoly; qcir print --stat" 2>/dev/null)
  b_cx=$(echo "$before" | awk '/^2-qubit/{print $3}')
  a_cx=$(echo "$after" | awk '/^2-qubit/{print $3}')
  bg=$(echo "$before" | awk '/^Clifford/{c=$3} /^T-family/{t=$3} END{print c+t+0}')
  ag=$(echo "$after" | awk '/^Clifford/{c=$3} /^T-family/{t=$3} END{print c+t+0}')
  echo "$name,$bg,$b_cx,$ag,$a_cx"
}

echo "circuit,gates_before,cnot_before,gates_after,cnot_after"
bench benchmark/qc/optimized/tof_3_pyzx.qc
bench benchmark/qc/optimized/vbe_adder_3_pyzx.qc
bench benchmark/qasm/spidernest_4_0.qasm
bench benchmark/SABRE/OLSQ/barenco_tof_4_after_heavy.qasm
bench benchmark/SABRE/OLSQ/vbe_adder_3_after_heavy.qasm
