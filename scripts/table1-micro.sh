#!/usr/bin/env bash
# Micro-benchmark: 6 small circuits, no Todd (faster path).
set -euo pipefail
QSYN="${QSYN_BIN:-/tmp/qsyn-build/qsyn}"
cd "${QSYN_ROOT:-/app/qsyn}"

CIRCUITS=(tof_3_pyzx.qc tof_4_pyzx.qc grover_5_pyzx.qc vbe_adder_3_pyzx.qc hwb6_pyzx.qc mod5_4_pyzx.qc)

run() {
  local c="$1" s="$2"
  timeout 90 "$QSYN" -q 2>/dev/null <<EOF || echo "TIMEOUT"
logger error
qcir read benchmark/qc/optimized/${c}
convert qc tableau
convert tableau qcir -r ${s}
qcir print --stat
EOF
}

parse_gates() {
  awk '/^Clifford/{c=$3} /^T-family/{t=$3} END{print c+t+0}'
}
parse_cx() { awk '/^2-qubit/{print $3}'; }

echo "circuit,mst_g,mst_cx,gray_g,gray_cx,pp_g,pp_cx"
tm=0; tg=0; tpm=0; tgm=0; tpg=0

for c in "${CIRCUITS[@]}"; do
  mo=$(run "$c" mst); go=$(run "$c" graysynth); po=$(run "$c" phasepoly)
  mg=$(echo "$mo"|parse_gates); mcx=$(echo "$mo"|parse_cx)
  gg=$(echo "$go"|parse_gates); gcx=$(echo "$go"|parse_cx)
  pg=$(echo "$po"|parse_gates); pcx=$(echo "$po"|parse_cx)
  echo "$c,$mg,$mcx,$gg,$gcx,$pg,$pcx"
  tm=$((tm+mg)); tg=$((tg+gg)); tpm=$((tpm+pg))
  tgm=$((tgm+mcx)); tpg=$((tpg+pcx))
done

echo ""
echo "TOTAL mst=$tm gray=$tg pp=$tpm"
echo "pp/mst=$((100*tpm/tm))% pp/gray=$((100*tpm/tg))%"
echo "CNOT pp/mst=$((100*tpg/tgm))%"
