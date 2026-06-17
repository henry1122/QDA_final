#!/usr/bin/env bash
# Fast benchmark on readable QASM circuits (OLSQ + spidernest).
set -euo pipefail
QSYN="${QSYN_BIN:-/tmp/qsyn-build/qsyn}"
cd "${QSYN_ROOT:-/app/qsyn}"

CIRCUITS=(
  benchmark/qasm/spidernest_4_0.qasm
  benchmark/SABRE/OLSQ/tof_4_after_heavy.qasm
  benchmark/SABRE/OLSQ/tof_5_after_heavy.qasm
  benchmark/SABRE/OLSQ/barenco_tof_4_after_heavy.qasm
  benchmark/SABRE/OLSQ/barenco_tof_5_after_heavy.qasm
  benchmark/SABRE/OLSQ/vbe_adder_3_after_heavy.qasm
  benchmark/SABRE/OLSQ/rc_adder_6_after_heavy.qasm
)

run() {
  local f="$1" s="$2"
  timeout 120 "$QSYN" -q 2>/dev/null <<EOF || echo "TIMEOUT"
logger error
qcir read ${f}
convert qc tableau
convert tableau qcir -r ${s}
qcir print --stat
EOF
}

parse_gates() { awk '/^Clifford/{c=$3} /^T-family/{t=$3} END{print c+t+0}'; }
parse_cx() { awk '/^2-qubit/{print $3}'; }

echo "circuit,mst_g,mst_cx,gray_g,gray_cx,pp_g,pp_cx"
tm=0; tg=0; tpm=0; tmcx=0; tgcx=0; tpcx=0; n=0

for f in "${CIRCUITS[@]}"; do
  [[ -f "$f" ]] || continue
  name=$(basename "$f" .qasm)
  mo=$(run "$f" mst); go=$(run "$f" graysynth); po=$(run "$f" phasepoly)
  mg=$(echo "$mo"|parse_gates); mcx=$(echo "$mo"|parse_cx)
  gg=$(echo "$go"|parse_gates); gcx=$(echo "$go"|parse_cx)
  pg=$(echo "$po"|parse_gates); pcx=$(echo "$po"|parse_cx)
  echo "$name,$mg,$mcx,$gg,$gcx,$pg,$pcx"
  if [[ "$mg" =~ ^[0-9]+$ ]]; then
    tm=$((tm+mg)); tg=$((tg+gg)); tpm=$((tpm+pg))
    tmcx=$((tmcx+mcx)); tgcx=$((tgcx+gcx)); tpcx=$((tpcx+pcx))
    n=$((n+1))
  fi
done

echo ""
echo "=== TOTAL ($n circuits) ==="
echo "MST  gates=$tm cnot=$tmcx"
echo "Gray gates=$tg cnot=$tgcx"
echo "PP   gates=$tpm cnot=$tpcx"
[[ $tm -gt 0 ]] && echo "pp/mst gates: $((100*tpm/tm))% (paper ~79%)"
[[ $tg -gt 0 ]] && echo "pp/gray gates: $((100*tpm/tg))% (paper ~70%)"
[[ $tmcx -gt 0 ]] && echo "CNOT reduction vs MST: $((100-100*tpcx/tmcx))% (paper ~28%)"
