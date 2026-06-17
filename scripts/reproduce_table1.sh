#!/usr/bin/env bash
# Full Table 1 reproduction: one qsyn process per circuit (crash-safe, timeout per circuit).
set -euo pipefail
ROOT="${QSYN_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
QSYN="${QSYN_BIN:-/tmp/qsyn-build/qsyn}"
OUT="${1:-$ROOT/doc/table1_results.md}"
TIMEOUT="${TABLE1_TIMEOUT:-600}"
cd "$ROOT"
mkdir -p "$(dirname "$OUT")"

CIRCUITS=(
  tof_3_pyzx.qc tof_4_pyzx.qc tof_5_pyzx.qc tof_10_pyzx.qc
  barenco_tof_3_pyzx.qc barenco_tof_4_pyzx.qc barenco_tof_5_pyzx.qc barenco_tof_10_pyzx.qc
  grover_5_pyzx.qc ham15-low_pyzx.qc ham15-med_pyzx.qc ham15-high_pyzx.qc
  adder_8_pyzx.qc Adder8_pyzx.qc Adder16_pyzx.qc Adder32_pyzx.qc Adder64_pyzx.qc
  vbe_adder_3_pyzx.qc rc_adder_6_pyzx.qc qcla_adder_10_pyzx.qc
  mod_red_21_pyzx.qc mod_mult_55_pyzx.qc mod_adder_1024_pyzx.qc mod5_4_pyzx.qc
  hwb6_pyzx.qc hwb8_pyzx.qc
)

CSV="$ROOT/doc/table1_results.csv"
PHASE_CSV="$ROOT/doc/table1_results_phase.csv"
echo "circuit,q,blk,Rz,pp,mst,gstair,gray,naive,todd_naive,pp_mst_pct,pp_gray_pct,ok,error" > "$CSV"
echo "# phase-only: CX+Rz inside synthesized phase-polynomial blocks (excludes H/Toffoli/T boundaries)" > "$PHASE_CSV"
echo "circuit,q,blk,Rz,pp,mst,gstair,gray,naive,todd_naive,pp_mst_pct,pp_gray_pct,ok,error" >> "$PHASE_CSV"

tot_pp=0 tot_mst=0 tot_gray=0 tot_gstair=0 tot_naive=0 n_ok=0
tot_pp_p=0 tot_mst_p=0 tot_gray_p=0

for c in "${CIRCUITS[@]}"; do
  echo "[$c] ..."
  if timeout "$TIMEOUT" "$QSYN" -q -c "logger error; benchmark table1 --format csv --circuits $c --output /tmp/row_one.csv --phase-output /tmp/row_phase_one.csv" > /tmp/row.err 2>&1; then
    tail -n +2 /tmp/row_one.csv >> "$CSV"
    tail -n 1 /tmp/row_phase_one.csv >> "$PHASE_CSV"
    IFS=',' read -r _ q blk rz pp mst gst gray naive _ pm pg ok err < <(tail -1 /tmp/row_one.csv)
    if [[ "$ok" == "1" ]]; then
      tot_pp=$((tot_pp+pp)); tot_mst=$((tot_mst+mst)); tot_gray=$((tot_gray+gray))
      tot_gstair=$((tot_gstair+gst)); tot_naive=$((tot_naive+naive)); n_ok=$((n_ok+1))
      IFS=',' read -r _ _ _ _ pp_p mst_p _ gray_p _ _ _ _ < <(tail -1 /tmp/row_phase_one.csv)
      tot_pp_p=$((tot_pp_p+pp_p)); tot_mst_p=$((tot_mst_p+mst_p)); tot_gray_p=$((tot_gray_p+gray_p))
    fi
  else
    ec=$?
    err="timeout or crash"
    [[ $ec -eq 124 ]] && err="timeout ${TIMEOUT}s"
    grep -q Assertion /tmp/row.err 2>/dev/null && err="assertion failed"
    echo "${c%.qc},,,,,,,,,,0,$err" >> "$CSV"
    echo "${c%.qc},,,,,,,,,,0,$err" >> "$PHASE_CSV"
  fi
done

python3 - "$CSV" "$OUT" "$tot_pp" "$tot_mst" "$tot_gstair" "$tot_gray" "$tot_naive" "$n_ok" <<'PY'
import csv, sys
csv_path, md_path, tp, tm, tg, ty, tn, n_ok = sys.argv[1:9]
tp,tm,tg,ty,tn,n_ok = map(int, [tp,tm,tg,ty,tn,n_ok])
rows = list(csv.DictReader(open(csv_path)))
with open(md_path, "w") as f:
    f.write("| circuit | q | blk | Rz | pp | mst+P | gstair+P | gray+P | naive+P | pp/mst | pp/gray |\n")
    f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in rows:
        if r.get("ok") != "1":
            f.write(f"| {r['circuit']} | ERR | | | | | | | | | |\n")
            continue
        pp,mst,gray = int(r['pp']), int(r['mst']), int(r['gray'])
        pm = f"{100*pp/mst:.1f}%" if mst else "-"
        pg = f"{100*pp/gray:.1f}%" if gray else "-"
        f.write(f"| {r['circuit']} | {r['q']} | {r['blk']} | {r['Rz']} | {pp} | {mst} | {r['gstair']} | {gray} | {r['naive']} | {pm} | {pg} |\n")
    f.write(f"\n| **TOTAL** ({n_ok} ok) | | | | **{tp}** | **{tm}** | **{tg}** | **{ty}** | **{tn}** | ")
    f.write(f"{100*tp//tm}% | {100*tp//ty}% |\n" if tm and ty else "| - | - |\n")
    f.write("\nPaper reference TOTAL: pp=3869, mst+P=4891 (79%), gray+P=5511 (70%)\n")
print(f"Wrote {md_path}")
PY
