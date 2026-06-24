#!/usr/bin/env bash
# Build benchmark-phasepoly and run pp+ vs TODD+PP on Table 1 (19 circuits).
set -euo pipefail

ROOT="${QSYN_ROOT:-/app/qsyn}"
BUILD="${QSYN_BUILD_DIR:-/tmp/qsyn-build}"
BENCH_DIR="$ROOT/benchmark/qc/optimized"
OUT_DIR="$ROOT/results"
BINARY="$BUILD/benchmark-phasepoly"

TABLE1=(
  tof_3_pyzx tof_4_pyzx tof_5_pyzx tof_10_pyzx
  barenco_tof_3_pyzx barenco_tof_4_pyzx barenco_tof_5_pyzx barenco_tof_10_pyzx
  grover_5_pyzx
  ham15-low_pyzx ham15-med_pyzx ham15-high_pyzx
  mod5_4_pyzx mod_mult_55_pyzx mod_red_21_pyzx
  qcla_com_7_pyzx qcla_mod_7_pyzx rc_adder_6_pyzx
  vbe_adder_3_pyzx
)

CIRCUITS=()
for name in "${TABLE1[@]}"; do
  f="$BENCH_DIR/${name}.qc"
  [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
  CIRCUITS+=("$f")
done

export CXXFLAGS="${CXXFLAGS:-} -Wno-error=uninitialized"
cmake -B "$BUILD" -S "$ROOT" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD" --target benchmark-phasepoly -j"$(nproc)"

mkdir -p "$OUT_DIR"
TS="$(date +%Y%m%d_%H%M%S)"
PP_OUT="$OUT_DIR/benchmark_${TS}_pp.txt"
TODD_OUT="$OUT_DIR/benchmark_${TS}_todd_pp.txt"

run_bench() {
  local label="$1"
  shift
  {
    echo "PhasePoly benchmark — $(date)"
    echo "Binary   : $BINARY"
    echo "Set      : table1 (19 circuits)"
    echo "Args     : $*"
    echo ""
    "$BINARY" "$@" "${CIRCUITS[@]}"
  } | tee "$label"
}

echo "=== pp+ (no block Todd) ==="
run_bench "$PP_OUT" --no-todd

echo ""
echo "=== TODD+PP (--block-todd) ==="
run_bench "$TODD_OUT" --no-todd --block-todd

echo ""
echo "pp+ results:     $PP_OUT"
echo "TODD+PP results: $TODD_OUT"
