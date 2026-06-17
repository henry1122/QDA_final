#!/usr/bin/env bash
# run_benchmarks.sh — build and run the PhasePoly benchmark on relevant circuits
#
# Usage:
#   ./scripts/run_benchmarks.sh                   # default: suitable circuits (no --set flag)
#   ./scripts/run_benchmarks.sh --set table1      # 19 circuits from PhasePoly paper Table 1
#   ./scripts/run_benchmarks.sh --set extended    # Table 1 + additional phase-poly-compatible circuits
#   ./scripts/run_benchmarks.sh --set all         # all 48 circuits (includes large/irrelevant ones)
#   ./scripts/run_benchmarks.sh --max-rz 20       # stricter per-block A* size limit
#   ./scripts/run_benchmarks.sh --no-todd         # skip Todd full-circuit comparison (faster)
#
# Circuit selection rationale:
#   "suitable" = "extended" (default):
#     Circuits with actual phase-polynomial structure (non-zero CNOT+Rz blocks) and all
#     blocks small enough for full PhasePoly A* within the default --max-rz=30 limit.
#     Excludes: pure-Clifford circuits (QFT8/16/32, QFTAdd*), circuits with very large
#     blocks that force A* fallback (gf2^*), and redundant scale-variants (Adder32/64).
#
# Results are printed to stdout and saved to results/benchmark_<timestamp>.txt

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BINARY="$REPO_ROOT/build/benchmark-phasepoly"
BENCH_DIR="$REPO_ROOT/benchmark/qc/optimized"
RESULTS_DIR="$REPO_ROOT/results"

# --- Circuit sets ---

# Table 1 from Chen et al. 2025 (PhasePoly paper)
TABLE1_CIRCUITS=(
    tof_3_pyzx tof_4_pyzx tof_5_pyzx tof_10_pyzx
    barenco_tof_3_pyzx barenco_tof_4_pyzx barenco_tof_5_pyzx barenco_tof_10_pyzx
    grover_5_pyzx
    ham15-low_pyzx ham15-med_pyzx ham15-high_pyzx
    mod5_4_pyzx mod_mult_55_pyzx mod_red_21_pyzx
    qcla_com_7_pyzx qcla_mod_7_pyzx rc_adder_6_pyzx
    vbe_adder_3_pyzx
)

# Additional circuits beyond Table 1 that are fully phase-poly-compatible:
# - have non-zero CNOT+Rz phase blocks
# - all blocks fit within the default max_rz=30 limit (no partial A*)
# - span a variety of circuit families and sizes
EXTRA_CIRCUITS=(
    Adder8_pyzx           # 23-qubit ripple-carry adder
    adder_8_pyzx          # alternate 24-qubit 8-bit adder
    csla_mux_3_original_pyzx  # carry-select look-ahead MUX
    csum_mux_9_corrected_pyzx # conditional sum MUX
    hwb6_pyzx             # hidden weighted bits, 6-bit
    mod_adder_1024_pyzx   # modular adder
    nth_prime6_pyzx       # nth prime oracle, 6-bit
    qcla_adder_10_pyzx    # carry-lookahead adder
    qft_4_pyzx            # QFT on 4 qubits (non-trivial phase structure)
)

# Circuits explicitly excluded and why:
#   QFT8/16/32, QFTAdd8/16/32 — 0 Rz blocks (pure Clifford, not phase-poly circuits)
#   gf2^* — large dense blocks (30+ Rz/block), PhasePoly A* must fall back for most blocks
#   Adder16/32/64 — redundant scaling of Adder8; included if you want a scaling study
#   hwb8, nth_prime8 — 1000+ blocks, very long runtime; use for scale experiments separately

# --- Parse flags ---
CIRCUIT_SET="extended"
BINARY_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --set)
            CIRCUIT_SET="$2"; shift 2 ;;
        --set=*)
            CIRCUIT_SET="${1#--set=}"; shift ;;
        *)
            BINARY_ARGS+=("$1"); shift ;;
    esac
done

# --- Build ---
echo "=== Building benchmark-phasepoly ==="
cmake --build "$REPO_ROOT/build" --target benchmark-phasepoly \
    -j"$(nproc 2>/dev/null || echo 4)" 2>&1
echo ""

# --- Collect circuit files ---
CIRCUIT_FILES=()

case "$CIRCUIT_SET" in
    table1)
        for name in "${TABLE1_CIRCUITS[@]}"; do
            f="$BENCH_DIR/${name}.qc"
            [[ -f "$f" ]] && CIRCUIT_FILES+=("$f") || echo "Warning: not found: $f" >&2
        done
        echo "=== Table 1 circuits (${#CIRCUIT_FILES[@]} files) ==="
        ;;
    extended|suitable|default)
        for name in "${TABLE1_CIRCUITS[@]}" "${EXTRA_CIRCUITS[@]}"; do
            f="$BENCH_DIR/${name}.qc"
            [[ -f "$f" ]] && CIRCUIT_FILES+=("$f") || echo "Warning: not found: $f" >&2
        done
        echo "=== Extended suitable set (${#CIRCUIT_FILES[@]} files: Table 1 + additional) ==="
        ;;
    all)
        mapfile -t CIRCUIT_FILES < <(find "$BENCH_DIR" -name "*.qc" | sort)
        echo "=== All circuits (${#CIRCUIT_FILES[@]} files, includes large/irrelevant ones) ==="
        ;;
    *)
        echo "Unknown --set value: $CIRCUIT_SET (expected: table1, extended, all)" >&2
        exit 1
        ;;
esac

if [[ ${#CIRCUIT_FILES[@]} -eq 0 ]]; then
    echo "No circuit files found." >&2; exit 1
fi

# --- Run ---
mkdir -p "$RESULTS_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTFILE="$RESULTS_DIR/benchmark_${TIMESTAMP}.txt"

{
    echo "PhasePoly benchmark — $(date)"
    echo "Binary   : $BINARY"
    echo "Set      : $CIRCUIT_SET (${#CIRCUIT_FILES[@]} circuits)"
    echo "Args     : ${BINARY_ARGS[*]:-<defaults>}"
    echo ""
    "$BINARY" "${BINARY_ARGS[@]}" "${CIRCUIT_FILES[@]}"
} 2>/dev/null | tee "$OUTFILE"

echo ""
echo "Results saved to: $OUTFILE"
