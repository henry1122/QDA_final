#!/usr/bin/env python3
"""
benchmark_phasepoly.py — Compare PhasePoly A* against qsyn baselines.

Usage:
    python scripts/benchmark_phasepoly.py [--qsyn PATH] [--bench DIR] [--out FILE]

For each .qc circuit in the benchmark directory the script runs qsyn's
QCir -> Tableau -> QCir pipeline with each rotation-synthesis strategy and
records the gate counts.  Results are printed as a TSV table and optionally
written to a CSV file.

Pipeline per (circuit, strategy):
    qcir read <file>
    convert qcir tableau
    convert tableau qcir --rotation <strategy>
    qcir print --statistics
    quit -f

Strategies compared: naive, graysynth, gstair, mst, phasepoly
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

STRATEGIES = ["naive", "graysynth", "gstair", "mst", "phasepoly"]

DOFILE_TEMPLATE = """\
qcir read {circuit}
convert qcir tableau
convert tableau qcir --rotation {strategy}
qcir print --stat
quit -f
"""


def make_dofile(circuit: Path, strategy: str) -> str:
    return DOFILE_TEMPLATE.format(circuit=circuit, strategy=strategy)


def parse_gate_counts(output: str) -> dict:
    """Parse the formatted output of `qcir print --stat`."""
    counts = {}
    for line in output.splitlines():
        # Match lines like: "2-qubit    : 15"  or  "T-family   : 28"
        m = re.match(r"^\s*([\w\-]+)\s*:\s*(\d+)", line)
        if m:
            key = m.group(1).lower().replace("-", "_")
            counts[key] = int(m.group(2))
    return counts


def run_qsyn(qsyn_bin: Path, dofile_content: str, timeout: int = 60) -> str:
    """Run qsyn with the given dofile content; return stdout."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".dof", delete=False) as f:
        f.write(dofile_content)
        fname = f.name
    try:
        result = subprocess.run(
            [str(qsyn_bin), "-f", fname],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    finally:
        Path(fname).unlink(missing_ok=True)


def benchmark(qsyn_bin: Path, bench_dir: Path, strategies: list[str]) -> list[dict]:
    circuits = sorted(bench_dir.glob("*.qc"))
    if not circuits:
        print(f"No .qc files found in {bench_dir}", file=sys.stderr)
        return []

    rows = []
    for circ in circuits:
        row = {"circuit": circ.stem}
        for strat in strategies:
            output = run_qsyn(qsyn_bin, make_dofile(circ, strat))
            counts = parse_gate_counts(output)
            row[f"{strat}_cx"]  = counts.get("2_qubit", -1)
            row[f"{strat}_t"]   = counts.get("t_family", -1)
        rows.append(row)
        print(f"  {circ.stem}", flush=True)
    return rows


def print_table(rows: list[dict], strategies: list[str]) -> None:
    if not rows:
        return

    # Header
    headers = ["circuit"] + [f"{s}_cx" for s in strategies] + [f"{s}_t" for s in strategies]
    col_w = {h: max(len(h), max(len(str(r.get(h, ""))) for r in rows)) for h in headers}

    def row_str(r: dict) -> str:
        return "  ".join(str(r.get(h, "")).rjust(col_w[h]) for h in headers)

    print("  ".join(h.rjust(col_w[h]) for h in headers))
    print("  ".join("-" * col_w[h] for h in headers))
    for r in rows:
        print(row_str(r))


def write_csv(rows: list[dict], strategies: list[str], path: Path) -> None:
    import csv
    headers = ["circuit"] + [f"{s}_cx" for s in strategies] + [f"{s}_t" for s in strategies]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        w.writerows({h: r.get(h, "") for h in headers} for r in rows)
    print(f"Results written to {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--qsyn",  default="./build/qsyn",          help="path to qsyn binary")
    parser.add_argument("--bench", default="./benchmark/qc/optimized", help="benchmark directory")
    parser.add_argument("--out",   default=None,                    help="write results to CSV")
    parser.add_argument("--strategies", nargs="+", default=STRATEGIES, choices=STRATEGIES,
                        help="strategies to compare")
    args = parser.parse_args()

    qsyn_bin  = Path(args.qsyn)
    bench_dir = Path(args.bench)

    if not qsyn_bin.is_file():
        sys.exit(f"qsyn binary not found: {qsyn_bin}")
    if not bench_dir.is_dir():
        sys.exit(f"Benchmark directory not found: {bench_dir}")

    print(f"Benchmarking {bench_dir} with strategies: {args.strategies}")
    rows = benchmark(qsyn_bin, bench_dir, args.strategies)

    print()
    print_table(rows, args.strategies)

    if args.out:
        write_csv(rows, args.strategies, Path(args.out))


if __name__ == "__main__":
    main()
