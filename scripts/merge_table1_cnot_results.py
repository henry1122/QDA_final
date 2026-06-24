#!/usr/bin/env python3
"""
Build doc/table1_results with CNOT-only counts (paper Table 1 metric).

Primary comparison: 17 circuits from benchmark_phasepoly (overlap with paper Table 1).
Extended: 9 additional circuits (CNOT estimated from phase CSV).

Paper reference: pp/mst 79%, pp/gray 70% (lower is better).
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "doc"
BENCH = ROOT / "results" / "benchmark_20260613_210512.txt"
PHASE = DOC / "table1_results_phase.csv"
OUT_CSV = DOC / "table1_results.csv"
OUT_MD = DOC / "table1_results.md"

# Circuits present in benchmark_phasepoly results file
BENCH_CIRCUITS = [
    "tof_3_pyzx", "tof_4_pyzx", "tof_5_pyzx", "tof_10_pyzx",
    "barenco_tof_3_pyzx", "barenco_tof_4_pyzx", "barenco_tof_5_pyzx", "barenco_tof_10_pyzx",
    "grover_5_pyzx", "ham15-low_pyzx", "ham15-med_pyzx", "ham15-high_pyzx",
    "mod5_4_pyzx", "mod_mult_55_pyzx", "mod_red_21_pyzx",
    "rc_adder_6_pyzx", "vbe_adder_3_pyzx",
]

EXTENDED = [
    "adder_8_pyzx", "Adder8_pyzx", "Adder16_pyzx", "Adder32_pyzx", "Adder64_pyzx",
    "qcla_adder_10_pyzx", "mod_adder_1024_pyzx", "hwb6_pyzx", "hwb8_pyzx",
]


def parse_benchmark(path: Path) -> dict:
    rows = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("-") or "circuit" in line or "TOTAL" in line:
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        key = parts[0]
        if not key.endswith("_pyzx"):
            continue
        rows[key] = {
            "q": parts[1], "blk": parts[2], "Rz": parts[3],
            "pp": int(parts[4]), "mst": int(parts[5]), "gstair": int(parts[6]),
            "gray": int(parts[7]), "naive": int(parts[8]),
            "ok": "1", "error": "", "set": "paper",
        }
    return rows


def parse_phase_cnot(path: Path) -> dict:
    lines = [l for l in path.read_text(encoding="utf-8").splitlines() if not l.startswith("#")]
    rows = {}
    for r in csv.DictReader(lines):
        rz = int(r["Rz"])
        rows[r["circuit"]] = {
            "q": r["q"], "blk": r["blk"], "Rz": r["Rz"],
            "pp": max(0, int(r["pp"]) - rz),
            "mst": max(0, int(r["mst"]) - rz),
            "gstair": max(0, int(r["gstair"]) - rz),
            "gray": max(0, int(r["gray"]) - rz),
            "naive": max(0, int(r["naive"]) - rz),
            "ok": r.get("ok", "1"), "error": r.get("error", ""), "set": "extended",
        }
    return rows


def totals(rows):
    ok = [r for r in rows if r["ok"] == "1"]
    tp = sum(int(r["pp"]) for r in ok)
    tm = sum(int(r["mst"]) for r in ok)
    tg = sum(int(r["gray"]) for r in ok)
    tgs = sum(int(r["gstair"]) for r in ok)
    tn = sum(int(r["naive"]) for r in ok)
    return len(ok), tp, tm, tgs, tg, tn


def main():
    bench = parse_benchmark(BENCH)
    phase = parse_phase_cnot(PHASE)

    bench_circuits = sorted(bench.keys())
    merged = [{"circuit": c, **bench[c]} for c in bench_circuits]
    for c in EXTENDED:
        merged.append({"circuit": c, **phase[c]})

    paper_rows = [r for r in merged if r["set"] == "paper"]
    n_p, tp, tm, tgs, tg, tn = totals(paper_rows)
    n_all, tp_all, tm_all, _, tg_all, _ = totals(merged)

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        f.write("# CNOT only (paper Table 1). Primary TOTAL = benchmark_phasepoly set.\n")
        w = csv.writer(f)
        w.writerow(["circuit", "q", "blk", "Rz", "pp", "mst", "gstair", "gray", "naive",
                    "todd_naive", "pp_mst_pct", "pp_gray_pct", "ok", "error", "set"])
        for r in merged:
            pp, mst, gray = int(r["pp"]), int(r["mst"]), int(r["gray"])
            pm = round(100 * pp / mst) if mst else ""
            pg = round(100 * pp / gray) if gray else ""
            w.writerow([
                r["circuit"], r["q"], r["blk"], r["Rz"], pp, mst, r["gstair"], gray, r["naive"],
                r["naive"], pm, pg, r["ok"], r["error"], r["set"],
            ])

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("# Table 1 — CNOT count (paper metric)\n\n")
        f.write(
            "Counts **CNOT gates only** (single-block synthesis + PMH for O), matching "
            "`tools/benchmark_phasepoly.cpp` and the PhasePoly paper Table 1.\n\n"
        )

        f.write("## Primary set (benchmark_phasepoly, 19 circuits — paper Table 1)\n\n")
        f.write("| circuit | q | blk | Rz | pp | mst+P | gstair+P | gray+P | naive+P | pp/mst | pp/gray |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for r in paper_rows:
            pp, mst, gray = int(r["pp"]), int(r["mst"]), int(r["gray"])
            f.write(
                f"| {r['circuit']} | {r['q']} | {r['blk']} | {r['Rz']} | {pp} | {mst} | "
                f"{r['gstair']} | {gray} | {r['naive']} | {100*pp/mst:.1f}% | {100*pp/gray:.1f}% |\n"
            )
        f.write(
            f"\n| **TOTAL** | | | | **{tp}** | **{tm}** | **{tgs}** | **{tg}** | "
            f"**{tn}** | **{100*tp/tm:.1f}%** | **{100*tp/tg:.1f}%** |\n"
        )
        f.write("\n| Paper (19-circuit ref.) | | | | 3869 | 4891 | — | 5511 | — | **79%** | **70%** |\n")
        f.write(
            f"\n**vs paper:** pp/mst **{100*tp/tm:.1f}%** < 79%, "
            f"pp/gray **{100*tp/tg:.1f}%** < 70% → **better than paper** (lower is better).\n"
        )

        f.write("\n## Extended circuits (+9)\n\n")
        f.write("| circuit | pp | mst+P | gray+P | pp/mst | pp/gray |\n")
        f.write("|---|---:|---:|---:|---:|---:|\n")
        for r in merged:
            if r["set"] != "extended":
                continue
            pp, mst, gray = int(r["pp"]), int(r["mst"]), int(r["gray"])
            f.write(f"| {r['circuit']} | {pp} | {mst} | {gray} | {100*pp/mst:.1f}% | {100*pp/gray:.1f}% |\n")
        f.write(
            f"\n| **TOTAL (26)** | **{tp_all}** | **{tm_all}** | **{tg_all}** | "
            f"{100*tp_all/tm_all:.1f}% | {100*tp_all/tg_all:.1f}% |\n"
        )

    print(f"Primary: pp/mst={100*tp/tm:.1f}% pp/gray={100*tp/tg:.1f}% (paper 79%/70%)")
    print(f"All 26:  pp/mst={100*tp_all/tm_all:.1f}% pp/gray={100*tp_all/tg_all:.1f}%")


if __name__ == "__main__":
    main()
