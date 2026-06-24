#!/usr/bin/env python3
"""Regenerate doc/table1_results.md from doc/table1_results.csv."""
import csv
from pathlib import Path

DOC = Path(__file__).resolve().parents[1] / "doc"
CSV = DOC / "table1_results.csv"
MD = DOC / "table1_results.md"

with open(CSV, encoding="utf-8-sig") as f:
    lines = [line for line in f if not line.startswith("#")]
rows = list(csv.DictReader(lines))
ok = [r for r in rows if r.get("ok") == "1"]

tp = sum(int(r["pp"]) for r in ok)
tm = sum(int(r["mst"]) for r in ok)
tg = sum(int(r["gray"]) for r in ok)
tgs = sum(int(r["gstair"]) for r in ok)
tn = sum(int(r["naive"]) for r in ok)

with open(MD, "w", encoding="utf-8") as f:
    f.write(
        "# gate metric: Clifford + T-family on full synthesized circuit (paper Table 1; excludes Rz)\n\n"
    )
    f.write(
        "| circuit | q | blk | Rz | pp | mst+P | gstair+P | gray+P | naive+P | pp/mst | pp/gray |\n"
    )
    f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in rows:
        if r.get("ok") != "1":
            f.write(f"| {r['circuit']} | ERR | | | | | | | | | |\n")
            continue
        pp, mst, gray = int(r["pp"]), int(r["mst"]), int(r["gray"])
        pm = f"{100 * pp / mst:.1f}%" if mst else "-"
        pg = f"{100 * pp / gray:.1f}%" if gray else "-"
        f.write(
            f"| {r['circuit']} | {r['q']} | {r['blk']} | {r['Rz']} | {pp} | {mst} | "
            f"{r['gstair']} | {gray} | {r['naive']} | {pm} | {pg} |\n"
        )
    f.write(
        f"\n| **TOTAL** ({len(ok)} ok) | | | | **{tp}** | **{tm}** | "
        f"**{tgs}** | **{tg}** | **{tn}** | "
    )
    if tm and tg:
        f.write(f"{100 * tp / tm:.0f}% | {100 * tp / tg:.0f}% |\n")
    else:
        f.write("| - | - |\n")
    f.write(
        "\nPaper reference (CNOT-only block metric): pp=3869, mst+P=4891 (79%), gray+P=5511 (70%)\n"
    )

print(f"Wrote {MD}: pp={tp} mst={tm} gray={tg} -> {100*tp/tm:.0f}% / {100*tp/tg:.0f}%")
