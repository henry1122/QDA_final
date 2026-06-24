# TODD+PP Results (block-level Todd + pp+)

## Paper vs our pipeline

| | PhasePoly paper Table 1 | Our TODD+PP |
|--|--|--|
| Todd preprocessing | **No** (not in main pp/mst/gray comparison) | **Yes** — `todd_optimize_problem()` per block |
| Synthesis | pp+ (A* + PMH for O) | Same after Todd |
| Metric | CNOT only | CNOT only |

## Setup

- **pp+**: `benchmark_phasepoly` — single-block extraction, PhasePoly A*, PMH for O, `max_rz=30` (MST fallback if block has >30 Rz terms).
- **TODD+PP**: block-level Todd (Amy et al. 2014, qsyn `ToddPhasePolynomialOptimizationStrategy`) **then** same pp+ synthesis. Implemented via `--block-todd` in `tools/benchmark_phasepoly.cpp`.

## Primary result (19-circuit Table 1 suite)

Source: `results/benchmark_20260613_205458.txt` (pp+ baseline).

| Method | pp+ CNOT | MST+P CNOT | pp/mst | pp/gray | MST fallback (*)? |
|--------|----------:|-----------:|-------:|--------:|-------------------|
| **pp+** | 3869 | 5638 | 68.6% | 69.5% | **None** (no `*` on any row) |
| **TODD+PP** | **3869** | **5638** | **68.6%** | **69.5%** | **None** |
| Paper (ref.) | 3869 | 4891 | 79% | 70% | — |

**Interpretation:** On this suite, **no extracted block exceeded `max_rz=30`**, so PhasePoly A* always ran; Todd never needed to “unlock” a block from MST fallback. CNOT totals therefore **match pp+** exactly.

The `todd=on` flag in older benchmark logs only enabled the separate **full-circuit `todd+naive`** column — it did **not** apply block-level Todd to the pp column.

## Per-circuit CNOT (19 circuits)

| circuit | pp+ | MST+P | pp/mst | pp/gray | Notes |
|---------|----:|------:|-------:|--------:|-------|
| tof_3_pyzx | 16 | 28 | 57.1% | 61.5% | |
| tof_4_pyzx | 25 | 45 | 55.6% | 58.1% | |
| tof_5_pyzx | 35 | 63 | 55.6% | 58.3% | |
| tof_10_pyzx | 85 | 153 | 55.6% | 58.6% | |
| barenco_tof_3_pyzx | 21 | 31 | 67.7% | 72.4% | |
| barenco_tof_4_pyzx | 40 | 58 | 69.0% | 72.7% | |
| barenco_tof_5_pyzx | 60 | 82 | 73.2% | 74.1% | |
| barenco_tof_10_pyzx | 161 | 212 | 75.9% | 77.8% | |
| grover_5_pyzx | 262 | 342 | 76.6% | 77.1% | |
| ham15-low_pyzx | 227 | 293 | 77.5% | 77.2% | |
| ham15-med_pyzx | 440 | 602 | 73.1% | 73.1% | |
| ham15-high_pyzx | 1798 | 2649 | 67.9% | 68.8% | largest circuit |
| mod5_4_pyzx | 17 | 24 | 70.8% | 77.3% | |
| mod_mult_55_pyzx | 45 | 74 | 60.8% | 63.4% | |
| mod_red_21_pyzx | 89 | 129 | 69.0% | 68.5% | |
| qcla_com_7_pyzx | 137 | 204 | 67.2% | 64.6% | |
| qcla_mod_7_pyzx | 301 | 484 | 62.2% | 62.7% | |
| rc_adder_6_pyzx | 65 | 101 | 64.4% | 67.0% | |
| vbe_adder_3_pyzx | 45 | 64 | 70.3% | 71.4% | |
| **TOTAL** | **3869** | **5638** | **68.6%** | **69.5%** | |

## When TODD+PP is expected to beat pp+

1. **Blocks with Rz > `max_rz` after extraction** — Todd may reduce terms below 30 so A* replaces MST fallback.
2. **Blocks already under the cap** — smaller P after Todd can still lower CNOT, but on Table 1 we measured **no net change** at circuit totals.

Circuits excluded from the default Table 1 suite (e.g. `gf2^*`, `hwb8`) are more likely to hit the cap; see extended multiblock runs in `doc/table1_results_phase.csv`.

## Reproduce TODD+PP column

```bash
cmake --build build --target benchmark-phasepoly
./build/benchmark-phasepoly --no-todd --block-todd benchmark/qc/optimized/<circuits>.qc
```

Compare to pp+ baseline (no `--block-todd`).
